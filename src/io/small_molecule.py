# Clipper plugin to UCSF ChimeraX
# Copyright (C) 2016-2019 Tristan Croll, University of Cambridge

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

# Note that this software makes use of modified versions of the Clipper, LibCCP4
# and MMDB libraries, as well as portions of the Intel Math Kernel Library. Each
# of these is redistributed under its own license terms.

'''
Headless extraction of small-molecule (core) CIF structures - primarily for the
Crystallography Open Database (COD) - into a plain dictionary for downstream
consumers (e.g. graph-ML force-field training). The model and its connectivity
come from ChimeraX's own corecif parser; the unit cell, space group and all
crystallographic quality metrics are recovered from the CIF (see crystal.py).

Most COD entries are model + metadata only (no reflections). Such entries are
returned with the structure-factor-derived fields (recomputed R-factor,
data/parameter ratio is still available from the refinement record) set to None;
the absence is reported, not treated as an error.
'''

from ..symmetry import crystal_symmetry_from_cif_file, _strip_su, _first_cif_field

# Nuclear (neutron) X-H bond lengths in Angstroms, keyed by the bonded heavy-atom
# element. Small-molecule X-ray refinement places H at the electron-density
# centroid, which sits ~0.1-0.2 A short of the nucleus along the bond. Extending
# to these values gives nuclear-normalized positions. Values from Allen & Bruno,
# Acta Cryst. (2010) B66, 380, with common fall-throughs.
_NUCLEAR_XH_LENGTHS = {
    'C': 1.083, 'N': 1.009, 'O': 0.983, 'B': 1.18,
    'S': 1.34, 'P': 1.40, 'SI': 1.48,
}
_DEFAULT_XH_LENGTH = 1.06

# Minimum number of measured reflections required to attempt the anisotropic-Gaussian +
# 20-term isotropic-spline scaling of Fcalc onto Fobs (scale_fcalc_to_fobs). That fit has
# ~27 parameters; with far fewer reflections it is singular and the degenerate Clipper
# ResolutionFn solve reads past its (empty) resolution bins - harmless garbage on Windows,
# a hard segfault on Linux. Some COD reflection files are truncated to a handful of rows
# (e.g. 2219385's _refln_ loop stops after a single line), so guard both consumers - the
# recomputed R-factor and the live map - and skip scaling below this floor.
_MIN_SCALING_REFLECTIONS = 50


def open_small_molecule_cif(session, path, file_name=None):
    '''Open a small-molecule CIF via ChimeraX's corecif parser. Returns the
    (single) AtomicStructure model, added to the session.

    The returned model has already had the corecif workarounds applied (see
    :func:`_prepare_corecif_model`): coordinates rebuilt in Clipper's orthogonal
    frame, and covalent connectivity repaired where corecif drops it on
    metal-coordinated atoms. Callers therefore never need to correct the frame or
    re-perceive bonds themselves.'''
    from chimerax.mmcif.corecif import open_corecif
    models, status = open_corecif(session, path, file_name=file_name)
    if not models:
        from chimerax.core.errors import UserError
        raise UserError('No small-molecule structure found in %r' % path)
    session.models.add(models)
    model = models[0]
    _prepare_corecif_model(session, model, path)
    return model


def extract_cod_structure(session, path, file_name=None, keep_model=False,
                          radiation='auto'):
    '''
    Open a small-molecule CIF and return the Mode-1 payload as a dict.

    Args:
        path: path to a (core) CIF file.
        keep_model: if False (default), the opened model is closed before
            returning (headless batch use). Set True to leave it in the session.
        radiation: 'xray', 'electron', or 'auto' (default; read from the CIF's
            _diffrn_radiation_probe) - selects the scattering factors for the
            recomputed R-factor. Recorded in the payload as 'radiation'.
    '''
    model = open_small_molecule_cif(session, path, file_name=file_name)
    try:
        payload = _build_payload(session, model, path, radiation=radiation)
    finally:
        if not keep_model:
            session.models.close([model])
    return payload


def _build_payload(session, model, path, radiation='auto'):
    from chimerax.mmcif import get_cif_tables
    refine_t, diffrn_t, reflns_t, cell_t = get_cif_tables(
        path, ['refine', 'diffrn', 'reflns', 'cell'])

    radiation = _resolve_radiation(radiation, path)
    cell, spacegroup, grid = crystal_symmetry_from_cif_file(path)

    atoms = model.atoms
    # Coordinates were rebuilt into Clipper's frame at open (see
    # open_small_molecule_cif / _prepare_corecif_model); read them straight off.
    coords = atoms.coords
    elements = atoms.element_names.tolist()
    occupancies = atoms.occupancies

    payload = {
        'source_file': path,
        'space_group': spacegroup.symbol_hm,
        'space_group_number': spacegroup.spacegroup_number,
        'cell': (cell.a, cell.b, cell.c,
                 cell.alpha_deg, cell.beta_deg, cell.gamma_deg),
        'elements': elements,
        'coordinates': coords,    # Angstroms, Clipper orthogonal frame from CIF fractionals
        'connectivity': _connectivity(model),
        'occupancies': occupancies,
        'disorder': _disorder_flags(atoms),
        'resolution': _resolution(diffrn_t, cell_t, reflns_t),
        'collection_temperature': _temperature(diffrn_t, cell_t),
        'data_parameter_ratio': _data_parameter_ratio(refine_t),
        'published_r_factor': _published_r(refine_t),
        'n_restraints': _int_field(refine_t, 'ls_number_restraints'),
        'hydrogen_treatment': _first_cif_field(refine_t, 'ls_hydrogen_treatment'),
        'hydrogens': _hydrogen_positions(model, coords),
        # Radiation used for the recomputed R-factor: 'xray' or 'electron' (micro-ED).
        'radiation': radiation,
        # Structure-factor-derived fields: populated only when reflections are
        # available (see _structure_factor_metrics). None otherwise.
        'recomputed_r_factor': None,
        'has_structure_factors': False,
    }

    sf = _structure_factor_metrics(session, model, path, cell, spacegroup, grid, radiation)
    if sf is not None:
        payload.update(sf)
    return payload


def _clipper_frame_coords(model, path, cell):
    '''
    Cartesian coordinates (Angstroms) for the model's atoms, rebuilt from the CIF
    fractional coordinates via Clipper's orthogonalisation, in the model's atom
    order (matched by atom label). This avoids a fractional->Cartesian error in
    ChimeraX's corecif parser that distorts coordinates for oblique unit cells
    (see CORECIF_BUG_REPORT.md): rebuilt coordinates reproduce the CIF's published
    bond distances essentially exactly. Falls back to the corecif coordinate for
    any atom whose label is not found in the CIF.
    '''
    import numpy
    from chimerax.mmcif import get_cif_tables
    from ..clipper_python import Coord_frac

    atoms = model.atoms
    coords = numpy.array(atoms.coords, numpy.double)
    at, _aniso = get_cif_tables(path, ['atom_site', 'atom_site_aniso'])
    if at is None or not at.has_field('fract_x') or not at.has_field('label'):
        return coords
    frac = {}
    for r in at.fields(('label', 'fract_x', 'fract_y', 'fract_z')):
        try:
            frac[r[0]] = (float(_strip_su(r[1])), float(_strip_su(r[2])), float(_strip_su(r[3])))
        except ValueError:
            pass
    for i, name in enumerate(atoms.names):
        uvw = frac.get(name)
        if uvw is not None:
            co = Coord_frac(*uvw).coord_orth(cell)
            coords[i] = (co.x, co.y, co.z)
    return coords


# Generous per-element covalent-degree ceilings for the connectivity repair's
# valence cap. Kept in step with the downstream valence-sanity check
# (garnet-isolde categorize._MAX_DEGREE) so a repaired atom is never then flagged
# over-valent. Metals are omitted: their connectivity stays as coordination
# pseudobonds (see repair_connectivity), so they get no covalent cap here.
_MAX_COVALENT_DEGREE = {
    'C': 4, 'N': 4, 'O': 3, 'H': 1, 'F': 1, 'B': 4, 'Si': 6, 'P': 6,
    'S': 6, 'Se': 6, 'Cl': 4, 'Br': 4, 'I': 7, 'Te': 6, 'As': 5,
}
_DEFAULT_MAX_COVALENT_DEGREE = 8


def _neighbour_pairs(coords, radius):
    '''Indices (i, j), i < j, of atom pairs within ``radius`` Angstroms.'''
    import numpy
    try:
        from scipy.spatial import cKDTree
        return cKDTree(coords).query_pairs(radius, output_type='ndarray')
    except Exception:
        n = len(coords)
        out = []
        for i in range(n):
            d = numpy.linalg.norm(coords[i + 1:] - coords[i], axis=1)
            for k in numpy.nonzero(d <= radius)[0]:
                out.append((i, i + 1 + int(k)))
        return out


def repair_connectivity(session, model, tolerance=0.4, min_bond=0.6):
    '''
    Repair covalent connectivity that ChimeraX's corecif parser drops on
    metal-coordinated atoms, and re-perceive covalent bonds corecif skipped.

    corecif's ``finished_parse`` runs distance-based bond perception ONLY on atoms
    that do not already carry a bond, then demotes metal bonds to coordination
    pseudobonds. When a CIF's ``_geom_bond`` loop tabulates only the metal
    coordination sphere - a common organometallic convention - each ligand atom
    receives its single bond to the metal from ``_geom_bond``, is therefore skipped
    by distance perception, and then loses that bond to the pseudobond conversion:
    it ends up with NO covalent bond at all (e.g. a cyclopentadienyl ring is left
    completely unconnected). See the ChimeraX corecif ``finished_parse`` /
    ``connect_residue_by_distance`` interaction (reported upstream 2026-07).

    This re-perceives covalent bonds among NON-METAL atoms using ChimeraX's own
    criterion (``Element.bond_length(e1, e2) + tolerance``), adding any missing
    bond, with an H-bonds-once rule and a per-element valence cap
    (``_MAX_COVALENT_DEGREE``) so disordered/overlapping structures cannot be
    over-bonded. Pairs closer than ``min_bond`` Angstroms are skipped: no real
    covalent bond is that short (shortest is H2 at 0.74 A, and H-H is never perceived
    here), so a near-zero separation means coincident atoms - substitutional/occupancy
    disorder or duplicate sites - which must not be bonded (doing so manufactures
    zero-length bonds). Metal-ligand connectivity is intentionally left as coordination
    pseudobonds (not forced into covalent bonds). Coordinates should already be in
    Clipper's frame (see :func:`_prepare_corecif_model`). Returns the number of
    bonds added.
    '''
    import numpy
    from chimerax.atomic import Element
    from chimerax.atomic.struct_edit import add_bond

    atoms = model.atoms
    n = len(atoms)
    if n < 2:
        return 0
    coords = numpy.asarray(atoms.coords, numpy.double)
    elements = [a.element for a in atoms]
    enames = [e.name for e in elements]
    is_metal = numpy.array([e.is_metal for e in elements], bool)
    znum = numpy.array([e.number for e in elements])

    idx = {a: i for i, a in enumerate(atoms)}
    existing = set()
    degree = numpy.zeros(n, int)
    for b in model.bonds:
        i, j = idx[b.atoms[0]], idx[b.atoms[1]]
        existing.add((min(i, j), max(i, j)))
        degree[i] += 1
        degree[j] += 1

    # Candidate non-metal pairs within a generous global cutoff, filtered per pair
    # by the element-specific bond length; shortest first so that under a valence
    # cap the most credible (shortest) bonds win.
    candidates = []
    for i, j in _neighbour_pairs(coords, 3.2):
        i, j = int(i), int(j)
        if is_metal[i] or is_metal[j]:
            continue
        if (i, j) in existing:
            continue
        bl = Element.bond_length(elements[i], elements[j])
        if bl == 0.0:
            continue
        d = float(numpy.linalg.norm(coords[i] - coords[j]))
        if d < min_bond or d > bl + tolerance:   # too short = coincident atoms, not a bond
            continue
        candidates.append((d, i, j))
    candidates.sort()

    def cap(i):
        return _MAX_COVALENT_DEGREE.get(enames[i], _DEFAULT_MAX_COVALENT_DEGREE)

    added = 0
    for d, i, j in candidates:
        if znum[i] <= 1 and degree[i] >= 1:      # H (or lighter) bonds only once
            continue
        if znum[j] <= 1 and degree[j] >= 1:
            continue
        if degree[i] >= cap(i) or degree[j] >= cap(j):
            continue
        add_bond(atoms[i], atoms[j])
        degree[i] += 1
        degree[j] += 1
        added += 1
    return added


def reassemble_symmetry_scattered_hydrogens(session, model, cell, spacegroup,
                                            search_cells=1, tolerance=0.45):
    '''
    Relocate hydrogens that corecif left orphaned because their bond to a heavy atom is
    symmetry-coded in the CIF's ``_geom_bond`` loop. corecif does not create
    symmetry-coded bonds, so such an H is placed at its deposited (asymmetric-unit)
    coordinate - which lies in a different ASU / unit cell from the heavy atom it
    actually bonds - and comes out with no covalent bond (a "scattered" hydrogen).

    For each still-unbonded H, this searches every space-group symop combined with
    lattice translations over a ``(2*search_cells+1)**3`` block of cells for the image
    that lands within X-H bonding distance (``Element.bond_length(H, X) + tolerance``)
    of a heavy atom. If EXACTLY ONE heavy partner is found (an unambiguous match), the H
    is moved onto that image and bonded to the partner; hydrogens with no match or an
    ambiguous (multi-partner) match are left untouched. The recovered operator matches
    the CIF's own ``_geom_bond`` symmetry code where one is present, and the search also
    handles entries with no ``_geom_bond`` loop.

    OPT-IN and MUTATING: unlike :func:`repair_connectivity` (which only adds bonds and
    runs at open), this MOVES atoms, so call it deliberately - e.g. when assembling clean
    covalent fragments for a small-molecule entry - not on every open. Restricted to
    hydrogen on purpose: a symmetry-related bond between heavy atoms is usually genuine
    extended structure (a coordination polymer), not a placement error. Run it AFTER
    :func:`repair_connectivity` so only genuinely scattered H remain. Returns the number
    of hydrogens relocated.
    '''
    import numpy
    from chimerax.atomic import Element
    from chimerax.atomic.struct_edit import add_bond
    from ..sym_realize import unit_cell_places

    atoms = model.atoms
    coords = numpy.asarray(atoms.coords, numpy.double)
    elements = [a.element for a in atoms]
    znum = numpy.array([e.number for e in elements])
    degree = numpy.zeros(len(atoms), int)
    idx = {a: i for i, a in enumerate(atoms)}
    for b in model.bonds:
        degree[idx[b.atoms[0]]] += 1
        degree[idx[b.atoms[1]]] += 1

    # Candidate parents are NON-METAL heavy atoms only: a hydride (H bonded solely to a
    # metal) must stay a metal-coordination interaction, not be forged into a covalent
    # H-metal bond, so it finds no partner here and is left untouched.
    is_metal = numpy.array([e.is_metal for e in elements])
    heavy = numpy.nonzero((znum > 1) & ~is_metal)[0]
    if len(heavy) == 0:
        return 0
    heavy_xyz = coords[heavy]
    Hn = Element.get_element('H')
    cut = numpy.array([Element.bond_length(Hn, elements[int(k)]) + tolerance for k in heavy])
    float_h = [i for i in range(len(atoms)) if znum[i] == 1 and degree[i] == 0]
    if not float_h:
        return 0

    n = 2 * int(search_cells) + 1
    places = unit_cell_places(cell, spacegroup, n, n, n,
                              origin=(-search_cells, -search_cells, -search_cells))

    relocated = 0
    for hi in float_h:
        h_xyz = coords[hi:hi + 1]
        best = None                     # (dist, target_xyz, heavy_local_index)
        partners = set()
        for P in places:
            t = P.transform_points(h_xyz)[0]
            d = numpy.linalg.norm(heavy_xyz - t, axis=1)
            k = int(numpy.argmin(d))
            if d[k] <= cut[k]:
                partners.add(k)
                if best is None or d[k] < best[0]:
                    best = (float(d[k]), t, k)
        if best is None or len(partners) != 1:
            continue                    # no match, or ambiguous -> leave it alone
        # clash guard: don't drop the H onto an existing atom (a duplicate)
        if numpy.min(numpy.linalg.norm(coords - best[1], axis=1)) < 0.3:
            continue
        atoms[hi].coord = best[1]
        add_bond(atoms[hi], atoms[int(heavy[best[2]])])
        relocated += 1
    return relocated


def _prepare_corecif_model(session, model, path, cell=None):
    '''
    Apply the corecif workarounds every freshly-opened small-molecule model needs,
    so outside callers get a corrected model straight from
    :func:`open_small_molecule_cif` and never touch the workarounds themselves:

      1. rebuild coordinates in Clipper's orthogonal frame (corecif
         mis-orthogonalises oblique cells - see :func:`_clipper_frame_coords`);
      2. repair covalent connectivity corecif drops on metal-coordinated atoms
         (see :func:`repair_connectivity`).

    Both steps are idempotent. A no-op when the CIF carries no crystal symmetry.
    '''
    if cell is None:
        try:
            cell, _sg, _grid = crystal_symmetry_from_cif_file(path)
        except Exception:
            return
    try:
        model.atoms.coords = _clipper_frame_coords(model, path, cell)
    except Exception as e:
        session.logger.warning('(CLIPPER) coordinate-frame correction failed for '
                               '%r: %s' % (path, e))
    try:
        repair_connectivity(session, model)
    except Exception as e:
        session.logger.warning('(CLIPPER) connectivity repair failed for %r: %s'
                               % (path, e))


def _connectivity(model):
    # Bond list as pairs of indices into model.atoms (ChimeraX atom order).
    atoms = model.atoms
    bond_atoms = model.bonds.atoms
    i0 = atoms.indices(bond_atoms[0])
    i1 = atoms.indices(bond_atoms[1])
    return list(zip(i0.tolist(), i1.tolist()))


def _disorder_flags(atoms):
    # An atom is flagged disordered if it carries an alternate-location id or a
    # partial occupancy.
    import numpy
    alt = atoms.alt_locs
    partial = atoms.occupancies < 1.0
    disordered = numpy.array([a not in ('', ' ', '\x00') for a in alt]) | partial
    return disordered


def _hydrogen_positions(model, coords):
    '''
    Return per-hydrogen positions:
        electron_density: the deposited coordinate (meaningful only where H were
            independently refined - see hydrogen_treatment flag).
        nuclear: the deposited coordinate extended along its bond to the standard
            nuclear (neutron) X-H distance for the bonded heavy atom.
    Both as {atom_index: (x,y,z)} dicts; only H bonded to exactly one heavy atom
    are normalizable (others are reported with electron_density only).
    '''
    import numpy
    atoms = model.atoms
    elements = atoms.element_names
    index_of = {a: i for i, a in enumerate(atoms)}
    h_idx = [i for i, e in enumerate(elements) if e == 'H']
    electron_density = {}
    nuclear = {}
    for i in h_idx:
        h_atom = atoms[i]
        electron_density[i] = tuple(float(x) for x in coords[i])
        heavy = [n for n in h_atom.neighbors if n.element.name != 'H']
        if len(heavy) != 1 or heavy[0] not in index_of:
            continue
        hv_xyz = coords[index_of[heavy[0]]]
        vec = coords[i] - hv_xyz
        d = numpy.linalg.norm(vec)
        if d < 1e-3:
            continue
        target = _NUCLEAR_XH_LENGTHS.get(heavy[0].element.name.upper(), _DEFAULT_XH_LENGTH)
        nuc = hv_xyz + vec / d * target
        nuclear[i] = tuple(float(x) for x in nuc)
    return {'electron_density': electron_density, 'nuclear': nuclear}


# ---- crystallographic quality metrics (read from the refinement record) ----

def _published_r(refine_t):
    # Prefer the all-data R; fall back to observed/gt subsets.
    for field in ('ls_R_factor_all', 'ls_R_factor_obs', 'ls_R_factor_gt'):
        v = _first_cif_field(refine_t, field)
        if v is not None:
            try:
                return float(_strip_su(v))
            except ValueError:
                pass
    return None


def _data_parameter_ratio(refine_t):
    nrefl = _int_field(refine_t, 'ls_number_reflns')
    nparm = _int_field(refine_t, 'ls_number_parameters')
    if nrefl and nparm:
        return nrefl / nparm
    return None


def _temperature(diffrn_t, cell_t):
    v = _first_cif_field(diffrn_t, 'ambient_temperature')
    if v is None:
        v = _first_cif_field(cell_t, 'measurement_temperature')
    if v is not None:
        try:
            return float(_strip_su(v))
        except ValueError:
            pass
    return None


def _resolution(diffrn_t, cell_t, reflns_t):
    # Reported resolution for the payload metadata. Only physically meaningful
    # (positive) values are returned; a non-positive wavelength/theta (COD's
    # '-1'/'0' sentinels) or a garbage d_min is reported as None ("unknown")
    # rather than a negative number poisoning downstream consumers.
    from math import sin, radians
    v = _first_cif_field(reflns_t, 'd_resolution_high')
    if v is not None:
        try:
            d = float(_strip_su(v))
            if d > 0:
                return d
        except ValueError:
            pass
    wl = _first_cif_field(diffrn_t, 'radiation_wavelength')
    theta = _first_cif_field(diffrn_t, 'reflns_theta_max')
    if theta is None:
        theta = _first_cif_field(cell_t, 'measurement_theta_max')
    try:
        if wl is not None and theta is not None:
            wl = float(_strip_su(wl))
            sin_theta = sin(radians(float(_strip_su(theta))))
            if wl > 0 and sin_theta > 0:
                return wl / (2.0 * sin_theta)
    except (ValueError, ZeroDivisionError):
        pass
    return None


def _int_field(table, field):
    v = _first_cif_field(table, field)
    if v is not None:
        try:
            return int(_strip_su(v))
        except ValueError:
            pass
    return None


def _find_reflection_table(path):
    '''
    Locate experimental reflections for a structure. COD stores them in a sibling
    CIF reflection file with a .hkl extension (itself a CIF carrying a _refln_
    loop); some CIFs embed the loop directly. Returns the CIFTable for 'refln',
    or None.
    '''
    import os
    from chimerax.mmcif import get_cif_tables
    candidates = []
    base, _ = os.path.splitext(path)
    candidates.append(base + '.hkl')
    candidates.append(path)               # embedded _refln_ in the model file
    for cand in candidates:
        if not os.path.exists(cand):
            continue
        tables = get_cif_tables(cand, ['refln'])
        refln = tables[0] if tables else None
        if refln is not None and refln.has_field('index_h') and refln.has_field('F_squared_meas'):
            return refln
    return None


def _radiation_enum(radiation):
    '''Map a radiation spec ('xray'/'electron') to the Clipper AtomShapeFn.RADIATION
    enum used by the SFcalc engines. Anything unrecognised is treated as X-ray.'''
    from ..clipper_python import AtomShapeFn
    if str(radiation).lower() in ('electron', 'ed', 'microed', '3ded', 'e'):
        return AtomShapeFn.ELECTRON
    return AtomShapeFn.XRAY


def _radiation_from_cif(path, default='xray'):
    '''Detect the experiment radiation from the CIF: 'electron' when
    _diffrn_radiation_probe (or _type) names electrons, else the default. Used to
    auto-select electron scattering factors for micro-ED / 3D-ED data.'''
    from chimerax.mmcif import get_cif_tables
    # get_cif_tables returns [] (not [empty_table]) when the file contains NONE of the
    # requested categories - which happens for entries with no _diffrn_radiation_* data.
    tables = get_cif_tables(path, ['diffrn_radiation'])
    dr = tables[0] if tables else None
    for field in ('probe', 'type'):
        v = _first_cif_field(dr, field)
        if v and 'electron' in v.lower():
            return 'electron'
    return default


def _resolve_radiation(radiation, path):
    '''Resolve 'auto' to the CIF-detected radiation; pass an explicit choice through.'''
    if radiation is None or str(radiation).lower() == 'auto':
        return _radiation_from_cif(path)
    return str(radiation).lower()


def read_small_molecule_fobs(path, cell, spacegroup):
    '''
    Read a small-molecule reflection file (a COD ``.hkl`` / CIF ``_refln_`` loop of
    ``F_squared_meas``) into the observed-data HKL_data objects a reciprocal
    differentiable target needs, in the small-molecule convention:

      * ``fobs`` (``HKL_data<F_sigF>``): observed amplitude ``Fo = sqrt(I)`` for
        ``I > 0`` (else 0), ``sigF = sigma(I) / (2 sqrt(I))`` where ``sigma(I)`` is
        present, else 1;
      * ``phi_fom`` (``HKL_data<Phi_fom>``) with ``fom() == 1`` — small-molecule data
        carry no sigma-A weighting, so the target reduces to a pure 1/sig^2 fit;
      * ``usage`` (``HKL_data<Flag>``) all-working — small-molecule datasets have no
        free set.

    Returns ``(hkl_info, fobs, phi_fom, usage)``. The resolution limit of the
    generated reflection list is taken from the measured reflections. NOTE: Clipper
    ``HKL_data`` holds a *non-owning* pointer to ``hkl_info``, so the caller MUST keep
    the returned ``hkl_info`` alive for as long as any of the HKL_data objects (or a
    target built from them) is in use — otherwise it dangles.

    This is the headless reflection half of the COD -> structure-factors ->
    differentiable-target pipeline; see
    :func:`chimerax.clipper.diff.crystal.small_molecule_ensemble_target`. It is a thin
    composition of :func:`_parse_reflection_file` (the file read) and
    :func:`fobs_from_arrays` (the pure numpy->Clipper transform); a caller that has
    cached the raw ``(hkl, fsq, sig)`` arrays can skip the file and call
    :func:`fobs_from_arrays` directly.
    '''
    hkl, fsq, sig = _parse_reflection_file(path)
    return fobs_from_arrays(hkl, fsq, sig, cell, spacegroup)


def _parse_reflection_file(path):
    '''
    Parse a small-molecule reflection file (COD ``.hkl`` / CIF ``_refln_`` loop) into
    the raw post-parse columns ``(hkl, fsq, sig)``: Miller indices ``int32[N,3]``,
    ``F_squared_meas`` (intensity I) ``float64[N]``, and ``sigma(F_squared)``
    ``float64[N]`` (1.0 where absent). No transform is applied — this is the
    file-reading half of :func:`read_small_molecule_fobs`, split out so a caller that
    has cached these arrays (e.g. the COD structure-factor cache) can rebuild the
    observed data via :func:`fobs_from_arrays` with no file access. The
    corrupt/short-loop guards live here, so cached arrays inherit the vetted output.
    Raises ``UserError`` when no reflections are present.
    '''
    import numpy
    from chimerax.core.errors import UserError

    refln = _find_reflection_table(path)
    if refln is None:
        raise UserError('No reflections (F_squared_meas) found for %r' % path)

    rows = refln.fields(('index_h', 'index_k', 'index_l',
                         'F_squared_meas', 'F_squared_sigma'),
                        allow_missing_fields=True)
    hkl = numpy.array([[int(r[0]), int(r[1]), int(r[2])] for r in rows], numpy.int32)
    fsq = numpy.array([float(_strip_su(r[3])) for r in rows], numpy.double)
    sig = numpy.array([float(_strip_su(r[4]))
                       if (len(r) > 4 and r[4] not in ('', '?', '.')) else 1.0
                       for r in rows], numpy.double)
    if len(hkl) == 0:
        raise UserError('Reflection loop for %r is empty' % path)
    return hkl, fsq, sig


def fobs_from_arrays(hkl, fsq, sig, cell, spacegroup):
    '''
    Build the observed-data HKL_data objects a reciprocal differentiable target needs
    from raw reflection arrays (the pure numpy->Clipper half of
    :func:`read_small_molecule_fobs`), in the small-molecule convention. Usable on
    cached arrays with no file access — the intended entry point for the COD
    structure-factor cache loader.

    Args:
        * hkl: Miller indices, integer array shaped ``(N, 3)``.
        * fsq: ``F_squared_meas`` (intensity I), ``(N,)``.
        * sig: ``sigma(F_squared)``, ``(N,)`` (use 1.0 where unknown).
        * cell / spacegroup: the crystal ``Cell`` / ``Spacegroup``.

    Returns ``(hkl_info, fobs, phi_fom, usage)``:
      * ``fobs`` (``HKL_data<F_sigF>``): ``Fo = sqrt(I)`` for ``I > 0`` (else 0),
        ``sigF = sigma(I) / (2 sqrt(I))``;
      * ``phi_fom`` (``HKL_data<Phi_fom>``) with ``fom() == 1`` (no sigma-A weighting);
      * ``usage`` (``HKL_data<Flag>``) all-working (small-molecule data have no free set).

    The resolution limit is taken from the highest-resolution supplied reflection.
    NOTE: Clipper ``HKL_data`` holds a *non-owning* pointer to ``hkl_info``, so the
    caller MUST keep the returned ``hkl_info`` alive for as long as any of the HKL_data
    objects (or a target built from them) is in use — otherwise it dangles (segfault).
    '''
    import numpy
    from chimerax.core.errors import UserError
    from .. import HKL_info, HKL_data_F_sigF, HKL_data_Phi_fom, HKL_data_Flag
    from ..clipper_python import HKL, Resolution

    hkls = numpy.ascontiguousarray(hkl, numpy.int32).reshape(-1, 3)
    fsq = numpy.ascontiguousarray(fsq, numpy.double).reshape(-1)
    sig = numpy.ascontiguousarray(sig, numpy.double).reshape(-1)
    n = len(hkls)
    if n == 0 or len(fsq) != n or len(sig) != n:
        raise UserError('fobs_from_arrays: hkl/fsq/sig length mismatch or empty '
                        '(%d/%d/%d)' % (n, len(fsq), len(sig)))

    valid = fsq > 0.0
    fo = numpy.zeros(n, numpy.double)
    sigf = numpy.ones(n, numpy.double)
    fo[valid] = numpy.sqrt(fsq[valid])
    sigf[valid] = sig[valid] / (2.0 * numpy.sqrt(fsq[valid]))

    # Resolution from the measured reflections, so the generated ASU covers them all.
    invresolsq = numpy.array([HKL(hkls[i].tolist()).invresolsq(cell) for i in range(n)])
    res = 1.0 / numpy.sqrt(invresolsq.max())
    hkl_info = HKL_info(spacegroup, cell, Resolution(res), True)

    fobs = HKL_data_F_sigF(hkl_info)
    fobs.set_data(hkls, numpy.stack([fo, sigf], axis=1).astype(numpy.float32))
    phi_fom = HKL_data_Phi_fom(hkl_info)
    phi_fom.set_data(hkls, numpy.stack(
        [numpy.zeros(n, numpy.float32), numpy.ones(n, numpy.float32)], axis=1))
    usage = HKL_data_Flag(hkl_info)
    usage.set_data(hkls, numpy.ones((n, 1), numpy.int32))
    return hkl_info, fobs, phi_fom, usage


def _structure_factor_metrics(session, model, path, cell, spacegroup, grid, radiation='xray'):
    '''
    If experimental structure factors are available (sibling .hkl reflection CIF,
    or an embedded _refln_ loop), compute a bulk-solvent-free R-factor of the
    model against the observed amplitudes and return the SF-derived fields.
    Returns None when no reflections are available (the common COD case).

    Small-molecule convention: Fo = sqrt(F_squared_meas); the headline R is
    R = sum||Fo| - |Fc_scaled|| / sum|Fo| over ALL reflections (a single R over all
    data, matching the live map and the published small-molecule R), with the
    I > 2 sigma observed-subset value reported alongside. Fcalc is placed on the Fo
    scale by the same anisotropic-Gaussian + isotropic-spline model the live map
    (SmallMoleculeXmapMgr) and the published R use — NOT a single linear scale, which
    piles the residual onto the heaviest scatterer and biases R low on metal-containing
    entries. (Not French-Wilson, a macromolecular technique that would not match the
    published small-molecule R.)
    '''
    refln = _find_reflection_table(path)
    if refln is None:
        return None

    import numpy
    rows = refln.fields(('index_h', 'index_k', 'index_l',
                         'F_squared_meas', 'F_squared_sigma'),
                        allow_missing_fields=True)
    hkls = numpy.array([[int(r[0]), int(r[1]), int(r[2])] for r in rows], numpy.int32)
    fsq = numpy.array([float(_strip_su(r[3])) for r in rows], numpy.double)
    # sigma(F^2) may be absent from some reflection files; default to 1.0 (the
    # I > 2 sigma observed cut then has no effect and R is reported over all data).
    sig = numpy.array([float(_strip_su(r[4]))
                       if (len(r) > 4 and r[4] not in ('', '?', '.')) else 1.0
                       for r in rows], numpy.double)

    # Too few measured reflections to support the aniso+spline scale (a truncated /
    # malformed COD _refln_ loop): skip it rather than feed a singular fit to Clipper's
    # ResolutionFn (see _MIN_SCALING_REFLECTIONS). Reflections exist, so has_sf is True.
    n_measured = int((fsq > 0).sum())
    if n_measured < _MIN_SCALING_REFLECTIONS:
        return {'has_structure_factors': True, 'recomputed_r_factor': None,
                'n_reflections_used': 0, 'n_reflections_measured': n_measured}

    from .. import HKL_info
    from ..clipper_python import HKL, Resolution, SFcalc_aniso_sum_float
    from ..clipper_python.data32 import HKL_data_F_phi_float, HKL_data_F_sigF_float
    from ..clipper_python.ext import scale_fcalc_to_fobs

    # Resolution limit from the measured reflections, so the generated ASU covers
    # every observation.
    invresolsq = numpy.array([HKL(hkls[i].tolist()).invresolsq(cell)
                              for i in range(len(hkls))])
    res = 1.0 / numpy.sqrt(invresolsq.max())
    hkl_info = HKL_info(spacegroup, cell, Resolution(res), True)

    # Observed amplitudes for the aniso+spline scaling: Fo = sqrt(I) (I > 0),
    # sigma(Fo) ~= sigma(I) / (2 sqrt(I)).
    valid = fsq > 0.0
    fo_amp = numpy.zeros(len(fsq), numpy.double)
    sig_amp = numpy.ones(len(fsq), numpy.double)
    fo_amp[valid] = numpy.sqrt(fsq[valid])
    sig_amp[valid] = sig[valid] / (2.0 * numpy.sqrt(fsq[valid]))
    fobs = HKL_data_F_sigF_float(hkl_info)
    fobs.set_data(hkls, numpy.stack([fo_amp, sig_amp], axis=1).astype(numpy.float32))

    atoms = _sfcalc_atom_list(path, cell, spacegroup, grid, radiation)
    if atoms is None:
        return {'has_structure_factors': True, 'recomputed_r_factor': None,
                'n_reflections_used': 0}

    # Exact direct summation of structure factors with NO bulk-solvent term
    # (correct for densely-packed small-molecule crystals; FFT is an unnecessary
    # approximation at this atom count). Build the functor then call it - the
    # "compute in constructor" shorthand does not fill fcalc through the bindings.
    fcalc = HKL_data_F_phi_float(hkl_info)
    SFcalc_aniso_sum_float(radiation=_radiation_enum(radiation))(fcalc, atoms)

    # Place Fcalc on the Fobs scale with the anisotropic-Gaussian + isotropic-spline
    # model - the SAME scaling the live small-molecule map (SmallMoleculeXmapMgr) and
    # the published R use, so headless R matches the GUI/published value. A single
    # linear scale k = sum(Fo Fc)/sum(Fc^2) piles the residual onto the heaviest
    # scatterer, biasing R low on metal-containing entries - the bug this replaces.
    scaled = HKL_data_F_phi_float(hkl_info)
    scale_fcalc_to_fobs(fcalc, fobs, scaled)

    # Gather Fo and scaled Fc by reference index (both indexed off hkl_info's ASU,
    # so no raw-HKL lookup that could miss the ASU). od.f = Fo, od.sigf = sigma(Fo),
    # scaled[ih.hkl].f = Fc on the Fo scale.
    fo_l, fc_l, sf_l = [], [], []
    ih = fobs.first_data
    while not ih.last():
        od = fobs[ih]
        if not od.missing and od.f > 0:
            cd = scaled[ih.hkl]
            if not cd.missing:
                fo_l.append(od.f); fc_l.append(cd.f); sf_l.append(od.sigf)
        fobs.next_data(ih)
    fo = numpy.array(fo_l); fc = numpy.array(fc_l); sf = numpy.array(sf_l)
    if not len(fo) or fc.sum() == 0:
        return {'has_structure_factors': True, 'recomputed_r_factor': None,
                'n_reflections_used': 0}

    # Observed set I > 2 sigma(I): with I = Fo^2 and sigma(I) = 2 Fo sigma(Fo) this is
    # Fo > 4 sigma(Fo). R uses the aniso+spline-scaled Fc (no manual scale factor) via the
    # shared R-factor core (reflection_tools.compute_r_factors) - the same routine the live
    # small-molecule map and the command log use, so the three surfaces cannot drift. No
    # epsilon weighting (small-molecule convention: unweighted, all reflections).
    from ..reflection_tools import compute_r_factors
    obs = fo > 4.0 * sf
    if not obs.any():
        obs = numpy.ones(len(fo), bool)
    rf = compute_r_factors(fo, fc, observed_mask=obs)
    # Headline recomputed_r_factor is over ALL reflections - matching the live map and the
    # published small-molecule R (a single R over all data) - so headless == GUI. The
    # I > 2 sigma observed value is retained under recomputed_r_factor_observed.
    return {
        'has_structure_factors': True,
        'recomputed_r_factor': float(rf.r_all),
        'recomputed_r_factor_all': float(rf.r_all),
        'recomputed_r_factor_observed': (float(rf.r_observed)
                                         if rf.r_observed is not None else None),
        'n_reflections_used': int(rf.n_all),
        'n_reflections_observed': int(rf.n_observed),
        'resolution_data': float(res),
    }


# Two-letter chemical element symbols, for parsing _atom_site_type_symbol values
# that may carry oxidation states or labels (e.g. "Cu2+", "O1-").
_TWO_LETTER_ELEMENTS = frozenset((
    'He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti',
    'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
    'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt',
    'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Th', 'Pa', 'Np', 'Pu',
))


def _element_from_type_symbol(symbol):
    leading = ''
    for c in symbol.strip().strip("'\""):
        if c.isalpha():
            leading += c
        else:
            break
    if len(leading) >= 2 and leading[:2].capitalize() in _TWO_LETTER_ELEMENTS:
        return leading[:2].capitalize()
    return leading[:1].upper() if leading else 'C'


def _sfcalc_atom_list(path, cell, spacegroup, grid, radiation='xray'):
    '''
    Build a Clipper Atom_list for structure-factor calculation directly from the
    CIF, entirely in CLIPPER's orthogonal frame. Built from the CIF (not the
    corecif AtomicStructure) so it is robust to disorder, where corecif may load a
    different atom set than the deposited atom_site list. Three things must be
    right for small-molecule crystals - none are handled by feeding ChimeraX
    Cartesian coordinates directly:

      * Coordinate frame: corecif orthogonalises with a different frame than
        Clipper, biasing F_calc; we rebuild Cartesian coordinates from the CIF
        fractional coordinates via Clipper's own Cell.
      * Special positions: Clipper's SFcalc applies every symop to every atom, so
        an atom on a special position is multiply counted; we divide its occupancy
        by the site multiplicity from Clipper's optimised Xmap.multiplicity.
      * Anisotropic ADPs: the CIF _atom_site_aniso_U values are referred to the
        reciprocal axes and must be scaled by a*_i a*_j before being interpreted as
        a Clipper U_aniso_frac, then transformed to the orthogonal frame.
    '''
    scaffold = sfcalc_scaffold(path, cell, spacegroup, grid, radiation)
    if scaffold is None:
        return None
    return atom_list_from_scaffold(scaffold)


def sfcalc_scaffold(path, cell, spacegroup, grid, radiation='xray'):
    '''
    Per-atom structure-factor inputs read once from the CIF, in Clipper's
    orthogonal frame: elements, Clipper-frame coordinates, special-position-
    corrected occupancies, isotropic U, and orthogonal anisotropic U. Also returns
    the atom labels (for aligning live model coordinates by label). The expensive
    parts (Coord_frac->orth, Xmap.multiplicity, ADP frame conversion) are done here
    so a live recalc only has to swap in fresh coordinates - see
    atom_list_from_scaffold.
    '''
    import numpy
    from chimerax.mmcif import get_cif_tables
    from ..clipper_python import Coord_frac, Xmap_double, U_aniso_frac

    at, an = get_cif_tables(path, ['atom_site', 'atom_site_aniso'])
    if at is None or not at.has_field('fract_x'):
        return None
    arows = at.fields(('label', 'type_symbol', 'fract_x', 'fract_y', 'fract_z',
                       'occupancy', 'U_iso_or_equiv'), allow_missing_fields=True)
    n = len(arows)
    if n == 0:
        return None

    aniso_map = {}
    if an is not None and an.has_field('U_11'):
        for r in an.fields(('label', 'U_11', 'U_22', 'U_33', 'U_12', 'U_13', 'U_23'),
                           allow_missing_fields=True):
            try:
                aniso_map[r[0]] = [float(_strip_su(x)) for x in r[1:7]]
            except (ValueError, IndexError):
                pass

    # CIF U_ij -> Clipper U_aniso_frac scaling (reciprocal-axis products).
    asx, bsx, csx = cell.a_star, cell.b_star, cell.c_star
    sc6 = numpy.array([asx * asx, bsx * bsx, csx * csx, asx * bsx, asx * csx, bsx * csx])

    labels = []
    elements = []
    coords = numpy.empty((n, 3), numpy.double)
    occupancies = numpy.empty(n, numpy.double)
    u_iso = numpy.zeros(n, numpy.double)
    u_aniso = numpy.ones((n, 6), numpy.double) * numpy.nan
    xm = Xmap_double(spacegroup, cell, grid)
    from ..scattering import clipper_species_from_type_symbol
    for i, r in enumerate(arows):
        labels.append(r[0])
        # Honour an ionic scattering type declared in the CIF type_symbol
        # (e.g. "O2-", "Ca2+"), validated against the table for this radiation
        # (X-ray = Waasmaier-Kirfel ions; electron = Peng-1998 ions, which add the
        # divergent Coulomb term in AtomShapeFn); neutral element otherwise.
        elements.append(clipper_species_from_type_symbol(
            r[1], _element_from_type_symbol(r[1]), radiation=radiation))
        cf = Coord_frac(float(_strip_su(r[2])), float(_strip_su(r[3])), float(_strip_su(r[4])))
        co = cf.coord_orth(cell)
        coords[i] = (co.x, co.y, co.z)
        occ = float(_strip_su(r[5])) if (len(r) > 5 and r[5] not in ('', '?', '.')) else 1.0
        occupancies[i] = occ / xm.multiplicity(cf.coord_grid(grid))
        if r[0] in aniso_map:
            uf = (numpy.array(aniso_map[r[0]]) * sc6).tolist()
            u_aniso[i] = U_aniso_frac(*uf).u_aniso_orth(cell).as_numpy()
        elif len(r) > 6 and r[6] not in ('', '?', '.'):
            u_iso[i] = float(_strip_su(r[6]))
        else:
            u_iso[i] = 0.05
    return {'labels': labels, 'elements': elements, 'coords': coords,
            'occupancies': occupancies, 'u_iso': u_iso, 'u_aniso': u_aniso}


def atom_list_from_scaffold(scaffold, coords=None):
    '''Build a Clipper Atom_list from a scaffold (see sfcalc_scaffold), optionally
    overriding the coordinates (Clipper orthogonal frame, n x 3) for a live recalc.'''
    from ..clipper_python import Atom_list
    xyz = scaffold['coords'] if coords is None else coords
    return Atom_list(scaffold['elements'], xyz, scaffold['occupancies'],
                     scaffold['u_iso'], scaffold['u_aniso'])


def register_clipper_atom_attributes(session):
    '''Register (once per session, idempotent) the per-atom attributes Clipper sets on
    small-molecule models. Both round-trip through .cxs sessions via ChimeraX's custom-
    attribute machinery; registering here means they exist before any session restore
    sets them.'''
    if getattr(session, '_clipper_atom_attrs_registered', False):
        return
    from chimerax.atomic import Atom
    Atom.register_attr(session, 'clipper_scattering_species', 'ChimeraX-Clipper', attr_type=str)
    Atom.register_attr(session, 'clipper_sf_exclude', 'ChimeraX-Clipper', attr_type=bool)
    session._clipper_atom_attrs_registered = True


def hydrate_small_molecule_model(session, model, path, cell, radiation='xray'):
    '''
    Give a freshly-opened corecif small-molecule model the crystallographic per-atom
    data corecif does not carry, so the live map can build its structure-factor inputs
    directly from the model (rather than from a frozen CIF scaffold):

      * isotropic B - corecif never parses _atom_site_U_iso_or_equiv, so Atom.bfactor
        is 0; set it (B = 8 pi^2 U), giving both the map input and a B the user can edit;
      * anisotropic U - corecif stores the raw reciprocal-fractional U_ij; overwrite with
        the orthogonal-frame tensor (a*-scaled, as sfcalc_scaffold computes it);
      * ionic scattering species - corecif discards the charge from _atom_site_type_symbol;
        store the Clipper species (honouring the CIF charge) on a custom attribute.

    Values are matched to atoms by CIF label; atoms whose label is absent keep defaults.
    '''
    import numpy
    from chimerax.mmcif import get_cif_tables
    from ..clipper_python import U_aniso_frac
    from ..scattering import clipper_species_from_type_symbol
    register_clipper_atom_attributes(session)

    at, an = get_cif_tables(path, ['atom_site', 'atom_site_aniso'])
    if at is None or not at.has_field('label'):
        return
    info = {}
    for r in at.fields(('label', 'type_symbol', 'U_iso_or_equiv'), allow_missing_fields=True):
        info[r[0]] = (r[1] if len(r) > 1 else '', r[2] if len(r) > 2 else '')
    aniso_raw = {}
    if an is not None and an.has_field('U_11'):
        for r in an.fields(('label', 'U_11', 'U_22', 'U_33', 'U_12', 'U_13', 'U_23'),
                           allow_missing_fields=True):
            try:
                aniso_raw[r[0]] = [float(_strip_su(x)) for x in r[1:7]]
            except (ValueError, IndexError):
                pass

    asx, bsx, csx = cell.a_star, cell.b_star, cell.c_star
    sc6 = numpy.array([asx*asx, bsx*bsx, csx*csx, asx*bsx, asx*csx, bsx*csx])
    b_from_u = 8.0 * numpy.pi ** 2

    atoms = model.atoms
    bfactors = numpy.array(atoms.bfactors, numpy.float32)
    aniso_atoms = []
    aniso_u6 = []
    for i, a in enumerate(atoms):
        ts, uiso = info.get(a.name, ('', ''))
        a.clipper_scattering_species = clipper_species_from_type_symbol(
            ts, a.element.name, radiation=radiation)
        if a.name in aniso_raw:
            uf = (numpy.array(aniso_raw[a.name]) * sc6).tolist()
            u6 = U_aniso_frac(*uf).u_aniso_orth(cell).as_numpy()  # [u11,u22,u33,u12,u13,u23]
            aniso_atoms.append(a)
            aniso_u6.append(u6)
            bfactors[i] = b_from_u * (u6[0] + u6[1] + u6[2]) / 3.0   # display B_eq
        else:
            u = None
            if uiso not in ('', '?', '.'):
                try:
                    u = float(_strip_su(uiso))
                except ValueError:
                    u = None
            bfactors[i] = b_from_u * (u if u is not None else 0.05)
    atoms.bfactors = bfactors
    if aniso_atoms:
        from chimerax.atomic import Atoms
        Atoms(aniso_atoms).aniso_u6 = numpy.array(aniso_u6, numpy.float32)


def show_cod_crystal(session, path, hkl_path=None, radiation='auto', fragments='rename'):
    '''
    Open a small-molecule (COD) CIF as a live crystal structure: the model in its
    unit cell (in Clipper's coordinate frame, so the density aligns and the corecif
    oblique-cell coordinate error is corrected), crystallographic symmetry, and -
    when reflections are available - a live 2Fo-Fc / Fo-Fc electron-density map
    computed by direct FFT with no bulk solvent, updating as the model changes.
    Returns the crystal SymmetryManager.

    radiation: 'xray', 'electron' (micro-ED / 3D-ED), or 'auto' (default) to read
    the experiment type from the CIF (_diffrn_radiation_probe), defaulting to X-ray.

    fragments: 'off' (leave the ASU as one UNL residue), 'rename' (default; split into
    named covalent-fragment residues) or 'complete' (also add symmetry-generated atoms
    to finish molecules split across a special position). See io.fragments.
    '''
    import os
    from ..symmetry import get_map_mgr
    radiation = _resolve_radiation(radiation, path)
    model = open_small_molecule_cif(session, path)
    cell, spacegroup, grid = crystal_symmetry_from_cif_file(path)
    # Coordinates are already in Clipper's frame (corrected at open); the corecif
    # oblique-cell distortion the old inline fix here handled is now internal.
    # Put the crystallographic per-atom data corecif omits (isotropic B, orthogonal
    # aniso U, ionic species) onto the model, so the live map reads it from the model.
    hydrate_small_molecule_model(session, model, path, cell, radiation)
    # Split into named fragments before the live map is set up (splitting reorders, and
    # 'complete' extends, model.atoms - and the live map reads atoms in model order).
    if fragments and fragments != 'off':
        from .fragments import split_fragments
        split_fragments(session, model, cell, spacegroup, grid, mode=fragments,
                         path=path, log=session.logger)
    mmgr = get_map_mgr(model, create=True)
    smd = _small_molecule_map_data(model, path, hkl_path, cell, spacegroup, grid, radiation)
    if smd is not None:
        from ..maps.xmapset import XmapSet
        xset = XmapSet(mmgr, small_molecule_data=smd)
        session.logger.info('(CLIPPER) %s: live maps using %s scattering factors.'
            % (os.path.basename(path), radiation))
        # Surface the model-vs-data R on the command line too (the macromolecular
        # 'clipper open' already logs Rwork/Rfree). This is the single R over all
        # reflections the live map computes through the shared R-factor core, so the
        # command log, the GUI status line and the headless metrics all agree.
        rwork = getattr(xset, 'rwork', None)
        if rwork:
            session.logger.info('(CLIPPER) %s: R = %.4f (all reflections).'
                % (os.path.basename(path), rwork))
    else:
        session.logger.info('(CLIPPER) No reflections found for %s; showing model '
            '+ symmetry only.' % os.path.basename(path))
    crystal_mgr = mmgr.crystal_mgr
    if crystal_mgr.id is None:
        session.models.add([crystal_mgr])
    return crystal_mgr


def _small_molecule_map_data(model, path, hkl_path, cell, spacegroup, grid, radiation='xray'):
    '''Assemble the small_molecule_data dict consumed by XmapSet: crystal definition,
    Fobs amplitudes aligned to a fresh HKL_info, and the live model. The live map builds
    its per-atom structure-factor inputs directly from the model each recompute (see
    SmallMoleculeXmapMgr), so no frozen scaffold or atom-index map is needed. Returns
    None if no reflections are available.'''
    import numpy
    from .. import HKL_info, HKL_data_F_sigF
    from ..clipper_python import HKL, Resolution

    refln = None
    if hkl_path:
        from chimerax.mmcif import get_cif_tables
        tbls = get_cif_tables(hkl_path, ['refln'])
        refln = tbls[0] if tbls else None
        if refln is not None and not refln.has_field('F_squared_meas'):
            refln = None
    if refln is None:
        refln = _find_reflection_table(path)
    if refln is None:
        return None

    rows = refln.fields(('index_h', 'index_k', 'index_l', 'F_squared_meas'),
                        allow_missing_fields=True)
    h = numpy.array([[int(r[0]), int(r[1]), int(r[2])] for r in rows], numpy.int32)
    fsq = numpy.array([float(_strip_su(r[3])) for r in rows], numpy.double)
    # Too few measured reflections to support the aniso+spline scale the live map runs
    # each recompute (see _MIN_SCALING_REFLECTIONS): treat as "no usable reflections" so
    # no live map is built, rather than let SmallMoleculeXmapMgr feed a singular fit to
    # Clipper (a hard crash on Linux). Also guards the max() below against an empty list.
    if int((fsq > 0).sum()) < _MIN_SCALING_REFLECTIONS:
        return None
    res = 1.0 / numpy.sqrt(max(HKL(h[i].tolist()).invresolsq(cell) for i in range(len(h))))
    hklinfo = HKL_info(spacegroup, cell, Resolution(res), True)
    fobs = HKL_data_F_sigF(hklinfo)
    fo = numpy.sqrt(numpy.clip(fsq, 0, None))
    fobs.set_data(h, numpy.stack([fo, numpy.ones_like(fo)], axis=1))

    return {'cell': cell, 'spacegroup': spacegroup, 'grid': grid, 'hklinfo': hklinfo,
            'resolution': res, 'fobs': fobs, 'structure': model, 'path': path,
            'radiation': radiation}
