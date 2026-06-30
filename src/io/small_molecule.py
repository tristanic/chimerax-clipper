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


def open_small_molecule_cif(session, path, file_name=None):
    '''Open a small-molecule CIF via ChimeraX's corecif parser. Returns the
    (single) AtomicStructure model, added to the session.'''
    from chimerax.mmcif.corecif import open_corecif
    models, status = open_corecif(session, path, file_name=file_name)
    if not models:
        from chimerax.core.errors import UserError
        raise UserError('No small-molecule structure found in %r' % path)
    session.models.add(models)
    return models[0]


def extract_cod_structure(session, path, file_name=None, keep_model=False):
    '''
    Open a small-molecule CIF and return the Mode-1 payload as a dict.

    Args:
        path: path to a (core) CIF file.
        keep_model: if False (default), the opened model is closed before
            returning (headless batch use). Set True to leave it in the session.
    '''
    model = open_small_molecule_cif(session, path, file_name=file_name)
    try:
        payload = _build_payload(session, model, path)
    finally:
        if not keep_model:
            session.models.close([model])
    return payload


def _build_payload(session, model, path):
    from chimerax.mmcif import get_cif_tables
    refine_t, diffrn_t, reflns_t, cell_t = get_cif_tables(
        path, ['refine', 'diffrn', 'reflns', 'cell'])

    cell, spacegroup, grid = crystal_symmetry_from_cif_file(path)

    atoms = model.atoms
    coords = _clipper_frame_coords(model, path, cell)
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
        # Structure-factor-derived fields: populated only when reflections are
        # available (see _structure_factor_metrics). None otherwise.
        'recomputed_r_factor': None,
        'has_structure_factors': False,
    }

    sf = _structure_factor_metrics(session, model, path, cell, spacegroup, grid)
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
    from math import sin, radians
    v = _first_cif_field(reflns_t, 'd_resolution_high')
    if v is not None:
        try:
            return float(_strip_su(v))
        except ValueError:
            pass
    wl = _first_cif_field(diffrn_t, 'radiation_wavelength')
    theta = _first_cif_field(diffrn_t, 'reflns_theta_max')
    if theta is None:
        theta = _first_cif_field(cell_t, 'measurement_theta_max')
    try:
        if wl is not None and theta is not None:
            return float(_strip_su(wl)) / (2.0 * sin(radians(float(_strip_su(theta)))))
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


def _structure_factor_metrics(session, model, path, cell, spacegroup, grid):
    '''
    If experimental structure factors are available (sibling .hkl reflection CIF,
    or an embedded _refln_ loop), compute a bulk-solvent-free R-factor of the
    model against the observed amplitudes and return the SF-derived fields.
    Returns None when no reflections are available (the common COD case).

    Small-molecule convention: Fo = sqrt(F_squared_meas) over the observed set
    (I > 2 sigma); R = sum||Fo| - k|Fc|| / sum|Fo|. (Not French-Wilson, which is a
    macromolecular technique and would not match the published small-molecule R.)
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

    from .. import HKL_info, HKL_data_I_sigI
    from ..clipper_python import HKL, Resolution, SFcalc_aniso_sum_double
    from ..clipper_python.data64 import HKL_data_F_phi_double

    # Resolution limit from the measured reflections, so the generated ASU covers
    # every observation.
    invresolsq = numpy.array([HKL(hkls[i].tolist()).invresolsq(cell)
                              for i in range(len(hkls))])
    res = 1.0 / numpy.sqrt(invresolsq.max())
    hkl_info = HKL_info(spacegroup, cell, Resolution(res), True)

    fobs = HKL_data_I_sigI(hkl_info)
    fobs.set_data(hkls, numpy.stack([fsq, sig], axis=1))

    atoms = _sfcalc_atom_list(path, cell, spacegroup, grid)
    if atoms is None:
        return {'has_structure_factors': True, 'recomputed_r_factor': None,
                'n_reflections_used': 0}

    # Exact direct summation of structure factors with NO bulk-solvent term
    # (correct for densely-packed small-molecule crystals; FFT is an unnecessary
    # approximation at this atom count). The functor must be built then called -
    # the shorthand "compute in constructor" form does not fill fcalc through the
    # bindings.
    fcalc = HKL_data_F_phi_double(hkl_info)
    sfcalc = SFcalc_aniso_sum_double()
    sfcalc(fcalc, atoms)

    fo, fc, fo_sig = [], [], []
    ih = fobs.first_data
    while not ih.last():
        od = fobs[ih]
        if not od.missing and od.i > 0:
            cd = fcalc[ih.hkl]
            if not cd.missing:
                fo.append(numpy.sqrt(od.i)); fc.append(cd.f); fo_sig.append(od.sigi)
        fobs.next_data(ih)
    fo = numpy.array(fo); fc = numpy.array(fc); fo_sig = numpy.array(fo_sig)
    if not len(fo) or fc.sum() == 0:
        return {'has_structure_factors': True, 'recomputed_r_factor': None,
                'n_reflections_used': 0}

    # Observed set: I > 2 sigma(I)  (fo_sig holds sigma(I), fo**2 = I).
    obs = fo ** 2 > 2.0 * fo_sig
    if not obs.any():
        obs = numpy.ones(len(fo), bool)
    k = numpy.sum(fo[obs] * fc[obs]) / numpy.sum(fc[obs] ** 2)
    r_obs = numpy.sum(numpy.abs(fo[obs] - k * fc[obs])) / numpy.sum(fo[obs])
    r_all = numpy.sum(numpy.abs(fo - k * fc)) / numpy.sum(fo)
    return {
        'has_structure_factors': True,
        'recomputed_r_factor': float(r_obs),
        'recomputed_r_factor_all': float(r_all),
        'n_reflections_used': int(obs.sum()),
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


def _sfcalc_atom_list(path, cell, spacegroup, grid):
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
    scaffold = sfcalc_scaffold(path, cell, spacegroup, grid)
    if scaffold is None:
        return None
    return atom_list_from_scaffold(scaffold)


def sfcalc_scaffold(path, cell, spacegroup, grid):
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
    for i, r in enumerate(arows):
        labels.append(r[0])
        elements.append(_element_from_type_symbol(r[1]))
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


def show_cod_crystal(session, path, hkl_path=None):
    '''
    Open a small-molecule (COD) CIF as a live crystal structure: the model in its
    unit cell (in Clipper's coordinate frame, so the density aligns and the corecif
    oblique-cell coordinate error is corrected), crystallographic symmetry, and -
    when reflections are available - a live 2Fo-Fc / Fo-Fc electron-density map
    computed by direct FFT with no bulk solvent, updating as the model changes.
    Returns the crystal SymmetryManager.
    '''
    import os
    from ..symmetry import get_map_mgr
    model = open_small_molecule_cif(session, path)
    cell, spacegroup, grid = crystal_symmetry_from_cif_file(path)
    model.atoms.coords = _clipper_frame_coords(model, path, cell)
    mmgr = get_map_mgr(model, create=True)
    smd = _small_molecule_map_data(model, path, hkl_path, cell, spacegroup, grid)
    if smd is not None:
        from ..maps.xmapset import XmapSet
        XmapSet(mmgr, small_molecule_data=smd)
    else:
        session.logger.info('(CLIPPER) No reflections found for %s; showing model '
            '+ symmetry only.' % os.path.basename(path))
    crystal_mgr = mmgr.crystal_mgr
    if crystal_mgr.id is None:
        session.models.add([crystal_mgr])
    return crystal_mgr


def _small_molecule_map_data(model, path, hkl_path, cell, spacegroup, grid):
    '''Assemble the small_molecule_data dict consumed by XmapSet: crystal definition,
    Fobs amplitudes aligned to a fresh HKL_info, the SF-calc scaffold, and a
    scaffold->model atom-index map (by label) for live coordinate updates. Returns
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
    res = 1.0 / numpy.sqrt(max(HKL(h[i].tolist()).invresolsq(cell) for i in range(len(h))))
    hklinfo = HKL_info(spacegroup, cell, Resolution(res), True)
    fobs = HKL_data_F_sigF(hklinfo)
    fo = numpy.sqrt(numpy.clip(fsq, 0, None))
    fobs.set_data(h, numpy.stack([fo, numpy.ones_like(fo)], axis=1))
    _, fo_d = fobs.data
    fo_aligned = fo_d[:, 0]

    scaffold = sfcalc_scaffold(path, cell, spacegroup, grid)
    name_to_idx = {n: i for i, n in enumerate(model.atoms.names)}
    s2m = numpy.array([name_to_idx.get(lbl, -1) for lbl in scaffold['labels']], int)
    return {'cell': cell, 'spacegroup': spacegroup, 'grid': grid, 'hklinfo': hklinfo,
            'resolution': res, 'scaffold': scaffold, 'fobs': fo_aligned,
            'structure': model, 'scaffold_to_model': s2m, 'path': path}
