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

'''
Split a small-molecule crystallographic asymmetric unit - which ChimeraX's corecif
parser delivers as a single "UNL" residue - into its individual covalent fragments,
each in its own Residue with a sensible name: CCD codes for the common simple species
(water, monatomic ions, small inorganic oxyanions) and LIG01, LIG02, ... otherwise.

Two crystallographic subtleties are handled:

  * A molecule can span a special position and be only partially present in the ASU,
    the rest supplied by symmetry (e.g. a water whose O sits on a 2-fold with only one
    of its H modelled looks like hydroxide in the ASU). Such fragments are identified
    - and, in mode 'complete', finished - via two complementary sources: Clipper's
    site multiplicity (for atoms sitting exactly ON a special position, e.g. the water
    O, whose H are usually absent from the CIF geometry tables) and the CIF's own
    _geom_bond_site_symmetry_2 records (which name every bond to a symmetry-equivalent
    atom, covering molecules split across a symmetry element that passes through a bond
    rather than an atom - e.g. a dimer across an inversion centre).
  * A molecule whose ASU pieces bond to each other only through a symmetry operator is
    reassembled into one fragment, again from the _geom_bond symmetry records.

Species naming uses a small CURATED STATIC TABLE only - this module imports neither
RDKit nor the chemsearch bundle at runtime (chemsearch is used only at development time
to build/verify the table). Metal atoms are always their own fragment (bonds to a metal
are never traversed), so coordinated water / ions / ligands separate out.
'''

# CCD 3-letter codes for monatomic ions, keyed by ChimeraX element name. A fragment
# qualifies only if it is a lone, bond-free heavy atom (so covalently bonded halogens
# are not mistaken for halide ions).
_MONATOMIC_ION = {
    'Li': 'LI', 'Na': 'NA', 'K': 'K', 'Rb': 'RB', 'Cs': 'CS',
    'Mg': 'MG', 'Ca': 'CA', 'Sr': 'SR', 'Ba': 'BA',
    'Mn': 'MN', 'Fe': 'FE', 'Co': 'CO', 'Ni': 'NI', 'Cu': 'CU', 'Zn': 'ZN',
    'Cd': 'CD', 'Hg': 'HG', 'Cr': 'CR', 'Al': 'AL',
    'F': 'F', 'Cl': 'CL', 'Br': 'BR', 'I': 'IOD',
}

# CCD codes for small inorganic species, keyed by their heavy-atom signature (the
# sorted (element, count) multiset, hydrogens ignored). Water / hydroxide / ammonium
# are single-heavy and handled separately (see _classify).
_INORGANIC = {
    (('O', 4), ('S', 1)): 'SO4',
    (('O', 4), ('P', 1)): 'PO4',
    (('N', 1), ('O', 3)): 'NO3',
    (('C', 1), ('O', 3)): 'CO3',
    (('Cl', 1), ('O', 4)): 'CLO4',
    (('B', 1), ('F', 4)): 'BF4',
    (('F', 6), ('P', 1)): 'PF6',
    (('N', 1), ('O', 2)): 'NO2',
    (('O', 3), ('S', 1)): 'SO3',
    (('Cl', 1), ('O', 3)): 'CLO3',
}

# Fragment categories, ordered for residue numbering: organics/unknown first, then
# named polyatomic species, then monatomic ions, then waters.
_CAT_ORGANIC = 0
_CAT_INORGANIC = 1
_CAT_ION = 2
_CAT_WATER = 3

# Position tolerance for treating a symmetry image as coincident with an existing atom.
_POS_TOL = 0.2

# A cross-symmetry _geom_bond entry longer than Element.bond_length(e1,e2) + this (A) is
# rejected as non-covalent (matches io.small_molecule.drop_implausibly_long_bonds so the
# cross-symmetry-reject and intra-ASU-drop thresholds agree).
_LONG_BOND_TOLERANCE = 0.5


def split_fragments(session, model, cell, spacegroup, grid, mode='rename', path=None,
                    log=None):
    '''
    Split model's asymmetric unit into named covalent-fragment residues.

    mode 'off'      -> no-op, returns None (status quo).
    mode 'rename'   -> split into fragments and name them; symmetry-aware identification
                       of special-position-split molecules; the atom set is UNCHANGED.
    mode 'complete' -> additionally add the symmetry-generated atoms that complete
                       special-position-split molecules, with corrected occupancies.

    path is the source CIF; when given, its _geom_bond_site_symmetry_2 records are used
    to reassemble molecules split across a symmetry operator (and, in 'complete' mode,
    to generate the missing symmetry copies). Everything runs headless.

    Must run BEFORE any live-map scaffold->model index map is built: it reorders (and,
    in 'complete' mode, extends) model.atoms. Returns a summary dict, or None for 'off'.
    '''
    if not mode or mode == 'off':
        return None
    if log is None:
        log = session.logger
    atoms = list(model.atoms)
    if not atoms:
        return None

    import numpy
    from ..clipper_util import site_multiplicities
    fracs = _frac_coords(atoms, cell)
    mults = site_multiplicities(numpy.array([a.coord for a in atoms], float),
                                cell, spacegroup, grid)
    special = {a: (int(mults[i]) > 1) for i, a in enumerate(atoms)}
    by_name = {a.name: a for a in atoms}

    # CIF-declared operators and the bonds that cross them (both optional).
    cif_symop_strings = _cif_symop_strings(path) if path else None
    xbonds = _cross_symmetry_bonds(path) if path else None
    if xbonds:
        xbonds = _reject_implausible_xbonds(xbonds, by_name, log, path)

    comps = _raw_components(model)
    comps = _merge_by_symmetry(comps, xbonds, by_name)

    # Classify, then order by category (organics, named species, ions, waters) with a
    # stable secondary key so LIG numbering follows the deposited atom order.
    classified = []
    for comp in comps:
        name, cat = _classify(comp, special)
        classified.append([comp, name, cat])
    classified.sort(key=lambda t: (t[2], min(a.serial_number for a in t[0])))

    # Reuse the original (UNL) residue for the first fragment - renamed in place - and
    # create new residues for the rest. This avoids leaving an empty residue behind:
    # ChimeraX's delete_residue only removes a residue by deleting its atoms, so it is a
    # no-op on an already-empty one.
    orig = list(model.residues)
    chain_id = orig[0].chain_id if orig else 'A'
    reuse = orig[0] if orig else None
    lig_n = 0
    made = []
    counts = {}
    for seq, (comp, name, cat) in enumerate(classified, start=1):
        if name is None:
            lig_n += 1
            name = 'LIG%02d' % lig_n
        if reuse is not None:
            res, reuse = reuse, None
            res.name = name
            res.number = seq
        else:
            res = model.new_residue(name, chain_id, seq)
        for a in comp:
            if a.residue is not res:
                a.residue.remove_atom(a)
                res.add_atom(a)
        made.append((res, comp))
        key = 'LIG' if name.startswith('LIG') else name
        counts[key] = counts.get(key, 0) + 1

    # Any other original residues emptied by the moves (unusual for corecif, which
    # delivers a single residue) - best-effort cleanup.
    for r in orig:
        if r is not None and not r.deleted and r.num_atoms == 0:
            model.delete_residue(r)

    n_added = 0
    if mode == 'complete':
        from .small_molecule import register_clipper_atom_attributes
        register_clipper_atom_attributes(session)
        n_added = _complete(model, made, special, fracs, cell,
                            cif_symop_strings, xbonds, by_name)
        # Completion wires bonds across the periodic boundary (a lattice-translation
        # _geom_bond partner, or a symmetry image that landed in an adjacent cell), which
        # would otherwise be drawn as a cell-spanning bond. Unwrap every molecule so each
        # bond is its minimum image; only integer lattice translations are applied, so all
        # structure factors - and hence maps, R factors and SF gradients - are invariant.
        _gather_molecules(model, cell)
        # corecif reads metal coordination from the CIF _geom_bond loop at open but ignores
        # the symmetry-coded entries, so a metal on a special position loses coordination to
        # the ligands supplied by symmetry (a metalloporphyrin on an inversion centre keeps
        # only half its metal-N bonds). Completion has now generated those ligand images, so
        # add the coordination pseudobonds the symmetry-coded _geom_bond entries name.
        _complete_metal_coordination(model, cell, cif_symop_strings, xbonds, by_name, fracs)
        # Gathering a molecule can translate a metal-coordinated atom to a lattice image
        # farther from its metal than the nearest one, leaving the coordination pseudobond
        # drawn across the cell; drop those (see _prune_wrapped_coordination).
        _prune_wrapped_coordination(model, cell)

    summary = {'mode': mode, 'n_fragments': len(made), 'n_added_atoms': n_added,
               'counts': counts}
    tally = ', '.join('%d×%s' % (n, k) for k, n in sorted(counts.items()))
    msg = '(CLIPPER) Split asymmetric unit into %d fragment(s): %s' % (len(made), tally)
    if mode == 'complete' and n_added:
        msg += '; added %d symmetry-generated atom(s) to complete special-position ' \
               'molecules (these are excluded from the structure-factor calculation ' \
               'and do not follow live edits of their source atoms)' % n_added
    log.info(msg)
    return summary


def _frac_coords(atoms, cell):
    '''Clipper fractional coordinate (Coord_frac) for each atom, keyed by atom.'''
    from ..clipper_python import Coord_orth
    out = {}
    for a in atoms:
        x, y, z = a.coord
        out[a] = Coord_orth(float(x), float(y), float(z)).coord_frac(cell)
    return out


def _raw_components(model):
    '''Connected components of the covalent-bond graph. A metal atom is always its own
    fragment: bonds to a metal are treated as fragment boundaries (never traversed), so
    coordinated water / ions / ligands separate out - matching ChimeraX's treatment of
    metal coordination as non-covalent.'''
    visited = set()
    comps = []
    for start in model.atoms:
        if start in visited:
            continue
        visited.add(start)
        comp = [start]
        if start.element.is_metal:
            comps.append(comp)
            continue
        stack = [start]
        while stack:
            a = stack.pop()
            for nb in a.neighbors:
                if nb in visited or nb.element.is_metal:
                    continue
                visited.add(nb)
                comp.append(nb)
                stack.append(nb)
        comps.append(comp)
    return comps


def _merge_by_symmetry(comps, xbonds, by_name):
    '''Merge components joined to each other only through a symmetry operator (e.g. two
    ASU pieces of one molecule), using the CIF's _geom_bond symmetry records. Bonds to a
    metal are ignored (a metal stays its own fragment).'''
    if not xbonds or len(comps) < 2:
        return comps
    comp_of = {}
    for i, c in enumerate(comps):
        for a in c:
            comp_of[a] = i
    parent = list(range(len(comps)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for l1, l2, _code in xbonds:
        a1, a2 = by_name.get(l1), by_name.get(l2)
        if a1 is None or a2 is None or a1.element.is_metal or a2.element.is_metal:
            continue
        i, j = comp_of.get(a1), comp_of.get(a2)
        if i is None or j is None or find(i) == find(j):
            continue
        parent[find(i)] = find(j)
    groups = {}
    for idx in range(len(comps)):
        groups.setdefault(find(idx), []).extend(comps[idx])
    return list(groups.values())


def _classify(comp, special):
    '''Return (residue_name, category). residue_name is None for the LIG fallback (the
    caller assigns LIG01, LIG02, ...).'''
    heavy = [a for a in comp if a.element.number > 1]
    n_h = sum(1 for a in comp if a.element.number == 1)
    if len(heavy) == 1:
        a = heavy[0]
        el = a.element.name
        if el == 'O':
            if n_h in (0, 2):
                return 'HOH', _CAT_WATER
            if n_h == 1:
                # One H on a 2-fold/mirror etc. is really water (symmetry supplies the
                # other H); on a general position it is named as-modelled (hydroxide),
                # leaving deeper is-it-really-water analysis to downstream tools.
                return ('HOH' if special.get(a) else 'OH'), _CAT_WATER
            return None, _CAT_ORGANIC
        if el == 'N' and n_h == 4:
            return 'NH4', _CAT_INORGANIC
        if n_h == 0 and el in _MONATOMIC_ION:
            return _MONATOMIC_ION[el], _CAT_ION
        return None, _CAT_ORGANIC
    name = _INORGANIC.get(_heavy_signature(heavy))
    if name is not None:
        return name, _CAT_INORGANIC
    return None, _CAT_ORGANIC


def _heavy_signature(heavy):
    from collections import Counter
    return tuple(sorted(Counter(a.element.name for a in heavy).items()))


# ---- symmetry completion (mode 'complete') ----

def _complete(model, made, special, fracs, cell, cif_symop_strings, xbonds, by_name):
    '''Add the symmetry-generated atoms that finish molecules split across a special
    position, and enforce unit occupancy on atoms lying ON a special position (PDB
    convention). Returns the number of atoms added.

    Completion is a bond-driven breadth-first search over the crystallographic bond graph
    (see _complete_fragment). Each atom is generated only as the endpoint of a bond being
    traversed - an in-model neighbour at the same rigid placement, a special-position
    atom's site-symmetry neighbour image, or a cross-symmetry _geom_bond partner at the
    COMPOSED placement P.op - so every generated atom is bonded by construction and no
    bondless atom is ever imaged. Traversing bonds (rather than applying whole-fragment
    operators once) is what lets completion (a) compose operators, reaching an atom that
    needs two symmetry steps - e.g. the fourth oxalate oxygen S5.S2(O6) in cod_2209193,
    which a single-operator pass cannot reach - and (b) image each bond with its OWN
    operator+lattice code, since one fragment's bonds can carry different codes (cod_2103691:
    the naphthalene ring closes via op2_656 while its sulfonate's third oxygen needs
    op2_657). Coincidence is tested by MINIMUM IMAGE, which both wires a deposited orphan
    whose real partner is a symmetry image and bounds the walk (images a lattice-vector
    apart fold together, so the orbit of a finite molecule is finite).'''
    for a, on_special in special.items():
        if on_special:
            a.occupancy = 1.0
    if not cif_symop_strings:
        return 0
    # Each operator as (rotation 3x3, translation 3) in fractional space; applying it is
    # then a plain rot.uvw + trn (no Clipper Symop needed - its constructor takes no
    # xyz string). CIF order is preserved, so index n-1 matches a _geom_bond code n.
    import numpy
    ops_rt = [_parse_xyz_op(s) for s in cif_symop_strings]

    # Cross-symmetry _geom_bond edges as (atom_1, atom_2, op_index, lattice): a bond
    # between atom_1 (identity placement) and op.atom_2. Metal endpoints are skipped (a
    # metal stays its own fragment); hydrogen endpoints are skipped because terminal H are
    # relocated by reassemble_symmetry_scattered_hydrogens (run before this) and water H
    # come through the site-symmetry branch, so wiring a symmetry-coded H bond here would
    # double it.
    edges = []
    if xbonds:
        for l1, l2, (n, dk, dl, dm) in xbonds:
            a1, a2 = by_name.get(l1), by_name.get(l2)
            if a1 is None or a2 is None:
                continue
            if a1.element.is_metal or a2.element.is_metal:
                continue
            if a1.element.number == 1 or a2.element.number == 1:
                continue
            if not (1 <= n <= len(ops_rt)):
                continue
            edges.append((a1, a2, n - 1, numpy.array([dk, dl, dm], float)))

    existing_names = set(model.atoms.names)
    n_added = 0
    for res, comp in made:
        n_added += _complete_fragment(res, comp, special, fracs, cell, ops_rt,
                                       edges, existing_names)
    return n_added


def _gather_molecules(model, cell):
    '''Unwrap every covalently-bonded molecule so each bond is its minimum image - i.e.
    make molecules whole across periodic boundaries. Completion wires bonds a molecule can
    have across the unit-cell boundary: a pure lattice-translation _geom_bond partner (a
    carboxylate whose second O is one cell over, cod_2223203/2223282) or a symmetry image
    that landed in an adjacent cell (a special-position molecule's generated half,
    cod_2103691); left alone these draw as cell-spanning bonds. A breadth-first walk of the
    covalent-bond graph translates each atom by the integer lattice vector that places it on
    the minimum image of its already-placed neighbour. ONLY integer lattice translations are
    applied, so every structure factor is exactly invariant (exp(2*pi*i h.(r+n)) ==
    exp(2*pi*i h.r)) - maps, R factors and SF gradients are untouched; this is purely a
    spatial gathering. Metals sit in their own component (their bonds were demoted to
    coordination pseudobonds, which are not covalent bonds), so they are left in place.'''
    import numpy
    from ..clipper_python import Coord_orth, Coord_frac
    atoms = model.atoms
    if len(atoms) == 0:
        return
    frac = {}
    for a in atoms:
        x, y, z = a.coord
        frac[a] = numpy.array(Coord_orth(float(x), float(y), float(z)).coord_frac(cell).uvw,
                              float)
    # Seed from the DEPOSITED atoms before the symmetry-generated images, so a real atom
    # anchors each molecule and only the generated copies are translated onto it. This
    # keeps deposited atoms at their deposited positions wherever a molecule has one to
    # anchor on - notably a metal-coordinated oxygen whose completed sulfonate image would
    # otherwise drag it (and its coordination bond) a cell away (cod_2103691).
    visited = set()
    for seed in sorted(atoms, key=lambda a: bool(getattr(a, 'clipper_sf_exclude', False))):
        if seed in visited:
            continue
        visited.add(seed)
        stack = [seed]
        while stack:
            a = stack.pop()
            fa = frac[a]
            for nb in a.neighbors:
                if nb in visited:
                    continue
                shift = numpy.round(frac[nb] - fa)
                if shift.any():
                    nf = frac[nb] - shift
                    frac[nb] = nf
                    co = Coord_frac(float(nf[0]), float(nf[1]), float(nf[2])).coord_orth(cell)
                    nb.coord = numpy.array([co.x, co.y, co.z])
                visited.add(nb)
                stack.append(nb)


def _complete_metal_coordination(model, cell, cif_symop_strings, xbonds, by_name, fracs):
    '''Add the metal-coordination pseudobonds to symmetry-generated ligand images that
    completion produced. corecif builds a metal's intra-ASU (identity) coordination from the
    CIF _geom_bond loop at open, but - like the symmetry-coded COVALENT bonds - ignores the
    symmetry-coded coordination entries, so a metal on a special position keeps coordination
    only to the ligands modelled in the ASU. The canonical case is a metalloporphyrin on an
    inversion centre (cod_1561255): _geom_bond lists 3 deposited Mg-N and 3 symmetry-coded
    ones, so the Mg comes out 3-coordinate when the completed macrocycle is 6-coordinate -
    fatal for the only sensible classical-MD model, which promotes the coordination to
    covalent. Each metal _geom_bond entry with a symmetry code names exactly one such bond -
    metal(identity) coordinating op(ligand) - and completion has already built that ligand
    image as covalent framework, so add the pseudobond to it (matched by MINIMUM-IMAGE
    position, since unwrapping may have translated the image by a lattice vector). This uses
    the deposited coordination partner + operator, so no distance criterion is guessed. A
    coordination image that lands across the cell is left to _prune_wrapped_coordination.'''
    if not cif_symop_strings or not xbonds:
        return
    import numpy
    from ..clipper_python import Coord_orth, Coord_frac
    ops_rt = [_parse_xyz_op(s) for s in cif_symop_strings]
    # (metal_atom, ligand_atom, rot, trn): the metal coordinates rot.ligand_frac + trn.
    metal_edges = []
    for l1, l2, (n, dk, dl, dm) in xbonds:
        a1, a2 = by_name.get(l1), by_name.get(l2)
        if a1 is None or a2 is None or not (1 <= n <= len(ops_rt)):
            continue
        if a1.element.is_metal == a2.element.is_metal:  # non-metal, or metal-metal
            continue
        if a1.element.number == 1 or a2.element.number == 1:
            continue
        rot, trn = ops_rt[n - 1]
        trn = trn + numpy.array([dk, dl, dm], float)
        if a1.element.is_metal:                         # a1 coordinates op(a2)
            metal_edges.append((a1, a2, rot, trn))
        else:                                           # a1(id)-op(a2) == a2(id)-op^-1(a1)
            irot, itrn = _invert_op(rot, trn)
            metal_edges.append((a2, a1, irot, itrn))
    if not metal_edges:
        return
    atoms = list(model.atoms)
    fr = {}
    for a in atoms:
        x, y, z = a.coord
        fr[a] = numpy.array(Coord_orth(float(x), float(y), float(z)).coord_frac(cell).uvw,
                            float)

    def find_at(target, elem):
        for a in atoms:
            if a.element.number != elem:
                continue
            d = fr[a] - target
            d = d - numpy.round(d)                      # minimum image
            co = Coord_frac(float(d[0]), float(d[1]), float(d[2])).coord_orth(cell)
            if (co.x * co.x + co.y * co.y + co.z * co.z) ** 0.5 < _POS_TOL:
                return a
        return None

    pbg = model.pseudobond_group('metal coordination bonds', create_type='per coordset')
    existing = set(frozenset(pb.atoms) for pb in pbg.pseudobonds)
    for metal, lig, rot, trn in metal_edges:
        img = find_at(rot.dot(numpy.array(fracs[lig].uvw)) + trn, lig.element.number)
        if img is None or img is metal:
            continue
        key = frozenset((metal, img))
        if key in existing:
            continue
        pbg.new_pseudobond(metal, img)
        existing.add(key)


def _prune_wrapped_coordination(model, cell):
    '''Delete metal-coordination pseudobonds left pointing across the unit cell by molecule
    unwrapping. _gather_molecules translates a molecule so its covalent bonds are whole; a
    metal-coordinated atom can thereby move to a lattice image farther from its metal than
    the nearest one, so the coordination pseudobond is drawn cell-spanning and coordinates a
    symmetry copy that is not a distinct atom in the gathered model. Drop any pseudobond
    whose current length exceeds its minimum-image distance (a nearer periodic image of the
    ligand exists) - so a genuine nearest-image coordination (e.g. a Ca...P at its deposited
    3.49 A) is kept, while a bond stretched to a wrapped copy is removed. This is unavoidable
    for a coordination polymer where a ligand bridges metals in adjacent cells (a sulfonyl O
    bound to one metal and, a cell over, coordinating another - cod_2223282): the covalent
    molecule is kept whole and the redundant cross-cell coordination display is dropped.'''
    import numpy
    from ..clipper_python import Coord_orth, Coord_frac
    def frac(a):
        x, y, z = a.coord
        return numpy.array(Coord_orth(float(x), float(y), float(z)).coord_frac(cell).uvw,
                           float)
    for pbg in list(model.pbg_map.values()):
        for pb in list(pbg.pseudobonds):
            a1, a2 = pb.atoms
            d = frac(a2) - frac(a1)
            dmin = d - numpy.round(d)
            co = Coord_frac(float(dmin[0]), float(dmin[1]), float(dmin[2])).coord_orth(cell)
            len_min = (co.x * co.x + co.y * co.y + co.z * co.z) ** 0.5
            if pb.length > len_min + _POS_TOL:
                pb.delete()


def _invert_op(rot, trn):
    '''Inverse of the affine fractional operator (rot, trn): (rot^-1, -rot^-1.trn).'''
    import numpy
    ir = numpy.linalg.inv(rot)
    return ir, -ir.dot(trn)


def _complete_fragment(res, comp, special, fracs, cell, ops_rt, edges, existing_names):
    '''Complete one fragment by bond-driven BFS. Returns the number of atoms added.'''
    import numpy
    from collections import deque
    from chimerax.atomic import struct_edit
    from ..clipper_python import Coord_frac

    comp_set = set(comp)
    # Partial-occupancy atoms are disorder components of a shared/averaged site; snapshot
    # the set up front (before any promotion below) so the skip is driven by the ORIGINAL
    # occupancy. A full-atom symmetry copy of one would clash with its own disorder partner
    # and break the fragment's charge/valence - the canonical case is a 0.5-occupancy
    # bridging proton on a 2-fold (cod_2215867): keep the single deposited atom, and
    # promote it to full occupancy so the completed molecule is a physical snapshot (one
    # whole shared proton, integer net charge) rather than two clashing half-atoms.
    disorder = {a for a in comp if a.occupancy < 0.99}
    promote = set()

    # A bridging proton - a hydrogen covalently bonded to two heavy atoms, i.e. a
    # short-strong O-H...O / N-H...O hydrogen bond the CIF _geom_bond loop lists as two bonds
    # (a carboxyl-carboxylate shared proton, cod_2017117) - fuses the two moieties it bridges
    # into one completion fragment. Imaging THROUGH it (propagating a symmetry operator from
    # one moiety across the bridge into the other) drags a spurious image of the proton and
    # the far moiety a whole cell away, wiring a ~17 A bond that no unwrap can resolve (the
    # bridged network is periodic). So the imaging BFS does not step onto a bridging proton:
    # its deposited bonds are left intact (the shared proton is kept, for downstream
    # tautomer / dual-bond force-field treatment), only its symmetry IMAGE is not generated.
    bridging = {a for a in comp if a.element.number == 1
                and sum(1 for nb in a.neighbors if nb.element.number > 1) >= 2}

    I3 = numpy.eye(3)
    Z3 = numpy.zeros(3)

    def frac_of(atom):
        return numpy.array(fracs[atom].uvw, float)

    # placed[i] = [element_number, fractional position, atom]; seeded with the ASU atoms.
    placed = [[a.element.number, frac_of(a), a] for a in comp]

    def cart_len(dfrac):
        co = Coord_frac(float(dfrac[0]), float(dfrac[1]), float(dfrac[2])).coord_orth(cell)
        return (co.x * co.x + co.y * co.y + co.z * co.z) ** 0.5

    def occupant(elem, fpos):
        for e, fp, at in placed:
            if e != elem:
                continue
            d = fpos - fp
            d = d - numpy.round(d)              # minimum image
            if cart_len(d) < _POS_TOL:
                return at
        return None

    added = []

    def get_or_create(source, gR, gt):
        '''Atom occupying source's image at placement (gR, gt), creating it if the position
        is not already occupied (minimum image). Returns (atom, is_new); (None, False) if a
        disorder atom's image would have been distinct (kept single, source promoted).'''
        fpos = gR.dot(frac_of(source)) + gt
        occ = occupant(source.element.number, fpos)
        if occ is not None:
            return occ, False
        if source in disorder:
            promote.add(source)
            return None, False
        co = Coord_frac(float(fpos[0]), float(fpos[1]), float(fpos[2])).coord_orth(cell)
        na = struct_edit.add_atom(_unique_atom_name(source.name, existing_names),
                                  source.element, res,
                                  numpy.array([co.x, co.y, co.z]),
                                  occupancy=source.occupancy, bfactor=source.bfactor,
                                  info_from=source)
        existing_names.add(na.name)
        # A crystallographic-symmetry image of an ASU atom; Clipper's SFcalc regenerates it
        # from the source, so keep it out of the structure-factor sum.
        na.clipper_sf_exclude = True
        placed.append([source.element.number, fpos, na])
        added.append(na)
        return na, True

    def link(atom, source, gR, gt, queue):
        partner, is_new = get_or_create(source, gR, gt)
        if partner is None:
            return
        if partner is not atom and partner not in atom.neighbors:
            struct_edit.add_bond(atom, partner)
        if is_new:
            queue.append((source, gR, gt, partner))

    def site_ops(atom):
        '''Non-identity operators (rot, translation) that fix atom's position modulo
        lattice - its site symmetry. The translation is folded so the fixed point maps
        exactly onto itself, so a neighbour imaged under it lands bonded.'''
        fa = frac_of(atom)
        out = []
        for idx, (rot, trn) in enumerate(ops_rt):
            fimg = rot.dot(fa) + trn
            shift = numpy.round(fimg - fa)
            if numpy.allclose(fimg - shift, fa, atol=1e-3):
                t = trn - shift
                if idx == 0 and numpy.allclose(t, Z3, atol=1e-9):
                    continue
                out.append((rot, t))
        return out

    # Bound the walk: a finite molecule has a finite minimum-image orbit, but a bad-data
    # extended framework / polymer would not terminate cleanly, so cap generation at a
    # generous multiple of the fragment's heavy-atom count (the runaway is dropped/flagged
    # downstream).
    hard_cap = max(200, 30 * sum(1 for a in comp if a.element.number > 1))

    queue = deque((a, I3, Z3, a) for a in comp)
    while queue:
        if len(added) > hard_cap:
            break
        source, gR, gt, atom = queue.popleft()
        # (1) in-model covalent neighbours, at the same rigid placement.
        for nb in source.neighbors:
            if not nb.element.is_metal and nb in comp_set and nb not in bridging:
                link(atom, nb, gR, gt, queue)
        # (1b) a special-position atom's neighbours are also imaged by its site symmetry
        #      (the water-H case: O on a 2-fold with only one H modelled).
        if special.get(source):
            for sR, st in site_ops(source):
                pR = gR.dot(sR)
                pt = gR.dot(st) + gt
                for nb in source.neighbors:
                    if not nb.element.is_metal and nb in comp_set and nb not in bridging:
                        link(atom, nb, pR, pt, queue)
        # (2) cross-symmetry _geom_bond edges, imaged with the edge's own operator+lattice
        #     at the COMPOSED placement.
        for a1, a2, opi, latt in edges:
            rot, trn = ops_rt[opi]
            trn = trn + latt
            if source is a1:                    # partner a2 at g.op
                pR = gR.dot(rot)
                pt = gR.dot(trn) + gt
                link(atom, a2, pR, pt, queue)
            elif source is a2:                  # partner a1 at g.op^-1
                irot, itrn = _invert_op(rot, trn)
                pR = gR.dot(irot)
                pt = gR.dot(itrn) + gt
                link(atom, a1, pR, pt, queue)

    for a in promote:
        a.occupancy = 1.0
    return len(added)


def _unique_atom_name(base, existing_names):
    '''A structure-unique atom name derived from base by suffixing a letter (then a
    number). Uniqueness across the whole structure matters: the live small-molecule map
    aligns its structure-factor scaffold to the model by atom name.'''
    import string
    for suf in string.ascii_uppercase:
        cand = base + suf
        if cand not in existing_names:
            return cand
    i = 1
    while (base + str(i)) in existing_names:
        i += 1
    return base + str(i)


# ---- CIF symmetry records ----

def _cif_symop_strings(path):
    '''The space group's symmetry operators as xyz strings, in CIF order (1-based, as
    referenced by _geom_bond_site_symmetry_2 codes). None if absent.'''
    from chimerax.mmcif import get_cif_tables
    # get_cif_tables returns [] (not padded) when the file has NONE of the requested
    # categories - guard the unpack for CIFs lacking any symmetry-operator table.
    tables = get_cif_tables(path, ['space_group_symop', 'symmetry'])
    if not tables:
        return None
    sgs_t, symm_t = tables
    for table, field in ((sgs_t, 'operation_xyz'), (symm_t, 'equiv_pos_as_xyz')):
        if table is not None and table.has_field(field):
            ops = [r[0].strip().strip("'\"").replace(' ', '')
                   for r in table.fields((field,))]
            ops = [o for o in ops if o]
            if ops:
                return ops
    return None


def _cross_symmetry_bonds(path):
    '''Bonds to a symmetry-equivalent atom, from the CIF _geom_bond loop: a list of
    (label_1, label_2, (n, dk, dl, dm), distance) where n is the 1-based symmetry-operator
    index, (dk, dl, dm) the accompanying lattice translation, and distance the deposited
    _geom_bond_distance in Angstroms (float, or None if the column is absent/unparseable -
    used to reject non-covalent close contacts). None if the loop is absent.'''
    from chimerax.mmcif import get_cif_tables
    # get_cif_tables returns [] (not [empty_table]) when the file has no _geom_bond loop.
    tables = get_cif_tables(path, ['geom_bond'])
    gb = tables[0] if tables else None
    if gb is None or not gb.has_field('atom_site_label_1'):
        return None
    rows = gb.fields(('atom_site_label_1', 'atom_site_label_2', 'site_symmetry_2',
                      'distance'), allow_missing_fields=True)
    out = []
    for r in rows:
        code = r[2].strip() if len(r) > 2 else '.'
        if not code or code in ('.', '?'):
            continue
        nklm = _parse_symcode(code)
        if nklm is None:
            continue
        dist = None
        if len(r) > 3:
            ds = r[3].strip()
            if ds and ds not in ('.', '?'):
                try:                        # strip the s.u. suffix, e.g. '2.8898(16)'
                    dist = float(ds.split('(')[0])
                except ValueError:
                    dist = None
        out.append((r[0], r[1], nklm, dist))
    return out


def _reject_implausible_xbonds(xbonds, by_name, log, path):
    '''Drop cross-symmetry _geom_bond entries too long to be covalent, and return the kept
    (label_1, label_2, (n,dk,dl,dm)) triples (distance stripped, so the downstream
    consumers are unchanged).

    Depositors routinely list secondary CLOSE CONTACTS - not bonds - in _geom_bond with a
    cross-symmetry code, especially in N/O-dense energetic materials (nitro / N-oxide
    O...O and O...N contacts at ~2.8 A; e.g. cod_2204168). Symmetry completion would take
    each such entry at face value and image it, wiring a physically impossible bond and
    generating a shell of spurious symmetry copies. An entry is rejected when its deposited
    distance exceeds ``Element.bond_length(e1, e2) + _LONG_BOND_TOLERANCE`` (the same bar
    io.small_molecule.drop_implausibly_long_bonds uses). Entries with no usable distance,
    or to/from a metal, are kept (no basis to reject). One warning summarises what was
    ignored. This is the cross-symmetry counterpart of
    io.small_molecule.drop_implausibly_long_bonds (which deletes the same kind of spurious
    bond among the ASU's own atoms at open time); both use the same bond_length + tolerance
    bar - here we decline to CREATE the bond, there it deletes one corecif already made.'''
    from chimerax.atomic import Element
    import os
    kept = []
    dropped = []
    for l1, l2, nklm, dist in xbonds:
        a1, a2 = by_name.get(l1), by_name.get(l2)
        if (dist is not None and a1 is not None and a2 is not None
                and not a1.element.is_metal and not a2.element.is_metal):
            bl = Element.bond_length(a1.element, a2.element)
            if bl and dist > bl + _LONG_BOND_TOLERANCE:
                dropped.append((l1, l2, dist, bl))
                continue
        kept.append((l1, l2, nklm))
    if dropped and log is not None:
        dropped.sort(key=lambda o: o[2] - o[3], reverse=True)   # most egregious first
        ex = '; '.join('%s-%s %.2f A (max ~%.2f)' % (n1, n2, d, bl + _LONG_BOND_TOLERANCE)
                       for n1, n2, d, bl in dropped[:3])
        name = os.path.basename(path) if path else '?'
        log.warning('(CLIPPER) %s: ignored %d cross-symmetry _geom_bond entr%s too long to '
            'be covalent (deposited non-bonding close contact(s) listed as bonds, e.g. '
            'O...O / O...N in N/O-dense structures); symmetry completion will not wire '
            'them: %s%s' % (name, len(dropped), 'y' if len(dropped) == 1 else 'ies', ex,
            ' ...' if len(dropped) > 3 else ''))
    return kept


def _parse_xyz_op(s):
    '''Parse a symmetry operator string ('x,y,z', '-x+1/2,y,-z-1/2', ...) into a
    (rotation 3x3, translation 3) pair of numpy arrays in fractional space.'''
    import numpy
    rot = numpy.zeros((3, 3))
    trn = numpy.zeros(3)
    axes = {'x': 0, 'y': 1, 'z': 2}
    for i, comp in enumerate(s.lower().replace(' ', '').split(',')[:3]):
        for term in comp.replace('-', '+-').split('+'):
            if not term:
                continue
            sign = 1.0
            if term[0] == '-':
                sign, term = -1.0, term[1:]
            if not term:
                continue
            if term[-1] in axes:
                coef = term[:-1]
                rot[i, axes[term[-1]]] += sign * (1.0 if coef == '' else _frac_to_float(coef))
            else:
                trn[i] += sign * _frac_to_float(term)
    return rot, trn


def _frac_to_float(t):
    if '/' in t:
        a, b = t.split('/')
        return float(a) / float(b)
    return float(t)


def _parse_symcode(code):
    '''CIF symmetry code 'n_klm' -> (n, k-5, l-5, m-5); bare 'n' -> (n, 0, 0, 0).'''
    code = code.strip()
    if '_' in code:
        n_s, t_s = code.split('_', 1)
        if not (n_s.isdigit() and len(t_s) == 3 and t_s.isdigit()):
            return None
        return (int(n_s), int(t_s[0]) - 5, int(t_s[1]) - 5, int(t_s[2]) - 5)
    if code.isdigit():
        return (int(code), 0, 0, 0)
    return None
