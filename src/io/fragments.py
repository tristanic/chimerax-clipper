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
    convention). Returns the number of atoms added.'''
    import numpy
    for a, on_special in special.items():
        if on_special:
            a.occupancy = 1.0
    if not cif_symop_strings:
        return 0
    # Each operator as (rotation 3x3, translation 3) in fractional space; applying it is
    # then a plain rot.uvw + trn (no Clipper Symop needed - its constructor takes no
    # xyz string). CIF order is preserved, so index n-1 matches a _geom_bond code n.
    ops_rt = [_parse_xyz_op(s) for s in cif_symop_strings]

    existing_names = set(model.atoms.names)
    n_added = 0
    for res, comp in made:
        comp_set = set(comp)
        ops = []   # each: (rot, translation, [(existing_atom, source_atom), ...])
        # (a) Site-symmetry of any atom sitting ON a special position generates the
        #     images of its (off-axis) neighbours - the water-H case. Fold in the shift
        #     that fixes the special atom exactly, so the image is the bonded one.
        anchor = next((a for a in comp if special.get(a)), None)
        if anchor is not None:
            fa = numpy.array(fracs[anchor].uvw)
            for rot, trn in ops_rt:
                fimg = rot.dot(fa) + trn
                shift = numpy.round(fimg - fa)
                if numpy.allclose(fimg - shift, fa, atol=1e-3):
                    ops.append((rot, trn - shift, []))
        # (b) Bonds to a symmetry copy (CIF _geom_bond) generate that copy - covers
        #     molecules split across an element passing through a bond (e.g. a dimer on
        #     an inversion centre, where no atom is on the special position).
        if xbonds:
            for l1, l2, (n, dk, dl, dm) in xbonds:
                a1, a2 = by_name.get(l1), by_name.get(l2)
                if a1 is None or a2 is None or a1.element.is_metal or a2.element.is_metal:
                    continue
                if a1 not in comp_set and a2 not in comp_set:
                    continue
                # Terminal hydrogen whose bond is symmetry-coded (deposited in a
                # neighbouring ASU): reassemble_symmetry_scattered_hydrogens (run before
                # this, under `fragments complete`) relocates the real H onto its parent.
                # Do NOT also generate a symmetry image of it here - that duplicates the H
                # and, because reassemble put the real atom at exactly the image position,
                # over-bonds it. Special-position water H are handled by branch (a) above.
                if a2.element.number == 1 or a1.element.number == 1:
                    continue
                if not (1 <= n <= len(ops_rt)):
                    continue
                # bond is a1 -- S(a2): applying S to the fragment places S(a2) bonded to a1.
                rot, trn = ops_rt[n - 1]
                ops.append((rot, trn + numpy.array([dk, dl, dm], float), [(a1, a2)]))
        n_added += _apply_ops(res, comp, ops, fracs, cell, existing_names)
    return n_added


def _apply_ops(res, comp, ops, fracs, cell, existing_names):
    import numpy
    from chimerax.atomic import struct_edit
    from ..clipper_python import Coord_frac
    added = []   # atoms generated for this fragment (dedup across ops)

    def occupied(xyz):
        for a in comp:
            if numpy.linalg.norm(numpy.array(a.coord) - xyz) < _POS_TOL:
                return a
        for na in added:
            if numpy.linalg.norm(numpy.array(na.coord) - xyz) < _POS_TOL:
                return na
        return None

    for rot, trn, links in ops:
        img = {}   # source atom -> atom occupying its image position (existing or new)
        for a in comp:
            uvw = rot.dot(numpy.array(fracs[a].uvw)) + trn
            co = Coord_frac(*uvw).coord_orth(cell)
            xyz = numpy.array([co.x, co.y, co.z])
            img[a] = [occupied(xyz), xyz]
        new_this = set()
        for a in comp:
            target, xyz = img[a]
            if target is not None:
                continue
            name = _unique_atom_name(a.name, existing_names)
            existing_names.add(name)
            na = struct_edit.add_atom(name, a.element, res, xyz,
                                      occupancy=a.occupancy, bfactor=a.bfactor,
                                      info_from=a)
            # This is a crystallographic-symmetry image of an ASU atom; Clipper's SFcalc
            # regenerates it from the source, so keep it out of the structure-factor sum.
            na.clipper_sf_exclude = True
            img[a][0] = na
            new_this.add(na)
            added.append(na)
        # Mirror the fragment's own connectivity onto the freshly generated atoms.
        for a in comp:
            ia = img[a][0]
            if ia is None:
                continue
            for nb in a.neighbors:
                if nb not in img:
                    continue
                inb = img[nb][0]
                if inb is None or inb is ia:
                    continue
                if (ia in new_this or inb in new_this) and inb not in ia.neighbors:
                    struct_edit.add_bond(ia, inb)
        # Wire the explicit cross-symmetry bond(s): original atom -> image of its partner.
        for orig, src in links:
            isrc = img.get(src, [None])[0]
            if isrc is not None and isrc is not orig and isrc not in orig.neighbors:
                struct_edit.add_bond(orig, isrc)
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
    (label_1, label_2, (n, dk, dl, dm)) where n is the 1-based symmetry-operator index
    and (dk, dl, dm) the accompanying lattice translation. None if the loop is absent.'''
    from chimerax.mmcif import get_cif_tables
    # get_cif_tables returns [] (not [empty_table]) when the file has no _geom_bond loop.
    tables = get_cif_tables(path, ['geom_bond'])
    gb = tables[0] if tables else None
    if gb is None or not gb.has_field('atom_site_label_1'):
        return None
    rows = gb.fields(('atom_site_label_1', 'atom_site_label_2', 'site_symmetry_2'),
                     allow_missing_fields=True)
    out = []
    for r in rows:
        code = r[2].strip() if len(r) > 2 else '.'
        if not code or code in ('.', '?'):
            continue
        nklm = _parse_symcode(code)
        if nklm is not None:
            out.append((r[0], r[1], nklm))
    return out


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
