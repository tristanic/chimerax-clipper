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
GUI-free crystallographic symmetry realisation.

The symmetry-copy machinery lives here rather than inside the ``clipper
symcopies`` command (or the GUI-only :class:`AtomicSymmetryModel`) so it can be
driven headless (``--nogui``) and reused by other tools. Nothing in this module
creates graphics, Drawings or touches ``session.main_view``: it operates purely
on a ChimeraX structure plus its crystallographic core (cell, spacegroup, grid,
Unit_Cell), all of which can be built from model metadata or a structure-factor
file without a running GUI.
'''

import numpy
from . import clipper_python

DEFAULT_SYM_SEARCH_FREQUENCY = 2


def sym_select_within(structure, cell, grid, unit_cell, atoms, cutoff,
        coords=None, whole_residues=True,
        sym_search_frequency=DEFAULT_SYM_SEARCH_FREQUENCY):
    '''
    Given a set of atoms, return a (atoms, symmetry matrices, sym_indices,
    symops) tuple giving all atoms and their symmetry operators within the
    given cutoff distance from any atom in the primary set. GUI-free core of
    the identically-named :class:`AtomicSymmetryModel` method (which delegates
    here).
    Args:
        structure:
            The ChimeraX AtomicStructure whose atoms define the ASU.
        cell, grid, unit_cell:
            The crystallographic core (Clipper Cell, Grid_sampling and
            Unit_Cell) for `structure`.
        atoms:
            The core set of atoms to be surrounded by the new selection.
        cutoff:
            Cutoff distance in Angstroms.
        coords (default None):
            Optionally, you can provide the coordinates for all atoms
            (useful for expanding a selection that already includes
            symmetry atoms).
        whole_residues (default true):
            Whether to expand the selections to whole_residues.
        sym_search_frequency (default 2):
            Grid sampling frequency for the symmetry-operator box search.
    '''
    if coords is None:
        coords = atoms.coords
    coords = coords.astype(numpy.float32)
    master_atoms = structure.atoms
    master_coords = master_atoms.coords.astype(numpy.float32)
    from .clipper_util import get_minmax_grid
    grid_minmax = get_minmax_grid(coords, cell, grid)
    from .crystal import calculate_grid_padding
    pad = calculate_grid_padding(cutoff, grid, cell)
    grid_minmax += numpy.array((-pad, pad))
    min_xyz = clipper_python.Coord_grid(grid_minmax[0]).coord_frac(grid).coord_orth(cell).xyz
    dim = grid_minmax[1]-grid_minmax[0]
    symops = unit_cell.all_symops_in_box(min_xyz, dim, True, sym_search_frequency)
    symmats = symops.all_matrices_orth(cell, '3x4')
    from chimerax.geometry import Place
    target = [(coords, Place().matrix.astype(numpy.float32))]
    search_list = []
    for i, s in enumerate(symmats):
        search_list.append((master_coords, s.astype(numpy.float32)))
    from chimerax.geometry import find_close_points_sets
    # We want to respond gracefully if the cutoff is zero, but
    # find_close_points_sets returns nothing if the cutoff is
    # *actually* zero. So just set it to a very small non-zero value.
    if cutoff == 0:
        cutoff = 1e-6
    i1, i2 = find_close_points_sets(search_list, target, cutoff)
    found_atoms = []
    sym_indices = []
    for i, (c, s) in enumerate(zip(i1, symmats)):
        if len(c):
            sel = master_atoms[c]
            if whole_residues:
                sel = sel.unique_residues.atoms
            found_atoms.append(sel)
            indices = numpy.empty(len(sel), numpy.uint8)
            indices[:] = i
            sym_indices.append(indices)
    if len(found_atoms) > 1:
        from chimerax.atomic import concatenate
        found_atoms = concatenate(found_atoms)
        sym_indices = numpy.concatenate(sym_indices)
    elif len(found_atoms) == 0:
        from chimerax.atomic import Atoms
        found_atoms = Atoms()
        sym_indices = numpy.array([], numpy.uint8)
    else:
        found_atoms = found_atoms[0]
        sym_indices = sym_indices[0]
    return (found_atoms, symmats, sym_indices, symops)


def places_from_matrices(matrices, indices, include_identity):
    '''
    Turn a set of orthogonal-space 3x4 symmetry matrices and a collection of
    indices into them into a list of ChimeraX Place objects. The identity
    (ASU) operator lives at index 0 (see Unit_Cell.all_symops_in_box) and is
    placed first when requested; all other operators follow in ascending
    index order so the result is deterministic.
    '''
    wanted = set(int(i) for i in indices)
    if include_identity:
        wanted.add(0)
    else:
        wanted.discard(0)
    ordered = ([0] if (include_identity and 0 in wanted) else []) \
        + sorted(i for i in wanted if i != 0)
    from chimerax.geometry import Place
    return [Place(matrix=matrices[i]) for i in ordered]


def unit_cell_places(cell, spacegroup, na, nb, nc, origin=(0, 0, 0)):
    '''
    Return the ChimeraX Place objects (orthogonal transforms in the structure's
    own coordinate frame) that generate an ``na`` x ``nb`` x ``nc`` block of
    complete unit cells: every space-group symmetry operator combined with each
    lattice translation in [0,na) x [0,nb) x [0,nc), offset by the integer
    ``origin`` cell. The identity operator (the untranslated ASU) comes first, so
    the result can be passed straight to :func:`realize_symmetry_copies`.

    Each copy is placed at ``symop(x)`` plus the whole-cell offset WITHOUT wrapping the
    symop's own translation, so a glide/inversion/screw copy can land outside [0,1) and
    the realised block is spatially scattered across roughly +/-1 cell. This is harmless
    (and intentional) for structure-factor use - Clipper applies the full space group to
    every atom regardless of position, so the SF is placement-invariant. For a compact
    packed block (e.g. a `clipper unitcells` display), pass the realised model through
    :func:`pack_copies_into_block`.
    '''
    from .clipper_python import RTop_fracs, Coord_frac
    oi, oj, ok = origin
    ops = RTop_fracs()
    nsym = spacegroup.num_symops
    for i in range(int(na)):
        for j in range(int(nb)):
            for k in range(int(nc)):
                offset = Coord_frac(oi + i, oj + j, ok + k)
                for s in range(nsym):
                    ops.append(spacegroup.symop(s), offset)
    matrices = ops.all_matrices_orth(cell, '3x4')
    from chimerax.geometry import Place
    return [Place(matrix=m) for m in matrices]


def pack_copies_into_block(box, cell, na, nb, nc, origin=(0, 0, 0)):
    '''
    Rigid-translate each symmetry copy of a realised box (from
    :func:`realize_symmetry_copies` on :func:`unit_cell_places`) by the integer
    lattice vector that brings its centroid into the ``na`` x ``nb`` x ``nc``
    unit-cell block at ``origin`` - turning the raw-symop placement (which scatters
    glide/inversion/screw copies across roughly +/-1 cell) into a compact "packed"
    periodic block for display.

    The shift is the SMALLEST integer lattice vector per axis that lands the copy's
    centroid in ``[origin, origin+dims)`` (so a single cell packs to ``[0,1)`` and a
    larger block keeps its tiling). Because every copy moves by a whole lattice
    vector this is structure-factor-invariant, and because it runs AFTER
    ``realize_symmetry_copies`` (special-position dedup already done) it cannot
    disturb connectivity or the dedup. The same shift is folded into the
    expansion's ``operators`` so ``operators``/``inverse_operators`` stay consistent
    with the packed coordinates. No-op if the box carries no per-atom provenance
    (older expansions). Returns the number of copies moved.
    '''
    import numpy
    from .clipper_python import Coord_orth, Coord_frac
    from chimerax.geometry import translation
    exp = getattr(box, 'clipper_sym_expansion', None)
    if exp is None or exp.operator_index is None:
        return 0
    opi = numpy.asarray(exp.operator_index)
    atoms = box.atoms
    coords = numpy.asarray(atoms.coords, numpy.double)
    frac = numpy.array([Coord_orth(float(x), float(y), float(z)).coord_frac(cell).uvw
                        for x, y, z in coords])
    org = numpy.asarray(origin, float)
    dims = numpy.asarray((na, nb, nc), float)
    o0 = Coord_orth(0., 0., 0.)
    moved = 0
    new_coords = coords.copy()
    for m in range(exp.n_asu):
        sel = (opi == m)
        if not sel.any():
            continue
        c = frac[sel].mean(axis=0)
        n = numpy.zeros(3, int)
        for a in range(3):
            lo, hi = org[a] - c[a], org[a] + dims[a] - c[a]   # valid shift n in [lo, hi)
            if lo <= 0 < hi:
                continue                                       # already in the block
            n[a] = int(numpy.ceil(lo)) if lo > 0 else int(numpy.floor(hi - 1e-9))
        if not n.any():
            continue
        on = Coord_frac(float(n[0]), float(n[1]), float(n[2])).coord_orth(cell)
        dv = numpy.array([on.x - o0.x, on.y - o0.y, on.z - o0.z])
        new_coords[sel] += dv
        exp.operators[m] = translation(dv) * exp.operators[m]
        moved += 1
    if moved:
        atoms.coords = new_coords
    return moved


class CrystalSymmetry:
    '''
    Lightweight, GUI-free holder for a structure's crystallographic core
    (cell, spacegroup, grid, Unit_Cell) with the one operation ``clipper
    symcopies`` needs: enumerating the symmetry operators (as ChimeraX Place
    objects, ASU first) whose copies approach within a cutoff of a selection.
    Build one with :func:`crystal_symmetry_for`.
    '''
    def __init__(self, structure, cell, spacegroup, grid, unit_cell, has_symmetry,
            sym_search_frequency=DEFAULT_SYM_SEARCH_FREQUENCY):
        self.structure = structure
        self.cell = cell
        self.spacegroup = spacegroup
        self.grid = grid
        self.unit_cell = unit_cell
        self.has_symmetry = has_symmetry
        self.sym_search_frequency = sym_search_frequency

    def sym_transforms_near(self, atoms, cutoff, include_identity=True):
        '''
        Return the symmetry operators (as a list of ChimeraX Place objects, in
        the structure's own coordinate frame) for every symmetry copy with at
        least one atom approaching within cutoff (Angstroms) of the given atoms.
        The identity/ASU transform is included first unless include_identity is
        False.
        '''
        found_atoms, symmats, sym_indices, symops = sym_select_within(
            self.structure, self.cell, self.grid, self.unit_cell, atoms, cutoff,
            sym_search_frequency=self.sym_search_frequency)
        return places_from_matrices(symmats, sym_indices, include_identity)


def crystal_symmetry_for(structure, refl_file=None):
    '''
    Build a :class:`CrystalSymmetry` for a ChimeraX structure without requiring
    the GUI SymmetryManager/AtomicSymmetryModel framework. The crystallographic
    core is resolved in priority order:

      1. If `structure` is already managed by a SymmetryManager, reuse its
         cached cell/spacegroup/grid/Unit_Cell (no rebuild).
      2. Else if `refl_file` (a .mtz/.cif structure-factor file) is given, take
         the cell and spacegroup from its header. Useful when the model itself
         carries no crystallographic metadata.
      3. Else build from the model's own metadata (PDB CRYST1 / mmCIF header).

    `CrystalSymmetry.has_symmetry` is False when no real crystallographic
    symmetry could be found (a plain P1 box), so callers can warn and skip.
    '''
    from .symmetry import (get_symmetry_handler, symmetry_from_model_metadata,
        Unit_Cell)
    sh = get_symmetry_handler(structure)
    if sh is not None:
        return CrystalSymmetry(structure, sh.cell, sh.spacegroup, sh.grid,
            sh.unit_cell, sh.has_symmetry)
    if refl_file is not None:
        from .clipper_mtz import load_hkl_data
        from . import Grid_sampling
        hklinfo = load_hkl_data(structure.session, refl_file,
            load_map_coeffs=False)[0]
        cell = hklinfo.cell
        spacegroup = hklinfo.spacegroup
        resolution = hklinfo.resolution
        grid = Grid_sampling(spacegroup, cell, resolution)
        unit_cell = Unit_Cell(structure.atoms, cell, spacegroup, grid)
        return CrystalSymmetry(structure, cell, spacegroup, grid, unit_cell,
            spacegroup.num_symops > 1)
    cell, spacegroup, grid, resolution, has_symmetry = \
        symmetry_from_model_metadata(structure)
    unit_cell = Unit_Cell(structure.atoms, cell, spacegroup, grid)
    return CrystalSymmetry(structure, cell, spacegroup, grid, unit_cell,
        has_symmetry)


class SymmetryExpansion:
    '''
    Book-keeping for an invertible symmetry expansion, attached to the combined
    model as `model.clipper_sym_expansion`. It records exactly which symmetry
    operator produced which output chains, so a caller can collapse the expanded
    model back onto the original ASU by applying each operator's inverse to its
    chains (see `collapse_to_asu`). This is what lets an MD engine relax a full
    P1 box and then fold the result back to a single ASU (at occupancy 1/n_asu)
    for a structure-factor calculation.
    '''
    def __init__(self, operators, chains_by_operator, deduplicated, source_name,
                 source_atom_index=None, operator_index=None):
        # operators[i] is the Place (structure-frame) applied to make copy i, with
        # the identity/ASU operator at index 0. chains_by_operator[i] is the list
        # of chain IDs that copy i contributed to the combined model.
        self.operators = list(operators)
        self.chains_by_operator = [list(c) for c in chains_by_operator]
        # deduplicated: records of special-position fragments dropped at expansion
        # time, each a dict {operator, chain_id, residue, multiplicity}.
        self.deduplicated = list(deduplicated)
        self.source_name = source_name
        # Flat per-box-atom provenance, aligned 1:1 with combined.atoms (same order):
        #   source_atom_index[a] = index of the SOURCE (ASU) atom that box atom a is a
        #       copy of - i.e. its unique-site id (0..n_src-1).
        #   operator_index[a]    = index into self.operators of the symmetry copy that
        #       produced box atom a.
        # Lets a caller count, per unique site, how many box copies exist (and how many
        # are "mobile") to build an occupancy-weighted overlay - see
        # n_copies_per_source_site. None on expansions built before this was recorded.
        self.source_atom_index = (None if source_atom_index is None
                                  else numpy.asarray(source_atom_index, dtype=int))
        self.operator_index = (None if operator_index is None
                               else numpy.asarray(operator_index, dtype=int))

    @property
    def n_asu(self):
        '''Number of ASU copies (symmetry operators) in the expansion.'''
        return len(self.operators)

    @property
    def n_copies_per_source_site(self):
        '''``(n_src,)`` count of box copies of each source (ASU) atom, from
        ``source_atom_index``. To make a site's box copies overlay to a target
        occupancy ``q``, give each copy ``q / n_copies_per_source_site[site]``; the
        uniform ``1/n_asu`` overlay is the special case where every site has ``n_asu``
        copies. For a soft/hard split, restrict the count to the mobile operators:
        ``numpy.bincount(source_atom_index[numpy.isin(operator_index, mobile_ops)])``
        gives the per-site mobile-copy count (the ``n_mobile`` divisor). ``None`` if
        provenance was not recorded.'''
        if self.source_atom_index is None:
            return None
        return numpy.bincount(self.source_atom_index)

    @property
    def inverse_operators(self):
        '''Place objects that undo each operator (map a copy back onto the ASU).'''
        return [op.inverse() for op in self.operators]

    @property
    def operator_by_chain(self):
        '''{chain_id: operator index} - the inverse of chains_by_operator.'''
        d = {}
        for i, chains in enumerate(self.chains_by_operator):
            for cid in chains:
                d[cid] = i
        return d


def _next_free_chain_id(used):
    '''Lowest chain ID (A, B, ... AA, ...) not already in the `used` set.'''
    from chimerax.atomic import next_chain_id
    cid = 'A'
    while cid in used:
        cid = next_chain_id(cid)
    return cid


def _assign_unique_chain_ids(copies):
    '''
    Give every copy globally-unique chain IDs (in place) so ChimeraX's `combine`
    never has to remap them - which both makes the {operator: chain_ids}
    provenance exact and avoids the per-chain "Remapping..." log flood. Copy 0
    (the identity/ASU) keeps its original non-blank IDs; every other chain is
    reassigned to the next free ID. Returns chains_by_copy: a list, per copy, of
    the chain IDs it now owns.

    Polymeric residues take their chain ID from their Chain (ChimeraX refuses a
    direct residue-level write), so those are renamed via the Chain object; free
    residues (waters, ions, ligands) are set through residue chain_ids.
    '''
    used = set()
    chains_by_copy = []
    for i, c in enumerate(copies):
        residues = c.residues
        orig_ids = sorted(set(residues.chain_ids))
        # Avoid this copy's own current IDs too, so a fresh ID can't collide with a
        # chain in this copy that hasn't been renamed yet (e.g. A -> B while B exists).
        avoid = set(used) | set(orig_ids)
        mapping = {}
        for cid in orig_ids:
            if i == 0 and cid and not cid.isspace() and cid not in used:
                new = cid
            else:
                new = _next_free_chain_id(avoid)
                avoid.add(new)
            mapping[cid] = new
            used.add(new)
        for chain in c.chains:
            new = mapping[chain.chain_id]
            if new != chain.chain_id:
                chain.chain_id = new
        free_mask = numpy.array([r.chain is None for r in residues], bool)
        if free_mask.any():
            free = residues[free_mask]
            free.chain_ids = [mapping[cid] for cid in free.chain_ids]
        chains_by_copy.append(sorted(set(mapping.values())))
    return chains_by_copy


# ChimeraX's `combine` drops per-atom custom attributes (unlike `copy`), so the
# physical ones worth keeping are re-applied to the merged model. Only
# `clipper_scattering_species` (the ionic scattering type, e.g. "O2-") is carried:
# it is a genuine per-atom property a later structure-factor calculation needs.
#
# `clipper_sf_exclude` is deliberately NOT preserved. It exists only so Clipper's
# single-ASU SFcalc skips symmetry-completed atoms it would otherwise regenerate;
# once we expand and de-duplicate, those atoms are explicit, distinct and correctly
# counted, so the flag is obsolete - and keeping it would wrongly drop real atoms
# from a SF calc on the expanded/collapsed model.
_CLIPPER_ATOM_ATTRS = ('clipper_scattering_species',)


def _propagate_clipper_attrs(copies, combined):
    '''
    Re-stamp the clipper per-atom custom attributes onto `combined` after a
    `combine`. Because every copy was given unique chain IDs, `combine` merges no
    chains, so combined.atoms is exactly the copies' atoms concatenated in order -
    letting us copy the attribute values positionally. Only attributes actually
    present on the source atoms are touched (proteins carry none, so skip).
    '''
    catoms = combined.atoms
    if len(catoms) != sum(c.num_atoms for c in copies):
        return  # concatenation assumption broken; do not risk mis-assigning
    for attr in _CLIPPER_ATOM_ATTRS:
        values = []
        present = False
        for c in copies:
            for a in c.atoms:
                v = getattr(a, attr, None)
                present = present or v is not None
                values.append(v)
        if not present:
            continue
        for a, v in zip(catoms, values):
            if v is not None:
                setattr(a, attr, v)


# ---------------------------------------------------------------------------
# TODO(ChimeraX aniso_u6 copy bug): AtomicStructure.copy() AND combine() both
# silently drop aniso_u6 (verified 2026-07-09; bfactor/occupancy survive, custom
# attrs survive copy() but not combine()). Reported upstream, and the ChimeraX team
# have PROMISED a fix (per TC, 2026-07-13) -- so REVISIT and REMOVE this workaround
# once that lands (check whether copy()/combine() preserve aniso_u6, then delete the
# source-based re-stamp). Because the ADPs never reach the copied/combined atoms, we
# re-stamp them from the SOURCE structure (rotated by each operator; see
# _propagate_aniso and the collapse_to_asu re-stamp). When the upstream fix lands,
# this whole source-based re-stamp can be replaced by a straight per-copy rotate
# (rotate the copy's own aniso_u6 in place, no surv_masks / source lookup needed).
# ---------------------------------------------------------------------------
def _rotate_aniso_array(u6, place):
    '''
    Rotate an (N,6) array of orthogonal `aniso_u6` ADPs by the operator `place`,
    i.e. U' = R U R^T, so ellipsoids stay consistent with coordinates the same
    operator moves. Uses Clipper's `U_aniso_orth.transform` (the ADP transform
    worked out for the "move symmetry fragment" tool). `place` here is always a
    `Place(matrix=...)` built from a Clipper orthogonal matrix, so `place.matrix`
    feeds straight into `RTop_orth`.
    '''
    from .clipper_python import U_aniso_orth, RTop_orth
    m = numpy.asarray(place.matrix, float)
    rtop = RTop_orth(numpy.ascontiguousarray(m[:, :3]),
                     numpy.ascontiguousarray(m[:, 3]))
    u6 = numpy.asarray(u6, float)
    out = numpy.empty_like(u6)
    for i in range(len(u6)):
        out[i] = U_aniso_orth(numpy.ascontiguousarray(u6[i])).transform(rtop).as_numpy()
    return out


def _propagate_aniso(structure, kept_places, surv_masks, combined):
    '''
    Re-stamp anisotropic ADPs onto the merged model, rotating each copy's by its
    operator (U' = R U R^T). Both `copy()` and `combine` drop `aniso_u6`, so the
    values are taken from the SOURCE structure (not the copies) and each copy's
    contribution is its surviving source atoms (`surv_masks[i]`, a bool mask over
    source atoms). Relies on the unique-chain-ID invariant that makes
    `combined.atoms` exactly the copies' surviving atoms concatenated in order.
    Per-atom presence is honoured (only atoms carrying an anisotropic U are written).
    '''
    src = structure.atoms
    src_has = numpy.asarray(src.has_aniso_u)
    if not src_has.any():
        return
    src_u6 = numpy.full((len(src), 6), numpy.nan)
    src_u6[src_has] = src.filter(src.has_aniso_u).aniso_u6

    has_parts, u6_parts = [], []
    for place, mask in zip(kept_places, surv_masks):
        block_has = src_has[mask]
        block_u6 = src_u6[mask].copy()
        present = numpy.asarray(block_has)
        if present.any():
            block_u6[present] = _rotate_aniso_array(block_u6[present], place)
        has_parts.append(block_has)
        u6_parts.append(block_u6)
    has_all = numpy.concatenate(has_parts)
    u6_all = numpy.concatenate(u6_parts)

    catoms = combined.atoms
    if len(catoms) != len(has_all):
        return
    catoms[has_all].aniso_u6 = u6_all[has_all].astype(numpy.float32)


def _redundant_fragments(structure_copy, matched_atoms, multiplicities, op_index):
    '''
    Among the connected molecules of a realised symmetry copy, find those that are
    special-position duplicates: every atom coincides with an atom already realised
    (`matched_atoms`), and the fragment sits on a special position (some atom has
    site multiplicity > 1). Coincidence of a whole fragment can only happen when a
    symmetry operator maps it onto an existing image - i.e. a special position - so
    unlike the old geometric heuristic this needs no molecule-size cap and no
    polymer exception (a metal on a crystallographic axis in a protein is handled
    too). `multiplicities` is the per-atom site multiplicity of the source ASU
    (aligned to structure_copy.atoms, which preserves source order).

    Returns (atoms_to_delete_or_None, dedup_records, straddle_warn):
      - dedup_records: one dict per dropped fragment residue.
      - straddle_warn: True if a partially-coincident fragment touching a special
        position was seen (the "molecule split across a special position but not
        symmetry-completed" case - imperfectly handled; caller should warn).
    '''
    from chimerax.atomic import concatenate
    all_atoms = structure_copy.atoms
    to_delete = []
    records = []
    straddle_warn = False
    for mol in structure_copy.molecules:
        n_match = len(mol.intersect(matched_atoms))
        if n_match == 0:
            continue
        idx = all_atoms.indices(mol)
        mol_mult = multiplicities[idx]
        special = mol_mult > 1
        if n_match == len(mol):
            if not special.any():
                # Fully coincident but no atom on a special position: an exact
                # overlap we cannot explain crystallographically - keep it rather
                # than silently deleting real atoms.
                continue
            to_delete.append(mol)
            mult = int(mol_mult.max())
            for r in mol.unique_residues:
                records.append({'operator': op_index, 'chain_id': r.chain_id,
                    'residue': '{} {}'.format(r.name, r.number),
                    'multiplicity': mult})
        elif special.any():
            # Some (not all) atoms coincide, and the fragment touches a special
            # position: the on-special atom is duplicated but its off-special
            # neighbours are not - the classic un-completed straddling molecule.
            straddle_warn = True
    if not to_delete:
        return None, records, straddle_warn
    merged = to_delete[0] if len(to_delete) == 1 else concatenate(to_delete)
    return merged, records, straddle_warn


def realize_symmetry_copies(session, structure, places, name=None,
        prune_special_positions=True, tolerance=0.5, multiplicities=None):
    '''
    Build real, whole-model copies of `structure` under each of `places` (a list
    of ChimeraX Place objects in the structure's own coordinate frame, with the
    identity/ASU operator first) and merge them into a single new AtomicStructure
    via ChimeraX's `combine` mechanism.

    Coordinates are never touched directly: `combine` bakes each source
    structure's scene_position into the merged result, so it is enough to make a
    plain copy per operator and set that copy's position to the operator. The
    merged model is then placed at the original's scene_position so it lands
    exactly where the crystallographic copies belong.

    Each copy is given globally-unique chain IDs before combining, so `combine`
    never remaps them. The combined model carries a `clipper_sym_expansion`
    (`SymmetryExpansion`) recording which operator produced which chains, making
    the expansion invertible (see `collapse_to_asu`).

    When `prune_special_positions` is True (default) and per-atom `multiplicities`
    are supplied (site multiplicity of the source ASU, e.g. from
    `clipper_util.site_multiplicities`, aligned to `structure.atoms`), fragment
    copies that fully coincide with an already-realised image at a special
    position are dropped - exactly the redundant `m-1` of every `m`-fold site, so
    a later collapse to the ASU reconstitutes the correct 1/m occupancy.
    Coincidence is tested within `tolerance` Angstroms, cumulatively.

    Returns the combined model (with `.clipper_sym_expansion` set), or None if no
    genuine symmetry copy survives.
    '''
    if not places:
        return None
    from chimerax.geometry import find_close_points

    if prune_special_positions and multiplicities is None:
        session.logger.warning('No site-multiplicity information supplied; '
            'special-position de-duplication skipped for {}.'.format(structure.name))
        prune_special_positions = False
    if multiplicities is not None:
        multiplicities = numpy.asarray(multiplicities)

    copies = []
    kept_places = []
    ref_coords = None
    dedup_records = []
    straddle_warn = False
    # Per kept copy, a boolean mask over the SOURCE atoms marking which survived
    # dedup - needed to re-stamp aniso_u6 from the source (copy()/combine drop it).
    surv_masks = []
    n_src = structure.num_atoms
    for i, place in enumerate(places):
        c = structure.copy()
        # Not yet in the session, so scene_position == position; `combine` reads
        # scene_position to bake in the transform.
        c.position = place
        surv = numpy.ones(n_src, bool)
        if prune_special_positions:
            world = place.transform_points(c.atoms.coords).astype(numpy.float32)
            if i == 0:
                # Identity/ASU: kept whole, seeds the coincidence reference set.
                ref_coords = world
            else:
                i1, i2 = find_close_points(world, ref_coords, tolerance)
                if len(i1):
                    to_delete, recs, sw = _redundant_fragments(
                        c, c.atoms[i1], multiplicities, i)
                    straddle_warn = straddle_warn or sw
                    dedup_records.extend(recs)
                    if to_delete is not None and len(to_delete):
                        # c.atoms is still 1:1 with the source here (pre-deletion).
                        surv[c.atoms.indices(to_delete)] = False
                        to_delete.delete()
                if not c.num_atoms:
                    c.delete()
                    continue
                ref_coords = numpy.concatenate([ref_coords,
                    place.transform_points(c.atoms.coords).astype(numpy.float32)])
        copies.append(c)
        kept_places.append(place)
        surv_masks.append(surv)

    if len(copies) <= 1:
        for c in copies:
            c.delete()
        return None

    chains_by_copy = _assign_unique_chain_ids(copies)

    from chimerax.atomic.cmd import combine_cmd
    if name is None:
        name = '{} symmetry copies'.format(structure.name)
    combined = combine_cmd(session, copies, close=False, name=name,
        add_to_session=False)
    _propagate_clipper_attrs(copies, combined)
    _propagate_aniso(structure, kept_places, surv_masks, combined)
    for c in copies:
        c.delete()
    combined.position = structure.scene_position
    session.models.add([combined])

    # Flat per-box-atom provenance, in combined.atoms order (= surviving source atoms of
    # each kept copy, concatenated in kept_places order - the same invariant
    # _propagate_aniso relies on). source_atom_index = which source ASU atom each box
    # atom copies; operator_index = which kept copy produced it.
    src_index_parts, op_index_parts = [], []
    for op_idx, mask in enumerate(surv_masks):
        surviving = numpy.nonzero(mask)[0]
        src_index_parts.append(surviving)
        op_index_parts.append(numpy.full(len(surviving), op_idx, dtype=int))
    source_atom_index = (numpy.concatenate(src_index_parts) if src_index_parts
                         else numpy.empty(0, int))
    operator_index = (numpy.concatenate(op_index_parts) if op_index_parts
                      else numpy.empty(0, int))

    combined.clipper_sym_expansion = SymmetryExpansion(
        kept_places, chains_by_copy, dedup_records, structure.name,
        source_atom_index=source_atom_index, operator_index=operator_index)

    if dedup_records:
        session.logger.info('De-duplicated {} special-position residue '
            'cop{} while expanding {}.'.format(len(dedup_records),
                'y' if len(dedup_records) == 1 else 'ies', structure.name))
    if straddle_warn:
        session.logger.warning('{} appears to contain molecules split across a '
            'special position that are not symmetry-complete; expansion may leave '
            'partial duplicates. Run "clipper fragments complete" first for exact '
            'special-position handling.'.format(structure.name))
    return combined


def collapse_to_asu(session, model, name=None, set_occupancies=True):
    '''
    Invert a symmetry expansion: fold every copy in `model` back onto the original
    ASU by applying the inverse of the operator that produced it, overlaying all
    `n_asu` images in the ASU frame. `model` must carry a `clipper_sym_expansion`
    (i.e. be the output of `realize_symmetry_copies` / `clipper symcopies` /
    `clipper unitcells`).

    Coordinates are taken as they currently stand, so a box whose atoms have moved
    (e.g. after an MD relaxation) collapses to the corresponding ensemble of ASU
    conformers. With `set_occupancies` True (default) every atom's occupancy is set
    to 1/n_asu, so a general-position atom's `n_asu` overlaid images sum to full
    occupancy and a special-position atom (present only `n_asu/m` times after
    de-duplication) sums to 1/m - exactly the occupancy Clipper's SFcalc expects.
    Anisotropic ADPs are rotated by the same inverse operator (re-stamped from the
    box, since `copy` drops `aniso_u6`); `clipper_scattering_species` rides along
    (preserved by `copy`).

    Returns the new collapsed model. This is a reference implementation; a
    differentiable pipeline can reproduce the same transform from
    `model.clipper_sym_expansion.operators`.
    '''
    from chimerax.core.errors import UserError
    exp = getattr(model, 'clipper_sym_expansion', None)
    if exp is None:
        raise UserError('Model #{} has no symmetry-expansion provenance to invert '
            '(clipper_sym_expansion). It must come from clipper symcopies / '
            'unitcells.'.format(model.id_string))

    # Capture the box ADPs up front: model.copy() drops aniso_u6 (see the
    # "ChimeraX aniso_u6 copy bug" note above _rotate_aniso_array), so they are
    # re-stamped (rotated by the inverse operator) onto the collapsed atoms below.
    box_atoms = model.atoms
    box_has = numpy.asarray(box_atoms.has_aniso_u)
    box_u6 = None
    if box_has.any():
        box_u6 = numpy.full((len(box_atoms), 6), numpy.nan)
        box_u6[box_has] = box_atoms.filter(box_atoms.has_aniso_u).aniso_u6

    collapsed = model.copy(name=name or '{} collapsed to ASU'.format(model.name))
    c_atoms = collapsed.atoms                       # 1:1 with box_atoms (order kept)
    obc = exp.operator_by_chain
    op_of_atom = numpy.array([obc.get(cid, 0) for cid in c_atoms.residues.chain_ids])
    for i, inv in enumerate(exp.inverse_operators):
        sel = op_of_atom == i
        if not sel.any():
            continue
        ai = c_atoms[sel]
        ai.coords = inv.transform_points(ai.coords)
        if box_u6 is not None:
            bh = box_has[sel]
            if bh.any():
                ai.filter(bh).aniso_u6 = _rotate_aniso_array(
                    box_u6[sel][bh], inv).astype(numpy.float32)
    if set_occupancies and collapsed.num_atoms:
        collapsed.atoms.occupancies = 1.0 / exp.n_asu
    collapsed.position = model.scene_position
    session.models.add([collapsed])
    return collapsed
