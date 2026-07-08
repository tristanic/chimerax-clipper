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

import contextlib
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


def _redundant_small_molecules(structure_copy, matched_atoms, max_prune_atoms):
    '''
    Among the connected molecules of a realised symmetry copy, find those that
    are safe to discard as special-position duplicates: small, self-contained
    (non-polymer) molecules - ions, waters, small solvent - every atom of which
    coincides with an atom already present (`matched_atoms`). Larger fragments
    or anything belonging to a polymer chain are kept even when they touch,
    because coincidence there is genuine crystallographic interface contact, not
    a duplicated special-position entity. Returns an Atoms of atoms to delete,
    or None.
    '''
    from chimerax.atomic import Residue, concatenate
    to_delete = []
    for mol in structure_copy.molecules:
        if len(mol) > max_prune_atoms:
            continue
        if not (mol.unique_residues.polymer_types == Residue.PT_NONE).all():
            continue
        # Fully coincident with pre-existing atoms?
        if len(mol.intersect(matched_atoms)) == len(mol):
            to_delete.append(mol)
    if not to_delete:
        return None
    if len(to_delete) == 1:
        return to_delete[0]
    return concatenate(to_delete)


@contextlib.contextmanager
def _suppress_chain_remap_log(session):
    '''
    Temporarily swallow the per-chain "Remapping chain ID ..." info messages that
    chimerax.atomic.cmd.combine_cmd emits. Every symmetry copy shares the ASU's
    chain IDs, so combining even a modest block of cells floods the log with one
    line per chain per copy. We install a highest-priority log that consumes only
    those messages (returning True + excludes_other_logs stops them reaching the
    real logs) and passes everything else through untouched.
    '''
    from chimerax.core.logger import HtmlLog
    class _Filter(HtmlLog):
        excludes_other_logs = True
        def log(self, level, msg, image_info, is_html):
            return level == self.LEVEL_INFO and 'Remapping chain ID' in msg
    f = _Filter()
    session.logger.add_log(f)
    try:
        yield
    finally:
        session.logger.remove_log(f)


def realize_symmetry_copies(session, structure, places, name=None,
        prune_special_positions=True, tolerance=0.5, max_prune_atoms=5,
        log_chain_remapping=False):
    '''
    Build real, whole-model copies of `structure` under each of `places` (a list
    of ChimeraX Place objects in the structure's own coordinate frame, with the
    identity/ASU operator first) and merge them into a single new
    AtomicStructure via ChimeraX's `combine` mechanism.

    Coordinates are never touched directly: `combine` bakes each source
    structure's scene_position into the merged result, so it is enough to make a
    plain copy per operator and set that copy's position to the operator. The
    merged model is then placed at the original's scene_position so it lands
    exactly where the crystallographic copies belong.

    When `prune_special_positions` is True (default), atoms on special positions
    that map back onto an already-realised position are removed from the copies,
    but only for small self-contained molecules (ions/water, up to
    `max_prune_atoms` atoms); larger fragments and polymer chains are kept (see
    `_redundant_small_molecules`). Coincidence is tested within `tolerance`
    Angstroms, cumulatively (each surviving copy is compared against the ASU and
    all copies accepted before it).

    Returns the combined model, or None if no genuine symmetry copy survives
    (e.g. every copy was a redundant special-position duplicate).

    `combine` remaps the (shared) chain IDs of every copy and logs one line per
    chain per copy; because that floods the log for large expansions it is
    suppressed by default. Set `log_chain_remapping` True to keep those messages
    (e.g. when the exact chain-ID assignment matters).
    '''
    if not places:
        return None
    from chimerax.geometry import find_close_points

    copies = []
    ref_coords = None
    n_pruned = 0
    for i, place in enumerate(places):
        c = structure.copy()
        # Not yet in the session, so scene_position == position; `combine` reads
        # scene_position to bake in the transform.
        c.position = place
        if i == 0:
            # Identity/ASU: kept whole, seeds the coincidence reference set.
            ref_coords = place.transform_points(c.atoms.coords).astype(numpy.float32)
            copies.append(c)
            continue
        if prune_special_positions:
            world = place.transform_points(c.atoms.coords).astype(numpy.float32)
            i1, i2 = find_close_points(world, ref_coords, tolerance)
            if len(i1):
                to_delete = _redundant_small_molecules(c, c.atoms[i1],
                    max_prune_atoms)
                if to_delete is not None and len(to_delete):
                    n_pruned += len(to_delete.unique_residues)
                    to_delete.delete()
            if not c.num_atoms:
                c.delete()
                continue
            ref_coords = numpy.concatenate([ref_coords,
                place.transform_points(c.atoms.coords).astype(numpy.float32)])
        copies.append(c)

    if len(copies) <= 1:
        for c in copies:
            c.delete()
        return None

    from chimerax.atomic.cmd import combine_cmd
    if name is None:
        name = '{} symmetry copies'.format(structure.name)
    suppressor = (contextlib.nullcontext() if log_chain_remapping
        else _suppress_chain_remap_log(session))
    with suppressor:
        combined = combine_cmd(session, copies, close=False, name=name,
            add_to_session=False)
    for c in copies:
        c.delete()
    combined.position = structure.scene_position
    session.models.add([combined])
    if n_pruned:
        session.logger.info('Pruned {} redundant special-position {} '
            '(ions/water) from the {} symmetry expansion.'.format(n_pruned,
                'molecule' if n_pruned == 1 else 'molecules', structure.name))
    return combined
