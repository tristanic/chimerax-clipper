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
ChimeraX-aware constructors for the differentiable crystallographic targets.

The numeric cores in :mod:`.state` speak only numpy (no ChimeraX-atomic
dependency); this module bridges a live ChimeraX model — specifically a
symmetry-expanded P1 box carrying a ``clipper_sym_expansion`` (from
``clipper unitcells`` / ``clipper symcopies`` / ``realize_symmetry_copies``) —
to the arrays those cores need.
'''

import math
import numpy

_B2U = 1.0 / (8.0 * math.pi * math.pi)


def ensemble_target_from_box(box_model, fobs, phi_fom, usage, *,
                             param_names=('X', 'Y', 'Z'), kind='amplitude',
                             radiation='xray', n_threads=1):
    '''
    Build an :class:`~chimerax.clipper.diff.state.EnsembleXrayTargetState` (the
    small-molecule ensemble crystallographic target) from a symmetry-expanded box.

    ``box_model`` is a P1 box produced by the symmetry-expansion engine — it must
    carry ``box_model.clipper_sym_expansion`` (a ``SymmetryExpansion``). The target
    scores the box (occupancy ``1/n_asu``) as an occupancy-weighted overlay against the
    observed data through a single-ASU (bulk-free) reciprocal SF calc. Clipper maps each
    box atom into the ASU internally, so the atoms are passed at their real box positions
    and the gradient propagates straight back to them — no explicit box→ASU fold, and
    ``collapse_to_asu`` is needed only for viewing. (Only ``n_asu`` is read from the
    expansion record here, to set the occupancy.)

    ``param_names`` selects which gradient terms the target returns (default coords
    only; include ``U11..U23`` to estimate anisotropic ADPs, etc.) — see
    :class:`~chimerax.clipper.diff.state.EnsembleXrayTargetState`.

    Args:
        * box_model: the expanded P1 box (``AtomicStructure`` with
          ``clipper_sym_expansion``).
        * fobs / phi_fom / usage: observed structure-factor data for the crystal —
          e.g. ``xm.f_obs`` / ``xm.weights`` / ``xm.usage_flags`` from the live
          small-molecule map manager. Small molecules supply ``fom()==1``.
        * param_names: gradient terms to return (default ``('X','Y','Z')``).
        * kind: ``'amplitude'`` (default) or ``'intensity'``.
        * radiation: ``'xray'`` (default) or ``'electron'`` — selects the
          scattering-factor identifiers.

    Returns the ready target; feed box coordinates (and, if refining them, box ADPs)
    to its ``value_and_gradient(...)`` or wrap it with
    :func:`chimerax.clipper.diff.targets.xray_loss`
    (or ``crystallographic_loss`` for the coords-only default).
    '''
    from ..scattering import ionic_scattering_names
    from .state import EnsembleXrayTargetState

    # Only n_asu is needed here (for the 1/n_asu occupancy overlay). A freshly-expanded
    # box carries the full SymmetryExpansion; a box restored from a session (see
    # diff.assembled_cache) has only the persisted scalar `clipper_n_asu` attribute (the
    # Python expansion object does not survive .cxs), so fall back to that.
    exp = getattr(box_model, 'clipper_sym_expansion', None)
    if exp is not None:
        n_asu = int(exp.n_asu)
    else:
        n_asu = getattr(box_model, 'clipper_n_asu', None)
        if n_asu is None:
            raise ValueError(
                'box_model #{} has neither a clipper_sym_expansion nor a clipper_n_asu '
                'attribute; build the box with clipper unitcells / symcopies (or '
                'realize_symmetry_copies), or restore one saved via '
                'chimerax.clipper.diff.assembled_cache.'.format(
                    getattr(box_model, 'id_string', '?')))
        n_asu = int(n_asu)

    atoms = box_model.atoms
    M = len(atoms)
    elements = ionic_scattering_names(atoms, radiation=radiation)

    # BOX-frame ADPs, passed straight through (Clipper folds each atom into the ASU).
    # Isotropic B → u_iso; anisotropic atoms carry aniso_u6 (in the box frame).
    u_iso = numpy.asarray(atoms.bfactors, numpy.double) * _B2U
    has_aniso = numpy.asarray(atoms.has_aniso_u)
    is_aniso = has_aniso.astype(numpy.uint8)
    u_aniso = numpy.zeros((M, 6), numpy.double)
    if has_aniso.any():
        u_aniso[has_aniso] = atoms.filter(has_aniso).aniso_u6

    return EnsembleXrayTargetState(
        elements, param_names=param_names, fobs=fobs, phi_fom=phi_fom,
        usage=usage, kind=kind, u_iso=u_iso, u_aniso=u_aniso, is_aniso=is_aniso,
        occupancy=1.0 / n_asu, n_threads=n_threads)


def _supercell_multiples(cell, n_cells, min_box_size):
    '''
    Resolve the ``(na, nb, nc)`` integer unit-cell repeats for the expansion box.

    ``min_box_size`` (Angstroms), when given, is the minimum **perpendicular width**
    (face-to-face distance) the box must have along each axis — the quantity a
    periodic-MD minimum-image / non-bonded cutoff constrains. The one-cell width
    along an axis is ``1 / cell.<x>_star`` (the (h00)/(0k0)/(00l) d-spacing), so the
    repeat count is ``ceil(min_box_size * cell.<x>_star)``. When both are given the
    element-wise max is used, so the result satisfies both. Every repeat is >= 1
    (the expansion needs at least one whole cell — a complete set of orbits).
    '''
    import math
    na, nb, nc = (n_cells if n_cells is not None else (1, 1, 1))
    na, nb, nc = int(na), int(nb), int(nc)
    if min_box_size is not None:
        if min_box_size <= 0:
            raise ValueError('min_box_size must be positive (got %r)' % (min_box_size,))
        na = max(na, int(math.ceil(min_box_size * cell.a_star)))
        nb = max(nb, int(math.ceil(min_box_size * cell.b_star)))
        nc = max(nc, int(math.ceil(min_box_size * cell.c_star)))
    return max(na, 1), max(nb, 1), max(nc, 1)


def small_molecule_ensemble_target(session, cif_path, hkl_path=None, *,
                                   param_names=('X', 'Y', 'Z'), kind='amplitude',
                                   radiation='auto', recover_scattered_hydrogens=True,
                                   complete_fragments=True,
                                   n_cells=None, min_box_size=None, n_threads=1):
    '''
    One-call, GUI-free builder: a small-molecule (COD) CIF plus its reflections ->
    a ready :class:`~chimerax.clipper.diff.state.EnsembleXrayTargetState` for
    differentiable structure-factor training. Returns ``(state, box_model)``.

    This performs the SAME crystallographic setup the turnkey ``clipper cod`` /
    ``clipper smallmol`` commands do — minus the live map manager, which is what
    makes those commands GUI-only — so it runs headless (``--nogui``). Feeding the
    pieces by hand is error-prone because several corecif corrections are easy to
    miss; this bakes them in, in the right order:

      1. open the CIF (corecif parser). ``open_small_molecule_cif`` internally
         rebuilds coordinates in Clipper's frame — correcting corecif's oblique-cell
         coordinate error (a NO-OP for orthorhombic cells but essential for
         monoclinic/triclinic, where skipping it silently corrupts the structure
         factors, not merely NaNs them) — and repairs covalent connectivity corecif
         drops on metal-coordinated atoms;
      2. **hydrate** the crystallographic per-atom data corecif omits — isotropic B,
         orthogonal-frame anisotropic U, and the ionic scattering species. Skipping
         this leaves ``Atom.bfactor == 0`` -> U = 0 -> infinitely sharp atoms ->
         **NaN gradients** (with a deceptively finite loss value);
      3. optionally recover hydrogens corecif scattered into the wrong ASU because
         their X-H bond was symmetry-coded in ``_geom_bond`` (``recover_scattered_
         hydrogens``; moves each orphaned H onto the correct symmetry image so it is
         present in the right cell for the SF sum);
      4. optionally complete molecules split across a special position
         (``complete_fragments``; required for correct SFs on such entries);
      5. read the reflections into fobs/phi_fom/usage (small-molecule convention:
         Fo = sqrt(I), fom = 1, all-working);
      6. expand the ASU to a full unit cell (with special-position de-duplication);
      7. build the ensemble target.

    Args:
        * cif_path: the model CIF.
        * hkl_path: reflections (COD ``.hkl`` / CIF ``_refln_`` loop). Default None
          -> the sibling ``.hkl`` of ``cif_path``.
        * param_names: gradient terms to return (default coords only; include
          ``U11..U23`` for anisotropic-ADP gradients, etc.).
        * kind: ``'amplitude'`` (default) or ``'intensity'``.
        * radiation: ``'xray'`` | ``'electron'`` | ``'auto'`` (default; read from the
          CIF's ``_diffrn_radiation_probe``).
        * recover_scattered_hydrogens: relocate symmetry-scattered hydrogens (orphaned
          by a symmetry-coded X-H bond) onto the correct image before expansion
          (default True). Mutating but restricted to unambiguously-matched H.
        * complete_fragments: complete special-position-straddling molecules before
          expansion (default True).
        * min_box_size: minimum box **perpendicular width** (Angstroms, face-to-face)
          required along each cell axis — the periodic-MD minimum-image / non-bonded
          cutoff constraint. The unit cell is repeated ``ceil(min_box_size /
          perp_width)`` times per axis to meet it. This is the usual real-world knob;
          leave ``n_cells`` unset and pass this.
        * n_cells: explicit ``(na, nb, nc)`` unit-cell repeats. Default (both unset):
          a single unit cell. If both ``n_cells`` and ``min_box_size`` are given, the
          per-axis maximum of the two is used.
        * n_threads: SF / gradient worker threads.

    The returned ``box_model`` defines the atom order the target expects
    coordinates in: feed ``box_model.atoms.coords`` (or a torch tensor of them) to
    ``state.value_and_gradient`` or to
    :func:`chimerax.clipper.diff.targets.xray_loss` /
    :func:`~chimerax.clipper.diff.targets.crystallographic_loss`. Both ``box_model``
    and the source model are added to the session; close them when done. The state
    keeps the reflection ``HKL_info`` alive (Clipper HKL_data holds a non-owning
    pointer to it).
    '''
    from chimerax.core.errors import UserError
    from ..symmetry import crystal_symmetry_from_cif_file
    from ..io.small_molecule import (open_small_molecule_cif,
        hydrate_small_molecule_model, _resolve_radiation,
        read_small_molecule_fobs, reassemble_symmetry_scattered_hydrogens)
    from ..sym_realize import unit_cell_places, realize_symmetry_copies
    from ..clipper_util import site_multiplicities

    radiation = _resolve_radiation(radiation, cif_path)
    # open corrects the oblique-cell coordinate frame + repairs connectivity internally.
    model = open_small_molecule_cif(session, cif_path)
    cell, spacegroup, grid = crystal_symmetry_from_cif_file(cif_path)
    if spacegroup.num_symops <= 1:
        raise ValueError('%r is P1 (no crystallographic symmetry); the ensemble '
                         'target needs a symmetry-expanded cell.' % cif_path)
    # fill the crystallographic per-atom data corecif omits (B / aniso U / ionic species).
    hydrate_small_molecule_model(session, model, cif_path, cell, radiation)
    # Bring home hydrogens corecif scattered into the wrong ASU (their bond to a heavy
    # atom was symmetry-coded in _geom_bond). Mutating, so opt-in; runs after open's
    # repair_connectivity (only genuinely orphaned H remain) and BEFORE fragment
    # completion / expansion, so the covalent fragments are assembled from correctly
    # placed atoms.
    if recover_scattered_hydrogens:
        n_h = reassemble_symmetry_scattered_hydrogens(session, model, cell, spacegroup)
        if n_h:
            session.logger.info('(CLIPPER) relocated %d symmetry-scattered '
                                'hydrogen(s) in %s' % (n_h, cif_path))
    if complete_fragments:
        from ..io.fragments import split_fragments
        split_fragments(session, model, cell, spacegroup, grid, mode='complete',
                        path=cif_path, log=session.logger)

    hkl_info, fobs, phi_fom, usage = read_small_molecule_fobs(
        hkl_path or cif_path, cell, spacegroup)

    na, nb, nc = _supercell_multiples(cell, n_cells, min_box_size)
    if (na, nb, nc) != (1, 1, 1):
        session.logger.info('(CLIPPER) expanding to a %dx%dx%d supercell '
            '(perp. widths %.1f/%.1f/%.1f A)' % (na, nb, nc,
            na / cell.a_star, nb / cell.b_star, nc / cell.c_star))
    mult = site_multiplicities(model.atoms.coords, cell, spacegroup, grid)
    places = unit_cell_places(cell, spacegroup, na, nb, nc)
    box = realize_symmetry_copies(session, model, places, multiplicities=mult)
    if box is None:
        raise UserError('Symmetry expansion of %r produced no copies.' % cif_path)

    state = ensemble_target_from_box(box, fobs, phi_fom, usage,
                                     param_names=param_names, kind=kind,
                                     radiation=radiation, n_threads=n_threads)
    # Keep the HKL_info alive for the lifetime of the state: Clipper HKL_data holds a
    # non-owning pointer to it, and the state's own refs cover fobs/phi_fom/usage but
    # not the HKL_info that owns their reflection list.
    state._hkl_info_keepalive = hkl_info
    return state, box
