# Clipper plugin to UCSF ChimeraX
# Copyright (C) 2016-2019 Tristan Croll, University of Cambridge
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

'''
Headless persist / fast-restore of an *assembled* small-molecule crystallographic
target (:class:`~chimerax.clipper.diff.state.EnsembleXrayTargetState`) for ML training.

The expensive step in building such a target is the crystallographic symmetry expansion
(:func:`~chimerax.clipper.sym_realize.realize_symmetry_copies`) that turns the ASU into
the periodic supercell ``box`` (an ``AtomicStructure``). Everything else — reading the
reflections (:func:`~chimerax.clipper.io.small_molecule.fobs_from_arrays`) and building
the state from the box (:func:`~chimerax.clipper.diff.crystal.ensemble_target_from_box`) —
is cheap.

The ``box`` is a plain ``AtomicStructure``, which round-trips losslessly through a native
ChimeraX session (``.cxs``) in ``--nogui`` — coordinates, bonds, occupancy, B-factors,
**anisotropic U**, and the registered per-atom ``clipper_scattering_species`` all survive.
So the fast path is: **save the expanded box in a session; restore it as stored data (no
re-expansion); rebuild the cheap half.** The live ``SymmetryExpansion`` Python object does
NOT survive the session, so the one scalar the target needs from it — ``n_asu`` — is
persisted as a registered structure attribute (along with the crystal definition), leaving
the restored box self-describing.

Typical use (consumer side, e.g. a geometry-store cache tier):

    # --- build + persist (cache miss) ---
    state, box = small_molecule_ensemble_target(session, cif, hkl_path=hkl, ...)
    stamp_box(session, box, cell, spacegroup, resolution, radiation=radiation)
    save_assembled_session(session, cache_path)      # native, --nogui-safe

    # --- fast restore (cache hit; fresh process) ---
    (box,) = load_assembled_structures(session, cache_path)   # no realize_symmetry_copies
    state, box = assembled_target_from_box(session, box, hkl, fsq, sig)

NOTE: a native session save writes the WHOLE session (there is no model-subset save).
Packaging — one combined session for the corpus, or one session per crystal built in an
isolated :class:`~chimerax.core.session.Session` at cache-miss time — is the caller's
choice. Do NOT isolate a box with ``copy()``/``combine`` first: those drop ``aniso_u6``.
'''

# Session-saved structure-level attributes that make a restored box self-describing.
_ASSEMBLED_ATTRS = (
    ('clipper_n_asu', int),
    ('clipper_cell', list),            # [a, b, c, alpha, beta, gamma]
    ('clipper_spacegroup_hall', str),
    ('clipper_resolution', float),
    ('clipper_radiation', str),
)


def register_assembled_attrs(session):
    '''Register (once per session, idempotent) the structure-level attributes stamped on an
    assembled box. Registration makes them round-trip through ``.cxs`` and be typed on
    restore; the values survive regardless, but registering keeps them first-class.'''
    if getattr(session, '_clipper_assembled_attrs_registered', False):
        return
    from chimerax.atomic import AtomicStructure
    for name, atype in _ASSEMBLED_ATTRS:
        AtomicStructure.register_attr(session, name, 'ChimeraX-Clipper', attr_type=atype)
    session._clipper_assembled_attrs_registered = True


def stamp_box(session, box, cell, spacegroup, resolution, radiation='xray', n_asu=None):
    '''Stamp the crystal definition + ``n_asu`` onto an expanded ``box`` so a session-restored
    copy can rebuild the target without re-expanding symmetry. ``n_asu`` defaults to
    ``box.clipper_sym_expansion.n_asu``. Returns ``box``.'''
    register_assembled_attrs(session)
    if n_asu is None:
        exp = getattr(box, 'clipper_sym_expansion', None)
        if exp is None:
            raise ValueError('box has no clipper_sym_expansion; pass n_asu explicitly.')
        n_asu = exp.n_asu
    box.clipper_n_asu = int(n_asu)
    box.clipper_cell = [float(x) for x in (list(cell.dim) + list(cell.angles_deg))]
    box.clipper_spacegroup_hall = spacegroup.symbol_hall
    box.clipper_resolution = float(resolution)
    box.clipper_radiation = str(radiation)
    return box


def save_assembled_session(session, path):
    '''Save the current session to ``path`` via ChimeraX's native session serializer
    (works in ``--nogui``). Stamp the assembled box(es) first (see :func:`stamp_box`) so
    the restore is self-contained. Saves the WHOLE session (no model-subset save).'''
    from chimerax.core.commands import run
    run(session, 'save "%s"' % path)


def load_assembled_structures(session, path):
    '''Open a session saved by :func:`save_assembled_session` (headless) and return the
    restored assembled box structures — those carrying ``clipper_n_asu``. No symmetry
    expansion is performed: the boxes are reloaded as stored atom data.'''
    from chimerax.core.commands import run
    from chimerax.atomic import AtomicStructure
    run(session, 'open "%s"' % path)
    return [m for m in session.models.list(type=AtomicStructure)
            if getattr(m, 'clipper_n_asu', None) is not None]


def _cell_spacegroup_from_box(box):
    from .. import Cell, Cell_descr, Spacegroup, Spgr_descr
    c = list(box.clipper_cell)
    cell = Cell(Cell_descr(*c[:3], *c[3:6]))
    spacegroup = Spacegroup(Spgr_descr(box.clipper_spacegroup_hall, Spgr_descr.Hall))
    return cell, spacegroup


def assembled_target_from_box(session, box, hkl, fsq, sig, *,
                              param_names=('X', 'Y', 'Z'), kind='amplitude', n_threads=1):
    '''Rebuild the :class:`~chimerax.clipper.diff.state.EnsembleXrayTargetState` for an
    already-expanded ``box`` (fresh or session-restored) WITHOUT re-running
    ``realize_symmetry_copies``.

    The raw reflection arrays ``(hkl, fsq, sig)`` come from the caller's cache (see
    :func:`~chimerax.clipper.io.small_molecule.fobs_from_arrays`); the crystal definition,
    ``n_asu`` and radiation are read from the box's stamped attributes (see
    :func:`stamp_box`). Returns ``(state, box)``.

    The returned state pins the rebuilt ``HKL_info`` (and the cell/spacegroup it references):
    Clipper ``HKL_data`` holds a *non-owning* pointer to its ``HKL_info``, so it must be kept
    alive for as long as the target is used.'''
    from ..io.small_molecule import fobs_from_arrays
    from .crystal import ensemble_target_from_box
    radiation = getattr(box, 'clipper_radiation', 'xray')
    cell, spacegroup = _cell_spacegroup_from_box(box)
    hkl_info, fobs, phi_fom, usage = fobs_from_arrays(hkl, fsq, sig, cell, spacegroup)
    state = ensemble_target_from_box(
        box, fobs, phi_fom, usage, param_names=param_names, kind=kind,
        radiation=radiation, n_threads=n_threads)
    # Non-owning HKL_data -> HKL_info pointer: keep the owner (and its cell/spacegroup) alive.
    state._restored_hkl_info = hkl_info
    state._restored_cell = cell
    state._restored_spacegroup = spacegroup
    return state, box
