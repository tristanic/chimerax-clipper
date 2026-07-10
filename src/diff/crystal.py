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
    folds the box onto a single ASU at occupancy ``1/n_asu`` and scores that
    occupancy-weighted overlay against the observed data through a single-ASU
    (bulk-free) reciprocal SF calc; the box→ASU fold is differentiable, so the
    coordinate gradient propagates straight back to the box atoms.

    The per-atom **inverse** operators come from ``clipper_sym_expansion`` (the same
    ones ``collapse_to_asu`` applies); the target folds box→ASU with them, rotating
    aniso-U tensors by ``U' = R U Rᵀ`` internally (matching ``collapse_to_asu``), so
    the box-frame model values are passed straight through here. Occupancy defaults to
    ``1/n_asu``.

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

    exp = getattr(box_model, 'clipper_sym_expansion', None)
    if exp is None:
        raise ValueError(
            'box_model #{} has no clipper_sym_expansion; build the box with '
            'clipper unitcells / symcopies (or realize_symmetry_copies) so the '
            'invertible expansion record is attached.'.format(
                getattr(box_model, 'id_string', '?')))

    atoms = box_model.atoms
    M = len(atoms)
    elements = ionic_scattering_names(atoms, radiation=radiation)

    # Per-atom operator index (which symmetry copy each atom belongs to), and the
    # per-atom inverse operator that folds it back onto the ASU.
    obc = exp.operator_by_chain
    op_of_atom = numpy.array(
        [obc.get(cid, 0) for cid in atoms.residues.chain_ids], dtype=int)
    inv_ops = exp.inverse_operators                    # list[Place]
    inv_mats = [numpy.asarray(p.matrix, numpy.double) for p in inv_ops]  # (3,4)
    inv_R = numpy.stack([inv_mats[o][:, :3] for o in op_of_atom])   # (M,3,3)
    inv_t = numpy.stack([inv_mats[o][:, 3] for o in op_of_atom])    # (M,3)

    # Default BOX-frame ADPs (the target rotates them into the ASU frame itself).
    # Isotropic B → u_iso; anisotropic atoms carry aniso_u6 (in the box frame).
    u_iso = numpy.asarray(atoms.bfactors, numpy.double) * _B2U
    has_aniso = numpy.asarray(atoms.has_aniso_u)
    is_aniso = has_aniso.astype(numpy.uint8)
    u_aniso = numpy.zeros((M, 6), numpy.double)
    if has_aniso.any():
        u_aniso[has_aniso] = atoms.filter(has_aniso).aniso_u6

    return EnsembleXrayTargetState(
        elements, inv_R, inv_t, param_names=param_names, fobs=fobs, phi_fom=phi_fom,
        usage=usage, kind=kind, u_iso=u_iso, u_aniso=u_aniso, is_aniso=is_aniso,
        occupancy=1.0 / exp.n_asu, n_threads=n_threads)
