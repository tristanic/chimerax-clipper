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
Differentiable crystallographic targets for machine-learning force-field
training.

Given a set of atoms in their periodic crystal context and experimental X-ray
data, this module computes a scalar crystallographic loss ``L`` together with
its analytic gradient ``dL/dp`` with respect to a runtime-selected subset of the
per-atom parameters ``{X, Y, Z, Uiso, Occ, U11, U22, U33, U12, U13, U23}``. The
gradient is exact (Agarwal 1978 real-space formulation), reusing Clipper's
``AtomShapeFn.rho_grad`` — the same machinery the B-factor/occupancy refiner
already relies on.

The numeric core (``state.XrayTargetState``) speaks only numpy and has **no
dependency on PyTorch** — it is usable on its own for finite-difference checks
or a non-torch optimiser. PyTorch is an *optional* dependency: it is needed only
for the ``targets`` autograd adapter, and is imported lazily (never at module
import time) via :func:`require_torch`. The ChimeraX-Clipper bundle therefore
installs and loads with no torch present.
'''


def require_torch():
    '''
    Return the imported ``torch`` module, or raise an informative error if it is
    not installed.

    PyTorch is deliberately **not** a dependency of the ChimeraX-Clipper bundle,
    so it may be absent. It is only needed for the :mod:`~chimerax.clipper.diff`
    autograd adapter (:mod:`~chimerax.clipper.diff.targets`); the numeric core
    in :mod:`~chimerax.clipper.diff.state` works without it.
    '''
    try:
        import torch
    except ImportError as e:
        raise ImportError(
            'PyTorch is required for chimerax.clipper.diff.targets but is not '
            'installed in this Python environment. Install a torch build '
            ' into ChimeraX\'s Python with "pip install torch"\n'
            'The numeric core (chimerax.clipper.diff.state) does not need torch.'
        ) from e
    return torch


def parameter_types():
    '''
    Return the mapping of parameter name -> ``AtomShapeFn.TYPE`` enum value for
    the eleven refinable per-atom parameters, in canonical order. Callers select
    a subset of these to request the corresponding gradient columns.
    '''
    from ..clipper_python import AtomShapeFn
    T = AtomShapeFn.TYPE
    return {
        'X': T.X, 'Y': T.Y, 'Z': T.Z,
        'Uiso': T.Uiso, 'Occ': T.Occ,
        'U11': T.U11, 'U22': T.U22, 'U33': T.U33,
        'U12': T.U12, 'U13': T.U13, 'U23': T.U23,
    }
