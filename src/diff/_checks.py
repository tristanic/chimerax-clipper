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
Finite-difference gradient checks — the primary correctness gate for every
differentiable target (torch-free).

Compares the analytic ``dL/dx`` returned by a ``*TargetState`` against a central
finite difference of its own ``L``. A target is trustworthy only once this
passes to ~1e-4 relative; nothing downstream (the torch adapter, a training run)
should be believed before it does.
'''

import numpy


def finite_difference_coord_gradient(state, coords, eps=1e-3):
    '''
    Central-difference estimate of ``dL/dx`` for a coordinate-differentiable
    target: perturb each Cartesian component of each atom by +/- ``eps`` (A) and
    difference ``state.value_and_gradient(...)[0]``.

    Returns an ``(N, 3)`` float64 array.
    '''
    coords = numpy.array(coords, numpy.double, copy=True)
    n = len(coords)
    num = numpy.empty((n, 3), numpy.double)
    for j in range(n):
        for k in range(3):
            saved = coords[j, k]
            coords[j, k] = saved + eps
            lp = state.value_and_gradient(coords)[0]
            coords[j, k] = saved - eps
            lm = state.value_and_gradient(coords)[0]
            coords[j, k] = saved
            num[j, k] = (lp - lm) / (2.0 * eps)
    return num


def check_coord_gradient(state, coords, eps=1e-3, rtol=1e-2):
    '''
    Compare analytic vs finite-difference coordinate gradients, using a metric
    that is meaningful for a float32 map.

    The pass criterion is **field-scaled**, not per-component: the largest
    absolute discrepancy across the whole gradient array must be within ``rtol``
    of the largest analytic gradient magnitude::

        passed = max|analytic - numeric| <= rtol * max|analytic|

    This is the right question to ask because (a) a central finite difference on
    a float32 map has a noise floor of roughly ``map_ULP / (2*eps)`` in gradient
    units — a few x1e-3 at ``eps=1e-3`` — so sub-1e-3 agreement is unattainable
    regardless of gradient correctness; and (b) atoms on special positions have
    components that are *exactly* zero by symmetry, for which a per-component
    relative error is meaningless (0/0). Comparing the worst error to the field
    scale sidesteps both. ``rtol`` defaults to 1e-2 (1%); real agreement for a
    correct analytic gradient is typically ~0.1%.

    The returned dict also reports, per component, how many "significant"
    components (|analytic| >= 5% of the field scale) agree to within ``rtol``,
    and the worst offender among them — this is the diagnostic that actually
    localises a real gradient bug (a genuine bug perturbs the *large* components,
    not just the near-zero ones).

    Returns a dict with ``passed``, ``field_scale``, ``max_abs_error``,
    ``rel_to_field`` (= max_abs_error / field_scale), ``n_significant``,
    ``n_significant_ok``, ``worst_significant`` ((atom, axis, analytic, numeric)
    or None), and the raw ``analytic`` / ``numeric`` arrays.
    '''
    coords = numpy.asarray(coords, numpy.double)
    _, analytic = state.value_and_gradient(coords)
    analytic = numpy.asarray(analytic, numpy.double)
    numeric = finite_difference_coord_gradient(state, coords, eps=eps)
    abs_err = numpy.abs(analytic - numeric)

    field_scale = float(numpy.abs(analytic).max()) if analytic.size else 0.0
    max_abs = float(abs_err.max()) if abs_err.size else 0.0
    rel_to_field = (max_abs / field_scale) if field_scale > 0 else 0.0

    # Per-component agreement, assessed only where the analytic gradient carries
    # real signal (>= 5% of the field scale). A correct gradient agrees here; a
    # bug shows up here first.
    sig_floor = 0.05 * field_scale
    significant = numpy.abs(analytic) >= sig_floor
    n_sig = int(significant.sum())
    worst = None
    n_sig_ok = n_sig
    if n_sig:
        sig_rel = numpy.where(significant, abs_err / numpy.maximum(numpy.abs(analytic), 1e-30), 0.0)
        n_sig_ok = int((sig_rel[significant] <= rtol).sum())
        flat = numpy.argmax(sig_rel)
        a, ax = divmod(int(flat), 3)
        worst = (a, 'XYZ'[ax], float(analytic[a, ax]), float(numeric[a, ax]))

    return {
        'passed': rel_to_field <= rtol,
        'field_scale': field_scale,
        'max_abs_error': max_abs,
        'rel_to_field': rel_to_field,
        'n_significant': n_sig,
        'n_significant_ok': n_sig_ok,
        'worst_significant': worst,
        'analytic': analytic,
        'numeric': numeric,
    }


# ---------------------------------------------------------------------------
# Multi-parameter checks for XrayTargetState (coords + Uiso + Occ + anisotropic U)
# ---------------------------------------------------------------------------

#: Default finite-difference step per parameter (Å for coords, Å² for U, unitless
#: for Occ). ADP steps are smaller because U enters the Gaussian exponent.
_DEFAULT_EPS = {
    'X': 1e-3, 'Y': 1e-3, 'Z': 1e-3,
    'Uiso': 1e-4, 'Occ': 1e-3,
    'U11': 1e-4, 'U22': 1e-4, 'U33': 1e-4,
    'U12': 1e-4, 'U13': 1e-4, 'U23': 1e-4,
}

#: Which mutable array (and column, or None for a 1-D array) each parameter
#: perturbs, for the finite difference.
_PARAM_SLOT = {
    'X': ('coords', 0), 'Y': ('coords', 1), 'Z': ('coords', 2),
    'Uiso': ('u_iso', None), 'Occ': ('occ', None),
    'U11': ('u_aniso', 0), 'U22': ('u_aniso', 1), 'U33': ('u_aniso', 2),
    'U12': ('u_aniso', 3), 'U13': ('u_aniso', 4), 'U23': ('u_aniso', 5),
}


def finite_difference_param_gradient(state, coords, u_iso=None, u_aniso=None,
                                     occ=None, is_aniso=None, eps_map=None):
    '''
    Central-difference estimate of ``dL/dp`` for every parameter in
    ``state.param_names``, for an :class:`~.state.XrayTargetState`.

    Every evaluation uses ``refresh_scale=False``, so the finite difference
    measures only the explicit parameter dependence — the caller is responsible
    for having primed (frozen) the scale first (see :func:`check_param_gradients`).
    Priming at a *reference* distinct from the evaluation point is important: a
    self-scaling target primed at the evaluation point makes occupancy / overall-B
    a null direction (the scale absorbs them), so those columns come out ~0.
    Returns an ``(N, P)`` array with columns in ``state.param_names`` order.
    '''
    n = state.n_atoms
    coords = numpy.array(coords, numpy.double, copy=True).reshape(n, 3)
    u_iso = (numpy.zeros(n) if u_iso is None
             else numpy.array(u_iso, numpy.double, copy=True).reshape(n))
    u_aniso = (numpy.zeros((n, 6)) if u_aniso is None
               else numpy.array(u_aniso, numpy.double, copy=True).reshape(n, 6))
    occ = (numpy.ones(n) if occ is None
           else numpy.array(occ, numpy.double, copy=True).reshape(n))
    arrays = {'coords': coords, 'u_iso': u_iso, 'u_aniso': u_aniso, 'occ': occ}
    eps_map = dict(_DEFAULT_EPS, **(eps_map or {}))

    def _L():
        return state.value_and_gradient(coords, u_iso, u_aniso, occ, is_aniso,
                                        refresh_scale=False)[0]

    P = len(state.param_names)
    num = numpy.zeros((n, P), numpy.double)
    for c, name in enumerate(state.param_names):
        key, col = _PARAM_SLOT[name]
        arr = arrays[key]
        eps = eps_map[name]
        for j in range(n):
            if col is None:
                saved = arr[j]
                arr[j] = saved + eps; lp = _L()
                arr[j] = saved - eps; lm = _L()
                arr[j] = saved
            else:
                saved = arr[j, col]
                arr[j, col] = saved + eps; lp = _L()
                arr[j, col] = saved - eps; lm = _L()
                arr[j, col] = saved
            num[j, c] = (lp - lm) / (2.0 * eps)
    return num


def check_param_gradients(state, coords, u_iso=None, u_aniso=None, occ=None,
                          is_aniso=None, eps_map=None, rtol=1e-2, prime=None):
    '''
    Compare analytic vs finite-difference gradients for an
    :class:`~.state.XrayTargetState`, **per parameter column** (the field-scaled
    metric of :func:`check_coord_gradient` applied independently to each column,
    because coordinate and ADP gradients differ by orders of magnitude in both
    units and scale).

    The scale factor is frozen for the whole check so analytic and finite
    difference are mutually consistent. Pass ``prime`` — a dict of reference
    parameters ``{'coords':..., 'u_iso':..., 'u_aniso':..., 'occ':...,
    'is_aniso':...}`` — to freeze the scale at that reference (recommended: the
    target-generating / deposited model). This avoids the scale-degeneracy that
    otherwise zeroes the occupancy / overall-B columns of a self-scaling target.
    If ``prime`` is None the scale is frozen at the evaluation point (fine for
    coordinate-only or fixed-scale targets).

    Returns a dict with an overall ``passed`` (all columns pass) and a
    ``per_param`` mapping name -> {passed, field_scale, max_abs_error,
    rel_to_field, worst_atom}. This is the check that validates the six
    anisotropic-U columns and x,y,z coordinate columns — in particular that the
    factor-2 (diagonal) / factor-4 (off-diagonal) ADP gradient scaling is applied
    exactly once (a wrong factor shows as a clean 2x error in that column).
    '''
    if prime is not None:
        state.value_and_gradient(
            prime['coords'], prime.get('u_iso'), prime.get('u_aniso'),
            prime.get('occ'), prime.get('is_aniso'), refresh_scale=True)
        _, analytic = state.value_and_gradient(
            coords, u_iso, u_aniso, occ, is_aniso, refresh_scale=False)
    else:
        # No reference given: prime + analytic at the evaluation point.
        _, analytic = state.value_and_gradient(
            coords, u_iso, u_aniso, occ, is_aniso, refresh_scale=True)
    analytic = numpy.asarray(analytic, numpy.double)
    numeric = finite_difference_param_gradient(
        state, coords, u_iso, u_aniso, occ, is_aniso, eps_map=eps_map)

    per_param = {}
    all_passed = True
    for c, name in enumerate(state.param_names):
        a = analytic[:, c]
        nmc = numeric[:, c]
        abs_err = numpy.abs(a - nmc)
        field = float(numpy.abs(a).max()) if a.size else 0.0
        max_abs = float(abs_err.max()) if abs_err.size else 0.0
        rel = (max_abs / field) if field > 0 else 0.0
        # If the whole column is ~0 (parameter inactive for these atoms), both
        # analytic and numeric are ~0 and the column trivially passes.
        passed = (rel <= rtol) or (field == 0.0 and max_abs <= 1e-9)
        worst_atom = int(numpy.argmax(abs_err)) if abs_err.size else -1
        per_param[name] = {
            'passed': passed, 'field_scale': field, 'max_abs_error': max_abs,
            'rel_to_field': rel, 'worst_atom': worst_atom,
        }
        all_passed = all_passed and passed
    return {'passed': all_passed, 'per_param': per_param,
            'analytic': analytic, 'numeric': numeric}
