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
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

'''
Finite-difference correctness gate for the differentiable X-ray target
(:mod:`chimerax.clipper.diff`).

This builds *synthetic* crystallographic reference data entirely from the live
Clipper API — a small unit cell, a handful of C/N/O atoms (mixed isotropic and
anisotropic), the structure factors those atoms produce, and the electron
density map they produce — then drives the finite-difference gate in
:mod:`chimerax.clipper.diff._checks` to prove that the analytic Agarwal gradients
returned by :class:`~chimerax.clipper.diff.state.XrayTargetState` (and the cubic-
interpolation gradient of :class:`~chimerax.clipper.diff.state.MapTargetState`)
match a central finite difference of the target's own value.

It is torch-free and needs no GUI/OpenGL, so it runs HEADLESS. There is no live
crystal, no MTZ file and no model — the reference "observations" ARE the
reference model's own Fcalc (reciprocal target) / density (real-space target),
so at the reference parameters the target is ~0 and any perturbation produces a
clean, physically meaningful gradient to difference against.

Coverage
--------
* :func:`~chimerax.clipper.diff._checks.check_coord_gradient` on a
  :class:`MapTargetState` (Tier-0, fixed-map, coordinate-only) — validates the
  ``Xmap.interp_cubic_grad_frac`` coordinate-gradient path against a genuine
  density field.
* :func:`~chimerax.clipper.diff._checks.check_param_gradients` on an
  :class:`XrayTargetState` in BOTH modes (real-space ``target_map`` and
  reciprocal ``fobs``), exercising all eleven per-atom parameter columns
  ``X,Y,Z,Uiso,Occ,U11,U22,U33,U12,U13,U23``. The isotropic atoms carry the
  ``Uiso`` signal; the anisotropic atoms carry the six ``U*`` signals; every
  column carries genuine analytic signal (non-trivial field scale), so a wrong
  ADP factor (the classic diagonal-2 / off-diagonal-4 bug) would show as a clean
  per-column error.
* Both are run for P1 and for the non-trivial monoclinic screw group P2_1, so
  the reciprocal symmetry-multiplicity (``epsilonc``) weighting is exercised.

The XrayTargetState coordinate columns are validated via the X/Y/Z columns of
``check_param_gradients`` (which supplies realistic ADPs) rather than by calling
``check_coord_gradient`` on the XrayTargetState directly: ``check_coord_gradient``
cannot pass a per-atom ``u_iso``, so it would drive the model with U=0 (infinitely
sharp atoms), whose density aliases on the grid and yields degenerate/NaN
gradients — a meaningless check. ``MapTargetState`` is the correct coordinate-only
partner for ``check_coord_gradient``.

How to run
----------
Headless, against an install that carries ``src/diff`` (the daily, from master)::

    "C:/Program Files/ChimeraX-Daily/bin/ChimeraX-console.exe" \
        --nogui --exit --script tests/test_diff_gradients.py

Prints a per-parameter table and a final ``OVERALL: PASS``/``FAIL`` line, and
exits non-zero on failure. The ``test_*`` functions are also plain asserts, so the
module can be driven by pytest from within a ChimeraX Python environment.
'''

import numpy

# (space-group Hermann-Mauguin symbol, cell a,b,c,alpha,beta,gamma) pairs.
# P1 (n_sym=1) is the baseline; P2_1 (n_sym=2, monoclinic screw axis, no special
# positions) exercises the reciprocal epsilonc symmetry weighting.
_SPACEGROUPS = [
    ('P 1',  (30.0, 32.0, 34.0, 90.0, 90.0, 90.0)),
    ('P 21', (28.0, 30.0, 32.0, 90.0, 102.0, 90.0)),
]

_ALL_PARAMS = ('X', 'Y', 'Z', 'Uiso', 'Occ',
               'U11', 'U22', 'U33', 'U12', 'U13', 'U23')


def _build_reference_crystal(sg_symbol, cell_params, res_limit=2.0):
    '''
    Construct a synthetic crystal + reference model and derive, purely from the
    reference model, everything the two XrayTargetState modes need:

      * ``fobs`` (``HKL_data<F_sigF>``) = |Fcalc| of the reference model, sigF=1;
      * ``phi_fom`` (``HKL_data<Phi_fom>``) with fom()=1 (pure 1/sig^2 weighting);
      * ``usage`` (``HKL_data<Flag>``) all-working;
      * ``xmap`` (``Xmap<ftype32>``) = FFT of the reference Fcalc (its density).

    Returns a dict of the crystal objects and the reference per-atom parameters.
    The ``HKL_info`` is returned in the dict and MUST be kept alive for as long as
    the reciprocal state is used (Clipper ``HKL_data`` holds a non-owning pointer
    to it).
    '''
    from chimerax.clipper import (HKL_info, HKL_data_F_sigF, HKL_data_Phi_fom,
                                  HKL_data_Flag)
    from chimerax.clipper.clipper_python import (
        Cell, Cell_descr, Spacegroup, Spgr_descr, Resolution, Grid_sampling,
        Atom_list, Coord_frac, AtomShapeFn, Xmap_float, SFcalc_aniso_sum_float)
    from chimerax.clipper.clipper_python.data32 import HKL_data_F_phi_float

    cell = Cell(Cell_descr(*cell_params))
    sg = Spacegroup(Spgr_descr(sg_symbol))
    res = Resolution(res_limit)
    hklinfo = HKL_info(sg, cell, res, True)          # generate=True -> reflection list
    grid = Grid_sampling(sg, cell, res)

    # Reference model: 5 atoms at general positions; atoms 0-2 isotropic,
    # atoms 3-4 anisotropic (so every parameter column below carries signal).
    elements = ['C', 'N', 'O', 'C', 'O']
    frac = numpy.array([
        [0.15, 0.20, 0.12],
        [0.35, 0.55, 0.30],
        [0.60, 0.25, 0.70],
        [0.80, 0.65, 0.45],
        [0.45, 0.85, 0.60],
    ])
    n = len(elements)
    coords = numpy.array([Coord_frac(*f).coord_orth(cell).xyz for f in frac], numpy.double)
    occ = numpy.ones(n)
    is_aniso = numpy.array([0, 0, 0, 1, 1], numpy.uint8)
    u_iso = numpy.array([0.30, 0.25, 0.35, 0.0, 0.0])
    u_aniso = numpy.zeros((n, 6))
    u_aniso[3] = [0.30, 0.35, 0.28, 0.03, 0.02, 0.010]   # positive-definite (Å²)
    u_aniso[4] = [0.25, 0.22, 0.30, -0.02, 0.01, 0.015]

    # Atom_list convention for SFcalc: anisotropic atoms carry their u_aniso;
    # isotropic atoms carry NaN in the aniso slots (and their u_iso).
    u_iso_al = u_iso.copy()
    u_aniso_al = numpy.full((n, 6), numpy.nan)
    for j in range(n):
        if is_aniso[j]:
            u_aniso_al[j] = u_aniso[j]
    atom_list = Atom_list(elements, coords, occ, u_iso_al, u_aniso_al)

    # Reference Fcalc (exact direct summation, X-ray factors, no bulk solvent).
    # NB: build the functor then call it — the "compute in constructor" shorthand
    # does not fill fcalc through the bindings.
    fcalc = HKL_data_F_phi_float(hklinfo)
    SFcalc_aniso_sum_float(radiation=AtomShapeFn.XRAY)(fcalc, atom_list)
    hkls_i, cdata = fcalc.data                 # cdata[:,0]=|Fc|, [:,1]=phi
    hkls_i = hkls_i.astype(numpy.int32)
    l = len(hkls_i)
    Fc = numpy.nan_to_num(cdata[:, 0]).astype(numpy.float32)

    # Synthetic observations = reference |Fcalc|, sigF=1, fom=1, all-working.
    fobs = HKL_data_F_sigF(hklinfo)
    fobs.set_data(hkls_i, numpy.stack([Fc, numpy.ones(l, numpy.float32)], axis=1))
    phi_fom = HKL_data_Phi_fom(hklinfo)
    phi_fom.set_data(hkls_i, numpy.stack(
        [numpy.zeros(l, numpy.float32), numpy.ones(l, numpy.float32)], axis=1))
    usage = HKL_data_Flag(hklinfo)
    usage.set_data(hkls_i, numpy.ones((l, 1), numpy.int32))

    # Reference density map (real-space target / MapTargetState field).
    xmap = Xmap_float(sg, cell, grid)
    xmap.fft_from(fcalc)

    # Synthetic measured intensities for the iobs (I-direct) target: |Fc|^2 shifted down
    # so the weakest ~15% are NEGATIVE, with a sigma that varies with I (as a SHELXL
    # effective sigma does). fobs_i is the matching amplitude set (Fo = sqrt(I) for I > 0,
    # else 0), which the evaluator still needs for its log-space scale stage.
    from chimerax.clipper import HKL_data_I_sigI
    i2 = Fc.astype(numpy.double) ** 2
    i_obs = i2 - numpy.percentile(i2, 15)
    s_obs = 0.05 * i2 + 0.1 * numpy.median(i2)
    iobs = HKL_data_I_sigI(hklinfo)
    iobs.set_data(hkls_i, numpy.stack([i_obs, s_obs], axis=1).astype(numpy.float32))
    pos = i_obs > 0
    fo_i = numpy.where(pos, numpy.sqrt(numpy.maximum(i_obs, 0)), 0.0)
    sf_i = numpy.where(pos, s_obs / (2 * numpy.sqrt(numpy.maximum(i_obs, 1e-30))), 1.0)
    fobs_i = HKL_data_F_sigF(hklinfo)
    fobs_i.set_data(hkls_i, numpy.stack([fo_i, sf_i], axis=1).astype(numpy.float32))

    return {
        'hklinfo': hklinfo, 'cell': cell, 'spacegroup': sg, 'grid': grid,
        'n_sym': sg.num_symops, 'n_refl': hklinfo.num_reflections,
        'elements': elements, 'coords': coords, 'occ': occ, 'is_aniso': is_aniso,
        'u_iso': u_iso, 'u_aniso': u_aniso,
        'fobs': fobs, 'phi_fom': phi_fom, 'usage': usage, 'xmap': xmap,
        'iobs': iobs, 'fobs_i': fobs_i, 'n_negative_i': int((i_obs < 0).sum()),
        'hkls': hkls_i, 'fc': Fc,
    }


def _evaluation_point(ref):
    '''A parameter set displaced from the reference — the point the gradient is
    checked at. Returns (coords, u_iso, u_aniso, occ, prime) where ``prime`` freezes
    the target scale at the reference model (avoids the occupancy / overall-B
    scale-degeneracy of a self-scaling target; see check_param_gradients).'''
    is_aniso = ref['is_aniso']
    coords = ref['coords'] + numpy.array([
        [0.12, -0.08, 0.10],
        [-0.10, 0.14, -0.06],
        [0.09, 0.07, -0.12],
        [-0.11, -0.05, 0.08],
        [0.06, 0.10, -0.09],
    ])
    u_iso = numpy.where(is_aniso == 0, ref['u_iso'] + 0.05, 0.0)
    occ = numpy.full(len(ref['elements']), 0.9)
    u_aniso = ref['u_aniso'] * 1.10
    prime = {'coords': ref['coords'], 'u_iso': ref['u_iso'],
             'u_aniso': ref['u_aniso'], 'occ': ref['occ'],
             'is_aniso': is_aniso}
    return coords, u_iso, u_aniso, occ, prime


def _map_coord_result(ref):
    from chimerax.clipper.diff.state import MapTargetState
    from chimerax.clipper.diff._checks import check_coord_gradient
    coords, _, _, _, _ = _evaluation_point(ref)
    state = MapTargetState(ref['xmap'])
    return check_coord_gradient(state, coords)


def _realspace_param_result(ref):
    from chimerax.clipper.diff.state import XrayTargetState
    from chimerax.clipper.diff._checks import check_param_gradients
    from chimerax.clipper.clipper_python import Coord_orth
    coords, u_iso, u_aniso, occ, prime = _evaluation_point(ref)
    state = XrayTargetState(ref['elements'], param_names=_ALL_PARAMS,
                            target_map=ref['xmap'],
                            target_origin=Coord_orth(0.0, 0.0, 0.0))
    return check_param_gradients(state, coords, u_iso=u_iso, u_aniso=u_aniso,
                                 occ=occ, is_aniso=ref['is_aniso'], prime=prime)


def _reciprocal_param_result(ref, kind='amplitude', direct_intensities=False):
    from chimerax.clipper.diff.state import XrayTargetState
    from chimerax.clipper.diff._checks import check_param_gradients
    coords, u_iso, u_aniso, occ, prime = _evaluation_point(ref)
    if direct_intensities:
        state = XrayTargetState(ref['elements'], param_names=_ALL_PARAMS,
                                fobs=ref['fobs_i'], phi_fom=ref['phi_fom'],
                                usage=ref['usage'], kind='intensity', iobs=ref['iobs'])
    else:
        state = XrayTargetState(ref['elements'], param_names=_ALL_PARAMS,
                                fobs=ref['fobs'], phi_fom=ref['phi_fom'],
                                usage=ref['usage'], kind=kind)
    return check_param_gradients(state, coords, u_iso=u_iso, u_aniso=u_aniso,
                                 occ=occ, is_aniso=ref['is_aniso'], prime=prime)


def _iobs_matches_derived(ref):
    '''Wiring check: iobs = (Fo^2, 2 Fo sigF) built from the reference fobs must give the
    legacy Fo-derived intensity target's value (to float32 rounding of Io; both scale
    fits then see the same, all-positive, intensities).'''
    from chimerax.clipper import HKL_data_I_sigI
    from chimerax.clipper.diff.state import XrayTargetState
    coords, u_iso, u_aniso, occ, _ = _evaluation_point(ref)
    fo = ref['fc'].astype(numpy.double)
    iobs = HKL_data_I_sigI(ref['hklinfo'])
    iobs.set_data(ref['hkls'], numpy.stack([fo ** 2, 2.0 * fo * 1.0], axis=1)
                  .astype(numpy.float32))
    common = dict(param_names=('X', 'Y', 'Z'), fobs=ref['fobs'], phi_fom=ref['phi_fom'],
                  usage=ref['usage'], kind='intensity')
    legacy = XrayTargetState(ref['elements'], **common)
    direct = XrayTargetState(ref['elements'], iobs=iobs, **common)
    args = (coords, u_iso, u_aniso, occ, ref['is_aniso'])
    l0, g0 = legacy.value_and_gradient(*args, refresh_scale=True)
    l1, g1 = direct.value_and_gradient(*args, refresh_scale=True)
    rel_l = abs(l1 - l0) / abs(l0)
    rel_g = numpy.abs(g1 - g0).max() / numpy.abs(g0).max()
    return rel_l, rel_g


# ---------------------------------------------------------------------------
# pytest-style entry points (also called by main()).
# ---------------------------------------------------------------------------

def test_map_coord_gradient():
    '''MapTargetState coordinate gradient == finite difference, P1 and P2_1.'''
    for sg, cp in _SPACEGROUPS:
        ref = _build_reference_crystal(sg, cp)
        r = _map_coord_result(ref)
        assert r['passed'], (
            '%s MapTargetState coord gradient FD check failed: '
            'rel_to_field=%.3e worst=%s' % (sg, r['rel_to_field'], r['worst_significant']))


def test_realspace_param_gradients():
    '''XrayTargetState (real-space) all 11 param gradients == FD, P1 and P2_1.'''
    for sg, cp in _SPACEGROUPS:
        ref = _build_reference_crystal(sg, cp)
        r = _realspace_param_result(ref)
        bad = {k: v for k, v in r['per_param'].items() if not v['passed']}
        assert r['passed'], '%s real-space param FD check failed: %s' % (sg, bad)


def test_reciprocal_param_gradients():
    '''XrayTargetState (reciprocal) all 11 param gradients == FD, P1 and P2_1.'''
    for sg, cp in _SPACEGROUPS:
        ref = _build_reference_crystal(sg, cp)
        r = _reciprocal_param_result(ref)
        bad = {k: v for k, v in r['per_param'].items() if not v['passed']}
        assert r['passed'], '%s reciprocal param FD check failed: %s' % (sg, bad)


def test_reciprocal_intensity_param_gradients():
    '''IntensityLS from Fo (Io = Fo^2): all 11 param gradients == FD, P1 and P2_1.'''
    for sg, cp in _SPACEGROUPS:
        ref = _build_reference_crystal(sg, cp)
        r = _reciprocal_param_result(ref, kind='intensity')
        bad = {k: v for k, v in r['per_param'].items() if not v['passed']}
        assert r['passed'], '%s intensity param FD check failed: %s' % (sg, bad)


def test_reciprocal_iobs_param_gradients():
    '''IntensityLS fed measured I/sigma(I) directly, with negative I and a varying sigma:
    all 11 param gradients == FD, P1 and P2_1.'''
    for sg, cp in _SPACEGROUPS:
        ref = _build_reference_crystal(sg, cp)
        assert ref['n_negative_i'] > 0
        r = _reciprocal_param_result(ref, direct_intensities=True)
        bad = {k: v for k, v in r['per_param'].items() if not v['passed']}
        assert r['passed'], '%s iobs param FD check failed: %s' % (sg, bad)


def test_iobs_matches_derived_intensities():
    '''iobs = (Fo^2, 2 Fo sigF) reproduces the Fo-derived intensity target.'''
    for sg, cp in _SPACEGROUPS:
        rel_l, rel_g = _iobs_matches_derived(_build_reference_crystal(sg, cp))
        assert rel_l < 1e-4 and rel_g < 1e-3, '%s: rel L %.2e, rel grad %.2e' % (sg, rel_l, rel_g)


# ---------------------------------------------------------------------------
# Headless runner (ChimeraX --script). Prints a table + OVERALL PASS/FAIL.
# ---------------------------------------------------------------------------

def _print_param_table(title, result):
    print('  %s: passed=%s' % (title, result['passed']))
    for name in _ALL_PARAMS:
        v = result['per_param'][name]
        print('     %-4s pass=%-5s rel=%.2e field=%.3g'
              % (name, v['passed'], v['rel_to_field'], v['field_scale']))


def main():
    all_ok = True
    for sg, cp in _SPACEGROUPS:
        ref = _build_reference_crystal(sg, cp)
        print('=== %s  (n_sym=%d, n_refl=%d) ==='
              % (sg, ref['n_sym'], ref['n_refl']))

        r_map = _map_coord_result(ref)
        print('  MapTargetState coord check: passed=%s rel=%.2e field=%.3g'
              % (r_map['passed'], r_map['rel_to_field'], r_map['field_scale']))
        all_ok = all_ok and r_map['passed']

        r_rs = _realspace_param_result(ref)
        _print_param_table('XrayTargetState real-space param check', r_rs)
        all_ok = all_ok and r_rs['passed']

        r_rc = _reciprocal_param_result(ref)
        _print_param_table('XrayTargetState reciprocal param check', r_rc)
        all_ok = all_ok and r_rc['passed']

        r_ri = _reciprocal_param_result(ref, kind='intensity')
        _print_param_table('XrayTargetState reciprocal intensity (Io = Fo^2) check', r_ri)
        all_ok = all_ok and r_ri['passed']

        r_io = _reciprocal_param_result(ref, direct_intensities=True)
        _print_param_table('XrayTargetState reciprocal iobs check (%d negative I)'
                           % ref['n_negative_i'], r_io)
        all_ok = all_ok and r_io['passed'] and ref['n_negative_i'] > 0

        rel_l, rel_g = _iobs_matches_derived(ref)
        ok = rel_l < 1e-4 and rel_g < 1e-3
        print('  iobs == Fo-derived intensities: passed=%s rel_L=%.2e rel_grad=%.2e'
              % (ok, rel_l, rel_g))
        all_ok = all_ok and ok
        print('')

    print('OVERALL: %s' % ('PASS' if all_ok else 'FAIL'))
    return all_ok


# ChimeraX `runscript` sets __name__ to 'ChimeraX_sandbox_N' (not '__main__');
# fire the runner for either, but stay quiet on a plain pytest import.
if __name__ == '__main__' or __name__.startswith('ChimeraX'):
    import sys
    sys.exit(0 if main() else 1)
