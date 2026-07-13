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
End-to-end headless test of the small-molecule COD -> structure-factors ->
differentiable-target pipeline via the supported one-call API
:func:`chimerax.clipper.diff.crystal.small_molecule_ensemble_target`.

Runs on the two bundled demo entries and asserts a FINITE, non-trivial gradient —
which specifically guards the three corecif corrections that silently break a
hand-assembled pipeline:

  * ``cod_2213867`` — organic, P b c a (orthorhombic, general positions): the
    straightforward case.
  * ``cod_1100908`` — a Cu(II) complex in C 2/c (monoclinic, metal on a 2-fold
    special position): exercises BOTH the oblique-cell coordinate fix (skip it and
    the monoclinic structure factors are silently wrong) AND special-position
    fragment completion (skip it and expansion leaves partial molecules -> NaN).

A B-factor of 0 (un-hydrated model) would give infinitely sharp atoms and a NaN
gradient with a deceptively finite loss, so ``numpy.isfinite(grad).all()`` is the
load-bearing assertion.

How to run (headless, against an install carrying src/diff + the demo data)::

    "C:/Program Files/ChimeraX-Daily/bin/ChimeraX-console.exe" \
        --nogui --exit --script tests/test_small_molecule_ensemble.py
'''

import os
import numpy


def _demo_dir():
    import chimerax.clipper
    return os.path.join(os.path.dirname(chimerax.clipper.__file__), 'demo')


def _check_entry(session, entry):
    '''Build the target for one demo entry and return the gradient result dict.'''
    from chimerax.core.commands import run
    from chimerax.clipper.diff.crystal import small_molecule_ensemble_target
    cif = os.path.join(_demo_dir(), entry + '.cif')
    state, box = small_molecule_ensemble_target(session, cif, param_names=('X', 'Y', 'Z'))
    coords = box.atoms.coords
    n_atoms = box.num_atoms                       # capture BEFORE close (box is deleted)
    state.refresh_scaling(coords)
    L, grad = state.value_and_gradient(coords)
    result = {'entry': entry, 'n_atoms': n_atoms, 'L': float(L),
              'finite': bool(numpy.isfinite(grad).all()),
              'grad_max': float(numpy.abs(grad).max()), 'shape': grad.shape}
    run(session, 'close session')
    return result


def _assert_result(r):
    assert r['finite'], '%s: gradient not finite (un-hydrated B=0? partial molecule?)' % r['entry']
    assert numpy.isfinite(r['L']), '%s: loss not finite' % r['entry']
    assert r['grad_max'] > 0.0, '%s: gradient is all-zero (no signal)' % r['entry']
    assert r['shape'][0] == r['n_atoms'] and r['shape'][1] == 3


def main(session):
    ok = True
    for entry in ('cod_2213867', 'cod_1100908'):
        r = _check_entry(session, entry)
        try:
            _assert_result(r)
            status = 'PASS'
        except AssertionError as e:
            status = 'FAIL (%s)' % e
            ok = False
        print('%-14s n_atoms=%3d L=%.6g finite=%s |grad|max=%.4g -> %s'
              % (r['entry'], r['n_atoms'], r['L'], r['finite'], r['grad_max'], status))
    print('OVERALL: %s' % ('PASS' if ok else 'FAIL'))
    return ok


# ChimeraX `runscript` sets __name__ to 'ChimeraX_sandbox_N' (not '__main__') and
# provides `session` as a global; fire the runner for either.
if __name__ == '__main__' or __name__.startswith('ChimeraX'):
    import sys
    sys.exit(0 if main(session) else 1)          # noqa: F821  (session is a script global)
