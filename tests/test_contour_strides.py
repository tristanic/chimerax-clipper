# Regression test: contouring a non-contiguous numpy view must give the same geometry
# as contouring a C-contiguous copy of the same logical array.
#
# "vol flip" hands Clipper a negative-stride view (FlipGrid.matrix() returns
# m[::-1,:,:]), and FastVolumeSurface.update_surface passes v.matrix().T straight to
# the contour engine without copying. Before contour.h's Stride was widened to
# int64_t, "Stride * Index" promoted to unsigned long on MSVC (where long is 32 bits),
# turning a negative stride into ~4.29e9 and segfaulting the session.
#
# Needs no session and no OpenGL. Run with:
#     run_chimerax.bat --nogui --exit --script tests/test_contour_strides.py

import numpy
from chimerax.clipper.contour_thread import Contour_Thread_Mgr

_TF = numpy.eye(4, dtype=numpy.float32)[:3]


def contour(data, level):
    cm = Contour_Thread_Mgr()
    cm.start_compute(data, level, 1.0, _TF, _TF, False, True)
    while not cm.ready():
        pass
    return cm.get_result()


def flip(a, axes):
    for ax in axes:
        a = {'z': a[::-1, :, :], 'y': a[:, ::-1, :], 'x': a[:, :, ::-1]}[ax]
    return a


def main():
    rng = numpy.random.default_rng(0)
    # Native ChimeraX layout is (nz, ny, nx); update_surface passes the .T of it.
    m = rng.random((40, 50, 60)).astype(numpy.float32)
    level = 0.5

    cases = [('unflipped', flip(m, '').T)]
    for axes in ('z', 'y', 'x', 'xy', 'xyz'):
        cases.append(('flip ' + axes, flip(m, axes).T))
    # matrix_slice() produces strided views for step > 1 even with no flip.
    cases.append(('step 2 + flip z', flip(m[::2, :, :], 'z').T))
    cases.append(('step 2', m[::2, :, :].T))

    for name, view in cases:
        ref = numpy.ascontiguousarray(view)      # the pre-9555d9e input
        got = contour(view, level)
        want = contour(ref, level)
        for arrays, label in zip(zip(got, want), ('vertices', 'triangles', 'normals')):
            a, b = arrays
            assert numpy.array_equal(a, b), \
                f'{name}: {label} differ (strides {view.strides}, shape {view.shape})'
        assert len(got[0]), f'{name}: contoured nothing - test is not exercising anything'
        print(f'  ok  {name:16s} strides={view.strides} '
              f'vertices={len(got[0])} triangles={len(got[1])}')

    print('test_contour_strides: PASS')


main()
