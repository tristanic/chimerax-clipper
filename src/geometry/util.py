# @Author: Tristan Croll <tic20>
# @Date:   17-Jan-2020
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 31-Jan-2020
# @License: Lesser GNU Public License version 3.0 (see LICENSE.md)
# @Copyright: 2016-2019 Tristan Croll


def rotation_axis_and_angle_3x3(rot33):
    import numpy
    from scipy.spatial.transform import Rotation as R
    r = R.from_dcm(rot33).as_rotvec()
    from math import degrees
    rnorm = numpy.linalg.norm(r)
    angle = degrees(rnorm)
    if rnorm > 0:
        r = r/rnorm
    return (r, angle)



def rotation_screw_axis_and_angle_affine(place):
    import numpy
    rot33 = place.axes().T
    central_axis, angle = place.rotation_axis_and_angle()
    if numpy.isclose(angle, 0):
        return (None, None, None, None)

    translation = place.translation()
    screw_component = numpy.dot(translation, central_axis)*central_axis
    perpendicular_component = translation-screw_component
    A = numpy.identity(4)
    A[:3,:3] = rot33
    A[:3,3] = perpendicular_component

    eigvals, eigvecs = numpy.linalg.eig(A)

    unit_eigvals = numpy.argwhere(numpy.logical_and(
        numpy.isreal(eigvals), numpy.isclose(eigvals.real, 1))
    ).ravel()


    real_eigvecs = [eigvecs[:,i].T.real for i in unit_eigvals]
    for ev in real_eigvecs:
        if ev[3] != 0:
            const_point = (ev/ev[3])[:3]
        else:
            slope = ev[:3]

    return const_point, slope, angle, screw_component
