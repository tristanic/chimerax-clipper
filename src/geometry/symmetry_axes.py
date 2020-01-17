# @Author: Tristan Croll <tic20>
# @Date:   16-Jan-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 17-Jan-2020
# @License: Lesser GNU Public License version 3.0 (see LICENSE.md)
# @Copyright: 2016-2019 Tristan Croll

def unit_cell_corners(unit_cell):
    origin = unit_cell.origin.coord_frac(unit_cell.grid).coord_orth(unit_cell.cell).as_numpy()
    top_corner = unit_cell.top_corner.coord_frac(unit_cell.grid).coord_orth(unit_cell.cell).as_numpy()
    corners = [[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
    import numpy
    bounds = numpy.array([origin, top_corner])
    for point in corners:
        for i in range(3):
            point[i] = bounds[point[i], i]
    return numpy.array(corners)


def unit_cell_planes(unit_cell):
    cell_corners = unit_cell_corners(unit_cell)
    xmin = [0,2,1,3]
    xmax = [4,5,6,7]
    ymin = [0,1,4,5]
    ymax = [2,6,3,7]
    zmin = [0,4,2,6]
    zmax = [1,3,5,7]
    import numpy
    plane_points = numpy.array([cell_corners[indices] for indices in (xmin, xmax, ymin, ymax, zmin, zmax)])
    plane_normals = [numpy.cross(plane[1]-plane[0], plane[2]-plane[0]) for plane in plane_points]
    plane_normals = numpy.array([normal / numpy.linalg.norm(normal) for normal in plane_normals])
    return plane_points, plane_normals



def crystallographic_symmetry_axis_and_screw_translation(place):
    from math import degrees
    from .util import rotation_screw_axis_and_angle_affine
    constant_point, axis_direction, angle, screw_component = rotation_screw_axis_and_angle_affine(place)
    if constant_point is None:
        return (None, None, None, None)

    angle = abs(round(degrees(angle)))
    if angle==0:
        # Pure translation
        return None
    if angle not in (60, 80, 120, 180):
        err_str = ('This transform has a rotation angle of {} degrees, '
            'which is not compatible with crystallographic symmetry. '
            'In real crystals, only rotations yielding 2, 3, 4 or 6-fold '
            'symmetry are possible.').format(angle)
    fold_symmetry = 360 // angle

    return constant_point, axis_direction, screw_component, fold_symmetry

def unit_cell_and_sym_axes(session, unit_cell):
    from chimerax.core.models import Model
    m = Model('Unit cell and symmetry', session)
    from chimerax.core.geometry import Places
    syms = Places(place_array = unit_cell.symops.all_matrices_orth(unit_cell.cell, '3x4'))
    for s in syms:
        d = sym_axis_drawing(s, unit_cell)
        if d is not None:
            m.add_drawing(d)
    return m



def sym_axis_drawing(place, unit_cell):
    constant_point, axis_direction, screw_component, fold_symmetry = \
        crystallographic_symmetry_axis_and_screw_translation(place)
    if constant_point is None:
        return None

    import numpy
    if numpy.linalg.norm(screw_component) > 0.1:
        return sym_axis_drawing_screw(unit_cell, place, constant_point, axis_direction, screw_component, fold_symmetry)
    return sym_axis_drawing_standard(unit_cell, place, constant_point, axis_direction, screw_component, fold_symmetry)

def sym_axis_drawing_screw(unit_cell, place, constant_point, axis_direction, screw_component, fold_symmetry):
    print('Screw axes not yet implemented!')
    return None

def sym_axis_drawing_standard(unit_cell, place, constant_point, axis_direction, screw_component, fold_symmetry):
    plane_points, plane_normals = unit_cell_planes(unit_cell)
    from chimerax.core.geometry import ray_segment
    clip_planes = [[pp[0],pn] for pp, pn in zip(plane_points, plane_normals)]
    cell_origin = plane_points[0][0]
    import numpy
    dist = numpy.linalg.norm(constant_point-cell_origin)
    constant_point -= axis_direction*dist
    f0, f1 = ray_segment(constant_point, axis_direction, clip_planes)

    from chimerax.core.graphics import Drawing
    d = Drawing('{}-fold symmetry axis'.format(fold_symmetry))
    d.set_geometry(*_drawing_geometry(fold_symmetry))
    from chimerax.core.geometry import Places, cylinder_rotations
    rot44 = numpy.empty([1,4,4], numpy.float32)
    xyz0 = numpy.array([constant_point+axis_direction*f0])
    xyz1 = numpy.array([constant_point + axis_direction*f1])
    cylinder_rotations(xyz0, xyz1, numpy.array([0.2]), rot44)
    rot44[:,3,:3] = 0.5*(xyz0+xyz1)
    d.positions = Places(opengl_array=rot44)


    # d.position = Place(axes=place.axes(), origin = constant_point + (f0+f1/2)*axis_direction)
    return d


def _drawing_geometry(fold_symmetry):
    from chimerax.surface.shapes import cylinder_geometry
    return cylinder_geometry()
