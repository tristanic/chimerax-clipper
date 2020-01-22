# @Author: Tristan Croll <tic20>
# @Date:   16-Jan-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 21-Jan-2020
# @License: Lesser GNU Public License version 3.0 (see LICENSE.md)
# @Copyright: 2016-2019 Tristan Croll

def unit_cell_corners(unit_cell, fractional_coords=False, pad=0):
    origin_grid = unit_cell.origin
    origin_frac = origin_grid.coord_frac(unit_cell.grid)
    from chimerax.clipper.clipper_python import Coord_frac
    top_frac = origin_frac+Coord_frac([1,1,1])
    origin = origin_frac.coord_orth(unit_cell.cell).as_numpy()
    top_corner = top_frac.coord_orth(unit_cell.cell).as_numpy()
    corners = [[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
    import numpy
    bounds = numpy.array([origin, top_corner])
    for i, point in enumerate(corners):
        if pad:
            point = numpy.array(point, numpy.float32)
            point[point==0]-=pad
            point[point==1]+=pad
        coord = (Coord_frac(point) + origin_frac)
        if not fractional_coords:
            coord = coord.coord_orth(unit_cell.cell)
        corners[i] = coord.as_numpy()
    return numpy.array(corners)

def unit_cell_box_drawing(unit_cell_corners, corner_radius=1.0, origin_radius=2.0, cylinder_radius=0.3):
    from chimerax.core.graphics import Drawing
    md = Drawing('Unit cell box')
    cd = Drawing('Corners')
    ed = Drawing('Edges')
    md.add_drawing(cd)
    md.add_drawing(ed)
    from chimerax.surface.shapes import sphere_geometry2, cylinder_geometry
    cd.set_geometry(*sphere_geometry2(80))
    ed.set_geometry(*cylinder_geometry())
    from chimerax.core.geometry import Places, cylinder_rotations
    import numpy
    shift_and_scale = numpy.empty((8,4), numpy.float32)
    shift_and_scale[:,:3] = unit_cell_corners
    shift_and_scale[1:,3] = corner_radius
    shift_and_scale[0,3] = origin_radius
    cd.positions=Places(shift_and_scale=shift_and_scale)
    edges = [[0,1],[0,2],[0,4],[1,3],[1,5],[2,3],[2,6],[3,7],[4,5],[4,6],[5,7],[6,7]]
    xyz=[[],[]]
    for e in edges:
        for i in range(2):
            xyz[i].append(unit_cell_corners[e[i]])
    xyz0, xyz1 = [numpy.array(x) for x in xyz]
    rot44 = numpy.empty((12,4,4), numpy.float32)
    cylinder_rotations(xyz0, xyz1, numpy.ones(12, numpy.float32)*cylinder_radius, rot44)
    rot44[:,3,:3] = 0.5*(xyz0+xyz1)
    ed.positions = Places(opengl_array=rot44)
    return md

def unit_cell_planes(unit_cell, fractional_coords=False, pad=0):
    cell_corners = unit_cell_corners(unit_cell, fractional_coords=fractional_coords, pad=pad)
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

def plane_intersections(origin, direction, plane_points, plane_normals):
    import numpy
    intersections = []
    for pp, pn in zip(plane_points, plane_normals):
        dn = numpy.dot(direction, pn)
        if dn != 0:
            pon = numpy.dot(origin-pp, pn)
            f = -pon/dn
            intersections.append((origin+f*direction))
    return intersections

def plane_intersections_in_unit_cell(axis_defs, plane_points, plane_normals):
    '''
    *** REQUIRES FRACTIONAL COORDINATES ***

    Finds all unique segments of the given axes passing through the unit cell.
    Returns a paired list of coordinates representing the start and end point of
    each axis segment.
    '''
    unique_intersections = set()
    import numpy
    from fractions import Fraction
    for (origin, direction) in axis_defs:
        #print('Current direction: {}'.format(numpy.round(direction, 4)))
        for u in [0]:# (-1, 0, 1):
            for v in [0]: #(-1, 0, 1):
                for w in [0]: #(-1, 0, 1):
                    current_origin = origin + [u,v,w]
                    #print('Current origin: {}'.format(numpy.round(current_origin, 4)))
                    intersections = plane_intersections(current_origin, direction, plane_points, plane_normals)
                    #print('Intersections: {}'.format([numpy.round(i, 4) for i in intersections]))
                    f_intersections = [numpy.array([float(Fraction(c).limit_denominator(24)) for c in coord]) for coord in intersections]
                    #print('Massaged intersections: {}'.format([numpy.round(f, 4) for f in f_intersections]))
                    in_box = []
                    for intersection in f_intersections:
                        if numpy.all(numpy.logical_and(intersection>=0, intersection<=1)):
                            in_box.append(intersection)
                    if len(in_box) > 0:
                        uis = frozenset([tuple(i) for i in in_box])
                        if len(uis) > 2:
                            print('WARNING: Impossible number of intersections detected for the one line! Values are {}; original: {}'.format(in_box, intersections))
                            print('Current direction: {}'.format(numpy.round(direction, 4)))
                            print('Current origin: {}'.format(numpy.round(current_origin, 4)))
                            continue
                        if len(uis) == 2:
                            unique_intersections.add(uis)
    if len(unique_intersections):
        uvw0, uvw1 = numpy.array([numpy.array([numpy.array(coord) for coord in s]) for s in unique_intersections]).transpose([1,0,2])
        #print('\n'.join([','.join([str(u0), str(u1)]) for (u0, u1) in zip(uvw0, uvw1)]))
    else:
        print('No intersections found!')
        return None, None
    return uvw0, uvw1





def crystallographic_symmetry_axis_and_screw_translation(place, unit_cell, fractional_coords=False):
    from math import degrees
    from .util import rotation_screw_axis_and_angle_affine
    constant_point, axis_direction, angle, screw_component = rotation_screw_axis_and_angle_affine(place)
    if constant_point is None:
        return (None, None, None, None)

    angle = abs(round(angle))
    if angle==0:
        # Pure translation
        return None
    if angle not in (60, 80, 120, 180):
        err_str = ('This transform has a rotation angle of {} degrees, '
            'which is not compatible with crystallographic symmetry. '
            'In real crystals, only rotations yielding 2, 3, 4 or 6-fold '
            'symmetry are possible.').format(angle)
        raise RuntimeError(err_str)
    fold_symmetry = 360 // angle

    if fractional_coords:
        from chimerax.clipper.clipper_python import Coord_orth
        cell = unit_cell.cell
        constant_point = Coord_orth(constant_point).coord_frac(cell).as_numpy()
        #print('Direction in orthographic coords: {}'.format(axis_direction))
        axis_direction = Coord_orth(axis_direction*100).coord_frac(cell).as_numpy()
        #print('Direction in fractional coords: {}'.format(axis_direction))
        screw_component = Coord_orth(screw_component).coord_frac(cell).as_numpy()

    return constant_point, axis_direction, screw_component, fold_symmetry

def unit_cell_and_sym_axes(session, unit_cell):
    from chimerax.core.models import Model
    m = Model('Unit cell and symmetry', session)
    from chimerax.core.geometry import Places
    from collections import defaultdict
    import numpy
    axis_defs = defaultdict(lambda: list())

    from chimerax.clipper.clipper_python import Coord_frac
    box_origin_frac = unit_cell.origin.coord_frac(unit_cell.grid)-Coord_frac([1,1,1])
    box_origin_xyz = box_origin_frac.coord_orth(unit_cell.cell).as_numpy()
    box_size = unit_cell.grid.dim*3

    syms = Places(place_array = unit_cell.all_symops_in_box(box_origin_xyz, box_size).all_matrices_orth(unit_cell.cell, '3x4'))

    #syms = Places(place_array = unit_cell.symops.all_matrices_orth(unit_cell.cell, '3x4'))

    for s in syms:
        constant_point, axis_direction, screw_component, fold_symmetry = \
            crystallographic_symmetry_axis_and_screw_translation(s, unit_cell, fractional_coords=True)
        if constant_point is None:
            continue
        from fractions import Fraction
        screw_component = Fraction(numpy.linalg.norm(screw_component)).limit_denominator(12) # LCM of 1,2,3,4 and 6
        axis_defs[(fold_symmetry, float(screw_component))].append([constant_point, axis_direction])

    plane_points, plane_normals = unit_cell_planes(unit_cell, fractional_coords=True, pad=0.001)
    plane_points = [pp[0] for pp in plane_points]

    for (fold_symmetry, screw_component), axes in axis_defs.items():
        uvw0, uvw1 = plane_intersections_in_unit_cell(axes, plane_points, plane_normals)
        if uvw0 is None:
            continue
        from chimerax.clipper.clipper_python import Coord_frac
        axyz0 = numpy.array([Coord_frac(uvw).coord_orth(unit_cell.cell).as_numpy() for uvw in uvw0])
        axyz1 = numpy.array([Coord_frac(uvw).coord_orth(unit_cell.cell).as_numpy() for uvw in uvw1])
        d = sym_axis_drawing(fold_symmetry, screw_component, axyz0, axyz1)
        if d is not None:
            m.add_drawing(d)

    bd = unit_cell_box_drawing(unit_cell_corners(unit_cell))
    m.add_drawing(bd)

    return m


def sym_axis_drawing(fold_symmetry, screw_component, axyz0, axyz1):
    if screw_component > 0:
        return sym_axis_drawing_screw(fold_symmetry, screw_component, axyz0, axyz1)
    return sym_axis_drawing_standard(fold_symmetry, axyz0, axyz1)

def sym_axis_drawing_screw(fold_symmetry, screw_component, axyz0, axyz1):
    print('Screw axes not yet implemented!')
    return None

_symmetry_colors = {
    2:  [[255,255,255,255],[255, 128, 128, 255]],
    3:  [[255,255,255,255],[128, 255, 255, 255],[255,128,255,255]],
    4:  [[255,255,255,255],[128, 255, 128, 255],[255,255,255,255],[128,255,128,255]],
    6:  [[255,255,255,255],[128,128,255,255],[255, 255, 128, 255],[255,255,255,255],[128,128,255,255],[255, 255, 128, 255]],
}
def sym_axis_drawing_standard(fold_symmetry, axyz0, axyz1, radius=0.3):
    import numpy
    n = len(axyz0)
    radii = numpy.ones(n, numpy.float32)*radius
    from chimerax.core.geometry import Places, cylinder_rotations, rotation
    rot44 = numpy.empty([n,4,4], numpy.float32)
    cylinder_rotations(axyz0, axyz1, radii, rot44)
    rot44[:,3,:3] = 0.5*(axyz0+axyz1)


    # d.clip_cap=True
    vertices = []
    normals = []
    triangles = []
    colors = []
    angle = 360/fold_symmetry
    v, n, t = cylindrical_wedge_geometry(angle=angle, nc=3)
    colors_base = numpy.ones((len(v), 4), numpy.uint8)
    # colors_1 = numpy.array([_symmetry_colors[fold_symmetry]]*len(v), numpy.uint8)
    # colors_2 = numpy.ones(colors_1.shape, numpy.uint8) * 255
    for i in range(fold_symmetry):
        colors_base[:,:] = _symmetry_colors[fold_symmetry][i]
        r = rotation((0,0,1), angle*i)
        vertices.append(r*v)
        normals.append(r.apply_without_translation(n))
        triangles.append(t + len(v)*i)
        colors.append(colors_base.copy())
    vertices = numpy.concatenate(vertices)
    normals = numpy.concatenate(normals)
    triangles = numpy.concatenate(triangles)
    colors = numpy.concatenate(colors)



    from chimerax.core.graphics import Drawing
    d = Drawing('{}-fold symmetry axis'.format(fold_symmetry))
    d.set_geometry(vertices, normals, triangles)
    d.set_vertex_colors(colors)
    #d.clip_cap = True

    d.positions = Places(opengl_array=rot44)

    # d.position = Place(axes=place.axes(), origin = constant_point + (f0+f1/2)*axis_direction)
    return d




def _drawing_geometry(fold_symmetry):
    from chimerax.surface.shapes import cylinder_geometry
    return cylinder_geometry()


def cylindrical_wedge_geometry(radius=1, height=1, angle=60, nz=2, nc=10, caps=True):
    '''
    Return vertex, normal vector and triangle arrays for cylinder geometry
    with specified radius and height centered at the origin.
    '''
    varray, narray, tarray = unit_cylindrical_wedge_geometry(nz, nc, angle)
    varray[:,0] *= radius
    varray[:,1] *= radius
    varray[:,2] *= height

    if not caps:
        return varray, narray, tarray

    nc += 1
    # Duplicate end rings to make sharp crease at cap edge.
    vc = varray.shape[0]
    varray.resize((vc+2*nc+2,3))
    narray.resize((vc+2*nc+2,3))
    varray[vc:vc+nc,:] = varray[0:nc,:] # Copy circle
    varray[vc+nc,:] = (0,0,-0.5*height) # Center of circle
    varray[vc+nc+1:vc+2*nc+1,:] = varray[vc-nc:vc,:] # Copy circle
    varray[vc+2*nc+1,:] = (0,0,0.5*height) # Center of circle
    narray[vc:vc+nc+1,:] = (0,0,-1)
    narray[vc+nc+1:,:] = (0,0,1)

    tc = tarray.shape[0]
    tarray.resize((tc+2*nc,3))
    for i in range(nc):
        tarray[tc+i,:] = (vc+nc,vc+(i+1)%nc,vc+i)
        tarray[tc+nc+i,:] = (vc+2*nc+1,vc+nc+1+i,vc+nc+1+(i+1)%nc)

    return varray, narray, tarray


def unit_cylindrical_wedge_geometry(nz, nc, angle=60):
    nc += 1 # One extra face to form point of wedge
    from numpy import empty, float32, arange, cos, sin, int32, pi, radians
    vc = nz*nc
    tc = (nz-1)*nc*2
    varray = empty((vc, 3), float32)
    narray = empty((vc, 3), float32)
    tarray = empty((tc, 3), int32)

    # Calculate vertices

    v = varray.reshape((nz, nc, 3))
    angles = radians(angle)/(nc-2)*arange(nc-1)

    v[:,:-1,0] = cos(angles)
    v[:,:-1,1] = sin(angles)
    v[:,-1,:2] = 0
    for z in range(nz):
        v[z,:,2] = z/(nz-1) - 0.5

    # Set normals
    narray[:,:] = varray
    narray[:,2] = 0

    # TODO: normals are wrong for wedge faces - not a huge issue at present
    # since these will normally be hidden anyway. Still, would be best to
    # correct them.

    # Create triangles
    t = tarray.reshape((nz-1, nc, 6))
    c = arange(nc)
    c1 = (c+1)%nc
    t[:,:,0] = t[1::2,:,3] = c
    t[::2,:,1] = t[::2,:,3] = t[1::2,:,1] = c1
    t[::2,:,4] = t[1::2,:,2] = t[1::2,:,4] = c1+nc
    t[::2,:,2] = t[:,:,5] = c+nc
    for z in range(1,nz-1):
        t[z,:,:] += z*nc

    return varray, narray, tarray
