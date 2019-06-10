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

# Note that this software makes use of modified versions of the Clipper, LibCCP4
# and MMDB libraries, as well as portions of the Intel Math Kernel Library. Each
# of these is redistributed under its own license terms.

import numpy

def atom_list_from_sel(atom_list):
    '''
    Takes a ChimeraX Atoms object, and creates a Clipper Atoms_list object
    from the relevant atom properties.
    '''
    n = len(atom_list)
    elements = atom_list.element_names.tolist()
    coords = atom_list.coords
    occupancies = atom_list.occupancies
    import numpy
    from math import pi
    u_iso = numpy.sqrt(atom_list.bfactors/(8*pi**2))
    # Have to be careful here. Atoms.aniso_u6 returns None if any atom in
    # the array has no aniso_u6 entry.

    u_aniso = numpy.ones([n,6],numpy.double)*numpy.nan
    u_aniso[atom_list.has_aniso_u] = atom_list.filter(atom_list.has_aniso_u).aniso_u6
    from .clipper_python import Atom_list
    clipper_atom_list = Atom_list(elements, coords, occupancies, u_iso, u_aniso)
    return clipper_atom_list

def make_unit_cell(model, xmap, draw = True):
    from chimerax.core.geometry import Place, Places
    cell = xmap.cell()
    atoms = model.atoms
    clipper_atoms = atom_list_from_sel(atoms)
    coord = model.bounds().center()
    from .clipper_python import Coord_orth
    coord_orth = Coord_orth(coord)
    coord_frac = coord_orth.coord_frac(cell)
    sg = xmap.spacegroup()
    unit_cell_frac_symops = xmap.unit_cell_symops(coord_frac, clipper_atoms)
    if draw:
        uc_places = []
        for op in unit_cell_frac_symops.symops():
            op_orth = op.rtop_orth(cell)
            uc_places.append(Place(matrix=op_orth.matrix()[0:3,:]))
        ucp = Places(uc_places)
        model.positions = ucp
    return unit_cell_frac_symops


def draw_box(min_corner, max_corner, name='box'):
    from chimerax.core.geometry import Place, Places
    from chimerax.core.models import Drawing, Model
    from chimerax.surface.shapes import sphere_geometry
    d = Drawing('corners')
    m = Model(name, session)
    #d.vertices, d.normals, d.triangles = sphere_geometry(80)
    d.set_geometry(*sphere_geometry(80))
    minmax = [min_corner, max_corner]
    dp = []
    base_color = numpy.array([255,255,255,255])
    color_increment = numpy.array([0,-32,-32,0])
    colors = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                dp.append(Place(origin=[minmax[i][0],minmax[j][1],minmax[k][2]]))
                colors.append(base_color-color_increment*(i+j+k));
    d.positions = Places(dp)
    d.colors = colors
    m.add_drawing(d)
    session.models.add([m])
    return d

def draw_asu(xmap):
    from chimerax.core.geometry import Place, Places
    from chimerax.core.models import Drawing, Model
    from chimerax.surface.shapes import sphere_geometry
    d = Drawing('asu corners')
    m = Model('asu box', session)
    #d.vertices, d.normals, d.triangles = sphere_geometry(80)
    d.set_geometry(*sphere_geometry(80))
    asu = xmap.grid_asu()
    grid = xmap.grid_sampling()
    cell = xmap.cell()
    minmax = [asu.min().coord_frac(grid).coord_orth(cell).xyz, asu.max().coord_frac(grid).coord_orth(cell).xyz]
    dp = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                dp.append(Place(origin=[minmax[i][0],minmax[j][1],minmax[k][2]]))
    d.positions = Places(dp)
    m.add_drawing(d)
    session.models.add([m])
    return d

def box_corners(origin_xyz, size_xyz):
    ret = []
    minmax = [origin_xyz, size_xyz]
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ret.append([minmax[i][0],minmax[j][1], minmax[k][2]])
    return ret
