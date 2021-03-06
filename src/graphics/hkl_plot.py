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

from chimerax.graphics import Drawing
from chimerax.core.models import Model

def test_hkl_plot(session, datafile, scale=1.0):
    from ..clipper_mtz import ReflectionDataContainer
    rdc = ReflectionDataContainer(session, datafile)
    session.models.add([rdc])
    fsigf = rdc.experimental_data.datasets['FOBS, SIGFOBS'].data.data
    hkls = fsigf[0]
    data = fsigf[1]
    plot = HKLPlot3D('fsigf', session, hkls, data, dim_scale=scale, highlight_negatives=True)
    session.models.add([plot])
    return plot

class HKLPlot3D(Model):
    SESSION_SAVE=False
    def __init__(self, name, session, hkls, vals, scale_to=None, dim_scale=2.0, highlight_negatives=True):
        super().__init__(name, session)
        self._values = vals
        if len(vals.shape) > 1:
            # Only want the real component
            real_vals = vals[:,0].copy()
        else:
            real_vals = vals.copy()
        # Remove missing (NaN) values
        import numpy
        mask = numpy.logical_not(numpy.isnan(real_vals))
        hkls = hkls[mask]
        real_vals = real_vals[mask]
        # Normalise values
        if scale_to is None:
            real_vals /= numpy.abs(real_vals).max()
        else:
            scale_to = scale_to.copy()
            if len(scale_to.shape) > 1:
                scale_to = scale_to[:,0]
            scale_to = scale_to[numpy.logical_not(numpy.isnan(scale_to))]
            real_vals /= numpy.abs(scale_to).max()
        self._axes = _HKL_Axes(axis_length=dim_scale*(hkls.max()-hkls.min()))
        self._data_plot = _HKLPlot3D('data', hkls, real_vals, dim_scale=dim_scale, highlight_negatives=highlight_negatives)
        self.add_drawing(self._axes)
        self.add_drawing(self._data_plot)

class _HKLPlot3D(Drawing):
    def __init__(self, name, hkls, vals, dim_scale = 2.0, highlight_negatives = True):
        super().__init__(name)
        from chimerax.surface.shapes import sphere_geometry2
        self.set_geometry(*sphere_geometry2(80))
        import numpy
        from chimerax.geometry import Places, Place, identity
        id_axes = identity().axes()
        abs_vals = numpy.abs(vals)
        import numpy
        place_array = numpy.zeros((len(hkls),3,4))
        place_array[:,:,3] = hkls*dim_scale
        place_array[:,:,:3] = id_axes
        for i in range(3):
            # Volume scales with value, so radius scales with cube root of value
            place_array[:,i,i] *=abs_vals**(1/3)

        # positions =[Place(origin=hkl*dim_scale, axes=id_axes*max(rval, 0))
        #     for (hkl, rval) in zip(hkls, vals)
        # ]
        if highlight_negatives:
            import numpy
            negatives = numpy.argwhere(vals < 0).flatten()
            place_array[negatives, :, :3] = id_axes
            # for ni in negatives:
            #     positions[ni] = Place(origin=positions[ni].origin(), axes=id_axes)

        positions = Places(place_array=place_array)

        self.set_positions(Places(positions))
        from chimerax.core.colors import BuiltinColormaps
        cmap = BuiltinColormaps['rdylbu'].linear_range(-1,1)
        colors = cmap.interpolated_rgba8(vals)
        if highlight_negatives:
            colors[vals<0] = [255,0,0,255]
        self.set_colors(colors)


class _HKL_Axes(Drawing):
    # skip_bounds = True
    def __init__(self, axis_length = 20.0, axis_radius = 0.05,
                 axis_colors = [(255,0,0,255),(0,255,0,255),(0,0,255,255)]):
        # self._center = None
        Drawing.__init__(self, 'HKL Axes')
        # self.pickable = False    # Don't set depth in frontCenter mode.
        # self.no_cofr = True	# Don't include in cofr calculation.
        self._create_geometry(axis_length, axis_radius, axis_colors)

    def _create_geometry(self, axis_length, axis_radius, axis_colors):
        self.set_size(axis_length, axis_radius)
        self.set_colors(axis_colors)

    def set_size(self, axis_length, axis_radius):
        from chimerax.surface.shapes import cylinder_geometry, cone_geometry
        vaz, naz, taz = cylinder_geometry(radius = axis_radius, height = axis_length)
        vcz, ncz, tcz = cone_geometry(radius = axis_radius * 2, height = axis_length * 0.2,
                                        caps = True)
        from chimerax.geometry import Place
        vct = Place(origin = (0,0,axis_length/2))
        vcz = vct.transform_points(vcz)
        nv = len(vaz)
        tcz = tcz + nv
        from numpy import concatenate
        vaz = concatenate((vaz, vcz))
        naz = concatenate((naz, ncz))
        taz = concatenate((taz, tcz))
        nv = len(vaz)
        tx = Place(axes = [[0,0,1],[0,-1,0],[1,0,0]])
        vax, nax, tax = tx.transform_points(vaz), tx.transform_vectors(naz), taz.copy() + nv
        ty = Place(axes = [[1,0,0],[0,0,-1],[0,1,0]])
        vay, nay, tay = ty.transform_points(vaz), ty.transform_vectors(naz), taz.copy() + 2*nv

        vc = self.vertex_colors
        self.set_geometry(concatenate((vax,vay,vaz)),
                          concatenate((nax,nay,naz)),
                          concatenate((tax,tay,taz)))
        self.vertex_colors = vc		# Setting geometry clears vertex colors

    def set_colors(self, axis_colors):
        # Axis colors red = x, green = y, blue = z
        from numpy import concatenate, empty, uint8
        nv = len(self.vertices)//3
        cax, cay, caz = empty((nv,4), uint8), empty((nv,4), uint8), empty((nv,4), uint8)
        cax[:], cay[:], caz[:] = axis_colors
        self.vertex_colors = concatenate((cax,cay,caz))
