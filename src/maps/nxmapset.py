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
from .mapset_base import MapSetBase

class NXmapSet(MapSetBase):
    '''
    Manages real-space maps. The most important difference between these and
    crystallographic maps is that there is no guarantee that two maps will have
    the same grid (i.e. voxel size and angles).
    '''
    def __init__(self, manager, name):
        super().__init__(manager, name)
        manager.add([self])

    def add_nxmap_handler_from_file(self, filename, is_difference_map=False,
        color=None, style=None, contour=None):
        from chimerax.map_data import open_file
        grid_data = open_file(filename)[0]
        from chimerax.map.volume import volume_from_grid_data
        h = self.add_nxmap_handler_from_volume(
            volume_from_grid_data(grid_data, self.session),
            is_difference_map=is_difference_map,
            color=color, style=style, contour=contour
        )
        self._mgr.rezone_needed()
        return h


    def add_nxmap_handler_from_volume(self, volume,
        is_difference_map=False,
        color=None, style=None, contour=None):
        h = NXmapHandler(self, volume)
        if self.spotlight_mode:
            corners = _find_box_corners(self.box_center, self.display_radius,
                h.data.xyz_to_ijk_transform)
            h.expand_to_cover_coords(corners, 15)
        self.add([h])
        self.set_nxmap_display_style(h, is_difference_map=is_difference_map,
            color=color, style=style, contour=contour)
        self.crystal_mgr.normalize_scene_positions()
        self._mgr.rezone_needed()
        return h

    def set_nxmap_display_style(self, nxmap_handler, is_difference_map=False,
        color=None, style=None, contour=None):
        if style is None:
            style='mesh'
        if is_difference_map and color is not None and len(color) != 2:
            err_string = '''
            ERROR: For a difference map you need to define colours for
            both positive and negative contours, as:
            [[r,g,b,a],[r,g,b,a]] in order [positive, negative].
            '''
            raise TypeError(err_string)
        if color is None:
            if is_difference_map:
                color = self.DEFAULT_DIFF_MAP_COLORS
            elif style == 'mesh':
                color = [self.DEFAULT_MESH_MAP_COLOR]
            else:
                color = [self.DEFAULT_SOLID_MAP_COLOR]
        if contour is None:
            from ..util import guess_suitable_contour
            if is_difference_map:
                # There's no easy way to define this. For now just start at 6 sigma
                pcontour = nxmap_handler.mean_sd_rms()[1]*6
                # pcontour = guess_suitable_contour(nxmap_handler, self.structure,
                #     atom_radius_scale = 0.25)
                contour = numpy.array([-pcontour, pcontour])
            else:
                contour = numpy.array([guess_suitable_contour(nxmap_handler, self.structure)])
        nxmap_handler.set_parameters(**{'cap_faces': False,
                                  'surface_levels': contour,
                                  'show_outline_box': False,
                                  'surface_colors': color,
                                  'square_mesh': False})
        nxmap_handler.set_display_style(style)



    def expand_to_cover_coords(self, coords, padding):
        for v in self:
            v.expand_to_cover_coords(coords, padding)

    def take_snapshot(self, session, flags):
        from chimerax.core.models import Model
        data = {
            'model state':          Model.take_snapshot(self, session, flags),
            'manager':              self._mgr,
        }
        from .. import CLIPPER_STATE_VERSION
        data['version']=CLIPPER_STATE_VERSION
        return data

    @staticmethod
    def restore_snapshot(session, data):
        from chimerax.core.models import Model
        nxs = NXmapSet(data['manager'], '')
        Model.set_state_from_snapshot(nxs, session, data['model state'])
        nxs.master_map_mgr.rezone_needed()
        return nxs



from .map_handler_base import MapHandlerBase
class NXmapHandler(MapHandlerBase):
    '''
    Real-space equivalent to XmapHandler_Static. Doesn't actually use any of
    the clipper engine, but provides a unified interface.
    '''
    #SESSION_SAVE=False
    def __init__(self, mapset, volume=None, data=None, name=None,
            is_difference_map=False):
        '''
        Takes ownership of the data from an existing Volume object.
        The input volume will be closed.
        '''
        if volume is None and data is None:
            raise RuntimeError('Must provide either a Volume or a GridData '
                'instance!')
        if volume is not None and data is not None:
            raise RuntimeError('Cannot provide both volume and data')
        session = mapset.session
        if volume is not None:
            data = volume.data
            name = volume.name
            if volume in session.models.list():
                session.models.close([volume])
            else:
                volume.delete()
        else:
            if name is None:
                name = "Unnamed map"
        super().__init__(mapset, name, data,
            is_difference_map=is_difference_map)


    def _box_changed_cb(self, name, params):
        self.update_mask()

    def _box_moved_cb(self, name, params):
        self.update_mask()

    @MapHandlerBase.is_difference_map.setter
    def is_difference_map(self, flag):
        if flag != self._is_difference_map:
            self._is_difference_map = flag
            if len(self.surfaces):
                style = self.surfaces[0].display_style
            else:
                style = 'mesh'
            self.mapset.set_nxmap_display_style(self, is_difference_map=flag, style=style)
            from .mask_handler import ZoneMask
            for s in self.surfaces:
                s.display=True
                ZoneMask(s, self.mapset.master_map_mgr.zone_mgr, None)


    property
    def stats(self):
        return self.mean_sd_rms()

    @property
    def sigma(self):
        return self.mean_sd_rms()[1]

    def update_mask(self):
        if not self.display:
            return
        corners = _find_box_corners(self.center(), self.display_radius, self.data.xyz_to_ijk_transform)
        self.new_region(ijk_min=corners[0], ijk_max=corners[1], ijk_step=[1,1,1],
            adjust_step=True)

    def expand_to_cover_coords(self, coords, padding):
        self.new_region(*self.bounding_region(coords, padding=padding, step=[1,1,1]),
            adjust_step=True)


    def take_snapshot(self, session, flags):
        from chimerax.map import Volume
        from chimerax.core.models import Model
        # from chimerax.map.session import state_from_map
        data = {
            'volume state': Volume.take_snapshot(self, session, flags),
            # 'original volume': Volume.take_snapshot(self, session, flags),
            # 'model state': Model.take_snapshot(self, session, flags),
            # 'volume state': state_from_map(self),
            'is difference map': self._is_difference_map,
            'mapset': self._mapset,
        }
        from .. import CLIPPER_STATE_VERSION
        data['version']=CLIPPER_STATE_VERSION
        return data

    @staticmethod
    def restore_snapshot(session, data):
        ov_data = data['volume state']
        grid_data = ov_data['grid data state'].grid_data
        if grid_data is None:
            return None
        v = NXmapHandler(data['mapset'], data=grid_data, is_difference_map=data['is difference map'])
        from chimerax.core.models import Model
        Model.set_state_from_snapshot(v, session, ov_data['model state'])
        from chimerax.map.session import set_map_state
        set_map_state(ov_data['volume state'], v)
        v._drawings_need_update()
        return v

_corners = numpy.array([(x,y,z) for x in (-1,1) for y in (-1,1) for z in (-1,1)])
def _find_box_corners(center, radius, xyz_to_ijk_transform):
    corners = xyz_to_ijk_transform*(center+radius*_corners)
    return (numpy.min(corners, axis=0), numpy.max(corners, axis=0))
