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
from chimerax.core.models import Model

class MapSetBase(Model):
    '''
    Base class for XmapSet_Live, XmapSet_Static and NXmapSet. Provides basic
    methods for visualisation, masking etc.
    '''
    # Default contour levels and colours for generic maps. Override in the
    # derived class if you wish

    STANDARD_LOW_CONTOUR = numpy.array([1.0])
    STANDARD_HIGH_CONTOUR = numpy.array([2.5])
    STANDARD_DIFFERENCE_MAP_CONTOURS = numpy.array([-3.0, 3.0])

    DEFAULT_MESH_MAP_COLOR = [0,1.0,1.0,1.0] # Solid cyan
    DEFAULT_SOLID_MAP_COLOR = [0,1.0,1.0,0.4] # Transparent cyan
    DEFAULT_DIFF_MAP_COLORS = [[1.0,0,0,1.0],[0,1.0,0,1.0]] #Solid red and green

    def __init__(self, manager, name):
        super().__init__(name, manager.session)
        self._mgr = manager
        self._base_name = name

        if not hasattr(self, 'triggers'):
            from chimerax.core.triggerset import TriggerSet
            self.triggers = TriggerSet()

        trigger_names = (
            'map box changed',
            'map box moved',
        )
        for t in trigger_names:
            self.triggers.add_trigger(t)

        mh = self._mgr_handlers = []
        mh.append((manager,
            manager.triggers.add_handler('spotlight moved',
                self._box_moved_cb)))
        mh.append((manager,
            manager.triggers.add_handler('spotlight changed',
                self._box_changed_cb)))
        mh.append((manager,
            manager.triggers.add_handler('cover coords',
                self._cover_coords_cb)))

    # @property
    # def triggers(self):
    #     return self._triggers

    @property
    def base_name(self):
        return self._base_name
    
    @base_name.setter
    def base_name(self, name):
        self._base_name = name
        from chimerax.core.models import MODEL_NAME_CHANGED

    @property
    def master_map_mgr(self):
        return self._mgr

    @property
    def box_center(self):
        return self.master_map_mgr.box_center

    @property
    def crystal_mgr(self):
        return self.master_map_mgr.crystal_mgr

    @property
    def structure(self):
        return self.crystal_mgr.structure

    @property
    def hklinfo(self):
        return self.crystal_mgr.hklinfo

    @property
    def spacegroup(self):
        return self.crystal_mgr.spacegroup

    @property
    def cell(self):
        return self.crystal_mgr.cell

    @property
    def grid(self):
        return self.crystal_mgr.grid

    @property
    def display_radius(self):
        '''Get/set the radius (in Angstroms) of the live map display sphere.'''
        return self.master_map_mgr.spotlight_radius

    def __getitem__(self, name_or_index):
        '''Get one of the child maps by name or index.'''
        if type(name_or_index) == str:
            for m in self.child_models():
                if m.name == name_or_index:
                    return m
            raise KeyError('No map with the name "{}"!'.format(name_or_index))
        else:
            return self.child_models()[name_or_index]

    @property
    def all_maps(self):
        from chimerax.map import Volume
        return [v for v in self.child_models() if isinstance(v, Volume)]

    @property
    def spotlight_mode(self):
        return self.master_map_mgr.spotlight_mode

    @spotlight_mode.setter
    def spotlight_mode(self, switch):
        raise NotImplementedError('Spotlight mode can only be enabled/disabled '
            'via the master symmetry manager!')

    @property
    def spotlight_center(self):
        return self.master_map_mgr.spotlight_center

    @spotlight_center.setter
    def spotlight_center(self, *_):
        raise NotImplementedError('Spotlight centre can only be changed '
            'via the master symmetry manager!')


    def expand_to_cover_coords(self, coords, padding):
        raise NotImplementedError('Function not defined in the base class!')

    # Callbacks

    def _cover_coords_cb(self, trigger_name, data):
        coords, padding = data
        self.expand_to_cover_coords(coords, padding)

    def _box_changed_cb(self, trigger_name, data):
        '''
        By default, just re-fires with the same data. Override in derived class
        if more complex handling is needed
        '''
        self.triggers.activate_trigger('map box changed', data)

    def _box_moved_cb(self, trigger_name, data):
        '''
        By default, just re-fires with the same data. Override in derived class
        if more complex handling is needed
        '''
        self.triggers.activate_trigger('map box moved', data)

    def delete(self):
        for (mgr, h) in self._mgr_handlers:
            try:
                mgr.triggers.remove_handler(h)
            except:
                continue
        super().delete()
