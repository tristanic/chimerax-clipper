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

from chimerax.core.models import Model

class SurfaceZone:
    '''
    Add this as a property to a Volume object to provide it with the
    necessary information to update its triangle mask after re-contouring.
    '''
    def __init__(self, session, distance, atoms = None, coords = None):
        '''
        Args:
          distance (float in Angstroms):
            distance from points to which the map will be masked
          atoms:
            Atoms to mask to (coordinates will be updated upon re-masking)
          coords:
            (x,y,z) coordinates to mask to (will not be updated upon
            re-masking).

          Set both atoms and coords to None to disable automatic re-masking.
        '''
        self.session = session
        self.update(distance, atoms, coords)

    def update(self, distance, atoms = None, coords = None):
        self.distance = distance
        self.atoms = atoms
        self.coords = coords


    @property
    def all_coords(self):
        if self.atoms is not None:
            import numpy
            if self.coords is not None:
                return numpy.concatenate(self.atoms.coords, self.coords)
            return self.atoms.coords
        return self.coords





    # from chimerax.surface import zone
    # for m in models:
    #     for s in m.surfaces:
    #         spoints = s.position.inverse(is_orthonormal=True) * points
    #         zone.surface_zone(s, spoints, distance, auto_update=True)


class MapMgr(Model):
    '''
    Top-level manager for all maps associated with a model.
    '''
    DEFAULT_MAX_VOXELS = 1500000

    # SESSION_SAVE=False
    def __init__(self, crystal_manager, spotlight_radius=12, default_oversampling_rate=2.0,
            auto_add = True):
        cm = self._mgr = crystal_manager
        super().__init__('Map Manager', cm.session)
        self._default_oversampling_rate=default_oversampling_rate

        self.contour_mgr = ContourResultMgr(self.session)
        self._zone_mgr = None

        if not hasattr(self, 'triggers'):
            from chimerax.core.triggerset import TriggerSet
            self.triggers = TriggerSet()

        trigger_names = (
            # Deprecated
            'map box changed',
            # Ask each MapSet to expand its volumes to cover an arbitrary set
            # of coordinates as efficiently as possible.
            'cover coords',
            # Change the radius of the "spotlight" sphere. It is up to each
            # MapSet to work out how to accommodate it
            'spotlight changed',
            'spotlight moved',    # Just changed the centre of the box
        )
        for t in trigger_names:
            self.triggers.add_trigger(t)

        self._max_voxels_for_live_remask = self.DEFAULT_MAX_VOXELS

        # Handler for live box update
        self._box_update_handler = None

        # Is the map box moving with the centre of rotation?
        self._spotlight_center = None

        self._initialize_zone_mgr()

        # Radius of the sphere in which the map will be displayed when
        # in live-scrolling mode
        self.spotlight_radius = spotlight_radius


        if self.spotlight_mode:
            self._start_spotlight_mode()

        # self.display=False
        self._rezone_pending = False
        # Apply the surface mask

        mh = self._mgr_handlers = []
        mh.append((cm, cm.triggers.add_handler('mode changed',
            self._spotlight_mode_changed_cb)))



        # self.session.triggers.add_handler('frame drawn', self._first_init_cb)
        if auto_add:
            cm.add([self])


    # def added_to_session(self, session):
    #     super().added_to_session(session)
    #     session.triggers.add_handler('frame drawn', self._first_init_cb)

    def _initialize_zone_mgr(self):
        if self._zone_mgr is None:
            from .mask_handler import ZoneMgr
            coords = [self.spotlight_center]
            self._zone_mgr = ZoneMgr(self.session, 1.5,
                5)

    @property
    def zone_mgr(self):
        return self._zone_mgr

    @property
    def xmapsets(self):
        '''
        Sets of crystallographic maps associated with this model. Each XmapSet
        handles the set of maps derived from a single crystallographic dataset.
        '''
        from .xmapset import XmapSet
        return [m for m in self.child_models() if isinstance(m, XmapSet)]

    @property
    def nxmapset(self):
        '''
        Handler for all real-space (non-crystallographic) maps associated with
        this model.
        '''
        from .nxmapset import NXmapSet
        for m in self.child_models():
            if isinstance(m, NXmapSet):
                return m
        return NXmapSet(self, 'Non-crystallographic maps')

    @property
    def all_xtal_maps(self):
        from .map_handler_base import XmapHandlerBase
        return [m for m in self.all_models() if isinstance(m, XmapHandlerBase)]

    @property
    def all_non_xtal_maps(self):
        from .nxmapset import NXmapHandler
        return [m for m in self.all_models() if isinstance(m, NXmapHandler)]

    @property
    def all_maps(self):
        from chimerax.map import Volume
        return [m for m in self.all_models() if isinstance(m, Volume)]

    @property
    def all_surfaces(self):
        return [s for m in self.all_maps for s in m.surfaces]

    @property
    def crystal_mgr(self):
        return self._mgr

    @property
    def box_center(self):
        return self.crystal_mgr.spotlight_center

    @property
    def structure(self):
        return self.crystal_mgr.structure

    # @property
    # def triggers(self):
    #     return self._triggers

    @property
    def spacegroup(self):
        return self.crystal_mgr.spacegroup

    @property
    def cell(self):
        return self.crystal_mgr.cell

    @property
    def spotlight_radius(self):
        '''Get/set the radius (in Angstroms) of the live map display sphere.'''
        return self._spotlight_radius

    @spotlight_radius.setter
    def spotlight_radius(self, radius):
        import numpy
        self._spotlight_radius = radius
        center = self.crystal_mgr.spotlight_center
        self.triggers.activate_trigger('spotlight changed',
            (center, radius)
        )
        # self._surface_zone.update(radius, coords = numpy.array([center]))
        self._zone_mgr.radius = radius
        self._zone_mgr.pad = radius
        # self._reapply_zone()

    @property
    def box_params(self):
        return (self._box_corner_xyz, self._box_corner_grid, self._box_dimensions)


    @property
    def spotlight_mode(self):
        '''
        Is live map scrolling turned on? Can only be changed via the master
        symmetry manager.
        '''
        return self.crystal_mgr.spotlight_mode

    @spotlight_mode.setter
    def spotlight_mode(self, switch):
        raise NotImplementedError(
            'Mode can only be changed via the master symmetry manager!')

    @property
    def spotlight_center(self):
        '''
        Current (x,y,z) position of the centre of the "spotlight". Read-only.
        '''
        return self.crystal_mgr.spotlight_center

    @property
    def last_covered_selection(self):
        '''
        Last set of coordinates covered by
        `self.crystal_mgr.isolate_and_cover_selection()`. Read-only.
        '''
        return self.crystal_mgr.last_covered_selection

    @property
    def display(self):
        return super().display

    @display.setter
    def display(self, switch):
        Model.display.fset(self, switch)
        if switch:
            if self.spotlight_mode:
                self._start_spotlight_mode()
            self._reapply_zone()

    def add_xmapset_from_mtz(self, mtzfile, oversampling_rate=None,
            auto_choose_reflection_data=True, auto_choose_free_flags=True):
        if oversampling_rate is None:
            oversampling_rate = self._default_oversampling_rate
        return self.add_xmapset_from_file(mtzfile, oversampling_rate,
            auto_choose_reflection_data, auto_choose_free_flags)

    def add_xmapset_from_file(self, sffile, oversampling_rate=None,
            auto_choose_reflection_data=True, auto_choose_free_flags=True):
        if oversampling_rate is None:
            oversampling_rate = self._default_oversampling_rate
        from ..clipper_mtz import ReflectionDataContainer
        mtzdata = ReflectionDataContainer(self.session, sffile,
            shannon_rate = oversampling_rate,
            auto_choose_reflection_data=auto_choose_reflection_data,
            auto_choose_free_flags=auto_choose_free_flags)
        cm = self.crystal_mgr
        if not cm.has_symmetry:
            self.session.logger.info('(CLIPPER) NOTE: No symmetry information found '
                'in model. Using symmetry from MTZ file.')
            cm = self.crystal_mgr
            cm.add_symmetry_info(mtzdata.cell, mtzdata.spacegroup, mtzdata.grid_sampling,
                mtzdata.resolution)

        elif not self.symmetry_matches(mtzdata):
            old_dim = ', '.join([f'{d:.3f}' for d in cm.cell.dim])
            new_dim = ', '.join([f'{d:.3f}' for d in mtzdata.cell.dim])
            old_angles = ', '.join([f'{d:.1f}' for d in cm.cell.angles_deg])
            new_angles = ', '.join([f'{d:.1f}' for d in mtzdata.cell.angles_deg])

            info_str = (f'Old spacegroup: {cm.spacegroup.symbol_hm}\t New:{mtzdata.spacegroup.symbol_hm}\n'
                    f'Old cell dimensions: ({old_dim})\t New: ({new_dim})\n'
                    f'Old cell angles: ({old_angles})\t New: ({new_angles})')
            if len(self.xmapsets):
                raise RuntimeError('(CLIPPER) Symmetry info from new structure factor file does not match '
                    'symmetry info from existing data!\n'+info_str)
            else:
                self.session.logger.warning('(CLIPPER) WARNING: Symmetry info from new structure factor file does '
                    'not match the symmetry info currently associated with the model. Since no other '
                    'crystal datasets are open, the new symmetry will be used.\n'+info_str
                    )
                cm.add_symmetry_info(mtzdata.cell, mtzdata.spacegroup, mtzdata.grid_sampling,
                mtzdata.resolution, override=True)
        from .xmapset import XmapSet
        return XmapSet(self, mtzdata)

    def symmetry_matches(self, xtal_data):
        return (
            xtal_data.cell.equals(self.cell, 1.0)
            and xtal_data.spacegroup.spacegroup_number == self.spacegroup.spacegroup_number
        )


    def _spotlight_mode_changed_cb(self, *_):
        if self.spotlight_mode:
            self._start_spotlight_mode()
        else:
            self._stop_spotlight_mode()

    def _start_spotlight_mode(self):
        zm = self._zone_mgr
        zm.radius = zm.pad = self.spotlight_radius
        zm.allow_remask_on_coord_updates()
        # zm.coords = [self.spotlight_center]
        self.triggers.activate_trigger('spotlight changed',
            (self.spotlight_center, self.spotlight_radius)
        )
        if self._box_update_handler is None:
            self._box_update_handler = self.crystal_mgr.triggers.add_handler(
                'spotlight moved', self.update_spotlight)
            self.update_spotlight(None, self.spotlight_center)
        from chimerax.geometry import Places
        self.positions = Places()
        self._reapply_zone()

    def _stop_spotlight_mode(self):
        if self._box_update_handler is not None:
            self.crystal_mgr.triggers.remove_handler(self._box_update_handler)
            self._box_update_handler = None

    def cover_box(self, minmax):
        '''
        Set the map box to fill a volume encompassed by the provided minimum
        and maximum xyz coordinates. Automatically turns off live scrolling.
        '''
        self.spotlight_mode = False
        xyz_min, xyz_max = minmax
        center = (xyz_min + xyz_max)/2
        self.triggers.activate_trigger('map box changed',
            (center, xyz_min, xyz_max))

    def cover_atoms(self, atoms, transforms=None, transform_indices=None,
            mask_radius=3, extra_padding=12):
        '''
        Expand all maps to a region large enough to cover the atoms, plus
        mask_radius+extra_padding in every direction. If provided, transforms
        should be a Places object, and transform_indices a Numpy array giving
        the index of the Place to be used to transform the coordinates of each
        atom. Unlike cover_coords(), the mask will be periodically updated in
        response to atom movements.
        '''
        zm = self._zone_mgr
        zm.set_symmetry_map(atoms, transforms, transform_indices)
        zm.radius = mask_radius
        zm.pad = extra_padding
        self.triggers.activate_trigger('cover coords',
            (zm.coords, mask_radius+extra_padding))
        displayed_volumes = [v for v in self.all_maps if v.display]
        num_displayed_voxels = sum(v.region_matrix(v.region).size for v in displayed_volumes)
        if num_displayed_voxels > self._max_voxels_for_live_remask:
            zm.block_remask_on_coord_updates()
        else:
            zm.allow_remask_on_coord_updates()
        # self._reapply_zone()


    def cover_coords(self, coords, mask_radius=3, extra_padding=3):
        '''
        Expand all maps to a region large enough to cover the coords, plus
        mask_radius+extra_padding in every direction.
        '''
        self.triggers.activate_trigger('cover coords',
            (coords, mask_radius+extra_padding))
        zm = self._zone_mgr
        zm.coords = coords
        zm.radius = mask_radius
        zm.pad = extra_padding
        # self._surface_zone.update(mask_radius, coords=coords)
        # self._reapply_zone()

    def update_zone_mask(self, coords):
        self._zone_mgr.coords = coords

    def update_spotlight(self, trigger_name, new_center):
        '''
        Update the position of the "spotlight" to surround the current centre of
        rotation. If this manager is not currently displayed, then filling the
        volumes with data around the new position will be deferred unless
        force is set to True.
        '''
        if self.spotlight_mode:
            import numpy
            self.triggers.activate_trigger('spotlight moved',
                new_center)
            zm = self._zone_mgr
            zm.coords = numpy.array([new_center])
            # self._surface_zone.update(self.spotlight_radius, coords = numpy.array([new_center]))
            # self._reapply_zone()

    def rezone_needed(self):
        if not self._rezone_pending:
            self._rezone_pending=True
            self.session.triggers.add_handler('new frame', self._rezone_once_cb)


    # Callbacks

    def _first_init_cb(self, *_):
        self.display = True
        self.update_spotlight(_, self.spotlight_center)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


    def _rezone_once_cb(self, *_):
        self._reapply_zone()
        self._rezone_pending=False
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _reapply_zone(self):
        '''
        Reapply any surface zone applied to the volume after changing box
        position.
        '''
        from .mask_handler import ZoneMask
        for s in self.all_surfaces:
            asm = s.auto_remask_triangles
            if asm is None:
                ZoneMask(s, self._zone_mgr, None)
            else:
                asm()

    def delete(self):
        self._stop_spotlight_mode()
        for (mgr, h) in self._mgr_handlers:
            try:
                mgr.triggers.remove_handler(h)
            except:
                continue
        self.contour_mgr.delete()
        super().delete()

    def take_snapshot(self, session, flags):
        from chimerax.core.models import Model
        data = {
            'symmetry manager': self._mgr,
            'model state': Model.take_snapshot(self, session, flags)
        }
        from .. import CLIPPER_STATE_VERSION
        data['version']=CLIPPER_STATE_VERSION
        return data

    @staticmethod
    def restore_snapshot(session, data):
        from chimerax.core.models import Model
        sh = data['symmetry manager']
        if sh is None:
            return None
        mmgr = MapMgr(sh, auto_add=False)
        Model.set_state_from_snapshot(mmgr, session, data['model state'])
        session.triggers.add_handler('end restore session', mmgr._end_restore_session_cb)
        return mmgr

    def _end_restore_session_cb(self, *_):
        self._reapply_zone()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

class ContourResultMgr:
    def __init__(self, session):
        self.session = session
        self._pending = []
        self._handler = self.session.triggers.add_handler('new frame', self._new_frame_cb)

    
    def enqueue(self, surface, rendering_options):
        self._pending.append((surface, rendering_options))
    
    def _new_frame_cb(self, *_):
        if len(self._pending):
            s, ro = self._pending.pop(0)
            if not s.deleted:
                s._use_fast_thread_result(ro)
    
    def delete(self):
        try:
            self.session.triggers.remove_handler(self._handler)
        except Exception as e:
            self.session.logger.warning(f'Failed to clean up ContourResultMgr new frame handler with the following error:\n{str(e)}')

    
