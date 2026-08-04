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
        # Reason -> holder-count of temporary locks blocking oversampling-rate
        # changes. Held e.g. while ISOLDE runs a simulation against a live map:
        # the GPU keeps a fixed-shape density box that is refilled in place, so
        # re-gridding mid-simulation would corrupt it. See set_oversampling_rate.
        self._oversampling_locks = {}

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
    def radiation_type(self):
        '''Scattering regime ('xray'/'electron') of the crystallographic dataset
        managed here, or None if there is no crystallographic (live) xmapset.'''
        for xs in self.xmapsets:
            rt = getattr(xs, 'radiation_type', None)
            if rt is not None:
                return rt
        return None

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
            auto_choose_reflection_data=True, auto_choose_free_flags=True,
            radiation='auto'):
        if oversampling_rate is None:
            oversampling_rate = self._default_oversampling_rate
        return self.add_xmapset_from_file(mtzfile, oversampling_rate,
            auto_choose_reflection_data, auto_choose_free_flags,
            radiation=radiation)

    def add_xmapset_from_file(self, sffile, oversampling_rate=None,
            auto_choose_reflection_data=True, auto_choose_free_flags=True,
            free_flag_label=None, free_flags_file=None,
            fsigf_name=None, map_columns=None, browse=False, radiation='auto'):
        # Resolve 'auto' here, the point every entry funnels through (the `clipper
        # open` command, drag-and-drop, and core ChimeraX's `open <id>
        # structureFactors true` fetch, which calls this method directly). 'auto'
        # reads the associated model's mmCIF _exptl.method, so electron-diffraction
        # entries get electron scattering factors without any extra argument.
        if str(radiation).lower() not in ('xray', 'electron'):
            from ..cmd import _resolve_macro_radiation
            radiation = _resolve_macro_radiation(radiation, sffile, self.structure,
                logger=self.session.logger)
        if radiation == 'electron':
            self.session.logger.info('(CLIPPER) Using electron scattering factors '
                '(micro-ED / 3D-ED) for structure-factor calculation and R-factors.')
        if oversampling_rate is None:
            oversampling_rate = self._default_oversampling_rate
        from ..clipper_mtz import ReflectionDataContainer
        # When browsing, keep every experimental dataset so the user can choose
        # among them; otherwise the automatic single-dataset choice is applied at
        # load time.
        mtzdata = ReflectionDataContainer(self.session, sffile,
            shannon_rate = oversampling_rate,
            free_flag_label=free_flag_label,
            auto_choose_reflection_data=(auto_choose_reflection_data and not browse),
            auto_choose_free_flags=auto_choose_free_flags)
        # Optionally pull the free-R flags from a separate reflection file.
        if free_flags_file is not None:
            from ..clipper_mtz import adopt_free_flags_from_file
            adopt_free_flags_from_file(self.session, mtzdata, free_flags_file)
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
        # Optional interactive browser: let the user pick the experimental
        # dataset, the free-R flags, and exactly which maps to open. Strictly
        # opt-in — the default (browse=False) path is unchanged and silent.
        if browse and self.session.ui.is_gui and not self.session.in_script:
            from ..ui.mtz_browser import run_mtz_browser
            selection = run_mtz_browser(self.session, mtzdata)
            if selection is None:
                # User cancelled: fall back to fully automatic behaviour.
                pass
            else:
                fsigf_name = selection.get('fsigf_name', fsigf_name)
                map_columns = selection.get('map_columns', map_columns)
        from .xmapset import XmapSet
        return XmapSet(self, mtzdata, fsigf_name=fsigf_name, map_columns=map_columns,
            radiation=radiation)

    def symmetry_matches(self, xtal_data):
        return (
            xtal_data.cell.equals(self.cell, 1.0)
            and xtal_data.spacegroup.spacegroup_number == self.spacegroup.spacegroup_number
        )

    @property
    def oversampling_change_blocked(self):
        'True while one or more holders are blocking oversampling-rate changes.'
        return bool(self._oversampling_locks)

    @property
    def oversampling_lock_reasons(self):
        'The set of reasons currently blocking oversampling-rate changes.'
        return set(self._oversampling_locks)

    def block_oversampling_changes(self, reason='a running task'):
        '''
        Temporarily prevent set_oversampling_rate() from changing the grid. Call
        this (e.g. from ISOLDE at the start of an interactive simulation, whose GPU
        density box has a fixed shape that re-gridding would corrupt) and pair it
        with allow_oversampling_changes(reason) when done. Nesting is supported:
        the block is released only when every matching allow_oversampling_changes
        call has been made.
        For guaranteed release, prefer the oversampling_changes_blocked() context
        manager.
        '''
        self._oversampling_locks[reason] = self._oversampling_locks.get(reason, 0) + 1

    def allow_oversampling_changes(self, reason='a running task'):
        'Release a block previously taken with block_oversampling_changes(reason).'
        n = self._oversampling_locks.get(reason, 0)
        if n <= 1:
            self._oversampling_locks.pop(reason, None)
        else:
            self._oversampling_locks[reason] = n - 1

    def oversampling_changes_blocked(self, reason='a running task'):
        '''
        Context manager form of block/allow_oversampling_changes(), releasing the
        block even if the enclosed code raises::

            with map_mgr.oversampling_changes_blocked('ISOLDE simulation'):
                ... run simulation ...
        '''
        from contextlib import contextmanager
        @contextmanager
        def _cm():
            self.block_oversampling_changes(reason)
            try:
                yield
            finally:
                self.allow_oversampling_changes(reason)
        return _cm()

    @property
    def oversampling_rate(self):
        '''
        The real-space oversampling (Shannon) rate currently in force. Reads the
        rate actually applied to the reflection data (kept uniform across
        datasets by set_oversampling_rate); if this manager has no
        crystallographic reflection data, returns the default rate used for
        newly-opened datasets. Assigning to this property is equivalent to
        calling set_oversampling_rate().
        '''
        for xs in self.xmapsets:
            cd = xs.crystal_data
            if cd is not None:
                return cd.shannon_rate
        return self._default_oversampling_rate

    @oversampling_rate.setter
    def oversampling_rate(self, rate):
        self.set_oversampling_rate(rate)

    def set_oversampling_rate(self, rate):
        '''
        Change the oversampling (Shannon) rate for all crystallographic maps
        managed here, at runtime. This rebuilds the map grid, re-grids every live
        and static map (re-FFT from cached coefficients - no Fcalc regeneration,
        the reflection data is untouched), rebuilds the symmetry Unit_Cell to match
        the new grid, and refreshes the display. A higher rate gives a finer
        real-space grid (smoother contours) at the cost of memory and FFT time.

        Refuses if a temporary block is in force (see block_oversampling_changes),
        e.g. while an interactive simulation is running against a live map.
        '''
        from chimerax.core.errors import UserError
        if self.oversampling_change_blocked:
            reasons = ', '.join(sorted(self._oversampling_locks)) or 'a running task'
            raise UserError('The oversampling rate cannot be changed right now '
                'because it is temporarily locked ({}). This typically happens '
                'while an interactive simulation is running against a live map '
                '(changing the grid would corrupt the fixed-size density box held '
                'on the GPU). Try again once it finishes.'.format(reasons))
        rate = float(rate)
        if rate <= 0:
            raise UserError('Oversampling (Shannon) rate must be positive!')
        new_grid = None
        n_regridded = 0
        for xs in self.xmapsets:
            cd = xs.crystal_data
            if cd is None:
                # e.g. a small-molecule or volume-only xmapset: no reflection grid
                continue
            cd.set_shannon_rate(rate)
            new_grid = cd.grid_sampling
            xm = xs.live_xmap_mgr
            if xm is not None and hasattr(xm, 'set_grid_sampling'):
                xm.set_grid_sampling(new_grid)
            for h in xs.static_xmaps:
                h._regrid(new_grid)
            n_regridded += 1
        if new_grid is None:
            from chimerax.core.errors import UserError
            raise UserError('No crystallographic reflection data found to re-grid.')
        # Push the same grid into the symmetry manager and rebuild the Unit_Cell,
        # whose integer symmetry operators are tied to the grid sampling.
        cm = self.crystal_mgr
        cm.add_symmetry_info(cm.cell, cm.spacegroup, new_grid, cm.resolution,
            override=True)
        # New datasets opened later on this manager should use the new rate too.
        self._default_oversampling_rate = rate
        self._refresh_after_grid_change()
        self.session.logger.info(
            '(CLIPPER) Set oversampling rate to {:.2f} for {} crystallographic '
            'map set(s); new grid {}.'.format(rate, n_regridded, tuple(new_grid.dim)))

    def _refresh_after_grid_change(self):
        '''
        Force display geometry (box dimensions, which depend on grid spacing) to be
        recomputed after a grid change, and re-mask surfaces. The 'spotlight
        changed' trigger is the same one the spotlight-start path uses to build
        boxes from scratch.
        '''
        if self.spotlight_mode:
            self.triggers.activate_trigger('spotlight changed',
                (self.spotlight_center, self.spotlight_radius))
        self.rezone_needed()


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

        When transforms and transform_indices are omitted (or None), a single
        identity transform is used so that all atoms are treated as-is.
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
    # Per-frame wall-clock budget (seconds) for applying completed contour results
    # to the GUI/GPU. Replaces the old "exactly one surface per frame" throttle,
    # whose refresh latency scaled with the number of maps (N maps x M surfaces
    # took N*M frames to fully update). We still bound the per-frame cost so a
    # burst of completed surfaces can't stall the frame, but we drain as many as
    # fit in the budget (always at least one, to guarantee forward progress).
    frame_budget = 0.008

    def __init__(self, session):
        self.session = session
        self._pending = []
        self._handler = self.session.triggers.add_handler('new frame', self._new_frame_cb)


    def enqueue(self, surface, rendering_options):
        self._pending.append((surface, rendering_options))

    def _new_frame_cb(self, *_):
        if not self._pending:
            return
        from time import perf_counter
        start = perf_counter()
        # Results are only enqueued once their contour thread is ready (see the
        # delayed_reaction in FastVolumeSurface.update_surface), so applying one
        # never blocks on get_result(). Drain within the frame budget; always
        # apply at least one so progress is guaranteed.
        while self._pending:
            s, ro = self._pending.pop(0)
            if not s.deleted:
                s._use_fast_thread_result(ro)
            if perf_counter() - start >= self.frame_budget:
                break
    
    def delete(self):
        try:
            self.session.triggers.remove_handler(self._handler)
        except Exception as e:
            self.session.logger.warning(f'Failed to clean up ContourResultMgr new frame handler with the following error:\n{str(e)}')

    
