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

from chimerax.map import Volume

class ZoneMgr:
    report_timing=False
    def __init__(self, session, grid_step, radius,
            atoms=None, transforms=None, transform_indices=None, coords=None, pad=None, interpolation_threshold=0.75):
        if pad is None:
            pad = radius
        self.structure = None
        self._symmetry_map = {}
        self._atoms = atoms
        if atoms is not None:
            if len(atoms) == 0:
                from chimerax.core.errors import UserError
                raise UserError('Attempting to apply a zone mask, but no atoms selected!')
            self.structure = self._unique_structure(atoms)
        self._structure_change_handler = None
        self.session = session
        self._step = grid_step
        self._radius = radius
        self._coords = coords
        self._pad = pad
        self.threshold=interpolation_threshold
        self._mask = None
        self._update_needed = False
        self._resize_box = True
        from chimerax.core.triggerset import TriggerSet
        self.triggers = TriggerSet()
        self.triggers.add_trigger('atom coords updated')

    @property
    def coords(self):
        if self._atoms is not None:
            if not len(self._symmetry_map):
                return self._atoms.coords
            else:
                import numpy
                coords = numpy.concatenate([tf*(atoms.coords) for atoms, tf in self._symmetry_map.values()])
                # print('Transformed coords: {}'.format(coords)')
                return coords
        return self._coords

    def block_remask_on_coord_updates(self):
        if not self.triggers.is_trigger_blocked('atom coords updated'):
            self.triggers.manual_block('atom coords updated')

    def allow_remask_on_coord_updates(self):
        if self.triggers.is_trigger_blocked('atom coords updated'):
            self.triggers.manual_release('atom coords updated')

    @coords.setter
    def coords(self, coords):
        self._symmetry_map.clear()
        self.stop_tracking_changes()
        if self.coords is not None:
            prev_coord_len = len(self.coords)
        else:
            prev_coord_len = None
        self._atoms = None
        self.structure = None
        self._coords = coords
        if len(coords) == 1 and prev_coord_len==1:
            self.update_needed(resize_box=False)
        else:
            self.update_needed(resize_box=True)

        # Clear all handlers
        self.triggers.delete_trigger('atom coords updated')
        self.triggers.add_trigger('atom coords updated')

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        self._coords = None
        self._atoms = atoms
        self._symmetry_map.clear()
        self.stop_tracking_changes()
        self.structure = self._unique_structure(atoms)
        self.start_tracking_changes()
        self._transforms = None
        self._transform_indices = None
        self.start_tracking_changes()
        self.update_needed(resize_box=True)

    @property
    def symmetry_map(self):
        return self._symmetry_map

    def set_symmetry_map(self, atoms, transforms, transform_indices):
        self._coords=None
        self._atoms = atoms
        import numpy
        unique_indices = numpy.unique(transform_indices)
        self._symmetry_map.clear()
        for i in unique_indices:
            mask = (transform_indices==i)
            self._symmetry_map[i] = (atoms[mask],transforms[i])
        self.stop_tracking_changes()
        self.structure = self._unique_structure(atoms)
        self.start_tracking_changes()
        self.update_needed(resize_box=True)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, radius):
        self._radius = radius
        self.update_needed(resize_box=True)

    @property
    def grid_step(self):
        return self._step

    @grid_step.setter
    def grid_step(self, step):
        self._step = step
        self.update_needed(resize_box=True)

    @property
    def pad(self):
        return self._pad

    @pad.setter
    def pad(self, pad):
        self._pad = pad
        self.update_needed(resize_box=True)

    def update_needed(self, resize_box=False):
        self._update_needed = True
        if not self._resize_box:
            self._resize_box = resize_box


    @property
    def mask(self):
        if self._mask is None:
            coords = self.coords
            if coords is None:
                raise TypeError('Must provide either atoms or coordinates!')
            self._mask = VolumeMask(self.session, coords, self.grid_step, self.radius, pad=self.pad)
        return self._mask

    def _update_mask(self):
        if self.report_timing:    
            from time import perf_counter
            start_time = perf_counter()
        update_origin = (len(self.coords) == 1)
        self.mask.generate_mask(self.coords, self.radius,
            reuse_existing=not self._resize_box, update_origin=update_origin,
            step=self.grid_step, pad=self.pad)
        self._update_needed = False
        self._resize_box = False
        if self.report_timing:
            self.session.logger.info(f'Updating mask took {(perf_counter()-start_time)*1e3} ms.')

    def get_vertex_mask(self, vertices):
        if self._update_needed:
            self._update_mask()
        return (self.mask.interpolated_values(vertices) >= self.threshold)


    def _unique_structure(self, atoms):
        us = atoms.unique_structures
        from chimerax.core.errors import UserError
        if len(us) == 0:
            raise UserError('Tried to define a zone mask with no atoms selected!')
        if len(us) != 1:
            raise UserError('All atoms for zone mask must be from a single model!')
        return us[0]

    def start_tracking_changes(self):
        if self.structure is None or self._atoms is None:
            raise RuntimeError('Live zone mask updating is only valid for atoms!')
        self._structure_change_handler = self.structure.triggers.add_handler(
            'changes', self._model_changes_cb
        )

    def stop_tracking_changes(self):
        if self._structure_change_handler is not None:
            if self.structure is not None:
                self.structure.triggers.remove_handler(self._structure_change_handler)
            self._structure_change_handler = None

    def _model_changes_cb(self, trigger_name, changes):
        if self._atoms is None:
            self._structure_change_handler = None
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER
        if 'coord changed' in changes[1].atom_reasons():
            self.update_needed()
            self.triggers.activate_trigger('atom coords updated', None)

class VolumeMask(Volume):
    '''
    Manages a binary mask defining displayed regions for other surfaces.
    '''
    SESSION_SAVE=False
    def __init__(self, session, coords, step, radius, pad=0):
        darray = self._generate_data_array(coords, step, radius, pad=pad)
        super().__init__(session, darray)
        self.name = "mask"
        self.generate_mask(coords, radius, clear=False, pad=pad)

    def _generate_data_array(self, coords, step, radius, pad=0):
        import numpy
        if hasattr(step, '__len__'):
            step = min(step)

        step = max((step, radius/4))
        step = min(step, 1.5)
        step = numpy.ones(3) * step
        cmin = coords.min(axis=0)
        cmax = coords.max(axis=0)
        mincoor = cmin-radius-pad
        maxcoor = cmax+radius+pad
        dim = numpy.ceil((maxcoor-mincoor)/step).astype(int)
        data = numpy.zeros(dim, numpy.uint8)
        self._data_fill_target = data
        from chimerax.map_data import ArrayGridData
        darray = ArrayGridData(data.transpose(), origin=mincoor, step=step)
        return darray

    def generate_mask(self, coords, radius, clear=True, reuse_existing=True, update_origin=False, step=None, pad=0):
        from chimerax.clipper import _map_mask
        if step is None:
            step = self.data.step
        if not reuse_existing:
            clear=False
            self.replace_data(self._generate_data_array(coords, step, radius, pad=pad))
        elif update_origin:
            self.data.set_origin(coords.min(axis=0)-radius-pad)
        if clear:
            self.clear()
        origin, step = self.data_origin_and_step()
        _map_mask.generate_mask(self._data_fill_target, origin, step,
            self._data_fill_target.shape, self.data.ijk_to_xyz_transform.matrix,
            self.data.xyz_to_ijk_transform.matrix, coords, len(coords), radius)
        self.data.values_changed()
        self._update_needed = False

    def clear(self):
        self.data.array[:] = 0
        self.data.values_changed()

    def interpolated_gradients(self, *args, **kwargs):
        raise NotImplementedError('Gradients are not available for a binary mask!')



from chimerax.core.state import State

class ZoneMask(State):
    report_timing=False
    SESSION_SAVE=False
    def __init__(self, surface, mgr, max_components, time_per_remask = 0.5):
        self.surface = surface
        self.mgr = mgr
        self.max_components = max_components
        self._time_per_remask = time_per_remask
        from time import time
        self._last_remask_time = time()
        self.set_surface_mask()

    def __call__(self):
        self.set_surface_mask()

    def _coord_update_cb(self, *_):
        from time import time
        current_time = time()
        if current_time - self._last_remask_time > self._time_per_remask:
            self.set_surface_mask()
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER



    def set_surface_mask(self):
        if self.report_timing:
            from time import perf_counter
            start_time = perf_counter()
        surface = self.surface
        v, t = surface.vertices, surface.triangles
        if t is None:
            surface.auto_remask_triangles = self
            return

        import numpy
        mask = self.mgr.get_vertex_mask(v)
        tmask = numpy.logical_and(mask[t[:,0]], mask[t[:,1]])
        numpy.logical_and(tmask, mask[t[:,2]], tmask)
        surface.triangle_mask = tmask

        if self.max_components is not None:
            from chimerax.surface import dust
            dust.show_only_largest_blobs(surface, True, self.max_components)
        surface.auto_remask_triangles = self
        if self.mgr.atoms is not None:
            self.mgr.triggers.add_handler('atom coords updated', self._coord_update_cb)
        from time import time
        self._last_remask_time = time()
        if self.report_timing:
            self.surface.session.logger.info(f'Applying zone mask to surface #{self.surface.id_string} took {(perf_counter()-start_time)*1e3:.3f} ms.')

    def take_snapshot(self, session, flags):
        data = {
            'surface':          self.surface,
            'mask':             self.mask,
            'max_components':   self.max_components,
        }
        from .. import CLIPPER_STATE_VERSION
        data['version']=CLIPPER_STATE_VERSION
        return data


    @classmethod
    def restore_snapshot(cls, session, data):
        surf = data['surface']
        if surf is None:
            return None     # Surface to mask is gone
        c = cls(surf, data['mask'], data['max_components'])
        surf.auto_remask_triangles = c
        c.set_surface_mask()
        return c


# class MultiZoneMask(State):
#     def __init__(self, surfaces, mask_volume, distance, max_components):
#         self.surfaces = surfaces
#         self.mask_volume = mask_volume
#         self.distance = distance
#         self.max_components = max_components
#         for s in surfaces:
#             self.set_surface_mask(s)
#
#     def __call__(self, surface=None):
#         if surface is None:
#             for s in self.surfaces:
#                 self.set_surface_mask(s)
#         else:
#             self.set_surface_mask(surface)
#
#     def set_surface_mask(self, surface):
#         v, t = surface.vertices, surface.triangles
#         if t is None:
#             return
#
#         indices = (self.mask_volume.interpolated_values(v)==1)
#         nv = len(v)
#         import numpy
#         mask = numpy.zeros((nv,), bool)
#         numpy.put(mask, indices, 1)
#         tmask = numpy.logical_and(mask[t[:,0]], mask[t[:,1]])
#         numpy.logical_and(tmask, mask[t[:,2]], tmask)
#         surface.triangle_mask = tmask
#
#         if self.max_components is not None:
#             from chimerax.surface import dust
#             dust.show_only_largest_blobs(surface, True, self.max_components)
#
#         surface.multi_remask_triangles = self
#
# class AtomZoneMask(MultiZoneMask):
#     def __init__(self, session, surfaces, step, distance, max_components,
#             atoms=None, coords=None, pad=0):
#         self.atoms = atoms
#         self.coords = coords
#         self.pad=pad
#         mv = VolumeMask(session, atoms.coords, step, distance, pad=pad)
#         super().__init__(surfaces, mv, distance, max_components)
#         self._mask_recalc_needed = False
#         self.reuse_existing_volume = True
#
#     def recalc_needed(self):
#         self._mask_recalc_needed = True
#
#     def set_surface_mask(self, surface):
#         if self._mask_recalc_needed:
#             self._recalculate_mask()
#             self._mask_recalc_needed=False
#         super().set_surface_mask(surface)
#
#     def _recalculate_mask(self):
#         if self.atoms is not None:
#             coords = self.atoms.coords
#         else:
#             coords = self.coords
#         self.mask_volume.generate_mask(coords, self.distance,
#             reuse_existing = False, pad=self.pad)
