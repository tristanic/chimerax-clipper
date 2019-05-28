# @Author: Tristan Croll <tic20>
# @Date:   22-May-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 28-May-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll

from chimerax.map import Volume

class Zone_Mgr:
    def __init__(self, session, grid_step, radius,
            atoms=None, coords=None, pad=None, interpolation_threshold=0.75):
        if pad is None:
            pad = radius
        self.session = session
        self._step = grid_step
        self._radius = radius
        self._atoms = atoms
        self._coords = coords
        self._pad = pad
        self.threshold=interpolation_threshold
        self._mask = None
        self._update_needed = False
        self._resize_box = True

    @property
    def coords(self):
        if self._atoms is not None:
            return self.atoms.coords
        return self._coords

    @coords.setter
    def coords(self, coords):
        if self.coords is not None:
            prev_coord_len = len(self.coords)
        else:
            prev_coord_len = None
        self._atoms = None
        self._coords = coords
        if len(coords) == 1 and prev_coord_len==1:
            self.update_needed(resize_box=False)
        else:
            self.update_needed(resize_box=True)

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        self._coords = None
        self._atoms = atoms
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
        update_origin = (len(self.coords) == 1)
        self.mask.generate_mask(self.coords, self.radius,
            reuse_existing=not self._resize_box, update_origin=update_origin,
            step=self.grid_step, pad=self.pad)
        self._update_needed = False
        self._resize_box = False

    def get_vertex_mask(self, vertices):
        if self._update_needed:
            self._update_mask()
        return (self.mask.interpolated_values(vertices) >= self.threshold)


class VolumeMask(Volume):
    '''
    Manages a binary mask defining displayed regions for other surfaces.
    '''

    def __init__(self, session, coords, step, radius, pad=0):
        import numpy

        darray = self._generate_data_array(coords, step, radius, pad=pad)
        super().__init__(session, darray)
        self.name = "mask"
        self.generate_mask(coords, radius, clear=False, pad=pad)

    def _generate_data_array(self, coords, step, radius, pad=0):
        import numpy
        if hasattr(step, '__len__'):
            step = min(step)

        step = max((step, radius/4))
        step = min(step, 3)
        step = numpy.ones(3) * step
        cmin = coords.min(axis=0)
        cmax = coords.max(axis=0)
        mincoor = cmin-radius-pad
        maxcoor = cmax+radius+pad
        dim = numpy.ceil((maxcoor-mincoor)/step).astype(numpy.int)
        data = numpy.zeros(dim, numpy.uint8)
        self._data_fill_target = data
        from chimerax.map.data import ArrayGridData
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
    def __init__(self, surface, mgr, max_components):
        self.surface = surface
        self.mgr = mgr
        self.max_components = max_components
        self.set_surface_mask()

    def __call__(self):
        self.set_surface_mask()

    def set_surface_mask(self):
        surface = self.surface
        v, t = surface.vertices, surface.triangles
        if t is None:
            return

        import numpy
        indices = numpy.where(self.mgr.get_vertex_mask(v))[0]
        nv = len(v)
        mask = numpy.zeros((nv,), numpy.bool)
        numpy.put(mask, indices, 1)
        tmask = numpy.logical_and(mask[t[:,0]], mask[t[:,1]])
        numpy.logical_and(tmask, mask[t[:,2]], tmask)
        surface.triangle_mask = tmask

        if self.max_components is not None:
            from chimerax.surface import dust
            dust.show_only_largest_blobs(surface, True, self.max_components)

        surface.auto_remask_triangles = self

    def take_snapshot(self, session, flags):
        data = {
            'surface':          self.surface,
            'mask':             self.mask,
            'max_components':   self.max_components,
            'version':          1,
        }

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
#         mask = numpy.zeros((nv,), numpy.bool)
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
