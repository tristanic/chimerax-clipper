# @Author: Tristan Croll <tic20>
# @Date:   22-May-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 22-May-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll

from chimerax.map import Volume


class VolumeMask(Volume):
    '''
    Manages a binary mask defining displayed regions for other surfaces.
    '''



    def __init__(self, session, coords, step, radius):
        import numpy

        darray = self._generate_data_array(coords, step, radius)
        super().__init__(session, darray)
        self.name = "mask"
        self.generate_mask(coords, radius, clear=False)

    def _generate_data_array(self, coords, step, radius):
        import numpy
        if hasattr(step, '__len__'):
            step = min(step)

        step = numpy.ones(3) * max((step, radius/3))
        cmin = coords.min(axis=0)
        cmax = coords.max(axis=0)
        dim = numpy.ceil(((cmax+radius)-(cmin-radius))/step).astype(numpy.int)
        data = self._data_fill_target = numpy.zeros(dim, numpy.uint8)
        from chimerax.map.data import ArrayGridData
        darray = ArrayGridData(data.transpose(), origin=cmin-radius, step=step)
        return darray

    def generate_mask(self, coords, radius, clear=True):
        from chimerax.clipper import _map_mask
        origin, step = self.data_origin_and_step()
        _map_mask.generate_mask(self._data_fill_target, self.data.origin, step,
            self._data_fill_target.shape, self.data.ijk_to_xyz_transform.matrix,
            self.data.xyz_to_ijk_transform.matrix, coords, len(coords), radius)
        self.data.values_changed()

    def clear(self):
        self.data.array[:] = 0
        self.data.values_changed()
