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
from .mapset_base import MapSet_Base



class XmapSet_Box_Params:
    def __init__(self, origin_xyz=None, origin_grid=None, dim=None):
        self.origin_grid = origin_grid
        self.origin_xyz = origin_xyz
        self.dim = dim


    @property
    def origin_xyz(self):
        return self._origin_xyz

    @origin_xyz.setter
    def origin_xyz(self, origin):
        if origin is None:
            origin = [0,0,0]
        self._origin_xyz = numpy.array(origin, numpy.float32)
        # print("Set XYZ origin to {}".format(self._origin_xyz))

    @property
    def origin_grid(self):
        return self._origin_grid

    @origin_grid.setter
    def origin_grid(self, origin):
        if origin is None:
            origin = [0,0,0]
        self._origin_grid = numpy.array(origin, numpy.int)
        # print("Set grid origin to {}".format(self._origin_grid))

    @property
    def dim(self):
        return self._dim

    @dim.setter
    def dim(self, dim):
        if dim is None:
            dim = [1,1,1]
        self._dim = numpy.array(dim, numpy.int)
        # print("Set XYZ origin to {}".format(self._dim))


    def __iter__(self):
        return iter((self._origin_xyz, self._origin_grid, self._dim))

    def set(self, params):
        self.origin_xyz=params[0]
        self.origin_grid = params[1]
        self.dim = params[2]


_pad_base = numpy.array([-1,1], numpy.int)
class XmapSet(MapSet_Base):
    '''
    Handles the organisation, visualisation and re-calculation (where
    applicable) of crystallographic maps. Two general types of maps are
    supported: "static" maps generated from pre-calculated amplitudes and phases
    (e.g. as generated by your favourite refinement package); and "live" maps
    calculated on-the-fly directly from atomic coordinates and experimentally
    observed amplitude data. Each XmapSet handles the data from one crystal.
    Furthermore, live map recalculation requires a single Fobs/sigFobs dataset -
    if more than one is found, the user will be asked to choose.
    '''

    STANDARD_HIGH_CONTOUR=0.5
    STANDARD_LOW_CONTOUR = 0.65


    _default_live_xmap_params = {
        '2mFo-DFc': {'b_sharp': 0, 'is_difference_map': False, 'display': True},
        'mFo-DFc': {'b_sharp': 0, 'is_difference_map': True, 'display': True},
    }
    def __init__(self, manager, crystal_data,
        use_live_maps = True,
        use_static_maps = True,
        fsigf_name=None,
        bsharp_vals=None,
        exclude_free_reflections=False,
        fill_with_fcalc=False,
        exclude_missing_reflections=False,
        show_r_factors=True):
        '''
        Prepare the XmapSet and create all required maps.

        Args:
            * manager:
                - the master `Map_Mgr` object
            * crystal_data:
                - a ReflectionDataContainer object holding all the data for a
                  single crystal
            * use_live_maps:
                - if True, then the F/sigF data in crystal_data (if present)
                  will be used to generate maps that automatically update on
                  changes to the atomic model
            * use_static_maps:
                - if True, each F/phi dataset found in crystal_data will be
                  used to generate a static (unchanging) crystallographic map.
            * fsigf_name:
                - If the dataset has multiple F/sigF columns, setting this
                  value will explicitly choose the one to use. Otherwise, a
                  dialog will be launched allowing interactive selection.
            * bsharp_vals:
                - (only applicable when use_live_maps is True) a list of values
                  for sharpening/smoothing B-factors. For each value in the
                  list, a new map will be created from amplitudes scaled by the
                  B-factor value. Negative values lead to smoothing, positive
                  values to sharpening. If None, a B-factor will be chosen for
                  optimum viewing.
            * exclude_free_reflections
                - (only applicable when use_live_maps is True) if True, live
                  maps will use all measured reflections in their calculations.
                  NOTE: if the maps are to be used for fitting purposes, it is
                  strongly advised to set this argument to False.
            * fill_with_fcalc
                - (only applicable when use_live_maps is True) if True, missing
                  reflections will be approximated from the model for map
                  calculations. Note that this will lead to some level of model
                  bias (potentially serious of a large fraction of reflections
                  are missing).
            * exclude_missing_reflections
                - (only applicable when both use_live_maps,
                  exclude_free_reflections and fill_with_fcalc are true)
                  if True, only reflections that are excluded due to being part
                  of the free set will be replaced with Fcalc.
            * show_r_factors
                - (only applicable when use_live_maps is True) if True, R-work
                  and R-free will be reported to the status bar every time maps
                  are recalculated
        '''
        super().__init__(manager, 'Crystallographic maps ({})'.format(
            crystal_data.filename))

        session = self.session

        from .. import (
            HKL_data_F_sigF, HKL_data_F_sigF_ano,
            HKL_data_I_sigI, HKL_data_I_sigI_ano
        )
        self._known_experimental_datatypes = (
            HKL_data_F_sigF, HKL_data_F_sigF_ano,
            HKL_data_I_sigI, HKL_data_I_sigI_ano
        )

        trigger_names = (
            'maps recalculated',
        )
        for t in trigger_names:
            self.triggers.add_trigger(t)

        self._crystal_data = crystal_data
        self.add([crystal_data])


        # Since the Unit_Cell class does much of its work in grid coordinates,
        # each XmapSet object needs its own Unit_Cell instance (since there is
        # no guarantee that resolutions will be equal).
        from .. import Unit_Cell, atom_list_from_sel
        alist = atom_list_from_sel(self.structure.atoms)
        self._unit_cell = Unit_Cell(alist, self.cell, self.spacegroup,
            self.grid, 5)
        # Work out if we have experimental data, and launch the live data
        # manager if so.
        xm = self._live_xmap_mgr = None
        self._live_update = False
        self._recalc_needed = False
        self._model_changes_handler = None
        self._model_swap_handler = None
        self._delayed_recalc_handler = None
        self._show_r_factors = show_r_factors

        self._box_params = XmapSet_Box_Params()
        if self.spotlight_mode:
            self._box_changed_cb('init', (self.spotlight_center, self.display_radius))
        else:
            self.expand_to_cover_coords(self.master_map_mgr.last_covered_selection, 10)

        if use_live_maps:
            if fsigf_name is not None:
                fsigf = crystal_data.experimental_data.datasets.get(fsigf_name, None)
                if fsigf is None:
                    session.logger.info('WARNING: no experimental reflection data found with the '
                        'given label {}! Launching chooser dialog...'.format(fsigf_name))
                    fsigf_name = None

            if fsigf_name is None:
                choice = self._choose_reflection_data(crystal_data)
                if choice is not None:
                    fsigf_name, fsigf = choice
                else:
                    fsigf_name = fsigf = None

            if fsigf is not None:
                if fsigf.dtype in (HKL_data_I_sigI, HKL_data_I_sigI_ano):
                    session.logger.info('Reflection data provided as intensities. '
                        'Performing French & Wilson scaling to convert to amplitudes...')
                    from chimerax.clipper.reflection_tools import french_wilson_analytical
                    fsigf_data = french_wilson_analytical(fsigf.data)
                else:
                    fsigf_data = fsigf.data

                self._launch_live_xmap_mgr(crystal_data, fsigf_data)
                if bsharp_vals is None:
                    bsharp_vals = [viewing_recommended_bsharp(self.resolution)]
                if bsharp_vals:
                    if max(bsharp_vals) > 0:
                        map_style = 'mesh'
                        transparency = 0.0
                    else:
                        map_style='surface'
                        transparency=0.6
                map_params = self._default_live_xmap_params.copy()
                map_params['2mFo-DFc']['style'] = map_style
                map_params['2mFo-DFc']['transparency'] = transparency

                self._prepare_standard_live_maps(exclude_free_reflections,
                    fill_with_fcalc, exclude_missing_reflections,
                    map_params)
                from math import sqrt
                for b in bsharp_vals:
                    style = 'mesh'
                    transparency = 0.0
                    if b<0:
                        name_str = "2mFo-DFc_smooth_{:.0f}".format(-b)
                        contour = sqrt(self.resolution) / 3
                        #contour = self.STANDARD_LOW_CONTOUR
                    else:
                        name_str = "2mFo-DFc_sharp_{:.0f}".format(b)
                        if b == max(bsharp_vals):
                            style = 'surface'
                            transparency = 0.6
                            contour = sqrt(self.resolution) / 4
                            #contour = self.STANDARD_HIGH_CONTOUR
                    self.add_live_xmap(name_str,
                        b_sharp=b,
                        is_difference_map=False,
                        exclude_free_reflections=exclude_free_reflections,
                        fill_with_fcalc=fill_with_fcalc,
                        exclude_missing_reflections=exclude_missing_reflections,
                        style=style,
                        transparency=transparency,
                        contour=contour
                        )
        if self.live_xmap_mgr is not None:
            self.live_update=True
        if use_static_maps:
            cdata = crystal_data.calculated_data
            for dataset in cdata:
                xmap = self.add_static_xmap(dataset)
                if xmap.is_probably_fcalc_map():
                    xmap.display = False
        manager.add([self])


    @property
    def box_params(self):
        return self._box_params

    @property
    def live_xmap_mgr(self):
        return self._live_xmap_mgr

    @property
    def live_xmaps(self):
        return [m for m in self.child_models() if isinstance(m, XmapHandler_Live)]

    @property
    def static_xmaps(self):
        return [m for m in self.child_models() if isinstance(m, XmapHandler_Static)]

    @property
    def free_flags(self):
        return self._crystal_data.free_flags

    @property
    def show_r_factors(self):
        return self._show_r_factors

    @show_r_factors.setter
    def show_r_factors(self, flag):
        self._show_r_factors = flag

    @property
    def live_update(self):
        return self._live_update

    @live_update.setter
    def live_update(self, flag):
        if flag == self._live_update:
            return
        if flag:
            if self.live_xmap_mgr is None:
                self.session.logger.warning('Map set {} has no experimental '
                    'data! Live map recalculation is not possible.'.format(self.name))
            else:
                if self._model_changes_handler is None:
                    self._model_changes_handler = self.crystal_mgr.triggers.add_handler(
                        'atoms changed', self._model_changed_cb
                    )
                if self._model_swap_handler is None:
                    self._model_swap_handler = self.crystal_mgr.triggers.add_handler(
                        'model replaced', self._model_changed_cb
                    )
        else:
            mch = self._model_changes_handler
            msh = self._model_swap_handler
            for h in (mch, msh):
                if h is not None:
                    self.crystal_mgr.triggers.remove_handler(h)

            self._model_changes_handler = None
            self._model_swap_handler = None
        self._live_update = flag

    @property
    def rfree(self):
        if self.live_xmap_mgr is None:
            return None
        return self.live_xmap_mgr.rfree

    @property
    def rwork(self):
        if self.live_xmap_mgr is None:
            return None
        return self.live_xmap_mgr.rwork

    @property
    def experimental_data(self):
        return self._crystal_data.experimental_data.datasets

    @property
    def hklinfo(self):
        return self._crystal_data.hklinfo

    @property
    def cell(self):
        return self._crystal_data.cell

    @property
    def spacegroup(self):
        return self._crystal_data.spacegroup

    @property
    def resolution(self):
        return self._crystal_data.resolution.limit

    @property
    def grid(self):
        return self._crystal_data.grid_sampling

    @property
    def unit_cell(self):
        return self._unit_cell

    def _launch_live_xmap_mgr(self, crystal_data, f_sigf):
        from ..util import available_cores
        # The master C++ manager for handling all map operations
        from ..clipper_python.ext import Xtal_thread_mgr
        xm = self._live_xmap_mgr = Xtal_thread_mgr(self.hklinfo,
            crystal_data.free_flags.data, self.grid, f_sigf,
            num_threads=available_cores())
        # from ..util import atom_list_from_sel
        # ca = self._clipper_atoms = atom_list_from_sel(self.structure.atoms)
        atoms = self.structure.atoms
        xm.init(atoms.pointers)

    def _prepare_standard_live_maps(self, exclude_free_reflections,
        fill_with_fcalc, exclude_missing_reflections, map_params):
        xm = self.live_xmap_mgr
        if xm is None:
            return
        for name, params in map_params.items():
            self.add_live_xmap(name, **params,
                exclude_free_reflections=exclude_free_reflections,
                fill_with_fcalc=fill_with_fcalc,
                exclude_missing_reflections=exclude_missing_reflections
                )


    def set_xmap_display_style(self, xmap_handler, is_difference_map=False,
        color=None, style=None, transparency=0.0, contour=None):
        if style is None:
            style = 'mesh'
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
            if is_difference_map:
                contour = self.STANDARD_DIFFERENCE_MAP_CONTOURS
            else:
                from math import sqrt
                contour = sqrt(self.resolution) / 3
                # contour = self.STANDARD_LOW_CONTOUR


        if is_difference_map:
            contour = numpy.array(contour) * xmap_handler.sigma
        else:
            from ..util import guess_suitable_contour
            contour = [guess_suitable_contour(xmap_handler, self.structure, atom_radius_scale = contour)]

        xmap_handler.set_parameters(**{'cap_faces': False,
                                  'surface_levels': contour,
                                  'show_outline_box': False,
                                  'surface_colors': color,
                                  'square_mesh': False,
                                  'transparency_factor': transparency,
                                  })
        xmap_handler.set_display_style(style)

    def add_live_xmap(self, name,
        b_sharp=0,
        is_difference_map=False,
        exclude_missing_reflections=False,
        exclude_free_reflections=True,
        fill_with_fcalc=True,
        color=None,
        style=None,
        transparency=0.0,
        contour=None,
        display=True):
        '''
        Add a "live" crystallographic map, calculated from atomic positions and
        experimental x-ray reflection data. The map will be automatically
        recalculated in the background whenever changes are made to the atomic
        model.
        '''
        xm = self._live_xmap_mgr
        if xm is None:
            raise RuntimeError('This crystal dataset has no experimental amplitudes!')
        xm.add_xmap(name, b_sharp, is_difference_map=is_difference_map,
            exclude_missing_reflections=exclude_missing_reflections,
            exclude_free_reflections=exclude_free_reflections,
            fill_with_fcalc = fill_with_fcalc)
        new_handler = XmapHandler_Live(self, name,
            is_difference_map=is_difference_map)
        self.set_xmap_display_style(new_handler,
            is_difference_map=is_difference_map,
            color=color,
            style=style,
            transparency=transparency,
            contour=contour)
        self.add([new_handler])
        if display:
            new_handler.show()
            self.master_map_mgr.rezone_needed()
        else:
            new_handler.display=False
        return new_handler

    def add_static_xmap(self, dataset,
        is_difference_map=None,
        color=None,
        style=None,
        contour=None,
        display=True):
        '''
        Add a crystallographic map based on pre-calculated amplitudes and
        phases. This map will remain unchanged no matter what happens to the
        atomic model.
        '''
        data = dataset.data
        from .. import HKL_data_F_phi
        if not isinstance(data, HKL_data_F_phi):
            raise RuntimeError('Invalid data type: {}!'.format(type(data)))
        if is_difference_map is None:
            is_difference_map = dataset.is_difference_map
        new_handler = XmapHandler_Static(self, dataset.name, data,
            is_difference_map = is_difference_map)
        self.set_xmap_display_style(new_handler,
            is_difference_map=is_difference_map,
            color=color,
            style=style,
            contour=contour
            )
        self.add([new_handler])
        if display:
            new_handler.show()
            self.master_map_mgr.rezone_needed()
        else:
            new_handler.display=False
        return new_handler

    def save_mtz(self, filename, save_input_data=True, save_output_fobs=True,
            save_map_coeffs=False):
        from chimerax.clipper.io import mtz_write
        data = {}
        if save_output_fobs or save_map_coeffs:
            data['out'] = {}
        if save_input_data:
            data['in'] = {dname: d.data for dname, d in self.experimental_data.items()}
        if save_output_fobs:
            data['out']['[FOBS, SIGFOBS]'] = self.live_xmap_mgr.f_obs
        if save_map_coeffs:
            out = data['out']
            out['[FCALC, PHFCALC]'] = self.live_xmap_mgr.f_calc
            out['[2FOFC, PH2FOFC]'] = self.live_xmap_mgr.base_2fofc
            out['[FOFC, PHFOFC]'] = self.live_xmap_mgr.base_fofc
        mtz_write.write_mtz(filename, self.hklinfo, self.free_flags.data, data)

    def _choose_reflection_data(self, crystal_data):
        from .. import HKL_data_F_sigF
        exp_data = crystal_data.experimental_data
        reflections = {key: data for key, data in exp_data.datasets.items() if data.dtype in self._known_experimental_datatypes}
        # fsigfs = {key: data for key, data in exp_data.datasets.items() if data.dtype==HKL_data_F_sigF}
        if not len(reflections):
            return None
        elif len(reflections) == 1:
            ret = list(reflections.items())[0]
            return ret
        else:
            choice = self._reflection_data_chooser(list(reflections.keys()))
            if choice is None:
                return None
            return (choice, reflections[choice])

    def _reflection_data_chooser(self, possible_names,
        title = 'Choose the dataset to use for map calculations'):
        from PyQt5.QtWidgets import QInputDialog
        choice, ok_pressed = QInputDialog.getItem(self.session.ui.main_window,
            title, 'Label: ', possible_names, 0, False)
        if ok_pressed and choice:
            return choice
        return None

    def expand_to_cover_coords(self, coords, padding):
        cell = self.cell
        grid = self.grid
        pad = _calculate_grid_padding(padding, grid, cell)
        from ..clipper_util import get_minmax_grid
        box_bounds_grid = \
            get_minmax_grid(coords, cell, grid)
        box_bounds_grid[0] -= pad
        box_bounds_grid[1] += pad
        self.set_box_limits(box_bounds_grid, force_fill=True)

    def set_box_limits(self, minmax, force_fill = False):
        '''
        Set the map box to fill a volume encompassed by the provided minimum
        and maximum grid coordinates. Automatically turns off live scrolling.
        '''
        self.live_scrolling = False
        from .. import Coord_grid
        cmin = Coord_grid(minmax[0])
        cmin_xyz = cmin.coord_frac(self.grid).coord_orth(self.cell).xyz
        dim = (minmax[1]-minmax[0])
        self.box_params.set((cmin_xyz, minmax[0], dim))
        self.triggers.activate_trigger('map box changed', None)

    # Callbacks

    _map_impacting_changes = set ((
        "aniso_u changed",
        "bfactor changed",
        "coord changed",
        "coordset changed",
        "occupancy changed",
    ))
    def _model_changed_cb(self, trigger_name, changes):
        recalc_needed = False
        if trigger_name == 'model replaced':
            recalc_needed = True
        elif changes is not None:
            changes = changes[1]
            atom_changes = set(changes.atom_reasons()).intersection(self._map_impacting_changes)
            added = len(changes.created_atoms())
            deleted = changes.num_deleted_atoms()
            recalc_needed = (atom_changes or added or deleted)

        if recalc_needed:
            self._recalc_needed = True
            if self._delayed_recalc_handler is None:
                self._delayed_recalc_handler = self.session.triggers.add_handler(
                    'new frame', self._recalculate_maps_if_needed
                )

    def _recalculate_maps_if_needed(self, *_):
        xm = self.live_xmap_mgr
        if self._recalc_needed and not xm.thread_running:
            self.recalculate_all_maps(self.structure.atoms)
            if self._delayed_recalc_handler is not None:
                self.session.triggers.remove_handler(self._delayed_recalc_handler)
                self._delayed_recalc_handler = None

    def recalculate_all_maps(self, atoms):
        import ctypes
        # from .. import atom_list_from_sel
        from ..delayed_reaction import delayed_reaction
        xm = self.live_xmap_mgr
        # ca = self._clipper_atoms = atom_list_from_sel(atoms)
        delayed_reaction(self.session.triggers, 'new frame',
            xm.recalculate_all_maps, [atoms.pointers,],
            xm.ready,
            self._apply_new_maps, []
            )
        self._recalc_needed = False

    def _apply_new_maps(self):
        if self.deleted:
            return
        xm = self.live_xmap_mgr
        xm.apply_new_maps()
        if self.show_r_factors:
            rfactor_string = 'R-work: {:0.4f}  Rfree: {:0.4f}'.format(xm.rwork, xm.rfree)
            self.session.logger.status(rfactor_string, secondary=True)
        self.triggers.activate_trigger('maps recalculated', None)

    def delete(self):
        self.live_update = False
        self.stop_showing_r_factors()
        super().delete()

    # Callbacks

    def _cover_coords_cb(self, trigger_name, params):
        coords, padding = params
        self.expand_to_cover_coords(coords, padding)

    def _box_changed_cb(self, trigger_name, params):
        center, radius = params
        origin_grid, origin_xyz = _find_box_corner(center, self.cell, self.grid, radius=radius)
        dim = 2*_calculate_grid_padding(radius, self.grid, self.cell)
        self.box_params.set((origin_xyz, origin_grid, dim))
        self.triggers.activate_trigger('map box changed', None)

    def _box_moved_cb(self, trigger_name, center):
        bp = self.box_params
        bp.origin_grid, bp.origin_xyz = _find_box_corner(center, self.cell,
            self.grid, radius=self.display_radius)
        self.triggers.activate_trigger('map box moved', None)

    def start_showing_r_factors(self):
        if not self.live_xmap_mgr:
            self.session.logger.warning('Attempted to start live R-factor display, '
                'but no live crystallographic maps are loaded. Command ignored.')
            return
        if not hasattr(self, '_r_factor_label') or self._r_factor_label is None:
            from chimerax.label import Label
            self._r_factor_label = Label(self.session, self.name+'_r_factors',
                text='', xpos=0.025, ypos=0.025, bold=True, size=30)
            self._r_factor_label_handler = self.triggers.add_handler(
                'maps recalculated', self._update_r_factor_text
            )
        self._update_r_factor_text()

    def _update_r_factor_text(self, *_):
        lb = self._r_factor_label
        lxm = self.live_xmap_mgr
        if not lxm:
            self._r_factor_label_handler = None
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER
        lb.text = 'Rwork: {:0.3f} Rfree: {:0.3f}  '.format(self.rwork, self.rfree)
        lb.update_drawing()

    def stop_showing_r_factors(self, *_):
        if hasattr(self, '_r_factor_label_handler') and self._r_factor_label_handler is not None:
            self.triggers.remove_handler(self._r_factor_label_handler)
            self._r_factor_label_handler = None
        if hasattr(self, '_r_factor_label') and self._r_factor_label is not None:
            self._r_factor_label.delete()
            self._r_factor_label = None



from .map_handler_base import XmapHandler_Base

class XmapHandler_Static(XmapHandler_Base):
    '''
    An XmapHandler_Static is in effect a resizable window into a periodic
    crystallographic map. The actual map data (a clipper Xmap object) is
    calculated via fast Fourier transform from a pre-calculated set of
    amplitudes and phases (e.g. as provided by a refinement package), and filled
    into the XmapHandler_Static.data array as needed. Mothods are provided for
    tracking and filling a box around the centre of rotation, and static display
    of a given region.
    '''
    def __init__(self, mapset, name, f_phi_data,
        is_difference_map=False):
        '''
        Args:
            mapset:
                The XmapSet object this belongs to
            name:
                A descriptive name for this map
            f_phi_data:
                The Clipper HKL_Data_F_phi object the map is to be calculated
                from.
            is_difference_map:
                Is this a difference map?
        '''
        name = '(STATIC) '+name
        self._mapset = mapset
        from .. import Xmap

        xmap = self._xmap = Xmap(self.spacegroup, self.cell, self.grid)
        xmap.fft_from(f_phi_data)
        super().__init__(mapset, name, is_difference_map=is_difference_map)


    @property
    def xmap(self):
        return self._xmap

    @property
    def stats(self):
        if not hasattr(self, '_stats') or self._stats is None:
            from .. import Map_stats
            all = self._all_stats = Map_stats(self._xmap)
            self._stats = (all.mean, all.std_dev, all.std_dev)
        return self._stats

    def is_probably_fcalc_map(self):
        '''
        Attempt to guess if a map loaded from a file is Fcalc (that is, simply
        calculated from the atomic coordinates and properties). All we can really go
        on here is some simple heuristics.
        '''
        name = self.name.upper()
        if name.startswith('FC'):
            return True
        if 'CALC' in name:
            return True
        return False



class XmapHandler_Live(XmapHandler_Base):
    '''
    An XmapHandler_Live is in effect a resizable window into a periodic
    crystallographic map. The actual map data (a clipper Xmap object) is
    recalculated as needed from the combination of atomic coordinates and
    measured diffraction data, and filled into the XmapHandler_Live.data array
    as needed. Mothods are provided for live recalculation, tracking and filling
    a box around the centre of rotation, and static display of a given region.
    '''
    def __init__(self, mapset, name,
        is_difference_map=False):
        '''
        Args:
            mapset:
                The XmapSet object this belongs to
            name:
                A descriptive name for this map
            is_difference_map:
                Is this a difference map?
        '''
        self._map_name = name
        name = '(LIVE) '+name
        super().__init__(mapset, name, is_difference_map=is_difference_map)
        self._mgr_handlers.append(
            (mapset,
            mapset.triggers.add_handler('maps recalculated', self._map_recalc_cb
            ))
        )

    @property
    def xmap_mgr(self):
        return self.mapset.live_xmap_mgr

    @property
    def xmap(self):
        return self.xmap_mgr.get_xmap_ref(self._map_name)

    @property
    def stats(self):
        all = self.xmap_mgr.get_map_stats(self._map_name)
        return (all.mean, all.std_dev, all.std_dev)

    def _map_recalc_cb(self, name, *_):
        if self.deleted:
            return
        for s in self.surfaces:
            s._use_thread=True
        self._fill_volume_data(self._data_fill_target, self.box_params.origin_grid)
        self.data.values_changed()

    def delete(self):
        self.xmap_mgr.delete_xmap(self._map_name)
        super().delete()

def map_potential_recommended_bsharp(resolution):
    '''
    Return a recommended sharpening/smoothing B-factor for a given map to
    optimise its use as an MDFF potential. For now this is a simple linear
    function of resolution, passing through zero at 2.5 Angstroms (smoothing
    below, sharpening above). Improved methods will be developed over time.
    '''
    # smooth by 30 A**2 at 1.5A res; sharpen by 30A**2 at 3.5A res; 0 at 2.5A res
    return viewing_recommended_bsharp(resolution, crossover=2.5, low_res_base=35,
        high_res_base=20)

def viewing_recommended_bsharp(resolution, crossover=2.5, low_res_base=45,
        high_res_base=25):
    '''
    For viewing purposes it is also often useful to have a smoothed or
    sharpened visualisation of your map, but the optimal degree of sharpening
    for visualisation is not necessarily the same as that for MDFF.
    '''
    # smooth by 50 A**2 at 1.5A res; sharpen by 50 A**2 at 3.5A res, 0 at 2.5A res
    from math import sqrt
    if resolution < crossover:
        return -high_res_base*sqrt(crossover-resolution) + 2.5
    else:
        return low_res_base*sqrt(resolution-crossover) - 2.5



def _calculate_grid_padding(radius, grid, cell):
    '''
    Calculate the number of grid steps needed on each crystallographic axis
    in order to capture at least radius angstroms in x, y and z.
    '''
    corner_mask = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,0,0],[1,1,1]])
    corners = (corner_mask * radius).astype(numpy.float32)
    grid_upper = numpy.zeros([8,3], numpy.int)
    grid_lower = numpy.zeros([8,3], numpy.int)
    from .. import Coord_orth
    for i, c in enumerate(corners):
        co = Coord_orth(c)
        cm = co.coord_frac(cell).coord_map(grid).uvw
        grid_upper[i,:] = numpy.ceil(cm).astype(int)
        grid_lower[i,:] = numpy.floor(cm).astype(int)
    return grid_upper.max(axis=0) - grid_lower.min(axis=0)


def _find_box_corner(center, cell, grid, radius = 20):
    '''
    Find the bottom corner (i.e. the origin) of a rhombohedral box
    big enough to hold a sphere of the desired radius.
    '''
    from .. import Coord_frac, Coord_orth, Coord_grid
    radii_frac = Coord_frac(radius/cell.dim)
    center_frac = Coord_orth(center).coord_frac(cell)
    bottom_corner_grid = center_frac.coord_grid(grid) \
                - Coord_grid(_calculate_grid_padding(radius, grid, cell))
    bottom_corner_orth = bottom_corner_grid.coord_frac(grid).coord_orth(cell)
    return bottom_corner_grid.uvw, bottom_corner_orth.xyz

def _get_bounding_box(coords, grid, cell, padding=2):
    '''
    Find the minimum and maximum grid coordinates of a box which will
    encompass the given (x,y,z) coordinates plus padding (in Angstroms).
    '''
    from .. import Util, Coord_grid
    grid_pad = _calculate_grid_padding(padding, grid, cell)
    box_bounds_grid = Util.get_minmax_grid(coords, cell, grid)\
                        + numpy.array((-grid_pad, grid_pad))
    box_origin_grid = box_bounds_grid[0]
    box_origin_xyz = Coord_grid(box_origin_grid).coord_frac(grid).coord_orth(cell)
    dim = box_bounds_grid[1] - box_bounds_grid[0]
    return [box_origin_grid, box_origin_xyz, dim]
