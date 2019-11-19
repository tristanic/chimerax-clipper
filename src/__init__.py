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

from .main import *
# General objects
from .clipper_python import (
    Cell,
    Cell_descr,
    CCP4MTZfile,
    CIFfile,
    Coord_frac,
    Coord_grid,
    Coord_orth,
    Grid,
    Grid_range,
    Grid_sampling,
    HKL_info,
    Map_stats,
    Resolution,
    RTop_frac,
    RTop_orth,
    Spacegroup,
    Spgr_descr,
    Symop,
    # Symops,
    Unit_Cell,
    Util,
    Xmap_float as Xmap,
    )


# Singular forms of HKL data
from .clipper_python import Flag, Flag_bool
from .clipper_python.data32 import (
    ABCD_float as ABCD,
    D_sigD_float as D_sigD,
    E_sigE_float as E_sigE,
    F_phi_float as F_phi,
    F_sigF_float as F_sigF,
    F_sigF_ano_float as F_sigF_ano,
    I_sigI_float as I_sigI,
    I_sigI_ano_float as I_sigI_ano,
    Phi_fom_float as Phi_fom,
    )

# Array forms of HKL data
from .clipper_python import HKL_data_Flag, HKL_data_Flag_bool
from .clipper_python.data32 import (
    HKL_data_ABCD_float as HKL_data_ABCD,
    HKL_data_D_sigD_float as HKL_data_D_sigD,
    HKL_data_E_sigE_float as HKL_data_E_sigE,
    HKL_data_F_phi_float as HKL_data_F_phi,
    HKL_data_F_sigF_float as HKL_data_F_sigF,
    HKL_data_F_sigF_ano_float as HKL_data_F_sigF_ano,
    HKL_data_I_sigI_float as HKL_data_I_sigI,
    HKL_data_I_sigI_ano_float as HKL_data_I_sigI_ano,
    HKL_data_Phi_fom_float as HKL_data_Phi_fom,
    )

from .clipper_mtz import ReflectionDataContainer

from .symmetry import get_symmetry_handler, get_map_mgr


from chimerax.core.toolshed import BundleAPI
class _ClipperBundle(BundleAPI):
    from chimerax.core.commands import FloatArg
    from chimerax.atomic import StructureArg

    @staticmethod
    def initialize(session, bundle_info):
        from chimerax.clipper import cmd
        # cmd.register_mtz_file_format(session)

    # @staticmethod
    # def fetch_from_database(session, identifier, ignore_cache=False,
    #         database_name=None, format_name=None, **kw):
    #     from .io import fetch_cif
    #     fetchers = {
    #         'pdb':  fetch_cif.fetch_structure_factors,
    #         'pdbe': fetch_cif.fetch_structure_factors_pdbe,
    #         'pdbe_updated': fetch_cif.fetch_structure_factors_pdbe,
    #         'pdbj': fetch_cif.fetch_structure_factors_pdbj,
    #     }
    #     try:
    #         fetcher = fetchers[database_name]
    #     except KeyError:
    #         from chimerax.core.errors import UserError
    #         raise UserError("Unknown database for fetching structure factors: {} Known databases are: {}".format(
    #             database_name, ', '.join(fetchers.keys())
    #         ))
    #     return fetcher(session, identifier, ignore_cache=ignore_cache, **kw)

    @staticmethod
    def register_command(command_name, logger):
        # 'register_command' is lazily called when the command is referenced
        from chimerax.clipper import cmd
        if command_name == 'clipper':
            cmd.register_clipper_cmd(logger)

    @staticmethod
    def open_file(session, path, format_name, structure_model=None,
            over_sampling=2.0):
        if format_name in ('mtz', 'sfcif'):
            if structure_model is None:
                from chimerax.core.errors import UserError
                raise UserError('Must specify a structure model to associate with crystallographic data')
            from .cmd import open_structure_factors
            return open_structure_factors(session, path, structure_model=structure_model,
                over_sampling=over_sampling)

    @staticmethod
    def save_file(session, path, *, models=None, preserve_input=False,
            save_map_coeffs=False):
        from .cmd import save_structure_factors
        return save_structure_factors(session, path, models=models,
            preserve_input=preserve_input, save_map_coeffs=save_map_coeffs)

    # get_class() is used for saving/restoring tools, but Clipper doesn't define
    # any.
    @staticmethod
    def get_class(class_name):
        from .symmetry import SymmetryManager, AtomicSymmetryModel
        from .graphics.hkl_plot import HKLPlot3D
        from .maps import (
            XmapSet, MapMgr, XmapHandler_Static,
            XmapHandler_Live, NXmapSet, NXmapHandler
            )
        from .maps.map_handler_base import FastVolumeSurface
        from .maps.mask_handler import VolumeMask, ZoneMask
        from .clipper_mtz import (
            ReflectionDataContainer, ReflectionDataNode, ReflectionData,
            ReflectionDataFreeFlags, ReflectionDataExp, ReflectionDataCalc,
        )
        ct = {
            'SymmetryManager':     SymmetryManager,
            'AtomicSymmetryModel':  AtomicSymmetryModel,
            'HKLPlot3D':          HKLPlot3D,
            'XmapSet':              XmapSet,
            'MapMgr':              MapMgr,
            'XmapHandler_Static':   XmapHandler_Static,
            'XmapHandler_Live':     XmapHandler_Live,
            'NXmapSet':             NXmapSet,
            'NXmapHandler':         NXmapHandler,
            'FastVolumeSurface':    FastVolumeSurface,
            'VolumeMask':           VolumeMask,
            'ZoneMask':             ZoneMask,
            'ReflectionDataContainer': ReflectionDataContainer,
            'ReflectionDataNode':  ReflectionDataNode,
            'ReflectionData':       ReflectionData,
            'ReflectionDataFreeFlags': ReflectionDataFreeFlags,
            'ReflectionDataExp':   ReflectionDataExp,
            'ReflectionDataCalc':  ReflectionDataCalc
        }
        return ct.get(class_name)

bundle_api = _ClipperBundle()
