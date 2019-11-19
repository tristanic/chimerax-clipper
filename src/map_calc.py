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


def generate_fcalc(container, atoms, fsigfdata, target=None):
    '''
    Generate a crystallographic map from the given atoms and observed
    structure factors.

    Args:
        * atoms:
            - A :class:`chimerax.Atoms` object containing all atoms to include
              in the calculation.
        * fsigfdata:
            - A :class:`ReflectionDataExp` object holding the observed
              structure factors
        * target:
            - A :class:`ReflectionDataCalc` object, or None. If not None, any
              data in the existing
    '''
    session = fsigfdata.session
    fsigf = fsigfdata.data
    hkls = container.hklinfo
    from . import atom_list_from_sel
    clipper_atoms = atom_list_from_sel(atoms)

    if target is None:
        from .clipper_mtz import ReflectionDataCalc
        from .clipper_python.data64 import HKL_data_F_phi_double
        target = ReflectionDataCalc('Fcalc', session,
            HKL_data_F_phi_double(hkls), is_difference_map=False)
    fcalc = target.data
    from .clipper_python import SFcalc_obs_bulk_double
    SFcalc_obs_bulk_double(fcalc, fsigf, clipper_atoms)
    container.calculated_data.add([target])

def generate_map_coeffs(container, fsigf, fcalc, free_r_flags):
    from .clipper_python.data64 import HKL_data_F_phi_double, HKL_data_Phi_fom_double
    session = container.session
    hkls = container.hklinfo
    best_coeffs = HKL_data_F_phi_double(hkls)
    diff_coeffs = HKL_data_F_phi_double(hkls)
    phiw = HKL_data_Phi_fom_double(hkls)

    from .clipper_python import SFweight_spline_double
    SFweight_spline_double(best_coeffs, diff_coeffs, phiw, fsigf.data, fcalc.data, free_r_flags.data)

    ret = []
    from .clipper_mtz import ReflectionDataCalc
    ret.append(ReflectionDataCalc('2FOFCWT, PH2FOFCWT', session, best_coeffs, is_difference_map=False))
    ret.append(ReflectionDataCalc('FOFCWT, PHFOFCWT', session, diff_coeffs, is_difference_map=True))
    return ret
    # for b in sharpening_factors:
    #     if b == 0:
    #         ret.append(best_coeffs)
    #     else:












    # if sharpening is None:
    #     if target is not None and hasattr(target, 'sharpening'):
    #         sharpening = target.sharpening
    #     else:
    #         sharpening = 0
    #
    # if sharpening < 0:
    #     name_ext = " sharp {:0.1f}".format(sharpening)
    # elif sharpening > 0:
    #     name_ext = " smooth {:0.1f}".format(sharpening)
    # else:
    #     name_ext = ""
    # name_string = type + name_ext
