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

def generate_free_set(flag_array, data_array, free_frac=0.05, max_free=2000):
    '''
    Reassign the values in flag_array to generate a new random set of test
    reflections, following the CCP4 convention (that is, free reflections have
    a flag value of 0). Assigns a value to every element of the array,
    regardless of whether the HKL is observed in the dataset. Returns the number
    of free reflections assigned.

    Args:
        flag_array: a Clipper HKL_data_Flag object
        data_array: a Clipper HKL_data_{F_sigF, F_sigF_ano, I_sigI or I_sigI_ano} array
        free_frac: the target fraction of reflections to be marked free.
        max_free: the maximum number of free reflections (in large datasets we
            don't want to throw away *too* much data)
    '''
    nobs = data_array.num_obs
    ntot = flag_array.hkl_info.num_reflections
    free_frac = min(free_frac, max_free/nobs)
    random_range = int(1/free_frac)
    import numpy
    test_vals = numpy.random.randint(0, high=random_range, size=ntot)
    flags = test_vals.astype(numpy.double).reshape([ntot, 1])
    # print('N_obs: {} N_tot: {}, flags shape: {}, flag_array shape: {}'.format(
    #     nobs, ntot, flags.shape, flag_array.data[1].shape
    # ))
    flag_array.set_data(flag_array.data[0],flags)
    free_obs = 0
    ih = data_array.first_data
    while not ih.last():
        if flag_array[ih].data[0] == 0:
            free_obs += 1
        data_array.next_data(ih)
    return free_obs
