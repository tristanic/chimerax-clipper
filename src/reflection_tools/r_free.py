# @Author: Tristan Croll <tic20>
# @Date:   06-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 06-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll

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
    print('N_obs: {} N_tot: {}, flags shape: {}, flag_array shape: {}'.format(
        nobs, ntot, flags.shape, flag_array.data[1].shape
    ))
    flag_array.set_data(flag_array.data[0],flags)
    free_obs = 0
    ih = data_array.first_data
    while not ih.last():
        if flag_array[ih].data[0] == 0:
            free_obs += 1
        data_array.next_data(ih)
    return free_obs
