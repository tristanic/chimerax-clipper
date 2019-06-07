# @Author: Tristan Croll <tic20>
# @Date:   05-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 07-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll

def _first_column_types():
    from chimerax.clipper import (
        HKL_data_Flag,
        HKL_data_D_sigD,
        HKL_data_F_sigF,
        HKL_data_F_phi,
        HKL_data_F_sigF_ano,
        HKL_data_I_sigI,
        HKL_data_I_sigI_ano,
    )
    first_column_type = {
        'D':    HKL_data_D_sigD,
        'F':    {'Q': HKL_data_F_sigF, 'P': HKL_data_F_phi},
        'G':    HKL_data_F_sigF_ano,
        'J':    HKL_data_I_sigI,
        'K':    HKL_data_I_sigI_ano,
        'I':    HKL_data_Flag
    }

    expected_following_columns = {
        'D':    ['Q'],
        # 'F':    ['Q'],
        'G':    ['L','G','L'],
        'J':    ['Q'],
        'K':    ['M','K','M'],
        'I':    [],
    }



    return (first_column_type, expected_following_columns)


def _defaultdict_as_dict(d):
    '''
    Convert a recursively-defined defaultdict back to a normal dict for easier
    downstream handling.
    '''
    for k, v in d.items():
        if isinstance(v, dict):
            d[k] = _defaultdict_as_dict(v)
        return dict(d)

def _parse_column_path(cpath):
    path, dtype = cpath.split(' ')
    crystal, dataset, name = path.split('/')[1:]
    return (crystal, dataset, name, dtype)

def load_mtz_data(session, filename):
    from chimerax.clipper import CCP4MTZfile, HKL_info

    mtzin = CCP4MTZfile()
    hklinfo = HKL_info()
    mtzin.open_read(filename)
    mtzin.import_hkl_info(hklinfo)

    first_column_type, expected_following_columns = _first_column_types()

    # Get all column names and types
    column_paths = [str(p) for p in mtzin.column_paths]

    from collections import defaultdict
    def nested_default_dict():
        return defaultdict(lambda: nested_default_dict())
    crystals = defaultdict(nested_default_dict)

    i = 0
    while i < len(column_paths):
        current_path = column_paths[i]
        crystal, dataset, name, dtype = _parse_column_path(current_path)
        if dtype == 'H':
            # HKLs, already captured in hklinfo
            i += 1
            continue
        array_type = first_column_type.get(dtype, None)
        if array_type is None:
            session.logger.info('Discarding unrecognised/unsupported data array {}'.format(current_path))
            i += 1
            continue

        if dtype == 'F':
            ncryst, ndat, nname, ndtype = _parse_column_path(column_paths[i+1])
            real_array_type = array_type.get(ndtype, None)
            if real_array_type is None:
                warn_str = ('WARNING: found column "{}" suggesting one of the '
                'data types ({}), but the following column type "{}" does not '
                'match any of the expected types ({}). Importing of this '
                'column has been skipped.')
                warn_str = warn_str.format(current_path,
                    ', '.join([atype.__name__ for atype in array_type.values()]),
                    ndtype, ', '.join(array_type.keys()))
                session.logger.warning(warn_str)
                i += 1
                continue
            next_expected = [ndtype]
            array_type = real_array_type
        # print('Path {} looks like {}. Trying...'.format(current_path, array_type.__name__))
        else:
            next_expected = expected_following_columns[dtype]

        import_names = [name]
        next_found_names = []
        next_found_dtypes = []
        for edtype in next_expected:
            i += 1
            npath = column_paths[i]
            ncryst, ndat, nname, ndtype = _parse_column_path(npath)
            if ncryst != crystal or ndat != dataset:
                warn_str = ('WARNING: found column "{}" suggesting data of type '
                '{}, but reached the end of dataset {} before finding the '
                'expected matching columns. Make sure your MTZ file is laid out '
                'in a logical order (e.g. F,sigF; F+,sigF+,F-,sigF-)')
                warn_str.format(current_path, array_type.__name__,
                    '/'.join(crystal, dataset))
                session.logger.warn(warn_str)
                break
            next_found_names.append(nname)
            next_found_dtypes.append(ndtype)
            if edtype != ndtype:
                warn_str = ('WARNING: found column "{}" suggesting data of type '
                '{}, but the data in the following column(s) is not as '
                'expected. For this data type the following columns should have '
                'types {}; found {} ({}). This column has been ignored.')
                warn_str.format(current_path, array_type.__name__,
                    ','.join(next_expected), ','.join(next_found_dtypes),
                    ','.join(next_found_names))
                next_found_names = next_found_names[:-1]
                next_found_dtypes = next_found_dtypes[:-1]
                break
            # print('Found {} dtype {}'.format(nname, ndtype))
        if len(next_found_names) == len(next_expected):
            import_path = '/'.join([crystal, dataset])
            import_names = [name] + next_found_names
            dset_name = ', '.join(import_names)
            import_path += '/'+'[{}]'.format(dset_name)
            this_data = array_type(hklinfo)
            mtzin.import_hkl_data(this_data, import_path)
            crystals[crystal][dataset][dset_name] = this_data
        i += 1

    return hklinfo, _defaultdict_as_dict(crystals)
