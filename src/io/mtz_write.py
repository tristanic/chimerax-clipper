# @Author: Tristan Croll <tic20>
# @Date:   04-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 04-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll

def write_mtz(filename, hkl_info, free_flags, datasets):
    '''
    Write a set of experimental and/or calculated structure factors to a MTZ
    file.

    Args:
        filename: a string
        hkl_info: a Clipper HKL_info object
        datasets: a {dataset_name: {data_name: data}} dict of dicts where
            dataset_name is a unique name, data_name is either a single unique
            name or a string in the format "[col1, col2, ...]" providing
            explicit names for each column in data, and data is a
            Clipper HKL_data array type.
    '''
    from chimerax.clipper.clipper_python import CCP4MTZfile
    if not filename.lower().endswith('.mtz'):
        filename += '.mtz'
    outfile = CCP4MTZfile()
    outfile.open_write(filename)
    outfile.export_hkl_info(hkl_info)
    outfile.export_hkl_data(free_flags, '/*/*/[Free_Flags]')
    seen_types = []
    for dataset_name, dataset in datasets.items():
        # cpath = '/{}/{}/'.format(crystal_name, dataset_name) + '[{}]'
        cpath = '/*/*/{}'
        for data_name, data in dataset.items():
            if data_name[0]=='[' and data_name[-1]==']':
                if len(data_name.split(',')) != data.data_size():
                    _bail_out(outfile)
                    from chimerax.core.errors import UserError
                    err_string = ('Number of column labels in "{}" does not '
                    'match number of columns in data (should be {}). Writing '
                    'has been cancelled.')
                    raise UserError(err_string.format(data_name, data.data_size()))
                dpath = cpath.format(data_name)
            else:
                dtype = type(data)
                if dtype in seen_types:
                    _bail_out(outfile)
                    from chimerax.core.errors import UserError
                    raise UserError('Multiple datasets of the same type were '
                        'found without explicit column labels. Writing has '
                        'been cancelled. To write an MTZ file, either remove '
                        'the extra dataset(s) or provide explicit column names '
                        'in the format "[col1, col2, ...]".')
                dpath = cpath.format(dataset_name)
            print("Writing data to: '{}'".format(dpath))
            outfile.export_hkl_data(data, dpath)
    outfile.close_write()

def _bail_out(outfile):
    outfile.close_write()
    import os
    os.remove(outfile)
