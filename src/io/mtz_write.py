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
