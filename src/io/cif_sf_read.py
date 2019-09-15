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

_cif_sf_table_names = (
    'audit',                        # General information about file
    'cell',                         # Crystallographic cell information
    'diffrn',                       # General diffraction information
    'diffrn_radiation_wavelength',  # X-ray wavelength
    'entry',                        # PDB ID
    'exptl_crystal',                # Unique crystal identifier for multi-dataset cifs
    'reflns_scale',                 #
    'symmetry',                     # Space group identification
    'reflns',                       # General information about reflections
    'refln',                        # Actual reflection data
    )

_metadata_tables = ('audit', 'diffrn', 'diffrn_radiation_wavelength', 'entry',
    'exptl_crystal', 'reflns')

_required_columns = {
    'cell':       ('length_a', 'length_b', 'length_c',
                    'angle_alpha', 'angle_beta', 'angle_gamma'),
    'refln':      ('index_h', 'index_k', 'index_l')
}

def _space_group_identifiers():
    from chimerax.clipper.clipper_python import Spgr_descr
    return (
        ('Int_Tables_number', Spgr_descr.TYPE.Number),
        ('space_group_name_Hall', Spgr_descr.TYPE.Hall),
        ('space_group_name_H-M', Spgr_descr.TYPE.HM),
    )

def _expand_hm_symbol(symbol_hm):
    s = symbol_hm.split(' ')
    if len(s) != 2:
        return symbol_hm
    lattice, sym = s
    if sym in ('2', '21', 'm', 'c', '2/m', '21/m', '2/c', '21/c'):
        return ' '.join([lattice, '1', sym, '1'])
    return symbol_hm

_ccp4_symop_code_to_hm_symbol = {
    '1003': 'P 1 2 1',
    '1004': 'P 1 1 21',
    '1005': 'B 1 1 2',
    '1006': 'P 1 1 m',
    '1007': 'P 1 1 b',
    '1008': 'B 1 1 m',
    '1009': 'B 1 1 b',
    '1010': 'P 1 1 2/m',
    '1011': 'P 1 1 21/m',
    '1012': 'B 1 1 2/m',
    '1013': 'P 1 1 2/b',
    '1014': 'P 1 1 21/b',
    '1015': 'B 1 1 2/b',
    '1017': 'P 21 2 2',
    '1059': 'P m m n',
    '1146': 'R 3',
    '1148': 'R -3',
    '1155': 'R 3 2',
    '1160': 'R 3 m',
    '1161': 'R 3 c',
    '1166': 'R -3 m',
    '1167': 'R -3 c',
    '2005': 'A 1 2 1',
    '2014': 'P 1 21/n 1',
    '2017': 'P 2 21 2',
    '2018': 'P 21 2 21',
    '3018': 'P 2 21 21',
    '4005': 'I 1 2 1',
}

from .. import (
    HKL_data_ABCD,
    HKL_data_D_sigD,
    HKL_data_E_sigE,
    HKL_data_F_phi,
    HKL_data_F_sigF,
    HKL_data_F_sigF_ano,
    HKL_data_I_sigI,
    HKL_data_I_sigI_ano,
    HKL_data_Phi_fom,
    HKL_data_Flag,
    HKL_data_Flag_bool
)

_data_columns_to_data_types = {
    ('F_meas_au', 'F_meas_sigma_au'):                   (float, HKL_data_F_sigF),
    ('F_meas', 'F_meas_sigma'):                         (float, HKL_data_F_sigF),
    ('pdbx_F_plus', 'pdbx_F_plus_sigma',
     'pdbx_F_minus', 'pdbx_F_minus_sigma'):             (float, HKL_data_F_sigF_ano),

    ('pdbx_anom_difference',
     'pdbx_anom_difference_sigma'):                     (float, HKL_data_D_sigD),
    ('pdbx_anomalous_diff',
     'pdbx_anomalous_diff_sigma'):                      (float, HKL_data_D_sigD),

    ('intensity_meas', 'intensity_sigma'):              (float, HKL_data_I_sigI),

    ('pdbx_I_plus', 'pdbx_I_plus_sigma',
      'pdbx_I_minus', 'pdbx_I_minus_sigma'):            (float, HKL_data_I_sigI_ano),

    ('pdbx_HL_A_iso', 'pdbx_HL_B_iso',
     'pdbx_HL_C_iso', 'pdbx_HL_D_iso'):                 (float, HKL_data_ABCD),

    ('phase_calc', 'fom'):                              (float, HKL_data_Phi_fom),

    ('pdbx_r_free_flag'):                               (int, HKL_data_Flag),
    ########
    # CALCULATED STRUCTURE FACTORS
    ########

    ('F_calc', 'phase_calc'):                           (float, HKL_data_F_phi), # Fcalc, phiFcalc
    ('F_calc_au', 'phase_calc'):                        (float, HKL_data_F_phi), # Fcalc, phiFcalc
    ('pdbx_DELFWT', 'pdbx_DELPHWT'):                    (float, HKL_data_F_phi), # mFo-DFc
    ('pdbx_FWT', 'pdbx_PHWT'):                          (float, HKL_data_F_phi), # 2mFo-DFc
    ('pdbx_F_calc_part_solvent',
     'pdbx_phase_calc_part_solvent'):                   (float, HKL_data_F_phi), # Solvent contribution to calculated structure factors
    ('pdbx_F_calc_with_solvent',
     'pdbx_phase_calc_with_solvent'):                   (float, HKL_data_F_phi), # Calculated structure factors including bulk solvent

}

_anomalous_data_columns = (
    ('pdbx_F_plus', 'pdbx_F_plus_sigma', 'pdbx_F_minus', 'pdbx_F_minus_sigma'),
    ('pdbx_I_plus', 'pdbx_I_plus_sigma', 'pdbx_I_minus', 'pdbx_I_minus_sigma')
)


def load_cif_sf(filename, load_map_coeffs=True):
    '''
    Load a set of structure factors from a .cif file.
    '''
    from chimerax.atomic.mmcif import get_cif_tables
    table_list = get_cif_tables(filename, _cif_sf_table_names)
    return _parse_tables(table_list, load_map_coeffs=load_map_coeffs)

def _parse_tables(table_list, load_map_coeffs=True):
    tables = dict(zip(_cif_sf_table_names, table_list))
    metadata = {l: tables[l] for l in _metadata_tables}
    cell_info = tables['cell']
    try:
        cell_dim = [float(d) for d in cell_info.fields(_required_columns['cell'])[0]]
    except:
        raise TypeError('Could not read cell information from file!')
    symm = tables['symmetry']
    for id, dtype in _space_group_identifiers():
        if symm.has_field(id):
            spgr_descriptor = symm.fields((id,))[0][0]
            spgr_dtype = dtype
            from chimerax.clipper.clipper_python import Spgr_descr
            if dtype == Spgr_descr.TYPE.HM:
                spgr_descriptor = _expand_hm_symbol(spgr_descriptor)
            elif dtype == Spgr_descr.TYPE.Number:
                # Check if it's a non-standard CCP4 number which Clipper doesn't recognise
                symbol_hm = _ccp4_symop_code_to_hm_symbol.get(spgr_descriptor, None)
                if symbol_hm is not None:
                    spgr_descriptor = symbol_hm
                    dtype = Spgr_descr.TYPE.HM
            break
    else:
        raise TypeError('Could not read spacegroup information from file!')

    refln_table = tables['refln']
    hkls = _get_miller_indices(refln_table)

    from ..clipper_python import (
        Cell_descr, Cell, Spgr_descr,
        Spacegroup, Resolution, Grid_sampling
    )
    cell = Cell(Cell_descr(*cell_dim))
    spacegroup = Spacegroup(Spgr_descr(spgr_descriptor, spgr_dtype))

    refln_info = tables['reflns']
    # Resolution information is explicitly given in the 'reflns' table in a
    # surprisingly small minority of structure factor cif files, but it's nice when
    # we can get it.
    if refln_info.has_field('d_resolution_high'):
        res = float(refln_info.fields(('d_resolution_high',))[0][0])
    else:
        # If this turns out to be too slow, it could fairly easily be pushed
        # back to C++
        from ..clipper_python import HKL
        invresolsq_limit = max(HKL(hkl).invresolsq(cell) for hkl in hkls)
        from math import sqrt
        res = 1/sqrt(invresolsq_limit)
    print('Resolution: {}'.format(res))
    resolution = Resolution(res)

    from .. import HKL_info
    hkl_info = HKL_info(spacegroup, cell, resolution, True)

    # OK, now we have all our vital information. Time to find all the data
    data = _parse_status(refln_table, hkl_info, hkls)
    exp_names = []
    exp_data = []
    calc_names = []
    calc_data = []
    for column_names, type_spec in _data_columns_to_data_types.items():
        result = _cif_columns_to_clipper(refln_table, hkl_info, hkls, column_names, type_spec)
        if result is not None:
            if isinstance(result, HKL_data_F_phi) and load_map_coeffs:
                calc_names.append(', '.join(column_names))
                calc_data.append(result)
            else:
                exp_names.append(', '.join(column_names))
                exp_data.append(result)
            # data[column_names] = result

    return hkl_info, ('Free flags', data['Free set']), (exp_names, exp_data), (calc_names, calc_data)

    # return cell, spacegroup, resolution, hkl_info, data


def _get_miller_indices(table):
    import numpy
    headers = ('index_h', 'index_k', 'index_l')
    hkls = numpy.array([[int(i) for i in row] for row in table.fields(headers)],
                        dtype=numpy.int32)
    return hkls

def _cif_columns_to_clipper(table, hkl_info, hkls, field_names, type_spec):
    import numpy
    for field in field_names:
        if not table.has_field(field):
            return None
    scalar_type, array_type = type_spec
    if scalar_type == float:
        data = numpy.array(
            [[float(d) if d!='?' else numpy.nan for d in row] for row in table.fields(field_names)],
            numpy.double
        )
    elif scalar_type == int:
        data = numpy.array(
            [[int(d) for d in row] for row in table.fields(field_names)],
            numpy.double
        )
    if field_names in _anomalous_data_columns:
        # Clipper's anomalous data types have a fifth covariance column, but
        # this is not provided in the .cif file. Just set it to 1.
        padded_data = numpy.empty((len(data), 5), numpy.double)
        padded_data[:,:4] = data
        padded_data[:,4] = 1
        data = padded_data
    clipper_data = array_type(hkl_info)
    try:
        clipper_data.set_data(hkls, data)
    except RuntimeError as e:
        err_string = " Field names: {}; HKL array shape: {}; Data array shape: {}".format(
            field_names, hkls.shape, data.shape
        )
        raise RuntimeError(str(e)+err_string)
    return clipper_data

_refln_status_flags = {
    'systematically_absent':    '-',    # Doesn't exist in the given space group
    'unobserved':               '<',    # Too low to measure, but not systematically absent
    'free':                     'f',    # Include in the R-free set
    'cut_high':                 'h',    # Measured, but above the decided high resolution cutoff
    'cut_low':                  'l',    # Measured, but below the decided low resolution cutoff
    'unreliable':               'x',    # Unreliable measurement. Do not use.
    'default':                  'o'     # Normal
}



def _parse_status(table, hkl_info, hkls):
    import numpy
    status = table.fields(('status',))
    n = len(hkls)
    # Due to a quirk of Clipper's template instantiation, all HKL_data types
    # import/export as double regardless of their actual data type.
    free = numpy.zeros((n,1), numpy.double)
    unreliable = numpy.zeros((n,1), numpy.double)
    systematic_absences = numpy.zeros((n,1), numpy.double)
    for i, s in enumerate(status):
        s = s[0]
        if s == 'o':
            continue
        elif s == '-':
            systematic_absences[i] = 1
        elif s == 'f':
            free[i] = 1
        elif s == 'x':
            unreliable[i] = 1
    c_free = HKL_data_Flag(hkl_info)
    c_unreliable = HKL_data_Flag_bool(hkl_info)
    c_systematic_absences = HKL_data_Flag_bool(hkl_info)
    c_free.set_data(hkls, free)
    c_unreliable.set_data(hkls, unreliable)
    c_systematic_absences.set_data(hkls, systematic_absences)
    return {"Free set":     c_free,
            "Unreliable":   c_unreliable,
            "Systematic absences":  c_systematic_absences
            }
