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

from chimerax.core.models import Model

def calculate_voxel_size(resolution, shannon_rate):
    return resolution.limit/2/shannon_rate

def calculate_shannon_rate(resolution, voxel_size):
    return resolution.limit/(2*voxel_size)

class ReflectionDataContainer(Model):
    #SESSION_SAVE=False
    '''
    A container class to hold a set of reciprocal-space data, and
    defining the methods to access and use it. A sub-class of the
    ChimeraX Model class allowing it to be loaded into the model
    hierarchy making it easily visible to the user.
    '''
    def __init__(self, session, hklfile=None, shannon_rate = 2.0,
        free_flag_label = None, auto_choose_reflection_data=True):
        '''
        This class should hold the information that's common to all
        the data contained in its children (e.g. the HKLinfo object,
        the Cell, Spacegroup and Grid_sampling objects, the Unit_Cell,
        etc.
        '''
        Model.__init__(self, 'Reflection Data', session)
        self.filename = None
        self.shannon_rate = shannon_rate
        self._hklinfo = None
        self._grid_sampling = None

        if hklfile is not None:
            self._init_from_hkl_file(hklfile, free_flag_label, auto_choose_reflection_data)


    def _init_from_hkl_file(self, hklfile, free_flag_label, auto_choose_reflection_data):
        import os
        self.filename = os.path.abspath(hklfile)
        hklinfo, free, exp, calc = load_hkl_data(self.session, hklfile,
            free_flag_label=free_flag_label,
            auto_choose_reflection_data=auto_choose_reflection_data)
        self._hklinfo = hklinfo
        self._grid_sampling = None


        if free[0] is not None:
            free_flags = ReflectionDataFreeFlags(free[0], self.session, free[1])
            self.add([free_flags])

        if len(exp[0]):
            dsets = []
            for name, data in zip(*exp):
                    dsets.append(ReflectionDataExp(name, self.session, data))
            self.experimental_data.add(dsets)

        if len(calc[0]):
            dsets = []
            for name, data in zip(*calc):
                    dsets.append(ReflectionDataCalc(name, self.session, data))
            self.calculated_data.add(dsets)

    @property
    def path(self):
        return self.filename

    @property
    def free_flags(self):
        for c in self.child_models():
            if isinstance(c, ReflectionDataFreeFlags):
                return c
        return None

    @property
    def experimental_data(self):
        for c in self.child_models():
            if isinstance(c, ReflectionDataNode) and c.name=='Experimental':
                return c
        ed = ReflectionDataNode('Experimental', self.session)
        self.add([ed])
        return ed

    @property
    def calculated_data(self):
        for c in self.child_models():
            if isinstance(c, ReflectionDataNode) and c.name=='Calculated':
                return c
        cd = ReflectionDataNode('Calculated', self.session)
        self.add([cd])
        return cd

    @property
    def hklinfo(self):
        return self._hklinfo

    @property
    def cell(self):
        return self.hklinfo.cell

    @property
    def spacegroup(self):
        return self.hklinfo.spacegroup

    @property
    def resolution(self):
        return self.hklinfo.resolution

    @property
    def grid_sampling(self):
        if self._grid_sampling is None:
            from . import Grid_sampling
            self._grid_sampling = Grid_sampling(
                self.spacegroup, self.cell, self.resolution, self.shannon_rate)
        return self._grid_sampling

    def take_snapshot(self, session, flags):
        from chimerax.core.models import Model
        data = {
            'model state':          Model.take_snapshot(self, session, flags),
            'original filename':    self.filename,
            'hall symbol':          self.spacegroup.symbol_hall,
            'resolution':           self.resolution.limit,
            'cell dim':             self.cell.dim,
            'cell angles':          self.cell.angles_deg,
            'shannon rate':         self.shannon_rate
        }
        from . import CLIPPER_STATE_VERSION
        data['version']=CLIPPER_STATE_VERSION
        return data

    @staticmethod
    def restore_snapshot(session, data):
        from chimerax.core.models import Model
        rdc = ReflectionDataContainer(session)
        Model.set_state_from_snapshot(rdc, session, data['model state'])
        rdc.filename = data['original filename']
        from chimerax.clipper.clipper_python import (
            Resolution, Spgr_descr, Spacegroup, Cell_descr, Cell, HKL_info
        )
        res = Resolution(data['resolution'])
        cell = Cell(Cell_descr(*data['cell dim'], *data['cell angles']))
        spgr_descr = Spgr_descr(data['hall symbol'], Spgr_descr.Hall)
        rdc._hklinfo = HKL_info(Spacegroup(spgr_descr), cell, res, True)
        return rdc



class ReflectionDataNode(Model):
    '''
    Container class to hold a subset of reflection data within a
    ReflectionDataContainer tree. Typically the subset will be either
    'Experimental' or 'Calculated'.
    '''
    # def __init__(self, name, session):
    #     Model.__init__(self, name, session)
    #     # self.datasets = datasets
    #     for name, data in datasets.items():
    #         self.add([data])
    #SESSION_SAVE=False
    @property
    def datasets(self):
        return dict((m.name, m) for m in self.child_models())

    def __iter__(self):
        return iter(self.datasets.values())

    def __getitem__(self, key):
        return self.datasets[key]

    def take_snapshot(self, session, flags):
        from chimerax.core.models import Model
        data = {
            'model state':          Model.take_snapshot(self, session, flags),
        }
        from . import CLIPPER_STATE_VERSION
        data['version']=CLIPPER_STATE_VERSION
        return data

    @staticmethod
    def restore_snapshot(session, data):
        rdn = ReflectionDataNode('', session)
        Model.set_state_from_snapshot(rdn, session, data['model state'])
        return rdn


class ReflectionData(Model):
    '''
    Prototype for ReflectionDataExp and ReflectionDataCalc. Should
    contain methods common to both (e.g. drawing of reciprocal-space
    reflections).
    '''
    #SESSION_SAVE=False
    def __init__(self, name, session, data):
        '''
        Args:
            name:
                A descriptive name.
            session:
                The ChimeraX session.
            data:
                A Clipper HKL_data_Flag or HKL_data_Flag_bool object.
        '''

        Model.__init__(self, name, session)
        self._data = data

    @property
    def data(self):
        return self._data

    @property
    def dtype(self):
        return type(self.data)

    @property
    def container(self):
        return self.parent.parent

    def take_snapshot(self, session, flags):
        from chimerax.core.models import Model
        data = {
            'container':            self.container,
            'model state':          Model.take_snapshot(self, session, flags),
            'hkl data':             self.data.data,
            'hkl data type':        self.dtype,
        }
        from . import CLIPPER_STATE_VERSION
        data['version']=CLIPPER_STATE_VERSION
        return data

    @classmethod
    def restore_snapshot(cls, session, data):
        from chimerax.core.models import Model
        dtype = data['hkl data type']
        container = data['container']
        hklinfo = container.hklinfo
        hkldata = data['hkl data']
        hklinfo.add_hkl_list(hkldata[0])
        clipper_array = dtype(hklinfo)
        clipper_array.set_data(*hkldata)
        rd = cls('', session, clipper_array)
        Model.set_state_from_snapshot(rd, session, data['model state'])
        return rd






class ReflectionDataFreeFlags(ReflectionData):
    '''Holds the array of free flags.'''

    @property
    def container(self):
        return self.parent


class ReflectionDataExp(ReflectionData):
    '''
    Holds one set of experimental reflection data (e.g. F/sigF, I/sigI),
    etc.
    '''
    pass

class ReflectionDataCalc(ReflectionData):
    '''Holds one set of calculated reflections and phases.'''
    def __init__(self, name, session, data, is_difference_map = None):
        '''
        Args:
            name:
                A descriptive name.
            session:
                The ChimeraX session.
            data:
                A Clipper HKL_data_F_phi object
            is_difference_map(bool):
                If True, maps generated from this data will be displayed
                with both positive and negative contours. If None, the
                code will attempt to guess from the name
        '''
        ReflectionData.__init__(self, name, session, data)
        if is_difference_map is None:
            self.is_difference_map = self._guess_if_difference_map(name)
        else:
            self.is_difference_map = is_difference_map

    def _guess_if_difference_map(self, name):
        first_part = name.split(',')[0]
        identifier = first_part.split(' ')[-1]
        if identifier in ('F','FWT', 'FC', 'FCALC'):
            return False
        if '2' in identifier:
            return False
        return True


def data_labels_match(data_col, sigma_or_phase_col):
  '''
  A quick-and-dirty approach to check if MTZ column names belong together.
  Relies on the fact that under most standard naming the phi or sigma column
  name corresponds to the amplitude/intensity column plus a prefix. We'll
  pair together columns where the second one matches the first plus a prefix
  or a suffix, but not both. We also need to make sure that we don't
  accidentally pair (for example) FOFC with PH2FOFC, so we'll exclude
  any cases where the prefix/suffix contains a digit
  '''
  if data_col in sigma_or_phase_col:
    remainder = [s for s in sigma_or_phase_col.split(data_col) if s]
    if len(remainder) == 1:
      return not any(i.isdigit() for i in remainder[0])
  return False

def find_data_pairs(mtzin, temp_tree, first_type, second_type, data_type):
    '''
    Find all column pairs matching a particular signature and with
    matching labels, create the corresponding Clipper objects, and
    return the objects and their names as two arrays.

    Args:
        mtzin:
            A currently open clipper.CCP4MTZfile object
        temp_tree:
            A DataTree (or similar dict) mirroring the internal data
            structure of the MTZ file.
        first_type (char):
            a character defining the data type for the first column. See:
            http://www.ccp4.ac.uk/html/mtzformat.html#coltypes
            for a list of valid types.
        second_type (char):
            a character defining the data type for the second column.
        data_type:
            a pointer to the clipper class corresponding to this
            data type.
    '''
    data_sets = []
    data_set_names = []
    for crystal in temp_tree.keys():
        for dataset in temp_tree[crystal].keys():
            for iname in temp_tree[crystal][dataset][first_type].keys():
                for sname in temp_tree[crystal][dataset][second_type].keys():
                    if data_labels_match(iname, sname):
                        keyname = ', '.join(map(str, [iname, sname]))
                        this_set = data_type()
                        mtzin.import_hkl_data(this_set,
                            '/'+'/'.join(
                            map(str, [crystal, dataset, '[{}]'.format(keyname)])))
                        data_sets.append(this_set)
                        data_set_names.append(keyname)
    return (data_sets, data_set_names)

def find_free_set(mtzin, temp_tree, label = None):
    '''Find the free set, optionally given an explicit label to look for.'''
    possible_free_flags = []
    possible_free_flags_names = []
    from . import HKL_data_Flag
    for crystal in temp_tree.keys():
        for dataset in temp_tree[crystal].keys():
            # Find and store the free set
            for iname in temp_tree[crystal][dataset]['I'].keys():
                if label is not None:
                    if iname == label:
                        free = HKL_data_Flag()
                        mtzin.import_hkl_data(free,
                            '/'+'/'.join(map(str, [crystal, dataset, '[{}]'.format(iname)])))
                        return ([free],[iname])

                elif 'free' in iname.lower():
                    thisFree = HKL_data_Flag()
                    mtzin.import_hkl_data(thisFree,
                        '/'+'/'.join(map(str, [crystal, dataset, '[{}]'.format(iname)])))
                    possible_free_flags.append(thisFree)
                    possible_free_flags_names.append(iname)

    if label is not None:
        raise TypeError('The label "{}" does not appear to be in this dataset!'.format(label))
    return (possible_free_flags, possible_free_flags_names)


def load_hkl_data(session, filename, free_flag_label = None,
        auto_choose_rfree=True, auto_choose_reflection_data=True,
        load_map_coeffs=True):
    '''
    Load in an mtz file, create Clipper objects from the data, and
    return the tuple:
        ( HKLinfo,
         (free_flags_name, free_flags),
         (experimental_set_names, experimental_sets),
         (calculated_set_names, calculated_sets) )

    where HKLinfo and free_flags are Clipper objects, and
    experimental_sets and calculated_sets are arrays of Clipper objects.

    If free_flag_label is a string, the column with that name will be assigned
    as the free flags.
    '''
    import os
    if os.path.isfile(filename):
        hklfile = os.path.abspath(filename)
    else:
        raise FileNotFoundError('Invalid filename!')

    extension = os.path.splitext(filename)[1].lower()
    if extension not in ('.cif', '.mtz'):
        raise ValueError('Reflection file must be in either .cif or .mtz format!')
        return

    if extension == '.mtz':
        (hklinfo, free, expt, calc) = load_mtz_data(
            session, hklfile, free_flag_label = free_flag_label,
            load_map_coeffs = load_map_coeffs)

    elif extension in ('.cif', '.ent'):
        from .io.cif_sf_read import load_cif_sf
        (hklinfo, free, expt, calc) = load_cif_sf(hklfile,
            load_map_coeffs = load_map_coeffs)

    else:
        from chimerax.core.errors import UserError
        raise UserError('Unrecognised structure factor file format!')

    if free[1] is None:
        from chimerax.clipper import HKL_data_Flag
        free = ('R-free-flags', HKL_data_Flag(hklinfo))

    if len(expt[0]) > 1:
        if auto_choose_reflection_data:
            choice = _auto_choose_reflections(*expt)
            warn_str = ('WARNING: multiple experimental reflection datasets '
                'found:\n {} \n'
                'Automatically choosing "{}".')
            session.logger.warning(warn_str.format(',\n'.join(expt[0]), choice[0][0]))
            expt = choice

    if len(expt[0]):
        _regenerate_free_set_if_necessary(session, free[1], expt[1][0])

    hklinfo, free, expt, calc = _filter_out_missing_free_flags(hklinfo, free, expt, calc)
    return (hklinfo, free, expt, calc)

def _filter_out_missing_free_flags(hklinfo, free, expt, calc):
    ih = hklinfo.first
    flags = free[1]
    good = []
    bad_count=0
    while not ih.last():
        if not flags[ih].missing:
            good.append(ih.hkl.as_numpy())
        else:
            bad_count += 1
        ih.next()
    if bad_count:
        from chimerax.clipper import (HKL_info, HKL_data_Flag, HKL_data_F_sigF, HKL_data_F_phi)
        import numpy
        good = numpy.array(good)
        new_hklinfo = HKL_info(hklinfo.spacegroup, hklinfo.cell, hklinfo.resolution, False)
        new_hklinfo.add_hkl_list(good)
        new_free = (free[0], free[1].restrict_to(new_hklinfo))
        new_expt = [[],[]]
        for i in range(len(expt[0])):
            new_expt[0].append(expt[0][i])
            new_expt[1].append( expt[1][i].restrict_to(new_hklinfo))
        new_calc = [[],[]]
        for i in range(len(calc[0])):
            new_calc[0].append(calc[0][i])
            new_calc[1].append(calc[1][i].restrict_to(new_hklinfo))
        return new_hklinfo, new_free, new_expt, new_calc
    else:
        return hklinfo, free, expt, calc


    



def _auto_choose_reflections(names, datasets):
    '''
    Attempt to automatically choose the best experimental reflection dataset
    from a loaded file to use for generating maps. The logic at present is very
    simple: if a non-anomalous intensity dataset is present, use that. Otherwise,
    if a non-anomalous amplitude dataset is present, use that. Otherwise, move
    on to anomalous intensities, followed by anomalous amplitudes.
    '''
    from chimerax.clipper import (HKL_data_F_sigF, HKL_data_I_sigI,
        HKL_data_F_sigF_ano, HKL_data_I_sigI_ano)
    for dtype in (HKL_data_F_sigF, HKL_data_I_sigI, HKL_data_F_sigF_ano,
            HKL_data_I_sigI_ano
        ):
        for name, dataset in zip(names, datasets):
            if type(dataset) == dtype:
                return ([name], [dataset])
    raise RuntimeError('No suitable experimental data found!')


def _regenerate_free_set_if_necessary(session, flag_array, data_array, free_frac=0.05, max_free = 2000):
    unique_vals = set()
    ih = data_array.first_data
    while not ih.last():
        unique_vals.add(flag_array[ih].data[0])
        data_array.next_data(ih)
    if len(unique_vals)==1:
        from chimerax.clipper.reflection_tools.r_free import generate_free_set
        num_free = generate_free_set(flag_array, data_array, free_frac, max_free)
        warn_str = ('No free flags detected in this dataset! '
            'Automatically generated a new random set with {} free from {} '
            'observed reflections. You should save your data to a new MTZ file '
            'and use this for any future rebuilding/refinement.')
        session.logger.warning(warn_str.format(num_free, data_array.num_obs))


def load_mtz_data(session, filename, free_flag_label = None,
        auto_choose_rfree=True, load_map_coeffs=True):
    from chimerax.clipper import HKL_data_Flag, HKL_data_F_phi
    from .io import mtz_read
    hklinfo, crystal_dict = mtz_read.load_mtz_data(session, filename, load_map_coeffs=load_map_coeffs)
    if len(crystal_dict) == 2:
        crystal_dict = _merge_if_mini_mtz_format(crystal_dict)
    if len(crystal_dict) > 1:
        warn_str = ('WARNING: This MTZ file contains data from multiple crystals. '
            'Only the data from the first crystal will be used. If you wish to '
            'use the other data, please split your MTZ file into individual '
            'datasets (you can do this using tools from the PHENIX or CCP suites).'
        )
        session.logger.warning(warn_str)

    for crystal, datasets in crystal_dict.items():
        break
    all_flag_arrays = {}
    all_exp_data_arrays = {}
    all_calc_data_arrays = {}
    for dataset_name, dataset in datasets.items():
        for data_name, data_array in dataset.items():
            dname_ext = '({}) '.format(dataset_name) + data_name
            if isinstance(data_array, HKL_data_Flag):
                all_flag_arrays[dname_ext] = data_array
            elif isinstance(data_array, HKL_data_F_phi):
                all_calc_data_arrays[dname_ext] = data_array
            else:
                all_exp_data_arrays[dname_ext] = data_array
    flag_names, flag_arrays = (list(all_flag_arrays.keys()), list(all_flag_arrays.values()))
    free_flag_name = None
    if free_flag_label is not None:
        try:
            name_index = [name.lower() for name in flag_names].index(free_flag_label.lower())
            free_flag_name = flag_names[name_index]
        except ValueError:
            from chimerax.core.errors import UserError
            err_str = ('The specified free flag label {} was not found in this '
                'file. Possible names are:\n{}')
            raise UserError(err_str.format(free_flag_label, ',\n'.join(flag_names)))
    free_flags = None
    if free_flag_name is None and len(flag_names):
        if len(flag_names) == 1:
            free_flag_name = flag_names[0]
            free_flags = flag_arrays[0]
        else:
            if auto_choose_rfree:
                for i, name in enumerate(flag_names):
                    if 'free' in name.lower():
                        break
                else:
                    name = flag_names[0]
                warn_str = ('WARNING: found multiple possible R-free arrays: \n {} \n'
                    'Automatically choosing "{}". If this is incorrect, please either '
                    'provide a file with a single flag array, or load the file '
                    'again using the free_flag_label argument.')
                session.logger.warning(warn_str.format(',\n'.join(flag_names), name))
                free_flag_name = name
            else:
                free_flag_name = _r_free_chooser(session, flag_names)
        if free_flag_name is None:
            free_flags = None
        else:
            index = flag_names.index(free_flag_name)
            free_flags = flag_arrays[index]
    exp_data_names, exp_data = (list(all_exp_data_arrays.keys()), list(all_exp_data_arrays.values()))
    calc_data_names, calc_data = (list(all_calc_data_arrays.keys()), list(all_calc_data_arrays.values()))
    return (
        hklinfo,
        (free_flag_name, free_flags),
        (exp_data_names, exp_data),
        (calc_data_names, calc_data)
    )




def _merge_if_mini_mtz_format(crystal_dict):
    '''
    The CCP4-7.0 "mini-MTZ" format is designed such that each MTZ file holds
    only a single data type. That's quite sensible in many respects, but the
    somewhat annoying aspect to its design is that the actual data is placed
    in a different "crystal" path compared to the HKL indices. So we need to
    do a little checking to work out if that's the case here.
    '''
    print('MTZ file either contains multiple crystal datasets or is a mini-MTZ file. Checking...')
    new_dict = {'crystal': {'dataset': {}}}
    inner_dict = new_dict['crystal']['dataset']
    hklbase = crystal_dict.get('HKL_base', None)
    if hklbase is None:
        return crystal_dict
    base_datasets = hklbase.get('HKL_base', None)
    if base_datasets is None:
        return crystal_dict
    # A mini-MTZ will contain only the free set here
    if len(base_datasets) == 1:
        for key, dat in base_datasets.items():
            inner_dict[key] = dat
        cryst_datasets = crystal_dict[list(crystal_dict.keys())[1]]
        if len(cryst_datasets) != 1:
            return crystal_dict
        for dkey, dataset in cryst_datasets.items():
            if len(dataset) != 1:
                return crystal_dict
            for key, data in dataset.items():
                inner_dict[key] = data
        return new_dict
    return crystal_dict



def load_mtz_data_old(session, filename, free_flag_label = None):
    from . import CCP4MTZfile, HKL_info
    mtzin = CCP4MTZfile()
    hkl = HKL_info()
    mtzin.open_read(filename)
    mtzin.import_hkl_info(hkl, True)
    # Get all the column names and types
    column_labels = mtzin.column_paths
    # Sort the columns into groups, and organise into a temporary tree
    from .data_tree import DataTree
    temp_tree = DataTree()['Experiment']
    i = 0
    for l in column_labels:
        thisname, thistype = l.__str__().split(' ')
        crystal, dataset, name = thisname.split('/')[1:]
        # The h, k, l indices are already captured in hklinfo
        if thistype != 'H':
            temp_tree[crystal][dataset][thistype][name] = i
            i += 1

    # FIXME: will need checking here to see if the crystal name already
    # exists, and offer options to replace, rename or append

    # FIXME: ultimately we want to be able to handle any legal MTZ file,
    #        but for now we'll balk if we find more than one crystal
    if len(temp_tree.keys()) != 1:
        errstring =\
        '''
        At present ChimeraX-Clipper cannot handle MTZ files containing
        data from more than one crystal. We hope to fix this in a future
        release, but for now you can split your MTZ file using tools
        available in the Phenix or CCP4 suites. Note that you *can*
        load multiple MTZ files into the one Clipper_MTZ datastructure
        if you wish, using repeat calls to load_hkl_data. However, they
        will all be treated as having the same unit cell parameters. If
        in doubt, create a new Xtal_Project object for each file.
        '''
        raise RuntimeError(errstring)

    # Find experimental data sets
    experimental_sets = []
    experimental_set_names = []


    # Find I/sigI pairs and pull them in as
    # Clipper HKL_data_I_sigI objects
    from . import HKL_data_I_sigI
    from . import HKL_data_F_sigF
    from . import HKL_data_F_phi
    isigi, isigi_names = find_data_pairs(
        mtzin, temp_tree, 'J', 'Q', HKL_data_I_sigI)
    experimental_sets.extend(isigi)
    experimental_set_names.extend(isigi_names)


    # Find F/sigF pairs and pull them in as
    # Clipper HKL_data_F_sigF objects
    fsigf, fsigf_names = find_data_pairs(
        mtzin, temp_tree, 'F', 'Q', HKL_data_F_sigF)
    experimental_sets.extend(fsigf)
    experimental_set_names.extend(fsigf_names)

    calculated_sets = []
    calculated_set_names = []

    # Find amplitude/phase pairs and pull them in as
    # Clipper HKL_data_F_phi objects
    fphi, fphi_names = find_data_pairs(
        mtzin, temp_tree, 'F', 'P', HKL_data_F_phi)
    calculated_sets.extend(fphi)
    calculated_set_names.extend(fphi_names)

    free_flags = None
    free_flags_name = None

    possible_free_flag_names = []
    all_integer_column_names = []
    for crystal in temp_tree.keys():
        for dataset in temp_tree[crystal].keys():
            for iname in temp_tree[crystal][dataset]['I'].keys():
                if 'free' in iname.lower():
                    possible_free_flag_names.append(iname)
                all_integer_column_names.append(iname)

    # possible_free_flags, possible_free_flag_names = find_free_set(
    #     mtzin, temp_tree, free_flag_label)

    if (len(possible_free_flag_names) == 0 and len(experimental_sets)):
        all_integer_column_names = temp_tree[crystal][dataset]['I'].keys()
        if not len(all_integer_column_names) and free_flag_label !=-1:
            err_string = 'This MTZ file does not appear to contain any '\
            + 'columns suitable for use as a free set. Please generate '\
            + 'a suitable set of free flags using your favourite '\
            + 'package (I recommend PHENIX) and try again.'
            raise RuntimeError(err_string)
        else:
            possible_free_flag_names = all_integer_column_names
            if len(possible_free_flag_names) == 1:
                session.logger.info('WARNING: assuming column with label {}'
                + ' defines the free set.'.format(possible_free_flag_names[0]))
                #possible_freefrom ._flags = [temp_tree[crystal][dataset]['I'][possible_free_flag_names[0]]]
    if not possible_free_flag_names:
        free_flags_name = None
    elif len(possible_free_flag_names) > 1:
        free_flags_name = _r_free_chooser(session, possible_free_flag_names)
        if free_flags_name is None:
            if free_flag_label != -1:
                raise RuntimeError('No free flags chosen. Bailing out.')
    else:
        free_flags_name = possible_free_flag_names[0]

    if free_flags_name:
        free_flags = find_free_set(mtzin, temp_tree, label=free_flags_name)[0][0]
    else:
        free_flags = None

    mtzin.close_read()

    return ( hkl,
            (free_flags_name, free_flags),
            (experimental_set_names, experimental_sets),
            (calculated_set_names, calculated_sets) )

def _r_free_chooser(session, possible_names):
    from PyQt5.QtWidgets import QInputDialog
    choice, ok_pressed = QInputDialog.getItem(session.ui.main_window, 'Choose R-free column', 'Label: ', possible_names, 0, False)
    if ok_pressed and choice:
        return choice
    return None
