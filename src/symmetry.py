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
import sys, os, glob
import ctypes

from . import clipper_python

from chimerax.atomic import molc, structure
# from chimerax.atomic.molc import CFunctions, string, cptr, pyobject, \
#     set_c_pointer, pointer, size_t

CFunctions = molc.CFunctions
string = molc.string
cptr = molc.cptr
pyobject = molc.pyobject
set_c_pointer = molc.set_c_pointer
pointer = molc.pointer
size_t = molc.size_t

dpath = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(dpath, 'lib','lib_symmetry*'))[0]

_c_functions = CFunctions(os.path.splitext(libfile)[0])

_symmetry = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), libfile))
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function

HIDE_ISOLDE = 0x02

auto_reset_camera=True

BACKBONE_MODE_RIBBON=0
BACKBONE_MODE_CA_TRACE=1

_backbone_mode_descr = {
    BACKBONE_MODE_RIBBON: "ribbon",
    BACKBONE_MODE_CA_TRACE: "CA trace"
}

def space_group_hm_synonym(symbol_hm):
    synonyms = {
        'H 3':  'R 3 :H',
        'H -3': 'R -3 :H',
        'H 3 2': 'R 3 2 :H',
    }
    return synonyms.get(symbol_hm)

def _format_sym_tuple(result):
    from chimerax.atomic import ctypes_support as convert
    from chimerax.geometry import Places
    primary_atoms = convert.atoms(result[0])
    sym_atoms = convert.atoms(result[1])
    n_sym_atoms = len(sym_atoms)
    sym_coords = result[2].reshape((n_sym_atoms,3))
    atom_syms = result[3]
    sym_bonds = convert.bonds(result[4])
    nbonds = len(sym_bonds)
    bond_positions = Places(opengl_array=result[5].reshape((nbonds*2,4,4)))
    bond_syms = result[6]
    return (primary_atoms, sym_atoms, sym_coords, atom_syms, sym_bonds, bond_positions, bond_syms)


def sym_transforms_in_sphere(atoms, transforms, center, cutoff, visible_only = True):
    f = c_function('atom_and_bond_sym_transforms',
        args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t, ctypes.POINTER(ctypes.c_double), ctypes.c_double,
            ctypes.c_bool),
        ret = ctypes.py_object)
    natoms = len(atoms)
    tf = numpy.empty(transforms.shape, numpy.double)
    tf[:] = transforms
    c = numpy.empty(3, numpy.double)
    c[:] = center
    n_tf = len(transforms)
    result = f(atoms._c_pointers, natoms, pointer(tf), n_tf, pointer(c), cutoff, visible_only)
    return _format_sym_tuple(result)

def whole_residue_sym_sphere(residues, transforms, center, cutoff, visible_only = True):
    f = c_function('atom_and_bond_sym_transforms_by_residue',
        args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t, ctypes.POINTER(ctypes.c_double), ctypes.c_double,
            ctypes.c_bool),
            ret = ctypes.py_object)
    nres = len(residues)
    tf = numpy.empty(transforms.shape,numpy.double)
    tf[:] = transforms
    c = numpy.empty(3, numpy.double)
    c[:] = center
    n_tf = len(transforms)
    result = f(residues._c_pointers, nres, pointer(tf), n_tf, pointer(c), cutoff, visible_only)
    return _format_sym_tuple(result)

def atom_and_bond_sym_transforms_from_sym_atoms(atoms, symops, sym_indices, visible_only = True):
    f = c_function('atom_and_bond_sym_transforms_from_sym_atoms',
        args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint8),
            ctypes.c_size_t, ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t, ctypes.c_bool),
        ret=ctypes.py_object)
    natoms = len(atoms)
    nsym = len(sym_indices)
    tf = numpy.empty(symops.shape, numpy.double)
    tf[:] = symops
    si = numpy.empty(len(sym_indices), numpy.uint8)
    si[:] = sym_indices
    result = f(atoms._c_pointers, pointer(sym_indices), natoms, pointer(tf), nsym, visible_only)
    return _format_sym_tuple(result)


def sym_ribbons_in_sphere(tether_coords, transforms, center, cutoff):
    f = c_function('close_sym_ribbon_transforms',
        args=(ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
              ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
              ctypes.POINTER(ctypes.c_double), ctypes.c_double),
        ret=ctypes.py_object)
    tc = numpy.empty(tether_coords.shape, numpy.double)
    tc[:] = tether_coords
    tf = numpy.empty(transforms.shape, numpy.double)
    tf[:] = transforms
    c = numpy.empty(3, numpy.double)
    c[:] = center
    result = f(pointer(tc), len(tc), pointer(tf), len(tf), pointer(c), cutoff)
    return result

def symmetry_coords(atoms, sym_matrices, sym_indices):
    unique_indices = numpy.unique(sym_indices)
    coords = atoms.coords
    from chimerax.geometry import Place
    for i in unique_indices:
        mask = sym_indices == i
        tf = Place(matrix=sym_matrices[i])
        coords[mask] = tf.transform_points(coords[mask])
    return coords

from chimerax.core.models import Model, Drawing

def get_symmetry_handler(structure, create=False, auto_add_to_session=False):
    sh = _get_symmetry_handler(structure, create)
    if auto_add_to_session and sh is not None:
        if sh not in structure.session.models.list():
            structure.session.models.add([sh])
    return sh

def _get_symmetry_handler(structure, create=False):
    p = structure.parent
    if isinstance(p, SymmetryManager):
            return p
    if create:
        return SymmetryManager(structure.session, model=structure)
    return None

def get_all_symmetry_handlers(session):
    from chimerax.atomic import AtomicStructure
    handlers = []
    for m in session.models.list(type=AtomicStructure):
        sh = get_symmetry_handler(m, create=False)
        if sh is not None:
            handlers.append(sh)
    return handlers

def get_map_mgr(structure, create=False, auto_add_to_session=False):
    sh = get_symmetry_handler(structure, create=create, auto_add_to_session=auto_add_to_session)
    if sh is not None:
        return sh.map_mgr
    return None

def is_crystal_map(volume):
    from .maps.map_handler_base import XmapHandlerBase
    return isinstance(volume, XmapHandlerBase)

def is_managed_map(volume):
    from .maps.map_handler_base import MapHandlerBase
    return isinstance(volume, MapHandlerBase)

def symmetry_from_restore_data(model, data):
    from chimerax.clipper import (
        Resolution, Spgr_descr, Spacegroup, Cell_descr, Cell, Grid_sampling
    )
    has_symmetry = data['has symmetry']
    if not has_symmetry:
        return simple_p1_box(model)
    res = Resolution(data['resolution'])
    cell = Cell(Cell_descr(*data['cell dim'], *data['cell angles']))
    spgr_descr = Spgr_descr(data['hall symbol'], Spgr_descr.Hall)
    spacegroup = Spacegroup(spgr_descr)   
    grid_sampling = Grid_sampling(spacegroup, cell, res)
    return cell, spacegroup, grid_sampling, res, True 


def symmetry_from_model_metadata(model):
    '''
    Create crystallographic symmetry objects based on the information in the
    model's metadata (PDB CRYST1 card or mmCIF header info)
    '''
    if 'CRYST1' in model.metadata.keys():
        return symmetry_from_model_metadata_pdb(model)
    elif getattr(model, 'is_corecif', False):
        return symmetry_from_corecif(model)
    elif 'cell' in model.metadata.keys():
        return symmetry_from_model_metadata_mmcif(model)
    return simple_p1_box(model)
    # raise TypeError('Model does not appear to have symmetry information!')


def symmetry_from_corecif(model):
    '''
    Crystallographic symmetry for a small-molecule (corecif) model. corecif
    consumes the unit cell to place atoms and stores neither the cell nor the
    symmetry in the model metadata, so both are recovered from the source CIF
    (model.filename), building the space group from the explicit list of symmetry
    operators. Returns the 5-tuple (cell, spacegroup, grid, resolution,
    has_symmetry) expected by SymmetryManager.add_model.
    '''
    path = getattr(model, 'filename', None)
    if path is None:
        return simple_p1_box(model)
    try:
        cell, spacegroup, grid, res = _crystal_symmetry_from_cif(path)
    except Exception:
        return simple_p1_box(model)
    from .clipper_python import Resolution
    return cell, spacegroup, grid, Resolution(res), True


# Small-molecule CIF numbers commonly carry their estimated standard uncertainty in
# parentheses, e.g. "19.4943(15)"; strip it (and any CIF quoting) before float().
def _strip_su(value):
    value = value.strip().strip("'\"")
    paren = value.find('(')
    if paren != -1:
        value = value[:paren]
    return value


def _first_cif_field(table, field):
    if table is not None and table.has_field(field):
        v = table.fields((field,))[0][0].strip()
        if v not in ('', '?', '.'):
            return v
    return None


def crystal_symmetry_from_cif_file(path):
    '''Build Clipper Cell, Spacegroup and Grid_sampling from a (small-molecule /
    core) CIF, with the space group built from its explicit symmetry operators.'''
    cell, spacegroup, grid, res = _crystal_symmetry_from_cif(path)
    return cell, spacegroup, grid


def _crystal_symmetry_from_cif(path):
    from chimerax.mmcif import get_cif_tables
    from .clipper_python import (Cell_descr, Cell, Spgr_descr, Spacegroup,
        Resolution, Grid_sampling)
    cell_t, symm_t, sg_t, sgs_t, diffrn_t = get_cif_tables(
        path, ['cell', 'symmetry', 'space_group', 'space_group_symop', 'diffrn'])

    names = ('length_a', 'length_b', 'length_c',
             'angle_alpha', 'angle_beta', 'angle_gamma')
    try:
        row = cell_t.fields(names)[0]
        params = [float(_strip_su(v)) for v in row]
    except Exception:
        raise TypeError('Could not read unit cell from CIF!')
    cell = Cell(Cell_descr(*params))

    # Prefer the explicit symop list (covers non-standard settings and the full
    # range of centrosymmetric / glide / mirror groups); fall back to number/name.
    ops = []
    for table, field in ((sgs_t, 'operation_xyz'), (symm_t, 'equiv_pos_as_xyz')):
        if table is not None and table.has_field(field):
            rows = table.fields((field,))
            ops = [r[0].strip().strip("'\"").replace(' ', '') for r in rows]
            ops = [o for o in ops if o]
            if ops:
                break
    if ops:
        spacegroup = Spacegroup(Spgr_descr(';'.join(ops), Spgr_descr.TYPE.Symops))
    elif sg_t is not None and sg_t.has_field('IT_number') and \
            _first_cif_field(sg_t, 'IT_number'):
        spacegroup = Spacegroup(Spgr_descr(int(_first_cif_field(sg_t, 'IT_number'))))
    else:
        sym = (_first_cif_field(symm_t, 'space_group_name_Hall')
               or _first_cif_field(symm_t, 'space_group_name_H-M'))
        if sym is None:
            raise TypeError('Could not determine space group from CIF!')
        spacegroup = Spacegroup(Spgr_descr(sym))

    res = _resolution_from_cif(diffrn_t, cell_t)
    grid = Grid_sampling(spacegroup, cell, Resolution(res))
    return cell, spacegroup, grid, res


def _resolution_from_cif(diffrn_t, cell_t):
    # High-resolution limit d_min = lambda / (2 sin(theta_max)); theta in degrees.
    #
    # This resolution only sizes the Grid_sampling for the real-space map and the
    # special-position multiplicity lookup - it is NOT the resolution of the
    # structure-factor calculation (that comes from the reflections themselves). It
    # must nonetheless be a sane positive number: COD is machine-aggregated, so the
    # CIF fields are adversarial. A non-positive wavelength (the '-1'/'0' "not
    # recorded" sentinels) yields a non-positive d_min -> Grid_sampling returns a
    # negative/overflowed grid -> the Xmap allocation SEGFAULTS; an electron-
    # diffraction wavelength (~0.02 A) fed to the X-ray Bragg relation yields an
    # absurd ~0.016 A d_min -> a multi-billion-point grid -> OUT-OF-MEMORY crash.
    # Reject non-positive inputs to the default and clamp to a physical small-molecule
    # range so neither can happen.
    from math import sin, radians, isfinite
    default = 0.84  # sensible small-molecule default when collection limits are absent
    wl = _first_cif_field(diffrn_t, 'radiation_wavelength')
    theta = _first_cif_field(diffrn_t, 'reflns_theta_max')
    if theta is None:
        theta = _first_cif_field(cell_t, 'measurement_theta_max')
    res = default
    try:
        if wl is not None and theta is not None:
            wl = float(_strip_su(wl))
            sin_theta = sin(radians(float(_strip_su(theta))))
            if wl > 0 and sin_theta > 0:
                res = wl / (2.0 * sin_theta)
    except (ValueError, ZeroDivisionError):
        res = default
    if not isfinite(res) or res <= 0:
        res = default
    return min(max(res, 0.5), 5.0)

def simple_p1_box(model, resolution=3):
    '''
    If the model lacks crystallographic information or the available symmetry
    information is incompatible with a real crystal structure, just create
    a large, rectangular P1 box.
    '''
    coord_range = model.atoms.coords.max(axis=0) - model.atoms.coords.min(axis=0)
    abc = coord_range*3
    angles = [90,90,90]
    sym_str = 'P 1'
    from .clipper_python import (Cell_descr, Cell, Spgr_descr, Spacegroup,
        Resolution, Grid_sampling)
    cell = Cell(Cell_descr(*abc, *angles))
    spacegroup = Spacegroup(Spgr_descr(sym_str, Spgr_descr.TYPE.HM))
    res = Resolution(resolution)
    grid_sampling=Grid_sampling(spacegroup, cell, res)
    return cell, spacegroup, grid_sampling, res, False


def symmetry_from_model_metadata_mmcif(model):
    metadata = model.metadata
    res = None

    if 'em_3d_reconstruction' in metadata.keys():
        em_dict = dict((key.lower(), data.lower()) for (key,data) in zip(
            metadata['em_3d_reconstruction'][1:], metadata['em_3d_reconstruction data']
        ))
        res_str = em_dict.get('resolution', '?')
        if res_str != '?':
            try:
                res = float(res_str)
            except:
                pass

    if not res and 'refine' in metadata.keys():
        refine_dict = dict((key.lower(), data.lower()) for (key,data) in zip(
            metadata['refine'][1:], metadata['refine data']
        ))
        res_str = refine_dict.get('ls_d_res_high', '?')
        if res_str != '?':
            try:
                res = float(res_str)
            except:
                pass
    if not res:
            res = 3.0

    # Extract the real cell + space group only for crystallographic experiments.
    # Electron *crystallography* (micro-ED / SerialED / 3D-ED) and neutron
    # diffraction are genuine crystals with a real unit cell, so include them
    # alongside X-ray; single-particle cryo-EM (em_3d_reconstruction) and other
    # non-crystalline methods still fall back to a P1 bounding box.
    CRYSTALLOGRAPHIC_METHODS = (
        'X-RAY DIFFRACTION', 'ELECTRON CRYSTALLOGRAPHY', 'NEUTRON DIFFRACTION')
    exptl_data = metadata.get('exptl data', None)
    if exptl_data is None or not any(m in exptl_data for m in CRYSTALLOGRAPHIC_METHODS):
        return simple_p1_box(model, resolution=res)

    try:
        cell_headers = metadata['cell'][1:]
        cell_data = metadata['cell data']
        cell_dict = dict((key.lower(), data.lower()) for (key, data) in zip(cell_headers, cell_data))
        abc = [float(f) for f  in [cell_dict['length_a'], cell_dict['length_b'], cell_dict['length_c']]]
        angles = [float(f) for f in [cell_dict['angle_alpha'], cell_dict['angle_beta'], cell_dict['angle_gamma']]]
    except:
        raise TypeError('No cell information available!')

    try:
        spgr_headers = metadata['symmetry'][1:]
        spgr_data = metadata['symmetry data']
    except:
        raise TypeError('No space group headers in metadata!')

    from .clipper_python import Cell_descr, Cell, Spgr_descr, Spacegroup, Resolution, Grid_sampling

    spgr_dict = dict((key.lower(), data.lower()) for (key, data) in zip(spgr_headers, spgr_data))
    spgr_str = spgr_dict['space_group_name_h-m'].capitalize()
    if spgr_str != '?':
        if spgr_str.startswith('H'):
            xspgr_str = space_group_hm_synonym(spgr_str)
            if xspgr_str is None:
                print('WARNING: Hermann-Maguin space group symbol "{}" was not recognised. Ignoring.'.format(spgr_str))
                spgr_str = '?'
            else:
                spgr_str = xspgr_str
        if ':' in spgr_str:
            spgr_dtype = Spgr_descr.TYPE.XHM
        else:
            spgr_dtype = Spgr_descr.TYPE.HM
    if spgr_str =='?':
        spgr_str = spgr_dict['space_group_name_hall'].capitalize()
        if spgr_str != '?':
            spgr_dtype = Spgr_descr.TYPE.Hall

    if spgr_str == '?':
        spgr_str = spgr_dict['int_tables_number']
        if spgr_str != '?':
            spgr_dtype = Spgr_descr.TYPE.Number
        else:
            raise TypeError('No space group information available!')


    cell_descr = Cell_descr(*abc, *angles)
    cell = Cell(cell_descr)
    try:
        spgr_descr = Spgr_descr(spgr_str, spgr_dtype)
    except RuntimeError:
        raise RuntimeError('No spacegroup of name {}!'.format(spgr_str))
    spacegroup = Spacegroup(spgr_descr)
    resolution = Resolution(res)
    grid_sampling = Grid_sampling(spacegroup, cell, resolution)
    return cell, spacegroup, grid_sampling, resolution, True

def parse_symops_from_pdb_header(remarks):
    # '''
    # The spacegroup identifier tends to be the most unreliable part
    # of the CRYST1 card, so it's considered safer to let Clipper
    # infer it from the list of symmetry operators at remark 290. This
    # typically looks something like the following:
    #
    # REMARK 290      SYMOP   SYMMETRY
    # REMARK 290     NNNMMM   OPERATOR
    # REMARK 290       1555   X,Y,Z
    # REMARK 290       2555   -X,-Y,Z+1/2
    # REMARK 290       3555   -Y+1/2,X+1/2,Z+1/4
    # REMARK 290       4555   Y+1/2,-X+1/2,Z+3/4
    # REMARK 290       5555   -X+1/2,Y+1/2,-Z+1/4
    # REMARK 290       6555   X+1/2,-Y+1/2,-Z+3/4
    # REMARK 290       7555   Y,X,-Z
    # REMARK 290       8555   -Y,-X,-Z+1/2
    #
    # Clipper is able to initialise a Spacegroup object from a
    # string containing a semicolon-delimited list of the symop
    # descriptors in the SYMMETRY OPERATOR column, so we need to
    # parse those out.
    # '''
    # # Find the start of the REMARK 290 section
    for i, remark in enumerate(remarks):
        if 'REMARK 290' in remark:
            break;
    else:
        return None
    while 'NNNMMM' not in remarks[i]:
        i += 1
    i+= 1
    thisline = remarks[i]
    symstrings = []
    while 'REMARK 290' in thisline and 'X' in thisline:
        symstrings.append(thisline.split()[3])
        i+= 1
        thisline = remarks[i]
    return ';'.join(symstrings)

def symmetry_from_model_metadata_pdb(model):
    '''
    Generate Cell, Spacegroup and a default Grid_Sampling from the PDB
    CRYST1 card.
    '''
    logger = model.session.logger
    from chimerax.clipper.clipper_python import Spgr_descr
    symstr = None
    symstr_type = None
    metadata = model.metadata
    remarks = metadata.get('REMARK', None)
    if remarks is not None and len(remarks):
        try:
            symstr = parse_symops_from_pdb_header(remarks)
            symstr_type = Spgr_descr.TYPE.Symops
        except:
            symstr = None
    cryst1 = metadata.get('CRYST1', None)
    if cryst1 is not None:
        cryst1 = cryst1[0]
        try:
            abc = [float(cryst1[6:15]), float(cryst1[15:24]), float(cryst1[24:33])]
            angles = [float(cryst1[33:40]), float(cryst1[40:47]), float(cryst1[47:54])]
            if symstr is None:
                symstr = cryst1[55:66].capitalize().strip()
                if symstr.startswith('H'):
                    x_symstr = space_group_hm_synonym(symstr)
                    if x_symstr is None:
                        err_str = ('Unrecognised space group symbol: {}. '
                            'If you believe this to be in error, please report '
                            'this as a bug.')
                        raise TypeError(err_str.format(symstr))
                    symstr = x_symstr
                    symstr_type = Spgr_descr.TYPE.XHM
                else:
                    if ':' in symstr:
                        symstr_type = Spgr_descr.TYPE.XHM
                    else:
                        symstr_type = Spgr_descr.TYPE.HM
        except:
            symstr = None
    if symstr is None:
        logger.warning('Missing or corrupted symmetry information in the PDB file. '
            'This model will be treated as a cryo-EM model until associated '
            'with an MTZ file containing symmetry information.')
        return simple_p1_box(model)
    # zval = int(cryst1[67:71])

    res = 3.0
    try:
        i = 0

        '''
        Get the resolution. We need this to define a Grid_sampling
        for the unit cell (needed even in the absence of a map since
        atomic symmetry lookups are done with integerised symops for
        performance). We want to be as forgiving as possible at this
        stage - we'll use the official resolution if we find it, and
        set a default resolution if we don't. This will be overridden
        by the value from any mtz file that's loaded later.
        '''
        while 'REMARK   2' not in remarks[i]:
            i += 1
        # The first 'REMARK   2' line is blank by convention, and
        # resolution is on the second line
        i += 1
        line = remarks[i].split()
        if line[2] == 'RESOLUTION':
            res = float(line[3])
    except:
        pass

    # If this is an EM structure, the CRYST1 card should look something like:
    # CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1
    # If this is the case, just put the model into a big P1 box.

    import numpy
    if symstr.strip() == 'P 1' and numpy.allclose(abc, 1.0) and numpy.allclose(angles, 90):
        return simple_p1_box(model, resolution=res)

    # It turns out that at least one cryo-EM package writes the cell dimensions
    # as 0.000    0.000    0.000, leading to a divide-by-zero hard crash in
    # Clipper. Let's do a sanity check for absurdly short cell lengths here.

    if any(l < 1 for l in abc):
        warn_str = ('CRYST1 reports physically impossible cell dimensions of '
            '({})Å³. Symmetry information in PDB file will be ignored.').format(
                ' x '.join('{:.2f}'.format(l) for l in abc)
                )
        logger.warning(warn_str)
        return simple_p1_box(model)

    from .clipper_python import Cell_descr, Cell, Spgr_descr, Spacegroup, Resolution, Grid_sampling

    cell_descr = Cell_descr(*abc, *angles)
    cell = Cell(cell_descr)
    try:
        spgr_descr = Spgr_descr(symstr, symstr_type)
    except RuntimeError as e:
        from chimerax.core.errors import UserError
        raise UserError(f'Unrecognised space group string: {symstr}')
    spacegroup = Spacegroup(spgr_descr)
    resolution = Resolution(float(res))
    grid_sampling = Grid_sampling(spacegroup, cell, resolution)
    return cell, spacegroup, grid_sampling, resolution, True

def apply_scene_positions(model):
    '''
    If a model has a non-identity transform, apply it to the atomic coordinates
    and reset the transform to identity.
    '''
    p = model.scene_position
    if not p.is_identity():
        model.atoms.transform(p)
        # model.atoms.coords = p*model.atoms.coords
        from chimerax.geometry import Place
        model.position = Place()
        if isinstance(model.parent, SymmetryManager):
            model.parent.position = Place()



class Unit_Cell(clipper_python.Unit_Cell):
   def __init__(self, atoms, cell,
                 spacegroup, grid_sampling, padding = 10):
        '''
        __init__(self, ref, atom_list, cell, spacegroup, grid_sampling) -> Unit_Cell

        Note: internal calculations for finding symmetry equivalents are
        run using integer symmetry operations on grid coordinates for
        improved performance. This means that if you change the sampling
        rate of your maps, you will need to create a new Unit_Cell object
        to match.

        Args:
            atom_list (clipper.Atom_list):
                A Clipper Atom_list object containing your reference model.
            cell (clipper.Cell)
            spacegroup (clipper.Spacegroup)
            grid_sampling (clipper.Grid_Sampling)
            padding (int):
                An optional extra padding (in number of grid steps) to
                add around the reference model when defining the reference
                box for finding symmetry operations. In most cases this
                can be left as zero.
        '''
        from .util import atom_list_from_sel
        atom_list = atom_list_from_sel(atoms)
        super().__init__(atom_list, cell, spacegroup, grid_sampling, padding)

   def update_reference(self, atoms):
        '''
        Recompute the reference-model bounding box and symmetry-search caches from
        the current coordinates of `atoms` (a ChimeraX Atoms selection), reusing
        the padding supplied at construction. Cheap (dominated by the number of
        symmetry operators); call after the model is moved far enough that the
        cached reference no longer represents it (e.g. after shifting a fragment to
        a different ASU).
        '''
        # Pass raw coordinates rather than building a Clipper Atom_list: only the
        # bounding box and centroid are needed, so constructing N Atom objects
        # (element strings, occupancies, ADPs) would be wasted work - significant
        # at ribosome scale.
        super().update_reference_from_coords(atoms.coords)


class SymmetryManager(Model):
    '''
    Handles crystallographic symmetry and maps for an atomic model.
    '''
    # Shift of the centre of rotation required to trigger an update of the
    # spotlight
    SPOTLIGHT_UPDATE_THRESHOLD = 0.1
    ATOMIC_SYM_EXTRA_RADIUS = 3
    # SESSION_ENDURING=False
    # SESSION_SAVE=False
    def __init__(self, session, model=None, mtzfile=None,
        spotlight_mode = True, spotlight_radius=12,
        hydrogens='polar', ignore_model_symmetry=False,
        set_lighting_to_simple=True, debug=False):
        if model is not None:
            if isinstance(model.parent, SymmetryManager):
                raise RuntimeError('This model already has a symmetry manager!')
            name = model.name
        else:
            name = 'N/A'
        name = 'Data manager ({})'.format(name)
        #session = model.session
        if set_lighting_to_simple:
            from chimerax.std_commands import lighting
            lighting.lighting(session, preset='simple')


        super().__init__(name, session)
        self._last_box_center = session.view.center_of_rotation
        self._session_restore = False
        # Durable "just restored from a session" flag. Unlike _session_restore
        # (True only during the synchronous add_model in _end_restore_session_cb),
        # this stays True until after the deferred 'frame drawn' default-styling
        # callbacks have run, so restored per-atom colours / draw modes / display
        # bits / radii / cartoon are NOT overwritten with Clipper defaults.
        self._restored_skip_default_styling = False
        self._debug = debug
        self._stepper = None
        self._last_covered_selection = None
        self._hydrogen_mode = hydrogens

        if not hasattr(self, 'triggers'):
            from chimerax.core.triggerset import TriggerSet
            self.triggers = TriggerSet()

        trigger_names = (
            # Fires when SymmetryManager is in spotlight mode and the centre of
            # rotation moves by more than the minimum distance from its previous
            # location
            'spotlight moved',

            # Fires when spotlight mode turned on/off
            'mode changed',

            # Fires whenever the size and/or location of the map box needs to
            # change - i.e. when the spotlight radius is changed, or when
            # expanding to cover a specific selection.
            'map box changed',
            # Like 'map box changed' but for atomic symmetry (box parameters
            # may be different)
            'atom box changed',
            # Changed backbone mode from ribbon to CA trace or vice versa
            'backbone mode changed',
            # Simply re-fires the model's 'changes' trigger
            'atoms changed',
            # Replaced the atomic model with a new one. All structure change
            # handlers will need to be reapplied
            'model replaced',
            # Shifted some subset of atoms into a different ASU. Firing data
            # should contain the shifted atoms and the transform applied
            'sym shifted atoms',
            # Fires whenever the displayed set of symmetry ghosts is (re)drawn -
            # i.e. at the end of AtomicSymmetryModel.update_graphics (spotlight
            # move, focal-set display, coordinate change, ...). Firing data is the
            # AtomicSymmetryModel. Lets dependent bundles (e.g. ISOLDE, which
            # draws its own validation/restraint markup on the ghosts) refresh in
            # lockstep without polling. Registered on the manager (not the lazily-
            # recreated AtomicSymmetryModel) so subscriptions are durable.
            'symmetry display changed',
        )
        for t in trigger_names:
            self.triggers.add_trigger(t)




        self._update_handler = None
        self._structure_remove_handler = None
        self._swapping_model = False
        self.initialized=False
        if model is not None:
            self.add_model(model, ignore_model_symmetry=ignore_model_symmetry,
                spotlight_mode=spotlight_mode)

        if mtzfile is not None:
            if model is None:
                raise RuntimeError('If providing a structure factor file during '
                    'initialisation, an atomic model must also be provided.')
            mmgr = self.map_mgr
            mmgr.add_xmapset_from_mtz(mtzfile)


    def add_model(self, model, ignore_model_symmetry = False, spotlight_mode=True,
            spotlight_radius=12, session_restore_data = None, debug=False):
        self._structure_change_handler = model.triggers.add_handler(
            'changes', self._structure_change_cb)
        if self.structure is not None and not self._session_restore:
            raise RuntimeError('This SymmetryManager already has an atomic '
                'structure associated with it! To swap for a new one, use '
                'swap_model() instead.')
        self.name = 'Data manager ({})'.format(model.name)
        self._last_box_center = model.atoms.coords.mean(axis=0)

        if session_restore_data is not None:
            f = symmetry_from_restore_data
            args = [model, session_restore_data]
        elif ignore_model_symmetry:
            f = simple_p1_box
            args = []
        else:
            f = symmetry_from_model_metadata
            args=[model]
        self.cell, self.spacegroup, self.grid, self.resolution, self._has_symmetry = f(*args)
        mmgr = self.map_mgr
        uc = self._unit_cell = Unit_Cell(model.atoms, self.cell, self.spacegroup, self.grid)
        self.spotlight_mode = spotlight_mode
        # self._atomic_symmetry_model = AtomicSymmetryModel(self,
        #     radius = spotlight_radius, live = spotlight_mode, debug=debug)
        self.spotlight_radius=spotlight_radius
        if self._update_handler is None:
            self._update_handler = self.session.triggers.add_handler('new frame',
                self.update)
        id = model.id
        self._anisou_sanity_check(model.atoms)
        if model.parent != self:
            apply_scene_positions(model)
            self._transplant_model(model)
        # "Touch" the atomic_symmetry_model to ensure it's created
        self.atomic_symmetry_model

        self.initialized=True
        if self._structure_remove_handler is None:
            self._structure_remove_handler = self.session.triggers.add_handler(
                'remove models', self._structure_remove_cb
            )
        if not self._session_restore:
            self.set_default_atom_display(mode=self._hydrogen_mode)
        # NB: no metadata write-back here. On plain association Clipper *reads* the
        # model's symmetry (from its metadata), so there is nothing to rewrite, and
        # doing so would clobber fields like _exptl.method. Write-back happens only
        # where Clipper *changes* the symmetry: add_symmetry_info (new/different
        # dataset) and discard_model_symmetry (revert to P1).

    def _structure_change_cb(self, trigger_name, changes):
        if 'aniso_u changed' in changes[1].atom_reasons():
            self._anisou_sanity_check(changes[1].modified_atoms())
        from chimerax.atomic import get_triggers
        self._current_changes = changes
        get_triggers().add_handler('changes done', self._changes_done_cb)

    def _changes_done_cb(self, *_):
        self.triggers.activate_trigger('atoms changed', self._current_changes)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


    def _anisou_sanity_check(self, atoms):
        from chimerax.clipper.sanity_check import remove_invalid_anisou
        remove_invalid_anisou(self.session, atoms)

    def add_symmetry_info(self, cell, spacegroup, grid_sampling, resolution, override=False):
        if self.has_symmetry and not override:
            from chimerax.core.errors import UserError
            raise UserError('This manager already has crystallographic symmetry information! '
                "If you're sure, you can force Clipper to use the new symmetry by setting override=True")
        self.cell = cell
        self.spacegroup = spacegroup
        self.grid = grid_sampling
        self.resolution = resolution
        self._unit_cell = Unit_Cell(self.structure.atoms, cell, spacegroup, grid_sampling)
        self.has_symmetry = True
        self._write_symmetry_to_model_metadata()

    @staticmethod
    def _metadata_field(md, category, tag):
        '''Value of one field in an mmCIF metadata category, or None. ChimeraX
        stores a category as md[category] = [category, tag1, ...] and
        md[category+' data'] = [val1, ...] aligned to the tags.'''
        headers = md.get(category)
        values = md.get(category + ' data')
        if not headers or not values:
            return None
        try:
            idx = [h.lower() for h in headers[1:]].index(tag.lower())
        except ValueError:
            return None
        return values[idx] if idx < len(values) else None

    def _write_symmetry_to_model_metadata(self):
        '''Rewrite the model's *native* symmetry metadata to Clipper's current
        authoritative symmetry, so the model stays self-consistent and any later
        re-derivation (or re-save to PDB/mmCIF) is correct. This matters when a
        model's original metadata symmetry is irrelevant to its current use -- an
        old crystal model reused in a cryo-EM (non-crystallographic) context, or
        reused with a new crystal dataset of different cell/space group.

        Only the model's *native* representation is written (PDB CRYST1, or mmCIF
        cell/symmetry): ChimeraX's writers cross-convert between the two via
        metadata (mmcif_write builds cell/symmetry from CRYST1; the PDB writer
        builds CRYST1 from cell/symmetry), so writing both risks confusing their
        "came-from-X" heuristics. No-op during session restore.'''
        if getattr(self, '_session_restore', False):
            return
        model = self.structure
        if model is None:
            return
        set_entry = getattr(model, 'set_metadata_entry', None)
        if set_entry is None:
            return
        if self.has_symmetry:
            a, b, c = self.cell.dim
            al, be, ga = self.cell.angles_deg
            hm = self.spacegroup.symbol_hm
        else:
            # Conventional "no crystallographic symmetry": unit cubic P1 box.
            a = b = c = 1.0
            al = be = ga = 90.0
            hm = 'P 1'
        md = model.metadata
        if 'cell' in md and 'CRYST1' not in md:
            # mmCIF-origin.
            if not self.has_symmetry:
                # Non-crystallographic: drop cell/symmetry so re-derivation falls
                # through to simple_p1_box (cleaner than fighting the _exptl gate).
                set_entry('cell', None)
                set_entry('symmetry', None)
            else:
                fmt = lambda v: '{:.4f}'.format(v)
                set_entry('cell', ['cell', 'length_a', 'length_b', 'length_c',
                                   'angle_alpha', 'angle_beta', 'angle_gamma'])
                set_entry('cell data', [fmt(a), fmt(b), fmt(c), fmt(al), fmt(be), fmt(ga)])
                set_entry('symmetry', ['symmetry', 'space_group_name_H-M'])
                set_entry('symmetry data', [hm])
                # symmetry_from_model_metadata_mmcif only reads the cell if
                # _exptl.method names a crystallographic experiment. Preserve an
                # existing crystallographic method (e.g. ELECTRON CRYSTALLOGRAPHY);
                # only set one when promoting a previously non-crystallographic model.
                CRYST = ('X-RAY DIFFRACTION', 'ELECTRON CRYSTALLOGRAPHY', 'NEUTRON DIFFRACTION')
                cur = self._metadata_field(md, 'exptl', 'method')
                if cur is None or cur.upper() not in CRYST:
                    rt = str(getattr(self.map_mgr, 'radiation_type', None)).lower()
                    method = 'ELECTRON CRYSTALLOGRAPHY' if rt == 'electron' else 'X-RAY DIFFRACTION'
                    set_entry('exptl', ['exptl', 'method'])
                    set_entry('exptl data', [method])
        else:
            # PDB-origin (or no prior symmetry metadata): one raw fixed-width
            # CRYST1 record, using ChimeraX's exact pdb_chars.cpp format -- which
            # is also exactly what symmetry_from_model_metadata_pdb parses by column.
            line = 'CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d' % (
                a, b, c, al, be, ga, hm, 1)
            set_entry('CRYST1', [line])

    def added_to_session(self, session):
        # Temporary hack due to issue #18427 in ChimeraX 1.10. TODO: remove for 1.11
        self.structure._cpp_notify_position(self.structure.scene_position)
        super().added_to_session(session)
        if self.structure is not None and not self._restored_skip_default_styling:
            self.set_default_atom_display(mode=self._hydrogen_mode)
        from .mousemodes import initialize_clipper_mouse_modes
        initialize_clipper_mouse_modes(session)

    def _structure_remove_cb(self, trigger_name, removed_models):
        if self.structure is None and not self._swapping_model:
            print('Deleting atomic symmetry model...')
            if self.atomic_symmetry_model is not None:
                asm = self.atomic_symmetry_model
                self.session.models.remove([asm])
                asm.delete()


    def swap_model(self, new_model, keep_old = False):
        '''
        Swap out the current atomic model for a new one.
        '''
        self._swapping_model = True
        old_model = self.structure
        if old_model is None:
            self._swapping_model=False
            return self.add_model(new_model)
        old_id = old_model.id
        old_model.triggers.remove_handler(self._structure_change_handler)
        self.session.models.remove([old_model])
        self._transplant_model(new_model)
        # self.add([new_model])
        self.session.models.assign_id(new_model,old_id)
        if keep_old:
            self.session.models.add([old_model])
        else:
            old_model.delete()
        new_model.triggers.add_handler('changes', self._structure_change_cb)
        def redraw_cb(*_):
            self.set_default_atom_display(mode=self._hydrogen_mode)
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER
        self.session.triggers.add_handler('frame drawn', redraw_cb)
        self._anisou_sanity_check(new_model.atoms)
        self.triggers.activate_trigger('model replaced', new_model)
        self._swapping_model = False
        return new_model

    def swap_model_from_file(self, filename, keep_old=False):
        '''
        Load a model from file, and replace the current coordinates with the
        result.
        '''
        from chimerax.open_command.cmd import provider_open
        new_model = provider_open(self.session, [filename], _add_models=False)[0]
        return self.swap_model(new_model, keep_old=keep_old)


    @property
    def position(self):
        return Model.position.fget(self)

    @position.setter
    def position(self, pos):
        Model.position.fset(self, pos)
        self.normalize_scene_positions()

    def normalize_scene_positions(self):
        apply_scene_positions(self.structure)

    def _transplant_model(self, model):
        if model.id is None:
            self.add([model])
            return
        mlist = model.all_models()
        mlist.sort(key=lambda m: len(m.id), reverse=True)
        from collections import defaultdict
        parent_to_models = defaultdict(list)
        for m in mlist:
            if m.parent in mlist:
                parent_to_models[m.parent].append(m)
        self.session.models.remove([model])
        self.add([model])
        def recursively_add_children(m):
            for child in parent_to_models[m]:
                m.add([child])
                recursively_add_children(child)
        recursively_add_children(model)


    @property
    def hydrogen_display_mode(self):
        return self._hydrogens

    @hydrogen_display_mode.setter
    def hydrogen_display_mode(self, mode):
        valid_names = ('polar', 'all', 'none')
        if mode not in valid_names:
            raise TypeError('Invalid mode! Mode must be one of ({})'.format(
                ','.join(valid_names)
            ))
        atoms = self.structure.atoms
        hydrogens = atoms[atoms.element_names == 'H']
        hydrogens.displays = False
        if mode == 'polar':
            from .util import exclude_nonpolar_hydrogens
            hydrogens[exclude_nonpolar_hydrogens(hydrogens)].displays = True
        elif mode == 'all':
            hydrogens.displays=True
        self._hydrogens = mode


    @property
    def structure(self):
        '''The atomic model managed by this symmetry manager.'''
        from chimerax.atomic import AtomicStructure
        for c in self.child_models():
            if type(c)==AtomicStructure:
                return c
        return None

    @property
    def map_mgr(self):
        ''' Master manager handling all maps associated with this model.'''
        from .maps import MapMgr
        for m in self.child_models():
            if isinstance(m, MapMgr):
                return m
        return MapMgr(self)

    @property
    def has_symmetry(self):
        return self._has_symmetry

    @has_symmetry.setter
    def has_symmetry(self, flag):
        self._has_symmetry = flag

    @property
    def last_covered_selection(self):
        return self._last_covered_selection

    @last_covered_selection.setter
    def last_covered_selection(self, sel):
        from chimerax.atomic import Atoms
        if not isinstance(sel, Atoms) and sel is not None:
            raise TypeError("Selection must be an Atoms array or None!")
        self._last_covered_selection = sel

    @property
    def has_experimental_maps(self):
        return self.xmapset is not None

    @property
    def stepper(self):
        '''
        Provides methods for "stepping" back and forth through the
        model according to secondary structure. For example, each call
        to stepper.step_forward() (with default arguments) will return
        an atom selection corresponding to the next pair of defined
        secondary structure elements plus their flanking loops.
        '''
        if self._stepper is None:
            from .structurestepper import StructureStepper
            self._stepper = StructureStepper(self.session, self.structure)
        return self._stepper

    @property
    def spotlight_mode(self):
        if not hasattr(self, "_spotlight_mode"):
            self._spotlight_mode=False
        return self._spotlight_mode

    @spotlight_mode.setter
    def spotlight_mode(self, flag):
        # Always fire the trigger even if we're already in spotlight mode, to ensure things behave correctly after user "vol unzone" commands
        if flag and auto_reset_camera:
            from chimerax.std_commands import cofr, camera
            cofr.cofr(self.session, 'centerOfView', show_pivot=True)
            camera.camera(self.session, 'ortho')
        self._spotlight_mode = flag
        self.triggers.activate_trigger('mode changed', None)
        if self.atomic_symmetry_model is not None:
            self.atomic_symmetry_model.live_scrolling = flag
        self.update(force=True)

    @property
    def spotlight_center(self):
        return self._last_box_center

    @spotlight_center.setter
    def spotlight_center(self, *_):
        raise NotImplementedError(
            'Centre of the spotlight should not be directly set. Change it by '
            'changing the centre of rotation of the main view: '
            'session.main_view.center_of_rotation = new_center'
        )

    @property
    def atomic_symmetry_model(self):
        for m in self.child_models():
            if isinstance(m, AtomicSymmetryModel):
                return m
        if self.structure is not None:
            return AtomicSymmetryModel(self, self.spotlight_radius, live=self.spotlight_mode, debug=self._debug)
        return None

    @property
    def atomic_sym_radius(self):
        if self.atomic_symmetry_model is not None:
            return self.atomic_symmetry_model.spotlight_radius
        return None

    @atomic_sym_radius.setter
    def atomic_sym_radius(self, radius):
        if self.atomic_symmetry_model is not None:
            self.atomic_symmetry_model.spotlight_radius = radius

    @property
    def spotlight_radius(self):
        return self.map_mgr.spotlight_radius

    @spotlight_radius.setter
    def spotlight_radius(self, radius):
        self.map_mgr.spotlight_radius = radius
        self.atomic_sym_radius = (radius + self.ATOMIC_SYM_EXTRA_RADIUS)

    @property
    def unit_cell(self):
        return self._unit_cell

    def discard_model_symmetry(self):
        if len(self.map_mgr.xmapsets):
            raise RuntimeError('Cannot discard model symmetry while a crystal '
                'dataset is loaded!')
        self.cell, self.spacegroup, self.grid, self.resolution, self._has_symmetry = \
            simple_p1_box(self.structure)
        self._unit_cell = Unit_Cell(self.structure.atoms, self.cell, self.spacegroup, self.grid)
        self._write_symmetry_to_model_metadata()


    def set_default_atom_display(self, mode='polar'):
        model = self.structure
        atoms = model.atoms
        atoms.draw_modes = atoms.STICK_STYLE
        atoms.displays=True
        self.hydrogen_display_mode = mode

    def delete(self):
        if self._update_handler is not None:
            self.session.triggers.remove_handler(self._update_handler)
            self._update_handler = None
        if self._structure_remove_handler is not None:
            self.session.triggers.remove_handler(self._structure_remove_handler)
            self._structure_remove_handler = None
        super().delete()

    def update(self, *_, force=False):
        v = self.session.main_view
        cofr = v.center_of_rotation
        # center = self.scene_position.inverse(is_orthonormal=True)*cofr
        update_needed = False
        if self._last_box_center is None:
            update_needed = True
        elif numpy.linalg.norm(cofr-self._last_box_center) > self.SPOTLIGHT_UPDATE_THRESHOLD:
            update_needed = True
        if update_needed or force:
            self._last_box_center = cofr
            self.triggers.activate_trigger('spotlight moved', cofr)


    def isolate_and_cover_selection(self, atoms, include_surrounding_residues=5,
        show_context=5, mask_radius=3, extra_padding=0, hide_surrounds=True,
        focus=True, include_hydrogens=False):
        '''
        Expand the map(s) (if present) to cover a given atomic selection,  then
        mask them to within a given distance of said atoms to reduce visual
        clutter. Adjust the atomic visualisation to show only the selected
        atoms, plus an optional surrounding buffer zone.
        Args:

          * atoms (ChimeraX `Atoms` object):

            The main selection we're interested in. The existing selection will
            be expanded to include the whole residue for every selected atom.

          * include_surrounding_residues (float):

            Any residue with an atom coming within this radius of the primary
            selection (including symmetry atoms) will be added to the selection
            covered by the map. To cover only the primary selection, set this
            value to zero.

          * show_context (float):

            Any residue within an atom coming within this radius of the previous
            two selections will be displayed, but will not be considered for the
            map masking calculation.

          * mask_radius (float):

            Components of the map more than this distance from any eligible atom
            will be hidden.

          * extra_padding (float):

            Optionally, further pad the volume by this distance. The extra
            volume will be hidden, but available for calculations.

          * hide_surrounds (bool):

            If true, all residues outside the selection region will be hidden

          * focus (bool):
            If true, the camera will be moved to focus on the selection (only
            the atoms in the master model will be considered)

          * include_hydrogens (bool):

            Include hydrogens for the purposes of masking (this should almost
            never be necessary)
        '''
        self.spotlight_mode = False
        original_atoms = self.last_covered_selection = atoms
        atoms = original_atoms.unique_residues.atoms
        atoms = atoms[atoms.element_names != 'H']
        asm = self.atomic_symmetry_model
        main_set, context_set = asm.show_selection_in_context(atoms,
            include_surrounding_residues=include_surrounding_residues,
            show_context=show_context
        )
        from chimerax.geometry import Places
        atoms = main_set[0]
        transforms = Places(place_array=main_set[1])
        indices = main_set[2]
        if not include_hydrogens:
            mask = (atoms.element_names != 'H')
            atoms = atoms[mask]
            indices = indices[mask]
        # self._map_cover_data = main_set
        # coords = self._update_cover_coords(main_set)
        self.map_mgr.cover_atoms(atoms, transforms=transforms, transform_indices=indices,
            mask_radius=mask_radius,
            extra_padding=extra_padding)
        if not hide_surrounds:
            self.structure.atoms.hides &=~HIDE_ISOLDE
        if focus:
            focus_on_selection(self.session, self.session.main_view, atoms)

    _cube_pairs = numpy.array([[0,1], [0,2], [0,4], [1,3], [1,5], [2,3], [2,6], [3,7], [4,5], [4,6], [5,7], [6,7]], int)

    def _update_cover_coords(self, cover_set):
        from chimerax.geometry import Places
        transforms = Places(place_array=cover_set[1])
        indices = cover_set[2]
        coords = cover_set[0].coords
        import numpy
        for i in numpy.unique(indices):
            mask = indices==i
            coords[mask] = transforms[i]*coords[mask]
        return coords

    def draw_unit_cell_box(self, offset=None, cylinder_radius = 0.05):
        '''
        Draw a parallellepiped around one unit cell.
        '''
        m, d = _get_special_positions_model(self)
        model = self.structure
        uc = self.unit_cell

        if offset is None:
            offset = numpy.zeros(3)

        from chimerax.geometry import Place, Places
        positions = []
        colors = []
        rgba_edge = numpy.array([0,255,255,128],numpy.uint8)
        corners_frac = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],numpy.double) + offset\
                        + uc.origin.coord_frac(self.grid).uvw

        corners = numpy.array([clipper_python.Coord_frac(c).coord_orth(self.cell).xyz for c in corners_frac])
        from chimerax.surface.shapes import cylinder_geometry
        d.set_geometry(*cylinder_geometry())
        d.set_color(rgba_edge)
        from chimerax.geometry import cylinder_rotations, Places
        xyz0 = numpy.array([corners[p[0]] for p in self._cube_pairs])
        xyz1 = numpy.array([corners[p[1]] for p in self._cube_pairs])
        radii = numpy.ones(12)*cylinder_radius
        p = numpy.empty((12, 4, 4), numpy.float32)
        cylinder_rotations(xyz0, xyz1, radii, p)
        p[:,3,:3] = 0.5*(xyz0+xyz1)
        pl = Places(opengl_array = p)
        d.set_positions(pl)
        m.display = True




    def draw_unit_cell_and_special_positions(self, offset = None):
        '''
        Quick-and-dirty drawing mapping out the special positions
        (positions which map back to themselves by at least one
        non-unity symop) within one unit cell. A sphere will be drawn
        at each grid-point with non-unit multiplicity, and colour-coded
        according to multiplicity:

            * 2-fold: white
            * 3-fold: cyan
            * 4-fold: yellow
            * 6-fold: magenta

        Ultimately it would be nice to replace this with something more
        elegant, that masks and scrolls continuously along with the model/
        map visualisation.

        Args:

            * offset (1x3 numpy array, default = None):
                
                Optional (u,v,w) offset (in fractions of a unit cell axis)
        '''
        m, d = _get_special_positions_model(self)
        model = self.structure

        if offset is None:
            offset = numpy.array([0,0,0],int)

        xmap = self.xmapset[0].xmap
        uc = self.unit_cell
        spc = numpy.array(xmap.special_positions_unit_cell_xyz(uc, offset))
        from chimerax.surface.shapes import sphere_geometry2
        sphere = numpy.array(sphere_geometry2(80))
        sphere[0]*=0.25
        d.set_geometry(*sphere)
        #d.vertices, d.normals, d.triangles = sphere

        from chimerax.geometry import Place, Places
        positions = []
        colors = []
        rgba_corner = numpy.array([255,0,255,128],numpy.uint8)
        corners_frac = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],numpy.double) + offset\
                        + self.unit_cell.min.coord_frac(self.grid).uvw

        corners = []
        for c in corners_frac:
            co = clipper_python.Coord_frac(c).coord_orth(self.cell)
            positions.append(Place(axes=numpy.identity(3)*4, origin=co.xyz))
            colors.append(rgba_corner)

        if len(spc):

            coords = spc[:,0:3]
            multiplicity = spc[:,3].astype(int)
            scale_2fold = numpy.identity(3)
            scale_3fold = numpy.identity(3)* 1.5
            scale_4fold = numpy.identity(3)* 2
            scale_6fold = numpy.identity(3)* 3
            rgba_2fold = numpy.array([255,255,255,255],numpy.uint8)
            rgba_3fold = numpy.array([0,255,255,255],numpy.uint8)
            rgba_4fold = numpy.array([255,255,0,255],numpy.uint8)
            rgba_6fold = numpy.array([255,0,0,255],numpy.uint8)

            for coord, mult in zip(coords, multiplicity):
                if mult == 2:
                    positions.append(Place(axes=scale_2fold, origin=coord))
                    colors.append(rgba_2fold)
                elif mult == 3:
                    positions.append(Place(axes=scale_3fold, origin=coord))
                    colors.append(rgba_3fold)
                elif mult == 4:
                    positions.append(Place(axes=scale_4fold, origin=coord))
                    colors.append(rgba_4fold)
                elif mult == 6:
                    positions.append(Place(axes=scale_6fold, origin=coord))
                    colors.append(rgba_6fold)
            for c in corners:
                positions.append(Place(axes=scale_6fold, origin=c))
                colors.append(rgba_corner)

        d.set_positions(Places(positions))
        d.set_colors(numpy.array(colors, numpy.uint8))
        # m.add_drawing(d)
        # model.parent.add([m])
        m.display = True

    def take_snapshot(self, session, flags):
        from chimerax.core.models import Model
        data = {
            'model state': Model.take_snapshot(self, session, flags),
            'hall symbol': self.spacegroup.symbol_hall,
            'resolution': self.resolution.limit,
            'cell dim': self.cell.dim,
            'cell angles': self.cell.angles_deg,
            'has symmetry': self.has_symmetry
        }
        from . import CLIPPER_STATE_VERSION
        data['version']=CLIPPER_STATE_VERSION
        return data

    @staticmethod
    def restore_snapshot(session, data):
        sh = SymmetryManager(session)
        # Preserve the restored model's saved display: block default styling until
        # the deferred 'frame drawn' callbacks have run (cleared in _end_restore_session_cb).
        sh._restored_skip_default_styling = True
        Model.set_state_from_snapshot(sh, session, data['model state'])
        session.triggers.add_handler('end restore session', sh._end_restore_session_cb)
        sh._session_restore_data = data
        return sh

    def _end_restore_session_cb(self, *_):
        self._session_restore=True
        if self.structure is not None:
            data = self._session_restore_data
            if data['version'] < 2:
                data = None
            self.add_model(self.structure, session_restore_data=data)
        delattr(self, '_session_restore_data')
        self._session_restore=False
        # NB: _restored_skip_default_styling is NOT cleared here on a frame count
        # (that released it before the deferred styling callback fired). It is
        # cleared inside AtomicSymmetryModel._set_default_cartoon_cb, which is the
        # last deferred one-time styling event — so the flag lives exactly as long
        # as needed regardless of how many frames the callback takes to fire.
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


def _get_special_positions_model(m):
    for cm in m.child_models():
        if cm.name == "Special positions":
            return cm, cm.child_drawings()[0]
    from chimerax.core.models import Model, Drawing
    d = Drawing("Special positions")
    sm = Model("Special positions", m.session)
    sm.add_drawing(d)
    m.add([sm])
    return sm, d


class AtomicSymmetryModel(Model):
    '''
    Finds and draws local symmetry atoms for an atomic structure
    '''
    from .clipper_python.ext import element_radius_map
    _element_radii = element_radius_map()
    SESSION_SAVE=False
    def __init__(self, parent, radius = 15,
        dim_colors_to = 0.4, backbone_mode = 'CA trace', live = True, debug=False):
        self._debug = debug
        self._live_scrolling = False
        session = self.session = parent.session
        self._color_dim_factor = dim_colors_to

        self._sym_search_frequency = 2
        self._current_focal_set = None

        self.manager = parent

        super().__init__('Atomic symmetry', session)
        parent.add([self])
        atomic_structure = self.structure
        self._last_hides = atomic_structure.atoms.hides

        self._box_changed_handler = None
        self._center = session.main_view.center_of_rotation
        self._spotlight_radius = radius
        self._box_dim = numpy.array([radius*2, radius*2, radius*2], numpy.double)
        ad = self._atoms_drawing = SymAtomsDrawing('Symmetry atoms')
        self.add_drawing(ad)
        #from chimerax.core.atomic.structure import PickedBonds
        bd = self._bonds_drawing = SymBondsDrawing('Symmetry bonds', PickedSymBond, structure.PickedBonds)
        self.add_drawing(bd)
        rd = self._ribbon_drawing = SymRibbonDrawing('Symmetry ribbons',
            self, dim_colors_to)
        self.add_drawing(rd)
        self._current_atoms = atomic_structure.atoms
        self._model_changes_handler = self.manager.triggers.add_handler(
                                        'atoms changed', self._model_changed_cb)
        self._model_swap_handler = self.manager.triggers.add_handler(
                                        'model replaced', self._model_swap_cb)
        # When a fragment is shifted to a different ASU ("Move symmetry fragment
        # here"), the Unit_Cell's cached reference model no longer represents the
        # atoms. Mark it stale; the refresh is done lazily on the next box search.
        self._reference_dirty = False
        self._sym_shift_handler = self.manager.triggers.add_handler(
                                        'sym shifted atoms', self._sym_shifted_cb)
        self.live_scrolling = live
        self._save_session_handler = session.triggers.add_handler('begin save session',
            self._start_save_session_cb)
        # On session restore, keep the restored atom radii; only assign default
        # per-element radii for a genuine fresh association.
        if not getattr(self.manager, '_restored_skip_default_styling', False):
            self._assign_atom_radii(self.structure.atoms)

    def added_to_session(self, session):
        super().added_to_session(session)
        # Have to delay setting the drawing settings until after the frame is
        # drawn, otherwise they get overridden by ChimeraX
        session.triggers.add_handler('frame drawn', self._set_default_cartoon_cb)


    def _set_default_cartoon_cb(self, *_):
        from chimerax.core.triggerset import DEREGISTER
        if not hasattr(self, 'session'):
            # Most likely the model has been closed prior to drawing. Just deregister
            return DEREGISTER
        if getattr(self.manager, '_restored_skip_default_styling', False):
            # Session restore carries the saved cartoon/colour/display state;
            # don't overwrite it. This is the last deferred one-time styling event
            # of a restore, so releasing the flag here (rather than on a fragile
            # frame count) is exactly when the restore-styling window closes -- and
            # lets later genuine restyles (swap_model, user-added atoms) proceed.
            self.manager._restored_skip_default_styling = False
            return DEREGISTER
        from .util import set_to_default_cartoon
        set_to_default_cartoon(self.session, model = self.structure)
        return DEREGISTER


    @property
    def unit_cell(self):
        return self.manager.unit_cell

    @property
    def cell(self):
        return self.manager.cell

    @property
    def spacegroup(self):
        return self.manager.spacegroup

    @property
    def grid(self):
        return self.manager.grid

    @property
    def structure(self):
        return self.manager.structure

    @property
    def _box_corner_drawing(self):
        if not hasattr(self, "_bcd") or self._bcd is None:
            from chimerax.core.models import Drawing
            from chimerax.surface.shapes import sphere_geometry2
            d = self._bcd = Drawing("symmetry box corners")
            d.set_geometry(*sphere_geometry2(200))
            self.add_drawing(d)
        return self._bcd

    @property
    def _search_positions_drawing(self):
        if not hasattr(self, "_spd") or self._spd is None:
            from chimerax.core.models import Drawing
            from chimerax.surface.shapes import sphere_geometry2
            d = self._spd = Drawing("symmetry search positions")
            v, n, t = sphere_geometry2(200)
            from chimerax.geometry import scale
            v = scale(0.5).transform_points(v)
            d.set_geometry(v,n,t)
            self.add_drawing(d)
        return self._spd

    def show_selection_in_context(self, atoms, include_surrounding_residues=5,
            show_context=5):
        self._last_context_params = (atoms, include_surrounding_residues, show_context)
        main_set = self.sym_select_within(atoms, include_surrounding_residues)
        main_coords = symmetry_coords(*main_set[0:3])
        context_set = self._current_focal_set = self.sym_select_within(
            main_set[0], show_context, coords=main_coords
        )
        self.set_sym_display(context_set[3],
            *atom_and_bond_sym_transforms_from_sym_atoms(*context_set[0:3])
        )
        return main_set, context_set

    def sym_select_within(self, atoms, cutoff, coords=None, whole_residues = True):
        '''
        Given a set of atoms, return a (atoms, symmetry matrices, sym_indices,
        symops) tuple giving all atoms and their symmetry operators within the
        given cutoff distance from any atom in the primary set.
        Args:
            atoms:
                The core set of atoms to be surrounded by the new selection.
            cutoff:
                Cutoff distance in Angstroms.
            coords (default None):
                Optionally, you can provide the coordinates for all atoms
                (useful for expanding a selection that already includes
                symmetry atoms).
            whole_residues (default true):
                Whether to expand the selections to whole_residues.
        '''
        from .sym_realize import sym_select_within as _sym_select_within
        return _sym_select_within(self.structure, self.cell, self.grid,
            self.unit_cell, atoms, cutoff, coords=coords,
            whole_residues=whole_residues,
            sym_search_frequency=self._sym_search_frequency)

    def _places_from_matrices(self, matrices, indices, include_identity):
        '''
        Turn a set of orthogonal-space 3x4 symmetry matrices and a collection of
        indices into them into a list of ChimeraX Place objects. The identity
        (ASU) operator lives at index 0 (see Unit_Cell.all_symops_in_box) and is
        placed first when requested; all other operators follow in ascending
        index order so the result is deterministic.
        '''
        from .sym_realize import places_from_matrices
        return places_from_matrices(matrices, indices, include_identity)

    def currently_displayed_sym_transforms(self, include_identity=True):
        '''
        Return the symmetry operators (as a list of ChimeraX Place objects, in the
        master model's own coordinate frame) for every symmetry copy currently
        drawn in the spotlight - i.e. every copy whose ribbon is being displayed.
        The identity/ASU transform is included first unless include_identity is
        False. Returns an empty list if nothing is currently displayed.
        '''
        tfs = getattr(self, '_current_tfs', None)
        ribbon_syms = getattr(self, '_current_ribbon_syms', None)
        if tfs is None or ribbon_syms is None or not len(tfs):
            return []
        return self._places_from_matrices(tfs, ribbon_syms, include_identity)

    def drawn_ghosts_by_operator(self, include_identity=False):
        '''
        Return the symmetry ghosts currently drawn, grouped by operator, as a list
        of ``(Place, Atoms)`` tuples: for each operator with at least one ghost on
        screen, the operator (as a :class:`Place` in the master model's coordinate
        frame) and the master :class:`Atoms` whose ghost copies are drawn under it.
        The identity/ASU operator (the real atoms) is excluded unless
        ``include_identity`` is True. Returns an empty list if nothing is drawn.

        This is the public accessor for "what symmetry is on screen right now",
        intended for dependent bundles that draw their own decorations on the
        ghosts. It reflects exactly what :meth:`update_graphics` last drew; pair it
        with the manager's ``'symmetry display changed'`` trigger to stay in sync.
        '''
        sym_atoms = getattr(self, '_current_sym_atoms', None)
        atom_syms = getattr(self, '_current_atom_syms', None)
        tfs = getattr(self, '_current_tfs', None)
        if sym_atoms is None or atom_syms is None or tfs is None or not len(sym_atoms):
            return []
        from chimerax.geometry import Place
        out = []
        for idx in numpy.unique(atom_syms):
            idx = int(idx)
            if idx == 0 and not include_identity:
                continue
            out.append((Place(matrix=tfs[idx]), sym_atoms[atom_syms == idx]))
        return out

    def sym_transforms_near(self, atoms, cutoff, include_identity=True):
        '''
        Return the symmetry operators (as a list of ChimeraX Place objects, in the
        master model's own coordinate frame) for every symmetry copy with at least
        one atom approaching within cutoff (Angstroms) of the given atoms. The
        identity/ASU transform is included first unless include_identity is False.
        '''
        found_atoms, symmats, sym_indices, symops = self.sym_select_within(atoms, cutoff)
        return self._places_from_matrices(symmats, sym_indices, include_identity)

    def delete(self):
        bh = self._box_changed_handler
        mh = self._model_changes_handler
        sh = self._model_swap_handler
        ssh = self._sym_shift_handler
        for h in (bh, mh, sh, ssh):
            if h is not None:
                self.manager.triggers.remove_handler(h)
        if self.structure is not None:
            self.unhide_all_atoms()
        self.session.triggers.remove_handler(self._save_session_handler)
        super().delete()

    @property
    def backbone_mode(self):
        return _backbone_mode_descr[self._backbone_mode]

    @backbone_mode.setter
    def backbone_mode(self, mode):
        old_mode = self.backbone_mode
        if mode == "ribbon":
            self._backbone_mode = BACKBONE_MODE_RIBBON
        elif mode == "CA trace":
            self._backbone_mode = BACKBONE_MODE_CA_TRACE
        else:
            raise TypeError('Unrecognised mode! Should be one of "ribbon" or "CA trace"')
        if old_mode != mode:
            # 'backbone mode changed' is registered on the manager, not on this
            # model's own (Model-default) TriggerSet, so fire it there.
            self.manager.triggers.activate_trigger('backbone mode changed', mode)

    def unhide_all_atoms(self):
        self.structure.atoms.hides &= ~HIDE_ISOLDE

    @property
    def dim_factor(self):
        return self._color_dim_factor

    @dim_factor.setter
    def dim_factor(self, factor):
        self._color_dim_factor = factor
        self._ribbon_drawing.dim_factor = factor

    @property
    def live_scrolling(self):
        return self._live_scrolling

    @live_scrolling.setter
    def live_scrolling(self, flag):
        bh = self._box_changed_handler
        if flag and not self._live_scrolling:
            if bh is None:
                self._box_changed_handler = self.parent.triggers.add_handler('spotlight moved',
                    self._box_moved_cb)
                self._update_box()
                self.update_graphics()
        elif not flag and self._live_scrolling:
            from chimerax.atomic import concatenate
            res = whole_residue_sym_sphere(self.structure.residues, self._current_tfs, self._center, self._spotlight_radius, visible_only=False)
            ma = res[0]
            sa = res[1]
            fa = concatenate((ma, sa))
            fs = numpy.concatenate((numpy.zeros(len(ma), numpy.uint8), res[3]))
            self._current_focal_set = (fa, self._current_tfs, fs, self._current_symops)
            self.unhide_all_atoms()
            if bh is not None:
                self.parent.triggers.remove_handler(bh)
                self._box_changed_handler = None
        self._live_scrolling = flag

    @property
    def spotlight_radius(self):
        return self._spotlight_radius

    @spotlight_radius.setter
    def spotlight_radius(self, radius):
        self._spotlight_radius = radius
        if self.visible:
            self.update_graphics()

    @property
    def display(self):
        return super().display

    @display.setter
    def display(self, flag):
        if not flag:
            self.unhide_all_atoms()
        super().set_display(flag)
        if flag:
            if self.live_scrolling:
                self._center = self.session.main_view.center_of_rotation
                self._update_box()
            self.update_graphics()

    @property
    def _level_of_detail(self):
        from chimerax.atomic import structure_graphics_updater
        gu = structure_graphics_updater(self.session)
        return gu.level_of_detail


    def _box_moved_cb(self, trigger_name, center):
        if not self.visible:
            return
        self._center = center
        self._update_box()
        self.update_graphics()

    def _sym_shifted_cb(self, trigger_name, data):
        # A bonded group was moved to a different ASU. Flag the Unit_Cell reference
        # for refresh (done lazily in _update_box) and re-run the search now so the
        # symmetry display is correct immediately.
        self._reference_dirty = True
        if self.visible and self.live_scrolling:
            self._update_box()
            self.update_graphics()

    def _update_box(self):
        from .crystal import find_box_params
        # Lazily refresh the Unit_Cell reference if the model changed in a way that
        # invalidates it (set via _sym_shifted_cb). Only fires when actually dirty,
        # so ordinary per-frame box updates during a simulation stay free.
        if getattr(self, '_reference_dirty', False):
            self.unit_cell.update_reference(self.structure.atoms)
            self._reference_dirty = False
        center = self._center
        radius = self._spotlight_radius
        box_corner_grid, box_corner_xyz, grid_dim = find_box_params(center, self.cell, self.grid, radius)
        dim = self._box_dim
        dim[:] = radius*2
        symops = self._current_symops = self.unit_cell.all_symops_in_box(box_corner_xyz, grid_dim, True, self._sym_search_frequency)
        # Identity symop will always be the first in the list
        tfs = self._current_tfs = symops.all_matrices_orth(self.cell, '3x4')
        atoms = self.structure.atoms
        atoms.hides |=HIDE_ISOLDE
        self._current_master_atoms, self._current_sym_atoms, self._current_sym_atom_coords, \
        self._current_atom_syms, self._current_bonds, self._current_bond_tfs, \
        self._current_bond_syms = \
                whole_residue_sym_sphere(self.structure.residues, tfs, center, radius)
                #sym_transforms_in_sphere(atoms, tfs[first_symop:], center, radius)
        self._current_master_atoms.hides &= ~HIDE_ISOLDE
        cs = self._current_sym_atoms
        if len(cs):
            crs = self._current_ribbon_syms = numpy.unique(self._current_atom_syms)
            self._current_ribbon_syms = crs[crs!=0]
        else:
            self._current_ribbon_syms = sym_ribbons_in_sphere(self._ribbon_drawing._tether_coords, tfs, center, radius)
        if(self._debug):
            self._draw_box_corners(box_corner_grid, grid_dim)
            self._draw_search_positions(box_corner_grid, grid_dim, self._sym_search_frequency)

    def _draw_box_corners(self, box_corner_grid, grid_dim):
        from . import clipper
        corner_mask = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]])
        colors = numpy.ones((8,4), numpy.uint8)*50
        colors[:,3]=255
        colors[:,:3]+=(corner_mask*200).astype(numpy.uint8)
        from chimerax.geometry import Place, Places
        places=[]
        for c in corner_mask:
            places.append(Place(origin=(box_corner_grid+clipper.Coord_grid(grid_dim*c)).coord_frac(self.grid).coord_orth(self.cell).xyz))
        self._box_corner_drawing.positions = Places(places)
        self._box_corner_drawing.colors = colors

    def _draw_search_positions(self, box_corner_grid, grid_dim, search_frequency):
        ref_size = self.unit_cell.ref_box.nuvw
        from .clipper_python import Coord_grid
        box_max = box_corner_grid + Coord_grid(grid_dim)
        step_size = max(ref_size.min()//search_frequency, 1)
        from chimerax.geometry import Place, Places
        places = []
        u = box_corner_grid.u
        u_done = False
        while not u_done:
            v = box_corner_grid.v
            v_done = False
            if u == box_max.u:
                u_done = True
            while not v_done:
                w = box_corner_grid.w
                w_done = False
                if v == box_max.v:
                    v_done = True
                while not w_done:
                    if w == box_max.w:
                        w_done = True
                    p = Coord_grid(u, v, w).coord_frac(self.grid).coord_orth(self.cell)
                    places.append(Place(origin=p.xyz))
                    w = min(w+step_size, box_max.w)
                v = min(v+step_size, box_max.v)
            u = min(u+step_size, box_max.u)
        self._search_positions_drawing.positions = Places(places)

    def set_sym_display(self, symops, primary_atoms, sym_atoms, sym_coords, atom_sym_indices,
        sym_bonds, sym_bond_tfs, bond_sym_indices, update_hides=True):
        '''
        Manually define the set of symmetry atoms to be displayed.
        Args:
            primary_atoms:
                A ChimeraX Atoms object defining the atoms in the master model
                to be displayed.
            symops:
                A Clipper Symops object defining the symmetry transforms
            sym_atoms:
                A ChimeraX Atoms object listing the symmetry atoms to be displayed
            sym_coords:
                The symmetry coordinates corresponding to sym_atoms
            atom_sym_indices:
                A Numpy integer array of length equal to sym_atoms, giving the
                index of the corresponding symmetry operator in symops
            sym_bonds:
                A ChimeraX Bonds object listing the symmetry bonds to be
                displayed.
            sym_bond_tfs:
                A Places object providing the halfbond cylinder transforms to
                draw the sym_bonds in their correct positions
            bond_sym_indices:
                A Numpy integer array of length equal to sym_bonds, giving the
                index of the corresponding symmetry operator in symops
        '''
        self.live_scrolling = False
        if update_hides:
            all_atoms = self.structure.atoms
            all_atoms.hides |= HIDE_ISOLDE
            primary_atoms.hides &= ~HIDE_ISOLDE
        self._current_symops = symops
        self._current_tfs = symops.all_matrices_orth(self.cell, '3x4')
        self._current_master_atoms = primary_atoms
        self._current_sym_atoms = sym_atoms
        self._current_sym_atom_coords = sym_coords
        csym = self._current_atom_syms = atom_sym_indices
        self._current_bonds = sym_bonds
        self._current_bond_tfs = sym_bond_tfs
        self._current_bond_syms = bond_sym_indices
        self._current_ribbon_syms = numpy.unique(csym[csym!=0])
        self.update_graphics()

    def _model_swap_cb(self, trigger_name, new_model):
        print('Updating cartoon style...')
        self.spotlight_mode = True
        from .util import set_to_default_cartoon
        self._assign_atom_radii(new_model.atoms)
        set_to_default_cartoon(self.session, model = self.structure)

        # self.session.triggers.add_handler('frame drawn', self._set_default_cartoon_cb)

    def _assign_atom_radii(self, atoms):
        '''
        For large models, the recalculation of idatm_types triggered by addition
        of new atoms when radii are left automatic becomes prohibitive for
        editing the model. So we just set them to fixed values.
        '''
        erm = self._element_radii
        unique_elements = set(atoms.element_names)
        for e in unique_elements:
            radius = erm.get(e, 1)
            atoms[atoms.element_names==e].radii = radius

    def _set_new_atom_style(self, atoms):
        from chimerax.atomic import selected_atoms, Atom
        from chimerax.core.commands import run
        session=self.session
        atoms.draw_modes = Atom.STICK_STYLE
        residues = atoms.unique_residues
        residues.ribbon_displays=True
        residues.ribbon_hide_backbones=False
        current_sel = selected_atoms(session)
        from chimerax.std_commands.color import color
        from chimerax.core.objects import Objects
        objects = Objects(atoms=atoms)
        color(session, objects, color='bychain', halfbond=True)
        color(session, objects, color='byhetero', halfbond=True)
        from .util import nonpolar_hydrogens
        atoms[nonpolar_hydrogens(atoms)].displays=False
        current_sel.selected = True


    def _model_changed_cb(self, trigger_name, changes):
        changes = changes[1]
        num_atoms_changed = False
        update_needed = False
        ribbon_update_needed = False
        created_atoms = changes.created_atoms()
        # Skip default radii/style while a session restore is settling: any
        # atoms reported here then are the restored structure's own atoms, whose
        # saved radii/colours/draw-modes must be preserved. Genuine user-added
        # atoms (after the restore flag clears) are styled normally.
        if len(created_atoms) and not getattr(self.manager, '_restored_skip_default_styling', False):
            self._assign_atom_radii(created_atoms)
            self._set_new_atom_style(created_atoms)
        if len(created_atoms) or changes.num_deleted_atoms()>0:
            num_atoms_changed = True
            self._last_hides = self.structure.atoms.hides&~HIDE_ISOLDE
            update_needed = True
            ribbon_update_needed = True
        if not self.visible:
            return
        reasons = changes.atom_reasons()
        if 'coord changed' in reasons:
            ribbon_update_needed = True
            update_needed = True
        if 'display changed' in reasons:
            ribbon_update_needed = True
            update_needed = True
        if  'hide changed' in reasons:
            hides = self.structure.atoms.hides
            # Prevent repeated callback with every display update in spotlight mode
            current_hides = hides&~HIDE_ISOLDE
            if not update_needed and numpy.any(current_hides != self._last_hides):
                update_needed = True
            self._last_hides = current_hides
        if 'color changed' in reasons:
            ribbon_update_needed = True
            update_needed = True
        if 'ribbon_display changed' in changes.residue_reasons():
            ribbon_update_needed = True
        if (ribbon_update_needed):
            self._ribbon_drawing.delayed_rebuild(self.session)
        if (update_needed):
            if self.live_scrolling:
                self._update_box()
                self.update_graphics()
            else:
                if num_atoms_changed:
                    self.show_selection_in_context(*self._last_context_params)

                self._update_sym_coords()


    def _update_sym_coords(self):
        focal_set = self._current_focal_set
        # try:
        res = atom_and_bond_sym_transforms_from_sym_atoms(*focal_set[0:3])
        self.set_sym_display(focal_set[3], *res, update_hides=False)
        # except:
        #     from chimerax.atomic import Atoms, Bonds
        #     self._current_sym_atoms = Atoms()
        #     self._current_bonds = Bonds()
        #     self.update_graphics()


    def update_graphics(self):
        lod = self._level_of_detail
        self._update_atom_graphics(lod)
        self._update_bond_graphics(lod)
        self._update_ribbon_graphics()
        # Every path that changes the displayed ghost set converges here, and all
        # _current_* state is fresh at this point. Notify dependent bundles (e.g.
        # ISOLDE) so they can redraw markup on the ghosts. Fired on the durable
        # manager so subscriptions survive AtomicSymmetryModel recreation.
        self.manager.triggers.activate_trigger('symmetry display changed', self)

    def _update_atom_graphics(self, lod):
        ad = self._atoms_drawing
        ca = self._current_sym_atoms
        if not len(ca):
            ad.display=False
            return
        syms = self._current_atom_syms
        coords = self._current_sym_atom_coords
        # if not self.live_scrolling:
        #     mask = _atoms_only_hidden_by_clipper(ca)
        #     ca = ca[mask]
        #     syms = syms[mask]
        #     coords = coords[mask]
        ad.visible_atoms = ca
        lod.set_atom_sphere_geometry(ad)

        na = len(ca)
        if na > 0:
            ad.display = True
            xyzr = numpy.empty((na, 4), numpy.float32)
            xyzr[:,:3] = coords
            xyzr[:,3] = self.structure._atom_display_radii(ca)
            from chimerax.geometry import Places
            ad.positions = Places(shift_and_scale = xyzr)
            colors = ca.colors.astype(numpy.float32)
            colors[:,:3] *= self._color_dim_factor
            ad.colors = colors.astype(numpy.uint8)
        else:
            ad.display = False


    def _update_bond_graphics(self, lod):
            bd = self._bonds_drawing
            bonds = self._current_bonds
            if not len(bonds):
                bd.display = False
                return
            bsym = self._current_bond_syms
            b_tfs = self._current_bond_tfs
            # if not self.live_scrolling:
            #     mask = _bonds_only_hidden_by_clipper(bonds)
            #     bonds = bonds[mask]
            #     mask = numpy.concatenate((mask, mask))
            #     bsym = bsym[mask]
            #     from chimerax.geometry import Places
            #     b_tfs = Places(place_array=b_tfs.array()[mask])
            lod.set_bond_cylinder_geometry(bd)
            bd.visible_bonds = bonds
            nb = len(bonds)
            if nb > 0:
                bd.display = True
                bd.positions = b_tfs
                colors = bonds.half_colors.astype(numpy.float32)
                colors[:,:3] *= self._color_dim_factor
                bd.colors = colors.astype(numpy.uint8)
            else:
                bd.display = False

    def _start_save_session_cb(self, *_):
        self._session_save_hides = self.structure.atoms.hides
        self.unhide_all_atoms()
        self.session.triggers.add_handler('end save session', self._end_save_session_cb)

    def _end_save_session_cb(self, *_):
        self.structure.atoms.hides = self._session_save_hides
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _update_ribbon_graphics(self):
        rd = self._ribbon_drawing
        prd = self.structure._ribbons_drawing
        position_indices = self._current_ribbon_syms
        tfs = self._current_tfs[position_indices]
        from chimerax.geometry import Places
        rd.positions = Places(place_array=tfs)
        rd.display = prd.display
        if not rd.display:
            return
        rd.update()

#from chimerax.core.atomic.structure import AtomsDrawing
class SymAtomsDrawing(structure.AtomsDrawing):
    def first_intercept(self, mxyz1, mxyz2, exclude=None):
        if not self.display or self.visible_atoms is None or (exclude and exclude(self)):
            return None
        xyzr = self.positions.shift_and_scale_array()
        coords, radii = xyzr[:,:3], xyzr[:,3]

        from chimerax.geometry import closest_sphere_intercept
        fmin, anum = closest_sphere_intercept(coords, radii, mxyz1, mxyz2)
        if fmin is None:
            return None
        atom = self.visible_atoms[anum]
        atom_syms = self.parent._current_atom_syms
        sym = self.parent._current_symops[int(atom_syms[anum])]
        return PickedSymAtom(atom, fmin, sym)

    def planes_pick(self, planes, exclude=None):
        return []

    def bounds(self, positions=True):
        if not positions:
            return self.geometry_bounds()
        cpb = self._cached_position_bounds
        if cpb is not None:
            return cpb
        xyzr = self.positions.shift_and_scale_array()
        if xyzr is None:
            return self.geometry_bounds()
        coords, radii = xyzr[:, :3], xyzr[:,3]
        from chimerax.geometry import sphere_bounds
        b = sphere_bounds(coords, radii)
        self._cached_position_bounds = b

    def update_selection(self):
        pass

from chimerax.graphics import Pick

class PickedSymAtom(Pick):
    def __init__(self, atom, distance, sym):
        super().__init__(distance)
        self.atom = atom
        self.sym = sym

    def description(self):
        return '({}) {}'.format(self.sym.format_as_symop(), self.atom)

    def select(self, mode = 'add'):
        pass

#from chimerax.core.atomic.structure import BondsDrawing
class SymBondsDrawing(structure.BondsDrawing):
    def first_intercept(self, mxyz1, mxyz2, exclude=None):
        if not self.display or (exclude and exclude(self)):
            return None
        asm = self.parent
        bonds = getattr(asm, '_current_bonds', None)
        if bonds is None or not len(bonds):
            return None
        tfs = asm._current_tfs
        symops = asm._current_symops
        # _current_bond_syms has one entry per half-bond (2 x len(bonds)), laid out
        # blocked as [half0 for each bond, half1 for each bond]; both halves of a
        # bond share an operator, so the per-bond operator index is the first half.
        nb = len(bonds)
        bond_syms = numpy.asarray(asm._current_bond_syms)[:nb]
        from chimerax.geometry import Place, closest_cylinder_intercept
        best_bond = best_f = best_k = None
        # A ghost bond under operator k is the master bond transformed by the
        # operator, so intercepting the ray with the ghost is the same as
        # intercepting the *inverse-transformed* ray with the master bond (the ray
        # fraction is preserved by the rigid transform). We hit-test against the
        # master bond geometry directly (closest_cylinder_intercept, not
        # _bond_intercept) so a ghost stays pickable even when its master bond is
        # hidden by the spotlight - matching how SymAtomsDrawing picks.
        for k in numpy.unique(bond_syms):
            k = int(k)
            bset = bonds[bond_syms == k]
            if not len(bset):
                continue
            lxyz1, lxyz2 = Place(matrix=tfs[k]).inverse() * (mxyz1, mxyz2)
            a1, a2 = bset.atoms
            f, bnum = closest_cylinder_intercept(a1.coords, a2.coords, bset.radii,
                                                 lxyz1, lxyz2)
            if f is not None and (best_f is None or f < best_f):
                best_bond, best_f, best_k = bset[bnum], f, k
        if best_bond is None:
            return None
        return PickedSymBond(best_bond, best_f, symops[best_k])

    def planes_pick(self, mxyz1, mxyz2, exclude=None):
        return []

    def bounds(self, positions=True):
        #TODO: properly calculate bounds
        return self.geometry_bounds()

    def select(self, mode = 'add'):
        pass

    def update_selection(self):
        pass

class PickedSymBond(Pick):
    def __init__(self, bond, distance, sym):
        super().__init__(distance)
        self.bond = bond
        self.sym = sym
    def description(self):
        return '({}) {}'.format(self.sym.format_as_symop(), self.bond)

class SymRibbonDrawing(Drawing):
    pickable = False
    def __init__(self, name, manager, dim_factor):
        super().__init__(name)
        self._manager = manager
        # m = self._master = master_ribbon
        self._tether_coords = numpy.array([], numpy.double)
        _copy_ribbon_drawing(self.master_ribbon, self, dim_factor)
        self._dim_factor = dim_factor

    @property
    def structure(self):
        return self._manager.structure

    @property
    def master_ribbon(self):
        return self.structure._ribbons_drawing

    @property
    def dim_factor(self):
        return self._dim_factor

    @dim_factor.setter
    def dim_factor(self, factor):
        self._dim_factor = factor
        self.update()

    def first_intercept(self, mxyz1, mxyz2, exclude=None):
        return None

    def planes_pick(self, myxz1, mxyz2, exclude=None):
        return []

    def delayed_rebuild(self, session):
        self._rebuild_handler = session.triggers.add_handler('frame drawn', self.rebuild)

    def rebuild(self, *args):
        self.remove_all_drawings()
        _copy_ribbon_drawing(self.master_ribbon, self, self.dim_factor)
        self.find_tether_coords()
        if args and args[0] == "frame drawn":
            self._rebuild_handler = None
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER

    def find_tether_coords(self):
        tethers = []
        for d in self.child_drawings():
            if 'ribbon_tethers' in d.name:
                tethers.append(d.positions.array()[:,:,3])
        if len(tethers):
            self._tether_coords = numpy.concatenate(tethers)
        else:
            self._tether_coords = numpy.array([], numpy.double)



    def update(self):
        dim_factor = self.dim_factor
        sad, mad = self.all_drawings(), self.master_ribbon.all_drawings()
        if len(sad) != len(mad):
            self.rebuild()
            return
        for td, md in zip(sad, mad):
            td.set_geometry(md.vertices, md.normals, md.triangles)
            td.display = md.display
            vertex_colors = md.vertex_colors
            if vertex_colors is not None:
                vertex_colors = vertex_colors.astype(numpy.float32)
                vertex_colors[:,:3] *= dim_factor
                td.vertex_colors = vertex_colors.astype(numpy.uint8)
            else:
                colors = md.colors.astype(numpy.float32)
                colors[:,:3] *= dim_factor
                td.colors = colors
        self.find_tether_coords()


def _atoms_only_hidden_by_clipper(atoms):
    hides = atoms.hides
    return numpy.logical_and(atoms.displays,
        numpy.logical_not(numpy.logical_and(hides&HIDE_ISOLDE, hides&~HIDE_ISOLDE)))

def _bonds_only_hidden_by_clipper(bonds):
    atoms = bonds.atoms
    return numpy.logical_and(bonds.displays, numpy.logical_and(
        _atoms_only_hidden_by_clipper(atoms[0]),
        _atoms_only_hidden_by_clipper(atoms[1])
    ))


def focus_on_selection(session, view, atoms, clip = True):
    v = view
    pad = 5.0
    bounds = atoms.scene_bounds
    bounds.xyz_min = bounds.xyz_min - pad
    bounds.xyz_max = bounds.xyz_max + pad
    radius = bounds.radius() + pad
    cofr_method = v.center_of_rotation_method
    v.view_all(bounds)
    v.center_of_rotation = center = bounds.center()
    v.center_of_rotation_method = cofr_method
    cam = v.camera
    vd = cam.view_direction()
    if clip:
        cp = v.clip_planes
        cp.set_clip_position('near', center - radius*vd, v)
        cp.set_clip_position('far', center + radius*vd, v)
    session.selection.clear()
    atoms.selected=True


def _copy_ribbon_drawing(master_drawing, target_drawing, dim_factor):
    from chimerax.graphics import Drawing
    d = master_drawing
    t = target_drawing
    def recursively_add_drawings(fromd, tod, dim_factor):
        children = fromd.child_drawings()
        for c in children:
            nc = Drawing(c.name)
            v, n, t, p, d = c.vertices, c.normals, c.triangles, c.positions, c.display
            if v is None or len(v) == 0:
                continue
            nc.set_geometry(v, n, t)
            nc.positions, nc.display = p, d
            vc = c.vertex_colors
            if vc is not None:
                vc = vc.astype(numpy.float32)
                vc[:,:3] *= dim_factor
                nc.vertex_colors = vc.astype(numpy.uint8)
            else:
                col = c.colors.astype(numpy.float32)
                col[:, :3] *= dim_factor
                nc.colors = col.astype(numpy.uint8)
            tod.add_drawing(nc)
            recursively_add_drawings(c, nc, dim_factor)
    recursively_add_drawings(d, t, dim_factor)

def _get_ca_pbg(m):
    from chimerax.atomic import PseudobondGroup
    for cm in m.child_models():
        if isinstance(cm, PseudobondGroup) and m.name == "CA trace":
            return cm
    pbg = m.pseudobond_group("CA trace")
    pbg.dashes = 1
    return pbg

def create_ca_trace(m):
    pbg = _get_ca_pbg(m)
    pbg.clear()
    chains = m.polymers(missing_structure_treatment = m.PMS_NEVER_CONNECTS)
    for chain in chains:
        chain = chain[0]
        if not chain[0].polymer_type == Residue.PT_AMINO:
            continue
        for i in range(len(chain)-1):
            r1, r2 = chain[i:i+2]
            try:
                ca1 = r1.atoms[r1.atoms.names=='CA'][0]
                ca2 = r2.atoms[r2.atoms.names=='CA'][0]
                pb = pbg.new_pseudobond(ca1, ca2)
                pb.color = ca1.color
            except:
                continue
