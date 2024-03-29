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

from chimerax.core.errors import UserError

def open_structure_factors(session, path, structure_model = None,
        over_sampling=2.0, always_raise_errors=True, auto_free_flags=True):
    if structure_model is None:
        raise UserError('Reflection data must be associated with an atomic '
            'structure, provided via the structure_model argument.')
    from .symmetry import get_map_mgr
    mmgr = get_map_mgr(structure_model, create=True)
    try:
        xmapset = mmgr.add_xmapset_from_file(path, oversampling_rate=over_sampling,
            auto_choose_free_flags=auto_free_flags)
        log_str = 'Opened crystallographic dataset from {}\n'.format(path)
        if xmapset.experimental_data:
            log_str += 'Found experimental reflection data: \n'
            log_str += '\n'.join(['\t{}'.format(n) for n in xmapset.experimental_data.keys()])
            log_str += '\n'
            log_str += 'Rwork: {:.4f}; Rfree: {:.4f}\n'.format(
                xmapset.rwork, xmapset.rfree
            )
        log_str += 'Generated maps: \n{}\n'.format(
            '\n'.join(['\t{}'.format(m.name) for m in xmapset]))
        log_str += 'Any unwanted maps may be safely closed via the Model panel.'
        cmgr = mmgr.crystal_mgr
        if cmgr not in session.models.list():
            return [mmgr.crystal_mgr], log_str
        return [], log_str

    except RuntimeError as e:
        if always_raise_errors:
            raise e
        else:
            session.logger.warning(str(e))
            return None, None

def open_structure_factors_and_add_to_session(session, path, structure_model=None,
        over_sampling=2.0, always_raise_errors=False, auto_free_flags=True):
    models, log_str = open_structure_factors(session, path, structure_model,
        over_sampling, always_raise_errors, auto_free_flags)
    if models is not None:
        session.models.add(models)

def save_structure_factors(session, path, models=None, preserve_input=False,
        save_map_coeffs=False):
    if models is None:
        models = session.models.list()
    from .maps import XmapSet
    xmapsets = [m for m in models if isinstance(m, XmapSet)]
    if not xmapsets:
        raise UserError('You need to specify at least one crystallographic map set to save!')
    suffix=''
    import os
    fname, ext = os.path.splitext(path)
    if ext.lower() != '.mtz':
        if len(ext)<=1:
            ext='.mtz'
        else:
            ext += '.mtz'
    suffix=''
    for i, xmapset in enumerate(xmapsets):
        if len(xmapsets) > 1:
            suffix = '-{}'.format(i)
        save_path = fname+suffix+ext
        xmapset.save_mtz(save_path, save_input_data=preserve_input,
            save_output_fobs=True, save_map_coeffs=save_map_coeffs)




def spotlight(session, models=None, enable=True, radius=None):
    from chimerax.clipper.symmetry import get_symmetry_handler, SymmetryManager
    if models is None:
        from chimerax.atomic import AtomicStructure
        models = session.models.list(type=AtomicStructure)
        sym_mgrs = session.models.list(type=SymmetryManager)
        if len(sym_mgrs) == 0:
            from chimerax.core.errors import UserError
            raise UserError('No models are initialised for Clipper. Use the command "clipper spotlight #{model id(s)}" to initialise (a) specific model(s).')
        models = [m for m in models if m.parent not in sym_mgrs]
        models = models+sym_mgrs
        create = False
    else:
        create = True
    for m in models:
        if isinstance(m, SymmetryManager):
            sh = m
        else:
            sh = get_symmetry_handler(m, create=create, auto_add_to_session=True)
        if sh is not None:
            sh.spotlight_mode=enable
            if radius is not None:
                sh.spotlight_radius = radius



def associate_volumes(session, volumes, to_model=None):
    if to_model is None:
        from chimerax.core.errors import UserError
        raise UserError('The toModel argument must be provided!')
    from chimerax.clipper.symmetry import get_map_mgr
    mgr = get_map_mgr(to_model, create=True)
    for v in volumes:
        mgr.nxmapset.add_nxmap_handler_from_volume(v)
    if mgr.crystal_mgr not in session.models.list():
        session.models.add([mgr.crystal_mgr])

def isolate(session, atoms,
        surround_distance=0,
        context_distance=5,
        mask_radius=3,
        hide_surrounds=True,
        focus=False):
    from chimerax.clipper.symmetry import get_symmetry_handler
    us = atoms.unique_structures
    for s in us:
        sel = us.atoms.intersect(atoms)
        sh = get_symmetry_handler(s)
        if sh is not None:
            sh.isolate_and_cover_selection(sel,
                include_surrounding_residues = surround_distance,
                show_context = context_distance,
                mask_radius = mask_radius,
                hide_surrounds = hide_surrounds,
                focus = focus)

def init_environ(session, mouse_modes=True, cofr=True, camera=True, lighting=True):
    from chimerax.core.commands import run
    if mouse_modes:
        from .mousemodes import initialize_clipper_mouse_modes
        initialize_clipper_mouse_modes(session)
    if cofr:
        run(session, 'cofr center showPivot t')
    if camera:
        run(session, 'camera ortho')
    if lighting:
        run(session, 'lighting simple')

def reset_environ(session, mouse_modes=True, cofr=True, camera=True, clip_planes=True):
    from chimerax.core.commands import run
    if mouse_modes:
        run(session, ('mousemode left select control;' 
                      'mousemode left none control shift;'
                      'mousemode middle none control;'
                      'mousemode wheel zoom; mousemode right none shift;'
                      'mousemode wheel none control;'
                      'mousemode wheel none alt;'
                      'mousemode wheel none shift'
                    )
        ) 
    if cofr:
        run(session, 'cofr front showPivot f')
    if camera:
        run(session, 'camera mono')
    if clip_planes:
        run(session, 'clip off')

from chimerax.core.commands.atomspec import AtomSpecArg
class VolumesArg(AtomSpecArg):
    """Parse command models specifier"""
    name = "a models specifier"

    @classmethod
    def parse(cls, text, session):
        '''
        Returns only Volume objects (not subclasses)
        '''
        from chimerax.map import Volume
        aspec, text, rest = super().parse(text, session)
        models = aspec.evaluate(session).models
        volumes = [m for m in models if type(m) == Volume]
        return volumes, text, rest

class XmapSetsArg(AtomSpecArg):
    name = "a models specifier"

    @classmethod
    def parse(cls, text, session):
        '''
        Returns any XmapSet objects found in the selection
        '''
        from .maps import XmapSet
        aspec, text, rest = super().parse(text, session)
        models = aspec.evaluate(session).models
        xmapsets = [m for m in models if type(m) == XmapSet]
        return xmapsets, text, rest

class AtomicStructuresOrSymmetryMgrsArg(AtomSpecArg):
    name = "a models specifier"

    @classmethod
    def parse(cls, text, session):
        '''
        Returns any AtomicStructures or SymmetryManagers found in the selection.
        If an AtomicStructure is already managed by a SymmetryManager, only the
        SymmetryManager will be included
        '''
        from .symmetry import SymmetryManager
        from chimerax.atomic import AtomicStructure
        aspec, text, rest = super().parse(text, session)
        models = aspec.evaluate(session).models
        sym_mgrs = [m for m in models if isinstance(m, SymmetryManager)]
        structures = [m for m in models if type(m)==AtomicStructure and
            m.parent not in sym_mgrs]
        return structures+sym_mgrs, text, rest


def draw_symmetry(session, models):
    logger = session.logger
    from chimerax.clipper.symmetry import SymmetryManager
    from chimerax.clipper.geometry import unit_cell_and_sym_axes
    for m in models:
        if isinstance(m.parent, SymmetryManager):
            unit_cell = m.parent.unit_cell
            has_symmetry_manager=True
        else:
            has_symmetry_manager=False
            from chimerax.clipper.symmetry import symmetry_from_model_metadata, Unit_Cell
            try:
                cell, spacegroup, grid, resolution, has_symmetry = symmetry_from_model_metadata(m)
            except Exception as e:
                warn_str = ('Failed to read symmetry information for model {}. '
                'Attempt returned the following error message: \n{}').format(
                    m.id_string, str(e)
                )
                logger.warning(warn_str)
                continue
            if not has_symmetry:
                logger.info('Model {} has no symmetry information. Skipping.').format(m.id_string)
                continue
            unit_cell = Unit_Cell(m.atoms, cell, spacegroup, grid)
        sym_drawing = unit_cell_and_sym_axes(session, unit_cell)
        if has_symmetry_manager:
            # Add the drawing to the symmetry manager so it's bundled with the model
            m.parent.add([sym_drawing])
        else:
            # Just add it at the top level
            sym_drawing.name += '({})'.format(m.id_string)
            session.models.add([sym_drawing])


def set_contour_sensitivity(session, sensitivity):
    from .mousemodes import ContourSelectedVolume
    mm = [b.mode for b in session.ui.mouse_modes.bindings if isinstance(b.mode, ContourSelectedVolume)]
    if len(mm):
        mm = mm[0]
    else:
        from chimerax.core.errors import UserError
        raise UserError('Clipper is not yet initialised. Command ignored.')
    mm.sensitivity = sensitivity

def set_camera_auto(session, flag):
    from . import symmetry
    symmetry.auto_reset_camera = flag

def register_clipper_cmd(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        BoolArg, FloatArg,
        OpenFileNameArg, SaveFileNameArg,
        create_alias
        )
    from chimerax.atomic import StructuresArg, StructureArg, AtomsArg

    init_desc = CmdDesc(
        keyword=[
            ('mouse_modes', BoolArg),
            ('cofr', BoolArg),
            ('camera', BoolArg),
            ('lighting', BoolArg)
        ],
        synopsis = 'Initialise the Clipper environment (mouse modes, center of rotation, camera and/or lighting) without affecting existing models'
    )
    register('clipper init', init_desc, init_environ, logger=logger)

    reset_desc = CmdDesc(
        keyword=[
            ('mouse_modes', BoolArg),
            ('cofr', BoolArg),
            ('camera', BoolArg),
            ('clip_planes', BoolArg)
        ],
        synopsis='Reset mouse modes, cofr, camera, and/or clip planes to ChimeraX defaults'
    )
    register('clipper off', reset_desc, reset_environ, logger=logger)

    open_desc = CmdDesc(
        required=[
            ('path', OpenFileNameArg),
        ],
        keyword=[
            ('structure_model', StructureArg),
            ('over_sampling', FloatArg),
            ('auto_free_flags', BoolArg)
        ],
        synopsis='Open a structure factor .mtz or .cif file and generate maps for the given model'
    )
    register('clipper open', open_desc, open_structure_factors_and_add_to_session, logger=logger)

    save_desc = CmdDesc(
        required=[
            ('path', SaveFileNameArg),
        ],
        optional=[
            ('models', XmapSetsArg),
        ],
        keyword=[
            ('preserve_input', BoolArg),
            ('save_map_coeffs', BoolArg)
        ],
        synopsis='Save structure factors to .mtz format.'
    )
    register('clipper save', save_desc, save_structure_factors)

    spot_desc = CmdDesc(
        optional=[
            ('models', AtomicStructuresOrSymmetryMgrsArg),
        ],
        keyword=[
            ('radius', FloatArg),
        ],
        synopsis='Switch on/off "Scrolling sphere" visualisation with live atomic symmetry'
    )
    register('clipper spotlight', spot_desc, spotlight, logger=logger)

    vol_desc = CmdDesc(
        required=[
            ('volumes', VolumesArg),
        ],
        keyword=[
            ('to_model', StructureArg),
        ],
        synopsis='Have Clipper take control of the chosen volumes and associate them with the given model'
    )
    register('clipper associate', vol_desc, associate_volumes, logger=logger)

    isol_desc = CmdDesc(
        required=[('atoms', AtomsArg)],
        keyword=[
            ('surround_distance', FloatArg),
            ('context_distance', FloatArg),
            ('mask_radius', FloatArg),
            ('hide_surrounds', BoolArg),
            ('focus', BoolArg),
        ],
        synopsis=('Visually isolate the selected atoms from their surroundings, '
            'and mask their maps to their immediate vicinity. The selection '
            'covered by the map(s) will be expanded to include all residues '
            'approaching within surroundDistance of the given selection. Any '
            'residues approaching within contextDistance of the result will be '
            'displayed, but not covered by the map(s). If hideSurrounds is '
            'True, all other atoms will be hidden. If focus is True, the view '
            'will be reset to cover the visible atoms. To revert to the default '
            'viewing mode, use "clipper spotlight".')
    )
    register('clipper isolate', isol_desc, isolate, logger=logger)

    sym_desc = CmdDesc(
        required=[('models', StructuresArg)],
        synopsis=('Draw the unit cell and all symmetry axes for the given '
            'crystal structures. NOTE: this tool is still a work in progress.'
        )
    )
    register('clipper symmetry', sym_desc, draw_symmetry, logger=logger)

    set_contour_sensitivity_desc = CmdDesc(
        required=[('sensitivity', FloatArg),],
        synopsis=('Set the sensitivity of map contouring to the mouse scroll '
            'wheel. Each "tick" of the wheel will change the contour by '
            '(sensitivity * map standard deviation).')
    )
    register('clipper set contourSensitivity', set_contour_sensitivity_desc, set_contour_sensitivity)

    set_camera_auto_desc = CmdDesc(
        required=[('flag', BoolArg)],
        synopsis='Tell Clipper whether to automatically reset the camera to orthographic view when switching to spotlight mode'
    )
    register('clipper set autoCameraOrtho', set_camera_auto_desc, set_camera_auto, logger=logger)


def register_cview_cmd(logger):
    from chimerax.core.commands import create_alias
    create_alias('cview', 'view $*; cofr center showpivot true', logger=logger)
