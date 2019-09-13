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
        over_sampling=2.0):
    if structure_model is None:
        raise UserError('Reflection data must be associated with an atomic '
            'structure, provided via the structure_model argument.')
    from .symmetry import get_map_mgr
    mmgr = get_map_mgr(structure_model, create=True)
    try:
        xmapset = mmgr.add_xmapset_from_mtz(path, oversampling_rate=over_sampling)
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
        return [mmgr.crystal_mgr], log_str

    except RuntimeError as e:
        if always_raise_errors:
            raise e
        else:
            session.logger.warning(str(e))
            return None, None

def open_structure_factors_and_add_to_session(session, path, structure_model=None,
        over_sampling=2.0, always_raise_errors=False):
    models, log_str = open_structure_factors(session, path, structure_model,
        over_sampling, always_raise_errors)
    session.models.add(models)

def save_structure_factors(session, path, models=None, preserve_input=False,
        save_map_coeffs=False):
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
    for i, xmapset in enumerate(xmapsets):
        if i > 0:
            suffix = '-{}'.format(i)
        save_path = fname+suffix+ext
        xmapset.save_mtz(save_path, save_input_data=preserve_input,
            save_output_fobs=True, save_map_coeffs=save_map_coeffs)




def spotlight(session, models=None, enable=True, radius=None):
    from chimerax.clipper.symmetry import get_symmetry_handler
    if models is None:
        from chimerax.atomic import AtomicStructure
        models = session.models.list(type=AtomicStructure)
    for m in models:
        sh = get_symmetry_handler(m, create=True, auto_add_to_session=True)
        session.logger.info('Setting spotlight mode for model {} to {}'.format(
            m.id_string, enable
        ))
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
    session.models.add([mgr.crystal_mgr])

def isolate(session, atoms,
        surround_distance=5,
        context_distance=5,
        mask_radius=3,
        hide_surrounds=True,
        focus=False,
        include_symmetry=True):
    from chimerax.clipper.symmetry import get_symmetry_handler
    us = atoms.unique_structures
    for s in us:
        sel = us.atoms.intersect(atoms)
        sh = get_symmetry_handler(s, create=True)
        sh.isolate_and_cover_selection(sel,
            include_surrounding_residues = surround_distance,
            show_context = context_distance,
            mask_radius = mask_radius,
            hide_surrounds = hide_surrounds,
            focus = focus,
            include_symmetry = include_symmetry)




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


def register_clipper_cmd(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        BoolArg, FloatArg,
        OpenFileNameArg
        )
    from chimerax.atomic import StructuresArg, StructureArg, AtomsArg

    open_desc = CmdDesc(
        required=[
            ('path', OpenFileNameArg),
        ],
        keyword=[
            ('structure_model', StructureArg),
            ('over_sampling', FloatArg)
        ],
        synopsis='Open a structure factor .mtz or .cif file and generate maps for the given model'
    )
    register('clipper open', open_desc, open_structure_factors_and_add_to_session, logger=logger)

    spot_desc = CmdDesc(
        optional=[
            ('models', StructuresArg),
            ('enable', BoolArg)
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
            ('include_symmetry', BoolArg)
        ],
        synopsis=('Visually isolate the selected atoms from their surroundings, '
            'and mask their maps to their immediate vicinity. The selection '
            'covered by the map(s) will be expanded to include all residues '
            'approaching within surroundDistance of the given selection. Any '
            'residues approaching within contextDistance of the result will be '
            'displayed, but not covered by the map(s). If hideSurrounds is '
            'True, all other atoms will be hidden. If focus is True, the view '
            'will be reset to cover the visible atoms. If includeSymmetry is '
            'True, symmetry atoms will be included in the contextDistance '
            'calculation. To revert to the default viewing mode, use '
            '"clipper spotlight".')
    )
    register('clipper isolate', isol_desc, isolate, logger=logger)
