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

def _radiation_from_structure(structure_model):
    '''Radiation implied by the model's mmCIF _exptl.method, or None.

    Deposited PDB entries record the experiment as "ELECTRON CRYSTALLOGRAPHY" vs
    "X-RAY DIFFRACTION" in _exptl.method - the most reliable radiation signal, and
    the natural one for `open <pdbid> structureFactors true`. Uses ChimeraX's own
    metadata accessor rather than re-parsing the header.'''
    if structure_model is None:
        return None
    try:
        from chimerax.mmcif import get_mmcif_tables_from_metadata
        tbl = get_mmcif_tables_from_metadata(structure_model, ['exptl'])[0]
        if not tbl:
            return None
        method = tbl.fields(['method'], allow_missing_fields=True)[0][0]
    except Exception:
        return None
    if method and 'electron' in method.lower():
        return 'electron'
    if method and 'x-ray' in method.lower():
        return 'xray'
    return None


def _resolve_macro_radiation(radiation, path, structure_model=None, logger=None):
    '''Resolve the macromolecular `radiation` argument to 'xray' or 'electron'.

    Explicit choices pass through. 'auto' consults, in order of reliability, the
    model's mmCIF _exptl.method, then a positive electron signal from a
    structure-factor CIF's _diffrn_radiation_probe. When nothing declares the
    experiment (e.g. a PDB-format model with no mmCIF metadata, or an mmCIF with
    no _exptl loop) it defaults to X-ray and, if a logger is given, warns - so a
    silent wrong default on electron data is at least visible and overridable.'''
    r = (radiation or 'auto').lower()
    if r in ('xray', 'electron'):
        return r
    # Model _exptl.method is authoritative when present and recognised.
    from_model = _radiation_from_structure(structure_model)
    if from_model is not None:
        return from_model
    # Model didn't declare a recognised experiment; trust only a *positive*
    # electron probe from a structure-factor CIF (its 'xray' is a bare default).
    if str(path).lower().endswith('.cif'):
        try:
            from .io.small_molecule import _radiation_from_cif
            if _radiation_from_cif(path) == 'electron':
                return 'electron'
        except Exception:
            pass
    if logger is not None:
        logger.warning(
            '(CLIPPER) Could not determine the experiment type from the model '
            'metadata (no mmCIF _exptl.method); defaulting to X-ray scattering '
            'factors. If this is electron-diffraction (micro-ED / 3D-ED) data, '
            're-open with "radiation electron".')
    return 'xray'


def open_structure_factors(session, path, structure_model = None,
        over_sampling=2.0, always_raise_errors=True, auto_free_flags=True,
        browse=False, fsigf=None, map_columns=None, free_flags=None,
        free_flags_file=None, radiation='auto'):
    if structure_model is None:
        raise UserError('Reflection data must be associated with an atomic '
            'structure, provided via the structure_model argument.')
    from .symmetry import get_map_mgr
    mmgr = get_map_mgr(structure_model, create=True)
    try:
        # radiation ('auto'|'xray'|'electron') is resolved centrally in
        # add_xmapset_from_file ('auto' reads the model's mmCIF _exptl.method).
        xmapset = mmgr.add_xmapset_from_file(path, oversampling_rate=over_sampling,
            auto_choose_free_flags=auto_free_flags,
            free_flag_label=free_flags, free_flags_file=free_flags_file,
            fsigf_name=fsigf, map_columns=map_columns, browse=browse,
            radiation=radiation)
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
        over_sampling=2.0, always_raise_errors=False, auto_free_flags=True,
        browse=False, fsigf=None, map_columns=None, free_flags=None,
        free_flags_file=None, radiation='auto'):
    models, log_str = open_structure_factors(session, path, structure_model,
        over_sampling, always_raise_errors, auto_free_flags,
        browse=browse, fsigf=fsigf, map_columns=map_columns,
        free_flags=free_flags, free_flags_file=free_flags_file, radiation=radiation)
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

def debug_maps(session, models=None, enable=True):
    '''Add (or remove) live structure-factor component maps -- Fcalc, Fatoms,
    Fbulk and Fmask -- for the given crystallographic model(s), as a debugging
    aid. Each is a live map that FFTs a raw component rather than sigmaa-weighted
    coefficients, and is display-gated: created hidden, recomputed only while
    shown. Not saved in sessions.'''
    from chimerax.clipper.symmetry import get_map_mgr
    from chimerax.core.errors import UserError
    if models is None:
        models = session.models.list()
    mgrs = []
    seen = set()
    for m in models:
        try:
            mgr = get_map_mgr(m, create=False)
        except Exception:
            mgr = None
        if mgr is not None and id(mgr) not in seen:
            seen.add(id(mgr))
            mgrs.append(mgr)
    if not mgrs:
        raise UserError('No Clipper crystallographic maps found for the specified '
            'model(s). Open reflection data first with "clipper open".')
    n_sets = 0
    for mgr in mgrs:
        for xs in mgr.xmapsets:
            if xs.live_xmap_mgr is None:
                continue
            if enable:
                xs.add_debug_maps()
            else:
                xs.remove_debug_maps()
            n_sets += 1
    if enable:
        session.logger.info('Added Fcalc/Fatoms/Fbulk/Fmask component maps to {} '
            'crystallographic dataset(s). They are hidden by default; show one to '
            'view it (each is recomputed only while shown).'.format(n_sets))
    else:
        session.logger.info('Removed structure-factor component maps from {} '
            'dataset(s).'.format(n_sets))

def set_scaling(session, models=None, all_reflections=None, reflections_per_bin=None,
                num_bins=None, seed=None):
    '''Control how the live crystallographic maps fit the Fcalc->Fobs bulk-solvent
    scale. By default the scale is fit on a fixed, seeded subset of reflections (fast
    and reproducible, with a tiny fixed bias). "allReflections true" fits over all
    reflections instead -- exact and unbiased, at a larger (one-off) cost.
    reflectionsPerBin / numBins set the size of the seeded subset (defaults 500 and
    20); seed sets which subset is drawn (default 10061865 -- vary it to probe
    sensitivity). All three are ignored when allReflections is true. With no options,
    reports the current settings.'''
    from chimerax.clipper.symmetry import get_map_mgr
    from chimerax.core.errors import UserError
    if models is None:
        models = session.models.list()
    mgrs = []
    seen = set()
    for m in models:
        try:
            mgr = get_map_mgr(m, create=False)
        except Exception:
            mgr = None
        if mgr is not None and id(mgr) not in seen:
            seen.add(id(mgr))
            mgrs.append(mgr)
    if not mgrs:
        raise UserError('No Clipper crystallographic maps found for the specified '
            'model(s). Open reflection data first with "clipper open".')
    n_sets = 0
    for mgr in mgrs:
        for xs in mgr.xmapsets:
            if xs.live_xmap_mgr is None:
                continue
            if all_reflections is not None:
                xs.all_reflections = all_reflections
            if reflections_per_bin is not None:
                xs.scaling_reflections_per_bin = reflections_per_bin
            if num_bins is not None:
                xs.scaling_num_bins = num_bins
            if seed is not None:
                xs.scaling_seed = seed
            n_sets += 1
            if xs.all_reflections:
                how = 'all reflections (exact)'
            else:
                how = 'seeded subset (seed={}, {}/bin x {} bins)'.format(
                    xs.scaling_seed, xs.scaling_reflections_per_bin, xs.scaling_num_bins)
            session.logger.info(
                '(CLIPPER) {}: bulk-solvent scaling = {}'.format(xs.name, how))
    if n_sets == 0:
        raise UserError('No live crystallographic maps found for the specified '
            'model(s) (static maps have no scaling to control).')

def discard_symmetry(session, models):
    '''Discard a model's crystallographic symmetry (revert to a P1 box) and
    rewrite its symmetry metadata so this is preserved across sessions. Useful
    when a model's original symmetry is irrelevant to its current use (e.g. an
    old crystal structure reused as a starting model in a cryo-EM map).'''
    from chimerax.clipper.symmetry import get_symmetry_handler, SymmetryManager
    for m in models:
        sh = m if isinstance(m, SymmetryManager) else get_symmetry_handler(m)
        if sh is None:
            session.logger.warning('No Clipper symmetry manager for {}; skipping.'.format(m))
            continue
        sh.discard_model_symmetry()

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

def _choose_bsharp_targets(named_maps, broad_live_maps):
    '''
    Pure target-selection for "clipper bsharp" (no session/display dependency, so
    it is unit-testable). Returns (targets, locked): `targets` is the list of maps
    to change; `locked` is set to the single explicitly-named map when it exists
    but is non-adjustable (so the caller can raise a clear error rather than
    silently doing nothing).
    '''
    if len(named_maps) == 1 and not broad_live_maps:
        h = named_maps[0]
        if not h.bsharp_adjustable:
            return [], h
        return [h], None
    # Broad or multi-map scope: only the designated viewing map(s). Explicitly
    # named maps are folded in but still filtered to designated + adjustable.
    pool = list(named_maps) + list(broad_live_maps)
    return [h for h in pool if h.bsharp_adjustable and h.bsharp_designated], None

def set_map_bsharp(session, b_sharp, maps=None):
    '''
    Set the sharpening (positive) / smoothing (negative) B-factor applied to live
    crystallographic maps. Affected maps are recomputed immediately from their
    cached coefficients (no atom movement, no Fcalc regeneration).

    Targeting:
      * one explicitly-named adjustable map  -> that map is changed;
      * no map given, or a whole model/dataset, or several maps -> only the
        "designated" viewing map(s) are changed (the heuristic sharpened map,
        not the raw reference or difference maps);
      * difference maps and maps explicitly locked (e.g. ISOLDE's MDFF map) are
        never changed.
    '''
    from .maps import XmapSet
    from .maps.xmapset import XmapHandler_Live

    # Split the argument into individually-named maps vs broad scopes (whole
    # XmapSets / structure models), so we can distinguish "change this one map"
    # from "retune the designated map of this dataset".
    named_maps = []
    broad_sets = []
    if maps is None:
        broad_sets = list(session.models.list(type=XmapSet))
    else:
        for m in maps:
            if isinstance(m, XmapHandler_Live):
                named_maps.append(m)
            elif isinstance(m, XmapSet):
                broad_sets.append(m)
            else:
                # A structure/top-level model: find its associated XmapSet(s).
                broad_sets.extend(mm for mm in m.all_models()
                    if isinstance(mm, XmapSet))

    broad_live_maps = []
    for xs in broad_sets:
        broad_live_maps.extend(xs.live_xmaps)

    targets, locked = _choose_bsharp_targets(named_maps, broad_live_maps)
    if locked is not None:
        raise UserError('Map "{}" cannot be sharpened/smoothed (it is a '
            'difference map or a locked map such as an MDFF map).'.format(locked.name))

    if not targets:
        raise UserError('No adjustable live crystallographic maps found to apply '
            'the sharpening B-factor to. Specify a single live map to change it '
            'directly, or open a dataset with "clipper open" first.')
    for h in targets:
        h.b_sharp = b_sharp
    session.logger.info('Set sharpening B-factor to {:.1f} for {} live map(s): {}'.format(
        b_sharp, len(targets), ', '.join(h._map_name for h in targets)))

def set_oversampling(session, oversampling, models=None):
    '''
    Set the real-space map oversampling (Shannon) rate at runtime for one or more
    crystallographic map managers. Re-grids all live and static maps and rebuilds
    the symmetry Unit_Cell to match; the reflection data is unchanged.
    '''
    from .maps import XmapSet
    from .maps.map_mgr import MapMgr
    from .symmetry import get_map_mgr
    mmgrs = []
    def add(mm):
        if mm is not None and mm not in mmgrs:
            mmgrs.append(mm)
    if models is None:
        for mm in session.models.list(type=MapMgr):
            add(mm)
    else:
        for m in models:
            if isinstance(m, MapMgr):
                add(m)
            elif isinstance(m, XmapSet):
                add(m.master_map_mgr)
            else:
                add(get_map_mgr(m, create=False))
    if not mmgrs:
        raise UserError('No crystallographic map manager found. Open a dataset '
            'with "clipper open" first, or specify a model that has one.')
    for mm in mmgrs:
        mm.set_oversampling_rate(oversampling)

def open_small_molecule(session, path, hkl=None, radiation='auto', fragments='rename'):
    '''
    Open a small-molecule (e.g. COD) CIF as a live crystal structure: the model in
    its unit cell, crystallographic symmetry, and - when reflections are available
    (an `hkl` file, or a sibling .hkl) - live 2mFo-DFc / mFo-DFc electron-density
    maps that update as the model changes.

    radiation: 'xray', 'electron' (micro-ED / 3D-ED), or 'auto' (default; read from
    the CIF _diffrn_radiation_probe) - selects the scattering-factor table.

    fragments: 'off', 'rename' (default) or 'complete' - split the asymmetric unit
    into named covalent-fragment residues (see the `clipper fragments` command).
    'complete' additionally reassembles whole molecules via crystallographic symmetry:
    it rewraps symmetry-scattered hydrogens onto their parents and adds the atoms that
    finish molecules split across a special position.
    '''
    from .io.small_molecule import show_cod_crystal
    return show_cod_crystal(session, path, hkl_path=hkl, radiation=radiation,
                            fragments=fragments)


def fetch_cod_crystal(session, cod_id, ignore_cache=False, radiation='auto',
                      fragments='rename'):
    '''
    Fetch a structure (and, if deposited, its reflections) from the Crystallography
    Open Database by numeric COD ID, and open it as a live crystal structure.
    '''
    from chimerax.core.errors import UserError
    cod_id = str(cod_id).strip()
    if not cod_id.isdigit():
        raise UserError('COD identifiers are numeric; got "%s"' % cod_id)
    from chimerax.core.fetch import fetch_file
    base = 'https://www.crystallography.net/cod/%s' % cod_id
    cif = fetch_file(session, base + '.cif', 'COD %s' % cod_id,
        '%s.cif' % cod_id, 'COD', ignore_cache=ignore_cache)
    hkl = None
    try:
        hkl = fetch_file(session, base + '.hkl', 'COD %s reflections' % cod_id,
            '%s.hkl' % cod_id, 'COD', ignore_cache=ignore_cache)
    except Exception:
        session.logger.info('(CLIPPER) No reflections deposited for COD %s; '
            'showing model + symmetry only.' % cod_id)
    from .io.small_molecule import show_cod_crystal
    return show_cod_crystal(session, cif, hkl_path=hkl, radiation=radiation,
                            fragments=fragments)


def split_fragments_cmd(session, structures, mode='rename'):
    '''
    Split each (small-molecule / corecif) structure's asymmetric unit into named
    covalent-fragment residues: CCD codes for the common simple species (water,
    monatomic ions, small inorganic ions) and LIG01, LIG02, ... otherwise.

    mode 'off' leaves the model unchanged; 'rename' (default) splits and names;
    'complete' reassembles whole molecules via crystallographic symmetry - it rewraps
    symmetry-scattered hydrogens onto their parents and adds the symmetry-generated
    atoms that finish molecules split across a special position.
    '''
    from chimerax.core.errors import UserError
    from .symmetry import crystal_symmetry_from_cif_file
    from .io.fragments import split_fragments
    from .io.small_molecule import _prepare_corecif_model, warn_implausibly_long_bonds
    if structures is None or not len(structures):
        raise UserError('No structures specified.')
    if mode == 'off':
        return []
    results = []
    for m in structures:
        if not getattr(m, 'is_corecif', False):
            session.logger.warning('(CLIPPER) %s is not a small-molecule (corecif) '
                'model; skipping fragment split.' % m)
            continue
        # A live small-molecule map builds its structure-factor atom list fresh from
        # model.atoms on every recompute (no frozen scaffold), so splitting - even
        # 'complete', which adds atoms - is safe while a map is open: the map simply
        # reflects the new atoms/residues on its next update.
        path = getattr(m, 'filename', None)
        if path is None:
            session.logger.warning('(CLIPPER) %s has no source file; cannot recover '
                'crystal symmetry for fragment split; skipping.' % m)
            continue
        cell, spacegroup, grid = crystal_symmetry_from_cif_file(path)
        # Apply the corecif workarounds to this (generically-opened) model: rebuild
        # Clipper-frame coordinates - corecif mis-orthogonalises oblique cells, which
        # would misplace completed symmetry atoms and defeat the special-position test -
        # and repair covalent connectivity corecif drops on metal-coordinated atoms.
        # Both idempotent, so this is safe on a model opened via `clipper smallmol` too.
        _prepare_corecif_model(session, m, path, cell)
        # 'complete' reassembles molecules via symmetry; rewrap symmetry-scattered H onto
        # their parents first (matches the `clipper cod`/`smallmol` complete path).
        if mode == 'complete':
            from .io.small_molecule import reassemble_symmetry_scattered_hydrogens
            n = reassemble_symmetry_scattered_hydrogens(session, m, cell, spacegroup)
            if n:
                session.logger.info('(CLIPPER) %s: reassembled %d symmetry-scattered '
                    'hydrogen(s).' % (m, n))
        results.append(split_fragments(session, m, cell, spacegroup, grid, mode=mode,
                                       path=path))
        # Flag (not auto-fix) implausibly long bonds corecif copied from the CIF
        # _geom_bond loop (e.g. a 1,3 distance listed as a bond) - connectivity to review.
        import os
        warn_implausibly_long_bonds(session, m, os.path.basename(path))
    return results


def _site_multiplicities_for(sym, structure):
    '''Per-atom site multiplicity of a structure's ASU, for special-position dedup.'''
    from .clipper_util import site_multiplicities
    return site_multiplicities(structure.atoms.coords, sym.cell, sym.spacegroup,
        sym.grid)


def symmetry_copies_real(session, models=None, contacting=None, distance=5.0,
        name=None, focus=False, refl_file=None, prune_special_positions=True):
    '''
    Create real, whole-model copies of crystallographic symmetry mates and merge
    them (together with the original ASU) into a single new model. Works both in
    the GUI and headless (--nogui): the crystallographic core is obtained without
    building the GUI symmetry framework (see clipper.crystal_symmetry_for).

    Selection criteria:
      * contacting given -> every symmetry copy with an atom within `distance` of
        the `contacting` selection is realised (per structure the selection
        touches);
      * contacting omitted, GUI spotlight active -> every symmetry copy currently
        drawn in Clipper's spotlight (i.e. whose ribbon is displayed) is realised;
      * contacting omitted, no spotlight (e.g. headless) -> every symmetry copy
        within `distance` of the whole model is realised.

    `refl_file` (a .mtz/.cif structure-factor file) supplies the cell and
    spacegroup for models whose own metadata lacks them. When
    `prune_special_positions` is True (default), redundant copies of atoms sitting
    on special positions are removed using exact Clipper site multiplicity. The
    result carries a `clipper_sym_expansion` recording which operator produced
    which chains, so the expansion can be inverted (see clipper.collapse_to_asu).
    '''
    from .symmetry import get_symmetry_handler, get_all_symmetry_handlers
    from .sym_realize import crystal_symmetry_for, realize_symmetry_copies
    cutoff_mode = contacting is not None and len(contacting) > 0
    if cutoff_mode:
        structures = list(contacting.unique_structures)
        if models is not None:
            structures = [s for s in structures if s in models]
    elif models is not None:
        structures = list(models)
    else:
        structures = [sh.structure for sh in get_all_symmetry_handlers(session)]

    if not structures:
        raise UserError('No crystal structures found to make symmetry copies for. '
            'Specify one or more models, or a "contacting" atom selection.')

    results = []
    for s in structures:
        sym = crystal_symmetry_for(s, refl_file=refl_file)
        if not sym.has_symmetry:
            session.logger.warning('Model {} has no crystallographic symmetry; '
                'skipping.'.format(s.id_string))
            continue
        if cutoff_mode:
            sel = contacting.intersect(s.atoms)
            if not len(sel):
                continue
            places = sym.sym_transforms_near(sel, distance)
            descr = 'within {:.3g} A of the selection'.format(distance)
        else:
            # In the GUI, prefer the spotlight (what the user can see) when it is
            # showing symmetry copies; otherwise - and always when headless - fall
            # back to the whole-model neighbours. Only touch the (GUI-only)
            # AtomicSymmetryModel when a GUI is present.
            places = []
            if session.ui.is_gui:
                sh = get_symmetry_handler(s)
                asm = sh.atomic_symmetry_model if sh is not None else None
                if asm is not None:
                    places = asm.currently_displayed_sym_transforms()
            if len(places) > 1:
                descr = 'currently displayed'
            else:
                places = sym.sym_transforms_near(s.atoms, distance)
                descr = 'within {:.3g} A of the model'.format(distance)
        # places always leads with the identity/ASU operator, so >1 means there is
        # at least one genuine symmetry copy to realise.
        if len(places) <= 1:
            session.logger.warning('No symmetry copies {} for model {}.'.format(
                descr, s.id_string))
            continue
        mults = _site_multiplicities_for(sym, s) if prune_special_positions else None
        combined = realize_symmetry_copies(session, s, places, name=name,
            prune_special_positions=prune_special_positions, multiplicities=mults)
        if combined is None:
            session.logger.warning('All symmetry copies {} for model {} were '
                'redundant special-position duplicates; nothing added.'.format(
                    descr, s.id_string))
            continue
        session.logger.info('Created model {} combining the ASU of {} with its '
            'symmetry copies {}.'.format(combined.id_string, s.id_string, descr))
        results.append(combined)

    if focus and results and session.ui.is_gui:
        from chimerax.core.commands import run
        run(session, 'view {}'.format(' '.join('#'+m.id_string for m in results)))
    return results


def generate_unit_cells(session, models=None, cells=None, box=None, name=None,
        focus=False, refl_file=None, prune_special_positions=True):
    '''
    Generate one or more whole crystallographic unit cells as a single new,
    real atomic model. Each unit cell is the model's ASU replicated by every
    space-group symmetry operator; a block of cells adds the lattice translations.

    The size of the block is set by exactly one of:
      * cells (n, m, o) -> an n x m x o block of unit cells along the a, b and c
        cell axes;
      * box  (x, y, z)  -> the smallest block of whole unit cells that fully
        fills an x by y by z Angstrom region (each count rounded up), so the
        result comes out slightly larger than the requested box.
    If neither is given, a single unit cell is generated.

    Runs headless as well as in the GUI. `refl_file`, `name`,
    `prune_special_positions` and `focus` behave as for `clipper symcopies`, and
    the result carries an invertible `clipper_sym_expansion` (see
    clipper.collapse_to_asu) - so a P1 box can be relaxed and folded back to the
    ASU for structure-factor work.
    '''
    from math import ceil
    from .symmetry import get_all_symmetry_handlers
    from .sym_realize import (crystal_symmetry_for, realize_symmetry_copies,
        unit_cell_places, pack_copies_into_block)

    if cells is not None and box is not None:
        raise UserError('Specify only one of "cells" or "box", not both.')
    if cells is not None and any(c < 1 for c in cells):
        raise UserError('"cells" must be three positive integers, e.g. cells 2,2,2')
    if box is not None and any(x <= 0 for x in box):
        raise UserError('"box" must be three positive distances (A), e.g. box 100,100,100')

    if models is not None:
        structures = list(models)
    else:
        structures = [sh.structure for sh in get_all_symmetry_handlers(session)]
    if not structures:
        raise UserError('No structures to generate unit cells for. Specify one or '
            'more models (and a "reflFile" if they carry no crystallographic '
            'metadata).')

    results = []
    for s in structures:
        sym = crystal_symmetry_for(s, refl_file=refl_file)
        if not sym.has_symmetry:
            session.logger.warning('Model {} has no crystallographic symmetry; '
                'skipping.'.format(s.id_string))
            continue
        if box is not None:
            na, nb, nc = (max(1, int(ceil(x/d))) for x, d in zip(box, sym.cell.dim))
        elif cells is not None:
            na, nb, nc = (int(c) for c in cells)
        else:
            na = nb = nc = 1
        places = unit_cell_places(sym.cell, sym.spacegroup, na, nb, nc)
        if len(places) <= 1:
            session.logger.warning('A single {} unit cell is just the input model; '
                'nothing to generate for {}.'.format(sym.spacegroup.symbol_hm,
                    s.id_string))
            continue
        cell_name = name
        if cell_name is None:
            cell_name = '{} unit cells ({}x{}x{})'.format(s.name, na, nb, nc)
        mults = _site_multiplicities_for(sym, s) if prune_special_positions else None
        combined = realize_symmetry_copies(session, s, places, name=cell_name,
            prune_special_positions=prune_special_positions, multiplicities=mults)
        if combined is None:
            session.logger.warning('No atoms survived generating unit cells for '
                '{}.'.format(s.id_string))
            continue
        # unit_cell_places places each copy at symop(x)+cell-offset without wrapping the
        # symop's own translation, so the raw block scatters across ~+/-1 cell. Pack each
        # copy's centroid into the block for a tidy periodic display (SF-invariant integer
        # lattice shift, applied after special-position dedup).
        pack_copies_into_block(combined, sym.cell, na, nb, nc)
        session.logger.info('Created model {} containing a {}x{}x{} block of unit '
            'cells for {}.'.format(combined.id_string, na, nb, nc, s.id_string))
        results.append(combined)

    if focus and results and session.ui.is_gui:
        from chimerax.core.commands import run
        run(session, 'view {}'.format(' '.join('#'+m.id_string for m in results)))
    return results


def register_clipper_cmd(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        BoolArg, FloatArg, IntArg, StringArg, EnumOf,
        Int3Arg, Float3Arg,
        OpenFileNameArg, SaveFileNameArg,
        ModelsArg,
        create_alias
        )
    from chimerax.atomic import StructuresArg, StructureArg, AtomsArg
    RadiationArg = EnumOf(['auto', 'xray', 'electron'])
    FragmentsArg = EnumOf(['off', 'rename', 'complete'])

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

    from chimerax.core.commands import ListOf
    open_desc = CmdDesc(
        required=[
            ('path', OpenFileNameArg),
        ],
        keyword=[
            ('structure_model', StructureArg),
            ('over_sampling', FloatArg),
            ('auto_free_flags', BoolArg),
            ('browse', BoolArg),
            ('fsigf', StringArg),
            ('map_columns', ListOf(StringArg)),
            ('free_flags', StringArg),
            ('free_flags_file', OpenFileNameArg),
            ('radiation', RadiationArg),
        ],
        synopsis='Open a structure factor .mtz or .cif file and generate maps for the given model '
                 '(radiation xray|electron for micro-ED / 3D-ED; default xray)'
    )
    register('clipper open', open_desc, open_structure_factors_and_add_to_session, logger=logger)

    smallmol_desc = CmdDesc(
        required=[
            ('path', OpenFileNameArg),
        ],
        keyword=[
            ('hkl', OpenFileNameArg),
            ('radiation', RadiationArg),
            ('fragments', FragmentsArg),
        ],
        synopsis='Open a small-molecule CIF as a live crystal structure: model in the '
                 'unit cell, symmetry, and (with reflections) live electron-density maps '
                 '(radiation xray|electron|auto for micro-ED support; fragments '
                 'off|rename|complete to split the ASU into named residues)'
    )
    register('clipper smallmol', smallmol_desc, open_small_molecule, logger=logger)

    cod_desc = CmdDesc(
        required=[
            ('cod_id', StringArg),
        ],
        keyword=[
            ('ignore_cache', BoolArg),
            ('radiation', RadiationArg),
            ('fragments', FragmentsArg),
        ],
        synopsis='Fetch a structure from the Crystallography Open Database (by numeric '
                 'COD ID) and open it as a live crystal structure'
    )
    register('clipper cod', cod_desc, fetch_cod_crystal, logger=logger)

    fragments_desc = CmdDesc(
        required=[
            ('structures', StructuresArg),
        ],
        keyword=[
            ('mode', FragmentsArg),
        ],
        synopsis='Split a small-molecule asymmetric unit into named covalent-fragment '
                 'residues (mode off|rename|complete; complete reassembles molecules via '
                 'symmetry - rewraps scattered hydrogens + adds atoms finishing molecules '
                 'split across a special position)'
    )
    register('clipper fragments', fragments_desc, split_fragments_cmd, logger=logger)

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

    debugmaps_desc = CmdDesc(
        optional=[
            ('models', StructuresArg),
            ('enable', BoolArg),
        ],
        synopsis='Add/remove live structure-factor component maps (Fcalc, Fatoms, '
                 'Fbulk, Fmask) for the given model(s), as a debugging aid'
    )
    register('clipper debugmaps', debugmaps_desc, debug_maps, logger=logger)

    scaling_desc = CmdDesc(
        optional=[
            ('models', StructuresArg),
        ],
        keyword=[
            ('all_reflections', BoolArg),
            ('reflections_per_bin', IntArg),
            ('num_bins', IntArg),
            ('seed', IntArg),
        ],
        synopsis='Control live-map bulk-solvent scaling: allReflections (exact, fit '
                 'over all reflections) vs a fixed seeded subset, and the subset '
                 'size/seed'
    )
    register('clipper scaling', scaling_desc, set_scaling, logger=logger)

    discard_sym_desc = CmdDesc(
        required=[
            ('models', AtomicStructuresOrSymmetryMgrsArg),
        ],
        synopsis='Discard a model\'s crystallographic symmetry (revert to a P1 box) '
                 'and rewrite its metadata so the change persists across sessions. '
                 'Fails if a crystal dataset is currently loaded.'
    )
    register('clipper discardSymmetry', discard_sym_desc, discard_symmetry, logger=logger)

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

    symcopies_desc = CmdDesc(
        optional=[('models', StructuresArg)],
        keyword=[
            ('contacting', AtomsArg),
            ('distance', FloatArg),
            ('name', StringArg),
            ('focus', BoolArg),
            ('refl_file', OpenFileNameArg),
            ('prune_special_positions', BoolArg),
        ],
        synopsis=('Make real, whole-model copies of crystallographic symmetry '
            'mates and combine them with the original ASU into a single new '
            'model. With "contacting", realises every symmetry copy approaching '
            'within "distance" of the given atom selection; otherwise realises '
            'the copies currently displayed in the spotlight (GUI), or - with no '
            'spotlight, e.g. headless - every copy within "distance" of the whole '
            'model. "reflFile" (a .mtz/.cif structure-factor file) supplies the '
            'cell and spacegroup for models lacking crystallographic metadata. '
            'By default, redundant special-position copies of small molecules '
            '(ions/water) are pruned; set pruneSpecialPositions false to keep them.')
    )
    register('clipper symcopies', symcopies_desc, symmetry_copies_real, logger=logger)

    unitcells_desc = CmdDesc(
        optional=[('models', StructuresArg)],
        keyword=[
            ('cells', Int3Arg),
            ('box', Float3Arg),
            ('name', StringArg),
            ('focus', BoolArg),
            ('refl_file', OpenFileNameArg),
            ('prune_special_positions', BoolArg),
        ],
        synopsis=('Generate whole crystallographic unit cells as a single new '
            'model. With "cells n,m,o", make an n x m x o block of unit cells; '
            'with "box x,y,z", make the smallest block of whole cells that fills '
            'an x by y by z Angstrom region; with neither, one unit cell. '
            '"reflFile" (a .mtz/.cif structure-factor file) supplies the cell and '
            'spacegroup for models lacking crystallographic metadata. By default, '
            'redundant special-position copies of small molecules (ions/water) '
            'are pruned; set pruneSpecialPositions false to keep them.')
    )
    register('clipper unitcells', unitcells_desc, generate_unit_cells, logger=logger)

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

    bsharp_desc = CmdDesc(
        required=[('b_sharp', FloatArg)],
        optional=[('maps', ModelsArg)],
        synopsis='Set the sharpening (positive) / smoothing (negative) B-factor '
            'applied to live crystallographic maps. If no maps are given, applies '
            'to all live maps in the session. Maps are recomputed immediately '
            'without moving atoms.'
    )
    register('clipper bsharp', bsharp_desc, set_map_bsharp, logger=logger)

    set_oversampling_desc = CmdDesc(
        required=[('oversampling', FloatArg)],
        optional=[('models', ModelsArg)],
        synopsis='Set the real-space map oversampling (Shannon) rate for live and '
            'static crystallographic maps at runtime. Higher values give a finer '
            'grid (smoother contours) at the cost of memory and FFT time; typical '
            'values are 1.5-3.0. Applies to all crystallographic map managers if '
            'none are specified.'
    )
    register('clipper set oversampling', set_oversampling_desc, set_oversampling, logger=logger)


def register_cview_cmd(logger):
    from chimerax.core.commands import create_alias
    create_alias('cview', 'view $*; cofr center showpivot true', logger=logger)
