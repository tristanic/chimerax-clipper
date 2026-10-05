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

'''
Opt-in variants of the small-molecule differentiable structure-factor target.

The default target (:func:`chimerax.clipper.io.small_molecule.fobs_from_arrays`) fits
amplitudes ``Fo = sqrt(I)`` weighted by counting statistics alone. :class:`SFTarget`
selects any combination of:

  * **observations** -- ``'amplitude'`` (the default convention); ``'intensity_masked'``
    (intensities rebuilt as ``Fo^2``, so every ``I <= 0`` reflection drops out);
    ``'french_wilson'`` (amplitudes from the French & Wilson posterior); or
    ``'intensity'`` (the measured ``I`` and ``sigma(I)`` themselves, negative ``I``
    included, fitted as ``1/2 sum w (Io - s^2|Fc|^2)^2``);
  * **weights** -- ``'sigma'`` (``w = 1/sigma^2``) or ``'shelxl'``, the depositor's own
    ``w = 1/[sigma^2(Fo^2) + (aP)^2 + bP]``, ``P = (max(Fo^2,0) + 2Fc^2)/3``, frozen at
    build from the depositor's ``Fc^2``;
  * **extinction** -- ``'none'`` or ``'shelxl'``, the depositor's refined SHELXL EXTI
    correction ``Fc* = kFc[1 + 0.001 x Fc^2 lambda^3/sin(2 theta)]^-1/4``. A SHELXL
    ``.fcf`` (``_shelx_refln_list_code`` 4 or 6 -- nearly every COD ``.hkl``) already
    carries it: its Fo^2 are the measured intensities divided by the refined ``y(h)``
    (and its Fc^2 are not multiplied by it), so on such data the correction is in the
    observations already and ``'shelxl'`` adds nothing (flagged
    ``extinction_in_data``). It is applied only to data without that mark.

Every variant is built by one pipeline on the raw ``(hkl, fsq, sig)`` arrays: merge
equivalent rows; apply the frozen corrections in intensity space; transform to the
chosen observation space. The weights become an effective ``sigma`` (``w`` is then
exactly the depositor's weight) and the extinction factor ``y(h)``, where it applies,
divides ``I`` and ``sigma`` (at a frozen ``y`` that is identical, in value and gradient,
to applying it to the model -- and is what SHELXL itself writes to an ``.fcf``), so both
corrections compose with every observation space and reach the C++ evaluator as
ordinary observations.

The per-crystal inputs the corrections need (the depositor's ``Fc^2``, ``a``, ``b``,
``x``, the wavelength, and Clipper's own ``Fc`` at the deposited model) are gathered once
from the CIF by :func:`small_molecule_target_aux`, row-aligned with the raw arrays, so a
caller that caches the raw arrays can cache these beside them and rebuild the target with
no file access.
'''

import re

SF_TARGET_VERSION = 1

_OBSERVATIONS = {'amplitude': 'a', 'intensity_masked': 'im', 'french_wilson': 'fw',
                 'intensity': 'i'}
_WEIGHTS = ('sigma', 'shelxl')
_EXTINCTION = ('none', 'shelxl')

_NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'
# SHELXL's P = (max(Fo^2,0) + 2Fc^2)/3, in either order, once carets are dropped.
_FO2 = r'(?:max\(fo?2,0\)|fo2)'
_SHELXL_P = r'p=[(\[](?:%s\+2fc2|2fc2\+%s)[)\]]/3' % (_FO2, _FO2)


class SFTarget:
    '''
    Which small-molecule structure-factor target to build (see the module docstring).
    ``SFTarget()`` is the default target; any other choice is a different target, so a
    caller with a signed corpus must carry :meth:`token` in its signature. Every
    non-default target is defined on merged reflections (``merge_equivalents=True``).
    '''
    def __init__(self, observations='amplitude', weights='sigma', extinction='none'):
        if observations not in _OBSERVATIONS:
            raise ValueError('SFTarget: observations must be one of %s, not %r'
                             % (sorted(_OBSERVATIONS), observations))
        if weights not in _WEIGHTS:
            raise ValueError('SFTarget: weights must be one of %s, not %r'
                             % (list(_WEIGHTS), weights))
        if extinction not in _EXTINCTION:
            raise ValueError('SFTarget: extinction must be one of %s, not %r'
                             % (list(_EXTINCTION), extinction))
        self.observations = observations
        self.weights = weights
        self.extinction = extinction

    @property
    def is_default(self):
        return (self.observations == 'amplitude' and self.weights == 'sigma'
                and self.extinction == 'none')

    @property
    def kind(self):
        '''The evaluator's target form: ``'intensity'`` or ``'amplitude'``.'''
        return ('intensity' if self.observations in ('intensity', 'intensity_masked')
                else 'amplitude')

    @property
    def needs_aux(self):
        return self.weights != 'sigma' or self.extinction != 'none'

    def token(self):
        '''A signature token for this target: ``''`` for the default (so signatures that
        predate it are unchanged), else e.g. ``'|sft1-i-wshelxl-xshelxl'``. The leading
        number is :data:`SF_TARGET_VERSION`, bumped whenever a variant's definition
        changes. Equivalent-row merging is not part of the token (every non-default
        target merges).'''
        if self.is_default:
            return ''
        return '|sft%d-%s-w%s-x%s' % (SF_TARGET_VERSION, _OBSERVATIONS[self.observations],
                                       self.weights, self.extinction)

    def as_dict(self):
        return {'observations': self.observations, 'weights': self.weights,
                'extinction': self.extinction}

    def __eq__(self, other):
        return isinstance(other, SFTarget) and self.as_dict() == other.as_dict()

    def __hash__(self):
        return hash(tuple(sorted(self.as_dict().items())))

    def __repr__(self):
        return 'SFTarget(%r, %r, %r)' % (self.observations, self.weights, self.extinction)


def is_default_target(sf_target):
    return sf_target is None or sf_target.is_default


def check_target_args(sf_target, kind, merge_equivalents):
    '''Validate an entry point's ``sf_target`` against its legacy ``kind`` /
    ``merge_equivalents`` arguments; returns the ``kind`` to build.'''
    if is_default_target(sf_target):
        return kind
    if kind != 'amplitude':
        raise ValueError('pass the target form through sf_target, not kind=%r' % kind)
    if not merge_equivalents:
        raise ValueError('%r is defined on merged reflections: pass '
                         'merge_equivalents=True' % (sf_target,))
    return sf_target.kind


# ---------------------------------------------------------------------------------------
# The depositor's refinement record
# ---------------------------------------------------------------------------------------

def cif_number(value):
    '''The value of a CIF number such as ``0.097(8)`` or ``'x = 0.0066(3)'``, as a float;
    None for an absent or non-numeric value.'''
    if value is None:
        return None
    s = str(value).strip().strip('\'"')
    if s in ('', '?', '.'):
        return None
    m = re.search(_NUM, s)
    if not m:
        return None
    return float(m.group(0))


def _top_level_terms(s):
    out, depth, cur = [], 0, ''
    for ch in s:
        if ch in '([':
            depth += 1
        elif ch in ')]':
            depth -= 1
        if ch == '+' and depth == 0:
            out.append(cur)
            cur = ''
        else:
            cur += ch
    out.append(cur)
    return [t for t in out if t]


def parse_shelxl_weighting(text):
    '''
    The weighting scheme in ``_refine_ls_weighting_details`` (or an old CIF's
    ``_refine_ls_weighting_scheme`` text). Returns ``{'form', 'a', 'b', 'p_max0'}``:

      * ``'shelxl'`` -- ``w = 1/[sigma^2(Fo^2) + (aP)^2 + bP]`` (SHELXL's WGHT, and
        CRYSTALS' "modified Sheldrick"); an absent a or b term is 0;
      * ``'sigma_only'`` -- ``w = 1/sigma^2(Fo^2)``, a = b = 0;
      * ``'other'`` -- anything else (refinement on F, Chebyshev, extra terms, or a P
        that is not SHELXL's ``(Fo^2 + 2Fc^2)/3``), a = b = None;
      * ``'missing'`` -- no text.

    ``p_max0`` is True when the text writes P with ``max(Fo^2, 0)``, False when it writes
    ``(Fo^2 + 2Fc^2)/3``, None when it does not define P. SHELXL itself always uses
    ``max(Fo^2, 0)``, whichever the CIF prints.
    '''
    if text is None or not str(text).strip() or str(text).strip() in ('?', '.'):
        return {'form': 'missing', 'a': None, 'b': None, 'p_max0': None}
    s = re.sub(r'[\s*~]', '', str(text))
    p_max0 = (True if re.search(r'max\(F', s, re.I)
              else (False if re.search(r'P=', s) else None))
    other = {'form': 'other', 'a': None, 'b': None, 'p_max0': p_max0}
    if p_max0 is not None and not re.search(_SHELXL_P, s.replace('^', '').lower()):
        return other
    m = re.search(r'w=1/\[([^\]]*)\]', s)
    if not m:
        m2 = re.search(r'w=1/\\s\^2\^\(F(o)?\^2\^\)(?![+\w(])', s)
        return ({'form': 'sigma_only', 'a': 0.0, 'b': 0.0, 'p_max0': p_max0}
                if m2 else other)
    terms = _top_level_terms(m.group(1))
    if not terms or not re.fullmatch(r'\\s\^2\^\(F(o)?\^2\^\)', terms[0]):
        return other
    a = b = 0.0
    for t in terms[1:]:
        ma = re.fullmatch(r'\((' + _NUM + r')P\)\^?2\^?', t)
        mb = re.fullmatch(r'(' + _NUM + r')P', t)
        if ma:
            a = float(ma.group(1))
        elif mb:
            b = float(mb.group(1))
        else:
            return other
    form = 'sigma_only' if len(terms) == 1 else 'shelxl'
    return {'form': form, 'a': a, 'b': b, 'p_max0': p_max0}


def parse_extinction(method, expression, coef):
    '''
    The depositor's extinction model. Returns ``{'kind', 'x'}``: ``kind`` is ``'none'``
    (no refined coefficient, or a zero one), ``'shelxl'`` (a positive coefficient under
    SHELX's ``Fc* = kFc[1 + 0.001 x Fc^2 lambda^3/sin(2 theta)]^-1/4``, recognised by
    the program name or the expression) or ``'other'`` (a positive coefficient under any
    other model: Larson, Zachariasen, Becker-Coppens, ...). ``x`` is the coefficient.
    '''
    x = cif_number(coef)
    if x is None or x <= 0.0:
        return {'kind': 'none', 'x': x}
    txt = ' '.join(str(v) for v in (method, expression) if v).lower()
    compact = re.sub(r'[\s~*]', '', txt)
    if 'shelx' in txt or '0.001xfc^2^' in compact or '0.001xfc2' in compact:
        return {'kind': 'shelxl', 'x': x}
    return {'kind': 'other', 'x': x}


def read_small_molecule_refinement(cif_path, hkl_path=None):
    '''
    The depositor's refinement record for a small-molecule crystal, with the reflection
    rows read exactly as :func:`chimerax.clipper.io.small_molecule._parse_reflection_file`
    reads them (``hkl_path`` or, by default, the ``.hkl`` beside ``cif_path``), so its
    per-row arrays line up with a cache of that function's ``(hkl, fsq, sig)``.

    Returns a dict:
      * ``'hkl'``: the rows' Miller indices (to check the alignment);
      * ``'fc_sq_dep'``: the depositor's ``_refln_F_squared_calc`` per row, NaN where the
        file gives none;
      * ``'weighting'``: :func:`parse_shelxl_weighting` of the CIF's weighting text, plus
        that ``'text'``;
      * ``'extinction'``: :func:`parse_extinction` plus the ``'method'`` / ``'expression'``;
      * ``'wavelength'`` (A, or None), ``'f_squared_multiplier'`` (the reflection file's
        ``_shelx_f_squared_multiplier``, default 1), ``'shelx_refln_list_code'`` (4 or 6
        for a SHELXL ``.fcf``, whose Fo^2 already carry SHELXL's extinction correction;
        None otherwise), ``'structure_factor_coef'``;
      * ``'goof'``, ``'wr2'``, ``'n_parameters'``: the depositor's reported fit, for
        checking the weights.
    '''
    import numpy
    from chimerax.core.errors import UserError
    from chimerax.mmcif import get_cif_tables
    from ..symmetry import _first_cif_field
    from .small_molecule import _find_reflection_source, _reflection_columns

    sf_path = hkl_path or cif_path
    refl_file, refln = _find_reflection_source(sf_path)
    if refln is None:
        raise UserError('No reflections (F_squared_meas) found for %r' % sf_path)
    hkl, _, _, fc_sq = _reflection_columns(refln, extra=('F_squared_calc',))

    multiplier = 1.0
    shelx = get_cif_tables(refl_file, ['shelx'])
    shelx = shelx[0] if shelx else None
    v = cif_number(_first_cif_field(shelx, 'F_squared_multiplier'))
    if v is not None and v > 0:
        multiplier = v
    list_code = cif_number(_first_cif_field(shelx, 'refln_list_code'))

    ref, dif, drad = (get_cif_tables(cif_path, ['refine', 'diffrn', 'diffrn_radiation'])
                      + [None] * 3)[:3]
    text = _first_cif_field(ref, 'ls_weighting_details')
    if text is None:
        scheme = _first_cif_field(ref, 'ls_weighting_scheme')
        if scheme is not None and '1/' in scheme:     # old CIFs: the formula is here
            text = scheme
    weighting = parse_shelxl_weighting(text)
    weighting['text'] = text
    method = _first_cif_field(ref, 'ls_extinction_method')
    expression = _first_cif_field(ref, 'ls_extinction_expression')
    extinction = parse_extinction(method, expression,
                                  _first_cif_field(ref, 'ls_extinction_coef'))
    extinction['method'] = method
    extinction['expression'] = expression
    wavelength = cif_number(_first_cif_field(drad, 'wavelength'))
    if wavelength is None:
        wavelength = cif_number(_first_cif_field(dif, 'radiation_wavelength'))
    if wavelength is not None and wavelength <= 0:
        wavelength = None
    n_par = cif_number(_first_cif_field(ref, 'ls_number_parameters'))
    return {
        'hkl': hkl,
        'fc_sq_dep': fc_sq,
        'weighting': weighting,
        'extinction': extinction,
        'wavelength': wavelength,
        'f_squared_multiplier': multiplier,
        'shelx_refln_list_code': None if list_code is None else int(list_code),
        'structure_factor_coef': _first_cif_field(ref, 'ls_structure_factor_coef'),
        'goof': cif_number(_first_cif_field(ref, 'ls_goodness_of_fit_ref')),
        'wr2': cif_number(_first_cif_field(ref, 'ls_wR_factor_ref')),
        'n_parameters': None if n_par is None else int(n_par),
    }


def _lookup_by_orbit(hkl, ref_hkl, ref_values, spacegroup):
    '''Per row of ``hkl``, the value of ``ref_values`` at whichever ``ref_hkl`` row is
    symmetry- or Friedel-equivalent to it (NaN where none is).'''
    import numpy
    from .small_molecule import _orbit_keys
    nan = numpy.full(len(hkl), numpy.nan)
    if not len(hkl) or not len(ref_hkl):
        return nan
    ref_keys = _orbit_keys(ref_hkl, spacegroup)[1].max(axis=1)
    keys = _orbit_keys(hkl, spacegroup)[1].max(axis=1)
    order = numpy.argsort(ref_keys)
    pos = numpy.clip(numpy.searchsorted(ref_keys[order], keys), 0, len(order) - 1)
    hit = ref_keys[order][pos] == keys
    out = nan
    out[hit] = numpy.asarray(ref_values, numpy.double)[order][pos[hit]]
    return out


def deposited_fcalc_sq(cif_path, hkl, fsq, sig, cell, spacegroup, grid, radiation='xray'):
    '''
    Clipper's structure factors at the deposited model, per row of the raw reflection
    arrays: ``(fc_sq_abs, fc_sq_scaled)``. ``fc_sq_abs`` is ``|Fc|^2`` on the absolute
    scale (electrons per cell, as SHELXL's extinction expression uses it);
    ``fc_sq_scaled`` is ``(s(h)|Fc|)^2`` on the observed scale, with ``s(h)`` the same
    anisotropic-Gaussian x isotropic-spline scale ``recomputed_r_factor`` fits to the
    merged data. Exact direct summation over the same CIF-derived atom list
    ``recomputed_r_factor`` uses. NaN where a row has no reflection (e.g. a systematic
    absence) and everywhere if the model cannot be read.
    '''
    import numpy
    from .. import HKL_info
    from ..clipper_python import SFcalc_aniso_sum_float
    from ..clipper_python.data32 import HKL_data_F_phi_float, HKL_data_F_sigF_float
    from ..clipper_python.ext import scale_fcalc_to_fobs
    from .small_molecule import (merge_equivalent_reflections, _sfcalc_atom_list,
        _d_min, _padded_resolution, _amplitudes_from_intensities, _radiation_enum)

    hkl = numpy.ascontiguousarray(hkl, numpy.int32).reshape(-1, 3)
    nan = numpy.full(len(hkl), numpy.nan)
    mh, mi, ms, _ = merge_equivalent_reflections(hkl, fsq, sig, spacegroup)
    atoms = _sfcalc_atom_list(cif_path, cell, spacegroup, grid, radiation)
    if atoms is None or not len(mh):
        return nan, nan.copy()
    hkl_info = HKL_info(spacegroup, cell, _padded_resolution(_d_min(mh, cell)), True)
    fo, sigf = _amplitudes_from_intensities(mi, ms)
    fobs = HKL_data_F_sigF_float(hkl_info)
    fobs.set_data(mh, numpy.stack([fo, sigf], axis=1).astype(numpy.float32))
    fcalc = HKL_data_F_phi_float(hkl_info)
    SFcalc_aniso_sum_float(radiation=_radiation_enum(radiation))(fcalc, atoms)
    scaled = HKL_data_F_phi_float(hkl_info)
    scale_fcalc_to_fobs(fcalc, fobs, scaled)
    fh, fv = fcalc.data
    sh, sv = scaled.data
    fc_abs = _lookup_by_orbit(hkl, fh, numpy.asarray(fv, numpy.double)[:, 0] ** 2, spacegroup)
    fc_scl = _lookup_by_orbit(hkl, sh, numpy.asarray(sv, numpy.double)[:, 0] ** 2, spacegroup)
    return fc_abs, fc_scl


def small_molecule_target_aux(cif_path, hkl_path=None, *, radiation='auto'):
    '''
    Everything the non-default targets need beyond the raw ``(hkl, fsq, sig)``, gathered
    once from the files and row-aligned with
    ``small_molecule._parse_reflection_file(hkl_path or cif_path)``: the
    :func:`read_small_molecule_refinement` record plus ``'fc_sq_model_abs'`` and
    ``'fc_sq_model_scaled'`` (:func:`deposited_fcalc_sq`). Plain numpy arrays, floats,
    strings and dicts, so it can be stored beside a cache of the raw arrays and passed
    back as ``aux`` with no file access.
    '''
    import numpy
    from chimerax.core.errors import UserError
    from ..symmetry import crystal_symmetry_from_cif_file
    from .small_molecule import _parse_reflection_file, _resolve_radiation

    radiation = _resolve_radiation(radiation, cif_path)
    cell, spacegroup, grid = crystal_symmetry_from_cif_file(cif_path)
    aux = read_small_molecule_refinement(cif_path, hkl_path)
    hkl, fsq, sig = _parse_reflection_file(hkl_path or cif_path)
    if not numpy.array_equal(hkl, aux['hkl']):
        raise UserError('small_molecule_target_aux: reflection rows of %r do not line up'
                        % (hkl_path or cif_path,))
    aux['fc_sq_model_abs'], aux['fc_sq_model_scaled'] = deposited_fcalc_sq(
        cif_path, hkl, fsq, sig, cell, spacegroup, grid, radiation)
    aux['radiation'] = radiation
    return aux


# ---------------------------------------------------------------------------------------
# The pipeline
# ---------------------------------------------------------------------------------------

def _copies(fc_sq, fsq):
    '''True when a deposited ``F_squared_calc`` column merely repeats ``F_squared_meas``
    (some reflection files carry Fo^2 in both), so it holds no calculated values.'''
    import numpy
    both = numpy.isfinite(fc_sq) & numpy.isfinite(fsq)
    n = int(both.sum())
    return n > 0 and int((fc_sq[both] == fsq[both]).sum()) >= 0.95 * n


def _invresolsq(hkl, cell):
    import numpy
    from ..clipper_python import HKL
    return numpy.array([HKL(h.tolist()).invresolsq(cell) for h in hkl], numpy.double)


def shelxl_extinction_factor(hkl, cell, x, wavelength, fc_sq_abs):
    '''
    SHELXL's EXTI factor on intensities, ``y(h) = [1 + 0.001 x |Fc|^2 lambda^3 /
    sin(2 theta)]^-1/2`` (``Fc* = kFc y^1/2``), with ``|Fc|^2`` on the absolute scale
    (``fc_sq_abs``, e.g. :func:`deposited_fcalc_sq`'s first array) and ``sin(theta) =
    lambda/(2d)``. NaN where ``fc_sq_abs`` is.
    '''
    import numpy
    sin_t = 0.5 * wavelength * numpy.sqrt(_invresolsq(numpy.asarray(hkl).reshape(-1, 3), cell))
    sin_2t = 2.0 * sin_t * numpy.sqrt(numpy.clip(1.0 - sin_t ** 2, 0.0, None))
    with numpy.errstate(divide='ignore', invalid='ignore'):
        return (1.0 + 0.001 * x * numpy.asarray(fc_sq_abs, numpy.double)
                * wavelength ** 3 / sin_2t) ** -0.5


def target_observation_arrays(hkl, fsq, sig, cell, spacegroup, sf_target, aux=None):
    '''
    The observations a target installs, as plain numpy arrays (no Clipper HKL_data):
    steps 1-2 of the pipeline (merge; frozen weights and extinction) and the amplitude
    conversion. Use it to describe a target exactly (e.g. its weights' sum) without
    building one. ``aux`` is :func:`small_molecule_target_aux`'s record for the same raw
    rows; it is required when ``sf_target.needs_aux``.

    Returns a dict over the merged reflections:
      * ``'hkl'``, ``'multiplicity'``;
      * ``'i'``, ``'sig'``: the intensity and effective sigma the target fits (after the
        weights and the extinction factor), so ``'w' = 1/sig^2`` is the weight of an
        intensity-space residual;
      * ``'y'``: the extinction factor (1 where none applies);
      * ``'fo'``, ``'sigf'``: the same data in the default amplitude convention
        (``Fo = sqrt(I)``, ``sigF = sig/(2 sqrt(I))``; ``Fo = 0``, ``sigF = 1`` for
        ``I <= 0``);
      * ``'sum_w_obs2'``: ``sum w Io^2`` (intensity targets) or ``sum w Fo^2`` (amplitude
        targets; French-Wilson excluded), the denominator of ``wR^2``;
      * ``'flags'``: what was applied -- ``weights_applied``, ``weight_form``,
        ``p_source`` (``'depositor'`` / ``'model'`` / ``'mixed'``: where P's ``Fc^2``
        came from), ``n_weight_fallback`` (reflections left at sigma-only weight because
        no ``Fc^2`` was available), ``fc_sq_dep_copies_fo`` (the deposited ``Fc^2``
        column only repeats ``Fo^2``, so the model's is used), ``f_squared_multiplier``,
        ``extinction_applied``, ``extinction_in_data`` (SHELXL .fcf data, already
        corrected, so nothing further applied), ``extinction_kind``, and a ``reasons``
        list for every correction that was asked for but could not be (fully) applied.
    '''
    import numpy
    from chimerax.core.errors import UserError
    from .small_molecule import _equivalence_groups, _amplitudes_from_intensities

    if sf_target is None:
        sf_target = SFTarget()
    if sf_target.needs_aux and aux is None:
        raise ValueError('%r needs the aux record (small_molecule_target_aux)' % (sf_target,))
    hkl = numpy.ascontiguousarray(hkl, numpy.int32).reshape(-1, 3)
    fsq = numpy.ascontiguousarray(fsq, numpy.double).reshape(-1)
    sig = numpy.ascontiguousarray(sig, numpy.double).reshape(-1)
    n_raw = len(hkl)
    if n_raw == 0 or len(fsq) != n_raw or len(sig) != n_raw:
        raise UserError('target_observation_arrays: hkl/fsq/sig length mismatch or empty '
                        '(%d/%d/%d)' % (n_raw, len(fsq), len(sig)))
    if aux is not None:
        for key in ('fc_sq_dep', 'fc_sq_model_abs', 'fc_sq_model_scaled'):
            if key in aux and len(aux[key]) != n_raw:
                raise UserError('target_observation_arrays: aux %r has %d rows, the '
                                'reflections %d' % (key, len(aux[key]), n_raw))

    g = _equivalence_groups(hkl, fsq, sig, spacegroup)
    i_m, sig_m, mult = g.merge_inverse_variance(fsq)
    rep = g.rep[g.present]
    if len(rep) == 0:
        raise UserError('target_observation_arrays: no usable reflections')

    flags = {'weights_applied': False, 'weight_form': None, 'p_source': None,
             'n_weight_fallback': 0, 'fc_sq_dep_copies_fo': False,
             'f_squared_multiplier': None, 'extinction_applied': False,
             'extinction_in_data': False, 'extinction_kind': None, 'reasons': []}
    sig_used = sig_m.copy()
    if sf_target.weights == 'shelxl':
        wt = aux['weighting']
        flags['weight_form'] = wt['form']
        if wt['form'] in ('shelxl', 'sigma_only'):
            fc_dep = numpy.asarray(aux['fc_sq_dep'], numpy.double)
            if _copies(fc_dep, fsq):
                flags['fc_sq_dep_copies_fo'] = True
                flags['reasons'].append('the depositor\'s Fc^2 column is a copy of Fo^2: '
                                        'P uses the model Fc^2')
                fc_dep = numpy.full(n_raw, numpy.nan)
            fc2 = g.merge_mean(fc_dep)
            dep_ok = numpy.isfinite(fc2)
            fc2[~dep_ok] = g.merge_mean(aux['fc_sq_model_scaled'])[~dep_ok]
            ok = numpy.isfinite(fc2)
            n_dep, n_mod = int(dep_ok.sum()), int((ok & ~dep_ok).sum())
            flags['p_source'] = ('depositor' if n_mod == 0 else
                                 ('model' if n_dep == 0 else 'mixed'))
            # The reflection file's values are m x SHELXL's own (m rescales them to fit
            # the format), and only the b term of the weight is not scale-invariant.
            m = float(aux.get('f_squared_multiplier', 1.0))
            flags['f_squared_multiplier'] = m
            p = (numpy.maximum(i_m[ok], 0.0) + 2.0 * fc2[ok]) / 3.0
            sig_used[ok] = numpy.sqrt(sig_m[ok] ** 2 + (wt['a'] * p) ** 2 + wt['b'] * m * p)
            flags['n_weight_fallback'] = int((~ok).sum())
            flags['weights_applied'] = True
            if flags['n_weight_fallback']:
                flags['reasons'].append('%d reflection(s) without Fc^2 kept the sigma-only '
                                        'weight' % flags['n_weight_fallback'])
        else:
            flags['reasons'].append('weighting form %r is not SHELXL\'s: sigma-only weights'
                                    % wt['form'])

    y = numpy.ones(len(rep), numpy.double)
    if sf_target.extinction == 'shelxl':
        ext = aux['extinction']
        flags['extinction_kind'] = ext['kind']
        wl = aux.get('wavelength')
        if ext['kind'] == 'shelxl' and aux.get('shelx_refln_list_code') in (4, 6):
            # A SHELXL .fcf's Fo^2 are already divided by the refined y(h) (checked
            # against the raw HKLF data such files embed), so at a frozen y the target
            # on them is already the depositor's extinction-corrected one.
            flags['extinction_in_data'] = True
        elif ext['kind'] == 'shelxl' and wl:
            y = shelxl_extinction_factor(rep, cell, ext['x'], wl,
                                         g.merge_mean(aux['fc_sq_model_abs']))
            bad = ~(numpy.isfinite(y) & (y > 0))
            y[bad] = 1.0
            flags['extinction_applied'] = True
            if bad.any():
                flags['reasons'].append('%d reflection(s) without Fc kept y = 1'
                                        % int(bad.sum()))
        elif ext['kind'] == 'shelxl':
            flags['reasons'].append('no wavelength: no extinction correction')
        elif ext['kind'] == 'other':
            flags['reasons'].append('extinction model is not SHELXL\'s (%r): no correction'
                                    % (ext.get('method'),))

    i_used = i_m / y
    sig_used = sig_used / y
    fo, sigf = _amplitudes_from_intensities(i_used, sig_used)
    w = 1.0 / sig_used ** 2
    if sf_target.kind == 'intensity':
        use = (fo > 0) if sf_target.observations == 'intensity_masked' else slice(None)
        sum_w_obs2 = float(numpy.sum((w * i_used ** 2)[use]))
    elif sf_target.observations == 'french_wilson':
        sum_w_obs2 = float('nan')
    else:
        sum_w_obs2 = float(numpy.sum((fo / sigf) ** 2))
    return {'hkl': rep, 'multiplicity': mult, 'i': i_used, 'sig': sig_used, 'y': y,
            'w': w, 'fo': fo, 'sigf': sigf, 'sum_w_obs2': sum_w_obs2, 'flags': flags}


class SFObservations:
    '''
    A built target's observed data (:func:`observations_from_arrays`): ``fobs``,
    ``phi_fom``, ``usage``, ``iobs`` (``HKL_data_I_sigI`` for the ``'intensity'``
    target, else None), the evaluator ``kind``, the pipeline's ``arrays`` and its
    ``flags``. Clipper ``HKL_data`` holds a *non-owning* pointer to ``hkl_info``: keep
    this object (or ``hkl_info``) alive for as long as a target built from it is used.
    '''
    def __init__(self, hkl_info, fobs, iobs, phi_fom, usage, kind, arrays):
        self.hkl_info = hkl_info
        self.fobs = fobs
        self.iobs = iobs
        self.phi_fom = phi_fom
        self.usage = usage
        self.kind = kind
        self.arrays = arrays
        self.flags = arrays['flags']


def observations_from_arrays(hkl, fsq, sig, cell, spacegroup, sf_target, aux=None):
    '''
    Build the observed data for ``sf_target`` from raw reflection arrays (the cached
    ``(hkl, fsq, sig)`` and, when the target needs it, the matching
    :func:`small_molecule_target_aux` record). Always on merged reflections, with the
    resolution limit just beyond the highest-resolution one. Returns an
    :class:`SFObservations`; pass its ``fobs`` / ``phi_fom`` / ``usage`` / ``kind`` /
    ``iobs`` to :func:`chimerax.clipper.diff.crystal.ensemble_target_from_box`.

    French-Wilson needs enough reflections for its resolution-binned Wilson prior (it
    raises for fewer than 120).
    '''
    import numpy
    from .. import HKL_info, HKL_data_F_sigF, HKL_data_Phi_fom, HKL_data_Flag, HKL_data_I_sigI
    from .small_molecule import _d_min, _padded_resolution

    if sf_target is None:
        sf_target = SFTarget()
    arrays = target_observation_arrays(hkl, fsq, sig, cell, spacegroup, sf_target, aux)
    rep = arrays['hkl']
    n = len(rep)
    hkl_info = HKL_info(spacegroup, cell, _padded_resolution(_d_min(rep, cell)), True)

    def intensities():
        isigi = HKL_data_I_sigI(hkl_info)
        isigi.set_data(rep, numpy.stack([arrays['i'], arrays['sig']], axis=1)
                       .astype(numpy.float32))
        return isigi

    iobs = None
    if sf_target.observations == 'french_wilson':
        from ..reflection_tools import french_wilson_analytical
        fobs = french_wilson_analytical(intensities())
    else:
        fobs = HKL_data_F_sigF(hkl_info)
        fobs.set_data(rep, numpy.stack([arrays['fo'], arrays['sigf']], axis=1)
                      .astype(numpy.float32))
        if sf_target.observations == 'intensity':
            iobs = intensities()
    phi_fom = HKL_data_Phi_fom(hkl_info)
    phi_fom.set_data(rep, numpy.stack(
        [numpy.zeros(n, numpy.float32), numpy.ones(n, numpy.float32)], axis=1))
    usage = HKL_data_Flag(hkl_info)
    usage.set_data(rep, numpy.ones((n, 1), numpy.int32))
    return SFObservations(hkl_info, fobs, iobs, phi_fom, usage, sf_target.kind, arrays)
