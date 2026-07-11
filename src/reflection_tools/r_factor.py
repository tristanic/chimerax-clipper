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

# Single source of truth for the crystallographic residual R = sum w|Fo - Fc| /
# sum w|Fo|. Every surface that reports an R-factor - the small-molecule live map
# (SmallMoleculeXmapMgr), the clipper cod/smallmol command log, and the headless
# import metrics (_structure_factor_metrics) - routes through here so the numbers
# cannot drift apart. The optional epsilon weighting reproduces the definition the
# C++ macromolecular manager (Xtal_mgr_base::calculate_r_factors) uses, so a work/
# free split computed here matches that reference.

from collections import namedtuple

# r_*   : the R-factor over each reflection set (None where the set is empty/absent).
# n_*   : the number of reflections contributing to each.
RFactors = namedtuple('RFactors',
                      ('r_all', 'r_work', 'r_free', 'r_observed',
                       'n_all', 'n_work', 'n_free', 'n_observed'))


def compute_r_factors(fo, fc, *, free_flags=None, free_value=0,
                      epsilon=None, observed_mask=None):
    '''
    Compute the crystallographic R-factor R = sum w|Fo - Fc| / sum w|Fo| over the
    relevant reflection subsets, in one place.

    Args:
        fo, fc: aligned 1-D arrays of observed and *scaled* calculated amplitudes.
            Fc must already be on the Fo scale (this routine applies no scaling).
        free_flags: optional aligned integer array partitioning work/free
            reflections; entries equal to free_value are the free set. When None,
            r_work and r_free are None (no work/free split - the small-molecule case).
        free_value: the flag value marking a free reflection (CCP4 convention: 0).
        epsilon: optional aligned per-reflection multiplicity (from
            ih.hkl_class.epsilon). When given, each term is weighted w = 2/epsilon,
            reproducing the C++ manager's definition; when None, w = 1 (unweighted,
            the small-molecule convention).
        observed_mask: optional aligned boolean mask selecting the observed subset
            (e.g. Fo > 4 sigma, i.e. I > 2 sigma). When None, r_observed is None.

    Non-finite entries and Fo <= 0 are dropped before summing. Returns an RFactors.
    '''
    import numpy
    fo = numpy.asarray(fo, dtype=float)
    fc = numpy.asarray(fc, dtype=float)

    valid = numpy.isfinite(fo) & numpy.isfinite(fc) & (fo > 0)

    if epsilon is not None:
        eps = numpy.asarray(epsilon, dtype=float)
        # 2/epsilon mirrors Xtal_mgr_base; the factor 2 cancels in the ratio, so it
        # is the 1/epsilon relative down-weighting of high-multiplicity reflections
        # that matters. Guard against a zero/invalid epsilon.
        with numpy.errstate(divide='ignore', invalid='ignore'):
            w = numpy.where(eps > 0, 2.0 / eps, 0.0)
        valid = valid & numpy.isfinite(w)
    else:
        w = numpy.ones_like(fo)

    def _r(mask):
        mask = mask & valid
        n = int(mask.sum())
        if n == 0:
            return None, 0
        den = float((w[mask] * fo[mask]).sum())
        if den <= 0:
            return None, n
        num = float((w[mask] * numpy.abs(fo[mask] - fc[mask])).sum())
        return num / den, n

    all_mask = numpy.ones(len(fo), bool)
    r_all, n_all = _r(all_mask)

    if free_flags is not None:
        flags = numpy.asarray(free_flags)
        free_mask = (flags == free_value)
        r_work, n_work = _r(~free_mask)
        r_free, n_free = _r(free_mask)
    else:
        r_work = r_free = None
        n_work = n_free = 0

    if observed_mask is not None:
        r_observed, n_observed = _r(numpy.asarray(observed_mask, bool))
    else:
        r_observed, n_observed = None, 0

    return RFactors(r_all=r_all, r_work=r_work, r_free=r_free, r_observed=r_observed,
                    n_all=n_all, n_work=n_work, n_free=n_free, n_observed=n_observed)


def r_factor(fo, fc, *, subset='all', **kwargs):
    '''
    Convenience wrapper returning a single R-factor for the requested subset
    ('all', 'work', 'free' or 'observed'). Accepts the same keyword arguments as
    compute_r_factors. Returns None when that subset is empty or absent.
    '''
    result = compute_r_factors(fo, fc, **kwargs)
    try:
        return getattr(result, 'r_' + subset)
    except AttributeError:
        raise ValueError("subset must be one of 'all', 'work', 'free', 'observed'")
