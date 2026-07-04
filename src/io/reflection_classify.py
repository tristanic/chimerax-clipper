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

'''
Shared "column semantics" layer for reflection data.

Every crystallographic package (Phenix, CCP4/Refmac, Buster, SHELX, MR programs)
labels its MTZ/CIF columns differently, so deciding which amplitude pairs with
which sigma/phase, which column is the free-R set, and what a given map
coefficient actually *is* (2mFo-DFc vs mFo-DFc vs Fcalc) has historically been
done by scattered ad-hoc name checks. This module centralises that knowledge into
one place so the automatic loader and the (optional) GUI browser agree by
construction: both call the same recognisers and act on the same ranked result.

The label tables below are the exact, case-sensitive strings the major packages
write. Sources: CCP4 MTZ format + Refmac5 docs; Phenix refine/maps docs; Buster
MTZcolumns; gemmi hkl docs (rfree_column label set). Note the free-R *value*
convention (which integer means "in the test set": CCP4=0, Phenix=1, SHELX=-1) is
NOT handled here — it is already solved in C++ (Xtal_mgr_base::guess_free_flag_value).
This layer only decides *which column* is which.
'''

# --- map-coefficient classification -----------------------------------------

# Human-facing map-type constants.
MAP_2FOFC = '2mFo-DFc'
MAP_FOFC = 'mFo-DFc'
MAP_FCALC = 'Fcalc'
MAP_ANOM = 'anomalous difference'
MAP_UNKNOWN = 'map coefficients'


class MapCoeffInfo:
    '''Recognised semantics of one F/phi (map-coefficient) dataset.'''
    def __init__(self, name, columns, map_type, is_difference, human_label,
            is_fcalc):
        self.name = name
        self.columns = columns
        self.map_type = map_type
        self.is_difference = is_difference
        self.human_label = human_label
        self.is_fcalc = is_fcalc

    def __repr__(self):
        return '<MapCoeffInfo {!r} type={} diff={}>'.format(
            self.name, self.map_type, self.is_difference)


def _strip_producer_prefix(label):
    '''Drop the leading 'pdbx_' seen on structure-factor CIF column names.'''
    if label.startswith('pdbx_'):
        return label[len('pdbx_'):]
    return label


def classify_map_coefficient(amplitude_label):
    '''
    Classify a map-coefficient dataset from the *amplitude* column label.

    Returns a tuple (map_type, is_difference, human_label). The rules follow the
    de-facto conventions shared across packages:
      - 'FWT'  or a label beginning '2FOFC'         -> 2mFo-DFc  (best map)
      - 'DELFWT', or beginning 'FOFC' (no leading 2),
        or containing 'DEL'                          -> mFo-DFc   (difference)
      - Fcalc-family names (FC, FCALC, FMODEL, ...)  -> Fcalc
      - 'ANOM'                                        -> anomalous difference
    Matching is done case-insensitively after stripping any 'pdbx_' prefix, so
    both MTZ ('FWT', 'DELFWT', '2FOFCWT') and CIF ('pdbx_FWT', 'pdbx_DELFWT')
    forms are recognised.
    '''
    raw = amplitude_label.strip()
    core = _strip_producer_prefix(raw)
    up = core.upper()

    # Fcalc family first: total-model (incl. bulk solvent) and atoms-only forms.
    # These are the maps we want hidden-by-default (rarely what the user wants to
    # look at directly).
    fcalc_exact = {
        'FC', 'FCALC', 'FC_ALL', 'FC_ALL_LS', 'FMODEL',
        'F_CALC', 'F_CALC_AU', 'F-MODEL',
        'F_CALC_WITH_SOLVENT', 'F_CALC_PART_SOLVENT',
    }
    if up in fcalc_exact or up.startswith('FMODEL') or up.startswith('FC_ALL') \
            or up.startswith('F_CALC') or up.startswith('F-MODEL'):
        return (MAP_FCALC, False, 'Fcalc ({})'.format(raw))

    # 2mFo-DFc "best" map.
    if up == 'FWT' or up.startswith('2FOFC'):
        return (MAP_2FOFC, False, '2mFo-DFc ({})'.format(raw))

    # mFo-DFc difference map. 'FOFC...' without a leading '2', or any 'DEL' form.
    if up.startswith('DELFWT') or up.startswith('FOFC') or 'DEL' in up:
        return (MAP_FOFC, True, 'mFo-DFc ({})'.format(raw))

    # Anomalous difference map.
    if 'ANOM' in up:
        return (MAP_ANOM, True, 'anomalous difference ({})'.format(raw))

    return (MAP_UNKNOWN, None, raw)


def map_coeff_info(name):
    '''
    Build a MapCoeffInfo from a dataset name (e.g. '(dataset) FWT, PHWT' or
    'pdbx_FWT, pdbx_PHWT'). The amplitude label is the first column.
    '''
    columns = extract_column_labels(name)
    amp = columns[0] if columns else name
    map_type, is_difference, human = classify_map_coefficient(amp)
    return MapCoeffInfo(name, columns, map_type, is_difference, human,
        is_fcalc=(map_type == MAP_FCALC))


# --- free-R flag classification ----------------------------------------------

# Canonical free-R flag column labels, best first. gemmi's rfree_column() matches
# the first five; the rest are common variants seen in the wild. Compared
# case-insensitively.
_FREE_FLAG_LABELS = [
    'FreeR_flag', 'R-free-flags', 'FREE', 'RFREE', 'FREER',
    'FreeR', 'FreeRflag', 'TEST', 'pdbx_r_free_flag',
]


def free_flag_rank(label):
    '''
    Rank a candidate free-R flag column by its label. Lower is better; None means
    "not a recognised free-flag name" (the caller may still fall back to it as a
    last resort if it is the only integer column). An exact canonical match beats
    a mere 'free' substring.
    '''
    core = _strip_producer_prefix(label.strip())
    low = core.lower()
    for i, canonical in enumerate(_FREE_FLAG_LABELS):
        if low == _strip_producer_prefix(canonical).lower():
            return i
    if 'free' in low:
        return len(_FREE_FLAG_LABELS)
    return None


def rank_free_flags(names):
    '''
    Given the flag-column names present in a file, return them ordered
    best-candidate-first, restricted to those that look like free-R flags. If
    none match by name the list is empty (caller decides whether to offer all
    integer columns).
    '''
    scored = []
    for name in names:
        labels = extract_column_labels(name)
        # A flag dataset is a single integer column.
        label = labels[0] if labels else name
        rank = free_flag_rank(label)
        if rank is not None:
            scored.append((rank, name))
    scored.sort(key=lambda t: t[0])
    return [name for _, name in scored]


# --- experimental-data ranking -----------------------------------------------

# Canonical experimental amplitude/intensity base labels, best first. Used as a
# name bonus on top of the data-type priority.
_EXPERIMENTAL_LABELS = [
    'FP', 'F-obs-filtered', 'F-obs', 'FOBS', 'F',
    'IMEAN', 'I-obs', 'IOBS', 'I',
]


def _experimental_name_bonus(label):
    core = _strip_producer_prefix(label.strip())
    low = core.lower()
    for i, canonical in enumerate(_EXPERIMENTAL_LABELS):
        if low == canonical.lower():
            return i
    return len(_EXPERIMENTAL_LABELS)


def experimental_sort_key(name, dtype_priority):
    '''
    Sort key for choosing the experimental dataset to drive live maps. Primary
    key is the data-type priority (supplied by the caller, which knows the
    Clipper types); secondary key is the canonical-name bonus. Lower sorts first.
    '''
    labels = extract_column_labels(name)
    amp = labels[0] if labels else name
    return (dtype_priority, _experimental_name_bonus(amp))


# --- shared helpers ----------------------------------------------------------

def extract_column_labels(name):
    '''
    Recover the individual MTZ/CIF column labels from a dataset name. Handles the
    forms the loaders produce, including stacked qualifiers:
      '(dataset) FWT, PHWT'         -> ['FWT', 'PHWT']
      '(STATIC) (dataset) FWT, PHWT'-> ['FWT', 'PHWT']
      'pdbx_FWT, pdbx_PHWT'         -> ['pdbx_FWT', 'pdbx_PHWT']
    Any leading '(...)' qualifiers are stripped.
    '''
    s = name.strip()
    while s.startswith('('):
        close = s.find(')')
        if close == -1:
            break
        s = s[close + 1:].strip()
    return [part.strip() for part in s.split(',') if part.strip()]


def labels_match(data_col, sigma_or_phase_col):
    '''
    Quick-and-dirty check that two MTZ column names belong together: the second
    should equal the first plus a prefix or a suffix, but not both, and the
    added part must contain no digit (so we don't pair e.g. FOFC with PH2FOFC).
    This is the name-based fallback used when strict column *ordering* fails.
    (Historically lived in clipper_mtz.data_labels_match.)
    '''
    if data_col in sigma_or_phase_col:
        remainder = [s for s in sigma_or_phase_col.split(data_col) if s]
        if len(remainder) == 1:
            return not any(c.isdigit() for c in remainder[0])
    return False
