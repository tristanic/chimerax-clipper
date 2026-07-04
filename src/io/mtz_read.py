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

def _defaultdict_as_dict(d):
    '''
    Convert a recursively-defined defaultdict back to a normal dict for easier
    downstream handling.
    '''
    for k, v in d.items():
        if isinstance(v, dict):
            d[k] = _defaultdict_as_dict(v)
        return dict(d)

def _parse_column_path(cpath):
    path, dtype = cpath.split(' ')
    crystal, dataset, name = path.split('/')[1:]
    return (crystal, dataset, name, dtype)


def _plan_imports(columns, load_map_coeffs=True):
    '''
    Decide how to group a flat, ordered list of MTZ columns into importable
    datasets, WITHOUT touching Clipper (pure and unit-testable).

    Args:
        columns: ordered list of (crystal, dataset, name, dtype) tuples, with the
            Miller-index ('H') columns already removed.
        load_map_coeffs: if False, amplitude/phase (F/phi) map coefficients are
            recognised (so their columns aren't mis-paired) but omitted from the
            plan.

    Returns (plans, messages) where each plan is
    (crystal, dataset, [column names in order], kind) with kind one of
    'F_sigF', 'F_phi', 'D_sigD', 'F_sigF_ano', 'I_sigI', 'I_sigI_ano', 'Flag';
    and messages is a list of (level, text) for the caller to log.

    Adjacency of the expected following column type is the primary signal — this
    is how essentially every producer writes F/sigF and amplitude/phase pairs.
    When adjacency fails (a mis-ordered or non-adjacent sigma), we fall back to
    pairing the amplitude/intensity with a sigma column matched *by name*
    (SIG-prefix etc.) within the same crystal/dataset, so the amplitude is no
    longer silently dropped.
    '''
    from .reflection_classify import labels_match
    from collections import OrderedDict

    groups = OrderedDict()
    for idx, c in enumerate(columns):
        groups.setdefault((c[0], c[1]), []).append(idx)

    def nm(idx):
        return columns[idx][2]

    def dt(idx):
        return columns[idx][3]

    plans = []
    messages = []
    consumed = set()

    for (crystal, dataset), idxs in groups.items():
        def next_unconsumed(after_pos):
            for pos in range(after_pos + 1, len(idxs)):
                if idxs[pos] not in consumed:
                    return pos
            return None

        def find_sigma_by_name(base_name, want_type):
            for pos in range(len(idxs)):
                j = idxs[pos]
                if j in consumed:
                    continue
                if dt(j) == want_type and labels_match(base_name, nm(j)):
                    return pos
            return None

        def take_adjacent(start_pos, expected_types):
            positions = []
            p = start_pos
            for et in expected_types:
                p = next_unconsumed(p)
                if p is None or dt(idxs[p]) != et:
                    return None
                positions.append(p)
            return positions

        def emit(positions, kind):
            for p in positions:
                consumed.add(idxs[p])
            plans.append((crystal, dataset, [nm(idxs[p]) for p in positions], kind))

        for pos, idx in enumerate(idxs):
            if idx in consumed:
                continue
            d = dt(idx)
            if d == 'F':
                np = next_unconsumed(pos)
                if np is not None and dt(idxs[np]) == 'Q':
                    emit([pos, np], 'F_sigF')
                elif np is not None and dt(idxs[np]) == 'P':
                    if load_map_coeffs:
                        emit([pos, np], 'F_phi')
                    else:
                        consumed.update([idx, idxs[np]])
                else:
                    qpos = find_sigma_by_name(nm(idx), 'Q')
                    if qpos is not None:
                        emit([pos, qpos], 'F_sigF')
                        messages.append(('info',
                            'Paired amplitude column "{}" with non-adjacent sigma '
                            '"{}" by name.'.format(nm(idx), nm(idxs[qpos]))))
                    else:
                        ppos = find_sigma_by_name(nm(idx), 'P')
                        if ppos is not None:
                            if load_map_coeffs:
                                emit([pos, ppos], 'F_phi')
                            else:
                                consumed.update([idx, idxs[ppos]])
                        else:
                            consumed.add(idx)
                            messages.append(('warning',
                                'WARNING: amplitude column "/{}/{}/{}" could not be '
                                'matched to a sigma or phase column and was skipped.'
                                .format(crystal, dataset, nm(idx))))
            elif d in ('D', 'J'):
                kind = 'D_sigD' if d == 'D' else 'I_sigI'
                np = next_unconsumed(pos)
                if np is not None and dt(idxs[np]) == 'Q':
                    emit([pos, np], kind)
                else:
                    qpos = find_sigma_by_name(nm(idx), 'Q')
                    if qpos is not None:
                        emit([pos, qpos], kind)
                        messages.append(('info',
                            'Paired column "{}" with non-adjacent sigma "{}" by '
                            'name.'.format(nm(idx), nm(idxs[qpos]))))
                    else:
                        consumed.add(idx)
                        messages.append(('warning',
                            'WARNING: column "/{}/{}/{}" has no matching sigma '
                            'column and was skipped.'.format(crystal, dataset, nm(idx))))
            elif d == 'G':
                seq = take_adjacent(pos, ['L', 'G', 'L'])
                if seq is not None:
                    emit([pos] + seq, 'F_sigF_ano')
                else:
                    consumed.add(idx)
                    messages.append(('warning',
                        'WARNING: anomalous amplitude column "/{}/{}/{}" was not '
                        'followed by the expected sigF+/F-/sigF- columns and was '
                        'skipped.'.format(crystal, dataset, nm(idx))))
            elif d == 'K':
                seq = take_adjacent(pos, ['M', 'K', 'M'])
                if seq is not None:
                    emit([pos] + seq, 'I_sigI_ano')
                else:
                    consumed.add(idx)
                    messages.append(('warning',
                        'WARNING: anomalous intensity column "/{}/{}/{}" was not '
                        'followed by the expected sigI+/I-/sigI- columns and was '
                        'skipped.'.format(crystal, dataset, nm(idx))))
            elif d == 'I':
                emit([pos], 'Flag')
            else:
                consumed.add(idx)
                messages.append(('info',
                    'Discarding unrecognised/unsupported data array "/{}/{}/{} {}"'
                    .format(crystal, dataset, nm(idx), d)))

    return plans, messages


def load_mtz_data(session, filename, load_map_coeffs=True):
    from chimerax.clipper import (CCP4MTZfile, HKL_info,
        HKL_data_Flag, HKL_data_D_sigD, HKL_data_F_sigF, HKL_data_F_phi,
        HKL_data_F_sigF_ano, HKL_data_I_sigI, HKL_data_I_sigI_ano)

    kind_to_type = {
        'F_sigF':       HKL_data_F_sigF,
        'F_phi':        HKL_data_F_phi,
        'D_sigD':       HKL_data_D_sigD,
        'F_sigF_ano':   HKL_data_F_sigF_ano,
        'I_sigI':       HKL_data_I_sigI,
        'I_sigI_ano':   HKL_data_I_sigI_ano,
        'Flag':         HKL_data_Flag,
    }

    mtzin = CCP4MTZfile()
    hklinfo = HKL_info()
    mtzin.open_read(filename)
    mtzin.import_hkl_info(hklinfo, True)

    # Get all column names and types, dropping the Miller indices (already in
    # hklinfo).
    columns = []
    for cpath in (str(p) for p in mtzin.column_paths):
        crystal, dataset, name, dtype = _parse_column_path(cpath)
        if dtype == 'H':
            continue
        columns.append((crystal, dataset, name, dtype))

    plans, messages = _plan_imports(columns, load_map_coeffs=load_map_coeffs)
    for level, text in messages:
        getattr(session.logger, level)(text)

    from collections import defaultdict
    def nested_default_dict():
        return defaultdict(lambda: nested_default_dict())
    crystals = defaultdict(nested_default_dict)

    for crystal, dataset, names, kind in plans:
        dset_name = ', '.join(names)
        import_path = '{}/{}/[{}]'.format(crystal, dataset, dset_name)
        this_data = kind_to_type[kind](hklinfo)
        mtzin.import_hkl_data(this_data, import_path)
        crystals[crystal][dataset][dset_name] = this_data

    return hklinfo, _defaultdict_as_dict(crystals)
