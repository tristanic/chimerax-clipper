# Clipper plugin to UCSF ChimeraX
# Copyright (C) 2016-2019 Tristan Croll, University of Cambridge
#
# Offline regression tests for small-molecule (COD) CIF support: model + crystal
# context, structure-factor reading, and the bulk-solvent-free recomputed R-factor.
# Test data are bundled (no network needed):
#   cod_1100908  C2/c (#15), structure factors, Cu + water O on a 2-fold special
#                position - exercises special-position occupancy handling.
#   cod_2213867  Pbca (#61), structure factors - orthorhombic glide-plane group.
#   cod_2010010  Pbca (#61), model only (no .hkl) - the common COD case.
#
# Run inside ChimeraX, e.g.:
#   run_chimerax.bat --nogui --exit --script src/tests/test_small_molecule.py
# or call an individual test_*(session) from the ChimeraX command line.

import os

_DATA = os.path.abspath(os.path.dirname(__file__))


def _extract(session, name):
    from chimerax.clipper.io.small_molecule import extract_cod_structure
    return extract_cod_structure(session, os.path.join(_DATA, name))


def test_structure_factors_special_position(session):
    '''C2/c with a metal on a 2-fold: recomputed R must match the published R,
    which only works if special-position occupancies are correctly halved.'''
    p = _extract(session, 'cod_1100908.cif')
    assert p['space_group_number'] == 15, p['space_group']
    assert p['has_structure_factors']
    assert p['published_r_factor'] is not None
    # Recomputed bulk-solvent-free R reproduces the published R (~0.041). A wide
    # tolerance still catches the special-position / frame / ADP regressions, which
    # pushed this to 0.2-0.8 when broken.
    assert abs(p['recomputed_r_factor_all'] - p['published_r_factor']) < 0.03, \
        (p['recomputed_r_factor_all'], p['published_r_factor'])
    assert len(p['elements']) == 25


def test_structure_factors_pbca(session):
    '''Pbca (centrosymmetric, glide planes) structure factors.'''
    p = _extract(session, 'cod_2213867.cif')
    assert p['space_group_number'] == 61, p['space_group']
    assert p['has_structure_factors']
    assert abs(p['recomputed_r_factor_all'] - p['published_r_factor']) < 0.03, \
        (p['recomputed_r_factor_all'], p['published_r_factor'])


def test_model_only(session):
    '''Model-only entry (no reflections): structure-factor fields are None, but
    the published metrics still come through.'''
    p = _extract(session, 'cod_2010010.cif')
    assert p['space_group_number'] == 61, p['space_group']
    assert p['has_structure_factors'] is False
    assert p['recomputed_r_factor'] is None
    assert p['published_r_factor'] is not None
    assert p['hydrogen_treatment'] is not None
    assert p['collection_temperature'] is not None


def test_clipper_frame_geometry(session):
    '''Payload coordinates are rebuilt from the CIF fractionals via Clipper, and
    must reproduce the CIF's own published _geom_bond_distance values (the corecif
    oblique-cell coordinate workaround).'''
    import numpy
    from chimerax.mmcif import get_cif_tables
    from chimerax.clipper.io.small_molecule import _element_from_type_symbol
    cif = os.path.join(_DATA, 'cod_1100908.cif')
    p = _extract(session, cif)
    coords = numpy.asarray(p['coordinates'])
    # Map atom_site labels -> payload index via the atom order corecif preserves.
    at, _an, gb = get_cif_tables(cif, ['atom_site', 'atom_site_aniso', 'geom_bond'])
    labels = [r[0] for r in at.fields(('label',))]
    lab2i = {labels[i]: i for i in range(len(labels))}
    worst = 0.0
    n = 0
    for r in gb.fields(('atom_site_label_1', 'atom_site_label_2', 'distance',
                        'site_symmetry_2'), allow_missing_fields=True):
        if len(r) > 3 and r[3] not in ('.', ''):
            continue  # skip symmetry-generated bonds
        if r[0] in lab2i and r[1] in lab2i and len(coords) == len(labels):
            d = numpy.linalg.norm(coords[lab2i[r[0]]] - coords[lab2i[r[1]]])
            worst = max(worst, abs(d - float(r[2].split('(')[0])))
            n += 1
    assert n > 0
    assert worst < 0.01, worst   # Clipper-frame coords match to <0.001 A; corecif's were ~0.017


def test_live_map_engine(session):
    '''The live-map compute engine (SmallMoleculeXmapMgr): the 2Fo-Fc map must
    place atoms on positive density, and the anisotropic+spline scaling must give an
    R-work close to the published value (which a crude overall scale would not - it
    left a large heavy-atom residual). Exercised headlessly; the GUI display itself
    needs an OpenGL context.'''
    import numpy
    from chimerax.clipper.symmetry import crystal_symmetry_from_cif_file
    from chimerax.clipper.io.small_molecule import (open_small_molecule_cif,
        hydrate_small_molecule_model, _small_molecule_map_data)
    from chimerax.clipper.maps.small_molecule_map import SmallMoleculeXmapMgr
    from chimerax.clipper.clipper_python import Coord_orth, Map_stats
    cif = os.path.join(_DATA, 'cod_1100908.cif')   # C2/c, Cu on a 2-fold; published R 0.041
    model = open_small_molecule_cif(session, cif)   # coords already correct (corecif + open)
    try:
        cell, sg, grid = crystal_symmetry_from_cif_file(cif)
        hydrate_small_molecule_model(session, model, cif, cell, 'xray')
        smd = _small_molecule_map_data(model, cif, None, cell, sg, grid)
        mgr = SmallMoleculeXmapMgr(smd['hklinfo'], smd['cell'], smd['spacegroup'],
            smd['grid'], smd['fobs'], smd['structure'])
        mgr.add_xmap('2Fo-Fc', is_difference_map=False)
        mgr.add_xmap('Fo-Fc', is_difference_map=True)
        # R-work close to the published 0.041 (aniso+spline scaling).
        assert abs(mgr.rwork - 0.041) < 0.02, mgr.rwork
        # 2Fo-Fc places atoms on positive density.
        xm = mgr.get_xmap_ref('2Fo-Fc')
        st = Map_stats(xm)
        vals = numpy.array([xm.get_data(
            Coord_orth(float(x[0]), float(x[1]), float(x[2])).coord_frac(cell).coord_grid(
                xm.grid_sampling)) for x in model.atoms.coords])
        assert (vals > st.mean).mean() > 0.9, (vals > st.mean).mean()
        assert (vals.mean() - st.mean) / st.std_dev > 2.0
    finally:
        session.models.close([model])


def test_electron_scattering_selected(session):
    '''Electron scattering factors (micro-ED): selecting radiation='electron' must
    actually swap the scattering-factor table, so the recomputed R differs from the
    X-ray value. (On this X-ray dataset the electron R is far worse - the point is
    only that the electron table is genuinely used; correctness of the electron
    coefficients is validated separately against Int. Tab. Vol C.)'''
    from chimerax.clipper.io.small_molecule import extract_cod_structure
    cif = os.path.join(_DATA, 'cod_1100908.cif')
    px = extract_cod_structure(session, cif, radiation='xray')
    pe = extract_cod_structure(session, cif, radiation='electron')
    assert px['radiation'] == 'xray' and pe['radiation'] == 'electron'
    assert px['recomputed_r_factor'] is not None
    assert pe['recomputed_r_factor'] is not None
    # X-ray reproduces the published R; electron (wrong for X-ray data) differs.
    assert px['recomputed_r_factor_all'] < 0.06, px['recomputed_r_factor_all']
    assert abs(px['recomputed_r_factor'] - pe['recomputed_r_factor']) > 1e-2, \
        (px['recomputed_r_factor'], pe['recomputed_r_factor'])


def test_radiation_autodetect(session):
    '''radiation='auto' reads _diffrn_radiation_probe from the CIF; a standard X-ray
    entry (no electron probe) resolves to X-ray.'''
    from chimerax.clipper.io.small_molecule import _radiation_from_cif
    assert _radiation_from_cif(os.path.join(_DATA, 'cod_1100908.cif')) == 'xray'


def test_electron_map_engine(session):
    '''The live-map engine runs under electron radiation (thread + GIL release) and
    produces a finite R-work that differs from the X-ray one - i.e. the electron
    table reaches the FFT Fcalc path, not just the summation R-factor path.'''
    from chimerax.clipper.symmetry import crystal_symmetry_from_cif_file
    from chimerax.clipper.io.small_molecule import (open_small_molecule_cif,
        hydrate_small_molecule_model, _small_molecule_map_data)
    from chimerax.clipper.maps.small_molecule_map import SmallMoleculeXmapMgr
    cif = os.path.join(_DATA, 'cod_1100908.cif')
    model = open_small_molecule_cif(session, cif)   # coords already correct (corecif + open)
    try:
        cell, sg, grid = crystal_symmetry_from_cif_file(cif)
        def rwork(radiation):
            hydrate_small_molecule_model(session, model, cif, cell, radiation)
            smd = _small_molecule_map_data(model, cif, None, cell, sg, grid, radiation)
            mgr = SmallMoleculeXmapMgr(smd['hklinfo'], smd['cell'], smd['spacegroup'],
                smd['grid'], smd['fobs'], smd['structure'], radiation=smd['radiation'])
            mgr.add_xmap('2Fo-Fc', is_difference_map=False)
            return mgr.rwork
        rx, re = rwork('xray'), rwork('electron')
        assert rx > 0 and re > 0 and abs(rx - re) > 1e-2, (rx, re)
    finally:
        session.models.close([model])


def test_ionic_electron_factors(session):
    '''Ionic electron scattering factors (Peng 1998): the compiled AtomShapeFn must
    reproduce eq (5) = screened Gaussians + 0.023934*dZ/s^2 exactly, the Coulomb
    term must NOT leak into X-ray, and Peng-only ions (absent from the X-ray table,
    e.g. Cr4+) must resolve under electron radiation.'''
    import numpy
    from chimerax.clipper.clipper_python import AtomShapeFn, Coord_orth
    K = 0.023934
    # (dZ, a[5], b[5]) from Peng 1998 Table 1.
    peng = {
        'O2-':  (-2, [0.0421,0.210,0.852,1.82,1.17],   [0.0609,0.559,2.96,11.5,37.7]),
        'Cu2+': (+2, [0.224,0.544,0.970,0.727,0.182],  [0.145,0.933,2.69,7.11,19.4]),
    }
    s = numpy.array([0.05, 0.15, 0.3, 0.6, 1.0])
    for ion, (dz, a, b) in peng.items():
        asf = AtomShapeFn(Coord_orth(0, 0, 0), ion, 0.0, 1.0, AtomShapeFn.ELECTRON)
        got = numpy.array([asf.f(4.0 * x * x) for x in s])
        ref = sum(a[i] * numpy.exp(-b[i] * s * s) for i in range(5)) + K * dz / (s * s)
        assert numpy.max(numpy.abs(got - ref)) < 1e-6, (ion, got, ref)
    # Cr4+ is not in Clipper's X-ray table but is in Peng's electron table.
    cr = AtomShapeFn(Coord_orth(0, 0, 0), 'Cr4+', 0.0, 1.0, AtomShapeFn.ELECTRON)
    assert cr.f(4.0 * 0.3 * 0.3) > 0
    # X-ray must stay finite at low s (no ionic Coulomb 1/s^2 term).
    xr = AtomShapeFn(Coord_orth(0, 0, 0), 'Cu2+', 0.0, 1.0, AtomShapeFn.XRAY)
    assert xr.f(4.0 * 0.02 * 0.02) < 30.0


def _pbca_arrays():
    from chimerax.clipper.symmetry import crystal_symmetry_from_cif_file
    from chimerax.clipper.io.small_molecule import _parse_reflection_file
    path = os.path.join(_DATA, 'cod_2213867.cif')
    cell, sg, grid = crystal_symmetry_from_cif_file(path)
    hkl, fsq, sig = _parse_reflection_file(path)
    return cell, sg, hkl, fsq, sig


def _slot_values(fobs, hkl):
    '''(f, sigf) held in the slot of each hkl; NaN where the slot is empty.'''
    import numpy
    from chimerax.clipper.clipper_python import HKL
    out = numpy.full((len(hkl), 2), numpy.nan)
    for i, h in enumerate(hkl):
        d = fobs[HKL(h.tolist())]
        if not d.missing:
            out[i] = d.f, d.sigf
    return out


# NB: keep each fobs_from_arrays HKL_info bound for as long as its HKL_data is used - the
# HKL_data hold a non-owning pointer to it (unpacking it into `_` leaves them dangling).
def _filled(fobs):
    import numpy
    return int(numpy.isfinite(fobs.data[1][:, 0]).sum())


def test_merge_equivalents_inverse_variance(session):
    '''Symmetry- and Friedel-equivalent copies of a reflection merge into its one slot as
    the inverse-variance mean of I, with the propagated sigma.'''
    import numpy
    from chimerax.clipper.io.small_molecule import fobs_from_arrays
    cell, sg, hkl, fsq, sig = _pbca_arrays()
    rot = numpy.rint(sg.primitive_symop(1).rot.as_numpy()).astype(numpy.int32)
    n = 200
    sub = hkl[:n]
    i_rows = numpy.stack([fsq[:n], fsq[:n] * 1.3 + 5.0, fsq[:n] * 0.8 + 2.0])
    s_rows = numpy.stack([sig[:n], sig[:n] * 2.0, sig[:n] * 0.5])
    hkl_d = numpy.concatenate([hkl, sub @ rot, -sub])
    fsq_d = numpy.concatenate([fsq, i_rows[1], i_rows[2]])
    sig_d = numpy.concatenate([sig, s_rows[1], s_rows[2]])

    ref_info, ref, _, _ = fobs_from_arrays(hkl, fsq, sig, cell, sg, merge_equivalents=True)
    info, fobs, _, _ = fobs_from_arrays(hkl_d, fsq_d, sig_d, cell, sg, merge_equivalents=True)
    assert _filled(fobs) == _filled(ref) == len(hkl), (_filled(fobs), _filled(ref), len(hkl))

    w = 1.0 / s_rows**2
    i_m = (w * i_rows).sum(0) / w.sum(0)
    s_m = 1.0 / numpy.sqrt(w.sum(0))
    got = _slot_values(fobs, sub)
    pos = i_m > 0
    assert numpy.allclose(got[pos, 0], numpy.sqrt(i_m[pos]), rtol=1e-6)
    assert numpy.allclose(got[pos, 1], s_m[pos] / (2 * numpy.sqrt(i_m[pos])), rtol=1e-6)
    assert numpy.all(got[~pos, 0] == 0)
    # Reflections without copies are untouched.
    assert numpy.array_equal(_slot_values(fobs, hkl[n:]), _slot_values(ref, hkl[n:]))


def test_merge_friedel_noncentrosymmetric(session):
    '''In a non-centrosymmetric group, Friedel mates h and -h share one slot and merge.'''
    import numpy
    from chimerax.clipper import Spacegroup, Spgr_descr
    from chimerax.clipper.io.small_molecule import fobs_from_arrays
    cell, _, hkl, fsq, sig = _pbca_arrays()
    sg = Spacegroup(Spgr_descr('P 2ac 2ab', Spgr_descr.Hall))     # P 21 21 21
    n = 300
    sub = hkl[:n]
    i_plus, i_minus = fsq[:n] + 10.0, fsq[:n] * 1.2 + 10.0
    s_plus, s_minus = sig[:n], sig[:n] * 3.0
    info, fobs, _, _ = fobs_from_arrays(
        numpy.concatenate([sub, -sub]), numpy.concatenate([i_plus, i_minus]),
        numpy.concatenate([s_plus, s_minus]), cell, sg, merge_equivalents=True)
    assert _filled(fobs) == n, _filled(fobs)
    w_p, w_m = 1 / s_plus**2, 1 / s_minus**2
    i_m = (w_p * i_plus + w_m * i_minus) / (w_p + w_m)
    assert numpy.allclose(_slot_values(fobs, sub)[:, 0], numpy.sqrt(i_m), rtol=1e-6)


def test_merge_orbit_matches_clipper_find_sym(session):
    '''The numpy Laue-orbit grouping agrees with Clipper's own reflection-to-slot mapping
    across crystal systems (hexagonal/trigonal included, where a fractional rotation's
    transpose is not itself a group operator).'''
    import numpy
    from chimerax.clipper import (Cell, Cell_descr, Spacegroup, Spgr_descr, HKL_info,
                                  HKL_data_F_sigF)
    from chimerax.clipper.clipper_python import Resolution
    from chimerax.clipper.io.small_molecule import (merge_equivalent_reflections,
                                                    fobs_from_arrays)
    hexagonal = (9.0, 9.0, 13.0, 90, 90, 120)
    cases = [
        ('-P 1', (9.0, 10.0, 11.0, 80, 85, 95)),
        ('-C 2yc', (9.0, 10.0, 11.0, 90, 105, 90)),
        ('P 4nw 2abw', (9.0, 9.0, 13.0, 90, 90, 90)),
        ('-R 3', hexagonal),
        ('-P 3 2c', hexagonal),
        ('P 61', hexagonal),
        ('P 2ac 2ab 3', (11.0, 11.0, 11.0, 90, 90, 90)),
    ]
    rng = numpy.random.default_rng(10061865)
    for hall, dims in cases:
        cell = Cell(Cell_descr(*dims))
        sg = Spacegroup(Spgr_descr(hall, Spgr_descr.Hall))
        hi = HKL_info(sg, cell, Resolution(1.5), True)
        asu = HKL_data_F_sigF(hi).data[0]
        asu = asu[asu.any(axis=1)]      # the list includes F000, which is never measured
        n = len(asu)
        rots = [numpy.rint(sg.primitive_symop(i).rot.as_numpy()).astype(numpy.int32)
                for i in range(sg.num_primitive_symops)]
        ops = rng.integers(len(rots), size=n)
        signs = rng.choice([-1, 1], size=n)
        rows = numpy.array([s * (h @ rots[o]) for h, o, s in zip(asu, ops, signs)],
                           numpy.int32)
        fsq = numpy.arange(n) + 1.0
        hkl_m, _, _, mult = merge_equivalent_reflections(rows, fsq, numpy.ones(n), sg)
        assert len(hkl_m) == n and numpy.all(mult == 1), (hall, len(hkl_m), n)
        info, fobs, _, _ = fobs_from_arrays(rows, fsq, numpy.ones(n), cell, sg,
                                            merge_equivalents=True)
        got = _slot_values(fobs, asu)[:, 0]
        assert numpy.allclose(got, numpy.sqrt(fsq), rtol=1e-6), hall
        del fobs, info


def test_merge_drops_unusable_sigma(session):
    '''Rows with no usable sigma do not enter the merge; with no usable sigma anywhere,
    the merge falls back to unit weights. The legacy path treats an unknown sigma as 1.'''
    import numpy
    from chimerax.clipper.io.small_molecule import fobs_from_arrays
    cell, sg, hkl, fsq, sig = _pbca_arrays()
    a, b = hkl[0], hkl[1]
    rows = numpy.concatenate([[a, -a, b, -b], hkl[2:]])
    i_rows = numpy.concatenate([[100.0, 300.0, 50.0, 70.0], fsq[2:]])
    s_rows = numpy.concatenate([[0.0, 10.0, numpy.nan, numpy.nan], sig[2:]])
    info, fobs, _, _ = fobs_from_arrays(rows, i_rows, s_rows, cell, sg, merge_equivalents=True)
    got = _slot_values(fobs, numpy.array([a, b]))
    assert numpy.isclose(got[0, 0], numpy.sqrt(300.0), rtol=1e-6), got
    assert numpy.isclose(got[0, 1], 10.0 / (2 * numpy.sqrt(300.0)), rtol=1e-6), got
    assert numpy.isnan(got[1, 0]), got

    n = 50
    nan = numpy.full(2 * n, numpy.nan)
    i_two = numpy.concatenate([fsq[:n] + 1.0, fsq[:n] + 3.0])
    info2, fobs2, _, _ = fobs_from_arrays(numpy.concatenate([hkl[:n], -hkl[:n]]), i_two,
                                          nan, cell, sg, merge_equivalents=True)
    i_m = fsq[:n] + 2.0
    got = _slot_values(fobs2, hkl[:n])
    assert numpy.allclose(got[:, 0], numpy.sqrt(i_m), rtol=1e-6)
    assert numpy.allclose(got[:, 1], (1 / numpy.sqrt(2)) / (2 * numpy.sqrt(i_m)), rtol=1e-6)

    unknown_info, unknown, _, _ = fobs_from_arrays(
        hkl, fsq, numpy.full(len(hkl), numpy.nan), cell, sg)
    unit_info, unit, _, _ = fobs_from_arrays(hkl, fsq, numpy.ones(len(hkl)), cell, sg)
    assert numpy.array_equal(unknown.data[1], unit.data[1], equal_nan=True)


def test_boundary_reflection_kept(session):
    '''The reflection that defines d_min gets a slot (the resolution limit is padded).'''
    import numpy
    from chimerax.clipper.symmetry import crystal_symmetry_from_cif_file
    from chimerax.clipper.clipper_python import HKL
    from chimerax.clipper.io.small_molecule import _parse_reflection_file, fobs_from_arrays
    for name in ('cod_1100908.cif', 'cod_2213867.cif'):
        path = os.path.join(_DATA, name)
        cell, sg, grid = crystal_symmetry_from_cif_file(path)
        hkl, fsq, sig = _parse_reflection_file(path)
        info, fobs, _, _ = fobs_from_arrays(hkl, fsq, sig, cell, sg, merge_equivalents=True)
        assert _filled(fobs) == len(hkl), (name, _filled(fobs), len(hkl))
        s = numpy.array([HKL(h.tolist()).invresolsq(cell) for h in hkl])
        edge = hkl[s == s.max()]
        assert not numpy.isnan(_slot_values(fobs, edge)[:, 0]).any(), (name, edge)
        del fobs, info


def run_all(session):
    tests = [v for k, v in sorted(globals().items())
             if k.startswith('test_') and callable(v)]
    failures = []
    for t in tests:
        try:
            t(session)
            print('PASS %s' % t.__name__)
        except Exception as e:
            failures.append((t.__name__, e))
            print('FAIL %s: %r' % (t.__name__, e))
    print('\n%d/%d passed' % (len(tests) - len(failures), len(tests)))
    if failures:
        raise SystemExit(1)


# When run via `run_chimerax --script`, `session` is available as a global.
if 'session' in globals():
    run_all(session)
