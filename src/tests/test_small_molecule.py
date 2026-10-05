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
#   cod_2021365  P-1 (#2), structure factors, SHELXL EXTI x = 0.022, and the raw HKLF 4
#                data embedded as _shelx_hkl_file - the extinction convention check.
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


# ---------------------------------------------------------------------------------------
# Opt-in structure-factor target variants (io.sf_target): intensities, SHELXL weights,
# extinction. Weighting strings marked "corpus" are verbatim from COD CIFs (2104628,
# 2212865, 2208469, 2018789, 2007192, 2201204), as in garnet's B13 scan tests.
# ---------------------------------------------------------------------------------------

def test_parse_shelxl_weighting(session):
    from chimerax.clipper.io.sf_target import parse_shelxl_weighting
    shelxl = [
        (r"calc w=1/[\s^2^(Fo^2^)+(0.0260P)^2^+1.1430P] where P=(Fo^2^+2Fc^2^)/3", 0.0260, 1.1430, False),
        (r"w = 1/[\s^2^(F^2^) + (0.04P)^2^ + 0.14P], where P=[max(Fo^2^,0) + 2Fc^2^]/3", 0.04, 0.14, True),
        (r"calc w=1/[\s^2^(Fo^2^)+(0.0474P)^2^] where P=(Fo^2^+2Fc^2^)/3", 0.0474, 0.0, False),
        (r"calc w = 1/[\s^2^(Fo^2^)+(0.0446P)^2^+0.0000P] where P=(Fo^2^+2Fc^2^)/3", 0.0446, 0.0, False),
        (r"w = 1/[\s^2^(Fo^2^)+(0.1106P)^2^+0.4663P] where P = (Fo^2^+2Fc^2^)/3", 0.1106, 0.4663, False),
        (r"w=1/[\s^2^(F~o~^2^)+(0.0412*P)^2^+0.2563*P]", 0.0412, 0.2563, None),
        (r"w=1/[\s^2^(Fo^2^)+0.2563P] where P=(Fo^2^+2Fc^2^)/3", 0.0, 0.2563, False),
        # the three bundled CIFs (cod_2010010's is in the old _refine_ls_weighting_scheme)
        (r" w = 1/[\s^2^(Fo^2^)+(0.04P)^2^+2P] where P = (Fo^2^+2Fc^2^)/3", 0.04, 2.0, False),
        (r"calc w = 1/[\s^2^(Fo^2^)+(0.099P)^2^+0.48P] where P=(Fo^2^+2Fc^2^)/3", 0.099, 0.48, False),
        ("calc w = 1/[\\s^2^(F~o~^2^)+(0.0362P)^2^+1.0198P] where\nP = (F~o~^2^+2F~c~^2^)/3",
         0.0362, 1.0198, False),
        # SHELXL's P in its manual's order, and without carets
        (r"w=1/[\s^2^(Fo^2^)+(0.05P)^2^+0.3P] where P=[2Fc^2^+Max(Fo^2^,0)]/3", 0.05, 0.3, True),
        (r"w=1/[\s^2^(Fo^2^)+(0.05P)^2^+0.3P] where P=(Fo2+2Fc2)/3", 0.05, 0.3, False),
    ]
    for text, a, b, p_max0 in shelxl:
        w = parse_shelxl_weighting(text)
        assert w['form'] == 'shelxl', (text, w)
        assert abs(w['a'] - a) < 1e-12 and abs(w['b'] - b) < 1e-12, (text, w)
        assert w['p_max0'] is p_max0, (text, w)
    for text, form in [
            (None, 'missing'), ('', 'missing'), ('?', 'missing'),
            (r"w=1/[\s^2^(Fo^2^)]", 'sigma_only'),
            (r"calc w=1/\s^2^(Fo^2^)", 'sigma_only'),
            (r"w=1/[\s^2^(F)+0.0004F^2^]", 'other'),                     # refined on F
            ("Method, part 1, Chebychev polynomial, (Watkin 1994, Prince 1982) [weight] = "
             "1.0/[A~0~*T~0~(x)+A~1~*T~1~(x)]", 'other'),
            (r"w=1/[\s^2^(Fo^2^)+(0.04P)^2^+0.1P+0.002Fo^2^]", 'other'),  # an extra term
            (r"w=1/[\s^2^(Fo^2^)+(0.04P)^2^+0.1P] where P=0.5Fo^2^+0.5Fc^2^", 'other')]:
        w = parse_shelxl_weighting(text)
        assert w['form'] == form, (text, w)
        if form == 'other':
            assert w['a'] is None and w['b'] is None


def test_parse_extinction(session):
    from chimerax.clipper.io.sf_target import parse_extinction, cif_number
    expr = r"Fc^*^=kFc[1+0.001xFc^2^\l^3^/sin(2\q)]^-1/4^"
    assert parse_extinction('SHELXL97', expr, '0.013(2)') == {'kind': 'shelxl', 'x': 0.013}
    assert parse_extinction('(SHELXL2018; Sheldrick, 2015)', None, '0.0012(3)')['kind'] == 'shelxl'
    assert parse_extinction('<i>SHELXTL</i> (Sheldrick, 2008)', expr, '0.097(8)')['kind'] == 'shelxl'
    assert parse_extinction(None, expr, '0.05')['kind'] == 'shelxl'
    # cod_2010010: the coefficient carries its own "x =", the formula sits in the method
    assert parse_extinction("F~c~^*^ = kF~c~[1+0.001xF~c~^2^\\l^3^/sin(2\\q)]^-1/4^", None,
                            'x = 0.0066(3)') == {'kind': 'shelxl', 'x': 0.0066}
    assert parse_extinction('Becker & Coppens (1974)', None, '120')['kind'] == 'other'
    assert parse_extinction('none', None, None)['kind'] == 'none'
    assert parse_extinction('SHELXL97', expr, '0.0')['kind'] == 'none'
    assert cif_number('1.2e-3(4)') == 1.2e-3 and cif_number('?') is None


def test_sf_target_spec(session):
    from chimerax.clipper.io.sf_target import SFTarget, check_target_args
    assert SFTarget().token() == '' and SFTarget().is_default
    assert SFTarget('intensity', 'shelxl', 'shelxl').token() == '|sft1-i-wshelxl-xshelxl'
    assert SFTarget('intensity_masked').token() == '|sft1-im-wsigma-xnone'
    assert SFTarget('french_wilson').kind == 'amplitude'
    assert SFTarget('intensity_masked').kind == 'intensity'
    assert SFTarget('amplitude', 'shelxl') == SFTarget('amplitude', 'shelxl', 'none')
    assert len({t.token() for t in _all_targets()}) == 16
    for bad in (dict(observations='f2'), dict(weights='unit'), dict(extinction='becker')):
        try:
            SFTarget(**bad)
        except ValueError:
            continue
        raise AssertionError(bad)
    assert check_target_args(None, 'intensity', False) == 'intensity'
    assert check_target_args(SFTarget('intensity'), 'amplitude', True) == 'intensity'
    for args in ((SFTarget('intensity'), 'amplitude', False),
                 (SFTarget('intensity'), 'intensity', True)):
        try:
            check_target_args(*args)
        except ValueError:
            continue
        raise AssertionError(args)


def _all_targets():
    from chimerax.clipper.io.sf_target import SFTarget
    return [SFTarget(o, w, e) for o in ('amplitude', 'intensity_masked', 'french_wilson',
                                        'intensity')
            for w in ('sigma', 'shelxl') for e in ('none', 'shelxl')]


def test_refinement_record(session):
    '''The depositor's record for the bundled crystals, row-aligned with the raw arrays.'''
    import numpy
    from chimerax.clipper.io.small_molecule import _parse_reflection_file
    from chimerax.clipper.io.sf_target import read_small_molecule_refinement
    for name, a, b in (('cod_1100908.cif', 0.04, 2.0), ('cod_2213867.cif', 0.099, 0.48)):
        path = os.path.join(_DATA, name)
        rec = read_small_molecule_refinement(path)
        hkl, fsq, sig = _parse_reflection_file(path)
        assert numpy.array_equal(rec['hkl'], hkl), name
        assert len(rec['fc_sq_dep']) == len(hkl) and numpy.isfinite(rec['fc_sq_dep']).all()
        w = rec['weighting']
        assert w['form'] == 'shelxl' and (w['a'], w['b']) == (a, b), (name, w)
        assert rec['extinction']['kind'] == 'none', rec['extinction']
        assert abs(rec['wavelength'] - 0.71073) < 1e-9, rec['wavelength']
        assert rec['f_squared_multiplier'] == 1.0
        assert rec['shelx_refln_list_code'] == 4
        assert rec['structure_factor_coef'] == 'Fsqd'


def test_shelxl_weights_reproduce_goof(session):
    '''With the depositor's Fc^2 and the parsed SHELXL weights, the depositor's own wR2
    and goodness of fit come back - an end-to-end check of a, b, P and the effective
    sigma that does not involve Clipper's Fcalc at all.'''
    import numpy
    from chimerax.clipper.symmetry import crystal_symmetry_from_cif_file
    from chimerax.clipper.io.small_molecule import _parse_reflection_file
    from chimerax.clipper.io.sf_target import (SFTarget, read_small_molecule_refinement,
        target_observation_arrays, _lookup_by_orbit)
    for name in ('cod_1100908.cif', 'cod_2213867.cif'):
        path = os.path.join(_DATA, name)
        cell, sg, grid = crystal_symmetry_from_cif_file(path)
        hkl, fsq, sig = _parse_reflection_file(path)
        aux = read_small_molecule_refinement(path)
        aux['fc_sq_model_abs'] = aux['fc_sq_model_scaled'] = numpy.full(len(hkl), numpy.nan)
        arr = target_observation_arrays(hkl, fsq, sig, cell, sg,
                                        SFTarget('intensity', 'shelxl'), aux)
        assert arr['flags']['weights_applied'] and arr['flags']['p_source'] == 'depositor'
        fc2 = _lookup_by_orbit(arr['hkl'], hkl, aux['fc_sq_dep'], sg)
        r2 = arr['w'] * (arr['i'] - fc2) ** 2
        wr2 = numpy.sqrt(r2.sum() / (arr['w'] * arr['i'] ** 2).sum())
        goof = numpy.sqrt(r2.sum() / (len(r2) - aux['n_parameters']))
        assert abs(wr2 / aux['wr2'] - 1) < 0.03, (name, wr2, aux['wr2'])
        assert abs(goof / aux['goof'] - 1) < 0.03, (name, goof, aux['goof'])
        # and counting statistics alone do not
        r2s = (arr['i'] - fc2) ** 2 / target_observation_arrays(
            hkl, fsq, sig, cell, sg, SFTarget('intensity'))['sig'] ** 2
        assert numpy.sqrt(r2s.sum() / (len(r2) - aux['n_parameters'])) > 1.2 * aux['goof']


def _synthetic_aux(n, fc_dep, fc_abs, a=0.05, b=0.3, m=1.0, ext_kind='shelxl', x=0.02,
                   form='shelxl', wavelength=0.71073):
    import numpy
    return {'fc_sq_dep': numpy.asarray(fc_dep, float), 'fc_sq_model_abs': numpy.asarray(fc_abs, float),
            'fc_sq_model_scaled': numpy.full(n, 7.0),
            'weighting': {'form': form, 'a': a, 'b': b, 'p_max0': False},
            'extinction': {'kind': ext_kind, 'x': x, 'method': None},
            'wavelength': wavelength, 'f_squared_multiplier': m}


def test_weight_and_extinction_arithmetic(session):
    '''sigma_eff^2 = sigma^2 + (aP)^2 + b m P with P = (max(I,0) + 2Fc^2)/3; then I and
    sigma divide by y = [1 + 0.001 x |Fc|^2 lambda^3 / sin(2 theta)]^-1/2.'''
    import numpy
    from chimerax.clipper.clipper_python import HKL
    from chimerax.clipper.io.sf_target import (SFTarget, target_observation_arrays,
                                               _lookup_by_orbit)
    cell, sg, hkl, fsq, sig = _pbca_arrays()
    hkl, fsq, sig = hkl[:40], fsq[:40].copy(), sig[:40]
    fsq[:5] = -3.0                           # negative intensities enter P as max(I, 0)
    n = len(hkl)
    fc_dep = numpy.linspace(10.0, 500.0, n)
    fc_abs = numpy.linspace(100.0, 9000.0, n)
    a, b, m, x, wl = 0.05, 0.3, 2.5, 0.02, 0.71073
    aux = _synthetic_aux(n, fc_dep, fc_abs, a, b, m, x=x, wavelength=wl)
    arr = target_observation_arrays(hkl, fsq, sig, cell, sg,
                                    SFTarget('intensity', 'shelxl', 'shelxl'), aux)
    assert len(arr['hkl']) == n and numpy.all(arr['multiplicity'] == 1)   # already unique

    def merged(per_row):       # per raw row -> the merged (canonical) order
        return _lookup_by_orbit(arr['hkl'], hkl, per_row, sg)
    p = (numpy.maximum(fsq, 0) + 2 * fc_dep) / 3
    s_eff = numpy.sqrt(sig ** 2 + (a * p) ** 2 + b * m * p)
    s = numpy.array([HKL(h.tolist()).invresolsq(cell) for h in hkl])
    sin_t = wl * numpy.sqrt(s) / 2
    y = (1 + 0.001 * x * fc_abs * wl ** 3 / (2 * sin_t * numpy.sqrt(1 - sin_t ** 2))) ** -0.5
    assert numpy.allclose(arr['y'], merged(y), rtol=1e-12)
    assert numpy.allclose(arr['i'], merged(fsq / y), rtol=1e-12)
    assert numpy.allclose(arr['sig'], merged(s_eff / y), rtol=1e-12)
    assert numpy.allclose(arr['w'], merged((y / s_eff) ** 2), rtol=1e-12)
    assert abs(arr['sum_w_obs2'] - (arr['w'] * arr['i'] ** 2).sum()) < 1e-9 * arr['sum_w_obs2']
    # each correction alone
    w_only = target_observation_arrays(hkl, fsq, sig, cell, sg, SFTarget('intensity', 'shelxl'), aux)
    assert numpy.allclose(w_only['sig'], merged(s_eff), rtol=1e-12)
    assert numpy.all(w_only['y'] == 1)
    x_only = target_observation_arrays(hkl, fsq, sig, cell, sg,
                                       SFTarget('intensity', extinction='shelxl'), aux)
    assert numpy.allclose(x_only['sig'], merged(sig / y), rtol=1e-12)
    assert numpy.allclose(x_only['i'], merged(fsq / y), rtol=1e-12)
    # the default convention, untouched
    plain = target_observation_arrays(hkl, fsq, sig, cell, sg, SFTarget('intensity'))
    assert numpy.allclose(plain['i'], merged(fsq), rtol=1e-12)
    assert numpy.allclose(plain['sig'], merged(sig), rtol=1e-12)
    # a SHELXL .fcf already carries the correction: nothing further is applied
    aux['shelx_refln_list_code'] = 4
    fcf = target_observation_arrays(hkl, fsq, sig, cell, sg,
                                    SFTarget('intensity', extinction='shelxl'), aux)
    assert fcf['flags']['extinction_in_data'] and not fcf['flags']['extinction_applied']
    assert not fcf['flags']['reasons'] and numpy.all(fcf['y'] == 1)
    assert numpy.array_equal(fcf['i'], plain['i']) and numpy.array_equal(fcf['sig'], plain['sig'])


def test_correction_fallback_flags(session):
    '''Corrections that cannot be applied fall back (sigma-only weights, y = 1) and say so.'''
    import numpy
    from chimerax.clipper.io.sf_target import SFTarget, target_observation_arrays
    cell, sg, hkl, fsq, sig = _pbca_arrays()
    hkl, fsq, sig = hkl[:40], fsq[:40], sig[:40]
    n = len(hkl)
    t = SFTarget('intensity', 'shelxl', 'shelxl')
    nan = numpy.full(n, numpy.nan)
    f = target_observation_arrays(hkl, fsq, sig, cell, sg, t,
                                  _synthetic_aux(n, numpy.full(n, 50.0), numpy.full(n, 500.0)))['flags']
    assert f['p_source'] == 'depositor' and f['extinction_applied'] and not f['reasons'], f
    f = target_observation_arrays(hkl, fsq, sig, cell, sg, t, _synthetic_aux(n, nan, nan))['flags']
    assert f['p_source'] == 'model' and f['n_weight_fallback'] == 0, f       # model Fc^2 for P
    assert f['reasons'], f                                                   # but no Fc for y
    part = numpy.full(n, 50.0)
    part[:3] = numpy.nan
    f = target_observation_arrays(hkl, fsq, sig, cell, sg, t, _synthetic_aux(n, part, part))['flags']
    assert f['p_source'] == 'mixed', f
    aux = _synthetic_aux(n, part, part, form='other', ext_kind='other')
    arr = target_observation_arrays(hkl, fsq, sig, cell, sg, t, aux)
    assert not arr['flags']['weights_applied'] and not arr['flags']['extinction_applied']
    assert len(arr['flags']['reasons']) == 2, arr['flags']
    plain = target_observation_arrays(hkl, fsq, sig, cell, sg, SFTarget('intensity'))
    assert numpy.array_equal(arr['sig'], plain['sig']) and numpy.all(arr['y'] == 1)
    try:
        target_observation_arrays(hkl, fsq, sig, cell, sg, t, None)
    except ValueError:
        pass
    else:
        raise AssertionError('a correction without aux must refuse')


def _shelx_hkl_file(path):
    '''The raw HKLF 4 rows (h k l I sigma, 3I4 2F8.2) a CIF embeds as _shelx_hkl_file.'''
    import numpy
    from chimerax.mmcif import get_cif_tables
    text = get_cif_tables(path, ['shelx'])[0].fields(('hkl_file',))[0][0]
    hkl, i, s = [], [], []
    for line in text.splitlines():
        try:
            h = [int(line[0:4]), int(line[4:8]), int(line[8:12])]
            row = float(line[12:20]), float(line[20:28])
        except ValueError:
            continue
        if h == [0, 0, 0]:
            break
        hkl.append(h); i.append(row[0]); s.append(row[1])
    return numpy.array(hkl, numpy.int32), numpy.array(i), numpy.array(s)


def test_extinction_already_in_fcf_data(session):
    '''A SHELXL .fcf's Fo^2 carry the refined extinction correction (Fo^2 = I/y): against
    the raw intensities the same CIF embeds, Fo^2/I tracks 1/y(h) as Clipper computes it
    (slope ~1, and the right size on the strongest reflections) - which also validates
    shelxl_extinction_factor's x / lambda^3 / sin(2 theta) / absolute-|Fc|^2 convention.
    So the 'shelxl' extinction target adds nothing on such data.'''
    import numpy
    from chimerax.clipper.symmetry import crystal_symmetry_from_cif_file
    from chimerax.clipper.io.small_molecule import (_parse_reflection_file,
                                                    merge_equivalent_reflections)
    from chimerax.clipper.io.sf_target import (SFTarget, small_molecule_target_aux,
        target_observation_arrays, shelxl_extinction_factor, _lookup_by_orbit)
    path = os.path.join(_DATA, 'cod_2021365.cif')
    cell, sg, grid = crystal_symmetry_from_cif_file(path)
    hkl, fsq, sig = _parse_reflection_file(path)
    aux = small_molecule_target_aux(path)
    assert aux['extinction'] == dict(aux['extinction'], kind='shelxl', x=0.022)
    assert aux['shelx_refln_list_code'] == 4
    plain = target_observation_arrays(hkl, fsq, sig, cell, sg, SFTarget('intensity'))
    ext = target_observation_arrays(hkl, fsq, sig, cell, sg,
                                    SFTarget('intensity', extinction='shelxl'), aux)
    assert ext['flags']['extinction_in_data'] and not ext['flags']['extinction_applied']
    assert numpy.array_equal(ext['i'], plain['i']) and numpy.array_equal(ext['sig'], plain['sig'])

    y = shelxl_extinction_factor(plain['hkl'], cell, 0.022, aux['wavelength'],
                                 _lookup_by_orbit(plain['hkl'], hkl, aux['fc_sq_model_abs'], sg))
    assert 0.7 < numpy.nanmin(y) < 0.8, numpy.nanmin(y)
    rh, ri, rs = _shelx_hkl_file(path)
    mh, mi, ms, _ = merge_equivalent_reflections(rh, ri, rs, sg)
    raw = _lookup_by_orbit(plain['hkl'], mh, mi, sg)
    raw_sig = _lookup_by_orbit(plain['hkl'], mh, ms, sg)
    ok = numpy.isfinite(raw) & (raw > 20 * raw_sig) & (plain['i'] > 0) & numpy.isfinite(y)
    log_ratio = numpy.log(plain['i'][ok] / raw[ok])
    log_y = numpy.log(y[ok])
    log_ratio -= numpy.median(log_ratio[log_y > numpy.median(log_y)])   # overall scale
    slope = numpy.polyfit(-log_y, log_ratio, 1)[0]
    assert 0.9 < slope < 1.1, slope
    strong = log_y < numpy.percentile(log_y, 3)
    assert abs(log_ratio[strong].mean() + log_y[strong].mean()) < 0.02, (
        log_ratio[strong].mean(), -log_y[strong].mean())


def _scaffold_state(path, obs):
    import numpy
    from chimerax.clipper.symmetry import crystal_symmetry_from_cif_file
    from chimerax.clipper.io.small_molecule import sfcalc_scaffold
    from chimerax.clipper.diff.state import XrayTargetState
    cell, sg, grid = crystal_symmetry_from_cif_file(path)
    sc = sfcalc_scaffold(path, cell, sg, grid)
    is_aniso = numpy.isfinite(sc['u_aniso'][:, 0]).astype(numpy.uint8)
    u_aniso = numpy.nan_to_num(sc['u_aniso'])
    state = XrayTargetState(sc['elements'], fobs=obs.fobs, phi_fom=obs.phi_fom,
                            usage=obs.usage, kind=obs.kind, iobs=obs.iobs)
    return state, (sc['coords'], sc['u_iso'], u_aniso, sc['occupancies'], is_aniso)


def test_target_variants_build(session):
    '''Every one of the 16 target variants builds on a real crystal and gives a finite
    value and gradient at the deposited model; the intensity target keeps the I <= 0
    reflections (with their own sigma) that the masked one leaves out. (This crystal has
    no refined extinction, so each +extinction variant equals its -extinction twin.)'''
    import numpy
    from chimerax.clipper.symmetry import crystal_symmetry_from_cif_file
    from chimerax.clipper.io.small_molecule import _parse_reflection_file
    from chimerax.clipper.io.sf_target import small_molecule_target_aux, observations_from_arrays
    path = os.path.join(_DATA, 'cod_1100908.cif')
    cell, sg, grid = crystal_symmetry_from_cif_file(path)
    hkl, fsq, sig = _parse_reflection_file(path)
    aux = small_molecule_target_aux(path)
    assert numpy.isfinite(aux['fc_sq_model_abs']).all()
    values = {}
    for t in _all_targets():
        obs = observations_from_arrays(hkl, fsq, sig, cell, sg, t, aux)
        assert obs.kind == t.kind and (obs.iobs is not None) == (t.observations == 'intensity')
        state, params = _scaffold_state(path, obs)
        L, g = state.value_and_gradient(*params, refresh_scale=True)
        assert numpy.isfinite(L) and L > 0 and numpy.isfinite(g).all(), (t, L)
        values[(t.observations, t.weights, t.extinction)] = L
        if t.observations == 'intensity':
            i = obs.arrays['i']
            n_le0 = int((i <= 0).sum())
            assert n_le0 > 0, 'fixture has no I <= 0 reflections'
            got = obs.iobs.data[1]
            assert int(numpy.isfinite(got[:, 0]).sum()) == len(i)
            assert (got[numpy.isfinite(got[:, 0]), 1] > 0).all()
            # the I <= 0 rows carry signal: hiding them changes the loss
            from chimerax.clipper import HKL_data_I_sigI
            keep = i > 0
            hidden = HKL_data_I_sigI(obs.hkl_info)
            hidden.set_data(obs.arrays['hkl'][keep], numpy.stack(
                [i[keep], obs.arrays['sig'][keep]], axis=1).astype(numpy.float32))
            obs.iobs = hidden
            state2, _ = _scaffold_state(path, obs)
            L2, _ = state2.value_and_gradient(*params, refresh_scale=True)
            assert L2 < L, (t, L, L2)
        if t.extinction == 'shelxl':
            assert not obs.flags['extinction_applied'] and obs.flags['extinction_kind'] == 'none'
            assert L == values[(t.observations, t.weights, 'none')], t
    assert len(set(values.values())) == 8, values


def test_sf_target_turnkey_matches_cached(session):
    '''small_molecule_ensemble_target (aux read from the files) and assembled_target_from_box
    (aux passed in, as from a cache) build the identical non-default target.'''
    import numpy
    from chimerax.clipper.symmetry import crystal_symmetry_from_cif_file
    from chimerax.clipper.io.small_molecule import _parse_reflection_file
    from chimerax.clipper.io.sf_target import SFTarget, small_molecule_target_aux
    from chimerax.clipper.diff.crystal import small_molecule_ensemble_target
    from chimerax.clipper.diff.assembled_cache import stamp_box, assembled_target_from_box
    path = os.path.join(_DATA, 'cod_2213867.cif')
    cell, sg, grid = crystal_symmetry_from_cif_file(path)
    hkl, fsq, sig = _parse_reflection_file(path)
    t = SFTarget('intensity', 'shelxl', 'shelxl')
    state, box = small_molecule_ensemble_target(session, path, sf_target=t,
                                                merge_equivalents=True)
    try:
        stamp_box(session, box, cell, sg, 0.8)
        cached, _ = assembled_target_from_box(session, box, hkl, fsq, sig, sf_target=t,
                                              merge_equivalents=True,
                                              aux=small_molecule_target_aux(path))
        c = numpy.array(box.atoms.coords)
        l0, g0 = state.value_and_gradient(c, refresh_scale=True)
        l1, g1 = cached.value_and_gradient(c, refresh_scale=True)
        assert l0 == l1 and numpy.array_equal(g0, g1), (l0, l1)
        assert state.sf_observations.flags == cached.sf_observations.flags
        try:
            small_molecule_ensemble_target(session, path, sf_target=t)
        except ValueError:
            pass
        else:
            raise AssertionError('a non-default target without merging must refuse')
    finally:
        from chimerax.core.commands import run
        run(session, 'close')


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
