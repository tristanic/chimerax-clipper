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
        _clipper_frame_coords, hydrate_small_molecule_model, _small_molecule_map_data)
    from chimerax.clipper.maps.small_molecule_map import SmallMoleculeXmapMgr
    from chimerax.clipper.clipper_python import Coord_orth, Map_stats
    cif = os.path.join(_DATA, 'cod_1100908.cif')   # C2/c, Cu on a 2-fold; published R 0.041
    model = open_small_molecule_cif(session, cif)
    try:
        cell, sg, grid = crystal_symmetry_from_cif_file(cif)
        model.atoms.coords = _clipper_frame_coords(model, cif, cell)
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
        _clipper_frame_coords, hydrate_small_molecule_model, _small_molecule_map_data)
    from chimerax.clipper.maps.small_molecule_map import SmallMoleculeXmapMgr
    cif = os.path.join(_DATA, 'cod_1100908.cif')
    model = open_small_molecule_cif(session, cif)
    try:
        cell, sg, grid = crystal_symmetry_from_cif_file(cif)
        model.atoms.coords = _clipper_frame_coords(model, cif, cell)
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
