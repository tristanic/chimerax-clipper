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
    '''The live-map compute engine (SmallMoleculeXmapMgr): the 2mFo-DFc map must
    place atoms on positive density, and the anisotropic+spline scaling must give an
    R-work close to the published value (which a crude overall scale would not - it
    left a large heavy-atom residual). Exercised headlessly; the GUI display itself
    needs an OpenGL context.'''
    import numpy
    from chimerax.clipper.symmetry import crystal_symmetry_from_cif_file
    from chimerax.clipper.io.small_molecule import (open_small_molecule_cif,
        _clipper_frame_coords, _small_molecule_map_data)
    from chimerax.clipper.maps.small_molecule_map import SmallMoleculeXmapMgr
    from chimerax.clipper.clipper_python import Coord_orth, Map_stats
    cif = os.path.join(_DATA, 'cod_1100908.cif')   # C2/c, Cu on a 2-fold; published R 0.041
    model = open_small_molecule_cif(session, cif)
    try:
        cell, sg, grid = crystal_symmetry_from_cif_file(cif)
        model.atoms.coords = _clipper_frame_coords(model, cif, cell)
        smd = _small_molecule_map_data(model, cif, None, cell, sg, grid)
        mgr = SmallMoleculeXmapMgr(smd['hklinfo'], smd['cell'], smd['spacegroup'],
            smd['grid'], smd['scaffold'], smd['fobs'], structure=model,
            scaffold_to_model=smd['scaffold_to_model'])
        mgr.add_xmap('2mFo-DFc', is_difference_map=False)
        mgr.add_xmap('mFo-DFc', is_difference_map=True)
        # R-work close to the published 0.041 (aniso+spline scaling).
        assert abs(mgr.rwork - 0.041) < 0.02, mgr.rwork
        # 2mFo-DFc places atoms on positive density.
        xm = mgr.get_xmap_ref('2mFo-DFc')
        st = Map_stats(xm)
        vals = numpy.array([xm.get_data(
            Coord_orth(float(x[0]), float(x[1]), float(x[2])).coord_frac(cell).coord_grid(
                xm.grid_sampling)) for x in model.atoms.coords])
        assert (vals > st.mean).mean() > 0.9, (vals > st.mean).mean()
        assert (vals.mean() - st.mean) / st.std_dev > 2.0
    finally:
        session.models.close([model])


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
