"""Pure Python association/lifetime regressions; no ChimeraX/OpenGL required.

Run with python3 -m unittest discover -s tests -p test_shared_grid.py.
Native rendering, masking and session serialization still need the GUI checks
in tests/shared_grid_manual.md. The doubles below model only APIs used here.
"""
import importlib.util
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch


class Triggers:
    def add_trigger(self, name):
        pass

    def add_handler(self, name, cb):
        return cb

    def remove_handler(self, cb):
        pass

    def activate_trigger(self, name, data):
        pass


class Drawing:
    def delete(self):
        pass


class Model:
    def __init__(self, name, session):
        self.name, self.session = name, session
        self._name = name
        self._deleted = False
        self.children = []
        self.triggers = Triggers()

    def add(self, models):
        self.children.extend(models)
        self.session.models.open.extend(models)

    def child_models(self):
        return self.children

    def delete(self):
        # Exact upstream Model.delete flow (786694db5), including re-entry via
        # Models.close before _deleted is set; Drawing has no native resources.
        if self._deleted:
            raise RuntimeError('Model %s was deleted twice' % self._name)
        models = self.session.models
        if models.have_model(self):
            models.close([self])
            return
        self._deleted = True
        Drawing.delete(self)
        self.triggers.activate_trigger('deleted', self)
        self.session = None

    @property
    def deleted(self):
        return self._deleted


class Volume(Model):
    def __init__(self, session, data):
        super().__init__('density', session)
        self.data = data
        data.add_change_callback(self.data_changed_cb)
        self.scene_position = SimpleNamespace(is_identity=lambda: True)
        self.surfaces = []
        self.region = ((0, 0, 0), (99, 99, 99), (1, 1, 1))
        self.display = True
        self.update_requests = []

    def show(self, style=None, rendering_options=None, show=True):
        # Exact executable behavior of Volume.show from upstream ChimeraX
        # 786694db5 (comments/docstring omitted). This MUST keep its eager read:
        # a no-op double hid a full-grid allocation during shared construction.
        if style is not None:
            self.set_display_style(style)
        if rendering_options:
            self.rendering_options = rendering_options
        self._keep_displayed_data = self.displayed_matrices() if show else None
        if show:
            self._drawings_need_update()
        self.display = show

    def _drawings_need_update(self):
        self.update_requests.append(self.region)

    def displayed_matrices(self):
        return [self.data.matrix(region=self.region)]

    def bounding_region(self, coords, **kw):
        return ((1, 1, 1), (3, 3, 3), (1, 1, 1))

    def new_region(self, ijk_min, ijk_max, ijk_step, **kw):
        self.region = (ijk_min, ijk_max, ijk_step)
        self._drawings_need_update()

    def data_changed_cb(self, reason):
        pass

    def set_parameters(self, **kw):
        if 'surface_levels' in kw:
            self.surfaces = [SimpleNamespace(level=x) for x in kw['surface_levels']]

    def set_display_style(self, style):
        self.style = style

    def copy_settings_from(self, source, **kw):
        assert kw == dict(copy_region=False, copy_xform=False, copy_zone=False)
        self.state_at_appearance_copy = (self.display, self.region)
        self.set_parameters(surface_levels=[s.level for s in source.surfaces])
        self.style = source.style

    def delete(self):
        # Exact upstream Volume.delete flow (786694db5).
        d = self.data
        if d:
            d.clear_cache()
            d.remove_change_callback(self.data_changed_cb)
        self.data = None
        self._keep_displayed_data = None
        self.outline_box = None
        self.close_models()
        Model.delete(self)

    def close_models(self):
        pass

    def replace_data(self, data):
        self.data.remove_change_callback(self.data_changed_cb)
        self.data = data
        data.add_change_callback(self.data_changed_cb)


class Grid:
    def __init__(self):
        self.callbacks = []
        self.reads = self.clears = 0
        self.allocations = []
        self.cache = {}
        self.xyz_to_ijk_transform = object()

    def matrix(self, region=None):
        self.reads += 1
        # Emulate file-backed data: cache eviction can allocate a NEW buffer on
        # a later request while another handler still pins the original buffer.
        region = region or ((0, 0, 0), (99, 99, 99), (1, 1, 1))
        if region not in self.cache:
            self.cache[region] = object()
            self.allocations.append(region)
        return self.cache[region]

    def add_change_callback(self, cb):
        self.callbacks.append(cb)

    def remove_change_callback(self, cb):
        self.callbacks.remove(cb)

    def clear_cache(self):
        self.clears += 1
        self.cache.clear()


def module(name, **attrs):
    m = ModuleType(name)
    m.__dict__.update(attrs)
    return m


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


class AtomSpecArg:
    @classmethod
    def parse(cls, text, session):
        spec = SimpleNamespace(evaluate=lambda session:
            SimpleNamespace(models=session.selection))
        return spec, text, ''


class SharedGridTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modules = patch.dict(sys.modules, {
            'numpy': module('numpy', array=lambda x: x),
            'chimerax': module('chimerax'),
            'chimerax.core': module('chimerax.core'),
            'chimerax.core.models': module('models', Model=Model),
            'chimerax.core.errors': module('errors', UserError=ValueError),
            'chimerax.core.commands': module('commands'),
            'chimerax.core.commands.atomspec': module('atomspec', AtomSpecArg=AtomSpecArg),
            'chimerax.map': module('map', Volume=Volume),
            'chimerax.map.volume': module('volume', VolumeSurface=object,
                mean_sd_rms=lambda matrix: (2, 1, 2.25)),
            'chimerax.clipper': module('chimerax.clipper'),
            'chimerax.clipper.maps': module('chimerax.clipper.maps'),
        })
        cls.modules.start()
        src = Path(__file__).resolve().parents[1] / 'src'
        load('chimerax.clipper.maps.mapset_base', src / 'maps/mapset_base.py')
        cls.base = load('chimerax.clipper.maps.map_handler_base',
            src / 'maps/map_handler_base.py')
        cls.nx = load('chimerax.clipper.maps.nxmapset', src / 'maps/nxmapset.py')
        cls.cmd = load('chimerax.clipper.cmd', src / 'cmd.py')

    @classmethod
    def tearDownClass(cls):
        cls.modules.stop()

    def setUp(self):
        models = SimpleNamespace(open=[], list=lambda: models.open)
        models.have_model = lambda model: model in models.open
        def close(volumes):
            # Models.close removes registrations BEFORE invoking delete. This
            # makes direct registered-model deletion re-enter the handler.
            for v in volumes:
                if v in models.open:
                    models.open.remove(v)
            for v in volumes:
                if not v._deleted:
                    v.delete()
        models.close = close
        self.session = SimpleNamespace(models=models)
        self.grid = Grid()

    def mapset(self):
        crystal = SimpleNamespace(normalize_scene_positions=lambda: None)
        mgr = Model('manager', self.session)
        mgr.crystal_mgr = crystal
        mgr.spotlight_mode = False
        mgr.rezone_needed = lambda: None
        return self.nx.NXmapSet(mgr, 'maps')

    def source(self):
        maps = self.mapset()
        h = self.nx.NXmapHandler(maps, data=self.grid, name='shared',
            is_difference_map=True)
        h.set_parameters(surface_levels=[-2, 2])
        h.set_display_style('mesh')
        maps.add([h])
        return h

    def test_share_retains_source_without_reading_grid_and_is_idempotent(self):
        source = self.source()
        target = self.mapset()
        reads = self.grid.reads
        copied = target.add_nxmap_handler_from_volume(source)
        self.assertIs(copied.data, source.data)
        self.assertFalse(source.deleted)
        self.assertTrue(copied.is_difference_map)
        self.assertEqual(self.grid.reads, reads)
        self.assertEqual([s.level for s in copied.surfaces], [-2, 2])
        copied.surfaces[0].level = -4
        self.assertEqual(source.surfaces[0].level, -2)
        self.assertIs(target.add_nxmap_handler_from_volume(source), copied)
        self.assertEqual(len(target.all_maps), 1)

    def test_uncached_file_grid_does_not_eagerly_reread_full_map(self):
        source = self.source()
        sibling = self.mapset().add_nxmap_handler_from_volume(source)
        pinned = source._keep_displayed_data[0]
        sibling.delete()  # Volume.delete evicts shared grid cache.
        self.assertEqual(self.grid.cache, {})
        reads, allocations = self.grid.reads, len(self.grid.allocations)
        target = self.mapset()
        target._mgr.spotlight_mode = True
        target._mgr.spotlight_radius = 2
        target._mgr.box_center = (2, 2, 2)
        with patch.object(self.nx, '_find_box_corners',
                return_value=((1, 1, 1), (3, 3, 3))):
            copied = target.add_nxmap_handler_from_volume(source)
        self.assertEqual(self.grid.reads, reads)
        self.assertEqual(len(self.grid.allocations), allocations)
        self.assertIs(source._keep_displayed_data[0], pinned)
        spotlight = ((1, 1, 1), (3, 3, 3), (1, 1, 1))
        self.assertEqual(copied.state_at_appearance_copy, (False, spotlight))
        self.assertTrue(copied.display)
        copied.displayed_matrices()  # Eventual redraw may read the intended region.
        self.assertEqual(self.grid.allocations[allocations:], [spotlight])

    def test_close_either_order_unregisters_only_own_callbacks(self):
        for reverse in (False, True):
            source = self.source()
            other = self.mapset().add_nxmap_handler_from_volume(source)
            first, last = (other, source) if reverse else (source, other)
            first.delete()
            first.delete()  # A second deletion must not remove sibling callbacks.
            self.assertNotIn(first._update_mean_sd_rms_cb, self.grid.callbacks)
            self.assertIn(last._update_mean_sd_rms_cb, self.grid.callbacks)
            self.assertIs(last.data, self.grid)
            for cb in list(self.grid.callbacks):
                cb('values changed')
            last.delete()
            self.assertEqual(self.grid.callbacks, [])

    def test_direct_registered_delete_reenters_before_deleted_and_preserves_sibling(self):
        source = self.source()
        sibling = self.mapset().add_nxmap_handler_from_volume(source)
        callbacks_before = len(self.grid.callbacks)
        entries = []
        original_delete = self.base.MapHandlerBase.delete
        def record_delete(handler):
            if handler is source:
                entries.append((handler.data, handler.deleted))
            return original_delete(handler)
        with patch.object(self.base.MapHandlerBase, 'delete', record_delete):
            source.delete()
        self.assertEqual(entries, [(self.grid, False), (None, False)])
        self.assertTrue(source.deleted)
        self.assertNotIn(source, self.session.models.open)
        self.assertIsNone(source.data)
        self.assertIs(sibling.data, self.grid)
        self.assertFalse(sibling.deleted)
        self.assertEqual(len(self.grid.callbacks), callbacks_before - 2)
        self.assertIn(sibling._update_mean_sd_rms_cb, self.grid.callbacks)
        sibling.displayed_matrices()
        for cb in list(self.grid.callbacks):
            cb('values changed')
        self.session.models.close([sibling])
        self.assertTrue(sibling.deleted)
        self.assertEqual(self.grid.callbacks, [])

    def test_trigger_cleanup_failure_still_unregisters_grid_callback(self):
        source = self.source()
        def missing_handler(cb):
            raise ValueError('manager already closed')
        source.mapset.triggers.remove_handler = missing_handler
        source.delete()
        self.assertEqual(self.grid.callbacks, [])

    def test_replacing_grid_detaches_statistics_callback_and_closes_cleanly(self):
        source = self.source()
        other = self.mapset().add_nxmap_handler_from_volume(source)
        new_grid = Grid()
        other.replace_data(new_grid)
        self.assertNotIn(other._update_mean_sd_rms_cb, self.grid.callbacks)
        self.assertIn(other._update_mean_sd_rms_cb, new_grid.callbacks)
        self.assertIn(source._update_mean_sd_rms_cb, self.grid.callbacks)
        other.delete()
        source.delete()
        self.assertEqual(self.grid.callbacks, [])
        self.assertEqual(new_grid.callbacks, [])

    def test_fresh_volume_keeps_replacement_behavior_and_difference_flag(self):
        source = Volume(self.session, self.grid)
        self.session.models.open.append(source)
        target = self.mapset()
        h = target.add_nxmap_handler_from_volume(source,
            is_difference_map=True, contour=[-2, 2])
        self.assertTrue(source.deleted)
        self.assertIs(h.data, self.grid)
        self.assertTrue(h.is_difference_map)
        self.assertEqual([s.level for s in h.surfaces], [-2, 2])

    def test_command_prevalidates_all_sources_before_creating_manager(self):
        source = self.source()
        transformed = Volume(self.session, self.grid)
        transformed.scene_position = SimpleNamespace(is_identity=lambda: False)
        called = []
        symmetry = module('symmetry', get_map_mgr=lambda *a, **kw: called.append(a))
        with patch.dict(sys.modules, {'chimerax.clipper.symmetry': symmetry}):
            with self.assertRaisesRegex(ValueError, 'absolute reference frame'):
                self.cmd.associate_volumes(self.session, [source, transformed], object())
        self.assertEqual(called, [])
        self.assertFalse(source.deleted)

    def test_rejects_crystallographic_and_transformed_sources(self):
        class CrystalVolume(Volume):
            pass
        with self.assertRaisesRegex(ValueError, 'Crystallographic'):
            self.mapset().add_nxmap_handler_from_volume(CrystalVolume(self.session, self.grid))
        source = self.source()
        source.scene_position = SimpleNamespace(is_identity=lambda: False)
        with self.assertRaisesRegex(ValueError, 'absolute reference frame'):
            self.mapset().add_nxmap_handler_from_volume(source)
        self.assertFalse(source.deleted)

    def test_command_parser_includes_managed_real_space_maps(self):
        source = self.source()
        ordinary = Volume(self.session, self.grid)
        class CrystalVolume(Volume):
            pass
        crystal = CrystalVolume(self.session, self.grid)
        self.session.selection = [source, ordinary, crystal, object()]
        volumes, _, _ = self.cmd.VolumesArg.parse('#1', self.session)
        self.assertEqual(volumes, [source, ordinary, crystal])

    def test_parser_to_command_rejects_empty_crystal_and_mixed_without_mutation(self):
        source = self.source()
        class CrystalVolume(Volume):
            pass
        crystal = CrystalVolume(self.session, self.grid)
        called = []
        symmetry = module('symmetry', get_map_mgr=lambda *a, **kw: called.append(a))
        for selection in ([], [object()], [crystal], [source, crystal]):
            with self.subTest(selection=selection):
                self.session.selection = selection
                volumes, _, _ = self.cmd.VolumesArg.parse('#1', self.session)
                with patch.dict(sys.modules, {'chimerax.clipper.symmetry': symmetry}):
                    with self.assertRaises(ValueError):
                        self.cmd.associate_volumes(self.session, volumes, object())
                self.assertEqual(called, [])
                self.assertFalse(source.deleted)
                self.assertFalse(crystal.deleted)


if __name__ == '__main__':
    unittest.main()
