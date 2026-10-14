"""A newly associated real-space map must be contoured before moving the view.

Run in ChimeraX with Clipper installed (no OpenGL context is needed)::

    ChimeraX --nogui --notools --exit --script tests/test_initial_nxmap_surface.py

The old initialization converted XYZ to grid indices twice, so a translated
origin or non-unit voxel size could leave an empty region until the first view
translation. Exercise real grid transforms, contour threads and zone masks.
"""

import time

import numpy

from chimerax.atomic import AtomicStructure, Element
from chimerax.map_data import ArrayGridData
from chimerax.map.volume import volume_from_grid_data
from chimerax.clipper import mousemodes, symmetry


def settle(session):
    # Match GUI deferred drawing rather than the synchronous in-script shortcut.
    # A contour is computed off-thread, queued on one frame, then applied on the
    # following frame. Allow a bounded interval for all of those stages.
    for frame in range(100):
        session.triggers.activate_trigger('new frame', frame)
        session.triggers.activate_trigger('graphics update', None)
        time.sleep(0.01)


def visible_triangles(volume):
    return sum(
        len(s.triangles) if s.triangle_mask is None else int(s.triangle_mask.sum())
        for s in volume.surfaces
        if s.display and s.triangles is not None
    )


def check_map(session, origin, step, spotlight=True):
    z, y, x = numpy.indices((33, 33, 33))
    # Keep the density's physical radius constant as the voxel size changes.
    distance_squared = sum(((axis - 16) * spacing)**2
                           for axis, spacing in zip((x, y, z), step))
    array = numpy.exp(-distance_squared / 32).astype(numpy.float32)
    data = ArrayGridData(array, origin=origin, step=step, name='gaussian')
    center = numpy.array(data.ijk_to_xyz((16, 16, 16)))

    model = AtomicStructure(session, name='initial surface regression')
    residue = model.new_residue('ALA', 'A', 1)
    atom = model.new_atom('CA', Element.get_element('C'))
    residue.add_atom(atom)
    atom.coord = center
    for name, offset in [('N', -1), ('C', 1)]:
        atom = model.new_atom(name, Element.get_element(name))
        residue.add_atom(atom)
        atom.coord = center + offset
    session.models.add([model])
    session.main_view.center_of_rotation = center
    volume = volume_from_grid_data(data, session)
    manager = symmetry.get_symmetry_handler(model, create=True,
                                           auto_add_to_session=True)
    manager.spotlight_mode = spotlight
    handler = manager.map_mgr.nxmapset.add_nxmap_handler_from_volume(
        volume, contour=[0.2])
    try:
        initial_region = tuple(tuple(r) for r in handler.region)
        settle(session)
        initial_triangles = visible_triangles(handler)
        assert handler.display and handler.parents_displayed

        # This explicitly reproduces the action that used to make the mesh
        # appear; it must no longer be necessary for a valid initial surface.
        session.main_view.center_of_rotation = center + (1., 0., 0.)
        settle(session)
        moved_triangles = visible_triangles(handler)
        print(f'origin={origin}, step={step}, spotlight={spotlight}: '
              f'initial region={initial_region}, '
              f'visible triangles before/after={initial_triangles}/{moved_triangles}')
        assert initial_triangles > 0, 'No visible surface before translating the view'
        assert moved_triangles > 0, 'Surface disappeared after translating the view'
        if not spotlight:
            assert initial_region == ((0, 0, 0), (32, 32, 32), (1, 1, 1))
    finally:
        session.models.close([manager])


def main(session):
    previous_script = session.in_script
    previous_camera_reset = symmetry.auto_reset_camera
    previous_mouse_setup = mousemodes.initialize_clipper_mouse_modes
    session.in_script = False
    # Camera reset and mouse binding require a GUI/OpenGL context. The actual
    # spotlight, region, contour and mask paths remain enabled below.
    symmetry.auto_reset_camera = False
    mousemodes.initialize_clipper_mouse_modes = lambda session: None
    try:
        for origin, step, spotlight in (
                ((100., 100., 100.), (1., 1., 1.), True),
                ((0., 0., 0.), (2., 3., 4.), True),
                ((-100., 50., 200.), (2., 3., 4.), True),
                ((0., 0., 0.), (1., 1., 1.), True),
                ((100., 100., 100.), (1., 1., 1.), False)):
            check_map(session, origin, step, spotlight)
    finally:
        symmetry.auto_reset_camera = previous_camera_reset
        mousemodes.initialize_clipper_mouse_modes = previous_mouse_setup
        session.in_script = previous_script
    print('test_initial_nxmap_surface: PASS')


main(session)
