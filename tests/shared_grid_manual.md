# Shared real-space maps: ChimeraX/ISOLDE verification

This change gives each Clipper system its own NXmapHandler, contours, display
region and zone/mask state, while sharing the exact GridData object. Grid values
and geometry (origin, step, rotation) remain shared. Editing those affects every
handler. It applies to ordinary real-space volumes and NXmapHandler sources;
crystallographic maps remain tied to their original system.

Constructing/associating a managed map makes no eager matrix request or density
array copy. Initial full-grid statistics are reused, avoiding the temporary
arrays required to recalculate them. The new handler stays hidden until its
target region and appearance are set; a subsequent redraw reads the intended
region as needed. Independent contour surfaces, region caches, masking and
ISOLDE fitting buffers can still use additional memory.

ChimeraX Volume.delete clears the grid's cache. Closing one shared handler can
therefore cause a later file read/cache refill in the other handler, but does
not destroy that handler's GridData. A survivor can still pin an old displayed
matrix after cache eviction. Redrawing a new handler can then allocate another
buffer, potentially the full grid if its display region covers the full map.
GridData identity does not guarantee one resident matrix under all cache and
display conditions. In-memory ArrayGridData retains its array.
The change deliberately uses ChimeraX's existing lifetime and session behavior;
it does not replace or monkeypatch cache eviction.

## GUI and fitting

Use a ChimeraX installation containing both modified bundles. The automated
Python tests use API doubles; they do not verify native rendering or simulation.

1. Open two copies of an atomic model and one real-space map. Associate the map
   with model A using ISOLDE's **From loaded volume** menu. Confirm the original
   standalone volume is replaced, as before.
2. Select model B. Open the menu and confirm A's managed map appears, labelled
   with A's model ID/name. Select it and check both systems retain their maps.
   The Non-crystallographic Map Settings controls should populate after the
   next frame. Enable MDFF for B, start/stop a short fitting simulation, then
   switch to A and repeat. Each should have a separate MDFF manager and controls.
3. In the Python shell, obtain the two NXmapHandlers (a and b) from
   `get_map_mgr(model).all_non_xtal_maps`, imported from
   `chimerax.clipper.symmetry`. Check `a is not b` and `a.data is b.data`.
   Set a breakpoint/spy on `grid.matrix` if verifying that constructing the new
   handler makes no eager matrix request. Include the scenario where a sibling
   was closed and evicted the cache while A still pins its displayed matrix.
   Check the first later display read uses B's intended spotlight region, rather
   than first requesting the entire map. Grid identity alone cannot establish
   the number of resident cached/displayed density arrays.
4. Change B's contour, color, surface style and spotlight/masking region. Confirm
   A's controls/region stay independent. Confirm a difference map keeps both
   signed contours, colors and its difference-map flag after reuse. Its MDFF
   setting should initially be disabled, following existing difference-map rules.
5. Reopen B's loaded-volume menu. Its grid must be absent. Associate A's map to B
   again from the command line (`clipper associate #<a> toModel #<modelB>`);
   confirm B still has exactly one handler for that grid.
6. Add a third system sharing the grid. On a fourth model's menu, confirm it is
   offered once, with the existing owners listed. Confirm live and static
   crystallographic map handlers are not offered.
7. Translate/rotate an ordinary volume, then try to associate it. Confirm the
   absolute-frame error appears before the volume closes or a system is created.
   Also test a transform inherited from a parent model.

## Closure and save/restore

1. Close A's map/system and continue contouring/fitting B. Repeat in a fresh
   session with B closed first. Check `grid.change_callbacks` contains no bound
   callbacks whose `__self__` is a closed handler. Change grid values and call
   `grid.values_changed()`; the survivor should update without stale callbacks.
2. Save a session containing A and B with their shared file-backed map. Close
   the session and restore it; verify independent appearance and fitting plus
   `a.data is b.data`. Repeat with maps included in the session.
3. Repeat using an in-memory `ArrayGridData` (no backing file), saving a session
   including maps. Verify the same grid identity and density values on restore,
   and survival of either closure order.

NXmapHandler snapshots already delegate GridData serialization to ChimeraX's
`Volume.take_snapshot` / `grid_data_state`, which deduplicates by grid identity.
No new session format or state version is needed. The round-trip and cache
behavior above must still be checked in the actual application.
