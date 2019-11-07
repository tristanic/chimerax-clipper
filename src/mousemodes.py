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

import numpy
from chimerax.mouse_modes import MouseMode, ZoomMouseMode as ZoomMouseMode_Base

def initialize_clipper_mouse_modes(session):
    initialize_zoom_mouse_modes(session)
    initialize_map_contour_mouse_modes(session)
    std_modes = session.ui.mouse_modes.modes
    from chimerax.mouse_modes.std_modes import TranslateMouseMode
    tmml = list(filter(lambda m: type(m)==TranslateMouseMode, std_modes))
    if len(tmml) != 1:
        # Just in case. There *should* only be one TranslateMouseMode in the list,
        # but if not let's play it safe and just create a fresh one
        tmm = TranslateMouseMode(session)
    else:
        tmm = tmml[0]
    session.ui.mouse_modes.bind_mouse_mode('left',['shift'], tmm)

def initialize_zoom_mouse_modes(session):
    z = ZoomMouseMode(session)
    c = ClipPlaneAdjuster(session, z)
    s = Z_Shift_CofR(session)
    session.ui.mouse_modes.bind_mouse_mode('right',['shift'], z) # Legacy mode
    session.ui.mouse_modes.bind_mouse_mode('wheel',[], z)
    session.ui.mouse_modes.bind_mouse_mode('wheel',['shift'], c)
    session.ui.mouse_modes.bind_mouse_mode('right',['control'], s)

def initialize_map_contour_mouse_modes(session):
    #z = ZoomMouseMode(session)
    s = SelectVolumeToContour(session)
    v = ContourSelectedVolume(session, s, True)
    session.ui.mouse_modes.bind_mouse_mode('wheel',['control'], s)
    # session.ui.mouse_modes.bind_mouse_mode('wheel',[], v)
    session.ui.mouse_modes.bind_mouse_mode('wheel',['alt'], v)

class ClipPlaneAdjuster(MouseMode):
    def __init__(self, session, zoom_mode):
        super().__init__(session)
        self._zoomer = zoom_mode

    def cleanup(self):
        pass

    def wheel(self, event):
        mult = 1-event.wheel_value()/30
        z = self._zoomer
        z.far_clip_multiplier *= mult
        z.near_clip_multiplier *= mult

class Z_Shift_CofR(MouseMode):
    def __init__(self, session):
        self._step_multiplier = 1
        super().__init__(session)
        from chimerax.core.graphics import Drawing
        d = self._drawing = Drawing('Depth Indicator')
        d.set_geometry(*self._drawing_geometry())
        d.display = False
        self.view.drawing.add_drawing(d)

    def cleanup(self):
        self.view.drawing.remove_drawing(self._drawing)

    def mouse_down(self, event):
        MouseMode.mouse_down(self, event)
        self._set_drawing_position()
        self._set_drawing_color()
        self._drawing.display=True

    def mouse_drag(self, event):
        dx, dy = self.mouse_motion(event)
        self.move_camera_and_cofr(-3*self.pixel_size()*dy)
        self._set_drawing_position()

    def mouse_up(self, event):
        self._drawing.display=False
        MouseMode.mouse_up(self, event)

    def wheel(self, event):
        d = event.wheel_value()
        psize = self.pixel_size()
        self.move_camera_and_cofr(-100*d*psize)

    def move_camera_and_cofr(self, dz):
        cofr = self.view.center_of_rotation
        cofr_method = self.view.center_of_rotation_method
        camera = self.view.camera
        cpos = self.camera_position.origin()
        vd = camera.view_direction()
        import numpy
        cc = cofr-cpos
        shift_vec = numpy.dot(cc/numpy.linalg.norm(cc), vd) *vd * dz
        from chimerax.core.geometry import translation
        t = translation(shift_vec)
        self.view.center_of_rotation += shift_vec
        self.view.center_of_rotation_method = cofr_method
        camera.set_position(t*camera.position)

    def _drawing_geometry(self):
        from chimerax.surface.shapes import box_geometry
        return box_geometry((-1,-1,-0.01), (1,1,0.01))

    def _set_drawing_position(self):
        from chimerax.core.geometry import Place, scale
        scale = 300*self.pixel_size()
        p = Place(axes=self.view.camera.position.axes()*scale, origin=self.view.center_of_rotation)
        self._drawing.position = p

    def _set_drawing_color(self):
        from chimerax.core.colors import contrast_with, Color
        import numpy
        color = numpy.array([0,0,0,0.5], numpy.float32)
        color[:3] = contrast_with(self.view.background_color)
        self._drawing.color = Color(color).uint8x4()



class ZoomMouseMode(ZoomMouseMode_Base):
    def __init__(self, session):
        super().__init__(session)
        self._far_clip_multiplier = 0.5
        self._near_clip_multiplier = 0.5

    @property
    def far_clip_multiplier(self):
        '''
        Multiplier applied to the camera-cofr distance to decide the position
        of the rear clipping plane. Clamped to the range (0..1)
        '''
        return self._far_clip_multiplier

    @far_clip_multiplier.setter
    def far_clip_multiplier(self, val):
        val = max(min(val, 1), 0)
        self._far_clip_multiplier = val
        # Zoom with a delta-z of zero to force redraw
        self.zoom(0)

    @property
    def near_clip_multiplier(self):
        '''
        Multiplier applied to the camera-cofr distance to decide the position
        of the near clipping plane. Clamped to the range (0..1)
        '''
        return self._near_clip_multiplier

    @near_clip_multiplier.setter
    def near_clip_multiplier(self, val):
        val = max(min(val, 1), 0)
        self._near_clip_multiplier = val
        # Zoom with a delta-z of zero to force redraw
        self.zoom(0)

    @property
    def far_clip(self):
        v = self.session.main_view
        cp = v.clip_planes
        fc = cp.find_plane('far')
        if fc is not None:
            return fc
        else:
            from chimerax.core.graphics.clipping import CameraClipPlane
            c = v.camera
            cofr = v.center_of_rotation
            cpos = c.position.origin()
            clip_point = cpos + (1+self.far_clip_multiplier)* (cofr-cpos)
            camera_point = c.position.inverse() * clip_point
            fc = CameraClipPlane('far', (0,0,1),
                                camera_point, v)
            # Put the near clip at the camera position for depth cueing
            v.clip_planes.add_plane(fc)
            return fc

    @property
    def near_clip(self):
        v = self.session.main_view
        cp = v.clip_planes
        nc = cp.find_plane('near')
        if nc is not None:
            return nc
        else:
            from chimerax.core.graphics.clipping import CameraClipPlane
            c = v.camera
            cofr = v.center_of_rotation
            cpos = c.position.origin()
            clip_point = cpos + (1-self.near_clip_multiplier) * (cofr-cpos)
            camera_point = c.position.inverse() * clip_point
            nc = CameraClipPlane('near', (0,0,-1),
                                camera_point, v)
            # Put the near clip at the camera position for depth cueing
            v.clip_planes.add_plane(nc)
            return nc

    def cleanup(self):
        pass

    def zoom(self, delta_z, stereo_scaling = False):
        v = self.view
        c = v.camera
        cofr = v.center_of_rotation
        if stereo_scaling and c.name == 'stereo':
            v.stereo_scaling(delta_z)
        if c.name == 'orthographic':
            import numpy
            c.field_width = max(c.field_width - delta_z, self.pixel_size())
            # TODO: Make camera field_width a property so it knows to redraw.
            from chimerax.core.geometry import place
            camera_to_cofr = cofr - c.position.origin()
            vd = c.view_direction()
            current_forward_distance = numpy.dot(camera_to_cofr, vd)
            new_view_distance = c.field_width * 0.7
            shift = current_forward_distance - new_view_distance
            new_origin = c.position.origin() + shift*vd
            new_pos = place.Place(axes = c.position.axes(), origin = new_origin)
            c.position = new_pos
            #distance = cofr-new_origin
            new_view_vec = new_view_distance*vd
            self.far_clip.plane_point = new_origin + new_view_vec*(1+self.far_clip_multiplier)
            self.near_clip.plane_point = new_origin + new_view_vec*(1-self.near_clip_multiplier)
            c.redraw_needed = True
        else:
            shift = c.position.transform_vectors((0, 0, delta_z))
            v.translate(shift)
            new_origin = c.position.origin()
            self.far_clip.plane_point = new_origin + (cofr-new_origin)*2

class SelectVolumeToContour(MouseMode):
    '''
    Designed to work together with ContourSelectedVolume.
    Each step of the mouse wheel increments through the currently
    loaded Volume objects, temporarily selecting the current pick to
    highlight it in the display. Stores a reference to the last selected
    volume, accessible via picked_volume. If the last selected volume
    has never been set or has been deleted, picked_volume defaults to
    the first Volume object in session.models.list().
    '''
    def __init__(self, session, cooldown_time = 0.15):
        super().__init__(session)
        self._last_picked_index = 0
        self._picked_volume = None
        self._deselect_handler = None
        from time import time
        self._last_event_time = time()
        self._cooldown_time = cooldown_time
        self._time_until_deselect = 5 # seconds
        self._deselect_start_time = 0
    def wheel(self, event):
        '''Select the next visible volume.'''
        # Simply acting on every event works well on a mouse, but horribly on
        # a trackpad - presumably because swiping on the trackpad generates a
        # much longer string of events with finer-grained fractional values:
        # either nothing happens, or a large number happen all at once, making
        # it nearly impossible to control exactly which map ends up selected.
        # Instead, we'll change the index by only +/- 1 per event, and impose
        # a cooldown period during which further events will be ignored.
        from time import time
        if time() - self._last_event_time < self._cooldown_time:
            return
        self._last_event_time = time()
        from math import copysign
        d = int(copysign(1, event.wheel_value()))
        # d = int(event.wheel_value())
        vol_list = self._get_vol_list()
        for v in vol_list:
            v.selected = False
        n = len(vol_list)
        if n == 0:
            return
        last = self._last_picked_index
        p = (last + d) % n
        sv = self._picked_volume = vol_list[p]
        sv.selected = True
        self._last_picked_index = p
        self.session.logger.status('Selected for contouring: {}'.format(sv.name))
        self._start_deselect_timer()
    def _get_vol_list(self):
        '''Get a list of currently visible volumes.'''
        from chimerax.map import Volume
        vlist = self.session.models.list(type=Volume)
        for i in reversed(range(len(vlist))):
            if not vlist[i].visible:
                vlist.pop(i)
        return vlist


    @property
    def picked_volume(self):
        try:
            vol_list = self._get_vol_list()
            if not len(vol_list):
                return None
            pv = self._picked_volume
            if pv is None:
                pv = self._picked_volume = vol_list[self._last_picked_index]
            elif pv not in vol_list or not pv.visible:
                pv = self._picked_volume = vol_list[0]
                self._last_picked_index = 0
        except IndexError:
            pv = self._picked_volume = vol_list[0]
            self._last_picked_index = 0
        return self._picked_volume

    def _start_deselect_timer(self):
        from time import time
        self._deselect_start_time = time()
        if self._deselect_handler is None:
            self._deselect_handler = self.session.triggers.add_handler(\
                                'new frame', self._deselect_on_timeout)

    def _deselect_on_timeout(self, *_):
        from time import time
        if time()- self._deselect_start_time > self._time_until_deselect:
            self.session.triggers.remove_handler(self._deselect_handler)
            self._deselect_handler = None
            self._picked_volume.selected = False

class ContourSelectedVolume(MouseMode):
    def __init__(self, session, selector, symmetrical=True):
        '''
        Modified volume contouring method which acts on a single volume
        at a time. By default, changes all contours towards/away from
        zero. If the volume has a surface_zone property set, it will be
        automatically masked back down after contouring.
        Args:
            session:
                The ChimeraX session.
            selector:
                A SelectVolumeToContour object used to define the current
                target volume.
            symmetrical:
                If True, scrolling up will adjust contours away from
                zero (that is, negative contours will get more negative).
                If False, all contours will be shifted in the same
                direction.
        '''
        super().__init__(session)
        # SelectVolumeToContour object telling us which volume to work on
        self.selector = selector
        self.symmetrical = symmetrical
        self.target_volume = None


    def wheel(self, event):
        d = event.wheel_value()
        v = self.selector.picked_volume
        if v is not None:
            self.target_volume = v
            if hasattr(v, 'sigma'):
                sd = v.sigma
            else:
                sd = v.mean_sd_rms()[1]
            step = d/30 * sd
            rep, levels = adjust_threshold_level(v, step, self.symmetrical)
            if rep != 'image':
                lsig = tuple(l/sd for l in levels)
                level_strings = []
                for l in levels:
                    if l > 0.1:
                        level_strings.append('{:.3f}'.format(l))
                    else:
                        level_strings.append('{:1.2e}'.format(l))
                lstr = ', '.join(level_strings)
                sstr = ', '.join(format(s, '.3f') for s in lsig)
                self.session.logger.status('Volume {} contour level(s): {} ({} sigma)'.format(v.name, lstr, sstr))




def adjust_threshold_level(m, step, sym):
    if m.image_shown:
        ms = m.matrix_value_statistics()
        new_levels = [(l+step,b) for l,b in m.image_levels]
        l,b = new_levels[-1]
        new_levels[-1] = (max(l,1.01*ms.maximum),b)
        m.set_parameters(image_levels = new_levels)
        return ('image', new_levels)
    else:
        #if sym and len(m.surface_levels) > 1:
        if sym and len(m.surfaces) > 1:
            #old_levels = m.surface_levels
            old_levels = [s.level for s in m.surfaces]
            new_levels = []
            for l in old_levels:
                if l < 0:
                    newl = l-step
                    if newl > 0:
                        newl = -(abs(step))
                    new_levels.append(newl)
                else:
                    newl = l+step
                    if newl < 0:
                        newl = abs(step)
                    new_levels.append(newl)
        else:
            old_levels = [s.level for s in m.surfaces]
            #new_levels = tuple(l+step for l in m.surface_levels)
            new_levels = tuple(l+step for l in old_levels)
        m.set_parameters(surface_levels = new_levels)
    return (None, new_levels)
