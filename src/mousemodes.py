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
from chimerax.mouse_modes import (
    MouseMode, 
    ZoomMouseMode as ZoomMouseMode_Base,
    RotateMouseMode as RotateMouseMode_Base,
    SelectContextMenuAction
)

def initialize_clipper_mouse_modes(session):
    initialize_zoom_mouse_modes(session)
    initialize_map_contour_mouse_modes(session)
    session.ui.mouse_modes.bind_mouse_mode('left', [], RotateMouseMode(session))

def initialize_zoom_mouse_modes(session):
    z = ZoomMouseMode(session)
    c = ClipPlaneAdjuster(session, z)
    s = Z_Shift_CofR(session)
    session.ui.mouse_modes.bind_mouse_mode('right',['shift'], z) # Legacy mode
    session.ui.mouse_modes.bind_mouse_mode('wheel',[], z)
    session.ui.mouse_modes.bind_mouse_mode('wheel',['shift'], c)
    session.ui.mouse_modes.bind_mouse_mode('middle',['control'], s)

def initialize_map_contour_mouse_modes(session):
    #z = ZoomMouseMode(session)
    s = SelectVolumeToContour(session)
    v = ContourSelectedVolume(session, s, True)
    session.ui.mouse_modes.bind_mouse_mode('wheel',['control'], s)
    # session.ui.mouse_modes.bind_mouse_mode('wheel',[], v)
    session.ui.mouse_modes.bind_mouse_mode('wheel',['alt'], v)

class RotateMouseMode(RotateMouseMode_Base):
    PAD_FRACTION = 0.3
    name = "clipper rotate"
    def mouse_double_click(self,event):
        x,y = event.position()
        pick = self.view.picked_object(x,y)
        atoms, sym = self._pick_atom_or_bond_atoms(pick)
        if atoms is not None:
            self.session.selection.clear()
            atoms.selected = True
            atoms.intra_bonds.selected=True
            self._center_on_picked_atoms(atoms, sym)
    
    def _pick_atom_or_bond_atoms(self, pick):
        if pick is None:
            return None, None
        atoms = None
        from chimerax.atomic import Atoms
        if hasattr(pick, 'atom'):
            atoms = Atoms([pick.atom])
        elif hasattr(pick, 'bond'):
            b = pick.bond
            atoms = Atoms(b.atoms)
        from .symmetry import PickedSymAtom, PickedSymBond
        if isinstance(pick, (PickedSymAtom, PickedSymBond)):
            from .util import rtop_frac_as_place
            sym = rtop_frac_as_place(pick.sym, atoms[0].structure.parent.cell)
        else:
            sym = None
        
        return atoms, sym
    
    def _center_on_picked_atoms(self, atoms, sym):
        v = self.session.view
        cofr = v.center_of_rotation
        cofr_method = v.center_of_rotation_method
        if sym is not None:
            v.move(sym)
            cofr = v.center_of_rotation = sym.inverse()*cofr

        center = atoms.coords.mean(axis=0)
        cd = v.camera.view_direction()
        shift = cofr-center
        from chimerax.geometry import translation
        v.move(translation(shift))
        v.center_of_rotation = center
        v.center_of_rotation_method = cofr_method


class ShiftToReferenceAsuMenuEntry(SelectContextMenuAction):
    dangerous=True
    def label(self, *args):
        return "Move symmetry fragment here"
    
    def criteria(self, session, event):
        from .symmetry import PickedSymAtom, PickedSymBond
        x,y = event.position()
        pick = session.view.picked_object(x,y)
        return isinstance(pick, (PickedSymAtom, PickedSymBond))
    
    def callback(self, session, event):
        x,y = event.position()
        pick = session.view.picked_object(x,y)
        if hasattr(pick, 'atom'):
            atom = pick.atom
        else:
            atom = pick.bond.atoms[0]
        from .util import rtop_frac_as_place
        m = atom.structure
        sym = rtop_frac_as_place(pick.sym, m.parent.cell)
        for g in m.bonded_groups():
            if atom in g:
                break
        g.transform(sym)
        m.parent.triggers.activate_trigger('sym shifted atoms', (g, sym))
        


    
        
        




class ClipPlaneAdjuster(MouseMode):
    name = "clipper clip adjust"
    def __init__(self, session, zoom_mode):
        super().__init__(session)
        self._zoomer = zoom_mode

    def cleanup(self):
        pass

    def wheel(self, event):
        incr = event.wheel_value()/120
        
        z = self._zoomer
        z.clip_multiplier -= incr
        z.clip_exponent *= (1-incr*20)

class Z_Shift_CofR(MouseMode):
    name = "clipper z shift"
    def __init__(self, session):
        self._step_multiplier = 1
        super().__init__(session)
        from chimerax.graphics import Drawing
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
        from chimerax.geometry import translation
        t = translation(shift_vec)
        self.view.center_of_rotation += shift_vec
        self.view.center_of_rotation_method = cofr_method
        camera.set_position(t*camera.position)

    def _drawing_geometry(self):
        from chimerax.surface.shapes import box_geometry
        return box_geometry((-1,-1,-0.01), (1,1,0.01))

    def _set_drawing_position(self):
        from chimerax.geometry import Place, scale
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
    name = "clipper zoom"
    MINIMUM_EXPONENT = 1e-4
    MINIMUM_MULTIPLIER = 1e-3
    def __init__(self, session):
        super().__init__(session)
        self._clip_multiplier = 0.1
        self._clip_exponent_multiplier = 1/50
        self.zoom(0)

    @property
    def clip_multiplier(self):
        '''
        Multiplier applied to the camera-cofr distance to decide the position
        of the rear clipping plane. Clamped to the range (0..1)
        '''
        return self._clip_multiplier

    @clip_multiplier.setter
    def clip_multiplier(self, val):
        val = max(min(val, 0.2), self.MINIMUM_MULTIPLIER)
        self._clip_multiplier = val
        # Zoom with a delta-z of zero to force redraw
        self.zoom(0)
    

    @property
    def clip_exponent(self):
        return self._clip_exponent_multiplier

    @clip_exponent.setter
    def clip_exponent(self, val):
        val = max(min(val, 0.1), self.MINIMUM_EXPONENT)
        self._clip_exponent_multiplier = val
        

    @property
    def far_clip(self):
        v = self.session.main_view
        cp = v.clip_planes
        fc = cp.find_plane('far')
        if fc is not None:
            return fc
        else:
            from chimerax.graphics.clipping import CameraClipPlane
            c = v.camera
            cofr = v.center_of_rotation
            cpos = c.position.origin()
            clip_point = cpos + (1+self.clip_multiplier)* (cofr-cpos)
            camera_point = c.position.inverse() * clip_point
            fc = CameraClipPlane('far', (0,0,1),
                                camera_point, v)
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
            from chimerax.graphics.clipping import CameraClipPlane
            c = v.camera
            cofr = v.center_of_rotation
            cpos = c.position.origin()
            clip_point = cpos + (1-self.clip_multiplier) * (cofr-cpos)
            camera_point = c.position.inverse() * clip_point
            nc = CameraClipPlane('near', (0,0,-1),
                                camera_point, v)
            v.clip_planes.add_plane(nc)
            return nc

    def cleanup(self):
        pass

    def adjusted_clip_multiplier(self, mult, dist):
        from math import exp
        return 1-(1-mult)*exp(-dist*self._clip_exponent_multiplier)
    
    def far_clip_point(self, origin, view_vec, dist):
        return origin + view_vec*(1+self.adjusted_clip_multiplier(self._clip_multiplier, dist))
    
    def near_clip_point(self, origin, view_vec, dist):
        return origin + view_vec*(1-self.adjusted_clip_multiplier(self._clip_multiplier, dist))



    def zoom(self, delta_z, stereo_scaling = False):
        import numpy
        v = self.view
        c = v.camera
        cofr = v.center_of_rotation
        if stereo_scaling and c.name == 'stereo':
            v.stereo_scaling(delta_z)
        if c.name == 'orthographic':
            c.field_width = max(c.field_width - delta_z, self.pixel_size())
            # TODO: Make camera field_width a property so it knows to redraw.
            from chimerax.geometry import place
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
            self.far_clip.plane_point = self.far_clip_point(new_origin, new_view_vec, new_view_distance)
            self.near_clip.plane_point = self.near_clip_point(new_origin, new_view_vec, new_view_distance)
            c.redraw_needed = True
        else:
            shift = c.position.transform_vector((0, 0, delta_z))
            v.translate(shift)
            new_origin = c.position.origin()
            vd = c.view_direction()
            distance = numpy.dot(cofr-new_origin, vd)
            self.far_clip.plane_point = self.far_clip_point(new_origin, vd*distance, distance)
            self.near_clip.plane_point = self.near_clip_point(new_origin, vd*distance, distance)

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
    name = "clipper contour select"
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
    name = "clipper contour adjust"
    def __init__(self, session, selector, symmetrical=True, sensitivity=0.02):
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
        self.sensitivity = sensitivity
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
            step = self.sensitivity * d * sd
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
