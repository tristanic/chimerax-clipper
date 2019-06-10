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

def delayed_reaction(triggerset, trigger_name, initiator_func, initiator_args, ready_test_func,
    final_func, final_func_args):
    '''
    Designed to work together with threaded (Python or C++) objects, to
    start a long-running threaded task and then automatically apply the
    result (in a later GUI update) when done. Can also be used to simply
    conveniently call the desired callback once on the next firing of the
    trigger, by setting ready_test_func to None.
    Args:
        triggerset: the triggerset providing the trigger (e.g. session.triggers)
        trigger_name: the name of the trigger (e.g. 'new frame')
        initiator_func: A handle to the function that kicks off the threaded
            process. Should not return anything.
        initiator_args: A tuple of arguments to be applied to initiator_func
        ready_test_func: Should return True when the threaded task is done,
            false otherwise. Set it to None to just run on the next firing
            of the trigger.
        final_func: Task to run once the thread is done.
        final_func_args: A tuple of arguments to be applied to final_func
            (e.g. to tell it what to do with the result)
    '''
    initiator_func(*initiator_args)
    class _cb:
        def __init__(self, triggerset, trigger_name, ready_test_func, final_func, final_func_args):
            self.tf = ready_test_func
            self.ff = final_func
            self.ff_args = final_func_args
            self.handler = triggerset.add_handler(trigger_name, self.callback)
        def callback(self, *_):
            if self.tf is None or self.tf():
                self.ff(*self.ff_args)
                from chimerax.core.triggerset import DEREGISTER
                return DEREGISTER
    cb = _cb(triggerset, trigger_name, ready_test_func, final_func, final_func_args)

class Delayed_Reaction_Tester:
    def __init__(self, frame_delay):
        self._counter = 0
        self._frame_delay = frame_delay

    def initiator_func(self):
        pass

    def ready_test_func(self):
        self._counter += 1
        if (self._counter >= self._frame_delay):
            return True
        return False

    def final_func(self, to_print):
        print(to_print)
