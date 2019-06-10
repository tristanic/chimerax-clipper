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
from .clipper_python import Util as _util

### Constants ###
from math import pi
eightpi2 = 8*pi*pi;

### Math functions ###

def atanh(x):
    ''' Hyperbolic tan '''
    from math import atanh
    return atanh(x)

def bfactor_to_uvalue(b_factor):
    ''' Convert isotropic B-factor to U-value '''
    return _util.b2u(b_factor)

def uvalue_to_bfactor(u_value):
    ''' Convert isotropic U-value to B-factor '''
    return _util.u2b(u_value)

def bessel_i0(x):
    ''' Modified Bessel function of the first kind '''
    return _util.bessel_i0(x)

def d2rad(angle):
    ''' Convert degrees to radians '''
    from math import radians
    return radians(angle)

def rad2d(angle):
    ''' Convert radians to degrees '''
    from math import degrees
    return degrees(angle)

def int_ceiling(fp_num):
    ''' Round the given number up and return as an integer. '''
    from math import ceil
    return ceil(fp_num)

def int_floor(fp_num):
    ''' Round the given number down and return as an integer. '''
    from math import floor
    return floor(fp_num)

def int_round(fp_num):
    ''' Round the given number to the nearest integer. '''
    return round(fp_num)

def invsim(x):
    ''' Inverse Sim function: I1(X)/I0(X). '''
    return _util.invsim(x)



### Useful utility functions ###
def get_minmax_grid(coords, cell, grid_sampling):
    '''
    Return a numpy array containing the (min, max) grid coordinates
    of a box encompassing the given coordinates.
    Args:
        coords ( [n*3 float]):
            An n*3 array of coordinates
        cell (clipper.Cell object)
        grid_sampling (clipper.Grid_sampling object)
    '''
    return _util.get_minmax_grid(coords, cell, grid_sampling)
