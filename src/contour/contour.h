// Marching cubes contouring code with permission from Tom Goddard, UCSF

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
// Note that this software makes use of modified versions of the Clipper, LibCCP4
// and MMDB libraries, as well as portions of the Intel Math Kernel Library. Each
// of these is redistributed under its own license terms.

// ----------------------------------------------------------------------------
// Compute an constant intensity surface from volume data.
//
#ifndef CONTOUR_HEADER_INCLUDED
#define CONTOUR_HEADER_INCLUDED

//
// The data values can be any of the standard C numeric types.
//
// The grid size array is in x, y, z order.
//
// The grid value for index (i0,i1,i2) where 0 <= ik < size[k] is
//
//	grid[i0*stride[0] + i1*stride[1] + i2*stride[2]]
//
// Returned vertex and triangle arrays should be freed with free_surface().
//

namespace Contour_Calculation
{
typedef unsigned int Index; // grid and edge indices, and surface vertex indices
typedef long Stride;	    // Array strides and pointer offsets, signed

class Contour_Surface
{
 public:
  virtual ~Contour_Surface() {};
  virtual Index vertex_count() = 0;
  virtual Index triangle_count() = 0;
  virtual void geometry(float *vertex_xyz, Index *triangle_vertex_indices) = 0;
  virtual void normals(float *normals) = 0;
};

template <class Data_Type>
Contour_Surface *surface(const Data_Type *grid,
			 const Index size[3], const Stride stride[3],
			 float threshold, bool cap_faces);
}

#include "contour.cpp"	// template implementation


#endif
