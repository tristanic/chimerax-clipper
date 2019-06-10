// PyBind11 Python bindings for Clipper
// Copyright (C) 2016-2019 Tristan Croll, University of Cambridge
//
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

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "type_conversions.h"
#include <clipper/clipper.h>
#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;



void declare_nx_operator(py::module& m)
{
    py::class_<NX_operator>(m, "NX_operator")
        .def(py::init<const Xmap_base&, const NXmap_base&, const RTop_orth&>())
        .def(py::init<const Cell&, const Grid_sampling&, const NXmap_base&, const RTop_orth&>())
        .def("coord_map", &NX_operator::coord_map)
        .def("coord_frac", &NX_operator::coord_frac)
        .def("nxmap_data_on_xmap_grid_linear_interp",
            [](const NX_operator& self, const NXmap<ftype32>& nxmap, const Coord_grid& c)
            {
                return self.template nxmap_data<Interp_linear, ftype64, NXmap<ftype32>>(nxmap, c);
            })
        .def("nxmap_data_on_xmap_grid_linear_interp",
            [](const NX_operator& self, const NXmap<ftype64>& nxmap, const Coord_grid& c)
            {
                return self.template nxmap_data<Interp_linear, ftype64, NXmap<ftype64>>(nxmap, c);
            })
        .def("nxmap_data_on_xmap_grid_cubic_interp",
            [](const NX_operator& self, const NXmap<ftype32>& nxmap, const Coord_grid& c)
            {
                return self.template nxmap_data<Interp_cubic, ftype64, NXmap<ftype32>>(nxmap, c);
            })
        .def("nxmap_data_on_xmap_grid_cubic_interp",
            [](const NX_operator& self, const NXmap<ftype64>& nxmap, const Coord_grid& c)
            {
                return self.template nxmap_data<Interp_cubic, ftype64, NXmap<ftype64>>(nxmap, c);
            })
        .def("xmap_data_on_nxmap_grid_linear_interp",
            [](const NX_operator& self, const Xmap<ftype32>& xmap, const Coord_grid& c)
            {
                return self.template xmap_data<Interp_linear, ftype64, Xmap<ftype32>>(xmap, c);
            })
        .def("xmap_data_on_nxmap_grid_linear_interp",
            [](const NX_operator& self, const Xmap<ftype64>& xmap, const Coord_grid& c)
            {
                return self.template xmap_data<Interp_linear, ftype64, Xmap<ftype64>>(xmap, c);
            })
        .def("xmap_data_on_nxmap_grid_cubic_interp",
            [](const NX_operator& self, const Xmap<ftype32>& xmap, const Coord_grid& c)
            {
                return self.template xmap_data<Interp_cubic, ftype64, Xmap<ftype32>>(xmap, c);
            })
        .def("xmap_data_on_nxmap_grid_linear_interp",
            [](const NX_operator& self, const Xmap<ftype64>& xmap, const Coord_grid& c)
            {
                return self.template xmap_data<Interp_cubic, ftype64, Xmap<ftype64>>(xmap, c);
            })
        ;
}


void init_nx_operator(py::module& m)
{
    declare_nx_operator(m);
} // init_nx_operator
