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

#include "type_conversions.h"
#include <clipper/clipper.h>

namespace py=pybind11;
using namespace clipper;

void declare_map_stats(py::module& m)
{
    py::class_<Map_stats>(m, "Map_stats")
        .def(py::init<Xmap<ftype32>>())
        .def(py::init<Xmap<ftype64>>())
        .def(py::init<NXmap<ftype32>>())
        .def(py::init<NXmap<ftype64>>())
        .def_property_readonly("mean", &Map_stats::mean)
        .def_property_readonly("std_dev", &Map_stats::std_dev)
        .def_property_readonly("min", &Map_stats::min)
        .def_property_readonly("max", &Map_stats::max)
        .def_property_readonly("range", &Map_stats::range)
        ;
}



void init_map_utils(py::module& m)
{
    declare_map_stats(m);
} // init_map_utils
