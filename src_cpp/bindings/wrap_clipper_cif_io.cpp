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
#include <pybind11/stl.h>

#include "type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-cif.h>


namespace py=pybind11;
using namespace clipper;

void declare_cif_io(py::module& m)
{
    py::class_<CIFfile> cif_file(m, "CIFfile");
    cif_file
        .def(py::init<>())
        .def("open_read", &CIFfile::open_read)
        .def("close_read", &CIFfile::close_read)
        .def_property_readonly("spacegroup", &CIFfile::spacegroup)
        .def_property_readonly("cell", &CIFfile::cell)
        .def_property_readonly("hkl_sampling", &CIFfile::hkl_sampling)
        .def_property_readonly("resolution", [](const CIFfile& self)
        {
            return self.resolution();
        })
        .def("import_hkl_info", &CIFfile::import_hkl_info)
        .def("import_hkl_data", &CIFfile::import_hkl_data)
        .def_property_readonly("contains_phases_predicate", &CIFfile::contains_phases_p)
        ;
} // declare_cif_file

void init_cif_io(py::module& m)
{
    declare_cif_io(m);
} // init_cif_file
