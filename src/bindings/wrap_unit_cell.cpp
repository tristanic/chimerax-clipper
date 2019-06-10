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
#include "numpy_helper.h"

#include "type_conversions.h"
#include <clipper/clipper.h>
#include "unit_cell.h"


namespace py=pybind11;
using namespace clipper;

const char* unit_cell_docstring =
"Describes the unit cell surrounding a reference coordinate (typically the "
"centroid of an atomic model). Provides a set of symmetry operators defined "
"with the modelled asu as the identity, along with the translations necessary "
"to bring symmetry copies back into the same unit cell. Also provides methods "
"to quickly find all the symmetry operators necessary to pack an arbitrary box "
"in space with copies of the atomic model.";

void declare_unit_cell(py::module& m)
{
    py::class_<Unit_Cell>(m, "Unit_Cell", unit_cell_docstring)
        //.def(py::init<>())
        .def(py::init<const Atom_list&, const Cell&,
                      const Spacegroup&, const Grid_sampling&, int>())
        .def_property_readonly("symops", &Unit_Cell::symops)
        .def_property_readonly("inv_symops", &Unit_Cell::inv_symops)
        .def_property("ref_coord",
            [](const Unit_Cell& self) {return self.ref_coord(); },
            [](Unit_Cell& self, const Coord_frac& new_ref)
                { self.set_ref_coord(new_ref); }
        )
        .def_property_readonly("cell", &Unit_Cell::cell)
        .def_property_readonly("spacegroup", &Unit_Cell::spacegroup)
        .def_property_readonly("grid", &Unit_Cell::grid)
        .def_property_readonly("ref_box", &Unit_Cell::ref_box)
        .def_property_readonly("origin", &Unit_Cell::min)
        .def_property_readonly("top_corner", &Unit_Cell::max)
        .def("update_reference_model_bounds", &Unit_Cell::update_reference_model_bounds)
        .def("all_symops_in_box",
            (RTop_fracs (Unit_Cell::*)(const Coord_orth& origin_xyz,
               const Vec3<int>& box_size_uvw,
               bool always_include_identity, int sample_frequency) const)
               &Unit_Cell::all_symops_in_box,
               py::arg("origin_xyz"), py::arg("box_size_uvw"),
               py::arg("always_include_identity")=true, py::arg("sample_frequency")=2
         )
         .def("all_symops_in_box",
         [](const Unit_Cell& self, py::array_t<ftype> origin_xyz, py::array_t<int> box_size_uvw,
            bool always_include_identity, int sample_frequency)
         {
             check_numpy_array_shape(origin_xyz, {3}, true);
             check_numpy_array_shape(box_size_uvw, {3}, true);
             ftype* optr = (ftype*)origin_xyz.request().ptr;
             int* bptr = (int*)box_size_uvw.request().ptr;
             auto orig = Coord_orth(optr[0], optr[1], optr[2]);
             auto box = Vec3<int>(bptr[0], bptr[1], bptr[2]);
             return self.all_symops_in_box(orig, box, always_include_identity, sample_frequency);
         },
         py::arg("origin_xyz"), py::arg("box_size_uvw"),
         py::arg("always_include_identity")=true, py::arg("sample_frequency")=2
         )
         ;
} // declare_unit_cell


void init_unit_cell(py::module& m)
{
    declare_unit_cell(m);
}
