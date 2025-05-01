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
#include <vector>

#include "type_conversions.h"
#include <clipper/clipper.h>

namespace py=pybind11;
using namespace clipper;

void init_hkl_info(py::module &m) {

    py::class_<HKL_info>(m, "HKL_info")
        .def(py::init<>())
        .def(py::init<const Spacegroup&, const Cell&, const Resolution&, const bool&>())
        .def("init", (void (HKL_info::*)(const Spacegroup&, const Cell&, const Resolution&, const bool&)) &HKL_info::init)
        .def("init", (void (HKL_info::*)(const Spacegroup&, const Cell&, const HKL_sampling&, const bool&)) &HKL_info::init)
        .def_property_readonly("is_null", &HKL_info::is_null)
        .def_property_readonly("cell", &HKL_info::cell)
        .def_property_readonly("spacegroup", &HKL_info::spacegroup)
        .def_property_readonly("hkl_sampling", &HKL_info::hkl_sampling)
        .def_property_readonly("resolution", &HKL_info::resolution)
        .def("generate_hkl_list", &HKL_info::generate_hkl_list)
        .def("add_hkl_list", &HKL_info::add_hkl_list)
        .def("add_hkl_list", [](HKL_info& self, py::array_t<int> hkls) {
            auto buf = hkls.request();
            if (buf.ndim != 2 || buf.shape[1] !=3)
                throw std::runtime_error("HKLs must be an nx3 Numpy array!");
            std::vector<HKL> vhkls;
            for (size_t i=0; i<buf.shape[0]; ++i)
            {
                vhkls.push_back(HKL(hkls.at(i,0), hkls.at(i,1), hkls.at(i,2)));
            }
            self.add_hkl_list(vhkls);
        })
        .def_property_readonly("num_reflections", &HKL_info::num_reflections)
        .def("hkl_of", &HKL_info::hkl_of)
        .def("index_of", &HKL_info::index_of)
        .def("invresolsq", &HKL_info::invresolsq)
        .def_property_readonly("invresolsq_range", &HKL_info::invresolsq_range)
        .def("hkl_class", &HKL_info::hkl_class)
        .def("find_sym", &HKL_info::find_sym)
        .def_property_readonly("first", &HKL_info::first);

    py::class_<HKL_info::HKL_reference_base>(m, "_HKL_reference_base")
        .def_property_readonly("base_hkl_info", &HKL_info::HKL_reference_base::base_hkl_info)
        .def_property_readonly("index", &HKL_info::HKL_reference_base::index)
        .def("invresolsq", (ftype (HKL_info::HKL_reference_base::*)(const HKL_data_base&) const)&HKL_info::HKL_reference_base::invresolsq)
        .def("invresolsq", (ftype (HKL_info::HKL_reference_base::*)() const)&HKL_info::HKL_reference_base::invresolsq)
        .def("last", &HKL_info::HKL_reference_base::last);

    py::class_<HKL_info::HKL_reference_index, HKL_info::HKL_reference_base>(m, "HKL_reference_index")
        .def(py::init<>())
        .def(py::init<const HKL_info&, const int&>())
        .def_property_readonly("hkl", &HKL_info::HKL_reference_index::hkl)
        .def_property_readonly("hkl_class", &HKL_info::HKL_reference_index::hkl_class)
        // avoid creating a new Python object every time we increment
        .def("next", [](HKL_info::HKL_reference_index& self) { self.next(); });

    py::class_<HKL_info::HKL_reference_coord, HKL_info::HKL_reference_base>(m, "HKL_reference_coord")
        .def(py::init<>())
        .def(py::init<const HKL_info&, const HKL&>())
        .def_property("hkl",
            &HKL_info::HKL_reference_coord::hkl,
            [](HKL_info::HKL_reference_coord& self, const HKL& hkl__) { self.set_hkl(hkl__); }
        )
        .def_property_readonly("sym", &HKL_info::HKL_reference_coord::sym)
        .def_property_readonly("friedel", &HKL_info::HKL_reference_coord::friedel)
        // avoid creating new Python objects when incrementing
        .def("next", [](HKL_info::HKL_reference_coord& self) { self.next(); })
        .def("next_h", [](HKL_info::HKL_reference_coord& self) { self.next_h(); })
        .def("next_k", [](HKL_info::HKL_reference_coord& self) { self.next_k(); })
        .def("next_l", [](HKL_info::HKL_reference_coord& self) { self.next_l(); })
        .def("prev_h", [](HKL_info::HKL_reference_coord& self) { self.prev_h(); })
        .def("prev_k", [](HKL_info::HKL_reference_coord& self) { self.prev_k(); })
        .def("prev_l", [](HKL_info::HKL_reference_coord& self) { self.prev_l(); })
        ;



} //init_hkl_info
