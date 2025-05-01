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

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>


#include <clipper/clipper.h>

#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;

template<class C, typename T>
py::class_<C> add_rtop_base_wrappings(py::class_<C>& derived)
{
    using Base=RTop<T>;
    derived
        //.def("inverse", [](const Base& self){ return self.inverse(); } )
        .def("equals", [](const Base& self, const Base& other, const T& tol) { return self.equals(other, tol); })
        .def_property("rot",
             [](const Base& self) { return self.rot(); },
             [](Base& self, const Mat33<T>& rot) { self.rot() = rot; }
        )
        .def_property("trn",
             [](const Base& self) { return self.trn(); },
             [](Base& self, const Vec3<T>& trn) { self.trn() = trn; }
        )
        //.def_static("identity", &Base::identity)
        //.def_static("null", &Base::null)
        .def("is_null",  [](const Base& self) { return self.is_null(); })

        /* Don't wrap base class format method, since this is overridden in some
           derived classes. It will have to be individually wrapped in each
           derived class.
         */
         /*
        .def("format", &Class::format)
        .def("__repr__", [&pyclass_name](const Class& self){
             return std::string("Clipper ")+pyclass_name+std::string(": ")+std::string(self.format());
         })
        .def("__str__", [](const Class& self) { return std::string(self.format()); })
        */
        .def("__mul__", [](const Base& self, const Vec3<T>& v){ return Cself*v; }, py::is_operator())
        .def("__mul__", [](const Base& self, const Base& other){ return C(self*other); }, py::is_operator())
        ;
        return derived;
} // init_rtop
