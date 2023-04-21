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
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "type_conversions.h"
#include <clipper/clipper.h>

#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;

void init_atomsf(py::module&m)
{
py::class_<AtomShapeFn> atomsf(m, "AtomShapeFn");
atomsf
    .def(py::init<>())
    .def(py::init<const Atom&>())
    .def(py::init<const Coord_orth&, const String&, const ftype, const ftype>())
    .def(py::init<const Coord_orth&, const String&, const U_aniso_orth&, const ftype>())
    .def("init", (void (AtomShapeFn::*)(const Atom&)) &AtomShapeFn::init)
    .def("init", (void (AtomShapeFn::*)(const Coord_orth&, const String&, const ftype, const ftype)) &AtomShapeFn::init)
    .def("init", (void (AtomShapeFn::*)(const Coord_orth&, const String&, const U_aniso_orth&, const ftype)) &AtomShapeFn::init)
    .def("f", (ftype (AtomShapeFn::*)(const Coord_reci_orth&) const) &AtomShapeFn::f)
    .def("f", (ftype (AtomShapeFn::*)(const ftype&) const) &AtomShapeFn::f)
    .def("rho", (ftype (AtomShapeFn::*)(const Coord_orth&) const) &AtomShapeFn::rho)
    .def("rho", (ftype (AtomShapeFn::*)(const ftype&) const) &AtomShapeFn::rho)
    .def("rho_grad", [] (const AtomShapeFn& self, const Coord_orth& xyz)
    {
        ftype rho;
        std::vector<ftype> grad;
        self.rho_grad(xyz, rho, grad);
        return py::make_tuple(rho, grad);
    })
    .def("rho_curv", [] (const AtomShapeFn& self, const Coord_orth& xyz)
    {
        ftype rho;
        std::vector<ftype> grad;
        Matrix<ftype> curv;
        self.rho_curv(xyz, rho, grad, curv);
        return py::make_tuple(rho, grad, curv);
    })
    .def("agarwal_params", &AtomShapeFn::agarwal_params)
    .def("set_agarwal_params", [](AtomShapeFn& self, const std::vector<AtomShapeFn::TYPE>& params) {self.agarwal_params()=params;})
    ;


py::enum_<AtomShapeFn::TYPE>(atomsf, "TYPE")
    .value("X", AtomShapeFn::TYPE::X)
    .value("Y", AtomShapeFn::TYPE::Y)
    .value("Z", AtomShapeFn::TYPE::Z)
    .value("Uiso", AtomShapeFn::TYPE::Uiso)
    .value("Occ", AtomShapeFn::TYPE::Occ)
    .value("U11", AtomShapeFn::TYPE::U11)
    .value("U22", AtomShapeFn::TYPE::U22)
    .value("U33", AtomShapeFn::TYPE::U33)
    .value("U12", AtomShapeFn::TYPE::U12)
    .value("U13", AtomShapeFn::TYPE::U13)
    .value("U23", AtomShapeFn::TYPE::U23)
    .export_values();


} // init_atomsf
