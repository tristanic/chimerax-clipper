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
#include <pybind11/stl.h>

#include "type_conversions.h"
#include "symops.h"

#include "numpy_helper.h"


namespace py=pybind11;
using namespace clipper;


int nrows_from_shape_string(const std::string& shape)
{
    if (shape == "4x4")
        return 4;
    else if (shape == "3x4")
        return 3;
    else
        throw std::logic_error("Shape should be one of \"4x4\" or \"3x4\"!");

}
void export_rot_trn_(const Mat33<ftype>& rot, const Vec3<ftype>& trn, int nrows, ftype* ptr)
{
    const static ftype lastrow[4] = {0,0,0,1};
    for (int j=0; j<nrows; ++j)
    {
        for (int k=0; k<4; ++k)
        {
            if (j==3)
                *ptr++ = lastrow[k];
            else
            {
                if (k==3)
                    *ptr++ = trn[j];
                else
                    *ptr++ = rot(j,k);
            }

        }
    }
}

void all_matrices_frac_(const RTop_fracs& self, int nrows, py::array_t<ftype>& target)
{
    int n = self.size();
    check_numpy_array_shape(target, {n, nrows, 4}, false);
    ftype* tptr = (ftype*)target.request().ptr;
    for (int i=0; i<n; ++i)
    {
        const RTop_frac& thisop = self[i];
        const Mat33<ftype>& rot = thisop.rot();
        const Vec3<ftype>& trn = thisop.trn();
        export_rot_trn_(rot, trn, nrows, tptr);
        tptr+= nrows*4;
    }
}

void all_matrices_orth_(const RTop_fracs& self, const Cell& cell, int nrows, py::array_t<ftype>& target)
{
    int n = self.size();
    check_numpy_array_shape(target, {n, nrows, 4}, false);
    ftype* tptr = (ftype*)target.request().ptr;
    for (int i=0; i<n; ++i)
    {
        const RTop_orth& thisop = self[i].rtop_orth(cell);
        const Mat33<ftype>& rot = thisop.rot();
        const Vec3<ftype>& trn = thisop.trn();
        export_rot_trn_(rot, trn, nrows, tptr);
        tptr+=nrows*4;
    }
}


namespace py=pybind11;
using namespace clipper;

void declare_symops(py::module& m)
{
    py::class_<RTop_fracs>(m, "RTop_fracs")
        .def(py::init<>())
        .def(py::init<std::vector<RTop_frac>&>())
        .def("__getitem__", &RTop_fracs::at)
        .def("__setitem__", &RTop_fracs::replace)
        .def("append", (void (RTop_fracs::*)(const Symop& op, const Coord_frac& offset)) &RTop_fracs::append)
        .def("append", (void (RTop_fracs::*)(const RTop_frac& op)) &RTop_fracs::append)
        .def("pop", &RTop_fracs::pop)
        .def("__len__", &RTop_fracs::size)
        .def("all_matrices_frac",
            [](const RTop_fracs& self, const std::string& shape) -> py::array_t<ftype>
            {
                auto nrows = nrows_from_shape_string(shape);
                py::array_t<ftype> target({(int)self.size(), nrows, 4});
                all_matrices_frac_(self, nrows, target);
                return target;
            })
        .def("all_matrices_frac",
            [](const RTop_fracs& self, const std::string& shape, py::array_t<ftype> target) -> void
            {
                auto nrows = nrows_from_shape_string(shape);
                all_matrices_frac_(self, nrows, target);
            })
            .def("all_matrices_orth",
                [](const RTop_fracs& self, const Cell& cell, const std::string& shape) -> py::array_t<ftype>
                {
                    auto nrows = nrows_from_shape_string(shape);
                    py::array_t<ftype> target({(int)self.size(), nrows, 4});
                    all_matrices_orth_(self, cell, nrows, target);
                    return target;
                })
            .def("all_matrices_orth",
                [](const RTop_fracs& self, const Cell& cell, const std::string& shape, py::array_t<ftype> target) -> void
                {
                    auto nrows = nrows_from_shape_string(shape);
                    all_matrices_orth_(self, cell, nrows, target);
                })
        ;
}

void init_symops(py::module& m)
{
    declare_symops(m);
}
