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
#include <pybind11/numpy.h>
#include "numpy_helper.h"

#include "type_conversions.h"
#include <clipper/clipper.h>

namespace py=pybind11;
using namespace clipper;

void declare_clipper_util(py::module &m)
{
    py::class_<Util>(m, "Util")
        .def_static("nan_double", &Util::nand)
        .def_static("nan_float", &Util::nanf)
        // .def_static("isnan_double", (bool (Util::*)(const ftype64)) &Util::isnan)
        // .def_static("isnan_float", (bool (Util::*)(const ftype32)) &Util::isnan)
        .def_static("sim", &Util::sim)
        .def_static("invsim", &Util::invsim)
        .def_static("sim_integ", &Util::sim_integ)
        .def_static("sim_deriv", &Util::sim_deriv)
        .def_static("sim_deriv_recur", &Util::sim_deriv_recur)
        .def_static("atanh", &Util::atanh)
        .def_static("bessel_i0", &Util::bessel_i0)
        .def_static("u2b", &Util::u2b)
        .def_static("uvalue_to_bfactor", &Util::u2b)
        .def_static("b2u", &Util::b2u)
        .def_static("bfactor_to_uvalue", &Util::b2u)
        .def_static("anom_mean_float", &Util::mean<ftype32>)
        .def_static("anom_mean_double", &Util::mean<ftype64>)
        .def_static("anom_mean_sig_float", &Util::sig_mean<ftype32>)
        .def_static("anom_mean_sig_double", &Util::sig_mean<ftype64>)
        // Extra goodies for Python
        .def_static("get_minmax_grid",
        [](py::array_t<ftype> coords, const Cell& cell, const Grid_sampling& grid)
            -> py::array_t<int>
        {
            auto cbuf = coords.request();
            int n = cbuf.shape[0];
            if (cbuf.shape[1] != 3)
                throw std::logic_error("Input should be an array of 3D coordinates!");
            ftype* cptr = (ftype*)cbuf.ptr;
            py::array_t<int> out({2,3});
            int* optr = (int*)out.request().ptr;
            Coord_grid ref_min = Coord_orth(cptr[0], cptr[1], cptr[2]).coord_frac(cell).coord_grid(grid);
            Coord_grid ref_max = ref_min;
            for (int i=3; i<n*3; i+=3)
            {
                Coord_grid thiscoord =
                    Coord_orth(cptr[i], cptr[i+1], cptr[i+2]).coord_frac(cell).coord_grid(grid);
                for (size_t j=0; j<3; ++j)
                {
                    if (thiscoord[j] < ref_min[j]) ref_min[j] = thiscoord[j];
                    else if (thiscoord[j] > ref_max[j]) ref_max[j] = thiscoord[j];
                }
            }
            for (size_t i=0; i<3; i++)
            {
                optr[i] = ref_min[i];
                optr[3+i] = ref_max[i];
            }
            return out;
        })
        ;

} // declare_clipper_util

void init_clipper_util(py::module &m)
{
    declare_clipper_util(m);
}
