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
#include <pybind11/numpy.h>

#include <clipper/clipper.h>

#include <vector>

namespace py=pybind11;
using ssize_t=py::ssize_t;
using namespace clipper;

// Numpy array that is guaranteed C-contiguous and of the right dtype by the time it
// reaches C++: pybind11 substitutes a converted copy if necessary.  Correct only for
// arrays we READ - never for an output array, where a silent copy would swallow
// everything we write into it (use check_c_contiguous_3d() for those).
template<typename T>
using c_array_in = py::array_t<T, py::array::c_style | py::array::forcecast>;

// Numpy array we WRITE into.  ExtraFlags=0 disables the forcecast that py::array_t
// applies by default: a dtype mismatch must raise, not silently hand us a converted
// copy whose contents are discarded when the call returns.
template<typename T>
using array_out = py::array_t<T, 0>;

// Check that an OUTPUT array is 3D and C-contiguous, and throw if not.  Arrays we
// fill by walking a bare pointer must be checked rather than coerced: declaring them
// py::array::c_style would let pybind11 hand us a contiguous copy instead, and the
// caller's array would silently never be written.  Takes py::array so it binds any
// py::array_t instantiation without triggering a conversion.
inline void check_c_contiguous_3d(const py::array& arr, const char* name)
{
    if (arr.ndim() != 3)
        throw std::runtime_error(std::string(name) + " must be a 3D array, but has "
            + std::to_string(arr.ndim()) + " dimension(s)!");
    if (!(arr.flags() & py::array::c_style))
        throw std::runtime_error(std::string(name) + " must be a C-contiguous array! "
            "It is written to in place, so a strided or transposed view cannot be used.");
}

// check that the given Numpy array matches expected dimensions, and throw an
// error if not. direction is true for incoming, false for outgoing.
template<typename T>
void check_numpy_array_shape(py::array_t<T> target, std::vector<int> dim, bool direction)
{
    auto buf = target.request();
    bool fail=false;
    if ((size_t)buf.ndim != dim.size())
        fail = true;
    else
        for (ssize_t i=0; i<buf.ndim; ++i)
            if (buf.shape[i] != dim[i])
                fail = true;
    if (fail) {
        auto shapet_txt = std::string("( ");
        auto shapeb_txt = std::string("( ");
        for (ssize_t i=0; i<buf.ndim; ++i)
            shapet_txt += std::to_string(buf.shape[i]) + " ";
        for (size_t i=0; i<dim.size(); ++i)
            shapeb_txt += std::to_string(dim[i]) + " ";
        shapet_txt += ")";
        shapeb_txt += ")";

        std::string message;

        auto msg = std::string("Array shape mismatch! ");
        if (direction)
            msg += "Input array shape is ";
        else
            msg += "Target array shape is ";
        msg += shapet_txt;
        msg += ", while expected shape is ";
        msg += shapeb_txt;
        msg += ".";
        throw std::runtime_error(msg.c_str());
    }
}

template<class C, typename T>
py::array_t<T> array_as_numpy_1d(const C& v, int n)
{
    py::array_t<T> ret(n);
    auto buf = ret.request();
    T* ptr = (T*)buf.ptr;
    for (int i=0; i<n; ++i)
        ptr[i] = v[i];
    return ret;
}

template<class C, typename T>
void array_as_numpy_1d(const C& v, int n, py::array_t<T> target)
{
    check_numpy_array_shape(target, {n}, true);
    auto buf = target.request();
    T* ptr = (T*)buf.ptr;
    for (int i=0; i<n; ++i)
        ptr[i] = v[i];
}


template <class C, typename T>
void fill_array_from_numpy_1d(C& v, int n, py::array_t<T> arr)
{
    check_numpy_array_shape(arr, {n}, true);
    auto buf = arr.request();
    T* ptr = (T*)buf.ptr;
    for (int i=0; i<n; ++i)
        v[i] = ptr[i];
}

template <class HKLdtype, class dtype>
py::array_t<dtype> hkl_data_export_numpy(const HKLdtype& self, const int& size)
{
    auto ret = py::array_t<dtype>(size);
    dtype* ptr = (dtype*)ret.request().ptr;
    self.data_export(ptr);
    return ret;
} // hkl_data_export_numpy

template <class HKLdtype, class dtype>
void hkl_data_import_numpy(HKLdtype& self, const int& size, py::array_t<dtype> vals)
{
    check_numpy_array_shape(vals, {size}, true);
    dtype* ptr = (dtype*)vals.request().ptr;
    self.data_import(ptr);
} // hkl_data_import_numpy

template <typename T>
std::unique_ptr<Vec3<T>> new_vec3_from_numpy(py::array_t<T> vals)
{
    check_numpy_array_shape(vals, {3}, true);
    return std::unique_ptr<Vec3<T>>(new Vec3<T>(vals.at(0), vals.at(1), vals.at(2)));
}

template <typename T>
std::unique_ptr<Mat33<T>> new_mat33_from_numpy(py::array_t<T> vals)
{
    check_numpy_array_shape(vals, {3,3}, true);
    auto r = vals.template unchecked<2>();
    return std::unique_ptr<Mat33<T>>(new Mat33<T>(
        r(0,0),r(0,1),r(0,2),r(1,0),r(1,1),r(1,2),r(2,0),r(2,1),r(2,2)
    ));
}
