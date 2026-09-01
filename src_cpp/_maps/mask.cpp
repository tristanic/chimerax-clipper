// ChimeraX-Clipper
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
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>
#include <algorithm>

namespace py=pybind11;
using ssize_t = py::ssize_t;

// Read-only arrays: let pybind11 hand us a C-contiguous array of the right dtype,
// converting if it has to.  Needed as well as convenient - callers legitimately pass
// a shape tuple for dim and a float32 Atoms.coords for coords.
template<typename T>
using c_array_in = py::array_t<T, py::array::c_style | py::array::forcecast>;

// Every array below is read through a bare pointer, so its element count is a bounds
// contract.  min_size is what the C++ actually dereferences.
static void check_size(const py::array& arr, ssize_t min_size, const char* name)
{
    if (arr.size() < min_size)
        throw std::runtime_error(std::string(name) + " must have at least "
            + std::to_string(min_size) + " elements, but has "
            + std::to_string(arr.size()) + "!");
}

// The mask is written through a flat C-order index (see stamp_atom_range), so it must
// be 3D, C-contiguous, and shaped exactly as dim says.  Checked rather than coerced:
// declaring it c_style would let pybind11 substitute a contiguous copy, and every
// voxel we stamp would be silently discarded when the call returns.
static void check_mask_array(const py::array& map, const size_t* dim)
{
    if (map.ndim() != 3)
        throw std::runtime_error("Mask array must be 3D!");
    if (!(map.flags() & py::array::c_style))
        throw std::runtime_error("Mask array must be C-contiguous! It is written to "
            "in place, so a strided or transposed view cannot be used.");
    for (ssize_t i=0; i<3; ++i)
        if ((size_t)map.shape(i) != dim[i])
            throw std::runtime_error("Mask array shape does not match the given dimensions!");
}

template<typename T>
void affine_transform(T* coord, T* tf, T* out )
{
    for (size_t i=0; i<3; ++i)
    {
        auto j=i*4;
        out[i] = coord[0]*tf[j] + coord[1]*tf[j+1] + coord[2]*tf[j+2] + tf[j+3];
    }
}

template<typename T>
void index_range(T* coord, T* radius, T* step, size_t* dim, size_t* minc, size_t* maxc)
{
    for (size_t i=0; i<3; ++i)
    {
        size_t ci_pl = (size_t)ceil(coord[i]+radius[i]);
        ssize_t ci_mi = (ssize_t)floor(coord[i]-radius[i]);
        ci_mi = ci_mi < 0 ? 0 : ci_mi;
        ci_pl = ci_pl > dim[i]-1 ? dim[i]-1 : ci_pl;
        minc[i] = ci_mi;
        maxc[i] = ci_pl;
    }
}

template<typename T>
T squared_distance(T* ref_coord, T* box_coord, T* ijk_to_xyz)
{
    // T tf_box[3];
    // affine_transform(box_coord, ijk_to_xyz, tf_box);
    T sqd = 0;
    for (size_t i=0; i<3; ++i)
        sqd += pow(ref_coord[i]-box_coord[i], 2);
    return sqd;
}


// Stamp the spheres for atoms [i0, i1) into the shared map. All writes set the
// same value (1), so this is safe to run concurrently over disjoint atom ranges
// even where their voxel boxes overlap (a single-byte store of 1 is idempotent;
// no locks needed). r_grid and sq_rad are precomputed and read-only here.
static void stamp_atom_range(
    uint8_t* map, double* step, size_t* dim, double* ijk_to_xyz, double* xyz_to_ijk,
    double* coords, size_t i0, size_t i1, double* r_grid, double sq_rad)
{
    size_t box_min[3], box_max[3];
    double transformed[3], bc[3];
    for (size_t i=i0; i<i1; ++i)
    {
        affine_transform(coords+3*i, xyz_to_ijk, transformed);
        index_range(transformed, r_grid, step, dim, box_min, box_max);
        for (size_t u=box_min[0]; u<=box_max[0]; ++u) {
            for (size_t v=box_min[1]; v<= box_max[1]; ++v) {
                for (size_t w=box_min[2]; w<= box_max[2]; ++w) {
                    bc[0]=u; bc[1]=v; bc[2]=w;
                    if (squared_distance(transformed, bc, ijk_to_xyz) < sq_rad)
                    {
                        map[w + dim[2] * ( v + dim[1] * u)] = 1;
                    }
                }
            }
        }
    }
}

void generate_mask(
    uint8_t* map, double* origin, double* step, size_t* dim,
    double* ijk_to_xyz, double* xyz_to_ijk, double* coords, size_t n, double radius,
    size_t num_threads)
{
    double r_xyz[3], r_grid[3];
    for (size_t i=0; i<3; ++i)
        r_xyz[i] = radius + origin[i];
    affine_transform(r_xyz, xyz_to_ijk, r_grid);
    double sq_rad = pow(radius/step[0],2); //radius*radius;

    // Small jobs aren't worth the thread overhead.
    if (num_threads <= 1 || n < 1024)
    {
        stamp_atom_range(map, step, dim, ijk_to_xyz, xyz_to_ijk, coords, 0, n, r_grid, sq_rad);
        return;
    }
    size_t nthreads = std::min(num_threads, (size_t)n);
    std::vector<std::thread> threads;
    threads.reserve(nthreads);
    size_t per = (n + nthreads - 1) / nthreads;
    for (size_t t=0; t<nthreads; ++t)
    {
        size_t i0 = t*per;
        size_t i1 = std::min(i0+per, n);
        if (i0 >= i1) break;
        threads.emplace_back(stamp_atom_range, map, step, dim, ijk_to_xyz, xyz_to_ijk,
            coords, i0, i1, r_grid, sq_rad);
    }
    for (auto& th: threads) th.join();
}

PYBIND11_MODULE(_map_mask, m){
    m.doc() = "Mask a map down to surround a set of coordinates.";
    m.def("generate_mask",
        [](py::array_t<uint8_t,0> map, c_array_in<double> origin,
            c_array_in<double> step, c_array_in<size_t> dim,
            c_array_in<double> ijk_to_xyz, c_array_in<double> xyz_to_ijk,
            c_array_in<double> coords, size_t n, double radius, size_t num_threads)
        {
            check_size(dim, 3, "dim");
            check_size(origin, 3, "origin");
            check_size(step, 3, "step");
            check_size(ijk_to_xyz, 12, "ijk_to_xyz");
            check_size(xyz_to_ijk, 12, "xyz_to_ijk");
            check_size(coords, 3*(ssize_t)n, "coords");
            auto dptr = static_cast<size_t*>(dim.request().ptr);
            check_mask_array(map, dptr);
            auto mptr = static_cast<uint8_t*>(map.request().ptr);
            auto optr = static_cast<double*>(origin.request().ptr);
            auto sptr = static_cast<double*>(step.request().ptr);
            auto itptr = static_cast<double*>(ijk_to_xyz.request().ptr);
            auto xtptr = static_cast<double*>(xyz_to_ijk.request().ptr);
            auto cptr = static_cast<double*>(coords.request().ptr);
            py::gil_scoped_release release;
            generate_mask(mptr, optr, sptr, dptr, itptr, xtptr, cptr, n, radius, num_threads);
        },
        py::arg("map"), py::arg("origin"), py::arg("step"), py::arg("dim"),
        py::arg("ijk_to_xyz"), py::arg("xyz_to_ijk"), py::arg("coords"),
        py::arg("n"), py::arg("radius"), py::arg("num_threads")=1);
}
