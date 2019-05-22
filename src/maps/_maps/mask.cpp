/**
 * @Author: Tristan Croll <tic20>
 * @Date:   22-May-2019
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 22-May-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>

namespace py=pybind11;

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


void generate_mask(
    uint8_t* map, double* origin, double* step, size_t* dim,
    double* ijk_to_xyz, double* xyz_to_ijk, double* coords, size_t n, double radius)
{
    size_t box_min[3], box_max[3];
    double transformed[3], bc[3];
    double r_xyz[3], r_grid[3];
    for (size_t i=0; i<3; ++i)
        r_xyz[i] = radius + origin[i];
    affine_transform(r_xyz, xyz_to_ijk, r_grid);
    double sq_rad = pow(radius/step[0],2); //radius*radius;
    for (size_t i=0; i<n; ++i)
    {
        affine_transform(coords+3*i, xyz_to_ijk, transformed);
        index_range(transformed, r_grid, step, dim, box_min, box_max);
        for (size_t u=box_min[0]; u<=box_max[0]; ++u) {
            for (size_t v=box_min[1]; v<= box_max[1]; ++v) {
                for (size_t w=box_min[2]; w<= box_max[2]; ++w) {
                    bc[0]=u; bc[1]=v; bc[2]=w;
                    if (squared_distance(transformed, bc, ijk_to_xyz) < sq_rad)
                    {
                        map[ w + dim[2] * ( v + dim[1] * u ) ] = 1;
                    }
                }
            }
        }
    }
}

PYBIND11_MODULE(_map_mask, m){
    m.doc() = "Mask a map down to surround a set of coordinates.";
    m.def("generate_mask",
        [](py::array_t<uint8_t,0> map, py::array_t<double> origin,
            py::array_t<double> step, py::array_t<size_t> dim,
            py::array_t<double> ijk_to_xyz, py::array_t<double> xyz_to_ijk,
            py::array_t<double> coords, size_t n, double radius)
        {
            auto mptr = static_cast<uint8_t*>(map.request().ptr);
            auto optr = static_cast<double*>(origin.request().ptr);
            auto sptr = static_cast<double*>(step.request().ptr);
            auto dptr = static_cast<size_t*>(dim.request().ptr);
            auto itptr = static_cast<double*>(ijk_to_xyz.request().ptr);
            auto xtptr = static_cast<double*>(xyz_to_ijk.request().ptr);
            auto cptr = static_cast<double*>(coords.request().ptr);
            generate_mask(mptr, optr, sptr, dptr, itptr, xtptr, cptr, n, radius);
        });
}
