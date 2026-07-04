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

#include <memory>
#include <future>
#include <iostream>

// #include "bindings/numpy_helper.h"

#include "contour.h"

using namespace Contour_Calculation;
namespace py=pybind11;

struct Contour_Geometry
{
public:
    std::vector<float>* vertex_xyz;
    std::vector<float>* normals;
    std::vector<int>* tv_indices;
    int vertex_count=0;
    int triangle_count=0;
    Contour_Geometry() {}
    Contour_Geometry(int vc, int tc)
        : vertex_count(vc), triangle_count(tc)
    {
        vertex_xyz = new std::vector<float>(vc*3);
        normals = new std::vector<float>(vc*3);
        tv_indices = new std::vector<int>(tc*3);
    }
};

void copy_transform(py::array_t<float> source, float dest[3][4])
{
    for (size_t i=0; i<3; ++i)
        for (size_t j=0; j<4; ++j)
            dest[i][j] = *(source.data(i,j));
}


/* HELPER FUNCTIONS */

template <typename T>
void transform_coord (T tf[3][4], T coord[3], T out[3])
{
    for (size_t i=0; i<3; ++i) {
        T* row = tf[i];
        out[i] = row[0]*coord[0] + row[1]*coord[1] + row[2]*coord[2] + row[3];
    }
}

// Transform coordinates in-place
template <typename T>
void transform_coords (T tf[3][4], T* coords, int n)
{
    T temp[3];
    for (int i=0; i<n; ++i)
    {
        transform_coord(tf, coords, temp);
        for (int j=0; j<3; ++j)
            *coords++ = temp[j];
    }
}

template <typename T>
T l2_norm_3d(T a[3])
{
    T accum = 0;
    for (int i = 0; i < 3; i++) {
        accum += a[i]*a[i];
    }
    return sqrt(accum);
}

// Normalize 3D vectors in-place
template <typename T>
void normalize_3d_vectors(T* vecs, int n)
{
    for (int i=0; i<n; ++i)
    {
        T norm = l2_norm_3d(vecs);
        // Avoid divide-by-zero
        norm = norm > 0 ? norm : 1;
        for (int j=0; j<3; ++j)
            *vecs++ /= norm;
    }
}

// ----------------------------------------------------------------------------
// Swap vertex 1 and 2 of each triangle.
//
void reverse_triangle_vertex_order(int* ta, int n)
{
  int s0=3, s1 = 1;
  // int s0 = triangles.stride(0), s1 = triangles.stride(1);
  for (int t = 0 ; t < n ; ++t)
    {
      int i1 = s0*t+s1, i2 = i1 + s1;
      int v1 = ta[i1], v2 = ta[i2];
      ta[i1] = v2;
      ta[i2] = v1;
    }
}

class Contour_Thread_Mgr
{
public:
    Contour_Thread_Mgr() {}
    void start_compute(py::array_t<float> data, float threshold, float det,
        py::array_t<float> vertex_transform, py::array_t<float> normal_transform,
        bool cap_faces=true, bool return_normals=false);
    bool ready() const { return ready_; }
    bool return_normals() const { return return_normals_; }
    Contour_Geometry get_result();
private:
    // We do NOT copy the volume data. Instead we keep a reference to the source
    // numpy array (co-ownership) and let the worker read it directly via data_ptr_.
    // This keeps the buffer alive even if the map is closed mid-contour (the map's
    // own reference can drop; ours holds it). The std::future's destructor blocks
    // until the worker finishes, and data_ref_ is declared before geom_ so it is
    // destroyed AFTER geom_ - i.e. it always outlives the worker's reads.
    // Tradeoff (option 2): if the box is refilled in place while a worker is
    // running (fast spotlight scroll, or unusually fast structure-factor recalc),
    // that surface can transiently "tear" for one frame; it self-corrects on the
    // next contour. No copy, no crash.
    py::array_t<float> data_ref_;
    float* data_ptr_ = nullptr;
    float threshold_;
    bool flip_triangles_; // if det < 0
    float vertex_transform_[3][4]; // Transform mapping vertices into model coordinates
    float normal_transform_[3][4]; // Transform mapping normals into model coordinates

    Stride stride_[3];
    Index size_[3];
    std::future<Contour_Geometry> geom_;
    bool working_=false;
    bool ready_=false;
    bool return_normals_=false;
    Contour_Geometry contour_surface_thread_(float data[], float threshold, bool cap_faces);
};

void Contour_Thread_Mgr::start_compute(py::array_t<float> data, float threshold,
    float det, py::array_t<float> vertex_transform, py::array_t<float> normal_transform,
    bool cap_faces, bool return_normals)
{
    if (working_)
        throw std::runtime_error("Contour thread is already running!");
    threshold_=threshold;
    // Reference (do NOT copy) the source data; the worker reads it directly via
    // the raw pointer, and the strides tell surface() how to index it (the map's
    // native (z,y,x) buffer presented as an (x,y,z) transposed view). Called with
    // the GIL held, so touching the numpy array here is safe.
    int stride_size = sizeof(float);
    for (size_t i=0; i<3; ++i)
    {
        // Numpy strides are per byte; surface() wants strides per float.
        stride_[i] = data.strides(i)/stride_size;
        size_[i] = data.shape(i);
    }
    data_ref_ = data;                             // co-own: keeps the buffer alive
    data_ptr_ = (float*)data_ref_.request().ptr;  // underlying (contiguous) buffer
    flip_triangles_ = det < 0;
    copy_transform(vertex_transform, vertex_transform_);
    copy_transform(normal_transform, normal_transform_);
    working_=true;
    return_normals_ = return_normals;
    ready_=false;
    try {
        geom_ = std::async(std::launch::async, &Contour_Thread_Mgr::contour_surface_thread_, this, data_ptr_, threshold, cap_faces);
    } catch (...) {
        working_ = false;
        throw;
    }
}

Contour_Geometry Contour_Thread_Mgr::contour_surface_thread_(float* data, float threshold, bool cap_faces)
{
    auto cptr = std::unique_ptr<Contour_Surface>(surface(data, size_, stride_, threshold, cap_faces));
    int vc = cptr->vertex_count();
    int tc = cptr->triangle_count();
    Contour_Geometry geom(vc, tc);
    cptr->geometry(geom.vertex_xyz->data(), reinterpret_cast<Index *>(geom.tv_indices->data()));
    if (return_normals_)
        cptr->normals(geom.normals->data());
    if (flip_triangles_)
        reverse_triangle_vertex_order(geom.tv_indices->data(), tc);

    // Transform coords and normals to model coordinates
    transform_coords(vertex_transform_, geom.vertex_xyz->data(), vc);
    if (return_normals_)
    {
        transform_coords(normal_transform_, geom.normals->data(), vc);
        // Set all normals to unit length;
        normalize_3d_vectors(geom.normals->data(), vc);

    }
    ready_=true;
    return geom;
}

Contour_Geometry
Contour_Thread_Mgr::get_result()
{
    if (!working_)
        throw std::runtime_error("No contour calculation has been started!");
    working_ = false;
    auto geom = geom_.get();   // blocks until the worker has finished reading data_ptr_
    // Release our reference to the source buffer now the worker is done. Called on
    // the GUI thread with the GIL held, so decref-ing the numpy array is safe.
    data_ref_ = py::array_t<float>();
    data_ptr_ = nullptr;
    return geom;
}

template<typename T>
void delete_when_done(void *data)
{
    std::vector<T>* d = reinterpret_cast<std::vector<T> *>(data);
    delete d;
}



PYBIND11_MODULE(contour_thread, m) {
    m.doc() = "Threaded contouring implementation";


    py::class_<Contour_Thread_Mgr>(m, "Contour_Thread_Mgr")
        .def(py::init<>())
        .def("start_compute", &Contour_Thread_Mgr::start_compute,
            py::arg("data"), py::arg("threshold"), py::arg("determinant"),
            py::arg("vertex_transform"), py::arg("normal_transform"),
            py::arg("cap_faces")=true,
            py::arg("return_normals")=false)
        .def("ready", &Contour_Thread_Mgr::ready)
        .def("get_result", [](Contour_Thread_Mgr& self) -> py::tuple
        {
            // Takes ownership of vertex_xyz, tv_indices and normals. They will
            // be automatically deleted when they go out of scope in Python.
            auto geom = self.get_result();
            const auto& vc = geom.vertex_count;
            const auto& tc = geom.triangle_count;
            py::array_t<float> va({vc, 3}, // shape
                 {3*4, 4}, // C-style contiguous strides for n*3 float
                 geom.vertex_xyz->data(),
                 py::capsule(geom.vertex_xyz, &delete_when_done<float>));
            py::array_t<int> ta({tc, 3}, // shape
                {3*4, 4}, // C-style contiguous strides for n*3 int32
                geom.tv_indices->data(),
                py::capsule(geom.tv_indices, &delete_when_done<int>));
            if (!self.return_normals())
                return py::make_tuple(va, ta);
            py::array_t<float> na({vc, 3}, // shape
                {3*4, 4}, // C-style contiguous strides for n*3 float
                geom.normals->data(),
                py::capsule(geom.normals, &delete_when_done<float>));
            return py::make_tuple(va, ta, na);
        })
        ;


}
