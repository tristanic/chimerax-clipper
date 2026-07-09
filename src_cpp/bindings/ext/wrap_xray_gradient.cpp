/**
 * @Author: Tristan Croll <tcroll@altoslabs.com>
 * Pybind11 bindings for XrayGradientEvaluator — synchronous crystallographic
 * value-and-gradient evaluation for end-to-end ML force-field training.
 *
 * The interface is numpy-only: it returns residual + per-atom gradients as a
 * plain float64 array for a PyTorch training loop, with NO torch dependency on
 * the C++ side (the thin autograd adapter lives in chimerax.clipper.diff).
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <clipper/clipper.h>
#include <clipper_ext/xray_gradient.h>

namespace py = pybind11;
namespace cx = clipper_cx;
using namespace clipper;
using namespace clipper::datatypes;

void init_xray_gradient(py::module& m)
{
    py::enum_<cx::XrayTargetKind>(m, "XrayTargetKind",
        "Reciprocal-space least-squares target form.\n"
        "  AmplitudeLS : 1/2 sum w (k|Fc| - m|Fo|)^2  (m = FOM; set m=1 for small\n"
        "                molecules via a phi_fom whose fom()=1; w = 1/sigma^2(Fo)).\n"
        "  IntensityLS : 1/2 sum w (Io - s|Fc|^2)^2, w = 1/sigma^2(Io) (SHELX-style;\n"
        "                Io = |Fo|^2, sigma(Io) ~ 2|Fo|sigma(F) when only F/sigF given).")
        .value("AmplitudeLS", cx::XrayTargetKind::AmplitudeLS)
        .value("IntensityLS", cx::XrayTargetKind::IntensityLS);

    py::class_<cx::XrayGradientEvaluator>(m, "XrayGradientEvaluator",
        "Synchronous crystallographic value-and-gradient evaluator.\n"
        "\n"
        "Seeds the fixed reference data (element list + observed structure factors,\n"
        "or a fixed target map) once at construction, then answers repeated\n"
        "value_and_gradient() calls that drive the mutable per-atom parameters from\n"
        "numpy arrays.  No L-BFGS loop, no background thread, no torch dependency.\n"
        "Gradients are returned for a runtime-selected subset of the eleven per-atom\n"
        "parameters {X,Y,Z,Uiso,Occ,U11,U22,U33,U12,U13,U23}, including coordinates.")

        // Reciprocal-space constructor (least-squares vs observed Fo / Io).
        .def(py::init([](std::vector<std::string>            elements,
                         const HKL_data<F_sigF<ftype32>>&    fobs,
                         const HKL_data<Phi_fom<ftype32>>&   phi_fom,
                         const HKL_data<Flag>&               usage,
                         cx::XrayTargetKind                  kind,
                         const HKL_data<F_phi<ftype32>>&     f_bulk,
                         int                                 n_threads,
                         bool                                threaded_density)
            {
                std::vector<String> els(elements.begin(), elements.end());
                return new cx::XrayGradientEvaluator(els, fobs, phi_fom, usage,
                                                     kind, f_bulk, n_threads,
                                                     threaded_density);
            }),
            py::arg("elements"),
            py::arg("fobs"),
            py::arg("phi_fom"),
            py::arg("usage"),
            py::arg("kind") = cx::XrayTargetKind::AmplitudeLS,
            py::arg("f_bulk") = HKL_data<F_phi<ftype32>>(),
            py::arg("n_threads") = 1,
            py::arg("threaded_density") = true)

        // Real-space constructor (least-squares vs a fixed target Xmap).
        .def(py::init([](std::vector<std::string> elements,
                         const Xmap<ftype32>&     target_map,
                         const Coord_orth&        target_origin,
                         int                      n_threads,
                         bool                     threaded_density)
            {
                std::vector<String> els(elements.begin(), elements.end());
                return new cx::XrayGradientEvaluator(els, target_map,
                                                     target_origin, n_threads,
                                                     threaded_density);
            }),
            py::arg("elements"),
            py::arg("target_map"),
            py::arg("target_origin"),
            py::arg("n_threads") = 1,
            py::arg("threaded_density") = true)

        .def_property_readonly("n_atoms", &cx::XrayGradientEvaluator::n_atoms)

        // One value + gradient evaluation.  Arrays are cast to contiguous float64
        // (uint8 for is_aniso).  `selected` is a list of AtomShapeFn.TYPE enum
        // values (their ints); the returned gradient columns follow that order.
        // Returns (float target, ndarray[N, P]).
        .def("value_and_gradient",
            [](cx::XrayGradientEvaluator& self,
               py::array_t<double,  py::array::c_style | py::array::forcecast> coords,
               py::array_t<double,  py::array::c_style | py::array::forcecast> u_iso,
               py::array_t<double,  py::array::c_style | py::array::forcecast> u_aniso,
               py::array_t<double,  py::array::c_style | py::array::forcecast> occ,
               py::array_t<uint8_t, py::array::c_style | py::array::forcecast> is_aniso,
               std::vector<int>                                                selected,
               bool                                                            refresh_scale) -> py::tuple
            {
                const int N = self.n_atoms();
                const int P = (int)selected.size();
                if (P == 0)
                    throw std::invalid_argument("value_and_gradient: empty parameter selection");
                if ((int)coords.size()   != 3*N) throw std::invalid_argument("coords must have shape (N,3)");
                if ((int)u_iso.size()    != N)   throw std::invalid_argument("u_iso must have length N");
                if ((int)u_aniso.size()  != 6*N) throw std::invalid_argument("u_aniso must have shape (N,6)");
                if ((int)occ.size()      != N)   throw std::invalid_argument("occ must have length N");
                if ((int)is_aniso.size() != N)   throw std::invalid_argument("is_aniso must have length N");

                std::vector<AtomShapeFn::TYPE> sel;
                sel.reserve(P);
                for (int i = 0; i < P; ++i)
                    sel.push_back((AtomShapeFn::TYPE)selected[i]);

                py::array_t<double> grad({N, P});
                double T = self.value_and_gradient(
                    coords.data(), u_iso.data(), u_aniso.data(), occ.data(),
                    is_aniso.data(), sel,
                    static_cast<double*>(grad.request().ptr), refresh_scale);
                return py::make_tuple(T, grad);
            },
            py::arg("coords"),
            py::arg("u_iso"),
            py::arg("u_aniso"),
            py::arg("occ"),
            py::arg("is_aniso"),
            py::arg("selected"),
            py::arg("refresh_scale") = false);
}
