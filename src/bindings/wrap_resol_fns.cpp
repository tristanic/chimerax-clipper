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
 // #include <pybind11/operators.h>
 // #include <pybind11/numpy.h>


 #include <clipper/clipper.h>

 #include "numpy_helper.h"

 namespace py=pybind11;
 using namespace clipper;

// from resol_fn.h

void declare_basisfn_base(py::module& m)
{
    py::class_<BasisFn_base, std::unique_ptr<BasisFn_base, py::nodelete>> basisfn_base(m, "_BasisFn_base");

    py::enum_<BasisFn_base::FNtype>(basisfn_base, "FNtype")
        .value("GENERAL", BasisFn_base::FNtype::GENERAL)
        .value("LINEAR", BasisFn_base::FNtype::LINEAR)
        ;

    basisfn_base
        .def("f", &BasisFn_base::f)
        .def_property_readonly("num_params", &BasisFn_base::num_params)
        ;
}

void declare_basisfn_fderiv(py::module& m)
{
    using Class = BasisFn_base::Fderiv;
    py::class_<Class>(m, "BasisFn_Fderiv")
        .def_property_readonly("f", [](const Class& self) { return self.f; })
        .def_property_readonly("df", [](const Class& self) { return self.df; })
        ;
}

void declare_targetfn_base(py::module& m)
{
    py::class_<TargetFn_base> targetfn_base(m, "_TargetFn_base");

    py::enum_<TargetFn_base::FNtype>(targetfn_base, "FNtype")
    .value("GENERAL", TargetFn_base::FNtype::GENERAL)
    .value("QUADRATIC", TargetFn_base::FNtype::QUADRATIC)
    ;
}

void declare_resolution_fn(py::module& m)
{
    py::class_<ResolutionFn>(m, "ResolutionFn")
        /*
        As of 26/10/2021, binding of STL containers passed by const ref appears to be a no-no - 
        the temporary vector produced by the implicit conversion machinery is out of scope by the 
        time it makes it to the function. Replaced with a pass by value alternative.
        */
        // .def(py::init<const HKL_info&, const BasisFn_base&, const TargetFn_base&,
        //         const std::vector<ftype>&, const ftype, const bool>(),
        //     py::arg("hkl_info"), py::arg("basis_fn"), py::arg("target_fn"),
        //     py::arg("params"), py::arg("damp") = 0.0, py::arg("debug") = false)
        .def(py::init([](const HKL_info& hkl_info, const BasisFn_base& basis_fn, const TargetFn_base& target_fn,
            std::vector<ftype> params, ftype damp, bool debug) {
                return std::unique_ptr<ResolutionFn>(new ResolutionFn(hkl_info, basis_fn, target_fn,
                    params, damp, debug));
            }),
            py::arg("hkl_info"), py::arg("basis_fn"), py::arg("target_fn"),
            py::arg("params"), py::arg("damp") = 0.0, py::arg("debug") = false)
        .def("f", &ResolutionFn::f)
        .def_property_readonly("params", &ResolutionFn::params)
        ;
}

void declare_resolution_fn_nonlinear(py::module& m)
{
    py::class_<ResolutionFn_nonlinear, ResolutionFn>(m, "ResolutionFn_nonlinear")
        /*
        As of 26/10/2021, binding of STL containers passed by const ref appears to be a no-no - 
        the temporary vector produced by the implicit conversion machinery is out of scope by the 
        time it makes it to the function. Replaced with a pass by value alternative.
        */
        // .def(py::init<const HKL_info&, const BasisFn_base&, const TargetFn_base&,
        //     const std::vector<ftype>&, const ftype, const bool >(),
        //     py::arg("hkl_info"), py::arg("basis_fn"), py::arg("target_fn"),
        //     py::arg("params"), py::arg("damp") = 0.0, py::arg("debug") = false)
        .def(py::init([](const HKL_info& hkl_info, const BasisFn_base& basis_fn, const TargetFn_base& target_fn,
            std::vector<ftype> params, ftype damp, bool debug) {
                return std::unique_ptr<ResolutionFn_nonlinear>(new ResolutionFn_nonlinear(hkl_info, basis_fn, target_fn,
                    params, damp, debug));
            }),
            py::arg("hkl_info"), py::arg("basis_fn"), py::arg("target_fn"),
            py::arg("params"), py::arg("damp") = 0.0, py::arg("debug") = false)
        ;
}

// from resol_basisfn.h
void declare_resolution_ordinal(py::module& m)
{
    py::class_<Resolution_ordinal, Generic_ordinal>(m, "Resolution_ordinal")
        .def(py::init([](const HKL_info& hklinfo, const ftype& power) 
        {
            auto ret = new Resolution_ordinal();
            ret->init(hklinfo, power);
            return std::unique_ptr<Resolution_ordinal>(ret);
        }))
        .def(py::init([](const HKL_data_base& hkldata, const ftype& power) 
        {
            auto ret = new Resolution_ordinal();
            ret->init(hkldata, power);
            return std::unique_ptr<Resolution_ordinal>(ret);
        }))
        .def(py::init([](const HKL_data_base& hkldata, const Cell& cell, const ftype& power) 
        {
            auto ret = new Resolution_ordinal();
            ret->init(hkldata, cell, power);
            return std::unique_ptr<Resolution_ordinal>(ret);
        }))
        .def("init", (void (Resolution_ordinal::*)(const HKL_info&, const ftype&)) &Resolution_ordinal::init)
        .def("init", (void (Resolution_ordinal::*)(const HKL_data_base&, const ftype&)) &Resolution_ordinal::init)
        .def("init", (void (Resolution_ordinal::*)(const HKL_data_base&, const Cell&, const ftype&)) &Resolution_ordinal::init)
        ;
}

template <class C, class B>
void add_basisfn_base_class_functions(py::class_<C,B>& pyclass)
{
    pyclass
        .def("f", &C::f)
        .def("fderiv", &C::fderiv)
        .def_property_readonly("type", &C::type)
        .def_property_readonly("num_diagonals", &C::num_diagonals)
        ;
}

void declare_basisfn_binner(py::module& m)
{
    py::class_<BasisFn_binner, BasisFn_base> basis_binner(m, "BasisFn_binner");
    basis_binner
        .def(py::init<const HKL_info&, const int&, const ftype>(),
            py::arg("hkl_info"), py::arg("num_bins"), py::arg("power")=1.0)
        .def(py::init<const HKL_data_base&, const int&, const ftype>(),
            py::arg("hkl_data"), py::arg("num_bins"), py::arg("power")=1.0)
        .def("f_s", &BasisFn_binner::f_s)
        .def("fderiv_s", &BasisFn_binner::fderiv_s)
        ;
    add_basisfn_base_class_functions<BasisFn_binner, BasisFn_base>(basis_binner);
}

void declare_basisfn_linear(py::module& m)
{
    py::class_<BasisFn_linear, BasisFn_base> basis_linear(m, "BasisFn_linear");
    basis_linear
        .def(py::init<const HKL_info&, const int&, const ftype>(),
            py::arg("hkl_info"), py::arg("num_bins"), py::arg("power")=1.0)
        .def(py::init<const HKL_data_base&, const int&, const ftype>(),
            py::arg("hkl_data"), py::arg("num_bins"), py::arg("power")=1.0)
        .def("f_s", &BasisFn_linear::f_s)
        .def("fderiv_s", &BasisFn_linear::fderiv_s)
        ;
    add_basisfn_base_class_functions<BasisFn_linear, BasisFn_base>(basis_linear);
}

void declare_basisfn_spline(py::module& m)
{
    py::class_<BasisFn_spline, BasisFn_base> basis_spline(m, "BasisFn_spline");
    basis_spline
        .def(py::init<const HKL_info&, const int&, const ftype>(),
            py::arg("hkl_info"), py::arg("num_bins"), py::arg("power")=1.0)
        .def(py::init<const HKL_data_base&, const int&, const ftype>(),
            py::arg("hkl_data"), py::arg("num_bins"), py::arg("power")=1.0)
        .def("f_s", &BasisFn_spline::f_s)
        .def("fderiv_s", &BasisFn_spline::fderiv_s)
        ;
    add_basisfn_base_class_functions<BasisFn_spline, BasisFn_base>(basis_spline);
}

void declare_basisfn_gaussian(py::module& m)
{
    py::class_<BasisFn_gaussian, BasisFn_base>(m, "BasisFn_gaussian")
        .def(py::init<>())
        .def("fderiv_s", &BasisFn_gaussian::fderiv_s)
        .def("fderiv", &BasisFn_gaussian::fderiv)
        .def("scale", &BasisFn_gaussian::scale)
        .def("u_iso", &BasisFn_gaussian::u_iso)
        ;
}

void declare_basisfn_aniso_gaussian(py::module& m)
{
    py::class_<BasisFn_aniso_gaussian, BasisFn_base>(m, "BasisFn_aniso_gaussian")
        .def(py::init<>())
        .def("fderiv_coord", &BasisFn_aniso_gaussian::fderiv_coord)
        .def("fderiv", &BasisFn_aniso_gaussian::fderiv)
        .def("scale", &BasisFn_aniso_gaussian::scale)
        .def("u_aniso_orth", &BasisFn_aniso_gaussian::u_aniso_orth)
        ;
}

void declare_basisfn_log_gaussian(py::module& m)
{
    py::class_<BasisFn_log_gaussian, BasisFn_base>(m, "BasisFn_log_gaussian")
        .def(py::init<>())
        .def("fderiv_s", &BasisFn_log_gaussian::fderiv_s)
        .def("fderiv", &BasisFn_log_gaussian::fderiv)
        .def_property_readonly("type", &BasisFn_log_gaussian::type)
        .def("scale", &BasisFn_log_gaussian::scale)
        .def("u_iso", &BasisFn_log_gaussian::u_iso)
        ;
}

void declare_basisfn_log_aniso_gaussian(py::module& m)
{
    py::class_<BasisFn_log_aniso_gaussian, BasisFn_base>(m, "BasisFn_log_aniso_gaussian")
        .def(py::init<>())
        .def("fderiv_coord", &BasisFn_log_aniso_gaussian::fderiv_coord)
        .def("fderiv", &BasisFn_log_aniso_gaussian::fderiv)
        .def_property_readonly("type", &BasisFn_log_aniso_gaussian::type)
        .def("scale", &BasisFn_log_aniso_gaussian::scale)
        .def("u_aniso_orth", &BasisFn_log_aniso_gaussian::u_aniso_orth)
        ;
}

// From resol_targetfn.h

template<class T>
void declare_targetfn_meanfnth(py::module& m, const std::string& suffix)
{
    std::string pyname = std::string("TargetFn_meanFnth_") + suffix;
    using Class = TargetFn_meanFnth<T>;
    py::class_<Class, TargetFn_base>(m, pyname.c_str())
        .def(py::init<const HKL_data<T>&, const ftype&>())
        .def("rderiv", &Class::rderiv)
        .def_property_readonly("type", &Class::type)
        ;
}

template<class T>
void declare_targetfn_meaninth(py::module& m, const std::string& suffix)
{
    std::string pyname = std::string("TargetFn_meanInth_") + suffix;
    using Class = TargetFn_meanInth<T>;
    py::class_<Class, TargetFn_base>(m, pyname.c_str())
        .def(py::init<const HKL_data<T>&, const ftype&>())
        .def("rderiv", &Class::rderiv)
        .def_property_readonly("type", &Class::type)
        ;
}

template<class T>
void declare_targetfn_meanenth(py::module& m, const std::string& suffix)
{
    std::string pyname = std::string("TargetFn_meanEnth_") + suffix;
    using Class = TargetFn_meanEnth<T>;
    py::class_<Class, TargetFn_base>(m, pyname.c_str())
        .def(py::init<const HKL_data<T>&, const ftype&>())
        .def("rderiv", &Class::rderiv)
        .def_property_readonly("type", &Class::type)
        ;
}

template<class T1, class T2>
void declare_targetfn_scalef1f2(py::module&m, const std::string& suffix)
{
    std::string pyname = std::string("TargetFn_scaleF1F2_") + suffix;
    using Class = TargetFn_scaleF1F2<T1,T2>;
    py::class_<Class, TargetFn_base>(m, pyname.c_str())
        .def(py::init<const HKL_data<T1>&, const HKL_data<T2>&>())
        .def("rderiv", &Class::rderiv)
        .def_property_readonly("type", &Class::type)
        ;
}

template<class T1, class T2>
void declare_targetfn_scalelogf1f2(py::module&m, const std::string& suffix)
{
    std::string pyname = std::string("TargetFn_scaleLogF1F2_") + suffix;
    using Class = TargetFn_scaleLogF1F2<T1,T2>;
    py::class_<Class, TargetFn_base>(m, pyname.c_str())
        .def(py::init<const HKL_data<T1>&, const HKL_data<T2>&>())
        .def("rderiv", &Class::rderiv)
        .def_property_readonly("type", &Class::type)
        ;
}

template<class T1, class T2>
void declare_targetfn_scalei1i2(py::module&m, const std::string& suffix)
{
    std::string pyname = std::string("TargetFn_scaleI1I2_") + suffix;
    using Class = TargetFn_scaleI1I2<T1,T2>;
    py::class_<Class, TargetFn_base>(m, pyname.c_str())
        .def(py::init<const HKL_data<T1>&, const HKL_data<T2>&>())
        .def("rderiv", &Class::rderiv)
        .def_property_readonly("type", &Class::type)
        ;
}

template<class T1, class T2>
void declare_targetfn_scalelogi1i2(py::module&m, const std::string& suffix)
{
    std::string pyname = std::string("TargetFn_scaleLogI1I2_") + suffix;
    using Class = TargetFn_scaleLogI1I2<T1,T2>;
    py::class_<Class, TargetFn_base>(m, pyname.c_str())
        .def(py::init<const HKL_data<T1>&, const HKL_data<T2>&>())
        .def("rderiv", &Class::rderiv)
        .def_property_readonly("type", &Class::type)
        ;
}

template<class T>
void declare_targetfn_scaleesq(py::module& m, const std::string& suffix)
{
    std::string pyname = std::string("TargetFn_scaleEsq_") + suffix;
    using Class = TargetFn_scaleEsq<T>;
    py::class_<Class, TargetFn_base>(m, pyname.c_str())
        .def(py::init<const HKL_data<T>&>())
        .def("rderiv", &Class::rderiv)
        .def_property_readonly("type", &Class::type)
        ;
}

template<class T>
void declare_targetfn_sigmaa_omegaa(py::module& m, const std::string& suffix)
{
    std::string pyname = std::string("TargetFn_sigmaa_omegaa_") + suffix;
    using Class = TargetFn_sigmaa_omegaa<T>;
    py::class_<Class, TargetFn_base>(m, pyname.c_str())
        .def(py::init<const HKL_data<T>&, const HKL_data<T>&>())
        .def("rderiv", &Class::rderiv)
        .def_static("sigmaa", &Class::sigmaa)
        ;
}



void init_resol_fn(py::module& m)
{
    declare_basisfn_base(m);
    declare_basisfn_fderiv(m);
    declare_targetfn_base(m);
    declare_resolution_fn(m);
    declare_resolution_fn_nonlinear(m);
    declare_resolution_ordinal(m);
    declare_basisfn_binner(m);
    declare_basisfn_linear(m);
    declare_basisfn_spline(m);
    declare_basisfn_gaussian(m);
    declare_basisfn_aniso_gaussian(m);
    declare_basisfn_log_gaussian(m);
    declare_basisfn_log_aniso_gaussian(m);

    declare_targetfn_meanfnth<data32::F_sigF>(m, "F_sigF_float");
    declare_targetfn_meanfnth<data32::F_sigF_ano>(m, "F_sigF_ano_float");
    declare_targetfn_meanfnth<data32::F_phi>(m, "F_phi_float");
    declare_targetfn_meaninth<data32::I_sigI>(m, "I_sigI_float");
    declare_targetfn_meaninth<data32::I_sigI_ano>(m, "I_sigI_ano_float");
    declare_targetfn_meanenth<data32::E_sigE>(m, "E_sigE_float");

    declare_targetfn_scalef1f2<data32::F_sigF,data32::F_sigF>(m, "F_sigF_float");
    declare_targetfn_scalelogf1f2<data32::F_sigF,data32::F_sigF>(m, "F_sigF_float");
    declare_targetfn_scalei1i2<data32::I_sigI,data32::I_sigI>(m, "I_sigI_float");
    declare_targetfn_scalelogi1i2<data32::I_sigI,data32::I_sigI>(m, "I_sigI_float");

    declare_targetfn_scaleesq<data32::E_sigE>(m, "E_sigE_float");

    declare_targetfn_sigmaa_omegaa<data32::E_sigE>(m, "E_sigE_float");



    declare_targetfn_meanfnth<data64::F_sigF>(m, "F_sigF_double");
    declare_targetfn_meanfnth<data64::F_sigF_ano>(m, "F_sigF_ano_double");
    declare_targetfn_meanfnth<data64::F_phi>(m, "F_phi_double");
    declare_targetfn_meaninth<data64::I_sigI>(m, "I_sigI_double");
    declare_targetfn_meaninth<data64::I_sigI_ano>(m, "I_sigI_ano_double");
    declare_targetfn_meanenth<data64::E_sigE>(m, "E_sigE_double");

    declare_targetfn_scalef1f2<data64::F_sigF,data64::F_sigF>(m, "F_sigF_double");
    declare_targetfn_scalelogf1f2<data64::F_sigF,data64::F_sigF>(m, "F_sigF_double");
    declare_targetfn_scalei1i2<data64::I_sigI,data64::I_sigI>(m, "I_sigI_double");
    declare_targetfn_scalelogi1i2<data64::I_sigI,data64::I_sigI>(m, "I_sigI_double");

    declare_targetfn_scaleesq<data64::E_sigE>(m, "E_sigE_double");

    declare_targetfn_sigmaa_omegaa<data64::E_sigE>(m, "E_sigE_double");


}
