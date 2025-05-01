#pragma once

#include <clipper/clipper.h>

using namespace clipper;
namespace clipper_cx
{

//template <typename T>
struct Packaged_ResolutionFn
{
    Packaged_ResolutionFn () {}
    //std::unique_ptr<TargetFn_scaleF1F2<f1, f2>> targetfn;
    std::unique_ptr<TargetFn_base> targetfn;
    BasisFn_aniso_gaussian aniso_basisfn;
    std::unique_ptr<ResolutionFn_nonlinear> resolution_fn;
};

template <typename T>
Packaged_ResolutionFn aniso_scale_fn(const HKL_data<F_phi<T>>& fcalc,
    const HKL_data<F_sigF<T>>& fobs, U_aniso_orth& uaniso, std::vector<ftype>& aniso_params)
{
    const auto& hkls = fcalc.hkl_info();
    Packaged_ResolutionFn ret;
    //ret.targetfn = std::unique_ptr<TargetFn_scaleF1F2<f1, f2>>(new TargetFn_scaleF1F2<f1, f2>(farr1, farr2));
    ret.targetfn = std::unique_ptr<TargetFn_base>(new TargetFn_scaleFobsFcalc<T>(fobs, fcalc));
    // First fit an isotropic Gaussian to improve stability of aniso step
    // BasisFn_gaussian iso_basisfn;
    // std::vector<ftype> iso_params = {1.0, 1.0};
    // ResolutionFn_nonlinear iso_rfn(hkls, iso_basisfn, *(ret.targetfn), iso_params);
    // iso_params = iso_rfn.params();
    // std::vector<ftype> aniso_params = {iso_params[0], iso_params[1], iso_params[1], iso_params[1], 0.0, 0.0, 0.0};
    ret.resolution_fn = std::unique_ptr<ResolutionFn_nonlinear>(
        new ResolutionFn_nonlinear(hkls, ret.aniso_basisfn, *(ret.targetfn), aniso_params, 10.0));
    aniso_params = ret.resolution_fn->params();
    uaniso = ret.aniso_basisfn.u_aniso_orth(ret.resolution_fn->params());
    return ret;
}

// Scale Fcalc to Fobs using an anisotropic Gaussian and then "polishing"
// with an isotropic spline
template <typename T>
void scale_fcalc_to_fobs(const HKL_data<F_phi<T>>& fcalc,
    const HKL_data<F_sigF<T>>& fobs, HKL_data<F_phi<T>>& scaled_fcalc,
    U_aniso_orth& uaniso, std::vector<ftype>& aniso_params, size_t n_spline_params = 20)
{
    const auto& hkls = fcalc.hkl_info();
    auto rfn_p = aniso_scale_fn(fcalc, fobs, uaniso, aniso_params);
    const auto& aniso_rfn = *(rfn_p.resolution_fn);
    for (auto ih = fobs.first(); !ih.last(); ih.next())
    {
        if (!fobs[ih].missing() && !fcalc[ih].missing())
        {
            scaled_fcalc[ih].f() = fcalc[ih].f()/sqrt(aniso_rfn.f(ih));
            scaled_fcalc[ih].phi() = fcalc[ih].phi();
        }
    }
    std::vector<ftype> iso_params(n_spline_params, 1.0);
    BasisFn_spline iso_basisfn(fobs, n_spline_params, 1.0);
    TargetFn_scaleF1F2<F_phi<T>, F_sigF<T>> iso_target(scaled_fcalc, fobs);
    ResolutionFn iso_rfn(hkls, iso_basisfn, iso_target, iso_params);
    for (auto ih=scaled_fcalc.first(); !ih.last(); ih.next())
    {
        if (!scaled_fcalc[ih].missing())
        {
            scaled_fcalc[ih].f() = scaled_fcalc[ih].f() * sqrt(iso_rfn.f(ih));
        }
    }
}

// Scale Fobs to Fcalc using an anisotropic Gaussian only
template <typename T>
void aniso_scale_fobs_to_fcalc(const HKL_data<F_phi<T>>& fcalc,
    const HKL_data<F_sigF<T>>& fobs, HKL_data<F_sigF<T>>& scaled_fobs,
    U_aniso_orth& uaniso, std::vector<ftype>& aniso_params)
{
    auto rfn_p = aniso_scale_fn(fcalc, fobs, uaniso, aniso_params);
    const auto& aniso_rfn = *(rfn_p.resolution_fn);
    for (auto ih = fobs.first(); !ih.last(); ih.next())
    {
        if (!fobs[ih].missing() && !fcalc[ih].missing())
        {
            auto scale = sqrt(aniso_rfn.f(ih));
            scaled_fobs[ih].f() = fobs[ih].f()*scale;
            scaled_fobs[ih].sigf() = fobs[ih].sigf()*scale;
        }
    }
    // std::cout << "Anisotropic scale params: ";
    // for (auto p: aniso_params)
    //     std::cout << p << ", ";
    // std::cout << std::endl;
}



} // namespace clipper_cx
