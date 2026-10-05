#pragma once

#include <clipper/clipper.h>

using namespace clipper;
namespace clipper_cx
{

//template <typename T>
struct Packaged_ResolutionFn
{
    Packaged_ResolutionFn () {}
    std::unique_ptr<TargetFn_base> targetfn;
    BasisFn_log_aniso_gaussian aniso_basisfn;
    std::unique_ptr<ResolutionFn> resolution_fn;
};

// Anisotropic scaling of Fcalc onto Fobs, fitted in LOG space.  The log-Gaussian
// basis + log target make this a linear least-squares problem, solved in a single
// closed-form pass by a plain ResolutionFn starting from zeros — stable across
// diverse datasets, with no isotropic pre-fit, no nonlinear convergence to babysit
// and no non-positive-definite failure mode (cf. Clipper's own SFscale_aniso).
// The fitted basis value fh = resolution_fn->f(ih) ~ log(Io/Ic), so the amplitude
// scale to bring Fc -> Fo is exp(0.5*fh) (and Fo -> Fc is exp(-0.5*fh)).
template <typename T>
Packaged_ResolutionFn aniso_scale_fn(const HKL_data<F_phi<T>>& fcalc,
    const HKL_data<F_sigF<T>>& fobs, U_aniso_orth& uaniso, std::vector<ftype>& aniso_params)
{
    const auto& hkls = fcalc.hkl_info();
    Packaged_ResolutionFn ret;
    ret.targetfn = std::unique_ptr<TargetFn_base>(
        new TargetFn_scaleLogF1F2<F_phi<T>, F_sigF<T>>(fcalc, fobs));
    aniso_params = std::vector<ftype>(7, 0.0);
    ret.resolution_fn = std::unique_ptr<ResolutionFn>(
        new ResolutionFn(hkls, ret.aniso_basisfn, *(ret.targetfn), aniso_params));
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
            // fh = aniso_rfn.f(ih) ~ log(Io/Ic); amplitude scale Fc -> Fo = exp(0.5*fh).
            scaled_fcalc[ih].f() = fcalc[ih].f() * exp(0.5 * aniso_rfn.f(ih));
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

// Isotropic-spline scaling target of |Fcalc|^2 onto a measured intensity I, where I
// keeps its sign. Same form as Clipper's TargetFn_scaleF1F2 (whose second set is
// |Fo|^2 and so can never be negative), so an intensity-space fit is not biased up
// by clamping the weak, negative measurements to zero.
template<class T1, class T2> class TargetFn_scaleF1I2 : public TargetFn_base
{
public:
    TargetFn_scaleF1I2(const HKL_data<T1>& hkl_data1_, const HKL_data<T2>& hkl_data2_)
        : hkl_data1(&hkl_data1_), hkl_data2(&hkl_data2_) {}
    Rderiv rderiv(const HKL_info::HKL_reference_index& ih, const ftype& fh) const
    {
        Rderiv result;
        const T1& ft1 = (*hkl_data1)[ih];
        const T2& ft2 = (*hkl_data2)[ih];
        if (!ft1.missing() && !ft2.missing()) {
            const ftype eps = ih.hkl_class().epsilon();
            const ftype f1 = pow(ft1.f(), 2) / eps;
            const ftype f2 = ft2.I() / eps;
            const ftype d = fh*f1 - f2;
            result.r = d * d / f1;
            result.dr = 2.0 * d;
            result.dr2 = 2.0 * f1;
        } else {
            result.r = result.dr = result.dr2 = 0.0;
        }
        return result;
    }
    FNtype type() const { return QUADRATIC; }
private:
    const HKL_data<T1>* hkl_data1;
    const HKL_data<T2>* hkl_data2;
};

// scale_fcalc_to_fobs for a target fed intensities directly. The anisotropic log-space
// stage is unchanged (on fobs it already uses only the Fo > 0, i.e. I > 0, rows, the
// only ones with a logarithm); the isotropic spline is then fitted against the signed
// I. Where the spline falls to or below zero (a shell with no net signal) the scale
// there is zero, rather than NaN.
template <typename T>
void scale_fcalc_to_fobs_iobs(const HKL_data<F_phi<T>>& fcalc,
    const HKL_data<F_sigF<T>>& fobs, const HKL_data<I_sigI<T>>& iobs,
    HKL_data<F_phi<T>>& scaled_fcalc, U_aniso_orth& uaniso,
    std::vector<ftype>& aniso_params, size_t n_spline_params = 20)
{
    const auto& hkls = fcalc.hkl_info();
    auto rfn_p = aniso_scale_fn(fcalc, fobs, uaniso, aniso_params);
    const auto& aniso_rfn = *(rfn_p.resolution_fn);
    for (auto ih = fobs.first(); !ih.last(); ih.next())
    {
        if (!fobs[ih].missing() && !fcalc[ih].missing())
        {
            scaled_fcalc[ih].f() = fcalc[ih].f() * exp(0.5 * aniso_rfn.f(ih));
            scaled_fcalc[ih].phi() = fcalc[ih].phi();
        }
    }
    std::vector<ftype> iso_params(n_spline_params, 1.0);
    BasisFn_spline iso_basisfn(fobs, n_spline_params, 1.0);
    TargetFn_scaleF1I2<F_phi<T>, I_sigI<T>> iso_target(scaled_fcalc, iobs);
    ResolutionFn iso_rfn(hkls, iso_basisfn, iso_target, iso_params);
    for (auto ih=scaled_fcalc.first(); !ih.last(); ih.next())
    {
        if (!scaled_fcalc[ih].missing())
        {
            const ftype fh = iso_rfn.f(ih);
            scaled_fcalc[ih].f() = scaled_fcalc[ih].f() * sqrt(fh > 0.0 ? fh : 0.0);
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
            // Inverse of the Fc -> Fo scale: Fo -> Fc = exp(-0.5*fh).
            auto scale = exp(-0.5 * aniso_rfn.f(ih));
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
