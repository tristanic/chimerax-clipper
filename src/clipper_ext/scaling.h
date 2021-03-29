#pragma once

#include <clipper/clipper.h>

namespace clipper_cx
{


template <typename T>
class TargetFn_scaleFcalcFobs: public TargetFn_base
{
public:
    //! constructor: takes the datalist against which to calc target
    TargetFn_scaleFcalcFobs( const HKL_data<F_phi<T>>& fcalc, const HKL_data<F_sigF<T>>& fobs );
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return QUADRATIC; }
  private:
    const HKL_data<F_phi<T>>* fcalc_;
    const HKL_data<F_sigF<T>>* fobs_;

};

template <typename T>
class TargetFn_scaleFobsFcalc: public TargetFn_base
{
public:
    //! constructor: takes the datalist against which to calc target
    TargetFn_scaleFobsFcalc( const HKL_data<F_sigF<T>>& fobs, const HKL_data<F_phi<T>>& fcalc );
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return QUADRATIC; }
  private:
    const HKL_data<F_sigF<T>>* fobs_;
    const HKL_data<F_phi<T>>* fcalc_;

};


template <class T>
T quick_r(const HKL_data<F_phi<T>>& fcalc, const HKL_data<F_sigF<T>>& fobs,
    const ResolutionFn& rfn)
{
    T r = 0.0;
    for ( auto ih = fobs.first(); !ih.last(); ih.next() )
    {
        if ( !fobs[ih].missing() ) {
            double eps = ih.hkl_class().epsilon();
            auto f1 = pow(fcalc[ih].f(), 2.0)/eps / rfn.f(ih);
            auto f2 = pow(fobs[ih].f(), 2.0)/eps;
            if (std::isnan(f1) || std::isnan(f2))
            {
                if (std::isnan(rfn.f(ih)))
                {
                    std::stringstream err_msg;
                    err_msg << "Scaling function value is NaN! Parameters are:" << std::endl;
                    for (const auto& p: rfn.params())
                        err_msg << p << ", ";
                    throw std::runtime_error(err_msg.str());
                }
                std::cerr << "NaN encountered at " << ih.hkl().format() << "!"
                  << "Fobs: " << f2 << ", fcalc: " << fcalc[ih].f() << ", scale: " << rfn.f(ih) << std::endl;
            }
            r += fabs(f1-f2);
        }
    }
    return r;

}

template <typename T>
std::vector<ftype> guess_initial_aniso_gaussian_params(
    const HKL_data<F_sigF<T>>& fobs, const HKL_data<F_phi<T>>& fcalc, T& r)
{
    // Fit an isotropic gaussian to the data, and use this to construct the
    // input to the anisotropic version.
    const auto& hkls = fobs.base_hkl_info();
    T sum_fobs_lo = 0, sum_fcalc_lo=0, sum_fobs_hi = 0, sum_fcalc_hi = 0;
    auto invressq_lim = pow(1/(hkls.resolution().limit()+0.5), 2);
    for (auto ih = fobs.first(); !ih.last(); ih.next())
    {
        if (!fobs[ih].missing() && !fcalc[ih].missing())
        {
          if ( ih.invresolsq() < 0.01 )
          {
            sum_fobs_lo += fobs[ih].f();
            sum_fcalc_lo += fcalc[ih].f();
          } else if (ih.invresolsq() > invressq_lim)
          {
            sum_fobs_hi += fobs[ih].f();
            sum_fcalc_hi += fcalc[ih].f();
          }
        }
    }

    auto p0 = log(pow(sum_fcalc_lo/sum_fobs_lo, 2));
    auto p1 = (p0-log(pow(sum_fcalc_hi/sum_fobs_hi, 2)))/((invressq_lim+hkls.resolution().invresolsq_limit())/2);

    std::vector<ftype> params = {p0, p1};

    std::cout << "Initial params: " << p0 << ", " << p1 << std::endl; // DELETEME

    BasisFn_gaussian basisfn;
    //TargetFn_scaleF1F2<F_phi<T>, F_sigF<T>> target(fcalc, fobs);
    TargetFn_scaleFobsFcalc<T> target(fobs, fcalc);
    ResolutionFn_nonlinear rfn(hkls, basisfn, target, params, 0.0);
    params = rfn.params();

    std::cout << "Fitted params: " << params[0] << ", " << params[1] << std::endl; // DELETEME
    // std::cerr << "Initial Gaussian params: " << params[0] << ", " << params[1] << std::endl;
    std::vector<ftype> output_params(7,0);
    output_params[0] = params[0];
    output_params[1] = params[1];
    output_params[2] = params[1];
    output_params[3] = params[1];
    BasisFn_aniso_gaussian aniso_basisfn;
    ResolutionFn_nonlinear aniso_rfn(hkls, aniso_basisfn, target, output_params, 0.0);
    try {
      r = quick_r(fcalc, fobs, aniso_rfn);
      output_params = aniso_rfn.params();
    } catch (std::runtime_error) {}
    std::cout << "Fitted aniso params: ";
    for (auto p: output_params)
      std::cout << p << ", ";
    std::cout << std::endl;

    return output_params;
}


template<typename T>
TargetFn_scaleFcalcFobs<T>::TargetFn_scaleFcalcFobs( const HKL_data<F_phi<T>>& fcalc, const HKL_data<F_sigF<T>>& fobs )
    : fcalc_(&fcalc), fobs_(&fobs)
{}

template<typename T>
TargetFn_base::Rderiv TargetFn_scaleFcalcFobs<T>::rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const
{
  Rderiv result;
  const F_phi<T>& fc = (*fcalc_)[ih];
  const F_sigF<T>& fo = (*fobs_)[ih];
  if ( !fc.missing() && !fo.missing() ) {
    const ftype eps = ih.hkl_class().epsilon();
    const ftype f1 = pow( fc.f(), 2 ) / eps;
    const ftype f2 = pow( fo.f(), 2 ) / eps;
    const ftype sigf = fo.sigf() * fo.f() / eps;
    const ftype d = (fh*f1 - f2)/sigf;
    result.r = d * d / f1;
    result.dr = 2.0 * d;
    result.dr2 = 2.0 * f1 / sigf;
  } else {
    result.r = result.dr = result.dr2 = 0.0;
  }
  return result;
}

template<typename T>
TargetFn_scaleFobsFcalc<T>::TargetFn_scaleFobsFcalc( const HKL_data<F_sigF<T>>& fobs, const HKL_data<F_phi<T>>& fcalc )
    : fobs_(&fobs), fcalc_(&fcalc)
{}


template<typename T>
TargetFn_base::Rderiv TargetFn_scaleFobsFcalc<T>::rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const
{
  Rderiv result;
  result.r = result.dr = result.dr2 = 0.0;
  const F_phi<T>& fc = (*fcalc_)[ih];
  const F_sigF<T>& fo = (*fobs_)[ih];
  if ( !fc.missing() && !fo.missing()) {
    const ftype eps = ih.hkl_class().epsilon();
    const ftype f1 = pow( fo.f(), 2 ) / eps;
    const ftype f2 = pow( fc.f(), 2 ) / eps;
    ftype w;
    if (fo.sigf() > 0)
      w = std::min(pow(fo.f()/fo.sigf(), 2), 5.0);
    else
      w = 1.0;
    // ftype sigf1 = 1.0; //fo.sigf() * fo.f() / eps + 1;
    const ftype d = fh*f1 - f2;
    result.r = w * d * d / f1;
    result.dr =  w * 2.0 * d;
    result.dr2 = 2.0 * w * f1;
  } else {
    result.r = result.dr = result.dr2 = 0.0;
  }
  return result;
}



} // namespace clipper_cx
