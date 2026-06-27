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

// (guess_initial_aniso_gaussian_params removed.)  The anisotropic scaling now
// fits in log space via BasisFn_log_aniso_gaussian + TargetFn_scaleLogF1F2 (a
// linear least-squares solve, see aniso_scale.h and the bulk-solvent optimiser),
// which is stable from zeros and needs no isotropic-then-anisotropic seeding.


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
    const ftype d = fh*f1 - f2;
    result.r = d * d / f1;
    result.dr =  2.0 * d;
    result.dr2 = 2.0 * f1;
  } else {
    result.r = result.dr = result.dr2 = 0.0;
  }
  return result;
}



} // namespace clipper_cx
