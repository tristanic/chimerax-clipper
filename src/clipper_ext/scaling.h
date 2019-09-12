#pragma once

#include <clipper/clipper.h>

namespace clipper_cx
{

template <typename T>
std::vector<ftype> guess_initial_aniso_gaussian_params(
    const HKL_data<F_sigF<T>>& fobs, const HKL_data<F_phi<T>>& fcalc)
{
    // Fit an isotropic gaussian to the data, and use this to construct the
    // input to the anisotropic version.
    const auto& hkls = fobs.base_hkl_info();
    std::vector<ftype> params = {1.0, 0.01};
    BasisFn_gaussian basisfn;
    TargetFn_scaleF1F2<F_phi<T>, F_sigF<T>> target(fcalc, fobs);
    ResolutionFn_nonlinear rfn(hkls, basisfn, target, params);
    params = rfn.params();
    // std::cerr << "Initial Gaussian params: " << params[0] << ", " << params[1] << std::endl;
    std::vector<ftype> output_params(7,0);
    output_params[0] = params[0];
    output_params[1] = params[1];
    output_params[2] = params[1];
    output_params[3] = params[1];
    return output_params;
}


} // namespace clipper_cx
