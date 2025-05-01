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

#include <clipper/clipper.h>

using namespace clipper;

namespace clipper_cx { namespace french_wilson {

static const double acentric_zf[71] = {
    0.423, 0.428, 0.432, 0.437, 0.442, 0.447, 0.453, 0.458, 0.464,  0.469,
    0.475, 0.482, 0.488, 0.495, 0.502, 0.509, 0.516, 0.524, 0.532,  0.540,
    0.549, 0.557, 0.567, 0.576, 0.586, 0.597, 0.608, 0.619, 0.631,  0.643,
    0.656, 0.670, 0.684, 0.699, 0.714, 0.730, 0.747, 0.765, 0.783,  0.802,
    0.822, 0.843, 0.865, 0.887, 0.911, 0.935, 0.960, 0.987, 1.014,  1.042,
    1.070, 1.100, 1.130, 1.161, 1.192, 1.224, 1.257, 1.289, 1.322,  1.355,
    1.388, 1.421, 1.454, 1.487, 1.519, 1.551, 1.583, 1.615, 1.646,  1.676,
    1.706
};

static const double acentric_sig_zf[71] = {
    0.216, 0.218, 0.220, 0.222, 0.224, 0.226, 0.229, 0.231, 0.234,  0.236,
    0.239, 0.241, 0.244, 0.247, 0.250, 0.253, 0.256, 0.259, 0.262,  0.266,
    0.269, 0.272, 0.276, 0.279, 0.283, 0.287, 0.291, 0.295, 0.298,  0.302,
    0.307, 0.311, 0.315, 0.319, 0.324, 0.328, 0.332, 0.337, 0.341,  0.345,
    0.349, 0.353, 0.357, 0.360, 0.364, 0.367, 0.369, 0.372, 0.374,  0.375,
    0.376, 0.377, 0.377, 0.377, 0.376, 0.374, 0.372, 0.369, 0.366,  0.362,
    0.358, 0.353, 0.348, 0.343, 0.338, 0.332, 0.327, 0.321, 0.315,  0.310,
    0.304
};

static const double centric_zf[81] = {
    0.269, 0.272, 0.276, 0.279, 0.282, 0.286, 0.289, 0.293, 0.297,  0.301,
    0.305, 0.309, 0.314, 0.318, 0.323, 0.328, 0.333, 0.339, 0.344,  0.350,
    0.356, 0.363, 0.370, 0.377, 0.384, 0.392, 0.400, 0.409, 0.418,  0.427,
    0.438, 0.448, 0.460, 0.471, 0.484, 0.498, 0.512, 0.527, 0.543,  0.560,
    0.578, 0.597, 0.618, 0.639, 0.662, 0.687, 0.713, 0.740, 0.769,  0.800,
    0.832, 0.866, 0.901, 0.938, 0.976, 1.016, 1.057, 1.098, 1.140,  1.183,
    1.227, 1.270, 1.313, 1.356, 1.398, 1.439, 1.480, 1.519, 1.558,  1.595,
    1.632, 1.667, 1.701, 1.735, 1.767, 1.799, 1.829, 1.859, 1.889,  1.917,
    1.945
};

static const double centric_sig_zf[81] = {
    0.203, 0.205, 0.207, 0.209, 0.211, 0.214, 0.216, 0.219, 0.222,  0.224,
    0.227, 0.230, 0.233, 0.236, 0.239, 0.243, 0.246, 0.250, 0.253,  0.257,
    0.261, 0.265, 0.269, 0.273, 0.278, 0.283, 0.288, 0.293, 0.298,  0.303,
    0.309, 0.314, 0.320, 0.327, 0.333, 0.340, 0.346, 0.353, 0.361,  0.368,
    0.375, 0.383, 0.390, 0.398, 0.405, 0.413, 0.420, 0.427, 0.433,  0.440,
    0.445, 0.450, 0.454, 0.457, 0.459, 0.460, 0.460, 0.458, 0.455,  0.451,
    0.445, 0.438, 0.431, 0.422, 0.412, 0.402, 0.392, 0.381, 0.370,  0.360,
    0.349, 0.339, 0.330, 0.321, 0.312, 0.304, 0.297, 0.290, 0.284,  0.278,
    0.272
};

    //! simple mean |I| target
    /*! This class implements the target function for calculating mean
      |I| as a function of position in reciprocal space. It
      includes the appropriate multiplicity correction, and so can be
      applied to any type with an 'i' member with the same dimensions as
      an |I|.
      */
template<class T>
class TargetFn_meanI : public TargetFn_base
{
public:
    //! constructor: takes the datalist against which to calc target, and power
    TargetFn_meanI( const HKL_data<T>& hkl_data) : hkl_data_(&hkl_data) {}
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return QUADRATIC; }
  private:
    const HKL_data<T>* hkl_data_;
};


const int MIN_REFLECTIONS_PER_BIN = 40;

template<typename T, template<typename TT> class Itype>
void french_wilson_acentric(const Itype<T>& isigi, F_sigF<T>& fsigf,
    const ftype& imean, const T rejection_cutoff=-4.0);

template<typename T, template<typename TT> class Itype>
void french_wilson_centric(const Itype<T>& isigi, F_sigF<T>& fsigf,
    const ftype& imean, const T rejection_cutoff=-4.0);

template<typename T, template<typename TT> class Itype>
void french_wilson(const HKL_data<Itype<T>>& isigi, HKL_data<F_sigF<T>>& fsigf,
    const T rejection_cutoff=-4.0, const int max_bins=60);





// IMPLEMENTATIONS

template<typename T>
T linear_interpolate(const double& lb, const double& ub, const T& frac)
{
    return frac*(ub-lb) + lb;
}

template<typename T> void sanity_check(const ftype& imean, const T& sigi)
{
    // negative sigmas are not possible
    if (sigi <= 0)
    {
        std::ostringstream msg;
        msg << "French and Wilson processing encountered a negative or zero sigI "
            << "value. Since negative sigmas are impossible, calculation has "
            << "been aborted. Please check your data.";
        throw std::runtime_error(msg.str());
    }
    if (imean <= 0)
    {
        std::ostringstream msg;
        msg << "Negative or zero mean intensity encountered in resolution bin "
            << "during French and Wilson processing. This indicates a severe "
            << "issue with your data. Processing has been cancelled.";
        throw std::runtime_error(msg.str());
    }
}

template<class T>
TargetFn_base::Rderiv TargetFn_meanI<T>::rderiv(const HKL_info::HKL_reference_index& ih, const ftype& fh) const
{
    Rderiv result;
    const HKL_data<T>& data = *hkl_data_;
    if (!data[ih].missing())
    {
        ftype d = fh - data[ih].I(); //(ftype(data[ih].I()) / sqrt(ih.hkl_class().epsilon()));
        result.r = d*d;
        result.dr = 2.0*d;
        result.dr2 = 2.0;
    } else {
        result.r = result.dr = result.dr2 = 0.0;
    }
    return result;
}

template<typename T, template<typename TT> class Itype>
void french_wilson_acentric(const Itype<T>& isigi, F_sigF<T>& fsigf,
    const ftype& imean, const T rejection_cutoff)
{
    const auto& i = isigi.I();
    const auto& sigi = isigi.sigI();
    sanity_check(imean, sigi);
    T i_over_sigi = i/sigi;
    T h = i_over_sigi - (sigi/imean);
    if (i_over_sigi < rejection_cutoff+0.3 || h < rejection_cutoff)
    {
        fsigf.set_null();
        return;
    }
    T f, sigf;

    if (h<3.0)
    {
        T root_sigi = sqrt(sigi);
        T point = 10.0*(h+4.0);
        int pt_1 = int(point);
        int pt_2 = pt_1+1;
        T delta = point-pt_1;
        f = linear_interpolate(acentric_zf[pt_1], acentric_zf[pt_2], delta) * root_sigi;
        sigf = linear_interpolate(acentric_sig_zf[pt_1], acentric_sig_zf[pt_2], delta) * root_sigi;
    } else {
        f = sqrt(h*sigi);
        sigf = 0.5*sigi/f;
    }
    fsigf.f() = f;
    fsigf.sigf() = sigf;
}


template<typename T, template<typename TT> class Itype>
void french_wilson_centric(const Itype<T>& isigi, F_sigF<T>& fsigf,
    const ftype& imean, const T rejection_cutoff)
{
    const auto& i = isigi.I();
    const auto& sigi = isigi.sigI();
    sanity_check(imean, sigi);
    T i_over_sigi = i/sigi;
    T h = i_over_sigi - (sigi/(2.0*imean));
    if (i_over_sigi < rejection_cutoff+0.3 || h < rejection_cutoff)
    {
        fsigf.set_null();
        return;
    }
    T f, sigf;
    if (h<4.0)
    {
        T root_sigi = sqrt(sigi);
        T point = 10.0*(h+4.0);
        int pt_1 = int(point);
        int pt_2 = pt_1 + 1;
        T delta = point-pt_1;
        f = linear_interpolate(centric_zf[pt_1], centric_zf[pt_2], delta) * root_sigi;
        sigf = linear_interpolate(centric_sig_zf[pt_1], centric_sig_zf[pt_2], delta) * root_sigi;
    } else {
        // as per CCTBX: French&Wilson + added x^6 term in Taylor series
        auto h_2 = 1.0/(h*h);
        auto h_4 = pow(h_2, 2.0);
        auto h_6 = h_2*h_4;
        // posterior of f
        auto post_f = sqrt(h)* (1.0 - (3.0/8.0)*h_2 - (87.0/128.0)*h_4 - (2889.0/1024.0)*h_6 );
        // posterior of sigf
        auto post_sigf = sqrt( h *( (1.0/4.0)*h_2 + (15.0/32.0)*h_4 + (273.0/128.0)*h_6) );

        f = post_f*sqrt(sigi);
        sigf = post_sigf * sqrt(sigi);
    }
    fsigf.f() = f;
    fsigf.sigf() = sigf;
}

template<typename T, template<typename TT> class Itype>
void french_wilson(const HKL_data<Itype<T>>& isigi, HKL_data<F_sigF<T>>& fsigf,
    const T rejection_cutoff, const int max_bins)
{
    const HKL_info& hkls = isigi.base_hkl_info();

    auto nrefl = hkls.num_reflections();
    // spline requires at least 3 bins
    if (nrefl < MIN_REFLECTIONS_PER_BIN*3)
    {
        std::ostringstream msg;
        msg << "Intensity dataset contains only " << nrefl
            << "reflections, which is insufficient to provide reliable "
            << "statistics for French & Wilson processing. Are you sure this "
            << "is a protein dataset?";
        throw std::runtime_error(msg.str());
    }
    int reflections_per_bin = std::max(nrefl/max_bins, MIN_REFLECTIONS_PER_BIN);
    int nbins = nrefl/reflections_per_bin;


    std::vector<ftype> params (nbins, 1.0);
    BasisFn_spline basisfn( hkls, nbins, 1.0 );
    TargetFn_meanI< Itype<T> > targetfn(isigi);
    ResolutionFn rfn (hkls, basisfn, targetfn, params);

    for (auto ih = isigi.first(); !ih.last(); ih.next() )
    {
        if (ih.hkl_class().centric() )
            french_wilson_centric( isigi[ih], fsigf[ih], rfn.f(ih) );
        else
            french_wilson_acentric( isigi[ih], fsigf[ih], rfn.f(ih) );
    }

}



}}
