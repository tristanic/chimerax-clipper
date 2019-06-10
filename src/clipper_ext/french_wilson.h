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

#include "french_wilson.cpp"

using namespace clipper;

namespace clipper_cx { namespace french_wilson {

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
