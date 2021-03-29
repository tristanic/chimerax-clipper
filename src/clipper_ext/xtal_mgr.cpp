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

#include "xtal_mgr.h"
#include "french_wilson.h"
#include "math_ext.h"
#include "scaling.h"
#include "aniso_scale.h"

#include <set>
#include <memory>
#include <algorithm>

#include <iostream> // debugging
#include <iomanip>

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{

Xmap_details::Xmap_details(const Xmap_details& other)
    : hkl_info_(other.hkl_info()), base_coeffs_(&other.base_coeffs()),
      b_sharp_(other.b_sharp()), is_difference_map_(other.is_difference_map()),
      exclude_missing_(other.exclude_missing()),
      exclude_freer_(other.exclude_free_reflections()),
      fill_(other.fill_with_fcalc())
{
    xmap_ = std::unique_ptr<Xmap<ftype32>>(new Xmap<ftype32>(other.xmap()));
    coeffs_ = std::unique_ptr<HKL_data<F_phi<ftype32>>>( new HKL_data<F_phi<ftype32>>(other.coeffs()) );
}


Xtal_mgr_base::Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
    const Grid_sampling& grid_sampling)
    : hklinfo_(hklinfo), free_flags_(free_flags), grid_sampling_(grid_sampling)
{
    cell_ = hklinfo.cell();
    fcalc_ = HKL_data<F_phi<ftype32>>(hklinfo_);
    base_2fofc_ = HKL_data<F_phi<ftype32>>(hklinfo_);
    base_fofc_ = HKL_data<F_phi<ftype32>>(hklinfo_);
    phi_fom_ = HKL_data<Phi_fom<ftype32>>(hklinfo_);
    usage_ = HKL_data<Flag>(hklinfo_);

} // Xtal_mgr_base


Xtal_mgr_base::Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
    const Grid_sampling& grid_sampling, const HKL_data<F_sigF<ftype32>>& fobs)
    : Xtal_mgr_base(hklinfo, free_flags, grid_sampling)
{
    fobs_original_ = fobs;
    fobs_ = HKL_data<F_sigF<ftype32>>(hklinfo_);
    remove_outliers(fobs_original_, fobs_, 750, OUTLIER_REJECTION_LIMIT);
    set_freeflag(guess_free_flag_value(free_flags_));
} // Xtal_mgr_base

Xtal_mgr_base::Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
    const Grid_sampling& grid_sampling, const HKL_data<F_sigF_ano<ftype32>>& fobs_anom)
    : Xtal_mgr_base(hklinfo, free_flags, grid_sampling)
{
    fobs_original_ = HKL_data<F_sigF<ftype32>>(hklinfo_);
    fobs_ = HKL_data<F_sigF<ftype32>>(hklinfo_);
    for (auto ih = fobs_anom.first(); !ih.last(); ih.next())
    {
        fobs_original_[ih].f() = fobs_anom[ih].f();
        fobs_original_[ih].sigf() = fobs_anom[ih].sigf();
    }
    remove_outliers(fobs_original_, fobs_, 750, OUTLIER_REJECTION_LIMIT);
    set_freeflag(guess_free_flag_value(free_flags_));
}

Xtal_mgr_base::Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
    const Grid_sampling& grid_sampling, const HKL_data<I_sigI<ftype32>>& iobs)
    : Xtal_mgr_base(hklinfo, free_flags, grid_sampling)
{
    fobs_original_ = HKL_data<F_sigF<ftype32>>(hklinfo_);
    fobs_ = HKL_data<F_sigF<ftype32>>(hklinfo_);
    french_wilson::french_wilson(iobs, fobs_original_);
    remove_outliers(fobs_original_, fobs_, 750, OUTLIER_REJECTION_LIMIT);
    set_freeflag(guess_free_flag_value(free_flags_));
}

Xtal_mgr_base::Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
    const Grid_sampling& grid_sampling, const HKL_data<I_sigI_ano<ftype32>>& iobs_ano)
    : Xtal_mgr_base(hklinfo, free_flags, grid_sampling)
{
    fobs_original_ = HKL_data<F_sigF<ftype32>>(hklinfo_);
    fobs_ = HKL_data<F_sigF<ftype32>>(hklinfo_);
    french_wilson::french_wilson(iobs_ano, fobs_original_);
    remove_outliers(fobs_original_, fobs_, 750, OUTLIER_REJECTION_LIMIT);
    set_freeflag(guess_free_flag_value(free_flags_));
}




int
Xtal_mgr_base::guess_free_flag_value(const HKL_data<Flag>& flags)
{
    HKL_info::HKL_reference_index ih;
    int f_min = 1e6, f_max=-1e6;
    std::set<int> flag_vals;
    for (ih=flags.first(); !ih.last(); ih.next())
    {
        const auto &f = flags[ih];
        if (!f.missing())
        {
            if (f.flag() < f_min) f_min = f.flag();
            else if (f.flag() > f_max) f_max = f.flag();
            flag_vals.insert(f.flag());
        }
    }
    if (flag_vals.size() == 1)
    {
        std::ostringstream err_str;
        err_str << "ERROR: all flags in the given free flags array have a " <<
            "single value (" << *(flag_vals.begin()) << "). Please retry " <<
            "with a reflection file containing a valid free set.";
        throw std::runtime_error(err_str.str());
    }
    if (f_max > 1) /* assume CCP4 */ return 0;
    if (f_min == -1 && flag_vals.size()==2) /* assume SHELX */ return -1;
    else
    {
        // Count the different values, and choose the one that makes most sense
        std::unordered_map<int, int> val_counts;
        for (const auto& v: flag_vals)
            val_counts[v] = 0;
        for (ih=flags.first(); !ih.last(); ih.next())
            val_counts[flags[ih].flag()] += 1;

        std::unordered_map<int, float> val_fracs;
        // Convert to fractions of total reflections
        for (const auto v: flag_vals)
        {
            std::cout << "Candidate R-free flag value: " << v << " Count: " << val_counts[v] << std::endl;
            val_fracs[v] = ((float)val_counts[v]) / flags.hkl_info().num_reflections();
        }

        // Weed out those that don't fit the criteria
        std::set<int> candidates;
        int last_candidate=-10;
        for (auto v: flag_vals)
        {
            if ((val_fracs[v] > 0.005 && val_fracs[v] < 0.15)
            || (val_fracs[v] < 0.15 && val_counts[v] > 2000))
            {
                candidates.insert(v);
                last_candidate = v;
            }
        }
        if (candidates.size() != 1) {
            throw std::runtime_error("Cannot determine the correct value for the "
            "free reflections in the given FreeR_flag array. Please convert your "
            "reflection file to CCP4 or PHENIX format using sftools or "
            "phenix.reflection_file_editor. ");
            return -1;
        }
        return last_candidate;
    }
} // guess_free_flag_value

void
Xtal_mgr_base::set_freeflag(int f)
{
    freeflag_ = f;
    HKL_info::HKL_reference_index ih;
    for (ih = usage_.first(); !ih.last(); ih.next())
    {
        auto f = free_flags_[ih];
        auto fobs = fobs_[ih];
        if (!f.missing() && !fobs.missing() && !(f.flag()==freeflag_))
            usage_[ih].flag() = SFweight_spline<ftype32>::BOTH;
        else
            usage_[ih].flag() = SFweight_spline<ftype32>::NONE;
    }
} // set_freeflag

void
Xtal_mgr_base::generate_fcalc(const Atom_list& atoms)
{
    // auto start_time = std::chrono::system_clock::now();
    bulk_solvent_calculator_(fcalc_, fobs_, atoms);
    // auto end_time = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed = end_time-start_time;
    // std::cout << "Complete Fcalc generation took " << elapsed.count() << " seconds." << std::endl;
    // guess_initial_gaussian_params_();
    fcalc_initialized_ = true;

    calculate_r_factors();
} // generate_fcalc

// Generate the standard set of map coefficients
void
Xtal_mgr_base::generate_base_map_coeffs()
{
    if (!fcalc_initialized())
        throw std::runtime_error("No Fcalc values have been calculated! Run "
            " generate_fcalc() on a suitable set of atoms first!");
    map_calculator_(base_2fofc_, base_fofc_, phi_fom_, scaled_fobs(), fcalc_, usage_);

    coeffs_initialized_=true;
} // generate_base_map_coeffs

void
Xtal_mgr_base::remove_outliers(const HKL_data<F_sigF<ftype32>>& f_sigf_in,
    HKL_data<F_sigF<ftype32>>& f_sigf_out, int reflections_per_bin,
    ftype high_p_cutoff, ftype beamstop_cutoff, ftype beamstop_d_min)
{
    ftype high_e_cutoff_centric = approx_inverse_erfc(high_p_cutoff) * sqrt(2.0);
    ftype high_e_cutoff_acentric = sqrt(log(1/high_p_cutoff));
    ftype beamstop_invrsq_cutoff = 1/pow(beamstop_d_min, 2.0);
    ftype beamstop_e_cutoff_centric = pow(beamstop_cutoff, 2.0)*M_PI/2.0;
    ftype beamstop_e_cutoff_acentric = -log(1.0-beamstop_cutoff);
    HKL_data<E_sigE<ftype32>> esige(hklinfo_);
    esige.compute(f_sigf_in, data32::Compute_EsigE_from_FsigF());
    int n_bins = std::max(hklinfo_.num_reflections() / reflections_per_bin, 1);
    //BasisFn_gaussian iso_basisfn;
    BasisFn_binner basisfn(hklinfo_, n_bins, 1.0);
    TargetFn_scaleEsq<E_sigE<ftype32>> targetfn(esige);
    std::vector<ftype> params = {(ftype)n_bins, 1.0};
    ResolutionFn rfn(hklinfo_, basisfn, targetfn, params);
    ftype32 e_norm;
    size_t outlier_count = 0;
    for (auto ih=esige.first(); !ih.last(); ih.next())
    {
        const auto& e = esige[ih];
        if (f_sigf_in[ih].missing())
        {
            f_sigf_out[ih].set_null();
            continue;
        }
        auto scale = sqrt(rfn.f(ih));
        e_norm = e.E() * scale;
        // e.sigE() *= scale;
        bool centric = ih.hkl_class().centric();
        if ((centric && e_norm >= high_e_cutoff_centric)
                || (!centric && e_norm >= high_e_cutoff_acentric))
        {
            f_sigf_out[ih].set_null();
            outlier_count++;
            continue;
        } else if (ih.invresolsq() < beamstop_invrsq_cutoff)
        {
            if ((centric && e_norm <= beamstop_e_cutoff_centric)
                || (!centric && e_norm <= beamstop_e_cutoff_acentric))
            {
                std::cout << "Removed probable beamstop shadow artefact at "
                    << ih.hkl().format() << std::endl;
                f_sigf_out[ih].set_null();
                outlier_count++;
                continue;
            }
        }
        f_sigf_out[ih] = f_sigf_in[ih];
    }
    if (outlier_count > 0)
    {
        std::cout << "Removed " << outlier_count << " outliers from reflection data." << std::endl;
    }

}


HKL_data<F_phi<ftype32>>
Xtal_mgr_base::scaled_fcalc()
{
    HKL_data<F_phi<ftype32>> sfcalc(hklinfo_);
    scale_fcalc_to_fobs(fcalc_, fobs_, sfcalc, aniso_ucryst_, aniso_scale_params_);
    // std::cout << "Overall Anisou scaling fcalc to fobs: " << aniso_ucryst_.format() << std::endl;
    return sfcalc;
    // int nparams = 20;
    // auto basisfn = choose_basisfn_(fcalc_, fobs_, nparams);
    // return scaled_fcalc_(fcalc_, fobs_, basisfn);

} // scaled_fcalc

HKL_data<F_phi<ftype32>>
Xtal_mgr_base::scaled_fcalc_(const HKL_data<F_phi<ftype32>>& fcalc,
    const HKL_data<F_sigF<ftype32>>& fobs, std::unique_ptr<BasisFn_base, BasisFn_Deleter>& basisfn)
{
    std::vector<ftype> params(basisfn->num_params());
    if (scaling_method_==ANISO_GAUSSIAN)
        params=aniso_scale_params_;
    HKL_data<F_phi<ftype32>> ret(hklinfo_);
    //TargetFn_scaleF1F2<F_phi<ftype32>, F_sigF<ftype32>> targetfn (fcalc_, fobs_);
    TargetFn_scaleF1F2<F_sigF<ftype32>, F_phi<ftype32>> targetfn (fobs_, fcalc_);
    ResolutionFn_nonlinear rfn(hklinfo_, *basisfn, targetfn, params, NONLINEAR_DAMP);
    params = rfn.params();
    if (scaling_method_==ANISO_GAUSSIAN)
        if (static_cast<BasisFn_aniso_gaussian*>(basisfn.get())->u_aniso_orth(params).det() <= 0.0)
            std::cerr << "WARNING: Overall anisotropic B-factor is non-positive definite when scaling fcalc to fobs!" << std::endl;
    // if (scaling_method_==ANISO_GAUSSIAN)
    //     std::cerr << "Anisotropic gaussian params: " << params[0] << ", " << params[1]
    //         << ", " << params[2] << ", " << params[3] << ", " << params[4]
    //         << ", " << params[5] << ", " << params[6] << std::endl;
    HKL_info::HKL_reference_index ih;
    for (ih=fcalc_.first(); !ih.last(); ih.next())
    {
        if (!fobs_[ih].missing())
        {
            ret[ih].f() = fcalc_[ih].f()/sqrt(rfn.f(ih));
            ret[ih].phi() = fcalc_[ih].phi();
        }
    }
    return ret;
}

HKL_data<F_sigF<ftype32>>
Xtal_mgr_base::scaled_fobs()
{
    HKL_data<F_sigF<ftype32>> sfobs(hklinfo_);
    aniso_scale_fobs_to_fcalc(fcalc_, fobs_, sfobs, aniso_ucryst_, aniso_scale_params_);
    // std::cout << "Overall Anisou scaling fcalc to fobs: " << aniso_ucryst_.format() << std::endl;
    return sfobs;


    // int nparams = 20;
    // auto basisfn = choose_basisfn_(fcalc_, fobs_, nparams);
    // return scaled_fobs_(fcalc_, fobs_, basisfn);

} // scaled_fcalc

HKL_data<F_sigF<ftype32>>
Xtal_mgr_base::scaled_fobs_(const HKL_data<F_phi<ftype32>>& fcalc,
    const HKL_data<F_sigF<ftype32>>& fobs, std::unique_ptr<BasisFn_base, BasisFn_Deleter>& basisfn)
{
    std::vector<ftype> params(basisfn->num_params());
    if (scaling_method_==ANISO_GAUSSIAN)
        params=aniso_scale_params_;
    HKL_data<F_sigF<ftype32>> ret(hklinfo_);
    TargetFn_scaleF1F2<F_sigF<ftype32>, F_phi<ftype32>> targetfn (fobs_, fcalc_);
    ResolutionFn_nonlinear rfn(hklinfo_, *basisfn, targetfn, params, NONLINEAR_DAMP);
    params = rfn.params();
    if (scaling_method_==ANISO_GAUSSIAN)
    if (static_cast<BasisFn_aniso_gaussian*>(basisfn.get())->u_aniso_orth(params).det() <= 0.0)
            std::cerr << "WARNING: Overall anisotropic B-factor is non-positive definite when scaling fobs to fcalc!" << std::endl;
    // if (scaling_method_==ANISO_GAUSSIAN)
    //     std::cerr << "Anisotropic gaussian params: " << params[0] << ", " << params[1]
    //         << ", " << params[2] << ", " << params[3] << ", " << params[4]
    //         << ", " << params[5] << ", " << params[6] << std::endl;
    HKL_info::HKL_reference_index ih;
    for (ih=fobs_.first(); !ih.last(); ih.next())
    {
        if (!fobs_[ih].missing())
        {
            ret[ih].f() = sqrt(rfn.f(ih))*fobs_[ih].f();
            ret[ih].sigf() = sqrt(rfn.f(ih))*fobs_[ih].sigf();
        }
    }
    return ret;
}


std::unique_ptr<BasisFn_base, BasisFn_Deleter> Xtal_mgr_base::choose_basisfn_(
    const HKL_data<F_phi<ftype32>>& fcalc, HKL_data<F_sigF<ftype32>> fobs,
    int nparams)
{
    if (scaling_method_ == NOT_CHOSEN)
    {
        std::cout << "Automatically detecting best scaling method..." << std::endl;
        std::vector<int> candidates = {SPLINE, ANISO_GAUSSIAN};
        std::vector<ftype> rfactors;
        for (auto c: candidates) {
            scaling_method_ = c;
            std::cout << scaling_method_names_.at(c);
            calculate_r_factors();
            std::cout << ": Rwork: " << std::setprecision(3) << rwork_
                      << " Rfree: " << std::setprecision(3) << rfree_ << std::endl;
            rfactors.push_back(rwork_);
        }
        auto best_index = std::distance(rfactors.begin(), std::min_element(rfactors.begin(), rfactors.end()));
        std::cout << "Chose " << scaling_method_names_.at(candidates[best_index])
                  << " for this dataset." << std::endl;
        scaling_method_ = candidates[best_index];
    }
    auto deleter = BasisFn_Deleter(scaling_method_);
    //std::unique_ptr<BasisFn_base> ret;
    switch (scaling_method_)
    {
        case SPLINE:         return std::unique_ptr<BasisFn_base, BasisFn_Deleter>(new BasisFn_spline(fobs, nparams, 1.0), deleter);
        case ANISO_GAUSSIAN: return std::unique_ptr<BasisFn_base, BasisFn_Deleter>(new BasisFn_aniso_gaussian(), deleter);
        default:
        {
            std::stringstream err;
            err << "Attempted to make a basis function with unknown type " << scaling_method_ << "!";
            throw std::runtime_error(err.str());
        }
    }
    return std::unique_ptr<BasisFn_base, BasisFn_Deleter>(nullptr, deleter);
}


std::vector<float>
Xtal_mgr_base::scaling_function()
{
    auto l = hklinfo_.num_reflections();
    std::vector<float> ret(l);

    int nparams = 20;
    std::vector<ftype> params(nparams, 0.0);
    auto basisfn = choose_basisfn_(fcalc_, fobs_, nparams);
    //BasisFn_spline basisfn(fobs_, nparams, 1.0);
    // BasisFn_aniso_gaussian basisfn;
    TargetFn_scaleF1F2<F_phi<ftype32>, F_sigF<ftype32>> targetfn (fcalc_, fobs_);
    ResolutionFn_nonlinear rfn(hklinfo_, *basisfn, targetfn, params, NONLINEAR_DAMP);
    HKL_info::HKL_reference_index ih;
    size_t i=0;
    for (ih = fcalc_.first(); !ih.last(); ih.next(), ++i)
        ret[i] = rfn.f(ih);
    return ret;
}


void
Xtal_mgr_base::calculate_r_factors()
{
    if (!fcalc_initialized())
        throw std::runtime_error("No Fcalc values have been calculated! Run "
            " generate_fcalc() on a suitable set of atoms first!");

    // int nparams = 7;
    // std::vector<ftype> params(nparams, 0.0);
    // //BasisFn_spline basisfn(hklinfo_, nparams, 1.0);
    // BasisFn_aniso_gaussian basisfn;
    // TargetFn_scaleF1F2<F_phi<ftype32>, F_sigF<ftype32>> targetfn (fcalc_, fobs_);
    // ResolutionFn_nonlinear rfn(hklinfo_, basisfn, targetfn, params);

    auto s_fcalc = scaled_fcalc();

    HKL_info::HKL_reference_index ih;
    // for standard rwork, rfree
    ftype sum_fwork=0, sum_ffree=0, sum_dwork=0, sum_dfree=0;
    // for sigma-weighted rwork, rfree
    ftype sum_wfwork2=0, sum_wffree2=0, sum_wdwork2=0, sum_wdfree2=0;
    for (ih=fcalc_.first(); !ih.last(); ih.next())
    {
        const auto& fo = fobs_[ih];
        const auto& fc = s_fcalc[ih];
        const auto& fflag = free_flags_[ih];
        if (!fo.missing() && !fc.missing())
        {
            ftype eps = ih.hkl_class().epsilon();
            ftype two_on_eps = 2.0/eps;
            // Shouldn't really have to do this, but the spline function can
            // dip below zero for sudden changes
            // float scale = rfn.f(ih);
            // if (scale < 0)
            // {
            //   std::cerr << "Below-zero scaling value of " << scale << " at " << ih.hkl().format() << "!" << std::endl;
            //   scale = 0;
            // }
            // auto scaled_fcalc = sqrt(scale)*fc.f();
            if (fflag.flag()==freeflag_) {
                sum_ffree+=two_on_eps*fo.f();
                sum_dfree+=two_on_eps*std::abs(fo.f()-fc.f());

                sum_wffree2 += 1/fo.sigf()*pow(two_on_eps*fo.f(), 2);
                sum_wdfree2 += 1/fo.sigf()*pow(two_on_eps*(fo.f()-fc.f()), 2);
            } else {
                sum_fwork+=two_on_eps*fo.f();
                sum_dwork+=two_on_eps*std::abs(fo.f()-fc.f());
                if (Util::is_nan(sum_fwork) || Util::is_nan(sum_dwork))
                {
                    std::cerr << "Invalid value encountered at " << ih.hkl().format() << "!" << std::endl;
                    std::cerr << "eps: " << eps << " two_on_eps: " << two_on_eps << std::endl;
                    std::cerr << "fobs: " << fo.f() << std::endl;
                    std::cerr << "fcalc: " << fcalc_[ih].f() << std::endl;
                    // std::cerr << "ResolutionFn value: " << rfn.f(ih) << std::endl;
                    std::cerr << "scaled_fcalc: " << fc.f() << std::endl;
                    std::cerr << "fsigf: " << fo.f() << ", " << fo.sigf() << std::endl;
                    throw std::runtime_error("Invalid value encountered in R-factor calculation!");
                }

                sum_wfwork2 += 1/fo.sigf()*pow(two_on_eps*fo.f(), 2);
                sum_wdwork2 += 1/fo.sigf()*pow(two_on_eps*(fo.f()-fc.f()), 2);
            }
        }
    }
    rfree_ = sum_dfree/sum_ffree;
    rwork_ = sum_dwork/sum_fwork;

    w_rfree_ = sqrt(sum_wdfree2/sum_wffree2);
    w_rwork_ = sqrt(sum_wdwork2/sum_wfwork2);
} // calculate_r_factors


void
Xtal_mgr_base::apply_b_factor_sharpening(HKL_data<F_phi<ftype32>>& coeffs,
    const ftype& bsharp)
{
    coeffs.compute(
        coeffs,
        Compute_scale_u_iso<F_phi<ftype32>>(
            1.0, Util::b2u(bsharp)));
} // apply_b_factor_sharpening

void
Xtal_mgr_base::add_xmap(const std::string& name,
    const ftype& bsharp, bool is_difference_map,
    bool exclude_missing_reflections,
    bool exclude_free_reflections, bool fill_with_fcalc, size_t num_threads)
{
    if (!coeffs_initialized())
        throw std::logic_error("You need to calculate base coefficients before "
        "creating any maps! Run init() on a suitable set of atoms first.");
    // Throw an error if a map with that name already exists
    auto it = maps_.find(name);
    if (it != maps_.end())
        throw std::logic_error("Each map must have a unique name!");
    if (is_difference_map)
        maps_.emplace(name, Xmap_details(hklinfo_, &base_fofc_, bsharp, grid_sampling_, is_difference_map,
                              exclude_missing_reflections, exclude_free_reflections, fill_with_fcalc));
    else
        maps_.emplace(name, Xmap_details(hklinfo_, &base_2fofc_, bsharp, grid_sampling_, is_difference_map,
                            exclude_missing_reflections, exclude_free_reflections, fill_with_fcalc));
    // Create a copy for thread-safe calculations

    recalculate_map(name, num_threads);
} // add_xmap

void
Xtal_mgr_base::recalculate_map(const std::string& name, size_t num_threads)
{
    try {
        recalculate_map(maps_.at(name), num_threads);
    } catch (std::out_of_range& exc) {
        std::stringstream msg;
        msg << "Could not find map with name " << name << "!";
        throw std::out_of_range(msg.str());
    }
} // recalculate_map

void
Xtal_mgr_base::recalculate_map(Xmap_details& xmd, size_t num_threads)
{
    if (xmd.exclude_free_reflections()) {
        if (xmd.fill_with_fcalc() && !xmd.is_difference_map())
            set_map_free_terms_to_dfc(xmd.base_coeffs(), xmd.coeffs());
        else
            set_map_free_terms_to_zero(xmd.base_coeffs(), xmd.coeffs());
    } else {
        xmd.coeffs() = xmd.base_coeffs();
    }
    if (xmd.exclude_missing())
        remove_missing_reflections_from_map_coeffs(xmd.coeffs(), fobs_);
    if (xmd.b_sharp() != 0.0)
        apply_b_factor_sharpening(xmd.coeffs(), xmd.b_sharp());
    xmd.xmap().fft_from(xmd.coeffs(), num_threads);
    // xmd.map_stats() = Map_stats(xmd.xmap());
    xmd.update_map_stats();
} // recalculate_map

void
Xtal_mgr_base::recalculate_all(const Atom_list& atoms)
{
    generate_fcalc(atoms);
    generate_base_map_coeffs();
    for (auto& it: maps_)
        recalculate_map(it.second);
} //recalculate_map


void
Xtal_mgr_base::set_map_free_terms_to_zero(const HKL_data<F_phi<ftype32>>& source,
    HKL_data<F_phi<ftype32>>& dest)
{
    if (!coeffs_initialized())
        throw std::runtime_error("Coefficients have not yet been calculated!");
    HKL_info::HKL_reference_index ih;
    for (ih = source.first(); !ih.last(); ih.next())
    {
        if (usage_[ih].missing() || usage_[ih].flag() == SFweight_spline<ftype32>::NONE)
        {
            if (!source[ih].missing())
                dest[ih].f() = 0;
            else
                dest[ih].set_null();
        } else {
            dest[ih] = source[ih];
        }
    }
} // set_map_free_terms_to_zero

void
Xtal_mgr_base::set_map_free_terms_to_dfc(const HKL_data<F_phi<ftype32>>& source,
    HKL_data<F_phi<ftype32>>& dest)
{
    if (!coeffs_initialized())
        throw std::runtime_error("Coefficients have not yet been calculated!");
    HKL_info::HKL_reference_index ih;
    HKL_data<Flag_bool> flag(source.hkl_info());
    for (ih=flag.first(); !ih.last(); ih.next())
        flag[ih].flag() = (!source[ih].missing() && (usage_[ih].flag()!=SFweight_base<ftype32>::NONE));
    const auto& param_s = map_calculator_.params_scale();
    // const auto& param_w = map_calculator_.params_error();
    BasisFn_spline basisfn( flag, param_s.size(), 1.0);
    for (ih=source.first(); !ih.last(); ih.next())
    {
        const auto& fpo = source[ih];
        if (usage_[ih].missing() || (usage_[ih].flag() == SFweight_spline<ftype32>::NONE))
        //if (flag[ih].flag())
        {
            if (!fpo.missing()) {
                auto s = basisfn.f_s( ih.invresolsq(), param_s);
                dest[ih].f() = s*fcalc_[ih].f();
            } else {
                dest[ih].set_null();
            }
        } else {
            dest[ih] = fpo;
        }
    }
} // set_map_free_terms_to_dfc

void
Xtal_mgr_base::remove_missing_reflections_from_map_coeffs(HKL_data<F_phi<ftype32>>& coeffs,
    const HKL_data<F_sigF<ftype32>>& f_sigf)
{
    if (!coeffs_initialized())
        throw std::runtime_error("Coefficients have not yet been calculated!");
    HKL_info::HKL_reference_index ih;
    for (ih=coeffs.first(); !ih.last(); ih.next())
    {
        if (f_sigf[ih].missing())
        {
            // std::cerr << "Setting " << ih.hkl().format() << " to null" << std::endl;
            coeffs[ih].set_null();
        }
    }


}



// THREADED IMPLEMENTATIONS
Xtal_thread_mgr::Xtal_thread_mgr(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
    const Grid_sampling& grid_sampling, const HKL_data<F_sigF<ftype32>>& fobs,
    const size_t num_threads)
    : num_threads_(num_threads)
    {
        mgr_ = std::unique_ptr<Xtal_mgr_base>(new Xtal_mgr_base(hklinfo, free_flags, grid_sampling, fobs));
        mgr_->bulk_solvent_calculator_.set_n_threads(num_threads);
    }

Xtal_thread_mgr::Xtal_thread_mgr(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
    const Grid_sampling& grid_sampling, const HKL_data<F_sigF_ano<ftype32>>& fobs_anom,
    const size_t num_threads)
    : num_threads_(num_threads)
    {
        mgr_ = std::unique_ptr<Xtal_mgr_base>(new Xtal_mgr_base(hklinfo, free_flags, grid_sampling, fobs_anom));
        mgr_->bulk_solvent_calculator_.set_n_threads(num_threads);
    }

Xtal_thread_mgr::Xtal_thread_mgr(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
    const Grid_sampling& grid_sampling, const HKL_data<I_sigI<ftype32>>& iobs,
    const size_t num_threads)
    : num_threads_(num_threads)
    {
        mgr_ = std::unique_ptr<Xtal_mgr_base>(new Xtal_mgr_base(hklinfo, free_flags, grid_sampling, iobs));
        mgr_->bulk_solvent_calculator_.set_n_threads(num_threads);
    }

Xtal_thread_mgr::Xtal_thread_mgr(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
    const Grid_sampling& grid_sampling, const HKL_data<I_sigI_ano<ftype32>>& iobs_anom,
    const size_t num_threads)
    : num_threads_(num_threads)
    {
        mgr_ = std::unique_ptr<Xtal_mgr_base>(new Xtal_mgr_base(hklinfo, free_flags, grid_sampling, iobs_anom));
        mgr_->bulk_solvent_calculator_.set_n_threads(num_threads);
    }


void
// Xtal_thread_mgr::recalculate_all(const Atom_list& atoms)
Xtal_thread_mgr::recalculate_all(std::vector<uintptr_t> cxatoms)
{
    deletion_guard();
    if (thread_running())
        throw std::runtime_error("Map recalculation already in progress! Run "
        "apply_new_maps() first.");
    auto cxa = std::vector<atomstruct::Atom*>();
    for (const auto& ptr: cxatoms)
        cxa.push_back(reinterpret_cast<atomstruct::Atom*>(ptr));
    //auto cxa = static_cast<atomstruct::Atom**>(cxatoms);
    atoms_ = bridge::clipper_atoms_from_cx_atoms_threaded(cxa.data(), cxa.size(), num_threads(), ignore_hydrogens());
    master_thread_result_ = std::async(std::launch::async,
        &Xtal_thread_mgr::recalculate_all_, this, atoms_);
} // recalculate_all

bool
Xtal_thread_mgr::recalculate_all_(const Atom_list& atoms)
{
    ready_ = false;
    try
    {
        mgr_->generate_fcalc(atoms);
        mgr_->generate_base_map_coeffs();
    } catch (...)
    {
        ready_ = true;
        std::cerr << "failed at generate_fcalc" << std::endl;
        std::rethrow_exception(std::current_exception());
    }

    auto map_names = mgr_->map_names();
    auto nmaps = map_names.size();
    size_t maps_per_thread = (size_t) ceil(( (float)(mgr_->n_maps())) /num_threads_);
    size_t threads_per_map = std::max(num_threads_ / nmaps * maps_per_thread, size_t(1));
    // Make copies of all maps to work on. Expensive, but necessary for thread
    // safety.
    // xmap_thread_results_.clear();
    // for (const auto& name: map_names)
    //     xmap_thread_results_.emplace(name, mgr_.maps_[name]);
    // auto start_time = std::chrono::system_clock::now();
    size_t n=0;
    std::vector<std::future<bool>> results;
    while (n<nmaps)
    {
        //std::cerr << "Recalculating maps " << n << " to " << std::min(n+maps_per_thread, nmaps)-1 << std::endl;
        results.push_back(std::async(std::launch::async,
            &Xtal_thread_mgr::recalculate_inner_, this, map_names,
            n, std::min(n+maps_per_thread, nmaps), threads_per_map));
        n+=maps_per_thread;
    }
    // Wait for all threads to finish
    std::vector<std::exception_ptr> exceptions;
    for (auto& r: results)
    {
        try {
            r.get();
        } catch (...)
        {
            std::cerr << "Caught exception in recalculate_all_" << std::endl;
            exceptions.push_back(std::current_exception());
        }
    }
    // auto end_time = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed = end_time-start_time;
    // std::cout << "FFTs took " << elapsed.count() << " seconds." << std::endl;

    ready_ = true;
    if (exceptions.size())
        std::rethrow_exception(exceptions[0]);
    return true;
} // recalculate_all_

bool
Xtal_thread_mgr::recalculate_inner_(const std::vector<std::string>& names,
    size_t i_min, size_t i_max, size_t threads)
{
    for (size_t i= i_min; i<i_max; ++i)
    {
        // std::cerr << "Recalculating map " << i << std::endl;
        mgr_->recalculate_map(xmap_thread_results_[names[i]], threads);
    }
    return true;
} // recalculate_inner_

void
Xtal_thread_mgr::apply_new_maps()
{
    deletion_guard();
    if (!thread_running())
        throw std::logic_error("You haven't generated any new maps yet! Call "
            "recalculate_all() first.");
    master_thread_result_.get();
    try {
        for (auto& it: xmap_thread_results_)
        {
            auto& source = it.second;
            auto& target = mgr_->maps_.at(it.first);
            source.coeffs_.swap(target.coeffs_);
            source.xmap_.swap(target.xmap_);
            source.map_stats_.swap(target.map_stats_);
        }
        // for (const auto& it: xmap_thread_results_)
        //     mgr_.maps_.at(it.first).xmap() = it.second.xmap();
        // xmap_thread_results_.clear();
        ready_ = false;
    } catch (std::out_of_range& err){
        std::stringstream msg;
        msg << "Something went wrong when applying new maps!";
        throw std::out_of_range(msg.str());
    }
} // apply_new_maps

void
Xtal_thread_mgr::set_freeflag(int f)
{
    deletion_guard();
    finalize_threads_if_necessary();
    mgr_->set_freeflag(f);
}

HKL_data<F_phi<ftype32>>
Xtal_thread_mgr::fcalc()
{
    deletion_guard();
    finalize_threads_if_necessary();
    return mgr_->fcalc();
}

HKL_data<F_phi<ftype32>>
Xtal_thread_mgr::scaled_fcalc()
{
    deletion_guard();
    finalize_threads_if_necessary();
    return mgr_->scaled_fcalc();
}

HKL_data<F_phi<ftype32>>
Xtal_thread_mgr::base_fofc()
{
    deletion_guard();
    finalize_threads_if_necessary();
    return mgr_->base_fofc();
}

HKL_data<F_phi<ftype32>>
Xtal_thread_mgr::base_2fofc()
{
    deletion_guard();
    finalize_threads_if_necessary();
    return mgr_->base_2fofc();
}

HKL_data<Phi_fom<ftype32>>
Xtal_thread_mgr::weights()
{
    deletion_guard();
    finalize_threads_if_necessary();
    return mgr_->weights();
}

void
Xtal_thread_mgr::init(std::vector<uintptr_t> cxatoms)
{
    deletion_guard();
    auto cxa = std::vector<atomstruct::Atom*>();
    for (const auto& ptr: cxatoms)
        cxa.push_back(reinterpret_cast<atomstruct::Atom*>(ptr));
    // auto cxa = static_cast<atomstruct::Atom**>(cxatoms);
    atoms_ = bridge::clipper_atoms_from_cx_atoms_threaded(cxa.data(), cxa.size(), num_threads(), ignore_hydrogens());
    init_(atoms_);
}

void
Xtal_thread_mgr::init_(const Atom_list& atoms)
{
    finalize_threads_if_necessary();
    mgr_->init(atoms);
}

void
Xtal_thread_mgr::add_xmap(const std::string& name,
    const ftype& bsharp, bool is_difference_map,
    bool exclude_missing_reflections,
    bool exclude_free_reflections, bool fill_with_fcalc)
{
    deletion_guard();
    finalize_threads_if_necessary();
    mgr_->add_xmap(name, bsharp, is_difference_map,
        exclude_missing_reflections, exclude_free_reflections, fill_with_fcalc, num_threads_);
    xmap_thread_results_.emplace(name, mgr_->maps_[name]);

}

void
Xtal_thread_mgr::delete_xmap(const std::string& name)
{
    deletion_guard();
    finalize_threads_if_necessary();
    mgr_->delete_xmap(name);
    xmap_thread_results_.erase(name);
}

void
Xtal_thread_mgr::deletion_guard() const
{
    if (mgr_==nullptr)
    {
        throw std::runtime_error("Xtal_thread_mgr C++ side has been deleted!");
    }
}

void
Xtal_thread_mgr::delete_all()
{
    deletion_guard();
    finalize_threads_if_necessary();
    atoms_.clear();
    xmap_thread_results_.clear();
    mgr_.reset();
}

} //namespace clipper_cx;
