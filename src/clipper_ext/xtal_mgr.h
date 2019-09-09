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

#pragma once
//#include <string>
#include <unordered_map>
#include <future>
#include <memory>
#include <chrono>


#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include "imex.h"
#include "sfcalc_obs_vdw.h"
#include "chimerax_bridge.h"

namespace clipper_cx { // ChimeraX extensions to clipper

using namespace clipper;
using namespace clipper::datatypes;

// Forward declaration
class Xtal_thread_mgr;

// Container for a crystallographic map and its coefficients, including the
// details needed to regenerate it. Just stores everything - calculations are
// handled by Xtal_mgr_base
struct CLIPPER_CX_IMEX Xmap_details
{
    friend class Xtal_thread_mgr;

public:
    Xmap_details() {} // null constructor
    Xmap_details(const Xmap_details& other);
    inline Xmap_details(const HKL_info& hklinfo, HKL_data<F_phi<ftype32>> const* base_coeffs, const ftype& b_sharp,
        const Grid_sampling& grid_sampling, bool is_difference_map=false,
        bool exclude_missing_reflections=false, bool exclude_free_reflections=true,
        bool fill_with_fcalc=true)
        : hkl_info_(hklinfo),
          base_coeffs_(base_coeffs), b_sharp_(b_sharp),
          is_difference_map_(is_difference_map),
          exclude_missing_(exclude_missing_reflections),
          exclude_freer_(exclude_free_reflections),
          fill_(fill_with_fcalc)
    {
        // initialise empty data objects
        coeffs_ = std::unique_ptr<HKL_data<F_phi<ftype32>>> (
            new HKL_data<F_phi<ftype32>>(hklinfo));
        xmap_ = std::unique_ptr<Xmap<ftype32>> (
            new Xmap<ftype32>(hklinfo.spacegroup(), hklinfo.cell(), grid_sampling));
    }
    Xmap_details operator= (const Xmap_details& other) { return Xmap_details(other); }

    inline const HKL_info& hkl_info() const { return hkl_info_; }
    inline const Xmap<ftype32>& xmap() const { return *xmap_; }
    inline Xmap<ftype32>& xmap() { return *xmap_; }
    inline const Map_stats& map_stats() const { return *map_stats_; }
    // inline Map_stats& map_stats() { return map_stats_; }
    void update_map_stats() { map_stats_ = std::unique_ptr<Map_stats> (new Map_stats(*xmap_)); }
    inline const HKL_data<F_phi<ftype32>>& coeffs() const { return *coeffs_; }
    inline HKL_data<F_phi<ftype32>>& coeffs() { return *coeffs_; }
    inline const HKL_data<F_phi<ftype32>>& base_coeffs() const { return *base_coeffs_; }
    inline const ftype& b_sharp() const { return b_sharp_; }
    inline bool is_difference_map() const { return is_difference_map_; }
    inline bool exclude_missing() const { return exclude_missing_; }
    inline bool exclude_free_reflections() const { return exclude_freer_; }
    inline bool fill_with_fcalc() const { return fill_; }


private:
    HKL_info hkl_info_;
    // Base coefficients before removing free reflections, sharpening etc..
    // These are common to all maps, so only a pointer is kept here. The actual
    // object is held by Xtal_mgr_base
    HKL_data<F_phi<ftype32>> const* base_coeffs_;
    // Final coefficients used to generate xmap_
    ftype b_sharp_=0;
    // If this is a difference map, fill_ is ignored
    bool is_difference_map_=false;
    // Exclude missing reflections from map? CAUTION: if false these will be
    // replaced by f_calc, which can lead to severe model bias for data with
    // a large proportion of missing reflections.
    bool exclude_missing_=false;
    // Exclude free reflections from map?
    bool exclude_freer_=true;
    // if exclude_freer_, replace free reflections with DFcalc?
    bool fill_=true;

    std::unique_ptr<HKL_data<F_phi<ftype32>>> coeffs_;
    // Sharpening/smoothing B-factor to apply
    // The real-space map

    std::unique_ptr<Xmap<ftype32>> xmap_;
    // Xmap<ftype32> xmap_;

    // max, min, mean, std_dev, range
    std::unique_ptr<Map_stats> map_stats_;
    // Map_stats map_stats_;

}; // class Xmap_details



class BasisFn_Deleter;


//! General manager class for handling operations on a set of crystallographic
//  maps and their data
class CLIPPER_CX_IMEX Xtal_mgr_base
{
    friend class Xtal_thread_mgr;
public:
    enum resol_fn {
        NOT_CHOSEN = 0,
        SPLINE,
        ANISO_GAUSSIAN,
    };

    //Xtal_mgr_base() {} // default constructor
    Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling, const HKL_data<F_sigF<ftype32>>& fobs);

    Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling, const HKL_data<F_sigF_ano<ftype32>>& fobs_anom);


    Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling, const HKL_data<I_sigI<ftype32>>& iobs);

    Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling, const HKL_data<I_sigI_ano<ftype32>>& iobs_anom);

    inline bool fcalc_initialized() const { return fcalc_initialized_; }
    inline bool coeffs_initialized() const { return coeffs_initialized_; }
    inline int freeflag() const { return freeflag_; }
    void set_freeflag(int f);
    inline const HKL_data<Flag>& flags() const { return free_flags_; }

    inline const ftype& rwork() { return rwork_; }
    inline const ftype& rfree() { return rfree_; }
    inline const ftype& weighted_rwork() { return w_rwork_; }
    inline const ftype& weighted_rfree() { return w_rfree_; }

    inline bool ignore_hydrogens() const { return ignore_hydrogens_; }
    inline void set_ignore_hydrogens(bool flag) { ignore_hydrogens_ = flag; }

    const Xmap_details& map_details(const std::string& name) const { return maps_.at(name); }

    size_t n_maps() const { return maps_.size(); }
    std::vector<std::string> map_names() const {
        std::vector<std::string> names;
        for (const auto& it: maps_)
            names.push_back(it.first);
        return names;
    }

    inline const ftype& bulk_frac()
    {
        if (!fcalc_initialized())
            throw std::runtime_error("Coefficients have not been calculated yet!");
        return bulk_solvent_calculator_.bulk_frac();
    }
    inline const ftype& bulk_scale()
    {
        if (!fcalc_initialized())
            throw std::runtime_error("Coefficients have not been calculated yet!");
        return bulk_solvent_calculator_.bulk_scale();
    }

    inline const HKL_data<F_sigF<ftype32>>& fobs() const
    {
        return fobs_;
    }
    inline const HKL_data<F_phi<ftype32>>& fcalc() const
    {
        if (!fcalc_initialized())
            throw std::runtime_error("Coefficients have not been calculated yet!");
        return fcalc_;
    }

    // Fcalc scaled to Fobs
    HKL_data<F_phi<ftype32>> scaled_fcalc();

    inline const HKL_data<Flag>& usage() const { return usage_; }
    inline const HKL_data<F_phi<ftype32>>& base_2fofc() const
    {
        if (!coeffs_initialized())
            throw std::runtime_error("Coefficients have not yet been calculated!");
        return base_2fofc_;
    }
    inline const HKL_data<F_phi<ftype32>>& base_fofc() const
    {
        if (!coeffs_initialized())
            throw std::runtime_error("Coefficients have not yet been calculated!");
        return base_fofc_;
    }

    inline const HKL_data<Phi_fom<ftype32>>& weights() const
    {
        if (!coeffs_initialized())
            throw std::runtime_error("Weights have not yet been calculated!");
        return phi_fom_;
    }

    int guess_free_flag_value(const HKL_data<Flag>& flags);

    // Regenerate Fcalc from a set of ChimeraX Atoms
    void generate_fcalc(const Atom_list& atoms);

    // Generate the standard set of map coefficients
    void generate_base_map_coeffs();

    // One-stop wrapper for to generate fcalcs and base map coefficients
    inline void init(const Atom_list& atoms)
    {
        generate_fcalc(atoms);
        generate_base_map_coeffs();
    }

    // Calculate Rwork and Rfree. Called automatically by generate_fcalc()
    void calculate_r_factors();

    // Apply in-place B-factor sharpening to a set of map coefficients
    void apply_b_factor_sharpening(HKL_data<F_phi<ftype32>>& coeffs, const ftype& bsharp);

    void add_xmap(const std::string& name,
        const ftype& bsharp, bool is_difference_map=false,
        bool exclude_missing_reflections=false,
        bool exclude_free_reflections=true, bool fill_with_fcalc=true);

    void delete_xmap(const std::string& name) { maps_.erase(name); }

    // Recalculate map coefficients and real-space map for one previously-stored
    // set of parameters (will throw std::out_of_range if nothing has been stored
    // under the given name)
    void recalculate_map(const std::string& name, size_t num_threads=1);

    // Recalculate a map in-place
    void recalculate_map(Xmap_details& xmd, size_t num_threads=1);

    // Generate fresh Fcalc, and regenerate all existing maps
    void recalculate_all(const Atom_list& atoms);

    inline const std::unordered_map<std::string, Xmap_details>& maps() const { return maps_; }

    //
    inline const Xmap<ftype32>& get_xmap(const std::string& name) const { return maps_.at(name).xmap(); }
    inline const Map_stats& get_map_stats(const std::string& name) { return maps_.at(name).map_stats(); }

    // Return the values of the current scaling function scaling fcalc to fobs
    // for each (H,K,L)
    std::vector<float> scaling_function();
    ~Xtal_mgr_base() {
#ifdef CLIPPER_VERBOSE
      std::cerr << "Shutting down Xtal_mgr_base" << std::endl;
#endif
    }


protected:
    // Safety
    // Have we generated any fcalcs yet?
    bool fcalc_initialized_=false;
    // Have we generated map coefficients yet?
    bool coeffs_initialized_=false;

private:
    const std::map<int, std::string> scaling_method_names_ = {
        {SPLINE, "Isotropic spline"},
        {ANISO_GAUSSIAN, "Anisotropic gaussian"}
    };
    //! Common portion of all constructors
    Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling);

    int scaling_method_ = NOT_CHOSEN;
    std::unique_ptr<BasisFn_base, BasisFn_Deleter> choose_basisfn_(
        const HKL_data<F_phi<ftype32>>& fcalc, HKL_data<F_sigF<ftype32>> fobs,
        int nparams);

    HKL_data<F_phi<ftype32>> scaled_fcalc_(
        const HKL_data<F_phi<ftype32>>& fcalc,
        const HKL_data<F_sigF<ftype32>>& fobs,
        std::unique_ptr<BasisFn_base, BasisFn_Deleter>& basisfn);


    // constants
    //static constexpr ftype ONE_1_ON_4_PI_SQUARED = 1/(4*M_PI*M_PI);
    int freeflag_ = 0;

    bool ignore_hydrogens_ = true;

    // Basic information
    HKL_info hklinfo_;
    Cell cell_;
    // Original free-r flags from the reflection data file
    HKL_data<Flag> free_flags_;
    // Flags defining which reflections to use in map calculations
    HKL_data<Flag> usage_;
    Grid_sampling grid_sampling_;
    // Observed structure factors
    HKL_data<F_sigF<ftype32>> fobs_;
    // Real-space maps

    // Base map coefficients before filling, sharpening etc.
    HKL_data<F_phi<ftype32>> base_2fofc_;
    HKL_data<F_phi<ftype32>> base_fofc_;

    // Weighting coefficients
    HKL_data<Phi_fom<ftype32>> phi_fom_;

    std::unordered_map<std::string, Xmap_details> maps_;
    // Most recent Fcalc values
    HKL_data<F_phi<ftype32>> fcalc_;

    SFcalc_obs_bulk_vdw<ftype32> bulk_solvent_calculator_ = SFcalc_obs_bulk_vdw<ftype32>();

    // SFcalc_obs_bulk<ftype32> bulk_solvent_calculator_ = SFcalc_obs_bulk<ftype32>();

    // Calculates sigmaa-weighted maximum-likelihood map coefficients
    SFweight_spline<ftype32> map_calculator_ = SFweight_spline<ftype32>();

    // "Straight" or standard rwork/rfree
    ftype rwork_ = 1;
    ftype rfree_ = 1;

    // sigma-weighted rwork/rfree
    ftype w_rwork_ = 1;
    ftype w_rfree_ = 1;

    void set_map_free_terms_to_zero(const HKL_data<F_phi<ftype32>>& source,
        HKL_data<F_phi<ftype32>>& dest);
    void set_map_free_terms_to_dfc(const HKL_data<F_phi<ftype32>>& source,
        HKL_data<F_phi<ftype32>>& dest);
    void remove_missing_reflections_from_map_coeffs(HKL_data<F_phi<ftype32>>& coeffs,
        const HKL_data<F_sigF<ftype32>>& f_sigf);


}; // class Xmap_mgr


//! Parallel handler for map calculations
/*! Uses std::async to push the task of (re)generating maps to a separate thread
    from the main program. Further splits up map calculations into individual
    threads if available.    for (const auto& name: map_names)
    {

    }

*/
class CLIPPER_CX_IMEX Xtal_thread_mgr
{
public:
    Xtal_thread_mgr(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling, const HKL_data<F_sigF<ftype32>>& fobs,
        const size_t num_threads = 1);

    Xtal_thread_mgr(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling, const HKL_data<F_sigF_ano<ftype32>>& fobs_ano,
        const size_t num_threads = 1);

    Xtal_thread_mgr(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling, const HKL_data<I_sigI<ftype32>>& iobs,
        const size_t num_threads = 1);

    Xtal_thread_mgr(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling, const HKL_data<I_sigI_ano<ftype32>>& iobs_anom,
        const size_t num_threads = 1);

    ~Xtal_thread_mgr() {
#ifdef CLIPPER_VERBOSE
      std::cerr << "Shutting down Xtal_thread_mgr" << std::endl;
#endif
    }


    inline size_t num_threads() const { return num_threads_; }
    inline void set_num_threads(size_t n)
    {
        deletion_guard();
        num_threads_=std::max(n, size_t(1));
        mgr_->bulk_solvent_calculator_.set_n_threads(num_threads_);
    }

    bool thread_running() const { return master_thread_result_.valid(); }
    bool ready() const { return ready_; }

    // Recalculate all maps in a separate thread (which will in turn spawn
    // worker threads up to the num_threads() limit). New maps will be stored
    // in xmap_thread_results_ until pushed to the manager by apply_new_maps();

    void recalculate_all(std::vector<uintptr_t> cxatoms);
    //void recalculate_all(const Atom_list& atoms);
    // Swap newly-created maps into the manager seen by Python. This will return
    // an error if recalculate_all() wasn't called first, or block until the
    // threads are finished. For best performance, wait until ready() returns
    // true before calling.
    void apply_new_maps();

    // It would be a Bad Thing(TM) to change parameters in the base manager
    // while a map recalculation is in progress. This function ensures that
    // any recalculations are done before making changes.
    void finalize_threads_if_necessary()
    {
        if (thread_running())
            apply_new_maps();
    }

    // Pass through calls to base manager functions in a thread-safe manner
    inline int freeflag() const { deletion_guard(); return mgr_->freeflag(); }
    void set_freeflag(int f);
    inline const HKL_data<Flag>& flags() const { deletion_guard(); return mgr_->flags(); }
    inline const HKL_data<Flag>& usage_flags() const { deletion_guard(); return mgr_->usage(); }
    inline const ftype& rwork() { deletion_guard(); return mgr_->rwork(); }
    inline const ftype& rfree() { deletion_guard(); return mgr_->rfree(); }
    inline const ftype& weighted_rwork() { deletion_guard(); return mgr_->weighted_rwork(); }
    inline const ftype& weighted_rfree() { deletion_guard(); return mgr_->weighted_rfree(); }

    inline const Xmap_details& map_details(const std::string& name) const { deletion_guard(); return mgr_->map_details(name); }
    inline size_t n_maps() const { deletion_guard(); return mgr_->n_maps(); }
    std::vector<std::string> map_names() const { deletion_guard(); return mgr_->map_names(); }

    inline const ftype& bulk_frac() { deletion_guard(); return mgr_->bulk_frac(); }
    inline const ftype& bulk_scale() { deletion_guard(); return mgr_->bulk_scale(); }

    inline const HKL_data<F_sigF<ftype32>>& fobs() const { deletion_guard(); return mgr_->fobs(); }

    inline bool ignore_hydrogens() const { return mgr_->ignore_hydrogens(); }
    inline void set_ignore_hydrogens(bool flag) { mgr_->set_ignore_hydrogens(flag); }

    // Finalise thread and return a copy
    HKL_data<F_phi<ftype32>> fcalc();

    HKL_data<F_phi<ftype32>> scaled_fcalc();

    // Finalise thread and return a copy
    HKL_data<F_phi<ftype32>> base_fofc();

    HKL_data<F_phi<ftype32>> base_2fofc();

    // Finalise thread and return a copy
    HKL_data<Phi_fom<ftype32>> weights();

    void init(std::vector<uintptr_t> cxatoms);

    void add_xmap(const std::string& name,
        const ftype& bsharp, bool is_difference_map=false,
        bool exclude_missing_reflections=false,
        bool exclude_free_reflections=true, bool fill_with_fcalc=true);

    void delete_xmap(const std::string& name);

    // Delete all "large memory" items (atoms, maps and structure factors).
    // Only really needed because Python's garbage collection doesn't seem to
    // reliably delete the Xtal_thread_mgr C++ object when the Python object
    // goes out of scope.)
    void delete_all();


    inline const Xmap<ftype32>& get_xmap(const std::string& name) { deletion_guard(); return mgr_->get_xmap(name); }

    inline const Map_stats& get_map_stats(const std::string& name) { deletion_guard(); return mgr_->get_map_stats(name); }

    // Return the values of the current scaling function scaling fcalc to fobs
    // for each (H,K,L)
    inline std::vector<float> scaling_function() const { deletion_guard(); return mgr_->scaling_function(); }


private:
    Atom_list atoms_; // Temporary Clipper atoms created from ChimeraX atoms
                      // at each iteration
    //Xtal_mgr_base mgr_;
    std::unique_ptr<Xtal_mgr_base> mgr_;
    size_t num_threads_;
    // std::future does not (yet) have an is_ready() function, so just have the
    // thread set a flag instead.
    bool ready_ = false;
    std::future<bool> master_thread_result_;
    std::unordered_map<std::string, Xmap_details> xmap_thread_results_;
    //std::vector<std::pair<std::string, Xmap_details>> xmap_thread_results_;

    // Safety check to prevent use (and nullptr dereference) after delete() has
    // been called.
    void deletion_guard() const;

    // Master thread function called by recalculate_all();
    bool recalculate_all_(const Atom_list& atoms);
    // Inner threads called by recalculate_all_();
    bool recalculate_inner_(const std::vector<std::string>& names, size_t i_min, size_t i_max, size_t threads=1);

    void init_(const Atom_list& atoms);

}; // Xtal_thread_mgr

class BasisFn_Deleter
{
public:
    BasisFn_Deleter(int fn_type)
        : fn_type_(fn_type) {}

    void operator()(BasisFn_base* ptr)
    {
        switch (fn_type_)
        {
            case Xtal_mgr_base::SPLINE: delete (BasisFn_spline*)ptr; break;
            case Xtal_mgr_base::ANISO_GAUSSIAN: delete (BasisFn_aniso_gaussian*)ptr; break;
            default:
                std::stringstream err;
                err << "Attempted to delete a BasisFn with unknown type id " << fn_type_ << "!";
                throw std::runtime_error(err.str());
        }
    }

private:
    int fn_type_ = 0;
};

} // namespace clipper_cx
