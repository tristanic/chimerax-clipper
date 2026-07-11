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

#include <unordered_map>
#include <cstdint>

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/core/atomsf.h>   // for AtomShapeFn::RADIATION (X-ray vs electron)

#include "imex.h"

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{

// Structure factor calculation using bulk solvent taking into account
// van der Waals radii
template<class T>
class SFcalc_obs_bulk_vdw : public SFcalc_obs_base<T>
{
public:
    // default constructor
    SFcalc_obs_bulk_vdw(std::vector<ftype>& params, const T tolerance, const size_t n_threads=8) : params_(&params), nthreads(n_threads), tolerance_frac_(tolerance) {}
    // constructor: shorthand for constructor+operator
    SFcalc_obs_bulk_vdw(HKL_data<F_phi<T> >& fphi, const HKL_data<F_sigF<T> >& fsig,
        const Atom_list& atoms, std::vector<ftype>& params, const T tolerance);
    bool operator() ( HKL_data<datatypes::F_phi<T> >& fphi,
            const HKL_data<datatypes::F_sigF<T> >& fsig, const Atom_list& atoms );
    const T& bulk_frac() const { return bulkfrc; }
    const T& bulk_scale() const { return bulkscl; }
    //! Bulk-solvent structure factors from the most recent operator() call.
    //! f_mask: FFT of the (smoothed) solvent mask.
    //! f_bulk: the additive bulk contribution k_sol·exp(-B_sol)·F_mask, so that
    //! F_total = F_atoms + f_bulk.  Both are empty until operator() has run.
    const HKL_data<F_phi<T>>& f_mask() const { return fmask_; }
    const HKL_data<F_phi<T>>& f_bulk() const { return fbulk_; }
    const size_t& n_threads() const { return nthreads; }
    void set_n_threads(size_t n) { nthreads=n; }
    //! If called, then the bulk solvent scale and B-factor will be re-optimised on the next run.
    void set_bulk_solvent_optimization_needed() { bulk_solvent_optimization_needed_ = true; }
    //! Select the scattering-factor table (X-ray vs electron for micro-ED). Read by
    //! the EDcalc built in each operator() call, so it must be set before the first
    //! (threaded) structure-factor calculation; it is never mutated concurrently.
    void set_radiation(AtomShapeFn::RADIATION r) { radiation_ = r; }
    AtomShapeFn::RADIATION radiation() const { return radiation_; }
    //! When true (default), the solvent mask is occupancy-weighted: a partial-
    //! occupancy atom excludes solvent only fractionally (protein content =
    //! min(1, sum of covering occupancies)).  When false, every non-zero-occupancy
    //! atom stamps a full binary sphere (original Jiang & Brunger behaviour).  Read
    //! by the EDcalc built in each operator() call and folded into mask_signature,
    //! so toggling it invalidates the mask cache on the next call.
    void set_occupancy_weighted(bool b) { occupancy_weighted_ = b; }
    bool occupancy_weighted() const { return occupancy_weighted_; }
    //! All-reflections scaling: when true, the bulk-solvent/scale fit uses ALL
    //! reflections instead of the seeded per-bin subset -- exact and unbiased, at
    //! the cost of a larger (one-off) fit. Default false: the fast, reproducible
    //! seeded subset used for live maps. The bin/seed knobs below are ignored when
    //! this is true. Invalidates the cached subset (the observed subset differs).
    void set_all_reflections(bool b) { all_reflections_ = b; scaling_subset_valid_ = false; }
    bool all_reflections() const { return all_reflections_; }
    //! Size of the reflection subset used for the (non-all-reflections) scale fit:
    //! reflections per resolution bin x number of bins. Defaults 500 x 20. Changing
    //! either invalidates the cached subset.
    void set_scaling_reflections_per_bin(size_t n) { scaling_refls_per_bin_ = n; scaling_subset_valid_ = false; }
    size_t scaling_reflections_per_bin() const { return scaling_refls_per_bin_; }
    void set_scaling_num_bins(size_t n) { scaling_num_bins_ = n; scaling_subset_valid_ = false; }
    size_t scaling_num_bins() const { return scaling_num_bins_; }
    //! Seed for the fixed reflection subset. The subset is drawn ONCE (seeded) and
    //! cached across recalculations, so the scale fit is reproducible; vary the seed
    //! to probe sensitivity to which reflections were chosen. Changing it invalidates
    //! the cached subset. Default 10061865.
    void set_scaling_seed(size_t s) { scaling_seed_ = s; scaling_subset_valid_ = false; }
    size_t scaling_seed() const { return scaling_seed_; }
private:
    AtomShapeFn::RADIATION radiation_ = AtomShapeFn::XRAY;
    bool occupancy_weighted_ = true;
    bool all_reflections_ = false;
    size_t scaling_refls_per_bin_ = 500;
    size_t scaling_num_bins_ = 20;
    size_t scaling_seed_ = 10061865;
    bool bulk_solvent_optimization_needed_ = true;
    //! True once a bulk-solvent solve has run, so subsequent solves can warm-start
    //! from the stored (bulkscl, bulk_u) rather than cold-starting.
    bool bulk_solvent_ever_optimized_ = false;
    std::vector<ftype> *const params_;
    size_t nthreads;
    T bulkfrc, bulkscl;
    T bulk_u;
    T tolerance_frac_;
    HKL_data<F_phi<T>> fmask_, fbulk_;
    //! Cache for the solvent-mask transform.  The mask (and hence fmask_, plus
    //! the mean solvent fraction bulkfrc) is a function of atom positions,
    //! elements, and -- when occupancy_weighted_ is on -- per-atom occupancy
    //! values (otherwise only zero-occupancy status).  It is invariant to
    //! B-factor changes.  When the mask-determining inputs are unchanged from the
    //! previous call (e.g. across B-factor refinement macrocycles, where
    //! coordinates and occupancies are fixed) the cached fmask_ is reused and a
    //! full EDcalc_mask + FFT is skipped -- roughly half the cost of operator().
    //! mask_signature_ fingerprints exactly those inputs; mask_cache_valid_
    //! guards the first call.  Note: with occupancy weighting on, an occupancy
    //! change now invalidates the mask (it genuinely alters the solvent content),
    //! so occupancy-refinement macrocycles rebuild the mask between cycles.
    bool mask_cache_valid_ = false;
    //! 64-bit FNV-1a fingerprint of the mask-determining inputs (atom count, and
    //! per-atom coordinates, element, and occupancy -- as a value when weighting
    //! is on, else a zero-occupancy flag).  A plain integer (not a heap-allocating
    //! std::string) so the member is trivially destructible and safe across the
    //! DLL boundary of this CLIPPER_CX_IMEX-exported class.
    uint64_t mask_signature_ = 0;
    //! Cache for the seeded reflection subset used by the (non-all-reflections)
    //! scale fit.  The selection depends only on the HKL set, the bin knobs, and
    //! scaling_seed_ -- none of which change during a refinement -- so the chosen
    //! HKL_info and the observed-F subset (fsig is invariant across recalcs) are
    //! built once and reused, sparing the HKL_info rebuild + obs-subset copy every
    //! fit.  Only the fcalc-dependent subsets are refreshed per call.  Any knob/seed
    //! setter, and the all_reflections toggle, clears scaling_subset_valid_; a
    //! change to the observed data re-derives via the manager, which drops the cache
    //! through set_all_reflections.
    bool scaling_subset_valid_ = false;
    HKL_info scaling_subset_hkls_;
    HKL_data<F_sigF<T>> scaling_subset_fsig_;
    //! Compute the mask-determining fingerprint for `atoms`.
    uint64_t mask_signature(const Atom_list& atoms) const;
}; // class SFcalc_obs_bulk_vdw


} // namespace clipper_cx
