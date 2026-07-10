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
private:
    AtomShapeFn::RADIATION radiation_ = AtomShapeFn::XRAY;
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
    //! elements and zero-occupancy status ONLY -- it is invariant to B-factor
    //! and occupancy-value changes.  When those inputs are unchanged from the
    //! previous call (e.g. across B-factor/occupancy refinement macrocycles,
    //! where coordinates are fixed) the cached fmask_ is reused and a full
    //! EDcalc_mask + FFT is skipped -- roughly half the cost of operator().
    //! mask_signature_ fingerprints exactly those inputs; mask_cache_valid_
    //! guards the first call.
    bool mask_cache_valid_ = false;
    //! 64-bit FNV-1a fingerprint of the mask-determining inputs (atom count, and
    //! per-atom coordinates, element and zero-occupancy flag).  A plain integer
    //! (not a heap-allocating std::string) so the member is trivially destructible
    //! and safe across the DLL boundary of this CLIPPER_CX_IMEX-exported class.
    uint64_t mask_signature_ = 0;
    //! Compute the mask-determining fingerprint for `atoms`.
    uint64_t mask_signature(const Atom_list& atoms) const;
}; // class SFcalc_obs_bulk_vdw


} // namespace clipper_cx
