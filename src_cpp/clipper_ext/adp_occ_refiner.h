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
#include <vector>
#include <future>
#include <memory>

#include <clipper/clipper.h>
#include "imex.h"
#include "adp_occ_types.h"
#include "chimerax_bridge.h"   // for AtomAltlocIndex

namespace clipper_cx {

using namespace clipper;
using namespace clipper::datatypes;

// Forward declaration — BFactorOccRefiner is fully defined in adp_occ_refiner.cpp
// to keep Eigen/LBFGSpp out of this header.
class BFactorOccRefiner;

//! Thread manager for background B-factor and occupancy refinement.
/*!
  Uses the Agarwal (1978) real-space gradient formulation.  The caller
  supplies a pre-computed "driving density" d(x) and an Atom_list; all data
  is deep-copied at launch time so no shared mutable state is accessed after
  launch() returns.

  For crystallographic ML refinement d(x) must be the Fo-DFc difference map
  (coefficients (m|Fo|-D|Fc|)*exp(i*phi_calc)), which is the negative of the
  ML gradient density G(h) = dT/dFcalc*(h).  The target
  T = -integral(d*rho_calc) then gives dT/dp = -sum_x d*d(rho)/dp
    = +sum_x G*d(rho)/dp, the correct Agarwal ML gradient.
  Do NOT pass the 2Fo-DFc map: its m*rho_obs component drives B-factors to
  the lower bound.

  Follows the same std::async + ready() polling pattern as Xtal_thread_mgr.

  Typical usage:
  @code
    BFactorOccRefinerThread mgr(cfg);
    mgr.launch(atoms, xtal_mgr.f_obs, xtal_mgr.weights, xtal_mgr.usage_flags);
    // ... poll mgr.ready() in the ChimeraX event loop ...
    if (mgr.ready()) {
        const auto& u_iso = mgr.refined_u_iso();  // Å², same order as atoms
        // write back to ChimeraX atoms via Python binding...
    }
  @endcode
*/
class CLIPPER_CX_IMEX BFactorOccRefinerThread {
public:
    //! Construct with a refinement configuration.
    explicit BFactorOccRefinerThread(const RefineConfig& cfg);

    //! Destructor — waits for any running thread to complete.
    ~BFactorOccRefinerThread();

    //! Launch background refinement against crystallographic data.
    /*!
      \param atoms    Clipper Atom_list (COPIED); B-factors/occupancies will be refined.
      \param mapping  Per-entry (Atom*, altloc, index) records built by
                      clipper_atoms_from_cx_atoms_with_map (COPIED); used by
                      apply_to_atoms() to write results back to ChimeraX.
      \param fobs     Observed amplitudes HKL_data<F_sigF> (COPIED).
      \param phi_fom  ML weights HKL_data<Phi_fom> containing figure-of-merit m_h
                      and best phases (COPIED); used as fixed sigma_A proxy.
      \param usage    Working-set flags HKL_data<Flag> (COPIED); free-set
                      reflections are excluded from both T and the driving density.

      All objects are deep-copied; the caller may modify them immediately
      after this call returns.
    */
    void launch(const Atom_list&                          atoms,
                const std::vector<bridge::AtomAltlocIndex>& mapping,
                const HKL_data<F_sigF<ftype32>>&          fobs,
                const HKL_data<Phi_fom<ftype32>>&         phi_fom,
                const HKL_data<Flag>&                     usage,
                bool                                      ignore_hydrogens = false,
                const HKL_data<F_phi<ftype32>>&           f_bulk
                    = HKL_data<F_phi<ftype32>>(),
                const BFactorRestraintSpec&               restraint_spec
                    = BFactorRestraintSpec());

    //! Write refined B-factors and occupancies back to ChimeraX atoms via the
    //! stored mapping.  Must be called from the main thread.  The caller must
    //! verify atom validity before calling (e.g. via the passive length check
    //! or the structural-changes hook in BFactorOccRefineManager).
    void apply_to_atoms();

    //! Real-space LS launch: refine against a fixed P1 Xmap target density.
    //! \param refined_atoms   Atoms being optimised (mapping used for write-back).
    //! \param mapping         Per-atom altloc mapping for refined_atoms.
    //! \param context_atoms   Additional atoms that contribute to the calculated
    //!                        density but are not optimised.
    //! \param target_map      Fixed P1 Xmap containing the target density.
    //! \param target_origin   Real-space origin of the P1 cell in the original
    //!                        orthogonal coordinate frame (Å).
    void launch_realspace(const Atom_list&                            refined_atoms,
                          const std::vector<bridge::AtomAltlocIndex>& mapping,
                          const Atom_list&                            context_atoms,
                          const Xmap<ftype32>&                        target_map,
                          const Coord_orth&                           target_origin,
                          bool                                        ignore_hydrogens = false,
                          const BFactorRestraintSpec&                 restraint_spec
                              = BFactorRestraintSpec(),
                          const BFactorTargetRestraintSpec&           target_spec
                              = BFactorTargetRestraintSpec());

    //! Compute R-work and R-free from the current refined parameters by running
    //! one additional EDcalc + FFT.  Only valid for the crystallographic path;
    //! returns {-1, -1} in real-space mode.  Blocks if the thread is running.
    std::pair<double, double> compute_rfactors();

    //! Standard R-work/R-free at the INITIAL parameters (captured at the start of
    //! refine()), using the same metric as compute_rfactors().  Lets the caller
    //! compare before→after like-for-like.  {-1, -1} in real-space mode / pre-run.
    //! Blocks if the thread is running.
    std::pair<double, double> initial_rfactors();

    //! Synchronous diagnostic: compute ρ_calc on the P1 grid without running
    //! the optimiser.  Atom coordinates are shifted by -target_origin (same as
    //! operator_realspace_()) before EDcalc, so any misregistration between the
    //! computed density and the target is immediately visible when both are
    //! displayed as ChimeraX Volumes at target_origin.
    //! Returns a heap-allocated Xmap whose ownership is transferred to the
    //! pybind11 caller via py::return_value_policy::take_ownership.  The
    //! returned Xmap has the same spacegroup/cell/grid as target_map and can
    //! be passed directly to _debug_show_p1_target() in Python.
    Xmap<ftype32>* compute_realspace_density(
        const Atom_list&                            refined_atoms,
        const std::vector<bridge::AtomAltlocIndex>& mapping,
        const Atom_list&                            context_atoms,
        const Xmap<ftype32>&                        target_map,
        const Coord_orth&                           target_origin);

    //! Signal the running refinement thread to stop at the next L-BFGS-B
    //! evaluation.  Thread-safe; may be called from any thread.
    void cancel();

    //! True if the background thread has been launched and has not yet finished.
    //! Uses !ready_ as the authoritative signal rather than future::valid()
    //! alone, which can be unreliable on some MSVC versions after get().
    bool thread_running() const { return thread_result_.valid() && !ready_; }

    //! True when the thread has finished and results are ready.
    bool ready() const { return ready_; }

    //! Number of atoms in the last-launched atom list.
    int n_atoms() const;

    //! Refined isotropic U values (Å²), in the same order as the input Atom_list.
    /*! Blocks if the thread is still running. */
    const std::vector<double>& refined_u_iso();

    //! Refined occupancies, in the same order as the input Atom_list.
    /*! Blocks if the thread is still running. */
    const std::vector<double>& refined_occ();

private:
    RefineConfig                          cfg_;
    std::vector<bridge::AtomAltlocIndex>  atom_mapping_;
    bool                                  ready_            = false;
    bool                                  ignore_hydrogens_ = false;
    std::future<bool>                     thread_result_;
    std::unique_ptr<BFactorOccRefiner>    refiner_;

    bool run_();

    //! Build the Clipper-indexed config the refiner actually runs on, from the
    //! ChimeraX-indexed cfg_ + the launch-time mapping: expands refine_occ to
    //! per-Clipper, resolves restraint specs to Clipper indices, and auto-derives
    //! occupancy groups (covalent fragment + altloc).
    RefineConfig build_effective_config_(
        const std::vector<bridge::AtomAltlocIndex>& mapping,
        int n_clipper_atoms,
        const BFactorRestraintSpec&       restraint_spec,
        const BFactorTargetRestraintSpec& target_spec);
};

} // namespace clipper_cx
