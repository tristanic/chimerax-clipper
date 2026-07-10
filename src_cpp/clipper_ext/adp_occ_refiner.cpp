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

#pragma warning(disable: 4251)
#include "adp_occ_refiner.h"
#include "edcalc_ext.h"

#include <algorithm>
#include <atomic>
#include <complex>
#include <numeric>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <string>

#include <atomstruct/Bond.h>   // Atom::bonds() / Bond::other_atom() for occ-group derivation

// Eigen and LBFGSpp are included only here, keeping them out of the public header.
#include <Eigen/Core>
#include <LBFGSB.h>

namespace clipper_cx {

using namespace clipper;
using namespace clipper::datatypes;

namespace {

// Barron's general robust loss (J.T. Barron, "A General and Adaptive Robust Loss
// Function", CVPR 2019 / arXiv 2017), evaluated in U-space.  Returns the energy
// and its first derivative dE/dr for a residual r, scale c (= sigma), shape
// alpha, and overall weight.  See adp_occ_types.h for the full description of the
// family.  The removable singularities at alpha = 2 (squared error), alpha = 0
// (Cauchy) and alpha -> -inf (Welsch) are handled by dedicated branches; every
// branch is ≈ weight·½(r/c)² near r = 0 (origin curvature weight/c²).
struct LossEval { double energy; double dE_dr; };

inline LossEval barron_loss(double r, double alpha, double c, double weight)
{
    const double inv_c2 = 1.0 / (c * c);
    const double x2     = r * r * inv_c2;   // (r/c)^2
    const double g0     = weight * r * inv_c2;  // common gradient prefactor weight·r/c²
    LossEval e;

    if (alpha <= -1e6) {
        // Welsch / Leclerc limit (alpha -> -inf)
        const double ex = std::exp(-0.5 * x2);
        e.energy = weight * (1.0 - ex);
        e.dE_dr  = g0 * ex;
    } else if (std::abs(alpha) < 1e-8) {
        // Cauchy / Lorentzian limit (alpha -> 0)
        e.energy = weight * std::log(0.5 * x2 + 1.0);
        e.dE_dr  = g0 / (0.5 * x2 + 1.0);
    } else if (std::abs(alpha - 2.0) < 1e-8) {
        // Squared-error / harmonic limit (alpha -> 2)
        e.energy = weight * 0.5 * x2;
        e.dE_dr  = g0;
    } else {
        // General case
        const double b = std::abs(alpha - 2.0);
        const double z = x2 / b + 1.0;
        e.energy = weight * (b / alpha) * (std::pow(z, 0.5 * alpha) - 1.0);
        e.dE_dr  = g0 * std::pow(z, 0.5 * alpha - 1.0);
    }
    return e;
}

} // anonymous namespace

// ============================================================================
// BFactorOccRefiner
// ============================================================================
//
// Implements the L-BFGS-B functor for isotropic B-factor and/or occupancy
// refinement using the Agarwal (1978) / Ten Eyck (1977) real-space approach.
//
// The target function and its gradient are recomputed from scratch at every
// L-BFGS-B evaluation, giving a fixed, self-consistent objective:
//
//   T(U, occ)  = ½ · Σ_{h ∈ working}  (k·|Fc(h)| − m_h·|Fo(h)|)²
//
// where k is an isotropic scale recomputed each evaluation (least-squares
// scale of Fc onto m·Fo over working reflections) and m_h = phi_fom.fom().
//
// The driving density and gradient are related by:
//
//   d(x)     = FFT{ k·(m_h·|Fo(h)| − k·|Fc(h)|) · exp(iφ_calc(h)) }
//   ∂T/∂p_j  = −Σ_x  d(x) · ∂ρ_j(x)/∂p_j
//
// Using a pre-computed, fixed d(x) (the old approach) is incorrect because
// d(x) depends on the current Fcalc: fixing it across many L-BFGS steps makes
// each step minimize a different target, causing oscillation between launches.
//
// Parameter vector layout (flat, for L-BFGS-B):
//   [U_iso_0, ..., U_iso_{N-1},               (if refine_b; one per atom)
//    occ_free_0, ...,                          (unconstrained occ representatives)
//    occ_group0_p0, ..., occ_group0_p{K-2},    (K-1 free params per K-entity sum group)
//    ...]
//
// Occupancy constraints — see adp_occ_types.h for full description.
//
// Strict working/free-set separation: usage_ flags exclude free-set reflections
// from both the target value and the driving density.

class BFactorOccRefiner {
public:
    //! Crystallographic ML constructor (per-iteration Fcalc recomputation).
    BFactorOccRefiner(const Atom_list&                  atoms,
                      const HKL_data<F_sigF<ftype32>>&  fobs,
                      const HKL_data<Phi_fom<ftype32>>& phi_fom,
                      const HKL_data<Flag>&             usage,
                      const RefineConfig&               cfg,
                      const HKL_data<F_phi<ftype32>>&   f_bulk
                          = HKL_data<F_phi<ftype32>>());

    //! Real-space LS constructor (fixed P1 target density).
    //! \param refined_atoms  Atoms whose parameters are optimised.
    //! \param context_atoms  Additional atoms that contribute to the
    //!                       calculated density but are not optimised.
    //! \param target_map     P1 Xmap containing the fixed target density.
    //!                       Must have the same spacegroup/cell/grid that
    //!                       density is to be calculated on.
    //! \param target_origin  Real-space origin of the P1 cell in the
    //!                       original orthogonal coordinate frame (Å).
    //!                       Atom positions are shifted by -target_origin
    //!                       before EDcalc so that they fall inside the P1 cell.
    BFactorOccRefiner(const Atom_list&     refined_atoms,
                      const Atom_list&     context_atoms,
                      const Xmap<ftype32>& target_map,
                      const Coord_orth&    target_origin,
                      const RefineConfig&  cfg);

    // L-BFGS-B functor interface.
    double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad);

    // Run the L-BFGS-B optimisation. Returns true on normal completion.
    bool refine();

    int n_atoms() const { return n_atoms_; }
    const std::vector<double>& refined_u_iso() const { return current_u_iso_; }
    const std::vector<double>& refined_occ()   const { return current_occ_; }

    //! Standard R-work/R-free at the INITIAL parameters, captured at the start of
    //! refine() (before optimisation).  Same metric as compute_rfactors(), so the
    //! before→after comparison in the manager's bail-out is exactly like-for-like.
    //! {-1, -1} in real-space mode or if refine() has not run.
    std::pair<double, double> initial_rfactors() const { return {rwork_before_, rfree_before_}; }

    //! Signal the refinement to stop at the next functor evaluation.
    void cancel() { cancel_.store(true, std::memory_order_relaxed); }

    std::pair<double, double> compute_rfactors();

private:
    std::atomic<bool> cancel_{false};
    bool              realspace_mode_ = false;

    // Fixed input data (copied at construction, constant through run)
    Atom_list                   atoms_;     // refined atoms
    Atom_list                   context_atoms_; // real-space only: fixed context
    // Crystallographic ML only:
    HKL_data<F_sigF<ftype32>>   fobs_;
    HKL_data<Phi_fom<ftype32>>  phi_fom_;
    HKL_data<Flag>              usage_;
    Cell                        cell_;
    Grid_sampling               grid_sampling_;
    // Crystallographic: bulk-solvent contribution F_bulk, supplied by the caller
    // (Xtal_mgr::f_bulk) and added to F_atoms(current U) in operator().
    HKL_data<F_phi<ftype32>>    fbulk_hkl_;
    bool                        has_bulk_solvent_ = false;
    // Standard R-work/R-free at the initial params, set at the start of refine().
    double                      rwork_before_ = -1.0;
    double                      rfree_before_ = -1.0;
    // Real-space only:
    Xmap<ftype32>               xmap_target_;   // fixed P1 target density
    Coord_orth                  target_origin_; // origin of P1 cell in original frame

    RefineConfig                cfg_;

    // Working space — reused each operator() call (single-threaded L-BFGS-B)
    Xmap<ftype32>               density_xmap_;   // EDcalc output; then d(x)
    HKL_data<F_phi<ftype32>>    fcalc_hkl_;      // Fc from forward FFT
    HKL_data<F_phi<ftype32>>    driving_hkl_;    // driving density coefficients
    EDcalc_aniso_thread<ftype32> edcalc_;

    // Derived from cfg_ and atoms_ during construction
    int  n_atoms_         = 0;
    int  n_params_        = 0;
    bool any_occ_refined_ = false;

    std::vector<int>               occ_representative_;
    std::vector<int>               free_occ_atoms_;
    std::vector<AtomShapeFn::TYPE> agarwal_types_;
    int uiso_grad_idx_ = -1;
    int  occ_grad_idx_ = -1;

    std::vector<double> current_u_iso_;
    std::vector<double> current_occ_;
    std::vector<double> atom_radii_;

    void init_derived();
    double operator_realspace_(const Eigen::VectorXd& x, Eigen::VectorXd& grad);
    void atoms_to_params(Eigen::VectorXd& x) const;
    void decode_params(const Eigen::VectorXd& x);
    void build_bounds(Eigen::VectorXd& lb, Eigen::VectorXd& ub) const;
    void assemble_occ_gradient(const std::vector<double>& raw_occ_grad,
                                Eigen::VectorXd& grad, int start_idx) const;
    bool group_is_refined(const OccConstraintGroup& grp) const
    {
        if (!any_occ_refined_ || grp.atom_indices.empty()) return false;
        int idx0 = grp.atom_indices[0];
        if (idx0 < 0 || idx0 >= n_atoms_)                  return false;
        return cfg_.refine_occ[occ_representative_[idx0]] != 0;
    }
};

// ----------------------------------------------------------------------------

BFactorOccRefiner::BFactorOccRefiner(const Atom_list&                  atoms,
                                     const HKL_data<F_sigF<ftype32>>&  fobs,
                                     const HKL_data<Phi_fom<ftype32>>& phi_fom,
                                     const HKL_data<Flag>&             usage,
                                     const RefineConfig&               cfg,
                                     const HKL_data<F_phi<ftype32>>&   f_bulk)
    : atoms_(atoms), fobs_(fobs), phi_fom_(phi_fom), usage_(usage), cfg_(cfg)
{
    // Bulk-solvent contribution F_bulk, supplied directly by the crystal manager
    // (Xtal_mgr::f_bulk).  A default-constructed HKL_data has num_obs()==0 and
    // signals "no bulk solvent".  Held fixed across the refinement (the atoms do
    // not move), so F_total = F_atoms(current U) + F_bulk at every step.
    if (f_bulk.num_obs() > 0) {
        fbulk_hkl_        = f_bulk;
        has_bulk_solvent_ = true;
    }
    init_derived();
}

BFactorOccRefiner::BFactorOccRefiner(const Atom_list&     refined_atoms,
                                     const Atom_list&     context_atoms,
                                     const Xmap<ftype32>& target_map,
                                     const Coord_orth&    target_origin,
                                     const RefineConfig&  cfg)
    : atoms_(refined_atoms), context_atoms_(context_atoms),
      target_origin_(target_origin), cfg_(cfg),
      realspace_mode_(true)
{
    // Copy target Xmap: Xmap copy-assignment is deleted so use init() + grid loop.
    xmap_target_.init(target_map.spacegroup(), target_map.cell(),
                      target_map.grid_sampling());
    for (Xmap<ftype32>::Map_reference_index ix = target_map.first();
         !ix.last(); ix.next())
        xmap_target_.set_data(ix.coord(), target_map[ix]);
    init_derived();
}

void BFactorOccRefiner::init_derived()
{
    n_atoms_ = (int)atoms_.size();

    current_u_iso_.resize(n_atoms_);
    current_occ_.resize(n_atoms_);
    atom_radii_.resize(n_atoms_);

    for (int j = 0; j < n_atoms_; ++j) {
        const Atom& a = atoms_[j];
        if (a.is_null()) {
            current_u_iso_[j] = Util::b2u(cfg_.b_min);
            current_occ_[j]   = 0.0;
            atom_radii_[j]    = 3.5;
            continue;
        }
        current_u_iso_[j] = std::max(a.u_iso(), Util::b2u(cfg_.b_min));
        current_occ_[j]   = a.occupancy();
        double u_clamped = std::min(std::max(a.u_iso(), 0.0), Util::b2u(200.0));
        atom_radii_[j] = std::max(1.5 + 3.0 * std::sqrt(u_clamped), 3.5);
    }

    // Real-space: bulk-shift the starting B-factors so the fragment mean matches
    // the context mean.  A newly-fitted fragment often carries a large systematic
    // B-factor offset relative to its surroundings; removing it here gives the
    // boundary restraints (and the optimiser) a near-correct starting point,
    // leaving only the genuine per-atom variation to refine.
    // (atoms_ and context_atoms_ are already H-filtered at construction when
    // ignore_hydrogens is set, so no element test is needed here.)
    if (realspace_mode_ && !context_atoms_.empty()) {
        double sum_frag = 0.0; int n_frag = 0;
        for (int j = 0; j < n_atoms_; ++j) {
            if (atoms_[j].is_null()) continue;
            sum_frag += current_u_iso_[j];
            ++n_frag;
        }
        double sum_ctx = 0.0; int n_ctx = 0;
        for (size_t c = 0; c < context_atoms_.size(); ++c) {
            const Atom& ca = context_atoms_[c];
            if (ca.is_null()) continue;
            sum_ctx += ca.u_iso();
            ++n_ctx;
        }
        if (n_frag > 0 && n_ctx > 0) {
            double shift = sum_ctx / n_ctx - sum_frag / n_frag;
            double u_min = Util::b2u(cfg_.b_min);
            double u_max = Util::b2u(cfg_.b_max);
            for (int j = 0; j < n_atoms_; ++j)
                current_u_iso_[j] = std::min(std::max(current_u_iso_[j] + shift,
                                                      u_min), u_max);
        }
    }

    // Working space — mode-dependent initialisation.
    // Xmap copy-assignment is deleted (unique_ptr<atomic_flag[]> member),
    // so use init() to populate the default-constructed member.
    if (realspace_mode_) {
        // Real-space: density_xmap_ shares grid with the P1 target map.
        density_xmap_.init(xmap_target_.spacegroup(), xmap_target_.cell(),
                           xmap_target_.grid_sampling());
        // fcalc_hkl_ / driving_hkl_ are not used in real-space mode.
    } else {
        // Crystallographic: derive grid from the HKL data.
        const HKL_info& hkls = fobs_.base_hkl_info();
        cell_          = fobs_.base_cell();
        grid_sampling_ = Grid_sampling(hkls.spacegroup(), cell_, hkls.resolution());
        density_xmap_.init(hkls.spacegroup(), cell_, grid_sampling_);
        fcalc_hkl_   = HKL_data<F_phi<ftype32>>(hkls, cell_);
        driving_hkl_ = HKL_data<F_phi<ftype32>>(hkls, cell_);
    }
    edcalc_ = EDcalc_aniso_thread<ftype32>((size_t)cfg_.n_threads);

    // Validate and derive any_occ_refined_
    if (!cfg_.refine_occ.empty()) {
        if ((int)cfg_.refine_occ.size() != n_atoms_)
            throw std::runtime_error(
                "BFactorOccRefiner: refine_occ array length ("
                + std::to_string(cfg_.refine_occ.size())
                + ") does not match atom count ("
                + std::to_string(n_atoms_) + ")");
        any_occ_refined_ = std::any_of(cfg_.refine_occ.begin(), cfg_.refine_occ.end(),
                                        [](uint8_t v){ return v != 0; });
    }

    // Representative map for EqualOccGroups
    occ_representative_.resize(n_atoms_);
    std::iota(occ_representative_.begin(), occ_representative_.end(), 0);
    for (const auto& grp : cfg_.equal_occ_groups) {
        if (grp.atom_indices.empty()) continue;
        int rep = grp.atom_indices[0];
        for (int idx : grp.atom_indices)
            occ_representative_[idx] = rep;
    }

    for (const auto& grp : cfg_.equal_occ_groups) {
        if (grp.atom_indices.empty()) continue;
        double shared = current_occ_[grp.atom_indices[0]];
        for (int idx : grp.atom_indices)
            current_occ_[idx] = shared;
    }

    std::set<int> grouped_rep_set;
    for (const auto& grp : cfg_.occ_groups)
        for (int idx : grp.atom_indices)
            grouped_rep_set.insert(occ_representative_[idx]);

    if (any_occ_refined_) {
        for (int j = 0; j < n_atoms_; ++j)
            if (cfg_.refine_occ[j]
                && occ_representative_[j] == j
                && grouped_rep_set.find(j) == grouped_rep_set.end())
                free_occ_atoms_.push_back(j);
    }

    if (cfg_.refine_b) {
        uiso_grad_idx_ = (int)agarwal_types_.size();
        agarwal_types_.push_back(AtomShapeFn::Uiso);
    }
    if (any_occ_refined_) {
        occ_grad_idx_ = (int)agarwal_types_.size();
        agarwal_types_.push_back(AtomShapeFn::Occ);
    }

    n_params_ = 0;
    if (cfg_.refine_b)
        n_params_ += n_atoms_;
    if (any_occ_refined_) {
        n_params_ += (int)free_occ_atoms_.size();
        for (const auto& grp : cfg_.occ_groups)
            if (group_is_refined(grp))
                n_params_ += (int)grp.atom_indices.size() - 1;
    }

    if (n_params_ == 0)
        throw std::runtime_error(
            "BFactorOccRefiner: no parameters to refine "
            "(refine_b is false and refine_occ is empty or all-zero)");
}

void BFactorOccRefiner::atoms_to_params(Eigen::VectorXd& x) const
{
    x.resize(n_params_);
    int idx = 0;
    if (cfg_.refine_b)
        for (int j = 0; j < n_atoms_; ++j)
            x[idx++] = current_u_iso_[j];
    if (any_occ_refined_) {
        for (int j : free_occ_atoms_)
            x[idx++] = current_occ_[j];
        for (const auto& grp : cfg_.occ_groups) {
            if (!group_is_refined(grp)) continue;
            for (size_t k = 0; k + 1 < grp.atom_indices.size(); ++k)
                x[idx++] = current_occ_[grp.atom_indices[k]];
        }
    }
}

void BFactorOccRefiner::decode_params(const Eigen::VectorXd& x)
{
    int idx = 0;
    if (cfg_.refine_b)
        for (int j = 0; j < n_atoms_; ++j)
            current_u_iso_[j] = x[idx++];
    if (any_occ_refined_) {
        for (int j : free_occ_atoms_)
            current_occ_[j] = x[idx++];
        for (const auto& grp : cfg_.occ_groups) {
            if (!group_is_refined(grp)) continue;
            double sum = 0.0;
            for (size_t k = 0; k + 1 < grp.atom_indices.size(); ++k) {
                double o = x[idx++];
                current_occ_[occ_representative_[grp.atom_indices[k]]] = o;
                sum += o;
            }
            current_occ_[occ_representative_[grp.atom_indices.back()]] = grp.total - sum;
        }
        for (const auto& grp : cfg_.equal_occ_groups) {
            if (grp.atom_indices.empty()) continue;
            double shared = current_occ_[grp.atom_indices[0]];
            for (int idx2 : grp.atom_indices)
                current_occ_[idx2] = shared;
        }
    }
}

void BFactorOccRefiner::build_bounds(Eigen::VectorXd& lb, Eigen::VectorXd& ub) const
{
    lb.resize(n_params_);
    ub.resize(n_params_);
    int idx = 0;
    if (cfg_.refine_b) {
        double u_min = Util::b2u(cfg_.b_min);
        double u_max = Util::b2u(cfg_.b_max);
        for (int j = 0; j < n_atoms_; ++j) {
            lb[idx] = u_min;
            ub[idx] = u_max;
            ++idx;
        }
    }
    if (any_occ_refined_) {
        for (size_t i = 0; i < free_occ_atoms_.size(); ++i) {
            lb[idx] = 0.0;
            ub[idx] = 1.0;
            ++idx;
        }
        for (const auto& grp : cfg_.occ_groups) {
            if (!group_is_refined(grp)) continue;
            for (size_t k = 0; k + 1 < grp.atom_indices.size(); ++k) {
                lb[idx] = 0.0;
                ub[idx] = grp.total;
                ++idx;
            }
        }
    }
}

void BFactorOccRefiner::assemble_occ_gradient(
    const std::vector<double>& raw_occ_grad,
    Eigen::VectorXd& grad, int start_idx) const
{
    int idx = start_idx;
    for (int j : free_occ_atoms_)
        grad[idx++] = raw_occ_grad[j];
    for (const auto& grp : cfg_.occ_groups) {
        if (!group_is_refined(grp)) continue;
        int last_rep = occ_representative_[grp.atom_indices.back()];
        for (size_t k = 0; k + 1 < grp.atom_indices.size(); ++k) {
            int rep_k = occ_representative_[grp.atom_indices[k]];
            grad[idx++] = raw_occ_grad[rep_k] - raw_occ_grad[last_rep];
        }
    }
}

double BFactorOccRefiner::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
{
    if (cancel_.load(std::memory_order_relaxed))
        throw std::runtime_error("BFactorOccRefiner: cancelled");

    if (realspace_mode_)
        return operator_realspace_(x, grad);

    decode_params(x);
    grad.setZero(n_params_);

    // -------------------------------------------------------------------------
    // 1. Build atom list with current U_iso and occupancy values.
    // -------------------------------------------------------------------------
    Atom_list atomu = atoms_;
    for (int j = 0; j < n_atoms_; ++j) {
        if (atomu[j].is_null()) continue;
        if (atomu[j].u_aniso_orth().is_null())
            atomu[j].set_u_iso((ftype)current_u_iso_[j]);
        atomu[j].set_occupancy((ftype)current_occ_[j]);
    }

    // -------------------------------------------------------------------------
    // 2. Compute electron density and forward-FFT to get Fcalc.
    // -------------------------------------------------------------------------
    density_xmap_ = ftype32(0);
    edcalc_(density_xmap_, atomu);
    density_xmap_.fft_to(fcalc_hkl_, cfg_.n_threads);

    // -------------------------------------------------------------------------
    // Add bulk-solvent contribution: F_total = F_atoms + F_bulk.
    // -------------------------------------------------------------------------
    if (has_bulk_solvent_) {
        for (HKL_info::HKL_reference_index ih = fobs_.first(); !ih.last(); ih.next()) {
            if (fcalc_hkl_[ih].missing() || fbulk_hkl_[ih].missing()) continue;
            std::complex<ftype32> f_total =
                (std::complex<ftype32>)fcalc_hkl_[ih]
              + (std::complex<ftype32>)fbulk_hkl_[ih];
            fcalc_hkl_.set_data(ih.hkl(), F_phi<ftype32>(f_total));
        }
    }

    // -------------------------------------------------------------------------
    // 3. Compute isotropic scale k = Σ(m|Fo||Fc|) / Σ|Fc|²  (working set only).
    //    usage_ flag convention (set in xtal_mgr): BOTH (!=0) = working,
    //    NONE (==0) = free.  The free set must be excluded here and from the
    //    target/gradient below, or it leaks into the fit (R-free ≈ R-work).
    // -------------------------------------------------------------------------
    ftype sum_cross = 0.0, sum_sq = 0.0;
    for (HKL_info::HKL_reference_index ih = fobs_.first(); !ih.last(); ih.next()) {
        if (fobs_[ih].missing() || usage_[ih].missing()
            || usage_[ih].flag() == 0) continue;
        ftype fc = fcalc_hkl_[ih].f();
        ftype fo = fobs_[ih].f();
        ftype m  = phi_fom_[ih].fom();
        sum_cross += m * fo * fc;
        sum_sq    += fc * fc;
    }
    ftype k = (sum_sq > ftype(1e-10)) ? sum_cross / sum_sq : ftype(1.0);

    // -------------------------------------------------------------------------
    // 4. Compute target T and driving density coefficients for all reflections.
    //    Working reflections: G(h) = k·(m|Fo| − k|Fc|)·exp(iφ_calc).
    //    Non-working (missing OR free, flag==0): set to null so fft_from treats
    //    them as zero — the free set contributes neither to T nor to the gradient.
    // -------------------------------------------------------------------------
    double T = 0.0;
    for (HKL_info::HKL_reference_index ih = fobs_.first(); !ih.last(); ih.next()) {
        if (fobs_[ih].missing() || usage_[ih].missing()
            || usage_[ih].flag() == 0) {
            driving_hkl_.set_data(ih.hkl(), F_phi<ftype32>());
            continue;
        }
        ftype fc_mag    = fcalc_hkl_[ih].f();
        ftype phi_c     = fcalc_hkl_[ih].phi();
        ftype fo        = fobs_[ih].f();
        ftype m         = phi_fom_[ih].fom();
        ftype fc_scaled = k * fc_mag;
        ftype residual  = fc_scaled - m * fo;
        T += ftype(0.5) * residual * residual;

        // d(h) = ε_c·k·(m|Fo| − k|Fc|)·exp(iφ_calc):  consistent with ∂T/∂p = −Σ d·∂ρ/∂p.
        // The per-reflection multiplicity weight ε_c = hkl_class().epsilonc()
        // (= 2·epsilon centric, epsilon acentric) compensates for fft_from spreading
        // each ASU reflection over its orbit (symmetry mates + Friedel) via set_hkl:
        // centric/axis reflections populate fewer distinct grid points, so their
        // single-count coefficient must be boosted or the back-FFT gradient is the
        // gradient of the WRONG target (T with centric reflections silently
        // under-weighted). T itself stays the plain unweighted ASU sum above; only
        // the driving coefficient carries ε_c, making this the exact gradient of T.
        // ε_c = 1 for every reflection in P1 (unchanged there). See xray_gradient.cpp.
        ftype eps_c = ih.hkl_class().epsilonc();
        std::complex<ftype32> coeff =
            ftype32(eps_c * k) * ftype32(m * fo - fc_scaled)
            * std::polar(ftype32(1.0f), ftype32(phi_c));
        driving_hkl_.set_data(ih.hkl(), F_phi<ftype32>(coeff));
    }

    // -------------------------------------------------------------------------
    // 5. Back-FFT driving coefficients → d(x) in real space.
    //    density_xmap_ is now the driving density, ready for rho_grad.
    // -------------------------------------------------------------------------
    density_xmap_.fft_from(driving_hkl_, cfg_.n_threads);

    // -------------------------------------------------------------------------
    // 6. Accumulate per-atom gradients from d(x) — parallelised over atoms.
    //
    // Each atom j writes only to raw_grad_u[j] and raw_grad_occ[j], so the
    // atom chunks are non-overlapping and require no locking.  density_xmap_
    // is read-only at this point.  rho_grad_vals is declared inside the
    // lambda so each worker owns its own thread-local copy.
    // -------------------------------------------------------------------------
    std::vector<double> raw_grad_u(n_atoms_, 0.0);
    std::vector<double> raw_grad_occ(n_atoms_, 0.0);

    const Cell&          cell = density_xmap_.cell();
    const Grid_sampling& grid = density_xmap_.grid_sampling();

    auto gradient_worker = [&](int start, int end) {
        std::vector<ftype>  rho_grad_vals;
        using MapRef = Xmap<ftype32>::Map_reference_coord;
        ftype rho;

        for (int j = start; j < end; ++j) {
            const Atom& a = atoms_[j];
            if (a.is_null()) continue;

            AtomShapeFn sf(a.coord_orth(), a.element(),
                           (ftype)current_u_iso_[j], (ftype)current_occ_[j]);
            sf.agarwal_params() = agarwal_types_;

            Grid_range gd(cell, grid, (ftype)atom_radii_[j]);
            Coord_grid cg = density_xmap_.coord_map(a.coord_orth()).coord_grid();
            Coord_grid g0 = cg + gd.min();
            Coord_grid g1 = cg + gd.max();

            MapRef iu(density_xmap_, g0);
            for (; iu.coord().u() <= g1.u(); iu.next_u()) {
                MapRef iv = iu;
                for (; iv.coord().v() <= g1.v(); iv.next_v()) {
                    MapRef iw = iv;
                    for (; iw.coord().w() <= g1.w(); iw.next_w()) {
                        sf.rho_grad(iw.coord_orth(), rho, rho_grad_vals);
                        double d = (double)density_xmap_[iw];
                        if (uiso_grad_idx_ >= 0)
                            raw_grad_u[j]   -= d * (double)rho_grad_vals[uiso_grad_idx_];
                        if (occ_grad_idx_  >= 0 && cfg_.refine_occ[j])
                            raw_grad_occ[j] -= d * (double)rho_grad_vals[occ_grad_idx_];
                    }
                }
            }
        }
    };

    int n_t   = std::max(1, cfg_.n_threads);
    int chunk = (n_atoms_ + n_t - 1) / n_t;   // ceiling division
    std::vector<std::future<void>> futures;
    futures.reserve(n_t);
    for (int t = 0; t < n_t; ++t) {
        int start = t * chunk;
        int end   = std::min(start + chunk, n_atoms_);
        if (start >= n_atoms_) break;
        futures.push_back(
            std::async(std::launch::async, gradient_worker, start, end));
    }
    for (auto& f : futures)
        f.get();

    // Sum EqualOccGroup member gradients onto their representatives.
    if (occ_grad_idx_ >= 0) {
        for (const auto& grp : cfg_.equal_occ_groups) {
            if (grp.atom_indices.empty()) continue;
            double sum = 0.0;
            for (int idx : grp.atom_indices)
                if (cfg_.refine_occ[idx])
                    sum += raw_grad_occ[idx];
            raw_grad_occ[grp.atom_indices[0]] = sum;
        }
    }

    // -------------------------------------------------------------------------
    // 7. B-factor restraints (Barron general loss, in U-space).
    // -------------------------------------------------------------------------
    for (const auto& r : cfg_.b_restraints) {
        double du = current_u_iso_[r.i] - current_u_iso_[r.j];
        LossEval e = barron_loss(du, r.alpha, r.sigma, r.weight);
        T += e.energy;
        raw_grad_u[r.i] += e.dE_dr;
        raw_grad_u[r.j] -= e.dE_dr;
    }

    // Assemble flat gradient vector.
    int idx = 0;
    if (cfg_.refine_b)
        for (int j = 0; j < n_atoms_; ++j)
            grad[idx++] = raw_grad_u[j];
    if (any_occ_refined_)
        assemble_occ_gradient(raw_grad_occ, grad, idx);

    return T;
}

double BFactorOccRefiner::operator_realspace_(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
{
    decode_params(x);
    grad.setZero(n_params_);

    // -------------------------------------------------------------------------
    // 1. Build combined atom list: context (fixed) + refined (current params).
    //    All positions are shifted to the P1 cell frame by subtracting the
    //    target_origin_ so that atoms lie within the P1 cell.
    // -------------------------------------------------------------------------
    Atom_list atomu_all;

    auto shift_atom = [&](const Atom& a, ftype u_iso, ftype occ) {
        Atom s = a;
        s.set_coord_orth(a.coord_orth() - target_origin_);
        if (s.u_aniso_orth().is_null())
            s.set_u_iso(u_iso);
        s.set_occupancy(occ);
        return s;
    };

    // Context atoms — fixed parameters
    for (size_t c = 0; c < context_atoms_.size(); ++c) {
        const Atom& ca = context_atoms_[c];
        if (!ca.is_null())
            atomu_all.push_back(shift_atom(ca, ca.u_iso(), ca.occupancy()));
    }

    // Refined atoms — current parameters
    for (int j = 0; j < n_atoms_; ++j) {
        const Atom& a = atoms_[j];
        if (!a.is_null())
            atomu_all.push_back(
                shift_atom(a, (ftype)current_u_iso_[j], (ftype)current_occ_[j]));
    }

    // -------------------------------------------------------------------------
    // 2. Compute ρ_calc on the P1 grid (same EDcalc machinery as crystallographic).
    // -------------------------------------------------------------------------
    density_xmap_ = ftype32(0);
    edcalc_(density_xmap_, atomu_all);

    // -------------------------------------------------------------------------
    // 3. Compute isotropic scale k = Σ(ρ_target·ρ_calc) / Σ ρ_calc²
    //    (analogous to the crystallographic scale of Fc onto m|Fo|).
    //
    //    Without this, the EDcalc density (e/Å³) and the cryo-EM map
    //    (arbitrary normalised units) are on completely different scales.
    //    A scale mismatch drives ρ_calc − ρ_target to a large constant offset
    //    at every atom centre, producing a gradient that saturates B-factors
    //    against the b_max bound.  With k recomputed every evaluation the
    //    target function ½Σ(k·ρ_calc − ρ_target)² is purely shape-based and
    //    converges to the correct B-factors regardless of map normalisation.
    // -------------------------------------------------------------------------
    double sum_cross = 0.0, sum_sq = 0.0;
    for (Xmap<ftype32>::Map_reference_index ix = density_xmap_.first();
         !ix.last(); ix.next())
    {
        double rc = density_xmap_[ix];
        double rt = xmap_target_[ix];
        sum_cross += rt * rc;
        sum_sq    += rc * rc;
    }
    double k = (sum_sq > 1e-20) ? sum_cross / sum_sq : 1.0;

    // -------------------------------------------------------------------------
    // 4. T = ½ Σ (k·ρ_calc − ρ_target)²;
    //    overwrite density_xmap_ with gradient map D(x) = k·(ρ_target − k·ρ_calc).
    //    gradient:  ∂T/∂p_j = Σ_x k(k·ρ_calc − ρ_target)·∂ρ_j/∂p_j
    //                        = −Σ_x D(x)·∂ρ_j/∂p_j   (same form as crystal path)
    // -------------------------------------------------------------------------
    double T = 0.0;
    for (Xmap<ftype32>::Map_reference_index ix = density_xmap_.first();
         !ix.last(); ix.next())
    {
        double rc   = density_xmap_[ix];          // ρ_calc (before overwrite)
        double rt   = xmap_target_[ix];
        double diff = k * rc - rt;                // k·ρ_calc − ρ_target
        T += 0.5 * diff * diff;
        density_xmap_[ix] = (ftype)(k * (rt - k * rc));  // D(x)
    }

    // -------------------------------------------------------------------------
    // 4. Per-atom gradients — refined atoms only, using shifted positions.
    // -------------------------------------------------------------------------
    std::vector<double> raw_grad_u(n_atoms_, 0.0);
    std::vector<double> raw_grad_occ(n_atoms_, 0.0);
    std::vector<ftype>  rho_grad_vals;

    const Cell&          cell = density_xmap_.cell();
    const Grid_sampling& grid = density_xmap_.grid_sampling();

    for (int j = 0; j < n_atoms_; ++j) {
        const Atom& a = atoms_[j];
        if (a.is_null()) continue;

        // Positions in the P1 cell frame
        Coord_orth pos_p1 = a.coord_orth() - target_origin_;

        AtomShapeFn sf(pos_p1, a.element(),
                       (ftype)current_u_iso_[j], (ftype)current_occ_[j]);
        sf.agarwal_params() = agarwal_types_;

        Grid_range gd(cell, grid, (ftype)atom_radii_[j]);
        Coord_grid cg = density_xmap_.coord_map(pos_p1).coord_grid();
        Coord_grid g0 = cg + gd.min();
        Coord_grid g1 = cg + gd.max();

        using MapRef = Xmap<ftype32>::Map_reference_coord;
        ftype rho;

        MapRef iu(density_xmap_, g0);
        for (; iu.coord().u() <= g1.u(); iu.next_u()) {
            MapRef iv = iu;
            for (; iv.coord().v() <= g1.v(); iv.next_v()) {
                MapRef iw = iv;
                for (; iw.coord().w() <= g1.w(); iw.next_w()) {
                    sf.rho_grad(iw.coord_orth(), rho, rho_grad_vals);
                    double d = (double)density_xmap_[iw];
                    if (uiso_grad_idx_ >= 0)
                        raw_grad_u[j]   -= d * (double)rho_grad_vals[uiso_grad_idx_];
                    if (occ_grad_idx_  >= 0 && cfg_.refine_occ[j])
                        raw_grad_occ[j] -= d * (double)rho_grad_vals[occ_grad_idx_];
                }
            }
        }
    }

    // Sum EqualOccGroup member gradients onto representatives.
    if (occ_grad_idx_ >= 0) {
        for (const auto& grp : cfg_.equal_occ_groups) {
            if (grp.atom_indices.empty()) continue;
            double sum = 0.0;
            for (int idx : grp.atom_indices)
                if (cfg_.refine_occ[idx])
                    sum += raw_grad_occ[idx];
            raw_grad_occ[grp.atom_indices[0]] = sum;
        }
    }

    // B-factor restraints (Barron general loss, U space) — identical to crystallographic.
    for (const auto& r : cfg_.b_restraints) {
        double du = current_u_iso_[r.i] - current_u_iso_[r.j];
        LossEval e = barron_loss(du, r.alpha, r.sigma, r.weight);
        T += e.energy;
        raw_grad_u[r.i] += e.dE_dr;
        raw_grad_u[r.j] -= e.dE_dr;
    }

    // One-sided Barron restraints toward fixed target U values, supplied by the
    // caller via RefineConfig.b_target_restraints.  Only atom i's gradient is
    // accumulated (the target is a fixed value).
    for (const auto& r : cfg_.b_target_restraints) {
        double du = current_u_iso_[r.i] - r.target_u;
        LossEval e = barron_loss(du, r.alpha, r.sigma, r.weight);
        T += e.energy;
        raw_grad_u[r.i] += e.dE_dr;
    }

    // Assemble flat gradient vector.
    int idx = 0;
    if (cfg_.refine_b)
        for (int j = 0; j < n_atoms_; ++j)
            grad[idx++] = raw_grad_u[j];
    if (any_occ_refined_)
        assemble_occ_gradient(raw_grad_occ, grad, idx);

    return T;
}

bool BFactorOccRefiner::refine()
{
    // F_bulk is supplied directly by the caller (Xtal_mgr::f_bulk) and stored in
    // the constructor — no reconstruction needed here.  fbulk_hkl_ is added to
    // F_atoms(current U) inside operator() when has_bulk_solvent_ is set.

    // Record the standard R-factors at the INITIAL parameters (current_u_iso_ still
    // hold the input values here), so the manager can compare before→after using the
    // same metric.  Crystallographic path only; compute_rfactors() returns {-1,-1}
    // in real-space mode.
    if (!realspace_mode_) {
        auto rf0 = compute_rfactors();
        rwork_before_ = rf0.first;
        rfree_before_ = rf0.second;
    }

    Eigen::VectorXd x, lb, ub;
    atoms_to_params(x);
    build_bounds(lb, ub);

    LBFGSpp::LBFGSBParam<double> param;
    param.epsilon        = cfg_.lbfgs_epsilon;
    param.max_iterations = cfg_.max_cycles;
    param.past           = cfg_.lbfgs_past;
    param.delta          = cfg_.lbfgs_delta;

    LBFGSpp::LBFGSBSolver<double> solver(param);
    double fx;
    try {
        solver.minimize(*this, x, fx, lb, ub);
    } catch (const std::runtime_error&) {
        // Iteration limit reached — use the current best parameters.
    }

    decode_params(x);
    return true;
}

// ============================================================================
// BFactorOccRefinerThread
// ============================================================================

BFactorOccRefinerThread::BFactorOccRefinerThread(const RefineConfig& cfg)
    : cfg_(cfg) {}

BFactorOccRefinerThread::~BFactorOccRefinerThread()
{
    if (thread_result_.valid())
        thread_result_.get();
}

void BFactorOccRefinerThread::launch(
    const Atom_list&                            atoms,
    const std::vector<bridge::AtomAltlocIndex>& mapping,
    const HKL_data<F_sigF<ftype32>>&            fobs,
    const HKL_data<Phi_fom<ftype32>>&           phi_fom,
    const HKL_data<Flag>&                       usage,
    bool                                        ignore_hydrogens,
    const HKL_data<F_phi<ftype32>>&             f_bulk,
    const BFactorRestraintSpec&                 restraint_spec)
{
    ready_             = false;
    ignore_hydrogens_  = ignore_hydrogens;
    atom_mapping_      = mapping;
    // Resolve ChimeraX-indexed restraint specs / refine_occ flags to the
    // Clipper Atom_list indexing the refiner uses (no target restraints in the
    // crystallographic path).
    RefineConfig eff = build_effective_config_(mapping, (int)atoms.size(),
                                               restraint_spec,
                                               BFactorTargetRestraintSpec());
    refiner_.reset(new BFactorOccRefiner(atoms, fobs, phi_fom, usage, eff,
                                          f_bulk));
    thread_result_ = std::async(std::launch::async,
                                &BFactorOccRefinerThread::run_, this);
}

void BFactorOccRefinerThread::apply_to_atoms()
{
    if (!refiner_) return;
    const auto& u_iso = refiner_->refined_u_iso();
    const auto& occ   = refiner_->refined_occ();
    // B = 8π² × U  (π² ≈ 9.8696)
    constexpr double u_to_b = 8.0 * 9.8696044010893586;

    for (const auto& e : atom_mapping_) {
        auto* a = e.cx_atom;
        int   i = e.index;
        // Save the atom's current altloc, switch to the target, set values,
        // then restore — prevents leaving atoms in a non-default altloc state.
        char saved_a = '\0';
        if (e.altloc != '\0') {
            saved_a = a->alt_loc();
            a->set_alt_loc(e.altloc);
        }
        float new_b   = (float)(u_iso[i] * u_to_b);
        float new_occ = occ.empty() ? 0.0f : (float)occ[i];
        a->set_bfactor(new_b);
        if (!occ.empty())
            a->set_occupancy(new_occ);
        if (e.altloc != '\0')
            a->set_alt_loc(saved_a);   // restore

        // Propagate B-factor and occupancy to bonded H atoms.
        // Only when ignore_hydrogens_ is set: if H atoms were included in the
        // optimisation they already have their own refined values in atom_mapping_.
        if (ignore_hydrogens_) {
            atomstruct::Structure* parent_struct = a->structure();
            for (auto* bond : a->bonds()) {
                auto* h = bond->other_atom(a);
                if (h->element().number() != 1)     continue;
                if (h->structure() != parent_struct) continue;  // skip symm copies
                if (e.altloc != '\0') {
                    // Only update the H altloc matching the parent's altloc.
                    if (h->alt_locs().find(e.altloc) == h->alt_locs().end())
                        continue;
                    char saved_h = h->alt_loc();
                    h->set_alt_loc(e.altloc);
                    h->set_bfactor(new_b);
                    if (!occ.empty()) h->set_occupancy(new_occ);
                    h->set_alt_loc(saved_h);   // restore
                } else {
                    // Single-conformer parent: update H's current (default) altloc.
                    h->set_bfactor(new_b);
                    if (!occ.empty()) h->set_occupancy(new_occ);
                }
            }
        }
    }
}

void BFactorOccRefinerThread::launch_realspace(
    const Atom_list&                            refined_atoms,
    const std::vector<bridge::AtomAltlocIndex>& mapping,
    const Atom_list&                            context_atoms,
    const Xmap<ftype32>&                        target_map,
    const Coord_orth&                           target_origin,
    bool                                        ignore_hydrogens,
    const BFactorRestraintSpec&                 restraint_spec,
    const BFactorTargetRestraintSpec&           target_spec)
{
    ready_            = false;
    ignore_hydrogens_ = ignore_hydrogens;
    atom_mapping_     = mapping;
    // Resolve ChimeraX-indexed restraint specs / refine_occ flags to the
    // Clipper Atom_list indexing the refiner uses.
    RefineConfig eff = build_effective_config_(mapping, (int)refined_atoms.size(),
                                               restraint_spec, target_spec);
    refiner_.reset(new BFactorOccRefiner(refined_atoms, context_atoms,
                                         target_map, target_origin, eff));
    thread_result_ = std::async(std::launch::async,
                                &BFactorOccRefinerThread::run_, this);
}

RefineConfig BFactorOccRefinerThread::build_effective_config_(
    const std::vector<bridge::AtomAltlocIndex>& mapping,
    int                                         n_clipper,
    const BFactorRestraintSpec&                 rspec,
    const BFactorTargetRestraintSpec&           tspec)
{
    RefineConfig eff = cfg_;

    // (cx_index, altloc) -> Clipper Atom_list index.
    std::map<std::pair<int, char>, int> resolver;
    for (const auto& e : mapping)
        resolver[std::make_pair(e.cx_index, e.altloc)] = e.index;

    auto altloc_char = [](const std::string& s) -> char {
        return s.empty() ? '\0' : s[0];
    };
    auto resolve = [&](int cx_idx, const std::string& al) -> int {
        auto it = resolver.find(std::make_pair(cx_idx, altloc_char(al)));
        if (it == resolver.end())
            throw std::out_of_range(
                "BFactorOccRefiner: restraint references ChimeraX atom index "
                + std::to_string(cx_idx) + " altloc '" + al + "', which has no "
                "corresponding refined Clipper atom (excluded by ignore_hydrogens, "
                "or a non-existent altloc).");
        return it->second;
    };

    // Expand the ChimeraX-indexed refine_occ flags to one entry per Clipper atom.
    const bool occ_requested = !cfg_.refine_occ.empty();
    if (occ_requested) {
        std::vector<uint8_t> expanded((size_t)n_clipper, 0);
        for (const auto& e : mapping)
            if (e.cx_index >= 0 && e.cx_index < (int)cfg_.refine_occ.size())
                expanded[e.index] = cfg_.refine_occ[e.cx_index];
        eff.refine_occ = std::move(expanded);
    }

    // Resolve pairwise restraints to Clipper indices.
    eff.b_restraints.clear();
    for (size_t k = 0; k < rspec.atoms1.size(); ++k) {
        BFactorRestraint r;
        r.i      = resolve(rspec.atoms1[k], rspec.altlocs1[k]);
        r.j      = resolve(rspec.atoms2[k], rspec.altlocs2[k]);
        r.weight = rspec.weights[k];
        r.sigma  = rspec.sigmas[k];
        r.alpha  = (k < rspec.alphas.size()) ? rspec.alphas[k] : 1.0;
        eff.b_restraints.push_back(r);
    }

    // Resolve one-sided target restraints to Clipper indices.
    eff.b_target_restraints.clear();
    for (size_t k = 0; k < tspec.atoms.size(); ++k) {
        BFactorTargetRestraint r;
        r.i        = resolve(tspec.atoms[k], tspec.altlocs[k]);
        r.target_u = tspec.target_us[k];
        r.weight   = tspec.weights[k];
        r.sigma    = tspec.sigmas[k];
        r.alpha    = (k < tspec.alphas.size()) ? tspec.alphas[k] : 1.0;
        eff.b_target_restraints.push_back(r);
    }

    // Auto-derive occupancy groups: atoms in one contiguous covalent fragment
    // that share an altloc get one shared occupancy (EqualOccGroup); the lettered
    // altlocs of a fragment sum to 1 (OccConstraintGroup).  Gated by refine_occ.
    eff.equal_occ_groups.clear();
    eff.occ_groups.clear();
    if (occ_requested) {
        // Flagged ChimeraX atoms -> their (Clipper index, altloc) entries.
        std::unordered_map<atomstruct::Atom*, std::vector<std::pair<int, char>>> flagged;
        std::unordered_set<atomstruct::Atom*> flagged_set;
        for (const auto& e : mapping) {
            if (e.cx_index < 0 || e.cx_index >= (int)cfg_.refine_occ.size()) continue;
            if (!cfg_.refine_occ[e.cx_index]) continue;
            flagged[e.cx_atom].push_back(std::make_pair(e.index, e.altloc));
            flagged_set.insert(e.cx_atom);
        }

        // Covalent fragments = connected components over bonds among flagged atoms.
        std::unordered_set<atomstruct::Atom*> visited;
        for (const auto& kv : flagged) {
            atomstruct::Atom* seed = kv.first;
            if (visited.count(seed)) continue;
            std::vector<atomstruct::Atom*> stack{seed}, component;
            visited.insert(seed);
            while (!stack.empty()) {
                atomstruct::Atom* a = stack.back(); stack.pop_back();
                component.push_back(a);
                for (auto* bond : a->bonds()) {
                    atomstruct::Atom* nb = bond->other_atom(a);
                    if (flagged_set.count(nb) && !visited.count(nb)) {
                        visited.insert(nb);
                        stack.push_back(nb);
                    }
                }
            }

            // Gather this fragment's Clipper indices grouped by altloc.
            std::map<char, std::vector<int>> by_altloc;
            for (atomstruct::Atom* a : component)
                for (const auto& ci_al : flagged[a])
                    by_altloc[ci_al.second].push_back(ci_al.first);

            std::vector<char> letters;
            for (const auto& kv2 : by_altloc)
                if (kv2.first != '\0') letters.push_back(kv2.first);

            if (letters.size() >= 2) {
                // Multi-conformer fragment: one EqualOccGroup per altloc, summing
                // to 1.  Shared ('\0') atoms are held at full occupancy (flag
                // cleared) and excluded from the sum.
                OccConstraintGroup og; og.total = 1.0;
                for (char c : letters) {
                    EqualOccGroup eg; eg.atom_indices = by_altloc[c];
                    eff.equal_occ_groups.push_back(eg);
                    og.atom_indices.push_back(eg.atom_indices.front());
                }
                eff.occ_groups.push_back(og);
                auto it0 = by_altloc.find('\0');
                if (it0 != by_altloc.end())
                    for (int ci : it0->second) eff.refine_occ[ci] = 0;
            } else {
                // Single altloc (a letter or '\0'): tie the fragment's atoms to one
                // shared free occupancy.  A lone atom needs no group — init_derived
                // treats an ungrouped flagged atom as a free occupancy parameter.
                std::vector<int> all;
                for (const auto& kv2 : by_altloc)
                    for (int ci : kv2.second) all.push_back(ci);
                if (all.size() >= 2) {
                    EqualOccGroup eg; eg.atom_indices = std::move(all);
                    eff.equal_occ_groups.push_back(eg);
                }
            }
        }
    }

    return eff;
}

std::pair<double, double> BFactorOccRefiner::compute_rfactors()
{
    if (realspace_mode_) return {-1.0, -1.0};

    // Build atom list with current refined parameters — same as operator().
    Atom_list atomu = atoms_;
    for (int j = 0; j < n_atoms_; ++j) {
        if (atomu[j].is_null()) continue;
        if (atomu[j].u_aniso_orth().is_null())
            atomu[j].set_u_iso((ftype)current_u_iso_[j]);
        atomu[j].set_occupancy((ftype)current_occ_[j]);
    }

    // EDcalc + FFT to obtain Fcalc (atomic only).
    density_xmap_ = ftype32(0);
    edcalc_(density_xmap_, atomu);
    density_xmap_.fft_to(fcalc_hkl_, cfg_.n_threads);

    // Add bulk-solvent contribution for a realistic R-factor estimate.
    if (has_bulk_solvent_) {
        for (HKL_info::HKL_reference_index ih = fobs_.first(); !ih.last(); ih.next()) {
            if (fcalc_hkl_[ih].missing() || fbulk_hkl_[ih].missing()) continue;
            std::complex<ftype32> f_total =
                (std::complex<ftype32>)fcalc_hkl_[ih]
              + (std::complex<ftype32>)fbulk_hkl_[ih];
            fcalc_hkl_.set_data(ih.hkl(), F_phi<ftype32>(f_total));
        }
    }

    // Isotropic scale k = Σ(|Fo||Fc|) / Σ|Fc|²  (working set only, flag != 0):
    // the least-squares scale of |Fc| onto |Fo|.  STANDARD (not FOM-weighted) so the
    // reported R-factors are the same quantity the manager/GUI reports and the
    // before/after bail-out comparison is like-for-like.  (FOM-weighting here used to
    // inflate the "R" by ~(1−⟨m⟩), spuriously tripping the bail-out on incomplete
    // models where ⟨m⟩ ≪ 1.)
    ftype sum_cross = 0.0, sum_sq = 0.0;
    for (HKL_info::HKL_reference_index ih = fobs_.first(); !ih.last(); ih.next()) {
        if (fobs_[ih].missing() || usage_[ih].missing()
            || usage_[ih].flag() == 0) continue;
        ftype fc = fcalc_hkl_[ih].f();
        ftype fo = fobs_[ih].f();
        sum_cross += fo * fc;
        sum_sq    += fc * fc;
    }
    ftype k = (sum_sq > ftype(1e-10)) ? sum_cross / sum_sq : ftype(1.0);

    // Standard R-work / R-free: Σ|k·|Fc| − |Fo|| / Σ|Fo|.
    // Convention: flag != 0 → working; flag == 0 → free.
    double rw_num = 0.0, rw_den = 0.0;
    double rf_num = 0.0, rf_den = 0.0;
    for (HKL_info::HKL_reference_index ih = fobs_.first(); !ih.last(); ih.next()) {
        if (fobs_[ih].missing() || usage_[ih].missing()) continue;
        ftype fc  = k * fcalc_hkl_[ih].f();
        ftype fo  = fobs_[ih].f();
        ftype res = std::abs(fc - fo);
        if (usage_[ih].flag() != 0) { rw_num += res; rw_den += fo; }
        else                         { rf_num += res; rf_den += fo; }
    }
    double rwork = (rw_den > 1e-10) ? rw_num / rw_den : -1.0;
    double rfree = (rf_den > 1e-10) ? rf_num / rf_den : -1.0;
    return {rwork, rfree};
}

std::pair<double, double> BFactorOccRefinerThread::compute_rfactors()
{
    if (thread_result_.valid()) thread_result_.get();   // block if running
    if (!refiner_) return {-1.0, -1.0};
    return refiner_->compute_rfactors();
}

std::pair<double, double> BFactorOccRefinerThread::initial_rfactors()
{
    if (thread_result_.valid()) thread_result_.get();   // block if running
    if (!refiner_) return {-1.0, -1.0};
    return refiner_->initial_rfactors();
}

Xmap<ftype32>* BFactorOccRefinerThread::compute_realspace_density(
    const Atom_list&                            refined_atoms,
    const std::vector<bridge::AtomAltlocIndex>& /*mapping*/,
    const Atom_list&                            context_atoms,
    const Xmap<ftype32>&                        target_map,
    const Coord_orth&                           target_origin)
{
    // Heap-allocate the result Xmap so ownership can be transferred to Python.
    Xmap<ftype32>* xmap = new Xmap<ftype32>();
    xmap->init(target_map.spacegroup(), target_map.cell(),
               target_map.grid_sampling());
    *xmap = ftype32(0);

    // Shift all atoms to the P1-cell frame — same transform as operator_realspace_().
    Atom_list atomu_all;
    auto shift_atom = [&](const Atom& a, ftype u_iso, ftype occ) {
        Atom s = a;
        s.set_coord_orth(a.coord_orth() - target_origin);
        if (s.u_aniso_orth().is_null())
            s.set_u_iso(u_iso);
        s.set_occupancy(occ);
        return s;
    };
    for (int c = 0; c < (int)context_atoms.size(); ++c) {
        const Atom& ca = context_atoms[c];
        if (!ca.is_null())
            atomu_all.push_back(shift_atom(ca, ca.u_iso(), ca.occupancy()));
    }
    for (int j = 0; j < (int)refined_atoms.size(); ++j) {
        const Atom& a = refined_atoms[j];
        if (!a.is_null())
            atomu_all.push_back(shift_atom(a, a.u_iso(), a.occupancy()));
    }

    EDcalc_aniso_thread<ftype32> edcalc((size_t)std::max(1, cfg_.n_threads));
    edcalc(*xmap, atomu_all);
    return xmap;   // ownership transferred to Python via take_ownership
}

void BFactorOccRefinerThread::cancel()
{
    if (refiner_)
        refiner_->cancel();
}

bool BFactorOccRefinerThread::run_()
{
    bool ok = refiner_->refine();
    ready_ = true;
    return ok;
}

int BFactorOccRefinerThread::n_atoms() const
{
    if (!refiner_) return 0;
    return refiner_->n_atoms();
}

const std::vector<double>& BFactorOccRefinerThread::refined_u_iso()
{
    if (thread_result_.valid())
        thread_result_.get();
    return refiner_->refined_u_iso();
}

const std::vector<double>& BFactorOccRefinerThread::refined_occ()
{
    if (thread_result_.valid())
        thread_result_.get();
    return refiner_->refined_occ();
}

} // namespace clipper_cx
