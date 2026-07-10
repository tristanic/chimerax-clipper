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

// C4251: exporting SFcalc_obs_bulk_vdw<T> (CLIPPER_CX_IMEX) with HKL_data members
// warns on MSVC; the members are correct (see CLAUDE.md gotcha #7).
#pragma warning(disable: 4251)

#include "sfcalc_obs_vdw.h"
#include "edcalc_ext.h"
#include "scaling.h"
#include "util.h"

#include <algorithm>
#include <cmath>

// Eigen + LBFGSpp for the joint 2-D bulk-solvent optimisation (same in-tree
// dependencies the ADP/occ refiner uses; sfcalc_obs_vdw.cpp compiles into the
// same clipper_cx library, so the include dirs are already present).
#include <Eigen/Core>
#include <LBFGSB.h>

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{

template<class dtype> class Compute_add_scaled_fphi
{
public:
    Compute_add_scaled_fphi(const dtype& scale) : scale_(scale) {}
    const F_phi<dtype> operator () (const HKL_info::HKL_reference_index& ih,
            const F_phi<dtype>& fphi1, const F_phi<dtype>& fphi2) const
    {
        clipper::datatypes::F_phi<dtype> fphi;
        if (!fphi1.missing() && !fphi2.missing() )
            fphi = F_phi<dtype>(std::complex<dtype>(fphi1) + scale_*std::complex<dtype>(fphi2));
        return fphi;
    }
private:
    dtype scale_;
};


// Smooth, NORMALISED least-squares residual on intensities — the differentiable
// analogue of quick_r (which sums |·| and so is non-smooth).  The scale is fitted
// in LOG space, so fh = rfn.f(ih) ~ log(Io/Ic) and the intensity scale bringing
// Fc² to Fo² is exp(fh).  The result is the dimensionless ratio
// Σ(Fc²·exp(fh) − Fo²)² / Σ(Fo²)² (an R²-like quantity, O(0.01)); the
// normalisation is essential so the L-BFGS-B gradient is well-scaled —
// un-normalised intensity residuals reach ~1e16 and blow the optimiser up.
template <class T>
double quick_lsq(const HKL_data<F_phi<T>>& fcalc, const HKL_data<F_sigF<T>>& fobs,
    const ResolutionFn& rfn)
{
    double num = 0.0, den = 0.0;
    for ( auto ih = fobs.first(); !ih.last(); ih.next() )
    {
        if ( !fobs[ih].missing() && !fcalc[ih].missing() ) {
            const double eps = ih.hkl_class().epsilon();
            const double fh  = rfn.f(ih);                   // log intensity scale
            if (!std::isfinite(fh)) continue;
            const double f1 = pow(fcalc[ih].f(), 2.0)/eps * std::exp(fh);
            const double f2 = pow(fobs[ih].f(),  2.0)/eps;
            if (!std::isfinite(f1) || !std::isfinite(f2)) continue;
            const double d  = f1 - f2;
            num += d * d;
            den += f2 * f2;
        }
    }
    if (!(den > 0.0)) return 1.0e6;     // no usable data → large finite penalty
    return num / den;
}

// L-BFGS-B objective for the bulk-solvent solve over x = [k_sol, mask_U] at a
// FIXED anisotropic-Gaussian scale (set via set_scale).  Keeping the scale fixed
// during the (k, U) optimisation removes the throw-prone nonlinear scale fit from
// the line search (the driver alternates scale-fit and (k, U)-optimise — block
// coordinate descent).  Objective = the normalised smooth LSQ (quick_lsq), with
// a finite-difference gradient; correct to first order near each iterate by the
// envelope theorem, exact at self-consistency.
template <class T>
class BulkSolventObjective
{
public:
    BulkSolventObjective(HKL_data<F_phi<T>>& fphi,
                         const HKL_data<F_phi<T>>& fphi_atom,
                         const HKL_data<F_phi<T>>& fphi_mask,
                         const HKL_data<F_sigF<T>>& fsig,
                         const HKL_info& hkls, const Cell& cell)
        : fphi_(fphi), fphi_atom_(fphi_atom), fphi_mask_(fphi_mask), fsig_(fsig),
          mask_final_(hkls, cell)
    {}

    //! Set the (already-fitted) scale held fixed during the next minimize().
    void set_scale(const ResolutionFn* rfn) { rfn_ = rfn; }

    //! Assemble F_total = F_atoms + k·exp(-ua)·F_mask into fphi_ at (k, ua).
    void assemble(double k, double ua)
    {
        mask_final_.compute(fphi_mask_,
            datatypes::Compute_scale_u_iso<datatypes::F_phi<T>>(1.0, -ua));
        fphi_.compute(fphi_atom_, mask_final_, Compute_add_scaled_fphi<T>((T)k));
    }

    double operator() (const Eigen::VectorXd& x, Eigen::VectorXd& grad)
    {
        const double k = x[0], ua = x[1];
        const double f0 = eval_(k, ua);
        const double dk  = 1.0e-3;
        const double dua = Util::b2u(1.0);
        grad[0] = (eval_(k+dk, ua) - eval_(k-dk, ua)) / (2.0*dk);
        grad[1] = (eval_(k, ua+dua) - eval_(k, ua-dua)) / (2.0*dua);
        if (!std::isfinite(f0) || !std::isfinite(grad[0]) || !std::isfinite(grad[1])) {
            grad[0] = grad[1] = 0.0;
            return 1.0e6;
        }
        return f0;
    }

private:
    double eval_(double k, double ua)
    {
        assemble(k, ua);
        return quick_lsq(fphi_, fsig_, *rfn_);
    }
    HKL_data<F_phi<T>>&        fphi_;
    const HKL_data<F_phi<T>>&  fphi_atom_;
    const HKL_data<F_phi<T>>&  fphi_mask_;
    const HKL_data<F_sigF<T>>& fsig_;
    HKL_data<F_phi<T>>         mask_final_;
    const ResolutionFn*        rfn_ = nullptr;
};

// Joint least-squares optimisation of (k_sol, mask_U), warm-started from
// (k_start, ua_start).  Alternates: (1) fit the aniso-Gaussian scale at the
// current (k, U); (2) bounded L-BFGS-B over (k, U) at that fixed scale.  Repeats
// to self-consistency.  Leaves `params` holding the final aniso-Gaussian scale.
template <class T>
void optimize_bulk_params(
        HKL_data<datatypes::F_phi<T>>& fphi,
        const HKL_data<datatypes::F_phi<T>>& fphi_atom,
        const HKL_data<datatypes::F_phi<T>>& fphi_mask,
        const HKL_data<datatypes::F_sigF<T>>& fsig,
        const HKL_info& hkls, const Cell& cell,
        std::vector<ftype>& params, ftype tolerance,
        T k_start, T ua_start, T& k_out, T& ua_out)
{
    BulkSolventObjective<T>    obj(fphi, fphi_atom, fphi_mask, fsig, hkls, cell);
    // Log-space anisotropic scale: log-Gaussian basis + log target → a LINEAR
    // least-squares fit solved in one closed-form pass by a plain ResolutionFn
    // (stable, no nonlinear convergence to throw mid line-search).  fh = rfn.f(ih)
    // ~ log(Io/Ic); quick_lsq applies exp(fh) as the intensity scale.
    BasisFn_log_aniso_gaussian basisfn;
    TargetFn_scaleLogF1F2<F_phi<T>, F_sigF<T>> targetfn(fphi, fsig);  // reads fphi (assembled below)
    if (params.size() != 7) params.assign(7, 0.0);
    (void)tolerance;   // unused: the linear ResolutionFn fit has no tolerance knob

    LBFGSpp::LBFGSBParam<double> p;
    p.epsilon        = 1e-3;
    p.max_iterations = 20;
    p.past           = 2;       // stop when the objective stops improving
    p.delta          = 1e-4;
    LBFGSpp::LBFGSBSolver<double> solver(p);

    Eigen::VectorXd x(2), lb(2), ub(2);
    lb[0] = 0.0;                 ub[0] = 2.0;                 // k_sol
    lb[1] = Util::b2u(-50.0);    ub[1] = Util::b2u(200.0);    // additional mask U
    x[0] = std::min(std::max((double)k_start,  lb[0]), ub[0]);
    x[1] = std::min(std::max((double)ua_start, lb[1]), ub[1]);

    for (int outer = 0; outer < 4; ++outer) {
        // (1) Fit the (log-space, linear) aniso scale at the current (k, U)
        //     (fphi assembled there).
        obj.assemble(x[0], x[1]);
        std::unique_ptr<ResolutionFn> rfn;
        try {
            rfn.reset(new ResolutionFn(hkls, basisfn, targetfn, params));
            params = rfn->params();
        } catch (const std::exception&) {
            if (!rfn) break;   // cannot fit a scale at all → keep current (k, U)
        }
        obj.set_scale(rfn.get());

        // (2) L-BFGS-B over (k, U) at the fixed scale.
        Eigen::VectorXd xprev = x;
        double fx = 0.0;
        try { solver.minimize(obj, x, fx, lb, ub); }
        catch (const std::exception&) { /* iteration cap — keep current x */ }

        if ((x - xprev).norm() < 1.0e-4) break;   // self-consistent
    }
    k_out  = (T)x[0];
    ua_out = (T)x[1];
}

template<class T>
bool SFcalc_obs_bulk_vdw<T>::operator() ( HKL_data<datatypes::F_phi<T> >& fphi,
    const HKL_data<datatypes::F_sigF<T> >& fsig, const Atom_list& atoms)
{
  // std::cout << "Starting bulk solvent calculation..." << std::endl << std::flush;
  double u_mask = Util::b2u( 50.0 );

  // increase the U values
  // Atom_list atomu = atoms;
  // U_aniso_orth uadd( u_atom ), u;
  // for ( int i = 0; i < atomu.size(); i++ ) if ( !atomu[i].is_null() ) {
  //   u = atomu[i].u_aniso_orth();
  //   if ( u.is_null() ) u = U_aniso_orth( atomu[i].u_iso() );
  //   atomu[i].set_u_aniso_orth( u + uadd );
  // }

  // now make the map for ed calcs
  const HKL_info&   hkls = fsig.base_hkl_info();
  const Spacegroup& spgr = hkls.spacegroup();
  const Cell&       cell = fsig.base_cell();
  HKL_data<datatypes::F_phi<T> >
    fphi_atom( hkls, cell ), fphi_mask( hkls, cell );
  const Grid_sampling grid( spgr, cell, hkls.resolution() );
  Xmap<float> xmap( spgr, cell, grid );

  // If any atom has an unrealistically low B-factor for this resolution,
  // apply a uniform additive shift ΔU to ALL atom U-values before the
  // electron density / FFT calculation.  A uniform shift preserves the
  // relative B-factor differences between atoms while ensuring the resulting
  // fphi_atom has the expected Gaussian resolution falloff.  In reciprocal
  // space this is equivalent to multiplying every Fcalc(h) by exp(-ΔU·s²),
  // a factor the BasisFn_aniso_gaussian scaler can absorb cleanly through its
  // isotropic parameter.
  //
  // Clamping individual atoms to the floor would instead destroy relative
  // B-factor information (two atoms at 0.5 and 8 Å² would both become 18 Å²)
  // and produce an irregular, per-atom correction the scaler cannot describe.
  //
  // Minimum credible U at resolution d_min: U_min = d_min² / (4π²),
  // i.e. the atom Gaussian width σ = d_min / 2π (Nyquist sampling criterion).
  const ftype d_min_calc  = hkls.resolution().limit();
  const ftype u_min_calc  = Util::b2u(2.0 * d_min_calc * d_min_calc);

  // Find the minimum isotropic U across all non-null, isotropic atoms.
  ftype u_min_model = u_min_calc; // sentinel: no shift needed unless below floor
  for (int i = 0; i < (int)atoms.size(); ++i) {
      if (atoms[i].is_null()) continue;
      if (atoms[i].u_aniso_orth().is_null())
          u_min_model = std::min(u_min_model, atoms[i].u_iso());
  }

  // Apply uniform shift only when needed; mask calculation uses original atoms.
  Atom_list atomu = atoms;
  const ftype u_delta = u_min_calc - u_min_model; // zero when no shift needed
  if (u_delta > 0.0) {
      for (int i = 0; i < (int)atomu.size(); ++i) {
          if (atomu[i].is_null()) continue;
          if (atomu[i].u_aniso_orth().is_null())
              atomu[i].set_u_iso(atomu[i].u_iso() + u_delta);
      }
  }

  // do ed calc from atomu (uniformly B-shifted for proper FFT sampling)
  EDcalc_aniso_thread<ftype32> edcalc(nthreads, radiation_);
  edcalc( xmap, atomu);

  xmap.fft_to( fphi_atom, nthreads );

  // do density calc from mask.
  //
  // The solvent mask (hence its transform fphi_mask, and the mean solvent
  // fraction bulkfrc) is a function of atom positions, elements and
  // zero-occupancy status ONLY -- invariant to B-factor and occupancy-value
  // changes.  When those mask-determining inputs are unchanged from the previous
  // call (e.g. across B-factor/occupancy refinement macrocycles, where
  // coordinates are fixed) reuse the cached transform fmask_ and skip a full
  // EDcalc_mask + FFT -- roughly half the cost of this routine.
  const uint64_t mask_sig = mask_signature(atoms);
  const bool mask_reused = mask_cache_valid_ && (mask_sig == mask_signature_);

  if ( !mask_reused )
  {
    auto emcalc = EDcalc_mask_vdw<ftype32>();
    emcalc.set_num_threads(nthreads);
    emcalc( xmap, atoms );
    for ( Xmap<ftype32>::Map_reference_index ix = xmap.first();
          !ix.last(); ix.next() )
      xmap[ix] = 1.0 - xmap[ix];

    // mean solvent fraction depends only on mask geometry -- compute it here,
    // while xmap holds the (inverted) mask, and cache it via the bulkfrc member.
    ftype64 w, s0 = 0.0, s1 = 0.0;
    for ( Xmap<ftype32>::Map_reference_index ix = xmap.first();
          !ix.last(); ix.next() ) {
      w = 1.0/ftype64( xmap.multiplicity( ix ) );
      s0 += w;
      s1 += w*xmap[ix];
    }
    bulkfrc = s1/s0;

    xmap.fft_to( fphi_mask, nthreads );
    fphi_mask.compute( fphi_mask, datatypes::Compute_scale_u_iso<datatypes::F_phi<T> >( 1.0, -u_mask ) );
    fphi_mask.set_data(HKL(0,0,0), datatypes::F_phi<T>());

    // Cache the smoothed solvent-mask transform and its fingerprint for reuse.
    fmask_ = fphi_mask;
    mask_signature_ = mask_sig;
    mask_cache_valid_ = true;
  }

  // Use the mask transform without copying it out of the persisted member:
  // reference the cached fmask_ on reuse, otherwise the freshly-computed local.
  const HKL_data<F_phi<T> >& mask = mask_reused ? fmask_ : fphi_mask;

  HKL_data<F_phi<T> > fphi_mask_final (hkls, cell);

  // set (0,0,0) term to null if it exists
  fphi_atom.set_data(HKL(0,0,0), datatypes::F_phi<T>());

  if (bulk_solvent_optimization_needed_)
  {
      HKL_info selected_hkls = select_random_reflections_in_bins(fsig, 500, 20);
      auto fphi_atom_temp = reflection_subset(fphi_atom, selected_hkls);
      auto fphi_mask_temp = reflection_subset(mask, selected_hkls);
      auto fsig_temp = reflection_subset(fsig, selected_hkls);
      auto fphi_temp = reflection_subset(fphi, selected_hkls);

      T sum_fobs=0, tolerance=0;
      for (auto ih=fsig_temp.first_data(); !ih.last(); fsig_temp.next_data(ih))
        sum_fobs += fsig_temp[ih].f();
      tolerance = tolerance_frac_*sum_fobs;

      T k_start, ua_start;
      if (!bulk_solvent_ever_optimized_) {
          // Cold start: the log-space aniso scale is a linear fit (no seeding
          // needed), so just start params at zeros and the (k_sol, U) at defaults.
          *params_ = std::vector<ftype>(7, 0.0);
          k_start = 0.35; ua_start = Util::b2u(0);
      } else {
          // Warm start (k_sol, U) from the previous solve.
          k_start = bulkscl; ua_start = bulk_u;
      }

      T k_opt, ua_opt;
      optimize_bulk_params<T>(fphi_temp, fphi_atom_temp, fphi_mask_temp, fsig_temp,
                              selected_hkls, cell, *params_, tolerance,
                              k_start, ua_start, k_opt, ua_opt);
      bulkscl = k_opt;
      bulk_u  = ua_opt;
      bulk_solvent_ever_optimized_ = true;
      std::cout << "Bulk solvent B-factor: " << Util::u2b(u_mask + bulk_u)
                << " scale: " << bulkscl << std::endl;
      bulk_solvent_optimization_needed_ = false;
  }
  fphi_mask_final.compute(mask,
      datatypes::Compute_scale_u_iso<datatypes::F_phi<T>>(1.0, -bulk_u));

  // Assemble F_total = F_atoms + k_sol·F_mask_final, capturing the bulk
  // contribution F_bulk = k_sol·F_mask_final as a first-class object.
  fbulk_ = HKL_data<F_phi<T> >(hkls, cell);
  for ( HKL_data<data32::F_phi>::HKL_reference_index ih = fphi.first();
        !ih.last(); ih.next() )
  {
    if (!fphi_mask_final[ih].missing())
        fbulk_.set_data(ih.hkl(),
            F_phi<T>(bulkscl * std::complex<T>(fphi_mask_final[ih])));
    fphi[ih] = std::complex<T>(fphi_atom[ih]) +
          bulkscl * std::complex<T>(fphi_mask_final[ih]);
  }

  // bulkfrc (mean solvent fraction) is computed alongside the mask above and
  // cached in the member, so it persists across mask-reuse calls.

  return true;
}

template <class T>
uint64_t SFcalc_obs_bulk_vdw<T>::mask_signature( const Atom_list& atoms ) const
{
  // 64-bit FNV-1a over exactly the inputs EDcalc_mask_vdw consumes: the atom
  // count, and per non-null atom its orthogonal coordinates, element name and a
  // zero-occupancy flag.  A moved atom, changed element, or occupancy crossing
  // zero changes the hash and so invalidates the cached mask.  Returned by value
  // as a plain integer (no heap allocation) so the cached member is safe across
  // this exported class's DLL boundary.
  uint64_t h = 1469598103934665603ull;            // FNV-1a offset basis
  auto mix = [&h](const void* p, size_t n) {
    const unsigned char* b = reinterpret_cast<const unsigned char*>(p);
    for ( size_t i = 0; i < n; ++i ) { h ^= b[i]; h *= 1099511628211ull; }
  };
  const size_t n = atoms.size();
  mix(&n, sizeof(n));
  for ( const auto& a : atoms )
  {
    const unsigned char null_flag = a.is_null() ? 1 : 0;
    mix(&null_flag, 1);
    if ( null_flag ) continue;
    const Coord_orth c = a.coord_orth();
    const double xyz[3] = { c.x(), c.y(), c.z() };
    mix(xyz, sizeof(xyz));
    const String& el = a.element();
    mix(el.c_str(), el.length());
    const unsigned char occ_zero = (a.occupancy() == 0.0) ? 1 : 0;
    mix(&occ_zero, 1);
  }
  return h;
}

// compile templates

template class CLIPPER_CX_IMEX SFcalc_obs_bulk_vdw<ftype32>;

template class CLIPPER_CX_IMEX SFcalc_obs_bulk_vdw<ftype64>;




} // namespace clipper_cx
