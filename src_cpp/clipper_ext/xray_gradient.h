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
#include <array>
#include <vector>
#include <string>
#include <memory>

#include <clipper/clipper.h>
#include "imex.h"

namespace clipper_cx {

using namespace clipper;
using namespace clipper::datatypes;

//! Reciprocal-space least-squares target form.
/*!
  AmplitudeLS reproduces the existing B-factor refiner's target,
    T = ½ Σ_h w_h (k·|Fc| − m·|Fo|)²,
  where m_h = phi_fom.fom() (the σ(A)/ML figure-of-merit weight for
  macromolecular data) and w_h = 1/σ²(Fo).  For a small-molecule crystal
  supply a phi_fom whose fom()=1 so m drops out, leaving a pure 1/σ²-weighted
  least-squares fit against the amplitudes (the phases used in the driving
  density always come from the current Fcalc, never from phi_fom).

  IntensityLS is the SHELX-style target for small-molecule data measured as
  intensities,
    T = ½ Σ_h w_h (Io − s·|Fc|²)²,   w_h = 1/σ²(Io),
  with s the intensity scale.  When only F/σF are available the intensity and
  its error are taken as Io = |Fo|², σ(Io) ≈ 2|Fo|σ(F).
*/
enum class XrayTargetKind {
    AmplitudeLS = 0,
    IntensityLS = 1
};

//! Shared Agarwal per-atom gradient accumulator.
/*!
  This is the single gradient kernel used by both the L-BFGS B-factor/occupancy
  refiner (BFactorOccRefiner) and the synchronous XrayGradientEvaluator.  Given a
  driving density d(x) already present on \p density_xmap and the current per-atom
  parameters, it accumulates the Agarwal (1978) real-space gradient

    raw_grad[j*P + c] += −Σ_x d(x) · ∂ρ_j(x)/∂p_c

  for every atom j and every selected parameter column c (0 ≤ c < P, P =
  types.size()).  Clipper works on the minimal ASU and reconstructs the full P1
  density from these atoms internally, so the integral runs over each ASU atom's
  own box exactly once — the symmetry expansion is Clipper's job (fft_from/fft_to),
  not this kernel's, and any reciprocal per-reflection symmetry weighting belongs
  in the driving coefficients that produced d(x).  Coordinate gradients (X,Y,Z)
  flow through the same AtomShapeFn::rho_grad call as the ADP/occupancy gradients,
  so requesting them costs nothing extra.

  Per atom, only the TYPEs that are meaningful for that atom are requested from
  rho_grad — an isotropic atom yields {X,Y,Z,Uiso,Occ} and an anisotropic atom
  yields {X,Y,Z,Occ,U11..U23}.  This avoids AtomShapeFn::rho_grad's two output
  hazards: the isotropic branch leaves U11..U23 columns untouched (stale), and
  the anisotropic branch never fills the Uiso slot (uninitialised).  Any output
  column that does not apply to an atom is simply left as the caller initialised
  it (accumulate nothing) — an anisotropic atom contributes 0 to a Uiso column,
  and an isotropic atom contributes 0 to a U11..U23 column.  No positive-
  definiteness is enforced here; the caller owns the physicality of u_aniso.

  \param density_xmap  Map holding the driving density d(x) (read-only here).
  \param positions     Per-atom orthogonal coordinates (already frame-shifted by
                       the caller if required; must match density_xmap's frame).
  \param elements      Per-atom element symbols.
  \param u_iso         Per-atom isotropic U (Å²); used where is_aniso[j]==0.
  \param u_aniso       Per-atom anisotropic U (Å², order u11,u22,u33,u12,u13,u23);
                       used where is_aniso[j]!=0.
  \param is_aniso      Per-atom flag (0 isotropic, non-zero anisotropic).
  \param occ           Per-atom occupancy.
  \param radii         Per-atom real-space truncation radius (Å) for the grid sum.
  \param types         Selected parameter columns, in the caller's output order.
  \param raw_grad      Output [N*P], row-major (atom-major); pre-sized and zeroed
                       by the caller.
  \param n_threads     Worker threads over the (independent) atom range.
*/
CLIPPER_CX_IMEX void accumulate_agarwal_gradient(
    const Xmap<ftype32>&                     density_xmap,
    const std::vector<Coord_orth>&           positions,
    const std::vector<String>&               elements,
    const std::vector<double>&               u_iso,
    const std::vector<std::array<double,6>>& u_aniso,
    const std::vector<uint8_t>&              is_aniso,
    const std::vector<double>&               occ,
    const std::vector<double>&               radii,
    const std::vector<AtomShapeFn::TYPE>&    types,
    std::vector<double>&                     raw_grad,
    int                                      n_threads);

//! Synchronous value-and-gradient evaluator for end-to-end ML training.
/*!
  Seeds the fixed reference data (element list + observed structure factors, or a
  fixed target map) ONCE, then answers repeated value_and_gradient() calls that
  drive the mutable per-atom parameters from plain numeric arrays.  There is no
  L-BFGS loop and no background thread: one call does one forward + gradient
  evaluation and returns immediately.  The interface speaks only C arrays /
  Clipper objects — it has no PyTorch dependency; the thin torch autograd adapter
  lives entirely on the Python side (chimerax.clipper.diff).

  Two modes, selected by constructor:
   - reciprocal-space least-squares against observed |Fo| (or Io), and
   - real-space least-squares against a fixed target Xmap.

  The overall scale (k for amplitudes, s for intensities in reciprocal mode; the
  map scale in real-space mode) is held FIXED within a single value_and_gradient
  call so that the returned value and gradient are mutually consistent.  Pass
  refresh_scale=true to re-estimate it from the current parameters before the
  evaluation (the training loop does this every few steps); pass false to reuse
  the cached scale — which is what a finite-difference gradient check must do so
  that the scale does not vary between the ±ε evaluations.
*/
class CLIPPER_CX_IMEX XrayGradientEvaluator {
public:
    //! Reciprocal-space constructor (least-squares vs observed structure factors).
    //! \param threaded_density  When true (default) the forward density ρ_calc is
    //!                          computed with the threaded EDcalc_aniso_thread
    //!                          (near-linear scaling to ~20 cores via n_threads);
    //!                          when false it falls back to the single-threaded
    //!                          companion builder. Both share the gradient kernel's
    //!                          density model, so results are identical.
    XrayGradientEvaluator(
        const std::vector<String>&        elements,
        const HKL_data<F_sigF<ftype32>>&  fobs,
        const HKL_data<Phi_fom<ftype32>>& phi_fom,
        const HKL_data<Flag>&             usage,
        XrayTargetKind                    kind,
        const HKL_data<F_phi<ftype32>>&   f_bulk    = HKL_data<F_phi<ftype32>>(),
        int                               n_threads = 1,
        bool                              threaded_density = true);

    //! Real-space constructor (least-squares vs a fixed target density).
    //! \param target_map     Fixed target Xmap.
    //! \param target_origin  Origin of the target cell in the caller's orthogonal
    //!                       frame (Å); atom coordinates are shifted by
    //!                       −target_origin before density calculation.
    //! \param threaded_density  See the reciprocal constructor.
    XrayGradientEvaluator(
        const std::vector<String>& elements,
        const Xmap<ftype32>&       target_map,
        const Coord_orth&          target_origin,
        int                        n_threads = 1,
        bool                       threaded_density = true);

    ~XrayGradientEvaluator();

    //! One value + gradient evaluation.  All array arguments have length
    //! N = elements.size() (coords: N*3, u_aniso: N*6, out_grad: N*P where
    //! P = selected.size()).  \p out_grad is caller-allocated and receives the
    //! per-atom gradient in row-major (atom-major) order, columns in the order of
    //! \p selected.  Returns the scalar target value T.
    double value_and_gradient(
        const double*                         coords,
        const double*                         u_iso,
        const double*                         u_aniso,
        const double*                         occ,
        const uint8_t*                        is_aniso,
        const std::vector<AtomShapeFn::TYPE>& selected,
        double*                               out_grad,
        bool                                  refresh_scale);

    //! Number of atoms fixed at construction.
    int n_atoms() const;

    //! Reciprocal (structure-factor) mode only. Runs the SAME forward density ->
    //! Fcalc -> scale as value_and_gradient (optionally re-fitting the scale via the
    //! shared scale_fcalc_to_fobs), then writes, per reflection in HKL_data index
    //! order (length n_reflections()), the observed amplitude Fo and the SCALED
    //! calculated amplitude s(h)*|Fc| into the two caller-allocated arrays. Off the
    //! measured/working set (or where Fcalc is missing) both entries are NaN. This
    //! lets the caller compute an R-factor through the shared
    //! reflection_tools.compute_r_factors, using the exact scaled Fcalc the loss sees
    //! (so the R is consistent with the live small-molecule map / recomputed_r_factor).
    void fobs_scaled_fcalc(
        const double* coords, const double* u_iso, const double* u_aniso,
        const double* occ, const uint8_t* is_aniso, bool refresh_scale,
        double* out_fo, double* out_scaled_fcalc);

    //! Number of reflections in the observed list (reciprocal mode; 0 in real-space
    //! mode). Sizes the arrays passed to fobs_scaled_fcalc.
    int n_reflections() const;

private:
    struct Impl;
    std::unique_ptr<Impl> p_;
};

} // namespace clipper_cx
