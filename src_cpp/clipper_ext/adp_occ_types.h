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
#include <cstdint>
#include <vector>
#include <string>

namespace clipper_cx {

//! A group of atoms that are all constrained to carry the same occupancy.
/*!
  Typically used to link all atoms of one alternate conformer: every member
  is covalently bonded to its counterpart in another conformer, so physically
  they must have identical occupancy.  The group contributes exactly one free
  parameter (the shared value); all other members are derived from it.
  The first atom in atom_indices is taken as the representative; its index is
  what OccConstraintGroup.atom_indices should reference when combining the
  two constraint types.
*/
struct EqualOccGroup {
    std::vector<int> atom_indices;  //!< Indices into the Atom_list being refined
};

//! A group of atoms (or EqualOccGroup representatives) whose occupancies sum to a given total.
/*!
  Represents e.g. the sum constraint occ_A + occ_B = 1.0 across two altlocs.
  When used together with EqualOccGroup, atom_indices should contain exactly
  one representative per EqualOccGroup (i.e. atom_indices[k] == the first
  element of the corresponding EqualOccGroup), so that each entry refers to
  one unique free parameter.
  The last entry in atom_indices is derived as total - sum(others).
*/
struct OccConstraintGroup {
    std::vector<int> atom_indices;  //!< Indices into the Atom_list being refined
    double total = 1.0;             //!< Occupancy sum constraint
};

//! Robust restraint between the isotropic U values of two atoms, using Barron's
//! general loss (J.T. Barron, "A General and Adaptive Robust Loss Function",
//! CVPR 2019 / arXiv 2017).
/*!
  With r = U_i - U_j (U = B / (8π²), in Å²), scale c = `sigma`, and shape α = `alpha`:

    E(r) = weight · (|α−2|/α) · ( ( (r/c)²/|α−2| + 1 )^(α/2) − 1 )      (general)

  with the removable-singularity limits handled separately:
    α = 2   : E = weight · ½(r/c)²                      (squared error / harmonic)
    α = 0   : E = weight · log(½(r/c)² + 1)             (Cauchy / Lorentzian)
    α → −∞  : E = weight · (1 − exp(−½(r/c)²))          (Welsch / Leclerc)

  The shape α smoothly interpolates the classic robust losses: α=2 harmonic,
  α=1 Charbonnier/pseudo-Huber (quadratic core, *linear* tail — bounded but
  non-vanishing influence), α=0 Cauchy, α=−2 Geman-McClure (redescending tail).
  For all finite α the loss is ≈ weight·½(r/c)² near r=0, so curvature at the
  origin is weight/c² regardless of α: c is a pure stiffness knob and α a pure
  tail-shape knob. α ∈ [0, 1] keeps a persistent restoring force on large
  differences (influence saturates rather than redescending to zero), which is
  generally better behaved than Geman-McClure at low resolution.

  weight : dimensionless overall scale (origin curvature = weight/c²).
  sigma  : scale c in U-space (Å²), the width of the quadratic core.
           Typical choice: sigma ≈ Util::b2u(target_delta_B), e.g.
           sigma = Util::b2u(5.0) ≈ 0.063 Å² for a 5 Å² B-factor tolerance.
  alpha  : shape parameter (default 1.0 = Charbonnier; −2 = Geman-McClure).
*/
//! One-sided Barron restraint: atom i toward a fixed target U value.
/*! Used to harmonise RSR-refined atoms with nearby context atoms.
    Energy / gradient are identical to BFactorRestraint (with r = U_i − target_u)
    but only atom i's gradient is accumulated (the target is fixed, not a refined
    parameter). */
struct BFactorTargetRestraint {
    int    i        = 0;
    double target_u = 0.0;   //!< Fixed target U_iso (Å²) from context atom
    double weight   = 1.0;
    double sigma    = 0.063; //!< Barron scale c (Å²), same default as BFactorRestraint
    double alpha    = 1.0;   //!< Barron shape (1 = Charbonnier, −2 = Geman-McClure)
};

struct BFactorRestraint {
    int    i      = 0;
    int    j      = 0;
    double weight = 1.0;
    double sigma  = 0.063;   //!< Barron scale c (Å², U-space); ≈ 5 Å² in B-space
    double alpha  = 1.0;     //!< Barron shape (1 = Charbonnier, −2 = Geman-McClure)
};

//! Caller-supplied pairwise B-factor restraint specification, expressed in terms
//! of the INPUT ChimeraX atom array: parallel arrays of atom-array indices plus
//! an explicit altloc per endpoint ("" = no altloc).  Resolved to Clipper-indexed
//! BFactorRestraints at launch.  The index/altloc/sigma/weight arrays have the
//! same length (= number of restraints); `alphas` may be empty (each restraint
//! then defaults to α=1) or the same length.  The client decides the altloc pairing.
struct BFactorRestraintSpec {
    std::vector<int>         atoms1;
    std::vector<std::string> altlocs1;
    std::vector<int>         atoms2;
    std::vector<std::string> altlocs2;
    std::vector<double>      sigmas;
    std::vector<double>      weights;
    std::vector<double>      alphas;   //!< Barron shape per restraint (empty = all α=1)
};

//! Caller-supplied one-sided B-factor restraint specification (toward a fixed
//! target U value).  ChimeraX-atom-indexed + altloc; resolved at launch.
//! `alphas` may be empty (each restraint then defaults to α=1) or the same length.
struct BFactorTargetRestraintSpec {
    std::vector<int>         atoms;
    std::vector<std::string> altlocs;
    std::vector<double>      target_us;
    std::vector<double>      sigmas;
    std::vector<double>      weights;
    std::vector<double>      alphas;   //!< Barron shape per restraint (empty = all α=1)
};

//! Configuration for B-factor and/or occupancy refinement.
struct RefineConfig {
    bool   refine_b        = true;
    //! Per-atom occupancy refinement flags (uint8_t: 0 = fixed, 1 = refine).
    //! Must be either empty (no occ refinement) or exactly n_atoms long.
    //! All atoms in an OccConstraintGroup must carry the same flag value;
    //! this is enforced at the Python level.
    std::vector<uint8_t> refine_occ;
    bool   use_curvature   = false;    //!< Use second derivatives for preconditioning (reserved)
    int    max_cycles      = 100;
    double b_min           = 0.5;      //!< Å², lower bound on isotropic B-factors
    double b_max           = 999.0;    //!< Å², upper bound on isotropic B-factors
    double lbfgs_epsilon   = 1e-5;     //!< Gradient norm convergence threshold
    int    lbfgs_past      = 3;        //!< Past iterations for delta-f stopping criterion
    double lbfgs_delta     = 0.0;      //!< Delta-f threshold (0 = disabled)
    int    n_threads       = 1;        //!< Thread count for EDcalc and FFT steps
    //! Use electron scattering factors (micro-ED / 3D-ED and cryo-EM potential
    //! maps) instead of X-ray for the density and Agarwal gradients. Default
    //! false (X-ray) keeps existing refinement byte-identical. Stored as a bool
    //! rather than AtomShapeFn::RADIATION to keep this header free of the clipper
    //! core include; converted at each AtomShapeFn/EDcalc construction site.
    bool   use_electron_scattering = false;

    //! Equal-occupancy groups (must be set before occ_groups).
    //! Each group contributes one free parameter. Representatives (first
    //! member of each group) should be referenced by occ_groups.
    std::vector<EqualOccGroup>      equal_occ_groups;
    std::vector<OccConstraintGroup> occ_groups;
    std::vector<BFactorRestraint>   b_restraints;
    std::vector<BFactorTargetRestraint> b_target_restraints;
};

} // namespace clipper_cx
