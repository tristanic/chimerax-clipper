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

//! Geman-McClure robust restraint between the isotropic U values of two atoms.
/*!
  Energy:   E = weight * r² / (σ² + r²)
  Gradient: dE/dr = weight * 2r·σ² / (σ² + r²)²
  where r = U_i - U_j  (U = B / (8π²), in Å²).

  The function is quadratic near r = 0 (curvature 2·weight/σ²) and plateaus
  toward `weight` as |r| → ∞, down-weighting outlier differences robustly.
  At |r| = σ, the energy is exactly weight/2.

  weight : dimensionless overall scale.
  sigma  : scale parameter in U-space (Å²); the half-maximum point.
           Typical choice: sigma ≈ Util::b2u(target_delta_B), e.g.
           sigma = Util::b2u(5.0) ≈ 0.063 Å² for a 5 Å² B-factor tolerance.
*/
//! One-sided Geman-McClure restraint: atom i toward a fixed target U value.
/*! Used to harmonise RSR-refined atoms with nearby context atoms.
    Energy / gradient are identical to BFactorRestraint but only atom i's
    gradient is accumulated (the target is fixed, not a refined parameter). */
struct BFactorTargetRestraint {
    int    i        = 0;
    double target_u = 0.0;   //!< Fixed target U_iso (Å²) from context atom
    double weight   = 1.0;
    double sigma    = 0.063; //!< Geman-McClure scale (Å²), same default as BFactorRestraint
};

struct BFactorRestraint {
    int    i      = 0;
    int    j      = 0;
    double weight = 1.0;
    double sigma  = 0.063;   //!< Geman-McClure scale (Å², U-space); ≈ 5 Å² in B-space
};

//! Caller-supplied pairwise B-factor restraint specification, expressed in terms
//! of the INPUT ChimeraX atom array: parallel arrays of atom-array indices plus
//! an explicit altloc per endpoint ("" = no altloc).  Resolved to Clipper-indexed
//! BFactorRestraints at launch.  All six arrays have the same length (= number
//! of restraints).  The client decides the altloc pairing.
struct BFactorRestraintSpec {
    std::vector<int>         atoms1;
    std::vector<std::string> altlocs1;
    std::vector<int>         atoms2;
    std::vector<std::string> altlocs2;
    std::vector<double>      sigmas;
    std::vector<double>      weights;
};

//! Caller-supplied one-sided B-factor restraint specification (toward a fixed
//! target U value).  ChimeraX-atom-indexed + altloc; resolved at launch.
struct BFactorTargetRestraintSpec {
    std::vector<int>         atoms;
    std::vector<std::string> altlocs;
    std::vector<double>      target_us;
    std::vector<double>      sigmas;
    std::vector<double>      weights;
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

    //! Equal-occupancy groups (must be set before occ_groups).
    //! Each group contributes one free parameter. Representatives (first
    //! member of each group) should be referenced by occ_groups.
    std::vector<EqualOccGroup>          equal_occ_groups;
    std::vector<OccConstraintGroup>     occ_groups;
    std::vector<BFactorRestraint>       b_restraints;
    std::vector<BFactorTargetRestraint> b_target_restraints;
};

} // namespace clipper_cx
