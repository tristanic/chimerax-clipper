// PyBind11 Python bindings for Clipper
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
#include <cmath>
#include <iostream>

#include <clipper/clipper.h>
#include "symops.h"
namespace clipper
{

class Unit_Cell
{
public:
    ~Unit_Cell() {};

    Unit_Cell(const Atom_list& atoms, const Cell& cell,
        const Spacegroup& sg, const Grid_sampling& grid, int padding=0);

    inline const RTop_fracs& symops() const { return symops_; } // Getter
    inline const RTop_fracs& inv_symops() const { return inv_symops_; } // Getter

    inline const Coord_frac& ref_coord() const { return ref_; } // Getter
    inline void set_ref_coord(const Coord_frac& new_ref) {
        ref_ = new_ref;
        recalculate_symops_();
    }

    inline const Cell& cell() const { return cell_; }
    inline const Spacegroup& spacegroup() const { return sg_; }
    inline const Grid_sampling& grid() const { return grid_; }

    inline const Grid_range& ref_box() const { return reference_model_bounds_; }

    inline const Coord_grid& min() const { return cell_origin_; }
    inline Coord_grid max() const {
        return cell_origin_ + Coord_grid(grid_.nu()-1, grid_.nv()-1, grid_.nw()-1); }

    RTop_fracs all_symops_in_box(const Grid_range& range, bool always_include_identity=false, int sample_frequency=2) const;
    RTop_fracs all_symops_in_box(const Coord_orth& origin_xyz, const Vec3<int>& box_size_uvw,
        bool always_include_identity=false, int sample_frequency=2) const;

    /* Recompute the reference-model bounding box, packing symops and per-symop
     * transformed bounding boxes from the current atomic coordinates. Cheap
     * enough (dominated by num_symops) to call whenever the model has moved -
     * see the "dirty flag" pattern on the Python side. Reuses the padding
     * supplied at construction.
     */
    void update_reference(const Atom_list& atoms);

    /* Same as update_reference(), but reads a flat row-major array of N
     * orthographic coordinates (N x 3) instead of a Clipper Atom_list. Avoids
     * constructing N heavyweight Clipper::Atom objects (element strings,
     * occupancies, ADPs) just to derive a bounding box - the dominant cost of the
     * atom-list path - so it is cheap even for very large models.
     */
    void update_reference_from_coords(const double* xyz, size_t n);

    /* Find the minimum and maximum grid coordinates of a box encompassing the atoms,
     * pad it by padding grid units either side, and save the result as a Grid_range.
     * Also refreshes ref_ (centroid) and cell_origin_. Callers that need the search
     * caches consistent should use update_reference() instead (which also rebuilds
     * the symop / bounding-box caches).
     */
    void update_reference_model_bounds(const Atom_list& atoms, const int& padding);


private:

    Cell cell_;
    Spacegroup sg_;
    Grid_sampling grid_;
    int padding_ = 0;
    Coord_frac ref_;  // Reference coordinate (model centroid), in fractional coords.
    RTop_fracs symops_; // RTops mapping our reference coordinate to equivalent positions in the same unit cell
    RTop_fracs inv_symops_; // For convenience, the inverse symops are also saved.
    Coord_grid cell_origin_; // The origin of the unit cell containing our reference model
    Grid_range reference_model_bounds_; // Smallest (padded) grid range encompassing our atomic model

    /* For the box search: for each spacegroup symop i, the axis-aligned bounding
     * box (in FRACTIONAL coordinates) of the reference model transformed by that
     * symop. Recomputed whenever the reference model changes. Because crystallographic
     * symops are integer matrices in the lattice basis, this envelope is exact for
     * axis-permuting space groups and a safe (superset) over-approximation for
     * trigonal/hexagonal - either way it never drops a needed operator.
     */
    std::vector<Coord_frac> sym_frac_min_;
    std::vector<Coord_frac> sym_frac_max_;

    // If the reference coordinate or model changes, re-derive the symop arrays and
    // the per-symop transformed bounding boxes.
    void recalculate_symops_();
    void recalculate_sym_bounds_();
}; // class Unit_Cell


// Implementations

Unit_Cell::Unit_Cell(const Atom_list& atoms, const Cell& cell,
    const Spacegroup& sg, const Grid_sampling& grid, int padding)
    : cell_(cell),  sg_(sg), grid_(grid), padding_(padding)
{
    update_reference_model_bounds(atoms, padding);
    recalculate_symops_();
}

void
Unit_Cell::update_reference(const Atom_list& atoms)
{
    update_reference_model_bounds(atoms, padding_);
    recalculate_symops_();
}

void
Unit_Cell::update_reference_from_coords(const double* xyz, size_t n)
{
    if (n == 0) return;
    Coord_orth c0(xyz[0], xyz[1], xyz[2]);
    Coord_grid ref_min = c0.coord_frac(cell_).coord_grid(grid_);
    Coord_grid ref_max = ref_min;
    Coord_grid pad(padding_, padding_, padding_);
    Coord_orth ref_orth(0,0,0);
    for (size_t a=0; a<n; ++a) {
        Coord_orth co(xyz[a*3], xyz[a*3+1], xyz[a*3+2]);
        ref_orth += co;
        Coord_grid g = co.coord_frac(cell_).coord_grid(grid_);
        for (int i=0; i<3; ++i) {
            if (g[i] < ref_min[i]) ref_min[i] = g[i];
            if (g[i] > ref_max[i]) ref_max[i] = g[i];
        }
    }
    for (int i=0; i<3; ++i)
        ref_orth[i] /= double(n);
    ref_ = ref_orth.coord_frac(cell_);
    Coord_grid grid_ref = ref_.coord_grid(grid_);
    cell_origin_ = grid_ref - grid_ref.unit(grid_);
    reference_model_bounds_ = Grid_range( ref_min-pad, ref_max+pad );
    recalculate_symops_();
}

void
Unit_Cell::recalculate_symops_()
{
    // For each symop in this spacegroup, find the whole-unit-cell translation
    // needed to bring the transformed reference coordinate back into this unit
    // cell. Kept because other code (unit-cell / symmetry-axis drawing) consumes
    // symops()/inv_symops().
    Coord_frac origin_frac = cell_origin_.coord_frac(grid_);
    symops_ = RTop_fracs();
    inv_symops_ = RTop_fracs();
    for (int i=0; i<sg_.num_symops(); ++i)
    {
        const Symop& this_symop = sg_.symop(i);
        Coord_frac tc = this_symop*ref_;
        Coord_frac cell_offset = tc - tc.lattice_copy_unit() - origin_frac;
        symops_.append(this_symop, -cell_offset);
    }
    for (size_t i=0; i<symops_.size(); ++i)
        inv_symops_.append(symops_[i].inverse());

    recalculate_sym_bounds_();
}

void
Unit_Cell::recalculate_sym_bounds_()
{
    // The 8 corners of the (padded) reference-model box, in fractional coords.
    const Coord_grid& gmin = reference_model_bounds_.min();
    const Coord_grid& gmax = reference_model_bounds_.max();
    Coord_frac corners[8];
    int n_corner = 0;
    for (int a=0; a<2; ++a) for (int b=0; b<2; ++b) for (int c=0; c<2; ++c)
    {
        Coord_grid cg(a ? gmax.u() : gmin.u(),
                      b ? gmax.v() : gmin.v(),
                      c ? gmax.w() : gmin.w());
        corners[n_corner++] = cg.coord_frac(grid_);
    }

    int n = sg_.num_symops();
    sym_frac_min_.assign(n, Coord_frac());
    sym_frac_max_.assign(n, Coord_frac());
    for (int i=0; i<n; ++i)
    {
        const Symop& s = sg_.symop(i);
        Coord_frac tmin, tmax;
        for (int k=0; k<8; ++k)
        {
            Coord_frac tc = s*corners[k];
            if (k==0) { tmin = tc; tmax = tc; }
            else for (int a=0; a<3; ++a)
            {
                if (tc[a] < tmin[a]) tmin[a] = tc[a];
                if (tc[a] > tmax[a]) tmax[a] = tc[a];
            }
        }
        sym_frac_min_[i] = tmin;
        sym_frac_max_[i] = tmax;
    }
}

void
Unit_Cell::update_reference_model_bounds(const Atom_list& atoms, const int& padding)
{
    padding_ = padding;
    Coord_grid ref_min = atoms[0].coord_orth().coord_frac(cell_).coord_grid(grid_);
    Coord_grid ref_max = ref_min;
    Coord_grid pad(padding,padding,padding);
    Coord_orth ref_orth(0,0,0);
    for (const auto& atom: atoms) {
      ref_orth += atom.coord_orth();
      Coord_grid thiscoord = atom.coord_orth().coord_frac(cell_).coord_grid(grid_);
      for (size_t i = 0; i < 3; i++) {
        if (thiscoord[i] < ref_min[i]) ref_min[i] = thiscoord[i];
        if (thiscoord[i] > ref_max[i]) ref_max[i] = thiscoord[i];
      }
    }
    for (int i=0; i<3; ++i)
        ref_orth[i] /= atoms.size();
    ref_ = ref_orth.coord_frac(cell_);
    Coord_grid grid_ref = ref_.coord_grid(grid_);
    cell_origin_ = grid_ref - grid_ref.unit(grid_);

    reference_model_bounds_ = Grid_range( ref_min-pad, ref_max+pad );
}

RTop_fracs
Unit_Cell::all_symops_in_box(const Coord_orth& origin_xyz, const Vec3<int>& box_size_uvw,
    bool always_include_identity, int sample_frequency) const
{
    Coord_grid grid_min = origin_xyz.coord_frac(cell_).coord_grid(grid_);
    Coord_grid grid_max = grid_min + Coord_grid(box_size_uvw);
    Grid_range range(grid_min, grid_max);
    return all_symops_in_box(range, always_include_identity, sample_frequency);
}

RTop_fracs
Unit_Cell::all_symops_in_box(const Grid_range& range, bool always_include_identity, int /*sample_frequency (deprecated, ignored)*/) const
{
    // Exact search: for each spacegroup symop, find every integer lattice
    // translation for which the symop-transformed reference-model bounding box
    // overlaps the query box, and emit that operator. Both boxes are axis-aligned
    // in fractional space, so per-axis overlap reduces to an integer interval on
    // the translation - no sampling, so no missed operators.
    RTop_fracs ret;

    // Grid -> fractional is a per-axis diagonal scaling, so the query box's grid
    // min/max map directly to its fractional min/max corners.
    Coord_frac qmin = range.min().coord_frac(grid_);
    Coord_frac qmax = range.max().coord_frac(grid_);

    // Downstream display code (atom_and_bond_sym_transforms_by_residue) REQUIRES
    // the identity operator to be transforms[0] - it iterates from index 1,
    // treating index 0 as the untransformed master model. So collect the
    // non-identity operators separately and put identity first.
    bool identity_overlaps = false;
    std::vector<RTop_frac> others;
    const int n = sg_.num_symops();
    for (int i=0; i<n; ++i)
    {
        const Coord_frac& smin = sym_frac_min_[i];
        const Coord_frac& smax = sym_frac_max_[i];
        int tlo[3], thi[3];
        bool empty = false;
        for (int a=0; a<3; ++a)
        {
            // (smin+t, smax+t) overlaps (qmin, qmax)  <=>  t in [qmin-smax, qmax-smin]
            tlo[a] = (int)std::ceil(qmin[a] - smax[a]);
            thi[a] = (int)std::floor(qmax[a] - smin[a]);
            if (tlo[a] > thi[a]) { empty = true; break; }
        }
        if (empty) continue;
        for (int tu=tlo[0]; tu<=thi[0]; ++tu)
          for (int tv=tlo[1]; tv<=thi[1]; ++tv)
            for (int tw=tlo[2]; tw<=thi[2]; ++tw)
            {
                if (i==0 && tu==0 && tv==0 && tw==0)
                {
                    identity_overlaps = true;
                    continue; // identity is added first, below
                }
                RTop_frac op = sg_.symop(i);
                op.trn() += Coord_frac(ftype(tu), ftype(tv), ftype(tw));
                others.push_back(op);
            }
    }

    if (identity_overlaps || always_include_identity)
        ret.append(RTop_frac::identity());
    for (const auto& op: others)
        ret.append(op);

    return ret;
}


} // namespace clipper
