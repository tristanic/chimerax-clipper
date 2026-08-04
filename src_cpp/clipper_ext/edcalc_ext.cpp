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

 #include <future>
 #include <chrono>
 #include <algorithm>
 #include <vector>

#include "edcalc_ext.h"
#include <clipper/core/atomsf.h>

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{

template <class T>
bool EDcalc_mask_vdw<T>::operator() (Xmap<T>& xmap, const Atom_list& atoms ) const
{
    bool result;
    if (n_threads_ > 1)
        result = edcalc_threaded_(xmap, atoms);
    else
        result = edcalc_single_thread_(xmap, atoms);
    return result;
}


template <class T>
bool EDcalc_mask_vdw<T>::edcalc_single_thread_ (Xmap<T>& xmap, const Atom_list& atoms ) const
{
    T zero = 0.0;
    xmap = zero; // Set all values to zero
    const Cell& cell = xmap.cell();
    const Grid_sampling& grid = xmap.grid_sampling();

    Coord_orth xyz;
    Coord_grid g0, g1;
    Grid_range gd(cell, grid, grid_radius_);
    ftype atom_radius;
    typename Xmap<T>::Map_reference_coord i0, iu, iv, iw;
    // Partial-occupancy atoms (0 < occ < 1) are deferred to a second pass that
    // accumulates fractional protein content; full-occupancy atoms stamp a solid
    // sphere here as before.  With occupancy weighting off, every non-zero atom is
    // treated as full and partial_atoms stays empty.
    std::vector<size_t> partial_atoms;
    // Step 1: set all grid points within (atom radius + probe radius) to 1.
    for (size_t ai=0; ai<atoms.size(); ++ai)
    {
        const Atom& a = atoms[ai];
        if (a.is_null()) continue;
        if (ignore_zero_occ_atoms_ && a.occupancy() == 0) continue;
        if (occupancy_weighted_ && a.occupancy() < 1.0) { partial_atoms.push_back(ai); continue; }
        try {
            atom_radius = clipper_cx::data::vdw_radius(a.element());
        } catch (...) {
            std::stringstream str;
            str << "Could not find van der Waals radius for atom name " << a.element();
            throw std::runtime_error(str.str());
        }
        xyz = a.coord_orth();
        g0 = xmap.coord_map( xyz ).coord_grid() + gd.min();
        g1 = xmap.coord_map( xyz ).coord_grid() + gd.max();
        i0 = typename Xmap<T>::Map_reference_coord( xmap, g0 );
        for (iu = i0; iu.coord().u() < g1.u(); iu.next_u() ) {
            for (iv=iu; iv.coord().v() < g1.v(); iv.next_v() ) {
                for (iw=iv; iw.coord().w() < g1.w(); iw.next_w() ) {
                    if ( (xyz-iw.coord_orth()).lengthsq() < pow(atom_radius+probe_radius_, 2) )
                        xmap[iw] = 1.0;
                }
            }
        }
    }
    // Step 1b: accumulate fractional protein content for partial-occupancy atoms.
    // protein_content = min(1, sum of covering occupancies).  Single-threaded, so a
    // plain read-modify-write is race-free; voxels already saturated to 1.0 by a
    // full-occupancy atom are left untouched.
    for (size_t ai: partial_atoms)
    {
        const Atom& a = atoms[ai];
        const T occ = (T)a.occupancy();
        atom_radius = clipper_cx::data::vdw_radius(a.element());
        xyz = a.coord_orth();
        g0 = xmap.coord_map( xyz ).coord_grid() + gd.min();
        g1 = xmap.coord_map( xyz ).coord_grid() + gd.max();
        i0 = typename Xmap<T>::Map_reference_coord( xmap, g0 );
        for (iu = i0; iu.coord().u() < g1.u(); iu.next_u() ) {
            for (iv=iu; iv.coord().v() < g1.v(); iv.next_v() ) {
                for (iw=iv; iw.coord().w() < g1.w(); iw.next_w() ) {
                    if ( xmap[iw] < 1.0 &&
                         (xyz-iw.coord_orth()).lengthsq() < pow(atom_radius+probe_radius_, 2) )
                        xmap[iw] = std::min(T(1.0), xmap[iw] + occ);
                }
            }
        }
    }
    // Step 2: "shrink-wrap" the solvent mask back to the model atoms. For every
    // grid point with a non-zero value, check if it is within shrink_radius_ of
    // a grid point with a value of zero. If so, set its value to zero.
    gd = Grid_range(cell, grid, shrink_radius_);
    Xmap<T> unshrunk(xmap);
    typename Xmap<T>::Map_reference_index ix;
    for (ix = unshrunk.first(); !ix.last(); ix.next()) {
        [&] {
            if (unshrunk[ix] != zero)
            {
                xyz = ix.coord_orth();
                g0 = xmap.coord_map( xyz ).coord_grid() + gd.min();
                g1 = xmap.coord_map( xyz ).coord_grid() + gd.max();
                i0 = typename Xmap<T>::Map_reference_coord( xmap, g0 );
                for (iu = i0; iu.coord().u() < g1.u(); iu.next_u() ) {
                    for (iv=iu; iv.coord().v() < g1.v(); iv.next_v() ) {
                        for (iw=iv; iw.coord().w() < g1.w(); iw.next_w() ) {
                            if (unshrunk[iw] == zero) {
                                if ( (xyz-iw.coord_orth()).lengthsq() < pow(shrink_radius_, 2) ) {
                                    //std::cerr << "Filling in value at " << iw.coord().format() << std::endl;
                                    xmap[ix] = 0.0;
                                    return;
                                }
                            }
                        }
                    }
                }

            }
        } ();
    }
    return true;
}



template <class T>
bool EDcalc_mask_vdw<T>::edcalc_threaded_ (Xmap<T>& xmap, const Atom_list& atoms ) const
{
    // Note: thread-safety interlocks are *not* required here, since we never
    // read from and write to the same map within any thread.
    T zero = 0.0;
    xmap = zero; // Set all values to zero
    const Cell& cell = xmap.cell();
    const Grid_sampling& grid = xmap.grid_sampling();

    // Partial-occupancy atoms (0 < occ < 1) are deferred to Step 1b; full-occupancy
    // atoms stamp a solid sphere in Step 1 as before.  With occupancy weighting off
    // the list stays empty and Step 1b is skipped entirely (byte-identical, same
    // speed as the original binary mask).
    std::vector<size_t> partial_atoms;
    if (occupancy_weighted_)
        for (size_t ai=0; ai<atoms.size(); ++ai)
        {
            const auto& a = atoms[ai];
            if (a.is_null()) continue;
            if (ignore_zero_occ_atoms_ && a.occupancy() == 0) continue;
            if (a.occupancy() < 1.0) partial_atoms.push_back(ai);
        }

    // Step 1: set all grid points within (atom radius + probe radius) to 1.
    std::vector<std::future<void>> thread_results;
    size_t atoms_per_thread = atoms.size() / n_threads_ + 1;
    size_t start=0, end;
    Grid_range gd(cell, grid, grid_radius_);
    bool izo = ignore_zero_occ_atoms_;
    bool ow = occupancy_weighted_;
    for (size_t i=0; i< n_threads_; ++i)
    {
        end = std::min(start+atoms_per_thread, atoms.size());
        thread_results.push_back(std::async(std::launch::async,
            [gd, izo, ow](Xmap<T>& xmap, const Atom_list& atoms, const ftype& probe_radius, size_t start, size_t end)
            {
                Coord_orth xyz;
                Coord_grid cg, g0, g1;
                typename Xmap<T>::Map_reference_coord i0, iu, iv, iw;
                for (size_t j=start; j<end; ++j)
                {
                    const auto& a = atoms[j];
                    if (a.is_null()) continue;
                    if (izo && a.occupancy() == 0) continue;
                    if (ow && a.occupancy() < 1.0) continue;  // deferred to Step 1b
                    const auto atom_radius = clipper_cx::data::vdw_radius(a.element());
                    xyz = a.coord_orth();
                    cg = xmap.coord_map(xyz).coord_grid();
                    g0 = cg + gd.min();
                    g1 = cg + gd.max();
                    i0 = typename Xmap<T>::Map_reference_coord( xmap, g0 );
                    for (iu=i0; iu.coord().u() < g1.u(); iu.next_u() ) {
                        for (iv=iu; iv.coord().v() < g1.v(); iv.next_v() ) {
                            for (iw=iv; iw.coord().w() < g1.w(); iw.next_w() ) {
                                if ( (xyz-iw.coord_orth()).lengthsq() < pow(atom_radius+probe_radius, 2) )
                                    xmap[iw] = 1.0;
                            }
                        }
                    }
                }
            },
            std::ref(xmap), std::cref(atoms), probe_radius_, start, end
        ));
        start += atoms_per_thread;
    }
    for (auto& r: thread_results)
        r.get();

    // Step 1b: accumulate fractional protein content for partial-occupancy atoms,
    // protein_content = min(1, sum of covering occupancies).  Uses Xmap's per-grid-
    // point lock (test_and_set) with the retry-later pattern from EDcalc_aniso_thread:
    // on a lock collision the voxel is deferred to a later sweep rather than spin-
    // waiting (fast because collisions are rare).  Safe because Step 1 has fully
    // joined -- no lock-free writer is live -- and every writer here takes the lock;
    // the capped-additive rule is order-independent, so the result is deterministic.
    if (!partial_atoms.empty())
    {
        thread_results.clear();
        const size_t n_partial = partial_atoms.size();
        size_t partial_per_thread = n_partial / n_threads_ + 1;
        start = 0;
        for (size_t i=0; i<n_threads_ && start < n_partial; ++i)
        {
            end = std::min(start+partial_per_thread, n_partial);
            thread_results.push_back(std::async(std::launch::async,
                [gd, &partial_atoms](Xmap<T>& xmap, const Atom_list& atoms, const ftype& probe_radius, size_t start, size_t end)
                {
                    Coord_orth xyz;
                    Coord_grid cg, g0, g1;
                    typename Xmap<T>::Map_reference_coord i0, iu, iv, iw;
                    std::vector<typename Xmap<T>::Map_reference_coord> in_range, remaining;
                    for (size_t k=start; k<end; ++k)
                    {
                        const auto& a = atoms[partial_atoms[k]];
                        const T occ = (T)a.occupancy();
                        const auto atom_radius = clipper_cx::data::vdw_radius(a.element());
                        const auto r2 = pow(atom_radius+probe_radius, 2);
                        xyz = a.coord_orth();
                        cg = xmap.coord_map(xyz).coord_grid();
                        g0 = cg + gd.min();
                        g1 = cg + gd.max();
                        i0 = typename Xmap<T>::Map_reference_coord( xmap, g0 );
                        in_range.clear();
                        remaining.clear();
                        for (iu=i0; iu.coord().u() < g1.u(); iu.next_u() )
                            for (iv=iu; iv.coord().v() < g1.v(); iv.next_v() )
                                for (iw=iv; iw.coord().w() < g1.w(); iw.next_w() )
                                    if ( (xyz-iw.coord_orth()).lengthsq() < r2 )
                                        in_range.push_back(iw);
                        while (!in_range.empty())
                        {
                            for (const auto& ix: in_range)
                            {
                                if (xmap.lock_element(ix)) { remaining.push_back(ix); continue; }
                                if (xmap[ix] < 1.0)
                                    xmap[ix] = std::min(T(1.0), xmap[ix] + occ);
                                xmap.unlock_element(ix);
                            }
                            in_range.swap(remaining);
                            remaining.clear();
                        }
                    }
                },
                std::ref(xmap), std::cref(atoms), probe_radius_, start, end
            ));
            start += partial_per_thread;
        }
        for (auto& r: thread_results)
            r.get();
    }

    // Step 2: "shrink-wrap" the solvent mask back to the model atoms. For every
    // grid point with a non-zero value, check if it is within shrink_radius_ of
    // a grid point with a value of zero. If so, set its value to zero.

    thread_results.clear();
    gd = Grid_range(cell, grid, shrink_radius_);
    Xmap<T> unshrunk(xmap);

    size_t points_per_thread = (size_t)unshrunk.unique_points() / n_threads_ + 1;
    start = 0;
    for (size_t i=0; i<n_threads_; ++i)
    {
        thread_results.push_back(std::async(std::launch::async,
            [start, points_per_thread, gd, zero](const Xmap<T>& unshrunk, Xmap<T>& xmap, const ftype& shrink_radius)
            {
                auto ix = typename Xmap<T>::Map_reference_index(unshrunk, (int)start);
                typename Xmap<T>::Map_reference_coord iu, iv, iw;
                Coord_orth xyz;
                Coord_grid cg, g0, g1;
                typename Xmap<T>::Map_reference_coord i0;
                for (size_t j=0; j<points_per_thread && !ix.last(); ix.next(), ++j)
                [&]{
                    if (unshrunk[ix] == zero)
                        return;
                    xyz = ix.coord_orth();
                    cg = xmap.coord_map(xyz).coord_grid();
                    g0 = cg+gd.min();
                    g1 = cg+gd.max();
                    i0 = typename Xmap<T>::Map_reference_coord(xmap, g0);
                    for (iu=i0; iu.coord().u() < g1.u(); iu.next_u() ) {
                        for (iv=iu; iv.coord().v() < g1.v(); iv.next_v() ) {
                            for (iw=iv; iw.coord().w() < g1.w(); iw.next_w() ) {
                                if (unshrunk[iw] == zero) {
                                    if ((xyz-iw.coord_orth()).lengthsq() < pow(shrink_radius, 2) ) {
                                        xmap[ix] = zero;
                                        return;
                                    }
                                }
                            }
                        }
                    }
                }();
            },
            std::cref(unshrunk), std::ref(xmap), shrink_radius_
        ));
        start += points_per_thread;
    }
    for (auto& r: thread_results)
        r.get();

    return true;
}


template<class T>
bool EDcalc_mask_vdw<T>::operator() (NXmap<T>& nxmap, const Atom_list& atoms) const
{
    throw std::logic_error("Not implemented");
    return false;
}

template<class T>
bool EDcalc_aniso_thread<T>::operator() (Xmap<T>& xmap, const Atom_list& atoms) const
{
    xmap = 0.0;
    std::vector<std::future<bool>> results(n_threads_);
    size_t atoms_per_thread = atoms.size() / n_threads_;
    size_t start=0, end;
    for (size_t i=0; i<n_threads_; ++i)
    {
        if (i<n_threads_-1)
            end = start+atoms_per_thread;
        else
            end = atoms.size();
        results[i] = std::async(std::launch::async,
            &EDcalc_aniso_thread<T>::edcalc_xmap_thread_, this,
            std::ref(xmap), std::cref(atoms), start, end);
        start += atoms_per_thread;
    }
    // Wait for all threads to finish
    for (auto&r: results)
        r.get();

    for (const auto& spos: xmap.special_positions()) {
        xmap.set_data(spos.first, spos.second*xmap.get_data(spos.first));
    }
    return true;
}

template <class T>
bool EDcalc_aniso_thread<T>::edcalc_xmap_thread_(Xmap<T>& xmap,
    const Atom_list& atoms, size_t start, size_t end) const
{
    const Cell& cell = xmap.cell();
    const Grid_sampling& grid = xmap.grid_sampling();

    Coord_orth xyz;
    Coord_grid g0, g1;
    Grid_range gd; // ( cell, grid, radius_ );
    typename Xmap<T>::Map_reference_coord i0, iu, iv, iw;
    std::vector<typename Xmap<T>::Map_reference_coord> grid_ref_coords, remaining_coords;
    for (size_t i=start; i<end; ++i)
    {
        const Atom& a = atoms[i];
        if (a.is_null())
            continue;
        gd = Grid_range(cell, grid, cutoff_radius(a));
        grid_ref_coords.clear();
        remaining_coords.clear();
        U_aniso_orth u (a.u_aniso_orth());
        if ( u.is_null() ) u = U_aniso_orth(a.u_iso());
        AtomShapeFn sf (a.coord_orth(), a.element(), u, a.occupancy(), radiation_);
        auto cg = xmap.coord_map(a.coord_orth()).coord_grid();
        g0 = cg+gd.min();
        g1 = cg+gd.max();
        i0 = typename Xmap<T>::Map_reference_coord(xmap, g0);
        for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u())
            for (iv=iu; iv.coord().v() <= g1.v(); iv.next_v())
                for (iw=iv; iw.coord().w() <= g1.w(); iw.next_w())
                    remaining_coords.push_back(iw);
        while (!remaining_coords.empty())
        {
            grid_ref_coords.swap(remaining_coords);
            for (const auto& ix: grid_ref_coords)
            {
                // xmap.lock_element(ix) returns true if ix is locked by another thread.
                if (xmap.lock_element(ix))
                {
                    remaining_coords.push_back(ix);
                    continue;
                }
                xmap[ix] += sf.rho( ix.coord_orth());
                xmap.unlock_element(ix);
            }
            grid_ref_coords.clear();
        }
    }
    return true;
}


template class CLIPPER_CX_IMEX EDcalc_mask_vdw<ftype32>;
template class CLIPPER_CX_IMEX EDcalc_mask_vdw<ftype64>;

template class CLIPPER_CX_IMEX EDcalc_aniso_thread<ftype32>;
template class CLIPPER_CX_IMEX EDcalc_aniso_thread<ftype64>;


} // namespace clipper_cx
