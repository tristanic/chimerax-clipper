/**
 * @Author: Tristan Croll <tic20>
 * @Date:   31-Aug-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 21-May-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */

 #include <future>
 #include <chrono>

#include "edcalc_ext.h"
#include <clipper/core/atomsf.h>

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{

template <class T>
bool EDcalc_mask_vdw<T>::operator() (Xmap<T>& xmap, const Atom_list& atoms ) const
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
    // Step 1: set all grid points within (atom radius + probe radius) to 1.
    for (const Atom &a: atoms) if ( !a.is_null() )
    {
        try {
            atom_radius = clipper_cx::data::vdw_radii.at(a.element().c_str());
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
    // Step 2: "shrink-wrap" the solvent mask back to the model atoms. For every
    // grid point with a non-zero value, check if it is within shrink_radius_ of
    // a grid point with a value of zero. If so, set its value to zero.
    gd = Grid_range(cell, grid, shrink_radius_);
    Xmap<T> unshrunk(xmap);
    typename Xmap<T>::Map_reference_index ix;
    for (ix = unshrunk.first(); !ix.last(); ix.next()) {
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
                                xmap[iw] = 0.0;

                            }
                        }
                    }
                }
            }

        }
    }
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
    auto start_time = std::chrono::steady_clock::now();
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
    auto end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end_time-start_time;
    std::cout << "EDcalc with " << n_threads_ << " threads took " << elapsed.count() << " seconds." << std::endl;
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
        AtomShapeFn sf (a.coord_orth(), a.element(), u, a.occupancy());
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
