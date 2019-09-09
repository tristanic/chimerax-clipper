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

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include "imex.h"
#include "vdw.h"

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{

//! Bulk solvent mask taking into account van der Waals radii
/*! As per Jiang & Brunger (1994), J Mol Biol 243: 100-115 */
template <class T>
class EDcalc_mask_vdw: public EDcalc_base<T>
{
public:
    // In my empirical tests, these slightly smaller probe/shrink radii seem better than the published ones
    EDcalc_mask_vdw( const ftype grid_radius = 3.0,
                     const ftype probe_radius = 0.9, // 0.6,
                     const ftype shrink_radius = 1.1, // 0.7,
                     const size_t n_threads = 1)
        : grid_radius_(grid_radius), probe_radius_(probe_radius),
          shrink_radius_(shrink_radius), n_threads_(n_threads)
    {}
    bool operator() (Xmap<T>& xmap, const Atom_list& atoms) const;
    bool operator() (NXmap<T>& nxmap, const Atom_list& atoms) const;
    void set_num_threads(size_t n) { n_threads_ = n; }

private:
    bool edcalc_single_thread_(Xmap<T>& xmap, const Atom_list& atoms) const;
    bool edcalc_threaded_(Xmap<T>& xmap, const Atom_list& atoms) const;
    const ftype grid_radius_;
    const ftype probe_radius_;
    const ftype shrink_radius_;
    size_t n_threads_;
}; // class EDcalc_mask_vdw

template<class T> class EDcalc_aniso_thread : public EDcalc_base<T> {
public:
    EDcalc_aniso_thread( const size_t n_threads = 2 ) : n_threads_(n_threads) {}
    bool operator() (Xmap<T>& xmap, const Atom_list& atoms) const;
    bool operator() ( NXmap<T>& nxmap, const Atom_list& atoms) const {return false;}
private:
    size_t n_threads_;
    bool edcalc_xmap_thread_(Xmap<T>& xmap, const Atom_list& atoms, size_t start, size_t end) const;
    bool edcalc_nxmap_thread_(NXmap<T>& nxmap, const Atom_list& atoms, size_t start, size_t end) const {return false;}
    // Super-simple heuristic, giving the radius at which the density value should
    // be less than 1% of its maximum.
    T cutoff_radius( const Atom& a) const
    {
        T atom_radius = data::vdw_radii.at(a.element().c_str());
        return std::max(atom_radius * (0.4 + 1.5 * pow(a.u_iso(), 0.5)), 3.0);
    }
};

} // namespace clipper_cx
