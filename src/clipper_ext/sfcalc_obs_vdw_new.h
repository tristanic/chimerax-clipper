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

#include <unordered_map>

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

#include "imex.h"

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{

// Structure factor calculation using bulk solvent taking into account
// van der Waals radii
template<class T>
class SFcalc_obs_bulk_vdw : public SFcalc_obs_base<T>
{
public:
    // default constructor
    SFcalc_obs_bulk_vdw(const int n_params=12, const size_t n_threads=8) : nparams(n_params), nthreads(n_threads) {}
    // constructor: shorthand for constructor+operator
    SFcalc_obs_bulk_vdw(HKL_data<F_phi<T> >& fphi, const HKL_data<F_sigF<T> >& fsig,
        const Atom_list& atoms, const int n_params = 12);
    bool operator() ( HKL_data<datatypes::F_phi<T> >& fphi,
            const HKL_data<datatypes::F_sigF<T> >& fsig, const Atom_list& atoms );
    const T& bulk_frac() { return bulkfrc; }
    const T& bulk_scale() { return bulkscl; }
    const size_t& n_threads() const { return nthreads; }
    void set_n_threads(size_t n) { nthreads=n; }
private:
    bool initialized_ = false;
    bool optimize_bulk_solvent = false;
    int nparams;
    size_t nthreads;
    T bulkfrc, bulkscl;
    T bulk_u;
}; // class SFcalc_obs_bulk_vdw


} // namespace clipper_cx
