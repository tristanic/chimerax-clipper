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
#include <utility>
#include <vector>
#include <atomstruct/Atom.h>
#include <atomstruct/Residue.h>
#include <clipper/clipper.h>
#include "imex.h"

namespace clipper_cx {

namespace bridge {

// Template helper — must remain in the header.
template <class ... Types>
clipper::Atom cl_atom_from_cx_atom(atomstruct::Atom* cxatom, Types ... args)
{
    clipper::Atom clatom;
    clatom.set_element(cxatom->element().name());
    const auto& cxcoord = cxatom->coord(args...);
    clatom.set_coord_orth(clipper::Coord_orth(cxcoord[0], cxcoord[1], cxcoord[2]));
    clatom.set_occupancy(cxatom->occupancy(args...));
    clatom.set_u_iso(clipper::Util::b2u(cxatom->bfactor(args...)));

    auto occ = cxatom->occupancy(args...);
    auto u_iso = clipper::Util::b2u(cxatom->bfactor(args...));
    clipper::U_aniso_orth uani;
    if (cxatom->has_aniso_u(args...))
    {
        const auto& cxau = *(cxatom->aniso_u(args...));
        // ChimeraX C++ layer stores aniso_u in row-major order
        uani = clipper::U_aniso_orth(cxau[0],cxau[3],cxau[5],cxau[1],cxau[2],cxau[4]);
    } else {
        uani = clipper::U_aniso_orth(clipper::U_aniso_orth::null());
    }
    clatom.set_u_aniso_orth(uani);
    return clatom;
}

// Non-template conversion functions — defined in chimerax_bridge.cpp and
// exported from clipper_cx so that multiple pybind11 translation units can
// call them without each TU emitting its own duplicate COMDAT symbol.
CLIPPER_CX_IMEX clipper::Atom_list
clipper_atoms_from_cx_atoms(atomstruct::Atom** cxatoms, size_t n, bool ignore_hydrogens);

CLIPPER_CX_IMEX clipper::Atom_list
clipper_atoms_from_cx_atoms_threaded(atomstruct::Atom** cxatoms, size_t n,
                                     size_t n_threads, bool ignore_hydrogens);

//! Maps one entry in the Clipper Atom_list back to its origin in ChimeraX.
//! altloc == '\0' means a single-conformer atom (no altloc).
//! index    : position in the Clipper Atom_list (one per altloc).
//! cx_index : position in the input ChimeraX atom array this entry came from
//!            (lets callers express restraints/flags in input-array terms).
struct CLIPPER_CX_IMEX AtomAltlocIndex {
    atomstruct::Atom* cx_atom;
    char              altloc;
    int               index;
    int               cx_index;
};

//! Like clipper_atoms_from_cx_atoms but also returns a per-entry mapping so
//! that refined values can be written back to the correct (Atom, altloc) pair.
CLIPPER_CX_IMEX
std::pair<clipper::Atom_list, std::vector<AtomAltlocIndex>>
clipper_atoms_from_cx_atoms_with_map(
    atomstruct::Atom** cxatoms, size_t n, bool ignore_hydrogens);

} // namespace bridge
} // namespace clipper_cx
