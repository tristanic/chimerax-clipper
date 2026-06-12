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

#include "chimerax_bridge.h"

#include <future>
#include <atomic>
#include <algorithm>

namespace clipper_cx {
namespace bridge {

clipper::Atom_list
clipper_atoms_from_cx_atoms(atomstruct::Atom** cxatoms, size_t n, bool ignore_hydrogens)
{
    auto al = clipper::Atom_list();
    for (size_t i = 0; i < n; ++i)
    {
        auto cxa = cxatoms[i];
        if (ignore_hydrogens && cxa->element().number() == 1) continue;
        if (cxa->element().number() == 0) continue;

        const auto& altlocs = cxa->alt_locs();
        if (altlocs.size())
        {
            for (const auto& altloc: altlocs)
                al.push_back(cl_atom_from_cx_atom<char>(cxa, altloc));
        } else {
            al.push_back(cl_atom_from_cx_atom<>(cxa));
        }
    }
    return al;
}

clipper::Atom_list
clipper_atoms_from_cx_atoms_threaded(atomstruct::Atom** cxatoms, size_t n,
                                     size_t n_threads, bool ignore_hydrogens)
{
    const size_t min_threaded_size   = 10000;
    const size_t min_atoms_per_thread = 4000;
    if (n_threads == 1 || n < min_threaded_size)
        return clipper_atoms_from_cx_atoms(cxatoms, n, ignore_hydrogens);

    size_t atoms_per_thread = std::max(min_atoms_per_thread, n / n_threads + 1);
    std::vector<std::future<clipper::Atom_list>> results;
    std::atomic<size_t> atom_count(0);
    size_t start = 0, end = 0;

    for (size_t i = 0; i < n_threads && end < n; ++i)
    {
        end = std::min(n, start + atoms_per_thread);
        results.push_back(std::async(std::launch::async,
            [&atom_count, ignore_hydrogens](atomstruct::Atom** atoms, size_t s, size_t e)
            {
                auto al = clipper::Atom_list();
                size_t counter = 0;
                for (size_t j = s; j < e; ++j)
                {
                    auto cxa = atoms[j];
                    if (ignore_hydrogens && cxa->element().number() == 1) continue;
                    if (cxa->element().number() == 0) continue;
                    auto altlocs = cxa->alt_locs();
                    if (altlocs.size())
                    {
                        for (const auto& altloc: altlocs)
                        {
                            al.push_back(cl_atom_from_cx_atom<char>(cxa, altloc));
                            counter++;
                        }
                    } else {
                        al.push_back(cl_atom_from_cx_atom<>(cxa));
                        counter++;
                    }
                }
                atom_count += counter;
                return al;
            },
            cxatoms, start, end));
        start += atoms_per_thread;
    }

    auto final_al = clipper::Atom_list();
    // Collect in reverse so the (smaller) last thread finishes first while
    // the larger earlier threads are still running.
    for (auto it = results.rbegin(); it != results.rend(); ++it)
    {
        auto al = it->get();
        std::move(std::begin(al), std::end(al), std::back_inserter(final_al));
    }
    return final_al;
}

std::pair<clipper::Atom_list, std::vector<AtomAltlocIndex>>
clipper_atoms_from_cx_atoms_with_map(
    atomstruct::Atom** cxatoms, size_t n, bool ignore_hydrogens)
{
    clipper::Atom_list al;
    std::vector<AtomAltlocIndex> mapping;
    int idx = 0;
    for (size_t i = 0; i < n; ++i)
    {
        auto cxa = cxatoms[i];
        if (ignore_hydrogens && cxa->element().number() == 1) continue;
        if (cxa->element().number() == 0)                     continue;

        const auto& altlocs = cxa->alt_locs();
        if (altlocs.size())
        {
            for (const auto& altloc : altlocs)
            {
                al.push_back(cl_atom_from_cx_atom<char>(cxa, altloc));
                mapping.push_back({cxa, altloc, idx++, (int)i});
            }
        } else {
            al.push_back(cl_atom_from_cx_atom<>(cxa));
            mapping.push_back({cxa, '\0', idx++, (int)i});
        }
    }
    return {std::move(al), std::move(mapping)};
}

} // namespace bridge
} // namespace clipper_cx
