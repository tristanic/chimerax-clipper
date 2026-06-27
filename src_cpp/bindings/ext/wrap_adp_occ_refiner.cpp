/**
 * @Author: Tristan Croll <tcroll@altoslabs.com>
 * Pybind11 bindings for BFactorOccRefinerThread — maximum-likelihood
 * isotropic B-factor and occupancy refinement.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <clipper/clipper.h>
#include <clipper_ext/adp_occ_refiner.h>
// chimerax_bridge.h defines clipper_atoms_from_cx_atoms as inline, so including
// it here is safe alongside wrap_xtal_mgr.cpp in the same DSO.
#include <clipper_ext/chimerax_bridge.h>

namespace py = pybind11;
namespace cx = clipper_cx;
using namespace clipper;
using namespace clipper::datatypes;

// ---------------------------------------------------------------------------

// EqualOccGroup / OccConstraintGroup / BFactorRestraint / BFactorTargetRestraint
// are engine-internal (Clipper-indexed) types only.  They are not exposed to
// Python: callers supply restraints as ChimeraX-indexed arrays to launch() /
// launch_realspace() and occupancy groups are auto-derived in C++.

void declare_refine_config(py::module& m)
{
    using Class = cx::RefineConfig;
    py::class_<Class>(m, "RefineConfig",
        "Configuration for ML B-factor and/or occupancy refinement.\n"
        "\n"
        "Holds the scalar refinement settings (refine_b, max_cycles, b_min/b_max,\n"
        "L-BFGS-B tolerances, n_threads) plus the per-atom occupancy-refinement\n"
        "flags `refine_occ`.\n"
        "\n"
        "Everything atom-specific is expressed in terms of the INPUT ChimeraX\n"
        "atom array (positions in the atoms passed to launch/launch_realspace),\n"
        "and C++ resolves it to the internal Clipper Atom_list at launch:\n"
        "  - `refine_occ`: per-ChimeraX-atom flags (length = input array size).\n"
        "  - B-factor restraints: passed as ChimeraX-index + altloc arrays to\n"
        "    launch()/launch_realspace() (not stored here).\n"
        "  - occupancy groups: auto-derived from covalent-fragment + altloc\n"
        "    topology when refine_occ is set (not built by the caller).")
        .def(py::init<>())
        .def_readwrite("refine_b",        &Class::refine_b,
            "bool: refine isotropic B-factors (default True)")
        // refine_occ is a per-atom uint8 array. Accept numpy bool/uint8 arrays
        // or Python sequences from the Python side.
        .def_property("refine_occ",
            [](const Class& self) -> py::array_t<uint8_t> {
                py::array_t<uint8_t> arr((py::ssize_t)self.refine_occ.size());
                auto buf = arr.request();
                std::copy(self.refine_occ.begin(), self.refine_occ.end(),
                          static_cast<uint8_t*>(buf.ptr));
                return arr;
            },
            [](Class& self, py::array_t<uint8_t> arr) {
                auto buf = arr.request();
                self.refine_occ.resize((size_t)buf.size);
                std::copy(static_cast<uint8_t*>(buf.ptr),
                          static_cast<uint8_t*>(buf.ptr) + buf.size,
                          self.refine_occ.begin());
            },
            "numpy.uint8 array: per-atom occupancy-refinement flags (0=fixed, "
            "1=refine), indexed by position in the ChimeraX atom array passed to "
            "launch()/launch_realspace() (length = that array's size).  Empty "
            "(default) means no occupancy refinement.  When set, occupancy groups "
            "are auto-derived in C++ from covalent-fragment + altloc topology.")
        .def_readwrite("use_curvature",    &Class::use_curvature,
            "bool: use second derivatives for preconditioning (reserved)")
        .def_readwrite("max_cycles",       &Class::max_cycles,
            "int: maximum L-BFGS-B iterations")
        .def_readwrite("b_min",            &Class::b_min,
            "float: lower bound on isotropic B-factors (A^2)")
        .def_readwrite("b_max",            &Class::b_max,
            "float: upper bound on isotropic B-factors (A^2)")
        .def_readwrite("lbfgs_epsilon",    &Class::lbfgs_epsilon,
            "float: gradient-norm convergence threshold")
        .def_readwrite("lbfgs_past",       &Class::lbfgs_past,
            "int: number of past iterations for the delta-f stopping criterion")
        .def_readwrite("lbfgs_delta",      &Class::lbfgs_delta,
            "float: relative delta-f convergence threshold (0 disables)")
        .def_readwrite("n_threads",        &Class::n_threads,
            "int: thread count for EDcalc and FFT steps")
        // b_restraints / b_target_restraints / occupancy groups are no longer
        // caller-set: restraints are supplied (ChimeraX-indexed + altloc) as
        // launch() / launch_realspace() arguments and occupancy groups are
        // auto-derived.  These RefineConfig fields remain as the engine-internal
        // Clipper-indexed representation, populated at launch.
        ;
}

// ---------------------------------------------------------------------------
// Restraint specs cross the binding as a single Python tuple per group (rather
// than many parallel-array arguments) so each bound method's argument count stays
// well under MSVC's pybind11 template-instantiation ceiling (C1202 with the v141
// toolset).  None → no restraints; a trailing `alphas` element is optional
// (omitted → the engine defaults each restraint to the Charbonnier shape α=1).
//   pairwise tuple : (atoms1, altlocs1, atoms2, altlocs2, sigmas, ks[, alphas])
//   target tuple   : (atoms,  altlocs,  target_us,        sigmas, ks[, alphas])
static cx::BFactorRestraintSpec parse_pairwise_restraints(py::object o)
{
    cx::BFactorRestraintSpec s;
    if (o.is_none()) return s;
    py::sequence t = o.cast<py::sequence>();
    s.atoms1   = t[0].cast<std::vector<int>>();
    s.altlocs1 = t[1].cast<std::vector<std::string>>();
    s.atoms2   = t[2].cast<std::vector<int>>();
    s.altlocs2 = t[3].cast<std::vector<std::string>>();
    s.sigmas   = t[4].cast<std::vector<double>>();
    s.weights  = t[5].cast<std::vector<double>>();
    if (py::len(t) > 6) s.alphas = t[6].cast<std::vector<double>>();
    return s;
}

static cx::BFactorTargetRestraintSpec parse_target_restraints(py::object o)
{
    cx::BFactorTargetRestraintSpec s;
    if (o.is_none()) return s;
    py::sequence t = o.cast<py::sequence>();
    s.atoms     = t[0].cast<std::vector<int>>();
    s.altlocs   = t[1].cast<std::vector<std::string>>();
    s.target_us = t[2].cast<std::vector<double>>();
    s.sigmas    = t[3].cast<std::vector<double>>();
    s.weights   = t[4].cast<std::vector<double>>();
    if (py::len(t) > 5) s.alphas = t[5].cast<std::vector<double>>();
    return s;
}

void declare_adp_occ_refiner(py::module& m)
{
    using Class = cx::BFactorOccRefinerThread;
    py::class_<Class>(m, "BFactorOccRefinerThread")
        .def(py::init<const cx::RefineConfig&>(),
             py::arg("config"))

        // Launch refinement.  Uses clipper_atoms_from_cx_atoms_with_map so that
        // per-altloc results can be written back correctly by apply_to_atoms().
        // f_bulk: optional bulk-solvent contribution (Xtal_mgr::f_bulk) from the
        // live map manager.  When provided it is held fixed and added to
        // F_atoms(current U) at each step, so the gradient correctly accounts for
        // the bulk-solvent contribution.
        // Pairwise B-factor restraints are supplied as a single Python tuple
        // `restraints` (or None) — see parse_pairwise_restraints above for its
        // shape.  Indices reference the input ChimeraX atom array, each endpoint
        // qualified by an altloc string ("" = no altloc); C++ resolves them to
        // Clipper indices.
        .def("launch",
            [](Class& self,
               std::vector<uintptr_t>                   cx_ptrs,
               const HKL_data<F_sigF<ftype32>>&         fobs,
               const HKL_data<Phi_fom<ftype32>>&        phi_fom,
               const HKL_data<Flag>&                    usage,
               bool                                     ignore_hydrogens,
               const HKL_data<F_phi<ftype32>>&          f_bulk,
               py::object                               restraints)
            {
                auto result =
                    clipper_cx::bridge::clipper_atoms_from_cx_atoms_with_map(
                        (atomstruct::Atom**)cx_ptrs.data(),
                        cx_ptrs.size(), ignore_hydrogens);
                cx::BFactorRestraintSpec rspec = parse_pairwise_restraints(restraints);
                self.launch(result.first, result.second, fobs, phi_fom, usage,
                            ignore_hydrogens, f_bulk, rspec);
            },
            py::arg("atoms"),
            py::arg("fobs"),
            py::arg("phi_fom"),
            py::arg("usage"),
            py::arg("ignore_hydrogens") = false,
            py::arg("f_bulk") = HKL_data<F_phi<ftype32>>(),
            py::arg("restraints")  = py::none())

        // Write refined B-factors and occupancies back to ChimeraX atoms via
        // the stored (Atom*, altloc, index) mapping.  Must be called from the
        // main thread after verifying atom validity.
        .def("apply_to_atoms", &Class::apply_to_atoms)

        // Real-space LS launch: refine against a fixed P1 Xmap target density.
        // context_ptrs: additional atoms contributing to ρ_calc but not refined.
        // target_map: P1 Xmap (built Python-side from the ChimeraX Volume subregion).
        // target_origin: Coord_orth giving the P1 cell origin in the original frame.
        // Restraints are supplied as single Python tuples (or None): `restraints`
        // (pairwise) and `target_restraints` (one-sided) — see the parse helpers
        // above for their shapes.  Indices reference the input refined_atoms array,
        // each endpoint qualified by an altloc string ("" = no altloc); C++
        // resolves them to Clipper indices.
        .def("launch_realspace",
            [](Class& self,
               std::vector<uintptr_t>    refined_ptrs,
               std::vector<uintptr_t>    context_ptrs,
               const Xmap<ftype32>&      target_map,
               const Coord_orth&         target_origin,
               bool                      ignore_hydrogens,
               py::object                restraints,
               py::object                target_restraints)
            {
                auto refined_result =
                    clipper_cx::bridge::clipper_atoms_from_cx_atoms_with_map(
                        (atomstruct::Atom**)refined_ptrs.data(),
                        refined_ptrs.size(), ignore_hydrogens);
                auto context_atoms = clipper_cx::bridge::clipper_atoms_from_cx_atoms(
                    (atomstruct::Atom**)context_ptrs.data(),
                    context_ptrs.size(), ignore_hydrogens);
                cx::BFactorRestraintSpec rspec = parse_pairwise_restraints(restraints);
                cx::BFactorTargetRestraintSpec tspec = parse_target_restraints(target_restraints);
                self.launch_realspace(refined_result.first, refined_result.second,
                                      context_atoms, target_map, target_origin,
                                      ignore_hydrogens, rspec, tspec);
            },
            py::arg("refined_atoms"),
            py::arg("context_atoms"),
            py::arg("target_map"),
            py::arg("target_origin"),
            py::arg("ignore_hydrogens") = false,
            py::arg("restraints")        = py::none(),
            py::arg("target_restraints") = py::none())

        // Compute R-work and R-free from the current refined parameters.
        // Runs one extra EDcalc + FFT; only meaningful for the crystallographic
        // path (returns (-1, -1) in real-space mode).  Blocks until the
        // refinement thread has finished.
        .def("compute_rfactors",
            [](Class& self) -> py::tuple {
                auto r = self.compute_rfactors();
                return py::make_tuple(r.first, r.second);
            })

        // Standard R-work/R-free at the initial (pre-refinement) parameters — same
        // metric as compute_rfactors(), for a like-for-like before/after comparison.
        .def("initial_rfactors",
            [](Class& self) -> py::tuple {
                auto r = self.initial_rfactors();
                return py::make_tuple(r.first, r.second);
            })

        // Synchronous diagnostic: compute ρ_calc on the P1 grid (no optimisation).
        // Returns a Python Xmap_float object (same grid as target_map) so the
        // caller can pass it directly to _debug_show_p1_target() alongside the
        // target density for a side-by-side registration check.
        .def("compute_realspace_density",
            [](Class& self,
               std::vector<uintptr_t>    refined_ptrs,
               std::vector<uintptr_t>    context_ptrs,
               const Xmap<ftype32>&      target_map,
               const Coord_orth&         target_origin,
               bool                      ignore_hydrogens) -> Xmap<ftype32>*
            {
                auto refined_result =
                    clipper_cx::bridge::clipper_atoms_from_cx_atoms_with_map(
                        (atomstruct::Atom**)refined_ptrs.data(),
                        refined_ptrs.size(), ignore_hydrogens);
                auto context_atoms =
                    clipper_cx::bridge::clipper_atoms_from_cx_atoms(
                        (atomstruct::Atom**)context_ptrs.data(),
                        context_ptrs.size(), ignore_hydrogens);
                return self.compute_realspace_density(
                    refined_result.first, refined_result.second,
                    context_atoms, target_map, target_origin);
            },
            py::arg("refined_atoms"),
            py::arg("context_atoms"),
            py::arg("target_map"),
            py::arg("target_origin"),
            py::arg("ignore_hydrogens") = false,
            py::return_value_policy::take_ownership)

        // Signal the running L-BFGS-B thread to stop at its next evaluation.
        .def("cancel", &Class::cancel)

        .def_property_readonly("thread_running", &Class::thread_running)
        .def("ready",          &Class::ready)
        .def_property_readonly("n_atoms",        &Class::n_atoms)

        // Return refined isotropic U values (Å²) as a 1-D float64 numpy array.
        // Blocks if the refinement thread has not yet finished.
        .def_property_readonly("refined_u_iso",
            [](Class& self) -> py::array_t<double> {
                const auto& v = self.refined_u_iso();
                py::array_t<double> arr((py::ssize_t)v.size());
                auto buf = arr.request();
                std::copy(v.begin(), v.end(),
                          static_cast<double*>(buf.ptr));
                return arr;
            })

        // Return refined occupancies as a 1-D float64 numpy array.
        .def_property_readonly("refined_occ",
            [](Class& self) -> py::array_t<double> {
                const auto& v = self.refined_occ();
                py::array_t<double> arr((py::ssize_t)v.size());
                auto buf = arr.request();
                std::copy(v.begin(), v.end(),
                          static_cast<double*>(buf.ptr));
                return arr;
            })
        ;
}

void init_adp_occ_refiner(py::module& m)
{
    declare_refine_config(m);
    declare_adp_occ_refiner(m);
}
