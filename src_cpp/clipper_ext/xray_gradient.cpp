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

#ifdef _MSC_VER
#pragma warning(disable: 4251)  // STL members of the exported XrayGradientEvaluator
#endif

#include "xray_gradient.h"
#include "edcalc_ext.h"
#include "aniso_scale.h"
#include "vdw.h"
#include <clipper/contrib/sfcalc.h>   // SFcalc_aniso_sum (exact direct-summation Fcalc)
#include <limits>

#include <algorithm>
#include <cmath>
#include <complex>
#include <future>
#include <stdexcept>

namespace clipper_cx {

using namespace clipper;
using namespace clipper::datatypes;

// ============================================================================
// Shared Agarwal per-atom gradient kernel
// ============================================================================

void accumulate_agarwal_gradient(
    const Xmap<ftype32>&                     density_xmap,
    const std::vector<Coord_orth>&           positions,
    const std::vector<String>&               elements,
    const std::vector<double>&               u_iso,
    const std::vector<std::array<double,6>>& u_aniso,
    const std::vector<uint8_t>&              is_aniso,
    const std::vector<double>&               occ,
    const std::vector<double>&               radii,
    const std::vector<AtomShapeFn::TYPE>&    types,
    std::vector<double>&                     raw_grad,
    int                                      n_threads)
{
    const int n_atoms = (int)positions.size();
    const int P       = (int)types.size();
    if (P == 0 || n_atoms == 0) return;

    // Clipper works on the minimal ASU: it reconstructs the full P1 density from
    // these atoms internally, so the gradient integrates the driving density over
    // each ASU atom's own box exactly once — no explicit symmetry-mate summation
    // (that would double-count what Clipper's fft_from/fft_to already expand). The
    // reciprocal per-reflection symmetry weighting is applied in the driving
    // coefficients (see target_and_driving_reciprocal), not here.
    //
    // Split the requested columns into the subset valid for an isotropic atom
    // ({X,Y,Z,Uiso,Occ}: enum value < U11) and the subset valid for an
    // anisotropic atom ({X,Y,Z,Occ,U11..U23}: everything except Uiso).  For each
    // subset keep the parallel list of output-column indices so results map back
    // to the caller's column order.  This is what avoids rho_grad's stale-Uaniso
    // (isotropic branch) and uninitialised-Uiso (anisotropic branch) hazards.
    std::vector<AtomShapeFn::TYPE> iso_types, aniso_types;
    std::vector<int>               iso_cols,  aniso_cols;
    iso_types.reserve(P); aniso_types.reserve(P);
    iso_cols.reserve(P);  aniso_cols.reserve(P);
    for (int c = 0; c < P; ++c) {
        const int t = (int)types[c];
        if (t < (int)AtomShapeFn::U11) { iso_types.push_back(types[c]);  iso_cols.push_back(c); }
        if (t != (int)AtomShapeFn::Uiso) { aniso_types.push_back(types[c]); aniso_cols.push_back(c); }
    }

    const Cell&          cell = density_xmap.cell();
    const Grid_sampling& grid = density_xmap.grid_sampling();

    auto worker = [&](int start, int end) {
        std::vector<ftype> vals;
        using MapRef = Xmap<ftype32>::Map_reference_coord;
        ftype rho;
        for (int j = start; j < end; ++j) {
            const bool aniso = is_aniso[j] != 0;
            const std::vector<AtomShapeFn::TYPE>& jtypes = aniso ? aniso_types : iso_types;
            const std::vector<int>&               jcols  = aniso ? aniso_cols  : iso_cols;
            const int jP = (int)jtypes.size();
            if (jP == 0) continue;

            AtomShapeFn sf;
            if (aniso) {
                const std::array<double,6>& u = u_aniso[j];
                sf = AtomShapeFn(positions[j], elements[j],
                                 U_aniso_orth(u[0], u[1], u[2], u[3], u[4], u[5]),
                                 (ftype)occ[j]);
            } else {
                sf = AtomShapeFn(positions[j], elements[j],
                                 (ftype)u_iso[j], (ftype)occ[j]);
            }
            sf.agarwal_params() = jtypes;

            Grid_range gd(cell, grid, (ftype)radii[j]);
            Coord_grid cg = density_xmap.coord_map(positions[j]).coord_grid();
            Coord_grid g0 = cg + gd.min();
            Coord_grid g1 = cg + gd.max();

            double* out = &raw_grad[(size_t)j * P];
            MapRef iu(density_xmap, g0);
            for (; iu.coord().u() <= g1.u(); iu.next_u()) {
                MapRef iv = iu;
                for (; iv.coord().v() <= g1.v(); iv.next_v()) {
                    MapRef iw = iv;
                    for (; iw.coord().w() <= g1.w(); iw.next_w()) {
                        sf.rho_grad(iw.coord_orth(), rho, vals);
                        const double d = (double)density_xmap[iw];
                        for (int c = 0; c < jP; ++c)
                            out[jcols[c]] -= d * (double)vals[c];
                    }
                }
            }
        }
    };

    const int n_t   = std::max(1, n_threads);
    const int chunk = (n_atoms + n_t - 1) / n_t;
    std::vector<std::future<void>> futures;
    futures.reserve(n_t);
    for (int t = 0; t < n_t; ++t) {
        int s = t * chunk;
        int e = std::min(s + chunk, n_atoms);
        if (s >= n_atoms) break;
        futures.push_back(std::async(std::launch::async, worker, s, e));
    }
    for (auto& f : futures) f.get();
}

// Companion forward-density builder: fills density_xmap (zeroed here) with
// Σ_atoms AtomShapeFn::rho over the SAME per-atom boxes (same radii, same iso/
// aniso ctor selection) that accumulate_agarwal_gradient differentiates.  This
// guarantees the value's density model and the gradient's density model are
// identical — the correctness requirement a finite-difference check exposes.
// Single-threaded: overlapping atom boxes write shared grid points, so this
// avoids the races (and the density discrepancy) of EDcalc_aniso_thread.
static void accumulate_model_density(
    Xmap<ftype32>&                           density_xmap,
    const std::vector<Coord_orth>&           positions,
    const std::vector<String>&               elements,
    const std::vector<double>&               u_iso,
    const std::vector<std::array<double,6>>& u_aniso,
    const std::vector<uint8_t>&              is_aniso,
    const std::vector<double>&               occ,
    const std::vector<double>&               radii)
{
    density_xmap = ftype32(0);
    const int n_atoms = (int)positions.size();
    const Cell&          cell = density_xmap.cell();
    const Grid_sampling& grid = density_xmap.grid_sampling();
    using MapRef = Xmap<ftype32>::Map_reference_coord;
    for (int j = 0; j < n_atoms; ++j) {
        AtomShapeFn sf;
        if (is_aniso[j] != 0) {
            const std::array<double,6>& u = u_aniso[j];
            sf = AtomShapeFn(positions[j], elements[j],
                             U_aniso_orth(u[0], u[1], u[2], u[3], u[4], u[5]),
                             (ftype)occ[j]);
        } else {
            sf = AtomShapeFn(positions[j], elements[j], (ftype)u_iso[j], (ftype)occ[j]);
        }
        Grid_range gd(cell, grid, (ftype)radii[j]);
        Coord_grid cg = density_xmap.coord_map(positions[j]).coord_grid();
        Coord_grid g0 = cg + gd.min();
        Coord_grid g1 = cg + gd.max();
        MapRef iu(density_xmap, g0);
        for (; iu.coord().u() <= g1.u(); iu.next_u()) {
            MapRef iv = iu;
            for (; iv.coord().v() <= g1.v(); iv.next_v()) {
                MapRef iw = iv;
                for (; iw.coord().w() <= g1.w(); iw.next_w())
                    density_xmap[iw] += (ftype)sf.rho(iw.coord_orth());
            }
        }
    }
}

// ============================================================================
// XrayGradientEvaluator
// ============================================================================

struct XrayGradientEvaluator::Impl {
    bool                         realspace = false;
    int                          n_threads = 1;
    std::vector<String>          elements;

    // Reciprocal-space fixed data
    HKL_data<F_sigF<ftype32>>    fobs;
    HKL_data<Phi_fom<ftype32>>   phi_fom;
    HKL_data<Flag>               usage;
    HKL_data<F_phi<ftype32>>     fbulk;
    bool                         has_bulk = false;
    XrayTargetKind               kind     = XrayTargetKind::AmplitudeLS;
    Cell                         cell;
    Grid_sampling                grid_sampling;

    // Real-space fixed data
    Xmap<ftype32>                xmap_target;
    Coord_orth                   target_origin;

    // Working space
    Xmap<ftype32>                density_xmap;
    HKL_data<F_phi<ftype32>>     fcalc_hkl;
    HKL_data<F_phi<ftype32>>     driving_hkl;
    EDcalc_aniso_thread<ftype32> edcalc;
    bool                         threaded_density = true;

    // Cached scale, frozen within a forward/backward and re-fit on refresh_scale:
    //  - reciprocal: per-reflection amplitude scale s(h) (aniso Gaussian x iso
    //    spline, Fc->Fo), indexed by HKL_reference_index.index(); intensity uses s².
    //  - real-space: a single isotropic map scale.
    bool                         has_scale = false;
    double                       scale     = 1.0;   // real-space map scale
    std::vector<double>          refl_scale_;        // reciprocal per-reflection s(h)

    // ---- reciprocal ctor ----
    Impl(const std::vector<String>&        el,
         const HKL_data<F_sigF<ftype32>>&  fo,
         const HKL_data<Phi_fom<ftype32>>& pf,
         const HKL_data<Flag>&             us,
         XrayTargetKind                    kd,
         const HKL_data<F_phi<ftype32>>&   fb,
         int                               nt,
         bool                              threaded)
        : realspace(false), n_threads(std::max(1, nt)), elements(el),
          fobs(fo), phi_fom(pf), usage(us), kind(kd), threaded_density(threaded)
    {
        // A default-constructed (unlinked) HKL_data has parent_hkl_info == NULL;
        // is_null() detects that safely, so guard before calling num_obs() (which
        // walks the parent). Without the is_null() guard the empty-bulk default
        // dereferences null.
        if (!fb.is_null() && fb.num_obs() > 0) { fbulk = fb; has_bulk = true; }
        const HKL_info& hkls = fobs.base_hkl_info();
        cell          = fobs.base_cell();
        grid_sampling = Grid_sampling(hkls.spacegroup(), cell, hkls.resolution());
        density_xmap.init(hkls.spacegroup(), cell, grid_sampling);
        fcalc_hkl   = HKL_data<F_phi<ftype32>>(hkls, cell);
        driving_hkl = HKL_data<F_phi<ftype32>>(hkls, cell);
        edcalc = EDcalc_aniso_thread<ftype32>((size_t)n_threads);
    }

    // ---- real-space ctor ----
    Impl(const std::vector<String>& el,
         const Xmap<ftype32>&       target,
         const Coord_orth&          origin,
         int                        nt,
         bool                       threaded)
        : realspace(true), n_threads(std::max(1, nt)), elements(el),
          threaded_density(threaded), target_origin(origin)
    {
        xmap_target.init(target.spacegroup(), target.cell(), target.grid_sampling());
        for (Xmap<ftype32>::Map_reference_index ix = target.first(); !ix.last(); ix.next())
            xmap_target.set_data(ix.coord(), target[ix]);
        density_xmap.init(target.spacegroup(), target.cell(), target.grid_sampling());
        edcalc = EDcalc_aniso_thread<ftype32>((size_t)n_threads);
    }

    int n_atoms() const { return (int)elements.size(); }

    // Fill the parallel per-atom arrays the density builder and gradient kernel
    // consume (positions frame-shifted by −target_origin in real-space mode).
    void build_atoms(const double* coords, const double* u_iso,
                     const double* u_aniso, const double* occ,
                     const uint8_t* is_aniso,
                     std::vector<Coord_orth>&           positions,
                     std::vector<double>&               u_iso_v,
                     std::vector<std::array<double,6>>& u_aniso_v,
                     std::vector<uint8_t>&              is_aniso_v,
                     std::vector<double>&               occ_v,
                     std::vector<double>&               radii_v) const
    {
        const int N = n_atoms();
        const Coord_orth shift = realspace ? target_origin : Coord_orth(0,0,0);
        positions.resize(N); u_iso_v.resize(N); u_aniso_v.resize(N);
        is_aniso_v.resize(N); occ_v.resize(N); radii_v.resize(N);
        for (int j = 0; j < N; ++j) {
            positions[j] = Coord_orth(coords[3*j], coords[3*j+1], coords[3*j+2]) - shift;
            occ_v[j]     = occ[j];
            is_aniso_v[j] = is_aniso[j];
            double u_eff;
            if (is_aniso[j] != 0) {
                for (int k = 0; k < 6; ++k) u_aniso_v[j][k] = u_aniso[6*j+k];
                u_iso_v[j] = 0.0;
                u_eff = (u_aniso_v[j][0] + u_aniso_v[j][1] + u_aniso_v[j][2]) / 3.0;
            } else {
                u_iso_v[j] = u_iso[j];
                for (int k = 0; k < 6; ++k) u_aniso_v[j][k] = 0.0;
                u_eff = u_iso[j];
            }
            // Per-atom truncation radius matching EDcalc_aniso_thread::cutoff_radius,
            // so the density builder and the gradient integrate over one support.
            const double u_pos = std::max(u_eff, 0.0);
            double vdw;
            try { vdw = data::vdw_radius(std::string(elements[j].c_str())); }
            catch (...) { vdw = 2.0; }
            radii_v[j] = std::max(vdw * (0.4 + 1.5 * std::sqrt(u_pos)), 3.0);
        }
    }

    // Build a Clipper Atom_list from the per-atom arrays, for the threaded
    // EDcalc_aniso_thread forward-density path. Aniso atoms carry both the aniso
    // tensor (density shape) and an equivalent u_iso = trace/3 (only used by
    // EDcalc's cutoff-radius heuristic, matching build_atoms' radii). Iso atoms
    // set u_iso only; their u_aniso_orth stays null (Atom() now nulls it) so
    // EDcalc's iso/aniso fallback selects the isotropic form correctly.
    Atom_list make_atom_list(const std::vector<Coord_orth>&           positions,
                             const std::vector<double>&               u_iso_v,
                             const std::vector<std::array<double,6>>& u_aniso_v,
                             const std::vector<uint8_t>&              is_aniso_v,
                             const std::vector<double>&               occ_v) const
    {
        const int N = n_atoms();
        Atom_list al; al.reserve(N);
        for (int j = 0; j < N; ++j) {
            Atom a;
            a.set_element(elements[j]);
            a.set_coord_orth(positions[j]);
            a.set_occupancy((ftype)occ_v[j]);
            if (is_aniso_v[j] != 0) {
                const std::array<double,6>& u = u_aniso_v[j];
                a.set_u_aniso_orth(U_aniso_orth(u[0], u[1], u[2], u[3], u[4], u[5]));
                a.set_u_iso((ftype)std::max((u[0] + u[1] + u[2]) / 3.0, 0.0));
            } else {
                a.set_u_iso((ftype)u_iso_v[j]);
            }
            al.push_back(a);
        }
        return al;
    }

    // Weight and observed-intensity helpers for the reciprocal target.
    static bool working(const HKL_data<Flag>& usage,
                        const HKL_info::HKL_reference_index& ih)
    {
        return !(usage[ih].missing() || usage[ih].flag() == 0);
    }

    // Per-reflection symmetry-multiplicity weight applied to the driving
    // coefficient (NOT the target). fft_from spreads each ASU reflection over its
    // orbit (num_primops mates + Friedel) via set_hkl; centric reflections (whose
    // Friedel mate coincides with a symmetry mate) and axis reflections (epsilon>1)
    // populate fewer distinct grid points, so their single-count driving
    // coefficient must be boosted to keep the FFT-assembled gradient equal to the
    // exact gradient of the ASU-summed target. epsilonc() = 2*epsilon (centric),
    // epsilon (acentric) is exactly that factor — it is 1 for every reflection in
    // P1 (so P1 is unchanged) and reproduces P-1's empirically-verified factor 2.
    // Finite-difference validated to <1% for both target kinds across
    // P1/P-1/P2/P2_1/P2_12_12_1/C2/C2/c.
    static double eps_weight(const HKL_info::HKL_reference_index& ih)
    {
        return ih.hkl_class().epsilonc();
    }

    // Fit the per-reflection amplitude scale s(h) (Fc->Fo) with the codebase's
    // robust anisotropic-Gaussian x isotropic-spline scaling (scale_fcalc_to_fobs,
    // all reflections -> deterministic, so gradients are reproducible). fcalc_hkl
    // must be current. Cached in refl_scale_ and held fixed until the next refit.
    void refit_scale_reciprocal()
    {
        const HKL_info& hkls = fobs.base_hkl_info();
        HKL_data<F_phi<ftype32>> scaled(hkls, cell);
        U_aniso_orth uaniso;
        std::vector<ftype> aniso_params;
        scale_fcalc_to_fobs<ftype32>(fcalc_hkl, fobs, scaled, uaniso, aniso_params);
        refl_scale_.assign((size_t)hkls.num_reflections(), 1.0);
        for (HKL_info::HKL_reference_index ih = fobs.first(); !ih.last(); ih.next()) {
            if (fcalc_hkl[ih].missing() || scaled[ih].missing()) continue;
            const double fc = fcalc_hkl[ih].f();
            // s(h) is a smooth function of h; recover it as scaled/unscaled (the
            // near-zero-|Fc| reflections it can't resolve contribute ~0 anyway).
            refl_scale_[ih.index()] = (fc > 1e-6) ? (scaled[ih].f() / fc) : 1.0;
        }
        has_scale = true;
    }

    // Fills driving_hkl and returns the target value T (reciprocal mode). Uses the
    // cached per-reflection scale s(h)=refl_scale_[index] (amplitude); intensity
    // uses s(h)². Everywhere the old scalar overall-k appeared, s(h) replaces it.
    double target_and_driving_reciprocal()
    {
        double T = 0.0;
        for (HKL_info::HKL_reference_index ih = fobs.first(); !ih.last(); ih.next()) {
            if (fobs[ih].missing() || !working(usage, ih)) {
                driving_hkl.set_data(ih.hkl(), F_phi<ftype32>());
                continue;
            }
            const double fc  = fcalc_hkl[ih].f();
            const double phi = fcalc_hkl[ih].phi();
            const double fo  = fobs[ih].f();
            const double sf  = fobs[ih].sigf();
            if (sf <= 0.0) { driving_hkl.set_data(ih.hkl(), F_phi<ftype32>()); continue; }
            const double s = refl_scale_[ih.index()];
            // Symmetry-multiplicity weight for the driving coefficient only; the
            // target T remains the plain ASU sum (each unique reflection once).
            const double ew = eps_weight(ih);

            double coeff_mag;
            if (kind == XrayTargetKind::AmplitudeLS) {
                const double m = phi_fom[ih].fom();
                const double w = 1.0 / (sf * sf);
                const double resid = s * fc - m * fo;
                T += 0.5 * w * resid * resid;
                // G(h) = ε_c·s(h)·w·(m|Fo| − s(h)|Fc|)·exp(iφc)
                coeff_mag = ew * s * w * (m * fo - s * fc);
            } else {
                const double Io   = fo * fo;
                const double sigI = 2.0 * fo * sf;
                if (sigI <= 0.0) { driving_hkl.set_data(ih.hkl(), F_phi<ftype32>()); continue; }
                const double w   = 1.0 / (sigI * sigI);
                const double s2  = s * s;                 // intensity scale = s(h)²
                const double fc2 = fc * fc;
                const double resid = s2 * fc2 - Io;
                T += 0.5 * w * resid * resid;
                // G(h) = ε_c·2·s²·w·(Io − s²|Fc|²)·|Fc|·exp(iφc)
                coeff_mag = ew * 2.0 * s2 * w * (Io - s2 * fc2) * fc;
            }
            std::complex<ftype32> coeff =
                ftype32(coeff_mag) * std::polar(ftype32(1.0f), ftype32(phi));
            driving_hkl.set_data(ih.hkl(), F_phi<ftype32>(coeff));
        }
        return T;
    }

    int n_reflections() const
    {
        return realspace ? 0 : fobs.base_hkl_info().num_reflections();
    }

    // Forward Fcalc -> scale, then emit (Fo, s·|Fc|) per reflection in HKL_data index
    // order. Off the measured/working set (or where Fcalc is missing) the pair is NaN,
    // so the caller masks exactly as the live small-molecule map does before handing to
    // reflection_tools.compute_r_factors. Reciprocal (fobs) mode only. FORWARD-ONLY (no
    // gradient) — the gradient path (evaluate()) is untouched and always FFT.
    //
    // `use_summation`: false -> the FFT Fcalc the loss/live map use (EDcalc density ->
    // fft_to); true -> EXACT direct summation (SFcalc_aniso_sum), which matches
    // io/small_molecule recomputed_r_factor and removes the FFT grid approximation that
    // costs ~0.01 in R on heavy-scatterer crystals. Both branches feed the SAME atoms
    // (build_atoms) through the SAME scale_fcalc_to_fobs, so they differ by ONE thing:
    // how Fcalc is computed. XRAY throughout, matching the evaluator's density model.
    //
    // SIDE-EFFECT-FREE: fits a LOCAL scale (does NOT touch refl_scale_/has_scale, which
    // the value_and_gradient loss path caches and reuses between frozen-scale steps), so
    // a metric call can be interleaved with training without repointing the loss scale.
    void fobs_scaled_fcalc(const double* coords, const double* u_iso,
                           const double* u_aniso, const double* occ,
                           const uint8_t* is_aniso, bool use_summation,
                           double* out_fo, double* out_sfc)
    {
        if (realspace)
            throw std::runtime_error(
                "fobs_scaled_fcalc: reciprocal (structure-factor) mode only");
        std::vector<Coord_orth>           positions;
        std::vector<double>               u_iso_v, occ_v, radii_v;
        std::vector<std::array<double,6>> u_aniso_v;
        std::vector<uint8_t>              is_aniso_v;
        build_atoms(coords, u_iso, u_aniso, occ, is_aniso, positions,
                    u_iso_v, u_aniso_v, is_aniso_v, occ_v, radii_v);
        Atom_list atoms = make_atom_list(positions, u_iso_v, u_aniso_v, is_aniso_v, occ_v);
        if (use_summation) {
            SFcalc_aniso_sum<ftype32> sfc;   // XRAY (matches the density path's AtomShapeFn)
            sfc(fcalc_hkl, atoms);
        } else {
            density_xmap = ftype32(0);
            edcalc(density_xmap, atoms);
            density_xmap.fft_to(fcalc_hkl, n_threads);
        }
        if (has_bulk) {
            for (HKL_info::HKL_reference_index ih = fobs.first(); !ih.last(); ih.next()) {
                if (fcalc_hkl[ih].missing() || fbulk[ih].missing()) continue;
                std::complex<ftype32> ft =
                    (std::complex<ftype32>)fcalc_hkl[ih] + (std::complex<ftype32>)fbulk[ih];
                fcalc_hkl.set_data(ih.hkl(), F_phi<ftype32>(ft));
            }
        }
        // Local aniso-Gaussian x iso-spline scale (the SAME scale_fcalc_to_fobs that
        // recomputed_r_factor and the live map use); scaled[ih].f() == s(h)·|Fc|.
        HKL_data<F_phi<ftype32>> scaled(fobs.base_hkl_info(), cell);
        U_aniso_orth      scale_u;
        std::vector<ftype> scale_params;
        scale_fcalc_to_fobs<ftype32>(fcalc_hkl, fobs, scaled, scale_u, scale_params);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (HKL_info::HKL_reference_index ih = fobs.first(); !ih.last(); ih.next()) {
            const size_t i = (size_t)ih.index();
            const bool measured = !fobs[ih].missing() && working(usage, ih)
                                  && fobs[ih].sigf() > 0.0;
            if (!measured || scaled[ih].missing()) {
                out_fo[i] = nan; out_sfc[i] = nan; continue;
            }
            out_fo[i]  = fobs[ih].f();
            out_sfc[i] = scaled[ih].f();
        }
    }

    double evaluate(const double* coords, const double* u_iso,
                    const double* u_aniso, const double* occ,
                    const uint8_t* is_aniso,
                    const std::vector<AtomShapeFn::TYPE>& selected,
                    double* out_grad, bool refresh_scale)
    {
        const int N = n_atoms();
        const int P = (int)selected.size();
        if (P == 0)
            throw std::runtime_error("XrayGradientEvaluator: empty parameter selection");

        std::vector<Coord_orth>           positions;
        std::vector<double>               u_iso_v, occ_v, radii_v;
        std::vector<std::array<double,6>> u_aniso_v;
        std::vector<uint8_t>              is_aniso_v;
        build_atoms(coords, u_iso, u_aniso, occ, is_aniso, positions,
                    u_iso_v, u_aniso_v, is_aniso_v, occ_v, radii_v);

        // ρ_calc on the working grid, from the SAME density model the gradient
        // kernel differentiates (AtomShapeFn.rho over the per-atom boxes). Default:
        // the threaded EDcalc_aniso_thread (near-linear scaling to ~20 cores);
        // fallback: the single-threaded companion. Both produce identical density
        // (same AtomShapeFn, same cutoff radius), so the gradient stays consistent.
        if (threaded_density) {
            Atom_list atoms = make_atom_list(positions, u_iso_v, u_aniso_v,
                                             is_aniso_v, occ_v);
            density_xmap = ftype32(0);
            edcalc(density_xmap, atoms);
        } else {
            accumulate_model_density(density_xmap, positions, elements, u_iso_v,
                                     u_aniso_v, is_aniso_v, occ_v, radii_v);
        }

        double T = 0.0;
        if (!realspace) {
            density_xmap.fft_to(fcalc_hkl, n_threads);
            if (has_bulk) {
                for (HKL_info::HKL_reference_index ih = fobs.first(); !ih.last(); ih.next()) {
                    if (fcalc_hkl[ih].missing() || fbulk[ih].missing()) continue;
                    std::complex<ftype32> ft =
                        (std::complex<ftype32>)fcalc_hkl[ih] + (std::complex<ftype32>)fbulk[ih];
                    fcalc_hkl.set_data(ih.hkl(), F_phi<ftype32>(ft));
                }
            }
            if (refresh_scale || !has_scale) refit_scale_reciprocal();
            T = target_and_driving_reciprocal();
            density_xmap.fft_from(driving_hkl, n_threads);   // density_xmap := d(x)
        } else {
            // Isotropic map scale k = Σ(ρ_t·ρ_c)/Σρ_c²  (frozen when refresh_scale=false).
            if (refresh_scale || !has_scale) {
                double sc = 0.0, sq = 0.0;
                for (Xmap<ftype32>::Map_reference_index ix = density_xmap.first(); !ix.last(); ix.next()) {
                    const double rc = density_xmap[ix];
                    sc += xmap_target[ix] * rc;
                    sq += rc * rc;
                }
                scale = (sq > 1e-20) ? sc / sq : 1.0;
                has_scale = true;
            }
            const double k = scale;
            for (Xmap<ftype32>::Map_reference_index ix = density_xmap.first(); !ix.last(); ix.next()) {
                const double rc   = density_xmap[ix];
                const double rt   = xmap_target[ix];
                const double diff = k * rc - rt;
                T += 0.5 * diff * diff;
                density_xmap[ix] = (ftype)(k * (rt - k * rc));   // D(x)
            }
        }

        std::vector<double> raw_grad((size_t)N * P, 0.0);
        accumulate_agarwal_gradient(density_xmap, positions, elements, u_iso_v,
                                    u_aniso_v, is_aniso_v, occ_v, radii_v,
                                    selected, raw_grad, n_threads);

        // Reciprocal route: the gradient is assembled by back-FFT (fft_from) of the
        // reciprocal residual, whose Clipper normalization differs from the direct
        // reciprocal-space target T by the geometric factor V²/(2·N_grid):
        //   V/N   — real-space grid quadrature,
        //   1/V   — Clipper's fft_from convention (extra power of V),
        //   ½     — Friedel pairing (the ASU sum counts each |F| once, the real FFT
        //           distributes over ±h).
        // This factor is spacegroup-INDEPENDENT: the per-reflection symmetry
        // multiplicity is carried exactly by the epsilonc() weight on each driving
        // coefficient (see eps_weight / target_and_driving_reciprocal), so no n_sym
        // appears here (a former n_sym factor was a wrong stand-in for it). The value
        // T is the ASU reflection-sum least-squares as defined; its true gradient is
        // this factor times the (epsilonc-weighted) back-FFT Agarwal sum over the one
        // ASU. The real-space route needs no such factor. Validated by finite
        // difference to <1% for both target kinds across
        // P1/P-1/P2/P2_1/P2_12_12_1/C2/C2/c.
        if (!realspace) {
            const double V     = cell.volume();
            const double grad_scale = V * V / (2.0 * (double)grid_sampling.size());
            for (double& g : raw_grad) g *= grad_scale;
        }

        std::copy(raw_grad.begin(), raw_grad.end(), out_grad);
        return T;
    }
};

XrayGradientEvaluator::XrayGradientEvaluator(
    const std::vector<String>&        elements,
    const HKL_data<F_sigF<ftype32>>&  fobs,
    const HKL_data<Phi_fom<ftype32>>& phi_fom,
    const HKL_data<Flag>&             usage,
    XrayTargetKind                    kind,
    const HKL_data<F_phi<ftype32>>&   f_bulk,
    int                               n_threads,
    bool                              threaded_density)
    : p_(new Impl(elements, fobs, phi_fom, usage, kind, f_bulk, n_threads,
                  threaded_density))
{}

XrayGradientEvaluator::XrayGradientEvaluator(
    const std::vector<String>& elements,
    const Xmap<ftype32>&       target_map,
    const Coord_orth&          target_origin,
    int                        n_threads,
    bool                       threaded_density)
    : p_(new Impl(elements, target_map, target_origin, n_threads, threaded_density))
{}

XrayGradientEvaluator::~XrayGradientEvaluator() = default;

int XrayGradientEvaluator::n_atoms() const { return p_->n_atoms(); }

int XrayGradientEvaluator::n_reflections() const { return p_->n_reflections(); }

void XrayGradientEvaluator::fobs_scaled_fcalc(
    const double* coords, const double* u_iso, const double* u_aniso,
    const double* occ, const uint8_t* is_aniso, bool use_summation,
    double* out_fo, double* out_scaled_fcalc)
{
    p_->fobs_scaled_fcalc(coords, u_iso, u_aniso, occ, is_aniso, use_summation,
                          out_fo, out_scaled_fcalc);
}

double XrayGradientEvaluator::value_and_gradient(
    const double*                         coords,
    const double*                         u_iso,
    const double*                         u_aniso,
    const double*                         occ,
    const uint8_t*                        is_aniso,
    const std::vector<AtomShapeFn::TYPE>& selected,
    double*                               out_grad,
    bool                                  refresh_scale)
{
    return p_->evaluate(coords, u_iso, u_aniso, occ, is_aniso, selected,
                        out_grad, refresh_scale);
}

} // namespace clipper_cx
