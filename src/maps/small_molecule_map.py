# Clipper plugin to UCSF ChimeraX
# Copyright (C) 2016-2019 Tristan Croll, University of Cambridge
#
# Live electron-density maps for small-molecule (COD) crystals. Pure-Python
# manager that duck-types the subset of the C++ Xtal_thread_mgr interface used by
# XmapSet/XmapHandler_Live, so it drops into the existing live-map display,
# spotlight-masking and trigger machinery.
#
# Structure factors are computed by FFT (SFcalc_aniso_fft, rate 1.5) with NO bulk
# solvent - correct for densely-packed small-molecule crystals and, with the
# corrected atom list (special positions, Clipper frame, a*-scaled ADPs), as
# accurate as exact summation to within ~0.003 in R while being far faster. Map
# coefficients are the small-molecule standard 2Fo-Fc / Fo-Fc with phases from
# Fcalc (no sigma-A weighting, no R-free set). The heavy compute runs in a worker
# thread; the SFcalc and fft_from bindings release the GIL, so graphics stay
# responsive during recalculation.

import numpy


class SmallMoleculeXmapMgr:
    # FFT oversampling rate for Fcalc (1.5 is the accuracy/speed sweet spot for
    # small molecules; higher rates are slower with no meaningful accuracy gain).
    FFT_RATE = 1.5

    def __init__(self, hklinfo, cell, spacegroup, grid_sampling, fobs, structure,
                 radiation='xray'):
        '''
        Args:
            hklinfo, cell, spacegroup, grid_sampling: the crystal definition.
            radiation: 'xray' or 'electron' - selects the scattering-factor table
                for Fcalc (electron for micro-ED / 3D-ED). Fixed at construction and
                read only in the worker thread, so it is race-free.
            fobs: an HKL_data_F_sigF of observed amplitudes (Fo = sqrt of
                intensity), over hklinfo. Used directly by the anisotropic+spline
                scaling (scale_fcalc_to_fobs).
            structure: the live AtomicStructure. Its per-atom coordinates, occupancy,
                B-factor, anisotropic U and (via the clipper_scattering_species custom
                attribute) scattering species are read fresh on each recompute, so the
                map reflects every model edit. The model is expected to carry Clipper-
                frame coords and orthogonal ADPs (see io.small_molecule.hydrate_...).
        '''
        self.hklinfo = hklinfo
        self.cell = cell
        self.spacegroup = spacegroup
        self.grid_sampling = grid_sampling
        self._structure = structure
        self._radiation = radiation
        self._fobs = fobs            # HKL_data_F_sigF, for aniso+spline scaling
        self._fo = numpy.asarray(fobs.data[1][:, 0], numpy.double)  # aligned to hklinfo
        self._meas = numpy.isfinite(self._fo)
        self._maps = {}            # name -> dict(is_difference, current, pending)
        self._thread = None
        self._ready = False
        self._rwork = 0.0
        self._n_reflections = int(self._meas.sum())
        # Reused only for its O(1) site-multiplicity lookup (special-position occupancy).
        from ..clipper_python import Xmap_double
        self._mult_xmap = Xmap_double(spacegroup, cell, grid_sampling)

    # ---- interface used by XmapSet / XmapHandler_Live (duck-types Xtal_thread_mgr) ----

    @property
    def thread_running(self):
        return self._thread is not None and self._thread.is_alive()

    def ready(self):
        return self._ready and not self.thread_running

    @property
    def rwork(self):
        return self._rwork

    @property
    def rfree(self):
        # Small-molecule data carries no R-free set; report rwork for the status line.
        return self._rwork

    def add_xmap(self, name, b_sharp=0, is_difference_map=False, **kw):
        from ..clipper_python import Xmap_float
        self._maps[name] = {
            'is_difference': bool(is_difference_map),
            'current': Xmap_float(self.spacegroup, self.cell, self.grid_sampling),
            'pending': None,
        }
        # Compute an initial map so something shows immediately.
        self.recalculate_now()

    def delete_xmap(self, name):
        self._maps.pop(name, None)

    def get_xmap_ref(self, name):
        return self._maps[name]['current']

    def get_map_stats(self, name):
        from ..clipper_python import Map_stats
        return Map_stats(self._maps[name]['current'])

    def init(self, *args):
        # Parity with Xtal_thread_mgr.init(atoms.pointers); maps are computed when
        # added, so nothing extra is needed here.
        pass

    def bulk_solvent_optimization_needed(self):
        # No-op: small-molecule crystals are densely packed and have no bulk
        # solvent to (re-)optimise. Present for Xtal_thread_mgr interface parity
        # (XmapSet calls it whenever the model changes).
        pass

    def recalculate_all_maps(self, *args):
        '''Start a background recalculation. Coordinates are snapshotted on the
        calling (main) thread; only the GIL-releasing compute runs in the worker.'''
        if self.thread_running:
            return
        atoms = self._current_atom_list()
        self._ready = False
        import threading
        self._thread = threading.Thread(target=self._worker, args=(atoms,), daemon=True)
        self._thread.start()

    def apply_new_maps(self):
        if self.thread_running:
            self._thread.join()
        for m in self._maps.values():
            if m['pending'] is not None:
                m['current'] = m['pending']
                m['pending'] = None
        self._ready = False

    # ---- internals ----

    def _current_atom_list(self):
        '''Build the Clipper Atom_list fresh from the LIVE model each recompute, so the
        map reflects every edit (coordinates, occupancy, B-factor, ADP). Called on the
        main thread (all reads are plain numpy); the immutable result goes to the worker.'''
        from ..clipper_python import Atom_list, Coord_orth
        atoms = self._structure.atoms
        # Exclude symmetry-completed atoms (see io.fragments): Clipper's SFcalc generates
        # them from their ASU source by crystallographic symmetry, so including them would
        # double-count.
        keep = numpy.array([not getattr(a, 'clipper_sf_exclude', False) for a in atoms], bool)
        if not keep.all():
            atoms = atoms.filter(keep)
        n = len(atoms)
        xyz = numpy.array(atoms.coords, numpy.double)
        # Ionic-aware scattering species (falls back to the neutral element for atoms
        # without the attribute, e.g. user-added atoms).
        elements = [getattr(a, 'clipper_scattering_species', None) or a.element.name
                    for a in atoms]
        # Occupancy with live special-position correction (Clipper's SFcalc applies every
        # symop, so an atom on a special position is otherwise multiply counted).
        occ = numpy.array(atoms.occupancies, numpy.double)
        cell, grid, xm = self.cell, self.grid_sampling, self._mult_xmap
        for i in range(n):
            cf = Coord_orth(xyz[i, 0], xyz[i, 1], xyz[i, 2]).coord_frac(cell)
            m = xm.multiplicity(cf.coord_grid(grid))
            if m > 1:
                occ[i] /= m
        # Isotropic U from B (U = B / 8 pi^2); anisotropic U where present (NaN elsewhere,
        # the per-atom iso/aniso convention Atom_list expects).
        u_iso = numpy.array(atoms.bfactors, numpy.double) / (8.0 * numpy.pi ** 2)
        u_aniso = numpy.full((n, 6), numpy.nan, numpy.double)
        has_aniso = atoms.has_aniso_u
        if has_aniso.any():
            u_aniso[has_aniso] = atoms.filter(has_aniso).aniso_u6
        return Atom_list(elements, xyz, occ, u_iso, u_aniso)

    def recalculate_now(self):
        '''Synchronous recalculation (used for the initial map and headless tests).'''
        self._compute(self._current_atom_list())
        for m in self._maps.values():
            if m['pending'] is not None:
                m['current'] = m['pending']
                m['pending'] = None
        self._ready = False

    def _worker(self, atoms):
        try:
            self._compute(atoms)
        finally:
            self._ready = True

    def _compute(self, atoms):
        '''Fcalc (FFT, no bulk solvent) -> scale -> 2Fo-Fc / Fo-Fc coefficients ->
        real-space maps. Runs in the worker thread; SFcalc and fft_from release the
        GIL so graphics keep rendering.'''
        # Float (ftype32) throughout, to match the live-map display path: the
        # XmapHandler fills a float32 numpy array via Xmap.export_section_numpy,
        # which is typed to the Xmap's scalar - a double Xmap would make pybind
        # fill a discarded float64 copy, leaving the displayed map all zeros.
        # Float (ftype32) throughout, to match the live-map display path: the
        # XmapHandler fills a float32 numpy array via Xmap.export_section_numpy,
        # which is typed to the Xmap's scalar - a double Xmap would make pybind
        # fill a discarded float64 copy, leaving the displayed map all zeros.
        from ..clipper_python import SFcalc_aniso_fft_float, Xmap_float, AtomShapeFn
        from ..clipper_python.data32 import HKL_data_F_phi_float
        from ..clipper_python.ext import scale_fcalc_to_fobs

        radiation = (AtomShapeFn.ELECTRON if str(self._radiation).lower() == 'electron'
                     else AtomShapeFn.XRAY)
        fcalc = HKL_data_F_phi_float(self.hklinfo)
        sfcalc = SFcalc_aniso_fft_float(2.5, self.FFT_RATE, 0.0, radiation=radiation)
        sfcalc(fcalc, atoms)

        # Scale Fcalc onto Fobs with an anisotropic Gaussian polished by an isotropic
        # spline (the same refinement-grade scaling Xtal_mgr_base uses) - far better
        # than an overall+B fit, so the residual no longer piles onto the strongest
        # scatterer (heavy atoms).
        scaled = HKL_data_F_phi_float(self.hklinfo)
        scale_fcalc_to_fobs(fcalc, self._fobs, scaled)

        _, cdata = fcalc.data            # cdata[:,0]=Fc, [:,1]=phi (radians)
        Fc = cdata[:, 0]
        phi = cdata[:, 1]
        hkls_i, sdata = scaled.data      # sdata[:,0] = scaled Fcalc (~Fo)
        scF = sdata[:, 0]
        fo = self._fo

        fit = self._meas & numpy.isfinite(scF) & (scF > 0) & (fo > 0)
        self._rwork = float(numpy.sum(numpy.abs(fo[fit] - scF[fit]))
                            / numpy.sum(fo[fit])) if fit.any() else 0.0

        # Put Fobs on the ABSOLUTE (electrons / A^3) scale via the same per-reflection
        # scale S = scaled_Fcalc/Fcalc (so Fo/S ~ Fc); fft_from of electron-scale
        # coefficients then yields e/A^3, for absolute-level contouring.
        with numpy.errstate(invalid='ignore', divide='ignore'):
            S = numpy.where(Fc > 1e-6, scF / Fc, numpy.nan)
            fobs_e = fo / S
        valid = self._meas & numpy.isfinite(fobs_e)
        fobs_e = numpy.where(valid, fobs_e, 0.0)
        Fc_f = numpy.nan_to_num(Fc)
        phi_f = numpy.nan_to_num(phi)
        twofofc = numpy.where(valid, 2.0 * fobs_e - Fc_f, 0.0)
        fofc = numpy.where(valid, fobs_e - Fc_f, 0.0)
        hkls_i = hkls_i.astype(numpy.int32)

        for m in self._maps.values():
            amp = fofc if m['is_difference'] else twofofc
            coeffs = HKL_data_F_phi_float(self.hklinfo)
            coeffs.set_data(hkls_i, numpy.stack([amp, phi_f], axis=1))
            xmap = Xmap_float(self.spacegroup, self.cell, self.grid_sampling)
            xmap.fft_from(coeffs)
            m['pending'] = xmap
