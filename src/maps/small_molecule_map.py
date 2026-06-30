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

    def __init__(self, hklinfo, cell, spacegroup, grid_sampling, scaffold,
                 fobs_amplitudes, structure=None, scaffold_to_model=None):
        '''
        Args:
            hklinfo, cell, spacegroup, grid_sampling: the crystal definition.
            scaffold: per-atom structure-factor inputs (see
                io.small_molecule.sfcalc_scaffold) - elements, Clipper-frame
                coords, occupancies, U_iso, U_aniso, labels.
            fobs_amplitudes: numpy array of Fo (sqrt of intensity), aligned to
                hklinfo reflection order; NaN where unmeasured.
            structure: optional live AtomicStructure (Clipper-frame coords) whose
                current coordinates are used on each recalc. If None, the
                scaffold's static coordinates are used (headless / static use).
            scaffold_to_model: optional int array mapping each scaffold atom to a
                model-atom index (or -1), for live coordinate override by label.
        '''
        self.hklinfo = hklinfo
        self.cell = cell
        self.spacegroup = spacegroup
        self.grid_sampling = grid_sampling
        self._scaffold = scaffold
        self._structure = structure
        self._scaffold_to_model = scaffold_to_model
        self._fo = numpy.asarray(fobs_amplitudes, numpy.double)
        self._meas = numpy.isfinite(self._fo)
        self._maps = {}            # name -> dict(is_difference, current, pending)
        self._thread = None
        self._ready = False
        self._rwork = 0.0
        self._n_reflections = int(self._meas.sum())

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
        from ..clipper_python import Xmap_double
        self._maps[name] = {
            'is_difference': bool(is_difference_map),
            'current': Xmap_double(self.spacegroup, self.cell, self.grid_sampling),
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
        from ..io.small_molecule import atom_list_from_scaffold
        coords = self._scaffold['coords']
        if self._structure is not None and self._scaffold_to_model is not None:
            live = self._structure.atoms.coords
            coords = coords.copy()
            sel = self._scaffold_to_model >= 0
            coords[sel] = live[self._scaffold_to_model[sel]]
        return atom_list_from_scaffold(self._scaffold, coords)

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
        from ..clipper_python import HKL, Resolution, SFcalc_aniso_fft_double, Xmap_double
        from ..clipper_python.data64 import HKL_data_F_phi_double

        fcalc = HKL_data_F_phi_double(self.hklinfo)
        sfcalc = SFcalc_aniso_fft_double(2.5, self.FFT_RATE, 0.0)
        sfcalc(fcalc, atoms)
        hkls, cdata = fcalc.data          # cdata[:,0]=Fc, cdata[:,1]=phi (radians)
        Fc = cdata[:, 0]
        phi = cdata[:, 1]
        fo = self._fo
        meas = self._meas & numpy.isfinite(Fc) & (Fc > 0)
        # Scale fit needs strictly-positive Fo (negative intensities clip to Fo=0).
        fit = meas & (fo > 0)

        # Scale Fc onto Fo (overall + isotropic B): ln(Fo/Fc) = ln k - (B/4) s^2.
        kFc = Fc.copy()
        if fit.sum() > 10:
            ss = numpy.array([HKL([int(h[0]), int(h[1]), int(h[2])]).invresolsq(self.cell)
                              for h in hkls])
            slope, intercept = numpy.polyfit(ss[fit], numpy.log(fo[fit] / Fc[fit]), 1)
            kFc = numpy.exp(intercept + slope * ss) * Fc
        elif fit.sum():
            k = numpy.sum(fo[fit] * Fc[fit]) / numpy.sum(Fc[fit] ** 2)
            kFc = k * Fc

        # R-work (over measured reflections), for the status line.
        self._rwork = float(numpy.sum(numpy.abs(fo[fit] - kFc[fit]))
                            / numpy.sum(fo[fit])) if fit.any() else 0.0

        fo_f = numpy.where(self._meas, fo, 0.0)
        kFc_f = numpy.nan_to_num(kFc)
        phi_f = numpy.nan_to_num(phi)
        twofofc = numpy.where(self._meas, 2.0 * fo_f - kFc_f, 0.0)
        fofc = numpy.where(self._meas, fo_f - kFc_f, 0.0)
        hkls_i = hkls.astype(numpy.int32)

        for m in self._maps.values():
            amp = fofc if m['is_difference'] else twofofc
            coeffs = HKL_data_F_phi_double(self.hklinfo)
            coeffs.set_data(hkls_i, numpy.stack([amp, phi_f], axis=1))
            xmap = Xmap_double(self.spacegroup, self.cell, self.grid_sampling)
            xmap.fft_from(coeffs)
            m['pending'] = xmap
