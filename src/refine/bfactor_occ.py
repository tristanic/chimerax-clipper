# Clipper plugin to UCSF ChimeraX
# Copyright (C) 2016-2019 Tristan Croll, University of Cambridge
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

"""
Python manager for background ML B-factor and occupancy refinement.

Refinement target
-----------------
The C++ refiner minimises a self-consistent crystallographic target at every
L-BFGS-B evaluation:

    T(U, occ)  = ½ · Σ_{h ∈ working}  (k·|Fc(h;U,occ)| − m_h·|Fo(h)|)²

where k is an isotropic least-squares scale recomputed each evaluation and
m_h = phi_fom.fom() (the sigma_A figure of merit, fixed for one L-BFGS run).

The driving density and gradient follow the Agarwal (1978) formulation:

    d(x)    = FFT{ k·(m_h·|Fo| − k·|Fc|)·exp(iφ_calc) }
    ∂T/∂p_j = −Σ_x  d(x)·∂ρ_j(x)/∂p_j

Because both T and d(x) are recomputed from the current Fcalc at each step,
the optimizer always sees a consistent gradient of a fixed objective —
avoiding the oscillation that results from using a pre-computed, stale map.

For crystallographic data pass an XmapSet as map_set (NXmapSet not yet
implemented).  The crystallographic data (fobs, phi_fom, usage flags) are
fetched from the XmapSet's live_xmap_mgr at each launch() call, so updated
sigma_A weights from the most recent map recalculation are always used.
"""

import numpy


class BFactorOccRefineManager:
    """
    Manages one or more macro-cycles of B-factor and/or occupancy refinement
    for a set of ChimeraX atoms against a Clipper density map.

    Parameters
    ----------
    session : chimerax.core.session.Session
    symmetry_manager : SymmetryManager
        Top-level Clipper model providing the limiting resolution and acting
        as the parent container for the structure and all associated maps.
    map_set : XmapSet or NXmapSet
        The map set containing the gradient density map. For crystallographic
        data this should be an XmapSet; NXmapSet support is not yet
        implemented.
    config : RefineConfig or None
        Refinement configuration. A default (B-factor only, 100 cycles) is
        used if not supplied.
    b_min_resolution_factor : float
        Multiplier for the resolution-dependent B-factor floor.
        config.b_min is raised to max(config.b_min, factor * d_min²) where
        d_min comes from symmetry_manager.resolution.limit.  The default of
        2.0 corresponds to the Gaussian sampling criterion (atom width ≥
        d_min / 2π).  For low-resolution data (> ~2.5 Å) a value of 3–4
        is typically more conservative and gives better behaviour.
    log_level : None | 'info' | 'debug'
        None    — no logging (default).
        'info'  — logs a one-line completion message via session.logger.info
                  when results are applied.
        'debug' — same message plus wall-clock time for the refinement cycle
                  (measured from the background thread launch to result
                  application).

    Restraints and occupancy
    ------------------------
    The engine generates no restraints of its own, and everything atom-specific
    is expressed in terms of the **input atom array** (positions in the atoms
    passed to ``launch`` / ``launch_realspace``) plus an explicit altloc — the
    C++ layer resolves these to its internal Clipper indices at launch.  Callers
    never need to know the Clipper Atom_list order.

      - **B-factor restraints** are passed to ``launch`` / ``launch_realspace``
        via the ``b_restraints`` keyword: a 6-tuple of equal-length sequences
        ``(atoms1, altlocs1, atoms2, altlocs2, sigmas, ks)`` — or a 7-tuple with
        a trailing ``alphas`` — where atomsN are indices into the launch atom
        array, altlocsN are altloc strings (``''`` = no altloc) selecting the
        conformer of each endpoint, sigmas/ks are the per-restraint Barron scale
        (c) and weight, and alphas (optional) are the per-restraint Barron shape.
        The restraint uses Barron's general robust loss: alpha=1 (the default
        when omitted) is Charbonnier/pseudo-Huber, alpha=-2 reproduces the
        Geman-McClure family, alpha=2 is harmonic.  The caller decides exactly
        which (atom, altloc) pairs to restrain.
      - **One-sided target restraints** (real-space only) are passed via the
        ``b_target_restraints`` keyword: a 5-tuple
        ``(atoms, altlocs, target_us, sigmas, ks)`` (or a 6-tuple with a trailing
        ``alphas``) pulling each atom's U toward a fixed ``target_u`` (Å²) — e.g.
        a neighbouring context atom's U, for fragment/context harmonisation.
      - **Occupancy refinement**: set ``config.refine_occ`` to a per-input-atom
        flag array (length = number of atoms passed to launch).  Occupancy
        groups are then auto-derived in C++: atoms in one contiguous covalent
        fragment that share an altloc get one shared occupancy, and the lettered
        altlocs of a fragment are constrained to sum to 1.

    Example::

        # refine B-factors of a ligand, restrained toward its protein context
        atoms = ...                       # Atoms passed to launch_realspace
        mgr.launch_realspace(
            atoms, context_atoms=ctx, target_volume=vol,
            b_target_restraints=(
                [0, 1, 2],                 # atom indices into `atoms`
                ['', '', ''],              # altlocs ('' = no altloc)
                [0.32, 0.30, 0.35],        # target U values (Å²)
                [0.063, 0.063, 0.063],     # sigmas
                [5.0, 5.0, 5.0]))          # weights
    """

    def __init__(self, session, symmetry_manager, map_set,
                 config=None,
                 b_min_resolution_factor=2.0,
                 n_threads=None,
                 log_level=None,
                 rfree_tolerance=0.02,
                 gap_tolerance=0.07):
        self._session      = session
        self._sym_mgr      = symmetry_manager
        self._map_set      = map_set
        self._b_min_factor = b_min_resolution_factor
        self._log_level         = log_level   # None | 'info' | 'debug'
        self._launch_time       = None        # set by timing wrapper in launch()
        self._rwork_before      = None        # R-work captured at launch() time
        self._rfree_before      = None        # R-free captured at launch() time
        # Overfitting guards.  R-work is deliberately NOT guarded — a small R-work
        # rise accompanied by an R-free fall is a *good* result.  The two guards
        # target the two failure modes:
        #   rfree_tolerance — max allowed per-run INCREASE in R-free (catches a
        #     single bad cycle; an absolute R-free ceiling is meaningless since the
        #     healthy value is wholly resolution-dependent).
        #   gap_tolerance   — absolute CEILING on the R-free − R-work gap (catches
        #     slow cumulative overfitting that a per-run check can't see).  Enforced
        #     only when this run also *widens* the gap, so an already-overfit model
        #     is never blocked from improving.  A fixed (not resolution-scaled)
        #     ceiling keeps this engine-only; clients may pass a resolution-aware value.
        # Set either to None to disable that guard.
        self._rfree_tolerance   = rfree_tolerance
        self._gap_tolerance     = gap_tolerance
        self._n_atoms_at_launch = 0           # passive length guard
        self._evaluate_only     = False       # set per-launch; skips write-back + bail-out
        self._last_rfactors     = None        # (rwork_before, rfree_before, rwork_after, rfree_after)
        self._cancel_results    = False       # set by proactive safety hooks
        self._watch_handles     = []          # (trigger_obj, handler) to clean up
        self._atoms             = None
        self._target_structure  = None        # single structure being refined; set in launch*()

        from ..clipper_python.ext import RefineConfig
        self._config = config if config is not None else RefineConfig()
        if n_threads is None:
            from ..util import available_cores
            n_threads = available_cores()
        self._config.n_threads = n_threads
        self._apply_resolution_b_min()
        self._thread_mgr = None   # created fresh on every launch()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @property
    def config(self):
        """The RefineConfig used by this manager."""
        return self._config

    def reset_config(self, config, n_threads=None):
        """Replace the refinement configuration. Safe to call between cycles.
        The new config takes effect on the next launch()."""
        self._config = config
        if n_threads is not None:
            self._config.n_threads = n_threads
        self._apply_resolution_b_min()

    def launch(self, atoms, ignore_hydrogens=False,
               auto_tolerances=True,
               delta_tolerance_factor=1.0,
               epsilon_tolerance_factor=1.0,
               b_restraints=None,
               evaluate_only=False):
        """
        Start a refinement cycle in the background.

        The crystallographic data (Fobs, sigma_A weights, usage flags) and
        atom data are fetched at call time and deep-copied into the background
        thread; the caller may modify them freely once this returns.

        The sigma_A weights (FOM values) are fixed for the duration of one
        L-BFGS-B run.  Re-calling launch() after map recalculation provides
        the updated sigma_A values, analogous to REFMAC macro-cycles.

        Parameters
        ----------
        atoms : chimerax.atomic.Atoms
            The atoms whose B-factors and/or occupancies will be refined.
        ignore_hydrogens : bool
            If True, hydrogen atoms are excluded from the refinement.
        auto_tolerances : bool
            If True (default), scale lbfgs_epsilon and lbfgs_delta to the
            current dataset (atom count and resolution) before launching.
            Applies only to this launch; self._config is not permanently
            modified.
        delta_tolerance_factor : float
            Multiplier on the auto-scaled lbfgs_delta.
            < 1 → tighter (more iterations, more precise).
            > 1 → looser (fewer iterations, exits sooner).
            Default 1.0.  Ignored when auto_tolerances=False.
        epsilon_tolerance_factor : float
            Multiplier on the auto-scaled lbfgs_epsilon.  Same semantics.
            Default 1.0.  Ignored when auto_tolerances=False.
        evaluate_only : bool
            If True, refine and compute the before/after R-factors (available
            afterwards via last_rfactors / initial_rfactors / compute_rfactors)
            but DO NOT write the result back to the model and DO NOT apply the
            overfitting bail-out.  Intended for scoring trial settings (e.g. an
            ISOLDE restraint-weight scan) without altering the structure, so
            trials are independent.  Default False == normal write-back behaviour.
        """
        if self.thread_running:
            raise RuntimeError(
                'A refinement cycle is already running. '
                'Wait for ready() before launching a new cycle.')

        # This manager refines a single model at a time.  Capture that structure
        # now (while the atoms are valid) for the model-closed safety check.
        _structures = atoms.unique_structures
        if len(_structures) != 1:
            raise ValueError(
                'B-factor refinement requires all atoms to belong to a single '
                f'structure (got {len(_structures)}).')
        self._target_structure = _structures[0]

        fobs, phi_fom, usage, f_bulk = self._get_hkl_data()
        self._atoms = atoms
        # Crystallographic refinement uses the dataset's own scattering factors,
        # which the map manager already knows (X-ray, or electron for micro-ED).
        # The same choice drives both the C++ table (use_electron_scattering) and
        # the per-atom identifiers below, so ionic species are validated against —
        # and logged for — the table actually used.
        _radiation = ('electron'
                      if str(getattr(self._map_set, '_radiation', 'xray')).lower() == 'electron'
                      else 'xray')
        self._config.use_electron_scattering = (_radiation == 'electron')
        # Explicit per-launch confirmation of the scattering factors in use. Not
        # routed through the de-duplicated ionic-species log, so a user-initiated
        # refinement always reports its radiation even when the map open already
        # announced the same ionic assignments.
        self._session.logger.info(
            '(CLIPPER) B-factor/occupancy refinement using %s scattering factors.'
            % ('electron' if self._config.use_electron_scattering else 'X-ray'))
        _validate_occ_groups(atoms, self._config)

        # Optionally compute data-aware convergence tolerances and apply them
        # temporarily to self._config before the thread manager is created.
        # self._config is restored afterwards so the user's explicit settings
        # are not permanently overwritten.
        _effective_eps    = None
        _effective_delta  = None
        _saved_tolerances = None
        if auto_tolerances:
            try:
                d_min = self._sym_mgr.resolution.limit
                _effective_eps, _effective_delta = self._compute_smart_tolerances(
                    len(atoms), d_min,
                    delta_factor=delta_tolerance_factor,
                    epsilon_factor=epsilon_tolerance_factor)
                _saved_tolerances = (self._config.lbfgs_epsilon,
                                     self._config.lbfgs_delta,
                                     self._config.lbfgs_past)
                self._config.lbfgs_epsilon = _effective_eps
                self._config.lbfgs_delta   = _effective_delta
                if self._config.lbfgs_past < 1:
                    self._config.lbfgs_past = 3
            except Exception:
                _effective_eps = _effective_delta = None   # log as config values

        # Always create a fresh thread manager from the current config so that
        # any changes made to self._config since the last launch() are reflected.
        from ..clipper_python.ext import BFactorOccRefinerThread
        self._thread_mgr = BFactorOccRefinerThread(self._config)

        # Restore the user's tolerance settings after the C++ copy is taken.
        if _saved_tolerances is not None:
            (self._config.lbfgs_epsilon,
             self._config.lbfgs_delta,
             self._config.lbfgs_past) = _saved_tolerances

        # Capture R-factors before refinement for before/after comparison.
        try:
            self._rwork_before = self._map_set.rwork
            self._rfree_before = self._map_set.rfree
        except Exception:
            self._rwork_before = None
            self._rfree_before = None

        # Passive guard: record atom count; bail in _apply_results() if it changes.
        self._n_atoms_at_launch = len(atoms)
        self._evaluate_only     = evaluate_only
        self._last_rfactors     = None
        self._cancel_results    = False

        # Proactive guards: cancel the running thread and log immediately if
        # atoms are deleted or the model is closed before results are applied.
        # Both handlers share a closure so each can clean up the other.
        from chimerax.atomic import get_triggers
        from chimerax.core.models import REMOVE_MODELS
        _cancelled  = [False]
        _handles    = []          # filled below; captured by closure

        def _do_cancel(reason):
            if _cancelled[0]:
                return
            _cancelled[0]        = True
            self._cancel_results = True
            self._thread_mgr.cancel()       # signal C++ thread to stop
            for t, h in _handles:
                try: t.remove_handler(h)
                except Exception: pass
            _handles.clear()
            self._watch_handles.clear()
            self._session.logger.warning(
                f'B-factor/occ refinement cancelled: {reason}')

        def _on_changes(trigger_name, changes):
            from chimerax.core.triggerset import DEREGISTER
            if changes.num_deleted_atoms() > 0:
                _do_cancel('atoms were deleted')
                return DEREGISTER

        def _on_remove_models(trigger_name, models):
            from chimerax.core.triggerset import DEREGISTER
            if self._target_structure in models:
                _do_cancel('target model was closed')
                return DEREGISTER

        _at  = get_triggers()
        _h1  = _at.add_handler('changes', _on_changes)
        _h2  = self._session.triggers.add_handler(REMOVE_MODELS, _on_remove_models)
        _handles.extend([(_at, _h1), (self._session.triggers, _h2)])
        self._watch_handles = list(_handles)   # keep reference for cleanup in _apply_results

        if self._log_level == 'debug':
            self._log_startup(atoms, fobs, ignore_hydrogens,
                              effective_epsilon=_effective_eps,
                              effective_delta=_effective_delta)

        # Wrap the start function to record the launch timestamp.
        # The closure captures self so _launch_time is set on the right object.
        def _timed_launch(*args):
            if self._log_level == 'debug':
                from time import perf_counter
                self._launch_time = perf_counter()
            self._thread_mgr.launch(*args)

        # The optional pairwise restraint spec is passed straight through as a
        # single tuple (or None); C++ unpacks it (see parse_pairwise_restraints).
        from ..delayed_reaction import delayed_reaction
        from ..scattering import ionic_scattering_names
        elements = ionic_scattering_names(atoms, radiation=_radiation)
        delayed_reaction(
            self._session.triggers, 'new frame',
            _timed_launch,
            [atoms.pointers, fobs, phi_fom, usage, ignore_hydrogens, f_bulk,
             b_restraints, elements],
            self._thread_mgr.ready,
            self._apply_results, []
        )

    @property
    def thread_running(self):
        return self._thread_mgr is not None and self._thread_mgr.thread_running

    @property
    def ready(self):
        return self._thread_mgr is not None and self._thread_mgr.ready()

    def initial_rfactors(self):
        """Standard (R_work, R_free) at the INPUT parameters of the most recent
        launch().  Crystallographic path only; blocks until the thread finishes,
        so valid once ready() is True.  Returns (-1.0, -1.0) if no launch() has
        created a thread yet."""
        if self._thread_mgr is None:
            return (-1.0, -1.0)
        return self._thread_mgr.initial_rfactors()

    def compute_rfactors(self):
        """Standard (R_work, R_free) at the REFINED parameters of the most recent
        launch().  Crystallographic path only; blocks until the thread finishes,
        so valid once ready() is True.  Returns (-1.0, -1.0) if no launch() has
        created a thread yet."""
        if self._thread_mgr is None:
            return (-1.0, -1.0)
        return self._thread_mgr.compute_rfactors()

    @property
    def last_rfactors(self):
        """The (rwork_before, rfree_before, rwork_after, rfree_after) tuple from
        the most recent run, or None until a run has produced R-factors (i.e. an
        evaluate_only run, or a normal run with an overfitting tolerance set)."""
        return self._last_rfactors

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _log_startup(self, atoms, fobs, ignore_hydrogens,
                     effective_epsilon=None, effective_delta=None):
        """Log a full configuration summary at debug level before launch.

        effective_epsilon / effective_delta: when auto_tolerances is active
        these are the values actually used this launch (vs the stored config).
        """
        cfg = self._config
        n_atoms_cx = len(atoms)

        # Resolution
        try:
            d_min = self._sym_mgr.resolution.limit
            res_str = f'{d_min:.2f} Å'
        except Exception:
            res_str = 'unknown'

        # Reflection count from the HKL_info attached to fobs
        try:
            n_refl = fobs.base_hkl_info.num_reflections
        except Exception:
            n_refl = None

        # refine_occ summary
        import numpy as _np
        ro = cfg.refine_occ          # numpy uint8 array or empty
        if len(ro) == 0:
            occ_str = 'all fixed'
        else:
            n_refined = int(_np.count_nonzero(ro))
            occ_str = f'{n_refined} / {len(ro)} atoms flagged'

        rw = self._rwork_before
        rf = self._rfree_before
        r_str = (f'Rw {rw:.4f}  Rf {rf:.4f}'
                 if rw is not None else 'not yet computed')

        lines = [
            'B-factor/occ refinement starting:',
            f'  R-factors (before): {r_str}',
            f'  Atoms:              {n_atoms_cx}'
            + (' (hydrogens excluded)' if ignore_hydrogens else ''),
            f'  Resolution:       {res_str}',
            f'  Reflections:      {n_refl if n_refl is not None else "unknown"} total',
            '  RefineConfig:',
            f'    refine_b:         {cfg.refine_b}',
            f'    refine_occ:       {occ_str}',
            f'    b_min:            {cfg.b_min:.3g} Å²',
            f'    b_max:            {cfg.b_max:.3g} Å²',
            f'    max_cycles:       {cfg.max_cycles}',
            f'    n_threads:        {cfg.n_threads}',
            f'    lbfgs_epsilon:    '
            + (f'{effective_epsilon:.2e} (auto)' if effective_epsilon is not None
               else f'{cfg.lbfgs_epsilon:.3g}'),
            f'    lbfgs_past:       {cfg.lbfgs_past}',
            f'    lbfgs_delta:      '
            + (f'{effective_delta:.2e} (auto)' if effective_delta is not None
               else f'{cfg.lbfgs_delta:.3g}'),
            f'    use_curvature:    {cfg.use_curvature}',
        ]
        self._session.logger.info('\n'.join(lines))

    def launch_realspace(self, refined_atoms, context_atoms=None,
                         target_volume=None, ignore_hydrogens=False,
                         padding=6.0, taper_width=3.0,
                         auto_tolerances=True,
                         delta_tolerance_factor=1.0,
                         epsilon_tolerance_factor=1.0,
                         resolution=None,
                         whole_map=False,
                         b_restraints=None,
                         b_target_restraints=None,
                         use_electron_scattering=True):
        """
        Refine B-factors against a fixed real-space target density.

        Calls ``isolate_and_cover_selection`` to ensure the map data is loaded
        for the region of interest, extracts a P1 Xmap subregion, then runs the
        background L-BFGS-B thread against it.  A per-iteration least-squares
        scale factor is applied so the result is independent of the map's
        absolute normalisation.

        Parameters
        ----------
        refined_atoms : chimerax.atomic.Atoms
            Atoms whose B-factors/occupancies are optimised.
        context_atoms : chimerax.atomic.Atoms or None
            Additional atoms contributing to ρ_calc but not refined.
        target_volume : chimerax.map.Volume
            The ChimeraX Volume to refine against.
        ignore_hydrogens : bool
            Exclude H from EDcalc and refinement.
        padding : float
            Extra Å on all sides when extracting the map subregion.
        taper_width : float
            Width of the cosine taper zone at the box edges (Å).
        auto_tolerances : bool
            If True (default), scale lbfgs_epsilon and lbfgs_delta to the
            current dataset.
        delta_tolerance_factor : float
            Multiplier on the auto-scaled lbfgs_delta.  < 1 tightens
            convergence; > 1 loosens it.  Ignored when auto_tolerances=False.
        epsilon_tolerance_factor : float
            Multiplier on the auto-scaled lbfgs_epsilon.  Same semantics.
            Ignored when auto_tolerances=False.
        resolution : float or None
            Nominal resolution of the target map in Å.  Used to scale the
            lbfgs_delta tolerance when auto_tolerances=True.  If None (default),
            the effective resolution is estimated as 3 × max(voxel_step) — a
            conservative Nyquist estimate that works for typical cryo-EM
            oversampling factors.  Pass the known value when refining against
            a crystallographic map or a cryo-EM map with a specified resolution.
        whole_map : bool
            If True, refine against the full extent of target_volume rather than
            an atom-bounded subregion.  ``isolate_and_cover_selection`` is not
            called.  Only valid for non-crystallographic (e.g. cryo-EM) maps;
            raises RuntimeError for XmapHandlerBase subclasses.
        use_electron_scattering : bool
            If True (default), use electron scattering factors for the calculated
            density and Agarwal gradients.  Real-space targets are most often
            cryo-EM potential maps, for which electron factors are correct; pass
            False to refine against an X-ray-derived real-space map.
        """
        if self.thread_running:
            raise RuntimeError(
                'A refinement cycle is already running. '
                'Wait for ready() before launching a new cycle.')
        # Real-space maps are usually cryo-EM electrostatic-potential maps, so
        # electron scattering factors are the sensible default (overridable).
        self._config.use_electron_scattering = bool(use_electron_scattering)
        self._session.logger.info(
            '(CLIPPER) Real-space B-factor/occupancy refinement using %s scattering factors.'
            % ('electron' if self._config.use_electron_scattering else 'X-ray'))
        if target_volume is None:
            raise ValueError('target_volume must be provided for real-space launch.')
        if whole_map:
            from .maps.map_handler_base import XmapHandlerBase
            if isinstance(target_volume, XmapHandlerBase):
                raise RuntimeError(
                    'whole_map=True is not supported for crystallographic maps; '
                    'the live map covers only the current spotlight box.')

        # Identify the full atom set needed to define the map box.
        try:
            from chimerax.atomic import concatenate
            _all_atoms = (concatenate([refined_atoms, context_atoms])
                          if context_atoms is not None else refined_atoms)
        except Exception:
            self._session.logger.warning('Failed to concatenate refined and context atom sets ')
            _all_atoms = refined_atoms

        # This manager refines a single model against a single map.  Enforce
        # that all atoms (refined + context) belong to one structure, and that
        # both that structure and the target map live under this manager's
        # SymmetryManager.  Capture the structure for the model-closed check.
        _structures = _all_atoms.unique_structures
        if len(_structures) != 1:
            raise ValueError(
                'Real-space refinement requires all atoms (refined + context) '
                f'to belong to a single structure (got {len(_structures)}).')
        self._target_structure = _structures[0]
        if self._target_structure is not self._sym_mgr.structure:
            raise ValueError(
                'The atoms to refine are not part of this Clipper session.')
        _owner = target_volume
        while _owner is not None and _owner is not self._sym_mgr:
            _owner = _owner.parent
        if _owner is not self._sym_mgr:
            raise ValueError(
                'The target map and the model must belong to the same Clipper '
                'session (SymmetryManager).')

        # Optionally scale L-BFGS-B tolerances to the dataset.
        # d_min comes from the caller when known; otherwise estimated as
        # 3 × max(voxel_step) — conservative Nyquist for cryo-EM oversampling.
        _saved_tolerances_rs = None
        if auto_tolerances:
            try:
                import math
                if resolution is not None:
                    d_min_eff = float(resolution)
                else:
                    try:
                        d_min_eff = self._sym_mgr.resolution.limit
                    except Exception:
                        d_min_eff = 3.0 * max(target_volume.data.step)
                _eff_eps, _eff_delta = self._compute_smart_tolerances(
                    len(refined_atoms), d_min_eff,
                    delta_factor=delta_tolerance_factor,
                    epsilon_factor=epsilon_tolerance_factor)
                _saved_tolerances_rs = (self._config.lbfgs_epsilon,
                                        self._config.lbfgs_delta,
                                        self._config.lbfgs_past)
                self._config.lbfgs_epsilon = _eff_eps
                self._config.lbfgs_delta   = _eff_delta
                if self._config.lbfgs_past < 1:
                    self._config.lbfgs_past = 3
            except Exception:
                _saved_tolerances_rs = None

        def _do_launch():
            if whole_map:
                target_xmap, target_origin = _build_p1_target_xmap(
                    target_volume, atoms=None, taper_width=taper_width)
            else:
                target_xmap, target_origin = _build_p1_target_xmap(
                    target_volume, _all_atoms, padding=padding, taper_width=taper_width)

            self._atoms             = refined_atoms
            self._n_atoms_at_launch = len(refined_atoms)
            # _apply_results is shared with launch(); reset so a prior
            # crystallographic evaluate_only=True can't leak into this run.
            self._evaluate_only     = False
            self._last_rfactors     = None
            self._cancel_results    = False
            context_ptrs = (context_atoms.pointers.tolist()
                            if context_atoms is not None else [])

            # Safety hooks — same pattern as crystallographic launch.
            from chimerax.atomic import get_triggers
            from chimerax.core.models import REMOVE_MODELS
            _cancelled = [False]
            _handles   = []

            def _do_cancel(reason):
                if _cancelled[0]:
                    return
                _cancelled[0]        = True
                self._cancel_results = True
                self._thread_mgr.cancel()
                for t, h in _handles:
                    try: t.remove_handler(h)
                    except Exception: pass
                _handles.clear()
                self._watch_handles.clear()
                self._session.logger.warning(
                    f'B-factor/occ real-space refinement cancelled: {reason}')

            def _on_changes(trigger_name, changes):
                from chimerax.core.triggerset import DEREGISTER
                if changes.num_deleted_atoms() > 0:
                    _do_cancel('atoms were deleted')
                    return DEREGISTER

            def _on_remove_models(trigger_name, models):
                from chimerax.core.triggerset import DEREGISTER
                if self._target_structure in models:
                    _do_cancel('target model was closed')
                    return DEREGISTER

            _at = get_triggers()
            _h1 = _at.add_handler('changes', _on_changes)
            _h2 = self._session.triggers.add_handler(REMOVE_MODELS, _on_remove_models)
            _handles.extend([(_at, _h1), (self._session.triggers, _h2)])
            self._watch_handles = list(_handles)

            from ..clipper_python.ext import BFactorOccRefinerThread
            self._thread_mgr = BFactorOccRefinerThread(self._config)
            # Config has been deep-copied into the C++ thread manager; safe to restore.
            if _saved_tolerances_rs is not None:
                (self._config.lbfgs_epsilon,
                 self._config.lbfgs_delta,
                 self._config.lbfgs_past) = _saved_tolerances_rs

            # Unpack the optional restraint arrays (ChimeraX-index + altloc),
            # captured by the closure below.
            # Restraint specs pass straight through as single tuples (or None);
            # C++ unpacks them (see parse_pairwise_restraints / parse_target_restraints).
            from ..scattering import ionic_scattering_names
            _radiation = 'electron' if use_electron_scattering else 'xray'
            refined_elements = ionic_scattering_names(refined_atoms, radiation=_radiation)
            context_elements = (ionic_scattering_names(context_atoms, radiation=_radiation)
                                if context_atoms is not None else [])

            def _timed_launch(refined_ptrs, ctx_ptrs, xmap, origin, ignore_h):
                if self._log_level == 'debug':
                    from time import perf_counter
                    self._launch_time = perf_counter()
                self._thread_mgr.launch_realspace(
                    refined_ptrs, ctx_ptrs, xmap, origin, ignore_h,
                    restraints=b_restraints,
                    target_restraints=b_target_restraints,
                    refined_elements=refined_elements,
                    context_elements=context_elements)

            from ..delayed_reaction import delayed_reaction
            delayed_reaction(
                self._session.triggers, 'new frame',
                _timed_launch,
                [refined_atoms.pointers, context_ptrs, target_xmap, target_origin,
                 ignore_hydrogens],
                self._thread_mgr.ready,
                self._apply_results, []
            )

        if whole_map:
            # Whole-map mode: the volume data is already fully loaded; skip
            # the spotlight repositioning step entirely.
            _do_launch()
        else:
            # isolate_and_cover_selection() fills the map data synchronously;
            # only the subsequent contour rebuild for display is async.  The P1
            # target can therefore be extracted immediately after it returns.
            self._sym_mgr.isolate_and_cover_selection(
                _all_atoms,
                hide_surrounds=True,
                focus=False,
                extra_padding=padding)
            _do_launch()

    def _get_hkl_data(self):
        """
        Return (fobs, phi_fom, usage, f_bulk) from the live map manager.

        f_bulk is the bulk-solvent contribution (F_total = F_atoms + f_bulk),
        held fixed during refinement and added to the calculated F_atoms each
        cycle.  Fetched fresh on every launch() call so that the sigma_A weights
        (phi_fom.fom()) and bulk-solvent model are current.

        The live map fits its Fcalc->Fobs scale (and hence f_bulk AND the sigma_A
        weights, since both derive from the scaled Fcalc) on a small random
        reflection subset for speed, so those inputs jitter by ~1% run-to-run --
        which propagates into the refined B-factors.  Refinement should be
        reproducible, so if the map is not already deterministic we run ONE
        fully-deterministic recalculation (scale fit over all reflections), pull
        the resulting f_bulk / weights / fobs, then restore the map's original
        (fast, stochastic) setting.  The pulled values are copies, so the restore
        does not disturb them.

        Raises NotImplementedError for NXmapSet (no crystallographic data).
        Raises RuntimeError if no live map manager is available.
        """
        from ..maps.nxmapset import NXmapSet
        if isinstance(self._map_set, NXmapSet):
            raise NotImplementedError(
                'B-factor refinement against non-crystallographic (NXmap) data '
                'is not yet implemented.')

        xm = self._map_set.live_xmap_mgr
        if xm is None:
            raise RuntimeError(
                'No live map manager is available on the XmapSet. '
                'B-factor refinement requires a live (Xtal_thread_mgr) map set.')

        was_deterministic = xm.deterministic_scaling
        if not was_deterministic:
            xm.deterministic_scaling = True
        try:
            # The live map fits the bulk scale only when flagged (a change to
            # B-factors/occupancies/coords does so via XmapSet._model_changed_cb);
            # force a fresh fit here so f_bulk matches the current model, and run one
            # blocking recalc so it (and the sigma_A weights) are ready to read. In
            # deterministic mode this fit uses all reflections, cutting the ~1% run-to-
            # run scale jitter (a residual warm-start dependence of ~0.01% remains).
            xm.bulk_solvent_optimization_needed()
            self._recompute_maps_blocking(xm)
            return xm.f_obs, xm.weights, xm.usage_flags, xm.f_bulk
        finally:
            if not was_deterministic:
                # Restore the fast stochastic scaling for interactive live maps.
                xm.deterministic_scaling = False

    def _recompute_maps_blocking(self, xm):
        '''Drive the live manager's recalc/ready/apply cycle synchronously (the
        XmapSet's normal path defers it to a "new frame" trigger, which never fires in
        the middle of this call).'''
        import time
        from ..scattering import ionic_scattering_names
        atoms = self._map_set.structure.atoms
        elements = ionic_scattering_names(
            atoms, radiation=getattr(self._map_set, '_radiation', 'xray'))
        for _ in range(60000):
            if not xm.thread_running:
                break
            time.sleep(0.005)
        xm.recalculate_all_maps(atoms.pointers, elements)
        for _ in range(60000):
            if xm.ready():
                break
            time.sleep(0.005)
        xm.apply_new_maps()

    def _compute_smart_tolerances(self, n_atoms, d_min,
                                   delta_factor=1.0, epsilon_factor=1.0):
        """
        Return (lbfgs_epsilon, lbfgs_delta) scaled to the dataset.

        lbfgs_epsilon scales with √N_atoms: the gradient infinity-norm is
        accumulated over all atoms so its natural magnitude grows with
        atom count.  Reference calibration is 1 000 atoms.

        lbfgs_delta scales with d_min²: B-factor precision from X-ray data
        is bounded by resolution (σ(B) ∝ d_min²), so converging more
        tightly than the data supports wastes iterations.

        delta_factor / epsilon_factor multiply the respective computed values,
        allowing the caller to tighten (< 1) or loosen (> 1) convergence
        relative to the data-aware baseline.

        Both outputs are clamped to prevent pathological inputs producing
        absurd tolerances.
        """
        import math
        eps   = 1e-5 * math.sqrt(max(1, n_atoms) / 1000.0) * epsilon_factor
        delta = 1e-4 * (d_min / 1.5) ** 2 * delta_factor
        return max(1e-7, min(1e-3, eps)), max(1e-6, min(1e-2, delta))

    def _apply_resolution_b_min(self):
        """
        Raise config.b_min to the resolution-dependent floor.

        Floor = b_min_resolution_factor × d_min²  (Å²), where d_min comes
        from symmetry_manager.resolution.limit.  This prevents L-BFGS-B from
        driving B-factors below values that are physically implausible at the
        data resolution, which would also corrupt the per-iteration Fcalc
        computation (atoms too sharp → FFT aliasing on the grid).
        """
        try:
            d_min = self._sym_mgr.resolution.limit
        except Exception:
            return
        b_min_res = self._b_min_factor * d_min * d_min
        if b_min_res > self._config.b_min:
            self._config.b_min = b_min_res

    def _apply_results(self):
        """
        Called in the main thread once the background refinement is done.
        Delegates all write-back to the C++ apply_to_atoms() which handles
        per-altloc B-factor and occupancy assignment correctly.
        """
        # Deregister any live change-watchers — they've done their job.
        for t, h in self._watch_handles:
            try: t.remove_handler(h)
            except Exception: pass
        self._watch_handles.clear()

        if self._atoms is None:
            return

        # Proactive cancel: already logged by _do_cancel.
        if self._cancel_results:
            return

        # Passive guard: any atom deletion shrinks the Atoms collection.
        try:
            current_n = len(self._atoms)
        except Exception:
            return
        if current_n == 0:
            return
        if current_n != self._n_atoms_at_launch:
            self._session.logger.warning(
                f'B-factor/occ refinement: atom count changed from '
                f'{self._n_atoms_at_launch} to {current_n} during refinement; '
                'results discarded.')
            return

        mgr = self._thread_mgr

        # Measure before/after R-factors when needed — either for the overfitting
        # bail-out (a tolerance is set) or for evaluate_only scoring.  Compare with
        # the SAME metric on both sides: the refiner's own initial_rfactors()
        # (standard R at the input parameters) vs compute_rfactors() (standard R at
        # the refined parameters).  Using the refiner for both makes the delta the
        # true change, immune to the small scale-model offset between the refiner and
        # the live map manager — and avoids the old FOM-weighted-vs-standard mismatch
        # that spuriously tripped the guard on incomplete models.  Skipped entirely
        # (no extra EDcalc+FFT) when neither a tolerance nor evaluate_only applies,
        # preserving the original fast path.
        self._last_rfactors = None
        rwork_before = rfree_before = rwork_new = rfree_new = None
        if (self._evaluate_only
                or self._rfree_tolerance is not None
                or self._gap_tolerance is not None):
            try:
                rwork_before, rfree_before = mgr.initial_rfactors()
                rwork_new, rfree_new = mgr.compute_rfactors()
                self._last_rfactors = (rwork_before, rfree_before,
                                       rwork_new, rfree_new)
            except Exception:
                rwork_before = rfree_before = rwork_new = rfree_new = None

        # Overfitting bail-out — normal mode only.  Discard the result if it overfits
        # (R-free rises, or the R-free−R-work gap widens) beyond tolerance.  R-work is
        # NOT guarded — a small R-work rise with R-free falling is a desirable outcome.
        # In evaluate_only mode we deliberately skip this so EVERY trial is measured,
        # including overfit ones the caller wants to see and reject itself.
        if (not self._evaluate_only
                and (self._rfree_tolerance is not None or self._gap_tolerance is not None)
                and rfree_before is not None and rfree_new is not None
                and rfree_before > 0.0 and rfree_new > 0.0):
            delta_rfree = rfree_new - rfree_before
            gap_before  = rfree_before - rwork_before
            gap_after   = rfree_new - rwork_new
            # R-free: reject a per-run increase beyond tolerance.
            rfree_bad = (self._rfree_tolerance is not None
                         and delta_rfree > self._rfree_tolerance)
            # Gap: reject if it exceeds the absolute ceiling AND this run widened it
            # (so an already-overfit model can still be improved).
            gap_bad = (self._gap_tolerance is not None
                       and gap_after > self._gap_tolerance
                       and gap_after > gap_before)
            if rfree_bad or gap_bad:
                self._session.logger.warning(
                    f'B-factor refinement reverted (overfitting guard): '
                    f'ΔR-free {delta_rfree:+.4f} '
                    f'(tol {self._rfree_tolerance}); '
                    f'R-free−R-work gap {gap_before:.4f}→{gap_after:.4f} '
                    f'(ceiling {self._gap_tolerance}).  '
                    f'Rw {rwork_before:.4f}→{rwork_new:.4f}  '
                    f'Rf {rfree_before:.4f}→{rfree_new:.4f}')
                return   # discard without touching ChimeraX atoms

        # evaluate_only: results are recorded in self._last_rfactors; never write
        # back to the model, so trials are independent.
        if self._evaluate_only:
            return

        mgr.apply_to_atoms()   # C++ writes back per (Atom*, altloc) pair

        if self._log_level:
            n_t  = self._config.n_threads
            msg  = (f'B-factor/occ refinement complete '
                    f'({mgr.n_atoms} atoms, '
                    f'{n_t} thread{"s" if n_t != 1 else ""})')
            if self._rwork_before is not None:
                msg += (f' — Rw {self._rwork_before:.4f}'
                        f'  Rf {self._rfree_before:.4f} (before)')
            if self._log_level == 'debug' and self._launch_time is not None:
                from time import perf_counter
                elapsed = perf_counter() - self._launch_time
                msg += f', {elapsed:.2f} s'
            self._session.logger.info(msg)

            # Register a one-shot handler to log R-factors after the maps are
            # recalculated with the new B-factors.  Only meaningful for the
            # crystallographic path where R-factors were captured at launch.
            if self._rwork_before is not None:
                rw_before = self._rwork_before
                rf_before = self._rfree_before
                map_set   = self._map_set
                session   = self._session

                def _log_r_after(trigger_name, _data):
                    from chimerax.core.triggerset import DEREGISTER
                    try:
                        rw_after = map_set.rwork
                        rf_after = map_set.rfree
                        if rw_after is not None:
                            session.logger.info(
                                f'B-factor/occ R-factors after recalculation: '
                                f'Rw {rw_before:.4f} → {rw_after:.4f}  '
                                f'Rf {rf_before:.4f} → {rf_after:.4f}')
                    except Exception:
                        pass
                    return DEREGISTER

                map_set.triggers.add_handler('maps recalculated', _log_r_after)


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _validate_occ_groups(atoms, config):
    """
    Validate the occupancy-refinement configuration against the atom array.

    Occupancy constraint groups are auto-derived in C++ from covalent-fragment +
    altloc topology, so the only caller-supplied occupancy input is the
    `refine_occ` flag array.  It must be either empty (no occupancy refinement)
    or exactly one flag per atom in `atoms` (ChimeraX-indexed).

    Raises ValueError on a length mismatch.
    """
    refine_occ = config.refine_occ          # numpy uint8 array, or empty
    if len(refine_occ) == 0:
        return
    n_atoms = len(atoms)
    if len(refine_occ) != n_atoms:
        raise ValueError(
            f'refine_occ array length ({len(refine_occ)}) does not match '
            f'atom count ({n_atoms}).')


def _scene_coords_to_ijk(coords, volume):
    """
    Convert atom SCENE coordinates (Å) to volume grid indices (i,j,k).
    Accounts for a non-identity volume.position by going
    scene → local-data frame → grid index.
    """
    local = volume.position.inverse().transform_points(coords)
    return volume.data.xyz_to_ijk_transform.transform_points(local)


def _atom_ijk_bounds(atoms, volume, padding):
    """Return (ijk_min, ijk_max) integer voxel bounds covering atoms + padding."""
    import numpy
    step     = numpy.array(volume.data.step)  # (sx, sy, sz) in Å per voxel
    ijk      = _scene_coords_to_ijk(atoms.coords, volume)
    pad_vox  = padding / step
    lo = numpy.floor(ijk.min(axis=0) - pad_vox).astype(int)
    hi = numpy.ceil (ijk.max(axis=0) + pad_vox).astype(int)
    grid_max = numpy.array(volume.data.size) - 1   # (nx-1, ny-1, nz-1)
    return numpy.maximum(lo, 0), numpy.minimum(hi, grid_max)


def _apply_edge_taper(data, step, taper_width):
    """
    Apply a cosine taper over the outer taper_width Å on all six faces.

    data has shape (nz, ny, nx).  ChimeraX step = (sx, sy, sz) where
    sx corresponds to the nx (last) axis, sz to nz (first axis).
    The loop therefore matches numpy axis 0→z, 1→y, 2→x.
    """
    import numpy
    tapered = data.copy()
    nz, ny, nx = data.shape
    # Pair each numpy axis with its extent and physical step:
    #   numpy axis 0 = z-direction  (nz voxels, step[2]=sz)
    #   numpy axis 1 = y-direction  (ny voxels, step[1]=sy)
    #   numpy axis 2 = x-direction  (nx voxels, step[0]=sx)
    for np_axis, (n, ds) in enumerate([(nz, step[2]), (ny, step[1]), (nx, step[0])]):
        t_vox = max(1, int(round(taper_width / ds)))
        ramp  = 0.5 * (1.0 - numpy.cos(numpy.pi * numpy.arange(t_vox) / t_vox))
        shape = [1, 1, 1]; shape[np_axis] = t_vox
        # low-index face
        sl = [slice(None)] * 3; sl[np_axis] = slice(0, t_vox)
        tapered[tuple(sl)] *= ramp.reshape(shape)
        # high-index face
        sl[np_axis] = slice(n - t_vox, n)
        tapered[tuple(sl)] *= ramp[::-1].reshape(shape)
    return tapered


def _debug_show_p1_target(session, xmap, target_origin, name='P1 debug target'):
    """
    Read a Clipper P1 Xmap back into a numpy array and open it as a ChimeraX
    Volume placed at target_origin (scene coordinates).

    Works for both the extracted target density and the EDcalc diagnostic density::

        target_xmap, origin = _build_p1_target_xmap(vol, atoms)
        _debug_show_p1_target(session, target_xmap, origin, 'target')

        from chimerax.clipper.clipper_python.ext import BFactorOccRefinerThread, RefineConfig
        mgr = BFactorOccRefinerThread(RefineConfig())
        edcalc_xmap = mgr.compute_realspace_density(
            atoms.pointers, [], target_xmap, origin)
        _debug_show_p1_target(session, edcalc_xmap, origin, 'EDcalc')

    If both volumes overlay the source map, the extraction and shift paths are
    consistent.  Any spatial offset between them indicates a frame mismatch in
    the -target_origin shift applied inside compute_realspace_density.
    """
    import numpy
    from chimerax.map_data import ArrayGridData
    from chimerax.map import volume_from_grid_data

    # Use the same pattern as XmapHandler_Base._generate_data_array:
    # export_section_numpy fills shape (nu, nv, nw) in C-order (w fastest).
    # data.transpose() converts to (nw, nv, nu) = (nz, ny, nx) for ChimeraX.
    dim  = numpy.array(xmap.grid_sampling.dim)    # [nu, nv, nw]
    step = xmap.cell.dim / dim                    # [sx, sy, sz]
    data = numpy.empty(dim, dtype='float32')      # shape (nu, nv, nw)
    origin_xyz = numpy.array([float(target_origin[0]),
                               float(target_origin[1]),
                               float(target_origin[2])])
    gd = ArrayGridData(data.transpose(), origin=origin_xyz,
                       step=step, cell_angles=xmap.cell.angles_deg,
                       name=name)
    xmap.export_section_numpy(xmap.first.coord, data, 1)
    v = volume_from_grid_data(gd, session)
    v.set_display_style('mesh')
    return v



def _build_p1_target_xmap(volume, atoms=None, padding=6.0, taper_width=3.0):
    """
    Extract the density subregion covering `atoms` + `padding` Å, apply a
    cosine taper at the box edges, and copy the result into a Clipper P1 Xmap.

    When `atoms` is None the full volume extent is used (whole-map mode).

    Returns (xmap, target_origin) where target_origin is a Clipper Coord_orth
    in SCENE coordinates giving the P1 cell origin.  All transforms account
    for a non-identity volume.position.
    """
    import numpy
    from ..clipper_python import (
        Xmap_float as _Xmap, Spacegroup, Spgr_descr,
        Cell, Cell_descr, Grid_sampling, Coord_orth)

    if atoms is None:
        ijk_min = numpy.zeros(3, dtype=int)
        ijk_max = numpy.array(volume.data.size) - 1
    else:
        ijk_min, ijk_max = _atom_ijk_bounds(atoms, volume, padding)

    # Extract subregion at full resolution.
    region = (tuple(ijk_min.tolist()), tuple(ijk_max.tolist()), (1, 1, 1))
    data   = volume.region_matrix(region).astype('float32')
    data   = _apply_edge_taper(data, volume.data.step, taper_width)

    # data.shape = (nz, ny, nx); step = (sx, sy, sz).
    # Cell axes: a=x (nx*sx), b=y (ny*sy), c=z (nz*sz) — all orthogonal.
    nz, ny, nx = data.shape
    sx, sy, sz  = volume.data.step

    alpha, beta, gamma = volume.data.cell_angles   # degrees; (90,90,90) for orthogonal
    cell = Cell(Cell_descr(nx * sx, ny * sy, nz * sz, alpha, beta, gamma))
    sg   = Spacegroup(Spgr_descr('P 1'))
    grid = Grid_sampling(nx, ny, nz)   # (nu, nv, nw) = (nx, ny, nz)
    xmap = _Xmap(sg, cell, grid)

    # Populate the Xmap using the same convention as numpy_import_core_:
    # export_section_numpy / import_section_numpy expect shape (nu, nv, nw)
    # in C-order with w fastest — the "Clipper-native" layout.
    # region_matrix returns (nz, ny, nx) = (nw, nv, nu) in ChimeraX C-order.
    # Transposing gives (nx, ny, nz) = (nu, nv, nw) — exactly what
    # import_section_numpy requires.  ascontiguousarray forces a C-contiguous
    # copy because data.transpose() is a non-contiguous view, and
    # numpy_import_core_ reads sequentially through *vptr++ without strides.
    clipper_data = numpy.ascontiguousarray(data.transpose())
    xmap.import_section_numpy(xmap.first.coord, clipper_data)

    # P1 cell origin in SCENE coordinates.
    # Chain: ijk_min → local-data xyz → scene xyz.
    local_origin = volume.data.ijk_to_xyz_transform.transform_points(
        ijk_min.astype(float).reshape(1, 3))[0]
    scene_origin = volume.position.transform_points(local_origin.reshape(1, 3))[0]
    target_origin = Coord_orth(*scene_origin.tolist())

    return xmap, target_origin
