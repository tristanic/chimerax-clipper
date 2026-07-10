# Clipper plugin to UCSF ChimeraX
# Copyright (C) 2016-2019 Tristan Croll, University of Cambridge

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

'''
Numeric (torch-free) cores for the differentiable crystallographic targets.

Each ``*TargetState`` object holds the fixed experimental reference data once,
and exposes ``value_and_gradient(coords)`` returning ``(L, dL/dp)`` as plain
numpy. These are the objects the PyTorch adapters in :mod:`.targets` wrap; they
also stand alone for finite-difference checking (:mod:`._checks`) or use with a
non-torch optimiser.

Coordinates are always in the **orthogonal (Cartesian) Angstrom** frame — the
same frame ``Grad_orth`` gradients are returned in — so there is no fractional /
grid unit mixing at the module boundary.

Tier 0 (:class:`MapTargetState`) is implemented here: an MDFF-style target that
scores atoms against a *fixed* experimental map (e.g. a precomputed 2mFo-DFc
map). Because the target field is fixed, only coordinate gradients are
meaningful for it. The model-density (Tier 1) and reciprocal-space (Tier 2)
targets — which additionally yield Occ / Uiso / anisotropic-U gradients via
``AtomShapeFn.rho_grad`` — are added on top of this.
'''

import numpy


#: Canonical order of the eleven refinable per-atom parameters.
_CANONICAL_PARAMS = ('X', 'Y', 'Z', 'Uiso', 'Occ',
                     'U11', 'U22', 'U33', 'U12', 'U13', 'U23')


class MapTargetState:
    '''
    MDFF-style differentiable target against a fixed experimental map.

    The loss is ``L = -sum_j w_j * rho(x_j)`` — the negative weighted sum of the
    (cubic-interpolated) experimental density sampled at each atom position — so
    that descending ``L`` pulls atoms up the density gradient. The gradient
    ``dL/dx_j = -w_j * grad(rho)(x_j)`` is the analytic gradient of the cubic
    interpolant (Clipper ``Xmap.interp_cubic_grad_frac``), so it agrees with a
    finite-difference of ``interp_cubic`` to interpolation precision.

    The map is crystallographic, so density is evaluated in the full periodic
    context automatically (symmetry + cell wrapping are handled by ``Xmap``).

    Args:
        * xmap:
            - a Clipper ``Xmap`` (e.g. ``XmapHandler_*.xmap``) holding the fixed
              target density.
        * weights:
            - optional ``(N,)`` array of per-atom weights (default: all ones).
              Held fixed; a common choice is the atomic number or occupancy.
    '''

    #: Parameters this target produces gradients for. A fixed map has no
    #: dependence on the model's own Occ/Uiso, so only coordinates apply.
    param_names = ('X', 'Y', 'Z')

    def __init__(self, xmap, weights=None):
        self._xmap = xmap
        self._cell = xmap.cell
        self._weights = None if weights is None else numpy.asarray(weights, numpy.double)

    def value_and_gradient(self, coords):
        '''
        Args:
            * coords: ``(N, 3)`` orthogonal (Cartesian) coordinates in Angstrom.

        Returns:
            * ``(L, grad)`` where ``L`` is a Python float and ``grad`` is an
              ``(N, 3)`` float64 array of ``dL/dx`` in the same frame as
              ``coords``.
        '''
        coords = numpy.ascontiguousarray(coords, numpy.double)
        n = len(coords)
        if self._weights is not None and len(self._weights) != n:
            raise ValueError('weights length {} does not match {} atoms'.format(
                len(self._weights), n))
        from ..clipper_python import Coord_orth
        xmap = self._xmap
        cell = self._cell
        grad = numpy.empty((n, 3), numpy.double)
        total = 0.0
        for j in range(n):
            x, y, z = coords[j]
            cf = Coord_orth(float(x), float(y), float(z)).coord_frac(cell)
            val, gfrac = xmap.interp_cubic_grad_frac(cf)
            gorth = gfrac.grad_orth(cell).dxyz()
            w = 1.0 if self._weights is None else float(self._weights[j])
            total += w * val
            grad[j] = -w * gorth
        return -total, grad


class XrayTargetState:
    '''
    Model-based differentiable crystallographic target (Tier 1/2).

    Unlike :class:`MapTargetState` (which scores against a *fixed* map and so has
    gradients only for coordinates), this target recomputes the model density /
    structure factors each call, so it yields analytic gradients for a
    runtime-selected subset of the eleven per-atom parameters
    ``{X, Y, Z, Uiso, Occ, U11, U22, U33, U12, U13, U23}`` — coordinates, isotropic
    and anisotropic ADPs, and occupancy. The gradient is the exact Agarwal (1978)
    real-space gradient, computed by the C++ ``XrayGradientEvaluator`` (which reuses
    the same ``AtomShapeFn.rho_grad`` machinery as the B-factor/occupancy refiner).

    Two modes, selected by the constructor arguments:

    * **reciprocal** — least-squares against observed structure factors. Pass
      ``fobs`` (``HKL_data<F_sigF>``), ``phi_fom`` (``HKL_data<Phi_fom>`` — its
      figure-of-merit is the σ(A)/ML weight; supply ``fom()=1`` for a small-molecule
      crystal so it reduces to a pure 1/σ²-weighted fit), ``usage``
      (``HKL_data<Flag>``; non-zero flag = working reflection). ``kind='amplitude'``
      (default) minimises ½Σw(k|Fc|−m|Fo|)²; ``kind='intensity'`` minimises the
      SHELX-style ½Σw(Io−s|Fc|²)². Optional ``f_bulk`` adds a fixed bulk-solvent
      contribution (omit for small molecules).

    * **real-space** — least-squares against a fixed target ``Xmap``. Pass
      ``target_map`` (and optionally ``target_origin``, a ``Coord_orth``).

    The numeric core speaks only numpy and has **no PyTorch dependency**; the thin
    torch autograd adapter lives in :mod:`.targets`.

    Args:
        * elements: sequence of ``N`` element symbols (fixed for the lifetime of
          the state — the evaluator seeds them once).
        * param_names: ordered subset of the eleven parameter names; the returned
          gradient has one column per name, in this order (default
          ``('X','Y','Z')``).
        * n_threads: worker threads for EDcalc / FFT / the gradient sum.

    The scale factor (k or s in reciprocal mode; the map scale in real-space mode)
    is held **fixed within a call** so value and gradient are consistent. Call
    :meth:`refresh_scaling` between optimiser steps to re-estimate it.
    '''

    def __init__(self, elements, *, param_names=('X', 'Y', 'Z'),
                 fobs=None, phi_fom=None, usage=None, kind='amplitude',
                 f_bulk=None, target_map=None, target_origin=None, n_threads=1,
                 threaded_density=True):
        from ..clipper_python import AtomShapeFn, Coord_orth
        from ..clipper_python.ext import XrayGradientEvaluator, XrayTargetKind

        self._elements = [str(e) for e in elements]
        self._n = len(self._elements)
        self.param_names = tuple(param_names)
        bad = [p for p in self.param_names if p not in _CANONICAL_PARAMS]
        if bad:
            raise ValueError('unknown parameter name(s): {}'.format(bad))
        T = AtomShapeFn.TYPE
        self._selected = [int(getattr(T, p)) for p in self.param_names]

        if target_map is not None:
            if fobs is not None:
                raise ValueError('give either target_map (real-space) or fobs '
                                 '(reciprocal), not both')
            if target_origin is None:
                target_origin = Coord_orth(0.0, 0.0, 0.0)
            self._ev = XrayGradientEvaluator(
                elements=self._elements, target_map=target_map,
                target_origin=target_origin, n_threads=n_threads,
                threaded_density=threaded_density)
        elif fobs is not None:
            if phi_fom is None or usage is None:
                raise ValueError('reciprocal mode requires fobs, phi_fom and usage')
            kd = (XrayTargetKind.IntensityLS if kind == 'intensity'
                  else XrayTargetKind.AmplitudeLS)
            kwargs = dict(elements=self._elements, fobs=fobs, phi_fom=phi_fom,
                          usage=usage, kind=kd, n_threads=n_threads,
                          threaded_density=threaded_density)
            if f_bulk is not None:
                kwargs['f_bulk'] = f_bulk
            self._ev = XrayGradientEvaluator(**kwargs)
        else:
            raise ValueError('give either target_map (real-space) or fobs '
                             '(reciprocal)')
        # Hold references to the reference-data objects so they (and the data they
        # own) are not garbage-collected while this state is alive. NOTE: Clipper
        # HKL_data holds a *non-owning* pointer to its HKL_info and the binding has
        # no keep_alive, so the caller must ALSO keep whatever owns that HKL_info
        # (e.g. the reflection dataset / crystal manager) alive — otherwise
        # base_hkl_info() dangles and the next evaluation reads garbage.
        self._refs = (fobs, phi_fom, usage, f_bulk, target_map, target_origin)

    @property
    def n_atoms(self):
        return self._n

    def _prepare(self, coords, u_iso, u_aniso, occ, is_aniso):
        n = self._n
        coords = numpy.ascontiguousarray(coords, numpy.double).reshape(n, 3)
        if is_aniso is None:
            is_aniso = numpy.zeros(n, numpy.uint8)
        else:
            is_aniso = numpy.ascontiguousarray(is_aniso, numpy.uint8).reshape(n)
        u_iso = (numpy.zeros(n, numpy.double) if u_iso is None
                 else numpy.ascontiguousarray(u_iso, numpy.double).reshape(n))
        u_aniso = (numpy.zeros((n, 6), numpy.double) if u_aniso is None
                   else numpy.ascontiguousarray(u_aniso, numpy.double).reshape(n, 6))
        occ = (numpy.ones(n, numpy.double) if occ is None
               else numpy.ascontiguousarray(occ, numpy.double).reshape(n))
        return coords, u_iso, u_aniso, occ, is_aniso

    def value_and_gradient(self, coords, u_iso=None, u_aniso=None, occ=None,
                           is_aniso=None, refresh_scale=False):
        '''
        Args:
            * coords: ``(N, 3)`` orthogonal (Cartesian) coordinates (Å).
            * u_iso: ``(N,)`` isotropic U (Å²); used where ``is_aniso`` is 0.
            * u_aniso: ``(N, 6)`` anisotropic U (Å², order u11,u22,u33,u12,u13,u23);
              used where ``is_aniso`` is non-zero.
            * occ: ``(N,)`` occupancies (default all 1).
            * is_aniso: ``(N,)`` flags (0 isotropic, 1 anisotropic; default all 0).
            * refresh_scale: re-estimate the fixed scale from these parameters
              before evaluating (call between optimiser steps). Leave ``False`` for
              a finite-difference check so the scale is frozen across ±ε.

        Returns:
            * ``(L, grad)`` — ``L`` a Python float, ``grad`` an ``(N, P)`` float64
              array of ``dL/dp`` with columns in :attr:`param_names` order.
        '''
        c, ui, ua, oc, ia = self._prepare(coords, u_iso, u_aniso, occ, is_aniso)
        L, grad = self._ev.value_and_gradient(
            c, ui, ua, oc, ia, self._selected, bool(refresh_scale))
        return float(L), grad

    def refresh_scaling(self, coords, u_iso=None, u_aniso=None, occ=None,
                        is_aniso=None):
        '''
        Re-estimate the fixed scale factor from the given parameters (one extra
        EDcalc/FFT). Returns the target value at those parameters. Call every few
        optimiser steps; the fixed-scale approximation (holding it constant within
        a forward/backward) is the standard Agarwal/REFMAC practice.
        '''
        c, ui, ua, oc, ia = self._prepare(coords, u_iso, u_aniso, occ, is_aniso)
        L, _ = self._ev.value_and_gradient(
            c, ui, ua, oc, ia, self._selected, True)
        return float(L)


class SupercellXrayTargetState:
    '''
    Periodic-construct wrapper around :class:`XrayTargetState` for ML force-field
    training (Tier 2).

    A differentiable force-field construct must span at least one full unit cell
    (a supercell where one cell is smaller than the minimum-image radius), so it
    holds ``n_asu`` symmetry copies of the asymmetric unit which break exact
    crystallographic symmetry as the structure relaxes. Clipper, by contrast, works
    on ONE unique ASU and reconstructs the full P1 density from it internally — and
    the structure factor of a symmetry-expanded ASU does not depend on which
    representative you feed, so each ASU can be handed to the evaluator **as-is**,
    with no coordinate transforms.

    This class folds the construct into its ``n_asu`` ASUs, evaluates the single-ASU
    target on each, and scatters the per-ASU gradients back onto the construct atoms.
    The per-ASU results are **summed, not averaged** (``L = Σ_k T_k``): each ASU's own
    deviation is a distinct, valid training signal, and averaging would discard it. A
    general-position construct atom belongs to exactly one ASU and receives that ASU's
    gradient; a special-position atom shared between several ASUs accumulates (scatter-
    add) a contribution from each.

    Args:
        * asu_elements: element symbols of ONE asymmetric unit (length ``N_asu``); all
          ASUs are compositionally identical (they are symmetry copies).
        * groups: sequence of ``n_asu`` integer index arrays into the construct atom
          set, each of length ``N_asu`` and aligned by ASU slot — ``groups[k][i]`` is
          the construct index of ASU ``k``'s ``i``-th atom. A general atom appears in
          exactly one group; a special-position atom shared between ``p`` ASUs appears
          in ``p`` groups.
        * multiplicity: optional ``(M,)`` per-construct-atom site multiplicity (1 for a
          general position, ``m`` for an ``m``-fold special position). Occupancies fed
          to Clipper are divided by it — the ASU convention the EDcalc special-position
          correction expects (feed ``occ/m``; the C++ layer scales the density back up).
          Default: all ones (all general positions).
        * remaining keyword args (param_names, fobs, phi_fom, usage, kind, f_bulk,
          target_map, target_origin, n_threads, threaded_density): forwarded verbatim
          to the per-ASU :class:`XrayTargetState`.

    A separate :class:`XrayTargetState` is built per ASU so each holds its own frozen
    scale (independent samples); they share the same read-only reference data. The
    gradient routing matches :class:`XrayTargetState` (columns in ``param_names``
    order) at construct scale, so the :func:`chimerax.clipper.diff.targets.xray_loss`
    torch adapter works unchanged with a construct-sized coordinate tensor.

    NB: special-position **gradient** exactness (the interplay between the EDcalc
    ``×multiplicity`` density correction and the uncorrected Agarwal gradient kernel)
    is not yet numerically verified; general positions are the validated common case.
    '''

    def __init__(self, asu_elements, groups, *, multiplicity=None,
                 param_names=('X', 'Y', 'Z'), fobs=None, phi_fom=None, usage=None,
                 kind='amplitude', f_bulk=None, target_map=None, target_origin=None,
                 n_threads=1, threaded_density=True):
        self._groups = [numpy.ascontiguousarray(g, numpy.intp).reshape(-1)
                        for g in groups]
        if not self._groups:
            raise ValueError('groups must contain at least one ASU')
        self._n_asu = len(self._groups)
        nslot = len(asu_elements)
        for k, g in enumerate(self._groups):
            if g.shape[0] != nslot:
                raise ValueError(
                    'groups[{}] has {} atoms but the ASU has {}'.format(
                        k, g.shape[0], nslot))
        # Construct atom count inferred from the grouping (every construct atom
        # belongs to at least one ASU, since the construct IS the union of ASUs).
        self._n_construct = 1 + int(max(int(g.max()) for g in self._groups))

        # One evaluator per ASU: each keeps its own frozen scale for consistent
        # value/gradient within a call and correct finite-difference freezing. They
        # share the same read-only reference data objects (fobs/phi_fom/usage/…).
        common = dict(param_names=param_names, fobs=fobs, phi_fom=phi_fom,
                      usage=usage, kind=kind, f_bulk=f_bulk, target_map=target_map,
                      target_origin=target_origin, n_threads=n_threads,
                      threaded_density=threaded_density)
        self._states = [XrayTargetState(asu_elements, **common)
                        for _ in range(self._n_asu)]
        self.param_names = self._states[0].param_names
        self._occ_col = (self.param_names.index('Occ')
                         if 'Occ' in self.param_names else None)
        self._mult = (None if multiplicity is None
                      else numpy.ascontiguousarray(multiplicity, numpy.double).reshape(-1))

    @property
    def n_atoms(self):
        '''Number of atoms in the whole construct (the fold/unfold operates at
        this scale; gradients are returned as ``(n_atoms, P)``).'''
        return self._n_construct

    @property
    def n_asu(self):
        return self._n_asu

    @property
    def asu_n_atoms(self):
        return self._states[0].n_atoms

    def value_and_gradient(self, coords, u_iso=None, u_aniso=None, occ=None,
                           is_aniso=None, refresh_scale=False):
        '''
        Fold → evaluate per ASU → scatter-add unfold. Signature and column
        convention match :class:`XrayTargetState`, but at construct scale: ``coords``
        is ``(M, 3)`` over all ``M`` construct atoms and the returned gradient is
        ``(M, P)`` with columns in :attr:`param_names` order. ``L = Σ_k T_k``.
        '''
        M = self._n_construct
        coords = numpy.asarray(coords, numpy.double).reshape(M, 3)
        u_iso_f = None if u_iso is None else numpy.asarray(u_iso, numpy.double).reshape(M)
        u_aniso_f = None if u_aniso is None else numpy.asarray(u_aniso, numpy.double).reshape(M, 6)
        is_aniso_f = None if is_aniso is None else numpy.asarray(is_aniso, numpy.uint8).reshape(M)
        occ_f = (numpy.ones(M, numpy.double) if occ is None
                 else numpy.asarray(occ, numpy.double).reshape(M))

        P = len(self.param_names)
        grad = numpy.zeros((M, P), numpy.double)
        L_total = 0.0
        for k, g in enumerate(self._groups):
            occ_g = occ_f[g]
            if self._mult is not None:
                occ_g = occ_g / self._mult[g]
            Lk, gk = self._states[k].value_and_gradient(
                coords[g],
                None if u_iso_f is None else u_iso_f[g],
                None if u_aniso_f is None else u_aniso_f[g],
                occ_g,
                None if is_aniso_f is None else is_aniso_f[g],
                refresh_scale=refresh_scale)
            L_total += Lk
            # occ was fed as occ/mult, so dL/docc_construct = (dL/docc_fed)/mult.
            if self._occ_col is not None and self._mult is not None:
                gk = numpy.array(gk, copy=True)
                gk[:, self._occ_col] /= self._mult[g]
            # Scatter-add: += (not =) so a special atom shared across ASUs
            # accumulates its contribution from each group it appears in.
            numpy.add.at(grad, g, gk)
        return float(L_total), grad

    def refresh_scaling(self, coords, u_iso=None, u_aniso=None, occ=None,
                        is_aniso=None):
        '''Re-fit every ASU's scale at the given construct parameters; returns
        ``L = Σ_k T_k`` there.'''
        L, _ = self.value_and_gradient(coords, u_iso, u_aniso, occ, is_aniso,
                                       refresh_scale=True)
        return L


def _u6_congruence(R):
    '''6×6 matrix ``A`` mapping stored ADP vectors ``(u11,u22,u33,u12,u13,u23)`` under
    the operator whose orthogonal rotation part is ``R``, matching Clipper's
    ``U_aniso_orth.transform``: ``U' = rᵀ·U·r`` with ``r = R⁻¹`` (NB **not** ``R·U·Rᵀ``
    — these differ when ``R`` is not exactly orthogonal, as symop ``Place`` rotations
    in a non-orthogonal cell generally are not). This is the exact linear map the
    engine's expansion / ``collapse_to_asu`` apply, so folding the box ADPs with it
    reproduces the ASU-frame ADPs; its transpose unfolds the aniso-U gradient
    (``dL/dU_box = Aᵀ·dL/dU_asu``). Column ``k`` is the transformed ``k``-th basis
    tensor.'''
    r = numpy.linalg.inv(R)
    idx = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
    A = numpy.empty((6, 6), numpy.double)
    for k, (a, b) in enumerate(idx):
        E = numpy.zeros((3, 3), numpy.double)
        E[a, b] = 1.0
        E[b, a] = 1.0                      # symmetric basis tensor
        RER = r.T @ E @ r                  # r·U·rᵀ with r = R⁻¹  (Clipper convention)
        A[:, k] = (RER[0, 0], RER[1, 1], RER[2, 2], RER[0, 1], RER[0, 2], RER[1, 2])
    return A


class EnsembleXrayTargetState:
    '''
    Small-molecule **ensemble** crystallographic target for a periodic MD box.

    The small-molecule counterpart to :class:`SupercellXrayTargetState`. Where the
    macromolecular path scores each symmetry copy as its own full-occupancy ASU (so
    the bulk-solvent mask stays correct), the small-molecule regime has no bulk
    solvent, so the crystallographically-standard **occupancy-weighted overlay** is
    used: a full P1 box (``n_asu`` symmetry copies produced by ``clipper unitcells`` /
    ``symcopies``, carrying a ``SymmetryExpansion``) is folded onto a single ASU at
    occupancy ``1/n_asu`` and scored — in **one** SF calc — against the observed
    structure factors through :class:`XrayTargetState` (reciprocal, **no** ``f_bulk``).

    No explicit per-atom multiplicity is needed: exact special-position dedup at
    expansion leaves a multiplicity-``m`` site with only ``n_asu/m`` images, so at
    ``1/n_asu`` each they overlay to occ ``1/m`` — exactly what Clipper's EDcalc
    ``×multiplicity`` correction then restores to a full peak.

    **Selectable gradient terms.** ``param_names`` chooses any subset of the eleven
    per-atom parameters (default ``('X','Y','Z')`` for the "frozen snapshots" use;
    include ``U11..U23`` to estimate anisotropic ADPs, ``Uiso``/``Occ`` as needed).
    The box→ASU fold rotates coordinates (``asu = R·box + t``) *and* aniso tensors
    (``U_asu = R·U·Rᵀ``), so the returned gradient is unfolded back into the **box**
    frame per term: **X,Y,Z** by ``Rᵀ``, **U11..U23** by the congruence-transpose
    ``Aᵀ``, and **Uiso**/**Occ** unchanged (rotation-invariant). Because coordinates
    couple within their triple and aniso-U within its sextet, each of those groups
    must be requested whole (all three / all six) or not at all.

    Interface mirrors :class:`XrayTargetState` — ``value_and_gradient(box_coords,
    u_iso, u_aniso, occ, is_aniso, refresh_scale) -> (L, (M, P))`` over all ``M`` box
    atoms, columns in ``param_names`` order — so
    :func:`chimerax.clipper.diff.targets.xray_loss` wraps it directly (and, for the
    coords-only default, ``crystallographic_loss`` does too). All coordinate/ADP/occ
    inputs are in the **box** frame; the fold/unfold is internal. ADP and occupancy
    args left ``None`` fall back to builder-supplied defaults, so a coords-only caller
    can simply pass ``box_coords``.

    Args:
        * elements: ``M`` box-atom element symbols.
        * inv_rot / inv_trn: ``(M,3,3)`` / ``(M,3)`` per-atom **inverse** operator
          (box→ASU rigid fold); from ``exp.inverse_operators`` via the builder.
        * param_names: subset of the eleven parameters (see above).
        * fobs / phi_fom / usage: observed data (reciprocal; small molecules supply
          ``phi_fom.fom()==1``).
        * u_iso / u_aniso / is_aniso / occupancy: default **box-frame** ADPs / flags /
          occupancy (the builder passes the model's values and ``1/n_asu``).
    '''

    _ALL = ('X', 'Y', 'Z', 'Uiso', 'Occ', 'U11', 'U22', 'U33', 'U12', 'U13', 'U23')

    def __init__(self, elements, inv_rot, inv_trn, *, param_names=('X', 'Y', 'Z'),
                 fobs, phi_fom, usage, kind='amplitude', u_iso=None, u_aniso=None,
                 is_aniso=None, occupancy=None, n_threads=1, threaded_density=True):
        self.param_names = tuple(param_names)
        bad = [p for p in self.param_names if p not in self._ALL]
        if bad:
            raise ValueError('unknown parameter name(s): {}'.format(bad))
        xyz = [p for p in ('X', 'Y', 'Z') if p in self.param_names]
        if xyz and len(xyz) != 3:
            raise ValueError('X, Y, Z must be selected together (the fold rotates them)')
        u6 = [p for p in ('U11', 'U22', 'U33', 'U12', 'U13', 'U23') if p in self.param_names]
        if u6 and len(u6) != 6:
            raise ValueError('U11..U23 must be selected together (the fold rotates the tensor)')

        n = self._M = len(elements)
        self._R = numpy.ascontiguousarray(inv_rot, numpy.double).reshape(n, 3, 3)
        self._t = numpy.ascontiguousarray(inv_trn, numpy.double).reshape(n, 3)
        # Per-atom 6×6 congruence for the box→ASU aniso fold (U_asu = R·U·Rᵀ).
        self._A = numpy.stack([_u6_congruence(self._R[a]) for a in range(n)])
        self._state = XrayTargetState(
            elements, param_names=self.param_names, fobs=fobs, phi_fom=phi_fom,
            usage=usage, kind=kind, n_threads=n_threads,
            threaded_density=threaded_density)          # NO f_bulk (small molecule)
        # Default BOX-frame ADPs / occupancy (used when a call leaves them None).
        self._def_uiso = (None if u_iso is None
                          else numpy.ascontiguousarray(u_iso, numpy.double).reshape(n))
        self._def_uaniso = (None if u_aniso is None
                            else numpy.ascontiguousarray(u_aniso, numpy.double).reshape(n, 6))
        self._def_isan = (None if is_aniso is None
                          else numpy.ascontiguousarray(is_aniso, numpy.uint8).reshape(n))
        if occupancy is None:
            self._def_occ = None
        elif numpy.isscalar(occupancy):
            self._def_occ = numpy.full(n, float(occupancy), numpy.double)
        else:
            self._def_occ = numpy.ascontiguousarray(occupancy, numpy.double).reshape(n)
        col = {name: i for i, name in enumerate(self.param_names)}
        self._xyz_cols = [col[c] for c in ('X', 'Y', 'Z')] if xyz else None
        self._u6_cols = [col[c] for c in ('U11', 'U22', 'U33', 'U12', 'U13', 'U23')] if u6 else None

    @property
    def n_atoms(self):
        return self._M

    def value_and_gradient(self, box_coords, u_iso=None, u_aniso=None, occ=None,
                           is_aniso=None, refresh_scale=False):
        '''
        Fold box→ASU, score, unfold. ``box_coords`` and any ADP/occ arrays are in the
        **box** frame (whole box, ``M`` atoms); ADP/occ args left ``None`` use the
        builder-supplied defaults. Returns ``(L, grad)`` with ``grad`` ``(M, P)`` =
        ``dL/d(box params)``, columns in :attr:`param_names` order. Scale is fixed
        unless ``refresh_scale`` (call :meth:`refresh_scaling` between optimiser steps).
        '''
        n = self._M
        box = numpy.ascontiguousarray(box_coords, numpy.double).reshape(n, 3)
        b_uiso = self._def_uiso if u_iso is None else numpy.ascontiguousarray(u_iso, numpy.double).reshape(n)
        b_uan = self._def_uaniso if u_aniso is None else numpy.ascontiguousarray(u_aniso, numpy.double).reshape(n, 6)
        b_occ = self._def_occ if occ is None else numpy.ascontiguousarray(occ, numpy.double).reshape(n)
        b_isan = self._def_isan if is_aniso is None else numpy.ascontiguousarray(is_aniso, numpy.uint8).reshape(n)

        # Fold to the ASU frame: coords rigidly, aniso-U by the congruence; Uiso/Occ
        # are rotation-invariant and pass straight through.
        asu = numpy.einsum('aij,aj->ai', self._R, box) + self._t
        asu_uan = None if b_uan is None else numpy.einsum('akl,al->ak', self._A, b_uan)
        L, g_asu = self._state.value_and_gradient(
            asu, u_iso=b_uiso, u_aniso=asu_uan, occ=b_occ, is_aniso=b_isan,
            refresh_scale=refresh_scale)
        g_asu = numpy.ascontiguousarray(g_asu, numpy.double)

        # Unfold ASU→box: X,Y,Z by Rᵀ; U11..U23 by Aᵀ; Uiso/Occ unchanged.
        g_box = g_asu.copy()
        if self._xyz_cols is not None:
            g_box[:, self._xyz_cols] = numpy.einsum('aij,ai->aj', self._R, g_asu[:, self._xyz_cols])
        if self._u6_cols is not None:
            g_box[:, self._u6_cols] = numpy.einsum('akl,ak->al', self._A, g_asu[:, self._u6_cols])
        return float(L), g_box

    def refresh_scaling(self, box_coords, u_iso=None, u_aniso=None, occ=None,
                        is_aniso=None):
        '''Re-fit the scale at the given box parameters; returns ``L`` there.'''
        L, _ = self.value_and_gradient(box_coords, u_iso, u_aniso, occ, is_aniso,
                                       refresh_scale=True)
        return L
