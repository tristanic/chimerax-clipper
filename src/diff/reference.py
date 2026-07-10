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
Pure-Python reference ("oracle") for the real-space model-density target.

This is retained from the `diffmap` line of work. It is **functionally superseded**
by :class:`chimerax.clipper.diff.state.XrayTargetState` used in its real-space mode
(``target_map=...``), which computes the same real-space least-squares value and
gradient through the compiled C++ ``XrayGradientEvaluator`` for all eleven per-atom
parameters. It is kept for two reasons:

1. **Independent cross-check.** Being a from-scratch numpy implementation over
   ``AtomShapeFn.rho`` / ``AtomShapeFn.rho_grad`` (deliberately written for clarity,
   not speed), it is a *fourth* independent way to compute the real-space
   value+gradient — beyond the threaded EDcalc, the single-threaded companion, and
   the finite-difference checks — useful for catching bugs common to the C++ paths.

2. **Two ideas worth carrying into the C++ engine.**
   - ``mask='atoms'``: the least-squares is evaluated only over grid voxels within an
     atom's cutoff radius, deliberately excluding far-field density (bulk solvent in
     proteins; unmodelled symmetry-mate peaks when only the ASU is supplied). The C++
     real-space route currently fits the *whole* grid, so this masking is a genuine
     feature gap to port into ``xray_gradient.cpp`` (real-space branch of
     ``evaluate()``).
   - Periodic-grid export: the true periodic dimensions must be read from
     ``grid_sampling`` (nu, nv, nw) via ``export_section_numpy`` over ``[0, n)`` — NOT
     via ``export_numpy``, which historically returned the ASU one grid-plane short in
     each axis (fixed in the binding as of this branch; the workaround is retained here
     for older builds).

Not wired into the torch adapter or the package's public surface — it is a reference
implementation, exercised from tests / finite-difference cross-checks.
'''

import numpy


class ModelDensityTargetState:
    '''
    Tier 1: differentiable real-space least-squares target between a *model*
    density built from the current atoms and a *fixed* experimental map, yielding
    analytic gradients for a runtime-selected subset of the per-atom parameters via
    the Agarwal (1978) real-space formula.

    Unlike :class:`MapTargetState` (which samples a fixed field at atom positions
    and so gives coordinate gradients only), this target rebuilds the model's own
    density, so its loss depends on every atomic parameter it is asked about —
    occupancy and isotropic B here (anisotropic-U in a follow-up). That is what
    makes ``dL/dOcc`` and ``dL/dU`` meaningful.

    Framing: the atoms are the full explicit contents of a **P1 cell with periodic
    boundary conditions** (the molecular-dynamics convention), so every atom is an
    independent degree of freedom and no crystallographic symmetry operators enter
    the gradient. Crystallographic symmetry of the experimental data is already
    baked into the fixed target map.

    Loss (least squares over the whole grid)::

        L = 0.5 * sum_x (s * rho_model(x) - rho_exp(x))^2

    ``scale='auto'`` sets ``s`` to the optimal LS value each call and holds it fixed
    for the gradient (envelope theorem → its first-order contribution vanishes);
    pass a float to fix it (use ``scale=1.0`` for the cleanest finite-difference
    check, which then has no scale-derivative term).

    Reference ("oracle") implementation: ``rho_model`` is built from
    ``AtomShapeFn.rho`` on the same grid the gradient uses ``AtomShapeFn.rho_grad``,
    so value and gradient share an identical density model and agree with a finite
    difference to ~machine precision. Written for clarity, not speed — the C++
    ``xray_gradient`` extension supersedes it for production.

    ADP mode: pass ``u_iso`` (N,) for isotropic atoms, or ``u_aniso`` (N,6) for
    anisotropic atoms (components U11,U22,U33,U12,U13,U23 — the Clipper /
    ``AtomShapeFn.TYPE`` order). The selectable params follow the mode: isotropic
    exposes ``Uiso``; anisotropic exposes ``U11..U23`` (the six ADP gradients used
    for Hessian/ADP training). Exactly one of ``u_iso`` / ``u_aniso`` must be given.

    Args:
        * target_xmap: Clipper Xmap of the fixed experimental density (supplies the
          P1 grid, cell and rho_exp).
        * elements: length-N Clipper scattering-factor identifiers.
        * occupancies: (N,) occupancies.
        * u_iso: (N,) isotropic U (A^2)   [isotropic mode], or
        * u_aniso: (N,6) U11,U22,U33,U12,U13,U23 (A^2)   [anisotropic mode].
        * param_names: subset of the allowed set for the mode:
          isotropic  -> ('X','Y','Z','Uiso','Occ');
          anisotropic-> ('X','Y','Z','U11','U22','U33','U12','U13','U23','Occ').
        * radius: density cutoff (A); grid points within it of an atom contribute.
        * scale: 'auto' or a float.
    '''

    _ALLOWED_ISO = ('X', 'Y', 'Z', 'Uiso', 'Occ')
    _ALLOWED_ANISO = ('X', 'Y', 'Z', 'U11', 'U22', 'U33', 'U12', 'U13', 'U23', 'Occ')

    def __init__(self, target_xmap, elements, occupancies, u_iso=None, u_aniso=None,
                 param_names=('X', 'Y', 'Z'), radius=3.0, scale='auto', mask='atoms'):
        if mask not in ('atoms', 'cell'):
            raise ValueError("mask must be 'atoms' or 'cell'; got %r" % (mask,))
        self._mask_mode = mask
        if (u_iso is None) == (u_aniso is None):
            raise ValueError('provide exactly one of u_iso (isotropic) or '
                             'u_aniso (anisotropic)')
        self._aniso = u_aniso is not None
        allowed = self._ALLOWED_ANISO if self._aniso else self._ALLOWED_ISO
        bad = [p for p in param_names if p not in allowed]
        if bad:
            raise ValueError('%s mode supports params %s; got unsupported %r' % (
                'anisotropic' if self._aniso else 'isotropic', allowed, bad))
        self._elements = list(elements)
        self._occ = numpy.ascontiguousarray(occupancies, numpy.double)
        if self._aniso:
            self._uaniso = numpy.ascontiguousarray(u_aniso, numpy.double)
            if self._uaniso.shape[1:] != (6,):
                raise ValueError('u_aniso must be (N,6); got %r' % (self._uaniso.shape,))
        else:
            self._uiso = numpy.ascontiguousarray(u_iso, numpy.double)
        self.param_names = tuple(param_names)
        self._radius = float(radius)
        self._scale = scale
        self._cell = target_xmap.cell
        # True periodic grid dimensions from the sampling (NOT export_numpy, which
        # historically returned the ASU as grid_asu.max()-min() = one plane short in
        # each axis and would corrupt periodic wrapping for atoms near a cell edge).
        from ..clipper_python import Coord_frac, AtomShapeFn, Coord_grid
        gs = target_xmap.grid_sampling
        self._nu, self._nv, self._nw = int(gs.nu), int(gs.nv), int(gs.nw)
        # rho_exp on the full periodic grid (nu, nv, nw), read losslessly via
        # export_section_numpy over [0, n) in each axis. float64 for a clean check.
        self._rho_exp = numpy.asarray(target_xmap.export_section_numpy(
            Coord_grid(0, 0, 0), Coord_grid(self._nu, self._nv, self._nw)), numpy.double)
        cell = self._cell
        self._c_u = numpy.asarray(Coord_frac(1.0 / self._nu, 0.0, 0.0).coord_orth(cell).xyz, numpy.double)
        self._c_v = numpy.asarray(Coord_frac(0.0, 1.0 / self._nv, 0.0).coord_orth(cell).xyz, numpy.double)
        self._c_w = numpy.asarray(Coord_frac(0.0, 0.0, 1.0 / self._nw).coord_orth(cell).xyz, numpy.double)
        import math
        self._steps = (
            int(math.ceil(self._radius / numpy.linalg.norm(self._c_u))),
            int(math.ceil(self._radius / numpy.linalg.norm(self._c_v))),
            int(math.ceil(self._radius / numpy.linalg.norm(self._c_w))),
        )
        T = AtomShapeFn.TYPE
        self._types = [getattr(T, p) for p in self.param_names]

    def _atom_box(self, coord):
        '''Return (wrapped_indices [(iu,iv,iw)...], orth_coords [np(3)...]) for grid
        points within the cutoff radius of an atom at orthogonal ``coord``.'''
        from ..clipper_python import Coord_orth
        cf = Coord_orth(float(coord[0]), float(coord[1]), float(coord[2])).coord_frac(self._cell)
        fu, fv, fw = cf.uvw
        u0 = int(round(fu * self._nu)); v0 = int(round(fv * self._nv)); w0 = int(round(fw * self._nw))
        du, dv, dw = self._steps
        r2 = self._radius * self._radius
        cu, cv, cw = self._c_u, self._c_v, self._c_w
        idxs = []; orths = []
        for iu in range(u0 - du, u0 + du + 1):
            for iv in range(v0 - dv, v0 + dv + 1):
                bu = iu * cu + iv * cv  # partial
                for iw in range(w0 - dw, w0 + dw + 1):
                    orth = bu + iw * cw
                    dx = orth - coord
                    if dx.dot(dx) <= r2:
                        idxs.append((iu % self._nu, iv % self._nv, iw % self._nw))
                        orths.append(orth)
        return idxs, orths

    def _make_shape(self, coord, j):
        '''Construct the AtomShapeFn for atom j at orthogonal ``coord``, in the
        appropriate ADP mode, with the selected Agarwal parameters set.'''
        from ..clipper_python import AtomShapeFn, Coord_orth, U_aniso_orth
        co = Coord_orth(float(coord[0]), float(coord[1]), float(coord[2]))
        if self._aniso:
            u = self._uaniso[j]
            adp = U_aniso_orth(float(u[0]), float(u[1]), float(u[2]),
                               float(u[3]), float(u[4]), float(u[5]))
            asf = AtomShapeFn(co, self._elements[j], adp, float(self._occ[j]))
        else:
            asf = AtomShapeFn(co, self._elements[j], float(self._uiso[j]), float(self._occ[j]))
        asf.set_agarwal_params(self._types)
        return asf

    def value_and_gradient(self, coords):
        '''Args: coords (N,3) orthogonal Angstrom. Returns (L, grad) with grad shape
        (N, P), P = len(param_names), columns ordered as param_names.'''
        coords = numpy.ascontiguousarray(coords, numpy.double)
        n = len(coords)
        from ..clipper_python import Coord_orth
        rho_model = numpy.zeros((self._nu, self._nv, self._nw), numpy.double)
        mask = None if self._mask_mode == 'cell' else numpy.zeros(rho_model.shape, bool)
        shapes = []; boxes = []
        # Pass 1: build the model density (and cache each atom's box). In 'atoms'
        # mode, mark every box voxel: the loss is evaluated only there, excluding
        # far-field density (bulk solvent in proteins; unmodelled symmetry-mate
        # peaks when only the ASU is supplied).
        for j in range(n):
            asf = self._make_shape(coords[j], j)
            idxs, orths = self._atom_box(coords[j])
            for (iu, iv, iw), orth in zip(idxs, orths):
                rho_model[iu, iv, iw] += asf.rho(Coord_orth(float(orth[0]), float(orth[1]), float(orth[2])))
                if mask is not None:
                    mask[iu, iv, iw] = True
            shapes.append(asf); boxes.append((idxs, orths))
        # Scale (fit over the evaluated region) + residual.
        if mask is None:
            rm, re = rho_model, self._rho_exp
        else:
            rm, re = rho_model[mask], self._rho_exp[mask]
        if self._scale == 'auto':
            denom = float(numpy.sum(rm * rm))
            s = float(numpy.sum(rm * re) / denom) if denom > 0 else 1.0
        else:
            s = float(self._scale)
        d = s * rho_model - self._rho_exp
        L = 0.5 * float(numpy.sum((d if mask is None else d[mask]) ** 2))
        # Pass 2: Agarwal gradient. dL/dp_j = s * sum_x d(x) * d(rho_j)/dp_j.
        grad = numpy.zeros((n, len(self._types)), numpy.double)
        for j in range(n):
            asf = shapes[j]; idxs, orths = boxes[j]
            acc = numpy.zeros(len(self._types), numpy.double)
            for (iu, iv, iw), orth in zip(idxs, orths):
                _, g = asf.rho_grad(Coord_orth(float(orth[0]), float(orth[1]), float(orth[2])))
                acc += d[iu, iv, iw] * numpy.asarray(g, numpy.double)
            grad[j] = s * acc
        return L, grad
