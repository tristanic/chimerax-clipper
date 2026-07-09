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
PyTorch autograd adapters for the differentiable crystallographic targets.

These are deliberately thin: the numeric cores in :mod:`.state` return both the
loss value and its analytic gradient in one call, so the ``autograd.Function``
just stashes the gradient in ``forward`` and scales it by the upstream gradient
in ``backward``.

PyTorch is an **optional** dependency of ChimeraX-Clipper. It is imported lazily
here (via :func:`chimerax.clipper.diff.require_torch`) and never at package
import time, so the bundle loads without torch present. Import this module only
when you actually intend to build a torch graph.
'''


def _autograd_function():
    '''
    Build and return the ``torch.autograd.Function`` subclass. Deferred into a
    function so that merely importing this module does not require torch; the
    class is created on first use and cached on the module.
    '''
    from . import require_torch
    torch = require_torch()

    class _CrystallographicTarget(torch.autograd.Function):
        '''
        Generic adapter: ``forward(coords, state)`` returns the scalar loss;
        ``backward`` returns ``upstream * dL/dcoords``. ``state`` is any object
        with ``value_and_gradient(coords_ndarray) -> (float, (N,3) ndarray)``.
        '''
        @staticmethod
        def forward(ctx, coords, state):
            coords_np = coords.detach().cpu().numpy()
            value, grad = state.value_and_gradient(coords_np)
            grad_t = torch.as_tensor(grad, dtype=coords.dtype, device=coords.device)
            ctx.save_for_backward(grad_t)
            return coords.new_tensor(float(value))

        @staticmethod
        def backward(ctx, grad_output):
            (grad_t,) = ctx.saved_tensors
            # grad_output is a scalar tensor; state is not differentiable.
            return grad_output * grad_t, None

    return _CrystallographicTarget


_FUNCTION = None


def crystallographic_loss(coords, state):
    '''
    Differentiable crystallographic loss.

    Args:
        * coords: an ``(N, 3)`` torch tensor of orthogonal (Cartesian) Angstrom
          coordinates, ``requires_grad=True`` to receive gradients.
        * state: a target-state object from :mod:`.state` (e.g.
          :class:`~chimerax.clipper.diff.state.MapTargetState`) exposing
          ``value_and_gradient``.

    Returns:
        * a scalar torch tensor; call ``.backward()`` to propagate ``dL/dcoords``
          (and thence into whatever produced ``coords``).
    '''
    global _FUNCTION
    if _FUNCTION is None:
        _FUNCTION = _autograd_function()
    return _FUNCTION.apply(coords, state)


def _xray_autograd_function():
    '''
    Build the multi-parameter ``torch.autograd.Function`` for
    :class:`~chimerax.clipper.diff.state.XrayTargetState`. The state returns an
    ``(N, P)`` gradient with columns in ``state.param_names`` order; ``backward``
    routes each column back to the matching input tensor (coords / u_iso / u_aniso
    / occ). Deferred so importing this module never requires torch.
    '''
    from . import require_torch
    torch = require_torch()

    class _XrayTarget(torch.autograd.Function):
        @staticmethod
        def forward(ctx, state, coords, u_iso, u_aniso, occ, is_aniso, refresh_scale):
            def _np(t):
                return None if t is None else t.detach().cpu().numpy()
            value, grad = state.value_and_gradient(
                _np(coords), _np(u_iso), _np(u_aniso), _np(occ), _np(is_aniso),
                refresh_scale=refresh_scale)
            # Column index of each parameter name in the (N, P) gradient.
            cols = {name: i for i, name in enumerate(state.param_names)}
            dev, dt = coords.device, coords.dtype
            g = torch.as_tensor(grad, dtype=dt, device=dev)

            def _pick(names, ncol):
                # Assemble a (N, ncol) gradient block; missing params → 0.
                out = torch.zeros((g.shape[0], ncol), dtype=dt, device=dev)
                any_present = False
                for k, nm in enumerate(names):
                    if nm in cols:
                        out[:, k] = g[:, cols[nm]]
                        any_present = True
                return out if any_present else None

            ctx._coords_g = _pick(('X', 'Y', 'Z'), 3) if coords.requires_grad else None
            ctx._uiso_g = (_pick(('Uiso',), 1)
                           if (u_iso is not None and u_iso.requires_grad) else None)
            ctx._uaniso_g = (_pick(('U11', 'U22', 'U33', 'U12', 'U13', 'U23'), 6)
                             if (u_aniso is not None and u_aniso.requires_grad) else None)
            ctx._occ_g = (_pick(('Occ',), 1)
                          if (occ is not None and occ.requires_grad) else None)
            return coords.new_tensor(float(value))

        @staticmethod
        def backward(ctx, grad_output):
            def _scale(blk):
                if blk is None:
                    return None
                b = grad_output * blk
                return b.reshape(-1) if b.shape[1] == 1 else b
            # Order matches forward's args: state, coords, u_iso, u_aniso, occ,
            # is_aniso, refresh_scale. Non-tensor / non-diff args return None.
            return (None, _scale(ctx._coords_g), _scale(ctx._uiso_g),
                    _scale(ctx._uaniso_g), _scale(ctx._occ_g), None, None)

    return _XrayTarget


_XRAY_FUNCTION = None


def xray_loss(state, coords, u_iso=None, u_aniso=None, occ=None, is_aniso=None,
              refresh_scale=False):
    '''
    Differentiable crystallographic loss for an
    :class:`~chimerax.clipper.diff.state.XrayTargetState`.

    Pass the per-atom parameter tensors you want gradients for with
    ``requires_grad=True``; any of ``coords`` (``(N,3)``), ``u_iso`` (``(N,)``),
    ``u_aniso`` (``(N,6)``) and ``occ`` (``(N,)``) may be differentiable. The
    columns the state actually produces are governed by its ``param_names``; a
    requested tensor whose parameters are absent from ``param_names`` simply
    receives no gradient. ``is_aniso`` (``(N,)``) and ``refresh_scale`` are
    non-differentiable.

    Returns a scalar torch tensor; call ``.backward()`` to propagate the
    parameter gradients.
    '''
    global _XRAY_FUNCTION
    if _XRAY_FUNCTION is None:
        _XRAY_FUNCTION = _xray_autograd_function()
    return _XRAY_FUNCTION.apply(state, coords, u_iso, u_aniso, occ, is_aniso,
                                refresh_scale)
