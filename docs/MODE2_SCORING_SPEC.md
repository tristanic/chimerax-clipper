# Mode 2 — in-the-loop differentiable scoring (Phase 5 design note)

Status: **design only**, no implementation yet. This records the agreed approach
for the GARNET Phase-5 "score a relaxed model against experimental data, with
gradients" capability, per `garnet-isolde/docs/clipper_interface.md`.

## Recommendation: real-space density-fit, not reciprocal-space R

Two scoring routes were considered:

1. **Reciprocal-space R-factor** (what Mode 1 computes). This is a scalar read-out
   and is **not differentiable** as implemented — there is no gradient of R w.r.t.
   atomic coordinates wired up. Building one would mean differentiating the
   structure-factor sum and the scaling, then back-propagating to coordinates:
   substantial new code.

2. **Real-space density-fit (MDFF-style)** — score the model against an
   experimental real-space map by sampling the map (and its gradient) at atom
   positions. **ISOLDE already implements exactly this, differentiably**: its
   molecular-dynamics flexible-fitting engine applies per-atom forces derived from
   a map, i.e. it already computes ∂(fit score)/∂(atomic coordinates) analytically
   and at interactive rates. This is the GARNET-side preference and the right
   choice.

So Phase 5 should reuse ISOLDE's map-fit potential rather than build a new
reciprocal-space gradient.

## Mechanics

Per structure (one-time setup):
- Read model + reflections (Mode 1 machinery already does this).
- Compute experimental map coefficients and a real-space map. For small molecules,
  a 2mFo–DFc-style or simply an Fobs map can be generated with the existing Clipper
  `SFweight`/`Xmap` path (cf. `map_calc.generate_map_coeffs`). The map is computed
  **once** and cached for the structure.

In the training loop (per candidate coordinate set):
- Feed coordinates to ISOLDE's map-fit potential; read back the scalar fit score
  and the per-atom gradient. No re-FFT is needed — only map + gradient sampling at
  the atom positions, which is cheap.

## Gradients

- **Analytic, via ISOLDE** — no finite differences required. ISOLDE's MDFF engine
  is built around ∂(map score)/∂(xyz); that is its inner loop.
- (Clipper also has reciprocal-space Fcalc derivatives, but they are not wired for
  autograd; the real-space path is the ready one.)

## Performance / batching

- Map generation is once-per-structure and fast at small-molecule scale (small
  cells, tens of atoms).
- In-loop scoring is map+gradient sampling at atom positions — ISOLDE does this for
  thousands of atoms at interactive frame rates, so small-molecule scoring is
  effectively free.
- Structures are independent → embarrassingly parallel; precompute maps, then score.

## Why this also dissolves the supercell / symmetry question

GARNET's relaxed coordinates come from an n×n×n MD supercell whose symmetry is
broken. Real-space scoring compares coordinates against a map placed in real space,
so there is no need to reconcile P1 (broken-symmetry) coordinates with the
crystallographic space group in reciprocal space — the comparison is well-posed by
construction. (If a reciprocal-space R were ever wanted in-loop, that convention
would have to be settled then; the real-space route avoids it.)

## Open items for when Phase 5 is built

- Confirm the ISOLDE API entry point for "score coordinates against a map + return
  gradient" and whether it can be driven headlessly without the full ISOLDE UI.
- Decide map type (Fobs vs sigma-A weighted) and resolution/grid for small molecules.
- Define the batching/throughput contract with the GARNET training loop.
