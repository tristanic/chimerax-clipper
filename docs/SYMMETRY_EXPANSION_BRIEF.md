# Invertible symmetry expansion — capability brief

Status: **implemented on `master`** (headless-verified; commit pending at time of
writing). Audience: the *Clipper B-factor refinement extensions* agent preparing the
ML-training pipeline. This is the capability you were waiting on: build a full P1 box
for the MD engine and fold it back onto a single ASU for the (differentiable)
structure-factor calculation.

Related: `SYMMETRY_EXPANSION_BRIEF` complements the differentiable-targets work
(`src/diff/`) and the GARNET scoring note (`MODE2_SCORING_SPEC.md`).

## TL;DR

`clipper symcopies` / `clipper unitcells` now (a) run headless, (b) de-duplicate
special positions **exactly** (Clipper site multiplicity, not a geometric heuristic),
and (c) attach an **invertible** `clipper_sym_expansion` record to the result. A new
`collapse_to_asu()` inverts the expansion. Everything is importable from
`chimerax.clipper`.

## The contract you rely on (P1 box ⇄ ASU)

1. Expand an ASU to a full P1 box: `clipper unitcells #m cells n,m,o` (or `box x,y,z`).
   The box holds `n_asu` copies with whole molecules and sensible connectivity.
2. Relax the box in your MD engine (PBC / minimum-image).
3. Collapse back for the SF calc: apply each copy's **inverse** operator, set
   occupancy `1/n_asu`, overlay in the ASU frame. `collapse_to_asu()` does this;
   your torch pipeline can reproduce it differentiably from the recorded operators.

**Why the occupancy is automatically correct:** collapsing `n_asu` copies at
`1/n_asu` sums a general-position atom to `n_asu·(1/n_asu)=1`. A special-position
atom of multiplicity `m`, after exact dedup, survives as only `n_asu/m` distinct
images, so it collapses to `(n_asu/m)·(1/n_asu)=1/m` — exactly Clipper's SFcalc
convention. This holds *because* dedup removes exactly the redundant `m−1` images.

## API (all from `chimerax.clipper`)

```python
from chimerax.clipper import (
    crystal_symmetry_for,      # (structure, refl_file=None) -> CrystalSymmetry
    unit_cell_places,          # (cell, spacegroup, na, nb, nc, origin=(0,0,0)) -> [Place]
    realize_symmetry_copies,   # (session, structure, places, name=None,
                               #  prune_special_positions=True, tolerance=0.5,
                               #  multiplicities=None) -> combined model | None
    collapse_to_asu,           # (session, model, name=None, set_occupancies=True) -> model
    site_multiplicities,       # (coords, cell, spacegroup, grid) -> int array
    SymmetryExpansion,
)
```

`model.clipper_sym_expansion` is a `SymmetryExpansion`:

| field | meaning |
|---|---|
| `operators` | `list[Place]`, structure-frame, **identity first** |
| `n_asu` | number of copies (`len(operators)`) |
| `chains_by_operator` | `list[list[str]]` — output chain IDs each operator produced |
| `operator_by_chain` | `{chain_id: operator_index}` |
| `inverse_operators` | `[Place.inverse() …]` (fold a copy back onto the ASU) |
| `deduplicated` | `[{operator, chain_id, residue, multiplicity}, …]` collapses done |
| `source_name` | source model name |

### Build a box, then collapse (reference)

```python
box = run(session, "clipper unitcells #1 cells 2,2,2")[0]     # or use the Python API
exp = box.clipper_sym_expansion
# ... MD relaxes box.atoms.coords ...
asu = collapse_to_asu(session, box)      # ASU-frame overlay, occ = 1/exp.n_asu
```

Differentiable path: don't call `collapse_to_asu`; read `exp.operators` /
`exp.inverse_operators` and `exp.chains_by_operator`, and apply the per-chain rigid
transform in torch (it's a linear map, so gradients flow straight through).

## Guarantees / semantics you can lean on

- **B-factors survive; anisotropic ADPs are handled for you.** `bfactor` and
  `occupancy` are preserved through `combine` and `copy()`, so your per-atom
  B-factors propagate box⇄ASU with no extra work. **`aniso_u6` is different** —
  ChimeraX's `copy()` *and* `combine` both silently drop it, so the expansion
  re-stamps ADPs itself, **rotated by the operator** (`U' = R U Rᵀ`, via Clipper's
  `U_aniso_orth.transform`), and `collapse_to_asu` rotates them back by the inverse.
  So expanded/collapsed models carry correctly-oriented ADPs — but if you build your
  own copies with `AtomicStructure.copy()`/`combine`, be aware you must re-apply
  `aniso_u6` yourself.
- **Unique chain IDs.** Every copy gets collision-free chain IDs, so the operator⇄
  chain mapping is exact and the collapse is unambiguous. (No "Remapping…" churn.)
- **`clipper_scattering_species` is preserved** onto the box and collapsed model
  (ionic scattering factors survive), so a SF calc on the collapsed ASU is correct.
- **`clipper_sf_exclude` is deliberately NOT carried.** It is a *single-ASU*
  artifact (it tells Clipper's ASU SFcalc to skip symmetry-completed atoms it would
  regenerate). After explicit expansion + dedup those atoms are distinct and counted,
  so the flag is obsolete — **your SF calc on the collapsed model must count every
  atom; do not re-introduce an exclusion.**
- **Coordinate frames.** Box atoms live in the structure frame (operators baked in);
  `collapse_to_asu` maps each copy back into the ASU frame.

## Preconditions & caveats

- **Complete special-position-straddling molecules first.** For clean dedup *and* MD
  connectivity, run `clipper fragments #m mode complete` before expanding (e.g. a
  water on a 2-fold with one modelled H). Expansion warns if it detects an
  un-completed straddling molecule.
- **`clipper fragments` now works on a mapped model.** (Previously it refused; that
  guard protected a scaffold index that no longer exists.) You can complete
  fragments before or after generating maps — the live small-molecule map rebuilds
  from `model.atoms` and reflects the split on its next recompute.
- **Anisotropic ADPs are handled in both directions** (see the guarantees above) —
  `collapse_to_asu` rotates `aniso_u6` by the inverse operator, so it is correct for
  anisotropic as well as isotropic refinement.
- **Site multiplicity is grid-quantised** (from `Xmap.multiplicity` at the cell's
  grid sampling). Exact for atoms on the special element; that's what `crystal_
  symmetry_for` supplies.
- **Box fill is per-axis** `⌈x/a⌉` — approximate for non-orthogonal cells.

## Verified (headless)

Provenance partitions every chain exactly once (`n_asu` matches symop count);
expand→perturb→`collapse_to_asu` round-trip shrinks a box back to one ASU at occ
`1/n_asu`; exact dedup records Cu-on-a-2-fold at multiplicity 2; scattering species
preserved through expansion; `sf_exclude` correctly absent; no chain-remap logging;
`bfactor`/`occupancy` preserved through `combine`; anisotropic ADPs rotated in the
box and restored exactly on collapse (no loss, `n_asu`×source count).

## Open coordination points

1. **SFcalc convention on the collapsed model** — the collapse yields an ASU-frame
   overlay of `n_asu` conformers at `1/n_asu`. Confirm your SF path expects that
   (original space group, occupancy-weighted) rather than a single conformer.
