# ChimeraX `corecif` fractional→Cartesian bug (oblique cells)

**Where:** `ChimeraX/src/bundles/mmcif/mmcif_cpp/corecif.cpp`, `SmallMolecule::to_cartesian` (~line 550).

**Symptom:** Small-molecule (core) CIF structures opened with `open … format corecif`
have subtly wrong atomic coordinates whenever the unit cell has a non-90° angle.
Bond lengths deviate from the CIF's own published `_geom_bond_distance` values by
~0.006 Å mean / 0.017 Å max for a β≈91° monoclinic cell (COD 1100908), and the
error grows as the cell becomes more oblique (small β, or triclinic). Orthogonal
cells (α=β=γ=90°) are unaffected.

**Root cause:** `compute_cell_matrix()` stores the cell vectors **a, b, c as the
rows** of `cell[3][3]`:
```
cell[0] = (a, 0, 0)
cell[1] = (b·cosγ, b·sinγ, 0)
cell[2] = (c·cosβ, c·n2, c·√(sin²β − n2²))
```
The correct fractional→Cartesian transform is therefore
`xyz[j] = Σ_i fract[i]·cell[i][j]`. But `to_cartesian` writes:
```cpp
for (auto i = 0; i < 3; ++i)
    for (auto j = 0; j < 3; ++j)
        xyz[i] += fract_xyz[i] * cell[i][j];   // BUG: output index should be j
```
i.e. `xyz[i] = fract[i]·Σ_j cell[i][j]` — a row-sum written into the wrong Cartesian
axis. For a diagonal (orthogonal) `cell` this coincides with the correct result,
which is why the bug is invisible for orthogonal cells and was never caught.

**Fix:**
```cpp
for (auto i = 0; i < 3; ++i)
    for (auto j = 0; j < 3; ++j)
        xyz[j] += fract_xyz[i] * cell[i][j];
```

**Verification:** rebuilding Cartesian coordinates from the CIF fractionals via a
correct orthogonalization (Clipper `Coord_frac::coord_orth`) reproduces the CIF's
published `_geom_bond_distance` values to max 0.0007 Å (mean 0.0002 Å) on COD
1100908, vs corecif's 0.017 Å max (0.006 Å mean).
