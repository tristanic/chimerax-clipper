# ChimeraX-Clipper

A ChimeraX plugin providing macromolecular crystallographic data handling: electron density maps, structure factors, and crystallographic symmetry. It is the primary crystallographic infrastructure layer for [ISOLDE](https://isolde.cimr.cam.ac.uk/), which lives in a separate, tightly-coupled repository. Changes to public API here frequently require corresponding changes in ISOLDE.

Built on Kevin Cowtan's [Clipper library](http://www.ysbl.york.ac.uk/~cowtan/clipper/), wrapped via pybind11. Licensed LGPLv3+.

## Architecture

```
Python (src/)
  └─ chimerax.clipper package
       ├─ crystal.py, symmetry.py, map_calc.py   # core logic
       ├─ maps/          # map management & visualisation
       ├─ io/            # MTZ, CIF, CCP4 I/O
       ├─ reflection_tools/  # R-free, French-Wilson
       ├─ geometry/, graphics/, ui/
       └─ clipper_python  ← pybind11 bindings (below)

C++ (src_cpp/)
  ├─ bindings/           # pybind11 wrap_*.cpp files
  ├─ clipper_ext/        # ChimeraX-specific Clipper extensions
  ├─ contour/            # threaded contour generation
  ├─ _maps/, symmetry/, util/
  └─ deps/               # Clipper, libccp4, mmdb2 source trees

extern/
  ├─ pybind11            # git submodule
  └─ pocketfft           # git submodule (header-only FFT)
```

Compiled extensions (all pybind11):
- `clipper_python` — main Clipper bindings (~40 wrap_*.cpp files)
- `_symmetry`, `_util`, `_map_mask`, `contour_thread`

## Build

### Windows (primary platform)

`make_win.bat` self-initialises the Visual Studio 2022 build environment (the
v141 toolset that matches ChimeraX), so it can be run directly from a plain
`cmd` *or* PowerShell prompt — no need to pre-open a vcvars terminal:

```bat
# Install into the ChimeraX daily build
make_win.bat app-install

# Install into the stable ChimeraX release
make_win.bat release app-install

# Full rebuild from scratch (required after any C++ changes)
make_win.bat clean app-install
make_win.bat release clean app-install
```

ChimeraX paths assumed by `make_win.bat`:
- Daily:   `C:\Program Files\ChimeraX-Daily\bin\ChimeraX-console.exe`
- Release: `C:\Program Files\ChimeraX\bin\ChimeraX-console.exe`

### Linux / macOS

```sh
make app-install          # against daily ChimeraX build
make RELEASE=1 app-install  # against stable release
make clean app-install    # full rebuild
```

### When to rebuild from scratch

Always run `clean` before `app-install` after modifying any C++ source (`.cpp` or `.h` files under `src_cpp/`). Pure Python changes do not require a clean rebuild.

## Testing

There is no automated test suite used in practice. Testing is manual and interactive: install into ChimeraX, open the application, and exercise the relevant functionality. The `tests/` directory exists but is not part of the regular development cycle.

## Key files

| File | Purpose |
|------|---------|
| `src/symmetry.py` | Symmetry manager — largest and most central Python file (82 KB) |
| `src/crystal.py` | Crystal data structures (67 KB) |
| `src/maps/xmapset.py` | Real-space map set management |
| `src/io/mtz_read.py` | MTZ structure factor file reading |
| `src_cpp/clipper_ext/xtal_mgr.cpp` | Crystal manager C++ implementation |
| `src_cpp/bindings/clipper_pybind11.cpp` | Main pybind11 binding module entry point |
| `pyproject.toml` | Build configuration; defines all C extensions and their source lists |

## Code conventions

- **Match the existing style** — do not introduce new conventions (type hints, docstring formats, etc.) into files that don't already use them.
- **No comments unless the why is non-obvious** — the codebase is lightly commented by design.
- Domain-level crystallographic concepts (space groups, structure factors, HKL data, R-free, etc.) can be mentioned without explanation — brief one-liner context is welcome where it aids code reasoning.

## Critical gotchas

### 1. numpy/C++ array ownership
When passing numpy arrays into C++ extensions, never assume the C++ side owns the data. Lifetime of the underlying buffer must be managed on the Python side. This is a common source of subtle memory bugs in the pybind11 wrappers.

### 2. Do not duplicate ChimeraX functionality
If ChimeraX already provides something (atomic data access, file I/O, coordinate handling), use the ChimeraX API — don't re-implement it here. (Exception: the vendored MMDB2 source under `src_cpp/deps/mmdb2/` is retained because Clipper's `CIFFile` class depends on it and at least one third-party downstream uses `CIFFile` — so although ChimeraX-Clipper itself doesn't use MMDB2, it stays.)

### 3. Public API stability vs. ISOLDE
This plugin and ISOLDE are co-evolved. Removing or renaming any symbol that ISOLDE imports will break ISOLDE. When modifying the public interface (anything importable from `chimerax.clipper`), check whether ISOLDE uses it before removing it.

### 4. C++ standard is C++11
ChimeraX's build system compiles all C++ sources with `/std:c++11` (MSVC) or
`-std=c++11` (GCC/Clang). C++14/17 features that will silently compile on some
toolchains but hard-fail here include:
- Structured bindings: `auto [a, b] = f();` → use `auto p = f(); p.first; p.second`
- `if constexpr`, `std::string_view`, fold expressions, etc.

### 5. C++ changes require a clean rebuild
The build system does not always detect incremental C++ changes correctly. When in doubt, `clean` first.

### 6. `chimerax_bridge` — non-template functions live in the `.cpp`
The non-template `clipper_atoms_from_cx_atoms` family (incl.
`clipper_atoms_from_cx_atoms_with_map`) is **defined in
`src_cpp/clipper_ext/chimerax_bridge.cpp`** (compiled into `clipper_cx.dll`) and
declared with `CLIPPER_CX_IMEX` in `chimerax_bridge.h`. Every pybind11 translation
unit (e.g. `wrap_xtal_mgr.cpp`, `wrap_adp_occ_refiner.cpp`) therefore links to the
single exported definition — no ODR violations, no LNK2001, and a clean linker
exit code even with `/FORCE:MULTIPLE`. Only the *template* helper
`cl_atom_from_cx_atom` remains `inline` in the header (templates must be
header-visible). Do **not** move the non-template definitions back into the
header.

### 7. Windows DLL export warnings (C4251)
Exporting classes with private `std::future`, `std::unique_ptr`, or other STL members via
`CLIPPER_CX_IMEX` generates C4251 warnings on MSVC. These are suppressed with
`#pragma warning(disable: 4251)` at the top of the relevant `.cpp` file. The warnings are
cosmetic — the code is correct because the STL members are private implementation details.

## Dependencies

Runtime (key):
- `ChimeraX-Core ~1.11`, `ChimeraX-Atomic ~1.61`, `ChimeraX-AtomicLibrary ~14.2`
- `numpy ~1.26.4`

Bundled C++ (source-compiled):
- Clipper (`src_cpp/deps/clipper/`) — Kevin Cowtan's crystallography library
- libccp4 (`src_cpp/deps/libccp4/`) — CCP4 format support
- MMDB2 (`src_cpp/deps/mmdb2/`) — CCP4 macromolecular database library; used only by Clipper's `CIFFile` class (see gotcha #2)

Git submodules:
- `extern/pybind11` — C++/Python bindings
- `extern/pocketfft` — header-only FFT
- `extern/eigen` — Eigen3 (MPL2, header-only); **pin to ≥ 5.0.1** — Eigen master nightly has a
  regression where `has_unary_operator`/`has_binary_operator` SFINAE gives false positives on
  MSVC for `scalar_zero_op`, causing C2064 (`__forceinline` operator()() accepted with extra args)
- `extern/lbfgspp` — LBFGSpp L-BFGS-B (MIT, header-only); depends on Eigen

## Version bumping

Version is defined in `src/__init__.py`. The `pyproject.toml` reads it from there. Bump before building a Toolshed distribution.
