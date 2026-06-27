"""ChimeraX startup hook: parallelise bundle_builder's serial static-library compile.

ChimeraX's bundle builder compiles a bundle's static C *libraries* (for us: mmdb2,
ccp4, clipper, clipper_cx -- ~110 source files) one at a time on a single core,
via distutils' ``compiler.compile()`` in a plain loop. The pybind11 *extensions*
already build with ``-j``, but the library phase does not, and on a universal2
build every file compiles twice -- that serial phase dominates a macOS build.

This module is installed into the ChimeraX *user* site-packages alongside a tiny
``.pth`` file whose single line is ``import chimerax_clipper_parallel_build``.
ChimeraX runs ``site.addsitedir`` on the user site-packages during toolshed
start-up (even under ``--safemode``), which processes the ``.pth`` and imports
this module. We then monkeypatch the *live* ``bundle_builder`` so the
static-library phase compiles across all cores.

Design notes:
- Patches the live module, so it tracks whatever ``bundle_builder`` version the
  running ChimeraX ships -- no rebuild needed after a daily update.
- Scoped: only the static-library phase is parallelised; the extension phase
  (already parallel via setuptools' ``-j``) is left untouched.
- Robust: structural guards make it a silent no-op if the code shape ever
  changes, and the whole thing is wrapped so a failure can never break start-up.
- Set ``CXC_PARALLEL_BUILD_DEBUG=1`` to print a line when the parallel path runs.

Cross-platform notes:
- Parallel compile (``CCompiler.compile``) covers macOS and Linux, since
  ``UnixCCompiler`` inherits that method. Windows ``MSVCCompiler`` overrides
  ``compile`` with its own loop and is intentionally not parallelised here.
- Compiler caching is split by platform: on Unix the Makefile wires ccache via
  ``CC``/``CXX``; on Windows this module routes ``cl.exe`` through sccache by
  wrapping ``MSVCCompiler.spawn`` (guarded by sccache being on PATH).
"""


def _install_parallel_compile_hook():
    import os
    import distutils.ccompiler
    from multiprocessing.pool import ThreadPool
    import chimerax.bundle_builder.bundle_builder_toml as bbt

    Bundle = bbt.Bundle
    method_name = "_clear_distutils_dir_and_prep_srcdir"

    # Structural guards: only patch when the code is shaped as we expect.
    if not hasattr(Bundle, method_name):
        return False
    for internal in ("_setup_compile", "_get_cc_args", "_compile"):
        if not hasattr(distutils.ccompiler.CCompiler, internal):
            return False
    if getattr(Bundle, "_cxc_parallel_build_patched", False):
        return True  # already patched in this interpreter

    def _parallel_compile(self, sources, output_dir=None, macros=None,
                          include_dirs=None, debug=0, extra_preargs=None,
                          extra_postargs=None, depends=None):
        # Mirrors distutils.ccompiler.CCompiler.compile, but runs the per-file
        # compile across a thread pool. Reuses distutils' own internals so the
        # result is byte-for-byte identical to the serial path.
        macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
            output_dir, macros, include_dirs, sources, depends, extra_postargs)
        cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)

        def _one(obj):
            try:
                src, ext = build[obj]
            except KeyError:
                return
            self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

        ThreadPool(os.cpu_count() or 1).map(_one, objects)
        return objects

    original = getattr(Bundle, method_name)

    def _wrapped(self, build_exts=False):
        # Only the static-library phase (build_exts=True) is serial. Swap in the
        # parallel compiler.compile just for its duration, then restore so the
        # setuptools-driven extension phase is unaffected.
        if not build_exts:
            return original(self, build_exts)
        saved = distutils.ccompiler.CCompiler.compile
        distutils.ccompiler.CCompiler.compile = _parallel_compile
        if os.environ.get("CXC_PARALLEL_BUILD_DEBUG"):
            import sys
            sys.stderr.write(
                "[chimerax-clipper] parallel static-library compile active "
                "(%d workers)\n" % (os.cpu_count() or 1))
        try:
            return original(self, build_exts)
        finally:
            distutils.ccompiler.CCompiler.compile = saved

    Bundle._clear_distutils_dir_and_prep_srcdir = _wrapped
    Bundle._cxc_parallel_build_patched = True
    return True


def _install_compiler_cache_hook():
    # Windows only: route MSVC cl.exe compiles through sccache (the ccache
    # equivalent that supports MSVC). On Unix the Makefile wires ccache via
    # CC/CXX instead, so this is a no-op there -- the cl.exe check never matches.
    import shutil
    if not shutil.which("sccache"):
        return False
    try:
        from distutils._msvccompiler import MSVCCompiler
    except Exception:
        return False
    if getattr(MSVCCompiler, "_cxc_sccache_patched", False):
        return True

    original_spawn = MSVCCompiler.spawn

    def _sccache_spawn(self, cmd, *args, **kwargs):
        # Only wrap cl.exe (the compiler). Linking/lib steps and every Unix
        # command fall through untouched; sccache itself ignores non-compiles.
        if cmd and str(cmd[0]).lower().endswith("cl.exe"):
            cmd = ["sccache", *cmd]
        # Call through to MSVC's own spawn so its PATH/env setup and
        # long-command-line response-file fallback are preserved.
        return original_spawn(self, cmd, *args, **kwargs)

    MSVCCompiler.spawn = _sccache_spawn
    MSVCCompiler._cxc_sccache_patched = True
    return True


for _hook in (_install_parallel_compile_hook, _install_compiler_cache_hook):
    try:
        _hook()
    except Exception:
        # A build optimisation must never prevent ChimeraX from starting.
        pass
