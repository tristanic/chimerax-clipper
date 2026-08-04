# The "make" targets are:
# 	wheel: build a Python wheel in "dist" directory.
# 	app-install: build wheel (if needed) and install in ChimeraX.
# 	test: run ChimeraX
# 	debug: run ChimeraX with debugging flag set
# 	clean: remove files used in building wheel
# 	distclean: remove files used in building wheel and license file

# These parameters may be changed as needed.

# ChimeraX bundle names must start with "ChimeraX_"
# to avoid clashes with package names in pypi.python.org.
# When uploaded to the ChimeraX toolshed, the bundle
# will be displayed without the ChimeraX- prefix.
# BUNDLE_NAME = ChimeraX-Clipper
# BUNDLE_VERSION = 0.9.4-dev1
# ChimeraX bundles should only include packages
# that install as chimerax.package_name.
# General Python packages should be uploaded to
# pypi.python.org rather than the ChimeraX toolshed.
# PKG_NAME = chimerax.clipper

# Define where ChimeraX is installed.
OS = $(patsubst CYGWIN_NT%,CYGWIN_NT,$(shell uname -s))
# CHIMERAX_APP is the ChimeraX install folder

ifeq ($(OS),CYGWIN_NT)
ifndef RELEASE
# Windows
CHIMERAX_APP = "/c/Program Files/ChimeraX_Daily"
else
CHIMERAX_APP = "/c/Program Files/ChimeraX"
endif
endif

ifeq ($(OS),Darwin)
# Mac
ifndef RELEASE
CHIMERAX_APP = /Applications/ChimeraX_Daily.app
else
CHIMERAX_APP = /Applications/ChimeraX.app
endif
endif

ifeq ($(OS),Linux)
ifndef RELEASE
CHIMERAX_APP = chimerax-daily
else
CHIMERAX_APP = chimerax
endif
endif

# ==================================================================
# Theoretically, no changes are needed below this line

# Platform-dependent settings.  Should not need fixing.
# For Windows, we assume Cygwin is being used.
ifeq ($(OS),CYGWIN_NT)
CHIMERAX_EXE = $(CHIMERAX_APP)/bin/ChimeraX.exe
endif
ifeq ($(OS),Darwin)
CHIMERAX_EXE = $(CHIMERAX_APP)/Contents/bin/ChimeraX
export MACOSX_DEPLOYMENT_TARGET=10.13
# Route compilation through ccache when available. distutils' customize_compiler
# honours CC/CXX for both the (serial) static-library phase and the (parallel)
# extension phase, so unchanged vendored sources become cache hits on rebuild.
CCACHE := $(shell command -v ccache 2>/dev/null)
ifdef CCACHE
export CC := ccache clang
export CXX := ccache clang++
endif
endif
ifeq ($(OS),Linux)
CHIMERAX_EXE = $(CHIMERAX_APP)
# Same ccache routing as macOS. cc/c++ are the generic driver names present on
# essentially all Linux toolchains, so this works whether the system compiler is
# gcc or clang.
CCACHE := $(shell command -v ccache 2>/dev/null)
ifdef CCACHE
export CC := ccache cc
export CXX := ccache c++
endif
endif

# --- Per-lane test isolation -------------------------------------------------
# Route every ChimeraX launch through run_chimerax.sh so each worktree/lane
# installs into (and tests from) its own isolated ChimeraX user directory --
# see run_chimerax.sh and tools/_isolated_chimerax.py. The CHIMERAX_APP block
# above is retained for its side effects (e.g. MACOSX_DEPLOYMENT_TARGET); this
# last assignment of CHIMERAX_EXE wins. On Windows use make_win.bat instead.
REL_TOK = $(if $(RELEASE),release,)
CHIMERAX_EXE = $(dir $(lastword $(MAKEFILE_LIST)))run_chimerax.sh $(REL_TOK)

BUNDLE_BASE_NAME = $(subst ChimeraX-,,$(BUNDLE_NAME))
SOURCE = src
SRCS = $(SOURCE)/*.py


:DEFAULT_GOAL := wheel

#
# Actual make dependencies!
#

wheel $(WHEEL): pyproject.toml
	$(CHIMERAX_EXE) --nogui --safemode --cmd "devel build . ; exit"

install app-install:	$(WHEEL)
	$(CHIMERAX_EXE) --nogui --safemode --cmd "devel install . ; exit"

uninstall app-uninstall:	$(WHEEL)
	$(CHIMERAX_EXE) --nogui --safemode --cmd "toolshed uninstall $(BUNDLE_BASE_NAME) ; exit"

docs:
	$(CHIMERAX_EXE) -m sphinx docs/source src/docs/user

test:
	$(CHIMERAX_EXE)

debug:
	$(CHIMERAX_EXE) --debug

clean:
	$(CHIMERAX_EXE) --nogui --safemode --cmd "devel clean . ; exit"

# Install a user-space startup hook that parallelises ChimeraX's serial
# static-library compile loop (the macOS build bottleneck). It lives in the
# ChimeraX user site-packages -- no write access to the .app bundle needed --
# and patches the live bundle_builder at start-up, so it survives daily updates
# with no re-run. See tools/chimerax_clipper_parallel_build.py.
install-parallel-build-hook:
	CXC_HOOK_SRC="$(CURDIR)/tools/chimerax_clipper_parallel_build.py" \
	  $(CHIMERAX_EXE) --nogui --silent --exit \
	  --script "$(CURDIR)/tools/install_parallel_build_hook.py"

uninstall-parallel-build-hook:
	CXC_HOOK_UNINSTALL=1 \
	  $(CHIMERAX_EXE) --nogui --silent --exit \
	  --script "$(CURDIR)/tools/install_parallel_build_hook.py"

.PHONY: docs install-parallel-build-hook uninstall-parallel-build-hook
