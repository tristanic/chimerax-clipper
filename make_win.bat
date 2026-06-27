@echo off
setlocal EnableExtensions EnableDelayedExpansion

REM ChimeraX is built with the v141 (VS2017) toolset; initialise the matching
REM Visual Studio environment if it isn't already active.  Because a .bat always
REM runs under cmd.exe, this works whether invoked from cmd or PowerShell, so a
REM pre-opened vcvars terminal is no longer required.  Guarded on VCINSTALLDIR so
REM running from an existing vcvars terminal is a harmless no-op.
REM
REM Delayed expansion (!VSWHERE! / !VSPATH!) is used inside the if-block because
REM the resolved paths contain "(x86)" -- a literal ")" that would close the
REM block prematurely if expanded with %...% at parse time.
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not defined VCINSTALLDIR (
    if not exist "!VSWHERE!" (
        echo ERROR: vswhere.exe not found - is Visual Studio 2022 installed?
        exit /b 1
    )
    set "VSPATH="
    for /f "usebackq delims=" %%i in (`"!VSWHERE!" -latest -products * -property installationPath`) do set "VSPATH=%%i"
    if not defined VSPATH (
        echo ERROR: no Visual Studio installation found by vswhere.
        exit /b 1
    )
    call "!VSPATH!\VC\Auxiliary\Build\vcvars64.bat" -vcvars_ver=14.1
)
set DISTUTILS_USE_SDK=1
set MSSdk=1

REM ChimeraX launches go through run_chimerax.bat so each worktree/lane installs
REM into its own isolated ChimeraX user directory (see run_chimerax.bat). The
REM "release" token (stable vs daily install) is forwarded to it; "clean",
REM "app-install" and "install-hook" are consumed here.
set "RELEASE_TOK="
set DO_CLEAN=
set DO_INSTALL=
set DO_INSTALL_HOOK=
for %%A in (%*) do (
    if /i "%%A"=="release"      set "RELEASE_TOK=release"
    if /i "%%A"=="clean"        set DO_CLEAN=1
    if /i "%%A"=="app-install"  set DO_INSTALL=1
    if /i "%%A"=="install-hook" set DO_INSTALL_HOOK=1
)

REM Install the user-space build-speedup hook (routes cl.exe through sccache).
REM Re-run after a ChimeraX update; it installs into the ChimeraX user site dir.
REM Routed through run_chimerax.bat so it lands in the same (lane-isolated or
REM standard) user dir the build/install below will use.
if defined DO_INSTALL_HOOK (
    set "CXC_HOOK_SRC=%CD%\tools\chimerax_clipper_parallel_build.py"
    call "%~dp0run_chimerax.bat" %RELEASE_TOK% --nogui --silent --exit --script "%CD%\tools\install_parallel_build_hook.py"
)

if defined DO_CLEAN   call "%~dp0run_chimerax.bat" %RELEASE_TOK% --nogui --safemode --exit --cmd "devel clean ."
if defined DO_INSTALL call "%~dp0run_chimerax.bat" %RELEASE_TOK% --nogui --safemode --exit --cmd "devel install ."

endlocal
