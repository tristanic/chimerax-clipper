# chimerax-clipper
Clipper plugin to ChimeraX, for handling of crystallographic maps and symmetry

A modified and extended version of [Clipper-Python](https://github.com/clipper-python/clipper-python) based on [Kevin Cowtan's Clipper Library](http://www.ysbl.york.ac.uk/~cowtan/clipper/) for use as a plugin to [UCSF ChimeraX](http://preview.cgl.ucsf.edu/chimerax/). Primarily built for use with [ISOLDE](https://isolde.cimr.cam.ac.uk/), but other uses are welcomed.

## Building
First of all, make sure you need to! Release builds of ChimeraX-Clipper can be installed directly from ChimeraX's ToolShed (Tools/More Tools...), and development builds are released semi-regularly. If you still want to build your own:

- Download and install a compatible version of ChimeraX (see the "ChimeraX-Core" entry in the *Dependencies* section of https://github.com/tristanic/chimerax-clipper/blob/master/bundle_info.xml for the current version restrictions).
- Make sure you have the correct compiler: for Linux, you'll need to build using the same toolchain as was used to build ChimeraX itself (currently `gcc-toolset-10` in RHEL 8). A Singularity recipe for their build environment can be found at https://github.com/RBVI/ChimeraX/blob/develop/prereqs/linux_buildenv/rhel-8.def. For MacOS you'll need XCode, and for Windows you'll need Visual Studio (ideally 2019, but newer versions tend to be good about maintaining backward compatibility).
- In Linux or MacOS, setting the environment variable RELEASE=1 will cause the Makefile to use the `chimerax` executable in Linux (assumed to be on the path), and `/Applications/ChimeraX.app/Contents/bin/ChimeraX in MacOS`. If the RELEASE environment variable is not set, then it will instead use `chimerax-daily` in Linux, and `/Applications/ChimeraX_Daily.app/Contents/bin/ChimeraX` in MacOS.
```
make clean
make app-install
```
- Since the documentation is generated using Sphinx, it needs to be generated *after* you've already installed the bundle once. 
  To regenerate and install the documentation, run:
```
make docs
make app-install
```
- In windows, "make" is only available in a CygWin or MinGW environment. To build from the Windows console, use the make_win.bat
  and make_docs.bat batch files:
```
make_win {release} clean app-install
make_docs {release}
make_win {release} app-install
```
If "release" is included, then building will use `c:\Program Files\ChimeraX\bin\ChimeraX-console.exe`, otherwise it will use `c:\Program Files\ChimeraX-Daily\bin\ChimeraX-console.exe`.
  Restart ChimeraX before using the plugin.
