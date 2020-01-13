# chimerax-clipper
Clipper plugin to ChimeraX, for handling of crystallographic maps and symmetry

A modified and extended version of [Clipper-Python](https://github.com/clipper-python/clipper-python) based on [Kevin Cowtan's Clipper Library](http://www.ysbl.york.ac.uk/~cowtan/clipper/) for use as a plugin to [UCSF ChimeraX](http://preview.cgl.ucsf.edu/chimerax/). Primarily built for use with [ISOLDE](https://isolde.cimr.cam.ac.uk/), but other uses are welcomed.

## Building
First of all, make sure you need to! Release builds of ChimeraX-Clipper can be installed directly from ChimeraX's ToolShed (Tools/More Tools...), and development builds will start being made available in the same manner in the near future. If you still want to build your own:

- install [PyBind11](https://pybind11.readthedocs.io/en/stable/) - that is, make sure its headers are on your system path.
- Download and install a recent *daily build* of ChimeraX (ChimeraX-Clipper is *not* guaranteed to work with older versions). For Linux, either of the "RedHat 7" or "Generic Linux" build is recommended.
- Make sure you have the correct compiler: for Linux with either of the above builds, you'll need GCC 4.9 (available in RedHat devtoolset-3). For Mac, you'll need to install Xcode. For Windows, you'll need Visual Studio 2015). 
- In Linux or MacOS, assuming ChimeraX is installed in the default location (/opt/UCSF/ChimeraX-daily for Linux, /Applications/ChimeraX.app for MacOS) run:
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
  and make_docs.bat batch files (you may need to edit these to set the path to the ChimeraX executable):
```
make_win clean app-install
make_docs
make_win app-install
```
  Restart ChimeraX before using the plugin.
