ChimeraX-Clipper: electron density and symmetry tools for ChimeraX
==================================================================

Contents:

.. toctree::
   :maxdepth: 2

   Home <self>
   commands/clipper
   api/clipper


Introduction
------------

ChimeraX-Clipper makes use of a modified version of the `Clipper`_ library from
Kevin Cowtan to add support for crystallographic symmetry and structure factor
calculations to ChimeraX. In addition, it adds visualisation modes suited to
macromolecular model building, and a convenience API to provide a unified
environment handling both cryo-EM and crystallographic data. It is a core
library required by `ISOLDE`_. All source code is available on `GitHub`_ and
is licensed under the GNU Lesser General Public License v3.0.

.. _Clipper: http://www.ysbl.york.ac.uk/~cowtan/clipper/doc/
.. _ISOLDE: https://isolde.cimr.cam.ac.uk
.. _GitHub: https://github.com/tristanic/chimerax-clipper

ChimeraX-Clipper consists of three essential parts:

1. A near-complete wrapping of the core Clipper libraries (modified for
   multithreading in key routines) using `PyBind11`_. Documentation for this
   layer is currently minimal - but since the bindings are for the most part a
   direct reflection of the underlying C++ objects the `Clipper API documentation`_
   is quite useful.
2. A C++ extension layer (also wrapped using PyBind11) built on the core API,
   to implement key tasks which need to run in separate threads to avoid
   disruption to the ChimeraX GUI. In particular, this code handles the
   parallelised, background structure factor calculations used to generate new
   crystallographic maps on changes to the model. Documentation of this layer
   is also currently minimal.
3. A ChimeraX-specific Python layer designed to facilitate quick scripting. This
   component *is* (lightly) documented here.

.. _PyBind11: https://github.com/pybind/pybind11
.. _Clipper API documentation: http://www.ysbl.york.ac.uk/~cowtan/clipper/doc/annotated.html



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
