ChimeraX-Clipper small-molecule example data
============================================

These are example small-molecule structures from the Crystallography Open
Database (COD, https://www.crystallography.net), bundled for trying out the
`clipper smallmol` command offline. Each entry is a model (.cif, core-CIF
dialect) plus its experimental structure factors (.hkl, a CIF with a _refln_
loop of F_squared_meas):

  cod_1100908  C 1 2/c 1 (#15)  - a Cu complex with the metal on a 2-fold
                                  special position; published R = 0.041.
  cod_2213867  P b c a   (#61)  - an orthorhombic organic structure with
                                  glide planes.

Open one in ChimeraX (model + structure factors + live 2Fo-Fc / Fo-Fc maps)
with, e.g.:

  clipper smallmol <this-folder>/cod_1100908.cif hkl <this-folder>/cod_1100908.hkl

To print the path of this folder from the ChimeraX Python shell (Tools >
General > Shell):

  import os, chimerax.clipper
  print(os.path.join(os.path.dirname(chimerax.clipper.__file__), 'demo'))

Or fetch either entry (and any other COD entry with structure factors) directly
from the database:

  clipper cod 1100908
