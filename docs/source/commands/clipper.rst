ChimeraX-Clipper Commands
=========================

.. _`open`:

clipper open
------------

Syntax: clipper open *path* [**structureModel** *structure*]
[**overSampling** *number*]

Open a structure factor file in .mtz or .cif format, and generate maps for the
model specified with *structureModel*.

If the structure factor file contains pre-calculated amplitudes and phases, they
should appear in the file in strict (Amplitudes, Phases, Amplitudes, Phases)
order. ChimeraX-Clipper will create one map for each set it finds, and prepend
'(STATIC)' to its name. To close any unwanted maps, look in the ChimeraX Models
panel under "Data manager"/"Map Manager"/"Crystallographic maps", select the
maps you wish to remove and click "Close":

.. image:: images/close_maps.png

Experimentally-measured reflections may be provided as amplitudes (F/sigF),
intensities (I/sigI), or their anomalous counterparts (F+/sigF+/F-/sigF- or
I+/sigI+/I-/sigI-) and should be stored in the file in these orders. Only one
set of experimental data will be used, with intensities preferred. If
intensities are provided, they will be internally converted to amplitudes using
the analytical French & Wilson method of `Read et al.`_. Anomalous datasets will
be merged.

.. _Read et al.: https://journals.iucr.org/d/issues/2016/03/00/dz5382/index.html

When experimental reflections are found, three live-updating maps will be
created (with "(LIVE)" prepended to their names). The first of these is a
standard 2mFo-DFc map, while the second is a 2mFo-DFc map with a degree of
B-factor sharpening or smoothing applied. The level of sharpening or smoothing
depends on the resolution of the data: maps with resolutions worse than 2.5Å
will be sharpened, while higher-resolution maps will be smoothed. The sharper
of the two maps will be displayed as a transparent surface, the other as a
wireframe. Finally, a standard mFo-DFc map will be generated and displayed
with contours at ± 3 sigma.

.. _`associate`:

clipper associate
-----------------

Syntax: clipper associate *volumes* [toModel *structure*]

Have ChimeraX-Clipper take control of the chosen volumes and associate them with
the given model.

In order to work with ChimeraX-Clipper's visualisation modes, a volumetric map
(such as a cryo-EM map) must first be associated with an atomic model using this
command.

.. _`spotlight`:

clipper spotlight
-----------------

Syntax: clipper spotlight [*structures*] [**radius** *number*]

Initiate "spotlight mode" (a sphere of visible atoms and density following the
centre of rotation) for the given models, and optionally set the radius of the
sphere. If *structures* is not specified, the command will only apply to models
which are already initialised into the ChimeraX-Clipper data structure.

.. _`isolate`:

clipper isolate
---------------

Syntax: clipper isolate *atoms* [**surroundDistance** *number* (5.0)]
[**contextDistance** *number* (5.0)] [**maskRadius** *number* (3.0)]
[**hideSurrounds** *true/false* (true)] [**focus** *true/false* (false)]
[**includeSymmetry** *true/false* (true)]

Visually isolate the selected atoms from their surroundings, and mask all
associated maps to cover the immediate vicinity. The algorithm for deciding
shown and
