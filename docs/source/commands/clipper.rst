.. _clipper_commands:

Commands
========

The primary purpose of ChimeraX-Clipper is to facilitate the handling of
crystallographic maps and symmetry in ChimeraX. However, it's not restricted to
crystallographic data: real-space volumetric maps of any kind may be associated
with your model to provide a unified visualisation scheme
(see :ref:`associate`). Using Clipper you can open (:ref:`open`) and save
(:ref:`save`) structure factors in MTZ format, and explore the detailed fit of
you model to maps (see :ref:`spotlight` and :ref:`isolate`). Symmetry-related
molecules can be turned into real, editable copies with :ref:`symcopies`.
Small-molecule
crystals (e.g. from the Crystallography Open Database) can be opened with live
symmetry and maps using :ref:`smallmol` and :ref:`cod`. Where applicable,
a model initialised for ChimeraX-Clipper will provide a real-time display of
crystallographic symmetry-related molecules. The appearance of live maps can be
tuned on the fly: adjust their sharpening/smoothing B-factor with :ref:`bsharp`,
or change the real-space oversampling (grid fineness) with :ref:`set oversampling`.

.. _`open`:

clipper open
------------

*(NOTE: MTZ files - but not .cif files - may also be opened using the standard
ChimeraX open command, with otherwise identical syntax to that described below)*

Syntax: clipper open *path* [**structureModel** *structure*]
[**overSampling** *number*] [**radiation** *xray/electron/auto* (auto)]

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

The **radiation** keyword selects the scattering-factor table used to calculate
structure factors from the model:

* ``xray`` (default) - X-ray form factors (Waasmaier & Kirfel), for conventional
  X-ray crystallography;
* ``electron`` - electron form factors (Peng, *International Tables* Vol C), for
  **electron diffraction (micro-ED / 3D-ED)** data. Electrons scatter from the
  electrostatic potential rather than the electron density, so the live maps and
  the reported R-factors are then correct for an electron-diffraction experiment.
  Any genuine monatomic ions in the model (e.g. an ordered ``Na+``, ``Ca2+`` or
  ``Cl-``) use the corresponding Peng-1998 ionic electron factors where tabulated,
  falling back to the neutral element otherwise;
* ``auto`` (default) - infer the radiation, in order of reliability, from the
  associated model's mmCIF ``_exptl.method`` (``ELECTRON CRYSTALLOGRAPHY`` vs
  ``X-RAY DIFFRACTION``), then from a structure-factor CIF's
  ``_diffrn_radiation_probe``; anything undeclared (e.g. a bare MTZ with no model
  metadata) resolves to X-ray. In practice this means electron-diffraction entries
  are detected automatically while conventional X-ray data is unaffected.

Because the model header is the primary signal, structure factors fetched from
the PDB are handled correctly with no extra input - e.g. for a micro-ED entry,
``open 8xyz structureFactors true`` auto-selects electron factors. An explicit
``radiation xray`` or ``radiation electron`` always overrides the auto-detection.

.. _`smallmol`:

clipper smallmol
----------------

Syntax: clipper smallmol *path* [**hkl** *path*]
[**radiation** *xray/electron/auto* (auto)]
[**fragments** *off/rename/complete* (rename)]

Open a small-molecule CIF (the "core CIF" dialect used by, e.g., the
Crystallography Open Database) as a live crystal structure: the model in its
unit cell, crystallographic symmetry mates, and - when experimental reflections
are available - live-updating electron-density maps. This is the small-molecule
counterpart to :ref:`open` (which handles macromolecular MTZ/structure-factor
data); use *smallmol* for single-molecule crystals where the unit cell and
symmetry come from the model's own CIF.

*path* is the model CIF. Reflections may be supplied in either of two ways:

* explicitly, with the **hkl** keyword, pointing at a small-molecule
  structure-factor file (itself a CIF, conventionally with a ``.hkl``
  extension, containing a ``_refln_`` loop of squared structure factors
  ``F_squared_meas``); or
* implicitly - if **hkl** is omitted, a sibling file alongside *path* with the
  same stem and a ``.hkl`` extension is used automatically when present.

If no reflections are found the model is still opened with its unit cell and
symmetry, but no maps are created.

When reflections are present, two live maps are generated (with "(LIVE)"
prepended to their names) and recomputed as the model moves:

* a **2mFo-DFc** map, contoured over the model, and
* an **mFo-DFc** difference map.

These follow the small-molecule conventions rather than the macromolecular ones
used by :ref:`open`. Structure factors are calculated by FFT with **no bulk
solvent** (small-molecule crystals are densely packed) and scaled onto the
observed amplitudes with an anisotropic-Gaussian-plus-isotropic-spline scaling;
there is no σA weighting and no free set. The difference map is contoured at an
**absolute level in e/Å³** (not the ± 3 σ used for macromolecular maps): the
familiar "3 σ" rule of thumb is calibrated for protein maps and is misleading at
small-molecule data quality, where a residual peak is best judged against an
absolute electron-density scale.

The **radiation** keyword selects the scattering-factor table used for the
calculated structure factors:

* ``xray`` (X-ray form factors; Waasmaier & Kirfel, including ionic species
  declared in the CIF ``_atom_site_type_symbol`` such as ``Cu2+`` / ``O2-``);
* ``electron`` (electron form factors; Peng 1996/1998, *International Tables*
  Vol C), for **electron diffraction — micro-ED / 3D-ED** data. Electron scattering
  senses the electrostatic potential, so the calculated maps and R-factor are
  correct for electron-diffraction experiments. Ionic species declared in the CIF
  (``Cu2+``, ``O2-``, ...) use Peng-1998 ionic electron factors, including the
  divergent Coulomb term of the ionic charge (the R-factor uses the full ionic
  factor; live maps use the screened part, since the Coulomb term is a long-range
  potential with no localised density).
* ``auto`` (default) reads ``_diffrn_radiation_probe`` from the CIF and uses
  electron factors when it names electrons, otherwise X-ray.

The **fragments** keyword controls how the asymmetric unit is divided into
residues (a core CIF carries no per-species annotation, so ChimeraX's parser
otherwise delivers the whole ASU as a single ``UNL`` residue). See
:ref:`fragments` for the full description; the default is ``rename``.

.. _`cod`:

clipper cod
-----------

Syntax: clipper cod *id* [**ignoreCache** *true/false* (false)]
[**radiation** *xray/electron/auto* (auto)]
[**fragments** *off/rename/complete* (rename)]

Fetch a structure from the Crystallography Open Database
(https://www.crystallography.net) by its numeric COD *id*, then open it with
:ref:`smallmol`. Both the model (``<id>.cif``) and, when the entry provides
them, its structure factors (``<id>.hkl``) are downloaded; roughly 10% of COD
entries include structure factors, and only those will produce live maps.

Downloads are cached locally; pass **ignoreCache true** to force a fresh fetch.
The **radiation** and **fragments** keywords are passed through to
:ref:`smallmol` (X-ray/electron/auto radiation; ``off``/``rename``/``complete``
fragment splitting).

For example, ``clipper cod 1100908`` fetches the Cu complex used as a bundled
example (see :ref:`smallmol_examples`).

.. _`fragments`:

clipper fragments
-----------------

Syntax: clipper fragments *structures* [**mode** *off/rename/complete* (rename)]

Split a small-molecule asymmetric unit into its individual covalent fragments,
each in its own residue with a sensible name. A core CIF encodes no per-species
identity, so ChimeraX's parser delivers the whole ASU as one ``UNL`` residue;
this command (also available as the **fragments** keyword of :ref:`smallmol` /
:ref:`cod`, applied automatically on opening) gives water, ions and each distinct
molecule their own residues. It operates on already-open small-molecule
(core-CIF) *structures*.

Fragments are the connected components of the covalent-bond graph, with two
crystallographic refinements:

* **Metals stand alone.** A bond to a metal atom is treated as a fragment
  boundary, so a coordinated water, counter-ion or ligand is separated from the
  metal it coordinates (the metal becomes its own residue). This is applied
  uniformly, including for metals embedded in a larger ligand (heme, chlorophyll):
  the metal is always its own residue.
* **Symmetry-split molecules are recognised.** A molecule can be only partly
  present in the ASU because it straddles a symmetry element - e.g. a water whose
  oxygen sits on a 2-fold axis, with only one of its hydrogens modelled, looks
  like a hydroxide. Such molecules are identified from the site symmetry and the
  CIF's own ``_geom_bond_site_symmetry_2`` records (which name every bond to a
  symmetry-equivalent atom), and named correctly.

Naming uses a small built-in table of the common simple species - water
(``HOH``), monatomic ions (``NA``, ``CL``, ``CU``, ...) and small inorganic ions
(``SO4``, ``PO4``, ``NO3``, ``CO3``, ...); everything else becomes ``LIG01``,
``LIG02``, ... .

The **mode** keyword selects:

* ``off`` - do nothing (leave the single ``UNL`` residue).
* ``rename`` (default) - split and name the fragments; the atom set is
  unchanged. A symmetry-split molecule is named for what it really is (the
  special-position water becomes ``HOH``, not hydroxide) but its missing
  symmetry-equivalent atoms are not added.
* ``complete`` - reassemble whole molecules through crystallographic symmetry, in
  two complementary ways. First, **rewrap symmetry-scattered hydrogens**: a hydrogen
  whose bond to its parent is symmetry-coded in the CIF (``_geom_bond``) is deposited
  in a neighbouring asymmetric unit and comes out unbonded; it is moved back onto its
  parent (each such H maps unambiguously under one symmetry operator). Second, **add
  the symmetry-generated atoms that finish split molecules** (e.g. the second hydrogen
  of a special-position water), correcting occupancies (atoms on a special position
  are given unit occupancy). Both move or add atoms - hence ``complete`` only, not the
  light-touch ``rename`` default. The *added* atoms (not the rewrapped hydrogens, which
  are real ASU atoms) are excluded from the live structure-factor calculation and do
  not follow live edits of their source atoms.

The command may be run whether or not the model already has a live small-molecule
map: the map rebuilds its structure factors from ``model.atoms`` on each recompute,
so it simply reflects the renamed (and, in ``complete`` mode, added) atoms on its
next update.

Opening a small-molecule CIF (this command, ``clipper smallmol`` and ``clipper cod``)
also logs a warning for any covalent bond far longer than physically possible - almost
always a spurious ``_geom_bond`` entry in the deposited CIF (e.g. a 1,3 distance listed
as a bond) that the corecif parser copies verbatim. Such bonds are **flagged, not
removed** (a long contact can be intentional); review the connectivity yourself.

.. _`smallmol_examples`:

Example data
------------

A couple of small-molecule examples (model + structure factors) are installed
alongside the bundle so the above commands can be tried offline:

* ``cod_1100908`` - a Cu complex in C 1 2/c 1 (#15), with the metal on a 2-fold
  special position (published R = 0.041);
* ``cod_2213867`` - an organic structure in P b c a (#61).

To locate them, run this in the ChimeraX Python shell (Tools → General → Shell)::

    import os, chimerax.clipper
    print(os.path.join(os.path.dirname(chimerax.clipper.__file__), 'demo'))

then, for example::

    clipper smallmol <that-folder>/cod_1100908.cif hkl <that-folder>/cod_1100908.hkl

.. _`save`:

clipper save
------------

*(NOTE: This command is also available via the top-level ChimeraX "save" command
with identical behaviour and syntax)*

Syntax: clipper save *filename* [*models*]
[**preserveInput** *true/false* (false)]
[**saveMapCoeffs** *true/false* (false)]

Save one or more sets of reflection data to (a) MTZ file(s). If more than one
dataset is specified by *models*, they will be saved as a series of files
(*filename*-0.mtz, *filename*-1.mtz, etc.) numbered in the order they are found
in the ChimeraX model tree.

If *preserveInput* is true, then the original experimental data loaded from file
will be saved, with "in." prepended to the column labels. **(IMPORTANT NOTE: if
the data was originally loaded from .cif, ONLY the columns selected by Clipper
for map calculations will be passed through)**

If *saveMapCoeffs* is true, amplitudes and phases for Clipper's current
Fc, 2Fo-Fc and Fo-Fc maps will be saved.

Free flags and the Fobs/SigFobs arrays used by Clipper will always be saved.


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

.. _`debugmaps`:

clipper debugmaps
-----------------

Syntax: clipper debugmaps [*structures*] [*enable* *true/false* (true)]

Add (or, with *enable false*, remove) a set of live "debug" maps exposing the
individual structure-factor components that go into the bulk-solvent model, as a
diagnostic aid:

* **Fcalc (total)** - the full model structure factors (atoms + bulk solvent).
* **Fatoms** - the atoms-only contribution (Fcalc minus the bulk term).
* **Fbulk** - the additive bulk-solvent contribution (k_sol scaled).
* **Fmask** - the (smoothed) solvent-mask transform; viewing this shows the
  bulk-solvent region directly, e.g. to inspect occupancy-weighted masking.

Each is a live map that updates with the model, just like the standard maps. They
are created hidden and are *display-gated*: while hidden they are not transformed
(so they cost nothing), and each is recomputed only while shown. Debug maps are
not saved in sessions; re-run this command to recreate them.

.. _`scaling`:

clipper scaling
---------------

Syntax: clipper scaling [*structures*] [*allReflections* *true/false*] [*reflectionsPerBin* *N*] [*numBins* *N*] [*seed* *N*]

Control how the live crystallographic maps fit the Fcalc-to-Fobs bulk-solvent scale
for the given model(s). By default the scale is fit on a fixed, *seeded subset* of
reflections (500 per bin across 20 resolution bins) - fast, and *reproducible*: the
subset is drawn once and reused across recalculations, so the R-factor and maps no
longer jitter between recalculations. It carries a tiny fixed bias versus the full
fit, which is negligible for interactive maps.

* **allReflections true** fits the scale over *all* reflections instead: exact and
  unbiased, at a larger (one-off) cost. Used internally where the exact scale matters
  (e.g. the B-factor/occupancy refiner); rarely needed interactively now that the
  default is already reproducible.
* **reflectionsPerBin** / **numBins** set the size of the seeded subset (defaults 500
  and 20). Exposed for tuning and characterisation.
* **seed** sets which subset is drawn (default 10061865). Vary it to check that the
  R-factor and maps are insensitive to *which* reflections were chosen.
* The subset knobs (reflectionsPerBin, numBins, seed) are ignored when
  *allReflections* is true.

With no options, reports the current settings. Applies to live maps only.

.. _`isolate`:

clipper isolate
---------------

Syntax: clipper isolate *atoms* [**surroundDistance** *number* (0.0)]
[**contextDistance** *number* (5.0)] [**maskRadius** *number* (3.0)]
[**hideSurrounds** *true/false* (true)] [**focus** *true/false* (false)]
[**includeSymmetry** *true/false* (true)]

Visually isolate the selected atoms from their surroundings, and mask all
associated maps to cover the immediate vicinity. The algorithm for deciding
the final view is as follows:

1. The *atoms* selection is expanded to complete residues.
2. All residues with any atoms within *surroundDistance* of any atoms in the
   results from (1) are added to the selection. If *includeSymmetry* is true,
   symmetry atoms within the *surroundDistance* will also be included. The final
   selection at this stage will be covered by the map(s).
3. *(This step only has an effect if hideSurrounds is true)* All residues with
   any atoms within *contextDistance* of the results from (2) will be displayed,
   but not covered by the map(s). If *includeSymmetry* is true, symmetry atoms
   will be included.
4. All maps associated with the model are masked down to within *maskRadius* of
   the atoms specified in (2).
5. If *hideSurrounds* is true, all atoms not found in (1) or (2) will be hidden.
   Cartoon display will not be affected.
6. If *focus* is true, the view will be reset to centre on and encompass the
   covered region.

.. _`symcopies`:

clipper symcopies
-----------------

Syntax: clipper symcopies [*structures*]
[**contacting** *atoms*] [**distance** *number* (5.0)]
[**name** *string*] [**focus** *true/false* (false)]
[**reflFile** *filename*] [**pruneSpecialPositions** *true/false* (true)]

Make real, whole-model copies of crystallographic symmetry mates and merge them
- together with the original ASU - into a single new atomic model, using
ChimeraX's ``combine`` mechanism. Unlike the live symmetry display (which
shows only the atoms near the spotlight), every copy produced here is a complete
model, and the result is an ordinary structure you can edit, save or otherwise
manipulate independently of the crystal. This command runs both in the GUI and
headless (``--nogui``): it does not require the live symmetry display, so it can
be scripted for batch symmetry expansion.

Which symmetry copies are realised depends on whether **contacting** is given:

* **contacting** *given* - every symmetry copy with at least one atom approaching
  within *distance* Å of the *contacting* atom selection is realised, for each
  crystal structure the selection touches. (If *structures* is also given, only
  those structures are considered.)
* **contacting** *omitted* - if a live spotlight is displaying symmetry copies
  (GUI only), exactly those currently-drawn copies are realised. With no spotlight
  - for example headless - every symmetry copy within *distance* Å of the whole
  model is realised instead. This applies to each of *structures*, or - if
  *structures* is also omitted - to every Clipper-managed model in the session.

**reflFile** takes a structure-factor file (``.mtz`` or ``.cif``) and uses its
cell and spacegroup for models whose own metadata (PDB ``CRYST1`` / mmCIF header)
does not carry them; no maps are generated. Models with no usable crystallographic
symmetry are skipped with a warning.

By default (**pruneSpecialPositions true**) a fragment copy that maps back exactly
onto a molecule already present - i.e. one sitting on a special position - is
dropped as a redundant duplicate. Which atoms lie on a special position is
determined exactly from Clipper's crystallographic site multiplicity (not a
distance/size heuristic), so a metal on a crystallographic axis in a protein is
handled just as cleanly as an ion or water in a small-molecule cell. Dropping
exactly the redundant copies is also what makes the occupancy correct after an
inverse collapse (below). Set **pruneSpecialPositions false** to realise every copy
verbatim. Molecules split across a special position should first be completed with
:ref:`fragments` ``complete`` (also needed for sensible connectivity); the command
warns if it detects an un-completed one.

**Invertible expansion.** The new model carries a ``clipper_sym_expansion`` record
mapping each symmetry operator to the chain IDs it produced (with its inverse and
``n_asu``). Every copy is given unique chain IDs, so nothing is renamed on merge and
the mapping is exact. From the Python layer this lets a caller fold the expanded
model back onto the original ASU - ``from chimerax.clipper import collapse_to_asu``
applies each operator's inverse and sets occupancy ``1/n_asu`` - which is how a full
P1 box can be relaxed (e.g. by an MD engine) and then collapsed to a single ASU for a
structure-factor calculation.

The new model is placed exactly where the crystallographic copies sit relative
to the original, so - unlike the general-purpose core symmetry tools - the view
is left undisturbed. Pass **focus true** to reset the view to encompass the
result (GUI only). Use **name** to set the new model's name (the default is
"*<structure>* symmetry copies").

.. _`unitcells`:

clipper unitcells
-----------------

Syntax: clipper unitcells [*structures*]
[**cells** *n,m,o*] [**box** *x,y,z*]
[**name** *string*] [**focus** *true/false* (false)]
[**reflFile** *filename*] [**pruneSpecialPositions** *true/false* (true)]

Generate one or more complete crystallographic unit cells as a single new,
real atomic model. Each unit cell is the model's ASU replicated by every
space-group symmetry operator; a block of cells adds the corresponding lattice
translations. Like :ref:`symcopies`, this runs both in the GUI and headless
(``--nogui``).

The size of the block is set by exactly one of:

* **cells** *n,m,o* - an *n* × *m* × *o* block of whole unit cells along the
  **a**, **b** and **c** cell axes (e.g. ``cells 2,2,2`` for eight cells).
* **box** *x,y,z* - the smallest block of whole unit cells that completely fills
  an *x* × *y* × *z* Å region. Each axis count is rounded up
  (*n* = ⌈*x*/*a*⌉, etc.), so the result comes out slightly larger than the
  requested box. (For non-orthogonal cells the box is filled approximately, per
  cell axis.)
* *neither* - a single unit cell.

**reflFile**, **name**, **focus** and **pruneSpecialPositions** behave exactly as
for :ref:`symcopies`; in particular, special-position copies are pruned by default
using exact site multiplicity, and the result carries the same invertible
``clipper_sym_expansion`` record (so a packed box can be collapsed back to the ASU
with ``collapse_to_asu``). Models with no usable crystallographic symmetry are
skipped with a warning (a single unit cell of a ``P 1`` structure is just the input
model, so nothing is generated). For an MD-ready box, complete any
special-position-split molecules first with :ref:`fragments` ``complete``.

.. _`bsharp`:

clipper bsharp
--------------

Syntax: clipper bsharp *b-factor* [*maps*]

Apply a sharpening (positive *b-factor*) or smoothing (negative *b-factor*)
B-factor, in Å\ :sup:`2`, to live crystallographic maps. The affected maps are
recomputed immediately from their cached structure-factor coefficients - no
atoms are moved and no F\ :sub:`calc` recalculation is performed, so the change
is near-instantaneous. Contour levels are held fixed in units of sigma, so the
surface sharpens or smooths in place rather than jumping.

If *maps* is omitted, the change is applied to the "designated" viewing map of
every crystallographic dataset in the session - i.e. the heuristically
sharpened 2mFo-DFc map that :ref:`open` displays as a surface, not the raw
reference map or the difference map. The same designated-only behaviour applies
when *maps* names a whole model or dataset, or more than one map. To retune one
specific map, name that single map explicitly.

Difference maps, and maps that have been explicitly locked (for example the map
ISOLDE fits against during an interactive simulation), are never changed;
naming such a map explicitly reports an error rather than silently doing
nothing.

*(This command applies only to macromolecular live maps. Small-molecule maps
created with* :ref:`smallmol` *use a different scaling scheme and do not support
a runtime sharpening B-factor.)*

.. _`set oversampling`:

clipper set oversampling
------------------------

Syntax: clipper set oversampling *rate* [*models*]

Change the real-space oversampling (Shannon) *rate* of crystallographic maps at
runtime - the same parameter set at load time by the **overSampling** option of
:ref:`open`. The map grid is rebuilt and every live and static map is
re-transformed onto it (a re-FFT from the already-cached coefficients: the
experimental reflection data, R-free set and calculated structure factors are
all untouched), and the crystallographic symmetry search is rebuilt to match the
new grid.

Higher rates give a finer real-space grid and smoother-looking contours, at the
cost of memory and FFT time; typical values are 1.5-3.0 (the default at load is
2.0). If *models* is omitted, the new rate is applied to every crystallographic
map manager in the session; otherwise only to the managers belonging to the
specified models.

.. note::

   While an interactive ISOLDE simulation is running against a live map, the
   oversampling rate is temporarily locked and this command reports an error:
   the simulation holds a fixed-size density box on the GPU that re-gridding
   would corrupt. Try again once the simulation has finished.
