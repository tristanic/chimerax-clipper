ChimeraX-Clipper Python API
===========================

.. toctree::
   :maxdepth: 2

.. contents::
    :local:

Note: this documentation is for the ChimeraX-specific portion of the ChimeraX-Clipper API.
For most day-to-day work in ChimeraX you should only need to use the :ref:`clipper_commands`.
The vast majority of the lower-level API is a direct mapping of the underlying `Clipper C++ API`_ 
(with some modifications to support threading of performance-critical functions),
and your first stop should be that documentation. When in doubt, the PyBind11 source code for 
the bindings themselves may be found `here`_.

.. _Clipper C++ API: http://www.ysbl.york.ac.uk/~cowtan/clipper/doc/annotated.html
.. _here: https://github.com/tristanic/chimerax-clipper/tree/master/src/bindings

High-level Management
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: chimerax.clipper.symmetry

    .. autoclass:: SymmetryManager
        :members:

    .. autoclass::AtomicSymmetryModel
        :members:

    .. autoclass:: Unit_Cell
        :members:

Map Management
~~~~~~~~~~~~~~

.. automodule:: chimerax.clipper.maps

    .. autoclass:: MapMgr
        :members:
    
    .. autoclass:: XmapSet
        :members:
    
    .. autoclass:: XmapHandler_Static
        :members:
    
    .. autoclass:: XmapHandler_Live
        :members:
    
    .. autoclass:: NXmapSet
        :members:
    
    .. autoclass:: NXmapHandler
        :members:

File I/O
~~~~~~~~

.. py:function:: chimerax.clipper.clipper_mtz.load_hkl_data


Mouse Modes
~~~~~~~~~~~

.. automodule:: chimerax.clipper.mousemodes

    .. autoclass:: RotateMouseMode
        :members:
    
    .. autoclass:: ClipPlaneAdjuster
        :members:
    
    .. autoclass:: Z_Shift_CofR
        :members:
    
    .. autoclass:: ZoomMouseMode
        :members:
    
    .. autoclass:: SelectVolumeToContour
        :members:
    
    .. autoclass:: ContourSelectedVolume
        :members:
    
    .. autoclass:: ShiftToReferenceAsuMenuEntry
        :members:

Miscellaneous Utilities
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: chimerax.clipper.util
    :members:

