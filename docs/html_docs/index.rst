.. pyNastran documentation master file, created by
   sphinx-quickstart2 on Tue Aug 14 13:26:34 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyNastran's documentation for v0.8!
==============================================
The pyNastran software interfaces to Nastran's complicated input and output
files and provides a simplified interface to read/edit/write the various files.
The software is compatible currently being used on Windows, Linux, and Mac.

The **BDF reader/editor/writer** supports 278 cards including coordinate
systems. Card objects have methods to access data such as Mass, Area, etc.
The BDF writer writes a small field formatted file, but makes full use of
the 8-character Nastran field. The OpenMDAO BDF parametrization syntax
is also supported.

The **OP2 reader** supports static/transient results, which unless you
analyzing frequency response data should be good enough. It also supports
**F06 Writing** for most of the objects. Results include: displacement,
velocity, acceleration, temperature, eigenvectors, eigenvalues, SPC forces,
MPC forces, grid point forces, load vectors, applied loads, strain energy,
as well as stress and strain.

The **Python OP4** reader/writer supports reading ASCII/binary sparse and dense
matrices, and writing ASCII matrices..

A simple GUI has been developed that can view BDF models and display static/dynamic
stress/strain/displacement/eignevectors (they must be real!) results from the OP2.
Additionally, Cart3d, STL, and Panair are included.


.. toctree::
   :maxdepth: 4

   quick_start/index
   index/pyNastran
..   manual/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

