============
Installation
============

-------------------------
Installation From Release
-------------------------

pyNastran is an easy package to install once you have the required Python modules.  It's a pure Python package so you shouldn't have too many problems.  On the command line:

``pip install pyNastran -U``

That will install the minimum set of what you need to run pyNastran (so no GUI).  If you want the GUI, choose one Qt binding.  On the command line, type:

    * ``pip install pyNastran[gui_pyside6] -U``
    * ``pip install pyNastran[gui_pyqt6]   -U``
    * ``pip install pyNastran[gui_pyqt5]   -U``

Python
------
The current development version supports Python 3.11 - 3.15.

Packages
--------
pyNastran is tested against a range of package versions (lowest to highest based on availbility), so it should work.  The recommended set of packages are:

 * **Required**:

   * numpy > 1.25, < 3
   * scipy > 1.8, < 2
   * matplotlib > 3.6, < 4
   * cpylog >= 1.6.0
   * docopt-ng >= 0.9.0   **(required for command line tools)**
   * colorama             **(colored logging)**
   * qtpy

 * **Optional**:

   * pandas >= 1.4, < 3
   * h5py                 **(HDF5 support; installed on Python < 3.15)**
   * tables (PyTables)    **(included in the dev extra on Python < 3.15)**

 * **GUI**:

   * vtk >= 9.0, < 10
   * Qt **(pick one)**

      * PyQt5
      * PyQt6
      * PySide6

   * QScintilla >= 2.13.0 **(optional for fancy scripting; PyQt5/6 only)**
   * pygments >= 2.2.0 **(optional for fancy scripting; PyQt5/6 only)**
   * imageio < 3          **(animation support)**
   * pillow >= 10.1.0
   * natsort
   * python-docx
   * openpyxl

*****************************************************
Install Procedure - From Source
*****************************************************
Base functionality:

 * `64-bit Python <https://www.python.org/downloads/>`_

On the command line, navigate to the directory with ``models``.  Then type:

 * ``pip install .``

If you would like to be able to modify the checked out source and have it have an effect, add ``-e`` before the ``.``:

 * ``pip install -e .``

If you would like pandas and hdf5 support to the BDF/OP2, but no gui options.

 * On the command line:

    * ``pip install .[nogui]``

For **optional** GUI support, choose one Qt binding:

 * On the command line:

    * ``pip install .[gui_pyside6]``
    * ``pip install .[gui_pyqt6]``
    * ``pip install .[gui_pyqt5]``

If you want to just install the dependencies yourself:

 * On the command line:

    * ``pip install . --no-dependencies``

Documentation
=============
If you don't want to use build the docs, just use the docs on the web.

See `docs <https://pynastran-git.readthedocs.io/en/latest/>`_

If you want to build the docs yourslf:
 - https://pynastran-git.readthedocs.io/en/latest/installation/building_docs.html
