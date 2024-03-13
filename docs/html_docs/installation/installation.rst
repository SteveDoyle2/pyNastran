============
Installation
============

-------------------------
Installation From Release
-------------------------

pyNastran is an easy package to install once you have the required Python modules.  It's a pure Python package so you shouldn't have too many problems.  On the command line:

``pip install pyNastran -U``

That will install the minimum set of what you need to run pyNastran (so no GUI).  If you want the GUI, there are two options.  On the command line, type:

    * ``pip install pyNastran[gui_pyside2] -U``
    * ``pip install pyNastran[gui_pyqt5]   -U``

Python
------
The software is tested on Windows and Linux against Python 3.9 - 3.12.

The latest dev version supports Python 3.9 - 3.12.

Packages
--------
pyNastran is tested against a range of package versions (lowest to highest based on availbility), so it should work.  The recommended set of packages are:

 * **Required**:

   * numpy >= 1.14
   * scipy >= 1.0
   * cpylog >= 1.4.0
   * docopt-ng == 0.7.2   **(required for command line tools)**

 * **Optional**:

   * colorama >= 0.3.9    **(colored logging)**
   * pandas >= 0.25,<2.0
   * matplotlib >= 2.2.4  **(plotting)**
   * h5py >= 2.8.0        **(HDF5 support)**

 * **GUI**:

   * vtk >= 9.0, <9.4
   * qtpy >= 1.4.0
   * Qt **(pick one)**

     * PyQt5 >= 5.9.2
     * PySide2 >= 5.11.2
     * PyQt6   (currently disabled)
     * PySide6 (currently disabled)
   * QScintilla >= 2.13.0 **(optional for fancy scripting; PyQt5/6 only)**
   * pygments >= 2.2.0 **(optional for fancy scripting; PyQt5/6 only)**
   * imageio >= 2.4.1,<3  **(optional for animation support)**

*****************************************************
Install Procedure - From Sorce
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

For **optional** GUI support, there are two options:

 * On the command line:

    * ``pip install .[gui_pyside2]``
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

