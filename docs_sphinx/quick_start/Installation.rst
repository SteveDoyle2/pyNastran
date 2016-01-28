============
Installation
============
pyNastran is an easy package to install once you have the required Python 
modules.  It's a pure Python package so you shouldn't have too many problems.

Install Procedure - From a Release
==================================

Overview
========
 * Install Python
 * Install numpy, scipy, pandas, vtk, PyQt4 (Python 2.x)
 * Install numpy, scipy, pandas (Python 3.x)
 * Download pyNastran from Github
 * Navigate to pyNastran directory on the command line
 * Install the package


Download Python
---------------
v0.8 is tested against:
 - Python 2.7.8
 - Python 3.3
 - Python 3.4
 - Python 3.5


Options include:
 * `Anaconda Python <https://store.continuum.io/cshop/anaconda/>`_ (recommended)
 * `WinPython <http://winpython.sourceforge.net/>`_
 * Base `Python <https://www.python.org/downloads/>`_ with the unofficial 
   `Windows binaries <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_

Make sure to get 64-bit Python, so memory usage becomes a non-issue.  It
shouldn't matter too much as long as your packages versions (e.g. numpy/scipy)
are consistent.  **With either distribution, both Python and all 3rd party
packages will be installed.  With base Python, you need to do that yourself.**


Additional packages
-------------------
The following packages are required.

 * `scipy <http://scipy.org/>`_
 * `numpy 1.9+ <http://numpy.org/>`_
 * `vtk 5.x or 6.x <http://www.vtk.org/VTK/resources/software.html>`_ (for the GUI)
 * `wxPython <http://wxpython.org/download.php#stable>`_ (for BDF/OP2/OP4 popups)
 * `PyQt4 <http://www.riverbankcomputing.com/software/pyqt/download>`_ (for the GUI; BDF/OP2/OP4 popups in v0.7)
 * `pandas <http://pandas.pydata.org/>`_ (optional way to use the OP2 in iPython; v0.8)

PyNastran's package requirements are tested with packages no older than 1 year
at the time of release.  If you require an older version, try changing version
requirements in `setup.py`.  It shouldn't be that different, but pyNastran does
make use of numpy's "new" axis option in `numpy.linalg.norm` and there was a
major bug fix in Python 2.7.8, so buyer beware.

Download pyNastran
------------------

 * Click the downloads page and download the most recent `zip version
 <https://github.com/SteveDoyle2/pynastran/archive/master.zip>`_ or `clone
 <github-windows://openRepo/https://github.com/SteveDoyle2/pynastran>`_ 
 it using Git.  Using Git allows you to easily update to the latest dev version
 when you want to as well as push any commits of your own.


Install pyNastran
-----------------
 * Navigate to pyNastran directory on the command line.  The ``setup.py`` file
   should exist in the current directory.
 
 * Either run...
   1. Able to edit the source code and have the changes propogate

    ``
    >>> python setup.py develop
    ``

   2. Changes will not propogate

    ``
    >>> python setup.py install
    ``

If you don't want the gui, use ``setup_no_gui.py`` instead of ``setup.py``.

Install Procedure - From Source
===============================

Overview
========
 * Install Python
 * Install numpy, scipy, pandas, vtk, PyQt4 (Python 2.x)
 * Install numpy, scipy, pandas (Python 3.x)
 * Install Sphinx, GraphViz, alabaster (for documentation)

 * Install Git
 * Use Git to clone pyNastran and checkout the `master`
 * Clone the latest pyNastran from Github
 * Install pyNastran

Install extra Python packages
-----------------------------
Install Sphinx and alabaster

``
pip install Sphinx
pip install alabaster
``

Install Git
-----------

 * Download & install `Git <http://git-scm.com/downloads/guis/>`_
 * Optionally, download a GUI for Git
    * `TortoiseGit <https://code.google.com/p/tortoisegit/>`_ (recommended for Windows)


Cloning pyNastran
-----------------
To checkout a branch

``
>>> git.exe clone --branch v0.8 --progress -v "https://github.com/SteveDoyle2/pyNastran.git" "C:\\work\\pyNastran_v0.8"
``

Checkout/clone the code by typing

``
>>> git clone https://github.com/SteveDoyle2/pynastran
``

or using a GUI and cloning the project.

If this seems hard, try TortoiseGit.  It's much easier.

Install pyNastran
-----------------
see above

Documentation
-------------
Two options for documentation exist

Build Docs
^^^^^^^^^^
Navigate to `pyNastran/docs_sphinx` directory on the command line.

``
>>> make html
``

Use existing docs
^^^^^^^^^^^^^^^^^
Use the `web docs <http://pynastran-git.readthedocs.org/en/latest/>`_

