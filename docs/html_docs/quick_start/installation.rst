============
Installation
============

-------------------------
Installation From Release
-------------------------

pyNastran is an easy package to install once you have the required Python
modules.  It's a pure Python package so you shouldn't have too many problems.
Just type:

``pip install pyNastran``


Python
------
The software is tested against:
 - Python 2.7.14 (Windows/Linux)
 - Python 3.5 (Windows/Linux)
 - Python 3.6 (Linux)

Packages
--------
The suggested set of packages include:
 - numpy 1.14
 - scipy 1.0
 - docopt == 0.6.2
 - VTK 7/8
 - PyQt 4/5
 - pandas ???
 - matplotlib >= 2.1.2
 - imageio >= 2.2.0
 - typing >= 3.6.1
 - pathlib2 >= 2.2.0
 - scandir >= 1.4.0

***********************************************
Install Procedure - From Anaconda (recommended)
***********************************************
Base functionality:
 - `Anaconda Python <https://store.continuum.io/cshop/anaconda/>`_
 * ``conda install numpy``
 * ``conda install scipy``
 * ``conda install docopt``
 * ``conda install typing`` (Python 2.7)
 * ``conda install pathlib2`` (Python 2.7)
 * ``conda install scandir`` (Python 2.7)
 * ``conda install pandas`` (optional)
 * ``conda install h5py`` (optional)
 * ``conda install matplotlib`` (optional)
 * ``pip install pyNastran``

For gui support (optional; required for GUI):
 * Python 2.7: 
   * From `Windows binaries <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_:
     * Download VTK:
       - VTK-7.1.1-cp27-cp27m-win_amd64.whl (Python 2.7)
  On the command line:
    * ``conda install imageio``
    * ``conda install pyqt``
    * ``pip install VTK*.whl``
    * ``pip install qtpy``

 * Python 3.5, 3.6, 3.7:
   * On the command line:
    * ``conda install imageio``
    * ``conda install pyqt``
    * ``conda install vtk``
    * ``pip install qtpy``


************************************
Install Procedure - From Base Python
************************************
Base functionality:
 * `Windows <https://www.python.org/downloads/windows/>`_
   * Download Windows x86-64 MSI installer
 * Linux/Mac `Python <https://www.python.org/downloads/>`_
   - Make sure to get 64-bit Python.
 * On the command line:
 
   * ``pip install numpy``
   * ``pip install scipy``
   * ``pip install docopt``
   * ``conda install typing`` (Python 2.7)
   * ``conda install pathlib2`` (Python 2.7)
   * ``conda install scandir`` (Python 2.7)
   * ``pip install pandas`` (optional)
   * ``pip install h5py`` (optional)
   * ``pip install matplotlib`` (optional)
   * ``pip install pyNastran``


For gui support (optional; required for GUI):
 * From `Windows binaries <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_:
   * Download VTK:
     - VTK-7.1.1-cp27-cp27m-win_amd64.whl (Python 2.7)
     - VTK-8.1.1-cp35-cp35m-win_amd64.whl (Python 3.5)
     - VTK-8.1.1-cp36-cp36m-win_amd64.whl (Python 3.6)
     - VTK-8.1.1-cp37-cp37m-win_amd64.whl (Python 3.7)
    * Download PyQt:
     - PyQt4-4.11.4-cp27-cp27m-win_amd64.whl (Python 2.7)
     - PyQt4-4.11.4-cp35-cp35m-win_amd64.whl (Python 3.5)
     - PyQt4-4.11.4-cp36-cp36m-win_amd64.whl (Python 3.6)
  On the command line:
    * ``pip install imageio``
    * ``pip install VTK*.whl``
    * ``pip install PyQt4*.whl``
    * ``pip install qtpy``

***********************************************
Install Procedure - From WinPython (incomplete)
***********************************************
 * `WinPython <http://winpython.sourceforge.net/>`_


Use Web docs
------------
See `docs <http://http://pynastran.m4-engineering.com/>`_

-------------------------
Installation From Source
-------------------------

pyNastran is an easy package to install once you have the required Python
modules.  It's a pure Python package so you shouldn't have too many problems.

Installing from source is recommened if:
 - You want the most recent version (see installation.rst-master)
 - You want easier access to the source
 - You're on an air-gapped machine

Overview
========
 * Install Python (see :doc:`installation_release`)
   - skip the `pip install pyNastran` step
 * Install Sphinx, GraphViz, alabaster (for documentation)

 * Install Git
 * Clone pyNastran-master from Github
 * Install pyNastran

Install extra packages (for Doucmentation)
==========================================

Install `GraphViz  <https://www.graphviz.org/>`_

Install additional python packages

.. code-block:: console

  pip install Sphinx
  pip install alabaster
  pip install numpydoc

Install Git
===========

 * Download & install `Git <http://git-scm.com/>`_ (required)
 * Download a GUI for Git (optional)
    * `TortoiseGit <https://code.google.com/p/tortoisegit/>`_ (recommended for Windows)


Install pyNastran
=================
There are two ways to install the master (dev) version of pyNastran

 1. Download the most recent `zip version <https://github.com/SteveDoyle2/pynastran/archive/master.zip>`_

 2. Clone pyNastran (see below).  Using Git allows you to easily update to the
    latest dev version when you want to as well as push any commits of your own.

If you don't want the gui, use ``setup_no_gui.py`` instead of ``setup.py``.

.. code-block:: console

  >>> python setup.py install

or:

.. code-block:: console

  >>> python setup_no_gui.py install


Cloning pyNastran using TortoiseGit
===================================
Right-click in a folder and select ``Git Clone``.

.. image:: clone.png

Enter the above information.  If desired, click the branch box and and enter a branch name
and click ``OK``.

Cloning pyNastran Using Command Line
====================================
Checkout/clone the dev code by typing (preferred):

.. code-block:: console

  >>> git clone https://github.com/SteveDoyle2/pynastran


To checkout a branch

.. code-block:: console

  >>> git.exe clone --branch v1.0-dev --progress -v "https://github.com/SteveDoyle2/pyNastran.git" "C:\\work\\pyNastran_v1.0-dev"


Documentation
=============
Two options for documentation exist.

Build Docs
----------
Navigate to ``pyNastran/docs_sphinx`` directory on the command line.

.. code-block:: console

  >>> make html

Alternatively, see `docs <http://pynastran-git.readthedocs.org/en/latest/>`_

