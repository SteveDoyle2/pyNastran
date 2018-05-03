=========================
Installation From Release
=========================

pyNastran is an easy package to install once you have the required Python
modules.  It's a pure Python package so you shouldn't have too many problems.

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
 * `conda install numpy`
 * `conda install scipy`
 * `conda install docopt`
 * `conda install typing` (Python 2.7)
 * `conda install pathlib2` (Python 2.7)
 * `conda install scandir` (Python 2.7)
 * `conda install pandas` (optional)
 * `conda install h5py` (optional)
 * `conda install matplotlib` (optional)
 * `pip install pyNastran`

For gui support (optional; required for GUI):
 * Python 2.7: 
   * From `Windows binaries <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_:
     * Download VTK:
       - VTK-7.1.1-cp27-cp27m-win_amd64.whl (Python 2.7)
  On the command line:
    * `conda install imageio`
    * `conda install pyqt`
    * `pip install VTK*.whl`
    * `pip install qtpy`

 * Python 3.5, 3.6, 3.7:
   * On the command line:
    * `conda install imageio`
    * `conda install pyqt`
    * `conda install vtk`
    * `pip install qtpy`

************************************
Install Procedure - From Base Python
************************************
Base functionality:
 * https://www.python.org/downloads/windows/
   * Download Windows x86-64 MSI installer
 * Download 64-bit `Python <https://www.python.org/downloads/>`_
   - Make sure to get 64-bit Python.
 * On the command line:
   * `pip install numpy`
   * `pip install scipy`
   * `pip install docopt`
   * `conda install typing` (Python 2.7)
   * `conda install pathlib2` (Python 2.7)
   * `conda install scandir` (Python 2.7)
   * `pip install pandas` (optional)
   * `pip install h5py` (optional)
   * `pip install matplotlib` (optional)
   * `pip install pyNastran`


For gui support (optional; required for GUI):
 * From `Windows binaries <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_:
   * Download VTK:
     - VTK-7.1.1-cp27-cp27m-win_amd64.whl (Python 2.7)
     - VTK-8.1.0-cp35-cp35m-win_amd64.whl (Python 3.5)
     - VTK-8.1.0-cp36-cp36m-win_amd64.whl (Python 3.6)
     - VTK-8.1.0-cp37-cp37m-win_amd64.whl (Python 3.7)
    * Download PyQt:
     - PyQt4-4.11.4-cp27-cp27m-win_amd64.whl (Python 2.7)
     - PyQt4-4.11.4-cp35-cp35m-win_amd64.whl (Python 3.5)
     - PyQt4-4.11.4-cp36-cp36m-win_amd64.whl (Python 3.6)
  On the command line:
    * `pip install imageio`
    * `pip install VTK*.whl`
    * `pip install PyQt4*.whl`
    * `pip install qtpy`

***********************************************
Install Procedure - From WinPython (incomplete)
***********************************************
 * `WinPython <http://winpython.sourceforge.net/>`_


Use Web docs
------------
See <http://pynastran-git.readthedocs.org/en/latest/>`_

