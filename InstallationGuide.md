# Introduction #
  * pyNastran is a very easy package to install once you have the required Python modules.  It's a pure Python package so no one should really be having too many problems.
  * The OP4 reader is a bit harder to get working, but is an optional package that requires building of C code in addition to requiring Cython.  Assuming you can build one of the sample Cython examples, it'll work.
  * It is recommended that all users install pyNastran as developers.
# User Overview #
  * Install Python 2.6, 2.7, 3.1 or 3.2
  * Install numpy, scipy, vtk, wxPython (Python 2.x)
  * Install numpy, scipy, distribute (Python 3.x)
  * Install Cython for OP4 support
  * Download the latest pyNastran release
  * Extract the zip file
  * Navigate to pyNastran directory on the command line
  * run **python setup.py develop** to install the package
  * run **python setup.py test** to run the tests.  All should pass
  * run **make** in the **pyNastran/op4** directory

# Install Procedure - User #

## Download Python ##
The software has been tested with Python 2.7 and Python 3.2.

I use **Python(x,y)-2.7.2.1**, but it shouldn't matter too much as long as your versions of numpy/scipy are consistent.  Python 2.6 and 2.7 are both fine, but the code uses a few features that aren't available in Python 2.5.

  * [Python(x,y) ](http://code.google.com/p/pythonxy/wiki/Welcome)
  * [scipy 0.9.0  ](http://scipy.org/)
  * [numpy 1.6.1  ](http://numpy.org/)
  * [vtk 5.8.0](http://www.vtk.org/VTK/resources/software.html) (for the GUI)
  * [wxPython 2.8.12.0](http://wxpython.org/download.php#stable) (for the GUI)
  * [distribute\_setup.py](http://pypi.python.org/pypi/distribute) (for Python 3.x)
  * [matplotlib 1.1.0](http://matplotlib.sourceforge.net/) (version 0.5+, for plotting)
  * [Cython 0.15.1](http://cython.org/) (for OP4 support)

## Install Python C compiler & Cython ##
  * Microsoft Visual Studio 2008 ([follow this](http://wiki.tiker.net/PyCuda/Installation/Windows)).  The version of Visual Studio is important as Python is built against MSVC 2008.  Don't worry, it's free.
  * Linux/Mac: Install gcc
  * Open a command prompt / terminal
  * run **easy\_install cython==0.15.1** to install Cython

## Download Distribute\_Setup.py (Python 3.x ONLY) ##
  * Download [distribute\_setup.py](http://pypi.python.org/pypi/distribute) to replace setuptools.  pyNastran needs it to run setup.py.

## Download pyNastran ##
  * Click the downloads page and download the most recent [version](http://code.google.com/p/pynastran/downloads/list)

## Install pyNastran ##
  * Navigate to pyNastran directory on the command line
  * run **python setup.py develop** to install pyNastran.  Don't type **python setup.py install**
  * run **python setup.py test** to run the tests.  Make sure all of them pass.

## Install OP4 Reader ##
  * Navigate to pyNastran/op4 directory on the command line

### Windows ###
  * run **make** to install the OP4 code
  * if that doesn't work, try **python setup.py build\_ext --inplace**

Enjoy!

# Developer Overview #
  * Install Python
  * Install numpy, scipy, matplotlib, vtk, PyQt4 / wxPython
  * install Sphinx, GraphViz (for documentation
  * Install Cython (if you want OP4 support)
  * Use subversion to checkout pyNastran
  * Navigate to pyNastran directory on the command line
  * run **python setup.py develop** to install the package

  * Note: The GUI requires PyQt4, but the BDF, OP2, F06, and OP4 have support for either wxPython or PyQt4 to select a file explicitly.
  * [PyQt4](http://www.riverbankcomputing.com/software/pyqt/download)

# Install Procedure - Developer #

## Install Subversion ##
  * download install [Subversion](http://tortoisesvn.net/)
  * checkout the code (see the subversion website if you need help) and the [Source](http://code.google.com/p/pynastran/source/checkout)


## Create the Documentation ##
  * Navigate to pyNastran directory on the command line
  * run sphinx to generate the documentation