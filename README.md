
|  Version  | Docs  | Status |
| :--- 	  | :--- 	  | :--- 	  |
|  [![PyPi Version](https://img.shields.io/pypi/v/pynastran.svg)](https://pypi.python.org/pypi/pyNastran) | [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=v1.0)](http://pynastran-git.readthedocs.io/en/v1.0.0/?badge=v1.0.0) | [![Build Status](https://img.shields.io/travis/SteveDoyle2/pyNastran/v1.0.svg)](https://travis-ci.org/SteveDoyle2/pyNastran) [![Coverage Status](https://img.shields.io/coveralls/SteveDoyle2/pyNastran/v1.0.svg)](https://coveralls.io/github/SteveDoyle2/pyNastran?branch=v1.0) |
|   Master | [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=latest)](http://pynastran-git.readthedocs.io/en/latest/?badge=latest) | [![Linux Status](https://img.shields.io/travis/SteveDoyle2/pyNastran/master.svg)](https://travis-ci.org/SteveDoyle2/pyNastran) [![Windows Status](https://ci.appveyor.com/api/projects/status/1qau107h43mbgghi/branch/master?svg=true)](https://ci.appveyor.com/project/SteveDoyle2/pynastran) [![codecov](https://codecov.io/gh/SteveDoyle2/pyNastran/branch/master/graph/badge.svg)](https://codecov.io/gh/SteveDoyle2/pyNastran) | 

<!---
[![Coverage Status](https://img.shields.io/coveralls/SteveDoyle2/pyNastran/master.svg)](https://coveralls.io/github/SteveDoyle2/pyNastran?branch=master)
--->

<!---

[![Build status](https://ci.appveyor.com/api/projects/status/1qau107h43mbgghi/branch/master?svg=true)](https://ci.appveyor.com/project/SteveDoyle2/pynastran/branch/master)

[![Build status](https://ci.appveyor.com/api/projects/status/1qau107h43mbgghi?svg=true)](https://ci.appveyor.com/project/SteveDoyle2/pynastran)

[![Build Status](https://travis-ci.org/SteveDoyle2/pyNastran.png)](https://travis-ci.org/SteveDoyle2/pyNastran)
[![Coverage Status](https://coveralls.io/repos/github/SteveDoyle2/pyNastran/badge.svg?branch=master)](https://coveralls.io/github/SteveDoyle2/pyNastran?branch=master)
[![Requirements Status](https://img.shields.io/requires/github/SteveDoyle2/pyNastran/master.svg)](https://requires.io/github/SteveDoyle2/pyNastran/requirements/?branch=master)	|
--->

### v1.0.0

[Download GUI](https://sourceforge.net/projects/pynastran/files/?source=navbar) (latest is from 2017/5/25)

<!---
[Doumentation](http://pynastran-git.readthedocs.org/en/latest/index.html)
--->

Also, check out the:
  * [Discussion forum](http://groups.google.com/group/pynastran-discuss)
  * [Developer forum](http://groups.google.com/group/pynastran-dev)

for more detailed information.

<!--- this isn't setup... -->
<!--- http://stevedoyle2.github.io/pyNastran/ --->

# Overview

pyNastran is an interface library to the various Nastran file formats (BDF, OP2, OP4).  Using the BDF interface, you can read/edit/write Nastran geometry without worrying about field formatting.  Many checks are also performed to verify that your model is correct.  Using the OP2 interface, you can read very large result files very quckly and very efficiently.  Additionally, you can also extract a subset of the result data and write F06 result files.

Using the pyNastran GUI, you can read in Nastran models and quickly view results for a model.  While it's no FEMAP/Patran, it can replace many tasks that would otherwise require a commercial program.

![GUI](https://github.com/SteveDoyle2/pynastran/blob/master/pyNastran/gui/images/caero.png)


<!--- Ripped off meshio  --->
<!--- [![Build Status](https://travis-ci.org/SteveDoyle2/pyNastran.svg?branch=master)](https://travis-ci.org/SteveDoyle2/pyNastran)  --->
<!--- [![codecov.io](https://codecov.io/github/SteveDoyle2/pyNastran/coverage.svg?branch=master)](https://codecov.io/github/SteveDoyle2/pyNastran?branch=master)  --->

<!--- ## pyNastran v0.8.0 has NOT been released (8/21/2016)   --->
<!--- [Download pyNastran v0.8] (https://github.com/SteveDoyle2/pyNastran/releases)  --->

<!--- 

This should be hidden...

v1.1-progress

Sort of the same goals as v1.0.  Just trying to make it not buggy.  There are some new features though.

Programmatics
 - dropping VTK 5/6 suppoprt for the GUI

OP2:
 - HDF5 export support
 - pandas support for matrices
 - couple more results vectorized (e.g., complex strain energy, DMIG strain energy, some forces)

BDF:
 - currently 330 cards supported
 - cross-referencing is now simpler
   - ``*_ref`` attributes are cross-referenced
      - ``element.nodes`` is not cross-referenced
      - ``element.nodes_ref`` is cross-referenced
 - pickling to load your deck ~5x faster
 - decreased time required for Case Control Deck with large SETs and many load cases
- improved optimization checks

GUI:
 - animation menu is more intuitive
    - in gui animation
    - multiple animation profiles
    - wipe deformed shape button (currently broken)
    - link to legend change (not done)
    - no longer a subset of the legend menu (not done)
 - nominal geometry (useful for deflection plots)
 - single click menus

  --->
  
# News

### pyNastran v1.0.0 has been released (5/25/2017)
This is a major release.  The focus this time has been on robustness and testing.
Hopefully, it shows.  The software has also been relicensed to be **BSD-3**, which
is a more permissive license and is the same one that numpy, scipy, and
matplotlib use.

Unfortunately, the GUI is more complicated.
 - For open source projects : LGPL 2/3
 - For companies that pay a license to Riverbank : LGPL 2/3
 - For companies that don't pay a license fee : GPL 2/3

However, you may distribute an unmodified binary.

 - Programmatics:
   - Dropping Python 3.3 support
   - Adding Python 3.6 support

 - [OP2](http://pynastran-git.readthedocs.io/en/latest/quick_start/op2_demo.html)
   - preliminary random results (ATO/CRM/PSD/RMS/NO)
   - improved geometry support
   - MONPNT1/MONPNT3 reading
   - MATPOOL matrices

 - F06 scripts
   - added a preliminary flutter (SOL 145) parser
     - supports multiple subcases
     - PK and PKNL methods supported
     - `plot_Vg_Vf(...)`, `plot_Vg(...)`, `plot_root_locus(...)`
     - input/output units
     - mode switching not fixed yet

 - [GUI](http://pynastran-git.readthedocs.io/en/latest/quick_start/gui.html):
   - complex displacement support
   - animation support
   - vector results (real/complex)
      - SPC Forces, MPC Forces, Load Vector, Applied Load
        - minimal control presently
      - No Grid Point Forces (e.g., freebody loads, interface loads)
   - signficant speedups

 - [BDF](http://pynastran-git.readthedocs.io/en/latest/quick_start/bdf_demo.html)
   - 312 cards supported
   - faster node transforms using:

     ```python
     >>> icd_transform, icp_transform, xyz_cp, nid_cp_cd = get_displacement_index_xyz_cp_cd(dtype='float64, sort_ids=True)
     >>> xyz_cid0 = transform_xyzcp_to_xyz_cid(xyz_cp, icp_transform, cid=0, in_place=False)
     ```

   - simplified card adding
     ```python
     >>> model.add_grid(nid, xyz=[4.,5.,6.], comment='nid, cp, x, y, z')
     ```

- comments can now be created without worrying about `$` signs

     ```python
     >>> model.add_card(['GRID', 10, None, 4.0, 5.0, 6.0], comment='GRID comment\ngrid,nid,cp,x,y,z')
     ```
     ```
     $GRID comment
     $grid,nid,cp,x,y,z
     GRID,10,,4.0,5.0,6.0
     ```

   - unit conversion

- [GUI](http://pynastran-git.readthedocs.io/en/latest/quick_start/gui.html)
   - buttons for picking, rotation center, distance
   - PyQt5 support
   - QScintilla & pygments support for scripting code editor

- [Matlab](http://pynastran-git.readthedocs.io/en/latest/quick_start/matlab.html) integration
   - pyNastran works with Matlab 2014a+

### pyNastran v0.8.0 has been released (8/21/2016)

|  Version  | Docs  | Status |
| :--- 	  | :--- 	  | :--- 	  |
|  [v0.8.0](https://github.com/SteveDoyle2/pyNastran/releases) | [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=v0.8)](http://pynastran-git.readthedocs.io/en/v0.8.0/?badge=v0.8.0) | [![Build Status](https://img.shields.io/travis/SteveDoyle2/pyNastran/v0.8.svg)](https://travis-ci.org/SteveDoyle2/pyNastran) [![Coverage Status](https://img.shields.io/coveralls/SteveDoyle2/pyNastran/v0.8.svg)](https://coveralls.io/github/SteveDoyle2/pyNastran?branch=v0.8) |

<!---

[![Build status](https://ci.appveyor.com/api/projects/status/1qau107h43mbgghi/branch/master?svg=true)](https://ci.appveyor.com/project/SteveDoyle2/pynastran/branch/master)

[![Build status](https://ci.appveyor.com/api/projects/status/1qau107h43mbgghi?svg=true)](https://ci.appveyor.com/project/SteveDoyle2/pynastran)

[![Build Status](https://travis-ci.org/SteveDoyle2/pyNastran.png)](https://travis-ci.org/SteveDoyle2/pyNastran)
[![Coverage Status](https://coveralls.io/repos/github/SteveDoyle2/pyNastran/badge.svg?branch=master)](https://coveralls.io/github/SteveDoyle2/pyNastran?branch=master)
[![Requirements Status](https://img.shields.io/requires/github/SteveDoyle2/pyNastran/master.svg)](https://requires.io/github/SteveDoyle2/pyNastran/requirements/?branch=master)	|
--->

The following are some new features from the v0.8 release.

 - BDF
   - unicode support
   - 278 cards supported
   - simplifed method to add cards : `grid = GRID(nid, cp, xyz, cd, ps, seid, comment)`
   - `verify=False` in `read_bdf(...)` will disable checks to allow reading in bad decks
   - lots of optimization work
   - bdf equivalence, renumber, deck merging
   - element quality
 - OP2
   - ~500x faster than v0.7.2
     - non-vectorized OP2 option has been removed
   - matrix support
   - improved SOL 200 support
   - transform displacement/eigenvectors/spc/mpc/applied loads to global coordinate system
   - transform stresses/forces to material coordinate system
   - geometry can be read directly from op2 (not perfect, but when it works, it's much faster)
   - [Pandas](http://pandas.pydata.org/) DataFrame support for use in the [Jupyter/iPython](http://jupyter.org/index.html) Notebook
 - GUI
   - additional results
      - multiple OP2s
      - deflection plots
      - area, max interior angle, skew angle, aspect ratio, taper ratio results
   - SOL 200 support
   - improved legend
   - custom user (nodal/centroidal) results
   - aero models now support sideslip coordinate systems
 - OP4
  - fixed sparse ASCII BIGMAT bug


### pyNastran v0.7.2 has been Released (4/25/2015)

|  Version  | Docs  |
| :--- 	  | :--- 	  |
|  [v0.7.2](https://github.com/SteveDoyle2/pyNastran/releases) |  [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=v0.7.2)](http://pynastran-git.readthedocs.io/en/v0.7.2/?badge=v0.7.2) |

Highlights:
 * OP2
   * superelement support
   * vectorized support (uses much less memory; Element Forces not vectorized yet)
   * additional results (e.g. grid point weight, eigenvalues)
   * `PARAM,POST,-2` support
   * catching of most FATAL errors without needing to read the F06
 * BDF
   * 238 BDF cards
   * large field format and double precision writing
 * GUI
   * much improved GUI with transient support (real only), a results sidebar, logging, and scripting support
 * Other
   * additional readers/converters to/from various other formats (e.g. STL, Cart3d, Panair) as well as GUI support
   * autogenerated online documentation for pyNastran using [readthedocs](https://rwww.readthedocs.org) and [Sphinx](http://sphinx-doc.org/)

### Version 0.6.1 has been released (6/2013)
**Version 0.6** improves BDF reading.  The reader is more robust and also requires proper BDF field formatting (e.g. a integer field can't be a float).  Additionally, cards also have a comment() method.

Marcin Gąsiorek participated in the latest pyNastran under the European Space Agency's (ESA) "Summer of Code In Space" [SOCIS](http://sophia.estec.esa.int/socis2012/?q=node/13) program.  The program provides a stipend to students to work on open-source projects.
He did a great job of simplifying code and creating nicer documentation.
