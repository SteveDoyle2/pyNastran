
|  Version  | Docs  | Status |
| :--- 	  | :--- 	  | :--- 	  |
|  [![PyPi Version](https://img.shields.io/pypi/v/pynastran.svg)](https://pypi.python.org/pypi/pyNastran) | [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=v1.0)](http://pynastran-git.readthedocs.io/en/v1.0.0/?badge=v1.0.0) | [![Build Status](https://img.shields.io/travis/SteveDoyle2/pyNastran/v1.0.svg)](https://travis-ci.org/SteveDoyle2/pyNastran) [![Coverage Status](https://img.shields.io/coveralls/SteveDoyle2/pyNastran/v1.0.svg)](https://coveralls.io/github/SteveDoyle2/pyNastran?branch=v1.0) |
|   Master | [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=latest)](http://pynastran-git.readthedocs.io/en/latest/?badge=latest) | [![Linux Status](https://img.shields.io/travis/SteveDoyle2/pyNastran/master.svg)](https://travis-ci.org/SteveDoyle2/pyNastran) [![Windows Status](https://ci.appveyor.com/api/projects/status/1qau107h43mbgghi/branch/master?svg=true)](https://ci.appveyor.com/project/SteveDoyle2/pynastran) ![Coverage Status](https://coveralls.io/repos/github/SteveDoyle2/pyNastran/badge.svg?branch=master) | 



<!---
[![codecov](https://codecov.io/gh/SteveDoyle2/pyNastran/branch/master/graph/badge.svg)](https://codecov.io/gh/SteveDoyle2/pyNastran) 

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
  * [Discussion forum](http://groups.google.com/group/pynastran-discuss) (intended for questions about the latest release)
  * [Developer forum](http://groups.google.com/group/pynastran-dev) (intended for questions about the master branch)
  * [Docs](http://pynastran-git.readthedocs.io/en/latest/?badge=latest) (the version isn't quite right, but it's close)

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

It's been roughly a year and 100 tickets closed since the last version, so it's probably
time for another release!  I'm pleased to announce that the continuing problems 
of up-to-date documentation will hopefully be a thing of the past.  xxx 
has graciously offered to host the documentation.  The pyNastranGUI exe (previously 
hosted on Sourceforge) will also be hosted there.  Outside of that, it's the same 
open-source project and will still be on Github.

Regarding features, the focus has again been on robustness and testing.  There has been
a 10% increase in the testing coverage (the same as v0.8 to v1.0).

There are a few changes (mainly in the BDF) though.  The GUI now also supports
PyQt4, PyQt5, and Pyside with the same API, which helps reduces the restriction on
licensing.

Programmatics
 - Dropping Python 3.4 support
 - dropping VTK 5/6 suppoprt for the GUI

BDF:
 - currently 340 cards supported
 - cross-referencing is now more straightforward to new users (much of v1.0 works using the `_ref` option)
   - ``*_ref`` attributes are cross-referenced
      - ``element.nodes`` is not cross-referenced
      - ``element.nodes_ref`` is cross-referenced
 - pickling to reload your deck ~5x faster
 - decreased time required for Case Control Deck with large SETs and many load cases
 - improved optimization checks

OP2:
 - HDF5 export/import support
 - pandas support for matrices
 - couple more results vectorized (e.g., complex strain energy, DMIG strain energy, some forces)
 - grid_point_stressses supported (disabled since v0.7)
 - fixed sparse matrices being stored as dense matrices

GUI:
 - can now mix and match fringe/displacement/vector results (e.g., max principal stress shown on a displaced model)
 - improved animation menu
    - in gui animation
    - more animation profiles
 - bar profile visualzation
 - nominal geometry (useful for deflection plots)
 - improved optimization support
 - improved picking display
 - better PSHELL/PCOMP distinction

Known issues:
 - Transient Pandas Dataframes will fail in newer versions of numpy/pandas.  If anyone knows how to use a MultiIndex,
   this is probably pretty easy to fix.
 
--->
  
# News

### pyNastran v1.0.0 has been released (5/25/2017)
This is a major release.  The focus this time has been on robustness and testing.
Hopefully, it shows.  The software has also been relicensed to be **BSD-3**, which
is a more permissive license and is the same one that numpy, scipy, and
matplotlib use.

Unfortunately, the GUI is more complicated.
 - For open source projects : GPL 2/3
 - For companies that pay a license to Riverbank : proprietary
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
|  [v0.8.0](https://github.com/SteveDoyle2/pyNastran/releases) | [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=v0.8)](http://pynastran-git.readthedocs.io/en/v0.8.0/?badge=v0.8.0) | [![Build Status](https://img.shields.io/travis/SteveDoyle2/pyNastran/v0.8.svg)](https://travis-ci.org/SteveDoyle2/pyNastran)  |

<!---

[![Build status](https://ci.appveyor.com/api/projects/status/1qau107h43mbgghi/branch/master?svg=true)](https://ci.appveyor.com/project/SteveDoyle2/pynastran/branch/master)

[![Build status](https://ci.appveyor.com/api/projects/status/1qau107h43mbgghi?svg=true)](https://ci.appveyor.com/project/SteveDoyle2/pynastran)

[![Build Status](https://travis-ci.org/SteveDoyle2/pyNastran.png)](https://travis-ci.org/SteveDoyle2/pyNastran)
[![Requirements Status](https://img.shields.io/requires/github/SteveDoyle2/pyNastran/master.svg)](https://requires.io/github/SteveDoyle2/pyNastran/requirements/?branch=master)	|
--->

See [v0.8.0](https://github.com/SteveDoyle2/pyNastran/releases/tag/v0.8.0) for information regarding enhancements.

### pyNastran v0.7.2 has been Released (4/25/2015)

|  Version  | Docs  |
| :--- 	  | :--- 	  |
|  [v0.7.2](https://github.com/SteveDoyle2/pyNastran/releases) |  [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=v0.7.2)](http://pynastran-git.readthedocs.io/en/v0.7.2/?badge=v0.7.2) |

See [v0.7.2](https://github.com/SteveDoyle2/pyNastran/releases) for information regarding enhancements.

### Version 0.6.1 has been released (6/2013)
**Version 0.6** improves BDF reading.  The reader is more robust and also requires proper BDF field formatting (e.g. a integer field can't be a float).  Additionally, cards also have a comment() method.

Marcin Gąsiorek participated in the latest pyNastran under the European Space Agency's (ESA) "Summer of Code In Space" [SOCIS](http://sophia.estec.esa.int/socis2012/?q=node/13) program.  The program provides a stipend to students to work on open-source projects.
He did a great job of simplifying code and creating nicer documentation.
