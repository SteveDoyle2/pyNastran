
|  Version	| Download  | Docs  | Status |
| :-------:	| :--- 	  | :--- 	  | :--- 	  |
|  v0.7.2 	| [Download] (https://github.com/SteveDoyle2/pyNastran/releases) |  [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=v0.7.2)](http://pynastran-git.readthedocs.io/en/v0.7.2/?badge=v0.7.2) | 
|  v0.8.0 	| [![PyPi Version](https://img.shields.io/pypi/v/pynastran.svg)](https://pypi.python.org/pypi/pyNastran) | [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=v0.8)](http://pynastran-git.readthedocs.io/en/v0.8.0/?badge=v0.8.0) | [![Build Status](https://img.shields.io/travis/SteveDoyle2/pyNastran/v0.8.svg)](https://travis-ci.org/SteveDoyle2/pyNastran) [![Coverage Status](https://img.shields.io/coveralls/SteveDoyle2/pyNastran/v0.8.svg)](https://coveralls.io/github/SteveDoyle2/pyNastran?branch=v0.8) | 
|   Master	|        | [![Documentation Status](https://readthedocs.org/projects/pynastran-git/badge/?version=latest)](http://pynastran-git.readthedocs.io/en/latest/?badge=latest) | [![Linux Status](https://img.shields.io/travis/SteveDoyle2/pyNastran/master.svg)](https://travis-ci.org/SteveDoyle2/pyNastran) [![Windows Status](https://ci.appveyor.com/api/projects/status/1qau107h43mbgghi/branch/master?svg=true)](https://ci.appveyor.com/project/SteveDoyle2/pynastran) [![Coverage Status](https://img.shields.io/coveralls/SteveDoyle2/pyNastran/master.svg)](https://coveralls.io/github/SteveDoyle2/pyNastran?branch=master) | 

<!---

[![Build status](https://ci.appveyor.com/api/projects/status/1qau107h43mbgghi/branch/master?svg=true)](https://ci.appveyor.com/project/SteveDoyle2/pynastran/branch/master)

[![Build status](https://ci.appveyor.com/api/projects/status/1qau107h43mbgghi?svg=true)](https://ci.appveyor.com/project/SteveDoyle2/pynastran)

[![Build Status](https://travis-ci.org/SteveDoyle2/pyNastran.png)](https://travis-ci.org/SteveDoyle2/pyNastran)
[![Coverage Status](https://coveralls.io/repos/github/SteveDoyle2/pyNastran/badge.svg?branch=master)](https://coveralls.io/github/SteveDoyle2/pyNastran?branch=master)
[![Requirements Status](https://img.shields.io/requires/github/SteveDoyle2/pyNastran/master.svg)](https://requires.io/github/SteveDoyle2/pyNastran/requirements/?branch=master)	|
--->

### v0.8.0

[Download GUI](https://sourceforge.net/projects/pynastran/files/?source=navbar)

<!--- 
[Doumentation](http://pynastran-git.readthedocs.org/en/latest/index.html)
--->

Also, check out the:
  * [Wiki](https://github.com/SteveDoyle2/pynastran/wiki)
  * [Discussion forum](http://groups.google.com/group/pynastran-discuss)
  * [Developer forum](http://groups.google.com/group/pynastran-dev)

for more detailed information.

<!--- this isn't setup... -->
<!--- http://stevedoyle2.github.io/pyNastran/ --->

# Overview

pyNastran is an interface library to the various Nastran file formats (BDF, OP2, OP4).  Using the BDF interface, you can read/edit/write Nastran geometry without worrying about field formatting.  Many checks are also performed to verify that your model is correct.  Using the OP2 interface, you can read very large result files very quckly and very efficiently.  Additionally, you can also extract a subset of the result data and write F06 result files.

Using the pyNastran GUI, you can read in Nastran models and quickly view results for a model.  While it's no FEMAP/Patran, it can replace many tasks that would otherwise require a commercial program.

![GUI](https://github.com/SteveDoyle2/pynastran/blob/v0.8/pyNastran/gui/images/caero.png)


<!--- Ripped off meshio  --->
<!--- [![Build Status](https://travis-ci.org/SteveDoyle2/pyNastran.svg?branch=master)](https://travis-ci.org/SteveDoyle2/pyNastran)  --->
<!--- [![codecov.io](https://codecov.io/github/SteveDoyle2/pyNastran/coverage.svg?branch=master)](https://codecov.io/github/SteveDoyle2/pyNastran?branch=master)  --->

<!--- ## pyNastran v0.8.0 has NOT been released (8/21/2016)   --->
<!--- [Download pyNastran v0.8] (https://github.com/SteveDoyle2/pyNastran/releases)  --->

# News

### pyNastran v0.9 progress (12/1/2016)
 - [OP2](http://pynastran-git.readthedocs.io/en/latest/quick_start/op2_demo.html)
   - preliminary random results (ATO/CRM/PSD/RMS/NO)
   - improved geometry support
   - MONPNT1/MONPNT3 reading

 - F06 scripts
   - added a preliminary flutter (SOL 145) parser
     - supports multiple subcases
     - PK and PKNL methods supported
     - `plot_Vg_Vf(...)`, `plot_Vg(...)`, `plot_root_locus(...)`
     - input/output units
     - mode switching not fixed yet

 - [BDF](http://pynastran-git.readthedocs.io/en/latest/quick_start/bdf_demo.html)
   - 297 cards supported
   - faster node transforms using:

     ```python
     >>> icd_transform, icp_transform, xyz_cp, nid_cp_cd = get_displacement_index_xyz_cp_cd(dtype='float64, sort_ids=True)
     >>> xyz_cid0 = transform_xyzcp_to_xyz_cid(xyz_cp, icp_transform, cid=0, in_place=False)
     ```

   - simplified card adding (**card defaults are still a work in progress**)
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
   - buttons for picking, rotation center
   - PyQt5 support
   - QScintilla & pygments support for scripting code editor

- [Matlab](http://pynastran-git.readthedocs.io/en/latest/quick_start/matlab.html) integration
   - pyNastran works with Matlab 2014a+

### pyNastran v0.8.0 has been released (8/21/2016)

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

Highlights:
 * OP2
   * superelement support
   * vectorized support (uses much less memory; Element Forces not vectorized yet)
     - this is the standard in v0.8
   * additional results (e.g. grid point weight, eigenvalues)
   * `PARAM,POST,-2` support
   * catching of most FATAL errors without needing to read the F06
 * F06
   * removed
 * BDF
   * 238 BDF cards
   * large field format and double precision writing
 * GUI
   * much improved GUI with transient support (real only), a results sidebar, logging, and scripting support
 * Other
   * additional readers/converters to/from various other formats (e.g. STL, Cart3d, Panair) as well as GUI support
   * autogenerated online documentation for pyNastran using [readthedocs](https://rwww.readthedocs.org) and [Sphinx](http://sphinx-doc.org/)

Most op2 object were changed in order to eliminate errors, and be more consistent.  For example, `plateStress` has been replaced by `ctria3_stress`, `cquad4_stress`, `ctria6_stress`, etc.  Also, plate centroids now have a `node_id` of `0`.  This greatly simplifies F06 writing and vectorized data extraction.

### Version 0.6.1 has been released (6/2013)
**Version 0.6** improves BDF reading.  The reader is more robust and also requires proper BDF field formatting (e.g. a integer field can't be a float).  Additionally, cards also have a comment() method.

Marcin Gąsiorek participated in the latest pyNastran under the European Space Agency's (ESA) "Summer of Code In Space" [SOCIS](http://sophia.estec.esa.int/socis2012/?q=node/13) program.  The program provides a stipend to students to work on open-source projects.
He did a great job of simplifying code and creating nicer documentation.
