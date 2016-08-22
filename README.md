
### v0.8.0
[![PyPi Version](https://img.shields.io/pypi/v/pynastran.svg)](https://pypi.python.org/pypi/pyNastran)
[![PyPi Downloads](https://img.shields.io/pypi/dm/pynastran.svg)](https://pypi.python.org/pypi/pyNastran)

<!--- [![Coverage Status](https://img.shields.io/coveralls/SteveDoyle2/pyNastrane/master.svg)](https://coveralls.io/github/SteveDoyle2/pyNastran?branch=master)   --->

[Doumentation](http://pynastran-git.readthedocs.io/en/v0.8/)

[Download pyNastran] (https://github.com/SteveDoyle2/pyNastran/releases)

[Download GUI](https://sourceforge.net/projects/pynastran/files/?source=navbar)


### v0.7.2

<!--- [![Github Downloads](	https://img.shields.io/github/downloads/SteveDoyle2/pyNastran/latest/total.svg)  --->
<!--- [![All Downloads](https://img.shields.io/github/downloads/SteveDoyle2/pyNastran/total.svg)  --->
<!--- [![Total PyPi Downloads](https://img.shields.io/github/downloads/atom/atom/latest/pynastran.svg)]  --->
<!--- [![v0.7.2 Downloads](https://img.shields.io/github/downloads/atom/atom/v0.7.2/total.svg)]  --->

[Documentation] (http://pynastran-git.readthedocs.org/en/v0.7.2/index.html)

### Master

[![Build Status](https://travis-ci.org/SteveDoyle2/pyNastran.png)](https://travis-ci.org/SteveDoyle2/pyNastran)

[Doumentation](http://pynastran-git.readthedocs.org/en/latest/index.html)

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

### pyNastran v0.8.0 has been released (8/21/2016)

  [Download pyNastran v0.8] (https://github.com/SteveDoyle2/pyNastran/releases)
  
  [Download GUI](https://sourceforge.net/projects/pynastran/files/?source=navbar)

The following are some new features from the v0.8 release.

 - BDF
   - unicode support
   - 278 cards supported
   - simplifed method to add cards : `grid = GRID(nid, cp, xyz, cd, ps, eid, comment)`
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
   - imporoved legend
   - custom user (nodal/centroidal) results
   - aero models now support sideslip coordinate systems
 - OP4
  - fixed sparse ASCII BIGMAT bug


### pyNastran v0.7.2 has been Released (4/25/2015)

[Download pyNastran v0.7] (https://github.com/SteveDoyle2/pyNastran/releases)

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
