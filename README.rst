+--------------+-----------------------------------------------------------------------+-------------+--------------+
| **Version**  | **Docs**                                                              | **Status**  | **Coverage** |
+--------------+-----------------------------------------------------------------------+-------------+--------------+
| |PyPi13|_    | `Docs <https://pynastran-git.readthedocs.io/en/1.3/>`_                | |Travis13|_ | |Codecov13|_ |
+--------------+-----------------------------------------------------------------------+-------------+--------------+
| Main         | `Docs <http://pynastran-git.readthedocs.io/en/latest/?badge=latest>`_ | |TravisMa|_ | |CodecovMa|_ |
+--------------+-----------------------------------------------------------------------+-------------+--------------+

.. |PyPi13| image:: https://img.shields.io/pypi/v/pynastran.svg
.. _PyPi13: https://pypi.python.org/pypi/pyNastran

.. comments
   #-----------------------------------------------------------------

.. |Travis11| image:: https://img.shields.io/travis/SteveDoyle2/pyNastran/v1.1.svg
.. _Travis13: https://travis-ci.org/SteveDoyle2/pyNastran

.. |Travis12| image:: https://img.shields.io/travis/SteveDoyle2/pyNastran/1.2.svg
.. _Travis12: https://travis-ci.org/SteveDoyle2/pyNastran

.. |Travis13| image:: https://img.shields.io/travis/SteveDoyle2/pyNastran/1.3.svg
.. _Travis13: https://travis-ci.org/SteveDoyle2/pyNastran

.. |TravisMa| image:: https://img.shields.io/travis/SteveDoyle2/pyNastran/main.svg
.. _TravisMa: https://travis-ci.org/SteveDoyle2/pyNastran

.. comments
   #-----------------------------------------------------------------

.. |Codecov11| image:: https://img.shields.io/coveralls/SteveDoyle2/pyNastran/1.1.svg
.. _Codecov11: https://coveralls.io/github/SteveDoyle2/pyNastran?branch=1.2

.. |Codecov12| image:: https://img.shields.io/coveralls/SteveDoyle2/pyNastran/1.2.svg
.. _Codecov12: https://coveralls.io/github/SteveDoyle2/pyNastran?branch=1.2

.. |Codecov13| image:: https://codecov.io/gh/SteveDoyle2/pyNastran/branch/1.3/graph/badge.svg
.. _Codecov13: https://codecov.io/gh/SteveDoyle2/pyNastran/branch/1.3

.. |CodecovMa| image:: https://codecov.io/gh/SteveDoyle2/pyNastran/branch/main/graph/badge.svg
.. _CodecovMa: https://codecov.io/gh/SteveDoyle2/pyNastran/branch/main

.. comments
   #-----------------------------------------------------------------

v1.3.3
------

`Download GUI <https://sourceforge.net/projects/pynastran/files/?source=navbar>`_ (latest is from 2020/6/28)

Also, check out the:

* `Discussion forum <http://groups.google.com/group/pynastran-discuss>`_ (questions about the latest release)

* `Developer forum <http://groups.google.com/group/pynastran-dev>`_ (questions about the main branch)

* `Docs <https://pynastran-git.readthedocs.io/en/1.3/>`_

for more detailed information.

Overview
========

pyNastran is an interface library to the various Nastran file formats (BDF, OP2, OP4).
Using the BDF interface, you can read/edit/write Nastran geometry without worrying about
field formatting.  Many checks are also performed to verify that your model is correct.
Using the OP2 interface, you can read large result files quickly and efficiently.
Additionally, you can also extract a subset of the result data and write OP2/F06 result
files.  For a more detailed list of features, see:

* `Features <https://pynastran-git.readthedocs.io/en/1.3/quick_start/features.html#overview>`_

Using the pyNastran GUI, you can read in Nastran models and quickly view results for a model.
While it's no FEMAP/Patran, it can replace many tasks that would otherwise require a
commercial program.

.. image:: https://github.com/SteveDoyle2/pynastran/blob/main/pyNastran/gui/images/caero.png

News
====

Check out this blog post covering `Flutter Analysis in pyNastran <https://www.m4-engineering.com/flutter-analysis-with-pynastran/>`_

Release Notes
=============

pyNastran v1.3.3 has been released (2020/6/28)
----------------------------------------------

v1.3.3 (2020/6/28)
------------------
This is a bug fix only release outside of:

     .. code-block:: python

        >>> subcase.add_set_from_values(set_id, values)

which was overly complicated to do before.


Programmatics:

* Supports Python 3.7 and 3.8

* GUI is compatible with PyQt5 and PySide2 and VTK 7-9

* improved testing of old versions of packages

* updated nptyping requirements

* support for NX 2019.0, 2019.1

See pyNastran `v1.3.3 release notes <https://github.com/SteveDoyle2/pyNastran/releases/tag/v1.3.3>`_ for details on the 54 bug fixes.  Not too bad for 2.5 months!


v1.3.2 (2020/4/8)
-----------------

With Python 2 now officially dead, it's time for a new killer feature to get the last few people to switch.

There is now support for writing OP2 files!  They're difficult to create from scratch,
but modifying an existing one isn't difficult.  This includes geometry support.

In addition, many new OP2 results have been added.  Modern NX Nastran should work much better.

Programmatics:

* Supports Python 3.7 and 3.8

* Dropping Python 2.7 and 3.6 support

* GUI is compatible with PyQt5 and PySide2

* improved testing of old versions of packages

BDF:
 * enhancements

   * 405 cards supported (up from 373)

   * improved mesh utilities

     * bdf mirror, bdf convert, bdf export_caero_mesh

     * additional `bdf scale`, which takes arbitrary mass, length, time, pressure, or velocity scale factors (3 of 5) to apply a scaling law

   - added deck guessing with punch=None

OP2:
 - enhancements:

   - OP2 write support

   - fixed most pandas deprecation warnings/added support for pandas 0.25

   - much improved NX 2019.2 support

 - minor enhancements:

   - more hdf5 results

   - a few more cards supported in the OP2 geometry reader

   - preliminary NX 64-bit support

   - more MSC versions supported

   - sped up R1TABRG (optimization) reading

   - various geometry cards added

   - supports more PARAM reading

 - new results:

   - composite failure indicies (OEFIT)

   - sensitivity support (DSCMCOL)

   - Cambpell diagrams (CDDATA)

   - eigenvectors (BOPHIGF)

   - grid point forces SORT2 (OGPFB2)

   - stress/strain/force

     - random stress/strain (OESVM1/2, OESVM1C, OSTRVM1/2, OSTRVM2, OESXRM1C, OESXNO1, OESXNO1C)

     - real/complex stress/strain/force output for centroidal CQUADR/CTRIAR

     - complex CBEAR forces

     - real CFAST, CWELD forces

     - nonlinear cbush stress/strain/force

   - other

     - PSDs

     - optimization

       - convergence table

       - design variables

       - weight, displacement, stress, strain, force, composite stress, composite strain, fractional mass response

   - SOL 401/402 results:

     - eigenvalue_fluid (LAMAF, LAMAS)

     - eigenvectors (BOPHIGF, BOPHIGS)

     - temperature (OTEMP1)

     - solution set

       - results: displacement, velocity, acceleration, eigenvectors:

       - tables: OUXY1, OUXY2, OPHSA

GUI:

- enhancements:

  - partial custom results many now be loaded

  - animation now supports complex fringe

  - result case description now shows the mode/time/frequency

  - result case description now shows the min/max value as well as the location

  - map centroidal to nodal option

  - adding node/element highlight menu

  - adding node/element mark/label menu

  - result case description now shows the min/max value as well as the location

  - result case description now shows the mode/time/frequency

- minor enhancements:

  - added export result option to right-click menu

  - more keyboard shortcuts

- nastran:

  - geometry:

    - added preferences for geometry/results to speed up model loading

    - element & material coordinate systems

  - results

    - added acoustic displacements (OUG1F)

    - real/complex stress/strain/force results

      - plate by upper/lower, composite plates by ply, bars, rods, springs, cbush, cdamp

    - fractional mass response

F06:
 - KE support for plot_145

Many bug fixes and a few more details on features, can be found in the `Release Notes <https://github.com/SteveDoyle2/pyNastran/blob/1.3/releaseNotes.md>`_.


v1.3.1/v1.3.0 (2020/4/8)
------------------------

This result has been superseded by 1.3.2.  The PyPi page was fixed.

v1.2.1 (2019/5/24)
------------------

I keep saying to myself there's not much to add, but Nastran is huge.  Beyond HDF5
support in the BDF, I'm a huge fan of the new ability to keep track of which include file a
card came from and write it as a separate file.  It's limited in usefulness, but very handy
in certain cases.  There's a new (still preliminary) superelement capability.  I'm far
from a superelement expert, but it's probably useful.

The OP2 reader now supports SORT2 along with much improved random results reading.
If you're using 60+ GB OP2s, you probably have had issues with RAM usage in the past.
With the new ability to dump the OP2 directly to HDF5, this should not be an as much of
an issue.  It's not 100% implemented, so let me know if you need it for another result.

Regarding the GUI, there are also some new features.  Groups work a bit better, but aren't
quite perfect.  Logging has been dramatically sped up so the GUI loads faster and you can
load Nastran models even faster if you disable additional results (e.g., element quality).

Finally, Python 2.7 is end of life.  Numpy, scipy, and matplotlib have all dropped
Python 2.7 support.  It's time for pyNastran to as well.  The OP2 reader is 30% faster in
Python 3.6+ than Python 2.7, so it's not all bad!

Programmatics:
 - This is the last Python 2.7 release.
 - Dropping Python 3.5 support
 - Supports Python 2.7, 3.6-3.7
 - GUI is compatible with PyQt4/PyQt5 as well as PySide/PySide2
 - improved testing of old versions of packages

BDF:
 - 373 cards supported (up from 343)

 - added abiltity to write models to separate include files

     .. code-block:: python

        >>> model = BDF()
        >>> model.read_bdf(bdf_filename, save_file_structure=True)

        out_filenames = {
            bdf_filename : bdf_filename_new,
            include_filename : include_filename_new,
        }
        >>> model.write_bdfs(out_filenames, relative_dirname=None, is_windows=None)
        >>> ifile = model.grids[1].ifile

 - HDF5 import/export

     .. code-block:: python

        >>> model = read_bdf(bdf_filename)
        >>> model.export_hdf5_filename(hdf5_filename)
        >>> model_new = OP2()
        >>> model_new.load_hdf5_filename(hdf5_filename)

 - preliminary superelement support

     .. code-block:: python

        >>> model.read_bdf(bdf_filename)
        >>> model.superelement_models[1].nodes

OP2:
 - reorganization of random op2 results into op2.results.psd (or ato, no, crm, rms) to aide in finding data
 - reorganization of op2 class to reduce number of functions in the object.  This affects any custom table reading.
 - improved optimzation response reading
 - limited SORT2 support
 - fixed CD transformation bug for BOUGV1 and BOPHIG1 tables
 - Improved HDF5 export/import support (e.g., matrices, random results)

 - Can optionally save directly to HDF5 instead of numpy (limited).
 - Loading OP2s to an HDF5 file to decrease memory usage

    .. code-block:: python

       >>> op2_model = OP2()
       >>> op2_model.load_as_h5 = True
       >>> op2_model.read_op2(op2_filename)

OP2Geom:
 - HDF5 support
 - reading EQEXIN/S, GPT, GPDT, CSTM/S tables (recovery of nodes & coordinate with OP2Geom)
 - fixed theta/mcid reading for CTRIA3/CQUAD4
 - fixed CQUAD8 bug

GUI:
 - sped up HTML logging
 - much improved groups menu
 - options for Nastran in preferences menu to speed up loading/limit memory usage
 - pyNastran BDF pickle reading
 - pyNastran OP2 HDF5 reading (not MSC's format)
 - visualization when pickling nodes/elements
 - min/max labels
 - highlight menu
 - Patran-style colors
 - custom force vectors
 - AVL support


Known issues:
 - Transient Pandas Dataframes will fail for newer versions of numpy/pandas.

Older Releases
--------------
 - pyNastran `v1.2.0 release notes <https://github.com/SteveDoyle2/pyNastran/releases/tag/v.1.2.0>`_ (2019/5/21)

 - pyNastran `v1.1.0 release notes <https://github.com/SteveDoyle2/pyNastran/releases/tag/v1.1.0>`_ (2018/6/26)

 - pyNastran `v1.0.0 release notes <https://github.com/SteveDoyle2/pyNastran/releases/tag/v1.0.0>`_ (2017/5/25)

 - pyNastran `v0.8.0 release notes <https://github.com/SteveDoyle2/pyNastran/releases/tag/v0.8.0>`_ (2016/8/21)

 - pyNastran `v0.7.2 release notes <https://github.com/SteveDoyle2/pyNastran/releases/tag/v0.7.2>`_ (2015/4/25)
