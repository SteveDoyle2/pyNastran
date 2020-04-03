Information
===========
If you like it/want to help leave some feedback:
 - [Dicussion Page](https://groups.google.com/forum/#!forum/pynastran-discuss)
 - [Developer Dicussion](https://groups.google.com/forum/#!forum/pynastran-dev)

If you have a bug/want a new feature or card, leave some feedback on the [Issue Tracker](https://github.com/SteveDoyle2/pyNastran/issues)

=============
Release Notes
=============

------------------------
v1.3.0 released 2020/4/8
------------------------
With Python 2 now officially dead, it's time for a new feature to encourage people to switch.

There is now support for writing OP2 files!  They're difficult to create from scratch, 
but modifying an existing is easy.  This also includes geometry support.

In addition, many new OP2 results have been added.  Modern NX Nastran should work much better.

Programmatics:
 - Supports Python 3.7 and 3.8
 - Dropping Python 2.7 and 3.6 support
 - GUI is compatible with PyQt5 and PySide2
 - improved testing of old versions of packages
 
BDF:
 - enhancements
   - 405 cards supported (up from 373)
   - improved mesh utilities
      - bdf mirror, bdf convert, bdf export_caero_mesh
      - additional `bdf scale`, which takes arbitrary mass, length, time, pressure, or velocity scale factors (3 of 5) to apply a scaling law
   - added deck guessing with punch=None
   - better parser for test_bdf
   - PCOMPG now supports DVPRELx
   - added get_position for SPOINT
   - added inertia for PBARL/L
   - bdf mirror:
     - added PLOTEL, RROD, RBAR, RBAR1, RBE1, AELIST, AESURF support
   - bdf_equivalence:
     - method='new' no longer requires neq_max
     - use method='old' for old behavior (deprecated)
 - bug fixes:
 - new cards:
   - DTI besides UNITS are now supported
   - PBEAM3, TEMPB3
 - minor enhancements:
   - handling ;; in DEQATN
   - support for SURF cid=None in _get_forces_moments_array
   - added check for number of loads = nscale factors for LOAD
   - more PARAM validation
- changes:
   - warning on RBAR dof check instead of RuntimeError

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
      - XSOP2DIR
      - PSDs
      - optimization
        - convergence table
        - design variables
        - weight response
        - displacement response
        - stress response
        - strain response
        - force response
        - composite stress response
        - composite strain response
        - fractional mass response
    - SOL 401/402 results:
      - eigenvalue_fluid (LAMAF, LAMAS)
      - eigenvectors (BOPHIGF, BOPHIGS)
      - temperature (OTEMP1)
      - solution set
        - results: displacement, velocity, acceleration, eigenvectors:
        - tables: OUXY1, OUXY2, OPHSA
 - API:
   - reorganized output results to use op2_results object to simply interface
     - backwards compatibility for stress, strain, forces, strain energy
     - model.cbar_force -> model.op2_results.force.cbar_force
     - model.cbar_stress -> model.op2_results.stress.cbar_stress
     - model.cbar_strain -> model.op2_results.strain.cbar_strain
     - model.cbar_strain_energy -> model.op2_results.strain_energy.cbar_strain_energy
 - TODO: linear/nonlinear hyperelastic solids
 - TODO: stress transforms.  This is probably a bit of work.
 - TODO: preliminary NX 64-bit support
 - TODO: CD transforms for cylindrical/spherical displacement, velocity, acceleration, forces.  This shouldn't be terrible.

GUI:
 - enhancements:
    - Can now load custom fringe/displacement/force results with an incorrect
      number of nodes/elements.  It's assumed that the node/element id in the
      first column should be the same as the id for the model.  Thus, you just
      filter out extra nodes/elements or you set some blank nodes/elements.
       - For integer results, no masking is applied.
       - For float results, masking is applied and nan corresponds to no color.
   - animation now supports complex fringe
   - result case description now shows the min/max value as well as the location
   - result case description now shows the mode/time/frequency
   - map centroidal to nodal option
   - adding node/element highlight menu
   - adding node/element mark/label menu
   - result case description now shows the min/max value as well as the location
   - result case description now shows the mode/time/frequency
 - minor enhancements:
   - added export result option to right-click menu
   - legend supports unicode
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

Bug fixes:
 - BDF:
   - fixing bug in set_param
   - PBARL & PBEAML DBOX now has defaults for DIM5-DIM10
   - added AECOMPL, DEFORM to bdf attributes
   - added RHO to MAT10 updater
   - Better case control SET parser (can now handle 1,THRU,10,EXCEPT,7,8,
     which doesn't have spaces that makes parsing easier).  It's also orders
     of magnitude faster on large problems.
   - Fixed NDDL bug in MONPNT2
   - Added check for DVPRELx to prevent referencing property fields when another
     similar property is used (e.g., a DVPRELx references a PCOMPG, but the PCOMP
     field id numbers are used to update the card).
   - added some missing hdf5 object support
   - fixing DEQATN bug caused by using python builtin
   - fixed missing variable in mass properties error message
   - aero
     - fixed issue where CAERO1 didn't set the xyz1 and xyz4 points to arrays
       when loading an HDF5 file
     - caero2 now sets box ids
   - elements/properties
     - added I12, J to PBEAM updater; fixed error message
     - fixed Area calculation for PBARL/CHAN, HAT1, DBOX
     - fixed PBAR optimization treating I12 as I1; added J
   - convert better handles:
     - T1-T4 in quad shells
     - fixed bug in CGAP x/g0
   - bdf mirror:
     - fixed bug where max element id is less than max rigid_element id
     - fixed bug with shells that have blanks (e.g., CTRIA6, CQUAD8)
     - plane is now case insensitive and auto-sorted
     - CAERO1 now supports lspan/lchord and Cp
     - CAERO1 handles xy plane properly now (doesn't assume it's xz)

OP2:
 - geom bug fixes:
    - fixed reading bugs for QHBDY and SPCOFF
    - better identifing duplicate property id (PBCOMP -> PBEAM)
    - fixed CONVM error (it can be 6 or 7 fields, not only 6)
    - fixed ACCEL key
 - bug fixes:
    - R1TABRG now stores response_types as strings
      (instead of bytes for Python 3)
    - f06 writer now writes CDAMP3/4 names correctly

GUI:
 - bug fixes:
  - better argument handling
  - fixed coordinate system scaling bug
  - added check on highlight menu for model existance
  - fixed import for new version download menu
  - improved command line error message
  - fixed support for CAEROx models without elements
  - export_cases now supports integers

v1.2.1 (2019/5/24)
------------------
OP2:
 - fixed bug with OUGV1PAT table

v1.2.0 (2019/5/21)
------------------
Programmatics:
 - This is the last Python 2.7 release.
 - Dropping Python 3.5 support
 - Supports Python 2.7, 3.6-3.7
 - GUI is compatible with PyQt4/PyQt5 as well as PySide/PySide2
 - improved testing of old versions of packages

BDF:
 - 373 cards supported (up from 343)
 - added abiltity to write models to separate include files
     ```python
     >>> model = BDF()
     >>> model.read_bdf(bdf_filename, save_file_structure=True)

     out_filenames = {
         bdf_filename : bdf_filename_new,
         include_filename : include_filename_new,
     }
     >>> model.write_bdfs(out_filenames, relative_dirname=None, is_windows=None)
     >>> ifile = model.grids[1].ifile
     ```

 - HDF5 import/export
      ```python
     >>> model = read_bdf(bdf_filename)
     >>> model.export_hdf5_filename(hdf5_filename)
     >>> model_new = OP2()
     >>> model_new.load_hdf5_filename(hdf5_filename)
     ```

 - preliminary superelement support
     ```python
     >>> model.read_bdf(bdf_filename)
     >>> model.superelement_models[1].nodes
     ```


OP2:
 - reorganization of random op2 results into op2.results.psd (or ato, no, crm, rms) to aide in finding data
 - reorganization of op2 class to reduce number of functions in the object.  This affects any custom table reading.
 - improved optimzation response reading
 - limited SORT2 support
 - fixed CD transformation bug for BOUGV1 and BOPHIG1 tables
 - Improved HDF5 export/import support (e.g., matrices, random results)

 - Can optionally save directly to HDF5 instead of numpy (limited).
 - Loading OP2s to an HDF5 file to decrease memory usage
      ```python
     >>> op2_model = OP2()
     >>> op2_model.load_as_h5 = True
     >>> op2_model.read_op2(op2_filename)
     ```

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
 - Transient Pandas Dataframes will fail for newer versions of numpy/pandas.  If anyone knows how to use a MultiIndex,
   this is probably pretty easy to fix.

v1.1.0 (2018/6/26)
------------------
This is a major release.  As in the v1.0 release, the focus has been on testing.

The bigest change has been to cross-referencing.  All cross-referenced
objects now have an ``_ref`` attribute that is cross-referened, which
is ``None`` if the model is not cross referenced.  Objects maintain the
uncross-referenced form as well.  For a CQUAD4, instead of:
 ``model.elements[eid].nodes[0].nid``
the referenced object is:
 ``model.elements[eid].nodes_ref[0].nid``

For detail on the license for the various components, see the v1.0 release.

Quick Overview:
===============
Programmatics:
 - Added support for numpy 1.14
 - Dropping support for Python 3.4 (2.7, 3.5, 3.6 are supported)
 - Dropping support for VTK 5/6 (7/8 are supported)

BDF:
 - model may be pickled (model.save('model.obj') and model.load('model.obj')
 - simplified cross-referencing

OP2:
 - model may be pickled (model.save('model.obj') and model.load('model.obj')
 - Added support for exporting OP2 to HDF5 (uses pyNastran format, not MSC)
 - real sparse matrices take much less memory now; were being converted to dense matrices

GUI:
 - improved animations
 - improved labels


Detailed Overview:
==================
BDF:
  - New features:
   - model.get_reduced_loads(load_id)
   - model.get_reduced_dloads(dload_id)
   - card1 == card2 now supported
   - PBARL/PBEAML support the NX TUBE2 type
   - added FREQ3, FREQ5, CAERO5, PAERO5, MATT3, SPLINE3, RSSCON, OMIT1
   - more xref
   - added model.clear_attributes()
   - renumbering:
     - bdf_renumber now supports renumbers SETs and SPLINEx cards
     - SPOINTs/EPOINTs now use dictionaries to enable SPOINT/EPOINT renumbering
     - caero sub-panels ids are now renumbered
     - renumbering mapper object returned now
   - improved removed_unused card support
   - improved mirroring
   - read_bdf StringIO option now parses pyNastran header
   - subcase copying speedup (helps with SETs)
   - preliminary ZONA loading
   - added atmosphere2.make_flfacts_eas_sweep, make_flfacts_mach_sweep, and
     make_flfacts_alt_sweep with an EAS (equivalent airspeed) limiter
   - rotate_v_wa_wb for CBAR/CBEAM to determine element vectors

  - Bug fixes:
   - more add_card documentation (e.g., add_grid, add_ctria3)
   - fixed NSMADD card type (was SPCADD)
   - fixed CPLTSTN3 card type (was CTRIA3)
   - fixed CPLTSTS3 card type (was CTRIA3)
   - fixed shell MCIDs not renumbering
   - fixed FREQx renumbering crash
   - model may be pickled again (model.save('model.obj') and model.load('model.obj')
   - fixed PBARL/PBEAML DBOX error
   - fixed PositionWRT bug
   - fixed LOAD card messing up load ids after cross referencing
   - fixed TRIM default on aeqr (was 0.0/rigid; should be 1.0/elastic)
   - NLPCI now gets written when there are no other dynamic cards
   - fixed CBUSH cid=0 bug
   - fixed TABLED4/TABLEM4 stopping error

  - API changes:
   - model.Node(nid, allow_empty_nodes=False msg='') no longer supports
     allow_empty_nodes.  Use:
       model.EmptyNode(nid, msg='') instead for that
       model.Node(nid, msg='') is the new form
   - model.Nodes(nid, allow_empty_nodes=False msg='') no longer supports
     allow_empty_nodes.  Use:
       model.EmptyNodes(nid, msg='') instead for that
       model.Nodes(nid, msg='') is the new form
   - PCOMPG.validate() now checks that global ply ids are unique
   - xref_nodes_with_elements now creates a list instead of a set
     (fixes a Python 3.x bug)
   - get_MPCx_node_ids_c1 is now get_MPCx_node_ids_c1
     (was inconsistent with what it does)
   - get_MPCx_node_ids_c1 created
   - xref'd objects now use _ref globally
   - aestat.id is now aestat.aestat_id
   - aeparm.id is now aeparm.aeparm_id
   - model.add_aset1/aset (also bset/cset/qset/uset) now a consistent set of
     function arguments and call the same function.  The card will be created
     based on your data instead of necessarily what you asked for.
   - LOAD cards are now stored in model.load_combinations instead of model.loads

  - Known bugs:
   - dynamic loads cross-referencing is buggy;
     reject the cards if there is a problem
   - PBEAM defaults with ENDA are slightly incorrect.

OP2:
  - New features:
   - added model.set_additional_generalized_tables_to_read(tables) to
     create custom OP2 readers
   - added complex/average strain energy
   - save/load hdf5 support
   - EIGRL support

  - Bug fixes:
   - improved table skipping
   - fixed RealCShearForceArray f06 writing
   - fixed CEN/3, CEN/4 writing for RealPlateBilinearForceArray
   - improved geometry table reading
   - real sparse matrices take much less memory now; were being converted to dense matrices
   - added RBAR on NX vs. MSC
   - fixed RBE2 with alpha bug
   - fixed CREEP bug
   - fixed RBE3 bug
   - fixed PBCOMP bug

  - API changes:
    - xlsx exporter removed

  - Known bugs:
    - pandas fails on some decks (numpy<1.13 is fine)
    - a large number of PSOLIDs will crash the read_op2/read_op2_geom;
      use PARAM,OGEOM,NO
    - transform_gpforce_to_global doesn't work properly with cylindrical
      or spherical coordinate systems

GUI:
  - New features:
    - control surfaces now get labels (label size doesn't resize properly)
    - in-gui animation
    - delete secondary actor support
    - delete result cases support
    - improved view buttons
    - preferences menu
    - right click support on results sidebar to apply fringe/displacement/veactor results

  - Bug fixes:
    - displays control surfaces again (aesurf)
    - changing secondary actor color works again
    - fixing random crash
    - "Show/Hide CAERO panels" updates the "Edit Geometry Properties" menu
    - "Toggle CAERO Subpanels" updates the "Edit Geometry Properties" menu
    - "Toggle CONM2s" updates the "Edit Geometry Properties" menu
    - fixing Windows taskbar icon bug
    - fixing first launch bug
    - qscintilla works in pyqt5

  - Known bugs:
    - after animating a model from within the GUI, the mouse behavior changes

OP4 bug fixes:
 - fixed Python 3 bytes bug

Applictions:
 - removed due to excessively amount of unmaintained code


v1.0.0 (2017/5/25)
------------------
This is a major release.  The focus this time has been on robustness and testing.
Hopefully, it shows.  The software has also been relicensed to be BSD-3, which
is a more permissive license and is the same one that numpy, scipy, and
matplotlib use.

Unfortunately, the GUI is more restrictive and more complicated.
 - For open source projects : LGPL 2/3
 - For companies that pay a license to Riverbank : LGPL 2/3
 - For companies that don't pay a license fee : GPL 2/3
However, you may distribute an unmodified binary.

Quick Overview:
===============

Programmatics:
 - Dropping Python 3.3 support
 - Adding Python 3.6 support

BDF:
 - This version makes it much easier to programmatically create BDFs with the
   ``add_grid``, ``add_ctria3(...)``, etc. methods.  Almost every card has a
    method and many have documentation of every parameter.

OP2:
 - Added support for MATPOOL matrices and RANDOM OUG-style tables.
 - First release of read_op2_geom.  Much improved OP2 Geometry reading.

F06:
 - added a preliminary flutter (SOL 145) parser
   - supports multiple subcases
   - PK and PKNL methods supported
   - `plot_Vg_Vf(...)`, `plot_Vg(...)`, `plot_root_locus(...)`
   - input/output units
   - mode switching not fixed yet

GUI:
  - Animation plots
  - NaN colors
  - improved picking tolerancing
  - probe result by mouse picking
  - measure distance by mouse picking
  - rotation center by mouse picking
  - area picking nodes/elements
  - added basic support for Patran *.nod results
  - significant speedups

Detailed Overview:
==================
BDF:
 - added cards:
   - CTRAX3, CTRAX6, CQUADX4, CQUADX8, CIHEX2
   - PAERO4, MONPNT2, MONPNT3, ROTORD, ROTORG
   - CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8, PPLANE, MATHE, MATG
   - POINT, PBRSECT, PBMSECT, CBARAO
   - QVECT
   - DVGRID
   - TIC, TSTEP1
   - NSM, NSM1, NSML, NSML1, NSMADD
   - SEQGP

 - fixed:
    - update_card method:
      # update cards similar to Nastran's optimization
      # on grid 1234, update cp (field 2) to a value of 7)
    - model.update_card('GRID', 1234, 2, 7)
    - speed issue on CORD1x cards
    - get_loads now supports PLOAD4 cid/N123

  - other bugs:
    - fixed DIVERG bug (bad key when adding card)
    - fixed PFAST bug (bad xref)
    - really annoying DTABLE/TABLED1 bug (they were both called by DTable,
      but are very different)
    - fixed bug in SET3 card (multiple THRUs on single card)
    - write_path function works (Windows)
    - corrections to CBAR/CBEAM area formulas
    - CBAR/CBEAM MSCBLM0 is not MSCBLMO

 - minor changes
    - simplified comment writing per Jeff Lyon
      >>> card.comment = 'this is a comment'
      '$ this is a comment\n'
      >>> card.comment = 'this is a multi lined\n  - comment'
      card._comment =
      '$ this is a multi lined\n'
      '$  - comment\n'

 - API changes
   - TABLEDx are no longer in tables and instead tables_d
   - TABLEMx are no longer in tables and instead tables_m
   - card init methods now use nids instead of ga/gb so all cards are consistent
   - frequencies is now a dictionary of lists (not a single value)

 - bdf_renumber
   - removed execs from bdf_renumber, so Python 3 will work
   - fixed bug with handling Nones (caused by exec change)
   - fixed case control update bug

 - Aero:
   - handling of non-uniform lsrb/lrib on CAERO2
   - fixed unxref bug in AESURF
   - fixed xref bug in CSSCHD
   - fixed FLUTTER unxref bug
   - FLFACT can now write in wide field format


OP2:
 - added:
   - read_op2(include_results=None, exclude_results=None)
   - read_op2_geom(include_results=None, exclude_results=None)

   - MSC Nastran 2005R3B supported
   - Radioss supported
   - eigenvector_ROUGV1 now supported
   - random tables (e.g., PSD, ATO, RMS) for OUG-style tables
   - OUG1, OAG1 support

 - vectorized:
   - RealCBeamForceVUArray -> cbeam_force_vu
   - ComplexCBeamForceVUArray -> cbeam_force_vu

 - fixed:
   - MEFF, ASSIG, ASEPS tables now skipped properly
   - removed beta_transforms from gpforce
   - stress/strain/force transforms now support complex numbers

 - API changes:
   - real eigenvectors have eign; complex eigenvectors have eigi/eigr;
     thus a different variable and frequency formula is required
   -
GUI:
 - improved PyQt5 support
 - significantly faster reading

 - added:
   - marker/arrow (force) plots
   - animation menu

       Name   Subcase     Time    Scale Factor  Phase Angle
       -----  --------  --------  ------------  -----------
       Scale  Constant  Constant    Variable      Constant
       Freq.  Constant  Constant    Constant      Variable
       Time   Variable  Variable    Constant      Constant

   - NaN colors for invalid results (e.g., normal=NaN for a
     CTRIA3 that has 2 coincident points)
   - area picking for groups
   - custom displacement results
   - probe support by mouse picking
   - distance support by mouse picking
   - rotation center by mouse picking
   - Python 2 only: nicer Python console
   - support for PyQt5 (not quite perfect)

 - fixed:
   - display issue with CHEXA20 (element was ordered wrong)
   - element selection for quadratic solids no longer crashes
   - all legend values are now saved (e.g., the ones from cart3d)
   - endian issues

 - changes:
   - significant speedups
   - cycle_results now takes an integer instead of a string, so it works
      when you call it with more than 10 results  or have transient solutions


Known Issues
============
BDF:
 - CBEAMs wa/wb NSMs lead to incorrect mass properties
 - write_path function has problems on Linux (low priority; 1/10)

OP2:
 - SPC bug (component=7)
 - SOL 601 issues

GUI:
 - node Cd deflection issue for results (high priority, 4/10)
 - is this a solid element issue?
   - fix Min Interior angle (moderate priority; 1/10; see optistruct model - quad?)
   - fix Aspect Ratio (moderate priority; 1/10; see optistruct model - quad?)
 - issue with picking/focusing on a CAERO/CONM2 grid, which doesn't have a result
   it causes a segfault (high priority; 4/10)
 - fix bar Cd frame for cylindrical/spherical coords (low priority; 4/10)
   - see femap_exhaust
   - old : transform_displacements_to_global(i_transforms, beta_transforms)
   - new : transform_displacements_to_global(i_transforms, coords, xyz_cid0=None)
