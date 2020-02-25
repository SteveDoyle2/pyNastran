This document discusses the major features (and limitations) added in each
version of the software.  See the HTML documentation for more information.

If you like it/want to help leave some feedback on the
Dicussion Page:       https://groups.google.com/forum/#!forum/pynastran-discuss
Developer Dicussion:  https://groups.google.com/forum/#!forum/pynastran-dev

If you have a bug/want a new feature or card, leave some feedback on the
Issues List:     https://github.com/SteveDoyle2/pyNastran/issues

*** Not done....

--------------------------------Future--------------------------------
Overall Notes

Added BDF
   - *** SET2,PEALST,DMI,CONM1

Overall OP2
   - ***added support for SORT2

BDF Bug Fixes
BDF API Additions

Removed BDF

Overall OP2
   - added plotting support for "standard" results (***complex)
      - ***complex results may plot in real/imag or mag/phase or PSD
        (power spectral density)
      - ***no temperature, stress, element forces, energy


Overall F06
F06 Reader
F06 Writer

Solver
   - not released
   - rewriting pyNastran solver (SOL 101)
      - improved support for 3D
      - uses case control deck to identify loads/BCs
      - still uses dense matricies, but layout will be easier to add sparse
        matricies later
      - still need to build Kgg/Fg/SPC/MPC matrices, partition matrices,
        solve, calculate responses
   - added D matricies for PSHELL, PCOMP, PSOLID

Solver (New):
   - not released
   - rewriting pyNastran solver (SOL 101)
      - uses case control deck to identify loads/BCs
      - F06/OP2 output
      - SOL 101 - Static Displacement
         - outputs:
           - OLOAD Resultant
           - Grid Point Weight Generator (calculated; no table)
           - Displacement
           - Applied Load
           - SPC Forces
           - Stress
           - Strain
           - Element Force
           - MPC Forces (TODO)
           - Strain Energy (TODO)
         - cards:
           - GRID
             - CD support - verified for CORDxR
             - PS support - not verified
             - SEID support
           - CORDx
           - CELAS1/CELAS2
             - no SPOINTs
             - stress/strain/force output - verified
             - stiffness matrix - verified
           - CELAS3/CELAS4
             - no GRIDs
             - no stress/strain/force output
             - stiffness matrix - verified
           - CONM1/CONM2
             - no cid
             - GRID only
             - mass matrix - not verified
           - CROD/CONROD/CTUBE
             - mass matrix - verified
             - stiffness matrix - verified
             - output - verified
           - constraints:
             - SPC, SPC1
           - loads:
             - LOAD
             - FORCE/MOMENT
               - GRID only
             - SLOAD
               - SPOINT only
             - SPCD
           - CBAR
             - PBAR/PBARL
             - stiffness matrix - not verified
             - no mass matrix
             - no output
             - not verified
           - CBEAM
             - PBEAM (no PBEAML)
             - no warping
             - stiffness matrix - not verified
             - no mass matrix
             - no output
           - PSHELL/PCOMP
             - CTRIA3
               - no stiffness matrix
               - mass matrix - not verified
             - CQUAD4
               - mass matrix - not verified
               - stiffness matrix - not verified
               - no PSHELL mid3
               - no output

-----------------------v1.3   released 2020/?/??-----------------------
BDF:
 - bug fixes:
 - new cards:
   - DTI besides UNITS are now supported
   - PBEAM3, TEMPB3
 - enhancements:
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
 - minor enhancements:
   - handling ;; in DEQATN
   - support for SURF cid=None in _get_forces_moments_array
   - added check for number of loads = nscale factors for LOAD
   - more PARAM validation
- changes:
   - warning on RBAR dof check instead of RuntimeError

OP2:
 - enhancements:
   - more hdf5 results
   - a few more cards supported in the OP2 geometry reader
   - adding op2 writer (includes geometry support)
   - adding composite failure indicies (OEFIT)
   - improved sensitivity support (DSCMCOL)
   - preliminary 64-bit support

GUI:
 - bug fixes:
   - better argument handling
 - enhancements:
    - Can now load custom fringe/displacement/force results with an incorrect
      number of nodes/elements.  It's assumed that the node/element id in the
      first column should be the same as the id for the model.  Thus, you just
      filter out extra nodes/elements or you set some blank nodes/elements.
       - For integer results, no masking is applied.
       - For float results, masking is applied and nan corresponds to no color.
    - added mark_elements, mark_elements_by_case

-----------------------v1.2.2 released 2020/1/x-----------------------
BDF:
 - bug fixes:
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
     - fixing bug in set_param
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
  - added check on highlight menu for model existance
  - fixed import for new version download menu
  - improved command line error message
  - fixed support for CAEROx models without elements
  - export_cases now supports integers

-----------------------v1.2.1 released 2019/5/24-----------------------
OP2:
 - fixed bug with OUGV1PAT table

-----------------------v1.2.0 released 2019/5/21-----------------------
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

-----------------------v1.1.0 released 2018/6/26-----------------------
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
---------------
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
------------------
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


-----------------------v1.0.0 released 2017/5/25-----------------------
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
---------------

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
------------------
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
------------
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

-----------------------v0.8.0 released 2016/8/21-----------------------
Major changes:
 - BDF:
   - supports unicode
   - 278 total cards
   - lots of bug fixes
   - much better SOL 200/optimization support
 - OP2:
   - 500x faster
     - removed v0.6 style dictionary result objects
     - it's just too slow to maintain
   - much better NX Nastran support
   - much better SOL 200/optimization support
   - pandas objects can be made for iPython/Jupyter notebook
 - much improved testing
 - much improved GUI

BDF:
 - unicode is now supported
 - support for the following in the header:
   $ pyNastran: encoding=latin1
   $ pyNastran: punch=True
   $ pyNastran: version=nx
   $ pyNastran: skip_cards=PBEAM,CBEAM
 - self.set_as_msc() / self.set_as_nx()
 - added bdf.get_cards_by_card_types to get list of cards by card type
  - self._type_to_id_map
  - self._type_to_slot_map
  - self._slot_to_type_map
 - added:
   - cards
     - splines: SPLINE3
     - aecomps: AECOMP
     - pbushts: PBUSHT
     - pelasts: PELAST
     - pdampts: PDAMPT
     - epoints: EPOINT
     - transfer_functions: TF
     - loads: QVOL, GMLOAD
     - csschd: CSSCHD
     - usets : USET, USET1
     - coords: GMCORD
     - spcs: GMSPC
     - tempds : TEMPD
     - se_sets: SESUP
     - se_bsets: SEBSET/SEBSET1
     - se_csets: SECSET/SECSET1
     - se_qsets: SEQSET/SEQSET1
     - dtable : DTABLE
     - dequations : DEQATN
     - dphases : DPHASE

   - more xref
     - attempting to xref every card
   - supports reading executive/case control decks
   - renumbering
   - parsing of "SET 1 THRU 10 EXCEPT, 8, 11, 12" in case control deck

 - fixes:
   - PBEAM has been reworked
   - SUPORT1 dropping extra cards with same sid
   - fixed bug in SPCAX with init
   - fixed bug in CBUSH with GB of None
   - QHBDY handles missing nodes better
   - CFAST doesn't crash for blank nodes
   - PFAST handles cid=-1
   - renamed nlparm.id to nlparm.nlparm_id as it was obviously wrong
   - renamed nlpci.nlparm_id to nlpci.nlpci_id as it was obviously wrong
   - self.setSuper is now self.se_sets
   - RBE3 is ambiguous; no longer handling bad cards
   - model.material_ids
   - MPC writing bug
   - add_card and print_card_16 now handles float32
   - CORD1x referencing bug
   - many DMIx fixes

 - API changes:
   - cards are now much easier to call directly
     e.g. instead of passing in a BDFCard(), you can do:
         grid = GRID(nid, cp, x, y, z, cd, ps, seid)
     Not all cards have defaults fully defined.  The order of cards may change
     such that all optional fields are at the end in field number order. This
     is card dependent. Old approach can be done using
     `model.add_card(...)` or `GRID.add_card(card, comment)`.
   - case control now returns PRESSURE as DISPLACEMENT and not VECTOR
   - SPCD is now in self.loads, not self.spcs
   - TABDMP1 is now in self.tables_sdamping
   - SUPORT1 is now in self.suport1, not self.spcs
   - renamed spoints.spoints to spoints.points
   - deprecated echo_bdf
   - renamed bdf.aesurfs to bdf.aesurf

OP2:
 - automatically identifying NX vs MSC Nastran if header table found
 - write_f06 for all classes
 - using _eigenvalue header across all transient results
 - many more f06 table writing
 - element_names removed from objects (element_name is used instead)
 - SPOINTs in tables (e.g. displacement/spc forces) now write properly
 - added subtable_name to objects
 - improved catching of SORT1 vs. SORT2
 - improved catching of REAL vs. IMAG
 - better SOL200 support
 - read_op2 is ~500x faster
 - removed non-vectorized OP2
 - pandas is supported

OP2_Geom:
 - Geometry can be accessed directly from the OP2 without needing to
   read the BDF
 - the same BDF objects are still created
 - not all cards are supported
 - not all cards are even written by Nastran
 - using the BDF/OP2 separately is recommended, but when it works,
   it's much faster
 - doesn't support card limiting

OP2 Issues:
 - not all NX tables supported; some still cause crashes
 - PLOTEL support is likely still buggy
 - XYPLOT can cause crashes

GUI:
 - nodal/elemental results at the same time
 - added displacement deflection results
 - improved legend control
 - added sub-grids (e.g. SPCs) & sub-grid control menu
 - added label support with 'p'
 - csv nodal/centroidal result loading for every format
 - multiple OP2 files
 - added scripting dock
 - preliminary groups

-----------------------v0.7.2 released 2015/04/25-----------------------
Relase notes:
   - This is a bug fix release
   - improved documentation

BDF:
   - fixed bug with write_card on GRID card
   - fixing bug with bdf.get_material_id_to_property_ids_map
   - fixing deprecation messages; missed a few stack traces

BDF API:
   - these are internal methods
     - build_thru(packs, maxDV=None) is now build_thru(packs, max_dv=None)
     - build_thru_packs(packs, maxDV=None) is now build_thru_packs(packs, max_dv=None)
     - build_thru_float(packs, maxDV=None) is now build_thru_float(packs, max_dv=None)
     - build_table_lines(fields, nStart=1, nEnd=0) is now build_table_lines(fields, nstart=1, nend=0)
   - added get_points(...) function for CAERO5
   - subcases now support:
      ``if 'LOAD in subcase:``
      ``data = subcase['LOAD']``

OP2:
   - fixed bug with gridPointForces not being saved
   - fixed issue with test_op2 suppressing exceptions
   - fixed bug with ComplexSolidStrainArray being a stress result

OP2 API:
   - op2._found_results and op2._saved_results have been replaced by op2._results
     to prevent invalid results from being requested

GUI:
   - fixed nastran GUI loading issue (BDFs didn't show up in popup)

-----------------------v0.7.1 released 2015/04/17-----------------------
Relase notes:
   - removed deprecation crashes

-----------------------v0.7.0 released 2015/04/17-----------------------
Release notes:
   - much improved Python 3 support
   - added large field BDF writing
   - more OP2/F06 results (e.g. mass properties)
   - OP2 reader is faster, more robust, and catches FatalErrors
   - improved GUI to supports Cart3d, Panair, and Usm3D and log messages
   - making bdf/f06/op2 APIs more consistent

BDF:
   - added MAT11, ACCEL1, CAERO3/PAERO3, DMI, TABLED4/TABDMP1,
     BCTPARA, BCRPARA, BCTADD, BCTSET, BSURF, BSURFS
   - MATS1, MATTx, TABLED4, TABDMP1, CAERO3, PAERO3 supported
   - fixed major PBEAM bug
   - fixed bugs in PBUSHT, PDAMPT, PELST, TABRND1, PTUBE, EIGC
   - imporved SET/SPOINT writing
   - improved messages when cross-referencing
   - readBDF is now gone, replaced by read_bdf (was previously deprecated)
   - large field & double precision format writing now supported
     write_bdf method now has a is_double=True/False option
   - fixed bug with SETx writing.  Nastran limits a SETx card to one THRU
     per card, which is inconsistent with other cards that use THRU.
   - fixed GRID Position bug when using referenced cylindrical/spherical
     coordinate systems
   - fixed CONM2 offset bug for coordinate systems
   - added CAERO1/2 sub-element access methods for use in GUI
   - additional xref options
   - improved error handling

OP2:
   - rewrite of the OP2 allows for increased robustness and speed
   - Most fatal errors caught
     - BDF input errors not caught
   - added grid point weight table (see beam_modes demo)
   - OP2 no longer defines the filename; use read_op2.
     readOP2 is now gone, replaced by read_op2 (was previously deprecated)
   - renamed OP2 class names to be more descriptive
     (e.g. DisplacementObject -> RealDisplacement)
   - OP2 vectorized support.  It's slower, but it takes up much less
     memory (~5x).  Note: the F06 writing cannot be vectorized, so if
     you do vectorize your code, it may be much faster.  The API and data
     members are still in development and should be expected to change.
   - superelement support

F06:
   - removed; use the OP2

GUI:
 - added much improved qt GUI
   - results sidebar
   - controllable legend
   - vectorized support for reduced memory usage
   - color coded logging
   - scripting capability
 - removed wx GUI

Converters / Additional GUI Options
 - pyNastran's code base makes it easy to develop other useful tools
   that make use of common code.  As such, additional formats are supported
   in terms of readers/writers/converters/viewing, but are not a main focus.
   These include:
   - Cart3d
   - Panair
   - STL
   - Usm3D

Known issues:
=============
  - 1. Unicode
  - 2. doesn't support elements/properties with the same ID
  - 3. card replication
  - 4. Executive Control Deck INCLUDE files
    5. DMIG headers not defined at the top of the matrix
  - 6. cannot parse doubly included BDFs if the include file
      is at the end of file A

      Fails
      =========
      A.bdf
      ---------
      CEND
      BEGIN BULK
      GRID,2,,1.0, 0.0, 0.0
      INCLUDE B.bdf

      B.bdf
      ---------
      INCLUDE C.bdf

      C.bdf
      ---------
      GRID,1,,0.0, 0.0, 0.0

     This will work
     ==============
      A.bdf
      ---------
      CEND
      BEGIN BULK
      INCLUDE B.bdf
      GRID,2,,1.0, 0.0, 0.0

      B.bdf
      ---------
      INCLUDE C.bdf

      C.bdf
      ---------
      GRID,1,,0.0, 0.0, 0.0


-----------------------v0.6.1 released 2013/06/13-----------------------
BDF Bug Fixes
   - Parsing of "1.2345D+04".  It worked in Python 2.6, but not 2.7.
   - RBE2 card with alpha field.
   - active I12 field on PBAR card now disables AssertionError if K1 or K2 are
     blank
   - fixed bug in LSEQ card where table field was cross-referenced as a load
    (field #4)

OP2/F06 Bug Fixes
   - improved memory usage in real plates (stress/strain) when writing F06

-----------------------v0.6 released 2013/05/07-----------------------

Overall Notes
   - significant improvement in BDF reader
   - BDF reader now supports comments.  It's stored on the card._comment
     variable, but not all cards have one, so use the card.comment() method.
   - BDF reader verifies types on the card fields (e.g. float is not allowed
     on integer/blank field)
   - integrated readBDF_Punch into readBDF method with a punch=True flag
   - minor bug fixes to OP2 reader
   - minor bug fixes to F06 reader
   - more PEP-8

Overall BDF
   - significant improvement in BDF reader
   - BDF reader now supports comments.  It's stored on the card._comment
     variable, but not all cards have one, so use the card.comment() method.
   - BDF reader verifies types on the card fields (e.g. float is not allowed
     on integer/blank field)
   - integrated readBDF_Punch into readBDF method with a punch=True flag

Added BDF
   - NLPCI, PRAC2D, PRAC3D, PLPLANE
   - added _verify methods to many cards

-----------------------v0.5 released 2012/08/13-----------------------
Overall Notes
   - lots of new BDF cards
   - OP2 now supports element real/imaginary forces
   - OP2 is much, much better at reading frequency analysis results
   - new OP4 reader in pure Python.  Writer isn't fully functional yet.
   - improved code to better support Python 3
   - implementing PEP-8

BDF API Changes
   - moved files around; if you want your imports to work, import from bdf.py
   - moved storage location of EIGC, EIGP to cMethod
   - grav card are now accessed using Load(lid) and not Grav(lid)
   - grav card loadID renamed from sid to lid
   - LOAD cards now have a getReducedLoads() which gets the scale factors and
     load objects for all sub-loads (supports LOAD, FORCEx, MOMENTx, GRAV)
     not implemented across all load types, so having a LOAD card is still a
     good thing
   - LOAD objects now have a getReducedLoads() method that determines
     a list of scaleFactors and FORCE/MOMENT/GRAV/etc. cards

BDF Bug Fixes
   - SPOINTs were lost when cross referenced
   - EIGC, EIGP now stored in self.cMethods instead of self.methods
   - TABRNDG is now in randomTables
   - CREEP field g now reads/writes out
   - PLOAD4 bug with g1 vs g3/g4
   - EIGB keeps method (INV/SINV) now (not clear if there's a default)
   - QBDY3 cross referencing fixed
   - TSTEPNL fixed possible issue (kStep default) caused by unlisted METHOD
     field
   - SLOAD nodes were referencing loads instead of nodes
   - PFAST kr3 was writing out as kr2
   - PSHELL now supports mid2=-1
   - fixed PVISC, PGAP cross referencing (these dont require it)
   - CELASx and CDAMPx cards now allow empty fields / grounding
   - GRAV cards don't crash if the max N value == min N value

Added BDF
   - ASET, ASET1, BSET, BSET1, CSET, CSET1, QSET, QSET1
   - DMIG, DMIJ, DMIJI, DMIK
   - TABRNDG
   - TLOAD1, ACCEL1,RFORCE,PLOADX1
   - SPLINE4, SPLINE5
   - CBUSH1D, PBUSH1D
   - PRAC2D, PRAC3D

Removed BDF
   - CONM1

Overall OP2
   - moved files around; if you want your imports to work, import from op2.py
   - test_op2 now crashes if there isn't a .op2 in the filename
   - major changes to reader to support non-standard solutions
     reader should be much better at getting correct values
   - added majority of element forces
   - added majority of element heat fluxes
   - significantly more support for complex results
   - added plotting support for "standard" results (real)
      - displacement, velocity, acceleration,
        eigenvectors, spcForces, mpcForces,
        gridPointForces, appliedLoads, loadVectors, forceVectors

OP2 API Changes
   - renamed many of the objects to follow pep-8, shouldn't really affect
     much code
   - complex tables now have same format as the real tables
     all results are stored in real/imaginary format and not
     magnitude/phase, regardless of what was defined in the OP2

Overall F06 Reader/Writer
   - updated to work with OP2 updates

Added F06 Writer
   - many, but not all loads

Overall GUI
   - added Cart3d support
   - added Panair support

Added GUI
   - CQUADR, CTRIAR

GUI Bug Fixes


-----------------------v0.4.1 released Py2.x 2012/06/05-----------------------
Overall
   - Python 3.x works!

BDF Bug Fixes
   - coordinte system transpose bug (messed up most coordinate systems)
     affected CORD1R, CORD2R, CORD1C, CORD2C, CORD1S, CORD2S
   - Fluid Structure Interaction coordinate system of -1 is now allowed
   - nastran double precision format (e.g. 1.23456D+03) now supported
   - better rounding of floats when writing BDF
   - SET3 printing

Added BDF
   - coordinate systems in CONM2 are fully supported
   - GRAV card prints out nicer

Added OP2
   - GEOM1S, GEOM2S, GEOM3S, GEOM4S tables (for Radioss)

-----------------------v0.4 released 2012/04/23-----------------------
Overall Notes
   - new ASCII/binary OP4 reader for dense matrices
   - lots of BDF bug fixes
   - lots of new BDF cards
   - added support for OpenMDAO user-defined parameterization syntax %varname
   - new results in OP2
   - new F06 Reader/Writer
   - F06 Writer integrates with OP2 Reader

Overall OP4
   - ASCII/binary reader
   - sparse matrices not supported

Overall BDF
   - 180 total cards
   - added tab support (no mixed tabs, spaces, commas)
   - reworked format of bdf_crossReferencing.txt to be more like
     bdf_readWrite.txt
   - added test code that compares cards values from one run to the next
   - added support for THRU-BY on certain cards (e.g. QBDY3)
   - added new custom Exceptions to make errors messages clearer
   - added support for user defined optimization values (using %varName)
     values set with "bdf.setDynamicSyntax(dictOfVars)"

BDF API Changes
   - removeCards(['GRID','CQUAD4']) is now called disableCards(cards)
   - removed obsolete isPrintable and Is methods from BDF class
   - added Mid1(), Mid2(), Mid3(), Mid4() methods for PSHELL.
     Mid() now returns Mid1() or Mid2() depending on the form of the PSHELL
     Mass() calculated based on mid1 (or mid2 if mid1=None) per QRG
   - Tri-type elements are now TriShells instead of CTRIA3s (same for QUADs,
     QuadShells)
   - component constraints (e.g. 123456) on GRID, SPC, SPC1 and MPC are now
     strings
     for easy looping instead of a single integer (default='')
   - added Mids() for PCOMP/PCOMPG
   - added getNodeIDToElementIDsMap(), getPropertyIDToElementIDsMap(),
     getMaterialIDToPropertyIDsMap()
     limited 0D,1D element support (e.g. CBAR, CBEAM, CELAS1, etc.) for these
     methods

BDF Bug Fixes
   - printCard can now handle cards with an embedded line of Nones
       with data below.  This is a rare bug but affects the NLPARM
   - fixed minor bug with small values not printing that were slightly
     above the tolerance level
   - fixed bug in CONM2 with inertia terms
   - PCOMP can now have 1 ply per line or 2 plies per line
   - PCOMP doesnt lose last plies anymore
   - fixed bugs in CQUADR, CQUADX, CVISC,  CONROD, CGAP,  CTRIAX, PTUBE
   - fixed bugs in NLPARM, PLOAD2, PLOAD4, QBDY1,  QBDY3, FREQ2,  CHBDYE
   - fixed bugs in PMASS,  CSHEAR, CBAR,   CQUAD,  RBE1,  RBE2,   RBE3
   - fixed bugs in PBEAML, MAT1, MAT3, CTRIAX6, CREEP
   - CREEP not lost when writing BDF out
   - disabled PBEAML

Added BDF
   - enabled CQUAD, RBE3
   - EIGB, EIGC, EIGP, EIGR, EIGRL
   - CBUSH, PBUSH, PDAMPT, CFAST, PFAST, CONM1, CVISC, PVISC
   - TSTEP, TSTEPNL
   - DLOAD, SLOAD, RLOAD1, RLOAD2
   - CAERO2, PAERO2, AEFACT
   - CBEND, PBEND, CBEAM3
   - MATS1, MATHP
   - DVPREL2, DVMREL1, DOPTPRM, DLINK
   - DRESP2
   - SPOINT

Overall OP2
   - more results
   - integrates with F06 Writer
   - more embedded BDF reading
   - added setTransientTimes(times) to explictly extract results at desired
     times
     - OUG
     - OES
     - OEE
     - OQG

OP2 API Changes
   - more results in CBEAM stress/strain
   - complex displacement reworked
   - added writeF06 method to most objects (no strain energy,
     shear stress/strain, celas stress/strain)

Added OP2
   - LAMA (requires an ESE = ALL to create table)
     - eigenvalues (real)

   - OUG (DISPLACMENTS=ALL)
     - improved complex displacements
     - ***complex velocity
     - ***complex acceleration
   - OQG (SPCFORCES=ALL, MPCFORCES=ALL)
     - SPC Forces (static/transient)
     - MPC Forces (static/transient)

   - OES (STRESS=ALL, STRAIN=ALL)
     - stress/strain
       - CTRIAX6 (static/transient)
     - hyperelastic stress/strain (transient)
       - CQUAD4
     - nonlinear strain (transient)
       - CTRIA3, CQUAD4
       - CROD, CONROD, CTUBE
       - ***CBEAM
   - OGF (GPFORCES=ALL)
     - grid point forces (static/transient)
     - load vector (static/transient)

Overall F06 Reader/Writer
   - rewrote f06 reader
      - format1, real, SORT1 supported
      - limited transient support
      - F06 Reader not nearly as complete as the OP2 Reader
      - F06 field parsing needs more testing, so beware!
      - juat use the OP2 reader
   - F06 Writer integrates with OP2 reader
     - transient variable names may not be correct (e.g. 'lftsfq' mean
       'LOAD STEP', but it should write 'LOAD STEP')
   - The op2 or f06 readers feed the f06 writer or another solver/user's own
     results.
   - see f06_readWrite.txt for more info

Overall GUI
   - CSHEAR, CTRIAX6 support
   - "python setup.py install" should work now, but it's still not recommended

-----------------------v0.3 released 2012/02/08-----------------------
Overall Notes
   - coordinate systems!!!
   - all user BDF card API methods will have the first letter capitalized
     (not fully implemented)
   - added OP2 demo program that calculates solidStress & displacement margins
   - addd preliminary GUI (pyNastranGUI) that can view results

BDF API Changes
   - added simpler method to remove cards from the BDF
   - standardizing methods for element/property/material cards
     - Eid, Pid, Mid, Mass, NodeIDs, Centroid
     - Thickness, Area, Volume, Length
     - J, I11, I12, I22
   - FREQx  methods
   - CAERO1 methods
   - SETx   methods

BDF Bug Fixes
   - enabled PDAMP, CREEP
   - bug fix in NLPARM, MAT9, GRDSET
   - BEGIN BULK may now have multiple spaces in it
   - Continuation markers can have more than 1 character in them for CSV cards
   - Unincluded cards at the end of an INCLUDE file won't crash the BDF reader
   - CORD2Rs can have axes that arent ijk!
   - fixed bug in fieldWriter that caused certain values (e.g. -0.9999 and
     -0.0999) to be written out as (-.)  They were assumed to not round up to
     -1.0 and -0.1 respectively so the leading -0. could be replaced by (-.)
     to save field width.

Overall BDF
   - reorganized code files (main BDF class is in the same place)
   - improved accuracy of PBEAM methods (Area, Nsm, etc.) by integrating
     instead of using the first value
   - added cross referencing of SPLINE1
   - test_bdf only prints a reject message once per card type now
   - more Volume, Centroid methods

Added BDF
   - PBEAML
   - TABLED1,TABLED2,TABLED3,
     TABLEM1,TABLEM2,TABLEM3,TABLEM4,
     TABLES1,TABLEST,TABRND1
   - PAERO1
   - FREQ,FREQ1,FREQ2
   - CGAP, PGAP
   - PDAMP5
   - PSHEAR
   - AEPARAM,AELINK,AESTAT
   - SET1,SET3,SESET
   - CORD1R,CORD1C,CORD2C,CORD1S,CORD2S
   - CONM2 (cid!=0)
   - SPLINE2
   - TRIM

Added OP2
   - GEOM2
     - CGAP
   - OES
     - CSHEAR
     - CTRIAR, CTRIA6
   - OUG
     - velocity
     - accleration

OP2 API Changes
   - OUG
      - displacementObject now has self.translations (no more
        self.displacements) & self.rotations
        This is standard across displacement/velocity/acceleration

Overall OP2
   - reorganized code files (main OP2 class is in the same place)

OP2 Bug Fixes
   - deviceCode is now passed properly to result objects
   - nodeIDs now are correct for transient results

-----------------------v0.2.1 released 2011/1/2-----------------------
BDF Bug Fixes
   - fixed bug in fieldWriter that caused certain values (e.g. -0.9999 and
     -0.0999) to be written out as (-.)  They were assumed to not round up to
     -1.0 and -0.1 respectively so the leading -0. could be replaced by (-.)
     to save field width.

-----------------------v0.2 released 2011/12/20-----------------------

BDF API Changes
   - test_bdf can be called from the command line
   - renamed read method to readBDF (it got confusing when the op2 object can
     read the bdf b/c it needs to initialize the bdf cards for the geomX tables)
   - renamed write to writeBDF
   - renamed writeAsPatran to writeBDFAsPatran
   - thermal materials (MAT4, MAT5) pulled from self.materials and put into
     self.thermalMaterials new ThermalMaterial & StructuralMaterial methods
     Material method checks for either still
   - rigid elements (RBAR, RBE2, etc.) separated out of elements to avoid
     overwriting cards new RigidElement method
     Element method doesnt check for rigid elements anymore

Overall BDF
   - a few more cards
   - added extra argument to all cards to set data directly from the
     OP2 (no data type setting as the OP2 is binary).
   - max lines in executive & case control decks is now 600
   - methods for Elements (Area, Nsm, Thicknesss)


Overall OP2
   - test_op2 can be called from the command line
   - reading support ONLY
   - static solutions (e.g. SOL 101) results will be very stable
   - thermal solutions have significantly less support than
     structural solutions.
   - SORT2 is not supported...for example:
         use   DISPLACEMENT(PLOT,SORT1)=ALL
         not   DISPLACEMENT(PLOT,SORT2)=ALL
   - Many examples outside of SOL 101 have been tested and don't crash.
     This includes large SOL 200 models.  Verify your results first.
   - The user has the option to read the Geometry tables (GEOM1 (nodes/coords),
     GEOM2 (elements), GEOM3 (loads), GEOM4 (constraints),
     EPT/EPTS (properties), MPT/MPTS (materials) and can write the data to a
     BDF.  This is still a preliminary option.  Use the makeGeom argument.
     You still need to write the BDF.

BDF Bug Fixes
   - DCONSTR, DDVAL, DAREA weren't processed properly
   - fixed RBE1 bug
   - fixed RBE3 bug & disabled RBE3
   - AERO/AEROS cards split from self.aeros object to avoid overwriting entries
   - FLFACT reading when the card is short
   - PCOMP Thickness method doesnt crash when you want the total thickness
     Better defaults on PCOMP
   - fixed some fields on disabled CORDxx cards
   - CELASx cards had issues with proper order of entries
   - PELAS cards werent processed properly
   - files on remote drives dont cause crashes anymore
   - cend now also signifies the end of the executive control deck
   - begin BULK is valid now
   - fixed issue with the case control deck parsing on some uncommon cards

Added BDF
   - PBARL, PDAMP
   - CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5
   - AESTAT
   - DRESP1, DVPREL1, DRESP1

OP2 Tables
   - GEOM1/GEOM2/GEOM3/EPT/MPT (basically stores the bdf)
     Cards that are not stable are not supported.  The
     writeAsPatran method may be used, but there may still be a
     few bugs in some of the cards.  These tables are very stable.

   - OUG  table (displacement/temperature,eigenvalues
                 velocity,acceleration) is very robust
       - just verify you're getting the result you want, there are
         a lot of results in this table
   - OQG  table (spc/mpc forces) None
   - OGP  table (grid point forces) works for SOL 101, untested on others
     (shouldnt crash)
   - OSTR table (strain) - see OES table
   - OES  table (stress) - elements that arent supported are skipped
       - Be careful about how you request the data from the object as
         the variables are dynamic based on the Case Control Deck parameters
         selected.  For example:
             STRESS(PLOT,SORT1,FIBER,VONMISES)=ALL
                - composite stress oxx is in the fiber direction (not the
                  x axis)
                - isotropic stress elements have von mises stress (ovm)
             STRESS(PLOT,SORT1,FIBER,MAXSHEAR)=ALL
                - composite stress oxx is in the fiber direction (not the
                  x axis)
                - isotropic stress elements have max shear stress (maxShear)
       - All this information is contained within boolean variables such as
         isVonMises(), isMaxShear() (stress/strain),isFiber(), isCurvature()
         (stress/strain) which are located on the result objects.
       - Elements can have different outputs depending on the
         solution and user requested output options.  For example, rods can
         contain the margin of safety in tension and compression, but they
         don't always have that data.  You can check this by using the
         boolean isMargins.

   - OEE table (element energy) - stable, but results are questionable
         and difficult to verify.


-----------------------v0.1.1 released 2011/12/06 -----------------------

BDF Bug Fixes
- bug fix to elements.py to fix import error

-----------------------v0.1   released 2011/11/09 -----------------------

Overall BDF
   - an executive control deck, case control deck, and bulk data deck
     are REQUIRED.
   - A ENDDATA card is strongly recommened!
   - Mass/Stiffness Matricies not supported except in a few, rare cases.
   - superelements are NOT supported
   - loads need some TLC in v0.2, but there is minor support for them
   - see the HTML documenation for a list of object data and methods.
     The Nastran QRG (Quick Reference Guide) is another good reference.

Executive Control Deck
   - Include files NOT supported
   - Limited to LESS than 200 lines
   - Extract soltion and method (if sol=600)
   - Can update solution/method

Case Control Deck
   - Include files supported (dont cross the BEGIN BULK entry)
   - Limited to LESS than 200 lines
   - In general, if a card is supported, EVERY field can be read,
     accessed, and written.
   - Most cards forms supported
      1.  STRESS(PLOT,PUNCH) = ALL     #  name(options) = value  (STRESS-type)
      2.  STRESS = ALL                 #  name = value           (STRESS-type)
      3.  DISP(PLOT,PUNCH) = ALL       #  name(options) = value  (STRESS-type)
      4.  PARAM,FIXEDB,-1              #  name, options, value   (PARAM-type)
      5.  SET 1 = 1,2,3,4, etc.        #  name id = options      (SET-type)
   - Some forms not supported
      1.  AUTOSPC (PRINT, PUNCH, SID=100, EPS=1.E-6, MPC)=YES    (AUTOSPC=YES is supported)
      2.  MODALSE(ESORT=ASCEND,THRESH=.0001)= 100
      3.  MODESELECT (FLUID LMODES = 5)
      4.  OUTPUT(XYOUT)
      5.  TFL = 1, 25, 77
      6.  SURFACE 10 SET 9 NORMAL X3
      7.  VOLUME 21 SET 2
   - Methods to update global/local subcase values


BDF Data Deck (no cross referencing required)
   - include files supported
   - 8/16 character fields, csv format for cards is supported.
     DONT USE TABS!!!
   - Cards are written in 8-character fields and make use of Nastran-scientific
     notation to maximize precision.
   - write and writeAsPatran methods to write the bdf file after modifying it.
     The write method writes each 'section' (e.g. nodes/elements/properties)
     as isolated sections.  writeAsPatran intersperses properties and elements.
   - Duplicate node/element/property/etc. IDs are NOT allowed!

Bulk Data Deck Cross-Referencing
   - Entire deck is read.  If a card isnt supported, it is rejected.  Unparsed
     cards are clearly marked in the $REJECTS or $REJECT_LINES sections.
   - Cross-referencing is a feature that links data from one card to another.
     An example is a CQUAD4 has four nodes.  When cross-referenced, the element
     object's nodes become node objects.
   - Coordinate systems are a work in progress and difficult to test without
     good test cases.  Feel free to add some.  Therefore, the global coordinate
     system is safest and is the only one supported at the moment (CORD2R).
   - In order to determine element information that requires information from
     other cards, cross-referencing is REQUIRED.  The read method has the
     option to disable this.
   - When using constraints, all constraints (SPC,SPC1,SPCD,MPC) must be
     part of an SPCADD or MPCADD, unless there are no SPCADDs or MPCADDs.
     In other words, DONT mix and match!
   - Structural elements have a Mass/Area/Length/Thickness depending on type.
   - Volume ONLY supported for TET4 and TET10
   - Lots of methods for PCOMPs
   - CBEAM elements do are not complete for their area/mass/etc calcuations.
   - Moments of intertias of elements are not calculated.
   - Stiffness Matricies are incomplete
