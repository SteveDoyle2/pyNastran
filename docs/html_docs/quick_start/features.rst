==============================
Features
==============================

Overview
========
 - Python 2.7, 3.6-3.7
 - BSD-3 license
 - unicode support
 - importable from within Matlab

 - limited package requirements for BDF/OP2/F06

  - additional features available with more packages

    - BDF/OP2:

       - h5py for HDF5 input/output support
       - PyQt4/PyQt5/PySide/PySide2/wxpython for file loading popup

    - OP2:

      - pandas for results/matrices for use in the Jupyter Notebook

    - F06:

      - matplotlib support for plotting

    - GUI: range of choices

      - PyQt4/PyQt5/PySide/PySide2
      - VTK 7/8

    - logging using **cpylog**

      - colorama for console logging
      - HTML logging for Jupyter Notebook
      - no markup when piping output to a file
      - supports overwriting logger object with user-defined logger

BDF Reader/Writer
=================
 - Input/Output:

   - 473 cards supported including:

     - optimization
     - aero
     - thermal
     - superelements

   - small, large, double precision file reading/writing
   - pickling
   - HDF5 reading/writing
   - comments are stored
   - simplified card adding
     ```python
     >>> model.add_grid(nid, xyz=[4.,5.,6.], comment='nid, cp, x, y, z')
     ```
     ```
     $GRID comment
     $grid,nid,cp,x,y,z
     GRID,10,,4.0,5.0,6.0
     ```

 - methods:

   - loads summation
   - mass properties (including NSM)
   - nodal equivalencing
   - mesh quality

     - aspect ratio, taper ratio, skew, min/max interior angle
     - quad collapsing
     - element deletion
   - deck merging
   - renumber
   - unit conversion
   - cutting plane
   - visualization of material coordinate systems
   - mirroring
   - solid skinning
   - length, area, volume, mass breakdowns

- list of cards supported...

+------------------------+------------------------------------------------------------+
| Card Group             | Cards                                                      |
+------------------------+------------------------------------------------------------+
| MATS1                  | MATS1                                                      |
+------------------------+------------------------------------------------------------+
| MATT1                  | MATT1                                                      |
+------------------------+------------------------------------------------------------+
| MATT2                  | MATT2                                                      |
+------------------------+------------------------------------------------------------+
| MATT3                  | MATT3                                                      |
+------------------------+------------------------------------------------------------+
| MATT4                  | MATT4                                                      |
+------------------------+------------------------------------------------------------+
| MATT5                  | MATT5                                                      |
+------------------------+------------------------------------------------------------+
| MATT8                  | MATT8                                                      |
+------------------------+------------------------------------------------------------+
| MATT9                  | MATT9                                                      |
+------------------------+------------------------------------------------------------+
| aecomps                | AECOMP                                                     |
+------------------------+------------------------------------------------------------+
| aefacts                | AEFACT                                                     |
+------------------------+------------------------------------------------------------+
| aelinks                | AELINK                                                     |
+------------------------+------------------------------------------------------------+
| aelists                | AELIST                                                     |
+------------------------+------------------------------------------------------------+
| aeparams               | AEPARM                                                     |
+------------------------+------------------------------------------------------------+
| aero                   | AERO                                                       |
+------------------------+------------------------------------------------------------+
| aeros                  | AEROS                                                      |
+------------------------+------------------------------------------------------------+
| aestats                | AESTAT                                                     |
+------------------------+------------------------------------------------------------+
| aesurf                 | AESURF                                                     |
+------------------------+------------------------------------------------------------+
| aesurfs                | AESURFS                                                    |
+------------------------+------------------------------------------------------------+
| ao_element_flags       | CBARAO                                                     |
+------------------------+------------------------------------------------------------+
| asets                  | ASET, ASET1                                                |
+------------------------+------------------------------------------------------------+
| axic                   | AXIC                                                       |
+------------------------+------------------------------------------------------------+
| axif                   | AXIF                                                       |
+------------------------+------------------------------------------------------------+
| bconp                  | BCONP                                                      |
+------------------------+------------------------------------------------------------+
| bcrparas               | BCRPARA                                                    |
+------------------------+------------------------------------------------------------+
| bcs                    | CONV, CONVM, RADBC, RADM, TEMPBC                           |
+------------------------+------------------------------------------------------------+
| bctadds                | BCTADD                                                     |
+------------------------+------------------------------------------------------------+
| bctparas               | BCTPARA                                                    |
+------------------------+------------------------------------------------------------+
| bctsets                | BCTSET                                                     |
+------------------------+------------------------------------------------------------+
| blseg                  | BLSEG                                                      |
+------------------------+------------------------------------------------------------+
| bsets                  | BSET, BSET1                                                |
+------------------------+------------------------------------------------------------+
| bsurf                  | BSURF                                                      |
+------------------------+------------------------------------------------------------+
| bsurfs                 | BSURFS                                                     |
+------------------------+------------------------------------------------------------+
| cMethods               | EIGC, EIGP                                                 |
+------------------------+------------------------------------------------------------+
| caeros                 | CAERO1, CAERO2, CAERO3, CAERO4, CAERO5                     |
+------------------------+------------------------------------------------------------+
| convection_properties  | PCONV, PCONVM                                              |
+------------------------+------------------------------------------------------------+
| coords                 | CORD1C, CORD1R, CORD1S, CORD2C, CORD2R, CORD2S, GMCORD     |
+------------------------+------------------------------------------------------------+
| creep_materials        | CREEP                                                      |
+------------------------+------------------------------------------------------------+
| csets                  | CSET, CSET1                                                |
+------------------------+------------------------------------------------------------+
| csschds                | CSSCHD                                                     |
+------------------------+------------------------------------------------------------+
| csuper                 | CSUPER                                                     |
+------------------------+------------------------------------------------------------+
| csupext                | CSUPEXT                                                    |
+------------------------+------------------------------------------------------------+
| dareas                 | DAREA                                                      |
+------------------------+------------------------------------------------------------+
| dconstrs               | DCONADD, DCONSTR                                           |
+------------------------+------------------------------------------------------------+
| ddvals                 | DDVAL                                                      |
+------------------------+------------------------------------------------------------+
| delays                 | DELAY                                                      |
+------------------------+------------------------------------------------------------+
| dequations             | DEQATN                                                     |
+------------------------+------------------------------------------------------------+
| desvars                | DESVAR                                                     |
+------------------------+------------------------------------------------------------+
| divergs                | DIVERG                                                     |
+------------------------+------------------------------------------------------------+
| dlinks                 | DLINK                                                      |
+------------------------+------------------------------------------------------------+
| dload_entries          | ACSRCE, QVECT, RANDPS, RANDT1, RLOAD1, RLOAD2, TLOAD1,     |
+------------------------+------------------------------------------------------------+
|                        | TLOAD2                                                     |
+------------------------+------------------------------------------------------------+
| dloads                 | DLOAD                                                      |
+------------------------+------------------------------------------------------------+
| dmigs                  | DMIG                                                       |
+------------------------+------------------------------------------------------------+
| dmijis                 | DMIJI                                                      |
+------------------------+------------------------------------------------------------+
| dmijs                  | DMIJ                                                       |
+------------------------+------------------------------------------------------------+
| dmiks                  | DMIK                                                       |
+------------------------+------------------------------------------------------------+
| dmis                   | DMI                                                        |
+------------------------+------------------------------------------------------------+
| doptprm                | DOPTPRM                                                    |
+------------------------+------------------------------------------------------------+
| dphases                | DPHASE                                                     |
+------------------------+------------------------------------------------------------+
| dresps                 | DRESP1, DRESP2, DRESP3                                     |
+------------------------+------------------------------------------------------------+
| dscreen                | DSCREEN                                                    |
+------------------------+------------------------------------------------------------+
| dtable                 | DTABLE                                                     |
+------------------------+------------------------------------------------------------+
| dti                    | DTI                                                        |
+------------------------+------------------------------------------------------------+
| dvcrels                | DVCREL1, DVCREL2                                           |
+------------------------+------------------------------------------------------------+
| dvgrids                | DVGRID                                                     |
+------------------------+------------------------------------------------------------+
| dvmrels                | DVMREL1, DVMREL2                                           |
+------------------------+------------------------------------------------------------+
| dvprels                | DVPREL1, DVPREL2                                           |
+------------------------+------------------------------------------------------------+
| elements               | CBAR, CBEAM, CBEAM3, CBEND, CBUSH, CBUSH1D, CBUSH2D,       |
+------------------------+------------------------------------------------------------+
|                        | CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5, CELAS1, CELAS2,    |
+------------------------+------------------------------------------------------------+
|                        | CELAS3, CELAS4, CFAST, CGAP, CHBDYE, CHBDYG, CHBDYP,       |
+------------------------+------------------------------------------------------------+
|                        | CHEXA, CIHEX1, CIHEX2, CONROD, CPENTA, CPLSTN3, CPLSTN4,   |
+------------------------+------------------------------------------------------------+
|                        | CPLSTN6, CPLSTN8, CPYRAM, CQUAD, CQUAD4, CQUAD8, CQUADR,   |
+------------------------+------------------------------------------------------------+
|                        | CQUADX, CQUADX4, CQUADX8, CRAC2D, CRAC3D, CROD, CSHEAR,    |
+------------------------+------------------------------------------------------------+
|                        | CTETRA, CTRAX3, CTRAX6, CTRIA3, CTRIA6, CTRIAR, CTRIAX,    |
+------------------------+------------------------------------------------------------+
|                        | CTRIAX6, CTUBE, CVISC, GENEL                               |
+------------------------+------------------------------------------------------------+
| flfacts                | FLFACT                                                     |
+------------------------+------------------------------------------------------------+
| flutters               | FLUTTER                                                    |
+------------------------+------------------------------------------------------------+
| frequencies            | FREQ, FREQ1, FREQ2, FREQ3, FREQ4, FREQ5                    |
+------------------------+------------------------------------------------------------+
| grdset                 | GRDSET                                                     |
+------------------------+------------------------------------------------------------+
| gridb                  | GRIDB                                                      |
+------------------------+------------------------------------------------------------+
| gusts                  | GUST                                                       |
+------------------------+------------------------------------------------------------+
| hyperelastic_materials | MATHE, MATHP                                               |
+------------------------+------------------------------------------------------------+
| load_combinations      | LOAD, LSEQ                                                 |
+------------------------+------------------------------------------------------------+
| loads                  | ACCEL, ACCEL1, FORCE, FORCE1, FORCE2, GMLOAD, GRAV,        |
+------------------------+------------------------------------------------------------+
|                        | LOADCYN, MOMENT, MOMENT1, MOMENT2, PLOAD, PLOAD1, PLOAD2,  |
+------------------------+------------------------------------------------------------+
|                        | PLOAD4, PLOADX1, QBDY1, QBDY2, QBDY3, QHBDY, QVOL,         |
+------------------------+------------------------------------------------------------+
|                        | RFORCE, RFORCE1, SLOAD, SPCD, TEMP                         |
+------------------------+------------------------------------------------------------+
| masses                 | CMASS1, CMASS2, CMASS3, CMASS4, CONM1, CONM2               |
+------------------------+------------------------------------------------------------+
| materials              | MAT1, MAT10, MAT11, MAT2, MAT3, MAT3D, MAT8, MAT9, MATG    |
+------------------------+------------------------------------------------------------+
| methods                | EIGB, EIGR, EIGRL                                          |
+------------------------+------------------------------------------------------------+
| mkaeros                | MKAERO1, MKAERO2                                           |
+------------------------+------------------------------------------------------------+
| monitor_points         | MONPNT1, MONPNT2, MONPNT3                                  |
+------------------------+------------------------------------------------------------+
| mpcadds                | MPCADD                                                     |
+------------------------+------------------------------------------------------------+
| mpcs                   | MPC                                                        |
+------------------------+------------------------------------------------------------+
| nlparms                | NLPARM                                                     |
+------------------------+------------------------------------------------------------+
| nlpcis                 | NLPCI                                                      |
+------------------------+------------------------------------------------------------+
| nodes                  | EPOINT, GRID, SPOINT                                       |
+------------------------+------------------------------------------------------------+
| normals                | SNORM                                                      |
+------------------------+------------------------------------------------------------+
| nsmadds                | NSMADD                                                     |
+------------------------+------------------------------------------------------------+
| nsms                   | NSM, NSM1, NSML, NSML1                                     |
+------------------------+------------------------------------------------------------+
| nxstrats               | NXSTRAT                                                    |
+------------------------+------------------------------------------------------------+
| omits                  | OMIT1                                                      |
+------------------------+------------------------------------------------------------+
| paeros                 | PAERO1, PAERO2, PAERO3, PAERO4, PAERO5                     |
+------------------------+------------------------------------------------------------+
| params                 | PARAM                                                      |
+------------------------+------------------------------------------------------------+
| pbusht                 | PBUSHT                                                     |
+------------------------+------------------------------------------------------------+
| pdampt                 | PDAMPT                                                     |
+------------------------+------------------------------------------------------------+
| pelast                 | PELAST                                                     |
+------------------------+------------------------------------------------------------+
| phbdys                 | PHBDY                                                      |
+------------------------+------------------------------------------------------------+
| plotels                | PLOTEL                                                     |
+------------------------+------------------------------------------------------------+
| points                 | POINT                                                      |
+------------------------+------------------------------------------------------------+
| properties             | PBAR, PBARL, PBCOMP, PBEAM, PBEAML, PBEND, PBMSECT,        |
+------------------------+------------------------------------------------------------+
|                        | PBRSECT, PBUSH, PBUSH1D, PCOMP, PCOMPG, PCOMPS, PCONEAX,   |
+------------------------+------------------------------------------------------------+
|                        | PDAMP, PDAMP5, PELAS, PFAST, PGAP, PIHEX, PLPLANE,         |
+------------------------+------------------------------------------------------------+
|                        | PLSOLID, PPLANE, PRAC2D, PRAC3D, PROD, PSHEAR, PSHELL,     |
+------------------------+------------------------------------------------------------+
|                        | PSOLID, PTUBE, PVISC                                       |
+------------------------+------------------------------------------------------------+
| properties_mass        | PMASS                                                      |
+------------------------+------------------------------------------------------------+
| qsets                  | QSET, QSET1                                                |
+------------------------+------------------------------------------------------------+
| radcavs                | RADCAV, RADLST                                             |
+------------------------+------------------------------------------------------------+
| radmtx                 | RADMTX                                                     |
+------------------------+------------------------------------------------------------+
| radset                 | RADSET                                                     |
+------------------------+------------------------------------------------------------+
| random_tables          | TABRND1, TABRNDG                                           |
+------------------------+------------------------------------------------------------+
| rigid_elements         | RBAR, RBAR1, RBE1, RBE2, RBE3, RROD, RSPLINE, RSSCON       |
+------------------------+------------------------------------------------------------+
| ringaxs                | POINTAX, RINGAX                                            |
+------------------------+------------------------------------------------------------+
| ringfl                 | RINGFL                                                     |
+------------------------+------------------------------------------------------------+
| rotors                 | ROTORD, ROTORG                                             |
+------------------------+------------------------------------------------------------+
| se_bsets               | SEBSET, SEBSET1                                            |
+------------------------+------------------------------------------------------------+
| se_csets               | SECSET, SECSET1                                            |
+------------------------+------------------------------------------------------------+
| se_qsets               | SEQSET, SEQSET1                                            |
+------------------------+------------------------------------------------------------+
| se_sets                | SESET                                                      |
+------------------------+------------------------------------------------------------+
| se_suport              | SESUP                                                      |
+------------------------+------------------------------------------------------------+
| se_usets               | SEQSET1                                                    |
+------------------------+------------------------------------------------------------+
| sebndry                | SEBNDRY                                                    |
+------------------------+------------------------------------------------------------+
| sebulk                 | SEBULK                                                     |
+------------------------+------------------------------------------------------------+
| seconct                | SECONCT                                                    |
+------------------------+------------------------------------------------------------+
| seelt                  | SEELT                                                      |
+------------------------+------------------------------------------------------------+
| seexcld                | SEEXCLD                                                    |
+------------------------+------------------------------------------------------------+
| selabel                | SELABEL                                                    |
+------------------------+------------------------------------------------------------+
| seload                 | SELOAD                                                     |
+------------------------+------------------------------------------------------------+
| seloc                  | SELOC                                                      |
+------------------------+------------------------------------------------------------+
| sempln                 | SEMPLN                                                     |
+------------------------+------------------------------------------------------------+
| senqset                | SENQSET                                                    |
+------------------------+------------------------------------------------------------+
| seqgp                  | SEQGP                                                      |
+------------------------+------------------------------------------------------------+
| setree                 | SETREE                                                     |
+------------------------+------------------------------------------------------------+
| sets                   | SET1, SET3                                                 |
+------------------------+------------------------------------------------------------+
| spcadds                | SPCADD                                                     |
+------------------------+------------------------------------------------------------+
| spcoffs                | SPCOFF, SPCOFF1                                            |
+------------------------+------------------------------------------------------------+
| spcs                   | GMSPC, SPC, SPC1, SPCAX                                    |
+------------------------+------------------------------------------------------------+
| splines                | SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5, SPLINE6,      |
+------------------------+------------------------------------------------------------+
|                        | SPLINE7                                                    |
+------------------------+------------------------------------------------------------+
| suport                 | SUPORT                                                     |
+------------------------+------------------------------------------------------------+
| suport1                | SUPORT1                                                    |
+------------------------+------------------------------------------------------------+
| tables                 | TABLEH1, TABLEHT, TABLES1, TABLEST                         |
+------------------------+------------------------------------------------------------+
| tables_d               | TABLED1, TABLED2, TABLED3, TABLED4                         |
+------------------------+------------------------------------------------------------+
| tables_m               | TABLEM1, TABLEM2, TABLEM3, TABLEM4                         |
+------------------------+------------------------------------------------------------+
| tables_sdamping        | TABDMP1                                                    |
+------------------------+------------------------------------------------------------+
| tempds                 | TEMPD                                                      |
+------------------------+------------------------------------------------------------+
| thermal_materials      | MAT4, MAT5                                                 |
+------------------------+------------------------------------------------------------+
| tics                   | TIC                                                        |
+------------------------+------------------------------------------------------------+
| transfer_functions     | TF                                                         |
+------------------------+------------------------------------------------------------+
| trims                  | TRIM, TRIM2                                                |
+------------------------+------------------------------------------------------------+
| tstepnls               | TSTEP1, TSTEPNL                                            |
+------------------------+------------------------------------------------------------+
| tsteps                 | TSTEP                                                      |
+------------------------+------------------------------------------------------------+
| usets                  | USET, USET1                                                |
+------------------------+------------------------------------------------------------+
| view3ds                | VIEW3D                                                     |
+------------------------+------------------------------------------------------------+
| views                  | VIEW                                                       |
+------------------------+------------------------------------------------------------+

- Executive Control Deck
- System Control Deck
- Case Control Deck
- cross-referencing to simplify accessing data
   - *_ref attributes are cross-referenced
   - element.nodes is not cross-referenced
   - element.nodes_ref is cross-referenced
- safe cross-referencing for imperfect models
- optional error storage to get a list of all discovered errors as once
- model validation



OP4 Reader
==========
 - For matrices, the OP2 is preffered.  It's simply faster.
 - Types:

   - ASCII/binary
   - SMALL/BIG MAT format
   - Real/Complex
   - Sparse/Dense
   - Single/Double Precision
 - ASCII writer

OP2 Reader / F06 Writer
=======================
- Supported Nastran versions:

  - MSC Nastran
  - NX Nastran
  - Optistruct
  - Radioss
  - IMAT
  - Autodesk Nastran/Nastran-in-CAD

    - geometry not supported

- Input/Output:

  - Very fast OP2 reader (up to 500 MB/sec with an SSD)
  - Memory efficient
  - support directly loading into HDF5 for very large models
  - HDF5 export/import support for MATLAB integration
  - pandas support (results & matrices)
  - F06 writing
  - Most fatal errors caught (BDF input errors not caught)
  - geometry can be read directly from op2 (it's not perfect, but it's much faster)

- Operations:

  - transform displacement/eigenvectors/spc/mpc/applied loads to global coordinate system
  - transform stresses/forces to material coordinate system

- Supports:

  - superelements
  - optimization
  - mesh adaptivity
  - preload
  - shape optimization

OP2 Results
------------
- This is probably an incomplete list.  **Most** results are supported.
- **Basic Tables**

  - Types:
     - Displacment
     - Velocity
     - Acceleration
     - Eigenvectors
     - SPC/MPC Forces
     - Applied Loads
     - Load Vectors
     - Temperature
  - Real/Complex
  - Random; no NO (Number of Crossings) or RMS results

- **Stress/Strain**

  - Real/Complex
  - Random; no NO (Number of Crossings) or RMS results
  - Types:

     - Spring, Rod, Bar, Beam, Bushing, Gap, Shell, Solid

- **Forces**

  - Real/Complex
  - Types:

     - Loads: Spring, Rod, Bar, Beam, Bushing, Gap, Shell (Isotropic/Composite), Solid
     - Thermal Gradient/Flux: 1D, 2D, 3D

- **Grid Point Forces**

  - Real/Complex

- **Strain Energy**

  - Real/Complex
  - Types:

    - Spring, Rod, Bar, Beam, Bushing, Gap, Shell (Isotropic/Composite), Solid, Rigid, DMIG

- **Matrices**

  - Basic:

    - Real/Complex
    - Sparse/Dense
    - Single/Double Precision

  - MATPOOL:

    - Real/Complex
    - Sparse/Dense
    - Single/Double Precision

- Other

  - Eigenvalues

    - Modal, Buckling, Complex

  - Grid Point Weight
  - Monitor Points
  - Design Optimization:

    - Convergence History
    - **Limited** Design Responses:

      - Weight
      - Stress (Isotropic/Composite)
      - Strain (Isotropic/Composite)
      - Force
      - Flutter

F06 Plotter
===========
- flutter (SOL 145) parser

  - Supports:
     - multiple subcases
     - PK and PKNL methods

  - `plot_Vg_Vf(...)`, `plot_Vg(...)`, `plot_root_locus(...)`
  - input/output units

GUI
========
   [GUI](http://pynastran-git.readthedocs.io/en/latest/quick_start/gui.html)
 - buttons for picking, rotation center, distance, min/max
 - GUI Features:

   - Packages:

     - PyQt4/PyQt5
     - PySide/PySide2
     - QScintilla & pygments support for scripting code editor
   - color coded logging

   - legend menu

     - min/max control
     - number of labels/colors
     - additional color maps
     - legend position

   - animation menu

      - mix and match fringe/displacement/vector results (e.g., stress shown on a displaced model)
      - Real/Complex Results
          - Scale factor
          - Phase
          - Time
      - Multiple Animation Profiles
      - Where:

        - in GUI
        - exported gif

   - node/element highlighting
   - element groups
   - high resolution screenshots
   - nodal/centroidal picking
   - coordinate systems
   - results sidebar
   - custom user results

     - nodal fringe
     - centroidal fringe
     - deflection
     - nodal vector results (e.g., SPC forces)
   - preferences menu

Nastran Specific Features
-------------------------
- multiple OP2s
- deflection plots
- SOL 200 support

- geometry

  - all elements supported in BDF

- bar profile visualzation

  - 3D
  - dimensional vectors

- aero models

  - CAERO panels & subpanels
  - sideslip coordinate systems support

- mass elements
- plotting elements (e.g., PLOTEL)
- nominal geometry (useful for deflection plots)

Nastran Geometry Results
^^^^^^^^^^^^^^^^^^^^^^^^
- node id
- element id
- property id

  - PSHELL breakdown

    - thickness, ts/t, 12I/t^3
    - for each material:

      - material id
      - stiffnesses
      - is_isotropic

  - PCOMP breakdown

    - total thickness
    - for each layer:

      - thickness
      - material id
      - stiffnesses
      - is_isotropic

  - PSOLID breakdown

    - material id
    - stiffnesses
    - is_isotropic

 - loads
 - optimization

   - design regions
   - current value
   - lower/upper bounds

 - mesh quality:

   - area, min/max interior angle, skew angle, aspect ratio, taper ratio results

Nastran OP2 Results
^^^^^^^^^^^^^^^^^^^
- solution types:

  - analysis types:

    - static
    - modal
    - frequency response
    - load step

  - additional model complexity

    - optimization
    - preload

- result quantities:

  - displacement, velocity, acceleration, eigenvectors
  - SPC/MPC forces
  - applied loads
  - temperature
  - stress/strain
  - strain energy
  - limited element forces
  - thermal gradient/flux

Converters / Additional GUI Options
-----------------------------------
pyNastran's code base makes it easy to develop other useful tools
that make use of common code.  As such, additional formats are supported
in terms of readers/writers/converters/viewing, but are not a main focus.

These include:

- AFLR
- AVL
- Cart3d
- Panair
- OpenFOAM
- S/HABP
- LAWGS
- FAST
- STL
- SU2
- Tetgen
- Tecplot
- Usm3d
- Abaqus
