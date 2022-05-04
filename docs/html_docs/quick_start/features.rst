==============================
Features
==============================

Overview
========
 - Python 3.7-3.8
 - BSD-3 license
 - unicode support
 - importable from within Matlab
 - limited package requirements for BDF, OP2, and F06
 - additional features available with more packages

   - BDF/OP2:

      - h5py for HDF5 input/output support
      - PyQt5/PySide2/wxpython for file loading popup
   - OP2:

     - pandas for results/matrices for use in the Jupyter Notebook
   - F06:

     - matplotlib support for plotting
   - GUI: range of choices

     - PyQt5/PySide2
     - VTK 7-9
   - logging using **cpylog**

     - colorama for console logging
     - HTML logging for Jupyter Notebook
     - no markup when piping output to a file
     - supports overwriting logger object with user-defined logger

BDF Reader/Writer
=================
 - Input/Output:

   - 402 cards supported including:

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
   - solid skinning (free faces)
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
| acmodl                 | ACMODL                                                     |
+------------------------+------------------------------------------------------------+
| aecomps                | AECOMP, AECOMPL                                            |
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
| bfric                  | BFRIC                                                      |
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
| cyax                   | CYAX                                                       |
+------------------------+------------------------------------------------------------+
| cyjoin                 | CYJOIN                                                     |
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
| dmi                    | DMI                                                        |
+------------------------+------------------------------------------------------------+
| dmiax                  | DMIAX                                                      |
+------------------------+------------------------------------------------------------+
| dmig                   | DMIG                                                       |
+------------------------+------------------------------------------------------------+
| dmij                   | DMIJ                                                       |
+------------------------+------------------------------------------------------------+
| dmiji                  | DMIJI                                                      |
+------------------------+------------------------------------------------------------+
| dmik                   | DMIK                                                       |
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
+------------------------+------------------------------------------------------------+
| feedge                 | FEEDGE                                                     |
+------------------------+------------------------------------------------------------+
| feface                 | FEFACE                                                     |
+------------------------+------------------------------------------------------------+
| flfacts                | FLFACT                                                     |
+------------------------+------------------------------------------------------------+
| flutters               | FLUTTER                                                    |
+------------------------+------------------------------------------------------------+
| frequencies            | FREQ, FREQ1, FREQ2, FREQ3, FREQ4, FREQ5                    |
+------------------------+------------------------------------------------------------+
| gmcurv                 | GMCURV                                                     |
+------------------------+------------------------------------------------------------+
| gmsurf                 | GMSURF                                                     |
+------------------------+------------------------------------------------------------+
| grdset                 | GRDSET                                                     |
+------------------------+------------------------------------------------------------+
| gridb                  | GRIDB                                                      |
+------------------------+------------------------------------------------------------+
| gusts                  | GUST                                                       |
+------------------------+------------------------------------------------------------+
| hyperelastic_materials | MATHE, MATHP                                               |
+------------------------+------------------------------------------------------------+
| load_combinations      | CLOAD, LOAD, LSEQ                                          |
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
| modtrak                | MODTRAK                                                    |
+------------------------+------------------------------------------------------------+
| monitor_points         | MONDSP1, MONPNT1, MONPNT2, MONPNT3                         |
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
| omits                  | OMIT, OMIT1                                                |
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
| pset                   | PSET                                                       |
+------------------------+------------------------------------------------------------+
| pval                   | PVAL                                                       |
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
| release                | RELEASE                                                    |
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
| topvar                 | TOPVAR                                                     |
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

   - ``*_ref`` attributes are cross-referenced
   - element.nodes is not cross-referenced
   - element.nodes_ref is cross-referenced

- safe cross-referencing for imperfect models
- optional error storage to get a list of all discovered errors as once
- model validation



OP4 Reader/Writer
=================
 - For matrices, the OP2 is preferred.  It's simply faster.
 - Types:

   - ASCII/binary
   - SMALL/BIG MAT format
   - Real/Complex
   - Sparse/Dense
   - Single/Double Precision
 - ASCII writer

OP2 Reader / OP2 Writer / F06 Writer
====================================
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
  - OP2/F06 writing
  - Most fatal errors caught (geometry input errors not caught)
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

Main OP2 Results
----------------
The main op2 results can be accessed directly from the op2 object
(e.g., model.displacements, model.celas1_stress).

 - OUG - displacement, temperatures, eigenvectors, velocity, acceleration

  - displacements
  - velocities
  - accelerations
  - displacements_scaled
  - temperatures
  - eigenvectors

 - OQG - spc/mpc forces

  - spc_forces
  - spc_forces_v
  - spc_forces_scaled_response_spectra_nrl
  - mpc_forces
  - mpc_forces_RAQCONS
  - mpc_forces_RAQEATC
  - thermal_gradient_and_flux

 - OGF - grid point forces

  - grid_point_forces

 - OPG - summation of loads for each element

  - load_vectors
  - load_vectors_v
  - thermal_load_vectors
  - applied_loads
  - force_vectors

 - OES/OSTR

  - 0d - CELASx stress/strain

   - celas1_stress
   - celas2_stress
   - celas3_stress
   - celas4_stress
   - celas1_strain
   - celas2_strain
   - celas3_strain
   - celas4_strain

  - isotropic CROD/CONROD/CTUBE stress/strain

   - crod_stress
   - conrod_stress
   - ctube_stress
   - crod_strain
   - conrod_strain
   - ctube_strain

  - isotropic CBAR stress/strain

   - cbar_stress
   - cbar_strain
   - cbar_stress_10nodes
   - cbar_strain_10nodes

  - isotropic CBEAM stress/strain

   - cbeam_stress
   - cbeam_strain
   - nonlinear_cbeam_stress

  - CBEND

   - cbend_stress
   - cbend_strain

  - isotropic CTRIA3/CQUAD4 stress

   - ctria3_stress
   - ctriar_stress
   - ctria6_stress
   - cquadr_stress
   - cquad4_stress
   - cquad8_stress

  - isotropic CTRIA3/CQUAD4 strain

   - ctria3_strain
   - ctriar_strain
   - ctria6_strain
   - cquadr_strain
   - cquad4_strain
   - cquad8_strain

  - isotropic CTETRA/CHEXA/CPENTA stress/strain

   - ctetra_stress
   - chexa_stress
   - cpenta_stress
   - ctetra_strain
   - chexa_strain
   - cpenta_strain

  - CSHEAR stress/strain

   - cshear_stress
   - cshear_strain

  - GAPNL 86

   - nonlinear_cgap_stress

  - CBUSH 226

   - nonlinear_cbush_stress
   - cbush1d_stress_strain
   - nonlinear_cbush1d_stress_strain
   - cplstn3_stress
   - cplstn4_stress
   - cplstn6_stress
   - cplstn8_stress
   - cplsts3_stress
   - cplsts4_stress
   - cplsts6_stress
   - cplsts8_stress
   - cplstn3_strain
   - cplstn4_strain
   - cplstn6_strain
   - cplstn8_strain
   - cplsts3_strain
   - cplsts4_strain
   - cplsts6_strain
   - cplsts8_strain

  - CTRIAX6

   - ctriax_stress
   - ctriax_strain
   - cbush_stress
   - cbush_strain

  - nonlinear CROD/CONROD/CTUBE stress

   - nonlinear_crod_stress
   - nonlinear_crod_strain
   - nonlinear_ctube_stress
   - nonlinear_ctube_strain
   - nonlinear_conrod_stress
   - nonlinear_conrod_strain

  - CEALS1 224, CELAS3 225

   - nonlinear_celas1_stress
   - nonlinear_celas3_stress

  - composite CTRIA3/CQUAD4 stress

   - cquad4_composite_stress
   - cquad8_composite_stress
   - cquadr_composite_stress
   - ctria3_composite_stress
   - ctria6_composite_stress
   - ctriar_composite_stress
   - cquad4_composite_strain
   - cquad8_composite_strain
   - cquadr_composite_strain
   - ctria3_composite_strain
   - ctria6_composite_strain
   - ctriar_composite_strain

 - OESNLXR - CTRIA3/CQUAD4 stress

  - nonlinear_cquad4_stress
  - nonlinear_ctria3_stress
  - nonlinear_cquad4_strain
  - nonlinear_ctria3_strain
  - hyperelastic_cquad4_strain

 - OESNLXR - solids

   - nonlinear_ctetra_stress_strain
   - nonlinear_cpenta_stress_strain
   - nonlinear_chexa_stress_strain

 - PVT

  - params

 - LAMA

  - eigenvalues

 - HISADD

  - convergence_history

 - R1TABRG

  -response1_table

 - OEF - Forces

  - 0-d

   - celas1_force
   - celas2_force
   - celas3_force
   - celas4_force
   - cvisc_force
   - coneax_force
   - cdamp1_force
   - cdamp2_force
   - cdamp3_force
   - cdamp4_force
   - cgap_force

  - rod

   - crod_force
   - conrod_force
   - ctube_force

 - bar/beam

  - cbar_force
  - cbar_force_abs
  - cbar_force_srss
  - cbar_force_nrl
  - cbar_force_10nodes
  - cbeam_force
  - cbeam_force_vu
  - cbush_force
  - cbend_force

 - shell

  - cquad4_force
  - cquad8_force
  - cquadr_force
  - ctria3_force
  - ctria6_force
  - ctriar_force
  - cshear_force

 - solid

  - chexa_pressure_force
  - cpenta_pressure_force
  - ctetra_pressure_force
  - vu_quad_force
  - vu_tria_force

 - OEF - Fluxes

  - conv_thermal_load
  - chbdye_thermal_load
  - chbdye_thermal_load_flux
  - chbdyg_thermal_load
  - chbdyg_thermal_load_flux
  - chbdyp_thermal_load
  - chbdyp_thermal_load_flux

  - thermalLoad_1D

   - crod_thermal_load
   - crod_thermal_load_flux
   - cbeam_thermal_load
   - cbeam_thermal_load_flux
   - ctube_thermal_load
   - ctube_thermal_load_flux
   - conrod_thermal_load
   - conrod_thermal_load_flux
   - cbar_thermal_load
   - cbar_thermal_load_flux
   - cbend_thermal_load
   - cbend_thermal_load_flux

  - thermalLoad_2D_3D

   - cquad4_thermal_load
   - cquad4_thermal_load_flux
   - ctriax6_thermal_load
   - ctriax6_thermal_load_flux
   - cquad8_thermal_load
   - cquad8_thermal_load_flux
   - ctria3_thermal_load
   - ctria3_thermal_load_flux
   - ctria6_thermal_load
   - ctria6_thermal_load_flux
   - ctetra_thermal_load
   - ctetra_thermal_load_flux
   - chexa_thermal_load
   - chexa_thermal_load_flux
   - cpenta_thermal_load
   - cpenta_thermal_load_flux
   - thermalLoad_VU
   - thermalLoad_VU_3D
   - vu_beam_thermal_load

 - OEFIT - Failure Indices

  - cquad4_composite_force_failure_indicies
  - cquad8_composite_force_failure_indicies
  - ctria3_composite_force_failure_indicies
  - ctria6_composite_force_failure_indicies

 - OGS1 - Grid Point Stresses

  - grid_point_surface_stresses
  - grid_point_stresses_volume_direct
  - grid_point_stresses_volume_principal
  - grid_point_stress_discontinuities

 - OEE - Strain Energy Density

  - cquad4_strain_energy
  - cquad8_strain_energy
  - cquadr_strain_energy
  - cquadx_strain_energy
  - ctria3_strain_energy
  - ctria6_strain_energy
  - ctriar_strain_energy
  - ctriax_strain_energy
  - ctriax6_strain_energy
  - cshear_strain_energy
  - ctetra_strain_energy
  - cpenta_strain_energy
  - chexa_strain_energy
  - cpyram_strain_energy
  - crod_strain_energy
  - ctube_strain_energy
  - conrod_strain_energy
  - cbar_strain_energy
  - cbeam_strain_energy
  - cgap_strain_energy
  - cbush_strain_energy
  - celas1_strain_energy
  - celas2_strain_energy
  - celas3_strain_energy
  - celas4_strain_energy
  - cdum8_strain_energy
  - cbend_strain_energy
  - dmig_strain_energy
  - genel_strain_energy
  - conm2_strain_energy
  - rbe1_strain_energy
  - rbe3_strain_energy

 - unused?

  - displacement_scaled_response_spectra_nrl
  - displacement_scaled_response_spectra_abs
  - displacement_scaled_response_spectra_srss
  - velocity_scaled_response_spectra_abs
  - acceleration_scaled_response_spectra_nrl
  - acceleration_scaled_response_spectra_abs

OP2.Results()
-------------

The OP2.Results() are accessed using model.results. as a prefix
(e.g., model.results.modal_contribution.celas1_stress).

 - eqexin
 - gpdt
 - bgpdt
 - ato # AutoCorrelationObjects()           - see below
 - psd # PowerSpectralDensityObjects()      - see below
 - rms # RootMeansSquareObjects()           - see below
 - no  # NumberOfCrossingsObjects()         - see below
 - crm # CumulativeRootMeansSquareObjects() - see below
 - modal_contribution

   - celas1_stress
   - celas2_stress
   - celas3_stress
   - celas4_stress
   - celas1_strain
   - celas2_strain
   - celas3_strain
   - celas4_strain
   - crod_stress
   - conrod_stress
   - ctube_stress
   - crod_strain
   - conrod_strain
   - ctube_strain
   - ctetra_stress
   - cpenta_stress
   - chexa_stress
   - ctetra_strain
   - cpenta_strain
   - chexa_strain
   - cbar_stress
   - cbar_strain
   - cbeam_stress
   - cbeam_strain
   - ctria3_stress
   - ctria6_stress
   - cquad4_stress
   - cquad8_stress
   - cquadr_stress
   - ctriar_stress
   - ctria3_strain
   - ctria6_strain
   - cquad4_strain
   - cquad8_strain
   - cquadr_strain
   - ctriar_strain
   - cquad4_composite_stress
   - cquad8_composite_stress
   - cquadr_composite_stress
   - ctria3_composite_stress
   - ctria6_composite_stress
   - ctriar_composite_stress
   - cquad4_composite_strain
   - cquad8_composite_strain
   - cquadr_composite_strain
   - ctria3_composite_strain
   - ctria6_composite_strain
   - ctriar_composite_strain
   - cshear_stress
   - cshear_strain
   - cshear_force
   - cbush_stress
   - cbush_strain
 - strength_ratio

   - cquad4_composite_stress
   - cquad8_composite_stress
   - cquadr_composite_stress
   - ctria3_composite_stress
   - ctria6_composite_stress
   - ctriar_composite_stress
   - cquad4_composite_strain
   - cquad8_composite_strain
   - cquadr_composite_strain
   - ctria3_composite_strain
   - ctria6_composite_strain
   - ctriar_composite_strain
 - ROUGV1  # relative disp/vel/acc/eigenvectors

   - displacements
   - velocities
   - accelerations
   - eigenvectors
 - RADEFFM

   - eigenvectors
 - RADCONS

   - eigenvectors
 - RAFCONS

   - cbar_force
   - cquad4_force
   - cbush_force
 - RASCONS

   - ctetra_stress
   - cpenta_stress
   - chexa_stress
   - ctetra_strain
   - cpenta_strain
   - chexa_strain
   - ctria3_stress
   - ctria6_stress
   - cquad4_stress
   - cquad8_stress
   - cquadr_stress
   - ctriar_stress
   - ctria3_strain
   - ctria6_strain
   - cquad4_strain
   - cquad8_strain
   - cquadr_strain
   - ctriar_strain
 - RAECONS

   - ctria3_strain
   - cquad4_strain
   - chexa_strain
 - RAGCONS

   - grid_point_forces
 - RAPCONS

   - cquad4_composite_stress
   - cquad8_composite_stress
   - cquadr_composite_stress
   - ctria3_composite_stress
   - ctria6_composite_stress
   - ctriar_composite_stress
 - RANCONS

   - cbar_strain_energy
   - cbush_strain_energy
   - chexa_strain_energy
   - ctria3_strain_energy
   - cquad4_strain_energy
 - RADEATC

   - eigenvectors
 - RAFEATC

   - cbar_force
   - cquad4_force
   - cbush_force
 - RASEATC

   - chexa_stress
   - cquad4_stress
 - RAEEATC

   - chexa_strain
   - ctria3_strain
   - cquad4_strain
 - RAGEATC

   - grid_point_forces
 - RAPEATC

   - cquad4_composite_stress
   - cquad8_composite_stress
   - cquadr_composite_stress
   - ctria3_composite_stress
   - ctria6_composite_stress
   - ctriar_composite_stress
 - RANEATC

   - cbar_strain_energy
   - cbush_strain_energy
   - chexa_strain_energy
   - ctria3_strain_energy
   - cquad4_strain_energy

All of these results have the same outputs (shown under model.results.crm).
For example, model.results.ato.displacements, model.results.crm.displacements.

 - ato # AutoCorrelationObjects()
 - psd # PowerSpectralDensityObjects()
 - rms # RootMeansSquareObjects()
 - no  # NumberOfCrossingsObjects()
 - crm # CumulativeRootMeansSquareObjects()

   - displacements
   - velocities
   - accelerations
   - load_vectors
   - spc_forces
   - mpc_forces
   - crod_force
   - conrod_force
   - ctube_force
   - cbar_force
   - cbeam_force
   - cbush_stress
   - cbush_strain
   - crod_stress
   - conrod_stress
   - ctube_stress
   - cbar_stress
   - cbeam_stress
   - crod_strain
   - conrod_strain
   - ctube_strain
   - cbar_strain
   - cbeam_strain
   - ctetra_strain
   - cpenta_strain
   - chexa_strain
   - ctetra_stress
   - cpenta_stress
   - chexa_stress
   - celas1_stress
   - celas2_stress
   - celas3_stress
   - celas4_stress
   - celas1_strain
   - celas2_strain
   - celas3_strain
   - celas4_strain
   - celas1_force
   - celas2_force
   - celas3_force
   - celas4_force
   - ctria3_force
   - ctria6_force
   - ctriar_force
   - cquad4_force
   - cquad8_force
   - cquadr_force
   - ctria3_stress
   - ctria6_stress
   - cquad4_stress
   - cquad8_stress
   - cquadr_stress
   - ctriar_stress
   - ctria3_strain
   - ctria6_strain
   - cquad4_strain
   - cquad8_strain
   - cquadr_strain
   - ctriar_strain
   - cbend_stress
   - cbend_strain
   - cbend_force
   - cshear_stress
   - cshear_strain
   - cshear_force
   - cbush_force
   - cdamp1_force
   - cdamp2_force
   - cdamp3_force
   - cdamp4_force
   - cvisc_force
   - cquad4_composite_stress
   - cquad8_composite_stress
   - cquadr_composite_stress
   - ctria3_composite_stress
   - ctria6_composite_stress
   - ctriar_composite_stress
   - cquad4_composite_strain
   - cquad8_composite_strain
   - cquadr_composite_strain
   - ctria3_composite_strain
   - ctria6_composite_strain
   - ctriar_composite_strain

Matrices with explicit methods
------------------------------
These are simply accessor methods to various matrices.  For example,
``model.total_effective_mass_matrix`` is the same as ``self.matrices['EFMFSMS']``.

 - total_effective_mass_matrix (EFMFSMS)
 - effective_mass_matrix (EFMASSS)
 - rigid_body_mass_matrix (RBMASS)
 - modal_effective_mass_fraction (EFMFACS)
 - modal_participation_factors (MPFACS)
 - modal_effective_mass (MEFMASS)
 - modal_effective_weight (MEFWTS)

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
