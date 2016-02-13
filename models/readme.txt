Released Example Problems
-------------------------
- bar3truss
- solid_bending
- sol_101_elements
- iSat
- aero

All Examples
------------
bar3truss
 - taken from OpenMDAO
 - https://github.com/OpenMDAO-Plugins/nastranwrapper
 - demonstrates OpenMDAO syntax

solid_bending
 - written by MSC.Nastran
 - BDF reader - SOL 101
 - OP2 reader - displacement, stress
 - works in GUI

plate_py
 - written by FEMAP
 - BDF reader - SOL 103
 - transient (modal) eigenvectors
 - no problems

beam_modes
 - written by FEMAP
 - BDF reader - SOL 103
 - transient (modal) eigenvectors
 - no problems

iSat
 - works in GUI
 - transient (modal) eigenvectors/stress

plate
 - written by Patran 2012 64-Bit
 - it has two versions:
     plate_formatted.bdf
       ASSIGN OUTPUT2 = 'plate.op2', UNIT = 12, FORM = FORMATTED
     plate_unformatted.bdf
       $ASSIGN OUTPUT2 = 'plate.op2', UNIT = 12, FORM = FORMATTED
  - only the unformatted op2 works in the op2 reader


contact
 - bdf
   - tests contact cards (BCTPARA, PCTSET, BSURFS, BSURF)
   - CTETRA10s with missing nodes
   - CPYRAM solids not supported in OP2
   - MATT1
   - shows a badly formatted BCTABLE
 - op2
   - skips CONTACT table, but reads stresses and displacements

aero
 - bdf: shows off aero support (CAERO1/SPLINE5)
 - op2: incomplete, but working test problem

pload4
 - various PLOAD4 test problems 

sol200
 - 