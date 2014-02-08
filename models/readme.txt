bar3truss
 - taken from OpenMDAO
 - https://github.com/OpenMDAO-Plugins/nastranwrapper
 - demonstrates OpenMDAO syntax
 
solid_bending
 - written by MSC.Nastran
 - BDF reader - SOL 101
 - OP2 reader - displacement, stress
 - tested in Windows GUI
 - tested using the test_bdf and test_op2 codes

plate_py
 - written by FEMAP
 - BDF reader - SOL 103
 - OP2 reader - modes on a plate
 - no results in GUI
 - tested using the test_bdf and test_op2 codes

beam_modes
 = written by FEMAP
 - BDF reader - SOL 103
 - OP2 reader - has problems with OP2 written by Linux (while on Windows), Windows OP2 is OK
   - note the bug is caused by PARAM,POST,-2
   - if you use PARAM,POST,-1 the op2 works fine
 - no results in GUI
 - tested using the test_bdf and test_op2 codes

plate
 - written using Patran 2012 64-Bit
 - it has two versions:
     plate_formatted.bdf
       ASSIGN OUTPUT2 = 'plate.op2', UNIT = 12, FORM = FORMATTED
     plate_unformatted.bdf
       $ASSIGN OUTPUT2 = 'plate.op2', UNIT = 12, FORM = FORMATTED
  - only the unormatted op2 works in the op2 reader
 - tested using the bdf_unit_tests, f06_unit_tests, and op2_unit_tests codes


contact
 - op2
   - skips CONTACT table, but reads stresses and displacements