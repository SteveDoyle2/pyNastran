solidBending
 - written by MSC.Nastran
 - BDF reader - SOL 101
 - OP2 reader - displacement, stress
 - tested in Windows GUI

plate_py
 - written by FEMAP
 - BDF reader - SOL 103
 - OP2 reader - modes on a plate
 - no results in GUI

beam_modes
 = written by FEMAP
 - BDF reader - SOL 103
 - OP2 reader - has problems with OP2 written by Linux (while on Windows), Windows OP2 is OK
   - note the bug is caused by PARAM,POST,-2
   - if you use PARAM,POST,-1 the op2 works fine
 - no results in GUI

