#ifdef SWIG

%module bdf
%{
/* Includes the header in the wrapper code */
#include "coord.h"
#include "grid.h"
%}
#endif
 

/* Parse the header file to generate wrappers */
%include "coord.h"
%include "grid.h"
