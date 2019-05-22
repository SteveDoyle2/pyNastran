$ AERODYNAMIC MODULE
$
$*******************************************************************************
$ WING SECTION
$
$ AERODYNAMIC COORDINATE SYSTEM
CORD2R,2,0,0.,0.,0.,0.,0.,-1.,+CORD200
+CORD200,-1.,0.,-1.
$
$ AERODYNAMIC PARAMETERS
AERO,2,1.,4.,1.225, 1
$
$ AERODYNAMIC PANNELS
CAERO1,601,701,2,20,10, , ,1,+CAERO10
+CAERO10,-2.,0.,0.,5.7,-1.11,12.7,0.,2.5
$
$ AERODYNAMIC PROPERTIES
PAERO1,701     $ SPECIFIES TO USE THE DLM
$
$*******************************************************************************
$ TAIL SECTION
$
$ COORDINATE SYSTEM FOR TAIL SPLINE
CORD2R,3,0,-13.,0.,-0.5,-13.,0.,-1.,+CORD201
+CORD201,-14.,0.,0.
$
$ AERODYNAMIC PANNELS
CAERO1,901,701,3,7,4, , ,1,+CAERO11
+CAERO11,-2.,0.,0.,3.5,-1.5,5.,0.,2. 
