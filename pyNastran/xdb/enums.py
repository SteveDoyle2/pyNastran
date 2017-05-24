from enum import Enum

class SolutionType(Enum):

    ZERO = 0 #Old data recovery (refer to prototype for "best guess" logic)
    STATIC = 1 #Static Data Recovery
    TRANSIENT = 2 #Transient Response Data Recovery
    FREQUENCY = 3 #Frequency Response Data Recovery
    REAL_EIGENVALUE = 4 #Real Eigenvalue Data Recovery
    NONLINEAR_STATIC = 5 #Nonlinear Statics Data Recovery
    BUCKLING = 6 #Buckling Data Recovery
    COMPLEX_EIGENVALUE = 7 #Complex Eigenvalue Data Recovery
    AERODYNAMIC = 8 #Aerodynamic Analysis
    RESPONSE_SPECTRA = 21 #Scaled Response Spectra Analysis - ABS