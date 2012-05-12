from math import radians,sin,cos,atan2,degrees

def polarToRealImag(mag,phase):
    """
    Converts magnitude-phase to real-imaginary
    so all complex results are consistent
    @param mag magnitude c^2
    @param phase phase angle phi (degrees; theta)
    @retval realValue the real component a of a+bi
    @retval imagValue the imaginary component b of a+bi
    """
    realValue = mag*cos(radians(phase)) # is phase in degrees/radians?
    imagValue = mag*sin(radians(phase))
    return complex(realValue,imagValue)

def realImagToMagPhase(realImag):
    """returns the magnitude and phase (degrees) of a complex number"""
    return abs(realImag),angle(realImag)

def angle(realImag):
    """returns the phase angle in degrees"""
    return degrees(atan2(realImag.imag,realImag.real))
