from math import radians,sin,cos

def polarToRealImag(mag,phase):
    """
    Converts magnitude-phase to real-imaginary
    so all complex results are consistent
    @param mag magnitude c^2
    @param phase phase angle phi (degrees; theta)
    @retval realValue the real component a of a+bi
    @retval imagValue the imaginary component b of a+bi
    """
    realValue = mag*cos(radians(phase))
    imagValue = mag*sin(radians(phase))
    return (realValue,imagValue)

