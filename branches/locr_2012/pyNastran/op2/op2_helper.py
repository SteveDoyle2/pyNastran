from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from numpy import radians, abs, angle  # ,sin, cos
#from math import radians
from cmath import rect  # polar


def polarToRealImag(mag, phase):
    """
    Converts magnitude-phase to real-imaginary
    so all complex results are consistent
    @param mag magnitude c^2
    @param phase phase angle phi (degrees; theta)
    @retval realValue the real component a of a+bi
    @retval imagValue the imaginary component b of a+bi
    """
    return rect(mag, radians(phase))
    #realValue = mag*cos(radians(phase)) # phase in degrees
    #imagValue = mag*sin(radians(phase))
    #return complex(realValue, imagValue)


def realImagToMagPhase(realImag):
    """returns the magnitude and phase (degrees) of a complex number"""
    return abs(realImag), angle(realImag, deg=True)
