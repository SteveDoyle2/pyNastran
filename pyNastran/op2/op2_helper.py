from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
"""
defines:
 - polar_to_real_imag
 - real_imag_to_mag_phase
"""

from numpy import radians, abs, angle  # ,sin, cos
from cmath import rect  # polar


def polar_to_real_imag(mag, phase):
    """
    Converts magnitude-phase to real-imaginary
    so all complex results are consistent

    Parameters
    ----------
    mag : float
        magnitude c^2
    phase : float
        phase angle phi (degrees; theta)

    Returns
    -------
    real_value : float
        the real component a of a+bi
    imag_value : float
        the imaginary component b of a+bi
    """
    return rect(mag, radians(phase))

def real_imag_to_mag_phase(real_imag):
    """returns the magnitude and phase (degrees) of a complex number"""
    return abs(real_imag), angle(real_imag, deg=True)

# def realImagToMagPhase(real_imag):
    # """returns the magnitude and phase (degrees) of a complex number"""
    # return abs(real_imag), angle(real_imag, deg=True)
