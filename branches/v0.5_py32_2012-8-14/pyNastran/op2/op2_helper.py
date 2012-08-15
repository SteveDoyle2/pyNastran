## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from numpy import radians, abs, angle #,sin, cos
#from math import radians
from cmath import rect #polar

def polarToRealImag(mag, phase):
    """
    Converts magnitude-phase to real-imaginary
    so all complex results are consistent
    @param mag magnitude c^2
    @param phase phase angle phi (degrees; theta)
    @retval realValue the real component a of a+bi
    @retval imagValue the imaginary component b of a+bi
    """
    return rect(mag,radians(phase))
    #realValue = mag*cos(radians(phase)) # phase in degrees
    #imagValue = mag*sin(radians(phase))
    #return complex(realValue, imagValue)

def realImagToMagPhase(realImag):
    """returns the magnitude and phase (degrees) of a complex number"""
    return abs(realImag), angle(realImag, deg=True)
