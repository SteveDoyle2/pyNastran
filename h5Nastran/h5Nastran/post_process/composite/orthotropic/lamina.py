from __future__ import print_function, absolute_import

from math import radians, sin, cos

from .material import OrthotropicMaterial


class OrthotropicLamina(object):
    def __init__(self):
        self.material = OrthotropicMaterial()
        self.orientation = 0.  # degrees
        self.thickness = 0.
        self.q11 = 0.
        self.q22 = 0.
        self.q12 = 0.
        self.q66 = 0.
        self.q16 = 0.
        self.q26 = 0.

    def update(self):
        self.material.update()

        theta = radians(self.orientation)

        c = cos(theta)
        s = sin(theta)
        c2 = c * c
        c3 = c2 * c
        c4 = c3 * c
        s2 = s * s
        s3 = s2 * s
        s4 = s3 * s

        c11 = self.material.c11
        c12 = self.material.c12
        c22 = self.material.c22
        c66 = self.material.c66

        self.q11 = c11 * c4 + 2. * (c12 + 2. * c66) * s2 * c2 + c22 * s4
        self.q22 = c11 * s4 + 2. * (c12 + 2. * c66) * s2 * c2 + c22 * c4
        self.q12 = (c11 + c22 - 4. * c66) * s2 * c2 + c12 * (s4 + c4)
        self.q66 = (c11 + c22 - 2. * c12 - 2. * c66) * s2 * c2 + c66 * (s4 + c4)
        self.q16 = (c11 - c12 - 2. * c66) * s * c3 + (c12 - c22 + 2. * c66) * s3 * c
        self.q26 = (c11 - c12 - 2. * c66) * s3 * c + (c12 - c22 + 2. * c66) * s * c3
