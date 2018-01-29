from __future__ import print_function, absolute_import

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.materials import MAT1, MAT8


class OrthotropicMaterial(object):
    def __init__(self):
        self.mid = 0
        self.material_name = ''
        self.material_type = ''
        self.material_product = ''
        self.env_condition = ''
        self.application = ''

        self.e11 = 0.
        self.e22 = 0.
        self.g12 = 0.
        self.nu12 = 0.
        self.nu21 = 0.
        self.nu32 = 0.

        self.e1t = 0.
        self.e1c = 0.
        self.e2t = 0.
        self.e2c = 0.
        self.e12 = 0.
        self.e3t = 0.
        self.e3c = 0.
        self.e23 = 0.

        self.f1t = 0.
        self.f1c = 0.
        self.f2t = 0.
        self.f2c = 0.
        self.f12 = 0.
        self.f3t = 0.
        self.f3c = 0.
        self.f23 = 0.

        self.c11 = 0.
        self.c22 = 0.
        self.c12 = 0.
        self.c66 = 0.

    def update(self):
        assert self.e11 > 0.
        assert self.e22 > 0.
        assert self.g12 > 0.
        # assert self.nu12 > 0.  # nastran supports nu12 == 0.

        self.nu21 = self.e22 * self.nu12 / self.e11

        m = 1. - self.nu12 * self.nu21

        assert m > 0.

        self.c11 = self.e11 / m
        self.c22 = self.e22 / m
        self.c12 = self.nu12 * self.c22
        self.c66 = self.g12

    def from_bdf(self, bdf, mid):
        # type: (BDF, int) -> None
        
        mat = bdf.materials[mid]
        card_type = mat.type
        
        if card_type == 'MAT1':
            self.from_mat1(mat)
        elif card_type == 'MAT8':
            self.from_mat8(mat)
        else:
            raise ValueError('Unsupported property %s!' % card_type)
            
    def from_mat1(self, mat):
        # type: (MAT1) -> None

        self.e11 = mat.e
        self.e22 = mat.e
        self.g12 = mat.g
        self.nu12 = mat.nu
        self.nu21 = mat.nu
        self.nu32 = mat.nu

        self.mid = mat.mid
        
    def from_mat8(self, mat):
        # type: (MAT8) -> None
        
        self.e11 = mat.e11
        self.e22 = mat.e22
        self.g12 = mat.g12
        self.nu12 = mat.nu12

        self.mid = mat.mid
