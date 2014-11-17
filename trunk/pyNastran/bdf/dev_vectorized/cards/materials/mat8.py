"""
All material cards are defined in this file.  This includes:

 * CREEP
 * MAT1 (isotropic solid/shell)
 * MAT2 (anisotropic)
 * MAT3 (linear orthotropic)
 * MAT4 (thermal)
 * MAT5 (thermal)
 * MAT8 (orthotropic shell)
 * MAT9 (anisotropic solid)
 * MAT10 (fluid element)
 * MATHP (hyperelastic)

All cards are Material objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import zeros, array, where, asarray, argsort, searchsorted

from pyNastran.bdf.fieldWriter import set_blank_if_default
#from pyNastran.bdf.cards.baseCard import BaseCard, Material
#from pyNastran.bdf.cards.tables import Table
from .mat1 import Material
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank,
    string, string_or_blank, blank)
from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16


class MAT8(Material):
    """
    Defines the material property for an orthotropic material for isoparametric
    shell elements.

    +-----+-----+-----+------+------+-----+-----+-----+-----+
    |  1  |  2  |  3  |  4   |  5   |  6  |  7  |  8  |  9  |
    +=====+=====+=====+======+======+=====+=====+=====+=====+
    |MAT8 | MID | E1  |  E2  | NU12 | G12 | G1Z | G2Z | RHO |
    +-----+-----+-----+------+------+-----+-----+-----+-----+
    |     | A1  |  A2 | TREF |  Xt  |  Xc |  Yt |  Yc |  S  |
    +-----+-----+-----+------+------+-----+-----+-----+-----+
    |     | GE1 | F12 | STRN |      |     |     |     |     |
    +-----+-----+-----+------+------+-----+-----+-----+-----+
    """
    type = 'MAT8'
    #_field_map = {
        #1: 'mid', 2:'e11', 3:'e22', 4:'nu12', 5: 'g12', 6:'g1z', 7:'g2z',
        #8: 'rho', 9:'a1', 10:'a2', 11:'TRef', 12:'Xt', 13:'Xc', 14:'Yt',
        #15:'Yc', 16: 'S', 17:'ge', 18:'F12', 19:'strn',
    #}

    def __init__(self, model):
        Material.__init__(self, model)
        self.i = 0
        self.material_id = None
        #self.model.log.info('start MAT8')

    def allocate(self, ncards):
        float_fmt = self.model.float
        self.material_id = zeros(ncards, dtype='int32')
        self.e11 = zeros(ncards, dtype=float_fmt)
        self.e22 = zeros(ncards, dtype=float_fmt)
        self.nu12 = zeros(ncards, dtype=float_fmt)

        self.g12 = zeros(ncards, dtype=float_fmt)
        self.g1z = zeros(ncards, dtype=float_fmt)
        self.g2z = zeros(ncards, dtype=float_fmt)
        self.rho = zeros(ncards, dtype=float_fmt)
        self.a1 = zeros(ncards, dtype=float_fmt)
        self.a2 = zeros(ncards, dtype=float_fmt)
        self.TRef = zeros(ncards, dtype=float_fmt)
        self.Xt = zeros(ncards, dtype=float_fmt)
        self.Xc = zeros(ncards, dtype=float_fmt)
        self.Yt = zeros(ncards, dtype=float_fmt)
        self.Yc = zeros(ncards, dtype=float_fmt)
        self.S = zeros(ncards, dtype=float_fmt)
        self.ge = zeros(ncards, dtype=float_fmt)
        self.F12 = zeros(ncards, dtype=float_fmt)
        self.strn = zeros(ncards, dtype=float_fmt)

    def add(self, card=None, data=None, comment=''):
        i = self.i

        self.mats8 = None
        self.matt8 = None
        #if comment:
            #self._comment = comment
        self.material_id[i] = integer(card, 1, 'material_id')
        self.e11[i] = double(card, 2, 'E11')    #: ..todo:: is this the correct default
        self.e22[i] = double(card, 3, 'E22')    #: ..todo:: is this the correct default
        self.nu12[i] = double(card, 4, 'nu12')  #: ..todo:: is this the correct default

        self.g12[i] = double_or_blank(card, 5, 'g12', 0.0)
        self.g1z[i] = double_or_blank(card, 6, 'g1z', 1e8)
        self.g2z[i] = double_or_blank(card, 7, 'g2z', 1e8)
        self.rho[i] = double_or_blank(card, 8, 'rho', 0.0)
        self.a1[i] = double_or_blank(card, 9, 'a1', 0.0)
        self.a2[i] = double_or_blank(card, 10, 'a2', 0.0)
        self.TRef[i] = double_or_blank(card, 11, 'TRef', 0.0)
        self.Xt[i] = double_or_blank(card, 12, 'Xt', 0.0)
        self.Xc[i] = double_or_blank(card, 13, 'Xc', self.Xt[i])
        self.Yt[i] = double_or_blank(card, 14, 'Yt', 0.0)
        self.Yc[i] = double_or_blank(card, 15, 'Yc', self.Yt[i])
        self.S[i] = double_or_blank(card, 16, 'S', 0.0)
        self.ge[i] = double_or_blank(card, 17, 'ge', 0.0)
        self.F12[i] = double_or_blank(card, 18, 'F12', 0.0)
        self.strn[i] = double_or_blank(card, 19, 'strn', 0.0)
        assert len(card) <= 20, 'len(MAT8 card) = %i' % len(card)
        self.i += 1

    def build(self):
        self.n = self.i
        if self.n:
            i = argsort(self.material_id)
            self.model.log.info('MAT8.n = %s' % self.n)
            self.material_id = self.material_id[i]
            self.e11 = self.e11[i]
            self.e22 = self.e22[i]
            self.nu12 = self.nu12[i]

            self.g12 = self.g12[i]
            self.g1z = self.g1z[i]
            self.g2z = self.g2z[i]
            self.rho = self.rho[i]
            self.a1 = self.a1[i]
            self.a2 = self.a2[i]
            self.TRef = self.TRef[i]
            self.Xt = self.Xt[i]
            self.Xc = self.Xc[i]
            self.Yt = self.Yt[i]
            self.Yc = self.Yc[i]
            self.S = self.S[i]
            self.ge = self.ge[i]
            self.F12 = self.F12[i]
            self.strn = self.strn[i]

    def get_density_by_material_id(self, material_id=None):
        i = self.get_index_by_material_id(material_id)
        density = self.rho[i]
        return density

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the MAT8 object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        mid = self.Mid()
        E11 = self.E11()
        E22 = self.E22()
        nu12 = self.Nu12()
        G12 = self.G12()
        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(E11, float), 'E11=%r' % E11
        assert isinstance(E22, float), 'E11=%r' % E11
        assert isinstance(G12, float), 'G12=%r' % G12
        assert isinstance(nu12, float), 'nu12=%r' % nu12

    def E11(self):
        return self.e11

    def E22(self):
        return self.e22

    def Nu12(self):
        return self.nu12

    def G12(self):
        return self.g12

    def D(self):
        """
        .. todo:: what about G1z and G2z
        """
        E11 = self.E11()
        E22 = self.E22()
        nu12 = self.Nu12()
        G12 = self.G12()

        D = zeros((3, 3), dtype='float32')
        mu = 1. - nu12 * nu12 * E11 / E22    # not necessary b/c they're equal
        D[0, 0] = E11 / mu
        D[1, 1] = E22 / mu
        D[0, 1] = nu12 * D[0, 0]
        D[1, 0] = D[0, 1]
        D[2, 2] = G12
        return D

    def rawFields(self):
        list_fields = ['MAT8', self.mid, self.e11, self.e22, self.nu12, self.g12,
                       self.g1z, self.g2z, self.rho, self.a1, self.a2, self.TRef,
                       self.Xt, self.Xc, self.Yt, self.Yc, self.S, self.ge,
                       self.F12, self.strn]
        return list_fields

    def reprFields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the MAT8 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        return list_fields

    def write_bdf(self, f, size=8, is_double=False, material_id=None):
        if self.n:
            for mid, e11, e22, nu12, g12, g1z, g2z, rho, a1, a2, TRef, Xt, Yt, Xc, Yc, S, ge, F12, strn in zip(
                self.material_id, self.e11, self.e22, self.nu12, self.g12, self.g1z, self.g2z, self.rho,
                self.a1, self.a2, self.TRef,
                self.Xt, self.Yt, self.Xc, self.Yc, self.S, self.ge, self.F12, self.strn):
                G12 = set_blank_if_default(g12, 0.)
                G1z = set_blank_if_default(g1z, 1e8)
                G2z = set_blank_if_default(g2z, 1e8)

                rho = set_blank_if_default(rho, 0.0)
                a1 = set_blank_if_default(a1, 0.0)
                a2 = set_blank_if_default(a2, 0.0)
                TRef = set_blank_if_default(TRef, 0.0)

                Xt = set_blank_if_default(Xt, 0.)
                Yt = set_blank_if_default(Yt, 0.)

                Xc = set_blank_if_default(Xc, Xt)
                Yc = set_blank_if_default(Yc, Yt)

                S = set_blank_if_default(S, 0.0)
                ge = set_blank_if_default(ge, 0.0)
                F12 = set_blank_if_default(F12, 0.0)
                strn = set_blank_if_default(strn, 0.0)

                list_fields = ['MAT8', mid, e11, e22, nu12, G12, G1z,
                          G2z, rho, a1, a2, TRef, Xt, Xc, Yt, Yc, S, ge, F12, strn]
                if size == 8:
                    #self.comment()
                    f.write(print_card_8(list_fields))
                else:
                    f.write(print_card_16(list_fields))

    def slice_by_index(self, i):
        i = asarray(i)
        #self.model.log.debug('i = %s' % i)
        obj = MAT8(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.material_id = self.material_id[i]
        obj.e11 = self.e11[i]
        obj.e22 = self.e22[i]
        obj.nu12 = self.nu12[i]

        obj.g12 = self.g12[i]
        obj.g1z = self.g1z[i]
        obj.g2z = self.g2z[i]
        obj.rho = self.rho[i]
        obj.a1 = self.a1[i]
        obj.a2 = self.a2[i]
        obj.TRef = self.TRef[i]
        obj.Xt = self.Xt[i]
        obj.Xc = self.Xc[i]
        obj.Yt = self.Yt[i]
        obj.Yc = self.Yc[i]
        obj.S = self.S[i]
        obj.ge = self.ge[i]
        obj.F12 = self.F12[i]
        obj.strn = self.strn[i]
        return obj