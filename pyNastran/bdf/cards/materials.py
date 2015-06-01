# pylint: disable=C0103,C0111,C0302,R0902,R0904,R0914,E1101,W0612,E0602
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
from six import integer_types
from numpy import zeros, array

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard, Material
from pyNastran.bdf.cards.bdf_tables import Table
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string, string_or_blank, blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class IsotropicMaterial(Material):
    """Isotropic Material Class"""
    def __init__(self, card, data):
        Material.__init__(self, card, data)


class OrthotropicMaterial(Material):
    """Orthotropic Material Class"""
    def __init__(self, card, data):
        Material.__init__(self, card, data)

class AnisotropicMaterial(Material):
    """Anisotropic Material Class"""
    def __init__(self, card, data):
        Material.__init__(self, card, data)


class ThermalMaterial(Material):
    """Thermal Material Class"""
    def __init__(self, card, data):
        Material.__init__(self, card, data)


class HyperelasticMaterial(Material):
    """Hyperelastic Material Class"""
    def __init__(self, card, data):
        Material.__init__(self, card, data)


class CREEP(Material):
    type = 'CREEP'

    def __init__(self, card=None, data=None, comment=''):
        Material.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.mid = integer(card, 1, 'mid')
            self.T0 = double_or_blank(card, 2, 'T0', 0.0)
            self.exp = double_or_blank(card, 3, 'exp', 1e-9)
            self.form = string_or_blank(card, 4, 'form') # blank?
            self.tidkp = integer_or_blank(card, 5, 'tidkp') # blank?
            self.tidcp = integer_or_blank(card, 6, 'tidcp') # blank?
            self.tidcs = integer_or_blank(card, 7, 'tidcs') # blank?
            self.thresh = double_or_blank(card, 8, 'thresh', 1e-5)
            self.Type = integer_or_blank(card, 9, 'Type') # 111, 112, 121, 122, 211, 212, 221, 222, 300 (or blank?)
            self.a = double_or_blank(card, 10, 'a')
            self.b = double_or_blank(card, 11, 'b')
            self.c = double_or_blank(card, 12, 'c')
            self.d = double_or_blank(card, 13, 'd')
            self.e = double_or_blank(card, 14, 'e')
            self.f = double_or_blank(card, 15, 'f')
            self.g = double_or_blank(card, 16, 'g')
            assert len(card) <= 17, 'len(CREEP card) = %i' % len(card)
        else:
            self.mid = data[0]
            self.T0 = data[1]
            self.exp = data[2]
            self.form = data[3]
            self.tidkp = data[4]
            self.tidcp = data[5]
            self.tidcs = data[6]
            self.thresh = data[7]
            self.Type = data[8]
            self.a = data[9]
            self.b = data[10]
            self.c = data[11]
            self.d = data[12]
            self.e = data[13]
            self.f = data[14]
            self.g = data[15]

    def cross_reference(self, model):
        msg = ' which is required by CREEP pid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)

    def Mid(self):  # links up to MAT1, MAT2, MAT9 or same mid
        if isinstance(self.mid, integer_types):
            return self.mid
        return self.mid.mid

    def raw_fields(self):
        list_fields = ['CREEP', self.Mid(), self.T0, self.exp, self.form,
                       self.tidkp, self.tidcp, self.tidcs, self.thresh, self.Type,
                       self.a, self.b, self.c, self.d, self.e, self.f, self.g]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the CREEP object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        thresh = set_blank_if_default(self.thresh, 1e-5)
        exp = set_blank_if_default(self.exp, 4.1e-9)
        T0 = set_blank_if_default(self.T0, 0.0)
        list_fields = ['CREEP', self.Mid(), T0, exp, self.form, self.tidkp,
                       self.tidcp, self.tidcs, thresh, self.Type,
                       self.a, self.b, self.c, self.d, self.e, self.f, self.g]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class MAT1(IsotropicMaterial):
    """
    Defines the material properties for linear isotropic materials.

    +-----+-----+-----+-----+-------+-----+------+------+-----+
    |  1  |  2  | 3   | 4   |   5   |  6  |  7   |  8   |  9  |
    +=====+=====+=====+=====+=======+=====+======+======+=====+
    |MAT1 | MID |  E  |  G  |  NU   | RHO |  A   | TREF | GE  |
    +-----+-----+-----+-----+-------+-----+------+------+-----+
    |     | ST  | SC  | SS  | MCSID |     |      |      |     |
    +-----+-----+-----+-----+-------+-----+------+------+-----+
    """
    type = 'MAT1'
    _field_map = {
        1: 'mid', 2:'e', 3:'g', 4:'nu', 5: 'rho', 6:'a', 7:'TRef', 8:'ge',
        9: 'St', 10:'Sc', 11:'Ss', 12:'Mcsid',
    }

    def __init__(self, card=None, data=None, comment=''):
        IsotropicMaterial.__init__(self, card, data)
        self.mats1 = None
        self.matt1 = None
        if comment:
            self._comment = comment
        if card:
            self.mid = integer(card, 1, 'mid')
            self.set_E_G_nu(card)
            self.rho = double_or_blank(card, 5, 'rho', 0.)
            self.a = double_or_blank(card, 6, 'a', 0.0)
            self.TRef = double_or_blank(card, 7, 'TRef', 0.0)
            self.ge = double_or_blank(card, 8, 'ge', 0.0)
            self.St = double_or_blank(card, 9, 'St', 0.0)
            self.Sc = double_or_blank(card, 10, 'Sc', 0.0)
            self.Ss = double_or_blank(card, 11, 'Ss', 0.0)
            self.Mcsid = integer_or_blank(card, 12, 'Mcsid', 0)
            assert len(card) <= 13, 'len(MAT1 card) = %i' % len(card)
        else:
            self.mid = data[0]
            self.e = data[1]
            self.g = data[2]
            self.nu = data[3]
            self.rho = data[4]
            self.a = data[5]
            self.TRef = data[6]
            self.ge = data[7]
            self.St = data[8]
            self.Sc = data[9]
            self.Ss = data[10]
            self.Mcsid = data[11]

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the MAT1 object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        mid = self.Mid()
        E = self.E()
        G = self.G()
        nu = self.Nu()
        assert isinstance(mid, int), 'mid=%r' % mid
        if xref:
            if [self.matt1, self.mats1] == [None, None]:
                assert isinstance(E, float), 'E=%r' % E
                assert isinstance(G, float), 'G=%r' % G
                assert isinstance(nu, float), 'nu=%r' % nu

    def D(self):
        E11 = self.E()
        E22 = E11
        nu12 = self.Nu()
        G12 = self.G()

        D = zeros((3, 3))
        #D = zeros((6,6))
        mu = 1. - nu12 * nu12  # *E11/E22    # not necessary b/c they're equal
        D[0, 0] = E11 / mu
        D[1, 1] = E22 / mu
        D[0, 1] = nu12 * D[0, 0]
        D[1, 0] = D[0, 1]
        D[2, 2] = G12
        return D

    def G(self):
        return self.g

    def E(self):
        return self.e

    def Nu(self):
        return self.nu

    def Rho(self):
        return self.rho

    def get_density(self):
        return self.rho

    def E_stress(self, stress):
        if self.mats1 is not None:
            E = self.matt1.E(self.e, stress)
        else:
            E = self.e
        return E

    def E_temperature(self, temperature):
        if self.matt1 is not None:
            E = self.matt1.E(self.e, temperature)
        else:
            E = self.e
        return E

    def set_E_G_nu(self, card):
        r"""
        \f[ G = \frac{E}{2 (1+\nu)} \f]
        """
        E = double_or_blank(card, 2, 'E')
        G = double_or_blank(card, 3, 'G')
        nu = double_or_blank(card, 4, 'nu')

        if G is None and E is None:  # no E,G
            raise RuntimeError('G=%s E=%s cannot both be None' % (G, E))
        elif E is not None and G is not None and nu is not None:
            pass
        elif E is not None and nu is not None:
            G = E / 2. / (1 + nu)
        elif G is not None and nu is not None:
            E = 2 * (1 + nu) * G
        elif G is not None and E is not None:
            nu = E / (2 * G) - 1.
        elif G is None and nu is None:
            G = 0.0
            nu = 0.0
        elif E is None and nu is None:
            E = 0.0
            nu = 0.0
        else:
            msg = 'G=%s E=%s nu=%s' % (G, E, nu)
            raise RuntimeError(msg)
        self.e = E
        self.g = G
        self.nu = nu

    def _write_calculix(self, elementSet='ELSetDummyMat'):
        # default value - same properties for all values
        temperature = self.TRef
        msg = '*ELASTIC,TYPE=ISO,ELSET=%s\n' % (elementSet)
        msg += '** E,NU,TEMPERATURE\n'
        msg += '%s,%s,%s\n' % (self.e, self.nu, temperature)

        if self.rho > 0.:
            msg += '*DENSITY\n'
            msg += '%s\n' % (self.rho)
        if self.a > 0:
            msg += '*EXPANSION,TYPE=ISO,ZERO=%s\n' % (self.TRef)
            msg += '** ALPHA,ALPHA*TREF\n'
            msg += '%s,%s\n\n' % (self.a, self.a * self.TRef)
        return msg

    def write_code_aster(self):
        msg = 'M%s = DEFI_MATRIAU(ELAS=_F(E=%g, # MAT1 mid=%s\n' % (
            self.mid, self.e, self.mid)
        #msg  = 'M%s = DEFI_MATRIAU(ELAS=_F( # MAT1\n' %(self.mid)
        #msg += '                       E  =%g,\n'  %(self.e)
        msg += '                       NU =%g,\n' % (self.nu)
        msg += '                       RHO=%g),);\n' % (self.rho)
        return msg

    def cross_reference(self, model):
        msg = ' which is required by MAT1 mid=%s' % self.mid
        #self.Mcsid = model.Coord(self.Mcsid, msg=msg)  # used only for PARAM,CURVPLOT
        if self.mid in model.MATS1:
            self.mats1 = model.MATS1[self.mid]  # not using a method...
        if self.mid in model.MATT1:
            self.matt1 = model.MATT1[self.mid]  # not using a method...

    def raw_fields(self):
        list_fields = ['MAT1', self.mid, self.e, self.g, self.nu, self.rho, self.a,
                       self.TRef, self.ge, self.St, self.Sc, self.Ss, self.Mcsid]
        return list_fields

    def getG_default(self):
        if self.g == 0.0 or self.nu == 0.0:
            G = self.g
        else:
            #G_default = self.e/2./(1+self.nu)
            if self.e is None:
                G = None
            else:
                G = self.e / 2. / (1 + self.nu)
        #print("MAT1 - self.e=%s self.nu=%s self.g=%s Gdef=%s G=%s"
        #      % (self.e, self.nu,self.g, G_default, G))
        return G

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the MAT1 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        Gdefault = self.getG_default()
        G = set_blank_if_default(self.g, Gdefault)

        rho = set_blank_if_default(self.rho, 0.)
        a = set_blank_if_default(self.a, 0.)
        TRef = set_blank_if_default(self.TRef, 0.)
        ge = set_blank_if_default(self.ge, 0.)
        St = set_blank_if_default(self.St, 0.)
        Sc = set_blank_if_default(self.Sc, 0.)
        Ss = set_blank_if_default(self.Ss, 0.)
        Mcsid = set_blank_if_default(self.Mcsid, 0)

        list_fields = ['MAT1', self.mid, self.e, G, self.nu, rho, a, TRef, ge,
                       St, Sc, Ss, Mcsid]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class MAT2(AnisotropicMaterial):
    """
    Defines the material properties for linear anisotropic materials for
    two-dimensional elements.

    +-----+-------+-----+-----+------+-----+------+-----+-----+
    |  1  |   2   |  3  |  4  |  5   |  6  |  7   | 8   |  9  |
    +=====+=======+=====+=====+======+=====+======+=====+=====+
    |MAT2 |  MID  | G11 | G12 | G13  | G22 | G23  | G33 | RHO |
    +-----+-------+-----+-----+------+-----+------+-----+-----+
    |     |  A1   | A2  | A3  | TREF | GE  |  ST  | SC  | SS  |
    +-----+-------+-----+-----+------+-----+------+-----+-----+
    |     | MCSID |     |     |      |     |      |     |     |
    +-----+-------+-----+-----+------+-----+------+-----+-----+
    """
    type = 'MAT2'
    _field_map = {
        1: 'mid', 2:'G11', 3:'G12', 4:'G13', 5: 'G22', 6:'G23', 7:'G33',
        8:'rho', 9:'a1', 10:'a2', 11:'a3', 12:'TRef', 13:'ge',
        14: 'St', 15:'Sc', 16:'Ss', 17:'Mcsid',
    }

    def __init__(self, card=None, data=None, comment=''):
        AnisotropicMaterial.__init__(self, card, data)
        self.matt2 = None
        if comment:
            self._comment = comment
        if card:
            self.mid = integer(card, 1, 'mid')
            self.G11 = double_or_blank(card, 2, 'G11', 0.0)
            self.G12 = double_or_blank(card, 3, 'G12', 0.0)
            self.G13 = double_or_blank(card, 4, 'G13', 0.0)
            self.G22 = double_or_blank(card, 5, 'G22', 0.0)
            self.G23 = double_or_blank(card, 6, 'G23', 0.0)
            self.G33 = double_or_blank(card, 7, 'G33', 0.0)

            self.rho = double_or_blank(card, 8, 'rho', 0.0)
            self.a1 = double_or_blank(card, 9, 'a1') # blank?
            self.a2 = double_or_blank(card, 10, 'a2') # blank?
            self.a3 = double_or_blank(card, 11, 'a3') # blank?
            self.TRef = double_or_blank(card, 12, 'TRef', 0.0)
            self.ge = double_or_blank(card, 13, 'ge', 0.0)
            self.St = double_or_blank(card, 14, 'St') # or blank?
            self.Sc = double_or_blank(card, 15, 'Sc') # or blank?
            self.Ss = double_or_blank(card, 16, 'Ss') # or blank?
            self.Mcsid = integer_or_blank(card, 17, 'Mcsid')
            assert len(card) <= 18, 'len(MAT2 card) = %i' % len(card)
        else:
            self.mid = data[0]
            self.G11 = data[1]
            self.G12 = data[2]
            self.G13 = data[3]
            self.G22 = data[4]
            self.G23 = data[5]
            self.G33 = data[6]

            self.rho = data[7]
            self.a1 = data[8]
            self.a2 = data[9]
            self.a3 = data[10]
            self.TRef = data[11]
            self.ge = data[12]
            self.St = data[13]
            self.Sc = data[14]
            self.Ss = data[15]
            self.Mcsid = data[16]

    def get_density(self):
        return self.rho

    def cross_reference(self, model):
        msg = ' which is required by MAT2 mid=%s' % self.mid
        if self.mid in model.MATT2:
            self.matt2 = model.MATT2[self.mid]  # not using a method...

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the MAT2 object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        pass

    def Dsolid(self):
        """
        Eq 9.4.7 in Finite Element Method using Matlab
        """
        D = zeros((6, 6))
        E = self.E()
        nu12 = self.nu12
        nu = nu12

        mu = 1. - nu12 * nu12  # *E11/E22    # not necessary b/c they're equal
        Emu = E / mu
        D[0, 0] = Emu  # E/(1-nu^2)
        D[1, 1] = Emu
        D[2, 2] = Emu
        D[0, 1] = nu * Emu  # nu*E/(1-nu^2)

        # nu*E/(1-nu^2)
        D[1, 2] = D[2, 1] = D[0, 2] = D[2, 0] = D[1, 0] = D[0, 1]

        # (1.-nu)/2.*E/(1-nu^2)
        D[3, 3] = (1. - nu) * 0.5 * Emu

        # (1.-nu)/2.*E/(1-nu^2)
        D[5, 5] = D[4, 4] = D[3, 3]

    def Dplate(self):
        """
        Eq 9.1.6 in Finite Element Method using Matlab
        """
        E = self.E()
        nu12 = self.Nu()
        nu = nu12
        #G12 = self.G()

        D = zeros((3, 3))
        mu = 1. - nu12 * nu12  # *E11/E22    # not necessary b/c they're equal
        Emu = E / mu
        D[0, 0] = Emu
        D[1, 1] = Emu
        D[0, 1] = nu * Emu
        D[1, 0] = D[0, 1]
        D[2, 2] = 1. - nu / 2. * Emu
        #D[4,4] =      #: .. todo:: verify
        #D[5,5] = G22
        #D[6,6] = G33
        return D

    def write_calculix(self):
        raise NotImplementedError(self.type)
        #msg = '*ELASTIC,TYPE=ORTHO\n'
        #temperature = 0.  # default value - same properties for all values
        #msg += '%s,%s,%s\n' % (self.e, self.nu, temperature)
        #D = Dplate
        #D1111 = D[0, 0]
        #D1122 = 0.
        #D2222 = D[1, 1]
        #D1133 = D[0, 2]
        #D2233 = D[1, 2]
        #D3333 = D[2, 2]
        #D1212 = D[0, 1]
        #D1313 = D[0, 2]
        #msg += '%s,%s,%s,%s,%s,%s,%s,%s\n\n' % (
            #D1111, D1122, D2222, D1133, D2233, D3333, D1212, D1313)

        ##G23
        #temperature = self.TRef
        #msg = '*ELASTIC,TYPE=ENGINEERING CONSTANTS  ** MAT2,mid=%s\n' % (
            #self.mid)
        #msg += '** E1,E2,E3,NU12,NU13,NU23,G12,G13\n'
        #msg += '** G23,TEMPERATURE\n'
        #msg += '%s,%s,%s,%s,%s,%s,%s,%s\n' % (
            #e1, e2, e3, nu12, nu13, nu23, g12, g13)
        #msg += '%s,%s\n' % (G23, temperature)
        #if self.rho > 0.:
            #msg += '*DENSITY\n'
            #msg += '%s\n' % (self.rho)
        #if self.a > 0:
            #msg += '*EXPANSION,TYPE=ISO,ZERO=%s\n' % (self.TRef)
            #msg += '** ALPHA,ALPHA*TREF\n'
            #msg += '%s,%s\n\n' % (self.a, self.a * self.TRef)
        #return msg

    def raw_fields(self):
        list_fields = ['MAT2', self.mid, self.G11, self.G12, self.G13, self.G22,
                       self.G23, self.G33, self.rho, self.a1, self.a2, self.a3,
                       self.TRef, self.ge, self.St, self.Sc, self.Ss,
                       self.Mcsid]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the MAT2 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        G11 = set_blank_if_default(self.G11, 0.0)
        G12 = set_blank_if_default(self.G12, 0.0)
        G13 = set_blank_if_default(self.G13, 0.0)
        G22 = set_blank_if_default(self.G22, 0.0)
        G23 = set_blank_if_default(self.G23, 0.0)
        G33 = set_blank_if_default(self.G33, 0.0)
        rho = set_blank_if_default(self.rho, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        list_fields = ['MAT2', self.mid, G11, G12, G13, G22, G23, G33, rho,
                       self.a1, self.a2, self.a3, TRef, ge,
                       self.St, self.Sc, self.Ss, self.Mcsid]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class MAT3(OrthotropicMaterial):
    """
    Defines the material properties for linear orthotropic materials used by
    the CTRIAX6 element entry.

    +-----+-----+----+-----+----+-------+-------+------+-----+
    |  1  |  2  |  3 |  4  | 5  |   6   |   7   |  8   |  9  |
    +=====+=====+====+=====+====+=======+=======+======+=====+
    |MAT3 | MID | EX | ETH | EZ | NUXTH | NUTHZ | NUZX | RHO |
    +-----+-----+----+-----+----+-------+-------+------+-----+
    |     |     |    | GZX | AX |  ATH  |  AZ   | TREF | GE  |
    +-----+-----+----+-----+----+-------+-------+------+-----+
    """
    type = 'MAT3'
    _field_map = {
        1: 'mid', 2:'ex', 3:'eth', 4:'ez', 5: 'nuxth', 6:'nuthz', 7:'nuzx',
        8:'rho', 11:'gzx', 12:'ax', 13:'ath', 14:'az', 15:'TRef',
        16: 'ge',
    }

    def __init__(self, card=None, data=None, comment=''):
        OrthotropicMaterial.__init__(self, card, data)
        self.mats3 = None
        self.matt3 = None
        if comment:
            self._comment = comment
        if card:
            self.mid = integer(card, 1, 'mid')
            self.ex = double(card, 2, 'ex')
            self.eth = double(card, 3, 'eth')
            self.ez = double(card, 4, 'ez')
            self.nuxth = double(card, 5, 'nuxth')
            self.nuthz = double(card, 6, 'nuthz')
            self.nuzx = double(card, 7, 'nuzx')
            self.rho = double_or_blank(card, 8, 'rho', 0.0)

            self.gzx = double_or_blank(card, 11, 'gzx')
            self.ax = double_or_blank(card, 12, 'ax', 0.0)
            self.ath = double_or_blank(card, 13, 'ath', 0.0)
            self.az = double_or_blank(card, 14, 'az', 0.0)
            self.TRef = double_or_blank(card, 15, 'TRef', 0.0)
            self.ge = double_or_blank(card, 16, 'ge', 0.0)
            assert len(card) <= 17, 'len(MAT3 card) = %i' % len(card)
        else:
            self.mid = data[0]
            self.ex = data[1]
            self.eth = data[2]
            self.ez = data[3]
            self.nuxth = data[4]
            self.nuthz = data[5]
            self.nuzx = data[6]

            self.rho = data[7]
            self.gzx = data[8]
            self.ax = data[9]
            self.ath = data[10]
            self.az = data[11]
            self.TRef = data[12]
            self.ge = data[13]

    def get_density(self):
        return self.rho

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the MAT1 object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        mid = self.Mid()
        assert isinstance(mid, int), 'mid=%r' % mid
        if xref:
            if [self.mats3, self.matt3] == [None, None]:
                pass

    def cross_reference(self, model):
        #msg = ' which is required by MAT3 mid=%s' % self.mid
        if self.mid in model.MATT3:
            self.matt3 = model.MATT3[self.mid]  # not using a method...

    def raw_fields(self):
        list_fields = ['MAT3', self.mid, self.ex, self.eth, self.ez, self.nuxth,
                       self.nuthz, self.nuzx, self.rho, None, None, self.gzx,
                       self.ax, self.ath, self.az, self.TRef, self.ge]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the MAT3 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        ax = set_blank_if_default(self.ax, 0.0)
        ath = set_blank_if_default(self.ath, 0.0)
        az = set_blank_if_default(self.az, 0.0)
        rho = set_blank_if_default(self.rho, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        list_fields = ['MAT3', self.mid, self.ex, self.eth, self.ez, self.nuxth,
                       self.nuthz, self.nuzx, rho, None, None, self.gzx,
                       ax, ath, az, TRef, ge]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class MAT4(ThermalMaterial):
    """
    Defines the constant or temperature-dependent thermal material properties
    for conductivity, heat capacity, density, dynamic viscosity, heat
    generation, reference enthalpy, and latent heat associated with a
    single-phase change.

    +-----+-----+--------+------+-----+----+-----+------+---------+
    |  1  |  2  |   3    |   4  |  5  | 6  |  7  |  8   |    9    |
    +=====+=====+========+======+=====+====+=====+======+=========+
    |MAT4 | MID |   K    |  CP  | RHO | MU |  H  | HGEN | REFENTH |
    +-----+-----+--------+------+-----+----+-----+------+---------+
    |     | TCH | TDELTA | QLAT |     |    |     |      |         |
    +-----+-----+--------+------+-----+----+-----+------+---------+
    """
    type = 'MAT4'
    _field_map = {
        1: 'mid', 2:'k', 3:'cp', 4:'rho', 5: 'mu', 6:'H', 7:'hgen',
        8:'refEnthalpy', 9:'tch', 10:'tdelta', 11:'qlat',
    }

    def __init__(self, card=None, data=None, comment=''):
        ThermalMaterial.__init__(self, card, data)
        self.matt4 = None
        if comment:
            self._comment = comment
        if card:
            self.mid = integer(card, 1, 'mid')
            self.k = double_or_blank(card, 2, 'k')
            self.cp = double_or_blank(card, 3, 'cp', 0.0)
            self.rho = double_or_blank(card, 4, 'rho', 1.0)
            self.H = double_or_blank(card, 5, 'H')
            self.mu = double_or_blank(card, 6, 'mu')
            self.hgen = double_or_blank(card, 7, 'hgen', 1.0)
            self.refEnthalpy = double_or_blank(card, 8, 'refEnthalpy')
            self.tch = double_or_blank(card, 9, 'tch')
            self.tdelta = double_or_blank(card, 10, 'tdelta')
            self.qlat = double_or_blank(card, 11, 'qlat')
            assert len(card) <= 12, 'len(MAT4 card) = %i' % len(card)
        else:
            self.mid = data[0]
            self.k = data[1]
            self.cp = data[2]
            self.rho = data[3]
            self.H = data[4]
            self.mu = data[5]
            self.hgen = data[6]
            self.refEnthalpy = data[7]
            self.tch = data[8]
            self.tdelta = data[9]
            self.qlat = data[10]

    def get_density(self):
        return self.rho

    def cross_reference(self, model):
        #msg = ' which is required by MAT4 mid=%s' % self.mid
        if self.mid in model.MATT4:
            self.matt4 = model.MATT4[self.mid]  # not using a method...

    def raw_fields(self):
        list_fields = ['MAT4', self.mid, self.k, self.cp, self.rho, self.H, self.mu,
                       self.hgen, self.refEnthalpy, self.tch, self.tdelta,
                       self.qlat]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the MAT4 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        rho = set_blank_if_default(self.rho, 1.0)
        hgen = set_blank_if_default(self.hgen, 1.0)
        cp = set_blank_if_default(self.cp, 0.0)
        list_fields = ['MAT4', self.mid, self.k, cp, rho, self.H, self.mu, hgen,
                       self.refEnthalpy, self.tch, self.tdelta, self.qlat]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class MAT5(ThermalMaterial):  # also AnisotropicMaterial
    """
    Defines the thermal material properties for anisotropic materials.

    +-----+-----+-------+-----+-----+-----+-----+-----+----+
    |  1  |  2  |   3   |  4  |  5  |  6  |  7  |  8  | 9  |
    +=====+=====+=======+=====+=====+=====+=====+=====+====+
    |MAT5 | MID |  KXX  | KXY | KXZ | KYY | KYZ | KZZ | CP |
    +-----+-----+-------+-----+-----+-----+-----+-----+----+
    |     | RHO |  HGEN |     |     |     |     |     |    |
    +-----+-----+-------+-----+-----+-----+-----+-----+----+
    """
    type = 'MAT5'
    _field_map = {
        1: 'mid', 2:'kxx', 3:'kxy', 4:'kxz', 5: 'kyy', 6:'kyz', 7:'kzz',
    }

    def __init__(self, card=None, data=None, comment=''):
        ThermalMaterial.__init__(self, card, data)
        self.matt5 = None
        if comment:
            self._comment = comment
        if card:
            self.mid = integer(card, 1, 'mid')
            #: Thermal conductivity (assumed default=0.0)
            self.kxx = double_or_blank(card, 2, 'kxx', 0.0)
            self.kxy = double_or_blank(card, 3, 'kxy', 0.0)
            self.kxz = double_or_blank(card, 4, 'kxz', 0.0)
            self.kyy = double_or_blank(card, 5, 'kyy', 0.0)
            self.kyz = double_or_blank(card, 6, 'kyz', 0.0)
            self.kzz = double_or_blank(card, 7, 'kzz', 0.0)

            self.cp = double_or_blank(card, 8, 'cp', 0.0)
            self.rho = double_or_blank(card, 9, 'rho', 1.0)
            self.hgen = double_or_blank(card, 10, 'hgen', 1.0)
            assert len(card) <= 11, 'len(MAT5 card) = %i' % len(card)
        else:
            self.mid = data[0]
            self.kxx = data[1]
            self.kxy = data[2]
            self.kxz = data[3]
            self.kyy = data[4]
            self.kyz = data[5]
            self.kzz = data[6]
            self.cp = data[7]
            self.rho = data[8]
            self.hgen = data[9]


    def cross_reference(self, model):
        #msg = ' which is required by MAT5 mid=%s' % self.mid
        if self.mid in model.MATT5:
            self.matt5 = model.MATT5[self.mid]  # not using a method...

    def get_density(self):
        return self.rho

    def K(self):
        """
        thermal conductivity matrix
        """
        k = array([[self.kxx, self.kxy, self.kxz],
                   [self.kxy, self.kyy, self.kyz],
                   [self.kxz, self.kyz, self.kzz]])
        return k

    def raw_fields(self):
        list_fields = ['MAT5', self.mid, self.kxx, self.kxy, self.kxz, self.kyy,
                       self.kyz, self.kzz, self.cp, self.rho, self.hgen]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the MAT5 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        kxx = set_blank_if_default(self.kxx, 0.0)
        kyy = set_blank_if_default(self.kyy, 0.0)
        kzz = set_blank_if_default(self.kzz, 0.0)
        kxy = set_blank_if_default(self.kxy, 0.0)
        kyz = set_blank_if_default(self.kyz, 0.0)
        kxz = set_blank_if_default(self.kxz, 0.0)

        rho = set_blank_if_default(self.rho, 1.0)
        hgen = set_blank_if_default(self.hgen, 1.0)
        cp = set_blank_if_default(self.cp, 0.0)
        list_fields = ['MAT5', self.mid, kxx, kxy, kxz, kyy, kyz, kzz, cp, rho,
                       hgen]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class MAT8(OrthotropicMaterial):
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
    _field_map = {
        1: 'mid', 2:'e11', 3:'e22', 4:'nu12', 5: 'g12', 6:'g1z', 7:'g2z',
        8: 'rho', 9:'a1', 10:'a2', 11:'TRef', 12:'Xt', 13:'Xc', 14:'Yt',
        15:'Yc', 16: 'S', 17:'ge', 18:'F12', 19:'strn',
    }

    def __init__(self, card=None, data=None, comment=''):
        OrthotropicMaterial.__init__(self, card, data)
        self.mats8 = None
        self.matt8 = None
        if comment:
            self._comment = comment
        if card:
            self.mid = integer(card, 1, 'mid')
            self.e11 = double(card, 2, 'E11')    #: .. todo:: is this the correct default
            self.e22 = double(card, 3, 'E22')    #: .. todo:: is this the correct default
            self.nu12 = double(card, 4, 'nu12')  #: .. todo:: is this the correct default

            self.g12 = double_or_blank(card, 5, 'g12', 0.0)
            self.g1z = double_or_blank(card, 6, 'g1z', 1e8)
            self.g2z = double_or_blank(card, 7, 'g2z', 1e8)
            self.rho = double_or_blank(card, 8, 'rho', 0.0)
            self.a1 = double_or_blank(card, 9, 'a1', 0.0)
            self.a2 = double_or_blank(card, 10, 'a2', 0.0)
            self.TRef = double_or_blank(card, 11, 'TRef', 0.0)
            self.Xt = double_or_blank(card, 12, 'Xt', 0.0)
            self.Xc = double_or_blank(card, 13, 'Xc', self.Xt)
            self.Yt = double_or_blank(card, 14, 'Yt', 0.0)
            self.Yc = double_or_blank(card, 15, 'Yc', self.Yt)
            self.S = double_or_blank(card, 16, 'S', 0.0)
            self.ge = double_or_blank(card, 17, 'ge', 0.0)
            self.F12 = double_or_blank(card, 18, 'F12', 0.0)
            self.strn = double_or_blank(card, 19, 'strn', 0.0)
            assert len(card) <= 20, 'len(MAT8 card) = %i' % len(card)
        else:
            self.mid = data[0]
            self.e11 = data[1]
            self.e22 = data[2]
            self.nu12 = data[3]

            self.g12 = data[4]
            self.g1z = data[5]
            self.g2z = data[6]
            self.rho = data[7]
            self.a1 = data[8]
            self.a2 = data[9]
            self.TRef = data[10]
            self.Xt = data[11]
            self.Xc = data[12]
            self.Yt = data[13]
            self.Yc = data[14]
            self.S = data[15]
            self.ge = data[16]
            self.F12 = data[17]
            self.strn = data[18]

    def cross_reference(self, model):
        #msg = ' which is required by MATT8 mid=%s' % self.mid
        if self.mid in model.MATT8:
            self.matt8 = model.MATT8[self.mid]  # not using a method...

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

    def raw_fields(self):
        list_fields = ['MAT8', self.mid, self.e11, self.e22, self.nu12, self.g12,
                       self.g1z, self.g2z, self.rho, self.a1, self.a2, self.TRef,
                       self.Xt, self.Xc, self.Yt, self.Yc, self.S, self.ge,
                       self.F12, self.strn]
        return list_fields

    def get_density(self):
        return self.rho

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the MAT8 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        G12 = set_blank_if_default(self.g12, 0.)
        G1z = set_blank_if_default(self.g1z, 1e8)
        G2z = set_blank_if_default(self.g2z, 1e8)

        rho = set_blank_if_default(self.rho, 0.0)
        a1 = set_blank_if_default(self.a1, 0.0)
        a2 = set_blank_if_default(self.a2, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)

        Xt = set_blank_if_default(self.Xt, 0.)
        Yt = set_blank_if_default(self.Yt, 0.)

        Xc = set_blank_if_default(self.Xc, self.Xt)
        Yc = set_blank_if_default(self.Yc, self.Yt)

        S = set_blank_if_default(self.S, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        F12 = set_blank_if_default(self.F12, 0.0)
        strn = set_blank_if_default(self.strn, 0.0)

        list_fields = ['MAT8', self.mid, self.e11, self.e22, self.nu12, G12, G1z,
                       G2z, rho, a1, a2, TRef, Xt, Xc, Yt, Yc, S, ge, F12, strn]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)

class MAT9(AnisotropicMaterial):
    """
    Defines the material properties for linear, temperature-independent,
    anisotropic materials for solid isoparametric elements (see PSOLID entry
    description).

    +-----+-----+-----+-----+-----+-----+------+-----+-----+
    |  1  |  2  | 3   | 4   |  5  |  6  |  7   | 8   |  9  |
    +=====+=====+=====+=====+=====+=====+======+=====+=====+
    |MAT9 | MID | G11 | G12 | G13 | G14 | G15  | G16 | G22 |
    +-----+-----+-----+-----+-----+-----+------+-----+-----+
    |     | G23 | G24 | G25 | G26 | G33 | G34  | G35 | G36 |
    +-----+-----+-----+-----+-----+-----+------+-----+-----+
    |     | G44 | G45 | G46 | G55 | G56 | G66  | RHO | A1  |
    +-----+-----+-----+-----+-----+-----+------+-----+-----+
    |     | A2  | A3  | A4  | A5  | A6  | TREF | GE  |     |
    +-----+-----+-----+-----+-----+-----+------+-----+-----+
    """
    type = 'MAT9'
    _field_map = {
        1: 'mid',
    }

    def __init__(self, card=None, data=None, comment=''):
        AnisotropicMaterial.__init__(self, card, data)
        self.matt9 = None
        if comment:
            self._comment = comment
        if card:
            #: Material ID
            self.mid = integer(card, 1, 'mid')
            self.G11 = double_or_blank(card, 2, 'G11', 0.0)
            self.G12 = double_or_blank(card, 3, 'G12', 0.0)
            self.G13 = double_or_blank(card, 4, 'G13', 0.0)
            self.G14 = double_or_blank(card, 5, 'G14', 0.0)
            self.G15 = double_or_blank(card, 6, 'G15', 0.0)
            self.G16 = double_or_blank(card, 7, 'G16', 0.0)
            self.G22 = double_or_blank(card, 8, 'G22', 0.0)
            self.G23 = double_or_blank(card, 9, 'G23', 0.0)
            self.G24 = double_or_blank(card, 10, 'G24', 0.0)
            self.G25 = double_or_blank(card, 11, 'G25', 0.0)
            self.G26 = double_or_blank(card, 12, 'G26', 0.0)
            self.G33 = double_or_blank(card, 13, 'G33', 0.0)
            self.G34 = double_or_blank(card, 14, 'G34', 0.0)
            self.G35 = double_or_blank(card, 15, 'G35', 0.0)
            self.G36 = double_or_blank(card, 16, 'G36', 0.0)
            self.G44 = double_or_blank(card, 17, 'G44', 0.0)
            self.G45 = double_or_blank(card, 18, 'G45', 0.0)
            self.G46 = double_or_blank(card, 19, 'G46', 0.0)
            self.G55 = double_or_blank(card, 20, 'G55', 0.0)
            self.G56 = double_or_blank(card, 21, 'G56', 0.0)
            self.G66 = double_or_blank(card, 22, 'G66', 0.0)
            self.rho = double_or_blank(card, 23, 'rho', 0.0)
            self.A = [double_or_blank(card, 24, 'A1', 0.0),
                      double_or_blank(card, 25, 'A2', 0.0),
                      double_or_blank(card, 26, 'A3', 0.0),
                      double_or_blank(card, 27, 'A4', 0.0),
                      double_or_blank(card, 28, 'A5', 0.0),
                      double_or_blank(card, 29, 'A6', 0.0)]
            self.TRef = double_or_blank(card, 30, 'TRef', 0.0)
            self.ge = double_or_blank(card, 31, 'ge', 0.0)
            assert len(card) <= 32, 'len(MAT9 card) = %i' % len(card)
        else:
            self.mid = data[0]
            self.G11 = data[1][0]
            self.G12 = data[1][1]
            self.G13 = data[1][2]
            self.G14 = data[1][3]
            self.G15 = data[1][4]
            self.G16 = data[1][5]
            self.G22 = data[1][6]
            self.G23 = data[1][7]
            self.G24 = data[1][8]
            self.G25 = data[1][9]
            self.G26 = data[1][10]
            self.G33 = data[1][11]
            self.G34 = data[1][12]
            self.G35 = data[1][13]
            self.G36 = data[1][14]
            self.G44 = data[1][15]
            self.G45 = data[1][16]
            self.G46 = data[1][17]
            self.G55 = data[1][18]
            self.G56 = data[1][19]
            self.G66 = data[1][20]
            self.rho = data[2]
            self.A = data[3]
            self.TRef = data[4]
            self.ge = data[5]

        assert len(self.A) == 6

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the MAT9 object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        mid = self.Mid()
        #E11 = self.E11()
        #E22 = self.E22()
        #nu12 = self.Nu12()
        #G12 = self.G12()
        assert isinstance(mid, int), 'mid=%r' % mid
        #assert isinstance(E11, float), 'E11=%r' % E11
        #assert isinstance(E22, float), 'E11=%r' % E11
        #assert isinstance(G12, float), 'G12=%r' % G12
        #assert isinstance(nu12, float), 'nu12=%r' % nu12

    def D(self):
        D = array(
            [[self.G11, self.G12, self.G13, self.G14, self.G15, self.G16],
             [self.G12, self.G22, self.G23, self.G24, self.G25, self.G26],
             [self.G13, self.G23, self.G33, self.G34, self.G35, self.G36],
             [self.G14, self.G24, self.G34, self.G44, self.G45, self.G46],
             [self.G15, self.G25, self.G35, self.G45, self.G55, self.G56],
             [self.G16, self.G26, self.G36, self.G46, self.G56, self.G66]])
        return D

    def raw_fields(self):
        list_fields = (['MAT9', self.mid, self.G11, self.G12, self.G13, self.G14,
                        self.G15, self.G16, self.G22, self.G23, self.G24, self.G25,
                        self.G26, self.G33, self.G34, self.G35, self.G36, self.G44,
                        self.G45, self.G46, self.G55, self.G56, self.G66, self.rho]
                       + self.A + [self.TRef, self.ge])
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the MAT9 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        A = []
        for a in self.A:
            a = set_blank_if_default(a, 0.0)
            A.append(a)

        rho = set_blank_if_default(self.rho, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        list_fields = (['MAT9', self.mid, self.G11, self.G12, self.G13, self.G14,
                        self.G15, self.G16, self.G22, self.G23, self.G24, self.G25,
                        self.G26, self.G33, self.G34, self.G35, self.G36, self.G44,
                        self.G45, self.G46, self.G55, self.G56, self.G66, rho]
                       + A + [TRef, ge])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class MAT10(Material):
    """
    Defines material properties for fluid elements in coupled fluid-structural
    analysis.

    +------+-----+------+-----+-----+-----+-----+-----+-----+
    |  1   |  2  |  3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +======+=====+======+=====+=====+=====+=====+=====+=====+
    |MAT10 | MID | BULK | RHO |  C  | GE  |     |     |     |
    +------+-----+------+-----+-----+-----+-----+-----+-----+
    """
    type = 'MAT10'
    _field_map = {
        1: 'mid', 2:'bulk', 3:'rho', 4:'c', 5:'ge',
    }

    def __init__(self, card=None, data=None, comment=''):
        Material.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.mid = integer(card, 1, 'mid')
            self.getBulkRhoC(card)
            self.ge = double_or_blank(card, 5, 'ge', 0.0)
            assert len(card) <= 6, 'len(MAT10 card) = %i' % len(card)
        else:
            self.mid = data[0]
            self.bulk = data[1]
            self.rho = data[2]
            self.c = data[3]
            self.ge = data[4]

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the MAT10 object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        mid = self.Mid()
        bulk = self.bulk
        rho = self.rho
        c = self.c
        ge = self.ge

        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(bulk, float), 'bulk=%r' % bulk
        assert isinstance(rho, float), 'rho=%r' % rho
        assert isinstance(c, float), 'c=%r' % c
        assert isinstance(ge, float), 'ge=%r' % ge

    def getBulkRhoC(self, card):
        r"""
        .. math:: bulk = c^2 \rho
        """
        bulk = double_or_blank(card, 2, 'bulk')
        rho = double_or_blank(card, 3, 'rho')
        c = double_or_blank(card, 4, 'c')

        if c is not None:
            if rho is not None:
                bulk = c ** 2. * rho
            elif bulk is not None:
                rho = bulk / c ** 2.
            else:
                msg = 'c is the only card defined on tbe MAT10'
                raise RuntimeError(msg)
        elif bulk is not None:
            if rho is not None:
                c = (bulk / rho) ** 0.5
            else:
                msg = 'c, bulk, and rho are all undefined on tbe MAT10'
                raise RuntimeError(msg)
        else:
            msg = 'c, bulk, and rho are all undefined on tbe MAT10'
            raise RuntimeError(msg)

        self.bulk = bulk
        self.rho = rho
        self.c = c

    def raw_fields(self):
        list_fields = ['MAT10', self.mid, self.bulk, self.rho, self.c, self.ge]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the MAT10 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        ge = set_blank_if_default(self.ge, 0.0)
        list_fields = ['MAT10', self.mid, self.bulk, self.rho, self.c, ge]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class MAT11(Material):
    """
    Defines the material properties for a 3D orthotropic material for
    isoparametric solid elements.

    +------+-----+-----+-----+----+------+------+------+-----+
    |  1   |  2  |  3  |  4  |  5 |   6  |  7   |  8   |  9  |
    +======+=====+=====+=====+====+======+======+======+=====+
    |MAT11 | MID |  E1 | E2  | E3 | NU12 | NU13 | NU23 | G12 |
    +------+-----+-----+-----+----+------+------+------+-----+
    |      | G13 | G23 | RHO | A1 |  A2  |  A3  | TREF | GE  |
    +------+-----+-----+-----+----+------+------+------+-----+
    """
    type = 'MAT11'

    _field_map = {
        1: 'mid', 2:'e1', 3:'e2', 4:'e3', 5: 'nu12', 6:'nu13', 7:'nu23',
        8: 'g12', 9:'g13', 10:'g23', 11:'rho', 12:'a1', 13:'a2', 14:'a3',
        15:'TRef', 16: 'ge',
    }
    def __init__(self, card=None, data=None, comment=''):
        Material.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.mid = integer(card, 1, 'mid')
            self.e1 = double(card, 2, 'E1')
            self.e2 = double(card, 3, 'E2')
            self.e3 = double(card, 4, 'E3')

            self.nu12 = double(card, 5, 'nu12')
            self.nu13 = double(card, 6, 'nu13')
            self.nu23 = double(card, 7, 'nu23')

            self.g12 = double(card, 8, 'g12')
            self.g13 = double(card, 9, 'g13')
            self.g23 = double(card, 10, 'g23')

            self.rho = double_or_blank(card, 11, 'rho', 0.0)
            self.a1 = double_or_blank(card, 12, 'a1', 0.0)
            self.a2 = double_or_blank(card, 13, 'a2', 0.0)
            self.a3 = double_or_blank(card, 14, 'a3', 0.0)

            self.TRef = double_or_blank(card, 15, 'TRef', 0.0)
            self.ge = double_or_blank(card, 16, 'ge', 0.0)
            assert len(card) <= 17, 'len(MAT11 card) = %i' % len(card)
        else:
            self.mid = data[0]
            self.e1 = data[1]
            self.e2 = data[2]
            self.e3 = data[3]
            self.nu12 = data[4]
            self.nu13 = data[5]
            self.nu23 = data[6]
            self.g12 = data[7]
            self.g13 = data[8]
            self.g23 = data[9]
            self.rho = data[10]
            self.a1 = data[11]
            self.a2 = data[12]
            self.a3 = data[13]
            self.TRef = data[14]
            self.ge = data[15]
        msg = 'MAT11 mid=%s does not have ' % self.mid
        assert self.e1 is not None, msg + 'E1 defined'
        assert self.e2 is not None, msg + 'E2 defined'
        assert self.e3 is not None, msg + 'E3 defined'
        assert self.g12 is not None, msg + 'G12 defined'
        assert self.g13 is not None, msg + 'G13 defined'
        assert self.g23 is not None, msg + 'G23 defined'
        assert self.nu12 is not None, msg + 'NU12 defined'
        assert self.nu13 is not None, msg + 'NU13 defined'
        assert self.nu23 is not None, msg + 'NU23 defined'

    def raw_fields(self):
        list_fields = ['MAT11', self.mid, self.e1, self.e2, self.e3, self.nu12,
                       self.nu13, self.nu23, self.g12, self.g13, self.g23, self.rho, self.a1,
                       self.a2, self.a3, self.TRef, self.ge]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the MAT11 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        a1 = set_blank_if_default(self.a1, 0.0)
        a2 = set_blank_if_default(self.a2, 0.0)
        a3 = set_blank_if_default(self.a3, 0.0)

        TRef = set_blank_if_default(self.TRef, 0.0)
        rho = set_blank_if_default(self.rho, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)

        list_fields = ['MAT11', self.mid, self.e1, self.e2, self.e3, self.nu12,
                       self.nu13, self.nu23, self.g12, self.g13, self.g23, rho, a1,
                       a2, a3, TRef, ge]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class MATHP(HyperelasticMaterial):
    type = 'MATHP'

    def __init__(self, card=None, data=None, comment=''):
        HyperelasticMaterial.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.mid = integer(card, 1, 'mid')
            self.a10 = double_or_blank(card, 2, 'a10', 0.)
            self.a01 = double_or_blank(card, 3, 'a01', 0.)
            self.d1 = double_or_blank(card, 4, 'd1', (self.a10 + self.a01) * 1000)
            self.rho = double_or_blank(card, 5, 'rho', 0.)
            self.av = double_or_blank(card, 6, 'av', 0.)
            self.TRef = double_or_blank(card, 7, 'TRef', 0.)
            self.ge = double_or_blank(card, 8, 'ge', 0.)

            self.na = integer_or_blank(card, 10, 'na', 1)
            self.nd = integer_or_blank(card, 11, 'nd', 1)

            self.a20 = double_or_blank(card, 17, 'a20', 0.)
            self.a11 = double_or_blank(card, 18, 'a11', 0.)
            self.a02 = double_or_blank(card, 19, 'a02', 0.)
            self.d2 = double_or_blank(card, 20, 'd2', 0.)

            self.a30 = double_or_blank(card, 25, 'a30', 0.)
            self.a21 = double_or_blank(card, 26, 'a21', 0.)
            self.a12 = double_or_blank(card, 27, 'a12', 0.)
            self.a03 = double_or_blank(card, 28, 'a03', 0.)
            self.d3 = double_or_blank(card, 29, 'd3', 0.)

            self.a40 = double_or_blank(card, 33, 'a40', 0.)
            self.a31 = double_or_blank(card, 34, 'a31', 0.)
            self.a22 = double_or_blank(card, 35, 'a22', 0.)
            self.a13 = double_or_blank(card, 36, 'a13', 0.)
            self.a04 = double_or_blank(card, 37, 'a04', 0.)
            self.d4 = double_or_blank(card, 38, 'd4', 0.)

            self.a50 = double_or_blank(card, 41, 'a50', 0.)
            self.a41 = double_or_blank(card, 42, 'a41', 0.)
            self.a32 = double_or_blank(card, 43, 'a32', 0.)
            self.a23 = double_or_blank(card, 44, 'a23', 0.)
            self.a14 = double_or_blank(card, 45, 'a14', 0.)
            self.a05 = double_or_blank(card, 46, 'a05', 0.)
            self.d5 = double_or_blank(card, 47, 'd5', 0.)

            self.tab1 = integer_or_blank(card, 49, 'tab1')
            self.tab2 = integer_or_blank(card, 50, 'tab2')
            self.tab3 = integer_or_blank(card, 51, 'tab3')
            self.tab4 = integer_or_blank(card, 52, 'tab4')
            self.tabd = integer_or_blank(card, 56, 'tabd')
            assert len(card) <= 57, 'len(MATHP card) = %i' % len(card)
        else:
            main = data[0]
            (mid, a10, a01, d1, rho, av, alpha, tref, ge, sf, na, nd, kp,
             a20, a11, a02, d2,
             a30, a21, a12, a03, d3,
             a40, a31, a22, a13, a04, d4,
             a50, a41, a32, a23, a14, a05, d5,
             continueFlag) = main

            self.mid = mid
            self.a10 = a10
            self.a01 = a01
            self.d1 = d1
            self.rho = rho
            self.av = av
            self.TRef = tref
            self.ge = ge

            self.na = na
            self.nd = nd

            self.a20 = a20
            self.a11 = a11
            self.a02 = a02
            self.d2 = d2

            self.a30 = a30
            self.a21 = a21
            self.a12 = a12
            self.a03 = a03
            self.d3 = d3

            self.a40 = a40
            self.a31 = a31
            self.a22 = a22
            self.a13 = a13
            self.a04 = a04
            self.d4 = d4

            self.a50 = a50
            self.a41 = a41
            self.a32 = a32
            self.a23 = a23
            self.a14 = a14
            self.a05 = a05
            self.d5 = d5

            if continueFlag:
                (tab1, tab2, tab3, tab4, x1, x2, x3, tab5) = data[1]
            else:
                tab1 = None
                tab2 = None
                tab3 = None
                tab4 = None
                tab5 = None

            self.tab1 = tab1
            self.tab2 = tab2
            self.tab3 = tab3
            self.tab4 = tab4
            self.tabd = tab5

    def raw_fields(self):
        list_fields = ['MATHP', self.mid, self.a10, self.a01, self.d1, self.rho,
                       self.av, self.TRef, self.ge,
                       None, self.na, self.nd, None, None, None, None, None,
                       self.a20, self.a11, self.a02, self.d2, None, None, None,
                       None,
                       self.a30, self.a21, self.a12, self.a03, self.d3, None,
                       None, None,
                       self.a40, self.a31, self.a22, self.a13, self.a04, self.d4,
                       None, None,
                       self.a50, self.a41, self.a32, self.a23, self.a14, self.a05,
                       self.d5, None,
                       self.tab1, self.tab2, self.tab4, None, None, None, self.tabd]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the MATHP object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        av = set_blank_if_default(self.av, 0.0)
        na = set_blank_if_default(self.na, 0.0)
        nd = set_blank_if_default(self.nd, 0.0)

        a01 = set_blank_if_default(self.a01, 0.0)
        a10 = set_blank_if_default(self.a10, 0.0)
        d1 = set_blank_if_default(self.d1, 1000 * (self.a01 + self.a10))

        a20 = set_blank_if_default(self.a20, 0.0)
        a11 = set_blank_if_default(self.a11, 0.0)
        a02 = set_blank_if_default(self.a02, 0.0)
        d2 = set_blank_if_default(self.d2, 0.0)

        a30 = set_blank_if_default(self.a30, 0.0)
        a12 = set_blank_if_default(self.a12, 0.0)
        a21 = set_blank_if_default(self.a21, 0.0)
        a03 = set_blank_if_default(self.a03, 0.0)
        d3 = set_blank_if_default(self.d3, 0.0)

        a40 = set_blank_if_default(self.a40, 0.0)
        a31 = set_blank_if_default(self.a31, 0.0)
        a22 = set_blank_if_default(self.a22, 0.0)
        a13 = set_blank_if_default(self.a13, 0.0)
        a04 = set_blank_if_default(self.a04, 0.0)
        d4 = set_blank_if_default(self.d4, 0.0)

        a50 = set_blank_if_default(self.a50, 0.0)
        a41 = set_blank_if_default(self.a41, 0.0)
        a32 = set_blank_if_default(self.a32, 0.0)
        a23 = set_blank_if_default(self.a23, 0.0)
        a14 = set_blank_if_default(self.a14, 0.0)
        a05 = set_blank_if_default(self.a05, 0.0)
        d5 = set_blank_if_default(self.d5, 0.0)

        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        list_fields = ['MATHP', self.mid, a10, a01, d1, self.rho, av, TRef, ge,
                       None, na, nd, None, None, None, None, None,
                       a20, a11, a02, d2, None, None, None, None,
                       a30, a21, a12, a03, d3, None, None, None,
                       a40, a31, a22, a13, a04, d4, None, None,
                       a50, a41, a32, a23, a14, a05, d5, None,
                       self.tab1, self.tab2, self.tab3, self.tab4,
                       None, None, None, self.tabd]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class EQUIV(Material):
    type = 'EQUIV'

    def __init__(self, card=None, data=None, comment=''):
        Material.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Identification number of a MAT1, MAT2, or MAT9 entry.
            self.mid = integer(card, 1, 'mid')
            self.field2 = integer(card, 2, 'field2')
            self.field3 = integer(card, 3, 'field3')
            self.field4 = blank(card, 4, 'field4')

            self.field5 = integer(card, 5, 'field5')
            self.field6 = integer(card, 6, 'field6')
            self.field7 = integer(card, 7, 'field7')
            #[u'EQUIV', 1, 106, 306, None, 1, 106, 306]
            #[u'EQUIV', 2, 100, 104, None, 1, 0, 4]
            assert len(card) <= 8, 'len(EQUIV card)=%i card=%s' % (len(card), card)
        else:
            raise NotImplementedError(data)

    def raw_fields(self):
        list_fields = ['EQUIV', self.Mid(), self.field2, self.field3,
                       self.field4, self.field5, self.field6, self.field7]
        return list_fields
