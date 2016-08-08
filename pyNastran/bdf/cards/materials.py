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
from numpy import zeros, array

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Material
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string_or_blank, blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class IsotropicMaterial(Material):
    """Isotropic Material Class"""
    def __init__(self):
        Material.__init__(self)


class OrthotropicMaterial(Material):
    """Orthotropic Material Class"""
    def __init__(self):
        Material.__init__(self)

class AnisotropicMaterial(Material):
    """Anisotropic Material Class"""
    def __init__(self):
        Material.__init__(self)


class ThermalMaterial(Material):
    """Thermal Material Class"""
    def __init__(self):
        Material.__init__(self)


class HyperelasticMaterial(Material):
    """Hyperelastic Material Class"""
    def __init__(self):
        Material.__init__(self)


class CREEP(Material):
    type = 'CREEP'

    def __init__(self, mid, T0, exp, form, tidkp, tidcp, tidcs, thresh, Type,
                 a, b, c, d, e, f, g, comment=''):
        Material.__init__(self)
        if comment:
            self._comment = comment
        self.mid = mid
        self.T0 = T0
        self.exp = exp
        self.form = form
        self.tidkp = tidkp
        self.tidcp = tidcp
        self.tidcs = tidcs
        self.thresh = thresh
        self.Type = Type
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f
        self.g = g

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        T0 = double_or_blank(card, 2, 'T0', 0.0)
        exp = double_or_blank(card, 3, 'exp', 1e-9)
        form = string_or_blank(card, 4, 'form') # blank?
        tidkp = integer_or_blank(card, 5, 'tidkp') # blank?
        tidcp = integer_or_blank(card, 6, 'tidcp') # blank?
        tidcs = integer_or_blank(card, 7, 'tidcs') # blank?
        thresh = double_or_blank(card, 8, 'thresh', 1e-5)
        Type = integer_or_blank(card, 9, 'Type')
        # 111, 112, 121, 122, 211, 212, 221, 222, 300 (or blank?)
        a = double_or_blank(card, 10, 'a')
        b = double_or_blank(card, 11, 'b')
        c = double_or_blank(card, 12, 'c')
        d = double_or_blank(card, 13, 'd')
        e = double_or_blank(card, 14, 'e')
        f = double_or_blank(card, 15, 'f')
        g = double_or_blank(card, 16, 'g')
        assert len(card) <= 17, 'len(CREEP card) = %i\ncard=%s' % (len(card), card)
        return CREEP(mid, T0, exp, form, tidkp, tidcp, tidcs, thresh, Type,
                     a, b, c, d, e, f, g, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        mid = data[0]
        T0 = data[1]
        exp = data[2]
        form = data[3]
        tidkp = data[4]
        tidcp = data[5]
        tidcs = data[6]
        thresh = data[7]
        Type = data[8]
        a = data[9]
        b = data[10]
        c = data[11]
        d = data[12]
        e = data[13]
        f = data[14]
        g = data[15]
        return CREEP(mid, T0, exp, form, tidkp, tidcp, tidcs, thresh, Type,
                     a, b, c, d, e, f, g, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CREEP pid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid

    def uncross_reference(self):
        self.mid = self.Mid()
        del self.mid_ref

    def Mid(self):  # links up to MAT1, MAT2, MAT9 or same mid
        if isinstance(self.mid, integer_types):
            return self.mid
        return self.mid_ref.mid

    def raw_fields(self):
        list_fields = ['CREEP', self.Mid(), self.T0, self.exp, self.form,
                       self.tidkp, self.tidcp, self.tidcs, self.thresh, self.Type,
                       self.a, self.b, self.c, self.d, self.e, self.f, self.g]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : [varies, ...]
            the fields that define the card
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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


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

    def __init__(self, mid, E, G, nu,
                 rho=0.0, a=0.0, TRef=0.0, ge=0.0,
                 St=0.0, Sc=0.0, Ss=0.0, Mcsid=0, comment=''):
        IsotropicMaterial.__init__(self)
        self.mats1 = None
        self.matt1 = None
        if comment:
            self._comment = comment
        self.mid = mid
        self.e = E
        self.g = G
        self.nu = nu
        self.rho = rho
        self.a = a
        self.TRef = TRef
        self.ge = ge
        self.St = St
        self.Sc = Sc
        self.Ss = Ss
        self.Mcsid = Mcsid

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        E, G, nu = cls.set_E_G_nu(card)

        rho = double_or_blank(card, 5, 'rho', 0.)
        a = double_or_blank(card, 6, 'a', 0.0)
        TRef = double_or_blank(card, 7, 'TRef', 0.0)
        ge = double_or_blank(card, 8, 'ge', 0.0)
        St = double_or_blank(card, 9, 'St', 0.0)
        Sc = double_or_blank(card, 10, 'Sc', 0.0)
        Ss = double_or_blank(card, 11, 'Ss', 0.0)
        Mcsid = integer_or_blank(card, 12, 'Mcsid', 0)
        assert len(card) <= 13, 'len(MAT1 card) = %i\ncard=%s' % (len(card), card)
        return MAT1(mid, E, G, nu, rho, a, TRef, ge,
                    St, Sc, Ss, Mcsid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        mid = data[0]
        e = data[1]
        g = data[2]
        nu = data[3]
        rho = data[4]
        a = data[5]
        TRef = data[6]
        ge = data[7]
        St = data[8]
        Sc = data[9]
        Ss = data[10]
        Mcsid = data[11]
        return MAT1(mid, e, g, nu, rho, a, TRef, ge,
                    St, Sc, Ss, Mcsid, comment=comment)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced
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
            E = self.matt1_ref.E(self.e, temperature)
        else:
            E = self.e
        return E

    @classmethod
    def set_E_G_nu(cls, card):
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
        return E, G, nu

    def _write_calculix(self, element_set='ELSetDummyMat'):
        # default value - same properties for all values
        temperature = self.TRef
        msg = '*ELASTIC,TYPE=ISO,ELSET=%s\n' % (element_set)
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
        pre = 'M%s = DEFI_MATRIAU(ELAS=_F(' % self.mid
        spaces = ' ' * len(pre)
        msg = '%sE=%g, # MAT1 mid=%s\n' % (pre, self.e, self.mid)
        #msg  = 'M%s = DEFI_MATRIAU(ELAS=_F( # MAT1\n' %(self.mid)
        #msg += spaces + 'E  =%g,\n'  %(self.e)
        msg += spaces + 'NU=%g,\n' % (self.nu)
        if '.' in '%g' % self.rho:
            msg += spaces + 'RHO=%g));\n' % (self.rho)
        else:
            msg += spaces + 'RHO=%.1f));\n' % (self.rho)
        return msg

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by MAT1 mid=%s' % self.mid
        #self.Mcsid = model.Coord(self.Mcsid, msg=msg)  # used only for PARAM,CURVPLOT
        if self.mid in model.MATS1:
            self.mats1 = model.MATS1[self.mid]  # not using a method...
            self.mats1_ref = self.mats1
        if self.mid in model.MATT1:
            self.matt1 = model.MATT1[self.mid]  # not using a method...
            self.matt1_ref = self.matt1

    def uncross_reference(self):
        self.mats1 = self.Mats1()
        self.matt1 = self.Matt1()
        del self.mats1_ref, self.matt1_Ref

    def Mats1(self):
        return self.mats1

    def Matt1(self):
        return self.matt1

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

        Returns
        -------
        fields : [varies, ...]
            the fields that define the card
        """
        Gdefault = self.getG_default()
        G = set_blank_if_default(self.g, Gdefault)

        rho = set_blank_if_default(self.rho, 0.)
        a = set_blank_if_default(self.a, 0.)
        TRef = set_blank_if_default(self.TRef, 0.)
        ge = set_blank_if_default(self.ge, 0.)

        if [self.St, self.Sc, self.Ss, self.Mcsid] == [0., 0., 0., 0]:
            list_fields = ['MAT1', self.mid, self.e, G, self.nu, rho, a, TRef, ge]
        else:
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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


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

    def __init__(self, mid, G11, G12, G13, G22, G23, G33,
                 rho, a1, a2, a3, TRef, ge, St, Sc, Ss, Mcsid, comment=''):
        AnisotropicMaterial.__init__(self)
        self.matt2 = None
        if comment:
            self._comment = comment
        self.mid = mid
        self.G11 = G11
        self.G12 = G12
        self.G13 = G13
        self.G22 = G22
        self.G23 = G23
        self.G33 = G33
        self.rho = rho
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.TRef = TRef
        self.ge = ge
        self.St = St
        self.Sc = Sc
        self.Ss = Ss
        self.Mcsid = Mcsid

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        G11 = double_or_blank(card, 2, 'G11', 0.0)
        G12 = double_or_blank(card, 3, 'G12', 0.0)
        G13 = double_or_blank(card, 4, 'G13', 0.0)
        G22 = double_or_blank(card, 5, 'G22', 0.0)
        G23 = double_or_blank(card, 6, 'G23', 0.0)
        G33 = double_or_blank(card, 7, 'G33', 0.0)

        rho = double_or_blank(card, 8, 'rho', 0.0)
        a1 = double_or_blank(card, 9, 'a1') # blank?
        a2 = double_or_blank(card, 10, 'a2') # blank?
        a3 = double_or_blank(card, 11, 'a3') # blank?
        TRef = double_or_blank(card, 12, 'TRef', 0.0)
        ge = double_or_blank(card, 13, 'ge', 0.0)
        St = double_or_blank(card, 14, 'St') # or blank?
        Sc = double_or_blank(card, 15, 'Sc') # or blank?
        Ss = double_or_blank(card, 16, 'Ss') # or blank?
        Mcsid = integer_or_blank(card, 17, 'Mcsid')
        assert len(card) <= 18, 'len(MAT2 card) = %i\ncard=%s' % (len(card), card)
        return MAT2(mid, G11, G12, G13, G22, G23, G33,
                    rho, a1, a2, a3, TRef, ge, St, Sc, Ss, Mcsid,
                    comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        mid = data[0]
        G11 = data[1]
        G12 = data[2]
        G13 = data[3]
        G22 = data[4]
        G23 = data[5]
        G33 = data[6]

        rho = data[7]
        a1 = data[8]
        a2 = data[9]
        a3 = data[10]
        TRef = data[11]
        ge = data[12]
        St = data[13]
        Sc = data[14]
        Ss = data[15]
        Mcsid = data[16]
        return MAT2(mid, G11, G12, G13, G22, G23, G33,
                    rho, a1, a2, a3, TRef, ge, St, Sc, Ss, Mcsid,
                    comment=comment)

    def get_density(self):
        return self.rho

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by MAT2 mid=%s' % self.mid
        if self.mid in model.MATT2:
            self.matt2 = model.MATT2[self.mid]  # not using a method...
            self.matt2_ref = self.matt2

    def uncross_reference(self):
        self.matt2 = self.Matt2()
        del self.matt2_Ref

    def _verify(self, xref):
        """
        Verifies all methods for this object work

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

        Returns
        -------
        fields : [varies, ...]
            the fields that define the card
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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


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

    def __init__(self, mid, ex, eth, ez, nuxth, nuthz, nuzx, rho, gzx,
                 ax, ath, az, TRef, ge, comment=''):
        OrthotropicMaterial.__init__(self)
        self.mats3 = None
        self.matt3 = None
        if comment:
            self._comment = comment
        self.mid = mid
        self.ex = ex
        self.eth = eth
        self.ez = ez
        self.nuxth = nuxth
        self.nuthz = nuthz
        self.nuzx = nuzx
        self.rho = rho
        self.gzx = gzx
        self.ax = ax
        self.ath = ath
        self.az = az
        self.TRef = TRef
        self.ge = ge

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        ex = double(card, 2, 'ex')
        eth = double(card, 3, 'eth')
        ez = double(card, 4, 'ez')
        nuxth = double(card, 5, 'nuxth')
        nuthz = double(card, 6, 'nuthz')
        nuzx = double(card, 7, 'nuzx')
        rho = double_or_blank(card, 8, 'rho', 0.0)

        gzx = double_or_blank(card, 11, 'gzx')
        ax = double_or_blank(card, 12, 'ax', 0.0)
        ath = double_or_blank(card, 13, 'ath', 0.0)
        az = double_or_blank(card, 14, 'az', 0.0)
        TRef = double_or_blank(card, 15, 'TRef', 0.0)
        ge = double_or_blank(card, 16, 'ge', 0.0)
        assert len(card) <= 17, 'len(MAT3 card) = %i\ncard=%s' % (len(card), card)
        return MAT3(mid, ex, eth, ez, nuxth, nuthz, nuzx, rho, gzx, ax, ath, az,
                    TRef, ge, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        mid = data[0]
        ex = data[1]
        eth = data[2]
        ez = data[3]
        nuxth = data[4]
        nuthz = data[5]
        nuzx = data[6]

        rho = data[7]
        gzx = data[8]
        ax = data[9]
        ath = data[10]
        az = data[11]
        TRef = data[12]
        ge = data[13]
        return MAT3(mid, ex, eth, ez, nuxth, nuthz, nuzx, rho, gzx,
                    ax, ath, az, TRef, ge, comment=comment)

    def get_density(self):
        return self.rho

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        mid = self.Mid()
        assert isinstance(mid, int), 'mid=%r' % mid
        if xref:
            if [self.mats3, self.matt3] == [None, None]:
                pass

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        #msg = ' which is required by MAT3 mid=%s' % self.mid
        if self.mid in model.MATT3:
            self.matt3 = model.MATT3[self.mid]  # not using a method...
            self.matt3_ref = self.matt3

    def raw_fields(self):
        list_fields = ['MAT3', self.mid, self.ex, self.eth, self.ez, self.nuxth,
                       self.nuthz, self.nuzx, self.rho, None, None, self.gzx,
                       self.ax, self.ath, self.az, self.TRef, self.ge]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : [varies, ...]
            the fields that define the card
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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


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

    def __init__(self, mid, k, cp, rho, H, mu, hgen, refEnthalpy, tch, tdelta, qlat, comment=''):
        ThermalMaterial.__init__(self)
        self.matt4 = None
        if comment:
            self._comment = comment
        self.mid = mid
        self.k = k
        self.cp = cp
        self.rho = rho
        self.H = H
        self.mu = mu
        self.hgen = hgen
        self.refEnthalpy = refEnthalpy
        self.tch = tch
        self.tdelta = tdelta
        self.qlat = qlat

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        k = double_or_blank(card, 2, 'k')
        cp = double_or_blank(card, 3, 'cp', 0.0)
        rho = double_or_blank(card, 4, 'rho', 1.0)
        H = double_or_blank(card, 5, 'H')
        mu = double_or_blank(card, 6, 'mu')
        hgen = double_or_blank(card, 7, 'hgen', 1.0)
        refEnthalpy = double_or_blank(card, 8, 'refEnthalpy')
        tch = double_or_blank(card, 9, 'tch')
        tdelta = double_or_blank(card, 10, 'tdelta')
        qlat = double_or_blank(card, 11, 'qlat')
        assert len(card) <= 12, 'len(MAT4 card) = %i\ncard=%s' % (len(card), card)
        return MAT4(mid, k, cp, rho, H, mu, hgen, refEnthalpy, tch, tdelta,
                    qlat, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        mid = data[0]
        k = data[1]
        cp = data[2]
        rho = data[3]
        H = data[4]
        mu = data[5]
        hgen = data[6]
        refEnthalpy = data[7]
        tch = data[8]
        tdelta = data[9]
        qlat = data[10]
        return MAT4(mid, k, cp, rho, H, mu, hgen, refEnthalpy, tch, tdelta,
                    qlat, comment=comment)

    def get_density(self):
        return self.rho

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        #msg = ' which is required by MAT4 mid=%s' % self.mid
        if self.mid in model.MATT4:
            self.matt4 = model.MATT4[self.mid]  # not using a method...
            self.matt4_ref = self.matt4

    def raw_fields(self):
        list_fields = ['MAT4', self.mid, self.k, self.cp, self.rho, self.H, self.mu,
                       self.hgen, self.refEnthalpy, self.tch, self.tdelta,
                       self.qlat]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


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

    def __init__(self, mid, kxx=0., kxy=0., kxz=0., kyy=0., kyz=0., kzz=0.,
                 cp=0., rho=1., hgen=1., comment=''):
        ThermalMaterial.__init__(self)
        self.matt5 = None
        if comment:
            self._comment = comment
        self.mid = mid
        #: Thermal conductivity (assumed default=0.0)
        self.mid = mid
        self.kxx = kxx
        self.kxy = kxy
        self.kxz = kxz
        self.kyy = kyy
        self.kyz = kyz
        self.kzz = kzz

        self.cp = cp
        self.rho = rho
        self.hgen = hgen

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        kxx = double_or_blank(card, 2, 'kxx', 0.0)
        kxy = double_or_blank(card, 3, 'kxy', 0.0)
        kxz = double_or_blank(card, 4, 'kxz', 0.0)
        kyy = double_or_blank(card, 5, 'kyy', 0.0)
        kyz = double_or_blank(card, 6, 'kyz', 0.0)
        kzz = double_or_blank(card, 7, 'kzz', 0.0)

        cp = double_or_blank(card, 8, 'cp', 0.0)
        rho = double_or_blank(card, 9, 'rho', 1.0)
        hgen = double_or_blank(card, 10, 'hgen', 1.0)
        assert len(card) <= 11, 'len(MAT5 card) = %i\ncard=%s' % (len(card), card)
        return MAT5(mid, kxx, kxy, kxz, kyy, kyz, kzz,
                    cp, rho, hgen, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        mid = data[0]
        kxx = data[1]
        kxy = data[2]
        kxz = data[3]
        kyy = data[4]
        kyz = data[5]
        kzz = data[6]
        cp = data[7]
        rho = data[8]
        hgen = data[9]
        return MAT5(mid, kxx, kxz, kyy, kyz, kzz,
                    cp, rho, hgen, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        #msg = ' which is required by MAT5 mid=%s' % self.mid
        if self.mid in model.MATT5:
            self.matt5 = model.MATT5[self.mid]  # not using a method...
            self.matt5_ref = self.matt5

    def uncross_reference(self):
        self.matt5 = self.Matt5()
        del self.matt5_ref, self.matt5_Ref

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

        Returns
        -------
        fields : [varies, ...]
            the fields that define the card
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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


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

    def __init__(self, mid, e11, e22, nu12, g12, g1z, g2z, rho, a1, a2, TRef,
                 Xt, Xc, Yt, Yc, S, ge, F12, strn, comment=''):
        OrthotropicMaterial.__init__(self)
        if comment:
            self._comment = comment
        self.mats8 = None
        self.matt8 = None

        self.mid = mid
        self.e11 = e11
        self.e22 = e22
        # this default was tested with a complicated model (Master_model_TAXI) using NX Nastran
        # it is not defined in the QRG, but will work
        self.nu12 = nu12
        self.g12 = g12
        self.g1z = g1z
        self.g2z = g2z
        self.rho = rho
        self.a1 = a1
        self.a2 = a2
        self.TRef = TRef
        self.Xt = Xt
        self.Xc = Xc
        self.Yt = Yt
        self.Yc = Yc
        self.S = S
        self.ge = ge
        self.F12 = F12
        self.strn = strn

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        e11 = double(card, 2, 'E11')    #: .. todo:: is this the correct default
        e22 = double(card, 3, 'E22')    #: .. todo:: is this the correct default

        nu12 = double_or_blank(card, 4, 'nu12', 0.0)

        g12 = double_or_blank(card, 5, 'g12', 0.0)
        g1z = double_or_blank(card, 6, 'g1z', 1e8)
        g2z = double_or_blank(card, 7, 'g2z', 1e8)
        rho = double_or_blank(card, 8, 'rho', 0.0)
        a1 = double_or_blank(card, 9, 'a1', 0.0)
        a2 = double_or_blank(card, 10, 'a2', 0.0)
        TRef = double_or_blank(card, 11, 'TRef', 0.0)
        Xt = double_or_blank(card, 12, 'Xt', 0.0)
        Xc = double_or_blank(card, 13, 'Xc', Xt)
        Yt = double_or_blank(card, 14, 'Yt', 0.0)
        Yc = double_or_blank(card, 15, 'Yc', Yt)
        S = double_or_blank(card, 16, 'S', 0.0)
        ge = double_or_blank(card, 17, 'ge', 0.0)
        F12 = double_or_blank(card, 18, 'F12', 0.0)
        strn = double_or_blank(card, 19, 'strn', 0.0)
        assert len(card) <= 20, 'len(MAT8 card) = %i\ncard=%s' % (len(card), card)
        return MAT8(mid, e11, e22, nu12, g12, g1z, g2z, rho, a1, a2, TRef,
                    Xt, Xc, Yt, Yc, S, ge, F12, strn, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        mid = data[0]
        e11 = data[1]
        e22 = data[2]
        nu12 = data[3]

        g12 = data[4]
        g1z = data[5]
        g2z = data[6]
        rho = data[7]
        a1 = data[8]
        a2 = data[9]
        TRef = data[10]
        Xt = data[11]
        Xc = data[12]
        Yt = data[13]
        Yc = data[14]
        S = data[15]
        ge = data[16]
        F12 = data[17]
        strn = data[18]
        return MAT8(mid, e11, e22, nu12, g12, g1z, g2z, rho, a1, a2, TRef,
                    Xt, Xc, Yt, Yc, S, ge, F12, strn,
                    comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        #msg = ' which is required by MATT8 mid=%s' % self.mid
        if self.mid in model.MATT8:
            self.matt8 = model.MATT8[self.mid]  # not using a method...
            self.matt8_ref = self.matt8

    def uncross_reference(self):
        self.matt8 = self.Matt8()
        del self.matt8_ref

    def Matt8(self):
        return self.matt8

    def _verify(self, xref):
        """
        Verifies all methods for this object work

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

        Returns
        -------
        fields : [varies, ...]
            the fields that define the card
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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

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

    def __init__(self, mid, G11, G12, G13, G14, G15, G16, G22, G23, G24,
                 G25, G26, G33, G34, G35, G36, G44, G45, G46, G55, G56, G66,
                 rho, A, TRef, ge, comment=''):
        AnisotropicMaterial.__init__(self)
        self.matt9 = None
        if comment:
            self._comment = comment
        #: Material ID
        self.mid = mid
        self.G11 = G11
        self.G12 = G12
        self.G13 = G13
        self.G14 = G14
        self.G15 = G15
        self.G16 = G16
        self.G22 = G22
        self.G23 = G23
        self.G24 = G24
        self.G25 = G25
        self.G26 = G26
        self.G33 = G33
        self.G34 = G34
        self.G35 = G35
        self.G36 = G36
        self.G44 = G44
        self.G45 = G45
        self.G46 = G46
        self.G55 = G55
        self.G56 = G56
        self.G66 = G66
        self.rho = rho
        self.A = A
        self.TRef = TRef
        self.ge = ge
        assert len(self.A) == 6, A

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        G11 = double_or_blank(card, 2, 'G11', 0.0)
        G12 = double_or_blank(card, 3, 'G12', 0.0)
        G13 = double_or_blank(card, 4, 'G13', 0.0)
        G14 = double_or_blank(card, 5, 'G14', 0.0)
        G15 = double_or_blank(card, 6, 'G15', 0.0)
        G16 = double_or_blank(card, 7, 'G16', 0.0)
        G22 = double_or_blank(card, 8, 'G22', 0.0)
        G23 = double_or_blank(card, 9, 'G23', 0.0)
        G24 = double_or_blank(card, 10, 'G24', 0.0)
        G25 = double_or_blank(card, 11, 'G25', 0.0)
        G26 = double_or_blank(card, 12, 'G26', 0.0)
        G33 = double_or_blank(card, 13, 'G33', 0.0)
        G34 = double_or_blank(card, 14, 'G34', 0.0)
        G35 = double_or_blank(card, 15, 'G35', 0.0)
        G36 = double_or_blank(card, 16, 'G36', 0.0)
        G44 = double_or_blank(card, 17, 'G44', 0.0)
        G45 = double_or_blank(card, 18, 'G45', 0.0)
        G46 = double_or_blank(card, 19, 'G46', 0.0)
        G55 = double_or_blank(card, 20, 'G55', 0.0)
        G56 = double_or_blank(card, 21, 'G56', 0.0)
        G66 = double_or_blank(card, 22, 'G66', 0.0)
        rho = double_or_blank(card, 23, 'rho', 0.0)
        A = [double_or_blank(card, 24, 'A1', 0.0),
             double_or_blank(card, 25, 'A2', 0.0),
             double_or_blank(card, 26, 'A3', 0.0),
             double_or_blank(card, 27, 'A4', 0.0),
             double_or_blank(card, 28, 'A5', 0.0),
             double_or_blank(card, 29, 'A6', 0.0)]
        TRef = double_or_blank(card, 30, 'TRef', 0.0)
        ge = double_or_blank(card, 31, 'ge', 0.0)
        assert len(card) <= 32, 'len(MAT9 card) = %i\ncard=%s' % (len(card), card)
        return MAT9(mid, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25,
                    G26, G33, G34, G35, G36, G44, G45, G46,
                    G55, G56, G66, rho, A, TRef, ge,
                    comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        mid = data[0]
        G11 = data[1][0]
        G12 = data[1][1]
        G13 = data[1][2]
        G14 = data[1][3]
        G15 = data[1][4]
        G16 = data[1][5]
        G22 = data[1][6]
        G23 = data[1][7]
        G24 = data[1][8]
        G25 = data[1][9]
        G26 = data[1][10]
        G33 = data[1][11]
        G34 = data[1][12]
        G35 = data[1][13]
        G36 = data[1][14]
        G44 = data[1][15]
        G45 = data[1][16]
        G46 = data[1][17]
        G55 = data[1][18]
        G56 = data[1][19]
        G66 = data[1][20]
        for gi in data[1]:
            assert isinstance(gi, float), data[1]
        rho = data[2]
        A = data[3]
        tref = data[4]
        ge = data[5]
        assert isinstance(rho, float), rho
        assert isinstance(tref, float), tref
        assert isinstance(ge, float), ge
        for ai in A:
            assert isinstance(ai, float), A
        return MAT9(mid, G11, G12, G13, G14, G15, G16, G22, G23, G24, G25,
                    G26, G33, G34, G35, G36, G44, G45, G46,
                    G55, G56, G66, rho, A, tref, ge,
                    comment=comment)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

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

        Returns
        -------
        fields : [varies, ...]
            the fields that define the card
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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


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

    def __init__(self, mid, bulk, rho, c, ge, comment=''):
        Material.__init__(self)
        if comment:
            self._comment = comment
        self.mid = mid
        self.bulk = bulk
        self.rho = rho
        self.c = c
        self.ge = ge

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        bulk, rho, c = _mat10_get_bulk_rho_c(card)
        ge = double_or_blank(card, 5, 'ge', 0.0)
        assert len(card) <= 6, 'len(MAT10 card) = %i\ncard=%s' % (len(card), card)
        return MAT10(mid, bulk, rho, c, ge, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        mid = data[0]
        bulk = data[1]
        rho = data[2]
        c = data[3]
        ge = data[4]
        return MAT10(mid, bulk, rho, c, ge, comment=comment)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

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

    def raw_fields(self):
        list_fields = ['MAT10', self.mid, self.bulk, self.rho, self.c, self.ge]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : [varies, ...]
            the fields that define the card
        """
        ge = set_blank_if_default(self.ge, 0.0)
        list_fields = ['MAT10', self.mid, self.bulk, self.rho, self.c, ge]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

def _mat10_get_bulk_rho_c(card):
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

    return bulk, rho, c

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
    def __init__(self, mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho,
                 a1, a2, a3, TRef, ge, comment=''):
        Material.__init__(self)
        if comment:
            self._comment = comment
        self.mid = mid
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3

        self.nu12 = nu12
        self.nu13 = nu13
        self.nu23 = nu23

        self.g12 = g12
        self.g13 = g13
        self.g23 = g23

        self.rho = rho
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3

        self.TRef = TRef
        self.ge = ge
        self._validate_input()

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        e1 = double(card, 2, 'E1')
        e2 = double(card, 3, 'E2')
        e3 = double(card, 4, 'E3')

        nu12 = double(card, 5, 'nu12')
        nu13 = double(card, 6, 'nu13')
        nu23 = double(card, 7, 'nu23')

        g12 = double(card, 8, 'g12')
        g13 = double(card, 9, 'g13')
        g23 = double(card, 10, 'g23')

        rho = double_or_blank(card, 11, 'rho', 0.0)
        a1 = double_or_blank(card, 12, 'a1', 0.0)
        a2 = double_or_blank(card, 13, 'a2', 0.0)
        a3 = double_or_blank(card, 14, 'a3', 0.0)

        tref = double_or_blank(card, 15, 'TRef', 0.0)
        ge = double_or_blank(card, 16, 'ge', 0.0)
        assert len(card) <= 17, 'len(MAT11 card) = %i\ncard=%s' % (len(card), card)
        return MAT11(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho,
                     a1, a2, a3, tref, ge, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        mid = data[0]
        e1 = data[1]
        e2 = data[2]
        e3 = data[3]
        nu12 = data[4]
        nu13 = data[5]
        nu23 = data[6]
        g12 = data[7]
        g13 = data[8]
        g23 = data[9]
        rho = data[10]
        a1 = data[11]
        a2 = data[12]
        a3 = data[13]
        tref = data[14]
        ge = data[15]
        return MAT11(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, rho,
                     a1, a2, a3, tref, ge, comment=comment)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced
        """
        pass

    def _validate_input(self):
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

        Returns
        -------
        fields : [varies, ...]
            the fields that define the card
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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class MATHP(HyperelasticMaterial):
    type = 'MATHP'

    def __init__(self, mid, a10, a01, d1, rho, av, TRef, ge, na, nd,
                 a20, a11, a02, d2,
                 a30, a21, a12, a03, d3,
                 a40, a31, a22, a13, a04, d4,
                 a50, a41, a32, a23, a14, a05, d5,
                 tab1, tab2, tab3, tab4, tabd, comment=''):
        HyperelasticMaterial.__init__(self)
        if comment:
            self._comment = comment
        self.mid = mid
        self.a10 = a10
        self.a01 = a01
        self.d1 = d1
        self.rho = rho
        self.av = av
        self.TRef = TRef
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

        self.tab1 = tab1
        self.tab2 = tab2
        self.tab3 = tab3
        self.tab4 = tab4
        self.tabd = tabd

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        a10 = double_or_blank(card, 2, 'a10', 0.)
        a01 = double_or_blank(card, 3, 'a01', 0.)
        d1 = double_or_blank(card, 4, 'd1', (a10 + a01) * 1000)
        rho = double_or_blank(card, 5, 'rho', 0.)
        av = double_or_blank(card, 6, 'av', 0.)
        tref = double_or_blank(card, 7, 'TRef', 0.)
        ge = double_or_blank(card, 8, 'ge', 0.)

        na = integer_or_blank(card, 10, 'na', 1)
        nd = integer_or_blank(card, 11, 'nd', 1)

        a20 = double_or_blank(card, 17, 'a20', 0.)
        a11 = double_or_blank(card, 18, 'a11', 0.)
        a02 = double_or_blank(card, 19, 'a02', 0.)
        d2 = double_or_blank(card, 20, 'd2', 0.)

        a30 = double_or_blank(card, 25, 'a30', 0.)
        a21 = double_or_blank(card, 26, 'a21', 0.)
        a12 = double_or_blank(card, 27, 'a12', 0.)
        a03 = double_or_blank(card, 28, 'a03', 0.)
        d3 = double_or_blank(card, 29, 'd3', 0.)

        a40 = double_or_blank(card, 33, 'a40', 0.)
        a31 = double_or_blank(card, 34, 'a31', 0.)
        a22 = double_or_blank(card, 35, 'a22', 0.)
        a13 = double_or_blank(card, 36, 'a13', 0.)
        a04 = double_or_blank(card, 37, 'a04', 0.)
        d4 = double_or_blank(card, 38, 'd4', 0.)

        a50 = double_or_blank(card, 41, 'a50', 0.)
        a41 = double_or_blank(card, 42, 'a41', 0.)
        a32 = double_or_blank(card, 43, 'a32', 0.)
        a23 = double_or_blank(card, 44, 'a23', 0.)
        a14 = double_or_blank(card, 45, 'a14', 0.)
        a05 = double_or_blank(card, 46, 'a05', 0.)
        d5 = double_or_blank(card, 47, 'd5', 0.)

        tab1 = integer_or_blank(card, 49, 'tab1')
        tab2 = integer_or_blank(card, 50, 'tab2')
        tab3 = integer_or_blank(card, 51, 'tab3')
        tab4 = integer_or_blank(card, 52, 'tab4')
        tabd = integer_or_blank(card, 56, 'tabd')
        assert len(card) <= 57, 'len(MATHP card) = %i\ncard=%s' % (len(card), card)
        return MATHP(mid, a10, a01, d1, rho, av, tref, ge, na, nd, a20, a11,
                     a02, d2, a30, a21, a12, a03, d3, a40,
                     a31, a22, a13, a04, d4, a50, a41,
                     a32, a23, a14, a05, d5, tab1, tab2,
                     tab3, tab4, tabd, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        if comment:
            self._comment = comment
        main = data[0]
        (mid, a10, a01, d1, rho, av, alpha, tref, ge, sf, na, nd, kp,
         a20, a11, a02, d2,
         a30, a21, a12, a03, d3,
         a40, a31, a22, a13, a04, d4,
         a50, a41, a32, a23, a14, a05, d5,
         continue_flag) = main

        if continue_flag:
            (tab1, tab2, tab3, tab4, x1, x2, x3, tab5) = data[1]
        else:
            tab1 = None
            tab2 = None
            tab3 = None
            tab4 = None
            tab5 = None

        return MATHP(mid, a10, a01, d1, rho, av, tref, ge, na, nd, a20, a11,
                     a02, d2, a30, a21, a12, a03, d3, a40,
                     a31, a22, a13, a04, d4, a50, a41,
                     a32, a23, a14, a05, d5, tab1, tab2,
                     tab3, tab4, tabd, comment=comment)

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

        Returns
        -------
        fields : [varies, ...]
            the fields that define the card
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
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class EQUIV(Material):
    type = 'EQUIV'

    def __init__(self, mid, field2, field3, field4, field5, field6, field7, comment=''):
        Material.__init__(self)
        if comment:
            self._comment = comment
        #: Identification number of a MAT1, MAT2, or MAT9 entry.
        self.mid = mid
        self.field2 = field2
        self.field3 = field3
        self.field4 = field4
        self.field5 = field5
        self.field6 = field6
        self.field7 = field7

    @classmethod
    def add_card(cls, card, comment=''):
        mid = integer(card, 1, 'mid')
        field2 = integer(card, 2, 'field2')
        field3 = integer(card, 3, 'field3')
        field4 = blank(card, 4, 'field4')

        field5 = integer(card, 5, 'field5')
        field6 = integer(card, 6, 'field6')
        field7 = integer(card, 7, 'field7')
        #[u'EQUIV', 1, 106, 306, None, 1, 106, 306]
        #[u'EQUIV', 2, 100, 104, None, 1, 0, 4]
        assert len(card) <= 8, 'len(EQUIV card)=%i card=%s' % (len(card), card)
        return EQUIV(mid, field2, field3, field4, field5, field6, field7,
                     comment=comment)

    def raw_fields(self):
        list_fields = ['EQUIV', self.Mid(), self.field2, self.field3,
                       self.field4, self.field5, self.field6, self.field7]
        return list_fields
