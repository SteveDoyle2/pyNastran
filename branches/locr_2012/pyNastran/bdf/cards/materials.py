# pylint: disable=C0103,R0902,R0904,R0914,E1101,W0612,E0602
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import zeros, array

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard, Material
from pyNastran.bdf.cards.tables import Table


class IsotropicMaterial(Material):
    """Isotropic Material Class"""
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

    def __init__(self, card=None, data=None):
        Material.__init__(self, card, data)
        if card:
            self.mid = card.field(1)
            self.T0 = card.field(2, 0.0)
            self.exp = card.field(3, 1e-9)
            self.form = card.field(4)
            self.tidkp = card.field(5)
            self.tidcp = card.field(6)
            self.tidcs = card.field(7)
            self.thresh = card.field(8, 1e-5)
            self.Type = card.field(9)
            self.a = card.field(10)
            self.b = card.field(11)
            self.c = card.field(12)
            self.d = card.field(13)
            self.e = card.field(14)
            self.f = card.field(15)
            self.g = card.field(16)
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
        ###

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def Mid(self):  # links up to MAT1, MAT2, MAT9 or same mid
        if isinstance(self.mid, int):
            return self.mid
        return self.mid.mid

    def rawFields(self):
        fields = ['CREEP', self.Mid(), self.T0, self.exp, self.form, self.tidkp, self.tidcp, self.tidcs, self.thresh,
                  self.Type, self.a, self.b, self.c, self.d, self.e, self.f, self.g]
        return fields

    def reprFields(self):
        thresh = set_blank_if_default(self.thresh, 1e-5)
        exp = set_blank_if_default(self.exp, 4.1e-9)
        T0 = set_blank_if_default(self.T0, 0.0)
        fields = ['CREEP', self.Mid(), T0, exp, self.form, self.tidkp, self.tidcp, self.tidcs, thresh,
                  self.Type, self.a, self.b, self.c, self.d, self.e, self.f, self.g]
        return fields


class MAT1(Material):
    """
    Defines the material properties for linear isotropic materials.
    MAT1     1      1.03+7  3.9615+6.3      .098
    """
    type = 'MAT1'

    def __init__(self, card=None, data=None):
        Material.__init__(self, card, data)

        if card:
            self.mid = card.field(1)
            self.set_E_G_nu(card)
            self.rho = card.field(5, 0.)
            self.a = card.field(6, 0.0)
            self.TRef = card.field(7, 0.0)
            self.ge = card.field(8, 0.0)
            self.St = card.field(9, 0.0)
            self.Sc = card.field(10, 0.0)
            self.Ss = card.field(11, 0.0)
            self.Mcsid = card.field(12, 0)
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

    def set_E_G_nu(self, card):
        r"""
        \f[ \large G = \frac{E}{2 (1+\nu)} \f]
        """
        #self.E  = card.field(2)
        #self.G  = card.field(3)
        #self.nu = card.field(4)
        #return

        E = card.field(2)
        G = card.field(3)
        nu = card.field(4)

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
        #print "mid = ",self.mid
        #print "E  = ",E
        #print "G  = ",G
        #print "nu = ",nu
        #print ""

    def writeCalculix(self, elementSet='ELSetDummyMat'):
        temperature = self.TRef  # default value - same properties for all values
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

    def writeCodeAster(self):
        msg = 'M%s = DEFI_MATRIAU(ELAS=_F(E=%g, # MAT1 mid=%s\n' % (
            self.mid, self.e, self.mid)
        #msg  = 'M%s = DEFI_MATRIAU(ELAS=_F( # MAT1\n' %(self.mid)
        #msg += '                       E  =%g,\n'  %(self.e)
        msg += '                       NU =%g,\n' % (self.nu)
        msg += '                       RHO=%g),);\n' % (self.rho)
        return msg

    #def cross_reference(self, model):
        #self.Mcsid = model.Coord(self.Mcsid)  # used only for PARAM,CURVPLOT
        #pass

    def rawFields(self):
        fields = ['MAT1', self.mid, self.e, self.g, self.nu, self.rho, self.a, self.TRef, self.ge,
                  self.St, self.Sc, self.Ss, self.Mcsid]
        return fields

    def getG_default(self):
        if self.g == 0.0 or self.nu == 0.0:
            G = self.g
        else:
            #G_default = self.e/2./(1+self.nu)
            G = self.e / 2. / (1 + self.nu)
        #print "MAT1 - self.e=%s self.nu=%s self.g=%s Gdef=%s G=%s" %(self.e, self.nu,self.g, G_default, G)
        return G

    def reprFields(self):
        G = self.getG_default()

        rho = set_blank_if_default(self.rho, 0.)
        a = set_blank_if_default(self.a, 0.)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.)
        St = set_blank_if_default(self.St, 0.)
        Sc = set_blank_if_default(self.Sc, 0.)
        Ss = set_blank_if_default(self.Ss, 0.)
        Mcsid = set_blank_if_default(self.Mcsid, 0)

        fields = ['MAT1', self.mid, self.e, G, self.nu, rho, a, TRef, ge,
                  St, Sc, Ss, Mcsid]
        return fields


class MAT2(AnisotropicMaterial):
    """
    Defines the material properties for linear anisotropic materials for
    two-dimensional elements.

    MAT2 MID G11 G12 G13 G22 G23 G33 RHO
    A1 A2 A3 TREF GE ST SC SS
    MCSID
    """
    type = 'MAT2'

    def __init__(self, card=None, data=None):
        AnisotropicMaterial.__init__(self, card, data)

        if card:
            self.mid = card.field(1)
            self.G11 = card.field(2, 0.0)
            self.G12 = card.field(3, 0.0)
            self.G13 = card.field(4, 0.0)
            self.G22 = card.field(5, 0.0)
            self.G23 = card.field(6, 0.0)
            self.G33 = card.field(7, 0.0)

            self.rho = card.field(8, 0.)
            self.a1 = card.field(9)
            self.a2 = card.field(10)
            self.a3 = card.field(11)
            self.TRef = card.field(12, 0.0)
            self.ge = card.field(13, 0.0)
            self.St = card.field(14)
            self.Sc = card.field(15)
            self.Ss = card.field(16)
            self.Mcsid = card.field(17)
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
        ###

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
        G12 = self.G()

        D = zeros((3, 3))
        mu = 1. - nu12 * nu12  # *E11/E22    # not necessary b/c they're equal
        Emu = E / mu
        D[0, 0] = Emu
        D[1, 1] = Emu
        D[0, 1] = nu * Emu
        D[1, 0] = D[0, 1]
        D[2, 2] = 1. - nu / 2. * Emu
        #D[4,4] =      ## @ todo verify
        #D[5,5] = G22
        #D[6,6] = G33
        return D

    def writeCalculix(self):
        raise NotImplementedError(self.type)
        msg = '*ELASTIC,TYPE=ORTHO\n'
        temperature = 0.  # default value - same properties for all values
        msg += '%s,%s,%s\n' % (self.e, self.nu, temperature)
        D = Dplate
        D1111 = D[0, 0]
        D1122 = 0.
        D2222 = D[1, 1]
        D1133 = D[0, 2]
        D2233 = D[1, 2]
        D3333 = D[2, 2]
        D1212 = D[0, 1]
        D1313 = D[0, 2]
        msg += '%s,%s,%s,%s,%s,%s,%s,%s\n\n' % (
            D1111, D1122, D2222, D1133, D2233, D3333, D1212, D1313)

        G23
        temperature = self.TRef
        msg = '*ELASTIC,TYPE=ENGINEERING CONSTANTS  ** MAT2,mid=%s\n' % (
            self.mid)
        msg += '** E1,E2,E3,NU12,NU13,NU23,G12,G13\n'
        msg += '** G23,TEMPERATURE\n'
        msg += '%s,%s,%s,%s,%s,%s,%s,%s\n' % (
            e1, e2, e3, nu12, nu13, nu23, g12, g13)
        msg += '%s,%s\n' % (G23, temperature)
        if self.rho > 0.:
            msg += '*DENSITY\n'
            msg += '%s\n' % (self.rho)
        if self.a > 0:
            msg += '*EXPANSION,TYPE=ISO,ZERO=%s\n' % (self.TRef)
            msg += '** ALPHA,ALPHA*TREF\n'
            msg += '%s,%s\n\n' % (self.a, self.a * self.TRef)
        return msg

    def rawFields(self):
        fields = ['MAT2', self.mid, self.G11, self.G12, self.G13, self.G22, self.G23, self.G33, self.rho,
                  self.a1, self.a2, self.a3, self.TRef, self.ge, self.St, self.Sc, self.Ss,
                  self.Mcsid]
        return fields

    def reprFields(self):
        G11 = set_blank_if_default(self.G11, 0.0)
        G12 = set_blank_if_default(self.G12, 0.0)
        G13 = set_blank_if_default(self.G13, 0.0)
        G22 = set_blank_if_default(self.G22, 0.0)
        G23 = set_blank_if_default(self.G23, 0.0)
        G33 = set_blank_if_default(self.G33, 0.0)
        rho = set_blank_if_default(self.rho, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        fields = ['MAT2', self.mid, G11, G12, G13, G22, G23, G33, rho,
                  self.a1, self.a2, self.a3, TRef, ge, self.St, self.Sc, self.Ss,
                  self.Mcsid]
        return fields


class MAT3(AnisotropicMaterial):
    """
    Defines the material properties for linear orthotropic materials used by the
    CTRIAX6 element entry.
    MAT3 MID EX  ETH EZ  NUXTH NUTHZ NUZX RHO
    -    -   GZX AX  ATH AZ TREF GE
    """
    type = 'MAT3'

    def __init__(self, card=None, data=None):
        AnisotropicMaterial.__init__(self, card, data)
        if card:
            self.mid = card.field(1)
            self.ex = card.field(2)
            self.eth = card.field(3)
            self.ez = card.field(4)
            self.nuxth = card.field(5)
            self.nuthz = card.field(6)
            self.nuzx = card.field(7)
            self.rho = card.field(8, 0.)

            self.gzx = card.field(11)
            self.ax = card.field(12)
            self.ath = card.field(13)
            self.az = card.field(14)
            self.TRef = card.field(15, 0.0)
            self.ge = card.field(16, 0.0)
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

    def rawFields(self):
        fields = ['MAT3', self.mid, self.ex, self.eth, self.ez, self.nuxth, self.nuthz, self.nuzx, self.rho,
                  None, None, self.gzx, self.ax, self.ath, self.az, self.TRef, self.ge]
        return fields

    def reprFields(self):
        rho = set_blank_if_default(self.rho, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        fields = ['MAT3', self.mid, self.ex, self.eth, self.ez, self.nuxth, self.nuthz, self.nuzx, rho,
                  None, None, self.gzx, self.ax, self.ath, self.az, TRef, ge]
        return fields


class MAT4(ThermalMaterial):
    """
    Defines the constant or temperature-dependent thermal material properties
    for conductivity, heat capacity, density, dynamic viscosity, heat
    generation, reference enthalpy, and latent heat associated with a
    single-phase change.

    MAT4 MID K CP H HGEN REFENTH
    TCH TDELTA QLAT
    """
    type = 'MAT4'

    def __init__(self, card=None, data=None):
        ThermalMaterial.__init__(self, card, data)
        if card:
            self.mid = card.field(1)
            self.k = card.field(2)
            self.cp = card.field(3, 0.0)
            self.rho = card.field(4, 1.0)
            self.H = card.field(5)
            self.mu = card.field(6)
            self.hgen = card.field(7, 1.0)
            self.refEnthalpy = card.field(8)
            self.tch = card.field(9)
            self.tdelta = card.field(10)
            self.qlat = card.field(11)
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
        ###

    def rawFields(self):
        fields = ['MAT4', self.mid, self.k, self.cp, self.rho, self.H, self.mu, self.hgen, self.refEnthalpy,
                          self.tch, self.tdelta, self.qlat]
        return fields

    def reprFields(self):
        rho = set_blank_if_default(self.rho, 1.0)
        hgen = set_blank_if_default(self.hgen, 1.0)
        cp = set_blank_if_default(self.cp, 0.0)
        fields = ['MAT4', self.mid, self.k, cp, rho, self.H, self.mu, hgen, self.refEnthalpy,
                          self.tch, self.tdelta, self.qlat]
        return fields


class MAT5(ThermalMaterial):  # also AnisotropicMaterial
    """
    Defines the thermal material properties for anisotropic materials.

    MAT5 MID KXX KXY KXZ KYY KYZ KZZ CP
    RHO HGEN
    """
    type = 'MAT5'

    def __init__(self, card=None, data=None):
        ThermalMaterial.__init__(self, card, data)
        if card:
            self.mid = card.field(1)
            ## Thermal conductivity (assumed default=0.0)
            self.kxx = card.field(2, 0.)
            self.kxy = card.field(3, 0.)
            self.kxz = card.field(4, 0.)
            self.kyy = card.field(5, 0.)
            self.kyz = card.field(6, 0.)
            self.kzz = card.field(7, 0.)

            self.cp = card.field(8, 0.0)
            self.rho = card.field(9, 1.0)
            self.hgen = card.field(10, 1.0)
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
        ###

    def K(self):
        """
        thermal conductivity matrix
        """
        k = array([[self.kxx, self.kxy, self.kxz],
                   [self.kxy, self.kyy, self.kyz],
                   [self.kxz, self.kyz, self.kzz]])
        return k

    def rawFields(self):
        fields = ['MAT5', self.mid, self.kxx, self.kxy, self.kxz, self.kyy, self.kyz, self.kzz, self.cp,
                          self.rho, self.hgen]
        return fields

    def reprFields(self):
        kxx = set_blank_if_default(self.kxx, 0.0)
        kyy = set_blank_if_default(self.kyy, 0.0)
        kzz = set_blank_if_default(self.kzz, 0.0)
        kxy = set_blank_if_default(self.kxy, 0.0)
        kyz = set_blank_if_default(self.kyz, 0.0)
        kxz = set_blank_if_default(self.kxz, 0.0)

        rho = set_blank_if_default(self.rho, 1.0)
        hgen = set_blank_if_default(self.hgen, 1.0)
        cp = set_blank_if_default(self.cp, 0.0)
        fields = ['MAT5', self.mid, kxx, kxy, kxz, kyy, kyz, kzz, cp,
                  rho, hgen]
        return fields


class MAT8(AnisotropicMaterial):
    """
    Defines the material property for an orthotropic material for isoparametric
    shell elements.
    MAT8          10  1.25+7  9.75+6     .28  1.11+7                   2.4-2
    """
    type = 'MAT8'

    def __init__(self, card=None, data=None):
        AnisotropicMaterial.__init__(self, card, data)
        if card:
            self.mid = card.field(1)
            self.e11 = card.field(2)  ## @todo is this the correct default
            self.e22 = card.field(3)  ## @todo is this the correct default
            self.nu12 = card.field(4)  ## @todo is this the correct default

            self.g12 = card.field(5, 0.0)
            self.g1z = card.field(6, 1e8)
            self.g2z = card.field(7, 1e8)
            self.rho = card.field(8, 0.0)
            self.a1 = card.field(9, 0.0)
            self.a2 = card.field(10, 0.0)
            self.TRef = card.field(11, 0.0)
            self.Xt = card.field(12, 0.0)
            self.Xc = card.field(13, self.Xt)
            self.Yt = card.field(14, 0.0)
            self.Yc = card.field(15, self.Yt)
            self.S = card.field(16, 0.0)
            self.ge = card.field(17, 0.0)
            self.F12 = card.field(18, 0.0)
            self.strn = card.field(19, 0.0)
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
        @todo what about G1z and G2z
        """
        E11 = self.E11()
        E22 = self.E22()
        nu12 = self.Nu12()
        G12 = self.G12()

        D = zeros((3, 3))
        mu = 1. - nu12 * nu12 * E11 / E22    # not necessary b/c they're equal
        D[0, 0] = E11 / mu
        D[1, 1] = E22 / mu
        D[0, 1] = nu12 * D[0, 0]
        D[1, 0] = D[0, 1]
        D[2, 2] = G12

        return D

    def rawFields(self):
        fields = ['MAT8', self.mid, self.e11, self.e22, self.nu12, self.g12, self.g1z, self.g2z, self.rho,
                  self.a1, self.a2, self.TRef, self.Xt, self.Xc, self.Yt, self.Yc, self.S,
                  self.ge, self.F12, self.strn]
        return fields

    def reprFields(self):
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

        fields = ['MAT8', self.mid, self.e11, self.e22, self.nu12, G12, G1z, G2z, rho,
                  a1, a2, TRef, Xt, Xc, Yt, Yc, S,
                  ge, F12, strn]
        return fields


class MAT9(AnisotropicMaterial):
    """
    Defines the material properties for linear, temperature-independent,
    anisotropic materials for solid isoparametric elements (see PSOLID entry
    description).

    MAT9 MID G11 G12 G13 G14 G15 G16 G22
    G23 G24 G25 G26 G33 G34 G35 G36
    G44 G45 G46 G55 G56 G66 RHO A1
    A2 A3 A4 A5 A6 TREF GE
    """
    type = 'MAT9'

    def __init__(self, card=None, data=None):
        AnisotropicMaterial.__init__(self, card, data)
        if card:
            self.mid = card.field(1)
            self.G11 = card.field(2, 0.0)
            self.G12 = card.field(3, 0.0)
            self.G13 = card.field(4, 0.0)
            self.G14 = card.field(5, 0.0)
            self.G15 = card.field(6, 0.0)
            self.G16 = card.field(7, 0.0)
            self.G22 = card.field(8, 0.0)
            self.G23 = card.field(9, 0.0)
            self.G24 = card.field(10, 0.0)
            self.G25 = card.field(11, 0.0)
            self.G26 = card.field(12, 0.0)
            self.G33 = card.field(13, 0.0)
            self.G34 = card.field(14, 0.0)
            self.G35 = card.field(15, 0.0)
            self.G36 = card.field(16, 0.0)
            self.G44 = card.field(17, 0.0)
            self.G45 = card.field(18, 0.0)
            self.G46 = card.field(19, 0.0)
            self.G55 = card.field(20, 0.0)
            self.G56 = card.field(21, 0.0)
            self.G66 = card.field(22, 0.0)
            self.rho = card.field(23, 0.0)
            self.A = card.fields(24, 30, [0.] * 6)
            self.TRef = card.field(30, 0.0)
            self.ge = card.field(31, 0.0)
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
        ###
        assert len(self.A) == 6

    def D(self):
        D = array(
            [[self.G11, self.G12, self.G13, self.G14, self.G15, self.G16],
             [self.G12, self.G22, self.G23, self.G24, self.G25, self.G26],
             [self.G13, self.G23, self.G33, self.G34, self.G35, self.G36],
             [self.G14, self.G24, self.G34, self.G44, self.G45, self.G46],
             [self.G15, self.G25, self.G35, self.G45, self.G55, self.G56],
             [self.G16, self.G26, self.G36, self.G46, self.G56, self.G66]])
        return D

    def rawFields(self):
        fields = ['MAT9', self.mid, self.G11, self.G12, self.G13, self.G14, self.G15, self.G16, self.G22,
                          self.G23, self.G24, self.G25, self.G26, self.G33, self.G34, self.G35, self.G36,
                          self.G44, self.G45, self.G46, self.G55, self.G56, self.G66, self.rho] + self.A + [self.TRef, self.ge]
        return fields

    def reprFields(self):
        A = []
        for a in self.A:
            a = set_blank_if_default(a, 0.0)
            A.append(a)

        rho = set_blank_if_default(self.rho, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        fields = ['MAT9', self.mid, self.G11, self.G12, self.G13, self.G14, self.G15, self.G16, self.G22,
                          self.G23, self.G24, self.G25, self.G26, self.G33, self.G34, self.G35, self.G36,
                          self.G44, self.G45, self.G46, self.G55, self.G56, self.G66, rho] + A + [TRef, ge]
        return fields


class MAT10(Material):
    """
    Defines material properties for fluid elements in coupled fluid-structural
    analysis.
    MAT10 MID BULK RHO C GE
    """
    type = 'MAT10'

    def __init__(self, card=None, data=None):
        Material.__init__(self, card, data)
        if card:
            self.mid = card.field(1)
            self.getBulkRhoC(card)
            self.ge = card.field(5, 0.0)
        else:
            self.mid = data[0]
            self.bulk = data[1]
            self.rho = data[2]
            self.c = data[3]
            self.ge = data[4]
        ###

    def getBulkRhoC(self, card):
        r"""
        \f[ \large bulk = c^2 \rho \f]
        """
        bulk = card.field(2)
        rho = card.field(3)
        c = card.field(4)

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
            ###
        else:
            msg = 'c, bulk, and rho are all undefined on tbe MAT10'
            raise RuntimeError(msg)
        ###
        self.bulk = bulk
        self.rho = rho
        self.c = c

    def rawFields(self):
        fields = ['MAT10', self.mid, self.bulk, self.rho, self.c, self.ge]
        return fields

    def reprFields(self):
        ge = set_blank_if_default(self.ge, 0.0)
        fields = ['MAT10', self.mid, self.bulk, self.rho, self.c, ge]
        return fields


class MATHP(HyperelasticMaterial):
    def __init__(self, card=None, data=None):
        HyperelasticMaterial.__init__(self, card, data)
        if card:
            self.mid = card.field(1)
            self.a10 = card.field(2, 0.)
            self.a01 = card.field(3, 0.)
            self.d1 = card.field(4, (self.a10 + self.a01) * 1000)
            self.rho = card.field(5, 0.)
            self.av = card.field(6, 0.)
            self.TRef = card.field(7, 0.)
            self.ge = card.field(8, 0.)

            self.na = card.field(10, 1)
            self.nd = card.field(11, 1)

            self.a20 = card.field(17, 0.)
            self.a11 = card.field(18, 0.)
            self.a02 = card.field(19, 0.)
            self.d2 = card.field(20, 0.)

            self.a30 = card.field(25, 0.)
            self.a21 = card.field(26, 0.)
            self.a12 = card.field(27, 0.)
            self.a03 = card.field(28, 0.)
            self.d3 = card.field(29, 0.)

            self.a40 = card.field(33, 0.)
            self.a31 = card.field(34, 0.)
            self.a22 = card.field(35, 0.)
            self.a13 = card.field(36, 0.)
            self.a04 = card.field(37, 0.)
            self.d4 = card.field(38, 0.)

            self.a50 = card.field(41, 0.)
            self.a41 = card.field(42)
            self.a32 = card.field(43)
            self.a23 = card.field(44)
            self.a14 = card.field(45)
            self.a05 = card.field(46)
            self.d5 = card.field(47, 0.)

            self.tab1 = card.field(49)
            self.tab2 = card.field(50)
            self.tab3 = card.field(51)
            self.tab4 = card.field(52)
            self.tabd = card.field(56)
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
            ###
            self.tab1 = tab1
            self.tab2 = tab2
            self.tab3 = tab3
            self.tab4 = tab4
            self.tabd = tab5

    def rawFields(self):
        fields = ['MATHP', self.mid, self.a10, self.a01, self.d1, self.rho, self.av, self.TRef, self.ge,
                  None, self.na, self.nd, None, None, None, None, None,
                  self.a20, self.a11, self.a02, self.d2, None, None, None, None,
                  self.a30, self.a21, self.a12, self.a03, self.d3, None, None, None,
                  self.a40, self.a31, self.a22, self.a13, self.a04, self.d4, None, None,
                  self.a50, self.a41, self.a32, self.a23, self.a14, self.a05, self.d5, None,
                  self.tab1, self.tab2, self.tab4, None, None, None, self.tabd]
        return fields

    def reprFields(self):
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
        fields = ['MATHP', self.mid, a10, a01, d1, self.rho, av, TRef, ge,
                          None, na, nd, None, None, None, None, None,
                          a20, a11, a02, d2, None, None, None, None,
                          a30, a21, a12, a03, d3, None, None, None,
                          a40, a31, a22, a13, a04, d4, None, None,
                          a50, a41, a32, a23, a14, a05, d5, None,
                          self.tab1, self.tab2, self.tab3, self.tab4, None, None, None, self.tabd]
        return fields


class MaterialDependence(BaseCard):
    def __init__(self, card, data):
        pass


class MATS1(MaterialDependence):
    """
    Specifies stress-dependent material properties for use in applications
    involving nonlinear materials. This entry is used if a MAT1, MAT2 or MAT9
    entry is specified with the same MID in a nonlinear solution sequence
    (SOLs 106 and 129).
    """
    type = 'MATS1'

    def __init__(self, card=None, data=None):
        MaterialDependence.__init__(self, card, data)
        if card:
            ## Identification number of a MAT1, MAT2, or MAT9 entry.
            self.mid = card.field(1)
            ## Identification number of a TABLES1 or TABLEST entry. If H is
            ## given, then this field must be blank.
            self.tid = card.field(2)
            ## Type of material nonlinearity. ('NLELAST' for nonlinear elastic
            ## or 'PLASTIC' for elastoplastic.)
            self.Type = card.field(3)
            ## Work hardening slope (slope of stress versus plastic strain) in
            ## units of stress. For elastic-perfectly plastic cases, H=0.0.
            ## For more than a single slope in the plastic range, the
            ## stress-strain data must be supplied on a TABLES1 entry referenced
            ## by TID, and this field must be blank
            self.h = card.field(4)
            ## Yield function criterion, selected by one of the following values
            ## (1) Von Mises (2) Tresca (3) Mohr-Coulomb (4) Drucker-Prager
            self.yf = card.field(5, 1)
            ## Hardening Rule, selected by one of the following values
            ## (Integer): (1) Isotropic (Default) (2) Kinematic
            ## (3) Combined isotropic and kinematic hardening
            self.hr = card.field(6, 1)
            ## Initial yield point
            self.limit1 = card.field(7)
            ## Internal friction angle, measured in degrees, for the
            ## Mohr-Coulomb and Drucker-Prager yield criteria
            self.limit2 = card.field(8)
        else:
            (mid, tid, Type, h, yf, hr, limit1, limit2) = data
            self.mid = mid
            self.tid = tid
            if Type == 1:
                self.Type = 'NLELAST'
            elif Type == 2:
                self.Type = 'PLASTIC'
            else:
                raise RuntimeError('invalid Type:  Type=%s; must be 1=NLELAST or 2=PLASTIC' % (Type))
            self.h = h
            self.yf = yf
            self.hr = hr
            self.limit1 = limit1
            self.limit2 = limit2
        ###

    def Yf(self):
        d = {1: 'VonMises', 2: 'Tresca', 3: 'MohrCoulomb', 4: 'Drucker-Prager'}
        return d[self.yf]

    def Hf(self):
        d = {1: 'Isotropic', 2: 'Kinematic', 3: 'Combined'}
        return d[self.hr]

    def E(self, strain=None):
        """
        Gets E (Young's Modulus) for a given strain
        @param self the object pointer
        @param strain the strain (None -> linear E value)
        @retval E (Young's Modulus)
        """
        msg = "E (Young's Modulus) not implemented for MATS1"
        raise NotImplementedError(msg)
        if self.tid:
            E = self.tid.Value(strain)
        ###
        return E

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)
        if self.tid:  # then self.h is used
            self.tid = model.Table(self.tid)
        ###

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        return self.mid.mid

    def Tid(self):
        if isinstance(self.tid, Table):
            return self.tid.tid
        return self.tid

    def rawFields(self):
        fields = ['MATS1', self.Mid(), self.Tid(), self.Type,
                  self.h, self.yf, self.hr, self.limit1, self.limit2]
        return fields

    def reprFields(self):
        return self.rawFields()
