#import sys
#from numpy import array,cross,dot
#from numpy import array
#from numpy.linalg import norm

# my code
from baseCard import BaseCard
#from thermal.materials import *

class Material(BaseCard):
    """Base Material Class"""
    def __init__(self,card,data):
        pass
        #self.type = card[0]

    def isSameCard(self,mat):
        if self.type!=mat.type:  return False
        fields1 = self.rawFields()
        fields2 = mat.rawFields()
        return self.isSameFields(fields1,fields2)

    def crossReference(self,model):
        pass

    def Mid(self):
        return self.mid

class IsotropicMaterial(Material):
    """Isotropic Material Class"""
    def __init__(self,card,data):
        Material.__init__(self,card,data)

class AnisotropicMaterial(Material):
    """Anisotropic Material Class"""
    def __init__(self,card,data):
        Material.__init__(self,card,data)

class ThermalMaterial(Material):
    """Thermal Material Class"""
    def __init__(self,card,data):
        Material.__init__(self,card,data) 


class CREEP(Material):
    type = 'CREEP'
    def __init__(self,card=None,data=None):
        Material.__init__(self,card,data)
        if card:
            self.mid    = card.field(1)
            self.T0     = card.field(2,0.0)
            self.exp    = card.field(3,1e-9)
            self.form   = card.field(4)
            self.tidkp  = card.field(5)
            self.tidcp  = card.field(6)
            self.tidcs  = card.field(7)
            self.thresh = card.field(8,1e-5)
            self.Type = card.field(9)
            self.a = card.field(10)
            self.b = card.field(11)
            self.c = card.field(12)
            self.d = card.field(13)
            self.e = card.field(14)
            self.f = card.field(15)
        else:
            self.mid    = data[0]
            self.T0     = data[1]
            self.exp    = data[2]
            self.form   = data[3]
            self.tidkp  = data[4]
            self.tidcp  = data[5]
            self.tidcs  = data[6]
            self.thresh = data[7]
            self.Type   = data[8]
            self.a = data[9]
            self.b = data[10]
            self.c = data[11]
            self.d = data[12]
            self.e = data[13]
            self.f = data[14]
        ###
    
    def crossReference(self,model):
        self.mid = model.Material(self.mid)

    def Mid(self): # links up to MAT1, MAT2, MAT9 or same mid
        if isinstance(self.mid):
            return self.mid
        return self.mid.mid

    def rawFields(self):
        fields = ['CREEP',self.Mid(),self.T0,self.exp,self.form,self.tidkp,self.tidcp,self.tidcs,self.thresh,
        self.Type,self.a,self.b,self.c,self.d,self.e,self.f]
        return fields

    def reprFields(self):
        thresh = self.setBlankIfDefault(self.thresh,1e-5)
        exp    = self.setBlankIfDefault(self.exp,4.1e-9)
        T0     = self.setBlankIfDefault(self.T0,0.0)
        fields = ['CREEP',self.Mid(),T0,exp,self.form,self.tidkp,self.tidcp,self.tidcs,thresh,
        self.Type,self.a,self.b,self.c,self.d,self.e,self.f]
        return fields

class MAT1(Material):
    """
    Defines the material properties for linear isotropic materials.
    MAT1     1      1.03+7  3.9615+6.3      .098
    """
    type = 'MAT1'
    def __init__(self,card=None,data=None):
        Material.__init__(self,card,data)

        if card:
            self.mid  = card.field(1)
            self.set_E_G_nu(card)
            self.rho  = card.field(5,1e-8)
            self.a    = card.field(6,0.0)
            self.TRef = card.field(7,0.0)
            self.ge   = card.field(8,0.0)
            self.St   = card.field(9,0.0)
            self.Sc   = card.field(10,0.0)
            self.Ss   = card.field(11,0.0)
            self.Mcsid = card.field(12,0)
        else:
            self.mid   = data[0]
            self.E     = data[1]
            self.G     = data[2]
            self.nu    = data[3]
            self.rho   = data[4]
            self.a     = data[5]
            self.TRef  = data[6]
            self.ge    = data[7]
            self.St    = data[8]
            self.Sc    = data[9]
            self.Ss    = data[10]
            self.Mcsid = data[11]

    #def G(self):
    #    return self.G
    
    #def E(self):
    #    return self.E
    
    #def nu(self):
    #    return self.nu

    def set_E_G_nu(self,card):
        """
        \f[ \large G = \frac{E}{2 (1+\nu)} \f]
        """
        #self.E  = card.field(2)
        #self.G  = card.field(3)
        #self.nu = card.field(4)
        #return

        E  = card.field(2)
        G  = card.field(3)
        nu = card.field(4)

        if G is None and E is None:
            raise RuntimeError('G=%s E=%s cannot both be None' %(G,E))
        if E  is not None and nu is not None:
            G = E/2./(1+nu)
        elif G is not None and nu is not None:
            E = 2*(1+nu)*G
        elif G is not None and E is not None:
            nu = E/(2*G)-1.
        elif G is None and nu is None:
            G  = 0.0
            nu = 0.0
        elif E is None and nu is None:
            E  = 0.0
            nu = 0.0
        else:
            msg = 'G=%s E=%s nu=%s' %(G,E,nu)
            raise RuntimeError(msg)
        self.E = E
        self.G = G
        self.nu = nu
        #print "E  = ",E
        #print "G  = ",G
        #print "nu = ",nu

    def writeCodeAster(self):
        msg  = 'M%s = DEFI_MATRIAU(ELAS=_F(E=%g, # MAT1\n' %(self.mid,self.E)
        #msg  = 'M%s = DEFI_MATRIAU(ELAS=_F( # MAT1\n' %(self.mid)
        #msg += '                       E  =%g,\n'  %(self.E)
        msg += '                       NU =%g,\n'  %(self.nu)
        msg += '                       RHO=%g),);\n'  %(self.rho)
        return msg

    def rawFields(self):
        fields = ['MAT1',self.mid,self.E,self.G,self.nu,self.rho,self.a,self.TRef,self.ge,
                  self.St,self.Sc,self.Ss,self.Mcsid]
        return fields

    def reprFields(self):
        #print "MAT1 - self.E=%s self.nu=%s" %(self.E,self.nu)
        G_default = self.E/2./(1+self.nu)
        G    = self.setBlankIfDefault(self.G,G_default)

        rho  = self.setBlankIfDefault(self.rho,1e-8)
        a    = self.setBlankIfDefault(self.a,0.)
        TRef = self.setBlankIfDefault(self.TRef,0.0)
        ge   = self.setBlankIfDefault(self.ge,0.)
        St   = self.setBlankIfDefault(self.St,0.)
        Sc   = self.setBlankIfDefault(self.Sc,0.)
        Ss   = self.setBlankIfDefault(self.Ss,0.)
        Mcsid = self.setBlankIfDefault(self.Mcsid,0)

        fields = ['MAT1',self.mid,self.E,G,self.nu,rho,a,TRef,ge,
                  St,Sc,Ss,Mcsid]
        return fields

class MAT2(AnisotropicMaterial):
    """
    Defines the material properties for linear anisotropic materials for two-dimensional
    elements.

    MAT2 MID G11 G12 G13 G22 G23 G33 RHO
    A1 A2 A3 TREF GE ST SC SS
    MCSID
    """
    type = 'MAT2'
    def __init__(self,card=None,data=None):
        AnisotropicMaterial.__init__(self,card,data)
        
        if card:
            self.mid  = card.field(1)
            self.G11  = card.field(2,0.0)
            self.G12  = card.field(3,0.0)
            self.G13  = card.field(4,0.0)
            self.G22  = card.field(5,0.0)
            self.G23  = card.field(6,0.0)
            self.G33  = card.field(7,0.0)

            self.rho  = card.field(8)
            self.a1   = card.field(9)
            self.a2   = card.field(10)
            self.a3   = card.field(11)
            self.TRef = card.field(12,0.0)
            self.ge   = card.field(13)
            self.St   = card.field(14)
            self.Sc   = card.field(15)
            self.Ss   = card.field(16)
            self.Mcsid = card.field(17)
        else:
            self.mid  = data[0]
            self.G11  = data[1]
            self.G12  = data[2]
            self.G13  = data[3]
            self.G22  = data[4]
            self.G23  = data[5]
            self.G33  = data[6]

            self.rho  = data[7]
            self.a1   = data[8]
            self.a2   = data[9]
            self.a3   = data[10]
            self.TRef = data[11]
            self.ge   = data[12]
            self.St   = data[13]
            self.Sc   = data[14]
            self.Ss   = data[15]
            self.Mcsid = data[16]
        ###

    def rawFields(self):
        fields = ['MAT2',self.mid,self.G11,self.G12,self.G13,self.G22,self.G23,self.G33,self.rho,
                  self.a1,self.a2,self.a3,self.TRef,self.ge,self.St,self.Sc,self.Ss,
                  self.Mcsid]
        return fields

    def reprFields(self):
        G11 = self.setBlankIfDefault(self.G11,0.0)
        G12 = self.setBlankIfDefault(self.G12,0.0)
        G13 = self.setBlankIfDefault(self.G13,0.0)
        G22 = self.setBlankIfDefault(self.G22,0.0)
        G23 = self.setBlankIfDefault(self.G23,0.0)
        G33 = self.setBlankIfDefault(self.G33,0.0)
        TRef = self.setBlankIfDefault(self.TRef,0.0)
        fields = ['MAT2',self.mid,G11,G12,G13,G22,G23,G33,self.rho,
                  self.a1,self.a2,self.a3,TRef,self.ge,self.St,self.Sc,self.Ss,
                  self.Mcsid]
        return fields

class MAT3(AnisotropicMaterial):
    """
    Defines the material properties for linear orthotropic materials used by the CTRIAX6 element entry.
    MAT3 MID EX  ETH EZ  NUXTH NUTHZ NUZX RHO
    -    -   GZX AX  ATH AZ TREF GE
    """
    type = 'MAT3'
    def __init__(self,card=None,data=None):
        AnisotropicMaterial.__init__(self,card,data)
        if card:
            self.mid   = card.field(1)
            self.ex    = card.field(2)
            self.eth   = card.field(3)
            self.ez    = card.field(4)
            self.nuxth = card.field(5)
            self.nuthz = card.field(6)
            self.nuzx  = card.field(7)

            self.rho  = card.field(8)
            self.gzx  = card.field(9)
            self.ax   = card.field(10)
            self.ath  = card.field(11)
            self.az   = card.field(12)
            self.TRef = card.field(13,0.0)
            self.ge   = card.field(14)
        else:
            self.mid   = data[0]
            self.ex    = data[1]
            self.eth   = data[2]
            self.ez    = data[3]
            self.nuxth = data[4]
            self.nuthz = data[5]
            self.nuzx  = data[6]

            self.rho  = data[7]
            self.gzx  = data[8]
            self.ax   = data[9]
            self.ath  = data[10]
            self.az   = data[11]
            self.TRef = data[12]
            self.ge   = data[13]

    def rawFields(self):
        fields = ['MAT3',self.mid,self.ex,self.eth,self.ez,self.nuxth,self.nuthz,self.nuzx,self.rho,
                         None,None,self.gzx,self.ax,self.ath,self.az,self.TRef,self.ge]
        return fields

    def reprFields(self):
        TRef = self.setBlankIfDefault(self.TRef,0.0)
        fields = ['MAT3',self.mid,self.ex,self.eth,self.ez,self.nuxth,self.nuthz,self.nuzx,self.rho,
                         None,None,self.gzx,self.ax,self.ath,self.az,TRef,self.ge]
        return fields

class MAT4(ThermalMaterial):
    """
    Defines the constant or temperature-dependent thermal material properties for
    conductivity, heat capacity, density, dynamic viscosity, heat generation, reference
    enthalpy, and latent heat associated with a single-phase change.

    MAT4 MID K CP H HGEN REFENTH
    TCH TDELTA QLAT
    """
    type = 'MAT4'
    def __init__(self,card=None,data=None):
        ThermalMaterial.__init__(self,card,data)
        if card:
            self.mid  = card.field(1)
            self.k    = card.field(2)
            self.cp   = card.field(3,0.0)
            self.rho  = card.field(4,1.0)
            self.H    = card.field(5)
            self.mu   = card.field(6)
            self.hgen = card.field(7,1.0)
            self.refEnthalpy = card.field(8)
            self.tch    = card.field(9)
            self.tdelta = card.field(10)
            self.qlat   = card.field(11)
        else:
            self.mid  = data[0]
            self.k    = data[1]
            self.cp   = data[2]
            self.rho  = data[3]
            self.H    = data[4]
            self.mu   = data[5]
            self.hgen = data[6]
            self.refEnthalpy = data[7]
            self.tch    = data[8]
            self.tdelta = data[9]
            self.qlat   = data[10]
        ###

    def rawFields(self):
        fields = ['MAT4',self.mid,self.k,self.cp,self.rho,self.H,self.mu,self.hgen,self.refEnthalpy,
                         self.tch,self.tdelta,self.qlat]
        return fields

    def reprFields(self):
        rho  = self.setBlankIfDefault(self.rho,1.0)
        hgen = self.setBlankIfDefault(self.hgen,1.0)
        cp   = self.setBlankIfDefault(self.cp,0.0)
        fields = ['MAT4',self.mid,self.k,cp,rho,self.H,self.mu,hgen,self.refEnthalpy,
                         self.tch,self.tdelta,self.qlat]
        return fields

class MAT5(ThermalMaterial):  # also AnisotropicMaterial
    """
    Defines the thermal material properties for anisotropic materials.

    MAT5 MID KXX KXY KXZ KYY KYZ KZZ CP
    RHO HGEN
    """
    type = 'MAT5'
    def __init__(self,card=None,data=None):
        ThermalMaterial.__init__(self,card,data)
        if card:
            self.mid  = card.field(1)
            self.kxx  = card.field(2)
            self.kxy  = card.field(3)
            self.kxz  = card.field(4)
            self.kyy  = card.field(5)
            self.kyz  = card.field(6)
            self.kzz  = card.field(7)
            self.cp   = card.field(3,0.0)
            self.rho  = card.field(4,1.0)
            self.hgen = card.field(7,1.0)
        else:
            self.mid  = data[0]
            self.kxx  = data[1]
            self.kxy  = data[2]
            self.kxz  = data[3]
            self.kyy  = data[4]
            self.kyz  = data[5]
            self.kzz  = data[6]
            self.cp   = data[7]
            self.rho  = data[8]
            self.hgen = data[9]
        ###

    def rawFields(self):
        fields = ['MAT5',self.mid,self.kxx,self.kxy,self.kxz,self.kyy,self.kyz,self.kzz,self.cp,
                         self.rho,self.hgen]
        return fields

    def reprFields(self):
        rho  = self.setBlankIfDefault(self.rho, 1.0)
        hgen = self.setBlankIfDefault(self.hgen,1.0)
        cp   = self.setBlankIfDefault(self.cp,  0.0)
        fields = ['MAT5',self.mid,self.kxx,self.kxy,self.kxz,self.kyy,self.kyz,self.kzz,cp,
                         rho,hgen]
        return fields

class MAT8(AnisotropicMaterial):
    """
    Defines the material property for an orthotropic material for isoparametric shell
    elements.
    MAT8          10  1.25+7  9.75+6     .28  1.11+7                   2.4-2"""
    type = 'MAT8'
    def __init__(self,card=None,data=None):
        AnisotropicMaterial.__init__(self,card,data)
        if card:
            self.mid  = card.field(1)
            self.E11  = card.field(2)
            self.E22  = card.field(3)
            self.nu12 = card.field(4)

            self.G12  = card.field(5,0.0)
            self.G1z  = card.field(6,1e8)
            self.G2z  = card.field(7,1e8)
            self.rho  = card.field(8,0.0)
            self.a1   = card.field(9,0.0)
            self.a2   = card.field(10,0.0)
            self.TRef = card.field(11,0.0)
            self.Xt   = card.field(12,0.0)
            self.Xc   = card.field(13,self.Xt)
            self.Yt   = card.field(14,0.0)
            self.Yc   = card.field(15,self.Yt)
            self.S    = card.field(16,0.0)
            self.ge   = card.field(17,0.0)
            self.F12  = card.field(18,0.0)
            self.strn = card.field(19,0.0)
        else:
            self.mid  = data[0]
            self.E11  = data[1]
            self.E22  = data[2]
            self.nu12 = data[3]

            self.G12  = data[4]
            self.G1z  = data[5]
            self.G2z  = data[6]
            self.rho  = data[7]
            self.a1   = data[8]
            self.a2   = data[9]
            self.TRef = data[10]
            self.Xt   = data[11]
            self.Xc   = data[12]
            self.Yt   = data[13]
            self.Yc   = data[14]
            self.S    = data[15]
            self.ge   = data[16]
            self.F12  = data[17]
            self.strn = data[18]

    def rawFields(self):
        fields = ['MAT8',self.mid,self.E11,self.E22,self.nu12,self.G12,self.G1z,self.G2z,self.rho,
                  self.a1,self.a2,self.TRef,self.Xt,self.Xc,self.Yt,self.Yc,self.S,
                  self.ge,self.F12,self.strn]
        return fields

    def reprFields(self):
        G12  = self.setBlankIfDefault(self.G12,0.)
        G1z  = self.setBlankIfDefault(self.G1z,1e8)
        G2z  = self.setBlankIfDefault(self.G2z,1e8)

        rho  = self.setBlankIfDefault(self.rho, 0.0)
        a1   = self.setBlankIfDefault(self.a1,  0.0)
        a2   = self.setBlankIfDefault(self.a2,  0.0)
        TRef = self.setBlankIfDefault(self.TRef,0.0)

        Xt   = self.setBlankIfDefault(self.Xt, 0.)
        Yt   = self.setBlankIfDefault(self.Yt, 0.)

        Xc   = self.setBlankIfDefault(self.Xc, self.Xt)
        Yc   = self.setBlankIfDefault(self.Yc, self.Yt)
        
        S    = self.setBlankIfDefault(self.S,   0.0)
        ge   = self.setBlankIfDefault(self.ge,  0.0)
        F12  = self.setBlankIfDefault(self.F12, 0.0)
        strn = self.setBlankIfDefault(self.strn,0.0)
        
        fields = ['MAT8',self.mid,self.E11,self.E22,self.nu12,G12,G1z,G2z,rho,
                  a1,a2,TRef,Xt,Xc,Yt,Yc,S,
                  ge,F12,strn]
        return fields

class MAT9(AnisotropicMaterial):
    """
    Defines the material properties for linear, temperature-independent, anisotropic
    materials for solid isoparametric elements (see PSOLID entry description).

    MAT9 MID G11 G12 G13 G14 G15 G16 G22
    G23 G24 G25 G26 G33 G34 G35 G36
    G44 G45 G46 G55 G56 G66 RHO A1
    A2 A3 A4 A5 A6 TREF GE
    """
    type = 'MAT9'
    def __init__(self,card=None,data=None):
        AnisotropicMaterial.__init__(self,card,data)
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
            self.G24 = card.field(10,0.0)
            self.G25 = card.field(11,0.0)
            self.G26 = card.field(12,0.0)
            self.G33 = card.field(13,0.0)
            self.G34 = card.field(14,0.0)
            self.G35 = card.field(15,0.0)
            self.G36 = card.field(16,0.0)
            self.G44 = card.field(17,0.0)
            self.G45 = card.field(18,0.0)
            self.G46 = card.field(19,0.0)
            self.G55 = card.field(20,0.0)
            self.G56 = card.field(21,0.0)
            self.G66 = card.field(22,0.0)
            self.rho = card.field(23,0.0)
            self.A   = card.fields(24,30,[0.]*6)
            self.TRef = card.field(30,0.0)
            self.ge   = card.field(31)
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
            self.A   = data[3]
            self.TRef =data[4]
            self.ge   =data[5]
        ###
        assert len(self.A)==6
            
    def rawFields(self):
        fields = ['MAT9',self.mid, self.G11, self.G12, self.G13, self.G14, self.G15, self.G16, self.G22,
                         self.G23, self.G24, self.G25, self.G26, self.G33, self.G34, self.G35, self.G36,
                         self.G44, self.G45, self.G46, self.G55, self.G56, self.G66, self.rho]+self.A+[self.TRef,self.ge]
        return fields

    def reprFields(self):
        A = []
        for a in self.A:
            a = self.setBlankIfDefault(a, 0.0)
            A.append(a)
        TRef = self.setBlankIfDefault(self.TRef,0.0)
        ge = self.ge
        
        fields = ['MAT9',self.mid, self.G11, self.G12, self.G13, self.G14, self.G15, self.G16, self.G22,
                         self.G23, self.G24, self.G25, self.G26, self.G33, self.G34, self.G35, self.G36,
                         self.G44, self.G45, self.G46, self.G55, self.G56, self.G66, self.rho]+A+[TRef,ge]
        return fields

class MAT10(Material):
    """
    Defines material properties for fluid elements in coupled fluid-structural analysis.
    MAT10 MID BULK RHO C GE
    """
    type = 'MAT10'
    def __init__(self,card=None,data=None):
        Material.__init__(self,card,data)
        if card:
            self.mid  = card.field(1)
            (self.bulk,self.rho,self.c) = self.getBulkRhoC(card)
            self.ge = card.field(5)
        else:
            self.mid  = data[0]
            self.bulk = data[1]
            self.rho  = data[2]
            self.c    = data[3]
            self.ge   = data[4]
        ###

    def getBulkRhoC(self,card):
        """
        \f[ \large bulk = c^2 \rho \f]
        """
        bulk = card.field(2)
        rho  = card.field(3)
        c    = card.field(4)

        if c is not None:
            if rho is not None:
                bulk = c**2.*rho
            elif bulk is not None:
                rho = bulk/c**2.
            else:
                raise RuntimeError('c is the only card defined on tbe MAT10')
        elif bulk is not None:
            if rho is not None:
                c = (bulk/rho)**0.5
            else:
                raise RuntimeError('c, bulk, and rho are all undefined on tbe MAT10')
            ###
        else:
            raise RuntimeError('c, bulk, and rho are all undefined on tbe MAT10')
        ###
        return(bulk,rho,c)
            
    def rawFields(self):
        fields = ['MAT10',self.mid,self.bulk,self.rho,self.c,self.ge]
        return fields

