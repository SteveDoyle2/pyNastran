#import sys
#from numpy import array,cross,dot
#from numpy import array
#from numpy.linalg import norm

# my code
from baseCard import BaseCard

class Material(BaseCard):
    def __init__(self,card):
        #self.type = card[0]
        self.mid  = card.field(1)
        
    def __repr__(self):
        fields = [self.type,self.mid]
        return self.printCard(fields)

class CREEP(Material):
    type = 'CREEP'
    def __init__(self,card):
        Material.__init__(self,card)
        self.id   = card.field(2) # mid
        self.T0   = card.field(3,0.0)
        self.exp  = card.field(4,1e-9)
        self.form = card.field(5)
        self.tidkp  = card.field(6)
        self.tidcp  = card.field(7)
        self.tidcs  = card.field(8)
        self.thresh = card.field(9,1e-5)
        self.Type = card.field(10)
        self.a = card.field(11)
        self.b = card.field(12)
        self.c = card.field(13)
        self.d = card.field(14)
        self.e = card.field(15)
        self.f = card.field(16)

    def __repr__(self):
        thresh = self.setBlankIfDefault(self.thresh,1e-5)
        exp    = self.setBlankIfDefault(self.G1z,4.1e-9)
        T0     = self.setBlankIfDefault(self.T0,0.0)
        fields = [self.type,self.id,T0,exp,self.form,self.tidkp,self.tidcp,self.tidcs,thresh,
        self.Type,self.a,self.b,self.c,self.d,self.e,self.f]
        return self.printCard(fields)

class MAT1(Material):
    """MAT1     1      1.03+7  3.9615+6.3      .098"""
    type = 'MAT1'
    def __init__(self,card):
        Material.__init__(self,card) #mid

        self.set_E_G_nu(card)
        self.rho  = card.field(5)
        self.a    = card.field(6)
        self.TRef = card.field(7,0.0)
        self.ge = card.field(8)
        self.St = card.field(9)
        self.Sc = card.field(10)
        self.Ss = card.field(11)
        self.Mcsid = card.field(12)


    def set_E_G_nu(self,card):
        self.E  = card.field(2)
        self.G  = card.field(3)
        self.nu = card.field(4)
        #return

        E  = card.field(2)
        G  = card.field(3)
        nu = card.field(4)

        assert G is not None and E is not None,'G=%s E=%s cannot both be None' %(G,E)
        if E and nu:
            print "G1 = ",G
            G = E/2./(1+nu)
            print "G2 = ",G
        elif G and nu:
            print "E1 = ",E
            E = 2*(1+nu)*G
            print "E2 = ",E
        elif G is None and nu is None:
            G  = 0.0
            nu = 0.0
        elif E is None and nu is None:
            E  = 0.0
            nu = 0.0
        self.E = E
        self.G = G
        self.nu = nu

    def __repr__(self):
        TRef = self.setBlankIfDefault(self.TRef,0.0)
        
        G_default = self.E/2./(1+self.nu)
        G    = self.setBlankIfDefault(self.G,G_default)
        #G = self.G
        fields = [self.type,self.mid,self.E,G,self.nu,self.rho,self.a,TRef,self.ge,
                  self.St,self.Sc,self.Ss,self.Mcsid]        
        return self.printCard(fields)

class MAT2(Material):
    """
    MAT2 MID G11 G12 G13 G22 G23 G33 RHO
    A1 A2 A3 TREF GE ST SC SS
    MCSID
    """
    type = 'MAT2'
    def __init__(self,card):
        Material.__init__(self,card) # mid
        
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

    def __repr__(self):
        G11 = self.setBlankIfDefault(self.G11,0.0)
        G12 = self.setBlankIfDefault(self.G12,0.0)
        G13 = self.setBlankIfDefault(self.G13,0.0)
        G22 = self.setBlankIfDefault(self.G22,0.0)
        G23 = self.setBlankIfDefault(self.G23,0.0)
        G33 = self.setBlankIfDefault(self.G33,0.0)
        Tref = self.setBlankIfDefault(self.TRef,0.0)
        fields = [self.type,self.mid,G11,G12,G13,G22,G23,G33,self.rho,
                  self.a1,self.a2,self.a3,Tref,self.ge,self.st,self.sc,self.ss,
                  self.Mcsid]
        return self.printCard(fields)

class MAT8(Material):
    """MAT8          10  1.25+7  9.75+6     .28  1.11+7                   2.4-2"""
    type = 'MAT8'
    def __init__(self,card):
        Material.__init__(self,card)
        self.E11  = card.field(2)
        self.E22  = card.field(3)
        self.nu12 = card.field(4)

        self.G12  = card.field(5,0.0)
        self.G1z  = card.field(6,1e8)
        self.G2z  = card.field(7,1e8)
        self.rho  = card.field(8)
        self.a1   = card.field(9,0.0)
        self.a2   = card.field(10,0.0)
        self.TRef = card.field(11,0.0)
        self.Xt   = card.field(12)
        self.Xc   = card.field(13,self.Xt)
        self.Yt   = card.field(14)
        self.Yc   = card.field(15,self.Yt)
        self.S    = card.field(16)
        self.ge   = card.field(17)
        self.F12  = card.field(18,0.0)
        self.strn = card.field(19)

    def __repr__(self):
        G1z  = self.setBlankIfDefault(self.G1z,1e8)
        G2z  = self.setBlankIfDefault(self.G2z,1e8)

        a1   = self.setBlankIfDefault(self.a1, 0.0)
        a2   = self.setBlankIfDefault(self.a2, 0.0)
        Tref = self.setBlankIfDefault(self.TRef,0.0)
        Xc   = self.setBlankIfDefault(self.Xc, self.Xt)
        Yc   = self.setBlankIfDefault(self.Yc, self.Yt)
        F12  = self.setBlankIfDefault(self.F12,0.0)
        
        fields = [self.type,self.mid,self.E11,self.E22,self.nu12,self.G12,G1z,G2z,self.rho,
                  a1,a2,Tref,self.Xt,Xc,self.Yt,Yc,self.S,
                  self.ge,F12,self.strn]
        #print "hi...\n"
        #print fields
        return self.printCard(fields)
