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
        fields = ['CREEP',self.id,T0,exp,self.form,self.tidkp,self.tidcp,self.tidcs,thresh,
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
        #self.E  = card.field(2)
        #self.G  = card.field(3)
        #self.nu = card.field(4)
        #return

        E  = card.field(2)
        G  = card.field(3)
        nu = card.field(4)

        if G is None and E is None:
            raise RuntimeError('G=%s E=%s cannot both be None' %(G,E))
        if E and nu:
            #print "G1 = ",G
            G = E/2./(1+nu)
            #print "G2 = ",G
        elif G and nu:
            #print "E1 = ",E
            E = 2*(1+nu)*G
            #print "E2 = ",E
        elif G is None and nu is None:
            G  = 0.0
            nu = 0.0
        elif E is None and nu is None:
            E  = 0.0
            nu = 0.0
        self.E = E
        self.G = G
        self.nu = nu
        #print "E  = ",E
        #print "G  = ",G
        #print "nu = ",nu

    def __repr__(self):
        TRef = self.setBlankIfDefault(self.TRef,0.0)
        
        G_default = self.E/2./(1+self.nu)
        G    = self.setBlankIfDefault(self.G,G_default)
        #G = self.G
        fields = ['MAT1',self.mid,self.E,G,self.nu,self.rho,self.a,TRef,self.ge,
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
        fields = ['MAT2',self.mid,G11,G12,G13,G22,G23,G33,self.rho,
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
        
        fields = ['MAT8',self.mid,self.E11,self.E22,self.nu12,self.G12,G1z,G2z,self.rho,
                  a1,a2,Tref,self.Xt,Xc,self.Yt,Yc,self.S,
                  self.ge,F12,self.strn]
        #print fields
        return self.printCard(fields)

class MAT9(Material):
    """
    MAT9 MID G11 G12 G13 G14 G15 G16 G22
    G23 G24 G25 G26 G33 G34 G35 G36
    G44 G45 G46 G55 G56 G66 RHO A1
    A2 A3 A4 A5 A6 TREF GE
    """
    type = 'MAT9'
    def __init__(self,card):
        Material.__init__(self,card)
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
        self.A = card.fields(24,30,[0.]*6)
        assert len(self.A)==6
        self.TRef = card.field(30,0.0)
        self.ge   = card.field(31)

    def __repr__(self):
        A = []
        for a in self.A:
            a = self.setBlankIfDefault(a, 0.0)
            A.append(a)
        TRef = self.setBlankIfDefault(self.TRef,0.0)
        
        fields = ['MAT9',self.mid, self.G11, self.G12, self.G13, self.G14, self.G15, self.G16, self.G22,
                         self.G23, self.G24, self.G25, self.G26, self.G33, self.G34, self.G35, self.G36,
                         self.G44, self.G45, self.G46, self.G55, self.G56, self.G66, self.rho]+A+[
                                                                                 TRef,self.ge]
        return self.printCard(fields)

class MAT10(Material):
    """MAT10 MID BULK RHO C GE"""
    type = 'MAT10'
    def __init__(self,card):
        Material.__init__(self,card)
        self.mid  = card.field(1)
        (self.bulk,self.rho.self.c) = self.getBulkRhoC(card)    
        self.ge = card.field(5)

    def getBulkRhoC(self,card):
        """
        #bulk = c^2*rho
        #c^2  = bulk/rho
        #rho  = bulk/c^2
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
            
    def __repr__(self):
        fields = ['MAT10',self.mid,self.bulk,self.rho,self.c,self.ge]
        return self.printCard(fields)

