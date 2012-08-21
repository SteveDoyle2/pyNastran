# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Property


class PFAST(Property):
    type = 'PFAST'

    def __init__(self, card=None, data=None):
        Property.__init__(self, card, data)
        ## Property ID
        self.pid = card.field(1)
        ## diameter of the fastener
        self.d = card.field(2)
        ## Specifies the element stiffness coordinate system
        self.mcid = card.field(3, -1)
        self.mflag = card.field(4, 0)  # 0-absolute 1-relative
        ## stiffness values in directions 1-3
        self.kt1 = card.field(5)
        self.kt2 = card.field(6)
        self.kt3 = card.field(7)
        ## Rotational stiffness values in directions 1-3
        self.kr1 = card.field(8, 0.0)
        self.kr2 = card.field(9, 0.0)
        self.kr3 = card.field(10, 0.0)
        ## Lumped mass of fastener
        self.mass = card.field(11, 0.0)
        ## Structural damping
        self.ge = card.field(12, 0.0)

    def cross_reference(self, model):
        self.mcid = model.Coord(self.mcid)

    def Mcid(self):
        if isinstance(self.mcid, int):
            return self.mcid
        return self.mcid.cid

    def Mass(self):
        return self.mass

    def rawFields(self):
        fields = ['PFAST', self.pid, self.d, self.Mcid(), self.mflag, self.kt1, self.kt2, self.kt3, self.kr1,
                          self.kr2, self.kr3, self.mass, self.ge]
        return fields

    def reprFields(self):
        mcid = set_blank_if_default(self.mcid, -1)
        mflag = set_blank_if_default(self.mflag, 0)
        kr1 = set_blank_if_default(self.kr1, 0.0)
        kr2 = set_blank_if_default(self.kr2, 0.0)
        kr3 = set_blank_if_default(self.kr3, 0.0)

        mass = set_blank_if_default(self.mass, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        fields = ['PFAST', self.pid, self.d, mcid, mflag, self.kt1, self.kt2, self.kt3, kr1,
                          kr2, kr3, mass, ge]
        return fields


class PGAP(Property):
    type = 'PGAP'

    def __init__(self, card=None, data=None):
        """
        Defines the properties of the gap element (CGAP entry).
        """
        Property.__init__(self, card, data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            ## initial gap opening
            self.u0 = card.field(2, 0.)
            ## preload
            self.f0 = card.field(3, 0.)
            ## axial stiffness of closed gap
            self.ka = card.field(4, 1.e8)
            ## axial stiffness of open gap
            self.kb = card.field(5, 1e-14 * self.ka)
            ## static friction coeff
            self.mu1 = card.field(7, 0.)
            ## transverse stiffness of closed gap
            self.kt = card.field(6, self.mu1 * self.ka)
            ## kinetic friction coeff
            self.mu2 = card.field(8, self.mu1)
            self.tmax = card.field(9, 0.)
            self.mar = card.field(10, 100.)
            self.trmin = card.field(11, 0.001)
        else:
            #(pid,u0,f0,ka,kb,kt,mu1,mu2,tmax,mar,trmin) = out
            self.pid = data[0]
            self.u0 = data[1]
            self.f0 = data[2]
            self.ka = data[3]
            self.kb = data[4]
            self.kt = data[5]
            self.mu1 = data[6]
            self.mu2 = data[7]
            self.tmax = data[8]
            self.mar = data[9]
            self.trmin = data[10]
        ###

    def cross_reference(self, model):
        pass

    def rawFields(self):
        fields = ['PGAP', self.pid, self.u0, self.f0, self.ka, self.kb, self.kt, self.mu1, self.mu2,
                  self.tmax, self.mar, self.trmin]
        return fields

    def reprFields(self):
        u0 = set_blank_if_default(self.u0, 0.)
        f0 = set_blank_if_default(self.f0, 0.)
        ka = set_blank_if_default(self.ka, 1.e8)
        kb = set_blank_if_default(self.kb, 1e-14 * self.ka)
        kt = set_blank_if_default(self.kt, self.mu1 * self.ka)
        mu1 = set_blank_if_default(self.mu1, 0.)
        mu2 = set_blank_if_default(self.mu2, self.mu1)
        tmax = set_blank_if_default(self.tmax, 0.)
        mar = set_blank_if_default(self.mar, 100.)
        trmin = set_blank_if_default(self.trmin, 0.001)

        fields = ['PGAP', self.pid, u0, f0, ka, kb, kt, mu1, mu2,
                  tmax, mar, trmin]
        return fields


class SolidProperty(Property):
    type = 'SolidProperty'

    def __init__(self, card, data):
        Property.__init__(self, card, data)

    def Rho(self):
        self.mid.rho


class PLSOLID(SolidProperty):
    """
    Defines a fully nonlinear (i.e., large strain and large rotation) hyperelastic solid
    element.
    PLSOLID PID MID STR
    PLSOLID 20 21
    """
    type = 'PLSOLID'

    def __init__(self, card=None, data=None):
        SolidProperty.__init__(self, card, data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            ## Material ID
            self.mid = card.field(2)
            self.ge = card.field(3)
            self.str = card.field(4, 'GRID')
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.ge = data[2]
            self.str = data[3]
        ###
        assert self.str in ['GRID', 'GAUS'], 'STR=|%s| doesnt have a valid stress/strain output value set\n' % (self.str)

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def rawFields(self):
        stressStrain = set_blank_if_default(self.str, 'GRID')
        fields = ['PLSOLID', self.pid, self.Mid(), stressStrain]
        return fields


class PSOLID(SolidProperty):
    """
    PSOLID PID MID CORDM IN STRESS ISOP FCTN
    PSOLID   1       1       0
    PSOLID 2 100 6 TWO GRID REDUCED
    """
    type = 'PSOLID'

    def __init__(self, card=None, data=None):
        SolidProperty.__init__(self, card, data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            ## Material ID
            self.mid = card.field(2)
            self.cordm = card.field(3, 0)
            self.integ = card.field(4)
            #validIntegration = ['THREE','TWO','FULL','BUBBLE',2,3,None,'REDUCED']
            self.stress = card.field(5)
            self.isop = card.field(6)
            self.fctn = card.field(7, 'SMECH')
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.cordm = data[2]
            self.integ = data[3]
            self.stress = data[4]
            self.isop = data[5]
            self.fctn = data[6]

            if self.fctn == 'SMEC':
                self.fctn = 'SMECH'
        ###

    def writeCalculix(self, elementSet=999):
        msg = '*SOLID SECTION,MATERIAL=M%s,ELSET=E_Mat%s\n' % (
            self.mid, elementSet)
        return msg

    def rawFields(self):
        fields = ['PSOLID', self.pid, self.Mid(), self.cordm,
                  self.integ, self.stress, self.isop, self.fctn]
        return fields

    def reprFields(self):
        cordm = set_blank_if_default(self.cordm, 0)
        fctn = set_blank_if_default(self.fctn, 'SMECH')
        fields = ['PSOLID', self.pid, self.Mid(), cordm,
                  self.integ, self.stress, self.isop, fctn]
        return fields


class CrackProperty(Property):
    def __init__(self, card, data):
        pass

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        return self.mid.mid


class PRAC2D(CrackProperty):
    """
    CRAC2D Element Property
    Defines the properties and stress evaluation techniques to be used with the CRAC2D
    structural element.
    """
    type = 'PRAC2D'

    def __init__(self, card=None, data=None):
        CrackProperty.__init__(self, card, data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            ## Material ID
            self.mid = card.field(2)
            self.thick = card.field(3)
            ## Plane strain or plane stress option. Use 0 for plane strain; 1 for plane
            ## stress. (Integer = 0 or 1)
            self.iPlane = card.field(4)
            assert self.iPlane in [0, 1], 'Invalid value for iPlane on PRAC2D, can only be 0,1 iPlane=|%s|' % (self.iPlane)
            ## Non-structural mass per unit area.(Real >= 0.0; Default = 0)
            self.nsm = card.field(5, 0.)
            ## Exponent used in the displacement field. See Remark 4. (Real; Default = 0.5)
            self.gamma = card.field(6, 0.5)
            ## Angle (in degrees) relative to the element x-axis along which stress
            ## intensity factors are to be calculated. See Remark 4. (Real; Default = 180.0)
            self.phi = card.field(7, 180.)

        else:
            raise NotImplementedError('not supported')
        ###

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)  # MAT1, MAT2, MAT8

    def rawFields(self):
        fields = ['PRAC2D', self.pid, self.Mid(), self.thick,
                  self.iPlane, self.nsm, self.gamma, self.phi]
        return fields

    def reprFields(self):
        nsm = set_blank_if_default(self.nsm, 0.)
        gamma = set_blank_if_default(self.gamma, 0.5)
        phi = set_blank_if_default(self.phi, 180.)
        fields = ['PRAC2D', self.pid, self.Mid(), self.thick,
                  self.iPlane, nsm, gamma, phi]
        return fields


class PRAC3D(CrackProperty):
    """
    CRAC3D Element Property
    Defines the properties of the CRAC3D structural element.
    """
    type = 'PRAC3D'

    def __init__(self, card=None, data=None):
        CrackProperty.__init__(self, card, data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            ## Material ID
            self.mid = card.field(2)
            ## Exponent used in the displacement field. See Remark 4.
            ## (Real; Default = 0.5)
            self.gamma = card.field(3, 0.5)
            ## Angle (in degrees) relative to the element x-axis along which
            ## stress intensity factors are to be calculated. See Remark 4.
            ## (Real; Default = 180.0)
            self.phi = card.field(4, 180.)

        else:
            raise NotImplementedError('not supported')
        ###

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)  # MAT1, MAT9

    def rawFields(self):
        fields = ['PRAC3D', self.pid, self.Mid(), self.gamma, self.phi]
        return fields

    def reprFields(self):
        gamma = set_blank_if_default(self.gamma, 0.5)
        phi = set_blank_if_default(self.phi, 180.)
        fields = ['PRAC3D', self.pid, self.Mid(), gamma, phi]
        return fields


class PCONEAX(Property):  # not done
    type = 'PCONEAX'

    def __init__(self, card=None, data=None):
        Property.__init__(self, card, data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            ## Material ID
            self.mid = card.field(2)
            self.group = card.field(3, 'MSCBMLO')
            self.Type = card.field(4)
            self.dim = []  # confusing entry...
        else:
            raise NotImplementedError('not supported')
        ###

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def rawFields(self):
        fields = ['PCONEAX', self.pid, self.Mid(), self.group, self.Type]
        raise NotImplementedError('not supported')
        return fields

