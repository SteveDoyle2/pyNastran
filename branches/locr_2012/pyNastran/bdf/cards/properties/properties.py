# pylint: disable=C0103,R0902,R0904,R0914,C0111
"""
All ungrouped properties are defined in this file.  This includes:
 * PFAST
 * PGAP
 * PLSOLID (SolidProperty)
 * PSOLID (SolidProperty)
 * PRAC2D (CrackProperty)
 * PRAC3D (CrackProperty)
 * PCONEAX (not done)
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Property
from pyNastran.bdf.format import (integer, integer_or_blank,
                                  double, double_or_blank,
                                  string_or_blank, integer_string_or_blank)

class PFAST(Property):
    type = 'PFAST'

    def __init__(self, card=None, data=None, comment=''):
        Property.__init__(self, card, data)
        if comment:
            self._comment = comment
        ## Property ID
        self.pid = integer(card, 1, 'pid')
        ## diameter of the fastener
        self.d = double(card, 2, 'd')
        assert self.d > 0
        ## Specifies the element stiffness coordinate system
        self.mcid = integer_or_blank(card, 3, 'mcid', -1)
        assert self.mcid >= -1
        self.mflag = integer_or_blank(card, 4, 'mflag', 0)  # 0-absolute 1-relative
        assert self.mflag in [0, 1]
        ## stiffness values in directions 1-3
        self.kt1 = double(card, 5, 'kt1')
        self.kt2 = double(card, 6, 'kt2')
        self.kt3 = double(card, 7, 'kt3')
        ## Rotational stiffness values in directions 1-3
        self.kr1 = double_or_blank(card, 8, 'kr1', 0.0)
        self.kr2 = double_or_blank(card, 9, 'kr2', 0.0)
        self.kr3 = double_or_blank(card, 10, 'kr3', 0.0)
        ## Lumped mass of fastener
        self.mass = double_or_blank(card, 11, 'mass', 0.0)
        ## Structural damping
        self.ge = double_or_blank(card, 12, 'ge', 0.0)

    def cross_reference(self, model):
        self.mcid = model.Coord(self.mcid)

    def Mcid(self):
        if isinstance(self.mcid, int):
            return self.mcid
        return self.mcid.cid

    def Mass(self):
        return self.mass

    def rawFields(self):
        fields = ['PFAST', self.pid, self.d, self.Mcid(), self.mflag, self.kt1,
                  self.kt2, self.kt3, self.kr1, self.kr2, self.kr3, self.mass,
                  self.ge]
        return fields

    def reprFields(self):
        mcid = set_blank_if_default(self.mcid, -1)
        mflag = set_blank_if_default(self.mflag, 0)
        kr1 = set_blank_if_default(self.kr1, 0.0)
        kr2 = set_blank_if_default(self.kr2, 0.0)
        kr3 = set_blank_if_default(self.kr3, 0.0)

        mass = set_blank_if_default(self.mass, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        fields = ['PFAST', self.pid, self.d, mcid, mflag, self.kt1, self.kt2,
                  self.kt3, kr1, kr2, kr3, mass, ge]
        return fields


class PGAP(Property):
    type = 'PGAP'

    def __init__(self, card=None, data=None, comment=''):
        """
        Defines the properties of the gap element (CGAP entry).
        """
        Property.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            ## Property ID
            self.pid = integer(card, 1, 'pid')
            ## initial gap opening
            self.u0 = double_or_blank(card, 2, 'u0', 0.)
            ## preload
            self.f0 = double_or_blank(card, 3, 'f0', 0.)
            ## axial stiffness of closed gap
            self.ka = double_or_blank(card, 4, 'ka', 1.e8)
            ## axial stiffness of open gap
            self.kb = double_or_blank(card, 5, 'kb', 1e-14 * self.ka)
            ## static friction coeff
            self.mu1 = double_or_blank(card, 7, 'mu1', 0.)
            ## transverse stiffness of closed gap
            self.kt = double_or_blank(card, 6, 'kt', self.mu1 * self.ka)
            ## kinetic friction coeff
            self.mu2 = double_or_blank(card, 8, 'mu2', self.mu1)
            self.tmax = double_or_blank(card, 9, 'tmax', 0.)
            self.mar = double_or_blank(card, 10, 'mar', 100.)
            self.trmin = double_or_blank(card, 11, 'trmin', 0.001)
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

    def cross_reference(self, model):
        pass

    def rawFields(self):
        fields = ['PGAP', self.pid, self.u0, self.f0, self.ka, self.kb,
                  self.kt, self.mu1, self.mu2, self.tmax, self.mar, self.trmin]
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
        return self.mid.rho


class PLSOLID(SolidProperty):
    """
    Defines a fully nonlinear (i.e., large strain and large rotation)
    hyperelastic solid element.
    PLSOLID PID MID STR
    PLSOLID 20 21
    """
    type = 'PLSOLID'

    def __init__(self, card=None, data=None, comment=''):
        SolidProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            ## Property ID
            self.pid = integer(card, 1, 'pid')
            ## Material ID
            self.mid = integer(card, 2, 'mid')
            self.ge = double_or_blank(card, 3, 'ge', 0.0)
            self.str = string_or_blank(card, 4, 'str', 'GRID')
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.ge = data[2]
            self.str = data[3]
        if self.str not in ['GRID', 'GAUS']:
            raise RuntimeError('STR=|%s| doesnt have a valid stress/strain '
                               'output value set\n' % self.str)

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

    def __init__(self, card=None, data=None, comment=''):
        SolidProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            ## Property ID
            self.pid = integer(card, 1, 'pid')
            ## Material ID
            self.mid = integer(card, 2, 'mid')
            self.cordm = integer_or_blank(card, 3, 'cordm', 0)
            self.integ = integer_string_or_blank(card, 4, 'integ')
            #validIntegration = ['THREE', 'TWO', 'FULL', 'BUBBLE',
            #                    2, 3, None, 'REDUCED']
            self.stress = integer_string_or_blank(card, 5, 'stress')
            self.isop = integer_string_or_blank(card, 6, 'isop')
            self.fctn = string_or_blank(card, 7, 'fctn', 'SMECH')
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

    def writeCalculix(self, elementSet=999):
        msg = '*SOLID SECTION,MATERIAL=M%s,ELSET=E_Mat%s\n' % (
            self.mid, elementSet)
        return msg

    def rawFields(self):
        fields = ['PSOLID', self.pid, self.Mid(), self.cordm, self.integ,
                  self.stress, self.isop, self.fctn]
        return fields

    def reprFields(self):
        cordm = set_blank_if_default(self.cordm, 0)
        fctn = set_blank_if_default(self.fctn, 'SMECH')
        fields = ['PSOLID', self.pid, self.Mid(), cordm, self.integ,
                  self.stress, self.isop, fctn]
        return fields


class CrackProperty(Property):
    def __init__(self, card, data):
        Property.__init__(self, card, data)

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        return self.mid.mid


class PRAC2D(CrackProperty):
    """
    CRAC2D Element Property
    Defines the properties and stress evaluation techniques to be used with
    the CRAC2D structural element.
    """
    type = 'PRAC2D'

    def __init__(self, card=None, data=None, comment=''):
        CrackProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            ## Property ID
            self.pid = integer(card, 1, 'pid')
            ## Material ID
            self.mid = integer(card, 2, 'mid')
            self.thick = double(card, 3, 'thick')
            ## Plane strain or plane stress option.
            ## Use 0 for plane strain; 1 for plane stress. (Integer = 0 or 1)
            self.iPlane = integer(card, 4, 'iPlane')
            if self.iPlane not in [0, 1]:
                raise RuntimeError('Invalid value for iPlane on PRAC2D, can '
                                   'only be 0,1 iPlane=|%s|' % self.iPlane)
            ## Non-structural mass per unit area.(Real >= 0.0; Default = 0)
            self.nsm = double_or_blank(card, 5, 'nsm', 0.)
            ## Exponent used in the displacement field. See Remark 4.
            ## (Real; Default = 0.5)
            self.gamma = double_or_blank(card, 6, 'gamma', 0.5)
            ## Angle (in degrees) relative to the element x-axis along which
            ## stress intensity factors are to be calculated. See Remark 4.
            ## (Real; Default = 180.0)
            self.phi = double_or_blank(card, 7, 'phi', 180.)

        else:
            raise NotImplementedError('not supported')

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

    def __init__(self, card=None, data=None, comment=''):
        CrackProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            ## Property ID
            self.pid = integer(card, 1, 'pid')
            ## Material ID
            self.mid = integer(card, 2, 'mid')
            ## Exponent used in the displacement field. See Remark 4.
            ## (Real; Default = 0.5)
            self.gamma = double_or_blank(card, 3, 'gamma', 0.5)
            ## Angle (in degrees) relative to the element x-axis along which
            ## stress intensity factors are to be calculated. See Remark 4.
            ## (Real; Default = 180.0)
            self.phi = double_or_blank(card, 4, 'gamma', 180.)

        else:
            raise NotImplementedError('not supported')

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


class PCONEAX(Property):
    type = 'PCONEAX'

    def __init__(self, card=None, data=None, comment=''):
        Property.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            ## Property ID
            self.pid = integer(card, 1, 'pid')
            ## Material ID
            self.mid1 = integer_or_blank(card, 2, 'mid1', 0)
            self.t1 = double_or_blank(card, 3, 't1')

            self.mid2 = integer(card, 4, 'mid2', 0)
            if self.mid2 > 0:
                self.i = double(card, 5, 'i')
                assert self.i > 0.0
            else:
                self.i = blank(card, 5, 'i')

            self.mid3 = integer(card, 6, 0)
            if self.mid3 > 0:
                self.t2 = double(card, 7, 't3')
                assert self.t2 > 0.0
            else:
                self.t2 = blank(card, 7, 't3')

            self.nsm = double(card, 8, 'nsm')
            self.z1 = double(card, 9, 'z1')
            self.z2 = double(card, 10, 'z2')
            self.phi = fields(double, card, 'phi', i=11, j=len(card))
        else:
            raise NotImplementedError('not supported')

    def cross_reference(self, model):
        if self.mid1 > 0:
            self.mid1 = model.Material(self.mid1, msg=msg)
        if self.mid2 > 0:
            self.mid2 = model.Material(self.mid2, msg=msg)
        if self.mid3 > 0:
            self.mid3 = model.Material(self.mid3, msg=msg)
        if self.mid4 > 0:
            self.mid4 = model.Material(self.mid4, msg=msg)

    def Mid1(self):
        if isinstance(self.mid1, Material):
            return self.mid1.mid
        return self.mid1

    def Mid2(self):
        if isinstance(self.mid2, Material):
            return self.mid2.mid
        return self.mid2

    def Mid3(self):
        if isinstance(self.mid3, Material):
            return self.mid3.mid
        return self.mid3

    def rawFields(self):
        fields = ['PCONEAX', self.pid, self.Mid1(), self.t1,
                  self.Mid2(), self.i, self.Mid3(), self.t2,
                  self.nsm, self.z1, self.z2] + self.phi
        return fields

    def reprFields(self):
        nsm = set_blank_if_default(self.nsm, 0.0)
        mid1 = set_blank_if_default(self.Mid1(), 0)
        mid2 = set_blank_if_default(self.Mid2(), 0)
        mid3 = set_blank_if_default(self.Mid3(), 0)
        i = set_blank_if_default(self.i, 0.0)
        t1 = set_blank_if_default(self.t1, 0.0)
        t2 = set_blank_if_default(self.t2, 0.0)
        fields = ['PCONEAX', self.pid, mid1, t1, mid2, i, mid3, t2,
                  nsm, self.z1, self.z2] + self.phi
        return fields
