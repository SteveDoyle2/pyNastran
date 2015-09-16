#pylint: disable=C0103,C0111,R0902
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
from six import integer_types
from six.moves import range

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Property, Material
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string_or_blank, integer_string_or_blank, blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class PFAST(Property):
    type = 'PFAST'
    _field_map = {
        1: 'pid', 2:'d', 3:'mcid', 4:'mflag',
        5:'kt1', 6:'kt2', 7:'kt3',
        8:'kr1', 9:'kr2', 10:'kr3',
        11:'mass', 12:'ge'
    }

    def __init__(self, card=None, data=None, comment=''):
        Property.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')

            #: diameter of the fastener
            self.d = double(card, 2, 'd')
            assert self.d > 0

            #: Specifies the element stiffness coordinate system
            self.mcid = integer_or_blank(card, 3, 'mcid', -1)
            assert self.mcid >= -1

            #: 0-absolute 1-relative
            self.mflag = integer_or_blank(card, 4, 'mflag', 0)
            assert self.mflag in [0, 1]

            #: stiffness values in directions 1-3
            self.kt1 = double(card, 5, 'kt1')
            self.kt2 = double(card, 6, 'kt2')
            self.kt3 = double(card, 7, 'kt3')

            #: Rotational stiffness values in directions 1-3
            self.kr1 = double_or_blank(card, 8, 'kr1', 0.0)
            self.kr2 = double_or_blank(card, 9, 'kr2', 0.0)
            self.kr3 = double_or_blank(card, 10, 'kr3', 0.0)

            #: Lumped mass of fastener
            self.mass = double_or_blank(card, 11, 'mass', 0.0)

            #: Structural damping
            self.ge = double_or_blank(card, 12, 'ge', 0.0)
            assert len(card) <= 13, 'len(PFAST card) = %i' % len(card)
        else:
            raise NotImplementedError(data)

    def cross_reference(self, model):
        msg = ' which is required by PFAST pid=%s' % self.pid
        if self.mcid != -1:
            self.mcid = model.Coord(self.Mcid(), msg)

    def Mcid(self):
        if isinstance(self.mcid, integer_types):
            return self.mcid
        return self.mcid.cid

    def Mass(self):
        return self.mass

    def raw_fields(self):
        fields = ['PFAST', self.pid, self.d, self.Mcid(), self.mflag, self.kt1,
                  self.kt2, self.kt3, self.kr1, self.kr2, self.kr3, self.mass,
                  self.ge]
        return fields

    def repr_fields(self):
        mcid = set_blank_if_default(self.Mcid(), -1)
        mflag = set_blank_if_default(self.mflag, 0)
        kr1 = set_blank_if_default(self.kr1, 0.0)
        kr2 = set_blank_if_default(self.kr2, 0.0)
        kr3 = set_blank_if_default(self.kr3, 0.0)

        mass = set_blank_if_default(self.mass, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        fields = ['PFAST', self.pid, self.d, mcid, mflag, self.kt1, self.kt2,
                  self.kt3, kr1, kr2, kr3, mass, ge]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PGAP(Property):
    type = 'PGAP'
    _field_map = {
        1: 'pid', 2:'u0', 3:'f0', 4:'ka', 5:'kb', 6:'kt', 7:'mu1',
        8:'mu2', 9:'tmax', 10:'mar', 11:'trmin',
    }

    def __init__(self, card=None, data=None, comment=''):
        """
        Defines the properties of the gap element (CGAP entry).
        """
        Property.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: initial gap opening
            self.u0 = double_or_blank(card, 2, 'u0', 0.)
            #: preload
            self.f0 = double_or_blank(card, 3, 'f0', 0.)
            #: axial stiffness of closed gap
            self.ka = double_or_blank(card, 4, 'ka', 1.e8)
            #: axial stiffness of open gap
            self.kb = double_or_blank(card, 5, 'kb', 1e-14 * self.ka)
            #: static friction coeff
            self.mu1 = double_or_blank(card, 7, 'mu1', 0.)
            #: transverse stiffness of closed gap
            self.kt = double_or_blank(card, 6, 'kt', self.mu1 * self.ka)
            #: kinetic friction coeff
            self.mu2 = double_or_blank(card, 8, 'mu2', self.mu1)
            self.tmax = double_or_blank(card, 9, 'tmax', 0.)
            self.mar = double_or_blank(card, 10, 'mar', 100.)
            self.trmin = double_or_blank(card, 11, 'trmin', 0.001)
            assert len(card) <= 12, 'len(PGAP card) = %i' % len(card)
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

    def _verify(self, xref=True):
        pid = self.Pid()
        assert isinstance(pid, int), 'pid=%r\n%s' % (pid, str(self))

    def cross_reference(self, model):
        pass

    def raw_fields(self):
        fields = ['PGAP', self.pid, self.u0, self.f0, self.ka, self.kb,
                  self.kt, self.mu1, self.mu2, self.tmax, self.mar, self.trmin]
        return fields

    def repr_fields(self):
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

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class SolidProperty(Property):
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
    _field_map = {
        1: 'pid', 2:'mid', 3:'str',
    }

    def __init__(self, card=None, data=None, comment=''):
        SolidProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: Material ID
            self.mid = integer(card, 2, 'mid')
            #: Location of stress and strain output
            self.str = string_or_blank(card, 3, 'str', 'GRID')
            assert len(card) <= 4, 'len(PLSOLID card) = %i' % len(card)
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.ge = data[2]
            self.str = data[3]

        if self.str == 'GAUS':
            self.str = 'GAUSS'
        if self.str not in ['GRID', 'GAUSS']:
            raise RuntimeError('STR="%s" doesnt have a valid stress/strain '
                               'output value set; valid=["GRID", "GAUSS"]\n'
                               % self.str)

    def cross_reference(self, model):
        msg = ' which is required by PLSOLID pid=%s' % self.pid
        self.mid = model.HyperelasticMaterial(self.mid, msg)

    def raw_fields(self):
        stress_strain = set_blank_if_default(self.str, 'GRID')
        fields = ['PLSOLID', self.pid, self.Mid(), stress_strain]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PSOLID(SolidProperty):
    """
    PSOLID PID MID CORDM IN STRESS ISOP FCTN
    PSOLID   1       1       0
    PSOLID 2 100 6 TWO GRID REDUCED
    """
    type = 'PSOLID'
    _field_map = {
        1: 'pid', 2:'mid', 3:'cordm', 4:'integ', 5:'stress',
        6:'isop', 7:'fctn',
    }

    def __init__(self, card=None, data=None, comment=''):
        SolidProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: Material ID
            self.mid = integer(card, 2, 'mid')
            self.cordm = integer_or_blank(card, 3, 'cordm', 0)
            self.integ = integer_string_or_blank(card, 4, 'integ')
            #validIntegration = ['THREE', 'TWO', 'FULL', 'BUBBLE',
            #                    2, 3, None, 'REDUCED']
            self.stress = integer_string_or_blank(card, 5, 'stress')
            self.isop = integer_string_or_blank(card, 6, 'isop')
            self.fctn = string_or_blank(card, 7, 'fctn', 'SMECH')
            assert len(card) <= 8, 'len(PSOLID card) = %i' % len(card)
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

    def E(self):
        return self.mid.E()

    def G(self):
        return self.mid.G()

    def Nu(self):
        return self.mid.Nu()

    def materials(self):
        return [self.mid]

    def _verify(self, xref=False):
        pid = self.Pid()
        mid = self.Mid()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mid, int), 'mid=%r' % mid

        if xref:
            if self.mid.type not in ['MAT1', 'MAT4', 'MAT5', 'MAT9', 'MAT10', 'MAT11']:
                msg = 'mid=%i self.mid.type=%s' % (mid, self.mid.type)
                raise TypeError(msg)

    def _write_calculix(self, elementSet=999):
        msg = '*SOLID SECTION,MATERIAL=M%s,ELSET=E_Mat%s\n' % (
            self.mid, elementSet)
        return msg

    def raw_fields(self):
        fields = ['PSOLID', self.pid, self.Mid(), self.cordm, self.integ,
                  self.stress, self.isop, self.fctn]
        return fields

    def repr_fields(self):
        cordm = set_blank_if_default(self.cordm, 0)
        fctn = set_blank_if_default(self.fctn, 'SMECH')
        fields = ['PSOLID', self.pid, self.Mid(), cordm, self.integ,
                  self.stress, self.isop, fctn]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        # this card has integers & strings, so it uses...
        return self.comment + print_card_8(card)

class CrackProperty(Property):
    def __init__(self, card, data):
        Property.__init__(self, card, data)

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        return self.mid.mid

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PRAC2D(CrackProperty):
    """
    CRAC2D Element Property
    Defines the properties and stress evaluation techniques to be used with
    the CRAC2D structural element.
    """
    type = 'PRAC2D'
    _field_map = {
        1: 'pid', 2:'mid', 3:'thick', 4:'iPlane', 5:'nsm', 6:'gamma', 7:'phi',
    }

    def __init__(self, card=None, data=None, comment=''):
        CrackProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: Material ID
            self.mid = integer(card, 2, 'mid')
            self.thick = double(card, 3, 'thick')
            #: Plane strain or plane stress option.
            #: Use 0 for plane strain; 1 for plane stress. (Integer = 0 or 1)
            self.iPlane = integer(card, 4, 'iPlane')
            if self.iPlane not in [0, 1]:
                raise RuntimeError('Invalid value for iPlane on PRAC2D, can '
                                   'only be 0,1 iPlane=|%s|' % self.iPlane)
            #: Non-structural mass per unit area.(Real >= 0.0; Default = 0)
            self.nsm = double_or_blank(card, 5, 'nsm', 0.)
            #: Exponent used in the displacement field. See Remark 4.
            #: (Real; Default = 0.5)
            self.gamma = double_or_blank(card, 6, 'gamma', 0.5)
            #: Angle (in degrees) relative to the element x-axis along which
            #: stress intensity factors are to be calculated. See Remark 4.
            #: (Real; Default = 180.0)
            self.phi = double_or_blank(card, 7, 'phi', 180.)
            assert len(card) <= 8, 'len(PRAC2D card) = %i' % len(card)

        else:
            raise NotImplementedError(data)

    def _verify(self, xref=True):
        pid = self.Pid()
        assert isinstance(pid, int)

    def cross_reference(self, model):
        msg = ' which is required by PRAC2D pid=%s' % self.pid
        self.mid = model.Material(self.mid, msg)  # MAT1, MAT2, MAT8

    def raw_fields(self):
        fields = ['PRAC2D', self.pid, self.Mid(), self.thick,
                  self.iPlane, self.nsm, self.gamma, self.phi]
        return fields

    def repr_fields(self):
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
    _field_map = {
        1: 'pid', 2:'mid', 3:'gamma', 4:'phi',
    }

    def __init__(self, card=None, data=None, comment=''):
        CrackProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: Material ID
            self.mid = integer(card, 2, 'mid')
            #: Exponent used in the displacement field. See Remark 4.
            #: (Real; Default = 0.5)
            self.gamma = double_or_blank(card, 3, 'gamma', 0.5)
            #: Angle (in degrees) relative to the element x-axis along which
            #: stress intensity factors are to be calculated. See Remark 4.
            #: (Real; Default = 180.0)
            self.phi = double_or_blank(card, 4, 'gamma', 180.)
            assert len(card) <= 5, 'len(PRAC3D card) = %i' % len(card)
        else:
            raise NotImplementedError(data)

    def _verify(self, xref=True):
        pid = self.Pid()
        assert isinstance(pid, int)

    def cross_reference(self, model):
        msg = ' which is required by PRAC3D pid=%s' % self.pid
        self.mid = model.Material(self.mid, msg)  # MAT1, MAT9

    def raw_fields(self):
        fields = ['PRAC3D', self.pid, self.Mid(), self.gamma, self.phi]
        return fields

    def repr_fields(self):
        gamma = set_blank_if_default(self.gamma, 0.5)
        phi = set_blank_if_default(self.phi, 180.)
        fields = ['PRAC3D', self.pid, self.Mid(), gamma, phi]
        return fields


class PCONEAX(Property):
    type = 'PCONEAX'
    _field_map = {
        1: 'pid', 2:'mid1', 3:'t1', 4:'mid2', 5:'i', 6:'mid3', 7:'t2',
        8: 'nsm', 9:'z1', 10:'z2',
    }
    def _update_field_helper(self, n, value):
        if n <= 0:
            msg = 'Field %r=%r is an invalid %s entry.' % (n, value, self.type)
            raise KeyError(msg)
        self.phi[n - 10] = value

    def __init__(self, card=None, data=None, comment=''):
        Property.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: Material ID
            self.mid1 = integer_or_blank(card, 2, 'mid1', 0)
            self.t1 = double_or_blank(card, 3, 't1')

            self.mid2 = integer_or_blank(card, 4, 'mid2', 0)
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

            j = 1
            self.phi = []
            for i in range(11, len(card)):
                phii = double(card, i, 'phi' % j)
                self.phi.append(phii)
                j += 1
        else:
            raise NotImplementedError(data)

    def cross_reference(self, model):
        msg = ' which is required by %s=%s' %(self.type, self.pid)
        if self.mid1 > 0:
            self.mid1 = model.Material(self.mid1, msg=msg)
        if self.mid2 > 0:
            self.mid2 = model.Material(self.mid2, msg=msg)
        if self.mid3 > 0:
            self.mid3 = model.Material(self.mid3, msg=msg)

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

    def raw_fields(self):
        fields = ['PCONEAX', self.pid, self.Mid1(), self.t1,
                  self.Mid2(), self.i, self.Mid3(), self.t2,
                  self.nsm, self.z1, self.z2] + self.phi
        return fields

    def repr_fields(self):
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

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
