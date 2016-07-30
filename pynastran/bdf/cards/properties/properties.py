#pylint: disable=C0103,C0111,R0902
"""
All ungrouped properties are defined in this file.  This includes:
 * PFAST
 * PGAP
 * PRAC2D (CrackProperty)
 * PRAC3D (CrackProperty)
 * PCONEAX (not done)
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import range

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Property, Material
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank,
    blank)
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

    def __init__(self, pid, d, mcid, mflag, kt1, kt2, kt3,
                 kr1, kr2, kr3, mass, ge, comment=''):
        Property.__init__(self)
        if comment:
            self._comment = comment
        #: Property ID
        self.pid = pid
        #: diameter of the fastener
        self.d = d
        #: Specifies the element stiffness coordinate system
        self.mcid = mcid
        #: 0-absolute 1-relative
        self.mflag = mflag

        #: stiffness values in directions 1-3
        self.kt1 = kt1
        self.kt2 = kt2
        self.kt3 = kt3

        #: Rotational stiffness values in directions 1-3
        self.kr1 = kr1
        self.kr2 = kr2
        self.kr3 = kr3
        #: Lumped mass of fastener
        self.mass = mass
        #: Structural damping
        self.ge = ge
        assert self.d > 0
        assert mflag in [0, 1]
        assert self.mcid >= -1

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')
        d = double(card, 2, 'd')
        mcid = integer_or_blank(card, 3, 'mcid', -1)
        mflag = integer_or_blank(card, 4, 'mflag', 0)

        kt1 = double(card, 5, 'kt1')
        kt2 = double(card, 6, 'kt2')
        kt3 = double(card, 7, 'kt3')

        kr1 = double_or_blank(card, 8, 'kr1', 0.0)
        kr2 = double_or_blank(card, 9, 'kr2', 0.0)
        kr3 = double_or_blank(card, 10, 'kr3', 0.0)
        mass = double_or_blank(card, 11, 'mass', 0.0)
        ge = double_or_blank(card, 12, 'ge', 0.0)
        assert len(card) <= 13, 'len(PFAST card) = %i\ncard=%s' % (len(card), card)
        return PFAST(pid, d, mcid, mflag, kt1, kt2, kt3,
                     kr1, kr2, kr3, mass, ge, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        (pid, d, mcid, mflag, kt1, kt2, kt3,
         kr1, kr2, kr3, mass, ge) = data
        return PFAST(pid, d, mcid, mflag, kt1, kt2, kt3,
                     kr1, kr2, kr3, mass, ge, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PFAST pid=%s' % self.pid
        if self.mcid != -1:
            self.mcid = model.Coord(self.Mcid(), msg)
            self.mcid_ref = self.mcid_ref

    def uncross_reference(self):
        self.mcid = self.Mcid()
        if self.mcid != -1:
            del self.mcid_ref

    def Mcid(self):
        if isinstance(self.mcid, integer_types):
            return self.mcid
        return self.mcid_ref.cid

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

    def __init__(self, pid, u0, f0, ka, kb, mu1, kt, mu2,
                 tmax, mar, trmin, comment=''):
        """
        Defines the properties of the gap element (CGAP entry).
        """
        Property.__init__(self)
        if comment:
            self._comment = comment

        #: Property ID
        self.pid = pid
        #: initial gap opening
        self.u0 = u0
        #: preload
        self.f0 = f0
        #: axial stiffness of closed gap
        self.ka = ka
        #: axial stiffness of open gap
        self.kb = kb
        #: static friction coeff
        self.mu1 = mu1
        #: transverse stiffness of closed gap
        self.kt = kt
        #: kinetic friction coeff
        self.mu2 = mu2
        self.tmax = tmax
        self.mar = mar
        self.trmin = trmin

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')
        u0 = double_or_blank(card, 2, 'u0', 0.)
        f0 = double_or_blank(card, 3, 'f0', 0.)
        ka = double_or_blank(card, 4, 'ka', 1.e8)
        kb = double_or_blank(card, 5, 'kb', 1e-14 * ka)
        mu1 = double_or_blank(card, 7, 'mu1', 0.)
        kt = double_or_blank(card, 6, 'kt', mu1 * ka)
        mu2 = double_or_blank(card, 8, 'mu2', mu1)
        tmax = double_or_blank(card, 9, 'tmax', 0.)
        mar = double_or_blank(card, 10, 'mar', 100.)
        trmin = double_or_blank(card, 11, 'trmin', 0.001)
        assert len(card) <= 12, 'len(PGAP card) = %i\ncard=%s' % (len(card), card)
        return PGAP(pid, u0, f0, ka, kb, mu1, kt, mu2, tmax, mar, trmin,
                    comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        pid = data[0]
        u0 = data[1]
        f0 = data[2]
        ka = data[3]
        kb = data[4]
        kt = data[5]
        mu1 = data[6]
        mu2 = data[7]
        tmax = data[8]
        mar = data[9]
        trmin = data[10]
        return PGAP(pid, u0, f0, ka, kb, mu1, kt, mu2, tmax, mar, trmin,
                    comment=comment)

    def _verify(self, xref=True):
        pid = self.Pid()
        assert isinstance(pid, int), 'pid=%r\n%s' % (pid, str(self))

    def cross_reference(self, model):
        pass

    def uncross_reference(self):
        pass

    def raw_fields(self):
        fields = ['PGAP', self.pid, self.u0, self.f0, self.ka, self.kb,
                  self.kt, self.mu1, self.mu2, self.tmax, self.mar, self.trmin]
        return fields

    def repr_fields(self):
        u0 = set_blank_if_default(self.u0, 0.)
        f0 = set_blank_if_default(self.f0, 0.)
        # ka doesn't have a default in MSC 2005r2
        #ka = set_blank_if_default(self.ka, 1.e8)
        kb = set_blank_if_default(self.kb, 1e-14 * self.ka)
        kt = set_blank_if_default(self.kt, self.mu1 * self.ka)
        mu1 = set_blank_if_default(self.mu1, 0.)
        mu2 = set_blank_if_default(self.mu2, self.mu1)
        tmax = set_blank_if_default(self.tmax, 0.)
        mar = set_blank_if_default(self.mar, 100.)
        trmin = set_blank_if_default(self.trmin, 0.001)

        fields = ['PGAP', self.pid, u0, f0, self.ka, kb, kt, mu1, mu2,
                  tmax, mar, trmin]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class CrackProperty(Property):
    def __init__(self):
        Property.__init__(self)

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        return self.mid_ref.mid

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

    def __init__(self, pid, mid, thick, iplane, nsm=0., gamma=0.5, phi=180.,
                 comment=''):
        CrackProperty.__init__(self)
        if comment:
            self._comment = comment
        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid
        self.thick = thick
        #: Plane strain or plane stress option.
        #: Use 0 for plane strain; 1 for plane stress. (Integer = 0 or 1)
        self.iplane = iplane
        #: Non-structural mass per unit area.(Real >= 0.0; Default = 0)
        self.nsm = nsm
        #: Exponent used in the displacement field. See Remark 4.
        #: (Real; Default = 0.5)
        self.gamma = gamma
        #: Angle (in degrees) relative to the element x-axis along which
        #: stress intensity factors are to be calculated. See Remark 4.
        #: (Real; Default = 180.0)
        self.phi = phi

        if iplane not in [0, 1]:
            raise RuntimeError('Invalid value for iPlane on PRAC2D, can '
                               'only be 0,1 iPlane=%r' % iplane)

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        thick = double(card, 3, 'thick')
        iplane = integer(card, 4, 'iplane')
        nsm = double_or_blank(card, 5, 'nsm', 0.)
        gamma = double_or_blank(card, 6, 'gamma', 0.5)
        phi = double_or_blank(card, 7, 'phi', 180.)
        assert len(card) <= 8, 'len(PRAC2D card) = %i\ncard=%s' % (len(card), card)
        return PRAC2D(pid, mid, thick, iplane, nsm, gamma, phi,
                      comment=comment)

    def _verify(self, xref=True):
        pid = self.Pid()
        assert isinstance(pid, int)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PRAC2D pid=%s' % self.pid
        self.mid = model.Material(self.mid, msg)  # MAT1, MAT2, MAT8
        self.mid_ref = self.mid

    def uncross_reference(self):
        self.mid = self.Mid()
        del self.mid_ref

    def raw_fields(self):
        fields = ['PRAC2D', self.pid, self.Mid(), self.thick,
                  self.iplane, self.nsm, self.gamma, self.phi]
        return fields

    def repr_fields(self):
        nsm = set_blank_if_default(self.nsm, 0.)
        gamma = set_blank_if_default(self.gamma, 0.5)
        phi = set_blank_if_default(self.phi, 180.)
        fields = ['PRAC2D', self.pid, self.Mid(), self.thick,
                  self.iplane, nsm, gamma, phi]
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

    def __init__(self, pid, mid, gamma=0.5, phi=180., comment=''):
        CrackProperty.__init__(self)
        if comment:
            self._comment = comment
        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid
        #: Exponent used in the displacement field. See Remark 4.
        #: (Real; Default = 0.5)
        self.gamma = gamma
        #: Angle (in degrees) relative to the element x-axis along which
        #: stress intensity factors are to be calculated. See Remark 4.
        #: (Real; Default = 180.0)
        self.phi = phi

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        gamma = double_or_blank(card, 3, 'gamma', 0.5)
        phi = double_or_blank(card, 4, 'gamma', 180.)
        assert len(card) <= 5, 'len(PRAC3D card) = %i\ncard=%s' % (len(card), card)
        return PRAC3D(pid, mid, gamma, phi, comment=comment)

    def _verify(self, xref=True):
        pid = self.Pid()
        assert isinstance(pid, int)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PRAC3D pid=%s' % self.pid
        self.mid = model.Material(self.mid, msg)  # MAT1, MAT9
        self.mid_ref = self.mid

    def uncross_reference(self):
        self.mid = self.Mid()
        del self.mid_ref

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

    def __init__(self, pid, mid1, t1, mid2, i, mid3, t2, nsm, z1, z2, phi, comment=''):
        Property.__init__(self)
        if comment:
            self._comment = comment
        self.pid = pid
        self.mid1 = mid1
        self.t1 = t1
        self.mid2 = mid2
        self.i = i
        self.mid3 = mid3
        self.t2 = t2
        self.nsm = nsm
        self.z1 = z1
        self.z2 = z2
        self.phi = phi

    @classmethod
    def add_card(cls, card, comment=''):
        #: Property ID
        pid = integer(card, 1, 'pid')
        #: Material ID
        mid1 = integer_or_blank(card, 2, 'mid1', 0)
        t1 = double_or_blank(card, 3, 't1')

        mid2 = integer_or_blank(card, 4, 'mid2', 0)
        if mid2 > 0:
            i = double(card, 5, 'i')
            assert i > 0.0
        else:
            i = blank(card, 5, 'i')

        mid3 = integer(card, 6, 0)
        if mid3 > 0:
            t2 = double(card, 7, 't3')
            assert t2 > 0.0
        else:
            t2 = blank(card, 7, 't3')

        nsm = double(card, 8, 'nsm')
        z1 = double(card, 9, 'z1')
        z2 = double(card, 10, 'z2')

        j = 1
        phi = []
        for i in range(11, len(card)):
            phii = double(card, i, 'phi' % j)
            phi.append(phii)
            j += 1
        return PCONEAX(pid, mid1, t1, mid2, i, mid3, t2, nsm, z1, z2, phi,
                       comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s=%s' %(self.type, self.pid)
        if self.mid1 > 0:
            self.mid1 = model.Material(self.mid1, msg=msg)
            self.mid1_ref = self.mid1
        if self.mid2 > 0:
            self.mid2 = model.Material(self.mid2, msg=msg)
            self.mid2_ref = self.mid2
        if self.mid3 > 0:
            self.mid3 = model.Material(self.mid3, msg=msg)
            self.mid3_ref = self.mid3

    def uncross_reference(self):
        self.mid1 = self.Mid1()
        self.mid2 = self.Mid2()
        self.mid3 = self.Mid3()
        del self.mid1_ref, self.mid2_ref, self.mid3_ref

    def Mid1(self):
        if isinstance(self.mid1, Material):
            return self.mid1_ref.mid
        return self.mid1

    def Mid2(self):
        if isinstance(self.mid2, Material):
            return self.mid2_ref.mid
        return self.mid2

    def Mid3(self):
        if isinstance(self.mid3, Material):
            return self.mid3_ref.mid
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
