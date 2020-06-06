#pylint: disable=C0103
"""
All ungrouped properties are defined in this file.  This includes:
 * PFAST
 * PGAP
 * PRAC2D (CrackProperty)
 * PRAC3D (CrackProperty)
 * PCONEAX (not done)

"""
from __future__ import annotations
from typing import TYPE_CHECKING
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class PFAST(Property):
    """
    +------+-----+-----+------+-------+---------+--------+--------+-----+
    |   1  |  2  |  3  |   4  |   5   |    6    |    7   |    8   |  9  |
    +======+=====+=====+======+=======+=========+========+========+=====+
    |PFAST | PID |  D  | MCID | MFLAG |   KT1   |   KT2  |   KT3  | KR1 |
    +------+-----+-----+------+-------+---------+--------+--------+-----+
    |      | KR2 | KR3 | MASS |   GE  |         |        |        |     |
    +------+-----+-----+------+-------+---------+--------+--------+-----+
    |PFAST |  7  | 1.1 | 70   |       | 100000. | 46000. | 12300. |     |
    +------+-----+-----+------+-------+---------+--------+--------+-----+
    """
    type = 'PFAST'
    _field_map = {
        1: 'pid', 2:'d', 3:'mcid', 4:'mflag',
        5:'kt1', 6:'kt2', 7:'kt3',
        8:'kr1', 9:'kr2', 10:'kr3',
        11:'mass', 12:'ge'
    }
    pname_fid_map = {
        'KT1' : 'kt1',
        'KT2' : 'kt2',
        'KT3' : 'kt3',
        'KR1' : 'kr1',
        'KR2' : 'kr2',
        'KR3' : 'kr3',
        'MASS' : 'mass',
    }

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        d = 0.1
        kt1 = 1.
        kt2 = 2.
        kt3 = 3.
        return PFAST(pid, d, kt1, kt2, kt3,
                     mcid=-1, mflag=0, kr1=0., kr2=0., kr3=0., mass=0., ge=0., comment='')

    def __init__(self, pid, d, kt1, kt2, kt3, mcid=-1, mflag=0,
                 kr1=0., kr2=0., kr3=0., mass=0., ge=0., comment=''):
        """
        Creates a PAST card

        Parameters
        ----------
        pid : int
            property id
        d : int
            diameter of the fastener
        kt1, kt2, kt3 : float
            stiffness values in directions 1-3
        mcid : int; default=01
            specifies the element stiffness coordinate system
        mflag : int; default=0
            0-absolute; 1-relative
        kr1, kr2, kr3 : float; default=0.0
            rotational stiffness values in directions 1-3
        mass : float; default=0.0
            lumped mass of the fastener
        ge : float; default=0.0
            structural damping
        comment : str; default=''
            a comment for the card
        """
        Property.__init__(self)
        if comment:
            self.comment = comment
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
        self.mcid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PFAST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        return PFAST(pid, d, kt1, kt2, kt3, mcid=mcid, mflag=mflag,
                     kr1=kr1, kr2=kr2, kr3=kr3, mass=mass, ge=ge,
                     comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PFAST card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        (pid, d, mcid, mflag, kt1, kt2, kt3,
         kr1, kr2, kr3, mass, ge) = data
        return PFAST(pid, d, kt1, kt2, kt3, mcid=mcid, mflag=mflag,
                     kr1=kr1, kr2=kr2, kr3=kr3, mass=mass, ge=ge,
                     comment=comment)

    def _verify(self, xref):
        assert isinstance(self.Mcid(), integer_types), self.Mcid()
        assert isinstance(self.Mass(), float), self.mass

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PFAST pid=%s' % self.pid
        if self.mcid != -1:
            self.mcid_ref = model.Coord(self.Mcid(), msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        """Removes cross-reference links"""
        self.mcid = self.Mcid()
        #if self.mcid != -1:
        self.mcid_ref = None

    def Mcid(self):
        if self.mcid_ref is None:
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PGAP(Property):
    """
    +------+------+-------+-------+------+------+------+------+------+
    |   1  |   2  |   3   |   4   |   5  |   6  |   7  |   8  |   9  |
    +======+======+=======+=======+======+======+======+======+======+
    | PGAP |  PID |   U0  |   F0  |  KA  |  KB  |  KT  |  MU1 |  MU2 |
    +------+------+-------+-------+------+------+------+------+------+
    |      | TMAX |  MAR  | TRMIN |      |      |      |      |      |
    +------+------+-------+-------+------+------+------+------+------+
    | PGAP |   2  | 0.025 |  2.5  | 1.E6 |      | 1.E6 | 0.25 | 0.25 |
    +------+------+-------+-------+------+------+------+------+------+
    """

    type = 'PGAP'
    _field_map = {
        1: 'pid', 2:'u0', 3:'f0', 4:'ka', 5:'kb', 6:'kt', 7:'mu1',
        8:'mu2', 9:'tmax', 10:'mar', 11:'trmin',
    }
    pname_fid_map = {
        5 : 'ka',
        6 : 'kb',
        7 : 'kt',
        8 : 'mu1',
        9 : 'mu2',
        10 : 'tmax',
        11 : 'mar',
        12 : 'trmin',
        'KA' : 'ka',
    }

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        return PGAP(pid, u0=0., f0=0., ka=1.e8, kb=None, mu1=0.,
                    kt=None, mu2=None, tmax=0., mar=100., trmin=0.001, comment='')

    def __init__(self, pid, u0=0., f0=0., ka=1.e8, kb=None, mu1=0.,
                 kt=None, mu2=None, tmax=0., mar=100., trmin=0.001,
                 comment=''):
        """
        Defines the properties of the gap element (CGAP entry).

        Parameters
        ----------
        pid : int
            property id for a CGAP
        u0 : float; default=0.
            Initial gap opening
        f0 : float; default=0.
            Preload
        ka : float; default=1.e8
            Axial stiffness for the closed gap
        kb : float; default=None -> 1e-14 * ka
            Axial stiffness for the open gap
        mu1 : float; default=0.
            Coefficient of static friction for the adaptive gap element
            or coefficient of friction in the y transverse direction
            for the nonadaptive gap element
        kt : float; default=None -> mu1*ka
            Transverse stiffness when the gap is closed
        mu2 : float; default=None -> mu1
            Coefficient of kinetic friction for the adaptive gap element
            or coefficient of friction in the z transverse direction
            for the nonadaptive gap element
        tmax : float; default=0.
            Maximum allowable penetration used in the adjustment of
            penalty values. The positive value activates the penalty
            value adjustment
        mar : float; default=100.
            Maximum allowable adjustment ratio for adaptive penalty
            values KA and KT
        trmin : float; default=0.001
            Fraction of TMAX defining the lower bound for the allowable
            penetration
        comment : str; default=''
            a comment for the card
        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        if kb is None:
            kb = 1e-14 * ka
        if kt is None:
            kt = mu1 * ka
        if mu2 is None:
            mu2 = mu1

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
        """
        Adds a PGAP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        """
        Adds a PGAP card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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

    def _verify(self, xref):
        pid = self.Pid()
        assert isinstance(pid, int), 'pid=%r\n%s' % (pid, str(self))

    def cross_reference(self, model: BDF) -> None:
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class CrackProperty(Property):
    def __init__(self):
        Property.__init__(self)

    def Mid(self) -> int:
        if self.mid_ref is not None:
            return self.mid_ref.mid
        return self.mid

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        thick = 0.1
        iplane = 1
        return PRAC2D(pid, mid, thick, iplane, nsm=0., gamma=0.5, phi=180., comment='')

    def __init__(self, pid, mid, thick, iplane, nsm=0., gamma=0.5, phi=180.,
                 comment=''):
        CrackProperty.__init__(self)
        if comment:
            self.comment = comment
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

        self.mid_ref = None

    def validate(self):
        if self.iplane not in [0, 1]:
            raise RuntimeError('Invalid value for iPlane on PRAC2D, can '
                               'only be 0,1 iPlane=%r' % self.iplane)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PRAC2D card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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

    def _verify(self, xref):
        pid = self.Pid()
        assert isinstance(pid, int)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PRAC2D pid=%s' % self.pid
        self.mid_ref = model.Material(self.mid, msg)  # MAT1, MAT2, MAT8

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

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

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        return PRAC3D(pid, mid, gamma=0.5, phi=180., comment='')

    def __init__(self, pid, mid, gamma=0.5, phi=180., comment=''):
        CrackProperty.__init__(self)
        if comment:
            self.comment = comment
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
        self.mid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PRAC3D card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        gamma = double_or_blank(card, 3, 'gamma', 0.5)
        phi = double_or_blank(card, 4, 'gamma', 180.)
        assert len(card) <= 5, 'len(PRAC3D card) = %i\ncard=%s' % (len(card), card)
        return PRAC3D(pid, mid, gamma, phi, comment=comment)

    def _verify(self, xref):
        pid = self.Pid()
        assert isinstance(pid, int)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PRAC3D pid=%s' % self.pid
        self.mid_ref = model.Material(self.mid, msg)  # MAT1, MAT9

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

    def raw_fields(self):
        fields = ['PRAC3D', self.pid, self.Mid(), self.gamma, self.phi]
        return fields

    def repr_fields(self):
        gamma = set_blank_if_default(self.gamma, 0.5)
        phi = set_blank_if_default(self.phi, 180.)
        fields = ['PRAC3D', self.pid, self.Mid(), gamma, phi]
        return fields
