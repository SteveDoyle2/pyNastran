# coding: utf-8
"""
All trim aero cards are defined in this file.  This includes:
 * AEROS
 * AESTAT
 * CSSCHD
 * DIVERG
 * TRIM

All cards are BaseCard objects.

"""
from __future__ import annotations
from itertools import count
from typing import TYPE_CHECKING

from pyNastran.bdf.cards.aero.dynamic_loads import Aero
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string,
    string_or_blank,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class AEROS(Aero):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.

    +-------+-------+-------+------+------+-------+-------+-------+
    |   1   |   2   |   3   |  4   |  5   |   6   |   7   |   8   |
    +=======+=======+=======+======+======+=======+=======+=======+
    | AEROS | ACSID | RCSID | REFC | REFB | REFS  | SYMXZ | SYMXY |
    +-------+-------+-------+------+------+-------+-------+-------+
    | AEROS |   10  |   20  | 10.  | 100. | 1000. |   1   |       |
    +-------+-------+-------+------+------+-------+-------+-------+
    """
    type = 'AEROS'
    _field_map = {
        1: 'acsid', 2:'rcsid', 3:'cRef', 4:'bRef', 5:'Sref',
        6:'symXZ', 7:'symXY',
    }
    _properties = ['is_anti_symmetric_xy', 'is_anti_symmetric_xz', 'is_symmetric_xy', 'is_symmetric_xz']

    @classmethod
    def _init_from_empty(cls):
        cref = 1.
        bref = 1.
        sref = 1.
        return AEROS(cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment='')

    def __init__(self, cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment=''):
        """
        Creates an AEROS card

        Parameters
        ----------
        cref : float
            the aerodynamic chord
        bref : float
            the wing span
            for a half model, this should be the full span
            for a full model, this should be the full span
        sref : float
            the wing area
            for a half model, this should be the half area
            for a full model, this should be the full area
        acsid : int; default=0
            aerodyanmic coordinate system
            defines the direction of the wind
        rcsid : int; default=0
            coordinate system for rigid body motions
        sym_xz : int; default=0
            xz symmetry flag (+1=symmetry; -1=antisymmetric)
        sym_xy : int; default=0
            xy symmetry flag (+1=symmetry; -1=antisymmetric)
        comment : str; default=''
            a comment for the card

        """
        Aero.__init__(self)
        if comment:
            self.comment = comment

        #: Aerodynamic coordinate system identification.
        self.acsid = acsid

        #: Reference coordinate system identification for rigid body motions.
        self.rcsid = rcsid

        #: Reference chord length
        self.cref = cref

        #: Reference span
        self.bref = bref

        #: Reference wing area
        self.sref = sref

        #: Symmetry key for the aero coordinate x-z plane. See Remark 6.
        #: (Integer = +1 for symmetry, 0 for no symmetry, and -1 for antisymmetry;
        #: Default = 0)
        self.sym_xz = sym_xz

        #: The symmetry key for the aero coordinate x-y plane can be used to
        #: simulate ground effects. (Integer = +1 for antisymmetry, 0 for no
        #: symmetry, and -1 for symmetry; Default = 0)
        self.sym_xy = sym_xy
        if self.acsid is None:
            self.acsid = 0
        if self.rcsid is None:
            self.rcsid = 0
        if self.sym_xz is None:
            self.sym_xz = 0
        if self.sym_xy is None:
            self.sym_xy = 0
        self.acsid_ref = None
        self.rcsid_ref = None

    def Acsid(self):
        """air velocity defined as moving into the +x direction"""
        try:
            return self.acsid_ref.cid
        except AttributeError:
            return self.acsid

    def Rcsid(self):
        """rigid body coordinate system"""
        try:
            return self.rcsid_ref.cid
        except AttributeError:
            return self.rcsid

    def validate(self):
        msg = ''
        if not isinstance(self.acsid, integer_types):
            msg += 'acsid=%s must be an integer; type=%s\n' % (self.acsid, type(self.acsid))
        if not isinstance(self.rcsid, integer_types):
            msg += 'rcsid=%s must be an integer; type=%s\n' % (self.rcsid, type(self.rcsid))
        if not isinstance(self.cref, float):
            msg += 'cref=%s must be an float; type=%s\n' % (self.cref, type(self.cref))
        if not isinstance(self.bref, float):
            msg += 'bref=%s must be an float; type=%s\n' % (self.bref, type(self.bref))
        if not isinstance(self.sref, float):
            msg += 'sref=%s must be an float; type=%s\n' % (self.sref, type(self.sref))
        if not isinstance(self.sym_xz, integer_types):
            msg += 'sym_xz=%s must be an integer; type=%s\n' % (self.sym_xz, type(self.sym_xz))
        if not isinstance(self.sym_xy, integer_types):
            msg += 'sym_xy=%s must be an integer; type=%s\n' % (self.sym_xy, type(self.sym_xy))
        if msg:
            raise TypeError('There are errors on the AEROS card:\n%s%s' % (msg, self))
        assert self.sym_xz in [-1, 0, 1], self.get_stats()
        assert self.sym_xy in [-1, 0, 1], self.get_stats()

    def cross_reference(self, model: BDF) -> None:
        """
        Cross refernece aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.

        """
        msg = ', which is required by AEROS'
        self.acsid_ref = model.Coord(self.acsid, msg=msg)
        self.rcsid_ref = model.Coord(self.rcsid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Safe cross refernece aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.

        """
        msg = ', which is required by AEROS'
        self.acsid_ref = model.safe_coord(self.acsid, None, xref_errors, msg=msg)
        self.rcsid_ref = model.safe_coord(self.rcsid, None, xref_errors, msg=msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AEROS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        acsid = integer_or_blank(card, 1, 'acsid', 0)
        rcsid = integer_or_blank(card, 2, 'rcsid', 0)
        cref = double(card, 3, 'cRef')
        bref = double(card, 4, 'bRef')
        sref = double(card, 5, 'Sref')
        sym_xz = integer_or_blank(card, 6, 'sym_xz', 0)
        sym_xy = integer_or_blank(card, 7, 'sym_xy', 0)
        assert len(card) <= 8, 'len(AEROS card) = %i\ncard=%s' % (len(card), card)
        return AEROS(cref, bref, sref, acsid, rcsid, sym_xz, sym_xy,
                     comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        acsid = data[0]
        rcsid = data[1]
        cref = data[2]
        bref = data[3]
        sref = data[4]
        sym_xz = data[5]
        sym_xy = data[6]
        assert len(data) == 7, 'data = %s' % data
        return AEROS(cref, bref, sref, acsid, rcsid, sym_xz, sym_xy,
                     comment=comment)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.acsid_ref = None
        self.rcsid_ref = None

    def update(self, maps):
        """
        maps = {
            'coord' : cid_map,
        }

        """
        cid_map = maps['coord']
        self.acsid = cid_map[self.acsid]
        self.rcsid = cid_map[self.rcsid]

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['AEROS', self.Acsid(), self.Rcsid(), self.cref,
                       self.bref, self.sref, self.sym_xz, self.sym_xy]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[varies]
          the fields that define the card

        """
        sym_xz = set_blank_if_default(self.sym_xz, 0)
        sym_xy = set_blank_if_default(self.sym_xy, 0)
        list_fields = ['AEROS', self.Acsid(), self.Rcsid(), self.cref,
                       self.bref, self.sref, sym_xz, sym_xy]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AESTAT(BaseCard):
    """
    Specifies rigid body motions to be used as trim variables in static
    aeroelasticity.

    +--------+------+--------+
    |    1   |   2  |    3   |
    +========+======+========+
    | AESTAT |  ID  | LABEL  |
    +--------+------+--------+
    | AESTAT | 5001 | ANGLEA |
    +--------+------+--------+
    """
    type = 'AESTAT'

    _field_map = {
        1: 'id', 2:'label',
    }
    @classmethod
    def _init_from_empty(cls):
        return AESTAT(1, 'name', comment='')

    def __init__(self, aestat_id, label, comment=''):
        """
        Creates an AESTAT card, which is a variable to be used in a TRIM analysis

        Parameters
        ----------
        id : int
            unique id
        label : str
            name for the id
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.aestat_id = aestat_id
        self.label = label

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AESTAT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        aestat_id = integer(card, 1, 'ID')
        label = string(card, 2, 'label')
        assert len(card) <= 3, 'len(AESTAT card) = %i\ncard=%s' % (len(card), card)
        return AESTAT(aestat_id, label, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        aestat_id = data[0]
        label = data[1]
        assert len(data) == 2, 'data = %s' % data
        return AESTAT(aestat_id, label, comment=comment)

    #def cross_reference(self, model: BDF) -> None:
        #pass

    #def uncross_reference(self) -> None:
        #pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/str]
            the fields that define the card

        """
        list_fields = ['AESTAT', self.aestat_id, self.label]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CSSCHD(Aero):
    """
    Defines a scheduled control surface deflection as a function of
    Mach number and angle of attack.

    +--------+-----+-------+--------+-------+-------+
    |    1   |  2  |   3   |   4    |   5   |   6   |
    +========+=====+=======+========+=======+=======+
    | CSSCHD | SlD | AESID | LALPHA | LMACH | LSCHD |
    +--------+-----+-------+--------+-------+-------+
    | CSSCHD |  5  |  50   |   12   |   15  |   25  |
    +--------+-----+-------+--------+-------+-------+
    """
    type = 'CSSCHD'
    _field_map = {
        1: 'sid', 2:'aesid', 3:'lalpha', 4:'lmach', 5:'lschd',
    }
    _properties = ['is_anti_symmetric_xy', 'is_anti_symmetric_xz', 'is_symmetric_xy', 'is_symmetric_xz'] ## TODO: remove these

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        aesid = 0
        lschd = 2
        return CSSCHD(sid, aesid, lschd, lalpha=None, lmach=None, comment='')

    def __init__(self, sid, aesid, lschd, lalpha=None, lmach=None, comment=''):
        """
        Creates an CSSCHD card, which defines a specified control surface
        deflection as a function of Mach and alpha (used in SOL 144/146).

        Parameters
        ----------
        sid : int
            the unique id
        aesid : int
            the control surface (AESURF) id
        lalpha : int; default=None
            the angle of attack profile (AEFACT) id
        lmach : int; default=None
            the mach profile (AEFACT) id
        lschd : int; default=None
            the control surface deflection profile (AEFACT) id
        comment : str; default=''
            a comment for the card

        """
        Aero.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.aesid = aesid
        self.lalpha = lalpha
        self.lmach = lmach
        self.lschd = lschd
        self.aesid_ref = None
        self.lalpha_ref = None
        self.lmach_ref = None
        self.lschd_ref = None

    def validate(self):
        if not(self.lalpha is None or isinstance(self.lalpha, integer_types)):
            raise TypeError('lalpha=%r must be an int or None' % self.lalpha)

        if not(self.lmach is None or isinstance(self.lmach, integer_types)):
            raise TypeError('lmach=%r must be an int or None' % self.lmach)

        if self.lalpha is None and self.lmach is None:
            msg = ('CSSCHD sid=%s; lalpha and lmach are both None'
                   ' (one must be an integer (AEFACT)\n%s' % (self.sid, str(self)))
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CSSCHD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        aesid = integer(card, 2, 'aesid')             # AESURF
        lalpha = integer_or_blank(card, 3, 'lAlpha')  # AEFACT
        lmach = integer_or_blank(card, 4, 'lMach')    # AEFACT
        lschd = integer(card, 5, 'lSchd')             # AEFACT
        assert len(card) <= 6, 'len(CSSCHD card) = %i\ncard=%s' % (len(card), card)
        return CSSCHD(sid, aesid, lalpha, lmach, lschd, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        aesid = data[1]   # AESURF
        lalpha = data[2]  # AEFACT
        lmach = data[3]   # AEFACT
        lschd = data[4]   # AEFACT
        return CSSCHD(sid, aesid, lalpha, lmach, lschd, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CSSCHD sid=%s' % self.sid
        self.aesid_ref = model.AESurf(self.aesid, msg=msg)
        self.lalpha_ref = model.AEFact(self.lalpha, msg=msg)
        self.lmach_ref = model.AEFact(self.lmach, msg=msg)
        self.lschd_ref = model.AEFact(self.lschd, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by CSSCHD sid=%s' % self.sid
        try:
            self.aesid_ref = model.AESurf(self.aesid, msg=msg)
        except KeyError:
            pass

        self.lalpha_ref = model.safe_aefact(self.lalpha, self.sid, xref_errors, msg=msg)
        self.lmach_ref = model.safe_aefact(self.lmach, self.sid, xref_errors, msg=msg)
        self.lschd_ref = model.safe_aefact(self.lschd, self.sid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.aesid = self.AESid()
        self.lalpha = self.LAlpha()
        self.lmach = self.LMach()
        self.lschd = self.LSchd()
        self.aesid_ref = None
        self.lalpha_ref = None
        self.lmach_ref = None
        self.lschd_ref = None

    def AESid(self):
        if self.aesid_ref is not None:
            return self.aesid_ref.aesid
        return self.aesid

    def LAlpha(self):
        if self.lalpha_ref is not None:
            return self.lalpha_ref.sid
        return self.lalpha

    def LMach(self):
        if self.lmach_ref is not None:
            return self.lmach_ref.sid
        return self.lmach

    def LSchd(self):
        if self.lschd_ref is not None:
            return self.lschd_ref.sid
        return self.lschd

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['CSSCHD', self.sid, self.AESid(), self.LAlpha(),
                       self.LMach(), self.LSchd()]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class DIVERG(BaseCard):
    """
    +--------+-----+--------+----+----+----+----+----+----+
    |   1    |  2  |   3    | 4  | 5  | 6  | 7  | 8  | 9  |
    +========+=====+========+====+====+====+====+====+====+
    | DIVERG | SID | NROOT  | M1 | M2 | M3 | M4 | M5 | M6 |
    +--------+-----+--------+----+----+----+----+----+----+
    |        |  M7 |  etc.  |    |    |    |    |    |    |
    +--------+-----+--------+----+----+----+----+----+----+

    Attributes
    ----------
    sid : int
        The name.
    nroots : int
        the number of roots
    machs : List[float, ..., float]
        list of Mach numbers

    """
    type = 'DIVERG'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        nroots = 10
        machs = [0.5, 0.75]
        return DIVERG(sid, nroots, machs, comment='')

    def __init__(self, sid, nroots, machs, comment=''):
        """
        Creates an DIVERG card, which is used in divergence
        analysis (SOL 144).

        Parameters
        ----------
        sid : int
            The name
        nroots : int
            the number of roots
        machs : List[float, ..., float]
            list of Mach numbers
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.nroots = nroots
        self.machs = machs

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DIVERG card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        nroots = integer(card, 2, 'nroot')
        j = 1
        machs = []
        for i in range(3, len(card)):
            mach = double(card, i, 'Mach_%i' % j)
            machs.append(mach)
            j += 1
        return DIVERG(sid, nroots, machs, comment=comment)

    #def cross_reference(self, model: BDF) -> None:
        #pass

    #def uncross_reference(self) -> None:
        #pass

    def raw_fields(self):
        list_fields = ['DIVERG', self.sid, self.nroots] + list(self.machs)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class TRIM(BaseCard):
    """
    Specifies constraints for aeroelastic trim variables.

    +------+--------+------+--------+--------+-----+--------+-----+----------+
    |   1  |   2    |   3  |    4   |    5   |  6  |    7   |  8  |     9    |
    +======+========+======+========+========+=====+========+=====+==========+
    | TRIM |   ID   | MACH |    Q   | LABEL1 | UX1 | LABEL2 | UX2 | IS_RIGID |
    +------+--------+------+--------+--------+-----+--------+-----+----------+
    |      | LABEL3 |  UX3 | LABEL4 |   UX4  | ... |        |     |          |
    +------+--------+------+--------+--------+-----+--------+-----+----------+
    """
    type = 'TRIM'
    _field_map = {
        1: 'sid', 2:'mach', 3:'q', 8:'aeqr',
    }

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the TRIM card

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : varies
            the value for the appropriate field

        """
        ni = 4
        for i in range(len(self.labels)):
            if n == ni:
                value = self.labels[i]
                return value

            if n + 1 == ni:
                value = self.uxs[i]
                return value

            #list_fields += [label, ux]
            if i == 1:
                #list_fields += [self.aeqr]
                ni += 1
        raise KeyError('Field %r is an invalid TRIM entry.' % n)

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the TRIM card

        Parameters
        ----------
        n : int
            the field number to update
        value : varies
            the value for the appropriate field

        """
        ni = 4
        for i in range(len(self.labels)):
            if n == ni:
                self.labels[i] = value
                return

            if n + 1 == ni:
                self.uxs[i] = value
                return

            if i == 1:
                ni += 1
        raise KeyError('Field %r=%r is an invalid TRIM entry.' % (n, value))

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        mach = 0.6
        q = 300.
        labels = ['ALPHA']
        uxs = [1.0]
        return TRIM(sid, mach, q, labels, uxs, aeqr=1.0, comment='')

    def __init__(self, sid, mach, q, labels, uxs, aeqr=1.0, comment=''):
        """
        Creates a TRIM card for a static aero (144) analysis.

        Parameters
        ----------
        sid : int
            the trim id; referenced by the Case Control TRIM field
        mach : float
            the mach number
        q : float
            dynamic pressure
        labels : List[str]
            names of the fixed variables
        uxs : List[float]
            values corresponding to labels
        aeqr : float
            0.0 : rigid trim analysis
            1.0 : elastic trim analysis (default)
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        #: Trim set identification number. (Integer > 0)
        self.sid = sid
        #: Mach number. (Real > 0.0 and != 1.0)
        self.mach = mach
        #: Dynamic pressure. (Real > 0.0)
        self.q = q

        #: The label identifying aerodynamic trim variables defined on an
        #: AESTAT or AESURF entry.
        self.labels = labels

        #: The magnitude of the aerodynamic extra point degree-of-freedom.
        #: (Real)
        self.uxs = uxs

        #: Flag to request a rigid trim analysis (Real > 0.0 and < 1.0;
        #: Default = 1.0. A value of 0.0 provides a rigid trim analysis.
        self.aeqr = aeqr

    def validate(self):
        assert self.mach >= 0.0, 'mach = %r' % self.mach
        assert self.mach != 1.0, 'mach = %r' % self.mach
        assert self.q > 0.0, 'q=%s' % self.q
        if len(set(self.labels)) != len(self.labels):
            msg = 'not all labels are unique; labels=%s' % str(self.labels)
            raise RuntimeError(msg)
        if len(self.labels) != len(self.uxs):
            msg = 'nlabels=%s != nux=%s; labels=%s uxs=%s' % (
                len(self.labels), len(self.uxs), str(self.labels), str(self.uxs))
            raise RuntimeError(msg)

    def verify_trim(self, suport, suport1, aestats, aeparms, aelinks, aesurf, xref=True):
        """
        Magic function that makes TRIM cards not frustrating.

        .. warning ::  This probably gets AELINKs/AEPARMs/AESURFSs wrong.

        **The TRIM equality**
        ndelta = (naestat + naesurf + naeparm) - (
               - (ntrim + ntrim_aesurf? + naelink + nsuport_dofs + nsuport1_dofs)
        ndelta = 0
        ntrim_aesurf is not included, but it might exist...

        **Steps to a TRIM analysis**
        1.  Define the number of independent control surfaces (naesurf)
            Make an AESURF for each.  Dual link the AESURFs if you can
            to avoid needing an AELINK (e.g., +roll is left aileron down,
            right aileron up).
            Horizontal Tail : name it DPITCH
            Vertical Tail   : name it DYAW
            Aileron         : name it DROLL
        2.  Create AELINKs if necessary.
        3.  Add the AESTAT variables.  Include one for each DOF the
            aircraft can move in the frame of the model
            (e.g., half/full model).
                Half model (2.5g pullup, abrupt pitch):
                  - 2d pitch/plunge, 1 control : URDD3, URDD5, PITCH, ANGLEA
                Full model (2.5g pullup, abrupt pitch):
                  - 3d pitch/plunge, 3 control : URDD3, URDD5, PITCH, ANGLEA, YAW (???)
        4.  Add the TRIM card to lock the variables that could theoretically move
            in the plane of the analysis that are known.
                Half model:
                   2.5g pullup   : lock URDD3=2.5, URDD5=0, PITCH=0
                                   solve for ANGLEA, DPITCH
                                   use DPITCH
                   abrupt pitch  : lock URDD3=1.0, URDD5=0, ANGLEA=5
                                   solve for PITCH, DPITCH
                                   use DPITCH
                Full model:
                   2.5g pullup   : lock URDD3=2.5, URDD4=0, URDD5=0,  PITCH=0, YAW=0,
                                   lock SIDES=0,  ROLL=0
                                   solve for ANGLEA, DPITCH
                                   use DPITCH, DYAW, DROLL
                                   TODO: probably wrong
                   30 degree yaw : lock URDD3=1.0, URDD4=0, ANGLEA=5, PITCH=0, YAW=30,
                                   lock DPITCH=0, ROLL=0
                                   solve for SIDES, URDD5
                                   use DPITCH, DYAW, DROLL
                                   TODO: probably wrong

        5.  Note that we could have simplified our full model AESTAT/TRIM
            cards (they can be the same as for a half model), but we'd
            like to be able to do multiple load cases in the same deck.

        6.  Add some SUPORT/SUPORT1 DOFs to ignore non-relevant motion in
            certain DOFs (e.g., z-motion).  Add enough to satisfy the TRIM
            equality.

        **Doesn't Consider**
         - AELINK
         - AEPARM
         - AESURFS

        +------------------------------------------------+
        |                 Default AESTATs                |
        +--------+---------+-----------------------------+
        | ANGLEA | ur (R2) | Angle of Attack             |
        | YAW    | ur (R3) | Yaw Rate                    |
        | SIDES  | ur (R3) | Angle of Sideslip           |
        +--------+---------+-----------------------------+
        | ROLL   | ůr (R1) | Roll Rate                   |
        | PITCH  | ůr (R2) | Pitch Rate                  |
        +--------+---------+-----------------------------+
        | URDD1  | ür (T1) | Longitudinal (See Remark 3) |
        | URDD2  | ür (T2) | Lateral                     |
        | URDD3  | ür (T3) | Vertical                    |
        | URDD4  | ür (R1) | Roll                        |
        | URDD5  | ür (R2) | Pitch                       |
        | URDD6  | ür (R3) | Yaw                         |
        +--------+---------+-----------------------------+
        """
        if xref:
            nsuport_dofs = 0
            nsuport1_dofs = 0
            suport_dofs = set()
            assert isinstance(suport, list), type(suport)
            for suporti in suport:
                for nid, cs in zip(suporti.node_ids, suporti.Cs):
                    for ci in cs:
                        #print('  SUPORT: nid=%r C=%r' % (nid, ci))
                        dof = (nid, ci)
                        if dof in suport_dofs:
                            msg = 'Duplicate DOF\n  dof=%s suport_dofs=%s' % (
                                str(dof), str(suport_dofs))
                            raise RuntimeError(msg)
                        suport_dofs.add(dof)
                        nsuport_dofs += 1

            suport_dof_msg2 = ''
            if suport1:
                #unused_conid = suport1.conid
                nids = suport1.node_ids
                suport_dof_msg = ''
                for nid, components in zip(nids, suport1.Cs):
                    for componenti in components:
                        dof = (nid, componenti)
                        suport_dof_msg += '    (%s, %s)\n' % (nid, componenti)
                        if dof in suport_dofs:
                            msg = 'dof=%s suport_dofs=%s' % (str(dof), str(suport_dofs))
                            raise RuntimeError(msg)
                        suport_dofs.add(dof)
                        nsuport1_dofs += 1
                suport_dof_msg2 = '\nsuport_dofs (nid, comp):\n%s\n' % suport_dof_msg.rstrip(',')

            aesurf_names = [aesurfi.label for aesurfi in aesurf.values()]
            aestat_labels = [aestat.label for aestat in aestats.values()]
            aeparm_labels = [aeparm.label for aeparm in aeparms.values()]
            naestat = len(aestat_labels)
            ntrim = len(self.labels)
            trim_aesurf_common = list(set(self.labels).intersection(set(aesurf_names)))
            trim_aesurf_common.sort()
            ntrim_aesurfs = len(trim_aesurf_common)
            naesurf = len(aesurf_names)
            naeparm = len(aeparm_labels)

            aelinksi = []
            if 0 in aelinks:
                aelinksi += [aelink.label for aelink in aelinks[0]]
            #if 'ALWAYS' in aelinks:
                #aelinksi += [aelink.label for aelink in aelinks['ALWAYS']]

            if self.sid in aelinks:
                aelinksi += [aelink.label for aelink in aelinks[self.sid]]
            naelink = len(aelinksi)


            ntrim_aesurf = 0
            labels = aestat_labels + aesurf_names + aeparm_labels
            msg = ''
            for label in self.labels:
                if label not in labels:
                    msg += 'TRIM label=%r is not defined\n' % label

                if label in aesurf_names:
                    #print('AESTAT/AESURF label = %r' % label)
                    ntrim_aesurf += 1
            if msg:
                msg += '\n aestat_labels=%s\n aeparm_labels=%s\n aesurf_names=%s\n%s' % (
                    aestat_labels, aeparm_labels, aesurf_names, str(self))
                raise RuntimeError(msg)

            # TODO: this doesn't work for multiple subcases
            #ntotal_suport_dofs = nsuport_dofs, nsuport1_dofs
            #ndelta = ntrim - nsuport_dofs - nsuport1_dofs - naesurf
            #if ndelta != 0:
                #msg = 'ntrim - nsuport_dofs - nsuport1_dofs - naesurf = ndelta = %s; ndelta != 0\n' % ndelta
                #msg += 'ntrim=%s nsuport_dofs=%s nsuport1_dofs=%s naesurfs=%s' % (
                    #ntrim, nsuport_dofs, nsuport1_dofs, naesurf)
                #raise RuntimeError(msg)

            #ndelta = (naestat + naesurf + naeparm + ntrim_aesurf) - (ntrim + naelink + nsuport_dofs + nsuport1_dofs)
            #if ndelta != 0:
                #msg = (
                    #'(naestat + naesurf + naeparm + ntrim_aesurf) - '
                    #'(ntrim + naelink + nsuport_dofs + nsuport1_dofs) = ndelta = %s; ndelta != 0\n'
                    #'naestat=%s naesurf=%s naeparm=%s ntrim_aesurfs=%s\n'
                    #'ntrim=%s naelink=%s nsuport_dofs=%s nsuport1_dofs=%s' % (
                        #ndelta,
                        #naestat, naesurf, naeparms, ntrim_aesurf,
                        #ntrim, naelink, nsuport_dofs, nsuport1_dofs))

            nplus = (naestat + naesurf + naeparm)
            nminus = ntrim + naelink + nsuport_dofs + nsuport1_dofs

            ndelta = nplus - nminus + 0*2*ntrim_aesurfs
            if ndelta != 0:
                #msg = (
                    #'(naestat + naesurf + naeparm) - (ntrim + ntrim_aesurf? + naelink + '
                    #'nsuport_dofs + nsuport1_dofs) = ndelta = %s; ndelta != 0\n'
                    #'naestat=%s naesurf=%s naeparm=%s ntrim=%s ntrim_aesurf=%s '
                    #'naelink=%s nsuport_dofs=%s nsuport1_dofs=%s\n' % (
                        #ndelta,
                        #naestat, naesurf, naeparm, ntrim, ntrim_aesurf,
                        #naelink, nsuport_dofs, nsuport1_dofs)
                #)
                msg = (
                    'Invalid trim state (ndelta != 0):\n'
                    '   (naestat + naesurf + naeparm + 0*2*ntrim_aesurf?) = (%s + %s + %s + 0*2*%s) = %s\n'
                    ' - (ntrim + naelink + nsuport_dofs + nsuport1_dofs) = (%s + %s + %s + %s) = %s\n'
                    '===================================================================\n'
                    '  ndelta = %s\n\n'
                    'Summary\n'
                    '-------\n'
                    '  +naestat = %s; %s\n'
                    '  +naesurf = %s; %s\n'
                    '  +naeparm = %s; %s\n'
                    '  +0*2*ntrim_aesurf? = %s -> 0; %s\n'
                    '  -ntrim = %s; %s\n'
                    '  -naelink = %s; %s\n'
                    '  -nsuport_dofs = %s\n'
                    '  -nsuport1_dofs = %s\n'
                    '%s\n\n' % (
                        naestat, naesurf, naeparm, ntrim_aesurf, nplus,
                        ntrim, naelink, nsuport_dofs, nsuport1_dofs, nminus,

                        ndelta,
                        naestat, aestat_labels,
                        naesurf, aesurf_names,
                        naeparm, aeparm_labels,
                        2*ntrim_aesurf, trim_aesurf_common,
                        ntrim, self.labels,
                        naelink, aelinksi,
                        nsuport_dofs, nsuport1_dofs, suport_dof_msg2)
                )
                msg += str(self)
                raise RuntimeError(msg)

    def cross_reference(self, model: BDF) -> None:
        pass
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        mach = double(card, 2, 'mach')
        q = double(card, 3, 'q')
        labels = []
        uxs = []

        label = string_or_blank(card, 4, 'label1')
        if label:
            ux = double(card, 5, 'ux1')
            uxs.append(ux)
            labels.append(label)

        label = string_or_blank(card, 6, 'label2')
        if label:
            ux = double(card, 7, 'ux1')
            uxs.append(ux)
            labels.append(label)
        aeqr = double_or_blank(card, 8, 'aeqr', 1.0)

        i = 9
        n = 3
        while i < len(card):
            label = string(card, i, 'label%i' % n)
            ux = double(card, i + 1, 'ux%i' % n)
            labels.append(label)
            uxs.append(ux)
            i += 2
            n += 1
        return TRIM(sid, mach, q, labels, uxs, aeqr, comment=comment)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['TRIM', self.sid, self.mach, self.q]
        nlabels = len(self.labels)
        assert nlabels > 0, self.labels
        for (i, label, ux) in zip(count(), self.labels, self.uxs):
            list_fields += [label, ux]
            if i == 1:
                list_fields += [self.aeqr]
        if nlabels == 1:
            list_fields += [None, None, self.aeqr]
        return list_fields

    def repr_fields(self):
        # fixes a Nastran bug
        aeqr = set_blank_if_default(self.aeqr, 1.0)

        list_fields = self.raw_fields()
        list_fields[8] = aeqr
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class TRIM2(TRIM):
    """
    Defines the state of the aerodynamic extra points for a trim analysis.
    All undefined extra points will be set to zero.

    +-------+--------+------+--------+--------+-----+--------+-----+----------+
    |   1   |   2    |   3  |    4   |    5   |  6  |    7   |  8  |     9    |
    +=======+========+======+========+========+=====+========+=====+==========+
    | TRIM2 |   ID   | MACH |    Q   |        |     |        |     | IS_RIGID |
    +-------+--------+------+--------+--------+-----+--------+-----+----------+
    |       | LABEL1 |  UX1 | LABEL2 |   UX2  | ... |        |     |          |
    +-------+--------+------+--------+--------+-----+--------+-----+----------+
    """
    type = 'TRIM2'
    _field_map = {
        1: 'sid', 2:'mach', 3:'q', 8:'aeqr',
    }

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        mach = 0.6
        q = 300.
        labels = ['ALPHA']
        uxs = [1.0]
        return TRIM2(sid, mach, q, labels, uxs, aeqr=1.0, comment='')

    def __init__(self, sid, mach, q, labels, uxs, aeqr=1.0, comment=''):
        TRIM.__init__(self, sid, mach, q, labels, uxs, aeqr=aeqr, comment=comment)

    @classmethod
    def add_card(cls, card, comment=''):  # TODO: not done...
        """
        Adds a TRIM2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        mach = double(card, 2, 'mach')
        q = double(card, 3, 'q')
        aeqr = double_or_blank(card, 8, 'aeqr', 1.0)

        i = 9
        n = 3
        labels = []
        uxs = []
        while i < len(card):
            label = string(card, i, 'label%i' % n)
            ux = double(card, i + 1, 'ux%i' % n)
            labels.append(label)
            uxs.append(ux)
            i += 2
        return TRIM2(sid, mach, q, labels, uxs, aeqr, comment=comment)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['TRIM2', self.sid, self.mach, self.q, None, None, None, None, self.aeqr]
        nlabels = len(self.labels)
        assert nlabels > 0, self.labels
        for label, ux in zip(self.labels, self.uxs):
            list_fields += [label, ux]
        return list_fields

    def repr_fields(self):
        # fixes a Nastran bug
        aeqr = set_blank_if_default(self.aeqr, 1.0)

        list_fields = self.raw_fields()
        list_fields[8] = aeqr
        return list_fields
