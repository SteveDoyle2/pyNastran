"""
All axisymmetric shell elements are defined in this file.  This includes:
 * AXIC
 * PRESAX

"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import (
    set_blank_if_default,
)
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, double, blank,
    string, string_or_blank, integer_or_string, fields,
)
from pyNastran.bdf.cards.base_card import BaseCard, Element

from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.cards.thermal.loads import ThermalLoad
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double


class AXIC(BaseCard):
    type = 'AXIC'
    def __init__(self, nharmonics, comment=''):
        if comment:
            self.comment = comment
        self.nharmonics = nharmonics

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a AXIC card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nharmonics = integer(card, 1, 'nharmonics')
        assert len(card) == 2, 'len(AXIC card) = %i\ncard=%s' % (len(card), card)
        return AXIC(nharmonics, comment=comment)

    def raw_fields(self):
        list_fields = ['AXIC', self.nharmonics]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        msg = self.comment + print_card_8(card)
        return msg

class AXIF(BaseCard):
    """
    AXIF CID G  DRHO DB NOSYM F
         N1  N2 N3   N4  N5   -etc.-
    """
    type = 'AXIF'
    def __init__(self, cid, g, drho, db, no_sym, f, n, comment=''):
        """
        cid : int
            Fluid coordinate system identification number. (Integer > 0)
        G : float
            Value of gravity for fluid elements in the axial direction. (Real)
        drho : float
            Default mass density for fluid elements. (Real > 0.0 or blank)
        db : float
            Default bulk modulus for fluid elements.
        no_sym : str
            Request for nonsymmetric (sine) terms of series.
            {YES, NO}
        F : str; default=None
            Flag specifying harmonics. (Blank if harmonic is specified, or Character:
            'NONE')
        Ni : List[int]
            Harmonic numbers for the solution, represented by an increasing
            sequence of integers. On continuation entries, without the 'THRU'
            option, blank fields are ignored. 'THRU' implies all numbers including
            upper and lower harmonics. (0 < Integer < 100, or Character: 'THRU',
            'STEP' or blank)
        """
        if comment:
            self.comment = comment
        self.cid = cid
        self.g = g
        self.drho = drho
        self.db = db
        self.no_sym = no_sym
        self.f = f
        self.n = n

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a AXIF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        cid = integer(card, 1, 'cid')
        g = double(card, 2, 'g')
        drho = double(card, 3, 'drho')
        db = double_or_blank(card, 4, 'db')
        no_sym = string(card, 5, 'no_sym')
        f = string_or_blank(card, 6, 'f')
        n = fields(integer_or_string, card, 'n', i=9, j=len(card))
        #cid : int
        #G : float
        #drho : float
        #db : float
        #no_sym : str
        #F : str
        #Ni : List[int]


        assert len(card) >= 7, 'len(AXIF card) = %i\ncard=%s' % (len(card), card)
        return AXIF(cid, g, drho, db, no_sym, f, n, comment=comment)

    def raw_fields(self):
        list_fields = ['AXIF', self.cid, self.g, self.drho, self.db, self.no_sym, self.f] + self.n
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        msg = self.comment + print_card_8(card)
        return msg

class POINTAX(BaseCard):
    """
    +---------+----+-----+------+
    |    1    |  2 |  3  |  4   |
    +=========+====+=====+======+
    | POINTAX | ID | RID | PHI  |
    +---------+----+-----+------+
    | POINTAX | 2  |  3  | 30.0 |
    +---------+----+-----+------+
    """
    type = 'POINTAX'
    def __init__(self, nid, ringax, phi, comment=''):
        """
        Creates a POINTAX card

        Parameters
        ----------
        nid : int
            Point identification number.
        ringax : int
            Identification number of a RINGAX entry.
        phi : float
            Azimuthal angle in degrees.
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        self.nid = nid
        self.ringax = ringax
        self.phi = phi

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RINGAX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nid = integer(card, 1, 'nid')
        ringax = integer(card, 2, 'ringax')
        phi = double(card, 3, 'phi')
        assert len(card) <= 4, 'len(POINTAX card) = %i\ncard=%s' % (len(card), card)
        return POINTAX(nid, ringax, phi, comment=comment)

    def raw_fields(self):
        list_fields = ['POINTAX', self.nid, self.ringax, self.phi]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        msg = self.comment + print_card_8(card)
        return msg


class RINGFL(BaseCard):
    type = 'RINGFL'

    #: allows the get_field method and update_field methods to be used
    #_field_map = {1: 'mid', 3:'R', 4:'z', 7:'ps'}

    def __init__(self, ringfl, xa, xb, comment=''):  # this card has missing fields
        # type: (int, float, float, Optional[str], str) -> None
        """
        Creates the RINGFL card
        """
        #Ring.__init__(self)
        if comment:
            self.comment = comment

        self.ringfl = ringfl
        self.xa = xa
        self.xb = xb

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        """
        Adds a RINGAX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        ioffset = 2 * icard
        j = icard + 1
        ringfl = integer(card, 1+ioffset, 'ringfl_%i' % j)
        xa = double(card, 2+ioffset, 'xa_%i' % j)
        xb = double_or_blank(card, 3+ioffset, 'xa_%i' % j)
        assert len(card) <= 9, 'len(RINGFL card) = %i\ncard=%s' % (len(card), card)
        return RINGFL(ringfl, xa, xb, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a RINGAX card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        nid = data[0]
        R = data[1]
        z = data[2]
        ps = data[3]
        assert len(data) == 4, data
        return RINGAX(nid, R, z, ps, comment=comment)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[varies]
            the fields that define the card

        """
        list_fields = ['RINGFL', self.ringfl, self.xa, self.xb]
        return list_fields

    def write_card(self, size=8, is_double=False):
        # (int, bool) -> str
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class RINGAX(BaseCard):
    """
    Defines a ring for conical shell problems.

    +-------+-----+-----+-----+----+-----+-----+------+
    |   1   |  2  |  3  |  4  |  5 |  6  |  7  |  8   |
    +=======+=====+=====+=====+====+=====+=====+======+
    |RINGAX | MID |     |  R  |  Z |     |     |  PS  |
    +-------+-----+-----+-----+----+-----+-----+------+
    """
    type = 'RINGAX'

    #: allows the get_field method and update_field methods to be used
    _field_map = {1: 'mid', 3:'R', 4:'z', 7:'ps'}

    def __init__(self, nid, R, z, ps=None, comment=''):  # this card has missing fields
        # type: (int, float, float, Optional[str], str) -> None
        """
        Creates the RINGAX card
        """
        #Ring.__init__(self)
        if comment:
            self.comment = comment

        #: Node ID
        self.nid = nid

        #: Radius
        self.R = R
        self.z = z

        #: local SPC constraint
        self.ps = ps

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RINGAX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nid = integer(card, 1, 'nid')
        blank(card, 2, 'blank')

        R = double(card, 3, 'R')
        z = double(card, 4, 'z')
        blank(card, 5, 'blank')
        blank(card, 6, 'blank')

        ps = integer_or_blank(card, 7, 'ps')
        assert len(card) <= 8, 'len(RINGAX card) = %i\ncard=%s' % (len(card), card)
        return RINGAX(nid, R, z, ps, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a RINGAX card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        nid = data[0]
        R = data[1]
        z = data[2]
        ps = data[3]
        assert len(data) == 4, data
        return RINGAX(nid, R, z, ps, comment=comment)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[varies]
            the fields that define the card

        """
        list_fields = ['RINGAX', self.nid, None, self.R, self.z, None,
                       None, self.ps]
        return list_fields

    def write_card(self, size=8, is_double=False):
        # (int, bool) -> str
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class CCONEAX(Element):
    """
    +---------+-----+-----+----+----+
    |    1    |  2  |  3  |  4 |  5 |
    +=========+=====+=====+====+====+
    | CCONEAX | EID | PID | N1 | N2 |
    +---------+-----+-----+----+----+
    """
    type = 'CCONEAX'
    _field_map = {
        1: 'eid', 2:'pid',
    }

    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 4:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, pid, rings, comment=''):
        """
        Creates a CCONEAX card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PCONEAX)
        nids : List[int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        Element.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        #self.prepare_node_ids(nids)
        self.rings = rings
        assert len(self.rings) == 2
        self.rings_ref = None
        self.pid_ref = None
        self.rings_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CCONEAX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2')]
        assert len(card) == 5, 'len(CROD card) = %i\ncard=%s' % (len(card), str(card))
        return CCONEAX(eid, pid, nids, comment=comment)

    @property
    def nodes(self):
        return self.rings

    @property
    def node_ids(self):
        return self.ring_ids

    @property
    def ring_ids(self):
        if self.rings_ref is not None:
            raise NotImplementedError(self.rings_ref)
        return self.rings

    def cross_reference(self, model):
        msg = ', which is required by CCONEAX eid=%s' % (self.eid)
        #self.rings_ref
        self.pid_ref = model.Property(self.pid, msg=msg)

    def uncross_reference(self):
        self.pid = self.Pid()
        self.pid_ref = None
        self.rings_ref = None

    def raw_fields(self):
        list_fields = ['CCONEAX', self.eid, self.Pid()] + self.ring_ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class PCONEAX(Property):
    """
    +---------+-------+--------+------+-------+-------+-------+-------+-------+
    |    1    |   2   |   3    |  4   |   5   |   6   |   7   |   8   |   9   |
    +=========+=======+========+======+=======+=======+=======+=======+=======+
    | PCONEAX |  ID   |  MID1  |  T1  | MID2  |   I   | MID3  |   T2  |  NSM  |
    +---------+-------+--------+------+-------+-------+-------+-------+-------+
    |         |  Z1   |   Z2   | PHIl | PHI2  | PHI3  | PHI4  | PHI5  | PHI6  |
    +---------+-------+--------+------+-------+-------+-------+-------+-------+
    |         | PHI7  |  PHI8  | PHI9 | PHI10 | PHI11 | PHI12 | PHI13 | PHI14 |
    +---------+-------+--------+------+-------+-------+-------+-------+-------+
    | PCONEAX |   2   |   4    | 1.0  |   6   | 16.3  |   8   |  2.1  |  0.5  |
    +---------+-------+--------+------+-------+-------+-------+-------+-------+
    |         | 0.001 | -0.002 | 23.6 | 42.9  |       |       |       |       |
    +---------+-------+--------+------+-------+-------+-------+-------+-------+
    """
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

    def __init__(self, pid, mid1, t1=None, mid2=0, i=None, mid3=None, t2=None,
                 nsm=None, z1=None, z2=None, phi=None, comment=''):
        """
        Creates a PCONEAX

        Parameters
        ----------
        pid : int
            PCONEAX property id for a CCONEAX.
        mid1 : int
            Membrane material id
        mid2 : int
            bending material id
        mid3 : int
            transverse shear material id
        t1 : float
            Membrane thickness. (Real > 0.0 if MID1 = 0)
        t2 : float
            Transverse shear thickness. (Real > 0.0 if MID3 = 0)
        I : float
            Moment of inertia per unit width.
        nsm : float
            Nonstructural mass per unit area.
        z1, z2 : float
            Fiber distances from the middle surface for stress recovery.
        phi : List[float]
            Azimuthal coordinates (in degrees) for stress recovery.
        comment : str; default=''
            a comment for the card

        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        if phi is None:
            phi = []
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
        self.mid1_ref = None
        self.mid2_ref = None
        self.mid3_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PCONEAX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
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

        mid3 = integer_or_blank(card, 6, 'mid3', 0)
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
        for icard in range(11, len(card)):
            phii = double(card, icard, 'phi%i' % icard)
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
        msg = ', which is required by PCONEAX=%s' %(self.pid)
        if self.mid1 > 0:
            self.mid1_ref = model.Material(self.mid1, msg=msg)
        if self.mid2 > 0:
            self.mid2_ref = model.Material(self.mid2, msg=msg)
        if self.mid3 > 0:
            self.mid3_ref = model.Material(self.mid3, msg=msg)

    def uncross_reference(self):
        self.mid1 = self.Mid1()
        self.mid2 = self.Mid2()
        self.mid3 = self.Mid3()
        self.mid1_ref = None
        self.mid2_ref = None
        self.mid3_ref = None

    @property
    def mid_ref(self):
        return self.mid1_ref

    def Mid1(self):
        if self.mid_ref is not None:
            return self.mid1_ref.mid
        return self.mid1

    def Mid2(self):
        if self.mid2_ref is not None:
            return self.mid2_ref.mid
        return self.mid2

    def Mid3(self):
        if self.mid3_ref is not None:
            return self.mid3_ref.mid
        return self.mid3

    def Mids(self):
        return [self.Mid1(), self.Mid2(), self.Mid3()]

    def raw_fields(self):
        list_fields = ['PCONEAX', self.pid, self.Mid1(), self.t1,
                       self.Mid2(), self.i, self.Mid3(), self.t2,
                       self.nsm, self.z1, self.z2] + self.phi
        return list_fields

    def repr_fields(self):
        #nsm = set_blank_if_default(self.nsm, 0.0)
        mid1 = set_blank_if_default(self.Mid1(), 0)
        mid2 = set_blank_if_default(self.Mid2(), 0)
        mid3 = set_blank_if_default(self.Mid3(), 0)
        i = set_blank_if_default(self.i, 0.0)
        t1 = set_blank_if_default(self.t1, 0.0)
        t2 = set_blank_if_default(self.t2, 0.0)
        list_fields = ['PCONEAX', self.pid, mid1, t1, mid2, i, mid3, t2,
                       self.nsm, self.z1, self.z2] + self.phi
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PRESAX(BaseCard):
    """
    +--------+-----+------+------+------+------+------+
    |   1    |  2  |   3  |  4   |  5   |   6  |   7  |
    +========+=====+======+======+======+======+======+
    | PRESAX | SID |   P  | RID1 | RID2 | PHI1 | PHI2 |
    +--------+-----+------+------+------+------+------+
    | PRESAX |  3  | 7.92 |  4   |   3  | 20.6 | 31.4 |
    +--------+-----+------+------+------+------+------+
    | PRESAX | 300 |  .1  |  2   |   1  | -90. | +90. |
    +--------+-----+------+------+------+------+------+
    """
    type = 'PRESAX'
    def __init__(self, sid, pressure, rid1, rid2, phi1=0., phi2=360., comment=''):
        if comment:
            self.comment = comment
        self.sid = sid
        self.pressure = pressure
        self.rid1 = rid1
        self.rid2 = rid2
        self.phi1 = phi1
        self.phi2 = phi2

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PRESAX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        pressure = double(card, 2, 'pressure')
        rid1 = integer(card, 3, 'rid1')
        rid2 = integer(card, 4, 'rid2')
        phi1 = double_or_blank(card, 5, 'phi1', 0.)
        phi2 = double_or_blank(card, 6, 'phi2', 360.)
        assert len(card) == 7, 'len(PRESAX card) = %i\ncard=%s' % (len(card), card)
        return PRESAX(sid, pressure, rid1, rid2, phi1, phi2, comment=comment)

    def cross_reference(self, model):
        pass

    def raw_fields(self):
        list_fields = ['PRESAX', self.sid, self.pressure, self.rid1, self.rid2,
                       self.phi1, self.phi2]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        msg = self.comment + print_card_8(card)
        return msg

class TEMPAX(ThermalLoad):
    """
    Defines temperature sets for conical shell problems.

    +--------+-----+------+-------+-------+-----+-------+------+----+
    |    1   |  2  |   3  |   4   |   5   |  6  |   7   |   8  |  9 |
    +========+=====+======+=======+=======+=====+=======+======+====+
    | TEMPAX | SID | RID1 |  PHI1 |  T1   | SID | RID2  | PHI2 | T2 |
    +--------+-----+------+-------+-------+-----+-------+------+----+
    | TEMPAX |  4  |   7  |  30.0 | 105.3 |     |       |      |    |
    +--------+-----+------+-------+-------+-----+-------+------+----+
    """
    type = 'TEMPAX'

    def __init__(self, sid, ring, phi, temperature, comment=''):
        """
        Creates a TEMPAX card

        Parameters
        ----------
        sid : int
            Load set identification number
        temperatures : dict[nid] : temperature
            nid : int
                node id
            temperature : float
                the nodal temperature
        comment : str; default=''
            a comment for the card

        """
        ThermalLoad.__init__(self)
        if comment:
            self.comment = comment
        #: Load set identification number. (Integer > 0)
        self.sid = sid

        self.ring = ring
        self.phi = phi
        self.temperature = temperature

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        """
        Adds a TEMPAX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            ???
        comment : str; default=''
            a comment for the card

        """
        istart = 1 + icard * 4
        sid = integer(card, istart, 'sid')
        ring = integer(card, istart + 1, 'ring' + str(icard))
        phi = double(card, istart + 2, 'phi' + str(icard))
        temperature = double(card, istart + 3, 'T' + str(icard))

        return TEMPAX(sid, ring, phi, temperature, comment=comment)

    def cross_reference(self, model):
        pass

    def safe_cross_reference(self, model, debug=True):
        pass

    def uncross_reference(self):
        pass

    def raw_fields(self):
        """Writes the TEMPAX card"""
        list_fields = ['TEMPAX', self.sid, self.ring, self.phi, self.temperature]
        return list_fields

    def repr_fields(self):
        """Writes the TEMP card"""
        return self.raw_fields()

    def get_loads(self):
        return [self]

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
