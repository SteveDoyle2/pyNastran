# pylint: disable=R0902,R0904,R0914
"""
All axisymmetric loads are defined in this file.  This includes:

 * FORCEAX
 * PLOADX1

"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.bdf.cards.thermal.loads import ThermalLoad
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class FORCEAX(BaseCard):
    """
    Defines a static concentrated force on a conical shell ring.

    +---------+-----+-----+-----+-----+-----+-----+-----+
    |    1    |  2  |  3  |  4  |  5  |  6  |  7  |  8  |
    +=========+=====+=====+=====+=====+=====+=====+=====+
    | FORCEAX | SID | RID | HID |  S  | FR  | FP  |  FZ |
    +---------+-----+-----+-----+-----+-----+-----+-----+
    | FORCEAX |  1  |  2  |  3  | 2.0 | 0.1 | 0.2 | 0.3 |
    +---------+-----+-----+-----+-----+-----+-----+-----+

    """
    type = 'FORCEAX'
    #_properties = ['node_ids', 'nodes']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        ring_id = 10
        hid = 0
        scale = 1.
        f_rtz = [0., 0., 1.]
        return FORCEAX(sid, ring_id, hid, scale, f_rtz, comment='')

    def __init__(self, sid, ring_id, hid, scale, f_rtz, comment=''):
        """
        Creates a PLOADX1 card, which defines surface traction for
        axisymmetric elements.

        Parameters
        ----------
        sid: int
            Load set identification number.
        ring_id : int
            RINGAX entry identification number.
        hid: int
            Harmonic identification number or a sequence of harmonics.
        scale : float
            Scale factor for the force.
        f_rtz
            Force components in r, φ, z directions.
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.ring_id = ring_id
        self.hid = hid
        self.scale = scale
        self.f_rtz = np.asarray(f_rtz)
        self.ring_id_ref = None

    def validate(self):
        assert isinstance(self.ring_id, int), 'ring_id=%r' % self.ring_id
        assert isinstance(self.hid, int), 'hid=%r' % self.hid

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FORCEAX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        ring_id = integer(card, 2, 'ring_id')
        hid = integer(card, 3, 'hid')
        scale = double(card, 4, 'scale')
        f_rtz = [
            double(card, 5, 'Fr'),
            double_or_blank(card, 6, 'Ft', 0.),
            double_or_blank(card, 7, 'Fz', 0.),
        ]
        assert len(card) <= 8, 'len(FORCEAX card) = %i\ncard=%s' % (len(card), card)
        return FORCEAX(sid, ring_id, hid, scale, f_rtz, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        pass
        #msg = ', which is required by PLOADX1 lid=%s' % self.sid
        #self.ring_id_ref = model.Element(self.eid, msg=msg)

    #@property
    #def node_ids(self):
        #return [self.Ga(), self.Gb()]

    #@property
    #def nodes(self):
        #return [self.ga, self.gb]

    #@property
    #def nodes_ref(self):
        #return [self.ga_ref, self.gb_ref]

    def safe_cross_reference(self, model, safe_coord):
        return self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.ring_id_ref = None


    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = [
            'FORCEAX', self.sid, self.ring_id, self.hid, self.scale] + list(self.f_rtz)
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class PLOADX1(BaseCard):
    """
    Pressure Load on Axisymmetric Element

    Defines surface traction to be used with the CQUADX, CTRIAX, and CTRIAX6
    axisymmetric element.

    +---------+-----+-----+----+----+----+----+-------+
    |    1    |  2  |  3  |  4 |  5 |  6 |  7 |   8   |
    +=========+=====+=====+====+====+====+====+=======+
    | PLOADX1 | SID | EID | PA | PB | GA | GB | THETA |
    +---------+-----+-----+----+----+----+----+-------+

    """
    type = 'PLOADX1'
    _properties = ['node_ids', 'nodes']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        eid = 1
        pa = 1.
        nids = [1, 2]
        return PLOADX1(sid, eid, pa, nids, pb=None, theta=0., comment='')

    def __init__(self, sid, eid, pa, nids, pb=None, theta=0., comment=''):
        """
        Creates a PLOADX1 card, which defines surface traction for
        axisymmetric elements.

        Parameters
        ----------
        sid : int
            load id
        eid : int
            element id (CQUADX, CTRIAX, or CTRIAX6)
        nids : List[int, int]
            Corner grid points.
            GA and GB are any two adjacent corner grid points of the element
        pa / pb : float / None
            Surface traction at grid point GA or GB
            pb : default is None -> pa
        theta : float; default=0.0
            Angle between surface traction and inward normal to the line
            segment.
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        if pb is None:
            pb = pa
        self.sid = sid
        self.eid = eid
        self.pa = pa
        self.pb = pb
        self.ga = nids[0]
        self.gb = nids[1]
        self.theta = theta
        self.eid_ref = None
        self.ga_ref = None
        self.gb_ref = None

    def validate(self):
        assert isinstance(self.ga, integer_types), 'ga=%r' % self.ga
        assert isinstance(self.gb, integer_types), 'gb=%r' % self.gb
        assert isinstance(self.pa, float), 'pa=%r' % self.pa
        assert isinstance(self.pb, float), 'pb=%r' % self.pb

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PLOADX1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2, 'eid')
        pa = double(card, 3, 'pa')
        pb = double_or_blank(card, 4, 'pb', pa)
        ga = integer(card, 5, 'ga')
        gb = integer(card, 6, 'gb')
        theta = double_or_blank(card, 7, 'theta', 0.)
        assert len(card) <= 8, 'len(PLOADX1 card) = %i\ncard=%s' % (len(card), card)
        nids = [ga, gb]
        return PLOADX1(sid, eid, pa, nids, pb=pb, theta=theta, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PLOADX1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        sid, eid, pa, pb, ga, gb, theta = data
        nids = [ga, gb]
        return PLOADX1(sid, eid, pa, nids, pb=pb, theta=theta, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PLOADX1 lid=%s' % self.sid
        self.eid_ref = model.Element(self.eid, msg=msg)
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    @property
    def nodes(self):
        return [self.ga, self.gb]

    @property
    def nodes_ref(self):
        return [self.ga_ref, self.gb_ref]

    def safe_cross_reference(self, model, safe_coord):
        return self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.eid = self.Eid()
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.eid_ref = None
        self.ga_ref = None
        self.gb_ref = None

    def Eid(self):
        if self.eid_ref is not None:
            return self.eid_ref.eid
        return self.eid

    def Ga(self):
        if self.ga_ref is not None:
            return self.ga_ref.nid
        return self.ga

    def Gb(self):
        if self.gb_ref is not None:
            return self.gb_ref.nid
        return self.gb

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = [
            'PLOADX1', self.sid, self.Eid(), self.pa, self.pb,
            self.Ga(), self.Gb(), self.theta]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
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

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        pressure = 1.
        rid1 = 1
        rid2 = 2
        return PRESAX(sid, pressure, rid1, rid2, phi1=0., phi2=360., comment='')

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

    def get_loads(self):
        return [self]

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        pass

    def raw_fields(self):
        list_fields = ['PRESAX', self.sid, self.pressure, self.rid1, self.rid2,
                       self.phi1, self.phi2]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        ring = 1.
        phi = 1.
        temperature = 1.
        return TEMPAX(sid, ring, phi, temperature, comment='')

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

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, debug=True):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
