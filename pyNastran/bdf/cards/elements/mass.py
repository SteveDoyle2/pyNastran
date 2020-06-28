# pylint: disable=C0103,R0902,R0904,R0914,C0111
"""
All mass elements are defined in this file.  This includes:

 * CMASS1
 * CMASS2
 * CMASS3
 * CMASS4
 * CONM1
 * CONM2

All mass elements are PointMassElement and Element objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.cards.base_card import Element
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank)
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def is_positive_semi_definite(A, tol=1e-8):
    """is the 3x3 matrix positive within tolerance"""
    vals = np.linalg.eigh(A)[0]
    return np.all(vals > -tol), vals

class PointMassElement(Element):
    def __init__(self):
        Element.__init__(self)

    def repr_fields(self) -> List[Optional[int, float, str]]:
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

    def write_card_16(self, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_16(card)

# class PointMass(BaseCard):
#     def __init__(self, card, data):
#         self.mass = None
#
#     def Mass(self):
#         return self.mass


class CMASS1(PointMassElement):
    """
    Defines a scalar mass element.

    +--------+-----+-----+----+----+----+----+
    |   1    |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CMASS1 | EID | PID | G1 | C1 | G2 | C2 |
    +--------+-----+-----+----+----+----+----+

    """
    type = 'CMASS1'
    _field_map = {
        1: 'eid', 2:'pid', 3:'g1', 4:'c1', 5:'g2', 6:'c2',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        nids = [1, 2]
        return CMASS1(eid, pid, nids, c1=0, c2=0, comment='')

    def __init__(self, eid, pid, nids, c1=0, c2=0, comment=''):
        # type: (int, int, [int, int], int, int, str) -> CMASS1
        """
        Creates a CMASS1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PMASS)
        nids : List[int, int]
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        PointMassElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.c1 = c1
        self.c2 = c2
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        self.nodes_ref = None
        self.pid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CMASS1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        g1 = integer(card, 3, 'g1')
        c1 = integer_or_blank(card, 4, 'c1')
        g2 = integer_or_blank(card, 5, 'g2')
        c2 = integer_or_blank(card, 6, 'c2')
        assert len(card) <= 7, 'len(CMASS1 card) = %i\ncard=%s' % (len(card), card)
        return CMASS1(eid, pid, [g1, g2], c1, c2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CMASS1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        pid = data[1]
        g1 = data[2]
        g2 = data[3]
        c1 = data[4]
        c2 = data[5]
        return CMASS1(eid, pid, [g1, g2], c1, c2, comment=comment)

    def Mass(self):
        return self.pid_ref.mass

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        mass = self.Mass()
        c1 = self.c1
        c2 = self.c2
        #self.nodes

        assert isinstance(eid, integer_types), 'eid=%r' % eid
        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(mass, float), 'mass=%r' % mass
        assert c1 is None or isinstance(c1, integer_types), 'c1=%r' % c1
        assert c2 is None or isinstance(c2, integer_types), 'c2=%r' % c2

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CMASS1 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.PropertyMass(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CMASS1 eid=%s' % self.eid
        self.nodes_ref, missing_nodes = model.safe_empty_nodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property_mass(self.pid, self.eid, xref_errors, msg=msg)
        if missing_nodes:
            model.log.warning(missing_nodes)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = [self.G1(), self.G2()]
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def G1(self):
        if self.nodes_ref is not None:
            return self.nodes_ref[0].nid
        return self.nodes[0]

    def G2(self):
        if self.nodes[1] in [0, None]:
            return 0
        if self.nodes_ref is not None:
            return self.nodes_ref[1].nid
        return self.nodes[1]

    def Centroid(self):
        """
        Centroid is assumed to be c=(g1+g2)/2.
        If g2 is blank, then the centroid is the location of g1.

        """
        factor = 0.
        p1 = np.array([0., 0., 0.])
        p2 = np.array([0., 0., 0.])
        if self.nodes_ref[0] is not None:
            p1 = self.nodes_ref[0].get_position()
            factor += 1.

        g2 = self.G2()
        if g2 not in [None, 0]:
            p2 = self.nodes_ref[1].get_position()
            factor += 1.
        c = (p1 + p2) / factor
        return c

    def center_of_mass(self):
        return self.Centroid()

    @property
    def node_ids(self):
        g1 = self.G1()
        g2 = self.G2()
        nodes = []
        if g1:
            nodes.append(g1)
        if g2:
            nodes.append(g2)
        return nodes

    def Pid(self):
        if self.pid_ref is None:
            return self.pid
        return self.pid_ref.pid

    def raw_fields(self):
        fields = ['CMASS1', self.eid, self.Pid(), self.G1(), self.c1,
                  self.G2(), self.c2]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class CMASS2(PointMassElement):
    """
    Defines a scalar mass element without reference to a property entry.

    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CMASS2 | EID |  M  | G1 | C1 | G2 | C2 |
    +--------+-----+-----+----+----+----+----+

    """
    type = 'CMASS2'
    _field_map = {
        1: 'eid', 2:'mass', 3:'g1', 4:'c1', 5:'g2', 6:'c2',
    }
    cp_name_map = {
        'M' : 'mass',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        mass = 1.
        nids = [1, 2]
        return CMASS2(eid, mass, nids, c1=0, c2=0, comment='')

    def __init__(self, eid, mass, nids, c1, c2, comment=''):
        """
        Creates a CMASS2 card

        Parameters
        ----------
        eid : int
            element id
        mass : float
            mass
        nids : List[int, int]
            node ids
        c1 / c2 : int; default=None
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        PointMassElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.mass = mass
        self.nodes = nids
        self.c1 = c1
        self.c2 = c2
        self.nodes_ref = None
        assert len(self.nodes) == 2, self.nodes

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the masses in a vectorized way"""
        #comments = []
        mass = []
        nodes = []
        components = []
        for eid in eids:
            element = model.masses[eid]
            #comments.append(element.comment)
            mass.append(element.mass)
            nodes.append([nid if nid is not None else 0 for nid in element.nodes])
            components.append([comp if comp is not None else 0
                               for comp in [element.c1, element.c2]])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('mass', data=mass)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('components', data=components)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CMASS2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        mass = double_or_blank(card, 2, 'mass', 0.)
        g1 = integer_or_blank(card, 3, 'g1')
        c1 = integer_or_blank(card, 4, 'c1')
        g2 = integer_or_blank(card, 5, 'g2')
        c2 = integer_or_blank(card, 6, 'c2')
        assert len(card) <= 7, 'len(CMASS2 card) = %i\ncard=%s' % (len(card), card)
        return CMASS2(eid, mass, [g1, g2], c1, c2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CMASS2 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        mass = data[1]
        g1 = data[2]
        g2 = data[3]
        c1 = data[4]
        c2 = data[5]
        #assert g1 > 0, data
        if g1 == 0:
            g1 = None
        else:
            assert g1 > 0, f'g1={g1}; g2={g2} c1={c1} c2={c2}'
        if g2 == 0:
            g2 = None
        else:
            assert g2 > 0, f'g2={g2}; g1={g1} c1={c1} c2={c2}'

        assert 0 <= c1 <= 123456, 'c1=%s data=%s' % (c1, data)
        assert 0 <= c2 <= 123456, 'c2=%s data=%s' % (c2, data)
        return CMASS2(eid, mass, [g1, g2], c1, c2, comment=comment)

    def validate(self):
        assert len(self.nodes) == 2, self.nodes

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        mass = self.Mass()
        c1 = self.c1
        c2 = self.c2
        #self.nodes

        assert isinstance(eid, integer_types), 'eid=%r' % eid
        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(mass, float), 'mass=%r' % mass
        assert c1 is None or isinstance(c1, integer_types), 'c1=%r' % c1
        assert c2 is None or isinstance(c2, integer_types), 'c2=%r' % c2

    #@property
    #def nodes(self):
        #return [self.g1, self.g2]

    @property
    def node_ids(self):
        g1 = self.G1()
        g2 = self.G2()
        nodes = []
        if g1:
            nodes.append(g1)
        if g2:
            nodes.append(g2)
        return nodes

    def Mass(self):
        return self.mass

    def Centroid(self):
        """
        Centroid is assumed to be c=(g1+g2)/2.
        If g2 is blank, then the centroid is the location of g1.
        """
        factor = 0.
        p1 = np.array([0., 0., 0.])
        p2 = np.array([0., 0., 0.])
        if self.nodes[0] is not None:
            p1 = self.nodes_ref[0].get_position()
            factor += 1.
        if self.nodes[1] is not None:
            p2 = self.nodes_ref[1].get_position()
            factor += 1.
        assert factor > 0., str(self)
        c = (p1 + p2) / factor
        return c

    def center_of_mass(self):
        return self.Centroid()

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CMASS2 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CMASS2 eid=%s' % self.eid
        self.nodes_ref, missing_nodes = model.safe_empty_nodes(self.node_ids, msg=msg)
        if missing_nodes:
            model.log.warning(missing_nodes)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = [self.G1(), self.G2()]
        self.nodes_ref = None

    def G1(self):
        if self.nodes_ref is not None and self.nodes_ref[0] is not None:
            return self.nodes_ref[0].nid
        return self.nodes[0]

    def G2(self):
        if self.nodes_ref is not None and self.nodes_ref[1] is not None:
            return self.nodes_ref[1].nid
        return self.nodes[1]

    def raw_fields(self):
        fields = ['CMASS2', self.eid, self.mass, self.G1(),
                  self.c1, self.G2(), self.c2]
        return fields

    def repr_fields(self):
        mass = set_blank_if_default(self.mass, 0.)
        fields = ['CMASS2', self.eid, mass, self.G1(), self.c1,
                  self.G2(), self.c2]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class CMASS3(PointMassElement):
    """
    Defines a scalar mass element that is connected only to scalar points.

    +--------+-----+-----+----+----+
    |   1    |  2  |   3 |  4 |  5 |
    +========+=====+=====+====+====+
    | CMASS3 | EID | PID | S1 | S2 |
    +--------+-----+-----+----+----+

    """
    type = 'CMASS3'
    _field_map = {
        1: 'eid', 2:'pid', 3:'s1', 4:'s2',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        nids = [1, 2]
        return CMASS3(eid, pid, nids, comment='')

    def __init__(self, eid, pid, nids, comment=''):
        """
        Creates a CMASS3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PMASS)
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        PointMassElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.nodes = nids
        assert self.nodes[0] != self.nodes[1]
        self.nodes_ref = None
        self.pid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CMASS3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        s1 = integer_or_blank(card, 3, 's1')
        s2 = integer_or_blank(card, 4, 's2')
        assert len(card) <= 5, 'len(CMASS3 card) = %i\ncard=%s' % (len(card), card)
        return CMASS3(eid, pid, [s1, s2], comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CMASS3 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        pid = data[1]
        s1 = data[2]
        s2 = data[3]
        return CMASS3(eid, pid, [s1, s2], comment=comment)

    def Mass(self):
        return self.pid_ref.mass

    def S1(self):
        if self.nodes_ref is not None:
            return self.nodes_ref[0].nid
        return self.nodes[0]

    def S2(self):
        if self.nodes_ref is not None:
            return self.nodes_ref[1].nid
        return self.nodes[1]

    @property
    def node_ids(self):
        return [self.S1(), self.S2()]

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CMASS3 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.PropertyMass(self.pid, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    #def Mass(self):
        #return self.mass

    def Centroid(self):
        return np.zeros(3)

    def center_of_mass(self):
        return self.Centroid()

    def raw_fields(self):
        fields = ['CMASS3', self.eid, self.Pid(), self.S1(), self.S2()]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class CMASS4(PointMassElement):
    """
    Defines a scalar mass element that is connected only to scalar points,
    without reference to a property entry

    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CMASS4 | EID |  M  | S1 | S2 |
    +--------+-----+-----+----+----+

    """
    type = 'CMASS4'
    #cp_name_map/update_by_cp_name
    cp_name_map = {
        'M' : 'mass',
    }
    _properties = ['node_ids']
    _field_map = {
        1: 'eid', 2:'mass', 3:'s1', 4:'s2',
    }

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        nids = [1, 2]
        return CMASS4(eid, pid, nids, comment='')

    def __init__(self, eid, mass, nids, comment=''):
        """
        Creates a CMASS4 card

        Parameters
        ----------
        eid : int
            element id
        mass : float
            SPOINT mass
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        PointMassElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.mass = mass
        self.nodes = nids
        assert self.nodes[0] != self.nodes[1]
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        ioffset = icard * 4
        eid = integer(card, 1 + ioffset, 'eid')
        mass = double_or_blank(card, 2 + ioffset, 'mass', 0.)
        s1 = integer(card, 3 + ioffset, 's1')
        s2 = integer_or_blank(card, 4 + ioffset, 's2', 0)
        assert len(card) <= 9, 'len(CMASS4 card) = %i\ncard=%s' % (len(card), card)
        return CMASS4(eid, mass, [s1, s2], comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CMASS4 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        mass = data[1]
        s1 = data[2]
        s2 = data[3]
        return CMASS4(eid, mass, [s1, s2], comment=comment)

    def Mass(self):
        return self.mass

    def Centroid(self):
        return np.zeros(3)

    def center_of_mass(self):
        return self.Centroid()

    @property
    def node_ids(self):
        return [self.S1(), self.S2()]

    def S1(self):
        if self.nodes_ref is not None and self.nodes_ref[0] is not None:
            return self.nodes_ref[0].nid
        return self.nodes[0]

    def S2(self):
        if self.nodes_ref is not None and self.nodes_ref[1] is not None:
            return self.nodes_ref[1].nid
        return self.nodes[1]

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CMASS4 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = [self.S1(), self.S2()]
        self.nodes_ref = None

    def raw_fields(self):
        fields = ['CMASS4', self.eid, self.mass] + self.node_ids
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class CONM1(PointMassElement):
    """
    Concentrated Mass Element Connection, General Form
    Defines a 6 x 6 symmetric mass matrix at a geometric grid point

    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    |    1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+=====+=====+=====+=====+=====+=====+=====+=====+
    |  CONM1 | EID |  G  | CID | M11 | M21 | M22 | M31 | M32 |
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    |        | M33 | M41 | M42 | M43 | M44 | M51 | M52 | M53 |
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    |        | M54 | M55 | M61 | M62 | M63 | M64 | M65 | M66 |
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+

    """
    type = 'CONM1'
    _field_map = {
        1: 'eid', 2:'nid', 3:'cid',
    }
    _properties = ['node_ids']
    def _update_field_helper(self, n, value):
        m = self.mass_matrix
        if n == 4:
            m[0, 0] = value
        elif n == 5:
            m[1, 0] = value
        elif n == 6:
            m[1, 1] = value
        elif n == 7:
            m[2, 0] = value
        elif n == 8:
            m[2, 1] = value
        elif n == 9:
            m[2, 2] = value
        elif n == 10:
            m[3, 0] = value
        elif n == 11:
            m[3, 1] = value
        elif n == 12:
            m[3, 2] = value
        elif n == 13:
            m[3, 3] = value
        elif n == 14:
            m[4, 0] = value
        elif n == 15:
            m[4, 1] = value
        elif n == 16:
            m[4, 2] = value
        elif n == 17:
            m[4, 3] = value
        elif n == 18:
            m[4, 4] = value
        elif n == 19:
            m[5, 0] = value
        elif n == 20:
            m[5, 1] = value
        elif n == 21:
            m[5, 2] = value
        elif n == 22:
            m[5, 3] = value
        elif n == 23:
            m[5, 4] = value
        elif n == 24:
            m[5, 5] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def update_by_cp_name(self, name, value):
        if name == 'M11':
            self.mass_matrix[0, 0] = value
        #elif name == 'M21':
            #self.mass_matrix[1, 0] = value
        elif name == 'M22':
            self.mass_matrix[1, 1] = value
        elif name == 'M33':
            self.mass_matrix[2, 2] = value
        elif name == 'M44':
            self.mass_matrix[3, 3] = value
        #elif name == 'X1':
            #self.X[0] = value
        else:
            raise NotImplementedError('element_type=%r has not implemented %r in cp_name_map' % (
                self.type, name))

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nid = 1
        mass_matrix = np.zeros((6, 6))
        return CONM1(eid, nid, mass_matrix, cid=0, comment='')

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        self.mass_matrix = np.asarray(self.mass_matrix)

    def __init__(self, eid, nid, mass_matrix, cid=0, comment=''):
        """
        Creates a CONM1 card

        Parameters
        ----------
        eid : int
            element id
        nid : int
            the node to put the mass matrix
        mass_matrix : (6, 6) float ndarray
            the 6x6 mass matrix, M
        cid : int; default=0
            the coordinate system for the mass matrix
        comment : str; default=''
            a comment for the card

        ::

          [M] = [M11 M21 M31 M41 M51 M61]
                [    M22 M32 M42 M52 M62]
                [        M33 M43 M53 M63]
                [            M44 M54 M64]
                [    Sym         M55 M65]
                [                    M66]

        """
        PointMassElement.__init__(self)
        if comment:
            self.comment = comment
        self.mass_matrix = mass_matrix
        self.eid = eid
        self.nid = nid
        self.cid = cid
        self.nid_ref = None
        self.cid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CONM2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        m = np.zeros((6, 6))
        eid = integer(card, 1, 'eid')
        nid = integer(card, 2, 'nid')
        cid = integer_or_blank(card, 3, 'cid', 0)

        m[0, 0] = double_or_blank(card, 4, 'M11', 0.)
        m[1, 0] = double_or_blank(card, 5, 'M21', 0.)
        m[1, 1] = double_or_blank(card, 6, 'M22', 0.)
        m[2, 0] = double_or_blank(card, 7, 'M31', 0.)
        m[2, 1] = double_or_blank(card, 8, 'M32', 0.)
        m[2, 2] = double_or_blank(card, 9, 'M33', 0.)
        m[3, 0] = double_or_blank(card, 10, 'M41', 0.)
        m[3, 1] = double_or_blank(card, 11, 'M42', 0.)
        m[3, 2] = double_or_blank(card, 12, 'M43', 0.)
        m[3, 3] = double_or_blank(card, 13, 'M44', 0.)
        m[4, 0] = double_or_blank(card, 14, 'M51', 0.)
        m[4, 1] = double_or_blank(card, 15, 'M52', 0.)
        m[4, 2] = double_or_blank(card, 16, 'M53', 0.)
        m[4, 3] = double_or_blank(card, 17, 'M54', 0.)
        m[4, 4] = double_or_blank(card, 18, 'M55', 0.)
        m[5, 0] = double_or_blank(card, 19, 'M61', 0.)
        m[5, 1] = double_or_blank(card, 20, 'M62', 0.)
        m[5, 2] = double_or_blank(card, 21, 'M63', 0.)
        m[5, 3] = double_or_blank(card, 22, 'M64', 0.)
        m[5, 4] = double_or_blank(card, 23, 'M65', 0.)
        m[5, 5] = double_or_blank(card, 24, 'M66', 0.)
        assert len(card) <= 25, 'len(CONM1 card) = %i\ncard=%s' % (len(card), card)
        return CONM1(eid, nid, m, cid=cid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CONM1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        m = np.zeros((6, 6))
        (eid, nid, cid, m1, m2a, m2b, m3a, m3b, m3c, m4a, m4b, m4c, m4d,
         m5a, m5b, m5c, m5d, m5e, m6a, m6b, m6c, m6d, m6e, m6f) = data
        m[0, 0] = m1   # M11
        m[1, 0] = m2a  # M21
        m[1, 1] = m2b  # M22
        m[2, 0] = m3a  # M31
        m[2, 1] = m3b  # M32
        m[2, 2] = m3c  # M33
        m[3, 0] = m4a  # M41
        m[3, 1] = m4b  # M42
        m[3, 2] = m4c  # M43
        m[3, 3] = m4d  # M44
        m[4, 0] = m5a  # M51
        m[4, 1] = m5b  # M52
        m[4, 2] = m5c  # M53
        m[4, 3] = m5d  # M54
        m[4, 4] = m5e  # M55
        m[5, 0] = m6a  # M61
        m[5, 1] = m6b  # M62
        m[5, 2] = m6c  # M63
        m[5, 3] = m6d  # M64
        m[5, 4] = m6e  # M65
        m[5, 5] = m6f  # M66
        return CONM1(eid, nid, m, cid=cid, comment=comment)

    def _verify(self, xref):
        eid = self.eid
        assert isinstance(eid, integer_types), 'eid=%r' % eid

    @staticmethod
    def Mass():
        return 0.0

    @staticmethod
    def Centroid():
        return np.zeros(3, dtype='float64')

    def center_of_mass(self):
        return self.Centroid()

    @property
    def node_ids(self):
        return [self.Nid()]

    def Nid(self):
        if self.nid_ref is not None:
            return self.nid_ref.nid
        return self.nid

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CONM1 eid=%s' % self.eid
        self.nid_ref = model.Node(self.nid, msg=msg)
        self.cid_ref = model.Coord(self.cid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CONM1 eid=%s' % self.eid
        self.nid_ref = model.Node(self.nid, msg=msg)
        self.cid_ref = model.safe_coord(self.cid, self.eid, xref_errors, msg='')

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nid = self.Nid()
        self.cid = self.Cid()
        self.nid_ref = None
        self.cid_ref = None

    def MassMatrix(self):
        return self.mass_matrix

    def raw_fields(self):
        cid = set_blank_if_default(self.Cid(), 0)
        nid = self.Nid()
        m = self.mass_matrix
        list_fields = [
            'CONM1', self.eid, nid, cid, m[0, 0], m[1, 0], m[1, 1],
            m[2, 0], m[2, 1], m[2, 2], m[3, 0], m[3, 1], m[3, 2],
            m[3, 3], m[4, 0], m[4, 1], m[4, 2], m[4, 3], m[4, 4],
            m[5, 0], m[5, 1], m[5, 2], m[5, 3], m[5, 4], m[5, 5]]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        list_fields2 = list_fields[0:4]
        for field in list_fields[4:]:
            val = set_blank_if_default(field, 0.)
            list_fields2.append(val)
        return list_fields2

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class CONM2(PointMassElement):
    """
    +-------+--------+-------+-------+---------+------+------+------+
    |   1   |    2   |    3  |   4   |    5    |  6   |  7   |   8  |
    +=======+========+=======+=======+=========+======+======+======+
    | CONM2 |   EID  |  NID  |  CID  |  MASS   |  X1  |  X2  |  X3  |
    +-------+--------+-------+-------+---------+------+------+------+
    |       |   I11  |  I21  |  I22  |   I31   |  I32 |  I33 |      |
    +-------+--------+-------+-------+---------+------+------+------+
    | CONM2 | 501274 | 11064 |       | 132.274 |      |      |      |
    +-------+--------+-------+-------+---------+------+------+------+

    """
    type = 'CONM2'
    _field_map = {
        1: 'eid', 2:'nid', 3:'cid', 4:'mass',
    }
    def update_by_cp_name(self, name, value):
        if name == 'M':
            self.mass = value
        elif name == 'X1':
            self.X[0] = value
        elif name == 'X2':
            self.X[1] = value
        elif name == 'X3':
            self.X[2] = value
        elif name == 'I11':
            #I11, I21, I22, I31, I32, I33 = I
            self.I[0] = value
        elif name == 'I21':
            self.I[1] = value
        elif name == 'I22':
            self.I[2] = value
        elif name == 'I31':
            self.I[3] = value
        elif name == 'I32':
            self.I[4] = value
        elif name == 'I33':
            self.I[5] = value
        else:
            raise NotImplementedError('element_type=%r has not implemented %r in cp_name_map' % (
                self.type, name))

    def _update_field_helper(self, n, value):
        if n == 5:
            self.X[0] = value
        elif n == 6:
            self.X[1] = value
        elif n == 7:
            self.X[2] = value

        elif n == 9:
            self.I[0] = value
        elif n == 10:
            self.I[1] = value
        elif n == 11:
            self.I[2] = value
        elif n == 12:
            self.I[3] = value
        elif n == 13:
            self.I[4] = value
        elif n == 14:
            self.I[5] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        unused_comments = []
        nid = []
        cid = []
        mass = []
        X = []
        I = []
        for eid in eids:
            element = model.masses[eid]
            #comments.append(element.comment)
            nid.append(element.nid)
            cid.append(element.cid)
            mass.append(element.mass)
            X.append(element.X)
            I.append(element.I)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('nid', data=nid)
        h5_file.create_dataset('cid', data=cid)
        h5_file.create_dataset('X', data=X)
        h5_file.create_dataset('I', data=I)
        h5_file.create_dataset('mass', data=mass)

        #self.eid = eid
        #self.nid = nid
        #self.cid = cid
        #self.mass = mass
        #self.X = np.asarray(X)
        #self.I = np.asarray(I)

    def __init__(self, eid, nid, mass, cid=0, X=None, I=None, comment=''):
        """
        Creates a CONM2 card

        Parameters
        ----------
        eid : int
           element id
        nid : int
           node id
        mass : float
           the mass of the CONM2
        cid : int; default=0
           coordinate frame of the offset (-1=absolute coordinates)
        X : (3, ) List[float]; default=None -> [0., 0., 0.]
            xyz offset vector relative to nid
        I : (6, ) List[float]; default=None -> [0., 0., 0., 0., 0., 0.]
            mass moment of inertia matrix about the CG
            I11, I21, I22, I31, I32, I33 = I
        comment : str; default=''
            a comment for the card

        """
        PointMassElement.__init__(self)
        if comment:
            self.comment = comment
        #: Element identification number. (0 < Integer < 100,000,000)
        self.eid = eid

        #: Grid point identification number. (Integer > 0)
        self.nid = nid

        #: Coordinate system identification number.
        #: For CID of -1; see X1, X2, X3 below.
        #: (Integer > -1; Default = 0)
        self.cid = cid

        #: Mass value. (Real)
        self.mass = mass

        if X is None:
            X = np.zeros(3)
        #: Offset distances from the grid point to the center of gravity of
        #: the mass in the coordinate system defined in field 4, unless
        #: CID = -1, in which case X1, X2, X3 are the coordinates, not
        #: offsets, of the center of gravity of the mass in the basic
        #: coordinate system. (Real)
        self.X = np.asarray(X)

        if I is None:
            I = np.zeros(6)
        #: Mass moments of inertia measured at the mass center of gravity in
        #: the coordinate system defined by field 4. If CID = -1, the basic
        #: coordinate system is implied. (Real)
        #: I11, I21, I22, I31, I32, I33 = I
        self.I = np.asarray(I)
        self.nid_ref = None
        self.cid_ref = None

    def validate(self):
        assert isinstance(self.cid, integer_types), self.cid
        assert isinstance(self.mass, float), self.mass
        assert self.mass >= 0., 'mass=%s' % self.mass

        assert isinstance(self.X, np.ndarray), self.X
        assert isinstance(self.I, np.ndarray), self.I
        assert len(self.X) == 3, self.X
        assert len(self.I) == 6, self.I

        I11, I12, I22, I13, I23, I33 = self.I
        I = np.array([
            [I11, I12, I13],
            [I12, I22, I23],
            [I13, I23, I33],
        ], dtype='float32')
        is_psd, eigi = is_positive_semi_definite(I)
        if not is_psd:
            msg = 'The eig(I) >= 0.\n'
            msg += 'I=\n%s\n' % str(I)
            msg += 'eigenvalues=%s' % str(eigi)
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CONM2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        nid = integer(card, 2, 'nid')
        cid = integer_or_blank(card, 3, 'cid', 0)
        mass = double_or_blank(card, 4, 'mass', 0.)

        X = [
            double_or_blank(card, 5, 'x1', 0.0),
            double_or_blank(card, 6, 'x2', 0.0),
            double_or_blank(card, 7, 'x3', 0.0),
        ]

        I = [
            double_or_blank(card, 9, 'I11', 0.0),
            double_or_blank(card, 10, 'I21', 0.0),
            double_or_blank(card, 11, 'I22', 0.0),
            double_or_blank(card, 12, 'I31', 0.0),
            double_or_blank(card, 13, 'I32', 0.0),
            double_or_blank(card, 14, 'I33', 0.0),
        ]
        assert len(card) <= 15, 'len(CONM2 card) = %i\ncard=%s' % (len(card), card)
        return CONM2(eid, nid, mass, cid=cid, X=X, I=I, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CONM2 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        nid = data[1]
        cid = data[2]
        mass = data[3]
        X = data[4:7]
        I = data[7:]
        return CONM2(eid, nid, mass, cid=cid, X=X, I=I, comment=comment)

    def _verify(self, xref):
        eid = self.eid
        nid = self.Nid()
        cid = self.Cid()
        mass = self.Mass()
        c = self.Centroid()

        assert isinstance(eid, integer_types), 'eid=%r' % eid
        assert isinstance(nid, integer_types), 'nid=%r' % nid
        assert isinstance(cid, integer_types), 'cid=%r' % cid
        assert isinstance(mass, float), 'mass=%r' % mass
        for i in range(3):
            assert isinstance(c[i], float), 'centroid[%i]=%r' % (i, c[i])

    def Mass(self):
        return self.mass

    def Inertia(self):
        """
        Returns the 3x3 inertia matrix
        .. warning:: doesnt handle offsets or coordinate systems

        """
        I = self.I
        A = [[I[0], -I[1], -I[3]],
             [-I[1], I[2], -I[4]],
             [-I[3], -I[4], I[5]]]
        if self.Cid() in [0, -1]:
            return A
        else:
            # transform to global
            #dx = self.cid_ref.transform_node_to_global(self.X)
            #matrix = self.cid_ref.beta()
            raise NotImplementedError('CONM2 intertia method for CID != 0 is not implemented.')
            #A2 = A * matrix
            #return A2  # correct for offset using dx???

    def offset(self, xyz_nid):
        cid = self.Cid()
        if cid == 0:
            # no transform needed
            X2 = xyz_nid + self.X
        elif cid == -1:
            # case X1, X2, X3 are the coordinates, not offsets, of the center of gravity of
            # the mass in the basic coordinate system.

            # 4. If CID = -1, offsets are internally computed as the difference between the grid
            # point location and X1, X2, X3.
            # this statement is not supported...
            return self.X
        else:
            # Offset distances from the grid point to the center of gravity of the mass
            # in the coordinate system

            # If CID > 0, then X1, X2, and X3 are defined by a local Cartesian system, even
            # if CID references a spherical or cylindrical coordinate system. This is similar
            # to the manner in which displacement coordinate systems are defined.
            # this statement is not supported...

            # convert self.X into the global frame
            x = self.cid_ref.transform_node_to_global(self.X)

            # self.X is an offset
            dx = x - self.cid_ref.origin

            # the actual position of the CONM2
            X2 = xyz_nid + dx
        return X2

    def Centroid(self):
        """
        This method seems way more complicated than it needs to be thanks
        to all these little caveats that don't seem to be supported.
        """
        cid = self.Cid()
        if cid == 0:
            # no transform needed
            X2 = self.nid_ref.get_position() + self.X
        elif cid == -1:
            # case X1, X2, X3 are the coordinates, not offsets, of the center of gravity of
            # the mass in the basic coordinate system.

            # 4. If CID = -1, offsets are internally computed as the difference between the grid
            # point location and X1, X2, X3.
            # this statement is not supported...
            return self.X
        else:
            # Offset distances from the grid point to the center of gravity of the mass
            # in the coordinate system

            # If CID > 0, then X1, X2, and X3 are defined by a local Cartesian system, even
            # if CID references a spherical or cylindrical coordinate system. This is similar
            # to the manner in which displacement coordinate systems are defined.
            # this statement is not supported...

            # convert self.X into the global frame
            x = self.cid_ref.transform_node_to_global(self.X)

            # self.X is an offset
            dx = x - self.cid_ref.origin

            # the actual position of the CONM2
            X2 = self.nid_ref.get_position() + dx
        return X2

    def Centroid_no_xref(self, model):
        """
        This method seems way more complicated than it needs to be thanks
        to all these little caveats that don't seem to be supported.

        """
        cid = self.cid
        nid_ref = model.Node(self.nid)
        if cid == 0:
            # no transform needed
            X2 = nid_ref.get_position() + self.X
        elif cid == -1:
            # case X1, X2, X3 are the coordinates, not offsets, of the center of gravity of
            # the mass in the basic coordinate system.

            # 4. If CID = -1, offsets are internally computed as the difference between the grid
            # point location and X1, X2, X3.
            # this statement is not supported...
            return self.X
        else:
            cid_ref = model.Coord(self.cid)
            # Offset distances from the grid point to the center of gravity of the mass
            # in the coordinate system

            # If CID > 0, then X1, X2, and X3 are defined by a local Cartesian system, even
            # if CID references a spherical or cylindrical coordinate system. This is similar
            # to the manner in which displacement coordinate systems are defined.
            # this statement is not supported...

            # convert self.X into the global frame
            x = cid_ref.transform_node_to_global(self.X)

            # self.X is an offset
            dx = x - cid_ref.origin

            # the actual position of the CONM2
            X2 = nid_ref.get_position() + dx
        return X2

    def center_of_mass(self):
        return self.Centroid()

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CONM2 eid=%s' % self.eid
        self.nid_ref = model.Node(self.nid, msg=msg)

        cid = self.Cid()
        if cid != -1:
            self.cid_ref = model.Coord(cid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CONM2 eid=%s' % self.eid
        self.nid_ref = model.Node(self.nid, msg=msg)
        cid = self.Cid()
        if cid != -1:
            self.cid_ref = model.safe_coord(cid, self.eid, xref_errors, msg='')

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nid = self.Nid()
        self.cid = self.Cid()
        self.nid_ref = None
        self.cid_ref = None

    @property
    def node_ids(self):
        return [self.Nid()]

    def Nid(self):
        if self.nid_ref is not None:
            return self.nid_ref.nid
        return self.nid

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def raw_fields(self):
        list_fields = (['CONM2', self.eid, self.Nid(), self.Cid(), self.mass] +
                       list(self.X) + [None] + list(self.I))
        return list_fields

    def repr_fields(self):
        I = []
        for i in self.I:
            if i == 0.:
                I.append(None)
            else:
                I.append(i)
        X = []
        for x in self.X:
            if x == 0.:
                X.append(None)
            else:
                X.append(x)

        cid = set_blank_if_default(self.Cid(), 0)
        list_fields = (['CONM2', self.eid, self.Nid(), cid, self.mass] + X +
                       [None] + I)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

    def write_card_16(self, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_16(card)
