# pylint: disable=C0103,R0902,R0904,R0914
"""
All spring elements are defined in this file.  This includes:

 * CELAS1
 * CELAS2
 * CELAS3
 * CELAS4

All spring elements are SpringElement and Element objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Element
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class SpringElement(Element):
    def __init__(self):
        Element.__init__(self)
        self.nodes = [None, None]

    def Mass(self):
        return 0.0

    def repr_fields(self):
        return self.raw_fields()

    def write_card_16(self, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_16(card)

class CELAS1(SpringElement):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CELAS1 | EID | PID | G1 | C1 | G2 | C2 |
    +--------+-----+-----+----+----+----+----+
    """
    type = 'CELAS1'
    _field_map = {
        1: 'eid', 2:'pid', 4:'c1', 6:'c2',
    }
    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 5:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, pid, nids, c1=0, c2=0, comment=''):
        """
        Creates a CELAS1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : List[int, int]
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card
        """
        SpringElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        #: property ID
        self.pid = pid
        #: component number
        self.c1 = c1
        self.c2 = c2
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        self.nodes_ref = None
        self.pid_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        components = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append([nid if nid is not None else 0 for nid in element.nodes])
            components.append([element.c1, element.c2])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('components', data=components)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CELAS1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'g1'), integer_or_blank(card, 5, 'g2', 0)]

        #: component number
        c1 = integer_or_blank(card, 4, 'c1', 0)
        c2 = integer_or_blank(card, 6, 'c2', 0)
        assert len(card) <= 7, 'len(CELAS1 card) = %i\ncard=%s' % (len(card), card)
        return CELAS1(eid, pid, nids, c1, c2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CELAS1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        pid = data[1]
        nids = [data[2], data[3]]
        c1 = data[4]
        c2 = data[5]
        return CELAS1(eid, pid, nids, c1, c2, comment=comment)

    def validate(self):
        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=%r %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=%r %s' % (self.c2, msg)
        assert len(self.nodes) == 2

    @property
    def node_ids(self):
        msg = ', which is required by CELAS1 eid=%s' % (self.eid)
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def _verify(self, xref):
        eid = self.eid
        node_ids = self.node_ids
        c1 = self.c1
        c2 = self.c2
        #ge = self.ge
        #s = self.s

        assert isinstance(eid, int), 'eid=%r' % eid
        assert isinstance(c1, int), 'c1=%r' % c1
        assert isinstance(c2, int), 'c2=%r' % c2
        #assert isinstance(ge, float), 'ge=%r' % ge
        #assert isinstance(s, float), 'ge=%r' % s
        if xref:
            k = self.K()
            assert self.pid_ref.type in ['PELAS'], self.pid_ref
            assert isinstance(k, float), 'k=%r' % k
            assert len(node_ids) == len(self.nodes)
            #for nodeID, node in zip(node_ids, self.nodes):
                #assert node.node.nid

    def K(self):
        return self.pid_ref.k

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CELAS1 eid=%s' % (self.eid)
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.Pid(), msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CELAS1 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def raw_fields(self):
        nodes = self.node_ids
        list_fields = ['CELAS1', self.eid, self.Pid(), nodes[0],
                       self.c1, nodes[1], self.c2]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CELAS2(SpringElement):
    """
    +--------+-----+-----+----+----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |  8 |  9 |
    +========+=====+=====+====+====+====+====+====+====+
    | CELAS2 | EID |  K  | G1 | C1 | G2 | C2 | GE | S  |
    +--------+-----+-----+----+----+----+----+----+----+
    """
    type = 'CELAS2'
    _field_map = {
        1: 'eid', 2:'k', 4:'c1', 6:'c2',
    }
    cp_name_map = {'K' : 'k', 'GE' : 'ge', 'S' : 's'}
    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 5:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, k, nids, c1=0, c2=0, ge=0., s=0., comment=''):
        """
        Creates a CELAS2 card

        Parameters
        ----------
        eid : int
            element id
        k : float
            spring stiffness
        nids : List[int, int]
            SPOINT ids
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        ge : int; default=0.0
            damping coefficient
        s : float; default=0.0
            stress coefficient
        comment : str; default=''
            a comment for the card
        """
        SpringElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        #: stiffness of the scalar spring
        self.k = k
        #: component number
        self.c1 = c1
        self.c2 = c2
        #: damping coefficient
        self.ge = ge
        #: stress coefficient
        self.s = s
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        self.nodes_ref = None
        self.pid_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        k = []
        ge = []
        s = []
        nodes = []
        components = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            k.append(element.k)
            ge.append(element.ge)
            s.append(element.s)
            nodes.append([nid if nid is not None else 0 for nid in element.nodes])
            components.append([element.c1, element.c2])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('K', data=k)
        h5_file.create_dataset('ge', data=ge)
        h5_file.create_dataset('s', data=s)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('components', data=components)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CELAS2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        k = double(card, 2, 'k')
        nids = [integer_or_blank(card, 3, 'g1', 0),
                integer_or_blank(card, 5, 'g2', 0)]
        c1 = integer_or_blank(card, 4, 'c1', 0)
        c2 = integer_or_blank(card, 6, 'c2', 0)
        ge = double_or_blank(card, 7, 'ge', 0.)
        s = double_or_blank(card, 8, 's', 0.)
        assert len(card) <= 9, 'len(CELAS2 card) = %i\ncard=%s' % (len(card), card)
        return CELAS2(eid, k, nids, c1, c2, ge, s, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CELAS2 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        k = data[1]
        nids = [data[2], data[3]]
        c1 = data[4]
        c2 = data[5]
        ge = data[6]
        s = data[7]
        return CELAS2(eid, k, nids, c1, c2, ge, s, comment=comment)

    def validate(self):
        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=%r %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=%r %s' % (self.c2, msg)
        if self.nodes[0] == self.nodes[1] and self.c1 == self.c2:
            msg = (
                'Thee nodes=%s must be unique or dofs=[%s, %s] must not be '
                'the same; CELAS2 eid=%s' % (self.nodes, self.c1, self.c2, self.eid))
            raise AssertionError(msg)

    @property
    def node_ids(self):
        msg = ', which is required by CELAS2 eid=%s' % (self.eid)
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CELAS2 eid=%s' % (self.eid)
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CELAS2 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def _verify(self, xref=True):
        eid = self.eid
        k = self.K()
        node_ids = self.node_ids
        c1 = self.c2
        c2 = self.c1
        ge = self.ge
        s = self.s

        assert isinstance(eid, int), 'eid=%r' % eid
        assert isinstance(c1, int), 'c1=%r' % c1
        assert isinstance(c2, int), 'c2=%r' % c2
        assert isinstance(k, float), 'k=%r' % k
        assert isinstance(ge, float), 'ge=%r' % ge
        assert isinstance(s, float), 'ge=%r' % s
        if xref:
            assert len(node_ids) == len(self.nodes)
            #for node_id, node in zip(node_id, self.nodes):
                #assert node.node.nid

    def K(self):
        return self.k

    def raw_fields(self):
        nodes = self.node_ids
        list_fields = ['CELAS2', self.eid, self.k, nodes[0], self.c1,
                       nodes[1], self.c2, self.ge, self.s]
        return list_fields

    def repr_fields(self):
        nodes = self.node_ids
        ge = set_blank_if_default(self.ge, 0.)
        s = set_blank_if_default(self.s, 0.)
        list_fields = ['CELAS2', self.eid, self.k, nodes[0], self.c1,
                       nodes[1], self.c2, ge, s]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CELAS3(SpringElement):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CELAS3 | EID | PID | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    type = 'CELAS3'
    _field_map = {
        1: 'eid', 2:'pid', #4:'s1', 6:'s2',
    }

    def __init__(self, eid, pid, nodes, comment=''):
        """
        Creates a CELAS3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card
        """
        SpringElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        #: property ID
        self.pid = pid
        #: Scalar point identification numbers
        self.nodes = nodes
        self.nodes_ref = None
        self.pid_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CELAS3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)

        s1 = integer_or_blank(card, 3, 's1', 0)
        s2 = integer_or_blank(card, 4, 's2', 0)
        assert len(card) <= 5, 'len(CELAS3 card) = %i\ncard=%s' % (len(card), card)
        return CELAS3(eid, pid, [s1, s2], comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CELAS3 card from the OP2

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
        return CELAS3(eid, pid, [s1, s2], comment=comment)

    def K(self):
        return self.pid_ref.k

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CELAS3 eid=%s' % (self.eid)
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.Pid(), msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CELAS3 eid=%s' % (self.eid)
        #self.nodes_ref = model.safe_empty_nodes(self.node_ids, msg=msg)
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    @property
    def node_ids(self):
        msg = ', which is required by CELAS3 eid=%s' % (self.eid)
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def get_edge_ids(self):
        return []

    def _verify(self, xref):
        eid = self.eid
        node_ids = self.node_ids
        s1 = self.nodes[0]
        s2 = self.nodes[1]
        #ge = self.ge
        #s = self.s

        assert isinstance(eid, int), 'eid=%r' % eid
        assert isinstance(s1, int), 's1=%r' % s1
        assert isinstance(s2, int), 's2=%r' % s2
        #assert isinstance(ge, float), 'ge=%r' % ge
        #assert isinstance(s, float), 'ge=%r' % s
        if xref:
            k = self.K()
            assert self.pid_ref.type in ['PELAS'], self.pid_ref
            assert isinstance(k, float), 'k=%r' % k
            assert len(node_ids) == len(self.nodes)
            #for nid, node in zip(node_ids, self.nodes):
                #assert node.node.nid

    def raw_fields(self):
        list_fields = ['CELAS3', self.eid, self.Pid()] + self.node_ids
        return list_fields

    #def repr_fields(self):
        #s1 = set_blank_if_default(self.s1,0)
        #s2 = set_blank_if_default(self.s2,0)
        #list_fields = ['CELAS3', self.eid, self.Pid(), s1, s2]
        #return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CELAS4(SpringElement):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CELAS4 | EID |  K  | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    type = 'CELAS4'
    _field_map = {
        1: 'eid', 2:'k', #4:'s1', 6:'s2',
    }
    cp_name_map = {'K': 'k',}

    def __init__(self, eid, k, nodes, comment=''):
        """
        Creates a CELAS4 card

        Parameters
        ----------
        eid : int
            element id
        k : float
            spring stiffness
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card
        """
        SpringElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        #: stiffness of the scalar spring
        self.k = k
        #: Scalar point identification numbers
        #self.nodes = nodes
        self.nodes = self.prepare_node_ids(nodes, allow_empty_nodes=True)
        self.nodes_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        k = []
        nodes = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            k.append(element.k)
            nodes.append([nid if nid is not None else 0 for nid in element.nodes])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('K', data=k)
        h5_file.create_dataset('nodes', data=nodes)

    def validate(self):
        if self.nodes[0] is not None and self.nodes[1] is not None:
            assert self.nodes[0] > 0 or self.nodes[1] > 0, 's1=%s s2=%s' % (self.nodes[0], self.nodes[1])

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CELAS4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        k = double(card, 2, 'k')
        s1 = integer_or_blank(card, 3, 's1', 0)
        s2 = integer_or_blank(card, 4, 's2', 0)
        assert len(card) <= 5, 'len(CELAS4 card) = %i\ncard=%s' % (len(card), card)
        return CELAS4(eid, k, [s1, s2], comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CELAS4 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        k = data[1]
        s1 = data[2]
        s2 = data[3]
        return CELAS4(eid, k, [s1, s2], comment=comment)

    def K(self):
        return self.k

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        msg = ', which is required by CELAS4 eid=%s' % (self.eid)
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def get_edge_ids(self):
        return []

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CELAS4 eid=%s' % (self.eid)
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CELAS4 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        #self.nodes_ref = model.safe_empty_nodes(self.node_ids, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def raw_fields(self):
        list_fields = ['CELAS4', self.eid, self.k] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)
