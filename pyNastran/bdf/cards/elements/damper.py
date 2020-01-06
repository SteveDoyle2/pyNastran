"""
All damper elements are defined in this file.  This includes:

 * CDAMP1
 * CDAMP2
 * CDAMP3
 * CDAMP4
 * CDAMP5
 * CVISC

All damper elements are DamperElement and Element objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import Element
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double)
from pyNastran.bdf.field_writer_8 import print_card_8
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class DamperElement(Element):
    def __init__(self):
        Element.__init__(self)

    #def Centroid(self):
        ## same as below, but we ignore the 2nd point it it's None
        ##p = (self.nodes_ref[1].get_position() + self.nodes_ref[0].get_position()) / 2.

        #p = self.nodes_ref[0].get_position()
        #if self.nodes_ref[1] is not None:
            #p += self.nodes_ref[1].get_position()
            #p /= 2.
        #return p

    #def center_of_mass(self):
        #return self.Centroid()

class LineDamper(DamperElement):
    def __init__(self):
        DamperElement.__init__(self)


class CDAMP1(LineDamper):
    type = 'CDAMP1'
    _field_map = {
        1: 'eid', 2:'pid', 'c1':4, 'c2':6,
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
        Creates a CDAMP1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP)
        nids : List[int, int]
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card
        """
        LineDamper.__init__(self)
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
            #components.append([comp if comp is not None else 0 for comp in [element.c1, element.c2]])
            components.append([element.c1, element.c2])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('components', data=components)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CDAMP1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer_or_blank(card, 3, 'g1', 0),
                integer_or_blank(card, 5, 'g2', 0)]

        #: component number
        c1 = integer_or_blank(card, 4, 'c1', 0)
        c2 = integer_or_blank(card, 6, 'c2', 0)
        assert len(card) <= 7, 'len(CDAMP1 card) = %i\ncard=%s' % (len(card), card)
        return CDAMP1(eid, pid, nids, c1, c2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CDAMP1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid, pid, g1, g2, c1, c2 = data
        nids = [g1, g2]
        return CDAMP1(eid, pid, nids, c1, c2, comment=comment)

    def validate(self):
        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=%r %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=%r %s' % (self.c2, msg)
        assert len(self.nodes) == 2

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, integer_types), 'nid%i is not an None/integer; nid=%s' %(i, nid)
        if xref:
            if self.pid_ref.type in 'PDAMP':
                b = self.B()
                assert isinstance(b, float)
            elif self.pid_ref.type in 'PDAMPT':
                pass
            else:
                raise NotImplementedError('pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type))

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        #return [nid if nid else None
                #for nid in self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)]
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def B(self):
        return self.pid_ref.b

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CDAMP1 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

        pid = self.pid
        if pid in model.properties:
            self.pid_ref = model.Property(pid, msg=msg)
        elif pid in model.pdampt:
            self.pid_ref = model.pdampt[pid]
        else:
            pids = model.properties.keys() + model.pdampt.keys()
            pids.sort()
            msg = ('pid=%i not found which is required by CDAMP1 eid=%i.  '
                   'Allowed Pids=%s' % (self.pid, self.eid, pids))
            raise KeyError(msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def raw_fields(self):
        nodes = self.node_ids
        fields = ['CDAMP1', self.eid, self.Pid(), nodes[0], self.c1,
                  nodes[1], self.c2]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CDAMP2(LineDamper):
    type = 'CDAMP2'
    _field_map = {
        1: 'eid', 2:'b', 'c1':4, 'c2':6,
    }
    cp_name_map = {'B' : 'b'}
    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 5:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, b, nids, c1=0, c2=0, comment=''):
        """
        Creates a CDAMP2 card

        Parameters
        ----------
        eid : int
            element id
        b : float
            damping
        nids : List[int, int]
            SPOINT ids
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card
        """
        LineDamper.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        #: Value of the scalar damper (Real)
        self.b = b
        #: component number
        self.c1 = c1
        self.c2 = c2

        # CDAMP2 do not have to be unique
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        self.nodes_ref = None
        self.pid = 0
        self.pid_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        b = []
        nodes = []
        components = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            b.append(element.b)
            nodes.append([nid if nid is not None else 0 for nid in element.nodes])
            components.append([element.c1, element.c2])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('B', data=b)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('components', data=components)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CDAMP2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        b = double(card, 2, 'b')
        nids = [integer_or_blank(card, 3, 'n1', 0),
                integer_or_blank(card, 5, 'n2', 0)]
        c1 = integer_or_blank(card, 4, 'c1', 0)
        c2 = integer_or_blank(card, 6, 'c2', 0)
        assert len(card) <= 7, 'len(CDAMP2 card) = %i\ncard=%s' % (len(card), card)
        return CDAMP2(eid, b, nids, c1, c2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CDAMP2 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        b = data[1]
        nids = [data[2], data[3]]
        c1 = data[4]
        c2 = data[5]
        return CDAMP2(eid, b, nids, c1, c2, comment=comment)

    def validate(self):
        assert len(self.nodes) == 2
        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=%r %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=%r %s' % (self.c2, msg)

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def get_edge_ids(self):
        node_ids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)
        if isinstance(node_ids[0], integer_types) and isinstance(node_ids[1], integer_types):
            return [tuple(sorted(node_ids))]
        return []

    def B(self):
        return self.b

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CDAMP2 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def _verify(self, xref):
        eid = self.eid
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, integer_types)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, integer_types), 'nid%i is not an integer/None; nid=%s' %(i, nid)

    def raw_fields(self):
        nodes = self.node_ids
        fields = ['CDAMP2', self.eid, self.b, nodes[0], self.c1,
                  nodes[1], self.c2]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CDAMP3(LineDamper):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CDAMP3 | EID | PID | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    type = 'CDAMP3'
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

    def __init__(self, eid, pid, nids, comment=''):
        """
        Creates a CDAMP3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP)
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        LineDamper.__init__(self)
        self.eid = eid
        self.pid = pid
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 2
        self.pid_ref = None
        self.nodes_ref = None

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
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CDAMP3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 's1', 0),
                integer_or_blank(card, 4, 's2', 0)]
        assert len(card) <= 5, 'len(CDAMP3 card) = %i\ncard=%s' % (len(card), card)
        return CDAMP3(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CDAMP3 card from the OP2

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
        return CDAMP3(eid, pid, nids, comment=comment)

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, integer_types), 'nid%i is not an integer/None; nid=%s' % (i, nid)
        if xref:
            assert self.pid_ref.type in ['PDAMP'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)

    def B(self):
        return self.pid_ref.b

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CDAMP3 eid=%s' % (self.eid)
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CDAMP3 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        #self.nodes_ref = model.safe_empty_nodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.pid_ref = None
        self.nodes_ref = None

    @property
    def node_ids(self):
        msg = ', which is required by CDAMP3 eid=%s' % (self.eid)
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def raw_fields(self):
        list_fields = ['CDAMP3', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CDAMP4(LineDamper):
    type = 'CDAMP4'
    _field_map = {
        1: 'eid', 2:'b',
    }
    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 4:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, b, nids, comment=''):
        """
        Creates a CDAMP4 card

        Parameters
        ----------
        eid : int
            element id
        b : float
            damping
        nids : List[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card
        """
        LineDamper.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.b = b
        self.nids = nids
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 2
        self.nodes_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        b = []
        nodes = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            b.append(element.b)
            nodes.append([nid if nid is not None else 0 for nid in element.nodes])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('B', data=b)
        h5_file.create_dataset('nodes', data=nodes)

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        ioffset = icard * 4
        eid = integer(card, 1 + ioffset, 'eid')
        #: Value of the scalar damper (Real)
        b = double(card, 2 + ioffset, 'b')
        nids = [
            integer_or_blank(card, 3 + ioffset, 'n1', 0),
            integer_or_blank(card, 4 + ioffset, 'n2', 0)
        ]
        assert len(card) <= 9, 'len(CDAMP4 card) = %i\ncard=%s' % (len(card), card)
        return CDAMP4(eid, b, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CDAMP4 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        b = data[1]
        nids = [data[2], data[3]]
        return CDAMP4(eid, b, nids, comment=comment)

    def _verify(self, xref):
        eid = self.eid
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, integer_types)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, integer_types), 'nid%i is not an integer/None; nid=%s' % (i, nid)

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        msg = ', which is required by CDAMP4 eid=%s' % (self.eid)
        nids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)
        return nids

    def B(self):
        return self.b

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CDAMP4 eid=%s' % (self.eid)
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        #msg = ', which is required by CDAMP4 eid=%s' % (self.eid)
        #self.nodes_ref = model.safe_empty_nodes(self.node_ids, msg=msg)
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def raw_fields(self):
        list_fields = ['CDAMP4', self.eid, self.b] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CDAMP5(LineDamper):
    """
    Defines a damping element that refers to a material property entry and connection to
    grid or scalar points.
    """
    type = 'CDAMP5'
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

    def __init__(self, eid, pid, nids, comment=''):
        """
        Creates a CDAMP5 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP5)
        nids : List[int, int]
            GRID/SPOINT ids
        comment : str; default=''
            a comment for the card
        """
        LineDamper.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        #: Property ID
        self.pid = pid
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 2
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
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CDAMP5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 'n1', 0),
                integer_or_blank(card, 4, 'n2', 0)]
        assert len(card) <= 5, 'len(CDAMP5 card) = %i\ncard=%s' % (len(card), card)
        return CDAMP5(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CDAMP5 card from the OP2

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
        return CDAMP5(eid, pid, nids, comment=comment)

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, integer_types), 'nid%i is not an integer/None; nid=%s' % (i, nid)
        if xref:
            assert self.pid_ref.type in ['PDAMP5'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            b = self.B()
            assert isinstance(b, float)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CDAMP5 eid=%s' % (self.eid)
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CDAMP5 eid=%s' % (self.eid)
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def B(self):
        return self.pid_ref.b

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def raw_fields(self):
        nodes = self.node_ids
        list_fields = ['CDAMP5', self.eid, self.Pid(), nodes[0], nodes[1]]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CVISC(LineDamper):
    """
    Viscous Damper Connection
    Defines a viscous damper element.

    +-------+-----+-----+----+----+
    |   1   |  2  |  3  | 4  | 5  |
    +=======+=====+=====+====+====+
    | CVISC | EID | PID | G1 | G2 |
    +-------+-----+-----+----+----+
    """
    type = 'CVISC'
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

    def __init__(self, eid, pid, nids, comment=''):
        """
        Creates a CVISC card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PVISC)
        nids : List[int, int]
            GRID ids
        comment : str; default=''
            a comment for the card
        """
        LineDamper.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 2
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
            nodes.append(element.nodes)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CVISC card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer_or_blank(card, 3, 'n1', 0),
                integer_or_blank(card, 4, 'n2', 0)]
        assert len(card) <= 5, 'len(CVISC card) = %i\ncard=%s' % (len(card), card)
        return CVISC(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CVISC card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        pid = data[1]
        nids = data[2:4]
        return CVISC(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CVISC eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CVISC eid=%s' % (self.eid)
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, integer_types), 'nid%i is not an integer/None; nid=%s' % (i, nid)
        if xref:
            assert self.pid_ref.type in ['PVISC'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)

    def B(self):
        return self.pid_ref.ce

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CVISC', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)
