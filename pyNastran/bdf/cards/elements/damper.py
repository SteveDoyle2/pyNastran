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
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.utils import integer_types
from pyNastran.bdf.cards.base_card import Element
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double)
from pyNastran.bdf.field_writer_8 import print_card_8


class DamperElement(Element):
    def __init__(self):
        Element.__init__(self)


class LineDamper(DamperElement):
    def __init__(self):
        DamperElement.__init__(self)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

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

    def __init__(self, eid, pid, nids, c1, c2, comment=''):
        LineDamper.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        self.c1 = c1
        self.c2 = c2
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        self._validate_input()

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 'g1', 0),
                integer_or_blank(card, 5, 'g2', 0)]

        #: component number
        c1 = integer_or_blank(card, 4, 'c1', 0)
        c2 = integer_or_blank(card, 6, 'c2', 0)
        assert len(card) <= 7, 'len(CDAMP1 card) = %i\ncard=%s' % (len(card), card)
        return CDAMP1(eid, pid, nids, c1, c2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid, pid, g1, g2, c1, c2 = data
        nids = [g1, g2]
        return CDAMP1(eid, pid, nids, c1, c2, comment=comment)

    def _validate_input(self):
        assert len(self.nodes) == 2
        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=%r %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=%r %s' % (self.c2, msg)

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, integer_types), 'nid%i is not an None/integer; nid=%s' %(i, nid)
        if xref:
            assert self.pid_ref.type in ['PDAMP'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)

    def Eid(self):
        return self.eid

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def _is_same_card(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid, self.Pid()] + self.node_ids + [self.c1, self.c2]
        fields2 = [elem.eid, elem.Pid()] + elem.node_ids + [elem.c1, elem.c2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self._is_same_fields(fields1, fields2)

    def B(self):
        return self.pid_ref.b

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CDAMP1 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes

        pid = self.pid
        if pid in model.properties:
            self.pid = model.Property(pid, msg=msg)
            self.pid_ref = self.pid
        elif pid in model.pdampt:
            self.pid = model.pdampt[pid]
            self.pid_ref = self.pid
        else:
            pids = model.properties.keys() + model.pdampt.keys()
            pids.sort()
            msg = ('pid=%i not found which is required by CDAMP1 eid=%i.  '
                   'Allowed Pids=%s' % (self.pid, self.eid, pids))
            raise KeyError(msg)

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def raw_fields(self):
        nodes = self.node_ids
        fields = ['CDAMP1', self.eid, self.Pid(), nodes[0], self.c1,
                  nodes[1], self.c2]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CDAMP2(LineDamper):
    type = 'CDAMP2'
    _field_map = {
        1: 'eid', 2:'b', 'c1':4, 'c2':6,
    }
    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 5:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, b, nids, c1, c2, comment=''):
        LineDamper.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        #: Value of the scalar damper (Real)
        self.b = b

        #: component number
        self.c1 = c1
        self.c2 = c2

        # CDAMP2 do not have to be unique
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        self._validate_input()

    @classmethod
    def add_card(cls, card, comment=''):
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
        eid = data[0]
        b = data[1]
        nids = [data[2], data[3]]
        c1 = data[4]
        c2 = data[5]
        return CDAMP2(eid, b, nids, c1, c2, comment=comment)

    def _validate_input(self):
        assert len(self.nodes) == 2
        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=%r %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=%r %s' % (self.c2, msg)

    def _verify(self, xref=True):
        eid = self.Eid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer/None; nid=%s' %(i, nid)

    def Eid(self):
        return self.eid

    def B(self):
        return self.b

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CDAMP2 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes

    def uncross_reference(self):
        self.nodes = self.node_ids
        del self.nodes_ref

    def get_edge_ids(self):
        node_ids = self._nodeIDs(allow_empty_nodes=True)
        if isinstance(node_ids[0], integer_types) and isinstance(node_ids[1], integer_types):
            return [tuple(sorted(node_ids))]
        return []

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def raw_fields(self):
        nodes = self.node_ids
        fields = ['CDAMP2', self.eid, self.b, nodes[0], self.c1,
                  nodes[1], self.c2]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CDAMP3(LineDamper):
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
        if comment:
            self._comment = comment
        LineDamper.__init__(self)
        self.eid = eid
        self.pid = pid
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 2

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 's1', 0),
                integer_or_blank(card, 4, 's2', 0)]
        assert len(card) <= 5, 'len(CDAMP3 card) = %i\ncard=%s' % (len(card), card)
        return CDAMP3(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = [data[2], data[3]]
        return CDAMP3(eid, pid, nids, comment=comment)

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer/None; nid=%s' % (i, nid)
        if xref:
            assert self.pid_ref.type in ['PDAMP'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)

    def Eid(self):
        return self.eid

    def B(self):
        return self.pid_ref.b

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CDAMP3', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def write_card(self, size=8, is_double=False):
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
        LineDamper.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.b = b
        self.nids = nids
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 2

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
        eid = data[0]
        b = data[1]
        nids = [data[2], data[3]]
        return CDAMP4(eid, b, nids, comment=comment)

    def _verify(self, xref=True):
        eid = self.Eid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer/None; nid=%s' % (i, nid)

    def Eid(self):
        return self.eid

    def B(self):
        return self.b

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CDAMP4', self.eid, self.b] + self.node_ids
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CDAMP5(LineDamper):
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
        LineDamper.__init__(self)
        if comment:
            self._comment = comment

        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 2

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 'n1', 0),
                integer_or_blank(card, 4, 'n2', 0)]
        assert len(card) <= 5, 'len(CDAMP5 card) = %i\ncard=%s' % (len(card), card)
        return CDAMP5(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = [data[2], data[3]]
        return CDAMP5(eid, pid, nids, comment=comment)

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer/None; nid=%s' % (i, nid)
        if xref:
            assert self.pid_ref.type in ['PDAMP5'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def Eid(self):
        return self.eid

    def B(self):
        return self.pid_ref.b

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def raw_fields(self):
        nodes = self.node_ids
        list_fields = ['CDAMP5', self.eid, self.Pid(), nodes[0], nodes[1]]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class CVISC(LineDamper):
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
        LineDamper.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 2

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer_or_blank(card, 3, 'n1', 0),
                integer_or_blank(card, 4, 'n2', 0)]
        assert len(card) <= 5, 'len(CVISC card) = %i\ncard=%s' % (len(card), card)
        return CVISC(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = data[2:4]
        return CVISC(eid, pid, nids, comment=comment)

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer/None; nid=%s' % (i, nid)
        if xref:
            assert self.pid_ref.type in ['PVISC'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)

    def Eid(self):
        return self.eid

    def B(self):
        return self.pid_ref.ce

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CVISC', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)
