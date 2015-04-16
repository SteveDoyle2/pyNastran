# pylint: disable=C0103,R0902,R0904,R0914,C0111
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

from pyNastran.bdf.cards.baseCard import Element
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
                                                    double)
from pyNastran.bdf.field_writer_8 import print_card_8


class DamperElement(Element):
    def __init__(self, card, data):
        Element.__init__(self, card, data)


class LineDamper(DamperElement):
    def __init__(self, card, data):
        DamperElement.__init__(self, card, data)

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)


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

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer(card, 2, 'pid')
            nids = [integer_or_blank(card, 3, 'g1', 0),
                    integer_or_blank(card, 5, 'g2', 0)]

            #: component number
            self.c1 = integer_or_blank(card, 4, 'c1', 0)
            self.c2 = integer_or_blank(card, 6, 'c2', 0)
            assert len(card) <= 7, 'len(CDAMP1 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = [data[2], data[4]]
            self.c1 = data[3]
            self.c2 = data[5]

        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 2
        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=|%s| %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=|%s| %s' % (self.c2, msg)

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
        if xref:
            assert self.pid.type in ['PDAMP'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)

    def Eid(self):
        return self.eid

    def nodeIDs(self):
        """deprecated"""
        return self.node_ids

    @property
    def node_ids(self):
        return [0 if nid is None else nid for nid in self._nodeIDs(allowEmptyNodes=True)]

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def _is_same_card(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid, self.Pid()] + self.node_ids + [self.c1, self.c2]
        fields2 = [elem.eid, elem.Pid()] + elem.nodeIDs() + [elem.c1, elem.c2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self._is_same_fields(fields1, fields2)

    def B(self):
        return self.pid.b

    def cross_reference(self, model):
        msg = ' which is required by CDAMP1 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        pid = self.pid
        if pid in model.properties:
            self.pid = model.Property(pid, msg=msg)
        elif pid in model.pdampt:
            self.pid = model.pdampt[pid]
        else:
            pids = model.properties.keys() + model.pdampt.keys()
            pids.sort()
            msg = 'pid=%i not found which is required by CDAMP1 eid=%i.  Allowed Pids=%s' % (self.pid, self.eid, pids)
            raise KeyError(msg)

    def raw_fields(self):
        nodes = self.node_ids
        fields = ['CDAMP1', self.eid, self.Pid(), nodes[0], self.c1,
                  nodes[1], self.c2]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment() + print_card_8(card)


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

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')

            #: Value of the scalar damper (Real)
            self.b = double(card, 2, 'b')
            nids = [integer_or_blank(card, 3, 'n1', 0),
                    integer_or_blank(card, 5, 'n2', 0)]

            #: component number
            self.c1 = integer_or_blank(card, 4, 'c1', 0)
            self.c2 = integer_or_blank(card, 6, 'c2', 0)
            assert len(card) <= 7, 'len(CDAMP2 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.b = data[1]
            nids = [data[2], data[4]]
            self.c1 = data[3]
            self.c2 = data[5]

        self.prepare_node_ids(nids, allow_empty_nodes=True)
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
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

    def Eid(self):
        return self.eid

    def B(self):
        return self.b

    def cross_reference(self, model):
        msg = ' which is required by CDAMP2 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)

    def nodeIDs(self):
        return self.node_ids

    @property
    def node_ids(self):
        return [0 if nid is None else nid for nid in self._nodeIDs(allowEmptyNodes=True)]

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        nodes = self.node_ids
        fields = ['CDAMP2', self.eid, self.b, nodes[0], self.c1,
                  nodes[1], self.c2]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment() + print_card_8(card)


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

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer(card, 2, 'pid')
            nids = [integer_or_blank(card, 3, 's1', 0),
                    integer_or_blank(card, 4, 's2', 0)]
            assert len(card) <= 5, 'len(CDAMP3 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = [data[2], data[3]]
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 2

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
        if xref:
            assert self.pid.type in ['PDAMP'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)

    def Eid(self):
        return self.eid

    def B(self):
        return self.pid.b

    def cross_reference(self, model):
        msg = ' which is required by CDAMP3 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def nodeIDs(self):
        """deprecated"""
        return self.node_ids

    @property
    def node_ids(self):
        return [0 if nid is None else nid for nid in self._nodeIDs(allowEmptyNodes=True) ]

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        nodes = self.node_ids
        list_fields = ['CDAMP3', self.eid, self.pid, nodes[0], nodes[1]]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment() + print_card_8(card)


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

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            #: Value of the scalar damper (Real)
            self.b = double(card, 2, 'b')
            nids = [integer_or_blank(card, 3, 'n1', 0),
                    integer_or_blank(card, 4, 'n2', 0)]
            assert len(card) <= 5, 'len(CDAMP4 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.b = data[1]
            nids = [data[2], data[3]]
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 2

    def _verify(self, xref=True):
        eid = self.Eid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

    def Eid(self):
        return self.eid

    def B(self):
        return self.b

    def cross_reference(self, model):
        msg = ' which is required by CDAMP4 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)

    def nodeIDs(self):
        """deprecated"""
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=True)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        nodes = self.node_ids
        list_fields = ['CDAMP4', self.eid, self.b, nodes[0], nodes[1]]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment() + print_card_8(card)


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

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            #: Property ID
            self.pid = integer(card, 2, 'pid')
            nids = [integer_or_blank(card, 3, 'n1', 0),
                    integer_or_blank(card, 4, 'n2', 0)]
            assert len(card) <= 5, 'len(CDAMP5 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = [data[2], data[3]]
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 2

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
        if xref:
            assert self.pid.type in ['PDAMP5'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)

    def cross_reference(self, model):
        msg = ''
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def Eid(self):
        return self.eid

    def B(self):
        return self.pid.b

    def nodeIDs(self):
        """deprecated"""
        return self.node_ids

    @property
    def node_ids(self):
        return [0 if nid is None else nid for nid in self._nodeIDs(allowEmptyNodes=True)]

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        nodes = self.node_ids
        list_fields = ['CDAMP5', self.eid, self.Pid(), nodes[0], nodes[1]]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment() + print_card_8(card)


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

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            nids = [integer_or_blank(card, 3, 'n1', 0),
                    integer_or_blank(card, 4, 'n2', 0)]
            assert len(card) <= 5, 'len(CVISC card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:4]
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 2

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        b = self.B()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        assert isinstance(b, float)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
        if xref:
            assert self.pid.type in ['PVISC'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)

    def Eid(self):
        return self.eid

    def B(self):
        return self.pid.ce

    def nodeIDs(self):
        """deprecated"""
        return self.node_ids

    @property
    def node_ids(self):
        return [0 if nid is None else nid for nid in self._nodeIDs(allowEmptyNodes=True)]

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = ['CVISC', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment() + print_card_8(card)
