# pylint: disable=C0103,R0902,R0904,R0914
"""
All spring elements are defined in this file.  This includes:

 * CELAS1
 * CELAS2
 * CELAS3
 * CELAS4

All spring elements are SpringElement and Element objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from numpy import array, zeros, dot, transpose
from numpy.linalg import norm

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Element
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
                                       double, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8


class SpringElement(Element):
    def __init__(self, card, data):
        Element.__init__(self, card, data)
        self.nodes = [None, None]

    def Eid(self):
        return self.eid

    def Centroid(self):
        p = (self.nodes[1].get_position() - self.nodes[0].get_position()) / 2.
        return p

    def Mass(self):
        return 0.0

    def repr_fields(self):
        return self.raw_fields()


class CELAS1(SpringElement):
    type = 'CELAS1'
    aster_type = 'CELAS1'
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

    def __init__(self, card=None, data=None, comment=''):
        SpringElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')

            #: property ID
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)

            nids = [integer(card, 3, 'g1'), integer_or_blank(card, 5, 'g2', 0)]
            #: component number
            self.c1 = integer_or_blank(card, 4, 'c1', 0)
            self.c2 = integer_or_blank(card, 6, 'c2', 0)
            assert len(card) <= 7, 'len(CELAS1 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = [data[2], data[3]]
            self.c1 = data[4]
            self.c2 = data[5]

        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=|%s| %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=|%s| %s' % (self.c2, msg)
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 2

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        return self._nodeIDs(allowEmptyNodes=True, msg=msg)

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def _verify(self, xref=True):
        eid = self.Eid()
        k = self.K()
        nodeIDs = self.node_ids
        c1 = self.c2
        c2 = self.c1
        #ge = self.ge
        #s = self.s

        assert isinstance(eid, int), 'eid=%r' % eid
        assert isinstance(c1, int), 'c1=%r' % c1
        assert isinstance(c2, int), 'c2=%r' % c2
        assert isinstance(k, float), 'k=%r' % k
        #assert isinstance(ge, float), 'ge=%r' % ge
        #assert isinstance(s, float), 'ge=%r' % s
        if xref:
            assert len(nodeIDs) == len(self.nodes)
            #for nodeID, node in zip(nodeIDs, self.nodes):
                #assert node.node.nid

    def _is_same_card(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid] + self.nodes + [self.pid, self.c1, self.c2]
        fields2 = [elem.eid] + elem.nodes + [elem.pid, elem.c1, elem.c2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self._is_same_fields(fields1, fields2)

    def K(self):
        return self.pid.k

    def cross_reference(self, model):
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.node_ids, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.Pid(), msg=msg)

    def uncross_reference(self, model):
        self.nodes = self.node_ids
        self.pid = self.Pid()

    def raw_fields(self):
        nodes = self.node_ids
        list_fields = ['CELAS1', self.eid, self.Pid(), nodes[0],
                       self.c1, nodes[1], self.c2]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CELAS2(SpringElement):
    type = 'CELAS2'
    aster_type = 'CELAS2'
    _field_map = {
        1: 'eid', 2:'k', 4:'c1', 6:'c2',
    }
    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 5:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        SpringElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')

            #: stiffness of the scalar spring
            self.k = double(card, 2, 'k')

            nids = [integer_or_blank(card, 3, 'g1', 0),
                    integer_or_blank(card, 5, 'g2', 0)]

            #: component number
            self.c1 = integer_or_blank(card, 4, 'c1', 0)
            self.c2 = integer_or_blank(card, 6, 'c2', 0)

            #: damping coefficient
            self.ge = double_or_blank(card, 7, 'ge', 0.)

            #: stress coefficient
            self.s = double_or_blank(card, 8, 's', 0.)
            assert len(card) <= 9, 'len(CELAS2 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.k = data[1]
            nids = [data[2], data[3]]
            self.c1 = data[4]
            self.c2 = data[5]
            self.ge = data[6]
            self.s = data[7]

        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=%r %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=%r %s' % (self.c2, msg)
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(set(self.nodes)) == 2, 'There are duplicate nodes=%s on CELAS2 eid=%s' % (self.nodes, self.eid)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        return self._nodeIDs(allowEmptyNodes=True, msg=msg)

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def cross_reference(self, model):
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.node_ids, allowEmptyNodes=True, msg=msg)

    def uncross_reference(self, model):
        self.nodes = self.node_ids

    def _verify(self, xref=True):
        eid = self.Eid()
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

    def _is_same_card(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid] + self.node_ids + [self.k, self.c1, self.c2]
        fields2 = [elem.eid] + elem.node_ids + [elem.k, elem.c1, elem.c2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self._is_same_fields(fields1, fields2)

    def K(self):
        return self.k

    def write_code_aster(self):
        nodes = self.node_ids
        msg = ''
        msg += 'DISCRET=_F( # CELAS2\n'
        if nodes[0]:
            msg += "     CARA='K_T_D_N'\n"
            msg += "     NOEUD=N%i,\n" % nodes[0]

        if nodes[1]:
            msg += "     CARA='K_T_D_L'\n"
            msg += "     NOEUD=N%i,\n" % nodes[1]
            msg += "     AMOR_HYST=%g # ge - damping\n" % self.ge
        msg += "     )\n"
        msg += "\n"

        if self.c1 == 1:
            msg += "VALE=(%g,0.,0.)\n" % self.k
        elif self.c1 == 2:
            msg += "VALE=(0.,%g,0.)\n" % self.k
        elif self.c1 == 2:
            msg += "VALE=(0.,0.,%g)\n" % self.k
        else:
            raise ValueError('unsupported value of c1=%s' % self.c1)
        return msg

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        return self._nodeIDs(allowEmptyNodes=True, msg=msg)

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

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

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CELAS3(SpringElement):
    type = 'CELAS3'
    aster_type = 'CELAS3'
    _field_map = {
        1: 'eid', 2:'pid', 4:'s1', 6:'s2',
    }

    def __init__(self, card=None, data=None, comment=''):
        SpringElement.__init__(self, card, data)

        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            #: property ID
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)

            #: Scalar point identification numbers
            self.s1 = integer_or_blank(card, 3, 's1', 0)
            self.s2 = integer_or_blank(card, 4, 's2', 0)
            assert len(card) <= 5, 'len(CELAS3 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.s1 = data[2]
            self.s2 = data[3]

    def _is_same_card(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid, self.pid, self.s1, self.s2]
        fields2 = [elem.eid, elem.pid, elem.s1, elem.s2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self._is_same_fields(fields1, fields2)

    def K(self):
        return self.pid.k

    def cross_reference(self, model):
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.node_ids, msg=msg)
        self.pid = model.Property(self.Pid(), msg=msg)

    def uncross_reference(self, model):
        self.nodes = self.node_ids
        self.pid = self.Pid()

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        return self._nodeIDs(allowEmptyNodes=True, msg=msg)

    def raw_fields(self):
        list_fields = ['CELAS3', self.eid, self.Pid(), self.s1, self.s2]
        return list_fields

    def get_edge_ids(self):
        return []

    #def repr_fields(self):
        #s1 = set_blank_if_default(self.s1,0)
        #s2 = set_blank_if_default(self.s2,0)
        #list_fields = ['CELAS3', self.eid, self.Pid(), s1, s2]
        #return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CELAS4(SpringElement):
    type = 'CELAS4'
    aster_type = 'CELAS4'
    _field_map = {
        1: 'eid', 2:'k', 4:'s1', 6:'s2',
    }

    def __init__(self, card=None, data=None, comment=''):
        SpringElement.__init__(self, card, data)

        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')

            #: stiffness of the scalar spring
            self.k = double(card, 2, 'k')

            #: Scalar point identification numbers
            self.s1 = integer_or_blank(card, 3, 's1', 0)
            self.s2 = integer_or_blank(card, 4, 's2', 0)
            assert self.s1 > 0 or self.s2 > 0, 's1=%s s2=%s' % (self.s1, self.s2)
            assert len(card) <= 5, 'len(CELAS4 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.k = data[1]
            self.s1 = data[2]
            self.s2 = data[3]

    def _is_same_card(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid, self.k, self.s1, self.s2]
        fields2 = [elem.eid, elem.k, elem.s1, elem.s2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self._is_same_fields(fields1, fields2)

    def K(self):
        return self.k

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        return self._nodeIDs(allowEmptyNodes=True, msg=msg)

    def get_edge_ids(self):
        return []

    def cross_reference(self, model):
        msg = ', which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.node_ids, allowEmptyNodes=True, msg=msg)

    def uncross_reference(self, model):
        self.nodes = self.node_ids

    def raw_fields(self):
        list_fields = ['CELAS4', self.eid, self.k, self.s1, self.s2]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)
