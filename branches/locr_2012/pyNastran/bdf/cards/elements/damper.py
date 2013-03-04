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
from pyNastran.bdf.format import integer, integer_or_blank, double

class DamperElement(Element):
    def __init__(self, card, data):
        Element.__init__(self, card, data)


class LineDamper(DamperElement):
    def __init__(self, card, data):
        DamperElement.__init__(self, card, data)

    def cross_reference(self, mesh):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = mesh.Nodes(self.nodes, msg=msg)
        self.pid = mesh.Property(self.pid, msg=msg)


class CDAMP1(LineDamper):
    type = 'CDAMP1'

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer(card, 2, 'pid')
            nids = [integer_or_blank(card, 3, 'g1', 0),
                    integer_or_blank(card, 5, 'g2', 0)]

            ## component number
            self.c1 = integer_or_blank(card, 4, 'c1', 0)
            self.c2 = integer_or_blank(card, 6, 'c2', 0)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = [data[2], data[4]]
            self.c1 = data[3]
            self.c2 = data[5]

        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 2
        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=|%s| %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=|%s| %s' % (self.c2, msg)

    def isSameCard(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid, self.Pid()] + self.nodeIDs() + [self.c1, self.c2]
        fields2 = [elem.eid, elem.Pid()] + elem.nodeIDs() + [elem.c1, elem.c2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def B(self):
        return self.pid.b

    def cross_reference(self, model):
        msg = ' which is required by CDAMP1 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        fields = ['CDAMP1', self.eid, self.Pid(), nodes[0], self.c1,
                  nodes[1], self.c2]
        return fields


class CDAMP2(LineDamper):
    type = 'CDAMP2'

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')

            ## Value of the scalar damper (Real)
            self.b = double(card, 2, 'b')
            nids = [integer_or_blank(card, 3, 'n1', 0),
                    integer_or_blank(card, 5, 'n2', 0)]

            ## component number
            self.c1 = integer_or_blank(card, 4, 'c1', 0)
            self.c2 = integer_or_blank(card, 6, 'c2', 0)
        else:
            self.eid = data[0]
            self.b = data[1]
            nids = [data[2], data[4]]
            self.c1 = data[3]
            self.c2 = data[5]

        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 2
        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=|%s| %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=|%s| %s' % (self.c2, msg)

    def B(self):
        return self.b

    def cross_reference(self, model):
        msg = ' which is required by CDAMP2 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)

    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        fields = ['CDAMP2', self.eid, self.b, nodes[0], self.c1,
                  nodes[1], self.c2]
        return fields


class CDAMP3(LineDamper):
    type = 'CDAMP3'

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer(card, 2, 'pid')
            nids = [integer_or_blank(card, 3, 'n1', 0),
                    integer_or_blank(card, 5, 'n2', 0)]
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = [data[2], data[3]]
        self.f(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 2

    def B(self):
        return self.pid.b

    def cross_reference(self, model):
        msg = ' which is required by CDAMP4 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        list_fields = ['CDAMP3', self.eid, self.pid, nodes[0], nodes[1]]
        return list_fields


class CDAMP4(LineDamper):
    type = 'CDAMP4'

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            ## Value of the scalar damper (Real)
            self.b = double(card, 2, 'b')
            nids = [integer_or_blank(card, 3, 'n1', 0),
                    integer_or_blank(card, 4, 'n2', 0)]
        else:
            self.eid = data[0]
            self.b = data[1]
            nids = [data[2], data[3]]
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 2

    def B(self):
        return self.b

    def cross_reference(self, model):
        msg = ' which is required by CDAMP4 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)

    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        list_fields = ['CDAMP4', self.eid, self.b, nodes[0], nodes[1]]
        return list_fields


class CDAMP5(LineDamper):
    type = 'CDAMP5'

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            ## Property ID
            self.pid = integer(card, 2, 'pid')
            nids = [integer_or_blank(card, 3, 'n1', 0),
                    integer_or_blank(card, 4, 'n2', 0)]
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = [data[2], data[3]]
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 2

    def cross_reference(self, model):
        msg = ''
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        list_fields = ['CDAMP5', self.eid, self.Pid(), nodes[0], nodes[1]]
        return list_fields


class CVISC(LineDamper):
    type = 'CVISC'

    def __init__(self, card=None, data=None, comment=''):
        LineDamper.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            nids = [integer_or_blank(card, 3, 'n1', 0),
                    integer_or_blank(card, 4, 'n2', 0)]
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:4]
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 2

    def B(self):
        return self.pid.ce

    def rawFields(self):
        list_fields = ['CVISC', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields

    def reprFields(self):
        return self.rawFields()