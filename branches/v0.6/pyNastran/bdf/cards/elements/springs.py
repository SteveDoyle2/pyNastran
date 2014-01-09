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

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Element
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
                                       double, double_or_blank)
from pyNastran.bdf.fieldWriter import print_card_8


class SpringElement(Element):
    def __init__(self, card, data):
        Element.__init__(self, card, data)
        self.nodes = [None, None]

    def Eid(self):
        return self.eid

    def Centroid(self):
        p = (self.nodes[1].Position() - self.nodes[0].Position()) / 2.
        return p

    def K(self):
        raise NotImplementedError('K not implemented in the %s class'
                                  % self.type)

    def Lambda(self, model, debug=False):
        """
          3d  [l,m,n,0,0,0]
              [0,0,0,l,m,n]
        """
        (n1, n2) = self.nodeIDs()
        p1 = model.Node(n1).Position()
        p2 = model.Node(n2).Position()
        v1 = p2 - p1
        v1 = v1 / norm(v1)
        (l, m, n) = v1
        Lambda = array(zeros((2, 6), 'float64'))

        Lambda[0, 0] = Lambda[1, 3] = l
        Lambda[0, 1] = Lambda[1, 4] = m
        Lambda[0, 2] = Lambda[1, 5] = n
        return Lambda

    def Mass(self):
        return 0.0

    def reprFields(self):
        return self.rawFields()

    def __repr__(self):
        return self.print_card(8)


class CELAS1(SpringElement):
    type = 'CELAS1'
    asterType = 'CELAS1'
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
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 2

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=True)

    def _verify(self, xref=False):
        eid = self.Eid()
        k = self.K()
        nodeIDs = self.nodeIDs()
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

    def isSameCard(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid] + self.nodes + [self.pid, self.c1, self.c2]
        fields2 = [elem.eid] + elem.nodes + [elem.pid, elem.c1, elem.c2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def K(self):
        return self.pid.k

    def cross_reference(self, model):
        msg = ' which is required by CELAS1 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def rawFields(self):
        nodes = self.nodeIDs()
        list_fields = ['CELAS1', self.eid, self.Pid(), nodes[0],
                  self.c1, nodes[1], self.c2]
        return list_fields

    def Stiffness(self, model, node_ids, index0s, fnorm):
        ki = self.K()
        k = ki * array([[1, -1,],
                        [-1, 1]])
        Lambda = self.Lambda(model)
        K = dot(dot(transpose(Lambda), k), Lambda)

        c1 = self.c1
        c2 = self.c2
        n1, n2 = node_ids
        delta1 = 0 if c1 in [1, 2, 3] else 3
        delta2 = 0 if c2 in [1, 2, 3] else 3

        nIJV = [
            (n1, 1 + delta1), (n1, 2 + delta1), (n1, 3 + delta1),
            (n2, 1 + delta2), (n2, 2 + delta2), (n2, 3 + delta2),
        ]
        dofs = nIJV
        return (K, dofs, nIJV)

    def displacement_stress(self, model, q, dofs, is3D=False):
        (n1, n2) = self.nodeIDs()
        Lambda = self.Lambda(model, debug=False)

        n11 = dofs[(n1, 1)]
        n21 = dofs[(n2, 1)]

        n12 = dofs[(n1, 2)]
        n22 = dofs[(n2, 2)]

        n13 = dofs[(n1, 3)]
        n23 = dofs[(n2, 3)]

        q_axial = array([
            q[n11], q[n12], q[n13],
            q[n21], q[n22], q[n23]
        ])
        u_axial = dot(array(Lambda), q_axial)
        du_axial = u_axial[0] - u_axial[1]

        s = self.pid.s
        ki = self.pid.k

        axial_strain = du_axial * s
        axial_force = ki * du_axial
        axial_stress = axial_force * s

        return (axial_strain, axial_stress, axial_force)

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + print_card_8(card)


class CELAS2(SpringElement):
    type = 'CELAS2'
    asterType = 'CELAS2'
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
            nids = [data[2], data[4]]
            self.c1 = data[3]
            self.c2 = data[5]
            self.ge = data[6]
            self.s = data[7]

        msg = 'on\n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % str(self)
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=|%s| %s' % (self.c1, msg)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=|%s| %s' % (self.c2, msg)
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(set(self.nodes)) == 2, 'There are duplicate nodes=%s on CELAS2 eid=%s' % (self.nodes, self.eid)

    def cross_reference(self, model):
        msg = ' which is required by CELAS2 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)

    def _verify(self, xref=False):
        eid = self.Eid()
        k = self.K()
        nodeIDs = self.nodeIDs()
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
            assert len(nodeIDs) == len(self.nodes)
            #for nodeID, node in zip(nodeIDs, self.nodes):
                #assert node.node.nid

    def isSameCard(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid] + self.nodes + [self.k, self.c1, self.c2]
        fields2 = [elem.eid] + elem.nodes + [elem.k, elem.c1, elem.c2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def K(self):
        return self.k

    def writeCodeAster(self):
        nodes = self.nodeIDs()
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
        return self._nodeIDs(allowEmptyNodes=True,
                             msg=str(['CELAS2', self.eid]))

    def rawFields(self):
        nodes = self.nodeIDs()
        list_fields = ['CELAS2', self.eid, self.k, nodes[0], self.c1,
                  nodes[1], self.c2, self.ge, self.s]
        return list_fields

    def reprFields(self):
        nodes = self.nodeIDs()
        ge = set_blank_if_default(self.ge, 0.)
        s = set_blank_if_default(self.s, 0.)
        list_fields = ['CELAS2', self.eid, self.k, nodes[0], self.c1,
                  nodes[1], self.c2, ge, s]
        return list_fields


    def Stiffness(self, model, node_ids, index0s, fnorm):
        ki = self.k
        k = ki * array([[1, -1,],
                        [-1, 1]])
        Lambda = self.Lambda(model)
        K = dot(dot(transpose(Lambda), k), Lambda)

        c1 = self.c1
        c2 = self.c2
        n1, n2 = node_ids
        delta1 = 0 if c1 in [1, 2, 3] else 3
        delta2 = 0 if c2 in [1, 2, 3] else 3

        nIJV = [
            (n1, 1 + delta1), (n1, 2 + delta1), (n1, 3 + delta1),
            (n2, 1 + delta2), (n2, 2 + delta2), (n2, 3 + delta2),
        ]
        dofs = nIJV
        return (K, dofs, nIJV)

    def displacement_stress(self, model, q, dofs):
        (n1, n2) = self.nodeIDs()
        Lambda = self.Lambda(model, debug=False)

        n11 = dofs[(n1, 1)]
        n21 = dofs[(n2, 1)]

        n12 = dofs[(n1, 2)]
        n22 = dofs[(n2, 2)]

        n13 = dofs[(n1, 3)]
        n23 = dofs[(n2, 3)]

        q_axial = array([
            q[n11], q[n12], q[n13],
            q[n21], q[n22], q[n23]
        ])
        u_axial = dot(Lambda, q_axial)
        du_axial = u_axial[0] - u_axial[1]

        s = self.s
        ki = self.k

        axial_strain = du_axial * s
        axial_force = ki * du_axial
        axial_stress = axial_force * s

        return (axial_strain, axial_stress, axial_force)

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + print_card_8(card)


class CELAS3(SpringElement):
    type = 'CELAS3'
    asterType = 'CELAS3'
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

    def isSameCard(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid, self.pid, self.s1, self.s2]
        fields2 = [elem.eid, elem.pid, elem.s1, elem.s2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def K(self):
        return self.pid.k

    def cross_reference(self, model):
        msg = ' which is required by CELAS3 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def rawFields(self):
        list_fields = ['CELAS3', self.eid, self.Pid(), self.s1, self.s2]
        return list_fields

    def Stiffness(self, model, node_ids, index0s, fnorm):
        ki = self.pid.k
        k = ki * array([[1, -1,],
                        [-1, 1]])

        c1 = self.c1
        c2 = self.c2
        n1, n2 = node_ids

        nIJV = [
            (n1, c1),
            (n2, c2),
        ]
        dofs = nIJV
        return (K, dofs, nIJV)

    def displacement_stress(self, model, q, dofs):
        (n1, n2) = self.nodeIDs()

        #print("**dofs =", dofs)
        n11 = dofs[(n1, self.c1)]
        n21 = dofs[(n2, self.c2)]

        u = array([
            q[n11],
            q[n21],
        ])
        du = u[0] - u[1]

        s = self.s
        axial_strain = du * s

        ki = self.pid.k
        axial_force = ki * du
        axial_stress = axial_force * s
        return (axial_strain, axial_stress, axial_force)

    #def reprFields(self):
        #s1 = set_blank_if_default(self.s1,0)
        #s2 = set_blank_if_default(self.s2,0)
        #list_fields = ['CELAS3',self.eid,self.Pid(),s1,s2]
        #return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + print_card_8(card)


class CELAS4(SpringElement):
    type = 'CELAS4'
    asterType = 'CELAS4'
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

    def isSameCard(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid, self.k, self.s1, self.s2]
        fields2 = [elem.eid, elem.k, elem.s1, elem.s2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def K(self):
        return self.k

    def cross_reference(self, model):
        msg = ' which is required by CELAS4 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)

    def rawFields(self):
        list_fields = ['CELAS4', self.eid, self.k, self.s1, self.s2]
        return list_fields

    def Stiffness(self, model, node_ids, index0s, fnorm):
        ki = self.k
        k = ki * array([[1, -1,],
                        [-1, 1]])

        c1 = self.c1
        c2 = self.c2
        n1, n2 = node_ids

        nIJV = [
            (n1, c1),
            (n2, c2),
        ]
        dofs = nIJV
        return (K, dofs, nIJV)

    def displacement_stress(self, model, q, dofs):
        (n1, n2) = self.nodeIDs()

        #print("**dofs =", dofs)
        n11 = dofs[(n1, self.c1)]
        n21 = dofs[(n2, self.c2)]

        u = array([
            q[n11],
            q[n21],
        ])
        du = u[0] - u[1]

        s = 0.0
        axial_strain = du * s

        ki = self.k
        axial_force = ki * du
        axial_stress = axial_force * s
        return (axial_strain, axial_stress, axial_force)

    #def reprFields(self):
        #s1 = set_blank_if_default(self.s1,0)
        #s2 = set_blank_if_default(self.s2,0)
        #list_fields = ['CELAS4',self.eid,self.Pid(),s1,s2]
        #return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + print_card_8(card)
