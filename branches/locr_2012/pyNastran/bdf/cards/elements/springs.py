# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from numpy import matrix, zeros, dot, transpose
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Element


class SpringElement(Element):
    def __init__(self, card, data):
        Element.__init__(self, card, data)

    def Length_noXref(self, n1=None, n2=None):
        r"""
        Returns the length of a bar/rod/beam element
        \f[ \large \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  } \f]
        @param self the object pointer
        @param n1 a Node object (default=None)
        @param n2 a Node object (default=None)
        @note
            if n1 AND n2 are both none (the default), then the model must
            be cross-referenced already
        """
        #print self.type
        L = norm(n1.Position() - n2.Position())
        return L

    def Centroid(self):
        p = (self.nodes[1].Position() - self.nodes[0].Position()) / 2.
        return p

    def K(self):
        raise NotImplementedError('K not implemented in the '
                                  '%s class' % (self.type))

    def Lambda(self, model):
        """
        2d  [l,m,0,0]
            [0,0,l,m]

        3d  [l,m,n,0,0,0]
            [0,0,0,l,m,n]
        """
        is3D = False
        #R = self.Rmatrix(model,is3D)

        (n1, n2) = self.nodeIDs()
        p1 = model.Node(n1).Position()
        p2 = model.Node(n2).Position()
        v1 = p2 - p1
        #print(v1)
        v1 = v1 / norm(v1)
        (l, m, n) = v1
        if is3D:
            Lambda = matrix(zeros((2, 6), 'd'))  # 3D
        else:
            Lambda = matrix(zeros((2, 4), 'd'))

        #print("R = \n",R)
        Lambda[0, 0] = Lambda[1, 2] = l
        Lambda[0, 1] = Lambda[1, 3] = m

        if is3D:
            Lambda[0, 2] = Lambda[1, 5] = n  # 3D
        #print("Lambda = \n",Lambda)
        return Lambda

    def Stiffness(self, model):
        ki = self.K()
        k = ki * matrix([[1, -1, ]
                         [-1, 1]])

        Lambda = self.Lambda(model)
        K = dot(dot(transpose(Lambda), k), Lambda)
        return K

    def Length(self):
        r"""
        Returns the length of a bar/rod/beam element
        \f[ \large \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  } \f]
        @param self the object pointer
        @note
            the model must be cross-referenced already
        """
        #print self.type
        return self.Length_noXref(self.nodes[1], self.nodes[0])

    def Mass(self):
        return 0.0

    def reprFields(self):
        return self.rawFields()

    def __repr__(self):
        fields = self.rawFields()
        return self.printCard(fields)


class CELAS1(SpringElement):
    type = 'CELAS1'
    asterType = 'CELAS1'

    def __init__(self, card=None, data=None):
        SpringElement.__init__(self, card, data)
        if card:
            self.eid = card.field(1)

            ## property ID
            self.pid = card.field(2, self.eid)

            nids = [card.field(3, 0), card.field(5, 0)]
            ## component number
            self.c1 = card.field(4, 0)
            self.c2 = card.field(6, 0)

        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = [data[2], data[3]]
            self.c1 = data[4]
            self.c2 = data[5]
        ###
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % (str(self), self.c1)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % (str(self), self.c2)
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 2

    def isSameCard(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid] + self.nodes + [self.pid, self.c1, self.c2]
        fields2 = [elem.eid] + self.nodes + [elem.pid, elem.c1, elem.c2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def K(self):
        return self.pid.k

    def cross_reference(self, model):
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True)
        self.pid = model.Property(self.pid)

    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True)
        fields = ['CELAS1', self.eid, self.Pid(), nodes[0],
                  self.c1, nodes[1], self.c2]
        return fields


class CELAS2(SpringElement):
    type = 'CELAS2'
    asterType = 'CELAS2'

    def __init__(self, card=None, data=None):
        SpringElement.__init__(self, card, data)

        if card:
            self.eid = card.field(1)

            ## stiffness of the scalar spring
            self.k = card.field(2)

            nids = [card.field(3, 0), card.field(5, 0)]

            ## component number
            self.c1 = card.field(4, 0)
            self.c2 = card.field(6, 0)

            ## damping coefficient
            self.ge = card.field(7, 0.)

            ## stress coefficient
            self.s = card.field(8, 0.)
        else:
            self.eid = data[0]
            self.k = data[1]
            nids = [data[2], data[4]]
            self.c1 = data[3]
            self.c2 = data[5]
            self.ge = data[6]
            self.s = data[7]
        ###
        assert self.c1 in [0, 1, 2, 3, 4, 5, 6], 'c1=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % (str(self), self.c1)
        assert self.c2 in [0, 1, 2, 3, 4, 5, 6], 'c2=|%s| on \n%s\n is invalid validComponents=[0,1,2,3,4,5,6]' % (str(self), self.c2)
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 2

    def isSameCard(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid] + self.nodes + [self.k, self.c1, self.c2]
        fields2 = [elem.eid] + self.nodes + [elem.k, elem.c1, elem.c2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def K(self):
        return self.k

    def cross_reference(self, model):
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True)
        #print("nodes = ",self.nodes)

    def writeCodeAster(self):
        nodes = self.nodeIDs()
        msg = ''
        msg += 'DISCRET=_F( # CELAS2\n'
        if nodes[0]:
            msg += "     CARA='K_T_D_N'\n"
            msg += "     NOEUD=N%i,\n" % (nodes[0])

        if nodes[1]:
            msg += "     CARA='K_T_D_L'\n"
            msg += "     NOEUD=N%i,\n" % (nodes[1])
            msg += "     AMOR_HYST=%g # ge - damping\n" % (self.ge)
        msg += "     )\n"
        msg += "\n"

        if self.c1 == 1:
            msg += "VALE=(%g,0.,0.)\n" % (self.k)
        elif self.c1 == 2:
            msg += "VALE=(0.,%g,0.)\n" % (self.k)
        elif self.c1 == 2:
            msg += "VALE=(0.,0.,%g)\n" % (self.k)
        else:
            raise ValueError('unsupported value of c1=%s' % (self.c1))
        ###
        return msg

    def rawFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True,
                             msg=str(['CELAS2', self.eid]))
        fields = ['CELAS2', self.eid, self.k, nodes[0], self.c1,
                  nodes[1], self.c2, self.ge, self.s]
        return fields

    def reprFields(self):
        nodes = self.nodeIDs(allowEmptyNodes=True,
                             msg=str(['CELAS2', self.eid]))
        ge = set_blank_if_default(self.ge, 0.)
        s = set_blank_if_default(self.s, 0.)
        fields = ['CELAS2', self.eid, self.k, nodes[0], self.c1,
                  nodes[1], self.c2, ge, s]
        return fields


class CELAS3(SpringElement):
    type = 'CELAS3'
    asterType = 'CELAS3'

    def __init__(self, card=None, data=None):
        SpringElement.__init__(self, card, data)

        if card:
            #nids = [card.field(3),card.field(5)]
            #self.prepareNodeIDs(nids)
            #assert len(self.nodes)==2

            self.eid = card.field(1)
            ## property ID
            self.pid = card.field(2, self.eid)

            ## Scalar point identification numbers
            self.s1 = card.field(3, 0)
            self.s2 = card.field(4, 0)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.s1 = data[2]
            self.s2 = data[3]
        ###

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
        self.nodes = model.Nodes(self.nodes)
        self.pid = model.Property(self.pid)

    def rawFields(self):
        fields = ['CELAS3', self.eid, self.Pid(), self.s1, self.s2]
        return fields

    #def reprFields(self):
        #s1 = set_blank_if_default(self.s1,0)
        #s2 = set_blank_if_default(self.s2,0)
        #fields = ['CELAS3',self.eid,self.Pid(),s1,s2]
        #return fields


class CELAS4(SpringElement):
    type = 'CELAS4'
    asterType = 'CELAS4'

    def __init__(self, card=None, data=None):
        SpringElement.__init__(self, card, data)

        if card:
            #nids = [card.field(3),card.field(5)]
            #self.prepareNodeIDs(nids)
            #assert len(self.nodes)==2

            self.eid = card.field(1)

            ## stiffness of the scalar spring
            self.k = card.field(2)

            ## Scalar point identification numbers
            self.s1 = card.field(3, 0)
            self.s2 = card.field(4, 0)
        else:
            self.eid = data[0]
            self.k = data[1]
            self.s1 = data[2]
            self.s2 = data[3]
        ###

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
        self.nodes = model.Nodes(self.nodes)

    def rawFields(self):
        fields = ['CELAS4', self.eid, self.k, self.s1, self.s2]
        return fields

    #def reprFields(self):
        #s1 = set_blank_if_default(self.s1,0)
        #s2 = set_blank_if_default(self.s2,0)
        #fields = ['CELAS4',self.eid,self.Pid(),s1,s2]
        #return fields
