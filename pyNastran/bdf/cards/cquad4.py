# pylint: disable=C0103,R0902,R0904,R0914,C0302
"""
All shell elements are defined in this file.  This includes:

 * CQUAD4
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from numpy import array, eye, cross, allclose, float32  # zeros,dot
from numpy.linalg import det  # inv

from pyNastran.bdf.fieldWriter import (set_blank_if_default,
                                       set_default_if_blank, print_card)
from pyNastran.bdf.cards.baseCard import Element
from pyNastran.utils.mathematics import Area, norm, centroid_triangle
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)

from numpy import array

class CQUAD4s(object):
    type = 'CQUAD4'
    asterType = 'QUAD4 # CQUAD4'
    calculixType = 'S4'

    def __init__(self):
        self._comments = {}
        self._ncquad4 = 0
        self._eidmap = {}

        #: Element ID
        self.eid = []

        #: Property ID
        self.pid = []

        self.nodes = []
        self.thetaMcid = []
        self.zOffset = []
        self.TFlag = []
        self.T = []

    def add_cquad4(self, card=None, data=None, comment=''):
        #QuadShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            eid = integer(card, 1, 'eid')
            pid = integer(card, 2, 'pid')
            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2'),
                    integer(card, 5, 'n3'),
                    integer(card, 6, 'n4')]
            thetaMcid = integer_double_or_blank(card, 7, 'thetaMcid', 0.0)
            zOffset = double_or_blank(card, 8, 'zOffset', 0.0)
            blank(card, 9, 'blank')
            TFlag = integer_or_blank(card, 10, 'TFlag', 0)
            T1 = double_or_blank(card, 11, 'T1', 1.0)
            T2 = double_or_blank(card, 12, 'T2', 1.0)
            T3 = double_or_blank(card, 13, 'T3', 1.0)
            T4 = double_or_blank(card, 14, 'T4', 1.0)
            assert len(card) <= 15, 'len(CQUAD4 card) = %i' % len(card)
        else:
            eid = data[0]
            pid = data[1]
            nids = data[2:6]

            thetaMcid = data[6]
            zOffset = data[7]
            TFlag = data[8]
            T1 = data[9]
            T2 = data[10]
            T3 = data[11]
            T4 = data[12]
            if T1 == -1.0:
                T1 = 1.0
            if T2 == -1.0:
                T2 = 1.0
            if T3 == -1.0:
                T3 = 1.0
            if T4 == -1.0:
                T4 = 1.0

        #self.prepareNodeIDs(nids)
        assert len(nids) == 4, 'CQUAD4'

        self._eidmap[eid] = self._ncquad4
        self.eid.append(eid)
        self.pid.append(pid)
        self.nodes.append(nids)
        self.thetaMcid.append(thetaMcid)
        self.zOffset.append(zOffset)
        self.TFlag.append(TFlag)
        self.T.append([T1, T2, T3, T4])

        self._ncquad4 += 1

    def cross_reference(self, model):
        msg = ' which is required by CQUAD4 eid=%s' % self.eid
        #self.nodes = model.Nodes(self.nodes, msg=msg)
        #self.pid = model.Property(self.pid, msg=msg)

    def Pid(self, eid):
        id = self._eidmap[eid]
        return self.pid[id]

    def _verify(self, model, xref):
        #eid = self.Eid()
        #pid = self.Pid()
        nids = self.nodeIDs(eid)

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i,nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref == 1:  # True
            assert self.pid.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            if not self.pid.type in ['PLPLANE']:
                t = self.Thickness()
                assert isinstance(t, float), 'thickness=%r' % t
                mass = self.Mass(model)
                assert isinstance(mass, float), 'mass=%r' % mass
            a,c,n = self.AreaCentroidNormal(model)
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float32), type(c[i])
                assert isinstance(n[i], float32), type(n[i])

    def flipNormal(self, eid):
        r"""
        ::
        
          1---2       1---4
          |   |  -->  |   |
          |   |       |   |
          4---3       2---3
        """
        id = self._eidmap[eid]
        (n1, n2, n3, n4) = self.nodes[id]
        self.nodes = [n1, n4, n3, n2]

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=False)

    def writeAsCTRIA3(self, newID):
        """
        triangle - 012
        triangle - 023
        """
        zOffset = set_blank_if_default(self.zOffset, 0.0)
        nodes1 = [self.nodes[0], self.nodes[1], self.nodes[2]]
        nodes2 = [self.nodes[0], self.nodes[2], self.nodes[3]]
        fields1 = ['CTRIA3', self.eid, self.Pid()] + nodes1 + [
            self.thetaMcid, zOffset]
        fields2 = ['CTRIA3', newID, self.Pid()] + nodes2 + [
            self.thetaMcid, zOffset]
        return self.print_card(fields1) + self.print_card(fields2)

    def rawFields(self):
        list_fields = ([self.type, self.eid, self.Pid()] + self.nodeIDs() +
                  [self.thetaMcid, self.zOffset, self.TFlag, self.T1, self.T2,
                   self.T3, self.T4])
        return list_fields

    def write(self, f=None, eids=None):
        msg = []
        size = 8
        for i in xrange(self._ncquad4):
            eid = self.eid[i]
            pid = self.pid[i]
            zOffset = self.zOffset[i]
            TFlag = self.TFlag[i]
            thetaMcid = self.thetaMcid[i]
            T1, T2, T3, T4 = self.T[i]

            zOffset   = set_blank_if_default(zOffset, 0.0)
            TFlag     = set_blank_if_default(TFlag, 0)
            thetaMcid = set_blank_if_default(thetaMcid, 0.0)

            T1 = set_blank_if_default(T1, 1.0)
            T2 = set_blank_if_default(T2, 1.0)
            T3 = set_blank_if_default(T3, 1.0)
            T4 = set_blank_if_default(T4, 1.0)

            fields = (['CQUAD4', eid, pid] + self.nodes[i] +
                      [thetaMcid, zOffset, None, TFlag, T1, T2, T3, T4])

            msg.append(print_card(fields))
        return ''.join(msg)
