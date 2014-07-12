# pylint: disable=C0103,R0902,R0904,R0914,C0302
"""
All shell elements are defined in this file.  This includes:

 * CQUAD4
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import array, eye, cross, allclose, float32  # zeros,dot
from numpy.linalg import det  # inv

from pyNastran.bdf.fieldWriter import (set_blank_if_default,
                                       set_default_if_blank, print_card)
from pyNastran.bdf.cards.baseCard import Element
from pyNastran.utils.mathematics import Area, norm, centroid_triangle
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)

from numpy import array

class CTRIA3s(object): # TriShell
    type = 'CTRIA3'
    asterType = 'TRIA3'
    calculixType = 'S3'

    def __init__(self):
        self._comments = {}
        self._nctria3 = 0
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

    def add_ctria3(self, card=None, data=None, comment=''):
        #TriShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            eid = integer(card, 1, 'eid')
            pid = integer(card, 2, 'pid')
            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2'),
                    integer(card, 5, 'n3')]

            thetaMcid = integer_double_or_blank(card, 6, 'thetaMcid', 0.0)
            zOffset = double_or_blank(card, 7, 'zOffset', 0.0)
            blank(card, 8, 'blank')
            blank(card, 9, 'blank')
            TFlag = integer_or_blank(card, 10, 'TFlag', 0)
            T1 = double_or_blank(card, 11, 'T1', 1.0)
            T2 = double_or_blank(card, 12, 'T2', 1.0)
            T3 = double_or_blank(card, 13, 'T3', 1.0)
            assert len(card) <= 14, 'len(CTRIA3 card) = %i' % len(card)
        else:
            eid = data[0]
            pid = data[1]
            nids = data[2:5]

            thetaMcid = data[5]
            zOffset = data[6]
            TFlag = data[7]
            T1 = data[8]
            T2 = data[9]
            T3 = data[10]
            if T1 == -1.0:
                T1 = 1.0
            if T2 == -1.0:
                T2 = 1.0
            if T3 == -1.0:
                T3 = 1.0

        #self.prepareNodeIDs(nids)
        assert len(self.nids) == 3
        self._eidmap[eid] = self._nctria3
        self.eid.append(eid)
        self.pid.append(pid)
        self.nodes.append(nids)
        self.thetaMcid.append(thetaMcid)
        self.zOffset.append(zOffset)
        self.TFlag.append(TFlag)
        self.T.append([T1, T2, T3])

        self._nctria3 += 1

    def cross_reference(self, model):
        msg = ' which is required by CTRIA3 eid=%s' % self.eid
        #self.nodes = model.Nodes(self.nodes, msg=msg)
        #self.pid = model.Property(self.pid, msg=msg)

    def Pid(self, eid):
        id = self._eidmap[eid]
        return self.pid[id]

    def _verify(self, model, xref):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.nodeIDs()

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

    #def Thickness(self):
        #if self.T1 + self.T2 + self.T3 > 0.0:
        #    if self.TFlag == 0:
        #        t = self.pid.Thickness()
        #        T1 = self.T1 / t
        #        T2 = self.T2 / t
        #        T3 = self.T3 / t
        #    else:
        #        T1 = self.T1
        #        T2 = self.T2
        #        T3 = self.T3
        #    t = (T1+T2+T3)/3.
        #else:
        #    t = self.pid.Thickness()
        #return t

    def flipNormal(self, eid):
        r"""
        ::

               1           1
              * *   -->   * *
             *   *       *   *
            2-----3     3-----2
        """
        id = self._eidmap[eid]
        (n1, n2, n3) = self.nodes[id]
        self.nodes[id] = [n1, n3, n2]

    def Interp(self, model, un):
        """
        Interpolation based on the area coordinates
        """
        (n0, n1, n2) = self.nodePositions(model)
        nc = (n0 + n1 + n2) / 3.

        a = n0 - nc
        b = n1 - nc
        c = n2 - nc

        tA1 = det(cross(b, c))   # 2*A1
        tA2 = det(cross(c, a))   # 2*A2
        tA3 = det(cross(a, b))   # 2*A3
        otA = 1. / (tA1 + tA2 + tA3)  # 1/2A

        S = array([tA1, tA2, tA3]) * otA  # Ai/A
        u = S * un
        return u

    def Jacob(self, model):
        (n0, n1, n2) = self.nodePositions(model)
        (nx0, ny0, nz0) = n0
        (nx1, ny1, nz1) = n1
        (nx2, ny2, nz2) = n2
        #J = matrix([n0,n1-n0,n2-n0])

        J = array([[nx0, nx1 - nx0, nx2 - nx0],
                   [ny0, ny1 - ny0, ny2 - ny0],
                   [nz0, nz1 - nz0, nz2 - nz0], ])
        #detJ = J.det()
        return J

    def nodeIDs(self, eid):
        id = self._eidmap[eid]
        return self.nids[id]

    #def rawFields(self):
        #list_fields = (['CTRIA3', self.eid, self.Pid()] + self.nodeIDs() +
                  #[self.thetaMcid, self.zOffset, None] + [None, self.TFlag,
                   #self.T1, self.T2, self.T3])
        #return list_fields

    def write(self, f=None, eids=None):
        msg = []
        size = 8
        for i in xrange(self._nctria3):
            eid = self.eid[i]
            pid = self.pid[i]
            zOffset = self.zOffset[i]
            TFlag = self.TFlag[i]
            thetaMcid = self.thetaMcid[i]
            T1, T2, T3 = self.T[i]

            zOffset   = set_blank_if_default(zOffset, 0.0)
            TFlag     = set_blank_if_default(TFlag, 0)
            thetaMcid = set_blank_if_default(thetaMcid, 0.0)

            T1 = set_blank_if_default(T1, 1.0)
            T2 = set_blank_if_default(T2, 1.0)
            T3 = set_blank_if_default(T3, 1.0)

            fields = (['CTRIA3', eid, pid] + self.nodeIDs(eid) +
                      [thetaMcid, zOffset, None] + [None, TFlag, T1, T2, T3])
            msg.append(print_card(fields))
        return ''.join(msg)
