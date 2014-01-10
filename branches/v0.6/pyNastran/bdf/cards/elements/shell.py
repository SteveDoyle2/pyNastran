# pylint: disable=C0103,R0902,R0904,R0914,C0302
"""
All shell elements are defined in this file.  This includes:

 * CTRIA3
 * CTRIA6
 * CTRIAX
 * CTRIAX6
 * CSHEAR
 * CQUAD
 * CQUAD4
 * CQUAD8
 * CQUADR
 * CQUADX

All tris are TriShell, ShellElement, and Element objects.
All quads are QuadShell, ShellElement, and Element objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from numpy import array, eye, cross, allclose, dot, transpose, zeros
from numpy.linalg import det, norm  # inv

from pyNastran.bdf.fieldWriter import (set_blank_if_default,
                                       set_default_if_blank)
from pyNastran.bdf.cards.baseCard import Element
from pyNastran.utils import list_print
from pyNastran.utils.mathematics import Area, norm, centroid_triangle
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)
from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16

def _triangle_area_centroid_normal(nodes):
    """
    Returns area,centroid,unitNormal

    :param nodes: list of three triangle vertices

    ::

      n = Normal = a x b
      Area = 1/2 * |a x b|
      V = <v1,v2,v3>
      |V| = sqrt(v1^0.5+v2^0.5+v3^0.5) = norm(V)

      Area = 0.5 * |n|
      unitNormal = n/|n|
    """
    (n0, n1, n2) = nodes
    vector = cross(n0 - n1, n0 - n2)
    length = norm(vector)
    normal = vector / length
    if not allclose(norm(normal), 1.):
        msg = ('function _triangle_area_centroid_normal, check...\n'
               'a = {0}\nb = {1}\nnormal = {2}\nlength = {3}\n'.format(
               n0 - n1, n0 - n2, normal, length))
        raise RuntimeError(msg)
    return (0.5 * length, (n0 + n1 + n2) / 3., normal)


def _normal(a, b):
    """Finds the unit normal vector of 2 vectors"""
    vector = cross(a, b)
    normal = vector / norm(vector)
    assert allclose(norm(normal), 1.)
    return normal


class ShellElement(Element):
    type = 'ShellElement'

    def __init__(self, card, data):
        Element.__init__(self, card, data)

    def Eid(self):
        return self.eid

    def Area(self):
        raise NotImplementedError('Area undefined for %s' % self.type)

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def mid(self):
        """
        Returns the material
        """
        return self.pid.mid()

    def Mid(self):
        """
        Returns the material ID
        """
        return self.pid.Mid()

    def Rho(self):
        """
        Returns the density
        """
        return self.pid.mid().rho

    def Nsm(self):
        """
        Returns the non-structural mass
        """
        return self.pid.Nsm()

    def MassPerArea(self):
        """
        Returns the mass per area
        """
        return self.pid.MassPerArea()

    def Mass(self):
        r"""
        .. math:: m = \frac{m}{A} A  \f]
        """
        return self.pid.MassPerArea() * self.Area()

    def flipNormal(self):
        raise NotImplementedError('flipNormal undefined for %s' % (self.type))


class TriShell(ShellElement):
    def __init__(self, card, data):
        ShellElement.__init__(self, card, data)

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def AreaCentroidNormal(self):
        """
        Returns area,centroid, normal as it's more efficient to do them
        together
        """
        (n0, n1, n2) = self.nodePositions()
        return _triangle_area_centroid_normal([n0, n1, n2])

    def Area(self):
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} (n_0-n_1) \times (n_0-n_2)
        """
        (n0, n1, n2) = self.nodePositions()
        a = n0 - n1
        b = n0 - n2
        area = Area(a, b)
        return area

    def Normal(self):
        r"""
        Get the normal vector, :math:`n`.

        .. math::
          n = \frac{(n_0-n_1) \times (n_0-n_2)}
             {\lvert (n_0-n_1) \times (n_0-n_2) \lvert}
        """
        (n0, n1, n2) = self.nodePositions()
        return _normal(n0 - n1, n0 - n2)

    def Centroid(self):
        r"""
        Get the centroid.

        .. math::
          CG = \frac{1}{3} (n_0+n_1+n_2)
        """
        (n0, n1, n2) = self.nodePositions()
        centroid = centroid_triangle(n0, n1, n2)
        return centroid

    def MassMatrix(self, isLumped=True):
        """
        6x6 mass matrix triangle
        http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/IFEM.Ch32.d/IFEM.Ch32.pdf
        """
        mass = self.Mass()  # rho*A*t
        if isLumped:  # lumped mass matrix
            Mass = mass / 3 * eye(6)
        else:  # consistent mass
            M = eye(6) * 2.
            M[2, 0] = M[4, 0] = M[0, 2] = M[0, 4] = 1.
            M[3, 1] = M[5, 1] = M[1, 3] = M[1, 5] = 1.
            M[4, 2] = M[2, 4] = 1.
            M[5, 3] = M[5, 3] = 1.
            Mass = mass / 12 * M
        return Mass


class CTRIA3(TriShell):
    type = 'CTRIA3'
    asterType = 'TRIA3'
    calculixType = 'S3'
    _field_map = {1: 'eid', 2:'pid', 6:'thetaMcid', 7:'zOffset', 10:'TFlag', 11:'T1', 12:'T2', 13:'T3'}

    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 4:
            self.nodes[1] = value
        elif n == 5:
            self.nodes[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        TriShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element ID
            self.eid = integer(card, 1, 'eid')
            #: Property ID
            self.pid = integer(card, 2, 'pid')

            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2'),
                    integer(card, 5, 'n3')]

            self.thetaMcid = integer_double_or_blank(card, 6, 'thetaMcid', 0.0)
            self.zOffset = double_or_blank(card, 7, 'zOffset', 0.0)
            blank(card, 8, 'blank')
            blank(card, 9, 'blank')

            self.TFlag = integer_or_blank(card, 10, 'TFlag', 0)
            self.T1 = double_or_blank(card, 11, 'T1', 1.0)
            self.T2 = double_or_blank(card, 12, 'T2', 1.0)
            self.T3 = double_or_blank(card, 13, 'T3', 1.0)
            assert len(card) <= 14, 'len(CTRIA3 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:5]

            self.thetaMcid = data[5]
            self.zOffset = data[6]
            self.TFlag = data[7]
            self.T1 = data[8]
            self.T2 = data[9]
            self.T3 = data[10]
            if self.T1 == -1.0:
                self.T1 = 1.0
            if self.T2 == -1.0:
                self.T2 = 1.0
            if self.T3 == -1.0:
                self.T3 = 1.0

        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 3

    def displacement_stress(self, model, q, dofs):
        n1, n2, n3 = self.nodes
        p1 = n1.Position()
        p2 = n2.Position()
        p3 = n3.Position()
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        x3, y3, z3 = p3

        v1 = p2 - p1
        v13 = p3 - p1

        v3 = cross(v1, v13)
        A = 0.5 * norm(v3)

        # constant strain triangle - PSHELL
        a1 = x2*y3 - x3*y2
        a2 = x3*y1 - x1*y3
        a3 = x1*y2 - x2*y1

        b1 = y2 - y3
        b2 = y3 - y1
        b3 = y1 - y2

        g1 = -x2 + x3
        g2 = -x3 + x1
        g3 = -x1 + x2
        B = array([
            [b1,  0, b2,  0, b3,  0],
            [0,  g1,  0, g2,  0, g3],
            [g1, b1, g2, b2, g3, b3],
        ])

        mat = self.pid.mid1
        nu = mat.Nu()
        E = mat.E() / (1-nu**2)
        G = mat.G()

        Ep = array([
            [E,      nu * E, 0.],
            [nu * E, E,      0.],
            [0.,     0.,      G],
        ])
        #===========
        # new stuff

        n1, n2, n3 = self.nodeIDs()
        n11 = dofs[(n1, 1)]
        n21 = dofs[(n2, 1)]
        n31 = dofs[(n3, 1)]

        n12 = dofs[(n1, 2)]
        n22 = dofs[(n2, 2)]
        n32 = dofs[(n3, 2)]

        q2 = array([
            q[n11], q[n12],  # 0, 1
            q[n21], q[n22],  # 2, 3
            q[n31], q[n32],  # 4, 5
        ])
        dx = q[2] - q[0]
        dy = q[4] - q[0]

        strain = dot(B, q2)
        stress = dot(Ep, dot(B, q2))
        F = stress * A

        return stress, strain, F

    def Stiffness(self, model, node_ids, index0s, gravLoad, is3D, fnorm):
        assert is3D == True

        print("node_ids =", node_ids)
        n1, n2, n3 = self.nodes
        p1 = n1.Position()
        p2 = n2.Position()
        p3 = n3.Position()

        v1 = p2 - p1
        v13 = p3 - p1

        v3 = cross(v1, v13)
        A = 0.5 * norm(v3)
        #A = 0.5 * norm(cross(v1, v13))

        mat = self.pid.mid1

        nu = mat.Nu()
        E = mat.E() / (1-nu**2)
        G = mat.G()

        x1, y1, z1 = p1
        x2, y2, z2 = p2
        x3, y3, z3 = p3

        # constant strain triangle - PSHELL
        a1 = x2*y3 - x3*y2
        a2 = x3*y1 - x1*y3
        a3 = x1*y2 - x2*y1

        b1 = y2 - y3
        b2 = y3 - y1
        b3 = y1 - y2

        g1 = -x2 + x3
        g2 = -x3 + x1
        g3 = -x1 + x2
        B = array([
            [b1,  0, b2,  0, b3,  0],
            [0,  g1,  0, g2,  0, g3],
            [g1, b1, g2, b2, g3, b3],
        ])
        Ep = array([
            [E,      nu * E, 0.],
            [nu * E, E,      0.],
            [0.,     0.,      G],
        ])

        t = self.Thickness()
        ki = t / (4 * A)
        K = ki * dot(transpose(B), dot(Ep, B))

        n1, n2, n3 = node_ids
        dofs = [
            (n1, 1), (n1, 2),
            (n2, 1), (n2, 2),
            (n3, 1), (n3, 2),
        ]
        nIJV = dofs
        #nIJV = []

        v1 /= norm(v1)
        v3 /= norm(v3)

        v2 = cross(v1, v3)
        R = array([v1, v2, v3])
        assert R.shape == (3,3), 'R=%s' % R

        Fg = []
        nGrav = []
        return (K, dofs, nIJV, Fg, nGrav)

    def cross_reference(self, model):
        msg = ' which is required by CTRIA3 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.nodeIDs()

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i,nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            if not self.pid.type in ['PLPLANE']:
                t = self.Thickness()
                assert isinstance(t, float), 'thickness=%r' % t
                mass = self.Mass()
                assert isinstance(mass, float), 'mass=%r' % mass
            a,c,n = self.AreaCentroidNormal()
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                assert isinstance(n[i], float)

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

    def flipNormal(self):
        """
        ::

               1           1
              * *   -->   * *
             *   *       *   *
            2-----3     3-----2
        """
        (n1, n2, n3) = self.nodes
        self.nodes = [n1, n3, n2]

    def Interp(self, un):
        """
        Interpolation based on the area coordinates
        """
        (n0, n1, n2) = self.nodePositions()
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

    def Jacob(self):
        (n0, n1, n2) = self.nodePositions()
        (nx0, ny0, nz0) = n0
        (nx1, ny1, nz1) = n1
        (nx2, ny2, nz2) = n2
        #J = matrix([n0,n1-n0,n2-n0])

        J = array([[nx0, nx1 - nx0, nx2 - nx0],
                   [ny0, ny1 - ny0, ny2 - ny0],
                   [nz0, nz1 - nz0, nz2 - nz0], ])
        #detJ = J.det()
        return J

    def getReprDefaults(self):
        zOffset = set_blank_if_default(self.zOffset, 0.0)
        TFlag = set_blank_if_default(self.TFlag, 0)
        thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        return (thetaMcid, zOffset, TFlag, T1, T2, T3)

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=False)

    def rawFields(self):
        list_fields = (['CTRIA3', self.eid, self.Pid()] + self.nodeIDs() +
                  [self.thetaMcid, self.zOffset, None] + [None, self.TFlag,
                   self.T1, self.T2, self.T3])
        return list_fields

    def reprFields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3) = self.getReprDefaults()
        list_fields = ([self.type, self.eid, self.Pid()] + self.nodeIDs() +
                  [thetaMcid, zOffset, None] + [None, TFlag, T1, T2, T3])
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        if size == 8 or len(card) == 6: # to last node
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)

class CTRIA6(TriShell):
    type = 'CTRIA6'
    asterType = 'TRIA6'
    calculixType = 'S6'

    def __init__(self, card=None, data=None, comment=''):
        TriShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element ID
            self.eid = integer(card, 1, 'eid')
            #: Property ID
            self.pid = integer(card, 2, 'pid')

            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2'),
                    integer(card, 5, 'n3'),
                    integer_or_blank(card, 6, 'n4', 0),
                    integer_or_blank(card, 7, 'n5', 0),
                    integer_or_blank(card, 8, 'n6', 0)]

            self.thetaMcid = integer_double_or_blank(card, 9, 'thetaMcid', 0.0)
            self.zOffset = double_or_blank(card, 10, 'zOffset', 0.0)

            self.T1 = double_or_blank(card, 11, 'T1', 1.0)
            self.T2 = double_or_blank(card, 12, 'T2', 1.0)
            self.T3 = double_or_blank(card, 13, 'T3', 1.0)
            self.TFlag = integer_or_blank(card, 14, 'TFlag', 0)
            assert len(card) <= 15, 'len(CTRIA6 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:8]
            self.thetaMcid = data[8]
            self.zOffset = data[8]
            self.T1 = data[9]
            self.T2 = data[10]
            self.T3 = data[11]
            self.TFlag = data[12]
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(nids) == 6, 'error on CTRIA6'

        #print self.thetaMcid
        #print card

        #print "self.xi = ",self.xi
        #raise

    def cross_reference(self, model):
        msg = ' which is required by CTRIA6 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.nodeIDs()

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            if not self.pid.type in ['PLPLANE']:
                t = self.Thickness()
                assert isinstance(t, float), 'thickness=%r' % t
                mass = self.Mass()
                assert isinstance(mass, float), 'mass=%r' % mass
            a,c,n = self.AreaCentroidNormal()
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                assert isinstance(n[i], float)

    def Thickness(self):
        """
        Returns the thickness, :math:`t`
        """
        return self.pid.Thickness()

    def AreaCentroidNormal(self):
        """
        Returns area, centroid, normal as it's more efficient to do them
        together
        """
        (n1, n2, n3, n4, n5, n6) = self.nodePositions()
        return _triangle_area_centroid_normal([n1, n2, n3])

    def Area(self):
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} (n_0-n_1) \times (n_0-n_2)
        """
        (n1, n2, n3, n4, n5, n6) = self.nodePositions()
        a = n1 - n2
        b = n1 - n3
        area = Area(a, b)
        return area

    def Normal(self):
        r"""
        Get the normal vector, :math:`n`.

        .. math::
          n = \frac{(n_0-n_1) \times (n_0-n_2)}{\lvert (n_0-n_1) \times (n_0-n_2) \lvert}
        """
        (n0, n1, n2) = self.nodePositions()[:3]
        return _normal(n0 - n1, n0 - n2)

    def Centroid(self):
        r"""
        Get the centroid.

        .. math::
          CG = \frac{1}{3} (n_1+n_2+n_3)
        """
        (n1, n2, n3, n4, n5, n6) = self.nodePositions()
        centroid = centroid_triangle(n1, n2, n3)
        return centroid

    def flipNormal(self):
        r"""
        ::

               1                1
               **               **
              *  *             *  *
             4    6   -->     6    4
            *      *         *      *
           2----5---3       3----5---2
        """
        (n1, n2, n3, n4, n5, n6) = self.nodes
        self.nodes = [n1, n3, n2, n6, n5, n4]

    def getReprDefaults(self):
        zOffset = set_blank_if_default(self.zOffset, 0.0)
        TFlag = set_blank_if_default(self.TFlag, 0)
        thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        return (thetaMcid, zOffset, TFlag, T1, T2, T3)

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=True)

    def rawFields(self):
        list_fields = (['CTRIA6', self.eid, self.Pid()] + self.nodeIDs() +
                  [self.thetaMcid, self.zOffset, None] + [None, self.TFlag,
                  self.T1, self.T2, self.T3])
        return list_fields

    def reprFields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3) = self.getReprDefaults()
        list_fields = (['CTRIA6', self.eid, self.Pid()] + self.nodeIDs() +
                  [thetaMcid, zOffset, None] + [None, TFlag, T1, T2, T3])
        return list_fields


class CTRIAR(TriShell):
    type = 'CTRIAR'
    def __init__(self, card=None, data=None, comment=''):
        TriShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = integer(card, 1, 'eid')
        #: Property ID
        self.pid = integer(card, 2, 'pid')

        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3')]

        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 3

        self.thetaMcid = integer_double_or_blank(card, 6, 'thetaMcid', 0.0)
        self.zOffset = double_or_blank(card, 7, 'zOffset', 0.0)
        blank(card, 8, 'blank')
        blank(card, 9, 'blank')
        blank(card, 10, 'blank')

        self.TFlag = integer_or_blank(card, 11, 'TFlag', 0)
        self.T1 = double_or_blank(card, 11, 'T1', 1.0)
        self.T2 = double_or_blank(card, 12, 'T2', 1.0)
        self.T3 = double_or_blank(card, 13, 'T3', 1.0)
        assert len(card) <= 14, 'len(CTRIAR card) = %i' % len(card)

    def cross_reference(self, model):
        msg = ' which is required by CTRIAR eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.nodeIDs()

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i,nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PSHELL', 'PCOMP'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            t = self.Thickness()
            a,c,n = self.AreaCentroidNormal()
            assert isinstance(t, float), 'thickness=%r' % t
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                assert isinstance(n[i], float)
            mass = self.Mass()
            assert isinstance(mass, float), 'mass=%r' % mass

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def flipNormal(self):
        r"""
        ::

               1           1
              * *   -->   * *
             *   *       *   *
            2-----3     3-----2
        """
        (n1, n2, n3) = self.nodes
        self.nodes = [n1, n3, n2]

    def getReprDefaults(self):
        zOffset = set_blank_if_default(self.zOffset, 0.0)
        TFlag = set_blank_if_default(self.TFlag, 0)
        thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        return (thetaMcid, zOffset, TFlag, T1, T2, T3)

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.nodeIDs()

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            # PSHELL/PCOMP
            assert self.pid.type in ['PSHELL'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            t = self.Thickness()
            a,c,n = self.AreaCentroidNormal()
            assert isinstance(t, float), 'thickness=%r' % t
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                #assert isinstance(n[i], float)
            mass = self.Mass()
            assert isinstance(mass, float), 'mass=%r' % mass

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=False)

    def rawFields(self):
        list_fields = [self.eid, self.Pid()] + self.nodeIDs() + [self.thetaMcid,
                  self.zOffset, self.TFlag, self.T1, self.T2, self.T3]
        return list_fields

    def reprFields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3) = self.getReprDefaults()
        list_fields = (['CTRIAR', self.eid, self.Pid()] + self.nodeIDs() +
                  [thetaMcid, zOffset, None, None, TFlag, T1, T2, T3])
        return list_fields


class CTRIAX(TriShell):
    type = 'CTRIAX'
    calculixType = 'CAX6'
    def __init__(self, card=None, data=None, comment=''):
        TriShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element ID
            self.eid = integer(card, 1, 'eid')
            #: Property ID
            self.pid = integer(card, 2, 'pid')

            nids = [integer_or_blank(card, 3, 'n1'),
                    integer_or_blank(card, 4, 'n2'),
                    integer_or_blank(card, 5, 'n3'),
                    integer_or_blank(card, 6, 'n4'),
                    integer_or_blank(card, 7, 'n5'),
                    integer_or_blank(card, 8, 'n6')]
            self.thetaMcid = integer_double_or_blank(card, 9, 'theta_mcsid', 0.0)
            assert len(card) <= 10, 'len(CTRIAX card) = %i' % len(card)
        else:
            raise NotImplementedError(data)
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(nids) == 6, 'error on CTRIAX'

    def cross_reference(self, model):
        msg = ' which is required by CTRIAX eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=True)

    def rawFields(self):
        nodeIDs = self.nodeIDs()
        list_fields = ['CTRIAX', self.eid, self.Pid()] + nodeIDs + [self.thetaMcid]
        return list_fields

    def reprFields(self):
        thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)
        nodeIDs = self.nodeIDs()
        list_fields = ['CTRIAX', self.eid, self.Pid()] + nodeIDs + [thetaMcid]
        return list_fields


class CTRIAX6(TriShell):
    """
    Nodes defined in a non-standard way
    ::
           5
          / \
         6   4
       /       \
      1----2----3
    """
    type = 'CTRIAX6'
    #calculixType = 'CAX6'
    def __init__(self, card=None, data=None, comment=''):
        TriShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element ID
            self.eid = integer(card, 1, 'eid')
            #: Material ID
            self.mid = integer(card, 2, 'mid')

            nids = [integer(card, 3, 'n1'),
                    integer_or_blank(card, 4, 'n2'),
                    integer(card, 5, 'n3'),
                    integer_or_blank(card, 6, 'n4'),
                    integer(card, 7, 'n5'),
                    integer_or_blank(card, 8, 'n6')]

            #: theta
            self.theta = double_or_blank(card, 9, 'theta', 0.0)
            assert len(card) <= 10, 'len(CTRIAX6 card) = %i' % len(card)
        else:
            raise NotImplementedError(data)
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(nids) == 6, 'error on CTRIAX6'

    def cross_reference(self, model):
        msg = ' which is required by CTRIAX6 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.mid = model.Material(self.mid)

    def _verify(self, xref=True):
        eid = self.Eid()
        #pid = self.Pid()
        nids = self.nodeIDs()
        assert self.pid == 0, 'pid = %s' % self.pid
        assert isinstance(eid, int)
        #assert isinstance(pid, int)
        for i,nid in enumerate(nids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer or blank; nid=%s' %(i, nid)

        if xref:
            assert self.mid.type in ['MAT1'], 'mid=%s self.mid.type=%s' % (mid, self.mid.type)
            #assert self.pid.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            #if not self.pid.type in ['PLPLANE']:
                #t = self.Thickness()
                #assert isinstance(t, float), 'thickness=%r' % t
                #mass = self.Mass()
                #assert isinstance(mass, float), 'mass=%r' % mass
            a,c,n = self.AreaCentroidNormal()
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                assert isinstance(n[i], float)

    def Pid(self):
        raise AttributeError("CTRIAX6 doesn't have a Property")

    def AreaCentroidNormal(self):
        """
        Returns area,centroid, normal as it's more efficient to do them
        together
        """
        (n0, n1, n2, n3, n4, n5) = self.nodePositions()
        return _triangle_area_centroid_normal([n0, n2, n4])

    def Area(self):
        r"""
        Get the normal vector.

        .. math:: A = \frac{1}{2} (n_0-n_1) times (n_0-n_2)
        """
        (n1, n2, n3, n4, n5, n6) = self.nodePositions()
        a = n1 - n3
        b = n1 - n5
        area = Area(a, b)
        return area

    def Thickness(self):
        """
        CTRIAX doesn't have a thickness because ???
        """
        raise AttributeError('CTRIAX6 does not have a thickness')

    def Nsm(self):
        raise AttributeError('CTRIAX6 does not have a non-structural mass')

    def MassPerArea(self):
        raise AttributeError('CTRIAX6 does not have a MassPerArea')

    def Mass(self):
        raise NotImplementedError('CTRIAX6 does not have a Mass method yet')

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        return self.mid.mid

    def flipNormal(self):
        r"""
        ::

               5               5
              / \             / \
             6   4   -->     6   4
           /       \       /       \
          1----2----3     1----2----3
        """
        (n1, n2, n3, n4, n5, n6) = self.nodes
        self.nodes = [n1, n6, n5, n4, n3, n2]

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=True)

    def rawFields(self):
        list_fields = (['CTRIAX6', self.eid, self.Mid(), self.Pid()] +
                  self.nodeIDs() +  [self.theta])
        return list_fields

    def reprFields(self):
        theta = set_default_if_blank(self.theta, 0.0)
        list_fields = ['CTRIAX6', self.eid, self.Mid()] + self.nodeIDs() + [theta]
        return list_fields


class QuadShell(ShellElement):
    def __init__(self, card, data):
        ShellElement.__init__(self, card, data)

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def Normal(self):
        (n1, n2, n3, n4) = self.nodePositions()
        return _normal(n1 - n3, n2 - n4)

    def AreaCentroidNormal(self):
        (area, centroid) = self.AreaCentroid()
        normal = self.Normal()
        return (area, centroid, normal)

    def AreaCentroid(self):
        """
        ::
          1-----2
          |    /|
          | A1/ |
          |  /  |
          |/ A2 |
          4-----3

        .. math:
            c = \frac{\sum(c_i A_i){\sum{A_i}}

         c = sum(ci*Ai)/sum(A)
         where:
           c=centroid
           A=area
        """
        (n1, n2, n3, n4) = self.nodePositions()
        a = n1 - n2
        b = n2 - n4
        area1 = Area(a, b)
        c1 = centroid_triangle(n1, n2, n4)

        a = n2 - n4
        b = n2 - n3
        area2 = Area(a, b)
        c2 = centroid_triangle(n2, n3, n4)

        area = area1 + area2
        centroid = (c1 * area1 + c2 * area2) / area
        return(area, centroid)

    def Centroid(self):
        #nodes = self.nodePositions()
        (area, centroid) = self.AreaCentroid()
        return centroid

    def Area(self):
        r"""
        \f[ A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert \f]
        where a and b are the quad's cross node point vectors
        """
        (n1, n2, n3, n4) = self.nodePositions()
        a = n1 - n3
        b = n2 - n4
        area = Area(a, b)
        return area

    def MassMatrix(self, isLumped=True, gauss=1):
        """
        6x6 mass matrix triangle
        http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/IFEM.Ch32.d/IFEM.Ch32.pdf
        """
        mass = self.Mass()  # rho*A*t
        if isLumped:  # lumped mass matrix
            Mass = mass / 3 * eye(6)
        else:  # consistent mass
            if gauss == 1:
                M = eye(8) * 1.  # 1x1 Gauss Rule
                M[2, 0] = M[3, 1] = M[4, 2] = M[5, 3] = M[6, 4] = M[7, 5] = 1.
                M[4, 0] = M[5, 1] = M[6, 2] = M[7, 3] = 1.
                M[6, 0] = M[7, 1] = 1.

                M[0, 2] = M[1, 3] = M[2, 4] = M[3, 5] = M[4, 6] = M[5, 7] = 1.
                M[0, 4] = M[1, 5] = M[2, 6] = M[3, 7] = 1.
                M[0, 6] = M[1, 7] = 1.
                Mass = mass / 32 * M
            if gauss == 2:
                M = eye(8) * 4.  # 2x2 Gauss Rule
                M[2, 0] = M[3, 1] = M[4, 2] = M[5, 3] = M[6, 4] = M[7, 5] = 2.
                M[4, 0] = M[5, 1] = M[6, 2] = M[7, 3] = 1.
                M[6, 0] = M[7, 1] = 2.

                M[0, 2] = M[1, 3] = M[2, 4] = M[3, 5] = M[4, 6] = M[5, 7] = 2.
                M[0, 4] = M[1, 5] = M[2, 6] = M[3, 7] = 1.
                M[0, 6] = M[1, 7] = 2.
                Mass = mass / 72 * M
        return Mass

    def flipNormal(self):
        r"""
        ::

          1---2       1---4
          |   |  -->  |   |
          |   |       |   |
          4---3       2---3
        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]

    def getReprDefaults(self):
        zOffset = set_blank_if_default(self.zOffset, 0.0)
        TFlag = set_blank_if_default(self.TFlag, 0)
        thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        T4 = set_blank_if_default(self.T4, 1.0)

        if 0:
            print("eid       = %s" % self.eid)
            print("nodes     = %s" % self.nodes)

            print("self.zOffset   = %s" % self.zOffset)
            print("self.TFlag     = %s" % self.TFlag)
            print("self.thetaMcid = %s" % self.thetaMcid)

            print("zOffset   = %s" % zOffset)
            print("TFlag     = %s" % TFlag)
            print("thetaMcid = %s" % thetaMcid)

            print("T1 = %s" % T1)
            print("T2 = %s" % T2)
            print("T3 = %s" % T3)
            print("T4 = %s\n" % T4)
        return (thetaMcid, zOffset, TFlag, T1, T2, T3, T4)


class CSHEAR(QuadShell):
    type = 'CSHEAR'
    calculixType = 'S4'
    def __init__(self, card=None, data=None, comment=''):
        QuadShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element ID
            self.eid = integer(card, 1, 'eid')
            #: Property ID
            self.pid = integer(card, 2, 'pid')
            nids = [integer_or_blank(card, 3, 'n1'),
                    integer_or_blank(card, 4, 'n2'),
                    integer_or_blank(card, 5, 'n3'),
                    integer_or_blank(card, 6, 'n4')]
            assert len(card) <= 7, 'len(CSHEAR card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 4

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=True)

    def cross_reference(self, model):
        msg = ' which is required by CSHEAR eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def Normal(self):
        (n1, n2, n3, n4) = self.nodePositions()
        return _normal(n1 - n3, n2 - n4)

    def AreaCentroidNormal(self):
        (area, centroid) = self.AreaCentroid()
        normal = self.Normal()
        return (area, centroid, normal)

    def AreaCentroid(self):
        """
        ::
          1-----2
          |    /|
          | A1/ |
          |  /  |
          |/ A2 |
          4-----3

          centroid
             c = sum(ci*Ai)/sum(A)
             where:
               c=centroid
               A=area
        """
        (n1, n2, n3, n4) = self.nodePositions()
        a = n1 - n2
        b = n2 - n4
        area1 = Area(a, b)
        c1 = centroid_triangle(n1, n2, n4)

        a = n2 - n4
        b = n2 - n3
        area2 = Area(a, b)
        c2 = centroid_triangle(n2, n3, n4)

        area = area1 + area2
        centroid = (c1 * area1 + c2 * area2) / area
        return(area, centroid)

    def Centroid(self):
        (area, centroid) = self.AreaCentroid()
        return centroid

    def Area(self):
        r"""
        \f[ A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert \f]
        where a and b are the quad's cross node point vectors
        """
        (n1, n2, n3, n4) = self.nodePositions()
        a = n1 - n3
        b = n2 - n4
        area = Area(a, b)
        return area

    def flipNormal(self):
        r"""
        ::

          1---2       1---4
          |   |  -->  |   |
          |   |       |   |
          4---3       2---3
        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=False)

    def rawFields(self):
        list_fields = ['CSHEAR', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def G(self):
        return self.pid.mid.G()

    def Thickness(self):
        return self.pid.t

    def Stiffness(self, model, node_ids, index0s, fnorm):
        n1, n2, n3, n4 = self.nodeIDs()

        (p1, p2, p3, p4) = self.nodePositions()
        v1 = p2 - p1
        v1b = p3 - p4

        v2 = p4 - p1
        v2b = p3 - p2
        # average of lengths & widths
        a = (norm(v1) + norm(v1b)) / 2.
        b = (norm(v2) + norm(v2b)) / 2.
        v1 /= a
        v2 /= b
        assert a > 0
        assert b > 0
        Lambda = array([v1, v2])
        assert Lambda.shape == (2, 3), 'shape=%s' % (Lambda.shape)
        Lambda2 = LambdaN(Lambda, 4) # n=4 points

        G = self.G()
        t = self.Thickness()

        ki = G * t / 4.0
        K = ki * array([
            [a/b, 1, a/b, -1, -a/b, -1, -a/b, 1],
            [1, b/a, 1, -b/a, -1, -b/a, -1, b/a],
            [a/b, 1, a/b, -1, -a/b, -1, -a/b, 1],  # row 0
            [-1, -b/a, -1, b/a, 1, b/a, 1, -b/a],  # -row 1
            [-a/b, -1, -a/b, 1, a/b, 1, a/b, -1],  # -row 0
            [-1, -b/a, -1, b/a, 1, b/a, 1, -b/a],  # -row 1
            [-a/b, -1, -a/b, 1, a/b, 1, a/b, -1],  # -row 0
            [1, b/a, 1, -b/a, -1, -b/a, -1, b/a],  # row 1
        ])
        K2 = dot(transpose(Lambda2), dot(K, Lambda2))

        dofs = []
        nIJV = [
            (n1, 1), (n1, 2), (n1, 3),
            (n2, 1), (n2, 2), (n2, 3),
            (n3, 1), (n3, 2), (n3, 3),
            (n4, 1), (n4, 2), (n4, 3),
        ]
        return (K2, dofs, nIJV)

    def displacement_stress(self, model, q, dofs):
        n1, n2, n3, n4 = self.nodeIDs()
        n = array([
           dofs[(n1, 1)],
           dofs[(n1, 2)],
           dofs[(n1, 3)],

           dofs[(n2, 1)],
           dofs[(n2, 2)],
           dofs[(n2, 3)],

           dofs[(n3, 1)],
           dofs[(n3, 2)],
           dofs[(n3, 3)],

           dofs[(n4, 1)],
           dofs[(n4, 2)],
           dofs[(n4, 3)]
        ])
        u = q[n]

        #=========
        (p1, p2, p3, p4) = self.nodePositions()
        v1 = p2 - p1
        v1b = p3 - p4

        v2 = p4 - p1
        v2b = p3 - p2
        # average of lengths & widths
        a = (norm(v1) + norm(v1b)) / 2.
        b = (norm(v2) + norm(v2b)) / 2.
        v1 /= a
        v2 /= b
        assert a > 0
        assert b > 0
        Lambda = array([v1, v2])
        assert Lambda.shape == (2, 3), 'shape=%s' % (Lambda.shape)
        Lambda2 = LambdaN(Lambda, 4) # n=4 points
        #=========
        G = self.G()
        t = self.Thickness()

        u2 = dot(Lambda2, u)

        ki = G * t / 4.
        K = ki * array([
            [a/b, 1, a/b, -1, -a/b, -1, -a/b, 1],
            [1, b/a, 1, -b/a, -1, -b/a, -1, b/a],
            [a/b, 1, a/b, -1, -a/b, -1, -a/b, 1],  # row 0
            [-1, -b/a, -1, b/a, 1, b/a, 1, -b/a],  # -row 1
            [-a/b, -1, -a/b, 1, a/b, 1, a/b, -1],  # -row 0
            [-1, -b/a, -1, b/a, 1, b/a, 1, -b/a],  # -row 1
            [-a/b, -1, -a/b, 1, a/b, 1, a/b, -1],  # -row 0
            [1, b/a, 1, -b/a, -1, -b/a, -1, b/a],  # row 1
        ])
        K2 = dot(transpose(Lambda2), dot(K, Lambda2))
        F = dot(K2, u)

        kit = G / (a * b * t)
        #Kt = kit * array([-a, -b, -a, b,  a, b,  a, -b])
        #tau_xy = Kt * u2

        Fx = F[ [0, 3, 6, 9] ]
        Fy = F[ [1, 4, 7, 10] ]

        # ^ b
        # |
        # 4----3
        # |    |
        # |    |
        # 1----2 ---> a
        F1x = Fx[0] + Fx[1]  # 1-2
        F2y = Fy[1] + Fy[2]  # 2-3
        F3x = Fx[2] + Fx[3]  # 3-4
        F4y = Fy[0] + Fy[3]  # 4-1

        tau41 = F1x / a
        tau12 = F2y / b
        tau23 = F3x / a
        tau34 = F4y / b

        tau14 = tau41

        ss = self.pid.mid.Ss  # shear stress
        exy = tau41 / G
        if ss:
            margin = tau41 / ss
        else:
            margin = 0.0
        shear12 = 0.0
        shear23 = 0.0

#        ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======
#        ELEMENT        F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1
#              ID       KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41
        #f14 = Fy[3] - Fy[0] # F4 - F1
        #f12 = Fx[1] - Fx[0] # F2 - F1

        #f21 = -f12
        #f23 = Fy[2] - Fy[1]  # F3 - F2

        #f32 = -f23
        #f34 = Fx[3] - Fx[2] # F4 - F3

        #f41 = -f14
        #f43 = Fx[2] - Fx[3] # F3 - F4

        #==========================
        f14 = Fy[0]
        f12 = Fx[0]

        f21 = Fx[1]
        f23 = Fy[2]

        f32 = Fy[2]
        f34 = Fx[2]

        f41 = Fy[3]
        f43 = Fx[3]

        kick1 = 0.0
        kick2 = 0.0
        kick3 = 0.0
        kick4 = 0.0

        F = [
            f41, f21, tau12, kick1,
            f12, f32, tau23, kick2,
            f23, f43, tau34, kick3,
            f14, f34, tau41, kick4,
             ]
        return ([tau14, tau14, margin], [exy, exy, margin/G], F)





def LambdaN(L, n):
    nx, ny = L.shape
    N = zeros((nx*n, ny*n), 'float64')
    for i in xrange(n):
        ix0 = i * nx
        iy0 = i * ny
        ix = (i + 1) * nx
        iy = (i + 1) * ny
        N[ix0:ix, iy0:iy] = L
    return N

class CQUAD4(QuadShell):
    type = 'CQUAD4'
    asterType = 'QUAD4 # CQUAD4'
    calculixType = 'S4'
    _field_map = {1: 'eid', 2:'pid', 7:'thetaMcid', 8:'zOffset', 10:'TFlag', 11:'T1', 12:'T2', 13:'T3'}

    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 4:
            self.nodes[1] = value
        elif n == 5:
            self.nodes[2] = value
        elif n == 6:
            self.nodes[3] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        QuadShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element ID
            self.eid = integer(card, 1, 'eid')
            #: Property ID
            self.pid = integer(card, 2, 'pid')
            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2'),
                    integer(card, 5, 'n3'),
                    integer(card, 6, 'n4')]
            self.thetaMcid = integer_double_or_blank(card, 7, 'thetaMcid', 0.0)
            self.zOffset = double_or_blank(card, 8, 'zOffset', 0.0)
            blank(card, 9, 'blank')
            self.TFlag = integer_or_blank(card, 10, 'TFlag', 0)
            self.T1 = double_or_blank(card, 11, 'T1', 1.0)
            self.T2 = double_or_blank(card, 12, 'T2', 1.0)
            self.T3 = double_or_blank(card, 13, 'T3', 1.0)
            self.T4 = double_or_blank(card, 14, 'T4', 1.0)
            assert len(card) <= 15, 'len(CQUAD4 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:6]

            self.thetaMcid = data[6]
            self.zOffset = data[7]
            self.TFlag = data[8]
            self.T1 = data[9]
            self.T2 = data[10]
            self.T3 = data[11]
            self.T4 = data[12]
            if self.T1 == -1.0:
                self.T1 = 1.0
            if self.T2 == -1.0:
                self.T2 = 1.0
            if self.T3 == -1.0:
                self.T3 = 1.0
            if self.T4 == -1.0:
                self.T4 = 1.0

        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 4, 'CQUAD4'

    def cross_reference(self, model):
        msg = ' which is required by CQUAD4 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.nodeIDs()
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i,nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            if not self.pid.type in ['PLPLANE']:
                t = self.Thickness()
                assert isinstance(t, float), 'thickness=%r' % t
                mass = self.Mass()
                assert isinstance(mass, float), 'mass=%r' % mass
            a,c,n = self.AreaCentroidNormal()
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                assert isinstance(n[i], float)

    def flipNormal(self):
        r"""
        ::

          1---2       1---4
          |   |  -->  |   |
          |   |       |   |
          4---3       2---3
        """
        (n1, n2, n3, n4) = self.nodes
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

    def reprFields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3,
            T4) = self.getReprDefaults()

        list_fields = (['CQUAD4', self.eid, self.Pid()] + self.nodeIDs() +
                  [thetaMcid, zOffset, None, TFlag, T1, T2, T3, T4])
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        if size == 8 or len(card) == 7: # to last node
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class CQUADR(QuadShell):
    type = 'CQUADR'
    #calculixType = 'CAX8'

    def __init__(self, card=None, data=None, comment=''):
        QuadShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element ID
            self.eid = integer(card, 1, 'eid')
            #: Property ID
            self.pid = integer(card, 2, 'pid')
            nids = [integer_or_blank(card, 3, 'n1'),
                    integer_or_blank(card, 4, 'n2'),
                    integer_or_blank(card, 5, 'n3'),
                    integer_or_blank(card, 6, 'n4')]

            self.thetaMcid = integer_double_or_blank(card, 7, 'thetaMcid', 0.0)
            self.zOffset = double_or_blank(card, 8, 'zOffset', 0.0)

            self.TFlag = integer_or_blank(card, 10, 'TFlag', 0)
            self.T1 = double_or_blank(card, 11, 'T1', 1.0)
            self.T2 = double_or_blank(card, 12, 'T2', 1.0)
            self.T3 = double_or_blank(card, 13, 'T3', 1.0)
            self.T4 = double_or_blank(card, 14, 'T4', 1.0)
            assert len(card) <= 15, 'len(CQUADR card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:6]

            self.thetaMcid = data[6]
            self.zOffset = data[7]
            self.TFlag = data[8]
            self.T1 = data[9]
            self.T2 = data[10]
            self.T3 = data[11]
            self.T4 = data[12]
            if self.T1 == -1.0:
                self.T1 = 1.0
            if self.T2 == -1.0:
                self.T2 = 1.0
            if self.T3 == -1.0:
                self.T3 = 1.0
            if self.T4 == -1.0:
                self.T4 = 1.0
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 4, 'CQUADR'

    def cross_reference(self, model):
        msg = ' which is required by CQUADR eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.nodeIDs()

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PSHELL', 'PCOMP'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            t = self.Thickness()
            a,c,n = self.AreaCentroidNormal()
            assert isinstance(t, float), 'thickness=%r' % t
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                #assert isinstance(n[i], float)
            mass = self.Mass()
            assert isinstance(mass, float), 'mass=%r' % mass

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def flipNormal(self):
        r"""
        ::

          1---2       1---4
          |   |  -->  |   |
          |   |       |   |
          4---3       2---3
        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=True)

    def rawFields(self):
        list_fields = (['CQUADR', self.eid, self.Pid()] + self.nodeIDs() +
                  [self.thetaMcid, self.zOffset, None, self.TFlag, self.T1,
                   self.T2, self.T3, self.T4])
        return list_fields

    def reprFields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3, T4) = self.getReprDefaults()
        list_fields = (['CQUADR', self.eid, self.Pid()] + self.nodeIDs() +
                  [thetaMcid, zOffset, None, TFlag, T1, T2, T3, T4])
        return list_fields


class CQUAD(QuadShell):
    type = 'CQUAD'

    def __init__(self, card=None, data=None, comment=''):
        QuadShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = integer(card, 1, 'eid')
        #: Property ID
        self.pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 'n1'),
                integer_or_blank(card, 4, 'n2'),
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5'),
                integer_or_blank(card, 8, 'n6'),
                integer_or_blank(card, 9, 'n7'),
                integer_or_blank(card, 10, 'n8'),
                integer_or_blank(card, 11, 'n9')]
        assert len(card) <= 12, 'len(CQUAD card) = %i' % len(card)
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 9

    def cross_reference(self, model):
        msg = ' which is required by CQUAD eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def flipNormal(self):
        r"""
        ::

          1--5--2       1--8--4
          |     |  -->  |     |
          8  9  6       5  9  7
          |     |       |     |
          4--7--3       2--6--3
        """
        (n1, n2, n3, n4, n5, n6, n7, n8, n9) = self.nodes
        self.nodes = [n1, n4, n3, n2, n8, n7, n6, n5, n9]
        assert len(self.nodes) == 9

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=True)

    def rawFields(self):
        list_fields = ['CQUAD', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields

    def reprFields(self):
        list_fields = ['CQUAD', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields


class CQUAD8(QuadShell):
    type = 'CQUAD8'
    asterType = 'QUAD8'

    def __init__(self, card=None, data=None, comment=''):
        QuadShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element ID
            self.eid = integer(card, 1, 'eid')
            #: Property ID
            self.pid = integer(card, 2, 'pid')
            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2'),
                    integer(card, 5, 'n3'),
                    integer(card, 6, 'n4'),
                    integer_or_blank(card, 7, 'n5', 0),
                    integer_or_blank(card, 8, 'n6', 0),
                    integer_or_blank(card, 9, 'n7', 0),
                    integer_or_blank(card, 10, 'n8', 0)]

            self.T1 = double_or_blank(card, 11, 'T1', 1.0)
            self.T2 = double_or_blank(card, 12, 'T2', 1.0)
            self.T3 = double_or_blank(card, 13, 'T3', 1.0)
            self.T4 = double_or_blank(card, 14, 'T4', 1.0)
            self.thetaMcid = integer_double_or_blank(card, 15, 'thetaMcid', 0.0)
            self.zOffset = double_or_blank(card, 16, 'zOffset', 0.0)
            self.TFlag = integer_or_blank(card, 17, 'TFlag', 0)
            assert len(card) <= 18, 'len(CQUAD4 card) = %i' % len(card)
        else:
            #print "CQUAD8 = ",data
            #(6401,
            #6400,
            #6401, 6402, 6405, 6403, 0, 0, 6404, 0,
            #-1.0, -1.0, -1.0, -1.0,
            #0.0, 0)
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:10]
            self.T1 = data[10]
            self.T2 = data[11]
            self.T3 = data[12]
            self.T4 = data[13]
            self.thetaMcid = data[14]
            self.zOffset = data[14]
            self.TFlag = data[15]

        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 8

    def cross_reference(self, model):
        msg = ' which is required by CQUAD8 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.nodeIDs()

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PSHELL', 'PCOMP'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            t = self.Thickness()
            a,c,n = self.AreaCentroidNormal()
            assert isinstance(t, float), 'thickness=%r' % t
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                #assert isinstance(n[i], float)
            mass = self.Mass()
            assert isinstance(mass, float), 'mass=%r' % mass

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def flipNormal(self):
        r"""
        ::

          1--5--2       1--8--4
          |     |  -->  |     |
          8     6       5     7
          |     |       |     |
          4--7--3       2--6--3
        """
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.nodes
        self.nodes = [n1, n4, n3, n2, n8, n7, n6, n5]

    def Normal(self):
        (n1, n2, n3, n4) = self.nodePositions()[:4]
        return _normal(n1 - n3, n2 - n4)

    def AreaCentroid(self):
        """
        ::

          1-----2
          |    /|
          | A1/ |
          |  /  |
          |/ A2 |
          4-----3

          centroid
             c = sum(ci*Ai)/sum(A)
             where:
               c=centroid
               A=area
        """
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.nodePositions()
        a = n1 - n2
        b = n2 - n4
        area1 = Area(a, b)
        c1 = centroid_triangle(n1, n2, n4)

        a = n2 - n4
        b = n2 - n3
        area2 = Area(a, b)
        c2 = centroid_triangle(n2, n3, n4)

        area = area1 + area2
        centroid = (c1 * area1 + c2 * area2) / area
        return(area, centroid)

    def Area(self):
        r"""
        \f[ A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert \f]
        where a and b are the quad's cross node point vectors
        """
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.nodePositions()
        a = n1 - n3
        b = n2 - n4
        area = Area(a, b)
        return area

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=True)

    def rawFields(self):
        list_fields = ['CQUAD8', self.eid, self.Pid()] + self.nodeIDs() + [
            self.T1, self.T2, self.T3, self.T4, self.thetaMcid, self.zOffset,
            self.TFlag]
        return list_fields

    def reprFields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3, T4) = self.getReprDefaults()
        list_fields = (['CQUAD8', self.eid, self.Pid()] + self.nodeIDs() + [
            T1, T2, T3, T4, thetaMcid, zOffset, TFlag])
        return list_fields


class CQUADX(QuadShell):
    type = 'CQUADX'
    calculixType = 'CAX8'

    def __init__(self, card=None, data=None, comment=''):
        QuadShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element ID
            self.eid = integer(card, 1, 'eid')
            #: Property ID
            self.pid = integer(card, 2, 'pid')
            nids = [integer_or_blank(card, 3, 'n1'),
                    integer_or_blank(card, 4, 'n2'),
                    integer_or_blank(card, 5, 'n3'),
                    integer_or_blank(card, 6, 'n4'),
                    integer_or_blank(card, 7, 'n5'),
                    integer_or_blank(card, 8, 'n6'),
                    integer_or_blank(card, 9, 'n7'),
                    integer_or_blank(card, 10, 'n8'),
                    integer_or_blank(card, 11, 'n9')]
            assert len(card) <= 12, 'len(CQUADX card) = %i' % len(card)
        else:
            raise NotImplementedError(data)
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 9

    def cross_reference(self, model):
        msg = ' which is required by CQUADX eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def flipNormal(self):
        r"""
        ::

          1--5--2       1--8--4
          |     |  -->  |     |
          8  9  6       5  9  7
          |     |       |     |
          4--7--3       2--6--3
        """
        (n1, n2, n3, n4, n5, n6, n7, n8, n9) = self.nodes
        self.nodes = [n1, n4, n3, n2, n8, n7, n6, n5, n9]

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=True)

    def rawFields(self):
        list_fields = ['CQUADX', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields

    def reprFields(self):
        return self.rawFields()