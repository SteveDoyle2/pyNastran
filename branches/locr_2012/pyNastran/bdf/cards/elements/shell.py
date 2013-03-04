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

from numpy import array, eye, cross, allclose  # zeros,dot
from numpy.linalg import det  # inv

from pyNastran.bdf.fieldWriter import (set_blank_if_default,
                                       set_default_if_blank)
from pyNastran.bdf.cards.baseCard import Element
from pyNastran.utils.mathematics import Area, norm, centroid_triangle
from pyNastran.bdf.format import (integer, integer_or_blank,
                                  double_or_blank, integer_double_or_blank)


def _triangle_area_centroid_normal(nodes):
    """
    Returns area,centroid,unitNormal

    @param nodes:
      list of three triangle vertices

    @code
    n = Normal = a x b
    Area   = 1/2 * |a x b|
    V = <v1,v2,v3>
    |V| = sqrt(v1^0.5+v2^0.5+v3^0.5) = norm(V)

    Area = 0.5 * |n|
    unitNormal = n/|n|
    @endcode
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
        \f[ \large  mass = \frac{mass}{area} area  \f]
        """
        return self.pid.MassPerArea() * self.Area()

    def flipNormal(self):
        raise NotImplementedError('flipNormal undefined for %s' % (self.type))

    def cross_reference(self, mesh):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = mesh.Nodes(self.nodes, msg=msg)
        self.pid = mesh.Property(self.pid, msg=msg)


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
        Returns the normal vector
        \f[ \large A = \frac{1}{2} (n_0-n_1) \times (n_0-n_2)  \f]
        """
        (n0, n1, n2) = self.nodePositions()
        a = n0 - n1
        b = n0 - n2
        area = Area(a, b)
        return area

    def Normal(self):
        r"""
        returns the normal vector
        \f[ \large a = (n_0-n_1) \times (n_0-n_2)  \f]
        \f[ \large n = \frac{n}{norm(N)}           \f]
        """
        (n0, n1, n2) = self.nodePositions()
        return _normal(n0 - n1, n0 - n2)

    def Centroid(self):
        r"""
        Returns the centroid
        \f[ \large CG = \frac{1}{3} (n_0+n_1+n_2)  \f]
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

    def __init__(self, card=None, data=None, comment=''):
        TriShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            ## element ID number
            self.eid = integer(card, 1, 'eid')
            self.pid = integer(card, 2, 'pid')

            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2'),
                    integer(card, 5, 'n3')]

            self.thetaMcid = integer_double_or_blank(card, 6, 'thetaMcid', 0.0)
            self.zOffset = double_or_blank(card, 7, 'zOffset', 0.0)

            self.TFlag = integer_or_blank(card, 10, 'TFlag', 0)
            self.T1 = double_or_blank(card, 11, 'T1', 1.0)
            self.T2 = double_or_blank(card, 12, 'T2', 1.0)
            self.T3 = double_or_blank(card, 13, 'T3', 1.0)
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

    def flipNormal(self):
        """
        @code
             1           1
            * *   -->   * *
           *   *       *   *
          2-----3     3-----2
        @endcode
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


class CTRIA6(TriShell):
    type = 'CTRIA6'
    asterType = 'TRIA6'
    calculixType = 'S6'

    def __init__(self, card=None, data=None, comment=''):
        TriShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            ## element ID number
            self.eid = integer(card, 1, 'eid')
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
        (n1, n2, n3, n4, n5, n6) = self.nodePositions()
        return _triangle_area_centroid_normal([n1, n2, n3])

    def Area(self):
        r"""
        Returns the normal vector
        \f[ \large A = \frac{1}{2} (n_0-n_1) \times (n_0-n_2)  \f]
        """
        (n1, n2, n3, n4, n5, n6) = self.nodePositions()
        a = n1 - n2
        b = n1 - n3
        area = Area(a, b)
        return area

    def Normal(self):
        r"""
        Returns the normal vector
        \f[ \large a = (n_0-n_1) \times (n_0-n_2)  \f]
        \f[ \large n = \frac{n}{norm(N)}           \f]
        """
        (n0, n1, n2) = self.nodePositions()[:3]
        return _normal(n0 - n1, n0 - n2)

    def Centroid(self):
        r"""
        Returns the centroid
        \f[ \large CG = \frac{1}{3} (n_1+n_2+n_3)  \f]
        """
        (n1, n2, n3, n4, n5, n6) = self.nodePositions()
        centroid = centroid_triangle(n1, n2, n3)
        return centroid

    def flipNormal(self):
        """
        @code
             1                1
             **               **
            *  *             *  *
           4    6   -->     6    4
          *      *         *      *
         2----5---3       3----5---2
        @endcode
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

    def rawFields(self):
        list_fields = (['CTRIA6', self.eid, self.Pid()] + self.nodeIDs(allowEmptyNodes=True) +
                  [self.thetaMcid, self.zOffset, None] + [None, self.TFlag,
                  self.T1, self.T2, self.T3])
        return list_fields

    def reprFields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3) = self.getReprDefaults()
        list_fields = (['CTRIA6', self.eid, self.Pid()] + self.nodeIDs(allowEmptyNodes=True) +
                  [thetaMcid, zOffset, None] + [None, TFlag, T1, T2, T3])
        return list_fields


class CTRIAR(TriShell):
    type = 'CTRIAR'
    def __init__(self, card=None, data=None, comment=''):
        TriShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        ## element ID number
        self.eid = integer(card, 1, 'eid')
        self.pid = integer(card, 2, 'pid')

        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3')]

        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 3

        self.thetaMcid = integer_double_or_blank(card, 6, 'thetaMcid', 0.0)
        self.zOffset = double_or_blank(card, 7, 'zOffset', 0.0)

        self.TFlag = integer_or_blank(card, 11, 'TFlag', 0)
        self.T1 = double_or_blank(card, 11, 'T1', 1.0)
        self.T2 = double_or_blank(card, 12, 'T2', 1.0)
        self.T3 = double_or_blank(card, 13, 'T3', 1.0)

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def flipNormal(self):
        """
        @code
             1           1
            * *   -->   * *
           *   *       *   *
          2-----3     3-----2
        @endcode
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
        ## element ID number
        self.eid = integer(card, 1, 'eid')

        nids = [integer_or_blank(card, 3, 'n1'),
                integer_or_blank(card, 4, 'n2'),
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5'),
                integer_or_blank(card, 8, 'n6')]
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(nids) == 6, 'error on CTRIAX'

    def rawFields(self):
        list_fields = ['CTRIAX', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields

    def reprFields(self):
        return self.rawFields()


class CTRIAX6(TriShell):
    """
    Nodes defined in a non-standard way
    @code
         5
        / \
       6   4
     /       \
    1----2----3
    @endcode
    """
    type = 'CTRIAX6'
    #calculixType = 'CAX6'
    def __init__(self, card=None, data=None, comment=''):
        TriShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        ## element ID number
        self.eid = integer(card, 1, 'eid')
        self.mid = integer(card, 2, 'mid')

        nids = [integer_or_blank(card, 3, 'n1'),
                integer_or_blank(card, 4, 'n2'),
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5'),
                integer_or_blank(card, 8, 'n6')]
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(nids) == 6, 'error on CTRIAX6'

        ## theta
        self.theta = double_or_blank(card, 10, 'theta', 0.0)

    def Area(self):
        r"""
        returns the normal vector
        \f[ \large A = \frac{1}{2} (n_0-n_1) times (n_0-n_2)  \f]
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

    def cross_reference(self, model):
        msg = ' which is required by CTRIAX6 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.mid = model.Material(self.mid)

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        return self.mid.mid

    def flipNormal(self):
        """
        @code
             5               5
            / \             / \
           6   4   -->     6   4
         /       \       /       \
        1----2----3     1----2----3
        @endcode
        """
        (n1, n2, n3, n4, n5, n6) = self.nodes
        self.nodes = [n1, n6, n5, n4, n3, n2]

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
        @code
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
        @endcode
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
        """
        @code
        1---2        1---4
        |   |  -->   |   |
        |   |        |   |
        4---3        2---3
        @endcode
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
            self.eid = integer(card, 1, 'eid')
            self.pid = integer(card, 2, 'pid')
            nids = [integer_or_blank(card, 3, 'n1'),
                    integer_or_blank(card, 4, 'n2'),
                    integer_or_blank(card, 5, 'n3'),
                    integer_or_blank(card, 6, 'n4')]
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 4

    def Normal(self):
        (n1, n2, n3, n4) = self.nodePositions()
        return _normal(n1 - n3, n2 - n4)

    def AreaCentroidNormal(self):
        (area, centroid) = self.AreaCentroid()
        normal = self.Normal()
        return (area, centroid, normal)

    def AreaCentroid(self):
        """
        @code
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
        @endcode
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
        """
        @code
        1---2        1---4
        |   |  -->   |   |
        |   |        |   |
        4---3        2---3
        @endcode
        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]

    def rawFields(self):
        list_fields = ['CSHEAR', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields

    def reprFields(self):
        return self.rawFields()


class CQUAD4(QuadShell):
    type = 'CQUAD4'
    asterType = 'QUAD4 # CQUAD4'
    calculixType = 'S4'

    def __init__(self, card=None, data=None, comment=''):
        QuadShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            ## element ID number
            self.eid = integer(card, 1, 'eid')
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

    def flipNormal(self):
        """
        @code
        1---2        1---4
        |   |  -->   |   |
        |   |        |   |
        4---3        2---3
        @endcode
        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]

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


class CQUADR(QuadShell):
    type = 'CQUADR'
    #calculixType = 'CAX8'

    def __init__(self, card=None, data=None, comment=''):
        QuadShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            ## element ID number
            self.eid = integer(card, 1, 'eid')
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

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def flipNormal(self):
        """
        @code
        1---2        1---4
        |   |  -->   |   |
        |   |        |   |
        4---3        2---3
        @endcode
        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]

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
        ## element ID number
        self.eid = integer(card, 1, 'eid')
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
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 9

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def flipNormal(self):
        """
        @code
        1--5--2        1--8--4
        |     |  -->   |     |
        8  9  6        5  9  7
        |     |        |     |
        4--7--3        2--6--3
        @endcode
        """
        (n1, n2, n3, n4, n5, n6, n7, n8, n9) = self.nodes
        self.nodes = [n1, n4, n3, n2, n8, n7, n6, n5, n9]
        assert len(self.nodes) == 9

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
            ## element ID number
            self.eid = integer(card, 1, 'eid')
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

    def cross_reference(self, mesh):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = mesh.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = mesh.Property(self.pid, msg=msg)

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def flipNormal(self):
        """
        1--5--2        1--8--4
        |     |  -->   |     |
        8     6        5     7
        |     |        |     |
        4--7--3        2--6--3
        """
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.nodes
        self.nodes = [n1, n4, n3, n2, n8, n7, n6, n5]

    def Normal(self):
        (n1, n2, n3, n4) = self.nodePositions()[:4]
        return _normal(n1 - n3, n2 - n4)

    def AreaCentroid(self):
        """
        @code
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
        @endcode
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

    def rawFields(self):
        list_fields = ['CQUAD8', self.eid, self.Pid()] + self.nodeIDs(allowEmptyNodes=True) + [
            self.T1, self.T2, self.T3, self.T4, self.thetaMcid, self.zOffset,
            self.TFlag]
        return list_fields

    def reprFields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3, T4) = self.getReprDefaults()
        list_fields = (['CQUAD8', self.eid, self.Pid()] + self.nodeIDs(allowEmptyNodes=True) + [
            T1, T2, T3, T4, thetaMcid, zOffset, TFlag])
        return list_fields


class CQUADX(QuadShell):
    type = 'CQUADX'
    calculixType = 'CAX8'

    def __init__(self, card=None, data=None, comment=''):
        QuadShell.__init__(self, card, data)
        if comment:
            self._comment = comment
        ## element ID number
        self.eid = integer(card, 1, 'eid')
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
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) == 9

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def flipNormal(self):
        """
        @code
        1--5--2        1--8--4
        |     |  -->   |     |
        8  9  6        5  9  7
        |     |        |     |
        4--7--3        2--6--3
        @endcode
        """
        (n1, n2, n3, n4, n5, n6, n7, n8, n9) = self.nodes
        self.nodes = [n1, n4, n3, n2, n8, n7, n6, n5, n9]

    def rawFields(self):
        list_fields = ['CQUADX', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields

    def reprFields(self):
        return self.rawFields()