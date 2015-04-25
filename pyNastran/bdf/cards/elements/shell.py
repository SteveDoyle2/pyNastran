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
from six.moves import range

from numpy import array, eye, cross, allclose
from numpy.linalg import det, norm  # inv

from pyNastran.bdf.deprecated import ShellElementDeprecated
from pyNastran.bdf.field_writer_8 import set_blank_if_default, set_default_if_blank
from pyNastran.bdf.cards.baseCard import Element
from pyNastran.utils.mathematics import Area, norm, centroid_triangle
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)
from pyNastran.bdf.field_writer_8 import print_card_8, print_field_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.utils import wipe_empty_fields

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


class ShellElement(Element, ShellElementDeprecated):
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

        .. todo:: possibly remove this
        """
        return self.pid.mid()

    def Mid(self):
        """
        Returns the material ID

        .. todo:: possibly remove this
        """
        return self.pid.Mid()

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
        A = self.Area()
        mpa = self.pid.MassPerArea()
        try:
            return mpa * A
        except TypeError:
            msg = 'mass/area=%s area=%s pidType=%s' % (mpa, A, self.pid.type)
            raise TypeError(msg)

    def flipNormal(self):
        raise NotImplementedError('flipNormal undefined for %s' % self.type)


class TriShell(ShellElement):
    def __init__(self, card, data):
        ShellElement.__init__(self, card, data)

    def get_edges(self):
        """
        Returns the edges
        """
        return [(self.nodes[0], self.nodes[1]),
                (self.nodes[1], self.nodes[2]),
                (self.nodes[2], self.nodes[0])]

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def AreaCentroidNormal(self):
        """
        Returns area,centroid, normal as it's more efficient to do them
        together

        Returns
        -------
        area : float
               the area
        centroid : (3,) array
               the centroid
        normal : (3,) array
               the normal vector
        """
        (n0, n1, n2) = self.nodePositions()
        return _triangle_area_centroid_normal([n0, n1, n2])

    def Area(self):
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} \lvert (n_0-n_1) \times (n_0-n_2) \rvert"""
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
        (n1, n2, n3) = self.nodePositions()
        try:
            n = _normal(n1 - n2, n1 - n3)
        except:
            msg = 'ERROR computing normal vector for eid=%i.\n' % self.eid
            msg += '  nid1=%i n1=%s\n' % (self.nodes[0].nid, n1)
            msg += '  nid2=%i n2=%s\n' % (self.nodes[1].nid, n2)
            msg += '  nid3=%i n3=%s\n' % (self.nodes[2].nid, n3)
            raise RuntimeError(msg)

        return n

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
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)

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

        self.prepare_node_ids(nids)
        assert len(self.nodes) == 3

    def cross_reference(self, model):
        msg = ' which is required by CTRIA3 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            if not self.pid.type in ['PLPLANE']:
                t = self.Thickness()
                assert isinstance(t, float), 'thickness=%r' % t
                mass = self.Mass()
                assert isinstance(mass, float), 'mass=%r' % mass
            a, c, n = self.AreaCentroidNormal()
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
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=False)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = (['CTRIA3', self.eid, self.Pid()] + self.node_ids +
                       [self.thetaMcid, self.zOffset, None] +
                       [None, self.TFlag, self.T1, self.T2, self.T3])
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3) = self.getReprDefaults()
        list_fields = (['CTRIA3', self.eid, self.Pid()] + self.node_ids +
                       [thetaMcid, zOffset, None] + [None, TFlag, T1, T2, T3])
        return list_fields

    def write_card(self, size=8, is_double=False):
        zOffset = set_blank_if_default(self.zOffset, 0.0)
        TFlag = set_blank_if_default(self.TFlag, 0)
        thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)

        #return self.write_card(size, double)
        nodes = self.node_ids
        row2_data = [thetaMcid, zOffset,
                     TFlag, T1, T2, T3]
        row2 = [print_field_8(field) for field in row2_data]
        data = [self.eid, self.Pid()] + nodes + row2
        msg = ('CTRIA3  %8i%8i%8i%8i%8i%8s%8s\n'
               '                %8s%8s%8s%8s\n' % tuple(data))
        return self.comment() + msg.rstrip() + '\n'

    #def write_card(self, size=8, is_double=False):
        #card = wipe_empty_fields(self.repr_fields())
        #if size == 8 or len(card) == 6: # to last node
            #msg = self.comment() + print_card_8(card)
        #else:
            #msg = self.comment() + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        #return msg

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
        self.prepare_node_ids(nids, allow_empty_nodes=True)
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
        nids = self.node_ids

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
            a, c, n = self.AreaCentroidNormal()
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

        .. math:: A = \frac{1}{2} (n_0-n_1) \times (n_0-n_2)"""
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
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=True)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = (['CTRIA6', self.eid, self.Pid()] + self.node_ids +
                       [self.thetaMcid, self.zOffset, None] + [None, self.TFlag,
                                                               self.T1, self.T2, self.T3])
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3) = self.getReprDefaults()
        list_fields = (['CTRIA6', self.eid, self.Pid()] + self.node_ids +
                       [thetaMcid, zOffset, None] + [None, TFlag, T1, T2, T3])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 8: # to last node
            msg = self.comment() + print_card_8(card)
        else:
            msg = self.comment() + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg


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

        self.prepare_node_ids(nids)
        assert len(self.nodes) == 3

        self.thetaMcid = integer_double_or_blank(card, 6, 'thetaMcid', 0.0)
        self.zOffset = double_or_blank(card, 7, 'zOffset', 0.0)
        blank(card, 8, 'blank')
        blank(card, 9, 'blank')
        blank(card, 10, 'blank')

        self.TFlag = integer_or_blank(card, 10, 'TFlag', 0)
        self.T1 = double_or_blank(card, 11, 'T1', 1.0)
        self.T2 = double_or_blank(card, 12, 'T2', 1.0)
        self.T3 = double_or_blank(card, 13, 'T3', 1.0)
        assert len(card) <= 14, 'len(CTRIAR card) = %i' % len(card)

    def cross_reference(self, model):
        msg = ' which is required by CTRIAR eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    #def _verify(self, xref=False):
        #eid = self.Eid()
        #pid = self.Pid()
        #nids = self.node_ids

        #assert isinstance(eid, int)
        #assert isinstance(pid, int)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        #if xref:
            #assert self.pid.type in ['PSHELL', 'PCOMP'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            #t = self.Thickness()
            #a,c,n = self.AreaCentroidNormal()
            #for i in range(3):
                #assert isinstance(c[i], float)
                #assert isinstance(n[i], float)

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
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            # PSHELL/PCOMP
            assert self.pid.type in ['PSHELL', 'PCOMP', 'PCOMPG'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            t = self.Thickness()
            a, c, n = self.AreaCentroidNormal()
            assert isinstance(t, float), 'thickness=%r' % t
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                #assert isinstance(n[i], float)
            mass = self.Mass()
            assert isinstance(mass, float), 'mass=%r' % mass

    def nodeIDs(self):
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=False)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = (['CTRIAR', self.eid, self.Pid()] + self.node_ids +
                       [self.thetaMcid, self.zOffset, self.TFlag,
                        self.T1, self.T2, self.T3])
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3) = self.getReprDefaults()
        list_fields = (['CTRIAR', self.eid, self.Pid()] + self.node_ids +
                       [thetaMcid, zOffset, None, None, TFlag, T1, T2, T3])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 5: # to last node
            msg = self.comment() + print_card_8(card)
        else:
            msg = self.comment() + print_card_16(card)
        return msg


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
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(nids) == 6, 'error on CTRIAX'

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            if i < 3:
                assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
            else:
                assert isinstance(nid, int) or nid is None, 'nid%i is not an integer or None nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PLPLANE'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            if not self.pid.type in ['PLPLANE']:
                t = self.Thickness()
                assert isinstance(t, float), 'thickness=%r' % t
                mass = self.Mass()
                assert isinstance(mass, float), 'mass=%r' % mass
            a, c, n = self.AreaCentroidNormal()
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                assert isinstance(n[i], float)

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

        .. math:: A = \frac{1}{2} \lvert (n_1-n_2) \times (n_1-n_3) \rvert"""
        (n1, n2, n3, n4, n5, n6) = self.nodePositions()
        a = n1 - n2
        b = n1 - n3
        area = Area(a, b)
        return area

    def cross_reference(self, model):
        msg = ' which is required by CTRIAX eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def nodeIDs(self):
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=True)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        nodeIDs = self.node_ids
        list_fields = ['CTRIAX', self.eid, self.Pid()] + nodeIDs + [self.thetaMcid]
        return list_fields

    def repr_fields(self):
        thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)
        nodeIDs = self.node_ids
        list_fields = ['CTRIAX', self.eid, self.Pid()] + nodeIDs + [thetaMcid]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 8: # to last node
            msg = self.comment() + print_card_8(card)
        else:
            msg = self.comment() + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg


class CTRIAX6(TriShell):
    r"""
    Nodes defined in a non-standard way::

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
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(nids) == 6, 'error on CTRIAX6'

    def cross_reference(self, model):
        msg = ' which is required by CTRIAX6 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.mid = model.Material(self.mid)

    def _verify(self, xref=True):
        eid = self.Eid()
        #pid = self.Pid()
        nids = self.node_ids
        assert self.pid == 0, 'pid = %s' % self.pid
        assert isinstance(eid, int)
        #assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer or blank; nid=%s' %(i, nid)

        if xref:
            assert self.mid.type in ['MAT1'], 'self.mid=%s self.mid.type=%s' % (self.mid, self.mid.type)
            #assert self.pid.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            #if not self.pid.type in ['PLPLANE']:
                #t = self.Thickness()
                #assert isinstance(t, float), 'thickness=%r' % t
                #mass = self.Mass()
                #assert isinstance(mass, float), 'mass=%r' % mass
            a, c, n = self.AreaCentroidNormal()
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                assert isinstance(n[i], float)

    def Pid(self):
        raise AttributeError("CTRIAX6 doesn't have a Property")

    def AreaCentroidNormal(self):
        """
        Returns area, centroid, normal as it's more efficient to do them
        together
        """
        (n0, n1, n2, n3, n4, n5) = self.nodePositions()
        return _triangle_area_centroid_normal([n0, n2, n4])

    def Area(self):
        r"""
        Get the normal vector.

        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_1-n_5) \rvert"""
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
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=True)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = (['CTRIAX6', self.eid, self.Mid(), self.Pid()] +
                       self.node_ids +  [self.theta])
        return list_fields

    def repr_fields(self):
        theta = set_default_if_blank(self.theta, 0.0)
        list_fields = ['CTRIAX6', self.eid, self.Mid()] + self.node_ids + [theta]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 8: # to last node
            msg = self.comment() + print_card_8(card)
        else:
            msg = self.comment() + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg


class QuadShell(ShellElement):
    def __init__(self, card, data):
        ShellElement.__init__(self, card, data)

    def get_edges(self):
        """
        Returns the edges
        """
        return [
            (self.nodes[0], self.nodes[1]),
            (self.nodes[1], self.nodes[2]),
            (self.nodes[2], self.nodes[3]),
            (self.nodes[3], self.nodes[0]),
        ]

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid.Thickness()

    def Normal(self):
        (n1, n2, n3, n4) = self.nodePositions()
        try:
            n = _normal(n1 - n3, n2 - n4)
        except:
            msg = 'ERROR computing normal vector for eid=%i.\n' % self.eid
            msg += '  nid1=%i n1=%s\n' % (self.nodes[0].nid, n1)
            msg += '  nid2=%i n2=%s\n' % (self.nodes[1].nid, n2)
            msg += '  nid3=%i n3=%s\n' % (self.nodes[2].nid, n3)
            msg += '  nid4=%i n4=%s\n' % (self.nodes[3].nid, n4)
            raise RuntimeError(msg)
        return n

    def AreaCentroidNormal(self):
        (area, centroid) = self.AreaCentroid()
        normal = self.Normal()
        return (area, centroid, normal)

    def AreaCentroid(self):
        r"""
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
        area = 0.5 * norm(cross(n3-n1, n4-n2))
        centroid = (n1 + n2 + n3 + n4) / 4.
        return(area, centroid)

    def Centroid(self):
        (n1, n2, n3, n4) = self.nodePositions()
        centroid = (n1 + n2 + n3 + n4) / 4.
        return centroid

    def Area(self):
        """
        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert
        where a and b are the quad's cross node point vectors"""
        (n1, n2, n3, n4) = self.nodePositions()
        area = 0.5 * norm(cross(n3-n1, n4-n2))
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

        #if 0:
            #print("eid       = %s" % self.eid)
            #print("nodes     = %s" % self.nodes)

            #print("self.zOffset   = %s" % self.zOffset)
            #print("self.TFlag     = %s" % self.TFlag)
            #print("self.thetaMcid = %s" % self.thetaMcid)

            #print("zOffset   = %s" % zOffset)
            #print("TFlag     = %s" % TFlag)
            #print("thetaMcid = %s" % thetaMcid)

            #print("T1 = %s" % T1)
            #print("T2 = %s" % T2)
            #print("T3 = %s" % T3)
            #print("T4 = %s\n" % T4)
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
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            nids = [integer_or_blank(card, 3, 'n1'),
                    integer_or_blank(card, 4, 'n2'),
                    integer_or_blank(card, 5, 'n3'),
                    integer_or_blank(card, 6, 'n4')]
            assert len(card) <= 7, 'len(CSHEAR card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 4

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
        r"""
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
        (area, centroid) = self.AreaCentroid()
        return centroid

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PSHEAR'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            if self.pid.type in ['PSHEAR']:
                t = self.Thickness()
                assert isinstance(t, float), 'thickness=%r' % t
                mass = self.Mass()
                assert isinstance(mass, float), 'mass=%r' % mass
            a, c, n = self.AreaCentroidNormal()
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                assert isinstance(n[i], float)

    def Area(self):
        r"""
        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert
        where a and b are the quad's cross node point vectors"""
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
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=False)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = ['CSHEAR', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        msg = self.comment() + print_card_8(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg

    def G(self):
        return self.pid.mid.G()

    def Thickness(self):
        return self.pid.t


class CQUAD4(QuadShell):
    type = 'CQUAD4'
    asterType = 'QUAD4 # CQUAD4'
    calculixType = 'S4'
    _field_map = {1: 'eid', 2:'pid', 7:'thetaMcid', 8:'zOffset',
                  10:'TFlag', 11:'T1', 12:'T2', 13:'T3'}

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
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
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

        self.prepare_node_ids(nids)
        assert len(self.nodes) == 4, 'CQUAD4'

    def cross_reference(self, model):
        msg = ' which is required by CQUAD4 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            if not self.pid.type in ['PLPLANE']:
                t = self.Thickness()
                assert isinstance(t, float), 'thickness=%r' % t
                mass = self.Mass()
                assert isinstance(mass, float), 'mass=%r' % mass
            a, c, n = self.AreaCentroidNormal()
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
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=False)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

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

    def raw_fields(self):
        list_fields = (['CQUAD4', self.eid, self.Pid()] + self.node_ids +
                       [self.thetaMcid, self.zOffset, self.TFlag, self.T1, self.T2,
                        self.T3, self.T4])
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3, T4) = self.getReprDefaults()
        list_fields = (['CQUAD4', self.eid, self.Pid()] + self.node_ids +
                       [thetaMcid, zOffset, None, TFlag, T1, T2, T3, T4])
        return list_fields

    def write_card(self, size=8, is_double=False):
        zOffset = set_blank_if_default(self.zOffset, 0.0)
        TFlag = set_blank_if_default(self.TFlag, 0)
        thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        T4 = set_blank_if_default(self.T3, 1.0)

        nodes = self.node_ids
        row2_data = [thetaMcid, zOffset,
                     TFlag, T1, T2, T3, T4]
        row2 = [print_field_8(field) for field in row2_data]
        data = [self.eid, self.Pid()] + nodes + row2
        msg = ('CQUAD4  %8i%8i%8i%8i%8i%8i%8s%8s\n'
               '                %8s%8s%8s%8s%8s\n' % tuple(data))
        return self.comment() + msg.rstrip() + '\n'

    #def write_card(self, size=8, is_double=False):
        #card = wipe_empty_fields(self.repr_fields())
        #if size == 8 or len(card) == 7: # to last node
            #msg = self.comment() + print_card_8(card)
        #else:
            #msg = self.comment() + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        #return msg


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
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 4, 'CQUADR'

    def cross_reference(self, model):
        msg = ' which is required by CQUADR eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PSHELL', 'PCOMP', 'PCOMPG'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            t = self.Thickness()
            a, c, n = self.AreaCentroidNormal()
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
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=True)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = (['CQUADR', self.eid, self.Pid()] + self.node_ids +
                       [self.thetaMcid, self.zOffset, None, self.TFlag, self.T1,
                        self.T2, self.T3, self.T4])
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3, T4) = self.getReprDefaults()
        list_fields = (['CQUADR', self.eid, self.Pid()] + self.node_ids +
                       [thetaMcid, zOffset, None, TFlag, T1, T2, T3, T4])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8 or len(card) == 7: # to last node
            msg = self.comment() + print_card_8(card)
        else:
            msg = self.comment() + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg


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
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5'),
                integer_or_blank(card, 8, 'n6'),
                integer_or_blank(card, 9, 'n7'),
                integer_or_blank(card, 10, 'n8'),
                integer_or_blank(card, 11, 'n9')]
        assert len(card) <= 12, 'len(CQUAD card) = %i' % len(card)
        self.prepare_node_ids(nids, allow_empty_nodes=True)
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
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=True)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = ['CQUAD', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def repr_fields(self):
        list_fields = ['CQUAD', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%8i' % node for node in nodes[4:]]
        data = [self.eid, self.Pid()] + nodes[:4] + nodes2
        msg = ('CQUAD   %8i%8i%8i%8i%8i%8i%8s%8s\n'  # 6 nodes
               '        %8s%8s%8s\n' % tuple(data))
        return self.comment() + msg.rstrip() + '\n'

    #def write_card(self, size=8, is_double=False):
        #card = self.repr_fields()
        #msg = self.comment() + print_card_8(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        #return msg


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

        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 8

    def cross_reference(self, model):
        msg = ' which is required by CQUAD8 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allowEmptyNodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid.type in ['PSHELL', 'PCOMP', 'PCOMPG'], 'pid=%i self.pid.type=%s' % (pid, self.pid.type)
            t = self.Thickness()
            a, c, n = self.AreaCentroidNormal()
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
        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert
        where a and b are the quad's cross node point vectors"""
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.nodePositions()
        a = n1 - n3
        b = n2 - n4
        area = Area(a, b)
        return area

    def nodeIDs(self):
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=True)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = ['CQUAD8', self.eid, self.Pid()] + self.node_ids + [
            self.T1, self.T2, self.T3, self.T4, self.thetaMcid, self.zOffset,
            self.TFlag]
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3, T4) = self.getReprDefaults()
        list_fields = (['CQUAD8', self.eid, self.Pid()] + self.node_ids + [
            T1, T2, T3, T4, thetaMcid, zOffset, TFlag])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8 or len(card) == 11: # to last node
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


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
        self.prepare_node_ids(nids, allow_empty_nodes=True)
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
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allowEmptyNodes=True)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = ['CQUADX', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes[:4]
        row2 = ['        ' if node is None else '%8i' % node for node in nodes[4:]]
        msg = ('CQUADX  %8i%8i%8i%8i%8i%8i%8s%8s\n'
               '        %8s%8s%8s' % tuple(data + row2))
        return self.comment() + msg.rstrip() + '\n'

