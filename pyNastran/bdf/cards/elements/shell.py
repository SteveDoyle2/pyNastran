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

import numpy as np
from numpy import cross, allclose
from numpy.linalg import norm

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, set_default_if_blank
from pyNastran.bdf.cards.base_card import Element
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, integer_double_or_blank, blank)
from pyNastran.bdf.field_writer_8 import print_card_8, print_field_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.utils import wipe_empty_fields

def _triangle_area_centroid_normal(nodes, card):
    """

    Parameters
    -------------
    nodes : list
        List of three triangle vertices.

    Returns
    --------
    area : float
        Area of triangle.
    centroid : ndarray
        Centroid of triangle.
    unit_normal : ndarray
        Unit normal of triangles.
    card : CTRIA3(), CTRIA6()
        the self parameter

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
    try:
        normal = vector / length
    except FloatingPointError as e:
        msg = e.strerror
        msg += '\nvector: %s ; length: %s' % (vector, length)
        raise RuntimeError(msg)

    if not allclose(norm(normal), 1.):
        msg = ('function _triangle_area_centroid_normal, check...\n'
               'a = {0}\nb = {1}\nnormal = {2}\nlength = {3}\n{4}'.format(
                   n0 - n1, n0 - n2, normal, length, str(card)))
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

    def __init__(self):
        Element.__init__(self)

    def Rho(self):
        """
        Returns the density
        """
        self.deprecated('Rho()', 'pid.mid().rho', '0.8')
        return self.pid_ref.mid().rho

    def Eid(self):
        return self.eid

    def Area(self):
        raise NotImplementedError('Area undefined for %s' % self.type)

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid_ref.Thickness()

    def mid(self):
        """
        Returns the material

        .. todo:: possibly remove this
        """
        return self.pid_ref.mid()

    def Mid(self):
        """
        Returns the material ID

        .. todo:: possibly remove this
        """
        return self.pid_ref.Mid()

    def Nsm(self):
        """
        Returns the non-structural mass
        """
        return self.pid_ref.Nsm()

    def MassPerArea(self):
        """
        Returns the mass per area
        """
        return self.pid_ref.MassPerArea()

    def Mass(self):
        r"""
        .. math:: m = \frac{m}{A} A  \f]
        """
        A = self.Area()
        mpa = self.pid_ref.MassPerArea()
        try:
            return mpa * A
        except TypeError:
            msg = 'mass/area=%s area=%s pidType=%s' % (mpa, A, self.pid_ref.type)
            raise TypeError(msg)

    def flipNormal(self):
        raise NotImplementedError('flipNormal undefined for %s' % self.type)


class TriShell(ShellElement):
    def __init__(self):
        ShellElement.__init__(self)

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [
            tuple(sorted([node_ids[0], node_ids[1]])),
            tuple(sorted([node_ids[1], node_ids[2]])),
            tuple(sorted([node_ids[2], node_ids[0]]))
        ]

    def get_edge_axes(self):
        n1, n2, n3 = self.nodes_ref
        g1 = n1.get_position()
        g2 = n2.get_position()
        g3 = n3.get_position()
        x = g2 - g1
        yprime = g3 - g1
        normal = cross(x, yprime)
        y = cross(normal, x)
        return x, y

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid_ref.Thickness()

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
        (n0, n1, n2) = self.get_node_positions()
        return _triangle_area_centroid_normal([n0, n1, n2], self)

    def get_area(self):
        return self.Area()

    def Area(self):
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} \lvert (n_0-n_1) \times (n_0-n_2) \rvert"""
        (n1, n2, n3) = self.get_node_positions()
        a = n1 - n2
        b = n1 - n3
        area = 0.5 * norm(cross(a, b))
        return area

    def Normal(self):
        r"""
        Get the normal vector, :math:`n`.

        .. math::
          n = \frac{(n_0-n_1) \times (n_0-n_2)}
             {\lvert (n_0-n_1) \times (n_0-n_2) \lvert}
        """
        n1, n2, n3 = self.get_node_positions()
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
        n1, n2, n3 = self.get_node_positions()
        centroid = (n1 + n2 + n3) / 3.
        return centroid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref


class CTRIA3(TriShell):
    """
    +--------+-------+-------+----+----+----+------------+---------+-----+
    |   1    |   2   |   3   |  4 |  5 |  6 |     7      |    8    |  9  |
    +========+=======+=======+=====+===+====+============+=========+=====+
    | CTRIA3 |  EID  |  PID  | N1 | N2 | N3 | THETA/MCID | ZOFFSET |     |
    +--------+-------+-------+----+----+----+------------+---------+-----+
    |        |       | TFLAG | T1 | T2 | T3 |            |         |     |
    +--------+-------+-------+----+----+----+------------+---------+-----+
    """
    type = 'CTRIA3'
    aster_type = 'TRIA3'
    calculixType = 'S3'
    _field_map = {
        1: 'eid', 2:'pid', 6:'thetaMcid', 7:'zOffset', 10:'TFlag',
        11:'T1', 12:'T2', 13:'T3'}

    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 4:
            self.nodes[1] = value
        elif n == 5:
            self.nodes[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, pid, nids, zOffset,
                 thetaMcid=0.0, TFlag=0, T1=1.0, T2=1.0, T3=1.0, comment=''):
        TriShell.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        assert len(nids) == 3, nids
        self.prepare_node_ids(nids)
        self.zOffset = zOffset
        self.thetaMcid = thetaMcid
        self.TFlag = TFlag
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 3

    def validate(self):
        assert len(set(self.nodes)) == 3, 'nodes=%s; n=%s\n%s' % (self.nodes, len(set(self.nodes)), str(self))

    @classmethod
    def add_op2_data(cls, data, comment=''):
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
        return CTRIA3(eid, pid, nids, zOffset, thetaMcid,
                      TFlag, T1, T2, T3, comment=comment)

    @classmethod
    def add_card(cls, card, comment=''):
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3')
        ]
        if len(card) > 5:
            thetaMcid = integer_double_or_blank(card, 6, 'thetaMcid', 0.0)
            zOffset = double_or_blank(card, 7, 'zOffset', 0.0)
            blank(card, 8, 'blank')
            blank(card, 9, 'blank')

            TFlag = integer_or_blank(card, 10, 'TFlag', 0)
            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            assert len(card) <= 14, 'len(CTRIA3 card) = %i\ncard=%s' % (len(card), card)
        else:
            thetaMcid = 0.0
            zOffset = 0.0
            TFlag = 0
            T1 = 1.0
            T2 = 1.0
            T3 = 1.0
        return CTRIA3(eid, pid, nids, zOffset, thetaMcid,
                      TFlag, T1, T2, T3, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CTRIA3 eid=%s' % self.eid
        self.nodes = model.Nodes(self.node_ids, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.Pid(), msg=msg)
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        edges = self.get_edge_ids()

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            if not self.pid_ref.type in ['PLPLANE']:
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
        #        t = self.pid_ref.Thickness()
        #        T1 = self.T1 / t
        #        T2 = self.T2 / t
        #        T3 = self.T3 / t
        #    else:
        #        T1 = self.T1
        #        T2 = self.T2
        #        T3 = self.T3
        #    t = (T1+T2+T3)/3.
        #else:
        #    t = self.pid_ref.Thickness()
        #return t

    def flipNormal(self):
        """
        Flips normal of element.

        ::

               1           1
              * *   -->   * *
             *   *       *   *
            2-----3     3-----2
        """
        (n1, n2, n3) = self.nodes
        self.nodes = [n1, n3, n2]

    def _get_repr_defaults(self):
        zOffset = set_blank_if_default(self.zOffset, 0.0)
        TFlag = set_blank_if_default(self.TFlag, 0)
        thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        return (thetaMcid, zOffset, TFlag, T1, T2, T3)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=False)

    def raw_fields(self):
        list_fields = (['CTRIA3', self.eid, self.Pid()] + self.node_ids +
                       [self.thetaMcid, self.zOffset, None] +
                       [None, self.TFlag, self.T1, self.T2, self.T3])
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3) = self._get_repr_defaults()
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
        return self.comment + msg.rstrip() + '\n'


class CTRIA6(TriShell):
    type = 'CTRIA6'
    aster_type = 'TRIA6'
    calculixType = 'S6'

    def __init__(self, eid, pid, nids, thetaMcid, zOffset,
                 TFlag, T1, T2, T3, comment=''):
        TriShell.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        self.thetaMcid = thetaMcid
        self.zOffset = zOffset
        self.TFlag = TFlag
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(nids) == 6, 'error on CTRIA6'

    @classmethod
    def add_card(cls, card, comment=''):
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer(card, 2, 'pid')

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4', 0),
            integer_or_blank(card, 7, 'n5', 0),
            integer_or_blank(card, 8, 'n6', 0)
        ]
        if len(card) > 9:
            thetaMcid = integer_double_or_blank(card, 9, 'thetaMcid', 0.0)
            zoffset = double_or_blank(card, 10, 'zOffset', 0.0)

            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            TFlag = integer_or_blank(card, 14, 'TFlag', 0)
            assert len(card) <= 15, 'len(CTRIA6 card) = %i\ncard=%s' % (len(card), card)
        else:
            thetaMcid = 0.0
            zoffset = 0.0
            T1 = 1.0
            T2 = 1.0
            T3 = 1.0
            TFlag = 0
        return CTRIA6(eid, pid, nids, thetaMcid, zoffset,
                      TFlag, T1, T2, T3, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = data[2:8]
        thetaMcid = data[8]
        zoffset = data[9]
        T1 = data[10]
        T2 = data[11]
        T3 = data[12]
        TFlag = data[13]
        assert isinstance(T1, float), data
        assert isinstance(T2, float), data
        assert isinstance(T3, float), data
        assert isinstance(TFlag, int), data
        #prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(nids) == 6, 'error on CTRIA6'
        return CTRIA6(eid, pid, nids, thetaMcid, zoffset,
                      TFlag, T1, T2, T3, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CTRIA6 eid=%s' % self.eid
        self.nodes = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.pid = model.Property(self.Pid(), msg=msg)
        self.nodes_ref = self.nodes
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        edges = self.get_edge_ids()

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert isinstance(nid, integer_types) or nid is None, 'nid%i is not an integer/None; nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            if not self.pid_ref.type in ['PLPLANE']:
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
        return self.pid_ref.Thickness()

    def AreaCentroidNormal(self):
        """
        Returns area, centroid, normal as it's more efficient to do them
        together
        """
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()
        return _triangle_area_centroid_normal([n1, n2, n3], self)

    def Area(self):
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} (n_0-n_1) \times (n_0-n_2)"""
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()
        a = n1 - n2
        b = n1 - n3
        area = 0.5 * norm(cross(a, b))
        return area

    def Normal(self):
        r"""
        Get the normal vector, :math:`n`.

        .. math::
          n = \frac{(n_0-n_1) \times (n_0-n_2)}{\lvert (n_0-n_1) \times (n_0-n_2) \lvert}
        """
        (n0, n1, n2) = self.get_node_positions()[:3]
        return _normal(n0 - n1, n0 - n2)

    def Centroid(self):
        r"""
        Get the centroid.

        .. math::
          CG = \frac{1}{3} (n_1+n_2+n_3)
        """
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()
        centroid = (n1 + n2 + n3) / 3.
        return centroid

    def flipNormal(self):
        r"""
        Flips normal of element.

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

    def _get_repr_defaults(self):
        zOffset = set_blank_if_default(self.zOffset, 0.0)
        assert isinstance(self.TFlag, int), self.TFlag
        TFlag = set_blank_if_default(self.TFlag, 0)
        thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        return (thetaMcid, zOffset, TFlag, T1, T2, T3)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = (['CTRIA6', self.eid, self.Pid()] + self.node_ids +
                       [self.thetaMcid, self.zOffset,
                        self.T1, self.T2, self.T3, self.TFlag,])
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3) = self._get_repr_defaults()
        list_fields = (['CTRIA6', self.eid, self.Pid()] + self.node_ids +
                       [thetaMcid, zOffset, T1, T2, T3, TFlag])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 8: # to last node
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg


class CTRIAR(TriShell):
    type = 'CTRIAR'
    def __init__(self, eid, pid, nids, thetaMcid, zOffset,
                 TFlag, T1, T2, T3, comment=''):
        TriShell.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid

        self.thetaMcid = thetaMcid
        self.zOffset = zOffset
        self.TFlag = TFlag
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.nodes = nids
        assert len(self.nodes) == 3

    def validate(self):
        self.validate_node_ids(allow_empty_nodes=False)

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3')
        ]

        thetaMcid = integer_double_or_blank(card, 6, 'thetaMcid', 0.0)
        zOffset = double_or_blank(card, 7, 'zOffset', 0.0)
        blank(card, 8, 'blank')
        blank(card, 9, 'blank')

        TFlag = integer_or_blank(card, 10, 'TFlag', 0)
        T1 = double_or_blank(card, 11, 'T1')
        T2 = double_or_blank(card, 12, 'T2')
        T3 = double_or_blank(card, 13, 'T3')
        assert len(card) <= 14, 'len(CTRIAR card) = %i\ncard=%s' % (len(card), card)
        return CTRIAR(eid, pid, nids, thetaMcid, zOffset,
                      TFlag, T1, T2, T3, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CTRIAR eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)
        self.nodes_ref = self.nodes
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    #def _verify(self, xref=False):
        #eid = self.Eid()
        #pid = self.Pid()
        #nids = self.node_ids

        #assert isinstance(eid, integer_types)
        #assert isinstance(pid, integer_types)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        #if xref:
            #assert self.pid_ref.type in ['PSHELL', 'PCOMP'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            #t = self.Thickness()
            #a,c,n = self.AreaCentroidNormal()
            #for i in range(3):
                #assert isinstance(c[i], float)
                #assert isinstance(n[i], float)

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid_ref.Thickness()

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

    def _get_repr_defaults(self):
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
        edges = self.get_edge_ids()

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            # PSHELL/PCOMP
            assert self.pid_ref.type in ['PSHELL', 'PCOMP', 'PCOMPG'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
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
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=False)

    def raw_fields(self):
        list_fields = (['CTRIAR', self.eid, self.Pid()] + self.node_ids +
                       [self.thetaMcid, self.zOffset, self.TFlag,
                        self.T1, self.T2, self.T3])
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3) = self._get_repr_defaults()
        list_fields = (['CTRIAR', self.eid, self.Pid()] + self.node_ids +
                       [thetaMcid, zOffset, None, None, TFlag, T1, T2, T3])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 5: # to last node
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        return msg


class CTRIAX(TriShell):
    """
    +--------+------------+-------+----+----+----+----+----+-----+
    |   1    |     2      |   3   |  4 |  5 |  6 | 7  |  8 |  9  |
    +========+============+=======+====+====+====+====+====+=====+
    | CTRIA3 |    EID     |  PID  | N1 | N2 | N3 | N4 | N5 | N6  |
    +--------+------------+-------+----+----+----+----+----+-----+
    |        | THETA/MCID |       |    |    |    |    |    |     |
    +--------+------------+-------+----+----+----+----+----+-----+
    """
    type = 'CTRIAX'
    calculixType = 'CAX6'
    def __init__(self, eid, pid, nids, thetaMcid, comment=''):
        TriShell.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.thetaMcid = thetaMcid
        self.nodes = nids
        assert len(nids) == 6, 'error on CTRIAX'

    def validate(self):
        self.validate_node_ids(allow_empty_nodes=True)

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')

        nids = [
            integer_or_blank(card, 3, 'n1'),
            integer_or_blank(card, 4, 'n2'),
            integer_or_blank(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4'),
            integer_or_blank(card, 7, 'n5'),
            integer_or_blank(card, 8, 'n6'),
            ]
        thetaMcid = integer_double_or_blank(card, 9, 'theta_mcsid', 0.0)
        assert len(card) <= 10, 'len(CTRIAX card) = %i\ncard=%s' % (len(card), card)
        return CTRIAX(eid, pid, nids, thetaMcid, comment=comment)

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        edges = self.get_edge_ids()

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            if i < 3:
                assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)
            else:
                assert isinstance(nid, integer_types) or nid is None, 'nid%i is not an integer or None nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PLPLANE'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            if not self.pid_ref.type in ['PLPLANE']:
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
        pass

    def AreaCentroidNormal(self):
        """
        Returns area, centroid, normal as it's more efficient to do them
        together
        """
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()
        return _triangle_area_centroid_normal([n1, n2, n3], self)

    def Area(self):
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} \lvert (n_1-n_2) \times (n_1-n_3) \rvert"""
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()
        a = n1 - n2
        b = n1 - n3
        area = 0.5 * norm(cross(a, b))
        return area

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CTRIAX eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)
        self.nodes_ref = self.nodes
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CTRIAX', self.eid, self.Pid()] + self.node_ids + [self.thetaMcid]
        return list_fields

    def repr_fields(self):
        thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)
        nodeIDs = self.node_ids
        list_fields = ['CTRIAX', self.eid, self.Pid()] + nodeIDs + [thetaMcid]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 8: # to last node
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg


class CTRIAX6(TriShell):
    """
    +--------+-------+-------+----+----+----+----+----+-----+
    |   1    |   2   |   3   |  4 |  5 |  6 |  7 |  8 |  9  |
    +========+=======+=======+=====+===+====+====+====+=====+
    | CTRIAX6 |  EID |  MID  | N1 | N2 | N3 | G4 | G5 | G6  |
    +--------+-------+-------+----+----+----+----+----+-----+
    |        |       | THETA |    |    |    |    |    |     |
    +--------+-------+-------+----+----+----+----+----+-----+

    Nodes are defined in a non-standard way::

           5
          / \
         6   4
       /       \
      1----2----3
    """
    type = 'CTRIAX6'
    #calculixType = 'CAX6'
    def __init__(self, eid, mid, nids, theta, comment=''):
        TriShell.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Material ID
        self.mid = mid
        #: theta
        self.theta = theta
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(nids) == 6, 'error on CTRIAX6'

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        mid = integer(card, 2, 'mid')

        nids = [
            integer(card, 3, 'n1'),
            integer_or_blank(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4'),
            integer(card, 7, 'n5'),
            integer_or_blank(card, 8, 'n6'),
        ]

        theta = double_or_blank(card, 9, 'theta', 0.0)
        assert len(card) <= 10, 'len(CTRIAX6 card) = %i\ncard=%s' % (len(card), card)
        return CTRIAX6(eid, mid, nids, theta, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CTRIAX6 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.mid = model.Material(self.mid)
        self.nodes_ref = self.nodes
        self.mid_ref = self.mid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.mid = self.Mid()
        del self.nodes_ref, self.mid_ref

    def _verify(self, xref=True):
        eid = self.Eid()
        nids = self.node_ids
        edges = self.get_edge_ids()

        assert self.pid == 0, 'pid = %s' % self.pid
        assert isinstance(eid, integer_types)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, integer_types), 'nid%i is not an integer or blank; nid=%s' %(i, nid)

        if xref:
            assert self.mid.type in ['MAT1', 'MAT3', 'MAT4'], 'self.mid=%s self.mid.type=%s' % (self.mid, self.mid.type)
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
        (n0, n1, n2, n3, n4, n5) = self.get_node_positions()
        return _triangle_area_centroid_normal([n0, n2, n4], self)

    def Area(self):
        r"""
        Get the normal vector.

        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_1-n_5) \rvert"""
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()
        a = n1 - n3
        b = n1 - n5
        area = 0.5 * norm(cross(a, b))
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
        if isinstance(self.mid, integer_types):
            return self.mid
        return self.mid_ref.mid

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
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        """
             5
            / \
           6   4
         /       \
        1----2----3
        """
        return self._nodeIDs(allow_empty_nodes=True)

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [
            tuple(sorted([node_ids[0], node_ids[2]])),
            tuple(sorted([node_ids[2], node_ids[4]])),
            tuple(sorted([node_ids[4], node_ids[0]]))
        ]

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
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg


class QuadShell(ShellElement):
    def __init__(self):
        ShellElement.__init__(self)

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [
            tuple(sorted([node_ids[0], node_ids[1]])),
            tuple(sorted([node_ids[1], node_ids[2]])),
            tuple(sorted([node_ids[2], node_ids[3]])),
            tuple(sorted([node_ids[3], node_ids[0]]))
        ]

    def get_edge_number_by_node_ids(self, n1, n2):
        edge_ids = self.get_edge_ids()
        edge = [n1, n2]
        edge.sort()
        tedge = tuple(edge)
        iedge = edge_ids.index(tedge)
        return iedge

    def get_edge_axes(self):
        n1, n2, n3, n4 = self.nodes_ref

        g1 = n1.get_position()
        g2 = n2.get_position()
        g3 = n3.get_position()
        g4 = n4.get_position()
        g12 = (g1 + g2) / 2.
        g23 = (g2 + g3) / 2.
        g34 = (g3 + g4) / 2.
        g14 = (g1 + g4) / 2.
        x = g23 - g14
        yprime = g34 - g12
        normal = cross(x, yprime)
        y = cross(normal, x)
        return x, y

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid_ref.Thickness()

    def Normal(self):
        (n1, n2, n3, n4) = self.get_node_positions()
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
        n1, n2, n3, n4 = self.get_node_positions()
        area = 0.5 * norm(cross(n3-n1, n4-n2))
        centroid = (n1 + n2 + n3 + n4) / 4.
        return(area, centroid)

    def Centroid(self):
        n1, n2, n3, n4 = self.get_node_positions()
        centroid = (n1 + n2 + n3 + n4) / 4.
        return centroid

    def get_area(self):
        return self.Area()

    def Area(self):
        """
        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert
        where a and b are the quad's cross node point vectors"""
        (n1, n2, n3, n4) = self.get_node_positions()
        area = 0.5 * norm(cross(n3-n1, n4-n2))
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

    def _get_repr_defaults(self):
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

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref


class CSHEAR(QuadShell):
    type = 'CSHEAR'
    calculixType = 'S4'
    def __init__(self, eid, pid, nids, comment=''):
        QuadShell.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 4

    def validate(self):
        assert len(set(self.nodes)) == 4, 'nodes=%s\n%s' % (self.nodes, str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer_or_blank(card, 3, 'n1'),
                integer_or_blank(card, 4, 'n2'),
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4')]
        assert len(card) <= 7, 'len(CSHEAR card) = %i\ncard=%s' % (len(card), card)
        return CSHEAR(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = data[2:]
        return CSHEAR(eid, pid, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CSHEAR eid=%s' % self.eid
        self.nodes = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.pid = model.Property(self.Pid(), msg=msg)
        self.nodes_ref = self.nodes
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def Normal(self):
        (n1, n2, n3, n4) = self.get_node_positions()
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
        (n1, n2, n3, n4) = self.get_node_positions()
        a = n1 - n2
        b = n2 - n4
        area1 = 0.5 * norm(cross(a, b))

        a = n2 - n4
        b = n2 - n3
        area2 = 0.5 * norm(cross(a, b))

        area = area1 + area2
        centroid = (n1 + n2 + n3 + n4) / 4.
        return(area, centroid)

    def Centroid(self):
        (n1, n2, n3, n4) = self.get_node_positions()
        centroid = (n1 + n2 + n3 + n4) / 4.
        return centroid

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        edges = self.get_edge_ids()

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PSHEAR'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            if self.pid_ref.type in ['PSHEAR']:
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
        (n1, n2, n3, n4) = self.get_node_positions()
        a = n1 - n3
        b = n2 - n4
        area = 0.5 * norm(cross(a, b))
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
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=False)

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
        msg = self.comment + print_card_8(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg

    def G(self):
        return self.pid_ref.mid_ref.G()

    def Thickness(self):
        return self.pid_ref.t


class CQUAD4(QuadShell):
    """
    +--------+-------+-------+----+----+----+----+------------+---------+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |     8      |    9    |
    +========+=======+=======+=====+===+====+====+============+=========+
    | CQUAD4 |  EID  |  PID  | N1 | N2 | N3 | N4 | THETA/MCID | ZOFFSET |
    +--------+-------+-------+----+----+----+----+------------+---------+
    |        |       | TFLAG | T1 | T2 | T3 | T4 |            |         |
    +--------+-------+-------+----+----+----+----+------------+---------+
    """
    type = 'CQUAD4'
    aster_type = 'QUAD4 # CQUAD4'
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

    def __init__(self, eid, pid, nids, thetaMcid=0.0, zOffset=0.,
                 TFlag=0, T1=1.0, T2=1.0, T3=1.0, T4=1.0, comment=''):
        QuadShell.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        assert len(nids) == 4, nids
        self.prepare_node_ids(nids)
        self.zOffset = zOffset
        self.thetaMcid = thetaMcid
        self.TFlag = TFlag
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.T4 = T4

    def validate(self):
        assert len(set(self.nodes)) == 4, 'nodes=%s\n%s' % (self.nodes, str(self))

    @classmethod
    def add_op2_data(cls, data, comment=''):
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
        return CQUAD4(eid, pid, nids, thetaMcid, zOffset,
                      TFlag, T1, T2, T3, T4, comment=comment)

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4')]
        if len(card) > 6:
            thetaMcid = integer_double_or_blank(card, 7, 'thetaMcid', 0.0)
            zOffset = double_or_blank(card, 8, 'zOffset', 0.0)
            blank(card, 9, 'blank')
            TFlag = integer_or_blank(card, 10, 'TFlag', 0)
            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            T4 = double_or_blank(card, 14, 'T4')
            assert len(card) <= 15, 'len(CQUAD4 card) = %i\ncard=%s' % (len(card), card)
        else:
            thetaMcid = 0.0
            zOffset = 0.0
            TFlag = 0
            T1 = 1.0
            T2 = 1.0
            T3 = 1.0
            T4 = 1.0
        return CQUAD4(eid, pid, nids, thetaMcid, zOffset,
                      TFlag, T1, T2, T3, T4, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CQUAD4 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def material_coordinate_system(self, normal=None, xyz1234=None):
        if normal is None:
            normal = self.Normal() # k = kmat

        if xyz1234 is None:
            xyz1 = self.nodes_ref[0].get_position()
            xyz2 = self.nodes_ref[1].get_position()
            xyz3 = self.nodes_ref[2].get_position()
            xyz4 = self.nodes_ref[3].get_position()
            #centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.
            #centroid = self.Centroid()
        else:
            #centroid = xyz1234.sum(axis=1)
            #assert len(centroid) == 3, centroid
            xyz1 = xyz1234[:, 0]
            xyz2 = xyz1234[:, 1]
            xyz3 = xyz1234[:, 2]
            xyz4 = xyz1234[:, 3]
        centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.

        if self.thetaMcid is None:
            raise NotImplementedError('thetaMcid=%r' % self.thetaMcid)
        if isinstance(self.thetaMcid, integer_types):
            i = self.thetaMcid_ref.i
            jmat = np.cross(normal, i) # k x i
            jmat /= np.linalg.norm(jmat)
            imat = np.cross(jmat, normal)
        elif isinstance(self.thetaMcid, float):
            raise NotImplementedError('thetaMcid=%r' % self.thetaMcid)
        return centroid, imat, jmat, normal

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        edges = self.get_edge_ids()
        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            if not self.pid_ref.type in ['PLPLANE']:
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
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=False)

    def writeAs_ctria3(self, newID):
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
                       [self.thetaMcid, self.zOffset, None, self.TFlag, self.T1, self.T2,
                        self.T3, self.T4])
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3, T4) = self._get_repr_defaults()
        list_fields = (['CQUAD4', self.eid, self.Pid()] + self.node_ids +
                       [thetaMcid, zOffset, None, TFlag, T1, T2, T3, T4])
        return list_fields

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids

        row2_data = [self.thetaMcid, self.zOffset,
                     self.TFlag, self.T1, self.T2, self.T3, self.T4]
        if row2_data == [0.0, 0.0, 0, 1.0, 1.0, 1.0, 1.0]:
            data = [self.eid, self.Pid()] + nodes
            msg = ('CQUAD4  %8i%8i%8i%8i%8i%8i\n' % tuple(data))
            return self.comment + msg
        else:
            thetaMcid = set_blank_if_default(self.thetaMcid, 0.0)
            zOffset = set_blank_if_default(self.zOffset, 0.0)
            TFlag = set_blank_if_default(self.TFlag, 0)
            T1 = set_blank_if_default(self.T1, 1.0)
            T2 = set_blank_if_default(self.T2, 1.0)
            T3 = set_blank_if_default(self.T3, 1.0)
            T4 = set_blank_if_default(self.T4, 1.0)

            row2_data = [thetaMcid, zOffset,
                         TFlag, T1, T2, T3, T4]
            row2 = [print_field_8(field) for field in row2_data]
            data = [self.eid, self.Pid()] + nodes + row2
            msg = ('CQUAD4  %8i%8i%8i%8i%8i%8i%8s%8s\n'
                   '                %8s%8s%8s%8s%8s\n' % tuple(data))
            return self.comment + msg.rstrip() + '\n'

    #def write_card(self, size=8, is_double=False):
        #card = wipe_empty_fields(self.repr_fields())
        #if size == 8 or len(card) == 7: # to last node
            #msg = self.comment + print_card_8(card)
        #else:
            #msg = self.comment + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        #return msg


class CQUADR(QuadShell):
    type = 'CQUADR'
    #calculixType = 'CAX8'

    def __init__(self, eid, pid, nids, thetaMcid, zOffset, TFlag,
                 T1, T2, T3, T4, comment=''):
        QuadShell.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.thetaMcid = thetaMcid
        self.zOffset = zOffset
        self.TFlag = TFlag
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.T4 = T4
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 4, 'CQUADR'

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 'n1'),
                integer_or_blank(card, 4, 'n2'),
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4')]

        thetaMcid = integer_double_or_blank(card, 7, 'thetaMcid', 0.0)
        zOffset = double_or_blank(card, 8, 'zOffset', 0.0)

        TFlag = integer_or_blank(card, 10, 'TFlag', 0)
        T1 = double_or_blank(card, 11, 'T1')
        T2 = double_or_blank(card, 12, 'T2')
        T3 = double_or_blank(card, 13, 'T3')
        T4 = double_or_blank(card, 14, 'T4')
        assert len(card) <= 15, 'len(CQUADR card) = %i\ncard=%s' % (len(card), card)
        return CQUADR(eid, pid, nids, thetaMcid, zOffset,
                      TFlag, T1, T2, T3, T4, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
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
        return CQUADR(eid, pid, nids, thetaMcid, zOffset,
                      TFlag, T1, T2, T3, T4, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CQUADR eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PSHELL', 'PCOMP', 'PCOMPG'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
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
        return self.pid_ref.Thickness()

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
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = (['CQUADR', self.eid, self.Pid()] + self.node_ids +
                       [self.thetaMcid, self.zOffset, None, self.TFlag, self.T1,
                        self.T2, self.T3, self.T4])
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3, T4) = self._get_repr_defaults()
        list_fields = (['CQUADR', self.eid, self.Pid()] + self.node_ids +
                       [thetaMcid, zOffset, None, TFlag, T1, T2, T3, T4])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8 or len(card) == 7: # to last node
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg


class CQUAD(QuadShell):
    type = 'CQUAD'

    def __init__(self, eid, pid, nids, comment=''):
        QuadShell.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 9

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5'),
                integer_or_blank(card, 8, 'n6'),
                integer_or_blank(card, 9, 'n7'),
                integer_or_blank(card, 10, 'n8'),
                integer_or_blank(card, 11, 'n9')]
        assert len(card) <= 12, 'len(CQUAD card) = %i\ncard=%s' % (len(card), card)
        return CQUAD(eid, pid, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CQUAD eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid_ref.Thickness()

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
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

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
        return self.comment + msg.rstrip() + '\n'

    #def write_card(self, size=8, is_double=False):
        #card = self.repr_fields()
        #msg = self.comment + print_card_8(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        #return msg


class CQUAD8(QuadShell):
    type = 'CQUAD8'
    aster_type = 'QUAD8'

    def __init__(self, eid, pid, nids, T1, T2, T3, T4, thetaMcid, zOffset, TFlag,
                 comment=''):
        QuadShell.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.T4 = T4
        self.TFlag = TFlag
        self.thetaMcid = thetaMcid
        self.zOffset = zOffset
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 8

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5', 0),
                integer_or_blank(card, 8, 'n6', 0),
                integer_or_blank(card, 9, 'n7', 0),
                integer_or_blank(card, 10, 'n8', 0)]
        if len(card) > 11:
            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            T4 = double_or_blank(card, 14, 'T4')
            thetaMcid = integer_double_or_blank(card, 15, 'thetaMcid', 0.0)
            zOffset = double_or_blank(card, 16, 'zOffset', 0.0)
            TFlag = integer_or_blank(card, 17, 'TFlag', 0)
            assert len(card) <= 18, 'len(CQUAD4 card) = %i\ncard=%s' % (len(card), card)
        else:
            thetaMcid = 0.0
            zOffset = 0.0
            T1 = 1.0
            T2 = 1.0
            T3 = 1.0
            T4 = 1.0
            TFlag = 0
        return CQUAD8(eid, pid, nids, T1, T2, T3, T4, thetaMcid, zOffset, TFlag,
                      comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        #print "CQUAD8 = ",data
        #(6401,
        #6400,
        #6401, 6402, 6405, 6403, 0, 0, 6404, 0,
        #-1.0, -1.0, -1.0, -1.0,
        #0.0, 0)
        eid = data[0]
        pid = data[1]
        nids = data[2:10]
        T1 = data[10]
        T2 = data[11]
        T3 = data[12]
        T4 = data[13]
        thetaMcid = data[14]
        zOffset = data[14]
        TFlag = data[15]
        return CQUAD8(eid, pid, nids, T1, T2, T3, T4, thetaMcid, zOffset, TFlag,
                      comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CQUAD8 eid=%s' % self.eid
        self.nodes = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.Pid(), msg=msg)
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        edges = self.get_edge_ids()

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert isinstance(nid, integer_types) or nid is None, 'nid%i is not an integer/None; nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
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
        return self.pid_ref.Thickness()

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
        (n1, n2, n3, n4) = self.get_node_positions()[:4]
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
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.get_node_positions()
        a = n1 - n2
        b = n2 - n4
        area1 = 0.5 * norm(cross(a, b))
        c1 = (n1 + n2 + n4) / 3.

        a = n2 - n4
        b = n2 - n3
        area2 = 0.5 * norm(cross(a, b))
        c2 = (n2 + n3 + n4) / 3.

        area = area1 + area2
        centroid = (c1 * area1 + c2 * area2) / area
        return(area, centroid)

    def Area(self):
        r"""
        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert
        where a and b are the quad's cross node point vectors"""
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.get_node_positions()
        a = n1 - n3
        b = n2 - n4
        area = 0.5 * norm(cross(a, b))
        return area

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CQUAD8', self.eid, self.Pid()] + self.node_ids + [
            self.T1, self.T2, self.T3, self.T4, self.thetaMcid, self.zOffset,
            self.TFlag]
        return list_fields

    def repr_fields(self):
        (thetaMcid, zOffset, TFlag, T1, T2, T3, T4) = self._get_repr_defaults()
        list_fields = (['CQUAD8', self.eid, self.Pid()] + self.node_ids + [
            T1, T2, T3, T4, thetaMcid, zOffset, TFlag])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8 or len(card) == 11: # to last node
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class CQUADX(QuadShell):
    type = 'CQUADX'
    calculixType = 'CAX8'

    def __init__(self, eid, pid, nids, comment=''):
        QuadShell.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 9

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer_or_blank(card, 3, 'n1'),
            integer_or_blank(card, 4, 'n2'),
            integer_or_blank(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4'),
            integer_or_blank(card, 7, 'n5'),
            integer_or_blank(card, 8, 'n6'),
            integer_or_blank(card, 9, 'n7'),
            integer_or_blank(card, 10, 'n8'),
            integer_or_blank(card, 11, 'n9')
        ]
        assert len(card) <= 12, 'len(CQUADX card) = %i\ncard=%s' % (len(card), card)
        return CQUADX(eid, pid, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CQUADX eid=%s' % self.eid
        self.nodes = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.Pid(), msg=msg)
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def Thickness(self):
        """
        Returns the thickness
        """
        return self.pid_ref.Thickness()

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
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced
        """
        pass

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
        return self.comment + msg.rstrip() + '\n'
