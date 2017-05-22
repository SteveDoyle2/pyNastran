## pylint: disable=C0103,R0902,R0904,R0914,C0302
"""
All axisymmetric shell elements are defined in this file.  This includes:
 * CTRAX3
 * CTRAX6
 * CTRIAX
 * CTRIAX6
 * CQUADX
 * CQUADX4
 * CQUADX8

All tris are TriShell, ShellElement, and Element objects.
All quads are QuadShell, ShellElement, and Element objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import range

from numpy import cross
from numpy.linalg import norm

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import (
    set_blank_if_default, set_default_if_blank,
    print_card_8, print_field_8)
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, integer_double_or_blank)
from pyNastran.bdf.cards.utils import wipe_empty_fields
from pyNastran.bdf.cards.elements.shell import TriShell, QuadShell, _triangle_area_centroid_normal, _normal
from pyNastran.bdf.cards.base_card import Element

__all__ = ['CTRAX3', 'CTRAX6', 'CTRIAX', 'CTRIAX6',
           'CQUADX', 'CQUADX4', 'CQUADX8']


class AxisymmetricElement(Element):
    def __init__(self):
        Element.__init__(self)


class AxisymmetricTri(AxisymmetricElement):
    def __init__(self):
        AxisymmetricElement.__init__(self)

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

    def Centroid(self):
        r"""
        Get the centroid.

        .. math::
          CG = \frac{1}{3} (n_0+n_1+n_2)
        """
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes[:3])
        centroid = (n1 + n2 + n3) / 3.
        return centroid

    def Mass(self):
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes[:3])
        return 0.

class CTRAX3(AxisymmetricTri):
    """
    +--------+------------+-------+----+----+----+-------+
    |   1    |     2      |   3   |  4 |  5 |  6 |   7   |
    +========+============+=======+====+====+====+=======+
    | CTRAX3 |    EID     |  PID  | N1 | N2 | N3 | THETA |
    +--------+------------+-------+----+----+----+-------+

    CTRAX3 is NX only!
    """
    type = 'CTRAX3'
    def __init__(self, eid, pid, nids, theta=0., comment=''):
        AxisymmetricTri.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.theta = theta
        self.nodes = nids
        assert len(nids) == 3, 'error on CTRAX3'

    def validate(self):
        self.validate_node_ids(allow_empty_nodes=True)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CTRAX3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')

        nids = [
            integer_or_blank(card, 3, 'n1'),
            integer_or_blank(card, 4, 'n2'),
            integer_or_blank(card, 5, 'n3'),
            ]
        theta = integer_double_or_blank(card, 6, 'theta', 0.0)
        assert len(card) <= 7, 'len(CTRAX3 card) = %i\ncard=%s' % (len(card), card)
        return CTRAX3(eid, pid, nids, theta=theta, comment=comment)

    def _verify(self, xref=True):
        eid = self.eid
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
            assert self.pid_ref.type in ['PSOLID', 'PLSOLID'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            mass = self.Mass()
            assert isinstance(mass, float), 'mass=%r' % mass
            c = self.Centroid()
            for i in range(3):
                assert isinstance(c[i], float)

    def flipNormal(self):
        pass

    def Mass(self):
        return 0.

    def AreaCentroidNormal(self):
        """
        Returns area, centroid, normal as it's more efficient to do them
        together
        """
        n1, n2, n3 = self.get_node_positions()
        return _triangle_area_centroid_normal([n1, n2, n3], self)

    def Area(self):
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} \lvert (n_1-n_2) \times (n_1-n_3) \rvert"""
        (n1, n2, n3) = self.get_node_positions()
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
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=False, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)
        self.nodes_ref = self.nodes
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=False)

    def raw_fields(self):
        list_fields = ['CTRAX3', self.eid, self.Pid()] + self.node_ids + [self.theta]
        return list_fields

    def repr_fields(self):
        theta = set_blank_if_default(self.theta, 0.0)
        nodeIDs = self.node_ids
        list_fields = ['CTRAX3', self.eid, self.Pid()] + nodeIDs + [theta]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 8: # to last node
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        return msg


class CTRAX6(AxisymmetricTri):
    """
    +--------+-------+-------+----+----+----+----+----+-----+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |  8 |  9  |
    +========+=======+=======+====+====+====+====+====+=====+
    | CTRAX6 |  EID  |  PID  | N1 | N2 | N3 | N4 | N5 | N6  |
    +--------+-------+-------+----+----+----+----+----+-----+
    |        | THETA |       |    |    |    |    |    |     |
    +--------+-------+-------+----+----+----+----+----+-----+

    Theta/Mcid is NX only!
    """
    type = 'CTRAX6'
    def __init__(self, eid, pid, nids, theta=0., comment=''):
        AxisymmetricTri.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.theta = theta
        self.nodes = nids
        assert len(nids) == 6, 'error on CTRAX6'

    def validate(self):
        self.validate_node_ids(allow_empty_nodes=True)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CTRAX6 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4'),
            integer_or_blank(card, 7, 'n5'),
            integer_or_blank(card, 8, 'n6'),
            ]
        theta = integer_double_or_blank(card, 9, 'theta', 0.0)
        assert len(card) <= 10, 'len(CTRAX6 card) = %i\ncard=%s' % (len(card), card)
        return CTRAX6(eid, pid, nids, theta=theta, comment=comment)

    def _verify(self, xref=True):
        eid = self.eid
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
            assert self.pid_ref.type in ['PSOLID', 'PLSOLID'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            mass = self.Mass()
            assert isinstance(mass, float), 'mass=%r' % mass
            c = self.Centroid()
            for i in range(3):
                assert isinstance(c[i], float)

    def flipNormal(self):
        pass

    def AreaCentroidNormal(self):
        """
        Returns area, centroid, normal as it's more efficient to do them
        together
        """
        (n1, n2, n3) = self.get_node_positions(nodes=self.nodes[:3])
        return _triangle_area_centroid_normal([n1, n2, n3], self)

    def Area(self):
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} \lvert (n_1-n_2) \times (n_1-n_3) \rvert"""
        (n1, n2, n3) = self.get_node_positions(nodes=self.nodes[:3])
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
        msg = ' which is required by CTRAX6 eid=%s' % self.eid
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)
        self.nodes_ref = self.nodes
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CTRAX6', self.eid, self.Pid()] + self.node_ids + [self.theta]
        return list_fields

    def repr_fields(self):
        theta = set_blank_if_default(self.theta, 0.0)
        nodeIDs = self.node_ids
        list_fields = ['CTRAX6', self.eid, self.Pid()] + nodeIDs + [theta]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 8: # to last node
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        return msg


class CTRIAX(TriShell):
    """
    +--------+------------+-------+----+----+----+----+----+-----+
    |   1    |     2      |   3   |  4 |  5 |  6 | 7  |  8 |  9  |
    +========+============+=======+====+====+====+====+====+=====+
    | CTRIAX |    EID     |  PID  | N1 | N2 | N3 | N4 | N5 | N6  |
    +--------+------------+-------+----+----+----+----+----+-----+
    |        | THETA/MCID |       |    |    |    |    |    |     |
    +--------+------------+-------+----+----+----+----+----+-----+

    Theta/Mcid is MSC only!
    """
    type = 'CTRIAX'
    def __init__(self, eid, pid, nids, theta_mcid=0., comment=''):
        TriShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID of a PLPLANE or PAXSYMH entry
        self.pid = pid
        self.theta_mcid = theta_mcid
        self.nodes = nids
        assert len(nids) == 6, 'error on CTRIAX'

    def validate(self):
        self.validate_node_ids(allow_empty_nodes=True)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CTRIAX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        theta_mcid = integer_double_or_blank(card, 9, 'theta_mcsid', 0.0)
        assert len(card) <= 10, 'len(CTRIAX card) = %i\ncard=%s' % (len(card), card)
        return CTRIAX(eid, pid, nids, theta_mcid=theta_mcid, comment=comment)

    def _verify(self, xref=True):
        eid = self.eid
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
        (n1, n2, n3) = self.get_node_positions(nodes=self.nodes[:3])
        return _triangle_area_centroid_normal([n1, n2, n3], self)

    def Area(self):
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} \lvert (n_1-n_2) \times (n_1-n_3) \rvert"""
        (n1, n2, n3) = self.get_node_positions(nodes=self.nodes[:3])
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

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CTRIAX', self.eid, self.Pid()] + self.node_ids + [self.theta_mcid]
        return list_fields

    def repr_fields(self):
        theta_mcid = set_blank_if_default(self.theta_mcid, 0.0)
        nodeIDs = self.node_ids
        list_fields = ['CTRIAX', self.eid, self.Pid()] + nodeIDs + [theta_mcid]
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
    +---------+-------+-------+----+----+----+----+----+-----+
    |    1    |   2   |   3   |  4 |  5 |  6 |  7 |  8 |  9  |
    +=========+=======+=======+=====+===+====+====+====+=====+
    | CTRIAX6 |  EID  |  MID  | N1 | N2 | N3 | G4 | G5 | G6  |
    +---------+-------+-------+----+----+----+----+----+-----+
    |         | THETA |       |    |    |    |    |    |     |
    +---------+-------+-------+----+----+----+----+----+-----+


    NX/MSC : Nodes are defined in a non-standard way::

           5
          / \
         6   4
       /       \
      1----2----3
    """
    type = 'CTRIAX6'
    pid = -53 # uses element type from OP2
    def __init__(self, eid, mid, nids, theta=0., comment=''):
        TriShell.__init__(self)
        if comment:
            self.comment = comment
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
        """
        Adds a CTRIAX6 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        return CTRIAX6(eid, mid, nids, theta=theta, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid, mid, n1, n2, n3, n4, n5, n6, theta, undef1, undef2 = data
        nids = [n1, n2, n3, n4, n5, n6]
        return CTRIAX6(eid, mid, nids, theta=theta, comment=comment)

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
        eid = self.eid
        nids = self.node_ids
        edges = self.get_edge_ids()

        assert self.pid == -53, 'pid = %s' % self.pid
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
        #raise AttributeError("CTRIAX6 doesn't have a Property")
        return self.pid

    def AreaCentroidNormal(self):
        """
        Returns area, centroid, normal as it's more efficient to do them
        together
        """
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()
        return _triangle_area_centroid_normal([n1, n3, n5], self)

    def Area(self):
        r"""
        Get the normal vector.

        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_1-n_5) \rvert"""
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()
        a = n1 - n3
        b = n1 - n5
        area = 0.5 * norm(cross(a, b))
        return area

    def Centroid(self):
        r"""
        Get the centroid.

        .. math::
          CG = \frac{1}{3} (n_0+n_1+n_2)
        """
        n1, n2, n3, n4, n5, n6 = self.get_node_positions()
        centroid = (n1 + n3 + n5) / 3.
        return centroid

    def Nsm(self):
        raise AttributeError('CTRIAX6 does not have a non-structural mass')

    def MassPerArea(self):
        raise AttributeError('CTRIAX6 does not have a MassPerArea')

    def Mass(self):
        #raise NotImplementedError('CTRIAX6 does not have a Mass method yet')
        return 0.

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
        return msg


class CQUADX(QuadShell):
    """
    Defines an axisymmetric quadrilateral element with up to nine grid
    points for use in fully nonlinear (i.e., large strain and large
    rotations) analysis or a linear harmonic or rotordynamic analysis.
    The element has between four and eight grid points

    +--------+-------+-------+----+------------+----+----+-----+-----+
    |   1    |   2   |   3   |  4 |      5     |  6 | 7  |  8  |  9  |
    +========+=======+=======+====+============+====+====+=====+=====+
    | CQUADX |  EID  |  PID  | N1 |     N2     | N3 | N4 |  G5 | G6  |
    +--------+-------+-------+----+------------+----+----+-----+-----+
    |        |  G7   |  G8   | G9 | THETA/MCID |    |    |     |     |
    +--------+-------+-------+----+------------+----+----+-----+-----+

    Theta/Mcid is MSC only!
    """
    type = 'CQUADX'
    def __init__(self, eid, pid, nids, theta_mcid=0., comment=''):
        QuadShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.theta_mcid = theta_mcid
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 9

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CQUADX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer(card, 6, 'n4'),
            integer_or_blank(card, 7, 'n5'),
            integer_or_blank(card, 8, 'n6'),
            integer_or_blank(card, 9, 'n7'),
            integer_or_blank(card, 10, 'n8'),
            integer_or_blank(card, 11, 'n9')
        ]
        theta_mcid = integer_double_or_blank(card, 12, 'theta/mcid', 0.)
        assert len(card) <= 12, 'len(CQUADX card) = %i\ncard=%s' % (len(card), card)
        return CQUADX(eid, pid, nids, theta_mcid=theta_mcid, comment=comment)

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
        theta_mcid = set_blank_if_default(self.theta_mcid, 0.0)
        row2 = ['        ' if node is None else '%8i' % node for node in nodes[4:]
                 ] + [print_field_8(theta_mcid)]
        msg = ('CQUADX  %8i%8i%8i%8i%8i%8i%8s%8s\n'
               '        %8s%8s%8s%s' % tuple(data + row2))
        return self.comment + msg.rstrip() + '\n'


class CQUADX4(QuadShell):
    """
    Defines an isoparametric and axisymmetric quadrilateral cross-section
    ring element for use in linear and fully nonlinear (i.e., large strain
    and large rotations) hyperelastic analysis.

    +---------+-------+-------+----+----+----+----+-------+
    |    1    |   2   |   3   |  4 |  5 |  6 | 7  |   8   |
    +=========+=======+=======+====+====+====+====+=======+
    | CQUADX4 |  EID  |  PID  | N1 | N2 | N3 | N4 | THETA |
    +---------+-------+-------+----+----+----+----+-------+

    CQUADX4 is an NX card only!
    """
    type = 'CQUADX4'

    def __init__(self, eid, pid, nids, theta=0., comment=''):
        QuadShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.theta = theta
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 4

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CQUADX4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer(card, 6, 'n4'),
        ]
        theta = integer_double_or_blank(card, 7, 'theta', 0.)
        assert len(card) <= 8, 'len(CQUADX4 card) = %i\ncard=%s' % (len(card), card)
        return CQUADX4(eid, pid, nids, theta=theta, comment=comment)

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

    def flipNormal(self):
        r"""
        ::

          1--5--2       1--8--4
          |     |  -->  |     |
          8  9  6       5  9  7
          |     |       |     |
          4--7--3       2--6--3
        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]

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
        list_fields = ['CQUADX4', self.eid, self.Pid()] + self.node_ids + [self.theta]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        data = ['CQUADX4', self.eid, self.Pid()] + nodes + [self.theta]
        return self.comment + print_card_8(data)

class CQUADX8(QuadShell):
    """
    Defines an isoparametric and axisymmetric quadrilateral cross-section
    ring element with midside nodes for use in linear and fully nonlinear
    (i.e., large strain and large rotations) hyperelastic analysis.

    +---------+-------+-------+-------+----+----+----+-----+-----+
    |    1    |   2   |   3   |   4   |  5 |  6 | 7  |  8  |  9  |
    +=========+=======+=======+=======+====+====+====+=====+=====+
    | CQUADX8 |  EID  |  PID  |  N1   | N2 | N3 | N4 |  G5 | G6  |
    +---------+-------+-------+-------+----+----+----+-----+-----+
    |         |  G7   |  G8   | THETA |    |    |    |     |     |
    +---------+-------+-------+-------+----+----+----+-----+-----+

    CQUADX8 is an NX card only!
    """
    type = 'CQUADX8'
    def __init__(self, eid, pid, nids, theta=0., comment=''):
        QuadShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.theta = theta
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 8

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CQUADX8 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer(card, 6, 'n4'),
            integer_or_blank(card, 7, 'n5'),
            integer_or_blank(card, 8, 'n6'),
            integer_or_blank(card, 9, 'n7'),
            integer_or_blank(card, 10, 'n8'),
        ]
        theta = integer_double_or_blank(card, 11, 'theta', 0.)
        assert len(card) <= 12, 'len(CQUADX8 card) = %i\ncard=%s' % (len(card), card)
        return CQUADX8(eid, pid, nids, theta=theta, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CQUADX8 eid=%s' % self.eid
        self.nodes = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.Pid(), msg=msg)
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def Normal(self):
        (n1, n2, n3, n4) = self.get_node_positions()[:4]
        return _normal(n1 - n3, n2 - n4)

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
        list_fields = ['CQUADX8', self.eid, self.Pid()] + self.node_ids + [self.theta]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes[:6]
        theta = set_blank_if_default(self.theta, 0.0)
        row2 = ['        ' if node is None else '%8i' % node for node in nodes[6:]
                 ] + [print_field_8(theta)]
        msg = ('CQUADX8 %8i%8i%8i%8i%8i%8i%8s%8s\n'
               '        %8s%8s%s' % tuple(data + row2))
        return self.comment + msg.rstrip() + '\n'
