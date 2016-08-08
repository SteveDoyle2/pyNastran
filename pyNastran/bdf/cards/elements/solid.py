# pylint: disable=R0902,R0904,R0914
"""
All solid elements are defined in this file.  This includes:

 * CHEXA8
 * CHEXA20
 * CPENTA6
 * CPENTA15
 * CTETRA4
 * CTETRA10

All solid elements are SolidElement and Element objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import range
from numpy import dot, cross
from numpy.linalg import norm

from pyNastran.bdf.cards.elements.elements import Element
from pyNastran.utils.mathematics import Area
from pyNastran.bdf.bdf_interface.assign_type import integer, integer_or_blank


_chexa_mapper = {
    (7, 5) : (7, 6, 5, 4),
    (5, 7) : (7, 6, 5, 4),
    (6, 4) : (7, 6, 5, 4),
    (4, 6) : (7, 6, 5, 4),

    (0, 2) : (0, 1, 2, 3),
    (2, 0) : (0, 1, 2, 3),
    (1, 3) : (0, 1, 2, 3),
    (3, 1) : (0, 1, 2, 3),

    (0, 7) : (0, 3, 7, 4),
    (7, 0) : (0, 3, 7, 4),
    (3, 4) : (0, 3, 7, 4),
    (4, 3) : (0, 3, 7, 4),

    (5, 2) : (5, 6, 2, 1),
    (2, 5) : (5, 6, 2, 1),
    (6, 1) : (5, 6, 2, 1),
    (1, 6) : (5, 6, 2, 1),

    (4, 1) : (4, 5, 1, 0),
    (1, 4) : (4, 5, 1, 0),
    (5, 0) : (4, 5, 1, 0),
    (0, 5) : (4, 5, 1, 0),

    (2, 7) : (2, 6, 7, 3),
    (7, 2) : (2, 6, 7, 3),
    (6, 3) : (2, 6, 7, 3),
    (3, 6) : (2, 6, 7, 3),
}
_chexa_faces = (
    (7, 6, 5, 4),
    (0, 1, 2, 3),
    (0, 3, 7, 4),
    (5, 6, 2, 1),
    (4, 5, 1, 0),
    (2, 6, 7, 3),
)

def volume4(n1, n2, n3, n4):
    r"""
    Gets the volume, :math:`V`, of the tetrahedron.

    .. math:: V = \frac{(a-d) \cdot \left( (b-d) \times (c-d) \right) }{6}
    """
    volume = -dot(n1 - n4, cross(n2 - n4, n3 - n4)) / 6.
    return volume


def area_centroid(n1, n2, n3, n4):
    """
    Gets the area, :math:`A`, and centroid of a quad.::

      1-----2
      |   / |
      | /   |
      4-----3
    """
    area1 = 0.5 * norm(cross(n1 - n2, n2 - n4))
    c1 = (n1 + n2 + n4) / 3.

    area2 = 0.5 * norm(cross(n2 - n4, n2 - n3))
    c2 = (n2 + n3 + n4) / 3.

    area = area1 + area2
    try:
        centroid = (c1 * area1 + c2 * area2) / area
    except FloatingPointError:
        msg = '\nc1=%r\narea1=%r\n' % (c1, area1)
        msg += 'c2=%r\narea2=%r' % (c2, area2)
        raise FloatingPointError(msg)
    return(area, centroid)


class SolidElement(Element):
    _field_map = {1: 'nid', 2:'pid'}

    def __init__(self):
        Element.__init__(self)

    def _update_field_helper(self, n, value):
        if n - 3 < len(self.nodes):
            self.nodes[n - 3] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def cross_reference(self, model):
        raise NotImplementedError('Element type=%r must implement cross_reference')

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref

    def Eid(self):
        return self.eid

    def E(self):
        return self.pid_ref.mid_ref.E()

    def G(self):
        return self.pid_ref.mid_ref.G()

    def Nu(self):
        return self.pid_ref.mid_ref.Nu()

    def Volume(self):
        """
        Base volume method that should be overwritten
        """
        pass

    def Mass(self):
        """
        Calculates the mass of the solid element
        Mass = Rho * Volume
        """
        return self.Rho() * self.Volume()

    def Mid(self):
        """
        Returns the material ID as an integer
        """
        return self.pid_ref.Mid()

    def Rho(self):
        """
        Returns the density
        """
        try:
            return self.pid_ref.Rho()
        except AttributeError:
            print("self.pid = %s" % (self.pid))
            #print("self.pid_ref.mid_ref = %s" % (str(self.pid_ref.mid_ref)))
            raise

    def _is_same_card(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid, self.Pid()] + self.nodes
        fields2 = [elem.eid, elem.Pid()] + elem.nodes
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self._is_same_fields(fields1, fields2)

    def raw_fields(self):
        list_fields = [self.type, self.eid, self.Pid()] + self.node_ids
        return list_fields

class CHEXA(object):
    def __init__(self, eid, pid, nids, comment=''):
        if len(nids) == 8:
            return CHEXA8(eid, pid, nids, comment=comment)
        elif len(nids) == 20:
            return CHEXA20(eid, pid, nids, comment=comment)
        else:
            msg = 'len(nids)=%s; expected 8 or 20; nids=%s' % (len(nids), nids)
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment):
        if card.nfields == 11:
            return CHEXA8.add_card(card, comment=comment)
        else:
            return CHEXA20.add_card(card, comment=comment)

class CTETRA(object):
    def __init__(self, eid, pid, nids, comment=''):
        if len(nids) == 4:
            return CTETRA4(eid, pid, nids, comment=comment)
        elif len(nids) == 10:
            return CTETRA10(eid, pid, nids, comment=comment)
        else:
            msg = 'len(nids)=%s; expected 8 or 20; nids=%s' % (len(nids), nids)
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment):
        if card.nfields == 7:
            return CTETRA4.add_card(card, comment=comment)
        else:
            return CTETRA10.add_card(card, comment=comment)

class CPENTA(object):
    def __init__(self, eid, pid, nids, comment=''):
        if len(nids) == 6:
            return CPENTA6(eid, pid, nids, comment=comment)
        elif len(nids) == 15:
            return CPENTA15(eid, pid, nids, comment=comment)
        else:
            msg = 'len(nids)=%s; expected 8 or 20; nids=%s' % (len(nids), nids)
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment):
        if card.nfields == 9:
            return CPENTA6.add_card(card, comment=comment)
        else:
            return CPENTA15.add_card(card, comment=comment)

class CPYRAM(object):
    def __init__(self, eid, pid, nids, comment=''):
        if len(nids) == 5:
            return CPYRAM5(eid, pid, nids, comment=comment)
        elif len(nids) == 13:
            return CPYRAM13(eid, pid, nids, comment=comment)
        else:
            msg = 'len(nids)=%s; expected 8 or 20; nids=%s' % (len(nids), nids)
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment):
        if card.nfields == 8:
            return CPYRAM5.add_card(card, comment=comment)
        else:
            return CPYRAM13.add_card(card, comment=comment)


class CHEXA8(SolidElement):
    """
    +-------+-----+-----+----+----+----+----+----+----+
    |   1   |  2  |  3  |  4 |  5 |  6 |  7 |  8 |  9 |
    +=======+=====+=====+====+====+====+====+====+====+
    | CHEXA | EID | PID | G1 | G2 | G3 | G4 | G5 | G6 |
    +-------+-----+-----+----+----+----+----+----+----+
    |       | G7  | G8  |    |    |    |    |    |    |
    +-------+-----+-----+----+----+----+----+----+----+
    """
    type = 'CHEXA'
    aster_type = 'HEXA8'
    calculixType = 'C3D8'

    def write_card(self, size=8, is_double=False):
        data = [self.eid, self.Pid()] + self.node_ids
        msg = ('CHEXA   %8i%8i%8i%8i%8i%8i%8i%8i\n'
               '        %8i%8i\n' % tuple(data))
        return self.comment + msg

    def write_card_16(self, is_double=False):
        data = [self.eid, self.Pid()] + self.node_ids
        msg = ('CHEXA*  %16i%16i%16i%16i\n'
               '*       %16i%16i%16i%16i\n'
               '*       %16i%16i\n' % tuple(data))
        return self.comment + msg

    def __init__(self, eid, pid, nids, comment=''):
        SolidElement.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 8

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'),
            integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'),
            integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'),
            integer(card, 8, 'nid6'),
            integer(card, 9, 'nid7'),
            integer(card, 10, 'nid8')
        ]
        assert len(card) == 11, 'len(CHEXA8 card) = %i\ncard=%s' % (len(card), card)
        return CHEXA8(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = data[2:]
        assert len(data) == 10, 'len(data)=%s data=%s' % (len(data), data)
        return CHEXA8(eid, pid, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=False, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(self.node_ids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
        if xref:
            centroid = self.Centroid()
            volume = self.Volume()
            assert isinstance(volume, float)
            for i in range(3):
                assert isinstance(centroid[i], float)

    def Centroid(self):
        """
        Averages the centroids at the two faces
        """
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.get_node_positions()
        c1 = area_centroid(n1, n2, n3, n4)[1]
        c2 = area_centroid(n5, n6, n7, n8)[1]
        centroid = (c1 + c2) / 2.
        return centroid

    def Volume(self):
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.get_node_positions()
        (area1, c1) = area_centroid(n1, n2, n3, n4)
        (area2, c2) = area_centroid(n5, n6, n7, n8)
        volume = (area1 + area2) / 2. * norm(c1 - c2)
        return abs(volume)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=False)

    def get_face(self, nid_opposite, nid):
        nids = self.node_ids[:8]
        return chexa_face(nid_opposite, nid, nids)

    def getFaceAreaCentroidNormal(self, nid, nid_opposite):
        """
        Parameters
        ----------
        nid : int
            G1 - a grid point on the corner of a face
        nid_opposite : int
            G3 - the grid point diagonally opposite of G1
        """
        nids = self.node_ids[:8]
        return chexa_face_area_centroid_normal(nid, nid_opposite, nids, self.nodes[:8])

    def get_edge_ids(self):
        """
        Return the edge IDs
        # top (5-6-7-8)
        # btm (1-2-3-4)
        # left (1-2-3-4)
        # right (5-6-7-8)
        # front (1-5-8-4)
        # back (2-6-7-3)
        """
        node_ids = self.node_ids
        return [
            # btm (1-2-3-4)
            tuple(sorted([node_ids[0], node_ids[1]])),
            tuple(sorted([node_ids[1], node_ids[2]])),
            tuple(sorted([node_ids[2], node_ids[3]])),
            tuple(sorted([node_ids[3], node_ids[0]])),

            # top (5-6-7-8)
            tuple(sorted([node_ids[4], node_ids[5]])),
            tuple(sorted([node_ids[5], node_ids[6]])),
            tuple(sorted([node_ids[6], node_ids[7]])),
            tuple(sorted([node_ids[7], node_ids[4]])),

            # up - (4-8, 3-7, 1-5, 2-6)
            tuple(sorted([node_ids[0], node_ids[4]])),
            tuple(sorted([node_ids[1], node_ids[5]])),
            tuple(sorted([node_ids[2], node_ids[6]])),
            tuple(sorted([node_ids[3], node_ids[7]])),
        ]


class CIHEX1(CHEXA8):
    type = 'CIHEX1'

    def write_card(self, size=8, is_double=False):
        data = [self.eid, self.Pid()] + self.node_ids
        msg = ('CIHEX1  %8i%8i%8i%8i%8i%8i%8i%8i\n'
               '        %8i%8i\n' % tuple(data))
        return self.comment + msg

    def write_card_16(self, is_double=False):
        data = [self.eid, self.Pid()] + self.node_ids
        msg = ('CIHEX1* %16i%16i%16i%16i\n'
               '*       %16i%16i%16i%16i\n'
               '*       %16i%16i\n' % tuple(data))
        return self.comment + msg

    def __init__(self, eid, pid, nids, comment=''):
        CHEXA8.__init__(self, eid, pid, nids, comment=comment)


class CHEXA20(SolidElement):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+-----+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +=======+=====+=====+=====+=====+=====+=====+=====+=====+
    | CHEXA | EID | PID | G1  | G2  | G3  | G4  | G5  | G6  |
    +-------+-----+-----+-----+-----+-----+-----+-----+-----+
    |       | G7  | G8  | G9  | G10 | G11 | G12 | G13 | G14 |
    +-------+-----+-----+-----+-----+-----+-----+-----+-----+
    |       | G15 | G16 | G17 | G18 | G19 | G20 |     |     |
    +-------+-----+-----+-----+-----+-----+-----+-----+-----+
    """
    type = 'CHEXA'
    aster_type = 'HEXA20'
    calculixType = 'C3D20'

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%8i' % node for node in nodes[8:]]

        data = [self.eid, self.Pid()] + nodes[:8] + nodes2
        msg = ('CHEXA   %8i%8i%8i%8i%8i%8i%8i%8i\n'
               '        %8i%8i%8s%8s%8s%8s%8s%8s\n'
               '        %8s%8s%8s%8s%8s%8s' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def write_card_16(self, is_double=False):
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%8i' % node for node in nodes[8:]]
        data = [self.eid, self.Pid()] + nodes[:8] + nodes2
        msg = ('CHEXA*  %16i%16i%16i%16i\n'
               '*       %16i%16i%16i%16i\n'
               '*       %16i%16i%16s%16s\n'
               '*       %16s%16s%16s%16s\n'
               '*       %16s%16s%16s%16s%16s%16s' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def __init__(self, eid, pid, nids, comment=''):
        SolidElement.__init__(self)

        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        msg = 'len(nids)=%s nids=%s' % (len(nids), nids)
        assert len(self.nodes) <= 20, msg

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'), integer(card, 8, 'nid6'),
            integer(card, 9, 'nid7'), integer(card, 10, 'nid8'),
            integer_or_blank(card, 11, 'nid9'),
            integer_or_blank(card, 12, 'nid10'),
            integer_or_blank(card, 13, 'nid11'),
            integer_or_blank(card, 14, 'nid12'),
            integer_or_blank(card, 15, 'nid13'),
            integer_or_blank(card, 16, 'nid14'),
            integer_or_blank(card, 17, 'nid15'),
            integer_or_blank(card, 18, 'nid16'),
            integer_or_blank(card, 19, 'nid17'),
            integer_or_blank(card, 20, 'nid18'),
            integer_or_blank(card, 21, 'nid19'),
            integer_or_blank(card, 22, 'nid20')
        ]
        assert len(card) <= 23, 'len(CHEXA20 card) = %i\ncard=%s' % (len(card), card)
        return CHEXA20(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = [d if d > 0 else None for d in data[2:]]
        return CHEXA20(eid, pid, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [
            # base
            tuple(sorted([node_ids[0], node_ids[1]])),
            tuple(sorted([node_ids[1], node_ids[2]])),
            tuple(sorted([node_ids[2], node_ids[3]])),
            tuple(sorted([node_ids[3], node_ids[0]])),

            # top
            tuple(sorted([node_ids[4], node_ids[5]])),
            tuple(sorted([node_ids[5], node_ids[6]])),
            tuple(sorted([node_ids[6], node_ids[7]])),
            tuple(sorted([node_ids[7], node_ids[4]])),

            # sides
            tuple(sorted([node_ids[0], node_ids[4]])),
            tuple(sorted([node_ids[1], node_ids[5]])),
            tuple(sorted([node_ids[2], node_ids[6]])),
            tuple(sorted([node_ids[3], node_ids[7]])),
        ]

    def get_face(self, nid_opposite, nid):
        nids = self.node_ids[:8]
        return chexa_face(nid_opposite, nid, nids)

    def getFaceAreaCentroidNormal(self, nid, nid_opposite):
        """
        Parameters
        ----------
        nid : int
            G1 - a grid point on the corner of a face
        nid_opposite : int
            G3 - the grid point diagonally opposite of G1
        """
        nids = self.node_ids[:8]
        return chexa_face_area_centroid_normal(nid, nid_opposite, nids, self.nodes[:8])

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        edges = self.get_edge_ids()
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(self.node_ids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer/blank; nid=%s' %(i, nid)
        if xref:
            centroid = self.Centroid()
            volume = self.Volume()
            assert isinstance(volume, float)
            for i in range(3):
                assert isinstance(centroid[i], float)

    def Centroid(self):
        """
        .. seealso:: CHEXA8.Centroid
        """
        (n1, n2, n3, n4, n5,
         n6, n7, n8) = self.get_node_positions()[:8]
        c1 = area_centroid(n1, n2, n3, n4)[1]
        c2 = area_centroid(n5, n6, n7, n8)[1]
        centroid = (c1 + c2) / 2.
        return centroid

    def Volume(self):
        """
        .. seealso:: CHEXA8.Volume
        """
        (n1, n2, n3, n4, n5,
         n6, n7, n8) = self.get_node_positions()[:8]
        (area1, c1) = area_centroid(n1, n2, n3, n4)
        (area2, c2) = area_centroid(n5, n6, n7, n8)
        volume = (area1 + area2) / 2. * norm(c1 - c2)
        return abs(volume)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)


class CPENTA6(SolidElement):
    r"""
    +--------+-----+-----+----+----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |  8 |  9 |
    +========+=====+=====+====+====+====+====+====+====+
    | CPENTA | EID | PID | G1 | G2 | G3 | G4 | G5 | G6 |
    +--------+-----+-----+----+----+----+----+----+----+

    ::
         3         6
        *----------*
       / \        / \
      / A \      / c \
      *---*-----*-----*
      1    2    4      5
      V = (A1+A2)/2  * norm(c1-c2)
      C = (c1-c2)/2
    """
    type = 'CPENTA'
    aster_type = 'PENTA6'
    calculixType = 'C3D6'

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes
        msg = 'CPENTA  %8i%8i%8i%8i%8i%8i%8i%8i\n' % tuple(data)
        return self.comment + msg

    def write_card_16(self, is_double=False):
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes
        msg = ('CPENTA  %16i%16i%16i%16i\n'
               '        %16i%16i%16i%16i\n' % tuple(data))
        return self.comment + msg

    def __init__(self, eid, pid, nids, comment=''):
        SolidElement.__init__(self)

        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 6

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'),
            integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'),
            integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'),
            integer(card, 8, 'nid6'),
        ]
        assert len(card) == 9, 'len(CPENTA6 card) = %i\ncard=%s' % (len(card), card)
        return CPENTA6(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = data[2:]
        assert len(data) == 8, 'len(data)=%s data=%s' % (len(data), data)
        return CPENTA6(eid, pid, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=False, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [
            # base
            tuple(sorted([node_ids[0], node_ids[1]])),
            tuple(sorted([node_ids[1], node_ids[2]])),
            tuple(sorted([node_ids[2], node_ids[0]])),

            # top
            tuple(sorted([node_ids[3], node_ids[4]])),
            tuple(sorted([node_ids[4], node_ids[5]])),
            tuple(sorted([node_ids[5], node_ids[3]])),

            # sides
            tuple(sorted([node_ids[0], node_ids[3]])),
            tuple(sorted([node_ids[1], node_ids[4]])),
            tuple(sorted([node_ids[2], node_ids[5]])),
        ]

    def get_face(self, nid, nid_opposite=None):
        nids = self.node_ids[:6]
        return cpenta_face(nid, nid_opposite, nids)

    def getFaceAreaCentroidNormal(self, nid, nid_opposite=None):
        nids = self.node_ids[:6]
        return cpenta_face_area_centroid_normal(nid, nid_opposite, nids, self.nodes[:6])

    def getFaceNodesAndArea(self, nid, nid_opposite):
        nids = self.node_ids[:6]
        indx1 = nids.index(nid)
        indx2 = nids.index(nid_opposite)

        #  offset so it's easier to map the nodes with the QRG
        pack = [indx1 + 1, indx2 + 1]
        pack.sort()
        mapper = {
            # reverse points away from the element
            [1, 2]: [1, 2, 3],  # close
            [2, 3]: [1, 2, 3],
            [1, 3]: [1, 2, 3],

            [4, 5]: [4, 5, 6],  # far-reverse
            [5, 6]: [4, 5, 6],
            [4, 6]: [4, 5, 6],

            [1, 5]: [1, 2, 5, 4],  # bottom
            [2, 4]: [1, 2, 5, 4],

            [1, 6]: [1, 3, 6, 4],  # left-reverse
            [3, 4]: [1, 3, 6, 4],

            [2, 6]: [2, 5, 6, 3],  # right
            [3, 5]: [2, 5, 6, 3],
        }

        pack2 = mapper[pack]
        if len(pack2) == 3:
            (n1, n2, n3) = pack2
            face_node_ids = [n1, n2, n3]
            n1i = nids.index(n1 - 1)
            n2i = nids.index(n2 - 1)
            n3i = nids.index(n3 - 1)
            p1 = self.nodes_ref[n1i].get_position()
            p2 = self.nodes_ref[n2i].get_position()
            p3 = self.nodes_ref[n3i].get_position()
            area = 0.5 * norm(cross(p1 - p2, p1 - p3))
        else:
            (n1, n2, n3, n4) = pack2
            n1i = nids.index(n1 - 1)
            n2i = nids.index(n2 - 1)
            n3i = nids.index(n3 - 1)
            n4i = nids.index(n4 - 1)
            face_node_ids = [n1, n2, n3, n4]
            p1 = self.nodes_ref[n1i].get_position()
            p2 = self.nodes_ref[n2i].get_position()
            p3 = self.nodes_ref[n3i].get_position()
            p4 = self.nodes_ref[n4i].get_position()
            area = 0.5 * norm(cross(p1 - p3, p2 - p4))
        return [face_node_ids, area]

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
        if xref:
            centroid = self.Centroid()
            volume = self.Volume()
            assert isinstance(volume, float)
            for i in range(3):
                assert isinstance(centroid[i], float)

    def Centroid(self):
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()
        c1 = (n1 + n2 + n3) / 3.
        c2 = (n4 + n5 + n6) / 3.
        centroid = (c1 + c2) / 2.
        return centroid

    def Volume(self):
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()
        area1 = 0.5 * norm(cross(n3 - n1, n2 - n1))
        area2 = 0.5 * norm(cross(n6 - n4, n5 - n4))
        c1 = (n1 + n2 + n3) / 3.
        c2 = (n4 + n5 + n6) / 3.
        volume = (area1 + area2) / 2. * norm(c1 - c2)
        return abs(volume)

    def raw_fields(self):
        list_fields = ['CPENTA', self.eid, self.Pid()] + self._nodeIDs(allow_empty_nodes=False)
        return list_fields

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=False)

def cpenta_face(nid, nid_opposite, nids):
    assert len(nids) == 6, nids
    indx1 = nids.index(nid)

    if nid_opposite is None:
        if indx1 in [0, 1, 2]:
            pack2 = tuple([2, 1, 0])
        elif indx1 in [3, 4, 5]:
            pack2 = tuple([3, 4, 5])
        else:
            raise RuntimeError(indx1)
        assert len(pack2) == 3, pack2
    else:
        indx2 = nids.index(nid_opposite)

        #  offset so it's easier to map the nodes with the QRG
        pack = tuple(sorted([indx1 + 1, indx2 + 1]))
        _cpenta_mapper = {
            # reverse points away from the element
            #(1, 2) : [1, 2, 3],  # close
            #(2, 3) : [1, 2, 3],
            #(1, 3) : [1, 2, 3],

            #(4, 5) : [4, 5, 6],  # far-reverse
            #(5, 6) : [4, 5, 6],
            #(4, 6) : [4, 5, 6],

            (1, 5) : [4, 5, 2, 1],  # bottom
            (2, 4) : [4, 5, 2, 1],

            (1, 6) : [1, 3, 6, 4],  # left-reverse
            (3, 4) : [1, 3, 6, 4],

            (2, 6) : [2, 5, 6, 3],  # right
            (3, 5) : [2, 5, 6, 3],
        }

        try:
            pack2 = _cpenta_mapper[pack]
        except KeyError:
            print('PLOAD4; remove a node')
            raise
        pack2 = [i - 1 for i in pack2]
    return pack2

def cpenta_face_area_centroid_normal(nid, nid_opposite, nids, nodes_ref):
    """
    Parameters
    ----------
    nid : int
        G1 - a grid point on the corner of a face
    nid_opposite : int / None
        G3 - the grid point diagonally opposite of G1
    """
    face = cpenta_face(nid, nid_opposite, nids)

    if nid_opposite is None:
        n1i, n2i, n3i = face
        p1 = nodes_ref[n1i].get_position()
        p2 = nodes_ref[n2i].get_position()
        p3 = nodes_ref[n3i].get_position()
        a = p3 - p1
        b = p2 - p1
        centroid = (p1 + p2 + p3) / 3.
    else:
        n1i, n2i, n3i, n4i = face
        p1 = nodes_ref[n1i].get_position()
        p2 = nodes_ref[n2i].get_position()
        p3 = nodes_ref[n3i].get_position()
        p4 = nodes_ref[n4i].get_position()
        a = p1 - p3
        b = p2 - p4
        centroid = (p1 + p2 + p3 + p4) / 4.
    normal = cross(a, b)
    n = norm(normal)
    area = 0.5 * n
    return face, area, centroid, normal / n

def chexa_face(nid_opposite, nid, nids):
    assert len(nids) == 8, nids
    g1i = nids.index(nid_opposite)
    g3i = nids.index(nid)

    for face in _chexa_faces:
        if g1i in face and g3i in face:
            found_face = face
    found_face = _chexa_mapper[tuple([g1i, g3i])]
    return found_face

def chexa_face_area_centroid_normal(nid, nid_opposite, nids, nodes_ref):
    """
    Parameters
    ----------
    nid : int
        G1 - a grid point on the corner of a face
    nid_opposite : int
        G3 - the grid point diagonally opposite of G1

    # top   (7-6-5-4)
    # btm   (0-1-2-3)
    # left  (0-3-7-4)
    # right (5-6-2-1)
    # front (4-5-1-0)
    # back  (2-6-7-3)
    """
    face = chexa_face(nid_opposite, nid, nids)
    nid1, nid2, nid3, nid4 = face
    n1 = nodes_ref[nid1].get_position()
    n2 = nodes_ref[nid2].get_position()
    n3 = nodes_ref[nid3].get_position()
    n4 = nodes_ref[nid4].get_position()

    axb = cross(n3 - n1, n4 - n2)
    areai = norm(axb)
    centroid = (n1 + n2 + n3 + n4) / 4.
    area = 0.5 * areai
    normal = axb / areai
    return face, area, centroid, normal


class CPENTA15(SolidElement):
    """
    +---------+-----+-----+----+-----+-----+-----+-----+-----+
    |    1    |  2  |  3  |  4 |  5  |  6  |  7  |  8  |  9  |
    +=========+=====+=====+====+=====+=====+=====+=====+=====+
    |  CPENTA | EID | PID | G1 | G2  | G3  | G4  | G5  | G6  |
    +---------+-----+-----+----+-----+-----+-----+-----+-----+
    |         | G7  | G8  | G9 | G10 | G11 | G12 | G13 | G14 |
    +---------+-----+-----+----+-----+-----+-----+-----+-----+
    |         | G15 |     |    |     |     |     |     |     |
    +---------+-----+-----+----+-----+-----+-----+-----+-----+
    """
    type = 'CPENTA'
    aster_type = 'PENTA15'
    calculixType = 'C3D15'

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%8i' % node for node in nodes[6:]]
        data = [self.eid, self.Pid()] + nodes[:6] + nodes2
        msg = ('CPENTA  %8i%8i%8i%8i%8i%8i%8i%8i\n'
               '        %8s%8s%8s%8s%8s%8s%8s%8s\n'
               '        %8s' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def write_card_16(self, is_double=False):
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%16i' % node for node in nodes[6:]]
        data = [self.eid, self.Pid()] + nodes[:6] + nodes2
        msg = ('CPENTA* %16i%16i%16i%16i\n'
               '*       %16s%16s%16s%16s\n'
               '*       %16s%16s%16s%16s\n'
               '*       %16s%16s%16s%16s\n'
               '*       %16s' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def __init__(self, eid, pid, nids, comment=''):
        SolidElement.__init__(self)

        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) <= 15

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'),
            integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'),
            integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'),
            integer(card, 8, 'nid6'),
            integer_or_blank(card, 9, 'nid7'),
            integer_or_blank(card, 10, 'nid8'),
            integer_or_blank(card, 11, 'nid9'),
            integer_or_blank(card, 12, 'nid10'),
            integer_or_blank(card, 13, 'nid11'),
            integer_or_blank(card, 14, 'nid12'),
            integer_or_blank(card, 15, 'nid13'),
            integer_or_blank(card, 16, 'nid14'),
            integer_or_blank(card, 17, 'nid15'),
        ]
        assert len(card) <= 18, 'len(CPENTA15 card) = %i\ncard=%s' % (len(card), card)
        return CPENTA15(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = [d if d > 0 else None for d in data[2:]]
        assert len(data) == 17, 'len(data)=%s data=%s' % (len(data), data)
        return CPENTA15(eid, pid, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def get_face(self, nid, nid_opposite):
        nids = self.node_ids[:6]
        return cpenta_face(nid_opposite, nid, nids)

    def getFaceAreaCentroidNormal(self, nid, nid_opposite=None):
        nids = self.node_ids[:6]
        return cpenta_face_area_centroid_normal(nid, nid_opposite, nids, self.nodes[:6])

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [
            # base
            tuple(sorted([node_ids[0], node_ids[1]])),
            tuple(sorted([node_ids[1], node_ids[2]])),
            tuple(sorted([node_ids[2], node_ids[0]])),

            # top
            tuple(sorted([node_ids[3], node_ids[4]])),
            tuple(sorted([node_ids[4], node_ids[5]])),
            tuple(sorted([node_ids[5], node_ids[3]])),

            # sides
            tuple(sorted([node_ids[0], node_ids[3]])),
            tuple(sorted([node_ids[1], node_ids[4]])),
            tuple(sorted([node_ids[2], node_ids[5]])),
        ]

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer/blank; nid=%s' %(i, nid)
        if xref:
            centroid = self.Centroid()
            volume = self.Volume()
            assert isinstance(volume, float)
            for i in range(3):
                assert isinstance(centroid[i], float)

    def Centroid(self):
        """
        .. seealso:: CPENTA6.Centroid
        """
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()[:6]
        c1 = (n1 + n2 + n3) / 3.
        c2 = (n4 + n5 + n6) / 3.
        centroid = (c1 - c2) / 2.
        return centroid

    def Volume(self):
        """
        .. seealso:: CPENTA6.Volume
        """
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()[:6]
        area1 = Area(n3 - n1, n2 - n1)
        area2 = Area(n6 - n4, n5 - n4)
        c1 = (n1 + n2 + n3) / 3.
        c2 = (n4 + n5 + n6) / 3.

        volume = (area1 + area2) / 2. * norm(c1 - c2)
        return abs(volume)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)


class CPYRAM5(SolidElement):
    """
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    |    1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+=====+=====+=====+=====+=====+=====+=====+=====+
    | CPYRAM | EID | PID | G1  | G2  | G3  | G4  | G5  |     |
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    """
    type = 'CPYRAM'
    #aster_type = 'CPYRAM5'
    #calculixType = 'C3D5'

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes
        msg = ('CPYRAM  %8i%8i%8i%8i%8i%8i%8i' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def write_card_16(self, is_double=False):
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes
        msg = ('CPYRAM  %16i%16i%16i%16i\n'
               '        %16i%16i%16i' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def __init__(self, eid, pid, nids, comment=''):
        SolidElement.__init__(self)

        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids)
        msg = 'len(nids)=%s nids=%s' % (len(nids), nids)
        assert len(self.nodes) <= 20, msg

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
                integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
                integer(card, 7, 'nid5')]
        assert len(card) == 8, 'len(CPYRAM5 1card) = %i\ncard=%s' % (len(card), card)
        return CPYRAM5(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = [d if d > 0 else None for d in data[2:]]
        return CPYRAM5(eid, pid, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [
            # base
            tuple(sorted([node_ids[0], node_ids[1]])),
            tuple(sorted([node_ids[1], node_ids[2]])),
            tuple(sorted([node_ids[2], node_ids[3]])),
            tuple(sorted([node_ids[3], node_ids[0]])),

            # sides
            tuple(sorted([node_ids[0], node_ids[4]])),
            tuple(sorted([node_ids[1], node_ids[4]])),
            tuple(sorted([node_ids[2], node_ids[4]])),
            tuple(sorted([node_ids[3], node_ids[4]])),
        ]

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer/blank; nid=%s' %(i, nid)
        if xref:
            centroid = self.Centroid()
            volume = self.Volume()
            assert isinstance(volume, float)
            for i in range(3):
                assert isinstance(centroid[i], float)

    def Centroid(self):
        """
        .. seealso:: CPYRAM5.Centroid
        """
        (n1, n2, n3, n4, n5) = self.get_node_positions()
        c1 = area_centroid(n1, n2, n3, n4)[1]
        centroid = (c1 + n5) / 2.
        return centroid

    def Volume(self):
        """
        .. seealso:: CPYRAM5.Volume
        """
        (n1, n2, n3, n4, n5) = self.get_node_positions()
        area1, c1 = area_centroid(n1, n2, n3, n4)
        volume = area1 / 3. * norm(c1 - n5)
        return abs(volume)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=False)


class CPYRAM13(SolidElement):
    """
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    |    1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+=====+=====+=====+=====+=====+=====+=====+=====+
    | CPYRAM | EID | PID | G1  | G2  | G3  | G4  | G5  | G6  |
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    |        | G7  | G8  | G9  | G10 | G11 | G12 |     |     |
    +--------+-----+-----+-----+-----+-----+-----+-----+-----+
    """
    type = 'CPYRAM'
    #aster_type = 'CPYRAM13'
    #calculixType = 'C3D13'

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%8i' % node for node in nodes[5:]]
        data = [self.eid, self.Pid()] + nodes[:5] + nodes2
        msg = ('CPYRAM  %8i%8i%8i%8i%8i%8i%8i%8s\n'
               '        %8s%8s%8s%8s%8s%8s%s' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def write_card_16(self, is_double=False):
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%16i' % node for node in nodes[5:]]
        data = [self.eid, self.Pid()] + nodes[:5] + nodes2
        msg = ('CPYRAM* %16i%16i%16i%16i\n'
               '*       %16i%16i%16i%16s\n'
               '*       %16s%16s%16s%16s\n'
               '*       %16s%16s%s\n' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def __init__(self, eid, pid, nids, comment=''):
        SolidElement.__init__(self)

        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        msg = 'len(nids)=%s nids=%s' % (len(nids), nids)
        assert len(self.nodes) <= 13, msg

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [
            integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'),
            integer_or_blank(card, 8, 'nid6'),
            integer_or_blank(card, 9, 'nid7'),
            integer_or_blank(card, 10, 'nid8'),
            integer_or_blank(card, 11, 'nid9'),
            integer_or_blank(card, 12, 'nid10'),
            integer_or_blank(card, 13, 'nid11'),
            integer_or_blank(card, 14, 'nid12'),
            integer_or_blank(card, 15, 'nid13')
        ]
        assert len(card) <= 16, 'len(CPYRAM13 1card) = %i\ncard=%s' % (len(card), card)
        return CPYRAM13(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = [d if d > 0 else None for d in data[2:]]
        return CPYRAM13(eid, pid, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [
            # base
            tuple(sorted([node_ids[0], node_ids[1]])),
            tuple(sorted([node_ids[1], node_ids[2]])),
            tuple(sorted([node_ids[2], node_ids[3]])),
            tuple(sorted([node_ids[3], node_ids[0]])),

            # sides
            tuple(sorted([node_ids[0], node_ids[4]])),
            tuple(sorted([node_ids[1], node_ids[4]])),
            tuple(sorted([node_ids[2], node_ids[4]])),
            tuple(sorted([node_ids[3], node_ids[4]])),
        ]

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer/blank; nid=%s' %(i, nid)
        if xref:
            centroid = self.Centroid()
            volume = self.Volume()
            assert isinstance(volume, float)
            for i in range(3):
                assert isinstance(centroid[i], float)

    def Centroid(self):
        """
        .. seealso:: CPYRAM5.Centroid
        """
        (n1, n2, n3, n4, n5,
         n6, n7, n8, n9, n10,
         n11, n12, n13) = self.get_node_positions()
        c1 = area_centroid(n1, n2, n3, n4)[1]
        centroid = (c1 + n5) / 2.
        return centroid

    def Volume(self):
        """
        .. seealso:: CPYRAM5.Volume
        """
        (n1, n2, n3, n4, n5,
         n6, n7, n8, n9, n10,
         n11, n12, n13) = self.get_node_positions()
        area1, c1 = area_centroid(n1, n2, n3, n4)
        volume = area1 / 2. * norm(c1 - n5)
        return abs(volume)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)


class CTETRA4(SolidElement):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CTETRA | EID | PID | G1 | G2 | G3 | G4 |
    +--------+-----+-----+----+----+----+----+
    """
    type = 'CTETRA'
    aster_type = 'TETRA4'
    calculixType = 'C3D4'

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes
        msg = 'CTETRA  %8i%8i%8i%8i%8i%8i\n' % tuple(data)
        return self.comment + msg

    def write_card_16(self, is_double=False):
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes
        msg = ('CTETRA  %16i%16i%16i%16i\n'
               '        %16i%16i\n' % tuple(data))
        return self.comment + msg

    def __init__(self, eid, pid, nids, comment=''):
        SolidElement.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids)
        assert len(self.nodes) == 4

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'nid1'),
                integer(card, 4, 'nid2'),
                integer(card, 5, 'nid3'),
                integer(card, 6, 'nid4'), ]
        assert len(card) == 7, 'len(CTETRA4 card) = %i\ncard=%s' % (len(card), card)
        return CTETRA4(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = data[2:]
        assert len(data) == 6, 'len(data)=%s data=%s' % (len(data), data)
        return CTETRA4(eid, pid, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=False, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)
        if xref:
            centroid = self.Centroid()
            volume = self.Volume()
            assert isinstance(volume, float)
            for i in range(3):
                assert isinstance(centroid[i], float)

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [
            # base
            tuple(sorted([node_ids[0], node_ids[1]])),
            tuple(sorted([node_ids[1], node_ids[2]])),
            tuple(sorted([node_ids[2], node_ids[0]])),

            # sides
            tuple(sorted([node_ids[0], node_ids[3]])),
            tuple(sorted([node_ids[1], node_ids[3]])),
            tuple(sorted([node_ids[2], node_ids[3]])),
        ]

    def Volume(self):
        (n1, n2, n3, n4) = self.get_node_positions()
        return volume4(n1, n2, n3, n4)

    def Centroid(self):
        (n1, n2, n3, n4) = self.get_node_positions()
        return (n1 + n2 + n3 + n4) / 4.

    def getFaceNodes(self, nid_opposite, nid=None):
        assert nid is None, nid
        nids = self.node_ids[:4]
        indx = nids.index(nid_opposite)
        nids.pop(indx)
        return nids

    def get_face(self, nid_opposite, nid):
        nids = self.node_ids[:6]
        return ctetra_face(nid_opposite, nid, nids)

    def getFaceAreaCentroidNormal(self, nid, nid_opposite):
        return ctetra_face_area_centroid_normal(nid, nid_opposite,
                                                self.node_ids, self.nodes)

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=False)

def ctetra_face(nid, nid_opposite, nids):
    assert len(nids) == 4, nids
    g1i = nids.index(nid)
    g4i = nids.index(nid_opposite)

    _ctetra_faces = (
        (3, 1, 0),
        (0, 1, 2),
        (3, 2, 1),
        (0, 2, 3),
    )
    for face in _ctetra_faces:
        if g1i in face and g4i not in face:
            found_face = face
    return found_face

def ctetra_face_area_centroid_normal(nid, nid_opposite, nids, nodes_ref):
    """
    Parameters
    ----------
    nid : int
        G1 - a grid point on the corner of a face
    nid_opposite : int
        G4 - a grid point not being loaded
    """
    face = ctetra_face(nid, nid_opposite, nids)
    nid1, nid2, nid3 = face
    n1 = nodes_ref[nid1].get_position()
    n2 = nodes_ref[nid2].get_position()
    n3 = nodes_ref[nid3].get_position()

    axb = cross(n2 - n1, n3 - n1)
    normi = norm(axb)
    centroid = (n1 + n2 + n3) / 3.
    area = 0.5 * normi
    assert area > 0, area
    normal = axb / normi
    return face, area, centroid, normal


class CTETRA10(SolidElement):
    """
    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    |    1   |  2  |  3  |  4  |  5  |  6  |  7 |  8  |  9  |
    +========+=====+=====+=====+=====+=====+====+=====+=====+
    | CTETRA | EID | PID | G1  | G2  | G3  | G4 | G5  | G6  |
    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    |        | G7  |  G8 | G9  | G10 |     |    |     |     |
    +--------+-----+-----+-----+-----+-----+----+-----+-----+

    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    | CTETRA | 1   | 1   | 239 | 229 | 516 | 99 | 335 | 103 |
    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    |        | 265 | 334 | 101 | 102 |     |    |     |     |
    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    """
    type = 'CTETRA'
    aster_type = 'TETRA10'
    calculixType = 'C3D10'

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%8i' % node for node in nodes[4:]]

        data = [self.eid, self.Pid()] + nodes[:4] + nodes2
        msg = ('CTETRA  %8i%8i%8i%8i%8i%8i%8s%8s\n'
               '        %8s%8s%8s%8s' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def write_card_16(self, is_double=False):
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%16i' % node for node in nodes[4:]]
        data = [self.eid, self.Pid()] + nodes[:4] + nodes2
        msg = ('CTETRA  %16i%16i%16i%16i\n'
               '        %16i%16i%16s%16s\n'
               '        %16s%16s%16s%16s' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def getFaceAreaCentroidNormal(self, nid_opposite, nid=None):
        return ctetra_face_area_centroid_normal(nid_opposite, nid,
                                                self.node_ids[:4], self.nodes[:4])

    def __init__(self, eid, pid, nids, comment=''):
        SolidElement.__init__(self)
        if comment:
            self._comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) <= 10

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'nid1'),
                integer(card, 4, 'nid2'),
                integer(card, 5, 'nid3'),
                integer(card, 6, 'nid4'),
                integer_or_blank(card, 7, 'nid5'),
                integer_or_blank(card, 8, 'nid6'),
                integer_or_blank(card, 9, 'nid7'),
                integer_or_blank(card, 10, 'nid8'),
                integer_or_blank(card, 11, 'nid9'),
                integer_or_blank(card, 12, 'nid10'), ]
        assert len(card) <= 13, 'len(CTETRA10 card) = %i\ncard=%s' % (len(card), card)
        return CTETRA10(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        pid = data[1]
        nids = [d if d > 0 else None for d in data[2:]]
        assert len(data) == 12, 'len(data)=%s data=%s' % (len(data), data)
        return CTETRA10(eid, pid, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [
            # base
            tuple(sorted([node_ids[0], node_ids[1]])),
            tuple(sorted([node_ids[1], node_ids[2]])),
            tuple(sorted([node_ids[2], node_ids[0]])),

            # sides
            tuple(sorted([node_ids[0], node_ids[3]])),
            tuple(sorted([node_ids[1], node_ids[3]])),
            tuple(sorted([node_ids[2], node_ids[3]])),
        ]

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        assert isinstance(eid, int)
        assert isinstance(pid, int)
        for i, nid in enumerate(nids):
            assert nid is None or isinstance(nid, int), 'nid%i is not an integer/blank; nid=%s' % (i, nid)
        if xref:
            centroid = self.Centroid()
            volume = self.Volume()
            assert isinstance(volume, float)
            for i in range(3):
                assert isinstance(centroid[i], float)

    #def N_10(self, g1, g2, g3, g4):
        #N1 = g1 * (2 * g1 - 1)
        #N2 = g2 * (2 * g2 - 1)
        #N3 = g3 * (2 * g3 - 1)
        #N4 = g4 * (2 * g4 - 1)
        #N5 = 4 * g1 * g2
        #N6 = 4 * g2 * g3
        #N7 = 4 * g3 * g1
        #N8 = 4 * g1 * g4
        #N9 = 4 * g2 * g4
        #N10 = 4 * g3 * g4
        #return (N1, N2, N3, N4, N5, N6, N7, N8, N9, N10)

    def Volume(self):
        """
        Gets the volume, :math:`V`, of the primary tetrahedron.

        .. seealso:: CTETRA4.Volume
        """
        (n1, n2, n3, n4) = self.get_node_positions()[:4]
        return volume4(n1, n2, n3, n4)

    def Centroid(self):
        """
        Gets the cenroid of the primary tetrahedron.

        .. seealso:: CTETRA4.Centroid
        """
        (n1, n2, n3, n4) = self.get_node_positions()[:4]
        return (n1 + n2 + n3 + n4) / 4.

    def getFaceNodes(self, nid_opposite, nid=None):
        nids = self.node_ids[:4]
        indx = nids.index(nid_opposite)
        nids.pop(indx)
        return nids

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

