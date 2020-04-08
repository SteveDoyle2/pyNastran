# pylint: disable=R0902,R0904,R0914
"""
All solid elements are defined in this file.  This includes:

 * CHEXA8
 * CHEXA20
 * CPENTA6
 * CPENTA15
 * CTETRA4
 * CTETRA10
 * CIHEX1
 * CIHEX2

All solid elements are SolidElement and Element objects.

"""
from __future__ import annotations
from typing import Tuple, Any, TYPE_CHECKING
import numpy as np
from numpy import dot, cross
from numpy.linalg import norm  # type: ignore

from pyNastran.bdf.cards.elements.elements import Element
from pyNastran.utils.mathematics import Area
from pyNastran.bdf.bdf_interface.assign_type import integer, integer_or_blank
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


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

def volume4(n1: Any, n2: Any, n3: Any, n4: Any) -> float:
    r"""
    Gets the volume, :math:`V`, of the tetrahedron.

    .. math:: V = \frac{(a-d) \cdot \left( (b-d) \times (c-d) \right) }{6}
    """
    volume = -dot(n1 - n4, cross(n2 - n4, n3 - n4)) / 6.
    return volume

def area_centroid(n1: Any, n2: Any, n3: Any, n4: Any) -> Tuple[float, float]:
    """
    Gets the area, :math:`A`, and centroid of a quad.::

      1-----2
      |   / |
      | /   |
      4-----3
    """
    area = 0.5 * norm(cross(n3 - n1, n4 - n2))
    centroid = (n1 + n2 + n3 + n4) / 4.
    return area, centroid


nnodes_map = {
    'CTETRA' : (4, 10),
    'CPENTA' : (6, 15),
    'CPYRAM' : (5, 13),
    'CHEXA' : (8, 20),
}
class SolidElement(Element):
    _field_map = {1: 'nid', 2:'pid'}
    _properties = ['faces']

    def __init__(self):
        Element.__init__(self)
        self.nodes_ref = None
        self.pid_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        nnodes = nnodes_map[cls.type]
        comments = []
        pids = []

        element0 = model.elements[eids[0]]
        nnodes0 = len(element0.nodes)
        nnodes_high_map = {
            4 : 10, 10 : 10, # CTETRA
            5 : 13, 13 : 13, # CYRAM
            6 : 15, 15 : 15, # CPENTA
            8 : 20, 20 : 20, # CHEXA
        }
        nnodes_low_map = {
            4 : 4, 10 : 4, # CTETRA
            5 : 5, 13 : 5, # CYRAM
            6 : 6, 15 : 6, # CPENTA
            8 : 8, 20 : 8, # CHEXA
        }
        neids = len(eids)
        nnodes = nnodes_high_map[nnodes0]
        nnodes_low = nnodes_low_map[nnodes0]

        nodes = np.zeros((neids, nnodes), dtype='int32')
        for i, eid in enumerate(eids):
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes[i, :len(element.nodes)] = [nid if nid is not None else 0 for nid in element.nodes]
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)

        if nodes[:, nnodes_low:].max() == 0:
            nodes = nodes[:, :nnodes_low]
        h5_file.create_dataset('nodes', data=nodes)

    def _update_field_helper(self, n, value):
        if n - 3 < len(self.nodes):
            self.nodes[n - 3] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def cross_reference(self, model: BDF) -> None:
        raise NotImplementedError('Element type=%r must implement cross_reference')

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def E(self) -> float:
        return self.pid_ref.mid_ref.E()

    def G(self) -> float:
        return self.pid_ref.mid_ref.G()

    def Nu(self) -> float:
       return self.pid_ref.mid_ref.Nu()

    def Volume(self) -> float:
        """
        Base volume method that should be overwritten
        """
        return 0.

    def Mass(self) -> float:
        """
        Calculates the mass of the solid element
        Mass = Rho * Volume
        """
        #print('  rho=%e volume=%e' % (self.Rho(), self.Volume()))
        return self.Rho() * self.Volume()

    def Mid(self) -> int:
        """
        Returns the material ID as an integer
        """
        return self.pid_ref.Mid()

    def Rho(self) -> float:
        """
        Returns the density
        """
        try:
            return self.pid_ref.Rho()
        except AttributeError:
            print("self.pid = %s" % (self.pid))
            #print("self.pid_ref.mid_ref = %s" % (str(self.pid_ref.mid_ref)))
            raise

    def get_face_area_centroid_normal(self, nid_opposite, nid=None):
        return self.get_face_area_centroid_normal(nid_opposite, nid)

    def raw_fields(self):
        list_fields = [self.type, self.eid, self.Pid()] + self.node_ids
        return list_fields

    def center_of_mass(self):
        return self.Centroid()


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
    def write_card(self, size: int=8, is_double: bool=False) -> str:
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
        """
        Creates a CHEXA8

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : List[int]
            node ids; n=8

        """
        SolidElement.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 8

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CHEXA8 card from ``BDF.add_card(...)``

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
        """
        Adds a CHEXA8 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        pid = data[1]
        nids = data[2:]
        assert len(data) == 10, 'len(data)=%s data=%s' % (len(data), data)
        return CHEXA8(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CHEXA eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CHEXA eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    @property
    def faces(self):
        """
        Gets the faces of the element

        Returns
        -------
        faces : Dict[int] = [face1, face2, ...]
            key = face number
            value = a list of nodes (integer pointers) as the values.

        .. note::  The order of the nodes are consistent with normals that point outwards
                   The face numbering is meaningless

        .. note::  The order of the nodes are consistent with ANSYS numbering; is this current?
        .. warning:: higher order element ids not verified with ANSYS; is this current?

        Examples
        --------
        >>> print(element.faces)
        """
        nodes = self.node_ids
        faces = {
            1 : [nodes[0], nodes[1], nodes[2], nodes[3]],
            2 : [nodes[0], nodes[1], nodes[5], nodes[4]],
            3 : [nodes[1], nodes[2], nodes[6], nodes[5]],
            4 : [nodes[2], nodes[3], nodes[7], nodes[6]],
            5 : [nodes[3], nodes[0], nodes[4], nodes[7]],
            6 : [nodes[4], nodes[5], nodes[6], nodes[7]],
        }
        return faces

    def material_coordinate_system(self, xyz=None):
        """http://www.ipes.dk/Files/Ipes/Filer/nastran_2016_doc_release.pdf"""
        #if normal is None:
            #normal = self.Normal() # k = kmat

        if xyz is None:
            x1 = self.nodes_ref[0].get_position()
            x2 = self.nodes_ref[1].get_position()
            x3 = self.nodes_ref[2].get_position()
            x4 = self.nodes_ref[3].get_position()
            x5 = self.nodes_ref[4].get_position()
            x6 = self.nodes_ref[5].get_position()
            x7 = self.nodes_ref[6].get_position()
            x8 = self.nodes_ref[7].get_position()
        else:
            x1 = xyz[:, 0]
            x2 = xyz[:, 1]
            x3 = xyz[:, 2]
            x4 = xyz[:, 3]
            x5 = xyz[:, 4]
            x6 = xyz[:, 5]
            x7 = xyz[:, 6]
            x8 = xyz[:, 7]

        #CORDM=-2
        centroid = (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8) / 8.
        xe = ((x2+x3+x6+x7) - (x1+x4+x8+x5)) / 4.
        xe /= np.linalg.norm(xe)
        v = ((x3+x7+x8+x4) - (x1+x2+x6+x5)) / 4
        ze = np.cross(xe, v)
        ze /= np.linalg.norm(ze)

        ye = np.cross(ze, xe)
        ye /= np.linalg.norm(ye)
        return centroid, xe, ye, ze

    def _verify(self, xref):
        eid = self.eid
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
        """Calculate the volume of the hex"""
        # https://www.osti.gov/servlets/purl/632793/
        #volume = (
            #det3(x7 - x0, x1 - x0, x3 - x5) +
            #det3(x7 - x0, x4 - x0, x5 - x6) +
            #det3(x7 - x0, x2 - x0, x6 - x3)
        #) / 6.
        #  swap points
        # x2 -> x3
        # x3 -> x2
        #
        # x6 -> x7
        # x7 -> x6
        #def det3(a, b, c):
            #return np.det(np.vstack(a, b, c))
        #volume = (
            #det3(x6 - x0, x1 - x0, x2 - x5) +
            #det3(x6 - x0, x4 - x0, x5 - x7) +
            #det3(x6 - x0, x3 - x0, x7 - x3)
        #) / 6.
        # add 1
        #volume = (
            #det3(x7 - x1, x2 - x1, x3 - x6) +
            #det3(x7 - x1, x5 - x1, x6 - x8) +
            #det3(x7 - x1, x4 - x1, x8 - x4)
        #) / 6.
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.get_node_positions()
        (area1, c1) = area_centroid(n1, n2, n3, n4)
        (area2, c2) = area_centroid(n5, n6, n7, n8)
        volume = (area1 + area2) / 2. * norm(c1 - c2)
        return abs(volume)

    @property
    def node_ids(self):
        nids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)
        return nids

    def get_face(self, nid_opposite, nid):
        nids = self.node_ids[:8]
        return chexa_face(nid_opposite, nid, nids)

    def get_face_area_centroid_normal(self, nid, nid_opposite):
        """
        Parameters
        ----------
        nid : int
            G1 - a grid point on the corner of a face
        nid_opposite : int
            G3 - the grid point diagonally opposite of G1
        """
        nids = self.node_ids[:8]
        return chexa_face_area_centroid_normal(nid, nid_opposite, nids, self.nodes_ref[:8])

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

    def flip_normal(self):  ## TODO verify
        """flips the element inside out"""
        # reverse the lower and upper quad faces
        n1, n2, n3, n4, n5, n6, n7, n8 = self.nodes
        self.nodes = [n1, n4, n3, n2,
                      n5, n8, n7, n6,]
        if self.nodes_ref is not None:
            n1_ref, n2_ref, n3_ref, n4_ref, n5_ref, n6_ref, n7_ref, n8_ref = self.nodes_ref
            self.nodes_ref = [
                n1_ref, n4_ref, n3_ref, n2_ref,
                n5_ref, n8_ref, n7_ref, n6_ref,]

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
    def write_card(self, size: int=8, is_double: bool=False) -> str:
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
        """
        Creates a CHEXA20

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : List[int]
            node ids; n=20
        """
        SolidElement.__init__(self)

        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid

        nnodes = len(nids)
        if nnodes < 20:
            nids.extend((20 - nnodes) * [None])
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        msg = 'len(nids)=%s nids=%s' % (len(nids), nids)
        assert len(self.nodes) == 20, msg

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CHEXA20 card from ``BDF.add_card(...)``

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
            integer_or_blank(card, 22, 'nid20'),
        ]
        assert len(card) <= 23, 'len(CHEXA20 card) = %i\ncard=%s' % (len(card), card)
        return CHEXA20(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CHEXA20 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        pid = data[1]
        nids = [d if d > 0 else None for d in data[2:]]
        return CHEXA20(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CHEXA eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CHEXA eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    @property
    def faces(self):
        """
        Gets the faces of the element

        Returns
        -------
        faces : Dict[int] = [face1, face2, ...]
            key = face number
            value = a list of nodes (integer pointers) as the values.

        .. note::  The order of the nodes are consistent with normals that point outwards
                   The face numbering is meaningless

        .. note::  The order of the nodes are consistent with ANSYS numbering; is this current?
        .. warning:: higher order element ids not verified with ANSYS; is this current?

        Examples
        --------
        >>> print(element.faces)
        """
        (n1, n2, n3, n4, n5, n6, n7, n8,
         n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20) = self.node_ids
        faces = {
            1 : [n1, n2, n3, n4, n9, n10, n11, n12],
            2 : [n1, n2, n6, n5, n9, n18, n13, n17],
            3 : [n2, n3, n7, n6, n10, n19, n14, n18],
            4 : [n3, n4, n8, n7, n11, n10, n15, n19],
            5 : [n4, n1, n5, n8, n12, n17, n16, n20],
            6 : [n5, n6, n7, n8, n13, n14, n15, n16],
        }
        return faces

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

    def get_face_area_centroid_normal(self, nid, nid_opposite):
        """
        Parameters
        ----------
        nid : int
            G1 - a grid point on the corner of a face
        nid_opposite : int
            G3 - the grid point diagonally opposite of G1
        """
        nids = self.node_ids[:8]
        return chexa_face_area_centroid_normal(nid, nid_opposite, nids, self.nodes_ref[:8])

    def _verify(self, xref):
        eid = self.eid
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

    @property
    def node_ids(self):
        nids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)
        return nids


class CIHEX1(CHEXA8):
    type = 'CIHEX1'

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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


class CIHEX2(CHEXA20):
    type = 'CIHEX2'

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%8i' % node for node in nodes[8:]]

        data = [self.eid, self.Pid()] + nodes[:8] + nodes2
        msg = ('CIHEX2  %8i%8i%8i%8i%8i%8i%8i%8i\n'
               '        %8i%8i%8s%8s%8s%8s%8s%8s\n'
               '        %8s%8s%8s%8s%8s%8s' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def write_card_16(self, is_double=False):
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%8i' % node for node in nodes[8:]]
        data = [self.eid, self.Pid()] + nodes[:8] + nodes2
        msg = ('CIHEX2* %16i%16i%16i%16i\n'
               '*       %16i%16i%16i%16i\n'
               '*       %16i%16i%16s%16s\n'
               '*       %16s%16s%16s%16s\n'
               '*       %16s%16s%16s%16s%16s%16s' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

    def __init__(self, eid, pid, nids, comment=''):
        CHEXA20.__init__(self, eid, pid, nids, comment=comment)


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
    def write_card(self, size: int=8, is_double: bool=False) -> str:
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
        """
        Creates a CPENTA6

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : List[int]
            node ids; n=6
        """
        SolidElement.__init__(self)

        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 6

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CPENTA6 card from ``BDF.add_card(...)``

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
        """
        Adds a CPENTA6 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        pid = data[1]
        nids = data[2:]
        assert len(data) == 8, 'len(data)=%s data=%s' % (len(data), data)
        return CPENTA6(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CPENTA eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CPENTA eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def material_coordinate_system(self, xyz=None):
        """http://www.ipes.dk/Files/Ipes/Filer/nastran_2016_doc_release.pdf"""
        #if normal is None:
            #normal = self.Normal() # k = kmat

        if xyz is None:
            x1 = self.nodes_ref[0].get_position()
            x2 = self.nodes_ref[1].get_position()
            x3 = self.nodes_ref[2].get_position()
            x4 = self.nodes_ref[3].get_position()
            x5 = self.nodes_ref[4].get_position()
            x6 = self.nodes_ref[5].get_position()
        else:
            x1 = xyz[:, 0]
            x2 = xyz[:, 1]
            x3 = xyz[:, 2]
            x4 = xyz[:, 3]
            x5 = xyz[:, 4]
            x6 = xyz[:, 5]

        #CORDM=-2
        centroid = self.Centroid()
        origin = (x1 + x4) / 2.
        xe = (x2 + x3 + x5 + x6) - origin
        xe /= np.linalg.norm(xe)
        v = ((x1 + x3 + x4 + x6) - (x1 + x2 + x4 + x5)) / 4.
        ze = np.cross(xe, v)
        ze /= np.linalg.norm(ze)

        ye = np.cross(ze, xe)
        ye /= np.linalg.norm(ye)
        return centroid, xe, ye, ze

    @property
    def faces(self):
        """
        Gets the faces of the element

        Returns
        -------
        faces : Dict[int] = [face1, face2, ...]
            key = face number
            value = a list of nodes (integer pointers) as the values.

        .. note::  The order of the nodes are consistent with normals that point outwards
                   The face numbering is meaningless

        .. note::  The order of the nodes are consistent with ANSYS numbering; is this current?
        .. warning:: higher order element ids not verified with ANSYS; is this current?

        Examples
        --------
        >>> print(element.faces)
        """
        nodes = self.node_ids
        faces = {
            1 : [nodes[0], nodes[1], nodes[2]],
            2 : [nodes[3], nodes[4], nodes[5]],
            3 : [nodes[0], nodes[1], nodes[4], nodes[3]],
            4 : [nodes[1], nodes[2], nodes[5], nodes[4]],
            5 : [nodes[2], nodes[0], nodes[3], nodes[5]],
        }
        return faces

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

    def get_face_area_centroid_normal(self, nid, nid_opposite=None):
        nids = self.node_ids[:6]
        return cpenta_face_area_centroid_normal(nid, nid_opposite, nids, self.nodes_ref[:6])

    def get_face_nodes_and_area(self, nid, nid_opposite):
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

    def _verify(self, xref):
        eid = self.eid
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
        """Calculate the volume of the penta"""
        (n1, n2, n3, n4, n5, n6) = self.get_node_positions()
        area1 = 0.5 * norm(cross(n3 - n1, n2 - n1))
        area2 = 0.5 * norm(cross(n6 - n4, n5 - n4))
        c1 = (n1 + n2 + n3) / 3.
        c2 = (n4 + n5 + n6) / 3.
        volume = (area1 + area2) / 2. * norm(c1 - c2)
        return abs(volume)
        #return volume4(n1, n2, n3, n4) + volume4(n2, n3, n4, n5) + volume4(n2, n4, n5, n6)

    def raw_fields(self):
        list_fields = ['CPENTA', self.eid, self.Pid()] + self.node_ids
        return list_fields

    @property
    def node_ids(self):
        nids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)
        return nids

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
    nodes_ref : List[GRID]
        the GRID objects

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
    def __init__(self, eid, pid, nids, comment=''):
        """
        Creates a CPENTA15

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : List[int]
            node ids; n=15
        """
        SolidElement.__init__(self)

        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        nnodes = len(nids)
        if nnodes < 15:
            nids.extend((15 - nnodes) * [None])
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 15

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CPENTA15 card from ``BDF.add_card(...)``

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
        """
        Adds a CPENTA15 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        pid = data[1]
        nids = [d if d > 0 else None for d in data[2:]]
        assert len(data) == 17, 'len(data)=%s data=%s' % (len(data), data)
        return CPENTA15(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CPENTA eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CPENTA eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    @property
    def faces(self):
        """
        Gets the faces of the element

        Returns
        -------
        faces : Dict[int] = [face1, face2, ...]
            key = face number
            value = a list of nodes (integer pointers) as the values.

        .. note::  The order of the nodes are consistent with normals that point outwards
                   The face numbering is meaningless

        .. note::  The order of the nodes are consistent with ANSYS numbering; is this current?
        .. warning:: higher order element ids not verified with ANSYS; is this current?

        Examples
        --------
        >>> print(element.faces)
        """
        n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15 = self.node_ids
        faces = {
            1 : [n1, n2, n3, n7, n8, n9],
            2 : [n4, n5, n6, n10, n11, n12],
            3 : [n1, n2, n5, n4, n7, n14, n10, n13],
            4 : [n2, n3, n6, n5, n8, n15, n11, n14],
            5 : [n3, n1, n4, n6, n9, n13, n12, n15],
        }
        return faces

    def get_face(self, nid, nid_opposite):
        nids = self.node_ids[:6]
        return cpenta_face(nid_opposite, nid, nids)

    def get_face_area_centroid_normal(self, nid, nid_opposite=None):
        nids = self.node_ids[:6]
        return cpenta_face_area_centroid_normal(nid, nid_opposite, nids, self.nodes_ref[:6])

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

    def _verify(self, xref):
        eid = self.eid
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
        centroid = (c1 + c2) / 2.
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

    @property
    def node_ids(self):
        nids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)
        return nids

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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


class CPYRAM5(SolidElement):
    """
    +--------+-----+-----+-----+-----+-----+-----+-----+
    |    1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |
    +========+=====+=====+=====+=====+=====+=====+=====+
    | CPYRAM | EID | PID | G1  | G2  | G3  | G4  | G5  |
    +--------+-----+-----+-----+-----+-----+-----+-----+
    """
    type = 'CPYRAM'
    def __init__(self, eid, pid, nids, comment=''):
        SolidElement.__init__(self)

        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.nodes = self.prepare_node_ids(nids)
        msg = 'len(nids)=%s nids=%s' % (len(nids), nids)
        assert len(self.nodes) <= 20, msg

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CPYRAM5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
                integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
                integer(card, 7, 'nid5')]
        assert len(card) == 8, 'len(CPYRAM5 1card) = %i\ncard=%s' % (len(card), card)
        return CPYRAM5(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CPYRAM5 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        pid = data[1]
        nids = [d if d > 0 else None for d in data[2:]]
        return CPYRAM5(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CPYRAM eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CPYRAM eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    @property
    def faces(self):
        """
        Gets the faces of the element

        Returns
        -------
        faces : Dict[int] = [face1, face2, ...]
            key = face number
            value = a list of nodes (integer pointers) as the values.

        .. note::  The order of the nodes are consistent with normals that point outwards
                   The face numbering is meaningless

        .. note::  The order of the nodes are consistent with ANSYS numbering; is this current?
        .. warning:: higher order element ids not verified with ANSYS; is this current?

        Examples
        --------
        >>> print(element.faces)
        """
        nodes = self.node_ids
        faces = {
            1 : [nodes[0], nodes[1], nodes[2], nodes[3]],
            2 : [nodes[0], nodes[1], nodes[4]],
            3 : [nodes[1], nodes[2], nodes[4]],
            4 : [nodes[2], nodes[3], nodes[4]],
            5 : [nodes[3], nodes[0], nodes[4]],
        }
        return faces

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

    def _verify(self, xref):
        eid = self.eid
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

        V = (l * w) * h / 3
        V = A * h / 3
        """
        (n1, n2, n3, n4, n5) = self.get_node_positions()
        area1, c1 = area_centroid(n1, n2, n3, n4)
        volume = area1 / 3. * norm(c1 - n5)
        return abs(volume)

    @property
    def node_ids(self):
        nids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)
        return nids

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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
    def __init__(self, eid, pid, nids, comment=''):
        SolidElement.__init__(self)

        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        nnodes = len(nids)
        if nnodes < 13:
            nids.extend((13 - nnodes) * [None])
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        msg = 'len(nids)=%s nids=%s' % (len(nids), nids)
        assert len(self.nodes) == 13, msg

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CPYRAM13 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        """
        Adds a CPYRAM13 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        pid = data[1]
        nids = [d if d > 0 else None for d in data[2:]]
        return CPYRAM13(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CPYRAM eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CPYRAM eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    @property
    def faces(self):
        """
        Gets the faces of the element

        Returns
        -------
        faces : Dict[int] = [face1, face2, ...]
            key = face number
            value = a list of nodes (integer pointers) as the values.

        .. note::  The order of the nodes are consistent with normals that point outwards
                   The face numbering is meaningless

        .. note::  The order of the nodes are consistent with ANSYS numbering; is this current?
        .. warning:: higher order element ids not verified with ANSYS; is this current?

        Examples
        --------
        >>> print(element.faces)
        """
        node_ids = self.node_ids
        n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13 = node_ids
        faces = {
            1 : [n1, n2, n3, n4, n6, n7, n8, n9],
            2 : [n1, n2, n5, n6, n11, n10],
            3 : [n2, n3, n5, n7, n12, n11],
            4 : [n3, n4, n5, n8, n13, n12],
            5 : [n4, n1, n5, n9, n10, n13],
        }
        return faces

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

    def _verify(self, xref):
        eid = self.eid
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

        V = (l * w) * h / 3
        V = A * h / 3
        """
        (n1, n2, n3, n4, n5,
         n6, n7, n8, n9, n10,
         n11, n12, n13) = self.get_node_positions()
        area1, c1 = area_centroid(n1, n2, n3, n4)
        volume = area1 / 3. * norm(c1 - n5)
        return abs(volume)

    @property
    def node_ids(self):
        nids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)
        return nids

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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


class CTETRA4(SolidElement):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CTETRA | EID | PID | G1 | G2 | G3 | G4 |
    +--------+-----+-----+----+----+----+----+
    """
    type = 'CTETRA'
    @property
    def faces(self):
        """
        Gets the faces of the element

        Returns
        -------
        faces : Dict[int] = [face1, face2, ...]
            key = face number
            value = a list of nodes (integer pointers) as the values.

        .. note::  The order of the nodes are consistent with normals that point outwards
                   The face numbering is meaningless

        Examples
        --------
        >>> print(element.faces)
        """
        nodes = self.node_ids
        faces = {
            1 : [nodes[0], nodes[1], nodes[3]],
            2 : [nodes[0], nodes[3], nodes[2]],
            3 : [nodes[1], nodes[2], nodes[3]],
            4 : [nodes[0], nodes[2], nodes[1]],
        }
        return faces

    @property
    def ansys_faces(self):
        """
        Gets the faces of the element

        Returns
        -------
        faces : Dict[int] = [face1, face2, ...]
            key = face number
            value = a list of nodes (integer pointers) as the values.

        .. note::  The order of the nodes are consistent with ANSYS numbering.
        .. warning:: higher order element ids not verified with ANSYS.

        Examples
        --------
        >>> print(element.faces)
        """
        nodes = self.node_ids
        faces = {
            1 : [nodes[0], nodes[1], nodes[2]],
            2 : [nodes[0], nodes[1], nodes[3]],
            3 : [nodes[1], nodes[2], nodes[3]],
            4 : [nodes[2], nodes[0], nodes[3]],
        }
        return faces

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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
        """
        Creates a CTETRA4

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : List[int]
            node ids; n=4
        comment : str; default=''
            a comment for the card
        """
        SolidElement.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 4

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CTETRA4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        """
        Adds a CTETRA4 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        pid = data[1]
        nids = data[2:]
        assert len(data) == 6, 'len(data)=%s data=%s' % (len(data), data)
        return CTETRA4(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CTETRA eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CTETRA eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def material_coordinate_system(self, xyz=None):
        """
        Returns
        -------
        centroid: (3,) float ndarray
           the centoid
        xe, ye, ze: (3,) float ndarray
            the element coordinate system

        """
        centroid, xe, ye, ze = _ctetra_element_coordinate_system(self, xyz=None)
        return centroid, xe, ye, ze

    def _verify(self, xref):
        eid = self.eid
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
        """Calculate the volume of the tet"""
        (n1, n2, n3, n4) = self.get_node_positions()
        return volume4(n1, n2, n3, n4)

    def Centroid(self):
        (n1, n2, n3, n4) = self.get_node_positions()
        return (n1 + n2 + n3 + n4) / 4.

    def get_face_nodes(self, nid_opposite, nid=None):
        assert nid is None, nid
        nids = self.node_ids[:4]
        indx = nids.index(nid_opposite)
        nids.pop(indx)
        return nids

    def get_face(self, nid_opposite, nid):
        nids = self.node_ids[:6]
        return ctetra_face(nid_opposite, nid, nids)

    def get_face_area_centroid_normal(self, nid, nid_opposite):
        return ctetra_face_area_centroid_normal(nid, nid_opposite,
                                                self.node_ids, self.nodes_ref)

    @property
    def node_ids(self):
        nids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)
        return nids

    def flip_normal(self):  ## TODO verify
        """flips the element inside out"""
        # flip n2 with n3
        n1, n2, n3, n4 = self.nodes
        self.nodes = [n1, n3, n2, n4]
        if self.nodes_ref is not None:
            n1_ref, n2_ref, n3_ref, n4_ref = self.nodes_ref
            self.nodes_ref = [n1_ref, n3_ref, n2_ref, n4_ref]

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
    | CTETRA | 1   | 1   | 239 | 229 | 516 | 99 | 335 | 103 |
    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    |        | 265 | 334 | 101 | 102 |     |    |     |     |
    +--------+-----+-----+-----+-----+-----+----+-----+-----+
    """
    type = 'CTETRA'
    def write_card(self, size: int=8, is_double: bool=False) -> str:
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

    def get_face_area_centroid_normal(self, nid_opposite, nid=None):
        return ctetra_face_area_centroid_normal(nid_opposite, nid,
                                                self.node_ids[:4], self.nodes_ref[:4])

    def __init__(self, eid, pid, nids, comment=''):
        """
        Creates a CTETRA10

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSOLID, PLSOLID)
        nids : List[int]
            node ids; n=10
        """
        SolidElement.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        nnodes = len(nids)
        if nnodes < 10:
            nids.extend((10 - nnodes) * [None])
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 10

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CTETRA10 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        """
        Adds a CTETRA10 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        pid = data[1]
        nids = [d if d > 0 else None for d in data[2:]]
        assert len(data) == 12, 'len(data)=%s data=%s' % (len(data), data)
        return CTETRA10(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CTETRA eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CTETRA eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def material_coordinate_system(self, xyz=None):
        """
        Returns
        -------
        centroid: (3,) float ndarray
           the centoid
        xe, ye, ze: (3,) float ndarray
            the element coordinate system

        """
        centroid, xe, ye, ze = _ctetra_element_coordinate_system(self, xyz=None)
        return centroid, xe, ye, ze

    @property
    def faces(self):
        """
        Gets the faces of the element

        Returns
        -------
        faces : Dict[int] = [face1, face2, ...]
            key = face number
            value = a list of nodes (integer pointers) as the values.

        .. note::  The order of the nodes are consistent with normals that point outwards
                   The face numbering is meaningless

        .. note::  The order of the nodes are consistent with ANSYS numbering; is this current?
        .. warning:: higher order element ids not verified with ANSYS; is this current?

        Examples
        --------
        >>> print(element.faces)
        """
        n1, n2, n3, n4, n5, n6, n7, n8, n9, n10 = self.node_ids
        faces = {
            1 : [n1, n2, n3, n5, n6, n7],  #More?
            2 : [n1, n2, n4, n5, n9, n8],
            3 : [n2, n3, n4, n6, n10, n9],
            4 : [n3, n1, n4, n7, n8, n10],
        }
        return faces

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

    def _verify(self, xref):
        eid = self.eid
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

    def get_face_nodes(self, nid_opposite, nid=None):
        nids = self.node_ids[:4]
        indx = nids.index(nid_opposite)
        nids.pop(indx)
        return nids

    @property
    def node_ids(self):
        nids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)
        return nids

def _ctetra_element_coordinate_system(element: Union[CTETRA4, CTETRA10], xyz=None):
    """
    Returns
    -------
    centroid: (3,) float ndarray
       the centoid
    xe, ye, ze: (3,) float ndarray
        the element coordinate system

    http://www.ipes.dk/Files/Ipes/Filer/nastran_2016_doc_release.pdf"""
    # this is the
    #if normal is None:
        #normal = element.Normal() # k = kmat

    if xyz is None:
        x1 = element.nodes_ref[0].get_position()
        x2 = element.nodes_ref[1].get_position()
        x3 = element.nodes_ref[2].get_position()
        x4 = element.nodes_ref[3].get_position()
    else:
        x1 = xyz[:, 0]
        x2 = xyz[:, 1]
        x3 = xyz[:, 2]
        x4 = xyz[:, 3]

    #CORDM=-2
    centroid = (x1 + x2 + x3 + x4) / 4.
    xe = (x2 + x3 + x4) / 3. - x1
    xe /= np.linalg.norm(xe)
    v = ((x1 + x3 + x4) - (x1 + x2 + x4)) / 3.
    ze = np.cross(xe, v)
    ze /= np.linalg.norm(ze)

    ye = np.cross(ze, xe)
    ye /= np.linalg.norm(ye)
    return centroid, xe, ye, ze

