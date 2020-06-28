## pylint: disable=C0103
"""
All shell elements are defined in this file.  This includes:

 * CTRIA3
 * CTRIA6

 * CSHEAR

 * CQUAD
 * CQUAD4
 * CQUAD8
 * CQUADR

 * CPLTSN3
 * CPLSTN4
 * CPLSTN6
 * CPLSTN8
 * SNORM

All tris are TriShell, ShellElement, and Element objects.
All quads are QuadShell, ShellElement, and Element objects.

"""
from __future__ import annotations
from typing import List, Tuple, Union, Any, TYPE_CHECKING

import numpy as np
from numpy import cross, allclose
from numpy.linalg import norm  # type: ignore

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_float_8
from pyNastran.bdf.cards.base_card import Element, BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, integer_double_or_blank, blank)
from pyNastran.bdf.field_writer_8 import print_card_8, print_field_8
from pyNastran.bdf.field_writer_16 import print_card_16, print_field_16
from pyNastran.bdf.cards.utils import wipe_empty_fields
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.nptyping import NDArray3float

__all__ = ['CTRIA3', 'CTRIA6', 'CSHEAR',
           'CQUAD', 'CQUAD4', 'CQUAD8', 'CQUADR',
           'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8',
           'CPLSTS3', 'CPLSTS4', 'CPLSTS6', 'CPLSTS8',
           '_triangle_area_centroid_normal', '_normal',
           'transform_shell_material_coordinate_system']

def _triangle_area_centroid_normal(nodes, card):
    """
    Gets the area, centroid and normal for a triangle.

    Parameters
    ----------
    nodes : List[np.ndarray]
        List of three triangle vertices.

    Returns
    -------
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
      unit_normal = n/|n|

    """
    (n1, n2, n3) = nodes
    vector = cross(n1 - n2, n1 - n3)
    length = norm(vector)
    try:
        normal = vector / length
    except FloatingPointError as error:
        #CTRIA3     20152     701   20174   20175   20176    8020
        #vector: [ 0.  0.  0.]; length: 0.0
            #[207.42750549, 0.0, -0.22425441]
            #[207.42750549, 0.0, 0.00631836]
            #[207.42750549, 0.0, 0.69803673]
        msg = error.message # strerror
        msg += '\nvector: %s; length: %s' % (vector, length)
        msg += '\n  %s\n  %s\n  %s' % (n1.tolist(), n2.tolist(), n3.tolist())
        raise RuntimeError(msg)

    if not allclose(norm(normal), 1.):
        msg = ('function _triangle_area_centroid_normal, check...\n'
               'a = {0}\nb = {1}\nnormal = {2}\nlength = {3}\n{4}'.format(
                   n1 - n2, n1 - n3, normal, length, str(card)))
        raise RuntimeError(msg)
    return (0.5 * length, (n1 + n2 + n3) / 3., normal)


def _normal(a, b):
    """Finds the unit normal vector of 2 vectors"""
    vector = cross(a, b)
    normal = vector / norm(vector)
    if not allclose(norm(normal), 1.):
        msg = ('function _normal, check...\n'
               'a = {0}\nb = {1}\nnormal = {2}\n'.format(
                   a, b, normal))
        raise RuntimeError(msg)
    return normal

def _normal4(n1, n2, n3, n4, card):
    """Finds the unit normal vector of 2 vectors"""
    a = n1 - n3
    b = n2 - n4
    vector = cross(a, b)
    normal = vector / norm(vector)
    if not allclose(norm(normal), 1.):
        msg = ('function _normal4, check...\n'
               f'a = {a}\nb = {b}\nnormal = {normal}\n{card}')
        raise RuntimeError(msg)
    return normal


class ShellElement(Element):
    type = 'ShellElement'

    def __init__(self):
        Element.__init__(self)

    def Theta_mcid(self):
        # () -> int
        if self.theta_mcid_ref is None:
            return self.theta_mcid
        return self.theta_mcid_ref.cid

    def _get_theta_mcid_repr(self):
        # : () -> str
        """
        set_blank_if_default doesn't distinguish between 0 and 0.0,
        so we fix it

        """
        theta_mcid = self.Theta_mcid()
        if isinstance(theta_mcid, float):
            theta_mcid = set_blank_if_default(theta_mcid, 0.0)
        return theta_mcid

    def Thickness(self) -> float:
        """Returns the thickness"""
        return self.pid_ref.Thickness()

    #def Volume(self):
        #"""
        #Returns the volume
        #"""
        #return self.Thickness() * self.Area()

    @property
    def material_ids(self) -> List[int]:
        """
        Returns the material

        .. todo:: possibly remove this

        """
        return self.pid_ref.material_ids

    def mid(self) -> Any:
        """
        Returns the material

        .. todo:: possibly remove this

        """
        return self.pid_ref.mid()

    def Mid(self) -> int:
        """
        Returns the material ID

        .. todo:: possibly remove this

        """
        return self.pid_ref.Mid()

    def Nsm(self) -> float:
        """Returns the non-structural mass"""
        return self.pid_ref.Nsm()

    def MassPerArea(self) -> float:
        """Returns the mass per area"""
        tscales = self.get_thickness_scale()
        return self.pid_ref.MassPerArea(tflag=self.tflag, tscales=tscales)

    def Mass(self) -> float:
        r"""
        .. math:: m = \frac{m}{A} A  \f]

        """
        A = self.Area()
        tscales = self.get_thickness_scale()
        try:
            mpa = self.pid_ref.MassPerArea(tflag=self.tflag, tscales=tscales)
        except TypeError:
            print(self.pid_ref)
            raise

        try:
            return mpa * A
        except TypeError:
            msg = 'mass/area=%s area=%s prop_type=%s' % (mpa, A, self.pid_ref.type)
            raise TypeError(msg)

    def Mass_breakdown(self) -> float:
        r"""
        .. math:: m = \frac{m}{A} A  \f]

        """
        A = self.Area()
        mpa = self.pid_ref.MassPerArea_structure()
        nsm = self.pid_ref.nsm
        try:
            return mpa * A, nsm * A
        except TypeError:
            msg = 'mass/area=%s area=%s prop_type=%s' % (mpa, A, self.pid_ref.type)
            raise TypeError(msg)

    def Mass_no_xref(self, model):
        r"""
        .. math:: m = \frac{m}{A} A  \f]

        """
        A = self.Area_no_xref(model)
        pid_ref = model.Property(self.pid)
        tavg_scale = self.get_thickness_scale()
        mpa = pid_ref.MassPerArea_no_xref(model, tavg_scale)
        try:
            return mpa * A
        except TypeError:
            msg = 'mass/area=%s area=%s prop_type=%s' % (mpa, A, self.pid_ref.type)
            raise TypeError(msg)

    #def flip_normal(self):
        #raise NotImplementedError('flip_normal undefined for %s' % self.type)


class TriShell(ShellElement):
    def __init__(self):
        ShellElement.__init__(self)
        self.nodes_ref = None  # type: Optional[List[Any]]
        self.pid_ref = None  # type: Optional[Any]

    def get_edge_ids(self):
        """Return the edge IDs"""
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


    def Thickness(self) -> float:
        """Returns the thickness"""
        tscales = self.get_thickness_scale()
        return self.pid_ref.Thickness(tflag=self.tflag, tscales=tscales)
        #return self.pid_ref.Thickness()

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
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
        return _triangle_area_centroid_normal([n1, n2, n3], self)

    def get_area(self) -> float:
        """see ``TriShell.Area()``"""
        return self.Area()

    def Area(self) -> float:
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} \lvert (n_0-n_1) \times (n_0-n_2) \rvert"""
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
        a = n1 - n2
        b = n1 - n3
        area = 0.5 * norm(cross(a, b))
        return area

    def Normal(self):
        # () -> np.ndarray
        r"""
        Get the normal vector, :math:`n`.

        .. math::
          n = \frac{(n_0-n_1) \times (n_0-n_2)}
             {\lvert (n_0-n_1) \times (n_0-n_2) \lvert}

        """
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
        try:
            n = _normal(n1 - n2, n1 - n3)
        except:
            msg = 'ERROR computing normal vector for eid=%i.\n' % self.eid
            msg += '  nid1=%i n1=%s\n' % (self.nodes_ref[0].nid, n1)
            msg += '  nid2=%i n2=%s\n' % (self.nodes_ref[1].nid, n2)
            msg += '  nid3=%i n3=%s\n' % (self.nodes_ref[2].nid, n3)
            raise RuntimeError(msg)
        return n

    def Centroid(self):
        # () -> np.ndarray
        r"""
        Get the centroid.

        .. math::
          CG = \frac{1}{3} (n_0+n_1+n_2)

        """
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
        centroid = (n1 + n2 + n3) / 3.
        return centroid

    def center_of_mass(self):
        # () -> np.ndarray
        return self.Centroid()

    #def uncross_reference(self):
        #self.nodes = self.node_ids
        #self.pid = self.Pid()
        #del self.nodes_ref, self.pid_ref

    def material_coordinate_system(self, normal=None, xyz123=None):
        """
        Determines the material coordinate system

        Parameters
        ----------
        normal (3, ) float ndarray
            the unit normal vector
        xyz123 (3, 3) float ndarray
            the xyz coordinates

        Returns
        -------
        centroid (3, ) float ndarray
            the centroid of the element
        imat (3, ) float ndarray
            the element unit i vector
        jmat (3, ) float ndarray
            the element unit j vector
        normal (3, ) float ndarray
            the unit normal vector

        .. todo:: rotate the coordinate system by the angle theta

        """
        if normal is None:
            normal = self.Normal() # k = kmat

        if xyz123 is None:
            xyz1 = self.nodes_ref[0].get_position()
            xyz2 = self.nodes_ref[1].get_position()
            xyz3 = self.nodes_ref[2].get_position()
            #centroid = (xyz1 + xyz2 + xyz3) / 3.
            #centroid = self.Centroid()
        else:
            #centroid = xyz1234.sum(axis=1)
            #assert len(centroid) == 3, centroid
            xyz1 = xyz123[:, 0]
            xyz2 = xyz123[:, 1]
            xyz3 = xyz123[:, 2]
        centroid = (xyz1 + xyz2 + xyz3) / 3.
        imat, jmat = _material_coordinate_system(self, normal, xyz1, xyz2)
        return centroid, imat, jmat, normal


def _material_coordinate_system(element, normal, xyz1, xyz2):
    """helper function for material_coordinate_system"""
    if element.theta_mcid is None:
        raise NotImplementedError('theta_mcid=%r' % element.theta_mcid)
    if isinstance(element.theta_mcid, integer_types):
        assert element.theta_mcid_ref is not None, f'mcid={element.theta_mcid} not found for\n{element}'
        i = element.theta_mcid_ref.i
        jmat = np.cross(normal, i) # k x i
        try:
            jmat /= np.linalg.norm(jmat)
        except FloatingPointError:
            raise ValueError(f'Cannot project i-axis onto element normal i={i} normal={normal}\n{element}')
        # we do an extra normalization here because
        # we had to project i onto the elemental plane
        # unlike in the next block
        imat = np.cross(jmat, normal)
    elif isinstance(element.theta_mcid, float):
        # TODO: rotate by the angle theta
        imat = xyz2 - xyz1
        imat /= np.linalg.norm(imat)
        jmat = np.cross(normal, imat) # k x i
        try:
            jmat /= np.linalg.norm(jmat)
        except FloatingPointError:
            raise ValueError(f'Cannot project i-axis onto element normal i={i} normal={normal}\n{element}')
    else:
        raise RuntimeError(element.theta_mcid)
    return imat, jmat


class CTRIA3(TriShell):
    """
    +--------+-------+-------+----+----+----+------------+---------+
    |   1    |   2   |   3   |  4 |  5 |  6 |     7      |    8    |
    +========+=======+=======+=====+===+====+============+=========+
    | CTRIA3 |  EID  |  PID  | N1 | N2 | N3 | THETA/MCID | ZOFFSET |
    +--------+-------+-------+----+----+----+------------+---------+
    |        |       | TFLAG | T1 | T2 | T3 |            |         |
    +--------+-------+-------+----+----+----+------------+---------+

    """
    type = 'CTRIA3'
    _field_map = {
        1: 'eid', 2:'pid', 6:'theta_mcid', 7:'zoffset', 10:'tflag',
        11:'T1', 12:'T2', 13:'T3'}
    cp_name_map = {
        'T1' : 'T1',
        'T2' : 'T2',
        'T3' : 'T3',
    }
    _properties = ['cp_name_map', '_field_map']

    def _update_field_helper(self, n, value):
        # (int, Any) -> None
        if n == 3:
            self.nodes[0] = value
        elif n == 4:
            self.nodes[1] = value
        elif n == 5:
            self.nodes[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        mcids = []
        thetas = []
        zoffsets = []
        #t123 = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
            if isinstance(element.theta_mcid, int):
                mcid = element.theta_mcid
                theta = 0.
            else:
                assert isinstance(element.theta_mcid, float), type(element.theta_mcid)
                mcid = -1
                theta = element.theta_mcid
            mcids.append(mcid)
            thetas.append(theta)
            zoffsets.append(element.zoffset)
            #t123.append([element.T1, element.T2, element.T3])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('mcid', data=mcids)
        h5_file.create_dataset('theta', data=thetas)
        h5_file.create_dataset('zoffset', data=zoffsets)
        #self.tflag = tflag

    def __init__(self, eid, pid, nids, zoffset=0., theta_mcid=0.0,
                 tflag=0, T1=None, T2=None, T3=None, comment=''):
        # (int, int, List[int], float, Union[int, float], int, Optional[float], Optional[float], Optional[float], str) -> None
        """
        Creates a CTRIA3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 : float; default=None
            If it is not supplied, then T1 through T3 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        TriShell.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        assert len(nids) == 3, nids
        self.nodes = self.prepare_node_ids(nids)
        self.zoffset = zoffset
        self.theta_mcid = theta_mcid
        self.tflag = tflag
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        assert len(self.nodes) == 3
        self.theta_mcid_ref = None  # type: Optional[Any]

    def validate(self):
        assert len(set(self.nodes)) == 3, 'nodes=%s; n=%s\n%s' % (self.nodes, len(set(self.nodes)), str(self))

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CTRIA3 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        pid = data[1]
        nids = data[2:5]

        theta_mcid = data[5]
        zoffset = data[6]
        tflag = data[7]
        T1 = data[8]
        T2 = data[9]
        T3 = data[10]
        if T1 == -1.0:
            T1 = 1.0
        if T2 == -1.0:
            T2 = 1.0
        if T3 == -1.0:
            T3 = 1.0
        assert tflag in [0, 1], data
        return CTRIA3(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, comment=comment)

    @classmethod
    def add_card(cls, card, comment=''):
        # type: (Any, str) -> CTRIA3
        """
        Adds a CTRIA3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
        ]
        if len(card) > 6:
            theta_mcid = integer_double_or_blank(card, 6, 'theta_mcid', 0.0)
            zoffset = double_or_blank(card, 7, 'zoffset', 0.0)
            blank(card, 8, 'blank')
            blank(card, 9, 'blank')

            tflag = integer_or_blank(card, 10, 'tflag', 0)
            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            assert len(card) <= 14, f'len(CTRIA3 card) = {len(card):d}\ncard={card}\n tflag={tflag} T123=[{T1}, {T2}, {T3}]'
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            tflag = 0
            T1 = None
            T2 = None
            T3 = None
        return CTRIA3(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CTRIA3 eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.Pid(), msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.Coord(self.theta_mcid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CTRIA3 eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.safe_coord(self.theta_mcid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.theta_mcid = self.Theta_mcid()
        self.pid_ref = None
        self.nodes_ref = None
        self.theta_mcid_ref = None

    def _verify(self, xref: bool) -> None:
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        unused_edges = self.get_edge_ids()

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
        #    if self.tflag == 0:
        #        t = self.pid_ref.Thickness()
        #        t1 = self.T1 / t
        #        t2 = self.T2 / t
        #        t3 = self.T3 / t
        #    else:
        #        t1 = self.T1
        #        t2 = self.T2
        #        t3 = self.T3
        #    t = (t1 + t2 + t3)/3.
        #else:
        #    t = self.pid_ref.Thickness()
        #return t

    def get_thickness_scale(self):
        return [self.T1, self.T2, self.T3]

    def flip_normal(self):
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
        if self.nodes_ref is not None:
            (n1, n2, n3) = self.nodes_ref
            self.nodes_ref = [n1, n3, n2]

    def _get_repr_defaults(self):
        zoffset = set_blank_if_default(self.zoffset, 0.0)
        tflag = set_blank_if_default(self.tflag, 0)
        theta_mcid = self._get_theta_mcid_repr()
        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        return theta_mcid, zoffset, tflag, T1, T2, T3

    @property
    def node_ids(self) -> List[int]:
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)

    def raw_fields(self):
        list_fields = (['CTRIA3', self.eid, self.Pid()] + self.node_ids +
                       [self.theta_mcid, self.zoffset, None] +
                       [None, self.tflag, self.T1, self.T2, self.T3])
        return list_fields

    def repr_fields(self):
        (theta_mcid, zoffset, tflag, T1, T2, T3) = self._get_repr_defaults()
        list_fields = (['CTRIA3', self.eid, self.Pid()] + self.node_ids +
                       [theta_mcid, zoffset, None] + [None, tflag, T1, T2, T3])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        zoffset = set_blank_if_default(self.zoffset, 0.0)
        tflag = set_blank_if_default(self.tflag, 0)
        theta_mcid = self._get_theta_mcid_repr()

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)

        #return self.write_card(size, double)
        nodes = self.node_ids
        row2_data = [theta_mcid, zoffset,
                     tflag, T1, T2, T3]
        row2 = [print_field_8(field) for field in row2_data]
        data = [self.eid, self.Pid()] + nodes + row2
        msg = ('CTRIA3  %8i%8i%8i%8i%8i%8s%8s\n'
               '                %8s%8s%8s%8s\n' % tuple(data))
        return self.comment + msg.rstrip() + '\n'

class CPLSTx3(TriShell):
    """
    +---------+-------+-------+----+----+----+-------+
    |    1    |   2   |   3   |  4 |  5 |  6 |   7   |
    +=========+=======+=======+=====+===+====+=======+
    | CPLSTN3 |  EID  |  PID  | N1 | N2 | N3 | THETA |
    +---------+-------+-------+----+----+----+-------+

    """
    type = 'CPLSTN3'
    _field_map = {1: 'eid', 2:'pid', 6:'theta', }

    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 4:
            self.nodes[1] = value
        elif n == 5:
            self.nodes[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, pid, nids, theta=0.0, comment=''):
        """NX specific card"""
        TriShell.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        assert len(nids) == 3, nids
        self.nodes = self.prepare_node_ids(nids)
        self.theta = theta
        assert len(self.nodes) == 3

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        thetas = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
            thetas.append(element.theta)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('theta', data=thetas)

    def validate(self):
        assert len(set(self.nodes)) == 3, 'nodes=%s; n=%s\n%s' % (self.nodes, len(set(self.nodes)), str(self))

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #eid = data[0]
        #pid = data[1]
        #nids = data[2:5]
        #theta = data[5]
        #return CPLSTN3(eid, pid, nids, theta, comment=comment)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CPLSTx3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
        ]
        if len(card) > 5:
            theta = double_or_blank(card, 6, 'theta', 0.0)
            assert len(card) <= 14, 'len(CPLSTN3 card) = %i\ncard=%s' % (len(card), card)
        else:
            theta = 0.0
        return cls(eid, pid, nids, theta, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid}'
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.Pid(), msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid}'
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        unused_edges = self.get_edge_ids()

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PPLANE', 'PLPLANE', 'PGPLSN'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)

            #if not self.pid_ref.type in ['PLPLANE']:
            #t = self.Thickness()
            #assert isinstance(t, float), 'thickness=%r' % t
            #mass = self.Mass()
            #assert isinstance(mass, float), 'mass=%r' % mass
            #a, c, n = self.AreaCentroidNormal()
            #assert isinstance(a, float), 'Area=%r' % a
            #for i in range(3):
                #assert isinstance(c[i], float)
                #assert isinstance(n[i], float)

    def flip_normal(self):
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
        if self.nodes_ref is not None:
            (n1, n2, n3) = self.nodes_ref
            self.nodes_ref = [n1, n3, n2]

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)

    def raw_fields(self):
        list_fields = ([self.type, self.eid, self.Pid()] + self.node_ids +
                       [self.theta])
        return list_fields

    def repr_fields(self):
        theta = set_blank_if_default(self.theta, 0.0)
        list_fields = ([self.type, self.eid, self.Pid()] + self.node_ids +
                       [theta])
        return list_fields


class CPLSTN3(CPLSTx3):
    type = 'CPLSTN3'

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes + [self.theta]
        msg = ('CPLSTN3 %8i%8i%8i%8i%8i%8s\n' % tuple(data))
        return self.comment + msg

class CPLSTS3(CPLSTx3):
    type = 'CPLSTS3'

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes + [self.theta]
        msg = ('CPLSTS3 %8i%8i%8i%8i%8i%8s\n' % tuple(data))
        return self.comment + msg

class CTRIA6(TriShell):
    """
    +--------+------------+---------+----+----+----+----+----+-----+
    |   1    |      2     |    3    |  4 |  5 |  6 | 7  | 8  |  9  |
    +========+============+=========+=====+===+====+====+====+=====+
    | CTRIA3 |    EID     |   PID   | N1 | N2 | N3 | N4 | N5 | N6  |
    +--------+------------+---------+----+----+----+----+----+-----+
    |        | THETA/MCID | ZOFFSET | T1 | T2 | T3 |    |    |     |
    +--------+------------+---------+----+----+----+----+----+-----+

    """
    type = 'CTRIA6'
    def __init__(self, eid, pid, nids, theta_mcid=0., zoffset=0., tflag=0,
                 T1=None, T2=None, T3=None, comment=''):
        """
        Creates a CTRIA6 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int/None, int/None, int/None]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 : float; default=None
            If it is not supplied, then T1 through T3 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        TriShell.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.theta_mcid = theta_mcid
        self.zoffset = zoffset
        self.tflag = tflag
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(nids) == 6, 'error on CTRIA6'
        self.theta_mcid_ref = None  # type: Optional[Any]

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        mcids = []
        thetas = []
        zoffsets = []
        #t123 = []

        #element0 = model.elements[eids[0]]
        neids = len(eids)
        nnodes = 6
        nodes = np.zeros((neids, nnodes), dtype='int32')

        for i, eid in enumerate(eids):
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes[i, :] = [nid if nid is not None else 0 for nid in element.nodes]
            if isinstance(element.theta_mcid, int):
                mcid = element.theta_mcid
                theta = 0.
            else:
                assert isinstance(element.theta_mcid, float), type(element.theta_mcid)
                mcid = -1
                theta = element.theta_mcid
            mcids.append(mcid)
            thetas.append(theta)
            zoffsets.append(element.zoffset)
            #t123.append([element.T1, element.T2, element.T3])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('mcid', data=mcids)
        h5_file.create_dataset('theta', data=thetas)
        h5_file.create_dataset('zoffset', data=zoffsets)
        #self.tflag = tflag

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CTRIA6 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
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
            integer_or_blank(card, 8, 'n6', 0),
        ]
        if len(card) > 9:
            theta_mcid = integer_double_or_blank(card, 9, 'theta_mcid', 0.0)
            zoffset = double_or_blank(card, 10, 'zoffset', 0.0)

            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            tflag = integer_or_blank(card, 14, 'tflag', 0)
            assert len(card) <= 15, 'len(CTRIA6 card) = %i\ncard=%s' % (len(card), card)
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            T1 = None
            T2 = None
            T3 = None
            tflag = 0
        return CTRIA6(eid, pid, nids, theta_mcid, zoffset,
                      tflag, T1, T2, T3, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CTRIA6 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        pid = data[1]
        nids = data[2:8]
        theta_mcid = data[8]
        zoffset = data[9]
        t1 = data[10]
        t2 = data[11]
        t3 = data[12]
        assert isinstance(t1, float), data
        assert isinstance(t2, float), data
        assert isinstance(t3, float), data

        if t1 == -1.:
            t1 = None
        if t2 == -1.:
            t2 = None
        if t3 == -1.:
            t3 = None
        tflag = data[13]
        assert isinstance(tflag, integer_types), data
        #prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(nids) == 6, 'error on CTRIA6'
        assert tflag in [0, 1], data
        return CTRIA6(eid, pid, nids, theta_mcid, zoffset,
                      tflag, t1, t2, t3, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CTRIA6 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.Pid(), msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.Coord(self.theta_mcid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CTRIA6 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.safe_coord(self.theta_mcid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.theta_mcid = self.Theta_mcid()
        self.pid_ref = None
        self.nodes_ref = None
        self.theta_mcid_ref = None

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        unused_edges = self.get_edge_ids()

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

    def get_thickness_scale(self):
        return [self.T1, self.T2, self.T3]

    def Thickness(self):
        """Returns the thickness, :math:`t`"""
        return self.pid_ref.Thickness()

    def AreaCentroidNormal(self):
        """
        Returns area, centroid, normal as it's more efficient to do them
        together

        """
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
        return _triangle_area_centroid_normal([n1, n2, n3], self)

    def Area(self):
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} (n_0-n_1) \times (n_0-n_2)"""
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
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
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
        return _normal(n1 - n2, n1 - n3)

    def Centroid(self):
        r"""
        Get the centroid.

        .. math::
          CG = \frac{1}{3} (n_1+n_2+n_3)

        """
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
        centroid = (n1 + n2 + n3) / 3.
        return centroid

    def center_of_mass(self):
        return self.Centroid()

    def flip_normal(self):
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
        if self.nodes_ref is not None:
            (n1, n2, n3, n4, n5, n6) = self.nodes_ref
            self.nodes_ref = [n1, n3, n2, n6, n5, n4]

    def _get_repr_defaults(self):
        zoffset = set_blank_if_default(self.zoffset, 0.0)
        assert isinstance(self.tflag, integer_types), self.tflag
        tflag = set_blank_if_default(self.tflag, 0)
        theta_mcid = self._get_theta_mcid_repr()

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        return theta_mcid, zoffset, tflag, T1, T2, T3

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = (['CTRIA6', self.eid, self.Pid()] + self.node_ids +
                       [self.theta_mcid, self.zoffset,
                        self.T1, self.T2, self.T3, self.tflag,])
        return list_fields

    def repr_fields(self):
        (theta_mcid, zoffset, tflag, T1, T2, T3) = self._get_repr_defaults()
        list_fields = (['CTRIA6', self.eid, self.Pid()] + self.node_ids +
                       [theta_mcid, zoffset, T1, T2, T3, tflag])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 8: # to last node
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg


class CTRIAR(TriShell):
    """
    +--------+-------+-------+----+----+----+------------+---------+
    |   1    |   2   |   3   |  4 |  5 |  6 |     7      |    8    |
    +========+=======+=======+=====+===+====+============+=========+
    | CTRIAR |  EID  |  PID  | N1 | N2 | N3 | THETA/MCID | ZOFFSET |
    +--------+-------+-------+----+----+----+------------+---------+
    |        |       | TFLAG | T1 | T2 | T3 |            |         |
    +--------+-------+-------+----+----+----+------------+---------+

    """
    type = 'CTRIAR'
    def __init__(self, eid, pid, nids, theta_mcid=0.0, zoffset=0.0,
                 tflag=0, T1=None, T2=None, T3=None, comment=''):
        """
        Creates a CTRIAR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 : float; default=None
            If it is not supplied, then T1 through T3 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        TriShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid

        self.theta_mcid = theta_mcid
        self.zoffset = zoffset
        self.tflag = tflag
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=False)
        assert len(self.nodes) == 3
        self.theta_mcid_ref = None  # type: Optional[Any]

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        mcids = []
        thetas = []
        zoffsets = []
        #t123 = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
            if isinstance(element.theta_mcid, int):
                mcid = element.theta_mcid
                theta = 0.
            else:
                assert isinstance(element.theta_mcid, float), type(element.theta_mcid)
                mcid = -1
                theta = element.theta_mcid
            mcids.append(mcid)
            thetas.append(theta)
            zoffsets.append(element.zoffset)
            #t123.append([element.T1, element.T2, element.T3])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('mcid', data=mcids)
        h5_file.create_dataset('theta', data=thetas)
        h5_file.create_dataset('zoffset', data=zoffsets)
        #self.tflag = tflag

    #def validate(self):
        #self.validate_node_ids(allow_empty_nodes=False)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CTRIAR card from ``BDF.add_card(...)``

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
        ]

        theta_mcid = integer_double_or_blank(card, 6, 'theta_mcid', 0.0)
        zoffset = double_or_blank(card, 7, 'zoffset', 0.0)
        blank(card, 8, 'blank')
        blank(card, 9, 'blank')

        tflag = integer_or_blank(card, 10, 'tflag', 0)
        T1 = double_or_blank(card, 11, 'T1')
        T2 = double_or_blank(card, 12, 'T2')
        T3 = double_or_blank(card, 13, 'T3')
        assert len(card) <= 14, 'len(CTRIAR card) = %i\ncard=%s' % (len(card), card)
        return CTRIAR(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CTRIAR card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        pid = data[1]
        nids = data[2:5]

        theta_mcid = data[5]
        zoffset = data[6]
        tflag = data[7]
        T1 = data[8]
        T2 = data[9]
        T3 = data[10]
        if T1 == -1.0:
            T1 = 1.0
        if T2 == -1.0:
            T2 = 1.0
        if T3 == -1.0:
            T3 = 1.0
        assert tflag in [0, 1], data
        return CTRIAR(eid, pid, nids, zoffset=zoffset, theta_mcid=theta_mcid,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CTRIAR eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.Coord(self.theta_mcid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CTRIAR eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.safe_coord(self.theta_mcid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.theta_mcid = self.Theta_mcid()
        self.pid_ref = None
        self.nodes_ref = None
        self.theta_mcid_ref = None

    def get_thickness_scale(self):
        return [self.T1, self.T2, self.T3]

    def Thickness(self):
        """Returns the thickness"""
        return self.pid_ref.Thickness()

    def flip_normal(self):
        r"""
        ::

               1           1
              * *   -->   * *
             *   *       *   *
            2-----3     3-----2

        """
        (n1, n2, n3) = self.nodes
        self.nodes = [n1, n3, n2]
        if self.nodes_ref is not None:
            (n1, n2, n3) = self.nodes_ref
            self.nodes_ref = [n1, n3, n2]

    def _get_repr_defaults(self):
        zoffset = set_blank_if_default(self.zoffset, 0.0)
        tflag = set_blank_if_default(self.tflag, 0)
        theta_mcid = self._get_theta_mcid_repr()

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        return (theta_mcid, zoffset, tflag, T1, T2, T3)


    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        unused_nids = self.node_ids
        unused_edges = self.get_edge_ids()

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            # PSHELL/PCOMP
            assert self.pid_ref.type in ['PSHELL', 'PCOMP', 'PCOMPG'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            t = self.Thickness()
            a, c, unused_normal = self.AreaCentroidNormal()
            assert isinstance(t, float), 'thickness=%r' % t
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                #assert isinstance(normal[i], float)
            mass = self.Mass()
            assert isinstance(mass, float), 'mass=%r' % mass

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)

    def raw_fields(self):
        list_fields = (['CTRIAR', self.eid, self.Pid()] + self.node_ids +
                       [self.theta_mcid, self.zoffset, self.tflag,
                        self.T1, self.T2, self.T3])
        return list_fields

    def repr_fields(self):
        (theta_mcid, zoffset, tflag, T1, T2, T3) = self._get_repr_defaults()
        list_fields = (['CTRIAR', self.eid, self.Pid()] + self.node_ids +
                       [theta_mcid, zoffset, None, None, tflag, T1, T2, T3])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 5: # to last node
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        return msg


class QuadShell(ShellElement):
    def __init__(self):
        ShellElement.__init__(self)
        self.nodes_ref = None  # type: Optional[List[Any]]
        self.pid_ref = None  # type: Optional[Any]

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

    def get_thickness_scale(self):
        return None

    def Thickness(self):
        """Returns the thickness"""
        tscales = self.get_thickness_scale()
        try:
            return self.pid_ref.Thickness(tflag=self.tflag, tscales=tscales)
        except TypeError:
            print(self.pid_ref)
            raise
        #return self.pid_ref.Thickness()

    def Normal(self):
        try:
            n1, n2, n3, n4 = self.get_node_positions(nodes=self.nodes_ref[:4])
        except ValueError:
            print(str(self))
            raise
        try:
            n = _normal(n1 - n3, n2 - n4)
        except:
            msg = 'ERROR computing normal vector for eid=%i.\n' % self.eid
            msg += '  nid1=%i n1=%s\n' % (self.nodes_ref[0].nid, n1)
            msg += '  nid2=%i n2=%s\n' % (self.nodes_ref[1].nid, n2)
            msg += '  nid3=%i n3=%s\n' % (self.nodes_ref[2].nid, n3)
            msg += '  nid4=%i n4=%s\n' % (self.nodes_ref[3].nid, n4)
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
        nodes_ref = self.nodes_ref[:4]
        n1, n2, n3, n4 = self.get_node_positions(nodes=nodes_ref)
        area = 0.5 * norm(cross(n3-n1, n4-n2))
        centroid = (n1 + n2 + n3 + n4) / 4.
        return(area, centroid)

    def Centroid(self):
        nodes_ref = self.nodes_ref[:4]
        n1, n2, n3, n4 = self.get_node_positions(nodes=nodes_ref)
        centroid = (n1 + n2 + n3 + n4) / 4.
        return centroid

    def Centroid_no_xref(self, model):
        nodes = self.nodes[:4]
        n1, n2, n3, n4 = self.get_node_positions_no_xref(model, nodes=nodes)
        centroid = (n1 + n2 + n3 + n4) / 4.
        return centroid

    def center_of_mass(self):
        return self.Centroid()

    def get_area(self):
        return self.Area()

    def Area(self):
        """
        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert
        where a and b are the quad's cross node point vectors"""
        (n1, n2, n3, n4) = self.get_node_positions(nodes=self.nodes_ref[:4])
        area = 0.5 * norm(cross(n3-n1, n4-n2))
        return area

    def Area_no_xref(self, model):
        """
        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert
        where a and b are the quad's cross node point vectors"""
        (n1, n2, n3, n4) = self.get_node_positions_no_xref(model, nodes=self.nodes[:4])
        area = 0.5 * norm(cross(n3-n1, n4-n2))
        return area

    def flip_normal(self):
        r"""
        ::

          1---2       1---4
          |   |  -->  |   |
          |   |       |   |
          4---3       2---3
        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]
        if self.nodes_ref is not None:
            (n1, n2, n3, n4) = self.nodes_ref
            self.nodes_ref = [n1, n4, n3, n2]

    def _get_repr_defaults(self):
        zoffset = set_blank_if_default(self.zoffset, 0.0)
        tflag = set_blank_if_default(self.tflag, 0)

        theta_mcid = self._get_theta_mcid_repr()
        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        T4 = set_blank_if_default(self.T4, 1.0)
        return (theta_mcid, zoffset, tflag, T1, T2, T3, T4)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def material_coordinate_system(self, normal=None, xyz1234=None):
        """
        Determines the material coordinate system

        Parameters
        ----------
        normal (3, ) float ndarray
            the unit normal vector
        xyz1234 (4, 3) float ndarray
            the xyz coordinates

        Returns
        -------
        centroid (3, ) float ndarray
            the centroid of the element
        imat (3, ) float ndarray
            the element unit i vector
        jmat (3, ) float ndarray
            the element unit j vector
        normal (3, ) float ndarray
            the unit normal vector

        .. todo:: rotate the coordinate system by the angle theta

        """
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

        imat, jmat = _material_coordinate_system(self, normal, xyz1, xyz2)
        return centroid, imat, jmat, normal


class CSHEAR(QuadShell):
    """
    +--------+-------+-------+----+----+----+----+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |
    +========+=======+=======+=====+===+====+====+
    | CSHEAR |  EID  |  PID  | N1 | N2 | N3 | N4 |
    +--------+-------+-------+----+----+----+----+

    """
    type = 'CSHEAR'
    def __init__(self, eid, pid, nids, comment=''):
        """
        Creates a CSHEAR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHEAR)
        nids : List[int, int, int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        QuadShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 4

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)

    def validate(self):
        assert len(set(self.nodes)) == 4, 'nodes=%s\n%s' % (self.nodes, str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CSHEAR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
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
        """
        Adds a CSHEAR card from the OP2

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
        return CSHEAR(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CSHEAR eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.Pid(), msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CSHEAR eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def Normal(self):
        (n1, n2, n3, n4) = self.get_node_positions()
        return _normal(n1 - n3, n2 - n4)

    def AreaCentroidNormal(self) -> Tuple[NDArray3float, NDArray3float, NDArray3float]:
        (area, centroid) = self.AreaCentroid()
        normal = self.Normal()
        return (area, centroid, normal)

    def AreaCentroid(self) -> Tuple[NDArray3float, NDArray3float]:
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

    def Centroid(self) -> NDArray3float:
        (n1, n2, n3, n4) = self.get_node_positions()
        centroid = (n1 + n2 + n3 + n4) / 4.
        return centroid

    def center_of_mass(self) -> NDArray3float:
        return self.Centroid()

    def _verify(self, xref: bool) -> None:
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        unused_edges = self.get_edge_ids()

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

    def Area(self) -> float:
        r"""
        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert
        where a and b are the quad's cross node point vectors"""
        (n1, n2, n3, n4) = self.get_node_positions()
        a = n1 - n3
        b = n2 - n4
        area = 0.5 * norm(cross(a, b))
        return area

    def flip_normal(self):
        r"""
        ::

          1---2       1---4
          |   |  -->  |   |
          |   |       |   |
          4---3       2---3

        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]
        if self.nodes_ref is not None:
            (n1, n2, n3, n4) = self.nodes_ref
            self.nodes_ref = [n1, n4, n3, n2]

    @property
    def node_ids(self) -> List[int]:
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self) -> List[Union[str, int]]:
        list_fields = ['CSHEAR', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def repr_fields(self) -> List[Union[str, int]]:
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        msg = self.comment + print_card_8(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg

    def G(self) -> float:
        return self.pid_ref.mid_ref.G()

    def Mass(self) -> float:
        r"""
        .. math:: m = \frac{m}{A} A  \f]
        """
        A = self.Area()
        mpa = self.pid_ref.MassPerArea()
        try:
            return mpa * A
        except TypeError:
            msg = 'mass/area=%s area=%s prop_type=%s' % (mpa, A, self.pid_ref.type)
            raise TypeError(msg)

    def Thickness(self) -> float:
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
    cp_name_map = {
        'T1' : 'T1',
        'T2' : 'T2',
        'T3' : 'T3',
        'T4' : 'T4',
    }
    _field_map = {1: 'eid', 2:'pid', 7:'theta_mcid', 8:'zoffset',
                  10:'tflag', 11:'T1', 12:'T2', 13:'T3'}
    _properties = ['cp_name_map', '_field_map']

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

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        mcids = []
        thetas = []
        zoffsets = []
        #t1234 = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
            if isinstance(element.theta_mcid, int):
                mcid = element.theta_mcid
                theta = 0.
            else:
                assert isinstance(element.theta_mcid, float), type(element.theta_mcid)
                mcid = -1
                theta = element.theta_mcid
            mcids.append(mcid)
            thetas.append(theta)
            zoffsets.append(element.zoffset)
            #t1234.append([element.T1, element.T2, element.T3, element.T4])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('mcid', data=mcids)
        h5_file.create_dataset('theta', data=thetas)
        h5_file.create_dataset('zoffset', data=zoffsets)
        #self.tflag = tflag

    def __init__(self, eid, pid, nids, theta_mcid=0.0, zoffset=0.,
                 tflag=0, T1=None, T2=None, T3=None, T4=None, comment=''):
        """
        Creates a CQUAD4 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 / T4 : float; default=None
            If it is not supplied, then T1 through T4 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        QuadShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        assert len(nids) == 4, nids
        self.nodes = self.prepare_node_ids(nids)
        self.zoffset = zoffset
        self.theta_mcid = theta_mcid
        self.tflag = tflag
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.T4 = T4
        self.theta_mcid_ref = None  # type: Optional[Any]

    def validate(self):
        assert len(set(self.nodes)) == 4, 'nodes=%s\n%s' % (self.nodes, str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CQUAD4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4'),]
        if len(card) > 7:
            theta_mcid = integer_double_or_blank(card, 7, 'theta_mcid', 0.0)
            zoffset = double_or_blank(card, 8, 'zoffset', 0.0)
            blank(card, 9, 'blank')
            tflag = integer_or_blank(card, 10, 'tflag', 0)
            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            T4 = double_or_blank(card, 14, 'T4')
            assert len(card) <= 15, 'len(CQUAD4 card) = %i\ncard=%s' % (len(card), card)
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            tflag = 0
            T1 = None
            T2 = None
            T3 = None
            T4 = None

        return CQUAD4(eid, pid, nids, theta_mcid, zoffset,
                      tflag, T1, T2, T3, T4, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CQUAD4 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        pid = data[1]
        nids = data[2:6]

        theta_mcid = data[6]
        zoffset = data[7]
        tflag = data[8]
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
        assert tflag in [0, 1], 'data=%s tflag=%s' % (data, tflag)
        return CQUAD4(eid, pid, nids, theta_mcid, zoffset,
                      tflag, T1, T2, T3, T4, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CQUAD4 eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.Coord(self.theta_mcid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CQUAD4 eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.safe_coord(self.theta_mcid, self.eid, xref_errors, msg=msg)

    #def x(self, eta, xi, xs):
        #"""Calculate the x-coordinate within the element.

        #Calculates the local xsect x-coordinate provided the desired master
        #coordinates eta and xi.

        #:Args:

        #- `eta (float)`: The eta coordinate in the master coordinate domain.*
        #- `xi (float)`: The xi coordinate in the master coordinate domain.*

        #:Returns:

        #- `x (float)`: The x-coordinate within the element.

        #.. note:: Xi and eta can both vary between -1 and 1 respectively.

        #per AeroComBAT
        #"""
        #return .25*(
            #xs[0]*(1.-xi)*(1.-eta) + xs[1]*(1.+xi)*(1.-eta)+
            #xs[2]*(1.+xi)*(1.+eta) + xs[3]*(1.-xi)*(1.+eta)
        #)

    #def y(self, eta, xi, ys):
        #"""Calculate the y-coordinate within the element.

        #Calculates the local xsect y-coordinate provided the desired master
        #coordinates eta and xi.

        #:Args:

        #- `eta (float)`: The eta coordinate in the master coordinate domain.*
        #- `xi (float)`: The xi coordinate in the master coordinate domain.*

        #:Returns:

        #- `y (float)': The y-coordinate within the element.

        #.. note:: Xi and eta can both vary between -1 and 1 respectively.

        #per AeroComBAT
        #"""
        #return .25*(
            #ys[0]*(1.-xi)*(1.-eta) + ys[1]*(1.+xi)*(1.-eta)+\
            #ys[2]*(1.+xi)*(1.+eta) + ys[3]*(1.-xi)*(1.+eta)
        #)

    #def Z(self, eta, xi, xs, ys):
        #"""Calculates transformation matrix relating stress to force-moments.

        #Intended primarily as a private method but left public, this method
        #calculates the transformation matrix that converts stresses to force
        #and moment resultants.

        #:Args:

        #- `eta (float)`: The eta coordinate in the master coordinate domain.*
        #- `xi (float)`: The xi coordinate in the master coordinate domain.*

        #:Returns:

        #- `Z (3x6 np.array[float])`: The stress-resutlant transformation array.

        #.. note:: Xi and eta can both vary between -1 and 1 respectively.

        #per AeroComBAT
        #"""
        #return np.array([
            #[1., 0, 0, 0, 0, -self.y(eta, xi, ys)],
            #[0, 1., 0, 0, 0, self.x(eta, xi, xs)],
            #[0, 0, 1., self.y(eta, xi, ys), -self.x(eta, xi, xs), 0]
        #])

    #def J(self, eta, xi):
        #"""Calculates the jacobian at a point in the element.

        #This method calculates the jacobian at a local point within the element
        #provided the master coordinates eta and xi.

        #:Args:

        #- `eta (float)`: The eta coordinate in the master coordinate domain.
        #- `xi (float)`: The xi coordinate in the master coordinate domain.

        #:Returns:

        #- `Jmat (3x3 np.array[float])`: The stress-resutlant transformation
            #array.

        #.. note:: Xi and eta can both vary between -1 and 1 respectively.

        #per AeroComBAT
        #"""
        #xs = self.xs
        #ys = self.ys
        #J11 = 0.25*(-xs[0]*(1-eta) + xs[1]*(1-eta) + xs[2]*(1+eta) - xs[3]*(1+eta))
        #J12 = 0.25*(-ys[0]*(1-eta) + ys[1]*(1-eta) + ys[2]*(1+eta) - ys[3]*(1+eta))
        #J21 = 0.25*(-xs[0]*(1-xi) - xs[1]*(1+xi) + xs[2]*(1+xi) + xs[3]*(1-xi))
        #J22 = 0.25*(-ys[0]*(1-xi) - ys[1]*(1+xi) + ys[2]*(1+xi) + ys[3]*(1-xi))
        #Jmat = np.array([
            #[J11, J12, 0.],
            #[J21, J22, 0.],
            #[0., 0., 1.]])
        #return Jmat

    #def _gauss(self):
        #"""
        #per AeroComBAT
        #"""
        #xyz = self.get_node_positions()
        #xs = xyz[:, 0]
        #ys = xyz[:, 0]

        ## Initialize coordinates for Guass Quadrature Integration
        #etas = np.array([-1,1]) * np.sqrt(3)/3
        #xis = np.array([-1,1]) * np.sqrt(3)/3

        ## Evaluate/sum the cross-section matricies at the Guass points
        #for k in range(0, np.size(xis)):
            #for l in range(0, np.size(etas)):
                ##Get Z Matrix
                #Zmat = self.Z(etas[l], xis[k], xs, ys)

                ##Get BN Matricies
                #Jmat = self.J(etas[l], xis[k], xs, ys)

                ##Get determinant of the Jacobian Matrix
                #Jdet = abs(np.linalg.det(Jmat))
                #Jmatinv = np.linalg.inv(Jmat)
                #Bxi = np.zeros((6,3))
                #Beta = np.zeros((6,3))
                #Bxi[0,0] = Bxi[2,1] = Bxi[3,2] = Jmatinv[0,0]
                #Bxi[1,1] = Bxi[2,0] = Bxi[4,2] = Jmatinv[1,0]
                #Beta[0,0] = Beta[2,1] = Beta[3,2] = Jmatinv[0,1]
                #Beta[1,1] = Beta[2,0] = Beta[4,2] = Jmatinv[1,1]
                #BN = np.dot(Bxi,self.dNdxi(etas[l])) + np.dot(Beta,self.dNdeta(xis[k]))

                ##Get a few last minute matricies
                #S = np.zeros((6,3))
                #S[3,0] = 1.
                #S[4,1] = 1.
                #S[5,2] = 1.
                #SZ = S @ Zmat
                #Nmat = self.N(etas[l], xis[k])
                #SN = S @ Nmat

                ## Calculate the mass per unit length of the element
                #self.mass += self.rho * Jdet

                ##Add to Ae Matrix
                #self.Ae += SZ.T @ self.Q @ SZ * Jdet

                ##Add to Re Matrix
                #self.Re += BN.T @ self.Q @ SZ * Jdet

                ##Add to Ee Matrix
                #self.Ee += BN.T @ self.Q @ BN * Jdet

                ##Add to Ce Matrix
                #self.Ce += BN.T @ self.Q @ SN * Jdet

                ##Add to Le Matrix
                #self.Le += SN.T @ self.Q @ SZ * Jdet

                ##Add to Me Matrix
                #self.Me += SN.T @ self.Q @ SN * Jdet

    #@staticmethod
    #def N(eta, xi):
        #"""Generates the shape-function value weighting matrix.

        #Intended primarily as a private method but left public, this method
        #generates the weighting matrix used to interpolate values within the
        #element. This method however is mainly reserved for the cross-sectional
        #analysis process.

        #:Args:

        #- `eta (float)`: The eta coordinate in the master coordinate domain.*
        #- `xi (float)`: The xi coordinate in the master coordinate domain.*

        #:Returns:

        #- `Nmat (3x12 np.array[float])`: The shape-function value weighting
            #matrix.

        #Notes
        #-----
        #Xi and eta can both vary between -1 and 1 respectively.

        #per AeroComBAT
        #"""
        #Nmat = np.zeros([3,12])
        #N1 = .25*(1.-xi)*(1.-eta)
        #N2 = .25*(1.+xi)*(1.-eta)
        #N3 = .25*(1.+xi)*(1.+eta)
        #N4 = .25*(1.-xi)*(1.+eta)
        #Nmat[0,0] = Nmat[1,1] = Nmat[2,2] = N1
        #Nmat[0,3] = Nmat[1,4] = Nmat[2,5] = N2
        #Nmat[0,6] = Nmat[1,7] = Nmat[2,8] = N3
        #Nmat[0,9] = Nmat[1,10] = Nmat[2,11] = N4
        #return Nmat

    #@staticmethod
    #def dNdxi(eta):
        #"""Generates a gradient of the shape-function value weighting matrix.

        #Intended primarily as a private method but left public, this method
        #generates the gradient of the weighting matrix with respect to xi and
        #is used to interpolate values within the element. This method however
        #is mainly reserved for the cross-sectional analysis process.

        #:Args:

        #- `eta (float)`: The eta coordinate in the master coordinate domain.*
        #- `xi (float)`: The xi coordinate in the master coordinate domain.*

        #:Returns:

        #- `dNdxi_mat (3x12 np.array[float])`: The gradient of the shape-
            #function value weighting matrix with respect to xi.

        #.. note:: Xi and eta can both vary between -1 and 1 respectively.

        #per AeroComBAT
        #"""
        #dNdxi_mat = np.zeros([3,12])
        #dN1dxi = -.25*(1-eta)
        #dN2dxi = .25*(1-eta)
        #dN3dxi = .25*(1+eta)
        #dN4dxi = -.25*(1+eta)
        #dNdxi_mat[0,0] = dNdxi_mat[1,1] = dNdxi_mat[2,2] = dN1dxi
        #dNdxi_mat[0,3] = dNdxi_mat[1,4] = dNdxi_mat[2,5] = dN2dxi
        #dNdxi_mat[0,6] = dNdxi_mat[1,7] = dNdxi_mat[2,8] = dN3dxi
        #dNdxi_mat[0,9] = dNdxi_mat[1,10] = dNdxi_mat[2,11] = dN4dxi
        #return dNdxi_mat

    #@staticmethod
    #def dNdeta(xi):
        #"""Generates a gradient of the shape-function value weighting matrix.

        #Intended primarily as a private method but left public, this method
        #generates the gradient of the weighting matrix with respect to eta and
        #is used to interpolate values within the element. This method however
        #is mainly reserved for the cross-sectional analysis process.

        #:Args:

        #- `eta (float)`: The eta coordinate in the master coordinate domain.*
        #- `xi (float)`: The xi coordinate in the master coordinate domain.*

        #:Returns:

        #- `dNdeta_mat (3x12 np.array[float])`: The gradient of the shape-
            #function value weighting matrix with respect to eta.

        #.. note:: Xi and eta can both vary between -1 and 1 respectively.

        #per AeroComBAT
        #"""
        #dNdeta_mat = np.zeros([3,12])
        #dN1deta = -.25*(1-xi)
        #dN2deta = -.25*(1+xi)
        #dN3deta = .25*(1+xi)
        #dN4deta = .25*(1-xi)
        #dNdeta_mat[0,0] = dNdeta_mat[1,1] = dNdeta_mat[2,2] = dN1deta
        #dNdeta_mat[0,3] = dNdeta_mat[1,4] = dNdeta_mat[2,5] = dN2deta
        #dNdeta_mat[0,6] = dNdeta_mat[1,7] = dNdeta_mat[2,8] = dN3deta
        #dNdeta_mat[0,9] = dNdeta_mat[1,10] = dNdeta_mat[2,11] = dN4deta
        #return dNdeta_mat

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.theta_mcid = self.Theta_mcid()
        self.pid_ref = None
        self.nodes_ref = None
        self.theta_mcid_ref = None

    def split_to_ctria3(self, eida, eidb):
        # (int, int) -> (CTRIA3, CTRIA3)
        """
        Splits a CQUAD4 into two CTRIA3s

        .. todo:: doesn't consider theta_mcid if a float correctly (use an integer)
        .. todo:: doesn't optimize the orientation of the nodes yet...

        """
        n1, n2, n3, n4 = self.nodes
        nids = [n1, n2, n3]
        elementa = CTRIA3(eida, self.pid, nids, zoffset=self.zoffset, theta_mcid=self.theta_mcid,
                          tflag=self.tflag, T1=self.T1, T2=self.T2, T3=self.T3,
                          comment=self.comment)
        nids = [n3, n4, n1]
        elementb = CTRIA3(eidb, self.pid, nids, zoffset=self.zoffset, theta_mcid=self.theta_mcid,
                          tflag=self.tflag, T1=self.T1, T2=self.T2, T3=self.T3,
                          comment='')
        return elementa, elementb

    def get_thickness_scale(self):
        return [self.T1, self.T2, self.T3, self.T4]

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        unused_edges = self.get_edge_ids()
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

    def flip_normal(self):
        r"""
        ::

          1---2       1---4
          |   |  -->  |   |
          |   |       |   |
          4---3       2---3

        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]
        if self.nodes_ref is not None:
            (n1, n2, n3, n4) = self.nodes_ref
            self.nodes_ref = [n1, n4, n3, n2]

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)

    def write_as_ctria3(self, new_eid):
        """
        triangle - 012
        triangle - 023

        """
        zoffset = set_blank_if_default(self.zoffset, 0.0)
        nodes1 = [self.nodes[0], self.nodes[1], self.nodes[2]]
        nodes2 = [self.nodes[0], self.nodes[2], self.nodes[3]]
        fields1 = ['CTRIA3', self.eid, self.Pid()] + nodes1 + [
            self.theta_mcid, zoffset]
        fields2 = ['CTRIA3', new_eid, self.Pid()] + nodes2 + [
            self.theta_mcid, zoffset]
        return self.print_card(fields1) + self.print_card(fields2)

    def raw_fields(self):
        list_fields = (['CQUAD4', self.eid, self.Pid()] + self.node_ids +
                       [self.theta_mcid, self.zoffset, None, self.tflag, self.T1, self.T2,
                        self.T3, self.T4])
        return list_fields

    def repr_fields(self):
        (theta_mcid, zoffset, tflag, T1, T2, T3, T4) = self._get_repr_defaults()
        list_fields = (['CQUAD4', self.eid, self.Pid()] + self.node_ids +
                       [theta_mcid, zoffset, None, tflag, T1, T2, T3, T4])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids

        row2_data = [self.theta_mcid, self.zoffset,
                     self.tflag, self.T1, self.T2, self.T3, self.T4]
        if row2_data == [0.0, 0.0, 0, 1.0, 1.0, 1.0, 1.0]:
            data = [self.eid, self.Pid()] + nodes
            msg = ('CQUAD4  %8i%8i%8i%8i%8i%8i\n' % tuple(data))
            return self.comment + msg
        else:
            theta_mcid = self._get_theta_mcid_repr()
            zoffset = set_blank_if_default(self.zoffset, 0.0)
            tflag = set_blank_if_default(self.tflag, 0)
            T1 = set_blank_if_default(self.T1, 1.0)
            T2 = set_blank_if_default(self.T2, 1.0)
            T3 = set_blank_if_default(self.T3, 1.0)
            T4 = set_blank_if_default(self.T4, 1.0)

            row2_data = [theta_mcid, zoffset,
                         tflag, T1, T2, T3, T4]
            if size == 8:
                row2 = [print_field_8(field) for field in row2_data]
                data = [self.eid, self.Pid()] + nodes + row2
                msg = ('CQUAD4  %8i%8i%8i%8i%8i%8i%8s%8s\n'
                       '                %8s%8s%8s%8s%8s' % tuple(data))
                return self.comment + msg.rstrip('\n ') + '\n'
            else:
                row2 = [print_field_16(field) for field in row2_data]
                is_stripped = [field.strip() == '' for field in row2]
                if all(is_stripped[2:]): # tflag, t1234 are blank
                    data = [self.eid, self.Pid()] + nodes + row2[:2]
                    msg = ('CQUAD4* %16i%16i%16i%16i\n'
                           '*       %16i%16i%16s%16s\n'
                           % tuple(data))
                else:
                    data = [self.eid, self.Pid()] + nodes + row2
                    msg = ('CQUAD4* %16i%16i%16i%16i\n'
                           '*       %16i%16i%16s%16s\n'
                           '*                     %16s%16s%16s\n'
                           '*       %16s%16s\n'
                           % tuple(data))
                    return self.comment + msg.rstrip('*\n ') + '\n'
            return self.comment + msg


class CPLSTx4(QuadShell):
    """
    +---------+-------+-------+----+----+----+----+-------+
    |    1    |   2   |   3   |  4 |  5 |  6 | 7  |   8   |
    +=========+=======+=======+====+====+====+====+=======+
    | CPLSTN4 |  EID  |  PID  | N1 | N2 | N3 | N4 | THETA |
    +---------+-------+-------+----+----+----+----+-------+

    """
    _field_map = {1: 'eid', 2:'pid', 7:'theta'}

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

    def __init__(self, eid, pid, nids, theta=0.0, comment=''):
        QuadShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        assert len(nids) == 4, nids
        self.nodes = self.prepare_node_ids(nids)
        self.theta = theta

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        thetas = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
            thetas.append(element.theta)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('theta', data=thetas)

    def validate(self):
        assert len(set(self.nodes)) == 4, 'nodes=%s\n%s' % (self.nodes, str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CPLSTx4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4'),]

        theta = double_or_blank(card, 7, 'theta', 0.0)
        assert len(card) <= 8, 'len(%s card) = %i\ncard=%s' % (self.type, len(card), card)
        return cls(eid, pid, nids, theta, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid}'
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.Pid(), msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid}'
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        unused_edges = self.get_edge_ids()
        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PPLANE', 'PLPLANE', 'PGPLSN'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            #if not self.pid_ref.type in ['PLPLANE']:
                #t = self.Thickness()
                #assert isinstance(t, float), 'thickness=%r' % t
                #mass = self.Mass()
                #assert isinstance(mass, float), 'mass=%r' % mass
            #a, c, n = self.AreaCentroidNormal()
            #assert isinstance(a, float), 'Area=%r' % a
            #for i in range(3):
                #assert isinstance(c[i], float)
                #assert isinstance(n[i], float)

    def flip_normal(self):
        r"""
        ::

          1---2       1---4
          |   |  -->  |   |
          |   |       |   |
          4---3       2---3

        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]
        if self.nodes_ref is not None:
            (n1, n2, n3, n4) = self.nodes_ref
            self.nodes_ref = [n1, n4, n3, n2]

    def Mass(self) -> float:
        return 0.

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)

    def raw_fields(self):
        list_fields = ([self.type, self.eid, self.Pid()] + self.node_ids +
                       [self.theta])
        return list_fields

    def repr_fields(self):
        theta = set_blank_if_default(self.theta, 0.0)
        list_fields = ([self.type, self.eid, self.Pid()] + self.node_ids +
                       [theta])
        return list_fields


class CPLSTS4(CPLSTx4):
    """
    +---------+-------+-------+----+-------+----+-------+-------+------+
    |    1    |   2   |   3   |  4 |   5   |  6 |   7   |   8   |   9  |
    +=========+=======+=======+====+=======+====+=======+=======+======+
    | CPLSTS4 |  EID  |  PID  | N1 |   N2  | N3 |   N4  | THETA |      |
    +---------+-------+-------+----+-------+----+-------+-------+------+
    |         |       |       |    | TFLAG | T1 |   T2  |   T3  |  T4  |
    +---------+-------+-------+----+-------+----+-------+-------+------+

    """
    type = 'CPLSTS4'

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes + [print_float_8(self.theta)]
        msg = ('CPLSTS4 %8i%8i%8i%8i%8i%8i%8s\n' % tuple(data))
        return self.comment + msg


class CPLSTN4(CPLSTx4):
    type = 'CPLSTN4'

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes + [print_float_8(self.theta)]
        msg = ('CPLSTN4 %8i%8i%8i%8i%8i%8i%8s\n' % tuple(data))
        return self.comment + msg


class CPLSTx6(TriShell):

    def __init__(self, eid, pid, nids, theta=0., comment=''):
        TriShell.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.theta = theta
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(nids) == 6, f'error on {self.type}'

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        thetas = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
            thetas.append(element.theta)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('theta', data=thetas)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CPLSTx6 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
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
            integer_or_blank(card, 8, 'n6', 0),
        ]
        if len(card) > 9:
            theta = double_or_blank(card, 9, 'theta', 0.0)
            assert len(card) <= 15, 'len(%s card) = %i\ncard=%s' % (self.type, len(card), card)
        else:
            theta = 0.0
        return cls(eid, pid, nids, theta=theta, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid}'
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.Pid(), msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid}'
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def Mass(self) -> float:
        return 0.

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CPLSTN6 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        pid = data[1]
        nids = data[2:8]
        theta = data[8]
        assert len(nids) == 6, 'error on CPLSTN6'
        return cls(eid, pid, nids, theta, comment=comment)

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        edges = self.get_edge_ids()

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert isinstance(nid, integer_types) or nid is None, 'nid%i is not an integer/None; nid=%s' %(i, nid)

        #if xref:
            #assert self.pid_ref.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            #if not self.pid_ref.type in ['PLPLANE']:
                #t = self.Thickness()
                #assert isinstance(t, float), 'thickness=%r' % t
                #mass = self.Mass()
                #assert isinstance(mass, float), 'mass=%r' % mass
            #a, c, n = self.AreaCentroidNormal()
            #assert isinstance(a, float), 'Area=%r' % a
            #for i in range(3):
                #assert isinstance(c[i], float)
                #assert isinstance(n[i], float)

    def Thickness(self):
        """Returns the thickness, :math:`t`"""
        return self.pid_ref.Thickness()

    def AreaCentroidNormal(self):
        """
        Returns area, centroid, normal as it's more efficient to do them
        together

        """
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
        return _triangle_area_centroid_normal([n1, n2, n3], self)

    def Area(self):
        r"""
        Get the area, :math:`A`.

        .. math:: A = \frac{1}{2} (n_0-n_1) \times (n_0-n_2)"""
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
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
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
        return _normal(n1 - n2, n1 - n3)

    def Centroid(self):
        r"""
        Get the centroid.

        .. math::
          CG = \frac{1}{3} (n_1+n_2+n_3)

        """
        n1, n2, n3 = self.get_node_positions(nodes=self.nodes_ref[:3])
        centroid = (n1 + n2 + n3) / 3.
        return centroid

    def flip_normal(self):
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
        if self.nodes_ref is not None:
            (n1, n2, n3, n4, n5, n6) = self.nodes_ref
            self.nodes_ref = [n1, n3, n2, n6, n5, n4]

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ([self.type, self.eid, self.Pid()] + self.node_ids +
                       [self.theta])
        return list_fields

    def repr_fields(self):
        theta = set_blank_if_default(self.theta, 0.0)
        list_fields = ([self.type, self.eid, self.Pid()] + self.node_ids +
                       [theta])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            msg = self.comment + print_card_8(card)
        msg = self.comment + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg


class CPLSTS6(CPLSTx6):
    type = 'CPLSTS6'

class CPLSTN6(CPLSTx6):
    type = 'CPLSTN6'


class CPLSTx8(QuadShell):
    def __init__(self, eid, pid, nids, theta=0., comment=''):
        QuadShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.theta = theta
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 8

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        thetas = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
            thetas.append(element.theta)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('theta', data=thetas)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CPLSTx8 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5', 0),
                integer_or_blank(card, 8, 'n6', 0),
                integer_or_blank(card, 9, 'n7', 0),
                integer_or_blank(card, 10, 'n8', 0),]
        if len(card) > 11:
            theta = double_or_blank(card, 15, 'theta', 0.0)
            assert len(card) <= 18, 'len(%s card) = %i\ncard=%s' % (self.type, len(card), card)
        else:
            theta = 0.0
        return cls(eid, pid, nids, theta=theta, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #"""
        #Adds a CPLSTx8 card from the OP2

        #Parameters
        #----------
        #data : List[varies]
            #a list of fields defined in OP2 format
        #comment : str; default=''
            #a comment for the card
        #"""
        ##print("CQUAD8 = ",data)
        ##(6401,
        ##6400,
        ##6401, 6402, 6405, 6403, 0, 0, 6404, 0,
        ##-1.0, -1.0, -1.0, -1.0,
        ##0.0, 0)
        #eid = data[0]
        #pid = data[1]
        #nids = data[2:10]
        #theta = data[10]
        #return CPLSTN8(eid, pid, nids, T1, T2, T3, T4, theta, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid}'
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.Pid(), msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid}'
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        unused_edges = self.get_edge_ids()

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert isinstance(nid, integer_types) or nid is None, 'nid%i is not an integer/None; nid=%s' %(i, nid)

        #if xref:
            #assert self.pid_ref.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            #t = self.Thickness()
            #a, c, n = self.AreaCentroidNormal()
            #assert isinstance(t, float), 'thickness=%r' % t
            #assert isinstance(a, float), 'Area=%r' % a
            #for i in range(3):
                #assert isinstance(c[i], float)
                ##assert isinstance(n[i], float)
            #mass = self.Mass()
            #assert isinstance(mass, float), 'mass=%r' % mass

    def Thickness(self):
        """Returns the thickness"""
        return self.pid_ref.Thickness()

    def flip_normal(self):
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
        n1, n2, n3, n4 = self.get_node_positions(nodes=self.nodes_ref[:4])
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
        n1, n2, n3, n4 = self.get_node_positions(nodes=self.nodes[:4])
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
        n1, n2, n3, n4 = self.get_node_positions(nodes=self.nodes[:4])
        a = n1 - n3
        b = n2 - n4
        area = 0.5 * norm(cross(a, b))
        return area

    def Mass(self) -> float:
        return 0.

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = [self.type, self.eid, self.Pid()] + self.node_ids + [
            self.theta]
        return list_fields

    def repr_fields(self):
        theta = set_blank_if_default(self.theta, 0.0)
        list_fields = ([self.type, self.eid, self.Pid()] + self.node_ids + [
            theta])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8: # to last node
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class CPLSTS8(CPLSTx8):
    type = 'CPLSTS8'

class CPLSTN8(CPLSTx8):
    type = 'CPLSTN8'



class CQUADR(QuadShell):
    """
    +--------+-------+-------+----+----+----+----+------------+---------+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |     8      |    9    |
    +========+=======+=======+=====+===+====+====+============+=========+
    | CQUADR |  EID  |  PID  | N1 | N2 | N3 | N4 | THETA/MCID | ZOFFSET |
    +--------+-------+-------+----+----+----+----+------------+---------+
    |        |       | TFLAG | T1 | T2 | T3 | T4 |            |         |
    +--------+-------+-------+----+----+----+----+------------+---------+

    """
    type = 'CQUADR'

    def __init__(self, eid, pid, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                 T1=None, T2=None, T3=None, T4=None, comment=''):
        """
        Creates a CQUADR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 / T4 : float; default=None
            If it is not supplied, then T1 through T4 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        QuadShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.theta_mcid = theta_mcid
        self.zoffset = zoffset
        self.tflag = tflag
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.T4 = T4
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 4, 'CQUADR'
        self.theta_mcid_ref = None  # type: Optional[Any]

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        mcids = []
        thetas = []
        zoffsets = []
        #t1234 = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
            if isinstance(element.theta_mcid, int):
                mcid = element.theta_mcid
                theta = 0.
            else:
                assert isinstance(element.theta_mcid, float), type(element.theta_mcid)
                mcid = -1
                theta = element.theta_mcid
            mcids.append(mcid)
            thetas.append(theta)
            zoffsets.append(element.zoffset)
            #t1234.append([element.T1, element.T2, element.T3, element.T4])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('mcid', data=mcids)
        h5_file.create_dataset('theta', data=thetas)
        h5_file.create_dataset('zoffset', data=zoffsets)
        #self.tflag = tflag

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CQUADR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 'n1'),
                integer_or_blank(card, 4, 'n2'),
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4'),]

        theta_mcid = integer_double_or_blank(card, 7, 'theta_mcid', 0.0)
        zoffset = double_or_blank(card, 8, 'zoffset', 0.0)

        tflag = integer_or_blank(card, 10, 'tflag', 0)
        T1 = double_or_blank(card, 11, 'T1')
        T2 = double_or_blank(card, 12, 'T2')
        T3 = double_or_blank(card, 13, 'T3')
        T4 = double_or_blank(card, 14, 'T4')
        assert len(card) <= 15, 'len(CQUADR card) = %i\ncard=%s' % (len(card), card)
        return CQUADR(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, T4=T4, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CQUADR card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        pid = data[1]
        nids = data[2:6]

        theta_mcid = data[6]
        zoffset = data[7]
        tflag = data[8]
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
        assert tflag in [0, 1], data
        return CQUADR(eid, pid, nids, theta_mcid, zoffset,
                      tflag, T1, T2, T3, T4, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CQUADR eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.Coord(self.theta_mcid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CQUADR eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.safe_coord(self.theta_mcid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.theta_mcid = self.Theta_mcid()
        self.pid_ref = None
        self.nodes_ref = None
        self.theta_mcid_ref = None

    def get_thickness_scale(self):
        return [self.T1, self.T2, self.T3, self.T4]

    def Thickness(self):
        """Returns the thickness"""
        return self.pid_ref.Thickness()

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        unused_nids = self.node_ids

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        #for i,nid in enumerate(nids):
            #assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PSHELL', 'PCOMP', 'PCOMPG'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            t = self.Thickness()
            a, c, unused_n = self.AreaCentroidNormal()
            assert isinstance(t, float), 'thickness=%r' % t
            assert isinstance(a, float), 'Area=%r' % a
            for i in range(3):
                assert isinstance(c[i], float)
                #assert isinstance(n[i], float)
            mass = self.Mass()
            assert isinstance(mass, float), 'mass=%r' % mass

    def flip_normal(self):
        r"""
        ::

          1---2       1---4
          |   |  -->  |   |
          |   |       |   |
          4---3       2---3

        """
        (n1, n2, n3, n4) = self.nodes
        self.nodes = [n1, n4, n3, n2]
        if self.nodes_ref is not None:
            (n1, n2, n3, n4) = self.nodes_ref
            self.nodes_ref = [n1, n4, n3, n2]

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    @node_ids.setter
    def node_ids(self, value):
        raise ValueError("You cannot set node IDs like this...modify the node objects")

    def raw_fields(self):
        list_fields = (['CQUADR', self.eid, self.Pid()] + self.node_ids +
                       [self.theta_mcid, self.zoffset, None, self.tflag, self.T1,
                        self.T2, self.T3, self.T4])
        return list_fields

    def repr_fields(self):
        (theta_mcid, zoffset, tflag, T1, T2, T3, T4) = self._get_repr_defaults()
        list_fields = (['CQUADR', self.eid, self.Pid()] + self.node_ids +
                       [theta_mcid, zoffset, None, tflag, T1, T2, T3, T4])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8 or len(card) == 7: # to last node
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        #msg2 = self.write_card(size)
        #assert msg == msg2, '\n%s---\n%s\n%r\n%r' % (msg, msg2, msg, msg2)
        return msg


class CPLSTS3(TriShell):
    """
    +---------+-------+-------+----+----+----+-------+-------+-----+
    |    1    |   2   |   3   |  4 |  5 |  6 |   7   |   8   |  9  |
    +=========+=======+=======+=====+===+====+=======+=======+=====+
    | CPLSTS3 |  EID  |  PID  | N1 | N2 | N3 | THETA |       |     |
    +---------+-------+-------+----+----+----+-------+-------+-----+
    |         |       | TFLAG | T1 | T2 | T3 |       |       |     |
    +---------+-------+-------+----+----+----+-------+-------+-----+

    """
    type = 'CPLSTS3'
    _field_map = {
        1: 'eid', 2:'pid', 6:'theta', 10:'tflag',
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

    def __init__(self, eid, pid, nids,
                 theta=0.0, tflag=0, T1=1.0, T2=1.0, T3=1.0, comment=''):
        TriShell.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        assert len(nids) == 3, nids
        self.nodes = self.prepare_node_ids(nids)
        self.theta = theta
        self.tflag = tflag
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 3

    def validate(self):
        assert len(set(self.nodes)) == 3, 'nodes=%s; n=%s\n%s' % (self.nodes, len(set(self.nodes)), str(self))

    def Mass(self) -> float:
        return 0.

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #eid = data[0]
        #pid = data[1]
        #nids = data[2:5]

        #theta = data[5]
        #tflag = data[7]
        #T1 = data[8]
        #T2 = data[9]
        #T3 = data[10]
        #if T1 == -1.0:
            #T1 = 1.0
        #if T2 == -1.0:
            #T2 = 1.0
        #if T3 == -1.0:
            #T3 = 1.0
        #return CPLSTS3(eid, pid, nids, zoffset, theta,
                       #tflag, T1, T2, T3, comment=comment)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CPLSTS3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
        ]
        if len(card) > 5:
            theta = double_or_blank(card, 6, 'theta', 0.0)
            blank(card, 8, 'blank')
            blank(card, 9, 'blank')

            tflag = integer_or_blank(card, 10, 'tflag', 0)
            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            assert len(card) <= 14, 'len(CPLSTS3 card) = %i\ncard=%s' % (len(card), card)
        else:
            theta = 0.0
            tflag = 0
            T1 = 1.0
            T2 = 1.0
            T3 = 1.0
        return CPLSTS3(eid, pid, nids, theta,
                       tflag, T1, T2, T3, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CPLSTS3 eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CPLSTS3 eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        unused_edges = self.get_edge_ids()

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)
        for i, nid in enumerate(nids):
            assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        #if xref:
            #assert self.pid_ref.type in ['PSHELL', 'PCOMP', 'PCOMPG', 'PLPLANE'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            #if not self.pid_ref.type in ['PLPLANE']:
                #t = self.Thickness()
                #assert isinstance(t, float), 'thickness=%r' % t
                #mass = self.Mass()
                #assert isinstance(mass, float), 'mass=%r' % mass
            #a, c, n = self.AreaCentroidNormal()
            #assert isinstance(a, float), 'Area=%r' % a
            #for i in range(3):
                #assert isinstance(c[i], float)
                #assert isinstance(n[i], float)

    def flip_normal(self):
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
        if self.nodes_ref is not None:
            (n1, n2, n3) = self.nodes_ref
            self.nodes_ref = [n1, n3, n2]

    def _get_repr_defaults(self):
        tflag = set_blank_if_default(self.tflag, 0)
        theta = set_blank_if_default(self.theta, 0.0)

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)
        return (theta, tflag, T1, T2, T3)

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=False)

    def raw_fields(self):
        list_fields = (['CPLSTS3', self.eid, self.Pid()] + self.node_ids +
                       [self.theta, None, None] +
                       [None, self.tflag, self.T1, self.T2, self.T3])
        return list_fields

    def repr_fields(self):
        (theta, tflag, T1, T2, T3) = self._get_repr_defaults()
        list_fields = (['CPLSTS3', self.eid, self.Pid()] + self.node_ids +
                       [theta, None, None] + [None, tflag, T1, T2, T3])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        (theta, tflag, T1, T2, T3) = self._get_repr_defaults()

        T1 = set_blank_if_default(self.T1, 1.0)
        T2 = set_blank_if_default(self.T2, 1.0)
        T3 = set_blank_if_default(self.T3, 1.0)

        nodes = self.node_ids
        row2_data = [theta, '', tflag, T1, T2, T3]
        row2 = [print_field_8(field) for field in row2_data]
        data = [self.eid, self.Pid()] + nodes + row2
        msg = ('CPLSTS3 %8i%8i%8i%8i%8i%8s%8s\n'
               '                %8s%8s%8s%8s\n' % tuple(data))
        return self.comment + msg.rstrip() + '\n'


class CQUAD(QuadShell):
    """
    +-------+-------+-----+----+------------+----+----+----+----+
    |    1  |   2   |  3  |  4 |     5      |  6 |  7 | 8  |  9 |
    +=======+=======+=====+====+============+====+====+====+====+
    | CQUAD |  EID  | PID | G1 |     G2     | G3 | G4 | G5 | G6 |
    +-------+-------+-----+----+------------+----+----+----+----+
    |       |   G7  | G8  | G9 | THETA/MCID |    |    |    |    |
    +-------+-------+-----+----+------------+----+----+----+----+

    theta_mcid is an MSC specific variable

    """
    type = 'CQUAD'
    #tflag = 1

    def __init__(self, eid, pid, nids, theta_mcid=0., comment=''):
        """
        Creates a CQUAD card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int, int/None, int/None,
                    int/None, int/None, int/None]
            node ids
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        comment : str; default=''
            a comment for the card

        """
        QuadShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.theta_mcid = theta_mcid
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 9
        self.theta_mcid_ref = None  # type: Optional[Any]

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        mcids = []
        thetas = []
        neids = len(eids)
        nodes = np.zeros((neids, 9), dtype='int32')
        for i, eid in enumerate(eids):
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes[i, :] = [nid if nid is not None else 0 for nid in element.nodes]
            if isinstance(element.theta_mcid, int):
                mcid = element.theta_mcid
                theta = 0.
            else:
                assert isinstance(element.theta_mcid, float), type(element.theta_mcid)
                mcid = -1
                theta = element.theta_mcid
            mcids.append(mcid)
            thetas.append(theta)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('mcid', data=mcids)
        h5_file.create_dataset('theta', data=thetas)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CQUAD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
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
                integer_or_blank(card, 11, 'n9'),]
        theta_mcid = integer_double_or_blank(card, 12, 'theta_mcid', 0.)
        assert len(card) <= 13, 'len(CQUAD card) = %i\ncard=%s' % (len(card), card)
        return CQUAD(eid, pid, nids, theta_mcid=theta_mcid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CQUAD card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        pid = data[1]
        nids = data[2:11]
        if len(data) == 11:
            theta_mcid = 0. #  msc specific
        else:
            raise RuntimeError(f'theta_mcid is defined; data={data}')
        return CQUAD(eid, pid, nids, theta_mcid=theta_mcid, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CQUAD eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.Coord(self.theta_mcid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CQUAD eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.safe_coord(self.theta_mcid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.theta_mcid = self.Theta_mcid()
        self.pid_ref = None
        self.nodes_ref = None
        self.theta_mcid_ref = None

    def Area(self):
        r"""
        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert
        where a and b are the quad's cross node point vectors"""
        n1, n2, n3, n4 = self.get_node_positions(nodes=self.nodes_ref[:4])
        a = n1 - n3
        b = n2 - n4
        area = 0.5 * norm(cross(a, b))
        return area

    def Mass(self) -> float:
        r""".. math:: m = \frac{m}{A} A  \f]"""
        A = self.Area()
        mpa = self.pid_ref.MassPerArea()
        try:
            return mpa * A
        except TypeError:
            msg = 'mass/area=%s area=%s prop_type=%s' % (mpa, A, self.pid_ref.type)
            raise TypeError(msg)

    def Thickness(self) -> float:
        """Returns the thickness"""
        return self.pid_ref.Thickness()

    def flip_normal(self) -> None:
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
        if self.nodes_ref is not None:
            (n1, n2, n3, n4, n5, n6, n7, n8, n9) = self.nodes_ref
            self.nodes_ref = [n1, n4, n3, n2, n8, n7, n6, n5, n9]

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def _verify(self, xref):
        pass

    def raw_fields(self):
        list_fields = ['CQUAD', self.eid, self.Pid()] + self.node_ids + [self.theta_mcid]
        return list_fields

    def repr_fields(self):
        list_fields = ['CQUAD', self.eid, self.Pid()] + self.node_ids + [self.theta_mcid]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        nodes2 = ['' if node is None else '%8i' % node for node in nodes[4:]]
        theta_mcid = self.theta_mcid

        data = [self.eid, self.Pid()] + nodes[:4] + nodes2 + [theta_mcid]
        msg = ('CQUAD   %8i%8i%8i%8i%8i%8i%8s%8s\n'  # 6 nodes
               '        %8s%8s%8s%8s\n' % tuple(data))
        return self.comment + msg.rstrip() + '\n'


class CQUAD8(QuadShell):
    """
    +--------+-------+-----+----+----+----+----+------------+-------+
    |    1   |   2   |  3  |  4 |  5 |  6 |  7 |      8     |   9   |
    +========+=======+=====+====+====+====+====+============+=======+
    | CQUAD8 |  EID  | PID | G1 | G2 | G3 | G4 |     G5     |  G6   |
    +--------+-------+-----+----+----+----+----+------------+-------+
    |        |   G7  | G8  | T1 | T2 | T3 | T4 | THETA/MCID | ZOFFS |
    +--------+-------+-----+----+----+----+----+------------+-------+
    |        | TFLAG |     |    |    |    |    |            |       |
    +--------+-------+-----+----+----+----+----+------------+-------+

    """
    type = 'CQUAD8'
    def __init__(self, eid, pid, nids, theta_mcid=0., zoffset=0.,
                 tflag=0, T1=None, T2=None, T3=None, T4=None,
                 comment=''):
        """
        Creates a CQUAD8 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : List[int, int, int, int, int/None, int/None, int/None, int/None]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 / T4 : float; default=None
            If it is not supplied, then T1 through T4 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        QuadShell.__init__(self)
        if comment:
            self.comment = comment
        #: Element ID
        self.eid = eid
        #: Property ID
        self.pid = pid
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.T4 = T4
        self.tflag = tflag
        self.theta_mcid = theta_mcid
        self.zoffset = zoffset
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 8
        self.theta_mcid_ref = None  # type: Optional[Any]

    @classmethod
    def export_to_hdf5(cls, h5_file: Any, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        mcids = []
        thetas = []
        zoffsets = []
        #t1234 = []
        neids = len(eids)
        nodes = np.zeros((neids, 8), dtype='int32')
        for i, eid in enumerate(eids):
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes[i, :] = [nid if nid is not None else 0 for nid in element.nodes]
            if isinstance(element.theta_mcid, int):
                mcid = element.theta_mcid
                theta = 0.
            else:
                assert isinstance(element.theta_mcid, float), type(element.theta_mcid)
                mcid = -1
                theta = element.theta_mcid
            mcids.append(mcid)
            thetas.append(theta)
            zoffsets.append(element.zoffset)
            #t1234.append([element.T1, element.T2, element.T3, element.T4])
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('mcid', data=mcids)
        h5_file.create_dataset('theta', data=thetas)
        h5_file.create_dataset('zoffset', data=zoffsets)
        #self.tflag = tflag

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CQUAD8 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5', 0),
                integer_or_blank(card, 8, 'n6', 0),
                integer_or_blank(card, 9, 'n7', 0),
                integer_or_blank(card, 10, 'n8', 0),]
        if len(card) > 11:
            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            T4 = double_or_blank(card, 14, 'T4')
            theta_mcid = integer_double_or_blank(card, 15, 'theta_mcid', 0.0)
            zoffset = double_or_blank(card, 16, 'zoffset', 0.0)
            tflag = integer_or_blank(card, 17, 'tflag', 0)
            assert len(card) <= 18, 'len(CQUAD8 card) = %i\ncard=%s' % (len(card), card)
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            T1 = None
            T2 = None
            T3 = None
            T4 = None
            tflag = 0
        return CQUAD8(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, T4=T4,
                      comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CQUAD8 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        #print("CQUAD8 = ",data)
        #(6401,
        #6400,
        #6401, 6402, 6405, 6403, 0, 0, 6404, 0,
        #-1.0, -1.0, -1.0, -1.0,
        #0.0, 0)
        #(eid, pid, n1, n2, n3, n4, n5, n6, n7, n8,
         #t1, t2, t3, t4, theta, zoffs, tflag) = out
        eid = data[0]
        pid = data[1]
        nids = data[2:10]
        T1 = data[10]
        T2 = data[11]
        T3 = data[12]
        T4 = data[13]
        theta_mcid = data[14]
        zoffset = data[15]
        tflag = data[16]
        #assert isinstance(tflag, int), tflag  # None for v2001; int for post-2001
        if tflag is None: # TODO: temporary
            tflag = 0
        assert tflag in [0, 1], 'tflag=%s data=%s' % (tflag, data)
        return CQUAD8(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                      tflag=tflag, T1=T1, T2=T2, T3=T3, T4=T4, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CQUAD8 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.Pid(), msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.Coord(self.theta_mcid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CQUAD8 eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        if isinstance(self.theta_mcid, integer_types):
            self.theta_mcid_ref = model.safe_coord(self.theta_mcid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.theta_mcid = self.Theta_mcid()
        self.pid_ref = None
        self.nodes_ref = None
        self.theta_mcid_ref = None

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids
        unused_edges = self.get_edge_ids()

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

    def get_thickness_scale(self):
        return [self.T1, self.T2, self.T3, self.T4]

    def Thickness(self):
        """Returns the thickness"""
        return self.pid_ref.Thickness()

    def flip_normal(self):
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
        if self.nodes_ref is not None:
            (n1, n2, n3, n4, n5, n6, n7, n8) = self.nodes_ref
            self.nodes_ref = [n1, n4, n3, n2, n8, n7, n6, n5]

    def Normal(self):
        n1, n2, n3, n4 = self.get_node_positions(nodes=self.nodes_ref[:4])
        return _normal4(n1, n2, n3, n4, self)

    def AreaCentroid(self) -> Tuple[float, np.ndarray]:
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
        n1, n2, n3, n4 = self.get_node_positions(nodes=self.nodes_ref[:4])
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

    def Area(self) -> float:
        r"""
        .. math:: A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert
        where a and b are the quad's cross node point vectors"""
        n1, n2, n3, n4 = self.get_node_positions(nodes=self.nodes_ref[:4])
        a = n1 - n3
        b = n2 - n4
        area = 0.5 * norm(cross(a, b))
        return area

    @property
    def node_ids(self) -> List[Union[int, None]]:
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CQUAD8', self.eid, self.Pid()] + self.node_ids + [
            self.T1, self.T2, self.T3, self.T4, self.theta_mcid, self.zoffset,
            self.tflag]
        return list_fields

    def repr_fields(self):
        (theta_mcid, zoffset, tflag, T1, T2, T3, T4) = self._get_repr_defaults()
        list_fields = (['CQUAD8', self.eid, self.Pid()] + self.node_ids + [
            T1, T2, T3, T4, theta_mcid, zoffset, tflag])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8 or len(card) == 11: # to last node
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class SNORM(BaseCard):
    """
    +--------+-------+-------+----+-----+----+
    |   1    |   2   |   3   |  4 |  5  |  6 |
    +========+=======+=======+====+=====+====+
    |  SNORM |  GID  |  CID  | N1 |  N2 | N3 |
    +--------+-------+-------+----+-----+----+
    |  SNORM |   3   |   2   | 0. | -1. | 0. |
    +--------+-------+-------+----+-----+----+

    """
    type = 'SNORM'

    @classmethod
    def _init_from_empty(cls):
        nid = 1
        normal = [0.1, 0.4, 0.3]
        return SNORM(nid, normal, cid=0, comment='')

    def __init__(self, nid, normal, cid=0, comment=''):
        """
        Creates an SNORM card

        Parameters
        ----------
        nid : int
            node id
        cid : int
            coordinate system
        normal : List[float, float, float]
            normal vector
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.nid = nid
        self.cid = cid

        self.normal = np.asarray(normal)
        self.cid_ref = None  # type: Optional[Any]

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an SNORM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nid = integer(card, 1, 'nid')
        cid = integer_or_blank(card, 2, 'cid', default=0)

        normal = [
            double_or_blank(card, 3, 'n1', default=0.0),
            double_or_blank(card, 4, 'n2', default=0.0),
            double_or_blank(card, 5, 'n3', default=0.0),
        ]

        assert len(card) <= 6, 'len(SNORM card) = %i\ncard=%s' % (len(card), card)
        return SNORM(nid, normal, cid=cid, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SNORM nid=%s' % self.nid
        self.cid_ref = model.Coord(self.cid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SNORM nid=%s' % self.nid
        self.cid_ref = model.safe_coord(self.cid, self.nid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cid = self.Cid()
        self.cid_ref = None

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def raw_fields(self):
        list_fields = ['SNORM', self.nid, self.Cid()] + list(self.normal)
        return list_fields

    def repr_fields(self):
        list_fields = ['SNORM', self.nid, self.Cid()] + list(self.normal)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 5: # to last node
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        return msg

def transform_shell_material_coordinate_system(cids: List[int],
                                               iaxes, theta_mcid, normal, p1, p2,
                                               idtype='int32', fdtype='float64'):
    """calculate the material transformation matrix"""
    is_mcid = np.array([isinstance(val, integer_types) for val in theta_mcid], dtype='bool')
    nelements = len(theta_mcid)
    imcid = np.where(is_mcid)[0]
    itheta = np.where(~is_mcid)[0]
    ielem = np.arange(nelements)
    mcid = np.array(theta_mcid, dtype=idtype)[imcid]
    #thetad = np.array(theta_mcid, dtype='float64')[itheta]
    #print('  nmcid=%s ntheta=%s' % (len(mcid), len(thetad)))
    #print('mcid', mcid)
    #print('thetad', thetad)

    assert nelements > 0, nelements
    telem = np.zeros((nelements, 3, 3), dtype=fdtype)
    if len(imcid):
        icids = np.searchsorted(cids, mcid)
        iaxesi = iaxes[icids, :]
        jmat = np.cross(normal[imcid, :], iaxesi) # k x i
        ji = np.linalg.norm(jmat, axis=1)[:, np.newaxis]
        ipos = np.where(ji > 0.)[0]
        #izero = np.where(ji == 0.)[0]
        #jmat[izero] = np.nan
        jmat[ipos] /= ji[ipos]
        # we do an extra normalization here because
        # we had to project i onto the elemental plane
        # unlike in the next block
        imat = np.cross(jmat, normal[imcid, :])
        telem[imcid, 0, :] = imat
        telem[imcid, 1, :] = jmat

    if len(itheta):
        assert len(p1.shape) == 2, p1.shape
        assert len(p2.shape) == 2, p1.shape
        imat = p2[itheta, :] - p1[itheta, :]
        assert len(imat.shape) == 2, imat.shape
        ni = np.linalg.norm(imat, axis=1)
        imat /= ni[:, np.newaxis]
        jmat = np.cross(normal[itheta, :], imat) # k x i
        nj = np.linalg.norm(jmat, axis=1)
        jmat /= nj[:, np.newaxis]
        telem[itheta, 0, :] = imat
        telem[itheta, 1, :] = jmat
    # rotate imat and jmat into the -theta frame?
    # the next step is to rotate the [S] matrix into the global xyz frame

    #       [ imat ]
    #[Te] = [ jmat ]
    #       [normal]
    # [T] = [ti][ti^T]
    # [Sbar] = [T_inv][S][T_inv^T]
    # [Sg] = [Te^T][Sbar][Te]
    #
    telem[ielem, 2, :] = normal

    T = np.eye(3, dtype=fdtype)
    #K = T
    K = np.zeros((6, 6), dtype=fdtype)
    K[:3, :3] = T
    K2 = K[:3, 3:]
    K3 = K[3:, :3]
    K4 = K[3:, 3:]
    for i in range(3):
        i1 = (i + 1) % 3
        i2 = (i + 2) % 3
        for j in range(3):
            j1 = (j + 1) % 3
            j2 = (j + 2) % 3
            K2[i, j] = T[i, j1] * T[i, j2]
            K3[i, j] = T[i1, j] * T[i2, j]
            K4[i, j] = T[i1, j1] * T[i2, j2] + T[i1, j2] * T[i2, j1]
    #K2[0, 0] = K2[0, 1]
    #K2[0, 1] = K2[0, 2]
    #K2[0, 2] = K2[0, 0]
    return telem

