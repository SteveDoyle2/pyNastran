"""
All shell elements are defined in this file.  This includes:

 * CQUAD1
 * CTRSHL

All tris are TriShell, ShellElement, and Element objects.
All quads are QuadShell, ShellElement, and Element objects.

"""
from __future__ import annotations
from typing import Optional, Any, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import Property, Material
#from pyNastran.bdf.cards.optimization import break_word_by_trailing_integer
from pyNastran.bdf.cards.materials import get_mat_props_S
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank,
    blank,
)

from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.utils import wipe_empty_fields

from pyNastran.bdf.cards.elements.shell import TriShell, QuadShell
from pyNastran.bdf.cards.properties.shell import get_2d_plate_transform

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF  # , MAT1, MAT8, MAT9
    #from pyNastran.nptyping_interface import NDArray3float

class CQUAD1(QuadShell):
    """
    +--------+-------+-------+----+----+----+----+-------+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |    8  |
    +========+=======+=======+=====+===+====+====+=======+
    | CQUAD4 |  EID  |  PID  | N1 | N2 | N3 | N4 | THETA |
    +--------+-------+-------+----+----+----+----+-------+
    |        |       | TFLAG | T1 | T2 | T3 | T4 | 0.0   |
    +--------+-------+-------+----+----+----+----+-------+
    $       eid     pid     n1      n2      n3      n4      theta
    CQUAD1  1       23      1       2       13      12      .00
    $       pid                     ???     ts/t???                 ???
    PQUAD1  23                      8       .6666667                13.55715

    """
    type = 'CQUAD1'
    #cp_name_map = {
    #}
    _field_map = {1: 'eid', 2:'pid', 7:'theta'}
    _properties = ['_field_map']

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
        thetas = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)

            theta = element.theta
            assert isinstance(theta, float), type(theta)
            thetas.append(theta)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('theta', data=thetas)

    def __init__(self, eid, pid, nids, theta=0.0, comment=''):
        """
        Creates a CQUAD1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int]
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
        self.theta = theta
        self.theta_ref = None  # type: Optional[Any]

    def validate(self):
        assert len(set(self.nodes)) == 4, 'nodes=%s\n%s' % (self.nodes, str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CQUAD1 card from ``BDF.add_card(...)``

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
        assert len(card) <= 7, f'len(CQUAD1 card) = {len(card):d}\ncard={card}'
        return CQUAD1(eid, pid, nids, theta, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CQUAD1 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        raise NotImplementedError('CQUAD1')
        eid = data[0]
        pid = data[1]
        nids = data[2:6]

        theta = data[6]
        for nid in nids:
            assert nid > 0, nids
        return CQUAD1(eid, pid, nids, theta, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CQUAD1 eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors):
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

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.pid_ref = None
        self.nodes_ref = None
        self.theta_mcid_ref = None

    #def split_to_ctria3(self, eida, eidb):
        ## (int, int) -> (CTRIA3, CTRIA3)
        #"""
        #Splits a CQUAD4 into two CTRIA3s

        #.. todo:: doesn't consider theta_mcid if a float correctly (use an integer)
        #.. todo:: doesn't optimize the orientation of the nodes yet...

        #"""
        #n1, n2, n3, n4 = self.nodes
        #nids = [n1, n2, n3]
        #elementa = CTRIA3(eida, self.pid, nids, zoffset=self.zoffset, theta_mcid=self.theta_mcid,
                          #tflag=self.tflag, T1=self.T1, T2=self.T2, T3=self.T3,
                          #comment=self.comment)
        #nids = [n3, n4, n1]
        #elementb = CTRIA3(eidb, self.pid, nids, zoffset=self.zoffset, theta_mcid=self.theta_mcid,
                          #tflag=self.tflag, T1=self.T1, T2=self.T2, T3=self.T3,
                          #comment='')
        #return elementa, elementb

    #def get_thickness_scale(self):
        #return [self.T1, self.T2, self.T3, self.T4]

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
            assert self.pid_ref.type in ['PQUAD1'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
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

    #def write_as_ctria3(self, new_eid):
        #"""
        #triangle - 012
        #triangle - 023

        #"""
        #zoffset = set_blank_if_default(self.zoffset, 0.0)
        #nodes1 = [self.nodes[0], self.nodes[1], self.nodes[2]]
        #nodes2 = [self.nodes[0], self.nodes[2], self.nodes[3]]
        #fields1 = ['CTRIA3', self.eid, self.Pid()] + nodes1 + [
            #self.theta_mcid, zoffset]
        #fields2 = ['CTRIA3', new_eid, self.Pid()] + nodes2 + [
            #self.theta_mcid, zoffset]
        #return self.print_card(fields1) + self.print_card(fields2)

    def raw_fields(self) -> list[Any]:
        list_fields = (['CQUAD1', self.eid, self.Pid()] + self.node_ids +
                       [self.theta])
        return list_fields

    def repr_fields(self) -> list[Any]:
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        data = [self.eid, self.Pid()] + nodes + [self.theta]
        if size == 8:
            msg = ('CQUAD1  %8i%8i%8i%8i%8i%8i%8s\n' % tuple(data))
        else:
            msg = ('CQUAD1* %16i%16i%16i%16i\n'
                   '*       %16i%16i%16s\n' % tuple(data))
        return self.comment + msg.rstrip('\n ') + '\n'


class CTRSHL(TriShell):
    """
    +--------+-------+-------+----+----+----+------------+---------+
    |   1    |   2   |   3   |  4 |  5 |  6 |     7      |    8    |
    +========+=======+=======+=====+===+====+============+=========+
    | CTRIA3 |  EID  |  PID  | N1 | N2 | N3 | THETA/MCID | ZOFFSET |
    +--------+-------+-------+----+----+----+------------+---------+

    """
    type = 'CTRSHL'
    def __init__(self, eid: int, pid: int, nids: list[int],
                 theta: float,
                 comment: str='') -> None:
        """
        Creates a CTRIA3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        TriShell.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.theta = theta
        assert len(nids) == 6, nids
        self.nodes = self.prepare_node_ids(nids)
        assert len(self.nodes) == 6
        self.theta_ref = None  # type: Optional[Any]

    def validate(self):
        assert len(set(self.nodes)) == 6, 'nodes=%s; n=%s\n%s' % (self.nodes, len(set(self.nodes)), str(self))

    @classmethod
    def add_card(cls: Any, card: Any, comment: str=''):
        """
        Adds a CTRSHL card from ``BDF.add_card(...)``

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
            integer(card, 6, 'n4'),
            integer(card, 7, 'n5'),
            integer(card, 8, 'n6'),
        ]
        theta = double_or_blank(card, 9, 'theta', default=0.0)
        if len(card) > 6:
            assert len(card) <= 10, f'len(CTRSHL card) = {len(card):d}\ncard={card}'

        return CTRSHL(eid, pid, nids, theta=theta, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CTRSHL eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.Pid(), msg=msg)
        #if isinstance(self.theta_mcid, integer_types):
            #self.theta_mcid_ref = model.Coord(self.theta_mcid, msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CTRSHL eid=%s' % self.eid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        #if isinstance(self.theta_mcid, integer_types):
            #self.theta_mcid_ref = model.safe_coord(self.theta_mcid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.theta_mcid = self.Theta_mcid()
        self.pid_ref = None
        self.nodes_ref = None
        self.theta_mcid_ref = None

    def Theta_mcid(self) -> float:
        return self.theta

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CTRSHL', self.eid, self.Pid()] + self.node_ids + [self.theta]
        return list_fields

    def repr_fields(self):
        list_fields = ['CTRSHL', self.eid, self.Pid()] + self.node_ids + [self.theta]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = wipe_empty_fields(self.repr_fields())
        if size == 8 or len(card) == 8: # to last node
            msg = self.comment + print_card_8(card)
        else:
            msg = self.comment + print_card_16(card)
        return msg



class PQUAD1(Property):
    """
    'PQUAD1  23                      8       .6666667                13.55715'
    |   1    |   2   |   3   |   4   |   5   |   6   |   7   |   8   |   9   | 10  |
    |--------|-------|-------|-------|-------|-------|-------|-------|-------|-----|
    |PQUAD1  |  PID  |  MID1 |   T1  |  MID2 |   I   |  MID3 |   T3  |  NSM  |abc  |
    |+bc     |   Z1  |   Z2  |       |       |       |       |       |       |     |
    |PQUAD1  |   32  |   16  |  2.98 |   9   | 6.45  |  16   |  5.29 |  6.32 |WXYZ1|
    |+XYZ1   |  0.09 | -0.06 |       |       |       |       |       |       |     |
    """
    type = 'PQUAD1'
    _field_map = {
        1: 'pid', 2:'mid1', 3:'t', 4:'mid2', 5:'twelveIt3', 6:'mid3',
        7: 'tst', 8:'nsm',
        11:'z1', 12:'z2',
    }
    pname_fid_map = {
        # 1 based
        4 : 't', 'T' : 't',
        6 : 'twelveIt3', # no option
        8 : 'tst', #'T' : 't',
    }

    def __init__(self, pid, mid1=None, t_membrane=None, mid2=None, inertia=1.0,
                 mid3=None, t_shear=0.833333, nsm=0.0,
                 z1=None, z2=None, comment=''):
        """
        Creates a PSHELL card

        Parameters
        ----------
        pid : int
            property id
        mid1 : int
            Material id number for membrane (Integer > 0).
        t_membrane : float
            Membrane thickness (Real).
        mid2 : float
            Material id number for bending (Integer > 0).
        inertia : float
            Area moment of inertia per unit width.
        mid3 : int
            Material id number for transverse shear (Integer >= 0).
        t3 : float
            Transverse shear thickness.
        nsm : float
            Nonstructural mass per unit area.
        z1, z2 : float; default=None
            Fiber distances for stress computation, positive according to the
            right-hand sequence defined on the CQUAD1 card.

        pid : int
            property id
        mid1 : int; default=None
            defines membrane material
            defines element density (unless blank)
        mid2 : int; default=None
            defines bending material
            defines element density if mid1=None
        mid3 : int; default=None
            defines transverse shear material
            (only defined if mid2 > 0)
        mid4 : int; default=None
            defines membrane-bending coupling material
            (only defined if mid1 > 0 and mid2 > 0; can't be mid1/mid2)
        twelveIt3 : float; default=1.0
            Bending moment of inertia ratio, 12I/T^3. Ratio of the actual
            bending moment inertia of the shell, I, to the bending
            moment of inertia of a homogeneous shell, T^3/12. The default
            value is for a homogeneous shell.
        nsm : float; default=0.0
            non-structural mass per unit area
        z1 / z2 : float; default=None
            fiber distance location 1/2 for stress/strain calculations
            z1 default : -t/2 if thickness is defined
            z2 default : t/2 if thickness is defined
        comment : str; default=''
            a comment for the card

        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        #if mid2 == -1:
            #mid2 = None

        #: Property ID
        self.pid = pid
        self.mid1 = mid1
        #: Material identification number for bending
        #: -1 for plane strin
        self.mid2 = mid2
        self.mid3 = mid3

        #: thickness
        self.t_membrane = t_membrane
        self.t_shear = t_shear

        self.inertia = inertia

        #: Non-structural Mass
        self.nsm = nsm

        #if z1 is None and self.t is not None:
            #z1 = -self.t / 2.
        #if z2 is None and self.t is not None:
            #z2 = self.t / 2.

        self.z1 = z1
        self.z2 = z2

        #if self.t is not None:
            #assert self.t >= 0.0, 'PSHELL pid=%s Thickness=%s must be >= 0' % (self.pid, self.t)

        self.mid1_ref = None
        self.mid2_ref = None
        self.mid3_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, pids):
        """exports the properties in a vectorized way"""
        #comments = []
        mids = []
        npids = len(pids)
        assert npids > 0, pids

        z = np.full((npids, 2), np.nan, dtype=None, order='C')
        t_membrane = np.full(npids, np.nan, dtype=None, order='C')
        inertia = np.full(npids, np.nan, dtype=None, order='C')
        t_shear = np.full(npids, np.nan, dtype=None, order='C')
        nsm = []
        for i, pid in enumerate(pids):
            prop = model.properties[pid]
            #comments.append(prop.comment)
            midsi = [0 if mid is None else mid for mid in
                     [prop.mid1, prop.mid2, prop.mid3]]
            mids.append(list(midsi))
            z[i, :] = [prop.z1, prop.z2]
            t_membrane[i] = prop.t_membrane
            inertia[i] = prop.inertia
            t_shear[i] = prop.t_shear
            nsm.append(prop.nsm)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('mids', data=mids)
        #print('z =', z)
        #print('t =', t)
        h5_file.create_dataset('z', data=z)
        h5_file.create_dataset('t_membrane', data=t_membrane)
        h5_file.create_dataset('inertia', data=inertia)
        h5_file.create_dataset('t_shear', data=t_shear)
        h5_file.create_dataset('nsm', data=nsm)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PQUAD1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #|   1    |   2   |   3   |   4   |   5   |   6   |   7   |   8   |   9   | 10  |
        #|--------|-------|-------|-------|-------|-------|-------|-------|-------|-----|
        #|PQUAD1  |  PID  |  MID1 |   T1  |  MID2 |   I   |  MID3 |   T3  |  NSM  |abc  |
        #|+bc     |   Z1  |   Z2  |       |       |       |       |       |       |     |
        #|PQUAD1  |   32  |   16  |  2.98 |   9   | 6.45  |  16   |  5.29 |  6.32 |WXYZ1|
        #|+XYZ1   |  0.09 | -0.06 |       |       |       |       |       |       |     |

        pid = integer(card, 1, 'pid')
        mid1 = blank(card, 2, 'mid1')
        t_membrane = blank(card, 3, 't_membrane')

        mid2 = integer(card, 4, 'mid2')
        inertia = double(card, 5, 'inertia')  # poor name
        mid3 = blank(card, 6, 'mid3')
        t_shear = blank(card, 7, 't_shear')
        nsm = double(card, 8, 'nsm')

        if t_membrane is not None:
            t_over_2 = t_membrane / 2.
            z1 = double_or_blank(card, 9, 'z1', -t_over_2)
            z2 = double_or_blank(card, 10, 'z2', t_over_2)
        else:
            z1 = double_or_blank(card, 9, 'z1')
            z2 = double_or_blank(card, 10, 'z2')

        #if self.mid2 is None:
        #    assert self.mid3 is None
        #else: # mid2 is defined
        #    #print (self.mid2 = ", self.mid2)
        #    assert self.mid2 >= -1
        #    #assert self.mid3 >  0
        assert len(card) <= 11, f'len(PQUAD1 card) = {len(card):d}\ncard={card}'
        return PQUAD1(pid, mid1, t_membrane, mid2, inertia=inertia,
                      mid3=mid3, t_shear=t_shear, nsm=nsm,
                      z1=z1, z2=z2, comment=comment)

    def _verify(self, xref):
        pid = self.Pid()
        mid = self.Mid()
        mid1 = self.Mid1()
        mid2 = self.Mid2()
        mid3 = self.Mid3()

        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(mid, integer_types), 'mid=%r' % mid
        assert mid1 is None or isinstance(mid1, integer_types), 'mid1=%r' % mid1
        assert mid2 is None or isinstance(mid2, integer_types), 'mid2=%r' % mid2
        assert mid3 is None or isinstance(mid3, integer_types), 'mid3=%r' % mid3

        mids = [mid for mid in [self.mid1, self.mid2, self.mid3]
                if mid is not None]
        unused_material_ids = self.material_ids
        assert len(mids) > 0
        if xref:
            assert isinstance(self.mid_ref, Material), 'mid=%r' % self.mid_ref

            mids_ref = [self.mid1_ref, self.mid2_ref, self.mid3_ref]
            for i, mid_ref in enumerate(mids_ref):
                if mid_ref is None or mid_ref == 0:
                    continue
                if i == 1: # mid2
                    if isinstance(mid_ref, integer_types):
                        assert mid_ref == -1, mid_ref
                        continue
                assert isinstance(mid_ref, Material), 'mid_ref=%r' % mid_ref
                if mid_ref.type == 'MAT1':
                    E = mid_ref.E()
                    G = mid_ref.G()
                    nu = mid_ref.Nu()
                    rho = mid_ref.Rho()
                    assert isinstance(E, float), 'E=%r' % E
                    assert isinstance(G, float), 'G=%r' % G
                    assert isinstance(nu, float), 'nu=%r' % nu
                    assert isinstance(rho, float), 'rho=%r' % rho
                elif mid_ref.type in ['MAT2', 'MAT4', 'MAT5', 'MAT8']:
                    pass
                #elif mid_ref.type == 'MAT2':
                    #pass
                #elif mid_ref.type == 'MAT4':
                    #pass
                #elif mid_ref.type == 'MAT5':
                    #pass
                #elif mid_ref.type == 'MAT8':
                    #pass
                else:
                    raise NotImplementedError('PQUAD1: pid=%s mid_ref.type=%s' % (
                        self.pid, mid_ref.type))

            t = self.Thickness()
            nsm = self.Nsm()
            mpa = self.MassPerArea()
            assert isinstance(t, float), 't=%r' % t
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            assert isinstance(mpa, float), 'mass_per_area=%r' % mpa

    def get_z_locations(self):
        """returns the locations of the bottom and top surface of the shell"""
        z = np.array([self.z1, self.z2])
        return z

    def materials(self):
        """returns the material objects referenced by the shell"""
        materials = [self.mid1_ref, self.mid2_ref, self.mid3_ref]
        return materials

    @property
    def material_ids(self):
        """returns the material ids"""
        return [self.Mid1(), self.Mid2(), self.Mid3()]

    #@property
    #def mid(self):
        #raise RuntimeError('use self.mid1, self.mid2, self.mid3,')

    #@mid.setter
    #def mid(self, value):
        #raise RuntimeError('use self.mid1, self.mid2, self.mid3, or self.mid4')

    @property
    def mid_ref(self):
        """returns the material used for mass"""
        if self.mid1_ref is not None:
            return self.mid1_ref
        return self.mid2_ref

    def Mid(self):
        """returns the material id used for mass"""
        mid1 = self.Mid1()
        if mid1 is not None:
            return mid1
        return self.Mid2()

    def Mid1(self):
        """returns the extension material id"""
        if self.mid1_ref is not None:
            return self.mid1_ref.mid
        return self.mid1

    def Mid2(self):
        """returns the bending material id"""
        if self.mid2_ref is not None:
            return self.mid2_ref.mid
        return self.mid2

    def Mid3(self):
        if self.mid3_ref is not None:
            return self.mid3_ref.mid
        return self.mid3

    def Thickness(self, tflag=1, tscales=None):
        """returns the thickness of the element"""
        return self.t_membrane

    def Rho(self):
        """returns the material density"""
        return self.mid_ref.rho

    def Nsm(self):
        """returns the non-structural mass"""
        return self.nsm

    def MassPerArea(self, tflag=1, tscales=None):
        """
        Calculates mass per area.

        .. math:: \frac{m}{A} = nsm + \rho t"""
        mid_ref = self.mid_ref
        rho = mid_ref.Rho()  # fails if mid1=None and mid2=None

        thickness = self.Thickness()
        try:
            mass_per_area = self.nsm + rho * thickness
        except Exception:
            print(f'nsm={self.nsm} rho={rho} t_membrane={self.t_membrane}')
            raise
        return mass_per_area

    def MassPerArea_no_xref(self, model, tflag=1, tscales=None):
        """
        Calculates mass per area.

        .. math:: \frac{m}{A} = nsm + \rho t"""
        mid_ref = model.Material(self.Mid())
        rho = mid_ref.Rho()
        thickness = self.Thickness(tflag=tflag, tscales=tscales)
        try:
            mass_per_area = self.nsm + rho * thickness
        except Exception:
            print(f'nsm={self.nsm} rho={rho} t_membrane={self.t_membrane}')
            raise
        return mass_per_area

    def MassPerArea_structure(self):
        """
        Calculates mass per area without considering non-structural mass.

        .. math:: \frac{m}{A} = nsm + \rho t"""
        mid_ref = self.mid_ref
        rho = mid_ref.Rho()
        try:
            mass_per_area = rho * self.t
        except Exception:
            print(f'nsm={self.nsm} rho={rho} t_membrane={self.t_membrane}')
            raise
        return mass_per_area

    def get_Qbar_matrix(self, mid_ref, theta=0.):
        """theta must be in radians"""
        S2, unused_S3 = get_mat_props_S(mid_ref)
        T = get_2d_plate_transform(theta)
        #Tinv = np.linalg.inv(T)
        #Qbar = np.linalg.multi_dot([Tinv, Q, Tinv.T])
        Sbar = np.linalg.multi_dot([T.T, S2, T])
        Qbar = np.linalg.inv(Sbar)
        return Qbar

    def get_Sbar_matrix(self, mid_ref, theta=0.):
        """theta must be in radians"""
        # this is the inverse of Sbar
        S2, unused_S3 = get_mat_props_S(mid_ref)
        T = get_2d_plate_transform(theta)
        #Tinv = np.linalg.inv(T)
        Sbar = np.linalg.multi_dot([T.T, S2, T])
        return Sbar

    def get_individual_ABD_matrices(
            self, theta_offset: float=0.) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Gets the ABD matrix

        Parameters
        ----------
        theta_offset : float
            rotates the ABD matrix; measured in degrees

        http://www2.me.rochester.edu/courses/ME204/nx_help/index.html#uid:id503291
        Understanding Classical Lamination Theory
        """
        #mids = self.get_material_ids()
        thickness = self.Thickness()

        #z0 = self.z1
        #z1 = self.z2
        #zmean = (z0 + z1) / 2.
        #dz = z1 - z0
        #dzsquared = z1 ** 2 - z0 ** 2
        #zcubed = z1 ** 3 - z0 ** 3
        # A11 A12 A16
        # A12 A22 A26
        # A16 A26 A66
        A = np.zeros((3, 3), dtype='float64')
        B = np.zeros((3, 3), dtype='float64')
        D = np.zeros((3, 3), dtype='float64') # TODO: 2x2 matrix?
        z0 = self.z1
        z1 = self.z2

        if self.mid1_ref:
            Qbar1 = self.get_Qbar_matrix(self.mid1_ref, theta=0.)
            A += Qbar1 * thickness

        if self.mid2_ref:
            Qbar2 = self.get_Qbar_matrix(self.mid2_ref, theta=0.)
            D += Qbar2 * (z1 ** 3 - z0 ** 3) * self.twelveIt3
        #Qbar3 = self.get_Qbar_matrix(self.mid3_ref, theta=0.)
        #if self.mid4_ref:
            #Qbar4 = self.get_Qbar_matrix(self.mid4_ref, theta=0.)
            #B += Qbar4 * (z1 ** 2 - z0 ** 2)

        ts = self.tst * thickness

        # [N, M, Q].T =   [TG1, T^2 * G4, 0]              * [epsilon0]
                        # [T^2 * G4, T^3/12 * G2, 0]        [xi]
                        # [0, 0, Ts * G3]                   [gamma]

        #B += Qbar * thickness * zmean
        #D += Qbar * thickness * (z1i ** 3 - z0i ** 3)
        #N += Qbar * alpha * thickness
        #M += Qbar * alpha * thickness * zmean
        #B /= 2.
        #D /= 3.
        #M /= 2.
        #np.set_printoptions(linewidth=120, suppress=True)
        #print(ABD)
        return A, B, D

    def get_ABD_matrices(self, theta_offset=0.) -> np.ndarray:
        """
        Gets the ABD matrix

        Parameters
        ----------
        theta_offset : float
            rotates the ABD matrix; measured in degrees

        """
        A, B, D = self.get_individual_ABD_matrices(theta_offset=theta_offset)
        ABD = np.block([
            [A, B],
            [B, D],
        ])
        return ABD

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PQUAD1 pid=%s' % self.pid
        if self.mid1:
            self.mid1_ref = model.Material(self.mid1, msg)
        if self.mid2 and self.mid2 != -1:
            self.mid2_ref = model.Material(self.mid2, msg)
        if self.mid3:
            self.mid3_ref = model.Material(self.mid3, msg)
        if self.t_membrane is not None:
            z1 = abs(self.z1)
            z2 = abs(self.z2)
            t = self.t_membrane
            if not ((-1.5*t <= z1 <= 1.5*t) or (-1.5*t <= z2 <= 1.5*t)):
                msg = (f'PQUAD1 pid={self.pid} midsurface: z1={self.z1:g} z2={self.z2:g} t={t:g} '
                       f'not in range of -1.5t < zi < 1.5t')
                model.log.warning(msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid1 = self.Mid1()
        self.mid2 = self.Mid2()
        self.mid3 = self.Mid3()
        self.mid1_ref = None
        self.mid2_ref = None
        self.mid3_ref = None

    def raw_fields(self):
        list_fields = ['PQUAD1', self.pid, self.Mid1(), self.t_membrane, self.Mid2(),
                       self.inertia, self.Mid3(), self.t_shear, self.nsm, self.z1,
                       self.z2]
        return list_fields

    def repr_fields(self):
        #twelveIt3 = set_blank_if_default(self.twelveIt3, 1.0)
        nsm = set_blank_if_default(self.nsm, 0.0)
        #if self.t is not None:
            #t_over_2 = self.t / 2.
            #z1 = set_blank_if_default(self.z1, -t_over_2)
            #z2 = set_blank_if_default(self.z2, t_over_2)
        #else:
            #z1 = self.z1
            #z2 = self.z2

        mid1 = self.Mid1()
        mid2 = self.Mid2()
        mid3 = self.Mid3()
        mid1 = None if mid1 == 0 else mid1
        mid2 = None if mid2 == 0 else mid2
        mid3 = None if mid3 == 0 else mid3
        list_fields = ['PSHELL', self.pid, mid1, self.t_membrane, mid2,
                       self.inertia, mid3, self.t_shear, nsm, self.z1, self.z2]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
