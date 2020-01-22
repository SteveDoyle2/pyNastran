# pylint: disable=C0103,R0902,R0904,R0914,C0111
"""
All bush elements are defined in this file.  This includes:

 * CBUSH
 * CBUSH1D
 * CBUSH2D

All bush elements are BushElement and Element objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Element
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_double_or_blank, double_or_blank,
    string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class BushElement(Element):
    def __init__(self):
        self.cid = None
        Element.__init__(self)

    def Cid(self):
        if self.cid is None:
            return None
        elif isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def Mass(self):
        return 0.

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        return [tuple(sorted(self.node_ids))]

    #def Centroid(self):
        ## same as below, but we ignore the 2nd point it it's None
        #p = (self.nodes_ref[1].get_position() + self.nodes_ref[0].get_position()) / 2.

        ##p = self.nodes_ref[0].get_position()
        ##if self.nodes_ref[1] is not None:
            ##p += self.nodes_ref[1].get_position()
            ##p /= 2.
        #return p

    #def center_of_mass(self):
        #return self.Centroid()

class CBUSH(BushElement):
    """
    Generalized Spring-and-Damper Connection

    Defines a generalized spring-and-damper structural element that
    may be nonlinear or frequency dependent.

    +-------+-----+------+----+----+-------+----+----+-----+
    |   1   |  2  |  3   |  4 |  5 |   6   |  7 |  8 |  9  |
    +=======+=====+======+====+====+=======+====+====+=====+
    | CBUSH | EID | PID  | GA | GB | GO/X1 | X2 | X3 | CID |
    +-------+-----+------+----+----+-------+----+----+-----+
    |       |  S  | OCID | S1 | S2 |   S3  |    |    |     |
    +-------+-----+------+----+----+-------+----+----+-----+
    """
    type = 'CBUSH'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', 8:'cid', 9:'s', 10:'ocid'
    }
    _properties = ['_field_map', ]

    def update_by_cp_name(self, cp_name, value):
        #if isinstance(pname_fid, int):
            #self._update_field_helper(pname_fid, value)
        if cp_name == 'X1':
            self.x[0] = value
        elif cp_name == 'X2':
            self.x[1] = value
        elif cp_name == 'X3':
            self.x[2] = value

        elif cp_name == 'S1':
            self.si[0] = value
        elif cp_name == 'S2':
            self.si[1] = value
        elif cp_name == 'S3':
            self.si[2] = value

        elif cp_name == 'S':
            self.s = value
        else:
            raise NotImplementedError('element_type=%r has not implemented %r in update_by_cp_name' % (
                self.type, cp_name))

    def _update_field_helper(self, n, value):
        if n == 11:
            self.si[0] = value
        elif n == 12:
            self.si[1] = value
        elif n == 13:
            self.si[2] = value
        else:
            if self.g0 is not None:
                if n == 5:
                    self.g0 = value
                else:
                    raise KeyError('Field %r=%r is an invalid CBUSH entry.' % (n, value))
            else:
                if n == 5:
                    self.x[0] = value
                elif n == 6:
                    self.x[1] = value
                elif n == 7:
                    self.x[2] = value
                else:
                    raise KeyError('Field %r=%r is an invalid CBUSH entry.' % (n, value))

    def __init__(self, eid, pid, nids, x, g0, cid=None, s=0.5, ocid=-1, si=None, comment=''):
        """
        Creates a CBUSH card

        Parameters
        ----------
        eid : int
            Element id
        pid : int
            Property id (PBUSH)
        nids : List[int, int]
            node ids; connected grid points at ends A and B
            The nodes may be coincident, but then cid is required.
        x : List[float, float, float]; None
            List : the directional vector used to define the stiffnesses
                   or damping from the PBUSH card
            None : use g0
        g0 : int/None
            int : the directional vector used to define the stiffnesses
                  or damping from the PBUSH card
            None : use x
        cid : int; default=None
            Element coordinate system identification. A 0 means the basic
            coordinate system. If CID is blank, then the element coordinate
            system is determined from GO or Xi.
        s: float; default=0.5
            Location of spring damper (0 <= s <= 1.0)
        ocid : int; default=-1
            Coordinate system identification of spring-damper offset.
            (Integer > -1; Default = -1, which means the offset
            point lies on the line between GA and GB)
        si : List[float, float, float]; default=None
            Components of spring-damper offset in the OCID coordinate system
            if OCID > 0.
            None : [None, None, None]
        comment : str; default=''
            a comment for the card
        """
        BushElement.__init__(self)
        if comment:
            self.comment = comment
            #: if OCID > 0.
        if si is None:
            si = [None, None, None]
        if x is None:
            x = [None, None, None]

        self.eid = eid
        self.pid = pid
        self.nodes = nids
        self.x = x
        self.g0 = g0
        self.cid = cid
        self.s = s
        self.ocid = ocid
        self.si = si
        self.nodes_ref = None
        self.g0_ref = None
        self.pid_ref = None
        self.cid_ref = None
        self.ocid_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        x = []
        g0 = []
        cid = []
        s = []
        ocid = []
        si = []
        nan = np.full(3, np.nan)
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append([nid if nid is not None else 0 for nid in element.nodes])

            if element.cid is None:
                cid.append(-1)
                g0i = element.g0
                #print(g0i, element.x)
                if g0i is not None:
                    assert element.x[0] is None
                    x.append(nan)
                    g0.append(g0i)
                else:
                    #assert element.x[0] is not None, element.get_stats()
                    if element.x[0] is None:
                        x.append(nan)
                    else:
                        x.append(element.x)
                    g0.append(-1)
            else:
                cid.append(element.cid)
                g0i = element.g0
                if g0i is not None:
                    assert element.x[0] is None
                    x.append(nan)
                    g0.append(g0i)
                else:
                    if element.x[0] is None:
                        x.append(nan)
                    else:
                        x.append(element.x)
                    #assert element.x[0] is None, element.get_stats()
                    #x.append(nan)
                    g0.append(-1)

            s.append(element.s)
            ocid.append(element.ocid)
            if element.si[0] is None:
                si.append(nan)
            else:
                si.append(element.si)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('pid', data=pids)
        #print('x =', x)
        #print('g0 =', g0)
        #print('cid =', cid)
        h5_file.create_dataset('x', data=x)
        h5_file.create_dataset('g0', data=g0)
        h5_file.create_dataset('cid', data=cid)

        h5_file.create_dataset('s', data=s)
        h5_file.create_dataset('ocid', data=ocid)
        #print('si =', si)
        h5_file.create_dataset('si', data=si)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CBUSH card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer_or_blank(card, 4, 'gb')

        #: Element coordinate system identification. A 0 means the basic
        #: coordinate system. If CID is blank, then the element coordinate
        #: system is determined from GO or Xi.
        #: (default=blank=element-based)
        cid = integer_or_blank(card, 8, 'cid')

        x1_g0 = integer_double_or_blank(card, 5, 'x1_g0')
        if isinstance(x1_g0, integer_types):
            g0 = x1_g0
            x = None
        elif isinstance(x1_g0, float):
            g0 = None
            x1 = x1_g0
            x2 = double_or_blank(card, 6, 'x2', 0.0)
            x3 = double_or_blank(card, 7, 'x3', 0.0)
            x = [x1, x2, x3]
            if not isinstance(cid, integer_types):
                assert max(x) != min(x), 'x=%s' % x
        else:
            g0 = None
            x = [None, None, None]

        #: Location of spring damper (0 <= s <= 1.0)
        s = double_or_blank(card, 9, 's', 0.5)
        #: Coordinate system identification of spring-damper offset. See
        #: Remark 9. (Integer > -1; Default = -1, which means the offset
        #: point lies on the line between GA and GB
        ocid = integer_or_blank(card, 10, 'ocid', -1)
        #: Components of spring-damper offset in the OCID coordinate system
        #: if OCID > 0.
        si = [double_or_blank(card, 11, 's1'),
              double_or_blank(card, 12, 's2'),
              double_or_blank(card, 13, 's3')]
        assert len(card) <= 14, 'len(CBUSH card) = %i\ncard=%s' % (len(card), card)
        return CBUSH(eid, pid, [ga, gb], x, g0, cid=cid, s=s, ocid=ocid, si=si, comment=comment)

    @classmethod
    def add_op2_data(cls, data, f, comment=''):
        """
        Adds a CBUSH card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        ((eid, pid, ga, gb, cid, s, ocid, si), x, g0) = data
        return CBUSH(eid, pid, [ga, gb], x, g0, cid=cid, s=s, ocid=ocid, si=si, comment=comment)

    #@property
    #def nodes(self):
        #return [self.ga, self.gb]

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def _verify(self, xref):
        ga = self.Ga()
        gb = self.Gb()
        cid = self.Cid()
        ocid = self.OCid()
        pid = self.Pid()
        #si = self.si
        assert isinstance(ga, integer_types), 'ga=%r' % ga
        assert isinstance(gb, integer_types) or gb is None, 'gb=%r' % gb
        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(cid, integer_types) or cid is None, 'cid=%r' % cid
        assert isinstance(ocid, integer_types), 'ocid=%r' % ocid

    def Ga(self):
        if self.nodes_ref is not None:
            return self.nodes_ref[0].nid
        return self.nodes[0]

    def Gb(self):
        if self.nodes[1] in [0, None]:
            return 0
        if self.nodes_ref is not None:
            return self.nodes_ref[1].nid
        return self.nodes[1]

    def G0(self):
        if self.g0_ref is not None:
            return self.g0_ref.nid
        return self.g0

    def OCid(self):
        if self.ocid_ref is not None:
            return self.ocid_ref.cid
        return self.ocid

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CBUSH eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)
        if self.g0 is not None:
            self.g0_ref = model.Node(self.g0, msg=msg)
        if self.cid is not None:
            self.cid_ref = model.Coord(self.cid, msg=msg)
        if self.ocid is not None and self.ocid != -1:
            self.ocid_ref = model.Coord(self.ocid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CBUSH eid=%s' % self.eid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        if self.g0 is not None:
            self.g0_ref = model.EmptyNode(self.g0, msg=msg)
        if self.cid is not None:
            self.cid_ref = model.safe_coord(self.cid, self.eid, xref_errors, msg=msg)
        if self.ocid is not None and self.ocid != -1:
            self.ocid_ref = model.safe_coord(self.ocid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.ga = self.Ga()
        self.pid = self.Pid()
        self.cid = self.Cid()
        self.gb = self.Gb()
        self.g0 = self.G0()
        self.cid_ref = None
        self.ga_ref = None
        self.gb_ref = None
        self.g0_ref = None
        self.pid_ref = None
        self.ocid_ref = None

    def _get_x_g0(self):
        if self.g0 is not None:
            x = [self.G0(), None, None]
        else:
            x = self.x
        return x

    def raw_fields(self):
        x = self._get_x_g0()
        list_fields = (['CBUSH', self.eid, self.Pid(), self.Ga(), self.Gb()] + x +
                       [self.Cid(), self.s, self.ocid] + self.si)
        return list_fields

    def repr_fields(self):
        x = self._get_x_g0()

        ocid = set_blank_if_default(self.OCid(), -1)
        s = set_blank_if_default(self.s, 0.5)
        list_fields = (['CBUSH', self.eid, self.Pid(), self.Ga(), self.Gb()] +
                       x + [self.Cid(), s, ocid] + self.si)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CBUSH1D(BushElement):
    type = 'CBUSH1D'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', 5:'cid',
    }

    def __init__(self, eid, pid, nids, cid=None, comment=''):
        if comment:
            self.comment = comment
        BushElement.__init__(self)
        self.eid = eid
        self.pid = pid
        self.ga = nids[0]
        self.gb = nids[1]
        self.cid = cid
        self.ga_ref = None
        self.gb_ref = None
        self.cid_ref = None
        self.pid_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        cid = []
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append([nid if nid is not None else 0 for nid in element.nodes])
            cid.append(element.cid if element.cid is not None else -1)
        #h5_file.create_dataset('_comment', data=comments)
        #print('cid =', cid)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('cid', data=cid)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CBUSH1D card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer_or_blank(card, 4, 'gb')
        cid = integer_or_blank(card, 5, 'cid')
        assert len(card) <= 6, 'len(CBUSH1D card) = %i\ncard=%s' % (len(card), card)
        return CBUSH1D(eid, pid, [ga, gb], cid=cid, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #eid = data[0]
        #pid = data[1]
        #ga = data[2]
        #gb = data[3]
        #raise NotImplementedError(data)
        #return CBUSH1D(eid, pid, [ga, gb], cid, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CBUSH1D eid=%s' % self.eid
        self.ga_ref = model.Node(self.ga, msg=msg)
        if self.gb:
            self.gb_ref = model.Node(self.gb, msg=msg)
        self.pid_ref = model.Property(self.pid, msg=msg)
        if self.cid is not None:
            self.cid_ref = model.Coord(self.cid)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CBUSH1D eid=%s' % self.eid
        self.ga_ref = model.Node(self.ga, msg=msg)
        if self.gb:
            self.gb_ref = model.Node(self.gb, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        if self.cid is not None:
            self.cid_ref = model.safe_coord(self.cid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.cid = self.Cid()
        self.pid = self.Pid()
        self.ga_ref = None
        self.gb_ref = None
        self.cid_ref = None
        self.pid_ref = None

    def _verify(self, xref):
        ga = self.Ga()
        gb = self.Gb()
        cid = self.Cid()
        pid = self.Pid()
        assert isinstance(ga, integer_types), 'ga=%r' % ga
        assert isinstance(gb, integer_types) or gb is None, 'gb=%r' % gb
        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(cid, integer_types) or cid is None, 'cid=%r' % cid

    def Ga(self):
        if self.ga_ref is not None:
            return self.ga_ref.nid
        return self.ga

    def Gb(self):
        if self.gb_ref is not None:
            return self.gb_ref.nid
        return self.gb

    @property
    def nodes(self):
        return [self.ga, self.gb]

    @property
    def nodes_ref(self):
        return [self.ga_ref, self.gb_ref]

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def raw_fields(self):
        list_fields = ['CBUSH1D', self.eid, self.Pid(), self.Ga(), self.Gb(),
                       self.Cid()]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CBUSH2D(BushElement):
    """
    2-D Linear-Nonlinear Connection
    Defines the connectivity of a two-dimensional Linear-Nonlinear element.
    """
    type = 'CBUSH2D'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', 5:'cid', 6:'plane', 7:'sptid',
    }

    def __init__(self, eid, pid, nids, cid=0, plane='XY', sptid=None, comment=''):
        BushElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.ga = nids[0]
        self.gb = nids[1]
        self.cid = cid
        self.plane = plane
        self.sptid = sptid
        if self.plane not in ['XY', 'YZ', 'ZX']:
            msg = ("plane not in required list, plane=%r\n"
                   "expected planes = ['XY','YZ','ZX']" % self.plane)
            raise RuntimeError(msg)
        self.ga_ref = None
        self.gb_ref = None
        self.pid_ref = None
        self.cid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CBUSH2D card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid')
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')
        cid = integer_or_blank(card, 5, 'cid', 0)
        plane = string_or_blank(card, 6, 'plane', 'XY')
        sptid = integer_or_blank(card, 7, 'sptid')
        assert len(card) <= 8, 'len(CBUSH2D card) = %i\ncard=%s' % (len(card), card)
        return CBUSH2D(eid, pid, [ga, gb], cid, plane, sptid, comment=comment)

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids, encoding='ascii'):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        cid = []
        plane = []
        sptid = []
        #nan = np.full(3, np.nan)
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
            cid.append(element.cid)
            plane.append(element.plane.encode(encoding))
            sptidi = 0 if element.sptid is None else element.sptid
            sptid.append(sptidi)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('pid', data=pids)
        #print('x =', x)
        #print('g0 =', g0)
        #print('cid =', cid)
        h5_file.create_dataset('cid', data=cid)
        h5_file.create_dataset('plane', data=plane)
        #print('si =', si)
        h5_file.create_dataset('sptid', data=sptid)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #eid = data[0]
        #pid = data[1]
        #ga = data[2]
        #gb = data[3]
        #raise NotImplementedError(data)
        #return CBUSH2D(eid, pid, [ga, gb], cid, plane, sptid, comment=comment)

    def _verify(self, xref):
        ga = self.Ga()
        gb = self.Gb()
        cid = self.Cid()
        pid = self.Pid()
        plane = self.plane
        assert isinstance(ga, integer_types), 'ga=%r' % ga
        assert isinstance(gb, integer_types), 'gb=%r' % gb
        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(cid, integer_types), 'cid=%r' % cid
        assert self.plane in ['XY', 'YZ', 'ZX'], 'plane=%r' % plane

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        if isinstance(self.gb, integer_types):
            return self.gb
        return self.gb_ref.nid

    @property
    def nodes(self):
        return [self.ga, self.gb]

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CBUSH2D eid=%s' % self.eid
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        self.pid_ref = model.Property(self.pid)
        if self.cid is not None:
            self.cid_ref = model.Coord(self.cid, msg=msg)
        #if self.sptid is not None:
            #pass

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by CBUSH2D eid=%s' % self.eid
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)
        if self.cid is not None:
            self.cid_ref = model.safe_coord(self.cid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.cid = self.Cid()
        self.pid = self.Pid()
        self.ga_ref = None
        self.gb_ref = None
        self.cid_ref = None
        self.pid_ref = None

    def raw_fields(self):
        list_fields = ['CBUSH2D', self.eid, self.Pid(), self.Ga(), self.Gb(),
                       self.Cid(), self.plane, self.sptid]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)
