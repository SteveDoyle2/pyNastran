# pylint: disable=C0103
"""
All ungrouped elements are defined in this file.  This includes:

 * CFAST
 * CGAP
 * CRAC2D
 * CRAC3D
 * PLOTEL
 * GENEL

All ungrouped elements are Element objects.

"""
from __future__ import annotations
from typing import Tuple, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import (
    Element, BaseCard, break_word_by_trailing_integer, _node_ids)
from pyNastran.bdf.bdf_interface.assign_type import (
    fields, integer, integer_or_blank, integer_double_or_blank,
    double, double_or_blank, string)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class CFAST(Element):
    """defines the CFAST element"""
    type = 'CFAST'
    _properties = ['node_ids', 'nodes']
    _field_map = {
        1: 'eid', 2:'pid', 3:'Type', 4:'ida', 5:'idb', 6:'gs', 7:'ga', 8:'gb',
        9:'xs', 10:'ys', 11:'zs',
    }
    cp_name_map = {
    }

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        Type = 'Type'
        ida = 1
        idb = 2
        return CFAST(eid, pid, Type, ida, idb,
                     gs=None, ga=None, gb=None, xs=None, ys=None, zs=None, comment='')

    def __init__(self, eid, pid, Type, ida, idb, gs=None, ga=None, gb=None,
                 xs=None, ys=None, zs=None, comment=''):
        Element.__init__(self)
        if comment:
            self.comment = comment
        if pid is None:
            pid = eid
        self.eid = eid
        self.pid = pid
        self.Type = Type
        self.ida = ida
        self.idb = idb
        self.gs = gs
        self.ga = ga
        self.gb = gb
        self.xs = xs
        self.ys = ys
        self.zs = zs
        self.pid_ref = None
        self.gs_ref = None
        self.ga_ref = None
        self.gb_ref = None

    def validate(self):
        if self.Type not in ['PROP', 'ELEM']:
            msg = 'CFAST; eid=%s Type=%r must be in [PROP, ELEM]' % (self.eid, self.Type)
            raise TypeError(msg)

        gab_is_none = self.ga is None or self.gb is None
        xyz_is_none = self.xs is None or self.ys is None or self.zs is None
        if self.gs is None and gab_is_none and xyz_is_none:
            msg = ('CFAST; eid=%s; gs=%s is not an integer or\n'
                   '              [ga=%s, gb=%s] are not integers or \n'
                   '              [xs=%s, ys=%s, zs=%s] are not floats' % (
                       self.eid, self.gs, self.ga, self.gb, self.xs, self.ys, self.zs))
            raise ValueError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CFAST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        Type = string(card, 3, 'Type')
        ida = integer(card, 4, 'ida')
        idb = integer(card, 5, 'idb')
        gs = integer_or_blank(card, 6, 'gs')
        ga = integer_or_blank(card, 7, 'ga')
        gb = integer_or_blank(card, 8, 'gb')
        xs = double_or_blank(card, 9, 'xs')
        ys = double_or_blank(card, 10, 'ys')
        zs = double_or_blank(card, 11, 'zs')
        assert len(card) <= 12, 'len(CFAST card) = %i\ncard=%s' % (len(card), card)
        #if self.Type=='PROP': # PSHELL/PCOMP  ida & idb
        return CFAST(eid, pid, Type, ida, idb, gs=gs, ga=ga, gb=gb,
                     xs=xs, ys=ys, zs=zs, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CFAST eid=%s' % self.eid
        self.pid_ref = model.Property(self.Pid(), msg=msg)
        if self.gs:
            self.gs_ref = model.Node(self.Gs(), msg=msg)
        if self.ga:
            self.ga_ref = model.Node(self.Ga(), msg=msg)
        if self.gb:
            self.gb_ref = model.Node(self.Gb(), msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        self.gs = self.Gs()
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.pid_ref = None
        self.gs_ref = None
        self.ga_ref = None
        self.gb_ref = None

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

    def raw_fields(self):
        gs, xs, ys, zs = self._gs_xyz()
        list_fields = ['CFAST', self.eid, self.Pid(), self.Type, self.ida, self.idb,
                       gs, self.Ga(), self.Gb(), xs, ys, zs]
        return list_fields

    @property
    def nodes(self):
        """gets all the nodes used on the CFAST (Gs, Ga, Gb)"""
        return [self.gs, self.ga, self.gb]

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def _gs_xyz(self) -> Tuple[Optional[int], Optional[float], Optional[float], Optional[float]]:
        if self.gs is None:
            out = (self.gs, self.xs, self.ys, self.zs)
        else:
            out = (self.Gs(), None, None, None)
        return out

    def Gs(self):
        """Gets the GS node"""
        if isinstance(self.gs, integer_types):
            return self.gs
        elif self.gs is not None:
            return self.gs_ref.nid
        # If neither GS nor GA is specified, then (XS, YS, ZS) in basic must be specified.
        raise RuntimeError(f'Gs was not returned from CFAST\n{self.get_stats()}')

    def Ga(self):
        """Gets the GA node"""
        if isinstance(self.ga, integer_types) or self.ga is None:
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        """Gets the GB node"""
        if isinstance(self.gb, integer_types) or self.gb is None:
            return self.gb
        return self.gb_ref.nid

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced

        """
        pass

    @property
    def node_ids(self):
        """gets all the node ids used on the CFAST (Gs, Ga, Gb)"""
        return [self.Gs(), self.Ga(), self.Gb()]

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CGAP(Element):
    """
    +------+-----+-----+-----+-----+-----+-----+------+-----+
    |  1   |  2  |  3  |  4  |  5  |  6  |  7  |   8  |  9  |
    +======+=====+=====+=====+=====+=====+=====+======+=====+
    | CGAP | EID | PID | GA  | GB  | X1  | X2  |  X3  | CID |
    +------+-----+-----+-----+-----+-----+-----+------+-----+
    | CGAP | 17  |  2  | 110 | 112 | 5.2 | 0.3 | -6.1 |     |
    +------+-----+-----+-----+-----+-----+-----+------+-----+

    or

    +------+-----+-----+-----+-----+-----+-----+------+-----+
    |  1   |  2  |  3  |  4  |  5  |  6  |  7  |   8  |  9  |
    +======+=====+=====+=====+=====+=====+=====+======+=====+
    | CGAP | EID | PID | GA  | GB  | GO  |     |      | CID |
    +------+-----+-----+-----+-----+-----+-----+------+-----+
    | CGAP | 17  |  2  | 110 | 112 | 13  |     |      |     |
    +------+-----+-----+-----+-----+-----+-----+------+-----+

    """
    type = 'CGAP'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb',
    }
    def update_by_cp_name(self, cp_name, value):
        if isinstance(cp_name, int):
            #self._update_field_helper(cp_name, value)
            raise NotImplementedError('element_type=%r has not implemented %r in cp_name_map' % (
                self.type, cp_name))
        #elif cp_name == 'Z0':
            #self.z0 = value
        #elif cp_name == 'SB':
            #self.sb = value
        #elif cp_name == 'TREF':
            #self.tref = value
        #elif cp_name == 'GE':
            #self.ge = value
        elif cp_name.startswith('X'):
            word, num = break_word_by_trailing_integer(cp_name)
            num = int(num)
            if word == 'X':
                self.x[num - 1] = value
            else:
                raise RuntimeError('eid=%s cp_name=%r word=%s\n' % (self.eid, cp_name, word))
        else:
            raise NotImplementedError('element_type=%r has not implemented %r in cp_name_map' % (
                self.type, cp_name))

    def __init__(self, eid, pid, nids, x, g0, cid=None, comment=''):
        """
        Creates a CGAP card

        Parameters
        ----------
        eid : int
            Element ID
        pid : int
            Property ID (PGAP)
        nids : List[int, int]
            node ids; connected grid points at ends A and B
        x : List[float, float, float]
            Components of the orientation vector,
            from GA, in the displacement coordinate system at GA
        g0 : int
            GO Alternate method to supply the orientation vector using
            grid point GO. Direction of is from GA to GO
        cid : int; default=None
            Element coordinate system identification number.
            CID must be specified if GA and GB are coincident
            (distance from GA to GB < 10^-4)
        comment : str; default=''
            a comment for the card

        """
        Element.__init__(self)
        if comment:
            self.comment = comment
        if pid is None:
            pid = eid
        self.eid = eid
        self.pid = pid
        self.ga = nids[0]
        self.gb = nids[1]
        self.x = x
        self.g0 = g0
        self.cid = cid
        self.ga_ref = None
        self.gb_ref = None
        self.g0_ref = None
        self.cid_ref = None
        self.pid_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        #comments = []
        pids = []
        nodes = []
        x = []
        g0 = []
        cid = []
        nan = np.full(3, np.nan)
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)

            if element.cid is None:
                cid.append(-1)
                g0i = element.g0
                if g0i is not None:
                    #assert element.x[0] is None, element.get_stats()
                    x.append(nan)
                    g0.append(g0i)
                else:
                    if element.x[0] is None:
                        x.append(nan)
                    else:
                        x.append(element.x)
                    #assert element.x[0] is not None, element.get_stats()
                    #x.append(element.x)
                    g0.append(-1)
            else:
                cid.append(element.cid)
                g0i = element.g0
                if g0i is not None:
                    assert x[0] is None
                    x.append(nan)
                    g0.append(g0i)
                else:
                    if element.x[0] is None:
                        x.append(nan)
                    else:
                        x.append(element.x)
                    #assert element.x[0] is None, element.get_stats()
                    x.append(nan)
                    g0.append(-1)
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

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CGAP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer_or_blank(card, 3, 'ga')
        gb = integer_or_blank(card, 4, 'gb')
        x1_g0 = integer_double_or_blank(card, 5, 'x1_g0')
        cid = integer_or_blank(card, 8, 'cid')
        if isinstance(x1_g0, integer_types):
            g0 = x1_g0
            x = None
        elif isinstance(x1_g0, float):
            g0 = None
            x1 = x1_g0
            x2 = double_or_blank(card, 6, 'x2', 0.0)
            x3 = double_or_blank(card, 7, 'x3', 0.0)
            x = [x1, x2, x3]
        else:
            #raise RuntimeError('invalid CGAP...x1/g0 = %r' %(x1_g0))
            g0 = None
            x = [None, None, None]
        assert len(card) <= 9, 'len(CGAP card) = %i\ncard=%s' % (len(card), card)
        return CGAP(eid, pid, [ga, gb], x, g0, cid=cid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CGAP card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        pid = data[1]
        ga = data[2]
        gb = data[3]
        g0 = data[4]
        x1 = data[5]
        x2 = data[6]
        x3 = data[7]
        x = [x1, x2, x3]
        cid = data[8]
        if cid == -1:
            cid = None
        return CGAP(eid, pid, [ga, gb], x, g0, cid=cid, comment=comment)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced

        """
        cid = self.Cid()
        eid = self.eid
        pid = self.Pid()
        nids = self.node_ids

        assert cid is None or isinstance(cid, integer_types), 'cid=%r\n%s' % (cid, str(self))
        assert isinstance(eid, integer_types), 'eid=%r\n%s' % (eid, str(self))
        assert isinstance(pid, integer_types), 'pid=%r\n%s' % (pid, str(self))
        for i, nid in enumerate(nids):
            assert isinstance(nid, integer_types), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PGAP'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            if self.cid is not None and self.cid != 0:
                assert self.cid_ref.type in ['CORD1R', 'CORD1C', 'CORD1S', 'CORD2R', 'CORD2C',
                                             'CORD2S'], 'cid=%i self.cid.type=%s' % (
                                                 cid, self.cid_ref.type)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CGAP eid=%s' % self.eid
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        if self.g0:
            self.g0_ref = model.Node(self.g0, msg=msg)
            self.x = self.g0_ref.get_position()
        self.pid_ref = model.Property(self.pid, msg=msg)
        if self.cid:
            self.cid_ref = model.Coord(self.cid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CGAP eid=%s' % self.eid
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        if self.g0:
            self.g0_ref = model.Node(self.g0, msg=msg)
            self.x = self.g0_ref.get_position()
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

    @property
    def nodes(self):
        return [self.ga, self.gb]

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def Cid(self):
        if self.cid_ref is None:
            return self.cid
        return self.cid_ref.cid

    def Ga(self):
        if self.ga_ref is None:
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        if self.gb_ref is None:
            return self.gb
        return self.gb_ref.nid

    def G0(self):
        if isinstance(self.g0, integer_types):
            return self.g0
        return self.g0_ref.nid

    def raw_fields(self):
        if self.g0 is not None:
            x = [self.G0(), None, None]
        else:
            x = self.x
        list_fields = (['CGAP', self.eid, self.Pid(), self.Ga(), self.Gb()] + x +
                       [self.Cid()])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CrackElement(Element):
    type = 'Crack'

    def __init__(self):
        Element.__init__(self)
        self.eid = 0
        self.nodes_ref = None
        self.pid_ref = None

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s eid=%s' % (self. type, self.eid)
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
        msg = ', which is required by %s eid=%s' % (self. type, self.eid)
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None


class CRAC2D(CrackElement):
    type = 'CRAC2D'
    _properties = ['node_ids']
    _field_map = {
        1: 'eid', 2:'pid',
    }
    ## todo:: not done

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        nids = [1]
        return CRAC2D(eid, pid, nids, comment='')

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        if isinstance(self.nodes, np.ndarray):
            self.nodes = self.nodes.tolist()
        self.nodes = [None if np.isnan(nid) else nid
                      for nid in self.nodes]

    def __init__(self, eid, pid, nids, comment=''):
        CrackElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid

        nnodes = len(nids)
        if nnodes < 18:
            nids = nids + [None] * (18-nnodes)

        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 18, 'eid=%s nnodes=%s' % (self.eid, len(self.nodes))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CRAC2D card from ``BDF.add_card(...)``

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
            integer(card, 3, 'n1'), integer(card, 4, 'n2'),
            integer(card, 5, 'n3'), integer(card, 6, 'n4'),
            integer(card, 7, 'n5'), integer(card, 8, 'n6'),
            integer(card, 9, 'n7'), integer(card, 10, 'n8'),
            integer(card, 11, 'n9'), integer(card, 12, 'n10'),
            integer_or_blank(card, 13, 'n11'),
            integer_or_blank(card, 14, 'n12'),
            integer_or_blank(card, 15, 'n13'),
            integer_or_blank(card, 16, 'n14'),
            integer_or_blank(card, 17, 'n15'),
            integer_or_blank(card, 18, 'n16'),
            integer_or_blank(card, 19, 'n17'),
            integer_or_blank(card, 20, 'n18')
        ]
        assert len(card) <= 21, 'len(CRAC2D card) = %i\ncard=%s' % (len(card), card)
        return CRAC2D(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CRAC2D card from the OP2

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
        return CRAC2D(eid, pid, nids, comment=comment)

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        unused_nids = self.node_ids

        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)

    def get_edge_ids(self):
        return []

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CRAC2D', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CRAC3D(CrackElement):
    type = 'CRAC3D'
    _properties = ['node_ids']
    _field_map = {
        1: 'eid', 2:'pid',
    }
    ## todo:: not done

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        nids = [1]
        return CRAC3D(eid, pid, nids, comment='')

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        if isinstance(self.nodes, np.ndarray):
            self.nodes = self.nodes.tolist()
        self.nodes = [None if np.isnan(nid) else nid
                      for nid in self.nodes]

    def __init__(self, eid, pid, nids, comment=''):
        CrackElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        nnodes = len(nids)
        if nnodes < 64:
            nids = nids + [None] * (64-nnodes)
        self.nodes = self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 64, 'eid=%s nnodes=%s' % (self.eid, len(self.nodes))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CRAC3D card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        # required 1-10, 19-28
        # optional 11-18, 29-36, 37-64
        # all/none 37-46
        nids = fields(integer_or_blank, card, 'nid', 3, 67)  # cap at +3 = 67
        assert len(card) <= 67, 'len(CRAC3D card) = %i\ncard=%s' % (len(card), card)
        return CRAC3D(eid, pid, nids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a CRAC3D card from the OP2

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
        return CRAC3D(eid, pid, nids, comment=comment)

    def _verify(self, xref):
        eid = self.eid
        pid = self.Pid()
        unused_nids = self.node_ids
        assert isinstance(eid, integer_types)
        assert isinstance(pid, integer_types)

    def get_edge_ids(self):
        return []

    @property
    def node_ids(self):
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CRAC3D', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PLOTEL(BaseCard):
    """
    Defines a 1D dummy element used for plotting.

    This element is not used in the model during any of the solution
    phases of a problem. It is used to simplify plotting of
    structures with large numbers of colinear grid points, where the
    plotting of each grid point along with the elements connecting
    them would result in a confusing plot.

    +--------+-----+-----+-----+
    |   1    |  2  |  3  |  4  |
    +========+=====+=====+=====+
    | PLOTEL | EID | G1  | G2  |
    +--------+-----+-----+-----+

    """
    type = 'PLOTEL'
    _field_map = {
        1: 'eid', 3:'g1', 4:'g2',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nodes = [1, 2]
        return PLOTEL(eid, nodes, comment='')

    @classmethod
    def export_to_hdf5(cls, h5_file, model, encoding):
        """exports the elements in a vectorized way"""
        #comments = []
        nodes = []
        eids = list(model.plotels.keys())
        for eid in eids:
            element = model.plotels[eid]
            #comments.append(element.comment)
            nodes.append(element.nodes)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('nodes', data=nodes)

    def __init__(self, eid, nodes, comment=''):
        """
        Adds a PLOTEL card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : List[int, int]
            Unique GRID point IDs

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.nodes = nodes
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card, icard: int, comment=''):
        """
        Adds a PLOTEL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        offset = icard * 4
        eid = integer(card, 1+offset, 'eid')
        nodes = [
            integer(card, 2+offset, 'g1'),
            integer(card, 3+offset, 'g2'),
        ]
        #assert len(card) <= 4, 'len(PLOTEL card) = %i\ncard=%s' % (len(card), card)
        assert len(card) <= 8, 'len(PLOTEL card) = %i\ncard=%s' % (len(card), card)
        return PLOTEL(eid, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PLOTEL card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        nodes = [data[1], data[2]]
        return PLOTEL(eid, nodes, comment=comment)

    def _verify(self, xref):
        pass

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PLOTEL eid=%s' % self.eid
        self.nodes_ref = [
            model.Node(self.nodes[0], msg=msg),
            model.Node(self.nodes[1], msg=msg),
        ]

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self):
        node_idsi = self.nodes
        n1, n2 = node_idsi
        if not isinstance(n1, integer_types):
            node_idsi[0] = n1.Nid()
        if not isinstance(n2, integer_types):
            node_idsi[1] = n2.Nid()
        return node_idsi

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def raw_fields(self):
        list_fields = ['PLOTEL', self.eid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        msg = 'PLOTEL  %8i%8i%8i\n' % (self.eid, nodes[0], nodes[1])
        return self.comment + msg


class GENEL(BaseCard):
    """
    +-------+------+-----+------+------+------+------+-------+------+
    |   1   |   2  |  3  |   4  |   5  |   6  |   7  |   8   |   9  |
    +=======+======+=====+======+======+======+======+=======+======+
    | GENEL | EID  |     | UI1  | CI1  | UI2  | CI2  |  UI3  | CI3  |
    +-------+------+-----+------+------+------+------+-------+------+
    |       | UI4  | CI4 | UI5  | CI5  | etc. |      |       |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       | UD   |     | UD1  | CD1  | UD2  | CD2  | etc.  |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       | K/Z  |     | KZ11 | KZ21 | KZ31 | etc. | KZ22  | KZ32 |
    +-------+------+-----+------+------+------+------+-------+------+
    |       | etc. |     | KZ33 | KZ43 | etc. |      |       |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       | S    |     | S11  | S12  | etc. |  S21 |  etc. |      |
    +-------+------+-----+------+------+------+------+-------+------+

    +-------+------+-----+------+------+------+------+-------+------+
    | GENEL |  629 |     |  1   |  1   |  13  |  4   |   42  |   0  |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  24  |  2  |      |      |      |      |       |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  UD  |     |  6   |  2   |  33  |  0   |       |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  Z   | 1.0 | 2.0  | 3.0  | 4.0  | 5.0  |  6.0  | 7.0  |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  8.0 | 9.0 | 10.0 |      |      |      |       |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  S   | 1.5 | 2.5  | 3.5  | 4.5  | 5.5  |  6.5  | 7.5  |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  8.5 |     |      |      |      |      |       |      |
    +-------+------+-----+------+------+------+------+-------+------+

    """
    type = 'GENEL'
    #pid = 0
    _properties = ['node_ids', 'ul_nodes', 'ud_nodes', 'nodes']
    @classmethod
    def _init_from_empty(cls):
        eid = 1
        ul = None
        ud = None
        k = [1.]
        z = None
        return GENEL(eid, ul, ud, k, z, s=None, comment='')

    def __init__(self, eid, ul, ud, k, z, s=None, comment=''):
        """creates a GENEL card

        The required input is the {UL} list and the lower triangular
        portion of [K] or [Z].  Additional input may include the {UD}
        list and [S].  If [S] is input, must also be input.  If {UD} is
        input but [S] is omitted, [S] is internally calculated. In this
        case, {UD} must contain six and only six degrees-of freedom.
        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.ul = ul
        self.ud = ud
        if k is not None:
            self.k = np.asarray(k)
            self.z = None
        else:
            self.z = np.asarray(z)
            self.k = None

        if s is not None:
            s = np.asarray(s)
        self.s = s

        self.ul_nodes_ref = None
        self.ud_nodes_ref = None

    def _finalize_hdf5(self, encoding):
        self.ul = np.array(self.ul, dtype='int32')#.reshape(len(self.ul) // 2, 2)
        self.ud = np.array(self.ud, dtype='int32')#.reshape(len(self.ud) // 2, 2)

        if self.k is None or (isinstance(self.k, list) and len(self.k) == 0):
            self.k = None
        else:
            self.k = np.array(self.k)

        if self.z is None or (isinstance(self.z, list) and len(self.z) == 0):
            self.z = None
        else:
            self.z = np.array(self.z)

        if self.s is None or (isinstance(self.s, list) and len(self.s) == 0):
            self.s = None
        else:
            self.s = np.array(self.s)

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        card_fields = card.fields()
        ucard_fields = [field.upper() if field is not None else None
                        for field in card_fields]
        nfields = card.nfields

        ul = []
        ud = []
        n_ud = ucard_fields.count('UD')

        nk = ucard_fields.count('K')
        nz = ucard_fields.count('Z')
        ns = ucard_fields.count('S')
        assert n_ud in [0, 1], 'n_UD=%s fields=%s' % (n_ud, card_fields)

        assert nk in [0, 1], 'n_K=%s fields=%s' % (nz, card_fields)
        assert nz in [0, 1], 'n_Z=%s fields=%s' % (nk, card_fields)
        assert ns in [0, 1], 'n_S=%s fields=%s' % (ns, card_fields)

        i_ul = 3
        _ul_fields, unused_istop = _read_genel_fields_until_char_blank(ucard_fields, i_ul)
        for i, _ul in enumerate(_ul_fields):
            uli = integer(card, i + i_ul, 'UL_%i' % (i + 1))
            ul.append(uli)

        if n_ud:
            i_ud = ucard_fields.index('UD')
            if i_ud < nfields and ucard_fields[i_ud] == 'UD':
                assert ucard_fields[i_ud] == 'UD', fields
                _ud_fields, unused_istop = _read_genel_fields_until_char_blank(ucard_fields, i_ud+2)
                for i, _ud in enumerate(_ud_fields):
                    udi = integer(card, i + i_ud+2, 'UD_%i' % (i + 1))
                    ud.append(udi)

        k = None
        z = None
        s = None
        if nk:
            k = []
            ik = ucard_fields.index('K')
            assert ucard_fields[ik] == 'K', card_fields
            _k_fields, unused_istop = _read_genel_fields_until_char_blank(ucard_fields, ik+1)
            for i, _k in enumerate(_k_fields):
                ki = double(card, i + ik+1, 'K_%i' % (i + 1))
                k.append(ki)
            unused_nblanks = _get_genel_offset(nk)
            #kz = k

        if nz:
            z = []
            assert k is None, k
            iz = card_fields.index('Z')
            assert card_fields[iz] == 'Z', card_fields
            _z_fields, unused_istop = _read_genel_fields_until_char_blank(ucard_fields, iz+1)
            for i, _z in enumerate(_z_fields):
                zi = double(card, i + iz+1, 'Z_%i' % (i + 1))
                z.append(zi)
            unused_nblanks = _get_genel_offset(nz)
            #kz = z

        if ns:
            s = []
            i_s = ucard_fields.index('S')
            assert ucard_fields[i_s] == 'S', card_fields
            _s_fields, unused_istop = _read_genel_fields_until_char_blank(ucard_fields, i_s+1)
            for i, _s in enumerate(_s_fields):
                si = double(card, i + i_s+1, 'S_%i' % (i + 1))
                s.append(si)
            unused_nblanks = _get_genel_offset(ns)

        #---------------------------------
        ul = np.array(ul).reshape(len(ul) // 2, 2)
        ud = np.array(ud).reshape(len(ud) // 2, 2)
        return GENEL(eid, ul, ud, k, z, s, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by GENEL eid=%s' % self.eid
        self.ul_nodes_ref = model.Nodes(self.ul[:, 0], msg=msg)
        if len(self.ud):
            self.ud_nodes_ref = model.Nodes(self.ud[:, 0], msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        #msg = ', which is required by GENEL eid=%s' % self.eid
        self.cross_reference(model)
        #self.ga_ref = model.Node(self.ga, msg=msg)
        #self.gb_ref = model.Node(self.gb, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.ul[:, 0] = self.ul_nodes
        self.ud[:, 0] = self.ud_nodes
        self.ul_nodes_ref = None
        self.ud_nodes_ref = None

    #def center_of_mass(self):
        #return 0.0

    def raw_fields(self):
        # we add 2 to represent the GENEL,eid fields
        n_ul = self.ul.shape[0] * 2 + 2
        ul_nones = _get_genel_offset(n_ul) * [None]

        # same as UL for UD
        n_ud = self.ud.shape[0] * 2 + 2
        ud_nones = _get_genel_offset(n_ud) * [None]

        # we call this kz to simplify our life
        kz_char = 'K' if self.k is not None else 'Z'
        kz = self.k if self.k is not None else self.z

        # K/Z has a +1 instead of +2 because there is no blank after K/Z
        n_kz = len(kz) + 1
        kz_nones = _get_genel_offset(n_kz) * [None]


        ud_line = []
        if self.ud_nodes_ref is None:
            if len(self.ud):
                ud_nodes_dofs = self.ud.ravel().tolist()
                ud_line = ['UD', None] + self.ud.ravel().tolist() + ud_nones
        else:
            ud_nodes = self.ud_nodes
            ud_dofs = self.ud[:, 1]
            ud_nodes_dofs = []
            for ud_node, ud_dof in zip(ud_nodes, ud_dofs):
                ud_nodes_dofs.extend([ud_node, ud_dof])
            ud_line = ['UD', None] + ud_nodes_dofs + ud_nones

        #print('s = %r' % self.s, self.s is not None)
        s_line = []
        if self.s is not None and len(self.s):
            s_line = ['S'] + self.s.tolist()

        if self.ul_nodes_ref is None:
            ul_nodes_dofs = self.ul.ravel().tolist()
        else:
            ul_nodes = self.ul_nodes
            ul_dofs = self.ul[:, 1]
            ul_nodes_dofs = []
            for ul_node, ul_dof in zip(ul_nodes, ul_dofs):
                ul_nodes_dofs.extend([ul_node, ul_dof])
        ul_line = ul_nodes_dofs + ul_nones

        list_fields = ['GENEL', self.eid, None] + (
            ul_line +
            ud_line +
            [kz_char] + kz.tolist() + kz_nones +
            s_line
        )
        return list_fields

    @property
    def nodes(self):
        return self.node_ids

    @property
    def node_ids(self):
        nodes = self.ul[:, 0].tolist()
        if len(self.ud):
            nodes += self.ud[:, 0].tolist()
        return nodes

    @property
    def ul_nodes(self):
        """gets the {UL} nodes"""
        if self.ul_nodes_ref is None:
            return self.ul[:, 0]
        nodes = _node_ids(self, nodes=self.ul_nodes_ref, allow_empty_nodes=False)
        return nodes

    @property
    def ud_nodes(self):
        """gets the {UD} nodes"""
        if self.ud_nodes_ref is None:
            return self.ud[:, 0]
        nodes = _node_ids(self, nodes=self.ud_nodes_ref, allow_empty_nodes=False)
        return nodes

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

def _get_genel_offset(n_ul):
    """we add to to represent the GENEL,eid fields"""
    n_ul_leftover = n_ul % 8
    return 8 - n_ul_leftover


def _read_genel_fields_until_char_blank(card_fields, istart):
    """somewhat loose parser helper function for GENEL"""
    new_fields = []
    i = 0
    for i, field in enumerate(card_fields[istart:]):
        if field is None:
            break
        if field.upper() in ['UD', 'K', 'S', 'Z']:
            break
        new_fields.append(field)
    return new_fields, istart+i
