# pylint: disable=C0103
"""
All ungrouped elements are defined in this file.  This includes:

 * CFAST
 * CGAP
 * CRAC2D
 * CRAC3D
 * PLOTEL

All ungrouped elements are Element objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import Element, BaseCard, break_word_by_trailing_integer
from pyNastran.bdf.bdf_interface.assign_type import (
    fields, integer, integer_or_blank, integer_double_or_blank,
    double_or_blank, string)
from pyNastran.bdf.field_writer_8 import print_card_8


class CFAST(Element):
    type = 'CFAST'
    _field_map = {
        1: 'eid', 2:'pid', 3:'Type', 4:'ida', 5:'idb', 6:'gs', 7:'ga', 8:'gb',
        9:'xs', 10:'ys', 11:'zs',
    }
    cp_name_map = {
    }

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
        if(self.gs is None and
           (self.ga is None or self.gb is None) and
           (self.xs is None or self.ys is None or self.zs is None)):
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

    def cross_reference(self, model):
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

    def uncross_reference(self):
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
        list_fields = ['CFAST', self.eid, self.Pid(), self.Type, self.ida, self.idb,
                       self.Gs(), self.Ga(), self.Gb(), self.xs, self.ys, self.zs]
        return list_fields

    @property
    def nodes(self):
        return [self.gs, self.ga, self.gb]

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def Gs(self):
        if isinstance(self.gs, integer_types):
            return self.gs
        elif self.gs is not None:
            return self.gs_ref.nid

    def Ga(self):
        if isinstance(self.ga, integer_types) or self.ga is None:
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
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
        return [self.Gs(), self.Ga(), self.Gb()]

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
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
        self.cid_ref = None
        self.pid_ref = None

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

    def cross_reference(self, model):
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

    def uncross_reference(self):
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

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CrackElement(Element):
    type = 'Crack'

    def __init__(self):
        self.eid = 0
        self.nodes_ref = None
        self.pid_ref = None

    def cross_reference(self, model):
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

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        self.nodes_ref = None
        self.pid_ref = None


class CRAC2D(CrackElement):
    type = 'CRAC2D'
    _field_map = {
        1: 'eid', 2:'pid',
    }
    ## todo:: not done

    def __init__(self, eid, pid, nids, comment=''):
        CrackElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid

        nnodes = len(nids)
        if nnodes < 18:
            nids = nids + [None] * (18-nnodes)

        self.prepare_node_ids(nids, allow_empty_nodes=True)
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

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CRAC3D(CrackElement):
    type = 'CRAC3D'
    _field_map = {
        1: 'eid', 2:'pid',
    }
    ## todo:: not done

    def __init__(self, eid, pid, nids, comment=''):
        CrackElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        nnodes = len(nids)
        if nnodes < 64:
            nids = nids + [None] * (64-nnodes)
        self.prepare_node_ids(nids, allow_empty_nodes=True)
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

    def write_card(self, size=8, is_double=False):
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
    def add_card(cls, card, comment=''):
        """
        Adds a PLOTEL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        nodes = [
            integer(card, 2, 'g1'),
            integer(card, 3, 'g2'),
        ]
        assert len(card) <= 4, 'len(PLOTEL card) = %i\ncard=%s' % (len(card), card)
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

    def cross_reference(self, model):
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

    def uncross_reference(self):
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

    def write_card(self, size=8, is_double=False):
        nodes = self.node_ids
        msg = 'PLOTEL  %8i%8i%8i\n' % (self.eid, nodes[0], nodes[1])
        return self.comment + msg
