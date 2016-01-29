# pylint: disable=C0103,R0902,R0904,R0914,C0111
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

from pyNastran.utils import integer_types
from pyNastran.bdf.cards.base_card import Element, BaseCard
from pyNastran.bdf.bdfInterface.assign_type import (fields, integer, integer_or_blank,
    integer_double_or_blank, double_or_blank, string)  # double
from pyNastran.bdf.field_writer_8 import print_card_8


class CFAST(Element):
    type = 'CFAST'
    _field_map = {
        1: 'eid', 2:'pid', 3:'Type', 4:'ida', 5:'idb', 6:'gs', 7:'ga', 8:'gb',
        9:'xs', 10:'ys', 11:'zs',
    }

    def __init__(self):
        Element.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.eid = integer(card, 1, 'eid')
        self.pid = integer_or_blank(card, 2, 'pid', self.eid)
        self.Type = string(card, 3, 'Type')
        self.ida = integer(card, 4, 'ida')
        self.idb = integer(card, 5, 'idb')
        self.gs = integer_or_blank(card, 6, 'gs')
        self.ga = integer_or_blank(card, 7, 'ga')
        self.gb = integer_or_blank(card, 8, 'gb')
        self.xs = double_or_blank(card, 9, 'xs')
        self.ys = double_or_blank(card, 10, 'ys')
        self.zs = double_or_blank(card, 11, 'zs')
        assert len(card) <= 12, 'len(CFAST card) = %i' % len(card)
        #if self.Type=='PROP': # PSHELL/PCOMP  ida & idb

    def cross_reference(self, model):
        msg = ' which is required by CFAST eid=%s' % self.eid
        self.pid = model.Property(self.Pid(), msg=msg)
        self.pid_ref = self.pid
        if self.gs:
            self.gs = model.Node(self.Gs(), msg=msg)
            self.gs_ref = self.gs
        if self.ga:
            self.ga = model.Node(self.Ga(), msg=msg)
            self.ga_ref = self.ga
        if self.gb:
            self.gb = model.Node(self.Gb(), msg=msg)
            self.gb_ref = self.gb

    def uncross_reference(self):
        self.gs = self.Gs()
        self.ga = self.Ga()
        self.gb = self.Gb()
        del self.gs_ref, self.ga_ref, self.gb_ref

    def raw_fields(self):
        list_fields = ['CFAST', self.eid, self.Pid(), self.Type, self.ida, self.idb,
                       self.Gs(), self.Ga(), self.Gb(), self.xs, self.ys, self.zs]
        return list_fields

    @property
    def nodes(self):
        return [self.ga, self.gb]

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def Gs(self):
        if isinstance(self.gs, integer_types):
            return self.gs
        return self.gs_ref.nid

    def Ga(self):
        if isinstance(self.ga, integer_types) or self.ga is None:
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        if isinstance(self.gb, integer_types) or self.gb is None:
            return self.gb
        return self.gb_ref.nid

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CGAP(Element):
    type = 'CGAP'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb',
    }

    def __init__(self):
        """
        # .. todo:: not done...
        """
        Element.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.eid = integer(card, 1, 'eid')
        self.pid = integer_or_blank(card, 2, 'pid', self.eid)
        self.ga = integer_or_blank(card, 3, 'ga')
        self.gb = integer_or_blank(card, 4, 'gb')
        x1G0 = integer_double_or_blank(card, 5, 'x1_g0')
        if isinstance(x1G0, integer_types):
            self.g0 = x1G0
            self.x = None
            self.cid = None
        elif isinstance(x1G0, float):
            self.g0 = None
            x1 = x1G0
            x2 = double_or_blank(card, 6, 'x2', 0.0)
            x3 = double_or_blank(card, 7, 'x3', 0.0)
            self.x = [x1, x2, x3]
            self.cid = integer_or_blank(card, 8, 'cid', 0)
        else:
            #raise RuntimeError('invalid CGAP...x1/g0 = |%s|' %(x1G0))
            self.g0 = None
            self.x = [None, None, None]
            self.cid = None
        assert len(card) <= 9, 'len(CGAP card) = %i' % len(card)

    def add_op2_data(self, data, comment=''):
        if comment:
            self._comment = comment
        self.eid = data[0]
        self.pid = data[1]
        self.ga = data[2]
        self.gb = data[3]
        self.g0 = data[4]
        x1 = data[5]
        x2 = data[6]
        x3 = data[7]
        self.x = [x1, x2, x3]
        self.cid = data[8]

    def _verify(self, xref=True):
        cid = self.Cid()
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids

        assert cid is None or isinstance(cid, int), 'cid=%r\n%s' % (cid, str(self))
        assert isinstance(eid, int), 'eid=%r\n%s' % (eid, str(self))
        assert isinstance(pid, int), 'pid=%r\n%s' % (pid, str(self))
        for i, nid in enumerate(nids):
            assert isinstance(nid, int), 'nid%i is not an integer; nid=%s' %(i, nid)

        if xref:
            assert self.pid_ref.type in ['PGAP'], 'pid=%i self.pid_ref.type=%s' % (pid, self.pid_ref.type)
            if self.cid is not None and self.cid != 0:
                assert self.cid_ref.type in ['CORD1R', 'CORD1C', 'CORD1S', 'CORD2R', 'CORD2C',
                                         'CORD2S'], 'cid=%i self.cid.type=%s' % (cid, self.cid_ref.type)

    def cross_reference(self, model):
        msg = ' which is required by CGAP eid=%s' % self.Eid()
        self.ga = model.Node(self.Ga(), msg=msg)
        self.gb = model.Node(self.Gb(), msg=msg)
        self.ga_ref = self.ga
        self.gb_ref = self.gb
        if self.g0:
            self.g0 = model.Node(self.g0, msg=msg)
            self.x = self.g0.get_position()
            self.g0_ref = self.g0
        self.pid = model.Property(self.Pid(), msg=msg)
        self.pid_ref = self.pid
        if self.cid:
            self.cid = model.Coord(self.Cid(), msg=msg)
            self.cid_ref = self.cid

    def uncross_reference(self):
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.cid = self.Cid()
        self.pid = self.Pid()
        if self.cid:
            del self.cid_ref
        del self.ga_ref, self.gb_ref, self.pid_ref

    def Eid(self):
        return self.eid

    @property
    def nodes(self):
        return [self.ga, self.gb]

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def Cid(self):
        if isinstance(self.cid, integer_types) or self.cid is None:
            return self.cid
        return self.cid_ref.cid

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        if isinstance(self.gb, integer_types):
            return self.gb
        return self.gb_ref.nid

    def raw_fields(self):
        if self.g0 is not None:
            x = [self.g0, None, None]
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
        pass

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self. type, self.eid)
        self.nodes = model.Nodes(self.nodes, allow_empty_nodes=True, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)
        self.nodes_ref = self.nodes
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_ref


class CRAC2D(CrackElement):
    type = 'CRAC2D'
    _field_map = {
        1: 'eid', 2:'pid',
    }
    ## todo:: not done

    def __init__(self):
        CrackElement.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.eid = integer(card, 1, 'eid')
        self.pid = integer(card, 2, 'pid')
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
        assert len(card) <= 21, 'len(CRAC2D card) = %i' % len(card)
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 18

    def add_op2_data(self, data, comment=''):
        if comment:
            self._comment = comment
        self.eid = data[0]
        self.pid = data[1]
        nids = data[2:]
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 18

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids

        assert isinstance(eid, int)
        assert isinstance(pid, int)

    def Eid(self):
        return self.eid

    def get_edge_ids(self):
        return []

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

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

    def __init__(self):
        CrackElement.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.eid = integer(card, 1, 'eid')
        self.pid = integer(card, 2, 'pid')
        # required 1-10, 19-28
        # optional 11-18, 29-36, 37-64
        # all/none 37-46
        nids = fields(integer_or_blank, card, 'nid', 3, 67)  # cap at +3 = 67
        assert len(card) <= 67, 'len(CRAC3D card) = %i' % len(card)
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 64

    def add_op2_data(self, data, comment=''):
        self.eid = data[0]
        self.pid = data[1]
        nids = data[2:]
        self.prepare_node_ids(nids, allow_empty_nodes=True)
        assert len(self.nodes) == 64

    def Eid(self):
        return self.eid

    def _verify(self, xref=True):
        eid = self.Eid()
        pid = self.Pid()
        nids = self.node_ids
        assert isinstance(eid, int)
        assert isinstance(pid, int)

    def get_edge_ids(self):
        return []

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return self._nodeIDs(allow_empty_nodes=True)

    def raw_fields(self):
        list_fields = ['CRAC3D', self.eid, self.Pid()] + self.node_ids
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PLOTEL(BaseCard):
    type = 'PLOTEL'
    _field_map = {
        1: 'eid', 3:'g1', 4:'g2',
    }

    def __init__(self):
        """
        Defines a 1D dummy element used for plotting.
        +--------+-----+-----+-----+
        |   1    |  2  |  3  |  4  |
        +--------+-----+-----+-----+
        | PLOTEL | EID | G1  | G2  |
        +--------+-----+-----+-----+
        """
        BaseCard.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.eid = integer(card, 1, 'eid')
        self.nodes = [
            integer(card, 2, 'g1'),
            integer(card, 3, 'g2'),
            ]
        assert len(card) <= 4, 'len(CGAP card) = %i' % len(card)

    def add_op2_data(self, data, comment=''):
        if comment:
            self._comment = comment
        self.eid = data[0]
        self.pid = data[1]
        self.nodes = [data[2], data[3]]

    def _verify(self, xref=True):
        pass

    def cross_reference(self, model):
        msg = ' which is required by PLOTEL eid=%s' % self.Eid()
        node_ids = self.node_ids
        self.nodes = [
            model.Node(node_ids[0], msg=msg),
            model.Node(node_ids[1], msg=msg),
        ]
        self.nodes_ref = self.nodes

    def uncross_reference(self):
        self.nodes = self.node_ids
        del self.nodes_ref

    def uncross_reference(self):
        self.nodes = self.node_ids

    def Eid(self):
        return self.eid

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
        card = self.repr_fields()
        return self.comment + print_card_8(card)
