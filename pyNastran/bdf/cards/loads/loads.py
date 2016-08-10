# pylint: disable=R0902,R0904,R0914,W0231,R0201
"""
All static loads are defined in this file.  This includes:

 * LSEQ
 * DAREA
 * SLOAD
 * RFORCE
 * RANDPS
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import integer_types, PY3
from six.moves import zip, range

#from pyNastran.bdf.errors import CrossReferenceError
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard, _node_ids
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, components_or_blank,
    string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double


class Load(BaseCard):
    """defines the DefaultLoad class"""
    type = 'DefLoad'

    def __init__(self):
        self.cid = None
        self.nodes = None

    @property
    def node_ids(self):
        try:
            return self._nodeIDs()
        except:
            #raise
            raise RuntimeError('error processing nodes for \n%s' % str(self))

    def _nodeIDs(self, nodes=None):
        """returns nodeIDs for repr functions"""
        if not nodes:
            nodes = self.nodes
        if isinstance(nodes[0], integer_types):
            return [node for node in nodes]
        else:
            return [node.nid for node in nodes]


class LoadCombination(Load):  # LOAD, DLOAD
    def __init__(self, sid, scale, scale_factors, load_ids, comment=''):
        Load.__init__(self)
        if comment:
            self._comment = comment

        #: load ID
        self.sid = sid

        #: overall scale factor
        self.scale = scale

        #: individual scale factors (corresponds to load_ids)
        self.scale_factors = scale_factors

        #: individual load_ids (corresponds to scale_factors)
        self.load_ids = load_ids

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        scale = double(card, 2, 'scale')

        scale_factors = []
        load_ids = []

        # alternating of scale factor & load set ID
        nloads = len(card) - 3
        assert nloads % 2 == 0
        for i in range(nloads // 2):
            n = 2 * i + 3
            scale_factors.append(double(card, n, 'scale_factor'))
            load_ids.append(integer(card, n + 1, 'load_id'))
        return cls(sid, scale, scale_factors, load_ids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        scale = data[1]
        scale_factors = data[2]
        load_ids = data[3]
        assert len(data) == 4, '%s data=%s' % (cls.type, data)
        return cls(sid, scale, scale_factors, load_ids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        load_ids2 = []
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        for load_id in self.load_ids:
            assert load_id != self.sid, 'Type=%s sid=%s load_id=%s creates a recursion error' % (self.type, self.sid, load_id)
            load_id2 = model.Load(load_id, msg=msg)
            assert isinstance(load_id2, list), load_id2
            load_ids2.append(load_id2)
        self.load_ids = load_ids2

        self.load_ids_ref = self.load_ids

    def safe_cross_reference(self, model, debug=True):
        load_ids2 = []
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        for load_id in self.load_ids:
            try:
                load_id2 = model.Load(load_id, msg=msg)
            except KeyError:
                if debug:
                    msg = 'Couldnt find load_id=%i, which is required by %s=%s' % (
                        load_id, self.type, self.sid)
                    print(msg)
                continue
            load_ids2.append(load_id2)
        self.load_ids = load_ids2

    def LoadID(self, lid):
        if isinstance(lid, integer_types):
            return lid
        elif isinstance(lid, list):
            return lid[0].sid
        else:
            raise NotImplementedError(lid)

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        """
        .. note:: requires a cross referenced load
        """
        loads = []
        for all_loads in self.load_ids:
            assert not isinstance(all_loads, int), 'all_loads=%s\n%s' % (str(all_loads), str(self))
            for load in all_loads:
                try:
                    loads += load.get_loads()
                except RuntimeError:
                    if PY3:
                        print('recursion error on load=\n%s' % str(load))
                        raise
                    try:
                        msg = 'recursion error on load=\n%s' % str(load)
                    except RuntimeError:
                        msg = 'Recursion Error on load=%s Type=%s' % (load.sid, load.type)
                    raise RuntimeError(msg)
            #loads += self.ID  #: :: todo:  what does this mean, was uncommented
        return loads


class LSEQ(BaseCard):  # Requires LOADSET in case control deck
    """
    Defines a sequence of static load sets

    .. todo:: how does this work...
    """
    type = 'LSEQ'

    def __init__(self, sid, excite_id, lid, tid, comment=''):
        if comment:
            self._comment = comment
        self.sid = sid
        self.excite_id = excite_id
        self.lid = lid
        self.tid = tid

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        lid = integer(card, 3, 'lid')
        tid = integer_or_blank(card, 4, 'tid')
        assert len(card) <= 5, 'len(LSEQ card) = %i\ncard=%s' % (len(card), card)
        return LSEQ(sid, excite_id, lid, tid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        excite_id = data[1]
        lid = data[2]
        tid = data[3]
        return LSEQ(sid, excite_id, lid, tid, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        self.lid = model.Load(self.lid, msg=msg)
        self.lid_ref = self.lid
        #self.excite_id = model.Node(self.excite_id, msg=msg)
        if self.tid:
            self.tid = model.Table(self.tid, msg=msg)
            self.tid_ref = self.tid

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.lid = self.Lid()
        self.tid = self.Tid()

        if self.tid is not None:
            self.tid_ref
        del self.lid_ref

    def LoadID(self, lid):
        if isinstance(lid, integer_types):
            return lid
        elif isinstance(lid, list):
            return self.LoadID(lid[0])
        else:
            return lid.sid

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return self.lid

    def Lid(self):
        if isinstance(self.lid, integer_types):
            return self.lid
        else:
            return self.LoadID(self.lid)

    #@property
    #def node_id(self):
        #print('self.excite_id =', self.excite_id)
        #if isinstance(self.excite_id, integer_types):
            #return self.excite_id
        #return self.excite_id.nid

    def Tid(self):
        if self.tid is None:
            return None
        elif isinstance(self.tid, integer_types):
            return self.tid
        return self.tid_ref.tid

    def raw_fields(self):
        list_fields = ['LSEQ', self.sid, self.excite_id, self.Lid(), self.Tid()]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class LOADCYN(Load):
    type = 'LOADCYN'

    def __init__(self, sid, scale, segment_id, scales, load_ids, segment_type=None, comment=''):
        if comment:
            self._comment = comment
        self.sid = sid
        self.scale = scale
        self.segment_id = segment_id
        self.scales = scales
        self.load_ids = load_ids
        self.segment_type = segment_type

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        scale = double(card, 2, 'scale')
        segment_id = integer(card, 3, 'segment_id')
        segment_type = string_or_blank(card, 4, 'segment_type')

        scalei = double(card, 5, 'sid1')
        loadi = integer(card, 6, 'load1')
        scales = [scalei]
        load_ids = [loadi]

        scalei = double_or_blank(card, 7, 'sid2')
        if scalei is not None:
            loadi = double_or_blank(card, 8, 'load2')
            scales.append(scaali)
            load_ids.append(loadi)
        return LOADCYN(sid, scale, segment_id, scales, load_ids,
                       segment_type=segment_type, comment=comment)

    def get_loads(self):
        return [self]

    def cross_reference(self, model):
        pass

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def raw_fields(self):
        end = []
        for scale, load in zip(self.scales, self.load_ids):
            end += [scale, load]
        list_fields = ['LOADCYN', self.sid, self.scale, self.segment_id, self.segment_type
                       ] + end
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class DAREA(BaseCard):
    """
    Defines scale (area) factors for static and dynamic loads. In dynamic
    analysis, DAREA is used in conjunction with ACSRCE, RLOADi and TLOADi
    entries.

    RLOAD1 -> DAREA by SID

     +-------+-----+----+----+-----+----+----+------+
     | DAREA | SID | P1 | C1 |A1   | P2 | C2 | A2   |
     +-------+-----+----+----+-----+----+----+------+
     | DAREA | 3   | 6  | 2  | 8.2 | 15 | 1  | 10.1 |
     +-------+-----+----+----+-----+----+----+------+
    """
    type = 'DAREA'

    def __init__(self, sid, p, c, scale, comment=''):
        self.sid = sid
        self.p = p
        self.c = c
        self.scale = scale

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        noffset = 3 * icard
        sid = integer(card, 1, 'sid')
        p = integer(card, 2 + noffset, 'p')
        c = components_or_blank(card, 3 + noffset, 'c', 0)
        scale = double(card, 4 + noffset, 'scale')
        return DAREA(sid, p, c, scale, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        p = data[1]
        c = data[2]
        scale = data[3]
        assert len(data) == 4, 'data = %s' % data
        return DAREA(sid, p, c, scale, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s=%s' % (self.type, self.sid)
        self.p = model.Node(self.node_id, allow_empty_nodes=False, msg=msg)
        self.p_ref = self.p

    def safe_cross_reference(self, model, debug=True):
        msg = ', which is required by %s=%s' % (self.type, self.sid)
        self.p = model.Node(self.node_id, allow_empty_nodes=False, msg=msg)
        self.p_ref = self.p

    def uncross_reference(self):
        self.p = self.node_id
        del self.p_ref

    @property
    def node_id(self):
        if isinstance(self.p, integer_types):
            return self.p
        return self.p_ref.nid

    def raw_fields(self):
        list_fields = ['DAREA', self.sid, self.node_id, self.c, self.scale]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class TabularLoad(BaseCard):
    def __init__(self):
        pass


class SPCD(Load):
    """
    Defines an enforced displacement value for static analysis and an enforced
    motion value (displacement, velocity or acceleration) in dynamic analysis.

     +------+-----+-----+-----+------+----+----+----+
     |   1  |  2  |  3  |  4  |   5  |  6 | 7  |  8 |
     +======+=====+=====+=====+======+====+====+====+
     | SPCD | SID |  G1 | C1  |  D1  | G2 | C2 | D2 |
     +------+-----+-----+-----+------+----+----+----+
     | SPCD | 100 | 32  | 436 | -2.6 | 5  | 2  | .9 |
     +------+-----+-----+-----+------+----+----+----+
    """
    type = 'SPCD'

    def __init__(self, sid, gids, constraints, enforced, comment=''):
        if comment:
            self._comment = comment
        self.sid = sid
        self.gids = gids
        self.constraints = constraints
        self.enforced = enforced

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        if card.field(5) in [None, '']:
            gids = [integer(card, 2, 'G1'),]
            constraints = [components_or_blank(card, 3, 'C1', 0)]
            enforced = [double_or_blank(card, 4, 'D1', 0.0)]
        else:
            gids = [
                integer(card, 2, 'G1'),
                integer(card, 5, 'G2'),
            ]
            # :0 if scalar point 1-6 if grid
            constraints = [components_or_blank(card, 3, 'C1', 0),
                           components_or_blank(card, 6, 'C2', 0)]
            enforced = [double_or_blank(card, 4, 'D1', 0.0),
                        double_or_blank(card, 7, 'D2', 0.0)]
        return SPCD(sid, gids, constraints, enforced, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        gids = [data[1]]
        constraints = [data[2]]
        enforced = [data[3]]
        return SPCD(sid, gids, constraints, enforced, comment=comment)

    @property
    def node_ids(self):
        msg = ', which is required by %s=%s' % (self.type, self.sid)
        return _node_ids(self, nodes=self.gids, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by %s=%s' % (self.type, self.sid)
        self.gids = model.Nodes(self.gids, allow_empty_nodes=True, msg=msg)
        self.gids_ref = self.gids

    def safe_cross_reference(self, model, debug=True):
        msg = ', which is required by %s=%s' % (self.type, self.sid)
        self.gids = model.Nodes(self.gids, allow_empty_nodes=True, msg=msg)
        self.gids_ref = self.gids

    def uncross_reference(self):
        self.gids = self.node_ids
        del self.gids_ref

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def raw_fields(self):
        fields = ['SPCD', self.sid]
        for (gid, constraint, enforced) in zip(self.node_ids, self.constraints,
                                               self.enforced):
            fields += [gid, constraint, enforced]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class SLOAD(Load):
    """
    Static Scalar Load

    Defines concentrated static loads on scalar or grid points.

    .. note:: Can be used in statics OR dynamics.

    If Si refers to a grid point, the load is applied to component T1 of the
    displacement coordinate system (see the CD field on the GRID entry).
    """
    type = 'SLOAD'

    def __init__(self, sid, nids, mags, comment=''):
        if comment:
            self._comment = comment
        #: load ID
        self.sid = sid
        self.nids = nids
        self.mags = mags

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')

        nfields = len(card) - 2
        n = nfields // 2
        if nfields % 2 == 1:
            n += 1
            msg = 'Missing last magnitude on SLOAD card=%s' % card.fields()
            raise RuntimeError(msg)

        nids = []
        mags = []
        for i in range(n):
            j = 2 * i + 2
            nids.append(integer(card, j, 'nid' + str(i)))
            mags.append(double(card, j + 1, 'mag' + str(i)))
        return SLOAD(sid, nids, mags, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        (sid, nid, scale_factor) = data
        return SLOAD(sid, [nid], [scale_factor], comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s=%s' % (self.type, self.sid)
        for (i, nid) in enumerate(self.nids):
            self.nids[i] = model.Node(nid, msg=msg)
        self.nids_ref = self.nids

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.nids = self.node_ids
        del self.nids_ref

    def Nid(self, node):
        if isinstance(node, integer_types):
            return node
        return node.nid

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        """
        .. todo::  not done
        """
        return []

    @property
    def node_ids(self):
        return [self.Nid(nid) for nid in self.nids]

    def raw_fields(self):
        list_fields = ['SLOAD', self.sid]
        for (nid, mag) in zip(self.nids, self.mags):
            list_fields += [self.Nid(nid), mag]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class RFORCE(Load):
    type = 'RFORCE'

    def __init__(self, sid, nid, cid, scale, r1, r2, r3, method, racc,
                 mb, idrf, comment=''):
        if comment:
            self._comment = comment
        self.sid = sid
        self.nid = nid
        self.cid = cid
        self.scale = scale
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.method = method
        self.racc = racc
        self.mb = mb
        self.idrf = idrf

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        nid = integer_or_blank(card, 2, 'nid', 0)
        cid = integer_or_blank(card, 3, 'cid', 0)
        scale = double_or_blank(card, 4, 'scale', 1.)
        r1 = double_or_blank(card, 5, 'r1', 0.)
        r2 = double_or_blank(card, 6, 'r2', 0.)
        r3 = double_or_blank(card, 7, 'r3', 0.)
        method = integer_or_blank(card, 8, 'method', 1)
        racc = double_or_blank(card, 9, 'racc', 0.)
        mb = integer_or_blank(card, 10, 'mb', 0)
        idrf = integer_or_blank(card, 11, 'idrf', 0)
        assert len(card) <= 12, 'len(RFORCE card) = %i\ncard=%s' % (len(card), card)
        return RFORCE(sid, nid, cid, scale, r1, r2, r3, method, racc, mb,
                      idrf, comment=comment)

    #@classmethod
    #def add_op2_data(self, data, comment=''):
        #self.sid = data[0]
        #print("RFORCE = %s" % data)
        #raise NotImplementedError(data)
        #return RFORCE(sid, nid, cid, scale, r1, r2, r3, method, racc, mb,
                      #idrf, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by RFORCE sid=%s' % self.sid
        if self.nid > 0:
            self.nid = model.Node(self.nid, msg=msg)
            self.nid_ref = self.nid
        self.cid = model.Coord(self.cid, msg=msg)
        self.cid_ref = self.cid

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.nid = self.Nid()
        self.cid = self.Cid()
        if self.nid != 0:
            del self.nid_ref
        del self.cid_ref

    def Nid(self):
        if isinstance(self.nid, integer_types):
            return self.nid
        return self.nid_ref.nid

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = ['RFORCE', self.sid, self.Nid(), self.Cid(), self.scale,
                       self.r1, self.r2, self.r3, self.method, self.racc,
                       self.mb, self.idrf]
        return list_fields

    def repr_fields(self):
        #method = set_blank_if_default(self.method,1)
        racc = set_blank_if_default(self.racc, 0.)
        mb = set_blank_if_default(self.mb, 0)
        idrf = set_blank_if_default(self.idrf, 0)
        list_fields = ['RFORCE', self.sid, self.Nid(), self.Cid(), self.scale,
                       self.r1, self.r2, self.r3, self.method, racc,
                       mb, idrf]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class RFORCE1(Load):
    type = 'RFORCE1'

    def __init__(self, sid, nid, scale, r1, r2, r3, racc,
                 mb, group_id, cid=0, method=2, comment=''):
        if comment:
            self._comment = comment
        self.sid = sid
        self.nid = nid
        self.cid = cid
        self.scale = scale
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.method = method
        self.racc = racc
        self.mb = mb
        self.group_id = group_id

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        nid = integer_or_blank(card, 2, 'nid', 0)
        cid = integer_or_blank(card, 3, 'cid', 0)
        scale = double_or_blank(card, 4, 'scale', 1.)
        r1 = double_or_blank(card, 5, 'r1', 0.)
        r2 = double_or_blank(card, 6, 'r2', 0.)
        r3 = double_or_blank(card, 7, 'r3', 0.)
        method = integer_or_blank(card, 8, 'method', 1)
        racc = double_or_blank(card, 9, 'racc', 0.)
        mb = integer_or_blank(card, 10, 'mb', 0)
        group_id = integer_or_blank(card, 11, 'group_id', 0)
        assert len(card) <= 12, 'len(RFORCE card) = %i\ncard=%s' % (len(card), card)
        return RFORCE1(sid, nid, scale, r1, r2, r3, racc, mb,
                       group_id, cid=cid, method=method, comment=comment)

    def get_loads(self):
        return [self]

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by RFORCE1 sid=%s' % self.sid
        if self.nid > 0:
            self.nid = model.Node(self.nid, msg=msg)
            self.nid_ref = self.nid
        self.cid = model.Coord(self.cid, msg=msg)
        self.cid_ref = self.cid

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.nid = self.Nid()
        self.cid = self.Cid()
        if self.nid != 0:
            del self.nid_ref
        del self.cid_ref

    def Nid(self):
        if isinstance(self.nid, integer_types):
            return self.nid
        return self.nid_ref.nid

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def raw_fields(self):
        list_fields = ['RFORCE1', self.sid, self.Nid(), self.Cid(), self.scale,
                       self.r1, self.r2, self.r3, self.method, self.racc,
                       self.mb, self.group_id]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class RandomLoad(BaseCard):
    def __init__(self, card, data):
        pass


class RANDPS(RandomLoad):
    r"""
    Power Spectral Density Specification

    Defines load set power spectral density factors for use in random analysis
    having the frequency dependent form:

    .. math:: S_{jk}(F) = (X+iY)G(F)
    """
    type = 'RANDPS'

    def __init__(self, sid, j, k, x=0., y=0., tid=0, comment=''):
        if comment:
            self._comment = comment
        #: Random analysis set identification number. (Integer > 0)
        #: Defined by RANDOM in the Case Control Deck.
        self.sid = sid

        #: Subcase identification number of the excited load set.
        #: (Integer > 0)
        self.j = j

        #: Subcase identification number of the applied load set.
        #: (Integer >= 0; K >= J)
        self.k = k

        #: Components of the complex number. (Real)
        self.x = x
        self.y = y

        #: Identification number of a TABRNDi entry that defines G(F).
        self.tid = tid

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        j = integer(card, 2, 'j')
        k = integer(card, 3, 'k')
        x = double_or_blank(card, 4, 'x', 0.0)
        y = double_or_blank(card, 5, 'y', 0.0)
        tid = integer_or_blank(card, 6, 'tid', 0)
        assert len(card) <= 7, 'len(RANDPS card) = %i\ncard=%s' % (len(card), card)
        return RANDPS(sid, j, k, x, y, tid, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        if self.tid:
            msg = ' which is required by RANDPS sid=%s' % (self.sid)
            #self.tid = model.Table(self.tid, msg=msg)
            self.tid = model.RandomTable(self.tid, msg=msg)
            self.tid_ref = self.tid

    def safe_cross_reference(self, model):
        return self.cross_reference(model)

    def uncross_reference(self):
        self.tid = self.Tid()
        if self.tid is not None:
            del self.tid_ref

    def getLoads(self):
        self.deprecated('getLoads()', 'get_loads()', '0.8')
        return self.get_loads()

    def get_loads(self):
        return [self]

    def Tid(self):
        if self.tid == 0:
            return None
        elif isinstance(self.tid, integer_types):
            return self.tid
        return self.tid_ref.tid

    def raw_fields(self):
        list_fields = ['RANDPS', self.sid, self.j, self.k, self.x, self.y,
                       self.Tid()]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
