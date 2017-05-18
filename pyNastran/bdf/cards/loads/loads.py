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
import numpy as np

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
    """Common method for LOAD, DLOAD"""
    def __init__(self, sid, scale, scale_factors, load_ids, comment=''):
        """
        Common method for LOAD, DLOAD

        Parameters
        ----------
        sid : int
            load id
        scale : float
            overall scale factor
        scale_factors : List[float]
            individual scale factors (corresponds to load_ids)
        load_ids : List[int]
            individual load_ids (corresponds to scale_factors)
        comment : str; default=''
            a comment for the card
        """
        Load.__init__(self)
        if comment:
            self.comment = comment

        #: load ID
        self.sid = sid

        #: overall scale factor
        self.scale = scale

        #: individual scale factors (corresponds to load_ids)
        if isinstance(scale_factors, float):
            scale_factors = [scale_factors]
        self.scale_factors = scale_factors

        #: individual load_ids (corresponds to scale_factors)
        if isinstance(load_ids, int):
            load_ids = [load_ids]
        self.load_ids = load_ids
        assert 0 not in load_ids, self

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        scale = double(card, 2, 'scale')

        scale_factors = []
        load_ids = []

        # alternating of scale factor & load set ID
        nloads = len(card) - 3
        assert nloads % 2 == 0, 'card=%s' % card
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
        self.load_ids_ref = self.load_ids

    def LoadID(self, lid):
        if isinstance(lid, integer_types):
            return lid
        elif isinstance(lid, list):
            return lid[0].sid
        else:
            raise NotImplementedError(lid)

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
    +------+-----+----------+-----+-----+
    |   1  |  2  |     3    |  4  |  5  |
    +======+=====+==========+=====+=====+
    | LSEQ | SID | EXCITEID | LID | TID |
    +------+-----+----------+-----+-----+

    ACSRCE : If there is no LOADSET Case Control command, then EXCITEID
             may reference DAREA and SLOAD entries. If there is a LOADSET
             Case Control command, then EXCITEID may reference DAREA
             entries as well as SLOAD entries specified by the LID field
             in the selected LSEQ entry corresponding to EXCITEID.

    DAREA :  Refer to RLOAD1, RLOAD2, TLOAD1, TLOAD2, or ACSRCE entries
             for the formulas that define the scale factor Ai in dynamic
             analysis.

    DPHASE :

    SLOAD :  In the static solution sequences, the load set ID (SID) is
             selected by the Case Control command LOAD. In the dynamic
             solution sequences, SID must be referenced in the LID field
             of an LSEQ entry, which in turn must be selected by the Case
             Control command LOADSET.

    LSEQ LID : Load set identification number of a set of static load
               entries such as those referenced by the LOAD Case Control
               command.


    LSEQ,  SID, EXCITEID, LID, TID

    #--------------------------------------------------------------
    # F:\\Program Files\\Siemens\\NXNastran\\nxn10p1\\nxn10p1\\nast\\tpl\\cube_iter.dat

    DLOAD       1001     1.0     1.0   55212
    sid = 1001
    load_id = [55212] -> RLOAD2.SID

    RLOAD2,     SID, EXCITEID, DELAYID, DPHASEID,   TB,     TP,  TYPE
    RLOAD2     55212   55120              55122   55123   55124
    EXCITEID = 55120 -> DAREA.SID
    DPHASEID = 55122 -> DPHASE.SID

    DARA        SID      NID    COMP  SCALE
    DAREA      55120     913    3     9.9E+9
    SID = 55120 -> RLOAD2.SID

    DPHASE      SID     POINTID   C1    TH1
    DPHASE     55122     913       3   -90.0
    SID = 55122
    POINTID = 913 -> GRID.NID

    GRID       NID       X     Y     Z
    GRID       913      50.  0.19  -39.9
    """
    type = 'LSEQ'

    def __init__(self, sid, excite_id, lid, tid=None, comment=''):
        """
        Creates a LSEQ card

        Parameters
        ----------
        sid : int
            loadset id; LOADSET points to this
        excite_id : int
            set id assigned to this static load vector
        lid : int
            load set id of a set of static load entries;
            LOAD in the Case Control
        tid : int; default=None
            temperature set id of a set of thermal load entries;
            TEMP(LOAD) in the Case Control
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.sid = sid
        self.excite_id = excite_id
        self.lid = lid
        self.tid = tid

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a LSEQ card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        lid = integer(card, 3, 'lid')
        tid = integer_or_blank(card, 4, 'tid')
        assert len(card) <= 5, 'len(LSEQ card) = %i\ncard=%s' % (len(card), card)
        return LSEQ(sid, excite_id, lid, tid=None, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds an LSEQ card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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
        msg = ', which is required by LSEQ=%s' % (self.sid)
        self.lid = model.Load(self.lid, msg=msg)
        self.lid_ref = self.lid
        #self.excite_id = model.Node(self.excite_id, msg=msg)
        if self.tid:
            # TODO: temperature set, not a table?
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
            self.comment = comment
        self.sid = sid
        self.scale = scale
        self.segment_id = segment_id
        self.scales = scales
        self.load_ids = load_ids
        self.segment_type = segment_type

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a LOADCYN card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
    |   1   |  2  | 3  |  4 |  5  | 6  |  7 |  8   |
    +=======+=====+====+====+=====+====+====+======+
    | DAREA | SID | P1 | C1 |  A1 | P2 | C2 |  A2  |
    +-------+-----+----+----+-----+----+----+------+
    | DAREA |  3  | 6  | 2  | 8.2 | 15 | 1  | 10.1 |
    +-------+-----+----+----+-----+----+----+------+
    """
    type = 'DAREA'

    def __init__(self, sid, nodes, components, scales, comment=''):
        """
        Creates a DAREA card

        Parameters
        ----------
        sid : int
            darea id
        nodes : List[int]
            GRID, EPOINT, SPOINT id
        components : List[int]
            Component number. (0-6; 0-EPOINT/SPOINT; 1-6 GRID)
        scales : List[float]
            Scale (area) factor
        """
        self.sid = sid
        if isinstance(nodes, integer_types):
            nodes = [nodes]
        if isinstance(components, integer_types):
            components = [components]
        if isinstance(scales, float):
            scales = [scales]

        self.nodes = nodes
        self.components = components

        assert isinstance(components, list), 'components=%r' % components
        for component in components:
            assert 0 <= component <= 6, component
        self.scales = scales
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        noffset = 3 * icard
        sid = integer(card, 1, 'sid')
        nid = integer(card, 2 + noffset, 'p')
        component = int(components_or_blank(card, 3 + noffset, 'c', 0))
        scale = double(card, 4 + noffset, 'scale')
        return DAREA(sid, nid, component, scale, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a DAREA card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid = data[0]
        p = data[1]
        c = data[2]
        scale = data[3]
        assert len(data) == 4, 'data = %s' % data
        return DAREA(sid, p, c, scale, comment=comment)

    def add(self, darea):
        assert self.sid == darea.sid, 'sid=%s darea.sid=%s' % (self.sid, darea.sid)
        if darea.comment:
            if hasattr('_comment'):
                self._comment += darea.comment
            else:
                self._comment = darea.comment
        self.nodes += darea.nodes
        self.components += darea.components
        self.scales += darea.scales

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by DAREA=%s' % (self.sid)
        self.nodes_ref = model.Nodes(self.node_ids, allow_empty_nodes=False, msg=msg)

    def safe_cross_reference(self, model, debug=True):
        nids2 = []
        msg = ', which is required by DAREA=%s' % (self.sid)
        for nid in self.node_ids:
            try:
                nid2 = model.Node(nid, msg=msg)
            except KeyError:
                if debug:
                    msg = 'Couldnt find nid=%i, which is required by DAREA=%s' % (
                        nid, self.sid)
                    print(msg)
                continue
            nids2.append(nid2)
        self.nodes_ref = nids2

    def uncross_reference(self):
        self.nodes_ref = None

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        msg = ', which is required by DAREA=%s' % (self.sid)
        return _node_ids(self, nodes=self.nodes_ref, allow_empty_nodes=False, msg=msg)

    def raw_fields(self):
        for nid, comp, scale in zip(self.node_ids, self.components, self.scales):
            list_fields = ['DAREA', self.sid, nid, comp, scale]
        return list_fields

    def write_card(self, size=8, is_double=False):
        msg = self.comment
        for nid, comp, scale in zip(self.node_ids, self.components, self.scales):
            msg += print_card_8(['DAREA', self.sid, nid, comp, scale])
        return msg


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
            self.comment = comment
        self.sid = sid
        self.gids = gids
        self.constraints = constraints
        self.enforced = enforced

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPCD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        """
        Adds an SPCD card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid = data[0]
        gids = [data[1]]
        constraints = [data[2]]
        enforced = [data[3]]
        return SPCD(sid, gids, constraints, enforced, comment=comment)

    @property
    def node_ids(self):
        msg = ', which is required by SPCD=%s' % (self.sid)
        return _node_ids(self, nodes=self.gids, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by SPCD=%s' % (self.sid)
        self.gids = model.Nodes(self.gids, allow_empty_nodes=True, msg=msg)
        self.gids_ref = self.gids

    def safe_cross_reference(self, model, debug=True):
        msg = ', which is required by SPCD=%s' % (self.sid)
        self.gids = model.Nodes(self.gids, allow_empty_nodes=True, msg=msg)
        self.gids_ref = self.gids

    def uncross_reference(self):
        self.gids = self.node_ids
        del self.gids_ref

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
        """
        Creates an SLOAD (SPOINT load)

        Parameters
        ----------
        sid : int
            load id
        nids : int; List[int]
            the SPOINT ids
        mags : float; List[float]
            the SPOINT loads
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        if isinstance(nids, integer_types):
            nids = [nids]
        if isinstance(mags, float):
            mags = [mags]
        #: load ID
        self.sid = sid
        self.nids = nids
        self.mags = mags

    def validate(self):
        assert len(self.nids) == len(self.mags), 'len(nids)=%s len(mags)=%s' % (len(self.nids), len(self.mags))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SLOAD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
        """
        Adds an SLOAD card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
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
        msg = ' which is required by SLOAD=%s' % (self.sid)
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

    def __init__(self, sid, nid, cid, scale, r1, r2, r3, method=1, racc=0.,
                 mb=0, idrf=0, comment=''):
        """
        idrf doesn't exist in MSC 2005r2; exists in MSC 2016
        """
        if comment:
            self.comment = comment
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

    def validate(self):
        assert self.method in [1, 2], self.method

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RFORCE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a RFORCE card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid, nid, cid, a, r1, r2, r3, method, racc, mb = data
        scale = 1.0
        return RFORCE(sid, nid, cid, scale, r1, r2, r3, method=method, racc=racc, mb=mb,
                      idrf=0, comment=comment)

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

    @property
    def node_id(self):
        if isinstance(self.nid, integer_types):
            return self.nid
        return self.nid_ref.nid

    def Nid(self):
        return self.node_id

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = ['RFORCE', self.sid, self.node_id, self.Cid(), self.scale,
                       self.r1, self.r2, self.r3, self.method, self.racc,
                       self.mb, self.idrf]
        return list_fields

    def repr_fields(self):
        #method = set_blank_if_default(self.method,1)
        racc = set_blank_if_default(self.racc, 0.)
        mb = set_blank_if_default(self.mb, 0)
        idrf = set_blank_if_default(self.idrf, 0)
        list_fields = ['RFORCE', self.sid, self.node_id, self.Cid(), self.scale,
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
    """
    NX Nastran specific card

    +---------+------+----+---------+---+----+----+----+--------+
    |    1    |   2  | 3  |    4    | 5 |  6 |  7 |  8 |   9    |
    +=========+======+====+=========+===+====+====+====+========+
    | RFORCE1 | SID  | G  |   CID   | A | R1 | R2 | R3 | METHOD |
    +---------+------+----+---------+---+----+----+----+--------+
    |         | RACC | MB | GROUPID |   |    |    |    |        |
    +---------+------+----+---------+---+----+----+----+--------+
    """
    type = 'RFORCE1'

    def __init__(self, sid, nid, scale, group_id,
                 cid=0, r123=None, racc=0., mb=0, method=2, comment=''):
        """
        Defines the RFORCE1 card

        Parameters
        ----------
        sid : int
            load set id
        nid : int
            grid point through which the rotation vector acts
        scale : float
            scale factor of the angular velocity in revolutions/time
        r123 : List[float, float, float] / (3, ) float ndarray
            rectangular components of the rotation vector R that passes
            through point G
        racc : int; default=0.0
        mb : int; default=0
            Indicates whether the CID coordinate system is defined in the main
            Bulk Data Section (MB = -1) or the partitioned superelement Bulk
            Data Section (MB = 0). Coordinate systems referenced in the main
            Bulk Data Section are considered stationary with respect to the
            assembly basic coordinate system.
        group_id : int
            Group identification number. The GROUP entry referenced in the
            GROUPID field selects the grid points to which the load is applied.
        cid : int; default=0
            Coordinate system defining the components of the rotation vector.
        method : int; default=2
            Method used to compute centrifugal forces due to angular velocity.
        comment : str
        """
        if comment:
            self.comment = comment
        self.sid = sid
        self.nid = nid
        self.cid = cid
        self.scale = scale
        if r123 is None:
            self.r123 = np.array([1., 0., 0.])
        else:
            self.r123 = np.asarray(r123)
        self.method = method
        self.racc = racc
        self.mb = mb
        self.group_id = group_id

    def validate(self):
        if not np.linalg.norm(self.r123) > 0.:
            msg = 'r123=%s norm=%s' % (self.r123, np.linalg.norm(self.r123))
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RFORCE1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        nid = integer_or_blank(card, 2, 'nid', 0)
        cid = integer_or_blank(card, 3, 'cid', 0)
        scale = double_or_blank(card, 4, 'scale', 1.)
        r123 = [
            double_or_blank(card, 5, 'r1', 1.),
            double_or_blank(card, 6, 'r2', 0.),
            double_or_blank(card, 7, 'r3', 0.),
        ]
        method = integer_or_blank(card, 8, 'method', 1)
        racc = double_or_blank(card, 9, 'racc', 0.)
        mb = integer_or_blank(card, 10, 'mb', 0)
        group_id = integer_or_blank(card, 11, 'group_id', 0)
        assert len(card) <= 12, 'len(RFORCE1 card) = %i\ncard=%s' % (len(card), card)
        return RFORCE1(sid, nid, scale, cid=cid, r123=r123, racc=racc,
                       mb=mb, group_id=group_id, method=method, comment=comment)

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
        self.nid = self.node_id
        self.cid = self.Cid()
        if self.nid != 0:
            del self.nid_ref
        del self.cid_ref

    @property
    def node_id(self):
        if isinstance(self.nid, integer_types):
            return self.nid
        return self.nid_ref.nid

    def Nid(self):
        return self.node_id

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def raw_fields(self):
        list_fields = (['RFORCE1', self.sid, self.node_id, self.Cid(), self.scale]
                       + list(self.r123) + [self.method, self.racc,
                                            self.mb, self.group_id])
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
        """
        Creates a RANDPS card

        Parameters
        ----------
        sid : int
            random analysis set id
            defined by RANDOM in the case control deck
        j : int
            Subcase id of the excited load set
        k : int
            Subcase id of the applied load set
            k > j
        x / y : float; default=0.0
            Components of the complex number
        tid : int; default=0
            TABRNDi id that defines G(F)
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
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
        assert self.sid > 0, 'sid=%s\n%s' % (self.sid, self)

    def validate(self):
        assert self.k >= self.j, 'k=%s j=%s\n%s' % (self.k, self.j, self)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RANDPS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
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
