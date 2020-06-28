# pylint: disable=R0902,R0904,R0914,W0231,R0201
"""
All static loads are defined in this file.  This includes:

 * LSEQ
 * DAREA
 * SLOAD
 * RFORCE
 * RANDPS

"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

#from pyNastran.bdf.errors import CrossReferenceError
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard, _node_ids
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, components_or_blank,
    string, string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.utils.numpy_utils import integer_types, float_types
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class Load(BaseCard):
    """defines the DefaultLoad class"""
    type = 'DefLoad'

    def __init__(self):
        self.cid = None
        self.nodes = None

    @property
    def node_ids(self):
        """get the node ids"""
        try:
            return self._node_ids()
        except:
            #raise
            raise RuntimeError('error processing nodes for \n%s' % str(self))

    def _node_ids(self, nodes=None):
        """returns node ids for repr functions"""
        if not nodes:
            nodes = self.nodes
        if isinstance(nodes[0], integer_types):
            return [node for node in nodes]
        else:
            return [node.nid for node in nodes]


class LoadCombination(BaseCard):
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
        BaseCard.__init__(self)
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
        self.load_ids_ref = None

    def validate(self):
        msg = ''
        if not isinstance(self.scale, float_types):
            msg += 'scale=%s must be a float; type=%s\n' % (self.scale, type(self.scale))
        assert isinstance(self.scale_factors, list), self.scale_factors
        assert isinstance(self.load_ids, list), self.load_ids
        if len(self.scale_factors) != len(self.load_ids):
            msg += 'scale_factors=%s load_ids=%s\n' % (self.scale_factors, self.load_ids)
        if msg:
            raise IndexError(msg)
        for scalei, load_id in zip(self.scale_factors, self.get_load_ids()):
            assert isinstance(scalei, float_types), scalei
            assert isinstance(load_id, integer_types), load_id

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

        assert len(card) > 3, 'len(%s card) = %i\ncard=%s' % (cls.__name__, len(card), card)
        return cls(sid, scale, scale_factors, load_ids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        scale = data[1]
        scale_factors = data[2]
        load_ids = data[3]
        assert len(data) == 4, '%s data=%s' % (cls.type, data)
        return cls(sid, scale, scale_factors, load_ids, comment=comment)

    def LoadID(self, lid):
        if isinstance(lid, integer_types):
            return lid
        elif isinstance(lid, list):
            return lid[0].sid
        else:
            raise NotImplementedError(lid)

    def get_load_ids(self):
        """
        xref/non-xref way to get the load ids
        """
        if self.load_ids_ref is None:
            return self.load_ids
        load_ids = []
        supported_loads = [
            'FORCE', 'FORCE1', 'FORCE2', 'MOMENT', 'MOMENT1', 'MOMENT2',
            'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4', 'GRAV', 'SPCD', 'GMLOAD',
            'RLOAD1', 'RLOAD2', 'TLOAD1', 'TLOAD2', 'PLOADX1', 'LOAD',
            'RFORCE', 'RFORCE1', #'RFORCE2'
            'ACCEL', 'ACCEL1', 'SLOAD', 'ACSRCE',
        ]
        for loads in self.load_ids_ref:
            load_idsi = []
            for load in loads:
                if isinstance(load, integer_types):
                    load_ids.append(load)
                #elif load.type == 'LOAD':
                    #load_ids.append(load.sid)
                elif load.type in supported_loads:
                    load_idsi.append(load.sid)
                else:
                    msg = ('The get_load_ids method doesnt support %s cards.\n'
                           '%s' % (load.__class__.__name__, str(load)))
                    raise NotImplementedError(msg)

                load_idi = list(set(load_idsi))
                assert len(load_idi) == 1, load_idsi
            load_ids.append(load_idi[0])
        return load_ids

    def get_loads(self):
        """
        .. note:: requires a cross referenced load
        """
        loads = []
        for all_loads in self.load_ids_ref:
            assert not isinstance(all_loads, int), 'all_loads=%s\n%s' % (str(all_loads), str(self))
            for load in all_loads:
                try:
                    loads += load.get_loads()
                except RuntimeError:
                    print('recursion error on load=\n%s' % str(load))
                    raise
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

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        excite_id = 2
        lid = 3
        return LSEQ(sid, excite_id, lid, tid=None, comment='')

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
        self.lid_ref = None
        self.tid_ref = None

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
        load_id = integer_or_blank(card, 3, 'lid')
        temp_id = integer_or_blank(card, 4, 'tid')
        if load_id is None and temp_id is None:
            msg = 'LSEQ load_id/temp_id must not be None; load_id=%s temp_id=%s' % (load_id, temp_id)
            raise RuntimeError(msg)
        assert len(card) <= 5, 'len(LSEQ card) = %i\ncard=%s' % (len(card), card)
        return LSEQ(sid, excite_id, load_id, tid=temp_id, comment=comment)

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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by LSEQ=%s' % (self.sid)
        if self.lid is not None:
            self.lid_ref = model.Load(self.lid, consider_load_combinations=True, msg=msg)
        #self.excite_id = model.Node(self.excite_id, msg=msg)
        if self.tid:
            # TODO: temperature set, not a table?
            self.tid_ref = model.Load(self.tid, consider_load_combinations=True, msg=msg)
            #self.tid_ref = model.Table(self.tid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        return self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.lid = self.Lid()
        self.tid = self.Tid()
        self.lid_ref = None
        self.tid_ref = None

    def LoadID(self, lid: int) -> int:
        if isinstance(lid, list):
            sid = self.LoadID(lid[0])
        elif isinstance(lid, integer_types):
            sid = lid
        else:
            sid = lid.sid
        return sid

    def get_loads(self) -> Any:
        return self.lid_ref

    def Lid(self) -> int:
        if self.lid_ref is not None:
            return self.LoadID(self.lid_ref)
        return self.lid

    #@property
    #def node_id(self):
        #print('self.excite_id =', self.excite_id)
        #if isinstance(self.excite_id, integer_types):
            #return self.excite_id
        #return self.excite_id.nid

    #def Tid(self):
        #if self.tid_ref is not None:
            #return self.tid_ref.tid
        #return self.tid

    def Tid(self):
        if self.tid_ref is not None:
            return self.LoadID(self.tid_ref)
        return self.tid

    def raw_fields(self):
        list_fields = ['LSEQ', self.sid, self.excite_id, self.Lid(), self.Tid()]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class LOADCYN(Load):
    type = 'LOADCYN'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        scale = 1.
        segment_id = 2
        scales = [1., 2.]
        load_ids = [10, 20]
        return LOADCYN(sid, scale, segment_id, scales, load_ids, segment_type=None, comment='')

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

        scalei = double(card, 5, 'scale1')
        loadi = integer(card, 6, 'load1')
        scales = [scalei]
        load_ids = [loadi]

        scalei = double_or_blank(card, 7, 'scale2')
        if scalei is not None:
            loadi = double_or_blank(card, 8, 'load2')
            scales.append(scalei)
            load_ids.append(loadi)
        return LOADCYN(sid, scale, segment_id, scales, load_ids,
                       segment_type=segment_type, comment=comment)

    def get_loads(self):
        return [self]

    def cross_reference(self, model: BDF) -> None:
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def safe_cross_reference(self, model, xref_errors):
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class LOADCYH(BaseCard):
    """
    Harmonic Load Input for Cyclic Symmetry
    Defines the harmonic coefficients of a static or dynamic load for
    use in cyclic symmetry analysis.

    +---------+-----+----------+-----+-------+-----+------+------+------+
    |    1    |  2  |     3    |  4  |   5   |  6  |   7  |   8  |   9  |
    +=========+=====+==========+=====+=======+=====+======+======+======+
    | LOADCYH | SID |     S    | HID | HTYPE | S1  |  L1  |  S2  |  L2  |
    +---------+-----+----------+-----+-------+-----+------+------+------+
    """
    type = 'LOADCYH'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        scale = 1.
        hid = 0
        htype = 'C'
        scales = [1.]
        load_ids = [2]
        return LOADCYH(sid, scale, hid, htype, scales, load_ids, comment='')

    def __init__(self, sid, scale, hid, htype, scales, load_ids, comment=''):
        """
        Creates a LOADCYH card

        Parameters
        ----------
        sid : int
            loadset id; LOADSET points to this
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        self.sid = sid
        self.scale = scale
        self.hid = hid
        self.htype = htype
        self.scales = scales
        self.load_ids = load_ids
        assert htype in {'C', 'S', 'CSTAR', 'SSTAR', 'GRAV', 'RFORCE', None}, htype

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a LOADCYH card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        scale = double(card, 2, 's')
        hid = integer(card, 3, 'hid')
        htype = string_or_blank(card, 4, 'htype')
        scale1 = double(card, 5, 'scale1')
        load1 = integer_or_blank(card, 6, 'load1')
        scale2 = double_or_blank(card, 7, 'scale2')
        load2 = integer_or_blank(card, 8, 'load2')
        scales = []
        load_ids = []
        if load1 != 0:
            load_ids.append(load1)
            scales.append(scale1)
        if load2 != 0:
            load_ids.append(load2)
            scales.append(scale2)
        assert len(card) <= 7, 'len(LOADCYH card) = %i\ncard=%s' % (len(card), card)
        return LOADCYH(sid, scale, hid, htype, scales, load_ids, comment=comment)

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
        raise NotImplementedError()
        #sid = data[0]
        #excite_id = data[1]
        #lid = data[2]
        #tid = data[3]
        #return LSEQ(sid, excite_id, lid, tid, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by LOADCYH=%s' % (self.sid)

    def safe_cross_reference(self, model, xref_errors):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass
    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = ['LOADCYH', self.sid, self.scale, self.hid, self.htype]
        for scale, load_id in zip(self.scales, self.load_ids):
            list_fields += [scale, load_id]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        nodes = [1]
        components = [1]
        scales = [1.]
        return DAREA(sid, nodes, components, scales, comment='')

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
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
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
            assert 0 <= component <= 6, 'component=%r' % component
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
            if hasattr(self, '_comment'):
                self._comment += darea.comment
            else:
                self._comment = darea.comment
        self.nodes += darea.nodes
        self.components += darea.components
        self.scales += darea.scales

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by DAREA=%s' % (self.sid)
        self.nodes_ref = model.Nodes(self.node_ids, msg=msg)

    def safe_cross_reference(self, model, xref_errors, debug=True):
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

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        msg = self.comment
        for nid, comp, scale in zip(self.node_ids, self.components, self.scales):
            msg += print_card_8(['DAREA', self.sid, nid, comp, scale])
        return msg


class DynamicLoad(BaseCard):
    def __init__(self):
        pass


class SPCD(Load):
    """
    Defines an enforced displacement value for static analysis and an
    enforced motion value (displacement, velocity or acceleration) in
    dynamic analysis.

     +------+-----+-----+-----+------+----+----+----+
     |   1  |  2  |  3  |  4  |   5  |  6 | 7  |  8 |
     +======+=====+=====+=====+======+====+====+====+
     | SPCD | SID |  G1 | C1  |  D1  | G2 | C2 | D2 |
     +------+-----+-----+-----+------+----+----+----+
     | SPCD | 100 | 32  | 436 | -2.6 | 5  | 2  | .9 |
     +------+-----+-----+-----+------+----+----+----+
    """
    type = 'SPCD'
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        nodes = [1]
        components = ['1']
        enforced = 1.
        return SPCD(sid, nodes, components, enforced, comment='')

    def __init__(self, sid, nodes, components, enforced, comment=''):
        """
        Creates an SPCD card, which defines the degree of freedoms to be
        set during enforced motion

        Parameters
        ----------
        conid : int
            constraint id
        nodes : List[int]
            GRID/SPOINT ids
        components : List[str]
            the degree of freedoms to constrain (e.g., '1', '123')
        enforced : List[float]
            the constrained value for the given node (typically 0.0)
        comment : str; default=''
            a comment for the card

        .. note:: len(nodes) == len(components) == len(enforced)
        .. warning:: Non-zero enforced deflection requires an SPC/SPC1 as well.
                     Yes, you really want to constrain the deflection to 0.0
                     with an SPC1 card and then reset the deflection using an
                     SPCD card.
        """
        if comment:
            self.comment = comment
        self.sid = sid
        if isinstance(nodes, int):
            nodes = [nodes]
        if isinstance(components, str):
            components = [components]
        elif isinstance(components, int):
            components = [str(components)]
        if isinstance(enforced, float):
            enforced = [enforced]
        self.nodes = nodes
        self.components = components
        self.enforced = enforced
        self.nodes_ref = None
        assert isinstance(self.nodes, list), self.nodes
        assert isinstance(self.components, list), self.components
        assert isinstance(self.enforced, list), self.enforced


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
            nodes = [integer(card, 2, 'G1'),]
            components = [components_or_blank(card, 3, 'C1', 0)]
            enforced = [double_or_blank(card, 4, 'D1', 0.0)]
        else:
            nodes = [
                integer(card, 2, 'G1'),
                integer(card, 5, 'G2'),
            ]
            # :0 if scalar point 1-6 if grid
            components = [components_or_blank(card, 3, 'C1', 0),
                          components_or_blank(card, 6, 'C2', 0)]
            enforced = [double_or_blank(card, 4, 'D1', 0.0),
                        double_or_blank(card, 7, 'D2', 0.0)]
        return SPCD(sid, nodes, components, enforced, comment=comment)

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
        nodes = [data[1]]
        components = [data[2]]
        enforced = [data[3]]
        return SPCD(sid, nodes, components, enforced, comment=comment)

    @property
    def constraints(self):
        self.deprecated('constraints', 'components', '1.2')
        return self.components

    @constraints.setter
    def constraints(self, constraints):
        self.deprecated('constraints', 'components', '1.2')
        self.components = constraints

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        msg = ', which is required by SPCD=%s' % (self.sid)
        return _node_ids(self, nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SPCD=%s' % (self.sid)
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model, xref_errors, debug=True):
        msg = ', which is required by SPCD=%s' % (self.sid)
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def get_loads(self):
        return [self]

    def raw_fields(self):
        fields = ['SPCD', self.sid]
        for (nid, component, enforced) in zip(self.node_ids, self.components,
                                               self.enforced):
            fields += [nid, component, enforced]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        elif is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)

class DEFORM(Load):
    """
    Defines an enforced displacement value for static analysis.

     +--------+-----+-----+------+----+----+----+----+
     |    1   |  2  |  3  |   5  |  6 |  8 |  6 |  8 |
     +========+=====+=====+======+====+====+====+====+
     | DEFORM | SID |  E1 |  D1  | E2 | D2 | E3 | D3 |
     +--------+-----+-----+------+----+----+----+----+
     | DEFORM | 100 | 32  | -2.6 | 5  | .9 | 6  | .9 |
     +--------+-----+-----+------+----+----+----+----+
    """
    type = 'DEFORM'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        eid = 1
        deformation = 1.
        return DEFORM(sid, eid, deformation, comment='')

    def __init__(self, sid, eid, deformation, comment=''):
        """
        Creates an DEFORM card, which defines applied deformation on
        a 1D elemment.  Links to the DEFORM card in the case control
        deck.

        Parameters
        ----------
        sid : int
            load id
        eid : int
            CTUBE/CROD/CONROD/CBAR/CBEAM element id
        deformation : float
            the applied deformation
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        self.sid = sid
        self.eid = eid
        self.deformation = deformation
        self.eid_ref = None

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        """
        Adds a DEFORM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        offset = 2 * icard
        sid = integer(card, 1, 'sid')
        eid = integer(card, 2 + offset, 'eid%i' % (icard + 1))
        deformation = double(card, 3 + offset, 'D%i' % (icard + 1))
        return DEFORM(sid, eid, deformation, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds an DEFORM card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        sid = data[0]
        eid = data[1]
        deformation = data[2]
        return DEFORM(sid, eid, deformation, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by DEFORM=%s' % (self.sid)
        self.eid_ref = model.Element(self.eid, msg)

    def safe_cross_reference(self, model, xref_errors, debug=True):
        msg = ', which is required by DEFORM=%s' % (self.sid)
        self.eid_ref = model.safe_element(self.eid, self.sid, xref_errors, msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.eid = self.Eid()
        self.eid_ref = None

    def get_loads(self):
        return [self]

    def Eid(self):
        if self.eid_ref is None:
            return self.eid
        return self.eid_ref.eid

    def raw_fields(self):
        fields = ['DEFORM', self.sid, self.Eid(), self.deformation]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class SLOAD(Load):
    """
    Static Scalar Load
    Defines concentrated static loads on scalar or grid points.

    +-------+-----+----+-----+----+------+----+-------+
    |   1   |  2  | 3  |  4  |  5 |  6   |  7 |   8   |
    +=======+=====+====+=====+====+======+====+=======+
    | SLOAD | SID | S1 | F1  | S2 |  F2  | S3 |   F3  |
    +-------+-----+----+-----+----+------+----+-------+
    | SLOAD | 16  | 2  | 5.9 | 17 | -6.3 | 14 | -2.93 |
    +-------+-----+----+-----+----+------+----+-------+

    .. note:: Can be used in statics OR dynamics.

    If Si refers to a grid point, the load is applied to component T1 of the
    displacement coordinate system (see the CD field on the GRID entry).
    """
    type = 'SLOAD'
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        nodes = [1]
        mags = [1.]
        return SLOAD(sid, nodes, mags, comment='')

    def __init__(self, sid, nodes, mags, comment=''):
        """
        Creates an SLOAD (GRID/SPOINT load)

        Parameters
        ----------
        sid : int
            load id
        nodes : int; List[int]
            the GRID/SPOINT ids
        mags : float; List[float]
            the load magnitude
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        if isinstance(nodes, integer_types):
            nodes = [nodes]
        if isinstance(mags, float):
            mags = [mags]
        #: load ID
        self.sid = sid
        self.nodes = nodes
        self.mags = mags
        self.nodes_ref = None

    def validate(self):
        assert len(self.nodes) == len(self.mags), 'len(nodes)=%s len(mags)=%s' % (len(self.nodes), len(self.mags))

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
        ngroups = nfields // 2
        if nfields % 2 == 1:
            ngroups += 1
            msg = 'Missing last magnitude on SLOAD card=%s' % card.fields()
            raise RuntimeError(msg)

        nodes = []
        mags = []
        for i in range(ngroups):
            j = 2 * i + 2
            nodes.append(integer(card, j, 'nid' + str(i)))
            mags.append(double(card, j + 1, 'mag' + str(i)))
        return SLOAD(sid, nodes, mags, comment=comment)

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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SLOAD=%s' % (self.sid)
        self.nodes_ref = []
        for nid in self.nodes:
            self.nodes_ref.append(model.Node(nid, msg=msg))
        #self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        return self.cross_reference(model)
        #msg = ', which is required by SLOAD=%s' % (self.sid)
        #self.nodes_ref = model.safe_empty_nodes(self.nodes, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

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
        if self.nodes_ref is None:
            return self.nodes
        return [self.Nid(nid) for nid in self.nodes_ref]

    def raw_fields(self):
        list_fields = ['SLOAD', self.sid]
        for nid, mag in zip(self.node_ids, self.mags):
            list_fields += [nid, mag]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double: bool=False) -> str:
        card = self.raw_fields()
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class RFORCE(Load):
    type = 'RFORCE'
    _properties = ['node_id']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        nid = 1
        scale = 1.
        r123 = [1., 0., 1.]
        return RFORCE(sid, nid, scale, r123,
                      cid=0, method=1, racc=0., mb=0, idrf=0, comment='')

    def __init__(self, sid, nid, scale, r123, cid=0, method=1, racc=0.,
                 mb=0, idrf=0, comment=''):
        """
        idrf doesn't exist in MSC 2005r2; exists in MSC 2016

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
            through point G (R1**2+R2**2+R3**2 > 0 unless A and RACC are
            both zero).
        cid : int; default=0
            Coordinate system defining the components of the rotation vector.
        method : int; default=1
            Method used to compute centrifugal forces due to angular velocity.
        racc : int; default=0.0
            Scale factor of the angular acceleration in revolutions per
            unit time squared.
        mb : int; default=0
            Indicates whether the CID coordinate system is defined in the main
            Bulk Data Section (MB = -1) or the partitioned superelement Bulk
            Data Section (MB = 0). Coordinate systems referenced in the main
            Bulk Data Section are considered stationary with respect to the
            assembly basic coordinate system.
        idrf : int; default=0
            ID indicating to which portion of the structure this particular
            RFORCE entry applies. It is possible to have multiple RFORCE
            entries in the same subcase for SOL 600 to represent different
            portions of the structure with different rotational accelerations.
            IDRF corresponds to a SET3 entry specifying the elements with this
            acceleration. A BRKSQL entry may also be specified with a matching
            IDRF entry.
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.sid = sid
        self.nid = nid
        self.cid = cid
        self.scale = scale
        self.r123 = r123
        self.method = method
        self.racc = racc
        self.mb = mb
        self.idrf = idrf
        self.nid_ref = None
        self.cid_ref = None
        self.validate()

    def validate(self):
        assert self.method in [1, 2], self.method
        assert isinstance(self.r123, list), self.r123
        assert isinstance(self.scale, float), self.scale
        assert isinstance(self.cid, int), self.cid
        assert isinstance(self.mb, int), self.mb
        assert isinstance(self.idrf, int), self.idrf
        assert isinstance(self.racc, float), self.racc

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
        return RFORCE(sid, nid, scale, [r1, r2, r3],
                      cid=cid, method=method, racc=racc, mb=mb, idrf=idrf, comment=comment)

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
        return RFORCE(sid, nid, scale, [r1, r2, r3], cid=cid, method=method, racc=racc, mb=mb,
                      idrf=0, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by RFORCE sid=%s' % self.sid
        if self.nid > 0:
            self.nid_ref = model.Node(self.nid, msg=msg)
        self.cid_ref = model.Coord(self.cid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by RFORCE sid=%s' % self.sid
        if self.nid > 0:
            self.nid_ref = model.Node(self.nid, msg=msg)
        self.cid_ref = model.safe_coord(self.cid, self.sid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nid = self.Nid()
        self.cid = self.Cid()
        self.nid_ref = None
        self.cid_ref = None

    @property
    def node_id(self):
        if self.nid_ref is not None:
            return self.nid_ref.nid
        return self.nid

    def Nid(self):
        return self.node_id

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = (['RFORCE', self.sid, self.node_id, self.Cid(), self.scale] +
                       list(self.r123) + [self.method, self.racc, self.mb, self.idrf])
        return list_fields

    def repr_fields(self):
        #method = set_blank_if_default(self.method,1)
        racc = set_blank_if_default(self.racc, 0.)
        mb = set_blank_if_default(self.mb, 0)
        idrf = set_blank_if_default(self.idrf, 0)
        list_fields = (['RFORCE', self.sid, self.node_id, self.Cid(), self.scale] +
                       list(self.r123) + [self.method, racc, mb, idrf])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        elif is_double:
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
    _properties = ['node_id']

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        nid = 1
        scale = 1.
        group_id = 1
        return RFORCE1(sid, nid, scale, group_id,
                       cid=0, r123=None, racc=0., mb=0, method=2, comment='')

    def __init__(self, sid, nid, scale, group_id,
                 cid=0, r123=None, racc=0., mb=0, method=2, comment=''):
        """
        Creates an RFORCE1 card

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
            through point G (R1**2+R2**2+R3**2 > 0 unless A and RACC are
            both zero).
        racc : int; default=0.0
            Scale factor of the angular acceleration in revolutions per
            unit time squared.
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
        comment : str; default=''
            a comment for the card

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
        self.nid_ref = None
        self.cid_ref = None

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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by RFORCE1 sid=%s' % self.sid
        #if self.nid > 0:  # TODO: why was this every here?
        self.nid_ref = model.Node(self.nid, msg=msg)
        self.cid_ref = model.Coord(self.cid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by RFORCE1 sid=%s' % self.sid
        #if self.nid > 0:  # TODO: why was this every here?
        self.nid_ref = model.Node(self.nid, msg=msg)
        self.cid_ref = model.safe_coord(self.cid, self.sid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nid = self.node_id
        self.cid = self.Cid()
        self.nid_ref = None
        self.cid_ref = None

    @property
    def node_id(self):
        if self.nid_ref is not None:
            return self.nid_ref.nid
        return self.nid

    def Nid(self):
        return self.node_id

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def raw_fields(self):
        list_fields = (['RFORCE1', self.sid, self.node_id, self.Cid(), self.scale]
                       + list(self.r123) + [self.method, self.racc,
                                            self.mb, self.group_id])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
