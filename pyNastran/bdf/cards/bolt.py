from __future__ import annotations
#from math import log, exp
from typing import Optional, TYPE_CHECKING

#import numpy as np
#from numpy import unique, hstack

from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import (
    BaseCard, expand_thru_by, expand_thru,
)
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string,
    #string_or_blank, blank, fields, components_or_blank,
    #integer_string_or_blank, integer_or_double, #parse_components,
    #modal_components_or_blank,
)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class BOLT_MSC(BaseCard):
    """
    +-------+--------+-------+-------+------+------+------+------+------+
    |   1   |   2    |   3   |   4   |  5   |  6   |   7  |  8   |  9   |
    +=======+========+=======+=======+======+======+======+======+======+
    | BOLT  | ID     | GRIDC |       |      |      |      |      |      |
    +-------+--------+-------+-------+------+------+------+------+------+
    |       | TOP    | GT1   |  GT2  |  GT3 |  GT4 |  GT5 |  GT6 |  GT7 |
    +-------+--------+-------+-------+------+------+------+------+------+
    |       | GT8    | GT9   |  etc  |      |      |      |      |      |
    +-------+--------+-------+-------+------+------+------+------+------+
    |       | BOTTOM | GB1   |  GB2  |  GB3 |  GB4 |  GB5 |  GB6 |  GB7 |
    +-------+--------+-------+-------+------+------+------+------+------+
    |       | GB8    | GB9   |  etc  |      |      |      |      |      |
    +-------+--------+-------+-------+------+------+------+------+------+
    """
    type = 'BOLT'
    def __init__(self, bolt_id: int, gridc: int,
                 nids_top: Optional[list[int]]=None,
                 nids_btm: Optional[list[int]]=None,
                 comment: str=''):
        BaseCard.__init__(self)
        self.bolt_id = bolt_id
        self.gridc = gridc
        self.nids_top = nids_top
        self.nids_btm = nids_btm

        #self.nid_ref = None
        #self.nids_ref = None

    @classmethod
    def add_card(cls, card, comment: str=''):
        bolt_id = integer(card, 1, 'bolt_id')
        gridc = integer(card, 2, 'gridc')
        top = string(card, 9, 'top')
        itop = card.index('TOP')
        ibtm = card.index('BOTTOM')
        nids_top = []
        nids_btm = []
        if itop and ibtm:
            assert itop < ibtm, (itop, ibtm)
            fields_top = card[itop+1:ibtm]
            fields_btm = card[ibtm+1:]
            nids_top = expand_thru(fields_top)
            nids_btm = expand_thru(fields_btm)
        elif itop:
            nids_top = card[itop+1:]
        elif ibtm:
            nids_btm = card[ibtm+1:]
        else:
            raise RuntimeError((itop, ibtm))
        assert top == 'TOP', top
        return BOLT_MSC(bolt_id, gridc, nids_top, nids_btm, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def repr_fields(self) -> list[str]:
        card = ['BOLT', self.bolt_id, self.gridc, None, None, None, None, None, None]
        if len(self.nids_top):
            nextra = 8 - (len(self.nids_top) + 1) % 8
            extra = [None] * nextra
            card.extend(['TOP', ] + self.nids_top + extra)
        if len(self.nids_btm):
            card.extend(['BOTTOM', ] + self.nids_btm)
        return card

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            #if self.tid > MAX_INT:
                #return self.comment + print_card_16(card)
            return self.comment + print_card_8(card)
        #if is_double:
            #return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)

class BOLT(BaseCard):
    """
    +-------+------+-------+-------+------+------+------+------+------+
    |   1   |  2   |   3   |   4   |  5   |  6   |   7  |  8   |  9   |
    +=======+======+=======+=======+======+======+======+======+======+
    | BOLT  | BID  | ETYPE | EID1  | EID2 | EID3 | EID4 | EID5 | EID6 |
    +-------+------+-------+-------+------+------+------+------+------+
    |       | EID7 | THRU  | EID8  |  BY  | INC  |      |      |      |
    +-------+------+-------+-------+------+------+------+------+------+
    |       | etc  |       |       |      |      |      |      |      |
    +-------+------+-------+-------+------+------+------+------+------+
    | BOLT  | BID  | ETYPE | CSID  | IDIR | G1   | G2   | G3   | G4   |
    +-------+------+-------+-------+------+------+------+------+------+
    |       | G5   | THRU  | G6    |  BY  | INC  |      |      |      |
    +-------+------+-------+-------+------+------+------+------+------+
    |       | etc  |       |       |      |      |      |      |      |
    +-------+------+-------+-------+------+------+------+------+------+
    | BOLT  | BID  | ETYPE | CSID  | IDIR | GP   |      |      |      |
    +-------+------+-------+-------+------+------+------+------+------+
    |       | EID1 | EID2  | EID3  | EID4 | EID5 | EID6 | EID7 | EID8 |
    +-------+------+-------+-------+------+------+------+------+------+
    |       | EID9 | THRU  | EID10 |  BY  | INC  |      |      |      |
    +-------+------+-------+-------+------+------+------+------+------+
    |       | etc  |       |       |      |      |      |      |      |
    +-------+------+-------+-------+------+------+------+------+------+
    """
    type = 'BOLT'
    def __init__(self, bolt_id: int, element_type: int,
                 eids: Optional[list]=None,  # element_type=1
                 nids: Optional[list]=None,  # element_type=2
                 csid: Optional[int]=None,   # element_type=2/3
                 idir: Optional[int]=None,   # element_type=2/3
                 nid: Optional[int]=None,    # element_type=3
                 comment: str=''):
        """
        bolt_id : int
            Bolt identification number. (Integer > 0)
        ETYPE : int
            Element type. (Integer; No default)
            = 1 to model bolts with CBAR and CBEAM elements in SOLs 101,
              103, 105, 107 through 112, 401, 402, and 601.
            = 2 to model bolts with CHEXA, CPENTA and CTETRA elements in
              SOLs 101, 103, 105, 107 through 112, 401 and 402.
            = 3 to model bolts with CHEXA, CPENTA, CTETRA, CPLSTS3,
              CPLSTS4, CPLSTS6, and CPLSTS8 elements in SOLs 401 and 402,
              or CHEXA, CPENTA, CPYRAM, and CTETRA elements in SOL 601.
        BY Specifies an increment when using THRU option. (Character)
            INC Increment used with THRU option. (Integer; Default = 1)
            INC > 0 can be defined, for example ...106,THRU,126,BY,INC,2
            INC < 0 can be defined, for example ...126,THRU,106,BY,INC,-2

        EIDi Selects element identification numbers to include in the bolt preload
             calculation. See Remark 2 and SOL 601 Remark 1. (Integer > 0, or
             using THRU; EID7 < EID8 for THRU option; No default); element_type=1
        CSID : int; default=???
             Identification number of the coordinate system used to define the bolt
             axis. For the basic coordinate system, CSID = 0. (Integer >= 0; See
             Remark 7 for default behaviour.)
        IDIR : int; default=???
            Direction of bolt axis relative to CSID. (Integer; See Remark 7 for
            default behaviour.)
            = 1 for the X direction
            = 2 for the Y direction
            = 3 for the Z direction
            See SOL 401 Remark 3.
            See SOL 402 Remark 7.
            Gi Identification numbers of grid points that form a cross section through
            the bolt. See Remarks 3 and 4. (Integer ≥ 0; No default)

        GP : int; default=??? (FORM=3)
            For SOL 401, the identification number of the grid point where the bolt
            cross sectional area is calculated. See Remarks 3 and 5. See also
            SOL 401 Remark 4. (Integer > 0 or blank)
            For SOL 402, the identification number of the grid point where the bolt
            cross sectional area is calculated. See Remark 3. See also SOL 402
            Remarks 4 and 5. (Integer > 0; No default)
            For SOL 601, the identification number of the grid point where the bolt
            is split. See Remark 3 and SOL 601 Remarks 3 and 4. (Integer ≥
            0 or blank)
        """
        BaseCard.__init__(self)
        self.bolt_id = bolt_id
        self.element_type = element_type
        self.eids = eids
        self.nids = nids
        self.csid = csid
        self.idir = idir
        if self.element_type in {2, 3}:
            assert self.idir in {None, 0, 1, 2, 3}, idir
        self.nid = None  #  GP

        self.nid_ref = None
        self.eids_ref = None
        self.nids_ref = None

    @classmethod
    def add_card(cls, card, comment: str=''):
        bolt_id = integer(card, 1, 'bolt_id')
        element_type = integer(card, 2, 'element_type')
        assert element_type in {1, 2, 3}, element_type

        eids = []
        nids = []
        idir = None
        csid = None
        nid = None
        if element_type == 1:
            ifield0 = 3
            for ifield in range(6):
                eid = integer_or_blank(card, ifield0+ifield, f'eid{ifield}')
                if eid:
                    eids.append(eid)
            ifield0 += ifield
            assert len(eids), eids
            assert len(card) <= 9, card

        elif element_type == 2:
            csid = integer_or_blank(card, 3, 'csid')
            idir = integer(card, 4, 'idir')
            assert idir in {1, 2, 3}
            ifield0 = 5
            for ifield in range(6):
                nid = integer_or_blank(card, ifield0+ifield, f'nid{ifield}')
                if nid:
                    nids.append(nid)
            ifield0 += ifield
            assert len(nids), nids
            assert len(card) <= 11, card

        elif element_type == 2:
            csid = integer_or_blank(card, 3, 'csid', default=0)
            idir = integer_or_blank(card, 4, 'idir', default=0)
            assert idir in {0, 1, 2, 3}, idir
            assert len(card) <= 5, card
        elif element_type == 3:
            csid = integer_or_blank(card, 3, 'csid', default=0)
            idir = integer_or_blank(card, 4, 'idir', default=0)
            nid = integer_or_blank(card, 5, 'nid/gp', default=0)
            assert idir in {0, 1, 2, 3}, idir

            #ifield = 9
            #while ifield < len(card):
            fields = card[9:]
            eids = expand_thru_by(fields, set_fields=True, sort_fields=True, require_int=True, allow_blanks=False)
            assert len(card) <= 15, card
        else:
            raise NotImplementedError(element_type)
        return BOLT(bolt_id, element_type, eids=eids, nids=nids, csid=csid, idir=idir, nid=nid)

    def cross_reference(self, model: BDF) -> None:
        """xref eids/nids/nid"""
        msg = ', which is required by BOLT bolt_id=%s' % self.bolt_id
        if self.nid is not None:
            self.nid_ref = model.Node(self.nid, msg=msg)

        if self.eids is not None:
            self.eids_ref = model.Elements(self.eids, msg=msg)
        if self.nids is not None:
            self.nids_ref = model.Nodes(self.nids, msg=msg)

    def repr_fields(self) -> list[str]:
        card = ['BOLT', self.bolt_id, self.element_type]
        if self.element_type == 1:
            assert len(self.eids)
            card.extend(self.eids)
        elif self.element_type == 2:
            assert len(self.nids)
            assert self.idir is not None
            card.extend([self.csid, self.idir] + self.nids)
            assert card[4] is not None
        elif self.element_type == 3:
            card.extend([self.csid, self.idir, self.nid, None, None, None] + self.eids)
        else:
            raise NotImplementedError(self.element_type)
        return card

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            #if self.tid > MAX_INT:
                #return self.comment + print_card_16(card)
            return self.comment + print_card_8(card)
        #if is_double:
            #return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class BOLTSEQ(BaseCard):
    """
    BOLTSEQ SID
    S_NOi B_IDi NINCi
    """
    type = 'BOLTSEQ'
    def __init__(self,
                 sid: int,
                 s_nos: list[int],
                 b_ids: list[int],
                 n_incs: Optional[list[int]]=None,
                 comment: str=''):
        """
        SID : int
            Bolt preload set identification number (Integer; No default)
        S_NOi : int
            Sequence order number for the BOLTLD, BOLTFOR, and BOLTFRC
            ID’s to be applied. (Integer; No default)
        B_IDi : int
            SID of BOLTLD, BOLTFOR, or BOLTFRC bulk entries defining a bolt
             preload. (Integer; No default)
        NINCi : int; default=1
            Number of increments in which to ramp up the bolt preload defined in
            BOLTLD, BOLTFOR, or BOLTFRC entries. (Integer; Default = 1)
        """
        BaseCard.__init__(self)
        self.sid = sid
        if isinstance(s_nos, integer_types):
            s_nos = [s_nos]
        if isinstance(b_ids, integer_types):
            b_ids = [b_ids] * len(s_nos)
        assert len(s_nos) == len(b_ids)

        if isinstance(n_incs, integer_types):
            n_incs = [n_incs] * len(s_nos)
        self.s_nos = s_nos
        self.b_ids = b_ids
        self.n_incs = n_incs

    @classmethod
    def add_card(cls, card, comment: str=''):
        sid = integer(card, 1, 'sid')
        ifield = 9
        s_nos = []
        b_ids = []
        n_incs = []
        while ifield < len(card):
            s_no = integer(card, ifield, 's_no')
            b_id = integer(card, ifield+1, 'b_id')
            n_inc = integer_or_blank(card, ifield+2, 'n_inc', default=1)
            s_nos.append(s_no)
            b_ids.append(b_id)
            n_incs.append(n_inc)
            ifield += 8
        assert len(s_nos) >= 1, s_nos
        assert len(b_ids) >= 1, b_ids
        assert len(n_incs) >= 1, n_incs
        return BOLTSEQ(sid, s_nos, b_ids, n_incs=n_incs)

    def cross_reference(self, model: BDF) -> None:
        self.so_nos_ref = []
        self.b_ids_ref = []

        for bolt_id in self.s_nos:
            boltfor = None
            boltfrc = None
            boltld = None
            # BOLTLD, BOLTFOR, and BOLTFRC
            if bolt_id in model.boltld:
                boltld = model.boltld[bolt_id]
            if bolt_id in model.boltfor:
                boltfor = model.boltfor[bolt_id]
            if bolt_id in model.boltfrc:
                boltfrc = model.boltfrc[bolt_id]
            bolts = (boltld, boltfor, boltfrc)
            self.so_nos_ref.append(bolts)

        for bolt_id in self.b_ids:
            boltfor = None
            boltfrc = None
            boltld = None
            # BOLTLD, BOLTFOR, and BOLTFRC
            if bolt_id in model.boltld:
                boltld = model.boltld[bolt_id]
            if bolt_id in model.boltfor:
                boltfor = model.boltfor[bolt_id]
            if bolt_id in model.boltfrc:
                boltfrc = model.boltfrc[bolt_id]
            bolts = (boltld, boltfor, boltfrc)
            self.b_ids_ref.append(bolts)
        #Sequence order number for the BOLTLD, BOLTFOR, and BOLTFRC
        #IDs to be applied. (Integer; No default)
        #B_IDi SID of BOLTLD, BOLTFOR, or BOLTFRC bulk entries defining a bolt
        #preload. (Integer; No default)
        pass

    def repr_fields(self) -> list[str]:

        fields = ['BOLTSEQ', self.sid, None, None, None, None, None, None, None]
        for s_no, b_id, n_inc in zip(self.s_nos, self.b_ids, self.n_incs):
            fields.extend([s_no, b_id, n_inc, None, None, None, None, None])
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            #if self.tid > MAX_INT:
                #return self.comment + print_card_16(card)
            return self.comment + print_card_8(card)
        #if is_double:
            #return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class BOLTFOR(BaseCard):
    """
    BOLTFOR SID LOAD B1 B2 B3 B4 B5 B6
    B7 THRU B8
    B9 B10 -etc-
    """
    type = 'BOLTFOR'
    def __init__(self,
                 sid: int,
                 load_value: float,
                 bolt_ids: list[int],
                 comment: str=''):
        BaseCard.__init__(self)
        self.sid = sid
        self.load_value = load_value
        self.bolt_ids = list(set(bolt_ids))
        bolt_ids.sort()
        self.bolt_ids_ref = []
        assert len(bolt_ids) > 0, bolt_ids

    @classmethod
    def add_card(self, card, comment: str=''):
        sid = integer(card, 1, 'sid')
        load_value = double(card, 2, 'load')
        fields = card[3:9]

        bolt_ids = expand_thru_by(fields, set_fields=False, sort_fields=False, require_int=True, allow_blanks=True)
        ifield = 9
        while ifield < len(card):
            fields = card[ifield:ifield+8]
            bolt_ids += expand_thru_by(fields, set_fields=False, sort_fields=False, require_int=True, allow_blanks=True)
            ifield += 8
        assert len(bolt_ids) >= 1, card
        return BOLTFOR(sid, load_value, bolt_ids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        self.bolt_ids_ref = []
        for bolt_id in self.bolt_ids:
            boltld = None
            boltfor = None
            boltfrc = None
            boltseq = None

            # SOL 401: BOLTLD, BOLTFOR, BOLTFRC, or BOLTSEQ
            if bolt_id in model.boltld:
                boltld = model.boltld[bolt_id]
            if bolt_id in model.boltfor:
                boltfor = model.boltfor[bolt_id]
            if bolt_id in model.boltfrc:
                boltfrc = model.boltfrc[bolt_id]
            if bolt_id in model.boltseq:
                boltseq = model.boltseq[bolt_id]
            bolts = (boltld, boltfor, boltfrc, boltseq)
            self.bolt_ids_ref.append(bolts)

    def repr_fields(self) -> list[str]:
        fields = ['BOLTFOR', self.sid, self.load_value] + self.bolt_ids
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            #if self.tid > MAX_INT:
                #return self.comment + print_card_16(card)
            return self.comment + print_card_8(card)
        #if is_double:
            #return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)

BOLTFRC = BOLT
BOLTLD = BOLT
