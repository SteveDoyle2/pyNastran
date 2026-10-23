"""
Defines the following contact cards:

nx contact:
 - BCONP (see contact_common.py)
 - BLSEG
 - BCPARA (see contact_common.py)
 - BCRPARA
 - BCTPARA
 - BCTADD
 - BCTSET
 - BSURF
 - BSURFS
 - BFRIC

"""
from __future__ import annotations
import warnings
from itertools import zip_longest
from typing import Optional, Any, TYPE_CHECKING

from pyNastran.bdf.cards.base_card import BaseCard, expand_thru_by, _node_ids
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_string_or_blank, double_or_blank,
    integer_double_or_blank, string, string_or_blank, string_choice_or_blank,
    double)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class BFRIC(BaseCard):
    """
    Slideline Contact Friction
    Defines frictional properties between two bodies in contact.

    +-------+------+-------+-------+-------+
    |   1   |   2  |   3   |   4   |   5   |
    +=======+======+=======+=======+=======+
    | BFRIC | FID  | FSTIF |       |  MU1  |
    +-------+------+-------+-------+-------+

    """
    type = 'BFRIC'

    @classmethod
    def _init_from_empty(cls):
        friction_id = 1
        mu1 = 0.2
        return BFRIC(friction_id, mu1)

    def __init__(self, friction_id: int, mu1: float, fstiff=None, comment: str=''):
        """
        Creates a BFRIC card, which defines a frictional contact.

        Parameters
        ----------
        friction_id : int
            Friction identification number.
        mu1 : float
            Coefficient of static friction.
        fstiff : float; default=None
            Frictional stiffness in stick. See Remarks 2 and 3
            Default=automatically selected by the program.

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.friction_id = friction_id
        self.fstiff = fstiff
        self.mu1 = mu1

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        friction_id = integer(card, 1, 'friction_id')
        fstiff = double_or_blank(card, 2, 'fstiff')
        #
        mu1 = double(card, 4, 'mu1')
        return BFRIC(friction_id, mu1, fstiff=fstiff, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
        # msg = f', which is required for BSURFS={self.sid:d}'
        # model.log.warning(f'cross-reference not supported for {self.type}')

    def raw_fields(self):
        list_fields = [
            'BFRIC', self.friction_id, self.fstiff, None, self.mu1]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class BLSEG(BaseCard):
    """
    3D Contact Region Definition by Shell Elements (SOLs 101, 601 and 701)

    Defines a 3D contact region by shell element IDs.

    +=======+====+====+======+====+====+=====+====+====+
    |   1   |  2 |  3 |   4  |  5 |  6 |  7  |  8 |  9 |
    +-------+----+----+------+----+----+-----+----+----+
    | BLSEG | ID | G1 |  G2  | G3 | G4 | G5  | G6 | G7 |
    +-------+----+----+------+----+----+-----+----+----+
    | BLSEG | ID | G1 | THRU | G2 | BY | INC |    |    |
    +-------+----+----+------+----+----+-----+----+----+

    """
    type = 'BLSEG'

    @classmethod
    def _init_from_empty(cls):
        line_id = 1
        nodes = [1]
        return BLSEG(line_id, nodes)

    def __init__(self, line_id, nodes, comment: str=''):
        if comment:
            self.comment = comment

        self.line_id = line_id
        self.nodes = expand_thru_by(nodes)
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BLSEG card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        line_id = integer(card, 1, 'line_id')
        #: Number (float)
        nfields = card.nfields
        i = 2
        nodes = []
        while i < nfields:
            d = integer_string_or_blank(card, i, 'field_%s' % i)
            if d is not None:
                nodes.append(d)
            i += 1
        return BLSEG(line_id, nodes, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = f', which is required by BLSEG line_id={self.line_id}'
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)

    def uncross_reference(self) -> None:
        #msg = f', which is required by BLSEG line_id={self.line_id}'
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self):
        """returns nodeIDs for repr functions"""
        return _node_ids(self, nodes=self.nodes_ref, allow_empty_nodes=False, msg='')

    def raw_fields(self):
        list_fields = ['BLSEG', self.line_id] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class BSURF(BaseCard):
    """
    3D Contact Region Definition by Shell Elements (SOLs 101, 601 and 701)

    Defines a 3D contact region by shell element IDs.

    +-------+------+------+-------+-------+--------+------+------+------+
    |   1   |   2  |   3  |   4   |   5   |    6   |  7   |   8  |   9  |
    +=======+======+======+=======+=======+========+======+======+======+
    | BSURF |  ID  | EID1 | EID2  | EID3  |  EID4  | EID5 | EID6 | EID7 |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       | EID8 | EID9 | EID10 |  etc. |        |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    | BSURF |  ID  | EID1 |  THRU | EID2  |   BY   | INC  |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       | EID8 | EID9 | EID10 | EID11 |  etc.  |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       | EID8 | THRU | EID9  |  BY   |  INC   |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    | BSURF |  15  |  5   | THRU  |  21   |   BY   |  4   |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       |  27  |  30  |  32   |  33   |        |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       |  35  | THRU |  44   |       |        |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+
    |       |  67  |  68  |  70   |  85   |   92   |      |      |      |
    +-------+------+------+-------+-------+--------+------+------+------+

    """
    type = 'BSURF'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        eids = [1]
        return BSURF(sid, eids)

    def __init__(self, sid, eids, comment: str=''):
        if comment:
            self.comment = comment
        #: Set identification number. (Unique Integer > 0)
        self.sid = sid
        #: Element identification numbers of shell elements. (Integer > 0)
        self.eids = eids
        self.eids_ref = []

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BSURF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        #: Number (float)
        nfields = card.nfields
        i = 2
        eid_data = []
        while i < nfields:
            d = integer_string_or_blank(card, i, 'field_%s' % i)
            if d is not None:
                eid_data.append(d)
            i += 1
        eids = expand_thru_by(eid_data)
        return BSURF(sid, eids, comment=comment)

    def cross_reference(self, model: BDF):
        msg = f', which is required for BSURF={self.sid:d}'
        eids_ref = []
        for eid in self.eids:
            eid_ref = model.Element(eid, msg=msg)
            eids_ref.append(eid_ref)
        self.eids_ref = eids_ref

    def raw_fields(self):
        fields = ['BSURF', self.sid]
        return fields + list(self.eids)

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class BSURFS(BaseCard):
    """
    Defines a 3D contact region by the faces of the CHEXA, CPENTA or CTETRA
    elements.

    Notes
    -----
    1. The continuation field is optional.
    2. BSURFS is a collection of one or more element faces on solid elements.
       BSURFS defines a contact region which may act as a contact source
       (contactor) or target.
    3. The ID must be unique with respect to all other BSURFS, BSURF, and
       BCPROP entries.

    """
    type = 'BSURFS'

    @classmethod
    def _init_from_empty(cls):
        bsurfs_id = 1
        eids = [1]
        g1s = [1]
        g2s = [1]
        g3s = [1]
        return BSURFS(bsurfs_id, eids, g1s, g2s, g3s)

    def __init__(self, bsurfs_id: int,
                 eids: list[int],
                 g1s: list[int],
                 g2s: list[int],
                 g3s: list[int],
                 comment: str=''):
        if comment:
            self.comment = comment
        #: Identification number of a contact region. See Remarks 2 and 3.
        #: (Integer > 0)
        self.id = bsurfs_id

        #: Element identification numbers of solid elements. (Integer > 0)
        self.eids = eids

        #: Identification numbers of 3 corner grid points on the face (triangular
        #: or quadrilateral) of the solid element. (Integer > 0)
        self.g1s = g1s
        self.g2s = g2s
        self.g3s = g3s
        self.grids_ref = []

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BSURFS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        bsurfs_id = integer(card, 1, 'id')
        eids = []
        g1s = []
        g2s = []
        g3s = []

        n = card.nfields - 5
        i = 0
        j = 1
        while i < n:
            eid = integer(card, 5 + i, 'eid%s' % j)
            g1 = integer(card, 5 + i + 1, 'g3_%s' % j)
            g2 = integer(card, 5 + i + 2, 'g2_%s' % j)
            g3 = integer(card, 5 + i + 3, 'g1_%s' % j)
            j += 1
            i += 4
            eids.append(eid)
            g1s.append(g1)
            g2s.append(g2)
            g3s.append(g3)
        return BSURFS(bsurfs_id, eids, g1s, g2s, g3s, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = f', which is required for BSURFS={self.id:d}'
        # model.log.warning(f'cross-reference not supported for {self.type}')
        eids_ref = []
        grids_ref = []
        for eid in self.eids:
            elem = model.Element(eid, msg=msg)
            eids_ref.append(elem)
        for g1, g2, g3 in zip_longest(self.g1s, self.g2s, self.g3s):
            g1_ref = model.Node(g1, msg=msg)
            g2_ref = model.Node(g2, msg=msg)
            g3_ref = model.Node(g3, msg=msg)
            grids_ref.append((g1_ref, g2_ref, g3_ref))
        self.eids_ref = eids_ref
        self.grids_ref = grids_ref

    def grids(self) -> list[tuple[int, int, int]]:
        if len(self.grids_ref) == 0:
            return zip(self.g1s, self.g2s, self.g3s)
        nodes_out = []
        for g1_ref, g2_ref, g3_ref in self.grids_ref:
            nodes = (g1_ref.nid, g2_ref.nid, g3_ref.nid)
            nodes_out.append(nodes)
        return nodes_out

    def raw_fields(self):
        fields = ['BSURFS', self.id, None, None, None]
        for eid, (g1, g2, g3) in zip(self.eids,self.grids()):
            fields += [eid, g1, g2, g3]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class BCTSET(BaseCard):
    """
    3D Contact Set Definition (SOLs 101, 601 and 701 only)
    Defines contact pairs of a 3D contact set.

    +--------+-------+------+-------+-------+-------+-------+
    |   1    |   2   | 3    |  4    |   5   |   6   |   7   |
    +========+=======+======+=======+=======+=======+=======+
    | BCTSET | CSID  | SID1 | TID1  | FRIC1 | MIND1 | MAXD1 |
    +--------+-------+------+-------+-------+-------+-------+
    |        |       | SID2 | TID2  | FRIC2 | MIND2 | MAXD2 |
    +--------+-------+------+-------+-------+-------+-------+
    |        |  etc. |      |       |       |       |       |
    +--------+-------+------+-------+-------+-------+-------+

    """
    type = 'BCTSET'

    @classmethod
    def _init_from_empty(cls):
        csid = 1
        sids = [1]
        tids = [1]
        frictions = [0.01]
        min_distances = [0.1]
        max_distances = [1.]
        return BCTSET(csid, sids, tids, frictions, min_distances, max_distances, sol=101)

    def __init__(self, csid, sids, tids, frictions, min_distances, max_distances,
                 comment: str='', sol: int=101):
        if comment:
            self.comment = comment
        #: CSID Contact set identification number. (Integer > 0)
        self.csid = csid
        #: SIDi Source region (contactor) identification number for contact pair i.
        #: (Integer > 0)
        self.sids = sids

        #: TIDi Target region identification number for contact pair i. (Integer > 0)
        self.tids = tids

        #: FRICi Static coefficient of friction for contact pair i. (Real; Default=0.0)
        self.frictions = frictions

        #: MINDi Minimum search distance for contact. (Real) (Sol 101 only)
        self.min_distances = min_distances

        #: MAXDi Maximum search distance for contact. (Real) (Sol 101 only)
        self.max_distances = max_distances

    @classmethod
    def add_card(cls, card, comment: str='', sol: int=101):
        csid = integer(card, 1, 'csid')
        sids = []
        tids = []
        frictions = []
        min_distances = []
        max_distances = []

        nfields = card.nfields
        i = 2
        j = 1
        while i < nfields:
            sids.append(integer(card, i, 'sid%s' % j))
            tids.append(integer(card, i + 1, 'tid%s' % j))
            frictions.append(double_or_blank(card, i + 2, 'fric%s' % j, default=0.0))

            # we ignore what the NX QRG says about the min/max distance for SOL 401
            # if you don't, Nastran will crash
            #if sol == 101:
            min_distances.append(double_or_blank(card, i + 3, 'mind%s' % j, default=0.0))
            max_distances.append(double_or_blank(card, i + 4, 'maxd%s' % j, default=0.0))
            #else:
                #min_distances.append(None)
                #max_distances.append(None)

            i += 8
            j += 1
        return BCTSET(csid, sids, tids, frictions, min_distances,
                      max_distances, comment=comment,
                      sol=sol)

    def cross_reference(self, model: BDF) -> None:
        # msg = f', which is required for BSURF={self.csid:d}'
        model.log.warning(f'cross-reference not supported for {self.type}')
        #: SIDi Source region (contactor) identification number for contact pair i.
        #: (Integer > 0)

        sids_ref = []
        tids_ref = []
        for sid in self.sids:
            # BCPROP, BCPROPS, BSURF, or BSURFS
            sid_ref = contact_surface(model, self, self.csid, msg='target contact surface')
            if sid_ref is None:
                continue
            sids_ref.append(sid_ref)

        #: TIDi Target region identification number for contact pair i. (Integer > 0)
        for tid in self.tids:
            # BCPROP, BCPROPS, BSURF, or BSURFS
            tid_ref = contact_surface(model, self, self.csid, msg='target contact surface')
            if tid_ref is None:
                continue
            tids_ref.append(sid_ref)

    def raw_fields(self):
        fields = ['BCTSET', self.csid]
        for sid, tid, fric, mind, maxd in zip(self.sids, self.tids, self.frictions,
                                              self.min_distances, self.max_distances):
            fields += [sid, tid, fric, mind, maxd, None, None, None]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BCRPARA(BaseCard):
    """
    +---------+------+------+--------+------+-----+---+---+---+----+
    |    1    |   2  |   3  |   4    |   5  |  6  | 7 | 8 | 9 | 10 |
    +=========+======+======+========+======+=====+===+===+===+====+
    | BCRPARA | CRID | SURF | OFFSET | TYPE | GP  |   |   |   |    |
    +---------+------+------+--------+------+-----+---+---+---+----+
    """
    type = 'BCRPARA'

    @classmethod
    def _init_from_empty(cls):
        crid = 1
        return BCRPARA(crid, offset=None, surf='TOP', Type='FLEX', grid_point=0)

    def __init__(self, crid: int, offset: Optional[float]=None,
                 surf: str='TOP', Type: str='FLEX', grid_point: int=0,
                 comment: str=''):
        """
        Creates a BCRPARA card

        Parameters
        ----------
        crid : int
            CRID Contact region ID.
        offset : float; default=None
            Offset distance for the contact region (Real > 0.0).
            None : OFFSET value in BCTPARA entry
        surf : str; default='TOP'
            SURF Indicates the contact side. See Remark 1.  {'TOP', 'BOT'; )
        Type : str; default='FLEX'
            Indicates whether a contact region is a rigid surface if it
            is used as a target region. {'RIGID', 'FLEX'}.
            This is not supported for SOL 101.
        grid_point : int; default=0
            Control grid point for a target contact region with TYPE=RIGID
            or when the rigid-target algorithm is used.  The grid point
            may be used to control the motion of a rigid surface.
            (Integer > 0).  This is not supported for SOL 101.
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: CRID Contact region ID. (Integer > 0)
        self.crid = crid

        #: SURF Indicates the contact side. See Remark 1. (Character = "TOP" or
        #: "BOT"; Default = "TOP")
        self.surf = surf

        #: Offset distance for the contact region. See Remark 2. (Real > 0.0,
        #: Default =OFFSET value in BCTPARA entry)
        self.offset = offset

        #: Indicates whether a contact region is a rigid surface if it is used as a
        #: target region. See Remarks 3 and 4. (Character = "RIGID" or "FLEX",
        #: Default = "FLEX"). This is not supported for SOL 101.
        self.Type = Type

        #: Control grid point for a target contact region with TYPE=RIGID or
        #: when the rigid-target algorithm is used. The grid point may be
        #: used to control the motion of a rigid surface. (Integer > 0)
        #: This is not supported for SOL 101.
        self.grid_point = grid_point

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCRPARA card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        crid = integer(card, 1, 'crid')
        surf = string_or_blank(card, 2, 'surf', default='TOP')
        offset = double_or_blank(card, 3, 'offset', default=None)
        Type = string_or_blank(card, 4, 'type', default='FLEX')
        grid_point = integer_or_blank(card, 5, 'grid_point', default=0)
        return BCRPARA(crid, surf=surf, offset=offset, Type=Type,
                       grid_point=grid_point, comment=comment)

    def cross_reference(self, model: BDF):
        pass

    def raw_fields(self):
        fields = ['BCRPARA', self.crid, self.surf, self.offset, self.Type, self.grid_point]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BCTPARM(BaseCard):
    """
    Contact Parameters (SOLs 101, 103, 111, 112, and 401).
    Control parameters for the contact algorithm.

    +---------+--------+--------+--------+--------+--------+---------+--------+
    |    1    |   2    |    3   |   4    |   5    |   6    |    7    |    8   |
    +=========+========+========+========+========+========+=========+========+
    | BCTPARM | CSID   | Param1 | Value1 | Param2 | Value2 | Param3  | Value3 |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | Param4 | Value4 | Param5 | Value5 |  etc.  |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    | BCTPARM |   1    | PENN   |  10.0  |  PENT  |  0.5   |   CTOL  | 0.001  |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | SHLTHK |   1    |        |        |        |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+

    """
    type = 'BCTPARM'

    @classmethod
    def _init_from_empty(cls):
        csid = 1
        params = {'CSTIFF' : 1}
        return BCTPARM(csid, params)

    def _finalize_hdf5(self, encoding: str):
        """hdf5 helper function"""
        keys, values = self.params
        self.params = {key : value for key, value in zip(keys, values)}

    def __init__(self, csid: int,
                 params: dict[str, int | float],
                 comment: str=''):
        """
        Creates a BCTPARM card

        Parameters
        ----------
        csid : int
            Contact set ID. Parameters defined in this command apply to
            contact set CSID defined by a BCTSET entry. (Integer > 0)
        params : dict[key] : value
            the optional parameters
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: Contact set ID. Parameters defined in this command apply to
        #: contact set CSID defined by a BCTSET entry. (Integer > 0)
        self.csid = csid
        self.params = params

    def _finalize_hdf5(self, encoding: str) -> None:
        keys = self.params[0]
        values = self.params[1]
        self.params = {}
        for key, value in zip(keys, values):
            if isinstance(value, bytes):
                value = value.decode(encoding)
            self.params[key] = value

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCTPARM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        csid = integer(card, 1, 'csid')
        params = read_param_dict(card, cls._get_param_value_from_card)
        return BCTPARM(csid, params, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def _get_param_value_from_card(
            card: BDFCard,
            i: int, j: int) -> tuple[str, int | float]:
        param = string(card, i, 'param%s' % j)
        if param == 'TYPE' and 0:
            value = integer_or_blank(card, i+1, f'value{j:d}', default=0)
            assert value in [0, 1, 2], 'TYPE must be [0, 1, 2]; TYPE=%r' % value
        elif param == 'PENN':
            # PENN 10.0
            value = double(card, i+1, 'value%s' % j)
        elif param == 'PENT':
            # PENT 0.5
            value = double(card, i+1, 'value%s' % j)
        elif param == 'CTOL':
            # CTOL 10.0
            value = double(card, i+1, 'value%s' % j)
        elif param == 'SHLTHK':
            # SHLTHK 1
            value = integer(card, i+1, 'value%s' % j)
        # elif param == 'TYPE': # NX
        # value = string_or_blank(card, i, 'value%s' % j, default='FLEX').upper()
        # assert value in ['FLEX', 'RIGID', 'COATING'], 'TYPE must be [FLEX, RIGID, COATING.]; CSTIFF=%r' % value

        # elif param == 'NSIDE':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=1)
        # assert value in [1, 2], 'NSIDE must be [1, 2]; NSIDE=%r' % value
        # elif param == 'TBIRTH':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)
        # elif param == 'TDEATH':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)
        # elif param == 'INIPENE':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=0)
        # assert value in [0, 1, 2, 3], 'INIPENE must be [0, 1, 2]; INIPENE=%r' % value
        # elif param == 'PDEPTH':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)
        # elif param == 'SEGNORM':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=0)
        # assert value in [-1, 0, 1], 'SEGNORM must be [-1, 0, 1]; SEGNORM=%r' % value
        # elif param == 'OFFTYPE':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=0)
        # assert value in [0, 1, 2], 'OFFTYPE must be [0, 1, 2]; OFFTYPE=%r' % value
        # elif param == 'OFFSET':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)
        # elif param == 'TZPENE':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)

        # elif param == 'CSTIFF':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=0)
        # assert value in [0, 1], 'CSTIFF must be [0, 1]; CSTIFF=%r' % value
        # elif param == 'TIED':
        # value = integer_or_blank(card, i+1, 'value%s' % j, default=0)
        # assert value in [0, 1], 'TIED must be [0, 1]; TIED=%r' % value
        # elif param == 'TIEDTOL':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.0)
        # elif param == 'EXTFAC':
        # value = double_or_blank(card, i+1, 'value%s' % j, default=0.001)
        # assert 1.0E-6 <= value <= 0.1, 'EXTFAC must be 1.0E-6 < EXTFAC < 0.1; EXTFAC=%r' % value
        else:
            # FRICMOD, FPARA1/2/3/4/5, EPSN, EPST, CFACTOR1, PENETOL
            # NCMOD, TCMOD, RFORCE, LFORCE, RTPCHECK, RTPMAX, XTYPE
            # ...
            value = integer_double_or_blank(card, i+1, 'value%s' % j)
            assert value is not None, '%s%i must not be None' % (param, j)
        return param, value

    def raw_fields(self) -> list[str | int | float]:
        fields = ['BCTPARM', self.csid]
        fields.extend(params_listify(self.params))
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class BCTPARA(BaseCard):
    """
    Defines parameters for a surface-to-surface contact region.

    +---------+--------+--------+--------+--------+--------+---------+--------+
    |    1    |   2    |    3   |   4    |   5    |   6    |    7    |    8   |
    +=========+========+========+========+========+========+=========+========+
    | BCTPARA | CSID   | Param1 | Value1 | Param2 | Value2 | Param3  | Value3 |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | Param4 | Value4 | Param5 | Value5 |  etc.  |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    | BCTPARA |   1    | TYPE   |   0    | NSIDE  |   2    | SEGNORM |  -1    |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | CSTIFF |   1    | OFFSET | 0.015  |        |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+

    """
    type = 'BCTPARA'

    @classmethod
    def _init_from_empty(cls):
        csid = 1
        params = {'CSTIFF' : 1}
        return BCTPARA(csid, params)

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        keys, values = self.params
        self.params = {key : value for key, value in zip(keys, values)}

    def __init__(self, csid: int,
                 params: dict[str, Any], comment: str=''):
        """
        Creates a BCTPARA card

        Parameters
        ----------
        csid : int
            Contact set ID. Parameters defined in this command apply to
            contact set CSID defined by a BCTSET entry. (Integer > 0)
        params : dict[key] : value
            the optional parameters
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: Contact set ID. Parameters defined in this command apply to
        #: contact set CSID defined by a BCTSET entry. (Integer > 0)
        self.csid = csid
        self.params = params

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCTPARA card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        csid = integer(card, 1, 'csid')

        all_params = [
            'TYPE', 'NSIDE', 'TBIRTH', 'TDEATH', 'INIPENE',
            'PDEPTH', 'SEGNORM', 'OFFTYPE', 'OFFSET',
            'TZPENE', 'CSTIFF', 'TIED', 'TIEDTOL', 'EXTFAC']
        all_params.sort()
        params = read_param_dict(card, cls._get_param_value_from_card)
        for param in params:
            if param not in all_params:
                warnings.warn(f'param={param!r} is not allowed; allowed={all_params}\n{str(card)}')

        return BCTPARA(csid, params, comment=comment)

    @staticmethod
    def _get_param_value_from_card(card: BDFCard,
                                   i: int, j: int) -> tuple[str, int | float]:
        param = string(card, i, f'param{j}')
        # assert param in all_params, f'param={param!r} is not allowed; allowed={all_params}'
        i += 1
        if param == 'TYPE':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [0, 1, 2], 'TYPE must be [0, 1, 2]; TYPE=%r' % value
        #elif param == 'TYPE': # NX
            #value = string_or_blank(card, i, 'value%s' % j, 'FLEX').upper()
            #assert value in ['FLEX', 'RIGID', 'COATING'], 'TYPE must be [FLEX, RIGID, COATING.]; CSTIFF=%r' % value

        elif param == 'NSIDE':
            value = integer_or_blank(card, i, 'value%s' % j, default=1)
            assert value in [1, 2], 'NSIDE must be [1, 2]; NSIDE=%r' % value
        elif param == 'TBIRTH':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)
        elif param == 'TDEATH':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)
        elif param == 'INIPENE':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [0, 1, 2, 3], 'INIPENE must be [0, 1, 2]; INIPENE=%r' % value
        elif param == 'PDEPTH':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)
        elif param == 'SEGNORM':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [-1, 0, 1], 'SEGNORM must be [-1, 0, 1]; SEGNORM=%r' % value
        elif param == 'OFFTYPE':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [0, 1, 2], 'OFFTYPE must be [0, 1, 2]; OFFTYPE=%r' % value
        elif param == 'OFFSET':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)
        elif param == 'TZPENE':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)

        elif param == 'CSTIFF':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [0, 1], 'CSTIFF must be [0, 1]; CSTIFF=%r' % value
        elif param == 'TIED':
            value = integer_or_blank(card, i, 'value%s' % j, default=0)
            assert value in [0, 1], 'TIED must be [0, 1]; TIED=%r' % value
        elif param == 'TIEDTOL':
            value = double_or_blank(card, i, 'value%s' % j, default=0.0)
        elif param == 'EXTFAC':
            value = double_or_blank(card, i, 'value%s' % j, default=0.001)
            assert 1.0E-6 <= value <= 0.1, 'EXTFAC must be 1.0E-6 < EXTFAC < 0.1; EXTFAC=%r' % value
        else:
            # FRICMOD, FPARA1/2/3/4/5, EPSN, EPST, CFACTOR1, PENETOL
            # NCMOD, TCMOD, RFORCE, LFORCE, RTPCHECK, RTPMAX, XTYPE
            # ...
            value = integer_double_or_blank(card, i, 'value%s' % j)
            assert value is not None, '%s%i must not be None' % (param, j)
        return param, value

    def cross_reference(self, model: BDF) -> None:
        pass

    def raw_fields(self) -> list:
        fields = ['BCTPARA', self.csid]
        fields.extend(params_listify(self.params))
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BCTADD(BaseCard):
    """
    +--------+------+----+-------+----+----+----+----+----+
    |   1    |  2   | 3  |   4   |  5 | 6  |  7 | 8  |  9 |
    +========+======+====+=======+====+====+====+====+====+
    | BCTADD | CSID | SI |  S2   | S3 | S4 | S5 | S6 | S7 |
    +--------+------+----+-------+----+----+----+----+----+
    |        |  S8  | S9 |  etc. |    |    |    |    |    |
    +--------+------+----+-------+----+----+----+----+----+

    Remarks:
    1. To include several contact sets defined via BCTSET entries in a model,
       BCTADD must be used to combine the contact sets. CSID in BCTADD is
       then selected with the Case Control command BCSET.
    2. Si must be unique and may not be the identification of this or any other
       BCTADD entry.

    """
    type = 'BCTADD'

    @classmethod
    def _init_from_empty(cls):
        csid = 1
        contact_sets = [1, 2]
        return BCTADD(csid, contact_sets)

    def __init__(self, csid, contact_sets, comment: str=''):
        if comment:
            self.comment = comment
        #: Contact set identification number. (Integer > 0)
        self.csid = csid

        #: Identification numbers of contact sets defined via BCTSET entries.
        #: (Integer > 0)
        self.contact_sets = contact_sets
        self.contact_sets_ref = []

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BCTADD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        csid = integer(card, 1, 'csid')
        contact_sets = []

        i = 1
        j = 1
        while i < card.nfields:
            contact_set = integer(card, i, 'S%i' % j)
            contact_sets.append(contact_set)
            i += 1
            j += 1
        return BCTADD(csid, contact_sets, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = f', which is required for {self.type}={self.csid:d}'
        contact_sets_ref = []
        bctsets = list(model.bctsets)
        for contact_set in self.contact_sets:
            try:
                contact_set_ref = model.bctsets[contact_set]
            except KeyError:
                # print(model.get_bdf_stats(contact_set))
                model.log.error(f'Invalid cross-reference for {self.type}={self.csid}\n'
                                f'bctsets={bctsets}\n{self}')
                continue
            contact_sets_ref.append(contact_set_ref)
        self.contact_sets_ref = contact_sets_ref

    def raw_fields(self):
        fields = ['BCTADD'] + self.contact_sets
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


def params_listify(params: dict[str, int | float],
                   ) -> list[str | int | float]:
    """
    Annoying function for contact params.

    3 sets of params are allowed per line, but:
     - 1 blank on the first line
     - 2 blanon the 2nd/3rd/etc. line

    """
    i = 0
    fields = []
    for j, (key, value) in enumerate(sorted(params.items())):
        if j == 3:
            fields.append(None)
            i = 0
        elif j > 0 and i == 3:
            fields.append(None)
            fields.append(None)
            i = 0
        fields.append(key)
        fields.append(value)
        i += 1
    return fields


def read_param_dict(card: BDFCard,
                    read_func) -> dict[str, int | float]:
    """reading version of params_listify"""
    params = {}
    iline = 0
    ifield = 2
    iparam = 0
    while ifield < card.nfields:
        param, value = read_func(card, ifield, iparam+1)
        params[param] = value
        iparam += 1
        if iparam % 3 == 0:
            # every 3 params, reset the line
            ifield = 9 + 8 * iline
            iline += 1
        else:
            ifield += 2
    assert len(params) > 0, card
    return params


def contact_surface(model: BDF, card, sid: int, msg: str=''):
    # 3D: BCPROP, BCPROPS, BSURF, or BSURFS
    # 2D: BEDGE or BLSEG
    if sid in model.bsurf:
        sid_ref = model.bsurf[sid]
    elif sid in model.bsurfs:
        sid_ref = model.bsurfs[sid]
    # elif sid in model.bedge:
    #     sid_ref = model.bedge[sid]
    elif sid in model.blseg:
        sid_ref = model.blseg[sid]
    # elif sid in model.bcprop:
    #     sid_ref = model.bcprop[sid]
    # elif sid in model.bcprops:
    #     sid_ref = model.bcprops[sid]
    else:
        # print(model.get_bdf_stats())
        bsurf = list(model.bsurf)
        bsurfs = list(model.bsurfs)
        bedge = [] # list(model.bedge)
        blseg = list(model.blseg)
        model.log.error(f'missing {msg} for {card.type}={sid}\n'
                        f'bsurf: {bsurf}\nbsurfs: {bsurfs}\n'
                        f'bedge: {bedge}\nblseg: {blseg}\n{card}')
        # raise RuntimeError(f'missing source contact surface {sid}\nbsurf: {bsurf}\nbsurfs: {bsurfs}\n{self}')
        sid_ref = None
    return sid_ref
