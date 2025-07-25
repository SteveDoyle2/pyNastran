"""
All nodes are defined in this file.  This includes:

 * Node
   * XPoint
     * EPOINT
     * SPOINT
   * XPoints
     * EPOINTs
     * SPOINTs
   * GRID
   * GRDSET
   * GRIDB
 * POINT
 * Ring
   * RINGAX
 * SEQGP

All ungrouped elements are Node objects.

The EPOINT/SPOINT classes refer to a single EPOINT/SPOINT.  The
EPOINTs/SPOINTs classes are for multiple degrees of freedom
(e.g. an SPOINT card).

"""
from __future__ import annotations
from itertools import count
from typing import Optional, Any, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf import MAX_INT
from pyNastran.bdf.field_writer_8 import set_string8_blank_if_default
from pyNastran.bdf.field_writer_16 import set_string16_blank_if_default

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru, write_card
from pyNastran.bdf.cards.collpase_card import collapse_thru_packs
from pyNastran.bdf.bdf_interface.internal_get import coord_id
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, blank, integer_or_string,
    integer_or_double, components_or_blank, parse_components_or_blank)
from pyNastran.bdf.bdf_interface.assign_type_force import (
    force_double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8, print_int_card
from pyNastran.bdf.field_writer_16 import print_float_16, print_card_16
from pyNastran.bdf.field_writer_double import print_scientific_double, print_card_double

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.cards.coordinate_systems import Coord
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.nptyping_interface import NDArray3float
    from pyNastran.bdf.cards.base_card import Element
    #from pyNastran.bdf.bdf_interface.typing import Coord, Element


class SEQGP(BaseCard):
    """defines the SEQGP class"""
    type = 'SEQGP'

    @classmethod
    def _init_from_empty(cls):
        nids = 1
        seqids = [2, 3]
        return SEQGP(nids, seqids, comment='')

    def __init__(self, nids, seqids, comment=''):
        """
        Creates the SEQGP card

        Parameters
        ----------
        nid : int
           the node id
        seqid : int/float
           the superelement id
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        if isinstance(nids, integer_types):
            nids = [nids]
        if isinstance(seqids, integer_types):
            seqids = [seqids]

        self.nids = nids
        self.seqids = seqids

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SEQGP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        ncard = len(card) - 1
        assert len(card) > 1, 'len(SEQGP) = 1; card=%s' % card
        assert ncard % 2 == 0, card
        nids = []
        seqids = []
        for ifield in range(1, ncard, 2):
            nid = integer(card, ifield, 'nid')
            seqid = integer_or_double(card, ifield+1, 'seqid')
            nids.append(nid)
            seqids.append(seqid)
        return SEQGP(nids, seqids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def append(self, seqgp):
        self.nids += seqgp.nids
        self.seqids += seqgp.seqids

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a SEQGP card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        nids, seqids = data
        return SEQGP(nids, seqids, comment=comment)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['SEQGP']
        for nid, seqid in zip(self.nids, self.seqids):
            list_fields.append(nid)
            list_fields.append(seqid)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card

        Parameters
        ----------
        size : int; default=8
            unused
        is_double : bool; default=False
            unused

        """
        _print_card = print_card_8
        if max(self.nids) > MAX_INT:
            _print_card = print_card_16

        msg = self.comment
        list_fields = ['SEQGP']
        for i, nid, seqid in zip(count(), self.nids, self.seqids):
            if i % 4 == 0 and i > 0:
                msg += _print_card(list_fields)
                list_fields = ['SEQGP']
            list_fields.append(nid)
            list_fields.append(seqid)

        if len(list_fields) > 1:
            msg += _print_card(list_fields)
        return msg


class XPoint(BaseCard):
    """common class for EPOINT/SPOINT"""

    @classmethod
    def _init_from_empty(cls):
        nid = 1
        return cls(nid, comment='')

    def __init__(self, nid, comment):
        #Node.__init__(self)
        if comment:
            self.comment = comment
        self.nid = nid
        assert isinstance(nid, integer_types), nid

    @classmethod
    def _export_to_hdf5(cls, h5_file, model: BDF, nids: list[int]) -> None:
        """exports the nodes in a vectorized way"""
        #comments = []
        #for nid in nids:
            #node = model.nodes[nid]
            #comments.append(element.comment)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('nid', data=nids)

    @property
    def type(self):
        """dummy method for EPOINT/SPOINT classes"""
        raise NotImplementedError('This method should be overwritten by the parent class')

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        lists_fields = []
        if isinstance(self.nid, integer_types):
            list_fields = [self.type, self.nid]
            lists_fields.append(list_fields)
        else:
            singles, doubles = collapse_thru_packs(self.nid)
            if singles:
                list_fields = [self.type] + singles
            if doubles:
                for spoint_double in doubles:
                    list_fields = [self.type] + spoint_double
                    lists_fields.append(list_fields)
        return lists_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card

        Parameters
        ----------
        size : int; default=8
            unused
        is_double : bool; default=False
            unused

        """
        msg = self.comment
        lists_fields = self.repr_fields()
        for list_fields in lists_fields:
            if 'THRU' not in list_fields:
                msg += print_int_card(list_fields)
            else:
                msg += print_card_8(list_fields)
        return msg

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        pass


class SPOINT(XPoint):
    """defines the SPOINT class"""
    type = 'SPOINT'

    def __init__(self, nid, comment=''):
        """
        Creates the SPOINT card

        Parameters
        ----------
        nid : int
           the SPOINT id
        comment : str; default=''
            a comment for the card

        """
        XPoint.__init__(self, nid, comment)

    def get_position(self):
        return np.zeros(3)


class EPOINT(XPoint):
    """defines the EPOINT class"""
    type = 'EPOINT'

    def __init__(self, nid, comment=''):
        """
        Creates the EPOINT card

        Parameters
        ----------
        nid : int
           the EPOINT id
        comment : str; default=''
            a comment for the card

        """
        XPoint.__init__(self, nid, comment)


def write_xpoints(cardtype: str, points, comment: str='') -> str:
    """writes SPOINTs/EPOINTs"""
    msg = comment
    if isinstance(points, dict):
        point_ids = []
        for point_id, point in sorted(points.items()):
            point_ids.append(point_id)
            if point.comment:
                msg += point.comment
    else:
        point_ids = points
    lists_fields = compress_xpoints(cardtype, point_ids)
    for list_fields in lists_fields:
        if 'THRU' not in list_fields:
            msg += print_int_card(list_fields)
        else:
            if max(point_ids) > MAX_INT:
                msg += print_card_16(list_fields)
            else:
                msg += print_card_8(list_fields)
    return msg


def compress_xpoints(point_type: str, xpoints: list[int]) -> list[list[int]]:
    """
    Gets the SPOINTs/EPOINTs in sorted, short form.

      uncompressed: SPOINT,1,3,5
      compressed:   SPOINT,1,3,5

      uncompressed: SPOINT,1,2,3,4,5
      compressed:   SPOINT,1,THRU,5

      uncompressed: SPOINT,1,2,3,4,5,7
      compressed:   SPOINT,7
                    SPOINT,1,THRU,5

    point_type = 'SPOINT'
    spoints = [1, 2, 3, 4, 5]
    fields = compressed_xpoints(point_type, spoints)
    >>> fields
    ['SPOINT', 1, 'THRU', 5]

    """
    spoints = list(xpoints)
    spoints.sort()

    singles, doubles = collapse_thru_packs(spoints)

    lists_fields = []
    if singles:
        list_fields = [point_type] + singles
        lists_fields.append(list_fields)
    if doubles:
        for spoint_double in doubles:
            list_fields = [point_type] + spoint_double
            lists_fields.append(list_fields)
    return lists_fields


class XPoints(BaseCard):
    """common class for EPOINTs and SPOINTs"""

    @property
    def type(self):
        """dummy method for EPOINTs/SPOINTs classes"""
        raise NotImplementedError('This method should be overwritten by the parent class')

    @classmethod
    def _init_from_empty(cls):
        ids = [1]
        return cls(ids, comment='')

    def __init__(self, ids, comment=''):
        #Node.__init__(self)
        if comment:
            self.comment = comment
        if isinstance(ids, integer_types):
            ids = [ids]
        self.points = set(expand_thru(ids))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPOINT/EPOINT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        points = []
        for i in range(1, len(card)):
            field = integer_or_string(card, i, 'ID%i' % i)
            points.append(field)
        return cls(points, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a SPOINT/EPOINT card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        points = data
        assert isinstance(points, list), points
        assert isinstance(points[0], integer_types), points
        assert min(points) > 0, points
        return cls(points, comment=comment)

    def __len__(self):
        """
        Returns the number of degrees of freedom for the EPOINTs/SPOINTs class

        Returns
        -------
        ndofs : int
            the number of degrees of freedom

        """
        return len(self.points)

    def add_points(self, slist):
        """Adds more EPOINTs/SPOINTs to this object"""
        self.points = self.points.union(set(slist))

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        points = list(self.points)
        points.sort()
        return [self.type] + points

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card

        Parameters
        ----------
        size : int; default=8
            unused
        is_double : bool; default=False
            unused

        """
        lists_fields = compress_xpoints(self.type, self.points)
        msg = self.comment
        for list_fields in lists_fields:
            if 'THRU' not in list_fields:
                msg += print_int_card(list_fields)
            else:
                msg += print_card_8(list_fields)
        return msg


class SPOINTs(XPoints):
    """
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |   1    |  2  |  3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+=====+======+=====+=====+=====+=====+=====+=====+
    | SPOINT | ID1 | THRU | ID2 |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    | SPOINT | ID1 | ID1  | ID3 | ID4 | ID5 | ID6 | ID7 | ID8 |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |        | ID8 | etc. |     |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+

    """
    type = 'SPOINT'

    def __init__(self, ids, comment=''):
        """
        Creates the SPOINTs card that contains many SPOINTs

        Parameters
        ----------
        ids : list[int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        XPoints.__init__(self, ids, comment=comment)

    def create_spointi(self):
        """Creates individal SPOINT objects"""
        spoints = []
        for nid in self.points:
            spoint = SPOINT(nid)
            spoints.append(spoint)
        if hasattr(self, 'ifile'):
            for spoint in spoints:
                spoint.ifile = self.ifile
        return spoints


class EPOINTs(XPoints):
    """
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |   1    |  2  |  3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+=====+======+=====+=====+=====+=====+=====+=====+
    | EPOINT | ID1 | THRU | ID2 |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    | EPOINT | ID1 | ID1  | ID3 | ID4 | ID5 | ID6 | ID7 | ID8 |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |        | ID8 | etc. |     |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+

    """
    type = 'EPOINT'

    def __init__(self, ids, comment=''):
        """
        Creates the EPOINTs card that contains many EPOINTs

        Parameters
        ----------
        ids : list[int]
            EPOINT ids
        comment : str; default=''
            a comment for the card

        """
        XPoints.__init__(self, ids, comment=comment)

    def create_epointi(self):
        """Creates individal EPOINT objects"""
        points = []
        for nid in self.points:
            points.append(EPOINT(nid))
        return points


class GRDSET(BaseCard):
    """
    Defines default options for fields 3, 7, 8, and 9 of all GRID entries.

    +--------+-----+----+----+----+----+----+----+------+
    |    1   |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
    +========+=====+====+====+====+====+====+====+======+
    | GRDSET |     | CP |    |    |    | CD | PS | SEID |
    +--------+-----+----+----+----+----+----+----+------+

    """
    type = 'GRDSET'

    #: allows the get_field method and update_field methods to be used
    _field_map = {2: 'cp', 6: 'cd', 7: 'ps', 8: 'seid'}

    @classmethod
    def _init_from_empty(cls):
        cp = 0
        cd = 1
        ps = '34'
        seid = 0
        return GRDSET(cp, cd, ps, seid, comment='')

    def __init__(self, cp: int, cd: int, ps: str, seid: int, comment: str=''):
        """
        Creates the GRDSET card

        Parameters
        ----------
        cp : int; default=0
            the xyz coordinate frame
        cd : int; default=0
            the analysis coordinate frame
        ps : str; default=''
            Additional SPCs in the analysis coordinate frame (e.g. '123').
            This corresponds to DOF set ``SG``.
        seid : int; default=0
            superelement id
            TODO: how is this used by Nastran???
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: Output Coordinate System
        self.cp = cp

        #: Analysis coordinate system
        self.cd = cd

        #: Default SPC constraint on undefined nodes
        self.ps = ps

        #: Superelement ID
        self.seid = seid

        self.cp_ref = None
        self.cd_ref = None
        self.seid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a GRDSET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #: Grid point coordinate system
        blank(card, 1, 'blank')

        cp = integer_or_blank(card, 2, 'cp', 0)
        blank(card, 3, 'blank')
        blank(card, 4, 'blank')
        blank(card, 5, 'blank')
        cd = integer_or_blank(card, 6, 'cd', 0)

        ps = str(parse_components_or_blank(card, 7, 'ps', default=''))
        seid = integer_or_blank(card, 8, 'seid', 0)
        assert len(card) <= 9, f'len(GRDSET card) = {len(card):d}\ncard={card}'
        return GRDSET(cp, cd, ps, seid, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by the GRDSET'
        self.cp_ref = model.Coord(self.cp, msg=msg)
        self.cd_ref = model.Coord(self.cd, msg=msg)
        #self.seid = model.SuperElement(self.seid, msg)
        #self.seid_ref = self.seid

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cp = self.Cp()
        self.cd = self.Cd()
        self.cp_ref = None
        self.cd_ref = None
        self.seid_ref = None

    def Cd(self):
        """
        Gets the output coordinate system

        Returns
        -------
        cd : int
            the output coordinate system

        """
        return coord_id(self.cd_ref, self.cd)

    def Cp(self) -> int:
        """
        Gets the analysis coordinate system

        Returns
        -------
        cp : int
            the analysis coordinate system

        """
        return coord_id(self.cp_ref, self.cp)

    def Ps(self):
        """
        Gets the GRID-based SPC

        Returns
        -------
        ps : str
            the GRID-based SPC

        """
        return self.ps

    def SEid(self):
        """
        Gets the Superelement ID

        Returns
        -------
        seid : int
            the Superelement ID

        """
        if self.seid_ref is None:
            return self.seid
        return self.seid_ref.seid

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref: bool
            has this model been cross referenced

        """
        cp = self.Cp()
        seid = self.SEid()
        cd = self.Cd()
        ps = self.Ps()
        assert isinstance(cp, integer_types), 'cp=%r' % cp
        assert isinstance(cd, integer_types), 'cd=%r' % cd
        assert isinstance(ps, str), 'ps=%r' % ps
        assert isinstance(seid, integer_types), 'seid=%r' % seid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['GRDSET', None, self.Cp(), None, None, None,
                       self.Cd(), self.ps, self.SEid()]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        cp = set_blank_if_default(self.Cp(), 0)
        cd = set_blank_if_default(self.Cd(), 0)
        ps = set_blank_if_default(self.ps, 0)
        seid = set_blank_if_default(self.SEid(), 0)
        list_fields = ['GRDSET', None, cp, None, None, None, cd, ps, seid]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card

        Parameters
        ----------
        size : int
            the size of the card (8/16)

        """
        card = self.repr_fields()
        return print_card_8(card)


class GRID(BaseCard):
    """
    +------+-----+----+----+----+----+----+----+------+
    |   1  |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
    +======+=====+====+====+====+====+====+====+======+
    | GRID | NID | CP | X1 | X2 | X3 | CD | PS | SEID |
    +------+-----+----+----+----+----+----+----+------+

    Attributes
    ----------
    nid : int
        node id
    xyz : float ndarray
        Raw location <:math:`x_1, x_2, x_3`>
    cp : int
        reference coordinate system
    cd : int
        analysis coordinate system
    ps : str
        nodal-based constraints
    seid : int
        superelement id
    cp_ref : Coord() or None
        cross-referenced cp
    cd_ref : Coord() or None
        cross-referenced cd

    Methods
    -------
    Nid()
        gets nid
    Cp()
        gets cp_ref.cid or cp depending on cross-referencing
    Cd()
        gets cd_ref.cid or cd depending on cross-referencing
    Ps()
        gets ps
    SEid()
        superelement id
    get_position()
        gets xyz in the global frame
    get_position_wrt(model, cid)
        gets xyz in a local frame
    cross_reference(model)
        cross-references the card
    uncross_reference()
        uncross-references the card
    set_position(model, xyz, cid=0, xref=True)
        updates the coordinate system

    Using the GRID object::

     model = read_bdf(bdf_filename)
     node = model.Node(nid)

     # gets the position of the node in the global frame
     node.get_position()
     node.get_position_wrt(model, cid=0)

     # gets the position of the node in a local frame
     node.get_position_wrt(model, cid=1)

     # change the location of the node
     node.set_position(model, array([1.,2.,3.]), cid=3)

    """
    type = 'GRID'

    #: allows the get_field method and update_field methods to be used
    _field_map = {1: 'nid', 2: 'cp', 6: 'cd', 7: 'ps', 8: 'seid'}

    def _get_field_helper(self, n: int):
        """
        Gets complicated parameters on the GRID card

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : float
            the value for the appropriate field

        """
        if n == 3:
            value = self.xyz[0]
        elif n == 4:
            value = self.xyz[1]
        elif n == 5:
            value = self.xyz[2]
        else:
            raise KeyError('Field %r is an invalid %s entry.' % (n, self.type))
        return value

    def _update_field_helper(self, n: int, value: Any):
        """
        Updates complicated parameters on the GRID card

        Parameters
        ----------
        n : int
            the field number to update
        value : float
            the value for the appropriate field

        """
        if n in {3, 'x1'}:
            self.xyz[0] = value
        elif n in {4, 'x2'}:
            self.xyz[1] = value
        elif n in {5, 'x3'}:
            self.xyz[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    @classmethod
    def export_to_hdf5(cls, h5_file, model, nids):
        """exports the nodes in a vectorized way"""
        comments = []
        cp = []
        xyz = []
        cd = []
        ps = []
        seid = []
        for nid in nids:
            node = model.nodes[nid]
            #comments.append(element.comment)
            xyz.append(node.xyz)
            cp.append(node.cp)
            cd.append(node.cd)
            psi = 0 if node.ps == '' else (int(node.ps))
            ps.append(psi)
            seid.append(node.seid)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('nid', data=nids)
        h5_file.create_dataset('xyz', data=xyz)
        h5_file.create_dataset('cp', data=cp)
        h5_file.create_dataset('cd', data=cd)
        h5_file.create_dataset('ps', data=ps)
        h5_file.create_dataset('seid', data=seid)

    def __init__(self, nid: int, xyz: np.ndarray | list[float],
                 cp: int=0, cd: int=0, ps: str='', seid: int=0,
                 comment: str='') -> None:
        """
        Creates the GRID card

        Parameters
        ----------
        nid : int
            node id
        cp : int; default=0
            the xyz coordinate frame
        xyz : (3, ) float ndarray; default=None -> [0., 0., 0.]
            the xyz/r-theta-z/rho-theta-phi values
        cd : int; default=0
            the analysis coordinate frame
        ps : str; default=''
            Additional SPCs in the analysis coordinate frame (e.g. '123').
            This corresponds to DOF set ``SG``.
        seid : int; default=0
            superelement id
            TODO: how is this used by Nastran???
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.nid = nid
        self.cp = cp
        if xyz is None:
            xyz = [0., 0., 0.]
        self.xyz = np.asarray(xyz, dtype='float64')
        assert self.xyz.size == 3, self.xyz.shape
        self.cd = cd
        self.ps = ps
        self.seid = seid
        self.cp_ref: Coord = None
        self.cd_ref: Coord = None
        self.elements_ref: list[Element] = None

    @classmethod
    def add_op2_data(cls, data: list[int | float], comment: str='') -> GRID:
        """
        Adds a GRID card from the OP2.

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        .. todo:: Currently unused, but is tested by test_nodes.py

        """
        nid = data[0]
        cp = data[1]
        xyz = data[2:5]
        cd = data[5]
        ps = data[6]
        seid = data[7]
        if ps == 0:
            ps = ''
        return GRID(nid, xyz, cp, cd, ps, seid, comment=comment)

    @classmethod
    def add_card(cls, card: BDFCard, comment: str='') -> GRID:
        """
        Adds a GRID card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nfields = len(card)
        #: Node ID
        nid = integer(card, 1, 'nid')

        #: Grid point coordinate system
        cp = integer_or_blank(card, 2, 'cp', default=0)

        #: node location in local frame
        xyz = [
            double_or_blank(card, 3, 'x1', default=0.),
            double_or_blank(card, 4, 'x2', default=0.),
            double_or_blank(card, 5, 'x3', default=0.)]

        if nfields > 6:
            #: Analysis coordinate system
            cd = integer_or_blank(card, 6, 'cd', default=0)

            #: SPC constraint
            ps = components_or_blank(card, 7, 'ps', default='')
            #u(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            seid = integer_or_blank(card, 8, 'seid', default=0)
            assert len(card) <= 9, f'len(GRID card) = {len(card):d}\ncard={card}'
        else:
            cd = 0
            ps = ''
            seid = 0
        return GRID(nid, xyz, cp, cd, ps, seid, comment=comment)

    @classmethod
    def add_card_lax(cls, card: BDFCard, comment: str='') -> GRID:
        """
        Adds a GRID card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nfields = len(card)
        #: Node ID
        nid = integer(card, 1, 'nid')

        #: Grid point coordinate system
        cp = integer_or_blank(card, 2, 'cp', default=0)

        #: node location in local frame
        xyz = [
            force_double_or_blank(card, 3, 'x1', default=0.0),
            force_double_or_blank(card, 4, 'x2', default=0.0),
            force_double_or_blank(card, 5, 'x3', default=0.0)]

        if nfields > 6:
            #: Analysis coordinate system
            cd = integer_or_blank(card, 6, 'cd', default=0)

            #: SPC constraint
            ps = components_or_blank(card, 7, 'ps', default='')
            #u(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            seid = integer_or_blank(card, 8, 'seid', default=0)
            assert len(card) <= 9, f'len(GRID card) = {len(card):d}\ncard={card}'
        else:
            cd = 0
            ps = ''
            seid = 0
        return GRID(nid, xyz, cp, cd, ps, seid, comment=comment)

    def validate(self) -> None:
        assert isinstance(self.cp, integer_types), f'cp={self.cp:d}'
        assert self.nid > 0, f'nid={self.nid:d}'
        assert self.cp >= 0, f'cp={self.cp:d}'
        assert self.cd >= -1, f'cd={self.cd:d}'
        assert self.seid >= 0, f'seid={self.seid:d}'
        assert len(self.xyz) == 3

    def Nid(self) -> int:
        """
        Gets the GRID ID

        Returns
        -------
        nid : int
            node ID

        """
        return self.nid

    def Ps(self) -> str:
        """
        Gets the GRID-based SPC

        Returns
        -------
        ps : str
            the GRID-based SPC

        """
        return self.ps

    def Cd(self) -> int:
        """
        Gets the output coordinate system

        Returns
        -------
        cd : int
            the output coordinate system

        """
        if self.cd_ref is None:
            return self.cd
        return self.cd_ref.cid

    def Cp(self) -> int:
        """
        Gets the analysis coordinate system

        Returns
        -------
        cp : int
            the analysis coordinate system

        """
        return coord_id(self.cp_ref, self.cp)

    def SEid(self) -> int:
        """
        Gets the Superelement ID

        Returns
        -------
        seid : int
            the Superelement ID

        """
        #if isinstance(self.seid, integer_types):
        return self.seid
        #return self.seid.seid

    def _verify(self, xref: bool) -> None:
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced

        """
        nid = self.Nid()
        cp = self.Cp()
        cd = self.Cd()
        xyz = self.xyz
        ps = self.Ps()
        seid = self.SEid()
        assert isinstance(xyz, np.ndarray), 'xyz=%r' % xyz
        assert isinstance(nid, integer_types), 'nid=%r' % nid
        assert isinstance(cp, integer_types), 'cp=%r' % cp
        assert isinstance(cd, integer_types), 'cd=%r' % cd
        assert isinstance(ps, str), 'ps=%r' % ps
        assert isinstance(seid, integer_types), 'seid=%r' % seid
        if xref:
            pos_xyz = self.get_position()
            assert isinstance(pos_xyz, np.ndarray), 'pos_xyz=%r' % pos_xyz

    def set_position(self, model: BDF, xyz: np.ndarray,
                     cid: int=0, xref: bool=True) -> None:
        # (Any, np.ndarray, int) -> None
        """
        Updates the GRID location

        Parameters
        ----------
        xyz : (3, ) float ndarray
            the location of the node.
        cp : int; default=0 (global)
            the analysis coordinate system
        xref : bool; default=True
            cross-references the coordinate system

        """
        self.xyz = xyz
        msg = ', which is required by GRID nid=%s' % self.nid
        self.cp = cid
        if xref:
            self.cp_ref = model.Coord(cid, msg=msg)

    def get_position_no_xref(self, model: Any) -> np.ndarray:
        # (Any) -> np.ndarray
        if self.cp == 0:
            return self.xyz
        assert isinstance(self.cp, integer_types), self.cp
        coord = model.Coord(self.cp)
        xyz = coord.transform_node_to_global_no_xref(self.xyz, model)
        return xyz

    def get_position(self) -> np.ndarray:
        """
        Gets the point in the global XYZ coordinate system.

        Returns
        -------
        xyz : (3, ) float ndarray
            the position of the GRID in the global coordinate system

        """
        try:
            xyz = self.cp_ref.transform_node_to_global(self.xyz)
        except AttributeError:
            if self.cp == 0:
                return self.xyz
            raise
        return xyz

    def get_position_assuming_rectangular(self) -> NDArray3float:
        """
        Gets the point in a coordinate system that has unit vectors
        in the referenced coordinate system, but is not transformed
        from a cylindrical/spherical system.  This is used by cards
        like CBAR/CBEAM for element offset vectors.

        Returns
        -------
        xyz : (3, ) float ndarray
            the position of the GRID in the global coordinate system

        """
        try:
            xyz = self.cp_ref.transform_node_to_global(self.xyz)
        except AttributeError:
            if self.cp == 0:
                return self.xyz
            raise
        return xyz

    def get_position_wrt_no_xref(self, model: BDF, cid: int) -> NDArray3float:
        """see get_position_wrt"""
        if cid == self.cp:  # same coordinate system
            return self.xyz
        msg = ', which is required by GRID nid=%s' % self.nid

        # converting the xyz point arbitrary->global
        cp_ref = model.Coord(self.cp, msg=msg)
        p = cp_ref.transform_node_to_global_no_xref(self.xyz, model)

        # a matrix global->local matrix is found
        coord_b = model.Coord(cid, msg=msg)
        xyz = coord_b.transform_node_to_local(p)
        return xyz

    def get_position_wrt(self, model: BDF, cid: int) -> np.ndarray:
        """
        Gets the location of the GRID which started in some arbitrary
        system and returns it in the desired coordinate system

        Parameters
        ----------
        model : BDF()
            the BDF object
        cid : int
            the desired coordinate ID

        Returns
        -------
        xyz : (3, ) float ndarray
            the position of the GRID in an arbitrary coordinate system

        """
        if cid == self.Cp():  # same coordinate system
            return self.xyz
        # a matrix global->local matrix is found
        msg = f', which is required by GRID nid={self.nid:d}'
        coord_b: Coord = model.Coord(cid, msg=msg)
        return self.get_position_wrt_coord_ref(coord_b)

    def get_position_wrt_coord_ref(self, coord_out: Coord) -> np.ndarray:
        """
        Gets the location of the GRID which started in some arbitrary
        system and returns it in the desired coordinate system

        Parameters
        ----------
        coord_out : Coord
            the desired coordinate system

        Returns
        -------
        xyz : (3, ) float ndarray
            the position of the GRID in an arbitrary coordinate system

        """
        # converting the xyz point arbitrary->global
        p = self.cp_ref.transform_node_to_global(self.xyz)

        xyz = coord_out.transform_node_to_local(p)
        return xyz

    def cross_reference(self, model: BDF, grdset: Optional[Any]=None) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        grdset : GRDSET / None; default=None
            a GRDSET if available (default=None)

        .. note::  The gridset object will only update the fields that
                   have not been set

        """
        if grdset:
            # update using a grdset object
            if not self.cp:
                self.cp_ref = grdset.cp_ref
            if not self.cd:
                self.cd = grdset.cd
                self.cd_ref = self.cd_ref
            if not self.ps:
                self.ps_ref = grdset.ps
            if not self.seid:
                self.seid_ref = grdset.seid
        msg = ', which is required by GRID nid=%s' % (self.nid)
        self.cp_ref = model.Coord(self.cp, msg=msg)
        if self.cd != -1:
            self.cd_ref = model.Coord(self.cd, msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors: dict[str, Any],
                             grdset: Optional[Any]=None) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        grdset : GRDSET / None; default=None
            a GRDSET if available (default=None)

        .. note::  The gridset object will only update the fields that
                   have not been set

        """
        if grdset:
            # update using a grdset object
            if not self.cp:
                self.cp_ref = grdset.cp_ref
            if not self.cd:
                self.cd = grdset.cd
                self.cd_ref = self.cd_ref
            if not self.ps:
                self.ps_ref = grdset.ps
            if not self.seid:
                self.seid_ref = grdset.seid
        msg = f', which is required by GRID nid={self.nid:d}'
        self.cp_ref = model.safe_coord(self.cp, self.nid, xref_errors, msg=msg)
        if self.cd != -1:
            self.cd_ref = model.safe_coord(self.cd, self.nid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cd_ref = None
        self.cp_ref = None
        self.elements_ref = None

    def raw_fields(self) -> list[Any]:
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[int/float/str]
            the fields that define the card

        """
        list_fields = ['GRID', self.nid, self.Cp()] + list(self.xyz) + \
                      [self.Cd(), self.ps, self.SEid()]
        return list_fields

    def repr_fields(self) -> list[Any]:
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[int/float/str]
            the fields that define the card

        """
        cp = set_blank_if_default(self.Cp(), 0)
        cd = set_blank_if_default(self.Cd(), 0)
        seid = set_blank_if_default(self.SEid(), 0)
        list_fields = ['GRID', self.nid, cp] + list(self.xyz) + [cd, self.ps,
                                                                 seid]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card

        Parameters
        ----------
        size : int; default=8
            the size of the card (8/16)
        is_double : bool; default=False
            should this card be written with double precision

        Returns
        -------
        msg : str
            the card as a string

        """
        if size == 8:
            return self.write_card_8()
        return self.write_card_16(is_double)

    def write_card_8(self) -> str:
        """Writes a GRID card in 8-field format"""
        xyz = self.xyz
        cp = self.Cp()
        cd = self.Cd()

        cps = set_string8_blank_if_default(cp, 0)
        if [cd, self.ps, self.seid] == [0, '', 0]:
            # default
            msg = 'GRID    %8i%8s%s%s%s\n' % (
                self.nid, cps,
                print_float_8(xyz[0]),
                print_float_8(xyz[1]),
                print_float_8(xyz[2]))
            #msg = 'GRID,%i,%s,%s,%s,%s\n' % (
                #self.nid, cps.strip(),
                #print_float_8(xyz[0]),
                #print_float_8(xyz[1]),
                #print_float_8(xyz[2]))
        else:
            cds = set_string8_blank_if_default(cd, 0)
            seid = set_string8_blank_if_default(self.SEid(), 0)
            msg = 'GRID    %8i%8s%s%s%s%s%8s%s\n' % (
                self.nid, cps,
                print_float_8(xyz[0]),
                print_float_8(xyz[1]),
                print_float_8(xyz[2]),
                cds, self.ps, seid)
        return self.comment + msg

    def write_card_16(self, is_double: bool=False) -> str:
        """Writes a GRID card in 16-field format"""
        xyz = self.xyz
        cp = set_string16_blank_if_default(self.Cp(), 0)
        cd = set_string16_blank_if_default(self.Cd(), 0)
        seid = set_string16_blank_if_default(self.SEid(), 0)

        if is_double:
            if [cd, self.ps, self.seid] == [0, '', 0]:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s\n' % (
                           self.nid,
                           cp,
                           print_scientific_double(xyz[0]),
                           print_scientific_double(xyz[1]),
                           print_scientific_double(xyz[2])))
            else:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s%16s%16s%16s\n' % (
                           self.nid,
                           cp,
                           print_scientific_double(xyz[0]),
                           print_scientific_double(xyz[1]),
                           print_scientific_double(xyz[2]),
                           cd, self.ps, seid))
        else:
            if [cd, self.ps, self.seid] == [0, '', 0]:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s\n' % (
                           self.nid,
                           cp,
                           print_float_16(xyz[0]),
                           print_float_16(xyz[1]),
                           print_float_16(xyz[2])))
            else:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s%16s%16s%16s\n' % (
                           self.nid,
                           cp,
                           print_float_16(xyz[0]),
                           print_float_16(xyz[1]),
                           print_float_16(xyz[2]),
                           cd, self.ps, seid))
        return self.comment + msg


class POINT(BaseCard):
    """
    +-------+-----+----+----+----+----+
    |   1   |  2  | 3  | 4  | 5  | 6  |
    +=======+=====+====+====+====+====+
    | POINT | NID | CP | X1 | X2 | X3 |
    +-------+-----+----+----+----+----+

    """
    type = 'POINT'
    _field_map = {1: 'nid', 2: 'cp'}

    def _get_field_helper(self, n: int) -> float:
        """
        Gets complicated parameters on the POINT card

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : varies
            the value for the appropriate field

        """
        if n == 3:
            value = self.xyz[0]
        elif n == 4:
            value = self.xyz[1]
        elif n == 5:
            value = self.xyz[2]
        else:
            raise KeyError('Field %r is an invalid %s entry.' % (n, self.type))
        return value

    def _update_field_helper(self, n: int, value: float) -> None:
        """
        Updates complicated parameters on the POINT card

        Parameters
        ----------
        n : int
            the field number to update
        value : varies
            the value for the appropriate field

        """
        if n == 3:
            self.xyz[0] = value
        elif n == 4:
            self.xyz[1] = value
        elif n == 5:
            self.xyz[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    @classmethod
    def _init_from_empty(cls):
        nid = 1
        xyz = [1., 2., 3.]
        return POINT(nid, xyz, cp=0, comment='')

    def __init__(self, nid: int,
                 xyz: np.ndarray | list[float],
                 cp: int=0, comment: str='') -> None:
        """
        Creates the POINT card

        Parameters
        ----------
        nid : int
            node id
        xyz : (3, ) float ndarray; default=None -> [0., 0., 0.]
            the xyz/r-theta-z/rho-theta-phi values
        cp : int; default=0
            coordinate system for the xyz location
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        if xyz is None:
            xyz = [0., 0., 0.]
        #Node.__init__(self)

        #: Node ID
        self.nid = nid

        #: Grid point coordinate system
        self.cp = cp

        #: node location in local frame
        self.xyz = np.asarray(xyz, dtype='float64')
        assert self.xyz.size == 3, self.xyz.shape
        self.cp_ref = None

    def validate(self) -> None:
        assert self.nid > 0, 'nid=%s' % self.nid
        assert self.cp >= 0, 'cp=%s' % self.cp
        assert len(self.xyz) == 3

    @classmethod
    def add_card(cls, card: BDFCard, comment: str='') -> POINT:
        """
        Adds a POINT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nid = integer(card, 1, 'nid')
        cp = integer_or_blank(card, 2, 'cp', default=0)

        xyz = np.array([
            double_or_blank(card, 3, 'x1', default=0.),
            double_or_blank(card, 4, 'x2', default=0.),
            double_or_blank(card, 5, 'x3', default=0.)], dtype='float64')

        assert len(card) <= 9, f'len(POINT card) = {len(card):d}\ncard={card}'
        return POINT(nid, xyz, cp=cp, comment=comment)

    @classmethod
    def add_op2_data(cls, data: list[int | float], comment: str='') -> POINT:
        """
        Adds a POINT card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        nid: int = data[0]
        cp: int = data[1]
        xyz: np.ndarray = np.array(data[2:5])
        return POINT(nid, xyz, cp=cp, comment=comment)

    def set_position(self, model: BDF, xyz: np.ndarray, cid: int=0) -> None:
        """
        Updates the POINT location

        Parameters
        ----------
        xyz : (3,) float ndarray
            the location of the node
        cp : int; default=0 (global)
            the analysis coordinate system

        """
        self.xyz = xyz
        msg = ', which is required by POINT nid=%s' % self.nid
        self.cp = model.Coord(cid, msg=msg)

    def get_position(self) -> NDArray3float:
        """
        Gets the point in the global XYZ coordinate system.

        Returns
        -------
        position : (3,) float ndarray
            the position of the POINT in the globaly coordinate system

        """
        p = self.cp_ref.transform_node_to_global(self.xyz)
        return p

    def get_position_wrt(self, model: BDF, cid: int) -> np.ndarray:
        """
        Gets the location of the POINT which started in some arbitrary
        system and returns it in the desired coordinate system

        Parameters
        ----------
        model : BDF()
            the BDF model object
        cid : int
            the desired coordinate ID

        Returns
        -------
        xyz : (3,) ndarray
            the position of the POINT in an arbitrary coordinate system

        """
        if cid == self.Cp():  # same coordinate system
            return self.xyz

        # converting the xyz point arbitrary->global
        p = self.cp_ref.transform_node_to_global(self.xyz)

        # a matrix global->local matrix is found
        msg = ', which is required by POINT nid=%s' % (self.nid)
        coord_b = model.Coord(cid, msg=msg)
        xyz = coord_b.transform_node_to_local(p)
        return xyz

    def Cp(self) -> int:
        """
        Gets the analysis coordinate system

        Returns
        -------
        cp : int
            the analysis coordinate system

        """
        return coord_id(self.cp_ref, self.cp)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        self.cp_ref = model.Coord(self.cp)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cp_ref = self.Cp()

    def raw_fields(self) -> list[str | int | float | None]:
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['POINT', self.nid, self.Cp()] + list(self.xyz)
        return list_fields

    def repr_fields(self) -> list[str | int | float]:
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        cp = set_blank_if_default(self.Cp(), 0)
        list_fields = ['POINT', self.nid, cp] + list(self.xyz)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card

        Parameters
        ----------
        size : int
            the size of the card (8/16)

        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
