"""
All superelements are defined in this file.  This includes:
 * CSUPER
 * CSUPEXT
 * SEBNDRY
 * SEBULK
 * SECONCT
 * SEELT
 * SEEXCLD
 * SELABEL
 * SELOC
 * SELOAD
 * SEMPLN
 * SENQSET
 * SETREE

"""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import (
    BaseCard, expand_thru #, _node_ids
)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_or_string,
    string, string_or_blank, double_or_blank, integer_string_or_blank,
    exact_string_or_blank,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class SEBNDRY(BaseCard):
    """
    Superelement Boundary-Point Definition

    Defines a list of grid points in a partitioned superelement for the
    automatic boundary search between a specified superelement or between all
    other superelements in the model.

    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    |    1    |   2   |   3   |   4   |   5   |   6   |   7   |   8   |   9   |
    +=========+=======+=======+=======+=======+=======+=======+=======+=======+
    | SEBNDRY | SEIDA | SEIDB | GIDA1 | GIDA2 | GIDA3 | GIDA4 | GIDA5 | GIDA6 |
    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    |         | GIDA7 | GIDA8 |  etc. |       |       |       |       |       |
    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    | SEBNDRY |  400  |   4   |   10  |   20  |   30  |   40  |       |       |
    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    | SEBNDRY |  400  |  ALL  |   10  |   20  |   30  |  THRU |   40  |       |
    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    """
    type = 'SEBNDRY'
    @classmethod
    def _init_from_empty(cls):
        seid_a = 1
        seid_b = 2
        ids = [10, 20, 30]
        return SEBNDRY(seid_a, seid_b, ids, comment='')

    def __init__(self, seid_a, seid_b, ids, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.seid_a = seid_a
        self.seid_b = seid_b

        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)
        self.ids_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid_a = integer(card, 1, 'seid_a')
        seid_b = integer_string_or_blank(card, 2, 'seid_b')
        ids = []
        i = 1
        nfields = len(card)
        for ifield in range(3, nfields):
            idi = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        assert len(card) >= 3, 'len(SEBNDRY card) = %i\ncard=%s' % (len(card), card)
        return SEBNDRY(seid_a, seid_b, ids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        list_fields = ['SEBNDRY', self.seid_a, self.seid_b] + self.ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class RELEASE(BaseCard):
    """
    Superelement Boundary Grid Point Release

    Defines degrees-of-freedom for superelement exterior grid points
    that are not connected to the superelement.

    +---------+------+------+------+------+------+------+------+------+
    |    1    |  2   |  3   |  4   |  5   |  6   |  7   |  8   |  9   |
    +=========+======+======+======+======+======+======+======+======+
    | RELEASE | SEID | COMP | GID1 | GID2 | GID3 | GID4 | GID5 | GID6 |
    +---------+------+------+------+------+------+------+------+------+
    |         | GID7 | GID8 | etc. |      |      |      |      |      |
    +---------+------+------+------+------+------+------+------+------+
    | RELEASE | 400  |  4   |  10  |  20  |  30  |  40  |      |      |
    +---------+------+------+------+------+------+------+------+------+
    | RELEASE | 400  | 156  |  30  | THRU |  40  |      |      |      |
    +---------+------+------+------+------+------+------+------+------+
    | RELEASE | 400  | 156  | ALL  |      |      |      |      |      |
    +---------+------+------+------+------+------+------+------+------+
    """
    type = 'RELEASE'
    @classmethod
    def _init_from_empty(cls):
        seid = 1
        comp = 1
        nids = [10, 20, 30]
        return SEBNDRY(seid, comp, nids, comment='')

    def __init__(self, seid, comp, nids, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.seid = seid
        self.comp = comp

        #:  Identifiers of grids points. (Integer > 0)
        self.nids = expand_thru(nids)
        self.nids_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        comp = integer(card, 2, 'comp')
        nids = []
        i = 3
        nfields = len(card)
        for ifield in range(3, nfields):
            idi = integer_or_string(card, ifield, 'ID%i' % i)
            nids.append(idi)
            i += 1
        assert len(card) >= 3, 'len(RELEASE card) = %i\ncard=%s' % (len(card), card)
        return RELEASE(seid, comp, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        list_fields = ['RELEASE', self.seid, self.comp] + self.nids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SEELT(BaseCard):
    """
    +-------+------+------+------+------+------+------+------+------+
    |   1   |   2  |   3  |   4  |   5  |   6  |   7  |   8  |   9  |
    +=======+======+======+======+======+======+======+======+======+
    | SEELT | SEID | EID1 | EID2 | EID3 | EID4 | EID5 | EID6 | EID7 |
    +-------+------+------+------+------+------+------+------+------+
    """
    type = 'SEELT'

    @classmethod
    def _init_from_empty(cls):
        seid = 10
        eids = [1, 2, 3]
        return SEELT(seid, eids, comment='')

    def __init__(self, seid, eids, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.seid = seid
        #:  Identifiers of grids points. (Integer > 0)
        self.eids = expand_thru(eids)
        self.eids_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        eids = []
        i = 1
        nfields = len(card)
        for ifield in range(2, nfields):
            eid = integer_string_or_blank(card, ifield, 'eid_%i' % i)
            if eid:
                i += 1
                eids.append(eid)
        assert len(card) <= 9, 'len(SEELT card) = %i\ncard=%s' % (len(card), card)
        return SEELT(seid, eids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SEELT seid=%s' % (self.seid)
        eids_ref = self._xref_elements_plotels(model, self.eids, msg=msg)
        self.eids_ref = eids_ref

    def _xref_elements_plotels(self, model, eids, msg=''):
        eids_ref = []
        missing_eids = []
        for eid in eids:
            if eid in model.elements:
                elem = model.elements[eid]
            elif eid in model.plotels:
                elem = model.plotels[eid]
            else:
                missing_eids.append(eid)
                continue
            eids_ref.append(elem)
        if missing_eids:
            raise KeyError('eids=%s not found%s' % (missing_eids, msg))
        return eids_ref

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        return self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.eids_ref = None

    def raw_fields(self):
        list_fields = ['SEELT', self.seid] + self.eids  ## TODO: xref
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SELOAD(BaseCard):
    """
    External Superelement Load Mapping to Residual

    Maps loads from an external superelement to a specified load set for
    the residual structure.

    +--------+-------+------+-------+
    |    1   |   2   |   3  |   4   |
    +========+=======+======+=======+
    | SELOAD | LIDS0 | SEID | LIDSE |
    +--------+-------+------+-------+
    | SELOAD | 10010 | 100  |  10   |
    +--------+-------+------+-------+
    """
    type = 'SELOC'

    @classmethod
    def _init_from_empty(cls):
        lid_s0 = 1
        seid = 2
        lid_se = 3
        return SELOAD(lid_s0, seid, lid_se, comment='')

    def __init__(self, lid_s0, seid, lid_se, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.lid_s0 = lid_s0
        self.seid = seid
        self.lid_se = lid_se

    @classmethod
    def add_card(cls, card, comment=''):
        lid_s0 = integer(card, 1, 'lid_s0')
        seid = integer(card, 2, 'seid')
        lid_se = integer(card, 3, 'lid_se')

        assert len(card) <= 4, 'len(SELOAD card) = %i\ncard=%s' % (len(card), card)
        return SELOAD(lid_s0, seid, lid_se, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def uncross_reference(self) -> None:
        pass

    def safe_cross_reference(self, model, xref_errors):
        pass

    def raw_fields(self):
        list_fields = ['SELOAD', self.lid_s0, self.seid, self.lid_se]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SEEXCLD(BaseCard):
    """
    Partitioned Superelement Exclusion

    Defines grids that will be excluded during the attachment of a
    partitioned superelement.

    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    |    1    |   2   |   3   |   4   |   5   |   6   |   7   |   8   |   9   |
    +=========+=======+=======+=======+=======+=======+=======+=======+=======+
    | SEEXCLD | SEIDA | SEIDB | GIDA1 | GIDA2 | GIDA3 | GIDA4 | GIDA5 | GIDA6 |
    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    |         | GIDA7 | GIDA8 |  etc. |       |       |       |       |       |
    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    """
    type = 'SEEXCLD'
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        seid_a = 1
        seid_b = 2
        nodes = [10, 20, 30]
        return SEEXCLD(seid_a, seid_b, nodes, comment='')

    def __init__(self, seid_a, seid_b, nodes, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.seid_a = seid_a
        self.seid_b = seid_b

        #:  Identifiers of grids points. (Integer > 0)
        self.nodes = expand_thru(nodes)
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid_a = integer(card, 1, 'seid_a')
        seid_b = integer_string_or_blank(card, 2, 'seid_b')
        nodes = []
        i = 1
        nfields = len(card)
        for ifield in range(3, nfields):
            nid = integer_string_or_blank(card, ifield, 'nid_%i' % i)
            if nid:
                i += 1
                nodes.append(nid)
        assert len(card) >= 3, 'len(SEEXCLD card) = %i\ncard=%s' % (len(card), card)
        return SEEXCLD(seid_a, seid_b, nodes, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    @property
    def node_ids(self):
        return self.nodes

    def raw_fields(self):
        list_fields = ['SEEXCLD', self.seid_a, self.seid_b, ] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SEMPLN(BaseCard):
    """
    Superelement Mirror Plane
    Defines a mirror plane for mirroring a partitioned superelement.

    +--------+------+-------+----+----+------+
    |    1   |   2  |    3  | 4  | 5  |   6  |
    +========+======+=======+====+====+======+
    | SEMPLN | SEID | PLANE | P1 | P2 |  P3  |
    +--------+------+-------+----+----+------+
    | SEMPLN | 110  | PLANE | 12 | 45 | 1125 |
    +--------+------+-------+----+----+------+
    """
    type = 'SEMPLN'
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        seid = 1
        p1 = 2
        p2 = 3
        p3 = 4
        return SEMPLN(seid, p1, p2, p3, comment='')

    def __init__(self, seid, p1, p2, p3, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.seid = seid
        self.nodes = [p1, p2, p3]
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'lid_s0')
        plane = string(card, 2, 'seid')
        p1 = integer(card, 3, 'p1')
        p2 = integer(card, 4, 'p2')
        p3 = integer(card, 5, 'p3')
        assert plane == 'PLANE', plane
        assert len(card) <= 6, 'len(SEMPLN card) = %i\ncard=%s' % (len(card), card)
        return SEMPLN(seid, p1, p2, p3, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SEMPLN seid=%s' % self.seid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SEMPLN seid=%s' % self.seid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self):
        return _node_ids(self, self.nodes, self.nodes_ref, allow_empty_nodes=False, msg='')

    def raw_fields(self):
        list_fields = ['SEMPLN', self.seid, 'PLANE'] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SELABEL(BaseCard):
    """
    Superelement Output Label

    Defines a label or name to be printed in the superelement output headings.

    +---------+------+---------------------------------+
    |    1    |   2  |                3                |
    +=========+======+=================================+
    | SELABEL | SEID | LABEL                           |
    +---------+------+---------------------------------+
    | SELABEL |  10  | LEFT REAR FENDER, MODEL XYZ2000 |
    +---------+------+---------------------------------+
    """
    type = 'SELABEL'

    @classmethod
    def _init_from_empty(cls):
        seid = 1
        label = 'LEFT REAR FENDER'
        return SELABEL(seid, label, comment='')

    def __init__(self, seid, label, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.seid = seid
        self.label = label

    def validate(self):
        assert isinstance(self.label, str), self.label

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        label = ''.join([exact_string_or_blank(card, ifield, 'label', '        ')
                         for ifield in range(2, len(card))])
        return SELABEL(seid, label, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        pass

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        return [self.write_card()]

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = 'SELABEL %-8s%s\n' % (self.seid, self.label)
        return self.comment + card


class SELOC(BaseCard):
    """
    Partitioned Superelement Location

    Defines a partitioned superelement relocation by listing three non-colinear points in
    the superelement and three corresponding points not belonging to the superelement.

    +-------+------+-----+-----+-----+------+-----+-----+
    |   1   |   2  |  3  |  4  |  5  |   6  |  7  |  8  |
    +=======+======+=====+=====+=====+======+=====+=====+
    | SELOC | SEID | PA1 | PA2 | PA3 |  PB1 | PB2 | PB3 |
    +-------+------+-----+-----+-----+------+-----+-----+
    | SELOC | 110  | 10  | 100 | 111 | 1010 | 112 | 30  |
    +-------+------+-----+-----+-----+------+-----+-----+

    """
    type = 'SELOC'
    _properties = ['nodes_0_ids', 'nodes_seid_ids']

    @classmethod
    def _init_from_empty(cls):
        seid = 1
        nodes_seid = [1, 2, 3]
        nodes0 = 42
        return SELOC(seid, nodes_seid, nodes0, comment='')

    def __init__(self, seid, nodes_seid, nodes0, comment=''):
        """
        Creates an SELOC card, which transforms the superelement SEID
        from PA to PB.  Basically, define two CORD1Rs.

        Parameters
        ----------
        seid : int
            the superelement to transform
        nodes_seid : List[int, int, int]
            the nodes in the superelement than define the resulting coordinate system
        nodes0 : List[int, int, int]
            the nodes in the superelement than define the starting coordinate system
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.seid = seid
        #:  Identifiers of grids points. (Integer > 0)
        self.nodes_0 = expand_thru(nodes0, set_fields=False, sort_fields=False)
        self.nodes_seid = expand_thru(nodes_seid, set_fields=False, sort_fields=False)
        self.nodes_0_ref = None
        self.nodes_seid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        nodes0 = []
        nodes_seid = []
        i = 1
        fields = card[9:]
        nfields = len(fields)
        assert nfields % 2 == 0, fields
        for ifield in [2, 3, 4]:
            nid_a = integer(card, ifield, 'nid_%i' % i)
            nodes_seid.append(nid_a)
        for ifield in [5, 6, 7]:
            nid_b = integer(card, ifield, 'nid_%i' % i)
            nodes0.append(nid_b)

        assert len(card) <= 8, 'len(SELOC card) = %i\ncard=%s' % (len(card), card)
        return SELOC(seid, nodes_seid, nodes0, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SELOC seid=%s' % (self.seid)

        #PA1-PA3 Three GRID entries in the PART that are to be used to move the PART. After moving,
        #these points will be coincident with PB1-PB3.
        #
        # Three GRID entries
        self.nodes_seid_ref = model.superelement_nodes(self.seid, self.nodes_seid, msg=msg)

        #PB1-PB3 Three points (either GRID or POINT entries) defined in the Main Bulk Data Section
        #that define where the PART should be.
        #
        # either GRID or POINT entries
        self.nodes_0_ref = model.get_point_grids(self.nodes_0, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SELOC seid=%s' % (self.seid)
        self.nodes_seid_ref = model.superelement_nodes(self.seid, self.nodes_seid, msg=msg)
        self.nodes_0_ref = model.get_point_grids(self.nodes_0, msg=msg)

    @property
    def nodes_seid_ids(self):
        return _node_ids(self, self.nodes_seid, self.nodes_seid_ref,
                         allow_empty_nodes=False, msg='')

    @property
    def nodes_0_ids(self):
        return _node_ids(self, self.nodes_0, self.nodes_0_ref, allow_empty_nodes=False, msg='')

    def transform(self, model, xyz_cid0):
        #if self.nodes_0_ref is None:
            #self.cross_reference(model)

        global_coord_ref = self.nodes_0_ref
        seid_coord_ref = self.nodes_seid_ref
        p123_0 = np.array([node.get_position() for node in global_coord_ref])
        p123_seid = np.array([node.get_position() for node in seid_coord_ref])
        #print('global_coord_ref:\n%s' % global_coord_ref)
        #print('seid_coord_ref:\n%s' % seid_coord_ref)
        #print('p123_seid:\n%s' % p123_seid)
        #print('p123_0:\n%s' % p123_0)
        cid = max(model.coords)
        coord_seid = model.add_cord2r(cid+1, p123_seid[0, :], p123_seid[1, :], p123_seid[2, :])
        coord_0 = model.add_cord2r(cid+2, p123_0[0, :], p123_0[1, :], p123_0[2, :])
        coord_0.setup()
        coord_seid.setup()
        #print('beta_seid:\n%s' % coord_seid.beta())
        #print('beta0:\n%s' % coord_0.beta())

        #print(coord_seid.get_stats())
        # TODO: coord xform:
        #   xform = coord0.T * coord_seid
        #   xform = coord_seid.T * coord0
        xform = coord_0.beta().T @ coord_seid.beta()
        #print('xform%i:\n%s' % (self.seid, xform))
        dorigin = p123_0[0, :] - p123_seid[0, :] # at least, I'm sure on this...
        del model.coords[cid + 1]
        del model.coords[cid + 2]

        # TODO: not 100% on this xform
        xyz_cid0 = xyz_cid0.dot(xform.T) + dorigin
        return xyz_cid0

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes_seid = self.nodes_seid_ids
        self.nodes_0 = self.nodes_0_ids
        self.nodes_0_ref = None
        self.nodes_seid_ref = None

    def raw_fields(self):
        list_fields = ['SELOC', self.seid] + list(self.nodes_seid_ids) + list(self.nodes_0_ids)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SETREE(BaseCard):
    """
    Superelement Tree Definition (Alternate Form of DTI,SETREE)

    Specifies superelement reduction order.

    +--------+-------+-------+-------+-------+-------+-------+-------+-------+
    |    1   |    2  |   3   |   4   |   5   |   6   |   7   |   8   |   9   |
    +========+=======+=======+=======+=======+=======+=======+=======+=======+
    | SETREE |  SEID | SEUP1 | SEUP2 | SEUP3 | SEUP4 | SEUP5 | SEUP6 | SEUP7 |
    +--------+-------+-------+-------+-------+-------+-------+-------+-------+
    |        | SEUP8 | SEUP9 |  etc. |       |       |       |       |       |
    +--------+-------+-------+-------+-------+-------+-------+-------+-------+
    | SETREE |  400  |   10  |   20  |   30  |   40  |       |       |       |
    +--------+-------+-------+-------+-------+-------+-------+-------+-------+
    """
    type = 'SETREE'

    @classmethod
    def _init_from_empty(cls):
        seid = 10
        superelements = [1, 2, 3]
        return SETREE(seid, superelements, comment='')

    def __init__(self, seid, superelements, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.seid = seid
        #:  Identifiers of grids points. (Integer > 0)
        self.superelements = expand_thru(superelements)
        self.superelements_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        superelements = []
        i = 1
        nfields = len(card)
        for ifield in range(2, nfields):
            superelement = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if superelement:
                i += 1
                superelements.append(superelement)
        assert len(card) >= 3, 'len(SETREE card) = %i\ncard=%s' % (len(card), card)
        return SETREE(seid, superelements, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SETREE seid=%s' % self.seid
        missing_superelements = []
        superelements_ref = []
        for super_id in self.superelements:
            if super_id in model.superelement_models:
                superelement = model.superelement_models[super_id]
            else:
                missing_superelements.append(super_id)
                continue
            superelements_ref.append(superelement)
        if missing_superelements:
            raise KeyError('cannot find superelements=%s%s' % (missing_superelements, msg))
        self.superelements_ref = superelements_ref

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
        pass

    def raw_fields(self):
        list_fields = ['SETREE', self.seid] + list(self.superelements)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class CSUPER(BaseCard):
    """
    Secondary Superelement Connection

    Defines the grid or scalar point connections for identical or mirror image
    superelements or superelements from an external source. These are all known as
    secondary superelements.

    +--------+------+------+------+-----+-----+-----=-----+-----+
    |    1   |   2  |   3  |   4  |  5  |  6  |  7  |  8  |  9  |
    +========+======+======+======+=====+=====+=====+=====+=====+
    | CSUPER | SSlD | PSID |  GP1 | GP2 | GP3 | GP4 | GP5 | GP6 |
    +--------+------+------+------+-----+-----+-----=-----+-----+
    |        |  GP7 |  GP8 | etc. |     |     |     |     |     |
    +--------+------+------+------+-----+-----+-----=-----+-----+
    """
    type = 'CSUPER'
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        seid = 1
        psid = 1
        nodes = [1, 2]
        return CSUPER(seid, psid, nodes, comment='')

    def __init__(self, seid, psid, nodes, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.seid = seid
        self.psid = psid
        #:  Identifiers of grids points. (Integer > 0)
        self.nodes = expand_thru(nodes)
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        psid = integer_or_blank(card, 2, 'psid', 0)
        nodes = []
        i = 1
        nfields = len(card)
        for ifield in range(3, nfields):
            nid = integer_string_or_blank(card, ifield, 'nid_%i' % i)
            if nid:
                i += 1
                nodes.append(nid)
        assert len(card) >= 3, 'len(CSUPER card) = %i\ncard=%s' % (len(card), card)
        return CSUPER(seid, psid, nodes, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CSUPER seid=%s' % self.seid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CSUPER seid=%s' % self.seid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self):
        return _node_ids(self, self.nodes, self.nodes_ref, allow_empty_nodes=False, msg='')

    def raw_fields(self):
        list_fields = ['CSUPER', self.seid, self.psid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CSUPEXT(BaseCard):
    """
    Superelement Exterior Point Definition

    Assigns exterior points to a superelement.

    +---------+------+-----+-----+-----+-----+-----+-----+-----+
    |    1    |   2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +=========+======+=====+=====+=====+=====+=====+=====+=====+
    | CSUPEXT | SEID | GP1 | GP2 | GP3 | GP4 | GP5 | GP6 | GP7 |
    +---------+------+-----+-----+-----+-----+-----+-----+-----+
    """
    type = 'CSUPEXT'
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        seid = 1
        nodes = [1]
        return CSUPEXT(seid, nodes, comment='')

    def __init__(self, seid, nodes, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.seid = seid
        #:  Identifiers of grids points. (Integer > 0)
        self.nodes = expand_thru(nodes)
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        nodes = []
        i = 1
        nfields = len(card)
        for ifield in range(2, nfields):
            nid = integer_string_or_blank(card, ifield, 'node_%i' % i)
            if nid:
                i += 1
                nodes.append(nid)
        assert len(card) <= 9, 'len(CSUPEXT card) = %i\ncard=%s' % (len(card), card)
        return CSUPEXT(seid, nodes, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CSUPEXT eid=%s' % self.seid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CSUPEXT eid=%s' % self.seid
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self):
        return _node_ids(self, self.nodes, self.nodes_ref, allow_empty_nodes=False, msg='')

    def raw_fields(self):
        list_fields = ['CSUPEXT', self.seid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SEBULK(BaseCard):
    """
    Partitional Superelement Connection

    Defines superelement boundary search options and a repeated,
    mirrored, or collector superelement.

    +--------+------+--------+-------+--------+--------+-----+--------+
    |    1   |   2  |    3   |    4  |    5   |   6    |  7  |    8   |
    +========+======+========+=======+========+========+=====+========+
    | SEBULK | SEID |  TYPE  | RSEID | METHOD |  TOL   | LOC | UNITNO |
    +--------+------+--------+-------+--------+--------+-----+--------+
    | SEBULK | 14   | REPEAT |   4   |  AUTO  | 1.0E-3 |     |        |
    +--------+------+--------+-------+--------+--------+-----+--------+

    """
    type = 'SEBULK'

    @classmethod
    def _init_from_empty(cls):
        seid = 1
        superelement_type = 'MIRROR'
        rseid = 42
        return SEBULK(seid, superelement_type, rseid,
                      method='AUTO', tol=1e-5, loc='YES', unitno=None, comment='')

    def __init__(self, seid, superelement_type, rseid,
                 method='AUTO', tol=1e-5, loc='YES', unitno=None,
                 comment=''):
        """
        Parameters
        ----------
        seid : int
            Partitioned superelement identification number.
        Type : str
            Superelement type.
            {PRIMARY, REPEAT, MIRROR, COLLCTR, EXTERNAL, EXTOP2}
        rseid : int; default=0
            Identification number of the reference superelement,
            used if TYPE = 'REPEAT' and 'MIRROR'.
        method : str; default='AUTO'
            Method to be used when searching for boundary grid points.
            {AUTO, MANUAL}
        tol : float; default=1e-5
            Location tolerance to be used when searching for boundary grid points.
        loc : str; default='YES'
            Coincident location check option for manual conection option.
            {YES, NO}
        unitno : int / None
            FORTRAN unit number for the OUTPUT2 file (applicable and
            meaningful only when TYPE='EXTOP2').

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.seid = seid
        self.superelement_type = superelement_type
        self.rseid = rseid
        self.method = method
        self.tol = tol
        self.loc = loc
        self.unitno = unitno

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        superelement_type = string(card, 2, 'superelement_type')
        rseid = integer_or_blank(card, 3, 'rseid', 0)
        method = string_or_blank(card, 4, 'method', 'AUTO')
        tol = double_or_blank(card, 5, 'tol', 1e-5)
        loc = string_or_blank(card, 6, 'loc', 'YES')
        unitno = integer_or_blank(card, 7, 'seid')
        assert len(card) <= 8, 'len(SEBULK card) = %i\ncard=%s' % (len(card), card)
        return SEBULK(seid, superelement_type, rseid, method=method, tol=tol,
                      loc=loc, unitno=unitno, comment=comment)

    def validate(self):
        assert self.superelement_type in ['PRIMARY', 'REPEAT', 'MIRROR', 'COLLCTR', 'EXTERNAL', 'EXTOP2', 'FRFOP2', 'MANUAL'], f'superelement_type={self.superelement_type}\n{self}'
        assert self.loc in ['YES', 'NO'], self.loc
        assert self.method in ['AUTO', 'MANUAL'], self.method

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        pass

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        list_fields = [
            'SEBULK', self.seid, self.superelement_type, self.rseid, self.method, self.tol,
            self.loc, self.unitno]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SECONCT(BaseCard):
    """
    Partitioned Superelement Boundary-Point Connection

    Explicitly defines grid and scalar point connection procedures for a
    partitioned superelement.

    +---------+-------+-------+--------+-------+-------+-------+------+------+
    |    1    |   2   |   3   |    4   |   5   |   6   |   7   |   8  |   9  |
    +=========+=======+=======+========+=======+=======+=======+======+======+
    | SECONCT | SEIDA | SEIDB |   TOL  |  LOC  |       |       |      |      |
    +---------+-------+-------+--------+-------+-------+-------+------+------+
    |         | GIDA1 | GIDB1 |  GIDA2 | GIDB2 | GIDA3 | GIDB3 | etc. | etc. |
    +---------+-------+-------+--------+-------+-------+-------+------+------+
    | SECONCT |   10  |   20  | 1.0E-4 |  YES  |       |       |      |      |
    +---------+-------+-------+--------+-------+-------+-------+------+------+
    |         |  1001 |  4001 |        |       |  2222 |  4444 |      |      |
    +---------+-------+-------+--------+-------+-------+-------+------+------+
    | SECONCT | SEIDA | SEIDB |   TOL  |  LOC  |       |       |      |      |
    +---------+-------+-------+--------+-------+-------+-------+------+------+
    |         | GIDA1 | THRU  |  GIDA2 | GIDB1 |  THRU | GIDB2 |      |      |
    +---------+-------+-------+--------+-------+-------+-------+------+------+
    | SECONCT |  10   |   20  |        |       |       |       |      |      |
    +---------+-------+-------+--------+-------+-------+-------+------+------+
    |         |  101  |  THRU |   110  |  201  |  THRU |  210  |      |      |
    +---------+-------+-------+--------+-------+-------+-------+------+------+
    """
    type = 'SECONCT'
    _properties = ['node_ids_a', 'node_ids_b']

    @classmethod
    def _init_from_empty(cls):
        seid_a = 1
        seid_b = 2
        tol = 0.1
        loc = 'YES'
        nodes_a = [10, 20, 30]
        nodes_b = [11, 21, 31]
        return SECONCT(seid_a, seid_b, tol, loc, nodes_a, nodes_b, comment='')

    def __init__(self, seid_a, seid_b, tol, loc, nodes_a, nodes_b, comment=''):
        """
        Parameters
        ----------
        SEIDA : int
            Partitioned superelement identification number.
        SEIDB : int
            Identification number of superelement for connection to SEIDA.
        TOL : float; default=1e-5
            Location tolerance to be used when searching for or checking boundary
            grid points.
        LOC : str; default='YES'
            Coincident location check option for manual connection.
            {YES, NO}
        GIDAi : int
            Identification number of a grid or scalar point in superelement SEIDA,
            which will be connected to GIDBi.
        GIDBi : int
            Identification number of a grid or scalar point in superelement SEIDB,
            which will be connected to GIDAi.

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.seid_a = seid_a
        self.seid_b = seid_b
        self.tol = tol
        self.loc = loc
        self.nodes_a = nodes_a
        self.nodes_b = nodes_b
        self.nodes_a_ref = None
        self.nodes_b_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid_a = integer(card, 1, 'seid_a')
        seid_b = integer(card, 2, 'seid_b')
        tol = double_or_blank(card, 3, 'tol', 1e-5)
        loc = string_or_blank(card, 4, 'loc', 'YES')
        fields = card[9:]
        if len(fields) < 2:
            assert len(card) >= 9, 'len(SECONCT card) = %i\ncard=%s' % (len(card), card)

        assert len(fields) % 2 == 0, 'card=%s\nfields=%s' % (card, fields)
        if 'THRU' in fields:
            raise NotImplementedError(fields)
            #start_a = integer(card, 9, 'start_a')
            #thru_a = string(card, 10, 'thru_a')
            #end_a = integer(card, 11, 'end_a')

            #start_b = integer(card, 12, 'start_b')
            #thru_b = string(card, 13, 'thru_b')
            #end_b = integer(card, 14, 'end_b')
            #assert thru_a == 'THRU', thru_a
            #assert thru_b == 'THRU', thru_b
            #nodes_a = list(range(start_a+1, end_a+1))
            #nodes_b = list(range(start_b+1, end_b+1))
            #print(nodes_a)
        else:
            nodes_a = []
            nodes_b = []
            inode = 1
            for ifield in range(0, len(fields), 2):
                node_a = integer_or_blank(card, 9+ifield, 'node_a%i' % inode)
                node_b = integer_or_blank(card, 9+ifield+1, 'node_b%i' % inode)
                if node_a is None and node_b is None:
                    continue
                assert node_a is not None, fields
                assert node_b is not None, fields
                nodes_a.append(node_a)
                nodes_b.append(node_b)
                inode += 1
        return SECONCT(seid_a, seid_b, tol, loc, nodes_a, nodes_b, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SECONCT seid_a=%s seid_b=%s' % (self.seid_a, self.seid_b)
        self.nodes_a_ref = model.superelement_nodes(self.seid_a, self.nodes_a, msg=msg)
        self.nodes_b_ref = model.superelement_nodes(self.seid_b, self.nodes_b, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SECONCT seid_a=%s seid_b=%s' % (self.seid_a, self.seid_b)
        self.nodes_a_ref = model.superelement_nodes(self.seid_a, self.nodes_a, msg=msg)
        self.nodes_b_ref = model.superelement_nodes(self.seid_b, self.nodes_b, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes_a = self.node_ids_a
        self.nodes_b = self.node_ids_b
        self.nodes_a_ref = None
        self.nodes_b_ref = None

    @property
    def node_ids_a(self):
        return _node_ids(self, self.nodes_a, self.nodes_a_ref, allow_empty_nodes=False, msg='')

    @property
    def node_ids_b(self):
        return _node_ids(self, self.nodes_b, self.nodes_b_ref, allow_empty_nodes=False, msg='')

    def raw_fields(self):
        list_fields = ['SECONCT', self.seid_a, self.seid_b, self.tol, self.loc,
                       None, None, None, None,]
        for (nid_a, nid_b) in zip(self.node_ids_a, self.node_ids_b):
            list_fields += [nid_a, nid_b]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SENQSET(BaseCard):
    """
    Superelement Internal Generalized Degree-of-Freedom

    Defines number of internally generated scalar points for superelement dynamic
    reduction.

    +---------+------+----+
    |    1    |   2  | 3  |
    +---------+------+----+
    | SENQSET | SEID | N  |
    +---------+------+----+
    | SENQSET | 110  | 45 |
    +---------+------+----+
    """
    type = 'SENQSET'

    @classmethod
    def _init_from_empty(cls):
        set_id = 1
        n = 45
        return SENQSET(set_id, n, comment='')

    def __init__(self, set_id, n, comment=''):
        """
        Parameters
        ----------
        set_id : int / str
            Partitioned superelement identification number.
            (Integer > 0 or Character='ALL')
        n : int; default=0
            Number of internally generated scalar points for dynamic
            reduction generalized coordinates (Integer > 0).
        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.set_id = set_id
        self.n = n

    @classmethod
    def add_card(cls, card, comment=''):
        set_id = integer_or_string(card, 1, 'set_id')
        n = integer_or_blank(card, 2, 'n', 0)
        assert len(card) <= 3, 'len(SENQSET card) = %i\ncard=%s' % (len(card), card)
        return SENQSET(set_id, n, comment=comment)

    def raw_fields(self):
        list_fields = ['SENQSET', self.set_id, self.n]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

def _node_ids(card, nodes, nodes_ref, allow_empty_nodes=False, msg=''):
    if nodes_ref is None:
        #nodes = card.nodes
        assert nodes is not None, card.__dict__
        return nodes

    try:
        if allow_empty_nodes:
            nodes2 = []
            for node in nodes_ref:
                if node == 0 or node is None:
                    nodes2.append(None)
                elif isinstance(node, integer_types):
                    nodes2.append(node)
                else:
                    nodes2.append(node.nid)
            assert nodes2 is not None, str(card)
            return nodes2
        else:
            try:
                node_ids = []
                for node in nodes_ref:
                    if isinstance(node, integer_types):
                        node_ids.append(node)
                    else:
                        node_ids.append(node.nid)

                #if isinstance(nodes[0], integer_types):
                    #node_ids = [node for node in nodes]
                #else:
                    #node_ids = [node.nid for node in nodes]
            except:
                print('type=%s nodes=%s allow_empty_nodes=%s\nmsg=%s' % (
                    card.type, nodes, allow_empty_nodes, msg))
                raise
            assert 0 not in node_ids, 'node_ids = %s' % node_ids
            assert node_ids is not None, str(card)
            return node_ids
    except:
        print('type=%s nodes=%s allow_empty_nodes=%s\nmsg=%s' % (
            card.type, nodes, allow_empty_nodes, msg))
        raise
    raise RuntimeError('huh...')
