from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double

from pyNastran.bdf.bdf_interface.assign_type import (
    string, string_or_blank, integer, integer_or_blank, double_or_blank,
    integer_or_double)
from pyNastran.bdf.cards.loads.static_loads import Load
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

# EQUIV
# DELETE
# ACTIVE
# RSPLNG

# CGEN
# SPCG
# SETG
# TEMPG
# ETEMP
# LIST
# EDGER
# INSECT
# GRIDMOD
# DISTORT

# TRICON
# CBARG
# PGEN

# EIDL
# EIDH

# GRIDMOD
# GRIDD
# GRIDG
# GRIDU
# EGRID

# CURVE
# EDGER

# MPOINT
# MDIM
# MLINE
# MAREA
# MTAB
# MFUN

# MESHOFF
# MESHON

# PLOTG
# PLOTG
# PLOTOPT

# PLOADG
class CGEN(BaseCard):
    """
    +------+--------+-----------+-----+----------+-----+-------------+------+------+
    |   1  |    2   |     3     | 4   |     5    |  6  |      7      |  8   |  9   |
    +======+========+===========+=====+==========+=====+=============+======+======+
    | CGEN |  TYPE  | FIELD_EID | PID | FIELD_ID | DIR | TH_GEOM_OPT | EIDL | EIDH |
    +------+--------+-----------+-----+----------+-----+-------------+------+------+
    | CGEN | TRIA3  |    550    | 78  |    25    |     |             |      |      |
    +------+--------+-----------+-----+----------+-----+-------------+------+------+
    | CGEN | TRIA6  |    520    | 78  |    26    |     |             |      |      |
    +------+--------+-----------+-----+----------+-----+-------------+------+------+
    | CGEN | QUAD4  |    450    | 145 |    33    |     |             |      |      |
    +------+--------+-----------+-----+----------+-----+-------------+------+------+
    | CGEN | HEXA8  |    610    | 57  |    41    |     |             |      |      |
    +------+--------+-----------+-----+----------+-----+-------------+------+------+
    | CGEN | HEXA20 |    620    | 57  |    42    |     |             |      |      |
    +------+--------+-----------+-----+----------+-----+-------------+------+------+
    """
    type = 'CGEN'

    def __init__(self, Type, field_eid, pid, field_id, th_geom_opt,
                 eidl, eidh, t_abcd=None, direction='L', comment=''):
        """
        Creates the CGEN card

        Parameters
        ----------
        Type : str
            LINE, TRIA, QUAD, HEX
        field_eid : int
           starting element id
        pid : int
           property id
        field_id : int
            GRIDG id
        cdir : str
            L, M, N
        th_geom_opt :
            TH : ???
                ???
            GEOM : ???
                only when Type = BEND
            OPT : ???
                ???
        eidl : int
            ???
        eidh : int
            ???
        """
        if comment:
            self.comment = comment
        self.Type = Type
        self.field_eid = field_eid
        self.pid = pid
        self.field_id = field_id
        self.direction = direction
        self.th_geom_opt = th_geom_opt
        self.eid = eidl
        self.eidh = eidh
        if t_abcd is None:
            t_abcd = []
        self.t_abcd = t_abcd
        self.nodes = []

    @property
    def node_ids(self):
        return []

    def validate(self):
        pass

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CGEN card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """

        Type = string(card, 1, 'Type')
        feid = integer(card, 2, 'feid')
        pid = integer(card, 3, 'pid')
        fid = integer(card, 4, 'fid')
        direction = string_or_blank(card, 5, 'direction', 'L')
        th = integer_or_blank(card, 6, 'th')
        eidl = integer(card, 7, 'eidl')
        eidh = integer(card, 8, 'eidh')
        t_abcd = [
            double_or_blank(card, 9, 'ta', 0.),
            double_or_blank(card, 10, 'tb', 0.),
            double_or_blank(card, 11, 'tc', 0.),
            double_or_blank(card, 12, 'td', 0.),
        ]
        return CGEN(Type, feid, pid, fid, th, eidl, eidh,
                    t_abcd=t_abcd, direction=direction,
                    comment=comment)

    def _verify(self, xref):
        pass

    def cross_reference(self, model: BDF) -> None:
        pass
    def safe_cross_reference(self, model):
        pass
    def uncross_reference(self, model):
        pass



    #def _validate_input(self):
        #assert self.nid > 0, 'nid=%s' % (self.nid)
        #assert self.cp >= 0, 'cp=%s' % (self.cp)
        #assert self.cd >= -1, 'cd=%s' % (self.cd)
        #assert self.seid >= 0, 'seid=%s' % (self.seid)
        #assert len(self.xyz) == 3

    #def _verify(self, xref):
        #"""
        #Verifies all methods for this object work

        #Parameters
        #----------
        #xref : bool
            #has this model been cross referenced
        #"""
        #pass

    #def cross_reference(self, model: BDF, grdset=None):
        #"""
        #Cross links the card so referenced cards can be extracted directly

        #Parameters
        #----------
        #model : BDF()
            #the BDF object
        #grdset : GRDSET / None; default=None
            #a GRDSET if available (default=None)

        #.. note::  The gridset object will only update the fields that
                   #have not been set
        #"""
        #if grdset:
            ## update using a gridset object
            #if not self.cp:
                #self.cp = grdset.cp
                #self.cp_ref = self.cp
            #if not self.cd:
                #self.cd = grdset.cd
                #self.cd_ref = self.cd
            #if not self.ps:
                #self.ps = grdset.ps
                #self.ps_ref = self.ps
            #if not self.seid:
                #self.seid = grdset.seid
                #self.seid_ref = self.seid
        #msg = ', which is required by CGEN nid=%s' % (self.nid)
        #self.cp = model.Coord(self.cp, msg=msg)
        #self.cp_ref = self.cp
        #if self.cd != -1:
            #self.cd = model.Coord(self.cd, msg=msg)
            #self.cd_ref = self.cd

    #def uncross_reference(self) -> None:
        #self.cp = self.Cp()
        #self.cd = self.Cd()
        #del self.cp_ref, self.cd_ref
        #if hasattr(self, 'elements'):
            #del self.elements

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[int/float/str]
            the fields that define the card
        """
        list_fields = [
            'CGEN', self.Type, self.field_eid, self.pid, self.field_id,
            self.direction, self.th_geom_opt, self.eid, self.eidh] + self.t_abcd
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[int/float/str]
            the fields that define the card
        """
        list_fields = [
            'CGEN', self.Type, self.field_eid, self.pid, self.field_id,
            self.direction, self.th_geom_opt, self.eid, self.eidh] + self.t_abcd
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class GMCORD(BaseCard):
    """defines the GMCOORD class

    GMCORD | CID | ENTITY | ID1 | ID2 |
    GMCORD | 101 | GMCURV |  26 |  44 |
    """
    type = 'GMCORD'

    @classmethod
    def _init_from_empty(cls):
        cid = 1
        entity = 'GMCURV'
        gm_ids = [34]
        return GMCORD(cid, entity, gm_ids, comment='')

    def __init__(self, cid, entity, gm_ids, comment=''):
        """Creates a GMCOORD"""
        if comment:
            self.comment = comment
        self.cid = cid
        self.entity = entity
        self.gm_ids = gm_ids

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a GMCORD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        cid = integer(card, 1, 'cid')
        entity = string(card, 2, 'entity')
        gm_ids = [
            integer(card, 3, 'GM_ID1'),
            integer_or_blank(card, 4, 'GM_ID2'),
        ]
        return GMCORD(cid, entity, gm_ids, comment=comment)

    @staticmethod
    def cross_reference(model: BDF) -> None:
        pass

    @staticmethod
    def uncross_reference() -> None:
        """Removes cross-reference links"""
        pass

    def setup(self):
        pass

    def raw_fields(self):
        list_fields = ['GMCORD', self.cid, self.entity] + self.gm_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class GMLOAD(Load):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and two grid points that determine the direction.

    """
    type = 'GMLOAD'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        normal = [0., 0., 1.]
        entity = 1
        entity_id = 1
        method = 1
        load_magnitudes = [1., 2.]
        return GMLOAD(sid, normal, entity, entity_id, method, load_magnitudes, cid=0, comment='')

    def __init__(self, sid, normal, entity, entity_id, method,
                 load_magnitudes, cid=0, comment=''):
        """Creates a GMLOAD object"""
        Load.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.cid = cid
        self.normal = normal
        self.entity = entity
        self.entity_id = entity_id
        self.method = method
        self.load_magnitudes = load_magnitudes
        self.cid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a GMLOAD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        cid = integer_or_blank(card, 2, 'cid', default=0)
        normal = np.array([
            double_or_blank(card, 3, 'N1', default=0.),
            double_or_blank(card, 4, 'N2', default=0.),
            double_or_blank(card, 5, 'N3', default=1.),
        ])
        entity = string(card, 6, 'entity')
        entity_id = integer(card, 7, 'entity_id')
        method = string(card, 8, 'method')

        load_magnitudes = []
        for i in range(9, len(card)):
            ifield = i - 8
            load_mag = integer_or_double(card, i, 'load_magnitude_%s' % ifield)
            load_magnitudes.append(load_mag)
        return GMLOAD(sid, normal, entity, entity_id, method,
                      load_magnitudes, cid=cid, comment=comment)

    #def DEquation(self):
        #if isinstance(self.dequation, int):
            #return self.dequation
        #return self.dequation.equation_id

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by GMLOAD sid=%s' % self.sid
        self.cid_ref = model.Coord(self.Cid(), msg=msg)
        #self.node = model.Node(self.node, msg=msg)
        #self.g1 = model.Node(self.g1, msg=msg)
        #self.g2 = model.Node(self.g2, msg=msg)
        #self.xyz = self.g2.get_position() - self.g1.get_position()
        #normalize(self, msg)

    def safe_cross_reference(self, model: BDF, xref_errors):
        msg = ', which is required by GMLOAD sid=%s' % self.sid
        self.cid_ref = model.safe_coord(self.Cid(), self.sid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cid = self.Cid()
        self.cid_ref = None

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    #def G1(self):
        #if isinstance(self.g1, (integer_types, float)):
            #return self.g1
        #return self.g1_ref.nid

    #def G2(self):
        #if isinstance(self.g2, (integer_types, float)):
            #return self.g2
        #return self.g2_ref.nid

    #def NodeID(self):
        #if isinstance(self.node, integer_types):
            #return self.node
        #return self.node_ref.nid

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = ['GMLOAD', self.sid, self.Cid()] + list(self.normal) + [
            self.entity, self.entity_id, self.method] + self.load_magnitudes
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


