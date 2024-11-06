from __future__ import annotations
from io import StringIO
from itertools import zip_longest
from collections import Counter, defaultdict
from typing import Callable, Optional, Any, TYPE_CHECKING
import numpy as np
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_scientific_16
from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.nodes import compress_xpoints
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru, _format_comment
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, components_or_blank,
    integer_or_string, blank)
from pyNastran.bdf.bdf_interface.assign_type_force import (
    lax_double_or_blank) # force_double_or_blank,

#from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_float, array_default_int)
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, parse_check, sort_duplicates)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import get_print_card_size, update_field_size
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check

from pyNastran.femutils.coord_transforms import (
    xyz_to_rtz_array,
    xyz_to_rtp_array, # xyz to xxx transforms
    rtz_to_xyz_array as cylindrical_to_rectangular,
    rtp_to_xyz_array as spherical_to_rectangular, # xxx to xyz transforms
    #rtz_to_rtp_array, rtp_to_rtz_array, # rtp/rtz and rtz/rtp transforms
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class XPOINT(VectorizedBaseCard):
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
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.ids = np.array([], dtype='int32')

    def add(self, ids: np.ndarray, comment: str='') -> int:
        if isinstance(ids, integer_types):
            ids = [ids]
        self.cards.append((ids, comment))
        self.n += len(ids)
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        ids = []
        for i in range(1, len(card)):
            field = integer_or_string(card, i, 'ID%d' % i)
            ids.append(field)
        self.cards.append((ids, comment))
        self.n += len(ids)
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        idtype = self.model.idtype
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        ids = []
        for icard, (idsi, comment) in enumerate(self.cards):
            idsi2 = expand_thru(idsi)
            assert min(idsi2) > 0, idsi
            ids.extend(idsi2)
            self.comment[icard] = comment

        self.ids = np.array(ids, dtype=idtype)
        self._sort()
        self.cards = []

    def _sort(self) -> None:
        iarg = np.argsort(self.ids)
        uarg = np.unique(iarg)
        #nvalues = len(self.node_id)
        if not np.array_equal(uarg, iarg):
            self.ids = self.ids[iarg]

    def __apply_slice__(self, xpoint: XPOINT, i: np.ndarray) -> None:
        xpoint.ids = self.ids[i]
        xpoint.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        bad_ids = self.ids[self.ids <= 0]
        assert self.ids.min() > 0, bad_ids

    @property
    def max_id(self) -> int:
        return self.ids.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        assert self.ids.min() > 0, self.ids[self.ids <= 0]
        print_card, size = get_print_card_size(size, self.max_id)

        #node_id = array_str(self.node_id, size=8)
        lists_fields = compress_xpoints(self.type, self.ids)
        for fields in lists_fields:
            msg = print_card(fields)
            bdf_file.write(msg)
        return

class SPOINT(XPOINT):
    _id_name = 'spoint_id'
    @property
    def spoint_id(self):
        return self.ids
    @spoint_id.setter
    def spoint_id(self, spoint_id: np.ndarray):
        self.ids = spoint_id

    #def slice_card(self, i: np.ndarray) -> SPOINT:
        #"""uses a node_index to extract SPOINTs"""
        #assert len(self.ids) > 0, self.ids
        #i = np.atleast_1d(np.asarray(i, dtype=self.ids.dtype))
        #i.sort()
        #spoint = SPOINT(self.model)
        #spoint.n = len(i)
        #spoint.ids = self.ids[i]
        #return spoint

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        pass
    def remove_unused(self, used_dict: [dict[str, np.ndarray]]) -> int:
        spoint_id = used_dict['spoint_id']
        return 0


class EPOINT(XPOINT):
    _id_name = 'epoint_id'
    @property
    def epoint_id(self) -> np.ndarray:
        return self.ids
    @epoint_id.setter
    def epoint_id(self, epoint_id: np.ndarray):
        self.ids = epoint_id

    #def slice_card(self, i: np.ndarray) -> EPOINT:
        #"""uses a node_index to extract SPOINTs"""
        #assert len(self.ids) > 0, self.ids
        #i = np.atleast_1d(np.asarray(i, dtype=self.ids.dtype))
        #i.sort()
        #epoint = EPOINT(self.model)
        #epoint.n = len(i)
        #epoint.ids = self.ids[i]
        #return epoint


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
    _field_map = {2:'cp', 6:'cd', 7:'ps', 8:'seid'}

    @classmethod
    def _init_from_empty(cls):
        cp = 0
        cd = 1
        ps = 0
        seid = 0
        return GRDSET(cp, cd, ps, seid, comment='')

    def __init__(self, cp: int=0, cd: int=0, ps: int=0, seid: int=0, comment: str=''):
        """
        Creates the GRDSET card

        Parameters
        ----------
        cp : int; default=0
            the xyz coordinate frame
        cd : int; default=0
            the analysis coordinate frame
        ps : int; default=0
            Additional SPCs in the analysis coordinate frame (e.g. 123).
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

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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

        cp = integer_or_blank(card, 2, 'cp', default=0)
        blank(card, 3, 'blank')
        blank(card, 4, 'blank')
        blank(card, 5, 'blank')
        cd = integer_or_blank(card, 6, 'cd', default=0)

        ps = integer_or_blank(card, 7, 'ps', default=0)
        seid = integer_or_blank(card, 8, 'seid', default=0)
        assert len(card) <= 9, f'len(GRDSET card) = {len(card):d}\ncard={card}'
        return GRDSET(cp, cd, ps, seid, comment=comment)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref: bool
            has this model been cross referenced

        """
        assert isinstance(self.cp, integer_types), 'cp=%r' % self.cp
        assert isinstance(self.cd, integer_types), 'cd=%r' % self.cd
        assert isinstance(self.ps, integer_types), 'ps=%r' % self.ps
        assert isinstance(self.seid, integer_types), 'seid=%r' % self.seid

    def raw_fields(self) -> list[Any]:
        list_fields = ['GRDSET', None, self.cp, None, None, None,
                       self.cd, self.ps, self.seid]
        return list_fields

    def write(self, size: int=8):
        cp = set_blank_if_default(self.cp, 0)
        cd = set_blank_if_default(self.cd, 0)
        ps = set_blank_if_default(self.ps, 0)
        seid = set_blank_if_default(self.SEid(), 0)
        list_fields = ['GRDSET', None, cp, None, None, None, cd, ps, seid]
        return print_card_8(list_fields)


class GRID(VectorizedBaseCard):
    _id_name = 'node_id'
    _show_attributes = ['node_id', 'xyz', 'cp', 'cd', 'ps', 'seid']
    #_remove_attributes = ['_xyz_cid0']
    def __init__(self, model: BDF):
        super().__init__(model)
        self._is_sorted = False
        self.clear()

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.node_id = np.array([], dtype='int32')
        self.cp = np.array([], dtype='int32')
        self.cd = np.array([], dtype='int32')
        self.ps = np.array([], dtype='int32')
        self.seid = np.array([], dtype='int32')
        self.xyz = np.zeros((0, 3), dtype='float64')
        self._xyz_cid0 = np.zeros((0, 3), dtype='float64')

    def add(self, nid: int, xyz: np.ndarray,
            cp: int=0, cd: int=0,
            ps: int=0, seid: int=0, comment: str='') -> int:
        ps = 0 if ps == '' else ps
        assert isinstance(ps, int), ps
        self.cards.append((nid, xyz, cp, cd, ps, seid, comment))
        if comment:
            self.comment[nid] = _format_comment(comment)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')
        nfields = len(card)
        #: Node ID
        nid = integer(card, 1, 'nid')

        GRID_CP_DEFAULT = GRID_CD_DEFAULT = GRID_PS_DEFAULT = GRID_SEID_DEFAULT = -100
        GRID_PS_DEFAULT_STR = str(GRID_PS_DEFAULT)

        #: Grid point coordinate system
        cp = integer_or_blank(card, 2, 'cp', default=GRID_CP_DEFAULT)

        fdouble_or_blank = lax_double_or_blank if self.model.is_lax_parser else double_or_blank

        #: node location in local frame
        xyz = [
            fdouble_or_blank(card, 3, 'x1', default=0.),
            fdouble_or_blank(card, 4, 'x2', default=0.),
            fdouble_or_blank(card, 5, 'x3', default=0.),
        ]

        if nfields > 6:
            #: Analysis coordinate system
            cd = integer_or_blank(card, 6, 'cd', default=GRID_CD_DEFAULT)

            #: SPC constraint
            ps_str = components_or_blank(card, 7, 'ps', default=GRID_PS_DEFAULT_STR)
            ps = int(ps_str)
            #u(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            seid = integer_or_blank(card, 8, 'seid', default=GRID_SEID_DEFAULT)
            assert len(card) <= 9, f'len(GRID card) = {len(card):d}\ncard={card}'
        else:
            cd = GRID_CD_DEFAULT
            ps = GRID_PS_DEFAULT
            seid = GRID_SEID_DEFAULT
            assert len(card) >= 2, f'len(GRID card) = {len(card):d}\ncard={card}'
        self.cards.append((nid, xyz, cp, cd, ps, seid, comment))
        if comment:
            self.comment[nid] = comment
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug('parse GRID')

        try:
            node_id, cp, cd, xyz, ps, seid, comment = self._setup(
                ncards, self.cards, 'int32', 'float64')
        except OverflowError:
            node_id, cp, cd, xyz, ps, seid, comment = self._setup(
                ncards, self.cards, 'int64', 'float64')

        self._save(node_id, cp, cd, xyz, ps, seid, comment)

        grdset = self.model.grdset
        if grdset is not None:
            cp_ = grdset.cp
            cd_ = grdset.cd
            ps_ = grdset.ps
            seid_ = grdset.seid
        else:
            cp_ = 0
            cd_ = 0
            ps_ = 0
            seid_ = 0

        GRID_CP_DEFAULT = GRID_CD_DEFAULT = GRID_PS_DEFAULT = GRID_SEID_DEFAULT = -100
        data_temp_default = [
            (self.cp, GRID_CP_DEFAULT, cp_),
            (self.cd, GRID_CD_DEFAULT, cd_),
            (self.ps, GRID_PS_DEFAULT, ps_),
            (self.seid, GRID_SEID_DEFAULT, seid_),
        ]
        for data, temp_value, default_value in data_temp_default:
            ibad = np.where(data == temp_value)[0]
            if len(ibad):
                data[ibad] = default_value

        assert self.cp.min() >= 0, self.cp.min()

        # -1 is for fluid
        # -3 is for ADAPT/PVAL/PINTS/SINT/GMBNDS
        assert self.cd.min() >= -3, self.cd.min()
        self.sort()
        icp0 = np.where(self.cp == 0)[0]
        self._xyz_cid0[icp0, :] = self.xyz[icp0, :]
        assert self.xyz.shape == self._xyz_cid0.shape
        self.cards = []

    def _setup(self, ncards: int, cards: list[Any],
               idtype: str, fdtype: str):
        node_id = np.zeros(ncards, dtype=idtype)
        cp = np.zeros(ncards, dtype=idtype)
        cd = np.zeros(ncards, dtype=idtype)
        xyz = np.zeros((ncards, 3), dtype=fdtype)
        ps = np.zeros(ncards, dtype=idtype)
        seid = np.zeros(ncards, dtype=idtype)

        comment = {}
        for icard, card in enumerate(cards):
            (nidi, xyzi, cpi, cdi, psi, seidi, commenti) = card
            node_id[icard] = nidi
            cp[icard] = cpi
            cd[icard] = cdi
            xyz[icard, :] = xyzi
            ps[icard] = psi
            seid[icard] = seidi
            if commenti:
                #comment[i] = commenti
                comment[nidi] = commenti
        return node_id, cp, cd, xyz, ps, seid, comment

    def _save(self,
              node_id: np.ndarray, cp: np.ndarray, cd: np.ndarray,
              xyz: np.ndarray, ps: np.ndarray, seid: np.ndarray,
              comment: Optional[dict[int, str]]=None) -> None:
        ncards = len(node_id)
        ncards_existing = len(self.node_id)

        _xyz_cid0 = np.full((ncards, 3), np.nan, dtype='float64')
        if ncards_existing != 0:
            node_id = np.hstack([self.node_id, node_id])
            cp = np.hstack([self.cp, cp])
            cd = np.hstack([self.cd, cd])
            xyz = np.vstack([self.xyz, xyz])
            ps = np.hstack([self.ps, ps])
            seid = np.hstack([self.seid, seid])
            _xyz_cid0 = np.vstack([self._xyz_cid0, _xyz_cid0])
        if comment:
            self.comment.update(comment)

        self.node_id = node_id
        self.cp = cp
        self.cd = cd
        self.xyz = xyz
        self.ps = ps
        self.seid = seid
        self._xyz_cid0 = _xyz_cid0
        self.n = len(node_id)
        #self.sort()
        #self.cards = []

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        #used_dict['node_id'].append(self.node_id)
        used_dict['coord_id'].append(self.cp)
        used_dict['coord_id'].append(self.cd)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        node_id = used_dict['node_id']
        ncards_removed = len(self.node_id) - len(node_id)
        if ncards_removed:
            self.slice_card_by_id(node_id, assume_sorted=True, sort_ids=False)
        return ncards_removed

    def convert(self, xyz_scale: float=1., **kwargs):
        self._xyz_cid0 *= xyz_scale
        if self.cp.max() == 0:
            self.xyz *= xyz_scale
            return
        coords = self.model.coord.slice_card_by_id(self.cp)
        is_xyz = (coords.icoord == 2) & (coords.coord_type == 'R')
        is_rtz = (coords.icoord == 2) & (coords.coord_type == 'C')
        is_rtp = (coords.icoord == 2) & (coords.coord_type == 'S')
        self.xyz[is_xyz, :] *= xyz_scale # xyz
        self.xyz[is_rtz, 0] *= xyz_scale # R
        self.xyz[is_rtz, 2] *= xyz_scale # z
        self.xyz[is_rtp, 0] *= xyz_scale # rho

    def has_ps(self) -> bool:
        return self.n > 0 and self.ps.max() > 0

    #def slice_by_node_id(self, node_id: np.ndarray) -> GRID:
        #inid = self._node_index(node_id)
        #return self.slice_card(inid)

    def slice_card_by_node_id(self, node_id: np.ndarray, sort_ids: bool=False) -> GRID:
        """uses a node_ids to extract GRIDs"""
        inid = self.index(node_id)
        #assert len(self.node_id) > 0, self.node_id
        #i = np.searchsorted(self.node_id, node_id)
        grid = self.slice_card_by_index(inid, sort_index=sort_ids)
        return grid

    def slice_card_by_index(self, i: np.ndarray, sort_index: bool=False) -> GRID:
        """uses a node_index to extract GRIDs"""
        assert self.xyz.shape == self._xyz_cid0.shape
        assert len(self.node_id) > 0, self.node_id
        i = np.atleast_1d(np.asarray(i, dtype=self.node_id.dtype))
        if sort_index:
            i.sort()
        grid = GRID(self.model)
        self.__apply_slice__(grid, i)
        return grid

    def __apply_slice__(self, grid: GRID, i: np.ndarray) -> None:
        self._slice_comment(grid, i)

        grid.n = len(i)
        grid._is_sorted = self._is_sorted
        grid.node_id = self.node_id[i]
        grid.cp = self.cp[i]
        grid.cd = self.cd[i]
        grid.ps = self.ps[i]
        grid.seid = self.seid[i]

        grid.xyz = self.xyz[i, :]
        grid._xyz_cid0 = self._xyz_cid0[i, :]
        assert grid.xyz.shape == grid._xyz_cid0.shape

    def get_cp_transform(self) -> dict[int, list[int]]:
        nid_cp_transform = defaultdict(list)
        icp = np.where(self.cp != 0)[0]
        for nid, cp in zip(self.node_id[icp], self.cp[icp]):
            nid_cp_transform[cp].append(nid)
        return dict(nid_cp_transform)

    def get_cd_transform(self) -> dict[int, list[int]]:
        nid_cd_transform = defaultdict(list)
        icd = np.where(self.cd != 0)[0]
        for nid, cd in zip(self.node_id[icd], self.cd[icd]):
            nid_cd_transform[cd].append(nid)
        return dict(nid_cd_transform)

    def get_nid_cp_cd_xyz_transforms(self):
        nnodes = len(self.node_id)
        nid_cp_cd = np.zeros((nnodes, 3), dtype=self.node_id.dtype)
        nid_cp_cd[:, 0] = self.node_id
        nid_cp_cd[:, 1] = self.cp
        nid_cp_cd[:, 2] = self.cd

        #for nid, node in sorted(model.nodes.items()):
            #cd = node.Cd()
            #cp = node.Cp()
            #nids_cp_transform[cp].append(nid)
            #nids_cd_transform[cd].append(nid)
            #nid_cp_cd[i, :] = [nid, cp, cd]
            #xyz_cp[i, :] = node.xyz
            #i += 1
        nids_cp_transform = self.get_cp_transform()
        nids_cd_transform = self.get_cd_transform()
        return nids_cp_transform, nids_cd_transform, nid_cp_cd, self.xyz

    @property
    def _coord_basic(self):
        coord_basic = self.cp.max() == 0 and self.cp.min() == 0 and self.cd.max() == 0 and self.cd.min() == 0
        return coord_basic

    @parse_check
    def write(self, size: int=8, is_double: bool=False,
              write_card_header: bool=False) -> str:
        stringio = StringIO()
        self.write_file(stringio, size=size, is_double=is_double, write_card_header=write_card_header)
        msg = stringio.getvalue()
        return msg

    @property
    def max_id(self) -> int:
        max_int = max(
            self.node_id.max(),
            self.cp.max(),
            self.cd.max(), )
        return max_int

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        size = update_field_size(self.max_id, size)
        if size == 8:
            self.write_file_8(bdf_file, write_card_header=write_card_header)
        elif is_double:
            self.write_file_double(bdf_file, write_card_header=write_card_header)
        else:
            self.write_file_16(bdf_file, write_card_header=write_card_header)
        return

    @parse_check
    def write_file_8(self, bdf_file: TextIOLike,
                     write_card_header: bool=False) -> None:
        coord_basic, extra_basic, is_basic = self._write_flags()
        node_id = array_str(self.node_id, size=8)
        xyzs = array_float(self.xyz, size=8)
        if is_basic:
            cps = ''
            for nid, xyz in zip_longest(node_id, xyzs):
                #x, y, z = xyz
                comment = self.comment.get(nid, '')
                bdf_file.write(comment)
                msg = 'GRID    %8s%8s%8s%8s%8s\n' % (
                    nid, cps, xyz[0], xyz[1], xyz[2])
                bdf_file.write(comment)
                bdf_file.write(msg)
        elif extra_basic: # cp=cd != 0, but seid=ps=0
            cps = array_default_int(self.cp, size=8, default=0)
            cds = array_default_int(self.cd, size=8, default=0)
            for nid, cp, cd, xyz in zip_longest(node_id, cps, cds, xyzs):
                #x, y, z = xyz
                #msg = print_card_8(['GRID', nid, cp, x, y, z, cd])
                comment = self.comment.get(nid, '')
                msg = 'GRID    %8s%8s%8s%8s%8s%8s\n' % (
                    nid, cp, xyz[0], xyz[1], xyz[2], cd)
                bdf_file.write(comment)
                bdf_file.write(msg)
        else:
            # ps=seid != 0
            cps = array_default_int(self.cp, size=8, default=0)
            cds = array_default_int(self.cd, size=8, default=0)
            pss = array_default_int(self.ps, size=8, default=0)
            seids = array_default_int(self.seid, size=8, default=0)
            for nid, cp, cd, xyz, ps, seid in zip_longest(node_id, cps, cds, xyzs, pss, seids):
                #msg = print_card_8(['GRID', nid, cp, x, y, z, cd])
                comment = self.comment.get(nid, '')
                msg = 'GRID    %8s%8s%8s%8s%8s%8s%8s%8s\n' % (
                    nid, cp, xyz[0], xyz[1], xyz[2], cd, ps, seid)
                bdf_file.write(comment)
                bdf_file.write(msg)
        return

    def _write_flags(self):
        coord_basic = self._coord_basic
        extra_basic = self.ps.max() == 0 and self.ps.min() == 0 and self.seid.max() == self.seid.min()
        is_basic = coord_basic and extra_basic
        return coord_basic, extra_basic, is_basic

    @parse_check
    def write_file_16(self, bdf_file: TextIOLike,
                      write_card_header: bool=False) -> None:
        return _write_grid_large(self, bdf_file, print_scientific_16,
                                 is_double=False, write_card_header=write_card_header)

    @parse_check
    def write_file_double(self, bdf_file: TextIOLike,
                          write_card_header: bool=False) -> None:
        return _write_grid_large(self, bdf_file, print_scientific_double,
                                 is_double=True, write_card_header=write_card_header)

    def sort2(self) -> None:  #  pragma: no cover
        r"""
        test_bdfv C:\MSC.Software\msc_nastran_docs_2020\tpl6\avl\ship_hull103.dat --skip_mass --skip_nominal

        large problem, but it's tricky...
        """
        # get the sort index
        ids = self.node_id

        isort = np.argsort(ids)
        #iarg_array = np.arange(len(iarg), dtype='int32')

        counter_counts = np.array(list(Counter(ids).values()))
        idup = np.where(counter_counts > 1)[0]

        ids_sorted = ids[isort]
        uids_sorted, iinv, counts = np.unique(ids, return_inverse=True, return_counts=True)

        if len(ids) == len(uids_sorted):
            # no duplicate nodes
            self.__apply_slice__(self, isort)
        else:
            # duplicate nodes
            iunique = np.where(counts==1)[0]
            iduplicate = np.where(counts>1)[0]

            dnode = ids_sorted[1:] - ids_sorted[:-1]
            np.ma.count()
        self._is_sorted = True

    def sort(self) -> None:
        """sorts the card by node_id"""
        sort_duplicates(self)
        self._is_sorted = True

    #def get_global_position_by_node_id(self, nids: np.ndarray) -> np.ndarray:
        #xyz_cid0 = self.xyz_cid0()
        #inid = np.searchsorted(self.node_id, nids)
        #return xyz_cid0[inid, :]

    def get_position_by_node_id(self, node_id: np.ndarray) -> np.ndarray:
        """returns the global position of the node"""
        xyz_cid0 = self._xyz_cid0
        inid = self.index(node_id)
        return xyz_cid0[inid, :]

    #def index(self, node_id: np.ndarray, safe: bool=False) -> np.ndarray:
        #assert len(self.node_id) > 0, self.node_id
        #node_id = np.atleast_1d(np.asarray(node_id, dtype=self.node_id.dtype))
        #inid = np.searchsorted(self.node_id, node_id)
        #if safe:
            #ibad = inid >= len(self.node_id)
            #if sum(ibad):
                ##self.model.log.error(f'bad nids; node_id={node_id[ibad]}')
                #raise RuntimeError(f'bad nids; node_id={node_id[ibad]}')
            #inids_leftover = inid[~ibad]
            #if len(inids_leftover):
                #actual_nids = self.node_id[inids_leftover]
                #assert np.array_equal(actual_nids, node_id)
        #return inid

    def xyz_cid0(self) -> np.ndarray:
        """not validated"""
        if not self._is_sorted:
            self.sort()
        xyz_cid0 = self._xyz_cid0
        assert xyz_cid0.shape[0] > 0, xyz_cid0.shape

        xyz = self.xyz
        #self.model.log.info(f'xyz = {xyz_cid0.shape}')
        assert xyz.shape == xyz_cid0.shape, f'xyz.shape={xyz.shape} xyz_cid0.shape={xyz_cid0.shape}'
        #icp = (self.cp > 0)
        #cps_to_transform = self.cp[icp]
        # TODO: filter?
        ucp = np.unique(self.cp)

        coord = self.model.coord
        coord_types = coord.coord_type
        coord_ids = coord.coord_id
        origins = coord.origin
        iaxes = coord.i
        jaxes = coord.j
        kaxes = coord.k
        for cp in ucp:
            if cp == 0:
                continue
            icp = np.where(self.cp == cp)[0]
            #nnodes = len(icp)

            icid = np.where(coord_ids == cp)[0]
            assert len(icid) == 1, f'Coord={cp} not found; referenced by CP field of GRID'
            #coord_type = np.asscalar(coord_types[icid])  #deprecated in numpy 1.16
            coord_type = coord_types[icid].item()

            #assert coord_type in {'R'}, coord_type
            assert coord_type in {'R', 'C', 'S'}, coord_type  # output coord
            origin = origins[icid, :]
            iaxis = iaxes[icid, :]
            jaxis = jaxes[icid, :]
            kaxis = kaxes[icid, :]
            beta = np.vstack([iaxis, jaxis, kaxis])
            assert origin.shape == (1, 3), origin.shape
            #print('--------------------------')
            #print(f'cp = {cp}')
            #print(f'beta = {beta}')
            #print(f'origin = {origin}')
            if 0:
                xyz_cid = xyz[icp, :] @ beta + origin
                if coord_type == 'R':
                    xyz_cid0i = xyz_cid
                elif coord_type == 'C':
                    xyz_cid0i = cylindrical_to_rectangular(xyz_cid)
                elif coord_type == 'S':
                    xyz_cid0i = spherical_to_rectangular(xyz_cid)
                else:
                    raise NotImplementedError(coord_type)
            else:
                # transform node from CORD2C into a rectangular frame
                # then transform it to the global
                if coord_type == 'R':
                    xyz_cidr = xyz[icp, :]
                elif coord_type == 'C':
                    xyz_cidr = cylindrical_to_rectangular(xyz[icp, :])
                elif coord_type == 'S':
                    xyz_cidr = spherical_to_rectangular(xyz[icp, :])
                else:
                    raise NotImplementedError(coord_type)
                xyz_cid0i = xyz_cidr @ beta + origin
            xyz_cid0[icp, :] = xyz_cid0i
        assert not np.any(np.isnan(xyz_cid0))
        assert xyz_cid0.shape[0] > 0, xyz_cid0.shape
        return xyz_cid0


class POINT(VectorizedBaseCard):
    """
    +-------+-----+----+----+----+----+
    |   1   |  2  | 3  | 4  | 5  | 6  |
    +=======+=====+====+====+====+====+
    | POINT | NID | CP | X1 | X2 | X3 |
    +-------+-----+----+----+----+----+

    """
    _id_name = 'point_id'
    def clear(self) -> None:
        self._is_sorted = False
        self.n = 0
        self.point_id = np.array([], dtype='int32')
        self.cp = np.array([], dtype='int32')
        self.xyz = np.zeros((0, 3), dtype='float64')
        self._xyz_cid0 = np.zeros((0, 3), dtype='float64')

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['point_id'].append(self.point_id)
        used_dict['coord_id'].append(self.cp)

    def convert(self, xyz_scale: float=1., **kwargs):
        self._xyz_cid0 *= xyz_scale
        if self.cp.max() == 0:
            self.xyz *= xyz_scale
            return
        coords = self.model.coord.slice_card_by_id(self.cp)
        is_xyz = (coords.icoord == 2) & (coords.coord_type == 'R')
        is_rtz = (coords.icoord == 2) & (coords.coord_type == 'C')
        is_rtp = (coords.icoord == 2) & (coords.coord_type == 'S')
        self.xyz[is_xyz, :] *= xyz_scale # xyz
        self.xyz[is_rtz, 0] *= xyz_scale # R
        self.xyz[is_rtz, 2] *= xyz_scale # z
        self.xyz[is_rtp, 0] *= xyz_scale # rho

    def add(self, nid: int, xyz: np.ndarray,
            cp: int=0, comment: str='') -> int:
        self.cards.append((nid, xyz, cp, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')
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
        self.cards.append((nid, xyz, cp, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug('parse POINT')

        point_id = np.zeros(ncards, dtype='int32')
        cp = np.zeros(ncards, dtype='int32')
        xyz = np.zeros((ncards, 3), dtype='float64')
        _xyz_cid0 = np.full((ncards, 3), np.nan, dtype='float64')
        comment = self.comment
        for icard, card in enumerate(self.cards):
            (nid, xyzi, cpi, commenti) = card
            point_id[icard] = nid
            cp[icard] = cpi
            xyz[icard, :] = xyzi
            comment[icard] = commenti

        self._save(point_id, cp, xyz, _xyz_cid0, comment)
        self.sort()
        icp0 = np.where(self.cp == 0)[0]
        self._xyz_cid0[icp0, :] = self.xyz[icp0, :]
        assert self.xyz.shape == self._xyz_cid0.shape
        self.cards = []

    def _save(self, point_id, cp, xyz, _xyz_cid0, comment=None):
        if len(self.point_id) > 0:
            asdf
        self.point_id = point_id
        self.cp = cp
        self.xyz = xyz
        self._xyz_cid0 = _xyz_cid0
        if comment:
            self.comment = comment

    #def slice_by_node_id(self, node_id: np.ndarray) -> GRID:
        #inid = self._node_index(node_id)
        #return self.slice_card(inid)

    def slice_card_by_point_id(self, node_id: np.ndarray) -> POINT:
        """uses a node_ids to extract POINTs"""
        inid = self.index(node_id)
        #assert len(self.node_id) > 0, self.node_id
        #i = np.searchsorted(self.node_id, node_id)
        point = self.slice_card_by_index(inid)
        return point

    def __apply_slice__(self, point: POINT, i: np.ndarray) -> None:
        point.point_id = self.point_id[i]
        point.cp = self.cp[i]
        point.xyz = self.xyz[i, :]
        point._xyz_cid0 = self._xyz_cid0[i, :]
        point.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        #nid = self.model.grid.node_id
        cids = self.model.coord.coord_id
        geom_check(self,
                   missing,
                   property_id=(cids, self.cp))

    @property
    def max_id(self) -> int:
        return self.point_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        point_ids = array_str(self.point_id, size=size)
        cps = array_default_int(self.cp, default=0, size=size)
        xyzs = array_float(self.xyz, size=size, is_double=False).tolist()
        for point_id, cp, xyz in zip(point_ids, cps, xyzs):
            list_fields = ['POINT', point_id, cp] + xyz
            bdf_file.write(print_card(list_fields))
        return

    #def sort(self) -> None:
        #i = np.argsort(self.point_id)
        #self.__apply_slice__(self, i)

    #def index(self, node_id: np.ndarray) -> np.ndarray:
        #assert len(self.node_id) > 0, self.node_id
        #node_id = np.atleast_1d(np.asarray(node_id, dtype=self.node_id.dtype))
        #inid = np.searchsorted(self.node_id, node_id)
        #return inid


def _write_grid_large(grid: GRID, bdf_file: TextIOLike,
                      print_scientific: Callable[float],
                      is_double: bool=False,
                      write_card_header: bool=False) -> None:
    coord_basic, extra_basic, is_basic = grid._write_flags()
    node_id = array_str(grid.node_id, size=16)
    xyzs = array_float(grid.xyz, size=16, is_double=is_double)
    if is_basic:
        for nid, xyz in zip_longest(node_id, xyzs):
            comment = grid.comment.get(nid, '')
            msg = ('GRID*   %16s                %16s%16s\n'
                   '*       %16s\n' % (
                       nid, xyz[0], xyz[1], xyz[2]))
            bdf_file.write(comment)
            bdf_file.write(msg)
    elif extra_basic: # cp=cd != 0, but seid=ps=0
        cps = array_default_int(grid.cp, size=16, default=0)
        cds = array_default_int(grid.cd, size=16, default=0)
        for nid, cp, cd, xyz in zip_longest(node_id, cps, cds, xyzs):
            comment = grid.comment.get(nid, '')
            msg = ('GRID*   %16s%16s%16s%16s\n'
                   '*       %16s%16s\n' % (
                       nid, cp, xyz[0], xyz[1], xyz[2], cd))
            bdf_file.write(comment)
            bdf_file.write(msg)
    else:
        # seid=ps != 0
        cps = array_default_int(grid.cp, size=16, default=0)
        cds = array_default_int(grid.cd, size=16, default=0)
        pss = array_default_int(grid.ps, size=16, default=0)
        seids = array_default_int(grid.seid, size=16, default=0)
        for nid, cp, cd, xyz, ps, seid in zip_longest(node_id, cps, cds, xyzs, pss, seids):
            comment = grid.comment.get(nid, '')
            msg = ('GRID*   %16s%16s%16s%16s\n'
                   '*       %16s%16s%16s%16s\n' % (
                       nid, cp, xyz[0], xyz[1], xyz[2], cd, ps, seid))
            bdf_file.write(comment)
            bdf_file.write(msg)
    return
