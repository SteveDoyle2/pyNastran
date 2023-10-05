from __future__ import annotations
from itertools import zip_longest
from collections import defaultdict
from typing import Any, TYPE_CHECKING
import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru, _format_comment
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, # components_or_blank,
    integer_string_or_blank,
    components_or_blank as fcomponents_or_blank)

#from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int
from pyNastran.dev.bdf_vectorized3.cards.base_card import VectorizedBaseCard, parse_node_check #get_print_card_8_16,
from pyNastran.dev.bdf_vectorized3.cards.write_utils import get_print_card, MAX_8_CHAR_INT

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


class ABCOQSET(VectorizedBaseCard):
    """
    Defines degrees-of-freedom in the analysis set (A-set).

    +-------+-----+-----+------+------+-----+-----+-----+-----+
    |  1    |  2  | 3   |  4   |  5   |  6  |  7  |  8  | 9   |
    +=======+=====+=====+======+======+=====+=====+=====+=====+
    | ASET  | ID1 | C1  | ID2  |  C2  | ID3 | C3  | ID4 | C4  |
    +-------+-----+-----+------+------+-----+-----+-----+-----+
    | ASET  | 16  |  2  |  23  | 3516 |  1  |  4  |     |     |
    +-------+-----+-----+------+------+-----+-----+-----+-----+

    +-------+-----+-----+------+------+-----+-----+-----+-----+
    |   1   |  2  |  3  |   4  |  5   |  6  |  7  |  8  |  9  |
    +=======+=====+=====+======+======+=====+=====+=====+=====+
    | ASET1 |  C  | ID1 |  ID2 | ID3  | ID4 | ID5 | ID6 | ID7 |
    +-------+-----+-----+------+------+-----+-----+-----+-----+
    |       | ID8 | ID9 |      |      |     |     |     |     |
    +-------+-----+-----+------+------+-----+-----+-----+-----+
    | ASET1 |  C  | ID1 | THRU | ID2  |     |     |     |     |
    +-------+-----+-----+------+------+-----+-----+-----+-----+
    """
    _id_name = 'node_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        #self._is_sorted = False
        self.component = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')

    def add_set(self, nids: list[int], components: list[int],
                 comment: str='') -> int:
        nnodes = len(nids)
        ncomp = len(components)
        assert nnodes == ncomp, (nnodes, ncomp)
        self.cards.append((nids, components, comment))
        #if comment:
            #self.comment[nid] = _format_comment(comment)
        self.n += nnodes
        return self.n

    def add_set1(self, component: int, nids: list[int],
                  comment: str='') -> int:
        nnodes = len(nids)
        assert isinstance(component, integer_types), component
        components = [component] * nnodes
        self.cards.append((nids, components, comment))
        #if comment:
            #self.comment[nid] = _format_comment(comment)
        self.n += nnodes
        return self.n

    def add_set1_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        components_str = fcomponents_or_blank(card, 1, 'components', '0')
        component = int(components_str)

        nfields = len(card)
        ids = []
        i = 1
        for ifield in range(2, nfields):
            idi = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        ids = expand_thru(ids, set_fields=True, sort_fields=True)
        components = [component] * len(ids)
        #return cls(ids, components, comment=comment)

        self.cards.append((ids, components, comment))
        #if comment:
            #self.comment[nid] = comment
        self.n += len(ids)
        return self.n

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        try:
            node_id, component = self._setup(ncards, self.cards, 'int32')
        except OverflowError:
            node_id, component = self._setup(ncards, self.cards, 'int64')
        self._save(node_id, component)

        self.sort()
        self.cards = []

    def _setup(self, ncards: int, cards: list[Any],
               idtype: str) -> tuple[np.ndarray, np.ndarray]:
        node_id = []
        component_list = []
        #comment = {}
        for i, card in enumerate(cards):
            (nidi, componenti, commenti) = card
            assert isinstance
            node_id.extend(nidi)
            component_list.extend(componenti)
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti
        node_id2 = np.array(node_id, dtype=idtype)
        component2 = np.array(component_list, dtype=idtype)
        return node_id2, component2

    def _save(self,
              node_id: np.ndarray,
              component: np.ndarray,
              comment: dict[int, str]=None) -> None:
        #ncards = len(node_id)
        ncards_existing = len(self.node_id)

        if ncards_existing != 0:
            node_id = np.hstack([self.node_id, node_id])
            component = np.hstack([self.component, component])
        #if comment:
            #self.comment.update(comment)
        self.node_id = node_id
        self.component = component
        self.n = len(node_id)
        #self.sort()
        #self.cards = []

    #def slice_by_node_id(self, node_id: np.ndarray) -> GRID:
        #inid = self._node_index(node_id)
        #return self.slice_card(inid)

    #def slice_card_by_node_id(self, node_id: np.ndarray) -> GRID:
        #"""uses a node_ids to extract GRIDs"""
        #inid = self.index(node_id)
        ##assert len(self.node_id) > 0, self.node_id
        ##i = np.searchsorted(self.node_id, node_id)
        #grid = self.slice_card_by_index(inid)
        #return grid

    #def slice_card_by_index(self, i: np.ndarray) -> GRID:
        #"""uses a node_index to extract GRIDs"""
        #assert self.xyz.shape == self._xyz_cid0.shape
        #assert len(self.node_id) > 0, self.node_id
        #i = np.atleast_1d(np.asarray(i, dtype=self.node_id.dtype))
        #i.sort()
        #grid = GRID(self.model)
        #self.__apply_slice__(grid, i)
        #return grid

    def __apply_slice__(self, grid: ASET, i: np.ndarray) -> None:
        self._slice_comment(grid, i)
        grid.n = len(i)
        grid.node_id = self.node_id[i]
        grid.component = self.component[i]

    @parse_node_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        max_int = self.node_id.max()
        #size = update_field_size(max_int, size)
        print_card = print_card_8

        comp_to_nids = defaultdict(list)
        for nid, comp in zip_longest(self.node_id, self.component):
            comp_to_nids[comp].append(nid)

        #bdf_file.write(comment)
        if self.type in {'ASET', 'BSET', 'CSET', 'QSET',
                         'ASET1', 'BSET1', 'CSET1', 'QSET1'}:
            class_name = self.type[0] + 'SET1'
        elif self.type in {'OMIT', 'OMIT1'}:
            class_name = 'OMIT1'
        else:
            raise NotImplementedError(self.type)

        for comp, nids in comp_to_nids.items():
            node_id = array_str(np.array(nids), size=size).tolist()
            list_fields = [class_name, comp, ] + node_id
            bdf_file.write(print_card(list_fields))
        return

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

class ASET(ABCOQSET):
    pass
class BSET(ABCOQSET):
    pass
class CSET(ABCOQSET):
    pass
class QSET(ABCOQSET):
    pass
class OMIT(ABCOQSET):
    pass

#BSET = CSET = QSET = OMIT = ASET
