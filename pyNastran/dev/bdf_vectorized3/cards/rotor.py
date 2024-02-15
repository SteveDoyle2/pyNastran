from __future__ import annotations
from itertools import zip_longest, count
from collections import defaultdict
from typing import Any, TYPE_CHECKING
import numpy as np
#from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import read_ids_thru, expand_thru, _format_comment
from pyNastran.bdf.cards.collpase_card import collapse_thru, collapse_thru_packs
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, # string, parse_components,
    #integer_or_string, fields,
    #integer_or_blank, double_or_blank, # components_or_blank,
    #integer_string_or_blank, # parse_components_or_blank,
    #components_or_blank as fcomponents_or_blank,
)

#from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, parse_check,
    hslice_by_idim, make_idim, remove_unused_duplicate)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    get_print_card_size, array_str, array_default_int)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class ROTORG(VectorizedBaseCard):
    """
    Rotor Line Model Grids

    Specifies grids that compose the rotor line model.

    +--------+-----------+-----+------+------+-----+-----+-----+-----+
    |    1   |     2     |  3  |   4  |  5   |  6  |  7  |  8  |  9  |
    +========+===========+=====+======+======+=====+=====+=====+=====+
    | ROTORG |  ROTOR_ID | ID1 |  ID2 | ID3  | ID4 | ID5 | ID6 | ID7 |
    +--------+-----------+-----+------+------+-----+-----+-----+-----+
    | ROTORG |    ID8    | ID9 |      |      |     |     |     |     |
    +--------+-----------+-----+------+------+-----+-----+-----+-----+
    """
    _id_name = 'rotor_id'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.rotor_id = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')

    def add(self, rotor_id: int, nids: list[int], comment: str='') -> int:
        """Creates a ROTORG card"""
        nids = expand_thru(nids, set_fields=True, sort_fields=False)
        nnodes = len(nids)
        self.cards.append((rotor_id, nids, comment))
        #if comment:
            #self.comment[nid] = _format_comment(comment)
        self.n += nnodes
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        rotor_id = integer(card, 1, 'rotor_id')
        ids = read_ids_thru(card, 2, base_str='node_id%d')
        nids = expand_thru(ids)

        self.cards.append((rotor_id, nids, comment))
        self.n += len(ids)
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        idtype = self.model.idtype
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        rotor_id = np.zeros(ncards, dtype='int32')
        nnode = np.zeros(ncards, dtype='int32')
        node_id_list = []
        #comment = {}
        for icard, card in enumerate(self.cards):
            (rotor_idi, nidi, commenti) = card
            assert isinstance(nidi, list), nidi
            rotor_id[icard] = rotor_idi
            nnodei = len(nidi)
            assert nnodei > 0, (seidi, nidi)
            nnode[icard] = nnodei
            node_id_list.extend(nidi)
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti

        node_id = np.array(node_id_list, dtype=idtype)
        self._save(rotor_id, node_id, nnode)
        self.sort()
        self.cards = []

    def _save(self,
              rotor_id: np.ndarray,
              node_id: np.ndarray,
              nnode: np.ndarray,
              comment: dict[int, str]=None) -> None:
        #ncards = len(node_id)
        ncards_existing = len(self.node_id)

        if ncards_existing != 0:
            rotor_id = np.hstack([self.rotor_id, rotor_id])
            nnode = np.hstack([self.nnode, nnode])
            node_id = np.hstack([self.node_id, node_id])
        #if comment:
            #self.comment.update(comment)

        assert len(rotor_id) == len(nnode)
        assert nnode.min() > 0, nnode
        self.rotor_id = rotor_id
        self.nnode = nnode
        self.node_id = node_id
        self.n = len(rotor_id)
        #self.cards = []

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.node_id
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['node_id'].append(self.node_id)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        pass

    @property
    def inode(self) -> np.ndarray:
        return make_idim(self.n, self.nnode)

    def __apply_slice__(self, rotor: ROTORG, i: np.ndarray) -> None:
        self._slice_comment(rotor, i)
        rotor.n = len(i)
        rotor.rotor_id = self.rotor_id[i]

        inode = self.inode
        rotor.node_id = hslice_by_idim(i, inode, self.node_id)
        rotor.nnode = self.nnode[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.node_id),)

    @property
    def max_id(self) -> int:
        return max(self.rotor_id.max(), self.node_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        rotor_to_nids = defaultdict(list)
        for rotor_idi, (inode0, inode1) in zip_longest(self.rotor_id, self.inode):
            nids = self.node_id[inode0:inode1].tolist()
            rotor_to_nids[rotor_idi].extend(nids)

        for rotor_idi, nids_list in rotor_to_nids.items():
            # TODO: collapse
            list_fields = ['ROTORG', rotor_idi] + nids_list
            bdf_file.write(print_card(list_fields))
        return
