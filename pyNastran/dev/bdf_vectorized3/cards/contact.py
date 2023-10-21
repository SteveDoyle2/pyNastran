from __future__ import annotations
from typing import Union, TYPE_CHECKING
from collections import defaultdict
from itertools import count
import numpy as np

from pyNastran.utils.numpy_utils import integer_types, float_types
#from pyNastran.bdf import MAX_INT
from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.cards.collpase_card import collapse_thru, collapse_thru_packs
from pyNastran.dev.bdf_vectorized3.cards.base_card import VectorizedBaseCard, hslice_by_idim, make_idim
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank,
    components_or_blank, parse_components,)
#from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16 # print_float_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.dev.bdf_vectorized3.cards.base_card import get_print_card_8_16
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    #array_default_str,
    array_str, array_default_int, update_field_size)
from pyNastran.dev.bdf_vectorized3.utils import print_card_8

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF




class ElementSet(VectorizedBaseCard):
    def clear(self):
        self.sid = np.array([], dtype='int32')
        self.nelement = np.array([], dtype='int32')
        self.element_id = np.array([], dtype='int32')

    def add(self, sid: list[int], element_ids: list[int],
            comment: str='') -> int:
        self.cards.append((sid, element_ids, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        sid = integer(card, 1, 'sid')
        nfields = len(card) - 1
        nrows = (nfields // 8)
        if (nfields % 8) > 0:
            nrows += 1

        ifield = 1
        element_ids_list = []
        for irow in range(nrows):
            if irow == 0:
                fieldsi = card[ifield+1:ifield+8]
            else:
                fieldsi = card[ifield:ifield+8]

            #print(f'{irow:d}/{nrows:d} = {fieldsi}')
            if 'THRU' in fieldsi:
                fieldsi = expand_thru(fieldsi, set_fields=True, sort_fields=False)
                element_ids_list += fieldsi
            else:
                for fieldi in fieldsi:
                    if fieldi is None:
                        continue
                    element_ids_list.append(fieldi)
            ifield += 8

        idtype = self.model.idtype
        element_ids = np.array(element_ids_list, dtype=idtype)
        assert len(card) > 2, f'len({self.type} card) = {len(card):d}\ncard={card}'
        assert 1100 not in element_ids
        self.cards.append((sid, element_ids, comment))
        self.n += 1
        return self.n

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        #ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        ncards = len(self.cards)
        idtype = self.model.idtype
        sid = np.zeros(ncards, dtype=idtype)
        nelement = np.zeros(ncards, dtype=idtype)
        element_ids_list = []
        #comment = {}
        for icard, card in enumerate(self.cards):
            (sidi, element_idsi, commenti) = card
            assert isinstance(sidi, int), sidi
            #assert isinstance(element_ids, (list), componenti
            sid[icard] = sidi
            nelement[icard] = len(element_idsi)
            element_ids_list.extend(element_idsi)
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti
        sid = np.array(sid, dtype=idtype)
        element_ids = np.array(element_ids_list, dtype=idtype)
        self._save(sid, nelement, element_ids, comment=None)
        #self.sort()
        self.cards = []

    def _save(self,
              sid: np.ndarray,
              nelement: np.ndarray,
              element_ids: np.ndarray,
              comment: dict[int, str]=None) -> None:
        ncards_existing = len(self.sid)
        if ncards_existing != 0:
            sid = np.hstack([self.sid, sid])
            nelement = np.hstack([self.nelement, nelement])
            element_ids = np.hstack([self.element_ids, element_ids])
        #if comment:
            #self.comment.update(comment)
        self.sid = sid
        self.nelement = nelement
        self.element_ids = element_ids
        self.n = len(self.sid)
        #self.sort()
        #self.cards = []

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['element_id'].append(self.element_id)

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

    def __apply_slice__(self, eset: ElementSet, i: np.ndarray) -> None:
        self._slice_comment(eset, i)
        eset.n = len(i)
        eset.sid = self.sid[i]
        eset.nelement = self.nelement[i]
        eset.element_id = self.element_id[i]

    @property
    def ielement(self) -> np.ndarray:
        return make_idim(self.n, self.nelement)

    @property
    def max_id(self) -> int:
        return self.element_ids.max()

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if self.n == 0:
            return
        size = update_field_size(self.max_id, size)
        print_card = get_print_card_8_16(size)

        set_to_eids = defaultdict(list)
        for sid, (ielement0, ielement1) in zip(self.sid, self.ielement):
            eids = self.element_ids[ielement0:ielement1].tolist()
            set_to_eids[sid].append(eids)

        for sid, eids_list in set_to_eids.items():
            eids = np.hstack(eids_list).tolist()
            singles, doubles = collapse_thru_packs(eids)
            #print('singles/doubles', singles, doubles)

            if len(singles):
                list_fields = [self.type, sid, ] + singles
                bdf_file.write(print_card(list_fields))

            for double in doubles:
                list_fields = [self.type, sid] + double
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

class BSURF(ElementSet):
    def geom_check(self, missing: dict[str, np.ndarray]):
        eid = self.model.shell_element_ids
        geom_check(self,
                   missing,
                   element_id=(eid, self.element_id),)


class BSURFS(ElementSet):
    def geom_check(self, missing: dict[str, np.ndarray]):
        eid = self.model.solid_element_ids
        geom_check(self,
                   missing,
                   element_id=(eid, self.element_id),)

