from __future__ import annotations
from typing import TYPE_CHECKING
from collections import defaultdict
#from itertools import count
import numpy as np

#from pyNastran.utils.numpy_utils import integer_types, float_types
#from pyNastran.bdf import MAX_INT
from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.cards.collpase_card import collapse_thru, collapse_thru_packs
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, hslice_by_idim, vslice_by_idim, make_idim)
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double,
    integer_or_blank, double_or_blank, string_or_blank,
    #components_or_blank, parse_components,
)
from pyNastran.bdf.bdf_interface.assign_type_force import lax_double_or_blank
#from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16 # print_float_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.dev.bdf_vectorized3.cards.base_card import get_print_card_8_16
from .constraints import ADD
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    #array_default_str,
    array_str, array_float,
    array_default_int, array_default_float,
    get_print_card_size)
from pyNastran.dev.bdf_vectorized3.utils import print_card_8

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


class ElementPropertyNodeSet(VectorizedBaseCard):
    def clear(self):
        self.sid = np.array([], dtype='int32')
        self.n_ids = np.array([], dtype='int32')
        self.ids = np.array([], dtype='int32')

    def add(self, sid: list[int], ids: list[int],
            comment: str='') -> int:
        self.cards.append((sid, ids, comment))
        self.n += 1
        return self.n - 1

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

        idtype = self.model.idtype
        sid = integer(card, 1, 'sid')
        nfields = len(card) - 1
        nrows = (nfields // 8)
        if (nfields % 8) > 0:
            nrows += 1

        ifield = 1
        ids_list = []
        for irow in range(nrows):
            if irow == 0:
                fieldsi = card[ifield+1:ifield+8]
            else:
                fieldsi = card[ifield:ifield+8]
            # need to make sure THRU is capitalized
            fieldsi = [field.upper() if isinstance(field, str) else field for field in fieldsi]

            #print(f'{irow:d}/{nrows:d} = {fieldsi}')
            if 'THRU' in fieldsi:
                fieldsi = expand_thru(fieldsi, set_fields=True, sort_fields=False)
                ids_list += fieldsi
            else:
                for fieldi in fieldsi:
                    if fieldi is None:
                        continue
                    ids_list.append(fieldi)
            ifield += 8

        ids = np.array(ids_list, dtype=idtype)
        assert len(card) > 2, f'len({self.type} card) = {len(card):d}\ncard={card}'
        self.cards.append((sid, ids, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        #ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        ncards = len(self.cards)
        idtype = self.model.idtype
        sid = np.zeros(ncards, dtype=idtype)
        n_ids = np.zeros(ncards, dtype=idtype)
        ids_list = []
        #comment = {}
        for icard, card in enumerate(self.cards):
            (sidi, idsi, commenti) = card
            assert isinstance(sidi, int), sidi
            #assert isinstance(element_ids, (list), componenti
            sid[icard] = sidi
            n_ids[icard] = len(idsi)
            ids_list.extend(idsi)
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti
        sid = np.array(sid, dtype=idtype)
        ids = np.array(ids_list, dtype=idtype)
        self._save(sid, n_ids, ids, comment=None)
        #self.sort()
        self.cards = []

    def _save(self,
              sid: np.ndarray,
              n_ids: np.ndarray,
              ids: np.ndarray,
              comment: dict[int, str]=None) -> None:
        ncards_existing = len(self.sid)
        if ncards_existing != 0:
            sid = np.hstack([self.sid, sid])
            n_ids = np.hstack([self.n_ids, n_ids])
            ids = np.hstack([self.ids, ids])
        #if comment:
            #self.comment.update(comment)
        self.sid = sid
        self.n_ids = n_ids
        self.ids = ids
        self.n = len(self.sid)
        #self.sort()
        #self.cards = []

    def __apply_slice__(self, eset: ElementPropertySet, i: np.ndarray) -> None:
        self._slice_comment(eset, i)
        eset.n = len(i)
        eset.sid = self.sid[i]
        eset.n_ids = self.n_ids[i]
        eset.ids = self.ids[i]
        raise RuntimeError(i)

    @property
    def max_id(self) -> int:
        return max(self.sid.max(), self.ids.max())

    @property
    def i_id(self) -> np.ndarray:
        return make_idim(self.n, self.n_ids)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if self.n == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)

        set_to_ids = defaultdict(list)
        for sid, (i0, i1) in zip(self.sid, self.i_id):
            ids0 = self.ids[i0:i1].tolist()
            set_to_ids[sid].append(ids0)

        for sid, ids_list in set_to_ids.items():
            ids = np.hstack(ids_list).tolist()
            singles, doubles = collapse_thru_packs(ids)
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


class PropertySet(ElementPropertyNodeSet):
    @property
    def nproperty(self) -> np.ndarray:
        return self.n_ids

    @nproperty.setter
    def nproperty(self, nproperty: np.ndarray) -> None:
        self.n_ids = nproperty

    @property
    def property_id(self) -> np.ndarray:
        return self.ids

    @property_id.setter
    def property_id(self, property_id: np.ndarray) -> None:
        self.ids = property_id

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['property_id'].append(self.property_id)

    @property
    def iproperty(self) -> np.ndarray:
        return self.i_id


class ElementSet(ElementPropertyNodeSet):
    @property
    def nelement(self) -> np.ndarray:
        return self.n_ids

    @nelement.setter
    def nelement(self, nelement: np.ndarray) -> None:
        self.n_ids = nelement

    @property
    def element_id(self) -> np.ndarray:
        return self.ids

    @element_id.setter
    def element_id(self, element_id: np.ndarray) -> None:
        self.ids = element_id

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['element_id'].append(self.element_id)

    @property
    def ielement(self) -> np.ndarray:
        return self.i_id


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


class BCPROP(PropertySet):
    def geom_check(self, missing: dict[str, np.ndarray]):
        pid = self.model.shell_property_ids
        geom_check(self,
                   missing,
                   property_id=(pid, self.property_id),)


class BCPROPS(PropertySet):
    def geom_check(self, missing: dict[str, np.ndarray]):
        pid = self.model.solid_property_ids
        geom_check(self,
                   missing,
                   property_id=(pid, self.property_id),)


class NodeSet(ElementPropertyNodeSet):
    @property
    def nnode(self) -> np.ndarray:
        return self.n_ids

    @nnode.setter
    def nnode(self, nnode: np.ndarray) -> None:
        self.n_ids = nnode

    @property
    def node_id(self) -> np.ndarray:
        return self.ids

    @node_id.setter
    def node_id(self, node_id: np.ndarray) -> None:
        self.ids = node_id

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['node_id'].append(self.node_id)

    @property
    def inode(self) -> np.ndarray:
        return self.i_id

class BLSEG(NodeSet):
    """Defines a glue or contact edge region or a curve for slideline contact via grid numbers."""
    def geom_check(self, missing: dict[str, np.ndarray]):
        nids = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nids, self.node_id),)


class BGSET(VectorizedBaseCard):
    """
    Defines glued contact pairs.

    +-------+------+------+------+---------+----+------+------+----+
    |   1   |  2   |  3   |   4  |    5    | 6  |  7   |   8  |  9 |
    +=======+======+======+======+=========+====+======+======+====+
    | BGSET | GSID | SID1 | TID1 | SDIST1  |    | EXT1 |      |    |
    +-------+------+------+------+---------+----+------+------+----+
    |       |      | SID2 | TID2 | SDIST2  |    | EXT2 |      |    |
    +-------+------+------+------+---------+----+------+------+----+
    """
    def clear(self):
        #self.sid = np.array([], dtype='int32')
        #self.nelement = np.array([], dtype='int32')
        #self.element_id = np.array([], dtype='int32')
        #: GSID Glue set identification number. (Integer > 0)
        self.glue_id = np.array([], dtype='int32')
        #: SIDi Source region (contactor) identification number for contact pair i.
        #: (Integer > 0)
        self.source_ids = np.array([], dtype='int32')

        #: TIDi Target region identification number for contact pair i. (Integer > 0)
        self.target_ids = np.array([], dtype='int32')

        #: SDISTi Search distance for glue regions (Real); (Default=10.0)
        self.search_distance = np.array([], dtype='float64')

        #: EXTi Extension factor for target region (SOLs 402 and 601 only).
        self.extension = np.array([], dtype='float64')

    def add(self, sid: list[int], element_ids: list[int],
            comment: str='') -> int:
        self.cards.append((sid, element_ids, comment))
        self.n += 1
        return self.n - 1

    #def remove_unused(self)
    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a BGSET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        glue_id = integer(card, 1, 'glue_id')
        sids = []
        tids = []
        sdists = []
        exts = []

        nfields = card.nfields
        i = 2
        j = 1
        while i < nfields:
            #SIDi Source region identification number for glue pair i. (Integer > 0)
            #TIDi Target region identification number for glue pair i. (Integer > 0)
            #SDISTi Search distance for glue regions (Real); (Default=10.0)
            #EXTi Extension factor for target region (SOLs 402 and 601 only).

            sids.append(integer(card, i, 'sid%s' % j))
            tids.append(integer(card, i + 1, 'tid%s' % j))
            sdists.append(double_or_blank(card, i + 2, 'fric%s' % j, default=0.0))
            #if sol == 101:
            exts.append(double_or_blank(card, i + 4, 'mind%s' % j, default=0.0))
            #else:
                #exts.append(None)
            i += 8
            j += 1
        #return BGSET(glue_id, sids, tids, sdists, exts,
                     #comment=comment, sol=sol)

        assert len(card) > 2, f'len({self.type} card) = {len(card):d}\ncard={card}'
        self.cards.append((glue_id, sids, tids, sdists, exts, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        #ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        ncards = len(self.cards)
        idtype = self.model.idtype
        fdtype = self.model.fdtype
        glue_id = np.zeros(ncards, dtype='int32')
        nsource = np.zeros(ncards, dtype='int32')
        source_ids = []
        target_ids = []
        search_distances = []
        extensions = []

        #comment = {}
        for icard, card in enumerate(self.cards):
            (glue_idi, source_idi, target_idi, search_distancei, extensioni, commenti) = card

            glue_id[icard] = glue_idi
            source_ids.extend(source_idi)
            target_ids.extend(target_idi)
            search_distances.extend(search_distancei)
            extensions.extend(source_idi)

            nsource[icard] = len(source_idi)
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti

        source_region = np.array(source_ids, dtype=idtype)
        target_id = np.array(target_ids, dtype=idtype)
        search_distance = np.array(search_distances, dtype=fdtype)
        extension = np.array(extensions, dtype=fdtype)
        self._save(glue_id, nsource, source_region, target_id, search_distance, extension, comment=None)
        #self.sort()
        self.cards = []

    def _save(self,
              glue_id,
              nsource,
              source_ids,
              target_ids,
              search_distance,
              extension: np.ndarray,
              comment: dict[int, str]=None) -> None:
        ncards_existing = len(self.glue_id)
        if ncards_existing != 0:
            asdf
            glue_id = np.hstack([self.glue_id, glue_id])
            #nelement = np.hstack([self.nelement, nelement])
            #element_ids = np.hstack([self.element_ids, element_ids])
        #if comment:
            #self.comment.update(comment)
        self.glue_id = glue_id
        self.nsource = nsource
        self.source_ids = source_ids
        self.target_ids = target_ids
        self.search_distance = search_distance
        self.extension = extension
        self.n = len(self.glue_id)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        contact_regions = np.hstack([self.source_ids, self.target_ids])
        used_dict['contact_set_id'].append(contact_regions)

    def convert(self, xyz_scale: float=1.0, **kwargs) -> None:
        self.search_distance *= xyz_scale

    def __apply_slice__(self, bgset: BGSET, i: np.ndarray) -> None:
        self._slice_comment(bgset, i)
        bgset.n = len(i)

        isource = self.isource
        bgset.glue_id = self.glue_id[i]
        bgset.nsource = self.nsource[i]
        bgset.source_ids = hslice_by_idim(i, isource, self.source_ids)
        bgset.target_ids = hslice_by_idim(i, isource, self.target_ids)
        bgset.search_distance = hslice_by_idim(i, isource, self.search_distance)
        bgset.extension = hslice_by_idim(i, isource, self.extension)

    @property
    def isource(self) -> np.ndarray:
        return make_idim(self.n, self.nsource)

    @property
    def max_id(self) -> int:
        return max(self.glue_id.max(), self.source_ids.max(), self.target_ids.max())

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if self.n == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)

        for glue_id, (isource0, isource1) in zip(self.glue_id, self.isource):
            target_ids = self.target_ids[isource0:isource1].tolist()
            source_ids = self.source_ids[isource0:isource1].tolist()
            extension = self.extension[isource0:isource1].tolist()
            search_distance = self.search_distance[isource0:isource1].tolist()
            list_fields = ['BGSET', glue_id]
            assert len(target_ids) > 0
            for source_id, target_id, ext, search_dist in zip(target_ids, source_ids, search_distance, extension):
                list_fields += [source_id, target_id, search_dist, None, ext, None, None, None]
            bdf_file.write(print_card(list_fields))
        return


class BCTSET(VectorizedBaseCard):
    """
    Defines contact pairs of a:
     - 2D (SOL 601 only)
     - 3D contact set (SOLs 101,103, 111, 112, 601 and 701)

    +--------+-------+------+-------+-------+-------+-------+-------+
    |   1    |   2   | 3    |  4    |   5   |   6   |   7   |   8   |
    +========+=======+======+=======+=======+=======+=======+=======+
    | BCTSET | CSID  | SID1 | TID1  | FRIC1 | MIND1 | MAXD1 |  DID  |
    +--------+-------+------+-------+-------+-------+-------+-------+
    |        |       | SID2 | TID2  | FRIC2 | MIND2 | MAXD2 |       |
    +--------+-------+------+-------+-------+-------+-------+-------+
    |        |  etc. |      |       |       |       |       |       |
    +--------+-------+------+-------+-------+-------+-------+-------+
    """
    def clear(self):
        #self.sid = np.array([], dtype='int32')
        #self.nelement = np.array([], dtype='int32')
        #self.element_id = np.array([], dtype='int32')
        #: GSID Contact set identification number. (Integer > 0)
        self.contact_id = np.array([], dtype='int32')

        # Identification number of a DESC bulk entry, which includes
        # a contact set description. Only supported by SOL 402.
        self.desc_id = np.array([], dtype='int32')

        #: SIDi Source region (contactor) identification number for contact pair i.
        #: (Integer > 0)
        self.source_ids = np.array([], dtype='int32')

        #: TIDi Target region identification number for contact pair i. (Integer > 0)
        self.target_ids = np.array([], dtype='int32')

        #: FRICi coefficient of friction for contact pair i. (Real); (Default=0.0)
        self.friction = np.array([], dtype='float64')

        # MINDi Minimum search distance for contact. (Real; Default=0.0)
        self.min_search_distance = np.array([], dtype='float64')

        #: SDISTi Search distance for glue regions (Real); (Default=0.0)
        self.max_search_distance = np.array([], dtype='float64')

    def add(self, contact_id: list[int],
            source_ids: list[int], target_ids: list[int],
            frictions: list[float],
            min_distances: list[float], max_distances: list[float],
            desc_id: int=0,
            comment: str='') -> int:
        """Adds a BCTSET card, which defines generalized contact"""
        self.cards.append((contact_id, source_ids, target_ids,
                           frictions,
                           min_distances, max_distances,
                           desc_id, comment))
        self.n += 1
        return self.n - 1

    #def remove_unused(self)
    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a BGSET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        contact_id = integer(card, 1, 'csid')
        desc_id = integer_or_blank(card, 7, 'desc_id')
        source_ids = []
        target_ids = []
        frictions = []
        min_distances = []
        max_distances = []

        nfields = card.nfields
        i = 2
        j = 1
        while i < nfields:
            source_ids.append(integer(card, i, 'sid%s' % j))
            target_ids.append(integer(card, i + 1, 'tid%s' % j))
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
        #return BCTSET(contact_id, source_ids, target_ids, frictions, min_distances,
                      #max_distances, comment=comment,
                      #sol=sol)

        assert len(card) > 2, f'len({self.type} card) = {len(card):d}\ncard={card}'
        self.cards.append((contact_id, source_ids, target_ids, frictions,
                           min_distances, max_distances,
                           desc_id, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        #ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        ncards = len(self.cards)
        idtype = self.model.idtype
        fdtype = self.model.fdtype
        contact_id = np.zeros(ncards, dtype='int32')
        desc_id = np.zeros(ncards, dtype='int32')
        nsource = np.zeros(ncards, dtype='int32')
        source_ids_list = []
        target_ids_list = []
        frictions_list = []
        min_search_distances_list = []
        max_search_distances_list = []

        #comment = {}
        for icard, card in enumerate(self.cards):
            (contact_idi, source_idsi, target_idsi, frictionsi,
             min_search_distancei, max_search_distancei,
             desc_idi, commenti) = card

            contact_id[icard] = contact_idi

            source_ids_list.extend(source_idsi)
            target_ids_list.extend(target_idsi)
            frictions_list.extend(frictionsi)
            min_search_distances_list.extend(min_search_distancei)
            max_search_distances_list.extend(max_search_distancei)

            nsource[icard] = len(source_idsi)
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti

        source_ids = np.array(source_ids_list, dtype=idtype)
        target_ids = np.array(target_ids_list, dtype=idtype)
        min_search_distances = np.array(min_search_distances_list, dtype=fdtype)
        max_search_distances = np.array(max_search_distances_list, dtype=fdtype)
        frictions = np.array(frictions_list, dtype=fdtype)
        self._save(contact_id, nsource, source_ids, target_ids,
                   frictions, min_search_distances, max_search_distances,
                   desc_id, comment=None)
        #self.sort()
        self.cards = []

    def _save(self,
              contact_id: np.ndarray,
              nsource: np.ndarray,
              source_ids: np.ndarray,
              target_ids: np.ndarray,
              frictions: np.ndarray,
              min_search_distances: np.ndarray,
              max_search_distances: np.ndarray,
              desc_id: np.ndarray,
              comment: dict[int, str]=None) -> None:
        ncards_existing = len(self.contact_id)
        if ncards_existing != 0:
            asdf
            contact_id = np.hstack([self.contact_id, contact_id])
            #nelement = np.hstack([self.nelement, nelement])
            #element_ids = np.hstack([self.element_ids, element_ids])
        #if comment:
            #self.comment.update(comment)
        self.contact_id = contact_id
        self.desc_id = desc_id

        self.nsource = nsource
        self.source_ids = source_ids
        self.target_ids = target_ids
        self.frictions = frictions
        self.min_search_distances = min_search_distances
        self.max_search_distances = max_search_distances
        self.n = len(self.contact_id)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        contact_regions = np.hstack([self.source_ids, self.target_ids])
        used_dict['contact_set_id'].append(contact_regions)

    def convert(self, xyz_scale: float=1.0, **kwargs) -> None:
        self.min_search_distances *= xyz_scale
        self.max_search_distances *= xyz_scale

    def __apply_slice__(self, bctset: BCTSET, i: np.ndarray) -> None:
        self._slice_comment(bctset, i)
        bctset.n = len(i)

        isource = self.isource
        bctset.contact_id = self.contact_id[i]
        bctset.desc_id = self.desc_id[i]
        bctset.nsource = self.nsource[i]
        bctset.source_ids = hslice_by_idim(i, isource, self.source_ids)
        bctset.target_ids = hslice_by_idim(i, isource, self.target_ids)
        bctset.min_search_distances = hslice_by_idim(i, isource, self.min_search_distances)
        bctset.max_search_distances = hslice_by_idim(i, isource, self.max_search_distances)
        bctset.frictions = hslice_by_idim(i, isource, self.frictions)

    @property
    def isource(self) -> np.ndarray:
        return make_idim(self.n, self.nsource)

    @property
    def max_id(self) -> int:
        return max(self.contact_id.max(), self.source_ids.max(), self.target_ids.max())

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if self.n == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)

        contact_ids = array_str(self.contact_id, size=size)
        desc_ids = array_default_int(self.desc_id, size=size)
        target_ids_ = array_str(self.target_ids, size=size).tolist()
        source_ids_ = array_str(self.source_ids, size=size).tolist()
        frictions_ = array_default_float(self.frictions, default=0.0, size=size).tolist()
        min_search_distances_ = array_default_float(self.min_search_distances, default=0.0, size=size).tolist()
        max_search_distances_ = array_default_float(self.max_search_distances, default=0.0, size=size).tolist()

        for contact_id, desc_id, (isource0, isource1) in zip(contact_ids, desc_ids, self.isource):
            target_ids = target_ids_[isource0:isource1]
            source_ids = source_ids_[isource0:isource1]
            frictions = frictions_[isource0:isource1]
            min_search_distances = min_search_distances_[isource0:isource1]
            max_search_distances = max_search_distances_[isource0:isource1]
            list_fields = ['BCTSET', contact_id]
            assert len(target_ids) > 0
            for (source_id, target_id, friction, dx_min,
                 dx_max) in zip(target_ids, source_ids, frictions,
                                min_search_distances, max_search_distances):
                list_fields += [source_id, target_id, friction, dx_min, dx_max, None, None, desc_id]
                desc_id = ''
            bdf_file.write(print_card(list_fields))
        return


class BGADD(ADD):
    _id_name = 'bgadd_id'
    @property
    def bgadd_id(self):
        return self.sid
    @bgadd_id.setter
    def bgadd_id(self, bgadd_id: np.ndarray):
        self.sid = bgadd_id

    @property
    def glue_ids(self):
        return self.sids
    @glue_ids.setter
    def glue_ids(self, glue_ids: np.ndarray):
        self.sids = glue_ids

    @property
    def nglue_ids(self):
        return self.nsids
    @nglue_ids.setter
    def nglue_ids(self):
        return self.nsids

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['glue_id'].append(self.glue_ids)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.bgadd_id) == 0:
            return
        if size == 8 and self.is_small_field:
            print_card = print_card_8
        else:
            print_card = print_card_16

        #self.get_reduced_spcs()
        bgadd_ids = array_str(self.bgadd_id, size=size)
        bgset_ids = array_str(self.glue_ids, size=size)
        for bgadd_id, idim in zip(bgadd_ids, self.idim):
            idim0, idim1 = idim
            bgset_idsi = bgset_ids[idim0:idim1].tolist()
            list_fields = ['BGADD', bgadd_id] + bgset_idsi
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        glue_id = np.unique(self.model.bgset.glue_id)
        #geom_check(self,
                   #missing,
                   #dconstr_id=(dconstr_id, self.dconstr_ids),
                   #)

class BCTADD(ADD):
    _id_name = 'bctadd_id'
    @property
    def bctadd_id(self):
        return self.sid
    @bctadd_id.setter
    def bctadd_id(self, bctadd_id: np.ndarray):
        self.sid = bctadd_id

    @property
    def contact_ids(self):
        return self.sids
    @contact_ids.setter
    def contact_ids(self, contact_ids: np.ndarray):
        self.sids = contact_ids

    @property
    def ncontact_ids(self):
        return self.nsids
    @ncontact_ids.setter
    def ncontact_ids(self):
        return self.nsids

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['contact_id'].append(self.contact_ids)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.bctadd_id) == 0:
            return
        if size == 8 and self.is_small_field:
            print_card = print_card_8
        else:
            print_card = print_card_16

        #self.get_reduced_spcs()
        bctadd_ids = array_str(self.bctadd_id, size=size)
        bctset_idsi = array_str(self.contact_ids, size=size)
        for bctadd_id, idim in zip(bctadd_ids, self.idim):
            idim0, idim1 = idim
            bctset_idsi = bctset_idsi[idim0:idim1].tolist()
            list_fields = ['BCTADD', bctadd_id] + bctset_idsi
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        contact_id = np.unique(self.model.bctset.contact_id)
        #geom_check(self,
                   #missing,
                   #dconstr_id=(dconstr_id, self.dconstr_ids),
                   #)


class BCONP(VectorizedBaseCard):
    """
    3D Contact Region Definition by Shell Elements (SOLs 101, 601 and 701)

    Defines a 3D contact region by shell element IDs.

    +-------+----+-------+--------+-----+------+--------+-------+-----+
    |   1   |  2 |   3   |   4    |  5  |   6  |   7    |   8   |  9  |
    +=======+====+=======+========+=====+======+========+=======+=====+
    | BCONP | ID | SLAVE | MASTER |     | SFAC | FRICID | PTYPE | CID |
    +-------+----+-------+--------+-----+------+--------+-------+-----+
    | BCONP | 95 |   10  |   15   |     |  1.0 |   33   |   1   |     |
    +-------+----+-------+--------+-----+------+--------+-------+-----+

    """
    _id_name = 'contact_id'
    def clear(self):
        self.contact_id = np.array([], dtype='int32')
        self.slave_id = np.array([], dtype='int32')
        self.master_id = np.array([], dtype='int32')
        self.sfac = np.array([], dtype='float64')
        self.friction_id = np.array([], dtype='int32')
        self.ptype = np.array([], dtype='int32')
        self.coord_id = np.array([], dtype='int32')

    def add(self, contact_id: int, slave: int, master: int,
            fric_id: int, sfac: float=1.0, ptype: int=1, cid: int=0,
            comment: str='') -> int:
        """Creates a BCONP card"""
        self.cards.append((contact_id, slave, master, sfac, fric_id, ptype, cid, comment))
        self.n += 1
        return self.n - 1

    #def remove_unused(self)
    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a BCONP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        fdouble_or_blank = lax_double_or_blank if self.model.is_lax_parser else double_or_blank

        contact_id = integer(card, 1, 'contact_id')
        slave = integer(card, 2, 'slave')
        master = integer(card, 3, 'master')
        sfac = fdouble_or_blank(card, 5, 'sfac', default=1.0)
        friction_id = integer_or_blank(card, 6, 'fric_id', default=0)
        ptype = integer_or_blank(card, 7, 'ptype', default=1)
        cid = integer_or_blank(card, 8, 'cid', default=0)
        #return BCONP(contact_id, slave, master, sfac, friction_id, ptype, cid,
                     #comment=comment)
        self.cards.append((contact_id, slave, master, sfac, friction_id, ptype, cid, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        #ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        ncards = len(self.cards)
        idtype = self.model.idtype
        fdtype = self.model.fdtype
        contact_id = np.zeros(ncards, dtype='int32')
        slave_id = np.zeros(ncards, dtype='int32')
        master_id = np.zeros(ncards, dtype='int32')
        sfac = np.zeros(ncards, dtype='float64')
        friction_id = np.zeros(ncards, dtype='int32')
        ptype = np.zeros(ncards, dtype='int32')
        coord_id = np.zeros(ncards, dtype='int32')

        #comment = {}
        for icard, card in enumerate(self.cards):
            (contact_idi, slave_idi, master_idi, sfaci, friction_idi, ptypei, cid, commenti) = card
            print(card)

            contact_id[icard] = contact_idi
            slave_id[icard] = slave_idi
            master_id[icard] = master_idi
            ptype[icard] = ptypei
            sfac[icard] = sfaci
            friction_id[icard] = friction_idi
            coord_id[icard] = cid
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti

        self._save(contact_id, slave_id, master_id, sfac, friction_id, ptype, coord_id)
        self.sort()
        self.cards = []

    def _save(self, contact_id, slave_id, master_id, sfac, friction_id, ptype, coord_id) -> None:
        ncards_existing = len(self.contact_id)
        if ncards_existing != 0:
            asdf
            contact_id = np.hstack([self.contact_id, contact_id])
        #if comment:
            #self.comment.update(comment)
        self.contact_id = contact_id
        self.slave_id = slave_id
        self.master_id = master_id
        self.sfac = sfac
        self.friction_id = friction_id
        self.ptype = ptype
        self.coord_id = coord_id
        self.n = len(self.contact_id)

    #def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        #asd
        #contact_regions = np.hstack([self.source_ids, self.target_ids])
        #used_dict['contact_set_id'].append(contact_regions)

    def convert(self, xyz_scale: float=1.0, **kwargs) -> None:
        pass
        #self.search_distance *= xyz_scale

    def __apply_slice__(self, bconp: BCONP, i: np.ndarray) -> None:
        #self._slice_comment(bconp, i)
        bconp.n = len(i)

        bconp.contact_id = self.contact_id[i]
        bconp.slave_id = self.slave_id[i]
        bconp.master_id = self.master_id[i]
        bconp.sfac = self.sfac[i]
        bconp.friction_id = self.friction_id[i]
        bconp.ptype = self.ptype[i]
        bconp.coord_id = self.coord_id[i]

    @property
    def max_id(self) -> int:
        return max(self.contact_id.max(), self.slave_id.max(),
                   self.master_id.max(), self.coord_id.max())

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if self.n == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)
        contact_ids = array_str(self.contact_id, size=size)
        slave_ids = array_str(self.slave_id, size=size)
        master_ids = array_str(self.master_id, size=size)
        sfacs = array_default_float(self.sfac, default=1., size=size, is_double=False)
        friction_ids = array_str(self.friction_id, size=size)
        ptypes = array_str(self.ptype, size=size)
        coord_ids = array_str(self.coord_id, size=size)

        for contact_id, slave_id, master_id, sfac, friction_id, ptype, cid in zip(
            contact_ids, slave_ids, master_ids, sfacs, friction_ids, ptypes, coord_ids):
            list_fields = [
                'BCONP', contact_id, slave_id, master_id, None, sfac,
                friction_id, ptype, cid]
            bdf_file.write(print_card(list_fields))
        return


class BFRIC(VectorizedBaseCard):
    """
    Slideline Contact Friction
    Defines frictional properties between two bodies in contact.

    +-------+------+-------+-------+-------+
    |   1   |   2  |   3   |   4   |   5   |
    +=======+======+=======+=======+=======+
    | BFRIC | FID  | FSTIF |       |  MU1  |
    +-------+------+-------+-------+-------+

    """
    _id_name = 'friction_id'
    def clear(self):
        self.friction_id = np.array([], dtype='int32')
        self.fstiff = np.array([], dtype='float64')
        self.mu1 = np.array([], dtype='float64')

    def add(self, friction_id: int, mu1: float, fstiff: float=None, comment: str='') -> int:
        """Creates a BFRIC card"""
        self.cards.append((friction_id, fstiff, mu1, comment))
        self.n += 1
        return self.n - 1

    #def remove_unused(self)
    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a BCONP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        fdouble_or_blank = lax_double_or_blank if self.model.is_lax_parser else double_or_blank

        friction_id = integer(card, 1, 'friction_id')
        fstiff = fdouble_or_blank(card, 2, 'fstiff', default=np.nan)  # default=program selected
        #
        mu1 = double(card, 4, 'mu1')
        #return BFRIC(friction_id, mu1, fstiff=fstiff, comment=comment)
        self.cards.append((friction_id, fstiff, mu1, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        #idtype = self.model.idtype
        #fdtype = self.model.fdtype
        friction_id = np.zeros(ncards, dtype='int32')
        fstiff = np.zeros(ncards, dtype='float64')
        mu1 = np.zeros(ncards, dtype='float64')

        #comment = {}
        for icard, card in enumerate(self.cards):
            (friction_idi, fstiffi, mu1i, commenti) = card

            friction_id[icard] = friction_idi
            fstiff[icard] = fstiffi
            mu1[icard] = mu1i
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti

        self._save(friction_id, fstiff, mu1)
        self.sort()
        self.cards = []

    def _save(self, friction_id, fstiff, mu1) -> None:
        ncards_existing = len(self.friction_id)
        if ncards_existing != 0:
            friction_id = np.hstack([self.friction_id, friction_id])
            fstiff = np.hstack([self.fstiff, fstiff])
            mu1 = np.hstack([self.mu1, mu1])
        #if comment:
            #self.comment.update(comment)
        self.friction_id = friction_id
        self.fstiff = fstiff
        self.mu1 = mu1
        self.n = len(self.friction_id)

    #def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        #asd
        #contact_regions = np.hstack([self.source_ids, self.target_ids])
        #used_dict['contact_set_id'].append(contact_regions)

    #def convert(self, xyz_scale: float=1.0, **kwargs) -> None:
        #pass
        #self.search_distance *= xyz_scale

    def __apply_slice__(self, bfric: BFRIC, i: np.ndarray) -> None:
        #self._slice_comment(bfric, i)
        bfric.n = len(i)
        bfric.friction_id = self.friction_id[i]
        bfric.fstiff = self.fstiff[i]
        bfric.mu1 = self.mu1[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        return self.friction_id.max()

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if self.n == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)
        friction_ids = array_str(self.friction_id, size=size)
        fstiffs = array_float(self.fstiff, size=size, is_double=False)
        mu1s = array_float(self.mu1, size=size, is_double=False)

        for friction_id, fstiff, mu1 in zip(friction_ids, fstiffs, mu1s):
            list_fields = [
                'BFRIC', friction_id, fstiff, None, mu1]
            bdf_file.write(print_card(list_fields))
        return


class BCRPARA(VectorizedBaseCard):
    """
    +---------+------+------+--------+------+-----+---+---+---+----+
    |    1    |   2  |   3  |   4    |   5  |  6  | 7 | 8 | 9 | 10 |
    +=========+======+======+========+======+=====+===+===+===+====+
    | BCRPARA | CRID | SURF | OFFSET | TYPE | GP  |   |   |   |    |
    +---------+------+------+--------+------+-----+---+---+---+----+
    """
    _id_name = 'contact_region_id'
    def clear(self):
        #: CRID Contact region ID. (Integer > 0)
        self.contact_region_id = np.array([], dtype='int32')

        #: SURF Indicates the contact side. See Remark 1. (Character = "TOP" or
        #: "BOT"; Default = "TOP")
        self.surf = np.array([], dtype='|U8')

        #: Offset distance for the contact region. See Remark 2. (Real > 0.0,
        #: Default =OFFSET value in BCTPARA entry)
        self.offset = np.array([], dtype='float64')

        #: Indicates whether a contact region is a rigid surface if it is used as a
        #: target region. See Remarks 3 and 4. (Character = "RIGID" or "FLEX",
        #: Default = "FLEX"). This is not supported for SOL 101.
        self.surface_type = np.array([], dtype='|U8')

        #: Control grid point for a target contact region with TYPE=RIGID or
        #: when the rigid-target algorithm is used. The grid point may be
        #: used to control the motion of a rigid surface. (Integer > 0)
        #: This is not supported for SOL 101.
        self.grid_point = np.array([], dtype='int32')

    def add(self, contact_region_id: int, offset: Optional[float]=None,
            surf: str='TOP', surface_type: str='FLEX', grid_point: int=0,
            comment: str='') -> int:
        """Creates a BFRIC card

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
        self.cards.append((contact_region_id, offset, surf, surface_type, grid_point, comment))
        self.n += 1
        return self.n - 1

    #def remove_unused(self)
    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a BCONP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        """
        Adds a BCRPARA card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        contact_region_id = integer(card, 1, 'crid')
        surf = string_or_blank(card, 2, 'surf', default='TOP')
        offset = double_or_blank(card, 3, 'offset', default=None)
        surface_type = string_or_blank(card, 4, 'type', default='FLEX')
        grid_point = integer_or_blank(card, 5, 'grid_point', default=0)
        #return BCRPARA(crid, surf=surf, offset=offset, Type=Type,
                       #grid_point=grid_point, comment=comment)
        self.cards.append((contact_region_id, offset, surf,
                           surface_type, grid_point, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        #idtype = self.model.idtype
        #fdtype = self.model.fdtype
        contact_region_id = np.zeros(ncards, dtype='int32')
        offset = np.zeros(ncards, dtype='float64')
        surf = np.zeros(ncards, dtype='|U8')
        surface_type = np.zeros(ncards, dtype='|U8')
        grid_point = np.zeros(ncards, dtype='int32')

        #comment = {}
        for icard, card in enumerate(self.cards):
            (contact_region_idi, offseti, surfi,
             surface_typei, grid_pointi, commenti) = card

            contact_region_id[icard] = contact_region_idi
            offset[icard] = offseti
            surf[icard] = surfi
            surface_type[icard] = surface_typei
            grid_point[icard] = grid_pointi
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti

        self._save(contact_region_id, offset, surf, surface_type, grid_point)
        self.sort()
        self.cards = []

    def _save(self, contact_region_id, offset, surf, surface_type, grid_point) -> None:
        ncards_existing = len(self.contact_region_id)
        if ncards_existing != 0:
            contact_region_id = np.hstack([self.contact_region_id, contact_region_id])
            offset = np.hstack([self.offset, offset])
            surf = np.hstack([self.surf, surf])
            surface_type = np.hstack([self.surface_type, surface_type])
            grid_point = np.hstack([self.grid_point, grid_point])
        #if comment:
            #self.comment.update(comment)
        self.contact_region_id = contact_region_id
        self.offset = offset
        self.surf = surf
        self.surface_type = surface_type
        self.grid_point = grid_point
        self.n = len(self.contact_region_id)

    #def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        #asd
        #contact_regions = np.hstack([self.source_ids, self.target_ids])
        #used_dict['contact_set_id'].append(contact_regions)

    #def convert(self, xyz_scale: float=1.0, **kwargs) -> None:
        #pass
        #self.search_distance *= xyz_scale

    def __apply_slice__(self, bcrpara: BCRPARA, i: np.ndarray) -> None:
        #self._slice_comment(bcrpara, i)
        bcrpara.n = len(i)
        bcrpara.contact_region_id = self.contact_region_id[i]
        bcrpara.surf = self.surf[i]
        bcrpara.offset = self.offset[i]
        bcrpara.surface_type = self.surface_type[i]
        bcrpara.grid_point = self.grid_point[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        nids = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nids, self.grid_point), filter_node0=True)

    @property
    def max_id(self) -> int:
        return max(self.contact_region_id.max(), self.grid_point.max())

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if self.n == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)
        contact_ids = array_str(self.contact_region_id, size=size)
        offsets = array_float(self.offset, size=size, is_double=False)
        grid_points = array_str(self.grid_point, size=size)

        for contact_id, surf, offset, surface_type, grid_point in zip(contact_ids, self.surf,
                                                                      offsets, self.surface_type, grid_points):
            list_fields = ['BCRPARA', contact_id, surf, offset, surface_type, grid_point]
            bdf_file.write(print_card(list_fields))
        return


class BEDGE(VectorizedBaseCard):
    """
    Defines a Region of Element Edges
    Defines a region of element edges for:
     - glue (BGSET entry)
     - contact (BCTSET entry)
     - fluid pressure load (PLOADFP entry)
     - wetted edges for a Co-simulation (CSMSET entry)

    +-------+------+-------+--------+--------+-------+------+--------+--------+
    |   1   |   2  |   3   |    4   |    5   |   6   |   7  |   8    |    9   |
    +=======+======+=======+========+========+=======+======+========+========+
    | BEDGE |  ID  |  EID1 | GID1,1 | GID1,2 |       | EID2 | GID2,1 | GID2,2 |
    +-------+------+-------+--------+--------+-------+------+--------+--------+
    |       |      |  EID3 | GID3,1 | GID3,2 | -etc- |      |        |        |
    +-------+------+-------+--------+--------+-------+------+--------+--------+

    """
    _id_name = 'bedge_id'
    def clear(self):
        self.bedge_id = np.array([], dtype='int32')
        self.elements = np.array([], dtype='int32')
        self.nelement = np.array([], dtype='int32')
        self.nodes = np.array([], dtype='int32')

    def add(self, bedge_id: int,
            eids: list[int],
            grids: list[tuple[int, int]],
            comment: str='') -> int:
        """Creates a BEDGE card"""
        self.cards.append((bedge_id, eids, grids, comment))
        self.n += 1
        return self.n - 1

    #def remove_unused(self)
    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a BCONP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #|   1   |   2  |   3   |    4   |    5   |   6   |   7  |   8    |    9   |
        #| BEDGE |  ID  |  EID1 | GID1,1 | GID1,2 |       | EID2 | GID2,1 | GID2,2 |
        #|       |      |  EID3 | GID3,1 | GID3,2 | -etc- |      |        |        |
        bedge_id = integer(card, 1, 'bedge_id')
        nfields = len(card) - 1

        nlines = nfields // 8
        if nfields % 8:
            nlines += 1

        ifield = 1
        eids = []
        nids = []
        neids = 1
        for iline in range(nlines):
            i0 = iline * 8 + ifield
            eid = integer_or_blank(card, i0, f'eid_{neids}', default=0)
            if eid:
                nid1 = integer(card, i0+1, f'nid1_{neids}')
                nid2 = integer(card, i0+2, f'nid2_{neids}')
                eids.append(eid)
                nids.append((nid1, nid2))
            neids += 1

            i0 += 4
            eid = integer_or_blank(card, i0, f'eid_{neids}', default=0)
            if eid:
                nid1 = integer(card, i0+1, f'nid1_{neids}')
                nid2 = integer(card, i0+2, f'nid2_{neids}')
                eids.append(eid)
                nids.append((nid1, nid2))
            neids += 1

        self.cards.append((bedge_id, eids, nids, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        #idtype = self.model.idtype
        #fdtype = self.model.fdtype
        bedge_id = np.zeros(ncards, dtype='int32')
        nelement = np.zeros(ncards, dtype='int32')

        #comment = {}
        element_list = []
        nodes_list = []
        for icard, card in enumerate(self.cards):
            (bedge_idi, eidsi, nidsi, commenti) = card
            bedge_id[icard] = bedge_idi
            nelement[icard] = len(eidsi)
            element_list.extend(eidsi)
            nodes_list.extend(nidsi)
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti

        element_id = np.array(element_list, dtype='int32')
        nodes = np.array(nodes_list, dtype='int32')
        assert nodes.ndim == 2, nodes.shape
        self._save(bedge_id, nelement, element_id, nodes)
        self.sort()
        self.cards = []

    def _save(self, bedge_id, nelement, element_id, nodes) -> None:
        ncards_existing = len(self.bedge_id)
        if ncards_existing != 0:
            bedge_id = np.hstack([self.bedge_id, bedge_id])
            nelement = np.hstack([self.nelement, nelement])
            element_id = np.hstack([self.element_id, element_id])
            nodes = np.hstack([self.nodes, nodes])
        #if comment:
            #self.comment.update(comment)
        self.bedge_id = bedge_id
        ielement = self.ielement
        self.nelement = nelement
        self.element_id = element_id
        self.nodes = nodes
        self.n = len(self.bedge_id)

    #def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        #asd
        #contact_regions = np.hstack([self.source_ids, self.target_ids])
        #used_dict['contact_set_id'].append(contact_regions)

    def convert(self, **kwargs) -> None:
        pass

    def __apply_slice__(self, bedge: BEDGE, i: np.ndarray) -> None:
        #self._slice_comment(bedge, i)
        bedge.bedge_id = self.bedge_id[i]
        ielement = self.ielement
        bedge.element_id = hslice_by_idim(i, ielement, self.element_id)
        bedge.nodes = vslice_by_idim(i, ielement, self.nodes)
        bedge.element_id = self.element_id[i]
        bedge.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        """
        Defines a region of element edges for:
         - glue (BGSET entry)
         - contact (BCTSET entry)
         - fluid pressure load (PLOADFP entry)
         - wetted edges for a Co-simulation (CSMSET entry)
        """
        nids = self.model.grid.node_id
        unids = np.unique(self.nodes.ravel())
        geom_check(self,
                   missing,
                   node=(nids, unids))

    @property
    def ielement(self) -> np.ndarray:
        return make_idim(self.n, self.nelement)

    @property
    def max_id(self) -> int:
        return max(self.bedge_id.max(), self.element_id.max(), self.nodes.max())

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if self.n == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)
        bedge_ids = array_str(self.bedge_id, size=size).tolist()
        element_id_ = array_str(self.element_id, size=size).tolist()
        nodes_ = array_str(self.nodes, size=size)

        for bedge_id, (ieid0, ieid1) in zip(bedge_ids, self.ielement):
            list_fields = ['BEDGE', bedge_id]
            element_id = element_id_[ieid0:ieid1]
            nodes = nodes_[ieid0:ieid1, :].tolist()
            for eid, (nid1, nid2) in zip(element_id, nodes):
                list_fields.extend([eid, nid1, nid2, None])
            bdf_file.write(print_card(list_fields))
        return
