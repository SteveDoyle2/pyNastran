from __future__ import annotations
from typing import Optional, TYPE_CHECKING
from collections import defaultdict
#from itertools import count
import numpy as np

#from pyNastran.utils.numpy_utils import integer_types, float_types
#from pyNastran.bdf import MAX_INT
from pyNastran.bdf.cards.base_card import expand_thru_by # expand_thru,
from pyNastran.bdf.cards.collpase_card import collapse_thru_packs # collapse_thru,
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, hslice_by_idim, vslice_by_idim, make_idim, parse_check)
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string, #blank,
    integer_or_blank, double_or_blank, string_or_blank,
    integer_double_or_blank, integer_string_or_blank,
    string_choice_or_blank,
    #components_or_blank, parse_components,
)
from pyNastran.bdf.bdf_interface.assign_type_force import lax_double_or_blank
#from pyNastran.bdf.field_writer_8 import print_card_8
#from pyNastran.bdf.field_writer_16 import print_card_16 # print_float_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.contact import _get_bcbody_section_values

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    remove_unused_primary, remove_unused_duplicate)
from .constraints import ADD
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    #array_default_str,
    array_str, array_float,
    array_default_int, array_default_float,
    array_float_nan, get_print_card_size)
#from pyNastran.dev.bdf_vectorized3.utils import print_card_8

from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF


class ElementPropertyNodeSet(VectorizedBaseCard):
    _id_name = 'sid'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
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
                fieldsi = expand_thru_by(fieldsi, set_fields=True, sort_fields=False)
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

    def __apply_slice__(self, eset: ElementPropertyNodeSet, i: np.ndarray) -> None:
        self._slice_comment(eset, i)
        eset.n = len(i)
        eset.sid = self.sid[i]

        i_id = self.i_id # [i, :]
        eset.ids = hslice_by_idim(i, i_id, self.ids)
        eset.n_ids = self.n_ids[i]

    @property
    def max_id(self) -> int:
        return max(self.sid.max(), self.ids.max())

    @property
    def i_id(self) -> np.ndarray:
        return make_idim(self.n, self.n_ids)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        contact_id = used_dict['contact_id']
        ncards_removed = remove_unused_duplicate(
            self, contact_id, self.sid, 'contact_id')
        return ncards_removed

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        set_to_ids = defaultdict(list)
        for sid, (i0, i1) in zip(self.sid, self.i_id):
            ids0 = self.ids[i0:i1].tolist()
            set_to_ids[sid].append(ids0)

        #bdf_file.write(f'$ {self.type}\n')
        for sid, ids_list in set_to_ids.items():
            ids = np.hstack(ids_list).tolist()
            singles, doubles = collapse_thru_packs(ids)
            #print(f'sid={sid} singles/doubles', len(singles), len(doubles))

            if len(singles):
                list_fields = [self.type, sid, ] + singles
                if len(doubles):
                    nfields = len(list_fields) - 1
                    nleftover = (nfields % 8)

                    if nleftover != 0:
                        list_fields.extend(['']*(8-nleftover))
                    for double in doubles:
                        list_fields += double + [''] * 5
                #print(print_card(list_fields).rstrip())
            elif len(doubles):
                assert len(singles) == 0, singles
                list_fields = [self.type, sid] + doubles[0] + ['', '', '', '']
                for double in doubles[1:]:
                    list_fields.extend(double + ['', '', '', '', ''])
                #print(print_card(list_fields).rstrip())
            else:  # pragma: no cover
                raise RuntimeError((singles, doubles))
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

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.node_id
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2



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
    #_id_name = 'bsurf_id'
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

class BOUTPUT(NodeSet):
    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.node_id),)


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
    _id_name = 'glue_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
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

        isource = self.isource
        bgset.glue_id = self.glue_id[i]
        bgset.nsource = self.nsource[i]
        bgset.source_ids = hslice_by_idim(i, isource, self.source_ids)
        bgset.target_ids = hslice_by_idim(i, isource, self.target_ids)
        bgset.search_distance = hslice_by_idim(i, isource, self.search_distance)
        bgset.extension = hslice_by_idim(i, isource, self.extension)
        bgset.n = len(i)

    @property
    def isource(self) -> np.ndarray:
        return make_idim(self.n, self.nsource)

    @property
    def max_id(self) -> int:
        return max(self.glue_id.max(), self.source_ids.max(), self.target_ids.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
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
    _id_name = 'contact_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
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

        isource = self.isource
        bctset.contact_id = self.contact_id[i]
        bctset.desc_id = self.desc_id[i]
        bctset.nsource = self.nsource[i]
        bctset.source_ids = hslice_by_idim(i, isource, self.source_ids)
        bctset.target_ids = hslice_by_idim(i, isource, self.target_ids)
        bctset.min_search_distances = hslice_by_idim(i, isource, self.min_search_distances)
        bctset.max_search_distances = hslice_by_idim(i, isource, self.max_search_distances)
        bctset.frictions = hslice_by_idim(i, isource, self.frictions)
        bctset.n = len(i)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        contact_id = used_dict['contact_id']
        ncards_removed = remove_unused_primary(
            self, contact_id, self.contact_id, 'contact_id')
        return ncards_removed

    @property
    def isource(self) -> np.ndarray:
        return make_idim(self.n, self.nsource)

    @property
    def max_id(self) -> int:
        return max(self.contact_id.max(), self.source_ids.max(), self.target_ids.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
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

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        contact_id = used_dict['contact_id']
        ncards_removed = remove_unused_primary(
            self, contact_id, self.bctadd_id, 'contact_id')
        return ncards_removed

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        #self.get_reduced_spcs()
        bctadd_ids = array_str(self.bctadd_id, size=size)
        bctset_ids = array_str(self.contact_ids, size=size).tolist()
        for bctadd_id, idim in zip(bctadd_ids, self.idim):
            idim0, idim1 = idim
            bctset_idsi = bctset_ids[idim0:idim1]
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
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
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

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['contact_id'].append(self.master_id)
        used_dict['contact_id'].append(self.slave_id)
        used_dict['bfric_id'].append(self.friction_id)
        used_dict['coord_id'].append(self.coord_id)

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
        #idtype = self.model.idtype
        #fdtype = self.model.fdtype
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
        bconp.contact_id = self.contact_id[i]
        bconp.slave_id = self.slave_id[i]
        bconp.master_id = self.master_id[i]
        bconp.sfac = self.sfac[i]
        bconp.friction_id = self.friction_id[i]
        bconp.ptype = self.ptype[i]
        bconp.coord_id = self.coord_id[i]
        bconp.n = len(i)

    @property
    def max_id(self) -> int:
        return max(self.contact_id.max(), self.slave_id.max(),
                   self.master_id.max(), self.coord_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
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
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
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
        bfric.friction_id = self.friction_id[i]
        bfric.fstiff = self.fstiff[i]
        bfric.mu1 = self.mu1[i]
        bfric.n = len(i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    @property
    def max_id(self) -> int:
        return self.friction_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        friction_ids = array_str(self.friction_id, size=size)
        fstiffs = array_float_nan(self.fstiff, size=size, is_double=False)
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
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
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
        contact_region_id : int
            CRID Contact region ID.
        offset : float; default=None
            Offset distance for the contact region (Real > 0.0).
            None : OFFSET value in BCTPARA entry
        surf : str; default='TOP'
            SURF Indicates the contact side. See Remark 1.  {'TOP', 'BOT'; )
        surface_type : str; default='FLEX'
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
        offset = double_or_blank(card, 3, 'offset', default=np.nan)
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
        bcrpara.contact_region_id = self.contact_region_id[i]
        bcrpara.surf = self.surf[i]
        bcrpara.offset = self.offset[i]
        bcrpara.surface_type = self.surface_type[i]
        bcrpara.grid_point = self.grid_point[i]
        bcrpara.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.grid_point)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        contact_id = used_dict['contact_id']
        ncards_removed = remove_unused_primary(
            self, contact_id, self.contact_region_id, 'contact_id')
        return ncards_removed

    def geom_check(self, missing: dict[str, np.ndarray]):
        nids = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nids, self.grid_point), filter_node0=True)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.grid_point
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @property
    def max_id(self) -> int:
        return max(self.contact_region_id.max(), self.grid_point.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        contact_ids = array_str(self.contact_region_id, size=size)
        offsets = array_float_nan(self.offset, size=size, is_double=False)
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
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
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

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['element_id'].append(self.element_id)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.nodes.ravel()
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

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

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
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


def bcbody_lines_to_card_rigid(lines: list[str]) -> tuple[BDFCard, BDFCard]:
    card_fields = []
    rigid_fields = []
    is_active_rigid = False
    for iline, line in enumerate(lines):
        #print(f'iline={iline} line={line!r}')
        nfields_expected = 9 if iline == 0 else 8
        line = line.expandtabs()

        #print(line)
        if '*' in line:
            #print(f'large {line!r}')
            # large field
            fields = _split_large_fields(line, iline)

            if fields[0].strip().upper().startswith('RIGID'):
                assert is_active_rigid is False, is_active_rigid
                assert is_active_rigid is not None, is_active_rigid
                rigid_fields.extend(fields)
                is_active_rigid = True
            elif is_active_rigid:
                assert is_active_rigid is not None, is_active_rigid
                rigid_fields.extend(fields)
                is_active_rigid = None
            else:
                card_fields.extend(fields)
            continue
        else:
            # small field
            if ',' in line:
                #print(f'comma {line!r}')
                fields = line.split(',')[:9]
                if iline > 0:
                    fields = fields[1:]

                    if len(fields) < nfields_expected:
                        nmissing = nfields_expected - len(fields)
                        fields.extend([''] * nmissing)

                    #print(f'  sline={fields}')
                    if fields[0].upper() == 'RIGID':
                        rigid_fields = fields
                        is_active_rigid = None
                    continue

                if len(fields) < nfields_expected:
                    nmissing = nfields_expected - len(fields)
                    fields.extend([''] * nmissing)
            else:
                #print(f'small {line!r}')
                if line[1:].strip().upper().startswith('RIGID'):
                    assert is_active_rigid is False, is_active_rigid
                    rigid_fields = _split_small_rigid_fields(line)
                    is_active_rigid = None
                    continue
                fields = _split_small_fields(line, iline)
        #print(f'iline={iline} fields={fields}; n={len(fields)}')
        assert len(fields) == nfields_expected, (fields, len(fields))
        card_fields.extend(fields)

    rigid_card = []
    if rigid_fields:
        rigid_fields2 = [
            '', rigid_fields[0], rigid_fields[1], rigid_fields[2],
            rigid_fields[3] + rigid_fields[4] + rigid_fields[5],
        ]
        rigid_card = BDFCard(rigid_fields2, has_none=False)

    card = BDFCard(card_fields, has_none=False)
    return card, rigid_card

def _split_large_fields(line: str, iline: int) -> list[str]:
    assert ',' not in line, line
    if iline == 0:
        fields = [
            line[0:8],
            line[8:24], line[24:40], line[40:56], line[56:72],
        ]
    else:
        fields = [
            line[8:24], line[24:40], line[40:56], line[56:72],
        ]
    return fields

def _split_small_fields(line: str, iline: int) -> list[str]:
    assert ',' not in line, line
    if iline == 0:
        fields = [
            line[0:8],
            line[8:16], line[16:24], line[24:32], line[32:40], line[40:48],
            line[48:56], line[56:64], line[64:72],
        ]
    else:
        fields = [
            line[8:16], line[16:24], line[24:32], line[32:40], line[40:48],
            line[48:56], line[56:64], line[64:72],
        ]
    return fields

def _split_small_rigid_fields(line: str) -> list[str]:
    if ',' in line:
        raise NotImplementedError(line)
    else:
        rigid_fields = [
            #line[0:8],
            line[8:16], line[16:24], line[24:32],
            line[32:40], line[40:48], line[48:56],
            #line[32:56],
        ]
    return rigid_fields


class BCBODY(VectorizedBaseCard):
    """
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |    1   |    2    |    3   |   4    |    5     |     6   |    7    |    8    |    9    |
    +========+=========+========+========+==========+=========+=========+=========+=========+
    | BCBODY | BID     | DIM    | BEHAV  |   BSID   |   ISTYP | FRIC    |  IDSPL  | CONTROL |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | NLOAD   | ANGVEL | DCOS1  |   DCOS2  |   DCOS3 | VELRB1  |  VELRB2 | VELRB3  |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | ADVANCE | SANGLE | COPTB  |   USER   |         |         |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | CTYPE   | ISMALL | ITYPE  |   IAUG   |  PENALT | AUGDIST |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | RIGID   | CGID   | NENT   | --- Rigid Body Name ---      |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | APPROV  | A      |  N1    | N2       |    N3   |    V1   |    V2   |    V3   |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | RTEMP   | G(temp)|  Tempr | T(Tempr) |         |         |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | SINK    | G(sink)|  Tsink | T(Tsink) |         |         |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | GROW    | GF1    |  GF2   |   GF3    | TAB-GF1 | TAB-GF2 | TAB-GF3 |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | HEAT    | CFILM  |  TSINK |   CHEAT  | TBODY   | HCV     | HNC     |  ITYPE  |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    |        | BNC     | EMISS  |  HBL   |          |         |         |         |         |
    +--------+---------+--------+--------+----------+---------+---------+---------+---------+
    """
    _id_name = 'bcbody_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.bcbody_id = np.array([], dtype='int32')

    #def add(self, bcbody: int,
            #eids: list[int],
            #grids: list[tuple[int, int]],
            #comment: str='') -> int:
        #"""Creates a BEDGE card"""
        #self.cards.append((bedge_id, eids, grids, comment))
        #self.n += 1
        #return self.n - 1

    #def remove_unused(self)
    def add_card(self, lines: list[str], comment: str='') -> int:
        """
        Adds a BCBODY card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        BID : int (4,1)
           Contact body identification number referenced by
           BCTABLE, BCHANGE, or BCMOVE. (Integer > 0; Required)
        DIM : str; default='3D'
           Dimension of body.
           DIM=2D planar body in x-y plane of the basic coordinate system,
                  composed of 2D elements or curves.
           DIM=3D any 3D body composed of rigid surfaces, shell elements or solid
                  elements.
        BEHAV (4,8)
           Behavior of curve or surface (Character; Default = DEFORM) DEFORM body is
           deformable, RIGID body is rigid, SYMM body is a symmetry body, ACOUS
           indicates an acoustic body, WORK indicates body is a workpiece, HEAT indicates
           body is a heat-rigid body. See Remark 3. for Rigid Bodies..
        BSID : int
            Identification number of a BSURF, BCBOX, BCPROP or BCMATL entry if
            BEHAV=DEFORM. (Integer > 0)
        ISTYP : int (4,3)
           Check of contact conditions. (Integer > 0; Default = 0)
        ISTYP : int
           is not supported in segment-to-segment contact.
           For a deformable body:
           =0 symmetric penetration, double sided contact.
           =1 unsymmetric penetration, single sided contact. (Integer > 0)
           =2 double-sided contact with automatic optimization of contact constraint
              equations (this option is known as “optimized contact”).
              Notes: single-sided contact (ISTYP=1) with the contact bodies arranged properly
              using the contact table frequently runs much faster than ISTYP=2.
              For a rigid body:
           =0 no symmetry condition on rigid body.
           =1 rigid body is a symmetry plane.
        FRIC : int/float (6,7)
            Friction coefficient. (Real > 0 or integer; Default = 0)
            If the value is an integer it represents the ID of a TABL3Di.
        IDSPL : int (4,5)
            Set IDSPL=1 to activate the SPLINE (analytical contact) option for a deformable
            body and for a rigid contact surface. Set it to zero or leave blank to not have
            analytical contact. (Integer; Default = 0)
        NLOAD : int or None
            Enter a positive number if "load controlled" and rotations are allowed (Integer). The
            positive number is the grid number where the moments or rotations are applied. The
            rotations are specified using SPCD at grid ID NLOAD and can be specified using dof's
            1-3 (for rotation about x, y, z respectively), or by dof's 4-6 (for rotation about x, y, z
            respectively).
            Note: This rotation takes the position of the grid point defined in CGID field as the
            center of rotation.
        ANGVEL : int/float; default=0.0
            Angular velocity or angular position about local axis through center of rotation. If the
            value is an integer it represents the ID of a TABLED1, TABLED2 or TABL3D, i.e., a
            time-dependent or multi-dimensional table; however, no log scales, only linear scales.
            (Real or Integer; Default = 0.0)
        DCOSi : int/float; default=0.0
            Components of direction cosine of local axis if ANGVEL is nonzero. If the value is an
            integer, it represents the ID of a TABLED1, TABLED2 or TABL3D, i.e., a time-dependent
            or multi-dimensional table; however, no log scales, only linear scales. (Real
            or Integer; Default=0.0) In 2D contact only DCOS3 is used and the Default is 1.0.
        VELRBi : int/float; default=0.0
            Translation velocity or final position (depending on the value of CONTROL) of rigid
            body at the grid point defined in CGID filed. For velocity control only, if the value is
            an integer, it represents the ID of TABLED1, TABLED2 or TABL3D, i.e., a time-dependent
            or multi-dimensional table; however, no log scales, only linear scales. Only
            VELRB1 and VELRB2 are used in 2D contact. (Real or Integer; Default = 0.0)

        """
        card, rigid_card = bcbody_lines_to_card_rigid(lines)

        contact_id = integer(card, 1, 'contact_id')
        field2 = card.field(2)  # string/None

        default = '3D'
        if field2 is None:
            dim = default
        else:
            field2 = field2.strip()
            if field2 in {'2', '3'}:
                dim = field2 + 'D'
            else:
                dim = string_choice_or_blank(card, 2, 'dim',
                                             ('2D', '3D'),
                                             default=default)

        behav = string_choice_or_blank(card, 3, 'behav',
                                       ('RIGID', 'DEFORM', 'SYMM', 'ACOUS', 'WORK', 'HEAT'),
                                       default='DEFORM')
        if behav == 'DEFORM':
            bsid = integer(card, 4, 'bsid')
        else:
            bsid = integer_double_or_blank(card, 4, 'bsid', default=0)

        istype = integer_or_blank(card, 5, 'istype', default=0)
        fric = integer_double_or_blank(card, 6, 'fric', default=0)
        idispl = integer_or_blank(card, 7, 'idispl', default=0)
        control = integer_or_blank(card, 8, 'control', default=0)

        # NLOAD   | ANGVEL | DCOS1  | DCOS2|  DCOS3 | VELRB1  | VELRB2 | VELRB3
        word_nload = integer_string_or_blank(card, 9, 'nload (int) / word (str)', default=None)
        i = 9
        if word_nload is None or isinstance(word_nload, int):
            nload = integer_or_blank(card, 9, 'nload', default=0) # made up
            ang_vel = double_or_blank(card, 10, 'ang_vel', default=0.0)
            dcos = [
                integer_double_or_blank(card, 11, 'dcos1', default=0.0),
                integer_double_or_blank(card, 12, 'dcos2', default=0.0),
                integer_double_or_blank(card, 13, 'dcos3', default=0.0),
            ]
            vel_rb = [
                integer_double_or_blank(card, 14, 'vel_rb1', default=0.0),
                integer_double_or_blank(card, 15, 'vel_rb2', default=0.0),
                integer_double_or_blank(card, 16, 'vel_rb3', default=0.0),
            ]
            i += 8
        else:
            nload = 0
            ang_vel = 0.0
            dcos = [0., 0., 0.]
            vel_rb = [0., 0., 0.]

        # advance
        sangle = 60.
        coptb = 0
        user = 0
        min_node = 0

        old_word = None
        while i < len(card):
            word = string_or_blank(card, i, 'word (str)', default=None)
            if word is None:
                raise RuntimeError(f'should be broken by {old_word}')

            #print('*', word)
            if word == 'ADVANCE':
                # TODO: USER is NX and MIN_NODE is MSC???
                #
                # | ADVANCE | SANGLE | COPTB | USER | MIDNO |
                sangle = double_or_blank(card, i+1, 'sangle', default=60.)
                coptb = integer_or_blank(card, i+2, 'coptb', default=0)
                user = integer_or_blank(card, i+3, 'user', default=0)
                min_node = integer_or_blank(card, i+4, 'min_node', default=0)
                # 'ADVANCE'
                #     The entries for this continuation line are for advanced options starting with
                #     MD Nastran R2.
                # SANGLE
                #     Threshold for automatic discontinuity detection in degrees. (Real; Default = 60.0)
                #     Used for SPLINE option in SOL 400 only. SANGLE is not used when IDSPL ≥ 0.
                # COPTB
                #     Flag to indicate how body surfaces may contact. See Remark 9. on the BCTABLE entry.
                #     (Integer; Default = 0)
                # MIDNOD
                #     Mid-side node projection flag. (Integer > 0; Default = 0)
                #     When MIDNOD > 0 and IDSPL 0, the mid-side grid of quadratic elements are
                #     projected onto the selected spline surfaces. This operation is performed before the
                #     contact process starts and it may change the location of grids in contact bodies. It may
                #     operate in combination with the initial stress-free contact.
                i += 8
            elif word == 'HEAT':
                # “HEAT”
                #    The entries of this continuation line(s) are for contact in heat transfer in a pure thermal
                #    analysis or in a coupled thermal/structural analysis. In a pure structural analysis they are
                #    ignored.
                # CFILM (9,1)/(10,1)
                #     Heat transfer coefficient (film) to environment. (Real or Integer, Default = 0.0) If Real,
                #     the value entered is the film coefficient. If Integer, the value entered is the ID of a
                #     TABLEM1 or TABLEM2 entry specifying the heat transfer coefficient vs temperature
                #     or a TABL3D entry specifying the film coefficient vs temperature and possibly other
                #     variables.
                # TSINK (9,2)/(10,2)
                #     Environment sink temperature. (Real or Integer, Default = 0.0). If Real, the value
                #     entered is the sink temperature. If Integer, the value entered is the ID of a TABLED1
                #     or TABLED2 entry specifying temperature vs time or a TABL3D entry specifying the
                #     sink temperature vs time and possibly other variables. When entered as a negative
                #     integer its absolute value is a scalar point identification number. If a scalar point is
                #     specified on this entry it need not be defined on an SPOINT entry.
                # CHEAT (9,3)/(10,3)
                #     Contact heat transfer coefficient. (Real or Integer; Default = 0.0). If Real, the value
                #     entered is the contact heat transfer coefficient. If Integer, the value entered is the ID of
                #     a TABLEM1 or TABLEM2 entry specifying the contact heat transfer coefficient vs
                #     temperature or a TABL3D entry specifying the contact heat transfer coefficient vs
                #     temperature and possibly other variables.
                # TBODY (9,4)/(10,4)
                #     Body temperature. (Real or Integer; Default = 0.0). If Real, the value entered is the body
                #     temperature. If Integer, the value entered is the ID of a TABLED1 or TABLED2 entry
                #     specifying the body temperature vs time or a TABL3D entry specifying the body
                #     temperature vs time and possibly other variables. When entered as a negative integer its
                #     absolute value is a scalar point identification number. If a scalar point is specified on
                #     this entry it need not be defined on an SPOINT entry.
                # HCV (9,5)/(10,5)
                #     Convection coefficient for near field behavior (Real or Integer; Default = 0.0). If Real
                #     the value entered is the near field convection coefficient. If Integer, the value entered is
                #     the ID of a TABLEM1 or TABLEM2 entry specifying the near field convection
                #     coefficient vs temperature or a TABL3D entry specifying the
                # HEAT CFILM TSINK CHEAT TBODY HCV HNC  ITYPE
                #      BNC   EMISS HBL   HNL   BNL HNLE BNLE
                #      HNCE  BNCE  CMB   CMS
                cfilm = integer_double_or_blank(card, i+1, 'cfilm', default=0.0)
                tsink = integer_double_or_blank(card, i+2, 'tsink', default=0.0)
                cheat = integer_double_or_blank(card, i+3, 'cheat', default=0.0)
                tbody = integer_double_or_blank(card, i+4, 'tbody', default=0.0)
                hcv = integer_double_or_blank(card, i+5, 'hcv', default=0.0)
                hnc = integer_double_or_blank(card, i+6, 'hnc', default=0.0)

                # no default...but this fails otherwise...
                #C:\MSC.Software\msc_nastran_runs\cpl_002_m1.dat
                itype = integer_or_blank(card, i+7, 'itype', default=0)
                i += 8

                bnc = double_or_blank(card, i+1, 'bnc', default=1.)
                emiss = double_or_blank(card, i+2, 'emiss', default=0.)
                hbl = double_or_blank(card, i+3, 'hbl', default=0.)
                hnl = integer_double_or_blank(card, i+4, 'hnl', default=0.)
                bnl = integer_double_or_blank(card, i+5, 'bnl', default=1.)
                hnle = integer_double_or_blank(card, i+6, 'hnle', default=0.)
                bnle = integer_double_or_blank(card, i+7, 'bnle', default=1.)
                i += 8

                hnce = integer_double_or_blank(card, i+1, 'hnce', default=0.)
                bnce = integer_double_or_blank(card, i+2, 'bnce', default=1.)
                cmb = double_or_blank(card, i+3, 'cmb', default=0.)
                cms = double_or_blank(card, i+4, 'cms', default=0.)
                i += 8
            elif word == 'GROW':
                #'GROW' GF1 GF2 GF3 TAB-GF1 TAB-GF2 TAB-GF3
                gf1 = double_or_blank(card, i+1, 'GF1', default=1.0)
                gf2 = double_or_blank(card, i+2, 'GF2', default=1.0)
                gf3 = double_or_blank(card, i+3, 'GF3', default=1.0)
                tab_gf1 = integer_or_blank(card, i+4, 'tab_GF1', default=0)
                tab_gf2 = integer_or_blank(card, i+5, 'tab_GF2', default=0)
                tab_gf3 = integer_or_blank(card, i+6, 'tab_GF3', default=0)
                #blank = blank(card, i+7, 'GROW blank')
                #print('grow values =', [gf1, gf2, gf3, tab_gf1, tab_gf2, tab_gf3])
                i += 8
                #GF1 GF2 GF3 TAB-GF1 TAB-GF2 TAB-GF3
            elif word == 'NURBS':
                i, values = _get_bcbody_section_values(card, i, word)
                #print('end of NURBS -> ', valuei)
            elif word in ('NURBS2', 'NURBS2D'):
                i, values = _get_bcbody_section_values(card, i, word)
                #print('end of NURBS -> ', valuei)
            elif word == 'PATCH3D':
                #'PATCH3D' NPATCH
                #          IDP G1 G2 G3 G4
                #          IDP G1 G2 G3 G4
                i, values = _get_bcbody_section_values(card, i, word)
            #elif word == 'RIGID':
                #i += 8
            elif word == 'BEZIER':
                i, values = _get_bcbody_section_values(card, i, word)
                #print(word, values)
            elif word == 'APPROV':
                i, values = _get_bcbody_section_values(card, i, word)
            else:  # pragma: no cover
                raise NotImplementedError(word)
            old_word = word

        if len(rigid_card):
            #CGID
            #(5,i) i=1,2,3
            #Grid point identification number defining the initial position of the reference point of
            #the rigid body or the point where a concentrated force or moment is applied
            #
            #NENT Number of geometric entities to describe this rigid surface. A rigid surface can be
            #described by multiple sets of patches, nurbs, etc. For example, if it takes 3 sets of
            #PATCH3D entries to describe a rigid surface, then set NENT=3.
            # (Integer > 0; Default=1)

            word = string(rigid_card, 1, 'RIGID')
            assert word == 'RIGID', word
            cgid = integer_or_blank(rigid_card, 2, 'cgid', default=-1)  #  TODO: made up default
            nent = integer_or_blank(rigid_card, 3, 'nent', default=1)
            rigid_body_name = string_or_blank(rigid_card, 4, 'rigid_body_name', default='dummy_name')
            #CGID NENT --- Rigid Body Name ---
        else:
            cgid = -1
            nent = 1
            rigid_body_name = ''
        #return BCBODY(contact_id, bsid,
                      #dim=dim, behav=behav, istype=istype,
                      #fric=fric, idispl=idispl,
                      #comment=comment)
        self.cards.append((contact_id, dim, behav, bsid, istype,
                           fric, idispl, control,
                           ('NLOAD', nload, dcos, ang_vel, vel_rb),
                           ('ADVANCE', sangle, coptb, user, min_node),
                           ('RIGID', cgid, nent, rigid_body_name),
                           comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        #idtype = self.model.idtype
        #fdtype = self.model.fdtype
        bcbody_id = np.zeros(ncards, dtype='int32')
        bsid = np.zeros(ncards, dtype='int32')
        dim = np.zeros(ncards, dtype='|U8')
        behavior = np.zeros(ncards, dtype='|U8')
        fric_int = np.zeros(ncards, dtype='float64')
        fric_float = np.zeros(ncards, dtype='float64')
        istype = np.zeros(ncards, dtype='int32')
        idispl = np.zeros(ncards, dtype='int32')
        control = np.zeros(ncards, dtype='int32')

        # line2
        nload = np.zeros(ncards, dtype='int32')
        ang_vel = np.zeros(ncards, dtype='float64')
        dcos = np.zeros((ncards, 3), dtype='float64')
        vel_rb = np.zeros((ncards, 3), dtype='float64')

        # ADVANCE
        sangle = np.zeros(ncards, dtype='float64')
        coptb = np.zeros(ncards, dtype='int32')
        user = np.zeros(ncards, dtype='int32')
        min_node = np.zeros(ncards, dtype='int32')

        # RIGID
        cgid = np.zeros(ncards, dtype='int32')
        nent = np.zeros(ncards, dtype='int32')
        rigid_body_name = np.zeros(ncards, dtype='|U24')

        #comment = {}
        #element_list = []
        #nodes_list = []
        for icard, card in enumerate(self.cards):
            (bcbody_idi,
             dimi, behavi, bsidi, istypei,
             frici, idispli, controli,
             nloadi, advancei, rigidi,
             commenti) = card
            assert isinstance(bsidi, int), bsidi
            bcbody_id[icard] = bcbody_idi
            bsid[icard] = bsidi
            dim[icard] = dimi
            behavior[icard] = behavi
            if isinstance(frici, float):
                fric_int[icard] = frici
            else:
                fric_float[icard] = frici
            istype[icard] = istypei
            idispl[icard] = idispli
            control[icard] = controli

            # nload
            (wprdo, nloadi, dcosi, ang_veli, vel_rbi) = nloadi
            nload[icard] = nloadi
            ang_vel[icard] = ang_veli
            dcos[icard, :] = dcosi
            vel_rb[icard, :] = vel_rbi

            # advance
            #('ADVANCE', sangle, coptb, user, min_node),
            (wordi, sanglei, coptbi, useri, min_nodei) = advancei
            sangle[icard] = sanglei
            coptb[icard] = coptbi
            user[icard] = useri
            min_node[icard] = min_nodei

            # rigid
            (wordi, cgidi, nenti, rigid_body_namei) = rigidi
            assert len(rigid_body_namei) < 24, rigid_body_namei
            cgid[icard] = cgidi
            nent[icard] = nenti
            rigid_body_name[icard] = rigid_body_namei
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti

        self._save(bcbody_id, dim, behavior, bsid, fric_int, fric_float, istype, idispl, control,
                   nload, ang_vel, dcos, vel_rb,
                   #advance
                   sangle, coptb, user, min_node,
                   #rigid
                   cgid, nent, rigid_body_name)
        self.sort()
        self.cards = []

    def _save(self, bcbody_id, dim, behavior, bsid, fric_int, fric_float, istype, idispl, control,
              nload, ang_vel, dcos, vel_rb,
              #advance
              sangle, coptb, user, min_node,
              #rigid
              cgid, nent, rigid_body_name) -> None:
        ncards_existing = len(self.bcbody_id)
        if ncards_existing != 0:
            bcbody_id = np.hstack([self.bcbody_id, bcbody_id])
            #nelement = np.hstack([self.nelement, nelement])
            #element_id = np.hstack([self.element_id, element_id])
            #nodes = np.hstack([self.nodes, nodes])
            asdf
        #if comment:
            #self.comment.update(comment)
        self.bcbody_id = bcbody_id
        self.dim = dim
        self.behavior = behavior
        self.bsid = bsid
        self.fric_int = fric_int
        self.fric_float = fric_float
        self.istype = istype
        self.idispl = idispl
        self.control = control

        self.nload = nload
        self.ang_vel = ang_vel
        self.dcos = dcos
        self.vel_rb = vel_rb

        # advance
        self.sangle = sangle
        self.coptb = coptb
        self.user = user
        self.min_node = min_node

        # rigid
        self.cgid = cgid
        self.nent = nent
        self.rigid_body_name = rigid_body_name
        self.n = len(self.bcbody_id)

    def __apply_slice__(self, bcbody: BCBODY, i: np.ndarray) -> None:
        #self._slice_comment(bedge, i)
        bcbody.bcbody_id = self.bcbody_id[i]
        bcbody.dim = self.dim[i]
        bcbody.behavior = self.behavior[i]
        bcbody.bsid = self.bsid[i]
        bcbody.fric_int = self.fric_int[i]
        bcbody.fric_float = self.fric_float[i]
        bcbody.istype = self.istype[i]
        bcbody.idispl = self.idispl[i]
        bcbody.control = self.control[i]

        bcbody.nload = self.nload[i]
        bcbody.ang_vel = self.ang_vel[i]
        bcbody.dcos = self.dcos[i, :]
        bcbody.vel_rb = self.vel_rb[i, :]

        # advance
        bcbody.sangle = self.sangle[i]
        bcbody.coptb = self.coptb[i]
        bcbody.user = self.user[i]
        bcbody.min_node = self.min_node[i]

        # rigid
        bcbody.cgid = self.cgid[i]
        bcbody.nent = self.nent[i]
        bcbody.rigid_body_name = self.rigid_body_name[i]
        bcbody.n = len(i)

    #def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        #used_dict['element_id'].append(self.element_id)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.nodes.ravel()
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    #def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        #asd
        #contact_regions = np.hstack([self.source_ids, self.target_ids])
        #used_dict['contact_set_id'].append(contact_regions)

    #def convert(self, **kwargs) -> None:
        #pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        """
        Defines a region of element edges for:
         - glue (BGSET entry)
         - contact (BCTSET entry)
         - fluid pressure load (PLOADFP entry)
         - wetted edges for a Co-simulation (CSMSET entry)
        """
        nids = self.model.grid.node_id
        #unids = np.unique(self.nodes.ravel())
        #geom_check(self,
                   #missing,
                   #node=(nids, unids))

    #@property
    #def ielement(self) -> np.ndarray:
        #return make_idim(self.n, self.nelement)

    @property
    def max_id(self) -> int:
        return self.bcbody_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        bcbody_ids = array_str(self.bcbody_id, size=size).tolist()

        istypes = array_default_int(self.istype, default=0)
        frics = array_default_float(self.fric_float, default=0.0)
        idispls = array_default_int(self.idispl, default=0)
        controls = array_default_int(self.control, default=0)

        nloads = array_str(self.nload, size=size).tolist()
        ang_vels = array_default_float(self.ang_vel, default=0., size=size, is_double=False)
        dcoss = array_default_float(self.dcos, default=0., size=size, is_double=False).tolist()
        vel_rbs = array_default_float(self.vel_rb, default=0., size=size, is_double=False).tolist()
        is_nloads = ~(
            (self.nload == 0) & (self.ang_vel == 0.0) &
            (self.dcos[:, 0]   == 0.0) & (self.dcos[:, 1]   == 0.0) & (self.dcos[:, 2]   == 0.0) &
            (self.vel_rb[:, 0] == 0.0) & (self.vel_rb[:, 1] == 0.0) & (self.vel_rb[:, 2] == 0.0)
        )
        # ADVANCE
        is_advances = ~((self.sangle == 60.0) & (self.coptb == 0) & (self.user == 0) & (self.min_node == 0))
        sangles = array_default_float(self.sangle, default=60.0, size=size, is_double=False)
        coptbs = array_default_int(self.coptb, default=0)
        users = array_default_int(self.user, default=0)
        min_nodes = array_default_int(self.min_node, default=0)

        # RIGID
        cgids = array_default_int(self.min_node, default=-1)

        for bcbody_id, dim, behav, bsid, istype, fric, idispl, control, \
            is_nload, nload, ang_vel, dcos, vel_rb, \
            is_advance, sangle, coptb, user, min_node, \
            cgid, nent, rigid_body_name in zip(
            bcbody_ids, self.dim, self.behavior, self.bsid, istypes, frics, idispls, controls,
            is_nloads, nloads, ang_vels, dcoss, vel_rbs,
            # advance
            is_advances, sangles, coptbs, users, min_nodes,
            # rigid
            cgids, self.nent, self.rigid_body_name):
            if is_nload:
                list_fields = [
                    'BCBODY', bcbody_id, dim, behav, bsid, istype, fric, idispl, control,
                    nload, ang_vel] + dcos + vel_rb
            else:
                list_fields = [
                    'BCBODY', bcbody_id, dim, behav, bsid, istype, fric, idispl, control]
            if is_advance:
                list_fields += ['ADVANCE', sangle, coptb, user, min_node, '', '']

            msgi = ''
            if rigid_body_name:
                msgi = f'        {"RIGID":<8s}{cgid:>8s}{nent:>8d}{rigid_body_name:>24s}\n'
                #list_fields += ['RIGID', cgid, nent, rigid_body_name]

            #list_fields = ['BEDGE', bedge_id]
            #element_id = element_id_[ieid0:ieid1]
            #nodes = nodes_[ieid0:ieid1, :].tolist()
            #for eid, (nid1, nid2) in zip(element_id, nodes):
                #list_fields.extend([eid, nid1, nid2, None])
            bdf_file.write(print_card(list_fields) + msgi)
        return


class BCBODY1(VectorizedBaseCard):
    """
    Format: (SOLs 101 and 400 only)
    +---------+-----+------+-----+-------+------+--------+-------+
    |    1    |  2  |   3  |  4  |   5   |   6  |   7    |   8   |
    +=========+=====+======+=====+=======+======+========+=======+
    | BCBODY1 | BID | BPID | DIM | BEHAV | BSID | BCRGID | BCGOUT|
    +---------+-----+------+-----+-------+------+--------+-------+

    Format: (SOL 700 only)
    +---------+-----+------+-----+-------+------+--------+-------+
    |    1    |  2  |   3  |  4  |   5   |   6  |   7    |   8   |
    +=========+=====+======+=====+=======+======+========+=======+
    | BCBODY1 | BID |      |     |       | BSID |        |       |
    +---------+-----+------+-----+-------+------+--------+-------+

    """
    _id_name = 'bcbody_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.bcbody_id = np.array([], dtype='int32')

    def add(self, bcbody: int,
            eids: list[int],
            grids: list[tuple[int, int]],
            comment: str='') -> int:
        """Creates a BCBODY1 card

        BID / bcbody_id : int
            Parameter identification number of a BCBDPRP entry. (Integer > 0 or blank) Ignored in
            SOL 700. See Remark 2.
        BPID : int
            Parameter identification number of a BCBDPRP entry. (Integer > 0 or blank) Ignored in
            SOL 700. See Remark 2.
        DIM : str; default='3D'
            Dimension of body. (Character; Default= 3D) Ignored in SOL 700.
            DIM=2D: planar body in x-y plane of the basic coordinate system, composed of 2D
                    elements or curves.
            DIM=3D: Any 3D body composed of rigid surfaces, shell elements or solid elements.
        BEHAV : str
            Behavior of curve or surface (Character; Default = DEFORM) Ignored in SOL 700.
            DEFORM body is deformable
            RIGID body is rigid (See Remark 1.)
            SYMM body is a symmetry rigid body
            HEAT indicates body is a heat-rigid body (See Remark 3..).
        """
        afsd
        #self.cards.append((bedge_id, eids, grids, comment))
        #self.n += 1
        #return self.n - 1

    #def remove_unused(self)
    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a BCBODY1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        BID : int (4,1)
           Contact body identification number referenced by
           BCTABLE, BCHANGE, or BCMOVE. (Integer > 0; Required)
        DIM : str; default='3D'
           Dimension of body.
           DIM=2D planar body in x-y plane of the basic coordinate system,
                  composed of 2D elements or curves.
           DIM=3D any 3D body composed of rigid surfaces, shell elements or solid
                  elements.
        BEHAV (4,8)
           Behavior of curve or surface (Character; Default = DEFORM) DEFORM body is
           deformable, RIGID body is rigid, SYMM body is a symmetry body, ACOUS
           indicates an acoustic body, WORK indicates body is a workpiece, HEAT indicates
           body is a heat-rigid body. See Remark 3. for Rigid Bodies..
        BSID : int
            Identification number of a BSURF, BCBOX, BCPROP or BCMATL entry if
            BEHAV=DEFORM. (Integer > 0)
        ISTYP : int (4,3)
           Check of contact conditions. (Integer > 0; Default = 0)
        ISTYP : int
           is not supported in segment-to-segment contact.
           For a deformable body:
           =0 symmetric penetration, double sided contact.
           =1 unsymmetric penetration, single sided contact. (Integer > 0)
           =2 double-sided contact with automatic optimization of contact constraint
              equations (this option is known as “optimized contact”).
              Notes: single-sided contact (ISTYP=1) with the contact bodies arranged properly
              using the contact table frequently runs much faster than ISTYP=2.
              For a rigid body:
           =0 no symmetry condition on rigid body.
           =1 rigid body is a symmetry plane.
        FRIC : int/float (6,7)
            Friction coefficient. (Real > 0 or integer; Default = 0)
            If the value is an integer it represents the ID of a TABL3Di.
        IDSPL : int (4,5)
            Set IDSPL=1 to activate the SPLINE (analytical contact) option for a deformable
            body and for a rigid contact surface. Set it to zero or leave blank to not have
            analytical contact. (Integer; Default = 0)
        NLOAD : int or None
            Enter a positive number if "load controlled" and rotations are allowed (Integer). The
            positive number is the grid number where the moments or rotations are applied. The
            rotations are specified using SPCD at grid ID NLOAD and can be specified using dof's
            1-3 (for rotation about x, y, z respectively), or by dof's 4-6 (for rotation about x, y, z
            respectively).
            Note: This rotation takes the position of the grid point defined in CGID field as the
            center of rotation.
        ANGVEL : int/float; default=0.0
            Angular velocity or angular position about local axis through center of rotation. If the
            value is an integer it represents the ID of a TABLED1, TABLED2 or TABL3D, i.e., a
            time-dependent or multi-dimensional table; however, no log scales, only linear scales.
            (Real or Integer; Default = 0.0)
        DCOSi : int/float; default=0.0
            Components of direction cosine of local axis if ANGVEL is nonzero. If the value is an
            integer, it represents the ID of a TABLED1, TABLED2 or TABL3D, i.e., a time-dependent
            or multi-dimensional table; however, no log scales, only linear scales. (Real
            or Integer; Default=0.0) In 2D contact only DCOS3 is used and the Default is 1.0.
        VELRBi : int/float; default=0.0
            Translation velocity or final position (depending on the value of CONTROL) of rigid
            body at the grid point defined in CGID filed. For velocity control only, if the value is
            an integer, it represents the ID of TABLED1, TABLED2 or TABL3D, i.e., a time-dependent
            or multi-dimensional table; however, no log scales, only linear scales. Only
            VELRB1 and VELRB2 are used in 2D contact. (Real or Integer; Default = 0.0)

        """


        #['BCBODY1', '1', '5001', '3D', 'DEFORM', '5']
        bcbody_id = integer(card, 1, 'bcbody_id')
        bpid = integer_or_blank(card, 2, 'bpid', default=0)
        #dimension = string_or_blank(card, 3, 'dim', default='3D')
        dimension = string_choice_or_blank(card, 3, 'dim',
                                           ('2D', '3D'),
                                           default='3D')

        behavior = string_or_blank(card, 4, 'behavior', default='DEFORM')

        assert behavior in {'DEFORM', 'RIGID', 'SYMM', 'HEAT'}, behavior

        #BSID : int
        #    For SOLs 101 and 400: Identification number of a BSURF or BCPROP entry if
        #    BEHAV=DEFORM or HEAT, or identification number of a BCRGSRF, BCPATCH,
        #    BCBZIER, BCNURB2, or BCNURBS entry if BEHAV=RIGID or SYMM (See Remark 4.).
        #    For SOL 700: Identification number (RBID) of a BSURF, BCBOX, BCPROP, BCMATL,
        #    BCSEG, BCGRID or BCELIPS entry. (Integer > 0).
        # BCRGID : int
        #    For SOLs 101 and 400: Identification number of a BCRIGID entry if BEHAV=RIGID or
        #    SYMM. Ignored in SOL 700. (Integer > 0)
        bsid = integer(card, 5, 'bsid')

        bc_rigid = integer_or_blank(card, 6, 'bc_rigid', default=0)

        #BCGOUT Reference point in basic coordinate system to calculate the
        #       global resultant contact force/moment (integer)
        #  -1 Origin
        #   0 Estimated centroid of deformable body (default)
        bcg_out = integer_or_blank(card, 7, 'bcg_out', default=0)
        assert len(card) < 9, card

        cardi = (bcbody_id, bpid, dimension, behavior, bsid, bc_rigid, bcg_out, comment)
        self.cards.append(cardi)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')

        #idtype = self.model.idtype
        #fdtype = self.model.fdtype
        bcbody_id = np.zeros(ncards, dtype='int32')
        bpid = np.zeros(ncards, dtype='int32')
        bsid = np.zeros(ncards, dtype='int32')
        dimension = np.zeros(ncards, dtype='|U8')
        behavior = np.zeros(ncards, dtype='|U8')
        bc_rigid = np.zeros(ncards, dtype='int32')
        bcg_out = np.zeros(ncards, dtype='int32')

        #comment = {}
        for icard, card in enumerate(self.cards):
            (bcbody_idi, bpidi, dimensioni, behaviori, bsidi, bc_rigidi, bcg_outi, comment) = card
            assert isinstance(bsidi, int), bsidi
            bcbody_id[icard] = bcbody_idi
            bpid[icard] = bpidi
            bsid[icard] = bsidi
            dimension[icard] = dimensioni
            behavior[icard] = behaviori
            bc_rigid[icard] = bc_rigidi
            bcg_out[icard] = bcg_outi
            #if commenti:
                #comment[i] = commenti
                #comment[nidi] = commenti

        self._save(bcbody_id, bpid, dimension, behavior, bsid, bc_rigid, bcg_out)
        self.sort()
        self.cards = []

    def _save(self, bcbody_id, bpid, dimension, behavior, bsid, bc_rigid, bcg_out) -> None:
        ncards_existing = len(self.bcbody_id)
        if ncards_existing != 0:
            bcbody_id = np.hstack([self.bcbody_id, bcbody_id])
            #nelement = np.hstack([self.nelement, nelement])
            #element_id = np.hstack([self.element_id, element_id])
            #nodes = np.hstack([self.nodes, nodes])
            asdf
        #if comment:
            #self.comment.update(comment)
        self.bcbody_id = bcbody_id
        self.bpid = bpid
        self.dimension = dimension
        self.behavior = behavior
        self.bsid = bsid
        self.bc_rigid = bc_rigid
        self.bcg_out = bcg_out
        self.n = len(self.bcbody_id)

    def __apply_slice__(self, bcbody: BCBODY1, i: np.ndarray) -> None:
        #self._slice_comment(bedge, i)
        bcbody.bcbody_id = self.bcbody_id[i]
        bcbody.bpid = self.bpid[i]
        bcbody.dimension = self.dimension[i]
        bcbody.behavior = self.behavior[i]
        bcbody.bsid = self.bsid[i]
        bcbody.bc_rigid = self.bc_rigid[i]
        bcbody.bcg_out = self.bcg_out[i]
        bcbody.n = len(i)

    #def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        #used_dict['element_id'].append(self.element_id)

    #def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        #"""helper for bdf_equivalence_nodes"""
        #nodes = self.nodes.ravel()
        #for i, nid1 in enumerate(nodes):
            #nid2 = nid_old_to_new.get(nid1, nid1)
            #nodes[i] = nid2

    #def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        #asd
        #contact_regions = np.hstack([self.source_ids, self.target_ids])
        #used_dict['contact_set_id'].append(contact_regions)

    #def convert(self, **kwargs) -> None:
        #pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        """
        Defines a region of element edges for:
         - glue (BGSET entry)
         - contact (BCTSET entry)
         - fluid pressure load (PLOADFP entry)
         - wetted edges for a Co-simulation (CSMSET entry)
        """
        nids = self.model.grid.node_id
        #unids = np.unique(self.nodes.ravel())
        #geom_check(self,
                   #missing,
                   #node=(nids, unids))

    #@property
    #def ielement(self) -> np.ndarray:
        #return make_idim(self.n, self.nelement)

    @property
    def max_id(self) -> int:
        return self.bcbody_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        bcbody_ids = array_str(self.bcbody_id, size=size).tolist()
        bsids = array_str(self.bsid, size=size).tolist()
        bpids = array_default_int(self.bpid, default=0, size=size).tolist()

        # RIGID
        bc_rigids = array_default_int(self.bc_rigid, default=0)
        bcg_outs = array_default_int(self.bcg_out, default=0, size=size).tolist()

        for bcbody_id, bpid, bsid, dim, behav, bc_rigid, bcg_out in zip(
            bcbody_ids, bpids, bsids, self.dimension, self.behavior, bc_rigids, bcg_outs):

            #+---------+-----+------+-----+-------+------+--------+-------+
            #|    1    |  2  |   3  |  4  |   5   |   6  |   7    |   8   |
            #+=========+=====+======+=====+=======+======+========+=======+
            #| BCBODY1 | BID | BPID | DIM | BEHAV | BSID | BCRGID | BCGOUT|
            #+---------+-----+------+-----+-------+------+--------+-------+
            list_fields = [
                'BCBODY1', bcbody_id, bpid, dim, behav, bsid, bc_rigid, bcg_out]
            bdf_file.write(print_card(list_fields))
        return
