from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string,
    integer_or_blank, double_or_blank, string_or_blank,
    integer_string_or_blank,
)
from pyNastran.utils.numpy_utils import integer_types

from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, hslice_by_idim, make_idim,
    parse_check,
)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_default_int, get_print_card_size,
)

if TYPE_CHECKING:
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.bdf.types import TextIOLike


class SETREE(VectorizedBaseCard):
    """
    Superelement Tree Definition

    +--------+------+-------+-------+-------+-------+-------+-------+-------+
    |   1    |   2  |   3   |   4   |   5   |   6   |   7   |   8   |   9   |
    +========+======+=======+=======+=======+=======+=======+=======+=======+
    | SETREE | SEID | SEUP1 | SEUP2 | SEUP3 | SEUP4 | SEUP5 | SEUP6 | SEUP7 |
    +--------+------+-------+-------+-------+-------+-------+-------+-------+
    """
    _id_name = 'seid'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.seid = np.array([], dtype='int32')
        self.nsuperelements = np.array([], dtype='int32')
        self.superelements = np.array([], dtype='int32')

    def add(self, seid: int, superelements: list[int],
            comment: str='') -> int:
        if isinstance(superelements, integer_types):
            superelements = [superelements]
        superelements = expand_thru(superelements, set_fields=False, sort_fields=False)
        self.cards.append((seid, superelements, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        seid = integer(card, 1, 'seid')
        superelements_raw = []
        nfields = len(card)
        for i in range(2, nfields):
            val = integer_string_or_blank(card, i, f'seup{i-1}')
            if val is not None:
                superelements_raw.append(val)
        assert len(superelements_raw) >= 1, f'SETREE card={card}'
        superelements = expand_thru(superelements_raw, set_fields=False, sort_fields=False)
        self.cards.append((seid, superelements, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        seid = np.zeros(ncards, dtype='int32')
        nsuperelements = np.zeros(ncards, dtype='int32')
        all_superelements = []
        for icard, card in enumerate(self.cards):
            (seidi, superelementsi, comment) = card
            seid[icard] = seidi
            nsuperelements[icard] = len(superelementsi)
            all_superelements.extend(superelementsi)
        superelements = np.array(all_superelements, dtype='int32')
        self._save(seid, nsuperelements, superelements)
        self.cards = []

    def _save(self, seid, nsuperelements, superelements):
        if len(self.seid) != 0:
            seid = np.hstack([self.seid, seid])
            nsuperelements = np.hstack([self.nsuperelements, nsuperelements])
            superelements = np.hstack([self.superelements, superelements])
        self.seid = seid
        self.nsuperelements = nsuperelements
        self.superelements = superelements
        self.n = len(seid)

    def __apply_slice__(self, card: SETREE, i: np.ndarray) -> None:
        card.n = len(i)
        card.seid = self.seid[i]
        idim = self.idim
        card.superelements = hslice_by_idim(i, idim, self.superelements)
        card.nsuperelements = self.nsuperelements[i]

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.nsuperelements)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        return max(self.seid.max(), self.superelements.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for seid, (idim0, idim1) in zip(self.seid, self.idim):
            superelements = self.superelements[idim0:idim1].tolist()
            list_fields = ['SETREE', seid] + superelements
            bdf_file.write(print_card(list_fields))


class SENQSET(VectorizedBaseCard):
    """
    Superelement Number of Internally Generated Scalar Points

    +---------+------+---+
    |    1    |   2  | 3 |
    +=========+======+===+
    | SENQSET | SEID | N |
    +---------+------+---+
    """
    _id_name = 'set_id'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.set_id = np.array([], dtype='int32')
        self.n_val = np.array([], dtype='int32')

    def add(self, set_id: int, n: int=0, comment: str='') -> int:
        self.cards.append((set_id, n, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        set_id = integer(card, 1, 'set_id')
        n = integer_or_blank(card, 2, 'n', default=0)
        assert len(card) <= 3, f'len(SENQSET card) = {len(card):d}\ncard={card}'
        self.cards.append((set_id, n, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        set_id = np.zeros(ncards, dtype='int32')
        n_val = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (set_idi, ni, comment) = card
            set_id[icard] = set_idi
            n_val[icard] = ni
        self._save(set_id, n_val)
        self.cards = []

    def _save(self, set_id, n_val):
        if len(self.set_id) != 0:
            set_id = np.hstack([self.set_id, set_id])
            n_val = np.hstack([self.n_val, n_val])
        self.set_id = set_id
        self.n_val = n_val
        self.n = len(set_id)

    def __apply_slice__(self, card: SENQSET, i: np.ndarray) -> None:
        card.n = len(i)
        card.set_id = self.set_id[i]
        card.n_val = self.n_val[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        return self.set_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for set_id, n in zip(self.set_id, self.n_val):
            list_fields = ['SENQSET', set_id, n]
            bdf_file.write(print_card(list_fields))


class SEBULK(VectorizedBaseCard):
    """
    Superelement Bulk Data Definition

    +--------+------+--------+-------+--------+------+------+--------+
    |   1    |   2  |   3    |   4   |   5    |  6   |   7  |   8    |
    +========+======+========+=======+========+======+======+========+
    | SEBULK | SEID |  TYPE  | RSEID | METHOD | TOL  | LOC  | UNITNO |
    +--------+------+--------+-------+--------+------+------+--------+
    """
    _id_name = 'seid'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.seid = np.array([], dtype='int32')
        self.superelement_type = np.array([], dtype='|U8')
        self.rseid = np.array([], dtype='int32')
        self.method = np.array([], dtype='|U8')
        self.tol = np.array([], dtype='float64')
        self.loc = np.array([], dtype='|U4')
        self.unitno = np.array([], dtype='int32')

    def add(self, seid: int, superelement_type: str, rseid: int=0,
            method: str='AUTO', tol: float=1e-5, loc: str='YES',
            unitno: int=0, comment: str='') -> int:
        self.cards.append((seid, superelement_type, rseid, method, tol, loc, unitno, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        seid = integer(card, 1, 'seid')
        superelement_type = string(card, 2, 'superelement_type')
        rseid = integer_or_blank(card, 3, 'rseid', default=0)
        method = string_or_blank(card, 4, 'method', default='AUTO')
        tol = double_or_blank(card, 5, 'tol', default=1e-5)
        loc = string_or_blank(card, 6, 'loc', default='YES')
        unitno = integer_or_blank(card, 7, 'unitno', default=0)
        assert len(card) <= 8, f'len(SEBULK card) = {len(card):d}\ncard={card}'
        self.cards.append((seid, superelement_type, rseid, method, tol, loc, unitno, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        seid = np.zeros(ncards, dtype='int32')
        superelement_type = np.full(ncards, '', dtype='|U8')
        rseid = np.zeros(ncards, dtype='int32')
        method = np.full(ncards, '', dtype='|U8')
        tol = np.zeros(ncards, dtype='float64')
        loc = np.full(ncards, '', dtype='|U4')
        unitno = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (seidi, se_typei, rseidi, methodi, toli, loci, unitnoi, comment) = card
            seid[icard] = seidi
            superelement_type[icard] = se_typei
            rseid[icard] = rseidi
            method[icard] = methodi
            tol[icard] = toli
            loc[icard] = loci
            unitno[icard] = unitnoi if unitnoi is not None else 0
        self._save(seid, superelement_type, rseid, method, tol, loc, unitno)
        self.cards = []

    def _save(self, seid, superelement_type, rseid, method, tol, loc, unitno):
        if len(self.seid) != 0:
            seid = np.hstack([self.seid, seid])
            superelement_type = np.hstack([self.superelement_type, superelement_type])
            rseid = np.hstack([self.rseid, rseid])
            method = np.hstack([self.method, method])
            tol = np.hstack([self.tol, tol])
            loc = np.hstack([self.loc, loc])
            unitno = np.hstack([self.unitno, unitno])
        self.seid = seid
        self.superelement_type = superelement_type
        self.rseid = rseid
        self.method = method
        self.tol = tol
        self.loc = loc
        self.unitno = unitno
        self.n = len(seid)

    def __apply_slice__(self, card: SEBULK, i: np.ndarray) -> None:
        card.n = len(i)
        card.seid = self.seid[i]
        card.superelement_type = self.superelement_type[i]
        card.rseid = self.rseid[i]
        card.method = self.method[i]
        card.tol = self.tol[i]
        card.loc = self.loc[i]
        card.unitno = self.unitno[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        return max(self.seid.max(), self.rseid.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for seid, se_type, rseid, method, tol, loc, unitno in zip(
                self.seid, self.superelement_type, self.rseid,
                self.method, self.tol, self.loc, self.unitno):
            unitno_val = unitno if unitno != 0 else None
            list_fields = ['SEBULK', seid, se_type, rseid, method, tol, loc, unitno_val]
            bdf_file.write(print_card(list_fields))


class SEBNDRY(VectorizedBaseCard):
    """
    Superelement Boundary Grid Point Definition

    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    |    1    |   2   |   3   |   4   |   5   |   6   |   7   |   8   |   9   |
    +=========+=======+=======+=======+=======+=======+=======+=======+=======+
    | SEBNDRY | SEIDA | SEIDB | GIDA1 | GIDA2 | GIDA3 | GIDA4 | GIDA5 | GIDA6 |
    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    """
    _id_name = 'seid_a'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.seid_a = np.array([], dtype='int32')
        self.seid_b = np.array([], dtype='int32')
        self.nids = np.array([], dtype='int32')
        self.nnodes = np.array([], dtype='int32')

    def add(self, seid_a: int, seid_b: int, ids: list[int],
            comment: str='') -> int:
        if isinstance(ids, integer_types):
            ids = [ids]
        ids = expand_thru(ids, set_fields=False, sort_fields=False)
        self.cards.append((seid_a, seid_b, ids, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        seid_a = integer(card, 1, 'seid_a')
        seid_b_raw = integer_string_or_blank(card, 2, 'seid_b')
        if seid_b_raw is None or seid_b_raw == 'ALL':
            seid_b = -1  # encode 'ALL' as -1
        else:
            seid_b = seid_b_raw
        nfields = len(card)
        ids_raw = []
        for i in range(3, nfields):
            val = integer_string_or_blank(card, i, f'gida{i-2}')
            if val is not None:
                ids_raw.append(val)
        assert len(ids_raw) >= 1, f'SEBNDRY card={card}'
        ids = expand_thru(ids_raw, set_fields=False, sort_fields=False)
        self.cards.append((seid_a, seid_b, ids, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        seid_a = np.zeros(ncards, dtype='int32')
        seid_b = np.zeros(ncards, dtype='int32')
        nnodes = np.zeros(ncards, dtype='int32')
        all_nids = []
        for icard, card in enumerate(self.cards):
            (seid_ai, seid_bi, idsi, comment) = card
            seid_a[icard] = seid_ai
            seid_b[icard] = seid_bi
            nnodes[icard] = len(idsi)
            all_nids.extend(idsi)
        nids = np.array(all_nids, dtype=idtype)
        self._save(seid_a, seid_b, nids, nnodes)
        self.cards = []

    def _save(self, seid_a, seid_b, nids, nnodes):
        if len(self.seid_a) != 0:
            seid_a = np.hstack([self.seid_a, seid_a])
            seid_b = np.hstack([self.seid_b, seid_b])
            nids = np.hstack([self.nids, nids])
            nnodes = np.hstack([self.nnodes, nnodes])
        self.seid_a = seid_a
        self.seid_b = seid_b
        self.nids = nids
        self.nnodes = nnodes
        self.n = len(seid_a)

    def __apply_slice__(self, card: SEBNDRY, i: np.ndarray) -> None:
        card.n = len(i)
        card.seid_a = self.seid_a[i]
        card.seid_b = self.seid_b[i]
        idim = self.idim
        card.nids = hslice_by_idim(i, idim, self.nids)
        card.nnodes = self.nnodes[i]

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.nnodes)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.nids)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self, missing, node=(nid, self.nids))

    @property
    def max_id(self) -> int:
        return max(self.seid_a.max(), self.nids.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for seid_a, seid_b, (idim0, idim1) in zip(
                self.seid_a, self.seid_b, self.idim):
            nids = self.nids[idim0:idim1].tolist()
            seid_b_val = 'ALL' if seid_b == -1 else seid_b
            list_fields = ['SEBNDRY', seid_a, seid_b_val] + nids
            bdf_file.write(print_card(list_fields))


class SECONCT(VectorizedBaseCard):
    """
    Superelement Boundary Connection

    +---------+-------+-------+------+------+------+------+------+------+
    |    1    |   2   |   3   |  4   |   5  |  6   |  7   |  8   |  9   |
    +=========+=======+=======+======+======+======+======+======+======+
    | SECONCT | SEIDA | SEIDB | TOL  | LOC  |      |      |      |      |
    +---------+-------+-------+------+------+------+------+------+------+
    |         | GIDA1 | GIDB1 | GIDA2| GIDB2| GIDA3| GIDB3|      |      |
    +---------+-------+-------+------+------+------+------+------+------+
    """
    _id_name = 'seid_a'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.seid_a = np.array([], dtype='int32')
        self.seid_b = np.array([], dtype='int32')
        self.tol = np.array([], dtype='float64')
        self.loc = np.array([], dtype='|U4')
        self.nodes_a = np.array([], dtype='int32')
        self.nodes_b = np.array([], dtype='int32')
        self.nnodes = np.array([], dtype='int32')

    def add(self, seid_a: int, seid_b: int, tol: float=1e-5,
            loc: str='YES', nodes_a: list[int]=None,
            nodes_b: list[int]=None, comment: str='') -> int:
        if nodes_a is None:
            nodes_a = []
        if nodes_b is None:
            nodes_b = []
        self.cards.append((seid_a, seid_b, tol, loc, nodes_a, nodes_b, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        seid_a = integer(card, 1, 'seid_a')
        seid_b = integer(card, 2, 'seid_b')
        tol = double_or_blank(card, 3, 'tol', default=1e-5)
        loc = string_or_blank(card, 4, 'loc', default='YES')
        # fields 5-8 are blank
        nfields = len(card)
        nodes_a = []
        nodes_b = []
        i = 9
        while i < nfields:
            nid_a = integer_or_blank(card, i, f'gida{len(nodes_a)+1}')
            nid_b = integer_or_blank(card, i + 1, f'gidb{len(nodes_b)+1}')
            if nid_a is not None and nid_b is not None:
                nodes_a.append(nid_a)
                nodes_b.append(nid_b)
            i += 2
        self.cards.append((seid_a, seid_b, tol, loc, nodes_a, nodes_b, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        seid_a = np.zeros(ncards, dtype='int32')
        seid_b = np.zeros(ncards, dtype='int32')
        tol = np.zeros(ncards, dtype='float64')
        loc = np.full(ncards, '', dtype='|U4')
        nnodes = np.zeros(ncards, dtype='int32')
        all_nodes_a = []
        all_nodes_b = []
        for icard, card in enumerate(self.cards):
            (seid_ai, seid_bi, toli, loci, nodes_ai, nodes_bi, comment) = card
            seid_a[icard] = seid_ai
            seid_b[icard] = seid_bi
            tol[icard] = toli
            loc[icard] = loci
            nnodes[icard] = len(nodes_ai)
            all_nodes_a.extend(nodes_ai)
            all_nodes_b.extend(nodes_bi)
        nodes_a = np.array(all_nodes_a, dtype=idtype)
        nodes_b = np.array(all_nodes_b, dtype=idtype)
        self._save(seid_a, seid_b, tol, loc, nodes_a, nodes_b, nnodes)
        self.cards = []

    def _save(self, seid_a, seid_b, tol, loc, nodes_a, nodes_b, nnodes):
        if len(self.seid_a) != 0:
            seid_a = np.hstack([self.seid_a, seid_a])
            seid_b = np.hstack([self.seid_b, seid_b])
            tol = np.hstack([self.tol, tol])
            loc = np.hstack([self.loc, loc])
            nodes_a = np.hstack([self.nodes_a, nodes_a])
            nodes_b = np.hstack([self.nodes_b, nodes_b])
            nnodes = np.hstack([self.nnodes, nnodes])
        self.seid_a = seid_a
        self.seid_b = seid_b
        self.tol = tol
        self.loc = loc
        self.nodes_a = nodes_a
        self.nodes_b = nodes_b
        self.nnodes = nnodes
        self.n = len(seid_a)

    def __apply_slice__(self, card: SECONCT, i: np.ndarray) -> None:
        card.n = len(i)
        card.seid_a = self.seid_a[i]
        card.seid_b = self.seid_b[i]
        card.tol = self.tol[i]
        card.loc = self.loc[i]
        idim = self.idim
        card.nodes_a = hslice_by_idim(i, idim, self.nodes_a)
        card.nodes_b = hslice_by_idim(i, idim, self.nodes_b)
        card.nnodes = self.nnodes[i]

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.nnodes)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.nodes_a)
        used_dict['node_id'].append(self.nodes_b)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        all_nodes = np.hstack([self.nodes_a, self.nodes_b])
        geom_check(self, missing, node=(nid, all_nodes))

    @property
    def max_id(self) -> int:
        max_val = self.seid_a.max()
        if len(self.nodes_a) > 0:
            max_val = max(max_val, self.nodes_a.max(), self.nodes_b.max())
        return max_val

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for seid_a, seid_b, tol, loc, (idim0, idim1) in zip(
                self.seid_a, self.seid_b, self.tol, self.loc, self.idim):
            nodes_a = self.nodes_a[idim0:idim1]
            nodes_b = self.nodes_b[idim0:idim1]
            list_fields = ['SECONCT', seid_a, seid_b, tol, loc,
                           None, None, None, None]
            for nid_a, nid_b in zip(nodes_a, nodes_b):
                list_fields += [nid_a, nid_b]
            bdf_file.write(print_card(list_fields))


class SEELT(VectorizedBaseCard):
    """
    Superelement Element Connection

    +-------+------+------+------+------+------+------+------+------+
    |   1   |   2  |   3  |   4  |   5  |   6  |   7  |   8  |   9  |
    +=======+======+======+======+======+======+======+======+======+
    | SEELT | SEID | EID1 | EID2 | EID3 | EID4 | EID5 | EID6 | EID7 |
    +-------+------+------+------+------+------+------+------+------+
    """
    _id_name = 'seid'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.seid = np.array([], dtype='int32')
        self.element_ids = np.array([], dtype='int32')
        self.neids = np.array([], dtype='int32')

    def add(self, seid: int, eids: list[int], comment: str='') -> int:
        if isinstance(eids, integer_types):
            eids = [eids]
        eids = expand_thru(eids, set_fields=False, sort_fields=False)
        self.cards.append((seid, eids, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        seid = integer(card, 1, 'seid')
        nfields = len(card)
        eids_raw = []
        for i in range(2, nfields):
            val = integer_string_or_blank(card, i, f'eid{i-1}')
            if val is not None:
                eids_raw.append(val)
        eids = expand_thru(eids_raw, set_fields=False, sort_fields=False)
        self.cards.append((seid, eids, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        seid = np.zeros(ncards, dtype='int32')
        neids = np.zeros(ncards, dtype='int32')
        all_eids = []
        for icard, card in enumerate(self.cards):
            (seidi, eidsi, comment) = card
            seid[icard] = seidi
            neids[icard] = len(eidsi)
            all_eids.extend(eidsi)
        element_ids = np.array(all_eids, dtype=idtype)
        self._save(seid, element_ids, neids)
        self.cards = []

    def _save(self, seid, element_ids, neids):
        if len(self.seid) != 0:
            seid = np.hstack([self.seid, seid])
            element_ids = np.hstack([self.element_ids, element_ids])
            neids = np.hstack([self.neids, neids])
        self.seid = seid
        self.element_ids = element_ids
        self.neids = neids
        self.n = len(seid)

    def __apply_slice__(self, card: SEELT, i: np.ndarray) -> None:
        card.n = len(i)
        card.seid = self.seid[i]
        idim = self.idim
        card.element_ids = hslice_by_idim(i, idim, self.element_ids)
        card.neids = self.neids[i]

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.neids)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['element_id'].append(self.element_ids)

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        return max(self.seid.max(), self.element_ids.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for seid, (idim0, idim1) in zip(self.seid, self.idim):
            eids = self.element_ids[idim0:idim1].tolist()
            list_fields = ['SEELT', seid] + eids
            bdf_file.write(print_card(list_fields))


class SELOC(VectorizedBaseCard):
    """
    Superelement Location

    +-------+------+-----+-----+-----+-----+-----+-----+
    |   1   |   2  |  3  |  4  |  5  |  6  |  7  |  8  |
    +=======+======+=====+=====+=====+=====+=====+=====+
    | SELOC | SEID | PA1 | PA2 | PA3 | PB1 | PB2 | PB3 |
    +-------+------+-----+-----+-----+-----+-----+-----+
    """
    _id_name = 'seid'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.seid = np.array([], dtype='int32')
        self.nodes_seid = np.zeros((0, 3), dtype='int32')
        self.nodes_0 = np.zeros((0, 3), dtype='int32')

    def add(self, seid: int, nodes_seid: list[int], nodes_0: list[int],
            comment: str='') -> int:
        self.cards.append((seid, nodes_seid, nodes_0, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        seid = integer(card, 1, 'seid')
        nodes_seid = [integer(card, 2, 'PA1'),
                      integer(card, 3, 'PA2'),
                      integer(card, 4, 'PA3')]
        nodes_0 = [integer(card, 5, 'PB1'),
                   integer(card, 6, 'PB2'),
                   integer(card, 7, 'PB3')]
        assert len(card) <= 8, f'len(SELOC card) = {len(card):d}\ncard={card}'
        self.cards.append((seid, nodes_seid, nodes_0, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        seid = np.zeros(ncards, dtype='int32')
        nodes_seid = np.zeros((ncards, 3), dtype=idtype)
        nodes_0 = np.zeros((ncards, 3), dtype=idtype)
        for icard, card in enumerate(self.cards):
            (seidi, nodes_seidi, nodes_0i, comment) = card
            seid[icard] = seidi
            nodes_seid[icard, :] = nodes_seidi
            nodes_0[icard, :] = nodes_0i
        self._save(seid, nodes_seid, nodes_0)
        self.cards = []

    def _save(self, seid, nodes_seid, nodes_0):
        if len(self.seid) != 0:
            seid = np.hstack([self.seid, seid])
            nodes_seid = np.vstack([self.nodes_seid, nodes_seid])
            nodes_0 = np.vstack([self.nodes_0, nodes_0])
        self.seid = seid
        self.nodes_seid = nodes_seid
        self.nodes_0 = nodes_0
        self.n = len(seid)

    def __apply_slice__(self, card: SELOC, i: np.ndarray) -> None:
        card.n = len(i)
        card.seid = self.seid[i]
        card.nodes_seid = self.nodes_seid[i, :]
        card.nodes_0 = self.nodes_0[i, :]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.nodes_seid.ravel())
        used_dict['node_id'].append(self.nodes_0.ravel())

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        all_nodes = np.hstack([self.nodes_seid.ravel(), self.nodes_0.ravel()])
        geom_check(self, missing, node=(nid, all_nodes))

    @property
    def max_id(self) -> int:
        return max(self.seid.max(), self.nodes_seid.max(), self.nodes_0.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for seid, nodes_s, nodes0 in zip(self.seid, self.nodes_seid, self.nodes_0):
            list_fields = ['SELOC', seid] + nodes_s.tolist() + nodes0.tolist()
            bdf_file.write(print_card(list_fields))


class SEMPLN(VectorizedBaseCard):
    """
    Superelement Mirror Plane

    +--------+------+-------+----+----+----+
    |   1    |   2  |   3   |  4 |  5 |  6 |
    +========+======+=======+====+====+====+
    | SEMPLN | SEID | PLANE | P1 | P2 | P3 |
    +--------+------+-------+----+----+----+
    """
    _id_name = 'seid'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.seid = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 3), dtype='int32')

    def add(self, seid: int, nodes: list[int], comment: str='') -> int:
        self.cards.append((seid, nodes, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        seid = integer(card, 1, 'seid')
        plane = string(card, 2, 'plane')
        assert plane == 'PLANE', f'SEMPLN plane={plane!r}'
        p1 = integer(card, 3, 'p1')
        p2 = integer(card, 4, 'p2')
        p3 = integer(card, 5, 'p3')
        assert len(card) <= 6, f'len(SEMPLN card) = {len(card):d}\ncard={card}'
        self.cards.append((seid, [p1, p2, p3], comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        seid = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 3), dtype=idtype)
        for icard, card in enumerate(self.cards):
            (seidi, nodesi, comment) = card
            seid[icard] = seidi
            nodes[icard, :] = nodesi
        self._save(seid, nodes)
        self.cards = []

    def _save(self, seid, nodes):
        if len(self.seid) != 0:
            seid = np.hstack([self.seid, seid])
            nodes = np.vstack([self.nodes, nodes])
        self.seid = seid
        self.nodes = nodes
        self.n = len(seid)

    def __apply_slice__(self, card: SEMPLN, i: np.ndarray) -> None:
        card.n = len(i)
        card.seid = self.seid[i]
        card.nodes = self.nodes[i, :]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.nodes.ravel())

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self, missing, node=(nid, self.nodes.ravel()))

    @property
    def max_id(self) -> int:
        return max(self.seid.max(), self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for seid, nodes in zip(self.seid, self.nodes):
            list_fields = ['SEMPLN', seid, 'PLANE'] + nodes.tolist()
            bdf_file.write(print_card(list_fields))


class SELABEL(VectorizedBaseCard):
    """
    Superelement Label

    +---------+------+-------------------------------+
    |    1    |   2  |  3-9                          |
    +=========+======+===============================+
    | SELABEL | SEID | LABEL                         |
    +---------+------+-------------------------------+
    """
    _id_name = 'seid'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.seid = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U56')

    def add(self, seid: int, label: str, comment: str='') -> int:
        self.cards.append((seid, label, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        seid = integer(card, 1, 'seid')
        # label spans remaining fields
        nfields = len(card)
        label_parts = []
        for i in range(2, nfields):
            field = card.field(i)
            if field is None:
                label_parts.append(' ' * 8)
            else:
                label_parts.append(str(field))
        label = ''.join(label_parts).rstrip()
        self.cards.append((seid, label, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        seid = np.zeros(ncards, dtype='int32')
        label = np.full(ncards, '', dtype='|U56')
        for icard, card in enumerate(self.cards):
            (seidi, labeli, comment) = card
            seid[icard] = seidi
            label[icard] = labeli
        self._save(seid, label)
        self.cards = []

    def _save(self, seid, label):
        if len(self.seid) != 0:
            seid = np.hstack([self.seid, seid])
            label = np.hstack([self.label, label])
        self.seid = seid
        self.label = label
        self.n = len(seid)

    def __apply_slice__(self, card: SELABEL, i: np.ndarray) -> None:
        card.n = len(i)
        card.seid = self.seid[i]
        card.label = self.label[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        return self.seid.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        for seid, label in zip(self.seid, self.label):
            bdf_file.write(f'SELABEL {seid:<8s}{label}\n' if isinstance(seid, str)
                           else f'SELABEL {seid:<8d}{label}\n')


class SEEXCLD(VectorizedBaseCard):
    """
    Superelement Exclusion

    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    |    1    |   2   |   3   |   4   |   5   |   6   |   7   |   8   |   9   |
    +=========+=======+=======+=======+=======+=======+=======+=======+=======+
    | SEEXCLD | SEIDA | SEIDB | GIDA1 | GIDA2 | GIDA3 | GIDA4 | GIDA5 | GIDA6 |
    +---------+-------+-------+-------+-------+-------+-------+-------+-------+
    """
    _id_name = 'seid_a'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.seid_a = np.array([], dtype='int32')
        self.seid_b = np.array([], dtype='int32')
        self.nodes = np.array([], dtype='int32')
        self.nnodes = np.array([], dtype='int32')

    def add(self, seid_a: int, seid_b: int, nodes: list[int],
            comment: str='') -> int:
        if isinstance(nodes, integer_types):
            nodes = [nodes]
        nodes = expand_thru(nodes, set_fields=False, sort_fields=False)
        self.cards.append((seid_a, seid_b, nodes, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        seid_a = integer(card, 1, 'seid_a')
        seid_b_raw = integer_string_or_blank(card, 2, 'seid_b')
        if seid_b_raw is None or seid_b_raw == 'ALL':
            seid_b = -1
        else:
            seid_b = seid_b_raw
        nfields = len(card)
        nodes_raw = []
        for i in range(3, nfields):
            val = integer_string_or_blank(card, i, f'gida{i-2}')
            if val is not None:
                nodes_raw.append(val)
        assert len(nodes_raw) >= 1, f'SEEXCLD card={card}'
        nodes = expand_thru(nodes_raw, set_fields=False, sort_fields=False)
        self.cards.append((seid_a, seid_b, nodes, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        seid_a = np.zeros(ncards, dtype='int32')
        seid_b = np.zeros(ncards, dtype='int32')
        nnodes = np.zeros(ncards, dtype='int32')
        all_nodes = []
        for icard, card in enumerate(self.cards):
            (seid_ai, seid_bi, nodesi, comment) = card
            seid_a[icard] = seid_ai
            seid_b[icard] = seid_bi
            nnodes[icard] = len(nodesi)
            all_nodes.extend(nodesi)
        nodes = np.array(all_nodes, dtype=idtype)
        self._save(seid_a, seid_b, nodes, nnodes)
        self.cards = []

    def _save(self, seid_a, seid_b, nodes, nnodes):
        if len(self.seid_a) != 0:
            seid_a = np.hstack([self.seid_a, seid_a])
            seid_b = np.hstack([self.seid_b, seid_b])
            nodes = np.hstack([self.nodes, nodes])
            nnodes = np.hstack([self.nnodes, nnodes])
        self.seid_a = seid_a
        self.seid_b = seid_b
        self.nodes = nodes
        self.nnodes = nnodes
        self.n = len(seid_a)

    def __apply_slice__(self, card: SEEXCLD, i: np.ndarray) -> None:
        card.n = len(i)
        card.seid_a = self.seid_a[i]
        card.seid_b = self.seid_b[i]
        idim = self.idim
        card.nodes = hslice_by_idim(i, idim, self.nodes)
        card.nnodes = self.nnodes[i]

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.nnodes)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.nodes)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self, missing, node=(nid, self.nodes))

    @property
    def max_id(self) -> int:
        max_val = self.seid_a.max()
        if len(self.nodes) > 0:
            max_val = max(max_val, self.nodes.max())
        return max_val

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for seid_a, seid_b, (idim0, idim1) in zip(
                self.seid_a, self.seid_b, self.idim):
            nodes = self.nodes[idim0:idim1].tolist()
            seid_b_val = 'ALL' if seid_b == -1 else seid_b
            list_fields = ['SEEXCLD', seid_a, seid_b_val] + nodes
            bdf_file.write(print_card(list_fields))


class CSUPER(VectorizedBaseCard):
    """
    Secondary Superelement Connection

    +--------+------+------+-----+-----+-----+-----+-----+-----+
    |   1    |   2  |   3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+======+======+=====+=====+=====+=====+=====+=====+
    | CSUPER | SSID | PSID | GP1 | GP2 | GP3 | GP4 | GP5 | GP6 |
    +--------+------+------+-----+-----+-----+-----+-----+-----+
    """
    _id_name = 'seid'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.seid = np.array([], dtype='int32')
        self.psid = np.array([], dtype='int32')
        self.nodes = np.array([], dtype='int32')
        self.nnodes = np.array([], dtype='int32')

    def add(self, seid: int, psid: int, nodes: list[int],
            comment: str='') -> int:
        if isinstance(nodes, integer_types):
            nodes = [nodes]
        nodes = expand_thru(nodes, set_fields=False, sort_fields=False)
        self.cards.append((seid, psid, nodes, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        seid = integer(card, 1, 'seid')
        psid = integer_or_blank(card, 2, 'psid', default=0)
        nfields = len(card)
        nodes_raw = []
        for i in range(3, nfields):
            val = integer_string_or_blank(card, i, f'gp{i-2}')
            if val is not None:
                nodes_raw.append(val)
        nodes = expand_thru(nodes_raw, set_fields=False, sort_fields=False)
        self.cards.append((seid, psid, nodes, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        seid = np.zeros(ncards, dtype='int32')
        psid = np.zeros(ncards, dtype='int32')
        nnodes = np.zeros(ncards, dtype='int32')
        all_nodes = []
        for icard, card in enumerate(self.cards):
            (seidi, psidi, nodesi, comment) = card
            seid[icard] = seidi
            psid[icard] = psidi
            nnodes[icard] = len(nodesi)
            all_nodes.extend(nodesi)
        nodes = np.array(all_nodes, dtype=idtype)
        self._save(seid, psid, nodes, nnodes)
        self.cards = []

    def _save(self, seid, psid, nodes, nnodes):
        if len(self.seid) != 0:
            seid = np.hstack([self.seid, seid])
            psid = np.hstack([self.psid, psid])
            nodes = np.hstack([self.nodes, nodes])
            nnodes = np.hstack([self.nnodes, nnodes])
        self.seid = seid
        self.psid = psid
        self.nodes = nodes
        self.nnodes = nnodes
        self.n = len(seid)

    def __apply_slice__(self, card: CSUPER, i: np.ndarray) -> None:
        card.n = len(i)
        card.seid = self.seid[i]
        card.psid = self.psid[i]
        idim = self.idim
        card.nodes = hslice_by_idim(i, idim, self.nodes)
        card.nnodes = self.nnodes[i]

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.nnodes)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.nodes)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self, missing, node=(nid, self.nodes))

    @property
    def max_id(self) -> int:
        max_val = self.seid.max()
        if len(self.nodes) > 0:
            max_val = max(max_val, self.nodes.max())
        return max_val

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for seid, psid, (idim0, idim1) in zip(self.seid, self.psid, self.idim):
            nodes = self.nodes[idim0:idim1].tolist()
            list_fields = ['CSUPER', seid, psid] + nodes
            bdf_file.write(print_card(list_fields))


class CSUPEXT(VectorizedBaseCard):
    """
    Superelement Exterior Point Definition

    +---------+------+-----+-----+-----+-----+-----+-----+-----+
    |    1    |   2  |  3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +=========+======+=====+=====+=====+=====+=====+=====+=====+
    | CSUPEXT | SEID | GP1 | GP2 | GP3 | GP4 | GP5 | GP6 | GP7 |
    +---------+------+-----+-----+-----+-----+-----+-----+-----+
    """
    _id_name = 'seid'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.seid = np.array([], dtype='int32')
        self.nodes = np.array([], dtype='int32')
        self.nnodes = np.array([], dtype='int32')

    def add(self, seid: int, nodes: list[int], comment: str='') -> int:
        if isinstance(nodes, integer_types):
            nodes = [nodes]
        nodes = expand_thru(nodes, set_fields=False, sort_fields=False)
        self.cards.append((seid, nodes, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        seid = integer(card, 1, 'seid')
        nfields = len(card)
        nodes_raw = []
        for i in range(2, nfields):
            val = integer_string_or_blank(card, i, f'gp{i-1}')
            if val is not None:
                nodes_raw.append(val)
        nodes = expand_thru(nodes_raw, set_fields=False, sort_fields=False)
        self.cards.append((seid, nodes, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        seid = np.zeros(ncards, dtype='int32')
        nnodes = np.zeros(ncards, dtype='int32')
        all_nodes = []
        for icard, card in enumerate(self.cards):
            (seidi, nodesi, comment) = card
            seid[icard] = seidi
            nnodes[icard] = len(nodesi)
            all_nodes.extend(nodesi)
        nodes = np.array(all_nodes, dtype=idtype)
        self._save(seid, nodes, nnodes)
        self.cards = []

    def _save(self, seid, nodes, nnodes):
        if len(self.seid) != 0:
            seid = np.hstack([self.seid, seid])
            nodes = np.hstack([self.nodes, nodes])
            nnodes = np.hstack([self.nnodes, nnodes])
        self.seid = seid
        self.nodes = nodes
        self.nnodes = nnodes
        self.n = len(seid)

    def __apply_slice__(self, card: CSUPEXT, i: np.ndarray) -> None:
        card.n = len(i)
        card.seid = self.seid[i]
        idim = self.idim
        card.nodes = hslice_by_idim(i, idim, self.nodes)
        card.nnodes = self.nnodes[i]

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.nnodes)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.nodes)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self, missing, node=(nid, self.nodes))

    @property
    def max_id(self) -> int:
        max_val = self.seid.max()
        if len(self.nodes) > 0:
            max_val = max(max_val, self.nodes.max())
        return max_val

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for seid, (idim0, idim1) in zip(self.seid, self.idim):
            nodes = self.nodes[idim0:idim1].tolist()
            list_fields = ['CSUPEXT', seid] + nodes
            bdf_file.write(print_card(list_fields))


class SELOAD(VectorizedBaseCard):
    """
    Superelement Load Transfer

    +--------+-------+------+-------+
    |    1   |   2   |   3  |   4   |
    +========+=======+======+=======+
    | SELOAD | LID_S0| SEID | LID_SE|
    +--------+-------+------+-------+
    """
    _id_name = 'lid_s0'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.lid_s0 = np.array([], dtype='int32')
        self.seid = np.array([], dtype='int32')
        self.lid_se = np.array([], dtype='int32')

    def add(self, lid_s0: int, seid: int, lid_se: int,
            comment: str='') -> int:
        self.cards.append((lid_s0, seid, lid_se, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        lid_s0 = integer(card, 1, 'lid_s0')
        seid = integer(card, 2, 'seid')
        lid_se = integer(card, 3, 'lid_se')
        assert len(card) <= 4, f'len(SELOAD card) = {len(card):d}\ncard={card}'
        self.cards.append((lid_s0, seid, lid_se, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        lid_s0 = np.zeros(ncards, dtype='int32')
        seid = np.zeros(ncards, dtype='int32')
        lid_se = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (lid_s0i, seidi, lid_sei, comment) = card
            lid_s0[icard] = lid_s0i
            seid[icard] = seidi
            lid_se[icard] = lid_sei
        self._save(lid_s0, seid, lid_se)
        self.cards = []

    def _save(self, lid_s0, seid, lid_se):
        if len(self.lid_s0) != 0:
            lid_s0 = np.hstack([self.lid_s0, lid_s0])
            seid = np.hstack([self.seid, seid])
            lid_se = np.hstack([self.lid_se, lid_se])
        self.lid_s0 = lid_s0
        self.seid = seid
        self.lid_se = lid_se
        self.n = len(lid_s0)

    def __apply_slice__(self, card: SELOAD, i: np.ndarray) -> None:
        card.n = len(i)
        card.lid_s0 = self.lid_s0[i]
        card.seid = self.seid[i]
        card.lid_se = self.lid_se[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        return max(self.lid_s0.max(), self.seid.max(), self.lid_se.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        for lid_s0, seid, lid_se in zip(self.lid_s0, self.seid, self.lid_se):
            list_fields = ['SELOAD', lid_s0, seid, lid_se]
            bdf_file.write(print_card(list_fields))
