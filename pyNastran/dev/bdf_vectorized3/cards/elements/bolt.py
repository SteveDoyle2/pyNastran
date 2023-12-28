from __future__ import annotations
from collections import defaultdict
from itertools import zip_longest
from typing import Union, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default # , set_string8_blank_if_default
from pyNastran.bdf.cards.base_card import (
    read_ids_thru, expand_thru, _format_comment)
#from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8, print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string,
    integer_or_blank, double_or_blank,
    components_or_blank, integer_string_or_blank, integer_double_or_blank,
    integer_or_string,
    modal_components_or_blank,
)
from pyNastran.bdf.bdf_interface.assign_type_force import force_double_or_blank
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, make_idim,
    hslice_by_idim, vslice_by_idim,
    remove_unused_primary, remove_unused_duplicate,
    parse_check,
)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_float,
    array_default_int, array_default_float,
    get_print_card_size)
from pyNastran.dev.bdf_vectorized3.cards.loads.static_loads import Combination
#from pyNastran.dev.bdf_vectorized3.utils import cast_int_array
#from .static_loads import get_loads_by_load_id, get_reduced_loads
from ..bdf_sets import split_set3_ids

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class BOLTLD(Combination):
    """
    +--------+-----+------+------+----+-----+----+----+----+
    |    1   |  2  |  3   |  4   | 5  |  6  | 7  | 8  | 9  |
    +========+=====+======+======+====+=====+====+====+====+
    | BOLTLD | SID |  S   |  S1  | L1 | S2  | L2 | S3 | L3 |
    +--------+-----+------+------+----+-----+----+----+----+
    |        | S4  |  L4  | etc. |    |     |    |    |    |
    +--------+-----+------+------+----+-----+----+----+----+
    | BOLTLD | 101 | -0.5 | 1.0  | 3  | 6.2 | 4  |    |    |
    +--------+-----+------+------+----+-----+----+----+----+

    """
    _id_name = 'bolt_id'
    @property
    def bolt_id(self) -> np.ndarray:
        return self._idi
    @property
    def nbolts(self) -> np.ndarray:
        return self._n_ids
    @property
    def bolt_ids(self) -> np.ndarray:
        return self._ids_data

    @bolt_id.setter
    def bolt_id(self, bolt_id: np.ndarray) -> None:
        self.__idi = bolt_id
    @nbolts.setter
    def nbolts(self, nbolts: np.ndarray) -> None:
        self._n_ids = nbolts
    @bolt_ids.setter
    def bolt_ids(self, bolt_ids: np.ndarray) -> None:
        self._ids_data = bolt_ids

    def add(self, bolt_id: int, scale: float,
            scale_factors: list[float], bolt_ids: list[int],
            comment: str='') -> int:
        """
        Creates a BOLTLD card

        Parameters
        ----------
        sid : int
            Load set identification number. See Remarks 1. and 4. (Integer > 0)
        scale : float
            Scale factor. See Remarks 2. and 8. (Real)
        Si : list[float]
            Scale factors. See Remarks 2., 7. and 8. (Real)
        load_ids : list[int]
            Load set identification numbers of RLOAD1, RLOAD2, TLOAD1,
            TLOAD2, and ACSRCE entries. See Remarks 3. and 7. (Integer > 0)
        comment : str; default=''
            a comment for the card

        """
        return super().add(
            sid, scale, scale_factors, bolt_ids, comment=comment)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['load_id'].append(self.load_ids)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        load_id = used_dict['load_id']
        ncards_removed = remove_unused_primary(
            self, load_id, self.load_id, 'load_id')
        return ncards_removed


class BOLT(VectorizedBaseCard):
    """
    MSC
    +--------+--------+-------+-------+------+------+------+------+------+
    |   1    |   2    |   3   |   4   |  5   |  6   |   7  |  8   |  9   |
    +========+========+=======+=======+======+======+======+======+======+
    |  BOLT  | ID     | GRIDC |       |      |      |      |      |      |
    +--------+--------+-------+-------+------+------+------+------+------+
    |        | TOP    | GT1   |  GT2  |  GT3 |  GT4 |  GT5 |  GT6 |  GT7 |
    +--------+--------+-------+-------+------+------+------+------+------+
    |        | GT8    | GT9   |  etc  |      |      |      |      |      |
    +--------+--------+-------+-------+------+------+------+------+------+
    |        | BOTTOM | GB1   |  GB2  |  GB3 |  GB4 |  GB5 |  GB6 |  GB7 |
    +--------+--------+-------+-------+------+------+------+------+------+
    |        | GB8    | GB9   |  etc  |      |      |      |      |      |
    +--------+--------+-------+-------+------+------+------+------+------+
    |  BOLT  |   100  | 1025  |       |      |      |      |      |      |
    +--------+--------+-------+-------+------+------+------+------+------+
    |        |   TOP  |  101  |  102  |  103 |  104 |  105 |      |      |
    +--------+--------+-------+-------+------+------+------+------+------+
    |        | BOTTOM |   1   |   2   |   3  |   4  |   5  |      |      |
    +--------+--------+-------+-------+------+------+------+------+------+

    NX
    ---
    +-------+------+--------+-------+------+------+------+------+------+
    |   1   |  2   |   3    |   4   |  5   |  6   |   7  |  8   |  9   |
    +=======+======+========+=======+======+======+======+======+======+
    | BOLT  | BID  | ETYPE1 | EID1  | EID2 | EID3 | EID4 | EID5 | EID6 |
    +-------+------+--------+-------+------+------+------+------+------+
    |       | EID7 | THRU   | EID8  |  BY  | INC  |      |      |      |
    +-------+------+--------+-------+------+------+------+------+------+
    |       | etc  |        |       |      |      |      |      |      |
    +-------+------+--------+-------+------+------+------+------+------+
    | BOLT  | BID  | ETYPE2 | CSID  | IDIR | G1   | G2   | G3   | G4   |
    +-------+------+--------+-------+------+------+------+------+------+
    |       | G5   | THRU   |  G6   |  BY  | INC  |      |      |      |
    +-------+------+--------+-------+------+------+------+------+------+
    |       | etc  |        |       |      |      |      |      |      |
    +-------+------+--------+-------+------+------+------+------+------+
    | BOLT  | BID  | ETYPE3 | CSID  | IDIR | GP   |      |      |      |
    +-------+------+--------+-------+------+------+------+------+------+
    |       | EID1 | EID2   | EID3  | EID4 | EID5 | EID6 | EID7 | EID8 |
    +-------+------+--------+-------+------+------+------+------+------+
    |       | EID9 | THRU   | EID10 |  BY  | INC  |      |      |      |
    +-------+------+--------+-------+------+------+------+------+------+
    |       | etc  |        |       |      |      |      |      |      |
    +-------+------+--------+-------+------+------+------+------+------+
    Type=1
    BOLT BID ETYPE EID1 EID2 EID3 EID4 EID5 EID6
    EID7 THRU EID8 BY INC
    -etc-
    BOLT 4 1 11

    Type=2
    BOLT BID ETYPE CSID IDIR G1 G2 G3 G4
    G5 THRU G6 BY INC
    -etc-

    Type=3
    BOLT BID ETYPE CSID IDIR GP
    EID1 EID2 EID3 EID4 EID5 EID6 EID7 EID8
    EID9 THRU EID10 BY INC
    -etc-
    """
    _id_name = 'bolt_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.bolt_id = np.array([], dtype='int32')

    #def add(self, sid: int, desc: str, ids: list[int], comment: str='') -> int:
        #self.cards.append((sid, desc, ids, comment))
        #self.n += 1
        #return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> None:
        """
        Adds a BOLT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        bolt_id = integer(card, 1, 'sid')

        top_bottom_flag = card.field(9)
        if isinstance(top_bottom_flag, str) and top_bottom_flag.upper() in {'TOP', 'BOTTOM'}:
            top_bottom_flag = top_bottom_flag.upper()
            # msc
            #BOLT ID GRIDC
            #TOP GT1 GT2 GT3 GT4 GT5 GT6 GT7
            #GT8 GT9 etc.
            #BOTTOM GB1 GB2 GB3 GB4 GB5 GB6 GB7
            #GB8 GB9 etc.
            gridc = integer(card, 2, 'gridc')

            i = 1
            ifield = 10
            top_bottom_dict = defaultdict(list)
            while ifield < len(card):
                value = integer_string_or_blank(card, ifield, f'value{i}', default=None)
                if value is not None:
                    if value in {'TOP', 'BOTTOM'}:
                        assert value not in top_bottom_dict
                        top_bottom_flag = value
                        i = 1
                        ifield += 1
                        continue
                    top_bottom_dict[top_bottom_flag].append(value)
                    i += 1
                ifield += 1
            card = ('msc', bolt_id, gridc, top_bottom_dict, comment)
        else:
            etype = integer(card, 2, 'etype')
            if etype == 1:
                eids = read_ids_thru(card, ifield0=3, base_str='EID%d')
            elif etype == 2:
                csid = integer_or_blank(card, 3, 'csid', default=0)
                idir = integer_or_blank(card, 4, 'idir', default=0)
                assert idir in {0, 1, 2, 3}, idir
                assert len(card) <= 5, card
            elif element_type == 3:
                csid = integer_or_blank(card, 3, 'csid', default=0)
                idir = integer_or_blank(card, 4, 'idir', default=0)
                nid = integer_or_blank(card, 5, 'nid/gp', default=0)
                assert idir in {0, 1, 2, 3}, idir

                eids = read_ids_thru(card, ifield0=9, base_str='EID%d')
                #fields = card[9:]
                #eids = expand_thru_by(fields, set_fields=True, sort_fields=True, require_int=True, allow_blanks=False)
                assert len(card) <= 15, card
            else:
                raise RuntimeError(etype)
            card = (bolt_id, etype, eids, comment)

        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        card = self.cards[0]
        if card[0] == 'msc':
            self.parse_cards_msc()
        else:
            assert card[0] == 'nx', card
            self.parse_cards_nx()

    def parse_cards_nx(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype

        bolt_id = np.zeros(ncards, dtype='int32')
        num_ids = np.zeros(ncards, dtype='int32')
        all_ids = []
        for icard, card in enumerate(self.cards):
            if card[0] == 'msc':
                asdf
            else:
                assert card[0] == 'nx', card
                bolt_idi, eidsi, comment = card

            bolt_id[icard] = bolt_idi
            ids2 = split_set3_ids(eidsi)
            num_ids[icard] = len(ids2)
            all_ids.extend(ids2)
        ids = np.array(all_ids, dtype=idtype)
        self._save(bolt_id, num_ids, ids)
        self.sort()
        self.cards = []

    def _save(self, bolt_id, num_ids, element_ids):
        if len(self.bolt_id) != 0:
            bolt_id = np.hstack([self.bolt_id, bolt_id])
            num_ids = np.hstack([self.num_ids, num_ids])
            element_ids = np.hstack([self.element_ids, element_ids])
        self.bolt_id = bolt_id
        self.num_ids = num_ids
        self.element_ids = element_ids
        self.n = len(bolt_id)

    def __apply_slice__(self, bolt: BOLT, i: np.ndarray) -> None:
        assert self.num_ids.sum() == len(self.element_ids)
        bolt.n = len(i)
        bolt.bolt_id = self.bolt_id[i]

        inid = self.iprop # [i, :]
        bolt.element_ids = hslice_by_idim(i, inid, self.element_ids)

        bolt.num_ids = self.num_ids[i]
        #assert isinstance(prop.ndim, np.ndarray), prop.ndim
        #assert prop.ndim.sum() == len(prop.dims), f'prop.ndim={prop.ndim} len(prop.dims)={len(prop.dims)}'

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def iprop(self) -> np.ndarray:
        return make_idim(self.n, self.num_ids)

    @property
    def max_id(self) -> int:
        return max(self.bolt_id.max(), self.element_ids.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        asdf
        print_card, size = get_print_card_size(size, self.max_id)

        set_id = array_str(self.set_id, size=size).tolist()
        ids_ = array_str(self.property_ids, size=size).tolist()
        for sid, property_type, iprop in zip(set_id, self.property_type, self.iprop):
            iprop0, iprop1 = iprop
            ids = ids_[iprop0:iprop1]

            list_fields = ['SET4', sid, 'PROP', property_type] + ids
            bdf_file.write(print_card(list_fields))
        return




class BOLTFOR(VectorizedBaseCard):
    _id_name = 'bolt_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.bolt_id = np.array([], dtype='int32')

    #def add(self, sid: int, desc: str, ids: list[int], comment: str='') -> int:
        #self.cards.append((sid, desc, ids, comment))
        #self.n += 1
        #return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> None:
        """
        Adds a BOLT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((bolt_id, eids, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype

        bolt_id = np.zeros(ncards, dtype='int32')
        num_ids = np.zeros(ncards, dtype='int32')
        all_ids = []
        for icard, card in enumerate(self.cards):
            sid, bolt_idi, eidsi, comment = card

            bolt_id[icard] = bolt_idi
            ids2 = split_set3_ids(eidsi)
            num_ids[icard] = len(ids2)
            all_ids.extend(ids2)
        ids = np.array(all_ids, dtype=idtype)
        self._save(bolt_id, num_ids, ids)
        self.sort()
        self.cards = []

    def _save(self, bolt_id, num_ids, element_ids):
        if len(self.set_id) != 0:
            bolt_id = np.hstack([self.bolt_id, bolt_id])
            num_ids = np.hstack([self.num_ids, num_ids])
            element_ids = np.hstack([self.element_ids, element_ids])
        self.bolt_id = bolt_id
        self.num_ids = num_ids
        self.element_ids = element_ids
        self.n = len(bolt_id)

    def __apply_slice__(self, bolt: BOLT, i: np.ndarray) -> None:
        assert self.num_ids.sum() == len(self.element_ids)
        bolt.n = len(i)
        bolt.bolt_id = self.bolt_id[i]

        inid = self.iprop # [i, :]
        bolt.element_ids = hslice_by_idim(i, inid, self.element_ids)

        bolt.num_ids = self.num_ids[i]
        #assert isinstance(prop.ndim, np.ndarray), prop.ndim
        #assert prop.ndim.sum() == len(prop.dims), f'prop.ndim={prop.ndim} len(prop.dims)={len(prop.dims)}'

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def iprop(self) -> np.ndarray:
        return make_idim(self.n, self.num_ids)

    @property
    def max_id(self) -> int:
        return max(self.bolt_id.max(), self.element_ids.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        asdf
        print_card, size = get_print_card_size(size, self.max_id)

        set_id = array_str(self.set_id, size=size).tolist()
        ids_ = array_str(self.property_ids, size=size).tolist()
        for sid, property_type, iprop in zip(set_id, self.property_type, self.iprop):
            iprop0, iprop1 = iprop
            ids = ids_[iprop0:iprop1]

            list_fields = ['SET4', sid, 'PROP', property_type] + ids
            bdf_file.write(print_card(list_fields))
        return


class BOLTFRC(VectorizedBaseCard):
    """
    BOLTFRC SID TYPE D LEN
    B1 B2 B3 B4 B5 B6 B7 B8
    B9 THRU B10
    B11 B12 -etc-
    """
    _id_name = 'bolt_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.bolt_id = np.array([], dtype='int32')

    #def add(self, sid: int, desc: str, ids: list[int], comment: str='') -> int:
        #self.cards.append((sid, desc, ids, comment))
        #self.n += 1
        #return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> None:
        """
        Adds a BOLT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #BOLTFRC SID TYPE D LEN
        #B1 B2 B3 B4 B5 B6 B7 B8
        #B9 THRU B10
        #B11 B12 -etc-
        bolt_id = integer(card, 1, 'bolt_id')
        bolt_type = string(card, 2, 'bolt_type') #  DISP, STRAIN, or LOAD
        preload = double(card, 3, 'D/preload') # preload dispacement/load/strain
        length = double_or_blank(card, 4, 'bolt_id', default=np.nan)
        bolt_ids = eids = read_ids_thru(card, ifield0=9, base_str='bolt%d')
        self.cards.append((bolt_id, bolt_type, preload, length, bolt_ids, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype

        bolt_id = np.zeros(ncards, dtype='int32')
        preload_type = np.zeros(ncards, dtype='|U8')
        preload = np.zeros(ncards, dtype='float64')
        length = np.zeros(ncards, dtype='float64')
        bolt_id = np.zeros(ncards, dtype='int32')
        bolt_id = np.zeros(ncards, dtype='int32')

        num_ids = np.zeros(ncards, dtype='int32')
        all_ids = []
        for icard, card in enumerate(self.cards):
            bolt_idi, preload_typei, preloadi, lengthi, bolt_idsi, comment = card

            bolt_id[icard] = bolt_idi
            preload_type[icard] = preload_typei
            preload[icard] = preloadi
            length[icard] = lengthi

            ids2 = split_set3_ids(bolt_idsi)
            num_ids[icard] = len(ids2)
            all_ids.extend(ids2)
        bolt_ids = np.array(all_ids, dtype=idtype)
        self._save(bolt_id, preload_type, preload, length, bolt_ids, num_ids)
        self.sort()
        self.cards = []

    def _save(self, bolt_id, preload_type, preload, length, bolt_ids, num_ids):
        if len(self.bolt_id) != 0:
            bolt_id = np.hstack([self.bolt_id, bolt_id])
            preload_type = np.hstack([self.preload_type, preload_type])
            preload = np.hstack([self.preload, preload])
            length = np.hstack([self.length, length])
            bolt_ids = np.hstack([self.bolt_ids, bolt_ids])
            num_ids = np.hstack([self.num_ids, num_ids])
        self.bolt_id = bolt_id
        self.preload_type = preload_type
        self.preload = preload
        self.length = length
        self.num_ids = num_ids
        self.bolt_ids = bolt_ids
        self.n = len(bolt_id)

    def __apply_slice__(self, bolt: BOLTFRC, i: np.ndarray) -> None:
        assert self.num_ids.sum() == len(self.bolt_ids)
        bolt.n = len(i)
        bolt.bolt_id = self.bolt_id[i]
        bolt.preload_type = self.preload_type[i]
        bolt.preload = self.preload[i]
        bolt.length = self.length[i]

        ibolt = self.ibolt # [i, :]
        bolt.bolt_ids = hslice_by_idim(i, ibolt, self.bolt_ids)

        bolt.num_ids = self.num_ids[i]
        #assert isinstance(prop.ndim, np.ndarray), prop.ndim
        #assert prop.ndim.sum() == len(prop.dims), f'prop.ndim={prop.ndim} len(prop.dims)={len(prop.dims)}'

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def ibolt(self) -> np.ndarray:
        return make_idim(self.n, self.num_ids)

    @property
    def max_id(self) -> int:
        return max(self.bolt_id.max(), self.bolt_ids.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        asdf
        print_card, size = get_print_card_size(size, self.max_id)

        bolt_ids = array_str(self.bolt_id, size=size).tolist()
        ids_ = array_str(self.bolt_ids, size=size).tolist()
        for bolt_id, preload_type, preload, length, ibolt in zip(bolt_ids, self.preload_type, self.preload,
                                                                 self.length, self.ibolt):
            ibolt0, ibolt1 = ibolt
            ids = ids_[ibolt0:ibolt1]

            list_fields = ['BOLTFRC', bolt_id, preload_type, preload, length, None, None, None, None] + ids
            bdf_file.write(print_card(list_fields))
        return



class BOLTSEQ(VectorizedBaseCard):
    """
    +---------+-------+-------+-------+
    |    1    |    2  |   3   |   4   |
    +=========+=======+=======+=======+
    | BOLTSEQ |  SID  |       |       |
    +---------+-------+-------+-------+
    |         | S_NO1 | B_ID1 | NINC1 |
    +---------+-------+-------+-------+
    |         | S_NO2 | B_ID2 | NINC2 |
    +---------+-------+-------+-------+
    |         | S_NO3 | B_ID3 | NINC3 |
    +---------+-------+-------+-------+

    NX card
    """
    _id_name = 'bolt_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.bolt_id = np.array([], dtype='int32')

    def add(self,
            bolt_id: int,
            s_nos: list[int],
            b_ids: list[int],
            n_incs: Optional[list[int]]=None,
            comment: str='') -> int:
        """
        SID : int
            Bolt preload set identification number (Integer; No default)
        S_NOi : int
            Sequence order number for the BOLTLD, BOLTFOR, and BOLTFRC
            ID’s to be applied. (Integer; No default)
        B_IDi : int
            SID of BOLTLD, BOLTFOR, or BOLTFRC bulk entries defining a bolt
             preload. (Integer; No default)
        NINCi : int; default=1
            Number of increments in which to ramp up the bolt preload defined in
            BOLTLD, BOLTFOR, or BOLTFRC entries. (Integer; Default = 1)
        """
        if isinstance(s_nos, integer_types):
            s_nos = [s_nos]
        if isinstance(b_ids, integer_types):
            b_ids = [b_ids] * len(s_nos)
        assert len(s_nos) == len(b_ids)

        if isinstance(n_incs, integer_types):
            n_incs = [n_incs] * len(s_nos)
        self.cards.append((bolt_id, s_nos, b_ids, n_incs, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a BOLT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        bolt_id = integer(card, 1, 'bolt_id')
        ifield = 9
        s_nos = []
        b_ids = []
        n_incs = []
        while ifield < len(card):
            s_no = integer(card, ifield, 's_no')
            b_id = integer(card, ifield+1, 'b_id')
            n_inc = integer_or_blank(card, ifield+2, 'n_inc', default=1)
            s_nos.append(s_no)
            b_ids.append(b_id)
            n_incs.append(n_inc)
            ifield += 8
        assert len(s_nos) >= 1, s_nos
        assert len(b_ids) >= 1, b_ids
        assert len(n_incs) >= 1, n_incs
        #return BOLTSEQ(sid, s_nos, b_ids, n_incs=n_incs)
        self.cards.append((bolt_id, s_nos, b_ids, n_incs, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype

        bolt_id = np.zeros(ncards, dtype='int32')
        nsteps = np.zeros(ncards, dtype='int32')


        num_steps = np.zeros(ncards, dtype='int32')
        all_s_nos = []
        all_b_ids = []
        all_n_incs = []
        for icard, card in enumerate(self.cards):
            bolt_idi, s_nosi, b_idsi, n_incsi, comment = card

            bolt_id[icard] = bolt_idi
            nsteps[icard] = len(s_nosi)

            all_s_nos.extend(s_nosi)
            all_b_ids.extend(b_idsi)
            all_n_incs.extend(n_incsi)

        s_nos = np.array(all_s_nos, dtype=idtype)
        b_ids = np.array(all_b_ids, dtype=idtype)
        n_incs = np.array(all_n_incs, dtype=idtype)

        self._save(bolt_id, s_nos, b_ids, n_incs, nsteps)
        self.sort()
        self.cards = []

    def _save(self, bolt_id, s_nos, b_ids, n_incs, nsteps):
        if len(self.bolt_id) != 0:
            bolt_id = np.hstack([self.bolt_id, bolt_id])
            asdf
            #preload_type = np.hstack([self.preload_type, preload_type])
            #preload = np.hstack([self.preload, preload])
            #length = np.hstack([self.length, length])
            #bolt_ids = np.hstack([self.bolt_ids, bolt_ids])
            #num_ids = np.hstack([self.num_ids, num_ids])
        self.bolt_id = bolt_id
        self.s_nos = s_nos
        self.b_ids = b_ids
        self.n_incs = n_incs
        self.nsteps = nsteps
        self.n = len(bolt_id)

    def __apply_slice__(self, bolt: BOLTSEQ, i: np.ndarray) -> None:
        assert self.nsteps.sum() == len(self.s_nos)
        assert self.nsteps.sum() == len(self.b_ids)
        assert self.nsteps.sum() == len(self.n_incs)
        bolt.n = len(i)
        bolt.bolt_id = self.bolt_id[i]

        ibolt = self.ibolt # [i, :]
        bolt.s_nos = hslice_by_idim(i, ibolt, self.s_nos)
        bolt.b_ids = hslice_by_idim(i, ibolt, self.b_ids)
        bolt.n_incs = hslice_by_idim(i, ibolt, self.n_incs)

        bolt.nsteps = self.nsteps[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def ibolt(self) -> np.ndarray:
        return make_idim(self.n, self.nsteps)

    @property
    def max_id(self) -> int:
        return max(self.bolt_id.max(), self.bolt_ids.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        bolt_ids = array_str(self.bolt_id, size=size).tolist()
        for bolt_id, ibolt in zip(bolt_ids, elf.ibolt):
            ibolt0, ibolt1 = ibolt
            s_nos = self.s_nos[ibolt0:ibolt1]
            b_ids = self.b_ids[ibolt0:ibolt1]
            n_incs = self.n_incs[ibolt0:ibolt1]

            list_fields = ['BOLTSEQ', bolt_id, None, None, None, None, None, None, None]
            for s_no, b_id, n_inc in zip(s_nos, b_ids, n_incs):
                list_fields.extend([s_no, b_id, n_inc, None, None, None, None, None])
            bdf_file.write(print_card(list_fields))
        return
