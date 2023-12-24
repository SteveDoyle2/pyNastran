from __future__ import annotations
from collections import defaultdict
from itertools import zip_longest
from typing import Union, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default # , set_string8_blank_if_default
from pyNastran.bdf.cards.base_card import expand_thru_by # expand_thru # BaseCard, expand_thru_by #  _node_ids,
#from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8, print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double,
    integer_or_blank, double_or_blank,
    components_or_blank, integer_string_or_blank, integer_double_or_blank,
    integer_or_string,
    modal_components_or_blank,
)
from pyNastran.bdf.bdf_interface.assign_type_force import force_double_or_blank
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.bdf.cards.loads.dloads import (
    fix_loadtype_tload1, fix_loadtype_tload2,
    fix_loadtype_rload1, # fix_loadtype_rload2,
)

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, make_idim,
    hslice_by_idim, vslice_by_idim,
    remove_unused_primary, remove_unused_duplicate,
    parse_load_check,
)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_float,
    array_default_int, array_default_float,
    get_print_card_size)
from pyNastran.dev.bdf_vectorized3.cards.loads.static_loads import LoadCombination
#from pyNastran.dev.bdf_vectorized3.utils import cast_int_array
#from .static_loads import get_loads_by_load_id, get_reduced_loads

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class DLOAD(LoadCombination):
    """
    +-------+-----+------+------+----+-----+----+----+----+
    |   1   |  2  |  3   |  4   | 5  |  6  | 7  | 8  | 9  |
    +=======+=====+======+======+====+=====+====+====+====+
    | DLOAD | SID |  S   |  S1  | L1 | S2  | L2 | S3 | L3 |
    +-------+-----+------+------+----+-----+----+----+----+
    |       | S4  |  L4  | etc. |    |     |    |    |    |
    +-------+-----+------+------+----+-----+----+----+----+
    | DLOAD | 101 | -0.5 | 1.0  | 3  | 6.2 | 4  |    |    |
    +-------+-----+------+------+----+-----+----+----+----+

    """
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')
        self.nloads = np.array([], dtype='int32')
        self.load_ids = np.array([], dtype='int32')
        self.scale_factors = np.array([], dtype='float64')

    def add(self, sid: int, scale: float,
            scale_factors: list[float], load_ids: list[int],
            comment: str='') -> int:
        """
        Creates a DLOAD card

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
            sid, scale, scale_factors, load_ids, comment=comment)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['load_id'].append(self.load_ids)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        load_id = used_dict['load_id']
        ncards_removed = remove_unused_primary(
            self, load_id, self.load_id, 'load_id')
        return ncards_removed

    def get_loads_by_load_id(self) -> dict[int, Loads]:
        model = self.model
        """"""
        #uload_ids = np.unique(self.load_ids)
        loads_by_load_id = defaultdict(list)

        for load in model.dynamic_load_cards:
            if load.type in {'DLOAD'}:
                continue
            #print(load)
            uload_idsi = np.unique(load.load_id)
            for uload_id in uload_idsi:
                i = np.where(uload_id == load.load_id)[0]
                if len(i) == 0:
                    continue
                loadi = load.slice_card_by_index(i)

                loads_by_load_id[uload_id].append(loadi)
        return dict(loads_by_load_id)

    def get_reduced_loads(self,
                          remove_missing_loads: bool=False,
                          filter_zero_scale_factors: bool=False,
                          stop_on_failure: bool=True) -> dict[int, Loads]:
        """
        Parameters
        ----------
        resolve_load_card : bool; default=False
            ???
        remove_missing_loads: bool; default=False
            LOAD cards can reference loads (e.g., GRAV) that don't exist
            Nastran sometimes ignores these loads leading to potentially incorrect results
        filter_zero_scale_factors: bool; default=False
            remove loads that are 0.0
        """
        reduced_loads = {}
        if self.n == 0:
            return reduced_loads

        stop_on_failure = True
        loads_by_load_id = self.get_loads_by_load_id()
        log = self.model.log
        for sid, global_scale, iload in zip(self.load_id, self.scale_factors, self.iload):
            reduced_loadsi = []
            #print(self.load_id, self.scale_factors, self.iload, iload)
            iload0, iload1 = iload
            if global_scale == 0. and filter_zero_scale_factors:
                continue
            scale_factors = global_scale * self.scale_factors[iload0:iload1]
            load_ids = self.load_ids[iload0:iload1]
            for (scale_factor, load_id) in zip_longest(scale_factors, load_ids):
                if scale_factor == 0. and filter_zero_scale_factors:
                    continue
                loads_found = loads_by_load_id[load_id]
                if len(loads_found) == 0:
                    msg = f'No referenced loads found for load_id={load_id} on DLOAD load_id={sid}'
                    log.error(msg)
                    if stop_on_failure:
                        raise RuntimeError(msg)
                reduced_loadsi.append((scale_factor, loads_found))
            reduced_loads[sid] = reduced_loadsi
        return reduced_loads


class DAREA(VectorizedBaseCard):
    """
    Defines scale (area) factors for static and dynamic loads. In dynamic
    analysis, DAREA is used in conjunction with ACSRCE, RLOADi and TLOADi
    entries.

    RLOAD1 -> DAREA by SID

    +-------+-----+----+----+-----+----+----+------+
    |   1   |  2  | 3  |  4 |  5  | 6  |  7 |  8   |
    +=======+=====+====+====+=====+====+====+======+
    | DAREA | SID | P1 | C1 |  A1 | P2 | C2 |  A2  |
    +-------+-----+----+----+-----+----+----+------+
    | DAREA |  3  | 6  | 2  | 8.2 | 15 | 1  | 10.1 |
    +-------+-----+----+----+-----+----+----+------+
    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')

    def slice_card_by_index(self, i: np.ndarray) -> DAREA:
        load = DAREA(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: DAREA, i: np.ndarray) -> None:
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.node_id = self.node_id[i]
        load.component = self.component[i]
        load.scale = self.scale[i]

    def add(self, sid: int, nid: int, component: str, scale: float,
            comment: str='') -> int:
        """
        Creates a DAREA card

        Parameters
        ----------
        sid : int
            darea id
        nid : int
            GRID, EPOINT, SPOINT id
        component : str
            Component number. (0-6; 0-EPOINT/SPOINT; 1-6 GRID)
        scale : float
            Scale (area) factor

        """
        self.cards.append((sid, nid, component, scale, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        # sid : int
        #     darea id
        # nodes : list[int]
        #     GRID, EPOINT, SPOINT id
        # components : list[int]
        #     Component number. (0-6; 0-EPOINT/SPOINT; 1-6 GRID)
        # scales : list[float]
        #     Scale (area) factor
        sid = integer(card, 1, 'sid')
        nid = integer(card, 2, 'p')
        component = int(components_or_blank(card, 3, 'c', default=0))
        scale = double(card, 4, 'scale')
        self.cards.append((sid, nid, component, scale, comment))
        self.n += 1

        #    0       1     2     3    4       5
        #['DAREA', '12', '102', '1', '1.0', '102', '4', '1.0']
        if card.field(5):  # node id
            comment = ''
            nid = integer(card, 5, 'p')
            component = int(components_or_blank(card, 6, 'c', default=0))
            scale = double(card, 7, 'scale')
            self.cards.append((sid, nid, component, scale, comment))
            self.n += 1
        assert len(card) <= 8, f'len(DAREA card) = {len(card):d}\ncard={card}'
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype

        #: Set identification number
        load_id = np.zeros(ncards, dtype='int32')


        #: Identification number of DAREA or SPCD entry set or a thermal load
        #: set (in heat transfer analysis) that defines {A}. (Integer > 0)
        node_id = np.zeros(ncards, dtype=idtype)
        component = np.zeros(ncards, dtype='int32')
        scale = np.zeros(ncards, dtype='float64')

        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (sid, nid, componenti, scalei, comment) = card
            load_id[icard] = sid
            node_id[icard] = nid
            component[icard] = componenti
            scale[icard] = scalei
        self._save(load_id, node_id, component, scale)
        assert len(self.load_id) == self.n
        self.cards = []

    def _save(self, load_id, node_id, component, scale):
        self.load_id = load_id
        self.node_id = node_id
        self.component = component
        self.scale = scale

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.node_id)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.node_id
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(), self.node_id.max())

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        #array_str, array_default_int
        load_ids = array_str(self.load_id, size=size)
        node_ids = array_str(self.node_id, size=size)
        components = array_default_int(self.component, default=0, size=size)
        for sid, nid, comp, scale in zip_longest(load_ids, node_ids, components, self.scale):
            list_fields = ['DAREA', sid, nid, comp, scale]
            bdf_file.write(print_card(list_fields))
        return


class TLOAD1(VectorizedBaseCard):
    r"""
    Transient Response Dynamic Excitation, Form 1

    Defines a time-dependent dynamic load or enforced motion of the form:

    .. math::
      \left\{ P(t) \right\} = \left\{ A \right\} \cdot F(t-\tau)

    for use in transient response analysis.

    MSC 20005.2
    +--------+-----+----------+-------+------+-----+-----+-----+
    |    1   |  2  |     3    |   4   |   5  |  6  |  7  |  8  |
    +========+=====+==========+=======+======+=====+=====+=====+
    | TLOAD1 | SID | EXCITEID | DELAY | TYPE | TID | US0 | VS0 |
    +--------+-----+----------+-------+------+-----+-----+-----+

    NX 11
    +--------+-----+----------+-------+------+-----+
    |    1   |  2  |     3    |   4   |   5  |  6  |
    +========+=====+==========+=======+======+=====+
    | TLOAD1 | SID | EXCITEID | DELAY | TYPE | TID |
    +--------+-----+----------+-------+------+-----+
    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        #: Set identification number
        self.load_id = np.array([], dtype='int32')

        #: Identification number of DAREA or SPCD entry set or a thermal load
        #: set (in heat transfer analysis) that defines {A}. (Integer > 0)
        self.excite_id = np.array([], dtype='int32')

        #: If it is a non-zero integer, it represents the
        #: identification number of DELAY Bulk Data entry that defines .
        #: If it is real, then it directly defines the value of that will
        #: be used for all degrees-of-freedom that are excited by this
        #: dynamic load entry.  See also Remark 9. (Integer >= 0,
        #: real or blank)
        self.delay_int = np.array([], dtype='int32')
        self.delay_float = np.array([], dtype='float64')

        #: Defines the type of the dynamic excitation. (LOAD,DISP, VELO, ACCE)
        self.load_type = np.array([], dtype='|U4')

        #: Identification number of TABLEDi entry that gives F(t). (Integer > 0)
        self.tabled_id = np.array([], dtype='int32')

        #: Factor for initial displacements of the enforced degrees-of-freedom.
        #: (Real; Default = 0.0)
        self.us0 = np.array([], dtype='float64')

        #: Factor for initial velocities of the enforced degrees-of-freedom.
        #: (Real; Default = 0.0)
        self.vs0 = np.array([], dtype='float64')

    def slice_card_by_index(self, i: np.ndarray) -> TLOAD1:
        load = TLOAD1(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: TLOAD1, i: np.ndarray) -> None:
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.excite_id = self.excite_id[i]

        load.delay_int = self.delay_int[i]
        load.delay_float = self.delay_float[i]
        load.load_type = self.load_type[i]
        load.tabled_id = self.tabled_id[i]
        load.us0 = self.us0[i]
        load.vs0 = self.vs0[i]

    def add(self, sid: int, excite_id: int, tabled_id: int, delay: int=0,
            load_type: str='LOAD', us0: float=0.0, vs0: float=0.0,
            comment: str='') -> int:
        """
        Creates a TLOAD1 card, which defines a load based on a table

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        tabled_id : int
            TABLEDi id that defines F(t) for all degrees of freedom in
            EXCITEID entry
            float : MSC not supported
        delay : int/float; default=0
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        load_type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        us0 : float; default=0.
            Factor for initial displacements of the enforced degrees-of-freedom
            MSC only
        vs0 : float; default=0.
            Factor for initial velocities of the enforced degrees-of-freedom
            MSC only
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, excite_id, delay, load_type, tabled_id, us0, vs0, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', default=0.0)
        load_type = integer_string_or_blank(card, 4, 'Type', default='LOAD')
        tabled_id = integer(card, 5, 'tid')
        us0 = double_or_blank(card, 6, 'us0', default=0.0)
        vs0 = double_or_blank(card, 7, 'vs0', default=0.0)

        assert len(card) <= 8, f'len(TLOAD1 card) = {len(card):d}\ncard={card}'
        #return TLOAD1(sid, excite_id, tid, delay=delay, Type=Type, us0=us0, vs0=vs0, comment=comment)
        load_type_str = fix_loadtype_tload1(load_type)
        self.cards.append((sid, excite_id, delay, load_type_str, tabled_id, us0, vs0, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)

        #: Set identification number
        load_id = np.zeros(ncards, dtype='int32')


        #: Identification number of DAREA or SPCD entry set or a thermal load
        #: set (in heat transfer analysis) that defines {A}. (Integer > 0)
        excite_id = np.zeros(ncards, dtype='int32')

        #: If it is a non-zero integer, it represents the
        #: identification number of DELAY Bulk Data entry that defines .
        #: If it is real, then it directly defines the value of that will
        #: be used for all degrees-of-freedom that are excited by this
        #: dynamic load entry.  See also Remark 9. (Integer >= 0,
        #: real or blank)
        delay_int = np.zeros(ncards, dtype='int32')
        delay_float = np.full(ncards, np.nan, dtype='float64')

        #: Defines the type of the dynamic excitation. (LOAD,DISP, VELO, ACCE)
        load_type = np.zeros(ncards, dtype='|U4')

        #: Identification number of TABLEDi entry that gives F(t). (Integer > 0)
        tabled_id = np.zeros(ncards, dtype='int32')

        #: Factor for initial displacements of the enforced degrees-of-freedom.
        #: (Real; Default = 0.0)
        us0 = np.zeros(ncards, dtype='float64')

        #: Factor for initial velocities of the enforced degrees-of-freedom.
        #: (Real; Default = 0.0)
        vs0 = np.zeros(ncards, dtype='float64')

        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (sid, excite_idi, delay, load_type_str, tabled_idi, us0i, vs0i, comment) = card
            load_id[icard] = sid
            excite_id[icard] = excite_idi
            _set_int_float(icard, delay_int, delay_float, delay)

            load_type[icard] = load_type_str
            tabled_id[icard] = tabled_idi
            us0[icard] = us0i
            vs0[icard] = vs0i
        self._save(load_id, excite_id, delay_int, delay_float, load_type, tabled_id, us0, vs0)
        assert len(self.load_id) == self.n
        self.cards = []

    def _save(self, load_id, excite_id, delay_int, delay_float, load_type, tabled_id, us0, vs0):
        if len(self.load_id) == 0:
            load_id = np.hstack([self.load_id, load_id])
            excite_id = np.hstack([self.excite_id, excite_id])
            delay_int = np.hstack([self.delay_int, delay_int])
            delay_float = np.hstack([self.delay_float, delay_float])
            load_type = np.hstack([self.load_type, load_type])
            tabled_id = np.hstack([self.tabled_id, tabled_id])
            us0 = np.hstack([self.us0, us0])
            vs0 = np.hstack([self.vs0, vs0])
        self.load_id = load_id
        self.excite_id = excite_id
        self.delay_int = delay_int
        self.delay_float = delay_float
        self.load_type = load_type
        self.tabled_id = tabled_id
        self.us0 = us0
        self.vs0 = vs0

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        # TODO: self.excite_id
        delay = self.delay_int[self.delay_int > 0]
        used_dict['delay_id'].append(delay)
        table = self.tabled_id
        tabled = table[table > 0]
        used_dict['tabled_id'].append(tabled)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        load_id = used_dict['load_id']
        ncards_removed = remove_unused_duplicate(
            self, load_id, self.load_id, 'load_id')
        return ncards_removed

    @property
    def max_id(self):
        return max(self.load_id.max(),
                   self.excite_id.max(),
                   self.delay_int.max(),
                   self.tabled_id.max())

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if self.n == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)

        #array_str, array_default_int
        load_ids = array_str(self.load_id, size=size)
        excite_ids = array_default_int(self.excite_id, default=0, size=size)
        delays = array_int_float(self.delay_int, self.delay_float, size=size, is_double=False)
        us0s = array_default_float(self.us0, default=0.0, size=size)
        vs0s = array_default_float(self.vs0, default=0.0, size=size)

        for sid, excite_id, delay, load_type, tabled_id, us0, \
            vs0 in zip_longest(load_ids, excite_ids,
                               delays, self.load_type,
                               self.tabled_id, us0s, vs0s):
            list_fields = ['TLOAD1', sid, excite_id, delay, load_type,
                           tabled_id, us0, vs0]
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        all_delay = np.unique(model.delay.delay_id)
        udelay = np.unique(self.delay_int)

        geom_check(
            self,
            missing,
            delay=(all_delay, udelay),
        )


class TLOAD2(VectorizedBaseCard):
    r"""
    Transient Response Dynamic Excitation, Form 1

    Defines a time-dependent dynamic load or enforced motion of the form:

    .. math::
    \left\{ P(t) \right\} = \left\{ A \right\} e^(C*t) cos(2 \pi f t + \phi)

      P(t) = 0                                            (t<T1+tau or t >  T2+tau)
      P(t) = {A} * t^b * e^(C*t) * cos(2*pi*f*t + phase)  (T1+tau <=   t <= T2+tau)

    for use in transient response analysis.

    MSC 2016.1
    +--------+-----+----------+-------+------+-----+-----+--------+---------+
    |    1   |  2  |     3    |   4   |   5  |  6  |  7  |    8   |    9    |
    +========+=====+==========+=======+======+=====+=====+========+=========+
    | TLOAD2 | SID | EXCITEID | DELAY | TYPE | T1  | T2  |  FREQ  |  PHASE  |
    +--------+-----+----------+-------+------+-----+-----+--------+---------+
    |        |  C  |     B    |  US0  |  VS0 |     |     |        |         |
    +--------+-----+----------+-------+------+-----+-----+--------+---------+

    NX 11
    +--------+-----+----------+-------+------+-----+-----+--------+---------+
    |    1   |  2  |     3    |   4   |   5  |  6  |  7  |    8   |    9    |
    +========+=====+==========+=======+======+=====+=====+========+=========+
    | TLOAD2 | SID | EXCITEID | DELAY | TYPE | T1  | T2  |  FREQ  |  PHASE  |
    +--------+-----+----------+-------+------+-----+-----+--------+---------+
    |        |  C  |     B    |       |      |     |     |        |         |
    +--------+-----+----------+-------+------+-----+-----+--------+---------+

    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> TLOAD2:
        #load = TLOAD2(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def add(self, sid: int, excite_id: int, delay: int=0,
            load_type: str='LOAD',
            T1: float=0., T2: Optional[float]=None,
            frequency: float=0., phase: float=0.,
            c: float=0., b: float=0.,
            us0: float=0., vs0: float=0.,
            comment: str='') -> int:
        """
        Creates a TLOAD2 card, which defines a exponential time load

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        delay : int/float; default=None
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        load_type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        T1 : float; default=0.
            time constant (t1 > 0.0)
            times below this are ignored
        T2 : float; default=None
            time constant (t2 > t1)
            times above this are ignored
        frequency : float; default=0.
            Frequency in cycles per unit time.
        phase : float; default=0.
            Phase angle in degrees.
        c : float; default=0.
            Exponential coefficient.
        b : float; default=0.
            Growth coefficient.
        us0 : float; default=0.
            Factor for initial displacements of the enforced degrees-of-freedom
            MSC only
        vs0 : float; default=0.
            Factor for initial velocities of the enforced degrees-of-freedom
            MSC only
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, excite_id, delay, load_type, [T1, T2],
                           frequency, phase, b, c, us0, vs0, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', default=0)
        load_type = integer_string_or_blank(card, 4, 'Type', default='LOAD')

        T1 = double_or_blank(card, 5, 'T1', default=0.0)
        T2 = double_or_blank(card, 6, 'T2', default=T1)
        frequency = double_or_blank(card, 7, 'frequency', default=0.)
        phase = double_or_blank(card, 8, 'phase', default=0.)
        c = double_or_blank(card, 9, 'c', default=0.)
        b = double_or_blank(card, 10, 'b', default=0.)
        us0 = double_or_blank(card, 11, 'us0', default=0.)
        vs0 = double_or_blank(card, 12, 'vs0', default=0.)
        assert len(card) <= 13, f'len(TLOAD2 card) = {len(card):d}\ncard={card}'
        load_type_str = fix_loadtype_tload2(load_type)
        self.cards.append((sid, excite_id, delay, load_type_str, [T1, T2],
                           frequency, phase, b, c, us0, vs0, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)

        #: Set identification number
        load_id = np.zeros(ncards, dtype='int32')

        #: Identification number of DAREA or SPCD entry set or a thermal load
        #: set (in heat transfer analysis) that defines {A}. (Integer > 0)
        excite_id = np.zeros(ncards, dtype='int32')

        #: If it is a non-zero integer, it represents the
        #: identification number of DELAY Bulk Data entry that defines .
        #: If it is real, then it directly defines the value of that will
        #: be used for all degrees-of-freedom that are excited by this
        #: dynamic load entry.  See also Remark 9. (Integer >= 0,
        #: real or blank)
        delay_int = np.zeros(ncards, dtype='int32')
        delay_float = np.zeros(ncards, dtype='float64')

        #: Defines the type of the dynamic excitation. (LOAD,DISP, VELO, ACCE)
        load_type = np.zeros(ncards, dtype='|U4')

        #: Identification number of TABLEDi entry that gives F(t). (Integer > 0)
        #tabled_id = np.zeros(ncards, dtype='int32')

        #: T1: Time constant. (Real >= 0.0)
        #: T2: Time constant. (Real; T2 > T1)
        T = np.zeros((ncards, 2), dtype='float64')

        #: Frequency in cycles per unit time. (Real >= 0.0; Default=0.0)
        frequency = np.zeros(ncards, dtype='float64')

        #: Phase angle in degrees. (Real; Default=0.0)
        phase = np.zeros(ncards, dtype='float64')

        #: Growth coefficient. (Real; Default=0.0)
        b = np.zeros(ncards, dtype='float64')

        #: Exponential coefficient. (Real; Default=0.0)
        c = np.zeros(ncards, dtype='float64')

        #: Factor for initial displacements of the enforced degrees-of-freedom.
        #: (Real; Default=0.0)
        us0 = np.zeros(ncards, dtype='float64')

        #: Factor for initial velocities of the enforced degrees-of-freedom.
        #: (Real; Default = 0.0)
        vs0 = np.zeros(ncards, dtype='float64')

        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (sid, excite_idi, delayi, load_type_str, time_constant,
             frequencyi, phasei, bi, ci, us0i, vs0i, comment) = card
            load_id[icard] = sid
            excite_id[icard] = excite_idi
            _set_int_float(icard, delay_int, delay_float, delayi)

            load_type[icard] = load_type_str

            #tabled_id[icard] = tabled_idi
            T[icard] = time_constant
            frequency[icard] = frequencyi
            phase[icard] = phasei
            b[icard] = bi
            c[icard] = ci
            us0[icard] = us0i
            vs0[icard] = vs0i
        self._save(load_id, excite_id, load_type, T,
                   delay_int, delay_float,
                   frequency, phase, b, c, us0, vs0)
        assert len(self.load_id) == self.n
        self.cards = []

    def _save(self, load_id, excite_id, load_type, T,
              delay_int, delay_float, frequency, phase, b, c, us0, vs0):
        self.load_id = load_id
        self.excite_id = excite_id
        self.load_type = load_type
        self.T = T
        self.frequency = frequency
        self.delay_int = delay_int
        self.delay_float = delay_float
        self.phase = phase
        self.b = b
        self.c = c
        self.us0 = us0
        self.vs0 = vs0

    def __apply_slice__(self, load: TLOAD2, i: np.ndarray) -> None:
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.excite_id = self.excite_id[i]

        load.T = self.T[i]
        load.frequency = self.frequency[i]
        load.delay_int = self.delay_int[i]
        load.delay_float = self.delay_float[i]

        load.phase = self.phase[i]
        load.b = self.b[i]
        load.c = self.c[i]
        load.us0 = self.us0[i]
        load.vs0 = self.vs0[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        # TODO: self.excite_id
        delay = self.delay_int[self.delay_int > 0]
        used_dict['delay_id'].append(delay)

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(),
                   self.excite_id.max(),
                   self.delay_int.max())

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        #array_str, array_default_int
        load_ids = array_str(self.load_id, size=size)
        excite_ids = array_default_int(self.excite_id, default=0, size=size)
        delays = array_int_float(self.delay_int, self.delay_float, size=size, is_double=False)

        for sid, excite_id, delay, load_type, T, \
            frequency, phase, b, c, \
            us0, vs0 in zip_longest(load_ids, excite_ids, delays,
                                    self.load_type, self.T,
                                    self.frequency, self.phase, self.b, self.c,
                                    self.us0, self.vs0):
            T1, T2 = T
            frequency = set_blank_if_default(frequency, 0.0)
            phase = set_blank_if_default(phase, 0.0)
            c = set_blank_if_default(c, 0.0)
            b = set_blank_if_default(b, 0.0)

            us0 = set_blank_if_default(us0, 0.0)
            vs0 = set_blank_if_default(vs0, 0.0)
            list_fields = ['TLOAD2', sid, excite_id, delay, load_type,
                           T1, T2, frequency, phase, c, b, us0, vs0]
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        all_delay = np.unique(model.delay.delay_id)
        udelay = np.unique(self.delay_int)

        geom_check(
            self,
            missing,
            delay=(all_delay, udelay),
        )


class RLOAD1(VectorizedBaseCard):
    r"""
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.

    .. math::
      \left\{ P(f)  \right\}  = \left\{A\right\} [ C(f)+iD(f)]
         e^{  i \left\{\theta - 2 \pi f \tau \right\} }

    +--------+-----+----------+-------+--------+----+----+------+
    |   1    |  2  |     3    |   4   |   5    |  6 |  7 |   8  |
    +========+=====+==========+=======+========+====+====+======+
    | RLOAD1 | SID | EXCITEID | DELAY | DPHASE | TC | TD | TYPE |
    +--------+-----+----------+-------+--------+----+----+------+
    | RLOAD1 |  5  |    3     |       |        | 1  |    |      |
    +--------+-----+----------+-------+--------+----+----+------+

    NX allows DELAY and DPHASE to be floats
    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> RLOAD1:
        #load = RLOAD1(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def __apply_slice__(self, load: RLOAD1, i: np.ndarray) -> None:
        #self.model.log.info(self.dphase_int)
        #self.model.log.info(i)
        assert len(i) > 0
        assert len(self.dphase_int) > 0
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.excite_id = self.excite_id[i]
        load.load_type = self.load_type[i]

        # ints
        load.dphase_int = self.dphase_int[i]
        load.delay_int = self.delay_int[i]
        load.tabled_c_int = self.tabled_c_int[i]
        load.tabled_d_int = self.tabled_d_int[i]

        # floats
        load.dphase_float = self.dphase_float[i]
        load.delay_float = self.delay_float[i]
        load.tabled_c_float = self.tabled_c_float[i]
        load.tabled_d_float = self.tabled_d_float[i]

    def add(self, sid: int, excite_id: int,
            delay: Union[int, float]=0,
            dphase: Union[int, float]=0,
            tc: Union[int, float]=0,
            td: Union[int, float]=0,
            load_type='LOAD', comment: str='') -> int:
        """
        Creates an RLOAD1 card, which defines a frequency-dependent load
        based on TABLEDs.

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        delay : int/float; default=None
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        dphase : int/float; default=None
            the dphase; if it's 0/blank there is no phase lag
            float : delay in units of time
            int : delay id
        tc : int/float; default=0
            TABLEDi id that defines C(f) for all degrees of freedom in
            EXCITEID entry
        td : int/float; default=0
            TABLEDi id that defines D(f) for all degrees of freedom in
            EXCITEID entry
        load_type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, excite_id, delay, dphase, tc, td, load_type, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', default=0)
        dphase = integer_double_or_blank(card, 4, 'dphase', default=0)
        tc = integer_double_or_blank(card, 5, 'tc', default=0)
        td = integer_double_or_blank(card, 6, 'td', default=0)

        load_type = integer_string_or_blank(card, 7, 'Type', default='LOAD')
        assert len(card) <= 8, f'len(RLOAD1 card) = {len(card):d}\ncard={card}'
        load_type_str = fix_loadtype_rload1(load_type)

        self.cards.append((sid, excite_id, delay, dphase, tc, td, load_type_str, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        #: Set identification number
        load_id = np.zeros(ncards, dtype='int32')

        #: Identification number of DAREA or SPCD entry set or a thermal load
        #: set (in heat transfer analysis) that defines {A}. (Integer > 0)
        excite_id = np.zeros(ncards, dtype='int32')

        #: If it is a non-zero integer, it represents the
        #: identification number of DELAY Bulk Data entry that defines .
        #: If it is real, then it directly defines the value of that will
        #: be used for all degrees-of-freedom that are excited by this
        #: dynamic load entry.  See also Remark 9. (Integer >= 0,
        #: real or blank)
        delay_int = np.zeros(ncards, dtype='int32')
        delay_float = np.full(ncards, np.nan, dtype='float64')

        #: Defines the type of the dynamic excitation. (LOAD,DISP, VELO, ACCE)
        load_type = np.zeros(ncards, dtype='|U4')

        dphase_int = np.zeros(ncards, dtype='int32')
        dphase_float = np.full(ncards, np.nan, dtype='float64')

        # tc : int/float; default=0
        #     TABLEDi id that defines C(f) for all degrees of freedom in EXCITEID entry
        # td : int/float; default=0
        #     TABLEDi id that defines D(f) for all degrees of freedom in EXCITEID entry
        tabled_c_int = np.zeros(ncards, dtype='int32')
        tabled_c_float = np.full(ncards, np.nan, dtype='float64')

        tabled_d_int = np.zeros(ncards, dtype='int32')
        tabled_d_float = np.full(ncards, np.nan, dtype='float64')

        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (sid, excite_idi, delay, dphase, tc, td, load_type_str, comment) = card
            load_id[icard] = sid
            excite_id[icard] = excite_idi
            _set_int_float(icard, delay_int, delay_float, delay)
            _set_int_float(icard, dphase_int, dphase_float, dphase)
            _set_int_float(icard, tabled_c_int, tabled_c_float, tc)
            _set_int_float(icard, tabled_d_int, tabled_d_float, td)

            load_type[icard] = load_type_str
        self._save(load_id, excite_id, load_type,
                   delay_int, delay_float,
                   dphase_int, dphase_float,
                   tabled_c_int, tabled_c_float,
                   tabled_d_int, tabled_d_float)
        assert len(self.load_id) == self.n
        self.cards = []

    def _save(self, load_id, excite_id, load_type,
              delay_int, delay_float,
              dphase_int, dphase_float,
              tabled_c_int, tabled_c_float,
              tabled_d_int, tabled_d_float) -> None:
        assert len(self.load_id) == 0
        self.load_id = load_id
        self.excite_id = excite_id
        self.load_type = load_type
        self.delay_int = delay_int
        self.delay_float = delay_float
        self.dphase_int = dphase_int
        self.dphase_float = dphase_float
        self.tabled_c_int = tabled_c_int
        self.tabled_c_float = tabled_c_float
        self.tabled_d_int = tabled_d_int
        self.tabled_d_float = tabled_d_float

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        # TODO: excite_id
        delay = self.delay_int[self.delay_int > 0]
        dphase = self.dphase_int[self.dphase_int > 0]
        used_dict['delay_id'].append(delay)
        used_dict['dphase_id'].append(dphase)
        for table in (self.tabled_c_int, self.tabled_d_int):
            tabled = table[table > 0]
            used_dict['tabled_id'].append(tabled)

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(), self.excite_id.max(),
                   self.delay_int.max(), self.dphase_int.max(),
                   self.tabled_c_int.max(), self.tabled_d_int.max())

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        #array_str, array_default_int
        load_ids = array_str(self.load_id, size=size)
        excite_ids = array_default_int(self.excite_id, default=0, size=size)
        delays = array_int_float(self.delay_int, self.delay_float, size=size, is_double=False)
        dphases = array_int_float(self.dphase_int, self.dphase_float, size=size, is_double=False)
        tcs = array_int_float(self.tabled_c_int, self.tabled_c_float, size=size, is_double=False)
        tds = array_int_float(self.tabled_d_int, self.tabled_d_float, size=size, is_double=False)

        for sid, excite_id, load_type, delay, dphase, \
            tc, td in zip_longest(load_ids, excite_ids, self.load_type,
                                  delays, dphases,
                                  tcs, tds):

            #dphase = _write_int_float(dphase_int, dphase_float)
            #delay = _write_int_float(delay_int, delay_float)
            #tc = _write_int_float(tc_int, tc_float)
            #td = _write_int_float(td_int, td_float)
            list_fields = ['RLOAD1', sid, excite_id, delay, dphase,
                           tc, td, load_type]
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        udelay = np.unique(self.delay_int)
        udphase = np.unique(self.dphase_int)
        all_delay = np.unique(model.delay.delay_id)
        all_dphase = np.unique(model.dphase.dphase_id)

        geom_check(
            self,
            missing,
            delay=(all_delay, udelay),
            dphase=(all_dphase, udphase),
        )


class RLOAD2(VectorizedBaseCard):
    r"""
    Defines a frequency-dependent dynamic load of the form
    for use in frequency response problems.

    .. math:: \left\{ P(f)  \right\}  = \left\{A\right\} * B(f)
        e^{  i \left\{ \phi(f) + \theta - 2 \pi f \tau \right\} }

    +--------+-----+----------+-------+--------+----+----+------+
    |   1    |  2  |     3    |   4   |    5   |  6 |  7 |  8   |
    +========+=====+==========+=======+========+====+====+======+
    | RLOAD2 | SID | EXCITEID | DELAY | DPHASE | TB | TP | TYPE |
    +--------+-----+----------+-------+--------+----+----+------+
    | RLOAD2 |  5  |    3     |       |        | 1  |    |      |
    +--------+-----+----------+-------+--------+----+----+------+

    NX allows DELAY and DPHASE to be floats
    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> RLOAD2:
        #load = RLOAD2(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def __apply_slice__(self, load: RLOAD2, i: np.ndarray) -> None:
        #self.model.log.info(self.dphase_int)
        #self.model.log.info(i)
        assert len(i) > 0
        assert len(self.dphase_int) > 0
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.excite_id = self.excite_id[i]
        load.load_type = self.load_type[i]

        # ints
        load.dphase_int = self.dphase_int[i]
        load.delay_int = self.delay_int[i]
        load.tabled_b_int = self.tabled_b_int[i]
        load.tabled_phi_int = self.tabled_phi_int[i]

        # floats
        load.dphase_float = self.dphase_float[i]
        load.delay_float = self.delay_float[i]
        load.tabled_b_float = self.tabled_b_float[i]
        load.tabled_phi_float = self.tabled_phi_float[i]

    def add(self, sid: int, excite_id: int,
            delay: Union[int, float]=0,
            dphase: Union[int, float]=0,
            tb: Union[int, float]=0,
            tphi: Union[int, float]=0,
            load_type: str='LOAD', comment: str='') -> int:
        """
        Creates an RLOAD2 card, which defines a frequency-dependent load
        based on TABLEDs.

        Parameters
        ----------
        sid : int
            load id
        excite_id : int
            node id where the load is applied
        delay : int/float; default=None
            the delay; if it's 0/blank there is no delay
            float : delay in units of time
            int : delay id
        dphase : int/float; default=None
            the dphase; if it's 0/blank there is no phase lag
            float : delay in units of time
            int : delay id
        tb : int/float; default=0
            TABLEDi id that defines B(f) for all degrees of freedom in
            EXCITEID entry
        tp : int/float; default=0
            TABLEDi id that defines phi(f) for all degrees of freedom in
            EXCITEID entry
        load_type : int/str; default='LOAD'
            the type of load
            0/LOAD
            1/DISP
            2/VELO
            3/ACCE
            4, 5, 6, 7, 12, 13 - MSC only
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, excite_id, delay, dphase, tb, tphi, load_type))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment:str='') -> None:
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        delay = integer_double_or_blank(card, 3, 'delay', default=0)
        dphase = integer_double_or_blank(card, 4, 'dphase', default=0)
        tb = integer_double_or_blank(card, 5, 'tb', default=0)
        tphi = integer_double_or_blank(card, 6, 'tp', default=0)
        load_type = integer_string_or_blank(card, 7, 'Type', default='LOAD')
        load_type_str = fix_loadtype_rload1(load_type)

        assert len(card) <= 8, f'len(RLOAD2 card) = {len(card):d}\ncard={card}'
        self.cards.append((sid, excite_id, delay, dphase, tb, tphi, load_type_str))
        self.n += 1
        return self.n

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        #excite_id : int
            #node id where the load is applied
        #delay : int/float; default=None
            #the delay; if it's 0/blank there is no delay
            #float : delay in units of time
            #int : delay id
        #dphase : int/float; default=None
            #the dphase; if it's 0/blank there is no phase lag
            #float : delay in units of time
            #int : delay id
        #tb : int/float; default=0
            #TABLEDi id that defines B(f) for all degrees of freedom in EXCITEID entry
        #tp : int/float; default=0
            #TABLEDi id that defines phi(f) for all degrees of freedom in EXCITEID entry
        #Type : int/str; default='LOAD'
            #the type of load
            #0/LOAD
            #1/DISP
            #2/VELO
            #3/ACCE
            #4, 5, 6, 7, 12, 13 - MSC only


        #: Set identification number
        load_id = np.zeros(ncards, dtype='int32')

        #: Identification number of DAREA or SPCD entry set or a thermal load
        #: set (in heat transfer analysis) that defines {A}. (Integer > 0)
        excite_id = np.zeros(ncards, dtype='int32')

        #: If it is a non-zero integer, it represents the
        #: identification number of DELAY Bulk Data entry that defines .
        #: If it is real, then it directly defines the value of that will
        #: be used for all degrees-of-freedom that are excited by this
        #: dynamic load entry.  See also Remark 9. (Integer >= 0,
        #: real or blank)
        delay_int = np.zeros(ncards, dtype='int32')
        delay_float = np.full(ncards, np.nan, dtype='float64')

        #: Defines the type of the dynamic excitation. (LOAD,DISP, VELO, ACCE)
        load_type = np.zeros(ncards, dtype='|U4')

        dphase_int = np.zeros(ncards, dtype='int32')
        dphase_float = np.full(ncards, np.nan, dtype='float64')

        # tc : int/float; default=0
        #     TABLEDi id that defines C(f) for all degrees of freedom in EXCITEID entry
        # td : int/float; default=0
        #     TABLEDi id that defines D(f) for all degrees of freedom in EXCITEID entry
        tabled_b_int = np.zeros(ncards, dtype='int32')
        tabled_phi_int = np.zeros(ncards, dtype='int32')

        tabled_b_float = np.full(ncards, np.nan, dtype='float64')
        tabled_phi_float = np.full(ncards, np.nan, dtype='float64')

        #self.tabled_b_id =
        #tabled_phi_id =

        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (sid, excite_idi, delay, dphase, tb, tphi, load_type_str) = card
            load_id[icard] = sid
            excite_id[icard] = excite_idi
            _set_int_float(icard, delay_int, delay_float, delay)
            _set_int_float(icard, dphase_int, dphase_float, dphase)
            _set_int_float(icard, tabled_b_int, tabled_b_float, tb)
            _set_int_float(icard, tabled_phi_int, tabled_phi_float, tphi)
            load_type[icard] = load_type_str

        self._save(load_id, excite_id, load_type,
                   delay_int, delay_float,
                   dphase_int, dphase_float,
                   tabled_b_int, tabled_b_float,
                   tabled_phi_int, tabled_phi_float)
        assert len(self.load_id) == self.n
        self.cards = []

    def _save(self, load_id, excite_id, load_type,
              delay_int, delay_float,
              dphase_int, dphase_float,
              tabled_b_int, tabled_b_float,
              tabled_phi_int, tabled_phi_float):
        self.load_id = load_id
        self.excite_id = excite_id
        self.load_type = load_type

        self.delay_int = delay_int
        self.delay_float = delay_float

        self.dphase_int = dphase_int
        self.dphase_float = dphase_float

        self.tabled_b_float = tabled_b_float
        self.tabled_b_int = tabled_b_int

        self.tabled_phi_int = tabled_phi_int
        self.tabled_phi_float = tabled_phi_float

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        # TODO: excite_id
        delay = self.delay_int[self.delay_int > 0]
        dphase = self.dphase_int[self.dphase_int > 0]
        used_dict['delay_id'].append(delay)
        used_dict['dphase_id'].append(dphase)
        for table in (self.tabled_b_int, self.tabled_phi_int):
            tabled = table[table > 0]
            used_dict['tabled_id'].append(tabled)

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(), self.excite_id.max(),
                   self.delay_int.max(), self.dphase_int.max(),
                   self.tabled_b_int.max(), self.tabled_phi_int.max())

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        #array_str, array_default_int
        load_ids = array_str(self.load_id, size=size)
        excite_ids = array_default_int(self.excite_id, default=0, size=size)

        delays = array_int_float(self.delay_int, self.delay_float, size=size, is_double=False)
        dphases = array_int_float(self.dphase_int, self.dphase_float, size=size, is_double=False)
        tbs = array_int_float(self.tabled_b_int, self.tabled_b_float, size=size, is_double=False)
        tphis = array_int_float(self.tabled_phi_int, self.tabled_phi_float, size=size, is_double=False)

        for sid, excite_id, load_type, delay, dphase, \
            tb, tphi in zip_longest(load_ids, excite_ids, self.load_type,
                                    delays, dphases,
                                    tbs, tphis):

            #dphase = _write_int_float(dphase_int, dphase_float)
            #delay = _write_int_float(delay_int, delay_float)
            #tb = _write_int_float(tb_int, tb_float)
            #tphi = _write_int_float(tphi_int, tphi_float)
            list_fields = ['RLOAD2', sid, excite_id, delay, dphase,
                           tb, tphi, load_type]
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        udelay = np.unique(self.delay_int)
        udphase = np.unique(self.dphase_int)

        all_delay = np.unique(model.delay.delay_id)
        all_dphase = np.unique(model.dphase.dphase_id)
        geom_check(
            self,
            missing,
            delay=(all_delay, udelay),
            dphase=(all_dphase, udphase),
        )


class LSEQ(VectorizedBaseCard):  # Requires LOADSET in case control deck
    """
    Defines a sequence of static load sets

    .. todo:: how does this work...
    +------+-----+----------+-----+-----+
    |   1  |  2  |     3    |  4  |  5  |
    +======+=====+==========+=====+=====+
    | LSEQ | SID | EXCITEID | LID | TID |
    +------+-----+----------+-----+-----+

    ACSRCE : If there is no LOADSET Case Control command, then EXCITEID
             may reference DAREA and SLOAD entries. If there is a LOADSET
             Case Control command, then EXCITEID may reference DAREA
             entries as well as SLOAD entries specified by the LID field
             in the selected LSEQ entry corresponding to EXCITEID.

    DAREA :  Refer to RLOAD1, RLOAD2, TLOAD1, TLOAD2, or ACSRCE entries
             for the formulas that define the scale factor Ai in dynamic
             analysis.

    DPHASE :

    SLOAD :  In the static solution sequences, the load set ID (SID) is
             selected by the Case Control command LOAD. In the dynamic
             solution sequences, SID must be referenced in the LID field
             of an LSEQ entry, which in turn must be selected by the Case
             Control command LOADSET.

    LSEQ LID : Load set identification number of a set of static load
               entries such as those referenced by the LOAD Case Control
               command.


    LSEQ,  SID, EXCITEID, LID, TID

    #--------------------------------------------------------------
    # F:\\Program Files\\Siemens\\NXNastran\\nxn10p1\\nxn10p1\\nast\\tpl\\cube_iter.dat

    DLOAD       1001     1.0     1.0   55212
    sid = 1001
    load_id = [55212] -> RLOAD2.SID

    RLOAD2,     SID, EXCITEID, DELAYID, DPHASEID,   TB,     TP,  TYPE
    RLOAD2     55212   55120              55122   55123   55124
    EXCITEID = 55120 -> DAREA.SID
    DPHASEID = 55122 -> DPHASE.SID

    DARA        SID      NID    COMP  SCALE
    DAREA      55120     913    3     9.9E+9
    SID = 55120 -> RLOAD2.SID

    DPHASE      SID     POINTID   C1    TH1
    DPHASE     55122     913       3   -90.0
    SID = 55122
    POINTID = 913 -> GRID.NID

    GRID       NID       X     Y     Z
    GRID       913      50.  0.19  -39.9
    """
    _id_name = 'lseq_id'
    def clear(self) -> None:
        self.n = 0
        self.lseq_id = np.array([], dtype='int32')
        self.excite_id = np.array([], dtype='int32')
        self.load_id = np.array([], dtype='int32')
        self.temperature_id = np.array([], dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> LSEQ:
        #load = LSEQ(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def __apply_slice__(self, load: LSEQ, i: np.ndarray) -> None:
        load.n = len(i)
        load.lseq_id = self.lseq_id[i]
        load.excite_id = self.excite_id[i]
        load.load_id = self.load_id[i]
        load.temp_id = self.temp_id[i]

    def add_card(self, card: BDFCard, comment: str=''):
        lseq_id = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id')
        load_id = integer_or_blank(card, 3, 'lid', default=0)
        temp_id = integer_or_blank(card, 4, 'tid', default=0)
        if max(load_id, temp_id) == 0:
            msg = 'LSEQ load_id/temp_id must not be None; load_id=%s temp_id=%s' % (load_id, temp_id)
            raise RuntimeError(msg)
        assert len(card) <= 5, f'len(LSEQ card) = {len(card):d}\ncard={card}'
        self.cards.append((lseq_id, excite_id, load_id, temp_id, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        lseq_id = np.zeros(ncards, dtype='int32')
        excite_id = np.zeros(ncards, dtype='int32')
        load_id = np.zeros(ncards, dtype='int32')
        temperature_id = np.zeros(ncards, dtype='int32')

        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (lseq_idi, excite_idi, load_idi, temp_idi, comment) = card
            lseq_id[icard] = lseq_idi
            excite_id[icard] = excite_idi
            load_id[icard] = load_idi
            temperature_id[icard] = temp_idi
        self._save(lseq_id, excite_id, load_id, temperature_id)
        assert len(self.load_id) == self.n, f'nloads={len(self.load_id)} n={self.n}'
        self.cards = []

    def _save(self, lseq_id, excite_id, load_id, temperature_id):
        if len(self.lseq_id) == 0:
            lseq_id = np.hstack([self.lseq_id, lseq_id])
            excite_id = np.hstack([self.excite_id, excite_id])
            load_id = np.hstack([self.load_id, load_id])
            temperature_id = np.hstack([self.temperature_id, temperature_id])
        nloads = len(lseq_id)
        self.lseq_id = lseq_id
        self.excite_id = excite_id
        self.load_id = load_id
        self.temp_id = temperature_id
        self.n = nloads

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(),
                   self.load_id.max(), self.excite_id.max(),
                   self.temp_id.max(), )

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        #self.get_loads()
        lseq_ids = array_str(self.lseq_id, size=size)
        load_ids = array_str(self.load_id, size=size)
        excite_ids = array_str(self.excite_id, size=size)
        load_ids = array_default_int(self.load_id, default=0, size=size)
        temp_ids = array_default_int(self.temp_id, default=0, size=size)
        for lseq_id, excite_id, load_id, temp_id in zip_longest(lseq_ids, excite_ids, load_ids, temp_ids):
            list_fields = ['LSEQ', lseq_id, excite_id, load_id, temp_id]
            bdf_file.write(print_card(list_fields))
        return

    #def get_loads_by_load_id(self) -> dict[int, Loads]:
        #return get_loads_by_load_id(self)

    #def get_reduced_loads(self,
                          #remove_missing_loads: bool=False,
                          #filter_zero_scale_factors: bool=False,
                          #stop_on_failure: bool=True) -> dict[int, Loads]:
        #return get_reduced_loads(
            #self, remove_missing_loads=remove_missing_loads,
            #filter_zero_scale_factors=filter_zero_scale_factors,
            #stop_on_failure=stop_on_failure)
    #def get_loads(self):
        #dynamic_loads = self.model.dynamic_loads
        #for card in


class TIC(VectorizedBaseCard):
    """
    Transient Initial Condition

    Defines values for the initial conditions of variables used in
    structural transient analysis. Both displacement and velocity
    values may be specified at independent degrees-of-freedom. This
    entry may not be used for heat transfer analysis.

    """
    _id_name = 'tic_id'
    def clear(self) -> None:
        self.n = 0
        self.tic_id = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')

    def add(self, sid: int, nid: int, components: int=0,
            u0: float=0.0, v0: float=0.0, comment: str='') -> int:
        assert isinstance(components, int), components
        self.cards.append((sid, nid, components, u0, v0, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        sid = integer(card, 1, 'sid')
        nid = integer(card, 2, 'G')
        comp = modal_components_or_blank(card, 3, 'C', default=0)
        u0 = double_or_blank(card, 4, 'U0', default=0.)
        v0 = double_or_blank(card, 5, 'V0', default=0.)

        self.cards.append((sid, nid, comp, u0, v0, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug('parse TIC')
        idtype = self.model.idtype
        tic_id = np.zeros(ncards, dtype='int32')
        node_id = np.zeros(ncards, dtype=idtype)
        component = np.zeros(ncards, dtype='int32')
        u0 = np.zeros(ncards, dtype='float64')
        v0 = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (sid, nid, comp, u0i, v0i, comment) = card
            tic_id[icard] = sid
            node_id[icard] = nid
            component[icard] = comp
            u0[icard] = u0i
            v0[icard] = v0i
            #self.comment[i] = comment
        self._save(tic_id, node_id, component, u0, v0)
        #self.sort()
        self.cards = []

    def _save(self, tic_id, node_id, component, u0, v0):
        assert len(self.tic_id) == 0, self.tic_id
        self.tic_id = tic_id
        self.node_id = node_id
        self.component = component
        self.u0 = u0
        self.v0 = v0

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.node_id)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.node_id
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    #def sort(self) -> None:
        #iarg = np.argsort(self.node_id)
        #uarg = np.unique(iarg)
        ##nvalues = len(self.node_id)
        #if not np.array_equal(uarg, iarg):
            #self.node_id = self.node_id[iarg]
            #self.cp = self.cp[iarg]
            #self.xyz = self.xyz[iarg, :]
            #self.ps = self.ps[iarg]
            #self.cd = self.cd[iarg]
            #self.seid = self.seid[iarg]
        #self._is_sorted = True

    #def index(self, node_id: np.ndarray) -> np.ndarray:
        #assert len(self.node_id) > 0, self.node_id
        #node_id = np.atleast_1d(np.asarray(node_id, dtype=self.node_id.dtype))
        #inid = np.searchsorted(self.node_id, node_id)
        #return inid

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        unode = np.unique(self.node_id)

        geom_check(
            self,
            missing,
            node=(model.grid.node_id, unode),
        )

    #def slice_by_node_id(self, node_id: np.ndarray) -> GRID:
        #inid = self._node_index(node_id)
        #return self.slice_card(inid)

    #def slice_card_by_node_id(self, node_id: np.ndarray) -> TIC:
        #"""uses a node_ids to extract GRIDs"""
        #inid = self.index(node_id)
        ##assert len(self.node_id) > 0, self.node_id
        ##i = np.searchsorted(self.node_id, node_id)
        #tic = self.slice_card_by_index(inid)
        #return grid

    def slice_card_by_index(self, i: np.ndarray) -> TIC:
        """uses a node_index to extract TICs"""
        #assert self.xyz.shape == self._xyz_cid0.shape
        assert len(self.node_id) > 0, self.node_id
        i = np.atleast_1d(np.asarray(i, dtype=self.node_id.dtype))
        i.sort()
        tic = TIC(self.model)
        self.__apply_slice__(tic, i)
        return tic

    def __apply_slice__(self, tic: TIC, i: np.ndarray) -> None:
        tic.n = len(i)
        #tic._is_sorted = self._is_sorted
        tic.node_id = self.node_id[i]
        tic.tic_id = self.tic_id[i]
        tic.component = self.component[i]
        tic.u0 = self.u0[i]
        tic.v0 = self.v0[i]

    @property
    def max_id(self) -> int:
        return max(self.tic_id.max(), self.node_id.max(),)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.tic_id) == 0:
            return

        print_card, size = get_print_card_size(size, self.max_id)
        tic_ids = array_str(self.tic_id, size=size)
        node_id = array_str(self.node_id, size=size)
        components = array_default_int(self.component, size=size, default=0)

        for sid, nid, comp, u0, v0 in zip_longest(tic_ids, node_id, components, self.u0, self.v0):
            list_fields = ['TIC', sid, nid, comp, u0, v0]
            bdf_file.write(print_card(list_fields))
        return


class TF(VectorizedBaseCard):
    """
    Defines a dynamic transfer function of the form:
        (B0 + B1 p + B2 *p2)*ud  sum(A0_i + A1_i*p + A2_i*p2)*ui = 0

    +----+-----+-----+------+------+------+--------+----+----+
    |  1 |  2  |  3  |   4  |   5  |   6  |    7   |  8 |  9 |
    +====+=====+=====+======+======+======+========+====+====+
    | TF | SID | GD  |  CD  |  B0  |  B1  |   B2   |    |    |
    +----+-----+-----+------+------+------+--------+----+----+
    |    | G_1 | C_1 | A0_1 | A1_1 | A2_1 |  etc.  |    |    |
    +----+-----+-----+------+------+------+--------+----+----+

    """
    _id_name = 'tf_id'
    def clear(self) -> None:
        self.n = 0
        self.tf_id = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')
        self.component = np.array([], dtype='int32')
        self.b = np.zeros((0, 3), dtype='float64')

        self.nodes = np.array([], dtype='int32')
        self.components = np.array([], dtype='int32')
        self.a = np.zeros((0, 3), dtype='float64')

    def add(self, tf_id: int,
            nid: int, component: int,
            b0: float, b1: float, b2: float,
            nids: list[int], components: list[int],
            a: list[tuple[float, float, float]],
            comment: str='') -> int:
        """Creates a TF card"""
        assert isinstance(component, int), component
        nnids = len(nids)
        assert nnids == len(components), f'nnids={nnids} ncomponents={len(components)}'
        assert nnids == len(a), f'nnids={nnids} na={len(a)}'

        self.cards.append((tf_id, nid, component, (b0, b1, b2),
                           a, nids, components, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a TF card from ``BDF.add_card(...)``

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
        nid0 = integer(card, 2, 'nid0')
        # component 0 means an SPOINT/EPOINT
        comp = components_or_blank(card, 3, 'components_0', default='0')
        b0 = double_or_blank(card, 4, 'b0', default=0.)
        b1 = double_or_blank(card, 5, 'b1', default=0.)
        b2 = double_or_blank(card, 6, 'b2', default=0.)

        nfields = len(card) - 9
        nrows = nfields // 8
        if nfields % 8 > 0:
            nrows += 1

        nids = []
        components = []
        a = []
        for irow in range(nrows):
            j = irow * 8 + 9
            #ifield = irow + 1
            nid = integer(card, j, 'grid_%i' % (irow + 1))
            component = components_or_blank(card, j + 1, 'components_%i' % (irow + 1), '0')
            a0 = double_or_blank(card, j + 2, 'a0_%i' % (irow + 1), 0.)
            a1 = double_or_blank(card, j + 3, 'a1_%i' % (irow + 1), 0.)
            a2 = double_or_blank(card, j + 4, 'a2_%i' % (irow + 1), 0.)
            nids.append(nid)
            components.append(component)
            a.append([a0, a1, a2])
        #return TF(sid, nid0, comp, [b0, b1, b2], nids, components, a,
                  #comment=comment)

        self.cards.append((sid, nid0, comp, (b0, b1, b2),
                           a, nids, components, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug(f'parse {self.type}')
        #idtype = self.model.idtype

        tf_id = np.zeros(ncards, dtype='int32')
        node_id = np.zeros(ncards, dtype='int32')
        component = np.zeros(ncards, dtype='int32')
        b = np.zeros((ncards, 3), dtype='float64')
        #a = np.zeros((ncards, 3), dtype='float64')

        a_list = []
        components_list = []
        nodes_list = []
        nnode = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (sid, nid0, comp, bi,
             ai, nidsi, componentsi, comment) = card
            tf_id[icard] = sid
            node_id[icard] = nid0
            component[icard] = comp
            b[icard] = bi

            for aii in ai:
                assert len(aii) == 3, ai
                a_list.append(aii)
            nnode[icard] = len(nidsi)
            nodes_list.extend(nidsi)
            components_list.extend(componentsi)
            #self.comment[icard] = comment
        a = np.array(a_list, dtype='float64')
        nodes = np.array(nodes_list, dtype='int32')
        components = np.array(components_list, dtype='int32')
        self._save(tf_id, node_id, component, b,
                   a, components, nodes, nnode)
        #self.sort()
        self.cards = []

    def _save(self, tf_id, node_id, component, b,
              a, components, nodes, nnode):
        assert len(self.tf_id) == 0, self.tf_id
        self.tf_id = tf_id
        self.node_id = node_id
        self.component = component
        self.b = b

        self.a = a
        self.components = components
        self.nodes = nodes
        self.nnode = nnode

    def __apply_slice__(self, tf: TF, i: np.ndarray) -> None:
        tf.n = len(i)
        #tf._is_sorted = self._is_sorted

        tf.tf_id = self.tf_id[i]
        tf.node_id = self.node_id[i]
        tf.component = self.component[i]
        tf.b = self.b[i, :]
        #asdf
        #tf.nnode = self.nnode[i]
        #if self.nnode.max() > 0:
        inode = self.inode
        tf.nodes = hslice_by_idim(i, inode, self.nodes)  # was vslice
        tf.components = hslice_by_idim(i, inode, self.components)  # was vslice
        tf.a = vslice_by_idim(i, inode, self.a)
        #tf.a = self.a[i, :]
        tf.nnode = self.nnode[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.node_id)
        used_dict['node_id'].append(self.nodes.ravel())

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        tf_id = used_dict['tf_id']
        ncards_removed = remove_unused_primary(
            self, tf_id, self.tf_id, 'tf_id')
        return ncards_removed

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        nodes = np.hstack([self.node_id, self.nodes.ravel()])
        unode = np.unique(nodes)

        geom_check(
            self,
            missing,
            node=(model.grid.node_id, unode),
        )

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.nodes.ravel()
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    @property
    def inode(self) -> np.ndarray:
        return make_idim(self.n, self.nnode)

    @property
    def max_id(self) -> int:
        tf_id = self.tf_id.max()
        node_id = self.node_id.max()
        nodes = 1
        if len(self.nodes):
            nodes = self.nodes.max()
        return max(tf_id, node_id, nodes)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.tf_id) == 0:
            return

        print_card, size = get_print_card_size(size, self.max_id)
        tf_ids = array_str(self.tf_id, size=size)
        node_id = array_str(self.node_id, size=size)
        #component = array_default_int(self.component, size=size, default=0)

        #self.a = a
        #self.self.components = components
        #self.self.nodes = nodes
        #self.self.nnode = nnode
        bb_ = array_float(self.b, size=size, is_double=False)
        component = array_str(self.component, size=size)

        nids_ = array_str(self.nodes, size=size)
        aa_ = array_float(self.b, size=size, is_double=False)
        components_ = array_str(self.components, size=size)
        for tf_id, nid, comp, (b0, b1, b2), (inode0, inode1) in zip_longest(tf_ids, node_id, component,
                                                                            bb_, self.inode):
            list_fields = ['TF', tf_id, nid, comp, b0, b1, b2, None, None]

            if inode0 != inode1:
                nidsi = nids_[inode0:inode1]
                componentsi = components_[inode0:inode1]
                aa = aa_[inode0:inode1, :]

                nidsi = self.nodes[inode0:inode1]
                componentsi = self.components[inode0:inode1]
                aa = self.a[inode0:inode1, :]
                for nid1, comp1, (a0, a1, a2) in zip(nidsi, componentsi, aa):
                    list_fields += [nid1, comp1, a0, a1, a2, None, None, None]

                #list_fields = ['TF', self.sid, self.nid0, self.c, self.b0, self.b1, self.b2, None, None]
                #for grid, c, (a0, a1, a2) in zip(self.nids, self.components, aa):
                    #list_fields += [grid, c, a0, a1, a2, None, None, None]
            bdf_file.write(print_card(list_fields))
        return


class DELAY(VectorizedBaseCard):
    """
    +-------+-----+-----------+-----+--------+------+-----+--------+
    |   1   |  2  |     3     |  4  |   5    |  6   |  7  |   8    |
    +=======+=====+===========+=====+========+======+=====+========+
    | DELAY | SID | POINT ID1 | C1  |   T1   | P2   | C2  |   T2   |
    +-------+-----+-----------+-----+--------+------+-----+--------+
    """
    #def __init__(self, model: BDF):
        #super().__init__(model)
        #self._is_sorted = False
    _id_name = 'delay_id'
    def clear(self) -> None:
        self.n = 0
        self.delay_id = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')

    #def slice_by_node_id(self, node_id: np.ndarray) -> GRID:
        #inid = self._node_index(node_id)
        #return self.slice_card(inid)

    #def slice_card_by_node_id(self, node_id: np.ndarray) -> TIC:
        #"""uses a node_ids to extract GRIDs"""
        #inid = self.index(node_id)
        ##assert len(self.node_id) > 0, self.node_id
        ##i = np.searchsorted(self.node_id, node_id)
        #tic = self.slice_card_by_index(inid)
        #return grid

    #def slice_card_by_index(self, i: np.ndarray) -> DELAY:
        #"""uses a node_index to extract DELAYs"""
        #assert len(self.node_id) > 0, self.node_id
        #i = np.atleast_1d(np.asarray(i, dtype=self.node_id.dtype))
        #i.sort()
        #delay = DELAY(self.model)
        #self.__apply_slice__(delay, i)
        #return delay

    def add(self, sid: int, nid: int, component: int=0,
            delay: float=0.0, comment: str='') -> int:
        if isinstance(nid, integer_types) and isinstance(component, integer_types):
            self.cards.append((sid, nid, component, delay, comment))
            self.n += 1
        else:
            assert isinstance(nid, list), nid
            assert isinstance(component, list), component
            assert isinstance(delay, list), delay
            for nidi, compi, delayi in zip(nid, component, delay):
                self.cards.append((sid, nidi, compi, delayi, comment))
                self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DELAY card from ``BDF.add_card(...)``

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
        node = integer(card, 2, 'node')
        component = integer_or_blank(card, 3, 'components', default=0)
        delay = double_or_blank(card, 4, 'delay')
        assert component in {0, 1, 2, 3, 4, 5, 6}, component
        self.cards.append((sid, node, component, delay, comment))
        self.n += 1

        if card.field(5):
            node = integer(card, 5, 'node')
            component = integer_or_blank(card, 6, 'components', default=0)
            delay = double_or_blank(card, 7, 'delay')
            assert component in {0, 1, 2, 3, 4, 5, 6}, component
            self.cards.append((sid, node, component, delay, ''))
            self.n += 1
        #return DELAY(sid, nodes, components, delays, comment=comment)
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug('parse DELAY')
        delay_id = np.zeros(ncards, dtype='int32')
        node_id = np.zeros(ncards, dtype='int32')
        component = np.zeros(ncards, dtype='int32')
        delay = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (delay_idi, nodei, componenti, delayi, comment) = card
            delay_id[icard] = delay_idi
            node_id[icard] = nodei
            component[icard] = componenti
            delay[icard] = delayi
            #self.comment[i] = comment
        self._save(delay_id, node_id, component, delay)
        #self.sort()
        self.cards = []

    def _save(self, delay_id, node_id, component, delay):
        assert len(self.delay_id) == 0, self.delay_id
        self.delay_id = delay_id
        self.node_id = node_id
        self.component = component
        self.delay = delay

    def __apply_slice__(self, delay: DELAY, i: np.ndarray) -> None:
        delay.n = len(i)
        delay.node_id = self.node_id[i]
        delay.delay_id = self.delay_id[i]
        delay.delay = self.delay[i]
        delay.component = self.component[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.node_id)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        delay_id = used_dict['delay_id']
        ncards_removed = remove_unused_duplicate(
            self, delay_id, self.delay_id, 'delay_id')
        return ncards_removed

    @property
    def max_id(self) -> int:
        return max(self.delay_id.max(), self.node_id.max(),)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.delay_id) == 0:
            return

        print_card, size = get_print_card_size(size, self.max_id)
        delay_ids = array_str(self.delay_id, size=size)
        node_id = array_str(self.node_id, size=size)
        components = array_default_int(self.component, size=size, default=0)

        for sid, nid, comp, delay in zip_longest(delay_ids, node_id, components, self.delay):
            list_fields = ['DELAY', sid, nid, comp, delay]
            bdf_file.write(print_card(list_fields))
        return

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.node_id
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    #def sort(self) -> None:
        #iarg = np.argsort(self.node_id)
        #uarg = np.unique(iarg)
        ##nvalues = len(self.node_id)
        #if not np.array_equal(uarg, iarg):
            #self.node_id = self.node_id[iarg]
            #self.cp = self.cp[iarg]
            #self.xyz = self.xyz[iarg, :]
            #self.ps = self.ps[iarg]
            #self.cd = self.cd[iarg]
            #self.seid = self.seid[iarg]
        #self._is_sorted = True

    #def index(self, node_id: np.ndarray) -> np.ndarray:
        #assert len(self.node_id) > 0, self.node_id
        #node_id = np.atleast_1d(np.asarray(node_id, dtype=self.node_id.dtype))
        #inid = np.searchsorted(self.node_id, node_id)
        #return inid

class DPHASE(VectorizedBaseCard):
    """
    Defines the phase lead term θ in the equation of the dynamic
    loading function.

    +--------+-----+-----------+-----+------+------+-----+-----+
    |   1    |  2  |     3     |  4  |  5   |  6   |  7  |  8  |
    +========+=====+===========+=====+======+======+=====+=====+
    | DPHASE | SID | POINT ID1 | C1  | TH1  |  P2  | C2  | TH2 |
    +--------+-----+-----------+-----+------+------+-----+-----+

    """
    #def __init__(self, model: BDF):
        #super().__init__(model)
        #self._is_sorted = False

    _id_name = 'dphase_id'
    def clear(self) -> None:
        self.n = 0
        self.dphase_id = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')
        self.component = np.array([], dtype='int32')

        # sid : int
        #     DPHASE id that is referenced by a RLOADx or ACSRCE card
        # nodes : list[int]
        #     list of nodes that see the delay
        #     len(nodes) = 1 or 2
        # components : list[int]
        #     the components corresponding to the nodes that see the delay
        #     len(nodes) = len(components)
        # phase_leads : list[float]
        #     Phase lead θ in degrees.
        #     len(nodes) = len(delays)
        # comment : str; default=''
        #     a comment for the card

    #def slice_by_node_id(self, node_id: np.ndarray) -> GRID:
        #inid = self._node_index(node_id)
        #return self.slice_card(inid)

    #def slice_card_by_node_id(self, node_id: np.ndarray) -> TIC:
        #"""uses a node_ids to extract GRIDs"""
        #inid = self.index(node_id)
        ##assert len(self.node_id) > 0, self.node_id
        ##i = np.searchsorted(self.node_id, node_id)
        #tic = self.slice_card_by_index(inid)
        #return grid

    #def slice_card_by_index(self, i: np.ndarray) -> DPHASE:
        #"""uses a node_index to extract DELAYs"""
        #assert len(self.node_id) > 0, self.node_id
        #i = np.atleast_1d(np.asarray(i, dtype=self.node_id.dtype))
        #i.sort()
        #dphase = DPHASE(self.model)
        #self.__apply_slice__(dphase, i)
        #return dphase

    def add(self, sid: int, nid: int, component: int=0,
            phase_lead: float=0.0, comment: str=''):
        if isinstance(nid, integer_types) and isinstance(component, integer_types):
            self.cards.append((sid, nid, component, phase_lead, comment))
            self.n += 1
        else:
            assert isinstance(nid, list), nid
            assert isinstance(component, list), component
            assert isinstance(phase_lead, list), phase_lead
            for nidi, compi, phase_leadi in zip(nid, component, phase_lead):
                self.cards.append((sid, nidi, compi, phase_leadi, comment))
                self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DPHASE card from ``BDF.add_card(...)``

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
        node = integer(card, 2, 'node0')
        component = integer_or_blank(card, 3, 'component0', default=0)
        phase_lead = double_or_blank(card, 4, 'phase_lead')
        assert component in {0, 1, 2, 3, 4, 5, 6}, component
        self.cards.append((sid, node, component, phase_lead, comment))
        self.n += 1

        if card.field(5):
            node = integer(card, 5, 'node1')
            component = integer_or_blank(card, 6, 'component1', default=0)
            phase_lead = double_or_blank(card, 7, 'phase_lead1')
            assert component in {0, 1, 2, 3, 4, 5, 6}, component
            self.cards.append((sid, node, component, phase_lead, ''))
            self.n += 1
        #return DPHASE(sid, nodes, components, phase_leads, comment=comment)
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug('parse DPHASE')
        idtype = self.model.idtype
        dphase_id = np.zeros(ncards, dtype='int32')
        node_id = np.zeros(ncards, dtype=idtype)
        component = np.zeros(ncards, dtype='int32')
        phase_lead = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (dphase_idi, nodei, componenti, phase_leadi, comment) = card
            dphase_id[icard] = dphase_idi
            node_id[icard] = nodei
            component[icard] = componenti
            phase_lead[icard] = phase_leadi
            #self.comment[i] = comment
        self._save(dphase_id, node_id, component, phase_lead)
        #self.sort()
        self.cards = []

    def _save(self, dphase_id, node_id, component, phase_lead):
        assert len(self.dphase_id) == 0, self.dphase_id
        self.dphase_id = dphase_id
        self.node_id = node_id
        self.component = component
        self.phase_lead = phase_lead

    def __apply_slice__(self, dphase: DPHASE, i: np.ndarray) -> None:
        dphase.n = len(i)
        #dphase._is_sorted = self._is_sorted
        dphase.dphase_id = self.dphase_id[i]
        dphase.node_id = self.node_id[i]
        dphase.component = self.component[i]
        dphase.phase_lead = self.phase_lead[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.node_id)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        dphase_id = used_dict['dphase_id']
        ncards_removed = remove_unused_duplicate(
            self, dphase_id, self.dphase_id, 'dphase_id')
        return ncards_removed

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.node_id
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    #def sort(self) -> None:
        #iarg = np.argsort(self.node_id)
        #uarg = np.unique(iarg)
        ##nvalues = len(self.node_id)
        #if not np.array_equal(uarg, iarg):
            #self.node_id = self.node_id[iarg]
            #self.cp = self.cp[iarg]
            #self.xyz = self.xyz[iarg, :]
            #self.ps = self.ps[iarg]
            #self.cd = self.cd[iarg]
            #self.seid = self.seid[iarg]
        #self._is_sorted = True

    #def index(self, node_id: np.ndarray) -> np.ndarray:
        #assert len(self.node_id) > 0, self.node_id
        #node_id = np.atleast_1d(np.asarray(node_id, dtype=self.node_id.dtype))
        #inid = np.searchsorted(self.node_id, node_id)
        #return inid

    @property
    def max_id(self) -> int:
        return max(self.dphase_id.max(), self.node_id.max(),)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dphase_id) == 0:
            return

        print_card, size = get_print_card_size(size, self.max_id)
        dphase_ids = array_str(self.dphase_id, size=size)
        node_id = array_str(self.node_id, size=size)
        components = array_default_int(self.component, size=size, default=0)

        for sid, nid, comp, phase_lead in zip_longest(dphase_ids, node_id, components, self.phase_lead):
            list_fields = ['DPHASE', sid, nid, comp, phase_lead]
            bdf_file.write(print_card(list_fields))
        return


class QVECT(VectorizedBaseCard):
    """
    Thermal Vector Flux Load

    Defines thermal vector flux from a distant source into a face of one
    or more CHBDYi boundary condition surface elements.

    +-------+------+------+-------+-----+---------+---------+---------+---------+
    |   1   |   2  |   3  |   4   |  5  |     6   |    7    |    8    |    9    |
    +=======+======+======+=======+=====+=========+=========+=========+=========+
    | QVECT | SID  |  Q0  | TSOUR | CE  | E1/TID1 | E2/TID2 | E3/TID3 | CNTRLND |
    +-------+------+------+-------+-----+---------+---------+---------+---------+
    |       | EID1 | EID2 |  etc. |     |         |         |         |         |
    +-------+------+------+-------+-----+---------+---------+---------+---------+

    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> QVECT:
        #load = QVECT(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def add(self, sid: int, q0: float, eids: list[int],
            t_source: float=None,
            ce: int=0,
            vector_tableds: list[Union[int, float]]=0.0,
            control_id: int=0, comment: str='') -> int:
        """
        Creates a QVECT card

        Parameters
        ----------
        sid : int
            Load set identification number. (Integer > 0)
        q0 : float; default=None
            Magnitude of thermal flux vector into face
        t_source : float; default=None
            Temperature of the radiant source
        ce : int; default=0
            Coordinate system identification number for thermal vector flux
        vector_tableds : list[int/float, int/float, int/float]
            vector : float; default=0.0
                directional cosines in coordinate system CE) of
                the thermal vector flux
            tabled : int
                TABLEDi entry identification numbers defining the
                components as a function of time
        control_id : int; default=0
            Control point
        eids : list[int] or THRU
            Element identification number of a CHBDYE, CHBDYG, or
            CHBDYP entry
        comment : str; default=''
            a comment for the card

        """

    def add_card(self, card: BDFCard, comment: str=''):
        sid = integer(card, 1, 'sid')
        q0 = double(card, 2, 'q0')
        t_source = double_or_blank(card, 3, 't_source')
        ce = integer_or_blank(card, 4, 'ce', default=0)
        vector_tableds = [
            integer_double_or_blank(card, 5, 'e1_tabled1', default=0.0),
            integer_double_or_blank(card, 6, 'e2_tabled2', default=0.0),
            integer_double_or_blank(card, 7, 'e3_tabled3', default=0.0),
        ]
        control_id = integer_or_blank(card, 8, 'control_id', default=0)

        i = 1
        eids = []
        for ifield in range(9, len(card)):
            eid = integer_or_string(card, ifield, 'eid_%d' % i)
            eids.append(eid)
            assert eid != 0, card
            i += 1
        elements = expand_thru_by(eids)
        self.cards.append((sid, q0, t_source, control_id, ce, vector_tableds, elements, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)

        #: Set identification number
        load_id = np.zeros(ncards, dtype='int32')
        control_id = np.zeros(ncards, dtype='int32')
        nelement = np.zeros(ncards, dtype='int32')
        q0 = np.full(ncards, np.nan, dtype='float64')
        t_source = np.full(ncards, np.nan, dtype='float64')
        ce = np.zeros(ncards, dtype='int32')
        vector = np.full((ncards, 3), np.nan, dtype='float64')
        tableds = np.full((ncards, 3), 0, dtype='int32')

        assert ncards > 0, ncards
        all_elements = []
        for icard, card in enumerate(self.cards):
            (sid, q0i, t_sourcei, control_idi, cei, vector_tabledsi, elementsi, comment) = card
            load_id[icard] = sid
            q0[icard] = q0i
            t_source[icard] = t_sourcei
            nelement[icard] = len(elementsi)
            control_id[icard] = control_idi
            ce[icard] = cei
            for i, vector_tabled in enumerate(vector_tabledsi):
                if isinstance(vector_tabled, int):
                    tableds[icard, i] = vector_tabled
                else:
                    vector[icard, i] = vector_tabled
            all_elements.extend(elementsi)

        elements = np.array(all_elements, dtype='int32')
        self._save(load_id, q0, t_source, control_id, ce, vector, tableds,
                   elements, nelement)
        assert len(self.load_id) == self.n
        self.cards = []

    def _save(self, load_id, q0, t_source, control_id, ce, vector, tableds,
              element, nelement):
        if len(self.load_id) != 0:
            adf
            self.load_id
        nloads = len(load_id)
        self.load_id = load_id
        self.q0 = q0
        self.t_source = t_source
        self.control_id = control_id
        self.ce = ce
        self.vector = vector
        self.tableds = tableds
        self.element = element
        assert len(nelement) > 0, nelement
        self.nelement = nelement
        self.n = nloads

    def __apply_slice__(self, load: QVECT, i: np.ndarray) -> None:
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.control_id = self.control_id[i]
        load.q0 = self.q0[i]
        load.t_source = self.t_source[i]
        load.ce = self.ce[i]
        load.vector = self.vector[i, :]
        load.tableds = self.tableds[i, :]

        load.element = hslice_by_idim(i, self.ielement, self.element)
        load.nelement = self.nelement[i]

    @property
    def ielement(self) -> np.ndarray:
        return make_idim(self.n, self.nelement)

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        all_node_ids = self.model.grid.node_id
        all_coord_ids = self.model.coord.coord_id
        uce = np.unique(self.ce)
        utabled = np.unique(self.tableds)
        uthermal_element = np.unique(self.element)
        unode = np.unique(self.control_id)

        geom_check(
            self,
            missing,
            node=(all_node_ids, unode),
            coord=(all_coord_ids, uce),
            #thermal_element=(all_thermal_elements, uthermal_element),
        )

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(),
                   self.control_id.max(),
                   self.element.max(),)

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        #array_str, array_default_int

        load_ids = array_str(self.load_id, size=size)
        control_id = array_default_int(self.control_id, size=size)
        for sid, q0, t_source, ce, vector, tableds, \
            control_id, (ieid0, ieid1) in zip_longest(load_ids, self.q0, self.t_source,
                                                      self.ce, self.vector, self.tableds,
                                                      control_id, self.ielement):
            element_ids = self.element[ieid0:ieid1].tolist()
            eids = collapse_thru_by(element_ids)
            vector_tableds = []
            for vectori, tabled in zip(vector, tableds):
                if tabled == 0:
                    vector_tableds.append(vectori)
                else:
                    vector_tableds.append(tabled)
            list_fields = [
                'QVECT', sid, q0, t_source, ce
                ] + vector_tableds + [control_id] + eids
            bdf_file.write(print_card(list_fields))
        return

    def sum_forces_moments(self) -> np.ndarray:
        nloads = len(self.load_id)
        force_moment = np.zeros((nloads, 6), dtype='float64')
        return force_moment


class ACSRCE(VectorizedBaseCard):
    r"""
    Defines acoustic source as a function of power vs. frequency.

    +--------+-----+----------+---------------+-----------------+-------+-----+---+
    |   1    |  2  |    3     |       4       |        5        |   6   |  7  | 8 |
    +========+=====+==========+===============+=================+=======+=====+===+
    | ACSRCE | SID | EXCITEID | DELAYI/DELAYR | DPHASEI/DPHASER | TP/RP | RHO | B |
    +--------+-----+----------+---------------+-----------------+-------+-----+---+

    ..math ::
      C = \sqrt(B ⁄ ρ)
      Source Strength = {A} * 1/(2πf)  * \sqrt( 8πC P(f) / ρ) ^ (ei(θ + 2πfτ))

    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')

    def add(self, sid: int, excite_id: int, rho: float, b: float,
            delay: Union[int, float]=0,
            dphase: Union[int, float]=0,
            power: Union[int, float]=0,
            comment: str='') -> int:
        """
        Creates an ACSRCE card

        Parameters
        ----------
        sid : int
            load set id number (referenced by DLOAD)
        excite_id : int
            Identification number of a DAREA or SLOAD entry that lists
            each degree of freedom to apply the excitation and the
            corresponding scale factor, A, for the excitation
        rho : float
            Density of the fluid
        b : float
            Bulk modulus of the fluid
        delay : int; default=0
            Time delay, τ.
        dphase : int / float; default=0
            the dphase; if it's 0/blank there is no phase lag
            float : delay in units of time
            int : delay id
        power : int; default=0
            Power as a function of frequency, P(f).
            float : value of P(f) used over all frequencies for all
                    degrees of freedom in EXCITEID entry.
            int : TABLEDi entry that defines P(f) for all degrees of
                  freedom in EXCITEID entry.
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, excite_id,
                           delay, dphase, power, rho, b, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a ACSRCE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        excite_id = integer(card, 2, 'excite_id') # DAREA, FBALOAD, SLOAD
        delay = integer_double_or_blank(card, 3, 'delay', default=0) # DELAY, FBADLAY
        dphase = integer_double_or_blank(card, 4, 'dphase', default=0) # DPHASE, FBAPHAS
        power = integer_double_or_blank(card, 5, 'power/tp/rp', default=0) # TABLEDi/power
        rho = double(card, 6, 'rho')
        b = double(card, 7, 'bulk modulus')

        assert len(card) <= 8, 'len(ACSRCE card) = %i\n%s' % (len(card), card)
        #return ACSRCE(sid, excite_id, rho, b,
                      #delay=delay, dphase=dphase, power=power, comment=comment)
        self.cards.append((sid, excite_id,
                           delay, dphase, power, rho, b, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)

        #: Set identification number
        load_id = np.zeros(ncards, dtype='int32')
        excite_id = np.zeros(ncards, dtype='int32')

        delay_float = np.zeros(ncards, dtype='float64')
        delay_int = np.zeros(ncards, dtype='int32')

        dphase_int = np.zeros(ncards, dtype='int32')
        dphase_float = np.zeros(ncards, dtype='float64')

        power_int = np.zeros(ncards, dtype='int32')
        power_float = np.zeros(ncards, dtype='float64')

        rho = np.zeros(ncards, dtype='float64')
        b = np.zeros(ncards, dtype='float64')

        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (sid, excite_idi,
             delay, dphase, power, rhoi, bi, comment) = card
            load_id[icard] = sid
            excite_id[icard] = excite_idi
            _set_int_float(icard, delay_int, delay_float, delay)
            _set_int_float(icard, dphase_int, dphase_float, dphase)
            _set_int_float(icard, delay_int, delay_float, delay)
            _set_int_float(icard, power_int, power_float, power)

            b[icard] = bi
            rho[icard] = rhoi

        self._save(load_id, excite_id, b, rho,
                   delay_int, delay_float,
                   dphase_int, dphase_float,
                   power_int, power_float)
        assert len(self.load_id) == self.n
        self.cards = []

    def _save(self, load_id, excite_id, b, rho,
              delay_int, delay_float,
              dphase_int, dphase_float,
              power_int, power_float,) -> None:
        if len(self.load_id) != 0:
            sadf
        self.load_id = load_id
        self.excite_id = excite_id
        self.b = b
        self.rho = rho

        self.power_int = power_int
        self.delay_int = delay_int
        self.dphase_int = dphase_int

        self.power_float = power_float
        self.delay_float = delay_float
        self.dphase_float = dphase_float

    def __apply_slice__(self, load: ACSRCE, i: np.ndarray) -> None:
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.excite_id = self.excite_id[i]
        load.b = self.b[i]
        load.rho = self.rho[i]

        load.power_int = self.power_int[i]
        load.delay_int = self.delay_int[i]
        load.dphase_int = self.dphase_int[i]

        load.power_float = self.power_float[i]
        load.delay_float = self.delay_float[i]
        load.dphase_float = self.dphase_float[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        #used_dict['???'].append(self.excite_id)
        used_dict['delay_id'].append(self.delay_int[self.delay_int > 0])
        used_dict['dphase_id'].append(self.dphase_int[self.dphase_int > 0])

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        load_id = used_dict['load_id']
        ncards_removed = remove_unused_duplicate(
            self, load_id, self.load_id, 'load_id')
        return ncards_removed

    def get_load_at_freq(self, freq: np.ndarray) -> np.ndarray:  # pragma: no cover
        r"""
        ..math ::
          C = \sqrt(B ⁄ ρ)
          Source_strength = {A} * 1/(2πf)  * \sqrt( 8πC P(f) / ρ) ^ (ei(θ + 2πfτ))
        """
        assert isinstance(freq, np.ndarray), freq

        C = np.sqrt(self.b / self.rho)
        ei = np.exp(1) * 1.j
        A = 0.0
        pi = np.pi

        if self.delay in [0, 0.]:
            tau = 0.
        else:
            #print('delay\n', self.delay_ref)
            tau = self.delay_ref.value

        Pf = self.power_ref.interpolate(freq)
        if self.dphase in [0, 0.]:
            theta = 0.
        else:
            #print('dphase\n', self.dphase_ref)
            theta = self.dphase_ref.interpolate(freq)

        omega = 2 * pi * freq[:, np.newaxis]
        strength = A / omega * np.sqrt(8*pi*C*Pf / self.rho) ** (ei*(theta + omega*tau))
        return strength

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(), self.excite_id.max(),
                   self.delay_int.max(), self.dphase_int.max(),
                   self.power_int.max())

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        load_ids = array_str(self.load_id, size=size)
        excite_ids = array_str(self.excite_id, size=size)

        delays = array_int_float(self.delay_int, self.delay_float, size=size, is_double=False)
        dphases = array_int_float(self.dphase_int, self.dphase_float, size=size, is_double=False)
        powers = array_int_float(self.power_int, self.power_float, size=size, is_double=False)
        bs = array_float(self.b, size=size, is_double=False)
        rhos = array_float(self.rho, size=size, is_double=False)        for sid, excite_id, delay, dphase, power, b, rho in zip_longest(load_ids, excite_ids,
                                                                        delays, dphases, powers,
                                                                        bs, rhos):
            list_fields = ['ACSRCE', sid, excite_id, delay, dphase,
                           power, rho, b]
            bdf_file.write(print_card(list_fields))
        return

class RANDPS(VectorizedBaseCard):
    r"""
    Power Spectral Density Specification

    Defines load set power spectral density factors for use in random analysis
    having the frequency dependent form:

    .. math:: S_{jk}(F) = (X+iY)G(F)
    """
    _id_name = 'load_id'
    def clear(self) -> None:
        self.n = 0
        self.load_id = np.array([], dtype='int32')

    #def slice_card_by_index(self, i: np.ndarray) -> RANDPS:
        #load = RANDPS(self.model)
        #self.__apply_slice__(load, i)
        #return load

    def add_card(self, card: BDFCard, comment: str=''):
        """
        Adds a RANDPS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        j = integer(card, 2, 'j')
        k = integer(card, 3, 'k')
        fdouble_or_blank = force_double_or_blank if self.model.is_lax_parser else double_or_blank
        x = fdouble_or_blank(card, 4, 'x', default=0.0)
        y = fdouble_or_blank(card, 5, 'y', default=0.0)
        tabrnd1_id = integer_or_blank(card, 6, 'tid', default=0)
        assert len(card) <= 7, f'len(RANDPS card) = {len(card):d}\ncard={card}'
        self.cards.append((sid, j, k, x, y, tabrnd1_id, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)

        # sid : int
        #     random analysis set id
        #     defined by RANDOM in the case control deck
        # j : int
        #     Subcase id of the excited load set
        # k : int
        #     Subcase id of the applied load set
        #     k > j
        # x / y : float; default=0.0
        #     Components of the complex number
        # tid : int; default=0
        #     TABRNDi id that defines G(F)
        # comment : str; default=''
        #     a comment for the card

        #: Set identification number
        load_id = np.zeros(ncards, dtype='int32')
        subcase_excite_load = np.zeros(ncards, dtype='int32')
        subcase_applied_load = np.zeros(ncards, dtype='int32')
        xy = np.zeros((ncards, 2), dtype='float64')
        tabrnd1_id = np.zeros(ncards, dtype='int32')

        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (sid, subcase_excite_loadi, subcase_applied_loadi,
             x, y, tabrnd1_idi, comment) = card
            load_id[icard] = sid
            subcase_excite_load[icard] = subcase_excite_loadi
            subcase_applied_load[icard] = subcase_applied_loadi
            xy[icard] = [x, y]
            tabrnd1_id[icard] = tabrnd1_idi

        self._save(load_id, subcase_excite_load, subcase_applied_load, xy, tabrnd1_id)
        assert len(self.load_id) == self.n
        self.cards = []

    def _save(self, load_id,
              subcase_excite_load, subcase_applied_load,
              xy, tabrnd1_id) -> None:
        if len(self.load_id) != 0:
            sadf
        self.load_id = load_id
        self.subcase_excite_load = subcase_excite_load
        self.subcase_applied_load = subcase_applied_load
        self.xy = xy
        self.tabrnd1_id = tabrnd1_id

    def __apply_slice__(self, load: RANDPS, i: np.ndarray) -> None:
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.subcase_excite_load = self.subcase_excite_load[i]
        load.subcase_applied_load = self.subcase_applied_load[i]
        load.xy = self.xy[i, :]
        load.tabrnd1_id = self.tabrnd1_id[i]

    @property
    def x(self) -> np.ndarray:
        return self.xy[:, 0]
    @property
    def y(self) -> np.ndarray:
        return self.xy[:, 1]

    @x.setter
    def x(self, x: np.ndarray) -> None:
        self.xy[:, 0] = x
    @y.setter
    def y(self, y: np.ndarray) -> None:
        self.xy[:, 1] = y

    @property
    def max_id(self) -> int:
        return max(self.load_id.max(),
                   self.subcase_excite_load.max(),
                   self.subcase_applied_load.max(),
                   self.tabrnd1_id.max(),)

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        #array_str, array_default_int

        load_ids = array_str(self.load_id, size=size)
        subcase_excite_loads = array_str(self.subcase_excite_load, size=size)
        subcase_applied_loads = array_str(self.subcase_applied_load, size=size)
        tabrnd1_id = array_default_int(self.tabrnd1_id, default=0, size=size)
        for sid, subcase_excite_load, subcase_applied_load, \
            x, y, tabrnd1_id in zip_longest(load_ids, subcase_excite_loads,
                                            subcase_applied_loads,
                                            self.x, self.y, tabrnd1_id):
            list_fields = ['RANDPS', sid, subcase_excite_load, subcase_applied_load,
                           x, y, tabrnd1_id]
            bdf_file.write(print_card(list_fields))
        return

    #def sum_forces_moments(self) -> np.ndarray:
        #nloads = len(self.load_id)
        #force_moment = np.zeros((nloads, 6), dtype='float64')
        #return force_moment


def _set_int_float(i: int, array_int: np.ndarray, array_float: np.ndarray, value: Union[int, float]) -> None:
    if isinstance(value, int):
        array_int[i] = value
    else:
        array_float[i] = value

def _write_int_float(value_int: int, value_float: float) -> Union[int, float]:
    if value_int > 0:
        value = value_int
    else:
        value = value_float
    return value

def array_int_float(delay_int: np.ndarray,
                    delay_float: np.ndarray,
                    size: int=8, is_double: bool=False) -> np.ndarray:
    delay_ids = array_str(delay_int, size=size)
    delay_floats = array_float(delay_float, size=size, is_double=is_double)
    idelay_int = (delay_int > 0)
    delay_floats[idelay_int] = delay_ids[idelay_int]
    return delay_floats
