from __future__ import annotations
from collections import defaultdict
from itertools import zip_longest
from typing import Union, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default # , set_string8_blank_if_default
#from pyNastran.bdf.cards.base_card import expand_thru_by # expand_thru # BaseCard, expand_thru_by #  _node_ids,
#from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8, print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
#from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double,
    integer_or_blank, double_or_blank,
    components_or_blank, integer_string_or_blank, integer_double_or_blank,
    #integer_or_string,
    modal_components_or_blank,
)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.bdf.cards.loads.dloads import (
    fix_loadtype_tload1, fix_loadtype_tload2,
    fix_loadtype_rload1, # fix_loadtype_rload2,
)

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, make_idim, # hslice_by_idim,
    parse_load_check, # get_print_card_8_16,
)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_float,
    array_default_int, array_default_float,
    get_print_card_size)
#from pyNastran.dev.bdf_vectorized3.utils import cast_int_array
#from .static_loads import get_loads_by_load_id, get_reduced_loads

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class DLOAD(VectorizedBaseCard):
    """
    +------+-----+------+------+----+-----+----+----+----+
    |   1  |  2  |  3   |  4   | 5  |  6  | 7  | 8  | 9  |
    +======+=====+======+======+====+=====+====+====+====+
    | LOAD | SID |  S   |  S1  | L1 | S2  | L2 | S3 | L3 |
    +------+-----+------+------+----+-----+----+----+----+
    |      | S4  |  L4  | etc. |    |     |    |    |    |
    +------+-----+------+------+----+-----+----+----+----+
    | LOAD | 101 | -0.5 | 1.0  | 3  | 6.2 | 4  |    |    |
    +------+-----+------+------+----+-----+----+----+----+

    """
    def clear(self) -> None:
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
        if isinstance(scale_factors, float) and isinstance(load_ids, int):
            scale_factors = [scale_factors]
            load_ids = [load_ids]
        elif isinstance(scale_factors, float):
            scale_factors = [scale_factors] * len(load_ids)
        else:
            load_ids = [load_ids] * len(scale_factors)

            scale_factors
        self.cards.append((sid, scale, scale_factors, load_ids, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> None:
        sid = integer(card, 1, 'sid')
        scale = double(card, 2, 'scale')

        scale_factors = []
        load_ids = []

        # alternating of scale factor & load set ID
        nload_fields = len(card) - 3
        assert nload_fields % 2 == 0, 'card=%s' % card
        for iload in range(nload_fields // 2):
            n = 2 * iload + 3
            scale_factors.append(double(card, n, 'scale_factor'))
            load_ids.append(integer(card, n + 1, 'load_id'))

        assert len(card) > 3, 'len(%s card) = %i\ncard=%s' % (self.type, len(card), card)
        self.cards.append((sid, scale, scale_factors, load_ids, comment))
        self.n += 1
        return self.n - 1

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        nloads = len(self.cards)
        if nloads == 0:
            return
        self.load_id = np.zeros(nloads, dtype='int32')
        self.scale = np.zeros(nloads, dtype='float64')
        self.nloads = np.zeros(nloads, dtype='int32')

        all_load_ids = []
        all_scale_factors = []
        assert nloads > 0, nloads
        for icard, card in enumerate(self.cards):
            (sid, scale, scale_factors, load_ids, comment) = card
            nloads_actual = len(scale_factors)

            self.load_id[icard] = sid
            self.scale[icard] = scale
            self.nloads[icard] = nloads_actual
            all_load_ids.extend(load_ids)
            all_scale_factors.extend(scale_factors)
        self.load_ids = np.array(all_load_ids, dtype='int32')
        self.scale_factors = np.array(all_scale_factors, dtype='float64')
        self.model.log.debug('load_id=%s scale=%s nloads=%s all_load_ids=%s all_scale_factors=%s' % (
            self.load_id, self.scale, self.nloads, all_load_ids, all_scale_factors))
        self.cards = []

    @property
    def iload(self) -> np.ndarray:
        return make_idim(self.n, self.nloads)

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        max_int = max(self.load_id.max(), self.load_ids.max())
        print_card, size = get_print_card_size(size, max_int)

        for sid, scale, iloadi in zip(self.load_id, self.scale, self.iload):
            iload0, iload1 = iloadi
            list_fields = ['DLOAD', sid, scale]
            scale_factors = self.scale_factors[iload0:iload1]
            load_ids = self.load_ids[iload0:iload1]
            for (scale_factor, load_id) in zip_longest(scale_factors, load_ids):
                list_fields += [scale_factor, load_id]
            #if len(load_ids) != len(scale_factors):
                #msg = 'nload_ids=%s nscale_factors=%s and arent the same\n' % (
                    #len(load_ids), len(scale_factors))
                #msg = 'load_ids=%s\n' % (load_ids)
                #msg += 'scale_factors=%s\n' % (scale_factors)
                #msg += print_card_8(list_fields)
                #msg += str(self.get_stats())
                #raise IndexError(msg)
            bdf_file.write(print_card(list_fields))
        return

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
    def __init__(self, model: BDF):
        super().__init__(model)
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

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        max_int = max(self.load_id.max(), self.node_id.max())
        print_card, size = get_print_card_size(size, max_int)

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
    def __init__(self, model: BDF):
        super().__init__(model)

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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.load_id = np.array([], dtype='int32')

    def slice_card_by_index(self, i: np.ndarray) -> TLOAD2:
        load = TLOAD2(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: TLOAD2, i: np.ndarray) -> None:
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.excite_id = self.excite_id[i]

        load.delay_int = self.delay_int[i]
        load.delay_float = self.delay_float[i]
        load.load_type = self.load_type[i]
        load.tabled_id = self.tabled_id[i]
        load.us0 = self.us0[i]
        load.vs0 = self.vs0[i]

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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.load_id = np.array([], dtype='int32')

    def slice_card_by_index(self, i: np.ndarray) -> RLOAD1:
        load = RLOAD1(self.model)
        self.__apply_slice__(load, i)
        return load

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
        tabled_d_int = np.full(ncards, np.nan, dtype='float64')

        tabled_c_float = np.zeros(ncards, dtype='int32')
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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.load_id = np.array([], dtype='int32')

    def slice_card_by_index(self, i: np.ndarray) -> RLOAD2:
        load = RLOAD2(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: RLOAD2, i: np.ndarray) -> None:
        self.model.log.info(self.dphase_int)
        self.model.log.info(i)
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

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.lseq_id = np.array([], dtype='int32')
        self.excite_id = np.array([], dtype='int32')
        self.load_id = np.array([], dtype='int32')
        self.temperature_id = np.array([], dtype='int32')

    def slice_card_by_index(self, i: np.ndarray) -> LSEQ:
        load = LSEQ(self.model)
        self.__apply_slice__(load, i)
        return load

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

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
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
    #def __init__(self, model: BDF):
        #super().__init__(model)
        #self._is_sorted = False
    def clear(self) -> None:
        self.tic_id = np.array([], dtype='int32')
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
        tic._is_sorted = self._is_sorted
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
    def parse_cards(self):
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug('parse TIC')
        tic_id = np.zeros(ncards, dtype='int32')
        node_id = np.zeros(ncards, dtype='int32')
        component = np.zeros(ncards, dtype='int32')
        u0 = np.zeros(ncards, dtype='float64')
        v0 = np.zeros(ncards, dtype='float64')

        for i, card in enumerate(self.cards):
            (sid, nid, comp, u0i, v0i, comment) = card
            tic_id[i] = sid
            node_id[i] = nid
            component[i] = comp
            u0[i] = u0i
            v0[i] = v0i
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
    def clear(self) -> None:
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

    def slice_card_by_index(self, i: np.ndarray) -> DELAY:
        """uses a node_index to extract DELAYs"""
        assert len(self.node_id) > 0, self.node_id
        i = np.atleast_1d(np.asarray(i, dtype=self.node_id.dtype))
        i.sort()
        delay = DELAY(self.model)
        self.__apply_slice__(delay, i)
        return delay

    def __apply_slice__(self, delay: DELAY, i: np.ndarray) -> None:
        delay.n = len(i)
        delay._is_sorted = self._is_sorted
        delay.node_id = self.node_id[i]
        delay.delay_id = self.delay_id[i]
        delay.delay = self.delay[i]
        delay.component = self.component[i]

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
    def parse_cards(self):
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug('parse DELAY')
        delay_id = np.zeros(ncards, dtype='int32')
        node_id = np.zeros(ncards, dtype='int32')
        component = np.zeros(ncards, dtype='int32')
        delay = np.zeros(ncards, dtype='float64')

        for i, card in enumerate(self.cards):
            (delay_idi, nodei, componenti, delayi, comment) = card
            delay_id[i] = delay_idi
            node_id[i] = nodei
            component[i] = componenti
            delay[i] = delayi
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

    def clear(self) -> None:
        self.dphase_id = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')

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

    def slice_card_by_index(self, i: np.ndarray) -> DPHASE:
        """uses a node_index to extract DELAYs"""
        assert len(self.node_id) > 0, self.node_id
        i = np.atleast_1d(np.asarray(i, dtype=self.node_id.dtype))
        i.sort()
        dphase = DPHASE(self.model)
        self.__apply_slice__(dphase, i)
        return dphase

    def __apply_slice__(self, dphase: DPHASE, i: np.ndarray) -> None:
        dphase.n = len(i)
        dphase._is_sorted = self._is_sorted
        dphase.dphase_id = self.dphase_id[i]
        dphase.node_id = self.node_id[i]
        dphase.component = self.component[i]
        dphase.phase_lead = self.phase_lead[i]

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
    def parse_cards(self):
        ncards = len(self.cards)
        if self.debug:
            self.model.log.debug('parse DPHASE')
        dphase_id = np.zeros(ncards, dtype='int32')
        node_id = np.zeros(ncards, dtype='int32')
        component = np.zeros(ncards, dtype='int32')
        phase_lead = np.zeros(ncards, dtype='float64')

        for i, card in enumerate(self.cards):
            (dphase_idi, nodei, componenti, phase_leadi, comment) = card
            dphase_id[i] = dphase_idi
            node_id[i] = nodei
            component[i] = componenti
            phase_lead[i] = phase_leadi
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

def array_int_float(delay_int: np.ndarray, delay_float: np.ndarray,
                    size: int=8, is_double: bool=False) -> np.ndarray:
    delay_ids = array_str(delay_int, size=size)
    delay_floats = array_float(delay_float, size=size, is_double=is_double)
    idelay_int = (delay_int > 0)
    delay_floats[idelay_int] = delay_ids[idelay_int]
    return delay_floats
