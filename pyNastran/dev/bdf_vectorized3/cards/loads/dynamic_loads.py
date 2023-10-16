from __future__ import annotations
#from collections import defaultdict
from itertools import zip_longest
from typing import Union, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.field_writer_8 import set_blank_if_default # , set_string8_blank_if_default
#from pyNastran.bdf.cards.base_card import expand_thru_by # expand_thru # BaseCard, expand_thru_by #  _node_ids,
#from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8, print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
#from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, # integer_or_blank,
    double_or_blank,
    components_or_blank, integer_string_or_blank, integer_double_or_blank,
    #integer_or_string, modal_components_or_blank,
)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.bdf.cards.loads.dloads import (
    fix_loadtype_tload1, fix_loadtype_tload2,
    #fix_loadtype_rload1, fix_loadtype_rload2,
)

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, # make_idim, hslice_by_idim,
    parse_load_check, # get_print_card_8_16,
)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_default_int, get_print_card_size)
#from pyNastran.dev.bdf_vectorized3.utils import cast_int_array
#from .static_loads import get_loads_by_load_id, get_reduced_loads

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


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
        return self.n

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
        return self.n

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
        return self.n

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
        return self.n

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

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        max_int = max(self.load_id.max(),
                      self.excite_id.max(),
                      self.delay_int.max(),
                      self.tabled_id.max())
        print_card, size = get_print_card_size(size, max_int)

        #array_str, array_default_int
        load_ids = array_str(self.load_id, size=size)
        excite_ids = array_default_int(self.excite_id, default=0, size=size)
        for sid, excite_id, delay_int, delay_float, \
            load_type, tabled_id, us0, vs0 in zip_longest(load_ids, excite_ids,
                                                          self.delay_int, self.delay_float,
                                                          self.load_type,
                                                          self.tabled_id, self.us0, self.vs0):
            us0 = set_blank_if_default(us0, 0.0)
            vs0 = set_blank_if_default(vs0, 0.0)
            delay = _write_int_float(delay_int, delay_float)
            list_fields = ['TLOAD1', sid, excite_id, delay, load_type,
                           tabled_id, us0, vs0]
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        udelay = np.unique(self.delay_int)

        geom_check(
            self,
            missing,
            #delay=(model.delay.delay_id, udelay),
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
        return self.n

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
        return self.n

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

    @parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        max_int = max(self.load_id.max(),
                      self.excite_id.max(),
                      self.delay_int.max())
        print_card, size = get_print_card_size(size, max_int)

        #array_str, array_default_int
        load_ids = array_str(self.load_id, size=size)
        excite_ids = array_default_int(self.excite_id, default=0, size=size)

        for sid, excite_id, delay_int, delay_float, \
            load_type, T, \
            frequency, phase, b, c, \
            us0, vs0 in zip_longest(load_ids, excite_ids, self.delay_int, self.delay_float,
                                    self.load_type, self.T,
                                    self.frequency, self.phase, self.b, self.c,
                                    self.us0, self.vs0):
            T1, T2 = T
            frequency = set_blank_if_default(frequency, 0.0)
            phase = set_blank_if_default(phase, 0.0)
            c = set_blank_if_default(c, 0.0)
            b = set_blank_if_default(b, 0.0)

            delay = _write_int_float(delay_int, delay_float)
            us0 = set_blank_if_default(us0, 0.0)
            vs0 = set_blank_if_default(vs0, 0.0)
            list_fields = ['TLOAD2', sid, excite_id, delay, load_type,
                           T1, T2, frequency, phase, c, b, us0, vs0]
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        udelay = np.unique(self.delay_int)

        geom_check(
            self,
            missing,
            #delay=(model.delay.delay_id, udelay),
        )

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
