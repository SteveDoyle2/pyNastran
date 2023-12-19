from __future__ import annotations
from itertools import zip_longest
from typing import Union, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types, float_types
#from pyNastran.bdf import MAX_INT
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    hslice_by_idim, remove_unused_primary) #, BaseCard, _node_ids, expand_thru
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, integer_or_string,
    string, string_or_blank, integer_double_string_or_blank,
    integer_string_or_blank, integer_double_or_blank,
    interpret_value)
from pyNastran.bdf.bdf_interface.assign_type_force import force_double, force_double_or_blank

from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.utils import build_table_lines

from pyNastran.bdf.cards.optimization import build_table_lines, parse_table_fields, DRESP2_PACK_LENGTH
from pyNastran.bdf.cards.base_card import expand_thru_by

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, make_idim, get_print_card_8_16,
    #hslice_by_idim,
)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_float, array_default_int, array_default_float, array_default_str)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.constraints import ADD

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

class DCONADD(ADD):
    _id_name = 'dconadd_id'
    @property
    def dconadd_id(self):
        return self.sid
    @dconadd_id.setter
    def dconadd_id(self, dconadd_id: np.ndarray):
        self.sid = dconadd_id

    @property
    def dconstr_ids(self):
        return self.sids
    @dconstr_ids.setter
    def dconstr_ids(self, dconstr_ids: np.ndarray):
        self.sids = dconstr_ids

    @property
    def ndconstrs(self):
        return self.nsids
    @ndconstrs.setter
    def ndconstrs(self):
        return self.nsids

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['dconstr_id'].append(self.dconstr_ids)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dconadd_id) == 0:
            return
        if size == 8 and self.is_small_field:
            print_card = print_card_8
        else:
            print_card = print_card_16

        #self.get_reduced_spcs()
        dconadd_ids = array_str(self.dconadd_id, size=size)
        dconstr_ids = array_str(self.dconstr_ids, size=size)
        for dconadd_id, idim in zip(dconadd_ids, self.idim):
            idim0, idim1 = idim
            dconstr_idsi = dconstr_ids[idim0:idim1].tolist()
            list_fields = ['DCONADD', dconadd_id] + dconstr_idsi
            bdf_file.write(print_card(list_fields))
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        dconstr_id = np.unique(self.model.dconstr.dconstr_id)
        geom_check(self,
                   missing,
                   dconstr_id=(dconstr_id, self.dconstr_ids),
                   )


class MODTRAK(VectorizedBaseCard):
    """
    MODTRAK SID LOWRNG HIGHRNG MTFILTER
    MODTRAK 100   1      26      0.80
    """
    _id_name = 'modtrak_id'
    def clear(self) -> None:
        self.n = 0
        self.modtrak_id = np.array([], dtype='int32')
        self.low_range = np.array([], dtype='int32')
        self.high_range = np.array([], dtype='int32')
        self.mt_filter = np.array([], dtype='float64')

    def add(self, modtrak_id: int, low_range: int, high_range: int,
            mt_filter: float, comment: str='') -> int:
        self.cards.append((modtrak_id, low_range, high_range, mt_filter, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        modtrak_id = integer(card, 1, 'sid')
        low_range = integer_or_blank(card, 2, 'low_range', default=0)
        high_range = integer(card, 3, 'high_range')
        mt_filter = double_or_blank(card, 4, 'mt_filter', default=0.9)
        #return MODTRAK(sid, low_range, high_range, mt_filter, comment=comment)
        assert len(card) <= 5, f'len(MODTRAK card) = {len(card):d}\ncard={card}'
        self.cards.append((modtrak_id, low_range, high_range, mt_filter, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        modtrak_id = np.zeros(ncards, dtype='int32')
        low_range = np.zeros(ncards, dtype='int32')
        high_range = np.zeros(ncards, dtype='int32')
        mt_filter = np.zeros(ncards, dtype='float64')
        for icard, card in enumerate(self.cards):
            modtrak_idi, low_rangei, high_rangei, mt_filteri, comment = card
            modtrak_id[icard] = modtrak_idi
            low_range[icard] = low_rangei
            high_range[icard] = high_rangei
            mt_filter[icard] = mt_filteri
        self._save(modtrak_id, low_range, high_range, mt_filter)
        self.cards = []

    def _save(self, modtrak_id, low_range, high_range, mt_filter):
        if len(self.modtrak_id) != 0:
            asfd
            modtrak_id = np.hstack([self.modtrak_id, modtrak_id])
        self.modtrak_id = modtrak_id
        self.low_range = low_range
        self.high_range = high_range
        self.mt_filter = mt_filter

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    #def geom_check(self, missing: dict[str, np.ndarray]):
        #pass

    #def index(self, desvar_id: np.ndarray) -> np.ndarray:
        #assert len(self.desvar_id) > 0, self.desvar_id
        #desvar_id = np.atleast_1d(np.asarray(desvar_id, dtype=self.desvar_id.dtype))
        #idesvar = np.searchsorted(self.desvar_id, desvar_id)
        #return idesvar

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.modtrak_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        modtrak_ids = array_str(self.modtrak_id, size=size)
        low_ranges = array_str(self.low_range, size=size)
        high_ranges = array_str(self.high_range, size=size)
        mt_filters = array_float(self.mt_filter, size=size, is_double=False)

        for modtrak_id, low_range, high_range, mt_filter in zip_longest(modtrak_ids, low_ranges, high_ranges, mt_filters):
            list_fields = ['MODTRAK', modtrak_id, low_range, high_range, mt_filter]
            bdf_file.write(print_card(list_fields))
        return


class DESVAR(VectorizedBaseCard):
    """
    +--------+-----+-------+-------+-----+-----+-------+-------+
    |    1   |  2  |    3  |   4   |  5  |  6  |    7  |   8   |
    +========+=====+=======+=======+=====+=====+=======+=======+
    | DESVAR | OID | LABEL | XINIT | XLB | XUB | DELXV | DDVAL |
    +--------+-----+-------+-------+-----+-----+-------+-------+
    """
    _id_name = 'desvar_id'
    def clear(self) -> None:
        self.n = 0
        self.desvar_id = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U8')
        self.xinit = np.array([], dtype='float64')
        self.xlb = np.array([], dtype='float64')
        self.xub = np.array([], dtype='float64')
        self.delx = np.array([], dtype='float64')
        self.ddval = np.array([], dtype='int32')

    def add(self, desvar_id: int, label: str, xinit: float,
            xlb: float=-1e20, xub: float=1e20,
            delx=None, ddval: Optional[int]=None,
            comment: str='') -> int:
        """
        Creates a DESVAR card

        Parameters
        ----------
        desvar_id : int
            design variable id
        label : str
            name of the design variable
        xinit : float
            the starting point value for the variable
        xlb : float; default=-1.e20
            the lower bound
        xub : float; default=1.e20
            the lower bound
        delx : float; default=1.e20
            fractional change allowed for design variables during
            approximate optimization
            NX  if blank : take from DOPTPRM; otherwise 1.0
            MSC if blank : take from DOPTPRM; otherwise 0.5
        ddval : int; default=None
            int : DDVAL id
                  allows you to set discrete values
            None : continuous
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((desvar_id, label, xinit, xlb, xub, delx, ddval, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        fdouble = force_double if self.model.is_lax_parser else double
        fdouble_or_blank = force_double_or_blank if self.model.is_lax_parser else double_or_blank

        desvar_id = integer(card, 1, 'desvar_id')
        label = string(card, 2, 'label')
        xinit = fdouble(card, 3, 'xinit')
        xlb = fdouble_or_blank(card, 4, 'xlb', -1e20)
        xub = fdouble_or_blank(card, 5, 'xub', 1e20)
        delx = fdouble_or_blank(card, 6, 'delx', default=np.nan)
        ddval = integer_or_blank(card, 7, 'ddval', default=0)
        assert len(card) <= 8, f'len(DESVAR card) = {len(card):d}\ncard={card}'
        self.cards.append((desvar_id, label, xinit, xlb, xub, delx, ddval, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        desvar_id = np.zeros(ncards, dtype='int32')
        #: user-defined name for printing purposes
        label = np.zeros(ncards, dtype='|U8')
        xinit = np.zeros(ncards, dtype='float64')
        xlb = np.zeros(ncards, dtype='float64')
        xub = np.zeros(ncards, dtype='float64')
        delx = np.zeros(ncards, dtype='float64')
        ddval = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            desvar_idi, labeli, xiniti, xlbi, xubi, delxi, ddvali, comment = card
            desvar_id[icard] = desvar_idi
            label[icard] = labeli
            xinit[icard] = xiniti
            xlb[icard] = xlbi
            xub[icard] = xubi
            delx[icard] = delxi
            if ddvali is None:
                ddvali = 0
            assert isinstance(ddvali, int), ddvali
            ddval[icard] = ddvali

        if self.model.apply_clip_to_desvar_range:
            i = (xinit < xlb)
            xinit[i] = xlb[i]

            i = (xinit > xub)
            xinit[i] = xub[i]

        #self.label = label
        ##xinit = np.clip(xinit, xlb, xub)
        #self.xinit = xinit
        #self.xlb = xlb
        #self.xub = xub
        #assert len(label) <= 8, f'desvar_id={desvar_id} label={label!r} must be less than 8 characters; length={len(label):d}'
        #assert xlb <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit >= xlb, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        ## controls change for a single optimization cycle
        ## taken from DOPTPRM if None; else default=1.
        #self.delx = delx
        ## DDVAL id if you want discrete values
        #self.ddval = ddval
        #self.ddval_ref = None
        ##assert ' ' not in label.rstrip(), self.get_stats()

        self._save(desvar_id, label, xinit, xlb, xub, delx, ddval)
        self.cards = []

    def _save(self, desvar_id, label, xinit, xlb, xub, delx, ddval):
        if len(self.desvar_id) != 0:
            desvar_id = np.hstack([self.desvar_id, desvar_id])
            label = np.hstack([self.label, label])
            xinit = np.hstack([self.xinit, xinit])
            xlb = np.hstack([self.xlb, xlb])
            xub = np.hstack([self.xub, xub])
            delx = np.hstack([self.delx, delx])
            ddval = np.hstack([self.ddval, ddval])
        self.desvar_id = desvar_id
        self.label = label
        self.xinit = xinit
        self.xlb = xlb
        self.xub = xub
        self.delx = delx
        self.ddval = ddval

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['ddval_id'].append(self.ddval[self.ddval > 0])

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        desvar_id = used_dict['desvar_id']
        ncards_removed = remove_unused_primary(self, desvar_id, self.desvar_id, 'desvar_id')
        return ncards_removed

    def index(self, desvar_id: np.ndarray) -> np.ndarray:
        assert len(self.desvar_id) > 0, self.desvar_id
        desvar_id = np.atleast_1d(np.asarray(desvar_id, dtype=self.desvar_id.dtype))
        idesvar = np.searchsorted(self.desvar_id, desvar_id)
        return idesvar

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.desvar_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        is_delx = not np.all(np.isnan(self.delx))
        is_ddval = not np.all(np.isnan(self.ddval))

        desvar_id = array_str(self.desvar_id, size=size)
        labels = self.label
        #xlbs = array_default_8(self.xlb, -1e20)
        #xubs = array_default_8(self.xub, 1e20)
        xlbs = array_float(self.xlb, size=size, is_double=is_double)
        xubs = array_float(self.xub, size=size, is_double=is_double)

        if is_delx and is_ddval:
            delxs = array_float(self.delx, size=size, is_double=is_double)
            for desvar_id, label, xinit, xlb, xub, delx, ddval in zip_longest(desvar_id, labels, self.xinit, xlbs, xubs,
                                                                              delxs, self.ddval):
                list_fields = ['DESVAR', desvar_id, label, xinit, xlb, xub,
                               delx, ddval]
                bdf_file.write(print_card(list_fields))
        elif is_delx:
            for desvar_id, label, xinit, xlb, xub, delx in zip_longest(desvar_id, labels, self.xinit, xlbs, xubs, delxs):
                list_fields = ['DESVAR', desvar_id, label, xinit, xlb, xub,
                               delx]
                bdf_file.write(print_card(list_fields))
        else:
            # is_ddval
            for desvar_id, label, xinit, xlb, xub, ddval in zip_longest(desvar_id, labels, self.xinit, xlbs, xubs, self.ddval):
                list_fields = ['DESVAR', desvar_id, label, xinit, xlb, xub,
                               None, ddval]
                bdf_file.write(print_card(list_fields))
        return


class DDVAL(VectorizedBaseCard):
    """
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    |   1   |  2  |   3   |   4   |   5   |   6   |   7   |   8   |   9   |
    +=======+=====+=======+=======+=======+=======+=======+=======+=======+
    | DDVAL | ID  | DVAL1 | DVAL2 | DVAL3 | DVAL4 | DVAL5 | DVAL6 | DVAL7 |
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    | DDVAL | ID  | DVAL1 | THRU  | DVAL2 | BY    |  INC  |       |       |
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+

    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    | DDVAL | 110 |  0.1  |  0.2  |  0.3  |  0.5  |  0.6  |  0.4  |       |
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    |       | .7  | THRU  |  1.0  |  BY   | 0.05  |       |       |       |
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    |       | 1.5 |  2.0  |       |       |       |       |       |       |
    +-------+-----+-------+-------+-------+-------+-------+-------+-------+
    """
    _id_name = 'ddval_id'
    def clear(self) -> None:
        self.n = 0
        self.ddval_id = np.array([], dtype='int32')
        self.nddval = np.array([], dtype='int32')
        self.value = np.array([], dtype='float64')

    def add(self, desvar_id: int, label: str, xinit: float,
            xlb: float=-1e20, xub: float=1e20,
            delx=None, ddval: Optional[int]=None,
            comment: str='') -> int:
        """
        Creates a DESVAR card

        Parameters
        ----------
        desvar_id : int
            design variable id
        label : str
            name of the design variable
        xinit : float
            the starting point value for the variable
        xlb : float; default=-1.e20
            the lower bound
        xub : float; default=1.e20
            the lower bound
        delx : float; default=1.e20
            fractional change allowed for design variables during
            approximate optimization
            NX  if blank : take from DOPTPRM; otherwise 1.0
            MSC if blank : take from DOPTPRM; otherwise 0.5
        ddval : int; default=None
            int : DDVAL id
                  allows you to set discrete values
            None : continuous
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((desvar_id, label, xinit, xlb, xub, delx, ddval, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        fdouble = force_double if self.model.is_lax_parser else double
        fdouble_or_blank = force_double_or_blank if self.model.is_lax_parser else double_or_blank

        """
        Adds a DDVAL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        oid = integer(card, 1, 'oid')
        n = 1
        ddvals = []
        for i in range(2, len(card)):
            ddval = integer_double_string_or_blank(card, i, 'DDVAL%s' % n)
            if ddval is not None:
                ddvals.append(ddval)
        #return DDVAL(oid, ddvals, comment=comment)
        self.cards.append((oid, ddvals, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ddval_id = np.zeros(ncards, dtype='int32')
        nddval = np.zeros(ncards, dtype='int32')
        all_values_list = []
        for icard, card in enumerate(self.cards):
            ddval_idi, values, comment = card
            ddval_id[icard] = ddval_idi

            if not self.model.is_lax_parser:
                for ddval in values:
                    assert not isinstance(ddval, integer_types), f'DDVALs id={ddval_idi} have integer fields={values}'
            values2 = expand_thru_by(values, require_int=False)
            values2.sort()
            all_values_list.extend(values2)
            nddval[icard] = len(values2)

        all_values = np.array(all_values_list, dtype='float64')
        self._save(ddval_id, nddval, all_values)
        self.sort()
        self.cards = []

    def _save(self, ddval_id, nddval, value):
        if len(self.ddval_id) != 0:
            ddval_id = np.hstack([self.ddval_id, ddval_id])
            nddval = np.hstack([self.nddval, nddval])
            value = np.hstack([self.value, value])
        self.ddval_id = ddval_id
        self.nddval = nddval
        self.value = value

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    #def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        #used_dict['ddval_id'].append(self.ddval[self.ddval > 0])

    #def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        #desvar_id = used_dict['desvar_id']
        #ncards_removed = remove_unused_primary(self, desvar_id, self.desvar_id, 'desvar_id')
        #return ncards_removed

    #def index(self, desvar_id: np.ndarray) -> np.ndarray:
        #assert len(self.desvar_id) > 0, self.desvar_id
        #desvar_id = np.atleast_1d(np.asarray(desvar_id, dtype=self.desvar_id.dtype))
        #idesvar = np.searchsorted(self.desvar_id, desvar_id)
        #return idesvar

    @property
    def iddval(self) -> np.ndarray:
        return make_idim(self.n, self.nddval)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.ddval_id) == 0:
            return
        print_card = get_print_card_8_16(size)

        ddval_ids = array_str(self.ddval_id, size=size)
        values = array_float(self.value, size=size, is_double=is_double).tolist()
        for ddval_id, (iddval0, iddval1) in zip_longest(ddval_ids, self.iddval):
            value = values[iddval0:iddval1]
            list_fields = ['DDVAL', ddval_id] + value
            bdf_file.write(print_card(list_fields))
        return


class DLINK(VectorizedBaseCard):
    """
    Multiple Design Variable Linking
    Relates one design variable to one or more other design variables.

    +-------+------+-------+--------+-------+------+----+------+----+
    |   1   |   2  |   3   |   4    |   5   |   6  |  7 |   8  | 9  |
    +=======+======+=======+========+=======+======+====+======+====+
    | DLINK |  ID  | DDVID |   C0   | CMULT | IDV1 | C1 | IDV2 | C2 |
    +-------+------+-------+--------+-------+------+----+------+----+
    |       | IDV3 |   C3  |  etc.  |       |      |    |      |    |
    +-------+------+-------+--------+-------+------+----+------+----+
    """
    _id_name = 'dlink_id'
    def clear(self) -> None:
        self.n = 0
        self.dlink_id = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U8')
        self.xinit = np.array([], dtype='float64')
        self.xlb = np.array([], dtype='float64')
        self.xub = np.array([], dtype='float64')
        self.delx = np.array([], dtype='float64')
        self.ddval = np.array([], dtype='float64')

    def add(self, dlink_id: int, dependent_desvar: int,
            independent_desvars: list[int],
            coeffs: list[float],
            c0: float=0.0, cmult: float=1.0, comment: str='') -> int:
        """
        Creates a DLINK card, which creates a variable that is a lienar
        ccombination of other design variables

        Parameters
        ----------
        dlink_id : int
            optimization id
        dependent_desvar : int
            the DESVAR to link
        independent_desvars : list[int]
            the DESVARs to combine
        coeffs : list[int]
            the linear combination coefficients
        c0 : float; default=0.0
            an offset
        cmult : float; default=1.0
            an scale factor
        comment : str; default=''
            a comment for the card

        """
        if isinstance(coeffs, float_types):
            coeffs = [coeffs]
        if isinstance(independent_desvars, integer_types):
            independent_desvars = [independent_desvars]
        self.cards.append((dlink_id, dependent_desvar, independent_desvars, coeffs,
                           c0, cmult, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DLINK card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        dlink_id = integer(card, 1, 'oid')
        dependent_desvar = integer(card, 2, 'dependent_desvar')
        c0 = double_or_blank(card, 3, 'c0', 0.)
        cmult = double_or_blank(card, 4, 'cmult', 1.)

        nfields = len(card) - 4
        n = nfields // 2
        independent_desvars = []
        coeffs = []

        for i in range(n):
            j = 2 * i + 5
            desvar = integer(card, j, 'independent_desvar_%d' % i)
            coeff = double(card, j + 1, 'coeff_%d' % i)
            independent_desvars.append(desvar)
            coeffs.append(coeff)
        #return DLINK(oid, dependent_desvar, independent_desvars, coeffs,
                     #c0=c0, cmult=cmult, comment=comment)
        self.cards.append((dlink_id, dependent_desvar, independent_desvars, coeffs,
                           c0, cmult, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        dlink_id = np.zeros(ncards, dtype='int32')
        dependent_desvar = np.zeros(ncards, dtype='int32')
        nindependent_desvars = np.zeros(ncards, dtype='int32')

        c0 = np.zeros(ncards, dtype='float64')
        cmult = np.zeros(ncards, dtype='float64')
        all_coefficients = []
        all_independent_desvars = []
        for icard, card in enumerate(self.cards):
            (dlink_idi, dependent_desvari, independent_desvarsi, coeffsi,
             c0i, cmulti, comment) = card
            dlink_id[icard] = dlink_idi
            dependent_desvar[icard] = dependent_desvari

            all_coefficients.extend(coeffsi)
            all_independent_desvars.extend(independent_desvarsi)
            nindependent_desvars[icard] = len(independent_desvarsi)
            #self.coefficents[icard] = coeffs
            c0[icard] = c0i
            cmult[icard] = cmulti

        coefficients = np.array(all_coefficients, dtype='float64')

        independent_desvars = np.array(all_independent_desvars, dtype='int32')
        self._save(dlink_id, dependent_desvar, c0, cmult,
                   nindependent_desvars, independent_desvars, coefficients)

    def _save(self, dlink_id, dependent_desvar,
              c0, cmult,
              nindependent_desvars, independent_desvars, coefficients):
        self.dlink_id = dlink_id
        self.dependent_desvar = dependent_desvar
        self.c0 = c0
        self.cmult = cmult

        self.nindependent_desvars = nindependent_desvars
        self.independent_desvars = independent_desvars
        self.coefficients = coefficients

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['desvar_id'].append(self.dependent_desvar)
        used_dict['desvar_id'].append(self.independent_desvars)

    @property
    def idesvar(self) -> np.ndarray:
        return make_idim(self.n, self.nindependent_desvars)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dlink_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        dlink_ids = array_str(self.dlink_id, size=size)
        dependent_desvar = array_str(self.dependent_desvar, size=size)
        c0s = array_default_float(self.c0, default=0.0, size=size, is_double=False)
        cmults = array_default_float(self.c0, default=1.0, size=size, is_double=False)
        for dlink_id, dep_desvar, (idesvar0, idesvar1), xinit, c0, cmult in zip_longest(
            dlink_ids, dependent_desvar, self.idesvar,
            self.xinit, c0s, cmults):
            #c0 = set_blank_if_default(c0, 0.)
            #cmult = set_blank_if_default(cmult, 1.)
            independent_desvars = self.independent_desvars[idesvar0 : idesvar1]
            coeffs = self.coefficients[idesvar0 : idesvar1]

            list_fields = ['DLINK', dlink_id, dep_desvar, c0, cmult]
            for (idv, ci) in zip(independent_desvars, coeffs):
                list_fields += [idv, ci]
            bdf_file.write(print_card(list_fields))
        return


class DVGRID(VectorizedBaseCard):
    """
    +--------+------+-----+-----+-------+----+----+----+
    |    1   |   2  |  3  |  4  |   5   |  6 |  7 |  8 |
    +========+======+=====+=====+=======+====+====+====+
    | DVGRID | DVID | GID | CID | COEFF | N1 | N2 | N3 |
    +--------+------+-----+-----+-------+----+----+----+
    """
    _id_name = 'desvar_id'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.n = 0
        self.desvar_id = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')
        self.coord_id = np.array([], dtype='int32')
        self.coefficient = np.array([], dtype='float64')
        self.dxyz = np.zeros((0, 3), dtype='float64')

    def add(self, desvar_id: int, nid: int, dxyz,
            cid: int=0, coeff: float=1.0, comment: str='') -> int:
        """
        Creates a DVGRID card

        Parameters
        ----------
        dvid : int
            DESVAR id
        nid : int
            GRID/POINT id
        dxyz : (3, ) float ndarray
            the amount to move the grid point
        cid : int; default=0
            Coordinate system for dxyz
        coeff : float; default=1.0
            the dxyz scale factor
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((desvar_id, nid, cid, coeff, dxyz, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DVGRID card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        desvar_id = integer(card, 1, 'dvid')
        nid = integer(card, 2, 'nid')
        cid = integer_or_blank(card, 3, 'cid', default=0)
        coeff = double_or_blank(card, 4, 'coeff', default=1.0)
        dxyz = [
            double_or_blank(card, 5, 'n1', default=0.),
            double_or_blank(card, 6, 'n2', default=0.),
            double_or_blank(card, 7, 'n3', default=0.),
        ]
        self.cards.append((desvar_id, nid, cid, coeff, dxyz, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        desvar_id = np.zeros(ncards, dtype='int32')
        node_id = np.zeros(ncards, dtype=idtype)
        coord_id = np.zeros(ncards, dtype='int32')
        coefficient = np.zeros(ncards, dtype='float64')
        dxyz = np.zeros((ncards, 3), dtype='float64')
        for icard, card in enumerate(self.cards):
            desvar_idi, nid, cid, coeff, dxyzi, comment = card
            desvar_id[icard] = desvar_idi
            node_id[icard] = nid
            coord_id[icard] = cid
            coefficient[icard] = coeff
            dxyz[icard, :] = dxyzi
        self._save(desvar_id, node_id, coord_id, coefficient, dxyz)
        self.cards = []

    def _save(self, desvar_id, node_id, coord_id,
              coefficient, dxyz) -> None:
        if len(self.desvar_id) != 0:
            asdf
        self.desvar_id = desvar_id
        self.node_id = node_id
        self.coord_id = coord_id
        self.coefficient = coefficient
        self.dxyz = dxyz

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['node_id'].append(self.node_id)
        used_dict['coord_id'].append(self.coord_id)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.node_id
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['node_id'].append(self.node_id)
        used_dict['coord_id'].append(self.coord_id)

    def convert(self, xyz_scale: float=1.0, **kwargs) -> None:
        """TODO: check the coordinate system..."""
        self.dxyz *= xyz_scale

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        cid = self.model.coord.coord_id
        geom_check(self,
                   missing,
                   node=(nid, self.node_id),
                   coord=(cid, self.coord_id))

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.desvar_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        desvar_ids = array_str(self.desvar_id, size=size)
        node_ids = array_str(self.node_id, size=size)
        coord_ids = array_default_int(self.coord_id, default=0, size=size)
        coefficients = array_float(self.coefficient, size=size, is_double=False)
        dxyzs = array_float(self.dxyz, size=size, is_double=False).tolist()
        for desvar_id, node_id, coord_id, coeff, dxyz in zip_longest(
            desvar_ids, node_ids, coord_ids, coefficients, dxyzs):
            list_fields = [
                'DVGRID', desvar_id, node_id, coord_id, coeff] + dxyz
            bdf_file.write(print_card(list_fields))
        return


FLOAT_RESPONSE_TYPES = {'DIVERG'}
class DRESP1(VectorizedBaseCard):
    """
    +--------+-------+---------+---------+--------+--------+-------+------+-------+
    |   1    |  2    |    3    |    4    |   5    |   6    |   7   |   8  |   9   |
    +========+=======+=========+=========+========+========+=======+======+=======+
    | DRESP1 |  OID  | LABEL   | RTYPE   | PTYPE  | REGION | ATTA  | ATTB | ATTI  |
    +--------+-------+---------+---------+--------+--------+-------+------+-------+
    | DRESP1 |  103  | STRESS2 |  STRESS | PSHELL |        |   9   |      |   3   |
    +--------+-------+---------+---------+--------+--------+-------+------+-------+
    | DRESP1 |  1S1  | CSTRAN3 | CSTRAIN |  PCOMP |        |   1   |  1   | 10000 |
    +--------+-------+---------+---------+--------+--------+-------+------+-------+
    """
    _id_name = 'dresp_id'
    def clear(self) -> None:
        self.n = 0
        self.dresp_id = np.array([], dtype='int32')

    def add(self, dresp_id: int, label: str,
            response_type: str, property_type: str, region: str,
            atta: Union[int, float, str, None],
            attb: Union[int, float, str, None],
            atti: list[Union[int, float, str]],
            validate: bool=True, comment: str='') -> int:
        """
        Creates a DRESP1 card.

        A DRESP1 is used to define a "simple" output result that may be
        optimized on.  A simple result is a result like stress, strain,
        force, displacement, eigenvalue, etc. for a node/element that
        may be found in a non-optimization case.

        Parameters
        ----------
        dresp_id : int
            response id
        label : str
            Name of the response
        response_type : str
            Response type
        property_type : str
            Element flag (PTYPE = 'ELEM'), or property entry name, or panel
            flag for ERP responses (PTYPE = 'PANEL' - See Remark 34), or
            RANDPS ID. Blank for grid point responses. 'ELEM' or property
            name used only with element type responses (stress, strain,
            force, etc.) to identify the relevant element IDs, or the property
            type and relevant property IDs.

            Must be {ELEM, PBAR, PSHELL, PCOMP, PANEL, etc.)
            PTYPE = RANDPS ID when RTYPE=PSDDISP, PSDVELO, or PSDACCL.
        region : str
            Region identifier for constraint screening
        atta : int / float / str / blank
            Response attribute
        attb : int / float / str / blank
            Response attribute
        atti : list[int / float / str]
            the response values to pull from
            list[int]:
                list of grid ids
                list of property ids
            list[str]
                'ALL'
        comment : str; default=''
            a comment for the card
        validate : bool; default=True
            should the card be validated when it's created

        Examples
        --------
        **Stress/PSHELL**

        >>> dresp_id = 103
        >>> label = 'resp1'
        >>> response_type = 'STRESS'
        >>> property_type = 'PSHELL'
        >>> pid = 3
        >>> atta = 9 # von mises upper surface stress
        >>> region = None
        >>> attb = None
        >>> atti = [pid]
        >>> DRESP1(dresp_id, label, response_type, property_type, region, atta, attb, atti)


        **Stress/PCOMP**

        >>> dresp_id = 104
        >>> label = 'resp2'
        >>> response_type = 'STRESS'
        >>> property_type = 'PCOMP'
        >>> pid = 3
        >>> layer = 4
        >>> atta = 9 # von mises upper surface stress
        >>> region = None
        >>> attb = layer
        >>> atti = [pid]
        >>> DRESP1(dresp_id, label, response_type, property_type, region, atta, attb, atti)

        """
        assert len(label) <= 8, label
        if isinstance(atti, (integer_types, str)):
            atti = [atti]
        self.cards.append((dresp_id, label, response_type, property_type, region,
                           atta, attb, atti, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DRESP1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        fdouble = force_double if self.model.is_lax_parser else double
        dresp_id = integer(card, 1, 'dresp_id')
        label = string(card, 2, 'label')
        #label = loose_string(card, 2, 'label')
        response_type = string(card, 3, 'rtype')

        # elem, pbar, pshell, etc. (ELEM flag or Prop Name)
        property_type = integer_string_or_blank(card, 4, 'ptype', default='')
        region = integer_or_blank(card, 5, 'region', default=-1)

        atta = integer_double_string_or_blank(card, 6, 'atta', default='')
        attb = integer_double_string_or_blank(card, 7, 'attb', default='')

        atti_list = []
        if response_type in FLOAT_RESPONSE_TYPES:
            #name   atta     attb  atti
            #----   -------- ----- -----------
            #DIVERG mode_num blank mach_number
            for i in range(8, len(card)):
                attii = fdouble(card, i, 'atti_%d' % (i + 1))
                atti_list.append(attii)
        else:
            for i in range(8, len(card)):
                #attii = integer_double_string_or_blank(card, i, 'atti_%d' % (i + 1))
                # DIVERG -> Mach Number
                attii = integer_string_or_blank(card, i, 'atti_%d' % (i + 1))
                atti_list.append(attii)
        #if len(atti_list) == 0:
            #['DRESP1', '10', 'WEIGHT', 'WEIGHT']
            #self.model.log.debug(card)
        #assert len(atti) > 0

        self.cards.append((dresp_id, label, response_type, property_type, region,
                           atta, attb, atti_list, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        dresp_id = np.zeros(ncards, dtype='int32')

        #: user-defined name for printing purposes
        label = np.zeros(ncards, dtype='|U8')
        response_type = np.zeros(ncards, dtype='|U8')

        property_type = np.zeros(ncards, dtype='|U8')
        region = np.zeros(ncards, dtype='int32')

        atta_type = np.zeros(ncards, dtype='|U1')
        atta_int = np.zeros(ncards, dtype='int32')
        atta_float = np.full(ncards, np.nan, dtype='float64')
        atta_str = np.zeros(ncards, dtype='|U8')

        attb_type = np.zeros(ncards, dtype='|U1')
        attb_int = np.zeros(ncards, dtype='int32')
        attb_float = np.full(ncards, np.nan, dtype='float64')
        attb_str = np.zeros(ncards, dtype='|U8')

        atti_type = np.full(ncards, 'i', dtype='|U1')
        atti_ints = []
        atti_strs = []
        atti_floats = []
        iatti = np.zeros((ncards, 2), dtype='int32')
        atti0 = 0
        for icard, card in enumerate(self.cards):
            (dresp_idi, labeli, response_typei, property_typei, regioni,
             attai, attbi, attii, comment) = card
            dresp_id[icard] = dresp_idi
            label[icard] = labeli
            response_type[icard] = response_typei
            property_type[icard] = property_typei

            if regioni is None:
                regioni = -1

            region[icard] = regioni
            if isinstance(attai, int):
                atta_type[icard] = 'i'
                atta_int[icard] = attai
            elif isinstance(attai, str):
                atta_type[icard] = 's'
                atta_str[icard] = attai
            else:
                atta_type[icard] = 'f'
                atta_float[icard] = attai

            if isinstance(attbi, int):
                attb_type[icard] = 'i'
                attb_int[icard] = attbi
            elif isinstance(attbi, str):
                attb_type[icard] = 's'
                attb_str[icard] = attbi
            else:
                attb_type[icard] = 'f'
                attb_float[icard] = attbi

            if response_typei in {'WEIGHT', 'VOLUME'} and attii in ([], ['ALL']):
                attii = [-1]

            natti = len(attii)
            atti1 = atti0 + natti
            if natti:
                if isinstance(attii[0], int):
                    assert response_typei not in FLOAT_RESPONSE_TYPES, response_typei
                    atti_ints.extend(attii)
                    atti_type[icard] = 'i'
                elif isinstance(attii[0], str):
                    assert response_typei in {'VOLUME'}, response_typei
                    atti_strs.extend(attii)
                    atti_type[icard] = 's'
                else:
                    assert isinstance(attii[0], float), attii
                    assert response_typei in FLOAT_RESPONSE_TYPES, response_typei
                    atti_floats.extend(attii)
                    atti_type[icard] = 'f'
            else:
                print(card)
                x = 1
            iatti[icard, :] = [atti0, atti1]  # TODO: change this to natti
            atti1 = atti0

        #if isinstance(attbi, int):
            #atti_type[icard] = 'i'
            #atti_int[icard] = attbi
        #elif isinstance(attbi, str):
            #atti_str[icard] = attii
        #else:
            #atti_float[icard] = attii


        idtype = self.model.idtype
        atti_float = np.array(atti_floats, dtype='float64')
        atti_int = np.array(atti_ints, dtype=idtype)
        atti_str = np.array(atti_strs, dtype='U8')
        ##xinit = np.clip(xinit, xlb, xub)

        self._save(dresp_id, label, response_type, property_type, region,
                   atta_type, atta_float, atta_int, atta_str,
                   attb_type, attb_float, attb_int, attb_str,
                   atti_type, atti_float, atti_int,
                   iatti, )
        #assert len(label) <= 8, f'desvar_id={desvar_id} label={label!r} must be less than 8 characters; length={len(label):d}'
        #assert xlb <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit >= xlb, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #x = 1
        self.cards = []

    def _save(self, dresp_id, label, response_type, property_type, region,
              atta_type, atta_float, atta_int, atta_str,
              attb_type, attb_float, attb_int, attb_str,
              atti_type, atti_float, atti_int, iatti,):
        if len(self.dresp_id) != 0:
            assdf

        self.dresp_id = dresp_id

        self.label = label
        self.response_type = response_type

        self.property_type = property_type
        self.region = region

        self.atta_type = atta_type
        self.atta_int = atta_int
        self.atta_float = atta_float
        self.atta_str = atta_str

        self.attb_type = attb_type
        self.attb_int = attb_int
        self.attb_float = attb_float
        self.attb_str = attb_str

        self.iatti = iatti
        self.atti_type = atti_type
        self.atti_float = atti_float
        self.atti_int = atti_int

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        skip_response_types = {
            'WEIGHT', 'COMPLIAN', 'DWEIGHT', 'FLUTTER', 'LAMA', 'EIGN', 'CEIG', 'STABDER',
            'STRESS', 'STRAIN', 'FORCE', 'CSTRESS', 'CSTRAIN', 'CFAILURE',
        }

        for dresp_id, label, response_type, property_type, region, \
            atta_type, atta_int, atta_float, atta_str, \
            attb_type, attb_int, attb_float, attb_str, \
            atti_type, iatti in zip_longest(
                self.dresp_id, self.label, self.response_type, self.property_type, self.region,
                self.atta_type, self.atta_int, self.atta_float, self.atta_str,
                self.attb_type, self.attb_int, self.attb_float, self.attb_str,
                self.atti_type, self.iatti):
            iatti0, iatti1 = iatti
            if atta_type == 'i':
                atta = atta_int
            elif atta_type == 'f':
                atta = atta_float
            else:
                atta = atta_str

            if attb_type == 'i':
                attb = attb_int
            elif attb_type == 'f':
                attb = attb_float
            else:
                attb = attb_str

            if atti_type == 'i':
                atti = self.atti_int[iatti0:iatti1]
            elif atti_type == 'f':
                atti = atti_float
                atti = atti_float
            else:
                asdf
                atti = atti_str

            if response_type == 'WEIGHT' and len(atti) and atti[0] == -1:
                atti_list = ['ALL']
                #print(self.response_type, self.atti)
                #iatti_all = (self.response_type == 'WEIGHT') & (self.atti == -1)
                #attis[iatti_all] = 'ALL'
            else:
                atti_list = atti.tolist()

            if response_type.upper() in skip_response_types:
                continue
            raise RuntimeError(response_type)
        #used_dict['node_id'].append(self.node_id)
        #used_dict['coord_id'].append(self.coord_id)

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        grid_flags = {
            'FRACCL', 'FRDISP', 'FRSPCF', 'FRVELO',
            'PRES',
            'PSDACCL', 'PSDDISP', 'PSDVELO',
            'RMSACCL', 'RMSDISP', 'RMSVELO', 'SPCFORCE',
            'TACCL', 'TDISP', 'TSPCF', 'TVELO'}
        for response_type, (iatti0, iatti1) in zip(self.response_type, self.iatti):
            if response_type in grid_flags:
                nodes = self.atti[iatti0:iatti1]

                for i, nid1 in enumerate(nodes):
                    nid2 = nid_old_to_new.get(nid1, nid1)
                    nodes[i] = nid2

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dresp_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        regions = array_default_int(self.region, default=-1, size=size)
        atta_ints = array_str(self.atta_int, size=size)
        attb_ints = array_str(self.attb_int, size=size)
        atti_ints = array_str(self.atti_int, size=size)
        atta_floats = array_float(self.atta_float, size=size, is_double=is_double)
        attb_floats = array_float(self.attb_float, size=size, is_double=is_double)
        atti_floats = array_float(self.atti_float, size=size, is_double=is_double)

        for dresp_id, label, response_type, property_type, region, \
            atta_type, atta_int, atta_float, atta_str, \
            attb_type, attb_int, attb_float, attb_str, \
            atti_type, iatti in zip_longest(
                self.dresp_id, self.label, self.response_type, self.property_type, regions,
                self.atta_type, atta_ints, atta_floats, self.atta_str,
                self.attb_type, attb_ints, attb_floats, self.attb_str,
                self.atti_type, self.iatti):
            iatti0, iatti1 = iatti
            if atta_type == 'i':
                atta = atta_int
            elif atta_type == 'f':
                atta = atta_float
            else:
                atta = atta_str

            if attb_type == 'i':
                attb = attb_int
            elif attb_type == 'f':
                attb = attb_float
            else:
                attb = attb_str

            if atti_type == 'i':
                atti_list = atti_ints[iatti0:iatti1].tolist()
                if response_type in {'WEIGHT', 'VOLUME'} and len(atti_list) == 1 and self.atti_int[iatti0] == -1:
                    atti_list = ['ALL']
            elif atti_type == 'f':
                atti_list = atti_floats[iatti0:iatti1].tolist()
            else:
                raise RuntimeError(atti_type)

            list_fields = ['DRESP1', dresp_id, label, response_type, property_type,
                           region, atta, attb] + atti_list
            bdf_file.write(print_card(list_fields))
        return


class DRESP2(VectorizedBaseCard):
    """
    Design Sensitivity Equation Response Quantities
    Defines equation responses that are used in the design, either as
    constraints or as an objective.

    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |   1    |    2    |    3   |     4     |    5   |    6   |    7   |    8   |    9   |
    +========+=========+========+===========+========+========+========+========+========+
    | DRESP2 |   ID    |  LABEL | EQID/FUNC | REGION | METHOD |   C1   |   C2   |   C3   |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DESVAR  | DVID1  |   DVID2   |  DVID3 |  DVID4 |  DVID5 |  DVID6 |  DVID7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DVID8  |   etc.    |        |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DTABLE  | LABL1  |  LABL2    |  LABL3 |  LABL4 |  LABL5 |  LABL6 |  LABL7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | LABL8  |  etc.     |        |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DRESP1  |  NR1   |    NR2    |   NR3  |   NR4  |   NR5  |   NR6  |  NR7   |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         |  NR8   |   etc.    |        |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DNODE   |   G1   |    C1     |   G2   |   C2   |   G3   |   C3   |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         |   G4   |    C4     |  etc.  |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVPREL1 | DPIP1  |   DPIP2   | DPIP3  | DPIP4  | DPIP5  | DPIP6  | DPIP7  |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DPIP8  |   DPIP9   |  etc.  |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVCREL1 | DCIC1  |   DCIC2   | DCIC3  | DCIC4  | DCIC5  | DCIC6  | DCIC7  |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DCIC8  |   DCIC9   |  etc.  |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVMREL1 | DMIM1  |   DMIM2   | DMIM3  | DMIM4  | DMIM5  | DMIM6  | DMIM7  |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DMIM8  |   DMIM9   |  etc.  |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVPREL2 | DPI2P1 |   DPI2P2  | DPI2P3 | DPI2P4 | DPI2P5 | DPI2P6 | DPI2P7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DPI2P8 |   DPI2P9  |  etc.  |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVCREL2 | DCI2C1 |   DCI2C2  | DCI2C3 | DCI2C4 | DCI2C5 | DCI2C6 | DCI2C7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DCI2C8 |   DCI2C9  |   etc. |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVMREL2 | DMI2M1 |   DMI2M2  | DMI2M3 | DMI2M4 | DMI2M5 | DMI2M6 | DMI2M7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DMI2M8 |   DMI2M9  |   etc. |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DRESP2  | NRR1   |   NRR2    |  NRR3  |  NRR4  |  NRR5  |  NRR6  |  NRR7  |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | NRR8   |   etc.    |        |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        | DVLREL1 | DLIL1  |   DLIL2   |  DLIL3 |  DLIL4 |  DLIL5 |  DLIL6 |  DLIL7 |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+
    |        |         | DLIL8  |   etc.    |        |        |        |        |        |
    +--------+---------+--------+-----------+--------+--------+--------+--------+--------+

    C1, C2, C3 are MSC specific
    """
    _id_name = 'dresp_id'
    def clear(self) -> None:
        self.n = 0
        self.dresp_id = np.array([], dtype='int32')

        #: user-defined name for printing purposes
        self.label = np.array([], dtype='|U8')
        self.dequation_id = np.array([], dtype='int32')
        self.dequation_str = np.array([], dtype='|U8')
        self.region = np.array([], dtype='int32')
        self.method = np.array([], dtype='|U8')
        self.c1 = np.array([], dtype='float64')
        self.c2 = np.array([], dtype='float64')
        self.c3 = np.array([], dtype='float64')

        self.nparams = np.array([], dtype='int32')
        self.param_type = np.array([], dtype='|U8')
        self.nparam_values = np.array([], dtype='int32')
        self.param_values = np.array([], dtype='int32')

    def add(self, dresp_id: int, label: str, dequation: int, region: int,
            params: dict[tuple[int, str], list[int]],
            method: str='MIN',
            c1: float=1., c2: float=0.005, c3: float=10.,
            validate: bool=True, comment: str='') -> int:
        """
        Creates a DRESP2 card.

        A DRESP2 is used to define a "complex" output result that may be
        optimized on.  A complex result is a result that uses:
          - simple (DRESP1) results
          - complex (DRESP2) results
          - default values (DTABLE)
          - DVCRELx values
          - DVMRELx values
          - DVPRELx values
          - DESVAR values
        Then, an equation (DEQATN) is used to formulate an output response.

        Parameters
        ----------
        dresp_id : int
            response id
        label : str
            Name of the response
        dequation : int
            DEQATN id
        region : str
            Region identifier for constraint screening
        params : dict[(index, card_type)] = values
            the storage table for the response function
            index : int
                a counter
            card_type : str
                the type of card to pull from
                DESVAR, DVPREL1, DRESP2, etc.
            values : list[int]
                the values for this response
        method : str; default=MIN
            flag used for FUNC=BETA/MATCH
            FUNC = BETA
                valid options are {MIN, MAX}
            FUNC = MATCH
                valid options are {LS, BETA}
        c1 / c2 / c3 : float; default=1. / 0.005 / 10.0
            constants for FUNC=BETA or FUNC=MATCH
        comment : str; default=''
            a comment for the card
        validate : bool; default=False
            should the card be validated when it's created

        params = {
           (0, 'DRESP1') = [10, 20],
           (1, 'DESVAR') = [30],
           (2, 'DRESP1') = [40],
        }

        """
        assert len(label) <= 8, label
        if region is None:
            region = -1
        if method is None:
            method = 'MIN'
        self.cards.append((
            dresp_id, label, dequation, region, params,
            method, c1, c2, c3, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DRESP2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        dresp_id = integer(card, 1, 'dresp_id')
        label = string(card, 2, 'label')
        dequation = integer_or_string(card, 3, 'dequation_id')
        region = integer_or_blank(card, 4, 'region', default=-1)
        method = string_or_blank(card, 5, 'method', default='MIN')

        # MSC 2005   Defaults: C1=100., C2=.005)
        # MSC 2016.1 Defaults: C1=1., C2=.005, C3=10.)
        c1 = double_or_blank(card, 6, 'c1', default=1.)
        c2 = double_or_blank(card, 7, 'c2', default=0.005)
        c3 = double_or_blank(card, 8, 'c3', default=10.)

        fields = [interpret_value(field) for field in card[9:]]

        # DRESP2, dresp_id,
        #         DRESP1, 10, 20
        #         DESVAR, 30
        #         DRESP1, 40
        # params = {
        #    (0, 'DRESP1') = [10, 20],
        #    (1, 'DESVAR') = [30],
        #    (2, 'DRESP1') = [40],
        # }
        params = parse_table_fields('DRESP2', card, fields)

        #print("--DRESP2 Params--")
        #for key, value_list in sorted(params.items()):
            #print("  key=%s value_list=%s" %(key, value_list))
        #return DRESP2(dresp_id, label, dequation, region, params,
                      #method, c1, c2, c3, comment=comment)


        # ------------------------
        self.cards.append((
            dresp_id, label, dequation, region, params,
            method, c1, c2, c3, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        dresp_id = np.zeros(ncards, dtype='int32')

        #: user-defined name for printing purposes
        label = np.zeros(ncards, dtype='|U8')

        dequation_id = np.full(ncards, -1, dtype='int32')
        dequation_str = np.zeros(ncards, dtype='|U8')
        region = np.zeros(ncards, dtype='int32')

        method = np.zeros(ncards, dtype='|U8')
        c1 = np.zeros(ncards, dtype='float64')
        c2 = np.zeros(ncards, dtype='float64')
        c3 = np.zeros(ncards, dtype='float64')

        nparams = np.zeros(ncards, dtype='int32')
        param_type_list = []
        param_values_list = []
        nparam_values_list = []
        for icard, card in enumerate(self.cards):
            (dresp_idi, labeli, dequationi, regioni, paramsi,
             methodi, c1i, c2i, c3i, comment) = card

            nparams[icard] = len(paramsi)
            for key, value in paramsi.items():
                param_type_list.append(key[1])
                if isinstance(value[0], integer_types):
                    assert key[1] != 'DNODE', key
                    pass
                elif isinstance(value[0], list):
                    assert key[1] == 'DNODE', key
                    assert isinstance(value[0][0], integer_types), value
                    value = np.array(value, dtype='int32').flatten().tolist()
                nparam_values_list.append(len(value))
                param_values_list.extend(value)
            #params = {
                #(0, 'DRESP1'): [42],
                #(1, 'DESVAR'): [12],
                #(3, 'DNODE'): [[100, 101],
                               #[1, 2]],
            #}

            if isinstance(dequationi, integer_types):
                dequation_id[icard] = dequationi
            else:
                dequation_str[icard] = dequationi

            dresp_id[icard] = dresp_idi
            label[icard] = labeli
            region[icard] = regioni
            method[icard] = methodi
            c1[icard] = c1i
            c2[icard] = c2i
            c3[icard] = c3i

        param_type = np.array(param_type_list)
        nparam_values = np.array(nparam_values_list)
        param_values = np.array(param_values_list)

        ##xinit = np.clip(xinit, xlb, xub)
        #assert len(label) <= 8, f'desvar_id={desvar_id} label={label!r} must be less than 8 characters; length={len(label):d}'
        #assert xlb <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit >= xlb, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'

        self._save(dresp_id, label, dequation_id, dequation_str,
                   region, method, c1, c2, c3,
                   nparams, param_type, nparam_values, param_values)
        self.cards = []

    def _save(self, dresp_id, label, dequation_id, dequation_str,
              region, method, c1, c2, c3,
              nparams, param_type, nparam_values, param_values):
        self.dresp_id = dresp_id

        #: user-defined name for printing purposes
        self.label = label
        #response_type = np.zeros(ncards, dtype='|U8')

        #self.property_type = np.zeros(ncards, dtype='|U8')
        #self.region = np.zeros(ncards, dtype='int32')

        self.dequation_id = dequation_id
        self.dequation_str = dequation_str
        self.region = region

        self.method = method
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.nparams = nparams
        self.param_type = param_type
        self.nparam_values = nparam_values
        self.param_values = param_values

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        iparam_value = 0

        deqatn_id = self.dequation_id[self.dequation_id > 0]
        used_dict['deqatn_id'].append(deqatn_id)

        for (dresp_id, label, deqatn_id, deqatn_str,
             region, method, iparam) in zip_longest(
            self.dresp_id, self.label,
            self.dequation_id, self.dequation_str,
            self.region, self.method,
            self.iparam):

            iparam0, iparam1 = iparam
            param_types = self.param_type[iparam0:iparam1]
            nparam_values = self.nparam_values[iparam0:iparam1]

            deqatn = deqatn_str if deqatn_id == -1 else deqatn_id
            list_fields = ['DRESP2', dresp_id, label, deqatn,
                           region, method,]

            for param_type, nvalues in zip(param_types, nparam_values):
                values_list2 = self.param_values[iparam_value:iparam_value + nvalues]
                #fields2 = [param_type] + values_list2.tolist()
                (i, j) = DRESP2_PACK_LENGTH[param_type]
                if param_type == 'DRESP1':
                    id_type = 'dresp_id'
                elif param_type == 'DESVAR':
                    id_type = 'desvar_id'
                elif param_type == 'DNODE':
                    id_type = 'desvar_id'
                    # should be (GRID, COMPONENT); only want GRID
                    values_list2 = values_list2[::2]
                else:  # pragma: no cover
                    raise NotImplementedError(param_type)

                used_dict[id_type].append(values_list2)
                #list_fields += build_table_lines(fields2, nstart=i, nend=j)
                iparam_value += nvalues

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        #self.model.log.warning('skipping DRESP2 nodal equivalence')
        iparam_value = 0
        for iparam in self.iparam:
            iparam0, iparam1 = iparam
            param_types = self.param_type[iparam0:iparam1]
            nparam_values = self.nparam_values[iparam0:iparam1]

            for param_type, nvalues in zip(param_types, nparam_values):
                if param_type in {'DRESP1', 'DESVAR', 'DVPREL1', 'DVPREL2', 'DVMREL1', 'DVMREL2'}:
                    continue
                values_list2 = self.param_values[iparam_value:iparam_value + nvalues]
                #print(values_list2)
                for i, nid1 in enumerate(values_list2):
                    if i % 2 == 0:
                        nid2 = nid_old_to_new.get(nid1, nid1)
                        #nodes[i] = nid2
                        values_list2[i] = nid2
                #print(values_list2)
                iparam_value += nvalues

    @property
    def iparam(self) -> np.ndarray:
        return make_idim(self.n, self.nparams)

    #def iparam_value(self, iparam: int) -> np.ndarray:
        #return make_idim(self.n, self.nparam_values[iparam])

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dresp_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        iparam_value = 0

        dresp_ids = array_str(self.dresp_id, size=size)
        regions = array_str(self.region, size=size)
        #pas = array_default_int(self.pa, default=0, size=size)
        methods = array_default_str(self.method, default='MIN', size=size)

        c1s = array_float(self.c1, size=size, is_double=False)
        c2s = array_float(self.c2, size=size, is_double=False)
        c3s = array_float(self.c3, size=size, is_double=False)
        for (dresp_id, label, deqatn_id, deqatn_str,
             region, method, c1, c2, c3, iparam) in zip_longest(
            dresp_ids, self.label,
            self.dequation_id, self.dequation_str,
            regions, methods,
            c1s, c2s, c3s, self.iparam):

            iparam0, iparam1 = iparam
            param_types = self.param_type[iparam0:iparam1]
            nparam_values = self.nparam_values[iparam0:iparam1]

            #method = set_blank_if_default(method, 'MIN')
            #c1 = None
            #c2 = None
            #c3 = None
            #if self.method == 'BETA':
                #if 'BETA' in self.func or 'MATCH' in  self.func:
                    ## MSC 2005   Defaults: C1=100., C2=.005)
                    ## MSC 2016.1 Defaults: C1=1., C2=.005, C3=10.)
                    #c1 = set_blank_if_default(c1, 100.)
                    #c2 = set_blank_if_default(c2, 0.005)
                    #c3 = set_blank_if_default(c3, 10.)

            deqatn = deqatn_str if deqatn_id == -1 else deqatn_id
            list_fields = ['DRESP2', dresp_id, label, deqatn,
                           region, method, c1, c2, c3]

            for param_type, nvalues in zip(param_types, nparam_values):
                values_list2 = self.param_values[iparam_value:iparam_value + nvalues]
                fields2 = [param_type] + values_list2.tolist()
                (i, j) = DRESP2_PACK_LENGTH[param_type]
                list_fields += build_table_lines(fields2, nstart=i, nend=j)
                iparam_value += nvalues

            bdf_file.write(print_card(list_fields))
        return

class DCONSTR(VectorizedBaseCard):
    """
    +---------+------+-----+------------+------------+-------+--------+
    |    1    |   2  |  3  |     4      |      5     |   6   |   7    |
    +=========+======+=====+============+============+=======+========+
    | DCONSTR | DCID | RID | LALLOW/LID | UALLOW/UID | LOWFQ | HIGHFQ |
    +---------+------+-----+------------+------------+-------+--------+
    | DCONSTR |  10  |  4  |    1.25    |            |       |        |
    +---------+------+-----+------------+------------+-------+--------+
    """
    _id_name = 'dconstr_id'
    def clear(self) -> None:
        self.n = 0
        self.dconstr_id = np.array([], dtype='int32')
        self.dresp_id = np.array([], dtype='int32')

        self.lower_allowable = np.array([], dtype='float64')
        self.lower_table = np.array([], dtype='int32')

        self.upper_allowable = np.array([], dtype='float64')
        self.upper_table = np.array([], dtype='int32')

        self.low_frequency = np.array([], dtype='float64')
        self.high_frequency = np.array([], dtype='float64')

    def add(self, dconstr_id: int, dresp_id: int,
            lid: float=-1.e20, uid: float=1.e20,
            lowfq: float=0.0, highfq: float=1.e20,
            comment: str='') -> int:
        """
        Creates a DCONSTR card

        Parameters
        ----------
        oid : int
            unique optimization id
        dresp_id : int
            DRESP1/2 id
        lid / uid=-1.e20 / 1.e20
            lower/upper bound
        lowfq / highfq : float; default=0. / 1.e20
            lower/upper end of the frequency range
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((dconstr_id, dresp_id, lid, uid, lowfq, highfq, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        dconstr_id = integer(card, 1, 'dconstr_id')
        dresp_id = integer(card, 2, 'dresp_id')
        lid = integer_double_or_blank(card, 3, 'lid', default=-1e20)
        uid = integer_double_or_blank(card, 4, 'uid', default=1e20)
        lowfq = double_or_blank(card, 5, 'lowfq', default=0.0)
        highfq = double_or_blank(card, 6, 'highfq', default=1e20)
        assert len(card) <= 7, f'len(DCONSTR card) = {len(card):d}\ncard={card}'
        self.cards.append((dconstr_id, dresp_id, lid, uid, lowfq, highfq, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        dconstr_id = np.zeros(ncards, dtype='int32')
        dresp_id = np.zeros(ncards, dtype='int32')

        lower_allowable = np.full(ncards, np.nan, dtype='float64')
        lower_table = np.full(ncards, 0, dtype='int32')

        upper_allowable = np.full(ncards, np.nan, dtype='float64')
        upper_table = np.full(ncards, 0, dtype='int32')

        low_frequency = np.zeros(ncards, dtype='float64')
        high_frequency = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (dconstr_idi, dresp_idi, lid, uid, lowfq, highfq, comment) = card
            dconstr_id[icard] = dconstr_idi
            dresp_id[icard] = dresp_idi
            if isinstance(lid, int):
                lower_table[icard] = lid
            else:
                lower_allowable[icard] = lid

            if isinstance(uid, int):
                upper_table[icard] = uid
            else:
                upper_allowable[icard] = uid

            low_frequency[icard] = lowfq
            high_frequency[icard] = highfq
        self._save(dconstr_id, dresp_id, lower_table, lower_allowable, upper_table, upper_allowable)

    def _save(self, dconstr_id, dresp_id, lower_table, lower_allowable, upper_table, upper_allowable):
        self.dconstr_id = dconstr_id
        self.dresp_id = dresp_id
        self.lower_table = lower_table
        self.lower_allowable = lower_allowable
        self.upper_table = upper_table
        self.upper_allowable = upper_allowable

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['tabled_id'].append(self.lower_table)
        used_dict['tabled_id'].append(self.upper_table)
        used_dict['dresp_id'].append(self.dresp_id)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dresp_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        dconstr_ids = array_str(self.dconstr_id, size=size)
        dresp_ids = array_str(self.dresp_id, size=size)
        for dconstr_id, dresp_id, lower_allowable, lower_table, upper_allowable, upper_table, \
            low_freq, high_freq in zip_longest(dconstr_ids, dresp_ids,
                                               self.lower_allowable, self.lower_table,
                                               self.upper_allowable, self.upper_table,
                                               self.low_frequency, self.high_frequency):
            lid = lower_table if lower_table > 0 else set_blank_if_default(lower_allowable, -1e20)
            uid = upper_table if upper_table > 0 else set_blank_if_default(upper_allowable, 1e20)

            lowfq = set_blank_if_default(low_freq, 0.0)
            highfq = set_blank_if_default(high_freq, 1e20)
            list_fields = ['DCONSTR', dconstr_id, dresp_id, lid, uid, lowfq, highfq]
            bdf_file.write(print_card(list_fields))
        return


class DVPREL1(VectorizedBaseCard):
    """
    +---------+--------+--------+--------+-----------+-------+--------+-----+
    |   1     |    2   |   3    |    4   |     5     |   6   |   7    |  8  |
    +=========+========+========+========+===========+=======+========+=====+
    | DVPREL1 |   ID   |  TYPE  |  PID   | PNAME/FID | PMIN  |  PMAX  |  C0 |
    +---------+--------+--------+--------+-----------+-------+--------+-----+
    |         | DVID1  | COEF1  | DVID2  |   COEF2   | DVID3 |  etc.  |     |
    +---------+--------+--------+--------+-----------+-------+--------+-----+
    | DVPREL1 | 200000 | PCOMP  | 2000   |     T2    |       |        |     |
    +---------+--------+--------+--------+-----------+-------+--------+-----+
    |         | 200000 |   1.0  |        |           |       |        |     |
    +---------+--------+--------+--------+-----------+-------+--------+-----+
    """
    _id_name = 'dvprel_id'
    def clear(self) -> None:
        self.n = 0
        self.dvprel_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.property_type = np.array([], dtype='|U8')
        self.property_name = np.array([], dtype='float64')
        self.field_num = np.array([], dtype='int32')
        self.p_min = np.array([], dtype='float64')
        self.p_max = np.array([], dtype='float64')
        self.c0 = np.array([], dtype='float64')

        self.ndesvar = np.array([], dtype='int32')
        self.desvar_id = np.array([], dtype='int32')
        self.coefficients = np.array([], dtype='float64')

    def add(self, oid: int, prop_type: str, pid: int, pname_fid: Union[int, str],
            desvar_ids: list[int],
            coeffs: list[float],
            p_min=None, p_max: float=1e20, c0: float=0.0,
            validate: bool=True, comment: str='') -> int:
        """
        Creates a DVPREL1 card

        Parameters
        ----------
        oid : int
            optimization id
        prop_type : str
            property card name (e.g., PSHELL)
        pid : int
            property id
        pname_fid : str/int
            optimization parameter as a pname (property name; T) or field number (fid)
        dvids : list[int]
            DESVAR ids
        coeffs : list[float]
            scale factors for DESVAR ids
        p_min : float; default=None
            minimum property value
        p_max : float; default=1e20
            maximum property value
        c0 : float; default=0.
            offset factor for the variable
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        """
        if isinstance(desvar_ids, integer_types):
            desvar_ids = [desvar_ids]
        if isinstance(coeffs, float_types):
            coeffs = [coeffs]
        assert len(desvar_ids) == len(coeffs), f'desvar_ids={desvar_ids} coeffs={coeffs}'
        card = (oid, prop_type, pid, pname_fid, desvar_ids, coeffs,
                p_min, p_max, c0,
                comment)
        assert oid > 0, oid
        assert pid > 0, pid
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DVPREL1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        oid = integer(card, 1, 'oid')
        prop_type = string(card, 2, 'prop_type')
        pid = integer(card, 3, 'pid')
        pname_fid = integer_or_string(card, 4, 'pName_FID')

        #: Minimum value allowed for this property.
        #: .. todo:: bad default (see DVMREL1)
        p_min = double_or_blank(card, 5, 'p_min', default=np.nan)
        p_max = double_or_blank(card, 6, 'p_max', default=1e20)
        c0 = double_or_blank(card, 7, 'c0', default=0.0)

        desvar_ids = []
        coeffs = []
        end_fields = [interpret_value(field) for field in card[9:]]

        nfields = len(end_fields) # - 1
        #if nfields % 2 == 1:
            #print('end_fields', end_fields)
            #end_fields.append(None)
            #nfields += 1

        i = 0
        for i in range(0, nfields, 2):
            desvar_id = end_fields[i]
            coeff = end_fields[i + 1]
            assert isinstance(desvar_id, integer_types), f'desvar_id={desvar_id} coeff={coeff}; card={card}'
            assert isinstance(coeff, float_types), f'desvar_id={desvar_id} coeff={coeff}; card={card}'
            desvar_ids.append(desvar_id)
            coeffs.append(coeff)

        if len(desvar_ids) != len(coeffs):
        #if nfields % 2 == 1:
            print(card)
            print("desvar_ids = %s" % (desvar_ids))
            print("coeffs = %s" % (coeffs))
            raise RuntimeError('invalid DVPREL1...')

        #return DVPREL1(oid, prop_type, pid, pname_fid, desvars, coeffs,
                       #p_min=p_min, p_max=p_max, c0=c0,
                       #comment=comment)
        card = (oid, prop_type, pid, pname_fid, desvar_ids, coeffs,
                p_min, p_max, c0,
                comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        dvprel_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        property_type = np.zeros(ncards, dtype='|U8')
        property_name = np.zeros(ncards, dtype='|U8')
        field_num = np.zeros(ncards, dtype='int32')
        p_min = np.zeros(ncards, dtype='float64')
        p_max = np.zeros(ncards, dtype='float64')
        c0 = np.zeros(ncards, dtype='float64')
        ndesvar = np.zeros(ncards, dtype='int32')

        all_desvars = []
        all_coeffs = []
        for icard, card in enumerate(self.cards):
            (oid, prop_type, pid, pname_fid, desvars, coeffs,
             p_mini, p_maxi, c0i,
             comment) = card

            dvprel_id[icard] = oid
            property_type[icard] = prop_type
            property_id[icard] = pid
            if isinstance(pname_fid, str):
                property_name[icard] = pname_fid
            else:
                field_num[icard] = pname_fid
            p_min[icard] = p_mini
            p_max[icard] = p_maxi
            c0[icard] = c0i

            ndesvar[icard] = len(desvars)
            all_desvars.extend(desvars)
            all_coeffs.extend(coeffs)

        try:
            desvar_id = np.array(all_desvars, dtype='int32')
        except TypeError:
            print(all_desvars)
            raise
        coefficients = np.array(all_coeffs, dtype='float64')
        self._save(dvprel_id, property_id, property_type,
                   property_name, field_num,
                   p_min, p_max, c0, ndesvar,
                   desvar_id, coefficients)
        self.sort()
        self.cards = []

    def _save(self, dvprel_id, property_id, property_type,
              property_name, field_num,
              p_min, p_max, c0, ndesvar,
              desvar_id, coefficients) -> None:
        if len(self.dvprel_id) != 0:
            dvprel_id = np.hstack([self.dvprel_id, dvprel_id])
            property_id = np.hstack([self.property_id, property_id])
            property_type = np.hstack([self.property_type, property_type])
            field_num = np.hstack([self.field_num, field_num])
            p_min = np.hstack([self.p_min, p_min])
            p_max = np.hstack([self.p_max, p_max])
            c0 = np.hstack([self.c0, c0])
            ndesvar = np.hstack([self.ndesvar, ndesvar])
            desvar_id = np.hstack([self.desvar_id, desvar_id])
            coefficients = np.hstack([self.coefficients, coefficients])
        self.dvprel_id = dvprel_id
        self.property_id = property_id
        self.property_type = property_type
        self.property_name = property_name
        self.field_num = field_num
        self.p_min = p_min
        self.p_max = p_max
        self.c0 = c0

        self.ndesvar = ndesvar
        self.desvar_id = desvar_id
        self.coefficients = coefficients
        assert self.property_id.min() > 0, self.property_id

    def __apply_slice__(self, opt: DVPREL1, i: np.ndarray) -> None:
        opt.dvprel_id = self.dvprel_id[i]
        opt.property_id = self.property_id[i]
        opt.property_type = self.property_type[i]
        opt.property_name = self.property_name[i]
        opt.field_num = self.field_num[i]
        opt.p_min = self.p_min[i]
        opt.p_max = self.p_max[i]
        opt.c0 = self.c0[i]

        idesvar = self.idesvar
        opt.desvar_id = hslice_by_idim(i, idesvar, self.desvar_id)
        opt.coefficients = hslice_by_idim(i, idesvar, self.coefficients)
        opt.ndesvar = self.ndesvar[i]
        opt.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['property_id'].append(self.property_id)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        dvprel_id = used_dict['dvprel_id']
        ncards_removed = remove_unused_primary(self, dvprel_id, self.dvprel_id, 'dvprel_id')
        return ncards_removed

    def geom_check(self, missing: dict[str, np.ndarray]) -> None:
        #ptype_to_pids = {}
        for ptype in np.unique(self.property_type):
            ptype_lower = ptype.lower()
            prop = getattr(self.model, ptype_lower)
            iptype = np.where(ptype == self.property_type)[0]
            #ptype_to_pids[ptype] = self.property_id[iptype]

            # TODO: add desvars
            geom_check(
                self,
                missing,
                property_id=(prop.property_id, self.property_id[iptype]),
            )

    @property
    def idesvar(self) -> np.ndarray:
        return make_idim(self.n, self.ndesvar)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dvprel_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        dvprel_ids = array_str(self.dvprel_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        desvar_ids = array_str(self.desvar_id, size=size)

        c0s = array_default_float(self.c0, default=0., size=size, is_double=False)
        for dvprel_id, pid, prop_type, \
            prop_name, field_num, \
            p_min, p_max, c0, idesvar in zip_longest(dvprel_ids, property_ids, self.property_type,
                                                  self.property_name, self.field_num,
                                                  self.p_max, self.p_min, c0s, self.idesvar):
            idesvar0, idesvar1 = idesvar
            desvars = desvar_ids[idesvar0:idesvar1]
            coeffs = self.coefficients[idesvar0:idesvar1]
            pname_fid = prop_name if prop_name else field_num
            p_max = set_blank_if_default(p_max, 1e20)
            #c0 = set_blank_if_default(c0, 0.)
            list_fields = ['DVPREL1', dvprel_id, prop_type, pid,
                           pname_fid, p_min, p_max, c0, None]
            for (dvid, coeff) in zip_longest(desvars, coeffs):
                list_fields.append(dvid)
                list_fields.append(coeff)
            bdf_file.write(print_card(list_fields))
        return


class DVPREL2(VectorizedBaseCard):
    """
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    |    1     |    2   |   3    |   4   |     5     |   6   |   7   |   8   |   9   |
    +==========+========+========+=======+===========+=======+=======+=======+=======+
    | DVPREL2  | ID     | TYPE   | PID   | PNAME/FID | PMIN  | PMAX  | EQID  |       |
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    |          | DESVAR | DVID1  | DVID2 |   DVID3   | DVID4 | DVID5 | DVID6 | DVID7 |
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    |          |        | DVID8  | etc.  |           |       |       |       |       |
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    |          | DTABLE | LABL1  | LABL2 |   LABL3   | LABL4 | LABL5 | LABL6 | LABL7 |
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    |          |        | LABL8  | etc.  |           |       |       |       |       |
    +----------+--------+--------+-------+-----------+-------+-------+-------+-------+
    """
    _id_name = 'dvprel_id'
    def clear(self) -> None:
        self.n = 0
        self.dvprel_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.property_type = np.array([], dtype='|U8')
        self.property_name = np.array([], dtype='float64')
        self.deqatn_id = np.array([], dtype='int32')
        self.field_num = np.array([], dtype='int32')
        self.p_min = np.array([], dtype='float64')
        self.p_max = np.array([], dtype='float64')


    def add(self, dvprel_id: int, prop_type: str, pid: int,
            pname_fid: Union[int, str], deqation: int,
            desvars: list[int]=None,
            labels: list[str]=None,
            p_min: Optional[float]=None, p_max: float=1.0e20,
            validate: bool=True, comment: str='') -> int:
        """
        Creates a DVPREL2 card

        Parameters
        ----------
        dvprel_id : int
            optimization id
        prop_type : str
            property card name (e.g., PSHELL)
        pid : int
            property id
        pname_fid : str/int
            optimization parameter as a pname (property name; T) or field number (fid)
        deqation : int
            DEQATN id
        dvids : list[int]; default=None
            DESVAR ids
        labels : list[str]; default=None
            DTABLE names
        p_min : float; default=None
            minimum property value
        p_max : float; default=1e20
            maximum property value
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        Notes
        -----
        either dvids or labels is required

        """
        if isinstance(desvars, integer_types):
            desvars = [desvars]

        # DTABLE
        if labels is None:
            labels = []
        elif isinstance(labels, str):
            labels = [labels]
        self.cards.append((dvprel_id, prop_type, pid, deqation, pname_fid,
                           p_min, p_max, desvars, labels, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DVPREL2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        dvprel_id = integer(card, 1, 'oid')
        prop_type = string(card, 2, 'prop_type')
        pid = integer(card, 3, 'pid')
        pname_fid = integer_or_string(card, 4, 'pName_FID')
        p_min = double_or_blank(card, 5, 'p_in', None)
        p_max = double_or_blank(card, 6, 'p_max', 1e20)
        dequation = integer_or_blank(card, 7, 'dequation') #: .. todo:: or blank?

        fields = [interpret_value(field) for field in card[9:]]
        ioffset = 9
        iend = len(fields) + ioffset

        #F:\work\pyNastran\examples\femap_examples\Support\nast\tpl\d200m20.dat
        #params = parse_table_fields('DRESP2', card, fields)
        #print(params)

        try:
            idesvar = fields.index('DESVAR') + ioffset
        except ValueError:
            idesvar = None

        try:
            idtable = fields.index('DTABLE') + ioffset
            #iDesMax  = idtable # the index to start parsing DESVAR
            ides_stop = idtable  # the index to stop  parsing DESVAR
        except ValueError:
            idtable = None
            ides_stop = iend

        desvars = []
        if idesvar:
            n = 1
            for i in range(10, ides_stop):
                dvid_name = 'DVID' + str(n)
                dvid = integer_or_blank(card, i, dvid_name)
                #print("%s = %s" % (dvid_name, dvid))
                if dvid:
                    assert dvid is not None
                    assert dvid != 'DESVAR'
                    desvars.append(dvid)
                    n += 1

        labels = []
        if idtable:
            n = 1
            for i in range(idtable + 1, iend):
                label_name = 'Label' + str(n)
                label = string(card, i, label_name)
                #print("%s = %s" % (label_name, label))
                if label:
                    assert label != 'DTABLE'
                    labels.append(label)

        #dvprel = DVPREL2(oid, prop_type, pid, pname_fid, dequation, dvids, labels,
                         #p_min=p_min, p_max=p_max, comment=comment)
        #if len(dvids) and len(labels) and idtable < idesvar:
            #raise SyntaxError('DESVARs must be defined before DTABLE\n%s' % str(dvprel))

        self.cards.append((dvprel_id, prop_type, pid, dequation, pname_fid,
                           p_min, p_max, desvars, labels, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        dvprel_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        property_type = np.zeros(ncards, dtype='|U8')
        property_name = np.zeros(ncards, dtype='|U8')
        field_num = np.zeros(ncards, dtype='int32')
        deqatn_id = np.zeros(ncards, dtype='int32')
        p_min = np.zeros(ncards, dtype='float64')
        p_max = np.zeros(ncards, dtype='float64')
        ndesvar = np.zeros(ncards, dtype='int32')
        ndtable = np.zeros(ncards, dtype='int32')

        all_desvars = []
        all_labels = []
        for icard, card in enumerate(self.cards):
            (dvprel_idi, prop_type, pid, dequation, pname_fid,
             p_mini, p_maxi, desvars, labels, comment) = card
            dvprel_id[icard] = dvprel_idi
            property_type[icard] = prop_type
            property_id[icard] = pid
            deqatn_id[icard] = dequation
            if isinstance(pname_fid, str):
                property_name[icard] = pname_fid
            else:
                field_num[icard] = pname_fid
            p_min[icard] = p_mini
            p_max[icard] = p_maxi

            ndesvar[icard] = len(desvars)
            ndtable[icard] = len(labels)
            all_desvars.extend(desvars)
            all_labels.extend(labels)

        desvar_ids = np.array(all_desvars, dtype='int32')
        labels = np.array(all_labels, dtype='|U8')
        self._save(dvprel_id, property_id, property_type, property_name, field_num,
                   deqatn_id, p_min, p_max, ndesvar, ndtable, desvar_ids, labels)
        self.cards = []

    def _save(self, dvprel_id, property_id, property_type, property_name, field_num,
              deqatn_id, p_min, p_max, ndesvar, ndtable, desvar_ids, labels):
        if len(self.dvprel_id) != 0:
            asdf
        self.dvprel_id = dvprel_id
        self.property_id = property_id
        self.property_type = property_type
        self.property_name = property_name
        self.field_num = field_num
        self.deqatn_id = deqatn_id
        self.p_min = p_min
        self.p_max = p_max
        #self.c0 = np.zeros(ncards, dtype='float64')
        self.ndesvar = ndesvar
        self.ndtable = ndtable
        self.desvar_ids = desvar_ids
        self.labels = labels

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['property_id'].append(self.property_id)

    def geom_check(self, missing: dict[str, np.ndarray]) -> None:
        #ptype_to_pids = {}
        for ptype in np.unique(self.property_type):
            ptype_lower = ptype.lower()
            prop = getattr(self.model, ptype_lower)

            iptype = np.where(ptype == self.property_type)[0]
            #ptype_to_pids[ptype] = self.property_id[iptype]

            # TODO: add desvars, dtable, deqatn
            geom_check(
                self,
                missing,
                property_id=(prop.property_id, self.property_id[iptype]),
            )

    @property
    def idesvar(self) -> np.ndarray:
        return make_idim(self.n, self.ndesvar)
    @property
    def ilabel(self) -> np.ndarray:
        return make_idim(self.n, self.ndtable)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dvprel_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        dvprel_ids = array_str(self.dvprel_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        desvar_ids = array_str(self.desvar_ids, size=size)

        for dvprel_id, pid, prop_type, \
            prop_name, field_num, deqatn_id, \
            p_min, p_max, idesvar, ilabel in zip_longest(dvprel_ids, property_ids, self.property_type,
                                                         self.property_name, self.field_num, self.deqatn_id,
                                                         self.p_max, self.p_min, self.idesvar, self.ilabel):
            idesvar0, idesvar1 = idesvar
            ilabel0, ilabel1 = ilabel
            desvars = desvar_ids[idesvar0:idesvar1]
            labels = self.labels[ilabel0:ilabel1]
            assert len(desvars) + len(labels) > 0, (desvars, labels)
            pname_fid = prop_name if prop_name else field_num
            p_max = set_blank_if_default(p_max, 1e20)

            list_fields = ['DVPREL2', dvprel_id, prop_type, pid,
                           pname_fid, p_min, p_max, deqatn_id, None]
            if len(desvars):
                fields2 = ['DESVAR'] + desvars.tolist()
                list_fields += build_table_lines(fields2, nstart=1, nend=0)
            if len(labels):
                fields2 = ['DTABLE'] + labels.tolist()
                list_fields += build_table_lines(fields2, nstart=1, nend=0)

            bdf_file.write(print_card(list_fields))
        return


class DVMREL1(VectorizedBaseCard):
    """
    Design Variable to Material Relation
    Defines the relation between a material property and design variables.

    +---------+-------+-------+-------+--------+-------+-------+--------+
    |    1    |   2   |   3   |   4   |    5   |   6   |   7   |    8   |
    +=========+=======+=======+=======+========+=======+=======+========+
    | DVMREL1 |  ID   | TYPE  |  MID  | MPNAME | MPMIN | MPMAX |   C0   |
    +---------+-------+-------+-------+--------+-------+-------+--------+
    |         | DVID1 | COEF1 | DVID2 | COEF2  | DVID3 | COEF3 |  etc.  |
    +---------+-------+-------+-------+--------+-------+-------+--------+
    """
    _id_name = 'dvmrel_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.dvmrel_id = np.array([], dtype='int32')

        self.material_id = np.array([], dtype='int32')
        self.material_type = np.array([], dtype='|U8')
        self.material_name = np.array([], dtype='float64')
        self.mp_min = np.array([], dtype='float64')
        self.mp_max = np.array([], dtype='float64')
        self.c0 = np.array([], dtype='float64')

        self.ndesvar = np.array([], dtype='int32')
        self.desvar_id = np.array([], dtype='int32')
        self.coefficients = np.array([], dtype='float64')

    def add(self, dvmrel_id: int, mat_type: str, mid: int, mp_name: str,
            desvar_ids: list[int], coeffs: list[float],
            mp_min: Optional[float]=None, mp_max: float=1e20,
            c0: float=0., validate: bool=False, comment: str=''):
        """
        Creates a DVMREL1 card

        Parameters
        ----------
        dvmrel_id : int
            optimization id
        mat_type : str
            material card name (e.g., MAT1)
        mid : int
            material id
        mp_name : str
            optimization parameter as a pname (material name; E)
        dvids : list[int]
            DESVAR ids
        coeffs : list[float]
            scale factors for DESVAR ids
        mp_min : float; default=None
            minimum material property value
        mp_max : float; default=1e20
            maximum material property value
        c0 : float; default=0.
            offset factor for the variable
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        """
        if isinstance(desvar_ids, integer_types):
            desvar_ids = [desvar_ids]
        if isinstance(coeffs, float_types):
            coeffs = [coeffs]
        assert len(desvar_ids) == len(coeffs), f'desvar_ids={desvar_ids} coeffs={coeffs}'
        card = (dvmrel_id, mat_type, mid, mp_name, desvar_ids, coeffs,
                mp_min, mp_max, c0, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DVMREL1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        dvmrel_id = integer(card, 1, 'dvmrel_id')
        mat_type = string(card, 2, 'mat_type')
        mid = integer(card, 3, 'mid')
        mp_name = string(card, 4, 'mp_name')
        #if self.mp_name in ['E', 'RHO', 'NU']:  positive values
            #self.mp_min = double_or_blank(card, 5, 'mpMin', 1e-15)
        #else: # negative
            #self.mp_min = double_or_blank(card, 5, 'mpMin', -1e-35)
        mp_min = double_or_blank(card, 5, 'mp_min', default=np.nan)  #: .. todo:: bad default
        mp_max = double_or_blank(card, 6, 'mp_max', default=1e20)
        c0 = double_or_blank(card, 7, 'c0', default=0.0)

        desvar_ids = []
        coeffs = []
        end_fields = [interpret_value(field) for field in card[9:]]
        nfields = len(end_fields) # - 1
        #if nfields % 2 == 1:
            #end_fields.append(None)
            #nfields += 1

        i = 0
        for i in range(0, nfields, 2):
            desvar_id = end_fields[i]
            coeff = end_fields[i + 1]
            assert isinstance(desvar_id, integer_types), card
            assert isinstance(coeff, float_types), card
            desvar_ids.append(desvar_id)
            coeffs.append(coeff)

        if len(desvar_ids) != len(coeffs):
        #if nfields % 2 == 1:
            print(card)
            print("desvar_ids = %s" % (desvar_ids))
            print("coeffs = %s" % (coeffs))
            raise RuntimeError('invalid DVMREL1...')
        #return DVMREL1(dvmrel_id, mat_type, mid, mp_name, desvar_ids, coeffs,
                       #mp_min=mp_min, mp_max=mp_max, c0=c0, comment=comment)
        card = (dvmrel_id, mat_type, mid, mp_name, desvar_ids, coeffs,
                mp_min, mp_max, c0, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        dvmrel_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        material_type = np.zeros(ncards, dtype='|U8')
        material_name = np.zeros(ncards, dtype='|U8')
        mp_min = np.zeros(ncards, dtype='float64')
        mp_max = np.zeros(ncards, dtype='float64')
        c0 = np.zeros(ncards, dtype='float64')
        ndesvar = np.zeros(ncards, dtype='int32')

        all_desvars = []
        all_coeffs = []
        for icard, card in enumerate(self.cards):
            (dvmrel_idi, mat_type, mid, mp_name, desvar_idsi, coeffsi,
             mp_mini, mp_maxi, c0i, comment) = card

            dvmrel_id[icard] = dvmrel_idi
            material_type[icard] = mat_type
            material_id[icard] = mid
            material_name[icard] = mp_name
            mp_min[icard] = mp_mini
            mp_max[icard] = mp_maxi
            c0[icard] = c0i

            ndesvar[icard] = len(desvar_idsi)
            all_desvars.extend(desvar_idsi)
            all_coeffs.extend(coeffsi)
        desvar_id = np.array(all_desvars, dtype='int32')
        coefficients = np.array(all_coeffs, dtype='float64')
        self._save(dvmrel_id, material_id, material_type, material_name,
                   mp_min, mp_max, c0, ndesvar, desvar_id, coefficients)
        self.cards = []

    def _save(self, dvmrel_id, material_id, material_type, material_name,
              mp_min, mp_max, c0, ndesvar, desvar_id, coefficients):
        if len(self.dvmrel_id) != 0:
            assdf
        self.dvmrel_id = dvmrel_id
        self.material_id = material_id
        self.material_type = material_type
        self.material_name = material_name
        self.mp_min = mp_min
        self.mp_max = mp_max
        self.c0 = c0
        self.ndesvar = ndesvar
        self.desvar_id = desvar_id
        self.coefficients = coefficients

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)

    def geom_check(self, missing: dict[str, np.ndarray]) -> None:
        #ptype_to_pids = {}
        for mtype in np.unique(self.material_type):
            mtype_lower = mtype.lower()
            mat = getattr(self.model, mtype_lower)
            imtype = np.where(mtype == self.material_type)[0]
            #ptype_to_pids[ptype] = self.property_id[iptype]

            # TODO: add desvars
            geom_check(
                self,
                missing,
                material_id=(mat.material_id, self.material_id[imtype]),
            )

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.ndesvar)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dvmrel_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        dvmrel_ids = array_str(self.dvmrel_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        desvar_ids = array_str(self.desvar_id, size=size)

        mp_mins = array_float(self.mp_min, size=size)
        mp_maxs = array_default_float(self.mp_max, default=1e20, size=size)
        c0s = array_default_float(self.c0, default=0., size=size)

        for dvmrel_id, mid, mat_type, mp_name, \
            mp_min, mp_max, c0, idim in zip_longest(dvmrel_ids, material_ids, self.material_type,
                                                    self.material_name,
                                                    mp_maxs, mp_mins, c0s, self.idim):
            idim0, idim1 = idim
            desvars = desvar_ids[idim0:idim1]
            coeffs = self.coefficients[idim0:idim1]
            #p_max = set_blank_if_default(p_max, 1e20)
            #c0 = set_blank_if_default(c0, 0.)
            list_fields = ['DVMREL1', dvmrel_id, mat_type, mid,
                           mp_name, mp_min, mp_max, c0, None]
            for (dvid, coeff) in zip_longest(desvars, coeffs):
                list_fields.append(dvid)
                list_fields.append(coeff)
            bdf_file.write(print_card(list_fields))
        return


class DVMREL2(VectorizedBaseCard):
    """
    Design Variable to Material Relation
    Defines the relation between a material property and design variables.

    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    |    1    |    2   |   3    |   4   |     5   |   6   |   7   |   8   |   9   |
    +=========+========+========+=======+=========+=======+=======+=======+=======+
    | DVMREL2 |   ID   | TYPE   |  MID  | MPNAME  | MPMIN | MPMAX | EQID  |       |
    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    |         | DESVAR | DVID1  | DVID2 | DVID3   | DVID4 | DVID5 | DVID6 | DVID7 |
    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    |         | DVID8  |  etc.  |       |         |       |       |       |       |
    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    |         | DTABLE | LABL1  | LABL2 | LABL3   | LABL4 | LABL5 | LABL6 | LABL7 |
    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    |         | LABL8  |  etc.  |       |         |       |       |       |       |
    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    """
    _id_name = 'dvmrel_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.dvmrel_id = np.array([], dtype='int32')

        self.material_id = np.array([], dtype='int32')
        self.material_type = np.array([], dtype='|U8')
        self.material_name = np.array([], dtype='float64')
        self.mp_min = np.array([], dtype='float64')
        self.mp_max = np.array([], dtype='float64')
        self.deqatn_id = np.array([], dtype='float64')

        self.ndesvar = np.array([], dtype='int32')
        self.desvar_id = np.array([], dtype='int32')
        self.coefficients = np.array([], dtype='float64')

    def add(self, dvmrel_id: int, mat_type: str, mid: int, mp_name: str,
            deqatn_id: int, desvar_ids: list[int], labels: list[str],
            mp_min: Optional[float]=None, mp_max: float=1e20,
            validate: bool=True, comment: str='') -> int:
        """
        Creates a DVMREL2 card

        Parameters
        ----------
        dvmrel_id : int
            optimization id
        mat_type : str
            material card name (e.g., MAT1)
        mid : int
            material id
        mp_name : str
            optimization parameter as a pname (material name; E)
        deqatn_id : int
            DEQATN id
        desvar_ids : list[int]; default=None
            DESVAR ids
        labels : list[str]; default=None
            DTABLE names
        mp_min : float; default=None
            minimum material property value
        mp_max : float; default=1e20
            maximum material property value
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        .. note:: either dvids or labels is required

        """
        #if isinstance(desvar_ids, integer_types):
            #desvar_ids = [desvar_ids]
        #if isinstance(coeffs, float_types):
            #coeffs = [coeffs]
        #assert len(desvar_ids) == len(coeffs), f'desvar_ids={desvar_ids} coeffs={coeffs}'
        #card = (dvmrel_id, mat_type, mid, mp_name, dequation, desvar_ids, coeffs,
                #mp_min, mp_max, deqatn_id, comment)
        card = (dvmrel_id, mat_type, mid, mp_name, deqatn_id, desvar_ids, labels,
                mp_min, mp_max, deqatn_id, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DVMREL1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        """
        Adds a DVMREL2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        dvmrel_id = integer(card, 1, 'dvmrel_id')
        mat_type = string(card, 2, 'mat_type')
        mid = integer(card, 3, 'mid')
        mp_name = string(card, 4, 'mp_name')
        #if self.mp_name in ['E', 'RHO', 'NU']:  positive values
            #self.mp_min = double_or_blank(card, 5, 'mpMin', 1e-15)
        #else: # negative
            #self.mp_min = double_or_blank(card, 5, 'mpMin', -1e-35)
        mp_min = double_or_blank(card, 5, 'mp_min', default=np.nan)  #: .. todo:: bad default
        mp_max = double_or_blank(card, 6, 'mp_max', default=1e20)
        deqatn_id = integer_or_blank(card, 7, 'deqatn_id') #: .. todo:: or blank?

        # --------------------------------------------------------------
        fields = [interpret_value(field) for field in card[9:]]
        ioffset = 9
        iend = len(fields) + ioffset

        try:
            idesvar = fields.index('DESVAR') + ioffset
        except ValueError:
            idesvar = None

        try:
            idtable = fields.index('DTABLE') + ioffset
            #iDesMax  = idtable # the index to start parsing DESVAR
            ides_stop = idtable  # the index to stop  parsing DESVAR
        except ValueError:
            idtable = None
            ides_stop = iend

        desvar_ids = []
        if idesvar:
            n = 1
            for i in range(10, ides_stop):
                desvar_id_name = 'DVID' + str(n)
                desvar_id = integer_or_blank(card, i, desvar_id_name)
                #print("%s = %s" % (dvid_name, desvar_id))
                if desvar_id:
                    assert desvar_id is not None
                    assert desvar_id != 'DESVAR'
                    desvar_ids.append(desvar_id)
                    n += 1

        labels = []
        if idtable:
            n = 1
            for i in range(idtable + 1, iend):
                label_name = 'Label' + str(n)
                label = string(card, i, label_name)
                #print("%s = %s" % (label_name, label))
                if label:
                    assert label != 'DTABLE'
                    labels.append(label)
        #return DVMREL2(oid, mat_type, mid, mp_name, dequation, dvids, labels,
                       #mp_min=mp_min, mp_max=mp_max, comment=comment)

        card = (dvmrel_id, mat_type, mid, mp_name, deqatn_id, desvar_ids, labels,
                mp_min, mp_max, deqatn_id, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        dvmrel_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        material_type = np.zeros(ncards, dtype='|U8')
        material_name = np.zeros(ncards, dtype='|U8')
        mp_min = np.zeros(ncards, dtype='float64')
        mp_max = np.zeros(ncards, dtype='float64')
        deqatn_id = np.zeros(ncards, dtype='int32')
        ndesvar = np.zeros(ncards, dtype='int32')
        nlabel = np.zeros(ncards, dtype='int32')

        all_desvars = []
        all_labels = []
        for icard, card in enumerate(self.cards):
            (dvmrel_idi, mat_type, mid, mp_name, deqatn_idi, desvar_idsi, labelsi,
             mp_mini, mp_maxi, c0i, comment) = card

            dvmrel_id[icard] = dvmrel_idi
            material_type[icard] = mat_type
            material_id[icard] = mid
            material_name[icard] = mp_name
            mp_min[icard] = mp_mini
            mp_max[icard] = mp_maxi
            deqatn_id[icard] = deqatn_idi
            if isinstance(desvar_idsi, integer_types):
                desvar_idsi = [desvar_idsi]
            if isinstance(labelsi, str):
                desvar_idsi = [labelsi]

            ndesvar[icard] = len(desvar_idsi)
            nlabel[icard] = len(labelsi)
            all_desvars.extend(desvar_idsi)
            all_labels.extend(labelsi)
        desvar_id = np.array(all_desvars, dtype='int32')
        labels = np.array(all_labels, dtype='|U8')
        self._save(dvmrel_id, material_id, material_type, material_name,
                   mp_min, mp_max, deqatn_id, ndesvar, desvar_id, nlabel, labels)
        self.cards = []

    def _save(self, dvmrel_id, material_id, material_type, material_name,
              mp_min, mp_max, deqatn_id, ndesvar, desvar_id, nlabel, labels):
        if len(self.dvmrel_id) != 0:
            assdf
        self.dvmrel_id = dvmrel_id
        self.material_id = material_id
        self.material_type = material_type
        self.material_name = material_name
        self.mp_min = mp_min
        self.mp_max = mp_max
        self.deqatn_id = deqatn_id
        self.ndesvar = ndesvar
        self.desvar_id = desvar_id
        self.nlabel = nlabel
        self.labels = labels

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)

    def geom_check(self, missing: dict[str, np.ndarray]) -> None:
        #ptype_to_pids = {}
        for mtype in np.unique(self.material_type):
            mtype_lower = mtype.lower()
            mat = getattr(self.model, mtype_lower)
            imtype = np.where(mtype == self.material_type)[0]
            #ptype_to_pids[ptype] = self.property_id[iptype]

            # TODO: add desvars, dtable, deqatn
            geom_check(
                self,
                missing,
                material_id=(mat.material_id, self.material_id[imtype]),
            )

    @property
    def idesvar(self) -> np.ndarray:
        return make_idim(self.n, self.ndesvar)
    @property
    def ilabel(self) -> np.ndarray:
        return make_idim(self.n, self.nlabel)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dvmrel_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        dvmrel_ids = array_str(self.dvmrel_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        desvar_ids = array_str(self.desvar_id, size=size)
        deqatn_ids = array_str(self.deqatn_id, size=size)

        mp_mins = array_float(self.mp_min, size=size)
        mp_maxs = array_default_float(self.mp_max, default=1e20, size=size)

        for dvmrel_id, mid, mat_type, mp_name, \
            mp_min, mp_max, deqatn_id, idim, ilabel in zip_longest(
                dvmrel_ids, material_ids,
                self.material_type, self.material_name,
                mp_maxs, mp_mins,
                deqatn_ids, self.idesvar, self.ilabel):
            idim0, idim1 = idim
            ilabel0, ilabel1 = ilabel
            desvars = desvar_ids[idim0:idim1]
            labels = self.labels[ilabel0:ilabel1]
            #p_max = set_blank_if_default(p_max, 1e20)
            #c0 = set_blank_if_default(c0, 0.)
            list_fields = ['DVMREL2', dvmrel_id, mat_type, mid,
                           mp_name, mp_min, mp_max, deqatn_id, None]
            if len(desvars):
                fields2 = ['DESVAR'] + desvars.tolist()
                list_fields += build_table_lines(fields2, nstart=1, nend=0)
            if len(labels):
                fields2 = ['DTABLE'] + labels.tolist()
                list_fields += build_table_lines(fields2, nstart=1, nend=0)

            bdf_file.write(print_card(list_fields))
        return


class DVCREL1(VectorizedBaseCard):
    """
    Design Variable to Connectivity Relation
    Defines the relation between a element and design variables.

    +---------+-------+-------+-------+--------+-------+-------+--------+
    |    1    |   2   |   3   |   4   |    5   |   6   |   7   |    8   |
    +=========+=======+=======+=======+========+=======+=======+========+
    | DVCREL1 |  ID   | TYPE  |  EID  | CPNAME | CPMIN | CPMAX |   C0   |
    +---------+-------+-------+-------+--------+-------+-------+--------+
    |         | DVID1 | COEF1 | DVID2 | COEF2  | DVID3 | COEF3 |  etc.  |
    +---------+-------+-------+-------+--------+-------+-------+--------+
    """
    _id_name = 'dvcrel_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.dvcrel_id = np.array([], dtype='int32')

        self.element_id = np.array([], dtype='int32')
        self.element_type = np.array([], dtype='|U8')
        self.cp_name = np.array([], dtype='|U8')
        self.cp_min = np.array([], dtype='float64')
        self.cp_max = np.array([], dtype='float64')
        self.c0 = np.array([], dtype='float64')

        self.ndesvar = np.array([], dtype='int32')
        self.desvar_id = np.array([], dtype='int32')
        self.coefficients = np.array([], dtype='float64')

    def add(self, dvcrel_id: int, element_type: str, eid: int, cp_name: str,
            desvar_ids: list[int], coeffs: list[float],
            cp_min: Optional[float]=None, cp_max: float=1e20,
            c0: float=0., validate: bool=False, comment: str=''):
        """
        Creates a DVCREL1 card

        Parameters
        ----------
        dvcrel_id : int
            optimization id
        element_type : str
            material card name (e.g., CONM2)
        eid : int
            material id
        cp_name : str
            optimization parameter as a pname (material name; X1)
        desvar_ids : list[int]
            DESVAR ids
        coeffs : list[float]
            scale factors for DESVAR ids
        cp_min : float; default=None
            minimum material property value
        cp_max : float; default=1e20
            maximum material property value
        c0 : float; default=0.
            offset factor for the variable
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        """
        if isinstance(desvar_ids, integer_types):
            desvar_ids = [desvar_ids]
        if isinstance(coeffs, float_types):
            coeffs = [coeffs]
        assert len(desvar_ids) == len(coeffs), f'desvar_ids={desvar_ids} coeffs={coeffs}'
        card = (dvcrel_id, element_type, eid, cp_name, desvar_ids, coeffs,
                cp_min, cp_max, c0, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DVCREL1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        dvcrel_id = integer(card, 1, 'dvcrel_id')
        element_type = string(card, 2, 'Type')
        eid = integer(card, 3, 'eid')
        cp_name = integer_or_string(card, 4, 'cp_name')

        cp_min = double_or_blank(card, 5, 'cp_min', default=None)
        cp_max = double_or_blank(card, 6, 'cp_max', default=1e20)
        c0 = double_or_blank(card, 7, 'c0', default=0.0)

        desvar_ids = []
        coeffs = []
        end_fields = [interpret_value(field) for field in card[9:]]
        nfields = len(end_fields) # - 1
        #if nfields % 2 == 1:
            #end_fields.append(None)
            #nfields += 1

        i = 0
        for i in range(0, nfields, 2):
            desvar_id = end_fields[i]
            coeff = end_fields[i + 1]
            assert isinstance(desvar_id, integer_types), card
            assert isinstance(coeff, float_types), card
            desvar_ids.append(desvar_id)
            coeffs.append(coeff)

        if len(desvar_ids) != len(coeffs): # nfields % 2 == 1:
            print(card)
            print("desvar_ids = %s" % (desvar_ids))
            print("coeffs = %s" % (coeffs))
            raise RuntimeError('invalid DVCREL1...')
        card = (dvcrel_id, element_type, eid, cp_name, desvar_ids, coeffs,
                cp_min, cp_max, c0, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        dvcrel_id = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        element_type = np.zeros(ncards, dtype='|U8')
        cp_name = np.zeros(ncards, dtype='|U8')
        cp_min = np.zeros(ncards, dtype='float64')
        cp_max = np.zeros(ncards, dtype='float64')
        c0 = np.zeros(ncards, dtype='float64')
        ndesvar = np.zeros(ncards, dtype='int32')

        all_desvars = []
        all_coeffs = []
        for icard, card in enumerate(self.cards):
            (dvmrel_idi, element_typei, eid, cp_namei,
             desvar_idsi, coeffsi,
             cp_mini, cp_maxi, c0i, comment) = card

            dvcrel_id[icard] = dvmrel_idi
            element_type[icard] = element_typei
            element_id[icard] = eid
            cp_name[icard] = cp_namei
            cp_min[icard] = cp_mini
            cp_max[icard] = cp_maxi
            c0[icard] = c0i

            ndesvar[icard] = len(desvar_idsi)
            all_desvars.extend(desvar_idsi)
            all_coeffs.extend(coeffsi)
        desvar_id = np.array(all_desvars, dtype='int32')
        coefficients = np.array(all_coeffs, dtype='float64')
        self._save(dvcrel_id, element_id, element_type, cp_name,
                   cp_min, cp_max, c0, ndesvar, desvar_id, coefficients)
        self.cards = []

    def _save(self, dvcrel_id, element_id, element_type, cp_name,
              cp_min, cp_max, c0, ndesvar, desvar_id, coefficients):
        if len(self.dvcrel_id) != 0:
            assdf
        self.dvcrel_id = dvcrel_id
        self.element_id = element_id
        self.element_type = element_type
        self.cp_name = cp_name
        self.cp_min = cp_min
        self.cp_max = cp_max
        self.c0 = c0
        self.ndesvar = ndesvar
        self.desvar_id = desvar_id
        self.coefficients = coefficients

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['element_id'].append(self.element_id)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        dvcrel_id = used_dict['dvcrel_id']
        ncards_removed = remove_unused_primary(self, dvcrel_id, self.dvcrel_id, 'dvcrel_id')
        return ncards_removed

    def geom_check(self, missing: dict[str, np.ndarray]) -> None:
        #ptype_to_pids = {}
        for etype in np.unique(self.element_type):
            etype_lower = etype.lower()
            elem = getattr(self.model, etype_lower)
            ietype = np.where(etype == self.element_type)[0]
            #ptype_to_pids[ptype] = self.property_id[iptype]

            # TODO: add desvars
            geom_check(
                self,
                missing,
                element_id=(elem.element_id, self.element_id[ietype]),
            )

    @property
    def idim(self) -> np.ndarray:
        return make_idim(self.n, self.ndesvar)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dvcrel_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        dvcrel_ids = array_str(self.dvcrel_id, size=size)
        element_ids = array_str(self.element_id, size=size)
        desvar_ids = array_str(self.desvar_id, size=size)
        cp_mins = array_float(self.cp_min, size=size)
        cp_maxs = array_default_float(self.cp_max, default=1e20, size=size)
        c0s = array_default_float(self.c0, default=0., size=size)

        #coeffs = array_float(self.cp_min, size=size)
        for dvcrel_id, eid, element_type, cp_name, \
            cp_min, cp_max, c0, idim in zip_longest(dvcrel_ids, element_ids, self.element_type,
                                                    self.cp_name,
                                                    cp_maxs, cp_mins, c0s, self.idim):
            idim0, idim1 = idim
            desvars = desvar_ids[idim0:idim1]
            coeffs = self.coefficients[idim0:idim1]
            #p_max = set_blank_if_default(p_max, 1e20)
            #c0 = set_blank_if_default(c0, 0.)
            list_fields = ['DVCREL1', dvcrel_id, element_type, eid,
                           cp_name, cp_min, cp_max, c0, None]
            for (dvid, coeff) in zip_longest(desvars, coeffs):
                list_fields.append(dvid)
                list_fields.append(coeff)
            bdf_file.write(print_card(list_fields))
        return


class DVCREL2(VectorizedBaseCard):
    """
    Design Variable to Material Relation
    Defines the relation between a material property and design variables.

    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    |    1    |    2   |   3    |   4   |     5   |   6   |   7   |   8   |   9   |
    +=========+========+========+=======+=========+=======+=======+=======+=======+
    | DVCREL2 |   ID   | TYPE   |  EID  | CPNAME  | CPMIN | CPMAX | EQID  |       |
    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    |         | DESVAR | DVID1  | DVID2 | DVID3   | DVID4 | DVID5 | DVID6 | DVID7 |
    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    |         | DVID8  |  etc.  |       |         |       |       |       |       |
    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    |         | DTABLE | LABL1  | LABL2 | LABL3   | LABL4 | LABL5 | LABL6 | LABL7 |
    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    |         | LABL8  |  etc.  |       |         |       |       |       |       |
    +---------+--------+--------+-------+---------+-------+-------+-------+-------+
    """
    _id_name = 'dvcrel_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.dvcrel_id = np.array([], dtype='int32')

        self.element_id = np.array([], dtype='int32')
        self.element_type = np.array([], dtype='|U8')
        self.cp_name = np.array([], dtype='float64')
        self.cp_min = np.array([], dtype='float64')
        self.cp_max = np.array([], dtype='float64')
        self.deqatn_id = np.array([], dtype='float64')

        self.ndesvar = np.array([], dtype='int32')
        self.desvar_id = np.array([], dtype='int32')
        self.labels = np.array([], dtype='|U8')

    def add(self, dvcrel_id: int, element_type: str, eid: int, cp_name: str,
            deqatn_id: int, desvar_ids: list[int], labels: list[str],
            cp_min: Optional[float]=None, cp_max: float=1e20,
            validate: bool=True, comment: str='') -> int:
        """
        Creates a DVCREL2 card

        Parameters
        ----------
        dvmrel_id : int
            optimization id
        element_type : str
            material card name (e.g., CONM2)
        eid : int
            material id
        mp_name : str
            optimization parameter as a pname (material name; X2)
        deqatn_id : int
            DEQATN id
        desvar_ids : list[int]; default=None
            DESVAR ids
        labels : list[str]; default=None
            DTABLE names
        cp_min : float; default=None
            minimum material property value
        cp_max : float; default=1e20
            maximum material property value
        validate : bool; default=False
            should the variable be validated
        comment : str; default=''
            a comment for the card

        .. note:: either dvids or labels is required

        """
        #if isinstance(desvar_ids, integer_types):
            #desvar_ids = [desvar_ids]
        #if isinstance(coeffs, float_types):
            #coeffs = [coeffs]
        #assert len(desvar_ids) == len(coeffs), f'desvar_ids={desvar_ids} coeffs={coeffs}'
        card = (dvcrel_id, element_type, eid, cp_name, deqatn_id, desvar_ids, labels,
                cp_min, cp_max, deqatn_id, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DVCREL2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        dvcrel_id = integer(card, 1, 'dvcrel_id')
        element_type = string(card, 2, 'element_type')
        eid = integer(card, 3, 'eid')
        cp_name = string(card, 4, 'cp_name')
        #if self.mp_name in ['E', 'RHO', 'NU']:  positive values
            #self.mp_min = double_or_blank(card, 5, 'mpMin', 1e-15)
        #else: # negative
            #self.mp_min = double_or_blank(card, 5, 'mpMin', -1e-35)
        cp_min = double_or_blank(card, 5, 'cp_min', default=np.nan)  #: .. todo:: bad default
        cp_max = double_or_blank(card, 6, 'cp_max', default=1e20)
        deqatn_id = integer_or_blank(card, 7, 'deqatn_id') #: .. todo:: or blank?

        # --------------------------------------------------------------
        fields = [interpret_value(field) for field in card[9:]]
        ioffset = 9
        iend = len(fields) + ioffset

        try:
            idesvar = fields.index('DESVAR') + ioffset
        except ValueError:
            idesvar = None

        try:
            idtable = fields.index('DTABLE') + ioffset
            #iDesMax  = idtable # the index to start parsing DESVAR
            ides_stop = idtable  # the index to stop  parsing DESVAR
        except ValueError:
            idtable = None
            ides_stop = iend

        desvar_ids = []
        if idesvar:
            n = 1
            for i in range(10, ides_stop):
                desvar_id_name = 'DVID' + str(n)
                desvar_id = integer_or_blank(card, i, desvar_id_name)
                #print("%s = %s" % (dvid_name, desvar_id))
                if desvar_id:
                    assert desvar_id is not None
                    assert desvar_id != 'DESVAR'
                    desvar_ids.append(desvar_id)
                    n += 1

        labels = []
        if idtable:
            n = 1
            for i in range(idtable + 1, iend):
                label_name = 'Label' + str(n)
                label = string(card, i, label_name)
                #print("%s = %s" % (label_name, label))
                if label:
                    assert label != 'DTABLE'
                    labels.append(label)
        #return DVMREL2(oid, mat_type, mid, mp_name, dequation, dvids, labels,
                       #mp_min=mp_min, mp_max=mp_max, comment=comment)

        card = (dvcrel_id, element_type, eid, cp_name, deqatn_id, desvar_ids, labels,
                cp_min, cp_max, deqatn_id, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        dvcrel_id = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        element_type = np.zeros(ncards, dtype='|U8')
        cp_name = np.zeros(ncards, dtype='|U8')
        cp_min = np.zeros(ncards, dtype='float64')
        cp_max = np.zeros(ncards, dtype='float64')
        deqatn_id = np.zeros(ncards, dtype='int32')
        ndesvar = np.zeros(ncards, dtype='int32')
        nlabel = np.zeros(ncards, dtype='int32')

        all_desvars = []
        all_labels = []
        for icard, card in enumerate(self.cards):
            (dvcrel_idi, element_typei, eid, cp_namei, deqatn_idi, desvar_idsi, labelsi,
             cp_mini, cp_maxi, deqatn_idi, comment) = card

            dvcrel_id[icard] = dvcrel_idi
            element_type[icard] = element_typei
            element_id[icard] = eid
            cp_name[icard] = cp_namei
            cp_min[icard] = cp_mini
            cp_max[icard] = cp_maxi
            deqatn_id[icard] = deqatn_idi
            if desvar_idsi is None:
                desvar_idsi = []
            elif isinstance(desvar_idsi, integer_types):
                desvar_idsi = [desvar_idsi]

            if labelsi is None:
                labelsi = []
            elif isinstance(labelsi, str):
                desvar_idsi = [labelsi]

            ndesvar[icard] = len(desvar_idsi)
            nlabel[icard] = len(labelsi)
            all_desvars.extend(desvar_idsi)
            all_labels.extend(labelsi)
        desvar_id = np.array(all_desvars, dtype='int32')
        labels = np.array(all_labels, dtype='|U8')
        self._save(dvcrel_id, element_id, element_type, cp_name,
                   cp_min, cp_max, deqatn_id, ndesvar, desvar_id, nlabel, labels)
        self.cards = []

    def _save(self, dvcrel_id, element_id, element_type, cp_name,
              cp_min, cp_max, deqatn_id, ndesvar, desvar_id, nlabel, labels):
        if len(self.dvcrel_id) != 0:
            assdf
        self.dvcrel_id = dvcrel_id
        self.element_id = element_id
        self.element_type = element_type
        self.cp_name = cp_name
        self.cp_min = cp_min
        self.cp_max = cp_max
        self.deqatn_id = deqatn_id
        self.ndesvar = ndesvar
        self.desvar_id = desvar_id
        self.nlabel = nlabel
        self.labels = labels

    def geom_check(self, missing: dict[str, np.ndarray]) -> None:
        #ptype_to_pids = {}
        for etype in np.unique(self.element_type):
            etype_lower = etype.lower()
            elem = getattr(self.model, etype_lower)
            ietype = np.where(etype == self.element_type)[0]
            #ptype_to_pids[ptype] = self.property_id[iptype]

            # TODO: add desvars, dtable, deqatn
            geom_check(
                self,
                missing,
                element_id=(elem.element_id, self.element_id[ietype]),
            )

    @property
    def idesvar(self) -> np.ndarray:
        return make_idim(self.n, self.ndesvar)
    @property
    def ilabel(self) -> np.ndarray:
        return make_idim(self.n, self.nlabel)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dvcrel_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        dvcrel_ids = array_str(self.dvcrel_id, size=size)
        element_ids = array_str(self.element_id, size=size)
        desvar_ids = array_str(self.desvar_id, size=size)

        cp_mins = array_float(self.cp_min, size=size)
        cp_maxs = array_default_float(self.cp_max, default=1e20, size=size)

        for dvcrel_id, eid, element_type, cp_name, \
            cp_min, cp_max, deqatn_id, idim, ilabel in zip_longest(
                dvcrel_ids, element_ids,
                self.element_type, self.cp_name,
                cp_maxs, cp_mins,
                self.deqatn_id, self.idesvar, self.ilabel):
            idim0, idim1 = idim
            ilabel0, ilabel1 = ilabel
            desvars = desvar_ids[idim0:idim1]
            labels = self.labels[ilabel0:ilabel1]
            #cp_max = set_blank_if_default(cp_max, 1e20)
            list_fields = ['DVCREL2', dvcrel_id, element_type, eid,
                           cp_name, cp_min, cp_max, deqatn_id, None]
            if len(desvars):
                fields2 = ['DESVAR'] + desvars.tolist()
                list_fields += build_table_lines(fields2, nstart=1, nend=0)
            if len(labels):
                fields2 = ['DTABLE'] + labels.tolist()
                list_fields += build_table_lines(fields2, nstart=1, nend=0)

            bdf_file.write(print_card(list_fields))
        return


class DSCREEN(VectorizedBaseCard):
    """
    +---------+-------+-------+------+
    |    1    |   2   |   3   |   4  |
    +=========+=======+=======+======+
    | DSCREEN | RTYPE |  TRS  | NSTR |
    +---------+-------+-------+------+
    | DSCREEN |  DISP | -0.3  | NSTR |
    +---------+-------+-------+------+

    """
    _id_name = 'dscreen_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.dscreen_id = np.array([], dtype='int32')
        self.response_type = np.array([], dtype='|U8')
        self.trs = np.array([], dtype='float64')
        self.nstr = np.array([], dtype='int32')

    def add(self, response_type: str, trs: float=-0.5, nstr: int=20,
            comment: str='') -> int:
        """
        Creates a DSCREEN object

        Parameters
        ----------
        response_type : str
            Response type for which the screening criteria apply
        trs : float
            Truncation threshold
        nstr : int
            Maximum number of constraints to be retained per region per
            load case
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((response_type, trs, nstr, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a DSCREEN card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        response_type = string(card, 1, 'rtype')
        trs = double_or_blank(card, 2, 'trs', default=-0.5)
        nstr = integer_or_blank(card, 3, 'nstr', default=20)
        assert len(card) <= 4, f'len(DSCREEN card) = {len(card):d}\ncard={card}'
        #return DSCREEN(response_type, trs=trs, nstr=nstr, comment=comment)
        self.cards.append((response_type, trs, nstr, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        #dscreen_id = np.zeros(ncards, dtype='int32')
        #: user-defined name for printing purposes
        response_type = np.zeros(ncards, dtype='|U8')
        trs = np.zeros(ncards, dtype='float64')
        nstr = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            response_typei, trsi, nstri, comment = card
            response_type[icard] = response_typei
            trs[icard] = trsi
            nstr[icard] = nstri

        self._save(response_type, trs, nstr)
        self.cards = []

    def _save(self, response_type, trs, nstr):
        if len(self.response_type) != 0:
            asdf
        #self.desvar_id = desvar_id
        self.response_type = response_type
        self.trs = trs
        self.nstr = nstr

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    #def geom_check(self, missing: dict[str, np.ndarray]):
        #pass

    #def index(self, desvar_id: np.ndarray) -> np.ndarray:
        #assert len(self.desvar_id) > 0, self.desvar_id
        #desvar_id = np.atleast_1d(np.asarray(desvar_id, dtype=self.desvar_id.dtype))
        #idesvar = np.searchsorted(self.desvar_id, desvar_id)
        #return idesvar

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.response_type) == 0:
            return
        print_card = get_print_card_8_16(size)

        trss = array_float(self.trs, size=size, is_double=is_double)
        nstrs = array_str(self.nstr, size=size)

        for response_typei, trsi, nstri in zip_longest(self.response_type, trss, nstrs):
            list_fields = ['DSCREEN', response_typei, trsi, nstri]
            bdf_file.write(print_card(list_fields))
        return
