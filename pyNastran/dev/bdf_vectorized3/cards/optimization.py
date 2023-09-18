from __future__ import annotations
from itertools import zip_longest
from typing import Union, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types, float_types
#from pyNastran.bdf import MAX_INT
#from pyNastran2.bdf.cards.base_card import BaseCard, _node_ids, expand_thru
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, integer_or_string,
    string, integer_double_string_or_blank,
    integer_string_or_blank, integer_double_or_blank, interpret_value)
from pyNastran.bdf.field_writer_8 import set_blank_if_default # , print_card_8
#from pyNastran.bdf.field_writer_16 import print_float_16, print_card_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.utils import build_table_lines

from pyNastran.dev.bdf_vectorized3.cards.base_card import VectorizedBaseCard, hslice_by_idim, make_idim, get_print_card_8_16
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class DESVAR(VectorizedBaseCard):
    """
    +--------+-----+-------+-------+-----+-----+-------+-------+
    |    1   |  2  |    3  |   4   |  5  |  6  |    7  |   8   |
    +========+=====+=======+=======+=====+=====+=======+=======+
    | DESVAR | OID | LABEL | XINIT | XLB | XUB | DELXV | DDVAL |
    +--------+-----+-------+-------+-----+-----+-------+-------+
    """
    def __init__(self, model: BDF):
        super().__init__(model)
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
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        desvar_id = integer(card, 1, 'desvar_id')
        label = string(card, 2, 'label')
        xinit = double(card, 3, 'xinit')
        xlb = double_or_blank(card, 4, 'xlb', -1e20)
        xub = double_or_blank(card, 5, 'xub', 1e20)
        delx = double_or_blank(card, 6, 'delx', default=np.nan)
        ddval = integer_or_blank(card, 7, 'ddval', default=0)
        assert len(card) <= 8, f'len(DESVAR card) = {len(card):d}\ncard={card}'
        self.cards.append((desvar_id, label, xinit, xlb, xub, delx, ddval, comment))
        self.n += 1
        return self.n

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards

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
        assert len(self.desvar_id) == 0
        self.desvar_id = desvar_id
        self.label = label
        self.xinit = xinit
        self.xlb = xlb
        self.xub = xub
        self.delx = delx
        self.ddval = ddval

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
        xlbs = self.xlb
        xubs = self.xub

        if is_delx and is_ddval:
            for desvar_id, label, xinit, xlb, xub, delx, ddval in zip_longest(desvar_id, labels, self.xinit, xlbs, xubs,
                                                                              self.delx, self.ddval):
                list_fields = ['DESVAR', desvar_id, label, xinit, xlb, xub,
                               delx, ddval]
                bdf_file.write(print_card(list_fields))
        elif is_delx:
            for desvar_id, label, xinit, xlb, xub, delx in zip_longest(desvar_id, labels, self.xinit, xlbs, xubs, self.delx):
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
    _id_name = 'desvar_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self.desvar_id = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U8')
        self.xinit = np.array([], dtype='float64')
        self.xlb = np.array([], dtype='float64')
        self.xub = np.array([], dtype='float64')
        self.delx = np.array([], dtype='float64')
        self.ddval = np.array([], dtype='float64')

    #def add(self, desvar_id: int, label: str, xinit: float,
            #xlb: float=-1e20, xub: float=1e20,
            #delx=None, ddval: Optional[int]=None,
            #comment: str=''):
        #"""
        #Creates a DESVAR card

        #Parameters
        #----------
        #desvar_id : int
            #design variable id
        #label : str
            #name of the design variable
        #xinit : float
            #the starting point value for the variable
        #xlb : float; default=-1.e20
            #the lower bound
        #xub : float; default=1.e20
            #the lower bound
        #delx : float; default=1.e20
            #fractional change allowed for design variables during
            #approximate optimization
            #NX  if blank : take from DOPTPRM; otherwise 1.0
            #MSC if blank : take from DOPTPRM; otherwise 0.5
        #ddval : int; default=None
            #int : DDVAL id
                  #allows you to set discrete values
            #None : continuous
        #comment : str; default=''
            #a comment for the card

        #"""
        #asdf
        ##self.desvar_id = np.hstack([self.desvar_id, desvar_id])
        ##self.label = np.hstack([self.label, label])
        ##self.xinit = np.hstack([self.xinit, xinit])
        ##self.xlb = np.hstack([self.xlb, xlb])
        ##self.xub = np.hstack([self.xub, xub])
        ##self.delx = np.hstack([self.delx, delx])
        ##self.ddval = np.hstack([self.ddval, ddval])

    def add(self, oid: int, dependent_desvar: int,
            independent_desvars: list[int],
            coeffs: list[float],
            c0: float=0.0, cmult: float=1.0, comment: str='') -> int:
        """
        Creates a DLINK card, which creates a variable that is a lienar
        ccombination of other design variables

        Parameters
        ----------
        oid : int
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
        self.cards.append((oid, dependent_desvar, independent_desvars, coeffs,
                           c0, cmult, comment))
        self.n += 1
        return self.n

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
        desvar_id = integer(card, 1, 'oid')
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
        self.cards.append((desvar_id, dependent_desvar, independent_desvars, coeffs,
                           c0, cmult, comment))
        self.n += 1
        return self.n

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards

        self.desvar_id = np.zeros(ncards, dtype='int32')
        self.dependent_desvar = np.zeros(ncards, dtype='int32')
        self.nindependent_desvars = np.zeros(ncards, dtype='int32')

        self.c0 = np.zeros(ncards, dtype='float64')
        self.cmult = np.zeros(ncards, dtype='float64')
        all_coefficents = []
        all_independent_desvars = []
        for icard, card in enumerate(self.cards):
            desvar_id, dependent_desvar, independent_desvars, coeffs, c0, cmult, comment = card
            self.desvar_id[icard] = desvar_id
            self.dependent_desvar[icard] = dependent_desvar

            all_coefficents.extend(coeffs)
            all_independent_desvars.extend(independent_desvars)
            self.nindependent_desvars[icard] = len(independent_desvars)
            #self.coefficents[icard] = coeffs
            self.c0[icard] = c0
            self.cmult[icard] = cmult
        self.coefficents = np.array(all_coefficents, dtype='float64')
        self.independent_desvars = np.array(all_independent_desvars, dtype='int32')

    @property
    def idesvar(self) -> np.ndarray:
        return make_idim(self.n, self.nindependent_desvars)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.desvar_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        desvar_ids = array_str(self.desvar_id, size=size)
        dependent_desvar = array_str(self.dependent_desvar, size=size)
        for desvar_id, dep_desvar, (idesvar0, idesvar1), xinit, c0, cmult in zip_longest(
            desvar_ids, dependent_desvar, self.idesvar, self.xinit, self.c0, self.cmult):
            c0 = set_blank_if_default(c0, 0.)
            cmult = set_blank_if_default(cmult, 1.)
            independent_desvars = self.independent_desvars[idesvar0 : idesvar1]
            coeffs = self.coefficents[idesvar0 : idesvar1]

            list_fields = ['DLINK', desvar_id, dep_desvar, c0, cmult]
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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.desvar_id = np.array([], dtype='int32')
        self.node_id = np.array([], dtype='int32')
        self.coord_id = np.array([], dtype='int32')
        self.coefficient = np.array([], dtype='float64')
        self.dxyz = np.zeros((0, 3), dtype='float64')

    #def add(self, desvar_id: int, label: str, xinit: float,
            #xlb: float=-1e20, xub: float=1e20,
            #delx=None, ddval: Optional[int]=None,
            #comment: str=''):
        #"""
        #Creates a DESVAR card

        #Parameters
        #----------
        #desvar_id : int
            #design variable id
        #label : str
            #name of the design variable
        #xinit : float
            #the starting point value for the variable
        #xlb : float; default=-1.e20
            #the lower bound
        #xub : float; default=1.e20
            #the lower bound
        #delx : float; default=1.e20
            #fractional change allowed for design variables during
            #approximate optimization
            #NX  if blank : take from DOPTPRM; otherwise 1.0
            #MSC if blank : take from DOPTPRM; otherwise 0.5
        #ddval : int; default=None
            #int : DDVAL id
                  #allows you to set discrete values
            #None : continuous
        #comment : str; default=''
            #a comment for the card

        #"""
        #self.desvar_id = np.hstack([self.desvar_id, desvar_id])
        #self.label = np.hstack([self.label, label])
        #self.xinit = np.hstack([self.xinit, xinit])
        #self.xlb = np.hstack([self.xlb, xlb])
        #self.xub = np.hstack([self.xub, xub])
        #self.delx = np.hstack([self.delx, delx])
        #self.ddval = np.hstack([self.ddval, ddval])

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
        return self.n

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards

        self.desvar_id = np.zeros(ncards, dtype='int32')
        self.node_id = np.zeros(ncards, dtype='int32')
        self.coord_id = np.zeros(ncards, dtype='int32')
        self.coefficient = np.zeros(ncards, dtype='float64')
        self.dxyz = np.zeros((ncards, 3), dtype='float64')
        for icard, card in enumerate(self.cards):
            desvar_id, nid, cid, coeff, dxyz, comment = card
            self.desvar_id[icard] = desvar_id
            self.node_id[icard] = nid
            self.coord_id[icard] = cid
            self.coefficient[icard] = coeff
            self.dxyz[icard, :] = dxyz
        self.cards = []

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.desvar_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        desvar_ids = array_str(self.desvar_id, size=size)
        node_ids = array_str(self.node_id, size=size)
        coord_ids = array_default_int(self.coord_id, default=0, size=size)
        for desvar_id, node_id, coord_id, coeff, dxyz in zip_longest(desvar_ids, node_ids, coord_ids,
                                                                     self.coefficient, self.dxyz.tolist()):
            list_fields = [
                'DVGRID', desvar_id, node_id, coord_id, coeff] + dxyz
            bdf_file.write(print_card(list_fields))
        return


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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.dresp_id = np.array([], dtype='int32')

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        self.dresp_id = np.zeros(ncards, dtype='int32')

        #: user-defined name for printing purposes
        self.label = np.zeros(ncards, dtype='|U8')
        self.response_type = np.zeros(ncards, dtype='|U8')

        self.property_type = np.zeros(ncards, dtype='|U8')
        self.region = np.zeros(ncards, dtype='int32')

        self.atta_type = np.zeros(ncards, dtype='|U1')
        self.atta_int = np.zeros(ncards, dtype='int32')
        self.atta_float = np.full(ncards, np.nan, dtype='float64')
        self.atta_str = np.zeros(ncards, dtype='|U8')

        self.attb_type = np.zeros(ncards, dtype='|U1')
        self.attb_int = np.zeros(ncards, dtype='int32')
        self.attb_float = np.full(ncards, np.nan, dtype='float64')
        self.attb_str = np.zeros(ncards, dtype='|U8')

        attis = []
        self.iatti = np.zeros((ncards, 2), dtype='int32')
        atti0 = 0
        for icard, card_comment in enumerate(self.cards):
            card, comment = card_comment

            dresp_id = integer(card, 1, 'dresp_id')
            label = string(card, 2, 'label')
            #label = loose_string(card, 2, 'label')
            response_type = string(card, 3, 'rtype')

            # elem, pbar, pshell, etc. (ELEM flag or Prop Name)
            property_type = integer_string_or_blank(card, 4, 'ptype', default='')
            region = integer_or_blank(card, 5, 'region', default=-1)

            atta = integer_double_string_or_blank(card, 6, 'atta', default='')
            attb = integer_double_string_or_blank(card, 7, 'attb', default='')

            atti = []
            for i in range(8, len(card)):
                #attii = integer_double_string_or_blank(card, i, 'atti_%d' % (i + 1))
                attii = integer_or_blank(card, i, 'atti_%d' % (i + 1))
                atti.append(attii)
            if len(atti) == 0:
                self.model.log.debug(card)
            #assert len(atti) > 0
            self.dresp_id[icard] = dresp_id
            self.label[icard] = label
            self.response_type[icard] = response_type
            self.property_type[icard] = property_type
            self.region[icard] = region
            if isinstance(atta, int):
                self.attb_type[icard] = 'i'
                self.atta_int[icard] = atta
            elif isinstance(atta, str):
                self.attb_type[icard] = 's'
                self.atta_str[icard] = atta
            else:
                self.attb_type[icard] = 'f'
                self.atta_float[icard] = atta

            if isinstance(attb, int):
                self.attb_type[icard] = 'i'
                self.attb_int[icard] = attb
            elif isinstance(attb, str):
                self.attb_type[icard] = 's'
                self.attb_str[icard] = attb
            else:
                self.attb_type[icard] = 'f'
                self.attb_float[icard] = attb

            natti = len(atti)
            atti1 = atti0 + natti
            attis.extend(atti)
            self.iatti[icard, :] = [atti0, atti1]
            atti1 = atti0
        self.atti = np.array(attis, dtype='int32')
        ##xinit = np.clip(xinit, xlb, xub)
        #assert len(label) <= 8, f'desvar_id={desvar_id} label={label!r} must be less than 8 characters; length={len(label):d}'
        #assert xlb <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit >= xlb, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #x = 1

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dresp_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        regions = array_default_int(self.region, default=-1, size=size)
        atta_ints = array_str(self.atta_int, size=size)
        attb_ints = array_str(self.attb_int, size=size)
        attis = array_str(self.atti, size=size)
        for dresp_id, label, response_type, property_type, region, \
            atta_type, atta_int, atta_float, atta_str, \
            attb_type, attb_int, attb_float, attb_str, iatti in zip_longest(
                self.dresp_id, self.label, self.response_type, self.property_type, regions,
                self.atta_type, atta_ints, self.atta_float, self.atta_str,
                self.attb_type, attb_ints, self.attb_float, self.attb_str, self.iatti):
            iatti0, iatti1 = iatti
            atti = attis[iatti0:iatti1]
            if atta_type == 'i':
                atta = atta_int
            elif atta_type == 'f':
                atta = atta_float
            else:
                atta = atta_str

            if atta_type == 'i':
                attb = attb_int
            elif atta_type == 'f':
                attb = attb_float
            else:
                attb = attb_str

            list_fields = ['DRESP1', dresp_id, label, response_type, property_type,
                           region, atta, attb] + atti.tolist()
            bdf_file.write(print_card(list_fields))
        return


class DRESP2(VectorizedBaseCard):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.dresp_id = np.array([], dtype='int32')

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        self.dresp_id = np.zeros(ncards, dtype='int32')

        #: user-defined name for printing purposes
        self.label = np.zeros(ncards, dtype='|U8')
        #self.response_type = np.zeros(ncards, dtype='|U8')

        #self.property_type = np.zeros(ncards, dtype='|U8')
        #self.region = np.zeros(ncards, dtype='int32')

        #self.atta_type = np.zeros(ncards, dtype='|U1')
        #self.atta_int = np.zeros(ncards, dtype='int32')
        #self.atta_float = np.full(ncards, np.nan, dtype='float64')
        #self.atta_str = np.zeros(ncards, dtype='|U8')

        #self.attb_type = np.zeros(ncards, dtype='|U1')
        #self.attb_int = np.zeros(ncards, dtype='int32')
        #self.attb_float = np.full(ncards, np.nan, dtype='float64')
        #self.attb_str = np.zeros(ncards, dtype='|U8')

        #attis = []
        #self.iatti = np.zeros((ncards, 2), dtype='int32')
        atti0 = 0
        for icard, card_comment in enumerate(self.cards):
            card, comment = card_comment

            dresp_id = integer(card, 1, 'dresp_id')
            label = string(card, 2, 'label')
            #label = loose_string(card, 2, 'label')
            #response_type = string(card, 3, 'rtype')

            ## elem, pbar, pshell, etc. (ELEM flag or Prop Name)
            #property_type = integer_string_or_blank(card, 4, 'ptype', default='')
            #region = integer_or_blank(card, 5, 'region', default=-1)

            #atta = integer_double_string_or_blank(card, 6, 'atta', default='')
            #attb = integer_double_string_or_blank(card, 7, 'attb', default='')

            #atti = []
            #for i in range(8, len(card)):
                ##attii = integer_double_string_or_blank(card, i, 'atti_%d' % (i + 1))
                #attii = integer_or_blank(card, i, 'atti_%d' % (i + 1))
                #atti.append(attii)
            #if len(atti) == 0:
                #print(card)
            #assert len(atti) > 0
            self.dresp_id[icard] = dresp_id
            self.label[icard] = label
            #self.response_type[icard] = response_type
            #self.property_type[icard] = property_type
            #self.region[icard] = region
            #if isinstance(atta, int):
                #self.attb_type[icard] = 'i'
                #self.atta_int[icard] = atta
            #elif isinstance(atta, str):
                #self.attb_type[icard] = 's'
                #self.atta_str[icard] = atta
            #else:
                #self.attb_type[icard] = 'f'
                #self.atta_float[icard] = atta

            #if isinstance(attb, int):
                #self.attb_type[icard] = 'i'
                #self.attb_int[icard] = attb
            #elif isinstance(attb, str):
                #self.attb_type[icard] = 's'
                #self.attb_str[icard] = attb
            #else:
                #self.attb_type[icard] = 'f'
                #self.attb_float[icard] = attb

            #natti = len(atti)
            #atti1 = atti0 + natti
            #attis.extend(atti)
            #self.iatti = [atti0, atti1]
            #atti1 = atti0
        ##xinit = np.clip(xinit, xlb, xub)
        #assert len(label) <= 8, f'desvar_id={desvar_id} label={label!r} must be less than 8 characters; length={len(label):d}'
        #assert xlb <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit >= xlb, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #x = 1

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.dresp_id) == 0:
            return

        print_card = get_print_card_8_16(size)
        for dresp_id, label in zip_longest(self.dresp_id, self.label):
            list_fields = ['DRESP2', dresp_id, label]
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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.dconstr_id = np.array([], dtype='int32')
        self.dresp_id = np.array([], dtype='int32')

        self.lower_allowable = np.array([], dtype='float64')
        self.lower_table = np.array([], dtype='int32')

        self.upper_allowable = np.array([], dtype='float64')
        self.upper_table = np.array([], dtype='int32')

        self.low_frequency = np.array([], dtype='float64')
        self.high_frequency = np.array([], dtype='float64')

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
        return self.n

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
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
    def __init__(self, model: BDF):
        super().__init__(model)
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
        self.cards.append(card)
        self.n += 1
        return self.n

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

        desvars = []
        coeffs = []
        end_fields = [interpret_value(field) for field in card[9:]]

        nfields = len(end_fields) - 1
        if nfields % 2 == 1:
            end_fields.append(None)
            nfields += 1

        i = 0
        for i in range(0, nfields, 2):
            desvars.append(end_fields[i])
            coeffs.append(end_fields[i + 1])
        if nfields % 2 == 1:
            print(card)
            print("desvars = %s" % (desvars))
            print("coeffs = %s" % (coeffs))
            raise RuntimeError('invalid DVPREL1...')

        #return DVPREL1(oid, prop_type, pid, pname_fid, desvars, coeffs,
                       #p_min=p_min, p_max=p_max, c0=c0,
                       #comment=comment)
        card = (oid, prop_type, pid, pname_fid, desvars, coeffs,
                p_min, p_max, c0,
                comment)
        self.cards.append(card)
        self.n += 1
        return self.n

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards

        self.dvprel_id = np.zeros(ncards, dtype='int32')
        self.property_id = np.zeros(ncards, dtype='int32')
        self.property_type = np.zeros(ncards, dtype='|U8')
        self.property_name = np.zeros(ncards, dtype='|U8')
        self.field_num = np.zeros(ncards, dtype='int32')
        self.p_min = np.zeros(ncards, dtype='float64')
        self.p_max = np.zeros(ncards, dtype='float64')
        self.c0 = np.zeros(ncards, dtype='float64')
        self.ndesvar = np.zeros(ncards, dtype='int32')

        all_desvars = []
        all_coeffs = []
        for icard, card in enumerate(self.cards):
            (oid, prop_type, pid, pname_fid, desvars, coeffs,
             p_min, p_max, c0,
             comment) = card

            self.dvprel_id[icard] = oid
            self.property_type[icard] = prop_type
            self.property_id[icard] = pid
            if isinstance(pname_fid, str):
                self.property_name[icard] = pname_fid
            else:
                self.field_num[icard] = pname_fid
            self.p_min[icard] = p_min
            self.p_max[icard] = p_max
            self.c0[icard] = c0

            self.ndesvar[icard] = len(desvars)
            all_desvars.extend(desvars)
            all_coeffs.extend(coeffs)
        self.desvar_id = np.array(all_desvars, dtype='int32')
        self.coefficients = np.array(all_coeffs, dtype='float64')

    @property
    def idim(self) -> np.ndarray:
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

        for dvprel_id, pid, prop_type, \
            prop_name, field_num, \
            p_min, p_max, c0, idim in zip_longest(dvprel_ids, property_ids, self.property_type,
                                                  self.property_name, self.field_num,
                                                  self.p_max, self.p_min, self.c0, self.idim):
            idim0, idim1 = idim
            desvars = desvar_ids[idim0:idim1]
            coeffs = self.coefficients[idim0:idim1]
            pname_fid = prop_name if prop_name else field_num
            p_max = set_blank_if_default(p_max, 1e20)
            c0 = set_blank_if_default(c0, 0.)
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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.dvprel_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.property_type = np.array([], dtype='|U8')
        self.property_name = np.array([], dtype='float64')
        self.deqatn_id = np.array([], dtype='int32')
        self.field_num = np.array([], dtype='int32')
        self.p_min = np.array([], dtype='float64')
        self.p_max = np.array([], dtype='float64')


    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards

        self.dvprel_id = np.zeros(ncards, dtype='int32')
        self.property_id = np.zeros(ncards, dtype='int32')
        self.property_type = np.zeros(ncards, dtype='|U8')
        self.property_name = np.zeros(ncards, dtype='|U8')
        self.field_num = np.zeros(ncards, dtype='int32')
        self.deqatn_id = np.zeros(ncards, dtype='int32')
        self.p_min = np.zeros(ncards, dtype='float64')
        self.p_max = np.zeros(ncards, dtype='float64')
        #self.c0 = np.zeros(ncards, dtype='float64')
        self.ndesvar = np.zeros(ncards, dtype='int32')
        self.ndtable = np.zeros(ncards, dtype='int32')

        all_desvars = []
        all_labels = []
        for icard, card_comment in enumerate(self.cards):
            card, comment = card_comment

            oid = integer(card, 1, 'oid')
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
            self.dvprel_id[icard] = oid
            self.property_type[icard] = prop_type
            self.property_id[icard] = pid
            self.deqatn_id[icard] = dequation
            if isinstance(pname_fid, str):
                self.property_name[icard] = pname_fid
            else:
                self.field_num[icard] = pname_fid
            self.p_min[icard] = p_min
            self.p_max[icard] = p_max
            #self.c0[icard] = c0

            self.ndesvar[icard] = len(desvars)
            self.ndtable[icard] = len(labels)
            all_desvars.extend(desvars)
            all_labels.extend(labels)
        self.desvar_id = np.array(all_desvars, dtype='int32')
        self.labels = np.array(all_labels, dtype='|U8')

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
        desvar_ids = array_str(self.desvar_id, size=size)

        for dvprel_id, pid, prop_type, \
            prop_name, field_num, deqatn_id, \
            p_min, p_max, idesvar, ilabel in zip_longest(dvprel_ids, property_ids, self.property_type,
                                                         self.property_name, self.field_num, self.deqatn_id,
                                                         self.p_max, self.p_min, self.idesvar, self.ilabel):
            idesvar0, idesvar1 = idesvar
            ilabel0, ilabel1 = ilabel
            desvars = desvar_ids[idesvar0:idesvar1]
            labels = self.labels[ilabel0:ilabel1]
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

DVCREL1 = DVCREL2 = DVPREL1
DVMREL1 = DVMREL2 = DVPREL1
