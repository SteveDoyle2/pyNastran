from __future__ import annotations
from itertools import zip_longest
from typing import Union, Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types, float_types
#from pyNastran.bdf import MAX_INT
#from pyNastran2.bdf.cards.base_card import BaseCard, _node_ids, expand_thru
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, integer_or_string,
    string, string_or_blank, integer_double_string_or_blank,
    integer_string_or_blank, integer_double_or_blank, interpret_value)
from pyNastran.bdf.field_writer_8 import set_blank_if_default # , print_card_8
#from pyNastran.bdf.field_writer_16 import print_float_16, print_card_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.cards.utils import build_table_lines

from pyNastran.bdf.cards.optimization import build_table_lines, parse_table_fields, DRESP2_PACK_LENGTH
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, make_idim, get_print_card_8_16,
    #hslice_by_idim,
)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int, array_default_str
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
        assert len(self.desvar_id) == 0
        self.desvar_id = desvar_id
        self.label = label
        self.xinit = xinit
        self.xlb = xlb
        self.xub = xub
        self.delx = delx
        self.ddval = ddval

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

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

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
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

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
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
        if isinstance(atti, integer_types):
            atti = [atti]
        self.cards.append((dresp_id, label, response_type, property_type, region,
                           atta, attb, atti, comment))
        self.n += 1
        return self.n

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
        for i in range(8, len(card)):
            #attii = integer_double_string_or_blank(card, i, 'atti_%d' % (i + 1))
            attii = integer_or_blank(card, i, 'atti_%d' % (i + 1))
            atti_list.append(attii)
        if len(atti_list) == 0:
            self.model.log.debug(card)
        #assert len(atti) > 0

        self.cards.append((dresp_id, label, response_type, property_type, region,
                           atta, attb, atti_list, comment))
        self.n += 1
        return self.n

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
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

        attis = []
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

            natti = len(attii)
            atti1 = atti0 + natti
            attis.extend(attii)
            iatti[icard, :] = [atti0, atti1]
            atti1 = atti0

        atti = np.array(attis, dtype='int32')
        ##xinit = np.clip(xinit, xlb, xub)

        self._save(dresp_id, label, response_type, property_type, region,
                   atta_type, atta_float, atta_int, atta_str,
                   attb_type, attb_float, attb_int, attb_str,
                   iatti, atti)
        #assert len(label) <= 8, f'desvar_id={desvar_id} label={label!r} must be less than 8 characters; length={len(label):d}'
        #assert xlb <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit >= xlb, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #x = 1
        self.cards = []

    def _save(self, dresp_id, label, response_type, property_type, region,
              atta_type, atta_float, atta_int, atta_str,
              attb_type, attb_float, attb_int, attb_str,
              iatti, atti):
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
        self.atti = atti

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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.dresp_id = np.array([], dtype='int32')

        #: user-defined name for printing purposes
        self.label = np.array([], dtype='|U8')
        self.dequation = np.array([], dtype='int32')
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
        return self.n

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
        return self.n

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        dresp_id = np.zeros(ncards, dtype='int32')

        #: user-defined name for printing purposes
        label = np.zeros(ncards, dtype='|U8')
        #response_type = np.zeros(ncards, dtype='|U8')

        #self.property_type = np.zeros(ncards, dtype='|U8')
        #self.region = np.zeros(ncards, dtype='int32')

        dequation = np.zeros(ncards, dtype='int32')
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

            assert isinstance(dequationi, integer_types), dequationi
            dresp_id[icard] = dresp_idi
            label[icard] = labeli
            dequation[icard] = dequationi
            region[icard] = regioni
            method[icard] = methodi
            c1[icard] = c1i
            c2[icard] = c2i
            c3[icard] = c3i
            #asf

        param_type = np.array(param_type_list)
        nparam_values = np.array(nparam_values_list)
        param_values = np.array(param_values_list)

        ##xinit = np.clip(xinit, xlb, xub)
        #assert len(label) <= 8, f'desvar_id={desvar_id} label={label!r} must be less than 8 characters; length={len(label):d}'
        #assert xlb <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit >= xlb, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #assert xinit <= xub, f'desvar_id={desvar_id:d} xlb={xlb} xub={xub}'
        #x = 1
        self._save(dresp_id, label, dequation, region, method, c1, c2, c3,
                   nparams, param_type, nparam_values, param_values)
        self.cards = []

    def _save(self, dresp_id, label, dequation, region, method, c1, c2, c3,
              nparams, param_type, nparam_values, param_values):
        self.dresp_id = dresp_id

        #: user-defined name for printing purposes
        self.label = label
        #response_type = np.zeros(ncards, dtype='|U8')

        #self.property_type = np.zeros(ncards, dtype='|U8')
        #self.region = np.zeros(ncards, dtype='int32')

        self.dequation = dequation
        self.region = region

        self.method = method
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.nparams = nparams
        self.param_type = param_type
        self.nparam_values = nparam_values
        self.param_values = param_values

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

        for (dresp_id, label, deqatn, region, method, c1, c2, c3, iparam) in zip_longest(
            dresp_ids, self.label, self.dequation, regions, methods,
            self.c1, self.c2, self.c3, self.iparam):

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

            if np.isnan(c3):
                c3 = ''

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

    def add(self, dconstr_id: int, dresp_id: int,
            lid: float=-1.e20, uid: float=1.e20,
            lowfq: float=0.0, highfq: float=1.e20, comment: str='') -> int:
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
        return self.n

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

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
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
        return self.n

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
        return self.n

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards

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
    def __init__(self, model: BDF):
        super().__init__(model)
        self.dvmrel_id = np.array([], dtype='int32')

        #self.dvmrel_id = np.zeros(ncards, dtype='int32')
        #self.material_id = np.zeros(ncards, dtype='int32')
        #self.material_type = np.zeros(ncards, dtype='|U8')
        #self.material_name = np.zeros(ncards, dtype='|U8')
        #self.mp_min = np.zeros(ncards, dtype='float64')
        #self.mp_max = np.zeros(ncards, dtype='float64')
        #self.c0 = np.zeros(ncards, dtype='float64')
        #self.ndesvar = np.zeros(ncards, dtype='int32')

        self.material_id = np.array([], dtype='int32')
        self.material_type = np.array([], dtype='|U8')
        self.material_name = np.array([], dtype='float64')
        ##self.field_num = np.array([], dtype='int32')
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
        oid : int
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
        return self.n

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
        mp_name = string(card, 4, 'mpName')
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
        nfields = len(end_fields) - 1
        if nfields % 2 == 1:
            end_fields.append(None)
            nfields += 1

        i = 0
        for i in range(0, nfields, 2):
            desvar_ids.append(end_fields[i])
            coeffs.append(end_fields[i + 1])
        if nfields % 2 == 1:
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
        return self.n

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

        for dvmrel_id, mid, mat_type, mp_name, \
            mp_min, mp_max, c0, idim in zip_longest(dvmrel_ids, material_ids, self.material_type,
                                                    self.material_name,
                                                    self.mp_max, self.mp_min, self.c0, self.idim):
            idim0, idim1 = idim
            desvars = desvar_ids[idim0:idim1]
            coeffs = self.coefficients[idim0:idim1]
            #p_max = set_blank_if_default(p_max, 1e20)
            #c0 = set_blank_if_default(c0, 0.)
            list_fields = ['DVCREL1', dvmrel_id, mat_type, mid,
                           mp_name, mp_min, mp_max, c0, None]
            for (dvid, coeff) in zip_longest(desvars, coeffs):
                list_fields.append(dvid)
                list_fields.append(coeff)
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
    def __init__(self, model: BDF):
        super().__init__(model)
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
        return self.n

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
        nfields = len(end_fields) - 1
        if nfields % 2 == 1:
            end_fields.append(None)
            nfields += 1

        i = 0
        for i in range(0, nfields, 2):
            desvar_ids.append(end_fields[i])
            coeffs.append(end_fields[i + 1])
        if nfields % 2 == 1:
            print(card)
            print("desvar_ids = %s" % (desvar_ids))
            print("coeffs = %s" % (coeffs))
            raise RuntimeError('invalid DVMREL1...')
        card = (dvcrel_id, element_type, eid, cp_name, desvar_ids, coeffs,
                cp_min, cp_max, c0, comment)
        self.cards.append(card)
        self.n += 1
        return self.n

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        dvcrel_id = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype='int32')
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

        for dvcrel_id, eid, element_type, cp_name, \
            cp_min, cp_max, c0, idim in zip_longest(dvcrel_ids, element_ids, self.element_type,
                                                    self.cp_name,
                                                    self.cp_max, self.cp_min, self.c0, self.idim):
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


DVCREL2 = DVPREL1
DVMREL2 = DVPREL1
