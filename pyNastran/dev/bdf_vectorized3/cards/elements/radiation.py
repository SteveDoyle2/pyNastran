from __future__ import annotations
from itertools import zip_longest
from typing import TYPE_CHECKING

import numpy as np
#from pyNastran.bdf.bdf_interface.assign_type_force import force_double, force_double_or_blank
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, # string,
    #blank,
    integer_or_blank, double_or_blank, string_or_blank,
    fields,
)
#from pyNastran.bdf.cards.base_card import BaseCard
#from pyNastran.bdf.cards.elements.bars import set_blank_if_default
#from pyNastran.bdf.cards.properties.bars import _bar_areaL # PBARL as pbarl, A_I1_I2_I12

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, parse_check,
    make_idim, hslice_by_idim, )
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_default_int, # array_default_float,
    array_float, get_print_card_size)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF


class RADCAV(VectorizedBaseCard):
    """
    Identifies the characteristics of each radiant enclosure.

    +--------+---------+--------+------+---------+--------+-------+------+--------+
    |    1   |    2    |    3   |  4   |     5   |    6   |   7   |   8  |    9   |
    +========+=========+========+========+=======+========+=======+======+========+
    | RADCAV | ICAVITY | ELEAMB | SHADOW | SCALE | PRTPCH | NFECI | RMAX |        |
    +--------+---------+--------+------+---------+--------+-------+------+--------+
    |        |  SET11  |  SET12 |  SET21 | SET22 |  SET31 | SET32 | etc. |        |
    +--------+---------+--------+------+---------+--------+-------+------+--------+
    | RADCAV |    1    |    1   |        |       |        |       | .99  |        |
    +--------+---------+--------+------+---------+--------+-------+------+--------+
    |        |    3    |    5   |    4   |   5   |    7   |   5   |      |        |
    +--------+---------+--------+------+---------+--------+-------+------+--------+

    """
    _id_name = 'icavity'
    def clear(self) -> None:
        self.icavity = np.array([], dtype='int32')
        self.ele_ambient = np.array([], dtype='int32')
        self.shadow = np.array([], dtype='|U8')
        self.scale = np.array([], dtype='float64')

        self.prtpch = np.array([], dtype='int32')
        self.nefci = np.array([], dtype='|U8')
        self.rmax = np.array([], dtype='float64')
        self.ncomp = np.array([], dtype='int32')
        self.nset = np.array([], dtype='int32')
        self.sets = np.array([], dtype='int32')
        self.n = 0

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a RADCAV card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        icavity = integer(card, 1, 'icavity')
        ele_ambient = integer_or_blank(card, 2, 'ele_amb', default=0)  # default made up
        shadow = string_or_blank(card, 3, 'shadow', default='YES')
        scale = double_or_blank(card, 4, 'scale', default=0.0)
        prtpch = integer_or_blank(card, 5, 'prtpch', default=0)    # default made up
        nefci = string_or_blank(card, 6, 'nefci', default='')      # default made up
        rmax = double_or_blank(card, 7, 'rmax', default=1.0)
        ncomp = integer_or_blank(card, 8, 'ncomp', default=32)

        sets = fields(integer, card, 'set', i=9, j=card.nfields)
        if not self.model.allow_empty_sets:
            assert len(sets) > 0, card
        #return RADCAV(icavity, sets, ele_amb=ele_amb,
                      #shadow=shadow, scale=scale, prtpch=prtpch,
                      #nefci=nefci, rmax=rmax, ncomp=ncomp, comment=comment)
        self.cards.append((icavity, sets, ele_ambient,
                           shadow, scale, prtpch,
                           nefci, rmax, ncomp, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        icavity = np.zeros(ncards, dtype='int32')
        ele_ambient = np.zeros(ncards, dtype='int32')
        shadow = np.zeros(ncards, dtype='|U8')
        scale = np.zeros(ncards, dtype='float64')

        prtpch = np.zeros(ncards, dtype='int32')
        nefci = np.zeros(ncards, dtype='|U8')
        rmax = np.zeros(ncards, dtype='float64')
        ncomp = np.zeros(ncards, dtype='int32')
        nset = np.zeros(ncards, dtype='int32')
        all_sets_list = []

        for icard, card in enumerate(self.cards):
            (icavityi, setsi, ele_ambi,
             shadowi, scalei, prtpchi,
             nefcii, rmaxi, ncompi, comment) = card
            icavity[icard] = icavityi
            ele_ambient[icard] = ele_ambi
            shadow[icard] = shadowi
            scale[icard] = scalei
            prtpch[icard] = prtpchi
            nefci[icard] = nefcii
            rmax[icard] = rmaxi
            ncomp[icard] = ncompi
            nset[icard] = len(setsi)
            all_sets_list.extend(setsi)
        all_sets = np.array(all_sets_list, dtype=idtype)
        self._save(icavity, ele_ambient, shadow, scale, prtpch,
                   nefci, rmax, ncomp, all_sets, nset)
        self.cards = []

    def _save(self, icavity, ele_ambient, shadow, scale, prtpch,
              nefci, rmax, ncomp, sets, nset) -> int:
        assert len(self.icavity) == 0
        self.icavity = icavity
        self.ele_ambient = ele_ambient
        self.shadow = shadow
        self.scale = scale
        self.prtpch = prtpch
        self.nefci = nefci
        self.rmax = rmax
        self.ncomp = ncomp
        self.sets = sets
        self.nset = nset
        self.n = len(icavity)

    def __apply_slice__(self, radcav: RADCAV, i: np.ndarray) -> None:
        radcav.icavity = self.icavity[i]
        radcav.ele_ambient = self.ele_ambient[i]
        radcav.shadow = self.shadow[i]
        radcav.scale = self.scale[i]
        radcav.prtpch = self.prtpch[i]
        radcav.nefci = self.nefci[i]
        radcav.rmax = self.rmax[i]
        radcav.ncomp = self.ncomp[i]

        iset = self.iset
        radcav.sets = hslice_by_idim(i, iset, self.sets)
        radcav.nset = self.nset[i]
        radcav.n = len(i)

    @property
    def iset(self) -> np.ndarray:
        return make_idim(self.n, self.nset)

    @property
    def max_id(self) -> int:
        set_max = 0 if len(self.sets) == 0 else self.sets.max()
        return max(self.icavity.max(), set_max)

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        icavitys = array_str(self.icavity, size=size)
        ele_ambients = array_str(self.ele_ambient, size=size)
        shadows = array_str(self.shadow, size=size)
        prtpchs = array_str(self.prtpch, size=size)
        ncomps = array_str(self.ncomp, size=size)
        for icavity, ele_ambient, shadow, scale, \
            prtpch, nefci, rmax, ncomp, (iset0, iset1) in zip_longest(icavitys, ele_ambients,
                                                                      shadows, self.scale,
                                                                      prtpchs, self.nefci, self.rmax,
                                                                      ncomps, self.iset):
            sets = self.sets[iset0:iset1].tolist()
            list_fields = ['RADCAV', icavity, ele_ambient, shadow, scale,
                           prtpch, nefci, rmax, ncomp] + sets
            bdf_file.write(print_card(list_fields))
        return


class RADLST(VectorizedBaseCard):
    """
    Identifies the characteristics of each radiant enclosure.

    +--------+---------+--------+------+-------+--------+-------+------+--------+
    |    1   |    2    |    3   |  4   |   5   |    6   |   7   |   8  |    9   |
    +========+=========+========+======+=======+========+=======+======+========+
    | RADLST | ICAVITY | MTXTYP | EID1 |  EID2 |  EID3  |  EID4 | EID5 |  EID6  |
    +--------+---------+--------+------+-------+--------+-------+------+--------+
    |        |   EID7  |  etc.  |      |       |        |       |      |        |
    +--------+---------+--------+------+-------+--------+-------+------+--------+
    | RADLST |    3    |    5   |  4   |   5   |    7   |   5   |      |        |
    +--------+---------+--------+------+-------+--------+-------+------+--------+

    """
    _id_name = 'icavity'
    def clear(self) -> None:
        self.icavity = np.array([], dtype='int32')
        self.ele_ambient = np.array([], dtype='int32')
        self.shadow = np.array([], dtype='|U8')
        self.scale = np.array([], dtype='float64')

        self.prtpch = np.array([], dtype='int32')
        self.nefci = np.array([], dtype='|U8')
        self.rmax = np.array([], dtype='float64')
        self.ncomp = np.array([], dtype='int32')
        self.nset = np.array([], dtype='int32')
        self.sets = np.array([], dtype='int32')
        self.n = 0

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a RADLST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        icavity = integer(card, 1, 'icavity')
        matrix_type = integer_or_blank(card, 2, 'matrix_type', default=1)

        eids = fields(integer, card, 'eid', i=3, j=card.nfields)

        if not self.model.allow_empty_sets:
            assert len(eids) > 0, card
        self.cards.append((icavity, matrix_type, eids, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        icavity = np.zeros(ncards, dtype='int32')
        nelement = np.zeros(ncards, dtype='int32')

        matrix_type = np.zeros(ncards, dtype='int32')
        all_eids_list = []

        for icard, card in enumerate(self.cards):
            (icavityi, matrix_typei, eids, comment) = card
            icavity[icard] = icavityi
            matrix_type[icard] = matrix_typei
            nelement[icard] = len(eids)
            all_eids_list.extend(eids)
        element_ids = np.array(all_eids_list, dtype=idtype)
        self._save(icavity, matrix_type, element_ids, nelement)
        self.cards = []

    def _save(self, icavity, matrix_type, element_ids, nelement) -> int:
        assert len(self.icavity) == 0
        self.icavity = icavity
        self.matrix_type = matrix_type
        self.element_ids = element_ids
        self.nelement = nelement
        self.n = len(icavity)

    def __apply_slice__(self, radlist: RADLST, i: np.ndarray) -> None:
        radlist.icavity = self.icavity[i]
        radlist.matrix_type = self.matrix_type[i]
        ielement = self.ielement
        radlist.element_ids = hslice_by_idim(i, ielement, self.element_ids)
        radlist.nelement = self.nelement[i]
        radlist.n = len(i)

    @property
    def ielement(self) -> np.ndarray:
        return make_idim(self.n, self.nelement)

    @property
    def max_id(self) -> int:
        eid_max = 0 if len(self.element_ids) == 0 else self.element_ids.max()
        return max(self.icavity.max(), eid_max)

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        icavitys = array_str(self.icavity, size=size)
        matrix_types = array_default_int(self.ele_ambient, default=1, size=size)
        element_ids = array_str(self.element_ids, size=size).tolist()
        for icavity, matrix_type, (ielement0, ielement1) in zip_longest(icavitys, matrix_types,
                                                                        self.ielement):
            eids = element_ids[ielement0:ielement1]
            list_fields = ['RADLST', icavity, matrix_type] + eids
            bdf_file.write(print_card(list_fields))
        return


class RADMTX(VectorizedBaseCard):
    """
    Provides the Fji=Aj*fji exchange factors for all the faces of a
    radiation enclosure specified in the corresponding RADLST entry.

    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |    1   |    2    |    3   |   4    |    5   |    6   |    7   |    8   |    9   |
    +========+=========+========+========+========+========+========+========+========+
    | RADMTX | ICAVITY | INDEX  |  Fi,j  | Fi+1,j | Fi+2,j | Fi+3,j | Fi+4,j | Fi+5,j |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    |        | Fi+6,j  |  etc.  |        |        |        |        |        |        |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+
    | RADMTX |    2    |    1   |  0.0   |  0.1   |   0.2  |   0.2  |   0.3  |   0.2  |
    +--------+---------+--------+--------+--------+--------+--------+--------+--------+

    """
    _id_name = 'icavity'
    def clear(self) -> None:
        self.icavity = np.array([], dtype='int32')
        self.ele_ambient = np.array([], dtype='int32')
        self.shadow = np.array([], dtype='|U8')
        self.scale = np.array([], dtype='float64')

        self.prtpch = np.array([], dtype='int32')
        self.nefci = np.array([], dtype='|U8')
        self.rmax = np.array([], dtype='float64')
        self.ncomp = np.array([], dtype='int32')
        self.nset = np.array([], dtype='int32')
        self.sets = np.array([], dtype='int32')
        self.n = 0

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a RADMTX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        icavity = integer(card, 1, 'icavity')
        index = integer(card, 2, 'index')
        exchange_factors = fields(double, card, 'eid', i=3, j=card.nfields)
        if not self.model.allow_empty_sets:
            assert len(exchange_factors) > 0, card
        #return RADMTX(icavity, index, exchange_factors, comment=comment)
        self.cards.append((icavity, index, exchange_factors, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        #idtype = self.model.idtype
        icavity = np.zeros(ncards, dtype='int32')
        nexchange_factors = np.zeros(ncards, dtype='int32')
        index = np.zeros(ncards, dtype='int32')

        exchange_factors_list = []
        for icard, card in enumerate(self.cards):
            (icavityi, indexi, exchange_factorsi, comment) = card
            icavity[icard] = icavityi
            index[icard] = indexi
            nexchange_factors[icard] = len(exchange_factors_list)
            exchange_factors_list.extend(exchange_factorsi)
        exchange_factors = np.array(exchange_factors_list, dtype='float64')
        self._save(icavity, index, exchange_factors, nexchange_factors)
        self.cards = []

    def _save(self, icavity, index, exchange_factors, nexchange_factors) -> int:
        assert len(self.icavity) == 0
        self.icavity = icavity
        self.index = index
        self.exchange_factors = exchange_factors
        self.nexchange_factors = nexchange_factors
        self.n = len(icavity)

    def __apply_slice__(self, radmtx: RADLST, i: np.ndarray) -> None:
        radmtx.icavity = self.icavity[i]
        radmtx.index = self.index[i]
        iexchange_factor = self.iexchange_factor
        radmtx.exchange_factors = hslice_by_idim(i, iexchange_factor, self.exchange_factors)
        radmtx.nexchange_factors = self.nexchange_factors[i]
        radmtx.n = len(i)

    @property
    def iexchange_factor(self) -> np.ndarray:
        return make_idim(self.n, self.nexchange_factors)

    @property
    def max_id(self) -> int:
        #eid_max = 0 if len(self.element_ids) == 0 else self.element_ids.max()
        return self.icavity.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        icavitys = array_str(self.icavity, size=size)
        indexs = array_str(self.index, size=size)
        exchange_factors = array_float(self.exchange_factors, size=size, is_double=False).tolist()
        for icavity, index, (iex0, iex1) in zip_longest(icavitys, indexs, self.iexchange_factor):
            exchange_factorsi = exchange_factors[iex0:iex1]
            list_fields = ['RADMTX', icavity, index] + exchange_factorsi
            bdf_file.write(print_card(list_fields))
        return


class VIEW(VectorizedBaseCard):
    """
    Defines radiation cavity and shadowing for radiation
    view factor calculations.

    +------+-------+---------+-------+----+----+--------+
    |   1  |   2   |    3    |   4   | 5  |  6 |    7   |
    +======+=======+=========+=======+====+====+========+
    | VIEW | IVIEW | ICAVITY | SHADE | NB | NG | DISLIN |
    +------+-------+---------+-------+----+----+--------+
    | VIEW |   1   |    1    | BOTH  | 2  | 3  |  0.25  |
    +------+-------+---------+-------+----+----+--------+

    """
    _id_name = 'icavity'
    def clear(self) -> None:
        self.iview = np.array([], dtype='int32')
        self.icavity = np.array([], dtype='int32')
        self.shade = np.array([], dtype='|U8')
        self.nbeta = np.array([], dtype='int32')
        self.ngamma = np.array([], dtype='int32')
        self.dislin = np.array([], dtype='float64')
        self.n = 0

    def add(self, iview: int, icavity: int, shade: str='BOTH',
                 nbeta: int=1, ngamma: int=1, dislin: float=0.0, comment: str='') -> int:
        """Creates a VIEW card"""
        self.cards.append((iview, icavity, shade, nbeta, ngamma,
                           dislin, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a VIEW card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        iview = integer(card, 1, 'iview')
        icavity = integer(card, 2, 'icavity')
        shade = string_or_blank(card, 3, 'shade', default='BOTH')
        nbeta = integer_or_blank(card, 4, 'nbeta', default=1)
        ngamma = integer_or_blank(card, 5, 'ngamma', default=1)
        dislin = double_or_blank(card, 6, 'dislin', default=0.0)
        #return VIEW(iview, icavity, shade=shade, nbeta=nbeta, ngamma=ngamma,
                    #dislin=dislin, comment=comment)

        self.cards.append((iview, icavity, shade, nbeta, ngamma,
                           dislin, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        iview = np.zeros(ncards, dtype='int32')
        icavity = np.zeros(ncards, dtype='int32')
        shade = np.zeros(ncards, dtype='|U8')
        nbeta = np.zeros(ncards, dtype='int32')
        ngamma = np.zeros(ncards, dtype='int32')
        dislin = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (iviewi, icavityi, shadei, nbetai, ngammai, dislini, comment) = card
            iview[icard] = iviewi
            icavity[icard] = icavityi
            shade[icard] = shadei
            nbeta[icard] = nbetai
            ngamma[icard] = ngammai
            dislin[icard] = dislini
        self._save(iview, icavity, shade, nbeta, ngamma, dislin)
        self.cards = []

    def _save(self, iview, icavity, shade, nbeta, ngamma, dislin) -> int:
        assert len(self.icavity) == 0
        self.iview = iview
        self.icavity = icavity
        self.shade = shade
        self.nbeta = nbeta
        self.ngamma = ngamma
        self.dislin = dislin
        self.n = len(iview)

    def __apply_slice__(self, view: VIEW, i: np.ndarray) -> None:
        view.iview = self.iview[i]
        view.icavity = self.icavity[i]
        view.shade = self.shade[i]
        view.nbeta = self.nbeta[i]
        view.ngamma = self.ngamma[i]
        view.dislin = self.dislin[i]
        view.n = len(i)

    @property
    def max_id(self) -> int:
        return max(self.iview.max(), self.icavity.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        iviews = array_str(self.iview, size=size)
        icavitys = array_str(self.icavity, size=size)
        nbetas = array_str(self.nbeta, size=size)
        ngammas = array_str(self.ngamma, size=size)
        for iview, icavity, shade, nbeta, ngamma, dislin in zip_longest(
            iviews, icavitys, self.shade, nbetas, ngammas, self.dislin):
            list_fields = ['VIEW', iview, icavity, shade, nbeta, ngamma, dislin]
            bdf_file.write(print_card(list_fields))
        return


class VIEW3D(VectorizedBaseCard):
    """
    View Factor Definition - Gaussian Integration Method

    Defines parameters to control and/or request the Gaussian Integration
    method of view factor calculation for a specified cavity.

    +--------+---------+------+------+------+------+--------+------+--------+
    |    1   |    2    |   3  |  4   |   5  |   6  |    7   |   8  |    9   |
    +========+=========+======+======+======+======+========+======+========+
    | VIEW3D | ICAVITY | GITB | GIPS | CIER | ETOL |  ZTOL  | WTOL | RADCHK |
    +--------+---------+------+------+------+------+--------+------+--------+
    | VIEW3D |    1    |   2  |   2  |   4  |      | 1.0E-6 |      |        |
    +--------+---------+------+------+------+------+--------+------+--------+

    """
    _id_name = 'icavity'
    def clear(self) -> None:
        self.icavity = np.array([], dtype='int32')
        self.gitb = np.array([], dtype='int32')
        self.gips = np.array([], dtype='int32')
        self.cier = np.array([], dtype='int32')
        self.error_tol = np.array([], dtype='float64')
        self.warp_tol = np.array([], dtype='float64')
        self.rad_check = np.array([], dtype='int32')
        self.n = 0

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a VIEW3D card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        icavity = integer(card, 1, 'icavity')
        gitb = integer_or_blank(card, 2, 'gitb', default=4)
        gips = integer_or_blank(card, 3, 'gips', default=4)
        cier = integer_or_blank(card, 4, 'cier', default=4)
        error_tol = double_or_blank(card, 5, 'error_tol', default=0.1)
        zero_tol = double_or_blank(card, 6, 'zero_tol', default=1e-10)
        warp_tol = double_or_blank(card, 7, 'warp_tol', default=0.01)
        rad_check = integer_or_blank(card, 8, 'rad_check', default=3)
        #return VIEW3D(icavity, gitb=gitb, gips=gips, cier=cier,
                      #error_tol=error_tol, zero_tol=zero_tol, warp_tol=warp_tol,
                      #rad_check=rad_check, comment=comment)

        self.cards.append((icavity, gitb, gips, cier,
                           error_tol, zero_tol, warp_tol,
                           rad_check, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        icavity = np.zeros(ncards, dtype='int32')
        gitb = np.zeros(ncards, dtype='int32')
        gips = np.zeros(ncards, dtype='int32')
        cier = np.zeros(ncards, dtype='int32')
        error_tol = np.zeros(ncards, dtype='float64')
        zero_tol = np.zeros(ncards, dtype='float64')
        warp_tol = np.zeros(ncards, dtype='float64')
        rad_check = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (icavityi, gitbi, gipsi, cieri,
             error_toli, zero_toli, warp_toli, rad_checki, comment) = card
            #iview[icard] = iviewi
            icavity[icard] = icavityi
            gitb[icard] = gitbi
            gips[icard] = gipsi
            cier[icard] = cieri
            error_tol[icard] = error_toli
            zero_tol[icard] = zero_toli
            warp_tol[icard] = warp_toli
            rad_check[icard] = rad_checki
        self._save(icavity, gitb, gips, cier, error_tol, zero_tol, warp_tol, rad_check)
        self.cards = []

    def _save(self, icavity, gitb, gips, cier, error_tol, zero_tol, warp_tol, rad_check) -> int:
        assert len(self.icavity) == 0
        #self.iview = iview
        self.icavity = icavity
        self.gitb = gitb
        self.gips = gips
        self.cier = cier
        self.error_tol = error_tol
        self.zero_tol = zero_tol
        self.warp_tol = warp_tol
        self.rad_check = rad_check
        self.n = len(icavity)

    def __apply_slice__(self, view: VIEW, i: np.ndarray) -> None:
        view.icavity = self.icavity[i]
        view.gitb = self.gitb[i]
        view.gips = self.gips[i]
        view.cier = self.cier[i]
        view.error_tol = self.error_tol[i]
        view.zero_tol = self.zero_tol[i]
        view.warp_tol = self.warp_tol[i]
        view.rad_check = self.rad_check[i]
        view.n = len(i)

    @property
    def max_id(self) -> int:
        return self.icavity.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        #iviews = array_str(self.iview, size=size)
        icavitys = array_str(self.icavity, size=size)
        gitbs = array_str(self.gitb, size=size)
        gips = array_str(self.gips, size=size)
        ciers = array_str(self.cier, size=size)
        for icavity, gitb, gips, cier, error_tol, zero_tol, warp_tol, rad_check in zip_longest(
            icavitys, gitbs, gips, ciers,
            self.error_tol, self.zero_tol, self.warp_tol, self.rad_check):
            list_fields = ['VIEW3D', icavity, gitb, gips, cier,
                           error_tol, zero_tol, warp_tol, rad_check]
            bdf_file.write(print_card(list_fields))
        return
