"""
Defines the following contact cards:

nx contact:
 - BCONP
 - BLSEG
 - BCPARA
 - BCRPARA
 - BCTPARA
 - BCTADD
 - BCTSET
 - BSURF
 - BSURFS
 - BFRIC

glue:
 - BGADD
 - BGSET

"""
from __future__ import annotations
from typing import TYPE_CHECKING

# from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank,
    integer_double_or_blank, string,)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.contact.nx_contact import read_param_dict, params_listify
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class BGADD(BaseCard):
    """
    +-------+------+----+-------+----+----+----+----+----+
    |   1   |  2   | 3  |   4   |  5 | 6  |  7 | 8  |  9 |
    +=======+======+====+=======+====+====+====+====+====+
    | BGADD | GSID | SI |  S2   | S3 | S4 | S5 | S6 | S7 |
    +-------+------+----+-------+----+----+----+----+----+
    |       |  S8  | S9 |  etc. |    |    |    |    |    |
    +-------+------+----+-------+----+----+----+----+----+

    """
    type = 'BGADD'

    @classmethod
    def _init_from_empty(cls):
        glue_id = 1
        contact_sets = [1, 2]
        return BGADD(glue_id, contact_sets)

    def __init__(self, glue_id, contact_sets, comment: str=''):
        if comment:
            self.comment = comment
        #: Glue identification number. (Integer > 0)
        self.glue_id = glue_id

        #: Identification numbers of contact sets defined via BCTSET entries.
        #: (Integer > 0)
        self.contact_sets = contact_sets

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BGADD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        glue_id = integer(card, 1, 'glue_id')
        contact_sets = []

        i = 1
        j = 1
        while i < card.nfields:
            contact_set = integer(card, i, 'S%i' % j)
            contact_sets.append(contact_set)
            i += 1
            j += 1
        return BGADD(glue_id, contact_sets, comment=comment)

    def raw_fields(self):
        fields = ['BGADD'] + self.contact_sets
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class BGPARM(BaseCard):
    """

    +---------+--------+--------+--------+--------+--------+---------+--------+
    |    1    |   2    |    3   |   4    |   5    |   6    |    7    |    8   |
    +=========+========+========+========+========+========+=========+========+
    | BGPARM  | CSID   | Param1 | Value1 | Param2 | Value2 | Param3  | Value3 |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | Param4 | Value4 | Param5 | Value5 |  etc.  |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    | BGPARM  |   1    | PENN   |  10.0  |  PENT  |  0.5   |   CTOL  | 0.001  |
    +---------+--------+--------+--------+--------+--------+---------+--------+
    |         | SHLTHK |   1    |        |        |        |         |        |
    +---------+--------+--------+--------+--------+--------+---------+--------+

    """
    type = 'BGPARM'

    @classmethod
    def _init_from_empty(cls):
        csid = 1
        params = {'REFINE' : 0}
        return BGPARM(csid, params)

    def _finalize_hdf5(self, encoding: str):
        """hdf5 helper function"""
        keys, values = self.params
        self.params = {key : value for key, value in zip(keys, values)}

    def __init__(self, csid: int,
                 params: dict[str, int | float],
                 comment: str=''):
        """
        Creates a BCTPARM card

        Parameters
        ----------
        csid : int
            Contact set ID. Parameters defined in this command apply to
            contact set CSID defined by a BGPARM entry. (Integer > 0)
        params : dict[key] : value
            the optional parameters
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment

        #: Contact set ID. Parameters defined in this command apply to
        #: contact set CSID defined by a BGPARM entry. (Integer > 0)
        self.csid = csid
        self.params = params

    def _finalize_hdf5(self, encoding: str) -> None:
        keys = self.params[0]
        values = self.params[1]
        self.params = {}
        for key, value in zip(keys, values):
            if isinstance(value, bytes):
                value = value.decode(encoding)
            self.params[key] = value

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BGPARM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        csid = integer(card, 1, 'csid')
        params = read_param_dict(card, cls._get_param_value_from_card)
        return BGPARM(csid, params, comment=comment)

    def _get_param_value_from_card(
            card: BDFCard,
            i: int, j: int) -> tuple[str, int | float]:
        param = string(card, i, 'param%s' % j)
        if param == 'REFINE':
            value = integer_or_blank(card, i+1, f'value{j:d}', default=0)
            assert value in [0], 'REFINE must be [0]; REFINE=%r' % value
        else:
            # FRICMOD, FPARA1/2/3/4/5, EPSN, EPST, CFACTOR1, PENETOL
            # NCMOD, TCMOD, RFORCE, LFORCE, RTPCHECK, RTPMAX, XTYPE
            # ...
            value = integer_double_or_blank(card, i+1, 'value%s' % j)
            assert value is not None, '%s%i must not be None' % (param, j)
        return param, value

    def raw_fields(self) -> list[str | int | float]:
        fields = ['BGPARM', self.csid]
        fields.extend(params_listify(self.params))
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class BGSET(BaseCard):
    """
    +-------+------+------+------+---------+----+------+------+----+
    |   1   |  2   |  3   |   4  |    5    | 6  |  7   |   8  |  9 |
    +=======+======+======+======+=========+====+======+======+====+
    | BGSET | GSID | SID1 | TID1 | SDIST1  |    | EXT1 |      |    |
    +-------+------+------+------+---------+----+------+------+----+
    |       |      | SID2 | TID2 | SDIST2  |    | EXT2 |      |    |
    +-------+------+------+------+---------+----+------+------+----+
    """
    type = 'BGSET'

    @classmethod
    def _init_from_empty(cls):
        glue_id = 1
        sids = [1]
        tids = [1]
        sdists = [0.01]
        exts = [1.]
        return BGSET(glue_id, sids, tids, sdists, exts, sol=101)

    def __init__(self, glue_id: int,
                 sids: list[int],
                 tids: list[int],
                 sdists: list[float],
                 exts: list[float],
                 comment: str='', sol: int=101):
        if comment:
            self.comment = comment
        #: GSID Glue set identification number. (Integer > 0)
        self.glue_id = glue_id
        #: SIDi Source region (contactor) identification number for contact pair i.
        #: (Integer > 0)
        self.sids = sids

        #: TIDi Target region identification number for contact pair i. (Integer > 0)
        self.tids = tids

        #: SDISTi Search distance for glue regions (Real); (Default=10.0)
        self.sdists = sdists

        #: EXTi Extension factor for target region (SOLs 402 and 601 only).
        self.exts = exts

    @classmethod
    def add_card(cls, card: BDFCard, comment: str='', sol: int=101):
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
        return BGSET(glue_id, sids, tids, sdists, exts,
                     comment=comment, sol=sol)

    def raw_fields(self):
        fields = ['BGSET', self.glue_id]
        for sid, tid, sdist, ext in zip(self.sids, self.tids, self.sdists,
                                        self.exts):
            fields += [sid, tid, sdist, None, ext, None, None, None]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
