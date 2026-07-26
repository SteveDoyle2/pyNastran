from __future__ import annotations
from itertools import count, zip_longest
from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank,
    double, double_or_blank,
    string, string_or_blank,)

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, make_idim, hslice_by_idim,
    parse_check, save_ifile_comment)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_float,
    array_default_int, get_print_card_size)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


class CSSCHD(VectorizedBaseCard):
    """
    Defines a scheduled control surface deflection as a function of
    Mach number and angle of attack.

    +--------+-----+-------+--------+-------+-------+
    |    1   |  2  |   3   |   4    |   5   |   6   |
    +========+=====+=======+========+=======+=======+
    | CSSCHD | SlD | AESID | LALPHA | LMACH | LSCHD |
    +--------+-----+-------+--------+-------+-------+
    | CSSCHD |  5  |  50   |   12   |   15  |   25  |
    +--------+-----+-------+--------+-------+-------+
    """
    _id_name = 'csschd_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.csschd_id = np.array([], dtype='int32')
        self.aesurf_id = np.array([], dtype='int32')
        self.lalpha = np.array([], dtype='int32')
        self.lmach = np.array([], dtype='int32')
        self.lschd = np.array([], dtype='int32')

    def add(self, sid: int, aesurf_id: int,
            lschd: int, lalpha: int=None, lmach: int=None,  # aefact
            ifile: int=0, comment: str='') -> int:
        """
        Creates an CSSCHD card, which defines a specified control surface
        deflection as a function of Mach and alpha (used in SOL 144/146).

        Parameters
        ----------
        sid : int
            the unique id
        aesid : int
            the control surface (AESURF) id
        lalpha : int; default=None
            the angle of attack profile (AEFACT) id
        lmach : int; default=None
            the mach profile (AEFACT) id
        lschd : int; default=None
            the control surface deflection profile (AEFACT) id
        comment : str; default=''
            a comment for the card

        """
        #assert lalpha is None or isinstance(lalpha, integer_types), lalpha
        #assert lmach is None or isinstance(lmach, integer_types), lmach
        #assert lschd is None or isinstance(lschd, integer_types), lschd
        card = (sid, aesurf_id, lalpha, lmach, lschd, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a CSSCHD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        aesurf_id = integer(card, 2, 'aesid')         # AESURF
        lalpha = integer_or_blank(card, 3, 'lAlpha')  # AEFACT
        lmach = integer_or_blank(card, 4, 'lMach')    # AEFACT
        lschd = integer(card, 5, 'lSchd')             # AEFACT
        assert len(card) <= 6, f'len(CSSCHD card) = {len(card):d}\ncard={card}'
        #return CSSCHD(sid, aesurf_id, lalpha, lmach, lschd, comment=comment)
        card = (sid, aesurf_id, lalpha, lmach, lschd, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        csschd_id = np.zeros(ncards, dtype='int32')
        aesurf_id = np.zeros(ncards, dtype='int32')
        lalpha = np.zeros(ncards, dtype='int32')
        lmach = np.zeros(ncards, dtype='int32')
        lschd = np.zeros(ncards, dtype='int32')
        comment = {}

        for icard, card in enumerate(self.cards):
            (sid, aesurf_idi, lalphai, lmachi, lschdi, ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[sid] = commenti
            csschd_id[icard] = sid
            aesurf_id[icard] = aesurf_idi
            lalpha[icard] = lalphai
            lmach[icard] = lmachi
            lschd[icard] = lschdi
        self._save(csschd_id, aesurf_id, lalpha, lmach, lschd,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, csschd_id, aesurf_id, lalpha, lmach, lschd,
              ifile=None, comment=None):
        ncards = len(csschd_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.csschd_id):
            ifile = np.hstack([self.ifile, ifile])
            csschd_id = np.hstack([self.csschd_id, csschd_id])
            aesurf_id = np.hstack([self.aesurf_id, aesurf_id])
            lalpha = np.hstack([self.lalpha, lalpha])
            lmach = np.hstack([self.lmach, lmach])
            lschd = np.hstack([self.lschd, lschd])
        save_ifile_comment(self, ifile, comment)
        self.csschd_id = csschd_id
        self.aesurf_id = aesurf_id
        self.lalpha = lalpha
        self.lmach = lmach
        self.lschd = lschd
        self.n = len(csschd_id)

    #def sort(self) -> None:
        #ueid = np.unique(self.csschd_id)
        #if np.array_equal(ueid, self.csschd_id):
            #return
        #i = np.argsort(self.csschd_id)
        #self.__apply_slice__(self, i)

    def __apply_slice__(self, load: CSSCHD, i: np.ndarray) -> None:
        load.n = len(i)
        load.ifile = self.ifile[i]
        load.csschd_id = self.csschd_id[i]
        load.aesurf_id = self.aesurf_id[i]
        load.lalpha = self.lalpha[i]
        load.lmach = self.lmach[i]
        load.lschd = self.lschd[i]

    def validate(self):
        #if not(self.lalpha is None or isinstance(self.lalpha, integer_types)):
            #raise TypeError('lalpha=%r must be an int or None' % self.lalpha)
        #if not(self.lmach is None or isinstance(self.lmach, integer_types)):
            #raise TypeError('lmach=%r must be an int or None' % self.lmach)
        nalpha = self.lalpha.size
        nmach = self.lmach.size
        if (nalpha == 0) or (nalpha != nmach):
            raise RuntimeError(f'CSSCHD csschd_id=%s; nalpha={nalpha} nmach={nmach}')

        ibad = (self.lalpha == 0) and (self.lmach == 0)
        if np.any(ibad):
            csschd = self.csschd_id[ibad]
            card = self.slice_card_by_index(ibad)
            msgi = card.write()
            msg = ('CSSCHD csschd_id=%s; lalpha and lmach are both 0'
                   ' (one must be an integer (AEFACT)\n%s' % (csschd, msgi))
            raise RuntimeError(msg)

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        uaesurf = np.unique(self.aesurf_id)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self, missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    @property
    def max_id(self) -> int:
        return max(self.csschd_id.max(), self.aesurf_id.max(),
                   self.lalpha.max(), self.lmach.max(), self.lschd.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        csschd_id_ = array_str(self.csschd_id, size=size)
        aesurf_id_ = array_str(self.aesurf_id, size=size)
        lalpha_ = array_default_int(self.lalpha, default=0, size=size)
        lmach_ = array_default_int(self.lmach, default=0, size=size)
        lschd_ = array_str(self.lschd, size=size)

        for csschd_id, aesurf_id, lalpha, lmach, lschd in zip(
                csschd_id_, aesurf_id_, lalpha_, lmach_, lschd_):

            list_fields = ['CSSCHD', csschd_id, aesurf_id, lalpha, lmach, lschd]
            bdf_file.write(print_card(list_fields))
        return


class DIVERG(VectorizedBaseCard):
    """
    +--------+-----+--------+----+----+----+----+----+----+
    |   1    |  2  |   3    | 4  | 5  | 6  | 7  | 8  | 9  |
    +========+=====+========+====+====+====+====+====+====+
    | DIVERG | SID | NROOT  | M1 | M2 | M3 | M4 | M5 | M6 |
    +--------+-----+--------+----+----+----+----+----+----+
    |        |  M7 |  etc.  |    |    |    |    |    |    |
    +--------+-----+--------+----+----+----+----+----+----+

    Attributes
    ----------
    sid : int
        The name.
    nroots : int
        the number of roots
    machs : list[float, ..., float]
        list of Mach numbers

    """
    _id_name = 'diverg_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.diverg_id = np.array([], dtype='int32')
        #self.aesurf_id = np.array([], dtype='int32')
        #self.lalpha = np.array([], dtype='int32')
        #self.lmach = np.array([], dtype='int32')
        #self.lschd = np.array([], dtype='int32')

    def add(self, diverg_id: int, nroots: int,
            machs: list[float],
            ifile: int=0, comment: str='') -> int:
        """
        Creates an DIVERG card, which is used in divergence
        analysis (SOL 144).

        Parameters
        ----------
        diverg_id : int
            The name
        nroots : int
            the number of roots
        machs : list[float, ..., float]
            list of Mach numbers
        comment : str; default=''
            a comment for the card

        """
        #assert lalpha is None or isinstance(lalpha, integer_types), lalpha
        #assert lmach is None or isinstance(lmach, integer_types), lmach
        #assert lschd is None or isinstance(lschd, integer_types), lschd
        card = (diverg_id, nroots, machs, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a DIVERG card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        nroots = integer(card, 2, 'nroot')
        j = 1
        machs = []
        for i in range(3, len(card)):
            mach = double(card, i, 'Mach_%i' % j)
            machs.append(mach)
            j += 1
        assert len(machs) > 0, card
        #return DIVERG(sid, nroots, machs, comment=comment)
        self.cards.append((sid, nroots, machs, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        diverg_id = np.zeros(ncards, dtype='int32')
        nroots = np.zeros(ncards, dtype='int32')
        nmach = np.zeros(ncards, dtype='int32')
        comment = {}

        machs_list = []
        for icard, card in enumerate(self.cards):
            (sid, nrootsi, machsi, ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[sid] = commenti
            diverg_id[icard] = sid
            nroots[icard] = nrootsi
            nmach[icard] = len(machsi)
            machs_list.extend(machsi)
        machs = np.array(machs_list, dtype='float64')
        self._save(diverg_id, nroots, machs, nmach,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, diverg_id, nroots, machs, nmach,
              ifile=None, comment=None):
        ncards = len(diverg_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.diverg_id):
            ifile = np.stack([self.ifile, ifile])
            asfd
        save_ifile_comment(self, ifile, comment)
        self.diverg_id = diverg_id
        self.nroots = nroots
        self.machs = machs
        self.nmach = nmach
        self.n = len(diverg_id)

    #def sort(self) -> None:
        #ueid = np.unique(self.csschd_id)
        #if np.array_equal(ueid, self.csschd_id):
            #return
        #i = np.argsort(self.csschd_id)
        #self.__apply_slice__(self, i)

    def __apply_slice__(self, load: DIVERG, i: np.ndarray) -> None:
        load.n = len(i)
        load.diverg_id = self.diverg_id[i]
        load.nroots = self.nroots[i]
        imach = self.imach
        load.machs = hslice_by_idim(i, imach, self.machs)
        load.nmach = self.nmach[i]

    @property
    def imach(self) -> np.ndarray:
        return make_idim(self.n, self.nmach)

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        #uaesurf = np.unique(self.aesurf_id)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self, missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    @property
    def max_id(self) -> int:
        return self.diverg_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        diverg_ids = array_str(self.diverg_id, size=size)
        nroots = array_str(self.nroots, size=size)
        machs = array_float(self.machs, size=size, is_double=False).tolist()

        for diverg_id, nroot, (imach0, imach1) in zip(diverg_ids, nroots, self.imach):
            mach = machs[imach0:imach1]
            list_fields = ['DIVERG', diverg_id, nroot] + mach
            bdf_file.write(print_card(list_fields))
        return


class TRIM(VectorizedBaseCard):
    """
    Specifies constraints for aeroelastic trim variables.

    +------+--------+------+--------+--------+-----+--------+-----+----------+
    |   1  |   2    |   3  |    4   |    5   |  6  |    7   |  8  |     9    |
    +======+========+======+========+========+=====+========+=====+==========+
    | TRIM |   ID   | MACH |    Q   | LABEL1 | UX1 | LABEL2 | UX2 | IS_RIGID |
    +------+--------+------+--------+--------+-----+--------+-----+----------+
    |      | LABEL3 |  UX3 | LABEL4 |   UX4  | ... |        |     |          |
    +------+--------+------+--------+--------+-----+--------+-----+----------+
    """
    _id_name = 'trim_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.trim_id = np.array([], dtype='int32')
        self.mach = np.array([], dtype='float64')
        self.q = np.array([], dtype='float64')
        self.aeqr = np.array([], dtype='float64')
        self.nlabel = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U8')
        self.ux = np.array([], dtype='float64')

    def add(self, sid: int, mach: float, q: float,
            labels: list[str], uxs: list[float], aeqr: float=1.0,
            trim_type: int=1,
            ifile: int=0, comment: str='') -> int:
        """
        Creates a TRIM/TRIM2 card for a static aero (144) analysis.

        Parameters
        ----------
        sid : int
            the trim id; referenced by the Case Control TRIM field
        mach : float
            the mach number
        q : float
            dynamic pressure
        labels : list[str]
            names of the fixed variables
        uxs : list[float]
            values corresponding to labels
        aeqr : float
            0.0 : rigid trim analysis
            1.0 : elastic trim analysis (default)
        trim_type : int
            1 : creates a TRIM
            2 : creates a TRIM2
        comment : str; default=''
            a comment for the card

        """
        assert len(labels) == len(uxs)
        card = (sid, mach, q, labels, uxs, aeqr, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        mach = double(card, 2, 'mach')
        q = double(card, 3, 'q')
        labels = []
        uxs = []

        label = string_or_blank(card, 4, 'label1')
        if label:
            ux = double(card, 5, 'ux1')
            uxs.append(ux)
            labels.append(label)

        label = string_or_blank(card, 6, 'label2')
        if label:
            ux = double(card, 7, 'ux1')
            uxs.append(ux)
            labels.append(label)
        aeqr = double_or_blank(card, 8, 'aeqr', default=1.0)

        i = 9
        n = 3
        while i < len(card):
            label = string(card, i, 'label%i' % n)
            ux = double(card, i + 1, 'ux%i' % n)
            labels.append(label)
            uxs.append(ux)
            i += 2
            n += 1
        #return TRIM(sid, mach, q, labels, uxs, aeqr, comment=comment)
        self.cards.append((sid, mach, q, labels, uxs, aeqr, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card2(self, card: BDFCard, ifile: int,
                  comment: str='') -> int:
        """
        Adds a TRIM2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        mach = double(card, 2, 'mach')
        q = double(card, 3, 'q')
        aeqr = double_or_blank(card, 8, 'aeqr', default=1.0)

        i = 9
        n = 3
        labels = []
        uxs = []
        while i < len(card):
            label = string(card, i, 'label%d' % n)
            ux = double(card, i + 1, 'ux%d' % n)
            labels.append(label)
            uxs.append(ux)
            i += 2
        #return TRIM2(sid, mach, q, labels, uxs, aeqr, comment=comment)
        self.cards.append((2, sid, mach, q, labels, uxs, aeqr, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        trim_id = np.zeros(ncards, dtype='int32')
        mach = np.zeros(ncards, dtype='float64')
        q = np.zeros(ncards, dtype='float64')
        aeqr = np.zeros(ncards, dtype='float64')
        nlabel = np.zeros(ncards, dtype='int32')
        label = []
        ux = []
        comment = {}
        for icard, card in enumerate(self.cards):
            (sid, machi, qi, labelsi, uxsi, aeqri, ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[sid] = commenti
            trim_id[icard] = sid
            mach[icard] = machi
            q[icard] = qi
            aeqr[icard] = aeqri
            nlabel[icard] = len(labelsi)
            label.extend(labelsi)
            ux.extend(uxsi)
        label = np.array(label, dtype='|U8')
        ux = np.array(ux, dtype='float64')
        self._save(trim_id, mach, q, aeqr, nlabel, label, ux,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, trim_id, mach, q, aeqr, nlabel, label, ux,
              ifile=None, comment=None):
        ncards = len(trim_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.trim_id):
            ifile = np.stack([self.ifile, ifile])
            afd
        save_ifile_comment(self, ifile, comment)
        self.trim_id = trim_id
        self.mach = mach
        self.q = q
        self.aeqr = aeqr
        self.nlabel = nlabel
        self.label = label
        self.ux = ux
        self.n = len(trim_id)

    def validate(self):
        assert self.mach.min() >= 0.0, 'mach = %s' % self.mach
        assert np.all(self.mach != 1.0), 'mach = %s' % self.mach
        assert self.q.min() > 0.0, 'q=%s' % self.q

        for trim_id, (ilabel0, ilabel1) in zip(self.trim_id, self.ilabel):
            labels = self.label[ilabel0:ilabel1]
            if len(set(labels)) != len(labels):
                msg = 'TRIM id=%d; not all labels are unique; labels=%s' % (trim_id, str(labels))
                raise RuntimeError(msg)
        #if len(self.labels) != len(self.uxs):
            #msg = 'nlabels=%s != nux=%s; labels=%s uxs=%s' % (
                #len(self.labels), len(self.uxs), str(self.labels), str(self.uxs))
            #raise RuntimeError(msg)

    @property
    def ilabel(self) -> np.ndarray:
        return make_idim(self.n, self.nlabel)

    def __apply_slice__(self, load: TRIM, i: np.ndarray) -> None:
        load.n = len(i)
        load.trim_id = self.trim_id[i]
        load.mach = self.mach[i]
        load.q = self.q[i]
        load.aeqr = self.aeqr[i]

        ilabel = self.ilabel
        load.label = hslice_by_idim(i, ilabel, self.label)
        load.ux = hslice_by_idim(i, ilabel, self.ux)
        load.nlabel = self.nlabel[i]

    def convert(self, pressure_scale: float=1.0, **kwargs) -> None:
        self.q *= pressure_scale

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        uaesurf = np.unique(self.trim_id)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    def verify_trim(self, suport1_id: int) -> None:
        """
        Magic function that makes TRIM cards not frustrating.

        .. warning ::  This probably gets AELINKs/AEPARMs/AESURFSs wrong.

        **The TRIM equality**
        ndelta = (naestat + naesurf + naeparm)
               - (ntrim + ntrim_aesurf? + naelink + nsuport_dofs + nsuport1_dofs)
        ndelta = 0
        ntrim_aesurf is not included, but it might exist...

        **Steps to a TRIM analysis**
        1.  Define the number of independent control surfaces (naesurf)
            Make an AESURF for each.  Dual link the AESURFs if you can
            to avoid needing an AELINK (e.g., +roll is left aileron down,
            right aileron up).
            Horizontal Tail : name it DPITCH
            Vertical Tail   : name it DYAW
            Aileron         : name it DROLL
        2.  Create AELINKs if necessary.
        3.  Add the AESTAT variables.  Include one for each DOF the
            aircraft can move in the frame of the model
            (e.g., half/full model).
                Half model (2.5g pullup, abrupt pitch):
                  - 2d pitch/plunge, 1 control : URDD3, URDD5, PITCH, ANGLEA
                Full model (2.5g pullup, abrupt pitch):
                  - 3d pitch/plunge, 3 control : URDD3, URDD5, PITCH, ANGLEA, YAW (???)
        4.  Add the TRIM card to lock the variables that could theoretically move
            in the plane of the analysis that are known.
                Half model:
                   2.5g pullup   : lock URDD3=2.5, URDD5=0, PITCH=0
                                   solve for ANGLEA, DPITCH
                                   use DPITCH
                   abrupt pitch  : lock URDD3=1.0, URDD5=0, ANGLEA=5
                                   solve for PITCH, DPITCH
                                   use DPITCH
                Full model:
                   2.5g pullup   : lock URDD3=2.5, URDD4=0, URDD5=0,  PITCH=0, YAW=0,
                                   lock SIDES=0,  ROLL=0
                                   solve for ANGLEA, DPITCH
                                   use DPITCH, DYAW, DROLL
                                   TODO: probably wrong
                   30 degree yaw : lock URDD3=1.0, URDD4=0, ANGLEA=5, PITCH=0, YAW=30,
                                   lock DPITCH=0, ROLL=0
                                   solve for SIDES, URDD5
                                   use DPITCH, DYAW, DROLL
                                   TODO: probably wrong

        5.  Note that we could have simplified our full model AESTAT/TRIM
            cards (they can be the same as for a half model), but we'd
            like to be able to do multiple load cases in the same deck.

        6.  Add some SUPORT/SUPORT1 DOFs to ignore non-relevant motion in
            certain DOFs (e.g., z-motion).  Add enough to satisfy the TRIM
            equality.

        **Doesn't Consider**
         - AELINK
         - AEPARM
         - AESURFS

        +------------------------------------------------+
        |                 Default AESTATs                |
        +--------+---------+-----------------------------+
        | ANGLEA | ur (R2) | Angle of Attack             |
        | YAW    | ur (R3) | Yaw Rate                    |
        | SIDES  | ur (R3) | Angle of Sideslip           |
        +--------+---------+-----------------------------+
        | ROLL   | ůr (R1) | Roll Rate                   |
        | PITCH  | ůr (R2) | Pitch Rate                  |
        +--------+---------+-----------------------------+
        | URDD1  | ür (T1) | Longitudinal (See Remark 3) |
        | URDD2  | ür (T2) | Lateral                     |
        | URDD3  | ür (T3) | Vertical                    |
        | URDD4  | ür (R1) | Roll                        |
        | URDD5  | ür (R2) | Pitch                       |
        | URDD6  | ür (R3) | Yaw                         |
        +--------+---------+-----------------------------+
        """
        #if not xref:
            #return
        #suport, suport1, aestats, aeparms, aelinks, aesurf, xref=True
        #suport = []
        #if 'SUPORT1' in subcase:
            #suport_id = subcase.get_int_parameter('SUPORT1')
            #suport1 = fem.suport1[suport_id]

        suport = self.model.suport

        suport_ids = np.unique(suport.suport_id)
        nsuport_dofs = 0
        nsuport1_dofs = 0
        suport_dofs = set()
        suport1_dofs = set()

        #for (inode0, inode1) in suport.inode:
            #nodes = suport.node[inode0:inode1]
            #components = suport.component[inode0:inode1].astype('|U8')
            #for nid, cs in zip(nodes, components):
                #for ci in str(cs):
                    ##print('  SUPORT: nid=%r C=%r' % (nid, ci))
                    #dof = (nid, ci)
                    #if dof in suport_dofs:
                        #msg = 'Duplicate DOF\n  dof=%s suport_dofs=%s' % (
                            #str(dof), str(suport_dofs))
                        #raise RuntimeError(msg)
                    #suport_dofs.add(dof)
                    #nsuport_dofs += 1

        suport_dof_msg2 = ''
        if 0 in suport_ids:
            suport_id = 0
            suport0 = suport.slice_card_by_id(suport_id) # , assume_sorted=True, sort_ids=False
            suport_dof_msg = ''
            for nid, component in zip(suport0.node_id, suport0.component):
                for componenti in str(component):
                    dof = (nid, componenti)
                    suport_dof_msg += '    (%s, %s)\n' % (nid, componenti)
                    if dof in suport_dofs:
                        msg = 'dof=%s suport_dofs=%s' % (str(dof), str(suport_dofs))
                        raise RuntimeError(msg)
                    suport_dofs.add(dof)
                    nsuport_dofs += 1
            suport_dof_msg2 = '\nsuport_dofs (nid, comp):\n%s\n' % suport_dof_msg.rstrip(',')

        if suport1_id > 0:
            suport1 = suport.slice_card_by_id(suport1_id) # , assume_sorted=True, sort_ids=False
            suport1_dof_msg = ''
            for nid, component in zip(suport1.node_id, suport1.component):
                for componenti in str(component):
                    dof = (nid, componenti)
                    suport1_dof_msg += '    (%s, %s)\n' % (nid, componenti)
                    if dof in suport1_dofs:
                        msg = 'dof=%s suport1_dofs=%s' % (str(dof), str(suport_dofs))
                        raise RuntimeError(msg)
                    suport_dofs.add(dof)
                    nsuport1_dofs += 1
            suport_dof_msg2 = '\nsuport_dofs (nid, comp):\n%s\n' % suport1_dof_msg.rstrip(',')

        aesurf_names = self.model.aesurf.label.tolist()
        aestat_labels = self.model.aestat.label.tolist()
        aeparm_labels = self.model.aeparm.label.tolist()
        print(self.get_stats())
        #aesurf_names = [aesurfi.label for aesurfi in aesurf.values()]
        #aestat_labels = [aestat.label for aestat in aestats.values()]
        #aeparm_labels = [aeparm.label for aeparm in aeparms.values()]
        naestat = len(aestat_labels)
        ntrim = len(self.label)
        trim_aesurf_common = list(set(self.label).intersection(set(aesurf_names)))
        trim_aesurf_common.sort()
        ntrim_aesurfs = len(trim_aesurf_common)
        naesurf = len(aesurf_names)
        naeparm = len(aeparm_labels)

        aelinksi = []
        trim_id = self.trim_id[0]
        (ilabel0, ilabel1) = self.ilabel[0]
        labels = self.label[ilabel0:ilabel1]

        i0 = np.where((self.model.aelink.aelink_id == 0) |
                       self.model.aelink.aelink_id == trim_id)[0]
        if len(i0):
            aelinksi = self.model.aelink.slice_card_by_index(i0)
        #if 'ALWAYS' in aelinks:
            #aelinksi += [aelink.label for aelink in aelinks['ALWAYS']]

        naelink = len(aelinksi)


        ntrim_aesurf = 0
        allowed_labels = aestat_labels + aesurf_names + aeparm_labels
        msg = ''
        for label in labels:
            if label not in allowed_labels:
                msg += 'TRIM label=%r is not defined\n' % label

            if label in aesurf_names:
                #print('AESTAT/AESURF label = %r' % label)
                ntrim_aesurf += 1
        if msg:
            msg += '\n aestat_labels=%s\n aeparm_labels=%s\n aesurf_names=%s\n%s' % (
                aestat_labels, aeparm_labels, aesurf_names, str(self))
            raise RuntimeError(msg)

        # TODO: this doesn't work for multiple subcases
        #ntotal_suport_dofs = nsuport_dofs, nsuport1_dofs
        #ndelta = ntrim - nsuport_dofs - nsuport1_dofs - naesurf
        #if ndelta != 0:
            #msg = 'ntrim - nsuport_dofs - nsuport1_dofs - naesurf = ndelta = %s; ndelta != 0\n' % ndelta
            #msg += 'ntrim=%s nsuport_dofs=%s nsuport1_dofs=%s naesurfs=%s' % (
                #ntrim, nsuport_dofs, nsuport1_dofs, naesurf)
            #raise RuntimeError(msg)

        #ndelta = (naestat + naesurf + naeparm + ntrim_aesurf) - (ntrim + naelink + nsuport_dofs + nsuport1_dofs)
        #if ndelta != 0:
            #msg = (
                #'(naestat + naesurf + naeparm + ntrim_aesurf) - '
                #'(ntrim + naelink + nsuport_dofs + nsuport1_dofs) = ndelta = %s; ndelta != 0\n'
                #'naestat=%s naesurf=%s naeparm=%s ntrim_aesurfs=%s\n'
                #'ntrim=%s naelink=%s nsuport_dofs=%s nsuport1_dofs=%s' % (
                    #ndelta,
                    #naestat, naesurf, naeparms, ntrim_aesurf,
                    #ntrim, naelink, nsuport_dofs, nsuport1_dofs))

        nplus = (naestat + naesurf + naeparm)
        nminus = ntrim + naelink + nsuport_dofs + nsuport1_dofs

        ndelta = nplus - nminus + 0*2*ntrim_aesurfs
        if ndelta != 0:
            #msg = (
                #'(naestat + naesurf + naeparm) - (ntrim + ntrim_aesurf? + naelink + '
                #'nsuport_dofs + nsuport1_dofs) = ndelta = %s; ndelta != 0\n'
                #'naestat=%s naesurf=%s naeparm=%s ntrim=%s ntrim_aesurf=%s '
                #'naelink=%s nsuport_dofs=%s nsuport1_dofs=%s\n' % (
                    #ndelta,
                    #naestat, naesurf, naeparm, ntrim, ntrim_aesurf,
                    #naelink, nsuport_dofs, nsuport1_dofs)
            #)
            msg = (
                'Invalid trim state (ndelta != 0):\n'
                f'   (naestat + naesurf + naeparm + 0*2*ntrim_aesurf?) = ({naestat} + {naesurf} + {naeparm} + 0*2*{ntrim_aesurf}) = {nplus}\n'
                f' - (ntrim + naelink + nsuport_dofs + nsuport1_dofs) = ({ntrim} + {naelink} + {nsuport_dofs} + {nsuport1_dofs}) = {nminus}\n'
                '===================================================================\n'
                f'  ndelta = {ndelta}\n\n'
                'Summary\n'
                '-------\n'
                f'  +naestat = {naestat}; {aestat_labels}\n'
                f'  +naesurf = {naesurf}; {aesurf_names}\n'
                f'  +naeparm = {naeparm}; {aeparm_labels}\n'
                f'  +0*2*ntrim_aesurf? = {2*ntrim_aesurf} -> 0; {trim_aesurf_common}\n'
                f'  -ntrim = {ntrim}; {labels}\n'
                f'  -naelink = {naelink}; {aelinksi}\n'
                f'  -nsuport_dofs = {nsuport_dofs}\n'
                f'  -nsuport1_dofs = {nsuport1_dofs}\n'
                f'{suport_dof_msg2}\n\n'
            )
            msg += str(self)
            raise RuntimeError(msg)

    @property
    def max_id(self) -> int:
        return self.trim_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        trim_ids = array_str(self.trim_id, size=size)

        for trim_id, mach, q, aeqr, (ilabel0, ilabel1) in zip_longest(
                trim_ids, self.mach, self.q, self.aeqr, self.ilabel):

            labels = self.label[ilabel0:ilabel1]
            uxs = self.ux[ilabel0:ilabel1]
            list_fields = ['TRIM', trim_id, mach, q]
            nlabels = len(labels)
            assert nlabels > 0, labels
            for (i, label, ux) in zip(count(), labels, uxs):
                list_fields += [label, ux]
                if i == 1:
                    list_fields += [aeqr]
            if nlabels == 1:
                list_fields += [None, None, aeqr]
            bdf_file.write(print_card(list_fields))
        return


class TRIM2(VectorizedBaseCard):
    """
    Defines the state of the aerodynamic extra points for a trim analysis.
    All undefined extra points will be set to zero.

    +-------+--------+------+--------+--------+-----+--------+-----+----------+
    |   1   |   2    |   3  |    4   |    5   |  6  |    7   |  8  |     9    |
    +=======+========+======+========+========+=====+========+=====+==========+
    | TRIM2 |   ID   | MACH |    Q   |        |     |        |     | IS_RIGID |
    +-------+--------+------+--------+--------+-----+--------+-----+----------+
    |       | LABEL1 |  UX1 | LABEL2 |   UX2  | ... |        |     |          |
    +-------+--------+------+--------+--------+-----+--------+-----+----------+
    """
    _id_name = 'trim_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.trim_id = np.array([], dtype='int32')
        self.mach = np.array([], dtype='float64')
        self.q = np.array([], dtype='float64')
        self.aeqr = np.array([], dtype='float64')
        self.nlabel = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U8')
        self.ux = np.array([], dtype='float64')

    def add(self, sid: int, mach: float, q: float,
            labels: list[str], uxs: list[float], aeqr: float=1.0,
            ifile: int=0, comment: str='') -> int:
        assert len(labels) == len(uxs)
        card = (sid, mach, q, labels, uxs, aeqr, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        mach = double(card, 2, 'mach')
        q = double(card, 3, 'q')
        aeqr = double_or_blank(card, 8, 'aeqr', default=1.0)

        i = 9
        n = 1
        labels = []
        uxs = []
        while i < len(card):
            label = string(card, i, 'label%d' % n)
            ux = double(card, i + 1, 'ux%d' % n)
            labels.append(label)
            uxs.append(ux)
            i += 2
            n += 1
        self.cards.append((sid, mach, q, labels, uxs, aeqr, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        trim_id = np.zeros(ncards, dtype='int32')
        mach = np.zeros(ncards, dtype='float64')
        q = np.zeros(ncards, dtype='float64')
        aeqr = np.zeros(ncards, dtype='float64')
        nlabel = np.zeros(ncards, dtype='int32')
        label = []
        ux = []
        comment = {}
        for icard, card in enumerate(self.cards):
            (sid, machi, qi, labelsi, uxsi, aeqri, ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[sid] = commenti
            trim_id[icard] = sid
            mach[icard] = machi
            q[icard] = qi
            aeqr[icard] = aeqri
            nlabel[icard] = len(labelsi)
            label.extend(labelsi)
            ux.extend(uxsi)
        label = np.array(label, dtype='|U8')
        ux = np.array(ux, dtype='float64')
        self._save(trim_id, mach, q, aeqr, nlabel, label, ux,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, trim_id, mach, q, aeqr, nlabel, label, ux,
              ifile=None, comment=None):
        ncards = len(trim_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.trim_id):
            trim_id = np.hstack([self.trim_id, trim_id])
            mach = np.hstack([self.mach, mach])
            q = np.hstack([self.q, q])
            aeqr = np.hstack([self.aeqr, aeqr])
            nlabel = np.hstack([self.nlabel, nlabel])
            label = np.hstack([self.label, label])
            ux = np.hstack([self.ux, ux])
        save_ifile_comment(self, ifile, comment)
        self.trim_id = trim_id
        self.mach = mach
        self.q = q
        self.aeqr = aeqr
        self.nlabel = nlabel
        self.label = label
        self.ux = ux
        self.n = len(trim_id)

    @property
    def ilabel(self) -> np.ndarray:
        return make_idim(self.n, self.nlabel)

    def __apply_slice__(self, load: TRIM2, i: np.ndarray) -> None:
        load.n = len(i)
        load.trim_id = self.trim_id[i]
        load.mach = self.mach[i]
        load.q = self.q[i]
        load.aeqr = self.aeqr[i]
        ilabel = self.ilabel
        load.label = hslice_by_idim(i, ilabel, self.label)
        load.ux = hslice_by_idim(i, ilabel, self.ux)
        load.nlabel = self.nlabel[i]

    def convert(self, pressure_scale: float=1.0, **kwargs) -> None:
        self.q *= pressure_scale

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        geom_check(self, missing)

    @property
    def max_id(self) -> int:
        return self.trim_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        trim_ids = array_str(self.trim_id, size=size)

        for trim_id, mach, q, aeqr, (ilabel0, ilabel1) in zip_longest(
                trim_ids, self.mach, self.q, self.aeqr, self.ilabel):

            labels = self.label[ilabel0:ilabel1]
            uxs = self.ux[ilabel0:ilabel1]
            list_fields = ['TRIM2', trim_id, mach, q, None, None, None, None, aeqr]
            nlabels = len(labels)
            assert nlabels > 0, labels
            for label, ux in zip(labels, uxs):
                list_fields += [label, ux]
            bdf_file.write(print_card(list_fields))
        return


class AEFORCE(VectorizedBaseCard):
    """
    Defines an aerodynamic force vector for static aeroelastic trim.

    +---------+------+-------+-------+------+------+-------+------+------+
    |    1    |   2  |   3   |   4   |   5  |   6  |   7   |   8  |   9  |
    +=========+======+=======+=======+======+======+=======+======+======+
    | AEFORCE | MACH | SYMXZ | SYMXY | UXID | MESH | FORCE | DMIK | PERQ |
    +---------+------+-------+-------+------+------+-------+------+------+
    """
    _id_name = 'ux_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.mach = np.array([], dtype='float64')
        self.sym_xz = np.array([], dtype='|U8')
        self.sym_xy = np.array([], dtype='|U8')
        self.ux_id = np.array([], dtype='int32')
        self.mesh = np.array([], dtype='|U8')
        self.force = np.array([], dtype='int32')
        self.dmik = np.array([], dtype='|U8')
        self.perq = np.array([], dtype='|U8')

    def add(self, mach: float, sym_xz: str, sym_xy: str, ux_id: int,
            mesh: str, force: int, dmik: str, perq: str='',
            ifile: int=0, comment: str='') -> int:
        card = (mach, sym_xz, sym_xy, ux_id, mesh, force, dmik, perq, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        mach = double(card, 1, 'mach')
        sym_xz = string_or_blank(card, 2, 'sym_xz', default='')
        sym_xy = string_or_blank(card, 3, 'sym_xy', default='')
        ux_id = integer(card, 4, 'ux_id')
        mesh = string_or_blank(card, 5, 'mesh', default='')
        force = integer_or_blank(card, 6, 'force', default=0)
        dmik = string_or_blank(card, 7, 'dmik', default='')
        perq = string_or_blank(card, 8, 'perq', default='')
        self.cards.append((mach, sym_xz, sym_xy, ux_id, mesh, force, dmik, perq, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        mach = np.zeros(ncards, dtype='float64')
        sym_xz = np.empty(ncards, dtype='|U8')
        sym_xy = np.empty(ncards, dtype='|U8')
        ux_id = np.zeros(ncards, dtype='int32')
        mesh = np.empty(ncards, dtype='|U8')
        force = np.zeros(ncards, dtype='int32')
        dmik = np.empty(ncards, dtype='|U8')
        perq = np.empty(ncards, dtype='|U8')
        comment = {}
        for icard, card in enumerate(self.cards):
            (machi, sym_xzi, sym_xyi, ux_idi, meshi, forcei, dmiki, perqi, ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[icard] = commenti
            mach[icard] = machi
            sym_xz[icard] = sym_xzi
            sym_xy[icard] = sym_xyi
            ux_id[icard] = ux_idi
            mesh[icard] = meshi
            force[icard] = forcei
            dmik[icard] = dmiki
            perq[icard] = perqi
        self._save(mach, sym_xz, sym_xy, ux_id, mesh, force, dmik, perq,
                   ifile=ifile, comment=comment)
        self.cards = []

    def _save(self, mach, sym_xz, sym_xy, ux_id, mesh, force, dmik, perq,
              ifile=None, comment=None):
        ncards = len(ux_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.ux_id):
            mach = np.hstack([self.mach, mach])
            sym_xz = np.hstack([self.sym_xz, sym_xz])
            sym_xy = np.hstack([self.sym_xy, sym_xy])
            ux_id = np.hstack([self.ux_id, ux_id])
            mesh = np.hstack([self.mesh, mesh])
            force = np.hstack([self.force, force])
            dmik = np.hstack([self.dmik, dmik])
            perq = np.hstack([self.perq, perq])
        save_ifile_comment(self, ifile, comment)
        self.mach = mach
        self.sym_xz = sym_xz
        self.sym_xy = sym_xy
        self.ux_id = ux_id
        self.mesh = mesh
        self.force = force
        self.dmik = dmik
        self.perq = perq
        self.n = len(ux_id)

    def __apply_slice__(self, load: AEFORCE, i: np.ndarray) -> None:
        load.n = len(i)
        load.mach = self.mach[i]
        load.sym_xz = self.sym_xz[i]
        load.sym_xy = self.sym_xy[i]
        load.ux_id = self.ux_id[i]
        load.mesh = self.mesh[i]
        load.force = self.force[i]
        load.dmik = self.dmik[i]
        load.perq = self.perq[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        geom_check(self, missing)

    @property
    def max_id(self) -> int:
        return self.ux_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        for mach, sym_xz, sym_xy, ux_id, mesh, force, dmik, perq in zip(
                self.mach, self.sym_xz, self.sym_xy, self.ux_id,
                self.mesh, self.force, self.dmik, self.perq):
            list_fields = ['AEFORCE', mach,
                           sym_xz if sym_xz else None,
                           sym_xy if sym_xy else None,
                           ux_id,
                           mesh if mesh else None,
                           force if force else None,
                           dmik if dmik else None,
                           perq if perq else None]
            bdf_file.write(print_card(list_fields))
        return


class UXVEC(VectorizedBaseCard):
    """
    Defines the state of the aerodynamic extra points for a trim analysis.
    All undefined extra points will be set to zero.

    +-------+--------+------+--------+--------+-----+--------+-----+----------+
    |   1   |   2    |   3  |    4   |    5   |  6  |    7   |  8  |     9    |
    +=======+========+======+========+========+=====+========+=====+==========+
    | UXVEC |   ID   |      |        |        |     |        |     |          |
    +-------+--------+------+--------+--------+-----+--------+-----+----------+
    |       | LABEL1 |  UX1 | LABEL2 |   UX2  | ... |        |     |          |
    +-------+--------+------+--------+--------+-----+--------+-----+----------+
    """
    _id_name = 'uxvec_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.uxvec_id = np.array([], dtype='int32')
        self.nlabel = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U8')
        self.ux = np.array([], dtype='float64')

    def add(self, sid: int, labels: list[str], uxs: list[float],
            ifile: int=0, comment: str='') -> int:
        assert len(labels) == len(uxs)
        card = (sid, labels, uxs, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        sid = integer(card, 1, 'sid')

        i = 9
        n = 1
        labels = []
        uxs = []
        while i < len(card):
            label = string(card, i, 'label%d' % n)
            ux = double(card, i + 1, 'ux%d' % n)
            labels.append(label)
            uxs.append(ux)
            i += 2
            n += 1
        self.cards.append((sid, labels, uxs, ifile, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        uxvec_id = np.zeros(ncards, dtype='int32')
        nlabel = np.zeros(ncards, dtype='int32')
        label = []
        ux = []
        comment = {}
        for icard, card in enumerate(self.cards):
            (sid, labelsi, uxsi, ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[sid] = commenti
            uxvec_id[icard] = sid
            nlabel[icard] = len(labelsi)
            label.extend(labelsi)
            ux.extend(uxsi)
        label = np.array(label, dtype='|U8')
        ux = np.array(ux, dtype='float64')
        self._save(uxvec_id, nlabel, label, ux,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, uxvec_id, nlabel, label, ux,
              ifile=None, comment=None):
        ncards = len(uxvec_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.uxvec_id):
            uxvec_id = np.hstack([self.uxvec_id, uxvec_id])
            nlabel = np.hstack([self.nlabel, nlabel])
            label = np.hstack([self.label, label])
            ux = np.hstack([self.ux, ux])
        save_ifile_comment(self, ifile, comment)
        self.uxvec_id = uxvec_id
        self.nlabel = nlabel
        self.label = label
        self.ux = ux
        self.n = len(uxvec_id)

    @property
    def ilabel(self) -> np.ndarray:
        return make_idim(self.n, self.nlabel)

    def __apply_slice__(self, load: UXVEC, i: np.ndarray) -> None:
        load.n = len(i)
        load.uxvec_id = self.uxvec_id[i]
        ilabel = self.ilabel
        load.label = hslice_by_idim(i, ilabel, self.label)
        load.ux = hslice_by_idim(i, ilabel, self.ux)
        load.nlabel = self.nlabel[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        geom_check(self, missing)

    @property
    def max_id(self) -> int:
        return self.uxvec_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        uxvec_ids = array_str(self.uxvec_id, size=size)

        for uxvec_id, (ilabel0, ilabel1) in zip_longest(
                uxvec_ids, self.ilabel):

            labels = self.label[ilabel0:ilabel1]
            uxs = self.ux[ilabel0:ilabel1]
            list_fields = ['UXVEC', uxvec_id, None, None, None, None, None, None, None]
            nlabels = len(labels)
            assert nlabels > 0, labels
            for label, ux in zip(labels, uxs):
                list_fields += [label, ux]
            bdf_file.write(print_card(list_fields))
        return

AEPRESS = None
AEDW = None
