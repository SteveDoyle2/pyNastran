from __future__ import annotations
from typing import Optional, TYPE_CHECKING

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, string,
    string_or_blank, string_multifield_or_blank,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class MKAEROZ(BaseCard):
    type = 'MKAEROZ'

    def __init__(self, sid: int, mach: float, flt_id: int, filename: str,
                 print_flag: int, freqs: list[float],
                 method: int=0, save: Optional[str]=None, comment: str=''):
        """
        Parameters
        ==========
        sid : int
            the MKAEROZ id
        mach : float
            the mach number for the TRIM solution
        save : str
            save the AIC data to the filename
            SAVE    save the AICs
            ACQUIRE load an AIC database
            ADD     append the new acids to the existing AIC database
            RESTART continue an analysis
        filename : str
            the length of the file must be at most 56 characters
        print_flag : int
            ???
        freqs : list[float]
            ???
        method : int
            ???
        save : ???
            ???
        comment : str; default=''
             a comment for the card
        """
        BaseCard.__init__(self)

        if comment:
            self.comment = comment
        self.sid = sid
        self.mach = mach
        self.method = method
        self.flt_id = flt_id
        self.save = save
        self.freqs = freqs
        self.filename = filename
        self.print_flag = print_flag

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        sid = integer(card, 1, 'IDMK')
        mach = double(card, 2, 'MACH')
        method = integer(card, 3, 'METHOD')
        flt_id = integer_or_blank(card, 4, 'IDFLT')
        save = string_or_blank(card, 5, 'SAVE')
        #print(f'card = {card}')
        filename = string_multifield_or_blank(card, (6, 7), 'filename', default='')
        print_flag = integer_or_blank(card, 8, 'PRINT_FLAG', 0)
        freqs = []
        ifreq = 1
        for ifield in range(9, len(card)):
            freq = double(card, ifield, 'FREQ%d' % ifreq)
            freqs.append(freq)
            ifreq += 1
        return MKAEROZ(sid, mach, flt_id, filename, print_flag, freqs,
                       method=method, save=save, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        return

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        filename_a = self.filename[:8]
        filename_b = self.filename[8:]
        list_fields = ['MKAEROZ', self.sid, self.mach, self.method, self.flt_id,
                       self.save, filename_a, filename_b, self.print_flag] + self.freqs
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class FLUTTER_ZONA(BaseCard):
    """
    Defines data needed to perform flutter, ASE, or a transient response analysis.

    +---------+-----+--------+------+------+-------+-------+-------------+------+
    |    1    |  2  |   3    |  4   |  5   |   6   |   7   |      8      |  9   |
    +=========+=====+========+======+======+=======+=======+=============+======+
    | FLUTTER | SID | METHOD | DENS | MACH | RFREQ | IMETH | NVALUE/OMAX | EPS  |
    +---------+-----+--------+------+------+-------+-------+-------------+------+
    | FLUTTER | 19  |   K    | 119  | 219  | 319   |   S   |      5      | 1.-4 |
    +---------+-----+--------+------+------+-------+-------+-------------+------+

    +---------+-------+-------+-------+-------+--------+-------+---------+--------+
    |    1    |   2   |   3   |   4   |   5   |    6   |   7   |    8    |    9   |
    +=========+=======+=======+=======+=======+========+=======+=========+========+
    | FLUTTER | SETID |  SYM  |  FIX  | NMODE | TABDMP | MLIST | CONMLST | NKSTEP |
    +---------+-------+-------+-------+-------+--------+-------+---------+--------+
    | FLUTTER |  100  | SYMM3 |   1   |   0   |  30    |  100  |    0    |   50   |
    +---------+-------+-------+-------+-------+--------+-------+---------+--------+

    """
    type = 'FLUTTER_ZONA'

    def __init__(self, sid: int, sym: str, fix: int, mlist: int, conmlst: int,
                 nmode: int=0, tabdmp: int=0, nkstep: int=25, comment: str=''):
        """
        Creates a FLUTTER card, which is required for a flutter (SOL 145)
        analysis.

        Parameters
        ----------
        sid : int
            Unique set identification number. (Integer > 0)
        sym : str
           Character string up to 8 characters with no embedded blanks.
           The first 4 characters can be either 'SYMM' (or 'SYM'), 'ANTI',
           or 'ASYM' that defines the boundary condition of the structural
           finite element model as well as the unsteady aerodynamics, where:
            - SYMM Symmetric boundary condition.
            - ANTI Antisymmetric boundary condition.
            - ASYM Asymmetric boundary condition.
          The last 4 characters are used to specify the interpolation scheme
          for the generalized aerodynamic matrices. They can be either:
           - blank for a cubic spline
           - L for a linear interpolation.
             (such as SYM = 'SYMML', 'ANTIL', or 'ASYML')
           - P for a second-order-polynomial interpolation.
             (such as SYM = 'SYMMP', 'ANTIP', or 'ASYMP')
           - integer for a hybrid cubic spline and linear interpolation scheme.
             (such as SYM = 'SYMM1', 'SYMM2', 'ANTI3', etc.)
           - (Default = SYMML)
        fix : int
           Identification number of a FIXHATM, FIXMATM, FIXMACH, or FIXMDEN
           bulk data card. (Integer > 0)
        nmode : int
            Number of structural modes used in the flutter analysis. (Integer >= 0)
            If NMODE = 0, all structural modes are used for flutter analysis
        tabdmp : int
            Identification number of a TABDMP1 bulk data card specifying modal damping as
           a function of natural frequency. (Integer ≥ 0)
        mlist : int
            Identification number of a SET1 or SETADD bulk data card
            specifying a list of normal modes to be omitted from the
            flutter analysis. (Integer >= 0)
        conmlst : int
            Identification number of a CONMLST bulk data card specifying
            a list of identification numbers of the CONM1 bulk data cards
            for mass perturbation. (Integer >= 0) (See Remark 8)
        nkstep : int; default=25
            Number of reduced frequency steps for the reduced-frequency-sweep
            technique of the g-Method flutter solution. (Integer >= 0)

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        sym = sym.upper()
        if sym == 'SYMM':
            sym = 'SYM'
        self.sid = sid
        self.sym = sym
        self.fix = fix
        self.nmode = nmode
        self.tabdmp = tabdmp
        self.mlist = mlist
        self.conmlst = conmlst
        self.nkstep = nkstep
        assert sym in {'SYM', 'ANTI', 'ASYM',
                       'SYMML', 'ANTIL', 'ASYMP',
                       'SYMMP', 'ANTIP',}, f'FLUTTER sid={self.sid} sym={sym!r}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a FLUTTER card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        sym = string(card, 2, 'sym')
        fix = integer(card, 3, 'fix')
        nmode = integer_or_blank(card, 4, 'nmode', default=0)
        tabdmp = integer_or_blank(card, 5, 'tabdmp', default=0)
        mlist = integer_or_blank(card, 6, 'mlist')
        conmlst = integer_or_blank(card, 7, 'conmlst')
        nkstep = integer_or_blank(card, 8, 'nkstep', default=25)
        assert len(card) <= 9, f'len(FLUTTER card) = {len(card):d}\ncard={card}'
        flutter = FLUTTER_ZONA(sid, sym, fix, mlist, conmlst,
                               nmode=nmode, tabdmp=tabdmp, nkstep=nkstep, comment=comment)
        return flutter

    def cross_reference(self, model: BDF) -> None:
        return
        #msg = ', which is required by SPLINE1 eid=%s' % self.eid
        #self.setg_ref = model.Set(self.setg, msg=msg)
        #self.setg_ref.cross_reference_set(model, 'Node', msg=msg)

        #self.panlst_ref = model.zona.panlsts[self.panlst]
        #self.panlst_ref.cross_reference(model)
        #self.aero_element_ids = self.panlst_ref.aero_element_ids

    def safe_cross_reference(self, model: BDF, xref_errors=None):
        return
        #msg = ', which is required by SPLINE1 eid=%s' % self.eid
        #try:
            #self.setg_ref = model.Set(self.setg, msg=msg)
            #self.setg_ref.safe_cross_reference(model, 'Node', msg=msg)
        #except KeyError:
            #model.log.warning('failed to find SETx set_id=%s%s; allowed_sets=%s' % (
                #self.setg, msg, np.unique(list(model.sets)))

        #try:
            #self.panlst_ref = model.zona.panlsts[self.panlst]
            #self.panlst_ref.safe_cross_reference(model, xref_errors)
            #self.aero_element_ids = self.panlst_ref.aero_element_ids
        #except KeyError:
            #pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        return
        #self.panlst_ref = None
        #self.setg_ref = None

    def convert_to_nastran(self, model: BDF):
        raise NotImplementedError()

    def raw_fields(self):
        list_fields = ['FLUTTER', self.sid, self.sym, self.fix, self.nmode,
                       self.tabdmp, self.mlist, self.conmlst, self.nkstep]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


