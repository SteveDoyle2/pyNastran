from __future__ import annotations
from itertools import count
from typing import Optional, TYPE_CHECKING

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string,
    string_or_blank, double_or_string, blank,
)
from pyNastran.bdf.cards.aero.aero import AELINK
from pyNastran.bdf.cards.aero.static_loads import TRIM
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class TRIM_ZONA(BaseCard):
    """
    Specifies constraints for aeroelastic trim variables.

    """
    type = 'TRIM_ZONA'
    _field_map = {
        1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    }

    def __init__(self, sid: int, mkaeroz: int, q: float,
                 cg: list[float], true_g: str,
                 nxyz: list[float], pqr: list[float], loadset: Optional[int],
                 labels: list[int], uxs: list[float], comment: str=''):
        """
        Creates a TRIM card for a static aero (144) analysis.

        Parameters
        ----------
        sid : int
            the trim id; referenced by the Case Control TRIM field
        q : float
            dynamic pressure
        true_g : str
            'TRUE': nxyz given in model units
            'G':    nxyz given in g's
        nxyz : list[float]
            g loading in xyz directions
        pqr : list[float]
            [roll_rate, pitch_rate, yaw_rate]
        loadset : int
            Identification number of a SET1 or SETADD bulk data card that
            specifies a set of identification numbers of TRIMFNC or
            TRIMADD bulk data card.  All values of the trim functions
            defined by the TRIMFNC or TRIMADD bulk data card are computed
            and printed out.
        labels : list[str]
            points to a TRIMVAR
            names of the fixed variables; TODO: why are these integers???
        uxs : list[float]
            values corresponding to labels
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        #: Trim set identification number. (Integer > 0)
        self.sid = sid
        self.mkaeroz = mkaeroz
        #: Dynamic pressure. (Real > 0.0)
        self.q = q

        self.cg = cg
        self.nxyz = nxyz
        self.true_g = true_g
        self.pqr = pqr
        self.loadset = loadset

        #: The label identifying aerodynamic trim variables defined on an
        #: AESTAT or AESURF entry.
        # points to a TRIMVAR
        self.labels = labels

        #: The magnitude of the aerodynamic extra point degree-of-freedom.
        #: (Real)
        self.uxs = uxs

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # TRIM TRIMID IDMK   QINF IDOBJ IDCONS RHOX RHOY RHOZ
        #      WTMASS WEIGHT IXX  IXY   IYY    IXZ  IYZ  IZZ
        #      TRNACC NX     NY   NZ    PDOT   QDOT RDOT LOADSET
        #      IDVAR1 VAL1 IDVAR2 VAL2
        sid = integer(card, 1, 'sid')
        mkaeroz = integer(card, 2, 'mkaeroz')
        qinf = double(card, 3, 'dynamic_pressure')
        # 5
        # 6
        cg = [
            double(card, 7, 'cg-x'),
            double(card, 8, 'cg-y'),
            double(card, 9, 'cg-z'),
        ]

        unused_wtmass = double(card, 9, 'wtmass')
        unused_weight = double(card, 10, 'weight')
        unused_inertia = [
            double(card, 11, 'Ixx'),
            double(card, 12, 'Ixy'),
            double(card, 13, 'Iyy'),
            double(card, 14, 'Ixz'),
            double(card, 15, 'Iyz'),
            double(card, 16, 'Izz'),
        ]
        #  TRUE/G  NX      NY      NZ      PDOT    QDOT    RDOT    LOADSET
        true_g = string(card, 17, 'TRUE/G')
        nx = double_or_string(card, 18, 'NX')
        ny = double_or_string(card, 19, 'NY')
        nz = double_or_string(card, 20, 'NZ')
        nxyz = [nx, ny, nz]

        p = double_or_string(card, 21, 'P')
        q = double_or_string(card, 22, 'Q')
        r = double_or_string(card, 23, 'R')
        pqr = [p, q, r]
        loadset = integer_or_blank(card, 24, 'loadset')

        labels = []
        uxs = []

        i = 25
        n = 1
        while i < len(card):
            label = integer(card, i, 'label%d' % n)
            ux = double_or_string(card, i + 1, 'ux%d' % n)
            if isinstance(ux, str):
                assert ux == 'FREE', 'ux=%r' % ux
            #print('  label=%s ux=%s' % (label, ux))
            labels.append(label)
            uxs.append(ux)
            i += 2
            n += 1
        assert len(card) >= 25, f'len(TRIM card) = {len(card):d}\ncard={card}'
        return TRIM_ZONA(sid, mkaeroz, qinf, cg, true_g, nxyz, pqr, loadset,
                         labels, uxs, comment=comment)

    def validate(self):
        assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

        assert isinstance(self.nxyz[0], float) or self.nxyz[0] in ['FREE', 'NONE'], 'nx=%r' % self.nxyz[0]
        assert isinstance(self.nxyz[1], float) or self.nxyz[1] in ['FREE', 'NONE'], 'ny=%r' % self.nxyz[1]
        assert isinstance(self.nxyz[2], float) or self.nxyz[2] in ['FREE', 'NONE'], 'nz=%r' % self.nxyz[2]

        assert isinstance(self.pqr[0], float) or self.pqr[0] in ['FREE', 'NONE'], 'p=%r' % self.pqr[0]
        assert isinstance(self.pqr[1], float) or self.pqr[1] in ['FREE', 'NONE'], 'q=%r' % self.pqr[1]
        assert isinstance(self.pqr[2], float) or self.pqr[2] in ['FREE', 'NONE'], 'r=%r' % self.pqr[2]

        assert self.q > 0.0, 'q=%s\n%s' % (self.q, str(self))
        if len(set(self.labels)) != len(self.labels):
            msg = 'not all labels are unique; labels=%s' % str(self.labels)
            raise RuntimeError(msg)
        if len(self.labels) != len(self.uxs):
            msg = 'nlabels=%s != nux=%s; labels=%s uxs=%s' % (
                len(self.labels), len(self.uxs), str(self.labels), str(self.uxs))
            raise RuntimeError(msg)

    def cross_reference(self, model: BDF) -> None:
        self.mkaeroz_ref = model.zona.mkaeroz[self.mkaeroz]
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model: BDF):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def convert_to_nastran(self, model: BDF):
        mkaeroz_id = self.mkaeroz
        mkaeroz = model.zona.mkaeroz[mkaeroz_id]
        #print(mkaeroz)
        mach = mkaeroz.mach
        labels = []
        uxs = []
        comment = str(self)
        for label_id, ux in zip(self.labels, self.uxs):
            if ux != 'FREE':
                trimvar = model.zona.trimvar[label_id]
                label = trimvar.label
                assert isinstance(label, str), 'label=%r' % label
                comment += str(trimvar)
                labels.append(label)
                uxs.append(ux)

        assert self.q is not None
        if self.q == 'NONE':
            self.q = 1.
        assert isinstance(self.q, float), str(self)
        trim = TRIM(self.sid, mach, self.q, labels, uxs,
                    aeqr=1.0, comment=comment)
        trim.validate()
        return trim

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        mach = 1.0
        aeqr = 0.0
        list_fields = ['TRIM', self.sid, mach, self.q]
        nlabels = len(self.labels)
        assert nlabels > 0, self.labels
        for (i, label, ux) in zip(count(), self.labels, self.uxs):
            list_fields += [label, ux]
            if i == 1:
                list_fields += [aeqr]
        if nlabels == 1:
            list_fields += [None, None, aeqr]
        return list_fields

    def repr_fields(self):
        # fixes a Nastran bug
        #aeqr = set_blank_if_default(self.aeqr, 1.0)
        aeqr = 0.

        list_fields = self.raw_fields()
        list_fields[8] = aeqr
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        return ''
        #card = self.repr_fields()
        #return self.comment + print_card_8(card)


class TRIMLNK(BaseCard):
    """
    Defines a set of coefficient and trim variable identification
    number pairs for trim variable linking.

    +=========+========+========+========+========+========+========+========+========+
    |    1    |    2   |    3   |    4   |   5   |    6    |    7   |   8    |    9   |
    +---------+--------+--------+--------+--------+--------+--------+--------+--------+
    | TRIMLNK | IDLINK |   SYM  | COEFF1 | IDVAR1 | COEFF2 | IDVAR2 | COEFF3 | IDVAR3 |
    +---------+--------+--------+--------+--------+--------+--------+--------+--------+
    |         | COEFF4 | IDVAR4 |  etc.  |        |        |        |        |        |
    +---------+--------+--------+--------+--------+--------+--------+--------+--------+

    """
    type = 'TRIMLNK'

    def __init__(self, link_id, sym, coeffs, var_ids, comment=''):
        """
        Creates a TRIMLNK card

        Parameters
        ----------
        link_id : int
            the TRIMLNK id
        sym : ???
            ???
        coeffs : ???
            ???
        var_ids : ???
            ???
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.link_id = link_id
        self.sym = sym
        self.coeffs = coeffs
        self.var_ids = var_ids
        assert sym in ['SYM', 'ASYM', 'ANTI'], sym

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIMLNK card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        link_id = integer(card, 1, 'var_id')
        sym = string_or_blank(card, 2, 'sym')

        nfields = len(card) - 3
        assert nfields % 2 == 0, card
        icoeff = 1
        coeffs = []
        var_ids = []
        for ifield in range(3, len(card), 2):
            coeff = double(card, ifield, 'coeff_%i' % icoeff)
            var_id = integer(card, ifield + 1, 'var_%i' % icoeff)
            coeffs.append(coeff)
            var_ids.append(var_id)
            icoeff += 1
        assert len(card) >= 5, f'len(TRIMLNK card) = {len(card):d}\ncard={card}'
        return TRIMLNK(link_id, sym, coeffs, var_ids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def convert_to_nastran(self, model):
        label = 'LNK_%s' % self.link_id
        trimvars = model.zona.trimvar

        comment = str(self)
        independent_labels = []
        for var_id in self.var_ids:
            trimvar = trimvars[var_id]
            label = trimvar.label
            comment += str(trimvar)
            independent_labels.append(label)

        Cis = self.coeffs
        aelink = AELINK(self.link_id, label, independent_labels, Cis, comment=comment)
        aelink.validate()
        return aelink

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['TRIMLNK', self.link_id, self.sym]
        for coeff, var in zip(self.coeffs, self.var_ids):
            list_fields.append(coeff)
            list_fields.append(var)
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class TRIMVAR(BaseCard):
    """
    Specifies a trim variable for static aeroelastic trim variables.

    """
    type = 'TRIMVAR'

    def __init__(self, var_id: int, label: str, lower: float, upper: float,
                 trimlnk: int, dmi: None, sym: int, initial: None,
                 dcd: float, dcy: float, dcl: float,
                 dcr: float, dcm: float, dcn: float, comment: str=''):
        """
        Creates a TRIMVAR card for a static aero (144) analysis.

        Parameters
        ----------
        var_id : int
            the trim id; referenced by the Case Control TRIM field
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.var_id = var_id
        self.label = label
        self.lower = lower
        self.upper = upper
        self.trimlnk = trimlnk
        self.dmi = dmi
        self.sym = sym
        self.initial = initial
        self.dcd = dcd
        self.dcy = dcy
        self.dcl = dcl
        self.dcr = dcr
        self.dcm = dcm
        self.dcn = dcn

    @classmethod
    def add_card(cls, card: BDFCard, comment:  str=''):
        """
        Adds a TRIMVAR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        var_id = integer(card, 1, 'var_id')
        label = string(card, 2, 'label')
        lower = double_or_blank(card, 3, 'lower')
        upper = double_or_blank(card, 4, 'upper')
        trimlnk = integer_or_blank(card, 5, 'TRIMLNK')
        dmi = blank(card, 6, 'DMI')
        sym = string_or_blank(card, 7, 'sym')
        initial = blank(card, 8, 'initial')
        dcd = double_or_blank(card, 9, 'DCD')
        dcy = double_or_blank(card, 10, 'DCY')
        dcl = double_or_blank(card, 11, 'DCL')
        dcr = double_or_blank(card, 12, 'DCR')
        dcm = double_or_blank(card, 13, 'DCM')
        dcn = double_or_blank(card, 14, 'DCN')
        return TRIMVAR(var_id, label, lower, upper, trimlnk, dmi, sym,
                       initial, dcd, dcy, dcl, dcr, dcm,
                       dcn, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model: BDF):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def convert_to_nastran(self, model: BDF):
        raise NotImplementedError()
        #mkaeroz_id = self.mkaeroz
        #mkaeroz = model.zona.mkaeroz[mkaeroz_id]
        #mach = mkaeroz.mach
        #labels = []
        #uxs = []
        #for label_id, ux in zip(self.labels, self.uxs):
            #if ux != 'FREE':
                #label = model.zona.trimvar[label_id]
                #labels.append(label)
                #uxs.append(ux)
        #trim = TRIM(self.sid, mach, self.q, labels, uxs,
                    #aeqr=1.0, comment=str(self))
        #trim.validate()
        #return trim

    def raw_fields(self) -> list:
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['TRIMVAR', self.var_id, self.label, self.lower, self.upper,
                       self.trimlnk, self.dmi, self.sym, self.initial,
                       self.dcd, self.dcy, self.dcl, self.dcr, self.dcm, self.dcn]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)
