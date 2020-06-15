"""
All superelements are defined in this file.  This includes:
 * CAABSF
 * CHACAB
 * CHACBR

"""
from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.bdf.cards.base_card import Element, Property
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (
    double,
    integer, integer_or_blank, # integer_or_string,
    string,
    string_or_blank, double_or_blank, # integer_string_or_blank,
    #exact_string_or_blank,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class CHACAB(Element):
    """
    Acoustic Absorber Element Connection
    Defines the acoustic absorber element in coupled fluid-structural analysis.

    | CHACAB | EID | PID | G1  | G2  | G3  | G4  | G5 | G6 |
    |        | G7  | G8  | G9  | G10 | G11 | G12 |    |    |
    |        |     |     | G17 | G18 | G19 | G20 |    |    |
    """
    type = 'CHACAB'

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        nodes = [1]
        return CHACAB(eid, pid, nodes, comment='')

    def __init__(self, eid, pid, nodes, comment=''):
        """
        EID Element identification number. (0 < Integer < 100,000,000)
        PID Property identification number of a PACABS entry. (Integer > 0)
        Gi Grid point identification numbers of connection points. (Integer > 0 or blank)
        """
        Element.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.nodes = self.prepare_node_ids(nodes, allow_empty_nodes=True)
        self.nodes_ref = None
        self.pid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CHACAB card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'), integer(card, 8, 'nid6'),
            integer(card, 9, 'nid7'), integer(card, 10, 'nid8'),
            integer_or_blank(card, 11, 'nid9'),
            integer_or_blank(card, 12, 'nid10'),
            integer_or_blank(card, 13, 'nid11'),
            integer_or_blank(card, 14, 'nid12'),
            None,
            None,
            None,
            None,
            #integer_or_blank(card, 15, 'nid13'),
            #integer_or_blank(card, 16, 'nid14'),
            #integer_or_blank(card, 17, 'nid15'),
            #integer_or_blank(card, 18, 'nid16'),
            integer_or_blank(card, 19, 'nid17'),
            integer_or_blank(card, 20, 'nid18'),
            integer_or_blank(card, 21, 'nid19'),
            integer_or_blank(card, 22, 'nid20'),
        ]
        assert len(card) <= 23, 'len(CHACAB card) = %i\ncard=%s' % (len(card), card)
        return CHACAB(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        #msg = f', which is required by CHACAB eid={self.eid}'
        #model.Property(self.pid, msg=msg)
        pass

    def safe_cross_reference(self, model: BDF, xref_error) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        self.pid_ref = None

    def Pid(self) -> int:
        if self.pid_ref is not None:
            return self.pid_ref.pid
        return self.pid

    @property
    def node_ids(self):
        nids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)
        return nids

    def raw_fields(self):
        nodes1 = self.nodes[:12]
        nodes2 = self.nodes[16:]
        assert len(nodes1) == 12, len(nodes1)
        assert len(nodes2) == 4, len(nodes2)
        return ['CHACAB', self.eid, self.pid, ] + nodes1 + [None, None, None, None] + nodes2

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)


class CHACBR(Element):
    """

    | CHACBR | EID | PID | G1  | G2  | G3  | G4  | G5 | G6 |
    |        | G7  | G8  | G9  | G10 | G11 | G12 |    |    |
    |        |     |     | G17 | G18 | G19 | G20 |    |    |
    """
    type = 'CHACBR'
    def __init__(self, eid, pid, nodes, comment=''):
        Element.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.nodes = self.prepare_node_ids(nodes, allow_empty_nodes=True)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CHACBR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'), integer(card, 8, 'nid6'),
            integer(card, 9, 'nid7'), integer(card, 10, 'nid8'),
            integer_or_blank(card, 11, 'nid9'),
            integer_or_blank(card, 12, 'nid10'),
            integer_or_blank(card, 13, 'nid11'),
            integer_or_blank(card, 14, 'nid12'),
            None,
            None,
            None,
            None,
            #integer_or_blank(card, 15, 'nid13'),
            #integer_or_blank(card, 16, 'nid14'),
            #integer_or_blank(card, 17, 'nid15'),
            #integer_or_blank(card, 18, 'nid16'),
            integer_or_blank(card, 19, 'nid17'),
            integer_or_blank(card, 20, 'nid18'),
            integer_or_blank(card, 21, 'nid19'),
            integer_or_blank(card, 22, 'nid20'),
        ]
        assert len(card) <= 23, 'len(CHACAB card) = %i\ncard=%s' % (len(card), card)
        return CHACAB(eid, pid, nids, comment=comment)

    def raw_fields(self):
        nodes1 = self.nodes[:12]
        nodes2 = self.nodes[16:]
        assert len(nodes1) == 12, len(nodes1)
        assert len(nodes2) == 4, len(nodes2)
        return ['CHACBR', self.eid, self.pid, ] + nodes1 + [None, None, None, None] + nodes2

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)


class CAABSF(Element):
    """
    Frequency-Dependent Acoustic Absorber Element

    Defines a frequency-dependent acoustic absorber element in coupled fluid-structural
    analysis.

    +--------+-------+-------+----+----+----+----+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |
    +========+=======+=======+=====+===+====+====+
    | CAABSF |  EID  |  PID  | N1 | N2 | N3 | N4 |
    +--------+-------+-------+----+----+----+----+
    """
    type = 'CAABSF'

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        pid = 1
        nodes = [1, 2, 3, 4]
        return CAABSF(eid, pid, nodes, comment='')

    def __init__(self, eid, pid, nodes, comment=''):
        """
        n1 - point source
        n1-n2 - line source
        n1-n3 - tri face
        n1-n4 - quad face
        """
        Element.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.pid = pid
        self.nodes_ref = None
        self.pid_ref = None
        self.nodes = self.prepare_node_ids(nodes, allow_empty_nodes=True)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CHACAB card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'),
            integer_or_blank(card, 4, 'nid2', 0),
            integer_or_blank(card, 5, 'nid3', 0),
            integer_or_blank(card, 6, 'nid4', 0),
        ]
        assert len(card) <= 7, 'len(CAABSF card) = %i\ncard=%s' % (len(card), card)
        return CAABSF(eid, pid, nids, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
    def safe_cross_reference(self, model: BDF, xref_error) -> None:
        pass
    def uncross_reference(self) -> None:
        pass

    @property
    def node_ids(self):
        nids = self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True)
        return nids

    def raw_fields(self):
        nodes1 = self.nodes
        assert len(nodes1) == 4, len(nodes1)
        return ['CAABSF', self.eid, self.pid, ] + self.nodes

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)

class PACBAR(Property):
    """
    +--------+-----+-------+--------+--------+--------+
    | PACBAR | PID | MBACK | MSEPTM | FRESON | KRESON |
    +--------+-----+-------+--------+--------+--------+
    """
    type = 'PACBAR'

    def __init__(self, pid, mback, mseptm, freson, kreson, comment=''):
        """
        Creates a PACBAR card

        Parameters
        ----------
        pid : int
            Property identification number. (Integer > 0)
        mback : float
            Mass per unit area of the backing material
        mseptm : float
            Mass per unit area of the septum material
        freson : float; default=None
            Resonant frequency of the sandwich construction in hertz.
        kreson : float; default=None
            Resonant stiffness of the sandwich construction.

        """
        super().__init__()
        if comment:
            self.comment = comment
        self.pid = pid
        self.mback = mback
        self.mseptm = mseptm
        self.freson = freson
        self.kreson = kreson

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PACBAR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mback = double(card, 2, 'mback')
        mseptm = double(card, 3, 'mseptm')
        freson = double_or_blank(card, 4, 'freson')
        kreson = double_or_blank(card, 5, 'kreson')
        assert len(card) <= 6, 'len(PACBAR card) = %i\ncard=%s' % (len(card), card)
        return PACBAR(pid, mback, mseptm, freson, kreson, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
    def safe_cross_reference(self, model: BDF, xref_error) -> None:
        pass
    def uncross_reference(self) -> None:
        pass

    def raw_fields(self):
        list_fields = [
            'PACBAR', self.pid, self.mback, self.mseptm, self.freson, self.kreson]
        return list_fields

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)

class PAABSF(Property):
    """
    +--------+-----+--------+--------+---+---+---+---+------+
    | PAABSF | PID | TZREID | TZIMID | S | A | B | K | RHOC |
    +--------+-----+--------+--------+---+---+---+---+------+
    """
    type = 'PAABSF'

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        return PAABSF(pid, tzreid=None, tzimid=None,
                      s=1.0, a=1.0, b=0.0, k=0.0, rhoc=1.0,
                      comment='')

    def __init__(self, pid, tzreid=None, tzimid=None,
                 s=1.0, a=1.0, b=0.0, k=0.0, rhoc=1.0,
                 comment=''):
        """
        Creates a PAABSF card

        Parameters
        ----------
        pid : int
            Property identification number of a CAABSF. (Integer > 0)
        TZREID : int; default=None
            TABLEDi id that defines the resistance as a function of
            frequency. The real part of the impedence.
        TZIMID : int; default=None
            TABLEDi id that defines the reactance as a function of
            frequency. The imaginary part of the impedance.
        S : float; default=1.0
            Impedance scale factor.
        A : float; default=1.0
            Area factor when 1 or 2 grid points are specified on the
            CAABSF entry.
        B : float; default=0.0
            Equivalent structural damping coefficient.
        K : float; default=0.0
            Equivalent structural stiffness coefficient.
        RHOC : float; default=1.0
            Constant used in data recovery for calculating an absorption
            coefficient. RHO is the media density, and C is the speed of
            sound in the media.
        """
        super().__init__()
        if comment:
            self.comment = comment
        self.pid = pid
        self.tzreid = tzreid
        self.tzimid = tzimid
        self.s = s
        self.a = a
        self.b = b
        self.k = k
        self.rhoc = rhoc

    def validate(self):
        assert not all([self.tzreid is None, self.tzimid is None,
                        self.b is None, self.k is None]), str(self)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PACABS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        tzreid = integer_or_blank(card, 2, 'tzreid')
        tzimid = integer_or_blank(card, 3, 'tzimid')
        s = double_or_blank(card, 4, 's', 1.0)
        a = double_or_blank(card, 5, 'a', 1.0)
        b = double_or_blank(card, 6, 'b', 0.0)
        k = double_or_blank(card, 7, 'k', 0.0)
        rhoc = double_or_blank(card, 8, 'rhoc', 1.0)
        assert len(card) <= 9, 'len(PAABSF card) = %i\ncard=%s' % (len(card), card)
        return PAABSF(pid, tzreid, tzimid, s, a, b, k, rhoc, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
    def safe_cross_reference(self, model: BDF, xref_error) -> None:
        pass
    def uncross_reference(self) -> None:
        pass

    def raw_fields(self):
        list_fields = [
            'PAABSF', self.pid, self.tzreid, self.tzimid,
            self.s, self.a, self.b, self.k, self.rhoc]
        return list_fields

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)

class PACABS(Element):
    """
    +--------+-----+-------+------+------+------+-------+-----+
    | PACABS | PID | SYNTH | TID1 | TID2 | TID3 | CUTFR |  B  |
    +--------+-----+-------+------+------+------+-------+-----+
    |        |  K  |   M   |      |      |      |       |     |
    +--------+-----+-------+------+------+------+-------+-----+
    """
    type = 'PACABS'

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        cutfr = 1.0
        b = 1.
        k = 1.
        m = 1.
        return PACABS(pid, cutfr, b, k, m)

    def __init__(self, pid, cutfr, b, k, m,
                 synth=True, tid_resistance=None, tid_reactance=None, tid_weight=None,
                 comment=''):
        """
        Creates a PACABS card

        Parameters
        ----------
        pid : int
            Property identification number.
        synth : bool; default=True
            Request the calculation of B, K, and M from the tables TIDi below
        tid_resistance : int; default=None
            Identification of the TABLEDi entry that defines the resistance.
        tid_reactance : int; default=None
            Identification of the TABLEDi entry that defines the reactance.
        tid_weight : int; default=None
            Identification of the TABLEDi entry that defines the weighting function.
        cutfr : float
            Cutoff frequency for tables referenced above. (Real > 0.0)
        B, K, M : float
            Equivalent damping, stiffness and mass values per unit area. (Real > 0.0)

        ..note:: tables are defined as a function of frequency in cycles/time
        """
        super().__init__()
        if comment:
            self.comment = comment
        self.pid = pid
        self.synth = synth
        assert isinstance(synth, bool), synth
        self.tid_resistance = tid_resistance
        self.tid_reactance = tid_reactance
        self.tid_weight = tid_weight
        self.cutfr = cutfr
        self.b = b
        self.k = k
        self.m = m

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PACABS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        synth = string_or_blank(card, 2, 'synth', 'YES')
        tid_resistance = integer(card, 3, 'tid_resistance')
        tid_reactance = integer(card, 4, 'tid_reactance')
        tid_weight = integer(card, 5, 'tid_weight')
        #
        cutfr = double(card, 7, 'cutfr')
        b = double_or_blank(card, 8, 'b')
        k = double_or_blank(card, 9, 'k')
        m = double_or_blank(card, 10, 'm')
        assert len(card) <= 8, 'len(PACABS card) = %i\ncard=%s' % (len(card), card)

        assert synth in ['YES', 'NO'], synth
        is_synth = synth == 'YES'
        return PACABS(pid, cutfr, b, k, m, synth=is_synth,
                      tid_resistance=tid_resistance, tid_reactance=tid_reactance, tid_weight=tid_weight,
                      comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
    def safe_cross_reference(self, model: BDF, xref_error) -> None:
        pass
    def uncross_reference(self) -> None:
        pass

    def raw_fields(self):
        synth = 'YES' if self.synth else 'NO'
        list_fields = [
            'PACABS', self.pid, synth, self.tid_resistance, self.tid_reactance, self.tid_weight,
            None, self.cutfr, self.b, self.k, self.m]
        return list_fields

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)

def is_msc(nastran_version: str):
    return nastran_version == 'msc'

class ACMODL(Element):
    """
    NX 2019.2
    +--------+-------+---------+----------+------+--------+--------+---------+----------+
    | ACMODL |       |  INFOR  |   FSET   | SSET | NORMAL |        | OVLPANG | SRCHUNIT |
    +--------+-------+---------+----------+------+--------+--------+---------+----------+
    |        | INTOL | AREAOP  |   CTYPE  |      |        |        |         |          |
    +--------+-------+---------+----------+------+--------+--------+---------+----------+


    MSC 2016.1
    +--------+-------+---------+----------+------+--------+--------+---------+----------+
    | ACMODL | INTER |  INFOR  |   FSET   | SSET | NORMAL | METHOD | SKNEPSG | DSKNEPS  |
    +--------+-------+---------+----------+------+--------+--------+---------+----------+
    |        | INTOL | ALLSSET | SRCHUNIT |      |        |        |         |          |
    +--------+-------+---------+----------+------+--------+--------+---------+----------+

    """
    type = 'ACMODL'

    @classmethod
    def _init_from_empty(cls):
        infor = 1
        fset = 1.0
        sset = 1.
        return ACMODL(infor, fset, sset)

    def __init__(self, infor: str, fset: int, sset: int,
                 normal=0.5, olvpang=60., search_unit='REL',
                 intol=0.2, area_op=0, ctype='STRONG',
                 method='BW',
                 sk_neps=0.5, dsk_neps=0.75, all_set='NO',
                 inter='DIFF',
                 nastran_version='nx', comment=''):
        """
        Creates a ACMODL card

        Parameters
        ----------
        INFOR : str
            Defines the meaning of the SID entered on the FSET and SSET fields.
            {ELEMENTS, PID, SET3}
        FSET : int; default=None
            Selects the ID of a SET1 or SET3 entry to define the fluid
            elements for the interface.  If the ID is entered, the
            corresponding fluid elements are considered.  If a negative
            sign is included in front of the ID, the corresponding fluid
            elements are excluded.  If blank, all fluid elements are considered.
        SSET : int; default=None
            Selects the ID of a SET1 or SET3 entry to define the structural
            elements for the interface.  If the ID is entered, the
            corresponding structural elements are considered.
            If a negative sign is included in front of the ID, the
            corresponding structural elements are excluded.
            If blank, all structural elements are considered.
        NORMAL : float; default=0.5
            Outward normal search distance to detect fluid-structure interface.
            SRCHUNIT=REL:
               NORMAL is a ratio of the height of the fluid box in the outward
               normal direction to the fluid surface to the maximum edge
               length of the fluid free face.
            SRCHUNIT=ABS:
               NORMAL is the outward search distance in the model/absolute units.
        OVLPANG :float; default=60.0
            Angular tolerance in degrees used to decide whether a fluid free face
            and a structural face can be considered as overlapping. If the angle
            between the normals of the fluid and structural faces exceeds this
           value, they cannot be coupled. (0.0 < Real ≤ 90.0; Default = 60.0)
        SRCHUNIT : str; default=REL
            Search units.
            ABS: then the model units are absolute.
            REL : then the relative model units are based on element size.
        INTOL : float; default=0.2
            Inward normal search distance to detect fluid-structure interface.
            SRCHUNIT = REL:
                INTOL is a ratio of the height of the fluid box in the
                inward normal direction to the fluid surface to the
                maximum edge length of the fluid free face.
            SRCHUNIT = ABS:
               then INTOL is the inward search distance in the
               model/absolute units.
        AREAOP : int; default=0
            Alternative fluid-structure coupling method selection.
            0 = The recommended method is used.
            1 = The RBE3 method is used.
        CTYPE : str; default='STRONG'
            Fluid-structure coupling type (only supported by new acoustics
            method introduced in NX Nastran 11).
            CTYPE = STRONG:
                two-way coupling is turned on.
            CTYPE = WEAK:
                one-way coupling is turned on. Here, the effect of
                the fluid on the structure is assumed to be negligible.
        DSKNEPS : float; default = 0.75
            Secondary fluid skin growth tolerance
        INTOL : float; default = 0.5
            Tolerance of inward normal.
        ALLSSET : str; default=NO
           NO then SSET structure is searched and coupled if found.
           YES then all the structure given by SSET is coupled.
        SRCHUNIT Search units ; str; default = REL
            ABS for absolute model units
            REL for relative model units based on element size

        """
        super().__init__()
        if comment:
            self.comment = comment
        self.nastran_version = nastran_version

        # NX
        self.infor = infor
        self.fset = fset
        self.sset = sset
        self.normal = normal
        self.olvpang = olvpang
        self.search_unit = search_unit
        self.intol = intol
        self.area_op = area_op
        self.ctype = ctype

        # MSC only
        self.method = method
        self.sk_neps = sk_neps
        self.dsk_neps = dsk_neps
        self.all_set = all_set
        self.inter = inter

    @classmethod
    def add_card(cls, card, nastran_version, comment=''):
        """
        Adds a ACMODL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        if is_msc(nastran_version):
            return cls.add_card_msc(card, comment=comment)
        return cls.add_card_nx(card, comment=comment)

    @classmethod
    def add_card_nx(cls, card, comment=''):
        #
        infor = string(card, 2, 'infor')
        fset = integer(card, 3, 'fset')
        sset = integer(card, 4, 'sset')
        normal = double_or_blank(card, 5, 'normal', 0.5)
        #
        olvpang = double_or_blank(card, 7, 'olvpang', 60.0)
        search_unit = string_or_blank(card, 8, 'search_unit', 'REL')
        intol = double_or_blank(card, 9, 'intol', 0.2)
        area_op = integer_or_blank(card, 10, 'area_op', 0)
        ctype = string_or_blank(card, 11, 'ctype', 'STRONG')
        assert len(card) <= 8, 'len(ACMODL card) = %i\ncard=%s' % (len(card), card)

        return ACMODL(infor, fset, sset,
                      normal=normal, olvpang=olvpang,
                      search_unit=search_unit, intol=intol,
                      area_op=area_op, ctype=ctype,
                      nastran_version='nx', comment=comment)

    @classmethod
    def add_card_msc(cls, card, comment=''):
        #print('MSC...ACMODL')
        inter = string_or_blank(card, 2, 'infor', 'DIFF')
        infor = string_or_blank(card, 2, 'infor', 'NONE')
        fset = integer_or_blank(card, 3, 'fset')
        sset = integer_or_blank(card, 4, 'sset')
        normal_default = 0.001 if inter == 'INDENT' else 1.0
        normal = double_or_blank(card, 5, 'normal', normal_default)
        method = string_or_blank(card, 6, 'method', 'BW') # BW/CP
        sk_neps = double_or_blank(card, 7, 'sk_neps', 0.5)
        dsk_neps = double_or_blank(card, 8, 'dsk_neps', 0.75)
        intol = double_or_blank(card, 9, 'intol', 0.2)
        all_set = string_or_blank(card, 10, 'all_set', 'NO') #
        search_unit = string_or_blank(card, 11, 'search_unit', 'REL')
        assert len(card) <= 12, 'len(ACMODL card) = %i\ncard=%s' % (len(card), card)
        return ACMODL(infor, fset, sset,
                      inter=inter, normal=normal, method=method,
                      sk_neps=sk_neps, dsk_neps=dsk_neps, all_set=all_set,
                      search_unit=search_unit, intol=intol,
                      nastran_version='msc',
                      comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
    def safe_cross_reference(self, model: BDF, xref_error) -> None:
        pass
    def uncross_reference(self) -> None:
        pass

    def raw_fields(self):
        if is_msc(self.nastran_version):
            list_fields = ['ACMODL', self.inter, self.infor,
                           self.fset, self.sset, self.normal, self.method,
                           self.sk_neps, self.dsk_neps, self.intol,
                           self.all_set, self.search_unit, ]
        else:
            list_fields = ['ACMODL', None, self.infor, self.fset, self.sset, self.normal, None, self.olvpang, self.search_unit,
                           self.intol, self.area_op, None, None, self.ctype]
            #list_fields = [
                #'ACMODL', self.pid, synth, self.tid_resistance, self.tid_reactance, self.tid_weight,
                #None, self.cutfr, self.b, self.k, self.m]
        return list_fields

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)
