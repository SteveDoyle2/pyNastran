"""
All superelements are defined in this file.  This includes:
 * CAABSF
 * CHACAB
 * CHACBR
 * PACBAR
 * PAABSF
 * ACMODL
 * PMIC
 * ACPLNW

"""
from __future__ import annotations
from typing import Union, TYPE_CHECKING

from pyNastran.bdf.cards.base_card import BaseCard, Element, Property
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (
    double,
    integer, integer_or_blank, # integer_or_string,
    string,
    string_or_blank, double_or_blank, # integer_string_or_blank,
    integer_double_or_blank,
    #exact_string_or_blank,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.bdf.bdf import BDF


class CHACAB(Element):
    """
    Acoustic Absorber Element Connection
    Defines the acoustic absorber element in coupled fluid-structural analysis.

    +--------+-----+-----+-----+-----+-----+-----+----+----+
    |   1    |  2  |  3  |  4  |  5  |  6  |  7  |  8 |  9 |
    +========+=====+=====+=====+=====+=====+=====+====+====+
    | CHACAB | EID | PID | G1  | G2  | G3  | G4  | G5 | G6 |
    +--------+-----+-----+-----+-----+-----+-----+----+----+
    |        | G7  | G8  | G9  | G10 | G11 | G12 |    |    |
    +--------+-----+-----+-----+-----+-----+-----+----+----+
    |        |     |     | G17 | G18 | G19 | G20 |    |    |
    +--------+-----+-----+-----+-----+-----+-----+----+----+
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
        assert len(card) <= 23, f'len(CHACAB card) = {len(card):d}\ncard={card}'
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
        assert len(card) <= 23, f'len(CHACAB card) = {len(card):d}\ncard={card}'
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
    +========+=======+=======+====+====+====+====+
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
        assert len(card) <= 7, f'len(CAABSF card) = {len(card):d}\ncard={card}'
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
        assert len(card) <= 6, f'len(PACBAR card) = {len(card):d}\ncard={card}'
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
    |    1   |  2  |    3   |    4   | 5 | 6 | 7 | 8 |   9  |
    +========+=====+========+========+===+===+===+===+=======
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

    def __init__(self, pid: int,
                 tzreid: int=0, tzimid: int=0,
                 s: float=1.0, a: float=1.0, b: float=0.0,
                 k: float=0.0, rhoc: float=1.0,
                 comment: str=''):
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
        Adds a PAABSF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        tzreid = integer_or_blank(card, 2, 'tzreid', default=0)
        tzimid = integer_or_blank(card, 3, 'tzimid', default=0)
        s = double_or_blank(card, 4, 's', default=1.0)
        a = double_or_blank(card, 5, 'a', default=1.0)
        b = double_or_blank(card, 6, 'b', default=0.0)
        k = double_or_blank(card, 7, 'k', default=0.0)
        rhoc = double_or_blank(card, 8, 'rhoc', default=1.0)
        assert len(card) <= 9, f'len(PAABSF card) = {len(card):d}\ncard={card}'
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
        assert len(card) <= 8, f'len(PACABS card) = {len(card):d}\ncard={card}'

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


class ACPLNW(BaseCard):
    type = 'ACPLNW'

    @classmethod
    def _init_from_empty(cls):
        infor = 1
        fset = 1.0
        sset = 1.
        return ACMODL(infor, fset, sset)

    def __init__(self, sid: int, form: str, scale: float,
                 real: Union[int, float],
                 imag: Union[int, float],
                 cid1: int, xyz: list[float],
                 cid2: int, nxyz: list[float], comment: str=''):
        """
        Creates a ACMODL card

        Parameters
        ----------

        """
        super().__init__()
        if comment:
            self.comment = comment

        self.sid = sid
        self.form = form
        self.scale = scale
        self.real = real
        self.imag = imag
        self.cid1 = cid1
        self.xyz = xyz
        self.cid2 = cid2
        self.nxyz = nxyz
        assert form in {'REAL', 'PHASE'}, form

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a ACPLNW card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        form = string_or_blank(card, 2, 'form', default='REAL')
        scale = double_or_blank(card, 3, 'scale', default=1.0)
        real = integer_double_or_blank(card, 4, 'real/mag', default=0.0)
        imag = integer_double_or_blank(card, 5, 'imag/phase', default=0.0)
        cid1 = integer(card, 9, 'cid1')
        xyz = [
            double_or_blank(card, 10, 'x', default=0.0),
            double_or_blank(card, 11, 'y', default=0.0),
            double_or_blank(card, 12, 'z', default=0.0),
        ]
        cid2 = integer(card, 13, 'cid2')
        nxyz = [
            double(card, 14, 'nx'),
            double(card, 15, 'ny'),
            double(card, 16, 'nz'),
        ]
        assert len(card) <= 17, f'len(ACMODL card) = {len(card):d}\ncard={card}'
        return ACPLNW(sid, form, scale,
                      real, imag,
                      cid1, xyz,
                      cid2, nxyz, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
    def safe_cross_reference(self, model: BDF, xref_error) -> None:
        pass
    def uncross_reference(self) -> None:
        pass

    def raw_fields(self):
        list_fields = ['ACPLNW', self.sid, self.form, self.scale,
                       self.real, self.imag, None, None, None,
                       self.cid1] + self.xyz + [self.cid2] + self.nxyz
        return list_fields

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)


class AMLREG(BaseCard):
    """
    +--------+---------+---------+--------------------------+
    | AMLREG |   RID   |   SID   | Name/Descriptor          |
    +--------+---------+---------+--------+--------+--------+
    |        | NLAYERS | RADSURF | INFID1 | INFID2 | INFID3 |
    +--------+---------+---------+--------+--------+--------+
    """
    type = 'AMLREG'

    @classmethod
    def _init_from_empty(cls):
        infor = 1
        fset = 1.0
        sset = 1.
        return AMLREG(infor, fset, sset)

    def __init__(self, rid: int, sid: int, name: str,
                 infid: list[int],
                 nlayers: int=5,
                 radsurf: str='AML', comment: str=''):
        """
        Creates a AMLREG card

        Parameters
        ----------
        rid: int
            AML region identification number
        sid : int
            Surface identification number
        name : str
            The Name/Descriptor is an optional character string
        infid : list[int]
            Identification number of an infinite plane.
            Up to three (3) infinite planes are considered for an AML
            region. The infinite planes are used when acoustic results
            are to be computed exterior to the AML region. (Integer
            >= 0; Default = 0) See Remark 5.
        nlayers : int; default=5
            Number of layers of extrusion to be formed by solver
        radsurf : str; default='AML'
           Radiation surface type
             1/AML:  the pressure and velocities on the AML boundary
                     are used to compute results in the far field.
             2/PHYB: the pressure and velocities on the physical
                     boundary (that is, all free fluid faces with the
                     exception of faces on the AML and the infinite
                     planes) are used to compute results in the far field.
             0/NONE: the region does not radiate.

        """
        super().__init__()
        if comment:
            self.comment = comment

        self.sid = sid
        self.rid = rid
        self.name = name
        self.infid = infid
        self.nlayers = nlayers
        self.radsurf = radsurf

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a ACPLNW card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        rid = integer(card, 1, 'rid')
        sid = integer(card, 2, 'sid')

        label_fields = [labeli for labeli in card[3:8] if labeli is not None]
        name = ''.join(label_fields).strip()
        assert len(name) <= 48, name

        nlayers = integer_or_blank(card, 9, 'nlayers', default=5)
        radsurf = string_or_blank(card, 10, 'radsurf', default='AML')
        infid = [
            integer(card, 11, 'infid1'),
            integer(card, 12, 'infid2'),
            integer(card, 13, 'infid3'),
        ]
        assert len(card) <= 14, f'len(AMLREG card) = {len(card):d}\ncard={card}'
        return AMLREG(rid, sid, name, infid,
                      nlayers=nlayers, radsurf=radsurf, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
    def safe_cross_reference(self, model: BDF, xref_error) -> None:
        pass
    def uncross_reference(self) -> None:
        pass

    def raw_fields(self):
        list_fields = ['AMLREG', self.rid, self.sid, self.name,
                       self.nlayers, self.radsurf] + self.infid
        return list_fields

    def write_card(self, size=8, is_double=False) -> str:
        infid = self.infid
        out = (
            f'AMLREG  {self.rid:<8d}{self.sid:<8d}{self.name:>s}\n'
            f'        {self.nlayers:<8d}{self.radsurf:8s}{infid[0]:8d}{infid[1]:8d}{infid[2]:8d}\n'
        )
        return out


class ACMODL(Element):
    """
    NX 2019.2
    +--------+-------+---------+----------+------+--------+--------+---------+----------+
    | ACMODL |       |  INFOR  |   FSET   | SSET | NORMAL |        | OVLPANG | SRCHUNIT |
    +--------+-------+---------+----------+------+--------+--------+---------+----------+
    |        | INTOL | AREAOP  |          |      | CTYPE  |        |         |          |
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
    def add_card(cls, card: BDFCard, nastran_version: str, comment: str=''):
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
        infor = string_or_blank(card, 2, 'infor', default='NONE')
        fset = integer_or_blank(card, 3, 'fset', default=None)
        sset = integer_or_blank(card, 4, 'sset', default=None)
        normal = double_or_blank(card, 5, 'normal', default=0.5)
        #
        olvpang = double_or_blank(card, 7, 'olvpang', default=60.0)
        search_unit = string_or_blank(card, 8, 'search_unit', default='REL')
        intol = double_or_blank(card, 9, 'intol', default=0.2)
        area_op = integer_or_blank(card, 10, 'area_op', default=0)
        ctype = string_or_blank(card, 13, 'ctype', default='STRONG')
        assert len(card) <= 14, f'len(ACMODL card) = {len(card):d}\ncard={card}'

        return ACMODL(infor, fset, sset,
                      normal=normal, olvpang=olvpang,
                      search_unit=search_unit, intol=intol,
                      area_op=area_op, ctype=ctype,
                      nastran_version='nx', comment=comment)

    @classmethod
    def add_card_msc(cls, card, comment=''):
        #print('MSC...ACMODL')
        inter = string_or_blank(card, 1, 'infor', default='DIFF')
        infor = string_or_blank(card, 2, 'infor', default='NONE')
        fset = integer_or_blank(card, 3, 'fset', default=None)
        sset = integer_or_blank(card, 4, 'sset', default=None)
        normal_default = 0.001 if inter == 'INDENT' else 1.0
        normal = double_or_blank(card, 5, 'normal', default=normal_default)
        method = string_or_blank(card, 6, 'method', 'BW') # BW/CP
        sk_neps = double_or_blank(card, 7, 'sk_neps', default=0.5)
        dsk_neps = double_or_blank(card, 8, 'dsk_neps', default=0.75)
        intol = double_or_blank(card, 9, 'intol', default=0.2)
        all_set = string_or_blank(card, 10, 'all_set', default='NO') #
        search_unit = string_or_blank(card, 11, 'search_unit', default='REL')
        assert len(card) <= 12, f'len(ACMODL card) = {len(card):d}\ncard={card}'
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
        infor = self.infor # None if self.infor == 'NONE' else self.infor
        fset = None if self.fset == 0 else self.fset
        sset = None if self.sset == 0 else self.sset
        if is_msc(self.nastran_version):
            list_fields = ['ACMODL', self.inter, infor,
                           fset, sset, self.normal, self.method,
                           self.sk_neps, self.dsk_neps, self.intol,
                           self.all_set, self.search_unit, ]
        else:
            list_fields = ['ACMODL', None, infor, fset, sset, self.normal, None, self.olvpang, self.search_unit,
                           self.intol, self.area_op, None, None, self.ctype]
            #list_fields = [
                #'ACMODL', self.pid, synth, self.tid_resistance, self.tid_reactance, self.tid_weight,
                #None, self.cutfr, self.b, self.k, self.m]
        return list_fields

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)


class PMIC(Property):
    type = 'PMIC'
    @classmethod
    def _init_from_empty(cls):
        pid = 1
        return PMIC(pid)

    def __init__(self, pid: int, comment: str=''):
        self.pid = pid
        self.comment = comment

    @classmethod
    def add_card(self, card: BDFCard, comment: str=''):
        pid = integer(card, 1, 'property_id')
        return PMIC(pid, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass
    def uncross_reference(self) -> None:
        return

    def raw_fields(self):
        return ['PMIC', self.pid]

    def write_card(self, size=8, is_double=False):
        fields = self.raw_fields()
        return print_card_8(fields)

class MATPOR(BaseCard):
    """

    model = CRAGG

    +--------+-----+-------+------+------+
    |   1    |  2  |   3   |   4  |  5   |
    +========+=====+=======+======+======+
    | MATPOR | MID | MODEL |  RHO |   C  |
    +--------+-----+-------+------+------+
    |        | RES |  POR  | TORT |      |
    +--------+-----+-------+------+------+

    model = JCA

    +--------+-----+-------+------+------+-------+-------+----+----+
    |   1    |  2  |   3   |   4  |  5   |   6   |   7   |  8 | 9  |
    +========+=====+=======+======+======+=======+=======+====+====+
    | MATPOR | MID | MODEL |  RHO |  C   | FRAME | GAMMA | PR | MU |
    +--------+-----+-------+------+------+-------+-------+----+----+
    |        | RES |  POR  | TORT | DENS |  L1   |       |    |    |
    +--------+-----+-------+------+------+-------+-------+----+----+
    """
    type = 'MATPOR'

    @classmethod
    def _init_from_empty(cls):
        mid = 1
        rho = 0.0
        c = 0.0
        resistivity = 0.0
        porosity = 0.0
        tortuosity = 0.0
        return cls.add_craggs(mid, rho, c, resistivity, porosity, tortuosity)

    def __init__(self, mid: int, model: str,
                 rho: float, c: float, resistivity: float,
                 porosity: float, tortuosity: float,
                 frame: str,
                 gamma: float, prandtl_number: float,
                 mu: float, L1, L2,
                 density: float=0.0, comment: str=''):
        """
        Creates a MATPOR card

        Parameters
        ----------
        """
        super().__init__()
        if comment:
            self.comment = comment

        self.mid = mid
        self.model = model
        self.rho = rho
        self.c = c
        self.resistivity = resistivity
        self.porosity = porosity
        self.tortuosity = tortuosity

        # JCA
        self.frame = frame
        self.gamma = gamma
        self.prandtl_number = prandtl_number
        self.mu = mu
        self.density = density
        self.L1 = L1
        self.L2 = L2

    @classmethod
    def add_craggs(self, mid: int,
                   rho: float, c: float, resistivity: float,
                   porosity: float, tortuosity: float,
                   comment: str=''):
        model = 'CRAGGS'
        frame = ''
        gamma = 0.0
        prandtl_number = 0.0
        mu = 0.0
        #density = 0.0
        L1 = 0.0
        L2 = 0.0
        return MATPOR(mid, model, rho, c, resistivity, porosity, tortuosity,
                      frame, gamma, prandtl_number, mu, L1, L2,
                      density=0.0, comment=comment)

    @classmethod
    def add_delmiki(self, mid: int,
                    rho: float, c: float, resistivity: float,
                    porosity: float, frame: str,
                    density: float=0.0, comment: str=''):
        model = 'CRAGGS'
        gamma = 0.0
        prandtl_number = 0.0
        mu = 0.0
        #density = 0.0
        L1 = 0.0
        L2 = 0.0
        return MATPOR(mid, model, rho, c, resistivity, porosity, tortuosity,
                      frame, gamma, prandtl_number, mu, L1, L2,
                      density=density, comment=comment)

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a MATPOR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #+--------+-----+-------+------+-----+--------+
        #| MATPOR | MID | MODEL |  RHO |  C  |        |
        #|        | RES |  POR  | TORT |     |        |
        #+--------+-----+-------+------+-----+--------+
        mid = integer(card, 1, 'material_id')
        model = string(card, 2, 'model')
        rho = double(card, 3, 'rho')
        c = double(card, 4, 'c')
        resistivity = double(card, 9, 'resistivity')
        porosity = double(card, 10, 'porosity')

        frame = ''
        gamma = 0.0
        prandlt_number = 0.0
        mu = 0.0
        density = 0.0
        L1 = 0.0
        L2 = 0.0
        tortuosity = 0.0
        if model == 'CRAGGS':
            tortuosity = double(card, 11, 'tortuosity')
            assert len(card) <= 12, f'len(MATPOR card) = {len(card):d}\ncard={card}'
        elif model == 'DELMIKI':
            # MATPOR MID MODEL RHO   C   FRAME
            #        RES POR        DENS
            frame = string(card, 5, 'frame')
            density = double_or_blank(card, 12, 'density', default=0.0)
            assert len(card) <= 13, f'len(MATPOR card) = {len(card):d}\ncard={card}'
        elif model == 'JCA':
            # MATPOR MID MODEL RHO  C    FRAME GAMMA PR MU
            #        RES  POR  TORT DENS   L1   L2
            frame = string(card, 5, 'frame')
            gamma = double(card, 6, 'gamma')
            prandlt_number = double(card, 7, 'prandlt_number')
            mu = double(card, 8, 'mu_dynamic_visc')
            density = double_or_blank(card, 12, 'density', default=0.0)
            tortuosity = double(card, 11, 'tortuosity')
            L1 = double(card, 13, 'L1')
            L2 = double(card, 14, 'L2')
            assert len(card) <= 15, f'len(MATPOR card) = {len(card):d}\ncard={card}'
        else:  # pragma: no cover
            raise RuntimeError(model)
        return MATPOR(mid, model, rho, c, resistivity, porosity, tortuosity,
                      frame, gamma, prandlt_number,
                      mu, L1, L2, density=density)

    def cross_reference(self, model: BDF) -> None:
        pass
    def safe_cross_reference(self, model: BDF, xref_error) -> None:
        pass
    def uncross_reference(self) -> None:
        pass

    def raw_fields(self):
        if self.model == 'CRAGGS':
            list_fields = [
                'MATPOR', self.mid, self.model, self.rho, self.c,
                None, None, None, None,
                self.resistivity, self.porosity, self.tortuosity]
        elif self.model == 'DELMIKI':
            list_fields = [
                'MATPOR', self.mid, self.model, self.rho, self.c,
                self.frame, self.gamma, None, None,
                self.resistivity, self.porosity, None,
                self.density]

        elif self.model == 'JCA':
            # MATPOR MID MODEL RHO  C    FRAME GAMMA PR MU
            #        RES  POR  TORT DENS   L1   L2
            list_fields = [
                'MATPOR', self.mid, self.model, self.rho, self.c,
                self.frame, self.gamma, self.prandtl_number, self.mu,
                self.resistivity, self.porosity, self.tortuosity,
                self.density, self.L1, self.L2]
        else:  # pragma: no cover
            raise NotImplementedError(self.model)
        return list_fields

    def write_card(self, size=8, is_double=False) -> str:
        fields = self.raw_fields()
        return print_card_8(fields)

