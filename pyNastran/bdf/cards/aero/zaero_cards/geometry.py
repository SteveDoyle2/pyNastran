from __future__ import annotations
import numpy as np
from typing import Optional, Any, TYPE_CHECKING


from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, string,
    integer_or_string, double,
)
from pyNastran.bdf.cards.aero.utils import (
    elements_from_quad,
    points_elements_from_quad_points,
    create_ellipse,
    #create_axisymmetric_body,
)
from pyNastran.bdf.cards.aero.aero import (
    Spline, CAERO1, CAERO2, PAERO2, AEFACT, AELIST, AESURF)
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    import matplotlib
    AxesSubplot = matplotlib.axes._subplots.AxesSubplot


class PANLST1(Spline):
    """
    Defines a set of aerodynamic boxes by the LABEL entry in CAERO7 or BODY7
    bulk data cards.

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE1 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE1 | 100  |       |       |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+

    +---------+-------+---------+------+------+------+----+-----+-------+
    |    1    |   2   |    3    |   4  |  5   |   6  |  7 |  8  |   9   |
    +=========+=======+=========+======+======+======+====+=====+=======+
    | PANLST1 | SETID | MACROID | BOX1 | BOX2 |      |    |     |       |
    +---------+-------+---------+------+------+------+----+-----+-------+
    | PANLST1 |  100  |   111   |  111 |  118 |      |    |     |       |
    +---------+-------+---------+------+------+------+----+-----+-------+

    PANLST1 is referred to by SPLINEi, ATTACH, LOADMOD, CPFACT, JETFRC, and/or
    AESURFZ bulk data card.
    """
    type = 'PANLST1'

    def __init__(self, eid: int, macro_id: int, box1: int, box2: int,
                 comment: str=''):
        """
        Creates a PANLST1 card

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.macro_id = macro_id  # points to CAERO7 / BODY7
        self.box1 = box1
        self.box2 = box2
        self.aero_element_ids = []
        self.caero_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PANLST3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        macro_id = integer(card, 2, 'macro_id')
        box1 = integer(card, 3, 'box1')
        box2 = integer(card, 4, 'box2')
        assert len(card) == 5, f'len(PANLST1 card) = {len(card):d}\ncard={card}'
        return PANLST1(eid, macro_id, box1, box2, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by PANLST1 eid=%s' % self.eid
        self.caero_ref = model.CAero(self.macro_id, msg=msg)
        self.aero_element_ids = np.arange(self.box1, self.box2)

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        self.caero_ref = None

    def raw_fields(self):
        list_fields = ['PANLST1', self.eid, self.macro_id, self.box1, self.box2]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PANLST2(Spline):
    """
    Defines a set of aerodynamic boxes by the LABEL entry in CAERO7 or BODY7
    bulk data cards.

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE1 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE1 | 100  |       |       |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+

    +---------+-------+---------+------+------+------+----+-----+-------+
    |    1    |   2   |    3    |   4  |  5   |   6  |  7 |  8  |   9   |
    +=========+=======+=========+======+======+======+====+=====+=======+
    | PANLST1 | SETID | MACROID | BOX1 | BOX2 |      |    |     |       |
    +---------+-------+---------+------+------+------+----+-----+-------+
    | PANLST1 |  100  |   111   |  111 |  118 |      |    |     |       |
    +---------+-------+---------+------+------+------+----+-----+-------+

    PANLST1 is referred to by SPLINEi, ATTACH, LOADMOD, CPFACT, JETFRC, and/or
    AESURFZ bulk data card.
    """
    type = 'PANLST2'

    def __init__(self, eid: int, macro_id: int, boxs: list[int],
                 comment: str=''):
        """
        Creates a PANLST2 card

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.2_Users_3rd_Ed.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.macro_id = macro_id  # points to CAERO7 / BODY7
        self.boxs = boxs
        self.aero_element_ids = []
        self.caero_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PANLST2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # SETID   MACROID BOX1    BOX2    ETC
        eid = integer(card, 1, 'eid')
        macro_id = integer(card, 2, 'macro_id')

        box1 = integer(card, 3, 'box1')
        thru = integer_or_string(card, 4, 'thru')
        if isinstance(thru, str):
            assert thru == 'THRU', thru
            box2 = integer(card, 5, 'box2')
            boxs = [box1, 'THRU', box2]
            assert len(card) == 6, f'len(PANLST2 card) = {len(card):d}\ncard={card}'
        else:
            boxs = [box1, thru]
            for ifield in range(5, len(card)):
                box2 = integer(card, ifield, 'box2')
                boxs.append(box2)
        return PANLST2(eid, macro_id, boxs, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by PANLST2 eid=%s' % self.eid
        self.caero_ref = model.CAero(self.macro_id, msg=msg)
        if 'THRU' in self.boxs:
            box1, thru, box2 = self.boxs
            # +1 leads to [box1, box2] instead of [box1, box2)
            boxs = np.arange(box1, box2+1)
        else:
            boxs = np.asarray(self.boxs, dtype='int32')
        self.aero_element_ids = boxs

    def uncross_reference(self) -> None:
        self.caero_ref = None

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def raw_fields(self):
        list_fields = ['PANLST2', self.eid, self.macro_id] + list(self.boxs)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PANLST3(Spline):
    """
    Defines a set of aerodynamic boxes by the LABEL entry in CAERO7 or BODY7
    bulk data cards.

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE1 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE1 | 100  |       |       |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+

    +---------+-------+--------+--------+--------+-----+----+-----+-------+
    |    1    |   2   |   3    |    4   |    5   |  6  |  7 |  8  |   9   |
    +=========+=======+========+========+========+=====+====+=====+=======+
    | PANLST3 | SETID | LABEL1 | LABEL2 | LABEL3 | etc |    |     |       |
    +---------+-------+--------+--------+--------+-----+----+-----+-------+
    | PANLST3 | 100   | WING   | HTAIL  |        |     |    |     |       |
    +---------+-------+--------+--------+--------+-----+----+-----+-------+

    PANLST3 is referred to by SPLINEi, ATTACH, LOADMOD, CPFACT, JETFRC, and/or
    AESURFZ bulk data card.

    """
    type = 'PANLST3'

    def __init__(self, eid: int, panel_groups, comment: str=''):
        """
        Creates a PANLST3 card

        Parameters
        ----------
        eid : int
            spline id
        comment : str; default=''
            a comment for the card

        """
        # https://www.zonatech.com/Documentation/ZAERO_9.3_Users_Full_Electronic.pdf
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.panel_groups = panel_groups  # points to CAERO7 / BODY7
        self.aero_element_ids = []
        self.caero_refs = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PANLST3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        group_id = 1
        panel_groups = []
        for ifield in range(2, len(card)):
            name = string(card, ifield, f'group_{group_id:d}')
            panel_groups.append(name)
        assert len(card) > 2, 'len(PANLST3 card)=%d; no panel_groups were defined\ncard=%s' % (len(card), card)
        return PANLST3(eid, panel_groups, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by PANLST3 eid=%d' % self.eid
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        caero_refs = []
        aero_element_ids = []
        for caero_label in self.panel_groups:
            try:
                caero_eid = model.zaero.caero_to_name_map[caero_label]
            except KeyError:
                keys = list(model.zaero.caero_to_name_map)
                keys.sort()
                caero_labels = [caero.label for caero in model.caeros.values()]
                caero_labels.sort()
                msg = (
                    f'PANLST3 eid={self.eid}: Missing label={caero_label!r}\n'
                    f' - keys        ={keys}\n'
                    f' - caero_labels={caero_labels}'
                )
                raise RuntimeError(msg)
            caero_ref = model.CAero(caero_eid, msg=msg)
            caero_refs.append(caero_ref)
            eid = caero_ref.eid
            npanels = caero_ref.npanels
            if npanels == 0:
                model.log.warning('skipping PANLST3 because there are 0 panels in:\n%r' % caero_ref)
                continue
            aero_element_ids2 = range(eid, eid + npanels)
            assert len(aero_element_ids2) == npanels, npanels
            aero_element_ids += aero_element_ids2
        self.caero_refs = caero_refs
        self.aero_element_ids = aero_element_ids

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        self.caero_refs = None
        # self.aero_element_ids = aero_element_ids

    def raw_fields(self):
        list_fields = ['PANLST3', self.eid] + self.panel_groups
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PAFOIL7(BaseCard):
    """
    Defines an aerodynamic body macroelement of a body-like component.
    Similar to Nastran's CAERO2.

    +---------+----+------+------+-------+------+------+-------+------+
    |   1     |  2 |   3  |  4   |   5   |   6  |   7  |   8   |  9   |
    +=========+====+======+======+=======+======+======+=======+======+
    | PAFOIL7 | ID | ITAX | ITHR | ICAMR | RADR | ITHT | ICAMT | RADT |
    +---------+----+------+------+-------+------+------+-------+------+
    | PAFOIL7 |  1 | -201 |  202 |  203  |  0.1 |  211 |  212  |  0.1 |
    +---------+----+------+------+-------+------+------+-------+------+

    """
    type = 'PAFOIL7'

    def __init__(self, pid: int, i_axial: int,
                 i_thickness_root: int, i_camber_root: int, le_radius_root: float,
                 i_thickness_tip: int, i_camber_tip: int, le_radius_tip: float,
                 comment: str=''):
        """
        Defines a BODY7 card, which defines a slender body
        (e.g., fuselage/wingtip tank).

        Parameters
        ----------
        pid : int
            PAFOIL7 identification number.
        i_axial : str
            Identification number of an AEFACT bulk data card used to
            specify the xcoordinate locations, in percentage of the
            chord length, where the thickness and camber are specified.
            ITAX can be a negative number (where ABS (ITAX) = AEFACT
            bulk data card identification number) to request linear
            interpolation.
        i_thickness_root / i_thickness_tip : int
            Identification number of an AEFACT bulk data card used to
            specify the half thickness of the airfoil at the wing
            root/tip.
        le_radius_root / le_radius_root: float
            Leading edge radius at the root/tip normalized by the
            root/tip chord.
        i_thickness_tip : int
            Identification number of an AEFACT bulk data card used to
            specify the half thickness at the wing tip.
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        self.pid = pid
        self.i_axial = i_axial

        self.i_thickness_root = i_thickness_root
        self.i_camber_root = i_camber_root
        self.le_radius_root = le_radius_root

        self.i_camber_tip = i_camber_tip
        self.le_radius_tip = le_radius_tip
        self.i_thickness_tip = i_thickness_tip

        self.i_thickness_root_ref = None
        self.i_camber_root_ref = None
        self.i_thickness_tip_ref = None
        self.i_camber_tip_ref = None
        self.i_axial_ref = None

    #@property
    #def cp(self) -> int:
        #return self.acoord
    #@property
    #def cp_ref(self):
        #return self.acoord_ref

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PAFOIL7 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """

        pid = integer(card, 1, 'pid')
        i_axial = integer(card, 2, 'i_axial')

        i_thickness_root = integer(card, 3, 'i_thickness_root')
        i_camber_root = integer(card, 4, 'i_camber_root')
        le_radius_root = double_or_blank(card, 5, 'le_radius_root')

        i_thickness_tip = integer(card, 6, 'i_thickness_tip')
        i_camber_tip = integer(card, 7, 'i_camber_tip')
        le_radius_tip = double_or_blank(card, 8, 'le_radius_tip')

        assert len(card) <= 9, f'len(PAFOIL7 card) = {len(card):d}\ncard={card}'
        return PAFOIL7(pid, i_axial,
                       i_thickness_root, i_camber_root, le_radius_root,
                       i_thickness_tip, i_camber_tip, le_radius_tip,
                       comment=comment)

    #def ACoord(self):
        #if self.acoord_ref is not None:
            #return self.acoord_ref.cid
        #return self.acoord

    #def Pid(self):
        #if self.pid_ref is not None:
            #return self.pid_ref.pid
        #return self.pid

    #def Lsb(self):  # AEFACT
        #if self.lsb_ref is not None:
            #return self.lsb_ref.sid
        #return self.lsb

    #def Lint(self):  # AEFACT
        #if self.lint_ref is not None:
            #return self.lint_ref.sid
        #return self.lint

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PAFOIL7 pid=%s' % self.pid
        self.i_axial_ref = model.AEFact(abs(self.i_axial), msg=msg)

        self.i_thickness_root_ref = model.AEFact(self.i_thickness_root, msg=msg)
        self.i_camber_root_ref = model.AEFact(self.i_camber_root, msg=msg)

        self.i_thickness_tip_ref = model.AEFact(self.i_thickness_tip, msg=msg)
        self.i_camber_tip_ref = model.AEFact(self.i_camber_tip, msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.i_thickness_root_ref = None
        self.i_camber_root_ref = None
        self.i_thickness_tip_ref = None
        self.i_camber_tip_ref = None

    def convert_to_nastran(self, model):
        """
        Should this be converted to a DMIG?

        +---------+----+------+------+-------+------+------+-------+------+
        |   1     |  2 |   3  |  4   |   5   |   6  |   7  |   8   |  9   |
        +=========+====+======+======+=======+======+======+=======+======+
        | PAFOIL7 | ID | ITAX | ITHR | ICAMR | RADR | ITHT | ICAMT | RADT |
        +---------+----+------+------+-------+------+------+-------+------+
        | PAFOIL7 |  1 | -201 |  202 |  203  |  0.1 |  211 |  212  |  0.1 |
        +---------+----+------+------+-------+------+------+-------+------+

        """
        raise NotImplementedError('PAFOIL7: convert_to_nastran')

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        #pid = integer(card, 1, 'pid')
        #i_axial = integer(card, 2, 'i_axial')

        #i_thickness_root = integer(card, 3, 'i_thickness_root')
        #i_camber_root = integer(card, 4, 'i_camber_root')
        #le_radius_root = double_or_blank(card, 5, 'le_radius_root')

        #i_thickness_tip = integer(card, 6, 'i_thickness_tip')
        #le_radius_tip = integer(card, 7, 'le_radius_tip')
        #i_camber_tip = double_or_blank(card, 8, 'i_camber_tip')

        list_fields = [
            'PAFOIL7', self.pid, self.i_axial,
            self.i_thickness_root, self.i_camber_root, self.le_radius_root,
            self.i_thickness_tip, self.i_camber_tip, self.le_radius_tip,
        ]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class BODY7(BaseCard):
    """
    Defines an aerodynamic body macroelement of a body-like component.
    Similar to Nastran's CAERO2.

    +--------+-----+-----+----+-----+------+-----+------+------+
    |    1   |  2  |  3  |  4 |  5  |   6  |   7 |   8  |  9   |
    +========+=====+=====+====+=====+======+=====+======+======+
    | CAERO2 | EID | PID | CP | NSB | NINT | LSB | LINT | IGID |
    +--------+-----+-----+----+-----+------+-----+------+------+
    |        | X1  |  Y1 | Z1 | X12 |      |     |      |      |
    +--------+-----+-----+----+-----+------+-----+------+------+

    +-------+---------+-------+---------+--------+------+---------+---------+---------+
    |   1   |    2    |   3   |     4   |    5   |   6  |     7   |    8    |    9    |
    +=======+=========+=======+=========+========+======+=========+=========+=========+
    | BODY7 |   BID   | LABEL | IPBODY7 | ACOORD | NSEG | IDMESH1 | IDMESH2 | IDMESH3 |
    +-------+---------+-------+---------+--------+------+---------+---------+---------+
    |       | IDMESH4 |  etc  |         |        |      |         |         |         |
    +-------+---------+-------+---------+--------+------+---------+---------+---------+
    | BODY7 |    4    |  BODY |    2    |    8   |   4  |    20   |    21   |    22   |
    +-------+---------+-------+---------+--------+------+---------+---------+---------+
    |       |    23   |       |         |        |      |         |         |         |
    +-------+---------+-------+---------+--------+------+---------+---------+---------+

    """
    type = 'BODY7'

    def __init__(self, eid: int, label: str, pid: int, nseg: int,
                 idmeshes: list[int], acoord: int=0, comment: str=''):
        """
        Defines a BODY7 card, which defines a slender body
        (e.g., fuselage/wingtip tank).

        Parameters
        ----------
        eid : int
            body id
        label : str
            An arbitrary character string used to define the body.
        pid  : int; default=0
            Identification number of PBODY7 bulk data card
            (specifying body wake and/or inlet aerodynamic boxes)
        nseg : int
            Number of body segments
        idmeshes : list[int]
            Identification number of SEGMESH bulk data card (specifying body segments).
        acoord : int; default=0
            Identification number of ACOORD bulk data card
            (specifying body center line location and orientation)
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        #: Element identification number
        self.eid = eid
        self.label = label

        #: Property identification number of a PAERO7 entry.
        self.pid = pid
        self.nseg = nseg
        self.idmeshes = idmeshes
        self.acoord = acoord
        self.pid_ref = None
        self.acoord_ref = None
        self.ascid_ref = None
        self.segmesh_refs = None

    #@property
    #def cp(self):
        #return self.acoord
    #@property
    #def cp_ref(self):
        #return self.acoord_ref

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a BODY7 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        label = string(card, 2, 'label')
        assert len(card) >= 3, f'len(BODY7 card) = {len(card):d}\ncard={card}'
        pid = integer_or_blank(card, 3, 'pid')
        acoord = integer_or_blank(card, 4, 'acoord', default=0)
        nseg = integer_or_blank(card, 5, 'nseg')

        idmeshes = []
        for i, ifield in enumerate(range(6, len(card))):
            segmesh = integer(card, ifield, 'idmesh_%i' % (i+1))
            idmeshes.append(segmesh)
        assert len(card) <= 13, f'len(BODY7 card) = {len(card):d}\ncard={card}'
        return BODY7(eid, label, pid, nseg, idmeshes, acoord=acoord, comment=comment)

    def ACoord(self):
        if self.acoord_ref is not None:
            return self.acoord_ref.cid
        return self.acoord

    def Pid(self):
        if self.pid_ref is not None:
            return self.pid_ref.pid
        return self.pid

    @property
    def nboxes(self):
        if self.nsb > 0:
            return self.nsb
        return len(self.lsb_ref.fractions) # AEFACT

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by BODY7 eid=%s' % self.eid
        self.segmesh_refs = []
        for segmesh_id in self.idmeshes:
            segmesh_ref = model.PAero(segmesh_id, msg=msg)  # links to SEGMESH/PAERO7
            self.segmesh_refs.append(segmesh_ref)

        #if self.pid is not None:
            #self.pid_ref = model.PAero(self.pid, msg=msg)  # links to PAERO7
        #self.acoord_ref = model.Coord(self.acoord, msg=msg)
        #if self.nsb == 0:
            #self.lsb_ref = model.AEFact(self.lsb, msg=msg)
        #if self.nint == 0:
            #self.lint_ref = model.AEFact(self.lint, msg=msg)
        if self.acoord is not None:
            self.acoord_ref = model.Coord(self.acoord, msg=msg)
        #self.ascid_ref = model.Acsid(msg=msg)
        self.ascid_ref = model.Coord(0, msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        self.acoord = self.ACoord()
        self.pid_ref = None
        self.acoord_ref = None

    def convert_to_nastran(self, model):
        """
        +--------+-----+-----+----+-----+------+-----+------+------+
        |    1   |  2  |  3  |  4 |  5  |   6  |   7 |   8  |  9   |
        +========+=====+=====+====+=====+======+=====+======+======+
        | CAERO2 | EID | PID | CP | NSB | NINT | LSB | LINT | IGID |
        +--------+-----+-----+----+-----+------+-----+------+------+
        |        | X1  |  Y1 | Z1 | X12 |      |     |      |      |
        +--------+-----+-----+----+-----+------+-----+------+------+

        +-------+---------+-------+---------+--------+------+---------+---------+---------+
        |   1   |    2    |   3   |     4   |    5   |   6  |     7   |    8    |    9    |
        +=======+=========+=======+=========+========+======+=========+=========+=========+
        | BODY7 |   BID   | LABEL | IPBODY7 | ACOORD | NSEG | IDMESH1 | IDMESH2 | IDMESH3 |
        +-------+---------+-------+---------+--------+------+---------+---------+---------+
        |       | IDMESH4 |  etc  |         |        |      |         |         |         |
        +-------+---------+-------+---------+--------+------+---------+---------+---------+
        | BODY7 |    4    |  BODY |    2    |    8   |   4  |    20   |    21   |    22   |
        +-------+---------+-------+---------+--------+------+---------+---------+---------+
        |       |    23   |       |         |        |      |         |         |         |
        +-------+---------+-------+---------+--------+------+---------+---------+---------+

        """
        pid = max(model.paeros) + 1000
        igroup = 1
        orient = 'ZY'
        cp = 0
        #width = 1

        #nsb  : AEFACT id for defining the location of the slender body elements
        #lsb  : AEFACT id for defining the location of interference elements
        #nint : Number of slender body elements
        #lint : Number of interference elements
        aefact_id = len(model.aefacts) + 1
        xs_id = aefact_id
        half_width_id = aefact_id + 1
        theta1_id = aefact_id + 2
        theta2_id = aefact_id + 3

        lsb = xs_id
        lint = xs_id

        #+---------+--------+-------+------+---------+-------+------+------+-----+
        #|    1    |    2   |   3   |   4  |    5    |   6   |   7  |   8  |  9  |
        #+=========+========+=======+======+=========+=======+======+======+=====+
        #| SEGMESH | IDMESH | NAXIS | NRAD | NOSERAD | IAXIS |      |      |     |
        #|         | ITYPE1 |   X1  | CAM1 |   YR1   |  ZR1  | IDY1 | IDZ1 |     |
        #|         | ITYPE2 |   X2  | CAM2 |   YR2   |  ZR2  | IDY2 | IDZ2 |     |
        #|         | ITYPE3 |   X3  | CAM3 |   YR3   |  ZR3  | IDY3 | IDZ3 |     |
        #+---------+--------+-------+------+---------+-------+------+------+-----+
        xpoints = []
        half_widths = []

        try:
            origin_x, origin_y, origin_z = self.acoord_ref.origin
        except AttributeError:  # pragma: no cover
            print(self.get_stats())
            raise
        #x_offset = origin_x + x
        #y_offset = origin_y + y
        #z_offset = origin_z + z
        nsegmesh = len(self.segmesh_refs)
        if nsegmesh == 0:
            raise RuntimeError('Number of SEGMESH  references on BODY7=0\n%s' % str(self))

        for isegmesh, segmesh in enumerate(self.segmesh_refs):
            itypes = segmesh.itypes

            #xs = segmesh.xs
            #idys_ref = segmesh.idys_ref
            #idzs_ref = segmesh.idzs_ref
            nitypes = len(itypes)
            idys_ref = [None] * nitypes if segmesh.idys_ref is None else segmesh.idys_ref
            idzs_ref = [None] * nitypes if segmesh.idzs_ref is None else segmesh.idzs_ref

            cambers = segmesh.cambers
            yrads = segmesh.ys
            zrads = segmesh.zs

            # what????
            if isegmesh in [0, nsegmesh - 1]:
                xs2 = segmesh.xs
                idys_ref2 = idys_ref
                idzs_ref2 = idzs_ref
            else:
                xs2 = segmesh.xs[1:]
                idys_ref2 = idys_ref[1:]
                idzs_ref2 = idzs_ref[1:]

            xpoints += xs2
            yz_mean = []

            thetas = self._get_thetas()
            for itype, camber, yrad, zrad, idy_ref, idz_ref in zip(
                    itypes, cambers, yrads, zrads, idys_ref2, idzs_ref2):
                out = self._get_body7_width_height_radius(
                    thetas, itype, camber, yrad, zrad, idy_ref, idz_ref)
                width, height, average_radius, ymeani, zmeani = out

                half_widths.append(average_radius)
                yz_mean.append([ymeani, zmeani])

            # I think you could area weight this and get a better mean...
            ymean, zmean = np.mean(yz_mean, axis=0)

        xpoints_local = [xi for xi in xpoints]
        assert len(half_widths) == len(xpoints_local)

        half_width = max(half_widths)
        AR = 1.0

        p1 = [origin_x, origin_y + ymean, origin_z + zmean]
        x12 = max(xpoints) - min(xpoints)
        dash = '-' * 80 + '\n'
        comment = dash
        comment += self.comment
        caero2 = CAERO2(self.eid, pid, igroup, p1, x12,
                        cp=cp,
                        nsb=0, lsb=lsb,
                        nint=0, lint=lint, comment=comment)

        #
        lrsb = half_width_id  # slender body
        lrib = half_width_id  # interference

        # theta arrays (AEFACTs)
        #lth1 = theta1_id
        #lth2 = theta2_id

        # 0-57 is excluded
        angles_body = [20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0,
                       200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0]
        angles_fin = [20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0,
                      200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0]

        aefact_xs = AEFACT(xs_id, xpoints_local, comment=dash+'Xs')
        aefact_width = AEFACT(half_width_id, half_widths, comment='half_widths')
        aefact_theta1 = AEFACT(theta1_id, angles_body, comment='angles_body')
        aefact_theta2 = AEFACT(theta2_id, angles_fin, comment='angles_fin')

        # which segments use theta1 array
        lth = [1, 10]  # nsegments] # t
        thi = [1]
        thn = [1]
        paero2 = PAERO2(pid, orient, half_width, AR, thi, thn,
                        lrsb=lrsb, lrib=lrib,
                        lth=lth, comment='')
        caero2.validate()
        paero2.validate()
        return caero2, paero2, aefact_xs, aefact_width, aefact_theta1, aefact_theta2

    def _get_body7_width_height_radius(self, thetas: np.ndarray,
                                       itype: int, camber: float,
                                       yrad: float, zrad: float,
                                       idy_ref, idz_ref) -> tuple[float, float, float, float, float]:
        if itype == 1:
            # Body of Revolution
            # Xi, CAMi, YRi
            radius = yrad
            aspect_ratio = 1.
            yz = create_ellipse(aspect_ratio, radius, thetas=thetas)
            ypoints = yz[:, 0]
            zpoints = camber + yz[:, 1]
        elif itype == 2:
            # Elliptical body
            height = zrad
            width = yrad
            aspect_ratio = height / width
            radius = height
            yz = create_ellipse(aspect_ratio, radius, thetas=thetas)
            ypoints = yz[:, 0]
            zpoints = yz[:, 1]

        elif itype == 3:
            # Arbitrary body using AEFACTss
            try:
                ypoints = idy_ref.fractions
                zpoints = idz_ref.fractions
            except AttributeError:  # pragma: no cover
                print('idy_ref = %s' % idy_ref)
                print('idz_ref = %s' % idz_ref)
                print(self.get_stats())
                raise
        else:  # pramga: no cover
            msg = f'Unsupported itype={itype} (must be 1/2/3)\n{str(self)}'
            raise NotImplementedError(msg)
        width = ypoints.max() - ypoints.min()
        height = zpoints.max() - zpoints.min()
        average_radius = (width + height) / 4.
        #elliptical_area = pi * width * height
        ymeani = ypoints.mean()
        zmeani = zpoints.mean()
        return width, height, average_radius, ymeani, zmeani

    def _get_nthetas(self) -> int:
        """gets the number of thetas for the body"""
        return self.segmesh_refs[0].nradial  # npoints
        #nthetas = 17
        #for itype, idy_ref, unused_idz_ref in zip(itypes, idys_ref2, idzs_ref2):
            #if itype == 3:
                #fractions = idy_ref.fractions
                #nthetas = len(fractions)
                #break
        #return nthetas

    def _get_thetas(self) -> np.ndarray:
        """gets the thetas for the body"""
        nthetas = self._get_nthetas()
        thetas = np.radians(np.linspace(0., 360., nthetas))
        return thetas

    def get_points(self) -> list[np.ndarray, np.ndarray]:
        """creates a 1D representation of the BODY7"""
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))

        #print("x12 = %s" % self.x12)
        #print("pcaero[%s] = %s" % (self.eid, [p1,p2]))
        return [p1, p2]

    @property
    def npanels(self) -> int:
        """gets the number of panels for the body"""
        if self.segmesh_refs is None:
            raise RuntimeError('xref the BODY7')
        nz = len(self.segmesh_refs)
        unused_segmesh = self.segmesh_refs[0]
        nthetas = self._get_nthetas()
        npanels = nz * (nthetas - 1)
        return npanels

    def get_points_elements_3d(self):
        """
        Gets the points/elements in 3d space as CQUAD4s
        The idea is that this is used by the GUI to display CAERO panels.

        TODO: doesn't support the aero coordinate system

        """
        #paero2 = self.pid_ref
        xyz = []
        element = []
        npoints = 0
        for segmesh in self.segmesh_refs:
            #print(segmesh)
            xyzi, elementi = self._get_points_elements_3di(segmesh)
            xyz.append(xyzi)
            element.append(elementi + npoints)
            npoints += xyzi.shape[0]

        xyzs = np.vstack(xyz)
        elements = np.vstack(element)
        assert xyzs is not None, str(self)
        assert elements is not None, str(self)

        return xyzs, elements

    def _get_points_elements_3di(self, segmesh: SEGMESH) -> tuple[np.ndarray, np.ndarray]:
        """
        points (nchord, nspan) float ndarray; might be backwards???
            the points
        elements (nquads, 4) int ndarray
            series of quad elements
            nquads = (nchord-1) * (nspan-1)
        """
        #lengths_y = []
        #lengths_z = []
        nx = segmesh.naxial
        ny = segmesh.nradial

        xs = []
        ys = []
        zs = []

        origin_x, origin_y, origin_z = self.acoord_ref.origin

        nthetas = segmesh.nradial
        thetas = np.radians(np.linspace(0., 360., nthetas))
        for itype, x, yrad, zrad, camber, idy_ref, idz_ref in zip(
                segmesh.itypes,
                segmesh.xs, segmesh.ys, segmesh.zs,
                segmesh.cambers,
                segmesh.idys_ref, segmesh.idzs_ref):

            xsi, ysi, zsi = self._get_xyzs_offset(
                origin_x, origin_y, origin_z, thetas,
                itype, x, yrad, zrad, camber, idy_ref, idz_ref)
            xs.append(xsi)
            ys.append(ysi)
            zs.append(zsi)

        xyz = np.vstack([
            np.hstack(xs),
            np.hstack(ys),
            np.hstack(zs),
        ]).T
        elements = elements_from_quad(nx, ny, dtype='int32')  # nx,ny are points
        return xyz, elements

    def _get_xyzs_offset(self, origin_x, origin_y, origin_z, thetas,
                         itype: int, x: float, yrad: float, zrad: float, camber: float,
                         idy_ref, idz_ref) -> tuple[list[float], np.ndarray, np.ndarray]:
        y = 0.
        z = 0.
        if itype == 1:
            # Body of Revolution
            # Xi, CAMi, YRi
            ## TODO: doesn't consider camber
            radius = yrad
            aspect_ratio = 1.
            yz = create_ellipse(aspect_ratio, radius, thetas=thetas)
            ypoints = yz[:, 0]
            zpoints = camber + yz[:, 1]
        elif itype == 2:
            # Elliptical body
            #  Xi, YRi, ZRi
            height = zrad
            width = yrad
            aspect_ratio = height / width
            radius = height
            yz = create_ellipse(aspect_ratio, radius, thetas=thetas)
            ypoints = yz[:, 0]
            zpoints = yz[:, 1]

        elif itype == 3:
            # Arbitrary body using AEFACTss
            #  Xi, IDYi, IDZi
            ypoints = idy_ref.fractions
            zpoints = idz_ref.fractions
            y = yrad
            z = zrad
        else:  # pramga: no cover
            msg = 'Unsupported itype=%s (must be 1/2/3)\n%s' % (itype, str(self))
            raise NotImplementedError(msg)

        assert len(ypoints) == len(zpoints), 'len(ypoints)=%s len(zpoints)=%s' % (len(ypoints), len(zpoints))
        nnodes = len(ypoints)

        x_offset = origin_x + x
        y_offset = origin_y + y
        z_offset = origin_z + z

        xsi = [x_offset] * nnodes
        ysi = y_offset + ypoints
        zsi = z_offset + zpoints
        return xsi, ysi, zsi

    #def set_points(self, points):
        #self.p1 = np.asarray(points[0])
        #p2 = np.asarray(points[1])
        #x12 = p2 - self.p1
        #self.x12 = x12[0]

    #def shift(self, dxyz):
        #"""shifts the aero panel"""
        #self.p1 += dxyz

    def raw_fields(self) -> list[Any]:
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        list_fields = ['BODY7', self.eid, self.label, self.Pid(), self.ACoord(),
                       self.nseg] + self.idmeshes
        return list_fields

    def repr_fields(self) -> list[Any]:
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SEGMESH(BaseCard):
    """
    Defines a grid system for a body segment; referenced by the BODY7 bulk data card.

    +---------+--------+-------+------+---------+-------+------+------+-----+
    |    1    |    2   |   3   |   4  |    5    |   6   |   7  |   8  |  9  |
    +=========+========+=======+======+=========+=======+======+======+=====+
    | SEGMESH | IDMESH | NAXIS | NRAD | NOSERAD | IAXIS |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         | ITYPE1 |   X1  | CAM1 |   YR1   |  ZR1  | IDY1 | IDZ1 |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         | ITYPE2 |   X2  | CAM2 |   YR2   |  ZR2  | IDY2 | IDZ2 |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         | ITYPE3 |   X3  | CAM3 |   YR3   |  ZR3  | IDY3 | IDZ3 |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    | SEGMESH |    2   |  3    |  6   |         |       |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         |    1   |  0.0  |  0.0 | 0.0     |       |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         |    1   |  1.0  |  0.0 | 0.5     |       |      |      |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+
    |         |    3   |  2.0  |      |         |       |  103 | 104  |     |
    +---------+--------+-------+------+---------+-------+------+------+-----+

    """
    type = 'SEGMESH'

    @property
    def pid(self) -> int:
        return self.segmesh_id

    @pid.setter
    def pid(self, segmesh_id: int) -> None:
        self.segmesh_id = segmesh_id

    def __init__(self, segmesh_id: int, naxial: int, nradial: int,
                 nose_radius: float, iaxis: int,
                 itypes: list[int],
                 xs: list[float], cambers: list[float],
                 ys: list[float], zs: list[float],
                 idys: list[int], idzs: list[int], comment: str=''):
        """
        Defines a SEGMESH card, which defines a cross-section for a PBODY7.

        Parameters
        ----------
        segmesh_id : int
            Body segment mesh identification number.
        naxial : int
            Number of axial stations (i.e., divisions) of the segment. (min=2).
        nradial : int
            Number of circumferential points of the segment (min=3).
        nose_radius : float
            Nose radius of blunt body.
            NOSERAD is active only if ZONA7U (Hypersonic Aerodynamic Method)
            is used (the METHOD entry of the MKAEROZ Bulk Data equals 2 or –2).
            Furthermore, NOSERAD is used only if the SEGMESH bulk data card is
            the first segment defined in the BODY7 bulk data card.
        iaxis : int
            The index of the axial station where the blunt nose ends.
            IAXIS is active only if ZONA7U (Hypersonic Aerodynamic
            Method) is used.
        itypes : int
            Type of input used to define the circumferential box cuts
            - 1 body of revolution
            - 2 elliptical body
            - 3 arbitrary body
        xs : list[float]
            X-location of the axial station; Xi must be in ascending
            order. (i.e., Xi+1 > Xi)
        cambers : list[float]
            Body camber at the Xi axial station. (Real)
        ys : list[float]
            Body cross-sectional radius if ITYPEi = 1 or the semi-axis length
            of the elliptical body parallel to the Y-axis if ITYPEi=2.
        zs : list[float]
            The semi-axis length of the elliptical body parallel to the Z-axis.
            Used only if ITYPEi=2.
        idys : int
            Identification number of AEFACT bulk data card that specifies
            nradial number of the Y-coordinate locations of the circumferential
            points at the Xi axial station. Use only if ITYPEi=3.
        idzs : int
            Identification number of AEFACT bulk data card that specifies
            nradial number of the Z-coordinate locations of the circumferential
            points at the Xi axial station. Use only if ITYPEi=3.
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        BaseCard.__init__(self)
        self.segmesh_id = segmesh_id
        self.naxial = naxial
        self.nradial = nradial
        self.nose_radius = nose_radius
        self.iaxis = iaxis
        self.itypes = itypes
        self.cambers = cambers
        self.xs = xs
        self.ys = ys
        self.zs = zs
        self.idys = idys
        self.idzs = idzs

        self.idys_ref = None
        self.idzs_ref = None
        self.pid_ref = None

    def validate(self):
        for i, itype in enumerate(self.itypes):
            assert itype in [1, 2, 3], 'itypes[%i]=%s is invalid; itypes=%s' % (i, itype, self.itypes)

        xi_old = self.xs[0]
        for i, xi in enumerate(self.xs[1:]):
            if xi <= xi_old:
                raise RuntimeError('xs=%s must be in ascending order\nx%i=%s x%i=%s (old)\n%s' % (
                    self.xs, i+2, xi, i+1, xi_old, str(self)))

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SEGMESH card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        segmesh_id = integer(card, 1, 'segmesh_id')
        naxial = integer(card, 2, 'naxial')
        nradial = integer(card, 3, 'nradial')
        nose_radius = double_or_blank(card, 4, 'nose_radius')
        iaxis = integer_or_blank(card, 5, 'iaxis')

        itypes = []
        xs = []
        ys = []
        zs = []
        cambers = []
        idys = []
        idzs = []
        assert len(card) >= 9, f'len(SEGMESH card) = {len(card):d}\ncard={card}'

        for counter, ifield in enumerate(range(9, len(card), 8)):
            itype = integer(card, ifield, 'itype%d' % (counter+1))
            x = double_or_blank(card, ifield+1, 'itype%d' % (counter+1), default=0.)
            camber = double_or_blank(card, ifield+2, 'camber%d' % (counter+1), default=0.)
            y = double_or_blank(card, ifield+3, 'y%d' % (counter+1), default=0.)
            z = double_or_blank(card, ifield+4, 'z%d' % (counter+1), default=0.)
            idy = integer_or_blank(card, ifield+5, 'idy%d' % (counter+1))
            idz = integer_or_blank(card, ifield+6, 'idz%d' % (counter+1))
            itypes.append(itype)
            xs.append(x)
            ys.append(y)
            zs.append(z)
            cambers.append(camber)
            idys.append(idy)
            idzs.append(idz)
        assert len(itypes) == naxial, f'naxial={naxial:%d} nradial={nradial:%d} len(itypes)={len(itypes):%d}'
        return SEGMESH(segmesh_id, naxial, nradial, nose_radius, iaxis,
                       itypes, xs, cambers, ys, zs, idys, idzs, comment=comment)

    def Cp(self) -> int:
        # return coord_id(self.cp_ref, self.cp)
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    def Pid(self) -> int:
        if self.pid_ref is not None:
            return self.pid_ref.pid
        return self.pid

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by SEGMESH eid=%d' % self.pid
        idys_ref = []
        idzs_ref = []
        for idy in self.idys:
            idy_ref = None
            if idy is not None and isinstance(idy, integer_types):
                idy_ref = model.AEFact(idy, msg=msg)
                assert len(idy_ref.fractions) > 2, 'idy_ref=%s' % idy_ref
            idys_ref.append(idy_ref)

        for idz in self.idzs:
            idz_ref = None
            if idz is not None and isinstance(idz, integer_types):
                idz_ref = model.AEFact(idz, msg=msg)
                assert len(idz_ref.fractions) > 2, 'idz_ref=%s' % idz_ref
            idzs_ref.append(idz_ref)
        self.idys_ref = idys_ref
        self.idzs_ref = idzs_ref
        #print(self.idys_ref)

    def safe_cross_reference(self, model: BDF, xref_errors):
        return self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        #self.cp = self.Cp()
        #self.idys = idys_ref
        #self.idzs = idzs_ref

        self.pid_ref = None
        #self.cp_ref = None
        self.idys_ref = None
        self.idzs_ref = None

    #def set_points(self, points):
        #self.p1 = np.asarray(points[0])
        #p2 = np.asarray(points[1])
        #x12 = p2 - self.p1
        #self.x12 = x12[0]

    def shift(self, dxyz: np.ndarray) -> None:
        """shifts the aero panel"""
        self.p1 += dxyz

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        list_fields = [
            'SEGMESH', self.segmesh_id, self.naxial, self.nradial, self.nose_radius,
            self.iaxis, None, None, None]
        for itype, x, camber, y, z, idy, idz in zip(self.itypes, self.xs, self.cambers,
                                                    self.ys, self.zs, self.idys, self.idzs):
            list_fields += [itype, x, camber, y, z, idy, idz, None]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO7(BaseCard):
    """
    Totally wrong...

    Defines an aerodynamic macro element (panel) in terms of two leading edge
    locations and side chords. This is used for Doublet-Lattice theory for
    subsonic aerodynamics and the ZONA51 theory for supersonic aerodynamics.

    +--------+-----+-------+--------+-------+--------+--------+--------+---------+
    |   1    |  2  |   3   |   4    |   5   |   6    |    7   |   8    |    9    |
    +========+=====+=======+========+=======+========+========+========+=========+
    | CAERO7 | WID | LABEL | ACOORD | NSPAN | NCHORD |  LSPAN |  ZTAIC | PAFOIL7 |
    +--------+-----+-------+--------+-------+--------+--------+--------+---------+
    |        | XRL |  YRL  |   ZRL  |  RCH  |  LRCHD | ATTCHR | ACORDR |         |
    +--------+-----+-------+--------+-------+--------+--------+--------+---------+
    |        | XTL |  YTL  |   ZTL  |  TCH  |  LTCHD | ATTCHT | ACORDT |         |
    +--------+-----+-------+--------+-------+--------+--------+--------+---------+

    ::

      1
      | \
      |   \
      |     \
      |      4
      |      |
      |      |
      2------3

    Attributes
    ----------
    eid : int
        element id
    p_airfoil : int
        int : PAFOIL7 ID
    p1 : (1, 3) ndarray float
        xyz location of point 1 (leading edge; inboard)
    p4 : (1, 3) ndarray float
        xyz location of point 4 (leading edge; outboard)
    x12 : float
        distance along the flow direction from node 1 to node 2; (typically x, root chord)
    x43 : float
        distance along the flow direction from node 4 to node 3; (typically x, tip chord)
    cp : int, Coord
        int : coordinate system
        Coord : Coordinate object (xref)
    nspan : int
        int > 0 : N spanwise boxes distributed evenly
        int = 0 : use lchord
    nchord : int
        int > 0 : N chordwise boxes distributed evenly
        int = 0 : use lchord
    lspan : int, AEFACT
        int > 0 : AEFACT reference for non-uniform nspan
        int = 0 : use nspan
        AEFACT : AEFACT object  (xref)
    lchord : int, AEFACT
        int > 0 : AEFACT reference for non-uniform nchord
        int = 0 : use nchord
        AEFACT : AEFACT object  (xref)
    comment : str; default=''
         a comment for the card

    """
    type = 'CAERO7'
    _field_map = {
        1: 'sid', 2: 'pid', 3: 'cp', 4: 'nspan', 5: 'nchord',
        6: 'lspan', 7: 'lchord', 8: 'igid', 12: 'x12', 16: 'x43',
    }

    def __init__(self, eid: int, label: str,
                 p1: np.ndarray, x12: float,
                 p4: np.ndarray, x43: float,
                 cp: int=0,
                 nspan: int=0, nchord: int=0,
                 lspan: int=0,
                 p_airfoil: Optional[int]=None,
                 ztaic: Optional[int]=None,
                 lchord_root: int=0, attach_root: int=0, achord_root: int=0,
                 lchord_tip: int=0, attach_tip: int=0, achord_tip: int=0,
                 comment: str=''):
        """
        Defines a CAERO1 card, which defines a simplified lifting surface
        (e.g., wing/tail).

        Parameters
        ----------
        eid : int
            element id
        p_airfoil : int
            int : PAFOIL7 ID
        p1 : (1, 3) ndarray float
            xyz location of point 1 (leading edge; inboard)
        p4 : (1, 3) ndarray float
            xyz location of point 4 (leading edge; outboard)
        x12 : float
            distance along the flow direction from node 1 to node 2; (typically x, root chord)
        x43 : float
            distance along the flow direction from node 4 to node 3; (typically x, tip chord)
        cp : int, Coord; default=0
            int : coordinate system
            Coord : Coordinate object (xref)
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        nchord : int; default=0
            int > 0 : N chordwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
            AEFACT : AEFACT object  (xref)
        lchord : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nchord
            int = 0 : use nchord
            AEFACT : AEFACT object  (xref)
        comment : str; default=''
             a comment for the card

        """
        BaseCard.__init__(self)
        if cp is None:
            cp = 0
        if nspan is None:
            nspan = 0
        if nchord is None:
            nchord = 0
        p1 = np.asarray(p1)
        p4 = np.asarray(p4)

        if comment:
            self.comment = comment
        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO2 entry.
        self.label = label

        #: Coordinate system for locating point 1.
        self.cp = cp
        self.nspan = nspan
        self.nchord = nchord
        self.lspan = lspan
        self.p_airfoil = p_airfoil
        self.p1 = p1
        self.x12 = x12
        self.p4 = p4
        self.x43 = x43
        self.ztaic = ztaic

        self.lchord_root = lchord_root
        self.attach_root = attach_root
        self.achord_root = achord_root

        self.lchord_tip = lchord_tip
        self.attach_tip = attach_tip
        self.achord_tip = achord_tip

        self.pid_ref = None
        self.cp_ref = None
        self.lchord_ref = None
        self.lspan_ref = None
        self.ascid_ref = None
        self.box_ids = None
        self.pafoil_ref = None
        #self._init_ids() #TODO: make this work here?

    def validate(self):
        msg = ''
        is_failed = False
        if self.nspan == 0 and self.lspan == 0:
            msg += 'NSPAN or LSPAN must be greater than 0; nspan=%r nlspan=%s\n' % (
                self.nspan, self.lspan)
            is_failed = True
        if self.nspan <= 0:
            msg += 'NSPAN must be greater than 0; nspan=%r\n' % (
                self.nspan)
            is_failed = True

        if self.nchord <= 0:
            msg += 'NCHORD must be greater than 0; nchord=%r\n' % (
                self.nchord)
            is_failed = True
        if is_failed:
            msg += str(self)
            raise ValueError(msg)
        assert len(self.p1) == 3, 'p1=%s' % self.p1
        assert len(self.p4) == 3, 'p4=%s' % self.p4
        assert self.nspan < 100, 'nspan=%s\n%s' % (self.nspan, str(self))
        assert self.nchord < 100, 'nchord=%s\n%s' % (self.nchord, str(self))

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a CAERO7 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        CAERO7  100101  RIBW1           2       25
                998.904  39.821 230.687 1298.159        310001
                1121.821 61.134 236.560 1175.243
        CAERO7  100201  RIBW2           2       25
                1121.821 61.134 236.560 1175.243
                1244.258 84.704 243.625 1052.805

        """
        eid = integer(card, 1, 'eid')
        name = string(card, 2, 'name')
        cp = integer_or_blank(card, 3, 'cp', default=0)
        nspan = integer_or_blank(card, 4, 'nspan', default=0)
        nchord = integer_or_blank(card, 5, 'nchord', default=0)
        lspan = integer_or_blank(card, 6, 'aefact_lchord', default=0)
        if lspan:
            lspan = 0
        ztaic = integer_or_blank(card, 7, 'ztaic')
        p_airfoil = integer_or_blank(card, 8, 'aefact')
        #assert cp == 0
        #igroup = integer(card, 8, 'igid')

        x1 = double_or_blank(card, 9, 'x1', default=0.0)
        y1 = double_or_blank(card, 10, 'y1', default=0.0)
        z1 = double_or_blank(card, 11, 'z1', default=0.0)
        p1 = np.array([x1, y1, z1])
        x12 = double_or_blank(card, 12, 'x12', default=0.)
        lchord_root = integer_or_blank(card, 13, 'lchord_root')
        attach_root = integer_or_blank(card, 14, 'attach_root')
        achord_root = integer_or_blank(card, 15, 'achord_root')

        x4 = double_or_blank(card, 17, 'x4', default=0.0)
        y4 = double_or_blank(card, 18, 'y4', default=0.0)
        z4 = double_or_blank(card, 19, 'z4', default=0.0)
        p4 = np.array([x4, y4, z4])
        x43 = double_or_blank(card, 20, 'x43', default=0.)
        lchord_tip = integer_or_blank(card, 21, 'lchord_tip')
        attach_tip = integer_or_blank(card, 22, 'attach_tip')
        achord_tip = integer_or_blank(card, 23, 'achord_tip')

        # print(card.write_card())
        assert len(card) <= 24, f'len(CAERO7 card) = {len(card):d}\ncard={card}'
        return CAERO7(eid, name, p1, x12, p4, x43,
                      cp=cp, nspan=nspan, nchord=nchord, lspan=lspan,
                      p_airfoil=p_airfoil, ztaic=ztaic,
                      lchord_root=lchord_root, attach_root=attach_root, achord_root=achord_root,
                      lchord_tip=lchord_tip, attach_tip=attach_tip, achord_tip=achord_tip,
                      comment=comment)

    def flip_normal(self):
        self.p1, self.p4 = self.p4, self.p1
        self.x12, self.x43 = self.x43, self.x12

    def _init_ids(self, dtype='int32'):
        """
        Fill `self.box_ids` with the sub-box ids. Shape is (nchord, nspan)

        """
        nchord, nspan = self.shape
        assert nchord >= 1, 'nchord=%s' % nchord
        assert nspan >= 1, 'nspan=%s' % nspan
        self.box_ids = np.zeros((nchord, nspan), dtype=dtype)

        npanels = nchord * nspan
        try:
            self.box_ids = np.arange(self.eid, self.eid + npanels,
                                     dtype=dtype).reshape(nspan, nchord) # .T
        except OverflowError:
            if dtype == 'int64':
                # we already tried int64
                msg = 'eid=%s nspan=%s nchord=%s' % (
                    self.eid, nspan, nchord)
                raise OverflowError(msg)
            self._init_ids(dtype='int64')

    def Cp(self):
        if self.cp_ref is not None:
            return self.cp_ref.cid
        return self.cp

    #def Pid(self):
        #if self.pid_ref is not None:
            #return self.pid_ref.pid
        #return self.pid

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CAERO7 eid=%s' % self.eid
        #self.pid_ref = model.PAero(self.pid, msg=msg)
        self.cp_ref = model.Coord(self.cp, msg=msg)
        self.ascid_ref = model.Acsid(msg=msg)

        #if self.nchord == 0:
            #assert isinstance(self.lchord, integer_types), self.lchord
            #self.lchord_ref = model.AEFact(self.lchord, msg)
        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan_ref = model.AEFact(self.lspan, msg)

        if self.p_airfoil:
            self.pafoil_ref = model.zaero.PAFOIL(self.p_airfoil, msg)
        self._init_ids()

    def safe_cross_reference(self, model: BDF, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CAERO1 eid=%s' % self.eid
        #try:
            #self.pid_ref = model.PAero(self.pid, msg=msg)
        #except KeyError:
            #pass

        self.cp_ref = model.safe_coord(self.cp, self.eid, xref_errors, msg=msg)
        self.ascid_ref = model.safe_acsid(msg=msg)

        #if self.nchord == 0:
            #assert isinstance(self.lchord, integer_types), self.lchord
            #self.lchord_ref = model.safe_aefact(self.lchord, self.eid, xref_errors, msg)

        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan_ref = model.safe_aefact(self.lspan, self.eid, xref_errors, msg)

        self._init_ids()

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        #self.pid = self.Pid()
        self.cp = self.Cp()
        self.lspan = self.get_LSpan()
        #self.pid_ref = None
        self.cp_ref = None
        self.lspan_ref = None
        self.ascid_ref = None

    def convert_to_nastran(self, pid: int=1):
        """
        +--------+-----+-----+----+-------+--------+--------+--------+------+
        |   1    |  2  |  3  | 4  |   5   |   6    |    7   |   8    |   9  |
        +========+=====+=====+====+=======+========+========+========+======+
        | CAERO1 | EID | PID | CP | NSPAN | NCHORD |  LSPAN | LCHORD | IGID |
        +--------+-----+-----+----+-------+--------+--------+--------+------+
        |        |  X1 | Y1  | Z1 |  X12  |   X4   |   Y4   |   Z4   | X43  |
        +--------+-----+-----+----+-------+--------+--------+--------+------+

        """
        igroup = 1
        caero = CAERO1(self.eid, pid, igroup, self.p1, self.x12,
                       self.p4, self.x43, cp=self.cp,
                       nspan=self.nspan, lspan=self.lspan,
                       nchord=self.nchord, lchord=0,
                       comment=self.comment)
        caero.validate()
        return caero

    @property
    def min_max_eid(self) -> list[int]:
        """
        Gets the min and max element ids of the CAERO card

        Returns
        -------
        min_max_eid : (2, ) list
            The [min_eid, max_eid]

        """
        nchord, nspan = self.shape
        return [self.eid, self.eid + nchord * nspan]

    def get_points(self):
        """
        Get the 4 corner points for the CAERO card

        Returns
        -------
        p1234 : (4, 3) list
             List of 4 corner points in the global frame

        """
        if self.cp_ref is None and self.cp == 0:
            p1 = self.p1
            p4 = self.p4
        else:
            p1 = self.cp_ref.transform_node_to_global(self.p1)
            p4 = self.cp_ref.transform_node_to_global(self.p4)

        if self.ascid_ref is None:
            # yes, this really does list + array addition
            p2 = p1 + np.array([self.x12, 0., 0.])
            p3 = p4 + np.array([self.x43, 0., 0.])
        else:
            p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))
            p3 = p4 + self.ascid_ref.transform_vector_to_global(np.array([self.x43, 0., 0.]))
        return [p1, p2, p3, p4]

    def get_box_index(self, box_id):
        """
        Get the index of ``self.box_ids`` that corresponds to the given box id.

        Parameters
        -----------
        box_id : int
            Box id to get the index of.

        Returns
        --------
        index : tuple
            Index of ``self.box_ids`` that corresponds to the given box id.

        """
        if box_id not in self.box_ids:
            self._box_id_error(box_id)
        index = np.where(self.box_ids == box_id)
        index = (index[0][0], index[1][0])
        return index

    def get_box_quarter_chord_center(self, box_id):
        """
        The the location of the quarter chord of the box along the centerline.

        Parameters
        -----------
        box_id : int
            Box id.

        Returns
        --------
        xyz_quarter_chord : ndarray
            Location of box quarter chord in global.

        """
        return self._get_box_x_chord_center(box_id, 0.25)

    def get_box_mid_chord_center(self, box_id):
        """
        The the location of the mid chord of the box along the centerline.

        Parameters
        -----------
        box_id : int
            Box id.

        Returns
        --------
        xyz_mid_chord : ndarray
            Location of box mid chord in global.

        """
        return self._get_box_x_chord_center(box_id, 0.5)

    def _get_box_x_chord_center(self, box_id, x_chord):
        """The the location of the x_chord of the box along the centerline."""
        raise NotImplementedError('CAERO7: _get_box_x_chord_center')
        #if self.lchord != 0 or self.lspan != 0:
            #raise NotImplementedError()
        #ichord, ispan = self.get_box_index(box_id)

        #le_vector = self.p4 - self.p1
        #delta_xyz = le_vector * ((ispan + 0.5)/self.nspan)
        #yz = delta_xyz[1:3] + self.p1[1:3]
        #chord = ((ispan + 0.5)/self.nspan) * (self.x43 - self.x12) + self.x12
        #x = (ichord + x_chord)/self.nchord * chord + self.p1[0] + delta_xyz[0]
        #return np.array([x, yz[0], yz[1]])

    def _box_id_error(self, box_id):
        """
        Raise box_id IndexError.
        """
        msg = '%i not in range of aero box ids\nRange: %i to %i' % (box_id, self.box_ids[0, 0],
                                                                    self.box_ids[-1, -1])
        raise IndexError(msg)

    @property
    def npanels(self):
        nchord, nspan = self.shape
        return nchord * nspan

    @property
    def shape(self):
        """returns (nelements_nchord, nelements_span)"""
        if self.nchord == 0:
            x = self.lchord_ref.fractions
            nchord = len(x) - 1
        else:
            nchord = self.nchord

        if self.nspan == 0:
            y = self.lspan_ref.fractions
            nspan = len(y) - 1
        else:
            nspan = self.nspan
        if nchord < 1 or nspan < 1:
            msg = 'CAERO7 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                self.eid, self.nchord, self.nspan, self.lchord, self.lspan)
            raise RuntimeError(msg)
        return nchord, nspan

    def get_panel_npoints_nelements(self):
        """
        Gets the number of sub-points and sub-elements for the CAERO card

        Returns
        -------
        npoints : int
            The number of nodes for the CAERO
        nelmements : int
            The number of elements for the CAERO

        """
        nchord, nspan = self.shape
        nelements = nchord * nspan
        npoints = (nchord + 1) * (nspan + 1)
        return npoints, nelements

    @property
    def xy(self):
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        if self.nchord == 0:
            x = self.lchord_ref.fractions
            nchord = len(x) - 1
        else:
            nchord = self.nchord
            x = np.linspace(0., 1., nchord + 1)

        if self.nspan == 0:
            y = self.lspan_ref.fractions
            nspan = len(y) - 1
        else:
            nspan = self.nspan
            y = np.linspace(0., 1., nspan + 1)

        if nchord < 1 or nspan < 1:
            msg = 'CAERO7 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                self.eid, self.nchord, self.nspan, self.lchord, self.lspan)
            raise RuntimeError(msg)
        return x, y

    def panel_points_elements(self):
        """
        Gets the sub-points and sub-elements for the CAERO card

        Returns
        -------
        points : (nnodes,3) ndarray of floats
            the array of points
        elements : (nelements,4) ndarray of integers
            the array of point ids

        """
        p1, p2, p3, p4 = self.get_points()
        x, y = self.xy
        # We're reordering the points so we get the node ids and element ids
        # to be consistent with Nastran.  This is only useful if you're plotting
        # aero panel forces
        #
        # this gives us chordwise panels and chordwise nodes
        return points_elements_from_quad_points(p1, p4, p3, p2, y, x, dtype='int32')

    def set_points(self, points):
        self.p1 = points[0]
        p2 = points[1]
        p3 = points[2]
        self.p4 = points[3]
        self.x12 = p2[0] - self.p1[0]
        self.x43 = p3[0] - self.p4[0]
        assert self.x12 >= 0., 'p1=%s p2=%s' % (self.p1, p2)
        assert self.x43 >= 0., 'p4=%s p3=%s' % (self.p4, p3)
        assert self.x12 > 0. or self.x43 > 0., 'points=%s' % str(points)

    def shift(self, dxyz):
        """shifts the aero panel"""
        self.p1 += dxyz
        self.p4 += dxyz

    def plot(self, ax: AxesSubplot) -> None:
        """plots the panels"""
        points, elements = self.panel_points_elements()
        for eid, elem in enumerate(elements[:, [0, 1, 2, 3, 0]]):
            pointsi = points[elem][:, [0, 1]]
            x = pointsi[:, 0]
            y = pointsi[:, 1]
            ax.plot(x, y, color='b')
            box_id = self.eid + eid
            centroid = (x[:-1].sum() / 4, y[:-1].sum() / 4)
            elem_name = f'e{box_id}'
            ax.annotate(elem_name, centroid, ha='center')

            for pid, point in zip(elem, pointsi):
                point_name = f'p{pid}'
                ax.annotate(point_name, point, ha='center')

    def get_LSpan(self) -> int:
        if self.lspan_ref is not None:
            return self.lspan_ref.sid
        return self.lspan

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
          the fields that define the card

        """
        cp = self.Cp()
        nspan = self.nspan
        nchord = self.nchord
        lspan = self.get_LSpan()
        list_fields = (
            ['CAERO7', self.eid, self.label, cp, nspan, nchord, lspan, self.ztaic,
             self.p_airfoil,] +
            list(self.p1) + [self.x12, self.lchord_root, self.attach_root, self.achord_root, None] +
            list(self.p4) + [self.x43, self.lchord_tip, self.attach_tip, self.achord_tip, None])
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : LIST
            The fields that define the card

        """
        cp = set_blank_if_default(self.Cp(), 0)
        nspan = set_blank_if_default(self.nspan, 0)
        nchord = set_blank_if_default(self.nchord, 0)
        #lchord = set_blank_if_default(self.get_LChord(), 0)
        lspan = set_blank_if_default(self.get_LSpan(), 0)
        #lchord = 0
        lspan = 0
        list_fields = (
            ['CAERO7', self.eid, self.label, cp, nspan, nchord, lspan, self.ztaic,
             self.p_airfoil,] +
            list(self.p1) + [self.x12, None, None, None, None] +
            list(self.p4) + [self.x43, None, None, None, None])
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PAFOIL8(BaseCard):
    """

    +---------+----+------+------+-------+------+------+-------+------+
    |   1     |  2 |   3  |  4   |   5   |   6  |   7  |   8   |  9   |
    +=========+====+======+======+=======+======+======+=======+======+
    | PAFOIL7 | ID | ITAX | ITHR | ICAMR | RADR | ITHT | ICAMT | RADT |
    +---------+----+------+------+-------+------+------+-------+------+
    | PAFOIL7 |  1 | -201 |  202 |  203  |  0.1 |  211 |  212  |  0.1 |
    +---------+----+------+------+-------+------+------+-------+------+

    """
    type = 'PAFOIL8'

    def __init__(self, pid: int, i_axial: int,
                 i_thickness_root: int, i_camber_root: int, le_radius_root: float,
                 i_thickness_tip: int, comment: str=''):
        """
        Defines a BODY7 card, which defines a slender body
        (e.g., fuselage/wingtip tank).

        Parameters
        ----------
        pid : int
            PAFOIL8 identification number.
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        self.pid = pid
        self.i_axial = i_axial

        self.i_thickness_root = i_thickness_root
        self.i_camber_root = i_camber_root
        self.le_radius_root = le_radius_root
        self.i_thickness_tip = i_thickness_tip

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PAFOIL8 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """

        # card = ['PAFOIL8', '4001',
        #         '1.0', '4002',
        #         '1.0', '4002', '0']
        pid = integer(card, 1, 'pid')

        i_axial = double(card, 2, 'i_axial')
        i_thickness_root = integer(card, 3, 'i_thickness_root')
        i_camber_root = double(card, 4, 'i_camber_root')
        le_radius_root = integer(card, 5, 'le_radius_root')
        i_thickness_tip = integer(card, 6, 'i_thickness_tip')

        assert len(card) == 7, f'len(PAFOIL8 card) = {len(card):d}\ncard={card}'
        return PAFOIL8(pid, i_axial,
                       i_thickness_root, i_camber_root, le_radius_root,
                       i_thickness_tip, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PAFOIL8 pid=%s' % self.pid

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def convert_to_nastran(self, model):
        raise NotImplementedError('PAFOIL8: convert_to_nastran')

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        list_fields = [
            'PAFOIL8', self.pid, self.i_axial,
            self.i_thickness_root, self.i_camber_root, self.le_radius_root,
            self.i_thickness_tip,
        ]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AESURFZ(BaseCard):
    """
    Specifies an aerodynamic control surface for aeroservoelastic, static
    aeroelastic/trim analysis, or the transient response analysis.

    +---------+--------+-------+-------+-------+--------+--------+--------+--------+
    |    1    |   2    |   3   |   4   |   5   |   6    |    7   |   8    |   9    |
    +=========+========+=======+=======+=======+========+========+========+========+
    | AESURFZ | LABEL  |  TYPE |  CID  |  SETK |  SETG  |  ACTID |        |        |
    +---------+--------+-------+-------+-------+--------+--------+--------+--------+
    | AESURFZ | RUDDER |  ASYM |   1   |   10  |   20   |    0   |        |        |
    +---------+--------+-------+-------+-------+--------+--------+--------+--------+
    """
    type = 'AESURFZ'

    @property
    def aesid(self) -> str:
        return self.label

    @property
    def alid1_ref(self):
        return None

    def __init__(self, label: str, surface_type: str, cid: int,
                 panlst: int, setg: int, actuator_tf: int,
                 comment: str=''):
        """
        Creates an AESURF card, which defines a control surface

        Parameters
        ----------
        label : str
            controller name
        surface_type : str
            defines the control surface type {SYM, ASYM}
        cid : int
            coordinate system id to define the hinge axis
        panlst : int
            aero panels defined by PANLST
        setg : int
           ???
        actuator_tf : int
            ???
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        #: Controller name.
        self.label = label
        self.surface_type = surface_type
        self.panlst = panlst
        self.setg = setg
        self.actuator_tf = actuator_tf

        #: Identification number of a rectangular coordinate system with a
        #: y-axis that defines the hinge line of the control surface
        #: component.
        self.cid = cid
        assert surface_type in ['SYM', 'ANTI', 'ANTISYM', 'ASYM'], f'surface_type={surface_type!r}'

        self.cid_ref = None
        self.panlst_ref = None
        self.aero_element_ids = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds an AESURF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        label = string(card, 1, 'label')
        surface_type = string(card, 2, 'TYPE')
        cid = integer(card, 3, 'CID')
        panlst = integer(card, 4, 'PANLST/SETK') # PANLST1, PANLST2, PANLST3
        setg = integer_or_blank(card, 5, 'SETG', default=0) # SET1, SETADD
        actuator_tf = integer_or_blank(card, 6, 'ACTID') # ACTU card
        assert len(card) <= 7, f'len(AESURFZ card) = {len(card):d}\ncard={card}'
        return AESURFZ(label, surface_type, cid, panlst, setg, actuator_tf, comment=comment)

    def Cid(self) -> int:
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def SetK(self):
        for panlst in self.panlst_ref:
            if panlst is not None:
                return panlst.eid
        return self.panlst

    #def aelist_id1(self):
        #if self.alid1_ref is not None:
            #return self.alid1_ref.sid
        #return self.alid1

    #def aelist_id2(self):
        #if self.alid2_ref is not None:
            #return self.alid2_ref.sid
        #return self.alid2

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        self.cid_ref = model.Coord(self.cid)
        self.panlst_ref, self.aero_element_ids = cross_reference_panlst(
            model, self.panlst)

    def safe_cross_reference(self, model: BDF, xref_errors):
        msg = ', which is required by AESURFZ aesid=%s' % self.aesid
        self.cid_ref = model.safe_coord(self.cid, self.aesid, xref_errors, msg=msg)
        self.panlst_ref, self.aero_element_ids = cross_reference_panlst(
            model, self.panlst)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.cid = self.Cid()
        self.cid_ref = None

        self.panlst = self.SetK()
        self.panlst_ref = None

    def convert_to_nastran(self, model, aesurf_id, aelist_id):
        """
        +--------+--------+-------+-------+-------+--------+--------+--------+--------+
        |    1   |   2    |   3   |   4   |   5   |   6    |    7   |   8    |   9    |
        +========+========+=======+=======+=======+========+========+========+========+
        | AESURF |   ID   | LABEL | CID1  | ALID1 |  CID2  | ALID2  |  EFF   |  LDW   |
        +--------+--------+-------+-------+-------+--------+--------+--------+--------+
        |        |  CREFC | CREFS | PLLIM | PULIM | HMLLIM | HMULIM | TQLLIM | TQULIM |
        +--------+--------+-------+-------+-------+--------+--------+--------+--------+

        +---------+--------+-------+-------+-------+--------+--------+
        |    1    |   2    |   3   |   4   |   5   |   6    |    7   |
        +=========+========+=======+=======+=======+========+========+
        | AESURFZ | LABEL  |  TYPE |  CID  |  SETK |  SETG  |  ACTID |
        +---------+--------+-------+-------+-------+--------+--------+
        | AESURFZ | RUDDER |  ASYM |   1   |   10  |   20   |    0   |
        +---------+--------+-------+-------+-------+--------+--------+
        """
        assert self.surface_type == 'ASYM', str(self)
        aelist = AELIST(aelist_id, self.aero_element_ids)
        aesurf = AESURF(aesurf_id, self.label, self.cid, aelist_id,
                        cid2=None, aelist_id2=None,
                        eff=1.0, ldw='LDW', crefc=1.0,
                        crefs=1.0, pllim=-np.pi/2.,
                        pulim=np.pi/2., hmllim=None,
                        hmulim=None, tqllim=None,
                        tqulim=None, comment=self.comment)
        aesurf.validate()
        aelist.validate()
        return aelist, aesurf

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fieldsreset_camera[int/float/str]
            the fields that define the card

        """
        list_fields = ['AESURFZ', self.label, self.surface_type, self.cid,
                       self.panlst, self.setg, self.actuator_tf]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[int/float/str]
            the fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        Writes the card with the specified width and precision

        Parameters
        ----------
        size : int (default=8)
            size of the field; {8, 16}
        is_double : bool (default=False)
            is this card double precision

        Returns
        -------
        msg : str
            the string representation of the card

        """
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AESLINK(BaseCard):
    r"""
    Defines relationships between or among AESTAT and AESURF entries, such
    that:

    .. math:: u^D + \Sigma_{i=1}^n C_i u_i^I = 0.0

    +---------+-------+-------+--------+------+-------+----+-------+----+
    |    1    |   2   |   3   |   4    |   5  |   6   |  7 |   8   |  9 |
    +=========+=======+=======+========+======+=======+====+=======+====+
    | AESLINK |  ID   | LABLD | LABL1  |  C1  | LABL2 | C2 | LABL3 | C3 |
    +---------+-------+-------+--------+------+-------+----+-------+----+
    |         | LABL4 |  C4   |  etc.  |      |       |    |       |    |
    +---------+-------+-------+--------+------+-------+----+-------+----+
    | AESLINK |  10   | INBDA |  OTBDA | -2.0 |       |    |       |    |
    +---------+-------+-------+--------+------+-------+----+-------+----+
    """
    type = 'AESLINK'

    # @classmethod
    # def _init_from_empty(cls):
    #     aelink_id = 1
    #     label = 'ELEV'
    #     independent_labels = ['ELEV1', 'ELEV2']
    #     linking_coefficients = [1., 2.]
    #     return AELINK(aelink_id, label, independent_labels, linking_coefficients, comment='')

    def __init__(self, label: int | str,
                 link_type: str,
                 actu_id: int, independent_labels: list[str],
                 linking_coefficients: list[float],
                 comment: str='') -> None:
        """
        Creates an AESLINK card, which defines an equation linking
        ??? and AESURFZ cards

        Parameters
        ----------
        label : int/str
            unique id
        actu_id : id
            ACTU card
        independent_labels : list[str, ..., str]
            name for the independent variables (AESTATs)
        linking_coefficients : list[float]
            linking coefficients
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.label = label
        self.link_type = link_type
        self.actu_id = actu_id

        #: defines the independent variable name (string)
        self.independent_labels = [independent_label.upper()
                                   for independent_label in independent_labels]

        #: linking coefficients (real)
        self.linking_coefficients = linking_coefficients

        self.dependent_label_ref = None
        self.independent_labels_ref = []

    @property
    def dependent_label(self) -> str:
        return self.label

    # def validate(self):
    #     errors = []
    #     if not self.label[0].isalpha():
    #         msgi = f'label={self.label!r} must start with a character'
    #         errors.append(msgi)
    #
    #     for label in self.independent_labels:
    #         if not label[0].isalpha():
    #             msgi = f' ind_label={label!r} must start with a character'
    #             errors.append(msgi)
    #
    #     if isinstance(self.aelink_id, integer_types) and self.aelink_id < 0:
    #         msgi = "aelink_id must be greater than or equal to 0 (or 'ALWAYS')"
    #         errors.append(msgi)
    #
    #     if len(self.independent_labels) != len(self.linking_coefficients):
    #         msgi = 'nlabels=%d nci=%d\nindependent_labels=%s linking_coefficients=%s\n%s' % (
    #             len(self.independent_labels), len(self.linking_coefficients),
    #             self.independent_labels, self.linking_coefficients, str(self))
    #         errors.append(msgi)
    #
    #     if len(self.independent_labels) == 0:
    #         msgi = 'nlabels=%d nci=%d\nindependent_labels=%s linking_coefficients=%s\n%s' % (
    #             len(self.independent_labels), len(self.linking_coefficients),
    #             self.independent_labels, self.linking_coefficients, str(self))
    #         errors.append(msgi)
    #
    #     if errors:
    #         msg = f'AESLINK id={self.aelink_id:d}\n -' + '\n - '.join(errors)
    #         raise RuntimeError(msg.rstrip('\n- '))

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds an AESLINK card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # AESLINK LABEL  TYPE    ACTID
        #         COEFF1 AESURF1 COEFF2 AESURF2 -etc-
        # AESLINK AES1 SYM  100
        #         1.0  AES2 0.5 AES3 0.3 AES4

        label = integer_or_string(card, 1, 'label')
        link_type = string(card, 2, 'link_type')
        actu_id = integer(card, 3, 'actu_id')
        assert link_type in {'SYM', 'ANTISYM', 'ASYM'}, f'link_type={link_type!r}'
        independent_labels = []
        linking_coefficients = []

        nfields = len(card)

        assert nfields % 2 == 1, f'nfields={nfields} must be odd'
        j = 1
        for ifield in range(9, nfields, 2):
            linking_coefficent = double(card, ifield, f'coeff{j}')
            independent_label = string(card, ifield+1, f'label{j}')
            independent_labels.append(independent_label)
            linking_coefficients.append(linking_coefficent)
        return AESLINK(label, link_type, actu_id, independent_labels, linking_coefficients,
                       comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        We're simply going to validate the labels
        Updating the labels on the AESURFs will NOT propogate to the AELINK
        """
        msg = ', which is required by:\n%s' % str(self)
        zaero = model.zaero
        self.actu_ref = zaero.actu[self.actu_id]


        aesurf_names = {aesurf.label for aesurf in model.aesurf.values()}
        #print(f'dependent_label={self.dependent_label}; aesurf_names={aesurf_names}')
        aestat_names = {aestat.label for aestat in model.aestats.values()}
        aparam_names = []
        is_aesurf = self.dependent_label in aesurf_names
        is_aeparam = False
        # if not (is_aesurf or is_aeparam):
        #     raise RuntimeError(f'dependent_label={self.dependent_label} is an AESURF and AEPARM\n{self}\n'
        #                        f'aesurf={list(model.aesurf.keys())} aeparam={list(model.aeparams.keys())}')

        if 0:
            if is_aesurf:
                self.dependent_label_ref = model.AESurf(self.dependent_label,
                                                        msg=f'dependent_label={self.dependent_label!r}; ' + msg)
            elif is_aeparam:
                self.dependent_label_ref = model.AEParam(self.dependent_label, msg=msg)
            else:
                asdf
                return
                raise RuntimeError(f'dependent_label={self.dependent_label} is not an AESURF or AEPARM\n{self}')

        self.independent_labels_ref = []
        for independent_label in self.independent_labels:
            is_aesurf = independent_label in aesurf_names
            is_aeparam = independent_label in aparam_names
            is_aestat = independent_label in aestat_names
            if not (is_aesurf or is_aeparam or is_aestat):
                raise RuntimeError(f'independent_label={independent_label} is an AESURF and AEPARM\n{self}\n'
                                   f'aesurf={list(model.aesurf.keys())} aeparam={list(model.aeparams.keys())}')
            elif is_aesurf:
                independent_aelink = model.AESurf(independent_label, msg=msg)
            #elif is_aeparam:
                #independent_aelink = model.AEParam(independent_label, msg=msg)
            elif is_aestat:
                independent_aelink = model.AEStat(independent_label, msg=msg)
            else:
                raise RuntimeError(f'independent_label={independent_label} is an AESURF and AEPARM\n{self}\n'
                                   f'aesurf={list(model.aesurf.keys())} aeparam={list(model.aeparams.keys())}')
            self.independent_labels_ref.append(independent_aelink)

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    #def uncross_reference(self) -> None:
        #"""Removes cross-reference links"""
        #pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        list_fields : list[int/float/str]
            the fields that define the card

        """
        list_fields = ['AESLINK', self.label, self.link_type, self.actu_id,
                       None, None, None, None, None]
        for (ivar, ival) in zip(self.independent_labels, self.linking_coefficients):
            list_fields += [ival, ivar]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        The writer method used by BDF.write_card()

        Parameters
        -----------
        size : int; default=8
            the size of the card (8/16)

        """
        card = self.raw_fields()
        return self.comment + print_card_8(card)


def cross_reference_panlst(model: BDF,
                           panlist_id: int) -> tuple[list[PANLST1 | PANLST2 | PANLST3],
                                                     np.ndarray]:
    panlst_ref = model.zaero.panlsts[panlist_id]
    aero_ids_list = []
    for panlst in panlst_ref:
        panlst.cross_reference(model)
        aero_ids_list.append(panlst.aero_element_ids)
    aero_element_ids = np.hstack(aero_ids_list)
    return panlst_ref, aero_element_ids
