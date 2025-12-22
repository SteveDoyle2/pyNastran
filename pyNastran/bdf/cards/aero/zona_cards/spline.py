from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.cards.aero.zona_cards.geometry import (
    PANLST1, PANLST2, PANLST3, cross_reference_panlst)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank,
    string_or_blank, blank,
    integer_string_or_blank,
)
from pyNastran.bdf.bdf_interface.assign_type_force import (
    force_double_or_blank,
)
from pyNastran.bdf.cards.aero.aero import Spline, SPLINE1
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class SPLINE1_ZAERO(Spline):
    """
    Defines an infinite plate spline method for displacements and loads
    transferal between CAERO7 macroelement and structural grid points.

    +---------+-------+-------+------+------+------+----+------+-------+
    |    1    |   2   |    3  |   4  |   5  |   6  |  7 |   8  |   9   |
    +=========+=======+=======+======+======+======+====+======+=======+
    | SPLINE1 | EID   | CAERO | BOX1 | BOX2 | SETG | DZ | METH | USAGE |
    +---------+-------+-------+------+------+------+----+------+-------+
    |         | NELEM | MELEM |      |      |      |    |      |       |
    +---------+-------+-------+------+------+------+----+------+-------+
    | SPLINE1 |   3   |  111  | 115  | 122  |  14  | 0. |      |       |
    +---------+-------+-------+------+------+------+----+------+-------+

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE1 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE1 | 100  |       |       |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+

    """
    type = 'SPLINE1_ZAERO'

    def __init__(self, eid: int, panlst: int, setg: int, model: str='', cp=None,
                 dz=None, eps=0.01, comment=''):
        """
        Creates a SPLINE1 card, which is useful for control surface
        constraints.

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
        self.model = model
        self.cp = cp
        self.panlst = panlst
        self.setg = setg
        self.dz = dz
        self.eps = eps
        self.panlst_ref = None
        self.setg_ref = None
        self.aero_element_ids = []

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SPLINE1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        model = string_or_blank(card, 2, 'model', default='')
        cp = blank(card, 3, 'cp')

        panlst = integer(card, 4, 'panlst/setk')
        setg = integer(card, 5, 'setg')
        dz = double_or_blank(card, 6, 'dz', default=0.0)
        eps = double_or_blank(card, 6, 'eps', default=0.01)
        return SPLINE1_ZAERO(eid, panlst, setg, model=model, cp=cp, dz=dz, eps=eps,
                             comment=comment)

    @classmethod
    def add_card_lax(cls, card, comment=''):
        """
        Adds a SPLINE1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        model = string_or_blank(card, 2, 'model')
        cp = blank(card, 3, 'cp')

        panlst = integer(card, 4, 'panlst/setk')
        setg = integer(card, 5, 'setg')
        dz = force_double_or_blank(card, 6, 'dz', default=0.0)
        eps = force_double_or_blank(card, 6, 'eps', default=0.01)
        return SPLINE1_ZAERO(eid, panlst, setg, model=model, cp=cp, dz=dz, eps=eps,
                             comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = f'SPLINE1 eid={self.eid}: setg is missing'
        self.setg_ref = cross_reference_set(model, self.setg, msg=msg)
        if self.setg_ref:
            msg = f', which is required by SPLINE1 eid={self.eid}'
            self.setg_ref.cross_reference_set(model, 'Node', msg=msg)

        self.panlst_ref, self.aero_element_ids = cross_reference_panlst(
            model, self.panlst)
        # model.log.info(f'SPLINE1={self.eid} model={self.model}: boxs={self.aero_element_ids}')

    def safe_cross_reference(self, model: BDF, xref_errors):
        msg = ', which is required by SPLINE1 eid=%s' % self.eid
        self.setg_ref = cross_reference_set(model, self.setg, msg=msg)

        try:
            self.setg_ref.safe_cross_reference(model, 'Node', msg=msg)
        except KeyError:
            model.log.warning('failed to find SETx set_id=%s%s; allowed_sets=%s' % (
                self.setg, msg, np.unique(list(model.sets.keys()))))

        try:
            self.panlst_ref, self.aero_element_ids = cross_reference_panlst(
                model, self.panlst)
            # model.log.info(f'SPLINE1={self.eid} model={self.model}: boxs={self.aero_element_ids}')
        except KeyError:
            pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.panlst_ref = None
        self.setg_ref = None

    def convert_to_nastran(self, model):
        """
        +---------+-------+-------+------+------+------+----+------+-------+
        |    1    |   2   |    3  |   4  |   5  |   6  |  7 |   8  |   9   |
        +=========+=======+=======+======+======+======+====+======+=======+
        | SPLINE1 | EID   | CAERO | BOX1 | BOX2 | SETG | DZ | METH | USAGE |
        +---------+-------+-------+------+------+------+----+------+-------+
        |         | NELEM | MELEM |      |      |      |    |      |       |
        +---------+-------+-------+------+------+------+----+------+-------+
        | SPLINE1 |   3   |  111  | 115  | 122  |  14  | 0. |      |       |
        +---------+-------+-------+------+------+------+----+------+-------+

        +---------+------+-------+-------+------+------+----+-----+-------+
        |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
        +=========+======+=======+=======+======+======+====+=====+=======+
        | SPLINE1 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
        +---------+------+-------+-------+------+------+----+-----+-------+
        | SPLINE1 | 100  |       |       |  1   |  10  | 0. |     |       |
        +---------+------+-------+-------+------+------+----+-----+-------+

        """
        #panlst = 100 # set_aero
        #return SPLINE1(self.eid, panlst, self.setg, model=None, cp=self.cp, dz=self.dz,
                       #eps=0.01, comment=self.comment)
        splines = []
        if not hasattr(self, '_comment'):
            self._comment = ''
        comment = '-' * 72 + '\n' #+ self._comment
        #self._comment = ''
        comment += str(self)
        for panlst in self.panlst_ref:
            for panel_groups in panlst.panel_groups:
                eid = model.zaero.caero_to_name_map[panel_groups]
                caero = model.caeros[eid]
                caero_id = eid
                box1 = caero.eid
                box2 = box1 + caero.npanels - 1
                assert caero.npanels > 0, caero
                #assert box1 > 0 and box2 > 0, 'box1=%s box2=%s' % (box1, box2)
                spline = SPLINE1(eid, caero_id, box1, box2, self.setg, dz=self.dz,
                                 method='IPS', usage='BOTH',
                                 nelements=10, melements=10, comment=comment)
                spline.validate()
                splines.append(spline)
                comment = ''
        return splines

    def raw_fields(self):
        list_fields = ['SPLINE1', self.eid, self.model, self.cp, self.panlst, self.setg,
                       self.dz, self.eps]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SPLINE2_ZAERO(Spline):
    """
    Defines an infinite plate spline method for displacements and loads
    transferal between CAERO7 macroelement and structural grid points.

    +---------+------+-------+------+------+----+-----+-------+-------+
    |    1    |  2   |   3   |  5   |   6  |  6 |  7  |   8   |   9   |
    +=========+======+=======+======+======+====+=====+=======+=======+
    | SPLINE2 | EID  | MODEL | SETK | SETG | DZ | EPS |  CP   | CURV  |
    +---------+------+-------+------+------+----+-----+-------+-------+
    | SPLINE2 | 100  |       |  1   |  10  | 0. |     |       |       |
    +---------+------+-------+------+------+----+-----+-------+-------+

    """
    type = 'SPLINE2_ZAERO'

    def __init__(self, eid: int, panlst: int,
                 setg: int, model: str='', dz=None,
                 eps: float=0.01, cp=None, curvature=None,
                 comment: str=''):
        """
        Creates a SPLINE2 card, which is useful for control surface
        constraints.

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
        self.model = model
        self.cp = cp
        self.panlst = panlst
        self.setg = setg
        self.dz = dz
        self.eps = eps
        self.curvature = curvature
        self.panlst_ref = None
        self.setg_ref = None
        self.aero_element_ids = []

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SPLINE2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        model = str(integer_string_or_blank(card, 2, 'model', default=''))
        panlst = integer(card, 3, 'panlst/setk')
        setg = integer(card, 4, 'setg')
        dz = double_or_blank(card, 5, 'dz', default=0.0)
        eps = double_or_blank(card, 6, 'eps', default=0.01)
        cp = integer_or_blank(card, 7, 'cp', default=0)
        curvature = double_or_blank(card, 8, 'curvature', default=1.0)
        return SPLINE2_ZAERO(eid, panlst, setg, model=model, cp=cp, dz=dz, eps=eps,
                             curvature=curvature, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by SPLINE1 eid=%s' % self.eid
        self.setg_ref = cross_reference_set(model, self.setg, msg=msg)
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        #self.caero_ref = model.CAero(self.caero, msg=msg)
        self.panlst_ref, self.aero_element_ids = cross_reference_panlst(model, self.panlst)

    def safe_cross_reference(self, model: BDF, xref_errors):
        try:
            msg = ', which is required by SPLINE1 eid=%s' % self.eid
            self.setg_ref = model.Set(self.setg, msg=msg)
            self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
        except Exception:
            pass
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        #self.caero_ref = model.CAero(self.caero, msg=msg)
        self.panlst_ref, self.aero_element_ids = cross_reference_panlst(model, self.panlst)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.panlst_ref = None
        self.setg_ref = None

    def raw_fields(self):
        list_fields = ['SPLINE2', self.eid, self.model, self.panlst, self.setg,
                       self.dz, self.eps, self.cp, self.curvature]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SPLINE3_ZAERO(Spline):
    """
    Defines a 3-D spline for the BODY7 and CAERO7 macroelement.

    +---------+------+-------+-------+------+------+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |   6  |  7 |  8  |   9   |
    +=========+======+=======+=======+======+======+====+=====+=======+
    | SPLINE3 | EID  | MODEL |  CP   | SETK | SETG | DZ | EPS |       |
    +---------+------+-------+-------+------+------+----+-----+-------+
    | SPLINE3 | 100  |       |  N/A  |  1   |  10  | 0. |     |       |
    +---------+------+-------+-------+------+------+----+-----+-------+

    """
    type = 'SPLINE3_ZAERO'

    def __init__(self, eid, panlst, setg, model=None, cp=None,
                 dz=None, eps=0.01, comment=''):
        """
        Creates a SPLINE3 card, which is useful for control surface
        constraints.

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
        self.model = model
        self.cp = cp
        self.panlst = panlst
        self.setg = setg
        self.dz = dz
        self.eps = eps
        self.panlst_ref = None
        self.setg_ref = None
        self.aero_element_ids = []

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SPLINE3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        model = blank(card, 2, 'model')
        cp = blank(card, 3, 'cp')

        panlst = integer(card, 4, 'panlst/setk')
        setg = integer(card, 5, 'setg')
        dz = blank(card, 6, 'dz')
        eps = double_or_blank(card, 6, 'eps', 0.01)
        return SPLINE3_ZAERO(eid, panlst, setg, model=model, cp=cp, dz=dz, eps=eps,
                             comment=comment)

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by SPLINE3 eid=%s' % self.eid
        self.setg_ref = model.Set(self.setg, msg=msg)
        self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
        #self.nodes_ref = model.Nodes(self.nodes, msg=msg)
        #self.caero_ref = model.CAero(self.caero, msg=msg)
        self.panlst_ref, self.aero_element_ids = cross_reference_panlst(
            model, self.panlst)

    def safe_cross_reference(self, model: BDF, xref_errors):
        msg = ', which is required by SPLINE3 eid=%s' % self.eid
        try:
            self.setg_ref = model.Set(self.setg, msg=msg)
            self.setg_ref.cross_reference_set(model, 'Node', msg=msg)
        except Exception:
            pass
        self.panlst_ref, self.aero_element_ids = cross_reference_panlst(
            model, self.panlst)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.panlst_ref = None
        self.setg_ref = None

    def convert_to_nastran(self, model):
        return []

    def raw_fields(self):
        """
        +---------+------+-------+-------+------+----+----+-----+-------+
        |    1    |  2   |   3   |   4   |  5   |  6 |  7 |  8  |   9   |
        +=========+======+=======+=======+======+====+====+=====+=======+
        | SPLINE3 | EID  | CAERO | BOXID | COMP | G1 | C1 | A1  | USAGE |
        +---------+------+-------+-------+------+----+----+-----+-------+
        |         |  G2  |  C2   |  A2   | ---- | G3 | C3 | A2  |  ---  |
        +---------+------+-------+-------+------+----+----+-----+-------+
        |         |  G4  |  C4   |  A4   | etc. |    |    |     |       |
        +---------+------+-------+-------+------+----+----+-----+-------+

        """
        list_fields = ['SPLINE3', self.eid, self.model, self.cp, self.panlst, self.setg,
                       self.dz, self.eps]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


def cross_reference_set(model: BDF,
                        set_id: int,
                        msg: str='',
                        which_msg: str='') -> tuple:
    zaero = model.zaero
    if set_id in zaero.setadd:
        set_ref = zaero.setadd[set_id]
    elif set_id in model.sets:
        set_ref = model.sets[set_id]
    else:
        setadd = list(zaero.setadd)
        sets = list(model.sets)
        setadd.sort()
        sets.sort()
        msg = (
            f'{msg}: is not [SETADD, SET1]\n'
            f' - setadd = {setadd}\n'
            f' - sets = {sets}'
        )
        model.log.warning(msg)
        set_ref = None
    return set_ref
