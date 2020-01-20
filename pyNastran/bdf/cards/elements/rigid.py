# pylint: disable=R0902,R0904,R0914
"""
All rigid elements are defined in this file.  This includes:

 * RBAR
 * RBAR1
 * RBE1
 * RBE2
 * RBE3
 * RSPLINE
 * RSSCON

All rigid elements are RigidElement and Element objects.

"""
from __future__ import annotations
import warnings
from itertools import count
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.cards.base_card import Element
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_double, integer_double_or_blank, integer_or_blank,
    double_or_blank, integer_double_or_string, parse_components, components_or_blank,
    blank, string)
from pyNastran.bdf.field_writer_16 import print_card_16
# from pyNastran.bdf.cards.utils import build_table_lines
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class RigidElement(Element):
    def cross_reference(self, model: BDF) -> None:
        pass

class RROD(RigidElement):
    """
    Rigid Pin-Ended Element Connection
    Defines a pin-ended element that is rigid in translation

    +------+-----+----+----+-----+-----+-------+
    |   1  |  2  | 3  | 4  |  5  |  6  |   7   |
    +======+=====+====+====+=====+=====+=======+
    | RROD | EID | GA | GB | CMA | CMB | ALPHA |
    +------+-----+----+----+-----+-----+-------+
    | RROD | 5   | 1  |  2 |     |     | 6.5-6 |
    +------+-----+----+----+-----+-----+-------+
    """
    type = 'RROD'
    _properties = ['dependent_nodes', 'independent_nodes']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nids = [1, 2]
        return RROD(eid, nids, cma=None, cmb=None, alpha=0.0, comment='')

    def __init__(self, eid, nids, cma=None, cmb=None, alpha=0.0, comment=''):
        """
        Creates a RROD element

        Parameters
        ----------
        eid : int
            element id
        nids : List[int, int]
            node ids; connected grid points at ends A and B
        cma / cmb : str; default=None
            dependent DOF
            must be in [None, '1', '2', '3']
            one of them must be None
        alpha : float; default=0.0
            coefficient of thermal expansion
        comment : str; default=''
            a comment for the card
        """
        RigidElement.__init__(self)
        if comment:
            self.comment = comment
        if cma == '0':
            cma = None
        if cmb == '0':
            cmb = None

        self.eid = eid
        self.nodes = nids
        self.cma = cma
        self.cmb = cmb
        self.alpha = alpha
        self.nodes_ref = None

    def validate(self):
        msg = ''
        if self.cma not in [None, '1', '2', '3']:
            msg += "  ga=%r; cma=%r must be in [None, '1', '2', '3']\n" % (
                self.nodes[0], self.cma)
        if self.cmb not in [None, '1', '2', '3']:
            msg += "  gb=%r; cmb=%r must be in [None, '1', '2', '3']\n" % (
                self.nodes[1], self.cmb)
        if self.cma is None and self.cmb is None:
            msg += 'A  ga=%s cma=%r; gb=%s cmb=%r; cma or cmb must be None (not both)' % (
                self.nodes[0], self.cma, self.nodes[1], self.cmb)
        elif self.cma is not None and self.cmb is not None:
            msg += 'D  ga=%s cma=%r; gb=%s cmb=%r; cma or cmb must be None (not both)' % (
                self.nodes[0], self.cma, self.nodes[1], self.cmb)

        if msg:
            raise RuntimeError('Invalid Dependent DOFs\n' + msg.rstrip())

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RROD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        ga = integer(card, 2, 'ga')
        gb = integer(card, 3, 'gb')
        cma = components_or_blank(card, 4, 'cma')
        cmb = components_or_blank(card, 5, 'cmb')
        alpha = double_or_blank(card, 6, 'alpha', 0.0)
        assert len(card) <= 7, 'len(RROD card) = %i\ncard=%s' % (len(card), card)
        return RROD(eid, [ga, gb], cma, cmb, alpha, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a RROD card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        ga = data[1]
        gb = data[2]
        cma = str(data[3])
        cmb = str(data[4])
        alpha = data[5]
        return RROD(eid, [ga, gb], cma, cmb, alpha, comment=comment)

    def Ga(self):
        if self.nodes_ref is None:
            return self.nodes[0]
        return self.nodes_ref[0].nid

    def Gb(self):
        if self.nodes_ref is None:
            return self.nodes[1]
        return self.nodes_ref[1].nid

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by RROD eid=%s' % (self.eid)
        self.nodes_ref = model.Nodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = [self.Ga(), self.Gb()]
        self.nodes_ref = None

    @property
    def independent_nodes(self):
        """gets the independent node ids"""
        if self.cma is None:
            return [self.Ga()]
        return [self.Gb()]

    @property
    def dependent_nodes(self):
        """gets the dependent node ids"""
        if self.cma is not None:
            return [self.Ga()]
        return [self.Gb()]

    def raw_fields(self):
        list_fields = ['RROD', self.eid, self.Ga(), self.Gb(),
                       self.cma, self.cmb, self.alpha]
        return list_fields

    def repr_fields(self):
        alpha = set_blank_if_default(self.alpha, 0.0)
        list_fields = ['RROD', self.eid, self.Ga(), self.Gb(),
                       self.cma, self.cmb, alpha]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class RBAR(RigidElement):
    """
    Defines a rigid bar with six degrees-of-freedom at each end.

    +------+-----+----+----+--------+-----+-----+-----+-------+
    |  1   |  2  |  3 |  4 |    5   |  6  |  7  |  8  |   9   |
    +======+=====+====+====+========+=====+=====+=====+=======+
    | RBAR | EID | GA | GB |  CNA   | CNB | CMA | CMB | ALPHA |
    +------+-----+----+----+--------+-----+-----+-----+-------+
    | RBAR |  5  | 1  |  2 | 123456 |     |     |     | 6.5-6 |
    +------+-----+----+----+--------+-----+-----+-----+-------+
    """
    type = 'RBAR'
    _properties = ['dependent_nodes', 'independent_nodes', 'nodes']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nids = [1, 2]
        cna = '123'
        cnb = '456'
        cma = None
        cmb = None
        return RBAR(eid, nids, cna, cnb, cma, cmb, alpha=0., comment='')

    def __init__(self, eid, nids, cna, cnb, cma, cmb, alpha=0., comment=''):
        """
        Creates a RBAR element

        Parameters
        ----------
        eid : int
            element id
        nids : List[int, int]
            node ids; connected grid points at ends A and B
        cna / cnb : str
            independent DOFs in '123456'
        cma / cmb : str
            dependent DOFs in '123456'
        alpha : float; default=0.0
            coefficient of thermal expansion
        comment : str; default=''
            a comment for the card

        """
        RigidElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.ga = nids[0]
        self.gb = nids[1]

        # If both CNA and CNB are blank, then CNA = 123456.
        if (cna, cnb) == ('', ''):
            cna = '123456'
            #cnb = ''
        elif (cna, cnb) == ('123456', '0'):
            cnb = ''

        self.cna = cna
        self.cnb = cnb

        #  If both CMA and CMB are zero or blank, all of the degrees-of-freedom
        #  not in CNA and CNB will be made dependent.
        if (cma, cmb) == ('', ''):
            for comp in '123456':
                if comp not in cna:
                    cma += comp
                if comp not in cnb:
                    cmb += comp

        self.cma = cma
        self.cmb = cmb
        self.alpha = alpha
        self.ga_ref = None
        self.gb_ref = None

    def validate(self):
        ncna = len(self.cna)
        ncma = len(self.cma)

        ncnb = len(self.cnb)
        ncmb = len(self.cmb)

        nindependent = ncna + ncnb
        ndependent = ncma + ncmb
        independent = self.cna + self.cnb
        dependent = self.cma + self.cmb

        independent_a = set()
        independent_b = set()
        dependent_a = set()
        dependent_b = set()
        msg = ''
        for comp in '123456':
            if comp in self.cna:
                independent_a.add(comp)
            if comp in self.cnb:
                independent_b.add(comp)
            if comp in self.cma:
                if comp in independent_a:
                    msg += 'dof=%s on node %s is both independent and dependent\n' % (self.ga)
                dependent_a.add(comp)
            if comp in self.cmb:
                if comp in independent_b:
                    msg += 'dof=%s on node %s is both independent and dependent\n' % (self.gb)
                dependent_b.add(comp)
        if msg:
            raise RuntimeError(msg + str(self))

        if 0:  # pragma: no cover
            msgi = ''
            msg1 = ''
            if nindependent != 6:
                msg1 = 'nindependent=%s; cna=%r (%s) cnb=%r (%s)\n' % (
                    nindependent, self.cna, ncna, self.cnb, ncnb)
                #raise RuntimeError(msg)
            for comp in '123456':
                if comp not in independent:
                    msgi += '  comp=%s is not independent\n' % (comp)
            if msgi:
                msg1 += 'cna=%r cnb=%r\n%s' % (self.cna, self.cnb, msgi)

            msgi2 = ''
            msg2 = ''
            if nindependent != 6:
                msg = 'ndependent=%s; cma=%r (%s) cmb=%r (%s)\n' % (
                    ndependent, self.cma, ncma, self.cmb, ncmb)
                raise RuntimeError(msg)
            for comp in '123456':
                if comp not in dependent:
                    msgi2 += '  comp=%s is not dependent\n' % (comp)
            if msgi2:
                msg2 = 'cma=%r cmb=%r\n%s' % (self.cma, self.cmb, msgi2)

            msg = msg1 + msg2
            if msg:
                raise RuntimeError(msg + str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RBAR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        ga = integer(card, 2, 'ga')
        gb = integer(card, 3, 'gb')
        cna = components_or_blank(card, 4, 'cna', '')
        cnb = components_or_blank(card, 5, 'cnb', '')
        cma = components_or_blank(card, 6, 'cma', '')
        cmb = components_or_blank(card, 7, 'cmb', '')
        alpha = double_or_blank(card, 8, 'alpha', 0.0)
        assert len(card) <= 9, 'len(RBAR card) = %i\ncard=%s' % (len(card), card)
        return RBAR(eid, [ga, gb], cna, cnb, cma, cmb, alpha, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a RBAR card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        ga = data[1]
        gb = data[2]
        cna = str(data[3])
        cnb = str(data[4])
        cma = str(data[5])
        cmb = str(data[6])
        alpha = data[7]
        return RBAR(eid, [ga, gb], cna, cnb, cma, cmb, alpha, comment=comment)

    # def convert_to_MPC(self, mpcID):
    #     """
    #     -Ai*ui + Aj*uj = 0
    #     where ui are the base DOFs (max=6)
    #     mpc sid   g1 c1 a1  g2 c2 a2
    #     rbe2 eid  gn cm g1  g2 g3 g4
    #     """
    #     raise NotImplementedError()
    #     #i = 0
    #     nCM = len(self.cm)
    #     Ai = nCM * len(self.Gmi) / len(self.gn)  # where nGN=1
    #
    #     card = ['MPC', mpcID]
    #     for cm in self.cm:  # the minus sign is applied to the base node
    #         card += [self.gn, cm, -Ai]
    #
    #     for gm in self.Gmi:
    #         for cm in self.cm:
    #             card += [gm, cm, Ai]
    #     return card

    def Ga(self):
        if self.ga_ref is not None:
            return self.ga_ref.nid
        return self.ga

    def Gb(self):
        if self.gb_ref is not None:
            return self.gb_ref.nid
        return self.gb

    @property
    def nodes(self):
        return [self.Ga(), self.Gb()]

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by RBAR eid=%s' % (self.eid)
        self.ga_ref = model.Node(self.Ga(), msg=msg)
        self.gb_ref = model.Node(self.Gb(), msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.ga_ref = None
        self.gb_ref = None

    @property
    def independent_nodes(self):
        """gets the independent node ids"""
        return [self.Ga()]

    @property
    def dependent_nodes(self):
        """gets the dependent node ids"""
        return [self.Gb()]

    def raw_fields(self):
        list_fields = ['RBAR', self.eid, self.Ga(), self.Gb(), self.cna,
                       self.cnb, self.cma, self.cmb, self.alpha]
        return list_fields

    def repr_fields(self):
        alpha = set_blank_if_default(self.alpha, 0.0)
        list_fields = ['RBAR', self.eid, self.Ga(), self.Gb(), self.cna, self.cnb,
                       self.cma, self.cmb, alpha]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class RBAR1(RigidElement):
    """
    +-------+-----+----+----+-----+-------+
    |   1   |  2  |  3 |  4 |  5  |   6   |
    +=======+=====+====+====+=====+=======+
    | RBAR1 | EID | GA | GB | CB  | ALPHA |
    +-------+-----+----+----+-----+-------+
    | RBAR1 | 5   |  1 |  2 | 123 | 6.5-6 |
    +-------+-----+----+----+-----+-------+
    """
    type = 'RBAR1'

    _properties = ['dependent_nodes', 'independent_nodes', 'nodes']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nids = [1, 2]
        cb = '123'
        return RBAR1(eid, nids, cb, alpha=0., comment='')

    def __init__(self, eid, nids, cb, alpha=0., comment=''):
        RigidElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.ga = nids[0]
        self.gb = nids[1]
        self.cb = cb
        self.alpha = alpha
        self.ga_ref = None
        self.gb_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RBAR1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        ga = integer(card, 2, 'ga')
        gb = integer(card, 3, 'gb')
        cb = components_or_blank(card, 4, 'cb')
        alpha = double_or_blank(card, 5, 'alpha', 0.0)
        assert len(card) <= 6, 'len(RBAR1 card) = %i\ncard=%s' % (len(card), card)
        return RBAR1(eid, [ga, gb], cb, alpha=alpha, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a RBAR1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        ga = data[1]
        gb = data[2]
        cb = data[3]
        alpha = data[4]
        return RBAR1(eid, [ga, gb], cb, alpha=alpha, comment=comment)

    def Ga(self):
        if self.ga_ref is not None:
            return self.ga_ref.nid
        return self.ga

    def Gb(self):
        if self.gb_ref is not None:
            return self.gb_ref.nid
        return self.gb

    @property
    def independent_nodes(self):
        """gets the independent node ids"""
        return [self.Ga()]

    @property
    def dependent_nodes(self):
        """gets the dependent node ids"""
        return [self.Gb()]

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by RBAR1 eid=%s' % (self.eid)
        self.ga_ref = model.Node(self.Ga(), msg=msg)
        self.gb_ref = model.Node(self.Gb(), msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.ga_ref = None
        self.gb_ref = None

    def raw_fields(self):
        list_fields = ['RBAR1', self.eid, self.Ga(), self.Gb(), self.cb, self.alpha]
        return list_fields

    def repr_fields(self):
        alpha = set_blank_if_default(self.alpha, 0.0)
        list_fields = ['RBAR1', self.eid, self.Ga(), self.Gb(), self.cb, alpha]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class RBE1(RigidElement):  # maybe not done, needs testing
    """
    +------+-----+-----+-----+-------+-----+-----+-----+
    |   1  |  2  |  3  |  4  |   5   |  6  |  7  |  8  |
    +======+=====+=====+=====+=======+=====+=====+=====+
    | RBE1 | EID | GN1 | CN1 |  GN2  | CN2 | GN3 | CN3 |
    +------+-----+-----+-----+-------+-----+-----+-----+
    |      |     | GN4 | CN4 |  GN5  | CN5 | GN6 | CN6 |
    +------+-----+-----+-----+-------+-----+-----+-----+
    |      | UM  | GM1 | CM1 |  GM2  | CM2 | GM3 | CM3 |
    +------+-----+-----+-----+-------+-----+-----+-----+
    |      | GM4 | CM4 | etc | ALPHA |     |     |     |
    +------+-----+-----+-----+-------+-----+-----+-----+

    +------+-----+-----+-----+-------+-----+-----+-----+
    | RBE1 | 59  | 59  | 123 |  60   | 456 |     |     |
    +------+-----+-----+-----+-------+-----+-----+-----+
    |      | UM  | 61  | 246 |       |     |     |     |
    +------+-----+-----+-----+-------+-----+-----+-----+
    """
    type = 'RBE1'
    _properties = ['Gmi_node_ids', 'Gni_node_ids', 'dependent_nodes', 'independent_nodes']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        Gni = [2]
        Cni = ['2']
        Gmi = [3]
        Cmi = ['4']
        return RBE1(eid, Gni, Cni, Gmi, Cmi, alpha=0., comment='')

    def __init__(self, eid, Gni, Cni, Gmi, Cmi, alpha=0., comment=''):
        """
        Creates an RBE1 element

        Parameters
        ----------
        eid : int
            element id
        Gni : List[int]
            independent node ids
        Cni : List[str]
            the independent components (e.g., '123456')
        Gmi : List[int]
            dependent node ids
        Cmi : List[str]
            the dependent components (e.g., '123456')
        alpha : float; default=0.
            thermal expansion coefficient
        comment : str; default=''
            a comment for the card
        """
        RigidElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.Gni = Gni
        self.Cni = Cni
        self.Gmi = Gmi
        self.Cmi = Cmi
        self.alpha = alpha
        assert len(Gmi) == len(Cmi), 'len(Gmi)=%s len(Cmi)=%s' % (len(Gmi), len(Cmi))
        assert len(Gmi) > 0, 'len(Gmi)=%s must be greater than 0' % (len(Gmi))
        self.Gni_ref = None
        self.Gmi_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RBE1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        ium = card.index('UM')
        if ium > 0:
            assert string(card, ium, 'UM') == 'UM', 'RBE1=%s must contain UM' % str(card)

        #assert isinstance(card[-1], str), 'card[-1]=%r type=%s' %(card[-1], type(card[-1]))
        alpha_last = integer_double_or_string(card, -1, 'alpha_last')
        if isinstance(alpha_last, float):
            alpha = alpha_last
            card.pop()  # remove the last field so len(card) will not include alpha
        else:
            alpha = 0.

        # loop till UM, no field9,field10
        n = 1
        i = 0
        offset = 2
        Gni = []
        Cni = []
        Gmi = []
        Cmi = []
        while offset + i < ium - 1:
            #print('field(%s) = %s' % (offset + i, card.field(offset + i)))
            gni = integer_or_blank(card, offset + i, 'gn%i' % n)
            cni = components_or_blank(card, offset + i + 1, 'cn%i' % n)

            if gni:
                #print("gni=%s cni=%s" % (gni ,cni))
                Gni.append(gni)
                Cni.append(cni)
                n += 1
            else:
                assert cni is None
            i += 2

        # loop till alpha, no field9,field10
        n = 1
        offset = ium + 1
        i = 0

        # dont grab alpha
        while offset + i < len(card):
            gmi = integer_or_blank(card, offset + i, 'gm%i' % n)
            cmi = components_or_blank(card, offset + i + 1, 'cm%i' % n)
            if gmi:
                Gmi.append(gmi)
                Cmi.append(cmi)
                n += 1
            else:
                assert cmi is None
            i += 2
        return RBE1(eid, Gni, Cni, Gmi, Cmi, alpha=alpha, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by RBE1 eid=%s' % (self.eid)
        self.Gni_ref = model.EmptyNodes(self.Gni, msg=msg)
        self.Gmi_ref = model.EmptyNodes(self.Gmi, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.Gni = self.Gni_node_ids
        self.Gmi = self.Gmi_node_ids
        self.Gni_ref = None
        self.Gmi_ref = None

    @property
    def Gni_node_ids(self):
        if self.Gni_ref is None:
            return self.Gni
        if len(self.Gni_ref) == 0:
            return []
        return self._node_ids(nodes=self.Gni_ref, allow_empty_nodes=True)

    @property
    def Gmi_node_ids(self):
        if self.Gmi_ref is None:
            return self.Gmi
        if len(self.Gmi_ref) == 0:
            return []
        return self._node_ids(nodes=self.Gmi_ref, allow_empty_nodes=True)

    @property
    def independent_nodes(self):
        """gets the independent node ids"""
        # checked
        return self.Gni_node_ids

    @property
    def dependent_nodes(self):
        """gets the dependent node ids"""
        # checked
        nodes = self.Gmi_node_ids
        return nodes

    def raw_fields(self):
        list_fields = [self.type, self.eid]
        for (i, gn, cn) in zip(count(), self.Gni_node_ids, self.Cni):
            list_fields += [gn, cn]
            if i > 0 and i % 3 == 0:
                list_fields += [None]

        nspaces = 8 - (len(list_fields) - 1) % 8  # puts UM/ALPHA onto next line
        if nspaces < 8:
            list_fields += [None] * nspaces

        # overly complicated loop to print the UM section
        list_fields += ['UM']
        j = 1
        for (i, gm, cm) in zip(count(), self.Gmi_node_ids, self.Cmi):
            list_fields += [gm, cm]
            if i > 0 and j % 3 == 0:
                list_fields += [None, None]
                j -= 3
            j += 1

        if self.alpha > 0.:  # handles default alpha value
            nspaces = 8 - (len(list_fields) - 1) % 8  # puts ALPHA onto next line
            if nspaces == 1:
                list_fields += [None, None]
            list_fields += [self.alpha]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class RBE2(RigidElement):
    """
    +-------+-----+-----+-----+------+-------+-----+-----+-----+
    |   1   |  2  |  3  |  4  |  5   |   6   |  7  |  8  |  9  |
    +=======+=====+=====+=====+======+=======+=====+=====+=====+
    |  RBE2 | EID | GN  | CM  | GM1  |  GM2  | GM3 | GM4 | GM5 |
    +-------+-----+-----+-----+------+-------+-----+-----+-----+
    |       | GM6 | GM7 | GM8 | etc. | ALPHA |     |     |     |
    +-------+-----+-----+-----+------+-------+-----+-----+-----+
    """
    type = 'RBE2'
    _field_map = {1: 'eid', 2:'gn', 3:'cm'}
    _properties = ['Gmi_node_ids', 'dependent_nodes', 'independent_nodes']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        gn = 1
        cm = '123'
        Gmi = [2, 3]
        return RBE2(eid, gn, cm, Gmi, alpha=0.0, comment='')

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the GRID card

        Parameters
        ----------
        n : int
            the field number to update
        value : int/float/str
            the value for the appropriate field
        """
        if 3 < n <= 3 + len(self.Gmi):
            self.Gmi[n - 4] = value
        elif n == 4 + len(self.Gmi):
            self.alpha = value
        else:
            raise KeyError('Field %r is an invalid %s entry.' % (n, self.type))
        return value

    def __init__(self, eid, gn, cm, Gmi, alpha=0.0, comment=''):
        """
        Creates an RBE2 element

        Parameters
        ----------
        eid : int
            element id
        gn : int
           Identification number of grid point to which all six independent
           degrees-of-freedom for the element are assigned.
        cm : str
            Component numbers of the dependent degrees-of-freedom in the
            global coordinate system at grid points GMi.
        Gmi : List[int]
            dependent nodes
        alpha : float; default=0.0
            ???
        """
        RigidElement.__init__(self)
        if comment:
            self.comment = comment
        #: Element identification number
        self.eid = eid

        #: Identification number of grid point to which all six independent
        #: degrees-of-freedom for the element are assigned. (Integer > 0)
        self.gn = gn

        #: Component numbers of the dependent degrees-of-freedom in the
        #: global coordinate system at grid points GMi. (Integers 1 through
        #: 6 with no embedded blanks.)
        self.cm = cm

        self.alpha = alpha

        #: Grid point identification numbers at which dependent
        #: degrees-of-freedom are assigned. (Integer > 0)
        if isinstance(Gmi, integer_types):
            Gmi = [Gmi]
        elif isinstance(Gmi, np.ndarray):
            Gmi = Gmi.tolist()
        self.Gmi = Gmi
        #self.nodes_ref = None
        self.Gmi_ref = None
        self.gn_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RBE2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        gn = integer(card, 2, 'gn')
        cm = components_or_blank(card, 3, 'cm')

        alpha = integer_or_double(card, len(card) - 1, 'alpha')
        if isinstance(alpha, float):
            # alpha is correct
            # the last field is not part of Gmi
            n = 1
        else:
            # the last field is part of Gmi
            n = 0
            alpha = 0.0

        j = 4
        Gmi = []
        for i in range(len(card) - 4 - n):
            gmi = integer(card, j + i, 'Gm%i' % (i + 1))
            Gmi.append(gmi)
        return RBE2(eid, gn, cm, Gmi, alpha, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a RBE2 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid = data[0]
        gn = data[1]
        cm = data[2]
        Gmi = data[3]
        alpha = data[4]
        #print("eid=%s gn=%s cm=%s Gmi=%s alpha=%s"
              #% (self.eid, self.gn, self.cm, self.Gmi, self.alpha))
        #raise NotImplementedError('RBE2 data...')
        assert len(Gmi) > 0, Gmi
        return RBE2(eid, gn, cm, Gmi, alpha, comment=comment)

    def update(self, maps):
        """
        Updates the card without xref
        """
        nid_map = maps['node']
        eid_map = maps['element']
        eid2 = eid_map[self.eid]
        gn2 = nid_map[self.gn]
        gm2 = [nid_map[nid] for nid in self.Gmi]
        self.eid = eid2
        self.gn = gn2
        self.Gmi = gm2

    def validate(self):
        assert self.gn is not None, 'gn=%s' % self.gn
        assert self.cm is not None, 'cm=%s' % self.cm
        self.gn = self.gn
        self.cm = str(self.cm)
        assert isinstance(self.alpha, float_types), 'alpha=%r type=%s' % (self.alpha, type(self.alpha))

    def convert_to_mpc(self, mpc_id):
        """
        .. math:: -A_i u_i + A_j u_j = 0

        where :math:`u_i` are the base DOFs (max=6)

         +------+------+----+----+-----+----+----+----+
         |   1  |   2  | 3  | 4  |  5  | 6  | 7  | 8  |
         +======+======+====+====+=====+====+====+====+
         | MPC  | sid  | g1 | c1 | a1  | g2 | c2 | a2 |
         +------+------+----+----+-----+----+----+----+
         | RBE2 | eid  | gn | cm | g1  | g2 | g3 | g4 |
         +------+------+----+----+-----+----+----+----+
        """
        n_cm = len(self.cm)
        Ai = n_cm * len(self.Gmi) / len(self.gn)  # where nGN=1

        card = ['MPC', mpc_id]
        for cm in self.cm:
            # the minus sign is applied to the base node
            card += [self.gn, cm, -Ai]

        for gm in self.Gmi:
            for cm in self.cm:
                card += [gm, cm, Ai]
        return card

    #def convert_to_RBE3(self):
        #raise NotImplementedError()
        #eid = self.eid
        #ref_node = self.gn
        #dof = self.cm
        #wf = 1.0
        #sDof = 123  # this is probably wrong...
        #boundary_nodes = self.Gmi

        ## this is to get the farthest nodes for the UM card
        #boundary_nodes.sort()
        #rbe3_nodes = boundary_nodes

        #rbe3 = ['RBE3', eid, ref_node, dof, wf, sDof] + rbe3_nodes
        #return rbe3

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by RBE2 eid=%s' % (self.eid)
        self.Gmi_ref = model.EmptyNodes(self.Gmi, msg=msg)
        self.gn_ref = model.Node(self.Gn(), msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by RBE2 eid=%s' % (self.eid)
        self.Gmi_ref, unused_missing_nodes = model.safe_empty_nodes(self.Gmi, msg=msg)
        self.gn_ref = model.Node(self.Gn(), msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.Gmi = self.Gmi_node_ids
        self.gn = self.Gn()
        self.Gmi_ref = None
        self.gn_ref = None

    def Gn(self):
        if self.gn_ref is not None:
            return self.gn_ref.nid
        return self.gn

    @property
    def Gmi_node_ids(self):
        if self.Gmi_ref is None or len(self.Gmi) == 0:
            return self.Gmi
        assert self.Gmi_ref is not None, self.Gmi

        # this lets us remove duplicate nodes when we xref
        non_unique_gmi_node_ids = self._node_ids(nodes=self.Gmi_ref, allow_empty_nodes=True)
        return np.unique(non_unique_gmi_node_ids).tolist()

    @property
    def independent_nodes(self):
        """gets the independent node ids"""
        nodes = [self.Gn()]
        return nodes

    @property
    def dependent_nodes(self):
        """gets the dependent node ids"""
        return self.Gmi_node_ids

    def raw_fields(self):
        list_fields = ['RBE2', self.eid, self.Gn(), self.cm] + self.Gmi_node_ids + [self.alpha]
        return list_fields

    def repr_fields(self):
        alpha = set_blank_if_default(self.alpha, 0.)
        list_fields = ['RBE2', self.eid, self.Gn(), self.cm] + self.Gmi_node_ids + [alpha]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_16(card)


class RBE3(RigidElement):
    """
    .. todo:: not done, needs testing badly

    +------+---------+---------+---------+------+--------+--------+------+--------+
    |   1  |    2    |    3    |    4    |  5   |    6   |    7   |   8  |    9   |
    +======+=========+=========+=========+======+========+========+======+========+
    | RBE3 |   EID   |         | REFGRID | REFC |  WT1   |   C1   | G1,1 |  G1,2  |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    |      |   G1,3  |   WT2   |   C2    | G2,1 |  G2,2  |  etc.  | WT3  |   C3   |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    |      |   G3,1  |   G3,2  |  etc.   | WT4  |  C4    |  G4,1  | G4,2 |  etc.  |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    |      |   'UM'  |   GM1   |   CM1   | GM2  |  CM2   |  GM3   | CM3  |        |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    |      |   GM4   |   CM4   |   GM5   | CM5  |  etc.  |        |      |        |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    |      | 'ALPHA' |   ALPHA |         |      |        |        |      |        |
    +------+---------+---------+---------+------+--------+--------+------+--------+
    """
    type = 'RBE3'
    _properties = ['wt_cg_groups', 'ref_grid_id', 'Gijs_node_ids',
                   'dependent_nodes', 'independent_nodes']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        refgrid = 1
        refc = '123'
        weights = [1.]
        comps = ['1']
        Gijs = [2]
        return RBE3(eid, refgrid, refc, weights, comps, Gijs,
                    Gmi=None, Cmi=None, alpha=0.0, comment='')

    def __init__(self, eid, refgrid, refc, weights, comps, Gijs,
                 Gmi=None, Cmi=None, alpha=0.0, comment=''):
        """
        Creates an RBE3 element

        Parameters
        ----------
        eid : int
            element id
        refgrid : int
            dependent node
        refc : str
            dependent components for refgrid???
        GiJs : List[int, ..., int]
            independent nodes
        comps : List[str, ..., str]
            independent components
        weights : List[float, ..., float]
            weights for the importance of the DOF
        Gmi : List[int, ..., int]; default=None -> []
            dependent nodes / UM Set
        Cmi : List[str, ..., str]; default=None -> []
            dependent components / UM Set
        alpha : float; default=0.0
            thermal expansion coefficient
        comment : str; default=''
            a comment for the card

        """
        RigidElement.__init__(self)
        if comment:
            self.comment = comment
        if Gmi is None:
            Gmi = []
        if Cmi is None:
            Cmi = []

        self.eid = eid
        self.refgrid = refgrid
        self.refc = refc
        self.refgrid_ref = None
        self.Gmi_ref = None
        self.Gijs_ref = None

        if not len(weights) == len(comps) and len(weights) == len(Gijs):
            msg = 'len(weights)=%s len(comps)=%s len(Gijs)=%s' % (
                len(weights), len(comps), len(Gijs))
            raise RuntimeError(msg)

        self.weights = weights
        self.comps = comps
        # allow for Gijs as a list or list of lists
        if isinstance(Gijs[0], integer_types):
            Gijs2 = []
            for Gij in Gijs:
                assert isinstance(Gij, integer_types), 'Gij=%s type=%s' % (Gij, type(Gij))
                Gijs2.append([Gij])
            self.Gijs = Gijs2
        else:
            # default
            self.Gijs = Gijs

        if not len(Gmi) == len(Cmi):
            msg = 'len(Gmi)=%s len(Cmi)=%s' % (len(Gmi), len(Cmi))
            raise RuntimeError(msg)
        self.Gmi = Gmi
        self.Cmi = Cmi

        self.alpha = alpha
        self.nodes_ref = None
        self.pid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RBE3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        blank(card, 2, 'blank')
        refgrid = integer(card, 3, 'refgrid')
        refc = components_or_blank(card, 4, 'refc')

        fields = [field.upper() if isinstance(field, str) else field for field in card[5:]]
        ioffset = 5
        iwt_max = len(fields) + ioffset
        try:
            ialpha = fields.index('ALPHA') + ioffset
            iwt_max = ialpha  # the index to start parsing UM
            ium_stop = ialpha  # the index to stop  parsing UM
        except ValueError:
            ialpha = None
            ium_stop = iwt_max

        try:
            ium = fields.index('UM') + ioffset
            iwt_max = ium
        except ValueError:
            ium = None

        i = ioffset
        n = 1
        weights = []
        comps = []
        Gijs = []
        while i < iwt_max:
            Gij = []
            wtname = 'wt' + str(n)
            wt = double_or_blank(card, i, wtname)
            if wt is not None:
                cname = 'c'+str(n)
                compi = components_or_blank(card, i + 1, cname)

                #print("%s=%s %s=%s" % (wtname, wt, cname, compi))
                i += 2
                gij = 0

                j = 0
                while isinstance(gij, int) and i < iwt_max:
                    j += 1
                    gij_name = 'g%s,%s' % (n, j)
                    gij = integer_double_or_blank(card, i, gij_name)
                    if isinstance(gij, float):
                        break
                    #print("%s = %s" % (gij_name, gij))
                    if gij is not None:
                        Gij.append(gij)
                    i += 1
                assert compi is not None
                assert len(Gij) > 0, Gij
                assert Gij[0] is not None, Gij
                weights.append(wt)
                comps.append(compi)
                Gijs.append(Gij)
                #print('----finished a group=%r----' % weight_cg_group)
            else:
                i += 1

        Gmi = []
        Cmi = []
        if ium:
            #print('UM = %s' % card.field(ium))  # UM
            i = ium + 1
            n = 1
            #print("i=%s iUmStop=%s" % (i, iUmStop))
            for j in range(i, ium_stop, 2):

                gm_name = 'gm' + str(n)
                cm_name = 'cm' + str(n)
                gmi = integer_or_blank(card, j, gm_name)
                if gmi is not None:
                    cmi = parse_components(card, j + 1, cm_name)
                    #print("gmi=%s cmi=%s" % (gmi, cmi))
                    Gmi.append(gmi)
                    Cmi.append(cmi)

        if ialpha:
            alpha = double_or_blank(card, ialpha + 1, 'alpha')
        else:
            #: thermal expansion coefficient
            alpha = 0.0
        return RBE3(eid, refgrid, refc, weights, comps, Gijs,
                    Gmi=Gmi, Cmi=Cmi, alpha=alpha, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a RBE3 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        eid, refgrid, refc, weights, comps, gijs, gmi, cmi, alpha = data
        return RBE3(eid, refgrid, refc, weights, comps, gijs,
                    Gmi=gmi, Cmi=cmi, alpha=alpha, comment=comment)

    @property
    def wt_cg_groups(self):
        wt_cg_groups = []
        for weight, comp, gijs in zip(self.weights, self.comps, self.Gijs):
            wt_cg_groups.append((weight, comp, gijs))
        return wt_cg_groups

    # def convert_to_mpc(self, mpc_id):
    #     """
    #     -Ai*ui + Aj*uj = 0
    #     where ui are the base DOFs (max=6)
    #     mpc sid   g1 c1 a1  g2 c2 a2
    #     rbe2 eid  gn cm g1  g2 g3 g4
    #     """
    #     raise NotImplementedError('this is the code for an RBE2...not RBE3')
    #     #i = 0
    #     nCM = len(self.cm)
    #     Ai = nCM * len(self.Gmi) / len(self.gn)  # where nGN=1
    #
    #     card = ['MPC', mpc_id]
    #     for cm in self.cm:  # the minus sign is applied to the base node
    #         card += [self.gn, cm, -Ai]
    #
    #     for gm in self.Gmi:
    #         for cm in self.cm:
    #             card += [gm, cm, Ai]
    #     return card

    @property
    def ref_grid_id(self):
        if self.refgrid_ref is not None:
            return self.refgrid_ref.nid
        return self.refgrid

    @property
    def Gmi_node_ids(self):
        if self.Gmi_ref is None:
            return self.Gmi
        if len(self.Gmi_ref) == 0:
            return []
        return self._node_ids(nodes=self.Gmi_ref, allow_empty_nodes=True)

    @property
    def Gijs_node_ids(self):
        if self.Gijs_ref is None:
            return self.Gijs
        Gijs = []
        for gij in self.Gijs_ref:
            gijs = self._node_ids(nodes=gij, allow_empty_nodes=True)
            Gijs.append(gijs)
        return Gijs

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by RBE3 eid=%s' % (self.eid)
        assert self.Gmi is not None
        self.Gmi_ref = model.EmptyNodes(self.Gmi, msg=msg)

        assert self.Gmi is not None
        self.refgrid_ref = model.Node(self.ref_grid_id, msg=msg)

        self.Gijs_ref = []
        for Gij in self.Gijs:
            self.Gijs_ref.append(model.EmptyNodes(Gij, msg=msg))

    def safe_cross_reference(self, model, debug=True):
        msg = ', which is required by RBE3 eid=%s' % (self.eid)
        assert self.Gmi is not None
        self.Gmi_ref, unused_missing_nodes = model.safe_empty_nodes(self.Gmi, msg=msg)

        assert self.Gmi_ref is not None
        self.refgrid_ref = model.Node(self.ref_grid_id, msg=msg)

        self.Gijs_ref = []
        for Gij in self.Gijs:
            nodes, msgi = model.safe_empty_nodes(Gij, msg=msg)
            if msgi:
                model.log.warning(msgi)
            self.Gijs_ref.append(nodes)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.Gijs = self.Gijs_node_ids
        self.Gmi = self.Gmi_node_ids
        self.refgrid = self.ref_grid_id
        self.Gijs_ref = None
        self.refgrid_ref = None
        self.Gmi_ref = None

        Gij = []
        for gij in self.Gijs:
            gij = self._node_ids(nodes=gij, allow_empty_nodes=True)
            Gij.append(gij)
        self.Gijs = Gij

    @property
    def independent_nodes(self):
        """
        gets the independent node ids
        TODO: not checked
        """
        nodes = []
        for gij in self.Gijs:
            giji = self._node_ids(nodes=gij, allow_empty_nodes=True)
            nodes += giji
        return nodes

    @property
    def dependent_nodes(self):
        """
        gets the dependent node ids
        TODO: not checked
        """
        nodes = [self.ref_grid_id]
        nodes += self.Gmi_node_ids
        return nodes

    def raw_fields(self):
        list_fields = ['RBE3', self.eid, None, self.ref_grid_id, self.refc]
        for (wt, ci, Gij) in zip(self.weights, self.comps, self.Gijs_node_ids):
            list_fields += [wt, ci] + Gij
        nspaces = 8 - (len(list_fields) - 1) % 8  # puts UM onto next line

        if nspaces < 8:
            list_fields += [None] * nspaces

        if self.Gmi:
            list_fields += ['UM']
            for (gmi, cmi) in zip(self.Gmi_node_ids, self.Cmi):
                list_fields += [gmi, cmi]

        nspaces = 8 - (len(list_fields) - 1) % 8  # puts ALPHA onto next line
        if nspaces < 8:
            list_fields += [None] * nspaces

        if self.alpha > 0.:  # handles the default value
            list_fields += ['ALPHA', self.alpha]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class RSPLINE(RigidElement):
    type = 'RSPLINE'
    _properties = ['dependent_nodes', 'independent_nodes']
    """
    Defines multipoint constraints for the interpolation of displacements
    at grid points.

    +---------+-----+-----+----+----+--------+----+----+----+
    |    1    |  2  |  3  |  4 |  5 |    6   |  7 |  8 |  9 |
    +=========+=====+=====+====+====+========+====+====+====+
    | RSPLINE | EID | D/L | G1 | G2 |   C2   | G3 | C3 | G4 |
    +---------+-----+-----+----+----+--------+----+----+----+
    |         |  C4 |  G5 | C5 | G6 |  etc.  |    |    |    |
    +---------+-----+-----+----+----+--------+----+----+----+
    """
    @classmethod
    def _init_from_empty(cls):
        eid = 1
        independent_nid = 2
        dependent_nids = [3, 4]
        dependent_components = ['4', '5']
        return RSPLINE(eid, independent_nid, dependent_nids, dependent_components,
                       diameter_ratio=0.1, comment='')

    def __init__(self, eid, independent_nid, dependent_nids, dependent_components,
                 diameter_ratio=0.1, comment=''):
        """
        Creates a RSPLINE card, which uses multipoint constraints for the
        interpolation of displacements at grid points

        Parameters
        ----------
        eid : int
            element id
        independent_nid : int
            the independent node id
        dependent_nids : List[int]
            the dependent node ids
        dependent_components : List[str]
            Components to be constrained
        diameter_ratio : float; default=0.1
            Ratio of the diameter of the elastic tube to the sum of the
            lengths of all segments
        comment : str; default=''
            a comment for the card
        """
        RigidElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.independent_nid = independent_nid
        self.dependent_nids = dependent_nids
        # Components to be constrained
        self.dependent_components = dependent_components

        # Ratio of the diameter of the elastic tube to the sum of the
        # lengths of all segments
        self.diameter_ratio = diameter_ratio

    def validate(self):
        assert len(self.dependent_nids) == len(self.dependent_components)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RSPLINE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        diameter_ratio = double_or_blank(card, 2, 'diameter_ratio', 0.1)
        nfields = len(card)
        #assert (nfields) % 2 == 1, 'nfields=%s card=%s'  % (nfields, card)
        #assert (nfields - 4) % 2 == 0, 'nfields=%s card=%s'  % (nfields, card)

        # blanks are allowed
        if (nfields - 4) % 2 != 0:
            nfields += 1

        dependent_nids = []
        dependent_components = []
        independent_nid = integer(card, 3, 'nid_1')
        j = 2
        for i in range(4, nfields, 2):
            nid = integer(card, i, 'nid_%s' % j)
            comp = components_or_blank(card, i+1, 'components_%i' % j, default='')
            dependent_nids.append(nid)
            dependent_components.append(comp)
            j += 1
        return RSPLINE(eid, independent_nid, dependent_nids, dependent_components,
                       diameter_ratio=diameter_ratio, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        return
        #msg = ', which is required by RSPLINE eid=%s' % (self.eid)
        #self.Gni = model.EmptyNodes(self.Gni, msg=msg)
        #self.Gmi = model.EmptyNodes(self.Gmi, msg=msg)
        #self.Gni_ref = self.Gni
        #self.Gmi_ref = self.Gmi

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass
        #self.Gni = self.Gni_node_ids
        #self.Gmi = self.Gmi_node_ids
        #del self.Gni_ref, self.Gmi_ref

    @property
    def independent_nodes(self):
        """gets the independent node ids"""
        # TODO: not quite right as it doesn't support blank entries
        return [self.independent_nid]

    @property
    def dependent_nodes(self):
        """gets the dependent node ids"""
        # TODO: not quite right as it doesn't support blank entries
        return self.dependent_nids

    def raw_fields(self):
        list_fields = [self.type, self.eid, self.diameter_ratio, self.independent_nid]
        for (gn, cn) in zip(self.dependent_nids, self.dependent_components):
            list_fields += [gn, cn]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class RSSCON(RigidElement):
    type = 'RSSCON'
    """
    Defines multipoint constraints to model clamped connections
    of shell-to-solid elements.

    +--------+------+------+-----+-----+-----+-----+-----+-----+
    |    1   |   2  |   3  |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+======+======+=====+=====+=====+=====+=====+=====+
    | RSSCON | RBID | TYPE | ES1 | EA1 | EB1 | ES2 | EA2 | EB2 |
    +--------+------+------+-----+-----+-----+-----+-----+-----+
    | RSSCON |  110 | GRID |  11 |  12 |  13 |  14 |  15 |  16 |
    +--------+------+------+-----+-----+-----+-----+-----+-----+
    | RSSCON |  111 | GRID |  31 |  74 |  75 |     |     |     |
    +--------+------+------+-----+-----+-----+-----+-----+-----+
    | RSSCON |  115 | ELEM | 311 | 741 |     |     |     |     |
    +--------+------+------+-----+-----+-----+-----+-----+-----+
    | RSSCON |  116 | INTC |  2  |  1  |  3  |     |     |     |
    +--------+------+------+-----+-----+-----+-----+-----+-----+
    """

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        rigid_type = 'GRID'
        return RSSCON(eid, rigid_type,
                      shell_eid=None, solid_eid=None, a_solid_grids=None,
                      b_solid_grids=None, shell_grids=None, comment='')

    def __init__(self, eid, rigid_type,
                 shell_eid=None, solid_eid=None,
                 a_solid_grids=None, b_solid_grids=None, shell_grids=None,
                 comment=''):
        """
        Creates an RSSCON card, which defines multipoint constraints to
        model clamped connections of shell-to-solid elements.

        Parameters
        ----------
        eid : int
            element id
        rigid_type : str
            GRID/ELEM
        shell/solid_eid : int; default=None
            the shell/solid element id (if rigid_type=ELEM)
        shell/solid_grids : List[int, int]; default=None
            the shell/solid node ids (if rigid_type=GRID)
        comment : str; default=''
            a comment for the card

        B----S----A
        """
        RigidElement.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.rigid_type = rigid_type
        if rigid_type == 'ELEM':
            self.shell_eid = shell_eid
            self.solid_eid = solid_eid
            self.a_solid_grids = None
            self.b_solid_grids = None
            self.shell_grids = None
        elif rigid_type == 'GRID':
            self.shell_eid = None
            self.solid_eid = None
            self.a_solid_grids = a_solid_grids
            self.b_solid_grids = b_solid_grids
            self.shell_grids = shell_grids
        elif rigid_type == 'INTC':
            self.shell_eid = None
            self.solid_eid = None
            self.shell_grids = shell_grids
            self.a_solid_grids = None
            self.b_solid_grids = None
        else:
            #| RSSCON |  116 | INTC |  2  |  1  |  3  |     |     |     |
            raise RuntimeError('rigid_type=%s and must be [ELEM, GRID]' % rigid_type)
        self.shell_eid_ref = None
        self.solid_eid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RSSCON card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        rigid_type = string(card, 2, 'rigid_type')

        if rigid_type == 'ELEM':
            a_solid_grids = None
            b_solid_grids = None
            shell_grids = None
            # ES1, EA1
            shell_eid = integer(card, 3, 'shell_eid')  # ES1
            solid_eid = integer(card, 4, 'solid_eid')  # EA1
            assert len(card) == 5, card
        elif rigid_type == 'GRID':
            shell_eid = None
            solid_eid = None
            # ES1, EA1, EB1
            shell_grids = [integer(card, 3, 'shell_nid_1')]  # ES1
            a_solid_grids = [integer(card, 4, 'a_solid_grid_1')]  # EA1
            b_solid_grids = [integer_or_blank(card, 5, 'b_solid_grid_1')]  # EB1

            shell_grids.append(integer_or_blank(card, 6, 'shell_nid_2'))  # ES2
            a_solid_grids.append(integer_or_blank(card, 7, 'a_solid_grid_2'))  # EA2
            b_solid_grids.append(integer_or_blank(card, 8, 'b_solid_grid_2'))  # EA2
            assert len(card) <= 9, card
        elif  rigid_type == 'INTC':
            shell_eid = None
            solid_eid = None
            shell_grids = [
                integer(card, 3, 'RSSCON INTC field 3'),
                integer(card, 4, 'RSSCON INTC field 4'),
                integer(card, 5, 'RSSCON INTC field 5'),
            ]
            a_solid_grids = None
            b_solid_grids = None
            assert len(card) == 6, card
        else:
            msg = 'RSSCON; eid=%s rigid_type=%s and must be [ELEM, GRID, INTC]' % (eid, rigid_type)
            raise RuntimeError(msg)
        return RSSCON(eid, rigid_type,
                      shell_eid=shell_eid, solid_eid=solid_eid,
                      a_solid_grids=a_solid_grids, b_solid_grids=b_solid_grids,
                      shell_grids=shell_grids,
                      comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        unused_msg = ', which is required by RSSCON eid=%s' % (self.eid)
        #if self.rigid_type == 'ELEM':
            #self.shell_eid_ref = model.Element(self.shell_eid, msg=msg)
            #self.solid_eid_ref = model.Element(self.shell_eid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.cross_reference(model)

    def EidShell(self):
        if self.shell_eid_ref is not None:
            return self.shell_eid_ref.eid
        return self.shell_eid

    def EidSolid(self):
        if self.solid_eid_ref is not None:
            return self.solid_eid_ref.eid
        return self.solid_eid

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.shell_eid = self.EidShell()
        self.solid_eid = self.EidSolid()
        self.shell_eid_ref = None
        self.solid_eid_ref = None

    @property
    def independent_nodes(self):
        """gets the independent node ids"""
        return []

    @property
    def dependent_nodes(self):
        """gets the dependent node ids"""
        return []

    def raw_fields(self):
        list_fields = ['RSSCON', self.eid, self.rigid_type]
        if self.rigid_type == 'ELEM':
            list_fields += [self.EidShell(), self.EidSolid()]
        elif self.rigid_type == 'GRID':
            for nid_shell, nid_a, nid_b in zip(self.shell_grids,
                                               self.a_solid_grids, self.b_solid_grids):
                list_fields += [nid_shell, nid_a, nid_b]
        else:
            raise RuntimeError('rigid_type=%s and must be [ELEM, GRID]' % self.rigid_type)
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
