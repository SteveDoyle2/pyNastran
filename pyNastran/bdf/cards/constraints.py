"""
All constraint cards are defined in this file.  This includes:

* Constraint
 * SUPORT
 * SUPORT1
 * SPC
 * SPC1
 * SPCAX
 * MPC
 * GMSPC
 * ConstraintADD
  * SPCADD
  * MPCADD

The ConstraintObject contain multiple constraints.

"""
from __future__ import annotations
from itertools import count
import warnings
from typing import TYPE_CHECKING

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf import MAX_INT
from pyNastran.bdf.cards.expand_card import remove_blanks_capitalize, set_list_print
from pyNastran.bdf.cards.base_card import BaseCard, _node_ids, expand_thru
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, parse_components,
    components_or_blank, string)
from pyNastran.bdf.field_writer_8 import print_card_8, print_float_8
from pyNastran.bdf.field_writer_16 import print_float_16, print_card_16
from pyNastran.bdf.field_writer_double import print_scientific_double
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class Constraint(BaseCard):
    """
    common class for:
     - SUPORT / SUPORT1 / SESUP
     - GMSPC
     - MPC
     - SPC / SPC1
     - SPCAX
     - SPCOFF / SPCOFF1

     """
    def __init__(self):
        pass

    def _node_ids(self, nodes=None, allow_empty_nodes=False, msg=''):
        """returns nodeIDs for repr functions"""
        return _node_ids(self, nodes, allow_empty_nodes, msg)


class SUPORT1(Constraint):
    """
    Defines determinate reaction degrees-of-freedom (r-set) in a free
    body-analysis.  SUPORT1 must be requested by the SUPORT1 Case
    Control command.

    +---------+-----+-----+----+-----+----+-----+----+
    |    1    |  2  |  3  |  4 |  5  | 6  |  7  | 8  |
    +=========+=====+=====+====+=====+====+=====+====+
    | SUPORT1 | SID | ID1 | C1 | ID2 | C2 | ID3 | C3 |
    +---------+-----+-----+----+-----+----+-----+----+
    | SUPORT1 |  1  |  2  | 23 |  4  | 15 |  5  |  0 |
    +---------+-----+-----+----+-----+----+-----+----+
    """
    type = 'SUPORT1'
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        conid = 1
        nodes = [1]
        components = ['123']
        return SUPORT1(conid, nodes, components)

    def __init__(self, conid, nodes, components, comment: str=''):
        """
        Creates a SUPORT card, which defines free-body reaction points.

        Parameters
        ----------
        conid : int
            Case Control SUPORT id
        nodes : list[int]
            the nodes to release
        components : list[str]
            components to support at each node
        comment : str; default=''
            a comment for the card

        """
        Constraint.__init__(self)
        if comment:
            self.comment = comment
        self.conid = conid
        self.nodes = nodes
        self.Cs = components
        assert len(self.nodes) > 0
        assert len(self.nodes) == len(self.Cs)
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SUPORT1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        conid = integer(card, 1, 'conid')  # really a support id sid

        nfields = len(card)
        assert len(card) > 2
        nterms = int((nfields - 1.) / 2.)
        n = 1
        nodes = []
        components = []
        for i in range(nterms):
            nstart = 2 + 2 * i
            nid = integer(card, nstart, 'ID%s' % n)
            comp = components_or_blank(card, nstart + 1, 'component%s' % n, '0')
            nodes.append(nid)
            components.append(comp)
            n += 1
        return SUPORT1(conid, nodes, components, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a SUPORT1 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        conid = data[0]
        assert (len(data) - 1) % 2 == 0, data
        nodes = []
        components = []
        for i in range(1, len(data), 2):
            nid = data[i]
            comp = data[i+1]
            nodes.append(nid)
            components.append(comp)
        return SUPORT1(conid, nodes, components, comment=comment)

    def add_suport1_to_set(self, suport1):
        assert self.conid == suport1.conid, 'SUPORT1 conid=%s new_conid=%s; they must be the same' % (self.conid, suport1.conid)
        comment = self.comment + suport1.comment
        if comment:
            self.comment = comment
        self.nodes += suport1.nodes
        self.Cs += suport1.Cs

    @property
    def node_ids(self) -> list[int]:
        msg = ', which is required by SUPORT1'
        if self.nodes_ref is None:
            return self.nodes
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SUPORT1'
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model: BDF, debug=True):
        nids2 = []
        msg = ', which is required by SUPORT1=%s' % self.conid
        for nid in self.nodes:
            try:
                nid2 = model.Node(nid, msg=msg)
            except KeyError:
                if debug:
                    msg = f'Couldnt find nid={nid}, which is required by SUPORT1={self.conid}'
                    model.log.warning(msg)
                continue
            nids2.append(nid2)
        self.nodes_ref = nids2

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def raw_fields(self):
        fields = ['SUPORT1', self.conid]
        for nid, comp in zip(self.node_ids, self.Cs):
            fields += [nid, comp]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class SUPORT(Constraint):
    """
    Defines determinate reaction degrees-of-freedom in a free body.

    +---------+-----+-----+-----+-----+-----+-----+-----+----+
    |    1    |  2  |  3  |  4  |  5  |  6  |  7  |  8  | 9  |
    +=========+=====+=====+=====+=====+=====+=====+=====+====+
    | SUPORT  | ID1 | C1  | ID2 |  C2 | ID3 | C3  | ID4 | C4 |
    +---------+-----+-----+-----+-----+-----+-----+-----+----+
    """
    type = 'SUPORT'
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        nodes = [1, 2]
        components = ['123', '456']
        return SUPORT(nodes, components)

    def __init__(self, nodes, components, comment: str=''):
        """
        Creates a SUPORT card, which defines free-body reaction points.
        This is always active.

        Parameters
        ----------
        nodes : list[int]
            the nodes to release
        components : list[str]
            components to support at each node
        comment : str; default=''
            a comment for the card

        """
        Constraint.__init__(self)
        if comment:
            self.comment = comment
        self.nodes = nodes
        self.Cs = []
        for compi in components:
            if isinstance(compi, integer_types):
                compi = str(compi)
            self.Cs.append(compi)
        self.nodes_ref = None

    def validate(self):
        assert len(self.nodes) > 0
        assert len(self.nodes) == len(self.Cs)

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SUPORT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nfields = len(card)
        assert len(card) > 1, card
        nterms = nfields // 2
        n = 1
        nodes = []
        components = []
        for i in range(nterms):
            nstart = 1 + 2 * i
            nid = integer(card, nstart, 'ID%s' % n)
            componentsi = components_or_blank(card, nstart + 1, 'component%s' % n, '0')
            nodes.append(nid)
            components.append(componentsi)
            n += 1
        return SUPORT(nodes, components, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a SUPORT card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        fields = data
        nodes = []
        components = []
        for i in range(0, len(fields), 2):
            nid, comp = fields[i:i+2]
            assert nid > 0, nid
            nodes.append(nid)
            components.append(comp)
        return SUPORT(nodes, components, comment=comment)

    @property
    def node_ids(self):
        msg = ', which is required by SUPORT'
        if self.nodes_ref is None:
            return self.nodes
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SUPORT'
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model: BDF, debug=True):
        nids2 = []
        msg = ', which is required by SUPORT'
        for nid in self.nodes:
            try:
                nid2 = model.Node(nid, msg=msg)
            except KeyError:
                if debug:
                    msg = 'Couldnt find nid=%i, which is required by SUPORT' % nid
                    model.log.warning(msg)
                continue
            nids2.append(nid2)
        self.nodes_ref = nids2

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def raw_fields(self):
        fields = [self.type]
        for nid, c in zip(self.node_ids, self.Cs):
            fields += [nid, c]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        if max(self.node_ids) > MAX_INT:
            return self.comment + print_card_16(card)
        return self.comment + print_card_8(card)


class SESUP(SUPORT):
    type = 'SESUP'

    @classmethod
    def _init_from_empty(cls):
        nodes = [1, 2]
        components = ['1', '2']
        return SESUP(nodes, components)

    def __init__(self, nodes, components, comment: str=''):
        SUPORT.__init__(self, nodes, components, comment=comment)


class MPC(Constraint):
    """
    Multipoint Constraint
    Defines a multipoint constraint equation of the form:
      sum(A_j * u_j) = 0

    where:
      uj represents degree-of-freedom Cj at grid or scalar point Gj.
      Aj represents the scale factor

    +-----+-----+----+----+-----+----+----+----+-----+
    |  1  |  2  |  3 |  4 |  5  |  6 |  7 |  8 |  9  |
    +=====+=====+====+====+=====+====+====+====+=====+
    | MPC | SID | G1 | C1 |  A1 | G2 | C2 | A2 |     |
    +-----+-----+----+----+-----+----+----+----+-----+
    |     |  G3 | C3 | A3 | ... |    |    |    |     |
    +-----+-----+----+----+-----+----+----+----+-----+
    """
    type = 'MPC'
    #'constraints', 'enforced', 'gids_ref', 'gids'
    _properties = ['node_ids', ]

    @classmethod
    def _init_from_empty(cls):
        conid = 1
        nodes = [1]
        components = ['1']
        coefficients = [1.]
        return MPC(conid, nodes, components, coefficients)

    def __init__(self, conid: int, nodes: list[int], components: list[str],
                 coefficients: list[float], comment: str=''):
        """
        Creates an MPC card

        Parameters
        ----------
        conid : int
            Case Control MPC id
        nodes : list[int]
            GRID/SPOINT ids
        components : list[str]
            the degree of freedoms to constrain (e.g., '1', '123')
        coefficients : list[float]
            the scaling coefficients

        """
        Constraint.__init__(self)
        if comment:
            self.comment = comment
        #: Set identification number. (Integer > 0)
        self.conid = conid

        #: Identification number of grid or scalar point. (Integer > 0)
        self.nodes = nodes

        #: Component number. (Any one of the Integers 1 through 6 for grid
        #: points; blank or zero for scalar points.)
        self.components = components

        #: Coefficient. (Real; Default = 0.0 except A1 must be nonzero.)
        self.coefficients = coefficients

        self.nodes_ref = None

    def validate(self):
        assert isinstance(self.nodes, list), type(self.nodes)
        assert isinstance(self.components, list), type(self.components)
        assert isinstance(self.coefficients, list), type(self.coefficients)
        assert len(self.nodes) == len(self.components)
        assert len(self.nodes) == len(self.coefficients)
        for nid, comp, coefficient in zip(self.nodes, self.components, self.coefficients):
            assert isinstance(nid, integer_types), self.nodes
            assert isinstance(comp, str), self.components
            assert isinstance(coefficient, float), self.coefficients

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds an MPC card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        conid = integer(card, 1, 'conid')
        nodes = []
        components = []
        coefficients = []

        nfields = len(card)

        i = 1
        for ifield in range(2, nfields, 8):
            nid = integer(card, ifield, 'G%i' % i)
            component = components_or_blank(card, ifield + 1, 'constraint%i' % i, '0')  # scalar point
            if i == 1:
                coefficient = double(card, ifield + 2, 'coefficient%i' % i)
                if coefficient == 0.0:
                    raise RuntimeError('coefficient1 must be nonzero; coefficient=%r' % coefficient)
            else:
                coefficient = double_or_blank(card, ifield + 2, 'coefficient%i' % i, 0.0)
            nodes.append(nid)
            components.append(component)
            coefficients.append(coefficient)
            i += 1

            if ifield + 4 > nfields and i != 2:
                # if G2 is empty (it's ifield+4 because nfields is length based
                # and not loop friendly)
                break
            nid = integer(card, ifield + 3, 'G%i' % i)
            component = components_or_blank(card, ifield + 4, 'constraint%i' % i, '0')  # scalar point
            coefficient = double_or_blank(card, ifield + 5, 'coefficient%i' % i, 0.0)
            nodes.append(nid)
            components.append(component)
            coefficients.append(coefficient)
            i += 1
        return MPC(conid, nodes, components, coefficients, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds an MPC card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        #msg = 'MPC has not implemented data parsing'
        conid = data[0]
        nodes = data[1]
        components = [str(component) for component in data[2]]
        enforced = data[3]
        return MPC(conid, nodes, components, enforced, comment=comment)

    @property
    def independent_dofs(self) -> tuple[list[int], list[int]]:
        """The first degree-of-freedom (G1, C1) in the sequence is defined to be the
        dependent degree-of-freedom. A dependent degree-of-freedom assigned by
        one MPC entry cannot be assigned dependent by another MPC entry or by a
        rigid element."""
        nodes = self.nodes[1:]
        components = self.components[1:]
        if isinstance(nodes, integer_types):
            return [nodes], [components]
        return nodes, components

    @property
    def dependent_dofs(self) -> tuple[list[int], list[int]]:
        """The first degree-of-freedom (G1, C1) in the sequence is defined to be the
        dependent degree-of-freedom. A dependent degree-of-freedom assigned by
        one MPC entry cannot be assigned dependent by another MPC entry or by a
        rigid element."""
        return [self.nodes[0]], [self.components[0]]

    @property
    def independent_nodes(self) -> list[int]:
        """The first degree-of-freedom (G1, C1) in the sequence is defined to be the
        dependent degree-of-freedom. A dependent degree-of-freedom assigned by
        one MPC entry cannot be assigned dependent by another MPC entry or by a
        rigid element."""
        nodes = self.nodes[1:]
        if isinstance(nodes, integer_types):
            return [nodes]
        return nodes

    @property
    def dependent_nodes(self) -> list[int]:
        """The first degree-of-freedom (G1, C1) in the sequence is defined to be the
        dependent degree-of-freedom. A dependent degree-of-freedom assigned by
        one MPC entry cannot be assigned dependent by another MPC entry or by a
        rigid element."""
        return [self.nodes[0]]

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        msg = ', which is required by MPC=%s' % self.conid
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MPC=%s' % self.conid
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model: BDF, debug=True):
        nids2 = []
        msg = ', which is required by SPC=%s' % self.conid
        for nid in self.nodes:
            try:
                nid2 = model.Node(nid, msg=msg)
            except KeyError:
                if debug:
                    msg = 'Couldnt find nid=%i, which is required by SPC=%s' % (
                        nid, self.conid)
                    model.log.warning(msg)
                continue
            nids2.append(nid2)
        self.nodes_ref = nids2

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def raw_fields(self):  # MPC
        fields = ['MPC', self.conid]
        for i, gid, component, coefficient in zip(count(), self.node_ids, self.components, self.coefficients):
            fields += [gid, component, coefficient]
            if i % 2 == 1 and i > 0:
                fields.append(None)
                fields.append(None)
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """see BaseCard.write_card``"""
        if size == 8:
            return self.write_card_8()
        return self.write_card_16(is_double)

    def write_card_8(self):
        msg = 'MPC     %8d' % self.conid
        grids, components, coefficients = self.node_ids, self.components, self.coefficients
        for i, grid, component, coefficient in zip(count(), grids, components, coefficients):
            msg += '%8d%8s%8s' % (grid, component, print_float_8(coefficient))
            if i % 2 == 1 and i > 0:
                msg += '\n%8s%8s' % ('', '')
        return self.comment + msg.rstrip() + '\n'

    def write_card_16(self, is_double=False):
        msg = 'MPC*    %16d' % self.conid
        grids, components, coefficients = self.node_ids, self.components, self.coefficients
        if is_double:
            for i, grid, component, coefficient in zip(count(), grids, components, coefficients):
                if i == 0:
                    msg += '%16d%16s%16s\n' % (
                        grid, component, print_scientific_double(coefficient))
                elif i % 2 == 1:
                    msg += '%-8s%16d%16s%16s\n' % (
                        '*', grid, component, print_scientific_double(coefficient))
                else:
                    msg += '%-8s%16s%16d%16s%16s\n' % (
                        '*', '', grid, component, print_scientific_double(coefficient))
        else:
            for i, grid, component, coefficient in zip(count(), grids, components, coefficients):
                if i == 0:
                    msg += '%16d%16s%16s\n' % (grid, component, print_float_16(coefficient))
                elif i % 2 == 1:
                    msg += '%-8s%16d%16s%16s\n' % (
                        '*', grid, component, print_float_16(coefficient))
                else:
                    msg += '%-8s%16s%16d%16s%16s\n' % (
                        '*', '', grid, component, print_float_16(coefficient))
        if i % 2 == 0:
            msg += '*'
        return self.comment + msg.rstrip() + '\n'


class SPC(Constraint):
    """
    Defines enforced displacement/temperature (static analysis)
    velocity/acceleration (dynamic analysis).

     +-----+-----+----+----+------+----+----+----+
     |  1  |  2  |  3 |  4 |   5  |  6 |  7 |  8 |
     +=====+=====+====+====+======+====+====+====+
     | SPC | SID | G1 | C1 |  D1  | G2 | C2 | D2 |
     +-----+-----+----+----+------+----+----+----+
     | SPC |  2  | 32 | 3  | -2.6 |  5 |    |    |
     +-----+-----+----+----+------+----+----+----+
    """
    type = 'SPC'
    _properties = ['node_ids', 'constraints', 'gids_ref', 'gids']

    @classmethod
    def _init_from_empty(cls):
        conid = 1
        nodes = [1, 2]
        components = ['123', '456']
        enforced = [0., 0.]
        return SPC(conid, nodes, components, enforced)

    def __init__(self, conid: int,
                 nodes: list[int],
                 components: list[str],
                 enforced: list[float], comment: str=''):
        """
        Creates an SPC card, which defines the degree of freedoms to be
        constrained

        Parameters
        ----------
        conid : int
            constraint id
        nodes : list[int]
            GRID/SPOINT ids
        components : list[str]
            the degree of freedoms to constrain (e.g., '1', '123')
        enforced : list[float]
            the constrained value for the given node (typically 0.0)
        comment : str; default=''
            a comment for the card

        .. note:: len(nodes) == len(components) == len(enforced)
        .. warning:: non-zero enforced deflection requires an SPCD as well

        """
        Constraint.__init__(self)
        if comment:
            self.comment = comment

        if isinstance(nodes, int):
            nodes = [nodes]

        if isinstance(components, str):
            components = [components]
        elif isinstance(components, int):
            components = [str(components)]

        if isinstance(enforced, float):
            enforced = [enforced]

        self.conid = conid
        self.nodes = nodes
        self.components = components
        self.enforced = enforced
        self.nodes_ref = None

    def validate(self):
        assert isinstance(self.nodes, list), self.nodes
        assert isinstance(self.components, list), self.components
        assert isinstance(self.enforced, list), self.enforced
        assert len(self.nodes) == len(self.components), 'len(self.nodes)=%s len(self.components)=%s' % (len(self.nodes), len(self.components))
        assert len(self.nodes) == len(self.enforced), 'len(self.nodes)=%s len(self.enforced)=%s' % (len(self.nodes), len(self.enforced))
        for nid, comp, enforcedi in zip(self.nodes, self.components, self.enforced):
            assert isinstance(nid, integer_types), self.nodes
            assert isinstance(comp, str), self.components
            assert isinstance(enforcedi, float), self.enforced

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds an SPC card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        conid = integer(card, 1, 'sid')
        if card.field(5) in [None, '']:
            nodes = [integer(card, 2, 'G1'),]
            components = [components_or_blank(card, 3, 'C1', '0')]
            enforced = [double_or_blank(card, 4, 'D1', 0.0)]
        else:
            nodes = [
                integer(card, 2, 'G1'),
                integer_or_blank(card, 5, 'G2'),
            ]
            # :0 if scalar point 1-6 if grid
            components = [components_or_blank(card, 3, 'C1', '0'),
                          components_or_blank(card, 6, 'C2', '0')]
            enforced = [double_or_blank(card, 4, 'D1', 0.0),
                        double_or_blank(card, 7, 'D2', 0.0)]
        return SPC(conid, nodes, components, enforced, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds an SPC card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        spc_id = data[0]
        nodes = [data[1]]
        components_str = sorted(str(data[2]))
        components = int(''.join(components_str))

        assert spc_id > 0, data
        for i, nid in enumerate(nodes):
            assert nodes[0] > 0, f'nodes={nodes} nodes[{i}]={nid}; data={data}'
        assert 0 <= components <= 123456, data
        enforced = [data[3]]
        assert spc_id > 0, data
        assert nodes[0] > 0, data
        components_str = str(components)
        assert len(components_str) <= 6, data
        components = [components_str]
        # if components[0] == 0:
        #     components[0] = 0
        # if components[0] == 16:
        #     components[0] = '16'
        # else:
        #     raise RuntimeError('SPC; components=%s data=%s' % (components, data))
        # assert 0 < components[0] > 1000, data
        return SPC(spc_id, nodes, components, enforced, comment=comment)

    @property
    def constraints(self):
        return self.components

    @constraints.setter
    def constraints(self, constraints):
        self.components = constraints

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        msg = f', which is required by SPC={self.conid:d}'
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by SPC={self.conid:d}'
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors, debug: bool=True):
        nids2 = []
        msg = ', which is required by SPC=%s' % self.conid
        for nid in self.node_ids:
            try:
                nid2 = model.Node(nid, msg=msg)
            except KeyError:
                if debug:
                    msg = 'Couldnt find nid=%i, which is required by SPC=%s' % (
                        nid, self.conid)
                    model.log.warning(msg)
                continue
            nids2.append(nid2)
        self.nodes_ref = nids2

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def raw_fields(self):
        fields = ['SPC', self.conid]
        for (node_id, constraint, enforced) in zip(self.node_ids, self.components,
                                                   self.enforced):
            fields += [node_id, constraint, enforced]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        try:
            out = self.comment + print_card_8(card)
        except RuntimeError:
            out = self.comment + print_card_16(card)
        return out


class GMSPC(Constraint):
    type = 'GMSPC'

    @classmethod
    def _init_from_empty(cls):
        conid = 1
        component = 2
        entity = 3
        entity_id = 4
        return GMSPC(conid, component, entity, entity_id)

    def __init__(self, conid, component, entity, entity_id, comment: str=''):
        Constraint.__init__(self)
        if comment:
            self.comment = comment
        self.conid = conid
        self.component = component
        self.entity = entity
        self.entity_id = entity_id

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a GMSPC card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        conid = integer(card, 1, 'sid')
        component = parse_components(card, 2, 'components')
        entity = string(card, 3, 'entity')
        entity_id = integer(card, 4, 'entity_id')
        return GMSPC(conid, component, entity, entity_id, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        raise NotImplementedError('GMSPC')

    def cross_reference(self, model: BDF) -> None:
        """TODO: xref"""
        #msg = ', which is required by GMSPC=%s' % (self.conid)
        pass

    def safe_cross_reference(self, model, xref_errors):
        """TODO: xref"""
        #msg = ', which is required by GMSPC=%s' % (self.conid)
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        fields = ['GMSPC', self.conid, self.component, self.entity, self.entity_id]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class SPCAX(Constraint):
    """
    Defines a set of single-point constraints or enforced displacements
    for conical shell coordinates.

    +-------+-----+-----+-----+----+-----+
    |   1   |  2  |  3  |  4  |  5 |  6  |
    +=======+=====+=====+=====+====+=====+
    | SPCAX | SID | RID | HID |  C |  D  |
    +-------+-----+-----+-----+----+-----+
    | SPCAX |  2  |  3  |  4  | 13 | 6.0 |
    +-------+-----+-----+-----+----+-----+
    """
    type = 'SPCAX'

    @classmethod
    def _init_from_empty(cls):
        conid = 1
        ringax = 2
        hid = 3
        component = 4
        enforced = 0.
        return SPCAX(conid, ringax, hid, component, enforced)

    def __init__(self, conid, ringax, hid, component, enforced, comment: str=''):
        """
        Creates an SPCAX card

        Parameters
        ----------
        conid : int
            Identification number of a single-point constraint set.
        ringax : int
            Ring identification number. See RINGAX entry.
        hid : int
            Harmonic identification number. (Integer >= 0)
        component : int
            Component identification number. (Any unique combination of the
            Integers 1 through 6.)
        enforced : float
            Enforced displacement value
        """
        # defines everything :) at least until cross-referencing methods are
        # implemented
        Constraint.__init__(self)
        if comment:
            self.comment = comment
        #: Identification number of a single-point constraint set.
        self.conid = conid
        #: Ring identification number. See RINGAX entry.
        self.ringax = ringax
        #: Harmonic identification number. (Integer >= 0)
        self.hid = hid
        #: Component identification number. (Any unique combination of the
        #: Integers 1 through 6.)
        self.component = component
        #: Enforced displacement value
        self.enforced = enforced

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SPCAX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        conid = integer(card, 1, 'conid')
        ringax = integer(card, 2, 'ringax')
        hid = integer(card, 3, 'hid')
        component = parse_components(card, 4, 'component')
        enforced = double(card, 5, 'enforced')
        return SPCAX(conid, ringax, hid, component, enforced, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment: str=''):
        #msg = '%s has not implemented data parsing' % cls.type
        #raise NotImplementedError(msg)

    def cross_reference(self, model: BDF) -> None:
        pass
        #msg = ', which is required by SPCAX=%s' % (self.conid)
        #self.ringax = model.ring[self.ringax]
        #self.hid = model.harmonic[self.hid]

    def safe_cross_reference(self, model, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        fields = ['SPCAX', self.conid, self.ringax, self.hid, self.component, self.enforced]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class SPC1(Constraint):
    """
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    |  1   |  2  |   3   |    4   |   5    |    6   |    7   |    8   |  9 |
    +======+=====+=======+========+========+========+========+========+====+
    | SPC1 | SID |   C   |   G1   |   G2   |   G3   |   G4   |   G5   | G6 |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    |      |  G7 |   G8  |   G9   |  etc.  |        |        |        |    |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    | SPC1 |  3  |  246  | 209075 | 209096 | 209512 | 209513 | 209516 |    |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    | SPC1 |  3  |   2   |   1    |   3    |   10   |   9    |   6    | 5  |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    |      |  2  |   8   |        |        |        |        |        |    |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    | SPC1 | SID |   C   |   G1   |  THRU  |   G2   |        |        |    |
    +------+-----+-------+--------+--------+--------+--------+--------+----+
    | SPC1 | 313 | 12456 |    6   |  THRU  |   32   |        |        |    |
    +------+-----+-------+--------+--------+--------+--------+--------+----+

    """
    type = 'SPC1'
    _properties = ['node_ids']  # 'constraints',

    @classmethod
    def _init_from_empty(cls):
        conid = 1
        components = '1'
        nodes = [1]
        return SPC1(conid, components, nodes)

    def __init__(self, conid: int, components: str, nodes: list[int], comment: str=''):
        """
        Creates an SPC1 card, which defines the degree of freedoms to be
        constrained to a value of 0.0

        Parameters
        ----------
        conid : int
            constraint id
        components : str
            the degree of freedoms to constrain (e.g., '1', '123')
        nodes : list[int]
            GRID/SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        Constraint.__init__(self)
        if comment:
            self.comment = comment
        if isinstance(nodes, integer_types):
            nodes = [nodes]
        if isinstance(components, int):
            components = str(components)
        self.conid = conid
        self.components = components
        self.nodes = expand_thru(nodes)
        self.nodes.sort()
        self.nodes_ref = None

    def validate(self):
        assert isinstance(self.nodes, list), 'nodes=%s\n%s' % (self.nodes, str(self))
        assert isinstance(self.components, str), 'components=%s\n%s' % (self.components, str(self))
        assert len(self.nodes) > 0, self.get_stats()

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SPC1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        conid = integer(card, 1, 'conid')
        components = components_or_blank(card, 2, 'components', '0')  # 246 = y; dx, dz dir
        #nodes = [node for node in card.fields(3) if node is not None]
        nodes = card.fields(3)
        return SPC1(conid, components, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds an SPC1 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        conid = data[0]
        components = str(data[1])
        nodes = data[2]
        if nodes[-1] == -1:
            nodes = nodes[:-1]
        assert conid > 0, data
        assert len(nodes) > 0, data
        for nid in nodes:
            assert nid > 0, data
        return SPC1(conid, components, nodes, comment=comment)

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        msg = ', which is required by SPC1; conid=%s' % self.conid
        return self._node_ids(self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SPC1; conid=%s' % self.conid
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)

    def safe_cross_reference(self, model: BDF, debug=True):
        nids2 = []
        missing_nids = []
        for nid in self.node_ids:
            try:
                nid2 = model.Node(nid)
            except KeyError:
                missing_nids.append(str(nid))
                continue
            nids2.append(nid2)
        if missing_nids and debug:
            model.log.warning("Couldn't find nids=[%s], which is required by SPC1=%s" % (
                ', '.join(missing_nids), self.conid))
        self.nodes_ref = nids2

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def raw_fields(self):
        fields = ['SPC1', self.conid, self.components] + self.node_ids
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class SPCOFF(Constraint):
    """
    +-----+----+----+----+----+----+----+----+----+
    |  1  |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 |
    +=====+====+====+====+====+====+====+====+====+
    | SPC | G1 | C1 | G2 | C2 | G3 | C3 | G4 | C4 |
    +-----+----+----+----+----+----+----+----+----+
    | SPC | 32 | 3  |  5 |    |    |    |    |    |
    +-----+----+----+----+----+----+----+----+----+

    """
    type = 'SPCOFF'

    @classmethod
    def _init_from_empty(cls):
        nodes = [1, 2]
        components = ['1', '2']
        return SPCOFF(nodes, components)

    def __init__(self, nodes, components, comment: str=''):
        Constraint.__init__(self)
        if comment:
            self.comment = comment
        self.nodes = nodes
        self.components = components
        self.nodes_ref = None

    def validate(self):
        assert isinstance(self.nodes, list), self.nodes
        assert isinstance(self.components, list), self.components
        assert len(self.nodes) == len(self.components), 'len(self.nodes)=%s len(self.components)=%s' % (len(self.nodes), len(self.components))
        for nid, comp in zip(self.nodes, self.components):
            assert isinstance(nid, integer_types), self.nodes
            assert isinstance(comp, str), self.components

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SPCOFF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nodes = []
        components = []
        nfields = len(card) - 1
        nconstraints = nfields // 2
        if nfields % 2 == 1:
            nconstraints += 1
        for counter in range(nconstraints):
            igrid = counter + 1
            ifield = counter * 2 + 1
            node = integer(card, ifield, 'G%i' % igrid)
            component = components_or_blank(card, ifield+1, 'C%i' % igrid, '0')
            nodes.append(node)
            components.append(component)
        assert len(card) > 1, f'len(SPCOFF card) = {len(card):d}\ncard={card}'
        return SPCOFF(nodes, components, comment=comment)

    @property
    def constraints(self):
        return self.components

    @constraints.setter
    def constraints(self, constraints):
        self.components = constraints

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        msg = ', which is required by SPCOFF'
        return self._node_ids(nodes=self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SPCOFF'
        self.nodes_ref = model.EmptyNodes(self.nodes, msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors, debug: bool=True):
        nids2 = []
        missing_nids = []
        for nid in self.node_ids:
            try:
                nid2 = model.Node(nid)
            except KeyError:
                missing_nids.append(str(nid))
                continue
            nids2.append(nid2)

        if missing_nids and debug:
            model.log.warning("Couldn't find nids=[%s], which is required by SPCOFF" % (
                ', '.join(missing_nids)))
        self.nodes_ref = nids2

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def raw_fields(self):
        fields = ['SPCOFF']
        for (gid, constraint) in zip(self.node_ids, self.components):
            fields += [gid, constraint]
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class SPCOFF1(Constraint):
    type = 'SPCOFF1'

    @classmethod
    def _init_from_empty(cls):
        components = '1'
        nodes = [1, 2]
        return SPCOFF1(components, nodes)

    def __init__(self, components, nodes, comment: str=''):
        Constraint.__init__(self)
        if comment:
            self.comment = comment
        self.components = components
        self.nodes = expand_thru(nodes)
        self.nodes.sort()
        self.nodes_ref = None

    def validate(self):
        assert isinstance(self.nodes, list), 'nodes=%s\n%s' % (self.nodes, str(self))
        assert isinstance(self.components, str), 'components=%s\n%s' % (self.components, str(self))

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SPCOFF1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        components = parse_components(card, 1, 'components')  # 246 = y; dx, dz dir
        nodes = card.fields(2)
        assert len(card) > 2, f'len(SPCOFF1 card) = {len(card):d}\ncard={card}'
        return cls(components, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds an SPCOFF1 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        components = str(data[0])
        nodes = data[1]
        if nodes[-1] == -1:
            nodes = nodes[:-1]
        for nid in nodes:
            assert nid > 0, data
        return SPCOFF1(components, nodes, comment=comment)

    @property
    def constraints(self):
        return self.components

    @constraints.setter
    def constraints(self, constraints):
        self.components = constraints

    @property
    def node_ids(self):
        if self.nodes_ref is None:
            return self.nodes
        msg = ', which is required by SPCOFF1'
        return self._node_ids(self.nodes_ref, allow_empty_nodes=True, msg=msg)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SPCOFF1'
        self.nodes_ref = model.EmptyNodes(self.node_ids, msg=msg)

    def safe_cross_reference(self, model: BDF, debug=True):
        nids2 = []
        missing_nids = []
        for nid in self.node_ids:
            try:
                nid2 = model.Node(nid)
            except KeyError:
                missing_nids.append(str(nid))
                continue
            nids2.append(nid2)

        if missing_nids and debug:
            model.log.warning("Couldn't find nids=[%s], which is required by SPCOFF1" % (
                ', '.join(missing_nids)))
        self.nodes_ref = nids2

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.nodes = self.node_ids
        self.nodes_ref = None

    def raw_fields(self):
        fields = ['SPCOFF1', self.components] + self.node_ids
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_16(card)


class ConstraintAdd(Constraint):
    """common class for SPCADD, MPCADD"""
    def __init__(self):
        Constraint.__init__(self)
        self.sets_ref = None


class SPCADD(ConstraintAdd):
    """
    Defines a single-point constraint set as a union of single-point constraint
    sets defined on SPC or SPC1 entries.

    +--------+----+----+-----+
    |    1   | 2  |  3 |  4  |
    +========+====+====+=====+
    | SPCADD | 2  |  1 |  3  |
    +--------+----+----+-----+
    """
    type = 'SPCADD'
    _properties = ['ids', 'spc_ids']

    @classmethod
    def _init_from_empty(cls):
        conid = 1
        sets = [1, 2]
        return SPCADD(conid, sets)

    def __init__(self, conid, sets, comment: str=''):
        ConstraintAdd.__init__(self)
        if comment:
            self.comment = comment
        self.conid = conid

        sets_clean = remove_blanks_capitalize(sets)
        if len(sets) != len(sets_clean):
            warnings.warn(f'NSMADD={conid} set length is changed\n'
                          f'old; n={len(sets)}: {set_list_print(sets)}\n'
                          f'new; n={len(sets_clean)}: {set_list_print(sets_clean)}')
        self.sets = expand_thru(sets_clean)
        self.sets.sort()

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a SPCADD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        conid = integer(card, 1, 'conid')
        sets = card.fields(2)
        return SPCADD(conid, sets, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds an SPCADD card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        conid = data[0]
        sets = data[1:].tolist()
        return SPCADD(conid, sets, comment=comment)

    @property
    def spc_ids(self):
        if self.sets_ref is None:
            return self.sets
        spc_ids = []
        for spc in self.sets_ref:
            if isinstance(spc, integer_types):
                spc_ids.append(spc)
            elif isinstance(spc, list):
                spc_ids.append(spc[0].conid)
            else:
                raise TypeError('type=%s; spc=\n%s' % (type(spc), spc))
        return spc_ids

    @property
    def ids(self):
        return self.spc_ids

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by SPCADD=%s' % self.conid
        self.sets_ref = []
        for spc_id in self.sets:
            self.sets_ref.append(model.SPC(spc_id, consider_spcadd=False, msg=msg))

    def safe_cross_reference(self, model: BDF, xref_errors, debug=True):
        self.sets_ref = []
        msg = ', which is required by SPCADD=%s' % self.conid
        for spc_id in self.sets:
            try:
                spc = model.SPC(spc_id, consider_spcadd=False, msg=msg)
            except KeyError:
                if debug:
                    msg = 'Couldnt find SPC=%i, which is required by SPCADD=%s' % (
                        spc_id, self.conid)
                    model.log.warning(msg)
                continue
            self.sets_ref.append(spc)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.sets = self.spc_ids
        self.sets_ref = []

    def raw_fields(self):
        fields = ['SPCADD', self.conid] + self.spc_ids
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_16(card)


class MPCADD(ConstraintAdd):
    r"""
    Defines a multipoint constraint equation of the form
    :math:`\Sigma_j A_j u_j =0` where :math:`u_j` represents
    degree-of-freedom :math:`C_j` at grid or scalar point :math:`G_j`.

    +--------+----+----+-----+
    |    1   | 2  |  3 |  4  |
    +========+====+====+=====+
    | MPCADD | 2  |  1 |  3  |
    +--------+----+----+-----+
    """
    type = 'MPCADD'
    _properties = ['ids', 'mpc_ids']

    @classmethod
    def _init_from_empty(cls):
        conid = 1
        sets = [1, 2]
        return MPCADD(conid, sets)

    def __init__(self, conid: int, sets: list[str | int], comment: str=''):
        ConstraintAdd.__init__(self)
        if comment:
            self.comment = comment
        self.conid = conid

        sets_clean = remove_blanks_capitalize(sets)
        if len(sets) != len(sets_clean):
            warnings.warn(f'NSMADD={conid} set length is changed\n'
                          f'old; n={len(sets)}: {set_list_print(sets)}\n'
                          f'new; n={len(sets_clean)}: {set_list_print(sets_clean)}')
        self.sets = expand_thru(sets_clean)
        self.sets.sort()

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a MPCADD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        conid = integer(card, 1, 'conid')
        sets = card.fields(2)
        return MPCADD(conid, sets, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds an MPCADD card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        conid = data[0]
        sets = data[1:].tolist()
        return MPCADD(conid, sets, comment=comment)

    @property
    def mpc_ids(self):
        if self.sets_ref is None:
            return self.sets
        mpc_ids = []
        for mpc in self.sets_ref:
            if isinstance(mpc, integer_types):
                mpc_ids.append(mpc)
            else:
                mpc_ids.append(mpc[0].conid)
        return mpc_ids

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MPCADD=%s' % self.conid
        self.sets_ref = []
        for mpc_id in self.sets:
            self.sets_ref.append(model.MPC(mpc_id, consider_mpcadd=False, msg=msg))

    def safe_cross_reference(self, model: BDF, debug=True):
        self.sets_ref = []
        msg = ', which is required by MPCADD=%s' % self.conid
        for mpc_id in self.sets:
            try:
                mpc = model.MPC(mpc_id, consider_mpcadd=False, msg=msg)
            except KeyError:
                if debug:
                    msg = 'Couldnt find MPC=%i, which is required by MPCADD=%s' % (
                        mpc_id, self.conid)
                    model.log.warning(msg)
                continue
            self.sets_ref.append(mpc)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.sets = self.mpc_ids
        self.sets_ref = []

    @property
    def ids(self):
        return self.mpc_ids

    def raw_fields(self):
        fields = ['MPCADD', self.conid] + self.mpc_ids
        return fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_16(card)
