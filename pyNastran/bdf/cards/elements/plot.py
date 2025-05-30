# pylint: disable=C0103
"""
All plot elements are defined in this file.  This includes:

 * PLOTEL
 * PLOTEL3
 * PLOTEL4
 * PLOTEL6
 * PLOTEL8
 * PLOTTET
 * PLOTPEN
 * PLOTPYR
 * PLOTHEX

All plot elements are Element objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import (
    BaseCard, _node_ids)
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank,
)
from pyNastran.bdf.field_writer_8 import print_card_8
# from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.cards.base_card import BDFCard


class PLOTEL(BaseCard):
    """
    Defines a 1D dummy element used for plotting.

    This element is not used in the model during any of the solution
    phases of a problem. It is used to simplify plotting of
    structures with large numbers of colinear grid points, where the
    plotting of each grid point along with the elements connecting
    them would result in a confusing plot.

    +--------+-----+-----+-----+
    |   1    |  2  |  3  |  4  |
    +========+=====+=====+=====+
    | PLOTEL | EID | G1  | G2  |
    +--------+-----+-----+-----+

    """
    type = 'PLOTEL'
    _field_map = {
        1: 'eid', 3: 'g1', 4: 'g2',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nodes = [1, 2]
        return PLOTEL(eid, nodes, comment='')

    @classmethod
    def export_to_hdf5(cls, h5_file, model: BDF, encoding: str):
        export_to_hdf5(h5_file, cls.type, model, encoding)

    def __init__(self, eid: int, nodes: list[int], comment: str=''):
        """
        Adds a PLOTEL card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : list[int, int]
            Unique GRID point IDs

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.nodes = nodes
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, icard: int, comment: str=''):
        """
        Adds a PLOTEL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard: int
            allows for two PLOTELs on a single line
        comment : str; default=''
            a comment for the card

        """
        offset = icard * 4
        eid = integer(card, 1+offset, 'eid')
        nodes = [
            integer(card, 2+offset, 'g1'),
            integer(card, 3+offset, 'g2'),
        ]
        #assert len(card) <= 4, f'len(PLOTEL card) = {len(card):d}\ncard={card}'
        assert len(card) <= 8, f'len(PLOTEL card) = {len(card):d}\ncard={card}'
        return PLOTEL(eid, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a PLOTEL card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        nodes = [data[1], data[2]]
        return PLOTEL(eid, nodes, comment=comment)

    def _verify(self, xref):
        pass

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PLOTEL eid=%s' % self.eid
        self.nodes_ref = [
            model.Node(self.nodes[0], msg=msg),
            model.Node(self.nodes[1], msg=msg),
        ]

    def safe_cross_reference(self, model: BDF, xref_errors):
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
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self) -> list[int]:
        if self.nodes_ref is None:
            return self.nodes
        node_idsi = self.nodes_ref
        n1, n2 = node_idsi
        nodes = [n1, n2]
        if not isinstance(n1, integer_types):
            nodes[0] = n1.Nid()
        if not isinstance(n2, integer_types):
            nodes[1] = n2.Nid()
        return nodes

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    def raw_fields(self):
        list_fields = ['PLOTEL', self.eid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        msg = 'PLOTEL  %8i%8i%8i\n' % (self.eid, nodes[0], nodes[1])
        return self.comment + msg


class PLOTEL3(BaseCard):
    """
    Defines a 2D dummy element used for plotting.
    """
    type = 'PLOTEL3'
    _field_map = {
        1: 'eid', 3: 'g1', 4: 'g2', 5: 'g3',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nodes = [1, 2, 3]
        return PLOTEL3(eid, nodes, comment='')

    @classmethod
    def export_to_hdf5(cls, h5_file, model: BDF, encoding: str):
        """exports the elements in a vectorized way"""
        export_to_hdf5(h5_file, cls.type, model, encoding)

    def __init__(self, eid: int, nodes: list[int], comment: str=''):
        """
        Adds a PLOTEL3 card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : list[int, int]
            Unique GRID point IDs

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.nodes = nodes
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PLOTEL3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        nodes = [
            integer(card, 2, 'g1'),
            integer(card, 3, 'g2'),
            integer(card, 4, 'g3'),
        ]
        assert len(card) <= 5, f'len(PLOTEL3 card) = {len(card):d}\ncard={card}'
        return PLOTEL3(eid, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a PLOTEL3 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        nodes = list(data[1:])
        return PLOTEL3(eid, nodes, comment=comment)

    def _verify(self, xref):
        pass

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid:d}'
        self.nodes_ref = [
            model.Node(self.nodes[0], msg=msg),
            model.Node(self.nodes[1], msg=msg),
            model.Node(self.nodes[2], msg=msg),
        ]

    def safe_cross_reference(self, model: BDF, xref_errors):
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
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self) -> list[int]:
        if self.nodes_ref is None:
            return self.nodes
        nodes = [nid_ref.nid for nid_ref in self.nodes_ref]
        return nodes

    # def get_edge_ids(self):
    #     return [tuple(sorted(self.node_ids))]

    def raw_fields(self):
        list_fields = ['PLOTEL3', self.eid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        fields = ['PLOTEL3', self.eid] + nodes
        msg = print_card_8(fields)
        return self.comment + msg


class PLOTEL4(BaseCard):
    """
    Defines a 2D dummy element used for plotting.
    """
    type = 'PLOTEL4'
    _field_map = {
        1: 'eid', 3: 'g1', 4: 'g2', 5: 'g3', 6: 'g4',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nodes = [1, 2, 3, 4]
        return PLOTEL4(eid, nodes, comment='')

    @classmethod
    def export_to_hdf5(cls, h5_file, model, encoding):
        """exports the elements in a vectorized way"""
        export_to_hdf5(h5_file, cls.type, model, encoding)

    def __init__(self, eid: int, nodes: list[int], comment: str=''):
        """
        Adds a PLOTEL4 card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : list[int, int, int, int]
            Unique GRID point IDs

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.nodes = nodes
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PLOTEL4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        nodes = [
            integer(card, 2, 'g1'),
            integer(card, 3, 'g2'),
            integer(card, 4, 'g3'),
            integer(card, 5, 'g4'),
        ]
        assert len(card) <= 6, f'len(PLOTEL4 card) = {len(card):d}\ncard={card}'
        return PLOTEL4(eid, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a PLOTEL4 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        nodes = list(data[1:])
        return PLOTEL4(eid, nodes, comment=comment)

    def _verify(self, xref):
        pass

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid:d}'
        self.nodes_ref = [
            model.Node(self.nodes[0], msg=msg),
            model.Node(self.nodes[1], msg=msg),
            model.Node(self.nodes[2], msg=msg),
            model.Node(self.nodes[3], msg=msg),
        ]

    def safe_cross_reference(self, model: BDF, xref_errors):
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
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self) -> list[int]:
        if self.nodes_ref is None:
            return self.nodes
        nodes = [nid_ref.nid for nid_ref in self.nodes_ref]
        return nodes

    # def get_edge_ids(self) -> list[tuple[int, int]]:
    #     return [tuple(sorted(self.node_ids))]

    def raw_fields(self):
        list_fields = ['PLOTEL4', self.eid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        fields = ['PLOTEL4', self.eid] + nodes
        msg = print_card_8(fields)
        return self.comment + msg


class PLOTEL6(BaseCard):
    """
    Defines a 2D dummy element used for plotting.
    """
    type = 'PLOTEL6'
    _field_map = {
        1: 'eid',
        3: 'g1', 4: 'g2', 5: 'g3',
        6: 'g4', 7: 'g5', 8: 'g6',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nodes = [1, 2, 3, 4, 5, 6]
        return PLOTEL6(eid, nodes, comment='')

    @classmethod
    def export_to_hdf5(cls, h5_file, model: BDF, encoding: str):
        """exports the elements in a vectorized way"""
        export_to_hdf5(h5_file, cls.type, model, encoding)

    def __init__(self, eid: int, nodes: list[int], comment: str=''):
        """
        Adds a PLOTEL6 card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : list[int, int]
            Unique GRID point IDs

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.nodes = nodes
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):

        """
        Adds a PLOTEL6 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        nodes = [
            integer(card, 2, 'g1'),
            integer(card, 3, 'g2'),
            integer(card, 4, 'g3'),

            integer(card, 5, 'g4'),
            integer(card, 6, 'g5'),
            integer(card, 7, 'g6'),
        ]
        assert len(card) <= 8, f'len(PLOTEL6 card) = {len(card):d}\ncard={card}'
        return PLOTEL6(eid, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a PLOTEL6 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        nodes = list(data[1:])
        return PLOTEL6(eid, nodes, comment=comment)

    def _verify(self, xref):
        pass

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid:d}'
        self.nodes_ref = [
            model.Node(self.nodes[0], msg=msg),
            model.Node(self.nodes[1], msg=msg),
            model.Node(self.nodes[2], msg=msg),
            model.Node(self.nodes[3], msg=msg),
            model.Node(self.nodes[4], msg=msg),
            model.Node(self.nodes[5], msg=msg),
        ]

    def safe_cross_reference(self, model: BDF, xref_errors):
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
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self) -> list[int]:
        if self.nodes_ref is None:
            return self.nodes
        nodes = [0 if nid_ref is None else nid_ref.nid
                 for nid_ref in self.nodes_ref]
        return nodes

    # def get_edge_ids(self):
    #     return [tuple(sorted(self.node_ids))]

    def raw_fields(self):
        list_fields = ['PLOTEL6', self.eid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        fields = ['PLOTEL6', self.eid] + nodes
        msg = print_card_8(fields)
        return self.comment + msg


class PLOTEL8(BaseCard):
    """
    Defines a 2D dummy element used for plotting.
    """
    type = 'PLOTEL8'
    _field_map = {
        1: 'eid',
        3: 'g1', 4: 'g2', 5: 'g3', 6: 'g4',
        7: 'g5', 8: 'g6', 9: 'g7', 10: 'g8',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nodes = [1, 2, 3, 4, 5, 6, 7, 8]
        return PLOTEL8(eid, nodes, comment='')

    @classmethod
    def export_to_hdf5(cls, h5_file, model: BDF, encoding):
        """exports the elements in a vectorized way"""
        export_to_hdf5(h5_file, cls.type, model, encoding)

    def __init__(self, eid: int, nodes: list[int], comment: str=''):
        """
        Adds a PLOTEL8 card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : list[int, int, int, int, int, int, int, int]
            Unique GRID point IDs

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.nodes = nodes
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PLOTEL8 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        nodes = [
            integer(card, 2, 'g1'),
            integer(card, 3, 'g2'),
            integer(card, 4, 'g3'),
            integer(card, 5, 'g4'),

            integer(card, 6, 'g5'),
            integer(card, 7, 'g6'),
            integer(card, 8, 'g7'),
            integer(card, 9, 'g8'),
        ]
        assert len(card) <= 10, f'len(PLOTEL8 card) = {len(card):d}\ncard={card}'
        return PLOTEL8(eid, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a PLOTEL8 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        nodes = list(data[1:])
        return PLOTEL8(eid, nodes, comment=comment)

    def _verify(self, xref):
        pass

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid:d}'
        self.nodes_ref = [
            model.Node(self.nodes[0], msg=msg),
            model.Node(self.nodes[1], msg=msg),
            model.Node(self.nodes[2], msg=msg),
            model.Node(self.nodes[3], msg=msg),

            model.Node(self.nodes[4], msg=msg),
            model.Node(self.nodes[5], msg=msg),
            model.Node(self.nodes[6], msg=msg),
            model.Node(self.nodes[7], msg=msg),
        ]

    def safe_cross_reference(self, model: BDF, xref_errors):
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
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self) -> list[int]:
        if self.nodes_ref is None:
            return self.nodes
        nodes = [0 if nid_ref is None else nid_ref.nid
                 for nid_ref in self.nodes_ref]
        return nodes

    # def get_edge_ids(self) -> list[tuple[int, int]]:
    #     return [tuple(sorted(self.node_ids))]

    def raw_fields(self):
        list_fields = ['PLOTEL8', self.eid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        fields = ['PLOTEL8', self.eid] + nodes
        msg = print_card_8(fields)
        return self.comment + msg


class PLOTTET(BaseCard):
    """
    Defines a 3D dummy element used for plotting.
    """
    type = 'PLOTTET'
    _field_map = {
        1: 'eid',
        3: 'g1', 4: 'g2', 5: 'g3', 6: 'g4', 7: 'g5',
        8: 'g6', 9: 'g7', 10: 'g8', 11: 'g9', 12: 'g10',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        return PLOTTET(eid, nodes, comment='')

    @classmethod
    def export_to_hdf5(cls, h5_file, model: BDF, encoding):
        """exports the elements in a vectorized way"""
        export_to_hdf5(h5_file, cls.type, model, encoding)

    def __init__(self, eid: int, nodes: list[int], comment: str=''):
        """
        Adds a PLOTTET card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : list[int, int, int, int, int,
                     int, int, int, int, int]
            Unique GRID point IDs

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        nnodes = len(nodes)
        if nnodes < 10:
            nodes.extend([0] * (10-nnodes))
        self.nodes = nodes
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PLOTTET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        nodes = [
            integer(card, 2, 'g1'),
            integer(card, 3, 'g2'),
            integer(card, 4, 'g3'),
            integer(card, 5, 'g4'),

            integer_or_blank(card, 6, 'g5', default=0),
            integer_or_blank(card, 7, 'g6', default=0),
            integer_or_blank(card, 8, 'g7', default=0),
            integer_or_blank(card, 9, 'g8', default=0),

            integer_or_blank(card, 10, 'g9', default=0),
            integer_or_blank(card, 11, 'g10', default=0),
        ]
        assert len(card) <= 12, f'len(PLOTTET card) = {len(card):d}\ncard={card}'
        return PLOTTET(eid, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a PLOTTET card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        nodes = list(data[1:])
        return PLOTTET(eid, nodes, comment=comment)

    def _verify(self, xref):
        pass

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid:d}'
        n1, n2, n3, n4, n5, n6, n7, n8, n9, n10 = self.nodes
        self.nodes_ref = [
            model.Node(n1, msg=msg),
            model.Node(n2, msg=msg),
            model.Node(n3, msg=msg),
            model.Node(n4, msg=msg),
            None if n5 == 0 else model.Node(n5, msg=msg),
            None if n6 == 0 else model.Node(n6, msg=msg),
            None if n7 == 0 else model.Node(n6, msg=msg),
            None if n8 == 0 else model.Node(n8, msg=msg),
            None if n9 == 0 else model.Node(n9, msg=msg),
            None if n10 == 0 else model.Node(n10, msg=msg),
        ]

    def safe_cross_reference(self, model: BDF, xref_errors):
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
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self) -> list[int]:
        if self.nodes_ref is None:
            return self.nodes
        nodes = [0 if nid_ref is None else nid_ref.nid
                 for nid_ref in self.nodes_ref]
        return nodes

    # def get_edge_ids(self) -> list[tuple[int, int]]:
    #     return [tuple(sorted(self.node_ids))]

    def raw_fields(self):
        list_fields = ['PLOTTET', self.eid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        fields = ['PLOTTET', self.eid] + nodes
        msg = print_card_8(fields)
        return self.comment + msg


class PLOTPYR(BaseCard):
    """
    Defines a 3D dummy element used for plotting.
    """
    type = 'PLOTPYR'
    _field_map = {
        1: 'eid',
        3: 'g1', 4: 'g2', 5: 'g3', 6: 'g4',
        7: 'g5', 8: 'g6', 9: 'g7', 10: 'g8',
        11: 'g9', 12: 'g10', 13: 'g11', 14: 'g12',
        15: 'g13',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
        return PLOTPYR(eid, nodes, comment='')

    @classmethod
    def export_to_hdf5(cls, h5_file, model: BDF, encoding):
        """exports the elements in a vectorized way"""
        export_to_hdf5(h5_file, cls.type, model, encoding)

    def __init__(self, eid: int, nodes: list[int], comment: str=''):
        """
        Adds a PLOTPYR card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : list[int, int, int, int, int,
                     int, int, int, int, int,
                     int, int, int]
            Unique GRID point IDs

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        nnodes = len(nodes)
        if nnodes < 13:
            nodes.extend([0] * (13-nnodes))
        self.nodes = nodes
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PLOTPYR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        nodes = [
            integer(card, 2, 'g1'),
            integer(card, 3, 'g2'),
            integer(card, 4, 'g3'),
            integer(card, 5, 'g4'),
            integer(card, 6, 'g5'),
            integer_or_blank(card, 7, 'g6', default=0),
            integer_or_blank(card, 8, 'g7', default=0),
            integer_or_blank(card, 9, 'g8', default=0),

            integer_or_blank(card, 10, 'g9', default=0),
            integer_or_blank(card, 11, 'g10', default=0),
            integer_or_blank(card, 12, 'g11', default=0),
            integer_or_blank(card, 13, 'g12', default=0),
            integer_or_blank(card, 14, 'g13', default=0),
        ]
        assert len(card) <= 15, f'len(PLOTPYR card) = {len(card):d}\ncard={card}'
        return PLOTPYR(eid, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a PLOTTET card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        nodes = list(data[1:])
        return PLOTPYR(eid, nodes, comment=comment)

    def _verify(self, xref):
        pass

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid:d}'
        (n1, n2, n3, n4, n5,
         n6, n7, n8, n9, n10,
         n11, n12, n13) = self.nodes
        self.nodes_ref = [
            model.Node(n1, msg=msg),
            model.Node(n2, msg=msg),
            model.Node(n3, msg=msg),
            model.Node(n4, msg=msg),
            model.Node(n5, msg=msg),
            None if n6 == 0 else model.Node(n6, msg=msg),
            None if n7 == 0 else model.Node(n6, msg=msg),
            None if n8 == 0 else model.Node(n8, msg=msg),
            None if n9 == 0 else model.Node(n9, msg=msg),
            None if n10 == 0 else model.Node(n10, msg=msg),
            None if n11 == 0 else model.Node(n11, msg=msg),
            None if n12 == 0 else model.Node(n12, msg=msg),
            None if n13 == 0 else model.Node(n13, msg=msg),
        ]

    def safe_cross_reference(self, model: BDF, xref_errors):
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
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self) -> list[int]:
        if self.nodes_ref is None:
            return self.nodes
        nodes = [0 if nid_ref is None else nid_ref.nid
                 for nid_ref in self.nodes_ref]
        return nodes

    # def get_edge_ids(self) -> list[tuple[int, int]]:
    #     return [tuple(sorted(self.node_ids))]

    def raw_fields(self):
        list_fields = ['PLOTPYR', self.eid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        fields = ['PLOTPYR', self.eid] + nodes
        msg = print_card_8(fields)
        return self.comment + msg


class PLOTPEN(BaseCard):
    """
    Defines a 3D dummy element used for plotting.
    """
    type = 'PLOTPEN'
    _field_map = {
        1: 'eid',
        3: 'g1', 4: 'g2', 5: 'g3', 6: 'g4',
        7: 'g5', 8: 'g6', 9: 'g7', 10: 'g8',
        11: 'g9', 12: 'g10', 13: 'g11', 14: 'g12',
        15: 'g13', 16: 'g14', 17: 'g15',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                 11, 12, 13, 14, 15]
        return PLOTPEN(eid, nodes, comment='')

    @classmethod
    def export_to_hdf5(cls, h5_file, model: BDF, encoding):
        """exports the elements in a vectorized way"""
        export_to_hdf5(h5_file, cls.type, model, encoding)

    def __init__(self, eid: int, nodes: list[int], comment: str=''):
        """
        Adds a PLOTPEN card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : list[int, int, int, int, int,
                     int, int, int, int, int,
                     int, int, int, int, int]
            Unique GRID point IDs

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        nnodes = len(nodes)
        if nnodes < 15:
            nodes.extend([0] * (15-nnodes))
        self.nodes = nodes
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PLOTPEN card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        nodes = [
            integer(card, 2, 'g1'),
            integer(card, 3, 'g2'),
            integer(card, 4, 'g3'),
            integer(card, 5, 'g4'),
            integer(card, 6, 'g5'),
            integer(card, 7, 'g6'),
            integer_or_blank(card, 8, 'g7', default=0),
            integer_or_blank(card, 9, 'g8', default=0),

            integer_or_blank(card, 10, 'g9', default=0),
            integer_or_blank(card, 11, 'g10', default=0),
            integer_or_blank(card, 12, 'g11', default=0),
            integer_or_blank(card, 13, 'g12', default=0),
            integer_or_blank(card, 14, 'g13', default=0),
            integer_or_blank(card, 15, 'g14', default=0),
            integer_or_blank(card, 16, 'g15', default=0),
        ]
        assert len(card) <= 17, f'len(PLOTPEN card) = {len(card):d}\ncard={card}'
        return PLOTPEN(eid, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a PLOTPEN card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        nodes = list(data[1:])
        return PLOTPEN(eid, nodes, comment=comment)

    def _verify(self, xref):
        pass

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid:d}'
        (n1, n2, n3, n4, n5,
         n6, n7, n8, n9, n10,
         n11, n12, n13, n14, n15) = self.nodes
        self.nodes_ref = [
            model.Node(n1, msg=msg),
            model.Node(n2, msg=msg),
            model.Node(n3, msg=msg),
            model.Node(n4, msg=msg),
            model.Node(n5, msg=msg),
            model.Node(n6, msg=msg),
            None if n7 == 0 else model.Node(n6, msg=msg),
            None if n8 == 0 else model.Node(n8, msg=msg),
            None if n9 == 0 else model.Node(n9, msg=msg),
            None if n10 == 0 else model.Node(n10, msg=msg),
            None if n11 == 0 else model.Node(n11, msg=msg),
            None if n12 == 0 else model.Node(n12, msg=msg),
            None if n13 == 0 else model.Node(n13, msg=msg),
            None if n14 == 0 else model.Node(n14, msg=msg),
            None if n15 == 0 else model.Node(n15, msg=msg),
        ]

    def safe_cross_reference(self, model: BDF, xref_errors):
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
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self) -> list[int]:
        if self.nodes_ref is None:
            return self.nodes
        nodes = [0 if nid_ref is None else nid_ref.nid
                 for nid_ref in self.nodes_ref]
        return nodes

    # def get_edge_ids(self) -> list[tuple[int, int]]:
    #     return [tuple(sorted(self.node_ids))]

    def raw_fields(self):
        list_fields = ['PLOTPEN', self.eid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        fields = ['PLOTPEN', self.eid] + nodes
        msg = print_card_8(fields)
        return self.comment + msg


class PLOTHEX(BaseCard):
    """
    Defines a 3D dummy element used for plotting.
    """
    type = 'PLOTHEX'
    _field_map = {
        1: 'eid',
        3: 'g1', 4: 'g2', 5: 'g3', 6: 'g4',
        7: 'g5', 8: 'g6', 9: 'g7', 10: 'g8',
        11: 'g9', 12: 'g10', 13: 'g11', 14: 'g12',
        15: 'g13', 16: 'g14', 17: 'g15',
        18: 'g16', 19: 'g17', 20: 'g18', 21: 'g19',
        22: 'g20',
    }
    _properties = ['node_ids']

    @classmethod
    def _init_from_empty(cls):
        eid = 1
        nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        return PLOTHEX(eid, nodes, comment='')

    @classmethod
    def export_to_hdf5(cls, h5_file, model: BDF, encoding):
        """exports the elements in a vectorized way"""
        export_to_hdf5(h5_file, cls.type, model, encoding)

    def __init__(self, eid: int, nodes: list[int], comment: str=''):
        """
        Adds a PLOTHEX card

        Parameters
        ----------
        eid : int
            Element ID
        nodes : list[int, int, int, int, int,
                     int, int, int, int, int,
                     int, int, int, int, int,
                     int, int, int, int, int]
            Unique GRID point IDs

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        nnodes = len(nodes)
        if nnodes < 20:
            nodes.extend([0] * (20-nnodes))
        self.nodes = nodes
        self.nodes_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PLOTHEX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        nodes = [
            integer(card, 2, 'g1'),
            integer(card, 3, 'g2'),
            integer(card, 4, 'g3'),
            integer(card, 5, 'g4'),
            integer(card, 6, 'g5'),
            integer(card, 7, 'g6'),
            integer(card, 8, 'g7'),
            integer(card, 9, 'g8'),

            integer_or_blank(card, 10, 'g9', default=0),
            integer_or_blank(card, 11, 'g10', default=0),
            integer_or_blank(card, 12, 'g11', default=0),
            integer_or_blank(card, 13, 'g12', default=0),
            integer_or_blank(card, 14, 'g13', default=0),
            integer_or_blank(card, 15, 'g14', default=0),
            integer_or_blank(card, 16, 'g15', default=0),
            integer_or_blank(card, 17, 'g16', default=0),
            integer_or_blank(card, 18, 'g17', default=0),
            integer_or_blank(card, 19, 'g18', default=0),
            integer_or_blank(card, 20, 'g19', default=0),
            integer_or_blank(card, 21, 'g20', default=0),
        ]
        assert len(card) <= 22, f'len(PLOTHEX card) = {len(card):d}\ncard={card}'
        return PLOTHEX(eid, nodes, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a PLOTHEX card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        eid = data[0]
        nodes = list(data[1:])
        return PLOTHEX(eid, nodes, comment=comment)

    def _verify(self, xref):
        pass

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = f', which is required by {self.type} eid={self.eid:d}'
        (n1, n2, n3, n4, n5,
         n6, n7, n8, n9, n10,
         n11, n12, n13, n14, n15,
         n16, n17, n18, n19, n20) = self.nodes
        self.nodes_ref = [
            model.Node(n1, msg=msg), model.Node(n2, msg=msg),
            model.Node(n3, msg=msg), model.Node(n4, msg=msg),
            model.Node(n5, msg=msg), model.Node(n6, msg=msg),
            model.Node(n7, msg=msg), model.Node(n8, msg=msg),
            None if n9 == 0 else model.Node(n9, msg=msg),
            None if n10 == 0 else model.Node(n10, msg=msg),
            None if n11 == 0 else model.Node(n11, msg=msg),
            None if n12 == 0 else model.Node(n12, msg=msg),
            None if n13 == 0 else model.Node(n13, msg=msg),
            None if n14 == 0 else model.Node(n14, msg=msg),
            None if n15 == 0 else model.Node(n15, msg=msg),
            None if n16 == 0 else model.Node(n16, msg=msg),
            None if n17 == 0 else model.Node(n17, msg=msg),
            None if n18 == 0 else model.Node(n18, msg=msg),
            None if n19 == 0 else model.Node(n19, msg=msg),
            None if n20 == 0 else model.Node(n20, msg=msg),
        ]

    def safe_cross_reference(self, model: BDF, xref_errors):
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
        self.nodes = self.node_ids
        self.nodes_ref = None

    @property
    def node_ids(self) -> list[int]:
        if self.nodes_ref is None:
            return self.nodes
        nodes = [0 if nid_ref is None else nid_ref.nid
                 for nid_ref in self.nodes_ref]
        return nodes

    # def get_edge_ids(self) -> list[tuple[int, int]]:
    #     return [tuple(sorted(self.node_ids))]

    def raw_fields(self):
        list_fields = ['PLOTHEX', self.eid] + self.node_ids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        nodes = self.node_ids
        fields = ['PLOTHEX', self.eid] + nodes
        msg = print_card_8(fields)
        return self.comment + msg


def export_to_hdf5(h5_file, plot_type: str,
                   model: BDF, encoding: str):
    """exports the elements in a vectorized way"""
    #comments = []
    nodes = []
    eids = list(model.plotels.keys())
    for eid in eids:
        element = model.plotels[eid]
        if element.type != plot_type:
            continue
        #comments.append(element.comment)
        nodes.append(element.nodes)
    #h5_file.create_dataset('_comment', data=comments)
    h5_file.create_dataset('eid', data=eids)
    h5_file.create_dataset('nodes', data=nodes)

PLOTELs = (PLOTEL | PLOTEL3 | PLOTEL4 |
           PLOTEL6 | PLOTEL8 |
           PLOTTET | PLOTPYR | PLOTPEN | PLOTHEX)
