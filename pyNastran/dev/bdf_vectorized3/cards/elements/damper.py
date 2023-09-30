from __future__ import annotations
from itertools import zip_longest
from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8, print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, integer_or_blank, double_or_blank,
    integer_double_or_blank)
from pyNastran.bdf.cards.elements.bars import set_blank_if_default

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, get_print_card_8_16,
    parse_element_check, parse_property_check)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class CDAMP1(Element):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CDAMP1 | EID | PID | G1 | C1 | G2 | C2 |
    +--------+-----+-----+----+----+----+----+
    """
    def add(self, eid: int, pid: int, nids: list[int], c1: int=0, c2: int=0,
            comment: str='') -> int:
        """
        Creates a CDAMP1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP)
        nids : list[int, int]
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nids, [c1, c2], comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        nids = [integer(card, 3, 'g1'),
                integer_or_blank(card, 5, 'g2', default=0)]

        #: component number
        c1 = integer_or_blank(card, 4, 'c1', default=0)
        c2 = integer_or_blank(card, 6, 'c2', default=0)
        assert len(card) <= 7, f'len(CDAMP1 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, [c1, c2], comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 2), dtype='int32')
        components = np.zeros((ncards, 2), dtype='int32')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, componentsi, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            components[icard, :] = componentsi
        self._save(element_id, property_id, nodes, components)
        self.cards = []

    def _save(self, element_id, property_id, nodes, components):
        nelements = len(element_id)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.components = components
        self.n = nelements

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no damper properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        components_ = array_default_int(self.components, default=0, size=size)
        for eid, pid, nodes, components in zip(element_id, property_id, nodes_, components_):
            list_fields = ['CDAMP1', eid, pid,
                           nodes[0], components[0],
                           nodes[1], components[1]]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pdamp]
                if prop.n > 0]

    #def mass(self) -> np.ndarray:
        #pid = self.property_id
        #mass_per_length = np.zeros(len(pid), dtype='float64')
        #for prop in self.allowed_properties:
            #i = np.searchsorted(prop.property_id, pid)

            ## filter properties that weren't found
            #ibad = (i == len(prop.property_id))
            #i[ibad] = 0

            ## make sure all properties we're setting were found
            #is_pid = (prop.property_id[i] == pid)
            #if len(is_pid) == 0:
                #continue

            ## we're at least using some properties
            #iprop = i[is_pid]
            #mass_per_lengthi = prop.mass_per_length()
            #mass_per_length[is_pid] = mass_per_lengthi[iprop]
        #length = self.length()
        #mass = mass_per_length * length
        #return mass

    #def length(self) -> np.ndarray:
        #length = line_length(self.model, self.nodes)
        #return length

    #def centroid(self) -> np.ndarray:
        #centroid = line_centroid(self.model, self.nodes)
        #return centroid


class CDAMP2(Element):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CDAMP2 | EID |  B  | G1 | C1 | G2 | C2 |
    +--------+-----+-----+----+----+----+----+
    """
    def add(self, eid: int, b: float, nids: list[int],
            c1: int=0, c2: int=0, comment: str='') -> int:
        """
        Creates a CDAMP2 card

        Parameters
        ----------
        eid : int
            element id
        b : float
            damping
        nids : list[int, int]
            SPOINT ids
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, b, nids, (c1, c2), comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        b = double(card, 2, 'b')
        nids = [integer_or_blank(card, 3, 'g1', default=0),
                integer_or_blank(card, 5, 'g2', default=0)]
        c1 = integer_or_blank(card, 4, 'c1', default=0)
        c2 = integer_or_blank(card, 6, 'c2', default=0)
        assert len(card) <= 7, f'len(CDAMP2 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, b, nids, (c1, c2), comment))
        self.n += 1
        return self.n

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        b = np.zeros(ncards, dtype='float64')
        nodes = np.zeros((ncards, 2), dtype='int32')
        components = np.zeros((ncards, 2), dtype='int32')

        for icard, card in enumerate(self.cards):
            eid, bi, nids, componentsi, comment = card
            element_id[icard] = eid
            b[icard] = bi
            nodes[icard, :] = nids
            components[icard, :] = componentsi
        self._save(element_id, nodes, components, b)
        self.cards = []

    def _save(self, element_id, nodes, components, b):
        nelements = len(element_id)
        self.element_id = element_id
        self.nodes = nodes
        self.components = components
        self.b = b
        self.n = nelements

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.nodes))

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        components_ = array_default_int(self.components, default=0, size=size)
        for eid, b, nodes, components in zip(element_id, self.b,
                                                    nodes_, components_):
            list_fields = ['CDAMP2', eid, b,
                           nodes[0], components[0],
                           nodes[1], components[1]]
            bdf_file.write(print_card(list_fields))
        return


class CDAMP3(Element):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CDAMP3 | EID | PID | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    def add(self, eid: int, pid: int, nids: list[int],
            comment: str='') -> int:
        """
        Creates a CDAMP3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP)
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)

        s1 = integer_or_blank(card, 3, 's1', default=0)
        s2 = integer_or_blank(card, 4, 's2', default=0)
        assert len(card) <= 5, f'len(CDAMP3 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, (s1, s2), comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        spoints = np.zeros((ncards, 2), dtype='int32')

        for icard, card in enumerate(self.cards):
            (eid, pid, spointsi, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            spoints[icard, :] = spointsi
        self._save(element_id, property_id, spoints)
        self.cards = []

    def _save(self, element_id, property_id, spoints):
        nelements = len(element_id)
        self.element_id = element_id
        self.property_id = property_id
        self.spoints = spoints
        self.n = nelements

    def geom_check(self, missing: dict[str, np.ndarray]):
        spoint = self.model.spoint
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no damper properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   spoint=(spoint, self.spoints),
                   property_id=(pids, self.property_id))

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        spoints_ = array_default_int(self.spoints, default=0, size=size)
        for eid, pid, spoints in zip(element_id, property_id, spoints_):
            msg = 'CDAMP3  %8s%8s%8s%8s\n' % (eid, pid, spoints[0], spoints[1])
            #list_fields = ['CDAMP3', eid, pid, spoints[0], spoints[1]]
            bdf_file.write(msg)
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pdamp]
                if prop.n > 0]


class CDAMP4(Element):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CDAMP4 | EID |  B  | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.spoints = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, b: float, nids: list[int],
            comment: str='') -> int:
        """
        Creates a CDAMP4 card

        Parameters
        ----------
        eid : int
            element id
        b : float
            damping
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, b, nids, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        b = double(card, 2, 'b')
        s1 = integer_or_blank(card, 3, 's1', default=0)
        s2 = integer_or_blank(card, 4, 's2', default=0)
        self.cards.append((eid, b, [s1, s2], comment))
        self.n += 1
        if card.field(5):
            eid = integer(card, 5, 'eid')
            b = double(card, 6, 'b')
            s1 = integer_or_blank(card, 7, 's1', default=0)
            s2 = integer_or_blank(card, 8, 's2', default=0)
            self.cards.append((eid, b, [s1, s2], comment))
            self.n += 1
        assert len(card) <= 9, f'len(CDAMP4 card) = {len(card):d}\ncard={card}'
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        b = np.zeros(ncards, dtype='float64')
        spoints = np.zeros((ncards, 2), dtype='int32')

        for icard, card in enumerate(self.cards):
            eid, bi, spointsi, comment = card
            element_id[icard] = eid
            b[icard] = bi
            spoints[icard, :] = spointsi
        self._save(element_id, b, spoints)
        self.cards = []

    def _save(self, element_id, b, spoints):
        assert len(self.element_id) == 0
        nelements = len(element_id)
        self.element_id = element_id
        self.b = b
        self.spoints = spoints
        self.n = nelements

    def geom_check(self, missing: dict[str, np.ndarray]):
        spoint = self.model.spoint
        geom_check(self,
                   missing,
                   spoint=(spoint, self.spoints))

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        spoints_ = array_default_int(self.spoints, default=0, size=size)
        for eid, b, spoints in zip(element_id, self.b, spoints_):
            list_fields = ['CDAMP4', eid, b, spoints[0], spoints[1]]
            bdf_file.write(print_card(list_fields))
        return


class CDAMP5(Element):
    """
    Defines a damping element that refers to a material property entry and connection to
    grid or scalar points.

    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CDAMP5 | EID | PID | N1 | N2 |
    +--------+-----+-----+----+----+
    """
    def add(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CDAMP5 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PDAMP5)
        nids : list[int, int]
            GRID/SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 'n1', 0),
                integer_or_blank(card, 4, 'n2', 0)]
        assert len(card) <= 5, f'len(CDAMP5 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 2), dtype='int32')

        for icard, card in enumerate(self.cards):
            eid, pid, nids, comment = card

            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes)
        self.cards = []

    def _save(self, element_id, property_id, nodes):
        assert len(self.element_id) == 0
        nelements = len(element_id)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.n = nelements

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        for eid, pid, nodes in zip(element_id, property_id, nodes_):
            list_fields = ['CDAMP5', eid, pid, nodes[0], nodes[1]]
            bdf_file.write(print_card(list_fields))
        return


class PDAMP(Property):
    """
    +-------+------+-----+------+----+------+----+------+----+
    |   1   |  2   |  3  |   4  | 5  |  6   |  7 |   8  |  9 |
    +=======+======+=====+======+====+======+====+======+====+
    | PDAMP | PID1 | B1  | PID2 | B2 | PID3 | B3 | PID4 | B4 |
    +-------+------+-----+------+----+------+----+------+----+
    | PDAMP |  1   | 2.0 |      |    |      |    |      |    |
    +-------+------+-----+------+----+------+----+------+----+
    """
    def add(self, pid: int, b: float, comment: str='') -> int:
        """
        Creates a PDAMP card

        Parameters
        ----------
        pid : int
            property id
        b : float
            viscous damping
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, b, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> list[int]:
        """adds a PDAMP"""
        ns = [self.n]
        pid = integer(card, 1, 'pid')
        b = double(card, 2, 'b')
        self.cards.append((pid, b, comment))
        self.n += 1

        if card.field(3):
            pid = integer(card, 3, 'pid2')
            b = double(card, 4, 'b2')
            self.cards.append((pid, b, comment))
            ns.append(self.n)
            self.n += 1

        if card.field(5):
            pid = integer(card, 5, 'pid3')
            b = double(card, 6, 'b3')
            self.cards.append((pid, b, comment))
            ns.append(self.n)
            self.n += 1

        if card.field(7):
            pid = integer(card, 7, 'pid4')
            b = double(card, 8, 'b4')
            self.cards.append((pid, b, comment))
            ns.append(self.n)
            self.n += 1

        assert len(card) <= 9, f'len(PDAMP card) = {len(card):d}\ncard={card}'
        return ns

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        b = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, bi, comment) = card
            property_id[icard] = pid
            b[icard] = bi
        self._save(property_id, b)
        self.cards = []

    def _save(self, property_id, b):
        nproperties = len(property_id)
        self.property_id = property_id
        self.b = b
        self.n = nproperties

    def validate(self) -> None:
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_id = array_str(self.property_id, size=size)
        for pid, b in zip(property_id, self.b):
            list_fields = ['PDAMP', pid, b]
            bdf_file.write(print_card(list_fields))
        return


class PDAMP5(Property):
    """
    +--------+------+-----+------+----+------+----+------+----+
    |    1   |  2   |  3  |   4  | 5  |  6   |  7 |   8  |  9 |
    +========+======+=====+======+====+======+====+======+====+
    | PDAMP5 | PID  | MID |  B   |    |      |    |      |    |
    +--------+------+-----+------+----+------+----+------+----+
    | PDAMP5 |  1   | 2   |  2.0 |    |      |    |      |    |
    +--------+------+-----+------+----+------+----+------+----+
    """
    def add(self, pid: int, mid: int, b: float, comment: str='') -> int:
        """
        Creates a PDAMP card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        b : float
            viscous damping
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, mid, b, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> list[int]:
        """adds a PDAMP"""
        ns = [self.n]
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        b = double(card, 3, 'b')
        self.cards.append((pid, mid, b, comment))
        self.n += 1

        assert len(card) <= 3, f'len(PDAMP5 card) = {len(card):d}\ncard={card}'
        return self.n

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        b = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, mid, bi, comment) = card
            property_id[icard] = pid
            material_id[icard] = mid
            b[icard] = bi
        self._save(property_id, material_id, b)
        self.cards = []

    def _save(self, property_id, material_id, b):
        nproperties = len(property_id)
        self.property_id = property_id
        self.material_id = material_id
        self.b = b
        self.n = nproperties

    def validate(self) -> None:
        return

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_id = array_str(self.property_id, size=size)
        material_id = array_str(self.material_id, size=size)
        for pid, mid, b in zip(property_id, material_id, self.b):
            list_fields = ['PDAMP5', pid, mid, b]
            bdf_file.write(print_card(list_fields))
        return


class PDAMPT(Property):
    def add(self, pid: int, tbid: int, comment: str='') -> int:
        """
        Creates a PDAMPT card

        Parameters
        ----------
        pid : int
            property id
        tbid : int
            TABLED1? id
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, tbid, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        tbid = integer_or_blank(card, 2, 'tbid', 0)
        assert len(card) <= 3, f'len(PDAMPT card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, tbid, comment))
        self.n += 1
        return self.n

    def parse_cards(self):
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return

        property_id = np.zeros(ncards, dtype='int32')

        #: Identification number of a TABLEDi entry that defines the
        #: damping force per-unit velocity versus frequency relationship
        table_b = np.zeros(ncards, dtype='int32')

        for icard, (pid, tbid, comment) in enumerate(self.cards):
            property_id[icard] = pid
            table_b[icard] = tbid
        self._save(property_id, table_b)
        self.cards = []

    def _save(self, property_id, table_b) -> None:
        assert len(self.property_id) == 0
        self.property_id = property_id
        self.table_b = table_b

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_id = array_str(self.property_id, size=size)
        table_bs = array_default_int(self.table_b, default=0, size=size)
        for pid, b in zip(property_id, table_bs):
            list_fields = ['PDAMPT', pid, b]
            bdf_file.write(print_card(list_fields))
        return


class CVISC(Element):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  6 |
    +========+=====+=====+====+====+
    | CVISC  | EID | PID | G1 | G2 |
    +--------+-----+-----+----+----+
    """
    def add(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CVISC card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PVISC)
        nids : list[int, int]
            GRID ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        nids = [integer(card, 3, 'g1'),
                integer_or_blank(card, 5, 'g2', default=0)]

        assert len(card) <= 6, f'len(CVISC card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 2), dtype='int32')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes)
        self.cards = []

    def _save(self, element_id, property_id, nodes):
        nelements = len(element_id)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.n = nelements

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no visc properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        for eid, pid, nodes in zip(element_id, property_id, nodes_):
            list_fields = ['CVISC', eid, pid, nodes[0], nodes[1]]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pvisc]
                if prop.n > 0]

    #def mass(self) -> np.ndarray:
        #pid = self.property_id
        #mass_per_length = np.zeros(len(pid), dtype='float64')
        #for prop in self.allowed_properties:
            #i = np.searchsorted(prop.property_id, pid)

            ## filter properties that weren't found
            #ibad = (i == len(prop.property_id))
            #i[ibad] = 0

            ## make sure all properties we're setting were found
            #is_pid = (prop.property_id[i] == pid)
            #if len(is_pid) == 0:
                #continue

            ## we're at least using some properties
            #iprop = i[is_pid]
            #mass_per_lengthi = prop.mass_per_length()
            #mass_per_length[is_pid] = mass_per_lengthi[iprop]
        #length = self.length()
        #mass = mass_per_length * length
        #return mass

    #def length(self) -> np.ndarray:
        #length = line_length(self.model, self.nodes)
        #return length

    #def centroid(self) -> np.ndarray:
        #centroid = line_centroid(self.model, self.nodes)
        #return centroid


class PVISC(Property):
    """
    Viscous Damping Element Property
    Defines properties of a one-dimensional viscous damping element (CVISC entry).

    +-------+------+-----+------+------+-----+-----+
    |   1   |  2   |  3  |  4   |   5  |  6  |  7  |
    +=======+======+=====+======+======+=====+=====+
    | PVISC | PID1 | CE1 | CR1  | PID2 | CE2 | CR2 |
    +-------+------+-----+------+------+-----+-----+
    | PVISC |  3   | 6.2 | 3.94 |      |     |     |
    +-------+------+-----+------+------+-----+-----+
    """
    def add(self, pid: int, ce: float, cr: float, comment: str='') -> int:
        self.cards.append((pid, ce, cr, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        ioffset = 0
        pid = integer(card, 1 + 4 * ioffset, 'pid')
        ce = double(card, 2 + 4 * ioffset, 'ce')
        cr = double_or_blank(card, 3 + 4 * ioffset, 'cr', default=0.)
        assert len(card) <= 8, f'len(PVISC card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, ce, cr, comment))
        self.n += 1

        if card.field(5):
            ioffset = 1
            pid = integer(card, 1 + 4 * ioffset, 'pid')
            ce = double(card, 2 + 4 * ioffset, 'ce')
            cr = double_or_blank(card, 3 + 4 * ioffset, 'cr', default=0.)
            assert len(card) <= 8, f'len(PVISC card) = {len(card):d}\ncard={card}'
            self.cards.append((pid, ce, cr, ''))
            self.n += 1
        return self.n

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        cr = np.zeros(ncards, dtype='float64')
        ce = np.zeros(ncards, dtype='float64')
        for icard, card in enumerate(self.cards):
            (pid, cei, cri, comment) = card
            property_id[icard] = pid
            cr[icard] = cri
            ce[icard] = cei
        self._save(property_id, cr, ce)
        self.cards = []


    def _save(self, property_id, cr, ce):
        self.property_id = property_id
        self.cr = cr
        self.ce = ce
        self.n = len(property_id)

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_id = array_str(self.property_id, size=size)
        for pid, ce, cr in zip_longest(property_id, self.ce, self.cr):
            cr = set_blank_if_default(cr, 0.)
            list_fields = ['PVISC', pid, ce, cr]
            bdf_file.write(print_card(list_fields))
        return


class CGAP(Element):
    """
    +------+-----+-----+-----+-----+-----+-----+------+-----+
    |  1   |  2  |  3  |  4  |  5  |  6  |  7  |   8  |  9  |
    +======+=====+=====+=====+=====+=====+=====+======+=====+
    | CGAP | EID | PID | GA  | GB  | X1  | X2  |  X3  | CID |
    +------+-----+-----+-----+-----+-----+-----+------+-----+
    | CGAP | 17  |  2  | 110 | 112 | 5.2 | 0.3 | -6.1 |     |
    +------+-----+-----+-----+-----+-----+-----+------+-----+

    or

    +------+-----+-----+-----+-----+-----+-----+------+-----+
    |  1   |  2  |  3  |  4  |  5  |  6  |  7  |   8  |  9  |
    +======+=====+=====+=====+=====+=====+=====+======+=====+
    | CGAP | EID | PID | GA  | GB  | GO  |     |      | CID |
    +------+-----+-----+-----+-----+-----+-----+------+-----+
    | CGAP | 17  |  2  | 110 | 112 | 13  |     |      |     |
    +------+-----+-----+-----+-----+-----+-----+------+-----+

    """
    def add(self, eid: int, pid: int, nids: list[int],
                 x: Optional[list[int]], g0: Optional[int],
                 cid: Optional[int]=None, comment: str='') -> int:
        self.cards.append((eid, pid, nids, x, g0, cid, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        ga = integer_or_blank(card, 3, 'ga', default=0)
        gb = integer_or_blank(card, 4, 'gb', default=0)
        x1_g0 = integer_double_or_blank(card, 5, 'x1_g0')
        cid = integer_or_blank(card, 8, 'cid', default=-1)

        if isinstance(x1_g0, integer_types):
            g0 = x1_g0
            x = None
        elif isinstance(x1_g0, float):
            g0 = None
            x1 = x1_g0
            x2 = double_or_blank(card, 6, 'x2', default=0.0)
            x3 = double_or_blank(card, 7, 'x3', default=0.0)
            x = [x1, x2, x3]
        else:
            #raise RuntimeError('invalid CGAP...x1/g0 = %r' %(x1_g0))
            g0 = -1
            x = [None, None, None]
        assert len(card) <= 9, f'len(CGAP card) = {len(card):d}\ncard={card}'

        self.cards.append((eid, pid, [ga, gb], x, g0, cid, comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        self.element_id = np.zeros(ncards, dtype='int32')
        self.property_id = np.zeros(ncards, dtype='int32')
        self.nodes = np.zeros((ncards, 2), dtype='int32')
        self.coord_id = np.zeros(ncards, dtype='int32')
        self.g0 = np.full(ncards, -1, dtype='int32')
        self.x = np.full((ncards, 3), np.nan, dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, x, g0, cid, comment) = card
            if g0 is None:
                self.x[icard, :] = x
            else:
                assert isinstance(g0, integer_types)
                self.g0[icard] = g0

            if cid is None:
                cid = -1

            self.element_id[icard] = eid
            self.property_id[icard] = pid
            self.nodes[icard, :] = nids
            self.coord_id[icard] = cid
        self.cards = []

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no PGAP properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        coord_ids = array_default_int(self.coord_id, default=-1, size=size)
        for eid, pid, nodes, g0, x, cid in zip_longest(element_id, property_id, nodes_, self.g0, self.x, coord_ids):
            ga, gb = nodes
            if g0 == -1:
                x1, x2, x3 = x # self.get_x_g0_defaults()
            else:
                x1 = g0
                x2 = ''
                x3 = ''

            list_fields = ['CGAP', eid, pid, ga, gb, x1, x2, x3, cid]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pgap]
                if prop.n > 0]

    #def mass(self) -> np.ndarray:
        #pid = self.property_id
        #mass_per_length = np.zeros(len(pid), dtype='float64')
        #for prop in self.allowed_properties:
            #i = np.searchsorted(prop.property_id, pid)

            ## filter properties that weren't found
            #ibad = (i == len(prop.property_id))
            #i[ibad] = 0

            ## make sure all properties we're setting were found
            #is_pid = (prop.property_id[i] == pid)
            #if len(is_pid) == 0:
                #continue

            ## we're at least using some properties
            #iprop = i[is_pid]
            #mass_per_lengthi = prop.mass_per_length()
            #mass_per_length[is_pid] = mass_per_lengthi[iprop]
        #length = self.length()
        #mass = mass_per_length * length
        #return mass

    #def length(self) -> np.ndarray:
        #length = line_length(self.model, self.nodes)
        #return length

    #def centroid(self) -> np.ndarray:
        #centroid = line_centroid(self.model, self.nodes)
        #return centroid

class PGAP(Property):
    """
    +------+------+-------+-------+------+------+------+------+------+
    |   1  |   2  |   3   |   4   |   5  |   6  |   7  |   8  |   9  |
    +======+======+=======+=======+======+======+======+======+======+
    | PGAP |  PID |   U0  |   F0  |  KA  |  KB  |  KT  |  MU1 |  MU2 |
    +------+------+-------+-------+------+------+------+------+------+
    |      | TMAX |  MAR  | TRMIN |      |      |      |      |      |
    +------+------+-------+-------+------+------+------+------+------+
    | PGAP |   2  | 0.025 |  2.5  | 1.E6 |      | 1.E6 | 0.25 | 0.25 |
    +------+------+-------+-------+------+------+------+------+------+
    """
    def add(self, pid: int, u0: float=0., f0: float=0.,
            ka: float=1.e8, kb: Optional[float]=None, mu1: float=0.,
            kt: Optional[float]=None, mu2: Optional[float]=None,
            tmax: float=0., mar: float=100., trmin: float=0.001,
            comment: str='') -> int:
        self.cards.append((pid, u0, f0, ka, kb, mu1, kt, mu2, tmax, mar, trmin, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        #assert isinstance(card, BDFCard), card
        pid = integer(card, 1, 'pid')
        u0 = double_or_blank(card, 2, 'u0', default=0.)
        f0 = double_or_blank(card, 3, 'f0', default=0.)
        ka = double_or_blank(card, 4, 'ka', default=1.e8)
        kb = double_or_blank(card, 5, 'kb', default=1e-14 * ka)
        mu1 = double_or_blank(card, 7, 'mu1', default=0.)
        kt = double_or_blank(card, 6, 'kt', default=mu1 * ka)
        mu2 = double_or_blank(card, 8, 'mu2', default=mu1)
        tmax = double_or_blank(card, 9, 'tmax', default=0.)
        mar = double_or_blank(card, 10, 'mar', default=100.)
        trmin = double_or_blank(card, 11, 'trmin', default=0.001)
        assert len(card) <= 12, f'len(PGAP card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, u0, f0, ka, kb, mu1, kt, mu2, tmax, mar, trmin, comment))
        #i = len(self.cards) - 1
        self.n += 1
        return self.n

    def parse_cards(self):
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return

        self.property_id = np.zeros(ncards, dtype='int32')
        self.cr = np.zeros(ncards, dtype='float64')
        self.u0 = np.zeros(ncards, dtype='float64')
        self.f0 = np.zeros(ncards, dtype='float64')
        self.ka = np.zeros(ncards, dtype='float64')
        self.kb = np.zeros(ncards, dtype='float64')
        self.kt = np.zeros(ncards, dtype='float64')
        self.mu1 = np.zeros(ncards, dtype='float64')
        self.mu2 = np.zeros(ncards, dtype='float64')
        self.tmax = np.zeros(ncards, dtype='float64')
        self.mar = np.zeros(ncards, dtype='float64')
        self.trmin = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, u0, f0, ka, kb, mu1, kt, mu2, tmax, mar, trmin, comment) = card
            self.property_id[icard] = pid
            self.u0[icard] = u0
            self.f0[icard] = f0
            self.ka[icard] = ka
            self.kb[icard] = kb
            self.kt[icard] = kt
            self.mu1[icard] = mu1
            self.mu2[icard] = mu2

            self.tmax[icard] = tmax
            self.mar[icard] = mar
            self.trmin[icard] = trmin
        self.cards = []

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_ids = array_str(self.property_id, size=size)
        for pid, u0, f0, ka, kb, kt, mu1, mu2, \
            tmax, mar, trmin in zip_longest(property_ids, self.u0, self.f0, self.ka, self.kb, self.kt,
                                            self.mu1, self.mu2, self.tmax, self.mar, self.trmin):
            u0 = set_blank_if_default(u0, 0.)
            f0 = set_blank_if_default(f0, 0.)
            # ka doesn't have a default in MSC 2005r2
            #ka = set_blank_if_default(ka, 1.e8)
            kb = set_blank_if_default(kb, 1e-14 * ka)
            kt = set_blank_if_default(kt, mu1 * ka)

            mu1 = set_blank_if_default(mu1, 0.)
            mu2 = set_blank_if_default(mu2, mu1)
            tmax = set_blank_if_default(tmax, 0.)
            mar = set_blank_if_default(mar, 100.)
            trmin = set_blank_if_default(trmin, 0.001)

            list_fields = ['PGAP', pid, u0, f0, ka, kb, kt, mu1, mu2,
                           tmax, mar, trmin]
            bdf_file.write(print_card(list_fields))
        return

