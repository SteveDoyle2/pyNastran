from __future__ import annotations
#from itertools import count
from typing import TYPE_CHECKING
import numpy as np

#from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8, print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, integer_or_blank, double_or_blank)
#from pyNastran.bdf.cards.elements.bars import set_blank_if_default

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, parse_check)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_default_int, get_print_card_size)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg
from .utils import get_mass_from_property
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class CMASS1(Element):
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 2), dtype='int32')
        self.components = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int],
            c1: int=0, c2: int=0, comment: str='') -> int:
        """
        Creates a CMASS1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PMASS)
        nids : list[int, int]
            node ids
        c1 / c2 : int; default=None
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nids, [c1, c2], comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CMASS1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        n1 = integer(card, 3, 'g1')
        c1 = integer_or_blank(card, 4, 'c1', default=0)
        n2 = integer_or_blank(card, 5, 'g2', default=0)
        c2 = integer_or_blank(card, 6, 'c2', default=0)
        assert len(card) <= 7, f'len(CMASS1 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, [n1, n2], [c1, c2], comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        components = np.zeros((ncards, 2), dtype='int32')

        for icard, card in enumerate(self.cards):
            (eidi, pidi, nidsi, componentsi, comment) = card
            element_id[icard] = eidi
            property_id[icard] = pidi
            nodes[icard, :] = nidsi
            components[icard, :] = componentsi
        self._save(element_id, property_id, nodes, components)
        self.cards = []

    def _save(self,
              element_id: np.ndarray,
              property_id: np.ndarray,
              nodes: np.ndarray, components: np.ndarray) -> None:
        assert len(self.element_id) == 0
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.components = components

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        nodes = self.nodes.ravel()
        nodes = nodes[nodes > 0]
        used_dict['property_id'].append(self.property_id)
        used_dict['node_id'].append(nodes)

    def __apply_slice__(self, elem: CMASS1, i: np.ndarray) -> None:
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.components = self.components[i, :]
        elem.n = len(i)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        components_ = array_default_int(self.components, default=0, size=size)
        for eid, pid, (n1, n2), (c1, c2) in zip(element_id, property_id, nodes_, components_):
            list_fields = ['CMASS1', eid, pid, n1, c1, n2, c2]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pmass]
                if prop.n > 0]

    def mass(self) -> np.ndarray:
        mass = get_mass_from_property(self.property_id, self.allowed_properties)
        return mass
    #def length(self) -> np.ndarray:
        #length = line_length(self.model, self.nodes)
        #return length

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no pmass properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    def centroid(self) -> np.ndarray:
        #ispoint = np.where(self.components.ravel() == 0)[0]
        #igrid = ~ispoint
        #nodes = self.nodes.ravel()[igrid]
        #nodes = nodes.reshape(len(nodes), 1)
        centroid = point_centroid(self.model, self.nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()


def point_centroid(model, nodes: np.ndarray) -> np.ndarray:
    #unids = np.unique(nodes.ravel())
    grid = model.grid # .slice_card_by_node_id(unids)

    xyz = grid.xyz_cid0()
    nid = grid.node_id
    nelement, nnodes = nodes.shape
    node_count = np.zeros(nelement, dtype='int32')
    centroid = np.zeros((nelement, 3), dtype='float64')
    inode = np.searchsorted(nid, nodes)
    exists = (nid[inode] == nodes)
    for i in range(nnodes):
        exist = exists[:, i]
        inodei = inode[exist, i]
        centroid[exist, :] = xyz[inodei, :]
        node_count[exist] += 1

    inode_count = (node_count > 0)
    centroid[inode_count, :] /= node_count[inode_count, np.newaxis]
    assert centroid.shape[0] == nodes.shape[0]
    return centroid


class CMASS2(Element):
    """
    Defines a scalar mass element without reference to a property entry.

    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CMASS2 | EID |  M  | G1 | C1 | G2 | C2 |
    +--------+-----+-----+----+----+----+----+

    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self._mass = np.array([], dtype='float64')
        self.nodes = np.zeros((0, 2), dtype='int32')
        self.components = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, mass: float, nids: list[int],
            c1: int, c2: int, comment: str='') -> int:
        """
        Creates a CMASS2 card

        Parameters
        ----------
        eid : int
            element id
        mass : float
            mass
        nids : list[int, int]
            node ids
        c1 / c2 : int; default=None
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, mass, nids, [c1, c2], comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        mass = double_or_blank(card, 2, 'mass', default=0.)
        n1 = integer_or_blank(card, 3, 'g1', default=0)
        c1 = integer_or_blank(card, 4, 'c1', default=0)
        n2 = integer_or_blank(card, 5, 'g2', default=0)
        c2 = integer_or_blank(card, 6, 'c2', default=0)
        assert len(card) <= 7, f'len(CMASS2 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, mass, [n1, n2], [c1, c2], comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        mass = np.zeros(ncards, dtype='float64')
        nodes = np.zeros((ncards, 2), dtype=idtype)
        components = np.zeros((ncards, 2), dtype='int32')

        for icard, card in enumerate(self.cards):
            (eidi, massi, nidsi, componentsi, comment) = card
            element_id[icard] = eidi
            mass[icard] = massi
            nodes[icard, :] = nidsi
            components[icard, :] = componentsi
        self._save(element_id, mass, nodes, components)
        self.cards = []

    def _save(self, element_id: np.ndarray,
              mass: np.ndarray,
              nodes: np.ndarray, components: np.ndarray) -> None:
        assert len(self.element_id) == 0
        self.element_id = element_id
        self._mass = mass
        self.nodes = nodes
        self.components = components

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        nodes = self.nodes.ravel()
        nodes = nodes[nodes > 0]
        used_dict['node_id'].append(nodes)

    def __apply_slice__(self, elem: CMASS2, i: np.ndarray) -> None:
        elem.element_id = self.element_id[i]
        elem._mass = self._mass[i]
        elem.nodes = self.nodes[i, :]
        elem.components = self.components[i, :]
        elem.n = len(i)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        components_ = array_default_int(self.components, default=0, size=size)
        for eid, mass, (n1, n2), (c1, c2) in zip(element_id, self._mass,
                                                nodes_, components_):
            list_fields = ['CMASS2', eid, mass, n1, c1, n2, c2]
            bdf_file.write(print_card(list_fields))
        return

    def mass(self) -> np.ndarray:
        return self._mass

    def centroid(self) -> np.ndarray:
        centroid = point_centroid(self.model, self.nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()


class CMASS3(Element):
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.spoints = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CMASS3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PMASS)
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, *nids, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CMASS3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        s1 = integer_or_blank(card, 3, 's1', default=0)
        s2 = integer_or_blank(card, 4, 's2', default=0)
        assert len(card) <= 5, f'len(CMASS3 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, s1, s2, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        spoints = np.zeros((ncards, 2), dtype=idtype)

        for icard, card in enumerate(self.cards):
            (eid, pid, s1, s2, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            spoints[icard, :] = [s1, s2]
        self._save(element_id, property_id, spoints)
        self.cards = []

    def _save(self, element_id: np.ndarray,
              property_id: np.ndarray,
              spoints: np.ndarray) -> None:
        assert len(self.element_id) == 0
        assert element_id.min() > 0, element_id
        assert spoints.min() >= 0, spoints
        self.element_id = element_id
        self.property_id = property_id
        self.spoints = spoints
        self.n = len(element_id)

    def __apply_slice__(self, elem: CMASS3, i: np.ndarray) -> None:
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.spoints = self.spoints[i, :]
        elem.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['property_id'].append(self.property_id)
        used_dict['spoint_id'].append(self.spoints.ravel())

    def mass(self) -> np.ndarray:
        return np.zeros(len(self.element_id), dtype='float64')

    def center_of_mass(self) -> np.ndarray:
        return np.zeros((len(self.element_id), 3), dtype='float64')

    def centroid(self) -> np.ndarray:
        return self.center_of_mass()

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.spoints.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        spoints_ = array_default_int(self.spoints, default=0, size=size)
        for eid, pid, (spoint1, spoint2) in zip(element_id, property_id, spoints_):
            msg = 'CMASS3  %8s%8s%8s%8s\n' % (eid, pid, spoint1, spoint2)
            bdf_file.write(msg)
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pmass]
                if prop.n > 0]


class CMASS4(Element):
    """
    Defines a scalar mass element that is connected only to scalar points,
    without reference to a property entry

    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CMASS4 | EID |  M  | S1 | S2 |
    +--------+-----+-----+----+----+

    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self._mass = np.array([], dtype='float64')
        self.spoints = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, mass: float, nids: list[int], comment: str='') -> int:
        """
        Creates a CMASS4 card

        Parameters
        ----------
        eid : int
            element id
        mass : float
            SPOINT mass
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, mass, *nids, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        mass = double(card, 2, 'mass')
        s1 = integer_or_blank(card, 3, 's1', default=0)
        s2 = integer_or_blank(card, 4, 's2', default=0)
        self.cards.append((eid, mass, s1, s2, comment))
        self.n += 1
        if card.field(5):
            eid = integer(card, 5, 'eid')
            mass = double(card, 6, 'mass')
            s1 = integer_or_blank(card, 7, 's1', default=0)
            s2 = integer_or_blank(card, 8, 's2', default=0)
            self.cards.append((eid, mass, s1, s2, comment))
            self.n += 1
        assert len(card) <= 9, f'len(CMASS4 card) = {len(card):d}\ncard={card}'
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        mass = np.zeros(ncards, dtype='float64')
        spoints = np.zeros((ncards, 2), dtype=idtype)

        for icard, card in enumerate(self.cards):
            (eid, massi, s1, s2, comment) = card
            element_id[icard] = eid
            mass[icard] = massi
            spoints[icard, :] = [s1, s2]
        self._save(element_id, mass, spoints)
        self.cards = []

    def _save(self, element_id, mass, spoints):
        assert len(self.element_id) == 0
        assert element_id.min() > 0, element_id
        assert spoints.min() >= 0, spoints
        self.element_id = element_id
        self._mass = mass
        self.spoints = spoints
        self.n = len(element_id)

    def __apply_slice__(self, elem: CMASS4, i: np.ndarray) -> None:
        elem.element_id = self.element_id[i]
        elem._mass = self._mass[i]
        elem.spoints = self.spoints[i, :]
        elem.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['spoint_id'].append(self.spoints.ravel())

    def mass(self) -> np.ndarray:
        return np.zeros(len(self.element_id), dtype='float64')

    def center_of_mass(self) -> np.ndarray:
        return np.zeros((len(self.element_id), 3), dtype='float64')

    def centroid(self) -> np.ndarray:
        return self.center_of_mass()

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.spoints.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        spoints_ = array_default_int(self.spoints, default=0, size=size)
        for eid, massi, spoints in zip(element_id, self._mass, spoints_):
            list_fields = ['CMASS4', eid, massi, spoints[0], spoints[1]]
            bdf_file.write(print_card(list_fields))
        return


class PMASS(Property):
    """
    Scalar Mass Property
    Specifies the mass value of a scalar mass element (CMASS1 or CMASS3 entries).

    +-------+------+------+------+------+------+----+------+----+
    |   1   |   2  |   3  |   4  |   5  |   6  |  7 |   8  |  9 |
    +=======+======+======+======+======+======+====+======+====+
    | PMASS | PID1 |  M1  | PID2 |  M2  | PID3 | M3 | PID4 | M4 |
    +-------+------+------+------+------+------+----+------+----+
    | PMASS |   7  | 4.29 |   6  | 13.2 |      |    |      |    |
    +-------+------+------+------+------+------+----+------+----+
    """
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self._mass = np.array([], dtype='float64')

    def add(self, pid: int, mass: float, comment: str='') -> int:
        """
        Creates an PMASS card, which defines a mass applied to a single DOF

        Parameters
        ----------
        pid : int
            Property id used by a CMASS1/CMASS3 card
        mass : float
            the mass to apply
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, mass, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        for icard, j in enumerate([1, 3, 5, 7]):
            if card.field(j):
                ioffset = icard * 2
                pid = integer(card, 1 + ioffset, 'pid')
                mass = double_or_blank(card, 2 + ioffset, 'mass', default=0.)
                self.cards.append((pid, mass, comment))
                comment = ''
                self.n += 1
        assert len(card) <= 9, f'len(PMASS card) = {len(card):d}\ncard={card}'
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        property_id = np.zeros(ncards, dtype=idtype)
        mass = np.zeros(ncards, dtype='float64')
        for icard, card in enumerate(self.cards):
            (pid, massi, comment) = card
            property_id[icard] = pid
            mass[icard] = massi
        self._save(property_id, mass)
        self.cards = []

    def _save(self, property_id: np.ndarray, mass: np.ndarray) -> None:
        assert len(self.property_id) == 0
        self.property_id = property_id
        self._mass = mass
        self.n = len(property_id)

    def __apply_slice__(self, prop: PMASS, i: np.ndarray) -> None:
        prop.property_id = self.property_id[i]
        prop._mass = self._mass[i]
        prop.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    @property
    def max_id(self) -> int:
        return self.property_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        property_id = array_str(self.property_id, size=size)
        for pid, mass in zip(property_id, self._mass):
            list_fields = ['PMASS', pid, mass]
            bdf_file.write(print_card(list_fields))
        return

    def mass(self) -> np.ndarray:
        return self._mass
