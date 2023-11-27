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
from pyNastran.bdf.bdf_interface.assign_type_force import force_double_or_blank
from pyNastran.bdf.cards.elements.bars import set_blank_if_default

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, get_print_card_8_16,
    parse_element_check, parse_property_check)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class CELAS1(Element):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CELAS1 | EID | PID | G1 | C1 | G2 | C2 |
    +--------+-----+-----+----+----+----+----+
    """
    @Element.clear_check
    def clear(self):
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 2), dtype='int32')
        self.components = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int],
            c1: int=0, c2: int=0, comment: str='') -> int:
        """
        Creates a CELAS1 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : list[int, int]
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nids, c1, c2, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        nids = [integer(card, 3, 'g1'),
                integer_or_blank(card, 5, 'g2', default=0)]

        #: component number
        c1 = integer_or_blank(card, 4, 'c1', default=0)
        c2 = integer_or_blank(card, 6, 'c2', default=0)
        assert len(card) <= 7, f'len(CELAS1 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, c1, c2, comment))
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
            (eid, pid, nids, c1, c2, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            components[icard, :] = [c1, c2]
        self._save(element_id, property_id, nodes, components)
        self.cards = []

    def _save(self, element_id, property_id, nodes, components):
        if len(self.element_id):
            raise NotImplementedError()
        nelements = len(element_id)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.components = components
        self.n = nelements

    def __apply_slice__(self, elem: CELAS1, i: np.ndarray) -> None:  # ignore[override]
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.components = self.components[i, :]
        elem.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        nodes = self.nodes.ravel()
        nodes = nodes[nodes > 0]
        used_dict['node_id'].append(nodes)
        used_dict['property_id'].append(self.property_id)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no spring properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes), filter_node0=True,
                   property_id=(pids, self.property_id))

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        components_ = array_default_int(self.components, default=0, size=size)
        for eid, pid, nodes, components in zip(element_id, property_id, nodes_, components_):
            list_fields = ['CELAS1', eid, pid,
                           nodes[0], components[0],
                           nodes[1], components[1]]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pelas]
                if prop.n > 0]


    #def length(self) -> np.ndarray:
        #length = line_length(self.model, self.nodes)
        #return length

    #def centroid(self) -> np.ndarray:
        #centroid = line_centroid(self.model, self.nodes)
        #return centroid


class CELAS2(Element):
    """
    +--------+-----+-----+----+----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |  8 |  9 |
    +========+=====+=====+====+====+====+====+====+====+
    | CELAS2 | EID |  K  | G1 | C1 | G2 | C2 | GE | S  |
    +--------+-----+-----+----+----+----+----+----+----+
    """
    @Element.clear_check
    def clear(self):
        self.element_id = np.array([], dtype='int32')
        self.k = np.array([], dtype='float64')
        self.nodes = np.zeros((0, 2), dtype='int32')
        self.components = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, k: float, nids: list[int],
            c1: int=0, c2: int=0, ge: float=0., s: float=0., comment: str='') -> int:
        """
        Creates a CELAS2 card

        Parameters
        ----------
        eid : int
            element id
        k : float
            spring stiffness
        nids : list[int, int]
            SPOINT ids
            node ids
        c1 / c2 : int; default=0
            DOF for nid1 / nid2
        ge : int; default=0.0
            damping coefficient
        s : float; default=0.0
            stress coefficient
        comment : str; default=''
            a comment for the card

        """
        nids = [0 if nid is None else nid
                for nid in nids]
        self.cards.append((eid, k, nids, c1, c2, ge, s, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        fdouble_or_blank = force_double_or_blank if self.model.is_lax_parser else double_or_blank

        eid = integer(card, 1, 'eid')
        k = double(card, 2, 'k')
        nids = [integer_or_blank(card, 3, 'g1', default=0),
                integer_or_blank(card, 5, 'g2', default=0)]
        c1 = integer_or_blank(card, 4, 'c1', default=0)
        c2 = integer_or_blank(card, 6, 'c2', default=0)
        ge = fdouble_or_blank(card, 7, 'ge', default=0.)
        s = fdouble_or_blank(card, 8, 's', default=0.)
        assert len(card) <= 9, f'len(CELAS2 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, k, nids, c1, c2, ge, s, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        k = np.zeros(ncards, dtype='float64')
        nodes = np.zeros((ncards, 2), dtype=idtype)
        components = np.zeros((ncards, 2), dtype='int32')
        ge = np.zeros(ncards, dtype='int32')
        s = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, ki, nids, c1, c2, gei, si, comment) = card
            element_id[icard] = eid
            k[icard] = ki
            nodes[icard, :] = nids
            components[icard, :] = [c1, c2]
            ge[icard] = gei
            s[icard] = si
        self._save(element_id, nodes, components, k, ge, s)
        self.cards = []

    def _save(self, element_id, nodes, components, k, ge, s):
        if len(self.element_id):
            raise NotImplementedError()
        nelements = len(element_id)
        self.element_id = element_id
        self.nodes = nodes
        self.components = components
        self.k = k
        self.ge = ge
        self.s = s
        self.n = nelements

    def __apply_slice__(self, elem: CELAS2, i: np.ndarray) -> None:  # ignore[override]
        elem.element_id = self.element_id[i]
        elem.nodes = self.nodes[i, :]
        elem.components = self.components[i, :]
        elem.k = self.k[i]
        elem.ge = self.ge[i]
        elem.s = self.s[i]
        elem.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        nodes = self.nodes.ravel()
        nodes = nodes[nodes > 0]
        used_dict['node_id'].append(nodes)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.nodes), filter_node0=True)

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        components_ = array_default_int(self.components, default=0, size=size)
        for eid, k, nodes, components, ge, s in zip(element_id, self.k,
                                                    nodes_, components_, self.ge, self.s):
            ge = set_blank_if_default(ge, 0.)
            s = set_blank_if_default(s, 0.)
            list_fields = ['CELAS2', eid, k,
                           nodes[0], components[0],
                           nodes[1], components[1], ge, s]
            bdf_file.write(print_card(list_fields))
        return


class CELAS3(Element):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CELAS3 | EID | PID | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    @Element.clear_check
    def clear(self):
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.spoints = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CELAS3 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PELAS)
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nids[0], nids[1], comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)

        s1 = integer_or_blank(card, 3, 's1', default=0)
        s2 = integer_or_blank(card, 4, 's2', default=0)
        assert len(card) <= 5, f'len(CELAS3 card) = {len(card):d}\ncard={card}'
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

    def _save(self, element_id, property_id, spoints):
        if len(self.element_id):
            raise NotImplementedError()
        nelements = len(element_id)
        self.element_id = element_id
        self.property_id = property_id
        self.spoints = spoints
        self.n = nelements

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        spoints = self.spoints.ravel()
        spoints = spoints[spoints > 0]
        used_dict['spoint_id'].append(spoints)
        used_dict['property_id'].append(self.property_id)

    def geom_check(self, missing: dict[str, np.ndarray]):
        spoint = self.model.spoint
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no spring properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   spoint=(spoint, self.spoints), filter_node0=True,
                   property_id=(pids, self.property_id))

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        spoints_ = array_default_int(self.spoints, default=0, size=size)
        for eid, pid, spoints in zip(element_id, property_id, spoints_):
            msg = 'CELAS3  %8s%8s%8s%8s\n' % (eid, pid, spoints[0], spoints[1])
            #list_fields = ['CELAS3', eid, pid, spoints[0], spoints[1]]
            bdf_file.write(msg)
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.pelas]
                if prop.n > 0]


class CELAS4(Element):
    """
    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CELAS4 | EID |  K  | S1 | S2 |
    +--------+-----+-----+----+----+
    """
    @Element.clear_check
    def clear(self):
        self.element_id = np.array([], dtype='int32')
        self.k = np.array([], dtype='float64')
        self.spoints = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, k: float, nids: list[int], comment: str='') -> int:
        """
        Creates a CELAS4 card

        Parameters
        ----------
        eid : int
            element id
        k : float
            spring stiffness
        nids : list[int, int]
            SPOINT ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, k, nids[0], nids[1], comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        k = double(card, 2, 'k')
        s1 = integer_or_blank(card, 3, 's1', default=0)
        s2 = integer_or_blank(card, 4, 's2', default=0)
        assert len(card) <= 5, f'len(CELAS4 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, k, s1, s2, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        k = np.zeros(ncards, dtype='float64')
        spoints = np.zeros((ncards, 2), dtype=idtype)

        for icard, card in enumerate(self.cards):
            (eid, ki, s1, s2, comment) = card
            element_id[icard] = eid
            k[icard] = ki
            spoints[icard, :] = [s1, s2]
        self._save(element_id, k, spoints)
        self.cards = []

    def _save(self, element_id, k, spoints):
        if len(self.element_id):
            raise NotImplementedError()
        nelements = len(element_id)
        self.element_id = element_id
        self.spoints = spoints
        self.k = k
        self.n = nelements

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        spoints = self.spoints.ravel()
        spoints = spoints[spoints > 0]
        used_dict['spoint_id'].append(spoints)


    def geom_check(self, missing: dict[str, np.ndarray]):
        spoint = self.model.spoint
        geom_check(self,
                   missing,
                   spoint=(spoint, self.spoints), filter_node0=True, )

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        spoints_ = array_default_int(self.spoints, default=0, size=size)
        for eid, k, spoints in zip(element_id, self.k, spoints_):
            list_fields = ['CELAS4', eid, k, spoints[0], spoints[1]]
            bdf_file.write(print_card(list_fields))
        return


class PELAS(Property):
    """
    Specifies the stiffness, damping coefficient, and stress coefficient of a
    scalar elastic (spring) element (CELAS1 or CELAS3 entry).
    """
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.k = np.array([], dtype='float64')
        self.ge = np.array([], dtype='float64')
        self.s = np.array([], dtype='float64')

    def add(self, pid: int, k: float, ge: float=0., s: float=0.,
            comment: str='') -> int:
        """
        Creates a PELAS card

        Parameters
        ----------
        pid : int
            property id
        k : float
            spring stiffness
        ge : int; default=0.0
            damping coefficient
        s : float; default=0.0
            stress coefficient
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, k, ge, s, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        fdouble_or_blank = force_double_or_blank if self.model.is_lax_parser else double_or_blank
        pid = integer(card, 1, 'pid')
        k = double(card, 2, 'k')
        ge = fdouble_or_blank(card, 3, 'ge', default=0.)
        s = fdouble_or_blank(card, 4, 's', default=0.)
        self.cards.append((pid, k, ge, s, comment))
        self.n += 1
        if card.field(5):
            comment = ''
            pid = integer(card, 5, 'pid')
            k = double(card, 6, 'k')
            ge = fdouble_or_blank(card, 7, 'ge', default=0.)
            s = fdouble_or_blank(card, 8, 's', default=0.)
            self.cards.append((pid, k, ge, s, comment))
            self.n += 1
        assert len(card) <= 7, f'len(PELAS card) = {len(card):d}\ncard={card}'
        return self.n - 1

        #i = len(self.cards) - 1
        #return i

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        k = np.zeros(ncards, dtype='float64')
        ge = np.zeros(ncards, dtype='float64')
        s = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, ki, gei, si, comment) = card
            property_id[icard] = pid
            k[icard] = ki
            ge[icard] = gei
            s[icard] = si
        self._save(property_id, k, ge, s)
        self.cards = []

    def _save(self, property_id, k, ge, s):
        if len(self.property_id):
            raise NotImplementedError()
        nproperties = len(property_id)
        self.property_id = property_id
        self.k = k
        self.ge = ge
        self.s = s
        self.n = nproperties

    def __apply_slice__(self, prop: PELAS, i: np.ndarray) -> None:  # ignore[override]
        prop.property_id = self.property_id[i]
        prop.k = self.k[i]
        prop.ge = self.ge[i]
        prop.s = self.s[i]
        prop.n = len(i)

    def validate(self) -> None:
        return

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['pelast_id'].append(self.property_id)

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_id = array_str(self.property_id, size=size)
        for pid, k, ge, s in zip(property_id, self.k, self.ge, self.s):
            ge = set_blank_if_default(ge, 0.)
            s = set_blank_if_default(s, 0.)
            list_fields = ['PELAS', pid, k, ge, s]
            bdf_file.write(print_card(list_fields))
        return


class PELAST(Property):
    """
    Specifies the stiffness, damping coefficient, and stress coefficient of a
    scalar elastic (spring) element (CELAS1 or CELAS3 entry).
    """
    @Property.clear_check
    def clear(self) -> None:
        self.property_id: np.array = np.array([], dtype='int32')
        self.table_k: np.array = np.array([], dtype='int32')
        self.table_ge: np.array = np.array([], dtype='int32')
        self.table_k_nonlinear: np.array = np.array([], dtype='int32')

    def add(self, pid: int, tkid: int=0, tgeid: int=0, tknid: int=0,
            comment: str='') -> int:
        """
        Creates a PELAST card

        Parameters
        ----------
        pid : int
            property id
        tkid : float
            TABLEDx that defines k vs. frequency
        tgeid : int; default=0
            TABLEDx that defines ge vs. frequency
        s : float; default=0.
            TABLEDx that defines force vs. displacement
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, tkid, tgeid, tknid, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        tkid = integer_or_blank(card, 2, 'tkid', default=0)
        tgeid = integer_or_blank(card, 3, 'tgeid', default=0)
        tknid = integer_or_blank(card, 4, 'tknid', default=0)
        assert len(card) <= 5, f'len(PELAST card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, tkid, tgeid, tknid, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')

        #: Identification number of a TABLEDi entry that defines the
        #: force per unit displacement vs. frequency relationship.
        #: (Integer > 0; Default = 0)
        table_k = np.zeros(ncards, dtype='int32')

        #: Identification number of a TABLEDi entry that defines the
        #: nondimensional structural damping coefficient vs. frequency
        #: relationship. (Integer > 0; Default = 0)
        table_ge = np.zeros(ncards, dtype='int32')

        #: Identification number of a TABELDi entry that defines the nonlinear
        #: force vs. displacement relationship. (Integer > 0; Default = 0)
        table_k_nonlinear = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (pid, tkid, tgeid, tknid, comment) = card
            property_id[icard] = pid
            table_k[icard] = tkid
            table_ge[icard] = tgeid
            table_k_nonlinear[icard] = tknid
        self._save(property_id, table_k, table_ge, table_k_nonlinear)
        self.cards = []

    def _save(self, property_id, table_k, table_ge, table_k_nonlinear):
        if len(self.property_id):
            raise NotImplementedError()
        self.property_id = property_id
        self.table_k = table_k
        self.table_ge = table_ge
        self.table_k_nonlinear = table_k_nonlinear

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_id = array_str(self.property_id, size=size)
        table_ks = array_default_int(self.table_k, default=0, size=size)
        table_ges = array_default_int(self.table_ge, default=0, size=size)
        table_kns = array_default_int(self.table_k_nonlinear, default=0, size=size)

        for pid, k, ge, kn in zip(property_id, table_ks, table_ges, table_kns):
            list_fields = ['PELAST', pid, k, ge, kn]
            bdf_file.write(print_card(list_fields))
        return
