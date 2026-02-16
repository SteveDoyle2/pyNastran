from __future__ import annotations
from itertools import zip_longest
from typing import Optional, Any, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, integer_or_blank, double_or_blank,
    integer_double_or_blank)
from pyNastran.bdf.cards.elements.bars import set_blank_if_default

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, parse_check, save_ifile_comment)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    get_print_card_size,
    array_str, array_float, array_float_nan,
    array_default_int, array_default_float, array_default_floats)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg
from .bar import get_bar_vector, safe_normalize, line_length

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class CDAMP1(Element):
    """
    +--------+-----+-----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |
    +========+=====+=====+====+====+====+====+
    | CDAMP1 | EID | PID | G1 | C1 | G2 | C2 |
    +--------+-----+-----+----+----+----+----+
    """
    _skip_equal_check = ['allowed_properties']
    @Element.clear_check
    def clear(self):
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 2), dtype='int32')
        self.components = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int], c1: int=0, c2: int=0,
            ifile: int=0, comment: str='') -> int:
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
        self.cards.append((eid, pid, nids, [c1, c2], ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        nids = [integer_or_blank(card, 3, 'g1', default=0),
                integer_or_blank(card, 5, 'g2', default=0)]

        #: component number
        c1 = integer_or_blank(card, 4, 'c1', default=0)
        c2 = integer_or_blank(card, 6, 'c2', default=0)
        assert len(card) <= 7, f'len(CDAMP1 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, [c1, c2], ifile, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        components = np.zeros((ncards, 2), dtype='int32')
        comment = {}

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, componentsi, ifilei, commenti) = card
            ifile[icard] = ifilei
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            components[icard, :] = componentsi
        self._save(element_id, property_id, nodes, components,
                   ifile=ifile)
        self.cards = []

    def _save(self, element_id, property_id, nodes, components,
              ifile=None, comment=None):
        if ifile is None:
            ncards = len(element_id)
            ifile = np.zeros(ncards, dtype='int32')
        save_ifile_comment(self, ifile, comment)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.components = components
        self.n = len(ifile)

    def __apply_slice__(self, elem: CDAMP1, i: np.ndarray) -> None:
        elem.ifile = self.ifile[i]
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.components = self.components[i]
        elem.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        nodes = self.nodes.ravel()
        nodes = nodes[nodes > 0]
        used_dict['property_id'].append(self.property_id)
        used_dict['node_id'].append(nodes)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no damper properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))
    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(),
                   self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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
    @Element.clear_check
    def clear(self):
        self.element_id = np.array([], dtype='int32')
        self.b = np.array([], dtype='float64')
        self.nodes = np.zeros((0, 2), dtype='int32')
        self.components = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, b: float, nids: list[int],
            c1: int=0, c2: int=0,
            ifile: int=0, comment: str='') -> int:
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
        self.cards.append((eid, b, nids, (c1, c2), ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        b = double(card, 2, 'b')
        nids = [integer_or_blank(card, 3, 'g1', default=0),
                integer_or_blank(card, 5, 'g2', default=0)]
        c1 = integer_or_blank(card, 4, 'c1', default=0)
        c2 = integer_or_blank(card, 6, 'c2', default=0)
        assert len(card) <= 7, f'len(CDAMP2 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, b, nids, (c1, c2), ifile, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        b = np.zeros(ncards, dtype='float64')
        nodes = np.zeros((ncards, 2), dtype=idtype)
        components = np.zeros((ncards, 2), dtype='int32')

        for icard, card in enumerate(self.cards):
            eid, bi, nids, componentsi, ifilei, comment = card
            ifile[icard] = ifilei
            element_id[icard] = eid
            b[icard] = bi
            nodes[icard, :] = nids
            components[icard, :] = componentsi
        self._save(element_id, nodes, components, b,
                   ifile=ifile)
        self.cards = []

    def _save(self, element_id, nodes, components, b,
              ifile=None, comment=None) -> None:
        ncards = len(element_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        save_ifile_comment(self, ifile, comment)
        self.element_id = element_id
        self.nodes = nodes
        self.components = components
        self.b = b
        self.n = len(ifile)

    def __apply_slice__(self, elem: CDAMP2, i: np.ndarray) -> None:
        elem.ifile = self.ifile[i]
        elem.element_id = self.element_id[i]
        elem.nodes = self.nodes[i, :]
        elem.components = self.components[i]
        elem.b = self.b[i]
        elem.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['node_id'].append(self.nodes.ravel())

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        geom_check(self,
                   missing,
                   node=(nid, self.nodes))
    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        bs = array_float(self.b, size=size, is_double=False)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        components_ = array_default_int(self.components, default=0, size=size)
        for eid, b, nodes, components in zip(element_id, bs,
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
    _skip_equal_check = ['allowed_properties']
    @Element.clear_check
    def clear(self):
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.spoints = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int],
            ifile: int=0, comment: str='') -> int:
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
        self.cards.append((eid, pid, nids, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)

        s1 = integer_or_blank(card, 3, 's1', default=0)
        s2 = integer_or_blank(card, 4, 's2', default=0)
        assert len(card) <= 5, f'len(CDAMP3 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, (s1, s2), ifile, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        spoints = np.zeros((ncards, 2), dtype=idtype)

        for icard, card in enumerate(self.cards):
            (eid, pid, spointsi, ifilei, comment) = card
            ifile[icard] = ifilei
            element_id[icard] = eid
            property_id[icard] = pid
            spoints[icard, :] = spointsi
        self._save(element_id, property_id, spoints,
                   ifile=ifile)
        self.cards = []

    def _save(self, element_id, property_id, spoints,
              ifile=None, comment=None) -> None:
        ncards = len(element_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        save_ifile_comment(self, ifile, comment)
        self.element_id = element_id
        self.property_id = property_id
        self.spoints = spoints
        self.n = len(ifile)

    def __apply_slice__(self, elem: CDAMP3, i: np.ndarray) -> None:
        elem.ifile = self.ifile[i]
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.spoints = self.spoints[i, :]
        elem.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['property_id'].append(self.property_id)
        used_dict['spoint_id'].append(self.spoints.ravel())

    def geom_check(self, missing: dict[str, np.ndarray]):
        spoint = self.model.spoint
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no damper properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   spoint=(spoint, self.spoints),
                   property_id=(pids, self.property_id))

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.spoints.max())

    @parse_check
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
    @Element.clear_check
    def clear(self):
        self.element_id = np.array([], dtype='int32')
        self.b = np.array([], dtype='float64')
        self.spoints = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, b: float, nids: list[int],
            ifile: int=0, comment: str='') -> int:
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
        self.cards.append((eid, b, nids, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        b = double(card, 2, 'b')
        s1 = integer_or_blank(card, 3, 's1', default=0)
        s2 = integer_or_blank(card, 4, 's2', default=0)
        self.cards.append((eid, b, [s1, s2], ifile, comment))
        self.n += 1
        if card.field(5):
            eid = integer(card, 5, 'eid')
            b = double(card, 6, 'b')
            s1 = integer_or_blank(card, 7, 's1', default=0)
            s2 = integer_or_blank(card, 8, 's2', default=0)
            self.cards.append((eid, b, [s1, s2], ifile, comment))
            self.n += 1
        assert len(card) <= 9, f'len(CDAMP4 card) = {len(card):d}\ncard={card}'
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        b = np.zeros(ncards, dtype='float64')
        spoints = np.zeros((ncards, 2), dtype=idtype)

        for icard, card in enumerate(self.cards):
            eid, bi, spointsi, ifilei, comment = card
            ifile[icard] = ifilei
            element_id[icard] = eid
            b[icard] = bi
            spoints[icard, :] = spointsi
        self._save(element_id, b, spoints, ifile=ifile)
        self.cards = []

    def _save(self, element_id, b, spoints,
              ifile=None, comment=None) -> None:
        ncards = len(element_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')

        if len(self.element_id):
            asdf
        save_ifile_comment(self, ifile, comment)
        self.element_id = element_id
        self.b = b
        self.spoints = spoints
        self.n = len(ifile)

    def __apply_slice__(self, elem: CDAMP4, i: np.ndarray) -> None:
        elem.ifile = self.ifile[i]
        elem.element_id = self.element_id[i]
        elem.b = self.b[i]
        elem.spoints = self.spoints[i, :]
        elem.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['spoint_id'].append(self.spoints.ravel())

    def geom_check(self, missing: dict[str, np.ndarray]):
        spoint = self.model.spoint
        geom_check(self,
                   missing,
                   spoint=(spoint, self.spoints))

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.spoints.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        bs = array_float(self.b, size=size, is_double=False)
        spoints_ = array_default_int(self.spoints, default=0, size=size)
        for eid, b, spoints in zip(element_id, bs, spoints_):
            list_fields = ['CDAMP4', eid, b, spoints[0], spoints[1]]
            bdf_file.write(print_card(list_fields))
        return


class CDAMP5(Element):
    """
    Defines a damping element that refers to a material property entry and
    connection to grid or scalar points.

    +--------+-----+-----+----+----+
    |    1   |  2  |  3  |  4 |  5 |
    +========+=====+=====+====+====+
    | CDAMP5 | EID | PID | N1 | N2 |
    +--------+-----+-----+----+----+
    """
    @Element.clear_check
    def clear(self):
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int],
            ifile: int=0, comment: str='') -> int:
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
        self.cards.append((eid, pid, nids, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 'n1', default=0),
                integer_or_blank(card, 4, 'n2', default=0)]
        assert len(card) <= 5, f'len(CDAMP5 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, ifile, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)

        for icard, card in enumerate(self.cards):
            eid, pid, nids, ifilei, comment = card
            ifile[icard] = ifilei
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes,
                   ifile=ifile)
        self.cards = []

    def _save(self, element_id, property_id, nodes,
              ifile=None, comment=None) -> None:
        ncards = len(element_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')

        if len(self.element_id):
            asdf
        save_ifile_comment(self, ifile, comment)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.n = len(ifile)

    def __apply_slice__(self, elem: CDAMP4, i: np.ndarray) -> None:
        elem.ifile = self.ifile[i]
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.n = len(i)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(),
                   self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.b = np.array([], dtype='float64')

    def add(self, pid: int, b: float,
            ifile: int=0, comment: str='') -> int:
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
        self.cards.append((pid, b, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> list[int]:
        """adds a PDAMP"""
        ns = [self.n]
        pid = integer(card, 1, 'pid')
        b = double(card, 2, 'b')
        self.cards.append((pid, b, ifile, comment))
        self.n += 1

        if card.field(3):
            pid = integer(card, 3, 'pid2')
            b = double(card, 4, 'b2')
            self.cards.append((pid, b, ifile, comment))
            ns.append(self.n)
            self.n += 1

        if card.field(5):
            pid = integer(card, 5, 'pid3')
            b = double(card, 6, 'b3')
            self.cards.append((pid, b, ifile, comment))
            ns.append(self.n)
            self.n += 1

        if card.field(7):
            pid = integer(card, 7, 'pid4')
            b = double(card, 8, 'b4')
            self.cards.append((pid, b, ifile, comment))
            ns.append(self.n)
            self.n += 1

        assert len(card) <= 9, f'len(PDAMP card) = {len(card):d}\ncard={card}'
        return ns

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype=idtype)
        b = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, bi, ifilei, commenti) = card
            ifile[icard] = ifilei
            property_id[icard] = pid
            b[icard] = bi
        self._save(property_id, b, ifile=ifile)
        self.cards = []

    def _save(self, property_id, b, ifile=None, comment=None):
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.property_id) != 0:
            ifile = np.hstack([self.ifile, ifile])
            property_id = np.hstack([self.property_id, property_id])
            b = np.hstack([self.b, b])
        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.b = b
        self.n = len(ifile)

    def __apply_slice__(self, prop: PDAMP, i: np.ndarray) -> None:
        prop.ifile = self.ifile[i]
        prop.property_id = self.property_id[i]
        prop.b = self.b[i]
        prop.n = len(i)

    def validate(self) -> None:
        return

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['pdampt_id'].append(self.property_id)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        return self.property_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.b = np.array([], dtype='float64')

    def add(self, pid: int, mid: int, b: float,
            ifile: int=0, comment: str='') -> int:
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
        self.cards.append((pid, mid, b, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> list[int]:
        """adds a PDAMP"""
        #ns = [self.n]
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        b = double(card, 3, 'b')
        self.cards.append((pid, mid, b, ifile, comment))
        self.n += 1

        assert len(card) <= 4, f'len(PDAMP5 card) = {len(card):d}\ncard={card}'
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype=idtype)
        material_id = np.zeros(ncards, dtype=idtype)
        b = np.zeros(ncards, dtype='float64')
        comment = {}

        for icard, card in enumerate(self.cards):
            (pid, mid, bi, ifilei, commenti) = card
            ifile[icard] = ifilei
            property_id[icard] = pid
            material_id[icard] = mid
            b[icard] = bi
        self._save(property_id, material_id, b, ifile=ifile)
        self.cards = []

    def _save(self, property_id, material_id, b,
              ifile=None, comment=None):
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.material_id = material_id
        self.b = b
        self.n = len(ifile)

    def validate(self) -> None:
        return

    @property
    def allowed_materials(self) -> list[Any]:
        return [mat for mat in [self.model.mat4, self.model.mat5]
                if mat.n > 0]

    def geom_check(self, missing: dict[str, np.ndarray]) -> None:
        mids = hstack_msg([mat.material_id for prop in self.allowed_materials],
                          msg=f'no thermal materials for {self.type}')
        mids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(mids, self.material_id))

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.material_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        property_id = array_str(self.property_id, size=size)
        material_id = array_str(self.material_id, size=size)
        for pid, mid, b in zip(property_id, material_id, self.b):
            list_fields = ['PDAMP5', pid, mid, b]
            bdf_file.write(print_card(list_fields))
        return


class PDAMPT(Property):
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.table_b = np.array([], dtype='int32')

    def add(self, pid: int, tbid: int,
            ifile: int=0, comment: str='') -> int:
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
        self.cards.append((pid, tbid, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        tbid = integer_or_blank(card, 2, 'tbid', 0)
        assert len(card) <= 3, f'len(PDAMPT card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, tbid, ifile, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')

        #: Identification number of a TABLEDi entry that defines the
        #: damping force per-unit velocity versus frequency relationship
        table_b = np.zeros(ncards, dtype='int32')
        comment = {}

        for icard, (pid, tbid, ifilei, commenti) in enumerate(self.cards):
            ifile[icard] = ifilei
            property_id[icard] = pid
            table_b[icard] = tbid
        self._save(property_id, table_b, ifile=ifile)
        self.cards = []

    def _save(self, property_id, table_b,
              ifile=None, comment=None) -> None:
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.property_id) != 0:
            asdf
        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.table_b = table_b
        self.n = len(ifile)

    def __apply_slice__(self, prop: PDAMPT, i: np.ndarray) -> None:
        prop.ifile = self.ifile[i]
        prop.property_id = self.property_id[i]
        prop.table_b = self.table_b[i]
        prop.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['tabled_id'].append(self.table_b)

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.table_b.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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
    @Element.clear_check
    def clear(self) -> None:
        self.ifile = np.array([], dtype='int32')
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int],
            ifile: int=0, comment: str='') -> int:
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
        self.cards.append((eid, pid, nids, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        nids = [integer(card, 3, 'g1'),
                integer_or_blank(card, 5, 'g2', default=0)]

        assert len(card) <= 6, f'len(CVISC card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, ifile, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, ifilei, commenti) = card
            ifile[icard] = ifilei
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes, ifile=ifile)
        self.cards = []

    def _save(self, element_id, property_id, nodes,
              ifile=None, comment=None):
        ncards = len(element_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.property_id) != 0:
            asdf
        assert ncards > 0, ncards
        save_ifile_comment(self, ifile, comment)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.n = len(ifile)

    def __apply_slice__(self, elem: CVISC, i: np.ndarray) -> None:
        elem.ifile = self.ifile[i]
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        nodes = self.nodes.ravel()
        nodes = nodes[nodes > 0]
        used_dict['property_id'].append(self.property_id)
        used_dict['node_id'].append(nodes)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no visc properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(),
                   self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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

    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.cr = np.array([], dtype='float64')
        self.ce = np.array([], dtype='float64')

    def add(self, pid: int, ce: float, cr: float,
            ifile: int=0, comment: str='') -> int:
        self.cards.append((pid, ce, cr, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        ioffset = 0
        pid = integer(card, 1 + 4 * ioffset, 'pid')
        ce = double(card, 2 + 4 * ioffset, 'ce')
        cr = double_or_blank(card, 3 + 4 * ioffset, 'cr', default=0.)
        assert len(card) <= 8, f'len(PVISC card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, ce, cr, ifile, comment))
        self.n += 1

        if card.field(5):
            ioffset = 1
            pid = integer(card, 1 + 4 * ioffset, 'pid')
            ce = double(card, 2 + 4 * ioffset, 'ce')
            cr = double_or_blank(card, 3 + 4 * ioffset, 'cr', default=0.)
            assert len(card) <= 8, f'len(PVISC card) = {len(card):d}\ncard={card}'
            self.cards.append((pid, ce, cr, ifile, ''))
            self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        cr = np.zeros(ncards, dtype='float64')
        ce = np.zeros(ncards, dtype='float64')
        for icard, card in enumerate(self.cards):
            (pid, cei, cri, ifilei, commenti) = card
            ifile[icard] = ifilei
            property_id[icard] = pid
            cr[icard] = cri
            ce[icard] = cei
        self._save(property_id, ce, cr, ifile=ifile)
        self.cards = []

    def _save(self, property_id, ce, cr,
              ifile=None, comment=None):
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.property_id):
            asdf
        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.ce = ce
        self.cr = cr
        self.n = len(ifile)

    def __apply_slice__(self, prop: PVISC, i: np.ndarray) -> None:
        prop.ifile = self.ifile[i]
        prop.property_id = self.property_id[i]
        prop.cr = self.cr[i]
        prop.ce = self.ce[i]
        prop.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass
    @property
    def max_id(self) -> int:
        return self.property_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        property_id = array_str(self.property_id, size=size)
        ces = array_float(self.ce, size=size, is_double=is_double)
        crs = array_default_float(self.ce, default=0., size=size, is_double=is_double)
        for pid, ce, cr in zip_longest(property_id, ces, crs):
            #cr = set_blank_if_default(cr, 0.)
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
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 2), dtype='int32')
        self.coord_id = np.array([], dtype='int32')
        self.g0 = np.array([], dtype='int32')
        self.x = np.zeros((0, 3), dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int],
            x: Optional[list[int]], g0: Optional[int],
            cid: Optional[int]=None,
            ifile: int=0, comment: str='') -> int:
        self.cards.append((eid, pid, nids, x, g0, cid, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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

        self.cards.append((eid, pid, [ga, gb], x, g0, cid, ifile, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        coord_id = np.zeros(ncards, dtype='int32')
        g0 = np.full(ncards, -1, dtype=idtype)
        x = np.full((ncards, 3), np.nan, dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, xi, g0i, cid, ifilei, commenti) = card
            if g0i is None or g0i == -1:
                x[icard, :] = xi
            else:
                assert isinstance(g0i, integer_types), g0i
                g0[icard] = g0i

            if cid is None:
                cid = -1

            ifile[icard] = ifilei
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            coord_id[icard] = cid
        self._save(element_id, property_id, nodes, coord_id, x, g0,
                   ifile=ifile)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, nodes, coord_id, x, g0,
              ifile=None, comment=None) -> None:
        ncards = len(element_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.element_id):
            raise RuntimeError(f'stacking of {self.type} is not supported')
        save_ifile_comment(self, ifile, comment)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.coord_id = coord_id
        self.x = x
        self.g0 = g0

    def convert(self, xyz_scale: float=1.0, **kwargs) -> None:
        """
        x is i
        """
        icoord = (self.coord_id != -1)
        ncoord = icoord.sum()
        is_x = self.is_x

        #xmax = self.x.max(axis=1)
        #xnan = np.isnan(xmax)
        nx = is_x.sum()
        if nx:
            ga = self.nodes[is_x, 0]
            grid_ga = self.model.grid.slice_card_by_id(ga, assume_sorted=True, sort_ids=False)
            cp = grid_ga.cp
            coords = self.model.coord.slice_card_by_id(cp, assume_sorted=True, sort_ids=False)
            xi = self.x[is_x, :].copy()
            ixyz = (coords.coord_type == 'R')
            irtz = (coords.coord_type == 'C')
            irtp = (coords.coord_type == 'S')
            xi[ixyz, :] *= xyz_scale
            xi[irtz, 0] *= xyz_scale
            xi[irtz, 2] *= xyz_scale
            xi[irtp, 0] *= xyz_scale
            self.x[is_x] = xi

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        coords = self.coord_id[self.coord_id >= 0]
        used_dict['property_id'].append(self.property_id)
        used_dict['node_id'].append(self.nodes.ravel())
        g0 = self.g0[self.g0 > 0]
        used_dict['node_id'].append(g0)
        used_dict['coord_id'].append(coords)

    def __apply_slice__(self, elem: CGAP, i: np.ndarray) -> None:
        elem.ifile = self.ifile[i]
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.g0 = self.g0[i]
        elem.x = self.x[i, :]
        elem.coord_id = self.coord_id[i]
        elem.n = len(i)

    def length(self) -> np.ndarray:
        length = line_length(self.model, self.nodes)
        return length

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no PGAP properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    @property
    def is_x(self) -> np.ndarray:
        return self.g0 == -1

    @property
    def is_g0(self) -> np.ndarray:
        return ~self.is_x

    def get_xyz(self) -> tuple[np.ndarray, np.ndarray]:
        #neids = len(self.element_id)
        grid = self.model.grid
        xyz = grid.xyz_cid0()
        nid = grid.node_id
        inode = np.searchsorted(nid, self.nodes)
        assert np.array_equal(nid[inode], self.nodes)
        in1 = inode[:, 0]
        in2 = inode[:, 1]
        xyz1 = xyz[in1, :]
        xyz2 = xyz[in2, :]
        return xyz1, xyz2

    def get_bar_vector(self, xyz1: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        idefault = (self.coord_id == -1)
        icoord = ~idefault

        i = np.full(self.x.shape, np.nan, self.x.dtype)
        j = np.full(self.x.shape, np.nan, self.x.dtype)
        k = np.full(self.x.shape, np.nan, self.x.dtype)
        xyz1, xyz2 = self.get_xyz()
        i_vector = xyz2 - xyz1
        nelement = self.n
        maxs = np.abs(i_vector).max(axis=1)
        assert len(maxs) == nelement, (len(maxs), nelement)

        if np.any(idefault):
            sub_i = i_vector[idefault, :]
            ihat = safe_normalize(sub_i)

            index = np.where(idefault)[0]
            sub_cgap = self.slice_card_by_index(index)
            v, cd = get_bar_vector(sub_cgap, xyz1[idefault, :])

            ki = np.cross(ihat, v, axis=1)
            khat = safe_normalize(ki)
            jhat = np.cross(khat, ihat, axis=1)
            i[idefault, :] = ihat
            j[idefault, :] = jhat
            k[idefault, :] = khat

        if np.any(icoord):
            #print('icoord =', icoord)
            coord_ids = self.coord_id[icoord]
            #print('coord_ids =', coord_ids)
            assert len(coord_ids) == icoord.sum()
            coords = self.model.coord.slice_card_by_id(coord_ids)
            ihat = coords.i
            jhat = coords.j
            khat = coords.k
            i[icoord, :] = ihat
            j[icoord, :] = jhat
            k[icoord, :] = khat
        #if 0:
            #v, cd = get_bar_vector(self, xyz1)
        return i, j, k

    def get_axes(self, xyz1: np.ndarray, xyz2: np.ndarray,
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                            np.ndarray, np.ndarray, np.ndarray]:
        log = self.model.log
        coords = self.model.coord
        #xyz1, xyz2 = self.get_xyz()

        neids = xyz1.shape[0]
        #i = xyz2 - xyz1
        #ihat_norm = np.linalg.norm(i, axis=1)
        #assert len(ihat_norm) == neids
        #if min(ihat_norm) == 0.:
            #msg = 'xyz1=%s xyz2=%s\n%s' % (xyz1, xyz2, self)
            #log.error(msg)
            #raise ValueError(msg)
        #i_offset = i / ihat_norm[:, np.newaxis]

        #log.info(f'x =\n{self.x}')
        #log.info(f'g0   = {self.g0}')
        ihat, yhat, zhat = self.get_bar_vector(xyz1)

        v = yhat
        wa = np.zeros(ihat.shape, ihat.dtype)
        wb = wa
        #ihat = xform[0, :]
        #yhat = xform[1, :]
        #zhat = xform[2, :]
        #wa, wb, _ihat, jhat, khat = out

        # we finally have the nodal coordaintes!!!! :)
        return v, ihat, yhat, zhat, wa, wb

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(),
                   self.nodes.max(), self.g0.max(), self.coord_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        g0s = array_default_int(self.g0, default=-1, size=size)
        coord_ids = array_default_int(self.coord_id, default=-1, size=size)
        xs = array_float_nan(self.x, size=size, is_double=is_double)
        for eid, pid, nodes, g0, x, cid in zip_longest(element_id, property_id, nodes_, g0s, xs, coord_ids):
            ga, gb = nodes
            if g0 == '':
                x1, x2, x3 = x
            else:
                x1 = g0
                x2 = ''
                x3 = ''

            list_fields = ['CGAP', eid, pid, ga, gb, x1, x2, x3, cid]
            assert all(['nan' not in field for field in list_fields]), list_fields
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

    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.u0 = np.array([], dtype='float64')
        self.f0 = np.array([], dtype='float64')
        self.ka = np.array([], dtype='float64')
        self.kb = np.array([], dtype='float64')
        self.kt = np.array([], dtype='float64')
        self.mu1 = np.array([], dtype='float64')
        self.mu2 = np.array([], dtype='float64')
        self.tmax = np.array([], dtype='float64')
        self.mar = np.array([], dtype='float64')
        self.trmin = np.array([], dtype='float64')

    def add(self, pid: int, u0: float=0., f0: float=0.,
            ka: float=1.e8, kb: Optional[float]=None, mu1: float=0.,
            kt: Optional[float]=None, mu2: Optional[float]=None,
            tmax: float=0., mar: float=100., trmin: float=0.001,
            ifile: int=0, comment: str='') -> int:
        self.cards.append((pid, u0, f0, ka, kb, mu1, kt, mu2, tmax, mar, trmin, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
        self.cards.append((pid, u0, f0, ka, kb, mu1, kt, mu2, tmax, mar, trmin, ifile, comment))
        #i = len(self.cards) - 1
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        u0 = np.zeros(ncards, dtype='float64')
        f0 = np.zeros(ncards, dtype='float64')
        ka = np.zeros(ncards, dtype='float64')
        kb = np.zeros(ncards, dtype='float64')
        kt = np.zeros(ncards, dtype='float64')
        mu1 = np.zeros(ncards, dtype='float64')
        mu2 = np.zeros(ncards, dtype='float64')
        tmax = np.zeros(ncards, dtype='float64')
        mar = np.zeros(ncards, dtype='float64')
        trmin = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, u0i, f0i, kai, kbi, mu1i, kti, mu2i, tmaxi, mari, trmini, ifilei, commenti) = card
            ifile[icard] = ifilei
            property_id[icard] = pid
            u0[icard] = u0i
            f0[icard] = f0i
            ka[icard] = kai
            kb[icard] = kbi
            kt[icard] = kti
            mu1[icard] = mu1i
            mu2[icard] = mu2i

            tmax[icard] = tmaxi
            mar[icard] = mari
            trmin[icard] = trmini
        self._save(property_id, u0, f0, ka, kb, kt, mu1, mu2, tmax, mar, trmin,
                   ifile=ifile)
        self.cards = []

    def _save(self, property_id, u0, f0, ka, kb, kt,
              mu1, mu2, tmax, mar, trmin,
              ifile=None, comment=None) -> None:
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.property_id) != 0:
            raise RuntimeError(f'stacking of {self.type} is not supported')
        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.u0 = u0
        self.f0 = f0
        self.ka = ka
        self.kb = kb
        self.kt = kt
        self.mu1 = mu1
        self.mu2 = mu2
        self.tmax = tmax
        self.mar = mar
        self.trmin = trmin
        self.n = len(property_id)

    def __apply_slice__(self, prop: PGAP, i: np.ndarray) -> None:
        prop.ifile = self.ifile[i]
        prop.property_id = self.property_id[i]
        prop.u0 = self.u0[i]
        prop.f0 = self.f0[i]
        prop.ka = self.ka[i]
        prop.kb = self.kb[i]
        prop.kt = self.kt[i]
        prop.mu1 = self.mu1[i]
        prop.mu2 = self.mu2[i]
        prop.tmax = self.tmax[i]
        prop.mar = self.mar[i]
        prop.trmin = self.trmin[i]
        prop.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        pass

    def convert(self, xyz_scale: float=1.0,
                force_scale: float=1.0,
                stiffness_scale: float=1.0, **kwargs) -> None:
        """
        u0 : float; default=0.
            Initial gap opening
        f0 : float; default=0.
            Preload
        ka : float; default=1.e8
            Axial stiffness for the closed gap
        kb : float; default=None -> 1e-14 * ka
            Axial stiffness for the open gap
        mu1 : float; default=0.
            Coefficient of static friction for the adaptive gap element
            or coefficient of friction in the y transverse direction
            for the nonadaptive gap element
        kt : float; default=None -> mu1*ka
            Transverse stiffness when the gap is closed
        mu2 : float; default=None -> mu1
            Coefficient of kinetic friction for the adaptive gap element
            or coefficient of friction in the z transverse direction
            for the nonadaptive gap element
        tmax : float; default=0.
            Maximum allowable penetration used in the adjustment of
            penalty values. The positive value activates the penalty
            value adjustment
        mar : float; default=100.
            Maximum allowable adjustment ratio for adaptive penalty
            values KA and KT
        trmin : float; default=0.001
            Fraction of TMAX defining the lower bound for the allowable
            penetration

        """
        self.u0 *= xyz_scale
        self.f0 *= force_scale
        self.ka *= stiffness_scale
        self.kb *= stiffness_scale
        self.tmax *= xyz_scale

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        return self.property_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        property_ids = array_str(self.property_id, size=size)
        u0s = array_default_float(self.u0, default=0.0, size=size, is_double=is_double)
        f0s = array_default_float(self.f0, default=0.0, size=size, is_double=is_double)
        mu1s = array_default_float(self.mu1, default=0.0, size=size, is_double=is_double)
        mars = array_default_float(self.mar, default=100.0, size=size, is_double=is_double)
        tmaxs = array_default_float(self.tmax, default=0.0, size=size, is_double=is_double)
        trmins = array_default_float(self.trmin, default=0.001, size=size, is_double=is_double)

        kas = array_default_float(self.ka, default=1e8, size=size, is_double=is_double)
        kb_default = 1e-14 * self.ka
        kbs = array_default_floats(self.kb, kb_default, size=size, is_double=is_double, nan_check=False)

        kt_default = self.mu1 * self.ka
        kts = array_default_floats(self.kt, kt_default, size=size, is_double=is_double, nan_check=False)

        for pid, u0, f0, ka, kb, kt, mu1, mu2, \
            tmax, mar, trmin in zip_longest(property_ids, u0s, f0s, kas, kbs, kts,
                                            mu1s, self.mu2, tmaxs, mars, trmins):
            #u0 = set_blank_if_default(u0, 0.)
            #f0 = set_blank_if_default(f0, 0.)

            # ka doesn't have a default in MSC 2005r2
            #ka = set_blank_if_default(ka, 1.e8)
            #kb = set_blank_if_default(kb, 1e-14 * ka)
            #kt = set_blank_if_default(kt, mu1 * ka)

            #mu1 = set_blank_if_default(mu1, 0.)
            mu2 = set_blank_if_default(mu2, mu1)
            #tmax = set_blank_if_default(tmax, 0.)
            #mar = set_blank_if_default(mar, 100.)
            #trmin = set_blank_if_default(trmin, 0.001)

            list_fields = ['PGAP', pid, u0, f0, ka, kb, kt, mu1, mu2,
                           tmax, mar, trmin]
            bdf_file.write(print_card(list_fields))
        return
