from __future__ import annotations
from itertools import zip_longest
from typing import Optional, Any, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, integer_or_blank, double_or_blank)
from pyNastran.bdf.cards.elements.bars import set_blank_if_default
#from pyNastran.bdf.cards.properties.bars import _bar_areaL # PBARL as pbarl, A_I1_I2_I12

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property,
    parse_check, searchsorted_filter, save_ifile_comment)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, get_print_card_size # , array_default_int
from .utils import get_density_from_material, basic_mass_material_id

from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.cards.materials import MAT1, MAT4, MAT5

def rod_materials(model: BDF) -> list[MAT1 | MAT4 | MAT5]:
    if model.is_thermal:
        return [model.mat4, model.mat5]
    return [model.mat1]


class CONROD(Element):
    """
    +--------+-----+-----+----+-----+---+---+---+-----+
    |   1    |  2  |  3  |  4 |  5  | 6 | 7 | 8 |  9  |
    +========+=====+=====+====+=====+===+===+===+=====+
    | CONROD | EID | N1  | N2 | MID | A | J | C | NSM |
    +--------+-----+-----+----+-----+---+---+---+-----+
    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 2), dtype='int32')
        self.A = np.array([], dtype='float64')
        self.J = np.array([], dtype='float64')
        self.c = np.array([], dtype='float64')
        self.nsm = np.array([], dtype='float64')

    def add(self, eid: int, mid: int, nodes: list[int],
            A: float=0.0, j: float=0.0, c: float=0.0, nsm: float=0.0,
            ifile: int=0, comment: str='') -> int:
        """
        Creates a CONROD card

        Parameters
        ----------
        eid : int
            element id
        mid : int
            material id
        nids : list[int, int]
            node ids
        A : float; default=0.
            area
        j : float; default=0.
            polar moment of inertia
        c : float; default=0.
            stress factor
        nsm : float; default=0.
            non-structural mass per unit length
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, nodes, mid, A, j, c, nsm, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        nodes = [integer(card, 2, 'n1'),
                 integer(card, 3, 'n2')]
        mid = integer(card, 4, 'mid')
        A = double_or_blank(card, 5, 'A', default=0.0)
        j = double_or_blank(card, 6, 'j', default=0.0)
        c = double_or_blank(card, 7, 'c', default=0.0)
        nsm = double_or_blank(card, 8, 'nsm', default=0.0)
        assert len(card) <= 9, 'len(CONROD card) = %i\ncard=%s' % (len(card), str(card))
        self.cards.append((eid, nodes, mid, A, j, c, nsm, ifile, comment))
        self.n += 1
        return self.n - 1

    def __apply_slice__(self, elem: CONROD, i: np.ndarray) -> None:  # ignore[override]
        elem.ifile = self.ifile[i]
        elem.element_id = self.element_id[i]
        elem.material_id = self.material_id[i]
        elem.nodes = self.nodes[i, :]
        elem.A = self.A[i]
        elem.J = self.J[i]
        elem.c = self.c[i]
        elem.nsm = self.nsm[i]
        elem.n = len(i)

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        material_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        A = np.zeros(ncards, dtype='float64')
        J = np.zeros(ncards, dtype='float64')
        c = np.zeros(ncards, dtype='float64')
        nsm = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, nodesi, mid, Ai, ji, ci, nsmi, ifilei, comment) = card
            ifile[icard] = ifilei
            element_id[icard] = eid
            material_id[icard] = mid
            nodes[icard, :] = nodesi
            A[icard] = Ai
            J[icard] = ji
            c[icard] = ci
            nsm[icard] = nsmi
        self._save(element_id, material_id, nodes, A, J, c, nsm)
        self.cards = []

    def _save(self, element_id, material_id, nodes, A, J, c, nsm,
              ifile=None, comment: Optional[dict[int, str]]=None) -> None:
        ncards = len(element_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.element_id) != 0:
            ifile = np.hstack([self.ifile, ifile])
            element_id = np.hstack([self.element_id, element_id])
            material_id = np.hstack([self.material_id, material_id])
            nodes = np.vstack([self.nodes, nodes])
            A = np.hstack([self.A, A])
            J = np.hstack([self.J, J])
            c = np.hstack([self.c, c])
            nsm = np.hstack([self.nsm, nsm])

        save_ifile_comment(self, ifile, comment)
        nelements = len(element_id)
        self.element_id = element_id
        self.material_id = material_id
        self.nodes = nodes
        self.A = A
        self.J = J
        self.c = c
        self.nsm = nsm
        self.n = nelements

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['element_id'].append(self.element_id)
        used_dict['material_id'].append(self.material_id)
        used_dict['node_id'].append(self.nodes.ravel())

    def convert(self, area_scale: float=1.,
                area_inertia_scale: float=1.0,
                nsm_per_length_scale: float=1.0, **kwargs):
        self.A *= area_scale
        self.J *= area_inertia_scale
        #self.c *= area_scale
        self.nsm *= nsm_per_length_scale

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        materials = [mat.material_id for mat in self.allowed_materials]
        mids = hstack_msg(materials,
                          msg=f'no conrod materials for {self.type}; material_id={self.material_id}')
        mids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   material_id=(mids, self.material_id))

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.material_id.max(), self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        for eid, mid, nodes, A, j, c, nsm in zip_longest(self.element_id, self.material_id, self.nodes,
                                                         self.A, self.J, self.c, self.nsm):
            n1, n2 = nodes
            j = set_blank_if_default(j, 0.0)
            c = set_blank_if_default(c, 0.0)
            nsm = set_blank_if_default(nsm, 0.0)
            list_fields = ['CONROD', eid, n1, n2, mid, A, j, c, nsm]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_materials(self) -> list[MAT1]:
        return [prop for prop in rod_materials(self.model) if prop.n > 0]

    def get_edge_ids(self) -> np.ndarray:
        edge_ids = _1d_edges(self.nodes)
        return edge_ids

    def mass_per_length(self) -> np.ndarray:
        mass_per_length = line_mid_mass_per_length(
            self.material_id, self.nsm, self.A,
            self.allowed_materials)
        return mass_per_length

    def nonstructural_mass(self) -> np.ndarray:
        return self.nsm

    def mass(self) -> np.ndarray:
        mass_per_length = line_mid_mass_per_length(
            self.material_id, self.nsm, self.A,
            self.allowed_materials)
        length = self.length()
        mass = mass_per_length * length
        return mass

    def volume(self) -> np.ndarray:
        length = self.length()
        volume = self.A * length
        return volume

    def length(self) -> np.ndarray:
        length = line_length(self.model, self.nodes)
        return length

    def centroid(self) -> np.ndarray:
        centroid = line_centroid(self.model, self.nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def area(self) -> np.ndarray:
        return self.A


class CROD(Element):
    """
    +------+-----+-----+----+----+
    |   1  |  2  |  3  |  4 |  5 |
    +======+=====+=====+====+====+
    | CROD | EID | PID | N1 | N2 |
    +------+-----+-----+----+----+
    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, pid: int, nodes: list[int],
            ifile: int=0, comment: str='') -> int:
        """
        Creates a CROD card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PROD)
        nids : list[int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nodes, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        nodes = [integer(card, 3, 'n1'),
                 integer(card, 4, 'n2')]
        assert len(card) == 5, 'len(CROD card) = %i\ncard=%s' % (len(card), str(card))
        self.cards.append((eid, pid, nodes, ifile, comment))
        self.n += 1
        return self.n - 1

    def __apply_slice__(self, elem: CROD, i: np.ndarray) -> None:  # ignore[override]
        elem.ifile = self.ifile[i]
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.n = len(i)

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        for icard, card in enumerate(self.cards):
            (eid, pid, nodesi, ifilei, comment) = card
            ifile[icard] = ifilei
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nodesi
        self._save(element_id, property_id, nodes)
        self.cards = []

    def _save(self, element_id, property_id, nodes,
              ifile=None, comment: Optional[dict[int, str]]=None) -> None:
        ncards = len(element_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.element_id) != 0:
            raise NotImplementedError()
        nelements = len(element_id)
        save_ifile_comment(self, ifile, comment)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.n = nelements

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        #used_dict['element_id'].append(self.element_id)
        used_dict['property_id'].append(self.property_id)
        used_dict['node_id'].append(self.nodes.ravel())

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no rod properties for {self.type}; property_id={self.property_id}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        for eid, pid, nodes in zip_longest(self.element_id, self.property_id, self.nodes):
            n1, n2 = nodes
            list_fields = ['CROD', eid, pid, n1, n2]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_properties(self):
        return [prop for prop in [self.model.prod]
                if prop.n > 0]

    def get_edge_ids(self) -> np.ndarray:
        edge_ids = _1d_edges(self.nodes)
        return edge_ids

    def mass_per_length(self) -> np.ndarray:
        mass_per_length = line_pid_mass_per_length(self.property_id, self.allowed_properties)
        return mass_per_length

    def nonstructural_mass(self) -> np.ndarray:
        nsm = line_pid_nsm_per_length(self.property_id, self.allowed_properties)
        return nsm

    def mass_material_id(self) -> np.ndarray:
        material_id = basic_mass_material_id(self.property_id, self.allowed_properties, 'CROD')
        return material_id

    def mass(self) -> np.ndarray:
        mass_per_length = line_pid_mass_per_length(self.property_id, self.allowed_properties)
        length = self.length()
        mass = mass_per_length * length
        return mass

    def area(self) -> np.ndarray:
        area = line_pid_area(self.property_id, self.allowed_properties, self.type)
        return area

    def length(self) -> np.ndarray:
        length = line_length(self.model, self.nodes)
        return length

    def centroid(self) -> np.ndarray:
        centroid = line_centroid(self.model, self.nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def volume(self) -> np.ndarray:
        length = self.length()
        volume = self.area() * length
        return volume


def _1d_edges(nodes: np.ndarray) -> np.ndarray:
    edge_ids = np.column_stack([nodes.min(axis=1), nodes.max(axis=1)])
    assert edge_ids.shape == nodes.shape
    return edge_ids


def tri_edges(nodes: np.ndarray) -> np.ndarray:
    nodes1 = nodes[:, [0, 1]]
    nodes2 = nodes[:, [1, 2]]
    nodes3 = nodes[:, [0, 2]]
    edge_ids1 = np.column_stack([nodes1.min(axis=1), nodes1.max(axis=1)])
    edge_ids2 = np.column_stack([nodes2.min(axis=1), nodes2.max(axis=1)])
    edge_ids3 = np.column_stack([nodes3.min(axis=1), nodes3.max(axis=1)])
    edge_ids = np.vstack([edge_ids1, edge_ids2, edge_ids3])
    return edge_ids

def quad_edges(nodes: np.ndarray) -> np.ndarray:
    nodes1 = nodes[:, [0, 1]]
    nodes2 = nodes[:, [1, 2]]
    nodes3 = nodes[:, [2, 3]]
    nodes4 = nodes[:, [0, 3]]
    edge_ids1 = np.column_stack([nodes1.min(axis=1), nodes1.max(axis=1)])
    edge_ids2 = np.column_stack([nodes2.min(axis=1), nodes2.max(axis=1)])
    edge_ids3 = np.column_stack([nodes3.min(axis=1), nodes3.max(axis=1)])
    edge_ids4 = np.column_stack([nodes4.min(axis=1), nodes4.max(axis=1)])
    edge_ids = np.vstack([edge_ids1, edge_ids2, edge_ids3, edge_ids4])
    return edge_ids


class PROD(Property):
    """
    +------+-----+-----+-----+-----+-----+-----+
    |  1   |  2  |  3  |  4  |  5  |  6  |  7  |
    +======+=====+=====+=====+=====+=====+=====+
    | PROD | PID | MID |  A  |  J  |  C  | NSM |
    +------+-----+-----+-----+-----+-----+-----+
    | PROD |  1  |  2  | 2.0 | 3.0 | 0.5 | 1.0 |
    +------+-----+-----+-----+-----+-----+-----+
    """
    _show_attributes = ['property_id', 'material_id', 'A', 'J', 'c', 'nsm']
    @Property.clear_check
    def clear(self) -> None:
        self.property_id: np.ndarray = np.array([], dtype='int32')
        self.material_id: np.ndarray = np.array([], dtype='int32')
        self.A: np.ndarray = np.array([], dtype='float64')
        self.J: np.ndarray = np.array([], dtype='float64')
        self.c: np.ndarray = np.array([], dtype='float64')
        self.nsm: np.ndarray = np.array([], dtype='float64')

    def add(self, pid: int, mid: int, A: float,
            j: float=0., c: float=0., nsm: float=0.,
            ifile: int=0, comment: str='') -> int:
        """
        Creates a PROD card

        Parameters
        ----------
        pid : int
           property id
        mid : int
           material id
        A : float
           area
        J : float; default=0.
           polar moment of inertia
        c : float; default=0.
           stress factor
        nsm : float; default=0.
           nonstructural mass per unit length
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, mid, A, j, c, nsm, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        A = double(card, 3, 'A')
        j = double_or_blank(card, 4, 'J', default=0.0)
        c = double_or_blank(card, 5, 'c', default=0.0)
        nsm = double_or_blank(card, 6, 'nsm', default=0.0)
        assert len(card) <= 7, f'len(PROD card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, mid, A, j, c, nsm, ifile, comment))
        self.n += 1
        return self.n - 1

    def __apply_slice__(self, prop: PROD, i: np.ndarray) -> None:  # ignore[override]
        prop.ifile = self.ifile[i]
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]
        prop.A = self.A[i]
        prop.J = self.J[i]
        prop.c = self.c[i]
        prop.nsm = self.nsm[i]
        prop.n = len(i)

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype=idtype)
        material_id = np.zeros(ncards, dtype=idtype)
        A = np.zeros(ncards, dtype='float64')
        J = np.zeros(ncards, dtype='float64')
        c = np.zeros(ncards, dtype='float64')
        nsm = np.zeros(ncards, dtype='float64')
        comment = {}

        for icard, card in enumerate(self.cards):
            (pid, mid, Ai, j, ci, nsmi, ifilei, commenti) = card
            ifile[icard] = ifilei
            property_id[icard] = pid
            material_id[icard] = mid
            A[icard] = Ai
            J[icard] = j
            c[icard] = ci
            nsm[icard] = nsmi
        self._save(property_id, material_id, A, J, c, nsm,
                   ifile=ifile)
        self.sort()
        self.cards = []

    def _save(self, property_id, material_id, A, J, c, nsm,
              ifile=None, comment=None):
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.property_id) != 0:
            raise NotImplementedError()
        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.material_id = material_id
        self.A = A
        self.J = J
        self.c = c
        self.nsm = nsm
        self.n = len(ifile)

    def convert(self, area_scale: float=1.,
                area_inertia_scale: float=1.0,
                nsm_per_length_scale: float=1.0, **kwargs):
        self.A *= area_scale
        self.J *= area_inertia_scale
        #self.c *= area_scale
        self.nsm *= nsm_per_length_scale

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          msg=f'no rod materials for {self.type}; material_id={self.material_id}')
        mids.sort()
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.material_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        for pid, mid, A, j, nsm, c in zip_longest(self.property_id, self.material_id,
                                                  self.A, self.J, self.nsm, self.c):
            j = set_blank_if_default(j, 0.0)
            c = set_blank_if_default(c, 0.0)
            nsm = set_blank_if_default(nsm, 0.0)
            list_fields = ['PROD', pid, mid, A, j, c, nsm]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_materials(self) -> list[Any]:
        return [mat for mat in rod_materials(self.model) if mat.n > 0]

    def mass_per_length(self) -> np.ndarray:
        return line_mid_mass_per_length(self.material_id, self.nsm, self.A,
                                        self.allowed_materials)

    def area(self) -> np.ndarray:
        return self.A


class CTUBE(Element):
    """
    +-------+-----+-----+----+----+
    |   1   |  2  |  3  |  4 |  5 |
    +=======+=====+=====+====+====+
    | CTUBE | EID | PID | N1 | N2 |
    +-------+-----+-----+----+----+
    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 2), dtype='int32')

    def add(self, eid: int, pid: int, nids: list[int],
            ifile: int=0, comment: str='') -> int:
        """
        Creates a CTUBE card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id
        nids : list[int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nids, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2')]
        assert len(card) == 5, 'len(CTUBE card) = %i\ncard=%s' % (len(card), str(card))
        self.cards.append((eid, pid, nids, ifile, comment))
        self.n += 1
        return self.n - 1

    def __apply_slice__(self, elem: CTUBE, i: np.ndarray) -> None:  # ignore[override]
        elem.ifile = self.ifile[i]
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]
        elem.n = len(i)

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        ifile = np.zeros(ncards, dtype='int32')
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 2), dtype=idtype)
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, ifilei, comment) = card
            ifile[icard] = ifilei
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes, ifile=ifile)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, nodes,
              ifile=None, comment: Optional[dict[int, str]]=None) -> None:
        ncards = len(element_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.element_id):
            raise NotImplementedError()
        nelements = len(element_id)
        save_ifile_comment(self, ifile, comment)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.n = nelements

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['element_id'].append(self.element_id)
        used_dict['property_id'].append(self.property_id)
        used_dict['node_id'].append(self.nodes.ravel())

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no tube properties for {self.type}; property_id={self.property_id}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        for eid, pid, nodes in zip_longest(self.element_id, self.property_id, self.nodes):
            n1, n2 = nodes
            list_fields = ['CTUBE', eid, pid, n1, n2]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def all_properties(self):
        return [self.model.ptube]

    @property
    def allowed_properties(self):
        all_properties = self.all_properties
        properties = [prop for prop in all_properties
                      if prop.n > 0]
        assert len(properties) > 0, all_properties
        return properties

    def get_edge_ids(self) -> np.ndarray:
        edge_ids = _1d_edges(self.nodes)
        return edge_ids

    def mass_per_length(self) -> np.ndarray:
        mass_per_length = line_pid_mass_per_length(self.property_id, self.allowed_properties)
        return mass_per_length

    def nonstructural_mass(self) -> np.ndarray:
        nsm_per_length = line_pid_nsm_per_length(self.property_id, self.allowed_properties)
        return nsm_per_length

    def mass_material_id(self) -> np.ndarray:
        material_id = basic_mass_material_id(self.property_id, self.allowed_properties, 'CTUBE')
        return material_id

    def mass(self) -> np.ndarray:
        mass_per_length = line_pid_mass_per_length(self.property_id, self.allowed_properties)
        length = self.length()
        mass = mass_per_length * length
        inan = np.isnan(mass)
        if np.any(inan):
            msg = f'{self.type} has nan mass'
            msg += f'element_id={self.element_id[inan]}'
            msg += f'property_id={self.property_id[inan]}\n'
            msg += f'length={length[inan]}\n'
            msg += f'mass_per_length={mass_per_length[inan]}\n'
            self.model.log.error(msg)
        return mass

    def length(self) -> np.ndarray:
        length = line_length(self.model, self.nodes)
        return length

    def centroid(self) -> np.ndarray:
        centroid = line_centroid(self.model, self.nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def area(self) -> np.ndarray:
        area = line_pid_area(self.property_id, self.allowed_properties, self.type)
        return area

    def volume(self) -> np.ndarray:
        length = self.length()
        volume = self.area() * length
        return volume


def line_pid_nsm_per_length(property_id: np.ndarray, allowed_properties: list[Any]):
    nsm_per_length = np.zeros(len(property_id), dtype='float64')
    for prop in allowed_properties:
        i_lookup, i_all = searchsorted_filter(prop.property_id, property_id, msg='')
        if len(i_lookup) == 0:
            continue

        # we're at least using some properties
        nsm_per_lengthi = prop.nsm # nonstructural_mass()
        nsm_per_length[i_lookup] = nsm_per_lengthi[i_all]
    return nsm_per_length

def line_pid_mass_per_length(property_id: np.ndarray, allowed_properties: list[Any]):
    mass_per_length = np.zeros(len(property_id), dtype='float64')
    for prop in allowed_properties:
        i_lookup, i_all = searchsorted_filter(prop.property_id, property_id, msg='')
        if len(i_lookup) == 0:
            continue

        # we're at least using some properties
        mass_per_lengthi = prop.mass_per_length()
        mass_per_length[i_lookup] = mass_per_lengthi[i_all]
    return mass_per_length

def line_pid_area(property_id: np.ndarray, allowed_properties: list[Any], card_type: str):
    area = np.zeros(len(property_id), dtype='float64')
    for prop in allowed_properties:
        try:
            i_lookup, i_all = searchsorted_filter(prop.property_id, property_id, msg=f'{prop.type} property_id')
        except RuntimeError:
            print(allowed_properties)
            raise
        if len(i_lookup) == 0:
            continue

        # we're at least using some properties
        areai = prop.area()
        area[i_lookup] = areai[i_all]
    return area


class PTUBE(Property):
    """
    +-------+------+-----+------+------+------+-----+
    |   1   |  2   |  3  |   4  |  5   |  6   |  7  |
    +=======+======+=====+======+======+======+=====+
    | PTUBE | PID  | MID |  OD  |   T  |  NSM | OD2 |
    +-------+------+-----+------+------+------+-----+
    | PTUBE |  2   |  6  | 6.29 | 0.25 |      |     |
    +-------+------+-----+------+------+------+-----+
    """
    _show_attributes = ['pid', 'mid', 'diameter', 't', 'nsm']
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.diameter = np.zeros((0, 2), dtype='float64')
        self.t = np.array([], dtype='float64')
        self.nsm = np.array([], dtype='float64')

    def add(self, pid: int, mid: int, OD1: float, t: Optional[float]=None,
            nsm: float=0., OD2: Optional[float]=None,
            ifile: int=0, comment: str='') -> int:
        """
        Adds a PTUBE card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        OD1 : float
            outer diameter at End A
        t : float; default=None -> OD1/2.
            thickness
        nsm : float; default=0.
            non-structural mass per unit length
        OD2 : float; default=None -> OD1
            outer diameter at End B
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, mid, OD1, OD2, t, nsm, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        OD1 = double(card, 3, 'OD1')
        t = double_or_blank(card, 4, 't', default=OD1/2)
        nsm = double_or_blank(card, 5, 'nsm', default=0.0)
        OD2 = double_or_blank(card, 6, 'OD2', default=OD1)
        assert len(card) <= 7, f'len(PTUBE card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, mid, OD1, OD2, t, nsm, ifile, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        diameter = np.zeros((ncards, 2), dtype='float64')
        t = np.zeros(ncards, dtype='float64')
        nsm = np.zeros(ncards, dtype='float64')
        comment = {}

        for icard, card in enumerate(self.cards):
            (pid, mid, OD1, OD2, ti, nsmi, ifilei, commenti) = card
            ifile[icard] = ifilei
            property_id[icard] = pid
            material_id[icard] = mid
            diameter[icard] = [OD1, OD2]
            t[icard] = ti
            nsm[icard] = nsmi
        self._save(property_id, material_id, diameter, t, nsm, ifile=ifile)
        self.sort()
        self.cards = []

    def _save(self, property_id, material_id, diameter, t, nsm,
              ifile=None, comment: Optional[dict[int, str]]=None) -> None:
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.property_id):
            raise NotImplementedError()
        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.material_id = material_id
        self.diameter = diameter
        self.t = t
        self.nsm = nsm
        assert self.diameter.ndim == 2, self.diameter.shape
        self.n = len(ifile)

    def __apply_slice__(self, prop: PTUBE, i: np.ndarray) -> None:  # ignore[override]
        assert self.diameter.ndim == 2, self.diameter.shape
        prop.ifile = self.ifile[i]
        prop.property_id = self.property_id[i]
        prop.diameter = self.diameter[i, :]
        prop.t = self.t[i]
        prop.nsm = self.nsm[i]
        prop.material_id = self.material_id[i]
        prop.n = len(i)
        assert prop.diameter.ndim == 2, prop.diameter.shape

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)

    #def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        #property_id = used_dict['property_id']
        #ncards_removed = len(self.property_id) - len(property_id)
        #if ncards_removed:
            #self.slice_card_by_id(property_id, assume_sorted=True, sort_ids=False)
        #return ncards_removed

    def convert(self, xyz_scale: float=1.,
                nsm_per_length_scale: float=1.0, **kwargs):
        self.t *= xyz_scale
        self.diameter *= xyz_scale
        self.nsm *= nsm_per_length_scale

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          msg=f'no rod materials for {self.type}; material_id={self.material_id}')
        mids.sort()
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.material_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        property_id = array_str(self.property_id, size=size)
        material_id = array_str(self.material_id, size=size)

        assert self.diameter.ndim == 2, self.diameter.shape
        for pid, mid, diameter, t, nsm in zip_longest(property_id, material_id,
                                                      self.diameter, self.t, self.nsm):
            OD1, OD2 = diameter
            ts = set_blank_if_default(t, OD1 / 2.)
            nsm = set_blank_if_default(nsm, 0.0)
            OD2s = set_blank_if_default(OD2, OD1)
            list_fields = ['PTUBE', pid, mid, OD1, ts, nsm, OD2s]

            bdf_file.write(print_card(list_fields))
        return

    def area(self) -> np.ndarray:
        t1 = self.t.copy()
        inan_t1 = np.isnan(t1)
        #if self.diameter.ndim == 1:
            #t1[inan_t1] = self.diameter[inan_t1]
            #return _tube_area(t1, self.diameter[inan_t1])

        diameter = self.diameter.copy()
        Dout1 = diameter[:, 0]
        Dout2 = diameter[:, 1]

        # fix D2
        inan_d2 = np.isnan(Dout2)
        Dout2[inan_d2] = Dout1[inan_d2]

        t1[inan_t1] = Dout1[inan_t1] / 2.

        if np.array_equal(Dout1, Dout2):
            A = _tube_area(t1, Dout1)
        else:
            t2 = self.t.copy()
            inan_t2 = np.isnan(t2)
            t2[inan_t2] = Dout2[inan_t2] / 2.
            A = (_tube_area(t1, Dout1) + _tube_area(t2, Dout2)) / 2.
        return A

    def J(self) -> np.ndarray:
        #if self.diameter.ndim == 1:
            #return self._areai(self.diameter)

        diameter = self.diameter.copy()
        Dout1 = diameter[:, 0]
        Dout2 = diameter[:, 1]
        inan = np.isnan(Dout2)
        Dout2[inan] = Dout1[inan]

        if np.array_equal(Dout1, Dout2):
            J = _polar(Dout1, self.t)
        else:
            # J = pi / 16 * (Do^4 - di^4) / Do
            Din1 = Dout1 - 2 * self.t
            Din2 = Dout2 - 2 * self.t
            J1 = np.pi / 16 * (Dout1**4 - Din1**4) / Dout1
            J2 = np.pi / 16 * (Dout2**4 - Din2**4) / Dout2
            #A = (self._areai(Dout1) + self._areai(Dout2)) / 2.
            J = (J1 + J2) / 2
        return J

    @property
    def all_materials(self):
        model = self.model
        if model.is_thermal:
            materials = [self.model.mat4]
        else:
            materials = [self.model.mat1]
        return materials

    @property
    def allowed_materials(self) -> list[Any]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def mass_per_length(self) -> np.ndarray:
        return line_mid_mass_per_length(self.material_id, self.nsm, self.area(),
                                        self.allowed_materials)

def _tube_area(t: np.ndarray, Dout: np.ndarray) -> np.ndarray:
    """Gets the Area of Section 1/2 of the CTUBE."""
    assert len(Dout) == len(t)
    A = np.zeros(len(t), dtype='float64')
    #Dout = self.diameter[:, 0]
    izero = np.where(t == 0.)[0]
    ipos = np.where(t != 0.)[0]
    if len(izero):
        A[izero] = np.pi / 4. * Dout[izero] ** 2
    if len(ipos):
        Din = Dout[ipos] - 2 * t[ipos]
        A[ipos] = np.pi / 4. * (Dout[ipos] * Dout[ipos] - Din[ipos] * Din[ipos])
    return A

def _polar(Dout: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Second polar moment of area
    https://en.wikipedia.org/wiki/Second_polar_moment_of_area

    """
    Din = Dout - 2 * t
    J = np.pi / 32 * (Dout ** 4 - Din ** 4)
    return J

def line_mid_mass_per_length(material_id: np.ndarray,
                             nsm: np.ndarray,
                             area: np.ndarray,
                             allowed_materials: list[Any]) -> np.ndarray:
    """
    calculates the mass per length for:
     - CONROD
     - PROD
     - PBAR
     - PBEAM

    """
    nproperties = len(material_id)
    rho = get_density_from_material(material_id, allowed_materials)
    mass_per_length = rho * area + nsm
    assert len(mass_per_length) == nproperties
    return mass_per_length

def line_vector_length(model: BDF, nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xyz = model.grid.xyz_cid0()
    nid = model.grid.node_id
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    line_vector = xyz2 - xyz1
    length = np.linalg.norm(line_vector, axis=1)
    assert len(length) == nodes.shape[0]
    return line_vector, length

def line_length_nan(model: BDF, nodes: np.ndarray, default_node: int=-1) -> np.ndarray:
    min_node = nodes.min(axis=1)
    assert len(min_node) == len(nodes)
    if min_node.min() == -1:
        length = np.full(min_node.shape, np.nan, dtype='float64')
        if min_node.max() > 0:
            inode = (min_node > 0)
            length[inode] = line_length(model, nodes[inode, :])
    else:
        length = line_length(model, nodes)
    return length

def line_length(model: BDF, nodes: np.ndarray) -> np.ndarray:
    return line_vector_length(model, nodes)[1]

def line_centroid_with_spoints(model: BDF, nodes: np.ndarray) -> np.ndarray:
    grid = model.grid
    xyz = grid.xyz_cid0().copy()
    nid = grid.node_id
    inode = np.searchsorted(nid, nodes)
    izero = (nodes == 0)
    inonzero = ~izero
    inode[izero] = -1

    assert np.array_equal(nid[inode][inonzero], nodes[inonzero])
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    xyz = np.vstack([xyz, [0., 0., 0.]])

    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    centroid = (xyz1 + xyz2) / 2.
    assert centroid.shape[0] == nodes.shape[0]
    return centroid

def line_centroid(model: BDF, nodes: np.ndarray) -> np.ndarray:
    grid = model.grid
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    centroid = (xyz1 + xyz2) / 2.
    assert centroid.shape[0] == nodes.shape[0]
    return centroid

def e_g_nu_from_property_id(property_id: np.ndarray,
                            allowed_properties: list[Property]) -> np.ndarray:
    npid = len(property_id)
    e_g_nu = np.full((npid, 3), np.nan, dtype='float64')
    for prop in allowed_properties:
        i_lookup, i_all = searchsorted_filter(prop.property_id, property_id, msg='')
        if len(i_lookup) == 0:
            continue

        # we're at least using some properties
        e_g_nui = prop.e_g_nu()
        e_g_nu[i_lookup] = e_g_nui[i_all]
    return e_g_nu
