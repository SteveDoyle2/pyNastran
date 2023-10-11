from __future__ import annotations
from typing import Optional, Any, TYPE_CHECKING
import numpy as np
#from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8, print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16 # , print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
#from pyNastran.bdf.cards.elements.bars import set_blank_if_default
#from pyNastran.bdf.cards.elements.solid import volume4
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, integer_or_blank, double_or_blank,
    integer_string_or_blank, string_or_blank)

from pyNastran.dev.bdf_vectorized3.utils import hstack_msg
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check, find_missing
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, make_idim, hslice_by_idim, get_print_card_8_16,
    parse_element_check, parse_property_check, ) # searchsorted_filter,
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    get_print_card, array_str, array_default_str, array_default_int, array_float_nan)
from .utils import get_density_from_material, get_density_from_property, basic_mass_material_id

from .solid_quality import chexa_quality, tetra_quality, penta_quality, pyram_quality, Quality
from .solid_utils import chexa_centroid
from .solid_volume import volume_ctetra, volume_cpenta, volume_cpyram, volume_chexa

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


class SolidElement(Element):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.nnode = 0
        self.nnode_base = 0
        self.nodes: np.ndarray = np.array([[]], dtype='int32')
        self.property_id: np.ndarray = np.array([], dtype='int32')

    def __apply_slice__(self, elem: SolidElement, i: np.ndarray) -> None:  # ignore[override]
        elem.n = len(i)
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.nodes = self.nodes[i, :]

    def add(self, eid: int, pid: int, nids: list[int], comment: str='') -> None:
        #nids2 = [0] * self.nodes.shape[1]
        nids2 = [0 if nid is None else nid
                 for nid in nids]
        nnodes_to_extend = self.nnode - len(nids)
        nids2.extend([0]*nnodes_to_extend)
        self.cards.append((eid, pid, nids2, comment))
        self.n += 1

    def volume(self) -> np.ndarray:
        raise RuntimeError(f'volume is not impemented for {self.type}')

    @property
    def allowed_properties(self) -> list[Any]:
        model = self.model
        allowed_properties = [model.psolid, model.plsolid, model.pcomps, model.pcompls]
        return [prop for prop in allowed_properties if prop.n > 0]

    def mass_breakdown(self) -> np.ndarray:
        """[rho, volume, mass]"""
        #material_id = get_material_from_property(self.property_id, self.allowed_properties)
        rho = get_density_from_property(self.property_id, self.allowed_properties)
        volume = self.volume()
        mass = rho * volume
        breakdown = np.column_stack([rho, volume, mass])
        return breakdown

    def mass(self) -> np.ndarray:
        #material_id = get_material_from_property(self.property_id, self.allowed_properties)
        rho = get_density_from_property(self.property_id, self.allowed_properties)
        if rho.max() == 0. and rho.min() == 0.:
            return np.zeros(len(rho), rho.dtype)
        mass = rho * self.volume()
        return mass

    def mass_material_id(self) -> np.ndarray:
        material_id = basic_mass_material_id(self.property_id, self.allowed_properties, 'CBAR')
        return material_id

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        nid = model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no solid properties for {self.type}')
        pids.sort()

        geom_check(self,
                   missing,
                   node=(nid, self.base_nodes),
                   property_id=(pids, self.property_id))
        geom_check(self,
                   missing,
                   node=(nid, self.midside_nodes), filter_node0=True)

    @property
    def base_nodes(self) -> np.ndarray:  # pragma: no cover
        raise NotImplementedError(f'add base_nodes to {self.type}')
    @property
    def midside_nodes(self) -> np.ndarray:  # pragma: no cover
        raise NotImplementedError(f'add midside_nodes to {self.type}')

    def _save(self,
             element_id: np.ndarray,
             property_id: np.ndarray,
             nodes: np.ndarray) -> None:
        ncards_existing = len(self.element_id)
        midside_nodes = nodes[:, self.nnode_base:]
        if self.model.filter_midside_nodes and midside_nodes.max() == 0:
            nodes = nodes[:, :self.nnode_base]

        if ncards_existing != 0:
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            nodes = np.vstack([self.nodes, nodes])
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.n = len(element_id)

        self.sort()
        self.cards = []

    #def write_8(self) -> str:
        #return self.write(size=8, is_double=False)

    #def write_16(self, is_double: bool=False) -> str:
        #return self.write(size=16, is_double=False)


class CTETRA(SolidElement):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.nodes = np.zeros((0, 10), dtype='int32')
        self.nnode = 10
        self.nnode_base = 4

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'nid1'),
                integer(card, 4, 'nid2'),
                integer(card, 5, 'nid3'),
                integer(card, 6, 'nid4'),
                integer_or_blank(card, 7, 'nid5', default=0),
                integer_or_blank(card, 8, 'nid6', default=0),
                integer_or_blank(card, 9, 'nid7', default=0),
                integer_or_blank(card, 10, 'nid8', default=0),
                integer_or_blank(card, 11, 'nid9', default=0),
                integer_or_blank(card, 12, 'nid10', default=0), ]
        assert len(card) <= 13, f'len(CTETRA card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 10), dtype=idtype)
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes)

    @property
    def base_nodes(self) -> np.ndarray:
        base_nodes = self.nodes[:, :4]
        return base_nodes

    @property
    def midside_nodes(self) -> np.ndarray:
        midside_nodes = self.nodes[:, 4:]
        if midside_nodes.shape[1] == 0:
            return midside_nodes
        assert midside_nodes.shape[1] == 6, midside_nodes.shape
        return midside_nodes

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        max_int = max(self.element_id.max(),
                      self.property_id.max(),
                      self.nodes.max())
        print_card = get_print_card(size, max_int)

        base_nodes = self.base_nodes
        midside_nodes = self.midside_nodes
        assert base_nodes.min() > 0, (base_nodes.min(), base_nodes.max())

        if midside_nodes.shape[1] == 0:
            for eid, pid, nodes in zip(self.element_id, self.property_id, base_nodes):
                msg = print_card(['CTETRA', eid, pid] + nodes.tolist())
                bdf_file.write(msg)
        else:
            for eid, pid, nodes in zip(self.element_id, self.property_id, self.nodes):
                msg = print_card(['CTETRA', eid, pid] + nodes.tolist())
                bdf_file.write(msg)
        return

    def volume(self) -> np.ndarray:
        xyz = self.model.grid.xyz_cid0()
        nid = self.model.grid.node_id
        nodes = self.base_nodes
        nelements = nodes.shape[0]

        inode = np.searchsorted(nid, nodes)
        n1 = xyz[inode[:, 0], :]
        n2 = xyz[inode[:, 1], :]
        n3 = xyz[inode[:, 2], :]
        n4 = xyz[inode[:, 3], :]
        volume = volume_ctetra(n1, n2, n3, n4)
        assert len(volume) == nelements
        #for n1i, n2i, n3i, n4i in zip(n1, n2, n3, n4):
            #v = volume4(n1i, n2i, n3i, n4i)
            #print(v)
        return volume

    def centroid(self) -> np.ndarray:
        xyz = self.model.grid.xyz_cid0()
        nid = self.model.grid.node_id
        nodes = self.base_nodes
        #nelements = nodes.shape[0]

        inode = np.searchsorted(nid, nodes)
        n1 = xyz[inode[:, 0], :]
        n2 = xyz[inode[:, 1], :]
        n3 = xyz[inode[:, 2], :]
        n4 = xyz[inode[:, 3], :]
        centroid = (n1 + n2 + n3 + n4) / 4
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        inode, umissing = find_missing(nid, self.nodes, 'nodes')
        if len(umissing):
            missing['node_id'] = umissing
        #x = 1

    def quality(self) -> Quality:
        out = tetra_quality(self)
        return out

CTETRA_FACE_MAPPER: dict[tuple[int, int], tuple[int, int, int]] = {
}
CPENTA_FACE_MAPPER = {
    # tri faces
    (0, -1): np.array([0, 1, 2]),  # close
    (1, -1): np.array([0, 1, 2]),
    (2, -1): np.array([0, 1, 2]),

    (3, -1): np.array([3, 4, 5]),  # far-reverse
    (4, -1): np.array([3, 4, 5]),
    (5, -1): np.array([3, 4, 5]),

    # reverse points away from the element
    (0, 1): np.array([0, 1, 2]),  # close
    (1, 2): np.array([0, 1, 2]),
    (0, 2): np.array([0, 1, 2]),

    (3, 4): np.array([3, 4, 5]),  # far-reverse
    (4, 5): np.array([3, 4, 5]),
    (3, 5): np.array([3, 4, 5]),

    (0, 4): np.array([0, 1, 4, 3]),  # bottom
    (1, 3): np.array([0, 1, 4, 3]),
    (4, 0): np.array([0, 1, 4, 3]),  # bottom (flipped)
    (3, 1): np.array([0, 1, 4, 3]),

    (0, 5): np.array([0, 2, 5, 3]),  # left-reverse
    (2, 3): np.array([0, 2, 5, 3]),
    (5, 0): np.array([0, 2, 5, 3]),  # left-reverse (flipped)
    (3, 2): np.array([0, 2, 5, 3]),

    (1, 5): np.array([1, 4, 5, 2]),  # right
    (2, 4): np.array([1, 4, 5, 2]),
    (5, 1): np.array([1, 4, 5, 2]),  # right (flipped)
    (4, 2): np.array([1, 4, 5, 2]),
}
CHEXA_FACE_MAPPER = {
    (7, 5) : np.array([7, 6, 5, 4]),
    (5, 7) : np.array([7, 6, 5, 4]),
    (6, 4) : np.array([7, 6, 5, 4]),
    (4, 6) : np.array([7, 6, 5, 4]),

    (0, 2) : np.array([0, 1, 2, 3]),
    (2, 0) : np.array([0, 1, 2, 3]),
    (1, 3) : np.array([0, 1, 2, 3]),
    (3, 1) : np.array([0, 1, 2, 3]),

    (0, 7) : np.array([0, 3, 7, 4]),
    (7, 0) : np.array([0, 3, 7, 4]),
    (3, 4) : np.array([0, 3, 7, 4]),
    (4, 3) : np.array([0, 3, 7, 4]),

    (5, 2) : np.array([5, 6, 2, 1]),
    (2, 5) : np.array([5, 6, 2, 1]),
    (6, 1) : np.array([5, 6, 2, 1]),
    (1, 6) : np.array([5, 6, 2, 1]),

    (4, 1) : np.array([4, 5, 1, 0]),
    (1, 4) : np.array([4, 5, 1, 0]),
    (5, 0) : np.array([4, 5, 1, 0]),
    (0, 5) : np.array([4, 5, 1, 0]),

    (2, 7) : np.array([2, 6, 7, 3]),
    (7, 2) : np.array([2, 6, 7, 3]),
    (6, 3) : np.array([2, 6, 7, 3]),
    (3, 6) : np.array([2, 6, 7, 3]),
}

class CPENTA(SolidElement):
    r"""
    +--------+-----+-----+----+----+----+----+----+----+
    |    1   |  2  |  3  |  4 |  5 |  6 |  7 |  8 |  9 |
    +========+=====+=====+====+====+====+====+====+====+
    | CPENTA | EID | PID | G1 | G2 | G3 | G4 | G5 | G6 |
    +--------+-----+-----+----+----+----+----+----+----+

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.nodes = np.zeros((0, 15), dtype='int32')
        self.nnode = 15
        self.nnode_base = 6

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'),
            integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'),
            integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'),
            integer(card, 8, 'nid6'),
            integer_or_blank(card, 9, 'nid7', default=0),
            integer_or_blank(card, 10, 'nid8', default=0),
            integer_or_blank(card, 11, 'nid9', default=0),
            integer_or_blank(card, 12, 'nid10', default=0),
            integer_or_blank(card, 13, 'nid11', default=0),
            integer_or_blank(card, 14, 'nid12', default=0),
            integer_or_blank(card, 15, 'nid13', default=0),
            integer_or_blank(card, 16, 'nid14', default=0),
            integer_or_blank(card, 17, 'nid15', default=0),
        ]
        assert len(card) <= 18, f'len(CPENTA card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 15), dtype=idtype)
        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes)

    @property
    def base_nodes(self) -> np.ndarray:
        base_nodes = self.nodes[:, :6]
        return base_nodes

    @property
    def midside_nodes(self) -> np.ndarray:
        midside_nodes = self.nodes[:, 6:]
        if midside_nodes.shape[1] == 0:
            return midside_nodes
        assert midside_nodes.shape[1] == 9, midside_nodes.shape
        return midside_nodes

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        max_int = max(self.element_id.max(),
                      self.property_id.max(),
                      self.nodes.max())
        print_card = get_print_card(size, max_int)

        base_nodes = self.base_nodes
        midside_nodes = self.midside_nodes
        assert base_nodes.min() > 0, (base_nodes.min(), base_nodes.max())

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        #card_name = 'CPENTA' if size == 8 else 'CPENTA*'
        base_nodes_str = array_str(base_nodes, size=size)
        if midside_nodes.shape[1] == 0:
            #if size == 8:
            for eid, pid, nodes in zip(element_ids, property_ids, base_nodes_str):
                #msg1 = (
                    #f'{card_name:<8s}{eid:>8s}{pid:>8s}{nodes[0]:>8s}'
                    #f'{nodes[1]:>8s}{nodes[2]:>8s}{nodes[3]:>8s}{nodes[4]:>8s}{nodes[5]:>8s}\n')
                msg2 = print_card(['CPENTA', eid, pid] + nodes.tolist())
                #assert msg1 == msg2
                bdf_file.write(msg2)
            #else:
                #for eid, pid, nodes in zip(element_ids, property_ids, base_nodes_str):
                    #msg1 = (
                        #f'{card_name:<8s}{eid:>16s}{pid:>16s}{nodes[0]:>16s}{nodes[1]:>16s}\n'
                        #f'{"*":<8s}{nodes[2]:>16s}{nodes[3]:>16s}{nodes[4]:>16s}{nodes[5]:>16s}\n')
                    #msg2 = print_card(['CPENTA', eid, pid] + nodes.tolist())
                    #assert msg1 == msg2
                    #lines.append(msg2)
        else:
            #midside_nodes_str = array_default_int(midside_nodes, default=0, size=size)
            nodes_str = array_default_int(self.nodes, default=0, size=size)
            for eid, pid, nodes in zip(element_ids, property_ids, nodes_str):
                msg = print_card(['CPENTA', eid, pid] + nodes.tolist())
                bdf_file.write(msg)
        return

    def volume(self):
        xyz = self.model.grid.xyz_cid0()
        nid = self.model.grid.node_id
        nodes = self.base_nodes

        inode = np.searchsorted(nid, nodes)
        n1 = xyz[inode[:, 0], :]
        n2 = xyz[inode[:, 1], :]
        n3 = xyz[inode[:, 2], :]
        n4 = xyz[inode[:, 3], :]
        n5 = xyz[inode[:, 4], :]
        n6 = xyz[inode[:, 5], :]
        volume = volume_cpenta(n1, n2, n3, n4, n5, n6)
        return volume

    def centroid(self) -> np.ndarray:
        xyz = self.model.grid.xyz_cid0()
        nid = self.model.grid.node_id
        nodes = self.base_nodes
        #nelements = nodes.shape[0]

        inode = np.searchsorted(nid, nodes)
        n1 = xyz[inode[:, 0], :]
        n2 = xyz[inode[:, 1], :]
        n3 = xyz[inode[:, 2], :]
        n4 = xyz[inode[:, 3], :]
        n5 = xyz[inode[:, 4], :]
        n6 = xyz[inode[:, 5], :]
        centroid = (n1 + n2 + n3 + n4 + n5 + n6) / 6
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def quality(self) -> Quality:
        out = penta_quality(self)
        return out


class CPYRAM(SolidElement):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.nodes = np.zeros((0, 13), dtype='int32')
        self.nnode = 13
        self.nnode_base = 5

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'),
            integer_or_blank(card, 8, 'nid6', default=0),
            integer_or_blank(card, 9, 'nid7', default=0),
            integer_or_blank(card, 10, 'nid8', default=0),
            integer_or_blank(card, 11, 'nid9', default=0),
            integer_or_blank(card, 12, 'nid10', default=0),
            integer_or_blank(card, 13, 'nid11', default=0),
            integer_or_blank(card, 14, 'nid12', default=0),
            integer_or_blank(card, 15, 'nid13', default=0)
        ]
        assert len(card) <= 16, f'len(CPYRAM13 1card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 13), dtype=idtype)
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes)

    @property
    def base_nodes(self) -> np.ndarray:
        base_nodes = self.nodes[:, :5]
        return base_nodes

    @property
    def midside_nodes(self) -> np.ndarray:
        midside_nodes = self.nodes[:, 5:]
        if midside_nodes.shape[1] == 0:
            return midside_nodes
        assert midside_nodes.shape[1] == 8, midside_nodes.shape
        return midside_nodes

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        max_int = max(self.element_id.max(),
                      self.property_id.max(),
                      self.nodes.max())
        print_card = get_print_card(size, max_int)

        base_nodes = self.base_nodes
        midside_nodes = self.midside_nodes
        assert midside_nodes.shape[1] in [0, 8], midside_nodes.shape
        assert base_nodes.min() > 0, (base_nodes.min(), base_nodes.max())

        if midside_nodes.shape[1] == 0:
            for eid, pid, nodes in zip(self.element_id, self.property_id, base_nodes):
                msg = print_card(['CPYRAM', eid, pid] + nodes.tolist())
                bdf_file.write(msg)
        else:
            for eid, pid, nodes in zip(self.element_id, self.property_id, self.nodes):
                msg = print_card(['CPYRAM', eid, pid] + nodes.tolist())
                bdf_file.write(msg)
        return

    def volume(self):
        xyz = self.model.grid.xyz_cid0()
        nid = self.model.grid.node_id
        nodes = self.base_nodes
        #nelements = nodes.shape[0]

        inode = np.searchsorted(nid, nodes)
        n1 = xyz[inode[:, 0], :]
        n2 = xyz[inode[:, 1], :]
        n3 = xyz[inode[:, 2], :]
        n4 = xyz[inode[:, 3], :]
        n5 = xyz[inode[:, 4], :]
        volume = volume_cpyram(n1, n2, n3, n4, n5)
        return volume

    def centroid(self) -> np.ndarray:
        xyz = self.model.grid.xyz_cid0()
        nid = self.model.grid.node_id
        nodes = self.base_nodes
        #nelements = nodes.shape[0]

        inode = np.searchsorted(nid, nodes)
        n1 = xyz[inode[:, 0], :]
        n2 = xyz[inode[:, 1], :]
        n3 = xyz[inode[:, 2], :]
        n4 = xyz[inode[:, 3], :]
        n5 = xyz[inode[:, 4], :]
        centroid = (n1 + n2 + n3 + n4 + n5) / 5
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def quality(self) -> Quality:
        out = pyram_quality(self)
        return out


class CHEXA(SolidElement):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.nodes = np.zeros((0, 20), dtype='int32')
        self.nnode_base = 8
        self.nnode = 20

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'), integer(card, 8, 'nid6'),
            integer(card, 9, 'nid7'), integer(card, 10, 'nid8'),
            integer_or_blank(card, 11, 'nid9', default=0),
            integer_or_blank(card, 12, 'nid10', default=0),
            integer_or_blank(card, 13, 'nid11', default=0),
            integer_or_blank(card, 14, 'nid12', default=0),
            integer_or_blank(card, 15, 'nid13', default=0),
            integer_or_blank(card, 16, 'nid14', default=0),
            integer_or_blank(card, 17, 'nid15', default=0),
            integer_or_blank(card, 18, 'nid16', default=0),
            integer_or_blank(card, 19, 'nid17', default=0),
            integer_or_blank(card, 20, 'nid18', default=0),
            integer_or_blank(card, 21, 'nid19', default=0),
            integer_or_blank(card, 22, 'nid20', default=0),
        ]
        assert len(card) <= 23, f'len(CHEXA20 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 20), dtype=idtype)
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes)

    @property
    def base_nodes(self) -> np.ndarray:
        base_nodes = self.nodes[:, :8]
        return base_nodes

    @property
    def midside_nodes(self) -> np.ndarray:
        midside_nodes = self.nodes[:, 8:]
        if midside_nodes.shape[1] == 0:
            return midside_nodes
        assert midside_nodes.shape[1] == 12, midside_nodes.shape
        return midside_nodes

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        max_int = max(self.element_id.max(),
                      self.property_id.max(),
                      self.nodes.max())
        print_card = get_print_card(size, max_int)

        base_nodes = self.base_nodes
        midside_nodes = self.midside_nodes
        assert base_nodes.min() > 0, (base_nodes.min(), base_nodes.max())

        if midside_nodes.shape[1] == 0:
            for eid, pid, nodes in zip(self.element_id, self.property_id, base_nodes):
                msg = print_card(['CHEXA', eid, pid] + nodes.tolist())
                bdf_file.write(msg)
        else:
            for eid, pid, nodes in zip(self.element_id, self.property_id, self.nodes):
                msg = print_card(['CHEXA', eid, pid] + nodes.tolist())
                bdf_file.write(msg)
        return

    def centroid(self) -> np.ndarray:
        centroid = chexa_centroid(self)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def volume(self) -> np.ndarray:
        xyz = self.model.grid.xyz_cid0()
        nid = self.model.grid.node_id
        nodes = self.base_nodes

        inode = np.searchsorted(nid, nodes)
        n1 = xyz[inode[:, 0], :]
        n2 = xyz[inode[:, 1], :]
        n3 = xyz[inode[:, 2], :]
        n4 = xyz[inode[:, 3], :]
        n5 = xyz[inode[:, 4], :]
        n6 = xyz[inode[:, 5], :]
        n7 = xyz[inode[:, 6], :]
        n8 = xyz[inode[:, 7], :]
        volume = volume_chexa(n1, n2, n3, n4, n5, n6, n7, n8)
        return volume

    def quality(self) -> Quality:
        out = chexa_quality(self)
        return out


class PSOLID(Property):
    """
    +--------+-----+-----+-------+-----+--------+---------+------+
    |    1   |  2  |  3  |   4   |  5  |    6   |    7    |   8  |
    +========+=====+=====+=======+=====+========+=========+======+
    | PSOLID | PID | MID | CORDM | IN  | STRESS |   ISOP  | FCTN |
    +--------+-----+-----+-------+-----+--------+---------+------+
    | PSOLID |  1  |     |   1   | 0   |        |         |      |
    +--------+-----+-----+-------+-----+--------+---------+------+
    | PSOLID |  2  | 100 |   6   | TWO |  GRID  | REDUCED |      |
    +--------+-----+-----+-------+-----+--------+---------+------+
    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.material_id = np.array([], dtype='int32')

    def slice_card_by_property_id(self, property_id: np.ndarray) -> PSOLID:
        """uses a node_ids to extract PSOLIDs"""
        iprop = self.index(property_id)
        prop = self.slice_card_by_index(iprop)
        return prop

    def __apply_slice__(self, prop: PSOLID, i: np.ndarray) -> None:  # ignore[override]
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]
        prop.coord_id = self.coord_id[i]
        prop.integ = self.integ[i]
        prop.stress = self.stress[i]
        prop.isop = self.isop[i]
        prop.fctn = self.fctn[i]
        prop.n = len(i)

    def add(self, pid: int, mid: int, cordm: int=0,
            integ: Optional[str|int]=None,
            stress: Optional[str]='',
            isop: Optional[str]='',
            fctn: str='SMECH', comment='') -> None:
        if integ is None:
            integ = ''
        if stress is None:
            stress = ''
        if isop is None:
            isop = ''
        self.cards.append((pid, mid, cordm, integ, stress, isop, fctn, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        cordm = integer_or_blank(card, 3, 'cordm', default=-1)
        integ = integer_string_or_blank(card, 4, 'integ', default='')
        stress = integer_string_or_blank(card, 5, 'stress', default='')
        isop = integer_string_or_blank(card, 6, 'isop', default='')
        fctn = string_or_blank(card, 7, 'fctn', default='SMECH')
        assert len(card) <= 8, f'len(PSOLID card) = {len(card):d}\ncard={card}'

        self.cards.append((pid, mid, cordm, integ, stress, isop, fctn, comment))
        self.n += 1
        return self.n

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        coord_id = np.zeros(ncards, dtype='int32')

        integ = np.full(ncards, '', dtype='|U8')
        stress = np.full(ncards, '', dtype='|U8')
        isop = np.full(ncards, '', dtype='|U8')
        fctn = np.full(ncards, '', dtype='|U8')

        for icard, card in enumerate(self.cards):
            (pid, mid, cordm, integi, stressi, isopi, fctni, comment) = card
            if integi == 0:
                integi = '0'
            elif integi == 1:
                integi = '1'
            elif integi == 2:
                integi = 'TWO'
            elif integi == 3:
                integi = 'THREE'

            if stressi == 0:
                stressi = '0' # TODO: not 100%
            elif stressi == 1:
                stressi = 'GAUSS'

            if isopi == 0:
                isopi = '0'
            if isopi == 1:
                isopi = 'FULL'

            property_id[icard] = pid
            material_id[icard] = mid
            coord_id[icard] = cordm

            assert isinstance(integi, str), integi
            assert isinstance(stressi, str), stressi
            assert isinstance(isopi, str), isopi
            assert isinstance(fctni, str), fctni
            integ[icard] = integi
            stress[icard] = stressi
            isop[icard] = isopi
            fctn[icard] = fctni
        self._save(property_id, material_id, coord_id, integ, stress, isop, fctn)
        self.sort()
        self.cards = []

    def _save(self, property_id, material_id, coord_id, integ, stress, isop, fctn):
        if len(self.property_id) != 0:
            property_id = np.hstack([self.property_id, property_id])
            material_id = np.hstack([self.material_id, material_id])
            coord_id = np.hstack([self.coord_id, coord_id])
            integ = np.hstack([self.integ, integ])
            stress = np.hstack([self.stress, stress])
            isop = np.hstack([self.isop, isop])
            fctn = np.hstack([self.fctn, fctn])

        nproperties = len(property_id)
        self.property_id = property_id
        self.material_id = material_id
        self.coord_id = coord_id

        #assert isinstance(integ, str), integ
        #assert isinstance(stress, str), stress
        #assert isinstance(isop, str), isop
        #assert isinstance(fctn, str), fctn
        self.integ = integ
        self.stress = stress
        self.isop = isop
        self.fctn = fctn
        self.n = nproperties

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        cid = model.coord.coord_id
        mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          msg=f'no solid materials for {self.type}')
        ucoord_id = np.unique(self.coord_id)
        ucoord_id = np.setdiff1d(ucoord_id, [-1])
        mids.sort()
        geom_check(self,
                   missing,
                   coord=(cid, ucoord_id),
                   material_id=(mids, self.material_id))

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        max_int = max(
            self.property_id.max(),
            self.material_id.max(),
            self.coord_id.max())
        print_card = get_print_card(size, max_int)

        property_id = array_str(self.property_id, size=size)
        material_id = array_str(self.material_id, size=size)
        coord_id = array_default_int(self.coord_id, default=-1, size=size)
        for pid, mid, cordm, integ, stress, isop, fctn in zip(property_id, material_id, coord_id,
                                                              self.integ, self.stress, self.isop, self.fctn):
            fctn = set_blank_if_default(fctn, 'SMECH')
            fields = ['PSOLID', pid, mid, cordm, integ, stress, isop, fctn]
            bdf_file.write(print_card(fields))
        return

    @property
    def allowed_materials(self) -> list[Any]:
        """TODO: what about SOL 200 or undefined case control decks?"""
        #MAT1, MAT3, MAT4, MAT5, MAT9, MAT10,
        #MAT10C, MATF10C, MAT11, or MATPOR
        model = self.model
        if model.is_thermal:
            all_materials = [model.mat4, model.mat5]
        else:
            all_materials = [model.mat1, model.mat9, model.mat10, model.mat11]
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def rho(self) -> np.ndarray:
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        return rho


class PLSOLID(Property):
    """
    Defines a fully nonlinear (i.e., large strain and large rotation)
    hyperelastic solid element.

    +---------+-----+-----+-----+
    |    1    |  2  |  3  |  4  |
    +=========+=====+=====+=====+
    | PLSOLID | PID | MID | STR |
    +---------+-----+-----+-----+
    | PLSOLID |  20 |  21 |     |
    +---------+-----+-----+-----+
    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.material_id = np.array([], dtype='int32')

    def __apply_slice__(self, prop: PLSOLID, i: np.ndarray) -> None:  # ignore[override]
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]
        prop.stress_strain = self.stress_strain[i]

    def add(self, pid, mid, stress_strain='GRID', ge=0., comment='') -> None:
        self.cards.append((pid, mid, stress_strain, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PLSOLID card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        stress_strain = string_or_blank(card, 3, 'stress_strain', 'GRID')
        assert len(card) <= 4, f'len(PLSOLID card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, mid, stress_strain, comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        self.property_id = np.zeros(ncards, dtype='int32')
        self.material_id = np.zeros(ncards, dtype='int32')
        self.stress_strain = np.zeros(ncards, dtype='|U4')

        for icard, card in enumerate(self.cards):
            (pid, mid, stress_strain, comment) = card
            assert isinstance(stress_strain, str), stress_strain

            self.property_id[icard] = pid
            self.material_id[icard] = mid
            self.stress_strain[icard] = stress_strain
        self.sort()
        self.cards = []

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          msg=f'no solid materials for {self.type}')
        mids.sort()
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_ids = array_str(self.property_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        stress_strains = array_default_str(self.stress_strain, default='GRID', size=size)
        for pid, mid, stress_strain in zip(property_ids, material_ids, stress_strains):
            fields = ['PLSOLID', pid, mid, stress_strain]
            bdf_file.write(print_card(fields))
        return

    @property
    def all_materials(self) -> list[Any]:
        model = self.model
        #if model.is_thermal:
            #materials = []
            #return [mat for mat in [model.mat4, model.mat5] if mat.n > 0]
        #else:
            #materials = []
            #return [mat for mat in [model.mat1, model.mat9, model.mat10, model.mat11] if mat.n > 0]

        #Identification number of a MATHP or MATHE entry.
        materials = [model.mathp, model.mathe]

        #3. For SOL 106, MID must refer to a MATHP entry.
        return materials

    @property
    def allowed_materials(self) -> list[Any]:
        """TODO: what about SOL 200 or undefined case control decks?"""
        #MAT1, MAT3, MAT4, MAT5, MAT9, MAT10,
        #MAT10C, MATF10C, MAT11, or MATPOR
        model = self.model
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={model.materials}'
        return materials

    def rho(self) -> np.ndarray:
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        return rho


class PCOMPS(Property):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.coord_id = np.array([], dtype='int32')
        self.psdir = np.array([], dtype='int32')
        self.sb = np.array([], dtype='float64')
        self.nb = np.array([], dtype='float64')
        self.tref = np.array([], dtype='float64')
        self.ge = np.array([], dtype='float64')

        self.global_ply_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.thickness = np.array([], dtype='float64')
        self.theta = np.array([], dtype='float64')

        self.failure_theory = np.array([], dtype='|U8')
        self.interlaminar_failure_theory = np.array([], dtype='|U8')
        self.sout = np.array([], dtype='|U8')

    def slice_card_by_property_id(self, property_id: np.ndarray) -> PCOMPS:
        """uses a node_ids to extract PCOMPSs"""
        iprop = self.index(property_id)
        prop = self.slice_card_by_index(iprop)
        return prop

    def __apply_slice__(self, prop: PCOMPS, i: np.ndarray) -> None:  # ignore[override]
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        #prop.material_id = self.material_id[i]

        prop.coord_id = self.coord_id[i]
        prop.psdir = self.psdir[i]
        prop.sb = self.sb[i]
        prop.nb = self.nb[i]
        prop.tref = self.tref[i]
        prop.ge = self.ge[i]

        #self.property_id = property_id
        #self.global_ply_id = global_ply_id
        #self.material_id = material_id
        #self.thickness = thickness
        #self.theta = theta
        #self.failure_theory = failure_theory
        #self.interlaminar_failure_theory = interlaminar_failure_theory
        #self.sout = sout

        iply = self.iply
        prop.global_ply_id = hslice_by_idim(i, iply, self.global_ply_id)
        prop.material_id = hslice_by_idim(i, iply, self.material_id)
        prop.thickness = hslice_by_idim(i, iply, self.thickness)
        prop.theta = hslice_by_idim(i, iply, self.theta)
        prop.failure_theory = hslice_by_idim(i, iply, self.failure_theory)
        prop.interlaminar_failure_theory = hslice_by_idim(i, iply, self.interlaminar_failure_theory)
        prop.sout = hslice_by_idim(i, iply, self.sout)

        prop.nply = self.nply[i]

    def add(self, pid: int,
            global_ply_ids: list[int],
            mids: list[int],
            thicknesses: list[float],
            thetas: list[float],
            cordm: int=0, psdir: int=13,
            sb: Optional[float]=None,
            nb: Optional[float]=None,
            tref: float=0.0, ge: float=0.0,
            failure_theories=None, interlaminar_failure_theories=None,
            souts=None, comment='') -> None:
        nb = nb if nb is not None else np.nan
        sb = sb if sb is not None else np.nan

        nplies = len(mids)
        if failure_theories is None:
            failure_theories = [''] * nplies
        if interlaminar_failure_theories is None:
            interlaminar_failure_theories = [''] * nplies
        if souts is None:
            souts = ['NO'] * nplies

        self.cards.append((pid, cordm, psdir, sb, nb, tref, ge,
                           global_ply_ids, mids, thicknesses, thetas,
                           failure_theories, interlaminar_failure_theories, souts))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PCOMPS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        cordm = integer_or_blank(card, 2, 'cordm', default=0)
        psdir = integer_or_blank(card, 3, 'psdir', default=13)
        sb = double_or_blank(card, 4, 'sb', default=np.nan)
        nb = double_or_blank(card, 5, 'nb', default=np.nan)
        tref = double_or_blank(card, 6, 'tref', default=0.0)
        ge = double_or_blank(card, 7, 'ge', default=0.0)

        nfields = len(card) - 1
        #nrows =
        ifield = 9
        global_ply_ids = []
        mids = []
        thicknesses = []
        thetas = []
        failure_theories = []
        interlaminar_failure_theories = []
        souts = []
        iply = 1
        while ifield < nfields:
            global_ply_id = integer(card, ifield, 'global_ply_id_%d' % iply)
            mid = integer(card, ifield + 1, 'mid_%d' % iply)
            t = double(card, ifield + 2, 'thickness_%d' % iply)
            theta = double(card, ifield + 3, 'theta_%d' % iply)
            ft = string_or_blank(card, ifield + 4, 'failure_theory_%d' % iply, default='')
            ift = string_or_blank(card, ifield + 5, 'interlaminar_failure_theory_%d' % iply, default='')
            sout = string_or_blank(card, ifield + 6, 'sout_%d' % iply, default='NO')
            global_ply_ids.append(global_ply_id)
            mids.append(mid)
            thicknesses.append(t)
            thetas.append(theta)
            failure_theories.append(ft)
            interlaminar_failure_theories.append(ift)
            assert ft in {'', 'PFA', 'HOFF', 'HILL'}, f'PCOMPS pid={pid} failure_theory={ft!r}'
            assert ift in {'', 'SB'}, f'PCOMPS pid={pid} interlaminar_failure_theory={ift!r}'
            souts.append(sout)
            iply += 1
            ifield += 8
        assert len(card) <= ifield, f'len(PCOMPS card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, cordm, psdir, sb, nb, tref, ge,
                           global_ply_ids, mids, thicknesses, thetas,
                           failure_theories, interlaminar_failure_theories, souts))
        self.n += 1
        return self.n

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        coord_id = np.zeros(ncards, dtype='int32')
        psdir = np.zeros(ncards, dtype='int32')
        sb = np.zeros(ncards, dtype='float64')
        nb = np.zeros(ncards, dtype='float64')
        tref = np.zeros(ncards, dtype='float64')
        ge = np.zeros(ncards, dtype='float64')

        nply = np.zeros(ncards, dtype='int32')
        all_global_ply_ids = []
        all_mids = []
        all_thicknesses = []
        all_thetas = []
        all_failure_theories = []
        all_interlaminar_failure_theories = []
        all_souts = []

        for icard, card in enumerate(self.cards):
            (pid, cordm, psdiri, sbi, nbi, trefi, gei,
             global_ply_ids, mids, thicknesses, thetas,
             failure_theories, interlaminar_failure_theories, souts) = card
            property_id[icard] = pid
            coord_id[icard] = cordm
            psdir[icard] = psdiri
            sb[icard] = sbi
            nb[icard] = nbi
            tref[icard] = trefi
            ge[icard] = gei
            nply[icard] = len(global_ply_ids)

            all_global_ply_ids.extend(global_ply_ids)
            all_mids.extend(mids)
            all_thicknesses.extend(thicknesses)
            all_thetas.extend(thetas)
            all_failure_theories.extend(failure_theories)
            all_interlaminar_failure_theories.extend(interlaminar_failure_theories)
            all_souts.extend(souts)
            for souti in souts:
                assert len(souti) <= 3, souti

        global_ply_id = np.array(all_global_ply_ids, dtype='int32')
        material_id = np.array(all_mids, dtype='int32')
        thickness = np.array(all_thicknesses, dtype='float64')
        theta = np.array(all_thetas, dtype='float64')
        failure_theory = np.array(all_failure_theories, dtype='|U8')
        interlaminar_failure_theory = np.array(all_interlaminar_failure_theories, dtype='|U8')
        sout = np.array(all_souts, dtype='|U3')

        self._save(property_id, global_ply_id, material_id, thickness, theta, nply,
                   failure_theory, interlaminar_failure_theory, sout,
                   sb, nb, psdir, tref, ge, coord_id)

        self.sort()
        self.cards = []

    def _save(self, property_id, global_ply_id, material_id, thickness, theta, nply,
              failure_theory, interlaminar_failure_theory, sout,
              sb, nb, psdir, tref, ge, coord_id):
        if len(self.property_id) != 0:
            raise NotImplementedError()
        self.property_id = property_id
        self.global_ply_id = global_ply_id
        self.material_id = material_id
        self.thickness = thickness
        self.theta = theta
        self.nply = nply

        self.failure_theory = failure_theory
        self.interlaminar_failure_theory = interlaminar_failure_theory
        self.sout = sout

        self.sb = sb
        self.nb = nb
        self.psdir = psdir
        self.tref = tref
        self.ge = ge
        self.coord_id = coord_id

    @property
    def iply(self) -> np.ndarray:
        return make_idim(self.n, self.nply)

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        for pid, cordm, psdir, sb, nb, tref, ge, iply in zip(self.property_id, self.coord_id,
                                                             self.psdir, self.sb, self.nb, self.tref, self.ge, self.iply):
            iply0, iply1 = iply

            if np.isnan(sb):
                sb = None
            if np.isnan(nb):
                nb = None
            list_fields = ['PCOMPS', pid, cordm, psdir, sb,
                           nb, tref, ge, None]

            global_ply_ids = self.global_ply_id[iply0:iply1]
            mids = self.material_id[iply0:iply1]
            thicknesses = self.thickness[iply0:iply1]
            thetas = self.theta[iply0:iply1]
            failure_theories = self.failure_theory[iply0:iply1]
            interlaminar_failure_theories = self.interlaminar_failure_theory[iply0:iply1]
            souts = self.sout[iply0:iply1]
            for glply, mid, t, theta, ft, ift, sout in zip(global_ply_ids,
                                                           mids, thicknesses, thetas,
                                                           failure_theories,
                                                           interlaminar_failure_theories,
                                                           souts):
                list_fields += [glply, mid, t, theta, ft, ift, sout, None]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_materials(self) -> list[Any]:
        """TODO: what about SOL 200 or undefined case control decks?"""
        model = self.model
        if model.is_thermal:
            all_materials = [model.mat4]
        else:
            all_materials = [model.mat1, model.mat10]
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def rho(self) -> np.ndarray:
        rho = get_density_from_material(self.material_id, self.allowed_materials)
        return rho


class PCOMPLS(Property):
    """
    +---------+------+--------+--------+--------+--------+
    | PCOMPLS |  PID | DIRECT |  CORDM |   SB   |  ANAL  |
    +---------+------+--------+--------+--------+--------+
    |         |  C8  |  BEH8  |  INT8  | BEH8H  | INT8H  |
    +---------+------+--------+--------+--------+--------+
    |         |  C20 |  BEH20 |  INT20 | BEH20H | INT20H |
    +---------+------+--------+--------+--------+--------+
    |         |  ID1 |  MID1  |   T1   | THETA1 |        |
    +---------+------+--------+--------+--------+--------+
    |         |  ID2 |  MID2  |   T2   | THETA2 |        |
    +---------+------+--------+--------+--------+--------+

    +---------+------+--------+--------+--------+--------+
    | PCOMPLS | 782  | 1      |        |        |        |
    +---------+------+--------+--------+--------+--------+
    |         | 1001 |   171  |   .3   |  12.3  |        |
    +---------+------+--------+--------+--------+--------+
    |         | 100  |   175  |   .7   |  77.7  |        |
    +---------+------+--------+--------+--------+--------+
    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.material_id = np.array([], dtype='int32')

    def slice_card_by_property_id(self, property_id: np.ndarray) -> PCOMPLS:
        """uses a node_ids to extract PCOMPLSs"""
        iprop = self.index(property_id)
        prop = self.slice_card_by_index(iprop)
        return prop

    #def __apply_slice__(self, prop: PCOMPLS, i: np.ndarray) -> None:
        #prop.n = len(i)
        #prop.property_id = self.property_id[i]
        ##prop.material_id = self.material_id[i]

        #prop.coord_id = self.coord_id[i]
        #prop.sb = self.sb[i]
        #prop.tref = self.tref[i]
        #prop.ge = self.ge[i]

        #iply = self.iply
        #prop.global_ply_id = hslice_by_idim(i, iply, self.global_ply_id)
        #prop.material_id = hslice_by_idim(i, iply, self.material_id)
        #prop.thickness = hslice_by_idim(i, iply, self.thickness)
        #prop.theta = hslice_by_idim(i, iply, self.theta)

        #prop.nply = self.nply[i]

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PCOMPLS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        direct = integer_or_blank(card, 2, 'direct', default=1)
        cordm = integer_or_blank(card, 3, 'cordm', default=0)
        sb = double_or_blank(card, 4, 'sb')
        analysis = string_or_blank(card, 5, 'tref', default='ISH')
        nfields = len(card) - 1
        #nrows =
        ifield = 9
        global_ply_ids = []
        mids = []
        thicknesses = []
        thetas = []
        iply = 1
        c8: list[str] = []
        c20: list[str] = []
        while ifield < nfields:
            value = card[ifield]
            """
            | PCOMPLS | PID | DIRECT | CORDM |   SB   |  ANAL  |
            |         | C8  |  BEH8  | INT8  | BEH8H  | INT8H  |
            |         | C20 |  BEH20 | INT20 | BEH20H | INT20H |
            |         | ID1 |  MID1  |   T1  | THETA1 |        |
            |         | ID2 |  MID2  |   T2  | THETA2 |        |
            """
            if value == 'C8':
                assert c8 == [], c8
                #['C8', 'SLCOMP', 'L', None, None, None, None, None]

                #Element structural behavior. See Remarks 4. and 7.
                #(Character default: SLCOMP for BEH8 and BEH20)
                behavior_structural8 = string_or_blank(card, ifield+1, 'behavior_structural8', default='SLCOMP')
                # (Character default: L for INT8, Q for INT20)
                integration_structural8 = string_or_blank(card, ifield+2, 'integration_structural8', default='L')

                # BEHiH Element heat behavior. See Remarks 4. and 8.
                # (Character Default: SLCOMP for BEH8H and BEH20H)
                behavior_heat8 = string_or_blank(card, ifield+3, 'behavior_heat8', default='SLCOMP')

                # INTiH Integration scheme. See Remark 9.
                # (Character Default: L for INT8H, Q for INT20H)
                integration_heat8 = string_or_blank(card, ifield+4, 'integration_heat8', default='L')
                c8 = [
                    behavior_structural8, integration_structural8,
                    behavior_heat8, integration_heat8,
                ]
                ifield += 8
                continue
            elif value == 'C20':
                assert c20 == [], c20
                #['C8', 'SLCOMP', 'L', None, None, None, None, None]

                #Element structural behavior. See Remarks 4. and 7.
                #(Character default: SLCOMP for BEH8 and BEH20)
                behavior_structural20 = string_or_blank(card, ifield+1, 'behavior_structural20', default='SLCOMP')
                # (Character default: L for INT8, Q for INT20)
                integration_structural20 = string_or_blank(card, ifield+2, 'integration_structural20', default='Q')

                # BEHiH Element heat behavior. See Remarks 4. and 8.
                # (Character Default: SLCOMP for BEH8H and BEH20H)
                behavior_heat20 = string_or_blank(card, ifield+3, 'behavior_heat20', default='SLCOMP')

                # INTiH Integration scheme. See Remark 9.
                # (Character Default: L for INT8H, Q for INT20H)
                integration_heat20 = string_or_blank(card, ifield+4, 'integration_heat20', default='Q')
                c8 = [
                    behavior_structural20, integration_structural20,
                    behavior_heat20, integration_heat20,
                ]
                ifield += 8
                continue
            global_ply_id = integer(card, ifield, 'global_ply_id_%d' % iply)
            mid = integer(card, ifield + 1, 'mid_%d' % iply)
            t = double(card, ifield + 2, 'thickness_%d' % iply)
            theta = double(card, ifield + 3, 'theta_%d' % iply)
            global_ply_ids.append(global_ply_id)
            mids.append(mid)
            thicknesses.append(t)
            thetas.append(theta)
            iply += 1
            ifield += 8
        assert len(card) <= ifield, f'len(PCOMPLS card) = {len(card):d}\ncard={card}'
        #return PCOMPLS(pid, global_ply_ids, mids, thicknesses, thetas,
                       #direct=direct, cordm=cordm, sb=sb, analysis=analysis,
                       #c8=c8, c20=c20,
                       #comment=comment)
        self.cards.append((pid, global_ply_ids, mids, thicknesses, thetas,
                           direct, cordm, sb, analysis,
                           c8, c20,
                           comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        property_id = np.zeros(ncards, dtype='int32')
        #material_id = np.zeros(ncards, dtype='int32')
        coord_id = np.zeros(ncards, dtype='int32')
        #psdir = np.zeros(ncards, dtype='int32')
        sb = np.zeros(ncards, dtype='float64')
        #nb = np.zeros(ncards, dtype='float64')
        #tref = np.zeros(ncards, dtype='float64')
        #ge = np.zeros(ncards, dtype='float64')
        analysis = np.zeros(ncards, dtype='|U4')
        direct = np.zeros(ncards, dtype='int32')

        self.nply = np.zeros(ncards, dtype='int32')
        all_global_ply_ids = []
        all_mids = []
        all_thicknesses = []
        all_thetas = []
        #all_failure_theories = []
        #all_interlaminar_failure_theories = []
        #all_souts = []

        c8_list = []
        c20_list = []
        for icard, card in enumerate(self.cards):
            (pid, global_ply_idi, midi, thicknessi, thetai,
             directi, cordmi, sbi, analysisi,
             c8i, c20i, comment) = card
            property_id[icard] = pid
            coord_id[icard] = cordmi
            sb[icard] = sbi
            analysis[icard] = analysisi
            assert len(analysisi) <= 4, analysisi
            direct[icard] = directi
            #print('8:', c8i)
            #print('20:', c20i)
            c8_list.extend(c8i)
            c20_list.extend(c20i)
            #self.nb[icard] = nb
            #self.tref[icard] = tref
            #self.ge[icard] = ge
            self.nply[icard] = len(global_ply_idi)

            all_global_ply_ids.extend(global_ply_idi)
            all_mids.extend(midi)
            all_thicknesses.extend(thicknessi)
            all_thetas.extend(thetai)
            #all_failure_theories.extend(failure_theories)
            #all_interlaminar_failure_theories.extend(interlaminar_failure_theories)
            #all_souts.extend(souts)

        global_ply_id = np.array(all_global_ply_ids, dtype='int32')
        material_id = np.array(all_mids, dtype='int32')
        thickness = np.array(all_thicknesses, dtype='float64')
        theta = np.array(all_thetas, dtype='float64')
        c8 = np.array(c8_list, dtype='|U8')
        c20 = np.array(c20_list, dtype='|U8')
        #failure_theory = np.array(all_failure_theories, dtype='|U8')
        #interlaminar_failure_theory = np.array(all_interlaminar_failure_theories, dtype='|U8')

        self._save(property_id, material_id, global_ply_id,
                   theta, thickness, coord_id,
                   c8, c20, direct, analysis, sb)

        #for sout in souts:
            #assert len(sout) <= 3, sout
        #self.sout = np.array(all_souts, dtype='|U3')

        self.sort()
        self.cards = []

    def _save(self, property_id, material_id, global_ply_id,
              theta, thickness, coord_id,
              c8, c20,
              direct, analysis, sb):
        if len(self.property_id) != 0:
            raise NotImplementedError()
        self.property_id = property_id
        self.material_id = material_id
        self.global_ply_id = global_ply_id
        self.theta = theta
        self.thickness = thickness
        self.coord_id = coord_id

        self.direct = direct
        self.analysis = analysis
        self.sb = sb

        self.c8 = c8
        self.c20 = c20
        #self.model.log.warning(f'saved PCOMPLS; material_id={self.material_id}')
        assert len(self.material_id) > 0, self.material_id

    @property
    def iply(self) -> np.ndarray:
        return make_idim(self.n, self.nply)

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        property_ids = array_str(self.property_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        coord_ids = array_default_int(self.coord_id, default=0, size=size)
        sbs = array_float_nan(self.sb, size=size, is_double=False)

        for pid, mid, cordm, direct, analysis, sb, iply in zip(property_ids, material_ids, coord_ids,
                                                               self.direct, self.analysis, sbs, self.iply):
            iply0, iply1 = iply
            global_ply_ids = self.global_ply_id[iply0:iply1]
            mids = self.material_id[iply0:iply1]
            thicknesses = self.thickness[iply0:iply1]
            thetas = self.theta[iply0:iply1]
            c8 = self.c8[iply0:iply1]
            c20 = self.c20[iply0:iply1]
            #failure_theories = self.failure_theory[iply0:iply1]
            #interlaminar_failure_theories = self.interlaminar_failure_theory[iply0:iply1]
            #souts = self.sout[iply0:iply1]
            list_fields = ['PCOMPLS', pid, direct, cordm, sb, analysis, None, None, None]
            print('c8 =', c8)
            print('c20 =', c20)
            if len(c8):
                list_fields += ['C8'] + c8.tolist() + [None, None, None]
            if len(c20):
                list_fields += ['C20'] + c20.tolist() + [None, None, None]

            for glply, mid, t, theta in zip(global_ply_ids,
                                            mids, thicknesses, thetas):
                list_fields += [glply, mid, t, theta, None, None, None, None]
            bdf_file.write(print_card(list_fields))
        return

    @property
    def allowed_materials(self) -> list[Any]:
        """TODO: what about SOL 200 or undefined case control decks?"""
        model = self.model
        #MAT1, MAT9, MATORT, MATHE, MATUSR or MATDIGI
        #if model.is_thermal:
            #all_materials = [model.mat4]
        #else:
        #all_materials = [model.mat1, model.mat10]
        all_materials = [model.mat9, model.matort]
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def rho(self) -> np.ndarray:
        print(self.material_id, self.model.matort.rho)
        rho = get_density_from_material(self.material_id, self.allowed_materials, debug=True)
        return rho

class CHACBR(SolidElement):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.nodes = np.zeros((0, 20), dtype='int32')
        self.nnode_base = 8
        self.nnode = 20

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'), integer(card, 8, 'nid6'),
            integer(card, 9, 'nid7'), integer(card, 10, 'nid8'),
            integer_or_blank(card, 11, 'nid9', default=0),
            integer_or_blank(card, 12, 'nid10', default=0),
            integer_or_blank(card, 13, 'nid11', default=0),
            integer_or_blank(card, 14, 'nid12', default=0),
            0, # integer_or_blank(card, 15, 'nid13', default=0),
            0, # integer_or_blank(card, 16, 'nid14', default=0),
            0, # integer_or_blank(card, 17, 'nid15', default=0),
            0, # integer_or_blank(card, 18, 'nid16', default=0),
            integer_or_blank(card, 19, 'nid17', default=0),
            integer_or_blank(card, 20, 'nid18', default=0),
            integer_or_blank(card, 21, 'nid19', default=0),
            integer_or_blank(card, 22, 'nid20', default=0),
        ]
        assert len(card) <= 23, f'len(CHACAB card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 20), dtype=idtype)
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes)

    @property
    def base_nodes(self) -> np.ndarray:
        base_nodes = self.nodes[:, :8]
        return base_nodes

    @property
    def midside_nodes(self) -> np.ndarray:
        midside_nodes = self.nodes[:, 8:]
        if midside_nodes.shape[1] == 0:
            return midside_nodes
        assert midside_nodes.shape[1] == 12, midside_nodes.shape
        return midside_nodes

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        base_nodes = self.base_nodes
        midside_nodes = self.midside_nodes
        assert base_nodes.min() > 0, (base_nodes.min(), base_nodes.max())

        if midside_nodes.shape[1] == 0:
            for eid, pid, nodes in zip(self.element_id, self.property_id, base_nodes):
                msg = print_card(['CHACBR', eid, pid] + nodes.tolist())
                bdf_file.write(msg)
        else:
            for eid, pid, nodes in zip(self.element_id, self.property_id, self.nodes):
                msg = print_card(['CHACBR', eid, pid] + nodes.tolist())
        return

    def centroid(self) -> np.ndarray:
        centroid = chexa_centroid(self)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def volume(self) -> np.ndarray:
        xyz = self.model.grid.xyz_cid0()
        nid = self.model.grid.node_id
        nodes = self.base_nodes

        inode = np.searchsorted(nid, nodes)
        n1 = xyz[inode[:, 0], :]
        n2 = xyz[inode[:, 1], :]
        n3 = xyz[inode[:, 2], :]
        n4 = xyz[inode[:, 3], :]
        n5 = xyz[inode[:, 4], :]
        n6 = xyz[inode[:, 5], :]
        n7 = xyz[inode[:, 6], :]
        n8 = xyz[inode[:, 7], :]
        volume = volume_chexa(n1, n2, n3, n4, n5, n6, n7, n8)
        return volume

    def mass(self) -> np.ndarray:
        #material_id = get_material_from_property(self.property_id, self.allowed_properties)
        #rho = get_density_from_property(self.property_id, self.allowed_properties)
        #if rho.max() == 0. and rho.min() == 0.:
        #return np.zeros(len(rho), rho.dtype)
        #mass = rho * self.volume()
        mass = np.zeros(len(self.element_id), dtype='float64')
        return mass

    def quality(self):
        out = chexa_quality(self)
        return out


class CHACAB(SolidElement):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.nodes = np.zeros((0, 20), dtype='int32')
        self.nnode_base = 8
        self.nnode = 20

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'nid1'), integer(card, 4, 'nid2'),
            integer(card, 5, 'nid3'), integer(card, 6, 'nid4'),
            integer(card, 7, 'nid5'), integer(card, 8, 'nid6'),
            integer(card, 9, 'nid7'), integer(card, 10, 'nid8'),
            integer_or_blank(card, 11, 'nid9', default=0),
            integer_or_blank(card, 12, 'nid10', default=0),
            integer_or_blank(card, 13, 'nid11', default=0),
            integer_or_blank(card, 14, 'nid12', default=0),
            0, # integer_or_blank(card, 15, 'nid13', default=0),
            0, # integer_or_blank(card, 16, 'nid14', default=0),
            0, # integer_or_blank(card, 17, 'nid15', default=0),
            0, # integer_or_blank(card, 18, 'nid16', default=0),
            integer_or_blank(card, 19, 'nid17', default=0),
            integer_or_blank(card, 20, 'nid18', default=0),
            integer_or_blank(card, 21, 'nid19', default=0),
            integer_or_blank(card, 22, 'nid20', default=0),
        ]
        assert len(card) <= 23, f'len(CHACAB card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 20), dtype=idtype)
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes)

    @property
    def base_nodes(self) -> np.ndarray:
        base_nodes = self.nodes[:, :8]
        return base_nodes

    @property
    def midside_nodes(self) -> np.ndarray:
        midside_nodes = self.nodes[:, 8:]
        if midside_nodes.shape[1] == 0:
            return midside_nodes
        assert midside_nodes.shape[1] == 12, midside_nodes.shape
        return midside_nodes

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        base_nodes = self.base_nodes
        midside_nodes = self.midside_nodes
        assert base_nodes.min() > 0, (base_nodes.min(), base_nodes.max())

        if midside_nodes.shape[1] == 0:
            for eid, pid, nodes in zip(self.element_id, self.property_id, base_nodes):
                msg = print_card(['CHACAB', eid, pid] + nodes.tolist())
                bdf_file.write(msg)
        else:
            for eid, pid, nodes in zip(self.element_id, self.property_id, self.nodes):
                msg = print_card(['CHACAB', eid, pid] + nodes.tolist())
                bdf_file.write(msg)
        return

    def centroid(self) -> np.ndarray:
        centroid = chexa_centroid(self)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def volume(self) -> np.ndarray:
        xyz = self.model.grid.xyz_cid0()
        nid = self.model.grid.node_id
        nodes = self.base_nodes

        inode = np.searchsorted(nid, nodes)
        n1 = xyz[inode[:, 0], :]
        n2 = xyz[inode[:, 1], :]
        n3 = xyz[inode[:, 2], :]
        n4 = xyz[inode[:, 3], :]
        n5 = xyz[inode[:, 4], :]
        n6 = xyz[inode[:, 5], :]
        n7 = xyz[inode[:, 6], :]
        n8 = xyz[inode[:, 7], :]
        volume = volume_chexa(n1, n2, n3, n4, n5, n6, n7, n8)
        return volume

    def mass(self) -> np.ndarray:
        #material_id = get_material_from_property(self.property_id, self.allowed_properties)
        #rho = get_density_from_property(self.property_id, self.allowed_properties)
        #if rho.max() == 0. and rho.min() == 0.:
        #return np.zeros(len(rho), rho.dtype)
        #mass = rho * self.volume()
        mass = np.zeros(len(self.element_id), dtype='float64')
        return mass

    def quality(self):
        out = chexa_quality(self)
        return out
