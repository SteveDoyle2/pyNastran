from __future__ import annotations
from typing import TYPE_CHECKING # Optional, Any,
import numpy as np
#from pyNastran.bdf.cards.elements.bars import set_blank_if_default
#from pyNastran.bdf.cards.elements.solid import volume4
#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.bdf_interface.assign_type_force import force_double_or_blank
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, integer_or_blank, double_or_blank,
    integer_string_or_blank, string_or_blank)

from pyNastran.dev.bdf_vectorized3.utils import hstack_msg
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check, find_missing
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, make_idim, hslice_by_idim,
    parse_check, ) # searchsorted_filter,
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_default_str, array_default_int,
    array_float, array_float_nan, get_print_card_size)
from .utils import get_density_from_material, get_density_from_property, basic_mass_material_id

#from .solid_quality import chexa_quality, tetra_quality, penta_quality, pyram_quality, Quality
#from .solid_utils import chexa_centroid
#from .solid_volume import volume_ctetra, volume_cpenta, volume_cpyram, volume_chexa
from .solid import SolidElement, chexa_quality, chexa_centroid, volume_chexa

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


class CHACBR(SolidElement):
    def clear(self):
        self.n = 0
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
        return self.n - 1

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

    @property
    def allowed_properties(self) -> list[Any]:
        model = self.model
        allowed_properties = [model.pacbar]
        return [prop for prop in allowed_properties if prop.n > 0]

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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
    def clear(self):
        self.n = 0
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
        return self.n - 1

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

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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

class PACBAR(Property):
    """
    +--------+-----+-------+--------+--------+--------+---------+------+
    |    1   |  2  |   3   |   4    |    5   |    6   |    7    |   8  |
    +========+=====+=======+========+========+========+=========+======+
    | PACBAR | PID | MBACK | MSEPTM | FRESON | KRESON |         |      |
    +--------+-----+-------+--------+--------+--------+---------+------+
    | PACBAR |  12 |  1.0  |  0.01  |  400.0 |        |         |      |
    +--------+-----+-------+--------+--------+--------+---------+------+

    PID Property identification number. (Integer > 0)
    MBACK Mass per unit area of the backing material. (Real > 0.0)
    MSEPTM Mass per unit area of the septum material. (Real > 0.0)
    FRESON Resonant frequency of the sandwich construction in hertz. (Real > 0.0 or blank)
    KRESON Resonant stiffness of the sandwich construction. (Real > 0.0 or blank)
    """
    _show_attributes = ['property_id', 'mass_backing', 'mass_septum',
                        'freq_resonant', 'k_resonant']
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.mass_backing = np.array([], dtype='float64')
        self.mass_septum = np.array([], dtype='float64')

        self.freq_resonant = np.array([], dtype='float64')
        self.k_resonant = np.array([], dtype='float64')

    def __apply_slice__(self, prop: PACBAR, i: np.ndarray) -> None:  # ignore[override]
        prop.property_id = self.property_id[i]
        prop.mass_backing = self.mass_backing[i]
        prop.mass_septum = self.mass_septum[i]
        prop.freq_resonant = self.freq_resonant[i]
        prop.k_resonant = self.k_resonant[i]
        prop.n = len(i)

    def add(self, pid: int, mass_backing: float, mass_septum: float,
            freq_resonant: float, k_resonant: float, comment: str='') -> int:
        self.cards.append((pid, mass_backing, mass_septum, freq_resonant, k_resonant, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        mass_backing = double(card, 2, 'mass_backing')
        mass_septum = double(card, 3, 'mass_septum')
        freq_resonant = double_or_blank(card, 4, 'freq_resonant', default=np.nan)
        k_resonant = double_or_blank(card, 5, 'k_resonant', default=np.nan)
        assert len(card) <= 6, f'len(PACBAR card) = {len(card):d}\ncard={card}'

        self.cards.append((pid, mass_backing, mass_septum, freq_resonant, k_resonant, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)

        property_id = np.zeros(ncards, dtype='int32')
        mass_backing = np.zeros(ncards, dtype='float64')
        mass_septum = np.zeros(ncards, dtype='float64')
        freq_resonant = np.zeros(ncards, dtype='float64')
        k_resonant = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, mass_backingi, mass_septumi, freq_resonanti, k_resonanti, comment) = card

            property_id[icard] = pid
            mass_backing[icard] = mass_backingi
            mass_septum[icard] = mass_septumi

            freq_resonant[icard] = freq_resonanti
            k_resonant[icard] = k_resonanti
        self._save(property_id, mass_backing, mass_septum, freq_resonant, k_resonant)
        self.sort()
        self.cards = []

    def _save(self, property_id, mass_backing, mass_septum, freq_resonant, k_resonant):
        if len(self.property_id) != 0:
            raise NotImplementedError()
            #property_id = np.hstack([self.property_id, property_id])
            #mass_backing = np.hstack([self.mass_backing, mass_backing])
            #mass_septum = np.hstack([self.mass_septum, mass_septum])
            #freq_resonant = np.hstack([self.freq_resonant, freq_resonant])
            #k_resonant = np.hstack([self.k_resonant, k_resonant])

        nproperties = len(property_id)
        self.property_id = property_id
        self.mass_backing = mass_backing
        self.mass_septum = mass_septum

        self.freq_resonant = freq_resonant
        self.k_resonant = k_resonant
        self.n = nproperties

    #def set_used(self, used_dict: [str, list[np.ndarray]]) -> None:
        #used_dict['material_id'].append(self.material_id)

    #def geom_check(self, missing: dict[str, np.ndarray]):
        #model = self.model
        #cid = model.coord.coord_id
        #mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          #msg=f'no solid materials for {self.type}')
        #ucoord_id = np.unique(self.coord_id)
        #ucoord_id = np.setdiff1d(ucoord_id, [-1])
        #mids.sort()
        #geom_check(self,
                   #missing,
                   #coord=(cid, ucoord_id),
                   #material_id=(mids, self.material_id))

    @property
    def max_id(self) -> int:
        return self.property_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        property_id = array_str(self.property_id, size=size)
        mass_backings = array_float(self.mass_backing, size=size, is_double=False)
        mass_septums = array_float(self.mass_septum, size=size, is_double=False)
        freq_resonants = array_float(self.freq_resonant, size=size, is_double=False)
        k_resonants = array_float_nan(self.k_resonant, size=size, is_double=False)
        for pid, mass_backing, mass_septum, freq_resonant, k_resonant in zip(
            property_id, mass_backings, mass_septums, freq_resonants, k_resonants):
            fields = ['PACBAR', pid, mass_backing, mass_septum, freq_resonant, k_resonant]
            bdf_file.write(print_card(fields))
        return

    #@property
    #def allowed_materials(self) -> list[Any]:
        #"""TODO: what about SOL 200 or undefined case control decks?"""
        ##MAT1, MAT3, MAT4, MAT5, MAT9, MAT10,
        ##MAT10C, MATF10C, MAT11, or MATPOR
        #model = self.model
        #if model.is_thermal:
            #all_materials = [model.mat4, model.mat5]
        #else:
            #all_materials = [model.mat1, model.mat9, model.mat10, model.mat11]
        #materials = [mat for mat in all_materials if mat.n > 0]
        #assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        #return materials

    #def rho(self) -> np.ndarray:
        #rho = get_density_from_material(self.material_id, self.allowed_materials)
        #return rho
