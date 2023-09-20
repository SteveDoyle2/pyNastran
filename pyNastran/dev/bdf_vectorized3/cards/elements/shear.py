from __future__ import annotations
#from abc import abstractmethod
from itertools import zip_longest
from typing import Any, TYPE_CHECKING

import numpy as np
from pyNastran.bdf.field_writer_8 import print_card_8, print_field_8
from pyNastran.bdf.field_writer_16 import print_card_16, print_field_16
from pyNastran.dev.bdf_vectorized3.cards.base_card import Element, Property, searchsorted_filter
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double,
    integer_or_blank, double_or_blank)
from pyNastran.bdf.cards.elements.bars import set_blank_if_default
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int
from .utils import get_density_from_material
from .shell import quad_area, quad_centroid
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


class CSHEAR(Element):
    """
    +--------+-------+-------+----+----+----+----+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |
    +========+=======+=======+=====+===+====+====+
    | CSHEAR |  EID  |  PID  | N1 | N2 | N3 | N4 |
    +--------+-------+-------+----+----+----+----+

    """
    def add(self, eid: int, pid: int, nids: list[int], comment: str='') -> int:
        """
        Creates a CSHEAR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHEAR)
        nids : list[int, int, int, int]
            node ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 4), dtype='int32')

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4'),]
        self.cards.append((eid, pid, nids, comment))
        self.n += 1
        return self.n

    def parse_cards(self):
        assert self.n >= 0, self.n
        if self.n == 0 or len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 4), dtype='int32')
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
        self._save(element_id, property_id, nodes)
        self.cards = []

    def _save(self, element_id, property_id, nodes):
        if len(self.element_id):
            raise NotImplementedError()
        nelements = len(element_id)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.n = nelements


    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no shear properties for {self.type}')
        pids.sort()
        geom_check(self,
                   missing,
                   node=(nid, self.nodes),
                   property_id=(pids, self.property_id))

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.element_id) == 0:
            return

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodess = array_str(self.nodes, size=size).tolist()
        for eid, pid, nodes in zip_longest(element_id, property_id, nodess):
            data = [eid, pid] + nodes
            if size == 8:
                msg = 'CSHEAR  %8s%8s%8s%8s%8s%8s\n' % tuple(data)
            else:
                msg = print_card_16(data)
            bdf_file.write(msg)
        return

    @property
    def allowed_properties(self) -> list[Any]:
        return [prop for prop in [self.model.pshear]
                if prop.n > 0]

    def area(self) -> np.ndarray:
        area = quad_area(self.model.grid, self.nodes)
        return area
    def centroid(self) -> np.ndarray:
        """centroid ignores density"""
        centroid = quad_centroid(self.model.grid, self.nodes)
        return centroid
    def center_of_mass(self) -> np.ndarray:
        """center_of_mass considers density"""
        return self.centroid()

    def volume(self) -> np.ndarray:
        A = self.area()
        t = shear_mass_per_area(self.property_id, self.allowed_properties)
        return A * t

    def mass_per_area(self) -> np.ndarray:
        mass_per_area = shear_mass_per_area(
            self.property_id, self.allowed_properties)
        return mass_per_area

    def mass(self) -> np.ndarray:
        mass_per_area = self.mass_per_area()
        area = self.area()
        mass = mass_per_area * area
        return mass


class PSHEAR(Property):
    """
    Defines the properties of a shear panel (CSHEAR entry).

    +--------+-----+-----+---+-----+----+----+
    |   1    |  2  |  3  | 4 |  5  |  6 |  7 |
    +========+=====+=====+===+=====+====+====+
    | PSHEAR | PID | MID | T | NSM | F1 | F2 |
    +--------+-----+-----+---+-----+----+----+
    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.material_id = np.array([], dtype='int32')
        self.t = np.array([], dtype='float64')
        self.nsm = np.array([], dtype='float64')
        self.f1 = np.array([], dtype='float64')
        self.f2 = np.array([], dtype='float64')

    def add(self, pid: int, mid: int, t: float, nsm: float=0.,
            f1: float=0., f2: float=0., comment: str='') -> int:
        """
        Creates a PSHEAR card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        t : float
            shear panel thickness
        nsm : float; default=0.
            nonstructural mass per unit length
        f1 : float; default=0.0
            Effectiveness factor for extensional stiffness along edges 1-2 and 3-4
        f2 : float; default=0.0
            Effectiveness factor for extensional stiffness along edges 2-3 and 1-4
        comment : str; default=''
            a comment for the card

        """
        #self.property_id = np.hstack([self.property_id, pid])
        #self.material_id = np.hstack([self.material_id, mid])
        #self.t = np.hstack([self.t, t])
        #self.nsm = np.hstack([self.nsm, nsm])
        #self.f1 = np.hstack([self.f1, f1])
        #self.f2 = np.hstack([self.f2, f2])
        #self.n += 1
        self.cards.append((pid, mid, t, nsm, f1, f2, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        t = double(card, 3, 't')
        nsm = double_or_blank(card, 4, 'nsm', default=0.0)
        f1 = double_or_blank(card, 5, 'f1', default=0.0)
        f2 = double_or_blank(card, 6, 'f2', default=0.0)
        assert len(card) <= 7, f'len(PSHEAR card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, mid, t, nsm, f1, f2, comment))
        self.n += 1
        return self.n

    def parse_cards(self):
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        t = np.zeros(ncards, dtype='float64')
        nsm = np.zeros(ncards, dtype='float64')
        f1 = np.zeros(ncards, dtype='float64')
        f2 = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            pid, mid, ti, nsmi, f1i, f2i, comment = card

            property_id[icard] = pid
            material_id[icard] = mid
            t[icard] = ti
            nsm[icard] = nsmi
            f1[icard] = f1i
            f2[icard] = f2i
        self._save(property_id, material_id, t, nsm, f1, f2)
        self.cards = []

    def _save(self, property_id, material_id, t, nsm, f1, f2):
        if len(self.property_id):
            raise NotImplementedError()
        self.property_id = property_id
        self.material_id = material_id
        self.t = t
        self.nsm = nsm
        self.f1 = f1
        self.f2 = f2

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          msg=f'no materials for {self.type}; {self.all_materials}')
        mids.sort()
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.property_id) == 0:
            return
        for pid, mid, t, nsm, f1, f2 in zip_longest(self.property_id, self.material_id, self.t, self.nsm, self.f1, self.f2):
            nsm = set_blank_if_default(nsm, 0.0)
            list_fields = ['PSHEAR', pid, mid, t, nsm,
                           f1, f2]
            msg = print_card_8(list_fields)
            bdf_file.write(msg)
        return

    @property
    def all_materials(self) -> list[Any]:
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[Any]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.materials}'
        return materials

    def total_thickness(self):
        return self.t

    def mass_per_area(self) -> np.ndarray:
        thickness = self.t
        nsm = self.nsm
        mid = self.material_id

        rho = get_density_from_material(mid, self.allowed_materials)
        mass_per_area = nsm + rho * thickness
        return mass_per_area

def shear_mass_per_area(property_id: np.ndarray,
                        allowed_properties: list[Any]) -> np.ndarray:
    mass_per_area = np.full(len(property_id), np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties
    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue
        mass_per_areai = prop.mass_per_area()
        mass_per_area[ilookup] = mass_per_areai[iall]
    return mass_per_area
