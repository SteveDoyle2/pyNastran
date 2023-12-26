from __future__ import annotations
from abc import abstractmethod
from itertools import zip_longest
from typing import Any, TYPE_CHECKING

import numpy as np
#from pyNastran.bdf.field_writer_8 import print_field_8 # print_card_8,
#from pyNastran.bdf.field_writer_16 import print_card_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, # string, # double,
    integer_or_blank, double_or_blank,
    #integer_double_or_blank, string_or_blank,
    blank,
    #integer_types,
)
#from pyNastran.bdf.cards.elements.bars import set_blank_if_default

from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property,
    #hslice_by_idim, make_idim, searchsorted_filter,
    parse_element_check)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_float,
    array_default_int, array_default_float, array_float_nan,
    get_print_card_size)
from .utils import get_density_from_material
from .shell import (
    tri_centroid, tri_area, # tri_area_centroid_normal, tri_quality_xyz, tri_quality_nodes,
    quad_area, quad_centroid, # quad_area_centroid_normal, quad_quality_nodes,
    shell_thickness, shell_mass_per_area,
)
from .shell_properties import nonlinear_thickness
from .shell_quality import tri_quality_nodes, quad_quality_nodes

from pyNastran.dev.bdf_vectorized3.utils import hstack_msg

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    #from pyNastran.bdf.cards.materials import MAT1, MAT8
    #from pyNastran.dev.bdf_vectorized3.cards.grid import GRID


class PlateStressElement(Element):
    @Element.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')

    def set_used(self, used_dict: [dict[str, list[np.ndarray]]]) -> None:
        used_dict['property_id'].append(self.property_id)
        #used_dict['material_id'].append(self.material_id)
        used_dict['node_id'].append(self.nodes.ravel())

    @property
    def all_properties(self) -> list[Any]:
        model = self.model
        all_props = [model.pplane]
        return all_props
    @property
    def allowed_properties(self) -> list[Any]:
        all_props = self.all_properties
        props = [prop for prop in all_props if prop.n > 0]
        assert len(props) > 0, f'{self.type}: all_props={all_props}'
        return props

    @abstractmethod
    def area(self):
        ...

    def total_thickness(self) -> np.ndarray:
        #print(self.tflag)
        #print(self.T)
        tflag = None
        T = None
        thickness = shell_thickness(
            self.model,
            tflag, T,
            self.property_id, self.allowed_properties)
        #thickness = shell_thickness(self.property_id, self.allowed_properties)
        inan = np.isnan(thickness)
        if np.any(inan):
            log = self.model.log
            pids = np.unique(self.property_id[inan])
            log.warning(f'eids={self.element_id[inan]} with pids={pids} has nan thickness')
            #allowed_properties = self.allowed_properties
            #iprops = self.get_allowed_property_index(allowed_properties)
            #print('iprops =', iprops)
            ##print(self.type, '.properties', self.property_id)
            #uprops = np.unique(iprops)
            #for iprop in uprops:
                ##print('pids in pshell', self.model.pshell.property_id)
                #i = np.where(iprop == iprops)[0]
                #pids = self.property_id[i]
                #print('iprop =', iprop, type(iprop))
                #print('pids =', pids, type(pids))
                #propcard = allowed_properties[iprop]
                #prop = propcard.slice_card_by_property_id(pids)
                #log.warning(prop.write(size=8))
        return thickness

    def mass(self) -> np.ndarray:
        A = self.area()
        #t = self.total_thickness()
        #mass = rho * A * t
        tflag = None
        T = None
        mass_per_area = shell_mass_per_area(
            self.model, tflag, T,
            self.property_id, self.allowed_properties)
        mass = mass_per_area * A
        return mass

    def volume(self) -> np.ndarray:
        A = self.area()
        t = self.total_thickness()
        volume = A * t
        inan = np.isnan(volume)
        if np.any(inan):
            msg = (f'{self.type} has nan volume; volume={volume[inan]}\n'
                   f'element_id={self.element_id[inan]}'
                   f'property_id={self.property_id[inan]}\n'
                   f'area={A[inan]}\n'
                   f't={t[inan]}')
            self.model.log.error(msg)
            raise RuntimeError(msg)
        return volume

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max())


class CPLSTS3(PlateStressElement):
    """
    +---------+-------+-----+----+-------+----+-------+-------+-----+
    |    1    |   2   |  3  | 4  |    5  |  6 |   7   |   8   |  9  |
    +=========+=======+=====+====+=======+====+=======+=======+=====+
    | CPLSTS3 |  EID  | PID | N1 |   N2  | N3 |       | THETA |     |
    +---------+-------+-----+----+-------+----+-------+-------+-----+
    |         |       |     |    | TFLAG | T1 |   T2  |   T3  |     |
    +---------+-------+-----+----+-------+----+-------+-------+-----+

    per NX 2019.2
    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 3), dtype='int32')
        self.tflag = np.array([], dtype='int32')
        self.T = np.zeros((0, 3), dtype='float64')
        self.theta = np.array([], dtype='float64')

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CPLSTS3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
        ]
        if len(card) > 5:
            blank(card, 6, 'blank')
            theta = double_or_blank(card, 7, 'theta', default=0.0)
            blank(card, 8, 'blank')
            blank(card, 9, 'blank')
            blank(card, 10, 'blank')
            blank(card, 11, 'blank')
            tflag = integer_or_blank(card, 12, 'tflag', default=0)
            T1 = double_or_blank(card, 13, 'T1')
            T2 = double_or_blank(card, 14, 'T2')
            T3 = double_or_blank(card, 15, 'T3')
            assert len(card) <= 16, f'len(CPLSTS3 card) = {len(card):d}\ncard={card}'
        else:
            theta = 0.0
            tflag = 0
            T1 = 1.0
            T2 = 1.0
            T3 = 1.0
        self.cards.append((eid, pid, nids, theta,
                           tflag, [T1, T2, T3],
                           comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 3), dtype=idtype)
        theta = np.zeros(ncards, dtype='float64')

        tflag = np.zeros(ncards, dtype='int32')
        T = np.zeros((ncards, 3), dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, tflagi, ti, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard] = nids
            theta[icard] = thetai
            tflag[icard] = tflagi
            T[icard] = ti
        self._save(element_id, property_id, nodes, theta, tflag, T)
        self.sort()
        self.cards = []

        #element._save(element_id, property_id, nodes, theta, tflag, T)

    def _save(self, element_id: np.ndarray, property_id: np.ndarray,
              nodes: np.ndarray, theta: np.ndarray,
              tflag: np.ndarray, T: np.ndarray) -> None:
        if len(self.element_id) != 0:
            raise NotImplementedError()
        nelements = len(element_id)
        assert nelements > 0, element_id
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta
        self.tflag = tflag
        self.T = T
        assert nodes.shape == (nelements, 3), nodes.shape
        assert T.shape == (nelements, 3), T.shape
        self.n = nelements

    def __apply_slice__(self, element: CPLSTS3, i: np.ndarray) -> None:  # ignore[override]
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.tflag = self.tflag[i]
        element.theta = self.theta[i]
        element.T = self.T[i, :]
        element.n = len(self.element_id)

    def area(self) -> np.ndarray:
        area = tri_area(self.model.grid, self.nodes)
        return area

    def centroid(self) -> np.ndarray:
        centroid = tri_centroid(self.model.grid, self.nodes)
        return centroid

    def quality(self):
        return tri_quality_nodes(self.model.grid, self.nodes)

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes = array_str(self.nodes, size=size)
        thetas = array_default_float(self.theta, default=0.0, size=size, is_double=False)
        Ts = array_float_nan(self.T, size=size, is_double=False)

        for eid, pid, (n1, n2, n3), theta, tflag, (t1, t2, t3) in zip_longest(
            element_id, property_id, nodes, thetas, self.tflag, Ts):
            list_fields = ['CPLSTS3', eid, pid, n1, n2, n3, None, theta, None,
                            None, None, None, tflag, t1, t2, t3]
            bdf_file.write(print_card(list_fields))
        return


class CPLSTS4(PlateStressElement):
    """
    +---------+-------+-------+----+-------+----+-------+-------+------+
    |    1    |   2   |   3   |  4 |   5   |  6 |   7   |   8   |   9  |
    +=========+=======+=======+====+=======+====+=======+=======+======+
    | CPLSTS4 |  EID  |  PID  | N1 |   N2  | N3 |   N4  | THETA |      |
    +---------+-------+-------+----+-------+----+-------+-------+------+
    |         |       |       |    | TFLAG | T1 |   T2  |   T3  |  T4  |
    +---------+-------+-------+----+-------+----+-------+-------+------+

    ['CPLSTS4', '1', '5', '17', '18', '19', '20', '0.0']
    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 4), dtype='int32')
        self.tflag = np.array([], dtype='int32')
        self.T = np.zeros((0, 4), dtype='float64')
        self.theta = np.array([], dtype='float64')

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CPLSTS4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer(card, 6, 'n4'),
        ]
        if len(card) > 6:
            theta = double_or_blank(card, 7, 'theta', default=0.0)
            blank(card, 8, 'blank')
            blank(card, 9, 'blank')

            tflag = integer_or_blank(card, 12, 'tflag', default=0)
            T1 = double_or_blank(card, 13, 'T1')
            T2 = double_or_blank(card, 14, 'T2')
            T3 = double_or_blank(card, 15, 'T3')
            T4 = double_or_blank(card, 16, 'T4')
            assert len(card) <= 17, f'len(CPLSTS4 card) = {len(card):d}\ncard={card}'
        else:
            theta = 0.0
            tflag = 0
            T1 = 1.0
            T2 = 1.0
            T3 = 1.0
            T4 = 1.0
        self.cards.append((eid, pid, nids, theta,
                           tflag, [T1, T2, T3, T4],
                           comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 4), dtype=idtype)
        theta = np.zeros(ncards, dtype='float64')

        tflag = np.zeros(ncards, dtype='int32')
        T = np.zeros((ncards, 4), dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, tflagi, ti, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard] = nids
            theta[icard] = thetai
            tflag[icard] = tflagi
            T[icard] = ti
        self._save(element_id, property_id, nodes, theta, tflag, T)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray, property_id: np.ndarray,
              nodes: np.ndarray, theta: np.ndarray,
              tflag: np.ndarray, T: np.ndarray) -> None:
        if len(self.element_id) != 0:
            raise NotImplementedError()
        nelements = len(element_id)
        assert nelements > 0, element_id
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta
        self.tflag = tflag
        self.T = T
        self.n = nelements

    def __apply_slice__(self, element: CPLSTS4, i: np.ndarray) -> None:  # ignore[override]
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.tflag = self.tflag[i]
        element.theta = self.theta[i]
        element.T = self.T[i, :]
        element.n = len(self.element_id)

    def area(self) -> np.ndarray:
        area = quad_area(self.model.grid, self.nodes)
        return area

    def centroid(self) -> np.ndarray:
        centroid = quad_centroid(self.model.grid, self.nodes)
        return centroid

    def quality(self):
        return quad_quality_nodes(self.model.grid, self.nodes)

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)
        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes = array_str(self.nodes, size=size)

        thetas = array_default_float(self.theta, default=0.0, size=size, is_double=False)
        Ts = array_float_nan(self.T, size=size, is_double=False)
        for eid, pid, (n1, n2, n3, n4), theta, tflag, (t1, t2, t3, t4) in zip_longest(
            element_id, property_id, nodes, thetas, self.tflag, Ts):
            list_fields = ['CPLSTS4', eid, pid, n1, n2, n3, n4, theta,
                            None, None, None, None, tflag, t1, t2, t3, t4]
            bdf_file.write(print_card(list_fields))
        return


class PPLANE(Property):
    """NX specific card"""
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.coord_id = np.array([], dtype='int32')
        self.stress_strain_output_location = np.array([], dtype='|U4')
        self.thickness = np.array([], dtype='float64')

    #def add(self, pid: int, mid1: int=None, t: float=None,
            #mid2: int=None, twelveIt3: float=1.0,
            #mid3: int=None, tst: float=0.833333, nsm: float=0.0,
            #z1: float=None, z2: float=None, mid4: int=None,
            #comment: str='') -> PSHELL:
        #"""
        #Creates a PSHELL card

        #Parameters
        #----------
        #pid : int
            #property id
        #mid1 : int; default=None
            #defines membrane material
            #defines element density (unless blank)
        #mid2 : int; default=None
            #defines bending material
            #defines element density if mid1=None
        #mid3 : int; default=None
            #defines transverse shear material
        #mid4 : int; default=None
            #defines membrane-bending coupling material
        #twelveIt3 : float; default=1.0
            #Bending moment of inertia ratio, 12I/T^3. Ratio of the actual
            #bending moment inertia of the shell, I, to the bending
            #moment of inertia of a homogeneous shell, T^3/12. The default
            #value is for a homogeneous shell.
        #nsm : float; default=0.0
            #non-structural mass per unit area
        #z1 / z2 : float; default=None
            #fiber distance location 1/2 for stress/strain calculations
            #z1 default : -t/2 if thickness is defined
            #z2 default : t/2 if thickness is defined
        #comment : str; default=''
            #a comment for the card

        #"""
        #if z1 is None and t is not None:
            #z1 = -t / 2.
        #if z2 is None and t is not None:
            #z2 = t / 2.
        #t = np.nan if t is None else t

        #self.cards.append((pid, mid1, t,
                           #mid2, twelveIt3, mid3, tst, nsm, z1, z2, mid4,
                           #comment))
        #self.n += 1

    def add(self, pid: int, mid: int, t: float=0.0, nsm: float=0.0,
            formulation_option: int=0, comment: str='') -> int:
        self.cards.append((pid, mid, t, nsm, formulation_option, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PPLANE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        if self.debug:
            self.model.log.debug(f'adding card {card}')
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')  # MATHE, MATHP
        t = double_or_blank(card, 3, 't', default=0.)
        nsm = double_or_blank(card, 4, 'nsm', default=0.)
        formulation_option = integer_or_blank(card, 5, 'formulation_option', default=0)
        self.cards.append((pid, mid, t, nsm, formulation_option, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        formulation_option = np.zeros(ncards, dtype='int32')
        thickness = np.zeros(ncards, dtype='float64')
        nsm = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, mid, t, nsmi, formulation_optioni, comment) = card
            property_id[icard] = pid
            material_id[icard] = mid
            thickness[icard] = t
            nsm[icard] = t
            formulation_option[icard] = formulation_optioni

        thickness[(thickness == 0.)] = np.nan
        self._save(property_id, material_id, thickness, nsm, formulation_option)
        self.sort()
        self.cards = []

    def _save(self, property_id, material_id, thickness, nsm, formulation_option):
        self.property_id = property_id
        self.material_id = material_id
        self.thickness = thickness
        self.nsm = nsm
        self.formulation_option = formulation_option

    def set_used(self, used_dict: [str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)

    def __apply_slice__(self, prop: PPLANE, i: np.ndarray) -> None:  # ignore[override]
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]
        prop.thickness = self.thickness[i]
        prop.nsm = self.nsm[i]
        prop.formulation_option = self.formulation_option[i]

    def write_file_8(self, bdf_file: TextIOLike,
                   write_card_header: bool=False) -> None:
        self.write_file(bdf_file, size=8, is_double=False,
                        write_card_header=write_card_header)

    def write_file_16(self, bdf_file: TextIOLike,
                      is_double: bool=False, write_card_header: bool=False) -> None:
        self.write_file(bdf_file, size=16, is_double=is_double,
                        write_card_header=write_card_header)


    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.material_id.max())

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.property_id) == 0:
            return
        print_card, size = get_print_card_size(size, self.max_id)

        pids = array_str(self.property_id, size=size)
        mids = array_str(self.material_id, size=size)
        thicknesses = array_float(self.thickness, size=size, is_double=is_double)
        for pid, mid, thickness, nsm, formulation_option in zip_longest(
            pids, mids, thicknesses, self.nsm, self.formulation_option):
            list_fields = ['PPLANE', pid, mid, thickness, nsm, formulation_option]
            msg = print_card(list_fields)
            bdf_file.write(msg)
        return

    @property
    def all_materials(self) -> list[Any]:
        #shell_materials(self.model)
        return [self.model.mat1]

    @property
    def allowed_materials(self) -> list[Any]:
        all_materials = self.all_materials
        materials = [mat for mat in all_materials if mat.n > 0]
        assert len(materials) > 0, f'{self.type}: all_allowed_materials={all_materials}\nall_materials={self.model.material_cards}'
        return materials

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          msg=f'no pplane materials for {self.type}')
        mids.sort()
        material_ids = np.unique(self.material_id.ravel())
        geom_check(self,
                   missing,
                   material_id=(mids, material_ids))

    def mass_per_area(self) -> np.ndarray:
        #self.model.log.warning(f'faking mass/area for {self.type}')
        thickness = self.thickness
        #thickness = np.zeros(len(self.property_id))
        #nsm = self.nsm
        #mid = self.material_id

        rho = get_density_from_material(self.material_id, self.allowed_materials)
        mass_per_area = rho * thickness
        return mass_per_area

    @property
    def all_properties(self) -> list[Any]:
        model = self.model
        all_props = [model.pshln2]
        return all_props
    @property
    def allowed_properties(self) -> list[Any]:
        all_props = self.all_properties
        props = [prop for prop in all_props if prop.n > 0]
        assert len(props) > 0, f'{self.type}: all_props={all_props}'
        return props

    def total_thickness(self) -> np.ndarray:
        is_all_nan = np.all(np.isnan(self.thickness))
        if is_all_nan:
            allowed_properties = self.allowed_properties
            thickness = nonlinear_thickness(self.property_id, allowed_properties)
        else:
            thickness = self.thickness
        return thickness


class PlateStrainElement(Element):
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')

    def set_used(self, used_dict: [dict[str, list[np.ndarray]]]) -> None:
        used_dict['property_id'].append(self.property_id)
        #used_dict['material_id'].append(self.material_id)
        used_dict['node_id'].append(self.nodes.ravel())

    def __apply_slice__(self, element: PlateStrainElement, i: np.ndarray) -> None:  # ignore[override]
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.theta = self.theta[i]
        element.n = len(self.element_id)

    @property
    def all_properties(self) -> list[Any]:
        model = self.model
        all_props = [model.pplane]
        return all_props
    @property
    def allowed_properties(self) -> list[Any]:
        all_props = self.all_properties
        props = [prop for prop in all_props if prop.n > 0]
        assert len(props) > 0, f'{self.type}: all_props={all_props}'
        return props

    @abstractmethod
    def area(self):
        ...

    def total_thickness(self) -> np.ndarray:
        #print(self.tflag)
        #print(self.T)
        T = None
        tflag = None
        thickness = shell_thickness(
            self.model, tflag, T,
            self.property_id, self.allowed_properties)
        inan = np.isnan(thickness)
        if np.any(inan):
            log = self.model.log
            pids = np.unique(self.property_id[inan])
            log.warning(f'eids={self.element_id[inan]} with pids={pids} has nan thickness')
            #allowed_properties = self.allowed_properties
            #iprops = self.get_allowed_property_index(allowed_properties)
            #print('iprops =', iprops)
            ##print(self.type, '.properties', self.property_id)
            #uprops = np.unique(iprops)
            #for iprop in uprops:
                ##print('pids in pshell', self.model.pshell.property_id)
                #i = np.where(iprop == iprops)[0]
                #pids = self.property_id[i]
                #print('iprop =', iprop, type(iprop))
                #print('pids =', pids, type(pids))
                #propcard = allowed_properties[iprop]
                #prop = propcard.slice_card_by_property_id(pids)
                #log.warning(prop.write(size=8))
        return thickness

    def mass(self) -> np.ndarray:
        A = self.area()
        #t = self.total_thickness()
        #mass = rho * A * t
        mass_per_area = shell_mass_per_area(
            self.model, None, None,
            self.property_id, self.allowed_properties)
        mass = mass_per_area * A
        return mass

    def volume(self) -> np.ndarray:
        A = self.area()
        t = self.total_thickness()
        volume = A * t
        inan = np.isnan(volume)
        if np.any(inan):
            msg = (f'{self.type} has nan volume; volume={volume[inan]}\n'
                   f'element_id={self.element_id[inan]}'
                   f'property_id={self.property_id[inan]}\n'
                   f'area={A[inan]}\n'
                   f't={t[inan]}')
            self.model.log.error(msg)
            raise RuntimeError(msg)
        return volume

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max())


class CPLSTN3(PlateStrainElement):
    """
    +---------+-------+-------+----+----+----+-------+-------+-----+
    |    1    |   2   |   3   |  4 |  5 |  6 |   7   |   8   |  9  |
    +=========+=======+=======+=====+===+====+=======+=======+=====+
    | CPLSTN3 |  EID  |  PID  | N1 | N2 | N3 | THETA |       |     |
    +---------+-------+-------+----+----+----+-------+-------+-----+

    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 3), dtype='int32')
        self.theta = np.array([], dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
            comment: str='') -> int:
        """Creates a CPLSTN3 card"""
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CPLSTS3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', default=eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
        ]
        theta = double_or_blank(card, 6, 'theta', default=0.0)
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 3), dtype=idtype)
        theta = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard] = nids
            theta[icard] = thetai
        self._save(element_id, property_id, nodes, theta)
        self.sort()
        self.cards = []

        #element._save(element_id, property_id, nodes, theta, tflag, T)

    def _save(self, element_id: np.ndarray, property_id: np.ndarray,
              nodes: np.ndarray, theta: np.ndarray) -> None:
        if len(self.element_id) != 0:
            raise NotImplementedError()
        nelements = len(element_id)
        assert nelements > 0, element_id
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta
        assert nodes.shape == (nelements, 3), nodes.shape
        self.n = nelements

    def area(self) -> np.ndarray:
        area = tri_area(self.model.grid, self.nodes)
        return area

    def centroid(self) -> np.ndarray:
        centroid = tri_centroid(self.model.grid, self.nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def quality(self):
        return tri_quality_nodes(self.model.grid, self.nodes)

    def flip_normal(self, i=None) -> None:
        """
        Flips normal of element.

        ::

               1           1            2
              * *   -->   * *    -->   * *
             *   *       *   *        *   *
            2-----3     3-----2      3-----1
            nominal     fast flip   perserve material orientation

        """
        if i is None:
            i = slice(len(self.element_id))
        #self.nodes[i, :] = self.nodes[i, [0, 2, 1]] # fast flip
        self.nodes[i, :] = self.nodes[i, [1, 0, 2]]  # preserve material orientation

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes = array_str(self.nodes, size=size)
        for eid, pid, (n1, n2, n3), theta in zip_longest(
            element_id, property_id, nodes, self.theta):
            list_fields = ['CPLSTN3', eid, pid, n1, n2, n3, theta]
            bdf_file.write(print_card(list_fields))
        return

class CPLSTN4(PlateStrainElement):
    """
    +---------+-------+-------+----+-------+----+-------+-------+------+
    |    1    |   2   |   3   |  4 |   5   |  6 |   7   |   8   |   9  |
    +=========+=======+=======+====+=======+====+=======+=======+======+
    | CPLSTN4 |  EID  |  PID  | N1 |   N2  | N3 |   N4  | THETA |      |
    +---------+-------+-------+----+-------+----+-------+-------+------+

    ['CPLSTN4', '1', '5', '17', '18', '19', '20', '0.0']
    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 4), dtype='int32')
        self.theta = np.array([], dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
            comment: str='') -> int:
        """Creates a CPLSTN4 card"""
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CPLSTS4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', default=eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer(card, 6, 'n4'),
        ]
        theta = double_or_blank(card, 7, 'theta', default=0.0)
        assert len(card) <= 8, f'len(CPLSTS4 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 4), dtype=idtype)
        theta = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard] = nids
            theta[icard] = thetai
        self._save(element_id, property_id, nodes, theta)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray, property_id: np.ndarray,
              nodes: np.ndarray, theta: np.ndarray) -> None:
        if len(self.element_id) != 0:
            raise NotImplementedError()
        nelements = len(element_id)
        assert nelements > 0, element_id
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta
        self.n = nelements

    def area(self) -> np.ndarray:
        area = quad_area(self.model.grid, self.nodes)
        return area

    def centroid(self) -> np.ndarray:
        centroid = quad_centroid(self.model.grid, self.nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def quality(self):
        return quad_quality_nodes(self.model.grid, self.nodes)

    def flip_normal(self, i=None):
        """
        1-2-3-4
        2-1-4-3
        """
        if i is None:
            i = slice(len(self.element_id))
        self.nodes[i, :] = self.nodes[i, [1, 0, 3, 2]]

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)


        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes = array_str(self.nodes, size=size)
        for eid, pid, (n1, n2, n3, n4), theta in zip_longest(
            element_id, property_id, nodes, self.theta):
            list_fields = ['CPLSTN4', eid, pid, n1, n2, n3, n4, theta]
            bdf_file.write(print_card(list_fields))
        return


class CPLSTN6(PlateStrainElement):
    """
    +---------+-------+-------+----+----+----+----+----+----+
    |    1    |   2   |   3   |  4 |  5 |  6 |  7 | 8  | 9  |
    +=========+=======+=======+=====+===+====+====+====+====+
    | CPLSTN6 |  EID  |  PID  | N1 | N2 | N3 | N4 | N5 | N6 |
    +---------+-------+-------+----+----+----+----+----+----+
    |         | THETA |       |    |    |    |    |    |    |
    +---------+-------+-------+----+----+----+----+----+----+

    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 6), dtype='int32')
        self.theta = np.array([], dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
            comment: str='') -> int:
        """Creates a CPLSTN6 card"""
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CPLSTN6 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', default=eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4', default=0),
            integer_or_blank(card, 7, 'n5', default=0),
            integer_or_blank(card, 8, 'n6', default=0),
        ]
        theta = double_or_blank(card, 9, 'theta', default=0.0)
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 6), dtype=idtype)
        theta = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard] = nids
            theta[icard] = thetai
        self._save(element_id, property_id, nodes, theta)
        self.sort()
        self.cards = []

        #element._save(element_id, property_id, nodes, theta, tflag, T)

    def _save(self, element_id: np.ndarray, property_id: np.ndarray,
              nodes: np.ndarray, theta: np.ndarray) -> None:
        if len(self.element_id) != 0:
            raise NotImplementedError()
        nelements = len(element_id)
        assert nelements > 0, element_id
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta
        assert nodes.shape == (nelements, 6), nodes.shape
        self.n = nelements

    def area(self) -> np.ndarray:
        area = tri_area(self.model.grid, self.base_nodes)
        return area

    def centroid(self) -> np.ndarray:
        centroid = tri_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def quality(self):
        return tri_quality_nodes(self.model.grid, self.base_nodes)

    @property
    def base_nodes(self):
        return self.nodes[:, :3]
    @property
    def midside_nodes(self):
        return self.nodes[:, 3:]

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes = array_str(self.nodes, size=size)
        for eid, pid, (n1, n2, n3, n4, n5, n6), theta in zip_longest(
            element_id, property_id, nodes, self.theta):
            list_fields = ['CPLSTN6', eid, pid, n1, n2, n3, n4, n5, n6,
                           theta]
            bdf_file.write(print_card(list_fields))
        return


class CPLSTN8(PlateStrainElement):
    """
    +---------+-------+-------+-------+----+----+----+----+----+
    |    1    |   2   |   3   |   4   |  5 |  6 |  7 | 8  | 9  |
    +=========+=======+=======+=======+===+====+====+====+====+
    | CPLSTN6 |  EID  |  PID  |  N1   | N2 | N3 | N4 | N5 | N6 |
    +---------+-------+-------+-------+----+----+----+----+----+
    |         |  N7   |  N8   | THETA |    |    |    |    |    |
    +---------+-------+-------+-------+----+----+----+----+----+

    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 8), dtype='int32')
        self.theta = np.array([], dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
            comment: str='') -> int:
        """Creates a CPLSTN8 card"""
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CPLSTN8 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', default=eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer(card, 6, 'n4'),
            integer_or_blank(card, 7, 'n5', default=0),
            integer_or_blank(card, 8, 'n6', default=0),
            integer_or_blank(card, 9, 'n7', default=0),
            integer_or_blank(card, 10, 'n8', default=0),
        ]
        theta = double_or_blank(card, 11, 'theta', default=0.0)
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 8), dtype=idtype)
        theta = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard] = nids
            theta[icard] = thetai
        self._save(element_id, property_id, nodes, theta)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray, property_id: np.ndarray,
              nodes: np.ndarray, theta: np.ndarray) -> None:
        if len(self.element_id) != 0:
            raise NotImplementedError()
        nelements = len(element_id)
        assert nelements > 0, element_id
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta
        assert nodes.shape == (nelements, 8), nodes.shape
        self.n = nelements

    def area(self) -> np.ndarray:
        area = quad_area(self.model.grid, self.base_nodes)
        return area

    def centroid(self) -> np.ndarray:
        centroid = quad_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def quality(self):
        return quad_quality_nodes(self.model.grid, self.base_nodes)

    @property
    def base_nodes(self):
        return self.nodes[:, :4]
    @property
    def midside_nodes(self):
        return self.nodes[:, 4:]

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes = array_str(self.nodes, size=size)
        for eid, pid, (n1, n2, n3, n4, n5, n6, n7, n8), theta in zip_longest(
            element_id, property_id, nodes, self.theta):
            list_fields = ['CPLSTN8', eid, pid, n1, n2, n3, n4, n5, n6, n7, n8,
                           theta]
            bdf_file.write(print_card(list_fields))
        return


class CPLSTS6(PlateStrainElement):
    """
    +---------+-------+-----+-------+-------+----+----+----+----+
    |    1    |   2   |  3  |   4   |    5  |  6 |  7 |  8 |  9 |
    +=========+=======+=====+=======+=======+====+====+====+====+
    | CPLSTS8 |  EID  | PID |  G1   |   G2  | G3 | G4 | G5 | G6 |
    +---------+-------+-----+-------+-------+----+----+----+----+
    |         |       |     | THETA | TFLAG | T1 | T2 | T3 |    |
    +---------+-------+-----+-------+-------+----+----+----+----+
    |         |       |     |       |       | T4 | T5 | T6 |    |
    +---------+-------+-----+-------+-------+----+----+----+----+

    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 6), dtype='int32')
        self.theta = np.array([], dtype='float64')
        self.tflag = np.array([], dtype='int32')
        self.thickness = np.zeros((0, 6), dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
            tflag: int=0, thickness=None, comment: str='') -> int:
        """Creates a CPLSTS6 card"""
        if thickness is None:
            thickness = np.full(6, np.nan, dtype='float64')
        self.cards.append((eid, pid, nids, theta, tflag, thickness, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CPLSTS6 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5', default=0),
                integer_or_blank(card, 8, 'n6', default=0),]
        # 9, 10
        if len(card) > 11:
            theta = double_or_blank(card, 11, 'theta', default=0.0)
            tflag = integer_or_blank(card, 12, 'tflag', default=0)
            thickness = [
                double_or_blank(card, 13, 't1', default=np.nan),
                double_or_blank(card, 14, 't2', default=np.nan),
                double_or_blank(card, 15, 't3', default=np.nan),
                #double_or_blank(card, 16, 't4', default=np.nan),
                #17, 18, 19, 20
                double_or_blank(card, 21, 't4', default=np.nan),
                double_or_blank(card, 22, 't5', default=np.nan),
                double_or_blank(card, 23, 't6', default=np.nan),
            ]
            assert len(card) <= 23, f'len(CPLSTS6 card) = {len(card):d}\ncard={card}'
        else:
            theta = 0.0
            tflag = 0
            thickness = np.full(6, np.nan, dtype='float64')
        self.cards.append((eid, pid, nids, theta, tflag, thickness, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 6), dtype=idtype)
        theta = np.zeros(ncards, dtype='float64')
        tflag = np.zeros(ncards, dtype='int32')
        thickness = np.zeros((ncards, 6), dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, tflagi, thicknessi, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard] = nids
            theta[icard] = thetai
            tflag[icard] = tflagi
            thickness[icard, :] = thicknessi
        self._save(element_id, property_id, nodes, theta, tflag, thickness)
        self.sort()
        self.cards = []

        #element._save(element_id, property_id, nodes, theta, tflag, T)

    def _save(self, element_id: np.ndarray, property_id: np.ndarray,
              nodes: np.ndarray, theta: np.ndarray, tflag: np.ndarray, thickness: np.ndarray) -> None:
        if len(self.element_id):
            raise NotImplementedError()
        nelements = len(element_id)
        assert nelements > 0, element_id
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta
        self.tflag = tflag
        self.thickness = thickness
        assert nodes.shape == (nelements, 6), nodes.shape
        self.n = nelements

    def area(self) -> np.ndarray:
        area = tri_area(self.model.grid, self.base_nodes)
        return area

    def centroid(self) -> np.ndarray:
        centroid = tri_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def quality(self):
        return tri_quality_nodes(self.model.grid, self.base_nodes)

    @property
    def base_nodes(self):
        return self.nodes[:, :3]
    @property
    def midside_nodes(self):
        return self.nodes[:, 3:]

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes = array_str(self.nodes, size=size)
        tflags = array_default_int(self.tflag, size=size, default=0)
        thicknesses = array_float_nan(self.thickness, size=size, is_double=False)
        for eid, pid, (n1, n2, n3, n4, n5, n6), theta, tflag, thickness in zip_longest(
            element_id, property_id, nodes, self.theta, tflags, thicknesses):
            list_fields = ['CPLSTS6', eid, pid, n1, n2, n3, n4, n5, n6, '', '',
                           theta, tflag]
            if not np.all(thickness == ''):
                t1, t2, t3, t4, t5, t6 = thickness
                list_fields.extend([t1, t2, t3, '', '', '', '', '', t4, t5, t6])
            bdf_file.write(print_card(list_fields))
        return


class CPLSTS8(PlateStrainElement):
    """
    +---------+-------+-----+-------+-------+----+----+----+----+
    |    1    |   2   |  3  |   4   |    5  |  6 |  7 |  8 |  9 |
    +=========+=======+=====+=======+=======+====+====+====+====+
    | CPLSTS8 |  EID  | PID |  G1   |   G2  | G3 | G4 | G5 | G6 |
    +---------+-------+-----+-------+-------+----+----+----+----+
    |         |   G7  | G8  | THETA | TFLAG | T1 | T2 | T3 | T4 |
    +---------+-------+-----+-------+-------+----+----+----+----+
    |         |       |     |       |       | T5 | T6 | T7 | T8 |
    +---------+-------+-----+-------+-------+----+----+----+----+

    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 8), dtype='int32')
        self.theta = np.array([], dtype='float64')
        self.tflag = np.array([], dtype='int32')
        self.thickness = np.zeros((0, 8), dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int], theta: float=0.0,
            tflag: int=0,
            thickness=None, comment: str='') -> int:
        """Creates a CPLSTS8 card"""
        if thickness is None:
            thickness = np.full(8, np.nan, dtype='float64')
        self.cards.append((eid, pid, nids, theta, tflag, thickness, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CPLSTS8 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4'),
                integer_or_blank(card, 7, 'n5', default=0),
                integer_or_blank(card, 8, 'n6', default=0),
                integer_or_blank(card, 9, 'n7', default=0),
                integer_or_blank(card, 10, 'n8', default=0),]
        if len(card) > 11:
            theta = double_or_blank(card, 11, 'theta', default=0.0)
            tflag = integer_or_blank(card, 12, 'tflag', default=0)
            thickness = [
                double_or_blank(card, 13, 't1', default=np.nan),
                double_or_blank(card, 14, 't2', default=np.nan),
                double_or_blank(card, 15, 't3', default=np.nan),
                double_or_blank(card, 16, 't4', default=np.nan),
                #17, 18, 19, 20
                double_or_blank(card, 21, 't5', default=np.nan),
                double_or_blank(card, 22, 't6', default=np.nan),
                double_or_blank(card, 23, 't7', default=np.nan),
                double_or_blank(card, 24, 't8', default=np.nan),
            ]
            assert len(card) <= 25, f'len(CPLSTS8 card) = {len(card):d}\ncard={card}'
        else:
            theta = 0.0
            tflag = 0
            thickness = np.full(8, np.nan, dtype='float64')
        self.cards.append((eid, pid, nids, theta, tflag, thickness, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        property_id = np.zeros(ncards, dtype=idtype)
        nodes = np.zeros((ncards, 8), dtype=idtype)
        theta = np.zeros(ncards, dtype='float64')
        tflag = np.zeros(ncards, dtype='int32')
        thickness = np.zeros((ncards, 8), dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, tflagi, thicknessi, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard] = nids
            theta[icard] = thetai
            tflag[icard] = tflagi
            thickness[icard] = thicknessi
        self._save(element_id, property_id, nodes, theta, tflag, thickness)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray, property_id: np.ndarray,
              nodes: np.ndarray, theta: np.ndarray,
              tflag: np.ndarray, thickness: np.ndarray) -> None:
        if len(self.element_id) != 0:
            raise NotImplementedError()
        nelements = len(element_id)
        assert nelements > 0, element_id
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta
        self.tflag = tflag
        self.thickness = thickness
        assert nodes.shape == (nelements, 8), nodes.shape
        self.n = nelements

    def area(self) -> np.ndarray:
        area = quad_area(self.model.grid, self.base_nodes)
        return area

    def centroid(self) -> np.ndarray:
        centroid = quad_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def quality(self):
        return quad_quality_nodes(self.model.grid, self.base_nodes)

    @property
    def base_nodes(self):
        return self.nodes[:, :4]
    @property
    def midside_nodes(self):
        return self.nodes[:, 4:]

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes = array_str(self.nodes, size=size)
        tflags = array_default_int(self.tflag, size=size, default=0)
        thicknesses = array_float_nan(self.thickness, size=size)
        for eid, pid, (n1, n2, n3, n4, n5, n6, n7, n8), theta, tflag, thickness, in zip_longest(
            element_id, property_id, nodes, self.theta, tflags, thicknesses):
            list_fields = ['CPLSTS8', eid, pid, n1, n2, n3, n4, n5, n6, n7, n8,
                           theta, tflag]
            if not np.all(thickness == ''):
                t1, t2, t3, t4, t5, t6, t7, t8 = thickness
                list_fields.extend([t1, t2, t3, t4, '', '', '', '', t5, t6, t7, t8])
            bdf_file.write(print_card(list_fields))
        return

