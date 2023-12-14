from __future__ import annotations
from abc import abstractmethod
from itertools import zip_longest
from typing import Any, TYPE_CHECKING

import numpy as np
#from pyNastran.bdf.field_writer_8 import print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, # string, # double,
    integer_or_blank, double_or_blank,
    integer_double_or_blank,
    #string_or_blank, blank, integer_types,
)
#from pyNastran.bdf.cards.elements.bars import set_blank_if_default

from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import Element, searchsorted_filter, parse_element_check # Property, hslice_by_idim, make_idim
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_float,
    array_default_int, # array_default_float,
    get_print_card_size)
from .utils import get_density_from_property, get_density_from_material, expanded_mass_material_id
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg
from .shell import (
    tri_centroid, tri_area,
    quad_centroid, quad_area,
    shell_thickness, shell_nonstructural_mass, shell_mass_per_area,
    _check_shell_mass,
    NUMPY_INTS, NUMPY_FLOATS)

if TYPE_CHECKING:
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    #from pyNastran.bdf.cards.materials import MAT1, MAT8
    #from pyNastran2.bdf.cards.grid import GRID


class AxisymmetricShellElement(Element):
    """
    CQUADX (MSC): PLPLANE or PAXSYMH or PLCOMP
    CQUADX4 (NX): PSOLID, PLSOLID, PMIC
    CQUADX8 (NX): PSOLID, PLSOLID
    """
    @Element.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.array([], dtype='int32')

    @property
    def all_properties(self) -> list[Any]:
        model = self.model
        if self.type in {'CQUADX4', 'CQUADX8', 'CTRAX3', 'CTRAX6'}:
            ## TODO: CQUADX4: PMIC not supported
            all_props = [model.psolid, model.plsolid]
        elif self.type == 'CTRIAX':
            #PLPLANE, PAXSYMH
            ## TODO: CTRIAX: PAXSYMH not supported
            all_props = [model.plplane]
        else:
            ## TODO: CQUADX: PLPLANE or PAXSYMH or PLCOMP
            assert self.type in 'CQUADX', self.type
            all_props = [model.plplane] # , model.plcomp
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

    #def volume(self) -> np.ndarray:
        #A = self.area()
        #t = self.total_thickness()
        #volume = A * t
        #inan = np.isnan(volume)
        #if np.any(inan):
            #msg = (f'{self.type} has nan volume; volume={volume[inan]}\n'
                   #f'element_id={self.element_id[inan]}'
                   #f'property_id={self.property_id[inan]}\n'
                   #f'area={A[inan]}\n'
                   #f't={t[inan]}')
            #self.model.log.error(msg)
            #raise RuntimeError(msg)
        #return volume

    def get_allowed_property_index(self, allowed_properties: list[Any]) -> np.ndarray:
        indexi = np.full(len(self.property_id), -1, dtype='int32')
        assert len(allowed_properties) > 0, allowed_properties
        for i, prop in enumerate(allowed_properties):
            ilookup, iall = searchsorted_filter(prop.property_id, self.property_id)
            if len(iall) == 0:
                continue
            indexi[ilookup] = i
        return indexi

    #def mass_per_area(self) -> np.ndarray:
        #mass_per_area = shell_mass_per_area(
            #self.property_id, self.allowed_properties)
        #return mass_per_area

    #def total_thickness(self) -> np.ndarray:
        #"""TODO: doesn't consider differential thickness"""
        #total_thickness = shell_total_thickness(
            #self.property_id, self.allowed_properties)
        #return total_thickness

    #def nsm_per_area(self) -> np.ndarray:
        #total_thickness = shell_nonstructural_mass(
            #self.property_id, self.allowed_properties)
        #return total_thickness

    def mass(self) -> np.ndarray:
        #material_id = get_material_from_property(self.property_id, self.allowed_properties)
        rho = get_density_from_property(self.property_id, self.allowed_properties)
        if rho.max() == 0. and rho.min() == 0.:
            return np.zeros(len(rho), rho.dtype)
        mass = rho * self.volume()

        inan = np.isnan(mass)
        if np.any(inan):
            #pids = np.unique(self.property_id[inan])
            msg = f'{self.type} has nan mass'
            msg += f'element_id={self.element_id[inan]}'
            msg += f'property_id={self.property_id[inan]}\n'
            #msg += f'area={area[inan]}\n'
            #t = self.total_thickness()
            #msg += f't={t[inan]}\n'
            self.model.log.warning(msg)
            #self.model.log.warning(f'eids={self.element_id[inan]} with pids={pids} has nan mass')
        return mass
        #tscales = self.get_thickness_scale()
        #try:
            #mpa = self.pid_ref.MassPerArea(tflag=self.tflag, tscales=tscales)
        #except TypeError:
            #print(self.pid_ref)
            #raise

        #if mpa == 0.0:
            #return 0.0

        #area = self.area()
        #try:
            #return mpa * A
        #except TypeError:
            #msg = 'mass/area=%s area=%s prop_type=%s' % (mpa, A, self.pid_ref.type)
            #raise TypeError(msg)

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no shell properties for {self.type}')
        #for prop in self.allowed_properties:
            #print(prop.write(size=8))
        assert len(pids) > 0, self.allowed_properties
        pids.sort()

        base_nodes = self.base_nodes
        midside_nodes = self.midside_nodes
        assert base_nodes is not None
        #print(self.base_nodes)
        geom_check(self,
                   missing,
                   node=(nid, base_nodes), filter_node0=True,
                   property_id=(pids, self.property_id))
        if midside_nodes is not None:
            geom_check(self,
                       missing,
                       node=(nid, midside_nodes), filter_node0=True)


class AxiShellElement(Element):
    """
    CTRIAX (MSC):  PLPLANE, PAXSYMH
    CTRIAX  (NX):  PLPLANE, PAXSYMH
    CTRIAX6 (NX):  uses mid

    CQUADX (MSC): PLPLANE or PAXSYMH or PLCOMP
    CQUADX4 (NX): PSOLID, PLSOLID, PMIC
    CQUADX8 (NX): PSOLID, PLSOLID
    """
    @Element.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.nodes = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        #self.tflag = np.array([], dtype='int32')
        #self.T = np.array([], dtype='int32')

    def check_types(self):
        assert self.element_id.dtype.name in NUMPY_INTS, self.element_id.dtype.name
        assert self.property_id.dtype.name in NUMPY_INTS, self.property_id.dtype.name
        #assert self.tflag.dtype.name in {'int8', 'int32', 'int64'}, self.tflag.dtype.name

    @property
    def all_properties(self) -> list[Any]:
        model = self.model
        if self.type in {'CQUADX4', 'CQUADX8'}:
            ## TODO: CQUADX4: PMIC not supported
            all_props = [model.psolid, model.plsolid]
        else:
            ## TODO: CQUADX: PLPLANE or PAXSYMH or PLCOMP
            assert self.type == 'CQUADX', self.type
            all_props = [model.plplane] # , model.plcomp
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

    def volume(self) -> np.ndarray:
        A = self.area()
        t = self.total_thickness()
        volume = A * t
        inan = np.isnan(volume)
        if np.any(inan):
            msg = (f'{self.type} has nan volume; volume={volume[inan]}\n'
                   f'element_id={self.element_id[inan]}'
                   f'property_id={self.property_id[inan]}\n')
            if not np.all(np.isfinite(A[inan])):
                msg += f'area={A[inan]}\n'
            if not np.all(np.isfinite(t[inan])):
                msg += f't={t[inan]}\n'
            if np.any(np.isnan(t[inan])):
                msg += (
                    f'tflag={self.tflag[inan]}\n'
                    f'T={self.T[inan, :]}')
            self.model.log.error(msg)
            raise RuntimeError(msg)
        return volume

    def get_allowed_property_index(self, allowed_properties: list[Any]) -> np.ndarray:
        indexi = np.full(len(self.property_id), -1, dtype='int32')
        assert len(allowed_properties) > 0, allowed_properties
        for i, prop in enumerate(allowed_properties):
            ilookup, iall = searchsorted_filter(prop.property_id, self.property_id)
            if len(iall) == 0:
                continue
            indexi[ilookup] = i
        return indexi

    def mass_material_id(self) -> np.ndarray:
        element_id, property_id, material_id = expanded_mass_material_id(
            self.element_id, self.property_id, self.allowed_properties)
        assert material_id.min() > 0, material_id
        return material_id

    def detailed_mass(self) -> np.ndarray:
        element_id, property_id, material_id = expanded_mass_material_id(
            self.element_id, self.property_id, self.allowed_properties)
        assert material_id.min() > 0, material_id
        return material_id

    def total_thickness(self) -> np.ndarray:
        #print(self.tflag)
        #print(self.T)
        thickness = shell_thickness(self.model,
                                    self.tflag, self.T,
                                    self.property_id, self.allowed_properties)
        inan = np.isnan(thickness)
        if np.any(inan):
            log = self.model.log
            pids = np.unique(self.property_id[inan])
            log.warning(f'eids={self.element_id[inan]} with pids={pids} has nan thickness')
            allowed_properties = self.allowed_properties
            iprops = self.get_allowed_property_index(allowed_properties)
            #print('iprops =', iprops)
            #print(self.type, '.properties', self.property_id)
            uprops = np.unique(iprops)
            for iprop in uprops:
                #print('pids in pshell', self.model.pshell.property_id)
                i = np.where(iprop == iprops)[0]
                pids = self.property_id[i]
                print('iprop =', iprop, type(iprop))
                print('pids =', pids, type(pids))
                propcard = allowed_properties[iprop]
                prop = propcard.slice_card_by_property_id(pids)
                log.warning(prop.write(size=8))

        thickness = shell_thickness(self.model,
                                    self.tflag, self.T,
                                    self.property_id, self.allowed_properties)
        inan = np.isnan(thickness)
        if inan.sum():
            self.model.log.error(thickness[inan])
        assert thickness.sum() > 0., thickness
        return thickness

    def mass_per_area(self) -> np.ndarray:
        nelement = len(self.element_id)
        assert nelement > 0, nelement
        mass_per_area = shell_mass_per_area(
            self.model, self.tflag, self.T,
            self.property_id, self.allowed_properties)
        assert len(mass_per_area) == nelement, mass_per_area
        return mass_per_area

    #def total_thickness(self) -> np.ndarray:
        #"""TODO: doesn't consider differential thickness"""
        #total_thickness = shell_total_thickness(
            #self.property_id, self.allowed_properties)
        #return total_thickness

    def nsm_per_area(self) -> np.ndarray:
        total_thickness = shell_nonstructural_mass(
            self.property_id, self.allowed_properties)
        return total_thickness

    def mass(self) -> np.ndarray:
        """TODO: doesn't consider differential thickness"""
        mass_per_area = self.mass_per_area()

        area = self.area()

        #print('mass_per_area =', mass_per_area)
        #print('area =', area)
        mass = mass_per_area * area
        #print('*mass =', mass)
        _check_shell_mass(self, mass, area)
        return mass
        #tscales = self.get_thickness_scale()
        #try:
            #mpa = self.pid_ref.MassPerArea(tflag=self.tflag, tscales=tscales)
        #except TypeError:
            #print(self.pid_ref)
            #raise

        #if mpa == 0.0:
            #return 0.0

        #area = self.area()
        #try:
            #return mpa * A
        #except TypeError:
            #msg = 'mass/area=%s area=%s prop_type=%s' % (mpa, A, self.pid_ref.type)
            #raise TypeError(msg)

    def set_from_op2(self, element_id, property_id, nodes, zoffset=None,
                     tflag=None, T=None, theta=None, mcid=None):
        #(eid, pid, n1, n2, n3, n4, n5, n6, theta, zoffs, t1, t2, t3) = out
        nelements = len(element_id)
        assert element_id.min() > 0, element_id
        assert property_id.min() > 0, property_id
        assert nodes.min() >= 0, nodes

        if mcid is not None:
            assert mcid.min() >= -1, nodes

        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes

        if zoffset is None:
            zoffset = np.full(nelements, np.nan, dtype='float64')
        assert zoffset is not None
        assert zoffset.dtype.name in NUMPY_FLOATS, zoffset.dtype.name
        self.zoffset = zoffset

        if theta is None:
            theta = np.full(nelements, np.nan, dtype='float64')
        self.theta = theta

        if mcid is None:
            mcid = np.full(nelements, -1, dtype=theta.dtype)
        self.mcid = mcid

        if tflag is None:
            tflag = np.zeros(nelements, dtype=element_id.dtype)
        else:
            utflag = np.unique(tflag)
            assert tflag.min() in {0, 1}, utflag
            assert tflag.max() in {0, 1}, utflag
        self.tflag = tflag

        nbase_nodes = self.base_nodes.shape[1]
        if T is None:
            T = np.zeros((nelements, nbase_nodes), dtype=theta.dtype)
        assert T.shape == (nelements, nbase_nodes), T.shape
        self.T = T

        self.n = nelements
        self.check_types()

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        pids = hstack_msg([prop.property_id for prop in self.allowed_properties],
                          msg=f'no shell properties for {self.type}')
        #for prop in self.allowed_properties:
            #print(prop.write(size=8))
        assert len(pids) > 0, self.allowed_properties
        pids.sort()

        base_nodes = self.base_nodes
        midside_nodes = self.midside_nodes
        assert base_nodes is not None
        #print(self.base_nodes)
        geom_check(self,
                   missing,
                   node=(nid, base_nodes), filter_node0=True,
                   property_id=(pids, self.property_id))
        if midside_nodes is not None:
            geom_check(self,
                       missing,
                       node=(nid, midside_nodes), filter_node0=True)


class CTRIAX(AxisymmetricShellElement):
    """
    +--------+------------+-------+----+----+----+----+----+-----+
    |   1    |     2      |   3   |  4 |  5 |  6 | 7  |  8 |  9  |
    +========+============+=======+====+====+====+====+====+=====+
    | CTRIAX |    EID     |  PID  | N1 | N2 | N3 | N4 | N5 | N6  |
    +--------+------------+-------+----+----+----+----+----+-----+
    |        | THETA/MCID |       |    |    |    |    |    |     |
    +--------+------------+-------+----+----+----+----+----+-----+

    Theta/Mcid is MSC only!
    """
    def add(self, eid: int, pid: int, nids: list[int],
            theta_mcid: int|float=0., comment: str='') -> int:
        """Creates a CTRIAX card"""
        self.cards.append((eid, pid, nids, theta_mcid, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4', default=0),
            integer_or_blank(card, 7, 'n5', default=0),
            integer_or_blank(card, 8, 'n6', default=0),
        ]
        theta_mcid = integer_double_or_blank(card, 9, 'theta_mcsid', default=0.0)
        assert len(card) <= 10, f'len(CTRIAX card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, theta_mcid, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 6), dtype='int32')
        mcid = np.full(ncards, -1, dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, theta_mcid, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :len(nids)] = nids
            if isinstance(theta_mcid, float):
                theta[icard] = theta_mcid
            else:
                mcid[icard] = theta_mcid
        self._save(element_id, property_id, nodes, theta, mcid)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray,
              property_id: np.ndarray,
              nodes: np.ndarray,
              theta: np.ndarray,
              mcid: np.ndarray) -> None:
        if len(self.element_id) != 0:
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            nodes = np.vstack([self.nodes, nodes])
            theta = np.hstack([self.theta, theta])
            mcid = np.hstack([self.mcid, mcid])
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta
        self.mcid = mcid

    def __apply_slice__(self, element: CTRIAX, i: np.ndarray) -> None:
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.theta = self.theta[i]
        element.mcid = self.mcid[i]
        element.n = len(self.element_id)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max(), self.mcid.max())

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size).tolist()
        mcids = array_default_int(self.mcid, default=-1, size=size)
        thetas = array_float(self.theta, size=size, is_double=False)
        imcid = (self.mcid == -1)
        thetas[imcid] = mcids[imcid]

        for eid, pid, nodes, theta_mcid in zip_longest(element_ids, property_ids, nodes_, thetas):
            list_fields = ['CTRIAX', eid, pid] + nodes + [theta_mcid]
            bdf_file.write(print_card(list_fields))
        return

    def area(self) -> np.ndarray:
        return tri_area(self.model.grid, self.base_nodes)

    def centroid(self) -> np.ndarray:
        centroid = tri_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    @property
    def base_nodes(self):
        return self.nodes[:, :3]
    @property
    def midside_nodes(self):
        return self.nodes[:, 3:]



class CTRIAX6(Element):
    """
    +---------+-------+-------+----+----+----+----+----+-----+
    |    1    |   2   |   3   |  4 |  5 |  6 |  7 |  8 |  9  |
    +=========+=======+=======+=====+===+====+====+====+=====+
    | CTRIAX6 |  EID  |  MID  | N1 | N2 | N3 | G4 | G5 | G6  |
    +---------+-------+-------+----+----+----+----+----+-----+
    |         | THETA |       |    |    |    |    |    |     |
    +---------+-------+-------+----+----+----+----+----+-----+

    NX/MSC : Nodes are defined in a non-standard way::

           5
          / \
         6   4
       /       \
      1----2----3
    """
    @Element.clear_check
    def clear(self) -> None:
        self.nodes = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')

    def add(self, eid: int, mid: int, nids: list[int], theta: float=0.,
            comment: str='') -> int:
        """Creates a CTRIAX6 card"""
        self.cards.append((eid, mid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CTRIAX6 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        mid = integer(card, 2, 'mid')

        nids = [
            integer(card, 3, 'n1'),
            integer_or_blank(card, 4, 'n2', default=0),
            integer(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4', default=0),
            integer(card, 7, 'n5'),
            integer_or_blank(card, 8, 'n6', default=0),
        ]
        theta = double_or_blank(card, 9, 'theta', default=0.0)
        assert len(card) <= 10, f'len(CTRIAX6 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, mid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 6), dtype='int32')
        #mcid = np.full(ncards, -1, dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, mid, nids, thetai, comment) = card
            element_id[icard] = eid
            material_id[icard] = mid
            nodes[icard, :] = nids
            theta[icard] = thetai
        self._save(element_id, material_id, nodes, theta)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray,
              material_id: np.ndarray,
              nodes: np.ndarray,
              theta: np.ndarray) -> None:
        if len(self.element_id) != 0:
            element_id = np.hstack([self.element_id, element_id])
            material_id = np.hstack([self.material_id, material_id])
            nodes = np.vstack([self.nodes, nodes])
            theta = np.hstack([self.theta, theta])
        self.element_id = element_id
        self.material_id = material_id
        self.nodes = nodes
        self.theta = theta

    def __apply_slice__(self, element: CTRIAX6, i: np.ndarray) -> None:
        assert len(self.material_id) > 0
        element.element_id = self.element_id[i]
        element.material_id = self.material_id[i]
        element.nodes = self.nodes[i, :]
        element.theta = self.theta[i]
        element.n = len(self.element_id)

    @property
    def tflag(self) -> np.ndarray:
        nelement = len(self.element_id)
        return np.zeros(nelement, dtype='int32')

    @property
    def T(self) -> np.ndarray:
        nelement = len(self.element_id)
        return np.full((nelement, 3), np.nan, dtype='float64')

    @property
    def all_materials(self) -> list[Any]:
        model = self.model
        all_mats = [model.mat1]
        return all_mats

    @property
    def allowed_materials(self) -> list[Any]:
        all_mats = self.all_materials
        mats = [mat for mat in all_mats if mat.n > 0]
        assert len(mats) > 0, all_mats
        return mats

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.nodes.max(), self.material_id.max())

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size).tolist()
        thetas = array_float(self.theta, size=size, is_double=is_double)
        for eid, mid, nodes, theta in zip_longest(element_ids, material_ids, nodes_, thetas):
            list_fields = ['CTRIAX6', eid, mid] + nodes + [theta]
            bdf_file.write(print_card(list_fields))
        return

    def area(self) -> np.ndarray:
        return tri_area(self.model.grid, self.base_nodes)

    def centroid(self) -> np.ndarray:
        centroid = tri_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    @property
    def base_nodes(self) -> np.ndarray:
        return self.nodes[:, [0, 2, 4]]
    @property
    def midside_nodes(self) -> np.ndarray:
        return self.nodes[:, [1, 3, 5]]

    def volume(self) -> np.ndarray:
        self.model.log.warning(f'faking {self.type} volume')
        nelement = len(self.element_id)
        volume = np.zeros(nelement, dtype='float64')
        return volume

    def mass(self) -> np.ndarray:
        volume = self.volume()
        unused_area = self.area()
        rho = get_density_from_material(
            self.material_id, self.allowed_materials, debug=False)
        #print('mass_per_area =', mass_per_area)
        #print('area =', area)
        mass = volume * rho
        #print('*mass =', mass)
        #_check_shell_mass(self, mass, area)
        return mass


class CQUADX(AxiShellElement):
    """
    Defines an axisymmetric quadrilateral element with up to nine grid
    points for use in fully nonlinear (i.e., large strain and large
    rotations) analysis or a linear harmonic or rotordynamic analysis.
    The element has between four and eight grid points

    +--------+-------+-------+----+------------+----+----+-----+-----+
    |   1    |   2   |   3   |  4 |      5     |  6 | 7  |  8  |  9  |
    +========+=======+=======+====+============+====+====+=====+=====+
    | CQUADX |  EID  |  PID  | N1 |     N2     | N3 | N4 |  G5 | G6  |
    +--------+-------+-------+----+------------+----+----+-----+-----+
    |        |  G7   |  G8   | G9 | THETA/MCID |    |    |     |     |
    +--------+-------+-------+----+------------+----+----+-----+-----+

    Theta/Mcid is MSC only!
    """
    def add(self, eid: int, pid: int, nids: list[int],
            theta_mcid: int|float=0., comment: str='') -> int:
        self.cards.append((eid, pid, nids, theta_mcid, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CQUADX card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer(card, 6, 'n4'),
            integer_or_blank(card, 7, 'n5', default=0),
            integer_or_blank(card, 8, 'n6', default=0),
            integer_or_blank(card, 9, 'n7', default=0),
            integer_or_blank(card, 10, 'n8', default=0),
            integer_or_blank(card, 11, 'n9', default=0),
        ]
        theta_mcid = integer_double_or_blank(card, 12, 'theta/mcid', default=0.)
        assert len(card) <= 13, f'len(CQUADX card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, theta_mcid, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 9), dtype='int32')
        mcid = np.full(ncards, -1, dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, theta_mcid, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            if isinstance(theta_mcid, float):
                theta[icard] = theta_mcid
            else:
                mcid[icard] = theta_mcid
        self._save(element_id, property_id, nodes, theta, mcid)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray,
              property_id: np.ndarray,
              nodes: np.ndarray,
              theta: np.ndarray,
              mcid: np.ndarray) -> None:
        if len(self.element_id) != 0:
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            nodes = np.vstack([self.nodes, nodes])
            mcid = np.hstack([self.mcid, mcid])
            theta = np.hstack([self.theta, theta])
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.mcid = mcid
        self.theta = theta

    def set_from_op2(self, element_id, property_id, nodes):
        nelements = len(element_id)
        assert element_id.min() > 0, element_id
        assert property_id.min() > 0, property_id
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.mcid = np.full(nelements, -1, dtype=element_id.dtype)
        self.theta = np.full(nelements, np.nan, dtype='float32')
        self.n = nelements

    def __apply_slice__(self, element: CQUADX, i: np.ndarray) -> None:
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.mcid = self.mcid[i]
        element.theta = self.theta[i]
        element.n = len(self.element_id)

    #@property
    #def tflag(self):
        #return None
    #@property
    #def T(self):
        #return None
    @property
    def all_properties(self) -> list[Any]:
        model = self.model
        #return [PLPLANE or PAXSYMH or PLCOMP]
        return [model.plplane]

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max(), self.mcid.max())

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.element_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size).tolist()
        mcids = array_default_int(self.mcid, default=-1, size=size)

        thetas = array_float(self.theta, size=size, is_double=False)
        imcid = (self.mcid == -1)
        thetas[imcid] = mcids[imcid]
        for eid, pid, nodes, theta_mcid in zip_longest(element_ids, property_ids, nodes_, thetas):
            list_fields = ['CQUADX', eid, pid] + nodes + [theta_mcid]
            bdf_file.write(print_card(list_fields))
        return

    def area(self) -> np.ndarray:
        return quad_area(self.model.grid, self.base_nodes)

    def centroid(self) -> np.ndarray:
        centroid = quad_centroid(self.model.grid, self.base_nodes)
        return centroid

    @property
    def base_nodes(self) -> np.ndarray:
        return self.nodes[:, :4]
    @property
    def midside_nodes(self) -> np.ndarray:
        return self.nodes[:, 4:]


class CQUADX4(AxisymmetricShellElement):
    """
    Defines an isoparametric and axisymmetric quadrilateral cross-section
    ring element for use in linear and fully nonlinear (i.e., large strain
    and large rotations) hyperelastic analysis.

    +---------+-------+-------+----+----+----+----+-------+
    |    1    |   2   |   3   |  4 |  5 |  6 | 7  |   8   |
    +=========+=======+=======+====+====+====+====+=======+
    | CQUADX4 |  EID  |  PID  | N1 | N2 | N3 | N4 | THETA |
    +---------+-------+-------+----+----+----+----+-------+

    CQUADX4 is an NX card only!
    """
    def add(self, eid: int, pid: int, nids: list[int],
            theta: float=0., comment: str='') -> int:
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CQUADX4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer(card, 6, 'n4'),
        ]
        theta = double_or_blank(card, 7, 'theta', default=0.)
        assert len(card) <= 8, f'len(CQUADX4 card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 4), dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            theta[icard] = thetai
        self._save(element_id, property_id, nodes, theta)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray,
              property_id: np.ndarray,
              nodes: np.ndarray,
              theta: np.ndarray) -> None:
        if len(self.element_id) != 0:
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            nodes = np.vstack([self.nodes, nodes])
            theta = np.hstack([self.theta, theta])
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta

    def __apply_slice__(self, element: CQUADX4, i: np.ndarray) -> None:
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.theta = self.theta[i]
        element.n = len(self.element_id)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max())

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size).tolist()
        thetas = array_float(self.theta, size=size, is_double=is_double)
        for eid, pid, nodes, theta in zip_longest(element_ids, property_ids, nodes_, thetas):
            list_fields = ['CQUADX4', eid, pid] + nodes + [theta]
            bdf_file.write(print_card(list_fields))
        return

    def area(self) -> np.ndarray:
        return quad_area(self.model.grid, self.base_nodes)

    def volume(self) -> np.ndarray:
        model = self.model
        model.log.warning(f'faking {self.type} volume')
        nelements, nnodes = self.nodes.shape
        assert nnodes == 4, nnodes
        nid = model.grid.node_id
        xyz = model.grid.xyz_cid0()
        inode = np.searchsorted(nid, self.nodes)
        assert np.array_equal(nid[inode], self.nodes)
        in1 = inode[:, 0]
        in2 = inode[:, 1]
        in3 = inode[:, 2]
        in4 = inode[:, 3]
        xyz1 = xyz[in1, :]
        xyz2 = xyz[in2, :]
        xyz3 = xyz[in3, :]
        xyz4 = xyz[in4, :]
        nelement = len(self.element_id)
        #area = 0.5 * norm(cross(n3-n1, n4-n2))
        volume = np.zeros(nelement, dtype='float64')
        return volume

    def centroid(self) -> np.ndarray:
        centroid = quad_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    @property
    def base_nodes(self) -> np.ndarray:
        return self.nodes
    @property
    def midside_nodes(self):
        return None


class CQUADX8(AxisymmetricShellElement):
    """
    Defines an isoparametric and axisymmetric quadrilateral cross-section
    ring element with midside nodes for use in linear and fully nonlinear
    (i.e., large strain and large rotations) hyperelastic analysis.

    +---------+-------+-------+-------+----+----+----+-----+-----+
    |    1    |   2   |   3   |   4   |  5 |  6 | 7  |  8  |  9  |
    +=========+=======+=======+=======+====+====+====+=====+=====+
    | CQUADX8 |  EID  |  PID  |  N1   | N2 | N3 | N4 |  G5 | G6  |
    +---------+-------+-------+-------+----+----+----+-----+-----+
    |         |  G7   |  G8   | THETA |    |    |    |     |     |
    +---------+-------+-------+-------+----+----+----+-----+-----+

    CQUADX8 is an NX card only!
    """
    def add(self, eid: int, pid: int, nids: list[int],
            theta: float=0., comment: str='') -> int:
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CQUADX8 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
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
        theta = double_or_blank(card, 12, 'theta', default=0.)
        assert len(card) <= 13, f'len(CQUADX card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 8), dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            theta[icard] = thetai
        self._save(element_id, property_id, nodes, theta)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray,
              property_id: np.ndarray,
              nodes: np.ndarray,
              theta: np.ndarray) -> None:
        if len(self.element_id) != 0:
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            nodes = np.vstack([self.nodes, nodes])
            theta = np.hstack([self.theta, theta])
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta

    def __apply_slice__(self, element: CQUADX8, i: np.ndarray) -> None:
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.theta = self.theta[i]
        element.n = len(self.element_id)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max())

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size).tolist()
        thetas = array_float(self.theta, size=size, is_double=is_double)
        for eid, pid, nodes, theta in zip_longest(element_ids, property_ids, nodes_, thetas):
            data = [eid, pid] + nodes + [theta]
            msg = print_card(data)
            bdf_file.write(msg)
        return

    def area(self) -> np.ndarray:
        return quad_area(self.model.grid, self.base_nodes)

    def centroid(self) -> np.ndarray:
        centroid = quad_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    @property
    def base_nodes(self):
        return self.nodes[:, :4]
    @property
    def midside_nodes(self):
        return self.nodes[:, 4:]


class CTRAX3(AxisymmetricShellElement):
    """
    +--------+------------+-------+----+----+----+-------+
    |   1    |     2      |   3   |  4 |  5 |  6 |   7   |
    +========+============+=======+====+====+====+=======+
    | CTRAX3 |    EID     |  PID  | N1 | N2 | N3 | THETA |
    +--------+------------+-------+----+----+----+-------+

    CTRAX3 is NX only!
    """
    #type = 'CTRAX3'

    def add(self, eid: int, pid: int, nids: list[int],
            theta: float=0., comment: str='') -> int:
        """Creates a CTRAX3 card"""
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')

        nids = [
            integer_or_blank(card, 3, 'n1'),
            integer_or_blank(card, 4, 'n2'),
            integer_or_blank(card, 5, 'n3'),
        ]
        theta = integer_double_or_blank(card, 6, 'theta', default=0.0)
        assert len(card) <= 7, f'len(CTRAX3 card) = {len(card):d}\ncard={card}'
        #return CTRAX3(eid, pid, nids, theta=theta, comment=comment)
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 3), dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            theta[icard] = thetai
        self._save(element_id, property_id, nodes, theta)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray,
              property_id: np.ndarray,
              nodes: np.ndarray,
              theta: np.ndarray) -> None:
        if len(self.element_id) != 0:
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            nodes = np.vstack([self.nodes, nodes])
            theta = np.hstack([self.theta, theta])
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta

    def __apply_slice__(self, element: CTRAX3, i: np.ndarray) -> None:
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.theta = self.theta[i]
        element.n = len(self.element_id)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max())

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size).tolist()
        thetas = array_float(self.theta, size=size, is_double=is_double)
        for eid, pid, nodes, theta in zip_longest(element_ids, property_ids, nodes_, thetas):
            list_fields = ['CTRAX3', eid, pid] + nodes + [theta]
            bdf_file.write(print_card(list_fields))
        return

    def area(self) -> np.ndarray:
        return tri_area(self.model.grid, self.base_nodes)

    def centroid(self) -> np.ndarray:
        centroid = tri_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    @property
    def base_nodes(self) -> np.ndarray:
        return self.nodes[:, :3]
    @property
    def midside_nodes(self) -> np.ndarray:
        return self.nodes[:, 3:]


class CTRAX6(AxisymmetricShellElement):
    """
    +--------+-------+-------+----+----+----+----+----+-----+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |  8 |  9  |
    +========+=======+=======+====+====+====+====+====+=====+
    | CTRAX6 |  EID  |  PID  | N1 | N2 | N3 | N4 | N5 | N6  |
    +--------+-------+-------+----+----+----+----+----+-----+
    |        | THETA |       |    |    |    |    |    |     |
    +--------+-------+-------+----+----+----+----+----+-----+

    Theta/Mcid is NX only!
    """
    def add(self, eid: int, pid: int, nids: list[int],
            theta: float=0., comment: str='') -> int:
        """Creates a CTRAX6 card"""
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CTRAX6 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4'),
            integer_or_blank(card, 7, 'n5'),
            integer_or_blank(card, 8, 'n6'),
            ]
        theta = integer_double_or_blank(card, 9, 'theta', 0.0)
        assert len(card) <= 10, f'len(CTRAX6 card) = {len(card):d}\ncard={card}'
        #return CTRAX6(eid, pid, nids, theta=theta, comment=comment)
        self.cards.append((eid, pid, nids, theta, comment))
        self.n += 1
        return self.n - 1

    @Element.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 6), dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, pid, nids, thetai, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            theta[icard] = thetai
        self._save(element_id, property_id, nodes, theta)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray,
              property_id: np.ndarray,
              nodes: np.ndarray,
              theta: np.ndarray) -> None:
        if len(self.element_id) != 0:
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            nodes = np.vstack([self.nodes, nodes])
            theta = np.hstack([self.theta, theta])
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.theta = theta

    def __apply_slice__(self, element: CTRAX6, i: np.ndarray) -> None:
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.theta = self.theta[i]
        element.n = len(self.element_id)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.property_id.max(), self.nodes.max())

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size).tolist()
        thetas = array_float(self.theta, size=size, is_double=is_double)
        for eid, pid, nodes, theta in zip_longest(element_ids, property_ids, nodes_, thetas):
            list_fields = ['CTRAX6', eid, pid] + nodes + [theta]
            bdf_file.write(print_card(list_fields))
        return

    def area(self) -> np.ndarray:
        return tri_area(self.model.grid, self.base_nodes)

    def centroid(self) -> np.ndarray:
        centroid = tri_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    @property
    def base_nodes(self) -> np.ndarray:
        return self.nodes[:, :3]
    @property
    def midside_nodes(self) -> np.ndarray:
        return self.nodes[:, 3:]
