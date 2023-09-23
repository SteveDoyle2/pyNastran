from __future__ import annotations
from abc import abstractmethod
from itertools import zip_longest
from typing import Union, Optional, Any, TYPE_CHECKING

import numpy as np
from pyNastran.bdf.field_writer_8 import print_field_8 # , print_card_8
from pyNastran.bdf.field_writer_16 import print_field_16, print_card_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, # string, # double,
    integer_or_blank, double_or_blank,
    integer_double_or_blank, blank, integer_types)
from pyNastran.bdf.cards.elements.bars import set_blank_if_default

from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, get_print_card_8_16,
    parse_element_check,
    #hslice_by_idim, make_idim,
    searchsorted_filter)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_default_int, # get_print_card,
    print_card_8_comment, print_card_16_comment)
from .utils import expanded_mass_material_id

from .shell_coords import element_coordinate_system, material_coordinate_system
from .shell_utils import (
    tri_area, tri_area_centroid_normal, tri_centroid,
    quad_area, quad_area_centroid_normal, quad_centroid,
    shell_mass_per_area, shell_mass_per_area_breakdown, shell_nonstructural_mass,
    shell_thickness,
)
from .shell_quality import tri_quality_nodes, quad_quality_nodes
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg

NUMPY_INTS = {'int32', 'int64'}
NUMPY_FLOATS = {'float32', 'float64'}


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    #from pyNastran.bdf.cards.materials import MAT1, MAT8
    #from pyNastran.dev.bdf_vectorized3.cards.grid import GRID


class ShellElement(Element):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.nodes = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.tflag = np.array([], dtype='int32')
        self.T = np.array([], dtype='int32')

    def check_types(self):
        assert self.element_id.dtype.name in NUMPY_INTS, self.element_id.dtype.name
        assert self.property_id.dtype.name in NUMPY_INTS, self.property_id.dtype.name
        assert self.tflag.dtype.name in {'int8', 'int32', 'int64'}, self.tflag.dtype.name

    @property
    def all_properties(self) -> list[Any]:
        model = self.model
        all_props = [
            model.pshell, model.pcomp,
            model.pcompg, model.plplane,
            ]  # shells
        return all_props

    @property
    def allowed_properties(self) -> list[Any]:
        all_props = self.all_properties
        props = [prop for prop in all_props if prop.n > 0]
        assert len(props) > 0, f'{self.type}: all_props={all_props}'
        return props

    def get_edge_axes(self) -> tuple[np.ndarray, np.ndarray]:
        grid = self.model.grid
        xyz = grid.xyz_cid0()
        nid = grid.node_id

        nnodes_per_element = self.nodes.shape[1]
        if nnodes_per_element in {3, 6}:
            inids = np.searchsorted(nid, self.nodes[:, :3])
            assert self.nodes.shape[1] == 3, inids.shape
            xyz1 = xyz[inids[:, 0]]
            xyz2 = xyz[inids[:, 1]]
            xyz3 = xyz[inids[:, 2]]
        elif nnodes_per_element in {4, 8}:
            inids = np.searchsorted(nid, self.nodes[:, :4])
            assert self.nodes.shape[1] == 4, inids.shape
            xyz1 = xyz[inids[:, 0], :]
            xyz2 = xyz[inids[:, 1], :]
            xyz3 = xyz[inids[:, 2], :]
            xyz4 = xyz[inids[:, 3], :]

            g12 = (xyz1 + xyz2) / 2.
            g23 = (xyz2 + xyz3) / 2.
            g34 = (xyz3 + xyz4) / 2.
            g14 = (xyz1 + xyz4) / 2.
            x = g23 - g14
            yprime = g34 - g12
            normal = np.cross(x, yprime)
            assert x.shape == normal.shape
            y = np.cross(normal, x)
            assert y.shape == normal.shape
        else:
            raise RuntimeError(nnodes_per_element)
        return x, y

    def material_coordinate_system(self) -> tuple[float,
                                                  np.ndarray, np.ndarray,
                                                  np.ndarray, np.ndarray]:
        """
        Determines the material coordinate system

        Parameters
        ----------
        normal (3, ) float ndarray
            the unit normal vector
        xyz1234 (4, 3) float ndarray
            the xyz coordinates

        Returns
        -------
        dxyz : float
            the mean length of the element
        centroid : (3, ) float ndarray
            the centroid of the element
        imat : (3, ) float ndarray
            the element unit i vector
        jmat : (3, ) float ndarray
            the element unit j vector
        normal : (3, ) float ndarray
            the unit normal vector

        .. todo:: rotate the coordinate system by the angle theta

        """
        nnodes_per_element = self.nodes.shape[1]
        if nnodes_per_element in {3, 6}:
            dxyz, centroid, normal, xyz1, xyz2 = self._dxyz_centroid_normal_xyz1_xyz2(ndim=3)
            imat, jmat = material_coordinate_system(self, normal, xyz1, xyz2)
        else:
            dxyz, centroid, normal, xyz1, xyz2 = self._dxyz_centroid_normal_xyz1_xyz2(ndim=4)
            imat, jmat = material_coordinate_system(self, normal, xyz1, xyz2)
        return dxyz, centroid, imat, jmat, normal

    def element_coordinate_system(self) -> tuple[float,
                                                 np.ndarray, np.ndarray,
                                                 np.ndarray, np.ndarray]:
        """
        Determines the element coordinate system

        Parameters
        ----------
        normal (3, ) float ndarray
            the unit normal vector
        xyz1234 (4, 3) float ndarray
            the xyz coordinates

        Returns
        -------
        dxyz : float
            the mean length of the element
        centroid : (3, ) float ndarray
            the centroid of the element
        imat : (3, ) float ndarray
            the element unit i vector
        jmat : (3, ) float ndarray
            the element unit j vector
        normal : (3, ) float ndarray
            the unit normal vector

        .. todo:: rotate the coordinate system by the angle theta

        """
        nnodes_per_element = self.nodes.shape[1]
        if nnodes_per_element in {3, 6}:
            dxyz, centroid, normal, xyz1, xyz2 = self._dxyz_centroid_normal_xyz1_xyz2(ndim=3)
            ielement, jelement = element_coordinate_system(self, normal, xyz1, xyz2)
        else:
            dxyz, centroid, normal, xyz1, xyz2 = self._dxyz_centroid_normal_xyz1_xyz2(ndim=4)
            ielement, jelement = element_coordinate_system(self, normal, xyz1, xyz2)
        return dxyz, centroid, ielement, jelement, normal

    def _dxyz_centroid_normal_xyz1_xyz2(self, ndim: int) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                                                  np.ndarray, np.ndarray]:
        normal = self.normal()  # k = kmat

        grid = self.model.grid
        xyz = grid.xyz_cid0()
        nid = grid.node_id

        inids = np.searchsorted(nid, self.nodes[:, :ndim])
        assert self.nodes.shape[1] == ndim, inids.shape
        xyz1 = xyz[inids[:, 0], :]
        xyz2 = xyz[inids[:, 1], :]
        xyz3 = xyz[inids[:, 2], :]
        if ndim == 3:
            centroid = (xyz1 + xyz2 + xyz3) / 3.

            # take the mean edge length to size the vectors in the GUI
            dxyz21 = np.linalg.norm(xyz2 - xyz1)
            dxyz32 = np.linalg.norm(xyz3 - xyz2)
            dxyz13 = np.linalg.norm(xyz1 - xyz3)
            dxyz = np.mean([dxyz21, dxyz32, dxyz13]) / 2.
        else:
            xyz4 = xyz[inids[:, 3], :]

            centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.

            # take the mean length to size the vectors in the GUI
            dxyz21 = np.linalg.norm(xyz2 - xyz1)
            dxyz32 = np.linalg.norm(xyz3 - xyz2)
            dxyz43 = np.linalg.norm(xyz4 - xyz3)
            dxyz14 = np.linalg.norm(xyz1 - xyz4)
            dxyz = np.mean([dxyz21, dxyz32, dxyz43, dxyz14]) / 2.
        return dxyz, centroid, normal, xyz1, xyz2

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

    def mass_breakdown(self) -> np.ndarray:
        """
        [area, nsm, rho, t, mass_per_area, mass]
        TODO: doesn't consider differential thickness
        """
        nelement = len(self.element_id)
        assert nelement > 0, nelement
        mass_per_area_breakdown = shell_mass_per_area_breakdown(
            self.model, self.tflag, self.T,
            self.property_id, self.allowed_properties)
        assert len(mass_per_area_breakdown) == nelement, mass_per_area_breakdown

        mass_per_area = mass_per_area_breakdown[:, -1]
        area = self.area()
        mass = mass_per_area * area

        _check_shell_mass(self, mass, area)
        breakdown = np.column_stack([area, mass_per_area_breakdown, mass])
        assert breakdown.shape[1] == 6, breakdown.shape
        return breakdown

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

    @parse_element_check
    def write_file(self, file_obj: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if size == 8:
            self.write_file_8(file_obj, write_card_header=write_card_header)
        else:
            self.write_file_16(file_obj, is_double=is_double, write_card_header=write_card_header)
        return


def _check_shell_mass(element: ShellElement, mass: np.ndarray, area=np.ndarray):
    inan = np.isnan(mass)
    if np.any(inan):
        #pids = np.unique(self.property_id[inan])
        msg = f'{element.type} has nan mass\n'
        msg += f'element_id={element.element_id[inan]}\n'
        msg += f'property_id={element.property_id[inan]}\n'
        msg += f'area={area[inan]}\n'
        mass_per_area = element.mass_per_area()
        t = element.total_thickness()
        msg += f'mass_per_area={mass_per_area[inan]}\n'
        msg += f't={t[inan]}\n'
        msg += f'all_properties={element.all_properties}\n'
        element.model.log.warning(msg)
        raise RuntimeError(msg)
        #self.model.log.warning(f'eids={self.element_id[inan]} with pids={pids} has nan mass')


class CTRIA3(ShellElement):
    """
    +--------+-------+-------+----+----+----+------------+---------+
    |   1    |   2   |   3   |  4 |  5 |  6 |     7      |    8    |
    +========+=======+=======+=====+===+====+============+=========+
    | CTRIA3 |  EID  |  PID  | N1 | N2 | N3 | THETA/MCID | ZOFFSET |
    +--------+-------+-------+----+----+----+------------+---------+
    |        |       | TFLAG | T1 | T2 | T3 |            |         |
    +--------+-------+-------+----+----+----+------------+---------+

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 3), dtype='int32')
        self.mcid = np.array([], dtype='int32')
        self.theta = np.array([], dtype='float64')
        self.zoffset = np.array([], dtype='float64')
        self.tflag = np.array([], dtype='int32')
        self.T = np.zeros((0, 3), dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int],
            theta_mcid: int|float=0.0, zoffset: float=0.,
            tflag: int=0, T1=None, T2=None, T3=None,
            comment: str='') -> int:
        self.cards.append(((eid, pid, nids,
                            theta_mcid, zoffset,
                            tflag, T1, T2, T3,
                            comment)))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')
        #: Element ID
        eid = integer(card, 1, 'eid')
        #: Property ID
        pid = integer_or_blank(card, 2, 'pid', default=eid)

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
        ]
        if len(card) > 6:
            theta_mcid = integer_double_or_blank(card, 6, 'theta_mcid', default=0.0)
            zoffset = double_or_blank(card, 7, 'zoffset', default=0.0)
            blank(card, 8, 'blank')
            blank(card, 9, 'blank')

            tflag = integer_or_blank(card, 10, 'tflag', default=0)
            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            assert len(card) <= 14, f'len(CTRIA3 card) = {len(card):d}\ncard={card}\n tflag={tflag} T123=[{T1}, {T2}, {T3}]'
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            tflag = 0
            T1 = None
            T2 = None
            T3 = None
        self.cards.append((eid, pid, nids,
                            theta_mcid, zoffset,
                            tflag, T1, T2, T3,
                            comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards

        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 3), dtype='int32')
        tflag = np.zeros(ncards, dtype='int8')
        mcid = np.full(ncards, -1, dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        zoffset = np.full(ncards, np.nan, dtype='float64')
        T = np.zeros((ncards, 3), dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, pid, nids,
             theta_mcid, zoffseti,
             tflagi, T1, T2, T3,
             comment) = card

            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            zoffset[icard] = zoffseti
            tflag[icard] = tflagi
            if isinstance(theta_mcid, float):
                theta[icard] = theta_mcid
            else:
                mcid[icard] = theta_mcid
            T[icard, :] = [T1, T2, T3]
        self._save(element_id, property_id, nodes,
                   zoffset, mcid, theta,
                   tflag, T)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, nodes,
              zoffset, mcid, theta,
              tflag, T):
        assert element_id.min() >= 0, element_id
        assert property_id.min() >= 0, property_id
        assert nodes.min() >= 0, nodes
        ncards_existing = len(self.element_id)
        if ncards_existing > 0:
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            nodes = np.vstack([self.nodes, nodes])
            mcid = np.hstack([self.mcid, mcid])
            theta = np.hstack([self.theta, theta])
            tflag = np.hstack([self.tflag, tflag])
            T = np.vstack([self.T, T])
            zoffset = np.hstack([self.zoffset, zoffset])
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.tflag = tflag
        self.mcid = mcid
        self.theta = theta
        self.zoffset = zoffset
        self.T = T

    def __apply_slice__(self, element: CTRIA3, i: np.ndarray) -> None:  # ignore[override]
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.tflag = self.tflag[i]
        element.mcid = self.mcid[i]
        element.theta = self.theta[i]
        element.zoffset = self.zoffset[i]
        element.T = self.T[i, :]
        element.n = len(self.element_id)

    def card_headers(self, size: int=8) -> list[str]:
        theta_mcid = 'th_mcid' if size == 8 else 'theta_mcid'
        headers = [
            'CTRIA3', 'eid', 'pid', 'node1', 'node2', 'node3',
            theta_mcid, 'zoffset', 'blank', 'blank', 'tflag', 'T1', 'T2', 'T3',
        ]
        return headers

    def write_file_8(self, bdf_file: TextIOLike,
                     write_card_header: bool=False) -> None:
        if len(self.element_id) == 0:
            return

        size = 8
        headers = self.card_headers()
        if write_card_header:
            bdf_file.write(print_card_8_comment(headers))
        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes_ = array_str(self.nodes, size=size)
        for eid, pid, nodes, theta, mcid, zoffset, tflag, T in zip_longest(element_id, property_id, nodes_, self.theta,
                                                                           self.mcid, self.zoffset, self.tflag, self.T):
            row1 = [eid, pid] + nodes.tolist()
            T1, T2, T3 = T
            if np.isnan(theta):
                theta_mcid = '%8d' % mcid
            else:
                theta_mcid = print_field_8(theta)

            row2_data0 = [theta_mcid, zoffset,  # actually part of line 1
                         tflag, T1, T2, T3]
            if row2_data0 == [0.0, 0.0, 0, 1.0, 1.0, 1.0]:
                msg = 'CTRIA3  %8s%8s%8s%8s%8s\n' % tuple(row1)
            else:
                zoffset = set_blank_if_default(zoffset, 0.0)
                tflag = set_blank_if_default(tflag, 0)
                #theta_mcid = self._get_theta_mcid_repr()

                T1 = set_blank_if_default(T1, 1.0)
                T2 = set_blank_if_default(T2, 1.0)
                T3 = set_blank_if_default(T3, 1.0)

                row2_data = [theta_mcid, zoffset, tflag, T1, T2, T3]
                row2 = [print_field_8(field) for field in row2_data]
                msg = ('CTRIA3  %8s%8s%8s%8s%8s%8s%8s\n'
                       '                %8s%8s%8s%8s\n' % tuple(row1 + row2)).rstrip(' \n') + '\n'
            bdf_file.write(msg)
        return

    def write_file_16(self, bdf_file: TextIOLike,
                      is_double: bool=False,
                      write_card_header: bool=False) -> None:
        if len(self.element_id) == 0:
            return
        size = 16
        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        nodes_ = array_str(self.nodes, size=size)
        headers = self.card_headers(size=size)
        if write_card_header:
            bdf_file.write(print_card_16_comment(headers))
        for eid, pid, nodes, theta, mcid, zoffset, tflag, T in zip_longest(
            element_id, property_id, nodes_, self.theta,
            self.mcid, self.zoffset, self.tflag, self.T):

            row1 = [eid, pid] + nodes.tolist()
            T1, T2, T3 = T
            if np.isnan(theta):
                theta_mcid = '%8d' % mcid
            else:
                theta_mcid = print_field_8(theta)

            row2_data0 = [theta_mcid, zoffset,  # actually part of line 1
                         tflag, T1, T2, T3]
            if row2_data0 == [0.0, 0.0, 0, 1.0, 1.0, 1.0]:
                msg = (
                    'CTRIA3* %16s%16s%16s%16s\n'
                    '*       %16s\n') % tuple(row1)
            else:
                zoffset = set_blank_if_default(zoffset, 0.0)
                tflag = set_blank_if_default(tflag, 0)
                #theta_mcid = self._get_theta_mcid_repr()

                T1 = set_blank_if_default(T1, 1.0)
                T2 = set_blank_if_default(T2, 1.0)
                T3 = set_blank_if_default(T3, 1.0)

                row2_data = [theta_mcid, zoffset, tflag, T1, T2, T3]
                row2 = [print_field_8(field) for field in row2_data]
                msg = ('CTRIA3  %8s%8s%8s%8s%8s%8s%8s\n'
                       '                %8s%8s%8s%8s\n' % tuple(row1 + row2)).rstrip(' \n') + '\n'
            bdf_file.write(msg)
        return

    def area(self):
        return tri_area(self.model.grid, self.nodes)

    def centroid(self) -> np.ndarray:
        """centroid ignores density"""
        centroid = tri_centroid(self.model.grid, self.nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        """center_of_mass considers density"""
        return self.centroid()

    def normal(self) -> np.ndarray:
        normal = self.area_centroid_normal()[2]
        return normal
    def area_centroid_normal(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        normal = tri_area_centroid_normal(self.model.grid, self.nodes)
        return normal

    @property
    def base_nodes(self) -> np.ndarray:
        return self.nodes
    @property
    def midside_nodes(self):
        return None

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

    def quality(self):
        return tri_quality_nodes(self.model.grid, self.nodes)


class CTRIAR(ShellElement):
    """
    +--------+-------+-------+----+----+----+------------+---------+
    |   1    |   2   |   3   |  4 |  5 |  6 |     7      |    8    |
    +========+=======+=======+=====+===+====+============+=========+
    | CTRIAR |  EID  |  PID  | N1 | N2 | N3 | THETA/MCID | ZOFFSET |
    +--------+-------+-------+----+----+----+------------+---------+
    |        |       | TFLAG | T1 | T2 | T3 |            |         |
    +--------+-------+-------+----+----+----+------------+---------+

    """

    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 3), dtype='int32')
        self.mcid = np.array([], dtype='int32')
        self.theta = np.array([], dtype='float64')
        self.zoffset = np.array([], dtype='float64')
        self.tflag = np.array([], dtype='int32')
        self.T = np.zeros((0, 3), dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int],
            theta_mcid: int|float=0.0, zoffset: float=0.,
            tflag: int=0, T1=None, T2=None, T3=None,
            comment: str='') -> int:
        """
        Creates a CTRIAR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 : float; default=None
            If it is not supplied, then T1 through T3 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((eid, pid, nids,
                           theta_mcid, zoffset,
                           tflag, [T1, T2, T3],
                           comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')
        #: Element ID
        #: Element ID
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
        ]

        theta_mcid = integer_double_or_blank(card, 6, 'theta_mcid', default=0.0)
        zoffset = double_or_blank(card, 7, 'zoffset', default=0.0)
        blank(card, 8, 'blank')
        blank(card, 9, 'blank')

        tflag = integer_or_blank(card, 10, 'tflag', default=0)
        T1 = double_or_blank(card, 11, 'T1')
        T2 = double_or_blank(card, 12, 'T2')
        T3 = double_or_blank(card, 13, 'T3')
        assert len(card) <= 14, f'len(CTRIAR card) = {len(card):d}\ncard={card}'

        card = (eid, pid, nids,
                theta_mcid, zoffset,
                tflag, [T1, T2, T3],
                comment)
        self.cards.append(card)
        self.n += 1
        return self.n

    def __apply_slice__(self, element: CTRIAR, i: np.ndarray) -> None:  # ignore[override]
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.tflag = self.tflag[i]
        element.mcid = self.mcid[i]
        element.theta = self.theta[i]
        element.zoffset = self.zoffset[i]
        element.T = self.T[i, :]
        element.n = len(self.element_id)

    def parse_cards(self) -> None:
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 3), dtype='int32')
        tflag = np.zeros(ncards, dtype='int8')
        mcid = np.full(ncards, -1, dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        zoffset = np.full(ncards, np.nan, dtype='float64')
        T = np.zeros((ncards, 3), dtype='float64')
        for icard, card in enumerate(self.cards):
            (eid, pid, nids,
             theta_mcid, zoffseti,
             tflagi, Ti,
             comment) = card

            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            zoffset[icard] = zoffseti
            tflag[icard] = tflagi
            if isinstance(theta_mcid, float):
                theta[icard] = theta_mcid
            else:
                mcid[icard] = theta_mcid
            T[icard, :] = Ti
        self._save(element_id, property_id, nodes,
                   zoffset, theta, mcid,
                   tflag, T)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, nodes,
              zoffset, theta, mcid,
              tflag, T):
        assert element_id.min() >= 0, element_id
        assert property_id.min() >= 0, property_id
        assert nodes.min() >= 0, nodes
        ncards_existing = len(self.element_id)
        if ncards_existing > 0:
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            nodes = np.vstack([self.nodes, nodes])
            mcid = np.hstack([self.mcid, mcid])
            theta = np.hstack([self.theta, theta])
            tflag = np.hstack([self.tflag, tflag])
            T = np.vstack([self.T, T])
            zoffset = np.hstack([self.zoffset, zoffset])
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes
        self.tflag = tflag
        self.mcid = mcid
        self.theta = theta
        self.zoffset = zoffset
        self.T = T

    @parse_element_check
    def write_file_8(self, bdf_file: TextIOLike,
                   write_card_header: bool=False) -> None:
        self.write_file(bdf_file, size=8, is_double=False,
                        write_card_header=write_card_header)

    @parse_element_check
    def write_file_16(self, bdf_file: TextIOLike,
                      is_double=False, write_card_header: bool=False) -> None:
        self.write_file(bdf_file, size=16, is_double=is_double,
                        write_card_header=write_card_header)

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        assert self.nodes.shape[1] == 3, self.nodes.shape
        print_card = get_print_card_8_16(size)

        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes = array_str(self.nodes, size=size)
        for eid, pid, nodes, theta, mcid, zoffset, tflag, T in zip_longest(
                element_ids, property_ids, nodes, self.theta,
                self.mcid, self.zoffset, self.tflag, self.T):
            if np.all(np.isnan(T)):
                T1 = T2 = T3 = None
            else:
                T1, T2, T3 = T
            if np.isnan(theta):
                theta_mcid = '%8d' % mcid
            else:
                theta_mcid = print_field_8(theta)

            #+--------+-------+-------+----+----+----+------------+---------+
            #|   1    |   2   |   3   |  4 |  5 |  6 |     7      |    8    |
            #+========+=======+=======+=====+===+====+============+=========+
            #| CTRIAR |  EID  |  PID  | N1 | N2 | N3 | THETA/MCID | ZOFFSET |
            #|        |       | TFLAG | T1 | T2 | T3 |            |         |
            #+--------+-------+-------+----+----+----+------------+---------+
            list_fields = ['CTRIAR', eid, pid] + nodes.tolist() + [
                theta_mcid, zoffset, None, None, tflag, T1, T2, T3]
            bdf_file.write(print_card(list_fields))
        return

    def area(self):
        return tri_area(self.model.grid, self.nodes)

    def centroid(self) -> np.ndarray:
        """centroid ignores density"""
        centroid = tri_centroid(self.model.grid, self.nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        """center_of_mass considers density"""
        return self.centroid()

    def normal(self) -> np.ndarray:
        normal = self.area_centroid_normal()[2]
        return normal

    def area_centroid_normal(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        normal = tri_area_centroid_normal(self.model.grid, self.nodes)
        return normal

    @property
    def base_nodes(self) -> np.ndarray:
        return self.nodes
    @property
    def midside_nodes(self):
        return None

    def flip_normal(self, i=None):
        """
        1-2-3
        3-1-2
        """
        if i is None:
            i = slice(len(self.element_id))
        self.nodes[i, :] = self.nodes[i, [2, 0, 1]]

    def quality(self):
        return tri_quality_nodes(self.model.grid, self.nodes)


class CQUAD4(ShellElement):
    """
    +--------+-------+-------+----+----+----+----+------------+---------+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |     8      |    9    |
    +========+=======+=======+=====+===+====+====+============+=========+
    | CQUAD4 |  EID  |  PID  | N1 | N2 | N3 | N4 | THETA/MCID | ZOFFSET |
    +--------+-------+-------+----+----+----+----+------------+---------+
    |        |       | TFLAG | T1 | T2 | T3 | T4 |            |         |
    +--------+-------+-------+----+----+----+----+------------+---------+

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 4), dtype='int32')
        self.mcid = np.array([], dtype='int32')
        self.theta = np.array([], dtype='float64')
        self.zoffset = np.array([], dtype='float64')
        self.tflag = np.array([], dtype='int32')
        self.T = np.zeros((0, 4), dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int],
            theta_mcid: int|float=0.0, zoffset: float=np.nan,
            tflag: int=0, T1=None, T2=None, T3=None, T4=None,
            comment: str='') -> int:
        self.cards.append((eid, pid, nids,
            theta_mcid, zoffset,
            tflag, T1, T2, T3, T4,
            comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', default=eid)
        nids = [integer(card, 3, 'n1'),
                integer(card, 4, 'n2'),
                integer(card, 5, 'n3'),
                integer(card, 6, 'n4'),]
        if len(card) > 7:
            theta_mcid = integer_double_or_blank(card, 7, 'theta_mcid', default=0.0)
            zoffset = double_or_blank(card, 8, 'zoffset', default=np.nan)
            blank(card, 9, 'blank')
            tflag = integer_or_blank(card, 10, 'tflag', default=0)
            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            T4 = double_or_blank(card, 14, 'T4')
            assert len(card) <= 15, f'len(CQUAD4 card) = {len(card):d}\ncard={card}'
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            tflag = 0
            T1 = None
            T2 = None
            T3 = None
            T4 = None

        self.cards.append((eid, pid, nids,
            theta_mcid, zoffset,
            tflag, T1, T2, T3, T4,
            comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        assert self.n >= 0, self.n
        if self.n == 0 or len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 4), dtype='int32')
        tflag = np.zeros(ncards, dtype='int8')
        mcid = np.full(ncards, -1, dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        zoffset = np.full(ncards, np.nan, dtype='float64')
        T = np.zeros((ncards, 4), dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids,
             theta_mcid, zoffseti,
             tflagi, T1, T2, T3, T4,
             comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            zoffset[icard] = zoffseti
            tflag[icard] = tflagi
            if isinstance(theta_mcid, float):
                theta[icard] = theta_mcid
            else:
                mcid[icard] = theta_mcid
            T[icard, :] = [T1, T2, T3, T4]
        self._save(element_id, property_id, nodes, zoffset, theta, mcid, tflag, T)
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray, property_id: np.ndarray, nodes: np.ndarray,
              zoffset=None, theta=None, mcid=None,
              tflag=None, T=None) -> None:
        _save_quad(self, element_id, property_id, nodes,
                   zoffset=zoffset, theta=theta, mcid=mcid, tflag=tflag, T=T)

    def __apply_slice__(self, element: CQUAD4, i: np.ndarray) -> None:  # ignore[override]
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.tflag = self.tflag[i]
        element.mcid = self.mcid[i]
        element.theta = self.theta[i]
        element.zoffset = self.zoffset[i]
        element.T = self.T[i, :]
        element.n = len(self.element_id)
        self.check_types()

    def check_types(self):
        super().check_types()
        assert self.T.dtype.name in NUMPY_FLOATS, self.T.dtype.name

    def _setup_write(self, size: int=8) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                                 np.ndarray, np.ndarray, np.ndarray]:
        self.check_types()
        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        #nodes_ = array_str(self.nodes, size=size)
        remove_tflag = (
            np.all(self.tflag == 0) and
            np.all(np.isnan(self.T))
        )
        mcids = array_default_int(self.mcid, default=-1, size=size)
        no_zoffset = np.all(np.isnan(self.zoffset))
        no_mcid = np.all(mcids == '')
        #CQUAD4    307517     105  247597  262585  262586  247591      -1     0.0
        return element_id, property_id, remove_tflag, no_zoffset, mcids, no_mcid

    def card_headers(self, size: int=8) -> list[str]:
        theta_mcid = 'th_mcid' if size == 8 else 'theta_mcid'
        headers = ['CQUAD4', 'eid', 'pid', 'node1', 'node2', 'node3', 'node4',
                   theta_mcid, 'zoffset', 'blank', 'tflag', 'T1', 'T2', 'T3', 'T4']
        return headers

    def write_file_8(self, bdf_file: TextIOLike,
                     write_card_header: bool=False) -> None:
        if len(self.element_id) == 0:
            return

        headers = self.card_headers()
        if write_card_header:
            bdf_file.write(print_card_8_comment(headers))
        element_id, property_id, remove_tflag, no_zoffset, mcids, no_mcid = self._setup_write()
        if remove_tflag:
            if no_zoffset and no_mcid:
                for eid, pid, nodes in zip_longest(element_id, property_id, self.nodes):
                    data = [eid, pid] + nodes.tolist()
                    msg = 'CQUAD4  %8s%8s%8d%8d%8d%8d\n' % tuple(data)
                    bdf_file.write(msg)
            elif no_zoffset:
                for eid, pid, nodes, theta, mcid in zip(element_id, property_id, self.nodes, self.theta, mcids):
                    data = [eid, pid, nodes[0], nodes[1], nodes[2], nodes[3], mcid]
                    msg = ('CQUAD4  %8s%8s%8d%8d%8d%8d%8s'  % tuple(data)).rstrip(' ') + '\n'
                    bdf_file.write(msg)
            elif no_mcid:
                for eid, pid, nodes, theta, zoffset in zip(element_id, property_id, self.nodes, self.theta, self.zoffset):
                    zoffset = '' if np.isnan(zoffset) else zoffset
                    data = [eid, pid, nodes[0], nodes[1], nodes[2], nodes[3], '', print_field_8(zoffset)]
                    msg = ('CQUAD4  %8s%8s%8d%8d%8d%8d%8s%8s'  % tuple(data)).rstrip(' ') + '\n'
                    bdf_file.write(msg)
            else:
                for eid, pid, nodes, theta, mcid, zoffset in zip(element_id, property_id, self.nodes, self.theta, mcids, self.zoffset):
                    zoffset = '' if np.isnan(zoffset) else zoffset
                    data = [eid, pid, nodes[0], nodes[1], nodes[2], nodes[3], mcid, print_field_8(zoffset)]
                    msg = ('CQUAD4  %8s%8s%8d%8d%8d%8d%8s%8s'  % tuple(data)).rstrip(' ') + '\n'
                    bdf_file.write(msg)
        else:
            for eid, pid, nodes, theta, mcid, zoffset, tflag, T in zip(element_id, property_id, self.nodes, self.theta,
                                                                       mcids, self.zoffset, self.tflag, self.T):
                zoffset = '' if np.isnan(zoffset) else zoffset
                T1, T2, T3, T4 = T
                if np.isnan(theta):
                    theta_mcid = '%8s' % mcid
                else:
                    theta_mcid = print_field_8(theta)

                row2_data = [theta_mcid, zoffset,  # actually part of line 1
                             tflag, T1, T2, T3, T4]
                if row2_data == [0.0, 0.0, 0, 1.0, 1.0, 1.0, 1.0]:
                    data = [eid, pid] + nodes.tolist()
                    msg = ('CQUAD4  %8s%8s%8d%8d%8d%8d\n' % tuple(data))
                    #return self.comment + msg
                else:
                    #theta_mcid = self._get_theta_mcid_repr()
                    zoffset = set_blank_if_default(zoffset, 0.0)
                    tflag = set_blank_if_default(tflag, 0)
                    T1 = set_blank_if_default(T1, 1.0)
                    T2 = set_blank_if_default(T2, 1.0)
                    T3 = set_blank_if_default(T3, 1.0)
                    T4 = set_blank_if_default(T4, 1.0)

                    row2_data = [theta_mcid, zoffset,
                                 tflag, T1, T2, T3, T4]
                    row2 = [print_field_8(field) for field in row2_data]
                    data = [eid, pid] + nodes.tolist() + row2
                    msg = ('CQUAD4  %8s%8s%8d%8d%8d%8d%8s%8s\n'
                           '                %8s%8s%8s%8s%8s\n' % tuple(data)).rstrip('\n ') + '\n'
                    #return self.comment + msg.rstrip('\n ') + '\n'
                bdf_file.write(msg)
        return

    def write_file_16(self, bdf_file: TextIOLike,
                      is_double: bool=False,
                      write_card_header: bool=False) -> None:
        if len(self.element_id) == 0:
            return
        print_card = print_card_16
        headers = self.card_headers(size=16)
        if write_card_header:
            bdf_file.write(print_card_16_comment(headers))
        element_id, property_id, remove_tflag, no_zoffset, mcids, no_mcid = self._setup_write()
        if remove_tflag:
            if no_zoffset and no_mcid:
                for eid, pid, nodes in zip_longest(element_id, property_id, self.nodes):
                    data = ['CQUAD4', eid, pid] + nodes.tolist()
                    bdf_file.write(print_card(data))
            elif no_zoffset:
                for eid, pid, nodes, theta, mcid in zip(element_id, property_id, self.nodes, self.theta, mcids):
                    data = ['CQUAD4', eid, pid, nodes[0], nodes[1], nodes[2], nodes[3], mcid]
                    bdf_file.write(print_card(data))
            elif no_mcid:
                for eid, pid, nodes, theta, zoffset in zip(element_id, property_id, self.nodes, self.theta, self.zoffset):
                    zoffset_str = '' if np.isnan(zoffset) else print_field_16(zoffset)
                    data = ['CQUAD4', eid, pid, nodes[0], nodes[1], nodes[2], nodes[3], '', zoffset]
                    bdf_file.write(print_card(data))
            else:
                for eid, pid, nodes, theta, mcid, zoffset in zip(element_id, property_id, self.nodes, self.theta, mcids, self.zoffset):
                    zoffset_str = '' if np.isnan(zoffset) else print_field_16(zoffset)
                    data = ['CQUAD4', eid, pid, nodes[0], nodes[1], nodes[2], nodes[3], mcid, zoffset_str]
                    bdf_file.write(print_card(data))
        else:
            for eid, pid, nodes, theta, mcid, zoffset, tflag, T in zip(element_id, property_id, self.nodes, self.theta,
                                                                       mcids, self.zoffset, self.tflag, self.T):
                #zoffset = '' if np.isnan(zoffset) else zoffset
                T1, T2, T3, T4 = T
                if np.isnan(theta):
                    theta_mcid = '%8s' % mcid
                else:
                    theta_mcid = print_field_8(theta)

                row2_data = [theta_mcid, zoffset,  # actually part of line 1
                             tflag, T1, T2, T3, T4]
                if row2_data == [0.0, 0.0, 0, 1.0, 1.0, 1.0, 1.0]:
                    data = [eid, pid] + nodes.tolist()
                    msg = ('CQUAD4* %16s%16s%16d%16d\n'
                           '*       %16d%16d\n' % tuple(data))
                    #return self.comment + msg
                else:
                    #theta_mcid = self._get_theta_mcid_repr()
                    zoffset = set_blank_if_default(zoffset, 0.0)
                    tflag = set_blank_if_default(tflag, 0)
                    T1 = set_blank_if_default(T1, 1.0)
                    T2 = set_blank_if_default(T2, 1.0)
                    T3 = set_blank_if_default(T3, 1.0)
                    T4 = set_blank_if_default(T4, 1.0)

                    row2_data = [theta_mcid, zoffset,
                                 tflag, T1, T2, T3, T4]
                    row2 = [print_field_16(field) for field in row2_data]
                    is_stripped = [field.strip() == '' for field in row2]
                    if all(is_stripped[2:]): # tflag, t1234 are blank
                        data = [eid, pid] + nodes.tolist() + row2[:2]
                        msg = ('CQUAD4* %16s%16s%16d%16d\n'
                               '*       %16d%16d%16s%16s\n'
                               % tuple(data))
                    else:
                        data = [eid, pid] + nodes.tolist() + row2
                        msg = ('CQUAD4* %16d%16d%16d%16d\n'
                               '*       %16d%16d%16s%16s\n'
                               '*                     %16s%16s%16s\n'
                               '*       %16s%16s\n'
                               % tuple(data)).rstrip('*\n ') + '\n'
                    #return self.comment + msg.rstrip('*\n ') + '\n'
                bdf_file.write(msg)
        return

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

    def normal(self) -> np.ndarray:
        normal = self.area_centroid_normal()[2]
        return normal

    def area_centroid_normal(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        normal = quad_area_centroid_normal(self.model.grid, self.nodes)
        return normal

    @property
    def base_nodes(self) -> np.ndarray:
        return self.nodes
    @property
    def midside_nodes(self):
        return None

    def flip_normal(self, i=None) -> None:
        r"""
        ::

          1---2       1---4         2---1
          |   |  -->  |   |   -->   |   |
          |   |       |   |         |   |
          4---3       2---3         3---4
          nominal     fast flip     preserves material orientation

        """
        if i is None:
            i = slice(len(self.element_id))
        #self.nodes[i, :] = self.nodes[i, [0, 3, 2, 1]] # fast flip
        self.nodes[i, :] = self.nodes[i, [1, 0, 3, 2]]  # preserve material orientation

    def quality(self):
        return quad_quality_nodes(self.model.grid, self.nodes)

def _set_shell(elem, eid: int, pid: int, nids: list[int],
               theta_mcid: Union[int, str], zoffset: float, tflag: int):
    elem.element_id = np.hstack([elem.element_id, eid])
    elem.property_id = np.hstack([elem.property_id, pid])
    elem.nodes = np.vstack([elem.nodes, nids])
    if isinstance(theta_mcid, integer_types):
        mcid = theta_mcid
        theta = np.nan
    else:
        mcid = -1
        theta = theta_mcid

    zoffset = np.nan if zoffset is None else zoffset
    elem.mcid = np.hstack([elem.mcid, mcid])
    elem.theta = np.hstack([elem.theta, theta])
    elem.zoffset = np.hstack([elem.zoffset, zoffset])
    elem.tflag = np.hstack([elem.tflag, tflag])


class CQUADR(ShellElement):
    """
    +--------+-------+-------+----+----+----+----+------------+---------+
    |   1    |   2   |   3   |  4 |  5 |  6 | 7  |     8      |    9    |
    +========+=======+=======+=====+===+====+====+============+=========+
    | CQUADR |  EID  |  PID  | N1 | N2 | N3 | N4 | THETA/MCID | ZOFFSET |
    +--------+-------+-------+----+----+----+----+------------+---------+
    |        |       | TFLAG | T1 | T2 | T3 | T4 |            |         |
    +--------+-------+-------+----+----+----+----+------------+---------+

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 4), dtype='int32')
        self.mcid = np.array([], dtype='int32')
        self.theta = np.array([], dtype='float64')
        self.zoffset = np.array([], dtype='float64')
        self.tflag = np.array([], dtype='int32')
        self.T = np.zeros((0, 4), dtype='float64')

    def __apply_slice__(self, element: CQUADR, i: np.ndarray) -> None:  # ignore[override]
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.tflag = self.tflag[i]
        element.mcid = self.mcid[i]
        element.theta = self.theta[i]
        element.zoffset = self.zoffset[i]
        element.T = self.T[i, :]
        element.n = len(self.element_id)

    def add(self, eid: int, pid: int, nids: list[int],
            theta_mcid: int|float=0.0, zoffset: float=0., tflag: int=0,
            T1=None, T2=None, T3=None, T4=None, comment: str='') -> int:
        """
        Creates a CQUADR card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 / T4 : float; default=None
            If it is not supplied, then T1 through T4 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        #_set_shell(self, eid, pid, nids, theta_mcid, zoffset, tflag)
        #self.element_id = np.hstack([self.element_id, eid])
        #self.property_id = np.hstack([self.property_id, pid])
        #self.nodes = np.vstack([self.nodes, nids])
        #if isinstance(theta_mcid, integer_types):
            #mcid = theta_mcid
            #theta = np.nan
        #else:
            #mcid = -1
            #theta = theta_mcid
        #self.mcid = np.hstack([self.mcid, mcid])
        #self.theta = np.hstack([self.theta, theta])
        #self.zoffset = np.hstack([self.zoffset, zoffset])
        #self.tflag = np.hstack([self.tflag, tflag])
        #self.T = np.vstack([self.T, np.array([T1, T2, T3, T4], dtype='float64')])
        self.cards.append((eid, pid, nids, theta_mcid, zoffset, tflag, [T1, T2, T3, T4]))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        nids = [integer_or_blank(card, 3, 'n1'),
                integer_or_blank(card, 4, 'n2'),
                integer_or_blank(card, 5, 'n3'),
                integer_or_blank(card, 6, 'n4'),]

        theta_mcid = integer_double_or_blank(card, 7, 'theta_mcid', default=0.0)
        zoffset = double_or_blank(card, 8, 'zoffset', default=0.0)

        tflag = integer_or_blank(card, 10, 'tflag', default=0)
        T1 = double_or_blank(card, 11, 'T1')
        T2 = double_or_blank(card, 12, 'T2')
        T3 = double_or_blank(card, 13, 'T3')
        T4 = double_or_blank(card, 14, 'T4')
        assert len(card) <= 15, f'len(CQUADR card) = {len(card):d}\ncard={card}'
        self.cards.append((eid, pid, nids, theta_mcid, zoffset, tflag, [T1, T2, T3, T4]))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        assert self.n >= 0, self.n
        if self.n == 0 or len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        self.element_id = np.zeros(ncards, dtype='int32')
        self.property_id = np.zeros(ncards, dtype='int32')
        self.nodes = np.zeros((ncards, 4), dtype='int32')
        self.tflag = np.zeros(ncards, dtype='int8')
        self.mcid = np.full(ncards, -1, dtype='int32')
        self.theta = np.full(ncards, np.nan, dtype='float64')
        self.zoffset = np.full(ncards, np.nan, dtype='float64')
        self.T = np.zeros((ncards, 4), dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, theta_mcid, zoffset, tflag, T) = card

            self.element_id[icard] = eid
            self.property_id[icard] = pid
            self.nodes[icard, :] = nids
            self.zoffset[icard] = zoffset
            self.tflag[icard] = tflag
            if isinstance(theta_mcid, float):
                self.theta[icard] = theta_mcid
            else:
                self.mcid[icard] = theta_mcid
            self.T[icard, :] = T
        self.sort()
        self.cards = []

    def _save(self, element_id: np.ndarray, property_id: np.ndarray, nodes: np.ndarray,
              zoffset=None, theta=None, mcid=None,
              tflag=None, T=None) -> None:
        _save_quad(self, element_id, property_id, nodes,
                   zoffset=zoffset, theta=theta, mcid=mcid, tflag=tflag, T=T)

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        #remove_tflag = (
            #np.all(self.tflag == 0) and
            #np.all(np.isnan(self.T))
        #)
        mcids = array_default_int(self.mcid, default=-1, size=size)
        #no_zoffset = np.all(np.isnan(self.zoffset))
        #no_mcid = np.all(mcids == '')
        #CQUAD4    307517     105  247597  262585  262586  247591      -1     0.0
        for eid, pid, nodes, theta, mcid, zoffset, tflag, T in zip_longest(
                self.element_id, self.property_id, self.nodes.tolist(), self.theta,
                mcids, self.zoffset, self.tflag, self.T):
            zoffset = '' if np.isnan(zoffset) else zoffset
            if np.isnan(theta):
                theta_mcid = '%8s' % mcid
            else:
                theta_mcid = print_field_8(theta)

            if np.all(np.isnan(T)):
                T1 = T2 = T3 = T4 = '' # , None, None, None
            else:
                T1, T2, T3, T4 = T

            list_fields = (['CQUADR', eid, pid] + nodes +
                           [theta_mcid, zoffset, None, tflag, T1, T2, T3, T4])
            bdf_file.write(print_card(list_fields))
        return

    def area(self) -> np.ndarray:
        area = quad_area(self.model.grid, self.nodes)
        return area

    def centroid(self) -> np.ndarray:
        """centroid ignores density"""
        centroid = quad_centroid(self.model.grid, self.nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        return self.centroid()

    def normal(self) -> np.ndarray:
        normal = self.area_centroid_normal()[2]
        return normal
    def area_centroid_normal(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        normal = quad_area_centroid_normal(self.model.grid, self.nodes)
        return normal

    @property
    def base_nodes(self) -> np.ndarray:
        return self.nodes
    @property
    def midside_nodes(self):
        return None

    def flip_normal(self, i=None):
        """
        1-2-3-4
        2-1-4-3
        """
        if i is None:
            i = slice(len(self.element_id))
        self.nodes[i, :] = self.nodes[i, [1, 0, 3, 2]]

    def quality(self):
        return quad_quality_nodes(self.model.grid, self.nodes)


class CTRIA6(ShellElement):
    """
    +--------+------------+---------+----+----+----+----+----+-----+
    |   1    |      2     |    3    |  4 |  5 |  6 | 7  | 8  |  9  |
    +========+============+=========+=====+===+====+====+====+=====+
    | CTRIA6 |    EID     |   PID   | N1 | N2 | N3 | N4 | N5 | N6  |
    +--------+------------+---------+----+----+----+----+----+-----+
    |        | THETA/MCID | ZOFFSET | T1 | T2 | T3 |    |    |     |
    +--------+------------+---------+----+----+----+----+----+-----+

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 6), dtype='int32')
        self.mcid = np.array([], dtype='int32')
        self.theta = np.array([], dtype='float64')
        self.zoffset = np.array([], dtype='float64')
        self.tflag = np.array([], dtype='int32')
        self.T = np.zeros((0, 3), dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int],
            theta_mcid: float=0., zoffset: float=0.,
            tflag: int=0, T1=None, T2=None, T3=None,
            comment: str='') -> int:
        """
        Creates a CTRIA6 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int/None, int/None, int/None]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 : float; default=None
            If it is not supplied, then T1 through T3 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        #self.element_id = np.hstack([self.element_id, eid])
        #self.property_id = np.hstack([self.property_id, pid])
        #self.nodes = np.vstack([self.nodes, nids])
        #if isinstance(theta_mcid, integer_types):
            #mcid = theta_mcid
            #theta = np.nan
        #else:
            #mcid = -1
            #theta = theta_mcid
        #self.mcid = np.hstack([self.mcid, mcid])
        #self.theta = np.hstack([self.theta, theta])
        #self.zoffset = np.hstack([self.zoffset, zoffset])
        #self.tflag = np.hstack([self.tflag, tflag])
        #self.T = np.vstack([self.T, [T1, T2, T3]])
        card = (eid, pid, nids, theta_mcid, zoffset,
                tflag, T1, T2, T3, comment)
        self.cards.append(card)
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CTRIA6 card from ``BDF.add_card(...)``

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
        pid = integer(card, 2, 'pid')

        nids = [
            integer(card, 3, 'n1'),
            integer(card, 4, 'n2'),
            integer(card, 5, 'n3'),
            integer_or_blank(card, 6, 'n4', default=0),
            integer_or_blank(card, 7, 'n5', default=0),
            integer_or_blank(card, 8, 'n6', default=0),
        ]
        if len(card) > 9:
            theta_mcid = integer_double_or_blank(card, 9, 'theta_mcid', default=0.0)
            zoffset = double_or_blank(card, 10, 'zoffset', default=0.0)

            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            tflag = integer_or_blank(card, 14, 'tflag', default=0)
            assert len(card) <= 15, f'len(CTRIA6 card) = {len(card):d}\ncard={card}'
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            T1 = None
            T2 = None
            T3 = None
            tflag = 0
        #return CTRIA6(eid, pid, nids, theta_mcid, zoffset,
                      #tflag, T1, T2, T3, comment=comment)
        card = (eid, pid, nids, theta_mcid, zoffset,
                tflag, T1, T2, T3, comment)
        self.cards.append(card)
        self.n += 1
        return self.n

    def __apply_slice__(self, element: CTRIA6, i: np.ndarray) -> None:  # ignore[override]
        assert element.type == 'CTRIA6'
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.tflag = self.tflag[i]
        element.mcid = self.mcid[i]
        element.theta = self.theta[i]
        element.zoffset = self.zoffset[i]
        element.T = self.T[i, :]
        element.n = len(self.element_id)

    def parse_cards(self) -> None:
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards

        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 6), dtype='int32')
        tflag = np.zeros(ncards, dtype='int8')
        mcid = np.full(ncards, -1, dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        zoffset = np.full(ncards, np.nan, dtype='float64')
        T = np.zeros((ncards, 3), dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, theta_mcid, zoffseti,
             tflagi, T1, T2, T3, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            zoffset[icard] = zoffseti
            tflag[icard] = tflagi
            if isinstance(theta_mcid, float):
                theta[icard] = theta_mcid
            else:
                mcid[icard] = theta_mcid
            T[icard, :] = [T1, T2, T3]
        self.sort()
        self.cards = []

        self._save(element_id, property_id, nodes,
                   zoffset, theta, mcid,
                   tflag, T)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, nodes,
              zoffset=None, theta=None, mcid=None,
              tflag=None, T=None):
        if len(self.element_id) != 0:
            raise NotImplementedError()
        assert element_id.min() >= 0, element_id
        assert property_id.min() >= 0, property_id
        assert nodes.min() >= 0, nodes
        nelements = len(element_id)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes

        if zoffset is None:
            zoffset = np.full(nelements, np.nan, dtype='float64')
        if theta is None:
            theta = np.full(nelements, 0., dtype='float64')
        if mcid is None:
            mcid = np.full(nelements, -1, dtype='int32')
        if tflag is None:
            tflag = np.zeros(nelements, dtype='int32')
        if T is None:
            T = np.full((nelements, 4), np.nan, dtype='float64')

        assert zoffset is not None
        self.zoffset = zoffset
        self.theta = theta
        self.mcid = mcid
        self.tflag = tflag
        self.T = T
        self.n = nelements

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)
        #remove_tflag = (
            #np.all(self.tflag == 0) and
            #np.all(np.isnan(self.T))
        #)
        element_ids = array_str(self.element_id, size=size)
        property_ids = array_str(self.property_id, size=size)
        nodes_ = array_default_int(self.nodes, default=0, size=size)
        mcids = array_default_int(self.mcid, default=-1, size=size)
        #no_zoffset = np.all(np.isnan(self.zoffset))
        #no_mcid = np.all(mcids == '')
        for eid, pid, nodes, theta, mcid, zoffset, tflag, T in zip_longest(
            element_ids, property_ids, nodes_, self.theta,
            mcids, self.zoffset, self.tflag, self.T):
            zoffset = '' if np.isnan(zoffset) else zoffset
            T1, T2, T3 = T
            if np.isnan(theta):
                theta_mcid = '%8s' % mcid
            else:
                theta_mcid = print_field_8(theta)

            zoffset = set_blank_if_default(zoffset, 0.0)
            tflag = set_blank_if_default(tflag, 0)
            T1 = set_blank_if_default(T1, 1.0)
            T2 = set_blank_if_default(T2, 1.0)
            T3 = set_blank_if_default(T3, 1.0)
            nodes2 = [None if node == 0 else node for node in nodes]
            list_fields = (['CTRIA6', eid, pid] + nodes2 +
                       [theta_mcid, zoffset, T1, T2, T3, tflag])
            bdf_file.write(print_card(list_fields))
        return

    def area(self):
        return tri_area(self.model.grid, self.base_nodes)

    def centroid(self) -> np.ndarray:
        """centroid ignores density"""
        centroid = tri_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        """center_of_mass considers density"""
        return self.centroid()

    def area_centroid_normal(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        normal = tri_area_centroid_normal(self.model.grid, self.base_nodes)
        return normal

    @property
    def base_nodes(self):
        return self.nodes[:, :3]
    @property
    def midside_nodes(self):
        return self.nodes[:, 3:]

    def flip_normal(self, i=None):
        r"""
        Flips normal of element.

        ::

               1                1
               **               **
              *  *             *  *
             4    6   -->     6    4
            *      *         *      *
           2----5---3       3----5---2

        """
        if i is None:
            i = slice(len(self.element_id))
        self.nodes[i, :] = self.nodes[i, [0, 2, 1, 5, 4, 3]]

    def quality(self):
        return tri_quality_nodes(self.model.grid, self.base_nodes)


class CQUAD8(ShellElement):
    """
    +--------+-------+-----+----+----+----+----+------------+-------+
    |    1   |   2   |  3  |  4 |  5 |  6 |  7 |      8     |   9   |
    +========+=======+=====+====+====+====+====+============+=======+
    | CQUAD8 |  EID  | PID | G1 | G2 | G3 | G4 |     G5     |  G6   |
    +--------+-------+-----+----+----+----+----+------------+-------+
    |        |   G7  | G8  | T1 | T2 | T3 | T4 | THETA/MCID | ZOFFS |
    +--------+-------+-----+----+----+----+----+------------+-------+
    |        | TFLAG |     |    |    |    |    |            |       |
    +--------+-------+-----+----+----+----+----+------------+-------+

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 8), dtype='int32')
        self.mcid = np.array([], dtype='int32')
        self.theta = np.array([], dtype='float64')
        self.zoffset = np.array([], dtype='float64')
        self.tflag = np.array([], dtype='int32')
        self.T = np.zeros((0, 4), dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int],
            theta_mcid: int|float=0.0, zoffset: float=0.,
            tflag: int=0, T1=None, T2=None, T3=None, T4=None,
            comment: str='') -> int:
        """
        Creates a CQUAD8 card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int, int/None, int/None, int/None, int/None]
            node ids
        zoffset : float; default=0.0
            Offset from the surface of grid points to the element reference
            plane.  Requires MID1 and MID2.
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        tflag : int; default=0
            0 : Ti are actual user specified thicknesses
            1 : Ti are fractions relative to the T value of the PSHELL
        T1 / T2 / T3 / T4 : float; default=None
            If it is not supplied, then T1 through T4 will be set equal
            to the value of T on the PSHELL entry.
        comment : str; default=''
            a comment for the card

        """
        #self.element_id = np.hstack([self.element_id, eid])
        #self.property_id = np.hstack([self.property_id, pid])
        #self.nodes = np.vstack([self.nodes, nids])
        #if isinstance(theta_mcid, integer_types):
            #mcid = theta_mcid
            #theta = np.nan
        #else:
            #mcid = -1
            #theta = theta_mcid
        #self.mcid = np.hstack([self.mcid, mcid])
        #self.theta = np.hstack([self.theta, theta])
        #self.zoffset = np.hstack([self.zoffset, zoffset])
        #self.tflag = np.hstack([self.tflag, tflag])
        #self.T = np.vstack([self.T, np.array([T1, T2, T3, T4], dtype='float64')])
        card = (eid, pid, nids, theta_mcid, zoffset,
                tflag, T1, T2, T3, T4, comment)
        self.cards.append(card)
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CQUAD8 card from ``BDF.add_card(...)``

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
            T1 = double_or_blank(card, 11, 'T1')
            T2 = double_or_blank(card, 12, 'T2')
            T3 = double_or_blank(card, 13, 'T3')
            T4 = double_or_blank(card, 14, 'T4')
            theta_mcid = integer_double_or_blank(card, 15, 'theta_mcid', default=0.0)
            zoffset = double_or_blank(card, 16, 'zoffset', default=0.0)
            tflag = integer_or_blank(card, 17, 'tflag', default=0)
            assert len(card) <= 18, f'len(CQUAD8 card) = {len(card):d}\ncard={card}'
        else:
            theta_mcid = 0.0
            zoffset = 0.0
            T1 = None
            T2 = None
            T3 = None
            T4 = None
            tflag = 0
        #return CQUAD8(eid, pid, nids, theta_mcid=theta_mcid, zoffset=zoffset,
                      #tflag=tflag, T1=T1, T2=T2, T3=T3, T4=T4,
                      #comment=comment)
        card = (eid, pid, nids, theta_mcid, zoffset,
                tflag, T1, T2, T3, T4, comment)
        self.cards.append(card)
        self.n += 1
        return self.n

    def __apply_slice__(self, element: CQUAD8, i: np.ndarray) -> None:  # ignore[override]
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.tflag = self.tflag[i]
        element.mcid = self.mcid[i]
        element.theta = self.theta[i]
        element.zoffset = self.zoffset[i]
        element.T = self.T[i, :]
        element.n = len(self.element_id)

    def parse_cards(self) -> None:
        assert self.n >= 0, self.n
        if self.n == 0 or len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nodes = np.zeros((ncards, 8), dtype='int32')
        tflag = np.zeros(ncards, dtype='int8')
        mcid = np.full(ncards, -1, dtype='int32')
        theta = np.full(ncards, np.nan, dtype='float64')
        zoffset = np.full(ncards, np.nan, dtype='float64')
        T = np.zeros((ncards, 4), dtype='float64')

        for icard, card in enumerate(self.cards):
            (eid, pid, nids, theta_mcid, zoffseti,
             tflagi, T1, T2, T3, T4, comment) = card
            #card, comment = card_comment

            #eid = integer(card, 1, 'eid')
            #pid = integer(card, 2, 'pid')
            #nids = [integer(card, 3, 'n1'),
                    #integer(card, 4, 'n2'),
                    #integer(card, 5, 'n3'),
                    #integer(card, 6, 'n4'),
                    #integer_or_blank(card, 7, 'n5', default=0),
                    #integer_or_blank(card, 8, 'n6', default=0),
                    #integer_or_blank(card, 9, 'n7', default=0),
                    #integer_or_blank(card, 10, 'n8', default=0),]
            #if len(card) > 11:
                #T1 = double_or_blank(card, 11, 'T1')
                #T2 = double_or_blank(card, 12, 'T2')
                #T3 = double_or_blank(card, 13, 'T3')
                #T4 = double_or_blank(card, 14, 'T4')
                #theta_mcid = integer_double_or_blank(card, 15, 'theta_mcid', default=0.0)
                #zoffseti = double_or_blank(card, 16, 'zoffset', default=0.0)
                #tflagi = integer_or_blank(card, 17, 'tflag', default=0)
                #assert len(card) <= 18, f'len(CQUAD8 card) = {len(card):d}\ncard={card}'
            #else:
                #theta_mcid = 0.0
                #zoffseti = 0.0
                #T1 = None
                #T2 = None
                #T3 = None
                #T4 = None
                #tflagi = 0

            element_id[icard] = eid
            property_id[icard] = pid
            nodes[icard, :] = nids
            zoffset[icard] = zoffseti
            tflag[icard] = tflagi
            if isinstance(theta_mcid, float):
                theta[icard] = theta_mcid
            else:
                mcid[icard] = theta_mcid
            T[icard, :] = [T1, T2, T3, T4]

        self._save(element_id, property_id, nodes,
                   zoffset, mcid, theta,
                   tflag, T)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, nodes,
              zoffset=None, theta=None, mcid=None,
              tflag=None, T=None):
        if len(self.element_id) != 0:
            raise NotImplementedError()
        assert element_id.min() >= 0, element_id
        assert property_id.min() >= 0, property_id
        assert nodes.min() >= 0, nodes
        nelements = len(element_id)
        self.element_id = element_id
        self.property_id = property_id
        self.nodes = nodes

        if zoffset is None:
            zoffset = np.full(nelements, np.nan, dtype='float64')
        if theta is None:
            theta = np.full(nelements, 0., dtype='float64')
        if mcid is None:
            mcid = np.full(nelements, -1, dtype='int32')
        if tflag is None:
            tflag = np.zeros(nelements, dtype='int32')
        if T is None:
            T = np.full((nelements, 4), np.nan, dtype='float64')

        assert zoffset is not None
        self.zoffset = zoffset
        self.theta = theta
        self.mcid = mcid
        self.tflag = tflag
        self.T = T
        self.n = nelements

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)
        #remove_tflag = (
            #np.all(self.tflag == 0) and
            #np.all(np.isnan(self.T))
        #)
        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        mcids = array_default_int(self.mcid, default=-1, size=size)
        #no_zoffset = np.all(np.isnan(self.zoffset))
        #no_mcid = np.all(mcids == '')
        #CQUAD4    307517     105  247597  262585  262586  247591      -1     0.0
        for eid, pid, nodes, theta, mcid, zoffset, tflag, T in zip_longest(
            element_id, property_id, self.nodes, self.theta,
            mcids, self.zoffset, self.tflag, self.T):
            zoffset = '' if np.isnan(zoffset) else zoffset
            T1, T2, T3, T4 = T
            if np.isnan(theta):
                theta_mcid = '%8s' % mcid
            else:
                theta_mcid = print_field_8(theta)

            zoffset = set_blank_if_default(zoffset, 0.0)
            tflag = set_blank_if_default(tflag, 0)
            T1 = set_blank_if_default(T1, 1.0)
            T2 = set_blank_if_default(T2, 1.0)
            T3 = set_blank_if_default(T3, 1.0)
            T4 = set_blank_if_default(T4, 1.0)
            nodes2 = [None if node == 0 else node for node in nodes]
            list_fields = ['CQUAD8', eid, pid] + nodes2 + [
                T1, T2, T3, T4, theta_mcid, zoffset, tflag]
            bdf_file.write(print_card(list_fields))
        return

    def area(self) -> np.ndarray:
        area = quad_area(self.model.grid, self.base_nodes)
        return area

    def centroid(self) -> np.ndarray:
        """centroid ignores density"""
        centroid = quad_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        """center_of_mass considers density"""
        return self.centroid()

    def area_centroid_normal(self) -> np.ndarray:
        normal = quad_area_centroid_normal(self.model.grid, self.base_nodes)
        return normal

    @property
    def base_nodes(self):
        return self.nodes[:, :4]
    @property
    def midside_nodes(self):
        return self.nodes[:, 4:]

    def flip_normal(self, i=None) -> None:
        r"""
        ::

          1--5--2       1--8--4
          |     |  -->  |     |
          8     6       5     7
          |     |       |     |
          4--7--3       2--6--3

        """
        if i is None:
            i = slice(len(self.element_id))
        self.nodes[i, :] = self.nodes[i, [0, 3, 2, 1, 7, 6, 5, 4]]

    def quality(self):
        return quad_quality_nodes(self.model.grid, self.base_nodes)


class CQUAD(ShellElement):
    """
    +-------+-------+-----+----+------------+----+----+----+----+
    |    1  |   2   |  3  |  4 |     5      |  6 |  7 | 8  |  9 |
    +=======+=======+=====+====+============+====+====+====+====+
    | CQUAD |  EID  | PID | G1 |     G2     | G3 | G4 | G5 | G6 |
    +-------+-------+-----+----+------------+----+----+----+----+
    |       |   G7  | G8  | G9 | THETA/MCID |    |    |    |    |
    +-------+-------+-----+----+------------+----+----+----+----+

    theta_mcid is an MSC specific variable

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 9), dtype='int32')
        self.mcid = np.array([], dtype='int32')
        self.theta = np.array([], dtype='float64')

    def add(self, eid: int, pid: int, nids: list[int],
            theta_mcid: int|float=0., comment: str='') -> int:
        """
        Creates a CQUAD card

        Parameters
        ----------
        eid : int
            element id
        pid : int
            property id (PSHELL/PCOMP/PCOMPG)
        nids : list[int, int, int, int, int/None, int/None,
                    int/None, int/None, int/None]
            node ids
        theta_mcid : float; default=0.0
            float : material coordinate system angle (theta) is defined
                    relative to the element coordinate system
            int : x-axis from material coordinate system angle defined by
                  mcid is projected onto the element
        comment : str; default=''
            a comment for the card

        """
        self.element_id = np.hstack([self.element_id, eid])
        self.property_id = np.hstack([self.property_id, pid])
        self.nodes = np.vstack([self.nodes, nids])
        if isinstance(theta_mcid, integer_types):
            mcid = theta_mcid
            theta = np.nan
        else:
            mcid = -1
            theta = theta_mcid

        self.mcid = np.hstack([self.mcid, mcid])
        self.theta = np.hstack([self.theta, theta])
        self.n += 1
        return self.n

    def __apply_slice__(self, element: CQUAD, i: np.ndarray) -> None:  # ignore[override]
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.mcid = self.mcid[i]
        element.theta = self.theta[i]
        element.n = len(self.element_id)

    def parse_cards(self) -> None:
        assert self.n >= 0, self.n
        if self.n == 0 or len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        self.element_id = np.zeros(ncards, dtype='int32')
        self.property_id = np.zeros(ncards, dtype='int32')
        self.nodes = np.zeros((ncards, 9), dtype='int32')
        self.mcid = np.full(ncards, -1, dtype='int32')
        self.theta = np.full(ncards, np.nan, dtype='float64')

        for icard, card_comment in enumerate(self.cards):
            card, comment = card_comment

            eid = integer(card, 1, 'eid')
            pid = integer(card, 2, 'pid')
            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2'),
                    integer(card, 5, 'n3'),
                    integer(card, 6, 'n4'),
                    integer_or_blank(card, 7, 'n5', default=0),
                    integer_or_blank(card, 8, 'n6', default=0),
                    integer_or_blank(card, 9, 'n7', default=0),
                    integer_or_blank(card, 10, 'n8', default=0),
                    integer_or_blank(card, 11, 'n9', default=0),]
            theta_mcid = integer_double_or_blank(card, 12, 'theta_mcid', default=0.)
            assert len(card) <= 13, f'len(CQUAD card) = {len(card):d}\ncard={card}'

            self.element_id[icard] = eid
            self.property_id[icard] = pid
            self.nodes[icard, :] = nids
            if isinstance(theta_mcid, float):
                self.theta[icard] = theta_mcid
            else:
                self.mcid[icard] = theta_mcid
        self.sort()
        self.cards = []

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:

        element_id = array_str(self.element_id, size=size)
        mcids = array_default_int(self.mcid, default=-1, size=size)
        for eid, pid, nodes, theta, mcid in zip_longest(element_id, self.property_id, self.nodes,
                                                        self.theta, mcids):
            if np.isnan(theta):
                theta_mcid = '%8s' % mcid
            else:
                theta_mcid = print_field_8(theta)

            nodes2 = ['' if node is None else '%8d' % node for node in nodes[4:]]

            data = [eid, pid] + nodes[:4].tolist() + nodes2 + [theta_mcid]
            msg = ('CQUAD   %8s%8i%8i%8i%8i%8i%8s%8s\n'  # 6 nodes
                   '        %8s%8s%8s%8s\n' % tuple(data))
            bdf_file.write(msg)
        return

    def area(self) -> np.ndarray:
        area = quad_area(self.model.grid, self.base_nodes)
        return area

    def centroid(self) -> np.ndarray:
        """centroid ignores density"""
        centroid = quad_centroid(self.model.grid, self.base_nodes)
        return centroid

    def center_of_mass(self) -> np.ndarray:
        """center_of_mass considers density"""
        return self.centroid()

    @property
    def base_nodes(self):
        return self.nodes[:, :4]
    @property
    def midside_nodes(self):
        return self.nodes[:, 4:]

    def flip_normal(self, i: Optional[np.ndarray]=None) -> None:
        r"""
        ::

          1--5--2       1--8--4
          |     |  -->  |     |
          8  9  6       5  9  7
          |     |       |     |
          4--7--3       2--6--3

        """
        if i is None:
            i = slice(len(self.element_id))
        self.nodes[i, :] = self.nodes[i, [0, 3, 2, 1, 7, 6, 5, 4, 8]]

    def quality(self):
        return quad_quality_nodes(self.model.grid, self.base_nodes)


def _save_quad(element: CQUAD4 | CQUADR,
               element_id: np.ndarray, property_id: np.ndarray, nodes: np.ndarray,
               zoffset=None, theta=None, mcid=None,
               tflag=None, T=None):
    assert element_id.min() >= 0, element_id
    assert property_id.min() >= 0, property_id
    assert nodes.min() >= 0, nodes
    nelements = len(element_id)
    element.element_id = element_id
    element.property_id = property_id
    element.nodes = nodes

    if zoffset is None:
        zoffset = np.full(nelements, np.nan, dtype='float64')
    if theta is None:
        theta = np.full(nelements, 0., dtype='float64')
    if mcid is None:
        mcid = np.full(nelements, -1, dtype='int32')
    if tflag is None:
        tflag = np.zeros(nelements, dtype='int32')
    if T is None:
        T = np.full((nelements, 4), np.nan, dtype='float64')

    assert zoffset is not None
    element.zoffset = zoffset
    element.theta = theta
    element.mcid = mcid
    element.tflag = tflag
    element.T = T
    element.n = nelements
