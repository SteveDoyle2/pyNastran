from __future__ import annotations
from abc import abstractmethod
from itertools import count, zip_longest
from typing import Union, Optional, Any, TYPE_CHECKING

import numpy as np
from pyNastran.bdf.field_writer_8 import print_field_8, print_card_8
from pyNastran.bdf.field_writer_16 import print_field_16, print_card_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, # string, # double,
    integer_or_blank, double_or_blank,
    integer_double_or_blank, string_or_blank, blank, integer_types)
from pyNastran.bdf.cards.elements.bars import set_blank_if_default
from pyNastran.bdf.cards.properties.shell import map_failure_theory_int

from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Element, Property, get_print_card_8_16,
    hslice_by_idim, make_idim, searchsorted_filter)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_default_int, get_print_card,
    print_card_8_comment, print_card_16_comment)
from .utils import get_density_from_material, get_density_from_property, expanded_mass_material_id

from .shell_coords import element_coordinate_system, material_coordinate_system
from .shell_utils import (
    tri_area, tri_area_centroid_normal, tri_centroid,
    quad_area, quad_area_centroid_normal, quad_centroid)
from .shell_quality import tri_quality_nodes, quad_quality_nodes
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg

NUMPY_INTS = {'int32', 'int64'}
NUMPY_FLOATS = {'float32', 'float64'}


if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.bdf.cards.materials import MAT1, MAT8
    #from pyNastran.dev.bdf_vectorized3.cards.grid import GRID


def shell_materials(model: BDF) -> list[Union[MAT1, MAT8]]:
    if model.is_thermal:
        return [model.mat4, model.mat5]
    return [model.mat1, model.mat2, model.mat8, model.mat9]


class PSHELL(Property):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')

        self.material_id = np.zeros((0, 4), dtype='int32')
        self.t = np.array([], dtype='float64')
        self.twelveIt3 = np.array([], dtype='float64')
        self.tst = np.array([], dtype='float64')
        self.nsm = np.array([], dtype='float64')
        self.z = np.zeros((0, 2), dtype='float64')

    def add(self, pid: int, mid1: int=None, t: float=None,
            mid2: int=None, twelveIt3: float=1.0,
            mid3: int=None, tst: float=0.833333, nsm: float=0.0,
            z1: float=None, z2: float=None, mid4: int=None,
            comment: str='') -> PSHELL:
        """
        Creates a PSHELL card

        Parameters
        ----------
        pid : int
            property id
        mid1 : int; default=None
            defines membrane material
            defines element density (unless blank)
        mid2 : int; default=None
            defines bending material
            defines element density if mid1=None
        mid3 : int; default=None
            defines transverse shear material
        mid4 : int; default=None
            defines membrane-bending coupling material
        twelveIt3 : float; default=1.0
            Bending moment of inertia ratio, 12I/T^3. Ratio of the actual
            bending moment inertia of the shell, I, to the bending
            moment of inertia of a homogeneous shell, T^3/12. The default
            value is for a homogeneous shell.
        nsm : float; default=0.0
            non-structural mass per unit area
        z1 / z2 : float; default=None
            fiber distance location 1/2 for stress/strain calculations
            z1 default : -t/2 if thickness is defined
            z2 default : t/2 if thickness is defined
        comment : str; default=''
            a comment for the card

        """
        if z1 is None and t is not None:
            z1 = -t / 2.
        if z2 is None and t is not None:
            z2 = t / 2.
        t = np.nan if t is None else t

        self.cards.append((pid, mid1, t,
                           mid2, twelveIt3, mid3, tst, nsm, z1, z2, mid4,
                           comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        pid = integer(card, 1, 'pid')
        mid1 = integer_or_blank(card, 2, 'mid1', default=-1)
        t = double_or_blank(card, 3, 't', default=np.nan)

        mid2 = integer_or_blank(card, 4, 'mid2', default=-1)
        twelveIt3 = double_or_blank(card, 5, '12*I/t^3', 1.0)  # poor name
        mid3 = integer_or_blank(card, 6, 'mid3', default=-1)
        tst = double_or_blank(card, 7, 'ts/t', 0.833333)
        nsm = double_or_blank(card, 8, 'nsm', 0.0)

        if np.isnan(t):
            z1 = double_or_blank(card, 9, 'z1', default=np.nan)
            z2 = double_or_blank(card, 10, 'z2', default=np.nan)
        else:
            t_over_2 = t / 2.
            z1 = double_or_blank(card, 9, 'z1', -t_over_2)
            z2 = double_or_blank(card, 10, 'z2', t_over_2)
        mid4 = integer_or_blank(card, 11, 'mid4', default=-1)
        assert len(card) <= 12, f'len(PSHELL card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, mid1, t,
                           mid2, twelveIt3, mid3, tst, nsm, z1, z2, mid4,
                           comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards

        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros((ncards, 4), dtype='int32')
        t = np.zeros(ncards, dtype='float64')
        twelveIt3 = np.zeros(ncards, dtype='float64')
        tst = np.zeros(ncards, dtype='float64')
        nsm = np.zeros(ncards, dtype='float64')
        z = np.zeros((ncards, 2), dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, mid1, ti,
             mid2, twelveIt3i, mid3, tsti, nsmi, z1, z2, mid4,
             comment) = card

            property_id[icard] = pid
            material_id[icard, :] = [mid if isinstance(mid, integer_types) else -1
                                     for mid in [mid1, mid2, mid3, mid4]]
            t[icard] = ti
            twelveIt3[icard] = twelveIt3i
            tst[icard] = tsti
            nsm[icard] = nsmi
            z[icard] = [z1, z2]
        self._save(property_id, material_id, t, twelveIt3, tst, nsm, z)
        self.sort()
        self.cards = []

    def _save(self, property_id, material_id, t, twelveIt3, tst, nsm, z) -> None:
        if len(self.property_id) != 0:
            raise NotImplementedError()
        self.property_id = property_id
        self.material_id = material_id
        self.t = t
        self.twelveIt3 = twelveIt3
        self.tst = tst
        self.nsm = nsm
        self.z = z
        self.n = len(property_id)

    def __apply_slice__(self, prop: PSHELL, i: np.ndarray) -> None:
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i, :]
        prop.t = self.t[i]
        prop.twelveIt3 = self.twelveIt3[i]
        prop.tst = self.tst[i]
        prop.nsm = self.nsm[i]
        prop.z = self.z[i, :]

    @property
    def z1(self) -> np.ndarray:
        return self.z[:, 0]
    @z1.setter
    def z1(self, z1: np.ndarray) -> None:
        self.z[:, 0] = z1

    @property
    def z2(self) -> np.ndarray:
        return self.z[:, 1]
    @z2.setter
    def z2(self, z2: np.ndarray) -> None:
        self.z[:, 1] = z2

    def _filter_pcomp_pshells(self) -> None:
        pshell_pids = self.property_id
        if len(pshell_pids) == 0:
            return
        model = self.model
        pcomp_pids_ = [propi.property_id for propi in [model.pcomp, model.pcompg]
                       if propi.n > 0]
        if not pcomp_pids_:
            return
        pcomp_pids = np.hstack(pcomp_pids_)
        pshells_pids_to_remove = np.intersect1d(pshell_pids, pcomp_pids)
        if not len(pshells_pids_to_remove):
            return
        model.log.warning(f'removing PSHELL property_ids={pshells_pids_to_remove}')
        pshells_to_keep = np.setdiff1d(pshell_pids, pshells_pids_to_remove)
        i = np.searchsorted(pshell_pids, pshells_to_keep)
        self.__apply_slice__(self, i)

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          msg=f'no shell materials for {self.type}')
        mids.sort()
        material_ids = np.unique(self.material_id.ravel())
        if -1 == material_ids[0]:
            material_ids = material_ids[1:]
        geom_check(self,
                   missing,
                   material_id=(mids, material_ids))

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
              write_card_header: bool=False) -> None:
        if len(self.property_id) == 0:
            return

        max_int = max(self.property_id.max(), self.material_id.max())
        print_card = get_print_card(size, max_int)

        for pid, mids, t, twelveIt3, tst, nsm, z in zip_longest(self.property_id, self.material_id, self.t,
                                                                self.twelveIt3, self.tst, self.nsm, self.z):
            mid1, mid2, mid3, mid4 = mids
            z1, z2 = z
            twelveIt3 = set_blank_if_default(twelveIt3, 1.0)
            tst = set_blank_if_default(tst, 0.833333)
            tst2 = set_blank_if_default(tst, 0.83333)
            if tst is None or tst2 is None:
                tst = None
            nsm = set_blank_if_default(nsm, 0.0)
            if t is not None:
                t_over_2 = t / 2.
                z1b = set_blank_if_default(z1, -t_over_2)
                z2b = set_blank_if_default(z2, t_over_2)
            else:
                z1b = z1
                z2b = z2

            mid1b = None if mid1 == -1 else mid1
            mid2b = None if mid2 == -1 else mid2
            mid3b = None if mid3 == -1 else mid3
            mid4b = None if mid4 == -1 else mid4
            list_fields = ['PSHELL', pid, mid1b, t, mid2b,
                           twelveIt3, mid3b, tst, nsm, z1b, z2b, mid4b]
            msg = print_card(list_fields)
            bdf_file.write(msg)
        return

    @property
    def allowed_materials(self) -> list[Any]:
        return [mat for mat in shell_materials(self.model) if mat.n]

    def expanded_mass_material_id(self) -> tuple[np.ndarray, np.ndarray]:
        mid = self.mass_material_id()
        return self.property_id, mid

    def mass_material_id(self) -> np.ndarray:
        mid1 = self.material_id[:, 0]
        mid2 = self.material_id[:, 1]
        assert mid1.dtype.name != 'object', mid1
        assert mid2.dtype.name != 'object', mid2

        imid_null = np.where(mid1 == -1)[0]
        mid = mid1.copy()
        mid[imid_null] = mid2[imid_null]
        return mid

    def mass_per_area(self) -> np.ndarray:
        thickness = self.t
        nsm = self.nsm
        mid = self.mass_material_id()
        rho = get_density_from_material(mid, self.allowed_materials)
        mass_per_area = nsm + rho * thickness
        return mass_per_area

    def nsm_rho_thickness(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        thickness = self.t
        nsm = self.nsm
        mid = self.mass_material_id()
        rho = get_density_from_material(mid, self.allowed_materials)
        return nsm, rho, thickness

    def density(self) -> np.ndarray:
        mid = self.mass_material_id()
        rho = get_density_from_material(mid, self.allowed_materials)
        return rho

    def total_thickness(self) -> np.ndarray:
        return self.t

    def get_individual_ABD_matrices(
            self, theta_offset: float=0.) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Gets the ABD matrix

        Parameters
        ----------
        theta_offset : float
            rotates the ABD matrix; measured in degrees

        http://www2.me.rochester.edu/courses/ME204/nx_help/index.html#uid:id503291
        Understanding Classical Lamination Theory
        """
        #mids = self.get_material_ids()
        thickness = self.total_thickness()

        #z0 = self.z1
        #z1 = self.z2
        #zmean = (z0 + z1) / 2.
        #dz = z1 - z0
        #dzsquared = z1 ** 2 - z0 ** 2
        #zcubed = z1 ** 3 - z0 ** 3
        # A11 A12 A16
        # A12 A22 A26
        # A16 A26 A66
        nprop = len(self.property_id)
        A = np.zeros((nprop, 3, 3), dtype='float64')
        B = np.zeros((nprop, 3, 3), dtype='float64')
        D = np.zeros((nprop, 3, 3), dtype='float64') # TODO: 2x2 matrix?

        theta = np.zeros(nprop, dtype='float64')
        z0 = self.z[:, 0]
        z1 = self.z[:, 1]

        #if self.mid1_ref:
        mid1 = self.material_id[:, 0]
        mid2 = self.material_id[:, 1]
        mid4 = self.material_id[:, 3]
        Qbar1 = self.get_Qbar_matrix(mid1, theta=theta)
        A += Qbar1 * thickness

        #if self.mid2_ref:
        Qbar2 = self.get_Qbar_matrix(mid2, theta=theta)
        D += Qbar2 * (z1 ** 3 - z0 ** 3) * self.twelveIt3
        #Qbar3 = self.get_Qbar_matrix(self.mid3_ref, theta=0.)
        #if self.mid4_ref:
        Qbar4 = self.get_Qbar_matrix(mid4, theta=theta)
        B += Qbar4 * (z1 ** 2 - z0 ** 2)

        unused_ts = self.tst * thickness

        # [N, M, Q].T =   [TG1, T^2 * G4, 0]              * [epsilon0]
                        # [T^2 * G4, T^3/12 * G2, 0]        [xi]
                        # [0, 0, Ts * G3]                   [gamma]

        #B += Qbar * thickness * zmean
        #D += Qbar * thickness * (z1i ** 3 - z0i ** 3)
        #N += Qbar * alpha * thickness
        #M += Qbar * alpha * thickness * zmean
        #B /= 2.
        #D /= 3.
        #M /= 2.
        #np.set_printoptions(linewidth=120, suppress=True)
        #print(ABD)
        return A, B, D

    def get_Qbar_matrix(self, mid, theta: float=0.):
        """theta must be in radians"""
        S2 = get_mat_s33(self.model, mid, self.allowed_materials)
        T = get_2d_plate_transform(theta)
        #Tt = np.transpose(T, axes=0)
        #Tinv = np.linalg.inv(T)
        #Qbar = np.linalg.multi_dot([Tinv, Q, Tinv.T])
        #ST = np.einsum('ijk,ilm->ijm', S2, T)
        #Sbarv = np.einsum('ijk,ilm->ikm', T, ST)
        Sbar = [np.linalg.multi_dot([Ti.T, Si, Ti]) for Ti, Si in zip(T, S2)]
        Qbar = [np.linalg.inv(Sbari) for Sbari in Sbar]
        #Ti = T[0, :, :]
        #Si = S2[0, :, :]
        #Sbar_test = np.linalg.multi_dot([Ti.T, Si, Ti])
        #Qbarv = np.linalg.inv(Sbarv)
        Qbarv = np.array(Qbar, dtype=S2.dtype)
        return Qbarv

    def get_Sbar_matrix(self, mid_ref, theta=0.):
        """theta must be in radians"""
        # this is the inverse of Sbar
        S2, unused_S3 = get_mat_props_S(mid_ref)
        T = get_2d_plate_transform(theta)
        #Tinv = np.linalg.inv(T)
        Sbar = np.linalg.multi_dot([T.T, S2, T])
        return Sbar

    def get_ABD_matrices(self, theta_offset=0.) -> np.ndarray:
        """
        Gets the ABD matrix

        Parameters
        ----------
        theta_offset : float
            rotates the ABD matrix; measured in degrees

        """
        A, B, D = self.get_individual_ABD_matrices(theta_offset=theta_offset)
        ABD = np.block([
            [A, B],
            [B, D],
        ])
        return ABD

    def get_Ainv_equivalent_pshell(self,
                                   imat_rotation_angle: float,
                                   thickness: float) -> tuple[float, float, float, float]:
        """imat_rotation_angle is in degrees...but is specified in radians and unused"""
        ABD = self.get_ABD_matrices(imat_rotation_angle)
        A = ABD[:3, :3]
        Ainv = np.linalg.inv(A)

        # equivalent compliance matrix
        S = thickness * Ainv
        #print(S[0, 0], S[1, 1], S[2, 2])

        # these are on the main diagonal
        mid_ref = self.mid1_ref
        Ex = mid_ref.e
        Ey = Ex
        Gxy = mid_ref.g

        #Ex = 1. / S[0, 0]
        #Ey = 1. / S[1, 1]
        #Gxy = 1. / S[2, 2]
        # S12 = -nu12 / E1 = -nu21/E2
        # nu12 = -S12 / E1
        nu_xy = -S[0, 1] / Ex
        return Ex, Ey, Gxy, nu_xy

def get_mat_s33(model: BDF, material_id: np.ndarray, allowed_materials: list[Any]):
    nmaterial = len(material_id)
    s33 = np.full((nmaterial, 3, 3), np.nan, dtype='float64')
    for mat in allowed_materials:
        ilookup, iall = searchsorted_filter(mat.material_id, material_id)
        if len(iall) == 0:
            continue
        s33i = mat.s33()
        s33[ilookup, :, :] = s33i[iall, :, :]
    return s33


class PLPLANE(Property):
    """
    Fully Nonlinear Plane Element Properties (SOL 601)
    Defines the properties of a fully nonlinear
    (i.e., large strain and large rotation)
    hyperelastic plane strain or axisymmetric element.

    +---------+-----+-----+-----+-----+
    |    1    |  2  |  3  |  4  |  5  |
    +=========+=====+=====+=====+=====+
    | PLPLANE | PID | MID | CID | STR |
    +---------+-----+-----+-----+-----+

    MSC

    +---------+-----+-----+-----+-----+---+
    |    1    |  2  |  3  |  4  |  5  | 6 |
    +=========+=====+=====+=====+=====+===+
    | PLPLANE | PID | MID | CID | STR | T |
    +---------+-----+-----+-----+-----+---+

    NX

    Referenced by:
     #- CQUAD, CQUAD4, CQUAD8, CQUADX, CTRIA3, CTRIA6, CTRIAX (MSC)
     #- CPLSTS3, CPLSTS4, CPLSTS6, CPLSTS8 entries (NX 10)

    NX:
    PLPLANE can be referenced by CQUAD, CQUAD4, CQUAD8, CTRIA3, and
    CTRIA6 entries in solutions 106 and 129. The CQUAD4, CQUAD8, CTRIA3,
    and CTRIA6 entries are treated as hyperelastic plane strain elements.
    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.coord_id = np.array([], dtype='int32')
        self.stress_strain_output_location = np.array([], dtype='|U4')
        self.thickness = np.array([], dtype='float64')

    def add(self, pid: int, mid: int, cid: int=0,
            stress_strain_output_location: str='GRID', thickness: float=np.nan,
            comment: str='') -> int:
        """Creates a PLPLANE card"""
        self.cards.append((pid, mid, cid, stress_strain_output_location, thickness, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PLPLANE card from ``BDF.add_card(...)``

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
        cid = integer_or_blank(card, 3, 'cid', default=0)
        stress_strain_output_location = string_or_blank(card, 4, 'str', default='GRID')

        # Thickness for plane stress elements (SOL 601 only). If T is blank or zero
        # then the thickness must be specified for Ti on the CPLSTS3, CPLSTS4,
        # CPLSTS6, and CPLSTS8 entries. T is ignored for plane strain elements.
        # (Real >= 0.0 or blank)
        thickness = double_or_blank(card, 5, 'thickness', default=np.nan)
        assert len(card) <= 6, f'len(PLPLANE card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, mid, cid, stress_strain_output_location, thickness, comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        assert self.n >= 0, self.n
        if len(self.cards) == 0:
            return
        ncards = len(self.cards)
        assert ncards > 0, ncards
        self.property_id = np.zeros(ncards, dtype='int32')
        self.material_id = np.zeros(ncards, dtype='int32')
        self.coord_id = np.zeros(ncards, dtype='int32')
        self.stress_strain_output_location = np.zeros(ncards, dtype='|U4')
        self.thickness = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, mid, cid, stress_strain_output_location, thickness, comment) = card
            assert len(stress_strain_output_location) <= 4, stress_strain_output_location

            self.property_id[icard] = pid
            self.material_id[icard] = mid
            self.coord_id[icard] = cid
            self.stress_strain_output_location[icard] = stress_strain_output_location
            self.thickness[icard] = thickness
        self.sort()
        self.cards = []

    def __apply_slice__(self, prop: PLPLANE, i: np.ndarray) -> None:
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]
        prop.coord_id = self.coord_id[i]
        prop.stress_strain_output_location = self.stress_strain_output_location[i]
        prop.thickness = self.thickness[i]

    #def write_file_8(self, bdf_file: TextIOLike,
    #               write_card_header: bool=False) -> None:
    #    self.write_file(bdf_file, size=8, is_double=False,
    #                    write_card_header=write_card_header)
#
#    def write_file_16(self, bdf_file: TextIOLike,
#                      is_double=False, write_card_header: bool=False) -> None:
#        self.write_file(bdf_file, size=16, is_double=is_double,
#                        write_card_header=write_card_header)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.property_id) == 0:
            return
        print_card = get_print_card_8_16(size)
        pids = array_str(self.property_id, size=size)
        mids = array_str(self.material_id, size=size)
        cids = array_str(self.coord_id, size=size)
        inan = np.isnan(self.thickness)

        if np.any(inan):
            for pid, mid, cid, stress_strain_output_location in zip_longest(
                pids[inan], mids[inan], cids[inan],
                self.stress_strain_output_location[inan]):
                list_fields = ['PLPLANE', pid, mid, cid,
                               stress_strain_output_location]
                msg = print_card(list_fields)
                bdf_file.write(msg)

        not_nan = ~inan
        if np.any(not_nan):
            for pid, mid, cid, stress_strain_output_location, t in zip_longest(
                pids[not_nan], mids[not_nan], cids[not_nan],
                self.stress_strain_output_location[not_nan],
                self.thickness[not_nan]):
                list_fields = ['PLPLANE', pid, mid, cid,
                               stress_strain_output_location, t]
                msg = print_card(list_fields)
                bdf_file.write(msg)
        return

    @property
    def allowed_materials(self) -> list[Any]:
        return []
        #return [mat for mat in shell_materials(self.model) if mat.n]

    def mass_per_area(self) -> np.ndarray:
        #thickness = self.thickness
        thickness = np.zeros(len(self.property_id))
        #nsm = self.nsm
        #mid = self.material_id

        #rho = get_density_from_material(mid, self.allowed_materials)
        #mass_per_area = rho * thickness
        #print(rho)
        #print(mass_per_area)
        mass_per_area = thickness
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
        assert len(props) > 0, all_props
        return props

    def total_thickness(self) -> np.ndarray:
        is_all_nan = np.all(np.isnan(self.thickness))
        if is_all_nan:
            thickness = nonlinear_thickness(self.property_id, self.allowed_properties)
        else:
            thickness = self.thickness
        return thickness


def nonlinear_thickness(property_id: np.ndarray,
                        allowed_properties: list[Any]) -> np.ndarray:
    thickness = np.full(len(property_id), np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties
    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue
        thicknessi = prop.thickness
        thickness[ilookup] = thicknessi[iall]
    return thickness


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

    def get_edge_axes(self):
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

    def _dxyz_centroid_normal_xyz1_xyz2(self, ndim: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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


def shell_thickness(model: BDF,
                    tflag: np.ndarray,
                    T: np.ndarray,
                    property_id: np.ndarray,
                    allowed_properties: list[Any]) -> np.ndarray:
    log = model.log
    thickness = np.full(len(property_id), np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties
    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue
        ti = prop.total_thickness()
        ti_all = ti[iall]

        # set the thickness, even if it's nan
        thickness[ilookup] = ti_all
        if prop.type != 'PSHELL':
            continue

        inan = np.isnan(ti_all)
        if not np.any(inan):
            continue

        #print('inan', inan)
        tflag_nan = tflag[ilookup]
        t_nan = T[ilookup, :]
        #print('tflag_nan', tflag_nan)
        #print('t_nan', t_nan)

        i0 = (tflag_nan == 0)
        i1 = ~i0
        if i0.sum():
            mean_thickness = t_nan[i0, :].mean(axis=1)
            ti_all[i0] = mean_thickness
        if i1.sum():
            scale = t_nan[i1, :].mean(axis=1)
            ti_all[i1] *= scale
            raise RuntimeError('tflag=1')
        #log.error(ti_all)
        # tflag=0: Thickness of element at grid points G1 through G4
        # TFLAG=1:Tthickness becomes a product of Ti and the thickness
        # on the PSHELL card. Ti is ignored for hyperelastic elements.
        # See Remark 6. (Real > 0.0 or blank. See Remark 4 for the default.)
        inan = np.isnan(ti_all)
        if np.any(inan):
            msg = (
                f'tflag={tflag_nan[inan]}\n'
                f'T={ti_all[inan, :]}')
            log.error(msg)
            raise RuntimeError(msg)
        thickness[ilookup] = ti_all

    inan = np.isnan(thickness)
    if np.any(inan):
        msg = f'Thickness has nan\nt={thickness}'
        log.error(msg)
        raise RuntimeError(msg)
    return thickness

def shell_total_thickness(property_id: np.ndarray,
                          allowed_properties: list[Any]) -> np.ndarray:
    thickness = np.full(len(property_id), np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties
    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue
        thicknessi = prop.total_thickness()
        thickness[ilookup] = thicknessi[iall]
    return thickness

def shell_nonstructural_mass(property_id: np.ndarray,
                             allowed_properties: list[Any]) -> np.ndarray:
    nsm = np.full(len(property_id), np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties
    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue
        nsmi = prop.nsm
        nsm[ilookup] = nsmi[iall]
    return nsm

def shell_mass_per_area(model: BDF,
                        tflag: np.ndarray,
                        T: np.ndarray,
                        property_id: np.ndarray,
                        allowed_properties: list[Any]) -> np.ndarray:
    nelement = len(property_id)
    assert nelement > 0, property_id
    mass_per_area = np.full(nelement, np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties

    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue

        if prop.type in {'PCOMP', 'PLPLANE'}:
            mass_per_areai = prop.mass_per_area()
            mass_per_areai_all = mass_per_areai[iall]
            mass_per_area[ilookup] = mass_per_areai_all
            continue

        assert prop.type == 'PSHELL', prop.type
        nsm, rho, ti = prop.nsm_rho_thickness()
        nsm_all = nsm[iall]
        rho_all = rho[iall]
        ti_all = ti[iall]
        mass_per_areai_all = nsm_all + rho_all * ti_all
        mass_per_area[ilookup] = mass_per_areai_all

        inan = np.isnan(ti_all)
        if not np.any(inan):
            continue

        #print('inan', inan)
        tflag_nan = tflag[ilookup]
        t_nan = T[ilookup, :]
        #print('tflag_nan', tflag_nan)
        #print('t_nan', t_nan)

        i0 = (tflag_nan == 0)
        i1 = ~i0
        if i0.sum():
            mean_thickness = t_nan[i0, :].mean(axis=1)
            ti_all[i0] = mean_thickness
        if i1.sum():
            scale = t_nan[i1, :].mean(axis=1)
            ti_all[i1] *= scale
            raise RuntimeError('tflag=1')
        #log.error(ti_all)
        # tflag=0: Thickness of element at grid points G1 through G4
        # TFLAG=1:Tthickness becomes a product of Ti and the thickness
        # on the PSHELL card. Ti is ignored for hyperelastic elements.
        # See Remark 6. (Real > 0.0 or blank. See Remark 4 for the default.)
        inan = np.isnan(ti_all)
        if np.any(inan):
            msg = (
                f'tflag={tflag_nan[inan]}\n'
                f'T={ti_all[inan, :]}')
            model.log.error(msg)
            raise RuntimeError(msg)

        mass_per_area_all = nsm_all + rho_all * ti_all
        #nsm_all = nsm[iall]
        #rho_all = rho[iall]
        #ti_all = ti[iall]
        mass_per_area[ilookup] = mass_per_area_all

        inan = np.isnan(ti_all)
        if np.any(inan):
            msg = f'Thickness has nan\nt={ti_all}'
            model.log.error(msg)
            raise RuntimeError(msg)
    assert nelement > 0, nelement
    assert len(mass_per_area) == nelement, mass_per_area
    return mass_per_area


def shell_mass_per_area_breakdown(model: BDF,
                                  tflag: np.ndarray,
                                  T: np.ndarray,
                                  property_id: np.ndarray,
                                  allowed_properties: list[Any]) -> np.ndarray:
    """
    PCOMP:    [nsm, nan, nan, mass_per_area]
    PLPLANE:  [nan, nan, nan, mass_per_area]
    PSHELL:   [nsm, rho, t,   mass_per_area]
    """
    nelement = len(property_id)
    assert nelement > 0, property_id
    mass_per_area_breakdown = np.full((nelement, 4), np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties

    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue

        if prop.type in {'PCOMP', 'PLPLANE'}:
            if prop.type == 'PCOMP':
                breakdowni = prop.mass_per_area_breakdown() # nsm, mass_per_area
                assert breakdowni.shape[1] == 2, breakdowni.shape
                mass_per_area_breakdown[ilookup, 0] = breakdowni[iall, 0]
                mass_per_area_breakdown[ilookup, 3] = breakdowni[iall, 1]
            elif prop.type == 'PLPLANE':
                mass_per_areai = prop.mass_per_area()
                mass_per_areai_all = mass_per_areai[iall]
                mass_per_area_breakdown[ilookup, 3] = mass_per_areai_all
            else:  ## pragma: no cover
                raise RuntimeError(prop.type)
            #mass_per_area_breakdown[ilookup, 0] = 0.
            #mass_per_area_breakdown[ilookup, 1] = 0.
            #mass_per_area_breakdown[ilookup, 2] = 0.
            continue

        assert prop.type == 'PSHELL', prop.type
        nsm, rho, ti = prop.nsm_rho_thickness()
        nsm_all = nsm[iall]
        rho_all = rho[iall]
        ti_all = ti[iall]
        mass_per_areai_all = nsm_all + rho_all * ti_all
        mass_per_area_breakdown[ilookup, 0] = mass_per_areai_all
        mass_per_area_breakdown[ilookup, 1] = rho_all
        mass_per_area_breakdown[ilookup, 2] = ti_all
        mass_per_area_breakdown[ilookup, 3] = mass_per_areai_all

        inan = np.isnan(ti_all)
        if not np.any(inan):
            continue

        #print('inan', inan)
        tflag_nan = tflag[ilookup]
        t_nan = T[ilookup, :]
        #print('tflag_nan', tflag_nan)
        #print('t_nan', t_nan)

        i0 = (tflag_nan == 0)
        i1 = ~i0
        if i0.sum():
            mean_thickness = t_nan[i0, :].mean(axis=1)
            ti_all[i0] = mean_thickness
        if i1.sum():
            scale = t_nan[i1, :].mean(axis=1)
            ti_all[i1] *= scale
            raise RuntimeError('tflag=1')
        #log.error(ti_all)
        # tflag=0: Thickness of element at grid points G1 through G4
        # TFLAG=1:Tthickness becomes a product of Ti and the thickness
        # on the PSHELL card. Ti is ignored for hyperelastic elements.
        # See Remark 6. (Real > 0.0 or blank. See Remark 4 for the default.)
        inan = np.isnan(ti_all)
        if np.any(inan):
            msg = (
                f'tflag={tflag_nan[inan]}\n'
                f'T={ti_all[inan, :]}')
            model.log.error(msg)
            raise RuntimeError(msg)

        mass_per_area_all = nsm_all + rho_all * ti_all
        #nsm_all = nsm[iall]
        #rho_all = rho[iall]
        #ti_all = ti[iall]
        mass_per_area_breakdown[ilookup, 0] = nsm_all
        mass_per_area_breakdown[ilookup, 1] = rho_all
        mass_per_area_breakdown[ilookup, 2] = ti_all
        mass_per_area_breakdown[ilookup, 3] = mass_per_area_all

        inan = np.isnan(ti_all)
        if np.any(inan):
            msg = f'Thickness has nan\nt={ti_all}'
            model.log.error(msg)
            raise RuntimeError(msg)
    assert nelement > 0, nelement
    assert len(mass_per_area_breakdown) == nelement, mass_per_area_breakdown
    return mass_per_area_breakdown


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
            comment: str=''):
        self.cards.append(((eid, pid, nids,
                            theta_mcid, zoffset,
                            tflag, T1, T2, T3,
                            comment)))
        self.n += 1

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

    def __apply_slice__(self, element: CTRIA3, i: np.ndarray) -> None:
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
            return ''
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
    def area_centroid_normal(self) -> np.ndarray:
        normal = tri_area_centroid_normal(self.model.grid, self.nodes)
        return normal

    @property
    def base_nodes(self):
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
            comment: str=''):
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

    def __apply_slice__(self, element: CTRIAR, i: np.ndarray) -> None:
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

    def write_file_8(self, bdf_file: TextIOLike,
                   write_card_header: bool=False) -> None:
        self.write_file(bdf_file, size=8, is_double=False,
                        write_card_header=write_card_header)

    def write_file_16(self, bdf_file: TextIOLike,
                      is_double=False, write_card_header: bool=False) -> None:
        self.write_file(bdf_file, size=16, is_double=is_double,
                        write_card_header=write_card_header)

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.element_id) == 0:
            return
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

    def area_centroid_normal(self) -> np.ndarray:
        normal = tri_area_centroid_normal(self.model.grid, self.nodes)
        return normal

    @property
    def base_nodes(self):
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
            comment: str=''):
        self.cards.append((eid, pid, nids,
            theta_mcid, zoffset,
            tflag, T1, T2, T3, T4,
            comment))
        self.n += 1

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

    def __apply_slice__(self, element: CQUAD4, i: np.ndarray) -> None:
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
    def area_centroid_normal(self) -> np.ndarray:
        normal = quad_area_centroid_normal(self.model.grid, self.nodes)
        return normal

    @property
    def base_nodes(self):
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

    def __apply_slice__(self, element: CQUADR, i: np.ndarray) -> None:
        element.element_id = self.element_id[i]
        element.property_id = self.property_id[i]
        element.nodes = self.nodes[i, :]
        element.tflag = self.tflag[i]
        element.mcid = self.mcid[i]
        element.theta = self.theta[i]
        element.zoffset = self.zoffset[i]
        element.T = self.T[i, :]
        element.n = len(self.element_id)

    #def add(self, eid: int, pid: int, nids: list[int],
            #theta_mcid: int|float=0.0, zoffset: float=0.,
            #tflag: int=0, T1=None, T2=None, T3=None, T4=None,
            #comment: str=''):
    def add(self, eid: int, pid: int, nids: list[int],
            theta_mcid: int|float=0.0, zoffset: float=0., tflag: int=0,
            T1=None, T2=None, T3=None, T4=None, comment: str='') -> CQUADR:
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

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> str:
        if len(self.element_id) == 0:
            return
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
    def area_centroid_normal(self) -> np.ndarray:
        normal = quad_area_centroid_normal(self.model.grid, self.nodes)
        return normal

    @property
    def base_nodes(self):
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
            comment: str='') -> CTRIA6:
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

    def add_card(self, card: BDFCard, comment: str=''):
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

    def __apply_slice__(self, element: CTRIA6, i: np.ndarray) -> None:
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

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.element_id) == 0:
            return
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

    def area_centroid_normal(self) -> np.ndarray:
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
            comment: str=''):
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

    def add_card(self, card: BDFCard, comment: str=''):
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

    def __apply_slice__(self, element: CQUAD8, i: np.ndarray) -> None:
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

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.element_id) == 0:
            return
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
            theta_mcid: int|float=0., comment: str='') -> CQUAD:
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

    def __apply_slice__(self, element: CQUAD, i: np.ndarray) -> None:
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

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
              write_card_header: bool=False) -> str:
        if len(self.element_id) == 0:
            return

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


class CompositeProperty(Property):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.material_id = np.array([], dtype='int32')
        self.thickness = np.array([], dtype='float64')
        self.theta = np.array([], dtype='float64')
        self.sout = np.array([], dtype='|U8')

        self.nsm = np.array([], dtype='float64')
        self.shear_bonding = np.array([], dtype='float64')
        self.failure_theory = np.array([], dtype='|U8')

        self.tref = np.array([], dtype='float64')
        self.lam = np.array([], dtype='|U3')
        self.z0 = np.array([], dtype='float64')

        self.nlayer = np.array([], dtype='int32')
        self.ge = np.array([], dtype='float64')
        self.group = np.array([], dtype='|U8')

    def sort(self):
        i = np.argsort(self.property_id)
        self.__apply_slice__(self, i)

    @abstractmethod
    def __apply_slice__(self, prop: Any, i: np.ndarray) -> None:  # pragma: no cover
        raise NotImplementedError(f'{self.type}: __apply_slice__')
    @abstractmethod
    def slice_card_by_property_id(self, property_id: np.ndarray):  # pragma: no cover
        raise NotImplementedError(f'{self.type}: slice_card_by_property_id')

    @property
    def ilayer(self) -> np.ndarray:
        return make_idim(self.n, self.nlayer)

    def geom_check(self, missing: dict[str, np.ndarray]):
        mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          msg=f'no shell materials for {self.type}')
        mids.sort()
        geom_check(self,
                   missing,
                   material_id=(mids, self.material_id))

    @property
    def allowed_materials(self) -> list[Any]:
        return [mat for mat in shell_materials(self.model) if mat.n]

    @property
    def symmetry_scale_factor(self):
        scale = np.ones(len(self.property_id), dtype='float64')
        isym = np.where(self.lam == 'SYM')
        scale[isym] = 2.
        return scale

    def expanded_mass_material_id(self) -> tuple[np.ndarray, np.ndarray]:
        property_ids = []
        material_ids = []

        idtype = self.property_id.dtype
        for pid, nsm, lam, z0, nlayer, ilayer in zip_longest(
                self.property_id, self.nsm, self.lam,
                self.z0, self.nlayer, self.ilayer):
            if lam not in {''}:
                msg = f'{self.type} property_id={pid:d} lam={lam!r} is not supported'
                raise NotImplementedError(msg)

            ilayer0, ilayer1 = ilayer
            property_idi = np.full(nlayer, pid, dtype=idtype)
            material_idsi = self.material_id[ilayer0 : ilayer1]
            property_ids.append(property_idi)
            material_ids.append(material_idsi)

            #thicknesses = self.thickness[ilayer0 : ilayer1]
            #_total_thickness = sum(thicknesses)

        property_id = np.hstack(property_ids)
        material_id = np.hstack(material_ids)
        return property_id, material_id

    def expanded_mass_material_id_lookup(self,
                                         element_ids_to_lookup: np.ndarray,
                                         property_ids_to_lookup: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        element_ids = []
        property_ids = []
        material_ids = []

        upids = np.unique(property_ids_to_lookup)
        idtype = self.property_id.dtype
        for pid, nsm, lam, z0, nlayer, ilayer in zip_longest(
                self.property_id, self.nsm, self.lam,
                self.z0, self.nlayer, self.ilayer):
            if lam not in {''}:
                msg = f'{self.type} property_id={pid:d} lam={lam!r} is not supported'
                raise NotImplementedError(msg)

            if pid not in upids:
                continue
            ipid = np.where(pid == property_ids_to_lookup)[0]
            npid = len(ipid)
            ilayer0, ilayer1 = ilayer
            element_idi = element_ids_to_lookup[ipid].reshape(npid, 1)
            #element_idi2 = np.hstack([element_idi]*nlayer).ravel()
            element_idi2 = np.tile(element_idi, nlayer).ravel()

            material_idsi = self.material_id[ilayer0 : ilayer1]
            #material_idsi2 = np.vstack([material_idsi]*npid).ravel()
            material_idsi2 = np.tile(material_idsi, npid)
            #if material_idsi.max() != material_idsi.min():
                #x = 1

            property_idi = np.full(nlayer*npid, pid, dtype=idtype)
            element_ids.append(element_idi2)
            property_ids.append(property_idi)
            material_ids.append(material_idsi2)

            #thicknesses = self.thickness[ilayer0 : ilayer1]
            #_total_thickness = sum(thicknesses)

        element_id = np.hstack(element_ids)
        property_id = np.hstack(property_ids)
        material_id = np.hstack(material_ids)
        return element_id, property_id, material_id

    def mass_per_area(self) -> np.ndarray:
        """"""
        nproperties = len(self.property_id)
        thickness = self.thickness
        symmetry_scale = self.symmetry_scale_factor
        nsm = self.nsm
        ilayers = self.ilayer

        #def _check_pid(pid):
            #if self.type != 'PCOMP':
                #return
            #ipid = np.where(self.property_id == pid)[0]
            #if len(ipid) == 0:
                #return
            #ipid = ipid.item()
            #scale = symmetry_scale[ipid]
            #ilayer = 1
            #idim0, idim1 = ilayers[ipid, :]

            #rho = get_density_from_material(self.material_id, self.allowed_materials, debug=True)
            #rho_t = rho * thickness
            #rhoi = rho[idim0:idim1]
            #thicknessi = thickness[idim0:idim1]
            #mass_per_area = np.zeros(nproperties, dtype='float64') + nsm

            #PCOMP          7
                           #1      .1      0.               1      .2      0.
                           #1      .3      0.               1      .4      0.
                           #1      .5      0.
            #elem.thicknesses
            #[0.1, 0.2, 0.3, 0.4, 0.5]

            #mass_per_area_pid = rho_t[idim0:idim1].sum()
            #pid = self.property_id[ipid]
            #mids_pid = self.material_id[idim0:idim1]
            #x = 1
        #_check_pid(2)
        #_check_pid(7)


        rho = get_density_from_material(self.material_id, self.allowed_materials, debug=False)
        rho_t = rho * thickness
        mass_per_area = np.zeros(nproperties, dtype='float64') + nsm
        for i, scale, ilayer in zip(count(), symmetry_scale, ilayers):
            idim0, idim1 = ilayer
            mass_per_area[i] += scale * rho_t[idim0:idim1].sum()
            #pid = self.property_id[i]
            #if pid == 2:
                #x = 1
        assert len(mass_per_area) == nproperties

        #pid_to_check = [2]
        #ipid2 = np.searchsorted(self.property_id, pid_to_check)
        #mpa = mass_per_area[ipid2]
        #pcard = self.slice_card_by_property_id(pid_to_check)
        #mids = pcard.material_id
        #assert np.allclose(mpa, 0.220467), f'mass/area={mpa}\n{pcard.write()}\nmids={mids}'
        return mass_per_area

    def mass_per_area_breakdown(self) -> np.ndarray:
        """[nsm, mass_per_area]"""
        nproperties = len(self.property_id)
        thickness = self.thickness
        symmetry_scale = self.symmetry_scale_factor
        nsm = self.nsm
        ilayers = self.ilayer

        rho = get_density_from_material(self.material_id, self.allowed_materials, debug=False)
        rho_t = rho * thickness
        mass_per_area = nsm.copy()

        for i, scale, ilayer in zip(count(), symmetry_scale, ilayers):
            idim0, idim1 = ilayer
            mass_per_area[i] += scale * rho_t[idim0:idim1].sum()
        assert len(mass_per_area) == nproperties
        mass_per_area_breakdown = np.column_stack([nsm, mass_per_area])
        return mass_per_area_breakdown

    def total_thickness(self) -> np.ndarray:
        nproperties = len(self.property_id)
        thickness = self.thickness
        symmetry_scale = self.symmetry_scale_factor
        ilayers = self.ilayer

        total_thickness = np.zeros(nproperties, dtype='float64')
        for i, scale, ilayer in zip(count(), symmetry_scale, ilayers):
            idim0, idim1 = ilayer
            total_thickness[i] += scale * thickness[idim0:idim1].sum()

        assert len(total_thickness) == nproperties
        return total_thickness

class PCOMP(CompositeProperty):
    """
    +-------+--------+--------+---------+-------+--------+--------+--------+-------+
    |   1   |    2   |    3   |    4    |   5   |   6    |   7    |    8   |   9   |
    +=======+========+========+=========+=======+========+========+========+=======+
    | PCOMP |   PID  |   Z0   |   NSM   |   SB  |   FT   |  TREF  |   GE   |  LAM  |
    +-------+--------+--------+---------+-------+--------+--------+--------+-------+
    |       |  MID1  |   T1   |  THETA1 | SOUT1 |  MID2  |  T2    | THETA2 | SOUT2 |
    +-------+--------+--------+---------+-------+--------+--------+--------+-------+
    |       |  MID3  |   T3   |  THETA3 | SOUT3 |  etc.  |        |        |       |
    +-------+--------+--------+---------+-------+--------+--------+--------+-------+

    +-------+--------+--------+---------+-------+--------+--------+--------+-------+
    | PCOMP | 701512 | 0.0+0  | 1.549-2 |       |        | 0.0+0  | 0.0+0  |  SYM  |
    +-------+--------+--------+---------+-------+--------+--------+--------+-------+
    |       | 300704 | 3.7-2  |  0.0+0  |  YES  | 300704 | 3.7-2  |   45.  |  YES  |
    +-------+--------+--------+---------+-------+--------+--------+--------+-------+
    |       | 300704 | 3.7-2  |  -45.   |  YES  | 300704 | 3.7-2  |   90.  |  YES  |
    +-------+--------+--------+---------+-------+--------+--------+--------+-------+
    |       | 300705 |   .5   |  0.0+0  |  YES  |        |        |        |       |
    +-------+--------+--------+---------+-------+--------+--------+--------+-------+
    """

    def add(self, pid: int, mids: list[int], thicknesses: list[float],
            thetas: Optional[list[float]]=None, souts: Optional[list[str]]=None,
            nsm: float=0., sb: float=0., ft: str='',
            tref: float=0., ge: float=0., lam: str='',
            z0=None, comment='') -> PCOMP:
        """
        Creates a PCOMP card

        Parameters
        ----------
        pid : int
            property id
        mids : list[int, ..., int]
            material ids for each ply
        thicknesses : list[float, ..., float]
            thicknesses for each ply
        thetas : list[float, ..., float]; default=None
            ply angle
            None : [0.] * nplies
        souts : list[str, ..., str]; default=None
            should the stress? be printed; {YES, NO}
            None : [NO] * nplies
        nsm : float; default=0.
            nonstructural mass per unit area
        sb : float; default=0.
            Allowable shear stress of the bonding material.
            Used by the failure theory
        ft : str; default='
            failure theory; {HILL, HOFF, TSAI, STRN, None}
        tref : float; default=0.
            reference temperature
        ge : float; default=0.
            structural damping
        lam : str; default=''
            symmetric flag; {SYM, MEM, BEND, SMEAR, SMCORE, ''}
            None : not symmmetric
        z0 : float; default=None
            Distance from the reference plane to the bottom surface
            None : -1/2 * total_thickness
        comment : str; default=''
            a comment for the card

        """
        if ft is None:
            ft = ''
        if lam is None:
            lam = ''

        if thetas is None:
            thetas = [0.] * len(mids)
        if souts is None:
            souts = ['YES'] * len(mids)

        assert len(souts) == len(mids)
        assert len(thetas) == len(mids)
        self.cards.append((pid, nsm, sb, ft, tref, ge, lam, z0,
                           mids, thicknesses, thetas, souts, comment))
        #self.property_id = np.hstack([self.property_id, pid])
        #self.lam = np.hstack([self.lam, lam])
        #if z0 is None:
            #total_thickness = sum(thicknesses)
            #if lam == 'SYM':
                #total_thickness *= 2
            #z0 = -0.5 * total_thickness

        #self.z0 = np.hstack([self.z0, z0])
        #self.ge = np.hstack([self.ge, ge])
        #self.shear_bonding = np.hstack([self.shear_bonding, sb])
        #self.nsm = np.hstack([self.nsm, nsm])
        #self.tref = np.hstack([self.tref, tref])
        #self.failure_theory = np.hstack([self.failure_theory, ft])

        #self.material_id = np.hstack([self.material_id, mids])
        #self.thickness = np.hstack([self.thickness, thicknesses])
        #self.theta = np.hstack([self.theta, thetas])
        #self.sout = np.hstack([self.sout, souts])
        #self.nlayer = np.hstack([self.nlayer, len(thicknesses)])
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        #print(card.write_card(size=16))

        pid = integer(card, 1, 'pid')

        # z0 is field 2 and is calculated at the end because we need the
        # thickness first
        #self.z0 = double_or_blank(card, 1, 'pid')

        nsm = double_or_blank(card, 3, 'nsm', default=0.0)
        shear_bonding = double_or_blank(card, 4, 'sb', default=0.0)
        ft = string_or_blank(card, 5, 'ft', default='')
        tref = double_or_blank(card, 6, 'tref', default=0.0)
        ge = double_or_blank(card, 7, 'ge', default=0.0)
        lam = string_or_blank(card, 8, 'lam', default='') # default=blank -> nothing
        assert len(lam) <= 6, lam  # MEM, SMEAR, SMCORE

        # -8 for the first 8 fields (1st line)
        nply_fields = card.nfields - 9

        # counting plies
        nmajor = nply_fields // 4
        nleftover = nply_fields % 4
        if nleftover:
            nmajor += 1
        nplies = nmajor

        mid_last = None
        thick_last = None
        #ply = None
        iply = 1

        # supports single ply per line
        mids = []
        thicknesses = []
        thetas = []
        souts = []
        for ioffset in range(9, 9 + nplies * 4, 4):
            actual = card.fields(ioffset, ioffset + 4)
            mid = integer_or_blank(card, ioffset, 'mid', default=mid_last)
            t = double_or_blank(card, ioffset + 1, 't', default=thick_last)
            theta = double_or_blank(card, ioffset + 2, 'theta', default=0.0)
            sout = string_or_blank(card, ioffset + 3, 'sout', default='NO')

            # if this card has 2 plies on the line
            if actual != [None, None, None, None]:
                mids.append(mid)
                thicknesses.append(t)
                thetas.append(theta)
                souts.append(sout)
                iply += 1
            mid_last = mid
            thick_last = t
        #print("nplies = %s" % nplies)

        #self.plies = []
        #if self.lam == 'SYM':
        #    if nplies % 2 == 1:  # 0th layer is the core layer
        #       # cut the thickness in half to make the ply have an
        #       # even number of plies, is there a better way???
        #       plies[0][1] = plies[0][1] / 2.
        #
        #    plies_lower = plies.reverse()
        #    self.plies = plies_lower + plies
        #    #print str(self)
        z0 = double_or_blank(card, 2, 'z0')

        self.cards.append((pid, nsm, shear_bonding, ft, tref, ge, lam, z0,
                           mids, thicknesses, thetas, souts, comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return

        property_id = np.zeros(ncards, dtype='int32')
        z0 = np.zeros(ncards, dtype='float64')
        nsm = np.zeros(ncards, dtype='float64')
        shear_bonding = np.zeros(ncards, dtype='float64')

        # 'HILL' for the Hill theory.
        # 'HOFF' for the Hoffman theory.
        # 'TSAI' for the Tsai-Wu theory.
        # 'STRN' for the Maximum Strain theory.
        failure_theory = np.zeros(ncards, dtype='|U8')
        tref = np.zeros(ncards, dtype='float64')
        ge = np.zeros(ncards, dtype='float64')
        lam = np.zeros(ncards, dtype='|U6')

        nlayer = np.zeros(ncards, dtype='int32')
        #group = np.full(ncards, '', dtype='|U8')

        mids_list = []
        thickness_list = []
        thetas_list = []
        sout_list = []
        for icard, card in enumerate(self.cards):
            (pid, nsmi, shear_bondingi, fti, trefi, gei, lami, z0i,
             mids, thicknesses, thetas, souts, comment) = card

            nsm[icard] = nsmi
            shear_bonding[icard] = shear_bondingi
            failure_theory[icard] = fti
            tref[icard] = trefi
            ge[icard] = gei
            lam[icard] = lami
            assert lami in {'', 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE'}, f'pid={pid} laminate={lami!r}'
            assert fti in {'', 'HILL', 'HOFF', 'STRN', 'TSAI', 'HFAI', 'HFAB', 'HTAP'}, f'pid={pid} failure_theory={fti!r}'
            nlayersi = len(mids)

            property_id[icard] = pid
            mids_list.extend(mids)
            thickness_list.extend(thicknesses)
            sout_list.extend(souts)
            thetas_list.extend(thetas)
            z0[icard] = z0i
            nlayer[icard] = nlayersi
            nsm[icard] = nsmi

        material_id = np.array(mids_list, dtype='int32')
        thickness = np.array(thickness_list, dtype='float64')
        sout = np.array(sout_list, dtype='|U4') # YES, NO, YESA
        theta = np.array(thetas_list, dtype='float64')
        assert len(mids_list) == nlayer.sum()
        self._save(property_id,
                   nlayer, material_id, thickness, sout, theta,
                   z0, nsm, shear_bonding, failure_theory, tref, ge, lam)
        self.sort()
        self.cards = []

    def _save(self, property_id,
              nlayer, material_id, thickness, sout, theta,
              z0, nsm, shear_bonding, failure_theory, tref, ge, lam):
        if len(self.property_id):
            property_id = np.hstack([self.property_id, property_id])

            nlayer = np.hstack([self.nlayer, nlayer])
            material_id = np.hstack([self.material_id, material_id])
            thickness = np.hstack([self.thickness, thickness])
            sout = np.hstack([self.sout, sout])
            theta = np.hstack([self.theta, theta])

            z0 = np.hstack([self.z0, z0])
            nsm = np.hstack([self.nsm, nsm])
            shear_bonding = np.hstack([self.shear_bonding, shear_bonding])
            failure_theory = np.hstack([self.failure_theory, failure_theory])
            tref = np.hstack([self.tref, tref])
            ge = np.hstack([self.ge, ge])
            lam = np.hstack([self.lam, lam])

        self.property_id = property_id

        self.nlayer = nlayer
        self.material_id = material_id
        self.thickness = thickness
        self.sout = sout
        self.theta = theta

        self.z0 = z0
        self.nsm = nsm
        self.shear_bonding = shear_bonding
        self.failure_theory = failure_theory
        self.tref = tref
        self.ge = ge
        self.lam = lam
        assert len(self.nlayer) == self.n

        inan = np.isnan(z0)
        if inan.sum():
            total_thickness = self.total_thickness()
            self.z0[inan] = -total_thickness[inan] / 2
        #x = 1

    def add_op2_data(self, data, comment=''):
        """
        Adds a PCOMP card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        #data_in = [
            #pid, z0, nsm, sb, ft, tref, ge,
            #is_symmetrical, Mid, T, Theta, Sout]
        pid = data[0]
        z0 = data[1]
        nsm = data[2]
        sb = data[3]
        ft_int = data[4]
        tref = data[5]
        ge = data[6]
        lam = data[7]
        Mid = data[8]
        T = data[9]
        Theta = data[10]
        Sout = data[11]

        if lam == 'NO':
            lam = None

        mids = []
        thicknesses = []
        thetas = []
        souts = []
        for (mid, t, theta, sout) in zip_longest(Mid, T, Theta, Sout):
            if sout == 0:
                sout = 'NO'
            elif sout == 1:
                sout = 'YES'
            #elif sout == 2:  #: .. todo:: what?!!
                #sout = 'YES'
            #elif sout == 3:  #: .. todo:: what?!!
                #sout = 'YES'
            else:
                raise RuntimeError(f'unsupported sout.  sout={sout!r} and must be 0 or 1.'
                                   f'\nPCOMP = {data}')
            mids.append(mid)
            thicknesses.append(t)
            thetas.append(theta)
            souts.append(sout)
            try:
                ft = map_failure_theory_int(ft_int)
            except NotImplementedError:
                raise RuntimeError(f'unsupported ft.  pid={pid} ft={ft_int!r}.'
                               f'\nPCOMP = {data}')
        self.add(pid, mids, thicknesses, thetas=thetas, souts=souts,
                 nsm=nsm, sb=sb, ft=ft, tref=tref, ge=ge, lam=lam, z0=z0)
        #return PCOMP(pid, mids, thicknesses, thetas, souts,
                     #nsm, sb, ft, tref, ge, lam, z0, validate=False, comment=comment)

    @property
    def nplies(self) -> np.ndarray:
        return self.nlayer
    @nplies.setter
    def nplies(self, nplies: np.ndarray) -> np.ndarray:
        self.nlayer = nplies

    @property
    def nplies_total(self) -> np.ndarray:
        return self.nlayer * (self.is_symmetrical + 1)
    #@nplies.setter
    #def nplies_total(self, nplies: np.ndarray) -> np.ndarray:
        #self.nlayer = nplies

    @property
    def is_symmetrical(self) -> np.ndarray:
        """
        Is the laminate symmetrical?

        Returns
        -------
        is_symmetrical : bool
            is the SYM flag active?

        """
        is_symmetric = (self.lam == 'SYM')
        return is_symmetric

    def slice_card_by_property_id(self, property_id: np.ndarray) -> PCOMP:
        """uses a node_ids to extract GRIDs"""
        iprop = self.index(property_id)
        #assert len(self.node_id) > 0, self.node_id
        #i = np.searchsorted(self.node_id, node_id)
        prop = self.slice_card_by_index(iprop)
        return prop

    def __apply_slice__(self, prop: PCOMP, i: np.ndarray) -> None:
        assert self.nlayer.sum() == len(self.thickness)
        prop.property_id = self.property_id[i]
        prop.z0 = self.z0[i]
        prop.nsm = self.nsm[i]
        prop.shear_bonding = self.shear_bonding[i]

        # 'HILL' for the Hill theory.
        # 'HOFF' for the Hoffman theory.
        # 'TSAI' for the Tsai-Wu theory.
        # 'STRN' for the Maximum Strain theory.
        prop.failure_theory = self.failure_theory[i]
        prop.tref = self.tref[i]
        prop.ge = self.ge[i]
        prop.lam = self.lam[i]
        #prop.group = self.group[i]

        ilayer = self.ilayer # [i, :]
        prop.material_id = hslice_by_idim(i, ilayer, self.material_id)
        prop.sout = hslice_by_idim(i, ilayer, self.sout)
        prop.theta = hslice_by_idim(i, ilayer, self.theta)
        prop.thickness = hslice_by_idim(i, ilayer, self.thickness)

        prop.nlayer = self.nlayer[i]
        assert prop.nlayer.sum() == len(prop.thickness), f'prop.nlayer={prop.nlayer} len(prop.thickness)={len(prop.thickness)}'
        prop.n = len(i)

    @property
    def sb(self):
        return self.shear_bonding
    @sb.setter
    def sb(self, sb):
        self.shear_bonding = sb

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.property_id) == 0:
            return
        print_card = get_print_card_8_16(size)
        assert self.failure_theory.dtype.name == 'str256', self.failure_theory.dtype.name
        for pid, nsm, shear_bonding, failure_theory, ge, tref, lam, z0, \
            nlayer, ilayer in zip_longest(self.property_id, self.nsm, self.shear_bonding,
                                          self.failure_theory, self.tref, self.ge, self.lam, self.z0,
                                          self.nlayer, self.ilayer):

            nsm2 = set_blank_if_default(nsm, 0.0)
            sb2 = set_blank_if_default(shear_bonding, 0.0)
            tref2 = set_blank_if_default(tref, 0.0)
            ge2 = set_blank_if_default(ge, 0.0)

            ilayer0, ilayer1 = ilayer
            material_ids = self.material_id[ilayer0 : ilayer1].tolist()
            thetas = self.theta[ilayer0 : ilayer1].tolist()
            thicknesses = self.thickness[ilayer0 : ilayer1].tolist()
            souts = self.sout[ilayer0 : ilayer1].tolist()

            _total_thickness = sum(thicknesses)
            z02 = set_blank_if_default(z0, -0.5 * _total_thickness)

            #ndim = self.valid_types[beam_type]
            #assert len(dim) == ndim, 'PCOMP ndim=%s len(dims)=%s' % (ndim, len(dim))
            list_fields = ['PCOMP', pid, z02, nsm2, sb2, failure_theory, tref2, ge2, lam]
            for (mid, t, theta, sout) in zip(material_ids, thicknesses,
                                             thetas, souts):
                #assert sout not in ['0', 0.0, 0]
                assert sout in {'YES', 'NO'}, sout
                #theta = set_blank_if_default(theta, 0.0)
                str_sout = set_blank_if_default(sout, 'NO')
                list_fields += [mid, t, theta, str_sout]

            bdf_file.write(print_card(list_fields))
        return

    def get_z_locations(self):
        nproperties = len(self.property_id)
        nplies = self.nplies_total
        nlayers_max = nplies.max() + 1
        shape = (nproperties, nlayers_max)

        z_locations = np.full(shape, np.nan, dtype='float64')
        for iprop, lam, z0, (ilayer0, ilayer1) in zip(count(), self.lam, self.z0, self.ilayer):
            if lam == '':
                thicknesses = self.thickness[ilayer0:ilayer1]
                csum = np.cumsum(thicknesses)
                nlayers = len(thicknesses)

                z0 = z0 + np.hstack([0., csum[:-1]])
                z1 = z0 + thicknesses
                z_locations[iprop, :nlayers] = z0
                z_locations[iprop, nlayers] = z1[-1]
            else:
                raise NotImplementedError(lam)
        return z_locations

    def get_individual_ABD_matrices(
            self, theta_offset: float=0.) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Gets the ABD matrix

        Parameters
        ----------
        theta_offset : float
            rotates the ABD matrix; measured in degrees

        http://www2.me.rochester.edu/courses/ME204/nx_help/index.html#uid:id503291
        Understanding Classical Lamination Theory
        """
        #mids = self.get_material_ids()
        thickness = self.total_thickness()

        #z0 = self.z1
        #z1 = self.z2
        #zmean = (z0 + z1) / 2.
        #dz = z1 - z0
        #dzsquared = z1 ** 2 - z0 ** 2
        #zcubed = z1 ** 3 - z0 ** 3
        # A11 A12 A16
        # A12 A22 A26
        # A16 A26 A66
        nprop = len(self.property_id)
        A = np.zeros((nprop, 3, 3), dtype='float64')
        B = np.zeros((nprop, 3, 3), dtype='float64')
        D = np.zeros((nprop, 3, 3), dtype='float64') # TODO: 2x2 matrix?

        #theta = np.zeros(nprop, dtype='float64')
        #z0 = self.z0
        #z1 = self.z[:, 1]
        #mean = (z0 + z1) / 2.

        #thicknesses = self.get_thicknesses()

        #if self.mid1_ref:

        #if self.is_symmetrical:
            #mids_ref = copy.deepcopy(self.mids_ref)
            #mids_ref += mids_ref[::-1]

        #assert len(mids) == len(mids_ref), f'mids={mids} ({len(mids)}) mids_ref:\n{mids_ref}; {len(mids_ref)}'
        for lam, z0, (ilayer0, ilayer1) in zip(self.lam, self.z0, self.ilayer):
            mids = self.material_id[ilayer0:ilayer1]
            thetas = self.theta[ilayer0:ilayer1]
            thicknesses = self.thickness[ilayer0:ilayer1]
            csum = np.cumsum(thicknesses)
            z0 = z0 + np.hstack([0., csum[:-1]])
            z1 = z0 + csum
            zmeans = (z0 + z1) / 2

            Qbar = self.get_Qbar_matrix(mids, thetas)
            for Qbari, thickness, zmean, z0i, z1i in zip(Qbar, thicknesses, zmeans, z0, z1):
                A += Qbari * thickness
                #B += Qbar * thickness * zmean
                #D += Qbar * thickness * (z1i ** 3 - z0i ** 3)
                B += Qbari * (z1i ** 2 - z0i ** 2)
                D += Qbari * (z1i ** 3 - z0i ** 3)

            # [N, M, Q].T =   [TG1, T^2 * G4, 0]              * [epsilon0]
                            # [T^2 * G4, T^3/12 * G2, 0]        [xi]
                            # [0, 0, Ts * G3]                   [gamma]

            #B += Qbar * thickness * zmean
            #D += Qbar * thickness * (z1i ** 3 - z0i ** 3)
            #N += Qbar * alpha * thickness
            #M += Qbar * alpha * thickness * zmean
            B /= 2.
            D /= 3.
        #M /= 2.
        #np.set_printoptions(linewidth=120, suppress=True)
        #print(ABD)
        return A, B, D

    def get_Qbar_matrix(self, mid, theta: float=0.):
        """theta must be in radians"""
        S2 = get_mat_s33(self.model, mid, self.allowed_materials)
        T = get_2d_plate_transform(theta)
        #Tt = np.transpose(T, axes=0)
        #Tinv = np.linalg.inv(T)
        #Qbar = np.linalg.multi_dot([Tinv, Q, Tinv.T])
        #ST = np.einsum('ijk,ilm->ijm', S2, T)
        #Sbarv = np.einsum('ijk,ilm->ikm', T, ST)
        Sbar = [np.linalg.multi_dot([Ti.T, Si, Ti]) for Ti, Si in zip(T, S2)]
        Qbar = [np.linalg.inv(Sbari) for Sbari in Sbar]
        #Ti = T[0, :, :]
        #Si = S2[0, :, :]
        #Sbar_test = np.linalg.multi_dot([Ti.T, Si, Ti])
        #Qbarv = np.linalg.inv(Sbarv)
        Qbarv = np.array(Qbar, dtype=S2.dtype)
        return Qbarv

    def get_Sbar_matrix(self, mid_ref, theta=0.):
        """theta must be in radians"""
        # this is the inverse of Sbar
        S2, unused_S3 = get_mat_props_S(mid_ref)
        T = get_2d_plate_transform(theta)
        #Tinv = np.linalg.inv(T)
        Sbar = np.linalg.multi_dot([T.T, S2, T])
        return Sbar

    def get_ABD_matrices(self, theta_offset=0.) -> np.ndarray:
        """
        Gets the ABD matrix

        Parameters
        ----------
        theta_offset : float
            rotates the ABD matrix; measured in degrees

        """
        A, B, D = self.get_individual_ABD_matrices(theta_offset=theta_offset)
        ABD = np.block([
            [A, B],
            [B, D],
        ])
        return ABD

    def get_Ainv_equivalent_pshell(self,
                                   imat_rotation_angle: float,
                                   thickness: float) -> tuple[float, float, float, float]:
        """imat_rotation_angle is in degrees...but is specified in radians and unused"""
        ABD = self.get_ABD_matrices(imat_rotation_angle)
        A = ABD[:, :3, :3]
        Ainv = np.linalg.inv(A)

        # equivalent compliance matrix
        S = thickness * Ainv
        #print(S[0, 0], S[1, 1], S[2, 2])

        # these are on the main diagonal
        Ex = 1. / S[0, 0]
        Ey = 1. / S[1, 1]
        Gxy = 1. / S[2, 2]
        #S12 = -nu12 / E1 = -nu21/E2
        #nu12 = -S12 / E1
        nu_xy = -S[0, 1] / Ex
        return Ex, Ey, Gxy, nu_xy


class PCOMPG(CompositeProperty):

    def add(self, pid, global_ply_ids, mids, thicknesses, thetas=None, souts=None,
                   nsm=0.0, sb=0.0, ft='', tref=0.0, ge=0.0, lam=None, z0=None,
                   comment='') -> PCOMPG:
        lam = lam if lam is not None else ''
        z0 = z0 if z0 is not None else np.nan

        nlayers = len(thicknesses)
        if thetas is None:
            thetas = [0.] * nlayers
        if souts is None:
            souts = ['NO'] * nlayers

        cardi = (
            pid, nsm, sb, ft, tref, ge, lam, z0,
            mids, thicknesses, thetas, souts, global_ply_ids,
            comment)
        self.cards.append(cardi)
        self.n += 1

    def __apply_slice__(self, prop: PCOMPG, i: np.ndarray) -> None:
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        prop.z0 = self.z0[i]
        prop.nsm = self.nsm[i]
        prop.shear_bonding = self.shear_bonding[i]

        # 'HILL' for the Hill theory.
        # 'HOFF' for the Hoffman theory.
        # 'TSAI' for the Tsai-Wu theory.
        # 'STRN' for the Maximum Strain theory.
        prop.failure_theory = self.failure_theory[i]
        prop.tref = self.tref[i]
        prop.ge = self.ge[i]
        prop.lam = self.lam[i]
        prop.group = self.group[i]

        ilayer = self.ilayer # [i, :]
        prop.material_id = hslice_by_idim(i, ilayer, self.material_id)
        prop.sout = hslice_by_idim(i, ilayer, self.sout)
        prop.theta = hslice_by_idim(i, ilayer, self.theta)
        prop.thickness = hslice_by_idim(i, ilayer, self.thickness)
        prop.global_ply_id = hslice_by_idim(i, ilayer, self.global_ply_id)
        prop.nlayer = self.nlayer[i]


    def add_card(self, card: BDFCard, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        # z0 will be calculated later
        nsm = double_or_blank(card, 3, 'nsm', default=0.0)
        sb = double_or_blank(card, 4, 'sb', default=0.0)
        ft = string_or_blank(card, 5, 'ft', default='')
        tref = double_or_blank(card, 6, 'tref', default=0.0)
        ge = double_or_blank(card, 7, 'ge', default=0.0)
        lam = string_or_blank(card, 8, 'lam', default='')
        assert len(lam) <= 2, lam

        fields = card.fields(9)

        #T = 0.  # thickness
        mid_last = None
        thick_last = None

        ioffset = 0
        mids = []
        thicknesses = []
        thetas = []
        souts = []
        global_ply_ids = []
        while ioffset < len(fields):
            global_ply_id = integer(card, 9 + ioffset, 'global_ply_id')
            mid = integer_or_blank(card, 9 + ioffset + 1, 'mid', default=mid_last)

            # can be blank 2nd time thru
            thickness = double_or_blank(card, 9 + ioffset + 2, 'thickness', default=thick_last)

            theta = double_or_blank(card, 9 + ioffset + 3, 'theta', default=0.0)
            sout = string_or_blank(card, 9 + ioffset + 4, 'sout', default='NO')
            #print('n=%s global_ply_id=%s mid=%s thickness=%s len=%s' %(
            #    n,global_ply_id,mid,thickness,len(fields)))

            mids.append(mid)
            thicknesses.append(thickness)
            thetas.append(theta)
            souts.append(sout)
            global_ply_ids.append(global_ply_id)

            assert mid is not None
            assert thickness is not None
            #assert isinstance(mid, integer_types), 'mid=%s' % mid
            assert isinstance(thickness, float), 'thickness=%s' % thickness
            mid_last = mid
            thick_last = thickness
            #T += thickness
            ioffset += 8
            #n += 1
        z0 = double_or_blank(card, 2, 'z0', default=np.nan)

        cardi = (
            pid, nsm, sb, ft, tref, ge, lam, z0,
            mids, thicknesses, thetas, souts, global_ply_ids,
            comment)
        self.cards.append(cardi)
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return

        self.property_id = np.zeros(ncards, dtype='int32')
        self.nsm = np.zeros(ncards, dtype='float64')
        self.shear_bonding = np.zeros(ncards, dtype='float64')
        self.failure_theory = np.zeros(ncards, dtype='|U8')

        self.tref = np.zeros(ncards, dtype='float64')
        self.lam = np.zeros(ncards, dtype='|U5')
        self.z0 = np.zeros(ncards, dtype='float64')

        self.nlayer = np.zeros(ncards, dtype='int32')
        self.ge = np.zeros(ncards, dtype='float64')
        self.group = np.full(ncards, '', dtype='|U8')

        global_ply_ids_list = []
        mids_list = []
        thickness_list = []
        thetas_list = []
        sout_list = []
        for icard, card in enumerate(self.cards):
            (pid, nsm, sb, ft, tref, ge, lam, z0,
             mids, thicknesses, thetas, souts, global_ply_ids, commment) = card

            self.failure_theory[icard] = ft
            self.shear_bonding[icard] = sb
            self.tref[icard] = tref
            self.ge[icard] = ge
            self.lam[icard] = lam

            nlayers = len(mids)
            self.property_id[icard] = pid
            global_ply_ids_list.extend(global_ply_ids)
            mids_list.extend(mids)
            thickness_list.extend(thicknesses)
            sout_list.extend(souts)
            thetas_list.extend(thetas)
            self.z0[icard] = z0
            self.nlayer[icard] = nlayers
            self.nsm[icard] = nsm

        self.global_ply_id = np.array(global_ply_ids_list, dtype='int32')
        self.material_id = np.array(mids_list, dtype='int32')
        self.thickness = np.array(thickness_list, dtype='float64')
        self.sout = np.array(sout_list, dtype='|U4') # YES, NO, YESA
        self.theta = np.array(thetas_list, dtype='float64')
        assert len(mids_list) == self.nlayer.sum()
        self.sort()
        self.cards = []

    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.property_id) == 0:
            return
        print_card = get_print_card_8_16(size)
        property_id = array_str(self.property_id, size=size)
        for pid, nsm, sb, failure_theory, ge, tref, lam, z0, \
            ilayer in zip_longest(property_id, self.nsm, self.shear_bonding,
                                  self.failure_theory, self.tref, self.ge, self.lam, self.z0,
                                  self.ilayer):
            nsm2 = set_blank_if_default(nsm, 0.0)
            sb2 = set_blank_if_default(sb, 0.0)
            tref2 = set_blank_if_default(tref, 0.0)
            ge2 = set_blank_if_default(ge, 0.0)

            ilayer0, ilayer1 = ilayer
            global_ids = self.global_ply_id[ilayer0 : ilayer1].tolist()
            material_ids = self.material_id[ilayer0 : ilayer1].tolist()
            thetas = self.theta[ilayer0 : ilayer1].tolist()
            thicknesses = self.thickness[ilayer0 : ilayer1].tolist()
            souts = self.sout[ilayer0 : ilayer1].tolist()

            _total_thickness = sum(thicknesses)
            z02 = set_blank_if_default(z0, -0.5 * _total_thickness)

            #ndim = self.valid_types[beam_type]
            #assert len(dim) == ndim, 'PCOMP ndim=%s len(dims)=%s' % (ndim, len(dim))
            list_fields = ['PCOMPG', pid, z02, nsm2, sb2, failure_theory, tref2, ge2, lam]
            for (global_id, mid, t, theta, sout) in zip(global_ids, material_ids, thicknesses,
                                                        thetas, souts):
                #theta = set_blank_if_default(theta, 0.0)
                str_sout = set_blank_if_default(sout, 'NO')
                list_fields += [global_id, mid, t, theta, str_sout, None, None, None]

            bdf_file.write(print_card(list_fields))
        return

def get_2d_plate_transform(theta: np.ndarray) -> np.ndarray:
    """theta must be in radians"""
    ntheta = len(theta)
    ct = np.cos(theta)
    st = np.sin(theta)
    ct2 = ct ** 2
    st2 = st ** 2
    cst = st * ct

    T126 = np.zeros((ntheta, 3, 3), dtype='float64')
    T126[:, 0, 0] = T126[:, 1, 1] = ct2
    T126[:, 0, 1] = T126[:, 0, 1] = st2
    T126[:, 0, 2] = 2 * cst
    T126[:, 1, 2] = -2 * cst
    T126[:, 2, 0] = -cst
    T126[:, 2, 1] = cst
    T126[:, 2, 2] = ct2 - st2
    #T126 = np.array([
        #[ct2, st2,    2 * cst],
        #[st2, ct2,   -2 * cst],
        #[-cst, cst, ct2 - st2],
    #])
    #T126inv = np.array([
        #[ct2,  st2,   -2 * cst],
        #[st2,  ct2,    2 * cst],
        #[cst, -cst, ct2 - st2],
    #])
    T = T126
    #T123456 = np.array([
        #[     ct2,     st2,   0., 0.,  0.,  -2 * cst],
        #[     st2,     ct2,   0., 0.,  0.,   2 * cst],
        #[      0.,      0.,   1., 0.,  0.,        0.],
        #[      0.,      0.,   0., ct, -st,        0.],
        #[      0.,      0.,   0., st,  ct,        0.],
        #[-ct * st, ct * st,   0., 0.,  0., ct2 - st2],
    #])
    #T123456 = np.array([
        #[     ct2,     st2,  -2 * cst],
        #[     st2,     ct2,   2 * cst],
        #[-ct * st, ct * st, ct2 - st2],
    #])
    return T

def get_mat_props_S(mid: np.ndarray):
    """
    Gets the material matrix [S] or [C] for plane strain

    [e] = [S][o]
    """
    mtype = mid_ref.type
    if mtype == 'MAT1':
        e = mid_ref.e
        g = mid_ref.g
        nu = mid_ref.nu
        # http://web.mit.edu/16.20/homepage/3_Constitutive/Constitutive_files/module_3_with_solutions.pdf
        # [e11, e22, 2*e12] = ei @  [o11, o22, o12]
        #[e] = [S][o]
        #[o] = [C][e]
        # eq 3.35 (2d)
        # eq 3.50 (3d)
        ei2 = np.array([
            [  1 / e, -nu / e,    0.],
            [-nu / e,   1 / e,    0.],
            [     0.,      0., 1 / g],
        ])
        #G = E / (2*(1 + nu))
        #1 / G = (2*(1 + nu)) / E
        nu2 = 2 * (1 + nu)
        ei3 = np.array([
            [1, -nu, -nu, 0., 0., 0.],
            [-nu, 1, -nu, 0., 0., 0.],
            [-nu, -nu, 1, 0., 0., 0.],
            [0., 0., 0., nu2, 0., 0.],
            [0., 0., 0., 0., nu2, 0.],
            [0., 0., 0., 0., 0., nu2],
        ]) / e

        #denom = e / (1 - nu ** 2)
        #C2 = np.array([
            #[1., -nu, 0.],
            #[nu, 1., 0.],
            #[0., 0., g / denom],
        #]) * denom

        #lambd = e * nu / (1 + nu) / (1 - 2 * nu)
        #lambda_2u = lambd + 2 * g
        #C3 = np.array([
            #[lambda_2u, lambd, lambd, 0., 0., 0.],
            #[lambd, lambda_2u, lambd, 0., 0., 0.],
            #[lambd, lambd, lambda_2u, 0., 0., 0.],
            #[0., 0., 0., g, 0., 0.],
            #[0., 0., 0., 0., g, 0.],
            #[0., 0., 0., 0., 0., g],
        #])

    elif mtype == 'MAT8':
        # orthotropic
        material = mid_ref
        # http://web.mit.edu/16.20/homepage/3_Constitutive/Constitutive_files/module_3_with_solutions.pdf
        # [e11, e22, 2*e12] = ei @  [o11, o22, o12]
        # eq 3.35 (2d)
        # eq 3.50 (3d)
        #ei2 = np.array([
            #[e, -nu / e, 0., ],
            #[-nu / e, e, 0., ],
            #[0., 0., 1/g],
        #])
        #G = E / (2*(1 + nu))
        #1 / G = (2*(1 + nu)) / E

        #  https://en.wikipedia.org/wiki/Orthotropic_material
        e1, e2 = material.e11, material.e22 # , material.e33
        e3 = 1.
        nu12 = material.nu12
        g12, g31, g23 = material.g12, material.g1z, material.g2z
        if g12 == 0.:
            g12 = 1.
        if g31 == 0.:
            g31 = 1.
        if g23 == 0.:
            g23 = 1.

        # nu21 * E1 = nu12 * E2
        nu13 = nu12 # assume; should fall out in calcs given e3=0
        nu23 = nu12 # assume; should fall out in calcs given e3=0
        nu21 = nu12 * e2 / e1
        nu31 = nu13 * e3 / e1
        nu32 = nu23 * e3 / e2
        ei2 = np.array([
            [    1/e1, -nu21/e2,    0.],
            [-nu12/e1,     1/e2,    0.],
            [      0.,       0., 1/g12],
        ])
        ei3 = np.array([
            [    1/e1, -nu21/e2, -nu31/e3,    0.,    0.,    0.],
            [-nu12/e1,     1/e2, -nu32/e3,    0.,    0.,    0.],
            [-nu13/e1, -nu23/e2,     1/e3,    0.,    0.,    0.],
            [      0.,       0.,       0., 1/g23,    0.,    0.],
            [      0.,       0.,       0.,    0., 1/g31,    0.],
            [      0.,       0.,       0.,    0.,    0., 1/g12],
        ])
        #denom = 1 - nu12 * nu21
        #C2 = np.array([
            #[e1, -nu21 * e1, 0.],
            #[nu12 * e2, e2, 0.],
            #[0., 0., g12 * denom],
        #]) / denom

    else:
        raise NotImplementedError(mid_ref.get_stats())
    return ei2, ei3

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
