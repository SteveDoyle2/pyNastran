from __future__ import annotations
from abc import abstractmethod
from itertools import count, zip_longest
from typing import Optional, Any, TYPE_CHECKING

import numpy as np
#from pyNastran.bdf.field_writer_8 import print_field_8, print_card_8
#from pyNastran.bdf.field_writer_16 import print_field_16, print_card_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, string, # double,
    integer_or_blank, double_or_blank,
    #integer_double_or_blank,
    string_or_blank, # blank,
    integer_types, float_types,
)
from pyNastran.bdf.cards.elements.bars import set_blank_if_default
from pyNastran.bdf.cards.properties.shell import map_failure_theory_int

from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    Property,
    hslice_by_idim, make_idim, searchsorted_filter,
    parse_check, save_ifile_comment,
    #vslice_by_idim,
)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str,
    array_default_int, array_default_float,
    get_print_card_size,
    #print_card_8_comment, print_card_16_comment,
)
from .utils import get_density_from_material # , expanded_mass_material_id

#from .shell_coords import element_coordinate_system, material_coordinate_system
#from .shell_utils import (
    #tri_area, tri_area_centroid_normal, tri_centroid,
    #quad_area, quad_area_centroid_normal, quad_centroid,
    #shell_mass_per_area, shell_mass_per_area_breakdown, shell_nonstructural_mass,
    #shell_thickness, shell_total_thickness,
#)
#from .shell_quality import tri_quality_nodes, quad_quality_nodes
from pyNastran.dev.bdf_vectorized3.utils import hstack_msg

NUMPY_INTS = {'int32', 'int64'}
NUMPY_FLOATS = {'float32', 'float64'}


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.cards.materials import MAT1, MAT8
    #from pyNastran.dev.bdf_vectorized3.cards.grid import GRID


def shell_materials(model: BDF) -> list[MAT1 | MAT8]:
    if model.is_thermal:
        return [model.mat4, model.mat5]
    return [model.mat1, model.mat2, model.mat8, model.mat9]


class PSHELL(Property):
    _show_attributes = ['property_id', 'material_id', 't',
                        'twelveIt3', 'tst', 'nsm', 'z']
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.zeros((0, 4), dtype='int32')
        self.t = np.array([], dtype='float64')
        self.twelveIt3 = np.array([], dtype='float64')
        self.tst = np.array([], dtype='float64')
        self.nsm = np.array([], dtype='float64')
        self.z = np.zeros((0, 2), dtype='float64')

    def add(self, pid: int, mid1: int=-1, t: float=0.0,
            mid2: int=-1, twelveIt3: float=1.0,
            mid3: int=-1, tst: float=0.833333, nsm: float=0.0,
            z1: Optional[float]=None, z2: Optional[float]=None, mid4: int=-1,
            comment: str='') -> int:
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
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        property_id = np.zeros(ncards, dtype=idtype)
        material_id = np.zeros((ncards, 4), dtype=idtype)
        t = np.zeros(ncards, dtype='float64')
        twelveIt3 = np.zeros(ncards, dtype='float64')
        tst = np.zeros(ncards, dtype='float64')
        nsm = np.zeros(ncards, dtype='float64')
        z = np.zeros((ncards, 2), dtype='float64')

        comment = {}
        for icard, card in enumerate(self.cards):
            (pid, mid1, ti,
             mid2, twelveIt3i, mid3, tsti, nsmi, z1, z2, mid4,
             commenti) = card

            property_id[icard] = pid
            material_id[icard, :] = [mid if isinstance(mid, integer_types) else -1
                                     for mid in [mid1, mid2, mid3, mid4]]
            t[icard] = ti
            twelveIt3[icard] = twelveIt3i
            tst[icard] = tsti
            nsm[icard] = nsmi
            z[icard] = [z1, z2]
            if commenti:
                comment[pid] = commenti
        self._save(property_id, material_id, t, twelveIt3, tst, nsm, z,
                   comment=comment)
        self.sort()
        self.cards = []

    def _save(self, property_id, material_id, t, twelveIt3, tst, nsm, z,
              ifile=None, comment: Optional[dict[int, str]]=None) -> None:
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.property_id) != 0:
            raise NotImplementedError()

        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.material_id = material_id
        self.t = t
        self.twelveIt3 = twelveIt3
        self.tst = tst
        self.nsm = nsm
        self.z = z
        self.n = ncards

    def set_used(self, used_dict: [str, list[np.ndarray]]) -> None:
        material_id = np.unique(self.material_id.flatten())
        material_id = material_id[material_id > 0]
        used_dict['material_id'].append(material_id)

    def convert(self, xyz_scale: float=1.0,
                nsm_per_area_scale: float=1.0, **kwargs):
        self.t *= xyz_scale
        self.z *= xyz_scale
        self.nsm *= nsm_per_area_scale

    def __apply_slice__(self, prop: PSHELL, i: np.ndarray) -> None:  # ignore[override]
        prop.n = len(i)
        self._slice_comment(prop, i)
        prop.ifile = self.ifile[i]
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
        mids = np.unique(mids)
        material_ids = np.unique(self.material_id.ravel())
        if -1 == material_ids[0]:
            material_ids = material_ids[1:]
        geom_check(self,
                   missing,
                   material_id=(mids, material_ids))

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.material_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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

    def get_Qbar_matrix(self, mid: int, theta: float=0.) -> np.ndarray:
        """theta must be in radians"""
        S2 = get_mat_s33(self.model, mid, self.allowed_materials)
        T = get_2d_plate_transform(theta)
        #Tt = np.transpose(T, axes=0)
        #Tinv = np.linalg.inv(T)
        #Qbar = np.linalg.multi_dot([Tinv, Q, Tinv.T])
        #ST = np.einsum('ijk,ilm->ijm', S2, T)
        #Sbarv = np.einsum('ijk,ilm->ikm', T, ST)
        Sbar_list = [np.linalg.multi_dot([Ti.T, Si, Ti]) for Ti, Si in zip(T, S2)]
        Qbar_list = [np.linalg.inv(Sbari) for Sbari in Sbar_list]
        #Ti = T[0, :, :]
        #Si = S2[0, :, :]
        #Sbar_test = np.linalg.multi_dot([Ti.T, Si, Ti])
        #Qbarv = np.linalg.inv(Sbarv)
        Qbar = np.array(Qbar_list, dtype=S2.dtype)
        return Qbar

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


class PAABSF(Property):
    """
    +--------+-----+--------+--------+---+---+---+---+------+
    |    1   |  2  |    3   |    4   | 5 | 6 | 7 | 8 |   9  |
    +========+=====+========+========+===+===+===+===+=======
    | PAABSF | PID | TZREID | TZIMID | S | A | B | K | RHOC |
    +--------+-----+--------+--------+---+---+---+---+------+
    """
    _show_attributes = ['property_id', 'table_reactance_real', 'table_reactance_imag',
                        's', 'a', 'b', 'k', 'rho_c']
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.table_reactance_real = np.array([], dtype='int32')
        self.table_reactance_imag = np.array([], dtype='int32')
        self.s = np.array([], dtype='float64')
        self.a = np.array([], dtype='float64')
        self.b = np.array([], dtype='float64')
        self.k = np.array([], dtype='float64')
        self.rho_c = np.array([], dtype='float64')

    def add(self, pid: int,
            table_reactance_real: int=0,
            table_reactance_imag: int=0,
            s: float=1.0, a: float=1.0, b: float=0.0,
            k: float=0.0, rhoc: float=1.0,
            comment: str=''):
        """
        Creates a PAABSF card

        Parameters
        ----------
        pid : int
            Property identification number of a CAABSF. (Integer > 0)
        TZREID : int; default=None
            TABLEDi id that defines the resistance as a function of
            frequency. The real part of the impedence.
        TZIMID : int; default=None
            TABLEDi id that defines the reactance as a function of
            frequency. The imaginary part of the impedance.
        S : float; default=1.0
            Impedance scale factor.
        A : float; default=1.0
            Area factor when 1 or 2 grid points are specified on the
            CAABSF entry.
        B : float; default=0.0
            Equivalent structural damping coefficient.
        K : float; default=0.0
            Equivalent structural stiffness coefficient.
        RHOC : float; default=1.0
            Constant used in data recovery for calculating an absorption
            coefficient. RHO is the media density, and C is the speed of
            sound in the media.

        """
        table_reactance_real = 0 if table_reactance_real is None else table_reactance_real
        table_reactance_imag = 0 if table_reactance_imag is None else table_reactance_imag
        self.cards.append((pid, table_reactance_real, table_reactance_imag, s, a, b, k, rhoc))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a PAABSF card from ``BDF.add_card(...)``

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
        table_reactance_real = integer_or_blank(card, 2, 'tzreid', default=0)
        table_reactance_imag = integer_or_blank(card, 3, 'tzimid', default=0)
        s = double_or_blank(card, 4, 's', default=1.0)
        a = double_or_blank(card, 5, 'a', default=1.0)
        b = double_or_blank(card, 6, 'b', default=0.0)
        k = double_or_blank(card, 7, 'k', default=0.0)
        rhoc = double_or_blank(card, 8, 'rhoc', default=1.0)
        assert len(card) <= 9, f'len(PAABSF card) = {len(card):d}\ncard={card}'
        #return PAABSF(pid, tzreid, tzimid, s, a, b, k, rhoc, comment=comment)
        self.cards.append((pid, table_reactance_real, table_reactance_imag, s, a, b, k, rhoc))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        table_reactance_real = np.zeros(ncards, dtype='int32')
        table_reactance_imag = np.zeros(ncards, dtype='int32')
        s = np.zeros(ncards, dtype='float64')
        a = np.zeros(ncards, dtype='float64')
        b = np.zeros(ncards, dtype='float64')
        k = np.zeros(ncards, dtype='float64')
        rho_c = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (pid, table_reactance_reali, table_reactance_imagi, si, ai, bi, ki, rhoci) = card

            property_id[icard] = pid
            table_reactance_real[icard] = table_reactance_reali
            table_reactance_imag[icard] = table_reactance_imagi
            s[icard] = si
            a[icard] = ai
            b[icard] = bi
            k[icard] = ki
            rho_c[icard] = rhoci
        self._save(property_id, table_reactance_real, table_reactance_imag, s, a, b, k, rho_c)
        self.sort()
        self.cards = []

    def _save(self, property_id,
              table_reactance_real, table_reactance_imag,
              s, a, b, k, rho_c,
              ifile=None, comment: Optional[dict[int, str]]=None) -> None:
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')

        if len(self.property_id) != 0:
            raise NotImplementedError()

        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.table_reactance_real = table_reactance_real
        self.table_reactance_imag = table_reactance_imag
        self.s = s
        self.a = a
        self.b = b
        self.k = k
        self.rho_c = rho_c
        self.n = ncards

    def set_used(self, used_dict: [str, list[np.ndarray]]) -> None:
        tableds = np.unique(np.hstack([
            self.table_reactance_real,
            self.table_reactance_imag]))
        tableds = tableds[tableds > 0]
        used_dict['tabled_id'].append(tableds)

    #def convert(self, xyz_scale: float=1.0,
                #nsm_per_area_scale: float=1.0, **kwargs):
        #self.t *= xyz_scale
        #self.z *= xyz_scale
        #self.nsm *= nsm_per_area_scale

    def __apply_slice__(self, prop: PAABSF, i: np.ndarray) -> None:  # ignore[override]
        prop.n = len(i)
        self._slice_comment(prop, i)
        prop.ifile = self.ifile[i]
        prop.property_id = self.property_id[i]
        prop.table_reactance_real = self.table_reactance_real[i]
        prop.table_reactance_imag = self.table_reactance_imag[i]
        prop.s = self.s[i]
        prop.a = self.a[i]
        prop.b = self.b[i]
        prop.k = self.k[i]
        prop.rho_c = self.rho_c[i]

    #def geom_check(self, missing: dict[str, np.ndarray]):
        #mids = hstack_msg([mat.material_id for mat in self.allowed_materials],
                          #msg=f'no shell materials for {self.type}')
        #mids = np.unique(mids)
        #material_ids = np.unique(self.material_id.ravel())
        #if -1 == material_ids[0]:
            #material_ids = material_ids[1:]
        #geom_check(self,
                   #missing,
                   #material_id=(mids, material_ids))

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(),
                   self.table_reactance_real.max(),
                   self.table_reactance_imag.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        for pid, table_real, table_imag, s, a, b, k, rhoc in zip_longest(self.property_id,
                                                                         self.table_reactance_real, self.table_reactance_real,
                                                                         self.s, self.a, self.b, self.k, self.rho_c):
            list_fields = [
                'PAABSF', pid, table_real, table_imag, s, a, b, k, rhoc]
            msg = print_card(list_fields)
            bdf_file.write(msg)
        return


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


class CompositeProperty(Property):
    def clear(self) -> None:
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

    #@abstractmethod
    #def __apply_slice__(self, prop: Any, i: np.ndarray) -> None:  # pragma: no cover
        #raise NotImplementedError(f'{self.type}: __apply_slice__')
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
        assert self.lam.size > 0, str(self)
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
        #if np.isnan(rho.max()):
            #inan = np.isnan(rho)
            #mid_nan = self.material_id[inan]
            #raise RuntimeError(f'material_id={mid_nan} has nan rho')

        rho_t = rho * thickness
        mass_per_area = np.zeros(nproperties, dtype='float64') + nsm
        for i, scale, ilayer in zip(count(), symmetry_scale, ilayers):
            idim0, idim1 = ilayer
            mass_per_areai = scale * rho_t[idim0:idim1].sum()
            #assert not np.isnan(mass_per_areai), (i, scale, ilayer, nsm, rho, thickness)
            mass_per_area[i] += mass_per_areai
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

        #if np.isnan(mass_per_area.max()):
            #inan = np.isnan(mass_per_area)
            #pid_nan = self.property_id[inan]
            #raise RuntimeError(f'property_id={pid_nan} has nan mass_per_area')
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
    _skip_equality_check = True  # assume unequal
    _show_attributes = [
        'property_id', 'z0', 'nsm', 'failure_theory', 'tref',
        'ge', 'lam', 'nlayer', 'material_id', 'thickness', 'sout', 'theta']
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.z0 = np.array([], dtype='float64')
        self.nsm = np.array([], dtype='float64')
        self.shear_bonding = np.array([], dtype='float64')

        # 'HILL' for the Hill theory.
        # 'HOFF' for the Hoffman theory.
        # 'TSAI' for the Tsai-Wu theory.
        # 'STRN' for the Maximum Strain theory.
        self.failure_theory = np.array([], dtype='|U8')
        self.tref = np.array([], dtype='float64')
        self.ge = np.array([], dtype='float64')
        self.lam = np.array([], dtype='|U8')

        self.nlayer = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.thickness = np.array([], dtype='float64')
        self.sout = np.array([], dtype='|U4') # YES, NO, YESA
        self.theta = np.array([], dtype='float64')

    def add(self, pid: int, mids: list[int], thicknesses: list[float],
            thetas: Optional[list[float]]=None, souts: Optional[list[str]]=None,
            nsm: float=0., sb: float=0., ft: str='',
            tref: float=0., ge: float=0., lam: str='',
            z0=None, ifile: int=0, comment: str='') -> int:
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

        nmids = len(mids)
        if isinstance(thicknesses, float_types):
            thicknesses = [thicknesses] * nmids

        if thetas is None:
            thetas = [0.] * nmids
        elif isinstance(thetas, float_types):
            thetas = [thetas] * nmids

        if souts is None:
            souts = ['YES'] * nmids
        elif isinstance(souts, str):
            souts = [souts] * nmids

        assert len(thicknesses) == nmids, thicknesses
        assert len(souts) == nmids, souts
        assert len(thetas) == nmids, thetas
        self.cards.append((pid, nsm, sb, ft, tref, ge, lam, z0,
                           mids, thicknesses, thetas, souts, ifile, comment))
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
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
                           mids, thicknesses, thetas, souts, ifile, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = 'int32'
        ifile = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype=idtype)
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
        comment = {}

        map_ft = {
            'HFAI': 'HFAIL',
            'HTAP': 'HTAPE',
            'HFAB': 'HFABR',
        }
        mids_list = []
        thickness_list = []
        thetas_list = []
        sout_list = []
        for icard, card in enumerate(self.cards):
            (pid, nsmi, shear_bondingi, fti, trefi, gei, lami, z0i,
             mids, thicknesses, thetas, souts, ifilei, commenti) = card

            nsm[icard] = nsmi
            shear_bonding[icard] = shear_bondingi
            failure_theory[icard] = fti
            tref[icard] = trefi
            ge[icard] = gei

            #'MEM' All plies must be specified, but only membrane terms (MID1 on the
            #      derived PSHELL entry) are computed.
            #'BEND' All plies must be specified, but only bending terms (MID2 on the derived
            #       PSHELL entry) are computed.
            #'SMEAR' All plies must be specified, stacking sequence is ignored MID1=MID2 on
            #        the derived PSHELL entry and MID3, MID4 and TS/T and 12I/T**3
            #        terms are set to zero).
            #'SMCORE' All plies must be specified, with the last ply specifying core properties and
            #         the previous plies specifying face sheet properties. The stiffness matrix is
            #         computed by placing half the face sheet thicknesses above the core and the
            #         other half below with the result that the laminate is symmetric about the
            #         mid-plane of the core. Stacking sequence is ignored in calculating the face
            #         sheet stiffness.
            lam[icard] = lami
            assert lami in {'', 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE'}, f'pid={pid} laminate={lami!r}'

            #'HILL' for the Hill theory.
            #'HOFF' for the Hoffman theory.
            #'TSAI' for the Tsai-Wu theory.
            #'STRN' for the Maximum Strain theory.
            #'HFAIL' for the Hashin failure criterion
            #'HTAPE' for the Hashin tape criterion
            #'HFABR' for the Hashin fabric criterion
            fti = map_ft.get(fti, fti)
            assert fti in {'', 'HILL', 'HOFF', 'STRN', 'TSAI', 'HFAIL', 'HFABR', 'HTAPE'}, f'pid={pid} failure_theory={fti!r}'
            nlayersi = len(mids)

            ifile[icard] = ifilei
            if commenti:
                comment[pid] = commenti
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
                   z0, nsm, shear_bonding, failure_theory, tref, ge, lam,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, property_id,
              nlayer, material_id, thickness, sout, theta,
              z0, nsm, shear_bonding, failure_theory, tref, ge, lam,
              ifile=None, comment: Optional[dict[int, str]]=None):
        assert isinstance(lam, np.ndarray), lam
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')

        if ncards:
            ifile = np.hstack([self.ifile, ifile])
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

        save_ifile_comment(self, ifile, comment)
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
        assert len(self.ifile) == len(self.property_id)

    def update_layers(self, property_id: np.ndarray,
                      thickness: Optional[np.ndarray]=None,
                      theta: Optional[np.ndarray]=None,
                      ilayer: Optional[np.ndarray]=None) -> None:
        assert thickness is not None or theta is not None, (thickness, theta)
        if isinstance(property_id, integer_types):
            property_id = np.array([property_id])
            if thickness is not None:
                thickness = np.array([thickness])
            if theta is not None:
                theta = np.array([theta])

        if theta is not None and thickness is not None:
            assert theta.shape == thickness.shape, (theta.shape, thickness.shape)

        uproperty_id = np.unique(property_id)
        assert len(property_id) == len(uproperty_id), f'property_id={property_id} uproperty_id={uproperty_id}'
        ipid = np.searchsorted(self.property_id, property_id)
        assert np.array_equal(self.property_id[ipid], property_id)

        if thickness is not None:
            assert thickness.ndim == 2, thickness.shape
            nlayer = self.nlayer[ipid]
            min_layers = nlayer.min()
            ndimensions = thickness.shape[1]
            if min_layers < ndimensions:
                imin = np.where(nlayer == min_layers)[0]
                pidi = property_id[imin]
                raise RuntimeError(f'too many dimensions for pid={pidi} '
                                   f'min(layers)={min_layers}; thickness.shape={thickness.shape}')

            ilayer_ = self.ilayer[ipid, :]
            for pid, (i0, i1), thicknessi in zip(property_id, ilayer_, thickness):
                thicknesses = self.thickness[i0:i1]
                if ilayer is None:
                    assert len(thicknesses) == len(thicknessi), (thicknesses, thicknessi)
                    thicknesses[:] = thicknessi
                else:
                    assert len(ilayer) == len(thicknessi), (ilayer, thicknessi)
                    thicknesses[ilayer] = thicknessi

        if theta is not None:
            assert theta.ndim == 2, theta.shape
            ilayer_ = self.ilayer[ipid, :]
            for pid, (i0, i1), thetai in zip(property_id, ilayer_, theta):
                thetas = self.theta[i0:i1]
                if ilayer is None:
                    assert len(thetas) == len(thicknessi), (thetas, thicknessi)
                    thetas[:] = thetai
                else:
                    assert len(ilayer) == len(thicknessi), (ilayer, thicknessi)
                    thetas[ilayer] = thetai

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['material_id'].append(self.material_id)

    def convert(self, xyz_scale: float=1.0,
                nsm_per_area_scale: float=1.0,
                temperature_scale: float=1.0, **kwargs):
        self.thickness *= xyz_scale
        self.z0 *= xyz_scale
        # self.thickness ## TODO: ???
        self.tref *= temperature_scale
        self.nsm *= nsm_per_area_scale

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
            except NotImplementedError:  # pragma: no cover
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
    def nplies(self, nplies: np.ndarray) -> None:
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
        print(f'iprop = {iprop}')
        assert len(self.ifile) == len(self.property_id)
        #assert len(self.node_id) > 0, self.node_id
        #i = np.searchsorted(self.node_id, node_id)
        prop = self.slice_card_by_index(iprop)
        assert len(prop.property_id) == len(prop.ifile)
        return prop

    def set_used(self, used_dict: [str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)

    def __apply_slice__(self, prop: PCOMP, i: np.ndarray) -> None:  # ignore[override]
        assert self.nlayer.sum() == len(self.thickness)
        assert len(self.ifile) == len(self.property_id)
        self._slice_comment(prop, i)
        prop.ifile = self.ifile[i]
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
        assert len(prop.ifile) == len(prop.property_id)

    @property
    def sb(self):
        return self.shear_bonding
    @sb.setter
    def sb(self, sb):
        self.shear_bonding = sb

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.material_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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
    _skip_equality_check = True  # assume unequal
    _show_attributes = [
        'property_id', 'z0', 'nsm', 'failure_theory', 'tref',
        'ge', 'lam', 'nlayer',
        'global_ply_id', 'material_id', 'thickness', 'sout', 'theta',
        'shear_bonding',
    ]
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.z0 = np.array([], dtype='float64')
        self.nsm = np.array([], dtype='float64')
        self.shear_bonding = np.array([], dtype='float64')

        # 'HILL' for the Hill theory.
        # 'HOFF' for the Hoffman theory.
        # 'TSAI' for the Tsai-Wu theory.
        # 'STRN' for the Maximum Strain theory.
        self.failure_theory = np.array([], dtype='|U8')
        self.tref = np.array([], dtype='float64')
        self.ge = np.array([], dtype='float64')
        self.lam = np.array([], dtype='|U6')

        self.nlayer = np.array([], dtype='int32')
        self.global_ply_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.thickness = np.array([], dtype='float64')
        self.sout = np.array([], dtype='|U4') # YES, NO, YESA
        self.theta = np.array([], dtype='float64')

    def add(self, pid: int,
            global_ply_ids: list[int], mids: list[int], thicknesses: list[float],
            thetas=None, souts=None,
            nsm: float=0.0, sb: float=0.0, ft: str='',
            tref: float=0.0, ge: float=0.0, lam: str='', z0=None,
            ifile: int=0, comment: str='') -> int:
        lam = lam if lam is not None else ''
        ft = ft if ft is not None else ''
        z0 = z0 if z0 is not None else np.nan

        nlayers = len(thicknesses)
        if thetas is None:
            thetas = [0.] * nlayers
        if souts is None:
            souts = ['NO'] * nlayers

        cardi = (
            pid, nsm, sb, ft, tref, ge, lam, z0,
            mids, thicknesses, thetas, souts, global_ply_ids,
            ifile, comment)
        self.cards.append(cardi)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        pid = integer(card, 1, 'pid')
        # z0 will be calculated later
        nsm = double_or_blank(card, 3, 'nsm', default=0.0)
        sb = double_or_blank(card, 4, 'sb', default=0.0)
        ft = string_or_blank(card, 5, 'ft', default='')
        tref = double_or_blank(card, 6, 'tref', default=0.0)
        ge = double_or_blank(card, 7, 'ge', default=0.0)
        lam = string_or_blank(card, 8, 'lam', default='')

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
            ifile, comment)
        self.cards.append(cardi)
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        nsm = np.zeros(ncards, dtype='float64')
        shear_bonding = np.zeros(ncards, dtype='float64')
        failure_theory = np.zeros(ncards, dtype='|U8')

        tref = np.zeros(ncards, dtype='float64')
        lam = np.zeros(ncards, dtype='|U5')
        z0 = np.zeros(ncards, dtype='float64')

        nlayer = np.zeros(ncards, dtype='int32')
        ge = np.zeros(ncards, dtype='float64')
        comment = {}

        global_ply_ids_list = []
        mids_list = []
        thickness_list = []
        thetas_list = []
        sout_list = []

        map_ft = {
            'HFAI': 'HFAIL',
            'HTAP': 'HTAPE',
            'HFAB': 'HFABR',
        }
        for icard, card in enumerate(self.cards):
            (pid, nsmi, sbi, fti, trefi, gei, lami, z0i,
             mids, thicknesses, thetas, souts, global_ply_ids,
             ifilei, commenti) = card

            ifile[icard] = ifilei
            if commenti:
                comment[pid] = commenti
            failure_theory[icard] = fti
            shear_bonding[icard] = sbi
            tref[icard] = trefi
            ge[icard] = gei
            lam[icard] = lami

            # SYM not listed in QRG
            assert lami in {'', 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE'}, f'pid={pid} laminate={lami!r}'

            #'HILL' for the Hill theory
            #'HOFF' for the Hoffman theory
            #'TSAI' for the Tsai-Wu theory
            #'STRN' for the Maximum Strain theory
            #'HFAIL' for the Hashin failure criterion
            #'HTAPE' for the Hashin tape criterion
            #'HFABR' for the Hashin fabric criterion
            fti = map_ft.get(fti, fti)
            assert fti in {'', 'HILL', 'HOFF', 'STRN', 'TSAI', 'HFAIL', 'HFABR', 'HTAPE'}, f'pid={pid} failure_theory={fti!r}'

            nlayeri = len(mids)
            property_id[icard] = pid
            global_ply_ids_list.extend(global_ply_ids)
            mids_list.extend(mids)
            thickness_list.extend(thicknesses)
            sout_list.extend(souts)
            thetas_list.extend(thetas)
            z0[icard] = z0i
            nlayer[icard] = nlayeri
            nsm[icard] = nsmi

        global_ply_id = np.array(global_ply_ids_list, dtype='int32')
        material_id = np.array(mids_list, dtype='int32')
        thickness = np.array(thickness_list, dtype='float64')
        sout = np.array(sout_list, dtype='|U4') # YES, NO, YESA
        theta = np.array(thetas_list, dtype='float64')
        self._save(property_id, nsm, shear_bonding, failure_theory, tref, ge, lam, z0,
                   nlayer, global_ply_id, material_id, thickness, sout, theta,
                   ifile=ifile, comment=comment)
        assert len(mids_list) == self.nlayer.sum()
        self.sort()
        self.cards = []

    def _save(self, property_id, nsm, shear_bonding, failure_theory, tref, ge, lam, z0,
              nlayer, global_ply_id, material_id, thickness, sout, theta,
              ifile=None, comment: Optional[dict[int, str]]=None) -> None:
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.property_id) != 0:
            raise RuntimeError(f'stacking of {self.type} is not supported')

        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.nsm = nsm
        self.shear_bonding = shear_bonding
        self.failure_theory = failure_theory

        self.tref = tref
        self.lam = lam
        self.z0 = z0
        self.ge = ge
        #self.group = group

        self.nlayer = nlayer
        self.global_ply_id = global_ply_id
        self.material_id = material_id
        self.thickness = thickness
        self.sout = sout
        self.theta = theta

    def __apply_slice__(self, prop: PCOMPG, i: np.ndarray) -> None:  # ignore[override]
        prop.n = len(i)
        self._slice_comment(prop, i)
        prop.ifile = self.ifile[i]
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

        ilayer = self.ilayer # [i, :]
        prop.material_id = hslice_by_idim(i, ilayer, self.material_id)
        prop.sout = hslice_by_idim(i, ilayer, self.sout)
        prop.theta = hslice_by_idim(i, ilayer, self.theta)
        prop.thickness = hslice_by_idim(i, ilayer, self.thickness)
        prop.global_ply_id = hslice_by_idim(i, ilayer, self.global_ply_id)
        prop.nlayer = self.nlayer[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['material_id'].append(self.material_id)

    def convert(self, xyz_scale: float=1.0,
                nsm_per_area_scale: float=1.0,
                temperature_scale: float=1.0, **kwargs):
        self.thickness *= xyz_scale
        self.z0 *= xyz_scale
        # self.thickness ## TODO: ???
        self.tref *= temperature_scale
        self.nsm *= nsm_per_area_scale

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.global_ply_id.max(),
                   self.material_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        property_id = array_str(self.property_id, size=size)

        nsms = array_default_float(self.nsm, default=0., size=size, is_double=False)
        shear_bondings = array_default_float(self.shear_bonding, default=0., size=size, is_double=False)
        trefs = array_default_float(self.tref, default=0., size=size, is_double=False)
        ges = array_default_float(self.ge, default=0., size=size, is_double=False)
        for pid, nsm, sb, failure_theory, ge, tref, lam, z0, \
            ilayer in zip_longest(property_id, nsms, shear_bondings,
                                  self.failure_theory, trefs, ges, self.lam, self.z0,
                                  self.ilayer):
            #nsm = set_blank_if_default(nsm, 0.0)
            #sb = set_blank_if_default(sb, 0.0)
            #tref = set_blank_if_default(tref, 0.0)
            #ge = set_blank_if_default(ge, 0.0)

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
            list_fields = ['PCOMPG', pid, z02, nsm, sb, failure_theory, tref, ge, lam]
            for (global_id, mid, t, theta, sout) in zip(global_ids, material_ids, thicknesses,
                                                        thetas, souts):
                #theta = set_blank_if_default(theta, 0.0)
                str_sout = set_blank_if_default(sout, 'NO')
                list_fields += [global_id, mid, t, theta, str_sout, None, None, None]

            bdf_file.write(print_card(list_fields))
        return


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
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.coord_id = np.array([], dtype='int32')
        self.stress_strain_output_location = np.array([], dtype='|U4')
        self.thickness = np.array([], dtype='float64')

    def add(self, pid: int, mid: int, cid: int=0,
            stress_strain_output_location: str='GRID', thickness: float=np.nan,
            ifile: int=0, comment: str='') -> int:
        """Creates a PLPLANE card"""
        self.cards.append((pid, mid, cid, stress_strain_output_location, thickness, ifile, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
        self.cards.append((pid, mid, cid, stress_strain_output_location, thickness, ifile, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        self.ifile = np.zeros(ncards, dtype='int32')
        self.property_id = np.zeros(ncards, dtype='int32')
        self.material_id = np.zeros(ncards, dtype='int32')
        self.coord_id = np.zeros(ncards, dtype='int32')
        self.stress_strain_output_location = np.zeros(ncards, dtype='|U4')
        self.thickness = np.zeros(ncards, dtype='float64')
        self.comment = {}

        for icard, card in enumerate(self.cards):
            (pid, mid, cid, stress_strain_output_location, thickness, ifilei, commenti) = card
            assert len(stress_strain_output_location) <= 4, stress_strain_output_location

            self.ifile[icard] = ifilei
            if commenti:
                self.comment[pid] = commenti
            self.property_id[icard] = pid
            self.material_id[icard] = mid
            self.coord_id[icard] = cid
            self.stress_strain_output_location[icard] = stress_strain_output_location
            self.thickness[icard] = thickness
        self.sort()
        self.cards = []

    def set_used(self, used_dict: [str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)
        used_dict['coord_id'].append(self.coord_id)

    def __apply_slice__(self, prop: PLPLANE, i: np.ndarray) -> None:  # ignore[override]
        prop.n = len(i)
        self._slice_comment(prop, i)
        prop.ifile = self.ifile[i]
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]
        prop.coord_id = self.coord_id[i]
        prop.stress_strain_output_location = self.stress_strain_output_location[i]
        prop.thickness = self.thickness[i]

    def convert(self, xyz_scale: float=1.0, **kwargs):
        self.thickness *= xyz_scale

    #def write_file_8(self, bdf_file: TextIOLike,
    #               write_card_header: bool=False) -> None:
    #    self.write_file(bdf_file, size=8, is_double=False,
    #                    write_card_header=write_card_header)
#
#    def write_file_16(self, bdf_file: TextIOLike,
#                      is_double=False, write_card_header: bool=False) -> None:
#        self.write_file(bdf_file, size=16, is_double=is_double,
#                        write_card_header=write_card_header)

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.material_id.max(), self.coord_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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


class PSHLN1(Property):
    """
    Fully Nonlinear Plane Element Properties (SOL 601)
    Defines the properties of a fully nonlinear
    (i.e., large strain and large rotation)
    hyperelastic plane strain or axisymmetric element.

    +--------+------+------+--------+----------+--------+
    |    1   |  2   |  3   |   4    |     5    |    6   |
    +========+======+======+========+==========+========+
    | PSHLN1 |  PID | MID1 |  MID2  | ANALYSIS |        |
    +--------+------+------+--------+----------+--------+
    |        | 'C3' | BEH3 |  INT3  |   BEH3H  | INT3H  |
    +--------+------+------+--------+----------+--------+
    |        | 'C4' | BEH4 |  INT4  |   BEH4H  | INT4H  |
    +--------+------+------+--------+----------+--------+
    |        | 'C6' | BEH6 |  INT6  |   BEH6H  | INT6H  |
    +--------+------+------+--------+----------+--------+
    |        | 'C8' | BEH8 |  INT8  |   BEH8H  | INT8H  |
    +--------+------+------+--------+----------+--------+
    MSC

    PSHLN1       100     100     100     ISH
                  C3    DCTN     LDK     DCT       L
                  C4     DCT       L     DCT       L
                  C8      MB       Q     DCT       Q
    """
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.zeros((0, 2), dtype='int32')
        #self.thickness = np.array([], dtype='float64')
        self.analysis = np.array([], dtype='|U8')

        self.beh = np.zeros((0, 4), dtype='|U8')
        self.beh_h = np.zeros((0, 4), dtype='|U8')
        self.integration = np.zeros((0, 4), dtype='|U8')
        self.integration_h = np.zeros((0, 4), dtype='|U8')
        #self.coord_id = np.array([], dtype='int32')
        #self.stress_strain_output_location = np.array([], dtype='|U4')
        #self.thickness = np.array([], dtype='float64')

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
        ##self.property_id = np.hstack([self.property_id, pid])
        ##mids = [mid if mid is not None else -1
                ##for mid in [mid1, mid2, mid3, mid4]]
        ##self.material_id = np.vstack([self.material_id, mids])
        ##self.t = np.hstack([self.t, t])
        ##self.twelveIt3 = np.hstack([self.twelveIt3, twelveIt3])
        ##self.tst = np.hstack([self.tst, tst])
        ##self.nsm = np.hstack([self.nsm, nsm])

        ##self.z = np.vstack([self.z, [z1, z2]])
        #self.n += 1

    def add(self, pid: int, mid1: int=0, mid2: int=0, analysis: str='ISH',
            behx=None, integration=None, behxh=None, integration_h=None,
            comment: str='') -> int:
        self.cards.append((pid, (mid1, mid2), analysis,
                           behx, integration, behxh, integration_h, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a PSHLN2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        #| PSHLN2  | PID  | MID  | DIRECT |    T   | ANALYSIS |
        #|         | 'C3' | BEH3 |  INT3  |  BEH3H | INT3H    |
        #|         | 'C4' | BEH4 |  INT4  |  BEH4H | INT4H    |
        #|         | 'C6' | BEH6 |  INT6  |  BEH6H | INT6H    |
        #|         | 'C8' | BEH8 |  INT8  |  BEH8H | INT8H    |

        pid = integer(card, 1, 'pid')
        mid1 = integer_or_blank(card, 2, 'mid1', default=0)  # MATHE, MATHP
        mid2 = integer_or_blank(card, 3, 'mid2', default=0)

        #: Analysis type.
        #  'IS' - Implicit structural elements are being referred to.
        #  'IH' - Implicit heat analysis elements are being referred to.
        #  'ISH' - Implicit structural and heat elements are being referred to. (Character Default ISH)
        analysis = string_or_blank(card, 4, 'analysis', default='ISH')
        i = 9
        beh3 = None
        int3 = None
        beh4 = None
        int4 = None
        beh6 = None
        int6 = None
        beh8 = None
        int8 = None

        beh3h = None
        int3h = None
        beh4h = None
        int4h = None
        beh6h = None
        int6h = None
        beh8h = None
        int8h = None
        while i < len(card):
            code = string(card, i, 'code')
            # C3 : applies to elements with three corner grids
            # C4 : applies to elements with four corner grids
            # C6 : applies to elements with three corner grids and three midside grids
            # C8 : applies to elements with four corner grids and four midside grids
            # BEHi Element structural behavior. See Remark 7. (Default: PLSTRN for BEH3, BEH4, BEH6, and BEH8)
            # INTi Integration scheme. See Remarks 7. and 10. (Default: L for INT3, INT4, Q for INT6 and INT8)
            # BEHiH Element heat behavior. (Default: PLSTRN for BEH3H, BEH4H, BEH6H, and BEH8H.)
            # INTiH Integration scheme.    (Default: L for INT3H, L for INT4H, Q for INT6H and INT8H)
            if code == 'C3':
                beh3 = string_or_blank(card, i+1, 'BEH3', default='DCTN')
                int3 = string_or_blank(card, i+2, 'INT3', default='LDK')
                beh3h = string_or_blank(card, i+3, 'BEH3H', default='DCT')
                int3h = string_or_blank(card, i+4, 'INT3H', default='L')
            elif code == 'C4':
                beh4 = string_or_blank(card, i+1, 'BEH4', default='DCT')
                int4 = string_or_blank(card, i+2, 'INT4', default='L')
                beh4h = string_or_blank(card, i+3, 'BEH4H', default='DCT')
                int4h = string_or_blank(card, i+4, 'INT4H', default='L')
            elif code == 'C6':
                beh6 = string_or_blank(card, i+1, 'BEH6', default='MB')
                int6 = string_or_blank(card, i+2, 'INT6', default='Q')
                beh6h = string_or_blank(card, i+3, 'BEH6H', default='MB')
                int6h = string_or_blank(card, i+4, 'INT6H', default='Q')
            elif code == 'C8':
                beh8 = string_or_blank(card, i+1, 'BEH8', default='DCT')
                int8 = string_or_blank(card, i+2, 'INT8', default='QRI')
                beh8h = string_or_blank(card, i+3, 'BEH8H', default='DCT')
                int8h = string_or_blank(card, i+4, 'INT8H', default='Q')
            else:
                raise NotImplementedError(f'PSHLN2 code={code!r}')
            i += 8
        behx = [beh3, beh4, beh6, beh8]
        behxh = [beh3h, beh4h, beh6h, beh8h]
        integration = [int3, int4, int6, int8]
        integration_h = [int3h, int4h, int6h, int8h]
        assert len(card) <= 38, f'len(PSHLN1 card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, (mid1, mid2), analysis,
                           behx, integration, behxh, integration_h, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros((ncards, 2), dtype='int32')
        analysis = np.zeros(ncards, dtype='|U8')
        beh = np.zeros((ncards, 4), dtype='|U8')
        beh_h = np.zeros((ncards, 4), dtype='|U8')
        integration = np.zeros((ncards, 4), dtype='|U8')
        integration_h = np.zeros((ncards, 4), dtype='|U8')

        for icard, card in enumerate(self.cards):
            (pid, mids, analysisi,
             behx, integrationx, behxh, integration_hx, comment) = card

            property_id[icard] = pid
            material_id[icard] = mids
            analysis[icard] = analysisi
            if behx is None:
                behx = [''] * 4
            if integrationx is None:
                integrationx = [''] * 4
            if behxh is None:
                behxh = [''] * 4
            if integration_hx is None:
                integration_hx = [''] * 4

            beh[icard] = [val if val is not None else ''
                          for val in behx]
            integration[icard] = [val if val is not None else ''
                                  for val in integrationx]
            beh_h[icard] = [val if val is not None else ''
                            for val in behxh]
            integration_h[icard] = [val if val is not None else ''
                                    for val in integration_hx]
        self._save(property_id, material_id, analysis, beh, beh_h, integration, integration_h)
        self.sort()
        self.cards = []

    def _save(self, property_id, material_id, analysis, beh, beh_h, integration, integration_h,
              ifile=None, comment: Optional[dict[int, str]]=None) -> None:
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.property_id) != 0:
            raise RuntimeError(f'stacking of {self.type} is not supported')
        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.material_id = material_id
        self.analysis = analysis
        self.beh = beh
        self.beh_h = beh_h
        self.integration = integration
        self.integration_h = integration_h

    def set_used(self, used_dict: [str, list[np.ndarray]]) -> None:
        material_id = self.material_id.flatten()
        material_id = material_id[material_id > 0]
        if len(material_id):
            used_dict['material_id'].append(material_id)

    def __apply_slice__(self, prop: PSHLN1, i: np.ndarray) -> None:
        prop.n = len(i)
        self._slice_comment(prop, i)
        prop.ifile = self.ifile[i]
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i, :]
        prop.analysis = self.analysis[i]

        prop.beh = self.beh[i, :]
        prop.beh_h = self.beh_h[i, :]
        prop.integration = self.integration[i, :]
        prop.integration_h = self.integration_h[i, :]

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.material_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        property_ids = array_str(self.property_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        #directs = array_default_int(self.direct, default=1, size=size)
        #inan = np.isnan(self.thickness)

        #self.property_id = np.zeros(ncards, dtype='int32')
        #self.material_id = np.zeros(ncards, dtype='int32')
        #self.direct = np.zeros(ncards, dtype='int32')
        #self.analysis = np.zeros(ncards, dtype='|U8')
        #self.thickness = np.zeros(ncards, dtype='float64')
        #self.beh = np.zeros((ncards, 2), dtype='|U8')
        #self.beh_h = np.zeros((ncards, 2), dtype='|U8')
        #self.integration = np.zeros((ncards, 2), dtype='|U8')
        #self.integration_h = np.zeros((ncards, 2), dtype='|U8')

        codes = ['C3', 'C4', 'C6', 'C8']
        for pid, (mid1, mid2), analysis, beh, integration, beh_h, integration_h in zip_longest(
            property_ids, material_ids, self.analysis,
            self.beh, self.integration, self.beh_h, self.integration_h):
            list_fields = ['PSHLN1', pid, mid1, mid2, analysis, None, None, None, None]
            #values = (beh, integration, beh_h, integration_h)
            for code, behx, intx, behxh, intxh in zip(codes, beh, integration, beh_h, integration_h):
                if behx == '' and intx == '' and behx == '' and intxh == '':
                    continue
                list_fields.extend([code, behx, intx, behxh, intxh, '', '', ''])
            bdf_file.write(print_card(list_fields))
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

    #def total_thickness(self) -> np.ndarray:
        #return self.thickness


class PSHLN2(Property):
    """
    Fully Nonlinear Plane Element Properties (SOL 601)
    Defines the properties of a fully nonlinear
    (i.e., large strain and large rotation)
    hyperelastic plane strain or axisymmetric element.

    +---------+------+------+--------+--------+----------+
    |    1    |  2   |  3   |   4    |    5   |     6    |
    +=========+======+======+========+========+==========+
    | PSHLN2  | PID  | MID  | DIRECT |    T   | ANALYSIS |
    +---------+------+------+--------+--------+----------+
    |         | 'C3' | BEH3 |  INT3  |  BEH3H | INT3H    |
    +---------+------+------+--------+--------+----------+
    |         | 'C4' | BEH4 |  INT4  |  BEH4H | INT4H    |
    +---------+------+------+--------+--------+----------+
    |         | 'C6' | BEH6 |  INT6  |  BEH6H | INT6H    |
    +---------+------+------+--------+--------+----------+
    |         | 'C8' | BEH8 |  INT8  |  BEH8H | INT8H    |
    +---------+------+------+--------+--------+----------+
    MSC

    """
    @Property.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.material_id = np.array([], dtype='int32')
        self.thickness = np.array([], dtype='float64')
        self.direct = np.array([], dtype='int32')
        self.analysis = np.array([], dtype='|U8')
        self.beh = np.zeros((0, 4), dtype='|U8')
        self.beh_h = np.zeros((0, 4), dtype='|U8')
        self.integration = np.zeros((0, 4), dtype='|U8')
        self.integration_h = np.zeros((0, 4), dtype='|U8')


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
        ##self.property_id = np.hstack([self.property_id, pid])
        ##mids = [mid if mid is not None else -1
                ##for mid in [mid1, mid2, mid3, mid4]]
        ##self.material_id = np.vstack([self.material_id, mids])
        ##self.t = np.hstack([self.t, t])
        ##self.twelveIt3 = np.hstack([self.twelveIt3, twelveIt3])
        ##self.tst = np.hstack([self.tst, tst])
        ##self.nsm = np.hstack([self.nsm, nsm])

        ##self.z = np.vstack([self.z, [z1, z2]])
        #self.n += 1

    def add(self, pid: int, mid: int, direct: int=1, thickness: float=1.0,
            analysis='ISH',
            behx: list[str]=None,
            integration: list[str]=None,
            behxh: list[str]=None,
            integration_h: list[str]=None,
            comment: str='') -> int:
        self.cards.append((pid, mid, direct, thickness, analysis,
                           behx, integration, behxh, integration_h, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a PSHLN2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        if self.debug:
            self.model.log.debug(f'adding card {card}')

        #| PSHLN2  | PID  | MID  | DIRECT |    T   | ANALYSIS |
        #|         | 'C3' | BEH3 |  INT3  |  BEH3H | INT3H    |
        #|         | 'C4' | BEH4 |  INT4  |  BEH4H | INT4H    |
        #|         | 'C6' | BEH6 |  INT6  |  BEH6H | INT6H    |
        #|         | 'C8' | BEH8 |  INT8  |  BEH8H | INT8H    |
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')  # MATHE, MATHP
        #: The layer direction for BEHi=COMPS or AXCOMP. (Integer=1/2; Default=1)
        direct = integer_or_blank(card, 3, 'direct', default=1)
        thickness = double_or_blank(card, 4, 'thickness', default=1.0)

        #: Analysis type.
        #  'IS' - Implicit structural elements are being referred to.
        #  'IH' - Implicit heat analysis elements are being referred to.
        #  'ISH' - Implicit structural and heat elements are being referred to. (Character Default ISH)
        analysis = string_or_blank(card, 5, 'analysis', default='ISH')
        i = 9
        beh3 = None
        int3 = None
        beh4 = None
        int4 = None
        beh6 = None
        int6 = None
        beh8 = None
        int8 = None

        beh3h = None
        int3h = None
        beh4h = None
        int4h = None
        beh6h = None
        int6h = None
        beh8h = None
        int8h = None
        while i < len(card):
            code = string(card, i, 'code')
            # C3 : applies to elements with three corner grids
            # C4 : applies to elements with four corner grids
            # C6 : applies to elements with three corner grids and three midside grids
            # C8 : applies to elements with four corner grids and four midside grids
            # BEHi Element structural behavior. See Remark 7. (Default: PLSTRN for BEH3, BEH4, BEH6, and BEH8)
            # INTi Integration scheme. See Remarks 7. and 10. (Default: L for INT3, INT4, Q for INT6 and INT8)
            # BEHiH Element heat behavior. (Default: PLSTRN for BEH3H, BEH4H, BEH6H, and BEH8H.)
            # INTiH Integration scheme.    (Default: L for INT3H, L for INT4H, Q for INT6H and INT8H)
            if code == 'C3':
                beh3 = string_or_blank(card, i+1, 'BEH3', default='PLSTRN')
                int3 = string_or_blank(card, i+2, 'INT3', default='L')
                beh3h = string_or_blank(card, i+3, 'BEH3H', default='PLSTRN')
                int3h = string_or_blank(card, i+4, 'INT3H', default='L')
            elif code == 'C4':
                beh4 = string_or_blank(card, i+1, 'BEH4', default='PLSTRN')
                int4 = string_or_blank(card, i+2, 'INT4', default='L')
                beh4h = string_or_blank(card, i+3, 'BEH4H', default='PLSTRN')
                int4h = string_or_blank(card, i+4, 'INT4H', default='L')
            elif code == 'C6':
                beh6 = string_or_blank(card, i+1, 'BEH6', default='PLSTRN')
                int6 = string_or_blank(card, i+2, 'INT6', default='Q')
                beh6h = string_or_blank(card, i+3, 'BEH6H', default='PLSTRN')
                int6h = string_or_blank(card, i+4, 'INT6H', default='Q')
            elif code == 'C8':
                beh8 = string_or_blank(card, i+1, 'BEH8', default='PLSTRN')
                int8 = string_or_blank(card, i+2, 'INT8', default='Q')
                beh8h = string_or_blank(card, i+3, 'BEH8H', default='PLSTRN')
                int8h = string_or_blank(card, i+4, 'INT8H', default='Q')
            else:
                raise NotImplementedError(f'PSHLN2 code={code!r}')
            i += 8
        behx = [beh3, beh4, beh6, beh8]
        behxh = [beh3h, beh4h, beh6h, beh8h]
        integration = [int3, int4, int6, int8]
        integration_h = [int3h, int4h, int6h, int8h]
        assert len(card) <= 38, f'len(PSHLN2 card) = {len(card):d}\ncard={card}'
        self.cards.append((pid, mid, direct, thickness, analysis,
                           behx, integration, behxh, integration_h, comment))
        self.n += 1
        return self.n - 1

    @Property.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        material_id = np.zeros(ncards, dtype='int32')
        direct = np.zeros(ncards, dtype='int32')
        analysis = np.zeros(ncards, dtype='|U8')
        thickness = np.zeros(ncards, dtype='float64')
        beh = np.zeros((ncards, 4), dtype='|U8')
        beh_h = np.zeros((ncards, 4), dtype='|U8')
        integration = np.zeros((ncards, 4), dtype='|U8')
        integration_h = np.zeros((ncards, 4), dtype='|U8')

        for icard, card in enumerate(self.cards):
            (pid, mid, directi, thicknessi, analysisi,
             behx, integrationi, behxh, integration_hi, comment) = card
            if behx is None:
                behx = [''] * 4
            if integrationi is None:
                integrationi = [''] * 4
            if behxh is None:
                behxh = [''] * 4
            if integration_hi is None:
                integration_hi = [''] * 4

            property_id[icard] = pid
            material_id[icard] = mid
            thickness[icard] = thicknessi
            direct[icard] = directi
            analysis[icard] = analysisi
            beh[icard] = [val if val is not None else ''
                          for val in behx]
            integration[icard] = [val if val is not None else ''
                                  for val in integrationi]
            beh_h[icard] = [val if val is not None else ''
                            for val in behxh]
            integration_h[icard] = [val if val is not None else ''
                                    for val in integration_hi]
        self._save(property_id, material_id,
                   thickness, direct, analysis,
                   beh, integration,
                   beh_h, integration_h)
        self.sort()
        #self.model.log.warning(f'PSHLN2 self.thickness={self.thickness}')
        self.cards = []

    def _save(self, property_id, material_id,
              thickness, direct, analysis,
              beh, integration,
              beh_h, integration_h,
              ifile=None, comment: Optional[dict[int, str]]=None) -> None:
        ncards = len(property_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if not len(self.property_id) == 0:
            raise RuntimeError(f'stacking of {self.type} is not supported')

        save_ifile_comment(self, ifile, comment)
        self.property_id = property_id
        self.material_id = material_id
        self.thickness = thickness
        self.direct = direct
        self.analysis = analysis

        self.beh = beh
        self.integration = integration
        self.beh_h = beh_h
        self.integration_h = integration_h

    def set_used(self, used_dict: [str, list[np.ndarray]]) -> None:
        used_dict['material_id'].append(self.material_id)
        #used_dict['coord_id'].append(self.coord_id)

    def __apply_slice__(self, prop: PLPLANE, i: np.ndarray) -> None:
        prop.n = len(i)
        self._slice_comment(prop, i)
        prop.ifile = self.ifile[i]
        prop.property_id = self.property_id[i]
        prop.material_id = self.material_id[i]
        prop.thickness = self.thickness[i]
        prop.direct = self.direct[i]
        prop.analysis = self.analysis[i]

        prop.beh = self.beh[i, :]
        prop.integration = self.integration[i, :]
        prop.beh_h = self.beh_h[i, :]
        prop.integration_h = self.integration_h[i, :]

    @property
    def max_id(self) -> int:
        return max(self.property_id.max(), self.material_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        property_ids = array_str(self.property_id, size=size)
        material_ids = array_str(self.material_id, size=size)
        directs = array_default_int(self.direct, default=1, size=size)
        #inan = np.isnan(self.thickness)

        #self.property_id = np.zeros(ncards, dtype='int32')
        #self.material_id = np.zeros(ncards, dtype='int32')
        #self.direct = np.zeros(ncards, dtype='int32')
        #self.analysis = np.zeros(ncards, dtype='|U8')
        #self.thickness = np.zeros(ncards, dtype='float64')
        #self.beh = np.zeros((ncards, 2), dtype='|U8')
        #self.beh_h = np.zeros((ncards, 2), dtype='|U8')
        #self.integration = np.zeros((ncards, 2), dtype='|U8')
        #self.integration_h = np.zeros((ncards, 2), dtype='|U8')

        codes = ['C3', 'C4', 'C6', 'C8']
        for pid, mid, direct, analysis, thickness, beh, integration, beh_h, integration_h in zip_longest(
            property_ids, material_ids, directs, self.analysis, self.thickness,
            self.beh, self.integration, self.beh_h, self.integration_h):
            list_fields = ['PSHLN2', pid, mid, direct, thickness, analysis, None, None, None]
            #values = (beh, integration, beh_h, integration_h)
            for code, behx, intx, behxh, intxh in zip(codes, beh, integration, beh_h, integration_h):
                if behx == '' and intx == '' and behx == '' and intxh == '':
                    continue
                list_fields.extend([code, behx, intx, behxh, intxh, '', '', ''])
            bdf_file.write(print_card(list_fields))
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

    #def total_thickness(self) -> np.ndarray:
        #return self.thickness


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
