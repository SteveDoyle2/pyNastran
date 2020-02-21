"""
All shell properties are defined in this file.  This includes:

 * PCOMP
 * PCOMPG
 * PLPLANE
 * PSHEAR
 * PSHELL
 * PPLANE

All shell properties are Property objects.

"""
from __future__ import annotations
import copy
from itertools import count
import warnings
from typing import List, Optional, Union, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Property, Material
from pyNastran.bdf.cards.optimization import break_word_by_trailing_integer
from pyNastran.bdf.cards.materials import get_mat_props_S
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, BDFCard, MAT1, MAT8, MAT9


class CompositeShellProperty(Property):
    """
    Common class for:
     - PCOMP
     - PCOMPG
    """
    def __init__(self):
        Property.__init__(self)
        self.mids = []
        self.thicknesses = []
        self.thetas = []
        self.souts = []
        self.z0 = 0.
        self.nsm = 0.
        self.tref = 0.
        self.ge = 0.
        self.sb = 0.
        self.ft = None
        self.lam = None
        self.mids_ref = None

    @property
    def TRef(self) -> float:
        return self.tref
    @TRef.setter
    def TRef(self, tref: float) -> None:
        self.tref = tref

    def MassPerArea(self, iply='all', method='nplies', tflag: int=1, tscales=None):
        return self.get_mass_per_area(iply, method)

    def MassPerArea_structure(self) -> float:
        rhos = [mat_ref.get_density() for mat_ref in self.mids_ref]
        return self.get_mass_per_area_structure(rhos)

    def Thickness(self, iply: Union[int, str]='all', tflag: int=1, tscales=None):
        return self.get_thickness(iply)

    def Nsm(self) -> float:
        return self.get_nonstructural_mass()

    def Rho(self, iply: int) -> float:
        return self.get_density(iply)

    def Theta(self, iply: int) -> float:
        return self.get_theta(iply)

    def sout(self, iply: int) -> str:
        return self.get_sout(iply)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        mids_ref = []
        for iply in range(len(self.thicknesses)):
            mid = self.mids[iply]
            msg = ', which is required by %s pid=%s iply=%s' % (self.type, self.pid, iply)
            mids_ref.append(model.Material(mid, msg))
        self.mids_ref = mids_ref

        z1 = self.z0
        t = self.Thickness()
        try:
            z2 = z1 + t
        except TypeError:
            msg = 'Type=%s z1=%s t=%s; z2=z1+t is undefined' % (self.type, z1, t)
            raise TypeError(msg)
        if not ((-1.5*t <= z1 <= 1.5*t) or (-1.5*t <= z2 <= 1.5*t)):
            msg = '%s pid=%s midsurface: z1=%s z2=%s t=%s not in range of -1.5t < zi < 1.5t' % (
                self.type, self.pid, z1, z2, t)
            model.log.warning(msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mids_ref = None

    @property
    def is_symmetrical(self) -> bool:
        """
        Is the laminate symmetrical?

        Returns
        -------
        is_symmetrical : bool
            is the SYM flag active?

        """
        if self.lam == 'SYM':
            return True
        return False

    def _adjust_ply_id(self, iply: int) -> int:
        """
        Gets the ply ID that's stored in **self.plies**.

        When a ply is not symmetric, this function returns the input iply.
        When a ply is symmetrical and the iply value is greater than the
        number of plies, we return the mirrored ply.  For the case of a
        symmetrical ply, the element will always have an even number of
        layers.

        Parameters
        ----------
        iply : int
            the ply ID

        Raises
        ------
         - IndexError if iply is invalid

        ::

          Case 1 (nplies=6, len(plies)=3, lam='SYM'):
              ply 2
              ply 1
              ply 0
              ------- sym
              ply 0 / 3
              ply 1 / 4
              ply 2 / 5
            Ask for ply 3, return ply 0
            Ask for ply 4, return ply 1
            Ask for ply 5, return ply 2

          Case 2 (nplies=5, len(plies)=5, lam='NO'):
              ply 5
              ply 4
              ply 3
              ply 1
              ply 0
            Ask for ply 3, return ply 1
            Ask for ply 4, return ply 2

        """
        if iply == 'all':
            return iply

        # iply = 5
        # nplies = 3 (sym=6)
        # iply2 = 'iply' or 'iply - nplies'
        # iply2 = 0 = 9
        # iply2 = 1 = 1
        # iply2 = 2 = 2
        # iply2 = 3 - 3 = 2
        # iply2 = 4 - 3 = 1
        # iply2 = 5 - 3 = 0
        nplies = len(self.thicknesses)
        if iply >= nplies:
            if iply < self.nplies:
                iply = nplies - iply - 1
            else:
                raise IndexError('invalid value for nplies=%s iply=%r (iply is 0-based)\n%s' % (
                    nplies, iply, str(self)))
        elif iply < 0:
            raise IndexError('invalid value for nplies=%s iply=%r (iply is 0-based)\n%s' % (
                nplies, iply, str(self)))
        return iply

    def get_material_id(self, iply: int) -> int:
        iply = self._adjust_ply_id(iply)
        mid = self.mids[iply]
        return mid

    def get_thickness(self, iply: Union[int, str]='all') -> float:
        """
        Gets the thickness of the :math:`i^{th}` ply.

        Parameters
        ----------
        iply : int/str; default='all'
            the string **'all'** (default) or the mass per area of
            the :math:`i^{th}` ply

        Returns
        -------
        thickness : float
            the thickness of the ply or plies

        """
        #nplies = len(self.thicknesses)
        if iply == 'all':  # get all layers
            thick = sum(self.thicknesses)
            if self.is_symmetrical:
                return thick * 2.
            return thick
        else:
            iply = self._adjust_ply_id(iply)
            thick = self.thicknesses[iply]
            return thick

    @property
    def nplies(self) -> int:
        r"""
        Gets the number of plies including the core.

            ::

              if Lam=SYM:
                returns nplies * 2   (even)
              else:
                returns nplies

            """
        nplies = len(self.thicknesses)
        if self.is_symmetrical:
            return nplies * 2
        return nplies

    def get_nonstructural_mass(self) -> float:
        """Gets the non-structural mass :math:`i^{th}` ply"""
        return self.nsm

    def Mids(self) -> List[int]:
        return self.material_ids

    @property
    def material_ids(self) -> List[int]:
        """
        Gets the material IDs of all the plies

        Returns
        -------
        mids : MATx
            the material IDs

        """
        mids = []
        for iply in range(self.nplies):
            mids.append(self.Mid(iply))
        return mids

    def get_density(self, iply) -> float:
        """
        Gets the density of the :math:`i^{th}` ply

        Parameters
        ----------
        iply : int
            the ply ID (starts from 0)

        """
        iply = self._adjust_ply_id(iply)
        mid_ref = self.mids_ref[iply]
        return mid_ref.rho

    def Mid(self, iply) -> int:
        """
        Gets the Material ID of the :math:`i^{th}` ply.

        Parameters
        ----------
        iply : int/str; default='all'
            the string **'all'** (default) or the mass per area of
            the :math:`i^{th}` ply

        Returns
        -------
        material_id : int
            the material id of the ith ply

        """
        iply = self._adjust_ply_id(iply)
        if self.mids_ref is not None:
            mid = self.mids_ref[iply].mid
        else:
            mid = self.mids[iply]
        return mid

    def Material(self, iply) -> Union[MAT1, MAT8, MAT9]:
        """
        Gets the material of the :math:`i^{th}` ply (not the ID unless
        it is not cross-referenced).

        Parameters
        ----------
        iply : int
            the ply ID (starts from 0)

        """
        iply = self._adjust_ply_id(iply)
        if self.mids_ref is not None:
            mid = self.mids_ref[iply]
        else:
            mid = self.mids[iply]
        return mid

    def get_theta(self, iply) -> float:
        """
        Gets the ply angle of the :math:`i^{th}` ply (not the ID)

        Parameters
        ----------
        iply : int
            the ply ID (starts from 0)

        """
        iply = self._adjust_ply_id(iply)
        theta = self.thetas[iply]
        return theta

    def get_sout(self, iply) -> str:
        """
        Gets the the flag identifying stress/strain outpur of the
        :math:`i^{th}` ply (not the ID).  default='NO'.

        Parameters
        ----------
        iply : int
            the ply ID (starts from 0)

        """
        iply = self._adjust_ply_id(iply)
        sout = self.souts[iply]
        return sout

    def get_material_ids(self):
        thickness = []
        for i in range(self.nplies):
            thick = self.get_material_id(i)
            thickness.append(thick)
        return np.array(thickness, dtype='int32')

    def get_thicknesses(self) -> np.ndarray:
        thickness = []
        for i in range(self.nplies):
            thick = self.get_thickness(i)
            thickness.append(thick)
        return np.array(thickness, dtype='float64')

    def get_thetas(self) -> np.ndarray:
        thetas = []
        for i in range(self.nplies):
            theta = self.get_theta(i)
            thetas.append(theta)
        return np.array(thetas, dtype='float64')

    def get_z_locations(self) -> np.ndarray:
        """
        Gets the z locations for the various plies.

        Parameters
        ----------
        iply : int
            the ply ID (starts from 0)

        Assume there are 2 plies, each of 1.0 thick, starting from :math:`z=0`.

        >>> pcomp.get_z_locations()
        [0., 1., 2.]

        """
        zi = self.z0
        z = [zi]
        for thick in self.get_thicknesses():
            zi += thick
            z.append(zi)
        return np.array(z)

    def get_mass_per_area(self, iply='all', method: str='nplies') -> float:
        r"""
        Gets the Mass/Area for the property.

        .. math:: \frac{m}{A} = \sum(\rho t) + nsm

        or

        .. math:: \frac{m}{A} - nsm = \sum(\rho t)

        and

        .. math:: \frac{m_i}{A} = rho_i t_i + nsm_i

        where :math:`nsm_i` is the non-structural mass of the
        :math:`i^{th}` ply

        Parameters
        ----------
        iply : str/int; default='all'
            str : the string 'all' (default) or the mass per area of
            the :math:`i^{th}` ply
        method : str
            the method to compute MassPerArea
            {nplies, rho*t, t}

           * **Case 1 (iply = all)**

             method has no effect because the total nsm is defined

           * **Case 2 (iply != all)**

             method **'nplies'** smear the nsm based on :math:`n_{plies}` (default)

             :math:`nsm_i = nsm / n_{plies}`  # smear based on nplies

           * **Case 3 (iply != all)**

             method **'rho*t'** smear the nsm based on the mass distribution

             .. math:: nsm_i = \rho_i t_i \frac{nsm}{\sum(\rho_i t_i)}

             .. math:: nsm_i = \rho_i t_i \frac{nsm}{\frac{m}{A} - nsm}

           * **Case 4 (iply != all)**

             method **'t'** smear the nsm based on the thickness distribution

             .. math:: nsm_i = t_i \frac{nsm}{\sum(t_i)}

        .. note:: final mass calculation will be done later

        """
        rhos = [mat_ref.get_density() for mat_ref in self.mids_ref]
        return self.get_mass_per_area_rho(rhos, iply, method)

    def get_mass_per_area_rho(self, rhos: List[float],
                              iply='all', method: str='nplies') -> float:
        r"""
        Gets the Mass/Area for the property.

        .. math:: \frac{m}{A} = \sum(\rho t) + nsm

        or

        .. math:: \frac{m}{A} - nsm = \sum(\rho t)

        and

        .. math:: \frac{m_i}{A} = rho_i t_i + nsm_i

        where :math:`nsm_i` is the non-structural mass of the
        :math:`i^{th}` ply

        Parameters
        ----------
        rhos : List[float]
            the densities of each ply
        iply : str/int; default='all'
            the mass per area of the :math:`i^{th}` ply
        method : str
            the method to compute MassPerArea
            {nplies, rho*t, t}

           * **Case 1 (iply = all)**

             method has no effect because the total nsm is defined

           * **Case 2 (iply != all)**

             method **'nplies'** smear the nsm based on :math:`n_{plies}` (default)

             :math:`nsm_i = nsm / n_{plies}`  # smear based on nplies

           * **Case 3 (iply != all)**

             method **'rho*t'** smear the nsm based on the mass distribution

             .. math:: nsm_i = \rho_i t_i \frac{nsm}{\sum(\rho_i t_i)}

             .. math:: nsm_i = \rho_i t_i \frac{nsm}{\frac{m}{A} - nsm}

           * **Case 4 (iply != all)**

             method **'t'** smear the nsm based on the thickness distribution

             .. math:: nsm_i = t_i \frac{nsm}{\sum(t_i)}

        .. note:: final mass calculation will be done later

        """
        assert method in ['nplies', 'rho*t', 't'], 'method=%r is invalid' % method
        nplies = len(self.thicknesses)
        iply = self._adjust_ply_id(iply)
        if iply == 'all':  # get all layers
            #mass_per_area_total = m/A = sum(rho*t) + nsm
            #mass_per_area_total = mpa-nsm = sum(rho*t)
            #(m/A)i = rho*t + nsmi
            # where nsmi has two methods
            mass_per_area = 0.
            #nplies = len(self.thicknesses)
            for iplyi in range(nplies):
                #rho = self.get_density(iplyi)
                rho = rhos[iplyi]
                t = self.thicknesses[iplyi]
                mass_per_area += rho * t

            if self.is_symmetrical:
                return 2. * mass_per_area + self.nsm
            return mass_per_area + self.nsm
        else:
            assert isinstance(iply, integer_types), 'iply must be an integer; iply=%r' % iply
            #rho = self.get_density(iply)
            rho = rhos[iply]
            t = self.thicknesses[iply]

            if method == 'nplies':
                # we divide by nplies b/c it's nsm per area and
                # we're working on a per ply basis
                # nsmi = nsmi/n  # smear based on nplies
                mass_per_area = rho * t + self.nsm / self.nplies
            elif method == 'rho*t':
                # assume you smear the nsm mass based on rho*t distribution
                #nsmi = rho*t / sum(rho*t) * nsm
                #rho*t + nsmi = rho*t + rho*t/(sum(rho*t) + nsm - nsm) * nsm
                #rho*t + nsmi = rho*t + rho*t/(mass_per_area_total - nsm) * nsm
                #             = rho*t * (1 + nsm/(mass_per_area_total-nsm))
                mass_per_area_total = self.get_mass_per_area_rho(rhos, iply='all', method='nplies')
                mass_per_area = rho * t * (1.0 + self.nsm / (mass_per_area_total - self.nsm))
            elif method == 't':
                # assume you smear the nsm mass based on t distribution
                #nsmi = t / sum(t) * nsm
                #rho*t + nsmi = rho*t + t/sum(t) * nsm
                #rho*t + nsmi = rho*t + t/thickness_total * nsm
                #             = t * (rho + nsm/thickness_total)
                thickness_total = self.get_thickness()
                mass_per_area = t * (rho + self.nsm / thickness_total)
            else:
                raise NotImplementedError('method=%r is not supported' % method)
            return mass_per_area

    def get_mass_per_area_structure(self, rhos: List[float]) -> float:
        r"""
        Gets the Mass/Area for the property structure only
        (doesn't consider nsm).

        .. math:: \frac{m}{A} = \sum(\rho t)

        Parameters
        ----------
        rhos : List[float]
            the densities of each ply

        """
        nplies = len(self.thicknesses)
        mass_per_area = 0.
        for iply in range(nplies):
            rho = rhos[iply]
            t = self.thicknesses[iply]
            mass_per_area += rho * t

        ksym = 2. if self.is_symmetrical else 1.
        return ksym * mass_per_area


class PCOMP(CompositeShellProperty):
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
    type = 'PCOMP'
    _field_map = {
        1: 'pid', 2: 'z0', 3:'nsm', 4:'sb', 5:'ft', 6:'tref', 7: 'ge', 8:'lam',
    }
    _properties = ['_field_map', 'plies', 'nplies', 'material_ids']
    def update_by_pname_fid(self, pname_fid: Union[str, int],
                            value: Union[int, float, str]) -> None:
        if isinstance(pname_fid, int):
            self._update_field_helper(pname_fid, value)
        elif pname_fid == 'Z0':
            self.z0 = value
        elif pname_fid == 'SB':
            self.sb = value
        elif pname_fid == 'TREF':
            self.tref = value
        elif pname_fid == 'GE':
            self.ge = value
        elif pname_fid.startswith(('T', 'THETA')):
            word, num = break_word_by_trailing_integer(pname_fid)
            num = int(num)
            if word == 'T':
                self.thicknesses[num - 1] = value
            elif word == 'THETA':
                self.thetas[num - 1] = value
            else:
                raise RuntimeError('pid=%s pname_fid=%r word=%s\n' % (self.pid, pname_fid, word))
        else:
            raise NotImplementedError('property_type=%r has not implemented %r in pname_map' % (
                self.type, pname_fid))

    def _update_field_helper(self, n: int, value) -> None:
        if n == 3:
            self.z0 = value
            return
        elif n == 4:
            self.nsm = value
            return
        elif n == 5:
            self.sb = value
            return
        elif n == 7:
            self.tref = value
            return
        elif n == 8:
            self.ge = value
            return

        assert n > 0, 'PCOMP pid=%s; negative indicies are not supported (pname_fid=%r)' % (self.pid, n)
        nnew = n - 10
        if nnew <= 0:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

        # the + 1 is to tell us we're on row 1 at minimum (not row 0)
        irow = n // 10
        irow_start = irow * 10 + 2
        offset = n - irow_start
        if offset < 0:
            raise RuntimeError('field=%s is invalid for the PCOMP' % n)

        irow_layer = irow - 1
        ilayer = irow_layer * 2 + offset // 4
        #print('\nn=%s nnew=%s irow_start=%s irow=%s offset=%s ilayer=%s' % (
            #n, nnew, irow_start, irow, offset, ilayer))
        assert ilayer >= 0
        iply = self._adjust_ply_id(ilayer)

        slot = offset % 4
        #print('iply=%s slot=%s' % (iply, slot))
        try:
            ply = self.plies[iply]
        except IndexError:
            msg = 'PCOMP pid=%s; n=%s nnew=%s ilayer=%s slot=%r\n' % (
                self.pid, n, nnew, ilayer, slot)
            msg += ('On PCOMP pid=%r, ply %i is not defined.  '
                    'iply_min=0; iply_max=%i' % (self.pid, ilayer, len(self.plies)))
            raise IndexError(msg)

        if slot == 0:
            self.material_ids[iply] = value
        elif slot == 1:
            self.thicknesses[iply] = value
        elif slot == 2:
            self.thetas[iply] = value
        elif slot == 3:
            self.souts[iply] = value
        else:
            raise RuntimeError('ilayer=%s iply=%s slot=%s in [mid, t, theta, sout]' % (
                ilayer, iply, slot))
        # ply = [mid, t, theta, sout]
        ply[slot] = value

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mids = [1]
        thicknesses = [1.]
        return PCOMP(pid, mids, thicknesses,
                     thetas=None, souts=None, nsm=0., sb=0., ft=None,
                     tref=0., ge=0., lam=None, z0=None, comment='')

    def __init__(self, pid: int,
                 mids: List[int], thicknesses: List[float],
                 thetas: Optional[List[float]]=None,
                 souts: Optional[List[str]]=None,
                 nsm: float=0., sb: float=0.,
                 ft: Optional[str]=None,
                 tref: float=0., ge: float=0.,
                 lam: Optional[str]=None, z0: Optional[float]=None,
                 comment: str=''):
        """
        Creates a PCOMP card

        pid : int
            property id
        mids : List[int, ..., int]
            material ids for each ply
        thicknesses : List[float, ..., float]
            thicknesses for each ply
        thetas : List[float, ..., float]; default=None
            ply angle
            None : [0.] * nplies
        souts : List[str, ..., str]; default=None
            should the stress? be printed; {YES, NO}
            None : [NO] * nplies
        nsm : float; default=0.
            nonstructural mass per unit area
        sb : float; default=0.
            Allowable shear stress of the bonding material.
            Used by the failure theory
        ft : str; default=None
            failure theory; {HILL, HOFF, TSAI, STRN, None}
        tref : float; default=0.
            reference temperature
        ge : float; default=0.
            structural damping
        lam : str; default=None
            symmetric flag; {SYM, MEM, BEND, SMEAR, SMCORE, None}
            None : not symmmetric
        z0 : float; default=None
            Distance from the reference plane to the bottom surface
            None : -1/2 * total_thickness
        comment : str; default=''
            a comment for the card

        """
        CompositeShellProperty.__init__(self)
        if comment:
            self.comment = comment

        nplies = len(mids)
        if thetas is None:
            thetas = [0.] * nplies
        if souts is None:
            souts = ['NO'] * nplies
        #: Property ID
        self.pid = pid

        #: Non-Structural Mass per unit Area
        self.nsm = nsm

        # Allowable shear stress of the bonding material
        # required if there is a failure theory
        self.sb = sb

        #: Failure Theory
        #:   ['HILL', 'HOFF', 'TSAI', 'STRN', None]
        self.ft = ft

        #: Reference Temperature (default=0.0)
        self.tref = tref
        self.ge = ge

        #: symmetric flag - default = No Symmetry (=None)
        self.lam = lam
        self.mids = mids
        self.thicknesses = thicknesses
        self.thetas = thetas
        self.souts = souts
        if z0 is None:
            z0 = -0.5 * self.Thickness()
        self.z0 = z0

    def validate(self):
        assert self.ft in ['HILL', 'HOFF', 'TSAI', 'STRN', 0.0, None], 'ft=%r' % self.ft

        # 'NO' is not an option!
        allowed_lam = [None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE']
        if self.lam not in allowed_lam:
            msg = 'lam=%r is invalid; allowed=[%s]' % (self.lam, ', '.join(allowed_lam))
            raise ValueError(msg)

        for iply, mid, thickness, theta, sout in  zip(count(), self.material_ids, self.thicknesses,
                                                      self.thetas, self.souts):
            if thickness <= 0.:
                ply = [mid, thickness, theta, sout]
                msg = ('thickness of PCOMP layer is invalid pid=%s'
                       ' iLayer=%s t=%s ply=[mid,t,theta,'
                       'sout]=%s' % (self.pid, iply, thickness, ply))
                warnings.warn(msg)
                continue
                #raise RuntimeError(msg)

        # this is a loose requirement
        #if self.ft in ['HILL', 'HOFF', 'TSAI', 'STRN'] and self.sb <= 0.:
            #msg = 'PCOMP pid=%s FT=%s sb=%s; sb must be greater than 0' % (
                #self.pid, self.ft, self.sb)
            #raise ValueError(msg)

    @classmethod
    def add_card(cls, card: BDFCard, comment=''):
        """
        Adds a PCOMP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')

        # z0 is field 2 and is calculated at the end because we need the
        # thickness first
        #self.z0 = double_or_blank(card, 1, 'pid')

        nsm = double_or_blank(card, 3, 'nsm', 0.0)
        sb = double_or_blank(card, 4, 'sb', 0.0)
        ft = string_or_blank(card, 5, 'ft')
        tref = double_or_blank(card, 6, 'tref', 0.0)
        ge = double_or_blank(card, 7, 'ge', 0.0)
        lam = string_or_blank(card, 8, 'lam') # default=blank -> nothing

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
        for i in range(9, 9 + nplies * 4, 4):
            actual = card.fields(i, i + 4)
            mid = integer_or_blank(card, i, 'mid', mid_last)
            t = double_or_blank(card, i + 1, 't', thick_last)
            theta = double_or_blank(card, i + 2, 'theta', 0.0)
            sout = string_or_blank(card, i + 3, 'sout', 'NO')

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
        return PCOMP(pid, mids, thicknesses, thetas, souts, nsm, sb, ft, tref, ge,
                     lam, z0, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PCOMP card from the OP2

        Parameters
        ----------
        data : List[varies]
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
        ft = data[4]
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
        for (mid, t, theta, sout) in zip(Mid, T, Theta, Sout):
            if sout == 0:
                sout = 'NO'
            elif sout == 1:
                sout = 'YES'
            #elif sout == 2:  #: .. todo:: what?!!
                #sout = 'YES'
            #elif sout == 3:  #: .. todo:: what?!!
                #sout = 'YES'
            else:
                raise RuntimeError('unsupported sout.  sout=%r and must be 0 or 1.'
                                   '\nPCOMP = %s' % (sout, data))
            mids.append(mid)
            thicknesses.append(t)
            thetas.append(theta)
            souts.append(sout)
        if ft == 0:
            ft = None
        elif ft == 1:
            ft = 'HILL'
        elif ft == 2:
            ft = 'HOFF'
        elif ft == 3:
            ft = 'TSAI'
        elif ft == 4:
            ft = 'STRN'
        else:
            raise RuntimeError('unsupported ft.  pid=%s ft=%r.'
                               '\nPCOMP = %s' % (pid, ft, data))
        return PCOMP(pid, mids, thicknesses, thetas, souts,
                     nsm, sb, ft, tref, ge, lam, z0, comment=comment)

    @property
    def plies(self):
        plies = []
        for mid, t, theta, sout in zip(self.material_ids, self.thicknesses,
                                       self.thetas, self.souts):
            plies.append([mid, t, theta, sout])
        return plies

    #@plies.setter
    #def plies(self, plies):
        #i = 0
        #for mid, t, theta, sout in zip(plies):
            #self.mids[i] = mid
            #self.thicknesses[i] = t
            #self.thetas[i] = theta
            #self.souts[i] = sout
            #i += 1

    def is_balanced_symmetric(self, debug=True):
        """assumes materials with different mids are unique"""
        # TODO: balanced is not done
        is_balanced = None
        is_symmetric_z = False
        is_symmetric_thetas = False
        thetad = self.get_thetas()
        if debug:  # pragma: no cover
            print(f'  theta (deg) = {thetad}; n={len(thetad)}')
            print(f'  nplies = {self.nplies}')

        z0, z1, zmeans = self.get_z0_z1_zmean()
        nz = len(zmeans)
        nhalf = nz // 2

        mids = self.get_material_ids()
        mlow = mids[:nhalf]
        mhigh = mids[-nhalf:][::-1]

        if self.is_symmetrical:
            is_symmetric_z = True
            #is_symmetric_materials = True
            is_symmetric_thetas = True
            if debug:  # pragma: no cover
                print('  symmetric flag')
        else:
            zlow = zmeans[:nhalf]
            zhigh = zmeans[-nhalf:][::-1]
            if debug:  # pragma: no cover
                print(f'zlow={zlow}; n={len(zlow)}')
                print(f'zhigh={zhigh}; n={len(zhigh)}')

            tlow = thetad[:nhalf]
            thigh = thetad[-nhalf:][::-1]

            # works for even and odd
            # if the z's are opposite and the thetas are the same
            # [-1, 0, 1] - symmetric z
            # [-1, 1] - symmetric z
            # [0, 45, 90, 90, 45, 0] - symmetric theta
            # [0, 45, 90, 45, 0] - symmetric theta
            if np.allclose(-zlow, zhigh):
                is_symmetric_z = True
            if np.allclose(tlow, thigh):
                is_symmetric_thetas = True

        # doesn't handle:
        #    mids = [2, 3]
        # where the materials are the same
        # (you need to look at the properties)

        if np.allclose(mlow, mhigh):
            is_symmetric_materials = True
        is_symmetric = (
            is_symmetric_z and
            is_symmetric_thetas and
            is_symmetric_materials
        )
        return is_balanced, is_symmetric

    def get_z0_z1_zmean(self):
        thicknesses = self.get_thicknesses()
        csum = np.cumsum(thicknesses)
        z0 = self.z0 + np.hstack([0., csum[:-1]])
        z1 = self.z0 + csum
        #dz = z1 - z0
        #dzsquared = z1 ** 2 - z0 ** 2
        #zcubed = z1 ** 3 - z0 ** 3
        zmeans = (z0 + z1) / 2.
        return z0, z1, zmeans

    def get_individual_ABD_matrices(
            self, theta_offset: float=0.) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Gets the ABD matrix

        Parameters
        ----------
        theta_offset : float
            rotates the ABD matrix; measured in degrees

        """
        mids = self.get_material_ids()
        thicknesses = self.get_thicknesses()
        thetad = self.get_thetas()
        theta = np.radians(thetad)
        assert len(mids) == len(thicknesses)
        assert len(mids) == len(thetad)
        z0, z1, zmeans = self.get_z0_z1_zmean()

        assert len(z0) == len(z1)
        # A11 A12 A16
        # A12 A22 A26
        # A16 A26 A66
        A = np.zeros((3, 3), dtype='float64')
        B = np.zeros((3, 3), dtype='float64')
        D = np.zeros((3, 3), dtype='float64')
        #Q = np.zeros((3, 3), dtype='float64')

        mids_ref = copy.deepcopy(self.mids_ref)
        assert mids_ref is not None, f'The following material hasnt been cross-referenced:\n{self}'
        if self.is_symmetrical:
            mids_ref += mids_ref[::-1]

        assert len(mids) == len(mids_ref), 'mids=%s (%s) mids_ref:\n%s; %s' % (mids, len(mids), mids_ref, len(mids_ref))
        for mid, mid_ref, thetai, thickness, zmean, z0i, z1i in zip(mids, mids_ref, theta,
                                                                    thicknesses, zmeans, z0, z1):
            Qbar = self.get_Qbar_matrix(mid_ref, thetai)
            A += Qbar * thickness
            #B += Qbar * thickness * zmean
            #D += Qbar * thickness * (z1i ** 3 - z0i ** 3)
            B += Qbar * (z1i ** 2 - z0i ** 2)
            D += Qbar * (z1i ** 3 - z0i ** 3)
            #N += Q * alpha * thickness
            #M += Q * alpha * thickness * zmean
        B /= 2.
        D /= 3.
        #M /= 2.
        return A, B, D

    def get_ABD_matrices(self, theta_offset: float=0.) -> np.ndarray:
        """
        Gets the ABD matrix

        Parameters
        ----------
        theta_offset : float
            rotates the ABD matrix; measured in degrees

        """
        A, B, D = self.get_individual_ABD_matrices(theta_offset)
        ABD = np.block([
            [A, B],
            [B, D],
        ])
        #np.set_printoptions(linewidth=120, suppress=True)
        #print(ABD)
        return ABD

    def get_Qbar_matrix(self, mid_ref, theta: float) -> np.ndarray:
        """theta must be in radians"""
        S2, unused_S3 = get_mat_props_S(mid_ref)
        T = get_2d_plate_transform(theta)

        #R = np.array([
            #[1., 0., 0.],
            #[0., 1., 0.],
            #[0., 0., 2.],
        #])
        Tinv = np.linalg.inv(T)
        Q = np.linalg.inv(S2)

        # [Qbar] = [T^-1][Q][T^-T]
        # [T^-T] = [R][T][R^-1] = [T^-1].T
        Qbar = np.linalg.multi_dot([Tinv, Q, Tinv.T])

        # [T.T] = [R][T^-1][R^-1]
        #Sbar = np.linalg.multi_dot([T.T, S2, T])
        #Qbar = np.linalg.inv(Sbar)
        return Qbar

    def get_Sbar_matrix(self, mid_ref, theta: float) -> np.ndarray:
        """theta must be in radians"""
        Qbar = self.get_Qbar_matrix(mid_ref, theta)
        Sbar = np.linalg.inv(Qbar)
        return Sbar

    def _verify(self, xref: bool) -> None:
        pid = self.Pid()
        is_sym = self.is_symmetrical
        nplies = self.nplies
        nsm = self.Nsm()
        mids = self.Mids()

        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(is_sym, bool), 'is_sym=%r' % is_sym
        assert isinstance(nplies, integer_types), 'nplies=%r' % nplies
        assert isinstance(nsm, float), 'nsm=%r' % nsm
        assert isinstance(mids, list), 'mids=%r' % mids

        t = self.Thickness()
        assert isinstance(t, float), 'thickness=%r' % t
        if xref:
            mpa = self.MassPerArea()
            assert isinstance(mpa, float), 'mass_per_area=%r' % mpa

        for iply in range(nplies):
            mid2 = self.Mid(iply)
            assert mids[iply] == mid2
            t = self.Thickness(iply)
            assert isinstance(t, float), 'thickness=%r' % t

            if xref:
                rho = self.Rho(iply)
                mpa = self.MassPerArea(iply)
                assert isinstance(rho, float), 'rho=%r' % rho
                assert isinstance(mpa, float), 'mass_per_area=%r' % mpa
                iplyi = self._adjust_ply_id(iply)
                mid = self.mids_ref[iplyi]
                assert mid.type in ['MAT1', 'MAT2', 'MAT8'], 'PCOMP: mid.type=%r' % mid.type

        for ply in self.plies:
            assert len(ply) == 4, ply

    def raw_fields(self):
        list_fields = ['PCOMP', self.pid, self.z0, self.nsm, self.sb, self.ft,
                       self.tref, self.ge, self.lam, ]
        for (mid, t, theta, sout) in zip(self.material_ids, self.thicknesses,
                                         self.thetas, self.souts):
            list_fields += [mid, t, theta, sout]
        return list_fields

    def repr_fields(self):
        nsm = set_blank_if_default(self.nsm, 0.0)
        sb = set_blank_if_default(self.sb, 0.0)
        tref = set_blank_if_default(self.tref, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        z0 = set_blank_if_default(self.z0, -0.5 * self.get_thickness())

        list_fields = ['PCOMP', self.pid, z0, nsm, sb, self.ft, tref, ge, self.lam]
        for (mid, t, theta, sout) in zip(self.material_ids, self.thicknesses,
                                         self.thetas, self.souts):
            #theta = set_blank_if_default(theta, 0.0)
            str_sout = set_blank_if_default(sout, 'NO')
            list_fields += [mid, t, theta, str_sout]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PCOMPG(CompositeShellProperty):
    """
    +--------+----------+------+-----+--------+--------+------+----+-----+
    |   1    |     2    |   3  |  4  |    5   |    6   |   7  | 8  |  9  |
    +========+==========+======+=====+========+========+======+====+=====+
    | PCOMPG |    PID   |  Z0  | NSM |   SB   |   FT   | TREF | GE | LAM |
    +--------+----------+------+-----+--------+--------+------+----+-----+
    |        |  GPLYID1 | MID1 |  T1 | THETA1 | SOUT1  |      |    |     |
    +--------+----------+------+-----+--------+--------+------+----+-----+
    |        |  GPLYID2 | MID2 |  T2 | THETA2 | SOUT2  |      |    |     |
    +--------+----------+------+-----+--------+--------+------+----+-----+
    """
    type = 'PCOMPG'
    _field_map = {
        1: 'pid', 2: 'z0', 3:'nsm', 4:'sb', 5:'ft', 6:'tref',
        7: 'ge', 8:'lam',
    }
    pname_fid_map = {
        # 1 based
        #4 : 't', 'T' : 't',
        #6 : 'twelveIt3', # no option
        #8 : 'tst', #'T' : 't',
    }
    _properties = ['_field_map', 'plies', 'nplies', 'material_ids']
    def update_by_pname_fid(self, pname_fid: Union[str, int], value) -> None:
        if isinstance(pname_fid, int):
            self._update_field_helper(pname_fid, value)
        else:
            raise NotImplementedError('PCOMPG: pname_fid=%r' % pname_fid)

    def _update_field_helper(self, n: int, value) -> None:
        nnew = n - 10
        if nnew <= 0:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

        ilayer = nnew // 10
        #print('ilayer=%s n=%s nnew=%s' % (ilayer, n, nnew))
        try:
            ply = self.plies[ilayer]
        except IndexError:
            msg = ('On PCOMPG pid=%r, ply %i is not defined.  '
                   'iply_min=0; iply_max=%i' % (self.pid, ilayer, len(self.plies)))
            raise IndexError(msg)

        #ply = [mid, thickness, theta, sout, global_ply_id]
        slot = nnew % 10
        #print('ply=', ply, nnew, slot)
        if slot == 4:
            self.thicknesses[ilayer] = value
        elif slot == 5:
            self.thetas[ilayer] = value
        else:
            raise NotImplementedError('On PCOMPG: cannot update n=%s; slot=%i\n'
                                      '  slot=2: Global Ply ID\n'
                                      '  slot=3: Material ID\n'
                                      '  slot=4: Thickness\n'
                                      '  slot=5: Theta\n'
                                      '  slot=6: SOUT\n' % (n, slot))
            #ply[slot] = value

        #for mid, t, theta, sout, global_ply_id in zip(self.mids, self.thicknesses,
                                                      #self.thetas, self.souts,
                                                      #self.global_ply_ids):

    @property
    def plies(self) -> List[Tuple[int, float, float, str, int]]:
        plies = []
        for mid, t, theta, sout, global_ply_id in zip(self.mids, self.thicknesses,
                                                      self.thetas, self.souts,
                                                      self.global_ply_ids):
            plies.append((mid, t, theta, sout, global_ply_id))
        return plies

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        global_ply_ids = [1]
        mids = [1]
        thicknesses = [1.]
        return PCOMPG(pid, global_ply_ids, mids, thicknesses,
                      thetas=None, souts=None, nsm=0.0, sb=0.0, ft=None,
                      tref=0.0, ge=0.0, lam=None, z0=None, comment='')

    def __init__(self, pid, global_ply_ids, mids, thicknesses, thetas=None, souts=None,
                 nsm=0.0, sb=0.0, ft=None, tref=0.0, ge=0.0, lam=None, z0=None, comment=''):
        """
        Creates a PCOMPG card

        Parameters
        ----------
        pid : int
            property id
        global_ply_ids : List[int]
            the ply id
        mids : List[int, ..., int]
            material ids for each ply
        thicknesses : List[float, ..., float]
            thicknesses for each ply
        thetas : List[float, ..., float]; default=None
            ply angle
            None : [0.] * nplies
        souts : List[str, ..., str]; default=None
            should the stress? be printed; {YES, NO}
            None : [NO] * nplies
        nsm : float; default=0.
            nonstructural mass per unit area
        sb : float; default=0.
            Allowable shear stress of the bonding material.
            Used by the failure theory
        ft : str; default=None
            failure theory; {HILL, HOFF, TSAI, STRN, None}
        tref : float; default=0.
            reference temperature
        ge : float; default=0.
            structural damping
        lam : str; default=None
            symmetric flag; {SYM, MEM, BEND, SMEAR, SMCORE, None}
            None : not symmmetric
        z0 : float; default=None
            Distance from the reference plane to the bottom surface
            None : -1/2 * total_thickness
        comment : str; default=''
            a comment for the card

        """
        CompositeShellProperty.__init__(self)
        if comment:
            self.comment = comment
        if thetas is None:
            nplies = len(mids)
            thetas = [0.] * nplies
        if souts is None:
            souts = ['NO'] * nplies

        #: Property ID
        self.pid = pid

        #: Non-Structural Mass per unit Area
        self.nsm = nsm
        self.sb = sb

        #: Failure Theory
        #:
        #:   ['HILL', 'HOFF', 'TSAI', 'STRN', None]
        self.ft = ft

        #: Reference Temperature (default=0.0)
        self.tref = tref
        self.ge = ge

        #: symmetric flag - default = No Symmetry (NO)
        self.lam = lam
        self.mids = mids
        self.thicknesses = thicknesses
        self.thetas = thetas
        self.global_ply_ids = global_ply_ids
        self.souts = souts
        if z0 is None:
            z0 = -0.5 * self.Thickness()
        self.z0 = z0

    def validate(self) -> None:
        assert isinstance(self.global_ply_ids, list), self.global_ply_ids
        sorted_global_ply_ids = sorted(self.global_ply_ids)
        if not np.array_equal(sorted_global_ply_ids, np.unique(self.global_ply_ids)):
            msg = 'PCOMPG pid=%s; global_ply_ids=%s must be unique' % (
                self.pid, self.global_ply_ids)
            raise ValueError(msg)

        #assert self.ft in ['HILL', 'HOFF', 'TSAI', 'STRN', 0.0, None], 'ft=%r' % self.ft # PCOMP
        assert self.ft in ['HILL', 'HOFF', 'TSAI', 'STRN', None], 'ft=%r' % self.ft

        # 'NO' is not an option!
        allowed_lam = [None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE']
        if self.lam not in allowed_lam:
            msg = 'lam=%r is invalid; allowed=[%s]' % (self.lam, ', '.join(
                str(lam) for lam in allowed_lam))
            raise ValueError(msg)
        for iply, sout in enumerate(self.souts):
            assert sout in ['YES', 'NO'], "iply=%s sout=%r; SOUT must be ['YES', 'NO']" % (iply, sout)

        # this is a loose requirement
        #if self.ft in ['HILL', 'HOFF', 'TSAI', 'STRN'] and self.sb <= 0.:
            #msg = 'PCOMP pid=%s FT=%s sb=%s; sb must be greater than 0' % (
                #self.pid, self.ft, self.sb)
            #raise ValueError(msg)

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a PCOMPG card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        # z0 will be calculated later
        nsm = double_or_blank(card, 3, 'nsm', 0.0)
        sb = double_or_blank(card, 4, 'sb', 0.0)
        ft = string_or_blank(card, 5, 'ft')
        tref = double_or_blank(card, 6, 'tref', 0.0)
        ge = double_or_blank(card, 7, 'ge', 0.0)
        lam = string_or_blank(card, 8, 'lam')

        fields = card.fields(9)

        #T = 0.  # thickness
        mid_last = None
        thick_last = None

        i = 0
        #n = 0
        mids = []
        thicknesses = []
        thetas = []
        souts = []
        global_ply_ids = []
        while i < len(fields):
            global_ply_id = integer(card, 9 + i, 'global_ply_id')
            mid = integer_or_blank(card, 9 + i + 1, 'mid', mid_last)

            # can be blank 2nd time thru
            thickness = double_or_blank(card, 9 + i + 2, 'thickness', thick_last)

            theta = double_or_blank(card, 9 + i + 3, 'theta', 0.0)
            sout = string_or_blank(card, 9 + i + 4, 'sout', 'NO')
            #print('n=%s global_ply_id=%s mid=%s thickness=%s len=%s' %(
            #    n,global_ply_id,mid,thickness,len(fields)))

            mids.append(mid)
            thicknesses.append(thickness)
            thetas.append(theta)
            souts.append(sout)
            global_ply_ids.append(global_ply_id)

            assert mid is not None
            assert thickness is not None
            assert isinstance(mid, integer_types), 'mid=%s' % mid
            assert isinstance(thickness, float), 'thickness=%s' % thickness
            mid_last = mid
            thick_last = thickness
            #T += thickness
            i += 8
            #n += 1
        z0 = double_or_blank(card, 2, 'z0')
        return PCOMPG(pid,
                      global_ply_ids, mids, thicknesses, thetas, souts,
                      nsm, sb, ft, tref, ge, lam, z0, comment=comment)

    def _verify(self, xref: bool) -> None:
        pid = self.Pid()
        is_sym = self.is_symmetrical
        nplies = self.nplies
        nsm = self.Nsm()
        mids = self.Mids()

        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(is_sym, bool), 'is_sym=%r' % is_sym
        assert isinstance(nplies, integer_types), 'nplies=%r' % nplies
        assert isinstance(nsm, float), 'nsm=%r' % nsm
        assert isinstance(mids, list), 'mids=%r' % mids

        t = self.Thickness()
        assert isinstance(t, float), 'thickness=%r' % t
        if xref:
            mpa = self.MassPerArea()
            assert isinstance(mpa, float), 'mass_per_area=%r' % mpa

        for iply in range(nplies):
            glply = self.GlobalPlyID(iply)
            mid2 = self.Mid(iply)
            assert mids[iply] == mid2
            t = self.Thickness(iply)
            assert isinstance(glply, integer_types), 'global_ply_id=%r' % glply
            assert isinstance(t, float), 'thickness=%r' % t
            if xref:
                rho = self.Rho(iply)
                mpa = self.MassPerArea(iply)
                assert isinstance(mpa, float), 'mass_per_area=%r' % mpa
                assert isinstance(rho, float), 'rho=%r' % rho
                mid = self.mids_ref[iply]
                assert mid.type in ['MAT1', 'MAT8'], 'PCOMPG: mid.type=%r' % mid.type

        #@plies.setter
        #def plies(self, plies):
            #i = 0
            #for mid, t, theta, sout, global_ply_id in zip(plies):
                #self.mids[i] = mid
                #self.thicknesses[i] = t
                #self.thetas[i] = theta
                #self.souts[i] = sout
                #self.global_ply_ids[i] = global_ply_id
                #i += 1

    def GlobalPlyID(self, iply: int) -> int:
        """returns the global ply id for the specified layer"""
        global_ply_id = self.global_ply_ids[iply]
        return global_ply_id

    def raw_fields(self):
        list_fields = [
            'PCOMPG', self.pid, self.z0, self.nsm, self.sb, self.ft,
            self.tref, self.ge, self.lam, ]
        for (mid, t, theta, sout, global_ply_id) in zip(
                self.material_ids, self.thicknesses, self.thetas, self.souts,
                self.global_ply_ids):
            list_fields += [global_ply_id, mid, t, theta, sout, None, None, None]
        return list_fields

    def repr_fields(self):
        nsm = set_blank_if_default(self.nsm, 0.0)
        sb = set_blank_if_default(self.sb, 0.0)
        tref = set_blank_if_default(self.tref, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        z0 = set_blank_if_default(self.z0, -0.5 * self.Thickness())

        list_fields = [
            'PCOMPG', self.pid, z0, nsm, sb, self.ft, tref, ge,
            self.lam]
        zipi = zip(self.material_ids, self.thicknesses, self.thetas, self.souts,
                   self.global_ply_ids)
        for (mid, t, theta, sout, global_ply_id) in zipi:
            sout = set_blank_if_default(sout, 'NO')
            list_fields += [global_ply_id, mid, t, theta, sout, None, None, None]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


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

    """
    type = 'PLPLANE'
    _field_map = {1: 'pid', 2:'mid', 6:'cid', 7:'str'}

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        return PLPLANE(pid, mid, cid=0, stress_strain_output_location='GRID', comment='')

    def __init__(self, pid, mid, cid=0, stress_strain_output_location='GRID', comment=''):
        """
        Creates a PLPLANE card, which defines the properties of a fully
        nonlinear (i.e., large strain and large rotation) hyperelastic
        plane strain or axisymmetric element.

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id; MATHP / MATHE
        cid : int; default=0
            ???
        stress_strain_output_location : str; default='GRID'
            ???
        comment : str; default=''
            a comment for the card

        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        #: Property ID
        self.pid = pid
        self.mid = mid
        self.cid = cid
        self.stress_strain_output_location = stress_strain_output_location
        self.mid_ref = None
        self.cid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PLPLANE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')  # MATHE, MATHP
        cid = integer_or_blank(card, 3, 'cid', 0)
        stress_strain_output_location = string_or_blank(card, 4, 'str', 'GRID')
        return PLPLANE(pid, mid, cid=cid,
                       stress_strain_output_location=stress_strain_output_location,
                       comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PLPLANE pid=%s' % self.pid
        self.mid_ref = model.HyperelasticMaterial(self.mid, msg=msg)
        self.cid_ref = model.Coord(self.cid, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.cid = self.Cid()
        self.mid_ref = None
        self.cid_ref = None

    def _verify(self, xref):
        unused_pid = self.Pid()
        unused_mid = self.Mid()
        unused_cid = self.Cid()
        #stress_strain_output_location = self.stress_strain_output_location
        if xref:
            assert self.mid_ref.type in ['MATHE', 'MATHP'], 'PLPLANE: mid_ref.type=%s' % self.mid_ref.type

    #def Pid(self):
        #return self.pid

    #def Thickness(self):
        #return 0.

    #def MassPerArea(self):
        #return 0.

    def Mid(self):
        """returns the material id"""
        if self.mid_ref is not None:
            return self.mid_ref.mid
        return self.mid

    def Cid(self):
        if self.cid_ref is not None:
            return self.cid_ref.cid
        return self.cid

    def raw_fields(self):
        list_fields = ['PLPLANE', self.pid, self.Mid(), self.Cid(),
                       self.stress_strain_output_location]
        return list_fields

    def repr_fields(self):
        list_fields = ['PLPLANE', self.pid, self.Mid(), self.Cid(),
                       self.stress_strain_output_location]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PPLANE(Property):
    type = 'PPLANE'
    _field_map = {1: 'pid', 2:'mid', 3:'t', 4:'nsm', 5:'formulation_option'}

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        return PPLANE(pid, mid, t=0., nsm=0., formulation_option=0, comment='')

    def __init__(self, pid, mid, t=0., nsm=0., formulation_option=0, comment=''):
        # type: (int, int, float, float, int, str) -> None
        """NX specific card"""
        Property.__init__(self)
        if comment:
            self.comment = comment
        #: Property ID
        self.pid = pid
        self.mid = mid
        self.t = t
        self.nsm = nsm
        self.formulation_option = formulation_option
        self.mid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PPLANE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')  # MATHE, MATHP
        t = double_or_blank(card, 3, 'cid', 0.)
        nsm = double_or_blank(card, 4, 'nsm', 0.)
        formulation_option = integer_or_blank(card, 5, 'formulation_option', 0)
        return PPLANE(pid, mid, t, nsm, formulation_option,
                      comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PPLANE pid=%s' % self.pid
        self.mid_ref = model.Material(self.mid, msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

    def _verify(self, xref):
        unused_pid = self.Pid()
        unused_mid = self.Mid()
        #stress_strain_output_location = self.stress_strain_output_location
        if xref:
            assert self.mid_ref.type in ['MAT1', 'MAT3', 'MAT4', 'MAT5', 'MAT8'], 'PPLANE: mid.type=%s' % self.mid.type

    #def Pid(self):
        #return self.pid

    def Thickness(self):
        """returns the thickness of the element"""
        return self.t

    def MassPerArea(self, tflag=1, tscales=None):
        """
        Calculates mass per area.

        .. math:: \frac{m}{A} = nsm + \rho t"""
        mid_ref = self.mid_ref
        rho = mid_ref.Rho()
        thickness = self.Thickness()  # Thickness(tflag=tflag, tscales=tscales)
        try:
            mass_per_area = self.nsm + rho * thickness
        except:
            print("nsm=%s rho=%s t=%s" % (self.nsm, rho, self.t))
            raise
        return mass_per_area

    def Mid(self):
        """returns the material id"""
        if self.mid_ref is not None:
            return self.mid_ref.mid
        return self.mid

    def raw_fields(self):
        list_fields = ['PPLANE', self.pid, self.Mid(), self.t, self.nsm, self.formulation_option]
        return list_fields

    def repr_fields(self):
        list_fields = ['PPLANE', self.pid, self.Mid(), self.t, self.nsm, self.formulation_option]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PSHEAR(Property):
    """
    Defines the properties of a shear panel (CSHEAR entry).

    +--------+-----+-----+---+-----+----+----+
    |   1    |  2  |  3  | 4 |  5  |  6 |  7 |
    +========+=====+=====+===+=====+====+====+
    | PSHEAR | PID | MID | T | NSM | F1 | F2 |
    +--------+-----+-----+---+-----+----+----+
    """
    type = 'PSHEAR'
    _field_map = {1: 'pid', 2:'mid', 3:'t', 4:'nsm', 5:'f1', 6:'f2'}
    pname_fid_map = {
        # 1 based
        4 : 't', 'T' : 't',
    }

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        t = 0.1
        return PSHEAR(pid, mid, t, nsm=0., f1=0., f2=0., comment='')

    def __init__(self, pid, mid, t, nsm=0., f1=0., f2=0., comment=''):
        # type: (int, int, float, float, float, float, str) -> None
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
        Property.__init__(self)
        if comment:
            self.comment = comment
        #: Property ID
        self.pid = pid
        self.t = t
        #: Material ID
        self.mid = mid
        self.nsm = nsm
        self.f1 = f1
        self.f2 = f2
        #assert self.f1 >= 0.0
        #assert self.f2 >= 0.0
        assert self.t > 0.0
        self.mid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PSHEAR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        t = double(card, 3, 't')
        nsm = double_or_blank(card, 4, 'nsm', 0.0)
        f1 = double_or_blank(card, 5, 'f1', 0.0)
        f2 = double_or_blank(card, 6, 'f2', 0.0)
        assert len(card) <= 7, 'len(PSHEAR card) = %i\ncard=%s' % (len(card), card)
        return PSHEAR(pid, mid, t, nsm=nsm, f1=f1, f2=f2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PSHEAR card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        pid = data[0]
        mid = data[1]
        t = data[2]
        nsm = data[3]
        f1 = data[4]
        f2 = data[5]
        assert isinstance(mid, integer_types), data
        return PSHEAR(pid, mid, t, nsm=nsm, f1=f1, f2=f2, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PSHEAR pid=%s' % self.pid
        self.mid_ref = model.Material(self.mid, msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

    def Thickness(self):
        """returns the thickness of the element"""
        return self.t

    def Rho(self):
        """returns the material density"""
        return self.mid_ref.Rho()

    def Mid(self):
        """returns the material id"""
        if self.mid_ref is not None:
            return self.mid_ref.mid
        return self.mid

    def MassPerArea(self, tflag=1, tscales=None):
        """
        Calculates mass per area.

        .. math::  \frac{m}{A} = nsm + \rho t"""
        rho = self.Rho()
        mass_per_area = self.nsm + rho * self.t
        return mass_per_area

    def _verify(self, xref):
        pid = self.Pid()
        midi = self.Mid()

        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(midi, integer_types), 'mid=%r' % midi

        if xref:
            assert isinstance(self.mid_ref, Material), 'mid=%r' % self.mid_ref
            assert self.mid_ref.type in ['MAT1', ], 'pid.type=%s mid.type=%s' % (self.type, self.mid_ref.type)
            if self.mid_ref.type == 'MAT1':
                E = self.mid_ref.E()
                G = self.mid_ref.G()
                nu = self.mid_ref.Nu()
                rho = self.mid_ref.Rho()
                assert isinstance(E, float), 'E=%r' % E
                assert isinstance(G, float), 'G=%r' % G
                assert isinstance(nu, float), 'nu=%r' % nu
                assert isinstance(rho, float), 'rho=%r' % rho

    def raw_fields(self):
        list_fields = ['PSHEAR', self.pid, self.Mid(), self.t, self.nsm,
                       self.f1, self.f2]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PSHELL(Property):
    """
    +--------+-------+------+--------+------+----------+------+------+---------+
    |   1    |   2   |   3  |    4   |  5   |    6     |   7  |  8   |    9    |
    +========+=======+======+========+======+==========+======+======+=========+
    | PSHELL |  PID  | MID1 |   T    | MID2 | 12I/T**3 | MID3 | TS/T |   NSM   |
    +--------+-------+------+--------+------+----------+------+------+---------+
    |        |  Z1   |  Z2  |  MID4  |      |          |      |      |         |
    +--------+-------+------+--------+------+----------+------+------+---------+
    | PSHELL | 41111 |  1   | 1.0000 |  1   |          |   1  |      | 0.02081 |
    +--------+-------+------+--------+------+----------+------+------+---------+
    """
    type = 'PSHELL'
    _field_map = {
        1: 'pid', 2:'mid1', 3:'t', 4:'mid2', 5:'twelveIt3', 6:'mid3',
        7: 'tst', 8:'nsm',
        11:'z1', 12:'z2',
    }
    pname_fid_map = {
        # 1 based
        4 : 't', 'T' : 't',
        6 : 'twelveIt3', # no option
        8 : 'tst', #'T' : 't',
    }

    def __init__(self, pid, mid1=None, t=None, mid2=None, twelveIt3=1.0,
                 mid3=None, tst=0.833333, nsm=0.0,
                 z1=None, z2=None, mid4=None, comment=''):
        # type: (int, Optional[int], float, Optional[int], float, Optional[int], float, float, Optional[float], Optional[float], Optional[int], str) -> None
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
            (only defined if mid2 > 0)
        mid4 : int; default=None
            defines membrane-bending coupling material
            (only defined if mid1 > 0 and mid2 > 0; can't be mid1/mid2)
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
        Property.__init__(self)
        if comment:
            self.comment = comment
        if mid2 == -1:
            mid2 = None

        #: Property ID
        self.pid = pid
        self.mid1 = mid1
        #: Material identification number for bending
        #: -1 for plane strin
        self.mid2 = mid2
        self.mid3 = mid3
        self.mid4 = mid4

        #: thickness
        self.t = t

        #: Scales the moment of interia of the element based on the
        #: moment of interia for a plate
        #:
        #: ..math:: I = \frac{12I}{t^3} I_{plate}
        self.twelveIt3 = twelveIt3

        #: Transverse shear thickness ratio, . Ratio of the shear thickness,
        #: ts/t, to the membrane thickness of the shell, t.
        #: The default value is for a homogeneous shell.
        self.tst = tst

        #: Non-structural Mass
        self.nsm = nsm

        if z1 is None and self.t is not None:
            z1 = -self.t / 2.
        if z2 is None and self.t is not None:
            z2 = self.t / 2.

        self.z1 = z1
        self.z2 = z2

        if self.t is not None:
            assert self.t >= 0.0, 'PSHELL pid=%s Thickness=%s must be >= 0' % (self.pid, self.t)

        self.mid1_ref = None
        self.mid2_ref = None
        self.mid3_ref = None
        self.mid4_ref = None

    @classmethod
    def export_to_hdf5(cls, h5_file, model, pids):
        """exports the properties in a vectorized way"""
        #comments = []
        mids = []
        npids = len(pids)
        assert npids > 0, pids

        z = np.full((npids, 2), np.nan, dtype=None, order='C')
        t = np.full(npids, np.nan, dtype=None, order='C')
        twelveIt3 = []
        tst = []
        nsm = []
        for i, pid in enumerate(pids):
            prop = model.properties[pid]
            #comments.append(prop.comment)
            midsi = [0 if mid is None else mid for mid in
                     [prop.mid1, prop.mid2, prop.mid3, prop.mid4]]
            mids.append(list(midsi))
            z[i, :] = [prop.z1, prop.z2]
            t[i] = prop.t
            twelveIt3.append(prop.twelveIt3)
            tst.append(prop.tst)
            nsm.append(prop.nsm)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('mids', data=mids)
        #print('z =', z)
        #print('t =', t)
        h5_file.create_dataset('z', data=z)
        h5_file.create_dataset('t', data=t)
        h5_file.create_dataset('twelveIt3', data=twelveIt3)
        h5_file.create_dataset('tst', data=tst)
        h5_file.create_dataset('nsm', data=nsm)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PSHELL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        mid1 = integer_or_blank(card, 2, 'mid1')
        t = double_or_blank(card, 3, 't')

        mid2 = integer_or_blank(card, 4, 'mid2')
        twelveIt3 = double_or_blank(card, 5, '12*I/t^3', 1.0)  # poor name
        mid3 = integer_or_blank(card, 6, 'mid3')
        tst = double_or_blank(card, 7, 'ts/t', 0.833333)
        nsm = double_or_blank(card, 8, 'nsm', 0.0)

        if t is not None:
            t_over_2 = t / 2.
            z1 = double_or_blank(card, 9, 'z1', -t_over_2)
            z2 = double_or_blank(card, 10, 'z2', t_over_2)
        else:
            z1 = double_or_blank(card, 9, 'z1')
            z2 = double_or_blank(card, 10, 'z2')
        mid4 = integer_or_blank(card, 11, 'mid4')

        #if self.mid2 is None:
        #    assert self.mid3 is None
        #else: # mid2 is defined
        #    #print (self.mid2 = ", self.mid2)
        #    assert self.mid2 >= -1
        #    #assert self.mid3 >  0

        #if self.mid1 is not None and self.mid2 is not None:
        #    assert self.mid4 == None
        assert len(card) <= 12, 'len(PSHELL card) = %i\ncard=%s' % (len(card), card)
        return PSHELL(pid, mid1, t, mid2, twelveIt3,
                      mid3, tst, nsm,
                      z1, z2, mid4, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PSHELL card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        pid = data[0]
        mid1 = data[1]
        t = data[2]
        mid2 = data[3]
        twelveIt3 = data[4]
        mid3 = data[5]
        tst = data[6]
        nsm = data[7]
        z1 = data[8]
        z2 = data[9]
        mid4 = data[10]
        #maxMid = max(self.mid1,self.mid2,self.mid3,self.mid4)
        #assert self.t > 0.0, ('the thickness must be defined on the PSHELL'
                              #' card (Ti field not supported)')
        return PSHELL(pid, mid1, t, mid2, twelveIt3,
                      mid3, tst, nsm,
                      z1, z2, mid4, comment=comment)

    def _verify(self, xref):
        pid = self.Pid()
        mid = self.Mid()
        mid1 = self.Mid1()
        mid2 = self.Mid2()
        mid3 = self.Mid3()
        mid4 = self.Mid4()

        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(mid, integer_types), 'mid=%r' % mid
        assert mid1 is None or isinstance(mid1, integer_types), 'mid1=%r' % mid1
        assert mid2 is None or isinstance(mid2, integer_types), 'mid2=%r' % mid2
        assert mid3 is None or isinstance(mid3, integer_types), 'mid3=%r' % mid3
        assert mid4 is None or isinstance(mid4, integer_types), 'mid4=%r' % mid4

        mids = [mid for mid in [self.mid1, self.mid2, self.mid3, self.mid4]
                if mid is not None]
        unused_material_ids = self.material_ids
        assert len(mids) > 0
        if xref:
            assert isinstance(self.mid_ref, Material), 'mid=%r' % self.mid_ref

            mids_ref = [self.mid1_ref, self.mid2_ref, self.mid3_ref, self.mid4_ref]
            for i, mid_ref in enumerate(mids_ref):
                if mid_ref is None or mid_ref == 0:
                    continue
                if i == 1: # mid2
                    if isinstance(mid_ref, integer_types):
                        assert mid_ref == -1, mid_ref
                        continue
                assert isinstance(mid_ref, Material), 'mid_ref=%r' % mid_ref
                if mid_ref.type == 'MAT1':
                    E = mid_ref.E()
                    G = mid_ref.G()
                    nu = mid_ref.Nu()
                    rho = mid_ref.Rho()
                    assert isinstance(E, float), 'E=%r' % E
                    assert isinstance(G, float), 'G=%r' % G
                    assert isinstance(nu, float), 'nu=%r' % nu
                    assert isinstance(rho, float), 'rho=%r' % rho
                elif mid_ref.type in ['MAT2', 'MAT4', 'MAT5', 'MAT8']:
                    pass
                #elif mid_ref.type == 'MAT2':
                    #pass
                #elif mid_ref.type == 'MAT4':
                    #pass
                #elif mid_ref.type == 'MAT5':
                    #pass
                #elif mid_ref.type == 'MAT8':
                    #pass
                else:
                    raise NotImplementedError('PSHELL: pid=%s mid_ref.type=%s' % (
                        self.pid, mid_ref.type))

            t = self.Thickness()
            nsm = self.Nsm()
            mpa = self.MassPerArea()
            assert isinstance(t, float), 't=%r' % t
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            assert isinstance(mpa, float), 'mass_per_area=%r' % mpa

    def get_z_locations(self):
        """returns the locations of the bottom and top surface of the shell"""
        z = np.array([self.z1, self.z2])
        return z

    def materials(self):
        """returns the material objects referenced by the shell"""
        materials = [self.mid1_ref, self.mid2_ref, self.mid3_ref, self.mid4_ref]
        return materials

    @property
    def material_ids(self):
        """returns the material ids"""
        return [self.Mid1(), self.Mid2(), self.Mid3(), self.Mid4()]

    #@property
    #def mid(self):
        #raise RuntimeError('use self.mid1, self.mid2, self.mid3, or self.mid4')

    #@mid.setter
    #def mid(self, value):
        #raise RuntimeError('use self.mid1, self.mid2, self.mid3, or self.mid4')

    @property
    def mid_ref(self):
        """returns the material used for mass"""
        if self.mid1_ref is not None:
            return self.mid1_ref
        return self.mid2_ref

    def Mid(self):
        """returns the material id used for mass"""
        mid1 = self.Mid1()
        if mid1 is not None:
            return mid1
        return self.Mid2()

    def Mid1(self):
        """returns the extension material id"""
        if self.mid1_ref is not None:
            return self.mid1_ref.mid
        return self.mid1

    def Mid2(self):
        """returns the bending material id"""
        if self.mid2_ref is not None:
            return self.mid2_ref.mid
        return self.mid2

    def Mid3(self):
        if self.mid3_ref is not None:
            return self.mid3_ref.mid
        return self.mid3

    def Mid4(self):
        if self.mid4_ref is not None:
            return self.mid4_ref.mid
        return self.mid4

    def Thickness(self, tflag=1, tscales=None):
        """returns the thickness of the element"""
        t0 = self.t
        if tscales is not None:
            nt = len(tscales)
            if tflag == 0: # absolute
                thickness = sum([ti if ti is not None else t0 for ti in tscales]) / nt
            elif tflag == 1: # relative
                thickness = sum([ti * t0 if ti is not None else t0 for ti in tscales]) / nt
            else:
                raise RuntimeError('tflag=%r and must be 0/1' % tflag)
            #print('t0 = %s' % t0)
            #print('  tscales = %s' % tscales)
            #print('  nt = %s' % nt)
            #print('  thickness = %s' % thickness)
            return thickness
        else:
            thickness = t0
        return thickness

    def Rho(self):
        """returns the material density"""
        return self.mid_ref.rho

    def Nsm(self):
        """returns the non-structural mass"""
        return self.nsm

    def MassPerArea(self, tflag=1, tscales=None):
        """
        Calculates mass per area.

        .. math:: \frac{m}{A} = nsm + \rho t"""
        mid_ref = self.mid_ref
        rho = mid_ref.Rho()  # fails if mid1=None and mid2=None

        thickness = self.Thickness(tflag=tflag, tscales=tscales)
        try:
            mass_per_area = self.nsm + rho * thickness
        except:
            print("nsm=%s rho=%s t=%s" % (self.nsm, rho, self.t))
            raise
        return mass_per_area

    def MassPerArea_no_xref(self, model, tflag=1, tscales=None):
        """
        Calculates mass per area.

        .. math:: \frac{m}{A} = nsm + \rho t"""
        mid_ref = model.Material(self.Mid())
        rho = mid_ref.Rho()
        thickness = self.Thickness(tflag=tflag, tscales=tscales)
        try:
            mass_per_area = self.nsm + rho * thickness
        except:
            print("nsm=%s rho=%s t=%s" % (self.nsm, rho, self.t))
            raise
        return mass_per_area

    def MassPerArea_structure(self):
        """
        Calculates mass per area without considering non-structural mass.

        .. math:: \frac{m}{A} = nsm + \rho t"""
        mid_ref = self.mid_ref
        rho = mid_ref.Rho()
        try:
            mass_per_area = rho * self.t
        except:
            print("nsm=%s rho=%s t=%s" % (self.nsm, rho, self.t))
            raise
        return mass_per_area

    def get_Qbar_matrix(self, mid_ref, theta=0.):
        """theta must be in radians"""
        S2, unused_S3 = get_mat_props_S(mid_ref)
        T = get_2d_plate_transform(theta)
        #Tinv = np.linalg.inv(T)
        #Qbar = np.linalg.multi_dot([Tinv, Q, Tinv.T])
        Sbar = np.linalg.multi_dot([T.T, S2, T])
        Qbar = np.linalg.inv(Sbar)
        return Qbar

    def get_Sbar_matrix(self, mid_ref, theta=0.):
        """theta must be in radians"""
        # this is the inverse of Sbar
        S2, unused_S3 = get_mat_props_S(mid_ref)
        T = get_2d_plate_transform(theta)
        #Tinv = np.linalg.inv(T)
        Sbar = np.linalg.multi_dot([T.T, S2, T])
        return Sbar

    def get_individual_ABD_matrices(
            self, theta_offset: float=0.) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
        thickness = self.Thickness()

        #z0 = self.z1
        #z1 = self.z2
        #zmean = (z0 + z1) / 2.
        #dz = z1 - z0
        #dzsquared = z1 ** 2 - z0 ** 2
        #zcubed = z1 ** 3 - z0 ** 3
        # A11 A12 A16
        # A12 A22 A26
        # A16 A26 A66
        A = np.zeros((3, 3), dtype='float64')
        B = np.zeros((3, 3), dtype='float64')
        D = np.zeros((3, 3), dtype='float64') # TODO: 2x2 matrix?
        z0 = self.z1
        z1 = self.z2

        if self.mid1_ref:
            Qbar1 = self.get_Qbar_matrix(self.mid1_ref, theta=0.)
            A += Qbar1 * thickness

        if self.mid2_ref:
            Qbar2 = self.get_Qbar_matrix(self.mid2_ref, theta=0.)
            D += Qbar2 * (z1 ** 3 - z0 ** 3) * self.twelveIt3
        #Qbar3 = self.get_Qbar_matrix(self.mid3_ref, theta=0.)
        if self.mid4_ref:
            Qbar4 = self.get_Qbar_matrix(self.mid4_ref, theta=0.)
            B += Qbar4 * (z1 ** 2 - z0 ** 2)

        ts = self.tst * thickness

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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PSHELL pid=%s' % self.pid
        if self.mid1:
            self.mid1_ref = model.Material(self.mid1, msg)
        if self.mid2 and self.mid2 != -1:
            self.mid2_ref = model.Material(self.mid2, msg)
        if self.mid3:
            self.mid3_ref = model.Material(self.mid3, msg)
        if self.mid4:
            self.mid4_ref = model.Material(self.mid4, msg)
        if self.t is not None:
            z1 = abs(self.z1)
            z2 = abs(self.z2)
            t = self.t
            if not ((-1.5*t <= z1 <= 1.5*t) or (-1.5*t <= z2 <= 1.5*t)):
                msg = (f'PSHELL pid={self.pid} midsurface: z1={self.z1:g} z2={self.z2:g} t={t:g} '
                       f'not in range of -1.5t < zi < 1.5t')
                model.log.warning(msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid1 = self.Mid1()
        self.mid2 = self.Mid2()
        self.mid3 = self.Mid3()
        self.mid4 = self.Mid4()
        self.mid1_ref = None
        self.mid2_ref = None
        self.mid3_ref = None
        self.mid4_ref = None

    def raw_fields(self):
        list_fields = ['PSHELL', self.pid, self.Mid1(), self.t, self.Mid2(),
                       self.twelveIt3, self.Mid3(), self.tst, self.nsm, self.z1,
                       self.z2, self.Mid4()]
        return list_fields

    def repr_fields(self):
        twelveIt3 = set_blank_if_default(self.twelveIt3, 1.0)
        tst = set_blank_if_default(self.tst, 0.833333)
        tst2 = set_blank_if_default(self.tst, 0.83333)
        if tst is None or tst2 is None:
            tst = None
        nsm = set_blank_if_default(self.nsm, 0.0)
        if self.t is not None:
            t_over_2 = self.t / 2.
            z1 = set_blank_if_default(self.z1, -t_over_2)
            z2 = set_blank_if_default(self.z2, t_over_2)
        else:
            z1 = self.z1
            z2 = self.z2

        mid1 = self.Mid1()
        mid2 = self.Mid2()
        mid3 = self.Mid3()
        mid4 = self.Mid4()
        mid1 = None if mid1 == 0 else mid1
        mid2 = None if mid2 == 0 else mid2
        mid3 = None if mid3 == 0 else mid3
        mid4 = None if mid4 == 0 else mid4
        list_fields = ['PSHELL', self.pid, mid1, self.t, mid2,
                       twelveIt3, mid3, tst, nsm, z1, z2, mid4]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

def get_2d_plate_transform(theta):
    """theta must be in radians"""
    ct = np.cos(theta)
    st = np.sin(theta)
    ct2 = ct ** 2
    st2 = st ** 2
    cst = st * ct

    T126 = np.array([
        [ct2, st2,    2 * cst],
        [st2, ct2,   -2 * cst],
        [-cst, cst, ct2 - st2],
    ])
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
