from __future__ import annotations
from typing import TYPE_CHECKING
from numpy import array

from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
                                                    #double_or_blank, integer_double_or_blank, blank, string_or_blank)
from pyNastran.bdf.bdf_interface.assign_type import integer, double_or_blank, string_or_blank, integer_or_blank
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer import print_card
from pyNastran.bdf.cards.base_card import _format_comment
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class BaseCard:
    def __init__(self):
        pass

    @property
    def comment(self):
        if hasattr(self, '_comment'):
            return '%s' % self._comment
        return ''

    @comment.setter
    def comment(self, new_comment):
        """sets a comment"""
        #comment = new_comment.rstrip()
        #self._comment = comment + '\n' if comment else ''
        self._comment = _format_comment(new_comment)

    def print_card(self, size=8):
        list_fields = self.repr_fields()
        return self.comment + print_card(list_fields, size=size)


    def __repr__(self):
        """
        Prints a card in the simplest way possible
        (default values are left blank).
        """
        try:
            return self.print_card()
        except:
            print('problem printing %s card' % self.type)
            fields = self.repr_fields()
            print("fields = ", fields)
            raise


class Property_i(BaseCard):
    def __init__(self):
        pass

    def Pid(self):
        """
        returns the property ID of an property

        Returns
        -------
        pid : int
            the Property ID
        """
        return self.pid

    def Mid(self):
        """
        returns the material ID of an element

        Returns
        -------
        mid : int
            the Material ID
        """
        if isinstance(self.mid, integer_types):
            return self.mid
        else:
            return self.mid_ref.mid

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s pid=%s' % (self.type, self.pid)
        self.mid = model.Material(self.mid, msg)
        self.mid_ref = self.mid

    def uncross_reference(self) -> None:
        self.mid = self.Mid()
        del self.mid_ref


class ShellProperty(Property_i):
    def __init__(self):
        Property_i.__init__(self)

class DeprecatedCompositeShellProperty:
    def MassPerArea(self, iply='all', method='nplies'):
        return self.get_mass_per_area(iply, method)

    def Thickness(self, iply='all'):
        return self.get_thickness(iply)

    def nPlies(self):
        return self.nplies

    def Nsm(self):
        return self.get_nonstructural_mass()

    def Rho(self, iply):
        return self.get_density(iply)

    def Theta(self, iply):
        return self.get_theta(iply)

    def sout(self, iply):
        return self.get_sout(iply)


class CompositeShellProperty(ShellProperty, DeprecatedCompositeShellProperty):
    def __init__(self):
        ShellProperty.__init__(self)
        self.nsm = 0.0

    def is_symmetrical(self):
        """
        Is the laminate symmetrical?

        :returns; True or False
        """
        if self.lam == 'SYM':
            return True
        return False

    def _adjust_ply_id(self, iply):
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

        nplies = len(self.plies)
        if iply >= nplies:
            if iply < self.nplies:
                iply = iply - nplies
            else:
                raise IndexError('invalid value for iply=%r' % iply)
        elif iply < 0:
            raise IndexError('invalid value for iply=%r' % iply)
        return iply

    def get_thickness(self, iply='all'):
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
        nplies = len(self.plies)
        if iply == 'all':  # get all layers
            t = 0.
            for iply in range(nplies):
                t += self.get_thickness(iply)

            if self.is_symmetrical:
                return t * 2.
            return t
        else:
            iply = self._adjust_ply_id(iply)
            t = self.plies[iply][1]
            return t

    @property
    def nplies(self):
        r"""
        Gets the number of plies including the core.

        ::

          if Lam=SYM:
            returns nplies*2   (even)
          else:
            returns nplies
        """
        nplies = len(self.plies)
        if self.is_symmetrical:
            return nplies * 2
        return nplies

    def get_nonstructural_mass(self):
        """
        Gets the non-structural mass :math:`i^{th}` ply
        """
        return self.nsm

    def Mid(self, iply):
        """
        Gets the Material ID of the :math:`i^{th}` ply.

        Parameters
        ----------
        iply : int
            the ply ID (starts from 0)

        Returns
        -------
        material_id : int
            the material id of the ith ply
        """
        iply = self._adjust_ply_id(iply)
        Mid = self.Material(iply)
        if isinstance(Mid, integer_types):
            return Mid
        return Mid.mid

    def get_material_ids(self):
        return self.Mids()

    def Mids(self):
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
            #theta = self.get_theta(iply)
            #sout = self.get_sout(iply)
        return mids

    def get_density(self, iply):
        """
        Gets the density of the :math:`i^{th}` ply

        Parameters
        ----------
        iply : int
            the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        mid = self.Material(iply)
        #print("rho =", mid.rho)
        return mid.rho

    def Material(self, iply):
        """
        Gets the material of the :math:`i^{th}` ply (not the ID unless
        it is not cross-referenced).

        Parameters
        ----------
        iply : int
            the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        Mid = self.plies[iply][0]
        return Mid

    def get_theta(self, iply):
        """
        Gets the ply angle of the :math:`i^{th}` ply (not the ID)

        Parameters
        ----------
        iply : int
            the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        Theta = self.plies[iply][2]
        return Theta

    def get_sout(self, iply):
        """
        Gets the the flag identifying stress/strain outpur of the
        :math:`i^{th}` ply (not the ID).  default='NO'.

        Parameters
        ----------
        iply : int
            the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        sout = self.plies[iply][3]
        return sout

    def get_z_locations(self):
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
        for i in range(self.nplies):
            t = self.get_thickness(i)
            zi += t
            z.append(zi)
        return array(z)

    def get_mass_per_area(self, iply='all', method='nplies'):
        r"""
        Gets the Mass/Area for the property.

        .. math:: \frac{m}{A} = \sum(\rho t) + nsm

        or

        .. math:: \frac{m}{A} - nsm = \sum(\rho t)

        and

        .. math:: \frac{m_i}{A} = rho_i t_i + nsm_i

        where :math:`nsm_i` is the non-structural mass of the
        :math:`i^{th}` ply

        :param iply:   the string 'all' (default) or the mass per area of
                       the :math:`i^{th}` ply
        :param method: the method to compute MassPerArea

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
        rhos = [ply[0].get_density() for ply in self.plies]
        return self.get_mass_per_area_rho(rhos, iply, method)

    def get_mass_per_area_rho(self, rhos, iply='all', method='nplies'):
        r"""
        Gets the Mass/Area for the property.

        .. math:: \frac{m}{A} = \sum(\rho t) + nsm

        or

        .. math:: \frac{m}{A} - nsm = \sum(\rho t)

        and

        .. math:: \frac{m_i}{A} = rho_i t_i + nsm_i

        where :math:`nsm_i` is the non-structural mass of the
        :math:`i^{th}` ply

        :param iply:   the string 'all' (default) or the mass per area of
                       the :math:`i^{th}` ply
        :param method: the method to compute MassPerArea

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
        nplies = len(self.plies)
        iply = self._adjust_ply_id(iply)
        if iply == 'all':  # get all layers
            #mass_per_area_total = m/A = sum(rho*t) + nsm
            #mass_per_area_total = mpa-nsm = sum(rho*t)
            #(m/A)i = rho*t + nsmi
            # where nsmi has two methods
            mass_per_area = 0.
            nplies = len(self.plies)
            for iply in range(nplies):
                #rho = self.get_density(iply)
                rho = rhos[iply]
                t = self.plies[iply][1]
                mass_per_area += rho * t

            if self.is_symmetrical:
                return 2. * mass_per_area + self.nsm
            return mass_per_area + self.nsm
        else:
            assert isinstance(iply, int), 'iply must be an integer; iply=%r' % iply
            #rho = self.get_density(iply)
            rho = rhos[iply]
            t = self.plies[iply][1]

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


class PCOMPi(CompositeShellProperty):
    """
    +-------+--------+--------+---------+------+--------+--------+-------+------+
    |   1   |    2   |    3   |    4    |  5   |    6   |    7   |   8   |   9  |
    +=======+========+========+=========+======+========+========+=======+======+
    | PCOMP |   PID  |   Z0   |   NSM   |  SB  |   FT   |  TREF  |  GE   | LAM  |
    +-------+--------+--------+---------+------+--------+--------+-------+------+
    | MID1  |   T1   | THETA1 |  SOUT1  | MID2 |   T2   | THETA2 | SOUT2 |      |
    +-------+--------+--------+---------+------+--------+--------+-------+------+
    | MID3  |   T3   | THETA3 |  SOUT3  | etc. |        |        |       |      |
    +-------+--------+--------+---------+------+--------+--------+-------+------+

    +-------+--------+--------+---------+------+--------+--------+-------+------+
    | PCOMP | 701512 | 0.0+0  | 1.549-2 |      |        | 0.0+0  | 0.0+0 | SYM  |
    +-------+--------+--------+---------+------+--------+--------+-------+------+
    |       | 300704 | 3.7-2  |  0.0+0  | YES  | 300704 | 3.7-2  |   45. | YES  |
    +-------+--------+--------+---------+------+--------+--------+-------+------+
    |       | 300704 | 3.7-2  |  -45.   | YES  | 300704 | 3.7-2  |   90. | YES  |
    +-------+--------+--------+---------+------+--------+--------+-------+------+
    |       | 300705 |   .5   |  0.0+0  | YES  |        |        |       |      |
    +-------+--------+--------+---------+------+--------+--------+-------+------+
    """
    type = 'PCOMP'

    def __init__(self, card=None, data=None, comment=''):  # not done, cleanup
        ShellProperty.__init__(self, card, data)

    @property
    def plies(self):
        plies = []
        for mid, t, theta, sout in zip(self.mids, self.thicknesses, self.thetas, self.souts):
            plies.append([mid, t, theta, sout])
        return plies

    def __init__(self, pid,
                 mids, thicknesses, thetas, souts,
                 nsm, sb, ft, tref, ge, lam, z0, comment=''):
        CompositeShellProperty.__init__(self)
        if comment:
            self.comment = comment

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
        #if lam is None:  # TODO: is NO an option?
            #lam = 'NO'
        self.lam = lam
        self.mids = mids
        self.thicknesses = thicknesses
        self.thetas = thetas
        self.souts = souts
        if z0 is None:
            z0 = -0.5 * self.Thickness()
        self.z0 = z0

        assert self.ft in ['HILL', 'HOFF', 'TSAI', 'STRN', 0.0, None], 'ft=%r' % self.ft
        # TODO: is NO an option?
        assert self.lam in [None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE'], 'lam=%r is invalid' % self.lam

    @classmethod
    def add_card(cls, card, comment=''):
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
        ply = None
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

            if t <= 0.:
                msg = ('thickness of PCOMP layer is invalid pid=%s'
                       ' iLayer=%s t=%s ply=[mid,t,theta,'
                       'sout]=%s' % (pid, iply, t, ply))
                raise RuntimeError(msg)

            # if this card has 2 plies on the line
            if actual != [None, None, None, None]:
                mids.append(mid)
                thicknesses.append(t)
                thetas.append(theta)
                souts.append(sout)
                iply += 1
            mid_last = mid
            thick_last = t
        #print("nplies = ",nplies)

        #self.plies = []
        #if self.lam == 'SYM':
        #    if nplies%2 == 1:  # 0th layer is the core layer
        #       # cut the thickness in half to make the ply have an
        #       # even number of plies, is there a better way???
        #       plies[0][1] = plies[0][1]/2.
        #
        #    pliesLower = plies.reverse()
        #    self.plies = pliesLower+plies
        #    #print str(self)
        z0 = double_or_blank(card, 2, 'z0')
        return PCOMPi(pid, mids, thicknesses, thetas, souts, nsm, sb, ft, tref, ge,
                      lam, z0, comment=comment)

    def update(self, pid_map, mid_map):
        """
        maps = {
            'node' : nid_map,
            'property' : pid_map,
        }
        """
        pid2 = pid_map[self.pid]
        mids2 = [mid_map[mid] if mid != 0 else 0 for mid in self.mids]
        self.pid = pid2
        self.mids = mids2

    @classmethod
    def add_op2_data(cls, data, comment=''):
        #data_in = [
            #pid, z0, nsm, sb, ft, Tref, ge,
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

        #ply = [mid,t,theta,sout]
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
        return PCOMPi(pid, mids, thicknesses, thetas, souts,
                      nsm, sb, ft, tref, ge, lam, z0, comment=comment)

    def raw_fields(self):
        list_fields = ['PCOMP', self.pid, self.z0, self.nsm, self.sb, self.ft,
                       self.tref, self.ge, self.lam, ]
        for (iply, ply) in enumerate(self.plies):
            (_mid, t, theta, sout) = ply
            mid = self.Mid(iply)
            list_fields += [mid, t, theta, sout]
        return list_fields

    def repr_fields(self):
        nsm = set_blank_if_default(self.nsm, 0.0)
        sb = set_blank_if_default(self.sb, 0.0)
        tref = set_blank_if_default(self.tref, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        z0 = set_blank_if_default(self.z0, -0.5 * self.get_thickness())

        list_fields = ['PCOMP', self.pid, z0, nsm, sb, self.ft, tref, ge, self.lam]
        for (iply, ply) in enumerate(self.plies):
            (_mid, t, theta, sout) = ply
            mid = self.Mid(iply)
            #theta = set_blank_if_default(theta,0.0)
            sout = set_blank_if_default(sout, 'NO')
            list_fields += [mid, t, theta, sout]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


