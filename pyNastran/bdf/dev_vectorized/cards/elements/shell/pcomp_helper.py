from __future__ import print_function

from numpy import array

from pyNastran.utils import integer_types
#from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
                                                    #double_or_blank, integer_double_or_blank, blank, string_or_blank)
from pyNastran.bdf.bdfInterface.assign_type import integer, double_or_blank, string_or_blank, integer_or_blank
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.fieldWriter import print_card


class BaseCard(object):
    def __init__(self):
        pass

    def comment(self):
        if hasattr(self, '_comment'):
            return '%s' % self._comment
        return ''

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
    def __init__(self, card, data):
        assert card is None or data is None

    def Pid(self):
        """
        returns the property ID of an property

        :returns pid: the Property ID
        :type pid:    int
        """
        return self.pid

    def Mid(self):
        """
        returns the material ID of an element

        :returns mid: the Material ID
        :type mid:    int
        """
        if isinstance(self.mid, integer_types):
            return self.mid
        else:
            return self.mid_ref.mid

    def cross_reference(self, model):
        msg = ' which is required by %s pid=%s' % (self.type, self.pid)
        self.mid = model.Material(self.mid, msg)
        self.mid_ref = self.mid

    def uncross_reference(self):
        self.mid = self.Mid()
        del self.mid_ref


class ShellProperty(Property_i):
    def __init__(self, card, data):
        Property_i.__init__(self, card, data)

class DeprecatedCompositeShellProperty(object):
    def MassPerArea(self, iply='all', method='nplies'):
        return self.get_mass_per_area(iply, method)

    def Thickness(self, iply='all'):
        return self.get_thickness(iply)

    def nPlies(self):
        return self.get_nplies()

    def Nsm(self):
        return self.get_nonstructural_mass()

    def isSymmetrical(self):
        return self.is_symmetrical()

    def Rho(self, iply):
        return self.get_density(iply)

    def Theta(self, iply):
        return self.get_theta(iply)

    def sout(self, iply):
        return self.get_sout(iply)


class CompositeShellProperty(ShellProperty, DeprecatedCompositeShellProperty):
    def __init__(self, card, data):
        ShellProperty.__init__(self, card, data)
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

        :param iply: the ply ID
        :raises: IndexError if iply is invalid

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
            if iply < self.get_nplies():
                iply = iply - nplies
            else:
                raise IndexError('invalid value for iply=%r' % iply)
        elif iply < 0:
            raise IndexError('invalid value for iply=%r' % iply)
        return iply

    def get_thickness(self, iply='all'):
        """
        Gets the thickness of the :math:`i^{th}` ply.

        :param iply: the string **'all'** (default) or the mass per area of
                     the :math:`i^{th}` ply
        """
        nplies = len(self.plies)
        if iply == 'all':  # get all layers
            t = 0.
            for iply in range(nplies):
                t += self.get_thickness(iply)

            if self.isSymmetrical():
                return t * 2.
            return t
        else:
            iply = self._adjust_ply_id(iply)
            t = self.plies[iply][1]
            return t

    def get_nplies(self):
        r"""
        Gets the number of plies including the core.

        ::

          if Lam=SYM:
            returns nPlies*2   (even)
          else:
            returns nPlies
        """
        nplies = len(self.plies)
        if self.is_symmetrical():
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

        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        Mid = self.Material(iply)
        if isinstance(Mid, int):
            return Mid
        return Mid.mid

    def get_material_ids(self):
        return self.Mids()

    def Mids(self):
        """
        Gets the material IDs of all the plies

        :returns mids: the material IDs
        """
        mids = []
        for iply in range(self.nPlies()):
            mids.append(self.Mid(iply))
            #theta = self.get_theta(iply)
            #sout = self.get_sout(iply)
        return mids

    def get_density(self, iply):
        """
        Gets the density of the :math:`i^{th}` ply

        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        mid = self.Material(iply)
        #print("rho =", mid.rho)
        return mid.rho

    def Material(self, iply):
        """
        Gets the material of the :math:`i^{th}` ply (not the ID unless
        it is not cross-referenced).

        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        Mid = self.plies[iply][0]
        return Mid

    def get_theta(self, iply):
        """
        Gets the ply angle of the :math:`i^{th}` ply (not the ID)

        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        Theta = self.plies[iply][2]
        return Theta

    def get_sout(self, iply):
        """
        Gets the the flag identifying stress/strain outpur of the
        :math:`i^{th}` ply (not the ID).  default='NO'.

        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        sout = self.plies[iply][3]
        return sout

    def get_z_locations(self):
        """
        Gets the z locations for the various plies.

        Assume there are 2 plies, each of 1.0 thick, starting from :math:`z=0`.

        >>> pcomp.get_z_locations()
        [0., 1., 2.]
        """
        zi = self.z0
        z = [zi]
        for i in range(self.get_nplies()):
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

            if self.is_symmetrical():
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
                mass_per_area = rho * t + self.nsm / self.get_nplies()
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
    ::

      PCOMP     701512   0.0+0 1.549-2                   0.0+0   0.0+0     SYM
                300704   3.7-2   0.0+0     YES  300704   3.7-2     45.     YES
                300704   3.7-2    -45.     YES  300704   3.7-2     90.     YES
                300705      .5   0.0+0     YES
    """
    type = 'PCOMP'

    def __init__(self, card=None, data=None, comment=''):  # not done, cleanup
        ShellProperty.__init__(self, card, data)

        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')

            # z0 is field 2 and is calculated at the end because we need the
            # thickness first
            #self.z0 = double_or_blank(card, 1, 'pid')

            #: Non-Structural Mass per unit Area
            self.nsm = double_or_blank(card, 3, 'nsm', 0.0)

            self.sb = double_or_blank(card, 4, 'sb', 0.0)
            #: Failure Theory
            #:
            #:   ['HILL', 'HOFF', 'TSAI', 'STRN', None]
            self.ft = string_or_blank(card, 5, 'ft')
            assert self.ft in ['HILL', 'HOFF', 'TSAI', 'STRN', None]

            #: Reference Temperature (default=0.0)
            self.TRef = double_or_blank(card, 6, 'TRef', 0.0)
            self.ge = double_or_blank(card, 7, 'ge', 0.0)

            #: symmetric flag - default = No Symmetry (NO)
            self.lam = string_or_blank(card, 8, 'lam', 'BLANK')
            assert self.lam in ['BLANK', 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE'], 'lam=%r is invalid' % self.lam

            # -8 for the first 8 fields (1st line)
            nply_fields = card.nfields - 9

            # counting plies
            nmajor = nply_fields // 4
            nleftover = nply_fields % 4
            if nleftover:
                nmajor += 1
            nplies = nmajor
            #print("nplies = ",nplies)

            plies = []
            mid_last = None
            thick_last = None
            ply = None
            iply = 1

            # supports single ply per line
            for i in range(9, 9 + nplies * 4, 4):
                actual = card.fields(i, i + 4)
                mid = integer_or_blank(card, i, 'mid', mid_last)
                t = double_or_blank(card, i + 1, 't', thick_last)
                theta = double_or_blank(card, i + 2, 'theta', 0.0)
                sout = string_or_blank(card, i + 3, 'sout', 'NO')

                if not t > 0.:
                    msg = ('thickness of PCOMP layer is invalid pid=%s'
                           ' iLayer=%s t=%s ply=[mid,t,theta,'
                           'sout]=%s' % (self.pid, iply, t, ply))
                    raise RuntimeError(msg)

                # if this card has 2 plies on the line
                if actual != [None, None, None, None]:
                    ply = [mid, t, theta, sout]
                    #print('ply =', ply)
                    plies.append(ply)
                    iply += 1
                mid_last = mid
                thick_last = t
            #print "nplies = ",nplies

            #: list of plies
            self.plies = plies

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
            self.z0 = double_or_blank(card, 2, 'z0', -0.5 * self.Thickness())
        else:
            self.pid = data[0]
            self.z0 = data[1]
            self.nsm = data[2]
            self.sb = data[3]
            self.ft = data[4]
            self.TRef = data[5]
            self.ge = data[6]
            self.lam = data[7]
            assert self.lam in ['SYM', 'NO'], "lam=%r and must be 'SYM'."
            Mid = data[8]
            T = data[9]
            Theta = data[10]
            Sout = data[11]

            self.plies = []
            #ply = [mid,t,theta,sout]
            for (mid, t, theta, sout) in zip(Mid, T, Theta, Sout):
                if sout == 0:
                    sout = 'NO'
                elif sout == 1:  #: .. todo:: not sure  0=NO,1=YES (most likely)
                    sout = 'YES'
                else:
                    raise RuntimeError('unsupported sout.  sout=%r and must be 0 or 1.'
                                       '\nPCOMP = %s' % (sout, data))
                self.plies.append([mid, t, theta, sout])

    def raw_fields(self):
        list_fields = ['PCOMP', self.pid, self.z0, self.nsm, self.sb, self.ft,
                       self.TRef, self.ge, self.lam, ]
        for (iply, ply) in enumerate(self.plies):
            (_mid, t, theta, sout) = ply
            mid = self.Mid(iply)
            list_fields += [mid, t, theta, sout]
        return list_fields

    def repr_fields(self):
        nsm = set_blank_if_default(self.nsm, 0.0)
        sb = set_blank_if_default(self.sb, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        z0 = set_blank_if_default(self.z0, -0.5 * self.get_thickness())

        list_fields = ['PCOMP', self.pid, z0, nsm, sb, self.ft, TRef, ge, self.lam]
        for (iply, ply) in enumerate(self.plies):
            (_mid, t, theta, sout) = ply
            mid = self.Mid(iply)
            #theta = set_blank_if_default(theta,0.0)
            sout = set_blank_if_default(sout, 'NO')
            list_fields += [mid, t, theta, sout]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


