"""
All shell properties are defined in this file.  This includes:

 * PCOMP
 * PCOMPG
 * PLPLANE
 * PSHEAR
 * PSHELL

All shell properties are ShellProperty and Property objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import zip, range
from numpy import array

from pyNastran.utils import integer_types
from pyNastran.bdf.deprecated import DeprecatedCompositeShellProperty
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Property, Material
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_8 import is_same

class ShellProperty(Property):
    def __init__(self):
        Property.__init__(self)


class CompositeShellProperty(ShellProperty, DeprecatedCompositeShellProperty):
    def __init__(self):
        ShellProperty.__init__(self)
        self.mids = []
        self.thicknesses = []
        self.thetas = []
        self.souts = []
        self.z0 = 0.
        self.nsm = 0.
        self.TRef = 0.
        self.ge = 0.
        self.sb = 0.
        self.ft = None
        self.lam = None

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        for iply in range(len(self.thicknesses)):
            mid = self.mids[iply]
            msg = ' which is required by %s pid=%s iply=%s' % (self.type, self.pid, iply)
            self.mids[iply] = model.Material(mid, msg)
        self.mids_ref = self.mids

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

    def uncross_reference(self):
        for iply in range(len(self.thicknesses)):
            mid = self.Mid(iply)
            self.mids[iply] = mid
        del self.mids_ref

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

        nplies = len(self.thicknesses)
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
        #nplies = len(self.thicknesses)
        if iply == 'all':  # get all layers
            thick = sum(self.thicknesses)
            if self.isSymmetrical():
                return thick * 2.
            return thick
        else:
            iply = self._adjust_ply_id(iply)
            thick = self.thicknesses[iply]
            return thick

    @property
    def nplies(self):
        r"""
            Gets the number of plies including the core.

            ::

              if Lam=SYM:
                returns nplies * 2   (even)
              else:
                returns nplies
            """
        nplies = len(self.thicknesses)
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
        Mid = self.Material(iply)
        if isinstance(Mid, integer_types):
            return Mid
        return Mid.mid

    def Mids(self):
        return self.material_ids

    @property
    def material_ids(self):
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

    def get_density(self, iply):
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

    def Material(self, iply):
        """
        Gets the material of the :math:`i^{th}` ply (not the ID unless
        it is not cross-referenced).

        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        mid = self.mids[iply]
        return mid

    def get_theta(self, iply):
        """
        Gets the ply angle of the :math:`i^{th}` ply (not the ID)

        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        theta = self.thetas[iply]
        return theta

    def get_sout(self, iply):
        """
        Gets the the flag identifying stress/strain outpur of the
        :math:`i^{th}` ply (not the ID).  default='NO'.

        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        sout = self.souts[iply]
        return sout

    def get_z_locations(self):
        """
        Gets the z locations for the various plies.

        :param iply: the ply ID (starts from 0)

        Assume there are 2 plies, each of 1.0 thick, starting from :math:`z=0`.

        >>> pcomp.get_z_locations()
        [0., 1., 2.]
        """
        zi = self.z0
        z = [zi]
        for i in range(self.nplies):
            thick = self.get_thickness(i)
            zi += thick
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
        rhos = [mat.get_density() for mat in self.mids]
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
        nplies = len(self.thicknesses)
        iply = self._adjust_ply_id(iply)
        if iply == 'all':  # get all layers
            #mass_per_area_total = m/A = sum(rho*t) + nsm
            #mass_per_area_total = mpa-nsm = sum(rho*t)
            #(m/A)i = rho*t + nsmi
            # where nsmi has two methods
            mass_per_area = 0.
            #nplies = len(self.thicknesses)
            for iply in range(nplies):
                #rho = self.get_density(iply)
                rho = rhos[iply]
                t = self.thicknesses[iply]
                mass_per_area += rho * t

            if self.is_symmetrical():
                return 2. * mass_per_area + self.nsm
            return mass_per_area + self.nsm
        else:
            assert isinstance(iply, int), 'iply must be an integer; iply=%r' % iply
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

    def _is_same_card(self, prop, debug=False):
        if self.type != prop.type:
            return False
        fields2 = [prop.nsm, prop.sb, prop.ft, prop.TRef, prop.ge, prop.lam]
        fields1 = [self.nsm, self.sb, self.ft, self.TRef, self.ge, self.lam]

        for ply in self.plies:
            fields1 += ply
        for ply in prop.plies:
            fields2 += ply

        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))

        for (field1, field2) in zip(fields1, fields2):
            if not is_same(field1, field2):
                return False
        return True


class PCOMP(CompositeShellProperty):
    """
    +-------+--------+--------+---------+------+--------+--------+-------+------+
    | PCOMP |   PID  |   Z0   |   NSM   | SB   |   FT   | TREF   | GE    | LAM  |
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
    _field_map = {
        1: 'pid', 2: 'z0', 3:'nsm', 4:'sb', 5:'ft', 6:'TRef',
        7: 'ge', 8:'lam',
    }

    def _update_field_helper(self, n, value):
        nnew = n - 9
        if nnew <= 0:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

        ilayer = nnew // 4
        try:
            ply = self.plies[ilayer]
        except IndexError:
            msg = ('On PCOMP pid=%r, ply %i is not defined.  '
                   'iply_min=0; iply_max=%i' % (self.pid, ilayer, len(self.plies)))
            raise IndexError(msg)

        # ply = [mid, t, theta, sout]
        slot = nnew % 4
        ply[slot] = value

    @property
    def plies(self):
        plies = []
        for mid, t, theta, sout in zip(self.mids, self.thicknesses, self.thetas, self.souts):
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

    def __init__(self, pid,
                 mids, thicknesses, thetas, souts,
                 nsm, sb, ft, TRef, ge, lam, z0, comment=''):
        CompositeShellProperty.__init__(self)
        if comment:
            self._comment = comment

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
        self.TRef = TRef
        self.ge = ge

        #: symmetric flag - default = No Symmetry (NO)
        self.lam = lam
        self.mids = mids
        self.thicknesses = thicknesses
        self.thetas = thetas
        self.souts = souts
        if z0 is None:
            z0 = -0.5 * self.Thickness()
        self.z0 = z0

        assert self.ft in ['HILL', 'HOFF', 'TSAI', 'STRN', 0.0, None], 'ft=%r' % self.ft
        assert self.lam in [None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE', 'NO'], 'lam=%r is invalid' % self.lam

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')

        # z0 is field 2 and is calculated at the end because we need the
        # thickness first
        #self.z0 = double_or_blank(card, 1, 'pid')

        nsm = double_or_blank(card, 3, 'nsm', 0.0)
        sb = double_or_blank(card, 4, 'sb', 0.0)
        ft = string_or_blank(card, 5, 'ft')
        TRef = double_or_blank(card, 6, 'TRef', 0.0)
        ge = double_or_blank(card, 7, 'ge', 0.0)
        lam = string_or_blank(card, 8, 'lam')

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
        #print "nplies = ",nplies

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
        return PCOMP(pid, mids, thicknesses, thetas, souts, nsm, sb, ft, TRef, ge,
                     lam, z0, comment=comment)

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
        TRef = data[5]
        ge = data[6]
        lam = data[7]
        Mid = data[8]
        T = data[9]
        Theta = data[10]
        Sout = data[11]

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
        return PCOMP(pid, mids, thicknesses, thetas, souts,
                     nsm, sb, ft, TRef, ge, lam, z0, comment=comment)

    def _verify(self, xref=False):
        pid = self.Pid()
        is_sym = self.isSymmetrical()
        nplies = self.nPlies()
        nsm = self.Nsm()
        mids = self.Mids()

        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(is_sym, bool), 'is_sym=%r' % is_sym
        assert isinstance(nplies, int), 'nplies=%r' % nplies
        assert isinstance(nsm, float), 'nsm=%r' % nsm
        assert isinstance(mids, list), 'mids=%r' % mids

        t = self.Thickness()
        mpa = self.MassPerArea()
        assert isinstance(t, float), 'thickness=%r' % t
        assert isinstance(mpa, float), 'mass_per_area=%r' % mpa
        for iply in range(nplies):
            mid2 = self.Mid(iply)
            assert mids[iply] == mid2
            t = self.Thickness(iply)
            rho = self.Rho(iply)
            mpa = self.MassPerArea(iply)
            assert isinstance(t, float), 'thickness=%r' % t
            assert isinstance(rho, float), 'rho=%r' % rho
            assert isinstance(mpa, float), 'mass_per_area=%r' % mpa
        for ply in self.plies:
            assert len(ply) == 4, ply

    def raw_fields(self):
        list_fields = ['PCOMP', self.pid, self.z0, self.nsm, self.sb, self.ft,
                       self.TRef, self.ge, self.lam, ]
        for (mid, t, theta, sout) in zip(self.material_ids, self.thicknesses,
                                         self.thetas, self.souts):
            list_fields += [mid, t, theta, sout]
        return list_fields

    def repr_fields(self):
        nsm = set_blank_if_default(self.nsm, 0.0)
        sb = set_blank_if_default(self.sb, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        z0 = set_blank_if_default(self.z0, -0.5 * self.get_thickness())

        list_fields = ['PCOMP', self.pid, z0, nsm, sb, self.ft, TRef, ge, self.lam]
        for (mid, t, theta, sout) in zip(self.material_ids, self.thicknesses,
                                         self.thetas, self.souts):
            #theta = set_blank_if_default(theta,0.0)
            str_sout = set_blank_if_default(sout, 'NO')
            list_fields += [mid, t, theta, str_sout]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PCOMPG(CompositeShellProperty):
    type = 'PCOMPG'
    _field_map = {
        1: 'pid', 2: 'z0', 3:'nsm', 4:'sb', 5:'ft', 6:'TRef',
        7: 'ge', 8:'lam',
    }

    def _update_field_helper(self, n, value):
        nnew = n - 9
        if nnew <= 0:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

        ilayer = nnew // 5
        try:
            ply = self.plies[ilayer]
        except IndexError:
            msg = ('On PCOMPG pid=%r, ply %i is not defined.  '
                   'iply_min=0; iply_max=%i' % (self.pid, ilayer, len(self.plies)))
            raise IndexError(msg)

        #ply = [mid, thickness, theta, sout, global_ply_id]
        slot = nnew % 5
        ply[slot] = value

    @property
    def plies(self):
        plies = []
        for mid, t, theta, sout, global_ply_id in zip(self.mids, self.thicknesses,
                                                      self.thetas, self.souts,
                                                      self.global_ply_ids):
            plies.append((mid, t, theta, sout, global_ply_id))
        return plies

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

    def __init__(self, pid,
                 global_ply_ids, mids, thicknesses, thetas, souts,
                 nsm, sb, ft, TRef, ge, lam, z0, comment=''):
        CompositeShellProperty.__init__(self)
        if comment:
            self._comment = comment

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
        self.TRef = TRef
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

        assert self.ft in ['HILL', 'HOFF', 'TSAI', 'STRN', None], 'ft=%r' % self.ft
        assert self.lam in [None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE'], 'lam=%r is invalid' % self.lam

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')
        # z0 will be calculated later
        nsm = double_or_blank(card, 3, 'nsm', 0.0)
        sb = double_or_blank(card, 4, 'sb', 0.0)
        ft = string_or_blank(card, 5, 'ft')
        TRef = double_or_blank(card, 6, 'TRef', 0.0)
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
            assert isinstance(mid, int), 'mid=%s' % mid
            assert isinstance(thickness, float), 'thickness=%s' % thickness
            mid_last = mid
            thick_last = thickness
            #T += thickness
            i += 8
            #n += 1
        z0 = double_or_blank(card, 2, 'z0')
        return PCOMPG(pid,
                      global_ply_ids, mids, thicknesses, thetas, souts,
                      nsm, sb, ft, TRef, ge, lam, z0, comment=comment)

    def _verify(self, xref=False):
        pid = self.Pid()
        is_sym = self.is_symmetrical()
        nplies = self.nplies
        nsm = self.Nsm()
        mids = self.Mids()

        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(is_sym, bool), 'is_sym=%r' % is_sym
        assert isinstance(nplies, int), 'nplies=%r' % nplies
        assert isinstance(nsm, float), 'nsm=%r' % nsm
        assert isinstance(mids, list), 'mids=%r' % mids

        t = self.Thickness()
        mpa = self.MassPerArea()
        assert isinstance(t, float), 'thickness=%r' % t
        assert isinstance(mpa, float), 'mass_per_area=%r' % mpa
        for iply in range(nplies):
            glply = self.GlobalPlyID(iply)
            mid2 = self.Mid(iply)
            assert mids[iply] == mid2
            t = self.Thickness(iply)
            rho = self.Rho(iply)
            mpa = self.MassPerArea(iply)
            assert isinstance(glply, int), 'global_ply_id=%r' % glply
            assert isinstance(t, float), 'thickness=%r' % t
            assert isinstance(rho, float), 'rho=%r' % rho
            assert isinstance(mpa, float), 'mass_per_area=%r' % mpa

    def GlobalPlyID(self, iply):
        global_ply_id = self.global_ply_ids[iply]
        return global_ply_id

    def raw_fields(self):
        list_fields = [
            'PCOMPG', self.pid, self.z0, self.nsm, self.sb, self.ft,
            self.TRef, self.ge, self.lam, ]
        zipi = zip(self.material_ids, self.thicknesses, self.thetas,
                   self.souts, self.global_ply_ids)
        for (mid, t, theta, sout, global_ply_id) in zipi:
            list_fields += [global_ply_id, mid, t, theta, sout, None, None, None]
        return list_fields

    def repr_fields(self):
        nsm = set_blank_if_default(self.nsm, 0.0)
        sb = set_blank_if_default(self.sb, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        z0 = set_blank_if_default(self.z0, -0.5 * self.Thickness())

        list_fields = [
            'PCOMPG', self.pid, z0, nsm, sb, self.ft, TRef, ge,
            self.lam]
        zipi = zip(self.material_ids, self.thicknesses, self.thetas, self.souts,
                   self.global_ply_ids)
        for (mid, t, theta, sout, global_ply_id) in zipi:
            sout = set_blank_if_default(sout, 'NO')
            list_fields += [global_ply_id, mid, t, theta, sout, None, None, None]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PLPLANE(ShellProperty):
    type = 'PLPLANE'
    _field_map = {1: 'pid', 2:'mid', 6:'cid', 7:'str'}

    def __init__(self, pid, mid, cid=0, stress_strain_output_location='GRID', comment=''):
        ShellProperty.__init__(self)
        if comment:
            self._comment = comment
        #: Property ID
        self.pid = pid
        self.mid = mid
        self.cid = cid
        self.stress_strain_output_location = stress_strain_output_location

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')  # MATHE, MATHP
        cid = integer_or_blank(card, 3, 'cid', 0)
        stress_strain_output_location = string_or_blank(card, 4, 'str', 'GRID')
        return PLPLANE(pid, mid, cid, stress_strain_output_location,
                       comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PLPLANE pid=%s' % self.pid
        self.mid = model.HyperelasticMaterial(self.mid, msg=msg)
        self.cid = model.Coord(self.cid, msg=msg)
        self.mid_ref = self.mid
        self.cid_ref = self.cid

    def _verify(self, xref=False):
        pid = self.Pid()
        mid = self.Mid()
        cid = self.Cid()
        #stress_strain_output_location = self.stress_strain_output_location
        if xref:
            assert self.mid.type in ['MATHE', 'MATHP'], 'mid.type=%s' % self.mid.type

    #def Pid(self):
        #return self.pid

    def Mid(self):
        if isinstance(self.mid, integer_types):
            return self.mid
        return self.mid_ref.mid

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def raw_fields(self):
        list_fields = ['PLPLANE', self.pid, self.Mid(), self.Cid(), self.stress_strain_output_location]
        return list_fields

    def repr_fields(self):
        list_fields = ['PLPLANE', self.pid, self.Mid(), self.Cid(), self.stress_strain_output_location]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PSHEAR(ShellProperty):
    type = 'PSHEAR'
    _field_map = {1: 'pid', 2:'mid', 3:'t', 4:'nsm', 5:'f1', 6:'f2'}

    def __init__(self, pid, t, mid, nsm, f1, f2, comment=''):
        """
        Defines the properties of a shear panel (CSHEAR entry).


        +--------+-----+-----+---+-----+----+----+
        |   1    |  2  |  3  | 4 |  5  |  6 |  7 |
        +========+=====+=====+===+=====+====+====+
        | PSHEAR | PID | MID | T | NSM | F1 | F2 |
        +--------+-----+-----+---+-----+----+----+
        """
        ShellProperty.__init__(self)
        if comment:
            self._comment = comment
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

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        t = double(card, 3, 't')
        nsm = double_or_blank(card, 4, 'nsm', 0.0)
        f1 = double_or_blank(card, 5, 'f1', 0.0)
        f2 = double_or_blank(card, 6, 'f2', 0.0)
        assert len(card) <= 7, 'len(PSHEAR card) = %i\ncard=%s' % (len(card), card)
        return PSHEAR(pid, t, mid, nsm, f1, f2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        pid = data[0]
        mid = data[1]
        t = data[2]
        nsm = data[3]
        f1 = data[4]
        f2 = data[5]
        assert isinstance(mid, integer_types), data
        return PSHEAR(pid, t, mid, nsm, f1, f2, comment=comment)

    def _is_same_card(self, prop, debug=False):
        if self.type != prop.type:
            return False
        fields1 = self.raw_fields()
        fields2 = prop.raw_fields()
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self._is_same_fields(fields1, fields2)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PSHEAR pid=%s' % self.pid
        self.mid = model.Material(self.mid, msg)
        self.mid_ref = self.mid

    def Rho(self):
        return self.mid_ref.Rho()

    def Mid(self):
        if isinstance(self.mid, integer_types):
            return self.mid
        return self.mid_ref.mid

    def MassPerArea(self):
        """
        Calculates mass per area.

        .. math::  \frac{m}{A} = nsm + \rho t"""
        rho = self.Rho()
        mass_per_area = self.nsm + rho * self.t
        return mass_per_area

    def _verify(self, xref=False):
        print('xref =', xref)
        pid = self.Pid()
        midi = self.Mid()

        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(midi, integer_types), 'mid=%r' % midi

        if xref:
            assert isinstance(self.mid_ref, Material), 'mid=%r' % self.mid_ref

            assert isinstance(self.mid_ref, Material), 'mid=%r' % self.mid_ref
            assert self.mid_ref.type in ['MAT1', ], 'pid.type=%s mid.type=%s' % (self.type, self.mid_ref.type)
            if self.mid.type == 'MAT1':
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

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PSHELL(ShellProperty):
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
        7: 'tst', 8:'nsm', 9:'z1', 10:'z2',
    }

    def __init__(self, pid, mid1=None, t=None, mid2=None, twelveIt3=1.0,
                 mid3=None, tst=0.833333, nsm=0.0,
                 z1=None, z2=None, mid4=None, comment=''):
        ShellProperty.__init__(self)
        if comment:
            self._comment = comment

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
            z2 = -self.t / 2.

        self.z1 = z1
        self.z2 = z2

        if self.t is not None:
            assert self.t >= 0.0, 'PSHELL pid=%s Thickness=%s must be >= 0' % (self.pid, self.t)

    @classmethod
    def add_card(cls, card, comment=''):
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

    def _verify(self, xref=False):
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
        material_ids = self.material_ids
        assert len(mids) > 0
        if xref:
            assert isinstance(self.mid, Material), 'mid=%r' % self.mid

            for i, mid in enumerate(mids):
                if mid == 0:
                    continue
                if i == 1: # mid2
                    if isinstance(mid, integer_types):
                        assert mid == -1, mid
                        continue
                assert isinstance(mid, Material), 'mid=%r' % mid
                assert mid.type in ['MAT1', 'MAT2', 'MAT4', 'MAT5', 'MAT8'], 'pid.type=%s mid.type=%s' % (self.type, mid.type)
                if mid.type == 'MAT1':
                    E = mid.E()
                    G = mid.G()
                    nu = mid.Nu()
                    rho = mid.Rho()
                    assert isinstance(E, float), 'E=%r' % E
                    assert isinstance(G, float), 'G=%r' % G
                    assert isinstance(nu, float), 'nu=%r' % nu
                    assert isinstance(rho, float), 'rho=%r' % rho
                #elif mid.type == 'MAT2':
                    #pass
                #elif mid.type == 'MAT4':
                    #pass
                #elif mid.type == 'MAT5':
                    #pass
                #elif mid.type == 'MAT8':
                    #pass
                #else:
                    #raise NotImplementedError('pid.type=%s mid.type=%s' % (self.type, mid.type))

            t = self.Thickness()
            nsm = self.Nsm()
            mpa = self.MassPerArea()
            assert isinstance(t, float), 't=%r' % t
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            assert isinstance(mpa, float), 'mass_per_area=%r' % mpa

    def get_z_locations(self):
        z = array([self.z1, self.z2])
        return z

    def materials(self):
        materials = []
        for i in range(1, 5):
            name = 'mid%i_ref' % i
            if hasattr(self, name):
                mat = getattr(self, name)
                materials.append(mat)
            else:
                materials.append(None)
        return materials

    @property
    def material_ids(self):
        return [self.Mid1(), self.Mid2(), self.Mid3(), self.Mid4()]

    #@property
    #def mid(self):
        #raise RuntimeError('use self.mid1, self.mid2, self.mid3, or self.mid4')

    #@mid.setter
    #def mid(self, value):
        #raise RuntimeError('use self.mid1, self.mid2, self.mid3, or self.mid4')

    @property
    def mid(self):
        if isinstance(self.mid1, Material):
            return self.mid1_ref
        return self.mid2_ref

    def Mid(self):
        if isinstance(self.mid1, Material):
            return self.mid1_ref.mid
        return self.Mid2()

    def Mid1(self):
        if isinstance(self.mid1, Material):
            return self.mid1_ref.mid
        return self.mid1

    def Mid2(self):
        if isinstance(self.mid2, Material):
            return self.mid2_ref.mid
        return self.mid2

    def Mid3(self):
        if isinstance(self.mid3, Material):
            return self.mid3_ref.mid
        return self.mid3

    def Mid4(self):
        if isinstance(self.mid4, Material):
            return self.mid4_ref.mid
        return self.mid4

    def Thickness(self):
        return self.t

    def Rho(self):
        return self.mid.rho

    def Nsm(self):
        return self.nsm

    def MassPerArea(self):
        """
        Calculates mass per area.

        .. math:: \frac{m}{A} = nsm + \rho t"""
        rho = self.Rho()
        try:
            mass_per_area = self.nsm + rho * self.t
        except:
            print("nsm=%s rho=%s t=%s" % (self.nsm, rho, self.t))
            raise
        return mass_per_area

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by PSHELL pid=%s' % self.pid
        if self.mid1:
            self.mid1 = model.Material(self.mid1, msg)
            self.mid1_ref = self.mid1
        if self.mid2 and self.mid2 != -1:
            self.mid2 = model.Material(self.mid2, msg)
            self.mid2_ref = self.mid2
        if self.mid3:
            self.mid3 = model.Material(self.mid3, msg)
            self.mid3_ref = self.mid3
        if self.mid4:
            self.mid4 = model.Material(self.mid4, msg)
            self.mid4_ref = self.mid4
        if self.t is not None:
            z1 = abs(self.z1)
            z2 = abs(self.z2)
            t = self.t
            if not ((-1.5*t <= z1 <= 1.5*t) or (-1.5*t <= z2 <= 1.5*t)):
                msg = 'PSHELL pid=%s midsurface: z1=%s z2=%s t=%s not in range of -1.5t < zi < 1.5t' % (
                    self.pid, self.z1, self.z2, t)
                model.log.warning(msg)

    def uncross_reference(self):
        self.mid1 = self.Mid1()
        self.mid2 = self.Mid2()
        self.mid3 = self.Mid3()
        self.mid4 = self.Mid4()

    def _write_calculix(self, marker='markerDummyProp',
                        element_set='ELsetDummyProp'):
        msg = '*SHELL SECTION,MATERIAL=M%s_%s,ELSET=%s,OFFSET=%s\n' % (
            marker, self.mid, element_set, self.z1)
        msg += '** THICKNESS\n'
        msg += '%s\n\n' % (self.t)
        return msg

    def write_code_aster(self):
        """
        * http://www.caelinux.org/wiki/index.php/Contrib:KeesWouters/shell/static
        * http://www.caelinux.org/wiki/index.php/Contrib:KeesWouters/platedynamics

        The angle_rep is a direction angle, use either angle(a,b) or
        vecteur(x,y,z)
        coque_ncou is the number of gauss nodes along the thickness, for
        linear analysis one node is sufficient.
        """
        msg = ''
        msg += "    COQUE=_F(GROUP_MA='P%s', # COQUE=PSHELL\n" % self.pid
        msg += "              EPAIS=%g, # EPAIS=thickness\n" % self.t
        msg += "              ANGL_REP=(0.,90.),  # ???\n"  #: .. todo:: what is this?
        #msg += "              VECTEUR=(1.0,0.0,0.0,)  #  local coordinate system\n"
        msg += "              EXCENTREMENT=%g,  # offset-Z1\n" % self.z1
        msg += "              COQUE_NCOU=1,  # Number of Integration Layers\n"
        msg += "              CARA=('NSM'), # ???\n"  #: .. todo:: check
        msg += "              VALE=(%g),),\n" % self.nsm
        return msg

    def raw_fields(self):
        list_fields = ['PSHELL', self.pid, self.Mid1(), self.t, self.Mid2(),
                       self.twelveIt3, self.Mid3(), self.tst, self.nsm, self.z1,
                       self.z2, self.Mid4()]
        return list_fields

    def repr_fields(self):
        twelveIt3 = set_blank_if_default(self.twelveIt3, 1.0)
        tst = set_blank_if_default(self.tst, 0.833333)
        tst2 = set_blank_if_default(self.tst, 0.83333)
        if tst or tst2 is None:
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

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
