#pylint: disable=C0103,C0111
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
#import sys
from itertools import izip
from numpy import array

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Property, Material
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank, double,
    double_or_blank, string_or_blank)
from pyNastran.bdf.fieldWriter16 import print_card_16

class ShellProperty(Property):
    def __init__(self, card, data):
        Property.__init__(self, card, data)

class CompositeShellProperty(ShellProperty):
    def __init__(self, card, data):
        ShellProperty.__init__(self, card, data)

    def cross_reference(self, model):
        """
        Links the Material IDs to the materials.

        :param self:  the PCOMP/PCOMPG object
        :param model: a BDF object
        """
        for iply in xrange(len(self.plies)):
            mid = self.plies[iply][0]
            self.plies[iply][0] = model.Material(mid)  # mid

    def isSymmetrical(self):
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

        :param self: the PCOMP object
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
            if iply < self.nPlies():
                iply = iply - nplies
            else:
                raise IndexError('invalid value for iply=%r' % iply)
        elif iply < 0:
            raise IndexError('invalid value for iply=%r' % iply)
        return iply

    def Thickness(self, iply='all'):
        """
        Gets the thickness of the :math:`i^{th}` ply.

        :param self: the PCOMP object
        :param iply: the string **'all'** (default) or the mass per area of
                     the :math:`i^{th}` ply
        """
        nplies = len(self.plies)
        if iply == 'all':  # get all layers
            t = 0.
            for iply in xrange(nplies):
                t += self.Thickness(iply)

            if self.isSymmetrical():
                return t * 2.
            return t
        else:
            iply = self._adjust_ply_id(iply)
            t = self.plies[iply][1]
            return t

    def nPlies(self):
        r"""
        Gets the number of plies including the core.

        ::

          if Lam=SYM:
            returns nPlies*2   (even)
          else:
            returns nPlies
        """
        nplies = len(self.plies)
        if self.isSymmetrical():
            return nplies * 2
        return nplies

    def Nsm(self):
        """
        Gets the non-structural mass :math:`i^{th}` ply

        :param self: the PCOMP/PCOMPG object
        """
        return self.nsm

    def Mid(self, iply):
        """
        Gets the Material ID of the :math:`i^{th}` ply.

        :param self: the PCOMP/PCOMPG object
        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        Mid = self.Material(iply)
        if isinstance(Mid, int):
            return Mid
        return Mid.mid

    def Mids(self):
        """
        Gets the material IDs of all the plies

        :param self: the PCOMP/PCOMPG object
        :returns mids: the material IDs
        """
        mids = []
        for iply in xrange(self.nPlies()):
            mids.append(self.Mid(iply))
            #theta = self.Theta(iply)
            #sout = self.sout(iply)
        return mids

    def Rho(self, iply):
        """
        Gets the density of the :math:`i^{th}` ply

        :param self: the PCOMP/PCOMPG object
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

        :param self: the PCOMP/PCOMPG object
        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        Mid = self.plies[iply][0]
        return Mid

    def Theta(self, iply):
        """
        Gets the ply angle of the :math:`i^{th}` ply (not the ID)

        :param self: the PCOMP/PCOMPG object
        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        Theta = self.plies[iply][2]
        return Theta

    def sout(self, iply):
        """
        Gets the the flag identifying stress/strain outpur of the
        :math:`i^{th}` ply (not the ID).  default='NO'.

        :param self: the PCOMP/PCOMPG object
        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        Sout = self.plies[iply][3]
        return Sout

    def get_z_locations(self):
        """
        Gets the z locations for the various plies.

        :param self: the PCOMP/PCOMPG object
        :param iply: the ply ID (starts from 0)

        Assume there are 2 plies, each of 1.0 thick, starting from :math:`z=0`.

        >>> pcomp.get_z_locations()
        [0., 1., 2.]
        """
        zi = self.z0
        z = [zi]
        for i in range(self.nPlies()):
            t = self.Thickness(i)
            zi += t
            z.append(zi)
        return array(z)

    def MassPerArea(self, iply='all', method='nplies'):
        r"""
        Gets the Mass/Area for the property.

        .. math:: \frac{m}{A} = \sum(\rho t) + nsm

        or

        .. math:: \frac{m}{A} - nsm = \sum(\rho t)

        and

        .. math:: \frac{m_i}{A} = rho_i t_i + nsm_i

        where :math:`nsm_i` is the non-structural mass of the
        :math:`i^{th}` ply

        :param self:   the PCOMP object
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
            #massPerAreaTotal = m/A = sum(rho*t) + nsm
            #massPerAreaTotal = mpa-nsm = sum(rho*t)
            #(m/A)i = rho*t + nsmi
            # where nsmi has two methods
            massPerArea = 0.
            nplies = len(self.plies)
            for iply in xrange(nplies):
                rho = self.Rho(iply)
                t = self.plies[iply][1]
                massPerArea += rho * t

            if self.isSymmetrical():
                return 2. * massPerArea + self.nsm
            return massPerArea + self.nsm
        else:
            assert isinstance(iply, int), 'iply must be an integer; iply=%r' % iply
            rho = self.Rho(iply)
            t = self.plies[iply][1]

            if method == 'nplies':
                # we divide by nplies b/c it's nsm per area and
                # we're working on a per ply basis
                # nsmi = nsmi/n  # smear based on nplies
                massPerArea = rho * t + self.nsm / self.nPlies()
            elif method == 'rho*t':
                # assume you smear the nsm mass based on rho*t distribution
                #nsmi = rho*t / sum(rho*t) * nsm
                #rho*t + nsmi = rho*t + rho*t/(sum(rho*t) + nsm - nsm) * nsm
                #rho*t + nsmi = rho*t + rho*t/(massPerAreaTotal - nsm) * nsm
                #             = rho*t * (1 + nsm/(massPerAreaTotal-nsm))
                massPerAreaTotal = self.MassPerArea()
                massPerArea = rho * t * (1.0 + self.nsm / (massPerAreaTotal-self.nsm))
            elif method == 't':
                # assume you smear the nsm mass based on t distribution
                #nsmi = t / sum(t) * nsm
                #rho*t + nsmi = rho*t + t/sum(t) * nsm
                #rho*t + nsmi = rho*t + t/thicknessTotal * nsm
                #             = t * (rho + nsm/thicknessTotal)
                thicknessTotal = self.Thickness()
                massPerArea = t * (rho + self.nsm / thicknessTotal)
            else:
                raise NotImplementedError('method=%r is not supported' % method)
            return massPerArea

    def isSameCard(self, prop, debug=False):
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

        for (field1, field2) in izip(fields1, fields2):
            if not self.isSame(field1, field2):
                return False
        return True


class PCOMP(CompositeShellProperty):
    """
    ::

      PCOMP     701512   0.0+0 1.549-2                   0.0+0   0.0+0     SYM
                300704   3.7-2   0.0+0     YES  300704   3.7-2     45.     YES
                300704   3.7-2    -45.     YES  300704   3.7-2     90.     YES
                300705      .5   0.0+0     YES
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
            msg = 'On PCOMP pid=%r, ply %i is not defined.  iply_min=0; iply_max=%i' % (self.pid, ilayer, len(self.plies))
            raise IndexError(msg)

        # ply = [mid, t, theta, sout]
        slot = nnew % 4
        ply[slot] = value

    def __init__(self, card=None, data=None, comment=''):  # not done, cleanup
        ShellProperty.__init__(self, card, data)

        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
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
            self.lam = string_or_blank(card, 8, 'lam')
            assert self.lam in [None, 'SYM', 'MEM', 'BEND', 'SMEAR', 'SMCORE'], 'lam=%r is invalid' % self.lam

            # -8 for the first 8 fields (1st line)
            nPlyFields = card.nFields() - 9

            # counting plies
            nMajor = nPlyFields // 4
            nLeftover = nPlyFields % 4
            if nLeftover:
                nMajor += 1
            nplies = nMajor
            #print("nplies = ",nplies)

            plies = []
            midLast = None
            tLast = None
            ply = None
            iply = 1

            # supports single ply per line
            for i in xrange(9, 9 + nplies * 4, 4):
                actual = card.fields(i, i + 4)
                mid = integer_or_blank(card, i, 'mid', midLast)
                t = double_or_blank(card, i + 1, 't', tLast)
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
                midLast = mid
                tLast = t
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
            for (mid, t, theta, sout) in izip(Mid, T, Theta, Sout):
                if sout == 0:
                    sout = 'NO'
                elif sout == 1:  #: .. todo:: not sure  0=NO,1=YES (most likely)
                    sout = 'YES'
                else:
                    raise RuntimeError('unsupported sout.  sout=%r and must be 0 or 1.'
                                       '\nPCOMP = %s' % (sout, data))
                self.plies.append([mid, t, theta, sout])

    def _verify(self, xref=False):
        pid = self.Pid()
        isSym = self.isSymmetrical()
        nplies = self.nPlies()
        nsm = self.Nsm()
        mids = self.Mids()

        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(isSym, bool), 'isSym=%r' % isSym
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

    def rawFields(self):
        list_fields = ['PCOMP', self.pid, self.z0, self.nsm, self.sb, self.ft,
                  self.TRef, self.ge, self.lam, ]
        for (iply, ply) in enumerate(self.plies):
            (_mid, t, theta, sout) = ply
            mid = self.Mid(iply)
            list_fields += [mid, t, theta, sout]
        return list_fields

    def reprFields(self):
        nsm = set_blank_if_default(self.nsm, 0.0)
        sb = set_blank_if_default(self.sb, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        z0 = set_blank_if_default(self.z0, -0.5 * self.Thickness())

        list_fields = ['PCOMP', self.pid, z0, nsm, sb, self.ft, TRef, ge, self.lam]
        for (iply, ply) in enumerate(self.plies):
            (_mid, t, theta, sout) = ply
            mid = self.Mid(iply)
            #theta = set_blank_if_default(theta,0.0)
            sout = set_blank_if_default(sout, 'NO')
            list_fields += [mid, t, theta, sout]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)


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
            msg = 'On PCOMPG pid=%r, ply %i is not defined.  iply_min=0; iply_max=%i' % (self.pid, ilayer, len(self.plies))
            raise IndexError(msg)

        #ply = [mid, thickness, theta, sout, gPlyID]
        slot = nnew % 5
        ply[slot] = value

    def __init__(self, card=None, data=None, comment=''):
        ShellProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.pid = integer(card, 1, 'pid')
            # z0 will be calculated later
            self.nsm = double_or_blank(card, 3, 'nsm', 0.0)
            self.sb = double_or_blank(card, 4, 'sb', 0.0)
            self.ft = string_or_blank(card, 5, 'ft')
            assert self.ft in ['HILL', 'HOFF', 'TSAI', 'STRN', None]
            self.TRef = double_or_blank(card, 6, 'TRef', 0.0)
            self.ge = double_or_blank(card, 7, 'ge', 0.0)
            self.lam = string_or_blank(card, 8, 'lam')
            fields = card.fields(9)

            T = 0.  # thickness
            midLast = None
            tLast = None
            self.plies = []

            i = 0
            #n = 0
            while i < len(fields):
                gPlyID = integer(card, 9 + i, 'gPlyID')
                mid = integer_or_blank(card, 9 + i + 1, 'mid', midLast)

                # can be blank 2nd time thru
                thickness = double_or_blank(card, 9 + i + 2, 'thickness', tLast)

                theta = double_or_blank(card, 9 + i + 3, 'theta', 0.0)
                sout = string_or_blank(card, 9 + i + 4, 'sout', 'NO')
                #print('n=%s gPlyID=%s mid=%s thickness=%s len=%s' %(
                #    n,gPlyID,mid,thickness,len(fields)))

                ply = [mid, thickness, theta, sout, gPlyID]
                #print("ply = %s" %(ply))
                self.plies.append(ply)
                #[mid,t,theta,sout] # PCOMP

                assert mid is not None
                assert thickness is not None
                assert isinstance(mid, int), 'mid=%s' % mid
                assert isinstance(thickness, float), 'thickness=%s' % thickness
                midLast = mid
                tLast = thickness
                T += thickness
                i += 8
                #n += 1
                self.z0 = double_or_blank(card, 2, 'z0', -0.5 * T)
        else:
            raise NotImplementedError('PCOMPG data')

    def _verify(self, xref=False):
        pid = self.Pid()
        isSym = self.isSymmetrical()
        nplies = self.nPlies()
        nsm = self.Nsm()
        mids = self.Mids()

        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(isSym, bool), 'isSym=%r' % isSym
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
        gPlyID = self.plies[iply][4]
        return gPlyID

    def rawFields(self):
        list_fields = ['PCOMPG', self.pid, self.z0, self.nsm, self.sb, self.ft,
                  self.TRef, self.ge, self.lam, ]
        for (iply, ply) in enumerate(self.plies):
            (_mid, t, theta, sout, gPlyID) = ply
            mid = self.Mid(iply)
            list_fields += [gPlyID, mid, t, theta, sout, None, None, None]
        return list_fields

    def reprFields(self):
        nsm = set_blank_if_default(self.nsm, 0.0)
        sb = set_blank_if_default(self.sb, 0.0)
        TRef = set_blank_if_default(self.TRef, 0.0)
        ge = set_blank_if_default(self.ge, 0.0)
        z0 = set_blank_if_default(self.z0, -0.5 * self.Thickness())

        list_fields = ['PCOMPG', self.pid, z0, nsm, sb, self.ft, TRef, ge,
                       self.lam]
        for (iply, ply) in enumerate(self.plies):
            (_mid, t, theta, sout, gPlyID) = ply
            mid = self.Mid(iply)
            #theta = set_blank_if_default(theta,0.0)
            sout = set_blank_if_default(sout, 'NO')
            list_fields += [gPlyID, mid, t, theta, sout, None, None, None]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)


class PLPLANE(ShellProperty):
    type = 'PLPLANE'
    _field_map = {1: 'pid', 2:'mid', 6:'cid', 7:'str'}

    def __init__(self, card=None, data=None, comment=''):
        ShellProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            self.mid = integer(card, 2, 'mid')  # MATHE, MATHP
            self.cid = integer_or_blank(card, 3, 'cid', 0)
            self.str = string_or_blank(card, 4, 'str', 'GRID')
        else:
            raise NotImplementedError(data)

    def cross_reference(self, model):
        msg = ' which is required by PLPLANE pid=%s' % self.pid
        self.mid = model.Material(self.mid, msg=msg)
        self.cid = model.Coord(self.cid, msg=msg)

    def _verify(self, xref=False):
        pid = self.Pid()
        mid = self.Mid()
        cid = self.Cid()
        str = self.str
        if xref:
            assert self.mid.type in ['MATHE', 'MATHP'], 'mid.type=%s' % self.mid.type

    #def Pid(self):
        #return self.pid

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        return self.mid.mid

    def Cid(self):
        if isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    def rawFields(self):
        list_fields = ['PLPLANE', self.pid, self.Mid(), self.Cid(), self.str]
        return list_fields

    def reprFields(self):
        list_fields = ['PLPLANE', self.pid, self.Mid(), self.Cid(), self.str]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + print_card_8(card)


class PSHEAR(ShellProperty):
    type = 'PSHEAR'
    _field_map = {1: 'pid', 2:'mid', 3:'t', 4:'nsm', 5:'f1', 6:'f2'}

    def __init__(self, card=None, data=None, comment=''):
        """
        Defines the properties of a shear panel (CSHEAR entry).
        PSHEAR PID MID T NSM F1 F2
        """
        ShellProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: Material ID
            self.mid = integer(card, 2, 'mid')
            self.t = double(card, 3, 't')
            self.nsm = double_or_blank(card, 4, 'nsm', 0.0)
            self.f1 = double_or_blank(card, 5, 'f1', 0.0)
            self.f2 = double_or_blank(card, 6, 'f2', 0.0)
            assert self.t > 0.0
            #assert self.f1 >= 0.0
            #assert self.f2 >= 0.0
            assert len(card) <= 7, 'len(PSHEAR card) = %i' % len(card)
        else:
            #(pid,mid,t,nsm,f1,f2) = out
            self.pid = data[0]
            self.mid = data[1]
            self.t = data[2]
            self.nsm = data[3]
            self.f1 = data[4]
            self.f2 = data[5]

    def isSameCard(self, prop, debug=False):
        if self.type != prop.type:
            return False
        fields1 = self.rawFields()
        fields2 = prop.rawFields()
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def rawFields(self):
        list_fields = ['PSHEAR', self.pid, self.Mid(), self.t, self.nsm,
                  self.f1, self.f2]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)


class PSHELL(ShellProperty):
    """
    ::

      PSHELL PID MID1 T MID2 12I/T**3 MID3 TS/T NSM
      Z1 Z2 MID4
      PSHELL   41111   1      1.0000   1               1               0.02081"""
    type = 'PSHELL'
    _field_map = {
        1: 'pid', 2:'mid1', 3:'t', 4:'mid2', 5:'twelveIt3', 6:'mid3',
        7: 'tst', 8:'nsm', 9:'z1', 10:'z2',
    }

    def __init__(self, card=None, data=None, comment=''):
        ShellProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            self.mid1 = integer_or_blank(card, 2, 'mid1')
            #: thickness
            self.t = double(card, 3, 't')

            #: Material identification number for bending
            self.mid2 = integer_or_blank(card, 4, 'mid2')
            #: Scales the moment of interia of the element based on the
            #: moment of interia for a plate
            #:
            #: ..math:: I = \frac{12I}{t^3} I_{plate}
            self.twelveIt3 = double_or_blank(card, 5, '12*I/t^3', 1.0)  # poor name
            self.mid3 = integer_or_blank(card, 6, 'mid3')
            self.tst = double_or_blank(card, 7, 'ts/t', 0.833333)
            #: Non-structural Mass
            self.nsm = double_or_blank(card, 8, 'nsm', 0.0)

            tOver2 = self.t / 2.
            self.z1 = double_or_blank(card, 9,  'z1', -tOver2)
            self.z2 = double_or_blank(card, 10, 'z2',  tOver2)
            self.mid4 = integer_or_blank(card, 11, 'mid4')

            #if self.mid2 is None:
            #    assert self.mid3 is None
            #else: # mid2 is defined
            #    #print "self.mid2 = ",self.mid2
            #    assert self.mid2 >= -1
            #    #assert self.mid3 >   0

            #if self.mid is not None and self.mid2 is not None:
            #    assert self.mid4==None
            assert len(card) <= 12, 'len(PSHELL card) = %i' % len(card)
        else:
            self.pid = data[0]
            self.mid1 = data[1]
            self.t = data[2]
            self.mid2 = data[3]
            self.twelveIt3 = data[4]
            self.mid3 = data[5]
            self.tst = data[6]
            self.nsm = data[7]
            self.z1 = data[8]
            self.z2 = data[9]
            self.mid4 = data[10]
            #maxMid = max(self.mid,self.mid2,self.mid3,self.mid4)

        assert self.t > 0.0, ('the thickness must be defined on the PSHELL'
                              ' card (Ti field not supported)')

    def _verify(self, xref=False):
        pid = self.Pid()
        mid = self.Mid()
        mid1 = self.Mid1()
        mid2 = self.Mid2()
        mid3 = self.Mid3()
        mid4 = self.Mid4()

        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mid, int), 'mid=%r' % mid
        assert mid1 is None or isinstance(mid1, int), 'mid1=%r' % mid1
        assert mid2 is None or isinstance(mid2, int), 'mid2=%r' % mid2
        assert mid3 is None or isinstance(mid3, int), 'mid3=%r' % mid3
        assert mid4 is None or isinstance(mid4, int), 'mid4=%r' % mid4

        mids = [mid for mid in [self.mid1, self.mid2, self.mid3, self.mid4]
                if mid is not None]
        assert len(mids) > 0
        if xref:
            assert isinstance(self.mid(), Material), 'mid=%r' % self.mid()

            for mid in mids:
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

    def mid(self):
        if isinstance(self.mid1, Material):
            return self.mid1
        return self.mid2

    def Mid(self):
        if isinstance(self.mid1, Material):
            return self.mid1.mid
        return self.Mid2()

    def Mid1(self):
        if isinstance(self.mid1, Material):
            return self.mid1.mid
        return self.mid1

    def Mid2(self):
        if isinstance(self.mid2, Material):
            return self.mid2.mid
        return self.mid2

    def Mid3(self):
        if isinstance(self.mid3, Material):
            return self.mid3.mid
        return self.mid3

    def Mid4(self):
        if isinstance(self.mid4, Material):
            return self.mid4.mid
        return self.mid4

    def Thickness(self):
        return self.t

    def Rho(self):
        return self.mid().rho

    def Nsm(self):
        return self.nsm

    def MassPerArea(self):
        """
        Calculates mass per area.

        .. math:: \frac{m}{A} = nsm + \rho t
        """
        try:
            massPerArea = self.nsm + self.Rho() * self.t
        except:
            print("nsm=%s rho=%s t=%s" % (self.nsm, self.Rho(), self.t))
            raise
        return massPerArea

    def D(self):
        return self.mid().Dplate()

    def cross_reference(self, model):
        if self.mid1:
            self.mid1 = model.Material(self.mid1)
        if self.mid2 and self.mid2 != -1:
            self.mid2 = model.Material(self.mid2)
        if self.mid3:
            self.mid3 = model.Material(self.mid3)
        if self.mid4:
            self.mid4 = model.Material(self.mid4)

    def writeCalculix(self, marker='markerDummyProp',
                      elementSet='ELsetDummyProp'):
        msg = '*SHELL SECTION,MATERIAL=M%s_%s,ELSET=%s,OFFSET=%s\n' % (
            marker, self.mid, elementSet, self.z1)
        msg += '** THICKNESS\n'
        msg += '%s\n\n' % (self.t)
        return msg

    def writeCodeAster(self):
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

    def rawFields(self):
        list_fields = ['PSHELL', self.pid, self.Mid1(), self.t, self.Mid2(),
                  self.twelveIt3, self.Mid3(), self.tst, self.nsm, self.z1,
                  self.z2, self.Mid4()]
        return list_fields

    def reprFields(self):
        twelveIt3 = set_blank_if_default(self.twelveIt3, 1.0)
        tst = set_blank_if_default(self.tst, 0.833333)
        nsm = set_blank_if_default(self.nsm, 0.0)

        tOver2 = self.t / 2.
        z1 = set_blank_if_default(self.z1, -tOver2)
        z2 = set_blank_if_default(self.z2, tOver2)

        list_fields = ['PSHELL', self.pid, self.Mid1(), self.t, self.Mid2(),
                  twelveIt3, self.Mid3(), tst, nsm, z1, z2, self.Mid4()]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        if size == 16:
            return self.comment() + print_card_16(card)
        return self.comment() + card_writer(card)