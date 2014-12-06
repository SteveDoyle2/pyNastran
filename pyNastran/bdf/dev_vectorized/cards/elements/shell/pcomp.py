from six import iteritems
from six.moves import StringIO, zip
from itertools import count

from numpy import (array, zeros, arange, searchsorted, where, unique,
                   hstack, concatenate, full, nan, asarray)

from pyNastran.bdf.dev_vectorized.utils import slice_to_iter
from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.fieldWriter import print_card_8, print_card
from pyNastran.bdf.fieldWriter16 import print_card_16
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank, string_or_blank)

from pyNastran.bdf.dev_vectorized.cards.elements.property import Property

class PCOMP(Property):
    type = 'PCOMP'
    def __len__(self):
        return self.n

    def allocate(self, ncards):
        self.n = 0

    def __init__(self, model):
        """
        Defines the PCOMP object.

        :param self: the PCOMP object
        :param model: the BDF object
        :param cards: the list of PCOMP cards
        """
        Property.__init__(self, model)
        #self.model = model
        del self.i
        self.properties = {}
        self.nplies = None
        #self.n = 0
        #self.i = 0
        #self._cards = []
        #self._comments = []

    def add(self, card, comment=''):
        prop = PCOMPi(card, comment=comment)
        self.properties[prop.pid] = prop
        #self._cards.append(card)
        #self._comments.append(comment)

    def get_thickness_by_property_id(self, property_id=None):
        if property_id is None:
            property_id = self.property_id
        elif isinstance(property_id, list):
            property_id = array(property_id, dtype='int32')
        elif isinstance(property_id, int):
            property_id = array([property_id], dtype='int32')
        #massPerArea = self.nsm + self.Rho() * self.t

        n = len(property_id)
        self.model.log.debug('nthickness = %s' % n)
        self.model.log.debug('property_id = %s %s' % (property_id, type(property_id)))
        thickness = zeros(n, dtype='float64')
        self.model.log.debug('thickness = %s %s' % (thickness, type(thickness)))
        upid = unique(property_id)
        for pid in upid:
            j = searchsorted(self.property_id, pid)
            i = where(pid == property_id)[0][0]
            t = 0.0
            for material_id, ti in zip(self.material_id[j, :], self.t[j, :]):
                if material_id == 0:
                    break
                t += ti
            self.model.log.debug('t = %s' % t)
            self.model.log.debug('i = %s' % i)
            thickness[i] = t
        return thickness

    #def get_nonstructural_mass_by_property_id(self, property_id=None):
        #props = self.get_properties(property_id)
        #d = array([prop.get_nonstructural_mass() for prop in props])
        #return d

    #def get_material_ids(self, property_id=None):
        #props = self.get_properties(property_id)
        #d = concatenate([prop.get_material_ids() for prop in props])
        #self.model.log.debug('mids PCOMP = %s' % d)
        #return d

    def get_mass_per_area_by_property_id(self, property_id=None):
        if property_id is None:
            property_id = self.property_id
        #massPerArea = self.nsm + self.Rho() * self.t

        n = len(property_id)
        mass_per_area = zeros(n, dtype='float64')
        upid = unique(property_id)
        for pid in upid:
            j = searchsorted(self.property_id, pid)
            i = where(pid == property_id)[0]
            mpa = self.nsm[i]
            for material_id, thickness in zip(self.material_id[j, :], self.t[j, :]):
                if material_id == 0:
                    break
                rho = self.model.materials.get_density_by_material_id(material_id)
                mpa += rho * thickness
            mass_per_area[i] = mpa
        return mass_per_area

    #def get_property_id_by_property_index(self):
        #asf

    #def get_property_id(self):
        #return self.properties.keys()

    #def get_properties(self, property_id=None):
        #props = []
        #if property_id is None:
            #property_id = self.property_id
        #for pid in self.properties:
            #prop = self.properties[pid]
            #props.append(prop)
        #return props

    def get_nplies_by_property_id(self, property_id=None):
        i = self.get_property_index_by_property_id(property_id)
        return self.get_nplies_by_property_index(i)

    def get_nplies_by_property_index(self, i=None):
        if self.nplies is None:
            raise RuntimeError('PCOMP.build() must be called')
        return self.nplies[i]

    def get_non_structural_mass_by_property_id(self, property_id=None):
        i = self.get_property_index_by_property_id(property_id)
        return self.get_non_structural_mass_by_property_index(i)

    def get_non_structural_mass_by_property_index(self, i=None):
        return self.nsm[i]

    def get_material_id_by_property_id_ply(self, property_id, j):
        i = self.get_property_index_by_property_id(property_id)
        if j < 0:
            msg = 'invalid jply=%s' % (j)
            raise IndexError(msg)
        mid = self.material_id[i, j]

        if mid.min() <= 0:
            msg = 'invalid jply=%s; material_id=%s' % (j, mid)
            raise IndexError(msg)
        return mid

    def get_thickness_by_property_id_ply(self, property_id, j):
        i = self.get_property_index_by_property_id(property_id)
        if j < 0:
            msg = 'invalid jply=%s' % (j)
            raise IndexError(msg)
        thickness = self.t[i, j]

        if thickness.min() <= 0.0:
            msg = 'invalid jply=%s; thickness=%s' % (j, thickness)
            raise IndexError(msg)
        return thickness

    def get_density_by_property_id_ply(self, property_id, j):
        i = self.get_property_index_by_property_id(property_id)
        if j < 0:
            msg = 'invalid jply=%s' % (j)
            raise IndexError(msg)

        mid = self.material_id[i, j]

        if mid.min() <= 0:
            msg = 'invalid jply=%s; material_id=%s' % (j, mid)
            raise IndexError(msg)

        density = self.model.materials.get_density_by_material_id(mid)
        #density = self.rho[i, j]

        #if density <= 0.0:
            #msg = 'invalid jply=%s; density=%s' % (j, density)
            #raise IndexError(msg)
        return density

    def get_theta_by_property_id_ply(self, property_id, j):
        i = self.get_property_index_by_property_id(property_id)
        if j < 0:
            msg = 'invalid jply=%s' % (j)
            raise IndexError(msg)

        theta = self.theta[i, j]
        return theta
        #if mid.min() <= 0:
            #msg = 'invalid jply=%s; material_id=%s' % (j, mid)
            #raise IndexError(msg)

    def get_sout_by_property_id_ply(self, property_id, j):
        i = self.get_property_index_by_property_id(property_id)
        if j < 0:
            msg = 'invalid jply=%s' % (j)
            raise IndexError(msg)

        sout = self.sout[i, j]
        return sout
        #if mid.min() <= 0:
            #msg = 'invalid jply=%s; material_id=%s' % (j, mid)
            #raise IndexError(msg)

    def get_material_ids_by_property_id(self, property_id=None):
        i = self.get_property_index_by_property_id(property_id)
        #if j < 0:
            #msg = 'invalid jply=%s' % (j)
            #raise IndexError(msg)
        jplies = self.nplies[i]
        mid = self.material_id[i, :jplies]
        self.model.log.debug('PCOMP.mid = %s' % mid)
        if mid.min() <= 0:
            msg = 'invalid jply=%s; mid=%s' % (j, mid)
            raise IndexError(msg)
        return mid

    def make_nplies(self, property_id):
        n = len(property_id)
        self.nplies = zeros(n, dtype='int32')
        for i, pid in enumerate(property_id):
            prop = self.properties[pid]
            npliesi = prop.get_nplies()
            self.nplies[i] = npliesi

    def is_symmetrical_by_property_index(self, i=None):
        if i is None:
            j = where(self.lam == 'SYM')[0]
            is_symmetrical = zeros(self.n, dtype='bool')
            is_symmetrical[j] = True
            return is_symmetrical

        n = len(i)
        is_symmetrical = zeros(n, dtype='bool')
        i = where(self.lam[i] == 'SYM')[0]
        is_symmetrical[i] = True
        return is_symmetrical

    def is_symmetrical_by_property_id(self, property_id=None):
        i = self.get_property_index_by_property_id(property_id)
        return self.is_symmetrical_by_property_index(i)

    def build(self):
        n = len(self.properties)
        self.n = n
        if self.n:
            float_fmt = self.model.float

            #: Property ID
            self.property_id = array(sorted(self.properties.keys()), dtype='int32')

            # number of plies
            self.make_nplies(self.property_id)
            self.model.log.debug('self.nplies = %s' % self.nplies)
            nplies = self.nplies.max()
            self.model.log.debug('nplies = %s' % nplies)

            #: Non-Structural Mass per unit Area
            self.nsm = zeros(n, dtype=float_fmt)

            self.sb = zeros(n, dtype=float_fmt)

            #: Failure Theory
            #:
            #:   ['HILL', 'HOFF', 'TSAI', 'STRN', '']
            self.ft = zeros((n, nplies), dtype='|S4') # 'HILL', 'HOFF', 'TSAI', 'STRN'

            #: Reference Temperature (default=0.0)
            self.TRef = zeros(n, dtype=float_fmt)
            self.ge = zeros(n, dtype=float_fmt)

            #: symmetric flag - default = No Symmetry (NO)
            self.lam = zeros(n, dtype='S8')

            self.material_id = zeros((n, nplies), dtype='int32')
            self.t = zeros((n, nplies), dtype=float_fmt)
            self.theta = zeros((n, nplies), dtype=float_fmt)
            self.sout = zeros((n, nplies), dtype='|S4') # YES, NO
            self.z0 = zeros(n, dtype=float_fmt)

            for i, (pid, prop) in enumerate(sorted(iteritems(self.properties))):
                self.nsm[i] = prop.nsm
                self.sb[i] = prop.sb
                self.ft[i] = prop.ft
                self.TRef[i] = prop.TRef
                self.ge[i] = prop.ge
                self.lam[i] = prop.lam
                self.z0[i] = prop.z0
                for iply, (mid, t, theta, sout) in zip(count(), prop.plies):
                    self.material_id[i, iply] = mid
                    self.t[i, iply] = t
                    self.theta[i, iply] = theta
                    self.sout[i, iply] = sout
            self.model.log.debug('PCOMP.material_id = %s' % self.material_id)
            return
            ncards = len(cards)
            for i, card in enumerate(cards):
                self.property_id[i] = integer(card, 1, 'pid')
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            unique_pids = unique(self.property_id)

            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PCOMP IDs...')
            self._cards = []
            self._comments = []
        else:
            self.property_id = array([], dtype='int32')
            #self.material_id = array([], dtype='int32')

    def write_bdf(self, f, size=8, property_id=None):
        if size == 8:
            for pid, pcomp in sorted(iteritems(self.properties)):
                f.write(pcomp.write_bdf(size, print_card_8))
        else:
            for pid, pcomp in sorted(iteritems(self.properties)):
                f.write(pcomp.write_bdf(size, print_card_16))

    #def slice_by_index(self, i):
        #i = asarray(i)
        #return self.__getitem__()

    #def __getitem_old__(self, property_id):
        #property_id, int_flag = slice_to_iter(property_id)
        #obj = PCOMP(self.model)

        #properties = {}
        #for pid in sorted(property_id):
            #properties[pid] = self.properties[pid]
        #obj.n = len(property_id)
        #obj.properties = properties
        #obj.property_id = sorted(self.properties.keys())
        ##obj._comments = obj._comments[index]
        ##obj.comments = obj.comments[index]
        #return obj

    def __getitem__(self, property_id):
        property_id, int_flag = slice_to_iter(property_id)
        i = searchsorted(self.property_id, property_id)
        return self.slice_by_index(i)

    def slice_by_index(self, i):
        i = asarray(i)
        obj = PCOMP(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.property_id = self.property_id[i]
        obj.nplies = self.nplies[i]
        obj.z0 = self.z0[i]
        obj.lam = self.lam[i]
        obj.ge = self.ge[i]
        obj.TRef = self.TRef[i]
        obj.sb = self.sb[i]
        obj.ft = self.ft[i]
        obj.nsm = self.nsm[i]

        obj.material_id = self.material_id[i, :]
        obj.t = self.t[i, :]
        obj.sout = self.sout[i, :]
        obj.theta = self.theta[i, :]

        obj.properties = {}
        for pid in obj.property_id:
            obj.properties[pid] = self.properties[pid]
        return obj

    #def __getitem__(self, property_id):
        #print('PCOMP.property_id = %s' % property_id)
        #property_id, int_flag = slice_to_iter(property_id)
        #properties = {}
        #for pid in sorted(property_id):
            #properties[pid] = self.properties[pid]
        ##print('intflag', int_flag)
        #return properties[pid] if int_flag else properties

    def __repr__(self):
        f = StringIO()
        f.write('<PCOMP object> n=%s\n' % self.n)
        self.write_bdf(f)
        return f.getvalue()

class BaseCard(object):
    def __init__(self):
        pass

    def comment(self):
        if hasattr(self, '_comment'):
            return '%s' % self._comment
        return ''

    def print_card(self, size=8):
        list_fields = self.repr_fields()
        return self.comment() + print_card(list_fields, size=size)


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

        :param self:  the Property pointer
        :returns pid: the Property ID
        :type pid:    int
        """
        return self.pid

    def Mid(self):
        """
        returns the material ID of an element

        :param self:  the Property pointer
        :returns mid: the Material ID
        :type mid:    int
        """
        if isinstance(self.mid, int):
            return self.mid
        else:
            return self.mid.mid

    def cross_reference(self, model):
        msg = ' which is required by %s pid=%s' % (self.type, self.pid)
        self.mid = model.Material(self.mid, msg)

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

        :param self: the PCOMP object
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

    def get_material_ids(self):
        return self.Mids()

    def Mids(self):
        """
        Gets the material IDs of all the plies

        :param self: the PCOMP/PCOMPG object
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

    def get_theta(self, iply):
        """
        Gets the ply angle of the :math:`i^{th}` ply (not the ID)

        :param self: the PCOMP/PCOMPG object
        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        Theta = self.plies[iply][2]
        return Theta

    def get_sout(self, iply):
        """
        Gets the the flag identifying stress/strain outpur of the
        :math:`i^{th}` ply (not the ID).  default='NO'.

        :param self: the PCOMP/PCOMPG object
        :param iply: the ply ID (starts from 0)
        """
        iply = self._adjust_ply_id(iply)
        sout = self.plies[iply][3]
        return sout

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
            for iply in range(nplies):
                #rho = self.get_density(iply)
                rho = rhos[iply]
                t = self.plies[iply][1]
                massPerArea += rho * t

            if self.is_symmetrical():
                return 2. * massPerArea + self.nsm
            return massPerArea + self.nsm
        else:
            assert isinstance(iply, int), 'iply must be an integer; iply=%r' % iply
            #rho = self.get_density(iply)
            rho = rhos[iply]
            t = self.plies[iply][1]

            if method == 'nplies':
                # we divide by nplies b/c it's nsm per area and
                # we're working on a per ply basis
                # nsmi = nsmi/n  # smear based on nplies
                massPerArea = rho * t + self.nsm / self.get_nplies()
            elif method == 'rho*t':
                # assume you smear the nsm mass based on rho*t distribution
                #nsmi = rho*t / sum(rho*t) * nsm
                #rho*t + nsmi = rho*t + rho*t/(sum(rho*t) + nsm - nsm) * nsm
                #rho*t + nsmi = rho*t + rho*t/(massPerAreaTotal - nsm) * nsm
                #             = rho*t * (1 + nsm/(massPerAreaTotal-nsm))
                massPerAreaTotal = self.get_mass_per_area_rho(rhos, iply='all', method='nplies')
                massPerArea = rho * t * (1.0 + self.nsm / (massPerAreaTotal - self.nsm))
            elif method == 't':
                # assume you smear the nsm mass based on t distribution
                #nsmi = t / sum(t) * nsm
                #rho*t + nsmi = rho*t + t/sum(t) * nsm
                #rho*t + nsmi = rho*t + t/thicknessTotal * nsm
                #             = t * (rho + nsm/thicknessTotal)
                thicknessTotal = self.get_thickness()
                massPerArea = t * (rho + self.nsm / thicknessTotal)
            else:
                raise NotImplementedError('method=%r is not supported' % method)
            return massPerArea


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
            for i in range(9, 9 + nplies * 4, 4):
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
            for (mid, t, theta, sout) in zip(Mid, T, Theta, Sout):
                if sout == 0:
                    sout = 'NO'
                elif sout == 1:  #: .. todo:: not sure  0=NO,1=YES (most likely)
                    sout = 'YES'
                else:
                    raise RuntimeError('unsupported sout.  sout=%r and must be 0 or 1.'
                                       '\nPCOMP = %s' % (sout, data))
                self.plies.append([mid, t, theta, sout])

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
        z0 = set_blank_if_default(self.z0, -0.5 * self.get_thickness())

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
        return self.comment() + print_card_8(card)
