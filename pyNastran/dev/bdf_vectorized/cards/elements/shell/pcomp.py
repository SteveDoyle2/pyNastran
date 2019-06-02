from io import StringIO
from itertools import count

from numpy import (array, zeros, searchsorted, where, unique,
                   hstack)

from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.dev.bdf_vectorized.utils import slice_to_iter
#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

from pyNastran.dev.bdf_vectorized.cards.elements.property import Property
from pyNastran.dev.bdf_vectorized.cards.elements.shell.pcomp_helper import PCOMPi

class PCOMP(Property):
    type = 'PCOMP'

    def allocate(self, ncards):
        self.n = 0

    def __init__(self, model):
        """
        Defines the PCOMP object.

        Parameters
        ----------
        model : BDF
           the BDF object
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

    def add_card(self, card, comment=''):
        prop = PCOMPi.add_card(card, comment=comment)
        self.properties[prop.pid] = prop
        #self._cards.append(card)
        #self._comments.append(comment)

    #def _get_base_parameter_by_property_id(self, property_id, name):
        #int_flag = True if isinstance(property_id, integer_types) else False

        #if property_id is None:
            #property_id = self.property_id
        #elif isinstance(property_id, list):
            #property_id = array(property_id, dtype='int32')
        #elif isinstance(property_id, integer_types):
            #property_id = array([property_id], dtype='int32')
        ##mass_per_area = self.nsm + self.Rho() * self.t

        #n = len(property_id)
        #thickness = zeros(n, dtype='float64')
        #upid = unique(property_id)
        #param = getattr(self, name)
        #for pid in upid:
            #j = searchsorted(self.property_id, pid)
            #i = where(pid == property_id)[0][0]
            #t = 0.0
            #for material_id, ti in zip(self.material_id[j, :], self.t[j, :]):
                #if material_id == 0:
                    #break
                #t += ti
            #self.model.log.debug('t[pid=%s] = %s' % (pid, t))
            #self.model.log.debug('i[pid=%s] = %s' % (pid, i))
            #thickness[i] = t
        #return thickness[0] if int_flag else thickness

    def get_thickness_by_property_id(self, property_id=None):
        int_flag = True if isinstance(property_id, integer_types) else False

        if property_id is None:
            property_id = self.property_id
        elif isinstance(property_id, list):
            property_id = array(property_id, dtype='int32')
        elif isinstance(property_id, integer_types):
            property_id = array([property_id], dtype='int32')
        #mass_per_area = self.nsm + self.Rho() * self.t

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
            #self.model.log.debug('t[pid=%s] = %s' % (pid, t))
            #self.model.log.debug('i[pid=%s] = %s' % (pid, i))
            sym = self.lam[i]
            if sym == 'SYM':
                t *= 2.0
            thickness[i] = t
        return thickness[0] if int_flag else thickness

    #def get_nonstructural_mass_by_property_id(self, property_id=None):
        #"""
        #Gets the nonstructural mass of the PHSELLs.

        #:param property_id: the property IDs to consider (default=None -> all)
        #"""
        #cow
        #int_flag = True if isinstance(property_id, integer_types) else False
        ##print('get_nonstructural_mass; pids = %s' % property_ids)
        #if property_id is None:
            #nsm = self.nsm
        #else:
            #i = self.get_property_index_by_property_id(property_id)
            ##self.model.log.debug('i = %s' % i)
            #nsm = self.nsm[i]
        #return nsm[0] if int_flag else nsm

    #def get_nonstructural_mass_by_property_id(self, property_id=None):
        #props = self.get_properties(property_id)
        #d = array([prop.get_nonstructural_mass() for prop in props])
        #return d

    #def get_material_ids(self, property_id=None):
        #props = self.get_properties(property_id)
        #d = concatenate([prop.get_material_ids() for prop in props])
        #self.model.log.debug('mids PCOMP = %s' % d)
        #return d

    def get_mass_per_area_by_property_id(self, property_id=None, method=None):
        int_flag = False
        if property_id is None:
            property_id = self.property_id
        elif isinstance(property_id, integer_types):
            int_flag = True #if isinstance(property_id, integer_types) else False
            property_id = [property_id]
            i = self.get_property_index_by_property_id([property_id])
        else:
            i = self.get_property_index_by_property_id(property_id)
        #mass_per_area = self.nsm + self.Rho() * self.t

        n = len(property_id)
        mass_per_area = zeros(n, dtype='float64')
        upid = unique(property_id)
        for pid in upid:
            j = searchsorted(self.property_id, pid)
            i = where(pid == property_id)[0]
            mpa = zeros(n, dtype='float64')
            for material_id, thickness in zip(self.material_id[j, :], self.t[j, :]):
                if material_id == 0:
                    break
                rho = self.model.materials.get_density_by_material_id(material_id)
                mpa += rho * thickness
            if self.lam[i] == 'SYM':
                mpa *= 2.0

            mass_per_area[i] = mpa
        mass_per_area += self.nsm[i]
        return mass_per_area[0] if int_flag else mass_per_area

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

    def set_nonstructural_mass_by_property_id(self, property_id, value):
        i = self.get_property_index_by_property_id(property_id)
        #print(i)
        #print('nsm = %s' % self.nsm)
        return self.set_nonstructural_mass_by_property_index(i, value)

    def get_nonstructural_mass_by_property_id(self, property_id=None):
        i = self.get_property_index_by_property_id(property_id)
        return self.get_nonstructural_mass_by_property_index(i)

    def set_nonstructural_mass_by_property_index(self, i, value):
        self.nsm[i] = value

    def get_nonstructural_mass_by_property_index(self, i=None):
        return self.nsm[i]

    def get_z_locations_by_property_id(self, property_id):
        #i = self.get_property_index_by_property_id(property_id)
        #return self.get_z_locations_by_property_index(i)

    #def get_z_locations_by_property_index(self, i=None):
        #for pid in upid:
            #j = searchsorted(self.property_id, pid)
            #i = where(pid == property_id)[0][0]
        assert isinstance(property_id, integer_types), type(property_id)
        i = 0
        j = searchsorted(self.property_id, property_id)
        #print('self.z0 = %s' % self.z0)
        t = 0.0
        z0 = self.z0[j]
        z = [z0]
        for material_id, ti in zip(self.material_id[j, :], self.t[j, :]):
            if material_id == 0:
                break
            t += ti
            z.append(z0 + t)
        z = array(z, dtype='float64')
        self.model.log.debug('t[pid=%s] = %s' % (property_id, t))
        self.model.log.debug('i[pid=%s] = %s' % (property_id, i))
        #print('t*** = ', t)
        #print('self.z0 = %s' % self.z0)
        return z

    def get_z0_by_property_id(self, property_id=None):
        i = self.get_property_index_by_property_id(property_id)
        return self.get_z0_by_property_index(i)

    def get_z0_by_property_index(self, i=None):
        return self.z0[i]

    def get_material_id_by_property_id_ply(self, property_id, jply):
        int_flag = True if isinstance(property_id, integer_types) else False
        i = self.get_property_index_by_property_id(property_id)
        jply2, nplies = self._adjust_ply_id(i, jply)
        mid = self.material_id[i, jply2]

        if mid.min() <= 0:
            msg = 'invalid jply=%s; mid=%s mats=%s' % (jply, mid, self.material_id[i, :])
            raise IndexError(msg)
        return mid[0] if int_flag else mid

    def _adjust_ply_id(self, i, jply):
        nplies = self.nplies[i]
        if jply < 0:
            msg = 'invalid jply=%s' % (jply)
            raise IndexError(msg)
        elif jply >= nplies:
            msg = 'invalid jply=%s' % (jply)
            raise IndexError(msg)
        else:
            if self.lam[i] == 'SYM':
                if jply >= nplies // 2:
                    jply2 = jply - nplies // 2
                    #print('valid ply %s -> %s' % (jply, jply2))
                else:
                    #print('valid ply', jply)
                    jply2 = jply
            else:
                jply2 = jply
        return jply2, nplies

    def get_thickness_by_property_id_ply(self, property_id, jply):
        int_flag = True if isinstance(property_id, integer_types) else False
        assert int_flag == True, property_id

        i = self.get_property_index_by_property_id(property_id)

        jply2, nplies = self._adjust_ply_id(i, jply)
        thickness = self.t[i, jply2]

        if thickness.min() <= 0.0:
            msg = 'invalid jply=%s; nplies=%s thickness=%s' % (jply, nplies, thickness)
            raise IndexError(msg)
        return thickness[0] if int_flag else thickness

    def get_density_by_property_id_ply(self, property_id, jply):
        int_flag = True if isinstance(property_id, integer_types) else False
        i = self.get_property_index_by_property_id(property_id)
        if jply < 0:
            msg = 'invalid jply=%s' % (jply)
            raise IndexError(msg)

        jply2, nplies = self._adjust_ply_id(i, jply)
        mid = self.material_id[i, jply2]

        if mid.min() <= 0:
            msg = 'invalid jply=%s; material_id=%s' % (jply, mid)
            raise IndexError(msg)

        density = self.model.materials.get_density_by_material_id(mid)
        #density = self.rho[i, j]

        #if density <= 0.0:
            #msg = 'invalid jply=%s; density=%s' % (j, density)
            #raise IndexError(msg)
        return density[0] if int_flag else density

    def get_theta_by_property_id_ply(self, property_id, jply):
        int_flag = True if isinstance(property_id, integer_types) else False
        i = self.get_property_index_by_property_id(property_id)

        jply2, nplies = self._adjust_ply_id(i, jply)
        theta = self.theta[i, jply2]
        return theta[0] if int_flag else theta
        #if mid.min() <= 0:
            #msg = 'invalid jply=%s; material_id=%s' % (j, mid)
            #raise IndexError(msg)

    def get_sout_by_property_id_ply(self, property_id, jply):
        int_flag = True if isinstance(property_id, integer_types) else False
        i = self.get_property_index_by_property_id(property_id)

        jply2, nplies = self._adjust_ply_id(i, jply)
        sout = self.sout[i, jply2]
        return sout[0] if int_flag else sout
        #if mid.min() <= 0:
            #msg = 'invalid jply=%s; material_id=%s' % (j, mid)
            #raise IndexError(msg)

    def get_mass_per_area_by_property_id_ply(self, property_id, jply, method='nplies'):
        int_flag = True if isinstance(property_id, integer_types) else False
        i = self.get_property_index_by_property_id(property_id)
        jply2, nplies = self._adjust_ply_id(i, jply)

        mid = self.material_id[i, jply2]
        thickness = self.t[i, jply2]
        if mid.min() <= 0:
            msg = 'invalid jply=%s; material_id=%s' % (jply, mid)
            raise IndexError(msg)
        if thickness.min() <= 0.0:
            msg = 'invalid jply=%s; thickness=%s' % (jply, thickness)
            raise IndexError(msg)

        density = self.model.materials.get_density_by_material_id(mid)

        # mass/A = t*rho #+ nsm
        mass_per_area = thickness * density
        nsm = self.nsm[i]
        if method == 'nplies':
            #nsm = self.nsm[i]
            nplies = self.get_nplies_by_property_index(i)
            delta = nsm / nplies
        elif method == 'rho*t':
            nplies = self.get_nplies_by_property_index(i)
            t = 0.
            rhot = []
            for j in range(nplies):
                ti = self.get_thickness_by_property_id_ply(property_id, j)
                rhoi = self.get_density_by_property_id_ply(property_id, j)
                t += ti
                rhot.append(rhoi * ti)
            delta = nsm * rhot[jply] / sum(rhot)
        elif method == 't':
            nplies = self.get_nplies_by_property_index(i)
            t = []
            for j in range(nplies):
                ti = self.get_thickness_by_property_id_ply(property_id, j)
                t.append(ti)
            delta = nsm * t[jply] / sum(t)
        else:
            raise NotImplementedError(method)
        mass_per_area += delta

        return mass_per_area[0] if int_flag else mass_per_area

    def get_material_ids_by_property_id(self, property_id=None):
        #int_flag = True if isinstance(property_id, integer_types) else False
        i = self.get_property_index_by_property_id(property_id)
        #if j < 0:
            #msg = 'invalid jply=%s' % (j)
            #raise IndexError(msg)

        jplies = self.nplies[i]
        jply2, nplies = self._adjust_ply_id(i, jplies - 1)

        mid = self.material_id[i, :jply2 + 1]
        is_sym = True
        self.model.log.debug('PCOMP.mid = %s' % mid)
        if mid.min() <= 0:
            msg = 'invalid jply=%s; mid=%s mats=%s' % (j, mid, self.material_id[i, :])
            raise IndexError(msg)
        return mid if not is_sym else hstack(list(mid) * 2)
        #return mid[0] if int_flag else mid

    def make_nplies(self, property_id):
        n = len(property_id)
        self.nplies = zeros(n, dtype='int32')
        for i, pid in enumerate(property_id):
            prop = self.properties[pid]
            npliesi = prop.nplies
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
            float_fmt = self.model.float_fmt

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
            self.ft = zeros((n, nplies), dtype='|U4') # 'HILL', 'HOFF', 'TSAI', 'STRN'

            #: Reference Temperature (default=0.0)
            self.tref = zeros(n, dtype=float_fmt)
            self.ge = zeros(n, dtype=float_fmt)

            #: symmetric flag - default = No Symmetry (NO)
            self.lam = zeros(n, dtype='|U8')

            self.material_id = zeros((n, nplies), dtype='int32')
            self.t = zeros((n, nplies), dtype=float_fmt)
            self.theta = zeros((n, nplies), dtype=float_fmt)
            self.sout = zeros((n, nplies), dtype='|U4') # YES, NO
            self.z0 = zeros(n, dtype=float_fmt)

            for i, (pid, prop) in enumerate(sorted(self.properties.items())):
                self.nsm[i] = prop.nsm
                self.sb[i] = prop.sb
                self.ft[i] = prop.ft
                self.tref[i] = prop.tref
                self.ge[i] = prop.ge
                self.lam[i] = prop.lam
                self.z0[i] = prop.z0
                for iply, (mid, t, theta, sout) in zip(count(), prop.plies):
                    self.material_id[i, iply] = mid
                    self.t[i, iply] = t
                    self.theta[i, iply] = theta
                    self.sout[i, iply] = sout
            #self.model.log.debug('PCOMP.material_id = %s' % self.material_id)

            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            unique_pids = unique(self.property_id)

            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PCOMP IDs...')
        else:
            self.property_id = array([], dtype='int32')
            #self.material_id = array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'node' : nid_map,
            'property' : pid_map,
        }
        """
        if self.n:
            pid_map = maps['property']
            mid_map = maps['material']
            properties2 = {}
            for i, (pid, mids) in enumerate(zip(self.property_id, self.material_id)):
                pid2 = pid_map[pid]
                mids2 = [mid_map[mid] if mid != 0 else 0 for mid in mids]
                prop = self.properties[pid]
                self.property_id[i] = pid2
                self.material_id[i, :] = mids2
                properties2[pid2] = prop
                prop.update(pid_map, mid_map)
            self.properties = properties2

    def write_card(self, bdf_file, size=8, is_double=False, property_id=None):
        if size == 8:
            for pid, pcomp in sorted(self.properties.items()):
                bdf_file.write(pcomp.write_card(size, print_card_8))
        else:
            for pid, pcomp in sorted(self.properties.items()):
                bdf_file.write(pcomp.write_card(size, print_card_16))

    def slice_by_index(self, i):
        i = self._validate_slice(i)
        obj = PCOMP(self.model)
        try:
            obj.n = len(i)
        except TypeError:
            msg = 'i=%s type(i)=%s shape=%s' % (i, type(i), str(i.shape))
            raise TypeError(msg)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.property_id = self.property_id[i]
        obj.nplies = self.nplies[i]
        obj.z0 = self.z0[i]
        obj.lam = self.lam[i]
        obj.ge = self.ge[i]
        obj.tref = self.tref[i]
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

    def __repr__(self):
        f = StringIO()
        f.write('<PCOMP object> n=%s\n' % self.n)
        self.write_card(f)
        return f.getvalue()

