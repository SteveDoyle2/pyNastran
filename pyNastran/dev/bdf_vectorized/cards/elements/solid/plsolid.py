from io import StringIO
from numpy import zeros, unique, searchsorted, array

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    string_or_blank)
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

from pyNastran.dev.bdf_vectorized.cards.elements.property import Property


class PLSOLID(Property):
    type = 'PLSOLID'
    def __init__(self, model):
        """
        Defines the PLSOLID object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        Property.__init__(self, model)

    def allocate(self, ncards):
        self.n = ncards
        #float_fmt = self.model.float_fmt
        #: Property ID
        self.property_id = zeros(ncards, dtype='int32')
        #: Material ID
        self.material_id = zeros(ncards, dtype='int32')
        #: Location of stress and strain output
        self.stress_strain = zeros(ncards, dtype='|U4')

    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i
        self.property_id[i] = integer(card, 1, 'pid')
        self.material_id[i] = integer(card, 2, 'mid')
        stress_strain = string_or_blank(card, 3, 'str', 'GRID')
        if stress_strain not in ('GRID', 'GAUS'):
            msg = 'STR="%s" doesnt have a valid stress/strain ' \
                  'output value set; valid=["GRID", "GAUS"]\n' \
                  % stress_strain
            #raise RuntimeError(msg)
        self.stress_strain[i] = stress_strain
        assert len(card) <= 4, 'len(PLSOLID card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        if self.n:
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_id = self.material_id[i]
            self.stress_strain = self.stress_strain[i]

            unique_pids = unique(self.property_id)
            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PLSOLID IDs...')
        else:
            self.property_id = array([], dtype='int32')
            self.material_id = array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'property' : pid_map,
            'material' : mid_map,
        }
        """
        if self.n:
            nid_map = maps['node']
            pid_map = maps['property']
            mid_map = maps['material']
            for i, pid, mids in enumerate(zip(self.property_id, self.material_ids)):
                self.property_id[i] = pid_map[pid]

    #def get_density_by_property_id(self, property_id=None):
        #if property_id is None:
            #property_id = self.property_id
        #elif isinstance(property_id, int):
            #property_id = array([property_id], dtype='int32')
        #property_id = asarray(property_id)
        #n = len(property_id)
        #rho = zeros(n, dtype='float64')
        #upid = unique(property_id)
        #for i, pid in enumerate(upid):
            ## get the location of pid id in global props
            #j = where(pid == self.property_id)[0]

            ## get the location of pid in the local props
            #k = where(pid == property_id)[0]

            #if len(j) == 0:
                #msg = 'pid=%s was not found in %s' % (upid, self.property_id)
                #raise ValueError(msg)
            #mid = self.material_id[j[0]]
            #rhoi = self.model.materials.get_density_by_material_id([mid])
            ##print('pid=%s rho[%s]=%s' % (pid, k, rhoi))
            #rho[k] = rhoi

        #if rho.min() == 0.0:
            #msg = 'property_id = %s\n' % property_id
            #msg += 'i = %s\n' % i
            #msg += 'rho = %s\n' % rho
            #msg += 'rhoj = %s' % rhoi
            #raise ValueError(msg)

        #assert rho.shape == (n, ), rho.shape
        #self.model.log.debug('rho = %s' % rho)
        #return rho

    def write_card(self, bdf_file, size=8, property_id=None):
        if self.n:
            #print("PSOLID.property_id =", self.property_id)
            for (pid, mid, stress) in zip(
                 self.property_id, self.material_id, self.stress_strain):
                if eid in self._comments:
                    bdf_file.write(self._comments[eid])
                stress_strain = set_blank_if_default(stress_strain, 'GRID')
                card = ['PLSOLID', pid, mid, stress_strain]
                bdf_file.write(print_card_8(card))

    def __getitem__(self, property_id):
        i = searchsorted(self.property_id, property_id)
        return self.slice_by_index(i)

    def slice_by_index(self, i):
        i = self._validate_slice(i)
        obj = PLSOLID(self.model)
        obj.n = len(i)
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.property_id = self.property_id[i]
        obj.material_id = self.material_id[i]
        obj.stress_strain = self.stress_strain[i, :]
        return obj

    def __repr__(self):
        f = StringIO()
        f.write('<PLSOLID object> n=%s\n' % self.n)
        self.write_card(f)
        return f.getvalue()
