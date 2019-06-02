from itertools import count
from numpy import array, zeros, arange, searchsorted, unique

from pyNastran.dev.bdf_vectorized.cards.elements.property import Property

from pyNastran.bdf.field_writer_8 import print_card_8 #, set_default_if_blank
from pyNastran.bdf.field_writer_16 import print_card_16
#from pyNastran.bdf.field_writer_double import print_card_double

from pyNastran.bdf.field_writer_8 import set_blank_if_default
#from pyNastran.bdf.field_writer_8 import set_string8_blank_if_default
#from pyNastran.bdf.field_writer_16 import set_string16_blank_if_default

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, string, double_or_blank, string_or_blank, double) # fields
#from pyNastran.dev.bdf_vectorized.utils import slice_to_iter


class PBARL(Property):
    type = 'PBARL'
    valid_types = {
        "ROD": 1,
        "TUBE": 2,
        "I": 6,
        "CHAN": 4,
        "T": 4,
        "BOX": 4,
        "BAR": 2,
        "CROSS": 4,
        "H": 4,
        "T1": 4,
        "I1": 4,
        "CHAN1": 4,
        "Z": 4,
        "CHAN2": 4,
        "T2": 4,
        "BOX1": 6,
        "HEXA": 3,
        "HAT": 4,
        "HAT1": 5,
        "DBOX": 10,  # was 12
    }  # for GROUP="MSCBML0"

    def __init__(self, model):
        """
        Defines the PCOMP object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        Property.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            self.model.log.debug('ncards PBARL = %s' % ncards)
            float_fmt = self.model.float_fmt

            #: Property ID
            self.property_id = zeros(ncards, dtype='int32')

            #: Material ID
            self.material_id = zeros(ncards, dtype='int32')

            self.group = zeros(ncards, dtype='|U8')

            #: Section Type (e.g. 'ROD', 'TUBE', 'I', 'H')
            self.Type = zeros(ncards, dtype='|U8')

            #: non-structural mass
            self.nsm = zeros(ncards, dtype=float_fmt)

            #: dimension list
            self.dim = {}

    def add_card(self, card, comment):
        i = self.i
        self.property_id[i] = integer(card, 1, 'property_id')
        self.material_id[i] = integer(card, 2, 'material_id')

        self.group[i] = string_or_blank(card, 3, 'group', 'MSCBMLO')
        Type = string(card, 4, 'Type')
        self.Type[i] = Type

        ndim = self.valid_types[Type]

        dim = []
        for idim in range(ndim):
            dimi = double(card, 9 + idim, 'dim%i' % (idim + 1))
            dim.append(dimi)

        assert len(dim) == ndim, 'PBARL ndim=%s len(dims)=%s' % (ndim, len(dim))
        assert None not in dim

        #: dimension list
        self.dim[i] = dim


        nsm = double_or_blank(card, 9 + ndim + 1, 'nsm', 0.0)
        self.nsm[i] = nsm
        assert isinstance(nsm, float), 'nsm=%r' % nsm

        if Type not in self.valid_types:
            msg = ('Invalid PBARL Type, Type=%s '
                   'valid_types=%s' % (Type, self.valid_types.keys()))
            raise RuntimeError(msg)

        if len(dim) != self.valid_types[Type]:
            msg = 'dim=%s len(dim)=%s Type=%s len(dimType)=%s' % (
                dim, len(dim), Type,
                self.valid_types[Type])
            raise RuntimeError(msg)

        self.i += 1

    def add_op2(self, data):
        i = self.i
        self.property_id[i] = data[0]
        self.material_id[i] = data[1]
        self.group[i] = data[2].strip()
        self.Type[i] = data[3].strip()
        self.dim[i] = list(data[4:-1])
        self.nsm[i] = data[-1]
        #print("group = %r" % self.group)
        #print("Type  = %r" % self.Type)
        #print("dim = ",self.dim)
        #print(str(self))
        #print("*PBARL = ",data)
        #raise NotImplementedError('not finished...')

        if Type not in self.valid_types:
            msg = ('Invalid PBARL Type, Type=%s '
                   'valid_types=%s' % (Type, self.valid_types.keys()))
            raise RuntimeError(msg)

        if len(dim) != self.valid_types[Type]:
            msg = 'dim=%s len(dim)=%s Type=%s len(dimType)=%s' % (
                dim, len(dim), Type,
                self.valid_types[Type])
            raise RuntimeError(msg)
        assert None not in dim
        self.i += 1

    def build(self):
        if self.n:
            unique_pids = unique(self.property_id)
            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PBARL IDs...')

            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_id = self.material_id[i]
            self.group = self.group[i]
            self.Type = self.Type[i]
            self.nsm = self.nsm[i]
            self.dim = {ii : self.dim[j] for ii, j in zip(count(), i)}
        else:
            self.property_id = array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'node' : nid_map,
            'property' : pid_map,
        }
        """
        if self.n:
            nid_map = maps['node']
            pid_map = maps['property']
            for i, pid in enumerate(self.property_id):
                try:
                    self.property_id[i] = pid_map[pid]
                except KeyError:
                    print('pid_map = %s' % pid_map)
                    raise

    #=========================================================================
    def get_index(self, property_ids):
        if isinstance(property_ids, int):
            property_ids = array([property_ids])
        if property_ids is None:
            return arange(self.n)

        indexs = searchsorted(self.property_id, property_ids)
        assert len(indexs) == len(property_ids), 'indexs=%s pids=%s' % (indexs, property_ids)
        return indexs

    def __getitem__(self, property_ids):
        """
        Allows for slicing:
         - elements[1:10]
         - elements[4]
         - elements[1:10:2]
         - elements[[1,2,5]]
         - elements[array([1,2,5])]
        """
        i = searchsorted(self.property_id, property_ids)
        return self.slice_by_index(i)

    #def __getitem__(self, property_ids):
        #property_ids, int_flag = slice_to_iter(property_ids)
        #obj = PBARL(self.model)

        #properties = {}
        #for pid in sorted(property_ids):
            #properties[pid] = self.properties[pid]
        #obj.n = len(property_ids)
        #obj.properties = properties
        #obj.property_id = sorted(self.properties.keys())
        ##obj._comments = obj._comments[index]
        ##obj.comments = obj.comments[index]
        #return obj

    #=========================================================================
    def write_card(self, bdf_file, size=8, property_id=None):
        #self.model.log.debug('PBARL.n = %s' % self.n)
        if self.n:
            if property_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.property_id, property_id)
            #self.model.log.debug('i = %s' % i)
            #cid = [cid if cid != 0 else '' for cid in self.coord_id]

            #group = set_blank_if_default(self.group, 'MSCBMLO')
            #list_fields = ['PBARL', self.pid, self.Mid(), group, self.Type, None,
                           #None, None, None] + self.dim + [self.nsm]

            #self.model.log.debug('*pbarl write pids=%s' % self.property_id)
            for (j, pid, mid, group, Type, nsm) in zip(count(), self.property_id[i],
                                                       self.material_id[i],
                                                       self.group[i], self.Type[i], self.nsm[i]):
                if pid in self._comments:
                    bdf_file.write(self._comments[pid])

                dim = self.dim[j]
                ndim = self.valid_types[Type]
                assert len(dim) == ndim, 'PBARL ndim=%s len(dims)=%s' % (ndim, len(self.dim))

                sgroup = set_blank_if_default(group, 'MSCBMLO')

                list_fields = ['PBARL', pid, mid, group, Type, None,
                               None, None, None] + dim + [nsm]

                #group = set_blank_if_default(self.group, 'MSCBMLO')
                #list_fields = ['PBARL', self.pid, self.Mid(), group, self.Type, None,
                               #None, None, None] + self.dim + [self.nsm]

                if size == 8:
                    bdf_file.write(print_card_8(list_fields))
                else:
                    bdf_file.write(print_card_16(list_fields))

    def slice_by_index(self, i):
        i = self._validate_slice(i)
        obj = PBARL(self.model)
        obj.n = len(i)
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.property_id = self.property_id[i]
        obj.material_id = self.material_id[i]
        obj.Type = self.Type[i]
        obj.group = self.group[i]
        obj.nsm = self.nsm[i]

        dim = {}
        j = 0
        for ii, dimi in self.dim.items():
            if ii in i:
                dim[j] = dimi
                j += 1
        obj.dim = dim

        return obj
