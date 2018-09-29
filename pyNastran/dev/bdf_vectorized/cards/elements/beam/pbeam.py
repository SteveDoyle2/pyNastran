from numpy import array, zeros, asarray

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.dev.bdf_vectorized.cards.elements.property import Property

from pyNastran.bdf.cards.properties.beam import PBEAM as vPBEAM
from pyNastran.dev.bdf_vectorized.utils import slice_to_iter


class PBEAM(Property):
    type = 'PBEAM'

    def __iter__(self):
        pids = self.property_id
        for pid in pids:
            yield pid

    def values(self):
        pids = self.property_id
        for pid in pids:
            yield self.__getitem__(pid)

    def items(self):
        pids = self.property_id
        for pid in pids:
            yield pid, self.__getitem__(pid)

    def __getitem__(self, property_id):
        property_id, int_flag = slice_to_iter(property_id)
        obj = PBEAM(self.model)

        properties = {}
        for pid in sorted(property_id):
            properties[pid] = self.properties[pid]
        obj.n = len(property_id)
        obj.properties = properties
        obj.property_id = sorted(self.properties.keys())
        #obj._comments = obj._comments[index]
        #obj.comments = obj.comments[index]
        return obj

    def slice_by_index(self, i):
        i = asarray(i)
        obj = PBEAM(self.model)
        raise NotImplementedError()
        return obj

    def __init__(self, model):
        """
        Defines the PBEAM object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.properties = {}
        self.model = model
        self.n = 0

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            pass

    def add_card(self, card, comment=''):
        prop = vPBEAM(card, comment=comment)
        self.properties[prop.pid] = prop

    def build(self):
        self.n = len(self.properties)
        self.property_id = array(sorted(self.properties.keys()), dtype='int32')

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
                self.property_id[i] = pid_map[pid]

    #=========================================================================
    def write_card(self, bdf_file, size=8, property_id=None):
        if size == 8:
            for pid, prop in sorted(self.properties.items()):
                bdf_file.write(prop.write_card(size, print_card_8))
        else:
            for pid, prop in sorted(self.properties.items()):
                bdf_file.write(prop.write_card(size, print_card_16))
