from io import StringIO
from numpy import array, zeros, unique

from pyNastran.dev.bdf_vectorized.utils import slice_to_iter

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import integer
from pyNastran.bdf.cards.properties.bush import PBUSH as vPBUSH

class PBUSH:
    type = 'PBUSH'

    def __init__(self, model):
        """
        Defines the PCOMP object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model
        self.properties = {}
        self.n = 0
        #self._cards = []
        #self._comments = []

    def add_card(self, card: BDFCard, comment: str=''):
        prop = vPBUSH(card, comment=comment)
        self.properties[prop.pid] = prop
        #self._cards.append(card)
        #self._comments.append(comment)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            pass

    def get_properties(self, property_id=None):
        props = []
        if property_id is None:
            property_id = self.property_id
        for pid in self.properties:
            prop = self.properties[pid]
            props.append(prop)
        return props

    def build(self):
        self.n = len(self.properties)
        self.property_id = array(sorted(self.properties.keys()), dtype='int32')
        return
        self.n = 0
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        #return
        if ncards:
            #: Property ID
            self.property_id = zeros(ncards, 'int32')

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

    def write_card(self, bdf_file, size=8, property_ids=None):
        if size == 8:
            for pid, pcomp in sorted(self.properties.items()):
                bdf_file.write(pcomp.write_card(size, print_card_8))
        else:
            for pid, pcomp in sorted(self.properties.items()):
                bdf_file.write(pcomp.write_card(size, print_card_16))

    #def slice_by_index(self, i):
        #i = self._validate_slice(i)
        #return self.__getitem__()

    def __getitem__(self, property_ids):
        property_ids, int_flag = slice_to_iter(property_ids)
        obj = PBUSH(self.model)

        properties = {}
        for pid in sorted(property_ids):
            properties[pid] = self.properties[pid]
        obj.n = len(property_ids)
        obj.properties = properties
        obj.property_id = sorted(self.properties.keys())
        #obj._comments = obj._comments[index]
        #obj.comments = obj.comments[index]
        return obj


    def __repr__(self):
        f = StringIO()
        f.write('<PBUSH object> n=%s\n' % self.n)
        self.write_card(f)
        return f.getvalue()
