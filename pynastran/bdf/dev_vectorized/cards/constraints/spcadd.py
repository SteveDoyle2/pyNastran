from six.moves import StringIO
#from itertools import count

#from numpy import array

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.dev_vectorized.bdf_interface.assign_type import integer

def get_spcadd_constraint(card):
    constraint_id = integer(card, 1, 'constraint_id')

    node_ids = card.fields(2)
    node_ids = expand_thru(node_ids)

    assert isinstance(constraint_id, int), constraint_id
    return constraint_id, node_ids


class SPCADD(object):
    """
    ::

      SPCADD SID S1 S2 S3 S4 S5 S6
      S7 S8 S9 -etc.-

      SPCADD SID S1 S2 S2
      SPCADD 313  1  2  3
    """
    type = 'SPCADD'
    def __init__(self, model):
        self.model = model
        self.spc_id = None
        self.spc_ids = []

    def add(self, spc_id, spc_ids, comment):
        #if comment:
            #self._comment = comment
        assert isinstance(spc_id, int), spc_id
        self.spc_id = spc_id
        self.spc_ids += spc_ids

    def build(self):
        self.spc_ids.sort()

    def write_card(self, f, size=8):
        card = ['SPCADD', self.spc_id] + self.spc_ids
        #print "card = ", card
        if size == 8:
            f.write(print_card_8(card))
        else:
            f.write(print_card_16(card))

    def __repr__(self):
        f = StringIO()
        self.write_card(f)
        return f.getvalue().rstrip()
