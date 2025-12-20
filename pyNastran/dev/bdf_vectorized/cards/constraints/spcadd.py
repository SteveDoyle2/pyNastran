from io import StringIO

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import integer

def get_spcadd_constraint(card):
    constraint_id = integer(card, 1, 'constraint_id')

    node_ids = card.fields(2)
    node_ids = expand_thru(node_ids)

    assert isinstance(constraint_id, int), constraint_id
    return constraint_id, node_ids


class SPCADD:
    """
    Defines a single-point constraint set as a union of single-point constraint
    sets defined on SPC or SPC1 entries.

    +--------+----+----+-----+
    |    1   | 2  |  3 |  4  |
    +========+====+====+=====+
    | SPCADD  | 2 | 1 |   3  |
    +---------+---+---+------+
    """
    type = 'SPCADD'
    def __init__(self, model):
        self.model = model
        self.spc_id = None
        self.spc_ids = []

    def add(self, spc_id, spc_ids, comment):
        #if comment:
            # self.comment = comment
        assert isinstance(spc_id, int), spc_id
        self.spc_id = spc_id
        self.spc_ids += spc_ids

    def build(self):
        self.spc_ids.sort()

    def write_card(self, bdf_file, size=8):
        card = ['SPCADD', self.spc_id] + self.spc_ids
        #print("card = ", card)
        if size == 8:
            bdf_file.write(print_card_8(card))
        else:
            bdf_file.write(print_card_16(card))

    def __repr__(self):
        f = StringIO()
        self.write_card(f)
        return f.getvalue().rstrip()
