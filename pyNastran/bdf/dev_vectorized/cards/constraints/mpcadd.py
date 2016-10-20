from six.moves import StringIO

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class MPCADD(object):
    """
    ::

      MPCADD SID G1 G2 G3 G4 G5 G6
      G7 G8 G9 -etc.-

      MPCADD SID G1 THRU G2
      MPCADD 313 6 THRU 32
    """
    def __init__(self, model):
        self.model = model
        self.constraint_id = None
        self.mpc_ids = []

    def add(self, constraint_id, mpc_ids, comment):
        #if comment:
            #self._comment = comment
        assert isinstance(constraint_id, int), constraint_id
        self.constraint_id = constraint_id
        self.mpc_ids += mpc_ids

    def build(self):
        self.mpc_ids.sort()

    def write_card(self, bdf_file, size=8):
        card = ['MPCADD', self.constraint_id] + self.mpc_ids
        if size == 8:
            bdf_file.write(print_card_8(card))
        else:
            bdf_file.write(print_card_16(card))

    def __repr__(self):
        f = StringIO()
        self.write_card(f)
        return f.getvalue()
