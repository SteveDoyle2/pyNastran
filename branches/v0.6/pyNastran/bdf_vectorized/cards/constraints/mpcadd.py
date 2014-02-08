import StringIO
#from collections import defaultdict
#from itertools import izip, count

#from numpy import array

from pyNastran.bdf.fieldWriter import print_card
#from pyNastran.bdf.cards.baseCard import expand_thru
#from pyNastran.bdf2.bdf_interface.assign_type import (integer, integer_or_blank,
    #double, double_or_blank, components, components_or_blank)

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

    def write_bdf(self, f, size=8):
        card = ['MPCADD', self.constraint_id] + self.mpc_ids
        f.write(print_card(card))

    def __repr__(self):
        f = StringIO.StringIO()
        self.write_bdf(f)
        return f.getvalue()