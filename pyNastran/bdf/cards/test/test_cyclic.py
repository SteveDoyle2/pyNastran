import copy
import unittest

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.utils import save_load_deck

#from pyNastran.bdf.field_writer_8 import print_card

class TestCyclic(unittest.TestCase):
    """tests cyclic cards"""

    def test_cyclic_01(self):
        """checks the CYAX, CYJOIN cards"""
        model = BDF(debug=False)
        nids = [1, 2]
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        cyax = model.add_cyax(nids, comment='cyax')

        side = 1
        coord = 'R'
        cyjoin = model.add_cyjoin(side, coord, nids, comment='cyjoin')

        cyax.raw_fields()
        cyjoin.raw_fields()
        model.validate()
        save_load_deck(model, run_save_load_hdf5=False)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
