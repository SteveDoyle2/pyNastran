import unittest

import pyNastran
from pyNastran.bdf.bdf import BDF

root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')

class TestDynamic(unittest.TestCase):
    """
    The cards tested are:
     * TSTEP
    """
    def test_tstep(self):
        """tests a TSTEP card"""
        model = BDF(debug=None)

        sid = 42
        n1 = n2 = 5
        dt1 = dt2 = 0.1
        no1 = no2 = 3
        card = ['TSTEP', sid,
                n1, dt1, no1, None, None, None, None, None,
                n2, dt2, no2]
        model.add_card(card, card[0], comment='tstep comment')
        model.validate()
        tstep = model.tsteps[42]
        tstep.raw_fields()
        tstep.write_card()
        tstep.write_card(size=16)

    def test_tstepnl(self):
        """tests a TSTEPNL card"""
        model = BDF(debug=None)
        card = ['TSTEPNL', 250, 100, .01, 1, 'ADAPT', 2, 10, 'PW',
                1.E-2, 1.E-3, 1.E-6, 2, 10, 2, .02, None,
                5, 5, 0, 0.75, 16.0, 0.1, 20.,]
        model.add_card(card, card[0], comment='tstepnl comment')
        model.validate()
        tstepnl = model.tstepnls[250]
        tstepnl.raw_fields()
        tstepnl.write_card()
        tstepnl.write_card(size=16)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
