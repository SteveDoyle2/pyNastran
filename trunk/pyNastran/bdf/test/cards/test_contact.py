import StringIO
from six.moves import zip
from math import pi, sqrt
from itertools import count
import copy
import unittest

from pyNastran.bdf.bdf import BDF

from pyNastran.bdf.fieldWriter import print_card

class TestContact(unittest.TestCase):

    def test_contact_01(self):
        model = BDF(debug=False)

        lines = [
            'BSURF          3       1       2       3       4       5       6       7',
            '               8       9      10      11      12      13      14      15',
            '              16      17      18      19      20      21      22      23',
        ]
        card = model.add_card(copy.deepcopy(lines), 'BSURF', is_list=False)
        out = model.bsurf[3].write_bdf(8, None)
        lines2 = out.split('\n')
        for i, (line, line2) in enumerate(zip(lines, lines2)):
            self.assertEqual(line, line2)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
