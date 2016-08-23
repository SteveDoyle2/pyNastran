import unittest

import pyNastran
from pyNastran.bdf.bdf import PARAM, BDFCard

root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')

class TestOther(unittest.TestCase):
    card = PARAM('NOCOMP', [-1])
    #print('%r' % card)
    assert str(card) == 'PARAM     NOCOMP      -1\n', '%r' % str(card)

    cardi = BDFCard(['PARAM', 'NOCOMP', -1])
    card = PARAM.add_card(cardi)
    assert str(card) == 'PARAM     NOCOMP      -1\n', '%r' % str(card)
    #print('%r' % card)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
