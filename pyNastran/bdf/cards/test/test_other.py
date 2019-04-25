import unittest

#import pyNastran
from pyNastran.bdf.bdf import BDF, PARAM, BDFCard

#ROOT_PATH = pyNastran.__path__[0]
#test_path = os.path.join(ROOT_PATH, 'bdf', 'cards', 'test')

class TestOther(unittest.TestCase):
    def test_param_01(self):
        card = PARAM('NOCOMP', [-1])
        #print('%r' % card)
        assert str(card) == 'PARAM     NOCOMP      -1\n', '%r' % str(card)

        cardi = BDFCard(['PARAM', 'NOCOMP', -1])
        card = PARAM.add_card(cardi)
        assert str(card) == 'PARAM     NOCOMP      -1\n', '%r' % str(card)
        #print('%r' % card)

    def test_param_02(self):
        model = BDF(debug=False)
        param = model.add_param('POST', -1, comment='param post')
        param.raw_fields()
        param.write_card(size=8, is_double=False)
        param.write_card(size=16, is_double=False)
        param.write_card(size=16, is_double=True)
        param.write_card_16(is_double=False)
        param.write_card_16(is_double=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
