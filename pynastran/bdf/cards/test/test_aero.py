import unittest

import pyNastran
from pyNastran.bdf.bdf import BDF

root_path = pyNastran.__path__[0]
#test_path = os.path.join(root_path, 'bdf', 'cards', 'test')

comment_bad = 'this is a bad comment'
comment_good = 'this is a good comment\n'
model = BDF(debug=None)
class TestAero(unittest.TestCase):
    """
    The Aero cards are:
     * AEFACT
     * AELINK
     * AELIST
     * AEPARM
     * AESTAT
     * AESURF / AESURFS
     * AERO / AEROS
     * CSSCHD
     * CAERO1 / CAERO2 / CAERO3 / CAERO4 / CAERO5
     * FLFACT
     * FLUTTER
     * GUST
     * MKAERO1 / MKAERO2
     * PAERO1 / PAERO2 / PAERO3
     * SPLINE1 / SPLINE2 / SPLINE4 / SPLINE5
    """
    def test_aefact_1(self):
        data = ['AEFACT', 97, .3, 0.7, 1.0]
        model.add_card(data, data[0], comment_bad, is_list=True)

        data = ['AEFACT', 97, .3, 0.7, 1.0]
        model.add_card(data, data[0], comment_bad, is_list=True)

        data = ['AEFACT', '98', '.3', '0.7', '1.0']
        model.add_card(data, data[0], comment_good, is_list=True)

        msg = 'this is a bad commentAEFACT        97      .3      .7      1.\n'
        aefact97 = model.aefacts[97]
        aefact98 = model.aefacts[98]
        self.assertTrue(all(aefact97.Di == [.3, .7, 1.0]))
        self.assertTrue(all(aefact98.Di == [.3, .7, 1.0]))

        out = aefact97.write_card(8, None)
        self.assertEqual(msg, out)

        msg = 'this is a good comment\nAEFACT        98      .3      .7      1.\n'
        out = aefact98.write_card(8, None)
        self.assertEqual(msg, out)

        data = ['AEFACT', 99, .3, 0.7, 1.0, None, 'cat']
        with self.assertRaises(SyntaxError):
            model.add_card(data, data[0], comment_good, is_list=True)

        data = ['AEFACT', 99, .3, 0.7, 1.0, 'cat']
        with self.assertRaises(SyntaxError):
            model.add_card(data, data[0], comment_good, is_list=True)

        data = ['AEFACT', 99, .3, 0.7, 1.0, 2]
        with self.assertRaises(SyntaxError):
            model.add_card(data, data[0], comment_good, is_list=True)

   # def test_aelink_1(self):
    def test_aelist_1(self):
        data = ['AELIST', 75, 1001, 'THRU', 1075, 1101, 'THRU', 1109, 1201, 1202]
        model.add_card(data, data[0], comment_bad, is_list=True)
        elements = list(range(1001, 1076)) + list(range(1101, 1110)) + [1201, 1202]
        aelist75 = model.aelists[75]
        #print(aelist.elements)
        #print(elements)
        self.assertTrue(elements == aelist75.elements)

        elements = list(range(1001, 1076)) + list(range(1101, 1110)) + [1108, 1202]
        data = ['AELIST', 76, 1001, 'THRU', 1075, 1101, 'THRU', 1109, 1108, 1202]
        model.add_card(data, data[0], comment_bad, is_list=True)
        aelist76 = model.aelists[76]
        #print(aelist76 .elements)
        #print(elements)
        self.assertFalse(elements == aelist76.elements)

        elements = list(set(elements))
        elements.sort()
        self.assertTrue(elements == aelist76.elements)

   # def test_aeparm_1(self):
   # def test_aestat_1(self):
   # def test_aesurf_1(self):
   # def test_aesurfs_1(self):

   # def test_aero_1(self):
   # def test_aeros_1(self):

   # def test_caero1_1(self):
   # def test_caero2_1(self):
   # def test_caero3_1(self):
   # def test_caero4_1(self):
   # def test_caero5_1(self):

   # def test_paero1_1(self):
   # def test_paero2_1(self):
   # def test_paero3_1(self):
   # def test_paero4_1(self):
   # def test_paero5_1(self):

   # def test_spline1_1(self):
   # def test_spline2_1(self):
   # def test_spline3_1(self):
   # def test_spline4_1(self):
   # def test_spline5_1(self):


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
