from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import unittest
import pyNastran
from pyNastran.bdf.cards.baseCard import collapse_thru, collapse_thru_by
from pyNastran.bdf.caseControlDeck import CaseControlDeck

testPath = pyNastran.__path__[0]
#print("testPath = %s" %(testPath))
from pyNastran.bdf.test.test_bdf import runBDF, run_all_files_in_folder


class Tester(unittest.TestCase):

    def runBDF(self, folder, bdfFilename, xref=False, cid=None,
               meshForm='combined', debug=False):
        cid = 0
        #xref = False
        runBDF(folder, bdfFilename, xref=xref, cid=cid, isFolder=True,
               meshForm=meshForm, debug=debug)

    def runAllFilesInFolder(self, folder, xref=False, cid=None, debug=False):
        run_all_files_in_folder(folder, xref=xref, cid=cid, debug=debug)


class BDF_Test(Tester):
    def test_bdf_01(self):
        bdfFilename = 'solidBending.bdf'
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        self.runBDF(folder, bdfFilename)
        self.runBDF(folder, bdfFilename, xref=True)

    def test_bdf_02(self):
        bdfFilename = 'plate_py.dat'
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        self.runBDF(folder, bdfFilename)
        self.runBDF(folder, bdfFilename, xref=True)

    def test_bdf_03(self):
        bdfFilename = 'beam_modes.dat'
        folder = os.path.abspath(os.path.join(testPath, '..', 'models'))
        self.runBDF(folder, bdfFilename)
        #self.runBDF(folder, bdfFilename, xref=True) ## PBEAML is not supported

    def test_bdf_04(self):
        bdfFilename = 'testA.bdf'
        folder = os.path.abspath(os.path.join(testPath, 'bdf', 'test', 'unit'))
        self.runBDF(folder, bdfFilename)
        #self.runBDF(folder, bdfFilename, xref=True) ## PBEAML is not supported

class BaseCard_Test(Tester):
    def test_base_card_01_collapse_thru(self):
        """
        tests collapse_thru method used by SETx cards
        """
        data = [1, 2, 3, 4, 5, 10]
        expected = [1, u'THRU', 5, 10]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 17]
        expected = [1, 3, 4, 5, 6, 17]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 7, 17]
        expected = [1, 3, 4, 'THRU', 7, 17]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [1, 3, 4, 6, 8, 10, 12, 14, 17]
        expected = [1, 3, 4, 'THRU', 14, 'BY', 2, 17]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 101]
        expected = [1, 3, 4, 5, 6, 8, 'THRU', 22, 'BY', 2, 101]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [1, 2, 3, 4, 5]
        expected = [1, 'THRU', 5]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [5]
        expected = [5]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [1,2,3,4,5, 7,9,11, 12,14,16]
        expected = [1, 'THRU', 5,
                    7, 9, 11,
                    12, 14, 16]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [1,2]
        expected = [1, 2]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [1,3,5,7,9,11]
        expected = [1, 'THRU', 11, 'BY', 2]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [1,2,3,4]
        expected = [1, 'THRU', 4]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [1,2,3]
        expected = [1, 2, 3]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

        data = [1, 2, 3, 4, 5, 6, 7, 8]
        expected = [1, 'THRU', 8]
        self.assertEquals(collapse_thru_by(data),expected,
                          collapse_thru_by(data))

class CaseControlTest(unittest.TestCase):
    def test_case_control_01(self):
        lines = ['SPC=2',
                 'MPC =3',
                 'STRESS= ALL',
                 'DISPLACEMENT(PLOT,PUNCH) = 8', ]

        deck = CaseControlDeck(lines)
        self.assertTrue(deck.has_parameter(0, 'SPC'))
        self.assertTrue(deck.has_parameter(0, 'sPC'))
        self.assertFalse(deck.has_parameter(0, 'JUNK'))
        #print("get_subcase_parameter(MPC) 3 = ", deck.get_subcase_parameter(
        #    0, 'MPC'))

        deck.add_parameter_to_global_subcase('GPFORCE = 7')

        deck.create_new_subcase(1)
        deck.create_new_subcase(2)

        deck.add_parameter_to_local_subcase(1, 'STRAIN = 7')

        out = deck.get_subcase_parameter(0, 'GPFORCE')

        deck.add_parameter_to_local_subcase(1, 'ANALYSIS = SAERO')
        deck.add_parameter_to_local_subcase(2, 'ANALYSIS = STATIC')
        out = deck.get_subcase_parameter(2, 'ANALYSIS')

        deck.add_parameter_to_local_subcase(1, 'SET 1 = 100')
        deck.add_parameter_to_local_subcase(1, 'SET 2 = 200')

        lines = ['DISPLACEMENT(PLOT,PUNCH) = 8',
                 'GPFORCE = 7',
                 'MPC = 3',
                 'SPC = 2',
                 'STRESS = ALL',
                 'SUBCASE 1',
                 '    ANALYSIS = SAERO',
                 '    STRAIN = 7',
                 '    SET 1 = 100',
                 '    SET 2 = 200',
                 'SUBCASE 2',
                 '    ANALYSIS = STATIC',]
        deck_string = '%s' %(deck)
        deck_lines = deck_string.strip().splitlines()
        self.assertEqual(lines, deck_lines)


if __name__ == '__main__':
    unittest.main()
