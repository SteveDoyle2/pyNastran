from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import unittest
from numpy import allclose, array
from numpy.linalg import norm

import pyNastran
from pyNastran.bdf.cards.baseCard import collapse_thru_by
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.caseControlDeck import CaseControlDeck

pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'bdf', 'test')
#print("testPath = %s" % testPath)
from pyNastran.bdf.test.test_bdf import run_bdf, run_all_files_in_folder


class Tester(unittest.TestCase):

    def run_bdf(self, folder, bdfFilename, xref=False, cid=None,
                meshForm='combined', debug=False):
        cid = 0
        #xref = False
        return run_bdf(folder, bdfFilename, xref=xref, cid=cid, isFolder=True,
                       meshForm=meshForm, debug=debug)

    def run_all_files_in_folder(self, folder, xref=False, cid=None, debug=False):
        run_all_files_in_folder(folder, xref=xref, cid=cid, debug=debug)


class TestBDF(Tester):
    def test_bdf_01(self):
        bdfFilename = os.path.join('solid_bending', 'solid_bending.bdf')
        folder = os.path.abspath(os.path.join(pkg_path, '..', 'models'))
        self.run_bdf(folder, bdfFilename)
        (fem1, fem2, diffCards2) = self.run_bdf(folder, bdfFilename, xref=True)

        for fem in [fem1, fem2]:
            assert len(fem.params) == 2, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 72, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 1, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 186, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 0, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 1, 'len(properties) = %i' % len(fem.properties)
        mass, cg, I = fem1.mass_properties()

        assert allclose(mass, 6.0), 'mass = %s' % mass
        cg_exact = array([0.5, 1., 1.5])
        for i, (cgi, cgie) in enumerate(zip(cg, cg_exact)):
            assert allclose(cgi, cgie), 'i=%s cg=%s' % (i, str(cg))

        self._compare_mass_cg_I(fem1)
        self._compare_mass_cg_I(fem1, reference_point=u'cg')
        mass, cg, I = fem1.mass_properties(reference_point='cg')

    def _compare_mass_cg_I(self, fem1, reference_point=None, sym_axis=None):
        num_cpus = 4
        mass1, cg1, I1 = fem1.mass_properties(reference_point=reference_point, sym_axis=sym_axis)
        mass2, cg2, I2 = fem1.mass_properties(reference_point=reference_point, sym_axis=sym_axis, num_cpus=num_cpus)

        assert mass1 == mass2, 'mass1=%s mass2=%s' % (mass1, mass2)
        assert allclose(norm((cg1 - cg2)**2), 0.0), 'cg1-cg2=%s' % (cg1 - cg2)
        assert allclose(norm((I1  -  I2)**2), 0.0), 'I1-I2=%s' % (I1 - I2)

    def test_bdf_02(self):
        bdfFilename = os.path.join('plate_py', 'plate_py.dat')
        folder = os.path.abspath(os.path.join(pkg_path, '..', 'models'))
        self.run_bdf(folder, bdfFilename)
        (fem1, fem2, diffCards2) = self.run_bdf(folder, bdfFilename, xref=True)

        for fem in [fem1, fem2]:
            assert len(fem.coords) == 3, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.params) == 6, 'len(params) = %i' % len(fem.params)
            assert len(fem.nodes) == 231, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 1, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 200, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 1, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 1, 'len(properties) = %i' % len(fem.properties)
        self._compare_mass_cg_I(fem1)
        self._compare_mass_cg_I(fem1, reference_point=u'cg')


    def test_bdf_04(self):
        bdfFilename = os.path.join('cbush', 'cbush.dat')
        folder = os.path.abspath(os.path.join(pkg_path, '..', 'models'))
        (fem1, fem2, diffCards2) = self.run_bdf(folder, bdfFilename)

        for fem in [fem1, fem2]:
            assert len(fem.params) == 6, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 2, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 0, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 1, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 0, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 1, 'len(properties) = %i' % len(fem.properties)  # PBEAML issue

        self._compare_mass_cg_I(fem1)
        self._compare_mass_cg_I(fem1, reference_point=u'cg')

        self.run_bdf(folder, bdfFilename, xref=True)

    def test_bdf_03(self):
        bdfFilename = os.path.join('beam_modes', 'beam_modes.dat')
        folder = os.path.abspath(os.path.join(pkg_path, '..', 'models'))
        (fem1, fem2, diffCards2) = self.run_bdf(folder, bdfFilename)

        for fem in [fem1, fem2]:
            assert len(fem.params) == 6, 'len(params) = %i' % len(fem.params)
            assert len(fem.coords) == 1, 'len(coords) = %i' % len(fem.coords)
            assert len(fem.nodes) == 12, 'len(nodes) = %i' % len(fem.nodes)
            assert len(fem.materials) == 1, 'len(materials) = %i' % len(fem.materials)
            assert len(fem.elements) == 11, 'len(elements) = %i' % len(fem.elements)
            assert len(fem.methods) == 1, 'len(methods) = %i' % len(fem.methods)
            assert len(fem.properties) == 3, 'len(properties) = %i' % len(fem.properties)  # PBEAML issue
        self._compare_mass_cg_I(fem1)
        #self._compare_mass_cg_I(fem1, reference_point=u'cg')

        #self.run_bdf(folder, bdfFilename, xref=True) # PBEAML is not supported

    def test_bdf_05(self):
        bdfFilename = 'testA.bdf'
        folder = os.path.abspath(os.path.join(pkg_path, 'bdf', 'test', 'unit'))
        self.run_bdf(folder, bdfFilename)
        #self.run_bdf(folder, bdfFilename, xref=True) # PBEAML is not supported

class BaseCard_Test(Tester):
    def test_base_card_01_collapse_thru(self):
        """
        tests collapse_thru method used by SETx cards
        """
        data = [1, 2, 3, 4, 5, 10]
        expected = [1, u'THRU', 5, 10]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 17]
        expected = [1, 3, 4, 5, 6, 17]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 7, 17]
        expected = [1, 3, 4, 'THRU', 7, 17]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [1, 3, 4, 6, 8, 10, 12, 14, 17]
        expected = [1, 3, 4, 'THRU', 14, 'BY', 2, 17]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [1, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 101]
        expected = [1, 3, 4, 5, 6, 8, 'THRU', 22, 'BY', 2, 101]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [1, 2, 3, 4, 5]
        expected = [1, 'THRU', 5]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [5]
        expected = [5]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [1,2,3,4,5, 7,9,11, 12,14,16]
        expected = [1, 'THRU', 5,
                    7, 9, 11,
                    12, 14, 16]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [1,2]
        expected = [1, 2]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [1,3,5,7,9,11]
        expected = [1, 'THRU', 11, 'BY', 2]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [1,2,3,4]
        expected = [1, 'THRU', 4]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [1,2,3]
        expected = [1, 2, 3]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

        data = [1, 2, 3, 4, 5, 6, 7, 8]
        expected = [1, 'THRU', 8]
        self.assertEqual(collapse_thru_by(data),expected,
                         collapse_thru_by(data))

class CaseControlTest(unittest.TestCase):
    def test_case_control_01(self):
        lines = ['SPC=2',
                 'MPC =3',
                 'STRESS= ALL',
                 'DISPLACEMENT(PLOT,PUNCH) = 8', ]

        deck = CaseControlDeck(lines)
        deck.write_begin_bulk = False
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
                 '    SET 1 = 100',
                 '    SET 2 = 200',
                 '    ANALYSIS = SAERO',
                 '    STRAIN = 7',
                 'SUBCASE 2',
                 '    ANALYSIS = STATIC',]
        deck_string = '%s' % deck
        deck_lines = deck_string.strip().splitlines()
        self.assertEqual(lines, deck_lines)

    def test_case_control_02(self):
        bdf_filename = os.path.join(test_path, 'unit', 'case_control.dat')
        bdf_filename2 = os.path.join(test_path, 'unit', 'case_control_out.dat')

        mesh = BDF(debug=True,log=None)
        mesh.readBDF(bdf_filename, includeDir=None, xref=True)
        str(mesh.caseControlDeck)

        mesh.caseControlDeck.create_new_subcase(1)

        #with self.assertRaises(AssertionError):
        str(mesh.caseControlDeck)
        subcase1 = mesh.caseControlDeck.subcases[1]
        str(subcase1)

        mesh.caseControlDeck.add_parameter_to_local_subcase(1, 'LOAD=1')
        str(mesh.caseControlDeck)

        mesh.caseControlDeck.create_new_subcase(2)
        mesh.caseControlDeck.add_parameter_to_local_subcase(2, 'LOAD=2')
        mesh.write_bdf(bdf_filename2)
        #print("---cc 3---\n%s" % str(mesh.caseControlDeck))

        f = open(bdf_filename2, 'r')
        lines = f.readlines()
        f.close()

        lines_expected = [
            '$EXECUTIVE CONTROL DECK',
            'SOL 101',
            'CEND',
            '$CASE CONTROL DECK',
            'TITLE = STATIC',
            'SUBCASE 1',
            '    LOAD = 1',
            'SUBCASE 2',
            '    LOAD = 2',
            'BEGIN BULK',
            '$PARAMS',
            'PARAM    AUTOSPC     YES',
            'PARAM     NOFISR       0',
            '$NODES',
            'GRID           1              0.      0.      0.',
            'ENDDATA',
        ]
        for line, line_expected in zip(lines, lines_expected):
            line = line.rstrip()
            msg = 'The lines are not the same...\n'
            msg += 'line     = %r\n' % line
            msg += 'expected = %r' % line_expected
            self.assertEqual(line, line_expected, msg)

    def test_case_control_03(self):
        lines = [
            'SUBCASE 1',
            '    ACCELERATION(PLOT,PRINT,PHASE) = ALL',
            '    DISPLACEMENT(PLOT,PRINT,PHASE) = ALL',
            '    DLOAD = 32',
            '    M2GG = 111',
            '    SET 88  = 5, 6, 7, 8, 9, 10 THRU 55 EXCEPT 15, 16, 77, 78, 79, '
            '100 THRU 300',
            '    SET 99  = 1 THRU 10',
            '    SET 105 = 1.009, 10.2, 13.4, 14.0, 15.0',
            '    SET 111 = MAAX1,MAAX2',
            '    SET 1001 = 101/T1, 501/T3, 991/R3',
            '    SET = ALL',
            '    SPC = 42',
            '    TSTEPNL = 22',
            '    VELOCITY(PLOT,PRINT,PHASE) = ALL',
            'BEGIN BULK',
            ]
        deck = CaseControlDeck(lines)
        deck.create_new_subcase(2)
        deck.add_parameter_to_local_subcase(2, 'SET 2 = 11,12,13,14,15,16,17,18,'
           '19,20,21,22,23,24,25,26,'
           '1000000000000000000000000000000000000000000000000000000,33')
        deck_msg = '%s' % deck
        #print('%r' % deck_lines)
        deck_lines = deck_msg.split('\n')

        lines_expected = [
            'SUBCASE 1',
            '    SET = ALL',
            '    SET 88 = 5, 6, 7, 8, 9, 10 THRU 55 EXCEPT 15, 16, 77, 78, 79,',
            '             100 THRU 300',
            '    SET 99 = 1 THRU 10',
            '    SET 105 = 1.009, 10.2, 13.4, 14.0, 15.0',
            '    SET 111 = MAAX1, MAAX2',
            '    SET 1001 = 101/T1, 501/T3, 991/R3',
            '    ACCELERATION(PLOT,PRINT,PHASE) = ALL',
            '    DISPLACEMENT(PLOT,PRINT,PHASE) = ALL',
            '    DLOAD = 32',
            '    M2GG = 111',
            '    SPC = 42',
            '    TSTEPNL = 22',
            '    VELOCITY(PLOT,PRINT,PHASE) = ALL',
            'SUBCASE 2',
            '    SET 2 = 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,',
            '            25, 26,',
            '            1000000000000000000000000000000000000000000000000000000,',
            '            33',
        ]
        for line, line_expected in zip(deck_lines, lines_expected):
            line = line.rstrip()
            msg = 'The lines are not the same...\n'
            msg += 'line     = %r\n' % line
            msg += 'expected = %r' % line_expected
            self.assertEqual(line, line_expected, msg)
        #print('%s' % deck)

    def test_case_control_04(self):
        lines_expected = [
            'ACCELERATION(PLOT,PRINT,PHASE) = ALL',
            'DISPLACEMENT(PLOT,PRINT,PHASE) = ALL',
            'BEGIN BULK',
        ]
        deck = CaseControlDeck(lines_expected)
        deck_msg = '%s' % deck
        deck_lines = deck_msg.split('\n')
        for line, line_expected in zip(deck_lines, lines_expected):
            line = line.rstrip()
            msg = 'The lines are not the same...\n'
            msg += 'line     = %r\n' % line
            msg += 'expected = %r' % line_expected
            self.assertEqual(line, line_expected, msg)
        #print('%s' % deck)


    def test_case_control_04(self):
        lines = [
            'TITLE= VIBRATION OF A BEAM.',
            'dsaprt=(end=sens)',
            'ECHO = UNSORT',
            'OLOAD = ALL',
            'DISP = ALL',
            'DESSUB = 2',
            'METHOD = 1',
            'ANALYSIS = MODES',
            'DESOBJ = 1',
            'SUBCASE 1',
            'DESSUB = 1',
            ' SUPER 1',
            'SUBCASE 2',
            'BEGIN BULK',
        ]
        lines_expected = [
            'ANALYSIS = MODES',
            'DESOBJ = 1',
            'DESSUB = 2',
            'DISPLACEMENT = ALL',
            'ECHO = UNSORT',
            'METHOD = 1',
            'OLOAD = ALL',
            'TITLE = VIBRATION OF A BEAM.',
            'dsaprt=(end=sens)',
            'SUBCASE 1',
            '    SUPER 1',
            '    DESSUB = 1',
            'SUBCASE 2',
            '    ANALYSIS = MODES',
            '    DESOBJ = 1',
            '    DESSUB = 2',
            '    DISPLACEMENT = ALL',
            '    ECHO = UNSORT',
            '    METHOD = 1',
            '    OLOAD = ALL',
            '    TITLE = VIBRATION OF A BEAM.',
            '    dsaprt=(end=sens)',
            'BEGIN BULK',
        ]
        deck = CaseControlDeck(lines)
        deck_msg = '%s' % deck
        #print('%s' % deck_msg)
        deck_lines = deck_msg.split('\n')
        for line, line_expected in zip(deck_lines, lines_expected):
            line = line.rstrip()
            msg = 'The lines are not the same...\n'
            msg += 'line     = %r\n' % line
            msg += 'expected = %r' % line_expected
            self.assertEqual(line, line_expected, msg)


if __name__ == '__main__':
    unittest.main()