import os
import unittest

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.caseControlDeck import CaseControlDeck

pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'bdf', 'test')


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
        #print("get_subcase_parameter(MPC) 3 = ", deck.get_subcase_parameter(0, 'MPC'))

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

        msg = '\n' + '\n'.join(deck_lines)
        self.assertEqual(lines, deck_lines, msg=msg)

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
            msg += '-------------\n--Actual--\n%s' % deck_msg
            msg += '-------------\n--Expected--\n%s' % '\n'.join(lines_expected)
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
            msg += '-------------\n--Actual--\n%s' % deck_msg
            msg += '-------------\n--Expected--\n%s' % '\n'.join(lines_expected)
            self.assertEqual(line, line_expected, msg)
        #print('%s' % deck)


    def test_case_control_05(self):
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
        #print "--------"
        #print deck_msg
        for line, line_expected in zip(deck_lines, lines_expected):
            line = line.rstrip()
            msg = 'The lines are not the same...\n'
            msg += 'line     = %r\n' % line
            msg += 'expected = %r\n' % line_expected
            msg += '-------------\n--Actual--\n%s' % deck_msg
            msg += '-------------\n--Expected--\n%s' % '\n'.join(lines_expected)
            self.assertEqual(line, line_expected, msg)

if __name__ == '__main__':
    unittest.main()