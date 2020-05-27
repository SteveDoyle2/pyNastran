"""tests the CaseControlDeck"""
import os
import unittest

import pyNastran
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.bdf.bdf_interface.subcase_utils import write_set, collapse_thru_packs

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(PKG_PATH, 'bdf', 'test')


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

        unused_out = deck.get_subcase_parameter(0, 'GPFORCE')

        deck.add_parameter_to_local_subcase(1, 'ANALYSIS = SAERO')
        deck.add_parameter_to_local_subcase(2, 'ANALYSIS = STATIC')
        unused_out = deck.get_subcase_parameter(2, 'ANALYSIS')

        deck.add_parameter_to_local_subcase(1, 'SET 1 = 100')
        deck.add_parameter_to_local_subcase(1, 'SET 2 = 200')

        lines_expected = [
            'DISPLACEMENT(PLOT,PUNCH) = 8',
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
        self.assertEqual(lines_expected, deck_lines, msg=msg)
        compare_lines(self, deck_lines, lines_expected, has_endline=False)

    def test_case_control_02(self):
        bdf_filename = os.path.join(TEST_PATH, 'unit', 'case_control.dat')
        bdf_filename2 = os.path.join(TEST_PATH, 'unit', 'case_control_out.dat')

        mesh = BDF(debug=False, log=None)
        mesh.read_bdf(bdf_filename, xref=True)
        str(mesh.case_control_deck)

        mesh.case_control_deck.create_new_subcase(1)

        #with self.assertRaises(AssertionError):
        str(mesh.case_control_deck)
        subcase1 = mesh.case_control_deck.subcases[1]
        str(subcase1)

        mesh.case_control_deck.add_parameter_to_local_subcase(1, 'LOAD=1')
        str(mesh.case_control_deck)

        mesh.case_control_deck.create_new_subcase(2)
        mesh.case_control_deck.add_parameter_to_local_subcase(2, 'LOAD=2')
        mesh.write_bdf(bdf_filename2)
        #print("---cc 3---\n%s" % str(mesh.case_control_deck))

        with open(bdf_filename2, 'r') as bdf_file:
            lines = bdf_file.readlines()

        lines_expected = [
            '$pyNastran: version=msc',
            '$pyNastran: punch=False',
            '$pyNastran: encoding=utf-8',
            '$pyNastran: nnodes=1',
            '$pyNastran: nelements=0',
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
        values = [
            11, 12, 13, 14, 15, 16, 17, 18,
            19, 20, 21, 22, 23, 24, 25, 26,
            1000000000000000000000000000000000000000000000000000000, 33]
        spaces = '    '
        set_id = 10

        singles, doubles = collapse_thru_packs(values)
        assert singles == [33, 1000000000000000000000000000000000000000000000000000000], singles
        assert doubles == [[11, 'THRU', 26]], doubles
        unused_msg = write_set(set_id, values, spaces)

        set_id = 11
        values = ['ALL']
        unused_msg = write_set(set_id, values, spaces)

    def test_case_control_04(self):
        seti = 'SET 88 = 5, 6, 7, 8, 9, 10 THRU 55 EXCEPT 15, 16, 77, 78, 79, 100 THRU 300'
        lines = []
        deck = CaseControlDeck(lines)
        deck.create_new_subcase(2)
        deck.add_parameter_to_local_subcase(2, seti)
        values, set_id = deck.get_subcase_parameter(2, 'SET 88')

        check = [
            (7, True),
            (13, True),
            (15, False),
            (16, False),
            (55, True),
            (77, True),
            (99, False),
            (150, True),
        ]
        for value, exists in check:
            if exists:
                assert value in values, 'value=%s should be in values=%s' % (value, values)
            else:
                assert value not in values, 'value=%s should NOT be in values=%s' % (value, values)

        unused_msg = write_set(set_id, values)

        singles, doubles = collapse_thru_packs(values)
        assert singles == [77, 78, 79], singles
        assert doubles == [[5, 'THRU', 14], [17, 'THRU', 55], [100, 'THRU', 300]], doubles


    def test_case_control_05(self):
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
            #'    SET = ALL',
            '    SET 20 = ALL',
            '    SPC = 42',
            '    TSTEPNL = 22',
            '    VELOCITY(PLOT,PRINT,PHASE) = ALL',
            'BEGIN BULK',
            ]
        deck = CaseControlDeck(lines)
        deck.create_new_subcase(2)
        deck.add_parameter_to_local_subcase(
            2, 'SET 2 = 11,12,13,14,15,16,17,18,'
            '19,20,21,22,23,24,25,26,'
            '1000000000000000000000000000000000000000000000000000000,33')
        deck_msg = '%s' % deck
        #print('%r' % deck_lines)
        deck_lines = deck_msg.split('\n')

        lines_expected = [
            'SUBCASE 1',
            #'    SET = ALL',
            '    SET 20 = ALL',
            '    SET 88 = 77, 78, 79, 5 THRU 14, 17 THRU 55, 100 THRU 300',
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
            '    SET 2 = 33,',
            '            1000000000000000000000000000000000000000000000000000000,',
            '            11 THRU 26',
        ]
        compare_lines(self, deck_lines, lines_expected, has_endline=False)

    def test_case_control_06(self):
        lines_expected = [
            'ACCELERATION(PLOT,PRINT,PHASE) = ALL',
            'DISPLACEMENT(PLOT,PRINT,PHASE) = ALL',
            'BEGIN BULK',
        ]
        deck = CaseControlDeck(lines_expected)
        deck_msg = '%s' % deck
        deck_lines = deck_msg.split('\n')
        compare_lines(self, deck_lines, lines_expected, has_endline=False)

    def test_case_control_07(self):
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
            'dsaprt=(end=sens)',
            'ECHO = UNSORT',
            'METHOD = 1',
            'OLOAD = ALL',
            'TITLE = VIBRATION OF A BEAM.',
            'SUBCASE 1',
            '    SUPER 1',
            '    DESSUB = 1',
            'SUBCASE 2',
            '    ANALYSIS = MODES',
            '    DESOBJ = 1',
            '    DESSUB = 2',
            '    DISPLACEMENT = ALL',
            '    dsaprt=(end=sens)',
            '    ECHO = UNSORT',
            '    METHOD = 1',
            '    OLOAD = ALL',
            '    TITLE = VIBRATION OF A BEAM.',
            'BEGIN BULK',
        ]
        deck = CaseControlDeck(lines)
        deck_msg = '%s' % deck
        #print('%s' % deck_msg)
        deck_lines = deck_msg.split('\n')
        #print "--------"
        #print deck_msg
        compare_lines(self, deck_lines, lines_expected, has_endline=False)

    def test_case_control_08(self):
        lines_expected = [
            '$pyNastran: version=msc\n',
            '$pyNastran: punch=True\n',
            '$pyNastran: encoding=utf-8\n',
            '$NODES\n',
            'GRID,100000,,43.91715,-29.,.8712984\n',
        ]
        bdf_filename = 'test7.bdf'
        bdf_filename2 = 'test7_bad.bdf'
        with open(bdf_filename, 'w', encoding='ascii') as bdf_file:
            for line in lines_expected:
                bdf_file.write(line)
        bdf = BDF(debug=False)
        bdf.read_bdf(bdf_filename)
        bdf.write_bdf(bdf_filename2)

        with open(bdf_filename, 'r', encoding='ascii') as bdf_file:
            lines = bdf_file.readlines()
            compare_lines(self, lines, lines_expected, has_endline=True)
        os.remove(bdf_filename)
        os.remove(bdf_filename2)

    def test_case_control_09(self):
        lines = [
            'GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),',
            'DATAREC=NO)=YES',
            'SUBCASE 1',
            '    DISPLACEMENT = ALL',
        ]
        lines_expected = [
            'GROUNDCHECK(PRINT, SET=(G,N,N+AUTOSPC,F,A), DATAREC=NO) = YES',
            'SUBCASE 1',
            '    DISPLACEMENT = ALL',
        ]
        deck = CaseControlDeck(lines)
        deck_msg = '%s' % deck
        #print('%s' % deck_msg)
        deck_lines = deck_msg.split('\n')
        compare_lines(self, deck_lines, lines_expected, has_endline=False)

    def test_case_control_10(self):
        lines = [
            'GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),',
            'THRESH=1e-2,DATAREC=NO)=YES',
            'SUBCASE 1',
            '    DISPLACEMENT = ALL',
        ]
        lines_expected = [
            'GROUNDCHECK(PRINT, SET=(G,N,N+AUTOSPC,F,A), THRESH=0.01, DATAREC=NO) = YES',
            'SUBCASE 1',
            '    DISPLACEMENT = ALL',
            'BEGIN BULK',
        ]
        deck = CaseControlDeck(lines)
        deck_msg = '%s' % deck
        #print('%s' % deck_msg)
        deck_lines = deck_msg.split('\n')
        compare_lines(self, deck_lines, lines_expected, has_endline=False)

    def test_case_control_11(self):
        """
        this test checks that subcase 3 uses its local definition of
        set 100 and subcase 4 uses the default definition
        """
        lines = [
            'SET 100 = 100',
            'DISP = 100',
            'SUBCASE 1',
            '     SPC = 1',
            '     LOAD = 1',
            'SUBCASE 2',
            '     SPC = 2',
            '     LOAD = 2',
            '     DISP = ALL',
            'SUBCASE 3',
            '     SET 100 = 100, 101',
            '     SPC = 3',
            '     LOAD = 3',
            '     DISP = 100',
            'SUBCASE 4',
            '     SPC = 3',
            '     LOAD = 3',
            '     DISP = 100',
        ]

        deck = CaseControlDeck(lines)

        default = deck.subcases[0]
        sc3 = deck.subcases[3]
        sc4 = deck.subcases[4]

        assert default.params['SET 100'] == [[100], 100, 'SET-type']
        assert sc3.params['SET 100'] == [[100, 101], 100, 'SET-type']
        assert sc4.params['SET 100'] == [[100], 100, 'SET-type']


    def test_groundcheck_1(self):
        """tests GROUNDCHECk"""
        lines = [
            'GROUNDCHECK(PRINT,',
            'THRESH=1e-2,DATAREC=NO)=YES',
            'SUBCASE 1',
            '    DISPLACEMENT = ALL',
        ]
        lines_expected = [
            'GROUNDCHECK(PRINT, THRESH=0.01, DATAREC=NO) = YES',
            'SUBCASE 1',
            '    DISPLACEMENT = ALL',
        ]
        deck = CaseControlDeck(lines)
        deck_msg = '%s' % deck
        #print('%s' % deck_msg)
        deck_lines = deck_msg.split('\n')
        compare_lines(self, deck_lines, lines_expected, has_endline=False)

    def test_groundcheck_2(self):
        lines = [
            'GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),THRESH=1e-2,DATAREC=NO,FAKE)=YES',
            'SUBCASE 1',
            '    DISPLACEMENT = ALL',
        ]
        with self.assertRaises(KeyError):  # FAKE is not a key
            unused_deck = CaseControlDeck(lines)

    def test_groundcheck_3(self):
        lines = [
            'GROUNDCHECK(SET=ALL)=YES',
            'SUBCASE 1',
            '    DISPLACEMENT = ALL',
        ]
        lines_expected = [
            'GROUNDCHECK(SET=ALL) = YES',
            'SUBCASE 1',
            '    DISPLACEMENT = ALL',
            'BEGIN BULK',
        ]
        deck = CaseControlDeck(lines)
        deck_msg = '%s' % deck
        #print('%s' % deck_msg)
        deck_lines = deck_msg.split('\n')
        compare_lines(self, deck_lines, lines_expected, has_endline=False)

        lines = [
            'GROUNDCHECK(PRINT, THRESH=0.01, DATAREC=NO, SET=ALL)=YES',
            'SUBCASE 1',
            '    DISPLACEMENT = ALL',
        ]
        lines_expected = [
            'GROUNDCHECK(PRINT, THRESH=0.01, DATAREC=NO, SET=ALL) = YES',
            'SUBCASE 1',
            '    DISPLACEMENT = ALL',
            'BEGIN BULK',
        ]
        deck = CaseControlDeck(lines)
        deck_msg = '%s' % deck
        #print('%s' % deck_msg)
        deck_lines = deck_msg.split('\n')
        compare_lines(self, deck_lines, lines_expected, has_endline=False)

    def test_weightcheck(self):
        weightchecks = [
            'WEIGHTCHECK=YES',
            'WEIGHTCHECK(GRID=12,SET=(G,N,A),MASS)=YES',
        ]
        for weightcheck in weightchecks:
            deck = CaseControlDeck([weightcheck])
            str(deck)

    def test_groundcheck(self):
        groundchecks = [
            'groundcheck=YES',
            'groundcheck(set=(g),datarec=yes,rthresh=.8)=yes',
        ]
        for groundcheck in groundchecks:
            deck = CaseControlDeck([groundcheck])
            str(deck)

    def test_sets(self):
        lines = [
            'SET 1 = 1,2, 21,THRU,25, 31,THRU,37',
            'SET 11 = 11, 12, 121 THRU 125, 131 THRU 137',
        ]
        lines_expected = [
            'SET 1 = 1, 2, 21 THRU 25, 31 THRU 37',
            'SET 11 = 11, 12, 121 THRU 125, 131 THRU 137',
            'BEGIN BULK',
        ]
        deck = CaseControlDeck(lines)
        deck_msg = '%s' % deck
        #print('%s' % deck_msg)
        deck_lines = deck_msg.split('\n')
        compare_lines(self, deck_lines, lines_expected, has_endline=False)

    def test_subcase_equals(self):
        lines = [
            'SET 100 = 100',
            'DISP = 100',
            'SUBCASE 1',
            '     SPC = 1',
            '     LOAD = 1',
            'SUBCASE= 2',
            '     SPC = 2',
            '     LOAD = 2',
            '     DISP = ALL',
            'SUBCASE =3',
            '     SET 100 = 100, 101',
            '     SPC = 3',
            '     LOAD = 3',
            '     DISP = 100',
            'SUBCASE=4',
            '     SPC = 3',
            '     LOAD = 3',
            '     DISP = 100',
        ]

        deck = CaseControlDeck(lines)
        unused_default = deck.subcases[0]
        unused_sc3 = deck.subcases[3]
        unused_sc4 = deck.subcases[4]

        #-------------------------------------
        lines = [
            'DISP = 100',
            'SUBCASEBAD',
            '     SPC = 3',
            '     LOAD = 3',
            '     DISP = 100',
        ]
        with self.assertRaises(SyntaxError):
            deck = CaseControlDeck(lines)

        #-------------------------------------
        lines = [
            'DISP = 100',
            'SUBCASE BAD',
            '     SPC = 3',
            '     LOAD = 3',
            '     DISP = 100',
        ]
        with self.assertRaises(ValueError):
            deck = CaseControlDeck(lines)

    def test_case_control_set(self):
        """tests adding SETs"""
        deck = CaseControlDeck([], log=None)
        subcase = deck.create_new_subcase(1)
        subcase.params['SET 42'] = [[1, 2, 3, 4, 5], 42, 'SET-type']
        set_str = write_set(41, [1, 2, 3]).rstrip()
        set_str2 = set_str.split('=', 1)[1]
        assert set_str2 == ' 1, 2, 3', set_str2
        subcase.add('SET 52', set_str2, 52, 'SET-type')
        subcase.add_set_from_values(100, [1, 2, 10])
        # subcase.add('SET', value, options, param_type)  # TODO: doesn't work

def compare_lines(self, lines, lines_expected, has_endline):
    i = 0
    for line, line_expected in zip(lines, lines_expected):
        line = line.rstrip()
        line_expected = line_expected.rstrip()
        msg = 'The lines are not the same...i=%s\n' % i
        msg += 'line     = %r\n' % line
        msg += 'expected = %r\n' % line_expected
        if has_endline:
            msg += '-------------\n--Actual--\n%s' % ''.join(lines)
            msg += '-------------\n--Expected--\n%s' % ''.join(lines_expected)
        else:
            msg += '-------------\n--Actual--\n%s' % '\n'.join(lines)
            msg += '-------------\n--Expected--\n%s' % '\n'.join(lines_expected)
        self.assertEqual(line, line_expected, msg)
        i += 1


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
