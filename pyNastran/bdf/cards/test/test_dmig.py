import unittest

import os
import pyNastran
from pyNastran.bdf.bdf import BDF, BDFCard, DMIG

from numpy import array, array_equal, sqrt, sin, cos, radians

root_path = pyNastran.__path__[0]
test_path = os.path.join(root_path, 'bdf', 'cards', 'test')

class TestDMIG(unittest.TestCase):

    def test_dmig_1(self):
        """
        Tests DMIG reading
        """
        model = BDF(debug=False)
        bdf_name = os.path.join(test_path, 'dmig.bdf')
        model.read_bdf(bdf_name, xref=False, punch=True)
        out = model.dmigs['REALS'].get_matrix(is_sparse=False)

        reals_actual, rows_reversed, cols_reversed = out
        #print "---reals_actual---\n", reals_actual
        #print "---out---\n", out

        reals_expected = [
            [1.0, 0.5, 0.25],
            [0.5, 2.0, 0.75],
            [0.25, 0.75, 3.0],
        ]
        a_matrix = model.dmigs['REALS']
        assert len(a_matrix.GCi) == 6, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 6, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)


        self.assertTrue(array_equal(reals_expected, reals_actual))

    def test_dmig_2(self):
        model = BDF(debug=False)
        bdf_name = os.path.join(test_path, 'dmig.bdf')
        model.read_bdf(bdf_name, xref=False, punch=True)

        out = model.dmigs['REAL'].get_matrix(is_sparse=False)
        REAL_actual, rows_reversed, cols_reversed = out
        #print "---REAL_actual---\n", REAL_actual
        REAL_expected = [
            [1.0, 0.5, 0.25],
            [0.0, 2.0, 0.75],
            [0.0, 0.0, 3.0],
        ]
        a_matrix = model.dmigs['REAL']
        assert len(a_matrix.GCi) == 6, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 6, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)

        self.assertTrue(array_equal(REAL_expected, REAL_actual))

        #model2 = BDF(debug=False)
        #bdf_name = os.path.join(test_path, 'include_dir', 'include.inc')
        #model2.read_bdf(bdf_name, xref=False, punch=True)

    def test_dmig_3(self):
        model = BDF(debug=False)
        bdf_name = os.path.join(test_path, 'dmig.bdf')
        model.read_bdf(bdf_name, xref=False, punch=True)
        out = model.dmigs['IMAG'].get_matrix(is_sparse=False)

        IMAG_actual, rows_reversed, cols_reversed = out
        #print "---IMAG_actual---\n", IMAG_actual
        IMAG_expected_real = [
            [1.0, 0.5, 0.25],
            [0.0, 2.0, 0.75],
            [0.0, 0.0, 3.0],
        ]
        IMAG_expected_imag = [
            [1.1, 0.51, 0.251],
            [0.0, 2.1, 0.751],
            [0.0, 0.0, 3.1],
        ]
        a_matrix = model.dmigs['IMAG']
        assert len(a_matrix.GCi) == 6, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 6, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)

        IMAG_expected = array(IMAG_expected_real) + array(IMAG_expected_imag)*1j
        self.assertTrue(array_equal(IMAG_expected, IMAG_actual))

    def test_dmig_4(self):
        model = BDF(debug=False)
        bdf_name = os.path.join(test_path, 'dmig.bdf')
        model.read_bdf(bdf_name, xref=False, punch=True)

        out = model.dmigs['IMAGS'].get_matrix(is_sparse=False)
        IMAGS_actual, rows_reversed, cols_reversed = out
        #print "---IMAGS_actual---\n", IMAGS_actual
        IMAGS_expected_real = [
            [1.0, 0.5, 0.25],
            [0.5, 2.0, 0.75],
            [0.25, 0.75, 3.0],
        ]
        IMAGS_expected_imag = [
            [1.1, 0.51, 0.251],
            [0.51, 2.1, 0.751],
            [0.251, 0.751, 3.1],
        ]
        a_matrix = model.dmigs['IMAGS']
        assert len(a_matrix.GCi) == 6, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 6, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)

        IMAGS_expected = array(IMAGS_expected_real) + array(IMAGS_expected_imag)*1j
        msg = '\n%s_actual\n%s\n\n----' % ('IMAGS', IMAGS_actual)
        msg += '\n%s_expected\n%s\n----' % ('IMAGS', IMAGS_expected)
        msg += '\n%s_delta\n%s\n----' % ('IMAGS', IMAGS_actual-IMAGS_expected)
        self.assertTrue(array_equal(IMAGS_expected, IMAGS_actual), msg)

    def test_dmig_5(self):
        model = BDF(debug=False)
        bdf_name = os.path.join(test_path, 'dmig.bdf')
        model.read_bdf(bdf_name, xref=False, punch=True)
        out = model.dmigs['POLE'].get_matrix(is_sparse=False)

        POLE_actual, rows_reversed, cols_reversed = out
        #print "---POLE_actual---\n", POLE_actual
        mag_expected = array([
            [1.0, 4.0, 5.0],
            [0.0, 2.0, 6.0],
            [0.0, 0.0, 3.0],
        ])

        a_matrix = model.dmigs['POLE']
        assert len(a_matrix.GCi) == 6, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 6, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)

        A_expected = mag_expected * cos(radians(45))
        B_expected = mag_expected * sin(radians(45))
        POLE_expected = A_expected + B_expected * 1j

        msg = '\n%s_actual\n%s\n\n----' % ('POLE', POLE_actual)
        msg += '\n%s_expected\n%s\n----' % ('POLE', POLE_expected)
        msg += '\n%s_delta\n%s\n----' % ('POLE', POLE_actual-POLE_expected)
        self.assertTrue(array_equal(POLE_expected, POLE_actual), msg)

    def test_dmig_06(self):
        lines = ['DMIG    ENFORCE 0       1       1       0']
        model = BDF(debug=False)
        card = model.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = DMIG(card)
        card.write_card(size, 'dummy')
        #card.rawFields()

    def test_dmig_07(self):
        cards = [
            ['DMIG, A, 0, 9, 1, 1,  ,    , 1'],
            ['DMIG, A, 1, 0,  , 2, 1, 1.0,'],
            ['DMIG, A, 1, 0,  , 2, 2, 1.0,'],
            ['DMIG, A, 1, 0,  , 2, 3, 1.0,'],
        ]
        model = BDF(debug=False)
        for card_lines in cards:
            model.add_card(card_lines, 'DMIG', is_list=False)
        model.fill_dmigs()

        a_matrix = model.dmigs['A']
        assert len(a_matrix.GCi) == 3, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 3, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)

    def test_dmig_08(self):
        cards = [
            ['DMIG, A, 1, 0,  , 2, 1, 1.0,'],
            ['DMIG, A, 0, 9, 1, 1,  ,    , 1'],
            ['DMIG, A, 1, 0,  , 2, 2, 1.0,'],
            ['DMIG, A, 1, 0,  , 2, 3, 1.0,'],
        ]
        model = BDF(debug=False)
        for card_lines in cards:
            model.add_card(card_lines, 'DMIG', is_list=False)
        model.fill_dmigs()

        a_matrix = model.dmigs['A']
        assert len(a_matrix.GCi) == 3, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 3, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)

    def test_dmig_09(self):
        cards = [
            ['DMIG, A, 0, 9, 1, 1,  ,    , 1'],
            ['DMIG, A, 1, ,  , 2, 1, 1.0,'],
            ['DMIG, A, 1, ,  , 2, 2, 1.0,'],
            ['DMIG, A, 1, ,  , 2, 3, 1.0,'],
        ]
        model = BDF(debug=False)
        for card_lines in cards:
            model.add_card(card_lines, 'DMIG', is_list=False)
        model.fill_dmigs()

        a_matrix = model.dmigs['A']
        assert len(a_matrix.GCi) == 3, 'len(GCi)=%s GCi=%s matrix=\n%s' % (len(a_matrix.GCi), a_matrix.GCi, a_matrix)
        assert len(a_matrix.GCj) == 3, 'len(GCj)=%s GCj=%s matrix=\n%s' % (len(a_matrix.GCj), a_matrix.GCj, a_matrix)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
