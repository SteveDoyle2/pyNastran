import unittest
import numpy as np
from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.cards.bdf_tables import (
    TABLED1, TABLED2, TABLED3, TABLED4,
    TABLEM1, TABLEM2, TABLEM3, TABLEM4,
    TABDMP1, #TABLES1, TABLEST, TABRND1, TABRNDG,
)
from pyNastran.bdf.cards.test.utils import save_load_deck

class TestTables(unittest.TestCase):

    def test_tabrndg(self):
        model = BDF(debug=False)
        tid = 101
        #1 : von Karman
        #2 : Dryden
        psd_type = 1
        LU = 2.17
        wg = 3.14
        table = model.add_tabrndg(tid, psd_type, LU, wg, comment='')
        str(table)
        table = model.add_tabrndg(tid+1, psd_type+1, LU, wg, comment='')
        str(table)
        save_load_deck(model, run_convert=False, run_quality=False, run_remove_unused=False)

    def test_tables1(self):
        model = BDF(debug=False)
        tid = 101
        #1 : von Karman
        #2 : Dryden
        Type = 1
        LU = 2.17
        wg = 3.14
        x = [1., 2., 3.]
        y = [10., 20., 30.]
        table = model.add_tables1(tid, x, y, Type=1, comment='tables1')
        str(table)
        tid = 110
        y = [10, 20, 30]
        table = model.add_tablest(tid, x, y, comment='tablest')
        str(table)
        save_load_deck(
            model, run_convert=False, run_quality=False,
            run_remove_unused=False,
            run_save_load_hdf5=False)

    def test_tabdmp1_01(self):
        model = BDF(debug=False)
        lines = ['TABDMP1,100,,,,,,,,+',
                 '+,1e-3,.02,200.,.02,ENDT',]
        card = model._process_card(lines)
        card = BDFCard(card)
        #print(card)
        table = TABDMP1.add_card(card, comment='table')
        table.raw_fields()
        str(table)

        table.write_card(size=8)
        table.write_card_16(is_double=False)
        table.write_card_16(is_double=True)

        assert table._comment == '$table\n', '%r' % table._comment
        assert table.comment == '$table\n', '%r' % table.comment
        fields = table.raw_fields()
        msg = table.write_card(size=8).rstrip()
        #print(msg)
        lines_expected = [
            '$table',
            'TABDMP1      100       G',
            '            .001     .02    200.     .02    ENDT']
            #'            1E-3    0.02   200.0    0.02    ENDT']
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

        msg = table.write_card(size=16).rstrip()
        #print(msg)
        lines_expected = [
            '$table',
            'TABDMP1*             100               G',
            '*',
            '*                   .001             .02            200.             .02',
            '*                   ENDT'
        ]
        lines_actual = [line.rstrip() for line in msg.rstrip().split('\n')]
        msg = '\n%r\n\n%r\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

        x = [1., 2., 3.]
        y = [10., 20., 30.]
        model.add_tabdmp1(100, x, y)
        save_load_deck(model, run_convert=False, run_quality=False, run_remove_unused=False)

    def test_tabled1(self):
        model = BDF(debug=False)
        lines = [
            'TABLED1, 32',
            ',-2.0, 6.0, 2.0, 5.0, 3.0, 5.6, ENDT',
        ]
        card = model._process_card(lines)
        card = BDFCard(card)
        table = TABLED1.add_card(card, comment='table')
        table.raw_fields()
        str(table)

        table.write_card(size=8)
        table.write_card_16(is_double=False)
        table.write_card_16(is_double=True)

        interp = table.interpolate(0.)
        assert np.allclose(interp, [5.5]), interp
        interp = table.interpolate([0., 5.])

        x = [1., 2., 3.]
        y = [10., 20., 30.]
        model.add_tabled1(101, x, y, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment='')
        model.add_tabled1(102, x, y, xaxis='LOG', yaxis='LOG', extrap=0, comment='')
        model.add_tabled1(103, x, y, xaxis='LINEAR', yaxis='LOG', extrap=1, comment='')
        save_load_deck(model, run_convert=False, run_quality=False, run_remove_unused=False)

    def test_tabled2(self):
        model = BDF(debug=False)
        lines = [
            'TABLED2, 15, -10.5',
            ',1.0, -4.5, 2.0, -4.2, 2.0, 2.8, 7.0, 6.5',
            ',SKIP, SKIP, 9.0, 6.5, ENDT',
        ]
        card = model._process_card(lines)
        card = BDFCard(card)
        table = TABLED2.add_card(card, comment='table')
        table.raw_fields()
        str(table)

        table.write_card(size=8)
        table.write_card_16(is_double=False)
        table.write_card_16(is_double=True)

        interp = table.interpolate(0.)
        assert np.allclose(interp, [-5.875]), interp
        interp = table.interpolate([0., 5.])

        x = [1., 2., 3.]
        y = [10., 20., 30.]
        x1 = 2.0
        model.add_tabled2(101, x1, x, y, extrap=0, comment='')
        x1 = 3.0
        model.add_tabled2(102, x1, x, y, extrap=1, comment='')
        save_load_deck(model, run_convert=False, run_quality=False, run_remove_unused=False)

    def test_tabled3(self):
        model = BDF(debug=False)
        lines = [
            'TABLED3, 62, 126.9, 30.0',
            ',2.9, 2.9, 3.6, 4.7, 5.2, 5.7, ENDT',
        ]
        card = model._process_card(lines)
        card = BDFCard(card)
        table = TABLED3.add_card(card, comment='table')
        table.raw_fields()
        str(table)

        #interp = table.interpolate(0.)
        #print('interp =', interp, type(interp))
        #assert np.allclose(interp, [5.5]), interp

        x = [1., 2., 3.]
        y = [10., 20., 30.]
        x1 = 2.0
        x2 = 3.0
        model.add_tabled3(101, x1, x2, x, y, extrap=0, comment='')
        x1 = 12.0
        x2 = 13.0
        model.add_tabled3(102, x1, x2, x, y, extrap=1, comment='')
        save_load_deck(model, run_convert=False, run_quality=False, run_remove_unused=False)

    def test_tabled4(self):
        model = BDF(debug=False)
        lines = [
            'TABLED4, 28, 0.0, 1.0, 0.0, 100.',
            ',2.91, -0.0329, 6.51-5, 0.0, -3.4-7, ENDT',
        ]
        card = model._process_card(lines)
        card = BDFCard(card)
        table = TABLED4.add_card(card)
        table.write_card(size=8)
        table.write_card_16(is_double=False)
        table.write_card_16(is_double=True)
        table.raw_fields()
        str(table)

        interp = table.interpolate(5.)
        assert np.allclose(interp, [2.746915]), interp
        with self.assertRaises(ValueError): # bug...
            interp = table.interpolate([0., 5.])

    #def test_tableht(self):
        #lines = [
            #'TABLEHT, 85',
            #'10.0, 101, 25.0, 102, 40.0, 110, ENDT',
        #]
        #card = model._process_card(lines)
        #card = BDFCard(card)
        #card2 = TABLEHT.add_card(card)
        x1 = 2.0
        x2 = 3.0
        x3 = 4.0
        x4 = 5.0
        a = [4., 5., 6.]
        model.add_tabled4(101, x1, x2, x3, x4, a, comment='tabled4')

        x1 = 12.0
        x2 = 13.0
        x3 = 14.0
        x4 = 15.0
        a = [40., 50., 60.]
        model.add_tabled4(102, x1, x2, x3, x4, a, comment='')
        save_load_deck(model, run_convert=False, run_quality=False, run_remove_unused=False)


    def test_tableh1(self):
        model = BDF(debug=False)
        lines = [
            'TABLEH1, 32',
            '-3.0, 6.9, 2.0, 5.6, 3.0, 5.6, ENDT',
        ]
        card = model._process_card(lines)
        card = BDFCard(card)
        table = TABLEM1.add_card(card, comment='table')
        table.raw_fields()
        str(table)

        table.write_card(size=8)
        table.write_card_16(is_double=False)
        table.write_card_16(is_double=True)

        x = [1., 2., 3.]
        y = [4., 5., 6.]
        model.add_tableh1(101, x, y, comment='tableh1')
        model.add_tableht(102, x, y, comment='tableht')
        save_load_deck(model)

    def test_tablem1(self):
        model = BDF(debug=False)
        lines = [
            'TABLEM1, 32',
            '-3.0, 6.9, 2.0, 5.6, 3.0, 5.6, ENDT',
        ]
        card = model._process_card(lines)
        card = BDFCard(card)
        table = TABLEM1.add_card(card, comment='table')
        table.raw_fields()
        str(table)

        x = [1., 2., 3.]
        y = [10., 20., 30.]
        model.add_tablem1(101, x, y, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment='')
        model.add_tablem1(102, x, y, xaxis='LOG', yaxis='LOG', extrap=0, comment='')
        model.add_tablem1(103, x, y, xaxis='LINEAR', yaxis='LOG', extrap=1, comment='')
        save_load_deck(model)

        #interp = table.interpolate(0.)
        #print('interp =', interp, type(interp))
        #assert np.allclose(interp, [5.5]), interp

    def test_tablem2(self):
        model = BDF(debug=False)
        lines = [
            'TABLEM2, 15, -10.5',
            ',1.0, -4.5, 2.0, -4.5, 2.0, 2.8, 7.0, 6.5',
            ',SKIP, SKIP, 9.0, 6.5, ENDT',
        ]
        card = model._process_card(lines)
        card = BDFCard(card)
        table = TABLEM2.add_card(card)
        table.raw_fields()
        str(table)
        table.write_card(size=8)
        table.write_card_16(is_double=False)
        table.write_card_16(is_double=True)

        x = [1., 2., 3.]
        y = [10., 20., 30.]
        x1 = 2.0
        model.add_tablem2(101, x1, x, y, extrap=0, comment='')
        x1 = 3.0
        model.add_tablem2(102, x1, x, y, extrap=1, comment='')
        save_load_deck(model, run_convert=False, run_quality=False, run_remove_unused=False)

        save_load_deck(model)

        #interp = table.interpolate(0.)
        #print('interp =', interp, type(interp))
        #assert np.allclose(interp, [5.5]), interp

    def test_tablem3(self):
        model = BDF(debug=False)
        lines = [
            'TABLEM3, 62, 126.9, 30.0',
            ',2.9, 2.9, 3.6, 4.7, 5.2, 5.7, ENDT',
        ]
        card = model._process_card(lines)
        card = BDFCard(card)
        table = TABLEM3.add_card(card, comment='table')
        table.raw_fields()
        str(table)

        table.write_card(size=8)
        table.write_card_16(is_double=False)
        table.write_card_16(is_double=True)

        x = [1., 2., 3.]
        y = [10., 20., 30.]
        x1 = 2.0
        x2 = 3.0
        model.add_tabled3(101, x1, x2, x, y, extrap=0, comment='')
        x1 = 12.0
        x2 = 13.0
        model.add_tabled3(102, x1, x2, x, y, extrap=1, comment='')
        save_load_deck(model, run_convert=False, run_quality=False, run_remove_unused=False)

        save_load_deck(model)

        #interp = table.interpolate(0.)
        #print('interp =', interp, type(interp))
        #assert np.allclose(interp, [5.5]), interp

    def test_tablem4(self):
        model = BDF(debug=False)
        lines = [
            'TABLEM4, 28, 0.0, 1.0, 0.0, 100.',
            ',2.91, -0.0329, 6.51-5, 0.0, -3.4-7, ENDT',
        ]
        card = model._process_card(lines)
        card = BDFCard(card)
        table = TABLEM4.add_card(card, comment='table')
        table.raw_fields()
        str(table)
        table.raw_fields()
        table.write_card(size=8)
        table.write_card_16(is_double=False)
        table.write_card_16(is_double=True)

        x = [1., 2., 3.]
        y = [10., 20., 30.]

        x1 = 2.0
        x2 = 3.0
        x3 = 4.0
        x4 = 5.0
        a = [4., 5., 6.]
        model.add_tablem4(101, x1, x2, x3, x4, a, comment='tabled4')

        x1 = 12.0
        x2 = 13.0
        x3 = 14.0
        x4 = 15.0
        a = [40., 50., 60.]
        model.add_tablem4(102, x1, x2, x3, x4, a, comment='')

        save_load_deck(model)

        #interp = table.interpolate(0.)
        #print('interp =', interp, type(interp))
        #assert np.allclose(interp, [5.5]), interp

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
