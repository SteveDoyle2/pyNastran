"""tests dynamic cards and dynamic load cards"""
import unittest

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.test.utils import save_load_deck


class TestMethods(unittest.TestCase):
    """
    The cards tested are:
     * TSTEP
    """
    def test_eigb(self):
        """tests an EIGB card"""
        model = BDF(debug=None)

        sid = 42
        L1 = 2.
        L2 = 42.
        nep = 1
        ndp = 31
        ndn = 3
        norm = 'cat'
        G = 31
        C = 4256
        method = 'INV'
        eigb = model.add_eigb(sid, method, L1, L2, nep, ndp, ndn, norm, G, C,
                              comment='eigb')
        eigb.raw_fields()

        sid = 43
        L1 = 2.
        L2 = -2.
        #nep = 1
        #ndp = 31
        #ndn = 3
        #norm = 'cat'
        #G = 31
        #C = 4256
        method = 'INV'
        with self.assertRaises(RuntimeError):
            eigb = model.add_eigb(sid, method, L1, L2, nep, ndp, ndn, norm, G, C,
                                  comment='eigb')

        model.pop_parse_errors()
        eigb.raw_fields()
        model.validate()
        save_load_deck(model)

    def test_eigp(self):
        """tests an EIGP card"""
        model = BDF(debug=False)

        sid = 15
        alpha1 = -5.2
        omega1 = 0.0
        m1 = 2

        alpha2 = 6.3
        omega2 = 5.5
        m2 = 3
        eigp = model.add_eigp(sid, alpha1, omega1, m1, alpha2, omega2, m2, comment='eigp')
        eigp.raw_fields()

        unused_eigp2 = model.CMethod(sid)
        model.validate()
        save_load_deck(model)

    def test_eigr(self):
        """tests an EIGR card"""
        model = BDF(debug=False)
        sid = 1
        nd = -42
        eigr = model.add_eigr(sid, method='LAN', f1=None, f2=None, ne=None, nd=nd,
                              norm='MASS', G=None, C=None, comment='eigr')

        sid = 2
        eigr = model.add_eigr(sid, method='SINV', f1=None, f2=None, ne=None, nd=None,
                              norm='MASS', G=None, C=None, comment='')

        sid = 3
        ne = -42
        eigr = model.add_eigr(sid, method='INV', f1=None, f2=None, ne=ne, nd=None,
                              norm='MASS', G=None, C=None, comment='')

        sid = 4
        eigr = model.add_eigr(sid, method='GIV', f1=None, f2=None, ne=None, nd=None,
                              norm='MASS', G=None, C=None, comment='')

        sid = 5
        eigr = model.add_eigr(sid, method='MGIV', f1=None, f2=None, ne=None, nd=None,
                              norm='MASS', G=None, C=None, comment='')

        sid = 6
        eigr = model.add_eigr(sid, method='HOU', f1=None, f2=None, ne=None, nd=None,
                              norm='MASS', G=None, C=None, comment='')

        sid = 7
        eigr = model.add_eigr(sid, method='MHOU', f1=None, f2=None, ne=None, nd=None,
                              norm='MASS', G=None, C=None, comment='')

        sid = 8
        nd = -6
        G = -42
        C = 123456
        eigr = model.add_eigr(sid, method='LAN', f1=None, f2=None, ne=None, nd=nd,
                              norm='POINT', G=G, C=C, comment='')

        sid = 9
        nd = -1
        eigr = model.add_eigr(sid, method='LAN', f1=None, f2=None, ne=None, nd=nd,
                              norm='MAX', G=None, C=None, comment='')

        model.pop_parse_errors()
        eigr.raw_fields()
        model.validate()
        save_load_deck(model)

    def test_eigrl(self):
        """tests an EIGRL card"""
        model = BDF(debug=False)
        sid = 1
        eigrl = model.add_eigrl(sid, v1=None, v2=None, nd=None, msglvl=0,
                                maxset=None, shfscl=None, norm=None,
                                options=None, values=None, comment='')
        eigrl.raw_fields()

        lines = [
            'EIGRL          2      0.    800.       6       4',
            '          NUMS=4',
        ]
        model.add_card_lines(lines, 'EIGRL', comment='eigrl')
        eigrl = model.methods[1]
        str(eigrl)

        model.pop_parse_errors()
        model.validate()
        save_load_deck(model)

    def test_eigc_1(self):
        """tests the EIGC"""
        model = BDF(debug=True, log=None, mode='msc')
        lines = [
            'EIGC    10      CLAN    MAX                     1.E-12',
            '                                                        20',
        ]
        model.add_card_lines(lines, 'EIGC', comment='', has_none=True)
        model.pop_parse_errors()
        eigc = model.cMethods[10]
        str(eigc)

    def test_eigc_2(self):
        """tests the EIGC"""
        model = BDF(debug=True, log=None, mode='msc')
        lines = [
            'EIGC    40      CLAN    MAX                     1.E-12',
            '        0.0     0.0                                     5',
            '        0.0      5.0                                    5',
            '        0.0     10.0                                    5',
            '        5.0      5.0                                    5',
            '        0.0     20.0                                    5',
            '        20.0    10.0                                    5',
            '        10.0    10.0                                    5',
        ]
        eigc = model.add_card_lines(lines, 'EIGC', comment='', has_none=True)
        model.pop_parse_errors()
        eigc = model.cMethods[40]
        str(eigc)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
