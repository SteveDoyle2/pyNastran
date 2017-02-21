from __future__ import print_function
import unittest
from six import StringIO
from pyNastran.bdf.bdf import read_bdf, BDF, CHBDYG, CaseControlDeck

class TestThermal(unittest.TestCase):
    def test_thermal_1(self):
        """tests various thermal cards"""
        model = BDF(debug=False)
        model.sol = 101
        lines = []
        model.case_control_deck = CaseControlDeck(lines, log=None)

        model.add_grid(11, xyz=[0., 0., 0.])
        model.add_grid(12, xyz=[1., 0., 0.])
        model.add_grid(13, xyz=[1., 1., 0.])
        model.add_grid(14, xyz=[0., 1., 0.])
        model.add_grid(15, xyz=[0., 2., 0.])

        eid = 1
        pid = 1
        mid = 1
        nodes = [11, 12, 13, 14]
        model.add_cquad4(eid, pid, nodes, theta_mcid=0.0, zoffset=0.,
                         TFlag=0,  T1=1.0, T2=1.0, T3=1.0, T4=1.0, comment='')
        model.add_pshell(pid, mid1=1, t=0.1)

        nids = [11, 12, 13, 15]
        pid = 2
        model.add_ctetra(eid, pid, nids)
        model.add_psolid(pid, mid)

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        eid = 2
        Type = 'AREA3'
        chbdyg = CHBDYG(eid, Type, nodes,
                        iview_front=0, ivew_back=0,
                        rad_mid_front=0, rad_mid_back=0,
                        comment='chbdyg')
        with self.assertRaises(ValueError):
            chbdyg.validate()

        Type = 'AREA4'
        chbdyg = model.add_chbdyg(eid, Type, nodes,
                         iview_front=0, ivew_back=0,
                         rad_mid_front=0, rad_mid_back=0,
                         comment='chbdyg')

        eid = 3
        eid2 = 4
        side = 1
        model.add_chbdye(eid, eid2, side,
                         iview_front=0, ivew_back=0,
                         rad_mid_front=0, rad_mid_back=0,
                         comment='chbdye')

        eid = 4
        g1 = 11
        g2 = 12
        pid = 10
        chbdyp = model.add_chbdyp(
            eid, pid, Type, g1, g2,
            g0=0, gmid=None, ce=0,
            iview_front=0, ivew_back=0,
            rad_mid_front=0, rad_mid_back=0,
            e1=None, e2=None, e3=None,
            comment='chbdyp')

        phbdy = model.add_phbdy(pid, af=None, d1=None, d2=None,
                                comment='phbdy')
        model.validate()
        model._verify_bdf(xref=False)
        model.cross_reference()
        model._verify_bdf(xref=True)

        bdf_filename = StringIO()
        bdf_filename2 = StringIO()
        bdf_filename3 = StringIO()

        model.write_bdf(bdf_filename, encoding=None, size=8,
                        is_double=False,
                        interspersed=False,
                        enddata=None, close=False)
        model.write_bdf(bdf_filename2, encoding=None, size=16,
                        is_double=False,
                        interspersed=False,
                        enddata=None, close=False)
        model.write_bdf(bdf_filename3, encoding=None, size=16,
                        is_double=True,
                        interspersed=False,
                        enddata=None, close=False)
        #print(bdf_filename.getvalue())

        bdf_filename2.seek(0)
        model2 = read_bdf(bdf_filename2)



if __name__ == '__main__':
    unittest.main()
