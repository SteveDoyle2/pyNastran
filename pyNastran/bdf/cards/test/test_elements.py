"""tests b-list elements"""
import unittest

from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.bdf import CGAP, PGAP, CDAMP1, CBUSH, CFAST
from pyNastran.bdf.cards.test.utils import save_load_deck



class TestElements(unittest.TestCase):

    def test_plotel_01(self):
        """tests a PLOTEL"""
        model = BDF(debug=False)
        eid = 9
        nodes = [10, 11]
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        plotel = model.add_plotel(eid, nodes, comment='plotel')
        plotel.write_card(size=8, is_double=False)
        plotel.write_card(size=16, is_double=False)
        plotel.write_card(size=16, is_double=True)
        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        save_load_deck(model, xref='standard', punch=True, run_remove_unused=True,
                       run_convert=True, run_renumber=False, run_mirror=True,
                       run_save_load=True)

    def test_cbush_01(self):
        """tests a CBUSH"""
        model = BDF(debug=False)
        lines = ['cbush,101,102,1,,,,,0']
        card = model._process_card(lines)
        card = BDFCard(card)

        size = 8
        elem = CBUSH.add_card(card)
        self.assertEqual(elem.eid, 101)
        self.assertEqual(elem.Pid(), 102)
        elem.write_card(size, 'dummy')
        elem.raw_fields()

    def test_cdamp1_01(self):
        """tests a CDAMP1"""
        model = BDF(debug=False)
        lines = ['CDAMP1, 2001, 20, 1001, 1']
        card = model._process_card(lines)
        card = BDFCard(card)

        size = 8
        elem = CDAMP1.add_card(card)
        self.assertEqual(elem.eid, 2001)
        self.assertEqual(elem.Pid(), 20)
        node_ids = elem.node_ids
        assert node_ids == [1001, None], node_ids
        elem.write_card(size, 'dummy')
        elem.raw_fields()

    def test_gap_01(self):
        """tests a CGAP/PGAP"""
        model = BDF(debug=False)
        lines = ['CGAP    899     90      21      99      0.      1.      0.      0']
        card = model._process_card(lines)
        card = BDFCard(card)

        cgap = CGAP.add_card(card, comment='cgap')
        node_ids = cgap.node_ids
        assert node_ids == [21, 99], node_ids
        self.assertEqual(cgap.eid, 899)
        self.assertEqual(cgap.Pid(), 90)
        cgap.write_card(size=8)
        cgap.raw_fields()

        lines = ['PGAP    90                      1.E+5']
        card = model._process_card(lines)
        card = BDFCard(card)

        pgap = PGAP.add_card(card, comment='pgap')
        pgap.write_card(size=8)
        pgap.write_card(size=16)
        self.assertEqual(pgap.Pid(), 90)
        pgap.raw_fields()

    def test_cfast(self):
        """tests a CFAST/PFAST"""
        model = BDF(debug=False)

        eid1 = 10
        pid = 11
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(eid1, pid, [1, 2, 3, 4])

        eid2 = 12
        model.add_grid(11, [0., 0., 1.])
        model.add_grid(12, [1., 0., 1.])
        model.add_grid(13, [1., 1., 1.])
        model.add_grid(14, [0., 1., 1.])
        model.add_cquad4(eid2, pid, [11, 12, 13, 14])

        mid = 13
        model.add_pshell(pid, mid1=mid, t=0.1)

        E = 1e7
        nu = 0.3
        G = E / (2. * (1. + nu))
        E = None
        mat1 = model.add_mat1(mid, E, G, nu, rho=0.2)

        eid = 14
        pid = 15
        Type = 'ELEM'
        ida = eid1
        idb = eid2
        cfast = model.add_cfast(eid, pid, Type, ida, idb,
                                gs=1, ga=None, gb=None,
                                xs=None, ys=None, zs=None, comment='cfast')


        Type = 'fake'
        pid2 = None
        cfast2 = CFAST(eid, pid2, Type, ida, idb,
                       gs=None, ga=None, gb=None,
                       xs=None, ys=None, zs=None, comment='')
        with self.assertRaises(TypeError):
            cfast2.validate()

        cfast2.Type = 'ELEM'
        with self.assertRaises(ValueError):
            cfast2.validate()
        cfast2.ga = 3
        cfast2.xs = 4.
        cfast2.ys = 4.
        cfast2.zs = 4.
        cfast2.validate()

        d = 1.0
        kt1 = 1.0
        kt2 = 1.0
        kt3 = 0.1
        pfast = model.add_pfast(pid, d, kt1, kt2, kt3, mcid=-1, mflag=0, kr1=0.,
                                kr2=0., kr3=0., mass=0., ge=0.,
                                comment='')
        model.validate()

        cfast.raw_fields()
        pfast.raw_fields()
        cfast.write_card()
        pfast.write_card()
        model._verify_bdf(xref=False)
        model.cross_reference()
        model.pop_xref_errors()
        model._verify_bdf(xref=True)
        cfast.raw_fields()
        pfast.raw_fields()
        cfast.write_card()
        pfast.write_card()

        model.uncross_reference()
        model.safe_cross_reference()
        model.mass_properties()
        save_load_deck(model, run_convert=False)

    def test_cbush1d(self):
        model = BDF(debug=False)

        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        nids = [2, 3]
        eid = 10
        pid = 100
        model.add_cbush1d(eid, pid, nids, cid=None, comment='cbush1d')
        model.add_pbush1d(pid, k=0., c=0., m=0., sa=0., se=0., optional_vars=None, comment='pbush1d')

        model.pop_parse_errors()
        model.cross_reference()
        save_load_deck(model)

    def test_cbush2d(self):
        model = BDF(debug=False)

        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        nids = [2, 3]
        eid = 10
        pid = 100
        cbush2d = model.add_cbush2d(eid, pid, nids, cid=0, plane='XY', sptid=None, comment='cbush2d')
        #model.add_pbush2d()

        #model.pop_parse_errors()
        #model.cross_reference()
        #save_load_deck(model)

    def test_crac2d(self):
        model = BDF(debug=False)

        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        model.add_grid(4, [1., 0., 0.])
        model.add_grid(5, [1., 0., 0.])
        model.add_grid(6, [1., 0., 0.])
        model.add_grid(7, [1., 0., 0.])
        model.add_grid(8, [1., 0., 0.])
        model.add_grid(9, [1., 0., 0.])
        model.add_grid(10, [1., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 0., 0.])
        model.add_grid(14, [1., 0., 0.])
        model.add_grid(15, [1., 0., 0.])
        model.add_grid(16, [1., 0., 0.])
        model.add_grid(17, [1., 0., 0.])
        model.add_grid(18, [1., 0., 0.])
        model.add_grid(19, [1., 0., 0.])
        nids = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18]
        eid = 10
        pid = 100
        mid = 1000
        thick = 20.
        iplane = 1
        crac2d = model.add_crac2d(eid, pid, nids, comment='crac2d')
        prac2d = model.add_prac2d(pid, mid, thick, iplane, nsm=0., gamma=0.5, phi=180., comment='')

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        #model.add_pbush2d()

        #model.pop_parse_errors()
        #model.cross_reference()
        save_load_deck(model, run_convert=False)

    def test_crac3d(self):
        model = BDF(debug=False)

        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        model.add_grid(4, [1., 0., 0.])
        model.add_grid(5, [1., 0., 0.])
        model.add_grid(6, [1., 0., 0.])
        model.add_grid(7, [1., 0., 0.])
        model.add_grid(8, [1., 0., 0.])
        model.add_grid(9, [1., 0., 0.])
        model.add_grid(10, [1., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 0., 0.])
        model.add_grid(14, [1., 0., 0.])
        model.add_grid(15, [1., 0., 0.])
        model.add_grid(16, [1., 0., 0.])
        model.add_grid(17, [1., 0., 0.])
        model.add_grid(18, [1., 0., 0.])
        model.add_grid(19, [1., 0., 0.])
        nids = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13, 14, 15, 16, 17, 18, 19]
        eid = 10
        pid = 100
        mid = 1000
        thick = 20.
        iplane = 1
        crac3d = model.add_crac3d(eid, pid, nids, comment='crac3d')
        prac3d = model.add_prac3d(pid, mid, gamma=0.5, phi=180., comment='crac3d')

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        #model.add_pbush2d()

        #model.pop_parse_errors()
        #model.cross_reference()
        save_load_deck(model, run_convert=False)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
