"""tests b-list elements"""
import unittest
import numpy as np

from cpylog import get_logger
from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.bdf import CGAP, PGAP, CBUSH, CFAST
from pyNastran.bdf.cards.test.utils import save_load_deck, read_write_op2_geom
from pyNastran.bdf.mesh_utils.mass_properties import mass_properties


class TestPlotElements(unittest.TestCase):

    def test_plotel_01(self):
        """tests a PLOTEL"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.set_error_storage(
            nparse_errors=0, stop_on_parsing_error=True,
            nxref_errors=0, stop_on_xref_error=True)
        eid = 9
        nodes = [10, 11]
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        plotel = model.add_plotel(eid, nodes, comment='plotel')
        plotel.raw_fields()
        plotel.write_card(size=8, is_double=False)
        plotel.write_card(size=16, is_double=False)
        plotel.write_card(size=16, is_double=True)
        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        save_load_deck(model, xref='standard', punch=True)

    def test_plotel_02(self):
        """tests a PLOTEL3/PLOTEL4"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.set_error_storage(
            nparse_errors=0, stop_on_parsing_error=True,
            nxref_errors=0, stop_on_xref_error=True)
        eid = 9
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        model.add_grid(12, [1., 1., 0.])
        model.add_grid(13, [0., 1., 0.])
        plotel3 = model.add_plotel3(eid, [10, 11, 12], comment='plotel')
        plotel4 = model.add_plotel4(eid+1, [10, 11, 12, 13], comment='plotel')
        plotel3.raw_fields()
        plotel3.write_card(size=8, is_double=False)
        plotel3.write_card(size=16, is_double=False)
        plotel3.write_card(size=16, is_double=True)

        plotel4.raw_fields()
        plotel4.write_card(size=8, is_double=False)
        plotel4.write_card(size=16, is_double=False)
        plotel4.write_card(size=16, is_double=True)
        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        save_load_deck(
            model, xref='standard', punch=True,
            run_op2_writer=False)

    def test_plotel_03(self):
        """tests a PLOTEL6/PLOTEL8"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.set_error_storage(
            nparse_errors=0, stop_on_parsing_error=True,
            nxref_errors=0, stop_on_xref_error=True)
        eid = 9
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        model.add_grid(12, [1., 1., 0.])
        model.add_grid(13, [0., 1., 0.])

        model.add_grid(20, [0.5, 0., 0.])
        model.add_grid(21, [1., 0.5, 0.])
        model.add_grid(22, [0.5, 1., 0.])
        model.add_grid(23, [0., 0.5, 0.])
        plotel6 = model.add_plotel6(eid, [10, 11, 12, 20,21, 22], comment='plotel')
        plotel8 = model.add_plotel8(eid+1, [10, 11, 12, 13,
                                            20, 21, 22, 23], comment='plotel')
        plotel6.raw_fields()
        plotel6.write_card(size=8, is_double=False)
        plotel6.write_card(size=16, is_double=False)
        plotel6.write_card(size=16, is_double=True)

        plotel8.raw_fields()
        plotel8.write_card(size=8, is_double=False)
        plotel8.write_card(size=16, is_double=False)
        plotel8.write_card(size=16, is_double=True)
        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        save_load_deck(
            model, xref='standard', punch=True,
            run_op2_writer=False)

    def test_plottet(self):
        """tests a PLOTTET"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.set_error_storage(
            nparse_errors=0, stop_on_parsing_error=True,
            nxref_errors=0, stop_on_xref_error=True)
        eid = 9
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        model.add_grid(12, [1., 1., 0.])
        model.add_grid(13, [0.5, 0.5, 1.])
        plottet = model.add_plottet(eid, [10, 11, 12, 13], comment='plotel')
        plottet.raw_fields()
        plottet.write_card(size=8, is_double=False)
        plottet.write_card(size=16, is_double=False)
        plottet.write_card(size=16, is_double=True)
        elem = plottet
        assert len(model._type_to_id_map[elem.type]) == 1, model._type_to_id_map

        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        save_load_deck(
            model, xref='standard', punch=True,
            run_op2_writer=False)

    def test_plotpyr(self):
        """tests a PLOTPYR"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.set_error_storage(
            nparse_errors=0, stop_on_parsing_error=True,
            nxref_errors=0, stop_on_xref_error=True)
        eid = 9
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(11, [1., 0., 0.])
        model.add_grid(12, [1., 1., 0.])
        model.add_grid(13, [1., 0., 0.])
        model.add_grid(20, [0., 0., 1.])
        plotpyr = model.add_plotpyr(eid, [10, 11, 12, 13, 20], comment='plotel')
        plotpyr.raw_fields()
        plotpyr.write_card(size=8, is_double=False)
        plotpyr.write_card(size=16, is_double=False)
        plotpyr.write_card(size=16, is_double=True)

        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        save_load_deck(
            model, xref='standard', punch=True,
            run_op2_writer=False)

    def test_plotpen(self):
        """tests a PLOTPEN"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.set_error_storage(
            nparse_errors=0, stop_on_parsing_error=True,
            nxref_errors=0, stop_on_xref_error=True)
        eid = 9
        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 1., 0.])

        model.add_grid(21, [0., 0., 1.])
        model.add_grid(22, [1., 0., 1.])
        model.add_grid(23, [1., 1., 1.])
        plotpen = model.add_plotpen(eid, [11, 12, 13, 21, 22, 23], comment='plotel')
        plotpen.raw_fields()
        plotpen.write_card(size=8, is_double=False)
        plotpen.write_card(size=16, is_double=False)
        plotpen.write_card(size=16, is_double=True)

        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        save_load_deck(
            model, xref='standard', punch=True,
            run_op2_writer=False)

    def test_plothex(self):
        """tests a PLOTHEX"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.set_error_storage(
            nparse_errors=0, stop_on_parsing_error=True,
            nxref_errors=0, stop_on_xref_error=True)
        eid = 9
        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 1., 0.])
        model.add_grid(14, [0., 1., 0.])

        model.add_grid(21, [0., 0., 1.])
        model.add_grid(22, [1., 0., 1.])
        model.add_grid(23, [1., 1., 1.])
        model.add_grid(24, [0., 1., 1.])
        plothex = model.add_plothex(eid, [11, 12, 13, 14, 21, 22, 23, 24], comment='plotel')
        plothex.raw_fields()
        plothex.write_card(size=8, is_double=False)
        plothex.write_card(size=16, is_double=False)
        plothex.write_card(size=16, is_double=True)

        model.cross_reference()
        model.uncross_reference()
        model.safe_cross_reference()
        save_load_deck(
            model, xref='standard', punch=True,
            run_op2_writer=False)

class TestElements(unittest.TestCase):
    def test_cbush_01(self):
        """tests a CBUSH"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        lines = ['cbush,101,102,1,,,,,0']
        card = model._process_card(lines)
        card = BDFCard(card)

        size = 8
        elem = CBUSH.add_card(card)
        self.assertEqual(elem.eid, 101)
        self.assertEqual(elem.Pid(), 102)
        elem.write_card(size, 'dummy')
        elem.raw_fields()

        pid = 101
        k_tables = [201]
        b_tables = [202]
        ge_tables = [203]
        kn_tables = [204]
        model.add_pbusht(pid, k_tables, b_tables, ge_tables, kn_tables, comment='pbusht')
        save_load_deck(model)

    def test_cbush_02(self):
        """tests a CGAP/PGAP"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid = 10
        pid = 11
        nodes = [1, None]
        g0 = None
        x = [1., 2., 3.]
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
        # assert cbush.ga == 1
        # assert cbush.gb is None
        model.add_cord2r(1, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.])
        cbush1 = model.add_cbush(
            eid, pid, nodes, x, g0,
            cid=1, s=0.5, ocid=1,
            si=[1., 2., 3.])
        cbush2 = model.add_cbush(eid+1, pid, [1, 2], x, g0)
        model.add_pbush(pid, k=[1.,], b=[2.])
        model.cross_reference()
        cbush1.Centroid()
        cbush2.Centroid()
        save_load_deck(model)

    def test_cbush_optistruct(self):
        """tests a CGAP/PGAP in optistruct"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid = 10
        pid = 11
        nodes = [1, None]
        g0 = None
        x = [1., 2., 3.]
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])
        # assert cbush.ga == 1
        # assert cbush.gb is None
        model.add_cord2r(1, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.])
        cbush1 = model.add_cbush(
            eid, pid, nodes, x, g0,
            cid=1, s=0.5, ocid=1,
            si=[1., 2., 3.])
        cbush2 = model.add_cbush(eid+1, pid, [1, 2], x, g0)
        model.add_pbush_optistruct(pid, k=[1.,], b=[2.])
        model.cross_reference()
        cbush1.Centroid()
        cbush2.Centroid()
        save_load_deck(model, run_op2_writer=False)

    def test_cgap_01(self):
        """tests a CGAP/PGAP"""
        log = get_logger(level='warning')
        model = BDF(log=log)
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
        model.elements[899] = cgap

        lines = ['PGAP    90                      1.E+5']
        card = model._process_card(lines)
        card = BDFCard(card)

        pgap = PGAP.add_card(card, comment='pgap')
        pgap.write_card(size=8)
        pgap.write_card(size=16)
        self.assertEqual(pgap.Pid(), 90)
        pgap.raw_fields()
        model.properties[90] = pgap

        model.add_grid(3, [-1., 0., 0.])
        model.add_grid(21, [0., 0., 0.])
        model.add_grid(99, [1., 0., 0.])
        eid = 100
        pid = 90
        nids = [21, 99]
        x = None
        g0 = 3
        cid = None
        cgap = model.add_cgap(eid, pid, nids,
                              x, g0, cid, comment='cgap')
        node_ids = cgap.node_ids
        assert node_ids == [21, 99], node_ids
        self.assertEqual(cgap.eid, 100)
        self.assertEqual(cgap.Pid(), 90)
        cgap.write_card(size=8)
        cgap.raw_fields()

        model.cross_reference()
        save_load_deck(model)

    def test_cweld_elemid_x(self):
        connectype = 'ELEMID'
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid = 10
        pid = 20
        gs = None
        ga = 32
        gb = 33
        mcid = 0
        shida = 101
        shidb = 102
        x = [1., 2., 3.]
        model.add_cweld(eid, pid, gs, connectype, ga, gb,
                        mcid=mcid, shida=shida, shidb=shidb,
                        x=x, comment='weld')
        mid = 201
        d = 0.1
        model.add_pweld(pid, mid, d, comment='pweld')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.add_grid(31, [0., 0., 0.])
        model.add_grid(32, [1., 0., 0.])
        model.add_grid(33, [2., 0., 0.])
        save_load_deck(model)

    def test_cweld_elpat(self):
        connectype = 'ELPAT'
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid = 10
        pid = 20
        gs = None
        ga = 32
        gb = 33
        mcid = 0
        shida = 101
        shidb = 102
        x = [1., 2., 3.]
        model.add_cweld(eid, pid, gs, connectype, ga, gb,
                        mcid=mcid, shida=shida, shidb=shidb,
                        x=x, comment='weld')
        mid = 201
        d = 0.1
        model.add_pweld(pid, mid, d, comment='pweld')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.add_grid(31, [0., 0., 0.])
        model.add_grid(32, [1., 0., 0.])
        model.add_grid(33, [2., 0., 0.])
        save_load_deck(model)

    def test_cweld_partpat(self):
        connectype = 'ELEMID'
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid = 10
        pid = 20
        gs = None
        ga = 32
        gb = 33
        mcid = 0
        shida = 101
        shidb = 102
        model.add_cweld(eid, pid, gs, connectype, ga, gb,
                        mcid=mcid, shida=shida, shidb=shidb,
                        x=None, comment='weld')
        mid = 201
        d = 0.1
        model.add_pweld(pid, mid, d, comment='pweld')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.add_grid(31, [0., 0., 0.])
        model.add_grid(32, [1., 0., 0.])
        model.add_grid(33, [2., 0., 0.])
        save_load_deck(model)

    def test_cweld_elemid_x_none(self):
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid = 10
        pid = 20
        gs = None
        ga = 32
        gb = 33
        mcid = 0
        shida = 101
        shidb = 102
        connectype = 'ELEMID'
        model.add_cweld(eid, pid, gs, connectype, ga, gb,
                        mcid, shida, shidb,
                        x=None, comment='weld')
        mid = 201
        d = 0.1
        model.add_pweld(pid, mid, d, comment='pweld')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        model.add_grid(31, [0., 0., 0.])
        model.add_grid(32, [1., 0., 0.])
        model.add_grid(33, [2., 0., 0.])
        save_load_deck(model)

    def test_cfast(self):
        """tests a CFAST/PFAST"""
        log = get_logger(level='warning')
        model = BDF(log=log)

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
        unused_mat1 = model.add_mat1(mid, E, G, nu, rho=0.2)

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
                                comment='pfast')
        pfastb = model.add_pfast(
            pid+1, d, kt1, kt2, kt3, mcid=0, mflag=0, kr1=0.,
            kr2=0., kr3=0., mass=0., ge=0.)
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
        mass_properties(model)
        #save_load_deck(model, run_convert=False)
        read_write_op2_geom(
            model, run_op2_writer=True, run_op2_reader=True,
            nastran_format='msc')
        save_load_deck(model)

    def test_cbush1d(self):
        model = BDF(debug=False)

        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        nids = [2, 3]
        eid = 10
        pid = 100
        model.add_cbush1d(eid, pid, nids, cid=None, comment='cbush1d')
        model.add_pbush1d(pid, k=0., c=0., m=0., sa=0., se=0., optional_vars=None,
                          comment='pbush1d')

        model.pop_parse_errors()
        model.cross_reference()
        save_load_deck(model, run_op2_reader=False)

    def test_cbush2d(self):
        log = get_logger(level='warning')
        model = BDF(log=log)

        model.add_grid(2, [0., 0., 0.])
        model.add_grid(3, [1., 0., 0.])
        nids = [2, 3]
        eid = 10
        pid = 100
        unused_cbush2d = model.add_cbush2d(eid, pid, nids, cid=0, plane='XY', sptid=None,
                                           comment='cbush2d')
        #model.add_pbush2d()

        #model.pop_parse_errors()
        #model.cross_reference()
        save_load_deck(model, run_convert=False, xref=False, run_renumber=False, run_test_bdf=False)

    def test_crac2d(self):
        log = get_logger(level='warning')
        model = BDF(log=log)

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
        nids = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        eid = 10
        pid = 100
        mid = 1000
        thick = 20.
        iplane = 1
        crac2d = model.add_crac2d(eid, pid, nids, comment='crac2d')
        prac2d = model.add_prac2d(pid, mid, thick, iplane, nsm=0., gamma=0.5, phi=180., comment='')
        crac2d.raw_fields()
        prac2d.raw_fields()

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        #model.add_pbush2d()

        #model.pop_parse_errors()
        #model.cross_reference()
        save_load_deck(model, run_convert=False)

    def test_crac3d(self):
        log = get_logger(level='warning')
        model = BDF(log=log)

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
        nids = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        eid = 10
        pid = 100
        mid = 1000
        unused_thick = 20.
        unused_iplane = 1
        crac3d = model.add_crac3d(eid, pid, nids, comment='crac3d')
        prac3d = model.add_prac3d(pid, mid, gamma=0.5, phi=180., comment='crac3d')
        crac3d.raw_fields()
        prac3d.raw_fields()

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)
        #model.add_pbush2d()

        #model.pop_parse_errors()
        #model.cross_reference()
        save_load_deck(model, run_convert=False)

    def test_genel_1(self):
        """tests a GENEL element"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        eid = 1

        model.add_grid(1, [0., 0., 0.])
        model.add_grid(13, [0., 0., 0.])
        model.add_grid(42, [0., 0., 0.])
        model.add_grid(24, [0., 0., 0.])
        model.add_grid(6, [0., 0., 0.])
        model.add_grid(33, [0., 0., 0.])
        ul = np.array([
            [1, 1],
            [13, 4],
            [42, 0],
            [24, 0],
        ], dtype='int32')
        ud = np.array([
            [6, 2],
            [33, 0],
        ], dtype='int32')

        #+-------+------+-----+------+------+------+------+-------+------+
        #| GENEL |  629 |     |  1   |  1   |  13  |  4   |   42  |   0  |
        #|       |  24  |  2  |      |      |      |      |       |      |
        #|       |  UD  |     |  6   |  2   |  33  |  0   |       |      |
        #|       |  Z   | 1.0 | 2.0  | 3.0  | 4.0  | 5.0  |  6.0  | 7.0  |
        #|       |  8.0 | 9.0 | 10.0 |      |      |      |       |      |
        #|       |  S   | 1.5 | 2.5  | 3.5  | 4.5  | 5.5  |  6.5  | 7.5  |
        #|       |  8.5 |     |      |      |      |      |       |      |
        #+-------+------+-----+------+------+------+------+-------+------+
        z = np.array([1., 2., 3., 4., 5., 6., 7., 8., 9., 10.])
        s = np.array([1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
        k = z
        model.add_genel_flexibility(eid, ul, ud, z, s)

        elem = model.add_genel_stiffness(10, ul, ud, z, k)
        str(elem)
        elem.eid = 11
        fields = elem.raw_fields()

        elem = model.add_card(fields, 'GENEL', comment='card', is_list=True, has_none=True)
        elemi = model.elements[11]
        str(elemi)
        #print('\n'+str(elem))

        elem = model.add_genel_stiffness(20, ul, ud, k)
        elem.eid = 21
        str(elem)
        fields = elem.raw_fields()
        elem = model.add_card(fields, 'GENEL', comment='card', is_list=True, has_none=True)
        elemi = model.elements[21]
        str(elemi)
        save_load_deck(model)

    def test_genel_2(self):
        """tests a GENEL element"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 0.])

        card_lines = [
            'genel,12001,,1,1,,,,,',
            ',ud,,2,1,,,,,',
            ',k,1.0+8,,,,,,,',
            ',s,-1.0+8',
        ]
        model.add_card(card_lines, 'GENEL', comment='', is_list=False, has_none=True)
        genel = model.elements[12001]
        genel.raw_fields()
        str(genel)
        assert len(genel.ul.ravel()) == 2, genel.ul
        save_load_deck(model, xref='standard', punch=True)


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
