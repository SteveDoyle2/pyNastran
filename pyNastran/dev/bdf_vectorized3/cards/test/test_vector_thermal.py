import unittest
from io import StringIO

from pyNastran.dev.bdf_vectorized3.bdf import read_bdf, BDF, CaseControlDeck, Subcase # CHBDYG,
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck
#from pyNastran.dev.bdf_vectorized3.mesh_utils.mirror_mesh import write_bdf_symmetric
from cpylog import SimpleLogger


class TestThermal(unittest.TestCase):
    def test_bdyor(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log, debug=False)
        surface_type = 'POINT'
        iview_front = None
        iview_back = None
        rad_mid_front = None
        rad_mid_back = None
        pid = 1
        g0 = None
        ce = 42
        card_lines = ['BDYOR', surface_type,
                      iview_front, iview_back,
                      rad_mid_front, rad_mid_back, None,
                      pid, g0,
                      ce]
        #model.add_card(card_lines, 'BDYOR')
        model.add_card_fields(card_lines, 'BDYOR', comment='', has_none=True)
        model.setup()

    def test_radcav(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log, debug=False)
        radcav = model.radcav
        card_fields = [
            'radcav', '1', None, 'NO', None, '0', None, None, '17', '3', '6', '5', '6', '5',
            '5', '6', '5', '6', '1', '6', '3', '1', '3', '3', '1', '1', '5', '3', '5', '1',
            '6', '3', '6', '3', '5', '3', '6']
        model.add_card(card_fields, 'RADCAV')
        model.setup()

    def test_radlst(self):
        log = SimpleLogger(level='warning')
        model = BDF(log=log, debug=False)
        radcav = model.radcav
        card_fields = ['radlst', '10', '1', '10', '20']
        model.add_card(card_fields, 'RADLST')
        model.setup()

    def test_thermal_1(self):
        """tests various thermal cards"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log, debug=False)
        model.sol = 101

        radm = model.radm
        radbc = model.radbc
        qvol = model.qvol
        qvect = model.qvect
        qhbdy = model.qhbdy
        qbdy1 = model.qbdy1
        qbdy2 = model.qbdy2
        qbdy3 = model.qbdy3
        view = model.view

        lines = [
            'SUBCASE 1',
            '  DISP(PLOT) = ALL',
            '  ANALYSIS = HEAT',
            'BEGIN BULK',
        ]
        model.case_control_deck = CaseControlDeck(lines, log=None)

        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 1., 0.])
        model.add_grid(14, [0., 1., 0.])
        model.add_grid(15, [0., 2., 0.])
        model.add_grid(16, [0., 3., 0.])
        model.add_grid(33, [0., 4., 0.])
        model.add_grid(57, [0., 5., 0.])
        model.add_grid(1000, [0., 6., 0.])
        model.add_grid(1001, [0., 7., 0.])

        eid = 1
        pid = 1
        mid = 1
        nodes = [11, 12, 13, 14]
        model.add_cquad4(eid, pid, nodes, theta_mcid=0.0, zoffset=0.,
                         tflag=0, T1=1.0, T2=1.0, T3=1.0, T4=1.0, comment='')
        model.add_pshell(pid, mid1=1, t=0.1)
        k = 1000.
        model.add_mat4(mid, k, cp=0.0, rho=1.0, H=None, mu=None, hgen=1.0,
                       ref_enthalpy=None, tch=None, tdelta=None, qlat=None, comment='')

        eid = 10
        nids = [11, 12, 13, 15]
        pid = 2
        model.add_ctetra(eid, pid, nids)
        model.add_psolid(pid, mid)

        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        RUN_THRMAL = False
        if RUN_THRMAL:
            eid = 2
            Type = 'AREA3'
            chbdygi = CHBDYG(eid, Type, nodes,
                             iview_front=0, iview_back=0,
                            rad_mid_front=0, rad_mid_back=0,
                            comment='chbdyg')
            with self.assertRaises(ValueError):
                chbdyg.validate()

        iview = 10
        Type = 'AREA4'
        chbdyg = model.add_chbdyg(eid, Type, nodes,
                                  iview_front=iview, iview_back=0,
                                  rad_mid_front=0, rad_mid_back=0,
                                  comment='chbdyg')
        #chbdyg.raw_fields()

        eid = 3
        eid2 = 4
        side = 1
        chbdye = model.add_chbdye(eid, eid2, side,
                                  iview_front=0, iview_back=0,
                                  rad_mid_front=0, rad_mid_back=0,
                                  comment='chbdye')
        #chbdye.raw_fields()

        eid = 4
        g1 = 11
        g2 = 12
        pid = 10
        # fails on AREA4 because op2 doesn't support it
        Type = 'LINE'
        chbdyp = model.add_chbdyp(
            eid, pid, Type, g1, g2,
            g0=0, gmid=None, ce=0,
            iview_front=0, iview_back=0,
            rad_mid_front=0, rad_mid_back=0,
            e1=None, e2=None, e3=None,
            comment='chbdyp')
        #chbdyp.raw_fields()

        phbdy = model.add_phbdy(pid, area_factor=None, d1=None, d2=None,
                                comment='phbdy')
        #phbdy.raw_fields()

        #---------------------------
        ta = 2
        ta1 = 2
        pconid = 11
        conv = model.add_conv(eid, pconid, ta, film_node=0, cntrlnd=0,
                              comment='conv')
        #conv.raw_fields()
        pconv = model.add_pconv(
            pconid, mid, form=0,
            exponent_free_convection=0.0,
            free_convection_type=0,
            table_id=None, chlen=None, gidin=None,
            ce=0, e1=None, e2=None, e3=None,
            comment='pconv')
        #pconv.raw_fields()

        pconid = 12
        convm = model.add_convm(eid, pconid, ta1, film_node=0, cntmdot=0,
                                ta2=None, mdot=1.0, comment='convm')
        #convm.raw_fields()
        coef = 0.023
        pconvm = model.add_pconvm(pconid, mid, coef, form=0, flag=0, expr=0.0,
                                  exppi=0.0, exppo=0.0, comment='pconvm')
        #pconvm.raw_fields()

        radmid = 42
        absorb = 0.2
        emissivity = 0.8
        radmi = model.add_radm(radmid, absorb, emissivity, comment='radm')
        #radm.raw_fields()

        famb = 100.
        nodamb = 33
        eids = [1]
        cntrlnd = 1000
        radbci = model.add_radbc(nodamb, famb, cntrlnd, eids, comment='radbc')
        #radbc.raw_fields()

        sid = 43
        qvol_val = 17.
        control_point = 1001
        elements = [1, 2]
        qvoli = model.add_qvol(sid, qvol_val, control_point, elements, comment='qvol')
        #qvol.raw_fields()

        q0 = 18.
        t_source = 19.
        eids = [2]
        qvecti = model.add_qvect(sid, q0, eids, t_source, ce=0,
                                 vector_tableds=None, control_id=0,
                                 comment='qvect')
        #qvecti.raw_fields()

        q0 = 15.8
        flag = 'POINT'
        grids = [-1]
        qhbdyi = model.add_qhbdy(sid, flag, q0, grids, area_factor=None, comment='qhbdy')
        #qhbdy.raw_fields()

        qflux = 20.
        eids = [1]
        qbdy1i = model.add_qbdy1(sid, qflux, eids, comment='qbdy1')
        #qbdy1i.raw_fields()

        eid = 1
        qfluxs = 12.
        qbdy2i = model.add_qbdy2(sid, eid, qfluxs, comment='qbdhy2')
        #qbdy2i.raw_fields()

        q0 = 14.
        cntrlnd = 57
        eids = [1, 2]
        qbdy3i = model.add_qbdy3(sid, q0, cntrlnd, eids, comment='qbdy3')
        #qbdy3i.raw_fields()

        icavity = 12
        model.add_view(iview, icavity, shade='BOTH', nbeta=1, ngamma=1, dislin=0.0, comment='')

        temperature = 13.3
        model.add_tempd(sid, temperature, comment='tempd')

        fields = ['TEMPD', 101, 1., 102, 2., 103, 3., 104, 4.]
        model.add_card(fields, 'TEMPD')

        temperatures = {
            15 : 37.,
            16 : 38.,
        }
        model.add_temp(sid, temperatures)
        model.setup()

        radm.write()
        radbc.write()
        qvol.write()
        qvect.write()
        qhbdy.write()
        qbdy1.write()
        qbdy2.write()
        qbdy3.write()
        view.write()
        save_load_deck(model)
        #-------------------------
        bdf_filename = StringIO()
        bdf_filename2 = StringIO()
        bdf_filename3 = StringIO()
        bdf_filename4 = StringIO()

        model.validate()
        model._verify_bdf(xref=False)
        model.write_bdf(bdf_filename, encoding=None, size=8,
                        is_double=False,
                        interspersed=False,
                        enddata=None, close=False)

        model.cross_reference()
        #model.pop_xref_errors()

        model._verify_bdf(xref=True)
        model.write_bdf(bdf_filename2, encoding=None, size=16,
                        is_double=False,
                        interspersed=False,
                        enddata=None, close=False)
        model.write_bdf(bdf_filename3, encoding=None, size=16,
                        is_double=True,
                        interspersed=False,
                        enddata=None, close=False)
        if RUN_THRMAL:
            write_bdf_symmetric(model, bdf_filename4, encoding=None, size=8,
                                is_double=False,
                                enddata=None, close=False, plane='xz')
        #model.cross_reference()

        #print(bdf_filename.getvalue())

        bdf_filename2.seek(0)
        model2 = read_bdf(bdf_filename2, xref=False, log=log, debug=False)
        #model2.safe_cross_reference()
        model2.setup()
        save_load_deck(model, punch=False, run_renumber=False, run_test_bdf=False)

    def test_thermal_2(self):
        """tests TABLEHT, TABLEH1"""
        model = BDF(debug=False, log=None, mode='msc')
        model.sol = 159
        model.case_control_deck = CaseControlDeck([], log=model.log)
        subcase: Subcase = model.case_control_deck.create_new_subcase(100)
        subcase.add('DISP', 'ALL', ['PLOT'], 'STRESS-type')
        subcase.add('TSTEP', 10, [], 'STRESS-type')
        sid = 10
        N = 40
        DT = 0.1
        NO = 40
        model.add_tstep(sid, N, DT, NO)

        tid = 101
        x = [1., 2., 3.]
        y = [10., 20., 30.]
        model.add_tableh1(tid, x, y, comment='tableh1')

        tid = 101
        model.add_tableh1(tid, x, y, comment='tableh1')

        tid = 110
        tableh1 = model.add_tableh1(tid, x, y, comment='tableh1')
        tableh1.raw_fields()

        tid_tableht = 85
        x = [10., 25., 40.]
        y = [101, 102, 110]

        #This table is referenced only by PCONV entries that define
        #free convection boundary condition properties.
        tableht = model.add_tableht(tid_tableht, x, y, comment='tableht')
        tableht.raw_fields()

        pconv_id = 100
        mid = None
        pconv = model.add_pconv(pconv_id, mid, form=0, exponent_free_convection=0.0,
                                free_convection_type=0,
                                table_id=tid_tableht, chlen=None, gidin=None,
                                ce=0, e1=None, e2=None, e3=None, comment='pconv')
        #pconv.raw_fields()

        # Every surface to which free convection is to be applied must
        # reference a PCONV entry. PCONV is referenced on the CONV Bulk Data entry.
        eid = 1
        ta = 1
        conv = model.add_conv(eid, pconv_id, ta, film_node=0, cntrlnd=0, comment='conv')
        #conv.raw_fields()

        conv = model.add_conv(2, pconv_id, ta, film_node=0, cntrlnd=0, comment='conv')
        conv = model.add_conv(3, pconv_id, ta, film_node=0, cntrlnd=0, comment='conv')

        # CHBDYG, CHBDYE, or CHBDYP surface element identification number.
        eid_fem = 1
        eid_conv = 1
        side = 3 # TODO: 1-6
        chbdye = model.add_chbdye(eid_fem, eid_conv, side, iview_front=0, iview_back=0,
                                  rad_mid_front=0, rad_mid_back=0, comment='chbdye')

        eid_fem = 2
        nodes = [1, 2, 3]
        surface_type = 'AREA3'
        chbdyg = model.add_chbdyg(eid_fem, surface_type, nodes,
                                  iview_front=0, iview_back=0,
                                  rad_mid_front=0, rad_mid_back=0, comment='chbdyg')

        eid_fem = 3
        pid_phybdy = 10
        g1 = 1
        g2 = None
        surface_type = 'POINT'
        chbdyp = model.add_chbdyp(eid_fem, pid_phybdy, surface_type, g1, g2,
                                  g0=0, gmid=None, ce=0,
                                  iview_front=0, iview_back=0,
                                  rad_mid_front=0, rad_mid_back=0,
                                  e1=None, e2=None, e3=None, comment='chbdyp')
        phbdy = model.add_phbdy(pid_phybdy, area_factor=None, d1=None, d2=None, comment='phbdy')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0., 1., 0.])

        #if RUN_THRMAL:
            #chbdye.raw_fields()
            #chbdyg.raw_fields()
            #chbdyp.raw_fields()

        model.validate()
        save_load_deck(model)

if __name__ == '__main__':   # pragma: no cover
    unittest.main()
