import unittest
from io import StringIO

from pyNastran.bdf.bdf import read_bdf, BDF, CHBDYG, CaseControlDeck
from pyNastran.bdf.cards.test.utils import save_load_deck
from pyNastran.bdf.mesh_utils.mirror_mesh import write_bdf_symmetric
from cpylog import SimpleLogger

class TestThermal(unittest.TestCase):
    def test_thermal_1(self):
        """tests various thermal cards"""
        log = SimpleLogger(level='warning')
        model = BDF(log=log, debug=False)
        model.sol = 101
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

        eid = 1
        pid = 1
        mid = 1
        nodes = [11, 12, 13, 14]
        model.add_cquad4(eid, pid, nodes, theta_mcid=0.0, zoffset=0.,
                         tflag=0, T1=1.0, T2=1.0, T3=1.0, T4=1.0, comment='')
        model.add_pshell(pid, mid1=1, t=0.1)

        eid = 10
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
                        iview_front=0, iview_back=0,
                        rad_mid_front=0, rad_mid_back=0,
                        comment='chbdyg')
        with self.assertRaises(ValueError):
            chbdyg.validate()

        Type = 'AREA4'
        chbdyg = model.add_chbdyg(eid, Type, nodes,
                                  iview_front=0, iview_back=0,
                                  rad_mid_front=0, rad_mid_back=0,
                                  comment='chbdyg')
        chbdyg.raw_fields()

        eid = 3
        eid2 = 4
        side = 1
        chbdye = model.add_chbdye(eid, eid2, side,
                                  iview_front=0, iview_back=0,
                                  rad_mid_front=0, rad_mid_back=0,
                                  comment='chbdye')
        chbdye.raw_fields()

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
        chbdyp.raw_fields()

        phbdy = model.add_phbdy(pid, af=None, d1=None, d2=None,
                                comment='phbdy')
        phbdy.raw_fields()

        #---------------------------
        ta = 2
        ta1 = 2
        pconid = 11
        conv = model.add_conv(eid, pconid, ta, film_node=0, cntrlnd=0,
                              comment='conv')
        conv.raw_fields()
        pconv = model.add_pconv(pconid, mid, form=0, expf=0.0, ftype=0,
                                tid=None, chlen=None, gidin=None,
                                ce=0, e1=None, e2=None, e3=None,
                                comment='pconv')
        pconv.raw_fields()

        pconid = 12
        convm = model.add_convm(eid, pconid, ta1, film_node=0, cntmdot=0,
                                ta2=None, mdot=1.0, comment='convm')
        convm.raw_fields()
        coef = 0.023
        pconvm = model.add_pconvm(pconid, mid, coef, form=0, flag=0, expr=0.0,
                                  exppi=0.0, exppo=0.0, comment='pconvm')
        pconvm.raw_fields()

        radmid = 42
        absorb = 0.2
        emissivity = 0.8
        radm = model.add_radm(radmid, absorb, emissivity, comment='radm')
        radm.raw_fields()

        famb = 100.
        nodamb = 33
        eids = [1]
        cntrlnd = 1000
        radbc = model.add_radbc(nodamb, famb, cntrlnd, eids, comment='radbc')
        radbc.raw_fields()

        sid = 43
        qvol = 17.
        control_point = 1001
        elements = [1, 2]
        qvol = model.add_qvol(sid, qvol, control_point, elements, comment='qvol')
        qvol.raw_fields()

        q0 = 18.
        t_source = 19.
        eids = [2]
        qvect = model.add_qvect(sid, q0, eids, t_source, ce=0,
                                vector_tableds=None, control_id=0,
                                comment='qvect')
        qvect.raw_fields()

        q0 = 15.8
        flag = 'POINT'
        grids = [-1]
        qhbdy = model.add_qhbdy(sid, flag, q0, grids, af=None, comment='qhbdy')
        qhbdy.raw_fields()

        qflux = 20.
        eids = [1]
        qbdy1 = model.add_qbdy1(sid, qflux, eids, comment='qbdy1')
        qbdy1.raw_fields()

        eid = 1
        qfluxs = 12.
        qbdy2 = model.add_qbdy2(sid, eid, qfluxs, comment='qbdhy2')
        qbdy2.raw_fields()

        q0 = 14.
        cntrlnd = 57
        eids = [1, 2]
        qbdy3 = model.add_qbdy3(sid, q0, cntrlnd, eids, comment='qbdy3')
        qbdy3.raw_fields()

        temperature = 13.3
        model.add_tempd(sid, temperature, comment='tempd')

        fields = ['TEMPD', 101, 1., 102, 2., 103, 3., 104, 4.]
        model.add_card(fields, 'TEMPD')


        temperatures = {
            15 : 37.,
            16 : 38.,
        }
        model.add_temp(sid, temperatures)
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
        model.pop_xref_errors()

        model._verify_bdf(xref=True)
        model.write_bdf(bdf_filename2, encoding=None, size=16,
                        is_double=False,
                        interspersed=False,
                        enddata=None, close=False)
        model.write_bdf(bdf_filename3, encoding=None, size=16,
                        is_double=True,
                        interspersed=False,
                        enddata=None, close=False)
        write_bdf_symmetric(model, bdf_filename4, encoding=None, size=8,
                            is_double=False,
                            enddata=None, close=False, plane='xz')
        #model.cross_reference()

        #print(bdf_filename.getvalue())

        bdf_filename2.seek(0)
        model2 = read_bdf(bdf_filename2, xref=False, log=log, debug=False)
        model2.safe_cross_reference()
        save_load_deck(model, punch=False, run_renumber=False, run_test_bdf=False)

    def test_thermal_2(self):
        """tests TABLEHT, TABLEH1"""
        model = BDF(debug=False, log=None, mode='msc')
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
        pconv = model.add_pconv(pconv_id, mid, form=0, expf=0.0, ftype=0,
                                tid=tid_tableht, chlen=None, gidin=None,
                                ce=0, e1=None, e2=None, e3=None, comment='pconv')
        pconv.raw_fields()

        # Every surface to which free convection is to be applied must
        # reference a PCONV entry. PCONV is referenced on the CONV Bulk Data entry.
        eid = 1
        ta = 1
        conv = model.add_conv(eid, pconv_id, ta, film_node=0, cntrlnd=0, comment='conv')
        conv.raw_fields()

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
        phbdy = model.add_phbdy(pid_phybdy, af=None, d1=None, d2=None, comment='phbdy')
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [0., 1., 0.])

        chbdye.raw_fields()
        chbdyg.raw_fields()
        chbdyp.raw_fields()

        model.validate()
        save_load_deck(model)



if __name__ == '__main__':   # pragma: no cover
    unittest.main()
