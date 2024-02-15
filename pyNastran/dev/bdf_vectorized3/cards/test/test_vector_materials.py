"""defines various Material card tests"""
import unittest
import numpy as np

from cpylog import get_logger
from pyNastran.dev.bdf_vectorized3.bdf import BDF, BDFCard #, MAT1, MAT8, MAT11
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck
from pyNastran.bdf.field_writer_8 import print_card_8


class TestMaterials(unittest.TestCase):
    """tests MAT1"""
    def test_mat1_01(self):
        """tests MAT1 initialization from a BDFCard"""
        #
        #MAT5           1    700.    300.    900.    400.    200.    600.     90.+
        #+             .1
        mid = 1
        E = 2e7
        G = 3e7
        nu = 0.4
        rho = 1e-50
        fields = ['MAT1', mid, E, G, nu, rho]

        card = BDFCard(fields)

        model = BDF(debug=False)
        mat1_id = model.mat1.add_card(card)
        mat1 = model.mat1
        mat1.parse_cards()
        self.assertEqual(mid, mat1.material_id)
        self.assertEqual(E, mat1.E)
        self.assertEqual(G, mat1.G)
        self.assertEqual(nu, mat1.nu)
        self.assertEqual(rho, mat1.rho)

        size = 8
        msg = mat1.write(size, 'dummy')
        self.assertEqual(msg,
                         'MAT1           1    2.+7    3.+7      .4   1.-50\n')

        size = 16
        expected = 'MAT1*                  1       20000000.       30000000.              .4\n*                  1.-50\n'
        actual = mat1.write(size, is_double=False)
        msg = 'actual=\n%s\nexpected=\n%s' % (actual, expected)
        self.assertEqual(actual, expected, msg)

    def test_mat1_02(self):
        """tests MAT1, MATT1"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        mat1 = model.mat1
        matt1 = model.matt1

        mid = 10
        #k = 1000.
        E = 3.0e7
        G = 4.0e6
        nu = 0.2
        mat1_id = model.add_mat1(mid, E, G, nu, comment='mat1')
        mat1.write(size=16, is_double=False)
        mat1.validate()

        e_table = 1
        g_table = 2
        nu_table = 3
        rho_table = 4
        a_table = 4
        ge_table = 4
        st_table = 4
        sc_table = 4
        ss_table = 4
        RUN_MATTX = False
        matt1_id = model.add_matt1(mid, e_table, g_table, nu_table, rho_table,
                                   a_table, ge_table, st_table, sc_table, ss_table,
                                   comment='matt1')
        #matt1.validate()

        if RUN_MATTX:
            x = np.linspace(1., 10.)
            y = np.sin(x) + 5.
            tablem1_id = model.add_tablem1(1, x, y, comment='tablem1')
            tablem1.write()

            x1 = 1.0
            tablem2_id = model.add_tablem2(2, x1, x, y, comment='tablem2')
            tablem2.write()

            x2 = 2.0
            tablem3_id = model.add_tablem3(3, x1, x2, x, y, comment='tablem3')
            tablem3.write()

            #x1 = 1.0
            #x2 = 2.0
            x3 = 3.0
            x4 = 4.0
            a = [5.0]
            tablem4_id = model.add_tablem4(4, x1, x2, x3, x4, a, comment='tablem4')
            tablem4.write()

        model.validate()
        model.cross_reference()
        matt1.write(size=16, is_double=False)

        save_load_deck(model)

    def _test_creep(self):
        """tests MAT1/CREEP"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        mid = 10
        #k = 1000.
        E = 3.0e7
        G = 4.0e6
        nu = 0.2
        mat1_id = model.add_mat1(mid, E, G, nu, comment='mat1')
        mat1 = model.mat1
        mat1.write(size=16, is_double=False)

        T0 = 42.
        exp = 1.
        form = 'cat'
        tidkp = 42
        tidcp = 43
        tidcs = 44
        thresh = 6.
        Type = 7
        a = 8.
        b = 9.
        c = 10.
        d = 11.
        e = 12.
        f = 13.
        g = 14.
        creep_id = model.add_creep(
            mid, T0, exp, form, tidkp, tidcp, tidcs, thresh, Type,
            a, b, c, d, e, f, g,
            comment='creep')
        model.pop_parse_errors()
        #creep.raw_fields()
        model.cross_reference()
        #model.uncross_reference()
        unused_model2 = save_load_deck(model)

    def test_mat2_01(self):
        """tests MAT2, MATT2"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        matt2 = model.matt2
        mid = 10
        G11 = G22 = G12 = G13 = G22 = G23 = G33 = 1.0
        #nuxth = nuthz = nuzx = 0.3
        mat2_D_id = model.add_mat2(mid, G11, G12, G13, G22, G23, G33, rho=0.,
                              a1=None, a2=None, a3=None, tref=0.,
                              ge=0., St=None, Sc=None, Ss=None,
                              mcsid=None, comment='mat2')
        #mat2.raw_fields()
        model.mat2.write(size=16, is_double=False)
        #mat2.validate()

        g11_table = 1
        g12_table = 2
        g13_table = 3
        g22_table = 4
        g23_table = 5
        g33_table = 6
        rho_table = 7
        a1_table = 8
        a2_table = 9
        a3_table = 10
        ge_table = 11
        st_table = 12
        sc_table = 13
        ss_table = 14
        RUN_MATTX = False
        matt2_id = model.add_matt2(mid, g11_table, g12_table, g13_table, g22_table,
                                   g23_table, g33_table, rho_table,
                                   a1_table, a2_table, a3_table,
                                   ge_table, st_table, sc_table, ss_table,
                                   comment='matt2')
        #matt2.validate()

        if RUN_MATTX:
            x = np.linspace(1., 10.)
            y = np.sin(x) + 5.
            for tid in [1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]:
                tablem1_id = model.add_tablem1(tid, x, y, comment='tablem1')
                tablem1.write()

                x1 = 1.0
                tablem2_id = model.add_tablem2(2, x1, x, y, comment='tablem2')
                tablem2.write()

                x2 = 2.0
                tablem3_id = model.add_tablem3(3, x1, x2, x, y, comment='tablem3')
                tablem3.write()

                #x1 = 1.0
                #x2 = 2.0
                x3 = 3.0
                x4 = 4.0
                a = [5.0]
                tablem4_id = model.add_tablem4(4, x1, x2, x3, x4, a, comment='tablem4')
                tablem4.write()

        model.validate()
        model.cross_reference()
        matt2.write(size=16, is_double=False)
        save_load_deck(model)

    def test_mat3_01(self):
        """tests MAT3"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        matt3 = model.matt3
        mat3 = model.mat3
        mid = 10
        ex = 3.0e7
        eth = 6.0e7
        ez = 6e4
        nuxth = nuthz = nuzx = 0.3
        mat3_id = model.add_mat3(mid, ex, eth, ez, nuxth, nuthz, nuzx, rho=0.0,
                                 gzx=None, ax=0.,
                                 ath=0., az=0.,
                                 tref=0., ge=0.,
                                 comment='mat3')
        mat3 = model.mat3
        mat3.write(size=16, is_double=False)
        mat3.rho
        #mat3.raw_fields()
        mat3.validate()

        matt3_id = model.add_matt3(
            mid, ex_table=1, eth_table=2, ez_table=3,
            nuth_table=4, nuxz_table=1, rho_table=1,
            gzx_table=1, ax_table=1, ath_table=1,
            az_table=1, ge_table=1, comment='matt3')
        #matt3.validate()
        matt3.write()

        x = np.linspace(1., 10.)
        y = np.sin(x) + 5.
        tablem1 = model.add_tablem1(1, x, y, comment='tablem1')
        tablem1.write_card()

        x1 = 1.0
        tablem2 = model.add_tablem2(2, x1, x, y, comment='tablem2')
        tablem2.write_card()

        x2 = 2.0
        tablem3 = model.add_tablem3(3, x1, x2, x, y, comment='tablem3')
        tablem3.write_card()

        #x1 = 1.0
        #x2 = 2.0
        x3 = 3.0
        x4 = 4.0
        a = [5.0]
        tablem4 = model.add_tablem4(4, x1, x2, x3, x4, a, comment='tablem4')
        tablem4.write_card()

        model.validate()
        model.pop_parse_errors()
        #model.cross_reference()
        #matt3.write(size=16, is_double=False)

        save_load_deck(model, run_remove_unused=False)

    def test_mat4_01(self):
        """tests MAT4, MATT4"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        matt4 = model.matt4
        mid = 10
        k = 1000.
        mat4_id = model.add_mat4(mid, k, cp=0.0, rho=1.0, H=None, mu=None,
                                 hgen=1.0, ref_enthalpy=None, tch=None, tdelta=None, qlat=None,
                                 comment='mat4')
        #mat4.raw_fields()
        mat4 = model.mat4
        mat4.write(size=16, is_double=False)
        mat4.validate()

        k_table = 1
        cp_table = 2
        h_table = 3
        mu_table = 4
        hgen_table = 3
        matt4_id = model.add_matt4(mid, k_table, cp_table, h_table, mu_table,
                                   hgen_table, comment='matt4')
        #matt4.validate()
        matt4.write()

        x = np.linspace(1., 10.)
        y = np.sin(x) + 5.
        tablem1 = model.add_tablem1(1, x, y, comment='tablem1')
        tablem1.write_card()

        x1 = 1.0
        tablem2 = model.add_tablem2(2, x1, x, y, comment='tablem2')
        tablem2.write_card()

        x2 = 2.0
        tablem3 = model.add_tablem3(3, x1, x2, x, y, comment='tablem3')
        tablem3.write_card()

        #x1 = 1.0
        #x2 = 2.0
        x3 = 3.0
        x4 = 4.0
        a = [5.0]
        tablem4 = model.add_tablem4(4, x1, x2, x3, x4, a, comment='tablem4')
        tablem4.write_card()

        model.validate()
        model.cross_reference()
        matt4.write(size=16, is_double=False)

        save_load_deck(model)

    def test_mat5_01(self):
        """tests MAT5, MATT5"""
        #
        #MAT5           1    700.    300.    900.    400.    200.    600.     90.+
        #+             .1
        log = get_logger(level='warning')
        model = BDF(log=log)
        matt5 = model.matt5
        mat5 = model.mat5

        mid = 10
        mat5_id = model.add_mat5(mid, kxx=0., kxy=0., kxz=0., kyy=0., kyz=0.,
                                 kzz=0., cp=0.,
                                 rho=1., hgen=1., comment='mat5')
        mat5 = model.mat5
        #model.parse_cards()

        mat5.parse_cards()

        mat5.k()
        #mat5.get_density()
        mat5.rho
        mat5.write(size=16, is_double=False)
        mat5.validate()

        kxx_table = 1
        kxy_table = 2
        kxz_table = 3
        kyy_table = 4
        kyz_table = 4
        kzz_table = 4
        cp_table = 4
        hgen_table = 4
        matt5_id = model.add_matt5(mid, kxx_table, kxy_table, kxz_table,
                                   kyy_table, kyz_table, kzz_table, cp_table, hgen_table,
                                   comment='matt5')
        #matt5.validate()
        matt5.write()

        x = np.linspace(1., 10.)
        y = np.sin(x) + 5.
        tablem1 = model.add_tablem1(1, x, y, comment='tablem1')
        tablem1.write_card()

        x1 = 1.0
        tablem2 = model.add_tablem2(2, x1, x, y, comment='tablem2')
        tablem2.write_card()

        x2 = 2.0
        tablem3 = model.add_tablem3(3, x1, x2, x, y, comment='tablem3')
        tablem3.write_card()

        #x1 = 1.0
        #x2 = 2.0
        x3 = 3.0
        x4 = 4.0
        a = [5.0]
        tablem4 = model.add_tablem4(4, x1, x2, x3, x4, a, comment='tablem4')
        tablem4.write_card()

        model.validate()
        model.cross_reference()
        matt5.write(size=16, is_double=False)

        save_load_deck(model, run_remove_unused=False)

    def test_mat8_01(self):  # should fail...
        """tests MAT8"""
        #lines = [  # fails???
        #    'MAT8*    4700007        1675.47         1675.47         .33             *   LHIG',
        #    '*   LHIG28.2            210000.         78000.                          *   LHIH',
        #    '*   LHIH1.32-5          1.32-5          75.             1.943           *   LHII',
        #    '*   LHII1.943           1.943           1.943           3.35',
        #]
        lines = [  # fails
            'MAT8*    4700010        2.83+6          1.14+6          .55             *   LHIJ',
            '*   LHIJ717000.         285194.         285194.                         *   LHIK',
            '*   LHIK9.17-6          2.606-5         70.                             *   LHIL',
            '*   LHIL',
        ]
        lines_expected = [
            'MAT8*            4700010        2830000.        1140000.             .55',
            '*                717000.         285194.         285194.',
            '*              .00000917       .00002606             70.',
            '*',
        ]

        log = get_logger(level='warning')
        model = BDF(log=log)
        card = model._process_card(lines)
        #print(print_card_8(card))
        cardi = BDFCard(card)

        mat8 = model.mat8
        card2 = mat8.add_card(cardi)

        #fields = card2.raw_fields()
        #msg = print_card_8(fields)
        size = 16
        msg = mat8.write()
        msg = mat8.write(size, 'dummy')

        lines_actual = msg.rstrip().split('\n')
        msg2 = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg2 += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        #print(msg2)
        self.assertEqual(len(lines_actual), len(lines_expected), msg2)
        for actual, expected in zip(lines_actual, lines_expected):
            msg = 'actual   = %r\n' % actual
            msg += 'expected = %r' % expected
            self.assertEqual(actual, expected, msg)

    def test_mat8_02(self):
        """tests MAT8, MATT8"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        mid = 10
        e11 = 3.0e7
        e22 = 6.0e7
        nu12 = 0.3
        mat8_id = model.add_mat8(mid, e11, e22, nu12)
        mat8 = model.mat8
        matt8 = model.matt8
        mat8.write(size=16, is_double=False)
        mat8.validate()

        RUN_MATTX = False
        matt8_id = model.add_matt8(
            mid, e1_table=1, e2_table=2, nu12_table=3,
            g12_table=4, g1z_table=1, g2z_table=1, rho_table=1,
            a1_table=1, a2_table=1,
            xt_table=1, xc_table=1, yt_table=1, yc_table=1,
            s_table=1, ge_table=1, f12_table=1, comment='matt8')
        matt8.validate()
        matt8.write()

        if RUN_MATTX:
            x = np.linspace(1., 10.)
            y = np.sin(x) + 5.
            tablem1_id = model.add_tablem1(1, x, y, comment='tablem1')
            tablem1.write()

            x1 = 1.0
            tablem2_id = model.add_tablem2(2, x1, x, y, comment='tablem2')
            tablem2.write()

            x2 = 2.0
            tablem3_id = model.add_tablem3(3, x1, x2, x, y, comment='tablem3')
            tablem3.write()

            #x1 = 1.0
            #x2 = 2.0
            x3 = 3.0
            x4 = 4.0
            a = [5.0]
            tablem4_id = model.add_tablem4(4, x1, x2, x3, x4, a, comment='tablem4')
            tablem4.write()

        model.validate()
        model.cross_reference()
        if RUN_MATTX:
            matt8.write(size=16, is_double=False)
        save_load_deck(model)

    def test_mat9(self):
        """tests MAT9"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        matt9 = model.matt9
        mid = 10
        #e11 = 3.0e7
        #e22 = 6.0e7
        #nu12 = 0.3
        mat9_id = model.add_mat9(mid, G11=0., G12=0., G13=0., G14=0., G15=0.,
                                 G16=0., G22=0., G23=0., G24=0.,
                                 G25=0., G26=0., G33=0., G34=0.,
                                 G35=0., G36=0., G44=0., G45=0.,
                                 G46=0., G55=0., G56=0., G66=0.,
                                 rho=0., A=None, tref=0., ge=0.,
                                 comment='mat9')
        mat9 = model.mat9
        mat9.write(size=16, is_double=False)
        mat9.rho
        mat9.validate()

        matt9_id = model.add_matt9(
            mid,
            g11_table=1, g12_table=1, g13_table=1, g14_table=1, g15_table=1, g16_table=1,
            g22_table=1, g23_table=1, g24_table=1, g25_table=1, g26_table=1,
            g33_table=1, g34_table=1, g35_table=1, g36_table=1,
            g44_table=1, g45_table=1, g46_table=1,
            g55_table=1, g56_table=1, g66_table=1,
            rho_table=2,
            a1_table=3, a2_table=3, a3_table=3, a4_table=3, a5_table=3, a6_table=3,
            ge_table=4, comment='')
        #matt9.validate()

        #if RUN_MATTX:
        model.validate()
        model.cross_reference()
        matt9.write()
        #matt8.write(size=16, is_double=False)
        save_load_deck(model)

    def test_mat11_01(self):
        """tests MAT11"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        lines = [
            'MAT11          1    1.+75000000. 700000.      .1     .13     .267000000.+',
            '+       9000000.3000000.      .1    1.-5    7.-6    8.-6     50.',
        ]
        lines_expected = [
            #'MAT11          1    1.+75000000. 700000..1000000.1300000.26000007000000.'
            #'        9000000.3000000.          .00001 .000007 .000008     50.'
            'MAT11          1    1.+75000000. 700000.      .1     .13     .267000000.',
            '        9000000.3000000.      .1  .00001 .000007 .000008     50.'
        ]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        mat11 = model.mat11
        mat_id = mat11.add_card(cardi)

        #fields = mat.raw_fields()
        #msg = print_card_8(fields)
        #f = StringIO.StringIO()
        size = 8
        msg = mat11.write(size, 'dummy')
        #msg = f.getvalue()
        #print(msg)

        lines_actual = msg.rstrip().split('\n')
        msg2 = '\n%s\n\n%s' % ('\n'.join(lines_expected), msg)
        msg2 += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        #print(msg)
        self.assertEqual(len(lines_actual), len(lines_expected), msg2)
        for actual, expected in zip(lines_actual, lines_expected):
            msg2 = '\nactual   = %r\n' % actual
            msg2 += 'expected =  %r' % expected
            self.assertEqual(actual, expected, msg2)

        save_load_deck(model, xref='standard', punch=True,
                       run_remove_unused=False)

    def _test_mat3d(self):
        """tests MAT3D"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        mid = 10
        e1 = 1.
        e2 = 2.
        e3 = 3.
        nu12 = 12.
        nu13 = 13.
        nu23 = 23.
        g12 = 112.
        g13 = 113.
        g23 = 123.
        mat3d_id = model.add_mat3d(
            mid,
            e1, e2, e3, nu12, nu13, nu23, g12, g13, g23,
            rho=1.0, comment='mat3d')
        #mat3d.raw_fields()
        mat3d.write(size=8, is_double=False)
        mat3d.write(size=16, is_double=False)
        mat3d.write(size=16, is_double=True)
        save_load_deck(model, xref='standard', punch=True,
                       run_remove_unused=False)
        #mat = MAT11(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23,
                    #rho=0.0, a1=0.0, a2=0.0, a3=0.0, tref=0.0, ge=0.0, comment='')

    def test_mats1(self):
        """tests MATS1"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        mid = 10
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu)

        tid = None
        Type = 'NLELAST'
        h = None
        hr = None
        yf = None
        limit1 = None
        limit2 = None
        unused_mats1_id = model.add_mats1(mid, tid, Type, h, hr, yf, limit1, limit2,
                                       comment='mats1')
        save_load_deck(model, xref='standard', punch=True, run_remove_unused=False)

    def _test_matdmg(self):
        """tests MATDMG"""
        log = get_logger(level='warning')
        model = BDF(log=log)

        mid = 10
        e1 = 1.3e8
        e2 = 7707000.0
        e3 = 7707000.0
        nu12 = 0.234
        nu13 = 0.234
        nu23 = 0.391
        g12 = 4154800.0
        g13 = 4154800.0
        g23 = 2770000.0
        model.add_mat11(mid, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, 0.0, comment="mat11")

        ppf_model = "UD"
        y012 = 22.5
        yc12 = 8526.4
        ys12 = 1918.0
        ys22 = 490.0
        y11limt = 15840.0
        y11limc = 15840.0
        ksit = None
        ksic = None
        b2 = 2.5
        b3 = 0.733
        a = 0.5
        litk = 20000.0
        bigk = 518810.0
        expn = 0.3648
        tau = 0.0025
        adel = 1.0
        plyuni = 0
        tid = None
        hbar = 0.0001
        dmax = 0.999
        pe = 0
        user = None
        r01 = None
        ds = None
        gic = None
        giic = None
        giiic = None

        model.add_matdmg(mid, ppf_model,
                         y012, yc12, ys12, ys22, y11limt, y11limc, ksit, ksic,
                         b2, b3, a, litk, bigk, expn, tau, adel,
                         plyuni, tid, hbar, dmax, pe,
                         user, r01, ds, gic, giic, giiic,
                         comment='matdmg')

        save_load_deck(model, xref='standard', punch=True, run_remove_unused=False)

    def test_multiple_materials(self):
        """tests multiple materials"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        E = 3.0e7
        G = None
        nu = 0.3
        mat1_id = model.add_mat1(1, E, G, nu)
        e11 = e22 = 3.0e7
        nu12 = 0.3
        model.add_mat8(8, e11, e22, nu12)

        model.add_mat4(4, 10.0)
        mat5_id = model.add_mat5(5)
        mat9_id = model.add_mat9(9, G11=0., G12=0., G13=0., G14=0., G15=0.,
                                 G16=0., G22=0.,
                                 G23=0., G24=0.,
                                 G25=0., G26=0.,
                                 G33=0., G34=0.,
                                 G35=0., G36=0.,
                                 G44=0., G45=0.,
                                 G46=0., G55=0.,
                                 G56=0., G66=0.,
                                 rho=0., A=None,
                                 tref=0., ge=0.,
                                 comment='mat9')

        bulk = 0.3
        rho = 0.2
        c = None
        model.add_mat10(10, bulk, rho, c)

        e1 = 1.
        e2 = 2.
        e3 = 3.
        nu13 = 0.3
        nu23 = 0.3
        g12 = 12.
        g13 = 13.
        g23 = 23.
        mat11_id = model.add_mat11(11, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23,
                                   rho=0.0, a1=0.0, a2=0.0,
                                   a3=0.0, tref=0.0, ge=0.0,
                                   comment='mat11')

        #mat1.raw_fields()
        #mat5.raw_fields()
        #mat9.raw_fields()
        #mat11.raw_fields()

        if 0:
            structural_material_ids = model.get_structural_material_ids()
            assert len(structural_material_ids) == 5, structural_material_ids

            thermal_material_ids = model.get_thermal_material_ids()
            assert len(thermal_material_ids) == 2, thermal_material_ids

            mats = model.Materials(1)
            assert len(mats) == 1, mats
            mats = model.Materials([1, 4, 5])
            assert len(mats) == 3, mats

            with self.assertRaises(KeyError):
                model.Material(-1)
            with self.assertRaises(KeyError):
                model.StructuralMaterial(-1)
            with self.assertRaises(KeyError):
                model.ThermalMaterial(-1)
        save_load_deck(model)

    def test_mathe(self):
        """tests the MATHE"""
        log = get_logger(level='warning')
        bdf_model = BDF(log=log)

        mid = 10
        bulk = 1.
        mus = []
        alphas = []
        betas = []
        mooney = []
        sussbat = []
        aboyce = []
        gent = []
        if 0:
            mat_model = 'JUNK'
            mat_mooney = bdf_model.add_mathe(mid, mat_model, bulk, mus, alphas, betas,
                                             mooney, sussbat, aboyce, gent,
                                             rho=0., texp=0., tref=0., ge=0., comment='mathe')
            with self.assertRaises(ValueError):
                mat_mooney.validate()

            #MOONEY, OGDEN, FOAM, ABOYCE, SUSSBAT, ABOYCE, GENT
            mat_mooney.mat_model = 'MOONEY'
            with self.assertRaises(AssertionError):
                mat_mooney.validate()
        #mooney = [
            #c10, c01,
            #c20, c11, c02,
            #c30, c21, c12, c03,
        #]
        mat_model = 'MOONEY'
        mooney = [1., 2., 3., 4., 5., 6., 7., 8., 9.]
        mat_mooney = bdf_model.add_mathe(
            mid, mat_model, bulk, mus, alphas, betas,
            mooney, sussbat, aboyce, gent,
            rho=0., texp=0., tref=0., ge=0., comment='mathe')

        #mat_mooney.validate()
        RUN_MAT = False
        if RUN_MAT:
            mat_model = 'OGDEN'
            mid = 11
            mat_ogden = bdf_model.add_mathe(mid, mat_model, bulk, mus, alphas, betas,
                                            mooney, sussbat, aboyce, gent,
                                            rho=0., texp=0., tref=0., ge=0., comment='mathe')
            mat_ogden.validate()

            mat_model = 'FOAM'
            mid = 12
            mat_foam = bdf_model.add_mathe(mid, mat_model, bulk, mus, alphas, betas,
                                           mooney, sussbat, aboyce, gent,
                                           rho=0., texp=0., tref=0., ge=0., comment='mathe')
            mat_foam.validate()


            mat_model = 'SUSSBAT'
            mid = 13

            sussbat = [1, 'PIG', 3.14] # #(tab1, sstype, relerr)
            mat_sussbat = bdf_model.add_mathe(mid, mat_model, bulk, mus, alphas, betas,
                                              mooney, sussbat, aboyce, gent,
                                              rho=0., texp=0., tref=0., ge=0., comment='mathe')
            mat_sussbat.validate()

            mat_model = 'ABOYCE'
            mid = 14
            aboyce = [1.1, 2.2] # (nkt, n)
            mat_aboyce = bdf_model.add_mathe(mid, mat_model, bulk, mus, alphas, betas,
                                             mooney, sussbat, aboyce, gent,
                                             rho=0., texp=0., tref=0., ge=0., comment='mathe')
            mat_aboyce.validate()

            model = 'GENT'
            mid = 15
            gent = [1.11, 2.22, 3.33, 4.44, 5.55, 6.66, 7.77, 8.88]  # (nkt, n, im, d1, d2, d3, d4, d5)
            mat_gent = bdf_model.add_mathe(mid, model, bulk, mus, alphas, betas,
                                           mooney, sussbat, aboyce, gent,
                                           rho=0., texp=0., tref=0., ge=0., comment='mathe')
            mat_gent.validate()

        #mat_mooney.raw_fields()
        #mat_foam.raw_fields()
        #mat_ogden.raw_fields()
        #mat_sussbat.raw_fields()
        #mat_aboyce.raw_fields()
        #mat_gent.raw_fields()
        bdf_model.setup()

        #bdf_model.pop_parse_errors()

        save_load_deck(bdf_model) # , run_remove_unused=False

    def test_mat10c(self):
        log = get_logger(level='warning')
        model = BDF(log=log)
        mid = 10
        model.add_mat10c(mid, form='REAL', rho_real=0.0, rho_imag=0.0, c_real=0.0, c_imag=0.0, comment='')
        model.setup()
        model.mat10c.write()
        save_load_deck(model, run_remove_unused=False)

    def test_matort(self):
        log = get_logger(level='warning')
        model = BDF(log=log)
        mid = 10
        E1 = 10.
        E2 = 20.
        E3 = 30.
        nu12 = 0.3
        nu23 = 0.2
        nu31 = 0.1
        G12 = 1.
        G23 = 2.
        G31 = 3.
        model.add_matort(mid, E1, E2, E3, nu12, nu23, nu31, G12, G23, G31,
                         rho=0.0, alpha1=0.0, alpha2=0.0, alpha3=0.0, tref=0.0, ge=0.0, comment='')
        model.setup()
        print(model.matort.write())
        save_load_deck(model)

    def test_mathp(self):
        log = get_logger(level='warning')
        model = BDF(log=log)
        mid = 10
        #E1 = 10.
        #E2 = 20.
        #E3 = 30.
        #nu12 = 0.3
        #nu23 = 0.2
        #nu31 = 0.1
        #G12 = 1.
        #G23 = 2.
        #G31 = 3.
        model.add_mathp(mid, a10=0., a01=0., d1=None, rho=0.,
                        av=0., tref=0., ge=0., na=1, nd=1,
                        a20=0., a11=0., a02=0., d2=0., a30=0., a21=0., a12=0., a03=0., d3=0.,
                        a40=0., a31=0., a22=0., a13=0., a04=0., d4=0., a50=0., a41=0.,
                        a32=0., a23=0., a14=0., a05=0., d5=0.,
                        tab1=None, tab2=None, tab3=None, tab4=None, tabd=None, comment='')
        print(model.mathp.write())
        model.setup()
        save_load_deck(model)

    def _test_matg(self):
        """tests the MATG"""
        log = get_logger(level='warning')
        model = BDF(log=log)
        mid = 10
        idmem = 2
        behav = 3
        tabld = 4
        tablu = [50, 51, 52, 53]
        yprs = 6.
        epl = 7.
        gpl = 8.
        matg_id = model.add_matg(mid, idmem, behav, tabld, tablu, yprs, epl, gpl,
                                 gap=0., tab_yprs=None, tab_epl=None, tab_gpl=None, tab_gap=None,
                                 comment='matg')
        matg.validate()
        #matg.raw_fields()
        save_load_deck(model, run_convert=False, run_op2_writer=False, run_remove_unused=False)

    def test_nxstrat(self):
        """tests the NXSTRAT"""
        params = {
            #'AUTO' : 1,
            #'MAXITE' : 30,
            'RTOL' : 0.005,
            #'ATSNEXT' : 3,
            'A' : 1,
            'B' : 2,
            'C' : 3,
            'D' : 4,
            'E' : 5,
            'F' : 6,
            'G' : 7,
            'H' : 8,
        }
        log = get_logger(level='warning')
        model = BDF(log=log)
        nxstrat_id = model.add_nxstrat(42, params)
        #nxstrat.raw_fields()
        save_load_deck(model) # , run_remove_unused=False

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
