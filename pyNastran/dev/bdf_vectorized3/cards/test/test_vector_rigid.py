import unittest
import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF, BDFCard # , RBE1, RBE2, RBE3, RROD
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck
from pyNastran.bdf.field_writer_8 import print_card_8

class TestRigid(unittest.TestCase):

    def test_rbe3_01(self):
        model = BDF(debug=False)
        lines = [
            'rbe3,6, ,3,123456,1.0,123456,41,4,+rbe3',
            '+rbe3,alpha,2.0e-4',
        ]
        card = model._process_card(lines)
        card = BDFCard(card)
        rbe = model.rbe3
        rbe_id = rbe.add_card(card)
        #fields = rbe.raw_fields()
        #msg = print_card_8(fields).rstrip()
        msg = rbe.write().rstrip()
        lines_expected = [
            'RBE3           6               3  123456      1.  123456      41       4',
            '           ALPHA   .0002'
        ]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

        if 0:
            dependent_nid_to_components = check_rbe(rbe)
            #print('dependent_nid_to_components = ', dependent_nid_to_components)
            assert dependent_nid_to_components == {3: '123456'}, dependent_nid_to_components

    def test_rbe3_02(self):
        """RBE3 Gmi/Cmi default"""
        model = BDF(debug=None, log=None, mode='msc')
        model.add_grid(1, [0., 0., 0])
        model.add_grid(4, [1., 0., 0])
        model.add_grid(5, [0., 1., 0])
        model.add_grid(6, [1., 1., 0])
        rbe3_id = model.add_rbe3(eid=1, refgrid=1, refc=1, weights=[.1, .5, 3.], comps=['123']*3,
                                 Gmi=None, Cmi=None, Gijs=[4, 5, 6])
        model.rbe3.write()
        save_load_deck(model)

    #-------------------------------------------------------------------------
    def test_rbe2_01(self):
        lines = [
            'RBE2      100045  166007  123456  117752  101899  117766  101898 117748',
            '+         117765  117764  117763  109821  117743  117744  117750 117751',
            '+         117745  117746  101902    1.-6',
        ]
        model = BDF(debug=False)
        card = model._process_card(lines)
        card = BDFCard(card)
        rbe2 = model.rbe2
        rbe_id = rbe2.add_card(card)
        rbe2.write(size=16)
        #fields = rbe.raw_fields()
        msg = rbe2.write(size=8).rstrip()
        #msg = print_card_8(fields).rstrip()
        #print(msg)
        lines_expected = [
            'RBE2      100045  166007  123456  117752  101899  117766  101898  117748',
            '          117765  117764  117763  109821  117743  117744  117750  117751',
            '          117745  117746  101902 .000001'
        ]
        lines_actual = msg.rstrip().split('\n')
        msg2 = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg2 += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg2)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg2)

        if 0:
            dependent_nid_to_components = check_rbe(rbe)
            expected = {117763: '123456', 117764: '123456', 117765: '123456', 117766: '123456',
                        101898: '123456', 101899: '123456', 101902: '123456', 117743: '123456',
                        117744: '123456', 117745: '123456', 117746: '123456', 117748: '123456',
                        117750: '123456', 117751: '123456', 117752: '123456', 109821: '123456'}
            assert dependent_nid_to_components == expected, dependent_nid_to_components

    def test_rbe2_02(self):
        lines = [
            'RBE2      100045  166007  123456  117752  101899  117766  101898 117748',
            '+         117765  117764  117763  109821  117743  117744  117750 117751',
            '+         117745  117746  101902    ',
        ]
        model = BDF(debug=False)
        card = model._process_card(lines)
        card = BDFCard(card)
        rbe2 = model.rbe2
        rbe2_id = rbe2.add_card(card)
        msg = rbe2.write().rstrip()

        eid = 2
        cm = 123
        gn = 10
        Gmi = [3, 4, 5]

        model.add_rbe2(eid, gn, cm, Gmi, alpha=1.0, tref=2.0, comment='cat', validate=False)
        model.parse_cards()

        #fields = rbe.raw_fields()
        #msg = print_card_8(fields).rstrip()
        lines_expected = [
            'RBE2      100045  166007  123456  117752  101899  117766  101898  117748',
            '          117765  117764  117763  109821  117743  117744  117750  117751',
            '          117745  117746  101902'
        ]
        lines_actual = msg.rstrip().split('\n')
        msg2 = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg2 += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg2)

        if 0:
            dependent_nid_to_components = check_rbe(rbe)
            expected = {117763: '123456', 117764: '123456', 117765: '123456', 117766: '123456',
                        101898: '123456', 101899: '123456', 101902: '123456', 117743: '123456',
                        117744: '123456', 117745: '123456', 117746: '123456', 117748: '123456',
                        117750: '123456', 117751: '123456', 117752: '123456', 109821: '123456'}
            assert dependent_nid_to_components == expected, dependent_nid_to_components

        if 0:
            model = BDF(debug=None, log=None, mode='msc')
            eid = rbe.eid
            gn = rbe.gn
            cm = rbe.cm
            Gmi = rbe.Gmi
            alpha = rbe.alpha
            model.add_rbe2(eid, gn, cm, Gmi, alpha=alpha, comment='rbe2')
        nids = [117752, 101899, 117766, 101898, 117748, 117765, 117764, 117763,
                109821, 117743, 117744, 117750, 117751, 117745, 117746, 101902,
                166007,
                3, 4, 5, 10]
        for nid in nids:
            model.add_grid(nid, [0., 0., 0.])
        save_load_deck(model)

    def test_rbe2_02b(self):
        """similar to test_rbe2_02, except with alpha"""
        model = BDF(debug=None, log=None, mode='msc')
        eid = 10
        gn = 100
        cm = '123'
        Gmi = [20, 21, 30]
        alpha = 1.
        model.add_rbe2(eid, gn, cm, Gmi, alpha=alpha, comment='rbe2')
        nids = [gn] + Gmi
        for nid in nids:
            model.add_grid(nid, [0., 0., 0.])
        save_load_deck(model)

    #-------------------------------------------------------------------------
    def test_rbe1_01(self):
        lines = [
            'RBE1    10201   10201   123     10202   456',
            '           UM   10201   456     10202   123',
        ]

        model = BDF(debug=False)
        card = model._process_card(lines)
        #print(print_card_8(card))
        card = BDFCard(card)
        rbe1 = model.rbe1
        rbe1_id = rbe1.add_card(card)
        #fields = rbe.raw_fields()
        #msg = print_card_8(fields).rstrip()
        msg = rbe1.write().rstrip()

        lines_expected = [
            'RBE1       10201   10201     123   10202     456',
            '              UM   10201     456   10202     123'
        ]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

        if 0:
            dependent_nid_to_components = check_rbe(rbe)
            assert dependent_nid_to_components == {10201: '456', 10202: '123'}, dependent_nid_to_components

    def test_rbe1_02(self):
        lines = [
            'RBE1        1001    1000  123456',
            '              UM    1002     123    1003     123    1004     123',
            '                    1005     123    1006     123    1008     123',
            '                    1009     123    1010     123    1011     123',
            '                    1012     123',
        ]
        model = BDF(debug=False)
        card = model._process_card(lines)
        #print(print_card_8(card))
        card = BDFCard(card)
        rbe1 = model.rbe1
        rbe1_id = rbe1.add_card(card)
        #fields = rbe.raw_fields()
        #msg = print_card_8(fields).rstrip()
        msg = rbe1.write().rstrip()

        lines_expected = [
            'RBE1        1001    1000  123456',
            '              UM    1002     123    1003     123    1004     123',
            '                    1005     123    1006     123    1008     123',
            '                    1009     123    1010     123    1011     123',
            '                    1012     123',
        ]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

        if 0:
            dependent_nid_to_components = check_rbe(rbe)
            assert dependent_nid_to_components == {1002: '123', 1003: '123', 1004: '123', 1005: '123', 1006: '123', 1008: '123', 1009: '123', 1010: '123', 1011: '123', 1012: '123'}, dependent_nid_to_components

        model = BDF(debug=None, log=None, mode='msc')
        if 0:
            eid = rbe.eid
            Gni = rbe.Gni
            Cni = rbe.Cni
            Gmi = rbe.Gmi
            Cmi = rbe.Cmi
            alpha = rbe.alpha
            model.add_rbe1(eid, Gni, Cni, Gmi, Cmi, alpha=alpha, comment='rbe1')
        nids = [1000, 1002, 1003, 1004, 1005, 1006, 1008, 1009, 1010, 1011, 1012]
        for nid in nids:
            model.add_grid(nid, [0., 0., 0.])
        save_load_deck(model)

    def test_rbe1_03(self):
        lines = [
            'rbe1,46,3,123456, , , , , ,+rbe46',
            '+rbe46,UM,4,123456,5,123456,2.0-6'
        ]
        model = BDF(debug=False)
        card = model._process_card(lines)
        card = BDFCard(card)
        rbe1 = model.rbe1
        rbe1_id = rbe1.add_card(card)
        #fields = rbe.raw_fields()
        #msg = print_card_8(fields).rstrip()
        msg = rbe1.write().rstrip()

        lines_expected = [
            'RBE1          46       3  123456',
            '              UM       4  123456       5  123456 .000002'
        ]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

        if 0:
            dependent_nid_to_components = check_rbe(rbe)
            assert dependent_nid_to_components == {4: '123456', 5: '123456'}, dependent_nid_to_components

    def _test_rsscon(self):
        model = BDF(debug=False)
        eid = 100
        shell_eid = 1
        solid_eid = 2
        rsscon = model.add_rsscon(
            eid, 'ELEM',
            shell_eid=shell_eid, solid_eid=solid_eid,
            a_solid_grids=None, b_solid_grids=None, shell_grids=None,
            comment='rsscon')

        eid = 101
        shell_grids = [31]
        a_solid_grids = [74]
        b_solid_grids = [75]
        model.add_rsscon(
            eid, 'GRID',
            shell_eid=None, solid_eid=None,
            a_solid_grids=a_solid_grids, b_solid_grids=b_solid_grids, shell_grids=shell_grids,
            comment='rsscon')

        eid = 102
        shell_grids = [11, 14]
        a_solid_grids = [12, 15]
        b_solid_grids = [13, 16]
        model.add_rsscon(
            eid, 'GRID',
            shell_eid=None, solid_eid=None,
            a_solid_grids=b_solid_grids, b_solid_grids=b_solid_grids, shell_grids=shell_grids,
            comment='rsscon')

        dependent_nid_to_components = check_rbe(rsscon)
        assert dependent_nid_to_components == {}, dependent_nid_to_components

        save_load_deck(model, punch=True)

    def test_rbar(self):
        """tests an RBAR"""
        model = BDF(debug=False, log=None, mode='msc')
        eid = 100
        nids = [10, 20]
        cna = '123'
        cnb = '12'
        cma = '456'
        cmb = '3456'
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(20, [0., 0., 0.])
        model.add_rbar(eid, nids, cna, cnb, cma, cmb, alpha=0., comment='rbar')
        save_load_deck(model)

    def test_rbar1(self):
        """tests an RBAR1"""
        model = BDF(debug=False, log=None, mode='msc')
        eid = 100
        nids = [10, 20]
        cb = '123'
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(20, [0., 0., 0.])
        rbar1_id = model.add_rbar1(eid, nids, cb, comment='rbar1')
        #rbar1.raw_fields()
        save_load_deck(model)

    def _test_rspline(self):
        """tests an RSPLINE"""
        model = BDF(debug=False)
        eid = 100
        independent_nid = 10
        dependent_nids = [20, 30]
        dependent_components = [4, 3]
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(20, [0., 0., 0.])
        model.add_grid(30, [0., 0., 0.])
        rspline_id = model.add_rspline(eid, independent_nid, dependent_nids,
                                       dependent_components,
                                       diameter_ratio=0.1,
                                       comment='rspline')
        card_lines = [
            #'$RSPLIN	eid 	D/L	G1	G2  	C2  	G3  	C3  	G4',
            'RSPLINE	17601		13	17601	123456	17603	13456	14'
        ]
        model.add_card(card_lines, 'RSPLINE', is_list=False)
        model.rspline.write(size=8)
        model.rspline.write(size=8)
        #rspline.raw_fields()
        save_load_deck(model)

    def test_rbar(self):
        """tests an RBAR"""
        model = BDF(debug=False)
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(20, [0., 0., 0.])
        eid = 2
        ga = 10
        gb = 20
        cna = '123456'
        cnb = ''
        cma = ''
        cmb = ''
        rbar = model.rbar
        rbar_id = model.add_rbar(eid, [ga, gb], cna, cnb, cma, cmb, alpha=0.,
                                 comment='rbar')
        rbar.write(size=16)
        rbar.write(size=8)
        #rbar.raw_fields()
        save_load_deck(model)

    def test_rrod(self):
        model = BDF(debug=False)
        model.add_grid(10, [0., 0., 0.])
        model.add_grid(20, [0., 0., 0.])
        eid = 2
        ga = 10
        gb = 20
        rrod = model.rrod
        rrod_id1 = model.add_rrod(eid, [ga, gb], cma='42', cmb='33')
        #with self.assertRaises(RuntimeError):
            #rrod_.validate()
        rrod_id2 = model.add_rrod(eid, [ga, gb], cma='3', cmb=None, alpha=0.0, comment='')
        rrod_id3 = model.add_rrod(eid, [ga, gb], cma=None, cmb=3, alpha=0.0, comment='')
        rrod_id3 = model.add_rrod(eid, [ga, gb], cma='', cmb=None, alpha=0.0, comment='')

        rrod.write(size=8)
        rrod.write(size=16)
        #rrod.raw_fields()
        save_load_deck(model, run_read_write=False)

    def test_rbe3_update(self):
        model = BDF(debug=False)
        eid = 1
        refgrid = 2
        refc = '123'
        weights = [1.0, 2.0]
        comps = ['123', '456']
        Gijs = [3, 4]
        rbe3_id = model.add_rbe3(eid, refgrid, refc, weights, comps, Gijs,
                                 Gmi=None, Cmi=None, alpha=0.0, comment='')
        rbe3 = model.rbe3
        model.parse_cards()
        assert rbe3.n == 1
        #print(rbe3)
        rbe3.refgrid = np.array([42])
        rbe3.write()
        rbe3.write_8()
        rbe3.write_16()
        #print(rbe3)
        #rbe3.get_field(4)
        #rbe3.update_field(2, 3)
        #self.assertRaises(IndexError):
        #rbe3.update_field(4, 6)
        save_load_deck(model)

def check_rbe(rbe):
    """simple RBE checks"""
    model = BDF(debug=None)
    model.rigid_elements[rbe.eid] = rbe
    node_ids = []
    for nid in rbe.independent_nodes + rbe.dependent_nodes:
        node_ids.append(nid)
    rigid_elements = model.get_rigid_elements_with_node_ids(node_ids)

    if rbe.type not in ['RSSCON']:
        # supported
        assert rbe.eid in rigid_elements, rbe
    else:
        # not supported
        assert rbe.eid not in rigid_elements, rbe
    dependent_nid_to_components = model.get_dependent_nid_to_components(mpc_id=None, stop_on_failure=True)
    return dependent_nid_to_components

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
