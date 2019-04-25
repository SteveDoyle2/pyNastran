import unittest
from pyNastran.bdf.bdf import BDF, BDFCard, RBE1, RBE2, RBE3
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.test.utils import save_load_deck

bdf = BDF(debug=False)
class TestRigid(unittest.TestCase):

    def test_rbe3_01(self):
        lines = [
            'rbe3,6, ,3,123456,1.0,123456,41,4,+rbe3',
            '+rbe3,alpha,2.0e-4',
        ]
        card = bdf._process_card(lines)
        card = BDFCard(card)
        rbe = RBE3.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()
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
        dependent_nid_to_components = check_rbe(rbe)
        #print('dependent_nid_to_components = ', dependent_nid_to_components)
        assert dependent_nid_to_components == {3: '123456'}, dependent_nid_to_components

    def test_rbe3_02(self):
        """RBE3 Gmi/Cmi default"""
        model = BDF()
        model.add_grid(1, [0.,0.,0])
        model.add_grid(4, [1.,0.,0])
        model.add_grid(5, [0.,1.,0])
        model.add_grid(6, [1.,1.,0])
        rbe3 = model.add_rbe3(eid=1, refgrid=1, refc=1, weights=[.1, .5, 3.], comps=['123']*3,
                              Gmi=None, Cmi=None, Gijs=[4, 5, 6])
        rbe3.write_card()

    #-------------------------------------------------------------------------
    def test_rbe2_01(self):
        lines = [
            'RBE2      100045  166007  123456  117752  101899  117766  101898 117748',
            '+         117765  117764  117763  109821  117743  117744  117750 117751',
            '+         117745  117746  101902    1.-6',
        ]
        card = bdf._process_card(lines)
        card = BDFCard(card)
        rbe = RBE2.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()
        #print(msg)
        lines_expected = [
            'RBE2      100045  166007  123456  117752  101899  117766  101898  117748',
            '          117765  117764  117763  109821  117743  117744  117750  117751',
            '          117745  117746  101902 .000001'
        ]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

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
        card = bdf._process_card(lines)
        card = BDFCard(card)
        rbe = RBE2.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()
        lines_expected = [
            'RBE2      100045  166007  123456  117752  101899  117766  101898  117748',
            '          117765  117764  117763  109821  117743  117744  117750  117751',
            '          117745  117746  101902      0.'
        ]
        lines_actual = msg.rstrip().split('\n')
        msg = '\n%s\n\n%s\n' % ('\n'.join(lines_expected), msg)
        msg += 'nlines_actual=%i nlines_expected=%i' % (len(lines_actual), len(lines_expected))
        self.assertEqual(len(lines_actual), len(lines_expected), msg)
        for actual, expected in zip(lines_actual, lines_expected):
            self.assertEqual(actual, expected, msg)

        dependent_nid_to_components = check_rbe(rbe)
        expected = {117763: '123456', 117764: '123456', 117765: '123456', 117766: '123456',
                    101898: '123456', 101899: '123456', 101902: '123456', 117743: '123456',
                    117744: '123456', 117745: '123456', 117746: '123456', 117748: '123456',
                    117750: '123456', 117751: '123456', 117752: '123456', 109821: '123456'}
        assert dependent_nid_to_components == expected, dependent_nid_to_components

    #-------------------------------------------------------------------------
    def test_rbe1_01(self):
        lines = [
            'RBE1    10201   10201   123     10202   456',
            '           UM   10201   456     10202   123',
        ]

        card = bdf._process_card(lines)
        #print(print_card_8(card))
        card = BDFCard(card)
        rbe = RBE1.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()

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
        card = bdf._process_card(lines)
        #print(print_card_8(card))
        card = BDFCard(card)
        rbe = RBE1.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()

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

        dependent_nid_to_components = check_rbe(rbe)
        assert dependent_nid_to_components == {1002: '123', 1003: '123', 1004: '123', 1005: '123', 1006: '123', 1008: '123', 1009: '123', 1010: '123', 1011: '123', 1012: '123'}, dependent_nid_to_components

    def test_rbe1_03(self):
        lines = [
            'rbe1,46,3,123456, , , , , ,+rbe46',
            '+rbe46,UM,4,123456,5,123456,2.0-6'
        ]
        card = bdf._process_card(lines)
        card = BDFCard(card)
        rbe = RBE1.add_card(card)
        fields = rbe.raw_fields()
        msg = print_card_8(fields).rstrip()

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


        dependent_nid_to_components = check_rbe(rbe)
        assert dependent_nid_to_components == {4: '123456', 5: '123456'}, dependent_nid_to_components

    def test_rsscon(self):
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

        save_load_deck(model, punch=True, run_save_load_hdf5=True)



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
