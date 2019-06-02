import unittest
from io import StringIO
from numpy import array, allclose

from pyNastran.dev.bdf_vectorized.bdf import BDF


class TestCoords(unittest.TestCase):
    def test_same(self):  # passes
        grids = [
            [1, 0, 0., 0., 1.],
            [2, 0, 0., 1., 0.],
            [3, 0, 1., 0., 0.],
            [4, 0, 1., 1., 1.],
            [5, 0, 1., 1., 0.],
        ]
        grids_expected = grids
        coords = []
        self._get_nodes(grids, grids_expected, coords)

    def test_shift(self):  # passes
        grids = [
            [1, 1, 0., 0., 1.],
            [2, 1, 0., 1., 0.],
            [3, 1, 1., 0., 0.],
            [4, 1, 1., 1., 1.],
            [5, 1, 1., 1., 0.],
        ]
        grids_expected = [
            [1, 1, 1., 1., 2.],
            [2, 1, 1., 2., 1.],
            [3, 1, 2., 1., 1.],
            [4, 1, 2., 2., 2.],
            [5, 1, 2., 2., 1.],
        ]

        coords = [
            # cid,rid, origin,      zaxis,     xaxis
            [1, 0, [1., 1., 1.], [1., 1., 2.], [2., 1., 1.]],
        ]
        self._get_nodes(grids, grids_expected, coords)

    def test_rotate(self):  # passes
        grids = [
            [1, 1, 0., 0., 1.],
            [2, 1, 0., 1., 0.],
            [3, 1, 1., 0., 0.],
            [4, 1, 1., 1., 1.],
            [5, 1, 1., 1., 0.],

        ]
        grids_expected = [
            #     y    z   x
            [1, 1., 1., 0., 0.],
            [2, 1., 0., -1., 0.],
            [3, 1., 0., 0., 1.],
            [4, 1., 1., -1., 1.],
            [5, 1., 0., -1., 1.],
        ]

        coords = [
            # cid, rid, origin,      zaxis,     xaxis
            [1, 0, [0., 0., 0.], [1., 0., 0.], [0., 0., 1.]],
        ]
        self._get_nodes(grids, grids_expected, coords)

    def test_rotate2(self):   # passes
        grids = [
            [1, 1, 0., 0., 1.],  # nid, cid, x,y,z
            [2, 1, 0., 1., 0.],
            [3, 1, 1., 0., 0.],
            [4, 1, 1., 1., 1.],
            [5, 1, 1., 1., 0.],
        ]
        grids_expected = [
            [1, 1, 0., 0., -1.],
            [2, 1, 0., -1., 0.],
            [3, 1, 1., 0., 0.],
            [4, 1, 1., -1., -1.],
            [5, 1, 1., -1., 0.],
        ]

        coords = [
            # cid, rid, origin,     zaxis        xaxis
            [1, 0, [0., 0., 0.], [0., 0., -1.], [1., 0., 0.]],
        ]
        self._get_nodes(grids, grids_expected, coords)

    def test_rotate3(self):  # passes
        grids = [
            [1, 1, 0., 0., 1.],
            [2, 1, 0., 1., 0.],
            [3, 1, 1., 0., 0.],
            [4, 1, 1., 1., 1.],
            [5, 1, 1., 1., 0.],
        ]
        grids_expected = [
            [1, 1, 0., 0., -1.],
            [2, 1, 0., 1., 0.],
            [3, 1, -1., 0., 0.],
            [4, 1, -1., 1., -1.],
            [5, 1, -1., 1., 0.],
        ]

        coords = [
            # cid, rid, origin,      zaxis          xaxis
            [1, 0, [0., 0., 0.], [0., 0., -1.], [-1., 0., 0.]],
        ]
        self._get_nodes(grids, grids_expected, coords)

    def test_rid_1(self):
        grids = [
            [1, 2, 10., 5., 3.],  # nid, cid, x,y,z
            [2, 3, 10., 5., 3.],
        ]
        grids_expected = [
            [1, 1, 11., 6., 4.],
            [2, 1, 11., 6., 4.],
        ]

        coords = [
            # cid, rid, origin,     zaxis        xaxis
            [1, 0, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.]],  # cid=1
            [2, 1, [1., 1., 1.], [1., 1., 2.], [2., 1., 1.]],  # cid=2
            #[2, 1, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.]],  # cid=2,equiv

            [3, 0, [1., 1., 1.], [1., 1., 2.], [2., 1., 1.]],  # cid=3
            #[3, 0, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.]],  # cid=3,equiv
        ]
        self._get_nodes(grids, grids_expected, coords)

    def test_cord1r_01(self):
        lines = ['cord1r,2,1,4,3']

        grids = [
            ['GRID', 4, 0, 0.0, 0., 0.],
            ['GRID', 3, 0, 0.0, 0., 1.],
            ['GRID', 1, 0, 0.0, 1., 0.],
        ]
        card_count = {
            'CORD1R' : 1,
            'GRID' : 3,
        }

        model = BDF(debug=False)
        #model.allocate(card_count)
        model.add_card(lines, 'CORD1R', is_list=False)
        for grid in grids:
            model.add_card(grid, 'GRID', is_list=True)
        model.build()

        size = 8
        coords = model.coords
        self.assertEqual(coords.get_cid_by_coord_id(2), 2)
        self.assertEqual(coords.get_rid_by_coord_id(2), 0)

        bdf_file = StringIO()
        coords.write_card(bdf_file, size=8, is_double=False)

    def test_cord2c_01(self):
        lines_a = [
            'CORD2C*                3               0              0.              0.',
            '*                     0.              0.              0.              1.*',
            '*                     1.              0.              1.'
        ]
        lines_b = [
            'CORD2R         4       3     10.      0.      5.     10.     90.      5.',
            '             10.      0.      6.'
        ]
        card_count = {
            'CORD2C' : 1,
            'CORD2R' : 1,
        }

        model = BDF(debug=False)
        cards = {
            'CORD2C' : ['', lines_a],
            'CORD2R' : ['', lines_b],
        }
        model.add_cards(cards, card_count)
        #model.allocate(cards, card_count)
        card = model.add_card(lines_a, 'CORD2C', is_list=False)
        card = model.add_card(lines_b, 'CORD2R', is_list=False)
        model.build()

        cord2r = model.coords.slice_by_coord_id(3)
        #print(type(cord2r))
        self.assertEqual(cord2r.get_cid_by_coord_id(), 3)
        self.assertEqual(cord2r.get_rid_by_coord_id(), 0)

        cord2r = model.coords.slice_by_coord_id(4)
        #print(type(cord2r))
        self.assertEqual(cord2r.get_cid_by_coord_id(), 4)
        self.assertEqual(cord2r.get_rid_by_coord_id(), 3)

        i = model.coords.get_coord_index_by_coord_id(4)
        T = model.coords.T[i, :, :].reshape(3, 3)
        Ti = T[0, :]
        Tj = T[1, :]
        Tk = T[2, :]
        msg = 'i=%s expected=(0., 0., 1.)'  % Ti
        self.assertTrue(allclose(Ti, array([0., 0., 1.])), msg)
        delta = Tj - array([1., 1., 0.]) / 2**0.5
        self.assertTrue(allclose(Tj, array([1., 1., 0.]) / 2**0.5), str(delta))
        delta = Tk - array([-1., 1., 0.]) / 2**0.5
        self.assertTrue(allclose(Tk, array([-1., 1., 0.]) / 2**0.5), str(delta))

    def test_grid_01(self):
        model = BDF(debug=False)
        card_lines = [
            #['CORD1R', 1, 1, 2, 3],  # fails on self.k
            ['GRID', 1, 0, 0., 0., 0.],
            ['GRID', 2, 0, 1., 0., 0.],
            ['GRID', 4, 0, 1., 2., 3.],
        ]
        #card_count = {
            #'GRID': 3,
        #}
        cards, card_count = model.add_cards_lines(card_lines)
        model.allocate(card_count, cards)
        for card in cards:
            model.add_card(card, card[0], comment='comment', is_list=True)
        model.build()
        #+------+-----+----+----+----+----+----+----+------+
        #|   0  |  1  | 2  | 3  | 4  | 5  |  6 | 7  |  8   |
        #+======+=====+====+====+====+====+====+====+======+
        #| GRID | NID | CP | X1 | X2 | X3 | CD | PS | SEID |
        #+------+-----+----+----+----+----+----+----+------+
        node = model.grid.slice_by_node_id(4)
        #self.assertEqual(node.get_field(1), 4)
        #self.assertEqual(node.get_field(2), 0)
        #self.assertEqual(node.get_field(3), 1.)
        #self.assertEqual(node.get_field(4), 2.)
        #self.assertEqual(node.get_field(5), 3.)

        #node.update_field(1, 5)
        #node.update_field(2, 6)
        #node.update_field(3, 7.)
        #node.update_field(4, 8.)
        #node.update_field(5, 9.)
        #with self.assertRaises(KeyError):
            #node.update_field(9, 'dummy')

        #self.assertEqual(node.get_field(1), 5)
        #self.assertEqual(node.get_field(2), 6)
        #self.assertEqual(node.get_field(3), 7.)
        #self.assertEqual(node.get_field(4), 8.)
        #self.assertEqual(node.get_field(5), 9.)
        #with self.assertRaises(KeyError):
            #node.get_field(9)

    def test_cord1r_02(self):
        model = BDF(debug=False)
        unused_card_count = {
            'CORD1R' : 1,
            'GRID' : 3,
        }
        #model.allocate(card_count)
        cards = [
            ['CORD1R', 1, 1, 2, 3],  # fails on self.k
            ['GRID', 1, 0, 0., 0., 0.],
            ['GRID', 2, 0, 1., 0., 0.],
            ['GRID', 3, 0, 1., 1., 0.],
        ]
        for card in cards:
            model.add_card(card, card[0], comment='comment\n', is_list=True)
        model.build()
        unused_c1 = model.coords.slice_by_coord_id(1)

    def test_cord2r_bad_01(self):
        model = BDF(debug=True)

        card_count = {
            'GRID' : 4,
            'CORD2R' : 3,
        }
        cards = {'GRID' : [], 'CORD2R' : []}
        grids = [
            ['GRID', 1, 0],
            ['GRID', 20, 0],
            ['GRID', 30, 0],
            ['GRID', 11, 5],
        ]
        for grid in grids:
            cards['GRID'].append(('', grid))

        coord_cards = [
            [
                'CORD2R', 1, 0, 0., 0., 0.,
                0., 0., 0.,
                0., 0., 0.],  # fails on self.k
            [
                'CORD2R', 2, 0, 0., 0., 0.,
                1., 0., 0.,
                0., 0., 0.],  # fails on normalize self.j
            [
                'CORD2R', 3, 0, 0., 0., 0.,
                1., 0., 0.,
                1., 1., 0.],  # passes
            [
                'CORD2R', 4, 0, 0., 1., 0.,
                1., 0., 0.,
                1., 1., 0.],  # passes
            [
                'CORD2R', 5, 4, 0., 1., 0.,
                1., 0., 0.,
                1., 1., 0.],  # passes
        ]
        for card in coord_cards:
            card_name = card[0]
            cid = card[1]
            if cid in [1, 2]:
                with self.assertRaises(RuntimeError):
                    cards[card_name].append(('', card))
            else:
                cards[card_name].append(('', card))
        model.add_cards(cards, card_count)
        model.build()

        # this runs because it's got rid=0
        coords = model.coords #.slice_by_coord_id(4)
        coords.transform_node_id_to_global_xyz(30) # nid
        coords.transform_node_id_to_global_xyz([30, 1, 11]) # [nid, nid, nid]
        coords.transform_node_id_to_local_by_coord_id([30, 1, 11], 4)  # [nid, ...], cp_goal
        coords.transform_node_id_to_local_by_coord_id([30, 1, 11], 0)  # [nid, ...], cp_goal

        xyz = [0., 0., 0.]
        coords.transform_xyz_to_global_by_coord_id(xyz, 4) # [xyz], cp_initial
        xyz = [
            [0., 0., 0.],
            [1., 1., 1.],
        ]
        coords.transform_xyz_to_global_by_coord_id(xyz, 4) # [xyz], cp_initial

        # global from global
        coords.transform_xyz_to_global_by_coord_id(xyz, 0) # [xyz], cp_initial

        # this doesn't run because rid != 0
        with self.assertRaises(RuntimeError):
            # cp=0 doesn't exist
            coords.transform_xyz_to_global_by_coord_id(xyz, 2)

    def test_cord2_rcs_01(self):
        """
        all points are located at <30,40,50>
        """
        model = BDF(debug=False)
        card_count = {
            'GRID' : 3,
            'CORD2R' : 1,
            'CORD2C' : 1,
            'CORD2S' : 1,
        }
        cards = []
        card_lines = [
            [
                #'$ Femap with NX Nastran Coordinate System 10 : rectangular defined in a rectangular',
                'CORD2R*               10               0             10.              5.',
                '*                     3.   10.3420201433   4.53015368961   3.81379768136*       ',
                '*          10.7198463104   5.68767171433   3.09449287122',],
            [
                #'$ Femap with NX Nastran Coordinate System 11 : cylindrical defined in rectangular',
                'CORD2C*               11               0              7.              3.',
                '*                     9.   7.64278760969   2.73799736977   9.71984631039*       ',
                '*          7.75440650673   3.37968226211   8.46454486422',],
            [
                #'$ Femap with NX Nastran Coordinate System 12 : spherical defined in rectangular',
                'CORD2S*               12               0             12.              8.',
                '*                     5.   12.6427876097   7.86697777844   5.75440650673*       ',
                '*          12.6634139482   8.58906867688   4.53861076379',],

            [
                'GRID*                 10              10   42.9066011565   34.2422137135',
                '*          28.6442730262               0',],
            [
                'GRID*                 11              11   48.8014631871   78.8394787869',
                '*          34.6037164304               0',],
            [
                'GRID*                 12              12   58.0775343829   44.7276544324',
                '*          75.7955331161               0',],
        ]
        cards, card_count = model.add_cards_lines(card_lines)
        model.allocate(card_count, cards)
        model.build()

        for nid in model.nodes:
            a = array([30., 40., 50.])
            b = model.Node(nid).get_position()
            self.assertTrue(allclose(array([30., 40., 50.]), model.Node(nid).get_position()), str(a - b))

    def test_cord2_rcs_02(self):
        """
        all points are located at <30,40,50>
        """
        model = BDF(debug=False)
        #card_count = {
            #'GRID' : 3,
            #'CORD2R' : 1,
            #'CORD2C' : 2,
            #'CORD2S' : 1,
        #}
        #model.allocate(card_count)
        card_lines = [
            [
                'CORD2C*                1               0              0.              0.',
                '*                     0.              0.              0.              1.*       ',
                '*                     1.              0.              1.',],
            [
                #'$ Femap with NX Nastran Coordinate System 20 : rectangular defined in cylindrical',
                'CORD2R*               20               1              7.             20.',
                '*                    -6.   7.07106781187   28.1301023542             -6.*       ',
                '*          7.70710678119             20.  -5.29289321881',],
            [
                #'$ Femap with NX Nastran Coordinate System 21 : cylindrical defined in cylindrical',
                'CORD2C*               21               1             15.            -30.',
                '*                    12.   14.6565766735  -30.3177805524   12.9355733712*       ',
                '*          14.6234241583  -26.4257323272   11.9304419665',],
            [
                #'$ Femap with NX Nastran Coordinate System 22 : spherical defined in cylindrical',
                'CORD2S*               22               1              5.            -75.',
                '*                    20.   5.66032384035  -82.9319986389   19.8502545865*       ',
                '*          4.88876051026  -73.8006653677   19.0116094889',],
            [
                'GRID*                 20              20   64.2559135157  -14.9400459772',
                '*          27.3271005317               0',],
            [
                'GRID*                 21              21   52.8328862418  -28.8729017195',
                '*           34.615939507               0',],
            [
                'GRID*                 22              22   61.1042111232   158.773483595',
                '*           -167.4951724               0',],
        ]
        cards, card_count = model.add_cards_lines(card_lines)
        model.allocate(card_count, cards)
        #for lines in cards:
            #card = model.add(lines)
            #model.add_card(card, card[0])
        model.build()
        for nid in model.nodes:
            a = array([30., 40., 50.])
            b = model.Node(nid).get_position()
            self.assertTrue(allclose(array([30., 40., 50.]), model.Node(nid).get_position()), str(a - b))

    def test_cord2_rcs_03(self):
        """
        all points are located at <30,40,50>
        """
        model = BDF(debug=False)
        #card_count = {
            #'GRID' : 3,
            #'CORD2R' : 1,
            #'CORD2C' : 1,
            #'CORD2S' : 2,
        #}

        card_lines = [
            [
                'CORD2S*                2               0              0.              0.',
                '*                     0.              0.              0.              1.*       ',
                '*                     1.              0.              1.',],
            [
                #'$ Femap with NX Nastran Coordinate System 30 : rectangular in spherical',
                'CORD2R*               30               2             14.             30.',
                '*                    70.    13.431863852   32.1458443949   75.2107442927*       ',
                '*          14.4583462334   33.4569982885   68.2297989286',],
            [
                #'$ Femap with NX Nastran Coordinate System 31 : cylindrical in spherical',
                'CORD2C*               31               2              3.             42.',
                '*                  -173.   2.86526881213   45.5425615252   159.180363517*       ',
                '*          3.65222385965   29.2536614627  -178.631312271',],
            [
                #'$ Femap with NX Nastran Coordinate System 32 : spherical in spherical',
                'CORD2S*               32               2             22.             14.',
                '*                    85.   22.1243073983   11.9537753718   77.9978191005*       ',
                '*          21.0997242967   13.1806120497   88.4824763008',],
            [
                'GRID*                 30              30   40.7437952957  -23.6254877994',
                '*           -33.09784854               0',],
            [
                'GRID*                 31              31   62.9378078196   15.9774797923',
                '*          31.0484428362               0',],
            [
                'GRID*                 32              32   53.8270847449   95.8215692632',
                '*          159.097767463               0',],
        ]
        cards, card_count = model.add_cards_lines(card_lines)
        model.allocate(card_count, cards)
        model.build()

        for nid in model.nodes:
            a = array([30., 40., 50.])
            b = model.Node(nid).get_position()
            self.assertTrue(allclose(array([30., 40., 50.]), model.Node(nid).get_position()), str(a - b))

    def test_cord1c_01(self):
        lines = ['cord1c,2,1,4,3']
        grids = [
            ['GRID', 4, 0, 0.0, 0., 0.],
            ['GRID', 3, 0, 0.0, 0., 1.],
            ['GRID', 1, 0, 0.0, 1., 0.],
        ]

        card_count = {
            'CORD1C' : 1,
            'GRID' : 3,
        }

        model = BDF(debug=False)
        model.allocate(card_count)
        model.add_card(lines, 'CORD1C', is_list=False)
        for grid in grids:
            model.add_card(grid, grid[0], is_list=True)
        model.build()

        size = 8
        bdf_file = StringIO()
        card = model.coords.slice_by_coord_id(2)
        self.assertEqual(card.get_cid_by_coord_id(), 2)
        self.assertEqual(card.get_rid_by_coord_id(), 0)
        card.write_card(bdf_file, size=8, is_double=False)

    def test_cord1s_01(self):
        cord1s = ['cord1s,2, 1,4,3']
        grids = [
            ['GRID', 4, 0, 0.0, 0., 0.],
            ['GRID', 3, 0, 0.0, 0., 1.],
            ['GRID', 1, 0, 0.0, 1., 0.],
        ]
        card_count = {
            'CORD1S' : 1,
            'GRID' : 3,
        }
        model = BDF(debug=False)
        cards = {
            'GRID' : [
                ('', grids[0]),
                ('', grids[1]),
                ('', grids[2]),
            ],
            'CORD1S' : [('', cord1s)]
        }
        model.add_cards(cards, card_count)
        model.build()

        size = 8
        bdf_file = StringIO()
        card = model.coords.slice_by_coord_id(2)
        self.assertEqual(card.get_cid_by_coord_id(), 2)
        self.assertEqual(card.get_rid_by_coord_id(), 0)
        card.write_card(bdf_file, size=8, is_double=False)
        #card.raw_fields()

    def test_cord2r_02(self):
        grid = ['GRID       20143       7 -9.31-4  .11841 .028296']
        coord = [
            'CORD2R         7           1.135 .089237  -.0676    .135 .089237  -.0676',
            '           1.135 .089237   .9324']

        model = BDF(debug=False)
        card_count = {
            'GRID' : 1,
            'CORD2R' : 1,
        }
        cards = {
            'CORD2R' : [('', coord)],
            'GRID' : [('', grid)],
        }
        model.add_cards(cards, card_count)
        model.build()

        g = model.grid.slice_by_node_id(20143)
        #xyz = g.get_position()
        xyz = model.coords.get_global_position_by_node_id(20143, g.cp[0])[0]

        # by running it through Patran...
        #GRID     20143          1.1067  .207647 -.068531
        expected = array([1.106704, .207647, -0.068531])
        diff = xyz - expected

        msg = '\nexpected=%s \nactual  =%s \ndiff    =%s' % (expected, xyz, diff)
        assert allclose(diff, 0.), msg
        coord = model.coords.slice_by_coord_id(7)
        T = coord.T[0, :, :]
        #self.assertTrue(array_equal(T, coord.beta_n(2)))

    def _get_nodes(self, grids, grids_expected, coords):
        model = BDF(debug=False)

        card_count = {
            'GRID' : len(grids),
            'CORD2R' : len(coords),
        }
        cards = {'GRID' : [], 'CORD2R' : []}
        for grid in grids:
            nid, cid, x, y, z = grid
            card = ['GRID', nid, cid, x, y, z]
            cards['GRID'].append(('', card))

        for coord in coords:
            cid, rid, x, y, z = coord
            card = ['CORD2R', cid, rid] + x + y + z
            cards['CORD2R'].append(('', card))
            #coordObj = model.coords.slice_by_coord_id(cid)
        model.add_cards(cards, card_count)
        model.build()

        for (i, grid) in enumerate(grids_expected):
            nid, cid, x, y, z = grid
            nodes = model.grid
            pos = nodes.get_position_by_node_id([nid])[0]
            n = array([x, y, z])
            msg = 'i=%s expected=%s actual=%s\n' % (i, n, pos)
            #print(msg)
            assert allclose(n, pos), msg


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
