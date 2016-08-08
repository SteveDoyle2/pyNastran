from __future__ import print_function
import unittest
from six import iteritems
from numpy import array, allclose, array_equal, cross

from pyNastran.bdf.bdf import BDF, BDFCard, CORD1R, CORD1C, CORD1S, CORD2R, CORD2C, CORD2S
from pyNastran.bdf.utils import TransformLoadWRT

bdf = BDF(debug=False)  # don't load this up with stuff
class TestCoords(unittest.TestCase):
    def test_same(self):
        grids = [
            [1, 0, 0., 0., 1.],
            [2, 0, 0., 1., 0.],
            [3, 0, 1., 0., 0.],
            [4, 0, 1., 1., 1.],
            [5, 0, 1., 1., 0.],
        ]
        grids_expected = grids
        coords = []
        self.get_nodes(grids, grids_expected, coords)

    def test_shift(self):
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

        coords = [  # cid,rid, origin,      zaxis,     xaxis
            [1, 0, [1., 1., 1.], [1., 1., 2.], [2., 1., 1.]],
        ]
        self.get_nodes(grids, grids_expected, coords)

    def test_rotate(self):
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
        self.get_nodes(grids, grids_expected, coords)

    def test_rotate2(self):
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
        self.get_nodes(grids, grids_expected, coords)

    def test_rotate3(self):
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
        self.get_nodes(grids, grids_expected, coords)

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
        self.get_nodes(grids, grids_expected, coords)

    def test_cord1r_01(self):
        lines = ['cord1r,2,1,4,3']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        coord = CORD1R.add_card(card)
        self.assertEqual(coord.Cid(), 2)
        self.assertEqual(coord.Rid(), 0)
        coord.write_card(size, 'dummy')
        coord.raw_fields()

    def test_cord2c_01(self):
        lines = [
            'CORD2C*                3               0              0.              0.',
            '*                     0.              0.              0.              1.*',
            '*                     1.              0.              1.'
        ]
        model = BDF(debug=False)
        card = model.process_card(lines)
        cardi = BDFCard(card)
        coord = CORD2C.add_card(cardi)
        model.add_coord(coord)

        lines = [
            'CORD2R         4       3     10.      0.      5.     10.     90.      5.',
            '             10.      0.      6.'
        ]
        card = model.process_card(lines)
        cardi = BDFCard(card)
        coord = CORD2R.add_card(cardi)
        model.add_coord(coord)
        model.cross_reference()

        cord2r = model.Coord(3)
        self.assertEqual(cord2r.Cid(), 3)
        self.assertEqual(cord2r.Rid(), 0)

        cord2r = model.Coord(4)
        self.assertEqual(cord2r.Cid(), 4)
        self.assertEqual(cord2r.Rid(), 3)

        self.assertTrue(allclose(cord2r.i, array([0., 0., 1.])))
        delta = cord2r.j - array([1., 1., 0.]) / 2**0.5
        self.assertTrue(allclose(cord2r.j, array([1., 1., 0.]) / 2**0.5), str(delta))
        delta = cord2r.k - array([-1., 1., 0.]) / 2**0.5
        self.assertTrue(allclose(cord2r.k, array([-1., 1., 0.]) / 2**0.5), str(delta))

    def test_grid_01(self):
        model = BDF(debug=False)
        cards = [
            #['CORD1R', 1, 1, 2, 3],  # fails on self.k
            ['GRID', 1, 0, 0., 0., 0.],
            ['GRID', 2, 0, 1., 0., 0.],
            ['GRID', 4, 0, 1., 2., 3.],
        ]
        for card in cards:
            model.add_card(card, card[0], comment='comment', is_list=True)

        #+------+-----+----+----+----+----+----+----+------+
        #|   0  |  1  | 2  | 3  | 4  | 5  |  6 | 7  |  8   |
        #+======+=====+====+====+====+====+====+====+======+
        #| GRID | NID | CP | X1 | X2 | X3 | CD | PS | SEID |
        #+------+-----+----+----+----+----+----+----+------+
        node = model.Node(4)
        self.assertEqual(node.get_field(1), 4)
        self.assertEqual(node.get_field(2), 0)
        self.assertEqual(node.get_field(3), 1.)
        self.assertEqual(node.get_field(4), 2.)
        self.assertEqual(node.get_field(5), 3.)

        node.update_field(1, 5)
        node.update_field(2, 6)
        node.update_field(3, 7.)
        node.update_field(4, 8.)
        node.update_field(5, 9.)
        with self.assertRaises(KeyError):
            node.update_field(9, 'dummy')

        self.assertEqual(node.get_field(1), 5)
        self.assertEqual(node.get_field(2), 6)
        self.assertEqual(node.get_field(3), 7.)
        self.assertEqual(node.get_field(4), 8.)
        self.assertEqual(node.get_field(5), 9.)
        with self.assertRaises(KeyError):
            node.get_field(9)

    def test_cord1_01(self):
        model = BDF(debug=False)
        cards = [
            ['CORD1R', 1, 1, 2, 3],  # fails on self.k
            ['GRID', 1, 0, 0., 0., 0.],
            ['GRID', 2, 0, 1., 0., 0.],
            ['GRID', 3, 0, 1., 1., 0.],
        ]
        for card in cards:
            model.add_card(card, card[0], comment='comment', is_list=True)
        c1 = model.Coord(1)
        self.assertEqual(c1.G1(), 1)
        self.assertEqual(c1.G2(), 2)
        self.assertEqual(c1.G3(), 3)

        model.cross_reference()
        self.assertEqual(c1.G1(), 1)
        self.assertEqual(c1.G2(), 2)
        self.assertEqual(c1.G3(), 3)

        self.assertEqual(c1.node_ids, [1, 2, 3])

    def test_cord2_bad_01(self):
        model = BDF(debug=False)
        cards = [
            ['CORD2R', 1, 0, 0., 0., 0.,
             0., 0., 0.,
             0., 0., 0.],  # fails on self.k
            ['CORD2R', 2, 0, 0., 0., 0.,
             1., 0., 0.,
             0., 0., 0.],  # fails on normalize self.j
            ['CORD2R', 3, 0, 0., 0., 0.,
             1., 0., 0.,
             1., 1., 0.],  # passes
            ['CORD2R', 4, 0, 0., 1., 0.,
             1., 0., 0.,
             1., 1., 0.],  # passes
            ['CORD2R', 5, 4, 0., 1., 0.,
             1., 0., 0.,
             1., 1., 0.],  # passes
        ]
        for card in cards:
            cid = card[1]
            if cid in [1, 2]:
                with self.assertRaises(RuntimeError):
                    model.add_card(card, card[0], is_list=True)
            else:
                model.add_card(card, card[0], is_list=True)

        # this runs because it's got rid=0
        cord4 = model.Coord(4)
        cord4.transform_node_to_global([0., 0., 0.])

        # this doesn't run because rid != 0
        cord5 = model.Coord(5)
        with self.assertRaises(RuntimeError):
            cord5.transform_node_to_global([0., 0., 0.])
        model.cross_reference()

    def test_cord2_rcs_01(self):
        """
        all points are located at <30,40,50>
        """
        model = BDF(debug=False)
        cards = [
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
        for lines in cards:
            card = model.process_card(lines)
            model.add_card(card, card[0])
        model.cross_reference()
        for nid in model.nodes:
            a = array([30., 40., 50.])
            b = model.Node(nid).get_position()
            self.assertTrue(allclose(array([30., 40., 50.]), model.Node(nid).get_position()), str(a - b))

    def test_cord2_rcs_02(self):
        """
        all points are located at <30,40,50>
        """
        model = BDF(debug=False)
        cards = [
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
        for lines in cards:
            card = model.process_card(lines)
            model.add_card(card, card[0])
        model.cross_reference()
        for nid in model.nodes:
            a = array([30., 40., 50.])
            b = model.Node(nid).get_position()
            self.assertTrue(allclose(array([30., 40., 50.]), model.Node(nid).get_position()), str(a - b))

    def test_cord2_rcs_03(self):
        """
        all points are located at <30,40,50>
        """
        model = BDF(debug=False)
        cards = [
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
        for lines in cards:
            card = model.process_card(lines)
            model.add_card(card, card[0])
        model.cross_reference()
        for nid in model.nodes:
            a = array([30., 40., 50.])
            b = model.Node(nid).get_position()
            self.assertTrue(allclose(array([30., 40., 50.]), model.Node(nid).get_position()), str(a - b))

    def test_cord1c_01(self):
        lines = ['cord1c,2,1,4,3']
        card = bdf.process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = CORD1C.add_card(cardi)
        self.assertEqual(card.Cid(), 2)
        self.assertEqual(card.Rid(), 0)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_cord1s_01(self):
        lines = ['cord1s,2,1,4,3']
        card = bdf.process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = CORD1S.add_card(cardi)
        self.assertEqual(card.Cid(), 2)
        self.assertEqual(card.Rid(), 0)
        card.write_card(size, 'dummy')
        card.raw_fields()

    def test_cord2r_02(self):
        grid = ['GRID       20143       7 -9.31-4  .11841 .028296']
        coord = [
            'CORD2R         7           1.135 .089237  -.0676    .135 .089237  -.0676',
            '           1.135 .089237   .9324'
        ]

        model = BDF(debug=False)
        card = model.process_card(grid)
        model.add_card(card, card[0])

        card = model.process_card(coord)
        model.add_card(card, card[0])
        model.cross_reference()
        coord = model.Coord(7)
        #print(coord.origin)
        #print(coord.i, coord.j, coord.k)

        g = model.Node(20143)
        #print(g.Position(debug=False))
        xyz = g.get_position()

        # by running it through Patran...
        #GRID     20143          1.1067  .207647 -.068531
        expected = array([1.106704, .207647, -0.068531])
        diff = xyz - expected

        msg = '\nexpected=%s \nactual  =%s \ndiff    =%s' % (expected, xyz, diff)
        assert allclose(diff, 0.), msg
        coord = model.Coord(7)
        coord.beta_n(1)
        coord.beta_n(2)
        coord.beta_n(3)
        coord.beta_n(6)
        with self.assertRaises(AttributeError):
            self.assertTrue(array_equal(coord.T(), coord.beta_n(2)))
        #with self.assertRaises(NotImplementedError):
            #self.assertTrue(array_equal(coord.T(), coord.beta_n(2)))

    def get_nodes(self, grids, grids_expected, coords):
        model = BDF(debug=False)

        for grid in grids:
            (nid, cid, x, y, z) = grid
            model.add_card(['GRID', nid, cid, x, y, z], 'GRID')

        for coord in coords:
            (cid, rid, x, y, z) = coord
            model.add_card(['CORD2R', cid, rid] + x + y + z, 'CORD2R')
            #coord_obj = model.Coord(cid)

        model.cross_reference()

        for (i, grid) in enumerate(grids_expected):
            (nid, cid, x, y, z) = grid
            node = model.Node(nid)
            pos = node.get_position()
            n = array([x, y, z])

            msg = 'i=%s expected=%s actual=%s\n' % (i, n, pos)
            msg += 'n=%s grid=\n%s' % (nid, node)
            coord = node.cp
            msg += 'n=%s coord=\n%s' % (node.nid, coord)
            while coord.rid:
                msg += 'n=%s rcoord=\n%s' % (node.nid, coord.rid)
                coord = coord.rid
            assert allclose(n, pos), msg


    def test_A(self):
        origin = array([0., 0., 0.])
        zaxis = array([0., 0., 1.])
        xzplane = array([1., 0., 0.])
        cid0 = CORD2R(cid=0, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane)
        Lx = 2.
        Ly = 0.
        Lz = 3.
        Fy = 1.
        origin = array([-Lx, 0., -Lz])
        z_axis = origin + array([0., 0., 1.])
        xz_plane = origin + array([1., 0., 1.])
        rid = 0
        data = [1, rid] + list(origin) + list(z_axis) + list(xz_plane)

        Fxyz = [0., -Fy, 0.]
        Mxyz = [0., 0., 0.]
        cid_new = CORD2R.add_op2_data(data=data)
        model = None

        Fxyz_local, Mxyz_local = TransformLoadWRT(Fxyz, Mxyz, cid0, cid_new,
                                                  model, is_cid_int=False)

        r = array([Lx, Ly, Lz])
        F = array([0., -Fy, 0.])
        M = cross(r, F)
        self.assertTrue(array_equal(Fxyz_local, F)), "expected=%s actual=%s" % (F, Fxyz_local)
        self.assertTrue(array_equal(Mxyz_local, cross(r, F))), "expected=%s actual=%s" % (M, Mxyz_local)

    def test_B(self):
        origin = array([0., 0., 0.])
        zaxis = array([0., 0., 1.])
        xzplane = array([1., 0., 0.])
        cid0 = CORD2R(cid=0, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane)

        Lx = 2.
        Ly = 3.
        Lz = 5.
        Fy = 1.5
        origin = array([-Lx, -Ly, -Lz])
        z_axis = origin + array([0., 0., 1.])
        xz_plane = origin + array([1., 0., 1.])
        rid = 0
        data = [1, rid] + list(origin) + list(z_axis) + list(xz_plane)

        Fxyz = [0., -Fy, 0.]
        Mxyz = [0., 0., 0.]
        cid_new = CORD2R.add_op2_data(data=data)
        model = None

        Fxyz_local, Mxyz_local = TransformLoadWRT(Fxyz, Mxyz, cid0, cid_new,
                                                  model, is_cid_int=False)
        r = array([Lx, Ly, Lz])
        F = array([0., -Fy, 0.])
        M = cross(r, F)
        self.assertTrue(array_equal(Fxyz_local, F)), "expected=%s actual=%s" % (F, Fxyz_local)
        self.assertTrue(array_equal(Mxyz_local, cross(r, F))), "expected=%s actual=%s" % (M, Mxyz_local)

    def test_coord_adding(self):
        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [1., 0., 0.]
        cid1 =CORD2R(cid=1, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane,
                     comment='')

        xaxis = [1., 0., 0.]
        yaxis = [0., 1., 0.]
        zaxis = [0., 0., 1.]
        xz_plane = [1., 0., 1.]
        yz_plane = [0., 1., 1.]
        xy_plane = [1., 1., 0.]
        # x-axis
        cid2 = CORD2R.add_axes(cid=2, rid=0, origin=origin, xaxis=xaxis, yaxis=None, zaxis=None,
                               xyplane=None, yzplane=None, xzplane=xz_plane)

        cid3 = CORD2R.add_axes(cid=2, rid=0, origin=origin, xaxis=xaxis, yaxis=None, zaxis=None,
                               xyplane=xy_plane, yzplane=None, xzplane=None)

        # y-axis
        cid4 = CORD2R.add_axes(cid=4, rid=0, origin=origin, xaxis=None, yaxis=yaxis, zaxis=None,
                               xyplane=xy_plane, yzplane=None, xzplane=None)

        cid5 = CORD2R.add_axes(cid=5, rid=0, origin=origin, xaxis=None, yaxis=yaxis, zaxis=None,
                               xyplane=None, yzplane=yz_plane, xzplane=None)

        # z-axis
        cid4 = CORD2R.add_axes(cid=4, rid=0, origin=origin, xaxis=None, yaxis=None, zaxis=zaxis,
                               xyplane=None, yzplane=None, xzplane=xz_plane)

        cid5 = CORD2R.add_axes(cid=5, rid=0, origin=origin, xaxis=None, yaxis=None, zaxis=zaxis,
                               xyplane=None, yzplane=yz_plane, xzplane=None)

        # ijk
        cid6 = CORD2R.add_ijk(cid=6, rid=0, origin=origin, i=xaxis, j=yaxis, k=None)
        cid7 = CORD2R.add_ijk(cid=7, rid=0, origin=origin, i=xaxis, j=None, k=zaxis)
        cid8 = CORD2R.add_ijk(cid=8, rid=0, origin=origin, i=None, j=yaxis, k=zaxis)
        #cid6.add_ijk(rid=0, origin=origin, i=None, j=None, k=None)

    def test_cord1_referencing_01(self):
        bulk_data_lines = [
            'CORD2C       300        253.345 171.174 197.495 242.9242270.3229205.2989+',
            '+       254.1607163.4128297.19',
            'CORD2R       301     3000.0     0.0     0.0     100.0   90.444676.406-13+',
            '+       100.0   -179.5557.514-13',
            'CORD2C       323        222.919 185.412 198.925 233.343986.26359191.1203+',
            '+       309.1199190.5056249.3577',
            'CORD2R        40        250.0   203.5575262.4969348.3264221.7756262.3345+',
            '+       250.0005202.6639162.5009',
            'CORD2R       111        253.345 171.174 197.495 263.764272.02488189.6917+',
            '+       153.8893160.7869196.6775',
            'CORD2C       112     1110.0     0.0     0.0     -7.88-14-2.58-14100.0   +',
            '+       100.0   -7.33-148.148-14',
            'CORD2C       114     1120.0     0.0     0.0     1.305-150.0     100.0   +',
            '+       100.0   -90.445 7.994-15',
            'CORD2R       116     1140.0     0.0     0.0     100.0   -90.0   -8.0-14 +',
            '+       100.0   180.0   9.859-14',
            'CORD2R       201        251.1275219.2556219.8197251.8979319.249 220.6651+',
            '+       253.9912220.0786119.8641',
            'CORD2R       202        251.4405218.3364184.3664256.3151317.7323174.5327+',
            '+       251.9847228.1553283.8817',
            'CORD2R       207        270.2522201.8565210.846 370.2484201.6459211.6994+',
            '+       271.1036200.9139110.8541',
            'CORD2R       208        230.5389215.4915162.7607234.0098310.4138131.4924+',
            '+       228.789 184.266867.77676',
            'CORD2R       511      400.0     0.0     0.0     -100.0  -1.7-13 -9.36-13+',
            '+       -9.19-131.074-12100.0',
            'CORD2R       321     3230.0     0.0     0.0     100.0   -3.26-141.363-13+',
            '+       100.0   90.0    1.11-13',
            'CORD2R        98        250.0   203.6961278.0626249.9841204.6729378.0579+',
            '+       315.2573127.9272278.8131',
            'CORD1R       932   23315   23310   22155',
            'GRID       22155     321-2.79-145.396   1.388-16     321       0',
            'GRID       23310        256.9914187.4238218.859      932       0',
            'GRID       23315     3210.0     0.0     0.0          932       0',
        ]
        model = BDF(debug=False)
        #model.echo = True
        #cards, card_count = model.get_bdf_cards(bulk_data_lines)
        cards, card_count = model.get_bdf_cards_dict(bulk_data_lines)
        model._parse_cards(cards, card_count)

        #print(model.card_count)
        assert model.card_count['CORD1R'] == 1, model.card_count
        assert model.card_count['CORD2C'] == 4, model.card_count
        assert model.card_count['CORD2R'] == 11, model.card_count
        assert model.card_count['GRID'] == 3, model.card_count
        model.cross_reference()
        for cid, coord in sorted(iteritems(model.coords)):
            assert coord.i is not None, coord

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
