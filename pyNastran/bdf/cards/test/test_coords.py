"""
defines:
 - TestCoords

"""
# pylint: disable=R0201,C0103
from io import StringIO
from copy import deepcopy
import unittest
import numpy as np
from numpy import array, allclose, array_equal, cross

from pyNastran.bdf.cards.coordinate_systems import (
    create_coords_along_line, get_nodes_along_axis_in_coords,
    define_coord_e123,
    CORD1R, CORD1C, CORD1S,
    CORD2R, CORD2C, #CORD2S,
    CORD3G)
from pyNastran.bdf.bdf import BDF, BDFCard
from pyNastran.bdf.utils import Position, PositionWRT, TransformLoadWRT
from pyNastran.bdf.cards.aero.utils import make_monpnt1s_from_cids
from pyNastran.bdf.cards.test.utils import save_load_deck

class TestCoords(unittest.TestCase):
    """tests the coordinate systems and their transforms"""
    def test_same(self):
        """simple coordinate equality test"""
        grids = [
            [1, 0, 0., 0., 1.],
            [2, 0, 0., 1., 0.],
            [3, 0, 1., 0., 0.],
            [4, 0, 1., 1., 1.],
            [5, 0, 1., 1., 0.],
        ]
        grids_expected = grids
        coords = []
        get_nodes(grids, grids_expected, coords)

    def test_shift(self):
        """simple coordinate test of origin shifting"""
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
        get_nodes(grids, grids_expected, coords)

    def test_rotate(self):
        """simple coordinate test of 90 degree rotations"""
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
        get_nodes(grids, grids_expected, coords)

    def test_rotate2(self):
        """simple coordinate test of 90 degree rotations"""
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
        get_nodes(grids, grids_expected, coords)

    def test_rotate3(self):
        """simple coordinate test of 90 degree rotations"""
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
        get_nodes(grids, grids_expected, coords)

    def test_rid_1(self):
        """simple coordinate test of a referenced coordinate system"""
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
        get_nodes(grids, grids_expected, coords)

    def test_cord1r_01(self):
        """simple CORD1R input/output test"""
        lines = ['cord1r,2,1,4,3']
        model = BDF(debug=False)
        card = model._process_card(lines)
        card = BDFCard(card)

        size = 8
        coord = CORD1R.add_card(card)
        self.assertEqual(coord.Cid(), 2)
        self.assertEqual(coord.Rid(), 0)
        coord.write_card(size, 'dummy')
        coord.raw_fields()
        make_tri(model)
        save_load_deck(model, run_renumber=False)

    def test_cord2c_01(self):
        """simple CORD2R/CORD2C input/output test"""
        lines = [
            'CORD2C*                3               0              0.              0.',
            '*                     0.              0.              0.              1.*',
            '*                     1.              0.              1.'
        ]
        model = BDF(debug=False)
        card = model._process_card(lines)
        cardi = BDFCard(card)
        cord2c = CORD2C.add_card(cardi)
        model._add_coord_object(cord2c)

        lines = [
            'CORD2R         4       3     10.      0.      5.     10.     90.      5.',
            '             10.      0.      6.'
        ]
        card = model._process_card(lines)
        cardi = BDFCard(card)
        cord2r = CORD2R.add_card(cardi)
        model._add_coord_object(cord2r)
        model.cross_reference()

        cord2r_b = model.Coord(3)
        self.assertEqual(cord2r_b.Cid(), 3)
        self.assertEqual(cord2r_b.Rid(), 0)

        cord2r_c = model.Coord(4)
        self.assertEqual(cord2r_c.Cid(), 4)
        self.assertEqual(cord2r_c.Rid(), 3)

        self.assertTrue(allclose(cord2r_c.i, array([0., 0., 1.])))
        delta = cord2r_c.j - array([1., 1., 0.]) / 2**0.5
        self.assertTrue(allclose(cord2r_c.j, array([1., 1., 0.]) / 2**0.5), str(delta))
        delta = cord2r_c.k - array([-1., 1., 0.]) / 2**0.5
        self.assertTrue(allclose(cord2r_c.k, array([-1., 1., 0.]) / 2**0.5), str(delta))


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
                #'$ Coordinate System 10 : rectangular defined in a rectangular',
                'CORD2R*               10               0             10.              5.',
                '*                     3.   10.3420201433   4.53015368961   3.81379768136*       ',
                '*          10.7198463104   5.68767171433   3.09449287122',],
            [
                #'$ Coordinate System 11 : cylindrical defined in rectangular',
                'CORD2C*               11               0              7.              3.',
                '*                     9.   7.64278760969   2.73799736977   9.71984631039*       ',
                '*          7.75440650673   3.37968226211   8.46454486422',],
            [
                #'$ Coordinate System 12 : spherical defined in rectangular',
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
            card = model._process_card(lines)
            model.add_card(card, card[0])

        unused_xyz_cid0b = model.get_xyz_in_coord_no_xref(cid=0, fdtype='float64')
        unused_xyz_cid0c = model.get_xyz_in_coord_no_xref(cid=12, fdtype='float64')
        model.cross_reference()

        xyz_cid0_actual = array([
            [30., 40., 50.],
            [30., 40., 50.],
            [30., 40., 50.],
        ], dtype='float64')
        for nid in model.nodes:
            node = model.Node(nid)
            a = array([30., 40., 50.])
            b = node.get_position()
            self.assertTrue(allclose(array([30., 40., 50.]),
                                     node.get_position()), str(a - b))

        xyz_cid0 = model.get_xyz_in_coord(cid=0, fdtype='float64')
        array_equal(xyz_cid0_actual, xyz_cid0)

        unused_icd_transform, icp_transform, xyz_cp, nid_cp_cd = model.get_displacement_index_xyz_cp_cd()
        xyz_cid0_xform = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nid_cp_cd[:, 0], icp_transform, cid=0)
        array_equal(xyz_cid0_actual, xyz_cid0_xform)
        assert array_equal(nid_cp_cd[:, 0], array([10, 11, 12]))

        unused_xyz_cid_10 = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nid_cp_cd[:, 0], icp_transform, cid=10)
        unused_xyz_cid_11 = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nid_cp_cd[:, 0], icp_transform, cid=11)
        unused_xyz_cid_12 = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nid_cp_cd[:, 0], icp_transform, cid=12)

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
                #'$ Coordinate System 20 : rectangular defined in cylindrical',
                'CORD2R*               20               1              7.             20.',
                '*                    -6.   7.07106781187   28.1301023542             -6.*       ',
                '*          7.70710678119             20.  -5.29289321881',],
            [
                #'$ Coordinate System 21 : cylindrical defined in cylindrical',
                'CORD2C*               21               1             15.            -30.',
                '*                    12.   14.6565766735  -30.3177805524   12.9355733712*       ',
                '*          14.6234241583  -26.4257323272   11.9304419665',],
            [
                #'$ Coordinate System 22 : spherical defined in cylindrical',
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
            card = model._process_card(lines)
            model.add_card(card, card[0])
        unused_xyz_cid0b = model.get_xyz_in_coord_no_xref(cid=0, fdtype='float64')
        unused_xyz_cid0c = model.get_xyz_in_coord_no_xref(cid=22, fdtype='float64')
        model.cross_reference()

        xyz_cid0_actual = array([
            [30., 40., 50.],
            [30., 40., 50.],
            [30., 40., 50.],
        ], dtype='float64')
        for nid in model.nodes:
            a = array([30., 40., 50.])
            b = model.Node(nid).get_position()
            self.assertTrue(allclose(array([30., 40., 50.]),
                                     model.Node(nid).get_position()), str(a - b))
        xyz_cid0 = model.get_xyz_in_coord(cid=0, fdtype='float64')
        array_equal(xyz_cid0_actual, xyz_cid0)

        unused_icd_transform, icp_transform, xyz_cp, nid_cp_cd = model.get_displacement_index_xyz_cp_cd()
        xyz_cid0_xform = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nid_cp_cd[:, 0], icp_transform, cid=0)
        array_equal(xyz_cid0_actual, xyz_cid0_xform)
        assert array_equal(nid_cp_cd[:, 0], array([20, 21, 22]))

        unused_xyz_cid_20 = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nid_cp_cd[:, 0], icp_transform, cid=20)
        unused_xyz_cid_21 = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nid_cp_cd[:, 0], icp_transform, cid=21)
        unused_xyz_cid_22 = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nid_cp_cd[:, 0], icp_transform, cid=22)


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
                #'$ Coordinate System 30 : rectangular in spherical',
                'CORD2R*               30               2             14.             30.',
                '*                    70.    13.431863852   32.1458443949   75.2107442927*       ',
                '*          14.4583462334   33.4569982885   68.2297989286',],
            [
                #'$ Coordinate System 31 : cylindrical in spherical',
                'CORD2C*               31               2              3.             42.',
                '*                  -173.   2.86526881213   45.5425615252   159.180363517*       ',
                '*          3.65222385965   29.2536614627  -178.631312271',],
            [
                #'$ Coordinate System 32 : spherical in spherical',
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
            card = model._process_card(lines)
            model.add_card(card, card[0])
        unused_xyz_cid0b = model.get_xyz_in_coord_no_xref(cid=0, fdtype='float64')
        unused_xyz_cid0c = model.get_xyz_in_coord_no_xref(cid=32, fdtype='float64')
        bdf_file = StringIO()
        model.write_bdf(bdf_file, close=False)
        bdf_file.seek(0)

        #-------------------------------------------------
        model.cross_reference()

        #-------------------------------------------------

        xyz_cid0_actual = array([
            [30., 40., 50.],
            [30., 40., 50.],
            [30., 40., 50.],
        ], dtype='float64')
        for nid in model.nodes:
            node = model.Node(nid)
            a = array([30., 40., 50.])
            b = node.get_position()
            self.assertTrue(allclose(array([30., 40., 50.]),
                                     node.get_position()), str(a - b))
        xyz_cid0 = model.get_xyz_in_coord(cid=0, fdtype='float64')
        assert np.allclose(xyz_cid0_actual, xyz_cid0), '%s' % (xyz_cid0_actual - xyz_cid0)

        out = model.get_displacement_index_xyz_cp_cd()
        unused_icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
        nids = nid_cp_cd[:, 0]
        xyz_cid0_xform = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nids, icp_transform, cid=0)
        array_equal(xyz_cid0_actual, xyz_cid0_xform)
        assert array_equal(nids, array([30, 31, 32]))

        for cid in [30, 31, 32]:
            unused_xyz_cid_a = model.transform_xyzcp_to_xyz_cid(
                xyz_cp, nids, icp_transform, cid=cid)
            #assert np.allclose(xyz_cid_a, xyz_cid_b), '%s' % np.isclose(xyz_cid_a, xyz_cid_b)

            #print(xyz_cid_a)
            #print(xyz_cid_b)
            #print(xyz_cid_a - xyz_cid_b)
            #print('-------------')
            #assert array_equal(xyz_cid_a, xyz_cid_b), 'error=%s'  % (
                #xyz_cid_a - xyz_cid_b)

        #---------------------------------------------
        xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nids, icp_transform,
            cid=0, atol=None)
        array_equal(xyz_cid0_actual, xyz_cid0)

        model.write_bdf(bdf_file, close=False)

        model3 = BDF(debug=False)
        origin = [14., 30., 70.]
        zaxis = [13.431863852, 32.1458443949, 75.2107442927]
        xzplane = [14.4583462334, 33.4569982885, 68.2297989286]
        cord2r = model3.add_cord2r(30, origin, zaxis, xzplane, rid=2, comment='')

        origin = [3., 42., -173.]
        zaxis = [2.86526881213, 45.5425615252, 159.180363517]
        xzplane = [3.65222385965, 29.2536614627, -178.631312271]
        cord2c = model3.add_cord2c(31, origin, zaxis, xzplane, rid=2, comment='')

        origin = [22., 14., 85.]
        zaxis = [22.1243073983, 11.9537753718, 77.9978191005]
        xzplane = [21.0997242967, 13.1806120497, 88.4824763008]
        cord2s = model3.add_cord2s(32, origin, zaxis, xzplane, rid=2, comment='')

        assert cord2r == model.coords[cord2r.cid], 'cord2r:\n%r\ncord2r[cid]:\n%r' % (str(cord2r), str(model.coords[cord2r.cid]))
        assert cord2c == model.coords[cord2c.cid], 'cord2c:\n%r\ncord2c[cid]:\n%r' % (str(cord2c), str(model.coords[cord2c.cid]))
        assert cord2s == model.coords[cord2s.cid], 'cord2s:\n%r\ncord2s[cid]:\n%r' % (str(cord2s), str(model.coords[cord2s.cid]))

    def test_cord1c_01(self):
        lines = ['cord1c,2,1,4,3']
        model = BDF(debug=False)
        card = model._process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = CORD1C.add_card(cardi)
        self.assertEqual(card.Cid(), 2)
        self.assertEqual(card.Rid(), 0)
        card.write_card(size, 'dummy')
        card.raw_fields()

        model = BDF(debug=False)
        cid = 2
        grid1, grid2, grid3 = 1, 4, 3
        coord = model.add_cord1c(cid, grid1, grid2, grid3, comment='cord1c')
        coord.comment = ''
        make_tri(model)

        assert coord == card, 'card:\n%r\ncoord:\n%r' % (str(coord), str(card))
        model.cross_reference()
        save_load_deck(model, run_renumber=False)

    def test_cord1s_01(self):
        lines = ['cord1s,2,1,4,3']
        model = BDF(debug=False)
        card = model._process_card(lines)
        cardi = BDFCard(card)

        size = 8
        card = CORD1S.add_card(cardi)
        self.assertEqual(card.Cid(), 2)
        self.assertEqual(card.Rid(), 0)
        card.write_card(size, 'dummy')
        card.raw_fields()

        model = BDF(debug=False)
        model.set_error_storage(nparse_errors=0, stop_on_parsing_error=True,
                                nxref_errors=0, stop_on_xref_error=True)

        cid = 2
        grid1, grid2, grid3 = 1, 4, 3
        coord = model.add_cord1s(cid, grid1, grid2, grid3, comment='cord1c')
        coord.comment = ''
        assert coord == card, 'card:\n%r\ncoord:\n%r' % (str(coord), str(card))

        make_tri(model)
        coord.cross_reference(model)
        model2 = deepcopy(model)
        model2.cross_reference()
        save_load_deck(model2, run_renumber=False)
        unused_cord2s = coord.to_cord2x(model, rid=0)

        model.pop_parse_errors()
        model.pop_xref_errors()
        model.coords[cid] = coord
        model.cross_reference()
        save_load_deck(model, run_renumber=False)

    def test_cord2r_02(self):
        grid = ['GRID       20143       7 -9.31-4  .11841 .028296']
        coord = [
            'CORD2R         7           1.135 .089237  -.0676    .135 .089237  -.0676',
            '           1.135 .089237   .9324'
        ]

        model = BDF(debug=False)
        card = model._process_card(grid)
        model.add_card(card, card[0])

        card = model._process_card(coord)
        model.add_card(card, card[0])
        model.cross_reference()
        coord = model.Coord(7)
        #print(coord.origin)
        #print(coord.i, coord.j, coord.k)

        node = model.Node(20143)
        xyzp1 = Position(node.xyz, node.cp, model)
        xyzp2 = Position(node.xyz, node.cp_ref, model)
        xyz = node.get_position()
        assert np.array_equal(xyz, xyzp1)
        assert np.array_equal(xyz, xyzp2)

        xyz_same = PositionWRT([1., 2., 3.], 100, 100, model)
        assert np.array_equal(xyz_same, [1., 2., 3.])

        xyz_wrt_p1 = PositionWRT(node.xyz, node.cp, 0, model)
        xyz_wrt_p2 = PositionWRT(node.xyz, node.cp_ref, 0, model)
        xyz_wrt = node.get_position_wrt(model, 0)
        assert np.array_equal(xyz, xyz_wrt_p1)
        assert np.array_equal(xyz, xyz_wrt_p2)
        assert np.array_equal(xyz, xyz_wrt)

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

        model2 = BDF(debug=False)
        cid = 7
        origin = [1.135, .089237, -.0676]
        zaxis = [.135, .089237, -.0676]
        xzplane = [1.135, .089237, .9324]
        coord2 = model2.add_cord2r(cid, origin, zaxis, xzplane, rid=0, comment='cord2r')
        coord2.comment = ''
        assert coord == coord2, 'coord:\n%r\ncoord2:\n%r' % (str(coord), str(coord2))

    def test_coord_xform_a(self):
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

        fxyz = [0., -Fy, 0.]
        mxyz = [0., 0., 0.]
        cid_new = CORD2R.add_op2_data(data=data)
        model = None

        fxyz_local, mxyz_local = TransformLoadWRT(fxyz, mxyz, cid0, cid_new,
                                                  model)

        r = array([Lx, Ly, Lz])
        F = array([0., -Fy, 0.])
        M = cross(r, F)
        self.assertTrue(array_equal(fxyz_local, F), 'expected=%s actual=%s' % (F, fxyz_local))
        self.assertTrue(array_equal(mxyz_local, cross(r, F)), 'expected=%s actual=%s' % (M, mxyz_local))

    def test_coord_xform_b(self):
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

        fxyz = [0., -Fy, 0.]
        mxyz = [0., 0., 0.]
        cid_new = CORD2R.add_op2_data(data=data)
        model = None

        fxyz_local, mxyz_local = TransformLoadWRT(fxyz, mxyz, cid0, cid_new,
                                                  model)
        r = array([Lx, Ly, Lz])
        F = array([0., -Fy, 0.])
        M = cross(r, F)
        self.assertTrue(array_equal(fxyz_local, F), 'expected=%s actual=%s' % (F, fxyz_local))
        self.assertTrue(array_equal(mxyz_local, M), 'expected=%s actual=%s' % (M, mxyz_local))

    def test_coord_adding(self):
        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [1., 0., 0.]
        unused_cid1 = CORD2R(cid=1, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane,
                             comment='cord2r')

        xaxis = [1., 0., 0.]
        yaxis = [0., 1., 0.]
        zaxis = [0., 0., 1.]
        xz_plane = [1., 0., 1.]
        yz_plane = [0., 1., 1.]
        xy_plane = [1., 1., 0.]
        # x-axis
        unused_cid2 = CORD2R.add_axes(cid=2, rid=0, origin=origin,
                                      xaxis=xaxis, yaxis=None, zaxis=None,
                                      xyplane=None, yzplane=None, xzplane=xz_plane)

        unused_cid3 = CORD2R.add_axes(cid=2, rid=0, origin=origin,
                                      xaxis=xaxis, yaxis=None, zaxis=None,
                                      xyplane=xy_plane, yzplane=None, xzplane=None)

        # y-axis
        unused_cid4 = CORD2R.add_axes(cid=4, rid=0, origin=origin,
                                      xaxis=None, yaxis=yaxis, zaxis=None,
                                      xyplane=xy_plane, yzplane=None, xzplane=None)

        unused_cid5 = CORD2R.add_axes(cid=5, rid=0, origin=origin,
                                      xaxis=None, yaxis=yaxis, zaxis=None,
                                      xyplane=None, yzplane=yz_plane, xzplane=None)

        # z-axis
        unused_cid4 = CORD2R.add_axes(cid=4, rid=0, origin=origin,
                                      xaxis=None, yaxis=None, zaxis=zaxis,
                                      xyplane=None, yzplane=None, xzplane=xz_plane)

        unused_cid5 = CORD2R.add_axes(cid=5, rid=0, origin=origin,
                                      xaxis=None, yaxis=None, zaxis=zaxis,
                                      xyplane=None, yzplane=yz_plane, xzplane=None)

        # ijk
        unused_cid6 = CORD2R.add_ijk(cid=6, rid=0, origin=origin, i=xaxis, j=yaxis, k=None)
        unused_cid7 = CORD2R.add_ijk(cid=7, rid=0, origin=origin, i=xaxis, j=None, k=zaxis)
        unused_cid8 = CORD2R.add_ijk(cid=8, rid=0, origin=origin, i=None, j=yaxis, k=zaxis)
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
        cards_list = []
        cards_dict, card_count = model.get_bdf_cards_dict(bulk_data_lines)
        model._parse_cards(cards_list, cards_dict, card_count)

        #print(model.card_count)
        assert model.card_count['CORD1R'] == 1, model.card_count
        assert model.card_count['CORD2C'] == 4, model.card_count
        assert model.card_count['CORD2R'] == 11, model.card_count
        assert model.card_count['GRID'] == 3, model.card_count
        model.cross_reference()
        for unused_cid, coord in sorted(model.coords.items()):
            assert coord.i is not None, coord

    def test_define_coords_from_axes(self):
        """define_coord_e123"""
        model = BDF(debug=False)
        cord2_type = 'CORD2R'
        cid = 1
        origin = [0., 0., 0.]
        unused_rid = 0
        xaxis = [1., 0., 0.]
        xzplane = [0., 0., 1.]
        define_coord_e123(model, cord2_type, cid, origin, rid=0,
                          xaxis=xaxis, yaxis=None, zaxis=None,
                          xyplane=None, yzplane=None, xzplane=xzplane, add=True)

        cid = 2
        xaxis = [1., 0., 0.]
        xyplane = [1., 0., 1.]
        define_coord_e123(model, cord2_type, cid, origin, rid=0,
                          xaxis=xaxis, yaxis=None, zaxis=None,
                          xyplane=xyplane, yzplane=None, xzplane=None, add=True)
        #yaxis = [0., 1., 0.]
        xyplane = [1., 1., 0.]
        define_coord_e123(model, cord2_type, cid, origin, rid=0,
                          xaxis=xaxis, yaxis=None, zaxis=None,
                          xyplane=xyplane, yzplane=None, xzplane=None, add=True)

    def test_add_coord_cards(self):
        """tests the ``add_card`` method"""
        model = BDF(debug=False)
        fields = ['CORD1R',
                  10, 1, 2, 3,
                  11, 7, 8, 9]
        model.add_card(fields, 'CORD1R')

        fields = ['CORD1R',
                  12, 1, 2, 3,
                  13, 7, 8, 9]
        model.add_card(fields, 'CORD1S')

        fields = ['CORD1R',
                  14, 1, 2, 3,
                  15, 7, 8, 9]
        model.add_card(fields, 'CORD1C')
        model.pop_parse_errors()
        #print(model.coords)
        model.pop_xref_errors()
        self.assertEqual(len(model.coords), 7)

    def test_transform(self):
        model = BDF(debug=False)
        log = model.log
        g1 = 1
        g2 = 2
        g3 = 3
        cord1r = model.add_cord1r(1, g1, g2, g3, comment='')
        cord1c = model.add_cord1c(2, g1, g2, g3, comment='')
        cord1s = model.add_cord1s(3, g1, g2, g3, comment='')

        origin = [0., 0., 0.]
        zaxis = [0., 0., 1.]
        xzplane = [1., 0., 0.]
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 0., 1.])
        model.add_grid(3, [1., 0., 1.])

        cord2r = model.add_cord2r(4, origin, zaxis, xzplane, rid=0, setup=True, comment='')
        cord2c = model.add_cord2c(5, origin, zaxis, xzplane, rid=0, setup=True, comment='')
        cord2s = model.add_cord2s(6, origin, zaxis, xzplane, rid=0, setup=True, comment='')

        cord1r.raw_fields()
        cord1c.raw_fields()
        cord1s.raw_fields()
        cord2r.raw_fields()
        cord2c.raw_fields()
        cord2s.raw_fields()

        g1 = 11
        g2 = 12
        g3 = 13
        model.add_grid(11, [0., 0., 0.], cp=1)
        model.add_grid(12, [0., 0., 1.], cp=1)
        model.add_grid(13, [1., 0., 1.], cp=1)

        cord1r = model.add_cord1r(101, g1, g2, g3, comment='')
        cord1c = model.add_cord1c(102, g1, g2, g3, comment='')
        cord1s = model.add_cord1s(103, g1, g2, g3, comment='')
        cord2r = model.add_cord2r(104, origin, zaxis, xzplane, rid=1, setup=True, comment='')
        cord2c = model.add_cord2c(105, origin, zaxis, xzplane, rid=1, setup=True, comment='')
        cord2s = model.add_cord2s(106, origin, zaxis, xzplane, rid=1, setup=True, comment='')

        model.cross_reference()
        xyz = [0., 0., 0.]
        p = [xyz]
        coord_to = cord2s
        for cid, coord in model.coords.items():
            if hasattr(coord, 'transform_node_to_global_array'):
                coord.transform_node_from_local_to_local_array(coord_to, xyz)
            else:
                log.warning(f'{coord.type} is missing transform_node_to_global_array')

            if hasattr(coord, 'coord_to_spherical'):
                coord.coord_to_spherical(xyz)
            else:
                log.warning(f'{coord.type} is missing coord_to_spherical')

            if hasattr(coord, 'coord_to_cylindrical'):
                coord.coord_to_cylindrical(xyz)
            else:
                log.warning(f'{coord.type} is missing coord_to_cylindrical')

            coord.global_to_basic(xyz)
            coord.transform_vector_to_global_array(p)
            coord.transform_node_to_global(xyz)
            coord.transform_vector_to_local(xyz)
            coord.coord_to_xyz_array(xyz)
            coord.global_to_local
            coord.local_to_global

    def test_gmcord(self):
        """tests GMCORD"""
        cid = 1
        entity = 'GMCURV'
        gm_ids = [3, 4]
        model = BDF(debug=False)
        gmcord = model.add_gmcord(cid, entity, gm_ids)
        gmcord.raw_fields()
        save_load_deck(model, run_convert=False)

    def test_cord3g(self):
        """tests the CORD3G card"""
        cid = 1
        #method_es = 'E313'
        method_es = 'E'
        method_int = 123
        form = 'EQN'
        thetas = [110, 111, 112]
        rid = 0

        cord3g_e = CORD3G(cid, method_es, method_int, form, thetas, rid,
                          comment='cord3g')
        fields = BDFCard(cord3g_e.raw_fields())
        cord3g_e.repr_fields()
        cord3g_e.add_card(fields)
        xyz = [0., 0., 0.]
        cord3g_e.coord3g_transform_to_global(xyz)

        method_es = 'S'
        cord3g_s = CORD3G(cid, method_es, method_int, form, thetas, rid,
                          comment='cord3g')
        fields = BDFCard(cord3g_s.raw_fields())
        cord3g_s.repr_fields()
        cord3g_s.add_card(fields)
        with self.assertRaises(NotImplementedError):  # TODO: add me
            cord3g_s.coord3g_transform_to_global(xyz)

    def test_create_coord_line(self):
        """tests creating a series of coordinate systems down an axis"""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [0., 100., 0.])
        eid = 10
        mid = 11
        nids = [1, 2]
        model.add_conrod(eid, mid, nids, A=1.0, j=0.0, c=0.0, nsm=0.0,
                         comment='')
        E = 3.0e7
        G = None
        nu = 0.3
        model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0, ge=0.0,
                       St=0.0, Sc=0.0, Ss=0.0, mcsid=0,
                       comment='')
        model.cross_reference()


        cid = 0
        unused_icd_transform, icp_transform, xyz_cp, nid_cp_cd = model.get_displacement_index_xyz_cp_cd(
            fdtype='float64', idtype='int32', sort_ids=True)
        nids = nid_cp_cd[:, 0]
        #unused_xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
            #xyz_cp, nids, icp_transform,
            #cid=cid, in_place=False, atol=1e-6)


        npoints = 50
        percents = np.linspace(0., 1., num=npoints, endpoint=True)
        p1 = np.array([50., -100, 0.])
        p2 = np.array([50., 100, 0.])

        cids = create_coords_along_line(model, p1, p2, percents, cid=0, axis=cid)
        cid_to_inids = get_nodes_along_axis_in_coords(
            model, nids, xyz_cp, icp_transform,
            cids)
        make_monpnt1s_from_cids(model, nids, cids, cid_to_inids)
        #model.write_bdf('spike.bdf')

def make_tri(model):
    model.add_grid(1, [0., 0., 0.])
    model.add_grid(3, [0., 0., 1.])
    model.add_grid(4, [1., 0., 1.])
    eid = 100
    pid = 10
    nids = [1, 3, 4]
    mid = 1000
    E = 3.0e7
    G = None
    nu = 0.3
    model.add_ctria3(eid, pid, nids)
    model.add_mat1(mid, E, G, nu)
    model.add_pshell(pid, mid1=mid, t=0.1)

def get_nodes(grids, grids_expected, coords):
    """
    Create each input grid/coord

    Loop over the expected grids, use the provided coordinate system
    and verify that the xyz location is correct.

    Parameters
    ----------
    grids : List[grid]
        the grids
        grid : List[nid, cp, x, y, z]
            the GRID fields
    coords : List[coord]
        coord : List[int cid, rid, origin, zaxis, xaxis]
            the coordinate system to add to the model
    grids_expected : List[grid]
        the expected grids
        grid : List[nid, cp, x, y, z]
            the GRID fields

    """
    model = BDF(debug=False)

    for grid in grids:
        (nid, cid, x, y, z) = grid
        model.add_card(['GRID', nid, cid, x, y, z], 'GRID')

    for coord in coords:
        (cid, rid, x, y, z) = coord
        model.add_card(['CORD2R', cid, rid] + x + y + z, 'CORD2R')
        #coord_obj = model.Coord(cid)

    model.cross_reference()
    save_load_deck(model, run_remove_unused=False, run_convert=False)

    for (i, grid) in enumerate(grids_expected):
        (nid, cid, x, y, z) = grid
        node = model.Node(nid)
        xyz_actual = node.get_position()
        xyz = array([x, y, z])

        msg = 'i=%s expected=%s actual=%s\n' % (i, xyz, xyz_actual)
        msg += 'n=%s grid=\n%s' % (nid, node)
        coord_ref = node.cp_ref
        msg += 'n=%s coord=\n%s' % (node.nid, coord_ref)

        if not allclose(xyz, xyz_actual):
            # TODO: this used to work, but the changing xref broke it somehow
            #       this block probably needs to be slightly updated
            while coord_ref.rid:
                msg += 'xyz=%s rcoord=\n%s' % (node.nid, coord_ref.rid)
                coord_ref = coord_ref.rid
            assert allclose(xyz, xyz_actual), msg


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
