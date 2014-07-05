from numpy import array, allclose
import unittest

from pyNastran.bdf.bdf import BDF, BDFCard, CORD1R, CORD1C, CORD1S

bdf = BDF()
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
        self.getNodes(grids, grids_expected, coords)

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

        coords = [  # cid,rid, origin,      zaxis,     xaxis
            [1, 0, [1., 1., 1.], [1., 1., 2.], [2., 1., 1.]],
        ]
        self.getNodes(grids, grids_expected, coords)

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
                     [1, 1., 1.,  0., 0.],
                     [2, 1., 0., -1., 0.],
                     [3, 1., 0.,  0., 1.],
                     [4, 1., 1., -1., 1.],
                     [5, 1., 0., -1., 1.],
                 ]

        coords = [  # cid, rid, origin,      zaxis,     xaxis
                   [1, 0, [0., 0., 0.], [1., 0., 0.], [0., 0., 1.]],
                 ]
        self.getNodes(grids, grids_expected, coords)

    def test_rotate2(self):   # passes
        grids = [
                     [1, 1, 0., 0., 1.],  # nid, cid, x,y,z
                     [2, 1, 0., 1., 0.],
                     [3, 1, 1., 0., 0.],
                     [4, 1, 1., 1., 1.],
                     [5, 1, 1., 1., 0.],
                 ]
        grids_expected = [
                     [1, 1, 0.,  0., -1.],
                     [2, 1, 0., -1.,  0.],
                     [3, 1, 1.,  0.,  0.],
                     [4, 1, 1., -1., -1.],
                     [5, 1, 1., -1.,  0.],
                 ]

        coords = [  # cid, rid, origin,     zaxis        xaxis
                   [1, 0, [0., 0., 0.], [0., 0., -1.], [1., 0., 0.]],
                 ]
        self.getNodes(grids, grids_expected, coords)

    def test_rotate3(self):  # passes
        #print('test_rotate3')
        grids = [
                     [1, 1, 0., 0., 1.],
                     [2, 1, 0., 1., 0.],
                     [3, 1, 1., 0., 0.],
                     [4, 1, 1., 1., 1.],
                     [5, 1, 1., 1., 0.],
                 ]
        grids_expected = [
                     [1, 1,  0., 0., -1.],
                     [2, 1,  0., 1.,  0.],
                     [3, 1, -1., 0.,  0.],
                     [4, 1, -1., 1., -1.],
                     [5, 1, -1., 1.,  0.],
                 ]

        coords = [  # cid, rid, origin,      zaxis          xaxis
                   [1, 0,      [0., 0., 0.], [0., 0., -1.], [-1., 0., 0.]],
                 ]
        self.getNodes(grids, grids_expected, coords)

    def test_rid_1(self):
        #print('test_rid_1')
        grids = [
                    [1, 2, 10., 5., 3.],  # nid, cid, x,y,z
                    [2, 3, 10., 5., 3.],
                 ]
        grids_expected = [
                    [1, 1, 11., 6., 4.],
                    [2, 1, 11., 6., 4.],
                 ]

        coords = [  # cid, rid, origin,     zaxis        xaxis
                   [1, 0, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.]],  # cid=1
                   [2, 1, [1., 1., 1.], [1., 1., 2.], [2., 1., 1.]],  # cid=2
                  #[2, 1, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.]],  # cid=2,equiv

                   [3, 0, [1., 1., 1.], [1., 1., 2.], [2., 1., 1.]],  # cid=3
                  #[4, 0, [0., 0., 0.], [0., 0., 1.], [1., 0., 0.]],  # cid=3,equiv
                 ]
        self.getNodes(grids, grids_expected, coords)

    def test_cord1r_01(self):
        lines = ['cord1r,2,1,4,3']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CORD1R(card)
        self.assertEquals(card.Cid(), 2)
        #self.assertEquals(card.Rid(), 1)
        card.write_bdf(size, 'dummy')
        card.rawFields()

    def test_cord1c_01(self):
        lines = ['cord1c,2,1,4,3']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CORD1C(card)
        self.assertEquals(card.Cid(), 2)
        #self.assertEquals(card.Rid(), 1)
        card.write_bdf(size, 'dummy')
        card.rawFields()

    def test_cord1s_01(self):
        lines = ['cord1s,2,1,4,3']
        card = bdf.process_card(lines)
        card = BDFCard(card)

        size = 8
        card = CORD1S(card)
        self.assertEquals(card.Cid(), 2)
        #self.assertEquals(card.Rid(), 1)
        card.write_bdf(size, 'dummy')
        card.rawFields()

    def test_cord2r_1(self):
        grid = ['GRID       20143       7 -9.31-4  .11841 .028296']
        coord = ['CORD2R         7           1.135 .089237  -.0676    .135 .089237  -.0676',
                 '           1.135 .089237   .9324']

        model = BDF()
        card = model.process_card(grid)
        model.add_card(card, card[0])

        card = model.process_card(coord)
        model.add_card(card, card[0])
        model.cross_reference()

        g = model.Node(20143)
        #print(g.Position(debug=False))

        # by running it through Patran...
        #GRID     20143          1.1067  .207647 -.068531
        diff = g.Position() - array([1.106704, .207647, -0.068531])

        msg = 'diff=%s' % diff
        assert allclose(diff, 0.), msg
        coord = model.Coord(7)
        coord.T()

    #def makeNodes(self, grids, coords):
        #grids2 = []

    def getNodes(self, grids, grids_expected, coords, debug=False):
        model = BDF(debug=False)

        for grid in grids:
            (nid, cid, x, y, z) = grid
            model.add_card(['GRID', nid, cid, x, y, z], 'GRID')
            gridObj = model.Node(nid)
            if debug:
                print(gridObj)

        for coord in coords:
            #print coord
            (cid, rid, x, y, z) = coord
            model.add_card(['CORD2R', cid, rid] + x + y + z, 'CORD2R')
            coordObj = model.Coord(cid)
            if debug:
                print(coordObj)

        model.cross_reference()

        for (i, grid) in enumerate(grids_expected):
            (nid, cid, x, y, z) = grid
            node = model.Node(nid)
            pos = node.Position()
            n = array([x, y, z])

            msg = 'i=%s expected=%s actual=%s\n' % (i, n, pos)
            msg += 'n=%s grid=\n%s' % (nid, node)
            coord = node.cp
            msg += 'n=%s coord=\n%s' % (node.nid, coord)
            while coord.rid:
                msg += 'n=%s rcoord=\n%s' % (node.nid, coord.rid)
                coord = coord.rid
            assert allclose(n, pos), msg

if __name__ == '__main__':
    unittest.main()
