from numpy import array, allclose
import unittest

from pyNastran.bdf.bdf import BDF


class TestCoords(unittest.TestCase):
    def test_same(self):  # passes
        grids = [
            [0, 0., 0., 1.],
            [0, 0., 1., 0.],
            [0, 1., 0., 0.],
            [0, 1., 1., 1.],
            [0, 1., 1., 0.],
        ]
        gridsExpected = grids
        coords = []
        self.getNodes(grids, gridsExpected, coords)

    def test_shift(self):  # passes
        grids = [
            [1, 0., 0., 1.],
            [1, 0., 1., 0.],
            [1, 1., 0., 0.],
            [1, 1., 1., 1.],
            [1, 1., 1., 0.],
        ]
        gridsExpected = [
            [1, 1., 1., 2.],
            [1, 1., 2., 1.],
            [1, 2., 1., 1.],
            [1, 2., 2., 2.],
            [1, 2., 2., 1.],
        ]

        coords = [  # rid origin,      zaxis,     xaxis
            [0, [1., 1., 1.], [1., 1., 2.], [2., 1., 1.]],
                 ]
        self.getNodes(grids, gridsExpected, coords)

    def test_rotate(self):  # passes
        grids = [
                     [1, 0., 0., 1.],
                     [1, 0., 1., 0.],
                     [1, 1., 0., 0.],
                     [1, 1., 1., 1.],
                     [1, 1., 1., 0.],
                 ]
        gridsExpected = [
                     #     y    z   x
                     [1., 1., 0, 0.],
                     [1., 0., -1, 0.],
                     [1., 0., 0, 1.],
                     [1., 1., -1, 1.],
                     [1., 0., -1, 1.],
                 ]

        coords = [  # rid origin,      zaxis,     xaxis
                   [0, [0., 0., 0.], [1., 0., 0.], [0., 0., 1.]],
                 ]
        self.getNodes(grids, gridsExpected, coords)

    def test_rotate2(self):   # passes
        grids = [
                     [1, 0., 0., 1.],
                     [1, 0., 1., 0.],
                     [1, 1., 0., 0.],
                     [1, 1., 1., 1.],
                     [1, 1., 1., 0.],
                 ]
        gridsExpected = [
                     [1, 0., 0., -1.],
                     [1, 0., -1., 0.],
                     [1, 1., 0., 0.],
                     [1, 1., -1., -1.],
                     [1, 1., -1., 0.],
                 ]

        coords = [  # rid origin,     zaxis        xaxis
                   [0, [0., 0., 0.], [0., 0., -1.], [1., 0., 0.]],
                 ]
        self.getNodes(grids, gridsExpected, coords)

    def test_rotate3(self):  # passes
        #print('test_rotate3')
        grids = [
                     [1, 0., 0., 1.],
                     [1, 0., 1., 0.],
                     [1, 1., 0., 0.],
                     [1, 1., 1., 1.],
                     [1, 1., 1., 0.],
                 ]
        gridsExpected = [
                     [1, 0., 0., -1.],
                     [1, 0., 1., 0.],
                     [1, -1., 0., 0.],
                     [1, -1., 1., -1.],
                     [1, -1., 1., 0.],
                 ]

        coords = [  # rid origin,     zaxis        xaxis
                   [0, [0., 0., 0.], [0., 0., -1.], [-1., 0., 0.]],
                 ]
        self.getNodes(grids, gridsExpected, coords)

    def off_test_rid_1(self):  # did i mess up the transform???
        #print('test_rid_1')
        grids = [
                     [2, 10., 5., 3.],  # cid, x,y,z
                    #[3,    10., 5.,  3.],
                 ]
        gridsExpected = [
                     ['x', 11., 6., 4.],  # ??? x,y,z
                    #['x',  11., 6.,  4.],
                 ]

        coords = [  # rid origin,     zaxis        xaxis
                   [0, [0., 0., 0.], [0., 0., -1.], [1., 0., 0.]],  # cid=1
                   [1, [1., 1., 1.], [1., 1., 2.], [2., 1., 1.]],  # cid=2
                  #[  1,  [0.,0.,0.], [0.,0., 1.], [1.,0.,0.]  ],  # cid=2,equiv

                   [0, [1., 1., 1.], [1., 1., 2.], [2., 1., 1.]],  # cid=3
                  #[  0,  [0.,0.,0.], [0.,0., 1.], [1.,0.,0.]  ],  # cid=3,equiv
                 ]
        self.getNodes(grids, gridsExpected, coords)

    def test_cord2r_1(self):
        grid = ['GRID       20143       7 -9.31-4  .11841 .028296']
        coord = ['CORD2R         7           1.135 .089237  -.0676    .135 .089237  -.0676',
                 '           1.135 .089237   .9324']

        mesh = BDF()
        card = mesh.processCard(grid)
        mesh.add_card(card, card[0])
        card = mesh.processCard(coord)
        mesh.add_card(card, card[0])
        mesh.cross_reference()

        g = mesh.Node(20143)
        #print(g.Position(debug=False))
        diff = g.Position() - array([1.106704, .207647, -0.068531])

        assert allclose(diff, 0.)

    def makeNodes(self, grids, coords):
        grids2 = []

    def getNodes(self, grids, gridsExpected, coords, debug=False):
        mesh = BDF(debug=False)

        for (nid, grid) in enumerate(grids):
            (cid, x, y, z) = grid
            mesh.add_card(['GRID', nid + 1, cid, x, y, z], 'GRID')
            gridObj = mesh.Node(nid + 1)
            if debug:
                print(gridObj)

        for (cid, coord) in enumerate(coords):
            #print coord
            (rid, x, y, z) = coord
            obj = mesh.add_card(
                ['CORD2R', cid + 1, rid] + x + y + z, 'CORD2R')
            coordObj = mesh.Coord(cid + 1)
            if debug:
                print(coordObj)

        mesh.cross_reference()

        for (i, grid) in enumerate(gridsExpected):
            (cid, x, y, z) = grid
            node = mesh.Node(i + 1)
            pos = node.Position()
            n = array([x, y, z])

            msg = 'expected=%s actual=%s\n' % (n, pos)
            msg += 'n=%s grid=\n%s' % (i + 1, node)
            coord = node.cp
            msg += 'n=%s coord=\n%s' % (node.nid, coord)
            while coord.rid:
                msg += 'n=%s rcoord=\n%s' % (node.nid, coord.rid)
                coord = coord.rid
            assert allclose(n, pos), msg
        ###

if __name__ == '__main__':
    unittest.main()
