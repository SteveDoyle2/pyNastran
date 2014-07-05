from numpy import array, allclose
import unittest

from pyNastran.bdf.bdf import BDF, BDFCard, CORD1R, CORD1C, CORD1S, CORD2R, CORD2C, CORD2S

bdf = BDF()  # don't load this up with stuff
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

    def test_cord2c_01(self):
        lines = [
            'CORD2C*                3               0              0.              0.',
            '*                     0.              0.              0.              1.*',
            '*                     1.              0.              1.'
        ]
        model = BDF()
        card = model.process_card(lines)
        card = BDFCard(card)
        card = CORD2C(card)
        model.add_coord(card)

        lines = [
            'CORD2R         4       3     10.      0.      5.     10.     90.      5.',
            '             10.      0.      6.'
        ]
        card = model.process_card(lines)
        card = BDFCard(card)
        card = CORD2R(card)
        model.add_coord(card)
        model.cross_reference()

        cord2r = model.Coord(4)

        #print('i = ', cord2r.i)
        #print('j = ', cord2r.j)
        #print('k = ', cord2r.k)
        self.assertTrue(allclose(cord2r.i, array([0., 0., 1.])))
        delta = cord2r.j - array([1., 1., 0.]) / 2**0.5
        self.assertTrue(allclose(cord2r.j, array([1., 1., 0.]) / 2**0.5), str(delta))
        delta = cord2r.k - array([-1., 1., 0.]) / 2**0.5
        self.assertTrue(allclose(cord2r.k, array([-1., 1., 0.]) / 2**0.5), str(delta))

    def test_cord2_rcs_01(self):
        model = BDF()
        cards = [
             [#'$ Femap with NX Nastran Coordinate System 10 : rectangular defined in a rectangular',
                 'CORD2R*               10               0             10.              5.',
                 '*                     3.   10.3420201433   4.53015368961   3.81379768136*',
                 '*          10.7198463104   5.68767171433   3.09449287122',],
             [#'$ Femap with NX Nastran Coordinate System 11 : cylindrical defined in rectangular',
                 'CORD2C*               11               0              7.              3.',
                 '*                     9.   7.64278760969   2.73799736977   9.71984631039*',
                 '*          7.75440650673   3.37968226211   8.46454486422',],
             [#'$ Femap with NX Nastran Coordinate System 12 : spherical defined in rectangular',
                 'CORD2S*               12               0             12.              8.',
                 '*                     5.   12.6427876097   7.86697777844   5.75440650673*',
                 '*          12.6634139482   8.58906867688   4.53861076379',],
             ['GRID*                 10              10   42.9066011565   34.2422137135',
                 '*          28.6442730262               0',],
             ['GRID*                 11              11   48.8014631871   78.8394787869',
                 '*          34.6037164304               0',],
             ['GRID*                 12              12   58.0775343829   44.7276544324',
                 '*          75.7955331161               0',],
            ]
        for lines in cards:
            card = model.process_card(lines)
            #print('card*, ', card, type(card))
            model.add_card(card, card[0])
        model.cross_reference()
        for nid in model.nodes:
            #print(model.Node(nid).Position())
            a = array([30.,40.,50.])
            b = model.Node(nid).Position()
            self.assertTrue(allclose(array([30.,40.,50.]), model.Node(nid).Position()), str(a-b))

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
