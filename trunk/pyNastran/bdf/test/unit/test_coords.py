
from numpy import array,allclose
import unittest

from pyNastran.bdf.bdf import BDF

class Tester(unittest.TestCase):
    def test_same(self):  # passes
        grids  = [
                     [0,  0.,  0.,  1.],
                     [0,  0.,  1.,  0.],
                     [0,  1.,  0.,  0.],
                     [0,  1.,  1.,  1.],
                     [0,  1.,  1.,  0.],
                 ]
        gridsExpected = grids
        coords = []
        self.getNodes(grids,gridsExpected,coords)

    def test_shift(self):  # passes
        grids  = [
                     [1,  0.,  0.,  1.],
                     [1,  0.,  1.,  0.],
                     [1,  1.,  0.,  0.],
                     [1,  1.,  1.,  1.],
                     [1,  1.,  1.,  0.],
                 ]
        gridsExpected  = [
                     [1,  1.,  1.,  2.],
                     [1,  1.,  2.,  1.],
                     [1,  2.,  1.,  1.],
                     [1,  2.,  2.,  2.],
                     [1,  2.,  2.,  1.],
                 ]

        coords = [   #rid origin,      zaxis,     xaxis
                   [  0,  [1.,1.,1.], [1.,1.,2.], [2.,1.,1.]   ],
                 ]
        self.getNodes(grids,gridsExpected,coords)


    def test_rotate(self):  # passes
        grids  = [
                     [1,  0.,  0.,  1.],
                     [1,  0.,  1.,  0.],
                     [1,  1.,  0.,  0.],
                     [1,  1.,  1.,  1.],
                     [1,  1.,  1.,  0.],
                 ]
        gridsExpected  = [
                     #     y    z   x
                     [1.,  1.,  0,  0.],
                     [1.,  0., -1,  0.],
                     [1.,  0.,  0,  1.],
                     [1.,  1., -1,  1.],
                     [1.,  0., -1,  1.],
                 ]

        coords = [   #rid origin,      zaxis,     xaxis
                   [  0,  [0.,0.,0.], [1.,0.,0.], [0.,0.,1.]   ],
                 ]
        self.getNodes(grids,gridsExpected,coords)


    def test_rotate2(self):   # passes
        grids  = [
                     [1,  0.,  0.,  1.],
                     [1,  0.,  1.,  0.],
                     [1,  1.,  0.,  0.],
                     [1,  1.,  1.,  1.],
                     [1,  1.,  1.,  0.],
                 ]
        gridsExpected  = [
                     [1,  0.,  0., -1.],
                     [1,  0., -1.,  0.],
                     [1,  1.,  0.,  0.],
                     [1,  1., -1., -1.],
                     [1,  1., -1.,  0.],
                 ]

        coords = [   #rid origin,     zaxis        xaxis
                   [  0,  [0.,0.,0.], [0.,0.,-1.], [1.,0.,0.]   ],
                 ]
        self.getNodes(grids,gridsExpected,coords)

    def Atest_rotate3(self):  # passes
        print 'test_rotate3'
        grids  = [
                     [1,  0.,  0.,  1.],
                     [1,  0.,  1.,  0.],
                     [1,  1.,  0.,  0.],
                     [1,  1.,  1.,  1.],
                     [1,  1.,  1.,  0.],
                 ]
        gridsExpected  = [
                     [1,   0.,  0., -1.],
                     [1,   0.,  1.,  0.],
                     [1,  -1.,  0.,  0.],
                     [1,  -1.,  1., -1.],
                     [1,  -1.,  1.,  0.],
                 ]

        coords = [   #rid origin,     zaxis        xaxis
                   [  0,  [0.,0.,0.], [0.,0.,-1.], [-1.,0.,0.]   ],
                 ]
        self.getNodes(grids,gridsExpected,coords)


    def test_rid_1(self):
        print 'test_rid_1'
        grids  = [
                     [2,    10., 5.,  3.],
                 ]
        gridsExpected  = [
                     ['x',  11., 6.,  4.],
                 ]

        coords = [   #rid origin,     zaxis        xaxis
                   [  0,  [0.,0.,0.], [0.,0.,-1.], [1.,0.,0.]  ],  # cid=1
                   [  1,  [1.,1.,1.], [1.,1., 2.], [2.,1.,1.]  ],  # cid=2
                 ]
        self.getNodes(grids,gridsExpected,coords)

    def test_cord2r_1(self):
        grid = ['GRID       20143       7 -9.31-4  .11841 .028296']
        coord = ['CORD2R         7           1.135 .089237  -.0676    .135 .089237  -.0676',
                 '           1.135 .089237   .9324']
        
        mesh = BDF()
        card = mesh.processCard(grid)
        mesh.addCard(card,card[0])
        card = mesh.processCard(coord)
        mesh.addCard(card,card[0])
        mesh.crossReference()

        g = mesh.Node(20143)
        print g.Position(debug=False)
        diff = g.Position() - array([1.106704,.207647,-0.068531])
        
        assert allclose(diff,0.)

    def makeNodes(self,grids,coords):
        grids2 = []
        
    def getNodes(self,grids,gridsExpected,coords,debug=False):
        mesh = BDF(debug=False)

        for nid,grid in enumerate(grids):
            (cid,x,y,z) = grid
            mesh.addCard(['GRID',nid+1,cid,x,y,z],'GRID')
            gridObj = mesh.Node(nid+1)
            if debug:
                print gridObj

        for cid,coord in enumerate(coords):
            #print coord
            (rid,x,y,z) = coord
            obj  = mesh.addCard(['CORD2R',cid+1,rid]+x+y+z, 'CORD2R')
            coordObj = mesh.Coord(cid+1)
            if debug:
                print coordObj

        mesh.crossReference()

        for i,grid in enumerate(gridsExpected):
            (cid,x,y,z) = grid
            node = mesh.Node(i+1)
            pos  = node.Position()
            n = array([x,y,z])
            
            msg = 'expected=%s actual=%s' %(n,pos)
            assert allclose(n,pos),msg
        ###

if __name__=='__main__':
    unittest.main()
