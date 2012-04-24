## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 

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
        print 'test_rid_2'
        grids  = [
                     [2,    10., 5.,  3.],
                 ]
        gridsExpected  = [
                     ['x',  11., 6.,  4.],
                 ]

        coords = [   #rid origin,     zaxis        xaxis
                   [  0,  [0.,0.,0.], [0.,0.,-1.], [1.,0.,0.]  ],  # cid=1
                   [  1,  [1.,1.,1.], [1.,1.,2.],  [2.,1., 1.]   ],  # cid=2
                 ]
        self.getNodes(grids,gridsExpected,coords)


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
        #print "hi"

if __name__=='__main__':
    unittest.main()
