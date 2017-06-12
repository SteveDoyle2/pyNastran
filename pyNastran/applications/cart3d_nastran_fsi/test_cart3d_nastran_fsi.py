from __future__ import print_function
import os
import unittest
import numpy as np
from pyNastran.applications.cart3d_nastran_fsi.math_functions import (
    triangle_area_centroid_normal, Centroid, Normal)
from pyNastran.applications.cart3d_nastran_fsi.map_deflections import Tet4, DeflectionReader
class Cart3dNastranFSI(unittest.TestCase):

    def test_acn(self):
        """tests triangle_area_centroid_normal"""
        n1 = np.array([0., 0., 0.])
        n2 = np.array([1., 1., 1.])
        n3 = np.array([1., 0., 0.])
        #n4 = np.array([5., 3., 0.])
        n5 = np.array([2., 0., 4.])

        n2 = np.array([0., 1., 0.])
        c1 = Centroid(n1, n2, n3)
        n = Normal(n5, n2)
        print("norm = ", n, np.linalg.norm(n))

        area, centroid, normal = triangle_area_centroid_normal([n1, n2, n3])
        print("area=%s centroid=%s normal=%s" % (area, centroid, normal))

    def test_tet(self):
        """tests a tet object"""
        b = [10., 0., 0.]
        a = [0., 10., 0.]
        c = [0., 0., 10.]
        d = [0., 0., 0.]
        m1 = [1., 1., 1.]
        m2 = [2., 2., 2.]
        tet = Tet4(a, b, c, d)
        print("volume = ", tet.volume())
        print("is_internal = \n%s\n" % (tet.is_internal_node(m1)))
        print("is_internal = \n%s"   % (tet.is_internal_node(m2)))

    def test_deflections(self):
        """reads an op2"""
        infilename = os.path.join('op2reader', 'solid_shell_bar.op2')
        deflections = {}
        op2 = DeflectionReader(infilename)
        #op2.print_displacement()
        #displacements = op2.convert_displacements()  # what is this...

        #for gridID, disp in sorted(displacements.items()):
           #print("gridID=%s disp=%s" % (gridID, disp))


if __name__ == '__main__': # pragma: no cover
    unittest.main()
