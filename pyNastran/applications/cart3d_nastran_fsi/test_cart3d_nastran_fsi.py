from __future__ import print_function
import os
import unittest
import numpy as np
from pyNastran.applications.cart3d_nastran_fsi.math_functions import (
    triangle_area_centroid_normal, Centroid, Normal)


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
        #print("norm = ", n, np.linalg.norm(n))

        area, centroid, normal = triangle_area_centroid_normal([n1, n2, n3])
        #print("area=%s centroid=%s normal=%s" % (area, centroid, normal))


if __name__ == '__main__': # pragma: no cover
    unittest.main()
