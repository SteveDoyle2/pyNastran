from __future__ import print_function
import os
import unittest

#from pyNastran.bdf.bdf import read_bdf
from pyNastran.bdf.patran_utils.read_patran import read_patran
import pyNastran

ROOTPATH = pyNastran.__path__[0]
MODELPATH = os.path.join(ROOTPATH, '..', 'models', 'patran_fmt')


class TestPatran(unittest.TestCase):
    def test_read_patran(self):
        """tests the nodal patran read function"""
        #bdf_filename = os.path.join(MODELPATH, '0012_20.bdf')
        node_filename = os.path.join(MODELPATH, 'normals.nod')
        data = read_patran(node_filename)
        nids = data['nids']
        assert nids.min() == 1, nids
        assert nids.max() == 2214, nids

if __name__ == '__main__':
    unittest.main()
