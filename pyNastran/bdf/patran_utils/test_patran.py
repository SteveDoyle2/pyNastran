from __future__ import print_function
import os
import unittest

from pyNastran.bdf.patran_utils.read_patran import read_patran
import pyNastran

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models', 'patran_fmt')


class TestPatran(unittest.TestCase):
    def test_read_patran(self):
        """tests the nodal patran read function"""
        #bdf_filename = os.path.join(MODEL_PATH, '0012_20.bdf')
        node_filename = os.path.join(MODEL_PATH, 'normals.nod')
        data = read_patran(node_filename, fdtype='float64', idtype='int32')
        nids = data['nids']
        assert nids.min() == 1, nids
        assert nids.max() == 2214, nids

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
