import os
import unittest
from numpy import allclose

import pyNastran
from pyNastran.converters.fast.fgrid_reader import FGridReader

PKG_PATH = pyNastran.__path__[0]
TEST_PATH = os.path.join(PKG_PATH, 'converters', 'fast')

class TestFast(unittest.TestCase):

    def test_fgrid_io_01(self):
        infile_name = os.path.join(TEST_PATH, 'flow_demo1', 'om6inviscid.fgrid')

        fgrid = FGridReader(log=None, debug=None)
        fgrid.read_fgrid(infile_name)

        #self.nodes = nodes
        #self.tris = tris
        #self.tets = tets
        assert len(fgrid.nodes) == 2800, 'nodes=%s' % len(fgrid.nodes)
        assert len(fgrid.tris) == 2004, 'tris=%s' % len(fgrid.tris)
        assert len(fgrid.tets) == 13576, 'tets=%s' % len(fgrid.tets)


if __name__ == '__main__':  # pragma: no cover
    import time
    time0 = time.time()
    unittest.main()
    print("dt = %s" % (time.time() - time0))
