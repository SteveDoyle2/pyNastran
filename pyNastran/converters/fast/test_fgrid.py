from six.moves import range
import os
from numpy import array_equal, allclose
import unittest

import pyNastran
from pyNastran.converters.fast.fgrid_reader import FGridReader
#from pyNastran.converters.cart3d.cart3d_to_nastran import cart3d_to_nastran_filename, cart3d_to_nastran_model

pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'converters', 'fast')

class TestFast(unittest.TestCase):

    def test_fgrid_io_01(self):
        infile_name = os.path.join(test_path, 'flow_demo1', 'om6inviscid.fgrid')

        fgrid = FGridReader(log=None, debug=False)
        fgrid.read_fgrid(infile_name)

        #self.nodes = nodes
        #self.tris = tris
        #self.tets = tets
        assert len(fgrid.nodes) == 2800, 'nodes=%s' % len(fgrid.nodes)
        assert len(fgrid.tris) == 2004, 'tris=%s' % len(fgrid.tris)
        assert len(fgrid.tets) == 13576, 'tets=%s' % len(fgrid.tets)


if __name__ == '__main__':  # pragma: no cover
    import time
    t0 = time.time()
    unittest.main()
    #test_1()
    #test_2()
    #test_3()
    print("dt = %s" % (time.time() - t0))
