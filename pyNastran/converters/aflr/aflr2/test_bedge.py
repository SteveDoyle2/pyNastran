from __future__ import print_function
import os
import unittest

import warnings
import numpy as np
warnings.simplefilter('always')
np.seterr(all='raise')

import pyNastran
from pyNastran.converters.aflr.aflr2.aflr2 import read_bedge

pkg_path = pyNastran.__path__[0]
test_path = os.path.join(pkg_path, 'converters', 'aflr', 'aflr2')

class TestBEdge(unittest.TestCase):

    def test_bedge_1(self):
        """tests the m3 model"""
        bedge_filename = os.path.join(test_path, 'm3.bedge')
        bdf_filename = os.path.join(test_path, 'm3.bdf')
        model = read_bedge(bedge_filename, beta_reverse=179.7, log=None, debug=False)
        model.write_nastran(bdf_filename)

        assert len(model.nodes) == 858, 'nodes=%s' % len(model.nodes)
        assert len(model.bars) == 858, 'nbars=%s' % len(model.bars)
        os.remove(bdf_filename)


def main():  # pragma: no cover
    import time
    t0 = time.time()
    unittest.main()
    print("dt = %s" % (time.time() - t0))

if __name__ == '__main__':  # pragma: no cover
    main()
