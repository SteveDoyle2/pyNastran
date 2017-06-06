from __future__ import print_function
import unittest
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.tables.geom.geom4 import read_rbe3s_from_idata_fdata

class TestOP2GeomUnit(unittest.TestCase):
    """lots of small OP2-Geom tests"""
    def test_rbe3(self):
        """
        data = [99           99 123456 1.0    123    44    45  48  49  -1    -3]
        data = [61           71 123456 1.0    123    70    75  77      -1    -3
                62           71 123456 1.0    123    58    59  72      -1    -3]
        data = [1001100 1001100 123456 1.0 123456 10011 10002          -1 -2 -3
                1002500 1002500 123456 1.0 123456 10025 10020          -1 -2 -3]
                eid     refg    refc   wt  c      g     ...
        """
        model = BDF()
        data = [99, 99, 123456, 1.0, 123, 44, 45, 48, 49, -1, -3]
        rbes = read_rbe3s_from_idata_fdata(
            model, np.array(data, dtype='int32'), np.array(data, dtype='float32'))
        assert len(rbes) == 1, rbes

        data = [61, 71, 123456, 1.0, 123, 70, 75, 77, -1, -3,
                62, 71, 123456, 1.0, 123, 58, 59, 72, -1, -3]
        rbes = read_rbe3s_from_idata_fdata(
            model, np.array(data, dtype='int32'), np.array(data, dtype='float32'))
        assert len(rbes) == 2, rbes

        data = [11, 10, 123456, 1.0, 456, 111, 100, -1, -2, -3,
                25, 50, 123456, 1.0, 456, 125, 103, -1, -2, -3]
        rbes = read_rbe3s_from_idata_fdata(
            model, np.array(data, dtype='int32'), np.array(data, dtype='float32'))
        assert len(rbes) == 0, rbes  ## TODO: not supported b/c the -2; should be 2


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
