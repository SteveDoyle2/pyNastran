#pylint: disable=R0201
from __future__ import print_function, unicode_literals
import unittest
import numpy as np
from pyNastran.op2.op2_geom import OP2Geom
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
        model = OP2Geom()
        data = [99, 99, 123456, 1.0, 123, 44, 45, 48, 49, -1, -3]
        rbes = read_rbe3s_from_idata_fdata(
            model, np.array(data, dtype='int32'), np.array(data, dtype='float32'))
        assert len(rbes) == 1, rbes

        rbes_expected = [
            ['RBE3', '99', '99', '123456', '1.', '123', '44', '45',
             '48', '49',]
        ]
        for rbe, rbe_expected in zip(rbes, rbes_expected):
            msg = 'rbe:\n%s\nexpected:\n%s' % (rbe.rstrip().split(), rbe_expected)
            assert rbe.rstrip().split() == rbe_expected, msg
        #--------------------------------------------------------------


        data = [61, 71, 123456, 1.0, 123, 70, 75, 77, -1, -3,
                62, 71, 123456, 1.0, 123, 58, 59, 72, -1, -3]
        rbes = read_rbe3s_from_idata_fdata(
            model, np.array(data, dtype='int32'), np.array(data, dtype='float32'))
        assert len(rbes) == 2, rbes
        rbes_expected = [
            ['RBE3', '61', '71', '123456', '1.', '123', '70', '75',
            '77'],
            ['RBE3', '62', '71', '123456', '1.', '123', '58', '59',
            '72'],
        ]
        for rbe, rbe_expected in zip(rbes, rbes_expected):
            msg = 'rbe:\n%s\nexpected:\n%s' % (rbe.rstrip().split(), rbe_expected)
            assert rbe.rstrip().split() == rbe_expected, msg

        #--------------------------------------------------------------
        data = [11, 10, 123456, 1.0, 456, 111, 100, -1, -2, -3,
                25, 50, 123456, 1.0, 456, 125, 103, -1, -2, -3]
        rbes = read_rbe3s_from_idata_fdata(
            model, np.array(data, dtype='int32'), np.array(data, dtype='float32'))
        assert len(rbes) == 2, rbes

        rbes_expected = [
            'RBE3          11              10  123456      1.     456     111     100',
            'RBE3          25              50  123456      1.     456     125     103',
        ]
        for rbe, rbe_expected in zip(rbes, rbes_expected):
            assert rbe.rstrip() == rbe_expected, rbe
        #--------------------------------------------------------------
        data = [407, 4, 123, 1.0, 123, 41201, 41210, 41212, 41221, -1,
                -0.25, 123, 41200, 41202, 41220, 41222, -1, -2, -3, 0.0,
                408, 4, 456, 1.0, 123, 41201, 41210, 41212, 41221, -1,
                1.0, 123, 41200, 41202, 41220, 41222, -1, -2, -3, 0.0]
        rbes = read_rbe3s_from_idata_fdata(
            model, np.array(data, dtype='int32'), np.array(data, dtype='float32'))
        assert len(rbes) == 2, rbes

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
