#pylint: disable=R0201
import struct
import unittest
import numpy as np
from pyNastran.op2.op2_geom import OP2Geom
from pyNastran.op2.tables.geom.geom4 import read_rbe3s_from_idata_fdata, ints_to_secset1s

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

    def test_group(self):
        """reads a GROUP card"""
        op2 = OP2Geom(make_geom=True, debug=False, log=None, debug_file=None, mode='msc')
        op2.op2_reader.factor = 1
        op2.idtype = 'int32'
        op2.idtype8 = 'int32'
        op2._uendian = '<'
        data = (17400, 174, 616,
                55, 0,
                -5, 90011, -1,
                -1,
                65, 0,
                -5, 90012, -1,
                -1,
                75, 0,
                -5, 90013, -1,
                -1)
        data_bytes1 = struct.pack(f'{len(data)}i', *data)
        op2.reader_edt._read_group(data_bytes1, 12)

        data2 = (
            17400, 174, 616,
            111, 5, b'THIS IS GROUP 111   ',
            -2, 5, b'THIS IS METADATA', -1,
            -5, 1, 0, 10, -1,
            -1)
        data_bytes2 = struct.pack(b'5i 20s 2i 16s 7i', *data2)
        op2.reader_edt._read_group(data_bytes2, 12)

        data3 = (17400, 174, 616,
                 55, 0,
                 -5, 90011, -1,
                 -1,
                 65, 0,
                 -5, 90012, -1,
                 -1,
                 75, 0,
                 -5, 90013, -1,
                 -1)
        data_bytes3 = struct.pack(f'{len(data)}i', *data3)
        op2.reader_edt._read_group(data_bytes3, 12)
        #print('done...')

    def test_uset1(self):
        """tests USET1"""
        op2 = OP2Geom(make_geom=True, debug=False, log=None, debug_file=None, mode='msc')
        op2.op2_reader.factor = 1
        op2.idtype = 'int32'
        op2.idtype8 = 'int32'
        op2.fdtype8 = 'float32'
        op2._uendian = '<'
        data = (2110, 21, 194,
                2.0, 123456, 0, 44, 45, 48, 49, -1)
        data_bytes = struct.pack(b'3i f 7i', *data)
        op2.reader_geom4._read_uset1(data_bytes, 12)


        data_bytes = struct.pack(b'3i f 7i', *data)
        op2.reader_geom4._read_uset1(data_bytes, 12)

    def test_nsml1_nx(self):
        """tests NSML1-NX"""
        op2 = OP2Geom(make_geom=True, debug=False, log=None, debug_file=None, mode='msc')
        op2.idtype8 = 'int32'
        op2.fdtype8 = 'float32'
        op2._uendian = '<'

        fmt1 = (b'<3i ' +
                b'i 8s f ' + b'2i8s2i ' * 4 + b' i ' +
                b'i 8s f ' + b'2i8s2i ' * 7 + b' i')
        data = (
            3701, 37, 995,
            1, b'ELEMENT ', 466.2,
            3, 249311, b'THRU    ', 250189, -1,
            3, 250656, b'THRU    ', 251905, -1,
            3, 270705, b'THRU    ', 275998, -1,
            3, 332687, b'THRU    ', 334734, -1,
            -2,

            2, b'ELEMENT ', 77.7,
            3, 225740, b'THRU    ', 227065, -1,
            3, 227602, b'THRU    ', 228898, -1,
            3, 229435, b'THRU    ', 230743, -1,
            3, 231280, b'THRU    ', 233789, -1,

            3, 233922, b'THRU    ', 235132, -1,
            3, 235265, b'THRU    ', 236463, -1,
            3, 338071, b'THRU    ', 341134, -1,
            -2)
        data_bytes = struct.pack(fmt1, *data)
        op2.reader_ept._read_nsml1_nx(data_bytes, 12)

    def test_seelt(self):
        op2 = OP2Geom(make_geom=True, debug=False, log=None, debug_file=None, mode='msc')
        op2.idtype8 = 'int32'
        #op2.fdtype8 = 'float32'
        #op2._uendian = '<'

        data = (
            7902, 79, 302,
            5, 24, -1,
            66, 662001, 662003, 662007, 662019, -1,
            67, 672001, 672008, 672015, -1,
            68, 682001, 682019, -1)
        fmt1 = b'i' * len(data)
        data_bytes = struct.pack(fmt1, *data)
        op2.reader_geom1._read_seelt(data_bytes, 12)

    def test_seset(self):
        op2 = OP2Geom(make_geom=True, debug=False, log=None, debug_file=None, mode='msc')
        op2.idtype8 = 'int32'
        data = (
            5601, 56, 296,
            1, 33, 34, 37, 38, -1,
            1, 93, -98, -1,
            7, 1, -8, -1)
        fmt1 = b'i' * len(data)
        data_bytes = struct.pack(fmt1, *data)
        op2.reader_geom1._read_seset(data_bytes, 12)

    def test_snorm(self):
        op2 = OP2Geom(make_geom=True, debug=False, log=None, debug_file=None, mode='msc')
        #op2.idtype8 = 'int32'
        op2._endian = b'<'
        op2.op2_reader.factor = 1
        data = (5678, 71, 475,
                1059, 10, 0.0, 0.2, 0.0,
                -1059, 101000001, 0.0, 0.2, 0.0,
                )
        fmt1 = b'3i ' + b'2i3f ' * 2
        data_bytes = struct.pack(fmt1, *data)
        op2.reader_geom1._read_snorm(data_bytes, 12)

    def test_read_grid_8(self):
        op2 = OP2Geom(make_geom=True, debug=False, log=None, debug_file=None, mode='msc')
        op2._endian = b'<'
        op2.op2_reader.factor = 1
        data = (
            4501, 45, 1,
            # nid, cp, x1, x2, x3, cd, ps, seid
            1, 0, 1., 0., 0., 1, 0, 0,
            2, 0, 2., 0., 0., 1, 0, 0,
            3, 0, 3., 0., 0., 1, 0, 0,
        )
        fmt = b'<3i ' + b'2i 3f 3i ' * 3
        data_bytes3 = struct.pack(fmt, *data)
        assert len(data_bytes3) == 108, len(data_bytes3)

        op2.table_name = b'GEOM1'
        op2.reader_geom1._read_grid(data_bytes3, 12)

        op2.table_name = b'GEOM1N'
        with self.assertRaises(AssertionError):
            op2.reader_geom1._read_grid_11(data_bytes3, 12)

        #-----------------------------------
    def test_read_grid_11(self):
        op2 = OP2Geom(make_geom=True, debug=False, log=None, debug_file=None, mode='msc')
        op2._endian = b'<'
        op2.op2_reader.factor = 1
        data = (
            4501, 45, 1,
            # nid, cp, x1, x2, x3, cd, ps, seid
            1, 0, 1., 0., 0., 1, 0, 0,
            2, 0, 2., 0., 0., 1, 0, 0,
            3, 0, 3., 0., 0., 1, 0, 0,
        )
        fmt = b'<3i ' + b'2i 3d 3i ' * 3
        data_bytes = struct.pack(fmt, *data)
        assert len(data_bytes) == 144, len(data_bytes)

        op2.table_name = b'GEOM1'
        with self.assertRaises(AssertionError):
            op2.reader_geom1._read_grid(data_bytes, 12)

        op2.table_name = b'GEOM1N'
        op2.reader_geom1._read_grid(data_bytes, 12)

    def test_extrn(self):
        """reads an EXTRN"""
        op2 = OP2Geom(make_geom=True, debug=False, log=None, debug_file=None, mode='msc')
        op2.idtype8 = 'int32'
        #op2.fdtype8 = 'float32'
        #op2._uendian = '<'

        data = (
            1627, 16, 463,
            1, 123456, 2, 123456, 3, 123456, 4, 123456, 100001, 1, 100002, 1, 100003, 1, 100004, 1,
            -1, -1)
        fmt1 = b'i' * len(data)
        data_bytes = struct.pack(fmt1, *data)
        op2.reader_geom1._read_extrn(data_bytes, 12)

    def test_ints_to_secset1s(self):
        ints_to_secset1s('SECSET1', [61,  123456,    1, 610101, 610124])
        ints_to_secset1s('SECSET1', [100,    123,    0,     41,     42,  -1,
                                     200,    123,    0,     35,     36,  -1])
        ints_to_secset1s('SECSET1', [ 1, 123456,    1,   1001,   1006,
                                      1, 123456,    0,   1066,     -1])

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
