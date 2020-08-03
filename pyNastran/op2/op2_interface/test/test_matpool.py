import unittest
import numpy as np
from pyNastran.op2.op2_interface.utils_matpool import find_all_dmigs_start_stop


class TestMatpool(unittest.TestCase):

    def test_matpool_1(self):
        ints = np.array([
            #0*         1          2  3      4   5    6  7  8  9      10*  11     12  13          14          15  16  17   18
            #name                  ?  shape  tin tout ?  ?  ?  colnid dof  rownid dof value
            538980418,  538976288, 0, 6,     2,  2,   0, 0, 0, 3,     1,   4,     1,  1717986918, 1072064102, -1, -1, -1, -1,
            #
            #19*        20         21 22     23  24   25 26 27 28     29   30     31 32           33          34  35  36  37
            540029250,  538976288, 0, 6,     1,  1,   0, 0, 0, 10402, 1,   10402, 3, -1717986918, 1071225241, -1, -1, -1, -1,
            #38*        39         40 41     42  43   44 45 46 47     47   49     50 51           52          53  54  55  56
            540029762,  538976288, 0, 6,     1,  1,   0, 0, 0, 30301, 1,   30301, 3, 0,           1071644672, -1, -1, -1, -1],
            dtype='int32')
        data = ints.tobytes()
        iminus1 = np.where(ints == -1)[0]
        #print(iminus1)
        #iminus1 = [15, 16, 17, 18,
                   #34, 35, 36, 37,
                   #53, 54, 55, 56]
        size = 4
        header_fmt = b'<8s 7i'
        istarts, istops, outs, kstarts, kstops = find_all_dmigs_start_stop(
            data, header_fmt, size, iminus1, debug=False)
        # good
        assert np.array_equal(istarts, [0, 19, 38])
        assert np.array_equal(istops, [15, 34, 53])
        #print(kstarts, kstops)
        assert np.array_equal(kstarts, [[9], [28], [47]]), kstarts
        assert np.array_equal(kstops, [[15], [34], [53]]), kstops
        #print(istarts, istops)

    def test_matpool_2(self):
        ints = np.array([
            #0*         1          2  3      4   5    6  7  8  9      10*  11     12  13          14          15  16  17   18
            #name                  ?  shape  tin tout ?  ?  ?  colnid dof  rownid dof value
            538980418,  538976288, 0, 6,     2,  2,   0, 0, 0, 3,     1,   4,     1,  1717986918, 1072064102, -1, -1, -1, -1],
            dtype='int32')
        data = ints.tobytes()
        iminus1 = np.where(ints == -1)[0]
        #print(iminus1)
        #iminus1 = [15, 16, 17, 18]
        size = 4
        header_fmt = b'<8s 7i'
        istarts, istops, outs, kstarts, kstops = find_all_dmigs_start_stop(
            data, header_fmt, size, iminus1, debug=False)
        # good
        assert np.array_equal(istarts, [0]), istarts
        assert np.array_equal(istops, [15]), istops
        #print(kstarts, kstops)
        assert np.array_equal(kstarts, [[9]]), kstarts
        assert np.array_equal(kstops, [[15]]), kstops
        #print(istarts, istops)

    def _test_matpool_64(self):
        """per C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\mnfexam32_0.op2"""
        ints = np.array([
            # 0                  1                    2  3  4  5  6  7  8
            2314885531441381707, 2314885530818453536, 0, 6, 1, 1, 0, 0, 1958,
            #col dof row dof
            #9   10 11   12 13  14  15
            200, 2, 200, 2, -1, -1, -1,
            #16
            200, 3, 200, 2, -1, 200, 3, 0, -1, -1,
            201, 1, 200, 2, -1, 200, 3, 0, 201, 1, 0, -1, -1,
            201, 2, 200, 2, -1, 200, 3, 0, 201, 1, 0, 201, 2, 0, -1, -1,
            202, 1, 200, 2, -1, 200, 3, 0, 201, 1, 0, 201, 2, 0, 202, 1, 0, -1, -1,
            207, 1, 207, 1, 4735090208478144015, -1, -1,
            207, 2, 207, 2, 4736007071530354485, -1, -1,
            207, 3, 207, 3, 4748537773998944216, -1, -1,
            207, 4, 207, 4, 4748798681640498386, -1, -1,
            207, 5, 207, 5, 4756815377936266833, -1, -1,
            207, 6, 207, 6, 4756889397520511640, -1, -1,
            243, 1, 243, 1, 4760613669112060268, -1, -1,
            243, 2, 243, 2, 4762603449244146308, -1, -1,
            243, 3, 243, 3, 4773200220112711257, -1, -1,
            243, 4, 243, 4, 4774286936441837097, -1, -1,
            243, 5, 243, 5, 4775661555386971815, -1, -1,
            243, 6, 243, 6, 4779329891370826714, -1, -1,
            100001, 1, 100001, 0, 4780310554305573060, -1, -1,
            100002, 1, 100002, 0, 4782399563524630904, -1, -1,
            100003, 1, 100003, 0, 4789460462940233755, -1, -1,
            100004, 1, 100004, 0, 4804544589863548165, -1, -1, -1, -1],
            dtype='int64')
        data = ints.tobytes()

        floats = np.frombuffer(data, dtype='float32')
        doubles = np.frombuffer(data, dtype='float64')
        #print(floats)
        #print(doubles)

        iminus1 = np.where(ints == -1)[0]
        #print(iminus1, print(len(ints)))

        header_fmt = b'<16s 7q'
        size = 8
        istarts, istops, outs, kstarts, kstops = find_all_dmigs_start_stop(
            data, header_fmt, size, iminus1, debug=False)
        # good
        assert np.array_equal(istarts, [0]), istarts
        assert np.array_equal(istops, [184]), istops
        #print(kstarts, kstops)
        assert np.array_equal(kstarts, [[9,  17, 24, 28, 34, 41, 47, 57, 63, 76, 83, 90, 97,  104, 111, 118, 125, 132, 139, 146, 153, 160, 167, 174, 181]]), kstarts
        assert np.array_equal(kstops,  [[13, 20, 24, 30, 37, 43, 53, 59, 72, 79, 86, 93, 100, 107, 114, 121, 128, 135, 142, 149, 156, 163, 170, 177, 184]]), kstops

if __name__ == '__main__':
    unittest.main()
