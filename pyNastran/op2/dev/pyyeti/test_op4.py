import numpy as np
import op4
from scipy.io import matlab
import os
from glob import glob
#from nose.tools import *

import unittest
from op2 import OP2

class TestOP4New(unittest.TestCase):

    def test_rdop4(self):
        matfile = 'nastran_op4_data/r_c_rc.mat'
        filenames = glob('nastran_op4_data/*.op4')
        o4 = op4.OP4()
        m = matlab.loadmat(matfile)
        for filename in filenames:
            if filename.find('badname') > -1:
                with assert_warns(RuntimeWarning) as cm:
                    dct = o4.dctload(filename)
                the_warning = str(cm.warning)
                assert 0 == the_warning.find('Output4 file has matrix '
                                             'name: 1mat')
                with assert_warns(RuntimeWarning) as cm:
                    names, mats, forms, mtypes = o4.listload(filename)
                the_warning = str(cm.warning)
                assert 0 == the_warning.find('Output4 file has matrix '
                                             'name: 1mat')
                with assert_warns(RuntimeWarning) as cm:
                    names2, sizes, forms2, mtypes2 = o4.dir(filename,
                                                            verbose=False)
                the_warning = str(cm.warning)
                assert 0 == the_warning.find('Output4 file has matrix '
                                             'name: 1mat')
            else:
                dct = o4.dctload(filename)
                names, mats, forms, mtypes = o4.listload(filename)
                names2, sizes, forms2, mtypes2 = o4.dir(filename,
                                                        verbose=False)
            assert sorted(dct.keys()) == sorted(names)
            assert names == names2
            assert forms == forms2
            assert mtypes == mtypes2
            for mat, sz in zip(mats, sizes):
                assert mat.shape == sz
            for nm in dct:
                if nm[-1] == 's':
                    matnm = nm[:-1]
                elif nm == '_1mat':
                    matnm = 'rmat'
                else:
                    matnm = nm
                assert np.allclose(m[matnm], dct[nm][0])
                pos = names.index(nm)
                assert np.allclose(m[matnm], mats[pos])
                assert dct[nm][1] == forms[pos]
                assert dct[nm][2] == mtypes[pos]
            nm2 = nm = 'rcmat'
            if filename.find('single') > -1:
                nm2 = 'rcmats'
            if filename.find('badname') > -1:
                with assert_warns(RuntimeWarning) as cm:
                    dct = o4.dctload(filename, nm2)
                    name, mat = o4.listload(filename, [nm2])[1:]
            else:
                dct = o4.dctload(filename, [nm2])
                name, mat = o4.listload(filename, nm2)[1:]
            assert np.allclose(m[nm], dct[nm2][0])
            assert np.allclose(m[nm], mat[0])


    def test_rdop4_zero_rowscutoff(self):
        matfile = 'nastran_op4_data/r_c_rc.mat'
        filenames = glob('nastran_op4_data/*.op4')
        o4 = op4.OP4()
        o4._rowsCutoff = 0
        m = matlab.loadmat(matfile)
        for filename in filenames:
            if filename.find('badname') > -1:
                with assert_warns(RuntimeWarning) as cm:
                    dct = o4.dctload(filename)
                the_warning = str(cm.warning)
                assert 0 == the_warning.find('Output4 file has matrix '
                                             'name: 1mat')
                with assert_warns(RuntimeWarning) as cm:
                    names, mats, forms, mtypes = o4.listload(filename)
                the_warning = str(cm.warning)
                assert 0 == the_warning.find('Output4 file has matrix '
                                             'name: 1mat')
                with assert_warns(RuntimeWarning) as cm:
                    names2, sizes, forms2, mtypes2 = o4.dir(filename,
                                                            verbose=False)
                the_warning = str(cm.warning)
                assert 0 == the_warning.find('Output4 file has matrix '
                                             'name: 1mat')
            else:
                dct = o4.dctload(filename)
                names, mats, forms, mtypes = o4.listload(filename)
                names2, sizes, forms2, mtypes2 = o4.dir(filename,
                                                        verbose=False)
            assert sorted(dct.keys()) == sorted(names)
            assert names == names2
            assert forms == forms2
            assert mtypes == mtypes2
            for mat, sz in zip(mats, sizes):
                assert mat.shape == sz
            for nm in dct:
                if nm[-1] == 's':
                    matnm = nm[:-1]
                elif nm == '_1mat':
                    matnm = 'rmat'
                else:
                    matnm = nm
                assert np.allclose(m[matnm], dct[nm][0])
                pos = names.index(nm)
                assert np.allclose(m[matnm], mats[pos])
                assert dct[nm][1] == forms[pos]
                assert dct[nm][2] == mtypes[pos]

            nm2 = nm = 'rcmat'
            if filename.find('single') > -1:
                nm2 = 'rcmats'
            if filename.find('badname') > -1:
                with assert_warns(RuntimeWarning) as cm:
                    dct = o4.dctload(filename, nm2)
                    name, mat = o4.listload(filename, [nm2])[1:]
            else:
                dct = o4.dctload(filename, [nm2])
                name, mat = o4.listload(filename, nm2)[1:]
            assert np.allclose(m[nm], dct[nm2][0])
            assert np.allclose(m[nm], mat[0])


    def test_rdop4_partb(self):
        filenames = glob('nastran_op4_data/*other')
        file1 = filenames[0]
        filenames = filenames[1:]
        o4 = op4.OP4()
        dct = o4.dctload(file1)
        for filename in filenames:
            dct2 = o4.dctload(filename)
            assert set(dct2.keys()) == set(dct.keys())
            for nm in dct:
                for j in range(3):
                    assert np.allclose(dct2[nm][j], dct[nm][j])


    def test_wtop4(self):
        matfile = 'nastran_op4_data/r_c_rc.mat'
        o4 = op4.OP4()
        m = matlab.loadmat(matfile)
        names = ['rmat', 'cmat', 'rcmat']
        mats = []
        wtdct = {}
        for nm in names:
            mats.append(m[nm])
            wtdct[nm] = m[nm]
        # write(filename, names, matrices=None,
        #       binary=True, digits=16, endian='')
        filenames = [
            ['nastran_op4_data/temp_ascii.op4', False, ''],
            ['nastran_op4_data/temp_le.op4', True, '<'],
            ['nastran_op4_data/temp_be.op4', True, '>'],
            ]
        for item in filenames:
            filename = item[0]
            binary = item[1]
            endian = item[2]
            o4.write(filename, names, mats,
                     binary=binary, endian=endian)
            names2, sizes, forms, mtypes = o4.dir(filename, verbose=False)
            assert names2 == names
            dct = o4.dctload(filename)
            for nm in dct:
                assert np.allclose(m[nm], dct[nm][0])
            o4.write(filename, wtdct,
                     binary=binary, endian=endian)
            dct = o4.dctload(filename)
            for nm in dct:
                assert np.allclose(m[nm], dct[nm][0])
        # clean up:
        for item in filenames:
            os.remove(item[0])


    def test_wtop4_single(self):
        matfile = 'nastran_op4_data/r_c_rc.mat'
        o4 = op4.OP4()
        m = matlab.loadmat(matfile)
        name = 'rmat'
        mat = m[name]
        # write(filename, names, matrices=None,
        #       binary=True, digits=16, endian='')
        filenames = [
            ['nastran_op4_data/temp_ascii.op4', False, ''],
            ]
        for item in filenames:
            filename = item[0]
            binary = item[1]
            endian = item[2]
            o4.write(filename, name, mat,
                     binary=binary, endian=endian)
            dct = o4.dctload(filename)
            for nm in dct:
                assert np.allclose(m[nm], dct[nm][0])
        # clean up:
        for item in filenames:
            os.remove(item[0])


    def test_wtop4_nonbigmat_binary(self):
        filenames = glob('nastran_op4_data/*.op4') +\
                    glob('nastran_op4_data/*.op4.other')
        o4 = op4.OP4()
        for name in filenames:
            if name.find('badname') != -1:
                continue
            data = o4.listload(name)
            o4.write('temp.op4', data[0], data[1], sparse='nonbigmat')
            data2 = o4.listload('temp.op4')
            assert data[0] == data2[0]
            for d1, d2 in zip(data[1], data2[1]):
                assert np.all(d1 == d2)
            os.remove('temp.op4')


    def test_wtop4_bigmat_binary(self):
        filenames = glob('nastran_op4_data/*.op4') +\
                    glob('nastran_op4_data/*.op4.other')
        o4 = op4.OP4()
        for name in filenames:
            if name.find('badname') != -1:
                continue
            data = o4.listload(name)
            o4.write('temp.op4', data[0], data[1], sparse='bigmat')
            data2 = o4.listload('temp.op4')
            assert data[0] == data2[0]
            for d1, d2 in zip(data[1], data2[1]):
                assert np.all(d1 == d2)
            os.remove('temp.op4')


    def test_wtop4_nonbigmat_ascii(self):
        filenames = glob('nastran_op4_data/*.op4') +\
                    glob('nastran_op4_data/*.op4.other')
        o4 = op4.OP4()
        for name in filenames:
            if name.find('badname') != -1:
                continue
            data = o4.listload(name)
            o4.write('temp.op4', data[0], data[1], sparse='nonbigmat',
                     binary=False)
            data2 = o4.listload('temp.op4')
            assert data[0] == data2[0]
            for d1, d2 in zip(data[1], data2[1]):
                assert np.all(d1 == d2)
            os.remove('temp.op4')


    def test_wtop4_bigmat_ascii(self):
        filenames = glob('nastran_op4_data/*.op4') +\
                    glob('nastran_op4_data/*.op4.other')
        o4 = op4.OP4()
        for name in filenames:
            if name.find('badname') != -1:
                continue
            data = o4.listload(name)
            o4.write('temp.op4', data[0], data[1], sparse='bigmat',
                     binary=False)
            data2 = o4.listload('temp.op4')
            assert data[0] == data2[0]
            for d1, d2 in zip(data[1], data2[1]):
                assert np.all(d1 == d2)
            os.remove('temp.op4')


    def test_non_float64(self):
        i8 = np.array([1, 2, 0, 4], np.int8)
        i16 = i8.astype(np.int16)
        i32 = i8.astype(np.int32)
        i64 = i8.astype(np.int64)
        f32 = i8.astype(np.float32)
        c32 = (f32 + 1j*f32).astype(np.complex64)
        o4 = op4.OP4()
        for mat in [i8, i16, i32, i64, f32, c32]:
            o4.write('temp.op4', dict(mat=mat))
            mat2 = o4.dctload('temp.op4', 'mat')['mat'][0]
            assert np.all(mat2 == mat)
            os.remove('temp.op4')


if __name__ == '__main__':
    unittest.main()
