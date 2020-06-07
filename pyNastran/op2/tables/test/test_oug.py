import os
import unittest
import numpy as np

from pyNastran.op2.tables.oug.oug_displacements import (
    RealDisplacementArray, ComplexDisplacementArray)

#from pyNastran.op2.tables.oug.oug_velocities import (
    #RealVelocityArray, ComplexVelocityArray)

#from pyNastran.op2.tables.oug.oug_accelerations import (
    #RealAccelerationArray, ComplexAccelerationArray)

#from pyNastran.op2.tables.oug.oug_temperatures import (
    #RealTemperatureArray)

from pyNastran.op2.tables.oug.oug_eigenvectors import (
    RealEigenvectorArray, ComplexEigenvectorArray,
)

#from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import RealThermalVelocityVectorArray


class TestOUG(unittest.TestCase):
    def test_real(self):
        """tests real displacements"""
        table_name = 'OUGV1'
        node_gridtype = np.zeros((10, 2), dtype='int32')
        node_gridtype[:, 0] = np.arange(1, 11)
        data = np.zeros((1, 10, 6), dtype='float32')
        isubcase = 1
        disp = RealDisplacementArray.add_static_case(
            table_name, node_gridtype, data, isubcase, is_sort1=True)
        disp.code_information()
        with open('disp.f06', 'w') as f06_file:
            disp.write_f06(f06_file, header=None, page_stamp='PAGE %s',
                           page_num=1, is_mag_phase=False, is_sort1=True)

        itable = 1
        new_result = False
        date = None
        with open('disp.op2', 'wb') as op2_file, open('disp.txt', 'w') as fascii:
            disp.write_op2(op2_file, fascii, itable, new_result,
                  date, is_mag_phase=False, endian=b'<')
        os.remove('disp.f06')
        os.remove('disp.op2')
        os.remove('disp.txt')

    def test_eigenvectors(self):
        """tests real displacements"""
        table_name = 'OUGV1'
        node_gridtype = np.zeros((10, 2), dtype='int32')
        node_gridtype[:, 0] = np.arange(1, 11)
        isubcase = 1

        nmodes = 3
        data = np.zeros((nmodes, 10, 6), dtype='float32')

        modes = np.arange(nmodes, dtype='int32')
        eigenvalues = np.zeros(nmodes, dtype='float32')
        mode_cycles = np.zeros(nmodes, dtype='float32')
        eig = RealEigenvectorArray.add_modal_case(
            table_name, node_gridtype, data, isubcase,
            modes, eigenvalues, mode_cycles, is_sort1=True)
        eig.code_information()
        with open('eig.f06', 'w') as f06_file:
            eig.write_f06(f06_file, header=None, page_stamp='PAGE %s',
                           page_num=1, is_mag_phase=False, is_sort1=True)

        itable = 1
        new_result = False
        date = None
        with open('eig.op2', 'wb') as op2_file, open('eig.txt', 'w') as fascii:
            eig.write_op2(op2_file, fascii, itable, new_result,
                  date, is_mag_phase=False, endian=b'>')
        os.remove('eig.f06')
        os.remove('eig.op2')
        os.remove('eig.txt')

    def test_transient(self):
        """tests transient displacements"""
        table_name = 'OUGV1'
        node_gridtype = np.zeros((10, 2), dtype='int32')
        node_gridtype[:, 0] = np.arange(1, 11, dtype='int32')
        isubcase = 1

        ntimes = 3
        data = np.zeros((ntimes, 10, 6), dtype='float32')

        times = np.linspace(0., 1., num=ntimes, dtype='float32')
        tdisp = RealDisplacementArray.add_transient_case(
            table_name, node_gridtype, data, isubcase,
            times, is_sort1=True)
        tdisp.code_information()
        header = []
        with open('tdisp.f06', 'w') as f06_file:
            tdisp.write_f06(f06_file, header=header, page_stamp='PAGE %s',
                           page_num=1, is_mag_phase=False, is_sort1=True)

        itable = 1
        new_result = False
        date = None
        with open('tdisp.op2', 'wb') as op2_file, open('tdisp.txt', 'w') as fascii:
            tdisp.write_op2(op2_file, fascii, itable, new_result,
                  date, is_mag_phase=False, endian=b'<')
        os.remove('tdisp.f06')
        os.remove('tdisp.op2')
        os.remove('tdisp.txt')

    def test_freq(self):
        """tests complex displacements"""
        table_name = 'OUGV1'
        node_gridtype = np.zeros((10, 2), dtype='int32')
        node_gridtype[:, 0] = np.arange(1, 11)
        isubcase = 1

        nfreqs = 3
        data = np.zeros((nfreqs, 10, 6), dtype='float32')

        freqs = np.zeros(nfreqs, dtype='float32')
        cdisp = ComplexDisplacementArray.add_freq_case(
            table_name, node_gridtype, data, isubcase,
            freqs, is_sort1=True)
        cdisp.code_information()
        header = ['', '', '',]
        with open('cdisp.f06', 'w') as f06_file:
            cdisp.write_f06(f06_file, header=header, page_stamp='PAGE %s',
                           page_num=1, is_mag_phase=False, is_sort1=True)

        itable = 1
        new_result = False
        date = None
        with open('cdisp.op2', 'wb') as op2_file, open('cdisp.txt', 'w') as fascii:
            cdisp.write_op2(op2_file, fascii, itable, new_result,
                  date, is_mag_phase=False, endian=b'<')
        os.remove('cdisp.f06')
        os.remove('cdisp.op2')
        os.remove('cdisp.txt')

if __name__ == '__main__':   # pragma: no cover
    unittest.main()
