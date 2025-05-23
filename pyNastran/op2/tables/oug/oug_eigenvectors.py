import numpy as np
from pyNastran.op2.result_objects.table_object import RealTableArray, ComplexTableArray
from pyNastran.f06.f06_formatting import write_floats_13e


class ComplexEigenvectorArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp: str='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []

        words = []
        #freq = self.eigrs[i]
        #freq = 0.0
        #words.append('%16s = %12E\n' % ('EIGENVALUE', freq))
        #words.append('%%16s = %%12E\n')
        has_cycle = False
        if hasattr(self, 'cycles'):
            has_cycle = True
        #if has_cycle:
        #words.append('%s = %s          C O M P L E X   E I G E N V E C T O R   N O. %s\n' % ('%16s', '%12E', '%10i'))
            #msg.append('%16s = %12E          C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' % ('CYCLES', self.mode_cycle, imode))
        #else:
        words.append('      COMPLEX EIGENVALUE = %s, %s\n' % ('%12E', '%12E'))
        words.append('                                       C O M P L E X   E I G E N V E C T O R   NO. %s\n' % '%10i')
        #msg.append('                                         C O M P L E X   E I G E N V E C T O R   N O . %10i\n \n' % (imode))

        #msg.append('      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n')
        #words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, is_mag_phase, is_sort1)


class RealEigenvectorArray(RealTableArray):
    def __init__(self, data_code, is_sort1: bool,
                 isubcase: int, dt, f06_flag: bool=False):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def scale(self, factors: np.ndarray) -> None:
        """
        generic scaling of a mode by a series of factors.
        Can generate a series of modes with a single call.
        """
        assert factors.dtype.name in ['complex64', 'complex128', 'float64', 'float32'], factors.dtype.name
        assert factors.ndim == 2, factors.shape
        nold_modes, nnew_modes = factors.shape
        self.modes = np.array(nnew_modes, dtype='int32')

        nold_modes_og, nnodes, six = self.data.shape
        assert nold_modes == nold_modes_og, (nold_modes, nold_modes_og)

        # checking the einsum
        data2 = np.zeros((nnew_modes, nnodes, 6), dtype=factors.dtype)
        for nmode in range(nnew_modes):
            for omode in range(nold_modes):
                tempi = factors[omode, nmode] * self.data[omode, :, :]
                data2[nmode, :, :] += tempi

        # new = factors * modes
        data = np.einsum('ij,ilm->jlm', factors, self.data)
        assert data.shape == (nnew_modes, nnodes, 6), data.shape
        assert np.allclose(data, data2)
        self.data = data

    def get_phi(self) -> np.ndarray:
        """
        gets the eigenvector matrix

        Returns
        -------
        phi : (ndof, nmodes)
            the eigenvector matrix

        TODO: doesn't consider SPOINTs/EPOINTs
        """
        nmodes, nnodes = self.data.shape[:2]
        ndof = nnodes * 6
        phi_transpose = self.data.reshape(nmodes, ndof)
        return phi_transpose.T

    @classmethod
    def phi_to_data(cls, phi: np.ndarray) -> np.ndarray:
        """(ndof, nmodes) -> (nmodes, nnodes, 6)"""
        ndof, nmodes = phi.shape
        nnodes = ndof // 6
        assert ndof % 6 == 0
        phi2 = phi.T  # nmodes, ndof
        data = phi2.reshape(nmodes, nnodes, 6)
        return data

    def set_phi(self, phi: np.ndarray) -> None:
        """(ndof, nmodes) -> (nmodes, nnodes, 6)"""
        self.data = self.phi_to_data(phi)

    def mac(self, obj) -> np.ndarray:
        """
        Parameters
        ----------
        obj : RealEigenvectorArray
            defines the phi2 matrix

        Returns
        -------
        mac_matrix: (nmodes1, nnodes2) float array
            the MAC matrix

        It is important to note that the mode shapes should be normalized
        before calculating the MAC. Also, the number of degrees of
        freedom (rows) in both mode shape matrices must be the same.

        """
        phi1 = self.get_phi()
        phi2 = self.get_phi()
        nmodes1, ndof1 = phi1.shape
        nmodes2, ndof2 = phi1.shape
        assert ndof1 == ndof2, (ndof1, ndof2)
        mac_matrix = np.zeros((nmodes1, nmodes2), dtype='float64')
        for i in range(nmodes1):
            for j in range(nmodes2):
                #mac_matrix[i, j]= np.abs(phia[:, i].T @ phib[:, j])**2 / ((phia[:, i].T @ phia[:, i]) * (phib[:, j].T @ phib[:, j]))
                mac_matrix[i, j] = np.abs(phi1[:, i].T @ phi2[:, j])**2 / ((phi1[:, i].T @ phi1[:, i]) * (phi2[:, j].T @ phi2[:, j]))
        return mac_matrix

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        #if self.nonlinear_factor not in (None, np.nan):
            #return self._write_f06_transient(header, page_stamp, page_num, f06_file,
                                             #is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        # modes get added
        words = '                                         R E A L   E I G E N V E C T O R   N O . %10i\n \n' \
                '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n'

        #if not len(header) >= 3:
            #header.append('')
        for itime in range(self.ntimes):
            node = self.node_gridtype[:, 0]
            gridtype = self.node_gridtype[:, 1]
            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            dt = self._times[itime]
            #if isinstance(dt, float):
                #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            #else:
                #header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            f06_file.write(''.join(header + [words % dt]))
            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
                sgridtype = self.recast_gridtype_as_string(gridtypei)
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                vals2 = write_floats_13e(vals)
                (dx, dy, dz, rx, ry, rz) = vals2
                f06_file.write(
                    '%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                        node_id, sgridtype, dx, dy, dz, rx, ry, rz))

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1
