from struct import pack
import inspect

import numpy as np
from numpy import zeros, allclose

from pyNastran.op2.result_objects.op2_objects import get_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import _eigenvalue_header, write_floats_13e #, get_key0


class RealFastArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]
        self.nelements = 0  # result specific

        #if is_sort1:
            #self.add_new_eid = self.add_new_eid_sort1
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        raise NotImplementedError()

    def build(self):
        """sizes the vectorized attributes of the RealShearArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        _times = zeros(self.ntimes, dtype=self.analysis_fmt)
        element = zeros(self.nelements, dtype=idtype)

        # [force_x, force_y, force_z, moment_x, moment_y, moment_z]
        data = zeros((self.ntimes, self.ntotal, 6), dtype=fdtype)

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element = group.create_dataset('element', data=element)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element = element
            self.data = data

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            #Mode                            1             2             3
            #Freq                 1.482246e-10  3.353940e-09  1.482246e-10
            #Eigenvalue          -8.673617e-19  4.440892e-16  8.673617e-19
            #Radians              9.313226e-10  2.107342e-08  9.313226e-10
            #ElementID Item
            #22        max_shear  8.050749e-13  5.871460e-07  2.035239e-12
            #         avg_shear -8.050749e-13  5.871460e-07  2.035239e-12
            #         margin     1.401298e-45  1.401298e-45  1.401298e-45
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            #Static     axial           SMa  torsion           SMt
            #ElementID
            #14           0.0  1.401298e-45      0.0  1.401298e-45
            #15           0.0  1.401298e-45      0.0  1.401298e-45
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (max_shear, avg_shear, margin) = t1
                        (max_shear2, avg_shear2, margin2) = t2
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s)\n  (%s, %s, %s)\n' % (
                                eid,
                                max_shear, avg_shear, margin,
                                max_shear2, avg_shear2, margin2)
                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, force_x, force_y, force_z, moment_x, moment_y, moment_z):
        """
                                  S T R E S S E S   I N   F A S T E N E R   E L E M E N T S   ( C F A S T )
        ELEMENT-ID         FORCE-X          FORCE-Y          FORCE-Z         MOMENT-X         MOMENT-Y         MOMENT-Z
               279      -1.485744E-01    -3.137333E-01    -6.343584E-01    -9.968021E-03     7.256226E-01     7.248363E-02

        """
        assert self.sort_method == 1, self
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [force_x, force_y, force_z, moment_x, moment_y, moment_z]
        self.ielement += 1

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append(f'  element.shape = {self.element.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append(f'  element type: {self.element_name}-{self.element_type}\n')
        msg += self.get_data_code()
        return msg

    def get_f06_header(self) -> list[str]:
        raise NotImplementedError('CWELD...')

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header()

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            force_x = self.data[itime, :, 0]
            force_y = self.data[itime, :, 1]
            force_z = self.data[itime, :, 2]
            moment_x = self.data[itime, :, 3]
            moment_y = self.data[itime, :, 4]
            moment_z = self.data[itime, :, 5]

            for eid, force_xi, force_yi, force_zi, moment_xi, moment_yi, moment_zi in zip(eids, force_x, force_y, force_z, moment_x, moment_y, moment_z):
                [force_xi, force_yi, force_zi, moment_xi, moment_yi, moment_zi] = write_floats_13e(
                    [force_xi, force_yi, force_zi, moment_xi, moment_yi, moment_zi])
                f06_file.write('      %8i   %13s  %13s  %13s  %13s  %13s  %s\n' % (
                    eid, force_xi, force_yi, force_zi, moment_xi, moment_yi, moment_zi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2_file, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian='>'):
        """writes an OP2"""
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        unused_eids = self.element

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        op2_ascii.write(f'  ntimes = {self.ntimes}\n')

        #eids = self.element
        eids_device = self.element * 10 + self.device_code

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if not self.is_sort1:
            raise NotImplementedError('SORT2')
        #struct1 = Struct(endian + b'i 6f')

        fdtype = self.data.dtype
        if self.size == fdtype.itemsize:
            pass
        else:
            # print(f'downcasting {self.class_name}...')
            #idtype = np.int32(1)
            fdtype = np.float32(1.0)

        # [eid, force_x, force_y, force_z, moment_x, moment_y, moment_z]
        data_out = np.empty((nelements, 7), dtype=fdtype)
        data_out[:, 0] = eids_device.view(fdtype)

        op2_ascii.write(f'nelements={nelements:d}\n')

        for itime in range(self.ntimes):
            #print('3, %s' % itable)
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            #print('stress itable = %s' % itable)
            itable -= 1
            #print('4, %s' % itable)
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')

            # [eid, max_shear, avg_shear, margin]
            data_out[:, 1:] = self.data[itime, :, :]
            op2_file.write(data_out)

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealFastStressArray(RealFastArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealFastArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['force_x', 'force_y', 'force_z', 'moment_x', 'moment_y', 'moment_z']
        return headers

    def get_f06_header(self) -> list[str]:
        msg = [
            '                           S T R E S S E S   I N   F A S T E N E R   E L E M E N T S   ( C F A S T )\n'
            ' \n'
            ' ELEMENT-ID         FORCE-X          FORCE-Y          FORCE-Z         MOMENT-X         MOMENT-Y         MOMENT-Z\n'
            #'        279      -1.485744E-01    -3.137333E-01    -6.343584E-01    -9.968021E-03     7.256226E-01     7.248363E-02\n'
        ]
        return msg

class RealFastStrainArray(RealFastArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealFastArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['disp_x', 'disp_y', 'disp_z', 'rotation_x', 'rotation_y', 'rotation_z']
        return headers

    def get_f06_header(self) -> list[str]:
        msg = [
            '                     S T R A I N (D I S P)   I N   F A S T E N E R   E L E M E N T S   ( C F A S T )\n'
            ' \n'
            ' ELEMENT-ID         DISP-X           DISP-Y           DISP-Z         ROTATION-X       ROTATION-Y       ROTATION-Z\n'
            #'        100      -1.554312E-14     1.059252E+00    -7.944389E-01    -4.062902E-08     9.630819E+01     1.284109E+02\n'
        ]
        return msg
