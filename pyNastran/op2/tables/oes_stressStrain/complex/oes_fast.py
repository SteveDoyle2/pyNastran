import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_complex_times_dtype
from pyNastran.op2.result_objects.utils_pandas import build_dataframe_transient_header, build_pandas_transient_elements
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_imag_floats_13e


class ComplexFastArray(OES_Object):
    """
    Common class for:
     - ComplexFastStressArray
     - ComplexFastStrainArray
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
        #self.element_node = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def _reset_indices(self) -> None:
        self.itotal = 0
        self.ielement = 0

    #def get_nnodes(self):
        #return get_nnodes(self)

    def build(self):
        """sizes the vectorized attributes of the ComplexShearArray"""
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.ntimes, self.nelements, self.ntotal, self.subtitle))
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        self.ntotal = self.nelements * nnodes# * 2
        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        idtype, cfdtype = get_complex_times_dtype(self.size)
        self._times = np.zeros(self.ntimes, dtype=self.analysis_fmt)
        #self.ntotal = self.nelements * nnodes

        self.element = np.zeros(self.nelements, dtype=idtype)

        # the number is messed up because of the offset for the element's properties
        if self.nelements != self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
                self.ntimes, self.nelements, nnodes, self.nelements * nnodes,
                self.ntotal)
            raise RuntimeError(msg)

        # stress: []
        # strain: [disp_x, disp_y, disp_z, rotation_x, rotation_y, rotation_z]
        self.data = np.zeros((self.ntimes, self.ntotal, 6), dtype=cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        #Mode                                          1                   2
        #EigenvalueReal                             -0.0                -0.0
        #EigenvalueImag                             -0.0                -0.0
        #Damping                                     0.0                 0.0
        #ElementID Item
        #22        max_shear  5.855954e-09+0.000000e+00j  0.000000+0.000000j
        #          avg_shear  5.855954e-09+0.000000e+00j  0.000000+0.000000j
        #import pandas as pd
        column_names, column_values = build_dataframe_transient_header(self)
        self.data_frame = build_pandas_transient_elements(
            self, column_values, column_names,
            self.headers, self.element, self.data)

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
                    for ieid, eid in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (tx1, ty1, unused_tz1, unused_rx1, unused_ry1, unused_rz1) = t1
                        (tx2, ty2, unused_tz2, unused_rx2, unused_ry2, unused_rz2) = t2
                        d = t1 - t2
                        if not np.allclose([tx1.real, tx1.imag, ty1.real, ty1.imag],
                                           [tx2.real, tx2.imag, ty2.real, ty2.imag], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s  (%s, %sj, %s, %sj)\n      (%s, %sj, %s, %sj)\n  dt12=(%s, %sj, %s, %sj)\n' % (
                                eid,
                                tx1.real, tx1.imag, ty1.real, ty1.imag,
                                tx2.real, tx2.imag, ty2.real, ty2.imag,
                                d[0].real, d[0].imag, d[1].real, d[1].imag,)
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
        """unvectorized method for adding SORT1 transient data"""
        assert self.sort_method == 1, self
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal] = [force_x, force_y, force_z, moment_x, moment_y, moment_z]
        self.element[self.itotal] = eid
        #self.ielement += 1
        self.itotal += 1

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.element.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i; table_name=%r\n' % (
                self.__class__.__name__, ntimes, nelements, nnodes, self.table_name))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i; table_name=%r\n' % (
                self.__class__.__name__, nelements, nnodes, self.table_name))
        msg.append('  data: [ntimes, nnodes, 2] where 2=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append(f'  element.shape = {self.element.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append('  %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        """
        '                      C O M P L E X  S T R E S S E S  I N  F A S T E N E R  E L E M E N T S   ( C F A S T )'
        '                                                          (REAL/IMAGINARY)'
        ' '
        ' ELEMENT-ID         FORCE-X          FORCE-Y          FORCE-Z         MOMENT-X         MOMENT-Y         MOMENT-Z'
        '        100       2.179551E+00     3.215275E+03    -2.411456E+03    -9.586007E-09     1.205728E+01     1.607637E+01'
        '                  1.893645E-03    -2.685570E+00     2.014177E+00     8.824217E-09    -1.007089E-02    -1.342785E-02'
        """
        #if header is None:
            #header = []
        #f.write(self.code_information())
        #return page_num
        msg_temp = _get_fast_msg(self.is_stress, is_mag_phase, is_sort1)
        ntimes = self.data.shape[0]
        eids = self.element
        if self.is_sort1:
            if is_sort1:
                for itime in range(ntimes):
                    dt = self._times[itime]

                    dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
                    header[1] = dt_line
                    msg = header + msg_temp
                    f06_file.write('\n'.join(msg))

                    force_x = self.data[itime, :, 0]
                    force_y = self.data[itime, :, 1]
                    force_z = self.data[itime, :, 2]
                    moment_x = self.data[itime, :, 3]
                    moment_y = self.data[itime, :, 4]
                    moment_z = self.data[itime, :, 5]
                    assert len(eids) == len(force_x)
                    assert len(force_x) > 0, force_x
                    for eid, force_xi, force_yi, force_zi, moment_xi, moment_yi, moment_zi in zip(eids, force_x, force_y, force_z, moment_x, moment_y, moment_z):
                        assert isinstance(eid, integer_types), 'eid=%s type=%s' % (eid, type(eid))
                        [rforce_x, iforce_x, rforce_y, iforce_y, rforce_z, iforce_z,
                         rmoment_x, imoment_x, rmoment_y, imoment_y, rmoment_z, imoment_z] = write_imag_floats_13e([force_xi, force_yi, force_zi, moment_xi, moment_yi, moment_zi], is_mag_phase)

                        #f.write('                      28                  0.0          /  0.0                           0.0          /  0.0\n')
                        f06_file.write(
                            '%24s  %13s  %13s  %13s  %13s  %13s  %s\n'
                            '%24s  %13s  %13s  %13s  %13s  %13s  %s\n' % (
                                eid, rforce_x, rforce_y, rforce_z, rmoment_x, rmoment_y, rmoment_z,
                                '', iforce_x, iforce_y, iforce_z, imoment_x, imoment_y, imoment_z))
                    f06_file.write(page_stamp % page_num)
                    page_num += 1
            else:
                # TODO: write in SORT2
                times = self._times
                raise NotImplementedError('Fast SORT2')
                #for ieid, eid in enumerate(eids):
                    #max_shear = self.data[:, ieid, 0].ravel()
                    #avg_shear = self.data[:, ieid, 1].ravel()
                    #for itime, max_sheari, avg_sheari in zip(times, max_shear, avg_shear):
                        #[rmax_shear, imax_shear, ravg_shear, iavg_shear
                         #] = write_imag_floats_13e([max_sheari, avg_sheari], is_mag_phase)

                        ##f06_file.write(
                            ##'   %6s   %-13s / %-13s     %-13s / %-13s\n' % (
                                ##eid, rmax_shear, imax_shear, ravg_shear, iavg_shear))
                        #f06_file.write(
                            #'%24s                 %-13s / %-13s                 %-13s / %-13s\n' % (
                                #eid, rmax_shear, imax_shear, ravg_shear, iavg_shear))
                    #f06_file.write(page_stamp % page_num)
                    #page_num += 1
        else:
            raise NotImplementedError('ComplexFastArray-sort2')
        return page_num - 1

    @property
    def headers(self) -> list[str]:
        return self._get_headers()

    def get_headers(self) -> list[str]:
        return self.headers

class ComplexFastStressArray(ComplexFastArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexFastArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self) -> list[str]:
        """
        '                      C O M P L E X  S T R E S S E S  I N  F A S T E N E R  E L E M E N T S   ( C F A S T )'
        '                                                          (REAL/IMAGINARY)'
        ' '
        ' ELEMENT-ID         FORCE-X          FORCE-Y          FORCE-Z         MOMENT-X         MOMENT-Y         MOMENT-Z'
        '        100       2.179551E+00     3.215275E+03    -2.411456E+03    -9.586007E-09     1.205728E+01     1.607637E+01'
        '                  1.893645E-03    -2.685570E+00     2.014177E+00     8.824217E-09    -1.007089E-02    -1.342785E-02'
        """
        return ['force_x', 'force_y', 'force_z', 'moment_x', 'moment_y', 'moment_z']

class ComplexFastStrainArray(ComplexFastArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexFastArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain, self.stress_bits

    def _get_headers(self):
        return ['disp_x', 'disp_y', 'disp_z', 'rotation_x', 'rotation_y', 'rotation_z']

def _get_fast_msg(is_stress: bool, is_mag_phase: bool, is_sort1: bool) -> list[str]:
    """
    '                      C O M P L E X  S T R E S S E S  I N  F A S T E N E R  E L E M E N T S   ( C F A S T )'
    '                                                          (REAL/IMAGINARY)'
    ' '
    ' ELEMENT-ID         FORCE-X          FORCE-Y          FORCE-Z         MOMENT-X         MOMENT-Y         MOMENT-Z'
    '        100       2.179551E+00     3.215275E+03    -2.411456E+03    -9.586007E-09     1.205728E+01     1.607637E+01'
    '                  1.893645E-03    -2.685570E+00     2.014177E+00     8.824217E-09    -1.007089E-02    -1.342785E-02'

    '              C O M P L E X  S T R A I N (D I S P)  I N  F A S T E N E R  E L E M E N T S   ( C F A S T )'
    '                                                          (REAL/IMAGINARY)'
    ' '
    ' ELEMENT-ID         DISP-X           DISP-Y           DISP-Z         ROTATION-X       ROTATION-Y       ROTATION-Z'
    '        100       2.638938E-04     1.068380E+00    -8.012853E-01    -1.814339E-07     9.713815E+01     1.295175E+02'
    '                 -1.073954E-05    -4.420905E-02     3.315679E-02     1.553854E-07    -4.019529E+00    -5.359372E+00'
    """
    out = []
    if is_stress:
        stress_strain = '                      C O M P L E X  S T R E S S E S  I N  F A S T E N E R  E L E M E N T S   ( C F A S T )\n'
        headers = ' ELEMENT-ID         FORCE-X          FORCE-Y          FORCE-Z         MOMENT-X         MOMENT-Y         MOMENT-Z\n'
    else:
        stress_strain = '              C O M P L E X  S T R A I N (D I S P)  I N  F A S T E N E R  E L E M E N T S   ( C F A S T )\n'
        headers = ' ELEMENT-ID         DISP-X           DISP-Y           DISP-Z         ROTATION-X       ROTATION-Y       ROTATION-Z\n'


    if is_mag_phase:
        raise NotImplementedError()
    else:
        mag_phase = '                                                          (REAL/IMAGINARY)\n'
    out = [stress_strain, mag_phase, ' \n', headers]
    return out
