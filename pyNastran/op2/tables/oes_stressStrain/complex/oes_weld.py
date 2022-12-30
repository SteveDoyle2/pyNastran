import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import get_complex_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_imag_floats_13e


class ComplexWeldArray(OES_Object):
    """
    Common class for:
     - ComplexWeldStressArray
     - ComplexWeldStrainArray
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
        self._times = np.zeros(self.ntimes, dtype=dtype)
        #self.ntotal = self.nelements * nnodes

        self.element = np.zeros(self.nelements, dtype=idtype)

        # the number is messed up because of the offset for the element's properties
        if self.nelements != self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
                self.ntimes, self.nelements, nnodes, self.nelements * nnodes,
                self.ntotal)
            raise RuntimeError(msg)

        # stress: [axial, max_a, min_a, max_b, min_b, max_shear, bearing]
        # strain: [axial, max_a, min_a, max_b, min_b, max_shear, bearing]
        self.data = np.zeros((self.ntimes, self.ntotal, 7), dtype=cfdtype)

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
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = self._build_pandas_transient_elements(
            column_values, column_names,
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

    def add_sort1(self, dt, eid, axial, max_a, min_a, max_b, min_b, max_shear, bearing):
        """unvectorized method for adding SORT1 transient data"""
        assert self.sort_method == 1, self
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal] = [axial, max_a, min_a, max_b, min_b, max_shear, bearing]
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
        msg_temp = _get_weld_msg(self.is_stress, is_mag_phase, is_sort1)
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

                    axial = self.data[itime, :, 0]
                    max_a = self.data[itime, :, 1]
                    min_a = self.data[itime, :, 2]
                    max_b = self.data[itime, :, 3]
                    min_b = self.data[itime, :, 4]
                    max_shear = self.data[itime, :, 5]
                    bearing = self.data[itime, :, 6]
                    assert len(eids) == len(force_x)
                    assert len(axial) > 0, axial
                    for eid, axiali, max_ai, min_ai, max_bi, min_bi, max_sheari, bearingi in zip(eids, axial, max_a, min_a, max_b, min_b, max_shear, bearing):
                        assert isinstance(eid, integer_types), 'eid=%s type=%s' % (eid, type(eid))
                        [raxial, iaxial,
                         rmax_a, imax_a, rmin_a, imin_a,
                         rmax_b, imax_b, rmin_b, imin_b,
                         rmax_shear, imax_shear, rbearing, ibearing, ] = write_imag_floats_13e([axiali, max_ai, min_ai, max_bi, min_bi, max_sheari, bearingi], is_mag_phase)

                        f06_file.write(
                            '%24s  %13s  %13s  %13s  %13s  %13s  %13s  %s\n'
                            '%24s  %13s  %13s  %13s  %13s  %13s  %13s  %s\n' % (
                                eid, raxial, rmax_a, rmin_a, rmax_b, rmin_b, rmax_shear, rbearing,
                                '',  iaxial, imax_a, imin_a, imax_b, imin_b, imax_shear, ibearing))
                    f06_file.write(page_stamp % page_num)
                    page_num += 1
            else:
                # TODO: write in SORT2
                times = self._times
                raise NotImplementedError('WELD SORT2')
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
            raise NotImplementedError('ComplexWeldArray-sort2')
        return page_num - 1

    @property
    def headers(self) -> list[str]:
        return self._get_headers()

    def get_headers(self) -> list[str]:
        return self.headers

class ComplexWeldStressArray(ComplexWeldArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexWeldArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self) -> list[str]:
        """
        '                        C O M P L E X   S T R E S S E S   I N   W E L D   E L E M E N T S   ( C W E L D P ) '
        '                                                          (REAL/IMAGINARY)'
        ' '
        '    ELEMENT          AXIAL         MAX  STRESS      MIN  STRESS      MAX  STRESS      MIN  STRESS        MAXIMUM          BEARING '
        '      ID             STRESS           END-A            END-A            END-B            END-B        SHEAR  STRESS       STRESS'
        '         30       1.146487E+00     1.278808E+03    -1.276515E+03     6.757888E+01    -6.528591E+01     1.298489E+02     3.134326E+01'
        '                 -1.373695E-01     1.532364E+02    -1.535111E+02     7.883755E+00    -8.158494E+00     1.559507E+01     3.761247E+00'
        """
        return ['axial', 'max_a', 'min_a', 'max_b', 'min_b', 'max_shear', 'bearing']

class ComplexWeldStrainArray(ComplexWeldArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexWeldArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain, self.stress_bits

    def _get_headers(self):
        return ['disp_x', 'disp_y', 'disp_z', 'rotation_x', 'rotation_y', 'rotation_z']

def _get_weld_msg(is_stress: bool, is_mag_phase: bool, is_sort1: bool) -> list[str]:
    """
    '                        C O M P L E X   S T R E S S E S   I N   W E L D   E L E M E N T S   ( C W E L D P ) '
    '                                                          (REAL/IMAGINARY)'
    ' '
    '    ELEMENT          AXIAL         MAX  STRESS      MIN  STRESS      MAX  STRESS      MIN  STRESS        MAXIMUM          BEARING '
    '      ID             STRESS           END-A            END-A            END-B            END-B        SHEAR  STRESS       STRESS'
    '         30       1.146487E+00     1.278808E+03    -1.276515E+03     6.757888E+01    -6.528591E+01     1.298489E+02     3.134326E+01'
    '                 -1.373695E-01     1.532364E+02    -1.535111E+02     7.883755E+00    -8.158494E+00     1.559507E+01     3.761247E+00'
    """
    out = []
    if is_stress:
        stress_strain = '                        C O M P L E X   S T R E S S E S   I N   W E L D   E L E M E N T S   ( C W E L D P ) \n'
        headers = (
            '    ELEMENT          AXIAL         MAX  STRESS      MIN  STRESS      MAX  STRESS      MIN  STRESS        MAXIMUM          BEARING \n'
            '      ID             STRESS           END-A            END-A            END-B            END-B        SHEAR  STRESS       STRESS\n')
    else:
        stress_strain = '              C O M P L E X  S T R A I N (D I S P)  I N  F A S T E N E R  E L E M E N T S   ( C F A S T )\n'
        #headers = ' ELEMENT-ID         DISP-X           DISP-Y           DISP-Z         ROTATION-X       ROTATION-Y       ROTATION-Z\n'
        raise RuntimeError('weld-strain')


    if is_mag_phase:
        raise NotImplementedError()
    else:
        mag_phase = '                                                          (REAL/IMAGINARY)\n'
    out = [stress_strain, mag_phase, ' \n', headers]
    return out
