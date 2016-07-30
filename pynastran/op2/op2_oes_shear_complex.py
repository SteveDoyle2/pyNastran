from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import get_key0
try:
    import pandas as pd
except ImportError:
    pass


class ComplexShearArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)   ## why???
        self.element_node = None
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        #self.itime = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_nnodes(self):
        return get_nnodes(self)

    def build(self):
        #print('data_code = %s' % self.data_code)
        if not hasattr(self, 'subtitle'):
            self.subtitle = self.data_code['subtitle']
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        self.ntotal = self.nelements * nnodes * 2
        #self.ntotal
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        #self.ntotal = self.nelements * nnodes

        self.element = zeros((self.nelements, 2), 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)

        # [max_shear, avg_shear]
        self.data = zeros((self.ntimes, self.ntotal, 2), 'complex64')

    def build_dataframe(self):
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = pd.Panel(self.data, items=column_values, major_axis=self.element, minor_axis=headers).to_frame()
        self.data_frame.columns.names = column_names
        self.data_frame.index.names = ['ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1():
                for itime in range(ntimes):
                    for ieid, eid in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (tx1, ty1, tz1, rx1, ry1, rz1) = t1
                        (tx2, ty2, tz2, rx2, ry2, rz2) = t2
                        d = t1 - t2
                        if not allclose([tx1.real, tx1.imag, ty1.real, ty1.imag],
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
                raise NotImplementedError(self.is_sort2())
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, max_shear, avg_shear):
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal] = [max_shear, avg_shear]
        self.element[self.itotal, :] = eid
        #self.ielement += 1
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.element_node.shape[0]
        #ntotal = self.ntotal
        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n' % (self.__class__.__name__, nelements, nnodes))
        msg.append('  data: [ntimes, nnodes, 2] where 2=[%s]\n' % str(', '.join(self._get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        raise NotImplementedError('ComplexShearArray')
        msg_temp, nnodes, is_bilinear = _get_plate_msg(self, is_mag_phase, is_sort1)

        ntimes = self.data.shape[0]
        eids = self.element
        if self.is_sort1():
            if is_sort1:
                for itime in range(ntimes):
                    dt = self._times[itime]

                    dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
                    header[1] = dt_line
                    msg = header + msg_temp
                    f.write('\n'.join(msg))

                    max_shear = self.data[itime, :, 0]
                    avg_shear = self.data[itime, :, 1]
                    for eid, max_sheari, avg_sheari in zip(eids, max_shear, avg_shear):
                        [rmax_shear, imax_shear, ravg_shear, iavg_shear
                         ,] = write_imag_floats_13e([max_sheari, avg_sheari], is_magnitude_phase)

                        f.write('   %6s   %-13s / %-13s     %-13s / %-13s\n' % (
                            eid, rmax_shear, imax_shear, ravg_shear, iavg_shear))
                    f.write(page_stamp % page_num)
                    page_num += 1
            else:
                times = self._times
                for ieid, eid in enumerate(eids):
                    max_shear = self.data[:, ieid, 0].ravel()
                    avg_shear = self.data[:, ieid, 1].ravel()
                    for itime, max_sheari, avg_sheari in zip(times, max_shear, avg_shear):
                        [rmax_shear, imax_shear, ravg_shear, iavg_shear
                         ] = write_imag_floats_13e([max_sheari, avg_sheari], is_magnitude_phase)

                        f.write('   %6s   %-13s / %-13s     %-13s / %-13s\n' % (
                            eid, rmax_shear, imax_shear, ravg_shear, iavg_shear))
                    f.write(page_stamp % page_num)
                    page_num += 1
        else:
            raise NotImplementedError('ComplexShearArray-sort2')
        return page_num - 1


class ComplexShearStressArray(ComplexShearArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def _get_headers(self):
        return ['max_shear', 'avg_shear']

class ComplexShearStrainArray(ComplexShearArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)
        assert self.is_strain(), self.stress_bits

    def _get_headers(self):
        return ['max_shear', 'avg_shear']
