from typing import List
import numpy as np
from numpy import zeros, float32, searchsorted, empty
from numpy import allclose, asarray, vstack

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.table_object import append_sort1_sort2
from pyNastran.op2.result_objects.op2_objects import BaseElement
from pyNastran.f06.f06_formatting import write_floats_13e, write_float_12e


class ElementTableArray(BaseElement):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.nonlinear_factor = np.nan
        self.table_name = None
        self.approach_code = None
        self.analysis_code = None
        BaseElement.__init__(self, data_code, isubcase, apply_data_code=True)  # no double inheritance
        self.is_sort1
        #self.dt = dt
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        self.ntotal = 0
        self.nelements = 0  # result specific

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        is_nan = (self.nonlinear_factor is not None and
                  np.isnan(self.nonlinear_factor) and
                  np.isnan(table.nonlinear_factor))
        if not is_nan:
            assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not is_nan:
            assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                self.element_name, self.element_type, self._times, table._times)
        if not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            for (eid, etype, eid2, etype2) in zip(self.element, self.element_data_type, table.element, table.element_data_type):
                msg += '(%s, %s)    (%s, %s)\n' % (eid, etype, eid2, etype2)
            print(msg)
            raise ValueError(msg)
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
                        (tx1, ty1, tz1, rx1, ry1, rz1) = t1
                        (tx2, ty2, tz2, rx2, ry2, rz2) = t2
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                tx1, ty1, tz1, rx1, ry1, rz1,
                                tx2, ty2, tz2, rx2, ry2, rz2)
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

    def combine(self, result, is_sort1=True):
        #print("combine; result=%s" % result)
        assert self.is_sort1 != result.is_sort1
        assert self.nonlinear_factor is not None
        assert result.nonlinear_factor is not None
        # self.ntimes += result.ntimes
        self.ntotal += result.data.shape[0]
        self.data = append_sort1_sort2(self.data, result.data)
        #print(self._times)
        #print(result._times)
        # self._times = hstack([self._times, result._times])
        self.node_gridtype = vstack([self.node_gridtype, result.node_gridtype])
        #print('%s' % ''.join(self.get_stats()))

    def _get_msgs(self, is_mag_phase):
        raise NotImplementedError()

    def data_type(self):
        raise NotImplementedError()

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]
        #ngrids = len(self.gridTypes)
        msg = []

        ntimesi, ntotal = self.data.shape[:2]
        ntimes = len(self._times)
        nelements = self.element.shape[0]

        nmajor = self.ntimes
        nminor = self.ntotal
        if self.is_sort1:
            assert nmajor == ntimes, 'ntimes=%s expected=%s' % (nmajor, ntimes)
            assert nminor == ntotal, 'ntotal=%s expected=%s' % (nminor, nelements)
        else:
            assert nmajor == nelements, 'nelements=%s expected=%s' % (nmajor, nelements)
            assert nminor == ntotal, 'ntotal=%s expected=%s' % (nminor, ntimes)

        msg.append('  isubcase = %s\n' % self.isubcase)
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n'
                       % (self.__class__.__name__, nelements))
        headers = ', '.join(self._get_headers())
        #msg.append('  data: [%s] shape=%s dtype=%s\n'
                   #% (headers, [int(i) for i in self.data.shape], self.data.dtype))
        msg.append('  data: [%s] shape=%s dtype=%s\n'
                   % (headers,
                      [int(i) for i in self.data.shape], self.data.dtype))
        msg += self.get_data_code()
        return msg

    @property
    def headers(self):
        return ['t1', 't2', 't3', 'r1', 'r2', 'r3']

    def _get_headers(self):
        return self.headers

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def build(self):
        """sizes the vectorized attributes of the ElementTableArray"""
        #print('_nelements=%s ntimes=%s sort1?=%s ntotal=%s -> _nelements=%s' % (
            #self._nelements, self.ntimes, self.is_sort1,
            #self.ntotal, self._nelements // self.ntimes))

        self.nelements //= self.ntimes
        self.itime = 0
        self.itotal = 0
        self.is_built = True

        if self.is_sort1:
            ntimes = self.ntimes
            nelements = self.ntotal
            nx = ntimes
            ny = self.ntotal
            #print("ntimes=%s nelements=%s" % (ntimes, nelements))
        if self.is_sort2:
            unused_ntotal = self.ntotal
            nelements = self.ntimes
            ntimes = self.ntotal
            nx = nelements
            ny = ntimes
            #print("ntotal=%s nelements=%s ntimes=%s" % (ntotal, nelements, ntimes))

        self._times = zeros(ntimes, dtype=self._times_dtype)
        #self.types = array(self.nelements, dtype='|S1')

        self.element = zeros(nelements, dtype='int32')
        self.element_data_type = empty(nelements, dtype='|U8')

        #[t1, t2, t3, r1, r2, r3]
        self.data = zeros((nx, ny, 6), self.data_type())

    def add_sort1(self, dt, eid, etype, v1, v2, v3, v4, v5, v6):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        # itotal - the node number
        # itime - the time/frequency step

        # the times/freqs
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        self.element_data_type[self.itotal] = etype
        self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]
        self.itotal += 1

    def add_sort2(self, dt, eid, etype, v1, v2, v3, v4, v5, v6):
        msg = "dt=%s eid=%s etype=%s v1=%s v2=%s v3=%s\n" % (dt, eid, etype, v1, v2, v3)
        msg += "                     v4=%s v5=%s v6=%s" % (v4, v5, v6)
        self._times[self.itotal] = dt

        if 1:  # this is needed for SORT1 tables
            inode = self.itime
            self.element[self.itime] = eid
            self.element_data_type[self.itime, :] = etype
            self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]
            # itotal - the node number
            # itime - the time/frequency step
        else:
            self.element[self.itime] = eid
            self.element_data_type[self.itime] = etype
            self.data[self.itotal, self.itime, :] = [v1, v2, v3, v4, v5, v6]
            # itotal - the time/frequency step
            # itime - the node number

        self.itotal += 1
        #self.itime += 1


class RealElementTableArray(ElementTableArray):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ElementTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    @property
    def nnodes_per_element(self) -> bool:
        return 1

    def data_type(self) -> str:
        return 'float32'

    #def spike():
        #import xlwings as xw
        #wb = xw.Workbook()  # Creates a connection with a new workbook
        #xw.Range('A1').value = 'Foo 1'
        #xw.Range('A1').value
        #'Foo 1'
        #xw.Range('A1').value = [['Foo 1', 'Foo 2', 'Foo 3'], [10.0, 20.0, 30.0]]
        #xw.Range('A1').table.value  # or: Range('A1:C2').value
        #[['Foo 1', 'Foo 2', 'Foo 3'], [10.0, 20.0, 30.0]]
        #xw.Sheet(1).name
        #'Sheet1'
        #chart = xw.Chart.add(source_data=xw.Range('A1').table)

    def _write_f06_block(self, words, header, page_stamp, page_num, f,
                         is_mag_phase=False, is_sort1=True):
        #words += self.getTableMarker()
        f.write(''.join(header + words))

        element = self.element
        element_type = self.element_data_type
        t1 = self.data[0, :, 0]
        t2 = self.data[0, :, 1]
        t3 = self.data[0, :, 2]
        r1 = self.data[0, :, 3]
        r2 = self.data[0, :, 4]
        r3 = self.data[0, :, 5]
        for element_id, etypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(element, element_type, t1, t2, t3, r1, r2, r3):
            vals = [t1i, t2i, t3i, r1i, r2i, r3i]
            (dx, dy, dz, rx, ry, rz) = write_floats_13e(vals)
            f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                element_id, etypei, dx, dy, dz, rx, ry, rz))
        f.write(page_stamp % page_num)
        return page_num

    def _write_sort1_as_sort2(self, f06_file, page_num, page_stamp, header, words):
        element = self.element
        element_type = self.element_data_type
        times = self._times

        node_id = 0  ## TODO: fix the node id
        for inode, (eid, etypei) in enumerate(zip(element, element_type)):
            t1 = self.data[:, inode, 0].ravel()
            t2 = self.data[:, inode, 1].ravel()
            t3 = self.data[:, inode, 2].ravel()
            r1 = self.data[:, inode, 3].ravel()
            r2 = self.data[:, inode, 4].ravel()
            r3 = self.data[:, inode, 5].ravel()

            header[1] = ' POINT-ID = %10i\n' % node_id
            f06_file.write(''.join(header + words))
            for dt, t1i, t2i, t3i, r1i, r2i, r3i in zip(times, t1, t2, t3, r1, r2, r3):
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                vals2 = write_floats_13e(vals)
                (dx, dy, dz, rx, ry, rz) = vals2
                f06_file.write('%14s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                    write_float_12e(dt), etypei, dx, dy, dz, rx, ry, rz))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort1_as_sort1(self, f06_file, page_num, page_stamp, header, words):
        element = self.element
        element_type = self.element_data_type
        times = self._times

        for itime in range(self.ntimes):
            dt = self._times[itime]
            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            if isinstance(dt, (float, float32)):
                header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            else:
                header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            f06_file.write(''.join(header + words))
            for element_id, etypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(element, element_type, t1, t2, t3, r1, r2, r3):
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                vals2 = write_floats_13e(vals)
                (dx, dy, dz, rx, ry, rz) = vals2
                f06_file.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                    element_id, etypei, dx, dy, dz, rx, ry, rz))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort2_as_sort2(self, f06_file, page_num, page_stamp, header, words):
        element = self.element
        element_type = self.element_data_type
        times = self._times
        for ieid, (element_id, etypei) in enumerate(zip(element, element_type)):
            t1 = self.data[ieid, :, 0]
            t2 = self.data[ieid, :, 1]
            t3 = self.data[ieid, :, 2]
            r1 = self.data[ieid, :, 3]
            r2 = self.data[ieid, :, 4]
            r3 = self.data[ieid, :, 5]

            header[1] = ' ELEMENT-ID = %10i\n' % element_id
            f06_file.write(''.join(header + words))
            for dt, t1i, t2i, t3i, r1i, r2i, r3i in zip(times, t1, t2, t3, r1, r2, r3):
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                vals2 = write_floats_13e(vals)
                (dx, dy, dz, rx, ry, rz) = vals2
                f06_file.write('%14s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                    write_float_12e(dt), etypei, dx, dy, dz, rx, ry, rz))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_f06_transient_block(self, words, header, page_stamp, page_num, f06_file,
                                   is_mag_phase=False, is_sort1=True):
        #words += self.getTableMarker()

        if not len(header) >= 3:
            header.append('')

        is_sort2 = not is_sort1
        if self.is_sort1 or self.nonlinear_factor in (None, np.nan):
            if is_sort2 and self.nonlinear_factor is not None:
                page_num = self._write_sort1_as_sort2(f06_file, page_num, page_stamp, header, words)
            else:
                page_num = self._write_sort1_as_sort1(f06_file, page_num, page_stamp, header, words)
        else:
            page_num = self._write_sort2_as_sort2(f06_file, page_num, page_stamp, header, words)
        return page_num - 1

    def extract_xyplot(self, element_ids, index):
        element_ids = asarray(element_ids, dtype='int32')
        i = index - 1
        assert index in [1, 2, 3, 4, 5, 6], index
        eids = self.element
        ieids = searchsorted(eids, element_ids)
        assert all(eids[ieids] == element_ids), 'nids=%s expected=%s; all=%s'  % (eids[ieids], element_ids, eids)
        return self.data[:, ieids, i]


#class ComplexElementTableArray(ElementTableArray):  # displacement style table
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #TableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    #def extract_xyplot(self, node_ids, index, index_str):
        #index_str = index_str.lower().strip()
        #if index_str in ['real', 'r']:
            #j = 1
        #elif index_str in ['imag', 'i']:
            #j = 2
        #elif index_str in ['mag', 'magnitude', 'm']:
            #j = 3
        #elif index_str in ['phase', 'p']:
            #j = 4
        #else:
            #raise ValueError('index_str=%r' % index_str)

        #node_ids = asarray(node_ids, dtype='int32')
        #i = index - 1
        #assert index in [1, 2, 3, 4, 5, 6,
                         #7, 8, 9, 10, 11, 12], index
        #nids = self.node_gridtype[:, 0]
        #inids = searchsorted(nids, node_ids)
        #assert all(nids[inids] == node_ids), 'nids=%s expected=%s; all=%s'  % (nids[inids], node_ids, nids)
        #if j == 1:
            ## real
            #return self.data[:, inids, i].real
        #elif j == 2:
            ## imag
            #return self.data[:, inids, i].imag
        #elif j == 3:
            ## mag
            #return np.abs(self.data[:, inids, i])
        #elif j == 4:
            ## phase
            #return angle(self.data[:, inids, i])
        #else:
            #raise RuntimeError()

    #@property
    #def is_real(self):
        #return False

    #@property
    #def is_complex(self):
        #return True

    #def data_type(self):
        #return 'complex64'

    ##def _write_f06_block(self, words, header, page_stamp, page_num, f06_file, is_mag_phase):
        ##self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, is_mag_phase, is_sort1)

    #def _write_f06_transient_block(self, words, header, page_stamp, page_num, f06_file,
                                   #is_mag_phase, is_sort1):
        #if is_mag_phase:
            #words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        #else:
            #words += ['                                                          (REAL/IMAGINARY)\n', ]

        #if not len(header) >= 3:
            #header.append('')

        #if self.is_sort1:
            #if is_sort1:
                #words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
                #page_num = self.write_sort1_as_sort1(f06_file, page_num, page_stamp, header, words, is_mag_phase)
            #else:
                #words += [' \n', '      FREQUENCY   TYPE          T1             T2             T3             R1             R2             R3\n']
                #page_num = self.write_sort1_as_sort2(f06_file, page_num, page_stamp, header, words, is_mag_phase)
        #else:
            #words += [' \n', '      FREQUENCY   TYPE          T1             T2             T3             R1             R2             R3\n']
            #page_num = self.write_sort2_as_sort2(f06_file, page_num, page_stamp, header, words, is_mag_phase)
        #return page_num - 1

    #def write_sort1_as_sort1(self, f06_file, page_num, page_stamp, header, words, is_mag_phase):
        #assert self.ntimes == len(self._times), 'ntimes=%s len(self._times)=%s' % (self.ntimes, self._times)
        #for itime, dt in enumerate(self._times):
            #node = self.node_gridtype[:, 0]
            #gridtype = self.node_gridtype[:, 1]
            #t1 = self.data[itime, :, 0]
            #t2 = self.data[itime, :, 1]
            #t3 = self.data[itime, :, 2]
            #r1 = self.data[itime, :, 3]
            #r2 = self.data[itime, :, 4]
            #r3 = self.data[itime, :, 5]

            #header[2] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            #f.write(''.join(header + words))
            #for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
                #sgridtype = self.recast_gridtype_as_string(gridtypei)
                #vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                #vals2 = write_imag_floats_13e(vals, is_mag_phase)
                #[dxr, dyr, dzr, rxr, ryr, rzr,
                 #dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                #if sgridtype == 'G':
                    #f06_file.write('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                              #'  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                  #node_id, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                  #'', '', dxi, dyi, dzi, rxi, ryi, rzi))
                #elif sgridtype == 'S':
                    #f06_file.write('0 %12i %6s     %-13s\n'
                              #'  %12s %6s     %-13s\n' % (node_id, sgridtype, dxr, '', '', dxi))
                #else:
                    #raise NotImplementedError(sgridtype)
            #f.write(page_stamp % page_num)
            #page_num += 1
        #return page_num

    #def write_sort1_as_sort2(self, f06_file, page_num, page_stamp, header, words, is_mag_phase):
        #node = self.node_gridtype[:, 0]
        #gridtype = self.node_gridtype[:, 1]

        #times = self._times
        #print(self.data.shape)
        #for inode, (node_id, gridtypei) in enumerate(zip(node, gridtype)):
            ## SORT1 is pretending to be SORT2
            #t1 = self.data[:, inode, 0].ravel()
            #t2 = self.data[:, inode, 1].ravel()
            #t3 = self.data[:, inode, 2].ravel()
            #r1 = self.data[:, inode, 3].ravel()
            #r2 = self.data[:, inode, 4].ravel()
            #r3 = self.data[:, inode, 5].ravel()
            #if len(r3) != len(times):
                #raise RuntimeError('len(d)=%s len(times)=%s' % (len(r3), len(times)))

            #header[2] = ' POINT-ID = %10i\n' % node_id
            #f.write(''.join(header + words))
            #for dt, t1i, t2i, t3i, r1i, r2i, r3i in zip(times, t1, t2, t3, r1, r2, r3):
                #sgridtype = self.recast_gridtype_as_string(gridtypei)
                #vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                #vals2 = write_imag_floats_13e(vals, is_mag_phase)
                #[dxr, dyr, dzr, rxr, ryr, rzr,
                 #dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                #sdt = write_float_12e(dt)
                ##if not is_all_zeros:
                #if sgridtype == 'G':
                    #f06_file.write('0 %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                              #'  %13s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                  #sdt, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                  #'', '', dxi, dyi, dzi, rxi, ryi, rzi))
                #elif sgridtype == 'S':
                    #f06_file.write('0 %12s %6s     %-13s\n'
                              #'  %12s %6s     %-13s\n' % (sdt, sgridtype, dxr, '', '', dxi))
                #else:
                    #msg = 'nid=%s dt=%s type=%s dx=%s dy=%s dz=%s rx=%s ry=%s rz=%s' % (
                        #node_id, dt, sgridtype, t1i, t2i, t3i, r1i, r2i, r3i)
                    #raise NotImplementedError(msg)
            #f06_file.write(page_stamp % page_num)
            #page_num += 1
        #return page_num

    #def write_sort2_as_sort2(self, f06_file, page_num, page_stamp, header, words, is_mag_phase):
        #node = self.node_gridtype[:, 0]
        #gridtype = self.node_gridtype[:, 1]

        #times = self._times
        #for inode, (node_id, gridtypei) in enumerate(zip(node, gridtype)):
            ## TODO: for SORT1 pretending to be SORT2
            ##t1 = self.data[:, inode, 0].ravel()
            #t1 = self.data[inode, :, 0]
            #t2 = self.data[inode, :, 1]
            #t3 = self.data[inode, :, 2]
            #r1 = self.data[inode, :, 3]
            #r2 = self.data[inode, :, 4]
            #r3 = self.data[inode, :, 5]
            #if len(r3) != len(times):
                #raise RuntimeError('len(d)=%s len(times)=%s' % (len(r3), len(times)))

            #header[2] = ' POINT-ID = %10i\n' % node_id
            #f.write(''.join(header + words))
            #for dt, t1i, t2i, t3i, r1i, r2i, r3i in zip(times, t1, t2, t3, r1, r2, r3):
                #sgridtype = self.recast_gridtype_as_string(gridtypei)
                #vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                #vals2 = write_imag_floats_13e(vals, is_mag_phase)
                #[dxr, dyr, dzr, rxr, ryr, rzr,
                 #dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                #sdt = write_float_12e(dt)
                ##if not is_all_zeros:
                #if sgridtype == 'G':
                    #f06_file.write('0 %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                              #'  %13s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                  #sdt, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                  #'', '', dxi, dyi, dzi, rxi, ryi, rzi))
                #elif sgridtype == 'S':
                    #f06_file.write('0 %12s %6s     %-13s\n'
                              #'  %12s %6s     %-13s\n' % (sdt, sgridtype, dxr, '', '', dxi))
                #else:
                    #msg = 'nid=%s dt=%s type=%s dx=%s dy=%s dz=%s rx=%s ry=%s rz=%s' % (
                        #node_id, dt, sgridtype, t1i, t2i, t3i, r1i, r2i, r3i)
                    #raise NotImplementedError(msg)
            #f06_file.write(page_stamp % page_num)
            #page_num += 1
        #return page_num

