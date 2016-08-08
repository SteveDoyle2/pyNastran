from __future__ import print_function, unicode_literals
from struct import Struct, pack
from six.moves import zip, range

import numpy as np
from numpy import zeros, float32, searchsorted, unique, where
from numpy import allclose, asarray, vstack, array_equal

from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.op2.result_objects.table_object import append_sort1_sort2
from pyNastran.f06.f06_formatting import write_floats_13e, write_float_12E
try:
    import pandas as pd
except ImportError:
    pass


class ScalarTableArray(ScalarObject):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.nonlinear_factor = None
        self.table_name = None
        self.approach_code = None
        self.analysis_code = None
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)  # no double inheritance
        self.is_sort1()
        self.is_sort2()
        #self.dt = dt

        #self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        self.ntotal = 0
        self._nnodes = 0  # result specific

    def __eq__(self, table):
        self._eq_header(table)
        assert self.is_sort1() == table.is_sort1()
        if not np.array_equal(self.node_gridtype, table.node_gridtype):
            assert self.node_gridtype.shape == table.node_gridtype.shape, 'shape=%s table.shape=%s' % (self.node_gridtype.shape, table.node_gridtype.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            for (nid, grid_type), (nid2, grid_type2) in zip(self.node_gridtype, table.node_gridtype):
                msg += '(%s, %s)    (%s, %s)\n' % (nid, grid_type, nid2, grid_type2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1():
                for itime in range(ntimes):
                    for inid, nid_gridtype, in enumerate(self.node_gridtype):
                        (nid, grid_type) = nid_gridtype
                        t1 = self.data[itime, inid, 0]
                        t2 = table.data[itime, inid, 0]
                        tx1 = t1[0]
                        tx2 = t2[0]
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '(%s, %s)\n  (%s, %s)\n' % (
                                nid, grid_type, tx1, tx2)
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

    def combine(self, result, is_sort1=True):
        #print("combine; result=%s" % result)
        assert self.is_sort1() != result.is_sort1()
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

    def get_stats(self, short=False):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]
        #ngrids = len(self.gridTypes)
        if short:
            return self._get_stats_short()
        msg = []

        ntimesi, ntotal = self.data.shape[:2]
        ntimes = len(self._times)
        nnodes = self.node_gridtype.shape[0]

        nmajor = self.ntimes
        nminor = self.ntotal
        if self.is_sort1():
            assert nmajor == ntimes, 'ntimes=%s expected=%s' % (nmajor, ntimes)
            assert nminor == ntotal, 'ntotal=%s expected=%s' % (nminor, nnodes)
        else:
            assert nmajor == nnodes, 'nnodes=%s expected=%s' % (nmajor, nnodes)
            assert nminor == ntotal, 'ntotal=%s expected=%s' % (nminor, ntimes)

        msg.append('  isubcase = %s\n' % self.isubcase)
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%s nnodes=%s\n'
                       % (self.__class__.__name__, ntimes, nnodes))
        else:
            msg.append('  type=%s nnodes=%s\n'
                       % (self.__class__.__name__, nnodes))
        headers = ', '.join(self._get_headers())
        #msg.append('  data: [%s] shape=%s dtype=%s\n'
                   #% (headers, [int(i) for i in self.data.shape], self.data.dtype))
        msg.append('  data: [%s] shape=%s dtype=%s\n'
                   % (headers,
                      [int(i) for i in self.data.shape], self.data.dtype))
        msg.append('  gridTypes\n  ')
        msg += self.get_data_code()
        return msg

    @property
    def headers(self):
        return ['t1']

    def _get_headers(self):
        return self.headers

    def get_headers(self):
        return self._get_headers()

    def _reset_indices(self):
        self.itotal = 0

    def build(self):
        #print('_nnodes=%s ntimes=%s sort1?=%s ntotal=%s -> _nnodes=%s' % (self._nnodes, self.ntimes, self.is_sort1(),
                                                                          #self.ntotal, self._nnodes // self.ntimes))
        if self.is_built:
            #print("resetting...")
            #self.itotal = 0
            return

        self._nnodes //= self.ntimes
        self.itime = 0
        self.itotal = 0
        self.is_built = True

        if self.is_sort1():
            ntimes = self.ntimes
            nnodes = self.ntotal
            nx = ntimes
            ny = self.ntotal
            #print("ntimes=%s nnodes=%s" % (ntimes, nnodes))
        if self.is_sort2():
            ntotal = self.ntotal
            nnodes = self.ntimes
            ntimes = self.ntotal
            nx = nnodes
            ny = ntimes
            #print("ntotal=%s nnodes=%s ntimes=%s" % (ntotal, nnodes, ntimes))

        self._times = zeros(ntimes, dtype=self._times_dtype)
        #self.types = array(self.nelements, dtype='|S1')

        self.node_gridtype = zeros((nnodes, 2), dtype='int32')

        #[t1]
        self.data = zeros((nx, ny, 1), self.data_type())

    def build_dataframe(self):
        headers = self.get_headers()
        node_gridtype = [self.node_gridtype[:, 0], self.gridtype_str]
        ugridtype_str = unique(self.gridtype_str)

        if self.nonlinear_factor is not None:
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=node_gridtype, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['NodeID', 'Type', 'Item']

            letter_dims = [
                ('G', 1),
                ('E', 1),
                ('S', 1),
                ('H', 1),
                ('L', 1),
            ]
            cat_keys = []
            for (letter, dim) in letter_dims:
                if letter not in ugridtype_str:
                    continue
                eig = self.data_frame.xs(letter, level=1)
                cat_keys.append(eig)
            self.data_frame = pd.concat(cat_keys)
        else:
            #self.data_frame = pd.Panel(self.data[0, :, :], major_axis=node_gridtype, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = ['Static']
            #self.data_frame.index.names = ['NodeID', 'Type', 'Item']

            df1 = pd.DataFrame(self.node_gridtype[:, 0])
            df1.columns = ['NodeID']
            df2 = pd.DataFrame(self.gridtype_str)
            df2.columns = ['Type']
            df3 = pd.DataFrame(self.data[0])
            df3.columns = headers
            self.data_frame = df1.join([df2, df3])
        #print(self.data_frame)

    def finalize(self):
        gridtypes = self.node_gridtype[:, 1]
        nnodes = len(gridtypes)
        self.gridtype_str = np.chararray((nnodes), unicode=True)
        ugridtypes = unique(gridtypes)
        for ugridtype in ugridtypes:
            i = where(gridtypes == ugridtype)
            self.gridtype_str[i] = self.recast_gridtype_as_string(ugridtype)

    def _write_xlsx(self, sheet, is_mag_phase=False):
        from xlwings import Range, Chart
        # 0.3.5 doesn't work, 0.5 does
        #from numpy import astype
        # print('xlsx_filename = %r' % xlsx_filename)
        #f = None
        #wb = Workbook()  # Creates a connection with a new workbook
        #wb.save(xlsx_filename)
        #Range('A1').value = 'Foo 1'
        #print(Range('A1').value)
        #'Foo 1'
        # Range('A1').value = xlsx_filename
        name = str(self.__class__.__name__)
        Range(sheet, 'A1').value = [name]
        Range(sheet, 'A2').value = ['Node', 'GridType'] + self.headers
        Range(sheet, 'A3').value = self.node_gridtype

        if self.is_real():
            Range(sheet, 'C3').value = self.data[0, :, :]
        else:
            pass
            #from numpy.core.defchararray import add as sadd
            #n, m = self.data[0, :, :].shape
            #nm = n * m
            #scomplex = array(['=complex('] * nm, dtype='|S10').reshape(n, m)
            #scomma = array([','] * nm, dtype='|S40').reshape(n, m)
            #sparen = array([')'] * nm, dtype='|S40').reshape(n, m)
            #data = sadd(
                #sadd(scomplex, self.data.real.astype('|S10')), # complex(5.
                #sadd(
                    #scomma, # ,
                    #sadd(self.data.imag.astype('|U10'), sparen), # 3j)
                #)
            #)
            #data = sadd(
                #scomplex,
                #self.data.real.astype('|S10'),
                #scomma,
                #self.data.imag.astype('|S10'),
                #sparen)
            #print(self.data.real)
            #Range(sheet, 'C3', atleast_2d=True).table.value = self.data.real
            #Range(sheet, 'C3').value = self.data.real
        #Range('C4').value = self.data[0, :, 0]
        #Range('D4').value = self.data[0, :, 1:]
        #print(Range('A1').table.value)  # or: Range('A1:C2').value
        #[['Foo 1', 'Foo 2', 'Foo 3'], [10.0, 20.0, 30.0]]
        #print(Sheet(1).name)
        #Sheet(isheet).name = 'displacements'
        #'Sheet1'
        #nrows = self.data.shape[1]
        #end_row = '%s' % (4 + nrows)
        #t1 = self.data[0, :, 0]
        #chart = Chart.add(source_data=Range('C4').value)
        #wb.save()
        # wb.save()


    def add(self, node_id, grid_type, v1):
        self.add_sort1(None, node_id, grid_type, v1)

    def add_sort1(self, dt, node_id, grid_type, v1):
        # itotal - the node number
        # itime - the time/frequency step

        # the times/freqs
        self._times[self.itime] = dt
        self.node_gridtype[self.itotal, :] = [node_id, grid_type]
        self.data[self.itime, self.itotal, 0] = v1
        self.itotal += 1

    def add_sort2(self, dt, node_id, grid_type, v1):
        self._times[self.itotal] = dt

        if 1:  # this is needed for SORT1 tables
            inode = self.itime
            self.node_gridtype[self.itime, :] = [node_id, grid_type]
            self.data[self.itime, self.itotal, 0] = v1
            # itotal - the node number
            # itime - the time/frequency step
        else:
            self.node_gridtype[self.itime, :] = [node_id, grid_type]
            self.data[self.itotal, self.itime, 0] = v1
            # itotal - the time/frequency step
            # itime - the node number

        self.itotal += 1
        #self.itime += 1

#def two_dee_string_add(string_lists):
    #string0 = string_lists[0]
    #n, m = string0.shape

    #s = []
    #for string_list in string_lists:
        #for string in string_list:
            #pass
    #return sumned

class RealScalarTableArray(ScalarTableArray):  # temperature style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def data_type(self):
        return 'float32'

    def _write_table_3(self, f, fascii, itable=-3, itime=0):
        import inspect
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        fascii.write('%s.write_table_3: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        f.write(pack('12i', *[4, itable, 4,
                              4, 1, 4,
                              4, 0, 4,
                              4, 146, 4,
                              ]))
        approach_code = self.approach_code
        table_code = self.table_code
        isubcase = self.isubcase
        random_code = self.random_code
        format_code = 1
        num_wide = self.num_wide
        acoustic_flag = 0
        thermal = 0
        title = b'%-128s' % bytes(self.title)
        subtitle = b'%-128s' % bytes(self.subtitle)
        label = b'%-128s' % bytes(self.label)
        ftable3 = b'50i 128s 128s 128s'
        oCode = 0
        if self.analysis_code == 1:
            lsdvmn = self.lsdvmn
        else:
            raise NotImplementedError(self.analysis_code)

        table3 = [
            approach_code, table_code, 0, isubcase, lsdvmn,
            0, 0, random_code, format_code, num_wide,
            oCode, acoustic_flag, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, thermal, thermal, 0,
            title, subtitle, label,
        ]

        n = 0
        for v in table3:
            if isinstance(v, (int, float)):
                n += 4
            else:
                n += len(v)
        assert n == 584, n
        data = [584] + table3 + [584]
        fmt = 'i' + ftable3 + 'i'
        #print(fmt)
        #f.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        fascii.write('%s header 3c = %s\n' % (self.table_name, data))
        f.write(pack(fmt, *data))

    def write_op2(self, f, fascii, itable, date, is_mag_phase=False, endian='>'):
        import inspect
        assert self.table_name in ['OUGV1', 'OQMG1', 'OQG1'], self.table_name

        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        fascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        #print('data_code =', self.data_code)
        if itable == -1:
            self._write_table_header(f, fascii, date)
            itable = -3

        if isinstance(self.nonlinear_factor, float):
            op2_format = endian + b'%sif' % (7 * self.ntimes)
            raise NotImplementedError()
        else:
            op2_format = endian + b'2i6f' * self.ntimes
        s = Struct(op2_format)

        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]
        #format_table4_1 = Struct(b(self._endian + '15i'))
        #format_table4_2 = Struct(b(self._endian + '3i'))

        # table 4 info
        #ntimes = self.data.shape[0]
        nnodes = self.data.shape[1]
        nnodes_device = self.node_gridtype[:, 0] * 10 + self.device_code

        #(2+6) => (node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i)
        ntotal = self.ntimes * nnodes * (2 + 6)

        #print('shape = %s' % str(self.data.shape))
        assert nnodes > 1, nnodes
        assert ntotal > 1, ntotal

        device_code = self.device_code
        fascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        for itime in range(self.ntimes):
            self._write_table_3(f, fascii, itable, itime)

            # record 4
            header = [4, -4, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4*ntotal]
            f.write(pack(b'%ii' % len(header), *header))
            fascii.write('r4 [4, 0, 4]\n')
            fascii.write('r4 [4, %s, 4]\n' % (itable-1))
            fascii.write('r4 [4, %i, 4]\n' % (4*ntotal))

            t1 = self.data[itime, :, 0]
            for node_id, gridtypei, t1i in zip(nnodes_device, gridtype, t1):
                data = [node_id, gridtypei, t1i, 0., 0., 0., 0., 0.]
                fascii.write('  nid, grid_type, dx, dy, dz, rx, ry, rz = %s\n' % data)
                f.write(s.pack(*data))

            itable -= 2
            header = [4 * ntotal,]
            f.write(pack(b'i', *header))
            fascii.write('footer = %s' % header)
        header = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4,
        ]
        f.write(pack(b'%ii' % len(header), *header))
        return itable

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

    def _write_f06_block(self, words, header, page_stamp, page_num, f, write_words,
                         is_mag_phase=False, is_sort1=True):
        if write_words:
            words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()
        f.write(''.join(header + words))

        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]
        t1 = self.data[0, :, 0]
        for node_id, gridtypei, t1i in zip(node, gridtype, t1):
            sgridtype = self.recast_gridtype_as_string(gridtypei)
            vals = [t1i]
            vals2 = write_floats_13e(vals)
            dx = vals2[0]
            f.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
        f.write(page_stamp % page_num)
        return page_num

    def _write_sort1_as_sort2(self, f, page_num, page_stamp, header, words):
        nodes = self.node_gridtype[:, 0]
        gridtypes = self.node_gridtype[:, 1]
        times = self._times

        for inode, (node_id, gridtypei) in enumerate(zip(nodes, gridtypes)):
            t1 = self.data[:, inode, 0].ravel()

            header[1] = ' POINT-ID = %10i\n' % node_id
            f.write(''.join(header + words))
            for dt, t1i in zip(times, t1):
                sgridtype = self.recast_gridtype_as_string(gridtypei)
                vals = [t1i]
                vals2 = write_floats_13e(vals)
                dx = vals2[0]
                if sgridtype == 'G':
                    f.write('%14s %6s     %s\n' % (write_float_12E(dt), sgridtype, dx))
                elif sgridtype == 'S':
                    f.write('%14s %6s     %s\n' % (node_id, sgridtype, dx))
                elif sgridtype == 'H':
                    f.write('%14s %6s     %s\n' % (write_float_12E(dt), sgridtype, dx))
                elif sgridtype == 'L':
                    f.write('%14s %6s     %s\n' % (write_float_12E(dt), sgridtype, dx))
                else:
                    raise NotImplementedError(sgridtype)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort1_as_sort1(self, f, page_num, page_stamp, header, words):
        nodes = self.node_gridtype[:, 0]
        gridtypes = self.node_gridtype[:, 1]
        #times = self._times

        for itime in range(self.ntimes):
            dt = self._times[itime]
            t1 = self.data[itime, :, 0]
            if isinstance(dt, (float, float32)):
                header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            else:
                header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            f.write(''.join(header + words))
            for node_id, gridtypei, t1i in zip(nodes, gridtypes, t1):
                sgridtype = self.recast_gridtype_as_string(gridtypei)
                vals = [t1i]
                vals2 = write_floats_13e(vals)
                dx = vals2[0]
                if sgridtype == 'G':
                    f.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
                elif sgridtype == 'S':
                    f.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
                elif sgridtype == 'H':
                    f.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
                elif sgridtype == 'L':
                    f.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
                else:
                    raise NotImplementedError(sgridtype)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort2_as_sort2(self, f, page_num, page_stamp, header, words):
        nodes = self.node_gridtype[:, 0]
        gridtypes = self.node_gridtype[:, 1]
        times = self._times
        for inode, (node_id, gridtypei) in enumerate(zip(nodes, gridtypes)):
            t1 = self.data[inode, :, 0]

            header[1] = ' POINT-ID = %10i\n' % node_id
            f.write(''.join(header + words))
            for dt, t1i in zip(times, t1):
                sgridtype = self.recast_gridtype_as_string(gridtypei)
                vals = [t1i]
                vals2 = write_floats_13e(vals)
                dx = vals2[0]
                if sgridtype == 'G':
                    f.write('%14s %6s     %s\n' % (write_float_12E(dt), sgridtype, dx))
                elif sgridtype == 'S':
                    f.write('%14s %6s     %s\n' % (node_id, sgridtype, dx))
                elif sgridtype == 'H':
                    f.write('%14s %6s     %s\n' % (write_float_12E(dt), sgridtype, dx))
                elif sgridtype == 'L':
                    f.write('%14s %6s     %s\n' % (write_float_12E(dt), sgridtype, dx))
                else:
                    raise NotImplementedError(sgridtype)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_f06_transient_block(self, words, header, page_stamp, page_num, f, write_words,
                                   is_mag_phase=False, is_sort1=True):
        if write_words:
            words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()

        if not len(header) >= 3:
            header.append('')

        is_sort2 = not is_sort1
        if self.is_sort1() or self.nonlinear_factor is None:
            if is_sort2 and self.nonlinear_factor is not None:
                page_num = self._write_sort1_as_sort2(f, page_num, page_stamp, header, words)
            else:
                page_num = self._write_sort1_as_sort1(f, page_num, page_stamp, header, words)
        else:
            page_num = self._write_sort2_as_sort2(f, page_num, page_stamp, header, words)
        return page_num - 1

    def extract_xyplot(self, node_ids, index):
        node_ids = asarray(node_ids, dtype='int32')
        i = index - 1
        assert index in [1, 2, 3, 4, 5, 6], index
        nids = self.node_gridtype[:, 0]
        inids = searchsorted(nids, node_ids)
        assert all(nids[inids] == node_ids), 'nids=%s expected=%s; all=%s'  % (nids[inids], node_ids, nids)
        return self.data[:, inids, i]


#class ComplexScalarTableArray(TableArray):  # displacement style table
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #ScalarTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
        #raise NotImplementedError()

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
            #return abs(self.data[:, inids, i])
        #elif j == 4:
            ## phase
            #return angle(self.data[:, inids, i])
        #else:
            #raise RuntimeError()

    #def is_real(self):
        #return False

    #def is_complex(self):
        #return True

    #def data_type(self):
        #return 'complex64'

    ##def _write_f06_block(self, words, header, page_stamp, page_num, f, is_mag_phase):
        ##self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_mag_phase, is_sort1)

    #def _write_f06_transient_block(self, words, header, page_stamp, page_num, f,
                                   #is_mag_phase, is_sort1):
        #if is_mag_phase:
            #words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        #else:
            #words += ['                                                          (REAL/IMAGINARY)\n', ]

        #if not len(header) >= 3:
            #header.append('')

        #if self.is_sort1():
            #if is_sort1:
                #words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
                #page_num = self.write_sort1_as_sort1(f, page_num, page_stamp, header, words, is_mag_phase)
            #else:
                #words += [' \n', '      FREQUENCY   TYPE          T1             T2             T3             R1             R2             R3\n']
                #page_num = self.write_sort1_as_sort2(f, page_num, page_stamp, header, words, is_mag_phase)
        #else:
            #words += [' \n', '      FREQUENCY   TYPE          T1             T2             T3             R1             R2             R3\n']
            #page_num = self.write_sort2_as_sort2(f, page_num, page_stamp, header, words, is_mag_phase)
        #return page_num - 1

    #def write_sort1_as_sort1(self, f, page_num, page_stamp, header, words, is_mag_phase):
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
                    #f.write('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                            #'  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                #node_id, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                #'', '', dxi, dyi, dzi, rxi, ryi, rzi))
                #elif sgridtype == 'S':
                    #f.write('0 %12i %6s     %-13s\n'
                            #'  %12s %6s     %-13s\n' % (node_id, sgridtype, dxr, '', '', dxi))
                #else:
                    #raise NotImplementedError(sgridtype)
            #f.write(page_stamp % page_num)
            #page_num += 1
        #return page_num

    #def write_sort1_as_sort2(self, f, page_num, page_stamp, header, words, is_mag_phase):
        #node = self.node_gridtype[:, 0]
        #gridtype = self.node_gridtype[:, 1]

        #times = self._times
        ## print(self.data.shape)
        #for inode, (node_id, gridtypei) in enumerate(zip(node, gridtype)):
            ## TODO: for SORT1 pretending to be SORT2
            ##t1 = self.data[:, inode, 0].ravel()
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
                #sdt = write_float_12E(dt)
                ##if not is_all_zeros:
                #if sgridtype == 'G':
                    #f.write('0 %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                            #'  %13s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                #sdt, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                #'', '', dxi, dyi, dzi, rxi, ryi, rzi))
                #elif sgridtype == 'S':
                    #f.write('0 %12s %6s     %-13s\n'
                            #'  %12s %6s     %-13s\n' % (sdt, sgridtype, dxr, '', '', dxi))
                #else:
                    #msg = 'nid=%s dt=%s type=%s dx=%s dy=%s dz=%s rx=%s ry=%s rz=%s' % (
                    #node_id, dt, sgridtype, t1i, t2i, t3i, r1i, r2i, r3i)
                    #raise NotImplementedError(msg)
            #f.write(page_stamp % page_num)
            #page_num += 1
        #return page_num

    #def write_sort2_as_sort2(self, f, page_num, page_stamp, header, words, is_mag_phase):
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
                #sdt = write_float_12E(dt)
                ##if not is_all_zeros:
                #if sgridtype == 'G':
                    #f.write('0 %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                            #'  %13s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                #sdt, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                #'', '', dxi, dyi, dzi, rxi, ryi, rzi))
                #elif sgridtype == 'S':
                    #f.write('0 %12s %6s     %-13s\n'
                            #'  %12s %6s     %-13s\n' % (sdt, sgridtype, dxr, '', '', dxi))
                #else:
                    #msg = 'nid=%s dt=%s type=%s dx=%s dy=%s dz=%s rx=%s ry=%s rz=%s' % (
                    #node_id, dt, sgridtype, t1i, t2i, t3i, r1i, r2i, r3i)
                    #raise NotImplementedError(msg)
            #f.write(page_stamp % page_num)
            #page_num += 1
        #return page_num


#class StaticArrayNode(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def node_ids(self):
        #return self.node_gridtype[:, 0]


#class StaticArrayElement(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def element_ids(self):
        #return self.node_gridtype[:, 0]


#class TimeArrayNodeSort1(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def times(self):
        #return self._times
    #@property
    #def node_ids(self):
        #return self.node_gridtype[:, 0]


#class TimeArrayElementSort1(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def times(self):
        #return self._times
    #@property
    #def element_ids(self):
        #return self.node_gridtype[:, 0]


#class TimeArrayNodeSort2(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def times(self):
        #return self._times
    #@property
    #def node_ids(self):
        #return self.node_gridtype[:, 0]


#class TimeArrayElementSort2(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def times(self):
        #return self._times
    #@property
    #def element_ids(self):
        #return self.node_gridtype[:, 0]


#class FrequencyArrayNodeSort2(ComplexTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def frequencies(self):
        #return self._times
    #@property
    #def node_ids(self):
        #return self.node_gridtype[:, 0]


#class FrequencyArrayElementSort2(ComplexTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def frequencies(self):
        #return self._times
    #@property
    #def node_ids(self):
        #return self.node_gridtype[:, 0]
