"""
defines:
 - ScalarTableObject

"""
from struct import Struct, pack
import warnings
from typing import List

import numpy as np
from numpy import zeros, searchsorted, unique, where, float32
from numpy import allclose, asarray, vstack

from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.op2.result_objects.table_object import append_sort1_sort2
from pyNastran.f06.f06_formatting import write_floats_13e, write_float_12e
from pyNastran.op2.op2_interface.write_utils import set_table3_field

float_types = (float, np.float32)
integer_types = (int, np.int32)


SORT2_TABLE_NAME_MAP = {
    'OUGV2' : 'OUGV1',
    'OPG2' : 'OPG1',
}
class ScalarTableArray(ScalarObject):  # displacement style table
    def __init__(self, data_code, unused_is_sort1, isubcase, unused_dt):
        self.nonlinear_factor = np.nan
        self.table_name = None
        self.approach_code = None
        self.analysis_code = None
        ScalarObject.__init__(self, data_code, isubcase,
                              apply_data_code=True)  # no double inheritance
        self.is_sort1
        self.is_sort2
        #self.dt = dt

        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        self.ntotal = 0
        self._nnodes = 0  # result specific

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.node_gridtype, table.node_gridtype):
            assert self.node_gridtype.shape == table.node_gridtype.shape, 'shape=%s table.shape=%s' % (self.node_gridtype.shape, table.node_gridtype.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'nid_gridtype:\n'
            msg += 'gridtype.shape=%s table.gridtype.shape=%s\n' % (str(self.node_gridtype.shape), str(table.node_gridtype.shape))
            for (nid, grid_type), (nid2, grid_type2) in zip(self.node_gridtype, table.node_gridtype):
                msg += '(%s, %s)    (%s, %s)\n' % (nid, grid_type, nid2, grid_type2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
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
                '<%s>; table_name=%r\n' % (self.__class__.__name__, self.table_name),
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]
        #ngrids = len(self.gridTypes)
        if short:
            return self._get_stats_short()
        msg = []

        unused_ntimesi, ntotal = self.data.shape[:2]
        ntimes = len(self._times)
        nnodes = self.node_gridtype.shape[0]

        nmajor = self.ntimes
        nminor = self.ntotal
        if self.is_sort1:
            assert nmajor == ntimes, 'ntimes=%s expected=%s' % (nmajor, ntimes)
            assert nminor == ntotal, 'ntotal=%s expected=%s' % (nminor, nnodes)
        else:
            if not nmajor == nnodes:
                msgi = 'nnodes=%s expected=%s' % (nmajor, nnodes)
                warnings.warn(msgi)
                msg.append('  WARNING: ' + msgi + '\n')
            assert nminor == ntotal, 'ntotal=%s expected=%s' % (nminor, ntimes)

        msg.append('  isubcase = %s\n' % self.isubcase)
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%s nnodes=%s, table_name=%s\n'
                       % (self.__class__.__name__, ntimes, nnodes, self.table_name))
        else:
            msg.append('  type=%s nnodes=%s, table_name=%s\n'
                       % (self.__class__.__name__, nnodes, self.table_name))
        headers = ', '.join(self._get_headers())
        #msg.append('  data: [%s] shape=%s dtype=%s\n'
                   #% (headers, [int(i) for i in self.data.shape], self.data.dtype))
        msg.append('  data: [%s] shape=%s dtype=%s\n'
                   % (headers,
                      [int(i) for i in self.data.shape], self.data.dtype))
        msg.append('  node_gridtype.shape = %s\n' % str(self.node_gridtype.shape).replace('L', ''))
        #msg.append('  gridTypes\n  ')
        msg += self.get_data_code()
        return msg

    @property
    def headers(self):
        return ['t1']

    def _get_headers(self):
        return self.headers

    def get_headers(self) -> List[str]:
        return self._get_headers()

    def _reset_indices(self):
        self.itotal = 0

    def build(self):
        """sizes the vectorized attributes of the ScalarTableArray"""
        #print('_nnodes=%s ntimes=%s sort1?=%s ntotal=%s -> _nnodes=%s' % (
            #self._nnodes, self.ntimes, self.is_sort1,
            #self.ntotal, self._nnodes // self.ntimes))

        self._nnodes //= self.ntimes
        self.itime = 0
        self.itotal = 0
        self.is_built = True

        if self.is_sort1:
            ntimes = self.ntimes
            nnodes = self.ntotal
            ntotal = self.ntotal
            nx = ntimes
            ny = nnodes
            #print("SORT1 ntimes=%s nnodes=%s" % (ntimes, nnodes))
        elif self.is_sort2:
            # flip this to sort1
            ntimes = self.ntotal
            nnodes = self.ntimes
            ntotal = nnodes
            nx = ntimes
            ny = nnodes
            #print("***SORT2 ntotal=%s nnodes=%s ntimes=%s" % (ntotal, nnodes, ntimes))
        else:
            raise RuntimeError('expected sort1/sort2\n%s' % self.code_information())
        self.build_data(ntimes, nnodes, ntotal, nx, ny, self._times_dtype)

    def build_data(self, ntimes, nnodes, ntotal, nx, ny, float_fmt):
        """actually performs the build step"""
        self.ntimes = ntimes
        self._nnodes = nnodes
        self.ntotal = ntotal

        _times = zeros(ntimes, dtype=float_fmt)
        #self.types = array(self.nelements, dtype='|S1')
        node_gridtype = zeros((nnodes, 2), dtype='int32')

        #[t1]
        data = zeros((nx, ny, 1), self.data_type())
        if self.load_as_h5:
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.node_gridtype = group.create_dataset('node_gridtype', data=node_gridtype)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.node_gridtype = node_gridtype
            self.data = data

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        node_gridtype = [self.node_gridtype[:, 0], self.gridtype_str]
        ugridtype_str = unique(self.gridtype_str)

        if self.nonlinear_factor not in (None, np.nan):
            #Time             0.0       10.0
            #NodeID Type Item
            #1      S    t1     0.0  0.006239
            #2      S    t1     0.0  0.006239
            #3      S    t1     0.0  0.006239
            #4      S    t1     0.0  0.006239
            #5      S    t1     0.0  0.006239
            #6      S    t1     0.0  0.006239
            #7      S    t1     0.0  0.006239
            #8      S    t1     0.0  0.006239
            #99     S    t1     0.0  5.000000
            #
            #Time             0.0       10.0
            #NodeID Type Item
            #1      S    t1     0.0  0.006239
            #2      S    t1     0.0  0.006239
            #3      S    t1     0.0  0.006239
            #4      S    t1     0.0  0.006239
            #5      S    t1     0.0  0.006239
            #6      S    t1     0.0  0.006239
            #7      S    t1     0.0  0.006239
            #8      S    t1     0.0  0.006239
            #99     S    t1     0.0  5.000000

            column_names, column_values = self._build_dataframe_transient_header()

            columns = pd.MultiIndex.from_arrays(column_values, names=column_names)

            nheaders = len(headers)
            def m_hstack(array, n):
                if n == 1:
                    return array
                return np.hstack([array]*n)

            ntimes, nnodes = self.data.shape[:2]
            node_gridtype_item = [
                m_hstack(self.node_gridtype[:, 0], nheaders),
                m_hstack(self.gridtype_str, nheaders),
                m_hstack(headers, nnodes)
            ]
            names = ['NodeID', 'Type', 'Item']
            index = pd.MultiIndex.from_arrays(node_gridtype_item, names=names)

            A = self.data.reshape(ntimes, nnodes*nheaders).T
            data_frame = pd.DataFrame(A, columns=columns, index=index)

            if 0:  # pragma: no cover
                # <=24.2
                data_frame = pd.Panel(self.data, items=column_values,
                                      major_axis=node_gridtype, minor_axis=headers).to_frame()
                data_frame.columns.names = column_names
                data_frame.index.names = ['NodeID', 'Type', 'Item']
                print(data_frame)

                letter_dims = [
                    ('G', 1),
                    ('E', 1),
                    ('S', 1),
                    ('H', 1),
                    ('L', 1),
                ]
                cat_keys = []
                for (letter, unused_dim) in letter_dims:
                    if letter not in ugridtype_str:
                        continue
                    eig = data_frame.xs(letter, level=1)
                    cat_keys.append(eig)
                data_frame = pd.concat(cat_keys)
                #print(data_frame)
        else:
            #self.data_frame = pd.Panel(self.data[0, :, :],
                                       #major_axis=node_gridtype, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = ['Static']
            #self.data_frame.index.names = ['NodeID', 'Type', 'Item']

            df1 = pd.DataFrame(self.node_gridtype[:, 0])
            df1.columns = ['NodeID']
            df2 = pd.DataFrame(self.gridtype_str)
            df2.columns = ['Type']
            df3 = pd.DataFrame(self.data[0])
            df3.columns = headers
            data_frame = df1.join([df2, df3])
        self.data_frame = data_frame
        #print(self.data_frame)

    def finalize(self):
        """
        Calls any OP2 objects that need to do any post matrix calcs
        """
        self.set_as_sort1()
        gridtypes = self.node_gridtype[:, 1]
        nnodes = len(gridtypes)
        self.gridtype_str = np.chararray((nnodes), unicode=True)
        ugridtypes = unique(gridtypes)
        for ugridtype in ugridtypes:
            i = where(gridtypes == ugridtype)
            self.gridtype_str[i] = self.recast_gridtype_as_string(ugridtype)
        #del self.itotal, self.itime

    def set_as_sort1(self):
        """changes the table into SORT1"""
        #if not self.table_name != 'OQMRMS1':
            #return
        if self.is_sort1:
            return
        #print('table_name=%r' % self.table_name)
        try:
            analysis_method = self.analysis_method
        except AttributeError:
            print(self.code_information())
            raise
        #print(self.get_stats())
        #print(self.node_gridtype)
        #print(self.data.shape)
        self.sort_method = 1
        self.sort_bits[1] = 0
        bit0, bit1, bit2 = self.sort_bits
        self.table_name = SORT2_TABLE_NAME_MAP[self.table_name]
        self.sort_code = bit0 + 2*bit1 + 4*bit2
        #print(self.code_information())
        assert self.is_sort1
        if analysis_method != 'N/A':
            self.data_names[0] = analysis_method
            #print(self.table_name_str, analysis_method, self._times)
            setattr(self, self.analysis_method + 's', self._times)
        del self.analysis_method

    def add_sort1(self, dt, node_id, grid_type, v1):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(node_id, int) and node_id > 0, 'dt=%s node_id=%s' % (dt, node_id)
        # itotal - the node number
        # itime - the time/frequency step

        # the times/freqs
        self._times[self.itime] = dt
        self.node_gridtype[self.itotal, :] = [node_id, grid_type]
        self.data[self.itime, self.itotal, 0] = v1
        self.itotal += 1

    def add_sort2(self, dt, node_id, grid_type, v1):
        """unvectorized method for adding SORT2 transient data"""
        # itotal - the time/frequency step
        # itime - the node number
        self._times[self.itotal] = dt
        self.node_gridtype[self.itime, :] = [node_id, grid_type]
        self.data[self.itotal, self.itime, 0] = v1
        self.itotal += 1
        #self.itime += 1


class RealScalarTableArray(ScalarTableArray):  # temperature style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False

    def data_type(self):
        return 'float32'

    def _write_table_3(self, op2_file, fascii, new_result, itable=-3, itime=0):
        import inspect
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        fascii.write('%s.write_table_3: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if new_result and itable != -3:
            header = [
                4, 146, 4,
            ]
        else:
            header = [
                4, itable, 4,
                4, 1, 4,
                4, 0, 4,
                4, 146, 4,
            ]
        op2_file.write(pack(b'%ii' % len(header), *header))
        fascii.write('table_3_header = %s\n' % header)

        approach_code = self.approach_code
        table_code = self.table_code
        isubcase = self.isubcase
        random_code = self.random_code
        format_code = 1
        num_wide = self.num_wide
        acoustic_flag = 0
        thermal = self.thermal
        title = b'%-128s' % self.title.encode('ascii')
        subtitle = b'%-128s' % self.subtitle.encode('ascii')
        label = b'%-128s' % self.label.encode('ascii')
        ftable3 = b'i' * 50 + b'128s 128s 128s'
        oCode = 0

        field6 = 0
        field7 = 0
        if self.analysis_code == 1:
            field5 = self.lsdvmns[itime]
        #elif self.analysis_code == 2:
            #field5 = self.modes[itime]
            #field6 = self.eigns[itime]
            #field7 = self.mode_cycles[itime]
            #assert isinstance(field6, float), type(field6)
            #assert isinstance(field7, float), type(field7)
            #ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            #ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
        #elif self.analysis_code == 5:
            #field5 = self.freqs[itime]
            #ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 6:
            if hasattr(self, 'dts'):
                field5 = self.dts[itime]
            else:
                field5 = self.times[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        #elif self.analysis_code == 7:  # pre-buckling
            #field5 = self.lsdvmns[itime] # load set number
        #elif self.analysis_code == 8:  # post-buckling
            #field5 = self.lsdvmns[itime] # load set number
            #field6 = self.eigns[itime]
            #ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
        #elif self.analysis_code == 9:  # complex eigenvalues
            #field5 = self.modes[itime]
            #field6 = self.eigns[itime]
            #field7 = self.eigis[itime]
            #ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            #ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
        elif self.analysis_code == 10:  # nonlinear statics
            field5 = self.lftsfqs[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5; load step
        #elif self.analysis_code == 11:  # old geometric nonlinear statics
            #field5 = self.lsdvmns[itime] # load set number
        else:
            raise NotImplementedError(self.analysis_code)

        table3 = [
            approach_code, table_code, 0, isubcase, field5,
            field6, field7, random_code, format_code, num_wide,
            oCode, acoustic_flag, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, thermal, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            title, subtitle, label,
        ]
        #self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False) # 9
        #self.thermal = self.add_data_parameter(data, 'thermal', b'i', 23, False)
        assert table3[22] == thermal

        n = 0
        for v in table3:
            if isinstance(v, (integer_types, float_types)):
                n += 4
            else:
                n += len(v)
        assert n == 584, n
        data = [584] + table3 + [584]
        fmt = b'i' + ftable3 + b'i'
        #print(fmt)
        #op2_file.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        fascii.write('%s header 3c = %s\n' % (self.table_name, data))
        op2_file.write(pack(fmt, *data))

    def write_op2(self, op2_file, fascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        assert self.table_name in ['OPG1', 'OUGV1', 'TOUGV1'], self.table_name  # 'OUGV1', 'OQMG1', 'OQG1'

        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        fascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2_file, fascii, date)
            itable = -3

        if self.is_sort1:
            s = Struct(endian + b'2i6f')
        else:
            raise NotImplementedError('SORT2')

        unused_node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]
        #format_table4_1 = Struct(self._endian + b'15i')
        #format_table4_2 = Struct(self._endian + b'3i')

        # table 4 info
        #ntimes = self.data.shape[0]
        nnodes = self.data.shape[1]
        nnodes_device = self.node_gridtype[:, 0] * 10 + self.device_code

        #(2+6) => (node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i)
        ntotal = nnodes * (2 + 6)

        #print('shape = %s' % str(self.data.shape))
        assert nnodes > 1, nnodes
        assert ntotal > 1, ntotal

        unused_device_code = self.device_code
        fascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        for itime in range(self.ntimes):
            self._write_table_3(op2_file, fascii, new_result, itable, itime)

            # record 4
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4*ntotal]
            op2_file.write(pack(b'%ii' % len(header), *header))
            fascii.write('r4 [4, 0, 4]\n')
            fascii.write('r4 [4, %s, 4]\n' % (itable))
            fascii.write('r4 [4, %i, 4]\n' % (4*ntotal))

            t1 = self.data[itime, :, 0]
            for node_id, gridtypei, t1i in zip(nnodes_device, gridtype, t1):
                data = [node_id, gridtypei, t1i, 0., 0., 0., 0., 0.]
                fascii.write('  nid, grid_type, dx, dy, dz, rx, ry, rz = %s\n' % data)
                op2_file.write(s.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack(b'i', *header))
            fascii.write('footer = %s\n' % header)
            new_result = False
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

    def _write_f06_block(self, words, header, page_stamp, page_num, f06_file, write_words,
                         is_mag_phase=False, is_sort1=True):
        if write_words:
            words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()
        f06_file.write(''.join(header + words))

        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]
        t1 = self.data[0, :, 0]
        for node_id, gridtypei, t1i in zip(node, gridtype, t1):
            sgridtype = self.recast_gridtype_as_string(gridtypei)
            vals = [t1i]
            vals2 = write_floats_13e(vals)
            dx = vals2[0]
            f06_file.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
        f06_file.write(page_stamp % page_num)
        return page_num

    def _write_sort1_as_sort2(self, f06_file, page_num, page_stamp, header, words):
        nodes = self.node_gridtype[:, 0]
        gridtypes = self.node_gridtype[:, 1]
        times = self._times

        for inode, (node_id, gridtypei) in enumerate(zip(nodes, gridtypes)):
            t1 = self.data[:, inode, 0].ravel()

            header[1] = ' POINT-ID = %10i\n' % node_id
            f06_file.write(''.join(header + words))
            for dt, t1i in zip(times, t1):
                sgridtype = self.recast_gridtype_as_string(gridtypei)
                vals = [t1i]
                vals2 = write_floats_13e(vals)
                dx = vals2[0]
                if sgridtype == 'G':
                    f06_file.write('%14s %6s     %s\n' % (write_float_12e(dt), sgridtype, dx))
                elif sgridtype == 'S':
                    f06_file.write('%14s %6s     %s\n' % (node_id, sgridtype, dx))
                elif sgridtype == 'H':
                    f06_file.write('%14s %6s     %s\n' % (write_float_12e(dt), sgridtype, dx))
                elif sgridtype == 'L':
                    f06_file.write('%14s %6s     %s\n' % (write_float_12e(dt), sgridtype, dx))
                else:
                    raise NotImplementedError(sgridtype)
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort1_as_sort1(self, f06_file, page_num, page_stamp, header, words):
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
            f06_file.write(''.join(header + words))
            for node_id, gridtypei, t1i in zip(nodes, gridtypes, t1):
                sgridtype = self.recast_gridtype_as_string(gridtypei)
                vals = [t1i]
                vals2 = write_floats_13e(vals)
                dx = vals2[0]
                if sgridtype == 'G':
                    f06_file.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
                elif sgridtype == 'S':
                    f06_file.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
                elif sgridtype == 'H':
                    f06_file.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
                elif sgridtype == 'L':
                    f06_file.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
                else:
                    raise NotImplementedError(sgridtype)
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort2_as_sort2(self, f06_file, page_num, page_stamp, header, words):
        nodes = self.node_gridtype[:, 0]
        gridtypes = self.node_gridtype[:, 1]
        times = self._times
        for inode, (node_id, gridtypei) in enumerate(zip(nodes, gridtypes)):
            t1 = self.data[inode, :, 0]

            header[1] = ' POINT-ID = %10i\n' % node_id
            f06_file.write(''.join(header + words))
            for dt, t1i in zip(times, t1):
                sgridtype = self.recast_gridtype_as_string(gridtypei)
                vals = [t1i]
                vals2 = write_floats_13e(vals)
                dx = vals2[0]
                if sgridtype == 'G':
                    f06_file.write('%14s %6s     %s\n' % (write_float_12e(dt), sgridtype, dx))
                elif sgridtype == 'S':
                    f06_file.write('%14s %6s     %s\n' % (node_id, sgridtype, dx))
                elif sgridtype == 'H':
                    f06_file.write('%14s %6s     %s\n' % (write_float_12e(dt), sgridtype, dx))
                elif sgridtype == 'L':
                    f06_file.write('%14s %6s     %s\n' % (write_float_12e(dt), sgridtype, dx))
                else:
                    raise NotImplementedError(sgridtype)
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_f06_transient_block(self, words, header, page_stamp, page_num, f06_file, write_words,
                                   is_mag_phase=False, is_sort1=True):
        if write_words:
            words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
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

    def extract_xyplot(self, node_ids, index):
        node_ids = asarray(node_ids, dtype='int32')
        i = index - 1
        assert index in [1, 2, 3, 4, 5, 6], index
        nids = self.node_gridtype[:, 0]
        inids = searchsorted(nids, node_ids)
        assert all(nids[inids] == node_ids), 'nids=%s expected=%s; all=%s'  % (nids[inids], node_ids, nids)
        return self.data[:, inids, i]
