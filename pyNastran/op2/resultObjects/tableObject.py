from __future__ import print_function
from six import iteritems
from six.moves import zip, range
from struct import Struct, pack

from numpy import array, zeros, abs, angle, float32, ndarray, searchsorted, asarray

from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.f06.f06_formatting import writeFloats13E, writeImagFloats13E, write_float_12E


class TableArray(ScalarObject):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.nonlinear_factor = None
        self.table_name = None
        self.analysis_code = None
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)  # no double inheritance
        #self.dt = dt

        #self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        self.ntotal = 0
        self._nnodes = 0  # result specific

    def _get_msgs(self, is_mag_phase):
        raise NotImplementedError()

    def data_type(self):
        raise NotImplementedError()

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
                    '  ntimes: %i\n' % self.ntimes,
                    '  ntotal: %i\n' % self.ntotal,
                    ]
        #ngrids = len(self.gridTypes)
        msg = []

        ntimes = len(self._times)
        #len(self.node_gridtype)
        nnodes = self.node_gridtype.shape[0]
        ntimes, ntotal = self.data.shape[:2]
        assert self.ntimes == ntimes, 'ntimes=%s expected=%s' % (self.ntimes, ntimes)
        assert self.ntotal == ntotal, 'ntotal=%s expected=%s' % (self.ntimes, ntimes)

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%s nnodes=%s\n'
                       % (self.__class__.__name__, ntimes, nnodes))
        else:
            msg.append('  type=%s nnodes=%s\n'
                       % (self.__class__.__name__, nnodes))
        msg.append('  data: [t1, t2, t3, r1, r2, r3] shape=%s dtype=%s\n'
                   % ([int(i) for i in self.data.shape], self.data.dtype))
        msg.append('  gridTypes\n  ')
        msg += self.get_data_code()
        return msg

    def _reset_indices(self):
        self.itotal = 0

    def build(self):
        if self.is_built:
            #print("resetting...")
            #self.itotal = 0
            return

        self._nnodes //= self.ntimes
        self.itime = 0
        self.itotal = 0
        self.is_built = True

        #print("ntimes=%s nnodes/elements=%s" % (self.ntimes, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        #self.types = array(self.nelements, dtype='|S1')

        self.node_gridtype = zeros((self.ntotal, 2), dtype='int32')

        #[t1, t2, t3, r1, r2, r3]
        self.data = zeros((self.ntimes, self.ntotal, 6), self.data_type())
        #print(str(self))

    def add(self, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        self.add_sort1(None, node_id, grid_type, v1, v2, v3, v4, v5, v6)

    def add_sort1(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        #print "dt=%s out=%s" %(dt,out)
        #if dt not in self.translations:
        #    self.add_new_transient(dt)
        msg = "node_id=%s v1=%s v2=%s v3=%s\n" % (node_id, v1, v2, v3)
        msg += "          v4=%s v5=%s v6=%s" % (v4, v5, v6)
        #print(msg)
        #assert node_id == 1575, msg
        assert -1 < node_id < 1000000000, msg
        assert isinstance(node_id, int), node_id
        #assert isinstance(node_id, int), msg
        #assert node_id not in self.translations[self.dt],'displacementObject - transient failure'

        #self.gridTypes[node_id] = self.recastGridType(grid_type)

        # [t1, t2, t3, r1, r2, r3]
        #print "%s node_gridtype[%s, :] = %s" % (self.__class__.__name__, self.itotal, [node_id, grid_type]),
        #print "%s data[%s, %s, :] = %s" % (self.__class__.__name__, self.itime, self.itotal, [v1, v2, v3, v4, v5, v6])

        # itotal - the node number
        # itime - the time/frequency step

        # the times/freqs
        self._times[self.itime] = dt
        self.node_gridtype[self.itotal, :] = [node_id, grid_type]
        self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]
        self.itotal += 1

    def add_sort2(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        #print "dt=%s out=%s" %(dt,out)
        #if dt not in self.translations:
        #    self.add_new_transient(dt)
        print('itime=%s itotal=%s' % (self.itotal, self.itime))
        msg = "dt=%s node_id=%s v1=%s v2=%s v3=%s\n" % (dt, node_id, v1, v2, v3)
        msg += "                    v4=%s v5=%s v6=%s" % (v4, v5, v6)
        #print(msg)
        #assert node_id == 1575, msg
        assert -1 < node_id < 1000000000, msg
        assert isinstance(node_id, int), node_id
        #assert isinstance(node_id, int), msg
        #assert node_id not in self.translations[self.dt],'displacementObject - transient failure'

        #self.gridTypes[node_id] = self.recastGridType(grid_type)

        # [t1, t2, t3, r1, r2, r3]
        #print "%s node_gridtype[%s, :] = %s" % (self.__class__.__name__, self.itotal, [node_id, grid_type]),
        #print "%s data[%s, %s, :] = %s" % (self.__class__.__name__, self.itime, self.itotal, [v1, v2, v3, v4, v5, v6])

        # the node IDs
        self._times[self.itime] = dt
        self.node_gridtype[self.itime, :] = [node_id, grid_type]
        self.data[self.itotal, self.itime, :] = [v1, v2, v3, v4, v5, v6]

        # itotal - the node number
        # itime - the time/frequency step
        self.itotal += 1
        #self.itime += 1


class RealTableArray(TableArray):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        TableArray.__init__(self, data_code, is_sort1, isubcase, dt)

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
        Title = '%-128s' % self.Title
        subtitle = '%-128s' % self.subtitle
        label = '%-128s' % self.label
        ftable3 = '50i 128s 128s 128s'
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
            Title, subtitle, label,
        ]

        n = 0
        for v in table3:
            if isinstance(v, int) or isinstance(v, float):
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
        #format_table4_1 = Struct(self._endian + b'15i')
        #format_table4_2 = Struct(self._endian + b'3i')

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
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(nnodes_device, gridtype, t1, t2, t3, r1, r2, r3):
                data = [node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i]
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

    def _write_f06_block(self, words, header, page_stamp, page_num, f, write_words=True):
        if write_words:
            words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()
        f.write(''.join(header + words))

        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]
        t1 = self.data[0, :, 0]
        t2 = self.data[0, :, 1]
        t3 = self.data[0, :, 2]
        r1 = self.data[0, :, 3]
        r2 = self.data[0, :, 4]
        r3 = self.data[0, :, 5]
        for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
            sgridtype = self.recast_gridtype_as_string(gridtypei)
            vals = [t1i, t2i, t3i, r1i, r2i, r3i]
            (vals2, is_all_zeros) = writeFloats13E(vals)
            (dx, dy, dz, rx, ry, rz) = vals2
            f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, sgridtype, dx, dy, dz, rx, ry, rz))
        f.write(page_stamp % page_num)
        return page_num

    def _write_f06_transient_block(self, words, header, page_stamp, page_num, f, write_words=True):
        if write_words:
            words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()

        if not len(header) >= 3:
            header.append('')

        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]
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
            f.write(''.join(header + words))
            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
                sgridtype = self.recast_gridtype_as_string(gridtypei)
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                (dx, dy, dz, rx, ry, rz) = vals2
                if sgridtype == 'G':
                    f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, sgridtype, dx, dy, dz, rx, ry, rz))
                elif sgridtype == 'S':
                    f.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
                elif sgridtype == 'H':
                    f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, sgridtype, dx, dy, dz, rx, ry, rz))
                elif sgridtype == 'L':
                    f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, sgridtype, dx, dy, dz, rx, ry, rz))
                else:
                    raise NotImplementedError(sgridtype)

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def extract_xyplot(self, node_ids, index):
        node_ids = asarray(node_ids, dtype='int32')
        i = index - 1
        assert index in [1, 2, 3, 4, 5, 6], index
        nids = self.node_gridtype[:, 0]
        inids = searchsorted(nids, node_ids)
        assert all(nids[inids] == node_ids), 'nids=%s expected=%s; all=%s'  % (nids[inids], node_ids, nids)
        return self.data[:, inids, i]


class ComplexTableArray(TableArray):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        TableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def extract_xyplot(self, node_ids, index, index_str):
        index_str = index_str.lower().strip()
        if index_str in ['real', 'r']:
            j = 1
        elif index_str in ['imag', 'i']:
            j = 2
        elif index_str in ['mag', 'magnitude', 'm']:
            j = 3
        elif index_str in ['phase', 'p']:
            j = 4
        else:
            raise ValueError('index_str=%r' % index_str)

        node_ids = asarray(node_ids, dtype='int32')
        i = index - 1
        assert index in [1, 2, 3, 4, 5, 6,
                         7, 8, 9, 10, 11, 12], index
        nids = self.node_gridtype[:, 0]
        inids = searchsorted(nids, node_ids)
        assert all(nids[inids] == node_ids), 'nids=%s expected=%s; all=%s'  % (nids[inids], node_ids, nids)
        if j == 1:
            # real
            return self.data[:, inids, i].real
        elif j == 2:
            # imag
            return self.data[:, inids, i].imag
        elif j == 3:
            # mag
            return abs(self.data[:, inids, i])
        elif j == 4:
            # phase
            return angle(self.data[:, inids, i])
        else:
            raise RuntimeError()

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def data_type(self):
        return 'complex64'

    #def _write_f06_block(self, words, header, page_stamp, page_num, f, is_mag_phase):
        #self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_mag_phase, is_sort2)

    def _write_f06_transient_block(self, words, header, page_stamp, page_num, f, is_mag_phase, is_sort2):
        if is_mag_phase:
            words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        else:
            words += ['                                                          (REAL/IMAGINARY)\n', ]

        words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()

        if not len(header) >= 3:
            header.append('')

        #is_sort1 = not is_sort2
        is_sort1 = True
        if is_sort1:
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
                header[2] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
                f.write(''.join(header + words))
                for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
                    sgridtype = self.recast_gridtype_as_string(gridtypei)
                    vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                    (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
                    [dxr, dyr, dzr, rxr, ryr, rzr,
                     dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                    if sgridtype == 'G':
                        f.write('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                                '  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                    node_id, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                    '', '', dxi, dyi, dzi, rxi, ryi, rzi))
                    elif sgridtype == 'S':
                        f.write('0 %12i %6s     %-s\n'
                                '  %12s %6s     %-s\n' % (node_id, sgridtype, dxr, '', '', dxi))
                    else:
                        raise NotImplementedError(sgridtype)
                f.write(page_stamp % page_num)
                page_num += 1
        else:
            node = self.node_gridtype[:, 0]
            gridtype = self.node_gridtype[:, 1]

            for inode, (node_id, gridtypei) in enumerate(zip(node, gridtype)):
                times = self._times
                t1 = self.data[:, inode, 0].ravel()
                t2 = self.data[:, inode, 1].ravel()
                t3 = self.data[:, inode, 2].ravel()
                r1 = self.data[:, inode, 3].ravel()
                r2 = self.data[:, inode, 4].ravel()
                r3 = self.data[:, inode, 5].ravel()
                if len(r3) != len(times):
                    raise RuntimeError('len(d)=%s len(times)=%s' % (len(r3), len(times)))

                header[2] = ' %s = %10i\n' % ('POINT-ID', node_id)
                f.write(''.join(header + words))
                for dt, t1i, t2i, t3i, r1i, r2i, r3i in zip(times, t1, t2, t3, r1, r2, r3):
                    sgridtype = self.recast_gridtype_as_string(gridtypei)
                    vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                    (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
                    [dxr, dyr, dzr, rxr, ryr, rzr,
                     dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                    sdt = write_float_12E(dt)
                    #if not is_all_zeros:
                    if sgridtype == 'G':
                        f.write('0 %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                                '  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                    sdt, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                    '', '', dxi, dyi, dzi, rxi, ryi, rzi))
                    elif sgridtype == 'S':
                        f.write('0 %12s %6s     %-s\n'
                                '  %12s %6s     %-s\n' % (sdt, sgridtype, dxr, '', '', dxi))
                    else:
                        raise NotImplementedError(sgridtype)
                f.write(page_stamp % page_num)
                page_num += 1
        return page_num - 1


class RealTableObject(ScalarObject):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.nonlinear_factor = None
        self.table_name = None
        self.analysis_code = None
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        self.gridTypes = {}
        self.translations = {}
        self.rotations = {}

        self.dt = dt

    def get_stats(self):
        ngrids = len(self.gridTypes)
        msg = []
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.translations)
            msg.append('  type=%s ntimes=%s ngrids=%s\n'
                       % (self.__class__.__name__, ntimes, ngrids))
        else:
            msg.append('  type=%s ngrids=%s\n'
                       % (self.__class__.__name__, ngrids))
        msg.append('  translations, rotations, gridTypes\n  ')
        msg += self.get_data_code()
        return msg

    def isImaginary(self):
        return False

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                node_id, grid_type, t1 = line[:3]
                if grid_type == 'G':
                    (t2, t3, r1, r2, r3) = line[3:]
                    self.gridTypes[node_id] = grid_type
                    self.translations[node_id] = array([t1, t2, t3], dtype='float32')
                    self.rotations[node_id] = array([r1, r2, r3], dtype='float32')
                elif grid_type == 'L':
                    # .. todo:: are L points single DOFs?
                    # we do know they're greater than the max value...
                    try:
                        (node_id, grid_type, t1, t2, t3, r1, r2, r3) = line
                    except ValueError:
                        print('L point line=%r' % line)
                        raise
                    self.gridTypes[node_id] = grid_type
                    self.translations[node_id] = array([t1, t2, t3], dtype='float32')
                    self.rotations[node_id] = array([r1, r2, r3], dtype='float32')
                elif grid_type == 'S':
                    t2 = t3 = r1 = r2 = r3 = 0.0
                    self.gridTypes[node_id] = grid_type
                    self.translations[node_id] = array([t1, t2, t3], dtype='float32')
                    self.rotations[node_id] = array([r1, r2, r3], dtype='float32')
                    for i, t1 in enumerate(line[3:]):
                        node_id2 = node_id + i + 1
                        self.gridTypes[node_id2] = grid_type
                        self.translations[node_id2] = array([t1, t2, t3], dtype='float32')
                        self.rotations[node_id2] = array([r1, r2, r3], dtype='float32')
                else:
                    raise NotImplementedError(line)
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.translations:
            self.update_dt(self.data_code, dt)

        for line in data:
            (node_id, grid_type, t1, t2, t3, r1, r2, r3) = line
            self.gridTypes[node_id] = grid_type
            self.translations[dt][node_id] = array([t1, t2, t3], dtype='float32')
            self.rotations[dt][node_id] = array([r1, r2, r3], dtype='float32')

    def update_dt(self, data_code, dt):
        self.data_code = data_code
        self.apply_data_code()
        if dt is not None:
            #print("updating %s...%s=%s  isubcase=%s"
            #      % (self.data_code['name'], self.data_code['name'],
            #         dt, self.isubcase))
            self.dt = dt
            self.add_new_transient(dt)

    def delete_transient(self, dt):
        del self.translations[dt]
        del self.rotations[dt]

    def get_transients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """initializes the transient variables"""
        self.dt = dt
        self.translations[dt] = {}
        self.rotations[dt] = {}

    def add(self, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        msg = "node_id=%s gridType=%s v1=%s v2=%s v3=%s" % (
            node_id, grid_type, v1, v2, v3)
        #print(msg)
        assert -1 < node_id < 1000000000, msg
        assert isinstance(node_id, int), msg
        #assert node_id not in self.translations,'displacementObject - static failure'

        self.gridTypes[node_id] = self.recast_gridtype_as_string(grid_type)
        self.translations[node_id] = array([v1, v2, v3], dtype='float32')  # dx,dy,dz
        self.rotations[node_id] = array([v4, v5, v6], dtype='float32')  # rx,ry,rz

    def add_sort1(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        msg = "node_id=%-8i v1=%s v2=%s v3=%s\n" % (node_id, v1, v2, v3)
        msg += "                 v4=%s v5=%s v6=%s" % (v4, v5, v6)
        #print(msg)
        if dt not in self.translations:
            self.add_new_transient(dt)
        assert -1 < node_id < 1000000000, msg
        assert dt is not None
        #assert isinstance(node_id,int),msg
        #assert node_id not in self.translations[self.dt],'displacementObject - transient failure'

        self.gridTypes[node_id] = self.recast_gridtype_as_string(grid_type)
        self.translations[dt][node_id] = array([v1, v2, v3], dtype='float32')  # dx,dy,dz
        self.rotations[dt][node_id] = array([v4, v5, v6], dtype='float32')  # rx,ry,rz

    def add_sort2(self, node_id, dt, grid_type, v1, v2, v3, v4, v5, v6):
        raise NotImplementedError('sort2')
        if dt not in self.translations:
            self.add_new_transient(dt)
        msg = "node_id=%s v1=%s v2=%s v3=%s\n" % (node_id, v1, v2, v3)
        msg += "          v4=%s v5=%s v6=%s" % (v4, v5, v6)
        msg = 'dt=%s node_id=%s' % (dt, node_id)
        #print msg
        assert 0 < node_id < 1000000000, msg    # remove
        assert dt is not None
        assert isinstance(node_id, int), msg  # remove
        assert -0.5 < dt, msg  # remove
        #assert node_id not in self.translations[self.dt],'displacementObject - transient failure'

        self.gridTypes[node_id] = self.recast_gridtype_as_string(grid_type)
        self.translations[dt][node_id] = array([v1, v2, v3], dtype='float32')  # dx,dy,dz
        self.rotations[dt][node_id] = array([v4, v5, v6], dtype='float32')  # rx,ry,rz

    def get_as_sort1(self):
        return (self.translations, self.rotations)

    def get_as_sort2(self):
        """returns translations and rotations in sort2 format"""
        translations2 = {}
        rotations2 = {}
        if self.dt is not None:
            for dt, translations in sorted(iteritems(self.translations)):
                node_ids = translations.keys()
                for node_id in node_ids:
                    translations2[node_id] = {}
                    rotations2[node_id] = {}

            for dt, translations in sorted(iteritems(self.translations)):
                for node_id, translation in sorted(iteritems(translations)):
                    rotation = self.rotations[dt][node_id]
                    translations2[node_id][dt] = translation
                    rotations2[node_id][dt] = rotation
        else:
            return (self.translations, self.rotations)
            #for node_id, translation in sorted(iteritems(self.translations)):
            #    rotation = self.rotations[node_id]
            #    translations2[node_id] = translation
            #    rotations2[node_id]    = rotation
        return (translations2, rotations2)

    def _write_f06_block(self, words, header, page_stamp, page_num=1, f=None):
        msg = words
        #assert f is not None # remove
        for node_id, translation in sorted(iteritems(self.translations)):
            rotation = self.rotations[node_id]
            grid_type = self.gridTypes[node_id]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, is_all_zeros) = writeFloats13E(vals)
            #if not is_all_zeros:
            (dx, dy, dz, rx, ry, rz) = vals2
            if grid_type == 'G':
                msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                           % (node_id, grid_type, dx, dy, dz, rx, ry, rz.rstrip()))
            elif grid_type == 'S':
                msg.append('%14i %6s     %s\n'
                           % (node_id, grid_type, dx.rstrip()))
            elif grid_type == 'H':
                msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                           % (node_id, grid_type, dx, dy, dz, rx, ry, rz.rstrip()))
            elif grid_type == 'L':
                msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                           % (node_id, grid_type, dx, dy, dz, rx, ry, rz.rstrip()))
            else:
                raise NotImplementedError(grid_type)
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient_block(self, words, header, page_stamp, page_num=1, f=None):
        msg = []
        #assert f is not None # remove
        for dt, translations in sorted(iteritems(self.translations)):
            if isinstance(dt, (float, float32)):  # fix
                #header[1] = ' %s = %10.4E float %s\n' % (self.data_code['name'], dt, self.analysis_code)
                header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            else:
                #header[1] = ' %s = %10i integer %s\n' % (self.data_code['name'], dt, self.analysis_code)
                header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            msg += header + words
            #if isinstance(translations, ndarray):
                # caused by duplicate table -> preload static solution before transient
                #F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_motion50t.op2
            for node_id, translation in sorted(iteritems(translations)):
                rotation = self.rotations[dt][node_id]
                grid_type = self.gridTypes[node_id]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation
                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                #if not is_all_zeros:
                [dx, dy, dz, rx, ry, rz] = vals2
                if grid_type == 'G':
                    msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, grid_type, dx, dy, dz, rx, ry, rz))
                elif grid_type == 'S':
                    msg.append('%14i %6s     %s\n' % (node_id, grid_type, dx))
                elif grid_type == 'L':
                    msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, grid_type, dx, dy, dz, rx, ry, rz))
                else:
                    raise NotImplementedError(grid_type)

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1

    def get_table_marker(self):
        if self.isATO():
            words = self.ATO_words()
        elif self.isCRM():
            words = self.CRM_words()
        elif self.isPSD():
            words = self.PSD_words()
        elif self.isRMS():
            words = self.RMS_words()
        elif self.isZERO():
            return self.ZERO_words()
        else:
            words = ['']
        return words

    def isATO(self):
        """Auto-Correlation Function"""
        if b'ATO' in self.table_name:
            return True
        return False

    def isCRM(self):
        """Correlated Root-Mean Square"""
        if b'CRM' in self.table_name:
            return True
        return False

    def isPSD(self):
        """Power Spectral Density"""
        if b'PSD' in self.table_name:
            return True
        return False

    def isRMS(self):
        """Root-Mean Square"""
        if b'RMS' in self.table_name:
            return True
        return False

    def isZERO(self):
        """Zero Crossings"""
        if b'NO' in self.table_name:
            return True
        return False

    def ATO_words(self):
        words = ['                                                 ( AUTO-CORRELATION FUNCTION )\n', ' \n']
        return words

    def CRM_words(self):
        words = ['                                               ( CUMULATIVE ROOT MEAN SQUARE )\n', ' \n']
        return words

    def PSD_words(self):
        words = ['                                             ( POWER SPECTRAL DENSITY FUNCTION )\n', ' \n']
        return words

    def RMS_words(self):
        words = ['                                                     ( ROOT MEAN SQUARE )\n', ' \n']
        return words

    def ZERO_words(self):
        words = ['                                                 ( NUMBER OF ZERO CROSSINGS )\n', ' \n']
        return words


class ComplexTableObject(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt, apply_data_code=True):
        self.nonlinear_factor = None
        self.table_name = None
        self.analysis_code = None
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=apply_data_code)
        self.gridTypes = {}
        self.translations = {}
        self.rotations = {}
        self.dt = dt

    def get_stats(self):
        ngrids = len(self.gridTypes)
        msg = []
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.translations)
            msg.append('  imaginary type=%s ntimes=%s ngrids=%s\n'
                       % (self.__class__.__name__, ntimes, ngrids))
        else:
            msg.append('  imaginary type=%s ngrids=%s\n'
                       % (self.__class__.__name__, ngrids))
        msg.append('  translations, rotations, gridTypes\n  ')
        msg += self.get_data_code()
        return msg

    def isImaginary(self):
        return True

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                try:
                    (node_id, grid_type, v1, v2, v3, v4, v5, v6) = line[2:]
                except:
                    print('line = %r' % line)
                    raise
                self.gridTypes[node_id] = grid_type
                self.translations[self.dt][node_id] = array([v1, v2, v3], dtype='complex64')  # dx,dy,dz
                self.rotations[self.dt][node_id] = array([v4, v5, v6], dtype='complex64')  # rx,ry,rz
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.translations:
            self.update_dt(self.data_code, dt)

        for line in data:
            try:
                (node_id, grid_type, v1, v2, v3, v4, v5, v6) = line
            except:
                print('line = %r' % line)
                raise
            #print "*dt=%s line=%s" % (self.dt, str(line))
            self.gridTypes[node_id] = grid_type
            self.translations[dt][node_id] = array([v1, v2, v3], dtype='complex64')  # dx,dy,dz
            self.rotations[dt][node_id] = array([v4, v5, v6], dtype='complex64')  # rx,ry,rz

    def add_complex_f06_data(self, data, transient):
        raise NotImplementedError()

    def update_dt(self, data_code, dt):
        self.data_code = data_code
        self.apply_data_code()
        if dt is not None:
            #print("updating %s...%s=%s  isubcase=%s" % (self.data_code['name'], self.data_code['name'], dt, self.isubcase))
            self.add_new_transient(dt)

    def delete_transient(self, dt):
        del self.translations[dt]
        del self.rotations[dt]

    def get_transients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """initializes the transient variables"""
        self.dt = dt
        self.translations[dt] = {}
        self.rotations[dt] = {}

    def add(self, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        #msg = "node_id=%s v1=%s v2=%s v3=%s" %(node_id, v1, v2, v3)
        #assert isinstance(node_id, int), node_id
        msg = "node_id=%s v1=%s v2=%s v3=%s\n" % (node_id, v1, v2, v3)
        msg += "          v4=%s v5=%s v6=%s" % (v4, v5, v6)
        #print(msg)
        assert 0 < node_id < 1000000000, msg  # -1
        assert -0.5 < dt, msg              # remove
        assert isinstance(node_id, int), msg
        #assert node_id not in self.translations,'complexDisplacementObject - static failure'

        self.gridTypes[node_id] = self.recast_gridtype_as_string(grid_type)
        self.translations[node_id] = array([v1, v2, v3], dtype='complex64')  # dx,dy,dz
        self.rotations[node_id] = array([v4, v5, v6], dtype='complex64')  # rx,ry,rz

    def add_sort1(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        #msg = "dt=%s node_id=%s v1=%s v2=%s v3=%s" %(dt, node_id, v1, v2, v3)
        #print(msg)
        if dt not in self.translations:
            self.add_new_transient(dt)
        #print(msg)
        if isinstance(node_id, float):
            msg = "node_id=%s v1=%s v2=%s v3=%s\n" % (node_id, v1, v2, v3)
            msg += "          v4=%s v5=%s v6=%s" % (v4, v5, v6)
            print(msg)

        if not 0 < node_id < 1000000000:
            msg = "node_id=%s v1=%s v2=%s v3=%s\n" % (node_id, v1, v2, v3)
            msg += "          v4=%s v5=%s v6=%s" % (v4, v5, v6)
            raise ValueError(msg)
        #assert -0.5 < dt, msg  # TODO: remove
        #assert isinstance(node_id, int), msg  # TODO: remove
        #assert node_id not in self.translations,'complexDisplacementObject - static failure'

        self.gridTypes[node_id] = self.recast_gridtype_as_string(grid_type)
        self.translations[dt][node_id] = array([v1, v2, v3], dtype='complex64')  # dx,dy,dz
        self.rotations[dt][node_id] = array([v4, v5, v6], dtype='complex64')  # rx,ry,rz

    def add_sort2(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        #[dt, grid_type, v1, v2, v3, v4, v5, v6] = data

        if dt not in self.translations:
            self.add_new_transient(dt)

        msg = "dt=%s node_id=%s v1=%s v2=%s v3=%s\n" % (dt, node_id, v1, v2, v3)
        msg += "                v4=%s v5=%s v6=%s" % (v4, v5, v6)
        #print msg
        assert 0 < node_id < 1000000000, msg  # -1
        assert -0.5 < dt, msg  # remove
        assert isinstance(node_id, int), msg

        self.gridTypes[node_id] = self.recast_gridtype_as_string(grid_type)
        self.translations[dt][node_id] = array([v1, v2, v3], dtype='complex64')  # dx,dy,dz
        self.rotations[dt][node_id] = array([v4, v5, v6], dtype='complex64')  # rx,ry,rz

    def _write_f06_block(self, words, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort2=False):
        raise RuntimeError('is this function used???')
        #words += self.getTableMarker()
        if is_mag_phase:
            words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        else:
            words += ['                                                          (REAL/IMAGINARY)\n', ]

        words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']

        msg = words
        for node_id, translation in sorted(iteritems(self.translations)):
            rotation = self.rotations[node_id]
            grid_type = self.gridTypes[node_id]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation

            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
            [dxr, dyr, dzr, rxr, ryr, rzr,
             dxi, dyi, dzi, rxi, ryi, rzi] = vals2
            if grid_type == 'G':
                f.write('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                        '  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                          % (node_id, grid_type, dxr, dyr, dzr, rxr, ryr, rzr,
                                        '', '', dxi, dyi, dzi, rxi, ryi, rzi))
            elif grid_type == 'S':
                f.write('0 %12i %6s     %-13s\n'
                        '  %12s %6s     %-13s\n' % (node_id, grid_type, dxr, '', '', dxi))
            else:
                raise NotImplementedError(grid_type)
        f.write(page_stamp % page_num)
        msg = ['']
        return page_num

    def _write_f06_transient_block(self, words, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort2=False):
        if is_mag_phase:
            words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        else:
            words += ['                                                          (REAL/IMAGINARY)\n', ]

        words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()

        msg = []
        if not len(header) >= 3:
            header.append('')

        for dt, translations in sorted(iteritems(self.translations)):
            header[2] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for node_id, translation in sorted(iteritems(translations)):
                rotation = self.rotations[dt][node_id]
                grid_type = self.gridTypes[node_id]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
                [dxr, dyr, dzr, rxr, ryr, rzr, dxi, dyi,
                    dzi, rxi, ryi, rzi] = vals2
                #if not is_all_zeros:
                if grid_type == 'G':
                    msg.append('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, grid_type, dxr, dyr, dzr, rxr, ryr, rzr))
                    msg.append('  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % ('', '', dxi, dyi, dzi, rxi, ryi, rzi))
                elif grid_type == 'S':
                    msg.append('0 %12i %6s     %s\n' % (node_id, grid_type, dxr))
                    msg.append('  %12s %6s     %s\n' % ('', '', dxi))
                else:
                    raise NotImplementedError(grid_type)
            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1
