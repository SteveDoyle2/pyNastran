from six import string_types, iteritems
from six.moves import zip, range
from struct import Struct, pack

from numpy import array, zeros, sqrt, abs, angle  # dot,

from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.f06.f06_formatting import writeFloats13E, writeImagFloats13E


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
        nnodes, two = self.node_gridtype.shape
        ntimes, ntotal, six = self.data.shape
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

    def add(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        self.add_sort1(dt, node_id, grid_type, v1, v2, v3, v4, v5, v6)

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
        #assert isinstance(nodeID, int), msg
        #assert nodeID not in self.translations[self.dt],'displacementObject - transient failure'

        #self.gridTypes[nodeID] = self.recastGridType(grid_type)

        # [t1, t2, t3, r1, r2, r3]
        #print "%s node_gridtype[%s, :] = %s" % (self.__class__.__name__, self.itotal, [node_id, grid_type]),
        #print "%s data[%s, %s, :] = %s" % (self.__class__.__name__, self.itime, self.itotal, [v1, v2, v3, v4, v5, v6])
        self._times[self.itime] = dt
        self.node_gridtype[self.itotal, :] = [node_id, grid_type]
        self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]
        self.itotal += 1


class RealTableArray(TableArray):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        TableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def data_type(self):
        return 'float32'

    def _write_op2_header(self, f, table_num, i):
        # record 3
        header1 = [  #4, 0, 4,
                  #4, -n, 4,
                  #4, 0, 4,

                  4, 146, 4,
                  4, 584, 4]
        header_format1 = b'6i '

        # end of record 3
        header2 = [4, 584, 4]
        header_format2 = b' 3i'

        #============
        header_format = b''  # the main block
        header = []
        for word in self.words:
            if word == '???':
                val = 0
            elif word in self.dataNames:
                val = getattr(self, word + 's')[i]  # self.times[i]
            else:
                val = getattr(self, word)  # self.time

            if isinstance(val, int):
                header_format += b'i'
            elif isinstance(val, float):
                header_format += b'f'
            elif isinstance(val, string_types) and word in ['Title', 'subtitle', 'label']:
                val = 4
                header_format += 'i'
                #val = '%128s' % val
                #header_format += b'128s'
            else:  # I don't think strings are allowed, but if they are, it's 4s
                raise RuntimeError('format(%r) = %r ???' % (word, val))
            header.append(val)

        #============
        Formats = header_format1 + header_format + header_format2
        headers = header1 + header + header2
        #print('Formats =', Formats)
        #print('headers =', headers)
        #f.write(Struct(header_format1).pack(*header1))
        #f.write(Struct(header_format2).pack(*header2))
        #f.write(Struct(header_format).pack(*header))

        f.write(Struct(Formats).pack(*headers))

    def write_op2(self, f, is_mag_phase=False):
        #print('data_code =', self.data_code)
        self._write_table_header(f)
        if isinstance(self.nonlinear_factor, float):
            op2_format = '%sif' % (7 * self.ntimes)
        else:
            op2_format = '2i6f' * self.ntimes
        s = Struct(op2_format)

        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]
        format_table4_1 = Struct(b'9i')
        format_table4_2 = Struct(b'3i')

        # table 4 info
        nnodes = self.data.shape[0]
        nnodes_device = self.node_gridtype[:, 0] * 10 + self.device_code

        #(2+6) => (node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i)
        ntotal = self.ntimes * nnodes * (2 + 6)

        table_num = 3
        for itime in range(self.ntimes):
            self._write_op2_header(f, table_num, itime)

            # record 4
            header = [4, 0, 4,
                      4, -table_num - 1, 4,
                      4, 4 * ntotal, 4]
            f.write(format_table4_1.pack(*header))

            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            i = 0

            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(nnodes_device, gridtype, t1, t2, t3, r1, r2, r3):
                vals = (node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i)
                #grid = nodeID*10+device_code
                f.write(s.pack(*vals))

            table_num -= 2
            header = [4, 4 * ntotal, 4]
            f.write(format_table4_2.pack(*header))
        #return n

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

            if isinstance(dt, float):
                header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            else:
                header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            f.write(''.join(header + words))
            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
                sgridtype = self.recast_gridtype_as_string(gridtypei)
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                (dx, dy, dz, rx, ry, rz) = vals2
                f.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, sgridtype, dx, dy, dz, rx, ry, rz))

            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexTableArray(TableArray):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        TableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def data_type(self):
        return 'complex64'

    #def _write_f06_block(self, words, header, page_stamp, page_num, f, is_mag_phase):
        #self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_mag_phase)

    def _write_f06_transient_block(self, words, header, page_stamp, page_num, f, is_mag_phase):
        if is_mag_phase:
            words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        else:
            words += ['                                                          (REAL/IMAGINARY)\n', ]

        words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()

        if not len(header) >= 3:
            header.append('')
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
                #if not is_all_zeros:
                f.write('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                        '  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                            node_id, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                        '', '', dxi, dyi, dzi, rxi, ryi, rzi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealTableObject(ScalarObject):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.nonlinear_factor = None
        self.table_name = None
        self.analysis_code = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.gridTypes = {}
        self.translations = {}
        self.rotations = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

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
            (nodeID, grid_type, t1, t2, t3, r1, r2, r3) = line
            self.gridTypes[nodeID] = grid_type
            self.translations[dt][nodeID] = array([t1, t2, t3], dtype='float32')
            self.rotations[dt][nodeID] = array([r1, r2, r3], dtype='float32')

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

    def add(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        msg = "node_id=%s gridType=%s v1=%s v2=%s v3=%s" % (
            node_id, grid_type, v1, v2, v3)
        #print(msg)
        assert -1 < node_id < 1000000000, msg
        assert isinstance(node_id, int), msg
        #assert nodeID not in self.translations,'displacementObject - static failure'

        self.gridTypes[node_id] = self.recast_gridtype_as_string(grid_type)
        self.translations[node_id] = array([v1, v2, v3], dtype='float32')  # dx,dy,dz
        self.rotations[node_id] = array([v4, v5, v6], dtype='float32')  # rx,ry,rz

    def add_sort1(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        msg = "node_id=%8i v1=%s v2=%s v3=%s\n" % (node_id, v1, v2, v3)
        msg += "                 v4=%s v5=%s v6=%s" % (v4, v5, v6)
        #print(msg)
        if dt not in self.translations:
            self.add_new_transient(dt)
        assert -1 < node_id < 1000000000, msg
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
                nodeIDs = translations.keys()
                for nodeID in nodeIDs:
                    translations2[nodeID] = {}
                    rotations2[nodeID] = {}

            for dt, translations in sorted(iteritems(self.translations)):
                for nodeID, translation in sorted(iteritems(translations)):
                    rotation = self.rotations[dt][nodeID]
                    translations2[nodeID][dt] = translation
                    rotations2[nodeID][dt] = rotation
        else:
            return (self.translations, self.rotations)
            #for nodeID,translation in sorted(iteritems(self.translations)):
            #    rotation = self.rotations[nodeID]
            #    translations2[nodeID] = translation
            #    rotations2[nodeID]    = rotation
        return (translations2, rotations2)

    def _write_f06_block(self, words, header, page_stamp, page_num=1, f=None):
        msg = words
        #assert f is not None # remove
        for nodeID, translation in sorted(iteritems(self.translations)):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, is_all_zeros) = writeFloats13E(vals)
            #if not is_all_zeros:
            (dx, dy, dz, rx, ry, rz) = vals2
            msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                    % (nodeID, grid_type, dx, dy, dz, rx, ry, rz.rstrip()))
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient_block(self, words, header, page_stamp, page_num=1, f=None):
        msg = []
        #assert f is not None # remove
        for dt, translations in sorted(iteritems(self.translations)):
            if isinstance(dt, float):  # fix
                #header[1] = ' %s = %10.4E float %s\n' % (self.data_code['name'], dt, self.analysis_code)
                header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            else:
                #header[1] = ' %s = %10i integer %s\n' % (self.data_code['name'], dt, self.analysis_code)
                header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            msg += header + words
            for nodeID, translation in sorted(iteritems(translations)):
                rotation = self.rotations[dt][nodeID]
                grid_type = self.gridTypes[nodeID]

                (dx, dy, dz) = translation
                (rx, ry, rz) = rotation
                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                #if not is_all_zeros:
                [dx, dy, dz, rx, ry, rz] = vals2
                msg.append('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (nodeID, grid_type, dx, dy, dz, rx, ry, rz))

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
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.nonlinear_factor = None
        self.table_name = None
        self.analysis_code = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.gridTypes = {}
        self.translations = {}
        self.rotations = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2

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
                    (nodeID, grid_type, v1, v2, v3, v4, v5, v6) = line[2:]
                except:
                    print('line = %r' % line)
                    raise
                self.gridTypes[nodeID] = grid_type
                self.translations[self.dt][nodeID] = array([v1, v2, v3], dtype='complex64')  # dx,dy,dz
                self.rotations[self.dt][nodeID] = array([v4, v5, v6], dtype='complex64')  # rx,ry,rz
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.translations:
            self.update_dt(self.data_code, dt)

        for line in data:
            try:
                (nodeID, grid_type, v1, v2, v3, v4, v5, v6) = line
            except:
                print('line = %r' % line)
                raise
            #print "*dt=%s line=%s" % (self.dt, str(line))
            self.gridTypes[nodeID] = grid_type
            self.translations[dt][nodeID] = array([v1, v2, v3], dtype='complex64')  # dx,dy,dz
            self.rotations[dt][nodeID] = array([v4, v5, v6], dtype='complex64')  # rx,ry,rz

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

    def add(self, dt, nodeID, grid_type, v1, v2, v3, v4, v5, v6):
        #msg = "dt=%s nodeID=%s v1=%s v2=%s v3=%s" %(dt,nodeID,v1,v2,v3)
        #assert isinstance(nodeID,int),nodeID
        msg = "nodeID=%s v1=%s v2=%s v3=%s\n" % (nodeID, v1, v2, v3)
        msg += "          v4=%s v5=%s v6=%s" % (v4, v5, v6)
        #print(msg)
        assert 0 < nodeID < 1000000000, msg  # -1
        assert -0.5 < dt, msg              # remove
        assert isinstance(nodeID, int), msg
        #assert nodeID not in self.translations,'complexDisplacementObject - static failure'

        self.gridTypes[nodeID] = self.recast_gridtype_as_string(grid_type)
        self.translations[nodeID] = array([v1, v2, v3], dtype='complex64')  # dx,dy,dz
        self.rotations[nodeID] = array([v4, v5, v6], dtype='complex64')  # rx,ry,rz

    def add_sort1(self, dt, nodeID, grid_type, v1, v2, v3, v4, v5, v6):
        #msg = "dt=%s nodeID=%s v1=%s v2=%s v3=%s" %(dt,nodeID,v1,v2,v3)
        #print msg
        if dt not in self.translations:
            self.add_new_transient(dt)
        #print msg
        if isinstance(nodeID, float):
            msg = "nodeID=%s v1=%s v2=%s v3=%s\n" % (nodeID, v1, v2, v3)
            msg += "          v4=%s v5=%s v6=%s" % (v4, v5, v6)
            print(msg)

        if not 0 < nodeID < 1000000000:
            msg = "nodeID=%s v1=%s v2=%s v3=%s\n" % (nodeID, v1, v2, v3)
            msg += "          v4=%s v5=%s v6=%s" % (v4, v5, v6)
            raise ValueError(msg)
        #assert -0.5 < dt, msg  # TODO: remove
        #assert isinstance(nodeID, int), msg  # TODO: remove
        #assert nodeID not in self.translations,'complexDisplacementObject - static failure'

        self.gridTypes[nodeID] = self.recast_gridtype_as_string(grid_type)
        self.translations[dt][nodeID] = array([v1, v2, v3], dtype='complex64')  # dx,dy,dz
        self.rotations[dt][nodeID] = array([v4, v5, v6], dtype='complex64')  # rx,ry,rz

    def add_sort2(self, nodeID, data):
        [dt, grid_type, v1, v2, v3, v4, v5, v6] = data

        if dt not in self.translations:
            self.add_new_transient(dt)

        msg = "dt=%s nodeID=%s v1=%s v2=%s v3=%s\n" % (dt, nodeID, v1, v2, v3)
        msg += "                v4=%s v5=%s v6=%s" % (v4, v5, v6)
        #print msg
        assert 0 < nodeID < 1000000000, msg  # -1
        assert -0.5 < dt, msg  # remove
        assert isinstance(nodeID, int), msg

        self.gridTypes[nodeID] = self.recast_gridtype_as_string(grid_type)
        self.translations[dt][nodeID] = array([v1, v2, v3], dtype='complex64')  # dx,dy,dz
        self.rotations[dt][nodeID] = array([v4, v5, v6], dtype='complex64')  # rx,ry,rz

    def _write_f06_block(self, words, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
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
            f.write('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                    '  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                      % (node_id, grid_type, dxr, dyr, dzr, rxr, ryr, rzr,
                                    '', '', dxi, dyi, dzi, rxi, ryi, rzi))

        f.write(page_stamp % page_num)
        msg = ['']
        return page_num

    def _write_f06_transient_block(self, words, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        if is_mag_phase:
            words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        else:
            words += ['                                                          (REAL/IMAGINARY)\n', ]

        words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()

        msg = []
        if not len(header) >= 3:
            header.append('')
        #assert f is not None
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
                msg.append('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (node_id, grid_type, dxr, dyr, dzr, rxr, ryr, rzr))
                msg.append('  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % ('', '', dxi, dyi, dzi, rxi, ryi, rzi))

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1
