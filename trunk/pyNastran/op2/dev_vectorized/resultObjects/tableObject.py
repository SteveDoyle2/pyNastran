#from struct import pack
import numpy

import pandas as pd
from numpy import array, zeros

from pyNastran import update_index
from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.f06.f06_formatting import writeFloats13E, writeImagFloats13E

class BaseTable(object):
    def __init__(self):
        self.shape = {}
        self._inode_start = None
        self._inode_end = None
        self.data = None

    def _preallocate(self, dt, nnodes):
        if self.shape is None:
            self._inode_start += nnodes
            self._inode_end += nnodes
        elif self.data is not None:
            #print "data =", self.data
            return self._inode_start, self._inode_end
        else:
            ndt, nnodes, dts = self._get_shape()
            #print "ndt=%s nnodes=%s dts=%s" % (ndt, nnodes, str(dts))

            data = {}
            columns = []
            indexs = []
            size_end = ndt * nnodes
            if dts[0] is not None:
                name = self.data_code['name']
                if isinstance(dt, int):
                    data[name] = pd.Series(zeros((size_end), dtype='int32'))
                else:
                    data[name] = pd.Series(zeros((size_end), dtype='float32'))
                columns.append(name)
                indexs = [name]
            data['node_id'] = pd.Series(zeros((size_end), dtype='int32'))
            indexs.append('node_id')

            #data['grid_type'] = pd.Series(zeros(ndt), dtype='int32'))
            #data['grid_type_str'] = pd.Series(zeros(nnodes), dtype='str'))
            float_complex = self._get_float_complex()

            #['T1', 'T2', 'T3', 'R1', 'R2', 'R3']
            headers = self._get_headers()
            for header in headers:
                data[header] = pd.Series(zeros((size_end), dtype=float_complex))

            self._size_start = 0
            self._size_end = size_end
            columns += ['node_id'] + headers

            self.data = pd.DataFrame(data, columns=columns)
            self._inode_start = 0
            self._inode_end = nnodes
        return self._inode_start, self._inode_end

    def _is_full(self, nnodes):
        self._size_start += nnodes
        if self._size_start == self._size_end:
            return True
        return False

    def _get_headers(self):
        return ['T1', 'T2', 'T3', 'R1', 'R2', 'R3']

    def _finalize(self):
        ndt, nnodes, dts = self._get_shape()

        mapper = {
            1 : 'G',
        }
        #grid_type_str = []
        #for grid_type in self.grid_type:
            #grid_type_str.append(mapper[grid_type])
        #self.grid_type_str = pd.Series(grid_type_str, dtype='str')

        if update_index:
            if dts[0] is not None:
                name = self.data_code['name']
                self.data = self.data.set_index([name, 'node_id'])
            else:
                self.data = self.data.set_index('node_id')

        #print "final\n", self.data
        del self._inode_start
        del self._inode_end

    def _increase_size(self, dt, nnodes):
        #self.shape += 1
        if dt in self.shape:  # default dictionary
            self.shape[dt] += nnodes
        else:
            self.shape[dt] = nnodes
            #self.n_nonlinear += 1

    def _get_shape(self):
        ndt = len(self.shape)
        dts = self.shape.keys()
        shape0 = dts[0]
        nnodes = self.shape[shape0]
        #print "ndt=%s nnodes=%s dts=%s" % (ndt, nnodes, str(dts))
        return ndt, nnodes, dts

    def get_stats(self):
        ndt, nnodes, dts = self._get_shape()
        msg = self._get_data_code()

        real_imag = 'real' if self.is_real() else 'imaginary'
        if dts[0] is not None:
            name = self.data_code['name']
            dtstring = name + ', '
            msg.append('  %s type=%s n%ss=%s nnodes=%s\n'
                       % (real_imag, self.__class__.__name__, name, ndt, nnodes))
        else:
            dtstring = ''
            msg.append('  %s type=%s nnodes=%s\n' % (real_imag, self.__class__.__name__, nnodes))
        headers = self._get_headers()
        msg.append('  element_data: index  : node_id\n')
        msg.append('              : results: node_type\n')
        msg.append('  data        : index  : %snode_id\n' % dtstring)
        msg.append('              : results: %s\n' % ', '.join(headers))
        return msg


class TableObject(BaseTable, ScalarObject):  # displacement style table
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        BaseTable.__init__(self)
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        self.nonlinear_factor = None
        self.table_name = None
        self.analysis_code = None
        self.nnonlinear = 0

        # new method, not finished
        self.nodeIDs_to_index = None
        self.gridTypes2 = {}
        self.dt = dt

    def _get_float_complex(self):
        return 'float32'

    def is_imaginary(self):
        return False

    def is_real(self):
        return True

    def add_array_f06_data(self, data, transient):
        nnodes = len(data)
        read_mode = 1
        nodeIDs_to_index = zeros(nnodes, dtype='int32')
        gridTypes = zeros(nnodes, dtype='string')
        translations = zeros((nnodes, 6), dtype='float32')

        for i, line in enumerate(data):
            (nodeID, grid_type, t1, t2, t3, r1, r2, r3) = line
            nodeIDs_to_index[i] = nodeID
            gridTypes[i] = grid_type
            translations[i, :] = array([t1, t2, t3, r1, r2, r3])

        if transient is None:
            self.add_array(None, nodeIDs_to_index, gridTypes, translations)
            #print('%s gridTypes=%s'  % (self.__class__.__name__, gridTypes))
            return
        else:
            (dtName, dt) = transient
            self.data_code['name'] = dtName
            if dt not in self.translations:
                self.update_dt(self.data_code, dt, read_mode)
            self.add_array_sort1(dt, nodeIDs_to_index, gridTypes, translations)

    def add_f06_data(self, data, transient):
        if as_array:
            return self.add_array_f06_data(data, transient)

        #raise RuntimeError('this should be commented out')
        if transient is None:
            for line in data:
                (nodeID, grid_type, t1, t2, t3, r1, r2, r3) = line
                self.gridTypes[nodeID] = grid_type
                self.translations[nodeID] = array([t1, t2, t3])
                self.rotations[nodeID] = array([r1, r2, r3])
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.translations:
            self.update_dt(self.data_code, dt, read_mode)

        for line in data:
            (nodeID, grid_type, t1, t2, t3, r1, r2, r3) = line
            self.gridTypes[nodeID] = grid_type
            self.translations[dt][nodeID] = array([t1, t2, t3])
            self.rotations[dt][nodeID] = array([r1, r2, r3])

    def _write_f06_block(self, words, header, pageStamp, f, pageNum=1):
        msg = words
        ndt, nnodes, dts = self._get_shape()

        #print self.data.to_string()
        headers = self._get_headers()
        assert headers == ['T1', 'T2', 'T3', 'R1', 'R2', 'R3'], headers

        ndata = len(self.data)
        for i in xrange(ndata):
            index = self.data.index[i]
            node_id = index
            grid_type = 'G' #self.gridTypes[nodeID]
            dx = self.data['T1'][index]
            dy = self.data['T2'][index]
            dz = self.data['T3'][index]
            rx = self.data['R1'][index]
            ry = self.data['R2'][index]
            rz = self.data['R3'][index]

            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, isAllZeros) = writeFloats13E(vals)
            if not isAllZeros:
                [dx, dy, dz, rx, ry, rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n'
                        % (node_id, grid_type, dx, dy, dz, rx, ry, rz.rstrip()))
            if i % 1000 == 0:
                f.write(''.join(msg))
                msg = ['']

        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient_block(self, words, header, pageStamp, f,
                                   pageNum=1):
        msg = []
        #assert f is not None # remove
        ndata = len(self.data)
        classname = self.__class__.__name__
        headers = self._get_headers()
        assert headers == ['T1', 'T2', 'T3', 'R1', 'R2', 'R3'], headers

        i = 0
        while i < ndata:
            index = self.data.index[i]
            (dt, node_id) = index
            dt_old = dt
            if isinstance(dt, float):  # @todo: fix
                #header[1] = ' %s = %10.4E float %s\n' % (self.data_code['name'], dt, self.analysis_code)
                header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            else:
                #header[1] = ' %s = %10i integer %s\n' % (self.data_code['name'], dt, self.analysis_code)
                header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            msg += header + words

            while dt == dt_old:
                index = self.data.index[i]
                (dt, node_id) = index

                data = self.data.ix[index]
                grid_type = 'G' #self.gridTypes[node_id]
                dx = data['T1']
                dy = data['T2']
                dz = data['T3']
                rx = data['R1']
                ry = data['R2']
                rz = data['R3']

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, isAllZeros) = writeFloats13E(vals)
                if not isAllZeros:
                    [dx, dy, dz, rx, ry, rz] = vals2
                    msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (node_id, grid_type, dx, dy, dz, rx, ry, rz.rstrip()))
                i += 1
                try:
                    dt = self.data.index[i+1][0]
                except IndexError:
                    break

                if i % 1000 == 0:
                    #print("class=%s i=%s ndata=%s" % (classname, i, ndata))
                    f.write(''.join(msg))
                    msg = ['']
            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        #print('returning...')
        return pageNum - 1

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
        if 'ATO' in self.table_name:
            return True
        return False

    def isCRM(self):
        """Correlated Root-Mean Square"""
        if 'CRM' in self.table_name:
            return True
        return False

    def isPSD(self):
        """Power Spectral Density"""
        if 'PSD' in self.table_name:
            return True
        return False

    def isRMS(self):
        """Root-Mean Square"""
        if 'RMS' in self.table_name:
            return True
        return False

    def isZERO(self):
        """Zero Crossings"""
        if 'NO' in self.table_name:
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


class ComplexTableObject(BaseTable, ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        BaseTable.__init__(self)
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        self.nonlinear_factor = None
        self.table_name = None
        self.analysis_code = None
        self.gridTypes = {}
        self.translations = {}
        self.rotations = {}
        #self._inode_start = None
        #self._inode_end = None
        self.dt = dt

    def _get_float_complex(self):
        return 'complex64'

    def is_imaginary(self):
        return True

    def is_real(self):
        return False

    def add_f06_data(self, data, transient):
        read_mode = 1
        if transient is None:
            for line in data:
                try:
                    (nodeID, grid_type, v1, v2, v3, v4, v5, v6) = line
                except:
                    print('line = %r' % line)
                    raise
                self.gridTypes[nodeID] = grid_type
                self.translations[self.dt][nodeID] = [v1, v2, v3]  # dx,dy,dz
                self.rotations[self.dt][nodeID] = [v4, v5, v6]  # rx,ry,rz
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.translations:
            self.update_dt(self.data_code, dt, read_mode)

        for line in data:
            try:
                (nodeID, grid_type, v1, v2, v3, v4, v5, v6) = line
            except:
                print('line = %r' % line)
                raise
            self.gridTypes[nodeID] = grid_type
            self.translations[self.dt][nodeID] = [v1, v2, v3]  # dx,dy,dz
            self.rotations[self.dt][nodeID] = [v4, v5, v6]  # rx,ry,rz

    def update_dt(self, data_code, dt, read_mode):
        self.data_code = data_code
        self.apply_data_code()

    def _write_f06_block(self, words, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        #words += self.getTableMarker()
        if is_mag_phase:
            words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        else:
            words += ['                                                          (REAL/IMAGINARY)\n', ]

        words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']

        msg = words
        for nodeID, translation in sorted(self.translations.iteritems()):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]
            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation

            vals = [dx, dy, dz, rx, ry, rz]
            (vals2, isAllZeros) = writeImagFloats13E(vals)
            [dxr, dyr, dzr, rxr, ryr, rzr, dxi, dyi, dzi, rxi,
                ryi, rzi] = vals2
            msg.append('0 %12i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (nodeID, grid_type, dxr, dyr, dzr, rxr, ryr, rzr.rstrip()))
            msg.append('  %12s %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %
                       ('', '', dxi, dyi, dzi, rxi, ryi, rzi.rstrip()))

        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient_block(self, words, header, pageStamp, f,
                                   pageNum=1, is_mag_phase=False):
        if is_mag_phase:
            words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        else:
            words += ['                                                          (REAL/IMAGINARY)\n', ]

        words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']

        msg = []
        #assert f is not None # remove
        ndata = len(self.data)
        classname = self.__class__.__name__
        i = 0
        while i < ndata:
            index = self.data.index[i]
            (dt, node_id) = index
            dt_old = dt
            if isinstance(dt, float):  # @todo: fix
                #header[1] = ' %s = %10.4E float %s\n' % (self.data_code['name'], dt, self.analysis_code)
                header[2] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            else:
                #header[1] = ' %s = %10i integer %s\n' % (self.data_code['name'], dt, self.analysis_code)
                header[2] = ' %s = %10i\n' % (self.data_code['name'], dt)
            msg += header + words

            while dt == dt_old:
                index = self.data.index[i]
                (dt, node_id) = index

                data = self.data.ix[index]
                grid_type = 'G' #self.gridTypes[node_id]
                dx = data['T1']
                dy = data['T2']
                dz = data['T3']
                rx = data['R1']
                ry = data['R2']
                rz = data['R3']

                vals = [dx, dy, dz, rx, ry, rz]
                (vals2, isAllZeros) = writeImagFloats13E(vals, is_mag_phase)
                [dxr, dyr, dzr, rxr, ryr, rzr,
                 dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                if not isAllZeros:
                    msg.append('0 %12i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % (node_id, grid_type, dxr, dyr, dzr, rxr, ryr, rzr.rstrip()))
                    msg.append('  %12s %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' % ('', '',            dxi, dyi, dzi, rxi, ryi, rzi.rstrip()))

                i += 1
                try:
                    dt = self.data.index[i+1][0]
                except IndexError:
                    break

                if i % 1000 == 0:
                    #print("class=%s i=%s ndata=%s" % (classname, i, ndata))
                    f.write(''.join(msg))
                    msg = ['']
            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        #print('returning...')
        return pageNum - 1