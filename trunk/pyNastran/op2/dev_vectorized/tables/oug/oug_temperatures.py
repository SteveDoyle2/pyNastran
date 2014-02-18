# pylint: disable=E1101
#import sys
#from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.op2.resultObjects.tableObject import RealTableObject, ComplexTableObject


class TemperatureObject(RealTableObject):  # approach_code=1, sort_code=0, thermal=1
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt, read_mode)
        #ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.gridTypes = {}
        #self.temperatures = {}
        #self.dt = dt

    def _get_stats(self):
        ngrids = len(self.gridTypes)

        msg = self._get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.temperatures)
            msg.append('  type=%s ntimes=%s ngrids=%s\n'
                       % (self.__class__.__name__, ntimes, ngrids))
        else:
            msg.append('  type=%s ngrids=%s\n' % (self.__class__.__name__,
                                                  ngrids))
        msg.append('  translations, rotations, gridTypes\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (gridID, grid_type) = line[0:2]
                temps = line[2:]
                for (i, temp) in enumerate(temps):
                    nodeID = gridID + i
                    self.gridTypes[nodeID] = grid_type
                    self.temperatures[nodeID] = temp
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.temperatures:
            self.update_dt(self.data_code, dt)
            self.isTransient = True

        for line in data:
            (gridID, grid_type) = line[0:2]
            temps = line[2:]
            for (i, temp) in enumerate(temps):
                nodeID = gridID + i
                self.gridTypes[nodeID] = grid_type
                self.temperatures[dt][nodeID] = temp

    def update_dt(self, data_code, dt):
        self.data_code = data_code
        self.apply_data_code()
        if dt is not None:
            self.log.debug("updating %s...%s=%s  isubcase=%s" % (self.data_code['name'], self.data_code['name'], dt, self.isubcase))
            self.dt = dt
            self.add_new_transient(dt)

    def delete_transient(self, dt):
        del self.temperatures[dt]

    def get_transients(self):
        k = self.temperatures.keys()
        k.sort()
        return k

   # def write_op2(self,block3,device_code=1):
   #     """
   #     creates the binary data for writing the table
   #     .. warning:: hasnt been tested...
   #     """
   #     msg = block3
   #     for nodeID,T in sorted(self.temperatures.iteritems()):
   #         grid = nodeID*10+device_code
   #         msg += pack('iffffff',grid,T,0,0,0,0,0)
   #     return msg
   #
   # def write_op2_transient(self,block3,device_code=1):
   #     """
   #     creates the binary data for writing the table
   #     .. warning:: hasnt been tested...
   #     .. warning:: dt slot needs to be fixed...
   #     """
   #     msg = ''
   #     for dt,temperatures in sorted(self.temperatures.iteritems()):
   #         XXX = 50 ## this isnt correct... .. todo:: update dt
   #         msg += block3[0:XXX] + pack('i',dt) + block3[XXX+4:]
   #         #msg += '%s = %g\n' %(self.data_code['name'],dt)
   #
   #         for nodeID,T in sorted(temperatures.iteritems()):
   #             grid = nodeID*10+device_code
   #             msg += pack('iffffff',grid,T,0,0,0,0,0)
   #     return msg

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        words += self.get_table_marker()
        return self._write_f06_block(words, header, pageStamp, f, pageNum)

    def _get_headers(self):
        return ['Temperature']

    def _write_f06_block(self, words, header, pageStamp, f, pageNum):
        ndata = len(self.data)
        ndt, nnodes, dts = self._get_shape()
        msg = ''
        i = 0
        index = self.data.index[i]
        (node_id_start) = index

        pack = []
        for inode in xrange(nnodes):
            node_id = self.data.index[i]

            grid_type = 'S'  # TODO: incorrect grid_type
            T = self.data['Temperature'][index]  # TODO: wrong label
            pack.append(T)

            try:
                next_node_id = self.data.index[i+1]
                next_grid_type = 'S'  # TODO: incorrect grid_type
            except IndexError:
                break # done with writing f06

            # TODO doesn't check grid type
            if len(pack) == 6 or (next_node_id != (node_id + 1)): # 6*2 # 6 spots, 2 values
                #print pack
                msg += write_pack(node_id_start, grid_type, pack)
                pack = []
                node_id_start = next_node_id
            if i % 1000:
                f.write(msg)
                msg = ''
            i += 1

        if len(pack): # 6*2 # 6 spots, 2 values
            msg += write_pack(node_id_start, grid_type, pack)
            pack = []
            f.write(msg)
            msg = ''
        msg += pageStamp + str(pageNum) + '\n'
        f.write(msg)
        return pageNum

    def _write_f06_transient_block(self, words, header, pageStamp, f, pageNum):
        ndata = len(self.data)
        ndt, nnodes, dts = self._get_shape()
        print('dts =', dts)
        msg = ''
        i = 0
        for dt in dts:
            index = self.data.index[i]
            (dt, node_id_start) = index
            if isinstance(dt, float):  # @todo: fix
                header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            else:
                header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            msg += ''.join(header + words)

            pack = []
            for inode in xrange(nnodes):
                index = self.data.index[i]
                (dt, node_id) = index
                grid_type = 'S'  # TODO: incorrect grid_type
                T = self.data['Temperature'][index]  # TODO: wrong label
                pack.append(T)

                try:
                    _next_dt, next_node_id = self.data.index[i+1]
                    next_grid_type = 'S'  # TODO: incorrect grid_type
                except IndexError:
                    break # done with writing f06

                # TODO doesn't check grid type
                if len(pack) == 6 or (next_node_id != (node_id + 1)): # 6*2 # 6 spots, 2 values
                    #print pack
                    msg += write_pack(node_id_start, grid_type, pack)
                    pack = []
                    node_id_start = next_node_id
                if i % 1000:
                    f.write(msg)
                    msg = ''
                i += 1

            if len(pack): # 6*2 # 6 spots, 2 values
                msg += write_pack(node_id_start, grid_type, pack)
                pack = []
                f.write(msg)
                msg = ''
            msg += pageStamp + str(pageNum) + '\n'
            pageNum += 1
        f.write(msg)
        return pageNum -  1

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1):
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, pageStamp, f, pageNum)


    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        msg = []
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(header, pageStamp, f, pageNum)

        return self._write_f06_block(words, header, pageStamp, f, pageNum)
        #msg += self.print_temp_lines(self.temperatures)
        #msg.append(pageStamp + str(pageNum) + '\n')
        #f.write(''.join(msg))
        #return pageNum  # static


def write_pack(node, grid_type, pack):
    msg = '      %8i   %4s    ' % (node, grid_type)
    for p in pack:
        msg += '  %10.6e' % p
    return msg + '\n'
