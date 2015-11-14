# pylint: disable=E1101
from six import iteritems
from pyNastran.op2.resultObjects.tableObject import RealTableArray  #, ComplexTableArray, TableObject, ComplexTableObject
from pyNastran.op2.resultObjects.op2_Objects import ScalarObject


class RealTemperatureArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        return self._write_f06_block(words, header, page_stamp, page_num, f, write_words=False)

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f, write_words=False)


class RealTemperature(ScalarObject):  # approach_code=1, sort_code=0, thermal=1
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.gridTypes = {}
        self.temperatures = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        ngrids = len(self.gridTypes)

        msg = self.get_data_code()
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
                    node_id = gridID + i
                    self.gridTypes[node_id] = grid_type
                    self.temperatures[node_id] = temp
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
                node_id = gridID + i
                self.gridTypes[node_id] = grid_type
                self.temperatures[dt][node_id] = temp

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

    def add_new_transient(self, dt):
        """initializes the transient variables"""
        self.temperatures[dt] = {}

    def add(self, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        # v2-v6 are 0
        assert 0 < node_id < 1000000000, 'node_id=%s' % node_id
        #assert nodeID not in self.temperatures

        grid_type = self.recast_gridtype_as_string(grid_type)
        self.gridTypes[node_id] = grid_type
        self.temperatures[node_id] = v1

    def add_sort1(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        # v2-v6 are 0
        if dt not in self.temperatures:
            self.add_new_transient(dt)

        assert 0 < node_id < 1000000000, 'node_id=%s' % node_id
        #assert node_id not in self.temperatures[self.dt]

        grid_type = self.recast_gridtype_as_string(grid_type)
        self.gridTypes[node_id] = grid_type
        self.temperatures[dt][node_id] = v1

    def add_sort2(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        # v2-v6 are 0
        if dt not in self.temperatures:
            self.add_new_transient(dt)

        assert 0 < node_id < 1000000000, 'node_id=%s' % node_id
        #assert node_id not in self.temperatures[self.dt]

        grid_type = self.recast_gridtype_as_string(grid_type)
        self.gridTypes[node_id] = grid_type
        self.temperatures[dt][node_id] = v1

   # def write_op2(self,block3,device_code=1):
   #     """
   #     creates the binary data for writing the table
   #     .. warning:: hasnt been tested...
   #     """
   #     msg = block3
   #     for nodeID,T in sorted(iteritems(self.temperatures)):
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
   #     for dt,temperatures in sorted(iteritems(self.temperatures)):
   #         XXX = 50 ## this isnt correct... .. todo:: update dt
   #         msg += block3[0:XXX] + pack('i',dt) + block3[XXX+4:]
   #         #msg += '%s = %g\n' %(self.data_code['name'],dt)
   #
   #         for nodeID,T in sorted(iteritems(temperatures)):
   #             grid = nodeID*10+device_code
   #             msg += pack('iffffff',grid,T,0,0,0,0,0)
   #     return msg

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        msg = []
        if self.nonlinear_factor is not None:
            for dt, temperatures in sorted(iteritems(self.temperatures)):
                dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
                header[2] = dtLine
                msg += header + words
                msg += self.print_temp_lines(temperatures)
                msg.append(page_stamp % page_num)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
            return page_num - 1  # transient

        msg += self.print_temp_lines(self.temperatures)
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num  # static

    def print_temp_lines(self, temperatures):
        msg = []
        ipack = []
        oldNodeID = -1
        oldGridType = None
        for node_id, T in sorted(iteritems(temperatures)):
            grid_type = self.gridTypes[node_id]

            if oldNodeID + 1 == node_id and grid_type == oldGridType:
                oldNodeID = node_id
                ipack.append(T)
            else:
                if oldNodeID > 0:
                    msg += self.print_pack(ipack)
                oldGridType = grid_type
                oldNodeID = node_id
                ipack = [node_id, grid_type, T]
        if ipack:
            msg += self.print_pack(ipack)
        return msg

    def print_pack(self, ipack):
        msg = []
        nID = ipack[0]
        gType = ipack[1]
        while len(ipack) > 8:
            nID = ipack[0]
            pack_out = ipack[:8]
            ipack = [nID + 6, gType] + ipack[8:]
            msg.append('      %8i   %4s      %10.6E   %10.6E   %10.6E   %10.6E   %10.6E   %10.6E\n' % (tuple(pack_out)))

        if ipack:
            fmt = '      %8i   %4s   ' + '   %10.6E' * (len(ipack) - 2) + '\n'
            out = fmt % (tuple(ipack))
            msg.append(out)
        return msg
