# pylint: disable=E1101
from six import iteritems
from pyNastran.op2.resultObjects.tableObject import RealTableArray  #, ComplexTableArray, TableObject, ComplexTableObject
from pyNastran.op2.resultObjects.op2_Objects import ScalarObject


class RealTemperatureArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, page_num, f)
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        return self._write_f06_block(words, header, pageStamp, page_num, f, write_words=False)

    def _write_f06_transient(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        return self._write_f06_transient_block(words, header, pageStamp, page_num, f, write_words=False)


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
            self.add = self.addSort2

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

    def add_new_transient(self, dt):
        """initializes the transient variables"""
        self.temperatures[dt] = {}

    def add(self, dt, nodeID, grid_type, v1, v2, v3, v4, v5, v6):
        # v2-v6 are 0
        assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)
        #assert nodeID not in self.temperatures

        grid_type = self.recast_gridtype_as_string(grid_type)
        self.gridTypes[nodeID] = grid_type
        self.temperatures[nodeID] = v1

    def add_sort1(self, dt, nodeID, grid_type, v1, v2, v3, v4, v5, v6):
        # v2-v6 are 0
        if dt not in self.temperatures:
            self.add_new_transient(dt)

        assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)
        #assert nodeID not in self.temperatures[self.dt]

        grid_type = self.recast_gridtype_as_string(grid_type)
        self.gridTypes[nodeID] = grid_type
        self.temperatures[dt][nodeID] = v1

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

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
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
                msg.append(pageStamp % page_num)
                f.write(''.join(msg))
                msg = ['']
                page_num += 1
            return page_num - 1  # transient

        msg += self.print_temp_lines(self.temperatures)
        msg.append(pageStamp % page_num)
        f.write(''.join(msg))
        return page_num  # static

    def print_temp_lines(self, temperatures):
        msg = []
        ipack = []
        oldNodeID = -1
        oldGridType = None
        for nodeID, T in sorted(iteritems(temperatures)):
            grid_type = self.gridTypes[nodeID]

            if oldNodeID + 1 == nodeID and grid_type == oldGridType:
                oldNodeID = nodeID
                ipack.append(T)
            else:
                if oldNodeID > 0:
                    msg += self.print_pack(ipack)
                oldGridType = grid_type
                oldNodeID = nodeID
                ipack = [nodeID, grid_type, T]
        if ipack:
            msg += self.print_pack(ipack)
        return msg

    def print_pack(self, ipack):
        msg = []
        nID = ipack[0]
        gType = ipack[1]
        while len(ipack) > 8:
            nID = ipack[0]
            packOut = ipack[:8]
            ipack = [nID + 6, gType] + ipack[8:]
            msg.append('      %8i   %4s      %10.6E   %10.6E   %10.6E   %10.6E   %10.6E   %10.6E\n' % (tuple(packOut)))

        if ipack:
            fmt = '      %8i   %4s   ' + '   %10.6E' * (len(ipack) - 2) + '\n'
            out = fmt % (tuple(ipack))
            msg.append(out)
        return msg
