from six import b
from six.moves import range
from struct import Struct
from pyNastran.op2.op2_common import OP2Common

class OGS(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_ogs1_3(self, data, ndata):
        three = self.parse_approach_code(data)
        self.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', 'dLoadID',
            'format_code', 'num_wide', 'o_code', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        self.parse_approach_code(data)
        #isubcase = self.get_values(data, 'i', 4)

        ## surface/volumeID
        self.ID = self.add_data_parameter(data, 'ID', 'i', 3, False)

        #: Reference coordinate system ID
        self.refid = self.add_data_parameter(data, 'refid', 'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        self.num_wide = self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## Stress/Strain code
        self.sCode = self.add_data_parameter(data, 'sCode', 'i', 11, False)

        ## Output Coordinate System
        self.oCoord = self.add_data_parameter(data, 'oCoord', 'i', 12, False)

        ## Axis Specification code
        self.axis = self.add_data_parameter(data, 'axis', 'i', 13, False)

        #: Normal Specification Code
        self.normal = self.add_data_parameter(data, 'normal', 'i', 14, False)

        self.fix_format_code()
        if not self.is_sort1():
            raise NotImplementedError('OGS sort2...')

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.mode_cycle = 0.0
            self.update_mode_cycle('mode_cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'mode_cycle'])
        #elif self.analysis_code == 3: # differential stiffness
        #elif self.analysis_code == 4: # differential stiffness
        #elif self.analysis_code == 5: # frequency
        elif self.analysis_code == 6:  # transient
            ## time step
            self.time = self.add_data_parameter(data, 'time', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['time'])
        #elif self.analysis_code == 7:  # pre-buckling
        #elif self.analysis_code == 8:  # post-buckling
        #elif self.analysis_code == 9:  # complex eigenvalues
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.lftsfq = self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.data_names = self.apply_data_code_value('data_names', ['lftsfq'])
        #elif self.analysis_code == 11:  # old geometric nonlinear statics
        #elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' % self.analysis_code)

        #print "*isubcase=%s" % (self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)

        #print self.code_information()
        if self.is_debug_file:
            self.binary_debug.write('  approach_code  = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode          = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase       = %r\n' % self.isubcase)
        self._read_title(data)
        self._write_debug_bits()

    def _read_ogs1_4(self, data, ndata):
        if self.read_mode == 1:
            return ndata

        if self.table_code == 26:
            # OGS1 - grid point stresses - surface
            assert self.table_name in [b'OGS1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_ogs1_table26(data, ndata)
        elif self.table_code == 27:
            # OGS1 - grid point stresses - volume direct
            assert self.table_name in [b'OGS1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_ogs1_table27(data, ndata)
        elif self.table_code == 28:
            # OGS1- grid point stresses - principal
            assert self.table_name in [b'OGS1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_ogs1_table28(data, ndata)
        elif self.table_code == 35:
            # OGS - Grid point stress discontinuities (plane strain)
            assert self.table_name in [b'OGS1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_ogs1_table35(data, ndata)
        else:
            raise NotImplementedError('table_code=%s table_name=%s' % (self.table_code, self.table_name))
        return n

    def _read_ogs1_table28(self, data, ndata):
        if self.num_wide == 15:
            pass
        else:
            raise NotImplementedError(self.num_wide)
        return ndata

    def _read_ogs1_table26(self, data, ndata):
        result_name = 'grid_point_stresses'
        if self.num_wide == 11:  # real/random
            #self.create_transient_object(self.gridPointStresses, GridPointStressesObject)
            n = self._read_ogs1_table26_numwide11(data, ndata)
        else:
            msg = 'only num_wide=11 is allowed  num_wide=%s' % self.num_wide
            raise RuntimeError(msg)
        return n

    def _read_ogs1_table26_numwide11(self, data, ndata):
        """surface stresses"""
        #dt = self.nonlinear_factor
        s = Struct(b'2i4s8f')

        n = 0
        nelements = ndata // 44  # 11*4
        for i in range(nelements):
            edata = data[n:n+44]
            out = s.unpack(edata)
            (ekey, eid, fiber, nx, ny, txy, angle, major, minor, tmax, ovm) = out
            nid = ekey // 10
            #fiber = fiber.decode('utf-8').strip()
            assert nid > 0, nid
            #self.obj.add(dt, nid, eid, fiber, nx, ny, txy,
            #             angle, major, minor, tmax, ovm)
            n += 44
        return n

    def _read_ogs1_table27(self, data, ndata):
        """OGS1 - grid point stresses - volume direct"""
        #is_sort1 = self.is_sort1()
        if self.num_wide == 9:  # real/random
            result_name = 'grid_point_volume_stresses'
            #self.create_transient_object(self.gridPointVolumeStresses, GridPointStressesVolumeObject)
            n = self._read_ogs1_table27_numwide9(data, ndata)
        else:
            msg = 'only num_wide=9 is allowed  num_wide=%s' % self.num_wide
            raise RuntimeError(msg)
        return n

    def _read_ogs1_table27_numwide9(self, data, ndata):
        """surface stresses"""
        s = Struct(b(self._endian + '2i7f'))
        n = 0
        nelements = ndata // 36  # 9*4
        for i in range(nelements):
            edata = data[n:n+36]
            out = s.unpack(edata)
            (ekey, nx, ny, nz, txy, tyz, txz, pressure, ovm) = out
            nid = ekey // 10
            assert nid > 0, nid
            #check_nid
            #self.obj.add(dt, nid, nx, ny, nz, txy, tyz, txz, pressure, ovm)
            n += 36
        return n

    def _read_ogs1_table35(self, data, ndata):
        """grid point stress discontinuities (plane stress/strain)"""
        result_name = 'grid_point_stresses'
        if self.num_wide == 6:
            #self.create_transient_object(self.gridPointStresses, GridPointStressesObject)
            n = self._read_ogs1_table35_numwide6(data, ndata)
        else:
            msg = 'only num_wide=11 is allowed  num_wide=%s' % self.num_wide
            raise RuntimeError(msg)
        return n

    def _read_ogs1_table35_numwide6(self, data, ndata):
        """grid point stress discontinuities (plane stress/strain)"""
        s = Struct(b(self._endian + 'i5f'))
        n = 0
        nelements = ndata // 24  # 6*4
        for i in range(nelements):
            out = s.unpack(data[n:n+24])
            (ekey, nx, ny, nz, txy, pressure) = out
            nid = ekey // 10
            assert nid > 0, nid
            #check_nid
            #self.obj.add(dt, nid, nx, ny, nz, txy, tyz, txz, pressure, ovm)
            n += 24
        return n
