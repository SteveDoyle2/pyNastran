from struct import Struct, unpack


class OGS(object):
    def __init__(self):
        pass

    def _read_ogs1_3(self, data):
        three = self.parse_approach_code(data)
        self.words = [
            'aCode',       'tCode',     '???',          'isubcase',
             '???',         '???',      '???',          'dLoadID'
             'format_code', 'num_wide', 'o_code',       '???',
            'acoustic_flag','???',      '???',          '???',
             '???',         '???',      '???',          '???',
             '???',         '???',      'thermal',      '???',
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

        #print "dLoadID(8)=%s format_code(9)=%s num_wide(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.format_code,self.num_wide,self.oCode,self.thermal)
        if not self.is_sort1():
            raise NotImplementedError('sort2...')

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics
            ## load set number
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.eign = self.add_data_parameter(data, 'eign', 'f', 6, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['mode', 'eign'])
        #elif self.analysis_code == 3: # differential stiffness
        #elif self.analysis_code == 4: # differential stiffness
        #elif self.analysis_code == 5: # frequency
        elif self.analysis_code == 6:  # transient
            ## time step
            self.time = self.add_data_parameter(data, 'time', 'f', 5)
            self.dataNames = self.apply_data_code_value('dataNames', ['time'])
        #elif self.analysis_code == 7:  # pre-buckling
        #elif self.analysis_code == 8:  # post-buckling
        #elif self.analysis_code == 9:  # complex eigenvalues
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.loadstep = self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.apply_data_code_value('dataNames', ['lftsfq'])
        #elif self.analysis_code == 11:  # old geometric nonlinear statics
        #elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' % self.analysis_code)

        #print "*isubcase=%s" % (self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)

        #print self.code_information()
        if self.debug:
            self.binary_debug.write('  aCode    = %r\n' % self.aCode)
            self.binary_debug.write('  tCode    = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase = %r\n' % self.isubcase)
        self.read_title(data)
        self.write_debug_bits()

    def _read_ogs1_4(self, data):
        if self.table_code == 26:  # OGS1 - grid point stresses - surface
            assert self.table_name in ['OGS1'], 'table_name=%s table_code=%s' % (
                self.table_name, self.table_code)
            self._read_og1s_table26(data)
        elif self.table_code == 27:  # OGS1 - grid point stresses - volume direct
            assert self.table_name in ['OGS1'], 'table_name=%s table_code=%s' % (
                self.table_name, self.table_code)
            self._read_og1s_table27(data)
        #elif self.table_code == 28:  # OGS1- grid point stresses - principal
            #assert self.table_name in ['OGS1'],'table_name=%s table_code=%s' %(self.table_name,self.table_code)
            #self.readOGS1_Data_table28(data)
            #self.not_implemented_or_skip()

        #elif self.table_code == 35:  # OGS - Grid point stress discontinuities (plane strain)
            #self.not_implemented_or_skip()
        else:
            raise NotImplementedError()

    def _read_og1s_table26(self, data):
        resultName = 'gridPointStresses'
        if self.num_wide == 11:  # real/random
            #self.create_transient_object(self.gridPointStresses, GridPointStressesObject)
            self.readOGS1_table26_numWide11(data)
        else:
            msg = 'only num_wide=11 is allowed  num_wide=%s' % (self.num_wide)
            raise RuntimeError(msg)

    def readOGS1_table26_numWide11(self, data):  # surface stresses
        #dt = self.nonlinear_factor
        format1 = b'2i4s8f'
        s = Struct(format1)

        n = 0
        nelements = len(data) // 44  # 11*4
        for i in xrange(nelements):
            edata = data[n:n+44]
            out = s.unpack(edata)
            (ekey, eid, fiber, nx, ny, txy, angle, major,
                minor, tmax, ovm) = out
            nid = (ekey - self.device_code) // 10
            #fiber = fiber.decode('utf-8').strip()
            check_nid
            #self.obj.add(dt, nid, eid, fiber, nx, ny, txy,
            #             angle, major, minor, tmax, ovm)

    def _read_og1s_table27(self, data):  # OGS1 - grid point stresses - volume direct
        #is_sort1 = self.is_sort1()
        print(self.code_information())
        if self.num_wide == 9:  # real/random
            resultName = 'gridPointVolumeStresses'
            #self.create_transient_object(self.gridPointVolumeStresses, GridPointStressesVolumeObject)
            self.readOGS1_table27_numWide9(data)
        else:
            msg = 'only num_wide=9 is allowed  num_wide=%s' % self.num_wide
            raise RuntimeError(msg)

    def readOGS1_table27_numWide9(self, data):  # surface stresses
        format1 = b'2i7f'
        s = Struct(format1)

        n = 0
        nelements = len(data) // 36  # 9*4
        for i in xrange(nelements):
            edata = data[n:n+36]
            out = unpack(format1, edata)
            (ekey, nx, ny, nz, txy, tyz, txz, pressure, ovm) = out
            nid = (ekey - self.device_code) // 10
            check_nid
            #self.obj.add(dt, nid, nx, ny, nz, txy, tyz, txz, pressure, ovm)
