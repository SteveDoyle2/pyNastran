from struct import unpack

from pyNastran.op2.op2_helper import polar_to_real_imag


class OQG(object):
    def __init__(self):
        pass


    def read_oqg1_3(self, data):
        three = self.parse_approach_code(data)
        self.words = [
            'aCode',       'tCode',    '???',           'isubcase',
             '???',         '???',      '???',          'random_code'
             'format_code', 'num_wide', '11',           '???',
            'acoustic_flag','???',      '???',          '???',
             '???',         '???',      '???',          '???',
             '???',         '???',      'thermal',      '???',
             '???', 'Title', 'subtitle', 'label']

        ## random code
        self.add_data_parameter( data, 'random_code', 'i', 8, False)

        ## format code
        self.add_data_parameter( data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## acoustic pressure flag
        self.add_data_parameter(data, 'acoustic_flag', 'f', 13, False)

        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.add_data_parameter(data, 'thermal', 'i', 23, False)

        if not self.is_sort1():
            raise NotImplementedError('sort2...')
        #assert self.isThermal()==False,self.thermal

        #self.print_block(data) # on
        ## assuming tCode=1
        if self.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode_cycle'])
        #elif self.analysis_code==3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code==4: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.add_data_parameter(data, 'freq', 'f', 5)
            self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.add_data_parameter(data, 'dt', 'f', 5)
            self.apply_data_code_value('dataNames', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue
            self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.apply_data_code_value('dataNames', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % (self.analysis_code)
            raise RuntimeError(msg)

        #print self.code_information()
        if self.debug:
            self.binary_debug.write('  aCode    = %r\n' % self.aCode)
            self.binary_debug.write('  tCode    = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase = %r\n' % self.isubcase)
        self.read_title(data)
        self.write_debug_bits()

    def write_debug_bits(self):
        if self.debug:
            msg = ''
            for i, param in enumerate(self.words):
                if param == '???':
                    param = 0
                msg += '%s, ' % param
                if i % 5 == 4:
                    msg += '\n             '
            if hasattr(self, 'format_code'):
                self.binary_debug.write('  sort_bits[0] = %i -> is_sort1 =%s\n' % (self.sort_bits[0], self.is_sort1() ))
                self.binary_debug.write('  sort_bits[1] = %i -> is_real  =%s vs real/imag\n' % (self.sort_bits[1], self.is_real()   ))
                self.binary_debug.write('  sort_bits[2] = %i -> is_random=%s vs mag/phase\n' % (self.sort_bits[2], self.is_random() ))
                if self.is_complex():
                    self.binary_debug.write('  format_code  = %i -> is_mag_phase=%s vs is_real_imag\n' % (self.format_code, self.is_mag_phase() ))
                else:
                    self.binary_debug.write('  format_code  = %i\n' % self.format_code)
            self.binary_debug.write('  recordi = [%s]\n\n' % msg)

    def read_oqg1_4(self, data):
        if self.table_code == 3:   # SPC Forces
            assert self.table_name in ['OQG1', 'OQGV1', 'OQP1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
        elif self.table_code == 39:  # MPC Forces
            assert self.table_name in ['OQMG1'], 'table_name=%s table_code=%s' % (
                self.table_name, self.table_code)
            self.read_mpc_forces(data)
        else:
            self.not_implemented_or_skip('bad OQG table')

    def read_spc_forces(self, data):
        result_name = 'SPC_forces'
        self.read_table(data, result_name, 'node')

    def read_mpc_forces(self, data):
        result_name = 'MPC_forces'
        self.read_table(data, result_name, 'node')

    def read_table(self, data, result_name, flag):
        if self.num_wide == 8:  # real/random
            if self.thermal == 0:
                #obj = self.create_transient_object(self.spcForces, SPCForcesObject)
                self.read_real_table(data, result_name, flag)
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()

    def read_real_table(self, data, result_name, flag):
        if self.debug:
            self.binary_debug.write('  read_real_table\n')
        assert flag in ['node', 'elem'], flag
        format1 = '2i6f' # 8

        ntotal = 32 # 8 * 4
        nnodes = len(data) // ntotal > 0
        assert len(data) % ntotal == 0

        n = 0
        for inode in xrange(nnodes):
            eData = data[n:n+32]
            out = unpack(format1, eData)
            (eid_device, gridType, tx, ty, tz, rx, ry, rz) = out
            eid = (eid_device - self.device_code) // 10
            #print "eType=%s" %(eType)

            dataIn = [eid, gridType, tx, ty, tz, rx, ry, rz]
            #self.obj.add(dt, dataIn)
            n += ntotal

    def read_complex_table(self, data, result_name, flag):
        if self.debug:
            self.binary_debug.write('  read_real_table\n')
        assert flag in ['node', 'elem'], flag

        format1 = '2i12f'
        is_magnitude_phase = self.is_magnitude_phase()

        n = 0
        ntotal = 56  # 14 * 4
        nnodes = len(data) // ntotal
        for inode in xrange(nnodes):
            edata = data[n:n+ntotal]

            out = unpack(format1, edata)
            if self.debug:
                self.binary_debug.write('read_complex_table - %s\n' % str(out))
            (eid_device, gridType, txr, tyr, tzr, rxr, ryr, rzr,
                                   txi, tyi, tzi, rxi, ryi, rzi) = out

            if is_magnitude_phase:
                tx = polar_to_real_imag(txr, txi)
                rx = polar_to_real_imag(rxr, rxi)
                ty = polar_to_real_imag(tyr, tyi)
                ry = polar_to_real_imag(ryr, ryi)
                tz = polar_to_real_imag(tzr, tzi)
                rz = polar_to_real_imag(rzr, rzi)
            else:
                tx = complex(txr, txi)
                rx = complex(rxr, rxi)
                ty = complex(tyr, tyi)
                ry = complex(ryr, ryi)
                tz = complex(tzr, tzi)
                rz = complex(rzr, rzi)

            eid = (eid_device - self.device_code) // 10

            dataIn = [eid, gridType, tx, ty, tz, rx, ry, rz]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            #self.obj.add(dt, dataIn)
            n += ntotal