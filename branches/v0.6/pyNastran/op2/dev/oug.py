from struct import unpack

class OUG(object):
    def __init__(self):
        pass


    def read_oug1_3(self, data):
        three = self.parse_approach_code(data)
        self.words = [
            'aCode',       'tCode',    '???',           'isubcase',
             '???',         '???',      '???',          'random_code'
             'format_code', 'num_wide', '???',          '???',
            'acoustic_flag','???',      '???',          '???',
             '???',         '???',      '???',          '???',
             '???',         '???',      'thermal',      '???',
             '???', 'Title', 'subtitle', 'label']

        ## random code
        self.add_data_parameter(data, 'random_code', 'i', 8, False)

        ## format code
        self.add_data_parameter(data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## acoustic pressure flag
        self.add_data_parameter(data, 'acoustic_flag', 'f', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.add_data_parameter(data, 'thermal', 'i', 23, False)

        if self.analysis_code == 1:   # statics / displacement / heat flux
            # load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            # mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            # real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.add_data_parameter(data, 'mode_cycle', 'i', 7, False)  # mode or cycle .. todo:: confused on the type - F1???
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode_cycle'])
        #elif self.analysis_code==3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code==4: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            # frequency
            self.add_data_parameter(data, 'freq', 'f', 5)
            self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            # time step
            self.add_data_parameter(data, 'dt', 'f', 5)
            self.apply_data_code_value('dataNames', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            # load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            # load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            # real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            # mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            # real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            # imaginary eigenvalue
            self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            # load step
            self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.apply_data_code_value('dataNames', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            # load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            # load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
            raise RuntimeError(msg)

        #print self.code_information()
        if self.debug:
            self.binary_debug.write('  aCode    = %r\n' % self.aCode)
            self.binary_debug.write('  tCode    = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase = %r\n' % self.isubcase)
        self.read_title(data)
        self.write_debug_bits()

    def read_oug1_4(self, data):
        if self.table_code == 1:   # Displacements
            assert self.table_name in ['OUG1', 'OUGV1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            self.read_displacement(data)
        elif self.table_code == 7:
            self.read_eigenvector(data)
        elif self.table_code == 10:
            self.read_velocity(data)
        elif self.table_code == 11:
            self.read_acceleration(data)
        else:
            self.not_implemented_or_skip('bad OUG table')

    def read_displacement(self, data):
        result_name = 'displacements'
        #real_obj = DisplacementObject
        #complex_obj = ComplexDisplacementObject
        #thermal_real_obj = TemperatureObject
        real_obj = None
        complex_obj = None
        thermal_real_obj = None
        self.read_oug_table(data, result_name, real_obj, complex_obj, thermal_real_obj, 'node')

    def read_velocity(self, data):
        result_name = 'velocity'
        #real_obj = VelocityObject
        #complex_obj = ComplexVelocityObject
        real_obj = None
        complex_obj = None
        thermal_real_obj = None
        self.read_oug_table(data, result_name, real_obj, complex_obj, thermal_real_obj, 'node')

    def read_acceleration(self, data):
        result_name = 'acceleration'
        #real_obj = AccelerationObject
        #complex_obj = ComplexAccelerationObject
        real_obj = None
        complex_obj = None
        thermal_real_obj = None
        self.read_oug_table(data, result_name, real_obj, complex_obj, thermal_real_obj, 'node')

    def read_eigenvector(self, data):
        result_name = 'eigenvectors'
        #real_obj = EigenVectorObject
        #complex_obj = ComplexEigenVectorObject
        real_obj = None
        complex_obj = None
        thermal_real_obj = None
        self.read_oug_table(data, result_name, real_obj, complex_obj, thermal_real_obj, 'node')