from struct import unpack

from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import(
    SPCForcesObject, ComplexSPCForcesObject)
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import(
    MPCForcesObject, ComplexMPCForcesObject)
from pyNastran.op2.tables.oqg_constraintForces.oqg_thermalGradientAndFlux import (
    TemperatureGradientAndFluxObject)

class OQG(object):
    """Table of spc/mpc forces/momenets"""

    def readTable_OQG(self):
        table3 = self.readTable_OQG_3
        table4Data = self.readOQG_Data
        self.read_results_table(table3, table4Data)
        self._delete_attributes_OQG()

    def _delete_attributes_OQG(self):
        #print self.obj
        params = ['lsdvm', 'mode', 'eigr', 'mode_cycle', 'freq', 'dt', 'lftsfq', 'thermal', 'rCode', 'fCode', 'num_wide', 'acousticFlag', 'thermal']
        self._delete_attributes(params)

    def readTable_OQG_3(self, iTable):  # iTable=-3
        buffer_words = self.get_marker()
        if self.make_op2_debug:
            self.op2Debug.write('buffer_words=%s\n' % (str(buffer_words)))
        #print "2-buffer_words = ",buffer_words,buffer_words*4,'\n'

        data = self.get_data(4)
        buffer_size, = unpack('i', data)
        data = self.get_data(4 * 50)
        #print self.print_block(data)

        (three) = self.parse_approach_code(data)

        ## random code
        self.add_data_parameter( data, 'randomCode', 'i', 8, False)

        ## format code
        self.add_data_parameter( data, 'format_code', 'i', 9, False)

        ## number of words per entry in record
        self.add_data_parameter(data, 'num_wide', 'i', 10, False)

        ## acoustic pressure flag
        self.add_data_parameter(data, 'acousticFlag', 'f', 13, False)

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
        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        #print "*isubcase=%s"%(self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)
        #print self.code_information()

        #self.print_block(data)
        self.read_title()

    def readOQG_Data(self):
        #tfsCode = [self.table_code,self.format_code,self.sort_code]

        #print "self.analysis_code=%s table_code(1)=%s thermal(23)=%g" %(self.analysis_code,self.table_code,self.thermal)
        assert self.thermal in [0, 1, 8], self.code_information()

        if   self.table_code == 3:   # SPC Forces
            assert self.table_name in ['OQG1', 'OQGV1', 'OQP1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            self.readOQG_Data_table3()
        elif self.table_code == 39:  # MPC Forces
            assert self.table_name in ['OQMG1'], 'table_name=%s table_code=%s' % (
                self.table_name, self.table_code)
            self.readOQG_Data_table3()
        else:
            self.not_implemented_or_skip('bad OQG table')
        #print self.obj

    def readOQG_Data_table3(self):  # SPC Forces
        #is_sort1 = self.is_sort1()
        #print(self.code_information())
        magPhase = self.is_magnitude_phase()
        if magPhase or self.num_wide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'spcForces'
                self.create_transient_object(self.spcForces,
                                           ComplexSPCForcesObject)
                self.handle_results_buffer(self.OUG_ComplexTable, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.num_wide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'spcForces'
                self.create_transient_object(
                    self.spcForces, SPCForcesObject)
                self.handle_results_buffer(self.OUG_RealTable, resultName)
            elif self.thermal == 1:
                resultName = 'thermalGradientAndFlux' #'finite element temperature gradients and fluxes'
                self.create_transient_object(
                    self.thermalGradientAndFlux, TemperatureGradientAndFluxObject)
                self.handle_results_buffer(self.OUG_RealTable, resultName)
            else:
                self.not_implemented_or_skip()
        else:
            self.not_implemented_or_skip('only num_wide=8 or 14 is allowed  num_wide=%s' % (self.num_wide))

        #if self.thermal not in [0,1]:
            #print self.obj
            #raise RuntimeError('check the printout for thermal...')

    def readOQG_Data_table39(self):  # MPC Forces
        #is_sort1 = self.is_sort1()
        if self.num_wide == 8:  # real/random
            if self.thermal == 0:
                resultName = 'mpcForces'
                self.create_transient_object(
                    self.mpcForces, MPCForcesObject)
                self.handle_results_buffer(self.OUG_RealTable, resultName)
            else:
                self.not_implemented_or_skip()
        elif self.num_wide == 14:  # real/imaginary or mag/phase
            if self.thermal == 0:
                resultName = 'mpcForces'
                self.create_transient_object(self.mpcForces,
                                           ComplexMPCForcesObject)
                self.handle_results_buffer(self.OUG_ComplexTable, resultName)
            else:
                self.not_implemented_or_skip()
        else:
            self.not_implemented_or_skip('only num_wide=8 or 14 is allowed  num_wide=%s' % (self.num_wide))
