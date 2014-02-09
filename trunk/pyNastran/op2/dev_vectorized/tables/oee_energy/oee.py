from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from struct import unpack

# pyNastran
from .oee_objects import StrainEnergyObject


class OEE(object):
    """Table of energy"""

    def readTable_OEE(self):
        table3 = self.read_table_OEE_3
        table4_data = self.read_OEE_data
        self.read_results_table(table3, table4_data)
        self._delete_attributes_OEE()

    def _delete_attributes_OEE(self):  # no thermal
        params = ['lsdvm', 'mode', 'eigr', 'freq', 'dt', 'lftsfq',
                  'format_code', 'num_wide']
        self._delete_attributes(params)

    def read_table_OEE_3(self, iTable):  # iTable=-3
        buffer_words = self.get_marker()
        if self.make_op2_debug:
            self.op2Debug.write('buffer_words=%s\n' % str(buffer_words))
        #print "2-buffer_words = ",buffer_words,buffer_words*4,'\n'

        data = self.get_data(4)
        buffer_size, = unpack(b'i', data)
        data = self.get_data(4 * 50)
        #print self.print_block(data)

        aCode = self.get_block_int_entry(data, 1)
        ## total energy of all elements in isubcase/mode
        self.eTotal = self._parse_approach_code(data)
        #print(self.print_section(100))
        element_name, = unpack(b'8s', data[24:32])
        #print("element_name = %s" %(element_name))
        try:
            element_name = element_name.decode('utf-8').strip()  # element name
        except UnicodeDecodeError:
            print("element_name = ", str(element_name))
            raise
        #print("element_name = %s" % element_name)
        if element_name.isalpha():
            self.data_code['element_name'] = element_name

        ## Load set or zero
        self._add_data_parameter(data, 'load_set', 'i', 8, False)
        ## format code
        self._add_data_parameter(data, 'format_code', 'i', 9, False)
        #: number of words per entry in record
        #: .. note:: is this needed for this table ???
        self._add_data_parameter(data, 'num_wide', 'i', 10, False)
        ## C
        self._add_data_parameter(data, 'cvalres', 'i', 11, False)

        ## Set identification number Number
        self._add_data_parameter(data, 'setID', 'i', 13, False)

        self._add_data_parameter(data, 'eigenReal', 'i', 14, False)
            ## Natural eigenvalue - real part
        self._add_data_parameter(data, 'eigenImag', 'i', 15, False)
            ## Natural eigenvalue - imaginary part
        self._add_data_parameter(
            data, 'freq', 'f', 16, False)  ## Natural frequency
        self._add_data_parameter(data, 'etotpos', 'f', 18)
            ## Total positive energy
        self._add_data_parameter(data, 'etotneg', 'f', 19, False)
            ## Total negative energy

        if not self.is_sort1():
            raise NotImplementedError('sort2...')

        #self.print_block(data) # on
        if self.analysis_code == 1:   # statics / displacement / heat flux
            #del self.data_code['nonlinear_factor']
            self._apply_data_code_value('dataNames', ['lsdvmn'])
            self._set_null_nonlinear_factor()
        elif self.analysis_code == 2:  # real eigenvalues
            self._add_data_parameter(data, 'mode', 'i', 5)  ## mode number
            self._apply_data_code_value('dataNames', ['mode'])
            #print "mode(5)=%s eigr(6)=%s mode_cycle(7)=%s" %(self.mode,self.eigr,self.mode_cycle)
        #elif self.analysis_code==3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code==4: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            self._add_data_parameter(data, 'freq2', 'f', 5)  ## frequency
            self._apply_data_code_value('dataNames', ['freq2'])
        elif self.analysis_code == 6:  # transient
            self._add_data_parameter(data, 'time', 'f', 5)  ## time step
            self._apply_data_code_value('dataNames', ['time'])
        #elif self.analysis_code==7: # pre-buckling
            #self.apply_data_code_value('dataNames',['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            self._add_data_parameter(data, 'mode', 'i', 5)  ## mode number
            self._apply_data_code_value('dataNames', ['mode'])
        elif self.analysis_code == 9:  # complex eigenvalues
            self._add_data_parameter(data, 'mode', 'i', 5)  ## mode number
            self._apply_data_code_value('dataNames', ['mode'])
        elif self.analysis_code == 10:  # nonlinear statics
            self._add_data_parameter(data, 'loadFactor', 'f', 5)  ## load factor
            self._apply_data_code_value('dataNames', ['loadFactor'])
        #elif self.analysis_code==11: # old geometric nonlinear statics
            #self.apply_data_code_value('dataNames',['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self._add_data_parameter(data, 'time', 'f', 5)  ## time step
            self._apply_data_code_value('dataNames', ['time'])
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               self.analysis_code)

        #print "*isubcase=%s element_name=|%s|" % (self.isubcase,self.element_name)
        #print "analysis_code=%s table_code=%s" % (self.analysis_code,self.table_code)
        #print self.code_information()

        #self.print_block(data)
        self.read_title()

    def read_OEE_data(self):
        #print "self.analysiscode=%s tablecode(1)=%s" %(self.analysiscode,self.tablecode)
        #tfsCode = [self.tablecode,self.formatcode,self.sortcode]

        if self.table_code == 18:
            assert self.table_name in ['ONRGY1', 'ONRGY2'], 'table_name=%s tablecode=%s' % (self.table_name, self.tablecode)
            self.read_strain_energy_table_18()
        else:
            self.not_implemented_or_skip('bad OEE table')
        #print str(self.obj)

    def read_strain_energy_table_18(self):  # real ???
        self.create_transient_object(self.strainEnergy, StrainEnergyObject)
        resultName = 'strainEnergy'
        name = resultName + ': Subcase %s' % self.isubcase

        if self.num_wide == 4:
            self.handle_results_buffer(self.OEE_Strain4, resultName, name)
        elif self.num_wide == 5:
            self.handle_results_buffer(self.OEE_Strain5, resultName, name)
        else:
            self.not_implemented_or_skip()
        #self.readMappedScalarsOut(debug=False) # handles dtMap, not correct...

    def OEE_Strain4(self, name):
        #device_code = self.device_code
        dt = self.nonlinear_factor

        (format1, extract) = self.getOUG_FormatStart()  # TODO change to OEE
        format1 += 'fff'
        format1 = bytes(format1)

        npoints = len(self.data) // 16
        if self.read_mode == 0 or name not in self._selected_names:
            if name not in self._result_names:
                self._result_names.append(name)
            # figure out the shape
            #self.obj._increase_size(dt, nnodes)
            iend = 16 * (npoints + 1)
            self.data = self.data[iend:]
            return
        else:  # read_mode = 1; # we know the shape so we can make a pandas matrix
            #index_data = (nodeIDs_to_index)
            #index_data = pd.MultiIndex.from_tuples(zip(data))

            #(inode_start, inode_end) = self.obj._preallocate(dt, nnodes)
            pass

        n0 = 0
        n1 = 16
        for i in xrange(npoints):
            eData = self.data[n0:n1]
            (eid, energy, percent, density) = unpack(format1, eData)
            eid2 = extract(eid, dt)
            #print("eType=%s" % eType)

            dataIn = [eid2, energy, percent, density]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, dataIn)
            n0 = n1
            n1 += 16
        self.data = self.data[:n0]
        #print('n0=%s n1=%s' % (n0, n1))
        #print self.strainEnergy

    def OEE_Strain5(self, name):
        #device_code = self.device_code
        dt = self.nonlinear_factor

        #(format1,extract) = self.getOUG_FormatStart()  # TODO change to OEE
        format1 = b'8s3f'

        npoints = len(self.data) // 20
        n0 = 0
        n1 = 20
        for i in xrange(npoints):
        #while len(self.data) >= 16:  # 5*4
            eData = self.data[n0:n1]

            (word, energy, percent, density) = unpack(format1, eData)
            #print "out = ",out
            word = word.strip()
            #print("eType=%s" % eType)

            dataIn = [word, energy, percent, density]
            #print "%s" %(self.get_element_type(self.element_type)),dataIn
            #eid = self.obj.add_new_eid(out)
            self.obj.add(dt, dataIn)
            n0 = n1
            n1 += 20
        self.data = self.data[:n0]
        #print self.strainEnergy