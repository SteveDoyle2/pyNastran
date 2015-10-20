#pylint: disable=C0326,C0301,C0103
from six import b
from six.moves import range
from struct import Struct, unpack

from pyNastran.op2.tables.oee_energy.oee_objects import RealStrainEnergy
from pyNastran.op2.op2_common import OP2Common

class ONR(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)
        self.words = None
        self.num_wide = None

    def _read_onr1_3(self, data, ndata):
        """
        reads ONRGY1 subtable 3
        """
        self.words = [
            'aCode',       'tCode',    'eTotal',        'isubcase',
             '???',         '???',      'element_name', 'load_set'
             'format_code', 'num_wide', 'cvalres',      'setID',
             'setID',       'eigenReal','eigenImag',     '???',
             'etotpos',     'etotneg',  '???',          '???',
             '???',         '???',      '???',      '???',
             '???', 'Title', 'subtitle', 'label']

        #aCode = self.get_block_int_entry(data, 1)

        ## total energy of all elements in isubcase/mode
        self.eTotal = self.parse_approach_code(data)
        if self.is_debug_file:
            self.binary_debug.flush()

        element_name, = unpack(b(self._endian + '8s'), data[24:32])
        print("element_name = %s" %(element_name))
        try:
            element_name = element_name.decode('utf-8').strip()  # element name
        except UnicodeDecodeError:
            print("element_name = ", str(element_name))
            #raise
        #print("element_name = %s" %(element_name))
        if element_name.isalpha():
            self.data_code['element_name'] = element_name
        else:
            self.data_code['element_name'] = 'UnicodeDecodeError???'

        #: Load set or zero
        self.load_set = self.add_data_parameter(data, 'load_set', 'i', 8, False)

        #: format code
        self.format_code = self.add_data_parameter(data, 'format_code', 'i', 9, False)

        #: number of words per entry in record
        #: .. note:: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', 'i', 10, False)
        ## C
        self.cvalres = self.add_data_parameter(data, 'cvalres', 'i', 11, False)

        #: Set identification number Number
        self.setID = self.add_data_parameter(data, 'setID', 'i', 13, False)

        #: Natural eigenvalue - real part
        self.eigenReal = self.add_data_parameter(data, 'eigenReal', 'i', 14, False)

        #: Natural eigenvalue - imaginary part
        self.eigenImag = self.add_data_parameter(data, 'eigenImag', 'i', 15, False)

        self.add_data_parameter(data, 'freq', 'f', 16, False)  ## Natural frequency

        #: Total positive energy
        self.etotpos = self.add_data_parameter(data, 'etotpos', 'f', 18)

        #: Total negative energy
        self.etotneg = self.add_data_parameter(data, 'etotneg', 'f', 19, False)

        if not self.is_sort1():
            raise NotImplementedError('sort2...')

        #self.print_block(data) # on
        if self.analysis_code == 1:   # statics / displacement / heat flux
            #del self.data_code['nonlinear_factor']
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.dataNames = self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)  ## mode number
            self.dataNames = self.apply_data_code_value('dataNames', ['mode'])
            #print "mode(5)=%s eigr(6)=%s mode_cycle(7)=%s" %(self.mode,self.eigr,self.mode_cycle)
        #elif self.analysis_code==3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code==4: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            self.freq2 = self.add_data_parameter(data, 'freq2', 'f', 5)  ## frequency
            self.dataNames = self.apply_data_code_value('dataNames', ['freq2'])
        elif self.analysis_code == 6:  # transient
            self.time = self.add_data_parameter(data, 'time', 'f', 5)  ## time step
            self.dataNames = self.apply_data_code_value('dataNames', ['time'])
        #elif self.analysis_code==7: # pre-buckling
            #self.dataNames = self.apply_data_code_value('dataNames',['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)  ## mode number
            self.dataNames = self.apply_data_code_value('dataNames', ['mode'])
        elif self.analysis_code == 9:  # complex eigenvalues
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)  ## mode number
            self.dataNames = self.apply_data_code_value('dataNames', ['mode'])
        elif self.analysis_code == 10:  # nonlinear statics
            self.loadFactor = self.add_data_parameter(data, 'loadFactor', 'f', 5)  ## load factor
            self.dataNames = self.apply_data_code_value('dataNames', ['loadFactor'])
        #elif self.analysis_code==11: # old geometric nonlinear statics
            #self.dataNames = self.apply_data_code_value('dataNames',['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.time = self.add_data_parameter(data, 'time', 'f', 5)  ## time step
            self.dataNames = self.apply_data_code_value('dataNames', ['time'])
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               self.analysis_code)

        self.fix_format_code()
        if self.is_debug_file:
            self.binary_debug.write('  approach_code  = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode          = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase       = %r\n' % self.isubcase)
        self._read_title(data)
        self._write_debug_bits()

    def _read_onr1_4(self, data, ndata):
        """
        reads ONRGY1 subtable 4
        """
        if self.read_mode == 1:
            return len(data)

        if self.table_code == 18:  # element strain energy
            assert self.table_name in [b'ONRGY1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_element_strain_energy(data, ndata)
        else:
            raise NotImplementedError(self.table_code)
        return n

    def _read_element_strain_energy(self, data, ndata):
        """
        table_code = 19
        """
        dt = self.nonlinear_factor
        n = 0
        result_name = 'strain_energy'
        if self.read_mode == 1:
            return len(data)
        self._results._found_result(result_name)

        if self.is_debug_file:
            self.binary_debug.write('cvalares = %s\n' % self.cvalres)
        if self.num_wide == 4:
            assert self.cvalres in [0, 1], self.cvalres
            self.create_transient_object(self.strain_energy, RealStrainEnergy)
            s = Struct(b(self._endian + 'i3f'))

            ntotal = 16
            nnodes = len(data) // ntotal
            for i in range(nnodes):
                edata = data[n:n+ntotal]

                out = s.unpack(edata)
                (eid_device, energy, percent, density) = out
                eid = (eid_device - self.device_code) // 10
                #print "eType=%s" % (eType)

                data_in = [eid, energy, percent, density]
                #print "%s" % (self.get_element_type(self.element_type)), data_in
                self.obj.add(dt, data_in)
                n += ntotal
        elif self.num_wide == 5:
            assert self.cvalres in [0, 1, 2], self.cvalres # 0??
            self.create_transient_object(self.strain_energy, RealStrainEnergy)  # why is this not different?
            ntotal = 20
            s = Struct(b(self._endian + '8s3f'))
            nnodes = len(data) // ntotal
            for i in range(nnodes):
                edata = data[n:n+20]
                out = s.unpack(edata)
                (word, energy, percent, density) = out
                #print "out = ",out
                word = word.strip()
                #print "eType=%s" % (eType)

                data_in = [word, energy, percent, density]
                #print "%s" %(self.get_element_type(self.element_type)), data_in
                #eid = self.obj.add_new_eid(out)
                self.obj.add(dt, data_in)
                n += ntotal
        elif self.num_wide == 6:  ## TODO: figure this out...
            self.create_transient_object(self.strain_energy, RealStrainEnergy)  # TODO: why is this not different?
            ntotal = 24
            s = Struct(b(self._endian + 'i8s3f'))
            nnodes = len(data) // ntotal
            for i in range(nnodes):
                edata = data[n:n+24]
                out = s.unpack(edata)
                (word, energy, percent, density) = out  # TODO: this has to be wrong...
                #print "out = ",out
                word = word.strip()
                #print "eType=%s" % (eType)

                data_in = [word, energy, percent, density]
                #print "%s" %(self.get_element_type(self.element_type)), data_in
                #eid = self.obj.add_new_eid(out)
                self.obj.add(dt, data_in)
                n += ntotal
        else:
            raise NotImplementedError('num_wide = %s' % self.num_wide)
        return n
