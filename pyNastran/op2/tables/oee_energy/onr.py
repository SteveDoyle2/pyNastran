#pylint: disable=C0326,C0301,C0103
from __future__ import print_function, unicode_literals
from six import b
from six.moves import range
from struct import Struct, unpack
from numpy import fromstring

from pyNastran.op2.tables.oee_energy.oee_objects import RealStrainEnergy, RealStrainEnergyArray
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

        #field_num = 6
        #datai = data[4 * (field_num - 1) : 4 * (field_num + 1)]
        #assert len(datai) == 8, len(datai)
        #print(4 * (field_num - 1), 4 * (field_num + 1))
        #element_name, = self.struct_8s.unpack(data[24:32])  # changed on 11/30/2015; was this for a long time...
        element_name, = self.struct_8s.unpack(data[20:28])
        #print("element_name = %s" % (element_name))
        try:
            element_name = element_name.decode('utf-8').strip()  # element name
        except UnicodeDecodeError:
            #self.log.warning("element_name = %s" % str(element_name))
            self.log.warning("element_name - UnicodeDecodeError")
            #self.show_data(data)
            raise
        if element_name.isalnum():
            self.data_code['element_name'] = element_name
        else:
            # print("element_name = %r" % (element_name))
            # print(type(element_name))
            self.data_code['element_name'] = 'UnicodeDecodeError???'
            self.log.warning('data[20:28]=%r instead of data[24:32]' % data[20:28])

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

        #if self.data_code['element_name'] == 'BAR':
            #pass
        #else:
            #raise NotImplementedError('element_name=%r' % (
                #self.data_code['element_name']))
        result_name = 'strain_energy'
        auto_return = False
        self._results._found_result(result_name)

        if self.is_debug_file:
            self.binary_debug.write('cvalares = %s\n' % self.cvalres)
        if self.num_wide == 4:
            assert self.cvalres in [0, 1], self.cvalres

            ntotal = 16
            nelements = ndata // ntotal
            slot = getattr(self, result_name)
            if 0:
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealStrainEnergyArray)

                obj = self.obj
                if auto_return:
                    if obj.dt_temp is None or obj.itime is None and obj.dt_temp == dt:
                        element_name = self.data_code['element_name']
                        if element_name in obj.element_name_count:
                            obj.element_name_count[element_name] += nelements
                        else:
                            obj.element_name_count[element_name] = nelements
                        obj.dt_temp = dt
                    return nelements * self.num_wide * 4
                itime = obj.itime // obj.nelement_types
            else:
                if self.read_mode == 1:
                    return ndata
                self.create_transient_object(self.strain_energy, RealStrainEnergy)
                obj = self.obj
            is_vectorized = False

            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, energy, percent, density]\n')
                self.binary_debug.write('  nelements=%i\n' % nelements)

            #self.element_names.add(self.data_code['element_name'])
            if self.use_vector and is_vectorized:
                n = nelements * 4 * self.num_wide
                itotal = obj.itotal2 # was obj.itotal
                ielement2 = obj.itotal2 + nelements # was obj.itotal
                itotal2 = ielement2
                print('itime=%s itotal=%s itotal2=%s' % (itime, itotal, itotal2))

                floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 4)
                obj._times[itime] = dt
                if obj.itime == 0:
                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, 4)
                    eids = ints[:, 0] // 10
                    print('eids =', eids)
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids

                #[energy, percent, density]
                obj.data[itime, itotal:itotal2, :] = floats[:, 1:]
                obj.itotal2 = itotal2 # was obj.itotal
                #obj.ielement = ielement2
            else:
                s = Struct(b(self._endian + 'i3f'))
                for i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = s.unpack(edata)
                    (eid_device, energy, percent, density) = out
                    eid = eid_device // 10
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))
                    self.obj.add_sort1(dt, eid, energy, percent, density)
                    n += ntotal
        elif self.num_wide == 5:
            if self.read_mode == 1:
                return ndata
            assert self.cvalres in [0, 1, 2], self.cvalres # 0??
            self.create_transient_object(self.strain_energy, RealStrainEnergy)  # why is this not different?
            ntotal = 20
            nnodes = ndata // ntotal

            obj = self.obj
            is_vectorized = False
            if self.use_vector and is_vectorized:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 5)
                obj._times[obj.itime] = dt

                strings = fromstring(data, dtype=self._endian + 'S4').reshape(nelements, 5)
                s = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])
                if obj.itime == 0:
                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, 5)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids
                    obj.element_type[obj.itime, itotal:itotal2, :] = s

                #[energy, percent, density]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:]
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(b(self._endian + '8s3f'))
                for i in range(nnodes):
                    edata = data[n:n+20]
                    out = s.unpack(edata)
                    (word, energy, percent, density) = out
                    word = word.strip()
                    #print "eType=%s" % (eType)
                    #print "%s" %(self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))
                    obj.add(dt, word, energy, percent, density)
                    n += ntotal
        elif self.num_wide == 6:  ## TODO: figure this out...
            if self.read_mode == 1:
                return ndata
            self.create_transient_object(self.strain_energy, RealStrainEnergy)  # TODO: why is this not different?
            ntotal = 24
            nnodes = ndata // ntotal

            obj = self.obj
            is_vectorized = False
            if self.use_vector and is_vectorized:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 5)
                obj._times[obj.itime] = dt

                strings = fromstring(data, dtype=self._endian + 'S4').reshape(nelements, 6)
                s = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])
                if obj.itime == 0:
                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, 6)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids
                    obj.element_type[obj.itime, itotal:itotal2, :] = s

                #[energy, percent, density]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 4:]
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(b(self._endian + 'i8s3f'))
                for i in range(nnodes):
                    edata = data[n:n+24]
                    out = s.unpack(edata)
                    (word, energy, percent, density) = out  # TODO: this has to be wrong...
                    word = word.strip()
                    #print "eType=%s" % (eType)
                    #print "%s" %(self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid(out)
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))
                    obj.add(dt, word, energy, percent, density)
                    n += ntotal
        else:
            raise NotImplementedError('num_wide = %s' % self.num_wide)
        return n
