#pylint: disable=C0326,C0301
from __future__ import print_function, unicode_literals
from six import b
from six.moves import range
from struct import Struct
from numpy import fromstring, array

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
            '???',         '???',      'element_name', 'load_set',
            'format_code', 'num_wide', 'cvalres',      'setID',
            'setID',       'eigenReal', 'eigenImag',     '???',
            'etotpos',     'etotneg',  '???',          '???',
            '???',         '???',      '???',      '???',
            '???', 'Title', 'subtitle', 'label']

        #aCode = self.get_block_int_entry(data, 1)

        ## total energy of all elements in isubcase/mode
        self.etotal = self.parse_approach_code(data)
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
        self.set_id = self.add_data_parameter(data, 'set_id', 'i', 13, False)

        #: Natural eigenvalue - real part
        self.eigen_real = self.add_data_parameter(data, 'eigen_real', 'i', 14, False)

        #: Natural eigenvalue - imaginary part
        self.eigen_imag = self.add_data_parameter(data, 'eigen_imag', 'i', 15, False)

        self.add_data_parameter(data, 'freq', 'f', 16, False)  ## Natural frequency

        #: Total positive energy
        self.etotpos = self.add_data_parameter(data, 'etotpos', 'f', 18)

        #: Total negative energy
        self.etotneg = self.add_data_parameter(data, 'etotneg', 'f', 19, False)

        if not self.is_sort1():
            raise NotImplementedError('sort2...')

        if self.analysis_code == 1:   # statics / displacement / heat flux
            #del self.data_code['nonlinear_factor']
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)  ## mode number
            #self.mode_cycle1 = self.add_data_parameter(data, 'mode', 'i', 7)
            #self.mode_cycle2 = self.add_data_parameter(data, 'mode', 'f', 7)
            #print('mode = ', self.mode)
            #print('mode_cycle1 = ', self.mode_cycle1)
            #print('mode_cycle2 = ', self.mode_cycle2)
            #self.show_data(data)
            #self.cycle = 0.
            #self.update_mode_cycle('cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode'])
            #print "mode(5)=%s eigr(6)=%s mode_cycle(7)=%s" %(self.mode,self.eigr,self.mode_cycle)
        #elif self.analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            self.freq2 = self.add_data_parameter(data, 'freq2', 'f', 5)  ## frequency
            self.data_names = self.apply_data_code_value('data_names', ['freq2'])
        elif self.analysis_code == 6:  # transient
            self.time = self.add_data_parameter(data, 'time', 'f', 5)  ## time step
            self.data_names = self.apply_data_code_value('data_names', ['time'])
        #elif self.analysis_code == 7: # pre-buckling
            #self.data_names = self.apply_data_code_value('data_names',['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)  ## mode number
            self.data_names = self.apply_data_code_value('data_names', ['mode'])
        elif self.analysis_code == 9:  # complex eigenvalues
            self.mode = self.add_data_parameter(data, 'mode', 'i', 5)  ## mode number
            self.data_names = self.apply_data_code_value('data_names', ['mode'])
        elif self.analysis_code == 10:  # nonlinear statics
            self.loadFactor = self.add_data_parameter(data, 'loadFactor', 'f', 5)  ## load factor
            self.data_names = self.apply_data_code_value('data_names', ['loadFactor'])
        #elif self.analysis_code == 11: # old geometric nonlinear statics
            #self.data_names = self.apply_data_code_value('data_names',['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.time = self.add_data_parameter(data, 'time', 'f', 5)  ## time step
            self.data_names = self.apply_data_code_value('data_names', ['time'])
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
            if self.table_name not in self.table_name in [b'ONRGY', b'ONRGY1']:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise NotImplementedError(msg)
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

        if self.data_code['element_name'] == 'BAR':
            result_name = 'cbar_strain_energy'
        elif self.data_code['element_name'] == 'BEAM':
            result_name = 'cbeam_strain_energy'
        elif self.data_code['element_name'] == 'BEND':
            result_name = 'cbend_strain_energy'

        elif self.data_code['element_name'] == 'ROD':
            result_name = 'crod_strain_energy'
        elif self.data_code['element_name'] == 'TUBE':
            result_name = 'ctube_strain_energy'
        elif self.data_code['element_name'] == 'CONROD':
            result_name = 'conrod_strain_energy'


        elif self.data_code['element_name'] in ['TRIA3', 'TRIAFD', 'TRIA3FD']:
            result_name = 'ctria3_strain_energy'
        elif self.data_code['element_name'] == 'TRIA6':
            result_name = 'ctria6_strain_energy'
        elif self.data_code['element_name'] == 'TRIAX6':
            result_name = 'ctriax6_strain_energy'
        elif self.data_code['element_name'] == 'TRIAR':
            result_name = 'ctriar_strain_energy'
        elif self.data_code['element_name'] in ['TRIAX3FD', 'TRIAXFD']:
            result_name = 'ctriax_strain_energy'


        elif self.data_code['element_name'] in ['QUAD4', 'QUADFD', 'QUAD4FD']:
            result_name = 'cquad4_strain_energy'
        elif self.data_code['element_name'] == 'QUAD8':
            result_name = 'cquad8_strain_energy'
        elif self.data_code['element_name'] == 'QUADR':
            result_name = 'cquadr_strain_energy'
        elif self.data_code['element_name'] in ['QUADXFD', 'QUADX4FD']:
            result_name = 'cquadx_strain_energy'
        elif self.data_code['element_name'] == 'SHEAR':
            result_name = 'cshear_strain_energy'

        elif self.data_code['element_name'] in ['HEXA', 'HEXAFD', 'HEXA8FD']:
            result_name = 'chexa_strain_energy'
        elif self.data_code['element_name'] in ['PENTA', 'PENTAFD', 'PENTA6FD']:
            result_name = 'cpenta_strain_energy'
        elif self.data_code['element_name'] in ['TETRA', 'TETRAFD', 'TETRA4FD']:
            result_name = 'ctetra_strain_energy'
        elif self.data_code['element_name'] in ['PYRAM']:
            result_name = 'cpyram_strain_energy'

        elif self.data_code['element_name'] == 'GAP':
            result_name = 'cgap_strain_energy'
        elif self.data_code['element_name'] == 'BUSH':
            result_name = 'cbush_strain_energy'

        elif self.data_code['element_name'] == 'ELAS1':
            result_name = 'celas1_strain_energy'
        elif self.data_code['element_name'] == 'ELAS2':
            result_name = 'celas2_strain_energy'
        elif self.data_code['element_name'] == 'ELAS3':
            result_name = 'celas3_strain_energy'
        elif self.data_code['element_name'] == 'ELAS4':
            result_name = 'celas4_strain_energy'

        elif self.data_code['element_name'] == 'DUM8':
            result_name = 'cdum8_strain_energy'
        elif self.data_code['element_name'] == 'DMIG':
            result_name = 'dmig_strain_energy'
        elif self.data_code['element_name'] == 'GENEL':
            result_name = 'genel_strain_energy'
        else:
            #result_name = 'chexa8fd_strain_energy'

            raise NotImplementedError('element_name=%r' % (
                self.data_code['element_name']))
        #result_name = 'strain_energy'
        slot = getattr(self, result_name)



        auto_return = False
        self._results._found_result(result_name)

        if self.is_debug_file:
            self.binary_debug.write('cvalares = %s\n' % self.cvalres)
        if self.num_wide == 4:
            assert self.cvalres in [0, 1], self.cvalres

            ntotal = 16
            nelements = ndata // ntotal
            if 1:
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealStrainEnergyArray)

                if auto_return:
                    #if obj.dt_temp is None or obj.itime is None and obj.dt_temp == dt:
                        #element_name = self.data_code['element_name']
                        #if element_name in obj.element_name_count:
                            #obj.element_name_count[element_name] += nelements
                        #else:
                            #obj.element_name_count[element_name] = nelements
                        #obj.dt_temp = dt
                    return nelements * self.num_wide * 4
                #itime = obj.itime #// obj.nelement_types
            else:
                if self.read_mode == 1:
                    return ndata
                self.create_transient_object(slot, RealStrainEnergy)

            obj = self.obj
            itime = obj.itime
            is_vectorized = True

            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, energy, percent, density]\n')
                self.binary_debug.write('  nelements=%i\n' % nelements)

            if self.use_vector and is_vectorized:
                n = nelements * 4 * self.num_wide
                ielement = obj.ielement
                ielement2 = obj.ielement + nelements
                itotal = obj.itotal
                itotal2 = obj.itotal + nelements * 4

                floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 4)
                obj._times[itime] = dt
                if obj.itime == 0:
                    ints = fromstring(data, dtype=self.idtype).reshape(nelements, 4)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[ielement:ielement2] = eids

                #[energy, percent, density]
                obj.data[itime, ielement:ielement2, :] = floats[:, 1:]
                obj.itotal2 = itotal2
                obj.ielement = ielement2
            else:
                struct1 = Struct(b(self._endian + 'i3f'))
                for i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = struct1.unpack(edata)
                    (eid_device, energy, percent, density) = out
                    eid = eid_device // 10
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))
                    self.obj.add_sort1(dt, eid, energy, percent, density)
                    n += ntotal
        elif self.num_wide == 5:
            assert self.cvalres in [0, 1, 2], self.cvalres # 0??
            ntotal = 20
            nnodes = ndata // ntotal

            if 0:
                if self.read_mode == 1:
                    return ndata
                self.create_transient_object(slot, RealStrainEnergy)  # why is this not different?
            else:
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealStrainEnergyArray)

            obj = self.obj
            is_vectorized = True
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
            ntotal = 24
            nnodes = ndata // ntotal
            if 0:
                if self.read_mode == 1:
                    return ndata
                self.create_transient_object(slot, RealStrainEnergy)  # TODO: why is this not different?
            else:
                auto_return, is_vectorized = self._create_oes_object4(
                    nelements, result_name, slot, RealStrainEnergyArray)

            obj = self.obj
            is_vectorized = True
            if self.use_vector and is_vectorized:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = fromstring(data, dtype=self.fdtype).reshape(nelements, 5)
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    strings = fromstring(data, dtype=self._endian + 'S4').reshape(nelements, 6)
                    s = array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])

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
                struct1 = Struct(b(self._endian + 'i8s3f'))
                for i in range(nnodes):
                    edata = data[n:n+24]
                    out = struct1.unpack(edata)
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
