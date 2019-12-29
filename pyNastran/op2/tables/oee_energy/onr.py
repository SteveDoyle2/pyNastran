#pylint: disable=C0326,C0301
from struct import Struct
import numpy as np

from pyNastran.op2.tables.oee_energy.oee_objects import RealStrainEnergyArray, ComplexStrainEnergyArray
from pyNastran.op2.op2_interface.op2_common import OP2Common
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block

RESULT_NAME_MAP = {
    'BAR' : 'cbar_strain_energy',
    'BEAM' : 'cbeam_strain_energy',
    'BEND' : 'cbend_strain_energy',

    'ROD' : 'crod_strain_energy',
    'TUBE' : 'ctube_strain_energy',
    'CONROD' : 'conrod_strain_energy',

    'TRIA3' : 'ctria3_strain_energy',
    'TRIAFD' : 'ctria3_strain_energy',
    'TRIA3FD' : 'ctria3_strain_energy',
    'TRIA6' : 'ctria6_strain_energy',
    'TRIAX6' : 'ctriax6_strain_energy',
    'TRIAR' : 'ctriar_strain_energy',
    'TRIAX3FD' : 'ctriax_strain_energy',
    'TRIAXFD' : 'ctriax_strain_energy',

    'QUAD4' : 'cquad4_strain_energy',
    'QUADFD' : 'cquad4_strain_energy',
    'QUAD4FD' : 'cquad4_strain_energy',
    'QUAD8' : 'cquad8_strain_energy',

    'QUADR' : 'cquadr_strain_energy',
    'QUADXFD' : 'cquadx_strain_energy',
    'QUADX4FD' : 'cquadx_strain_energy',
    'SHEAR' : 'cshear_strain_energy',

    'TETRA' : 'ctetra_strain_energy',
    'TETRAFD' : 'ctetra_strain_energy',
    'TETRA4FD' : 'ctetra_strain_energy',
    'PENTA' : 'cpenta_strain_energy',
    'PENTAFD' : 'cpenta_strain_energy',
    'PENTA6FD' : 'cpenta_strain_energy',
    'HEXA' : 'chexa_strain_energy',
    'HEXAFD' : 'chexa_strain_energy',
    'HEXA8FD' : 'chexa_strain_energy',
    'PYRAM' : 'cpyram_strain_energy',

    'GAP' : 'cgap_strain_energy',
    'BUSH' : 'cbush_strain_energy',
    'ELAS1' : 'celas1_strain_energy',
    'ELAS2' : 'celas2_strain_energy',
    'ELAS3' : 'celas3_strain_energy',
    'ELAS4' : 'celas4_strain_energy',

    'DUM8' : 'cdum8_strain_energy',
    'DMIG' : 'dmig_strain_energy',
    'GENEL' : 'genel_strain_energy',
    'CONM2' : 'conm2_strain_energy',
    'RBE1' : 'rbe1_strain_energy',
    'RBE3' : 'rbe3_strain_energy',
    'WELDC' : 'weldc_strain_energy',
}

class ONR(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)
        self.words = None
        self.num_wide = None

    def get_onr_prefix_postfix(self):
        """
        Creates the prefix/postfix that splits off ATO, CRM, PSD, nonlinear,
        etc. results.  We also fix some of the sort bits as typing:

            STRESS(PLOT,SORT1,RALL) = ALL

        will actually create the OESRMS2 table (depending on what else
        is in your case control).  However, it's in an OESATO2 table, so
        we know it's really SORT2.

        Also, if you're validating the sort_bit flags, *RMS2 and *NO2 are
        actually SORT1 tables.

        NX Case Control  Block         Description
        ===============  ==========    ===========
        NLSTRESS         OESNLXR       Nonlinear static stresses
        BOUTPUT          OESNLBR       Slideline stresses
        STRESS           OESNLXD       Nonlinear Transient Stresses
        STRESS           OES1C/OSTR1C  Ply stresses/strains
        STRESS           OES1X         Element stresses with intermediate (CBAR and CBEAM)
                                       station stresses and stresses on nonlinear elements
        STRESS           OES/OESVM     Element stresses (linear elements only)
        STRAIN           OSTR1         Element strains
        STRESS/STRAIN    DOES1/DOSTR1  Scaled Response Spectra
        MODCON           OSTRMC        Modal contributions
        """
        prefix = ''
        postfix = ''
        if self.table_name in [b'ONRGY1', b'ONRGY2', b'ONRGY']:
            prefix = 'strain_energy.'
        elif self.table_name in [b'RANEATC']: #, b'OSTRMS1C']:
            self.format_code = 1
            self.sort_bits[0] = 0 # real
            prefix = 'RANEATC.'
        elif self.table_name in [b'RANCONS']: #, b'OSTRMS1C']:
            self.format_code = 1
            self.sort_bits[0] = 0 # real
            prefix = 'RANCONS.'
        else:
            raise NotImplementedError(self.table_name)
        self.data_code['sort_bits'] = self.sort_bits
        self.data_code['nonlinear_factor'] = self.nonlinear_factor
        return prefix, postfix

    def _read_onr1_3(self, data, ndata):
        """
        reads ONRGY1 subtable 3
        """
        self._analysis_code_fmt = b'i'
        self.words = [
            'aCode',       'tCode',    'eTotal',       'isubcase',
            '???',         '???',      'element_name', 'load_set',
            'format_code', 'num_wide', 'cvalres',      'setID',
            'setID',       'eigenReal', 'eigenImag',   'rmssf',
            'etotpos',     'etotneg',  'thresh',       '???',
            '???',         '???',      '???',      '???',
            '???', 'Title', 'subtitle', 'label']

        #aCode = self.get_block_int_entry(data, 1)

        ## total energy of all elements in isubcase/mode
        self.etotal = self.parse_approach_code(data)
        if self.is_debug_file:
            self.binary_debug.flush()

        self._onr_element_name(data)

        #: Load set or zero
        self.load_set = self.add_data_parameter(data, 'load_set', b'i', 8, False)

        #: format code
        self.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        #: number of words per entry in record
        #: .. note:: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)
        ## C
        self.cvalres = self.add_data_parameter(data, 'cvalres', b'i', 11, False)

        #: Set identification number Number
        self.set_id = self.add_data_parameter(data, 'set_id', b'i', 13, False)

        #: Natural eigenvalue - real part
        self.eigen_real = self.add_data_parameter(data, 'eigen_real', b'i', 14, False)

        #: Natural eigenvalue - imaginary part
        self.eigen_imag = self.add_data_parameter(data, 'eigen_imag', b'i', 15, False)

        #: Natural frequency
        self.freq = self.add_data_parameter(data, 'freq', b'f', 16, False)

        #: RMS and CRMS scale factor - NX
        self.rmssf = self.add_data_parameter(data, 'rmssf', b'f', 17)

        #: Total positive energy
        self.etotpos = self.add_data_parameter(data, 'etotpos', b'f', 18)

        #: Total negative energy
        self.etotneg = self.add_data_parameter(data, 'etotneg', b'f', 19, False)

        #: Energy Threshold - NX
        self.thresh = self.add_data_parameter(data, 'thresh', b'f', 17)

        if not self.is_sort1:
            raise NotImplementedError('sort2...')

        if self.analysis_code == 1:   # statics / displacement / heat flux
            #del self.data_code['nonlinear_factor']
            self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)  ## mode number
            #self.mode_cycle1 = self.add_data_parameter(data, 'mode', b'i', 7)
            #self.mode_cycle2 = self.add_data_parameter(data, 'mode', b'f', 7)
            #print('mode = ', self.mode)
            #print('mode_cycle1 = ', self.mode_cycle1)
            #print('mode_cycle2 = ', self.mode_cycle2)
            #self.show_data(data)
            #self.cycle = 0.
            #self.update_mode_cycle('cycle')
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'freq'])
            #print("mode(5)=%s eign(6)=%s mode_cycle(7)=%s" % (
                #self.mode, self.eign, self.mode_cycle))
        #elif self.analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            self.freq2 = self.add_data_parameter(data, 'freq2', b'f', 5)  ## frequency
            self.data_names = self.apply_data_code_value('data_names', ['freq2'])
        elif self.analysis_code == 6:  # transient
            self.time = self.add_data_parameter(data, 'time', b'f', 5)  ## time step
            self.data_names = self.apply_data_code_value('data_names', ['time'])
        #elif self.analysis_code == 7: # pre-buckling
            #self.data_names = self.apply_data_code_value('data_names',['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)  ## mode number
            self.data_names = self.apply_data_code_value('data_names', ['mode'])
        elif self.analysis_code == 9:  # complex eigenvalues
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)  ## mode number
            self.eigr = self.eigen_real
            self.eigi = self.eigen_imag
            self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'eign'])
        elif self.analysis_code == 10:  # nonlinear statics
            self.loadFactor = self.add_data_parameter(data, 'loadFactor', b'f', 5)  ## load factor
            self.data_names = self.apply_data_code_value('data_names', ['loadFactor'])
        #elif self.analysis_code == 11: # old geometric nonlinear statics
            #self.data_names = self.apply_data_code_value('data_names',['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.time = self.add_data_parameter(data, 'time', b'f', 5)  ## time step
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

    def _onr_element_name(self, data):
        #field_num = 6
        #datai = data[4 * (field_num - 1) : 4 * (field_num + 1)]
        #assert len(datai) == 8, len(datai)
        #print(4 * (field_num - 1), 4 * (field_num + 1))
        #element_name, = self.struct_8s.unpack(data[24:32])  # changed on 11/30/2015; was this for a long time...

        #self.show_data(data[:28])
        if self.size == 4:
            element_name, = self.struct_8s.unpack(data[20:28])
        else:
            element_name, = self.struct_16s.unpack(data[40:56])
            element_name = reshape_bytes_block(element_name)

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
            #print("element_name = %r" % (element_name))
            self.data_code['element_name'] = 'UnicodeDecodeError???'
            self.log.warning('data[20:28]=%r instead of data[24:32]' % data[20:28])

    def _read_onr2_3(self, data, ndata):
        """reads the SORT2 version of table 4 (the data table)"""
        self.nonlinear_factor = np.nan
        self.is_table_1 = False
        self.is_table_2 = True
        unused_three = self.parse_approach_code(data)
        self.words = [
            'aCode',       'tCode',    'eTotal',       'isubcase',
            '???',         '???',      'element_name', 'load_set',
            'format_code', 'num_wide', 'cvalres',      'setID',
            'setID',       'eigenReal', 'eigenImag',   'rmssf',
            'etotpos',     'etotneg',  'thresh',       '???',
            '???',         '???',      '???',      '???',
            '???', 'Title', 'subtitle', 'label']

        self._onr_element_name(data)

        #: Load set or zero
        self.load_set = self.add_data_parameter(data, 'load_set', b'i', 8, False)

        #: format code
        self.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        #: number of words per entry in record
        #: .. note:: is this needed for this table ???
        self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)
        ## C
        self.cvalres = self.add_data_parameter(data, 'cvalres', b'i', 11, False)

        #: Set identification number Number
        self.set_id = self.add_data_parameter(data, 'set_id', b'i', 13, False)

        #: Natural eigenvalue - real part
        self.eigen_real = self.add_data_parameter(data, 'eigen_real', b'i', 14, False)

        #: Natural eigenvalue - imaginary part
        self.eigen_imag = self.add_data_parameter(data, 'eigen_imag', b'i', 15, False)

        #: Natural frequency
        self.freq = self.add_data_parameter(data, 'freq', b'f', 16, False)

        #: RMS and CRMS scale factor - NX
        self.rmssf = self.add_data_parameter(data, 'rmssf', b'f', 17)

        #: Total positive energy
        self.etotpos = self.add_data_parameter(data, 'etotpos', b'f', 18)

        #: Total negative energy
        self.etotneg = self.add_data_parameter(data, 'etotneg', b'f', 19, False)

        #: Energy Threshold - NX
        self.thresh = self.add_data_parameter(data, 'thresh', b'f', 17)

        self.element_id = self.add_data_parameter(data, 'node_id', b'i', 5, fix_device_code=True)
        #if self.analysis_code == 1:  # statics / displacement / heat flux
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            #self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            #self.setNullNonlinearFactor()

        if self.analysis_code == 1:  # static...because reasons.
            self._analysis_code_fmt = b'i'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'N/A')
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            self._analysis_code_fmt = b'i'
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', b'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names',
                                                         ['node_id', 'eigr', 'mode_cycle'])
            self.apply_data_code_value('analysis_method', 'mode')
        #elif self.analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #self.data_names = self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            #self.freq = self.add_data_parameter(data, 'freq', b'f', 5)
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'freq')
        elif self.analysis_code == 6:  # transient
            ## time step
            #self.dt = self.add_data_parameter(data, 'dt', b'f', 5)
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'dt')
        elif self.analysis_code == 7:  # pre-buckling
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self._analysis_code_fmt = b'i'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'lsdvmn')
        elif self.analysis_code == 8:  # post-buckling
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self._analysis_code_fmt = b'i'
            ## real eigenvalue
            self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr'])
            self.apply_data_code_value('analysis_method', 'eigr')
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            self._analysis_code_fmt = b'i'
            ## real eigenvalue
            #self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            self.eigi = self.add_data_parameter(data, 'eigi', b'f', 7, False)
            self.data_names = self.apply_data_code_value('data_names', ['node_id', 'eigr', 'eigi'])
            self.apply_data_code_value('analysis_method', 'mode')
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            #self.lftsfq = self.add_data_parameter(data, 'lftsfq', b'f', 5)
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'lftsfq')
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self._analysis_code_fmt = b'f'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
        elif self.analysis_code == 12:
            # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            self._analysis_code_fmt = b'i'
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            self.apply_data_code_value('analysis_method', 'lsdvmn')
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
            raise RuntimeError(msg)

        self.fix_format_code()
        if self.num_wide == 8:
            self.format_code = 1
            self.data_code['format_code'] = 1
        else:
            #self.fix_format_code()
            if self.format_code == 1:
                self.format_code = 2
                self.data_code['format_code'] = 2
            assert self.format_code in [2, 3], self.code_information()

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

        element_name = self.data_code['element_name']
        try:
            result_name = RESULT_NAME_MAP[element_name]
        except KeyError:
            raise NotImplementedError('element_name=%r' % (self.data_code['element_name']))
        prefix, postfix = self.get_onr_prefix_postfix()
        result_name = prefix + result_name + postfix
        #result_name = 'strain_energy'

        if self._results.is_not_saved(result_name):
            return ndata
        self._results._found_result(result_name)

        slot = self.get_result(result_name)

        #auto_return = False
        if self.is_debug_file:
            self.binary_debug.write('cvalares = %s\n' % self.cvalres)
        if self.format_code in [1, 2] and self.num_wide == 4:
            assert self.cvalres in [0, 1], self.cvalres

            assert self.num_wide == 4
            ntotal = 16 * self.factor  # 4*4=16
            nelements = ndata // ntotal
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
                return nelements * ntotal
            #itime = obj.itime #// obj.nelement_types

            obj = self.obj
            itime = obj.itime

            if self.is_debug_file:
                self.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                self.binary_debug.write('  #elementi = [eid_device, energy, percent, density]\n')
                self.binary_debug.write('  nelements=%i\n' % nelements)

            if self.use_vector and self.sort_method == 1: # and self.is_sort1:
                n = nelements * ntotal
                ielement = obj.ielement
                ielement2 = obj.ielement + nelements
                itotal = obj.itotal
                itotal2 = obj.itotal + nelements * 4

                floats = np.frombuffer(data, dtype=self.fdtype8).reshape(nelements, 4)
                obj._times[itime] = dt
                #if obj.itime == 0:
                ints = np.frombuffer(data, dtype=self.idtype8).reshape(nelements, 4)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, f'etype={self.element_name} isubtable={self.isubtable} eids.min()={eids.min()}'
                obj.element[itime, ielement:ielement2] = eids

                #[energy, percent, density]
                obj.data[itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal2 = itotal2
                obj.ielement = ielement2
            else:
                fmt = mapfmt(self._endian + self._analysis_code_fmt + b'3f', self.size)
                struct1 = Struct(fmt)
                for unused_i in range(nelements):
                    edata = data[n:n+ntotal]

                    out = struct1.unpack(edata)
                    (eid_device, energy, percent, density) = out
                    if self.sort_method == 1:
                        eid = eid_device // 10
                    else:
                        eid = self.nonlinear_factor
                        dt = eid_device
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))
                    self.obj.add_sort1(dt, eid, energy, percent, density)
                    n += ntotal
        elif self.format_code == 1 and self.num_wide == 5:
            assert self.cvalres in [0, 1, 2], self.cvalres # 0??
            ntotal = 20
            nnodes = ndata // ntotal
            nelements = nnodes

            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealStrainEnergyArray)
            if auto_return:
                return nelements * self.num_wide * 4

            obj = self.obj
            if self.use_vector:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = np.frombuffer(data, dtype=self.fdtype).reshape(nelements, 5).copy()
                obj._times[obj.itime] = dt

                strings = np.frombuffer(data, dtype=self._uendian + 'S4').reshape(nelements, 5)
                if obj.itime == 0:
                    ints = np.frombuffer(data, dtype=self.idtype).reshape(nelements, 5)
                    if obj.element_name == 'DMIG':
                        s = np.array([(s1+s2).decode('latin1').strip()
                                      for s1, s2 in zip(strings[:, 0], strings[:, 1])], dtype='|U8')
                        obj.element[itotal:itotal2] = s
                    else:
                        eids = ints[:, 0] // 10
                        assert eids.min() > 0, eids.min()

                        s = np.array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])
                        obj.element[itotal:itotal2] = eids
                        obj.element_type[obj.itime, itotal:itotal2, :] = s

                #[energy, percent, density]
                if obj.element_name == 'DMIG':
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 2:]
                else:
                    obj.data[obj.itime, itotal:itotal2, :] = floats[:, 3:]
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(self._endian + b'8s3f')
                for unused_i in range(nnodes):
                    edata = data[n:n+20]
                    out = s.unpack(edata)
                    (word, energy, percent, density) = out
                    word = word.strip()
                    if self.is_debug_file:
                        self.binary_debug.write('  eid/word=%r; %s\n' % (word, str(out)))
                    obj.add_sort1(dt, word, energy, percent, density)
                    n += ntotal
        elif self.format_code in [2, 3] and self.num_wide == 5:
            #ELEMENT-ID   STRAIN-ENERGY (MAG/PHASE)  PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
            #    5         2.027844E-10 /   0.0            1.2581            2.027844E-09
            ntotal = 20
            nelements = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, ComplexStrainEnergyArray)
            if auto_return:
                return nelements * self.num_wide * 4

            obj = self.obj
            if self.use_vector:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = np.frombuffer(data, dtype=self.fdtype).reshape(nelements, 5)
                obj._times[obj.itime] = dt

                #if obj.itime == 0:
                ints = np.frombuffer(data, dtype=self.idtype).reshape(nelements, 5)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids
                #obj.element_type[obj.itime, itotal:itotal2, :] = s

                #[energyr, energyi, percent, density]
                obj.element[obj.itime, itotal:itotal2] = eids
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                s = Struct(self._endian + b'i4f')
                for unused_i in range(nelements):
                    edata = data[n:n+20]
                    out = s.unpack(edata)
                    (eid_device, energyr, energyi, percent, density) = out
                    eid = eid_device // 10
                    #if is_magnitude_phase:
                        #energy = polar_to_real_imag(energyr, energyi)
                    #else:
                        #energy = complex(energyr, energyi)

                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))
                    obj.add_sort1(dt, eid, energyr, energyi, percent, density)
                    n += ntotal

        elif self.format_code == 1 and self.num_wide == 6:  ## TODO: figure this out...
            ntotal = 24
            nnodes = ndata // ntotal
            auto_return, is_vectorized = self._create_oes_object4(
                nelements, result_name, slot, RealStrainEnergyArray)

            if auto_return:
                return nelements * self.num_wide * 4

            obj = self.obj
            if self.use_vector:
                n = nelements * 4 * self.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = np.frombuffer(data, dtype=self.fdtype).reshape(nelements, 5)
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    strings = np.frombuffer(data, dtype=self._uendian + 'S4').reshape(nelements, 6)
                    s = np.array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])

                    ints = np.frombuffer(data, dtype=self.idtype).reshape(nelements, 6)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids
                    obj.element_type[obj.itime, itotal:itotal2, :] = s

                #[energy, percent, density]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 4:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                struct1 = Struct(self._endian + b'i8s3f')
                for unused_i in range(nnodes):
                    edata = data[n:n+24]
                    out = struct1.unpack(edata)
                    (word, energy, percent, density) = out  # TODO: this has to be wrong...
                    word = word.strip()
                    #print "eType=%s" % (eType)
                    #print "%s" %(self.get_element_type(self.element_type)), data_in
                    #eid = self.obj.add_new_eid_sort1(out)
                    if self.is_debug_file:
                        self.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))
                    obj.add_sort1(dt, word, energy, percent, density)
                    n += ntotal
        elif self.format_code in [2, 3] and self.num_wide == 4:
            #
            #  FREQUENCY =  1.000000E+01
            #                           E L E M E N T   S T R A I N   E N E R G I E S   ( A V E R A G E )
            #
            #            ELEMENT-TYPE = QUADR               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =   3.662188E+06
            #            SUBCASE               1            * TOTAL ENERGY OF ALL ELEMENTS IN SET       9 =   1.853189E+05
            #
            #                                ELEMENT-ID          STRAIN-ENERGY           PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY
            #                                         9          8.723258E+03                 0.2382              5.815505E+00
            #                                        10          7.815898E+03                 0.2134              5.210599E+00
            #                                        11          8.512115E+04                 2.3243              5.674743E+01
            #                                        12          5.200864E+04                 1.4202              3.467243E+01
            #
            #                    TYPE = QUADR    SUBTOTAL        1.536690E+05                 4.1961
            #
            #device_code   = 1   Print
            #analysis_code = 5   Frequency
            #table_code    = 18  ONRGY1-OEE - Element strain energy
            #format_code   = 2   Real/Imaginary
            #sort_method   = 1
            #sort_code     = 0
                #sort_bits   = (0, 0, 0)
                #data_format = 0   Real
                #sort_type   = 0   Sort1
                #is_random   = 0   Sorted Responses
            #random_code   = 0
            #s_code        = None ???
            #num_wide      = 4
            #isubcase      = 1
            #MSC Nastran
            raise NotImplementedError('onr')
        else:
            msg = self.code_information()
            return self._not_implemented_or_skip(data, ndata, msg)
            #raise NotImplementedError(self.code_information())
        return n
