#pylint: disable=C0326,C0301
from __future__ import annotations
from typing import TYPE_CHECKING
from struct import Struct
import numpy as np

from pyNastran.op2.tables.oee_energy.oee_objects import RealStrainEnergyArray, ComplexStrainEnergyArray
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

RESULT_NAME_MAP = {
    'BAR' : 'cbar_strain_energy',
    'BEAM' : 'cbeam_strain_energy',
    'BEND' : 'cbend_strain_energy',
    'BEAM3' : 'cbeam3_strain_energy',

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
    # TODO: this will probably be a problem someday...cquad8_nonlinear_strain_energy
    'QUAD8N' : 'cquad8_strain_energy',

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
    'PYRA' : 'cpyram_strain_energy',

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
    'WELDP' : 'cweld_strain_energy',
    'WELDC' : 'cweld_strain_energy',
    'WELD' : 'cweld_strain_energy',
    'FASTP' : 'cfast_strain_energy',
    'SEAMP' : 'cseam_strain_energy',
}

class ONR:
    def __init__(self, op2: OP2):
        self.op2 = op2
        #op2.words = None
        #op2.num_wide = None

    @property
    def size(self) -> int:
        return self.op2.size
    @property
    def factor(self) -> int:
        return self.op2.factor

    def get_onr_prefix_postfix(self) -> tuple[str, str]:
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
        op2 = self.op2
        prefix = ''
        postfix = ''
        if op2.table_name in [b'ONRGY1', b'ONRGY2', b'ONRGY']:
            prefix = 'strain_energy.'
        elif op2.table_name in [b'RANEATC']: #, b'OSTRMS1C']:
            op2.format_code = 1
            op2.sort_bits[0] = 0 # real
            prefix = 'RANEATC.'
        elif op2.table_name in [b'RANCONS']: #, b'OSTRMS1C']:
            op2.format_code = 1
            op2.sort_bits[0] = 0 # real
            prefix = 'RANCONS.'
        else:
            raise NotImplementedError(op2.table_name)
        op2.data_code['sort_bits'] = op2.sort_bits
        op2.data_code['nonlinear_factor'] = op2.nonlinear_factor
        return prefix, postfix

    def _read_onr1_3(self, data: bytes, ndata: int):
        """
        reads ONRGY1 subtable 3
        """
        op2 = self.op2
        op2._analysis_code_fmt = b'i'
        op2.words = [
            'aCode',       'tCode',    'eTotal',       'isubcase',
            '???',         '???',      'element_name', 'load_set',
            'format_code', 'num_wide', 'cvalres',      'setID',
            'setID',       'eigenReal', 'eigenImag',   'rmssf',
            'etotpos',     'etotneg',  'thresh',       '???',
            '???',         '???',      '???',      '???',
            '???', 'Title', 'subtitle', 'label']

        #aCode = self.get_block_int_entry(data, 1)

        ## total energy of all elements in isubcase/mode
        self.etotal = op2.parse_approach_code(data)
        if op2.is_debug_file:
            op2.binary_debug.flush()

        self._onr_element_name(data)

        #: Load set or zero
        op2.load_set = op2.add_data_parameter(data, 'load_set', b'i', 8, False)

        #: format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        #: number of words per entry in record
        #: .. note:: is this needed for this table ???
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)
        ## C
        op2.cvalres = op2.add_data_parameter(data, 'cvalres', b'i', 11, False)

        #: Set identification number Number
        op2.set_id = op2.add_data_parameter(data, 'set_id', b'i', 13, False)

        #: Natural eigenvalue - real part
        op2.eigen_real = op2.add_data_parameter(data, 'eigen_real', b'f', 14, False)

        #: Natural eigenvalue - imaginary part
        op2.eigen_imag = op2.add_data_parameter(data, 'eigen_imag', b'f', 15, False)

        #: Natural frequency
        op2.freq = op2.add_data_parameter(data, 'freq', b'f', 16, False)

        #: RMS and CRMS scale factor - NX
        op2.rmssf = op2.add_data_parameter(data, 'rmssf', b'f', 17)

        #: Total positive energy
        op2.etotpos = op2.add_data_parameter(data, 'etotpos', b'f', 18)

        #: Total negative energy
        op2.etotneg = op2.add_data_parameter(data, 'etotneg', b'f', 19, False)

        #: Energy Threshold - NX
        op2.thresh = op2.add_data_parameter(data, 'thresh', b'f', 17)

        if not op2.is_sort1:
            raise NotImplementedError('sort2...')

        if op2.analysis_code == 1:   # statics / displacement / heat flux
            #del op2.data_code['nonlinear_factor']
            op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['lsdvmn'])
            op2.setNullNonlinearFactor()
        elif op2.analysis_code == 2:  # real eigenvalues
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)  ## mode number
            #op2.mode_cycle1 = op2.add_data_parameter(data, 'mode', b'i', 7)
            #op2.mode_cycle2 = op2.add_data_parameter(data, 'mode', b'f', 7)
            #print('mode = ', op2.mode)
            #print('mode_cycle1 = ', op2.mode_cycle1)
            #print('mode_cycle2 = ', op2.mode_cycle2)
            #self.show_data(data)
            #op2.cycle = 0.
            #self.reader_oug.update_mode_cycle('cycle')
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'freq'])
            #print("mode(5)=%s eign(6)=%s mode_cycle(7)=%s" % (
                #op2.mode, self.eign, op2.mode_cycle))
        #elif op2.analysis_code == 3: # differential stiffness
            #op2.lsdvmn = self.get_values(data,'i',5) ## load set number
            #op2.data_code['lsdvmn'] = op2.lsdvmn
        #elif op2.analysis_code == 4: # differential stiffness
            #op2.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif op2.analysis_code == 5:   # frequency
            op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)  ## frequency
            op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        elif op2.analysis_code == 6:  # transient
            op2.time = op2.add_data_parameter(data, 'time', b'f', 5)  ## time step
            op2.data_names = op2.apply_data_code_value('data_names', ['time'])
        #elif op2.analysis_code == 7: # pre-buckling
            #op2.data_names = op2.apply_data_code_value('data_names',['lsdvmn'])
        elif op2.analysis_code == 8:  # post-buckling
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)  ## mode number
            op2.data_names = op2.apply_data_code_value('data_names', ['mode'])
        elif op2.analysis_code == 9:  # complex eigenvalues
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)  ## mode number
            op2.eigr = op2.eigen_real
            op2.eigi = op2.eigen_imag
            op2.data_code['eigr'] = op2.eigr
            op2.data_code['eigi'] = op2.eigi
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif op2.analysis_code == 10:  # nonlinear statics
            self.loadFactor = op2.add_data_parameter(data, 'loadFactor', b'f', 5)  ## load factor
            op2.data_names = op2.apply_data_code_value('data_names', ['loadFactor'])
        #elif op2.analysis_code == 11: # old geometric nonlinear statics
            #op2.data_names = op2.apply_data_code_value('data_names',['lsdvmn'])
        elif op2.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            op2.time = op2.add_data_parameter(data, 'time', b'f', 5)  ## time step
            op2.data_names = op2.apply_data_code_value('data_names', ['time'])
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               op2.analysis_code)

        op2.fix_format_code()
        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()

    def _onr_element_name(self, data: bytes) -> str:
        op2 = self.op2
        #field_num = 6
        #datai = data[4 * (field_num - 1) : 4 * (field_num + 1)]
        #assert len(datai) == 8, len(datai)
        #print(4 * (field_num - 1), 4 * (field_num + 1))
        #element_name, = op2.struct_8s.unpack(data[24:32])  # changed on 11/30/2015; was this for a long time...

        #self.show_data(data[:28])
        if self.size == 4:
            element_name, = op2.struct_8s.unpack(data[20:28])
        else:
            element_name, = op2.struct_16s.unpack(data[40:56])
            element_name = reshape_bytes_block(element_name)

        #print("element_name = %s" % (element_name))
        try:
            element_name = element_name.decode('utf-8').strip()  # element name
        except UnicodeDecodeError:
            #self.log.warning("element_name = %s" % str(element_name))
            op2.log.warning("element_name - UnicodeDecodeError")
            #self.show_data(data)
            raise
        if element_name.isalnum():
            op2.data_code['element_name'] = element_name
        else:
            #print("element_name = %r" % (element_name))
            op2.data_code['element_name'] = 'UnicodeDecodeError???'
            op2.log.warning('data[20:28]=%r instead of data[24:32]' % data[20:28])

    def _read_onr2_3(self, data: bytes, ndata: int):
        """reads the SORT2 version of table 4 (the data table)"""
        op2 = self.op2
        op2.nonlinear_factor = np.nan
        op2.is_table_1 = False
        op2.is_table_2 = True
        unused_three = op2.parse_approach_code(data)
        op2.words = [
            'aCode',       'tCode',    'eTotal',       'isubcase',
            '???',         '???',      'element_name', 'load_set',
            'format_code', 'num_wide', 'cvalres',      'setID',
            'setID',       'eigenReal', 'eigenImag',   'rmssf',
            'etotpos',     'etotneg',  'thresh',       '???',
            '???',         '???',      '???',      '???',
            '???', 'Title', 'subtitle', 'label']

        self._onr_element_name(data)

        #: Load set or zero
        op2.load_set = op2.add_data_parameter(data, 'load_set', b'i', 8, False)

        #: format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        #: number of words per entry in record
        #: .. note:: is this needed for this table ???
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)
        ## C
        op2.cvalres = op2.add_data_parameter(data, 'cvalres', b'i', 11, False)

        #: Set identification number Number
        op2.set_id = op2.add_data_parameter(data, 'set_id', b'i', 13, False)

        #: Natural eigenvalue - real part
        #op2.eigen_real = op2.add_data_parameter(data, 'eigen_real', b'f', 14, False)

        #: Natural eigenvalue - imaginary part
        #op2.eigen_imag = op2.add_data_parameter(data, 'eigen_imag', b'f', 15, False)

        #: Natural frequency
        op2.freq = op2.add_data_parameter(data, 'freq', b'f', 16, False)

        #: RMS and CRMS scale factor - NX
        op2.rmssf = op2.add_data_parameter(data, 'rmssf', b'f', 17)

        #: Total positive energy
        op2.etotpos = op2.add_data_parameter(data, 'etotpos', b'f', 18)

        #: Total negative energy
        op2.etotneg = op2.add_data_parameter(data, 'etotneg', b'f', 19, False)

        #: Energy Threshold - NX
        op2.thresh = op2.add_data_parameter(data, 'thresh', b'f', 17)

        op2.element_id = op2.add_data_parameter(data, 'node_id', b'i', 5, fix_device_code=True)
        #if op2.analysis_code == 1:  # statics / displacement / heat flux
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            #op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            #op2.setNullNonlinearFactor()

        if op2.analysis_code == 1:  # static...because reasons.
            op2._analysis_code_fmt = b'i'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'N/A')
        elif op2.analysis_code == 2:  # real eigenvalues
            ## mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            op2._analysis_code_fmt = b'i'
            ## real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            op2.mode_cycle = op2.add_data_parameter(data, 'mode_cycle', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names',
                                                       ['node_id', 'eigr', 'mode_cycle'])
            op2.apply_data_code_value('analysis_method', 'mode')
        #elif op2.analysis_code == 3: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #op2.data_names = op2.data_code['lsdvmn'] = op2.lsdvmn
        #elif op2.analysis_code == 4: # differential stiffness
            #op2.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        elif op2.analysis_code == 5:   # frequency
            ## frequency
            #op2.freq = op2.add_data_parameter(data, 'freq', b'f', 5)
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'freq')
        elif op2.analysis_code == 6:  # transient
            ## time step
            #op2.dt = op2.add_data_parameter(data, 'dt', b'f', 5)
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'dt')
        elif op2.analysis_code == 7:  # pre-buckling
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2._analysis_code_fmt = b'i'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'lsdvmn')
        elif op2.analysis_code == 8:  # post-buckling
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2._analysis_code_fmt = b'i'
            ## real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id', 'eigr'])
            op2.apply_data_code_value('analysis_method', 'eigr')
        elif op2.analysis_code == 9:  # complex eigenvalues
            ## mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            op2._analysis_code_fmt = b'i'
            ## real eigenvalue
            #op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id', 'eigr', 'eigi'])
            op2.apply_data_code_value('analysis_method', 'mode')
        elif op2.analysis_code == 10:  # nonlinear statics
            ## load step
            #self.lftsfq = op2.add_data_parameter(data, 'lftsfq', b'f', 5)
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'lftsfq')
        elif op2.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2._analysis_code_fmt = b'f'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
        elif op2.analysis_code == 12:
            # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            #op2.lsdvmn = op2.add_data_parameter(data, 'lsdvmn', b'i', 5)
            op2._analysis_code_fmt = b'i'
            op2.data_names = op2.apply_data_code_value('data_names', ['node_id'])
            op2.apply_data_code_value('analysis_method', 'lsdvmn')
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % op2.analysis_code
            raise RuntimeError(msg)

        op2.fix_format_code()
        if op2.num_wide == 8:
            op2.format_code = 1
            op2.data_code['format_code'] = 1
        else:
            #op2.fix_format_code()
            if op2.format_code == 1:
                op2.format_code = 2
                op2.data_code['format_code'] = 2
            assert op2.format_code in [2, 3], op2.code_information()

        if op2.is_debug_file:
            op2.binary_debug.write('  approach_code  = %r\n' % op2.approach_code)
            op2.binary_debug.write('  tCode          = %r\n' % op2.tCode)
            op2.binary_debug.write('  isubcase       = %r\n' % op2.isubcase)
        op2._read_title(data)
        op2._write_debug_bits()

    def _read_onr1_4(self, data: bytes, ndata: int):
        """
        reads ONRGY1 subtable 4
        """
        op2 = self.op2
        if op2.table_code == 18:  # element strain energy
            if op2.table_name not in op2.table_name in [b'ONRGY', b'ONRGY1']:
                msg = f'table_name={op2.table_name} table_code={op2.table_code}'
                raise NotImplementedError(msg)
            n = self._read_element_strain_energy(data, ndata)
        else:
            raise NotImplementedError(op2.table_code)
        return n

    def _read_element_strain_energy(self, data: bytes, ndata: int):
        """
        table_code = 19
        """
        op2 = self.op2
        dt = op2.nonlinear_factor
        n = 0

        element_name = op2.data_code['element_name']
        try:
            result_name = RESULT_NAME_MAP[element_name]
        except KeyError:
            raise NotImplementedError('element_name1=%r element_name=%r' % (element_name, op2.data_code['element_name']))
        prefix, postfix = self.get_onr_prefix_postfix()
        result_name = prefix + result_name + postfix
        #result_name = 'strain_energy'

        if op2._results.is_not_saved(result_name):
            return ndata
        op2._results._found_result(result_name)

        slot = op2.get_result(result_name)

        #auto_return = False
        if op2.is_debug_file:
            op2.binary_debug.write('cvalares = %s\n' % op2.cvalres)
        if op2.format_code in [1, 2] and op2.num_wide == 4:
            assert op2.cvalres in [0, 1], op2.cvalres

            assert op2.num_wide == 4
            ntotal = 16 * self.factor  # 4*4=16
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, RealStrainEnergyArray)

            if auto_return:
                #if obj.dt_temp is None or obj.itime is None and obj.dt_temp == dt:
                    #element_name = op2.data_code['element_name']
                    #if element_name in obj.element_name_count:
                        #obj.element_name_count[element_name] += nelements
                    #else:
                        #obj.element_name_count[element_name] = nelements
                    #obj.dt_temp = dt
                return nelements * ntotal
            #itime = obj.itime #// obj.nelement_types
            #op2.show_data(data, types='if')

            obj = op2.obj
            itime = obj.itime

            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, energy, percent, density]\n')
                op2.binary_debug.write('  nelements=%i\n' % nelements)

            if op2.is_optistruct:
                op2.use_vector = False

            if op2.use_vector and op2.sort_method == 1: # and op2.is_sort1:
                n = nelements * ntotal
                ielement = obj.ielement
                ielement2 = obj.ielement + nelements
                itotal = obj.itotal
                itotal2 = obj.itotal + nelements * 4

                floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 4)
                obj._times[itime] = dt
                #if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 4)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, f'etype={element_name} isubtable={op2.isubtable} eids.min()={eids.min()}'
                obj.element[itime, ielement:ielement2] = eids

                #[energy, percent, density]
                obj.data[itime, ielement:ielement2, :] = floats[:, 1:].copy()
                obj.itotal2 = itotal2
                obj.ielement = ielement2
            else:
                n = real_strain_energy_4(op2, data, op2.sort_method,
                                         self.size, n, ntotal, nelements, dt)
        elif op2.format_code == 1 and op2.num_wide == 5:
            assert op2.cvalres in [0, 1, 2], op2.cvalres # 0??
            ntotal = 20
            nnodes = ndata // ntotal
            nelements = nnodes

            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, RealStrainEnergyArray)
            if auto_return:
                return nelements * op2.num_wide * 4

            obj = op2.obj
            if op2.use_vector:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 5).copy()
                obj._times[obj.itime] = dt

                strings = np.frombuffer(data, dtype=op2._uendian + 'S4').reshape(nelements, 5)
                if obj.itime == 0:
                    ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 5)
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
                n = complex_strain_energy_4(op2, data, op2.sort_method,
                                            self.size, n, ntotal, nelements, dt)
        elif op2.format_code in [2, 3] and op2.num_wide == 5:
            #ELEMENT-ID   STRAIN-ENERGY (MAG/PHASE)  PERCENT OF TOTAL  STRAIN-ENERGY-DENSITY
            #    5         2.027844E-10 /   0.0            1.2581            2.027844E-09
            ntotal = 20
            nelements = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, ComplexStrainEnergyArray)
            if auto_return:
                return nelements * op2.num_wide * 4

            obj = op2.obj
            if op2.use_vector:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 5)
                obj._times[obj.itime] = dt

                #if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 5)
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
                s = Struct(op2._endian + b'i4f')
                for unused_i in range(nelements):
                    edata = data[n:n+20]
                    out = s.unpack(edata)
                    (eid_device, energyr, energyi, percent, density) = out
                    eid = eid_device // 10
                    #if is_magnitude_phase:
                        #energy = polar_to_real_imag(energyr, energyi)
                    #else:
                        #energy = complex(energyr, energyi)

                    if op2.is_debug_file:
                        op2.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))
                    obj.add_sort1(dt, eid, energyr, energyi, percent, density)
                    n += ntotal

        elif op2.format_code == 1 and op2.num_wide == 6:  ## TODO: figure this out...
            ntotal = 24
            nnodes = ndata // ntotal
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, RealStrainEnergyArray)

            if auto_return:
                return nelements * op2.num_wide * 4

            obj = op2.obj
            if op2.use_vector:
                n = nelements * 4 * op2.num_wide
                itotal = obj.ielement
                ielement2 = obj.itotal + nelements
                itotal2 = ielement2

                floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 5)
                obj._times[obj.itime] = dt

                if obj.itime == 0:
                    strings = np.frombuffer(data, dtype=op2._uendian + 'S4').reshape(nelements, 6)
                    s = np.array([s1+s2 for s1, s2 in zip(strings[:, 1], strings[:, 2])])

                    ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 6)
                    eids = ints[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    obj.element[itotal:itotal2] = eids
                    obj.element_type[obj.itime, itotal:itotal2, :] = s

                #[energy, percent, density]
                obj.data[obj.itime, itotal:itotal2, :] = floats[:, 4:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2
            else:
                struct1 = Struct(op2._endian + b'i8s3f')
                for unused_i in range(nnodes):
                    edata = data[n:n+24]
                    out = struct1.unpack(edata)
                    (word, energy, percent, density) = out  # TODO: this has to be wrong...
                    word = word.strip()
                    #print "eType=%s" % (eType)
                    #print "%s" %(self.get_element_type(self.element_type)), data_in
                    #eid = op2.obj.add_new_eid_sort1(out)
                    if op2.is_debug_file:
                        op2.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))
                    obj.add_sort1(dt, word, energy, percent, density)
                    n += ntotal
        elif op2.format_code in [2, 3] and op2.num_wide == 4:
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
            msg = op2.code_information()
            return op2._not_implemented_or_skip(data, ndata, msg)
            #raise NotImplementedError(op2.code_information())
        return n

def real_strain_energy_4(op2: OP2,
                         data: bytes,
                         sort_method: int,
                         size: int,
                         n: int,
                         ntotal: int,
                         nelements: int,
                         dt) -> int:
    """
    (eid_device       eid energy percent  density)
    (11                 1 0.0114 0.1983   0.01147) typical
    ( 0                 0      0      0        -1) optistruct - final
    (1000000000 100000000    sum    sum       NaN) nx/msc     - final

    """
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'3f', size)
    struct1 = Struct(fmt)
    obj = op2.obj  # type: RealStrainEnergyArray

    if op2.is_optistruct:
        fmt2 = mapfmt(op2._endian + op2._analysis_code_fmt + b'2f i', size)
        struct2 = Struct(fmt2)

        edata = data[n:n+ntotal]
        sum_energy = 0.
        sum_percent = 0.
        for unused_i in range(nelements-1):
            edata = data[n:n+ntotal]

            out = struct1.unpack(edata)
            (eid_device, energy, percent, density) = out
            if sort_method == 1:
                eid = eid_device // 10
            else:
                eid = op2.nonlinear_factor
                dt = eid_device
            #print(f'adding dt={dt:g} eid_device={eid_device} eid={eid} energy={energy:g} percent={percent:g} density={density:g}')
            if op2.is_debug_file:
                op2.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))
            sum_energy += energy
            sum_percent += percent
            obj.add_sort1(dt, eid, energy, percent, density)
            n += ntotal
            edata = data[n:n+ntotal]

        out = struct2.unpack(edata)
        (eid_device, energy, percent, density) = out
        assert eid_device == 0, eid_device
        if sort_method == 1:
            eid = eid_device // 10
        else:
            raise NotImplementedError(sort_method)
            #eid = op2.nonlinear_factor
            #dt = eid_device
        #print(f'adding dt={dt:g} eid_device={eid_device} eid={eid} energy={energy:g} percent={percent:g} density={density:g}')
        #if op2.is_debug_file:
            #op2.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))

        eid = 100000000
        obj.add_sort1(dt, eid, sum_energy, sum_percent, np.nan)
        n += ntotal
    else:
        for unused_i in range(nelements):
            edata = data[n:n+ntotal]
            out = struct1.unpack(edata)
            (eid_device, energy, percent, density) = out
            if sort_method == 1:
                eid = eid_device // 10
            else:
                eid = op2.nonlinear_factor
                dt = eid_device
            #print(f'adding dt={dt:g} eid_device={eid_device} eid={eid} energy={energy} percent={percent:g} density={density:g}')
            if op2.is_debug_file:
                op2.binary_debug.write('  eid=%i; %s\n' % (eid, str(out)))
            obj.add_sort1(dt, eid, energy, percent, density)
            n += ntotal
    return n

def complex_strain_energy_4(op2, data, sort_method,
                            size, n, ntotal, nnodes, dt):
    obj = op2.obj  # type: ComplexStrainEnergyArray
    s = Struct(op2._endian + b'8s3f')
    for unused_i in range(nnodes):
        edata = data[n:n+20]
        out = s.unpack(edata)
        (word, energy, percent, density) = out
        word = word.strip()
        if op2.is_debug_file:
            op2.binary_debug.write('  eid/word=%r; %s\n' % (word, str(out)))
        obj.add_sort1(dt, word, energy, percent, density)
        n += ntotal
    return n
