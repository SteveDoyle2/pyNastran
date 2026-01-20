from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

class OEFPK:
    """Defines OEFx table reading for element forces/heat flux"""
    def __init__(self, op2: OP2):
        self.op2 = op2

    def _create_oes_object4(self, *args, **kwargs):
        return self.op2._create_oes_object4(*args, **kwargs)

    #def get_oef_prefix_postfix(self):
        #"""
        #NX Case Control  Block         Description
        #===============  ==========    ===========
        #NLSTRESS         OESNLXR       Nonlinear static stresses
        #BOUTPUT          OESNLBR       Slideline stresses
        #STRESS           OESNLXD       Nonlinear Transient Stresses
        #STRESS           OES1C/OSTR1C  Ply stresses/strains
        #STRESS           OES1X         Element stresses with intermediate (CBAR and CBEAM)
                                       #station stresses and stresses on nonlinear elements
        #STRESS           OES/OESVM     Element stresses (linear elements only)
        #STRAIN           OSTR1         Element strains
        #STRESS/STRAIN    DOES1/DOSTR1  Scaled Response Spectra
        #MODCON           OSTRMC        Modal contributions
        #"""
        #op2 = self.op2
        #prefix = ''
        #postfix = ''
        #table_name_bytes = op2.table_name
        #assert isinstance(table_name_bytes, bytes), table_name_bytes
        #is_sort1 = table_name_bytes in SORT1_TABLES_BYTES
        #assert table_name_bytes in TABLES_BYTES, table_name_bytes

        #if table_name_bytes in [b'OEF1X', b'OEF1', b'OEF2']:
            #if op2.thermal == 0:
                #prefix = 'force.'
            #elif op2.thermal == 1:
                #prefix = 'thermal_load.'
            #else:
                #raise NotImplementedError(op2.code_information())
        #elif table_name_bytes in [b'HOEF1']:
            #postfix = '_flux'
        ##elif op2.table_name in ['OESNLXR']:
            ##prefix = 'sideline_'
        ##elif op2.table_name in ['OESNLXD', 'OESNL1X', 'OESNLXR']:
            ##prefix = 'nonlinear_'
        ##elif op2.table_name == 'OESNLBR':
            ##prefix = 'sideline_'
        ##elif op2.table_name == 'OESRT':
            ##prefix = 'strength_ratio.'
        ##elif op2.table_name in ['OESCP', 'OESTRCP']:
            ##pass # TODO: update
        #elif table_name_bytes in [b'OEFCRM1', b'OEFCRM2']:
            #assert op2.table_code in [4, 504], op2.code_information()
            #prefix = 'crm.'
            #op2.reader_oes._set_as_random()
        #elif table_name_bytes in [b'OEFPSD1', b'OEFPSD2']:
            #assert op2.table_code in [4, 604], op2.code_information()
            #op2.reader_oes._set_as_random()
            #prefix = 'psd.'
        #elif table_name_bytes in [b'OEFRMS1', b'OEFRMS2']:
            #assert op2.table_code in [4, 804], op2.code_information()
            #op2.reader_oes._set_as_random()
            #is_sort1 = True
            #op2._analysis_code_fmt = b'i'
            #prefix = 'rms.'
        #elif table_name_bytes in [b'OEFNO1', b'OEFNO2']:
            #assert op2.table_code in [4, 904], op2.code_information()
            #op2.reader_oes._set_as_random()
            #op2.sort_method = 1
            #op2.data_code['nonlinear_factor'] = None
            #op2._analysis_code_fmt = b'i'
            #assert op2.sort_method == 1, op2.code_information()
            #prefix = 'no.'
        #elif table_name_bytes in [b'OEFATO1', b'OEFATO2']:
            #assert op2.table_code in [4], op2.code_information()
            #prefix = 'ato.'

        #elif table_name_bytes in [b'RAFCONS']:
            #prefix = 'RAFCONS.'
        #elif table_name_bytes in [b'RAFEATC']:
            #prefix = 'RAFEATC.'
        #elif table_name_bytes in [b'DOEF1']:
            #assert op2.table_code in [4], op2.code_information()
            #prefix = shock_response_prefix(op2.thermal)
        #elif table_name_bytes in [b'OEFIT']:
            #assert op2.table_code in [25], op2.code_information()
            #prefix = 'failure_indices.'
            ##raise NotImplementedError(op2.code_information())
        #elif table_name_bytes in [b'OEFITSTN']: # composite failure indicies
            #assert op2.table_code in [25], op2.code_information()
            #prefix = 'failure_indices.'
        #else:
            #raise NotImplementedError('%r' % op2.table_name)

        #op2.sort_bits.is_sort1 = is_sort1 # sort2
        #op2.sort_method = 1 if is_sort1 else 2
        #return prefix, postfix

    def _oef_force_code(self):
        """
        Gets the numwide codes for the element to determine if
        the real or complex result should be found.
        The format and sort codes do not always give the right answer...
        """
        op2 = self.op2
        if op2.is_nx:
            real_mapper = NX_OEF_REAL_MAPPER
            imag_mapper = NX_OEF_IMAG_MAPPER
        else:
            real_mapper = MSC_OEF_REAL_MAPPER
            imag_mapper = MSC_OEF_IMAG_MAPPER


        try:
            real = real_mapper[op2.element_type]
        except KeyError:
            real = None

        try:
            imag = imag_mapper[op2.element_type]
        except KeyError:
            imag = None
        return real, imag

    def _read_oefpk1_3(self, data: bytes, ndata: int):
        """Table 3 parser for OEF1 table"""
        op2 = self.op2
        op2._analysis_code_fmt = b'i'
        op2._data_factor = 1
        op2.words = [
            'aCode', 'tCode', 'element_type', 'isubcase',
            '???', '???', '???', '???',
            'format_code', 'num_wide', 'o_code', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', 'Title', 'subtitle', 'label']

        op2.parse_approach_code(data)

        #: element type
        op2.element_type = op2.add_data_parameter(data, 'element_type', b'i', 3, False)

        # dynamic load set ID/random code
        #self.dLoadID = op2.add_data_parameter(data, 'dLoadID', b'i', 8, False)

        #: format code
        op2.format_code = op2.add_data_parameter(data, 'format_code', b'i', 9, False)

        #: number of words per entry in record
        #: .. note: is this needed for this table ???
        op2.num_wide = op2.add_data_parameter(data, 'num_wide', b'i', 10, False)

        #: undefined in DMAP...
        op2.o_code = op2.add_data_parameter(data, 'o_code', b'i', 11, False)

        #: thermal flag; 1 for heat ransfer, 0 otherwise
        op2.thermal = op2.add_data_parameter(data, 'thermal', b'i', 23, False)

        ## assuming tCode=1
        if op2.analysis_code == 1:   # statics
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5, False)  # load set ID number
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
            op2.setNullNonlinearFactor()
        elif op2.analysis_code == 2:  # normal modes/buckling (real eigenvalues)
            #: mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            #: eigenvalue
            op2.eign = op2.add_data_parameter(data, 'eign', b'f', 6, False)
            op2.cycle = 0.
            op2._op2_readers.reader_oug.update_mode_cycle('cycle')
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'cycle'])
            # TODO: mode_cycle is not defined?
            #op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])
        elif op2.analysis_code == 3:  # differential stiffness 0
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        elif op2.analysis_code == 4:  # differential stiffness 1
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        elif op2.analysis_code == 5:   # frequency
            self.freq = op2.add_data_parameter(data, 'freq', b'f', 5)  # frequency
            op2.data_names = op2.apply_data_code_value('data_names', ['freq'])
        elif op2.analysis_code == 6:  # transient
            self.time = op2.add_data_parameter(data, 'time', b'f', 5)  # time step
            op2.data_names = op2.apply_data_code_value('data_names', ['time'])
        elif op2.analysis_code == 7:  # pre-buckling
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            #op2.apply_data_code_value('data_names',['lsdvmn'])
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        elif op2.analysis_code == 8:  # post-buckling
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            #: real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID', 'eigr'])
        elif op2.analysis_code == 9:  # complex eigenvalues
            #: mode number
            op2.mode = op2.add_data_parameter(data, 'mode', b'i', 5)
            #: real eigenvalue
            op2.eigr = op2.add_data_parameter(data, 'eigr', b'f', 6, False)
            #: imaginary eigenvalue
            op2.eigi = op2.add_data_parameter(data, 'eigi', b'f', 7, False)
            op2.data_names = op2.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        elif op2.analysis_code == 10:  # nonlinear statics
            #: load step
            self.load_step = op2.add_data_parameter(data, 'load_step', b'f', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['load_step'])
        elif op2.analysis_code == 11:  # geometric nonlinear statics
            #: load set ID number
            op2.loadID = op2.add_data_parameter(data, 'loadID', b'i', 5)
            op2.data_names = op2.apply_data_code_value('data_names', ['loadID'])
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % str(op2.analysis_code)
            raise RuntimeError(msg)

        op2.fix_format_code()
        op2._parse_thermal_code()
        op2._op2_readers.reader_oef._set_force_stress_element_name()

        if op2.is_debug_file:
            op2.binary_debug.write('  %-14s = %r\n' % ('element_name', op2.element_name))
            op2.binary_debug.write('  %-14s = %r %s\n' % ('approach_code', op2.approach_code,
                                                           op2.approach_code_str(op2.approach_code)))
            op2.binary_debug.write('  %-14s = %r\n' % ('tCode', op2.tCode))
            op2.binary_debug.write('  %-14s = %r\n' % ('isubcase', op2.isubcase))


        op2._read_title(data)
        if op2.element_type not in op2.element_mapper:
            msg = 'element_type = %s' % op2.element_type
            return op2._not_implemented_or_skip(data, ndata, msg)
        op2._write_debug_bits()
        assert op2.num_wide != 146, op2.code_information()
        #print('OEF-%s' % op2.element_name)
        #self._check_result_type()

    def _read_oefpk1_4(self, data: bytes, ndata: int):
        if self.op2.table_code == 404:
            pass
        else:
            raise RuntimeError(self.op2.code_information())

        n = self.op2._op2_readers.reader_oef._read_oef1_loads(data, ndata)
        return n
