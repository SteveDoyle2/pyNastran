# pylint: disable=C0301,W0201
import copy
from struct import Struct, unpack
from typing import Tuple, Dict, Union, Any

import numpy as np

from pyNastran import is_release
from pyNastran.op2.errors import OverwriteTableError
from pyNastran.f06.f06_writer import F06Writer
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_interface.utils import (
    build_obj, get_superelement_adaptivity_index, update_subtitle_with_adaptivity_index,
    update_label2)
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block
from pyNastran.op2.op2_interface.op2_codes import (
    Op2Codes, get_scode_word, get_sort_method_from_table_name)

from pyNastran.op2.errors import SortCodeError, MultipleSolutionNotImplementedError

NX_TABLES = [
    501, 510, 511,
    601, 610, 611,
    701, 710, 711,
    801, 810, 811,
    901, 910, 911
]

class OP2Common(Op2Codes, F06Writer):
    def __init__(self):
        Op2Codes.__init__(self)
        F06Writer.__init__(self)

        #: flag for vectorization
        #: 0 - no vectorization
        #: 1 -   first pass
        #: 2 -   second pass
        self.read_mode = None

        # Cross valdation flag so we can write:
        #   >>> modelA = OP2()
        #   >>> modelA.read_op2(op2_filename)
        #   >>> modelB = OP2()
        #   >>> model.B.use_vector = False
        #   >>> modelB.read_op2(op2_filename)
        #   >>> assert model A == modelB
        # vectorized:
        #    False : uses range loops (testing)
        #    True : uses vectorization
        # non-vectorized:
        #    does nothing
        # -----------------
        self.use_vector = True

        # is a debug file being written to
        self.is_debug_file = False

        # how many optimization passes have there been
        self._count = 0

        #: the results
        self.result_names = set()
        #: bool
        self.is_vectorized = None
        self.combine = False

        #: the storage dictionary that is passed to OP2 objects (e.g. RealDisplacementArray)
        #: the key-value pairs are extracted and used to generate dynamic self
        #: variables for the OP2 objects
        #self.data_code = {
            #'_encoding' : self._encoding,
            #'load_as_h5' :
        #}

        #: current subcase ID
        #: non-transient (SOL101) cases have isubcase set to None
        #: transient (or frequency/modal) cases have isubcase set to a int/float value
        self.isubcase = None

        #: the corresponding piece to isubcase
        #: used only for SORT2 (not supported)
        self.ID = None

        #: should the op2 debugging file be written
        self.debug = False

        #: op2 debug file or None (for self.debug=False)
        self.binary_debug = None

        #: the list of "words" on a subtable 3
        self.words = []  # List[Optional[str]]

        #: The current table_name (e.g. OES1)
        #: None indicates no table_name has been read
        self.table_name = None

        # the date stamp used in the F06
        self.date = (1, 1, 2000)

        #: set of all the subcases that have been found
        #self.subcases = set()

        #: the list/set/tuple of times/modes/frequencies that should be read
        #: currently unused
        self.expected_times = None

        self._endian = None

        #self.show_table3_map = [
            ##'OUGV1',
            ##'OEF1X',
            ##'OES1X1',
        #]
        #self.show_table4_map = [
            ##'OUGV1',
            ##'OEF1X',
            ##'OES1X1',
        #]

        # sets the element mapper
        #self.get_element_type(33)

    def _device_code_(self) -> None:
        """
        0 -
        1 - PRINT
        2 - PLOT
        3 - PRINT, PLOT
        4 - PUNCH
        5 - PRINT, PUNCH
        6 - PLOT, PUNCH
        7 - PRINT, PLOT, PUNCH

        """
        pass

    def fix_format_code(self) -> None:
        """
        Nastran correctly calculates the proper defaults for the analysis
        based on the solution type and the the user's requests.  However,
        the user doesn't always set the values correctly.  When Nastran
        goes to write the output, it uses the original values, rather
        than the correct values that were used for analysis.

        In a SOL 101 case:
            STRESS(PLOT, SORT1, IMAG) = ALL

        the user has set an incorrect value (IMAG), which gets turned into
        a format code of 2, where:
          1 - real
          2 - real/imag (complex results)
          3 - mag/phase (complex results)

        This inconsistency causes problems with the parser.  Thus, based on
        the analysis_code (1 is like SOL 101, but really means static), we
        can switch mag/phase results to real static results.

        .. note:: A case of 4 is not used and is used below as a placeholder,
                  while a case of -1 is some bizarre unhandled,
                  undocumented case.

        """
        self._set_times_dtype()
        self.format_code_original = self.format_code
        if self.format_code == -1:
            if self.is_debug_file:
                self.op2_reader._write_ndata(self.binary_debug, 100)
            if self.table_name in [b'OESNLXR', b'OESNLBR', b'OESNLXD', b'OESNL1X', b'OESNLBR2']:
                assert self.format_code == -1, self.format_code
                self.format_code = 1
            else:
                raise RuntimeError(self.code_information())
            #return

        random_code = self.random_code if hasattr(self, 'random_code') else 0
        #if self.format_code == 0:
            #self.code_information()
        if random_code == 0:
            #self.log.debug(f'random_code = {random_code}')
            if self.analysis_code == 1:   # statics / displacement / heat flux
                assert self.format_code in [1, 3], self.code_information()
                self.format_code = 1
            elif self.analysis_code == 2:  # real eigenvalues
                assert self.format_code in [1, 3], self.code_information()
                self.format_code = 1
            #elif self.analysis_code==3: # differential stiffness
            #elif self.analysis_code==4: # differential stiffness
            elif self.analysis_code == 5:   # frequency
                assert self.format_code in [1, 2, 3], self.code_information()
                if self.format_code == 1:
                    self.format_code = 2
            elif self.analysis_code == 6:  # transient
                self.format_code = 1
                assert self.format_code in [1, 2, 3], self.code_information()
            elif self.analysis_code == 7:  # pre-buckling
                assert self.format_code in [1], self.code_information()
            elif self.analysis_code == 8:  # post-buckling
                assert self.format_code in [1, 2], self.code_information()
            elif self.analysis_code == 9:  # complex eigenvalues
                assert self.format_code in [1, 2, 3], self.code_information()
                if self.format_code == 1:
                    self.format_code = 2
            elif self.analysis_code == 10:  # nonlinear statics
                assert self.format_code in [1], self.code_information()
            elif self.analysis_code == 11:  # old geometric nonlinear statics
                assert self.format_code in [1], self.code_information()
            elif self.analysis_code == 12:
                # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
                assert self.format_code in [4], self.code_information() # invalid value
            else:
                msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
                raise RuntimeError(msg)
            self.data_code['format_code'] = self.format_code
        #assert self.format_code == 1, self.code_information()
        #if self.format_code != self.format_code_original:
            #print('self.format_code=%s orig=%s' % (self.format_code,
                                                   #self.format_code_original))

    def _set_times_dtype(self) -> None:
        self.data_code['_times_dtype'] = 'float32'
        if self.analysis_code == 1:   # statics / displacement / heat flux
            pass # static doesn't have a type
        elif self.analysis_code == 2:  # real eigenvalues
            pass
        #elif self.analysis_code==3: # differential stiffness
        #elif self.analysis_code==4: # differential stiffness
        elif self.analysis_code == 5:  # frequency
            pass
        elif self.analysis_code == 6:  # transient
            pass
        elif self.analysis_code == 7:  # pre-buckling
            pass
        elif self.analysis_code == 8:  # post-buckling
            pass
        elif self.analysis_code == 9:  # complex eigenvalues
            pass
        elif self.analysis_code == 10:  # nonlinear statics
            pass
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            pass
        elif self.analysis_code == 12:
            # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            pass
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % self.analysis_code
            raise RuntimeError(msg)

    def add_data_parameter(self, data: bytes, var_name: str,
                           Type: bytes, field_num: int,
                           apply_nonlinear_factor: bool=True,
                           fix_device_code: bool=False,
                           add_to_dict: bool=True):
        if self.size == 4:
            return self.add_data_parameter4(
                data, var_name, Type, field_num,
                apply_nonlinear_factor=apply_nonlinear_factor,
                fix_device_code=fix_device_code,
                add_to_dict=add_to_dict)
        return self.add_data_parameter8(
            data, var_name, Type, field_num,
            apply_nonlinear_factor=apply_nonlinear_factor,
            fix_device_code=fix_device_code,
            add_to_dict=add_to_dict)

    def add_data_parameter8(self, data: bytes,
                            var_name: str, var_type: bytes, field_num: int,
                            apply_nonlinear_factor: bool=True,
                            fix_device_code: bool=False,
                            add_to_dict: bool=True):
        assert len(data) == 1168, len(data)
        datai = data[8 * (field_num - 1) : 8 * (field_num)]
        assert len(datai) == 8, len(datai)
        if var_type == b'i':
            var_type = b'q'
        elif var_type == b'f':
            var_type = b'd'
        else:  # pragma: no cover
            raise NotImplementedError(var_type)
        value = self._unpack_data_parameter(
            datai, var_name, var_type, field_num,
            apply_nonlinear_factor=apply_nonlinear_factor,
            fix_device_code=fix_device_code,
            add_to_dict=add_to_dict)
        return value

    def add_data_parameter4(self, data: bytes, var_name: str,
                            Type: bytes, field_num: int,
                            apply_nonlinear_factor: bool=True,
                            fix_device_code: bool=False,
                            add_to_dict: bool=True):
        assert len(data) == 584, len(data)
        datai = data[4 * (field_num - 1) : 4 * (field_num)]
        assert len(datai) == 4, len(datai)
        #assert type(self._endian) == type(Type), 'endian=%r Type=%r' % (self._endian, Type)
        value = self._unpack_data_parameter(
            datai, var_name, Type, field_num,
            apply_nonlinear_factor=apply_nonlinear_factor,
            fix_device_code=fix_device_code,
            add_to_dict=add_to_dict)
        return value

    def _unpack_data_parameter(self, datai: bytes, var_name: str,
                               Type: bytes, field_num: int,
                               apply_nonlinear_factor=True,
                               fix_device_code: bool=False,
                               add_to_dict: bool=True):
        value, = unpack(self._endian + Type, datai)
        if fix_device_code:
            # was value = (value - self.device_code) // 10
            value = value // 10
        if self.is_debug_file:
            self.binary_debug.write('  %-14s = %r\n' % (var_name, value))
        #setattr(self, var_name, value)  # set the parameter to the local namespace

        if apply_nonlinear_factor:
            npvalue = _cast_nonlinear_factor(value)
            self.nonlinear_factor = npvalue
            #if self.table_name == b'OUGV2':
                #assert isinstance(self.nonlinear_factor, int), self.nonlinear_factor
            self.data_code['nonlinear_factor'] = npvalue
            self.data_code['name'] = var_name

        if add_to_dict:
            self.data_code[var_name] = value

        try:
            self.words[field_num - 1] = var_name
        except IndexError:
            msg = 'Trying to set word, but...len(words)=%s ifield=%s' % (len(self.words), field_num - 1)
            raise IndexError(msg)
        return value

    def apply_data_code_value(self, name: str, value: Union[int, float, str]) -> None:
        self.data_code[name] = value

    def setNullNonlinearFactor(self) -> None:
        """
        Initializes the nonlinear factor, which lets us know if
        this is a transient solution or not.

        """
        self.nonlinear_factor = np.nan #np.float32(None)
        self.data_code['nonlinear_factor'] = np.nan

    def _read_title_helper(self, data):
        if self.size == 4:
            assert len(data) == 584, len(data)
            # title_subtitle_label
            title, subtitle, label = unpack(self._endian + b'128s128s128s', data[200:])
        else:
            assert len(data) == 1168, len(data)
            # title_subtitle_label
            title, subtitle, label = unpack(self._endian + b'256s256s256s', data[400:])
            title = reshape_bytes_block(title)
            subtitle = reshape_bytes_block(subtitle)
            label = reshape_bytes_block(label)

        self.title = title.decode(self.encoding).strip()
        subtitle = subtitle.decode(self.encoding)
        subtitle_original = subtitle.strip()

        label = label.decode(self.encoding).strip()
        nlabel = 65
        label2 = label[nlabel:]
        try:
            label2 = update_label2(label2, self.isubcase)
        except AssertionError:
            pass

        assert len(label[:nlabel]) <= nlabel, f'len={len(label)} \nlabel     ={label!r} \nlabel[:{nlabel}]={label[:nlabel]!r}'
        assert len(label2) <= 55, f'len={len(label2)} label = {label!r}\nlabel[:{nlabel}]={label[:nlabel]!r}\nlabel2    ={label2!r}'
        # not done...
        # 65 + 55 = 120 < 128

        self.label = label
        self.pval_step = label2

        nsubtitle_break = 67
        adpativity_index = subtitle[nsubtitle_break:99]
        superelement = subtitle[99:].strip()

        #print('subtitle = %r' % subtitle)
        #print('aindex   = %r' % adpativity_index)
        #print('superele = %r' % superelement)

        subtitle = subtitle[:nsubtitle_break].strip()
        assert len(superelement) <= 26, f'len={len(superelement)} superelement={superelement!r}'
        superelement = superelement.strip()

        assert len(subtitle) <= 67, f'len={len(subtitle)} subtitle={subtitle!r}'
        superelement_adaptivity_index = get_superelement_adaptivity_index(subtitle, superelement)
        subtitle, superelement_adaptivity_index = update_subtitle_with_adaptivity_index(
            subtitle, superelement_adaptivity_index, adpativity_index)
        self.subtitle = subtitle
        self.superelement_adaptivity_index = superelement_adaptivity_index
        assert len(self.label) <= 124, f'len={len(self.label)} label={self.label!r}'

        #: the subtitle of the subcase
        self.data_code['subtitle'] = self.subtitle
        self.data_code['subtitle_original'] = subtitle_original

        # the sub-key
        self.data_code['pval_step'] = self.pval_step
        self.data_code['superelement_adaptivity_index'] = self.superelement_adaptivity_index

        #: the label of the subcase
        self.data_code['label'] = self.label
        self.data_code['title'] = self.title

        if self.is_debug_file:
            self.binary_debug.write(
                '  %-14s = %r\n' * 6 % (
                    'count', self._count,
                    'title', self.title,
                    'subtitle', self.subtitle,
                    'label', self.label,
                    'pval_step', self.pval_step,
                    'superelement_adaptivity_index', self.superelement_adaptivity_index))

    def _read_title(self, data):
        self._read_title_helper(data)

        if hasattr(self, 'isubcase'):
            if self.isubcase not in self.isubcase_name_map:
                # 100 from label
                # 20 from subtitle line
                # 'SUBCASE 2'
                #self.isubcase_name_map[isubcase] = [self.Subtitle, self.label]
                self.isubcase_name_map[self.isubcase] = [
                    self.subtitle, self.superelement_adaptivity_index,
                    self.analysis_code, self.label]
        else:
            raise  RuntimeError('isubcase is not defined')

    def _write_debug_bits(self):
        """
        s_code =  0 -> stress_bits = [0,0,0,0,0]
        s_code =  1 -> stress_bits = [0,0,0,0,1]
        s_code =  2 -> stress_bits = [0,0,0,1,0]
        s_code =  3 -> stress_bits = [0,0,0,1,1]
        etc.
        s_code = 32 -> stress_bits = [1,1,1,1,1]

        stress_bits[0] = 0 -> isMaxShear=True       isVonMises=False
        stress_bits[0] = 1 -> isMaxShear=False      isVonMises=True

        stress_bits[1] = 0 -> is_stress=True        is_strain=False
        stress_bits[2] = 0 -> isFiberCurvature=True isFiberDistance=False
        stress_bits[3] = 0 -> duplicate of Bit[1] (stress/strain)
        stress_bits[4] = 0 -> material coordinate system flag

        """
        if self.is_debug_file:
            msg = ''
            assert len(self.words) in [0, 28], f'table_name={self.table_name!r} len(self.words)={len(self.words)} words={self.words}'
            for i, param in enumerate(self.words):
                if param == 's_code':
                    try:
                        s_word = get_scode_word(self.s_code, self.stress_bits)
                    except AttributeError:
                        raise
                    self.binary_debug.write('  s_code         = %s -> %s\n' % (self.s_code, s_word))
                    self.binary_debug.write('    stress_bits[0] = %i -> is_von_mises    =%-5s vs is_max_shear\n' % (self.stress_bits[0], self.is_von_mises))
                    self.binary_debug.write('    stress_bits[1] = %i -> is_strain       =%-5s vs is_stress\n' % (self.stress_bits[1], self.is_strain))
                    self.binary_debug.write('    stress_bits[2] = %i -> strain_curvature=%-5s vs fiber_dist\n' % (self.stress_bits[2], self.is_curvature))
                    self.binary_debug.write('    stress_bits[3] = %i -> is_strain       =%-5s vs is_stress\n' % (self.stress_bits[3], self.is_strain))
                    self.binary_debug.write('    stress_bits[4] = %i -> material coordinate system flag=%s vs ???\n' % (self.stress_bits[4], self.stress_bits[4]))
                elif param == '???':
                    param = 0
                msg += '%s, ' % param
                if i % 5 == 4:
                    msg += '\n             '

            if hasattr(self, 'format_code'):
                try:
                    is_complex = self.is_complex
                except AssertionError:
                    self.binary_debug.write('\n  ERROR: cannot determine is_complex properly; '
                                            'check_sort_bits!!!\n')
                    is_complex = '???'

                try:
                    is_random = self.is_random
                except AssertionError:
                    is_random = '???'

                try:
                    is_sort1 = self.is_sort1
                except AssertionError:
                    is_sort1 = '???'

                try:
                    is_real = self.is_real
                except AssertionError:
                    is_real = '???'

                if is_complex:
                    msg = '\n  %-14s = %i -> is_mag_phase vs is_real_imag vs. is_random\n' % (
                        'format_code', self.format_code)
                    self.binary_debug.write(msg)
                else:
                    self.binary_debug.write('  %-14s = %i\n' % ('format_code', self.format_code))
                self.binary_debug.write('    sort_bits[0] = %i -> is_random=%s vs mag/phase\n' % (self.sort_bits[0], is_random))
                self.binary_debug.write('    sort_bits[1] = %i -> is_sort1 =%s vs sort2\n' % (self.sort_bits[1], is_sort1))
                self.binary_debug.write('    sort_bits[2] = %i -> is_real  =%s vs real/imag\n' % (self.sort_bits[2], is_real))

                try:
                    sort_method, is_real, is_random = self._table_specs()
                    self.binary_debug.write('    sort_method = %s\n' % sort_method)
                except AssertionError:
                    self.binary_debug.write('    sort_method = ???\n')

                if is_complex:
                    msg = '\n  %-14s = %i -> is_mag_phase vs is_real_imag vs. is_random\n' % (
                        'format_code', self.format_code)
                    self.binary_debug.write(msg)
                else:
                    self.binary_debug.write('  %-14s = %i\n' % ('format_code', self.format_code))

            self.binary_debug.write('  recordi = [%s]\n\n' % msg)

    def get_table_count(self) -> int:
        """identifiers superelements"""
        #{#b'PVT0': 1, b'CASECC': 1,
        #b'GPLS': 2, b'GPDTS': 2, b'EPTS': 2, b'MPTS': 2,
        #b'GEOM1S': 1, b'GEOM2S': 2, b'GEOM3S': 1,
        #b'BGPDTS': 1, b'EQEXINS': 1,}
        keys = [
            #b'PVT0':, b'CASECC', b'BOUGV1',
            b'GPDTS', b'BGPDTS', b'GPLS', b'EQEXINS',
            b'GEOM1S', b'GEOM2S', b'GEOM3S', b'GEOM4S',
            b'EPTS', b'MPTS', 'DITS',
        ]
        #print('self.table_count =', self.table_count)
        max_geom_id = max([self.table_count[key] for key in keys])
        return max_geom_id

    def _read_geom_4(self, mapper: Dict[Tuple[int, int, int], Any],
                     data: bytes, ndata: int) -> int:
        """
        Reads a geometry table

        TODO: Callable[[bytes, int]]
        """
        if self.read_mode == 1:
            return ndata
        if not self.make_geom:
            return ndata

        max_geom_id = self.get_table_count()
        if max_geom_id > 1:
            self.log.warning('superelement 2 not supported.  This may crash.')

        n = 0
        if self.size == 4:
            ngeom = 12
            structi = self.struct_3i
        else:
            ngeom = 24
            structi = self.struct_3q

        keys = structi.unpack(data[n:n+ngeom])
        n += ngeom
        if len(data) == ngeom:
            #print('*self.istream = %s' % self.istream)
            #print('self.isubtable = %s' % self.isubtable)
            #self.istream -= 1 ## TODO: removed because it doesn't exist???
            self.isubtable_old = self.isubtable
            return n

        #print('is_start_of_subtable=%s' % self.is_start_of_subtable)
        #print('self.istream = %s' % self.istream)
        #if hasattr(self, 'isubtable_old'):
            #print("self.isubtable_old=%r self.isubtable=%s" % (self.isubtable_old, self.isubtable))
        #else:
            #print("self.isubtable=%s" % (self.isubtable))
        if not hasattr(self, 'isubtable_old'):
            self.isubtable_old = 1
        elif self.isubtable_old > self.isubtable:
            self.isubtable_old = 1

        #self.binary_debug.write('isubtable=%s isubtable_old=%s\n' % (self.isubtable, self.isubtable_old))
        #ni = self.f.tell() - len(data) + 12
        #self.binary_debug.write('**:  f.tell()=%s; n=%s:%s\n\n' % (self.f.tell(), ni, self.n))

        #if 0:
            ## we're only going to use the keys if istream=0 (so the beginning of the record)
            #if self.istream == 0 and keys in mapper:
                #pass
            #elif self.isubtable_old == self.isubtable:
                ## we didn't increment the record, so we fix the n+=12 statement we called before
                ## then we toss the keys and use the old geom_keys
                #n = 0
                #keys = self.geom_keys
            #else:
                #msg = 'keys=%s not found - %s; istream=%s; isubtable=%s isubtable_old=%s\n mapper=%s' % (
                    #str(keys), self.table_name, self.istream, self.isubtable, self.isubtable_old,
                    #mapper.keys())
                #raise NotImplementedError(msg)

        try:
            name, func = mapper[keys]
        except KeyError:
            #raise KeyError('table_name=%s keys=%s' % (self.table_name_str, str(keys)))
            return n
        if self.is_debug_file:
            self.binary_debug.write('  found keys=%s -> name=%-6s - %s\n' % (
                str(keys), name, self.table_name))
        if self.debug:
            self.log.debug("  found keys=(%5s,%4s,%4s) name=%-6s - %s" % (
                keys[0], keys[1], keys[2], name, self.table_name))
        self.card_name = name
        n = func(data, n)  # gets all the grid/mat cards
        assert n is not None, name
        del self.card_name

        self.geom_keys = keys
        self.is_start_of_subtable = False
        self.isubtable_old = self.isubtable

        #assert n == len(data), 'n=%s len(data)=%s' % (n, len(data))
        return n

    def _fix_format_code(self, format_code=1):
        """
        Nastran can mess up the format code by using what the user specified,
        which may be wrong.

        For a SOL 101, if the user uses the following in their BDF:
            DISP(PLOT,PHASE)=ALL
        it's wrong, and should be:
            DISP(PLOT,REAL)=ALL

        """
        # we'll probably remove this later because we're fixing
        #it before we get to the object
        return
        # if self.format_code != format_code:
            # self.format_code = format_code
            # self.obj.format_code = format_code
            # self.obj.data_code['format_code'] = format_code

    def _read_random_table(self, data, ndata, result_name, storage_obj,
                           real_vector, node_elem,
                           random_code=None, is_cid=False):
        """
        Reads a real table (for random analysis)
        """
        assert self.format_code == 1, self.format_code
        assert self.num_wide == 8, self.num_wide
        is_vectorized = True

        #if self.format_code == 1 and self.num_wide == 8:  # real/random
        # real
        ntotal = 32 * self.factor
        nnodes = ndata // ntotal  # 8*4
        #self.log.debug('  create table_vector')
        auto_return = self._create_table_vector(
            result_name, nnodes, storage_obj, real_vector, is_cid=is_cid)
        if auto_return:
            return ndata
        #self.log.debug('  *create table_vector')

        #self._fix_format_code(format_code=1)
        is_sort1 = self.is_sort1  # uses the sort_bits

        if is_sort1:
            #self.log.debug('   sort1; table_name=%r' % self.table_name)
            if self.nonlinear_factor in (None, np.nan):
                n = self._read_real_table_static(
                    data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
            else:
                n = self._read_real_table_sort1(
                    data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
        else:
            #self.log.debug('   sort2; table_name=%r' % self.table_name)
            n = self._read_real_table_sort2(
                data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
        assert n is not None
        return n

    def _read_table_sort1_real(self, data, ndata, result_name, storage_obj,
                               real_vector, node_elem,
                               random_code=None, is_cid=False):
        """Reads a real table (for random analysis)"""
        assert self.format_code == 1, self.format_code
        assert self.num_wide == 8, self.num_wide
        is_vectorized = True

        #if self.format_code == 1 and self.num_wide == 8:  # real/random
        # real
        ntotal = 32 * self.factor
        nnodes = ndata // ntotal  # 8*4
        self.log.debug('  create table_vector')
        auto_return = self._create_table_vector(
            result_name, nnodes, storage_obj, real_vector, is_cid=is_cid)
        if auto_return:
            return ndata
        self.log.debug('  *create table_vector')

        #self._fix_format_code(format_code=1)
        self.log.debug('   sort1; table_name=%r' % self.table_name)
        if self.nonlinear_factor in (None, np.nan):
            n = self._read_real_table_static(
                data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
        else:
            n = self._read_real_table_sort1(
                data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
        assert n is not None
        return n

    def _read_table_vectorized(self, data, ndata, result_name, storage_obj,
                               real_vector, complex_vector,
                               node_elem, random_code=None, is_cid=False):
        """Reads a generalized real/complex SORT1/SORT2 table"""
        assert isinstance(result_name, str), 'result_name=%r' % result_name
        assert isinstance(storage_obj, dict), 'storage_obj=%r' % storage_obj
        #print('self.num_wide =', self.num_wide)
        #print('random...%s' % self.isRandomResponse())
        #if not self.isRandomResponse():
        is_vectorized = True
        factor = self.factor
        if self.format_code == 1 and self.num_wide == 8:  # real/random
            # real
            ntotal = 32 * factor
            nnodes = ndata // ntotal  # 8*4
            auto_return = self._create_table_vector(
                result_name, nnodes, storage_obj, real_vector, is_cid=is_cid)
            if auto_return:
                return ndata

            self._fix_format_code(format_code=1)
            if self.is_sort1:
                if self.nonlinear_factor in (None, np.nan):
                    n = self._read_real_table_static(
                        data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                else:
                    n = self._read_real_table_sort1(
                        data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
            else:
                n = self._read_real_table_sort2(
                    data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                #n = ndata
                #msg = self.code_information()
                #n = self._not_implemented_or_skip(data, ndata, msg)
        elif self.format_code in [2, 3] and self.num_wide == 14:  # real or real/imaginary or mag/phase
            # complex
            ntotal = 56 * factor
            nnodes = ndata // ntotal  # 14*4
            if self.is_debug_file:
                self.binary_debug.write('nnodes=%s' % nnodes)
            auto_return = self._create_table_vector(
                result_name, nnodes, storage_obj, complex_vector)
            if auto_return:
                return ndata
            if self.is_sort1:
                if self.is_magnitude_phase():
                    n = self._read_complex_table_sort1_mag(
                        data, is_vectorized, nnodes, result_name, node_elem)
                else:
                    n = self._read_complex_table_sort1_imag(
                        data, is_vectorized, nnodes, result_name, node_elem)
            else:
                if self.is_magnitude_phase():
                    n = self._read_complex_table_sort2_mag(
                        data, is_vectorized, nnodes, result_name, node_elem)
                else:
                    n = self._read_complex_table_sort2_imag(
                        data, is_vectorized, nnodes, result_name, node_elem)
                #msg = self.code_information()
                #n = self._not_implemented_or_skip(data, ndata, msg)
        else:
            #msg = 'COMPLEX/PHASE is included in:\n'
            #msg += '  DISP(PLOT)=ALL\n'
            #msg += '  but the result type is REAL\n'
            msg = self.code_information()
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_scalar_table_vectorized(self, data, ndata, result_name, storage_obj,
                                      real_vector, complex_vector,
                                      node_elem, random_code=None, is_cid=False):
        """
        Reads a table

        Parameters
        ----------
        data : bytes
            the data to read
        ndata : int
            the length of data
        result_name : str
            the name
        storage_obj : dict
            the slot for the result
        real_vector : RealTableArray()
            the result object if this is a real result
        complex_vector : ComplexTableArray()
            the result object if this is a complex result
        node_elem : str
            'node' or 'elem'
        random_code : int; default=None
            unused
        is_cid : bool; default=False
            unused

        Returns
        -------
        n : int
            the new position in the OP2

        >>> n = self._read_scalar_table_vectorized(
            data, ndata, result_name, storage_obj,
            RealTemperatureVectorArray, ComplexThermalLoadVectorArray,
            'node', random_code=self.random_code)

        """
        assert isinstance(result_name, str), 'result_name=%r' % result_name
        assert isinstance(storage_obj, dict), 'storage_obj=%r' % storage_obj
        #print('self.num_wide =', self.num_wide)
        #print('random...%s' % self.isRandomResponse())
        #if not self.isRandomResponse():
        is_vectorized = True
        factor = self.factor
        if self.format_code == 1 and self.num_wide == 8:  # real/random
            # real
            ntotal = 32 * factor
            nnodes = ndata // ntotal  # 8*4
            auto_return = self._create_table_vector(
                result_name, nnodes, storage_obj, real_vector, is_cid=is_cid)
            if auto_return:
                return ndata

            self._fix_format_code(format_code=1)
            if self.is_sort1:
                if self.nonlinear_factor in (None, np.nan):
                    n = self._read_real_scalar_table_static(
                        data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                else:
                    n = self._read_real_scalar_table_sort1(
                        data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
            else:
                n = self._read_real_scalar_table_sort2(
                    data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                #n = ndata
                #msg = self.code_information()
                #n = self._not_implemented_or_skip(data, ndata, msg)
        elif self.format_code in [2, 3] and self.num_wide == 14:  # real or real/imaginary or mag/phase
            raise NotImplementedError('real/imaginary or mag/phase')
            # complex
            ntotal = 56 * factor
            nnodes = ndata // ntotal  # 14*4
            if self.is_debug_file:
                self.binary_debug.write('nnodes=%s' % nnodes)
            auto_return = self._create_table_vector(
                result_name, nnodes, storage_obj, complex_vector)
            if auto_return:
                return ndata
            if self.is_sort1:
                if self.is_magnitude_phase():
                    n = self._read_complex_table_sort1_mag(data, is_vectorized, nnodes, result_name, node_elem)
                else:
                    n = self._read_complex_table_sort1_imag(data, is_vectorized, nnodes, result_name, node_elem)
            else:
                if self.is_magnitude_phase():
                    n = self._read_complex_table_sort2_mag(data, is_vectorized, nnodes, result_name, node_elem)
                else:
                    n = self._read_complex_table_sort2_imag(data, is_vectorized, nnodes, result_name, node_elem)
                #msg = self.code_information()
                #n = self._not_implemented_or_skip(data, ndata, msg)
        else:
            #msg = 'COMPLEX/PHASE is included in:\n'
            #msg += '  DISP(PLOT)=ALL\n'
            #msg += '  but the result type is REAL\n'
            msg = self.code_information()
            n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def function_code(self, value):
        """
        This is a new specification from NX that's really important and
        not in the MSC manual, even though they use it.

        ACODE,4=05 vs. ACODE=05  -> function code 4
        TCODE,1=02 vs. TCODE=02  -> function code 1

        """
        if self._function_code == 1:
            if value // 1000 in [2, 3, 6]:
                out = 2
            else:
                out = 1
        elif self._function_code == 2:
            out = value % 100
        elif self._function_code == 3:
            out = value % 1000
        elif self._function_code == 4:
            out = value // 10
        elif self._function_code == 5:
            out = value % 10
        #elif self._function_code == 6:
            #raise NotImplementedError(self.function_code)
        #elif self._function_code == 7:
            #raise NotImplementedError(self.function_code)
        else:
            raise NotImplementedError(self.function_code)
        return out

    def _read_real_scalar_table_static(self, data, is_vectorized, nnodes,
                                       unused_result_name, flag, is_cid=False):
        """
        With a static (e.g. SOL 101) result, reads a complex OUG-style
        table created by:
              DISP(PLOT,SORT1,REAL) = ALL

        """
        if self.is_debug_file:
            self.binary_debug.write('  _read_real_scalar_table_static\n')
        assert flag in ['node', 'elem'], flag
        dt = self.nonlinear_factor
        assert self.obj is not None
        obj = self.obj

        factor = self.factor
        ntotal = 32 * factor # 32=4*8
        if self.use_vector and is_vectorized:
            n = nnodes * ntotal
            itotal2 = obj.itotal + nnodes
            #print('ndata=%s n=%s nnodes=%s' % (ndata, n, nnodes))
            ints = np.frombuffer(data, dtype=self.idtype).reshape(nnodes, 8)
            floats = np.frombuffer(data, dtype=self.fdtype).reshape(nnodes, 8)
            obj._times[obj.itime] = dt
            #self.node_gridtype[self.itotal, :] = [node_id, grid_type]
            #self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]
            #obj.node_gridtype[obj.itotal:itotal2, :] = ints[:, 0:1]
            nids = ints[:, 0] // 10
            assert nids.min() > 0, nids.min()
            obj.node_gridtype[obj.itotal:itotal2, 0] = nids
            obj.node_gridtype[obj.itotal:itotal2, 1] = ints[:, 1].copy()
            obj.data[obj.itime, obj.itotal:itotal2, 0] = floats[:, 2].copy()
            if np.abs(floats[:, 1:]).max() != 0:
                msg = '%s is not a scalar result...do you have p-elements?\n' % (
                    obj.__class__.__name__)
                for icol in range(1, 6):
                    abs_max = np.abs(floats[:, icol]).max()
                    if abs_max != 0:
                        msg += 'itime=%s icol=%s max=%s min=%s\n' % (
                            obj.itime, icol, floats[:, icol].max(), floats[:, icol].min())
                self.log.warning(msg.rstrip())
                #raise ValueError(msg.rstrip())
            obj.itotal = itotal2
        else:
            dt = np.nan
            n = 0
            fmt = mapfmt(self._endian + b'2i6f', self.size)
            s = Struct(fmt)
            for unused_inode in range(nnodes):
                out = s.unpack(data[n:n+ntotal])
                eid_device, grid_type, tx = out[:3]
                eid = eid_device // 10
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
                obj.add_sort1(dt, eid, grid_type, tx)
                n += ntotal
        return n

    def _read_real_scalar_table_sort1(self, data, is_vectorized, nnodes,
                                      unused_result_name, flag, is_cid=False):
        """
        With a real transient result (e.g. SOL 109/159), reads a
        real OUG-style table created by:
              DISP(PLOT,SORT1,REAL) = ALL

        """
        #print('result_name=%s use_vector=%s is_vectorized=%s' % (
            #result_name, self.use_vector, is_vectorized))
        if self.is_debug_file:
            self.binary_debug.write('  _read_real_scalar_table_sort1\n')
        #assert flag in ['node', 'elem'], flag
        #assert self.obj is not None
        dt = self.nonlinear_factor
        obj = self.obj
        ntotal = 32 * self.factor # 32=4 * 8
        if self.use_vector and is_vectorized:
            itime = obj.itime
            n = nnodes * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nnodes

            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=self.idtype).reshape(nnodes, 8)
                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[itotal:itotal2, 0] = nids
                obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1].copy()

            floats = np.frombuffer(data, dtype=self.fdtype).reshape(nnodes, 8)
            obj.data[obj.itime, obj.itotal:itotal2, 0] = floats[:, 2].copy()
            assert np.abs(floats[:, 3:]).max() == 0, '%s is not a scalar result...' % obj.__class__.__name__
            obj._times[itime] = dt
            obj.itotal = itotal2
        else:
            n = 0
            assert nnodes > 0, nnodes
            fmt = mapfmt(self._endian + b'2i6f', self.size)
            s = Struct(fmt)
            for unused_inode in range(nnodes):
                out = s.unpack(data[n:n+ntotal])
                eid_device, grid_type, tx = out[:3]
                eid = eid_device // 10
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
                obj.add_sort1(dt, eid, grid_type, tx)
                n += ntotal
        return n

    def _read_real_scalar_table_sort2(self, data, is_vectorized, nnodes, result_name,
                                      flag, is_cid=False):
        """
        With a real transient result (e.g. SOL 109/159), reads a
        real OUG-style table created by:
              DISP(PLOT,SORT2,REAL) = ALL

        """
        if self.is_debug_file:
            self.binary_debug.write('  _read_real_scalar_table_sort2\n')
        assert flag in ['node', 'elem'], flag
        eid = self.nonlinear_factor
        #assert self.obj is not None

        obj = self.obj
        if self.use_vector and is_vectorized and 0:  # TODO: not done....
            itime = obj.itime
            n = nnodes * 4 * 8
            itotal = obj.itotal
            itotal2 = itotal + nnodes
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=self.idtype).reshape(nnodes, 8)
                #nids = ints[:, 0] // 10
                nids = np.ones(nnodes, dtype='int32') * eid
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[itotal:itotal2, 0] = nids
                obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1].copy()

            floats = np.frombuffer(data, dtype=self.fdtype).reshape(nnodes, 8).copy()
            obj._times[itime] = floats[:, 0]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 2]
            assert np.abs(floats[:, 3:]).max() == 0, '%s is not a scalar result...' % obj.__class__.__name__
            obj.itotal = itotal2
        else:
            n = 0
            assert nnodes > 0

            flag = 'freq/dt/mode'
            s = Struct(self._endian + self._analysis_code_fmt + b'i6f')
            assert eid > 0, self.code_information()
            for unused_inode in range(nnodes):
                edata = data[n:n+32]
                out = s.unpack(edata)
                (dt, grid_type, tx) = out[:3]
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, dt, str(out)))
                obj.add_sort2(dt, eid, grid_type, tx)
                n += 32
        return n

    def _read_real_table_static(self, data, is_vectorized, nnodes,
                                unused_result_name, flag, is_cid=False):
        """
        With a static (e.g. SOL 101) result, reads a complex OUG-style
        table created by:
              DISP(PLOT,SORT1,REAL) = ALL

        """
        if self.is_debug_file:
            self.binary_debug.write('  _read_real_table_static\n')
        assert flag in ['node', 'elem'], flag
        dt = self.nonlinear_factor
        assert self.obj is not None
        obj = self.obj

        ntotal = 32 * self.factor # 4 * 8
        if self.use_vector and is_vectorized:
            n = nnodes * ntotal
            itotal2 = obj.itotal + nnodes
            #print('ndata=%s n=%s nnodes=%s' % (ndata, n, nnodes))
            ints = np.frombuffer(data, dtype=self.idtype8).reshape(nnodes, 8)
            floats = np.frombuffer(data, dtype=self.fdtype8).reshape(nnodes, 8)
            obj._times[obj.itime] = dt
            #self.node_gridtype[self.itotal, :] = [node_id, grid_type]
            #self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]
            #obj.node_gridtype[obj.itotal:itotal2, :] = ints[:, 0:1]
            nids = ints[:, 0] // 10
            assert nids.min() > 0, nids.min()
            obj.node_gridtype[obj.itotal:itotal2, 0] = nids
            obj.node_gridtype[obj.itotal:itotal2, 1] = ints[:, 1].copy()
            obj.data[obj.itime, obj.itotal:itotal2, :] = floats[:, 2:].copy()
            obj.itotal = itotal2
        else:
            n = 0
            dt = np.nan
            fmt = mapfmt(self._endian + b'2i6f', self.size)
            s = Struct(fmt)
            for unused_inode in range(nnodes):
                out = s.unpack(data[n:n+ntotal])
                (eid_device, grid_type, tx, ty, tz, rx, ry, rz) = out
                eid = eid_device // 10
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
                obj.add_sort1(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
                n += ntotal
        return n

    def _read_real_table_sort1(self, data, is_vectorized, nnodes,
                               unused_result_name, flag, is_cid=False):
        """
        With a real transient result (e.g. SOL 109/159), reads a
        real OUG-style table created by:
              DISP(PLOT,SORT1,REAL) = ALL

        """
        #print('result_name=%s use_vector=%s is_vectorized=%s' % (
            #result_name, self.use_vector, is_vectorized))
        if self.is_debug_file:
            self.binary_debug.write('  _read_real_table_sort1\n')
        #assert flag in ['node', 'elem'], flag
        #assert self.obj is not None
        dt = self.nonlinear_factor
        obj = self.obj

        ntotal = 32 * self.factor  # 32=4*8
        if self.use_vector and is_vectorized:
            itime = obj.itime
            n = nnodes * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nnodes

            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=self.idtype8).reshape(nnodes, 8)

                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[itotal:itotal2, 0] = nids
                obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1].copy()

            floats = np.frombuffer(data, dtype=self.fdtype8).reshape(nnodes, 8)
            obj.data[obj.itime, obj.itotal:itotal2, :] = floats[:, 2:].copy()
            obj._times[itime] = dt
            obj.itotal = itotal2
        else:
            n = 0
            assert nnodes > 0, nnodes
            fmt = mapfmt(self._endian + b'2i6f', self.size)
            s = Struct(fmt)
            for unused_inode in range(nnodes):
                out = s.unpack(data[n:n+ntotal])
                (eid_device, grid_type, tx, ty, tz, rx, ry, rz) = out
                eid = eid_device // 10
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
                obj.add_sort1(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
                n += ntotal
        return n

    def _read_real_table_sort2(self, data, is_vectorized, nnodes, result_name, flag, is_cid=False):
        """
        With a real transient result (e.g. SOL 109/159), reads a
        real OUG-style table created by:
              DISP(PLOT,SORT2,REAL) = ALL

        """
        if self.is_debug_file:
            self.binary_debug.write('  _read_real_table_sort2\n')
        assert flag in ['node', 'elem'], flag
        nid = self.nonlinear_factor
        #assert self.obj is not None

        obj = self.obj
        ntotal = 32 * self.factor # 4*8
        if self.use_vector and is_vectorized:
            itime = obj.itime
            n = nnodes * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nnodes

            floats = np.frombuffer(data, dtype=self.fdtype8).reshape(nnodes, 8).copy()
            ints = np.frombuffer(data, dtype=self.idtype8).reshape(nnodes, 8)

            self._set_sort2_time(obj, self._analysis_code_fmt, ints, floats)
        #def _set_sort2_time(self, obj, self._analysis_code_fmt, ints, floats):
            #if obj.itime == 0:
                #if self._analysis_code_fmt == b'i':
                    #times = ints[:, 0]
                #else:
                    #assert self._analysis_code_fmt == b'f'
                    #times = floats[:, 0]
                #obj._times = times
            obj.node_gridtype[itime, 0] = nid
            obj.node_gridtype[itime, 1] = ints[0, 1].copy()
            obj.data[itotal:itotal2, obj.itime, :] = floats[:, 2:]
            obj.itotal = itotal2
        else:
            n = 0
            assert nnodes > 0
            flag = self.data_code['analysis_method']
            fmt = mapfmt(self._endian + self._analysis_code_fmt + b'i6f', self.size)
            structi = Struct(fmt)

            #psds = ('CRM2', 'NO2', 'PSD2', 'RMS2')
            #print('sort_method=%s' % self.sort_method)
            #if self.table_name_str.endswith(psds):
            for unused_inode in range(nnodes):
                edata = data[n:n+ntotal]
                out = structi.unpack(edata)
                (dt, grid_type, tx, ty, tz, rx, ry, rz) = out
                if self.is_debug_file:
                    self.binary_debug.write(
                        f'  nid={nid} {flag}={dt} ({type(dt)}); {str(out)}\n')
                obj.add_sort2(dt, nid, grid_type, tx, ty, tz, rx, ry, rz)
                n += ntotal
        #if self.table_name_str == 'OQMRMS1':
            #print(obj.node_gridtype)
            #print('------------')
        return n

    def _set_sort2_time(self, obj, analysis_code_fmt, ints, floats):
        if obj.itime == 0:
            if analysis_code_fmt == b'i':
                times = ints[:, 0]
            else:
                assert analysis_code_fmt == b'f'
                times = floats[:, 0]
            obj._times = times

    def _read_complex_table_sort1_mag(self, data, is_vectorized, nnodes, result_name, flag):
        if self.is_debug_file:
            self.binary_debug.write('  _read_complex_table_sort1_mag\n')

        #force_flux = self.get_force_flux(self.thermal)
        #disp_temp = self.get_disp_temp(self.thermal)
        #self.log.info('_read_complex_table_sort1_mag - %s; %s\n' % (
            #self.code, self.get_table_code_name(disp_temp, force_flux, stress_word='')))
        assert flag in ['node', 'elem'], flag
        dt = self.nonlinear_factor

        n = 0
        obj = self.obj
        ntotal = 56 * self.factor # 4 * 14
        if self.use_vector and is_vectorized:
            n = nnodes * ntotal
            itotal2 = obj.itotal + nnodes

            ints = np.frombuffer(data, dtype=self.idtype8).reshape(nnodes, 14)
            #print('ints[:, 0] =', ints[:, 0], ints[:, 0] // 10)
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=self.idtype8).reshape(nnodes, 14)
                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[obj.itotal:itotal2, 0] = nids
                obj.node_gridtype[obj.itotal:itotal2, 1] = ints[:, 1].copy()

            floats = np.frombuffer(data, dtype=self.fdtype8).reshape(nnodes, 14).copy()
            mag = floats[:, 2:8]
            phase = floats[:, 8:]
            rtheta = np.radians(phase)
            real_imag = mag * (np.cos(rtheta) + 1.j * np.sin(rtheta))
            #abs(real_imag), angle(real_imag, deg=True)

            obj._times[obj.itime] = dt
            obj.data[obj.itime, obj.itotal:itotal2, :] = real_imag
            obj.itotal = itotal2
        else:
            fmt = mapfmt(self._endian + b'2i12f', self.size)
            s = Struct(fmt)
            for unused_inode in range(nnodes):
                out = s.unpack(data[n:n+ntotal])
                (eid_device, grid_type, txr, tyr, tzr, rxr, ryr, rzr,
                 txi, tyi, tzi, rxi, ryi, rzi) = out
                eid = eid_device // 10
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i %s\n' % (flag, eid, str(out)))
                tx = polar_to_real_imag(txr, txi)
                ty = polar_to_real_imag(tyr, tyi)
                tz = polar_to_real_imag(tzr, tzi)
                rx = polar_to_real_imag(rxr, rxi)
                ry = polar_to_real_imag(ryr, ryi)
                rz = polar_to_real_imag(rzr, rzi)
                obj.add_sort1(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
                n += ntotal
        return n

    def _read_complex_table_sort1_imag(self, data, is_vectorized, nnodes,
                                       unused_result_name, flag):
        if self.is_debug_file:
            self.binary_debug.write('  _read_complex_table_sort1_imag\n')
        #assert flag in ['node', 'elem'], flag
        dt = self.nonlinear_factor
        obj = self.obj

        ntotal = 56 * self.factor # 4 * 14
        if self.use_vector and is_vectorized:
            n = nnodes * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nnodes

            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=self.idtype8).reshape(nnodes, 14)
                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                try:
                    obj.node_gridtype[itotal:itotal2, 0] = nids
                except ValueError:  # pragma: no cover
                    msg = f'nnids={len(nids)} itotal={itotal} itotal2={itotal2}'
                    print(obj.node_gridtype[:, 0].shape)
                    print(obj.node_gridtype[itotal:itotal2, 0].shape)
                    print(nids.shape)
                    raise ValueError(msg)
                obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1].copy()

            floats = np.frombuffer(data, dtype=self.fdtype8).reshape(nnodes, 14).copy()
            real = floats[:, 2:8]
            imag = floats[:, 8:]

            obj._times[obj.itime] = dt
            obj.data[obj.itime, itotal:itotal2, :] = real + 1.j * imag
            obj.itotal = itotal2
        else:
            n = 0
            fmt = mapfmt(self._endian + b'2i12f', self.size)
            s = Struct(fmt)

            assert self.obj is not None
            assert nnodes > 0
            for unused_inode in range(nnodes):
                out = s.unpack(data[n:n+ntotal])
                (eid_device, grid_type, txr, tyr, tzr, rxr, ryr, rzr,
                 txi, tyi, tzi, rxi, ryi, rzi) = out
                eid = eid_device // 10
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i %s\n' % (flag, eid, str(out)))
                tx = complex(txr, txi)
                ty = complex(tyr, tyi)
                tz = complex(tzr, tzi)
                rx = complex(rxr, rxi)
                ry = complex(ryr, ryi)
                rz = complex(rzr, rzi)
                obj.add_sort1(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
                n += ntotal
        return n

    def _check_id(self, eid_device, unused_flag, unused_out):
        """
        Somewhat risky method for calculating the eid because the device code
        is ignored.  However, this might be the actual way to parse the id.

        """
        #print('eid =', eid)
        #print('flag =', flag)
        eid2 = eid_device // 10
        return eid2
        # eid = (eid_device - self.device_code) // 10
        # if eid != eid2 or eid2 <= 0:
            # msg = 'eid_device=%s device_code=%s eid=%s eid2=%s\n\n' % (eid_device, self.device_code,
                                                                       # eid, eid2)
            # msg += 'The device code is set wrong, probably because you used:\n'
            # msg += "  '%s=ALL' instead of '%s(PLOT,PRINT,REAL)=ALL'"  % (bdf_name, bdf_name)
            # msg += '  %s=%i; %s\n' % (flag, eid, str(out))
            # msg += str(self.code_information())
            # raise DeviceCodeError(msg)
        # return eid2

    def get_oug2_flag(self):
        if self.analysis_code == 5:
            flag = 'freq'
            flag_type = '%.2f'
        else:
            raise RuntimeError(self.code_information())
        #flag = 'freq/dt/mode'
        return flag, flag_type

    def _read_complex_table_sort2_mag(self, data, is_vectorized, nnodes, result_name, flag):
        if self.is_debug_file:
            self.binary_debug.write('  _read_complex_table_sort2_mag\n')
        #self.log.info('_read_complex_table_sort2_mag')
        assert flag in ['node', 'elem'], flag
        flag, flag_type = self.get_oug2_flag()
        node_id = self.nonlinear_factor

        #ntotal = 56  # 14 * 4
        assert self.obj is not None
        assert nnodes > 0
        #assert ndata % ntotal == 0

        obj = self.obj
        ntotal = 56 * self.factor # 4 * 14
        if self.use_vector and is_vectorized:
            itime = obj.itime
            n = nnodes * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nnodes

            floats =np.frombuffer(data, dtype=self.fdtype8).reshape(nnodes, 14).copy()
            ints = np.frombuffer(data, dtype=self.idtype8).reshape(nnodes, 14)

            self._set_sort2_time(obj, self._analysis_code_fmt, ints, floats)
            obj.node_gridtype[itime, 0] = node_id
            obj.node_gridtype[itime, 1] = ints[0, 1].copy()

            mag = floats[:, 2:8]
            phase = floats[:, 8:]
            rtheta = np.radians(phase)
            real_imag = mag * (np.cos(rtheta) + 1.j * np.sin(rtheta))
            obj.data[itotal:itotal2, obj.itime, :] = real_imag
            obj.itotal = itotal2
        else:
            n = 0
            s = Struct(mapfmt(self._endian + self._analysis_code_fmt + b'i12f', self.size))
            binary_debug_fmt = '  %s=%s %%s\n' % (flag, flag_type)
            for unused_inode in range(nnodes):
                edata = data[n:n+ntotal]
                out = s.unpack(edata)
                (freq, grid_type, txr, tyr, tzr, rxr, ryr, rzr,
                 txi, tyi, tzi, rxi, ryi, rzi) = out

                if self.is_debug_file:
                    self.binary_debug.write(binary_debug_fmt % (freq, str(out)))
                tx = polar_to_real_imag(txr, txi)
                ty = polar_to_real_imag(tyr, tyi)
                tz = polar_to_real_imag(tzr, tzi)
                rx = polar_to_real_imag(rxr, rxi)
                ry = polar_to_real_imag(ryr, ryi)
                rz = polar_to_real_imag(rzr, rzi)
                obj.add_sort2(freq, node_id, grid_type, tx, ty, tz, rx, ry, rz)
                n += ntotal
        return n

    def _read_complex_table_sort2_imag(self, data, is_vectorized, nnodes, result_name, flag):
        """
        With a complex result (e.g. SOL 103/108), reads a complex OUG-style
        table created by:
          DISP(PLOT,SORT2,PHASE) = ALL
          DISP(PLOT,SORT2,IMAG) = ALL

        """
        if self.is_debug_file:
            self.binary_debug.write('  _read_complex_table_sort2_imag\n')
        #self.log.info('_read_complex_table_sort2_imag')

        assert flag in ['node', 'elem'], flag
        flag, flag_type = self.get_oug2_flag()
        node_id = self.nonlinear_factor

        obj = self.obj
        ntotal = 56 * self.factor # 4 * 14
        if self.use_vector and is_vectorized:
            itime = obj.itime
            n = nnodes * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nnodes

            floats = np.frombuffer(data, dtype=self.fdtype8).reshape(nnodes, 14).copy()
            ints = np.frombuffer(data, dtype=self.idtype8).reshape(nnodes, 14)

            self._set_sort2_time(obj, self._analysis_code_fmt, ints, floats)
            obj.node_gridtype[itime, 0] = node_id
            obj.node_gridtype[itime, 1] = ints[0, 1].copy()

            real = floats[:, 2:8]
            imag = floats[:, 8:]
            obj.data[itotal:itotal2, obj.itime, :] = real + 1.j * imag
            obj.itotal = itotal2
        else:
            n = 0
            #ntotal = 56  # 14 * 4
            fmt = mapfmt(self._endian + self._analysis_code_fmt + b'i12f', self.size)
            s = Struct(fmt)
            assert self.obj is not None
            assert nnodes > 0
            #assert ndata % ntotal == 0

            binary_debug_fmt = '  %s=%s %%s\n' % (flag, flag_type)

            for unused_inode in range(nnodes):
                edata = data[n:n+ntotal]
                out = s.unpack(edata)

                (freq, grid_type, txr, tyr, tzr, rxr, ryr, rzr,
                 txi, tyi, tzi, rxi, ryi, rzi) = out

                if self.is_debug_file:
                    self.binary_debug.write(binary_debug_fmt % (freq, str(out)))
                tx = complex(txr, txi)
                ty = complex(tyr, tyi)
                tz = complex(tzr, tzi)
                rx = complex(rxr, rxi)
                ry = complex(ryr, ryi)
                rz = complex(rzr, rzi)
                obj.add_sort2(freq, node_id, grid_type, tx, ty, tz, rx, ry, rz)
                n += ntotal
        return n

    def create_transient_object(self, result_name, storage_obj, class_obj,
                                is_cid=False, debug=False):
        """
        Creates a transient object (or None if the subcase should be skippied).

        Parameters
        ----------
        result_name : str
            the name of the dictionary to store the object in (e.g. 'displacements')
        storage_obj : dict
            the dictionary to store the object in (e.g. self.displacements)
        class_obj : object()
            the class object to instantiate
        debug : bool
            developer debug

        .. python ::

            result_name = 'displacements'
            slot = self.displacements
            class_obj = RealDisplacementArray
            self.create_transient_object(result_name, storage_obj, class_obj, is_cid=is_cid)

        .. note:: dt can also be load_step depending on the class

        """
        assert not isinstance(class_obj, str), 'class_obj=%r' % class_obj
        assert class_obj is not None, class_obj
        if debug:
            print("create Transient Object")
            print("***NF = %s" % self.nonlinear_factor)
        #if not hasattr(self, storageName):
            #attrs =  object_attributes(obj, mode="public")
            #msg = 'storage_obj=%r does not exist.\n' % storage_obj
            #msg += 'Attributes = [%s]' , ', %s'.join(attrs)
            #raise RuntimeError(msg)
        #storage_obj = getattr(self, storageName)
        #assert class_obj is not None, 'name=%r has no associated classObject' % storageName

        #self.log.debug('self.table_name=%s isubcase=%s subtitle=%r' % (
            #self.table_name, self.isubcase, self.subtitle.strip()))
        table_name = self.table_name.decode(self.encoding)
        self.data_code['table_name'] = table_name
        self.data_code['result_name'] = result_name
        self.data_code['_count'] = self._count
        if 'h5_file' in self.data_code:
            del self.data_code['h5_file']
        assert self.log is not None

        code = self._get_code()
        #print('code =', code)
        if hasattr(self, 'isubcase'):
            if self.code in storage_obj:
                self.obj = storage_obj[code]
                if self.nonlinear_factor not in (None, np.nan):
                    if self.obj.nonlinear_factor in (None, np.nan):
                        msg = (
                            'The %s object is flipping from a \n'
                            'static (e.g. preload) result to a transient/frequency based results\n'
                            '%s -> %s\n' % (
                                result_name, self.obj.nonlinear_factor, self.nonlinear_factor))
                        msg += ('code = (subcase=%s, analysis_code=%s, sort=%s, count=%s, '
                                'ogs=%s, superelement_adaptivity_index=%r pval_step=%r)\n' % tuple(code))
                        msg += '%s\n' % str(self.obj)
                        msg += '\nIf this isnt correct, check if the data code was applied on the object'
                        raise MultipleSolutionNotImplementedError(msg)
                try:
                    data_codei = copy.deepcopy(self.data_code)
                except:
                    print("self.data_code =", self.data_code)
                    raise
                assert 'table_name' in data_codei
                assert data_codei['table_name'] is not None, data_codei
                self.obj.update_data_code(data_codei)
                assert self.obj.table_name is not None, self.data_code
            else:
                #if 'element_name' in self.data_code:
                    #print('code not in object', self.data_code['element_name'])
                #else:
                    #print('code not in object')
                class_obj.is_cid = is_cid
                is_sort1 = self.is_sort1  # uses the sort_bits

                #self.load_as_h5 = data_code['load_as_h5']
                #del data_code['load_as_h5']
            ##if 'h5_file' in data_code:
                #self.h5_file = data_code['h5_file']
                #del data_code['h5_file']

                if self.data_code['load_as_h5']:
                    self.data_code['h5_file'] = self.op2_reader.h5_file
                #if self.data_code['table_name'] not in ['OUGV1']:
                    #print(self.data_code)
                self.obj = class_obj(self.data_code, is_sort1, self.isubcase, self.nonlinear_factor)
                assert self.obj.table_name is not None, self.obj.data_code
            storage_obj[code] = self.obj
            #assert self.obj.table_name is not None
        else:
            if code in storage_obj:
                self.obj = storage_obj[code]
            else:
                storage_obj[code] = self.obj
        assert self.obj.table_name is not None, self.data_code

    def _get_code(self):
        """
        The code is a the way you access something like self.displacements.
        Ideally, it's just the subcase id, but for things like optimization, it
        gets more complicated.  So we make it in the complicated way and simplify
        it later if we can.

        """
        ogs = 0
        if hasattr(self, 'ogs'):
            ogs = self.ogs
        #if self.binary_debug:
            #self.binary_debug.write(self.code_information(include_time=True))

        code = (self.isubcase, self.analysis_code, self._sort_method, self._count, ogs,
                self.superelement_adaptivity_index, self.pval_step)
        #code = (self.isubcase, self.analysis_code, self._sort_method, self._count,
                #self.superelement_adaptivity_index, self.table_name_str)
        #print('%r' % self.subtitle)
        self.code = code
        #self.log.debug('code = %s' % str(self.code))
        return self.code

    def _not_implemented_or_skip(self, unused_data, ndata, msg=''):
        """
        A simple pass loop for unsupported tables that can be hacked on
        to crash the program everywhere that uses it.

        """
        #msg = 'table_name=%s table_code=%s %s\n%s' % (
            #self.table_name, self.table_code, msg, self.code_information())
        #if any([card_name in msg for card_name in ['VUHEXA', 'VUPENTA', 'VUTETRA', 'VUQUAD']]):
            #return ndata
        #raise NotImplementedError(msg)
        #if self.table_name_str.startswith(('OSTR', 'OES', 'OEF')):
            #if self.element_type in [145, 146, 147, # VUHEXA, VUPENTA, VUTETRA,
                                     #189, 190, 191, # VUQUAD, VUTRIA, VUBEAM
                                     #69, # CBEND
                                     #]:
                #return ndata

        if is_release:
            if msg != self._last_comment:
                #print(self.code_information())
                if self.read_mode == 2:
                    if msg == '':
                        self.log.warning(self.code_information())
                    self.log.warning(msg)
                    #ddd

                    #if not('VUHEXA' in msg or 'VUPENTA' in msg or 'VUTETRA' in msg
                           #or 'Element Stress' in self.code_information()
                           #or 'Element Strain' in self.code_information()):
                        #aaa
                #if self.table_name in ['OEFPSD1', 'OEFPSD2',
                #                       'OEFCRM1', 'OEFCRM2',
                #                       'OEFNO1', 'OEFNO2',
                #                       'OEFRMS1', 'OEFRMS2']:
                #    pass
                #else:
                #    self.log.warning(self.code_information())

                self._last_comment = msg
            return ndata
        else:  # pragma: no cover
            msg = 'table_name=%s table_code=%s %s\n%s' % (
                self.table_name, self.table_code, msg, self.code_information())
            raise NotImplementedError(msg)

    @property
    def size(self) -> int:
        return self.op2_reader.size

    @property
    def factor(self):
        return self.op2_reader.factor

    def parse_approach_code(self, data: bytes) -> None:
        if self.size == 4:
            self.parse_approach_code4(data)
        else:
            self.parse_approach_code8(data)


    def parse_approach_code4(self, data: bytes) -> None:
        """
        Function  Formula                                                Manual
        ========  =======                                                ======
        1         if item_name/1000 in 2,3,6: 2; else 1                  if(item_name/1000 = 2,3,6) then return 2, else return 1
        2         item_name % 100                                        mod(item_name,100)
        3         item_name % 1000                                       mod(item_name,1000)
        4         item_name // 10                                        item_name/10
        5         item_name % 10                                         mod(item_name,10)
        6         if item_name != 8; 0; else 1 # ???                     if iand(item_name,8)<> then set to 0, else set to 1
        7         if item_name in [0,2]: 0; elif item_name in [1,3] 1    if item_name/1000 = 0 or 2, then set to 0; = 1 or 3, then set to 1

        TCODE,1=02 means:
          TCODE1 = 2
          TCODE1/1000 = 0
          TCODE = f1(TCODE1)

        """
        (approach_code, tCode, int3, isubcase) = unpack(self._endian + b'4i', data[:16])
        self._set_approach_code(approach_code, tCode, int3, isubcase)

    def parse_approach_code8(self, data: bytes) -> None:
        (approach_code, tCode, int3, isubcase) = unpack(self._endian + b'4q', data[:32])
        self._set_approach_code(approach_code, tCode, int3, isubcase)

    def _set_approach_code(self, approach_code: int, tCode: int, int3: int, isubcase: int) -> None:
        self.approach_code = approach_code
        self.tCode = tCode
        self.int3 = int3
        self.data_code['size'] = self.size
        self.data_code['is_msc'] = self.is_msc
        self.data_code['is_nasa95'] = self.is_nasa95

        if not hasattr(self, 'subtable_name'):
            self.data_code['subtable_name'] = self.subtable_name

        self.data_code['table_name'] = self.table_name
        self.data_code['approach_code'] = approach_code

        #: the local subcase ID
        self.isubcase = isubcase
        self.data_code['isubcase'] = self.isubcase
        #self.subcases.add(self.isubcase)  # set notation

        table_code = tCode % 1000
        if table_code in NX_TABLES:
            self.to_nx()

        #: the type of result being processed
        self.table_code = table_code
        self.data_code['table_code'] = self.table_code
        self.data_code['tCode'] = self.tCode

        #: used to create sort_bits
        self.sort_code = tCode // 1000
        #Sort 1 - SortCode=((TCODE//1000)+2)//2

        self.data_code['sort_code'] = self.sort_code
        #print('tCode=%s tCode%%1000=%-2s tCode//1000=%s' % (tCode, tCode%1000, tCode//1000))
        self.sort_method = _function1(tCode)
        self.data_code['sort_method'] = self.sort_method

        #: what type of data was saved from the run; used to parse the
        #: approach_code and grid_device.  device_code defines what options
        #: inside a result, STRESS(PLOT,PRINT), are used.
        self.device_code = approach_code % 10

        self.data_code['device_code'] = self.device_code
        assert self.device_code in [0, 1, 2, 3, 4, 5, 6, 7], self.device_code

        #: what solution was run (e.g. Static/Transient/Modal)
        self.analysis_code = (approach_code - self.device_code) // 10

        self.data_code['analysis_code'] = self.analysis_code

        #print(
            #f'parse_approach_code - approach_code={approach_code} tCode={tCode} int3={int3} isubcase={isubcase}\n'
            #f'                 so - analysis_code={self.analysis_code} device_code={self.device_code} '
            #f'table_code={self.table_code} sort_code={self.sort_code}\n'
        #)
        if 0:
            if self.device_code == 3:
                #sys.stderr.write('The op2 may be inconsistent...\n')
                #sys.stderr.write("  print and plot can cause bad results..."
                #                 "if there's a crash, try plot only\n")
                self.device_code = 1

                #self.log.info('The op2 may be inconsistent...')
                #self.log.info('  print and plot can cause bad results...'
                #              'if there's a crash, try plot only')
                self.data_code['device_code'] = self.device_code

        if self.is_debug_file:
            self.binary_debug.write('  %-14s = %r\n' % ('table_name', self.table_name))
            self.binary_debug.write('  %-14s = analysis_code * 10 + device_code\n' % 'approach_code')
            self.binary_debug.write('  %-14s = %r\n' % ('approach_code', self.approach_code))
            self.binary_debug.write('  %-14s = %r\n' % ('  device_code', self.device_code))
            self.binary_debug.write('  %-14s = %r\n' % ('  analysis_code', self.analysis_code))
            self.binary_debug.write('  %-14s = sort_code * 1000 + table_code\n' % ('tCode'))
            self.binary_debug.write('  %-14s = %r\n' % ('tCode', self.tCode))
            self.binary_debug.write('  %-14s = %r\n' % ('  table_code', self.table_code))
            self.binary_debug.write('  %-14s = %r\n' % ('  sort_code', self.sort_code))
        self._parse_sort_code()
        assert self.sort_code in [0, 1, 2, 3, 4, 5, 6], self.sort_code #self.code_information()

    def _parse_thermal_code(self):
        """
        +------------+---------------+
        |  thermal   | thermal_bits  |
        +============+===============+
        | 0          |  [0, 0, 0]    |
        +------------+---------------+
        | 1          |  [0, 0, 1]    |
        +------------+---------------+
        | 2          |  [0, 1, 0]    |
        +------------+---------------+
        | 3          |  [0, 1, 1]    |
        +------------+---------------+
        | ...        | ...           |
        +------------+---------------+
        | 7          |  [1, 1, 1, 1] |
        +------------+---------------+

        1 Thermal
        2 Scaled response spectra ABS
        3 Scaled response spectra SRSS
        4 Scaled response spectra NRL
        5 Scaled response spectra NRLO
        ::
          thermal =  0 -> thermal_bits = [0,0,0,0,0]  # no thermal
          thermal =  1 -> thermal_bits = [0,0,0,0,1]  # 1- thermal
          thermal =  2 -> thermal_bits = [0,0,0,1,0]  # 2 - Scaled response spectra ABS
          thermal =  3 -> thermal_bits = [0,0,0,1,1]
          thermal =  4 -> thermal_bits = [0,0,1,0,0]  # 3 - Scaled response spectra SRSS
          thermal =  5 -> thermal_bits = [0,0,1,0,1]
          thermal =  6 -> thermal_bits = [0,0,1,1,0]
          thermal =  7 -> thermal_bits = [0,0,1,1,1]

          thermal =  8 -> thermal_bits = [0,1,0,0,0]  # 4-Scaled response spectra NRL
          thermal =  9 -> thermal_bits = [0,1,0,0,1]  # NRL + thermal
          thermal = 10 -> thermal_bits = [0,1,0,1,0]
          thermal = 11 -> thermal_bits = [0,1,0,1,1]
          thermal = 12 -> thermal_bits = [0,1,1,0,0]
          thermal = 13 -> thermal_bits = [0,1,1,0,1]
          thermal = 14 -> thermal_bits = [0,1,1,1,0]
          thermal = 15 -> thermal_bits = [0,1,1,1,1]

          #------
          thermal = 16 -> thermal_bits = [1,0,0,0,0]  # 5 - Scaled response spectra NRLO
          thermal = 17 -> thermal_bits = [1,0,0,0,1]
          thermal = 18 -> thermal_bits = [1,0,0,1,0]
          thermal = 19 -> thermal_bits = [1,0,0,1,1]
          thermal = 20 -> thermal_bits = [1,0,1,0,0]
          thermal = 21 -> thermal_bits = [1,0,1,0,1]
          thermal = 22 -> thermal_bits = [1,0,1,1,0]
          thermal = 23 -> thermal_bits = [1,0,1,1,1]

          thermal = 24 -> thermal_bits = [1,1,0,0,0]
          thermal = 25 -> thermal_bits = [1,1,0,0,1]
          thermal = 26 -> thermal_bits = [1,1,0,1,0]
          thermal = 27 -> thermal_bits = [1,1,0,1,1]
          thermal = 28 -> thermal_bits = [1,1,1,0,0]
          thermal = 29 -> thermal_bits = [1,1,1,0,1]
          thermal = 30 -> thermal_bits = [1,1,1,1,0]
          thermal = 31 -> thermal_bits = [1,1,1,1,1]


          thermal_bits[4] = 0 -> thermal
          thermal_bits[3] = 0 -> ABS
          thermal_bits[2] = 0 -> SRSS
          thermal_bits[1] = 0 -> NRL
          thermal_bits[0] = 0 -> NRLO

        """
        bits = [0, 0, 0, 0, 0]
        thermal_code = self.thermal

        # Sort codes can range from 0 to 7, but most of the examples
        # are covered by these.  The ones that break are incredibly large.
        #if self.thermal not in [0, 1, 2, 3, 4, 5, 6, 7]:
            #msg = 'Invalid sort_code=%s' % (self.sort_code)
            #raise SortCodeError(msg)
            #if self.sort_code == 1145655:
                #return
        i = 4
        while thermal_code > 0:
            value = thermal_code % 2
            thermal_code = (thermal_code - value) // 2
            bits[i] = value
            i -= 1

        #: the bytes describe the Random information
        self.thermal_bits = bits
        self.data_code['thermal_bits'] = self.thermal_bits

    def _parse_sort_code(self):
        """
        +------------+------------+
        | sort_code  | sort_bits  |
        +============+============+
        | 0          | [0, 0, 0]  |
        +------------+------------+
        | 1          | [0, 0, 1]  |
        +------------+------------+
        | 2          | [0, 1, 0]  |
        +------------+------------+
        | 3          | [0, 1, 1]  |
        +------------+------------+
        | ...        | ...        |
        +------------+------------+
        | 7          | [1, 1, 1]  |
        +------------+------------+

        ::
          sort_code = 0 -> sort_bits = [0,0,0]  #         sort1, real
          sort_code = 1 -> sort_bits = [0,0,1]  #         sort1, complex
          sort_code = 2 -> sort_bits = [0,1,0]  #         sort2, real
          sort_code = 3 -> sort_bits = [0,1,1]  #         sort2, complex
          sort_code = 4 -> sort_bits = [1,0,0]  # random, sort1, real
          sort_code = 5 -> sort_bits = [1,0,1]  # random, sort1, real
          sort_code = 6 -> sort_bits = [1,1,0]  # random, sort2, real
          sort_code = 7 -> sort_bits = [1,1,1]  # random, sort2, complex
          # random, sort2, complex <- [1, 1, 1]

          sort_bits[0] = 0 -> isSorted=True isRandom=False
          sort_bits[1] = 0 -> is_sort1=True is_sort2=False
          sort_bits[2] = 0 -> isReal=True   isReal/Imaginary=False

        """
        bits = [0, 0, 0]
        sort_code = self.sort_code

        # Sort codes can range from 0 to 7, but most of the examples
        # are covered by these.  The ones that break are incredibly large.
        if self.sort_code not in [0, 1, 2, 3, 4, 5, 6, 7]:
            msg = 'Invalid sort_code=%s' % (self.sort_code)
            raise SortCodeError(msg)
        i = 2
        while sort_code > 0:
            value = sort_code % 2
            sort_code = (sort_code - value) // 2
            bits[i] = value
            i -= 1

        # fixing bit[1]
        bits[1] = 0 if self.is_table_1 else 1
        #: the bytes describe the SORT information
        self.sort_bits = bits
        self.data_code['sort_bits'] = self.sort_bits

    @property
    def _sort_method(self):
        try:
            sort_method, unused_is_real, unused_is_random = self._table_specs()
        except:
            sort_method = get_sort_method_from_table_name(self.table_name)
        #is_sort1 = self.table_name.endswith('1')
        #is_sort1 = self.is_sort1  # uses the sort_bits
        assert sort_method in [1, 2], 'sort_method=%r\n%s' % (sort_method, self.code_information())
        return sort_method

    @property
    def is_real(self):
        unused_sort_method, is_real, unused_is_random = self._table_specs()
        return is_real

    @property
    def is_complex(self):
        return not self.is_real

    @property
    def is_random(self):
        unused_sort_method, unused_is_real, is_random = self._table_specs()
        return is_random

    #def is_mag_phase(self):
        #assert self.format_code in [0, 1], self.format_code
        #return bool(self.format_code)

    def is_mag_phase(self):
        return self.is_magnitude_phase()

    def is_magnitude_phase(self):
        if self.format_code == 3:
            return True
        return False

    def debug3(self):
        return self.is_debug_file
        #if self.debug and self.table_name in self.show_table3_map:
            #return True
        #return False

    def debug4(self):
        return self.is_debug_file
        #if self.debug and self.table_name in self.show_table4_map:
            #return True
        #return False

    #@property
    #def is_stress(self):
        #if self.stress_bits[1] == 0:
            #return True
        #return False

    @property
    def is_curvature(self):
        if self.is_stress:
            curvature_flag = False
        else:
            # strain only
            curvature_flag = self.stress_bits[2] == 0
        if self.s_code in [10, 11, 20, 27]:
            assert curvature_flag, curvature_flag
            return True
        assert not curvature_flag, curvature_flag
        return False

    @property
    def is_fiber_distance(self):
        return not self.is_curvature

    @property
    def is_max_shear(self):
        return self.stress_bits[4] == 0

    @property
    def is_von_mises(self):
        return not self.is_max_shear

    @property
    def is_stress(self):
        return not self.is_strain

    @property
    def is_strain(self):
        if self.stress_bits[1] == 1:
            return True
        return False

    def _create_table_object(self, result_name, nnodes,
                             slot, slot_object, slot_vector, is_cid=False):
        assert isinstance(result_name, str), result_name
        assert isinstance(slot, dict), slot
        auto_return = False
        #print('%s nnodes=%s' % (result_name, nnodes))
        is_vectorized = self.is_vectorized
        if is_vectorized and slot_vector is None:
            is_vectorized = False

        if is_vectorized:
            if self.read_mode == 1:
                self.create_transient_object(result_name, slot, slot_vector, is_cid=is_cid)
                self.result_names.add(result_name)
                self.obj._nnodes += nnodes
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                self.obj = slot[self.code]
                #self.obj.update_data_code(self.data_code)
                build_obj(self.obj)
        else:  # not vectorized
            self.result_names.add(result_name)
            if self.read_mode == 1:
                auto_return = True
                return auto_return, is_vectorized
            # pass = 0/2
            self.create_transient_object(result_name, slot, slot_object)
        return auto_return, is_vectorized

    def _create_table_vector(self, result_name, nnodes,
                             slot, slot_vector, is_cid=False):
        assert isinstance(result_name, str), result_name
        assert isinstance(slot, dict), slot
        auto_return = False
        #print('%s nnodes=%s' % (result_name, nnodes))
        self.result_names.add(result_name)
        if self.read_mode == 1:
            self.create_transient_object(result_name, slot, slot_vector, is_cid=is_cid)
            #self.result_names.add(result_name)
            self.obj._nnodes += nnodes
            auto_return = True
        elif self.read_mode == 2:
            self.code = self._get_code()
            self.obj = slot[self.code]
            #self.obj.update_data_code(self.data_code)
            build_obj(self.obj)
        else:
            auto_return = True
        return auto_return

    def _create_oes_object4(self, nelements, result_name, slot, obj_vector):
        """
        Creates the self.obj parameter based on if this is vectorized or not.

        Parameters
        ----------
        nelements :  int
            the number of elements to preallocate for vectorization
        result_name : str
            the name of the dictionary to store the object in (e.g. 'displacements')
        slot : dict[(int, int, str)=obj
            the self dictionary that will be filled with a
            non-vectorized result
        obj_vector : OESArray
            a pointer to the vectorized class

        Returns
        -------
        auto_return : bool
            a flag indicating a return n should be called
        is_vectorized : bool
            True/False

        Since that's confusing, let's say we have real CTETRA stress data.
        We're going to fill self.ctetra_stress with the class
        RealSolidStressArray.  So we call:

        if self._is_vectorized(RealSolidStressArray):
            if self._results.is_not_saved(result_vector_name):
                return ndata
        else:
            if self._results.is_not_saved(result_name):
                return ndata

        auto_return, is_vectorized = self._create_oes_object4(self, nelements,
                            'ctetra_stress', self.ctetra_stress,
                            RealSolidStressArray)
        if auto_return:
            return nelements * ntotal

        """
        auto_return = False
        #is_vectorized = True
        is_vectorized = self._is_vectorized(obj_vector)
        #print("vectorized...read_mode=%s...%s; %s" % (self.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if self.read_mode == 1:
                #print('oes-self.nonlinear_factor =', self.nonlinear_factor)
                #print(self.data_code)
                self.create_transient_object(result_name, slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % obj.ntimes)
                self.result_names.add(result_name)
                #print('self.obj =', self.obj)
                self.obj.nelements += nelements
                assert self.obj.table_name is not None
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                #self.log.info("code = %s" % str(self.code))
                #print("code = %s" % str(self.code))

                try:
                    self.obj = slot[self.code] # if this is failing, you probably set obj_vector to None...
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (self.code, result_name)
                    msg += "There's probably an extra check for read_mode=1...%s" % result_name
                    self.log.error(msg)
                    raise
                if not self.obj.table_name == self.table_name.decode('utf-8'):
                    print(self.obj)
                    msg = 'obj.table_name=%s table_name=%s; this shouldnt happen for read_mode=2' %  (
                        self.obj.table_name, self.table_name)
                    raise OverwriteTableError(msg)

                #obj.update_data_code(self.data_code)
                try:
                    build_obj(self.obj)
                except AssertionError:
                    print(self.code)
                    print(self.code_information())
                    raise

            else:  # not vectorized
                auto_return = True
        else:
            auto_return = True

        assert is_vectorized, '%r is not vectorized; obj=%s' % (result_name, obj_vector)
        return auto_return, is_vectorized

    def _is_vectorized(self, obj_vector):
        """
        Checks to see if the data array has been vectorized

        Parameters
        ----------
        obj_vector:  the object to check
            (obj or None; None happens when vectorization hasn't been implemented)

        Returns
        -------
        is_vectorized : bool
            should the data object be vectorized

        .. note :: the Vectorized column refers to the setting given by the user
        """
        is_vectorized = False
        if self.is_vectorized:
            if obj_vector is not None:
                is_vectorized = True
        return is_vectorized

    def _set_structs(self, size):
        """
        defines common struct formats

        https://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html#arrays-dtypes-constructing

        """
        self.fdtype = np.dtype(self._uendian + 'f4')
        self.idtype = np.dtype(self._uendian + 'i4')
        self.double_dtype = np.dtype(self._uendian + 'd')
        self.long_dtype = np.dtype(self._uendian + 'i8')
        #self.idtype = np.dtype(self._uendian + 'i8')

        #self.sdtype = np.dtype(self._uendian + '4s')
        self.struct_i = Struct(self._endian + b'i')
        self.struct_8s = Struct(self._endian + b'8s')
        if size == 4:
            self.idtype8 = self.idtype
            self.fdtype8 = self.fdtype
            self.struct_3i = Struct(self._endian + b'3i')
            self.struct_2i = Struct(self._endian + b'ii')
            self.struct_8s_i = Struct(self._endian + b'8si')
            #self.op2_reader.read_block = self.op2_reader.read_block4
            #self.op2_reader.read_markers = self.op2_reader.read_markers4
        elif size == 8:
            self.fdtype8 = np.dtype(self._uendian + 'd')
            self.idtype8 = np.dtype(self._uendian + 'i8')
            #self.fdtype8 = np.dtype(self._uendian + 'f4')
            #self.idtype8 = np.dtype(self._uendian + 'i4')

            self.log.warning('64-bit precison is poorly supported')
            #raise NotImplementedError('64-bit precison is not supported')
            self.struct_q = Struct(self._endian + b'q')
            self.struct_16s = Struct(self._endian + b'16s')
            #self.struct_8s = Struct(self._endian + b'8s')
            self.struct_2q = Struct(self._endian + b'2q')
            self.struct_16s_q = Struct(self._endian + b'16sq')
            #self.op2_reader.read_block = self.op2_reader.read_block8
            #self.op2_reader.read_markers = self.op2_reader.read_markers8

            # geom
            self.struct_3q = Struct(self._endian + b'3q')
        else:
            NotImplementedError(size)
        self.op2_reader.size = size
        self.op2_reader.factor = size // 4

    def del_structs(self):
        """deepcopy(OP2) fails on Python 3.6 without doing this"""
        del self.fdtype, self.idtype, self.double_dtype, self.long_dtype
        if hasattr(self, 'struct_2q'):
            del self.struct_i, self.struct_2q, self.struct_16s, self.struct_16s_q
            del self.struct_q, self.struct_8s
            if hasattr(self, 'struct_2i'):
                del self.struct_2i
            if hasattr(self, 'struct_3q'): # geom
                del self.struct_3q

        elif hasattr(self, 'struct_2i'):
            del self.struct_i, self.struct_2i, self.struct_3i, self.struct_8s, self.struct_8s_i
        out = [outi for outi in self.object_attributes() if 'struct_' in outi]
        assert len(out) == 0, out

def _cast_nonlinear_factor(value):
    """h5py is picky about it's data types"""
    if isinstance(value, int):
        try:
            value = np.int32(value)
        except OverflowError:
            value = np.int64(value)
    elif isinstance(value, float):
        value = np.float32(value)
    elif isinstance(value, (np.int32, np.float32)):  # pragma: no cover
        pass
    else: # pragma: no cover
        raise NotImplementedError('value=%s type=%s' % (value, type(value)))
    return value

def _function1(value):
    """function1(value)"""
    if value // 1000 in [2, 3, 6]:
        return 2
    return 1

def _function2(value):
    """function2(value)"""
    return value % 100

def _function3(value):
    """function3(value)"""
    return value % 1000

def _function4(value):
    """function4(value)"""
    return value // 10

def _function5(value):
    """function5(value)"""
    return value % 10

def _function6(value):
    """weird..."""
    if value != 8:
        return 0
    return 1

def _function7(value):
    """function7(value)"""
    if value in [0, 2]:
        out = 0
    elif value in [1, 3]:
        out = 1
    else:
        raise RuntimeError(value)
    return out
