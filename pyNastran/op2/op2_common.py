from __future__ import print_function, unicode_literals
from six import string_types, b
from six.moves import range, zip
import copy
from struct import Struct, unpack

from numpy import radians, sin, cos, fromstring, ones, float32, dtype as npdtype

import numpy as np

from pyNastran import is_release
from pyNastran.f06.f06_writer import F06Writer
from pyNastran.op2.op2_codes import Op2Codes, get_scode_word
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.errors import SortCodeError, MultipleSolutionNotImplementedError # DeviceCodeError,

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
        #    False : uses range loops
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
        self.result_names = set([])
        #: bool
        self.is_vectorized = None
        self.combine = False

        #: the storage dictionary that is passed to OP2 objects (e.g. RealDisplacementArray)
        #: the key-value pairs are extracted and used to generate dynamic self
        #: variables for the OP2 objects
        self.data_code = {}

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
        self.words = []

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
        self.get_element_type(33)

    def _device_code_(self):
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

    def fix_format_code(self):
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

        Note that a case of 4 is not used and is used below as a placeholder,
        while a case of -1 is some bizarre unhandled, undocumented case.
        """
        self._set_times_dtype()
        self.format_code_original = self.format_code
        if self.format_code == -1:
            if self.is_debug_file:
                self.write_ndata(self.binary_debug, 100)
            if self.table_name in [b'OESNLXR', b'OESNLBR', b'OESNLXD', b'OESNL1X']:
                assert self.format_code == -1, self.format_code
                self.format_code = 1
            else:
                raise RuntimeError(self.code_information())
            #return

        random_code = self.random_code if hasattr(self, 'random_code') else 0
        if random_code == 0:
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
                assert self.format_code in [1, 2, 3], self.code_information()
                self.format_code = 1
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

    def _set_times_dtype(self):
        self.data_code['_times_dtype'] = 'float32'
        if self.analysis_code == 1:   # statics / displacement / heat flux
            pass # static doesn't have a type
        elif self.analysis_code == 2:  # real eigenvalues
            pass
        #elif self.analysis_code==3: # differential stiffness
        #elif self.analysis_code==4: # differential stiffness
        elif self.analysis_code == 5:   # frequency
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

    def add_data_parameter(self, data, var_name, Type, field_num,
                           apply_nonlinear_factor=True, fix_device_code=False, add_to_dict=True):
        datai = data[4 * (field_num - 1) : 4 * (field_num)]
        assert len(datai) == 4, len(datai)
        value, = unpack(b(self._endian + Type), datai)
        if fix_device_code:
            value = (value - self.device_code) // 10
        if self.is_debug_file:
            self.binary_debug.write('  %-14s = %r\n' % (var_name, value))
        #setattr(self, var_name, value)  # set the parameter to the local namespace

        if apply_nonlinear_factor:
            self.nonlinear_factor = value
            #if self.table_name == b'OUGV2':
                #assert isinstance(self.nonlinear_factor, int), self.nonlinear_factor
            self.data_code['nonlinear_factor'] = value
            self.data_code['name'] = var_name

        if add_to_dict:
            self.data_code[var_name] = value

        try:
            self.words[field_num - 1] = var_name
        except IndexError:
            msg = 'Trying to set word, but...len(words)=%s ifield=%s' % (len(self.words), field_num - 1)
            raise IndexError(msg)
        return value

    def apply_data_code_value(self, name, value):
        self.data_code[name] = value

    def setNullNonlinearFactor(self):
        """
        initializes the nonlinear factor, which lets us know if
        this is a transient solution or not
        """
        self.nonlinear_factor = None
        self.data_code['nonlinear_factor'] = None

    def _read_title_helper(self, data):
        assert len(data) == 584, len(data)
        # titleSubtitleLabel
        title, subtitle, label = unpack(b(self._endian + '128s128s128s'), data[200:])
        self.title = title.strip().decode(self.encoding).strip()
        self.subtitle = subtitle.strip().decode(self.encoding).strip()
        self.label = label.strip().decode(self.encoding).strip()
        #print('title = %r' % self.title)

        #: the subtitle of the subcase
        self.data_code['subtitle'] = self.subtitle

        #: the label of the subcase
        self.data_code['label'] = self.label
        self.data_code['title'] = self.title#.decode(self.encoding)

        if self.is_debug_file:
            self.binary_debug.write('  %-14s = %r\n' % ('count', self._count))
            self.binary_debug.write('  %-14s = %r\n' % ('title', self.title))
            self.binary_debug.write('  %-14s = %r\n' % ('subtitle', self.subtitle))
            self.binary_debug.write('  %-14s = %r\n' % ('label', self.label))

    def _read_title(self, data):
        self._read_title_helper(data)

        if hasattr(self, 'isubcase'):
            if self.isubcase not in self.iSubcaseNameMap:
                self.iSubcaseNameMap[self.isubcase] = [self.subtitle, self.analysis_code, self.label]
        else:
            raise  RuntimeError('isubcase is not defined')

        if hasattr(self, 'subtitle') and hasattr(self, 'label'):
            code = (self.isubcase, self.analysis_code, self.subtitle)
            if code not in self.labels:
                self.subtitles[self.isubcase].append(self.subtitle)
                self.labels[code] = self.label

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
            for i, param in enumerate(self.words):
                if param == 's_code':
                    s_word = get_scode_word(self.s_code, self.stress_bits)
                    self.binary_debug.write('  s_code         = %s -> %s\n' % (self.s_code, s_word))
                    self.binary_debug.write('    stress_bits[0] = %i -> is_von_mises    =%-5s vs is_max_shear\n' % (self.stress_bits[0], self.is_von_mises()))
                    self.binary_debug.write('    stress_bits[1] = %i -> is_strain       =%-5s vs is_stress\n' % (self.stress_bits[1], self.is_strain()))
                    self.binary_debug.write('    stress_bits[2] = %i -> strain_curvature=%-5s vs fiber_dist\n' % (self.stress_bits[2], self.is_curvature()))
                    self.binary_debug.write('    stress_bits[3] = %i -> is_strain       =%-5s vs is_stress\n' % (self.stress_bits[3], self.is_strain()))
                    self.binary_debug.write('    stress_bits[4] = %i -> material coordinate system flag=%s vs ???\n' % (self.stress_bits[4], self.stress_bits[4]))
                elif param == '???':
                    param = 0
                msg += '%s, ' % param
                if i % 5 == 4:
                    msg += '\n             '

            if hasattr(self, 'format_code'):
                if self.is_complex():
                    self.binary_debug.write('\n  %-14s = %i -> is_mag_phase vs is_real_imag vs. is_random\n' % ('format_code', self.format_code))
                else:
                    self.binary_debug.write('  %-14s = %i\n' % ('format_code', self.format_code))
                self.binary_debug.write('    sort_bits[0] = %i -> is_random=%s vs mag/phase\n' % (self.sort_bits[0], self.is_random()))
                self.binary_debug.write('    sort_bits[1] = %i -> is_sort1 =%s vs sort2\n' % (self.sort_bits[1], self.is_sort1()))
                self.binary_debug.write('    sort_bits[2] = %i -> is_real  =%s vs real/imag\n' % (self.sort_bits[2], self.is_real()))
                sort_method, is_real, is_random = self._table_specs()
                self.binary_debug.write('    sort_method = %s\n' % sort_method)

                if self.is_complex():
                    self.binary_debug.write('\n  %-14s = %i -> is_mag_phase vs is_real_imag vs. is_random\n' % ('format_code', self.format_code))
                else:
                    self.binary_debug.write('  %-14s = %i\n' % ('format_code', self.format_code))

            self.binary_debug.write('  recordi = [%s]\n\n' % msg)

    def _read_geom_4(self, mapper, data, ndata):
        if self.read_mode == 1:
            return ndata
        if not self.make_geom:
            return ndata
        n = 0
        keys = unpack(b'3i', data[n:n+12])
        n += 12
        if len(data) == 12:
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

        if 0:
            # we're only going to use the keys if istream=0 (so the beginning of the record)
            if self.istream == 0 and keys in mapper:
                pass
            elif self.isubtable_old == self.isubtable:
                # we didn't increment the record, so we fix the n+=12 statement we called before
                # then we toss the keys and use the old geom_keys
                n = 0
                keys = self.geom_keys
            else:
                msg = 'keys=%s not found - %s; istream=%s; isubtable=%s isubtable_old=%s\n mapper=%s' % (
                    str(keys), self.table_name, self.istream, self.isubtable, self.isubtable_old,
                    mapper.keys())
                raise NotImplementedError(msg)

        try:
            name, func = mapper[keys]
        except KeyError:
            return n
        if self.is_debug_file:
            self.binary_debug.write('  found keys=%s -> name=%-6s - %s\n' % (str(keys), name, self.table_name))
        if self.debug:
            self.log.debug("  found keys=(%5s,%4s,%4s) name=%-6s - %s" % (keys[0], keys[1], keys[2], name, self.table_name))

        n = func(data, n)  # gets all the grid/mat cards
        assert n is not None, name

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

    def _read_table(self, data, ndata, result_name, storage_obj,
                    real_obj, complex_obj,
                    real_vector, complex_vector,
                    node_elem, random_code=None, is_cid=False):
        raise RuntimeError('is this still used?')
        assert isinstance(result_name, string_types), 'result_name=%r' % result_name
        assert isinstance(storage_obj, dict), 'storage_obj=%r' % storage_obj
        #assert real_obj is None
        #assert complex_obj is None
        #assert thermal_real_obj is None

        #print('self.num_wide =', self.num_wide)
        #print('random...%s' % self.isRandomResponse())
        #if not self.isRandomResponse():
        if self.format_code == 1 and self.num_wide == 8:  # real/random
            # real_obj
            assert real_obj is not None
            nnodes = ndata // 32  # 8*4
            auto_return, is_vectorized = self._create_table_object(result_name, nnodes, storage_obj, real_obj, real_vector, is_cid=is_cid)
            if auto_return:
                return ndata

            self._fix_format_code(format_code=1)
            if self.is_sort1():
                if self.nonlinear_factor is None:
                    n = self._read_real_table_static(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                else:
                    n = self._read_real_table_sort1(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
            else:
                n = self._read_real_table_sort2(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                #n = ndata
                #msg = self.code_information()
                #n = self._not_implemented_or_skip(data, ndata, msg)
        elif self.format_code in [2, 3] and self.num_wide == 14:  # real or real/imaginary or mag/phase
            # complex_obj
            assert complex_obj is not None
            nnodes = ndata // 56  # 14*4
            if self.is_debug_file:
                self.binary_debug.write('nnodes=%s' % nnodes)
            auto_return, is_vectorized = self._create_table_object(result_name, nnodes, storage_obj, complex_obj, complex_vector)
            if auto_return:
                return ndata
            if self.is_sort1():
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
        #else:
        #msg = 'invalid random_code=%s num_wide=%s' % (random_code, self.num_wide)
        #n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_table_vectorized(self, data, ndata, result_name, storage_obj,
                               real_vector, complex_vector,
                               node_elem, random_code=None, is_cid=False):

        assert isinstance(result_name, string_types), 'result_name=%r' % result_name
        assert isinstance(storage_obj, dict), 'storage_obj=%r' % storage_obj
        #print('self.num_wide =', self.num_wide)
        #print('random...%s' % self.isRandomResponse())
        #if not self.isRandomResponse():
        is_vectorized = True
        if self.format_code == 1 and self.num_wide == 8:  # real/random
            # real
            nnodes = ndata // 32  # 8*4
            auto_return = self._create_table_vector(
                result_name, nnodes, storage_obj, real_vector, is_cid=is_cid)
            if auto_return:
                return ndata

            self._fix_format_code(format_code=1)
            if self.is_sort1():
                if self.nonlinear_factor is None:
                    n = self._read_real_table_static(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                else:
                    n = self._read_real_table_sort1(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
            else:
                n = self._read_real_table_sort2(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                #n = ndata
                #msg = self.code_information()
                #n = self._not_implemented_or_skip(data, ndata, msg)
        elif self.format_code in [2, 3] and self.num_wide == 14:  # real or real/imaginary or mag/phase
            # complex
            nnodes = ndata // 56  # 14*4
            if self.is_debug_file:
                self.binary_debug.write('nnodes=%s' % nnodes)
            auto_return = self._create_table_vector(
                result_name, nnodes, storage_obj, complex_vector)
            if auto_return:
                return ndata
            if self.is_sort1():
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
        #else:
        #msg = 'invalid random_code=%s num_wide=%s' % (random_code, self.num_wide)
        #n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_scalar_table_vectorized(self, data, ndata, result_name, storage_obj,
                               real_vector, complex_vector,
                               node_elem, random_code=None, is_cid=False):

        assert isinstance(result_name, string_types), 'result_name=%r' % result_name
        assert isinstance(storage_obj, dict), 'storage_obj=%r' % storage_obj
        #print('self.num_wide =', self.num_wide)
        #print('random...%s' % self.isRandomResponse())
        #if not self.isRandomResponse():
        is_vectorized = True
        if self.format_code == 1 and self.num_wide == 8:  # real/random
            # real
            nnodes = ndata // 32  # 8*4
            auto_return = self._create_table_vector(
                result_name, nnodes, storage_obj, real_vector, is_cid=is_cid)
            if auto_return:
                return ndata

            self._fix_format_code(format_code=1)
            if self.is_sort1():
                if self.nonlinear_factor is None:
                    n = self._read_real_scalar_table_static(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                else:
                    n = self._read_real_scalar_table_sort1(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
            else:
                n = self._read_real_scalar_table_sort2(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                #n = ndata
                #msg = self.code_information()
                #n = self._not_implemented_or_skip(data, ndata, msg)
        elif self.format_code in [2, 3] and self.num_wide == 14:  # real or real/imaginary or mag/phase
            raise NotImplementedError('real/imaginary or mag/phase')
            # complex
            nnodes = ndata // 56  # 14*4
            if self.is_debug_file:
                self.binary_debug.write('nnodes=%s' % nnodes)
            auto_return = self._create_table_vector(
                result_name, nnodes, storage_obj, complex_vector)
            if auto_return:
                return ndata
            if self.is_sort1():
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
        #else:
        #msg = 'invalid random_code=%s num_wide=%s' % (random_code, self.num_wide)
        #n = self._not_implemented_or_skip(data, ndata, msg)
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
                return 2
            return 1
        elif self._function_code == 2:
            return value % 100
        elif self._function_code == 3:
            return value % 1000
        elif self._function_code == 4:
            return value // 10
        elif self._function_code == 5:
            return value % 10
        #elif self._function_code == 6:
            #raise NotImplementedError(self.function_code)
        #elif self._function_code == 7:
            #raise NotImplementedError(self.function_code)
        raise NotImplementedError(self.function_code)

    def _read_real_scalar_table_static(self, data, is_vectorized, nnodes, result_name, flag, is_cid=False):
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

        if self.use_vector and is_vectorized:
            n = nnodes * 4 * 8
            itotal2 = obj.itotal + nnodes
            #print('ndata=%s n=%s nnodes=%s' % (ndata, n, nnodes))
            ints = fromstring(data, dtype=self.idtype).reshape(nnodes, 8)
            floats = fromstring(data, dtype=self.fdtype).reshape(nnodes, 8)
            obj._times[obj.itime] = dt
            #self.node_gridtype[self.itotal, :] = [node_id, grid_type]
            #self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]
            #obj.node_gridtype[obj.itotal:itotal2, :] = ints[:, 0:1]
            nids = ints[:, 0] // 10
            assert nids.min() > 0, nids.min()
            obj.node_gridtype[obj.itotal:itotal2, 0] = nids
            obj.node_gridtype[obj.itotal:itotal2, 1] = ints[:, 1]
            obj.data[obj.itime, obj.itotal:itotal2, 0] = floats[:, 2]
            assert np.abs(floats[:, 3:]).max() == 0, '%s is not a scalar result...' % obj.__class__.__name__
            obj.itotal = itotal2
        else:
            n = 0
            s = Struct(b(self._endian + '2i6f'))
            for inode in range(nnodes):
                out = s.unpack(data[n:n+32])
                eid_device, grid_type, tx = out[:3]
                eid = eid_device // 10
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
                #print(out)
                obj.add(eid, grid_type, tx)
                n += 32
        return n

    def _read_real_scalar_table_sort1(self, data, is_vectorized, nnodes, result_name, flag, is_cid=False):
        """
        With a real transient result (e.g. SOL 109/159), reads a
        real OUG-style table created by:
              DISP(PLOT,SORT1,REAL) = ALL
        """
        # print('result_name=%s use_vector=%s is_vectorized=%s' % (result_name, self.use_vector, is_vectorized))
        if self.is_debug_file:
            self.binary_debug.write('  _read_real_scalar_table_sort1\n')
        #assert flag in ['node', 'elem'], flag
        #assert self.obj is not None
        dt = self.nonlinear_factor
        obj = self.obj
        if self.use_vector and is_vectorized:
            itime = obj.itime
            n = nnodes * 4 * 8
            itotal = obj.itotal
            itotal2 = itotal + nnodes

            is_mask = False
            if obj.itime == 0:
                ints = fromstring(data, dtype=self.idtype).reshape(nnodes, 8)
                #ints = floats[:, :2].view('int32')
                #from numpy import array_equal
                #assert array_equal(ints, intsB)

                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[itotal:itotal2, 0] = nids
                obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1]
                #obj.Vn = ones(floats.shape, dtype=bool)
                #obj.Vn[itotal:itotal2, :2] = False
                #print(obj.Vn)
                if obj.nonlinear_factor is not None and is_mask:
                    raise NotImplementedError('masking1')
                    float_mask = np.arange(nnodes * 8, dtype=np.int32).reshape(nnodes, 8)[:, 2:]
                    obj.float_mask = float_mask

            if obj.nonlinear_factor is not None and is_mask:
                raise NotImplementedError('masking2')
                results = fromstring(data, dtype=self.fdtype)[obj.float_mask]
            else:
                floats = fromstring(data, dtype=self.fdtype).reshape(nnodes, 8)
                obj.data[obj.itime, obj.itotal:itotal2, 0] = floats[:, 2]
                assert np.abs(floats[:, 3:]).max() == 0, '%s is not a scalar result...' % obj.__class__.__name__
            obj._times[itime] = dt
            obj.itotal = itotal2
        else:
            n = 0
            assert nnodes > 0, nnodes
            s = Struct(b(self._endian + '2i6f'))
            for inode in range(nnodes):
                out = s.unpack(data[n:n+32])
                eid_device, grid_type, tx = out[:3]
                eid = eid_device // 10
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
                obj.add_sort1(dt, eid, grid_type, tx)
                n += 32
        return n

    def _read_real_scalar_table_sort2(self, data, is_vectorized, nnodes, result_name, flag, is_cid=False):
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
                ints = fromstring(data, dtype=self.idtype).reshape(nnodes, 8)
                #nids = ints[:, 0] // 10
                nids = ones(nnodes, dtype='int32') * eid
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[itotal:itotal2, 0] = nids
                obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1]

            floats = fromstring(data, dtype=self.fdtype).reshape(nnodes, 8)
            obj._times[itime] = floats[:, 0]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 2]
            assert np.abs(floats[:, 3:]).max() == 0, '%s is not a scalar result...' % obj.__class__.__name__
            obj.itotal = itotal2
        else:
            n = 0
            assert nnodes > 0

            flag = 'freq/dt/mode'
            s = Struct(b(self._endian + self._analysis_code_fmt + 'i6f'))
            assert eid > 0, self.code_information()
            for inode in range(nnodes):
                edata = data[n:n+32]
                out = s.unpack(edata)
                (dt, grid_type, tx) = out[:3]
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, dt, str(out)))
                obj.add_sort2(dt, eid, grid_type, tx)
                n += 32
        return n

    def _read_real_table_static(self, data, is_vectorized, nnodes, result_name, flag, is_cid=False):
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

        if self.use_vector and is_vectorized:
            n = nnodes * 4 * 8
            itotal2 = obj.itotal + nnodes
            #print('ndata=%s n=%s nnodes=%s' % (ndata, n, nnodes))
            ints = fromstring(data, dtype=self.idtype).reshape(nnodes, 8)
            floats = fromstring(data, dtype=self.fdtype).reshape(nnodes, 8)
            obj._times[obj.itime] = dt
            #self.node_gridtype[self.itotal, :] = [node_id, grid_type]
            #self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]
            #obj.node_gridtype[obj.itotal:itotal2, :] = ints[:, 0:1]
            nids = ints[:, 0] // 10
            assert nids.min() > 0, nids.min()
            obj.node_gridtype[obj.itotal:itotal2, 0] = nids
            obj.node_gridtype[obj.itotal:itotal2, 1] = ints[:, 1]
            obj.data[obj.itime, obj.itotal:itotal2, :] = floats[:, 2:]
            obj.itotal = itotal2
        else:
        #if 1:
            n = 0
            s = Struct(b(self._endian + '2i6f'))
            for inode in range(nnodes):
                out = s.unpack(data[n:n+32])
                (eid_device, grid_type, tx, ty, tz, rx, ry, rz) = out
                eid = eid_device // 10
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
                #print(out)
                obj.add(eid, grid_type, tx, ty, tz, rx, ry, rz)
                n += 32
        return n


    #@jit(str, str, str, boolean)
    #@jit('int64(str, str, str, boolean)')
    #@int_(str, str, str, boolean)
    #@autojit
    def _read_real_table_sort1(self, data, is_vectorized, nnodes, result_name, flag, is_cid=False):
        """
        With a real transient result (e.g. SOL 109/159), reads a
        real OUG-style table created by:
              DISP(PLOT,SORT1,REAL) = ALL
        """
        # print('result_name=%s use_vector=%s is_vectorized=%s' % (result_name, self.use_vector, is_vectorized))
        if self.is_debug_file:
            self.binary_debug.write('  _read_real_table_sort1\n')
        #assert flag in ['node', 'elem'], flag
        #assert self.obj is not None
        dt = self.nonlinear_factor
        obj = self.obj
        if self.use_vector and is_vectorized:
            itime = obj.itime
            n = nnodes * 4 * 8
            itotal = obj.itotal
            itotal2 = itotal + nnodes

            is_mask = False
            if obj.itime == 0:
                ints = fromstring(data, dtype=self.idtype).reshape(nnodes, 8)
                #ints = floats[:, :2].view('int32')
                #from numpy import array_equal
                #assert array_equal(ints, intsB)

                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[itotal:itotal2, 0] = nids
                obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1]
                #obj.Vn = ones(floats.shape, dtype=bool)
                #obj.Vn[itotal:itotal2, :2] = False
                #print(obj.Vn)
                if obj.nonlinear_factor is not None and is_mask:
                    float_mask = np.arange(nnodes * 8, dtype=np.int32).reshape(nnodes, 8)[:, 2:]
                    obj.float_mask = float_mask

            if obj.nonlinear_factor is not None and is_mask:
                results = fromstring(data, dtype=self.fdtype)[obj.float_mask]
            else:
                floats = fromstring(data, dtype=self.fdtype).reshape(nnodes, 8)
                obj.data[obj.itime, obj.itotal:itotal2, :] = floats[:, 2:]
            obj._times[itime] = dt
            obj.itotal = itotal2
        else:
            n = 0
            assert nnodes > 0, nnodes
            s = Struct(b(self._endian + '2i6f'))
            for inode in range(nnodes):
                out = s.unpack(data[n:n+32])
                (eid_device, grid_type, tx, ty, tz, rx, ry, rz) = out
                eid = eid_device // 10
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
                obj.add_sort1(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
                n += 32
        return n

    #@autojit
    def _read_real_table_sort2(self, data, is_vectorized, nnodes, result_name, flag, is_cid=False):
        """
        With a real transient result (e.g. SOL 109/159), reads a
        real OUG-style table created by:
              DISP(PLOT,SORT2,REAL) = ALL
        """
        if self.is_debug_file:
            self.binary_debug.write('  _read_real_table_sort2\n')
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
                ints = fromstring(data, dtype=self.idtype).reshape(nnodes, 8)
                #nids = ints[:, 0] // 10
                nids = ones(nnodes, dtype='int32') * eid
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[itotal:itotal2, 0] = nids
                obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1]

            floats = fromstring(data, dtype=self.fdtype).reshape(nnodes, 8)
            obj._times[itime] = floats[:, 0]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 2:]
            obj.itotal = itotal2
        else:
            n = 0
            assert nnodes > 0

            flag = 'freq/dt/mode'
            s = Struct(b(self._endian + self._analysis_code_fmt + 'i6f'))
            assert eid > 0, self.code_information()
            for inode in range(nnodes):
                edata = data[n:n+32]
                out = s.unpack(edata)
                (dt, grid_type, tx, ty, tz, rx, ry, rz) = out
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, dt, str(out)))
                obj.add_sort2(dt, eid, grid_type, tx, ty, tz, rx, ry, rz)
                n += 32
        return n

    def _read_complex_table_sort1_mag(self, data, is_vectorized, nnodes, result_name, flag):
        if self.is_debug_file:
            self.binary_debug.write('  _read_complex_table_sort1_mag\n')
        assert flag in ['node', 'elem'], flag
        dt = self.nonlinear_factor

        n = 0
        obj = self.obj
        s = Struct(b(self._endian + '2i12f'))

        if self.use_vector and is_vectorized:
            n = nnodes * 4 * 14
            itotal2 = obj.itotal + nnodes

            if obj.itime == 0:
                ints = fromstring(data, dtype=self.idtype).reshape(nnodes, 14)
                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[obj.itotal:itotal2, 0] = nids
                obj.node_gridtype[obj.itotal:itotal2, 1] = ints[:, 1]

            floats = fromstring(data, dtype=self.fdtype).reshape(nnodes, 14)
            mag = floats[:, 2:8]
            phase = floats[:, 8:]
            rtheta = radians(phase)
            real_imag = mag * (cos(rtheta) + 1.j * sin(rtheta))
            #assert mag.shape == phase.shape, 'mag.shape=%s phase.shape=%s' % (str(mag.shape), str(phase.shape))
            #abs(real_imag), angle(real_imag, deg=True)

            obj._times[obj.itime] = dt
            obj.data[obj.itime, obj.itotal:itotal2, :] = real_imag
            obj.itotal = itotal2
        else:
            for inode in range(nnodes):
                out = s.unpack(data[n:n+56])
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
                n += 56
        return n

    def _read_complex_table_sort1_imag(self, data, is_vectorized, nnodes, result_name, flag):
        if self.is_debug_file:
            self.binary_debug.write('  _read_complex_table_sort1_imag\n')
        #assert flag in ['node', 'elem'], flag
        dt = self.nonlinear_factor
        obj = self.obj

        if self.use_vector and is_vectorized:
            n = nnodes * 4 * 14
            itotal = obj.itotal
            itotal2 = itotal + nnodes

            if obj.itime == 0:
                ints = fromstring(data, dtype=self.idtype).reshape(nnodes, 14)
                #print(ints[:, :2])
                nids = ints[:, 0] // 10
                assert nids.min() > 0, nids.min()
                obj.node_gridtype[itotal:itotal2, 0] = nids
                obj.node_gridtype[itotal:itotal2, 1] = ints[:, 1]

            floats = fromstring(data, dtype=self.fdtype).reshape(nnodes, 14)
            real = floats[:, 2:8]
            imag = floats[:, 8:]
            #assert real.shape == imag.shape, 'real.shape=%s imag.shape=%s' % (str(real.shape), str(imag.shape))

            obj._times[obj.itime] = dt
            obj.data[obj.itime, itotal:itotal2, :] = real + 1.j * imag
            obj.itotal = itotal2
        else:
        #if 1:
            n = 0
            s = Struct(b(self._endian + '2i12f'))

            assert self.obj is not None
            assert nnodes > 0
            for inode in range(nnodes):
                out = s.unpack(data[n:n+56])
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
                n += 56
        return n

    def _check_id(self, eid_device, flag, bdf_name, out):
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
            self.binary_debug.write('  _read_complex_table\n')
        assert flag in ['node', 'elem'], flag
        flag, flag_type = self.get_oug2_flag()
        node_id = self.nonlinear_factor

        #ntotal = 56  # 14 * 4

        assert self.obj is not None
        assert nnodes > 0
        #assert ndata % ntotal == 0

        if self.use_vector and is_vectorized and 0:
            pass
        else:
            n = 0
            s = Struct(self._endian + self._analysis_code_fmt + 'i12f')
            binary_debug_fmt = '  %s=%s %%s\n' % (flag, flag_type)
            for inode in range(nnodes):
                edata = data[n:n+56]
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
                self.obj.add_sort2(freq, node_id, grid_type, tx, ty, tz, rx, ry, rz)
                n += 56
        return n

    def _read_complex_table_sort2_imag(self, data, is_vectorized, nnodes, result_name, flag):
        """
        With a complex result (e.g. SOL 103/108), reads a complex OUG-style
        table created by:
          DISP(PLOT,SORT2,PHASE) = ALL
          DISP(PLOT,SORT2,IMAG) = ALL
        """
        if self.is_debug_file:
            self.binary_debug.write('  _read_complex_table\n')
        assert flag in ['node', 'elem'], flag
        flag, flag_type = self.get_oug2_flag()
        node_id = self.nonlinear_factor

        if self.use_vector and is_vectorized and 0:
            pass
        else:
            n = 0
            #ntotal = 56  # 14 * 4
            s = Struct(self._endian + self._analysis_code_fmt + 'i12f')
            assert self.obj is not None
            assert nnodes > 0
            #assert ndata % ntotal == 0

            binary_debug_fmt = '  %s=%s %%s\n' % (flag, flag_type)

            for inode in range(nnodes):
                edata = data[n:n+56]
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
                self.obj.add_sort2(freq, node_id, grid_type, tx, ty, tz, rx, ry, rz)
                n += 56
        return n

    def create_transient_object(self, storage_obj, class_obj, is_cid=False, debug=False):
        """
        Creates a transient object (or None if the subcase should be skippied).

        Parameters
        ----------
        storageName : str
            the name of the dictionary to store the object in (e.g. 'displacements')
        class_obj : object()
            the class object to instantiate
        debug : bool
            developer debug

        .. python ::

            slot = self.displacements
            slot_vector = RealDisplacementArray
            self.create_transient_object(slot, slot_vector, is_cid=is_cid)

        .. note:: dt can also be load_step depending on the class
        """
        assert not isinstance(class_obj, string_types), 'class_obj=%r' % class_obj
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

        #self.log.debug('self.table_name=%s isubcase=%s subtitle=%r' % (self.table_name, self.isubcase, self.subtitle.strip()))
        self.data_code['table_name'] = self.table_name.decode(self.encoding)
        assert self.log is not None

        code = self._get_code()
        #print('code =', code)
        if hasattr(self, 'isubcase'):
            if self.code in storage_obj:
                self.obj = storage_obj[code]
                if self.nonlinear_factor is not None:
                    if self.obj.nonlinear_factor is None:
                        msg = 'The object is flipping from a static (e.g. preload)\n'
                        msg += 'result to a transient/frequency based results\n'
                        msg += '%s -> %s\n' % (self.obj.nonlinear_factor, self.nonlinear_factor)
                        msg += 'code = (subcase=%s, analysis_code=%s, sort=%s, count=%s, subtitle=%s)\n' % tuple(code)
                        msg += '%s\n' % str(self.obj)
                        msg += '\nIf this isnt correct, check if the data code was applied on the object'
                        raise MultipleSolutionNotImplementedError(msg)
                self.obj.update_data_code(copy.deepcopy(self.data_code))
            else:
                class_obj.is_cid = is_cid
                self.obj = class_obj(self.data_code, self.is_sort1(), self.isubcase, self.nonlinear_factor)
            storage_obj[code] = self.obj
        else:
            if code in storage_obj:
                self.obj = storage_obj[code]
            else:
                storage_obj[code] = self.obj

    def _get_code(self):
        code = self.isubcase
        #code = (self.isubcase, self.analysis_code, self._sort_method, self._count, self.subtitle)
        code = (self.isubcase, self.analysis_code, self._sort_method, self._count, self.subtitle)
        #print('%r' % self.subtitle)
        self.code = code
        #self.log.debug('code = %s' % str(self.code))
        return self.code

    def _not_implemented_or_skip(self, data, ndata, msg=''):
        """
        A simple pass loop for unsupported tables that can be hacked on
        to crash the program everywhere that uses it.
        """
        #msg = 'table_name=%s table_code=%s %s\n%s' % (
            #self.table_name, self.table_code, msg, self.code_information())
        #raise NotImplementedError(msg)
        if is_release:
            if msg != self._last_comment:
                #print(self.code_information())
                if self.read_mode == 2:
                    self.log.warning(msg)
                #self.log.error(self.code_information())
                #self.log.warning(self.code_information())
                self._last_comment = msg
            return ndata
        else:
            msg = 'table_name=%s table_code=%s %s\n%s' % (
                self.table_name, self.table_code, msg, self.code_information())
            raise NotImplementedError(msg)

    def parse_approach_code(self, data):
        (approach_code, tCode, int3, isubcase) = unpack(b(self._endian + '4i'), data[:16])
        self.approach_code = approach_code
        self.tCode = tCode
        self.int3 = int3
        self.data_code['is_msc'] = self.is_msc

        if not hasattr(self, 'subtable_name'):
            self.data_code['subtable_name'] = self.subtable_name

        self.data_code['table_name'] = self.table_name
        self.data_code['approach_code'] = approach_code

        #: the local subcase ID
        self.isubcase = isubcase
        self.data_code['isubcase'] = self.isubcase
        #self.subcases.add(self.isubcase)  # set notation

        #: the type of result being processed
        self.table_code = tCode % 1000
        self.data_code['table_code'] = self.table_code

        #: used to create sort_bits
        self.sort_code = tCode // 1000

        self.data_code['sort_code'] = self.sort_code

        #: what type of data was saved from the run; used to parse the
        #: approach_code and grid_device.  device_code defines what options
        #: inside a result, STRESS(PLOT,PRINT), are used.
        self.device_code = approach_code % 10

        self.data_code['device_code'] = self.device_code
        assert self.device_code in [0, 1, 2, 3, 4, 5, 6, 7], self.device_code

        #: what solution was run (e.g. Static/Transient/Modal)
        self.analysis_code = (approach_code - self.device_code) // 10
        self.data_code['analysis_code'] = self.analysis_code

        #print('parse_approach_code - approach_code=%s tCode=%s int3=%s isubcase=%s' % (approach_code, tCode, int3, isubcase))
        #print('                 so - analysis_code=%s device_code=%s table_code=%s sort_code=%s\n' % (self.analysis_code, self.device_code, self.table_code, self.sort_code))
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
        sort_method, is_real, is_random = self._table_specs()
        assert sort_method in [1, 2], sort_method
        return sort_method

    def is_real(self):
        sort_method, is_real, is_random = self._table_specs()
        return is_real

    def is_complex(self):
        sort_method, is_real, is_random = self._table_specs()
        return not is_real

    def is_random(self):
        sort_method, is_real, is_random = self._table_specs()
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
        # if self.debug and self.table_name in self.show_table3_map:
            # return True
        # return False

    def debug4(self):
        return self.is_debug_file
        # if self.debug and self.table_name in self.show_table4_map:
            # return True
        # return False

    #def is_stress(self):
        #if self.stress_bits[1] == 0:
            #return True
        #return False

    def is_curvature(self):
        if self.is_stress():
            curvature_flag = False
        else:
            # strain only
            curvature_flag = True if self.stress_bits[2] == 0 else False
        if self.s_code in [10, 11, 20, 27]:
            assert curvature_flag, curvature_flag
            return True
        assert not curvature_flag, curvature_flag
        return False

    def is_fiber_distance(self):
        return not self.is_curvature()

    def is_max_shear(self):
        return True if self.stress_bits[4] == 0 else False

    def is_von_mises(self):
        return not self.is_max_shear()

    def is_stress(self):
        return not self.is_strain()

    def is_strain(self):
        if self.stress_bits[1] == 1:
            return True
        return False

    def _create_table_object(self, result_name, nnodes,
                             slot, slot_object, slot_vector, is_cid=False):
        assert isinstance(result_name, string_types), result_name
        assert isinstance(slot, dict), slot
        auto_return = False
        #print('%s nnodes=%s' % (result_name, nnodes))
        is_vectorized = self.is_vectorized
        if is_vectorized and slot_vector is None:
            is_vectorized = False

        if is_vectorized:
            if self.read_mode == 1:
                self.create_transient_object(slot, slot_vector, is_cid=is_cid)
                self.result_names.add(result_name)
                self.obj._nnodes += nnodes
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                self.obj = slot[self.code]
                #self.obj.update_data_code(self.data_code)
                self.obj.build()
        else:  # not vectorized
            self.result_names.add(result_name)
            if self.read_mode == 1:
                auto_return = True
                return auto_return, is_vectorized
            # pass = 0/2
            self.create_transient_object(slot, slot_object)
        return auto_return, is_vectorized

    def _create_table_vector(self, result_name, nnodes,
                             slot, slot_vector, is_cid=False):
        assert isinstance(result_name, string_types), result_name
        assert isinstance(slot, dict), slot
        auto_return = False
        #print('%s nnodes=%s' % (result_name, nnodes))
        self.result_names.add(result_name)
        if self.read_mode == 1:
            self.create_transient_object(slot, slot_vector, is_cid=is_cid)
            self.result_names.add(result_name)
            self.obj._nnodes += nnodes
            auto_return = True
        elif self.read_mode == 2:
            self.code = self._get_code()
            self.obj = slot[self.code]
            #self.obj.update_data_code(self.data_code)
            self.obj.build()
        else:
            auto_return = True
        return auto_return

    def _create_oes_object3(self, nelements, result_name, slot, obj, obj_vector):
        """
        Creates the self.obj parameter based on if this is vectorized or not.

        Parameters
        ----------
        nelements :  int
            the number of elements to preallocate for vectorization
        result_name : str
            unused
        slot : dict[(int, int, str)=obj
            the self dictionary that will be filled with a
            non-vectorized result
        obj : OES
            a pointer to the non-vectorized class
        obj_vector : OESArray
            a pointer to the vectorized class

        Returns
        -------
        auto_return : bool
            a flag indicating a return n should be called
        is_vectorized : bool
            True/False

        Since that's confusing, let's say we have real CTETRA stress data.
        We're going to fill self.solidStress with the class
        RealSolidStress.  If it were vectorized, we'd fill
        self.ctetra_stress. with RealSolidStressArray.  So we call:

        if self._is_vectorized(RealSolidStressArray, self.ctetra_stress):
            if self._results.is_not_saved(result_vector_name):
                return ndata
        else:
            if self._results.is_not_saved(result_name):
                return ndata

        auto_return, is_vectorized = self._create_oes_object3(self, nelements,
                            'ctetra_stress', self.ctetra_stress,
                            RealSolidStress, RealSolidStressArray)
        if auto_return:
            return nelements * ntotal
        """
        auto_return = False
        is_vectorized = self._is_vectorized(obj_vector, slot)
        #print('is_vectorized=%s result_name=%r' % (is_vectorized, result_name))
        if is_vectorized:
            #print("vectorized...read_mode=%s...%s" % (self.read_mode, result_name))
            if self.read_mode == 1:
                self.create_transient_object(slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % self.obj.ntimes)
                self.result_names.add(result_name)
                #print('self.obj =', self.obj)
                self.obj.nelements += nelements
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                #self.log.info("code = %s" % str(self.code))

                # if this is failing, you probably set obj_vector to None...
                try:
                    self.obj = slot[self.code]
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (self.code, result_name)
                    msg += "There's probably an extra check for read_mode=1..."
                    self.log.error(msg)
                    raise
                #self.obj.update_data_code(self.data_code)
                self.obj.build()

        else:  # not vectorized
            self.code = self._get_code()
            self.result_names.add(result_name)
            #print("not vectorized...read_mode=%s...%s" % (self.read_mode, result_name))
            #self.log.info("code = %s" % str(self.code))
            if self.read_mode == 1:
                self.result_names.add(result_name)
                auto_return = True
            # pass = 0/2
            self.create_transient_object(slot, obj)

        if auto_return and self.read_mode == 2:
            raise RuntimeError('this should never happen...auto_return=True read_mode=2')
        return auto_return, is_vectorized

    def _create_oes_object4(self, nelements, result_name, slot, obj_vector):
        """same as _create_oes_object4 except it doesn't support unvectorized objects"""
        auto_return = False
        #is_vectorized = True
        is_vectorized = self._is_vectorized(obj_vector, slot)
        #print("vectorized...read_mode=%s...%s; %s" % (self.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if self.read_mode == 1:
                #print('oes-self.nonlinear_factor =', self.nonlinear_factor)
                #print(self.data_code)
                self.create_transient_object(slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % obj.ntimes)
                self.result_names.add(result_name)
                #print('self.obj =', self.obj)
                self.obj.nelements += nelements
                auto_return = True
            elif self.read_mode == 2:
                self.code = self._get_code()
                #self.log.info("code = %s" % str(self.code))
                #print("code = %s" % str(self.code))

                # if this is failing, you probably set obj_vector to None...
                try:
                    self.obj = slot[self.code]
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (self.code, result_name)
                    msg += "There's probably an extra check for read_mode=1...%s" % result_name
                    self.log.error(msg)
                    raise
                assert self.obj.table_name == self.table_name.decode('utf-8'), 'obj.table_name=%s table_name=%s' % (self.obj.table_name, self.table_name)

                #obj.update_data_code(self.data_code)
                self.obj.build()

            else:  # not vectorized
                auto_return = True
        else:
            auto_return = True
        return auto_return, is_vectorized

    def _set_structs(self):
        """defines common struct formats"""
        self.fdtype = npdtype(self._endian + 'f4')
        self.idtype = npdtype(self._endian + 'i4')
        #self.sdtype = npdtype(self._endian + '4s')
        self.struct_i = Struct(b(self._endian + 'i'))
        self.struct_8s = Struct(b(self._endian + '8s'))
        self.struct_2i = Struct(b(self._endian + 'ii'))
