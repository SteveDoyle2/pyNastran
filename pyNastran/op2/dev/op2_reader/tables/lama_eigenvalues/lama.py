"""
defines the LAMA class to read:
 - RealEigenvalues
 - ComplexEigenvalues
 - BucklingEigenvalues

from the OP2
"""
from __future__ import print_function
import copy
from struct import Struct, unpack
from six import string_types
from six.moves import range

#from pyNastran.op2.op2_interface.op2_common import OP2Common
from pyNastran.op2.tables.lama_eigenvalues.lama_objects import (
    RealEigenvalues, ComplexEigenvalues, BucklingEigenvalues)

class TableReader(object):
    def __init__(self, op2):
        self.op2 = op2

    @property
    def log(self):
        """interface to the op2 object"""
        return self.op2.log

    @property
    def encoding(self):
        """interface to the op2 object"""
        return self.op2.encoding

    @property
    def ntimes(self):
        """interface to the op2 object"""
        return self.op2.ntimes

    @property
    def _results(self):
        """interface to the op2 object"""
        return self.op2._results

    @property
    def is_vectorized(self):
        """interface to the op2 object"""
        return self.op2.is_vectorized

    @property
    def sort_bits(self):
        """interface to the op2 object"""
        return self.op2.sort_bits

    @property
    def table_name(self):
        """interface to the op2 object"""
        return self.op2.table_name

    @property
    def read_mode(self):
        """interface to the op2 object"""
        return self.op2.read_mode

    @property
    def is_debug_file(self):
        """interface to the op2 object"""
        return self.op2.is_debug_file

    @property
    def binary_debug(self):
        """interface to the op2 object"""
        return self.op2.binary_debug

    @property
    def _endian(self):
        """interface to the op2 object"""
        return self.op2._endian

    @property
    def data_code(self):
        """interface to the op2 object"""
        return self.op2.data_code

    def get_result(self, result_name):
        """interface to the op2 object"""
        return self.op2.get_result(result_name)

    @property
    def is_sort1(self):
        #is_sort1_table = self.is_sort1
        op2 = self.op2
        try:
            sort_method, is_real, is_random = op2._table_specs()
            return True if sort_method == 1 else False
        except AssertionError:
            table_name = self.table_name_str
            if table_name in SORT1_TABLES:
                is_sort1_table = True
            elif table_name in SORT2_TABLES:
                is_sort1_table = False
            else:
                try:
                    is_sort1_table = int(table_name[-1]) == 1
                except ValueError:
                    raise ValueError('is this SORT1/2?  table_name=%r' % table_name)
            return is_sort1_table

    def _read_title(self, data):
        op2 = self.op2
        op2._read_title_helper(data)

        if hasattr(op2, 'isubcase'):
            if op2.isubcase not in op2.isubcase_name_map:
                # 100 from label
                # 20 from subtitle line
                # 'SUBCASE 2'
                #op2.isubcase_name_map[isubcase] = [op2.subtitle, op2.label]
                op2.isubcase_name_map[op2.isubcase] = [
                    op2.subtitle, op2.superelement_adaptivity_index,
                    op2.analysis_code, op2.label]
        else:
            raise  RuntimeError('isubcase is not defined')

        if hasattr(op2, 'subtitle') and hasattr(op2, 'label'):
            ogs = 0
            if hasattr(op2, 'ogs'):
                ogs = op2.ogs
            unused_code = (op2.isubcase, op2.analysis_code,
                           op2.superelement_adaptivity_index,
                           op2.pval_step, ogs)
            #code = (op2.isubcase, op2.analysis_code,
                    #op2.superelement_adaptivity_index, op2.table_name_str)
            #print("code =", code)
            #if code not in op2.labels:
                #op2.subtitles[op2.isubcase].append(op2.subtitle)
                #op2.labels[code] = op2.label

    def parse_approach_code(self, data):
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
        op2 = self.op2
        (approach_code, tCode, int3, isubcase) = unpack(self._endian + b'4i', data[:16])
        op2.approach_code = approach_code
        op2.tCode = tCode
        op2.int3 = int3
        self.data_code['is_msc'] = op2.is_msc

        if not hasattr(op2, 'subtable_name'):
            self.data_code['subtable_name'] = op2.subtable_name

        self.data_code['table_name'] = op2.table_name
        self.data_code['approach_code'] = approach_code

        #: the local subcase ID
        op2.isubcase = isubcase
        self.data_code['isubcase'] = op2.isubcase
        #op2.subcases.add(op2.isubcase)  # set notation

        #: the type of result being processed
        op2.table_code = tCode % 1000
        self.data_code['table_code'] = op2.table_code
        self.data_code['tCode'] = op2.tCode

        #: used to create sort_bits
        op2.sort_code = tCode // 1000
        #Sort 1 - SortCode=((TCODE//1000)+2)//2

        self.data_code['sort_code'] = op2.sort_code
        #print('tCode=%s tCode%%1000=%-2s tCode//1000=%s' % (tCode, tCode%1000, tCode//1000))
        op2.sort_method = op2._function1(tCode)
        self.data_code['sort_method'] = op2.sort_method

        #: what type of data was saved from the run; used to parse the
        #: approach_code and grid_device.  device_code defines what options
        #: inside a result, STRESS(PLOT,PRINT), are used.
        op2.device_code = approach_code % 10

        self.data_code['device_code'] = op2.device_code
        assert op2.device_code in [0, 1, 2, 3, 4, 5, 6, 7], op2.device_code

        #: what solution was run (e.g. Static/Transient/Modal)
        op2.analysis_code = (approach_code - op2.device_code) // 10
        self.data_code['analysis_code'] = op2.analysis_code

        #print('parse_approach_code - approach_code=%s tCode=%s int3=%s isubcase=%s' % (approach_code, tCode, int3, isubcase))
        #print('                 so - analysis_code=%s device_code=%s table_code=%s sort_code=%s\n' % (op2.analysis_code, self.device_code, self.table_code, self.sort_code))
        if 0:  # pragma: no cover
            if op2.device_code == 3:
                #sys.stderr.write('The op2 may be inconsistent...\n')
                #sys.stderr.write("  print and plot can cause bad results..."
                #                 "if there's a crash, try plot only\n")
                op2.device_code = 1

                #self.log.info('The op2 may be inconsistent...')
                #self.log.info('  print and plot can cause bad results...'
                #              'if there's a crash, try plot only')
                self.data_code['device_code'] = op2.device_code

        if self.is_debug_file:
            self.binary_debug.write('  %-14s = %r\n' % ('table_name', op2.table_name))
            self.binary_debug.write('  %-14s = analysis_code * 10 + device_code\n' % 'approach_code')
            self.binary_debug.write('  %-14s = %r\n' % ('approach_code', op2.approach_code))
            self.binary_debug.write('  %-14s = %r\n' % ('  device_code', op2.device_code))
            self.binary_debug.write('  %-14s = %r\n' % ('  analysis_code', op2.analysis_code))
            self.binary_debug.write('  %-14s = sort_code * 1000 + table_code\n' % ('tCode'))
            self.binary_debug.write('  %-14s = %r\n' % ('tCode', op2.tCode))
            self.binary_debug.write('  %-14s = %r\n' % ('  table_code', op2.table_code))
            self.binary_debug.write('  %-14s = %r\n' % ('  sort_code', op2.sort_code))
        self._parse_sort_code()
        assert op2.sort_code in [0, 1, 2, 3, 4, 5, 6], op2.sort_code #op2.code_information()

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
        op2 = self.op2
        bits = [0, 0, 0]
        sort_code = op2.sort_code

        # Sort codes can range from 0 to 7, but most of the examples
        # are covered by these.  The ones that break are incredibly large.
        if op2.sort_code not in [0, 1, 2, 3, 4, 5, 6, 7]:
            msg = 'Invalid sort_code=%s' % (op2.sort_code)
            raise SortCodeError(msg)
        i = 2
        while sort_code > 0:
            value = sort_code % 2
            sort_code = (sort_code - value) // 2
            bits[i] = value
            i -= 1

        # fixing bit[1]
        bits[1] = 0 if op2.is_table_1 else 1
        #: the bytes describe the SORT information
        op2.sort_bits = bits
        self.data_code['sort_bits'] = op2.sort_bits

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

        .. note:: A case of 4 is not used and is used below as a placeholder,
                  while a case of -1 is some bizarre unhandled,
                  undocumented case.

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

    def is_mag_phase(self):
        return self.is_magnitude_phase()

    def is_magnitude_phase(self):
        if self.format_code == 3:
            return True
        return False

    def _create_oes_object4(self, nelements, result_name, slot, obj_vector):
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

        if self._is_vectorized(RealSolidStressArray, self.ctetra_stress):
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
        op2 = self.op2
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
                op2.result_names.add(result_name)
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
                if not self.obj.table_name == self.table_name.decode('utf-8'):
                    msg = 'obj.table_name=%s table_name=%s' % (self.obj.table_name, self.table_name)
                    raise TypeError(msg)

                #obj.update_data_code(self.data_code)
                self.obj.build()

            else:  # not vectorized
                auto_return = True
        else:
            auto_return = True

        assert is_vectorized, '%r is not vectorized; obj=%s' % (result_name, obj_vector)
        return auto_return, is_vectorized

    def _is_vectorized(self, obj_vector, slot_vector):
        """
        Checks to see if the data array has been vectorized

        Parameters
        ----------
        obj_vector:  the object to check
            (obj or None; None happens when vectorization hasn't been implemented)
        slot_vector: the dictionary to put the object in
            (dict or None; None happens when obj hasn't been implemented)

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
        op2 = self.op2
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

        #self.log.debug('self.table_name=%s isubcase=%s subtitle=%r' % (
            #self.table_name, self.isubcase, self.subtitle.strip()))
        self.data_code['table_name'] = self.table_name.decode(self.encoding)
        assert self.log is not None

        code = self._get_code()
        #print('code =', code)
        if hasattr(op2, 'isubcase'):
            if self.code in storage_obj:
                self.obj = storage_obj[code]
                if op2.nonlinear_factor is not None:
                    if self.obj.nonlinear_factor is None:
                        msg = 'The object is flipping from a static (e.g. preload)\n'
                        msg += 'result to a transient/frequency based results\n'
                        msg += '%s -> %s\n' % (self.obj.nonlinear_factor, self.nonlinear_factor)
                        msg += 'code = (subcase=%s, analysis_code=%s, sort=%s, count=%s, ogs=%s, superelement_adaptivity_index=%r pval_step=%r)\n' % tuple(code)
                        msg += '%s\n' % str(self.obj)
                        msg += '\nIf this isnt correct, check if the data code was applied on the object'
                        raise MultipleSolutionNotImplementedError(msg)
                self.obj.update_data_code(copy.deepcopy(self.data_code))
            else:
                print('making %s' % class_obj)
                class_obj.is_cid = is_cid
                is_sort1 = self.is_sort1  # uses the sort_bits

                self.obj = class_obj(self.data_code, is_sort1, op2.isubcase, op2.nonlinear_factor)
            storage_obj[code] = self.obj
        else:
            if code in storage_obj:
                self.obj = storage_obj[code]
            else:
                storage_obj[code] = self.obj

    def _get_code(self):
        op2 = self.op2
        code = op2.isubcase
        ogs = 0
        if hasattr(self, 'ogs'):
            ogs = self.ogs
        #if self.binary_debug:
            #self.binary_debug.write(self.code_information(include_time=True))

        code = (op2.isubcase, op2.analysis_code, op2._sort_method, op2._count, ogs,
                op2.superelement_adaptivity_index, op2.pval_step)
        #code = (self.isubcase, self.analysis_code, self._sort_method, self._count,
                #self.superelement_adaptivity_index, self.table_name_str)
        #print('%r' % self.subtitle)
        self.code = code
        #self.log.debug('code = %s' % str(self.code))
        return self.code


class LAMA(TableReader):
    def __init__(self, op2):
        TableReader.__init__(self, op2)

    #---------------------------------------------------------------------------
    def _read_complex_eigenvalue_3(self, data, ndata):
        """parses the Complex Eigenvalues Table 3 Data"""
        #raise NotImplementedError(self.table_name)
        op2 = self.op2
        op2.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        unused_three = self.parse_approach_code(data)
        op2.six = op2.add_data_parameter(data, 'six', b'i', 10, False)  # seven
        self._read_title(data)

    def _read_buckling_eigenvalue_3(self, data, ndata):
        """parses the Buckling Eigenvalues Table 3 Data"""
        #print(self.show_data(data))
        #self._read_title_helper(data)
        op2 = self.op2
        op2.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']
        #op2_reader.show_data(data)

        unused_three = op2.parse_approach_code(data)

        op2.seven = op2.add_data_parameter(data, 'seven', b'i', 10, False)  # seven
        #: residual vector augmentation flag
        op2.residual_flag = op2.add_data_parameter(data, 'residual_flag', b'i', 11, False)
        #: fluid modes flag
        op2.fluid_flag = op2.add_data_parameter(data, 'fluid_flag', b'i', 12, False)

        op2._read_title(data)

    def _read_complex_eigenvalue_4(self, data, ndata):
        """parses the Complex Eigenvalues Table 4 Data"""
        op2 = self.op2
        if op2.read_mode == 1:
            return ndata

        ntotal = 4 * 6
        nmodes = ndata // ntotal
        n = 0
        #assert op2.isubcase != 0, op2.isubcase
        clama = ComplexEigenvalues(op2.title, nmodes)
        op2.eigenvalues[op2.title] = clama
        #op2.eigenvalues[op2.isubcase] = lama
        structi = Struct(self._endian + b'ii4f')
        for i in range(nmodes):
            edata = data[n:n+ntotal]
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, order, eigr, eigc, freq, damping) = out # CLAMA
            #print('imode=%s order=%s eigr=%s eigc=%s freq=%s damping=%s' %
                  #(imode, order, eigr, eigc, freq, damping))
            clama.add_f06_line(out, i)
            n += ntotal
        assert n == ndata, 'clama length error'
        return n

    def _read_buckling_eigenvalue_4(self, data, ndata):
        """parses the Buckling Eigenvalues Table 4 Data"""
        # BLAMA - Buckling eigenvalue summary table
        # CLAMA - Complex eigenvalue summary table
        # LAMA - Normal modes eigenvalue summary table
        op2 = self.op2
        if op2.read_mode == 1:
            return ndata

        ntotal = 4 * 7
        nmodes = ndata // ntotal
        n = 0
        #assert op2.isubcase != 0, op2.isubcase
        blama = BucklingEigenvalues(op2.title, nmodes)
        op2.eigenvalues[op2.title] = blama
        #op2.eigenvalues[self.isubcase] = lama
        structi = Struct(self._endian + b'ii5f')
        for i in range(nmodes):
            edata = data[n:n+ntotal]
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, order, eigen, omega, freq, mass, stiff) = out # BLAMA??
            #(mode_num, extract_order, eigenvalue, radian, cycle, genM, genK) = line  # LAMA
            #(root_num, extract_order, eigr, eigi, cycle, damping) = data  # CLAMA
            blama.add_f06_line(out, i)
            n += ntotal
        return n

    def _read_real_eigenvalue_3(self, data, ndata):
        """parses the Real Eigenvalues Table 3 Data"""
        op2 = self.op2
        op2.words = [
            'aCode', 'tCode', '???', 'isubcase',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', '???', '???']
        #self.show_data(data)

        unused_three = op2.parse_approach_code(data)

        op2.seven = op2.add_data_parameter(data, 'seven', b'i', 10, False)  # seven
        ## residual vector augmentation flag
        op2.residual_flag = op2.add_data_parameter(data, 'residual_flag', b'i', 11, False)
        ## fluid modes flag
        op2.fluid_flag = op2.add_data_parameter(data, 'fluid_flag', b'i', 12, False)
        op2.title = None

        #print(op2.data_code)
        #op2.add_data_parameter(data,'format_code',  'i',9,False)   ## format code

        #: number of words per entry in record;
        #.. todo:: is this needed for this table ???
        #op2.add_data_parameter(data,'num_wide',     'i',10,False)

        #if self.analysis_code == 2: # sort2
            #op2.lsdvmn = op2.get_values(data,'i',5)

        #print("*isubcase=%s" % op2.isubcase)
        #print("analysis_code=%s table_code=%s thermal=%s" % (
            #op2.analysis_code, op2.table_code, op2.thermal))

        #op2.print_block(data)
        op2._read_title(data)

    def _read_real_eigenvalue_4(self, data, ndata):
        """parses the Real Eigenvalues Table 4 Data"""
        op2 = self.op2
        if self.read_mode == 1:
            return ndata
        #self.show_data(data)
        nmodes = ndata // 28
        n = 0
        ntotal = 28
        #assert op2.isubcase != 0, op2.isubcase
        lama = RealEigenvalues(op2.title, nmodes=nmodes)
        op2.eigenvalues[op2.title] = lama
        structi = Struct(self._endian + b'ii5f')
        for i in range(nmodes):
            edata = data[n:n+28]
            out = structi.unpack(edata)
            if self.is_debug_file:
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(imode, extract_order, eigenvalue, radian, cycle, gen_mass, gen_stiffness) = out
            lama.add_f06_line(out, i)
            n += ntotal
        return n
