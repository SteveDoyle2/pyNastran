from __future__ import print_function
import copy
from struct import Struct, unpack
from six import string_types


class GenericTableReader(object):
    def __init__(self, op2_reader, op2):
        self.op2 = op2
        self.op2_reader = op2_reader

    @property
    def idtype(self):
        return self.op2.idtype
    @property
    def fdtype(self):
        return self.op2.fdtype

    @property
    def thermal(self):
        """interface to the op2 object"""
        return self.op2_reader.thermal
    @property
    def sort_method(self):
        """interface to the op2 object"""
        return self.op2_reader.sort_method
    @property
    def nonlinear_factor(self):
        """interface to the op2 object"""
        return self.op2_reader.nonlinear_factor
    @property
    def use_vector(self):
        """interface to the op2 object"""
        return self.op2_reader.use_vector

    @property
    def obj(self):
        """interface to the op2 object"""
        return self.op2_reader.obj
    @obj.setter
    def obj(self):
        """interface to the op2 object"""
        raise RuntimeError()

    @property
    def log(self):
        """interface to the op2 object"""
        return self.op2.log

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

    #@property
    #def data_code(self):
        #"""interface to the op2 object"""
        #return self.op2.data_code

    def get_result(self, result_name):
        """interface to the op2 object"""
        return self.op2.get_result(result_name)

    #---------------------------------------------------------------------------
    def apply_data_code_value(self, name, value):
        self.op2_reader.data_code[name] = value

    def add_data_parameter(self, data, var_name, Type, field_num,
                           apply_nonlinear_factor=True, fix_device_code=False,
                           add_to_dict=True):
        op2_reader = self.op2_reader
        datai = data[4 * (field_num - 1) : 4 * (field_num)]
        assert len(datai) == 4, len(datai)
        #assert type(self._endian) == type(Type), 'endian=%r Type=%r' % (self._endian, Type)
        value, = unpack(self._endian + Type, datai)
        if fix_device_code:
            # was value = (value - self.device_code) // 10
            value = value // 10
        if self.is_debug_file:
            self.binary_debug.write('  %-14s = %r\n' % (var_name, value))
        #setattr(self, var_name, value)  # set the parameter to the local namespace

        if apply_nonlinear_factor:
            self.nonlinear_factor = value
            #if self.table_name == b'OUGV2':
                #assert isinstance(self.nonlinear_factor, int), self.nonlinear_factor
            op2_reader.data_code['nonlinear_factor'] = value
            op2_reader.data_code['name'] = var_name

        if add_to_dict:
            op2_reader.data_code[var_name] = value

        try:
            op2_reader.words[field_num - 1] = var_name
        except IndexError:
            msg = 'Trying to set word, but...len(words)=%s ifield=%s' % (len(self.words), field_num - 1)
            raise IndexError(msg)
        return value

    @property
    def is_sort1(self):
        #is_sort1_table = self.is_sort1
        op2_reader = self.op2_reader
        try:
            sort_method, is_real, is_random = op2_reader._table_specs()
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
        (approach_code, tCode, int3, isubcase) = unpack(self._endian + b'4i', data[:16])
        op2_reader = self.op2_reader
        op2_reader.approach_code = approach_code
        op2_reader.tCode = tCode
        op2_reader.int3 = int3

        data_code = op2_reader.data_code
        data_code['is_msc'] = self.op2.is_msc

        if not hasattr(op2_reader, 'subtable_name'):
            data_code['subtable_name'] = op2_reader.subtable_name


        data_code['table_name'] = op2_reader.table_name
        data_code['approach_code'] = approach_code

        #: the local subcase ID
        op2_reader.isubcase = isubcase
        data_code['isubcase'] = op2_reader.isubcase
        #op2.subcases.add(op2_reader.isubcase)  # set notation

        #: the type of result being processed
        op2_reader.table_code = tCode % 1000
        data_code['table_code'] = op2_reader.table_code
        data_code['tCode'] = op2_reader.tCode

        #: used to create sort_bits
        op2_reader.sort_code = tCode // 1000
        #Sort 1 - SortCode=((TCODE//1000)+2)//2

        data_code['sort_code'] = op2_reader.sort_code
        #print('tCode=%s tCode%%1000=%-2s tCode//1000=%s' % (tCode, tCode%1000, tCode//1000))
        op2_reader.sort_method = op2_reader._function1(tCode)
        data_code['sort_method'] = op2_reader.sort_method

        #: what type of data was saved from the run; used to parse the
        #: approach_code and grid_device.  device_code defines what options
        #: inside a result, STRESS(PLOT,PRINT), are used.
        op2_reader.device_code = approach_code % 10
        data_code['device_code'] = op2_reader.device_code
        assert op2_reader.device_code in [0, 1, 2, 3, 4, 5, 6, 7], op2_reader.device_code

        #: what solution was run (e.g. Static/Transient/Modal)
        op2_reader.analysis_code = (approach_code - op2_reader.device_code) // 10
        data_code['analysis_code'] = op2_reader.analysis_code

        #print('parse_approach_code - approach_code=%s tCode=%s int3=%s isubcase=%s' % (approach_code, tCode, int3, isubcase))
        #print('                 so - analysis_code=%s device_code=%s table_code=%s sort_code=%s\n' % (
            #op2_reader.analysis_code, op2_reader.device_code, op2_reader.table_code, op2_reader.sort_code))
        if 0:  # pragma: no cover
            if op2_reader.device_code == 3:
                #sys.stderr.write('The op2 may be inconsistent...\n')
                #sys.stderr.write("  print and plot can cause bad results..."
                #                 "if there's a crash, try plot only\n")
                op2_reader.device_code = 1

                #self.log.info('The op2 may be inconsistent...')
                #self.log.info('  print and plot can cause bad results...'
                #              'if there's a crash, try plot only')
                data_code['device_code'] = op2_reader.device_code

        if self.is_debug_file:
            self.binary_debug.write('  %-14s = %r\n' % ('table_name', op2_reader.table_name))
            self.binary_debug.write('  %-14s = analysis_code * 10 + device_code\n' % 'approach_code')
            self.binary_debug.write('  %-14s = %r\n' % ('approach_code', op2_reader.approach_code))
            self.binary_debug.write('  %-14s = %r\n' % ('  device_code', op2_reader.device_code))
            self.binary_debug.write('  %-14s = %r\n' % ('  analysis_code', op2_reader.analysis_code))
            self.binary_debug.write('  %-14s = sort_code * 1000 + table_code\n' % ('tCode'))
            self.binary_debug.write('  %-14s = %r\n' % ('tCode', op2_reader.tCode))
            self.binary_debug.write('  %-14s = %r\n' % ('  table_code', op2_reader.table_code))
            self.binary_debug.write('  %-14s = %r\n' % ('  sort_code', op2_reader.sort_code))
        self._parse_sort_code()
        assert op2_reader.sort_code in [0, 1, 2, 3, 4, 5, 6], op2_reader.sort_code #op2.code_information()

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
        op2_reader = self.op2_reader
        sort_code = op2_reader.sort_code

        # Sort codes can range from 0 to 7, but most of the examples
        # are covered by these.  The ones that break are incredibly large.
        if op2_reader.sort_code not in [0, 1, 2, 3, 4, 5, 6, 7]:
            msg = 'Invalid sort_code=%s' % (op2_reader.sort_code)
            raise SortCodeError(msg)
        i = 2
        while sort_code > 0:
            value = sort_code % 2
            sort_code = (sort_code - value) // 2
            bits[i] = value
            i -= 1

        # fixing bit[1]
        bits[1] = 0 if op2_reader.is_table_1 else 1
        #: the bytes describe the SORT information
        op2_reader.sort_bits = bits
        op2_reader.data_code['sort_bits'] = op2_reader.sort_bits

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
        op2_reader = self.op2_reader
        op2_reader.format_code_original = op2_reader.format_code
        if op2_reader.format_code == -1:
            if self.is_debug_file:
                self.write_ndata(self.binary_debug, 100)
            if op2_reader.table_name in [b'OESNLXR', b'OESNLBR', b'OESNLXD', b'OESNL1X']:
                assert op2_reader.format_code == -1, op2_reader.format_code
                op2_reader.format_code = 1
            else:
                raise RuntimeError(op2_reader.code_information())
            #return


        analysis_code = op2_reader.analysis_code
        random_code = op2_reader.random_code if hasattr(op2_reader, 'random_code') else 0
        format_code = op2_reader.format_code
        if random_code == 0:
            if analysis_code == 1:   # statics / displacement / heat flux
                assert format_code in [1, 3], op2_reader.code_information()
                format_code = 1
            elif analysis_code == 2:  # real eigenvalues
                assert format_code in [1, 3], op2_reader.code_information()
                format_code = 1
            #elif analysis_code == 3: # differential stiffness
            #elif analysis_code == 4: # differential stiffness
            elif analysis_code == 5:   # frequency
                assert format_code in [1, 2, 3], op2_reader.code_information()
                if format_code == 1:
                    format_code = 2
            elif analysis_code == 6:  # transient
                assert format_code in [1, 2, 3], op2_reader.code_information()
                format_code = 1
            elif analysis_code == 7:  # pre-buckling
                assert format_code in [1], op2_reader.code_information()
            elif analysis_code == 8:  # post-buckling
                assert format_code in [1, 2], op2_reader.code_information()
            elif analysis_code == 9:  # complex eigenvalues
                assert format_code in [1, 2, 3], op2_reader.code_information()
                if format_code == 1:
                    format_code = 2
            elif analysis_code == 10:  # nonlinear statics
                assert format_code in [1], op2_reader.code_information()
            elif analysis_code == 11:  # old geometric nonlinear statics
                assert format_code in [1], op2_reader.code_information()
            elif analysis_code == 12:
                # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
                assert format_code in [4], op2_reader.code_information() # invalid value
            else:
                msg = 'invalid analysis_code...analysis_code=%s' % analysis_code
                raise RuntimeError(msg)

            op2_reader.format_code = format_code
            op2_reader.data_code['format_code'] = format_code
        #assert self.format_code == 1, op2_reader.code_information()
        #if op2_reader.format_code != op2_reader.format_code_original:
            #print('op2_reader.format_code=%s orig=%s' % (
                #op2_reader.format_code, op2_reader.format_code_original))

    def _set_times_dtype(self):
        op2_reader = self.op2_reader
        op2_reader.data_code['_times_dtype'] = 'float32'
        analysis_code = op2_reader.analysis_code
        if analysis_code == 1:   # statics / displacement / heat flux
            pass # static doesn't have a type
        elif analysis_code == 2:  # real eigenvalues
            pass
        #elif self.analysis_code==3: # differential stiffness
        #elif self.analysis_code==4: # differential stiffness
        elif analysis_code == 5:   # frequency
            pass
        elif analysis_code == 6:  # transient
            pass
        elif analysis_code == 7:  # pre-buckling
            pass
        elif analysis_code == 8:  # post-buckling
            pass
        elif analysis_code == 9:  # complex eigenvalues
            pass
        elif analysis_code == 10:  # nonlinear statics
            pass
        elif analysis_code == 11:  # old geometric nonlinear statics
            pass
        elif analysis_code == 12:
            # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            pass
        else:
            msg = 'invalid analysis_code...analysis_code=%s' % analysis_code
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
        op2_reader = self.op2_reader
        bits = [0, 0, 0, 0, 0]
        thermal_code = op2_reader.thermal

        # Sort codes can range from 0 to 7, but most of the examples
        # are covered by these.  The ones that break are incredibly large.
        #if op2_reader.thermal not in [0, 1, 2, 3, 4, 5, 6, 7]:
            #msg = 'Invalid sort_code=%s' % (op2_reader.sort_code)
            #raise SortCodeError(msg)
            #if op2_reader.sort_code == 1145655:
                #return
        i = 4
        while thermal_code > 0:
            value = thermal_code % 2
            thermal_code = (thermal_code - value) // 2
            bits[i] = value
            i -= 1

        #: the bytes describe the Random information
        op2_reader.thermal_bits = bits
        op2_reader.data_code['thermal_bits'] = bits

    def is_mag_phase(self):
        return self.is_magnitude_phase()

    def is_magnitude_phase(self):
        if self.op2_reader.format_code == 3:
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
        op2_reader = self.op2_reader
        auto_return = False
        #is_vectorized = True
        is_vectorized = self._is_vectorized(obj_vector, slot)
        #print("vectorized...read_mode=%s...%s; %s" % (self.read_mode, result_name, is_vectorized))

        if is_vectorized:
            if self.read_mode == 1:
                #print('oes-self.nonlinear_factor =', self.nonlinear_factor)
                #print(self.data_code)
                op2_reader.create_transient_object(slot, obj_vector)
                #print("read_mode 1; ntimes=%s" % obj.ntimes)
                op2.result_names.add(result_name)
                #print('self.obj =', self.obj)
                op2_reader.obj.nelements += nelements
                auto_return = True
            elif self.read_mode == 2:
                op2_reader.code = op2_reader._get_code()
                #self.log.info("code = %s" % str(self.code))
                #print("code = %s" % str(self.code))

                # if this is failing, you probably set obj_vector to None...
                try:
                    op2_reader.obj = slot[op2_reader.code]
                except KeyError:
                    msg = 'Could not find key=%s in result=%r\n' % (op2_reader.code, result_name)
                    msg += "There's probably an extra check for read_mode=1...%s" % result_name
                    self.log.error(msg)
                    raise
                if not op2_reader.obj.table_name == op2_reader.table_name.decode('utf-8'):
                    msg = 'obj.table_name=%s table_name=%s' % (obj.table_name, op2_reader.table_name)
                    raise TypeError(msg)

                #obj.update_data_code(self.data_code)
                op2_reader.obj.build()
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



class TableReader(GenericTableReader):
    def __init__(self, op2_reader, op2):
        self.op2 = op2
        self.op2_reader = op2_reader
