#pylint: disable=C0301,C0111
from __future__ import print_function, unicode_literals
import copy
from itertools import count
from struct import pack
from six import text_type, binary_type, PY3, string_types
import numpy as np

from pyNastran import is_release
from pyNastran.utils import object_attributes, object_methods

#from pyNastran.utils import list_print
from pyNastran.op2.op2_interface.op2_codes import Op2Codes, get_sort_method_from_table_name
from pyNastran.op2.op2_interface.write_utils import write_table_header, export_to_hdf5


class BaseScalarObject(Op2Codes):
    """
    The base scalar class is used by:
     - RealEigenvalues
     - BucklingEigenvalues
     - ComplexEigenvalues
     - ScalarObject
    """
    def __init__(self):
        Op2Codes.__init__(self)
        self.is_built = False

        # init value
        #  None - static
        #  int/float - transient/freq/load step/eigenvalue
        self.nonlinear_factor = np.nan
        self.format_code = None
        self.sort_code = None
        self.table_code = None
        self.title = None
        self.subtitle = None
        self.label = None
        self.num_wide = None
        self.device_code = None
        self.table_name = None
        #self.ntimes = 0
        #self.ntotal = 0
        #assert isinstance(self.name, (text_type, binary_type)), 'name=%s type=%s' % (self.name, type(self.name))

    def object_attributes(self, mode='public', keys_to_skip=None):
        if keys_to_skip is None:
            keys_to_skip = []
        elif isinstance(keys_to_skip, string_types):
            keys_to_skip = [keys_to_skip]

        my_keys_to_skip = [
            'object_methods', 'object_attributes',
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def object_methods(self, mode='public', keys_to_skip=None):
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = []

        my_keys_to_skip = [
            'object_methods', 'object_attributes',
        ]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def __eq__(self, table):
        #raise NotImplementedError(str(self.get_stats()))
        return False

    def __ne__(self, table):
        return not self == table

    @property
    def class_name(self):
        return self.__class__.__name__

    def get_headers(self):  # pragma: no cover
        raise RuntimeError()

    def _get_stats_short(self):  # pragma: no cover
        raise NotImplementedError('_get_stats_short')

    def build_dataframe(self):  # pragma: no cover
        """creates a pandas dataframe"""
        print('build_dataframe is not implemented in %s' % self.__class__.__name__)

    def export_to_hdf5(self, group, log):
        """exports the object to HDF5 format"""
        export_to_hdf5(self, group, log)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient(header, page_stamp, page_num=page_num, f06_file=f06_file,
                                             is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        msg = 'write_f06 is not implemented in %s\n' % self.__class__.__name__
        f06_file.write(msg)
        print(msg[:-1])
        #raise NotImplementedError(msg)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f06_file=None,
                             is_mag_phase=False, is_sort1=True):
        msg = '_write_f06_transient is not implemented in %s\n' % self.__class__.__name__
        f06_file.write(msg)
        print(msg[:-1])
        #raise NotImplementedError(msg)
        return page_num

    def __repr__(self):
        return ''.join(self.get_stats())

    def get_stats(self, short=False):
        msg = 'get_stats is not implemented in %s\n' % self.__class__.__name__
        if not is_release:
            raise NotImplementedError(msg)
        return msg


class ScalarObject(BaseScalarObject):
    """
    Used by all vectorized objects including:
     - DisplacementArray
     - RodStressArray
    """
    def __init__(self, data_code, isubcase, apply_data_code=True):
        assert 'nonlinear_factor' in data_code, data_code

        # new terms...I think they're all valid...
        self.result_name = None
        self.approach_code = None
        self.analysis_code = None
        self.data = None
        self._times = None

        self.isubcase = None
        self.ogs = None
        self.pval_step = None
        self.name = None
        self.superelement_adaptivity_index = None
        self._count = None
        #--------------------------------

        BaseScalarObject.__init__(self)
        self.isubcase = isubcase

        #--------------------------------
        self.data_frame = None
        # the nonlinear factor; None=static; float=transient
        self.dt = None
        # the number of time steps
        self.ntimes = 0
        # the number of entries in a table for a single time step
        self.ntotal = 0
        # there are a few tables (e.g. GridPointForcesArray) that change
        # length from time=1 to time=2.  As such, we will track the
        # length of each time step
        self._ntotals = []

        self.load_as_h5 = False
        self.h5_file = None
        if 'load_as_h5' in data_code:
            self.load_as_h5 = data_code['load_as_h5']
            del data_code['load_as_h5']
        if 'h5_file' in data_code:
            self.h5_file = data_code['h5_file']
            del data_code['h5_file']

        self.data_code = copy.deepcopy(data_code)
        #self.table_name = self.data_code['table_name']
        # if data code isn't being applied and you don't have
        # parameters that were in data_code (e.g.
        # self.element_name/self.element_type), you need to define
        # your vectorized class as setting data_code
        #
        # see pyNastran.op2.oes_stressStrain.real.oes_bars for an example
        # it's really subtle...
        if apply_data_code:
            self.apply_data_code()
            self._set_data_members()
        #print(self.code_information())

    #def isImaginary(self):
        #return bool(self.sort_bits[1])

    def _get_stats_short(self):
        msg = []
        class_name = self.__class__.__name__
        if hasattr(self, 'data'):
            unused_data = self.data
            shape = [int(i) for i in self.data.shape]
            headers = self.get_headers()
            headers_str = str(', '.join(headers))
            msg.append('%s[%s]; %s; [%s]\n' % (
                class_name, self.isubcase, shape, headers_str))
        return msg

    def __eq__(self, table):
        self._eq_header(table)
        #raise NotImplementedError(str(self.get_stats()))
        return False

    def _eq_header(self, table):
        is_nan = (self.nonlinear_factor is not None and
                  np.isnan(self.nonlinear_factor) and
                  np.isnan(table.nonlinear_factor))
        if not is_nan:
            assert self.nonlinear_factor == table.nonlinear_factor, 'nonlinear_factor=%s table.nonlinear_factor=%s' % (self.nonlinear_factor, table.nonlinear_factor)
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (
            self.table_name, table.table_name)
        assert self.approach_code == table.approach_code

    def __getstate__(self):
        """we need to remove the saved functions"""
        state = self.__dict__.copy()
        if 'add' in state:
            del state['add']
        if 'add_new_eid' in state:
            del state['add_new_eid']
        if 'class_name' in state:
            del state['class_name']
        if 'nnodes_per_element' in state:
            del state['nnodes_per_element']
        if 'add_node' in state:
            del state['add_node']
        if 'add_eid' in state:
            del state['add_eid']
        if '_add' in state:
            del state['_add']
        #if '_add_new_eid_sort1' in state:
            #del state['_add_new_eid_sort1']
        #if '_add_new_node_sort1' in state:
            #del state['_add_new_node_sort1']
        if '_add_new_eid' in state:
            del state['_add_new_eid']
        if '_add_new_node' in state:
            del state['_add_new_node']
        if 'dataframe' in state:
            del state['dataframe']

        #for key, value in state.items():
            #if isinstance(value, (int, float, str, np.ndarray, list)) or value is None:
                #continue
            #print(' %s = %s' % (key, value))
        #print(state)
        return state

    def _get_result_group(self):
        """gets the h5 result group"""
        code = self._get_code()
        case_name = 'Subcase=%s' % str(code) # (self.isubcase)
        #print('self.h5_file =', self.h5_file)
        if case_name in self.h5_file:
            subcase_group = self.h5_file[case_name]
        else:
            subcase_group = self.h5_file.create_group(case_name)
        group = subcase_group.create_group(self.result_name)
        return group

    def _get_code(self):
        code = self.isubcase
        ogs = 0
        if hasattr(self, 'ogs'):
            ogs = self.ogs
        #if self.binary_debug:
            #self.binary_debug.write(self.code_information(include_time=True))

        code = (self.isubcase, self.analysis_code, self._sort_method(), self._count, ogs,
                self.superelement_adaptivity_index, self.pval_step)
        #code = (self.isubcase, self.analysis_code, self._sort_method, self._count,
                #self.superelement_adaptivity_index, self.table_name_str)
        #print('%r' % self.subtitle)
        #self.code = code
        #self.log.debug('code = %s' % str(self.code))
        return code

    #@property
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
    def dataframe(self):
        """alternate way to get the dataframe"""
        return self.data_frame

    def apply_data_code(self):
        if self.table_name is not None and self.table_name != self.data_code['table_name']:
            print(self.data_code)
            msg = 'old_table_name=%r new_table_name=%r' % (
                self.table_name, self.data_code['table_name'])
            raise RuntimeError(msg)
        for key, value in sorted(self.data_code.items()):
            if PY3 and isinstance(value, bytes):
                print("  key=%s value=%s; value is bytes" % (key, value))
            self.__setattr__(key, value)
            #print("  key=%s value=%s" %(key, value))
        #if self.table_name in [b'OES1X', b'OES1X1']:

    def get_data_code(self, prefix='  '):
        msg = ''
        if 'data_names' not in self.data_code:
            return ['']

        msg += '%ssort1\n' % prefix if self.is_sort1 else '%ssort2\n' % prefix
        for name in self.data_code['data_names']:
            if hasattr(self, name + 's'):
                vals = getattr(self, name + 's')
                name = name + 's'
                vals_array = np.array(vals)
            elif hasattr(self, name):
                vals = getattr(self, name)
                vals_array = np.array(vals)
            else:
                vals_array = '???'
            #msg.append('%s = [%s]\n' % (name, ', '.join(['%r' % val for val in vals])))
            msg += '%s%s = %s\n' % (prefix, name, vals_array)
        #print("***data_names =", self.data_names)
        return [msg]

    def get_unsteady_value(self):
        name = self.data_code['name']
        return self._get_var(name)

    def _get_var(self, name):
        return getattr(self, name)

    def _set_var(self, name, value):
        return self.__setattr__(name, value)

    def _start_data_member(self, var_name, value_name):
        if hasattr(self, var_name):
            return True
        elif hasattr(self, value_name):
            self._set_var(var_name, [])
            return True
        return False

    def _append_data_member(self, var_name, value_name):
        """
        this appends a data member to a variable that may or may not exist
        """
        has_list = self._start_data_member(var_name, value_name)
        if has_list:
            listA = self._get_var(var_name)
            if listA is not None:
                #print("has %s" % var_name)
                value = self._get_var(value_name)
                try:
                    n = len(listA)
                except:
                    print("listA = ", listA)
                    raise
                listA.append(value)
                assert len(listA) == n + 1

    def _set_data_members(self):
        if 'data_names' not in self.data_code:
            msg = ('No "transient" variable was set for %s ("data_names" '
                   'was not defined in self.data_code).\n' % self.table_name)
            raise NotImplementedError(msg + self.code_information())

        for name in self.data_code['data_names']:
            self._append_data_member(name + 's', name)

    def update_data_code(self, data_code):
        if not self.data_code or (data_code['nonlinear_factor'] != self.data_code['nonlinear_factor']):
            self.data_code = data_code
            self.apply_data_code()  # take all the parameters in data_code and make them attributes of the class
            self._set_data_members()  # set the transient variables
        #else:
            #print('NF_new=%r NF_old=%r' % (data_code['nonlinear_factor'], self.data_code['nonlinear_factor']))

    def print_data_members(self):
        """
        Prints out the "unique" vals of the case.

        Uses a provided list of data_code['data_names'] to set the values for
        each subcase.  Then populates a list of self.name+'s' (by using
        setattr) with the current value.  For example, if the variable name is
        'mode', we make self.modes.  Then to extract the values, we build a
        list of of the variables that were set like this and then loop over
        them to print their values.

        This way there is no dependency on one result type having ['mode'] and
        another result type having ['mode','eigr','eigi'].
        """
        key_vals = []
        for name in self.data_code['data_names']:
            vals = getattr(self, name + 's')
            key_vals.append(vals)
            #print("%ss = %s" % (name, vals))

        msg = ''
        for name in self.data_code['data_names']:
            msg += '%-10s ' % name
        msg += '\n'

        nmodes = len(key_vals[0])
        for i in range(nmodes):
            for vals in key_vals:
                msg += '%-10g ' % vals[i]
            msg += '\n'
        return msg + '\n'

    def recast_gridtype_as_string(self, grid_type):
        """
        converts a grid_type integer to a string

        Point type (per NX 10; OUG table; p.5-663):
        =1, GRID Point
        =2, Scalar Point
        =3, Extra Point
        =4, Modal
        =5, p-elements, 0-DOF
        -6, p-elements, number of DOF
        """
        if grid_type == 1:
            grid_type_str = 'G'  # GRID
        elif grid_type == 2:
            grid_type_str = 'S'  # SPOINT
        elif grid_type == 7:
            grid_type_str = 'L'  # RIGID POINT (e.g. RBE3)
        elif grid_type == 0:
            grid_type_str = 'H'  # SECTOR/HARMONIC/RING POINT
        else:
            raise RuntimeError('grid_type=%s' % grid_type)
        return grid_type_str

    def cast_grid_type(self, grid_type_str):
        """converts a grid_type string to an integer"""
        if grid_type_str == 'G':
            grid_type = 1  # GRID
        elif grid_type_str == 'S':
            grid_type = 2  # SPOINT
        elif grid_type_str == 'L':
            grid_type = 7  # RIGID POINT (e.g. RBE3)
        elif grid_type_str == 'H':
            grid_type = 0  # SECTOR/HARMONIC/RING POINT
        else:
            raise RuntimeError('grid_type=%r' % grid_type_str)
        return grid_type

    def update_dt(self, data_code, unused_dt):
        """
        This method is called if the object already exits and a new
        time/freq/load step is found
        """
        self.data_code = data_code
        self.apply_data_code()
        msg = 'update_dt not implemented in the %s class' % self.__class__.__name__
        raise RuntimeError(msg)
        #assert dt>=0.
        #print("updating dt...dt=%s" % dt)
        #if dt is not None:
            #self.dt = dt
            #self.add_new_transient()

    def _build_dataframe_transient_header(self):
        """builds the header for the Pandas DataFrame/table"""
        assert isinstance(self.name, (text_type, binary_type)), 'name=%s type=%s' % (self.name, type(self.name))
        #print('self.name = %r' % self.name)
        #name = self.name #data_code['name']
        times = self._times
        utimes = np.unique(times)
        if not len(times) == len(utimes):
            msg = 'WARNING : %s - times=%s unique_times=%s...assuming new values...new_times=%s' % (
                self.__class__.__name__, times, utimes, np.arange(len(times)))
            print(msg)
            #raise RuntimeError(msg)
            times = np.arange(len(times))
        column_names = []
        column_values = []

        data_names = self.data_code['data_names']
        for name in data_names:
            #if name == primary_name:
            #times = self.da
            times = np.array(getattr(self, name + 's'))
            if name == 'mode':
                column_names.append('Mode')
                column_values.append(times)

                #if freq not in data_names:
                #if name == 'freq':
                ##if hasattr(self, 'freqs'):
                    #column_names.append('Freq')
                    #column_values.append(self.freqs)
                #elif name == 'eigr':
                    #column_names.append('eigenvalue_real')
                    #column_values.append(self.eigrs)
                #elif hasattr(self, 'eigrs') and 0:
                    #try:
                        #abs_freqs = np.sqrt(np.abs(self.eigrs)) / (2 * np.pi)
                    #except FloatingPointError:
                        #msg = 'Cant analyze freq = sqrt(eig)/(2*pi)\neigr=%s\n' % (self.eigrs)
                        #abs_freqs = np.sqrt(np.abs(self.eigrs)) / (2 * np.pi)
                        #msg += 'freq = sqrt(abs(self.eigrs)) / (2 * np.pi)=%s' % abs_freqs
                        #raise FloatingPointError(msg)
                    #column_names.append('Freq')
                    #column_values.append(abs_freqs)
                #else:
                    #pass

                # Convert eigenvalues to frequencies
                # TODO: add damping header
            elif name in ['eign']:
                abs_freqs = np.sqrt(np.abs(self.eigns)) / (2 * np.pi)
                column_names.append('Freq')
                column_values.append(abs_freqs)
                column_names.append('Eigenvalue')
                column_values.append(times)
                column_names.append('Radians')
                column_values.append(abs_freqs * 2 * np.pi)

            elif name in ['eigr']:
                column_names.append('EigenvalueReal')
                column_values.append(times)

            elif name in ['eigi']:
                column_names.append('EigenvalueImag')
                column_values.append(times)
                eigr = np.array(self.eigrs)
                eigi = np.array(self.eigis)
                denom = np.sqrt(eigr ** 2 + eigi ** 2)
                damping = np.zeros(len(eigr), dtype=eigr.dtype)
                inonzero = np.where(denom != 0)[0]
                if len(inonzero):
                    damping[inonzero] = -eigr[inonzero] / denom[inonzero]
                column_names.append('Damping')
                column_values.append(damping)
                #calculate_damping
            elif name in ['mode_cycle']:
                continue
                #column_names.append('mode_cycle(Freq?)')
                #column_values.append(times)
            elif name in ['mode2']:
                continue
                #column_names.append('mode2(Freq?)')
                #column_values.append(times)
            elif name in ['cycle']:
                continue
                #column_names.append('Freq (Cycles/s)')
                #column_values.append(times)

            elif name in ['freq', 'freq2']:
                column_names.append('Freq')
                column_values.append(times)
            elif name in ['dt', 'time']:
                column_names.append('Time')
                column_values.append(times)
            elif name in ['lftsfq', 'lsdvmn', 'load_step', 'loadID', 'loadFactor', 'loadIDs']:
                column_names.append('LoadStep')
                column_values.append(times)
            elif name == 'node_id':
                column_names.append('NodeID')
                column_values.append(times)
            elif name == 'element_id':
                column_names.append('EleemntID')
                column_values.append(times)
            else:
                msg = 'build_dataframe; name=%r' % name
                print(msg)
                raise NotImplementedError(msg)
        assert len(column_names) > 0, column_names
        assert len(column_names) == len(column_values), 'names=%s values=%s' % (column_names, column_values)
        assert len(self.get_headers()) == self.data.shape[-1], 'headers=%s; n=%s\ndata.headers=%s' % (self.get_headers(), len(self.get_headers()), self.data.shape[-1])
        return column_names, column_values

    def _write_table_header(self, op2_file, fascii, date):
        table_name = '%-8s' % self.table_name # 'BOUGV1  '
        fascii.write('%s._write_table_header\n' % table_name)
        #get_nmarkers- [4, 0, 4]
        #marker = [4, 2, 4]
        #table_header = [8, 'BOUGV1  ', 8]
        write_table_header(op2_file, fascii, table_name)


        #read_markers -> [4, -1, 4]
        #get_nmarkers- [4, 0, 4]
        #read_record - marker = [4, 7, 4]
        #read_record - record = [28, recordi, 28]

        #write_markers(op2_file, fascii, '  %s header1a' % self.table_name, [-1, 0, 7])
        data_a = [4, -1, 4,]
        #data_a = []
        #data_b = [4, -1, 4,]
        data_c = [4, 7, 4,]
        data = data_a + data_c
        blank = ' ' * len(self.table_name)
        fascii.write('%s header1a_i = %s\n' % (self.table_name, data_a))
        #fascii.write('%s            = %s\n' % (blank, data_b))
        fascii.write('%s            = %s\n' % (blank, data_c))
        op2_file.write(pack('<6i', *data))

        table1_fmt = b'<9i'
        table1 = [
            28,
            102, 0, 0, 0, 512, 0, 0,
            28,
        ]
        fascii.write('%s header1b = %s\n' % (self.table_name, table1))
        op2_file.write(pack(table1_fmt, *table1))

        #recordi = [subtable_name, month, day, year, 0, 1]

        data = [
            4, -2, 4,
            4, 1, 4,
            4, 0, 4,
            4, 7, 4,
        ]
        fascii.write('%s header2a = %s\n' % (self.table_name, data))
        op2_file.write(pack(b'<12i', *data))

        month, day, year = date
        try:
            unused_subtable_name = self.subtable_name
        except AttributeError:
            #print('attrs =', self.object_attributes())
            #raise
            pass

        self.subtable_name = b'OUG1    '
        table2 = [
            28,  # 4i -> 13i
            # subtable,todays date 3/6/2014, 0, 1  ( year=year-2000)
            b'%-8s' % self.subtable_name, month, day, year - 2000, 0, 1,
            28,
            ]
        table2_format = 'i8s6i'
        fascii.write('%s header2b = %s\n' % (self.table_name, table2))
        op2_file.write(pack(table2_format, *table2))


class BaseElement(ScalarObject):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        #--------------------------------
        # TODO: remove ???
        #self.element_name = None
        #self.element = None
        #self.element_node = None
        #self.element_type = None
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=apply_data_code)

    def _eq_header(self, table):
        ScalarObject._eq_header(self, table)
        is_nan = (self.nonlinear_factor is not None and
                  np.isnan(self.nonlinear_factor) and
                  np.isnan(table.nonlinear_factor))
        if hasattr(self, 'element_name'):
            if self.nonlinear_factor not in (None, np.nan) and not is_nan:
                assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                    self.element_name, self.element_type, self._times, table._times)

        if hasattr(self, 'element') and not np.array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'shape=%s element.shape=%s' % (
                self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\nEid\n' % str(self.code_information())
            for eid1, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid1, eid2)
            print(msg)
            raise ValueError(msg)

        if hasattr(self, 'element_node') and not np.array_equal(self.element_node, table.element_node):
            if self.element_node.shape != table.element_node.shape:
                msg = 'self.element_node.shape=%s table.element_node.shape=%s' % (
                    self.element_node.shape, table.element_node.shape)

                print(self.element_node.tolist())
                print(table.element_node.tolist())
                raise ValueError(msg)

            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            for i, (eid1, nid1), (eid2, nid2) in zip(count(), self.element_node, table.element_node):
                msg += '%s : (%s, %s), (%s, %s)\n' % (i, eid1, nid1, eid2, nid2)
            print(msg)
            raise ValueError(msg)
