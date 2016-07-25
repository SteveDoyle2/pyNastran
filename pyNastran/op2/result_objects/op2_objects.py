#pylint: disable=C0301,C0111
from __future__ import print_function, unicode_literals
from six import text_type, binary_type, iteritems, PY3, string_types
from six.moves import range
import copy
from struct import pack
import numpy as np

from pyNastran import is_release
from pyNastran.op2.op2_codes import Op2Codes
from pyNastran.utils import object_attributes, object_methods

#from pyNastran.utils import list_print
from pyNastran.op2.write_utils import write_table_header

class BaseScalarObject(Op2Codes):
    def __init__(self):
        Op2Codes.__init__(self)
        self.is_built = False

        # init value
        #  None - static
        #  int/float - transient/freq/load step/eigenvalue
        self.nonlinear_factor = None
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
        self._eq_header(table)
        return False

    def __ne__(self, table):
        return not self == table

    def _eq_header(self, table):
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (
            self.table_name, table.table_name)
        assert self.approach_code == table.approach_code

        if hasattr(self, 'element_name'):
            if self.nonlinear_factor is not None:
                assert np.array_equal(self._times, table._times), 'ename=%s-%s times=%s table.times=%s' % (
                    self.element_name, self.element_type, self._times, table._times)

        if hasattr(self, 'element'):
            if not np.array_equal(self.element, table.element):
                assert self.element.shape == table.element.shape, 'shape=%s element.shape=%s' % (
                    self.element.shape, table.element.shape)
                msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
                msg += '%s\nEid\n' % str(self.code_information())
                for eid1, eid2 in zip(self.element, table.element):
                    msg += '%s, %s\n' % (eid1, eid2)
                print(msg)
                raise ValueError(msg)

        if hasattr(self, 'element_node'):
            if not np.array_equal(self.element_node, table.element_node):
                assert self.element_node.shape == table.element_node.shape, 'shape=%s element_node.shape=%s' % (
                    self.element_node.shape, table.element_node.shape)
                msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
                msg += '%s\n' % str(self.code_information())
                for (eid1, nid1), (eid2, nid2) in zip(self.element_node, table.element_node):
                    msg += '(%s, %s), (%s, %s)\n' % (eid1, nid1, eid2, eid2)
                print(msg)
                raise ValueError(msg)

    @property
    def class_name(self):
        return self.__class__.__name__

    def get_headers(self):
        raise RuntimeError()

    def _get_stats_short(self):
        msg = []
        class_name = self.__class__.__name__
        if hasattr(self, 'data'):
            data = self.data
            shape = [int(i) for i in self.data.shape]
            headers = self.get_headers()
            headers_str = str(', '.join(headers))
            msg.append('%s[%s]; %s; [%s]\n' % (
            class_name, self.isubcase, shape, headers_str))
        return msg

    def _build_dataframe_transient_header(self):
        """builds the header for the Pandas DataFrame/table"""
        assert isinstance(self.name, (text_type, binary_type)), 'name=%s type=%s' % (self.name, type(self.name))
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
            elif name in ['eigr']:
                column_names.append('EigenvalueReal')
                column_values.append(times)
                abs_freqs = np.sqrt(np.abs(self.eigrs)) / (2 * np.pi)
                column_names.append('Freq')
                column_values.append(abs_freqs)
                column_names.append('Radians')
                column_values.append(abs_freqs * 2 * np.pi)

            elif name in ['eigi']:
                column_names.append('EigenvalueImag')
                column_values.append(times)
                eigr = np.array(self.eigrs)
                eigi = np.array(self.eigis)
                damping = -eigr / np.sqrt(eigr ** 2 + eigi ** 2)
                column_names.append('Damping (-eigr/sqrt(eigr^2+eigi^2); check)')
                column_values.append(times)
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
            elif name in ['dt','time']:
                column_names.append('Time')
                column_values.append(times)
            elif name in ['lftsfq', 'lsdvmn', 'load_step', 'loadID', 'loadFactor', 'loadIDs']:
                column_names.append('LoadStep')
                column_values.append(times)
            elif name == 'node_id':
                column_names.append('NodeID')
                column_values.append(times)
            else:
                msg = 'build_dataframe; name=%r' % name
                print(msg)
                raise NotImplementedError(msg)
        assert len(column_names) > 0, column_names
        assert len(column_names) == len(column_values), 'names=%s values=%s' % (column_names, column_values)
        assert len(self.get_headers()) == self.data.shape[-1], 'headers=%s; n=%s\ndata.headers=%s' % (self.get_headers(), len(self.get_headers()), self.data.shape[-1])
        return column_names, column_values

    def build_dataframe(self):
        print('build_dataframe is not implemented in %s' % self.__class__.__name__)

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num=page_num, f=f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        msg = 'write_f06 is not implemented in %s\n' % self.__class__.__name__
        f.write(msg)
        print(msg[:-1])
        #raise NotImplementedError(msg)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg = '_write_f06_transient is not implemented in %s\n' % self.__class__.__name__
        f.write(msg)
        print(msg[:-1])
        #raise NotImplementedError(msg)
        return page_num

    def __repr__(self):
        return ''.join(self.get_stats())

    def get_stats(self):
        msg = 'get_stats is not implemented in %s\n' % self.__class__.__name__
        if not is_release:
            raise NotImplementedError(msg)
        return msg


class ScalarObject(BaseScalarObject):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        assert 'nonlinear_factor' in data_code, data_code
        BaseScalarObject.__init__(self)
        self.isubcase = isubcase
        self.isTransient = False
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

        self.data_code = copy.deepcopy(data_code)

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

    def apply_data_code(self):
        if self.table_name is not None and self.table_name != self.data_code['table_name']:
            print(self.data_code)
            raise RuntimeError('old_table_name=%r new_table_name=%r' % (self.table_name, self.data_code['table_name']))
        for key, value in sorted(iteritems(self.data_code)):
            if PY3 and isinstance(value, bytes):
                print("  key=%s value=%s; value is bytes" % (key, value))
            self.__setattr__(key, value)
            #print("  key=%s value=%s" %(key, value))
        #if self.table_name in [b'OES1X', b'OES1X1']:

    def get_data_code(self):
        msg = ''
        if 'data_names' not in self.data_code:
            return ['']

        msg += 'sort1\n  ' if self.is_sort1() else 'sort2\n  '
        for name in self.data_code['data_names']:
            if hasattr(self, name + 's'):
                vals = getattr(self, name + 's')
                name = name + 's'
            else:
                vals = getattr(self, name)
            #msg.append('%s = [%s]\n' % (name, ', '.join(['%r' % val for val in vals])))
            msg += '%s = %s\n  ' % (name, np.array(vals))
        #print("***data_names =", self.data_names)
        return [msg.rstrip(' ')]

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
            msg = 'No "transient" variable was set for %s ("data_names" was not defined in self.data_code).\n' % self.table_name
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
        """converts a grid_type integer to a string"""
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

    def update_dt(self, data_code, dt):
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
        if dt is not None:
            self.dt = dt
            self.add_new_transient()

    def _write_table_header(self, f, fascii, date):
        table_name = '%-8s' % self.table_name # 'BOUGV1  '
        fascii.write('%s._write_table_header\n' % table_name)
        #get_nmarkers- [4, 0, 4]
        #marker = [4, 2, 4]
        #table_header = [8, 'BOUGV1  ', 8]
        write_table_header(f, fascii, table_name)


        #read_markers -> [4, -1, 4]
        #get_nmarkers- [4, 0, 4]
        #read_record - marker = [4, 7, 4]
        #read_record - record = [28, recordi, 28]

        #write_markers(f, fascii, '  %s header1a' % self.table_name, [-1, 0, 7])
        data_a = [4, -1, 4,]
        #data_a = []
        #data_b = [4, -1, 4,]
        data_c = [4, 7, 4,]
        data = data_a + data_c
        blank = ' ' * len(self.table_name)
        fascii.write('%s header1a_i = %s\n' % (self.table_name, data_a))
        #fascii.write('%s            = %s\n' % (blank, data_b))
        fascii.write('%s            = %s\n' % (blank, data_c))
        f.write(pack('<6i', *data))

        table1_fmt = b'<9i'
        table1 = [
            28,
            1, 2, 3, 4, 5, 6, 7,
            28,
        ]
        fascii.write('%s header1b = %s\n' % (self.table_name, table1))
        f.write(pack(table1_fmt, *table1))

        #recordi = [subtable_name, month, day, year, 0, 1]

        data = [
            4, -2, 4,
            4, 1, 4,
            4, 0, 4,
            4, 7, 4,
        ]
        fascii.write('%s header2a = %s\n' % (self.table_name, data))
        print('data =', data, len(data))
        f.write(pack(b'<12i', *data))

        month, day, year = date
        try:
            subtable_name = self.subtable_name
        except AttributeError:
            print('attrs =', self.object_attributes())

        self.subtable_name = b'OUG1    '
        table2 = [
            28,  # 4i -> 13i
            b'%-8s' % self.subtable_name, month, day, year - 2000, 0, 1,   # subtable,todays date 3/6/2014, 0, 1  ( year=year-2000)
            28,
            ]
        table2_format = 'i8s6i'
        fascii.write('%s header2b = %s\n' % (self.table_name, table2))
        f.write(pack(table2_format, *table2))
