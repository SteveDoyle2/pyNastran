#pylint: disable=C0301,C0111
from __future__ import print_function
from six import  iteritems
from six.moves import range
import copy
from struct import pack
from numpy import array, sqrt, pi

from pyNastran import is_release
from pyNastran.op2.op2Codes import Op2Codes
#from pyNastran.utils import list_print
#from pyNastran.op2.write_utils import write_table_header

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
        self.Title = None
        self.subtitle = None
        self.label = None
        self.num_wide = None
        self.device_code = None
        self.table_name = None

    def __eq__(self, table):
        return True

    def __ne__(self, table):
        return not self == table

    def name(self):
        return self.__class__.__name__

    def get_headers(self):
        raise RuntimeError()

    def _build_dataframe_transient_header(self):
        """builds the header for the Pandas DataFrame/table"""
        name = self.name #data_code['name']
        times = self._times
        column_names = []
        column_values = []
        skip_names = [] # 'Time'
        if name == 'mode':
            column_names.append('Mode')
            column_names.append('Freq')

            # Convert eigenvalues to frequencies
            freq  = sqrt(self.eigrs) / (2 * pi)
            column_values.append(times)
            column_values.append(freq)
        elif name == 'freq':
            column_names.append('Freq')
            column_values.append(times)
        elif name == 'loadID':
            column_names.append('LoadID')
            column_values.append(times)
        else:
            raise NotImplementedError('build_dataframe; name=%r' % name)
        assert len(column_names) > 0, column_names
        assert len(column_names) == len(column_values), column_names
        assert len(self.get_headers()) == self.data.shape[-1], 'headers=%s; n=%s\ndata.headers=%s' % (self.get_headers(), len(self.get_headers()), self.data.shape[-1])
        return column_names, column_values

    def build_dataframe(self):
        print('build_dataframe is not implemented in %s' % self.__class__.__name__)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
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
            self.set_data_members()
        #print(self.code_information())

    def isImaginary(self):
        return bool(self.sort_bits[1])

    def apply_data_code(self):
        for key, value in sorted(iteritems(self.data_code)):
            self.__setattr__(key, value)
            #print("  key=%s value=%s" %(key, value))
        #if self.table_name in [b'OES1X', b'OES1X1']:

    def get_data_code(self):
        msg = []
        if 'dataNames' not in self.data_code:
            return []

        msg.append('sort1\n  ' if self.is_sort1() else 'sort2\n  ')
        for name in self.data_code['dataNames']:
            if hasattr(self, name + 's'):
                vals = getattr(self, name + 's')
                name = name + 's'
            else:
                vals = getattr(self, name)
            #msg.append('%s = [%s]\n' % (name, ', '.join(['%r' % val for val in vals])))
            msg.append('%s = %s\n' % (name, array(vals)))
        #print("***dataNames =", self.dataNames)
        return msg

    def getUnsteadyValue(self):
        name = self.data_code['name']
        return self.getVar(name)

    def getVar(self, name):
        return getattr(self, name)

    def set_var(self, name, value):
        return self.__setattr__(name, value)

    def start_data_member(self, var_name, value_name):
        if hasattr(self, var_name):
            return True
        elif hasattr(self, value_name):
            self.set_var(var_name, [])
            return True
        return False

    def append_data_member(self, var_name, value_name):
        """
        this appends a data member to a variable that may or may not exist
        """
        hasList = self.start_data_member(var_name, value_name)
        if hasList:
            listA = self.getVar(var_name)
            if listA is not None:
                #print("has %s" % var_name)
                value = self.getVar(value_name)
                try:
                    n = len(listA)
                except:
                    print("listA = ", listA)
                    raise
                listA.append(value)
                assert len(listA) == n + 1

    def set_data_members(self):
        if 'dataNames' not in self.data_code:
            msg = 'No "transient" variable was set for %s ("dataNames" was not defined in self.data_code).\n' % self.table_name
            raise NotImplementedError(msg + self.code_information())

        for name in self.data_code['dataNames']:
            self.append_data_member(name + 's', name)

    def update_data_code(self, data_code):
        if not self.data_code or (data_code['nonlinear_factor'] != self.data_code['nonlinear_factor']):
            self.data_code = data_code
            self.apply_data_code()  # take all the parameters in data_code and make them attributes of the class
            self.set_data_members()  # set the transient variables
        #else:
            #print('NF_new=%r NF_old=%r' % (data_code['nonlinear_factor'], self.data_code['nonlinear_factor']))

    def print_data_members(self):
        """
        Prints out the "unique" vals of the case.
        Uses a provided list of data_code['dataNames'] to set the values for
        each subcase.  Then populates a list of self.name+'s' (by using
        setattr) with the current value.  For example, if the variable name is
        'mode', we make self.modes.  Then to extract the values, we build a
        list of of the variables that were set like this and then loop over
        them to print their values.

        This way there is no dependency on one result type having ['mode'] and
        another result type having ['mode','eigr','eigi'].
        """
        key_vals = []
        for name in self.data_code['dataNames']:
            vals = getattr(self, name + 's')
            key_vals.append(vals)
            #print("%ss = %s" % (name, vals))

        msg = ''
        for name in self.data_code['dataNames']:
            msg += '%-10s ' % name
        msg += '\n'

        nModes = len(key_vals[0])
        for i in range(nModes):
            for vals in key_vals:
                msg += '%-10g ' % vals[i]
            msg += '\n'
        return msg + '\n'

    def recast_gridtype_as_string(self, grid_type):
        """converts a grid_type integer to a string"""
        if grid_type == 1:
            Type = 'G'  # GRID
        elif grid_type == 2:
            Type = 'S'  # SPOINT
        elif grid_type == 7:
            Type = 'L'  # RIGID POINT (e.g. RBE3)
        elif grid_type == 0:
            Type = 'H'  # SECTOR/HARMONIC/RING POINT
        else:
            raise RuntimeError('grid_type=%s' % grid_type)
        return Type

    def cast_grid_type(self, grid_type):
        """converts a grid_type string to an integer"""
        if grid_type == 'G':
            Type = 1  # GRID
        elif grid_type == 'S':
            Type = 2  # SPOINT
        elif grid_type == 'L':
            Type = 7  # RIGID POINT (e.g. RBE3)
        elif grid_type == 'H':
            Type = 0  # SECTOR/HARMONIC/RING POINT
        else:
            raise RuntimeError('grid_type=%r' % grid_type)
        return Type

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

        table1_fmt = '<9i'
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
        f.write(pack('12i', *data))

        month, day, year = date
        table2 = [
            28,  # 4i -> 13i
            b'%-8s' % self.subtable_name, month, day, year - 2000, 0, 1,   # subtable,todays date 3/6/2014, 0, 1  ( year=year-2000)
            28,
            ]
        table2_format = 'i8s6i'
        fascii.write('%s header2b = %s\n' % (self.table_name, table2))
        f.write(pack(table2_format, *table2))
