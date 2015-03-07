#pylint: disable=C0301,C0111
from __future__ import print_function
from six import  iteritems
from six.moves import range
import copy

from pyNastran.op2.op2Codes import Op2Codes
#from pyNastran.utils import list_print


class BaseScalarObject(Op2Codes):
    def __init__(self):
        Op2Codes.__init__(self)
        self.is_built = False

    def name(self):
        return self.__class__.__name__

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        msg = 'write_f06 is not implemented in %s\n' % self.__class__.__name__
        f.write(msg)
        print(msg[:-1])
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        msg = '_write_f06_transient is not implemented in %s\n' % self.__class__.__name__
        f.write(msg)
        print(msg[:-1])
        return page_num

    def __repr__(self):
        return ''.join(self.get_stats())

    def get_stats(self):
        msg = 'get_stats is not implemented in %s\n' % self.__class__.__name__
        return msg


class ScalarObject(BaseScalarObject):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        assert 'nonlinear_factor' in data_code, data_code
        BaseScalarObject.__init__(self)
        self.isubcase = isubcase
        self.isTransient = False
        self.dt = None
        self.ntimes = 0
        self.ntotal = 0
        self.data_code = copy.deepcopy(data_code)

        #print("***data_code =", self.data_code)
        if apply_data_code:
            #print("adding data code...")
            self.apply_data_code()
            self.set_data_members()
        #self.log.debug(self.code_information())

    def isImaginary(self):
        return bool(self.sort_bits[1])

    def apply_data_code(self):
        #self.log = self.data_code['log']
        for key, value in sorted(iteritems(self.data_code)):
            #if key is not 'log':
                self.__setattr__(key, value)
                #self.log.debug("  key=%s value=%s" %(key, value))
                #print "  key=%s value=%s" %(key, value)
        #self.log.debug("")

    def get_data_code(self):
        msg = []
        if 'dataNames' not in self.data_code:
            return []

        for name in self.data_code['dataNames']:
            #try:
                if hasattr(self, name + 's'):
                    vals = getattr(self, name + 's')
                    name = name + 's'
                else:
                    vals = getattr(self, name)
                #msg.append('%s = %s\n' % (name, list_print(vals)))
                msg.append('%s = [%s]\n' % (name, ', '.join(['%r' % val for val in vals])) )
            #except AttributeError:  # weird case...
                #pass
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

        if 'element_name' in self.data_code:
            self.element_names = [self.data_code['element_name']]

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