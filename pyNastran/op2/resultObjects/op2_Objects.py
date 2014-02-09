#pylint: disable=C0301,C0111
from __future__ import print_function
import copy

from pyNastran.op2.op2Codes import Op2Codes
from pyNastran.utils import list_print


class baseScalarObject(Op2Codes):
    def __init__(self):
        Op2Codes.__init__(self)

    def name(self):
        return self.__class__.__name__

    def write_f06(self, header, page_stamp, pageNum=1, f=None, is_mag_phase=False):
        msg = 'write_f06 is not implemented in %s\n' % self.__class__.__name__
        f.write(msg)
        return pageNum

    def _write_f06_transient(self, header, page_stamp, pageNum=1, f=None, is_mag_phase=False):
        msg = '_write_f06_transient is not implemented in %s\n' % self.__class__.__name__
        f.write(msg)
        return pageNum

    def __repr__(self):
        return self.get_stats()

    def get_stats(self):
        msg = 'get_stats is not implemented in %s\n' % self.__class__.__name__
        return msg


class scalarObject(baseScalarObject):
    def __init__(self, data_code, isubcase):
        assert 'nonlinear_factor' in data_code, data_code
        baseScalarObject.__init__(self)
        self.isubcase = isubcase
        self.isTransient = False
        self.dt = None
        self.data_code = copy.deepcopy(data_code)

        if 0:
            if 'Title' in data_code:
                self.Title = data_code['Title']
                del data_code['Title']
            #else:
                #raise RuntimeError('Title doesnt exist in data_code=%s' % str(self.data_code))

            if 'subtitle' in data_code:
                self.subtitle = data_code['subtitle']
                del data_code['subtitle']
            #else:
                #raise RuntimeError('subtitle doesnt exist in data_code=%s' % str(self.data_code))


            if 'isubcase' in data_code:
                self.subtitle = data_code['isubcase']
                del data_code['isubcase']
            else:
                raise RuntimeError('isubcase doesnt exist in data_code=%s' % str(self.data_code))

            if 'label' in data_code:
                self.label = data_code['label']
                del data_code['label']
            #else:
                #raise RuntimeError('label doesnt exist in data_code=%s' % str(self.data_code))

            if 'table_name' in data_code:
                self.table_name = data_code['table_name']
                del data_code['table_name']
            else:
                raise RuntimeError('table_name doesnt exist in data_code=%s' % str(self.data_code))

            if 'table_code' in data_code:
                self.table_code = data_code['table_code']
                del data_code['table_code']
            else:
                raise RuntimeError('table_code doesnt exist in data_code=%s' % str(self.data_code))

            if 'aCode' in data_code:
                self.aCode = data_code['aCode']
                del data_code['aCode']
            if 'tCode' in data_code:
                self.tCode = data_code['tCode']
                del data_code['tCode']
            if 's_code' in data_code:
                self.s_code = data_code['s_code']
                del data_code['s_code']


            if 'mode_cycle' in data_code:
                self.mode_cycle = data_code['mode_cycle']
                del data_code['mode_cycle']

            if 'eigr' in data_code:
                self.eigr = data_code['eigr']
                del data_code['eigr']
            if 'eigi' in data_code:
                self.eigr = data_code['eigi']
                del data_code['eigi']

            if 'dataNames' in data_code:
                self.dataNames = data_code['dataNames']
                del data_code['dataNames']
            else:
                raise RuntimeError('dataNames doesnt exist in data_code=%s' % str(self.data_code))

            if 'nonlinear_factor' in data_code:
                self.nonlinear_factor = data_code['nonlinear_factor']
                del data_code['nonlinear_factor']
            else:
                raise RuntimeError('nonlinear_factor doesnt exist in data_code=%s' % str(self.data_code))

            if 'acousticFlag' in data_code:
                self.acousticFlag = data_code['acousticFlag']
                del data_code['acousticFlag']
            if 'randomCode' in data_code:
                self.randomCode = data_code['randomCode']
                del data_code['randomCode']
            if 'o_code' in data_code:
                self.o_code = data_code['o_code']
                del data_code['o_code']

            if 'num_wide' in data_code:
                self.num_wide = data_code['num_wide']
                del data_code['num_wide']
            else:
                raise RuntimeError('num_wide doesnt exist in data_code=%s' % str(self.data_code))

            if 'format_code' in data_code:
                self.format_code = data_code['format_code']
                del data_code['format_code']
            else:
                raise RuntimeError('format_code doesnt exist in data_code=%s' % str(self.data_code))

            if 'analysis_code' in data_code:
                self.analysis_code = data_code['analysis_code']
                del data_code['analysis_code']
            else:
                raise RuntimeError('analysis_code doesnt exist in data_code=%s' % str(self.data_code))

            if 'sort_code' in data_code:
                self.sort_code = data_code['sort_code']
                del data_code['sort_code']

            if 'sort_bits' in data_code:
                self.sort_bits = data_code['sort_bits']
                del data_code['sort_bits']

            if 'device_code' in data_code:
                self.device_code = data_code['device_code']
                del data_code['device_code']

            if 'thermal' in data_code:
                self.thermal = data_code['thermal']
                del data_code['thermal']
            if 'stress_bits' in data_code:
                self.stress_bits = data_code['stress_bits']
                del data_code['stress_bits']

            if 'element_name' in data_code:
                self.element_name = data_code['element_name']
                del data_code['element_name']
            if 'element_type' in data_code:
                self.element_type = data_code['element_type']
                del data_code['element_type']

            if 'time' in data_code:
                self.time = data_code['time']
                self.times = []
                self.name = 'time'
                del data_code['time']
                del data_code['name']
            elif 'dt' in data_code:
                self.dt = data_code['dt']
                self.dts = []
                self.name = 'dt'
                del data_code['dt']
                del data_code['name']
            elif 'freq' in data_code:
                self.freq = data_code['freq']
                self.freqs = []
                self.name = 'freq'
                del data_code['freq']
                del data_code['name']
            elif 'mode' in data_code:
                self.mode = data_code['mode']
                self.modes = []
                self.name = 'mode'
                del data_code['mode']
                del data_code['name']
            elif 'lsdvmn' in data_code:
                self.lsdvmn = data_code['lsdvmn']
                self.lsdvmns = []
                self.name = 'lsdvmn'
                del data_code['lsdvmn']
            elif 'lftsfq' in data_code:
                self.lftsfq = data_code['lftsfq']
                self.lftsfqs = []
                self.name = 'lftsfq'
                del data_code['lftsfq']
            elif 'load_set' in data_code:
                self.load_set = data_code['load_set']
                self.load_sets = []
                self.name = 'load_set'
                del data_code['load_set']
            elif 'load_step' in data_code:
                self.load_step = data_code['load_step']
                self.load_steps = []
                self.name = 'load_step'
                del data_code['load_step']
            elif 'loadID' in data_code:
                self.loadID = data_code['loadID']
                self.loadIDs = []
                self.name = 'loadID'
                del data_code['loadID']
            else:
                raise RuntimeError('time/mode/freq/lsdvmn/dt/load_set/lftsfq/loadID/load_set/load_step doesnt exist in data_code=%s' % str(data_code))
            #del data_code['log']


            self.log = data_code['log']
            #del data_code['log']
            if not len(data_code) == 0:
                #print(data_code)
                pass

        self.apply_data_code()
        self.set_data_members()
        #self.log.debug(self.code_information())

    def isImaginary(self):
        return bool(self.sort_bits[1])

    def _write_matlab_args(self, name, isubcase, f):
        for key, value, in sorted(self.data_code.iteritems()):
            if key is not 'log':
                if isinstance(value, str):
                    value = "'%s'" % value
                    msg = 'fem.%s(%i).%s = %s;\n' % (name, isubcase, key, value)
                elif isinstance(value, list) and isinstance(value[0], str):
                    msgTemp = "','".join(value)
                    msg = "fem.%s(%i).%s = {'%s'};\n" % (
                        name, isubcase, key, msgTemp)

                elif value is None:
                    value = "'%s'" % value
                else:
                    msg = 'fem.%s(%i).%s = %s;\n' % (
                        name, isubcase, key, value)
                f.write(msg)

    def apply_data_code(self):
        self.log = self.data_code['log']
        for key, value in sorted(self.data_code.iteritems()):
            if key is not 'log':
                self.__setattr__(key, value)
                #self.log.debug("  key=%s value=%s" %(key, value))
                #print "  key=%s value=%s" %(key, value)
        #self.log.debug("")

    def get_data_code(self):
        msg = []
        if 'dataNames' not in self.data_code:
            return []

        for name in self.data_code['dataNames']:
            try:
                if hasattr(self, name + 's'):
                    vals = getattr(self, name + 's')
                    name = name + 's'
                else:
                    vals = getattr(self, name)
                msg.append('  %s = %s\n' % (name, list_print(vals)))
            except AttributeError:  # weird case...
                pass
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

    def append_data_member(self, varName, valueName):
        """
        this appends a data member to a variable that may or may not exist
        """
        hasList = self.start_data_member(varName, valueName)
        if hasList:
            listA = self.getVar(varName)
            if listA is not None:
                #print "has %s" %(varName)
                value = self.getVar(valueName)
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
            #print "name = ",name
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
        then to print their values.

        This way there is no dependency on one result type having ['mode'] and
        another result type having ['mode','eigr','eigi'].
        """
        keyVals = []
        for name in self.data_code['dataNames']:
            vals = getattr(self, name + 's')
            keyVals.append(vals)
            #print "%ss = %s" %(name,vals)

        msg = ''
        for name in self.data_code['dataNames']:
            msg += '%-10s ' % name
        msg += '\n'

        nModes = len(keyVals[0])
        for i in xrange(nModes):
            for vals in keyVals:
                msg += '%-10g ' % vals[i]
            msg += '\n'
        return msg + '\n'

    def recastGridType(self, grid_type):
        """converts a grid_type integer to a string"""
        if grid_type == 1:
            Type = 'G'  # GRID
        elif grid_type == 2:
            Type = 'S'  # SPOINT
        elif grid_type == 7:
            Type = 'L'  # RIGID POINT (e.g. RBE3)
        elif grid_type == 0:
            Type = 'H'      # SECTOR/HARMONIC/RING POINT
        else:
            raise RuntimeError('grid_type=%s' % grid_type)
        return Type

    def cast_grid_type(self, gridType):
        """converts a gridType integer to a string"""
        if gridType == 'G':
            Type = 1  # GRID
        elif gridType == 'S':
            Type = 2  # SPOINT
        elif gridType == 'L':
            Type = 7  # RIGID POINT (e.g. RBE3)
        elif gridType == 'H':
            Type = 0      # SECTOR/HARMONIC/RING POINT
        else:
            raise RuntimeError('gridType=%s' % gridType)
        return Type

    def update_dt(self, data_code, dt):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        self.data_code = data_code
        self.apply_data_code()
        msg = 'update_dt not implemented in the %s class' % self.__class__.__name__
        raise RuntimeError(msg)
        #assert dt>=0.
        #print "updating dt...dt=%s" %(dt)
        if dt is not None:
            self.dt = dt
            self.add_new_transient()