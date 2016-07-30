# pylint: disable=C0103
"""
Subcase creation/extraction class
"""
from __future__ import print_function
from six import string_types, iteritems, PY2, PY3
from numpy import ndarray

from pyNastran.utils import integer_types
from pyNastran.bdf.cards.base_card import deprecated
from pyNastran.bdf.cards.collpase_card import collapse_thru_packs
from pyNastran.bdf.bdf_interface.assign_type import interpret_value

int_cards = (
    # these are cards that look like:
    #    LOAD = 6
    'SPC', 'MPC', 'TRIM', 'FMETHOD', 'METHOD', 'LOAD',
    'SUPORT', 'SUPORT1', 'TEMPERATURE(INITIAL)', 'TEMPERATURE(LOAD)',
    'DLOAD', 'MFLUID', 'CLOAD', 'NLPARM', 'CMETHOD',
    'FREQUENCY', 'TSTEP', 'TSTEPNL', 'SDAMPING', 'DESOBJ',
    'TEMPERATURE(INIT)', 'RANDOM', 'DESSUB', 'ADAPT', 'MAXLINES',
    'TFL', 'DESGLB', 'SMETHOD', 'DYNRED', 'GUST', 'TEMPERATURE(MATE)',
    'OTIME', 'NONLINEAR', 'AUXM', 'IC', 'BC', 'OUTRCV', 'DIVERG',
    'DATAREC', 'TEMPERATURE(BOTH)', 'DEFORM', 'MODES', 'CASE',
    'SEDR', 'SELG', 'SEFINAL', 'SEKR', 'TEMPERATURE(ESTIMATE)',
    'GPSDCON', 'AUXMODEL',
    'MODTRAK', 'OFREQ', 'DRSPAN', 'OMODES', 'ADACT', 'SERESP', 'STATSUB',
    'CURVESYM', 'ELSDCON', 'CSSCHD', 'NSM', 'TSTRU', 'RANDVAR', ''
    'RGYRO', 'SELR', 'TEMPERATURE(ESTI)', 'RCROSS', 'SERE', 'SEMR',
)

plottable_types = (
    # these are types that look like:
    #    STRESS(PLOT,PRINT,PUNCH,SORT1) = ALL
    # they all support PLOT
    'STRESS', 'STRAIN', 'SPCFORCES', 'DISPLACEMENT', 'MPCFORCES', 'SVECTOR',
    'VELOCITY', 'ACCELERATION', 'FORCE', 'ESE', 'OLOAD', 'SEALL', 'GPFORCE',
    'GPSTRESS', 'GPSTRAIN', 'FLUX', 'AEROF', 'THERMAL', 'STRFIELD',
    'NOUTPUT', 'SEDV', 'APRES', 'HTFLOW', 'NLSTRESS', 'GPKE',
    'SACCELERATION', 'SDISPLACEMENT', 'SEMG', 'HARMONICS', 'PRESSURE', 'VUGRID',
    'ELSUM', 'SVELOCITY', 'STRFIELD REAL', 'SENSITY', 'MONITOR',
    'NLLOAD', 'GPSDCON', 'BOUTPUT',
)


class Subcase(object):
    """
    Subcase creation/extraction class
    """
    solCodeMap = {
        1: 101,
        21: 101,
        24: 101,
        26: 101,
        61: 101,
        64: 106,  # correct
        66: 106,  # correct
        68: 106,  # correct
        76: 101,
        99: 129,  # correct
        144: 101,  # correct
        187: 101,
    }

    def __init__(self, id=0):
        self.id = id
        self.params = {}
        self.sol = None
        #print("\n***adding subcase %s***" % self.id)

    def deprecated(self, old_name, new_name, deprecated_version):
        return deprecated(old_name, new_name, deprecated_version, levels=[0, 1, 2])

    def get_stress_code(self, key, options, value):
        """
        Method get_stress_code:

        .. note:: the individual element must take the stress_code and reduce
        it to what the element can return.  For example, for an isotropic
        CQUAD4 the fiber field doesnt mean anything.

        BAR - no von mises/fiber
        ISOTROPIC - no fiber

        .. todo:: how does the MATERIAL bit get turned on?  I'm assuming it's
                  element dependent...
        """
        stress_code = 0
        if 'VONMISES' in options:
            stress_code += 1
        if key == 'STRAIN':
            stress_code += 10  # 2+8=10 - fields 2 and 4
        if 'FIBER' in options:
            stress_code += 4
        #if 'MATERIAL' in options:
        #    stress_code += 16  material coord (1) vs element (0)
        return stress_code

    def add_op2_data(self, data_code, msg, log=None):
        """
        >>> self.subcase.add_op2_data(self.data_code, 'VECTOR')
        """
        #subtable_name = data_code['subtable_name']
        table_name = data_code['table_name']
        if PY2 and isinstance(table_name, str):
            # tablename is a byte string
            table_name = table_name = table_name.decode('latin1')
        elif PY3 and not isinstance(table_name, str):
            # table_name is a byte string
            table_name = table_name = table_name.decode('latin1')
        else:
            raise NotImplementedError('table_name=%r PY2=%s PY3=%s' % (table_name, PY2, PY3))

        table_code = data_code['table_code']
        sort_code = data_code['sort_code']
        device_code = data_code['device_code']
        #print(data_code)
        #print('table_name=%r table_code=%s sort_code=%r device_code=%r' % (
            #table_name, table_code, sort_code, device_code))
        table_name = table_name.strip()
        #if 'TITLE' in
        #print(data_code)
        options = []
        if data_code['title']:
            self.add('TITLE', data_code['title'], options, 'STRING-type')
        if data_code['subtitle']:
            self.add('SUBTITLE', data_code['subtitle'], options, 'STRING-type')
        if data_code['label']:
            self.add('LABEL', data_code['label'], options, 'STRING-type')

        if table_name in ['OUGV1', 'BOUGV1', 'OUGV2', 'OUG1']:
            if table_code == 1:
                thermal = data_code['thermal']
                if thermal == 0:
                    self.add('DISPLACEMENT', 'ALL', options, 'STRESS-type')
                elif thermal == 1:
                    self.add('ANALYSIS', 'HEAT', options, 'KEY-type')
                else:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                    #raise NotImplementedError(data_code)
            elif table_code == 7:
                self.add('VECTOR', 'ALL', options, 'STRESS-type')
            elif table_code == 10:
                self.add('VELOCITY', 'ALL', options, 'STRESS-type')
            elif table_code == 11:
                self.add('ACCELERATION', 'ALL', options, 'STRESS-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name == 'TOUGV1':
            thermal = data_code['thermal']
            if thermal == 1:
                self.add('ANALYSIS', 'HEAT', options, 'KEY-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name == 'ROUGV1':
            thermal = data_code['thermal']
            if thermal == 0:
                self.add('DISPLACEMENT', 'ALL', options, 'STRESS-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name == 'BOPHIG':
            if table_code == 7:
                self.add('ANALYSIS', 'HEAT', options, 'KEY-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name == 'OUPV1':
            if table_code == 1:
                self.add('SDISPLACEMENT', 'ALL', options, 'STRESS-type')
            elif table_code == 10:
                self.add('SVELOCITY', 'ALL', options, 'STRESS-type')
            elif table_code == 11:
                self.add('SACCELERATION', 'ALL', options, 'STRESS-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)

        elif table_name in ['OQG1', 'OQG2']:
            if table_code == 3:
                self.add('SPCFORCES', 'ALL', options, 'STRESS-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name in ['OEF1X', 'OEF1']:
            if table_code in [4]:
                self.add('FORCE', 'ALL', options, 'STRESS-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name in ['OEFIT']:
            if table_code in [25]:
                self.add('FORCE', 'ALL', options, 'STRESS-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name == 'OQMG1':
            if table_code in [3, 39]:
                self.add('MPCFORCES', 'ALL', options, 'STRESS-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name in ['OGPFB1']:
            if table_code == 19:
                self.add('GPFORCE', 'ALL', options, 'STRESS-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name in ['OES1', 'OES1X', 'OES1X1', 'OES1C', 'OESCP',
                            'OESNLXD', 'OESNLXR', 'OESNLBR', 'OESTRCP']:
            #assert data_code['is_stress_flag'] == True, data_code
            if table_code == 5:
                self.add('STRESS', 'ALL', options, 'STRESS-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name in ['OESRT']:
            #assert data_code['is_stress_flag'] == True, data_code
            if table_code == 25:
                self.add('STRESS', 'ALL', options, 'STRESS-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name in ['OSTR1X', 'OSTR1C']:
            assert data_code['is_strain_flag'] == True, data_code
            if table_code == 5:
                self.add('STRAIN', 'ALL', options, 'STRESS-type')
            else:
                if log is None:
                    print('Error calling subcase.add_op2_data...')
                    print(msg)
                    print(data_code)
                else:
                    self.log.error(msg)
                    self.log.error(data_code)
                #raise NotImplementedError(data_code)
        elif table_name in ['RADCONS', 'RADEFFM', 'RADEATC', 'RAPEATC', 'RAQEATC', 'RADCONS',
                            'RASEATC', 'RAFEATC', 'RAEEATC', 'RANEATC', 'RAGEATC',]:
            pass
        else:
            if log is None:
                print('Error calling subcase.add_op2_data...')
                print(msg)
                print(data_code)
            else:
                self.log.error(msg)
                self.log.error(data_code)
            #raise NotImplementedError(data_code)
        #print(self)

    def get_format_code(self, options, value):
        """
        Gets the format code that will be used by the op2 based on
        the options.

        Parameters
        ----------
        options : list[int/float/str]
            the options for a parameter
        value : int/float/str
            the value of the parameter

        .. todo::  not done...only supports REAL, IMAG, PHASE, not RANDOM
        """
        format_code = 0
        if 'REAL' in options:
            format_code += 1
        if 'IMAG' in options:
            format_code += 2
        if 'PHASE' in options:
            format_code += 4
        format_code = max(format_code, 1)
        return format_code

    def get_sort_code(self, options, value):
        """
        Gets the sort code of a given set of options and value

        Parameters
        ----------
        options : List[int/str]
            the options for a parameter
        value : int; str
            the value of the parameter
        """
        sort_code = 0
        if 'COMPLEX' in options:
            sort_code += 1
        if 'SORT2' in options:
            sort_code += 2
        if 'RANDOM' in options:
            sort_code += 4
        return sort_code

    def get_device_code(self, options, value):
        """
        Gets the device code of a given set of options and value

        Parameters
        ----------
        options : list[int/float/str]
            the options for a parameter
        value : int/float/str
            the value of the parameter

        Returns
        -------
        device_code : int
           The OP2 device code

           0 - No output
           1 - PRINT
           2 - PLOT
           3 - PRINT, PLOT
           4 - PUNCH
           5 - PRINT, PUNCH
           6 - PRINT, PLOT, PUNCH
        """
        device_code = 0
        if 'PRINT' in options:
            device_code += 1
        if 'PLOT' in options:
            device_code += 2
        if 'PUNCH' in options:
            device_code += 4
        device_code = max(device_code, 1)
        #if device_code==0:
        #    device_code=1  # PRINT
        return device_code

    def get_analysis_code(self, sol):
        """
        Maps the solution number to the OP2 analysis code.

         * 8 - post-buckling (maybe 7 depending on NLPARM???)

        # not important
         * 3/4 - differential stiffness (obsolete)
         * 11  - old geometric nonlinear statics
         * 12  - contran (???)

        .. todo:: verify
        """
        codes = {
            101 : 1,  # staics
            103 : 2,  # modes
            105 : 7,  # pre-buckling
            106 : 10, # nonlinear statics
            107 : 9,  # complex eigenvalues
            108 : 5,  # frequency
            111 : 5,
            112 : 6,
            114 : 1,
            115 : 2,
            116 : 7,
            118 : 5,
            129 : 6,  # nonlinear
            144 : 1,  # static aero
            145 : 1,
            146 : 1,  # flutter
            153 : 10,
            159 : 6,  # transient thermal
        }
        #print("sol=%s" % sol)
        approach_code = codes[sol]
        #print('approach_code = %s' % approach_code)
        return approach_code

    def get_table_code(self, sol, table_name, options):
        """
        Gets the table code of a given parameter.  For example, the
        DISPLACMENT(PLOT,POST)=ALL makes an OUGV1 table and stores the
        displacement.  This has an OP2 table code of 1, unless you're running a
        modal solution, in which case it makes an OUGV1 table of eigenvectors
        and has a table code of 7.

        Parameters
        ----------
        options : list[int/float/str]
            the options for a parameter
        value : int/float/str
            the value of the parameter

        Returns
        -------
        table_code : int
           the OP2 table_code
        """
        if table_name in ['VECTOR', 'PRESSURE']:
            table_name = 'DISPLACEMENT'  # equivalent tables...

        key = (sol, table_name)
        tables = {
            # SOL, table_name      table_code
            (101, 'ACCELERATION'): 11,
            (103, 'ACCELERATION'): 11,
            (106, 'ACCELERATION'): 11,
            (107, 'ACCELERATION'): 11,
            (108, 'ACCELERATION'): 11,
            (129, 'ACCELERATION'): 11,
            #(144, 'ACCELERATION'): 11,
            (145, 'ACCELERATION'): 11,
            (146, 'ACCELERATION'): 11,

            (101, 'DISPLACEMENT'): 1,
            (103, 'DISPLACEMENT'): 7,  # VECTOR
            (105, 'DISPLACEMENT'): 7,
            (106, 'DISPLACEMENT'): 1,
            (107, 'DISPLACEMENT'): 7,
            (108, 'DISPLACEMENT'): 1,
            (109, 'DISPLACEMENT'): 1,
            (111, 'DISPLACEMENT'): 7,
            (112, 'DISPLACEMENT'): 1,
            (129, 'DISPLACEMENT'): 7,
            #(144, 'DISPLACEMENT'): 1,
            (145, 'DISPLACEMENT'): 1,
            (146, 'DISPLACEMENT'): 1,

            (101, 'ESE'): 18,  # energy
            (103, 'ESE'): 18,  # energy
            (105, 'ESE'): 18,  # energy
            (106, 'ESE'): 18,  # energy
            (107, 'ESE'): 18,  # energy
            (108, 'ESE'): 18,  # energy
            (109, 'ESE'): 18,  # energy
            (110, 'ESE'): 18,  # energy
            (111, 'ESE'): 18,  # energy
            (112, 'ESE'): 18,  # energy
            (145, 'ESE'): 18,  # energy
            (146, 'ESE'): 18,  # energy

            (101, 'FORCE'): 3,  # ???
            (103, 'FORCE'): 3,  # ???
            (105, 'FORCE'): 3,  # ???
            (106, 'FORCE'): 3,  # ???
            (107, 'FORCE'): 4,  # ???
            (108, 'FORCE'): 3,  # ???
            (111, 'FORCE'): 3,  # ???
            (112, 'FORCE'): 3,  # ???
            (129, 'FORCE'): 3,  # ???
            (145, 'FORCE'): 3,  # ???
            (146, 'FORCE'): 3,  # ???

            (101, 'GPFORCE'): 19,
            (105, 'GPFORCE'): 19,
            (106, 'GPFORCE'): 19,
            (107, 'GPFORCE'): 19,
            (108, 'GPFORCE'): 19,
            (111, 'GPFORCE'): 19,
            (112, 'GPFORCE'): 19,
            (129, 'GPFORCE'): 19,
            (145, 'GPFORCE'): 19,
            (146, 'GPFORCE'): 19,

            (101, 'GPSTRESS'): 20,
            (105, 'GPSTRESS'): 20,
            (106, 'GPSTRESS'): 20,
            (107, 'GPSTRESS'): 20,
            (108, 'GPSTRESS'): 20,
            (111, 'GPSTRESS'): 20,
            (112, 'GPSTRESS'): 20,
            (129, 'GPSTRESS'): 20,
            (145, 'GPSTRESS'): 20,
            (146, 'GPSTRESS'): 20,

            (101, 'GPSTRAIN'): 21,
            (105, 'GPSTRAIN'): 21,
            (106, 'GPSTRAIN'): 21,
            (107, 'GPSTRAIN'): 21,
            (108, 'GPSTRAIN'): 21,
            (111, 'GPSTRAIN'): 21,
            (112, 'GPSTRAIN'): 21,
            (129, 'GPSTRAIN'): 21,
            (145, 'GPSTRAIN'): 21,
            (146, 'GPSTRAIN'): 21,

            (101, 'MPCFORCES'): 3,
            (103, 'MPCFORCES'): 3,
            (106, 'MPCFORCES'): 3,
            (108, 'MPCFORCES'): 3,
            (112, 'MPCFORCES'): 3,
            (129, 'MPCFORCES'): 3,
            #(144, 'MPCFORCES'): 3,
            (145, 'MPCFORCES'): 3,
            (146, 'MPCFORCES'): 3,

            (101, 'OLOAD'): 2,
            (103, 'OLOAD'): 2,
            (105, 'OLOAD'): 2,
            (106, 'OLOAD'): 2,
            (107, 'OLOAD'): 2,
            (108, 'OLOAD'): 2,
            (111, 'OLOAD'): 2,
            (112, 'OLOAD'): 2,
            (129, 'OLOAD'): 2,
            #(144, 'OLOAD'): 2,
            (145, 'OLOAD'): 2,
            (146, 'OLOAD'): 2,

            (101, 'SPCFORCES'): 3,
            (103, 'SPCFORCES'): 3,
            (105, 'SPCFORCES'): 3,
            (106, 'SPCFORCES'): 3,
            (107, 'SPCFORCES'): 3,
            (108, 'SPCFORCES'): 3,
            (110, 'SPCFORCES'): 3,
            (111, 'SPCFORCES'): 3,
            (112, 'SPCFORCES'): 3,
            (129, 'SPCFORCES'): 3,
            #(144, 'SPCFORCES'): 3,
            (145, 'SPCFORCES'): 3,
            (146, 'SPCFORCES'): 3,

            (101, 'STRAIN'): 5,  # 5/20/21 ???
            (105, 'STRAIN'): 5,
            (106, 'STRAIN'): 5,
            (107, 'STRAIN'): 5,
            (108, 'STRAIN'): 5,
            (110, 'STRAIN'): 5,
            (111, 'STRAIN'): 5,
            (112, 'STRAIN'): 5,
            (129, 'STRAIN'): 5,
            (145, 'STRAIN'): 5,
            (146, 'STRAIN'): 5,

            (101, 'STRESS'): 5,  # 5/20/21 ???
            (103, 'STRESS'): 5,
            (105, 'STRESS'): 5,
            (106, 'STRESS'): 5,
            (107, 'STRESS'): 5,
            (108, 'STRESS'): 5,
            (111, 'STRESS'): 5,
            (112, 'STRESS'): 5,
            (129, 'STRESS'): 5,
            (145, 'STRESS'): 5,
            (146, 'STRESS'): 5,

            (145, 'SVECTOR'): 14,

            (101, 'FLUX'): 4,
            (103, 'FLUX'): 4,
            (106, 'FLUX'): 4,
            (112, 'FLUX'): 4,
            (108, 'FLUX'): 4,
            (153, 'FLUX'): 4,
            (159, 'FLUX'): 4,

            (101, 'THERMAL'): 3,  # 3/4 ???
            (159, 'THERMAL'): 3,  # 3/4 ???

            (101, 'VELOCITY'): 10,
            (103, 'VELOCITY'): 10,
            (106, 'VELOCITY'): 10,
            (107, 'VELOCITY'): 10,
            (108, 'VELOCITY'): 10,
            (111, 'VELOCITY'): 10,
            (112, 'VELOCITY'): 10,
            (129, 'VELOCITY'): 10,
            #(144, 'VELOCITY'): 10,
            (145, 'VELOCITY'): 10,
            (146, 'VELOCITY'): 10,

            (101, 'VUGRID'): 10,
        }
        #print("key=%s" % str(key))
        if key not in tables:
            raise KeyError(key)
        table_code = tables[key]
        return table_code

    def __contains__(self, param_name):
        """
        Checks to see if a parameter name is in the subcase.

        Parameters
        ----------
        param_name : str
            the case control parameters to check for

        .. code-block:: python

          model = BDF()
          model.read_bdf(bdf_filename)
          case_control = model.case_control_deck
          subcase1 = case_control.subcases[1]
          if 'LOAD' in subcase1:
              print('found LOAD for subcase 1')
        """
        if param_name in self.params:
            return True
        return False

    def has_parameter(self, *param_names):
        """
        Checks to see if one or more parameter names are in the subcase.

        Parameters
        ----------
        param_names : str; List[str]
            the case control parameters to check for

        Returns
        -------
        exists : List[bool]
            do the parameters exist

        .. code-block:: python

          model = BDF()
          model.read_bdf(bdf_filename)
          case_control = model.case_control_deck
          subcase1 = case_control.subcases[1]
          if any(subcase1.has_parameter('LOAD', 'TEMPERATURE(LOAD)')):
              print('found LOAD for subcase 1')
        """
        exists = [True if param_name.upper() in self.params else False for param_name in param_names]
        return exists

    def __getitem__(self, param_name):
        """
        Gets the [value, options] for a subcase.

        Parameters
        ----------
        param_name : str
            the case control parameters to get

        Returns
        -------
        value : varies
            the value of the parameter
            'ALL' in STRESS(PLOT,PRINT) = ALL
        options : List[varies]
            the values in parentheses
            ['PLOT', 'PRINT'] in STRESS(PLOT,PRINT) = ALL

        .. code-block:: python

          model = BDF()
          model.read_bdf(bdf_filename)
          case_control = model.case_control_deck
          subcase1 = case_control.subcases[1]
          value, options = subcase1['LOAD']
        """
        return self.get_parameter(param_name)

    def suppress_output(self, suppress_to='PLOT'):
        """
        Replaces F06 printing with OP2 printing

        Converts:
            STRESS(PRINT,SORT1,REAL)
            FORCE(PRINT,PLOT,SORT1,REAL)

        to:
            STRESS(PLOT,SORT1,REAL)
            FORCE(PLOT,SORT1,REAL)

        .. warning:: needs more validation
        """
        for key, param in iteritems(self.params):
            value, options, param_type = param
            if key in int_cards or key in ('SUBTITLE', 'LABEL', 'TITLE', 'ECHO'):
                pass
            elif key in plottable_types:
                if suppress_to not in options:
                    param[1].append(suppress_to)
                if 'PRINT' in options:
                    param[1].remove('PRINT')
            else:
                raise NotImplementedError(key)

    def get_parameter(self, param_name, msg=''):  # possibly deprecate...
        """
        Gets the [value, options] for a subcase.

        Parameters
        ----------
        param_name : str
            the case control parameters to get

        Returns
        -------
        value : varies
            the value of the parameter
            'ALL' in STRESS(PLOT,PRINT) = ALL
        options : List[varies]
            the values in parentheses
            ['PLOT', 'PRINT'] in STRESS(PLOT,PRINT) = ALL

        .. code-block:: python

          model = BDF()
          model.read_bdf(bdf_filename)
          case_control = model.case_control_deck
          subcase1 = case_control.subcases[1]
          value, options = subcase1['LOAD']
        """
        param_name = update_param_name(param_name)
        if param_name not in self.params:
            raise KeyError('%s doesnt exist in subcase=%s in the case '
                           'control deck%s.' % (param_name, self.id, msg))
        value, options, param_type = self.params[param_name]
        #print('param_name=%r value=%s options=%s param_type=%r' % (param_name, value, options, param_type))
        return value, options  #, param_type

    def update_parameter_in_subcase(self, key, value, options, param_type):
        self.deprecated('update_parameter_in_subcase', 'update', '0.8')
        self.update(key, value, options, param_type)

    def add_parameter_to_subcase(self, key, value, options, param_type):
        self.deprecated('add_parameter_to_subcase', 'add', '0.8')
        self.add(key, value, options, param_type)

    def add(self, key, value, options, param_type):
        assert param_type in ['SET-type', 'CSV-type', 'SUBCASE-type', 'KEY-type', 'STRESS-type', 'STRING-type'], param_type
        self._add_data(key, value, options, param_type)

    def update(self, key, value, options, param_type):
        assert param_type in ['SET-type', 'CSV-type', 'SUBCASE-type', 'KEY-type', 'STRESS-type', 'STRING-type'], param_type
        assert key in self.params, 'key=%r is not in isubcase=%s' % (key, self.id)
        self._add_data(key, value, options, param_type)

    def _add_data(self, key, value, options, param_type):
        key = update_param_name(key)
        if key == 'ANALYSIS' and value == 'FLUT':
                value = 'FLUTTER'

        #print("adding isubcase=%s key=%r value=%r options=%r "
        #      "param_type=%r" %(self.id, key, value, options, param_type))
        if isinstance(value, string_types) and value.isdigit():
            value = int(value)

        (key, value, options) = self._simplify_data(key, value, options, param_type)
        self.params[key] = [value, options, param_type]

    def _simplify_data(self, key, value, options, param_type):
        if param_type == 'SET-type':
            #print("adding isubcase=%s key=%r value=%r options=%r "
                  #"param_type=%r" % (self.id, key, value, options, param_type))
            values2 = expand_thru_case_control(value)
            assert isinstance(values2, list), type(values2)
            if isinstance(options, list):
                msg = 'invalid type for options=%s value; expected an integer; got a list' % key
                raise TypeError(msg)
            options = int(options)
            return (key, values2, options)

        elif param_type == 'CSV-type':
            #print("adding isubcase=%s key=%r value=|%s| options=|%s| "
            #      "param_type=%s" %(self.id, key, value, options, param_type))
            if value.isdigit():  # PARAM,DBFIXED,-1
                value = value
        else:
            if 0:
                a = 'key=%r' % key
                b = 'value=%r' % value
                c = 'options=%r' % options
                d = 'param_type=%r' % param_type
                #print("_adding isubcase=%s %-18s %-12s %-12s %-12s" %(self.id, a, b, c, d))
            if isinstance(value, integer_types) or value is None:
                pass
            elif isinstance(value, (list, ndarray)):  # new???
                msg = 'invalid type for key=%s value; expected an integer; got a list' % key
                raise TypeError(msg)
            elif value.isdigit():  # STRESS = ALL
                value = value
            #else: pass
        return key, value, options

    def get_op2_data(self, sol, solmap_toValue):
        self.sol = sol
        label = 'SUBCASE %s' % (self.id)
        op2Params = {'isubcase': None, 'tables': [], 'analysisCodes': [],
                     'device_codes': [], 'sortCodes': [], 'tableCodes': [],
                     'label': label, 'subtitle': None, 'title': None,
                     'formatCodes': [], 'stressCodes': [], 'thermal': None}

        results = ['DISPLACEMENT', 'EKE', 'EDE', 'ELSDCON', 'ENTHALPY',
                   'EQUILIBRIUM', 'ESE', 'FLUX', 'FORCE', 'GPFORCE', 'GPKE',
                   'GPSDCON', 'GPSTRAIN', 'GPSTRESS', 'HOUTPUT', 'MODALKE',
                   'MODALSE', 'MPCFORCES', 'NLSTRESS', 'NOUTPUT', 'OLOAD',
                   'PFGRID', 'PFMODE', 'PFPANEL', 'RCROSS', 'RESVEC',
                   'SACCELERATION', 'SDISPACEMENT', 'SPCFORCES', 'STRAIN',
                   'STRESS', 'SVECTOR', 'SVELOCITY', 'THERMAL', 'VECTOR',
                   'VELOCITY', 'VUGRID', 'WEIGHTCHECK']

        # converts from solution 200 to solution 144
        if self.sol == 200 or 'ANALYSIS' in self.params:
            param = self.params['ANALYSIS']
            (value, options, paramType) = param

            sol = solmap_toValue[value.upper()]
            #print("***value=%s sol=%s" % (value, sol))
        else:  # leaves SOL the same
            sol = self.sol

        if sol in self.solCodeMap:  # reduces SOL 144 to SOL 101
            sol = self.solCodeMap[sol]

        for (key, param) in iteritems(self.params):
            key = key.upper()
            (value, options, paramType) = param
            #msg = ("  -key=|%s| value=|%s| options=%s paramType=|%s|"
            #    % (key, value, options, paramType))

        thermal = 0
        for (key, param) in iteritems(self.params):
            key = key.upper()
            (value, options, paramType) = param
            #msg = ("  *key=|%s| value=|%s| options=%s paramType=|%s|"
            #    % (key, value, options, paramType)
            #print(msg)
            #msg += self.printParam(key, param)
            if paramType == 'SUBCASE-type':
                op2Params['isubcase'].append(value)
            elif key in ['BEGIN', 'ECHO', 'ANALYSIS'] or 'SET' in key:
                pass
            elif key == 'TEMPERATURE':
                thermal = 1
            elif key in results:
                sort_code = self.get_sort_code(options, value)
                device_code = self.get_device_code(options, value)

                if key in ['STRESS', 'STRAIN']:
                    stressCode = self.get_stress_code(key, options, value)
                    op2Params['stressCodes'].append(stressCode)
                else:
                    op2Params['stressCodes'].append(0)

                format_code = self.get_format_code(options, value)
                table_code = self.get_table_code(sol, key, options)
                analysis_code = self.get_analysis_code(sol)

                approach_code = analysis_code * 10 + device_code
                tCode = table_code * 1000 + sort_code
                op2Params['tables'].append(key)

                op2Params['analysisCodes'].append(analysis_code)
                op2Params['approachCodes'].append(approach_code)
                op2Params['device_codes'].append(device_code)
                op2Params['formatCodes'].append(format_code)
                op2Params['sortCodes'].append(sort_code)
                op2Params['tableCodes'].append(table_code)
                op2Params['tCodes'].append(tCode)
                #analysisMethod = value

            #elif key in ['ADACT', 'ADAPT', 'AERCONFIG', 'TITLE', 'SUBTITLE',
            #             'LABEL', 'LOAD', 'SUPORT', 'SUPORT1', 'MPC', 'SPC',
            #            'TSTEPNL', 'NLPARM', 'TRIM', 'GUST', 'METHOD',
            #            'DESOBJ', 'DESSUB', 'FMETHOD', 'SEALL']:
            else:
                op2Params[key.lower()] = value

            #else:
            #    raise NotImplementedErrror('unsupported entry...\n%s' %(msg))

        op2Params['thermal'] = thermal

        #print("\nThe estimated results...")
        #for (key, value) in sorted(iteritems(op2Params)):
            #if value is not None:
                #print("   key=%r value=%r" % (key, value))

    def print_param(self, key, param):
        """
        Prints a single entry of the a subcase from the global or local
        subcase list.
        """
        msg = ''
        #msg += 'id=%s   ' %(self.id)
        (value, options, param_type) = param

        spaces = ''
        if self.id > 0:
            spaces = '    '

        #print('key=%s param=%s param_type=%s' % (key, param, param_type))
        if param_type == 'SUBCASE-type':
            if self.id > 0:
                msgi = 'SUBCASE %s\n' % (self.id)
                assert len(msgi) < 72, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
                msg += msgi
            #else:  global subcase ID=0 and is not printed
            #    pass
        elif param_type == 'KEY-type':
            #print (KEY-TYPE:  %r" % value)
            assert value is not None, param
            if ',' in value:
                sline = value.split(',')
                two_spaces = ',\n' + 2 * spaces
                msgi = spaces + two_spaces.join(sline) + '\n'
                assert len(msgi) < 68, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
                msg += msgi
            else:
                msgi = spaces + '%s\n' % value
                #assert len(msgi) < 68, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
                msg += msgi
        elif param_type == 'STRING-type':
            msgi = spaces + '%s = %s\n' % (key, value)
            if key not in ['TITLE', 'LABEL', 'SUBTITLE']:
                assert len(msgi) < 68, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
            msg += msgi
        elif param_type == 'CSV-type':
            msgi = spaces + '%s,%s,%s\n' % (key, value, options)
            assert len(msgi) < 68, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
            msg += msgi
        elif param_type == 'STRESS-type':
            str_options = ','.join(options)
            #print("str_options = %r" % (str_options))
            #print("STRESSTYPE key=%s value=%s options=%s"
            #    % (key, value, options))
            if value is None:
                val = ''
            else:
                val = ' = %s' % value
            if len(str_options) > 0:
                msgi = '%s(%s)%s\n' % (key, str_options, val)
                if len(msgi) > 64:
                    msgi = '%s(' % key
                    msg_done = ''
                    i = 0
                    while i < len(options):
                        option = options[i]
                        new_msgi = '%s,' % options[i]
                        if (len(msgi) + len(new_msgi)) < 64:
                            msgi += new_msgi
                        else:
                            msg_done += msgi + '\n'
                            msgi = spaces + new_msgi
                        i += 1
                    msg_done += msgi
                    msgi = ''
                    msg_done = msg_done.rstrip(' ,\n') + ')%s\n' % val
                    assert len(msgi) < 68, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
                    msg += msg_done
                else:
                    assert len(msgi) < 68, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
                    msg += msgi
            else:
                msgi = '%s%s\n' % (key, val)
                assert len(msgi) < 68, 'len(msg)=%s; msg=\n%s' % (len(msgi), msgi)
                msg += msgi
            msg = spaces + msg

        elif param_type == 'SET-type':
            #: .. todo:: collapse data...not written yet
            msg += write_set(value, options, spaces)
        else:
            # SET-type is not supported yet...
            msg = ('\nkey=%r param=%r\n'
                   'allowed_params=[SET-type, STRESS-type, STRING-type, SUBCASE-type, KEY-type]\n'
                   'CSV-type     -> PARAM,FIXEDB,-1\n'
                   'KEY-type     -> ???\n' # the catch all
                   'SET-type     -> SET 99 = 1234\n'
                   'SUBCASE-type -> ???\n'
                   'STRESS-type  -> DISP(PLOT, SORT1)=ALL\n'
                   '                STRESS(PLOT, SORT1)=ALL\n'
                   'STRING-type  -> LABEL = A label\n'
                   '                TITLE = A title\n'
                   '                SUBTITLE = A subtitle\n'
                   ''% (key, param))
            raise NotImplementedError(msg)

        #print("msg = %r" % (msg))
        return msg

    #def cross_reference(self, model):
        #"""
        #Method crossReference:

        #Parameters
        #----------
        #model : BDF()
            #the BDF object

        #.. note:: this is not integrated and probably never will be as it's
          #not really that necessary.  it's only really useful when running an
          #analysis.
        #"""
        #raise NotImplementedError()
        #print("keys = %s" % (sorted(self.params.keys())))
        #if 'LOAD' in self.params:
            #load_id = self.params['LOAD'][0]
            #load_obj = model.loads[load_id]
            #load_obj.cross_reference(model)
        #if 'SUPORT' in self.params:
            #pass
        #if 'MPC' in self.params:
            ##mpcID = self.params['MPC'][0]
            ##mpcObj = model.mpcs[mpcID]
            ##mpcObj.cross_reference(model)
            #pass
        #if 'SPC' in self.params:
            ##spcID = self.params['SPC'][0]
            ##print "SPC ID = ",spcID
            ##spcObj = model.spcObject
            ##spcObj.cross_reference(spcID, model)
            #pass
        #if 'TSTEPNL' in self.params:
            #tstepnl_id = self.params['TSTEPNL'][0]
            #tstepnl_obj = model.tstepnl[tstepnl_id]
            #tstepnl_obj.cross_reference(model)
        #if 'NLPARM' in self.params:
            #nlparm_id = self.params['NLPARM'][0]
            #nlparm_obj = model.nlparms[nlparm_id]
            #nlparm_obj.cross_reference(model)
        #if 'TRIM' in self.params:
            #trim_id = self.params['TRIM'][0]
            #trim_obj = model.trims[trim_id]
            #trim_obj.cross_reference(model)
        #if 'GUST' in self.params:
            #gust_id = self.params['GUST'][0]
            #gust_obj = model.gusts[gust_id]
            #gust_obj.cross_reference(model)
        #if 'DLOAD' in self.params:  # ???
            #pass

    def finish_subcase(self):
        """
        Removes the subcase parameter from the subcase to avoid printing it in
        a funny spot
        """
        if 'SUBCASE' in self.params:
            del self.params['SUBCASE']
        #print "self.params %s = %s" %(self.id,self.params)

    def write_subcase(self, subcase0):
        """
        Internal method to print a subcase

        Parameters
        ----------
        subcase0 : Subcase()
            the global Subcase object

        Returns
        -------
        msg : str
           the string of the current Subcase
        """
        if self.id == 0:
            msg = str(self)
        else:
            msg = 'SUBCASE %s\n' % self.id
            nparams = 0
            for (key, param) in self.subcase_sorted(self.params.items()):
                if key in subcase0.params and subcase0.params[key] == param:
                    pass  # dont write global subcase parameters
                else:
                    #print("key=%s param=%s" %(key, param))
                    (value, options, param_type) = param
                    #print("  *key=%r value=|%s| options=%s "
                          #"param_type=|%s|" % (key, value, options, param_type))
                    msg += self.print_param(key, param)
                    nparams += 1
                    #print ""
            if nparams == 0:
                for (key, param) in self.subcase_sorted(self.params.items()):
                    #print("key=%s param=%s" %(key, param))
                    (value, options, param_type) = param
                    #print("  *key=%r value=|%s| options=%s "
                          #"param_type=|%s|" % (key, value, options, param_type))
                    msg += self.print_param(key, param)
                    nparams += 1
                assert nparams > 0, 'No subcase paramters are defined for isubcase=%s...' % self.id

        return msg

    def subcase_sorted(self, lst):
        """
        Does a "smart" sort on the keys such that SET cards increment in
        numerical order.  Also puts the sets first.

        Parameters
        ----------
        lst : List[str]
            the list of subcase list objects (list_a)

        Returns
        -------
        list_b : List[str]
            the sorted version of list_a
        """
        # presort the list to put all the SET cards next to each other
        # instead of list_a.sort() as it allows lst to be any iterable
        lst = sorted(lst)

        i = 0  # needed in case the loop doesn't execute
        iset = None  # index of first SET card in the deck
        set_dict = {}
        list_before = []
        set_keys = []
        for (i, entry) in enumerate(lst):
            key = entry[0]
            if 'SET' in key[0:3]:
                if key == 'SET':  # handles "SET = ALL"
                    key = 0
                else:  # handles "SET 100 = 1,2,3"
                    sline = key.split(' ')
                    try:
                        key = int(sline[1])
                    except:
                        msg = 'error caclulating key; sline=%s' % sline
                        raise RuntimeError(msg)

                # store the integer ID and the SET-type list
                set_dict[key] = entry
                set_keys.append(key)
            else:
                # only store the entries before the SET cards
                list_before.append(entry)
                if iset:
                    break

        # grab the other entries
        list_after = lst[i + 1:]

        # write the SET cards in a sorted order
        set_list = []
        for key in sorted(set_keys):
            set_list.append(set_dict[key])

        # combine all the cards
        list_b = set_list + list_before + list_after
        return list_b

    def __repr__(self):
        """
        Prints out EVERY entry in the subcase.  Skips parameters already in
        the global subcase.

        .. note:: this function is only used for debugging.
        """
        #msg = "-------SUBCASE %s-------\n" %(self.id)
        msg = ''
        if self.id > 0:
            msg += 'SUBCASE %s\n' % self.id

        nparams = 0
        for key, param in self.subcase_sorted(iteritems(self.params)):
            (value, options, param_type) = param
            #print('key=%r value=%s options=%s' % (key, value, options))
            msg += self.print_param(key, param)
            nparams += 1
        if self.id > 0:
            assert nparams > 0, 'No subcase paramters are defined for isubcase=%s...' % self.id
        return msg

def update_param_name(param_name):
    """
    Takes an abbreviated name and expands it so the user can type DISP or
    DISPLACEMENT and get the same answer

    Parameters
    ----------
    param_name : str
        the parameter name to be standardized (e.g. DISP vs. DIPLACEMENT)

    .. todo:: not a complete list
    """
    if   param_name.startswith('ACCE'):
        param_name = 'ACCELERATION'
    elif param_name.startswith('DESO'):
        param_name = 'DESOBJ'
    elif param_name.startswith('DESS'):
        param_name = 'DESSUB'
    elif param_name.startswith('DISP'):
        param_name = 'DISPLACEMENT'
    elif param_name.startswith('EXPO'):
        param_name = 'EXPORTLID'
    elif param_name.startswith('ELFO'):
        param_name = 'FORCE'
    elif param_name.startswith('ELST'):
        param_name = 'STRESS' # or ELSTRESS
    elif param_name.startswith('FORC'):
        param_name = 'FORCE'
    elif param_name.startswith('FREQ'):
        param_name = 'FREQUENCY'
    elif param_name.startswith('GPFO'):
        param_name = 'GPFORCE'
    elif param_name == 'GPST':
        raise SyntaxError('invalid GPST stress/strain')
    elif param_name.startswith('HARMONIC'):
        param_name = 'HARMONICS'
    elif param_name.startswith('METH'):
        param_name = 'METHOD'
    elif param_name.startswith('MPCF'):
        param_name = 'MPCFORCES'
    elif param_name.startswith('NLPAR'):
        param_name = 'NLPARM'
    elif param_name.startswith('OLOA'):
        param_name = 'OLOAD'
    elif param_name.startswith('PRES'):
        param_name = 'PRESSURE'
    elif param_name.startswith('SDAMP'):
        param_name = 'SDAMPING'
    elif param_name.startswith('SDISP'):
        param_name = 'SDISPLACEMENT'
    elif param_name.startswith('SMETH'):
        param_name = 'SMETHOD'
    elif param_name.startswith('SPCF'):
        param_name = 'SPCFORCES'
    elif param_name.startswith('STRA'):
        param_name = 'STRAIN'
    elif param_name.startswith('STRE'):
        param_name = 'STRESS'
    elif param_name.startswith('SUBT'):
        param_name = 'SUBTITLE'
    elif param_name.startswith('SUPO'):
        param_name = 'SUPORT1'
    elif param_name.startswith('SVEC'):
        param_name = 'SVECTOR'
    elif param_name.startswith('SVELO'):
        param_name = 'SVELOCITY'
    elif param_name.startswith('THER'):
        param_name = 'THERMAL'
    elif param_name.startswith('VECT'):
        #param_name = 'PRESSURE' # or VECTOR
        param_name = 'DISPLACEMENT' # or VECTOR
    elif param_name.startswith('VELO'):
        param_name = 'VELOCITY'
    elif param_name.startswith('TITL'):
        param_name = 'TITLE'
    elif param_name.startswith('MAXLINE'):
        param_name = 'MAXLINES'
    elif param_name.startswith('LINE'):
        param_name = 'LINE'
    elif param_name.startswith('AXISYM'):
        param_name = 'AXISYMMETRIC'
    elif param_name.startswith('SUBSE'):
        param_name = 'SUBSEQ'
    elif param_name.startswith('XTIT'):
        param_name = 'XTITLE'
    elif param_name.startswith('YTIT'):
        param_name = 'YTITLE'
    elif param_name.startswith('SACCE'):
        param_name = 'SACCELERATION'
    elif param_name.startswith('GPSTRE'):
        param_name = 'GPSTRESS'
    elif param_name.startswith('GPSTR'):
        param_name = 'GPSTRAIN'
    elif param_name in ['DEFO', 'DEFOR']:
        param_name = 'DEFORM'
    elif param_name == 'TEMPERATURE(INIT)':
        param_name = 'TEMPERATURE(INITIAL)'

    #elif param_name.startswith('DFRE'):  param_name = 'D'

    # handled in caseControlDeck.py
    #elif param_name.startswith('TEMP'):  param_name = 'TEMPERATURE'
    #print '*param_name = ',param_name
    return param_name

def expand_thru_case_control(set_value):
    set_value2 = set([])
    add_mode = True
    imin = 0
    imax = 0
    #print('set_value = %r' % set_value)
    for ivalue in set_value:
        if isinstance(ivalue, integer_types):
            assert add_mode is True, add_mode
            set_value2.add(ivalue)
            continue
        ivalue = ivalue.strip()
        #print('  ivalue=%r; type=%s' % (ivalue, type(ivalue)))
        if '/' in ivalue:
            set_value2.add(ivalue)
        else:

            #if 'ALL' in ivalue:
                #msg = ('ALL is not supported on CaseControlDeck '
                       #'SET card\nvalue=%r\nset=%r' % (ivalue, set_value))
                #raise RuntimeError(msg)
            #elif 'EXCEPT' in ivalue:
                #msg = ('EXCEPT is not supported on CaseControlDeck '
                       #'SET card\nvalue=%r\nset=%r' % (ivalue, set_value))
                #raise RuntimeError(msg)

            ivalue = interpret_value(ivalue, card=str(set_value))
            if isinstance(ivalue, integer_types):
                #print('  isdigit')
                #ivalue = int(ivalue)
                if not imin < ivalue < imax:
                    add_mode = True
                    #print('  break...\n')
                if add_mode:
                    #print('  adding %s' % ivalue)
                    set_value2.add(ivalue)
                else:
                    #print('  removing %s' % ivalue)
                    set_value2.remove(ivalue)
                    imin = ivalue
            elif isinstance(ivalue, float):
                assert add_mode is True, add_mode
                set_value2.add(ivalue)
            elif isinstance(ivalue, string_types):
                #print('  not digit=%r' % set_value)
                if set_value is 'EXCLUDE':
                    msg = ('EXCLUDE is not supported on CaseControlDeck '
                           'SET card\n')
                    raise RuntimeError(msg)
                elif 'THRU' in ivalue:
                    svalues = ivalue.split()
                    #print('  svalues=%s' % svalues)
                    if len(svalues) == 3:
                        assert add_mode is True, add_mode
                        imin, thru, imax = svalues
                        assert thru == 'THRU', thru
                        imin = int(imin)
                        imax = int(imax)
                        assert imax > imin, 'imin=%s imax=%s' % (imin, imax)
                        for jthru in range(imin, imax + 1):
                            set_value2.add(jthru)

                    elif len(svalues) == 4:
                        imin, thru, imax, by_except = svalues
                        imin = int(imin)
                        imax = int(imax)
                        assert imax > imin, 'imin=%s imax=%s' % (imin, imax)
                        assert by_except == 'EXCEPT', by_except

                        for jthru in range(imin, imax + 1):
                            set_value2.add(jthru)
                        add_mode = False

                    elif len(svalues) == 5:
                        imin, thru, imax, by_except, increment_except = svalues
                        imin = int(imin)
                        imax = int(imax)
                        assert imax > imin, 'imin=%s imax=%s' % (imin, imax)
                        increment_except = int(increment_except)
                        if by_except == 'BY':
                            for jthru in range(imin, imax + 1, by_except):
                                set_value2.add(jthru)
                            add_mode = True
                        elif by_except == 'EXCEPT':
                            for jthru in range(imin, imax + 1):
                                if jthru == increment_except:
                                    continue
                                set_value2.add(jthru)
                            add_mode = False
                        else:
                            raise RuntimeError(ivalue)
                    else:
                        raise RuntimeError('expected data of the form: 10 THRU 20 or 10 THRU 20 BY 5; actual=\n%r' % ivalue)
                else:
                    assert add_mode is True, add_mode
                    set_value2.add(ivalue)

    list_values = list(set_value2)
    try:
        list_values.sort()
    except TypeError:
        msg = 'sort error: list_values=%s'  % (list_values)
        raise TypeError(msg)

    #print('end of expansion', list_values)
    return list_values


def write_set(value, options, spaces=''):
    """
    writes
    SET 80 = 3926, 3927, 3928, 4141, 4142, 4143, 4356, 4357, 4358, 4571,
         4572, 4573, 3323 THRU 3462, 3464 THRU 3603, 3605 THRU 3683,
         3910 THRU 3921, 4125 THRU 4136, 4340 THRU 4351

    Parameters
    ----------
    value : int
        the Set ID
    options : List[int]
        the Set values
    spaces : str; default=''
        indentation

    Returns
    -------
    msg : str
       the string of the set

    Example
    -------
    value = 80
    options = [1, 2, 3, 4, 5, 7]
    set = write_set(value, options, spaces='')
    print(set)
    >>> SET 80 = 1 THRU 5, 7
    """
    value.sort()
    starter = 'SET %s = ' % (options)
    msg2 = spaces + starter

    msg = ''
    nchars = len(msg2)
    is_valid = True
    for valuei in value:
        if not isinstance(valuei, integer_types):
            is_valid = False
            break

    if is_valid:
        singles, doubles = collapse_thru_packs(value)

        out_value = singles
        for double in doubles:
            assert len(double) == 3, double
            sdouble = '%i THRU %i' % (double[0], double[2])
            out_value.append(sdouble)
    else:
        out_value = value

    for i, out_valuei in enumerate(out_value):
        new_string = '%s, ' % out_valuei
        if len(msg2 + new_string) > 70:
            msg += msg2 + '\n'
            msg2 = ' ' * nchars + new_string
        else:
            msg2 += new_string
    return msg + msg2.rstrip(' \n,') + '\n'
