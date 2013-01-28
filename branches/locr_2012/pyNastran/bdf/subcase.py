# pylint: disable=C0103,R0201
"""
Subcase creation/extraction class
"""
#import sys


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
        #print "\n***adding subcase %s***" %(self.id)

    def get_stress_code(self, key, options, value):
        """
        Method get_stress_code:
        @note the individual element must take the stress_code and reduce it to
         what the element can return.  For example, for an isotropic CQUAD4
         the fiber field doesnt mean anything.

        BAR - no von mises/fiber
        ISOTROPIC - no fiber

        @todo how does the MATERIAL bit get turned on?  I'm assuming it's
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

    def get_format_code(self, options, value):
        """
        Gets the format code that will be used by the op2 based on
        the options
        @param self the Subcase object
        @param options the options for a parameter
        @param value the value of the parameter
        @todo not done...only supports REAL, IMAG, PHASE, not RANDOM
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
        @param self the Subcase object
        @param options the options for a parameter
        @param value the value of the parameter
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
        @param self the Subcase object
        @param options the options for a parameter
        @param value the value of the parameter
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
        8 - post-buckling (maybe 7 depending on NLPARM???)

        # not important
        3/4 - differential stiffness (obsolete)
        11  - old geometric nonlinear statics
        12  - contran (???)
        @todo verify
        """
        codes = {
            101: 1,  # staics
            103: 2,  # modes
            105: 7,  # pre-buckling
            106: 10,  # nonlinear statics
            107: 9,  # complex eigenvalues
            108: 5,  # frequency
            111: 5,
            112: 6,
            114: 1,
            115: 2,
            116: 7,
            118: 5,
            129: 6,  # nonlinear
            144: 1,  # static aero
            145: 1,
            146: 1,  # flutter
            153: 10,
            159: 6,  # transient thermal
        }
        #print "sol=%s" %(sol)
        approach_code = codes[sol]
        #print 'approach_code = ',approach_code
        return approach_code

    def get_table_code(self, sol, table_name, options):
        """
        Gets the table code of a given parameter.  For example, the
        DISPLACMENT(PLOT,POST)=ALL makes an OUGV1 table and stores the
        displacement.  This has an OP2 table code of 1, unless you're running a
        modal solution, in which case it makes an OUGV1 table of eigenvectors
        and has a table code of 7.

        @param self the Subcase object
        @param options the options for a parameter
        @param value the value of the parameter
        """
        if table_name in ['VECTOR', 'PRESSURE']:
            table_name = 'DISPLACEMENT'  # equivalent tables...

        key = (sol, table_name)
        tables = {  # SOL, table_name      table_code
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
        print("key=%s" % (str(key)))
        if key not in tables:
            raise KeyError(key)
        table_code = tables[key]
        return table_code

    def has_parameter(self, paramName):
        """
        Checks to see if a parameter name is in the subcase.
        @param self the Subcase object
        """
        if paramName in self.params:
            return True
        return False

    def get_parameter(self, paramName):
        """
        Gets the [value, options] for a subcase.
        @param self the Subcase object
        """
        paramName = self.update_param_name(paramName)
        if paramName not in self.params:
            raise KeyError('%s doesnt exist in subcase=%s in the case '
                           'control deck.' % (paramName, self.id))
        return self.params[paramName][0:2]

    def update_param_name(self, param_name):
        """
        Takes an abbreviated name and expands it so the user can type DISP or
        DISPLACEMENT and get the same answer
        @param self
          the Subcase object
        @param param_name
          the parameter name to be standardized (e.g. DISP vs. DIPLACEMENT)
        @todo not a complete list
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
        elif param_name.startswith('FORC'):
            param_name = 'FORCE'
        elif param_name.startswith('FREQ'):
            param_name = 'FREQUENCY'
        elif param_name.startswith('GPFO'):
            param_name = 'GPFORCE'
        elif param_name == 'GPST':
            raise SyntaxError('invalid GPST stress/strain')
        elif param_name.startswith('METH'):
            param_name = 'METHOD'
        elif param_name.startswith('MPCF'):
            param_name = 'MPCFORCES'
        elif param_name.startswith('OLOA'):
            param_name = 'OLOAD'
        elif param_name.startswith('PRES'):
            param_name = 'PRESSURE'
        elif param_name.startswith('SPCF'):
            param_name = 'SPCFORCES'
        elif param_name.startswith('STRA'):
            param_name = 'STRAIN'
        elif param_name.startswith('STRE'):
            param_name = 'STRESS'
        elif param_name.startswith('SUPO'):
            param_name = 'SUPORT1'
        elif param_name.startswith('SVEC'):
            param_name = 'SVECTOR'
        elif param_name.startswith('THER'):
            param_name = 'THERMAL'
        elif param_name.startswith('VECT'):
            param_name = 'VECTOR'
        elif param_name.startswith('VELO'):
            param_name = 'VELOCITY'

        #elif param_name.startswith('DFRE'):  param_name = 'D'

        # handled in caseControlDeck.py
        #elif param_name.startswith('TEMP'):  param_name = 'TEMPERATURE'
        #print '*param_name = ',param_name
        return param_name

    def _add_data(self, key, value, options, param_type):
        key = self.update_param_name(key)
        #print("adding isubcase=%s key=|%s| value=|%s| options=|%s| "
        #      "param_type=%s" %(self.id, key, value, options, param_type))
        if isinstance(value, unicode) and value.isdigit():
            value = int(value)

        (key, value, options) = self._simplify_data(key, value, options,
                                                    param_type)
        self.params[key] = [value, options, param_type]

    def _simplify_data(self, key, value, options, param_type):
        if param_type == 'SET-type':
            #print("adding isubcase=%s key=|%s| value=|%s| options=|%s| "
            #      "param_type=%s" %(self.id, key, value, options, param_type))
            values2 = []
            for (i, ivalue) in enumerate(value):
                ivalue = ivalue.strip()
                if ivalue.isdigit():
                    values2.append(int(ivalue))
                else:
                    if value is 'EXCLUDE':
                        msg = ('EXCLUDE is not supported on CaseControlDeck '
                               'SET card\n')
                        raise RuntimeError(msg)
                    values2.append(ivalue)

            # TODO expand values with THRU and EXCLUDE
            # TODO sort values
            # TODO collapse values when printing

            #print "values2 = ",values2
            options = int(options)
            return (key, values2, options)

        elif param_type == 'CSV-type':
            #print("adding isubcase=%s key=|%s| value=|%s| options=|%s| "
            #      "param_type=%s" %(self.id, key, value, options, param_type))
            if value.isdigit():  # PARAM,DBFIXED,-1
                value = value
        else:
            #a = 'key=|%s|'       % key
            #b = 'value=|%s|'     % value
            #c = 'options=|%s|'   % options
            #d = 'param_type=|%s|' % param_type
            #print("_adding isubcase=%s %-18s %-12s %-12s %-12s" %(self.id, a,
            #                                                      b, c, d))
            if isinstance(value, int) or value is None:
                pass
            elif value.isdigit():  # STRESS = ALL
                value = value
            #else: pass
        return (key, value, options)

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
            print("***value=%s sol=%s" % (value, sol))
        else:  # leaves SOL the same
            sol = self.sol

        if sol in self.solCodeMap:  # reduces SOL 144 to SOL 101
            sol = self.solCodeMap[sol]

        for (key, param) in self.params.iteritems():
            key = key.upper()
            (value, options, paramType) = param
            #msg = ("  -key=|%s| value=|%s| options=%s paramType=|%s|"
            #    % (key, value, options, paramType))

        thermal = 0
        for (key, param) in self.params.iteritems():
            key = key.upper()
            (value, options, paramType) = param
            #msg = ("  *key=|%s| value=|%s| options=%s paramType=|%s|"
            #    % (key, value, options, paramType)
            #print(msg)
            #msg += self.printParam(key, param, printBeginBulk=False)
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

        print("\nThe estimated results...")
        for (key, value) in sorted(op2Params.iteritems()):
            if value is not None:
                print("   key=|%s| value=|%s|" % (key, value))

    def print_param(self, key, param, printBeginBulk=True):
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

        if param_type == 'SUBCASE-type':
            if self.id > 0:
                msg += 'SUBCASE %s\n' % (self.id)
            #else:  global subcase ID=0 and is not printed
            #    pass
        elif param_type == 'KEY-type':
            #print "KEY-TYPE:  |%s|" %(value)
            assert value is not None, param
            msg += spaces + '%s\n' % (value)
        elif param_type == 'STRING-type':
            msg += spaces + '%s = %s\n' % (key, value)
        elif param_type == 'CSV-type':
            msg += spaces + '%s,%s,%s\n' % (key, value, options)
        elif param_type == 'STRESS-type':
            sOptions = ','.join(options)
            #print("sOptions = |%s|" %(sOptions))
            #print("STRESSTYPE key=%s value=%s options=%s"
            #    %(key, value, options))
            if len(sOptions) > 0:
                msg += '%s(%s) = %s\n' % (key, sOptions, value)
            else:
                msg += '%s = %s\n' % (key, value)
            msg = spaces + msg

        elif param_type == 'BEGIN_BULK-type':
            msg += '%s %s\n' % (key, value)
            if 'BEGIN BULK' not in msg:
                msg = spaces + msg
            elif printBeginBulk:
                pass
            else:
                msg = ''

        elif param_type == 'SET-type':
            # TODO collapse data...not written yet
            starter = 'SET %s = ' % (options)
            msg2 = spaces + starter
            nChars = len(msg2)

            i = 0
            while i < len(value):
                #print "type(value[i]) = ",type(value[i])
                newString = '%s, ' % (value[i])
                #print "newString[%i] = |%s|" %(i,newString)
                if len(msg2 + newString) > 70:
                    msg += msg2 + '\n'
                    msg2 = ' ' * nChars + newString
                else:
                    msg2 += newString
                i += 1

            msg += msg2.rstrip(' \n,') + '\n'
        else:
            # SET-type is not supported yet...
            raise NotImplementedError((key, param))

        #print "msg = |%r|" %(msg)
        return msg

    def cross_reference(self, model):
        """
        Method crossReference:
        @note
          this is not integrated and probably never will be as it's not
          really that necessary.  it's only really useful when running an
          analysis
        """
        print("keys = %s" % (sorted(self.params.keys())))
        if 'LOAD' in self.params:
            loadID = self.params['LOAD'][0]
            loadObj = model.loads[loadID]
            loadObj.cross_reference(model)
        if 'SUPORT' in self.params:
            pass
        if 'MPC' in self.params:
            #mpcID = self.params['MPC'][0]
            #mpcObj = model.mpcs[mpcID]
            #mpcObj.cross_reference(model)
            pass
        if 'SPC' in self.params:
            #spcID = self.params['SPC'][0]
            #print "SPC ID = ",spcID
            #spcObj = model.spcObject
            #spcObj.cross_reference(spcID, model)
            pass
        if 'TSTEPNL' in self.params:
            tstepnlID = self.params['TSTEPNL'][0]
            tstepnlObj = model.tstepnl[tstepnlID]
            tstepnlObj.cross_reference(model)
        if 'NLPARM' in self.params:
            nlparmID = self.params['NLPARM'][0]
            nlparmObj = model.nlparms[nlparmID]
            nlparmObj.cross_reference(model)
        if 'TRIM' in self.params:
            trimID = self.params['TRIM'][0]
            trimObj = model.trims[trimID]
            trimObj.cross_reference(model)
        if 'GUST' in self.params:
            gustID = self.params['GUST'][0]
            gustObj = model.gusts[gustID]
            gustObj.cross_reference(model)
        if 'DLOAD' in self.params:  # ???
            pass

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
        internal method to print a subcase
        @param self the subcase object
        @param subcase0 the global subcase
        """
        if self.id == 0:
            msg = str(self)
        else:
            msg = 'SUBCASE %s\n' % (self.id)
            for (key, param) in self.subcase_sorted(self.params.items()):
                if key in subcase0.params and subcase0.params[key] == param:
                    pass  # dont write global subcase parameters
                else:
                    if 'key' == 'BEGIN':
                        continue
                    else:
                        #print "key=%s param=%s" %(key, param)
                        (value, options, paramType) = param
                        #print("  *key=|%s| value=|%s| options=%s "
                        #      "paramType=|%s|" % (key, value, options,
                        msg += self.print_param(key, param,
                                                printBeginBulk=False)
                        #print ""

        ## self.id>0 and 'BEGIN' in self.params used to prevent printing of
        ## BEGIN BULK multiple times
        if self.id > 0 and 'BEGIN' in self.params:
            msg += self.print_param('BEGIN', self.params['BEGIN'])
        return msg

    def subcase_sorted(self, lst):
        """
        Does a "smart" sort on the keys such that SET cards increment in
        numerical order.  Also puts the sets first.
        @param self the Subcase object
        @param lst the list of subcase list objects
        @retval listB the sorted version of listA
        """
        # presort the list to put all the SET cards next to each other
        # instead of listA.sort() as it allows lst to be any iterable
        lst = sorted(lst) 

        i = 0  # needed in case the loop doesn't execute
        iSet = None  # index of first SET card in the deck
        setDict = {}
        listBefore = []
        setKeys = []
        for (i, entry) in enumerate(lst):
            key = entry[0]
            if 'SET' in key[0:3]:
                if key == 'SET':  # handles "SET = ALL"
                    key = 0
                else:  # handles "SET 100 = 1,2,3"
                    sline = key.split(' ')
                    key = int(sline[1])

                # store the integer ID and the SET-type list
                setDict[key] = entry
                setKeys.append(key)
            else:
                # only store the entries before the SET cards
                listBefore.append(entry)
                if iSet:
                    break

        # grab the other entries
        listAfter = lst[i + 1:]

        # write the SET cards in a sorted order
        setList = []
        for key in sorted(setKeys):
            setList.append(setDict[key])

        # combine all the cards
        listB = setList + listBefore + listAfter
        return listB

    def __repr__(self):
        """
        Prints out every entry in the subcase.  Skips parameters already in
        the global subcase.
        """
        #msg = "-------SUBCASE %s-------\n" %(self.id)
        msg = ''
        if self.id > 0:
            msg += 'SUBCASE %s\n' % (self.id)

        for (key, param) in self.subcase_sorted(self.params.items()):
            if 'key' == 'BEGIN':
                continue
            else:
                #print "key=%s param=%s" %(key,param)
                (value, options, paramType) = param
                #print("  ?*key=|%s| value=|%s| options=%s paramType=|%s|"
                #      %(key,value,options,paramType))
                msg += self.print_param(key, param, printBeginBulk=False)
                #print ""

        if self.id > 0 and 'BEGIN' in self.params:  # prevents 2 BEGIN BULKs
            msg += self.print_param('BEGIN', self.params['BEGIN'])
        return msg