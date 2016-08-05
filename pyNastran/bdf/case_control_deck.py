# pylint: disable=R0904,R0902,C0103
"""
CaseControlDeck parsing and extraction class
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems, itervalues
import sys
import copy

#from pyNastran.bdf import subcase
from pyNastran.bdf.subcase import Subcase, update_param_name
from pyNastran.utils.log import get_logger


    #def hasParameter(self, isubcase, param_name):
        #self.has_parameter(isubcase, param_name)

    #def get_subcase_parameter(self, isubcase, param_name):
    #def has_subcase(self, isubcase):
    #def create_new_subcase(self, isubcase):
    #def delete_subcase(self, isubcase):
    #def copy_subcase(self, i_from_subcase, i_to_subcase, overwrite_subcase=True):
    #def get_subcase_list(self):
    #def get_local_subcase_list(self):
    #def update_solution(self, isubcase, sol):
    #def add_parameter_to_global_subcase(self, param):
    #def add_parameter_to_local_subcase(self, isubcase, param):
    #def finish_subcases(self):
    #def convert_to_sol_200(self, model):

class CaseControlDeck(object):
    """
    CaseControlDeck parsing and extraction class
    """

    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        del state['log']
        return state

    def __init__(self, lines, log=None):
        """
        Parameters
        ----------
        lines : List[str]
            list of lines that represent the case control deck
            ending with BEGIN BULK
        log : log()
            a :mod: `logging` object
        """
        # pulls the logger from the BDF object
        self.log = get_logger(log, "debug")
        self.debug = False

        self.sol_200_map = {
            'STATICS' : 101,
            'STATIC' : 101,

            'MODES' : 103,
            'MODE' : 103,

            'BUCK' : 105,
            'BUCKLING' : 105,

            'DFREQ' : 108,
            'MFREQ' : 111,
            'SAERO' : 144,

            'FLUTTER' : 145,
            'FLUT' : 145,

            'DIVERGE' : 144,
            'DIVERG' : 145,

            # 'HEAT' : ,
            # 'STRUCTURE' : ,
            'NLSTATICS' : 400,
            'LNSTATICS' : 400,
            'MTRAN' : 112,
            'DCEIG' : 107,
        }
        # 'HEAT', 'ANALYSIS', 'MFREQ', 'STATICS', 'MODES', 'DFREQ',
        # 'MTRAN', 'BUCK', 'MCEIG', 'DCEIG', 'SAERO', 'NLSTATIC', 'NLSTAT',
        # 'STATIC', 'MTRANS', 'MODE', 'FLUTTER', 'DIVERG', 'NLTRAN', 'FLUT',
        #self.debug = True

        #: stores a single copy of 'BEGIN BULK' or 'BEGIN SUPER'
        self.reject_lines = []
        self.begin_bulk = ['BEGIN', 'BULK']

        # allows BEGIN BULK to be turned off
        self.write_begin_bulk = True
        self._begin_count = 0

        self.lines = lines
        self.subcases = {0: Subcase(id=0)}
        try:
            self._read(self.lines)
        except:
            print('\n'.join(self.lines))
            raise

    def suppress_output(self):
        """
        Replaces F06 printing with OP2 printing

        Converts:
            STRESS(PRINT,SORT1,REAL)
            FORCE(PRINT,PLOT,SORT1,REAL)

        to:
            STRESS(PLOT,SORT1,REAL)
            FORCE(PLOT,SORT1,REAL)

        .. warning:: most case control types are not supported
        """
        for isubcase, subcase in iteritems(self.subcases):
            # if isubcase == 0:
                # continue
            subcase.suppress_output()

    def has_parameter(self, isubcase, *param_names):
        """
        Checks to see if a parameter (e.g. STRESS) is defined in a certain
        subcase ID.

        Parameters
        ----------
        isubcase : int
            the subcase ID to check
        param_names : List[str]
            the parameter name to look for
        """
        if self.has_subcase(isubcase):
            return any(self.subcases[isubcase].has_parameter(*param_names))

    def get_subcase_parameter(self, isubcase, param_name):
        """
        Get the [value, options] of a subcase's parameter.  For example, for
        STRESS(PLOT,POST)=ALL, param_name=STRESS, value=ALL, options=['PLOT',
        'POST']


        Parameters
        ----------
        isubcase : int
            the subcase ID to check
        param_name : str
            the parameter name to look for
        """
        if self.has_subcase(isubcase):
            return self.subcases[isubcase].get_parameter(param_name.upper())
        msg = ('isubcase=%r does not exist...subcases=%s'
               % (isubcase, str(sorted(self.subcases.keys()))))
        raise RuntimeError(msg)

    def has_subcase(self, isubcase):
        """
        Checks to see if a subcase exists.

        Parameters
        ----------
        isubcase : int
            the subcase ID

        Returns
        -------
        val : bool
            does_subcase_exist (type = bool)
        """
        if isubcase in self.subcases:
            return True
        return False

    def create_new_subcase(self, isubcase):
        """
        Method create_new_subcase:

        Parameters
        ----------
        isubcase : int
            the subcase ID

        Returns
        -------
        subcase : Subcase()
            the new subcase

        Warning
        -------
        - be careful you dont add data to the global subcase
          after running this...is this True???
        """
        #print("creating subcase=%s" % isubcase)
        if self.has_subcase(isubcase):
            sys.stderr.write('subcase=%s already exists...skipping\n' %
                             isubcase)
            return self.subcases[isubcase]
        subcase = self.copy_subcase(i_from_subcase=0, i_to_subcase=isubcase,
                                    overwrite_subcase=True)
        #self.subcases[isubcase] = Subcase(id=isubcase)
        return subcase

    def delete_subcase(self, isubcase):
        """
        Deletes a subcase.

        :param isubcase: the Subcase to delete
        :type isubcase: int
        """
        if not self.has_subcase(isubcase):
            sys.stderr.write('subcase %s doesnt exist...skipping\n' % isubcase)
        del self.subcases[isubcase]

    def copy_subcase(self, i_from_subcase, i_to_subcase, overwrite_subcase=True):
        """
        Overwrites the parameters from one subcase to another.

        Parameters
        ----------
        i_from_subcase : int
            the Subcase to pull the data from
        i_to_subcase : int
            the Subcase to map the data to
        overwrite_subcase : bool; default=True
            NULLs i_to_subcase before copying i_from_subcase

        Returns
        -------
        subcase : Subcase()
            the new subcase
        """
        #print("copying subcase from=%s to=%s overwrite=%s" % (i_from_subcase, i_to_subcase, overwrite_subcase))
        if not self.has_subcase(i_from_subcase):
            msg = 'i_from_subcase=%r does not exist...subcases=%s' % (i_from_subcase, str(sorted(self.subcases.keys())))
            raise RuntimeError(msg)
        if overwrite_subcase:
            subcase_from = self.subcases[i_from_subcase]
            subcase_to = copy.deepcopy(subcase_from)
            subcase_to.id = i_to_subcase
            #for key, param in sorted(iteritems(subcase_from.params)):
                #print("going to copy key=%s param=%s" % (key, param))
            self.subcases[i_to_subcase] = subcase_to
        else:
            if not self.has_subcase(i_to_subcase):
                msg = ('i_from_subcase=%r does not exist...subcases=%s'
                       % (i_to_subcase, str(sorted(self.subcases.keys()))))
                raise RuntimeError(msg)
            subcase_to = self.subcases[i_to_subcase]

            for key, param in sorted(iteritems(subcase_to)):
                #print('copying key=%s param=%s' % (key, param))
                if key == 'BEGIN':
                    pass
                subcase_to[key] = copy.deepcopy(param)
        return subcase_to

    def get_subcase_list(self):
        """
        Gets the list of subcases including the global subcase ID (0)
        """
        return sorted(self.subcases.keys())

    def get_local_subcase_list(self):
        """
        Gets the list of subcases that aren't the global subcase ID
        """
        id_list = [id for id in self.subcases if id != 0]  # skip the global
        return sorted(id_list)

    def update_solution(self, isubcase, sol):
        """
        sol = STATICS, FLUTTER, MODAL, etc.

        Parameters
        ----------
        isubcase : int
            the subcase ID to update
        sol : str
            the solution type to change the solution to

        >>> print(bdf.case_control)
        SUBCASE 1
            DISP = ALL

        >>> bdf.case_control.update_solution(1, 'FLUTTER')
        >>> print(bdf.case_control)
        SUBCASE 1
            ANALYSIS FLUTTER
            DISP = ALL
        >>>
        """
        self.add_parameter_to_local_subcase(isubcase, 'ANALYSIS %s' % sol)

    def add_parameter_to_global_subcase(self, param):
        """
        Takes in a single-lined string and adds it to the global subcase.

        Parameters
        ----------
        param : str
            the variable to add

        .. note:: dont worry about overbounding the line

        >>> bdf = BDF()
        >>> bdf.read_bdf(bdf_filename)
        >>> bdf.case_control.add_parameter_to_global_subcase('DISP=ALL')
        >>> print(bdf.case_control)
        TITLE = DUMMY LINE
        DISP = ALL
        """
        (j, key, value, options, paramType) = self._parse_data_from_user(param)
        subcase_list = self.get_subcase_list()
        for isubcase in subcase_list:
            self._add_parameter_to_subcase(key, value, options, paramType,
                                           isubcase)

    def add_parameter_to_local_subcase(self, isubcase, param):
        """
        Takes in a single-lined string and adds it to a single Subcase.

        Parameters
        ----------
        isubcase : int
            the subcase ID to add
        param_name : List[str]
            the parameter name to add

        .. note::  dont worry about overbounding the line

        >>> bdf = BDF()
        >>> bdf.read_bdf(bdf_filename)
        >>> bdf.case_control.add_parameter_to_local_subcase(1, 'DISP=ALL')
        >>> print(bdf.case_control)
        TITLE = DUMMY LINE
        SUBCASE 1
            DISP = ALL
        >>>
        """
        (j, key, value, options, param_type) = self._parse_data_from_user(param)
        self._add_parameter_to_subcase(key, value, options, param_type,
                                       isubcase)

    def _parse_data_from_user(self, param):
        """
        Parses a case control line

        Parameters
        ----------
        param : str
            the variable to add
        """
        if '\n' in param or '\r' in param or '\t' in param:
            msg = 'doesnt support embedded endline/tab characters\n'
            msg += '  param = %r' % param
            raise SyntaxError(msg)
        #self.read([param])
        lines = _clean_lines(self, [param])
        (j, key, value, options, param_type) = self._parse_entry(lines)
        return (j, key, value, options, param_type)

    def _read(self, lines):
        """
        Reads the case control deck

        .. note::    supports comment lines
        .. warning:: doesnt check for 72 character width lines, but will
                     follow that when it's written out
        """
        isubcase = 0
        lines = _clean_lines(self, lines)
        self.output_lines = []
        i = 0
        is_output_lines = False
        while i < len(lines):
            line = lines[i] #[:72]
            #comment = lines[i][72:]

            lines2 = [line]
            while ',' in lines[i][-1]:
                i += 1
                lines2.append(lines[i])
                #comment = lines[i][72:]
            (j, key, value, options, param_type) = self._parse_entry(lines2)
            i += 1

            line_upper = line.upper()
            if key == 'BEGIN':
                if 'BULK' not in line_upper and 'SUPER' not in line_upper:
                    raise NotImplementedError('line=%r' % line)
                if self._begin_count == 1:
                    raise NotImplementedError('multiple BEGIN lines are defined...')
                self.begin_bulk = [key, value]
                self._begin_count += 1
                continue
            elif line_upper.startswith('OUTPUT'):
                is_output_lines = True
                #output_line = '%s(%s) = %s\n' % (key, options, value)
                key = 'OUTPUT'

                # OUTPUT(POST) -> POST
                post = line_upper.split('OUTPUT')[1].strip('( )')
                options = [post]
                value = None
                param_type = 'STRESS-type'

                isubcase = self._add_parameter_to_subcase(key, value, options,
                                                          param_type, isubcase)
                self.output_lines.append(line)
                continue
            #print("key=%-12r icase=%i value=%r options=%r param_type=%r" %(key,
            #    isubcase, value, options, param_type))
            isubcase = self._add_parameter_to_subcase(key, value, options,
                                                      param_type, isubcase)

        #print(str(self))
        self.finish_subcases()

    def _parse_entry(self, lines):
        r"""
        Internal method for parsing a card of the case control deck

        Parses a single case control deck card into 4 sections

        1.  paramName - obvious
        2.  Value     - still kind of obvious
        3.  options   - rarely used data
        4.  paramType - STRESS-type, SUBCASE-type, PARAM-type, SET-type, BEGIN_BULK-type

        It's easier with examples:

        paramType = SUBCASE-type
          SUBCASE 1              ->   paramName=SUBCASE  value=1            options=[]
        paramType = STRESS-type
          STRESS       = ALL     ->   paramName=STRESS    value=ALL         options=[]
          STRAIN(PLOT) = 5       ->   paramName=STRAIN    value=5           options=[PLOT]
          TITLE        = stuff   ->   paramName=TITLE     value=stuff       options=[]
        paramType = SET-type
          SET 1 = 10,20,30       ->   paramName=SET       value=[10,20,30]  options = 1
        paramType = BEGIN_BULK-type
          BEGIN BULK             ->   paramName=BEGIN     value=BULK        options = []
        paramType = CSV-type
          PARAM,FIXEDB,-1        ->   paramName=PARAM     value=FIXEDB      options = [-1]

        The paramType is the "macro" form of the data (similar to integer, float, string).
        The value is generally whats on the RHS of the equals sign (assuming it's there).
        Options are modifiers on the data.  Form things like the PARAM card or the SET card
        they arent as clear, but the paramType lets the program know how to format it
        when writing it out.

        Parameters
        ----------
        lines : List[str, str, ...]
            list of lines

        Returns
        -------
        paramName : str
            see brief
        value : List[...]
            see brief
        options : List[str/int/float]
            see brief
        paramType : str/int/float/List
            see brief
        """
        i = 0
        options = []
        value = None
        key = None
        param_type = None

        line = lines[i]
        #print(line)
        #print("*****lines = ", lines)

        equals_count = 0
        for letter in line:
            if letter == '=':
                equals_count += 1
        line_upper = line.upper()

        #print("line_upper = %r" % line)
        #print('  equals_count = %s' % equals_count)
        if line_upper.startswith('SUBCASE'):
            #print("line = %r" % line)
            line2 = line.replace('=', '')
            sline = line2.split()
            if len(sline) != 2:
                msg = "trying to parse %r..." % line
                raise RuntimeError(msg)
            (key, param_type) = sline
            key = key.upper()

            #print("key=%r isubcase=%r" % (key, isubcase))
            value = int(param_type)
            #self.isubcase = int(isubcase)
            param_type = 'SUBCASE-type'
            assert key.upper() == key, key

        elif line_upper.startswith(('LABEL', 'SUBT', 'TITL')):  # SUBTITLE/TITLE
            try:
                eindex = line.index('=')
            except:
                msg = "cannot find an = sign in LABEL/SUBTITLE/TITLE line\n"
                msg += "line = %r" % line_upper.strip()
                raise RuntimeError(msg)

            key = line_upper[0:eindex].strip()
            value = line[eindex + 1:].strip()
            options = []
            param_type = 'STRING-type'
        elif (line_upper.startswith('SET ') or line_upper.startswith('SETMC ')
              ) and equals_count == 1:
            # would have been caught by STRESS-type
            sline = line_upper.split('=')
            assert len(sline) == 2, sline

            key, value = sline
            try:
                (key, ID) = key.split()
            except:
                raise RuntimeError(key)
            key = key + ' ' + ID

            assert key.upper() == key, key
            options = int(ID)

            if self.debug:
                self.log.debug('SET-type key=%r ID=%r' % (key, ID))
            fivalues = value.rstrip(' ,').split(',')  # float/int values

            #: .. todo:: should be more efficient multiline reader...
            # read more lines....
            if line[-1].strip() == ',':
                i += 1
                #print("rawSETLine = %r" % (lines[i]))
                while 1:
                    if ',' == lines[i].strip()[-1]:
                        fivalues += lines[i][:-1].split(',')
                    else:  # last case
                        fivalues += lines[i].split(',')
                        #print("fivalues last = i=%s %r" % (i, lines[i]))
                        i += 1
                        break
                    i += 1
            #print("len(fivalues) = ",len(fivalues))
            value = fivalues
            param_type = 'SET-type'

        elif equals_count == 1:  # STRESS
            if '=' in line:
                (key, value) = line_upper.strip().split('=')
            else:
                msg = 'expected item of form "name = value"   line=%r' % line.strip()
                raise RuntimeError(msg)

            key = key.strip().upper()
            value = value.strip()
            if self.debug:
                self.log.debug("key=%r value=%r" % (key, value))
            param_type = 'STRESS-type'
            assert key.upper() == key, key

            if '(' in key:  # comma may be in line - STRESS-type
                #paramType = 'STRESS-type'
                sline = key.strip(')').split('(')
                key = sline[0]
                options = sline[1].split(',')

                # handle TEMPERATURE(INITIAL) and TEMPERATURE(LOAD) cards
                if key in ['TEMPERATURE', 'TEMP']:
                    option = options[0]
                    if option == '':
                        option = 'BOTH'
                    key = 'TEMPERATURE'
                    options = [option]
                #print("key=%r options=%s" %(key,options))

            #elif ' ' in key and ',' in value:  # SET-type
                #(key, ID) = key.split()
                #key = key + ' ' + ID

                #if self.debug:
                    #self.log.debug('SET-type key=%r ID=%r' % (key, ID))
                #fivalues = value.rstrip(' ,').split(',')  # float/int values

                ##: .. todo:: should be more efficient multiline reader...
                ## read more lines....
                #if line[-1].strip() == ',':
                    #i += 1
                    ##print("rawSETLine = %r" % (lines[i]))
                    #while 1:
                        #if ',' == lines[i].strip()[-1]:
                            #fivalues += lines[i][:-1].split(',')
                        #else:  # last case
                            #fivalues += lines[i].split(',')
                            ##print("fivalues last = i=%s %r" % (i, lines[i]))
                            #i += 1
                            #break
                        #i += 1
                ##print("len(fivalues) = ",len(fivalues))
                #value = fivalues

                #options = ID  # needed a place to put it...
                #param_type = 'SET-type'
            elif ',' in value:  # STRESS-type; special TITLE = stuffA,stuffB
                #print('A ??? line = ',line)
                #raise RuntimeError(line)
                pass
            else:  # STRESS-type; TITLE = stuff
                #print('B ??? line = ',line)
                if key in ['TEMPERATURE', 'TEMP']:
                    assert len(options) == 0, options
                    key = 'TEMPERATURE'
                    options = ['BOTH']

            key = update_param_name(key.strip().upper())
            verify_card(key, value, options, line)
            assert key.upper() == key, key
        elif equals_count > 2 and '(' in line and 'FLSPOUT' not in line:
            #GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES
            #print('****', lines)
            assert len(lines) == 1, lines
            line = lines[0]
            try:
                key, value_options = line.split('(', 1)
            except ValueError:
                msg = 'Expected a "(", but did not find one.\n'
                msg += 'Looking for something of the form:\n'
                msg += '   GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES\n'
                msg += '%r' % line
                raise ValueError(msg)

            try:
                options_paren, value = value_options.rsplit('=', 1)
            except ValueError:
                msg = 'Expected a "=", but did not find one.\n'
                msg += 'Looking for something of the form:\n'
                msg += '   GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES\n'
                msg += 'value_options=%r\n' % value_options
                msg += '%r' % line
                raise ValueError(msg)
            options_paren = options_paren.strip()

            value = value.strip()
            if value.isdigit():
                value = int(value)
            if not options_paren.endswith(')'):
                raise RuntimeError(line)
            str_options = options_paren[:-1]

            if '(' in str_options:
                options_start = []
                options_end = []
                icomma = str_options.index(',')
                iparen = str_options.index('(')
                #print('icomma=%s iparen=%s' % (icomma, iparen))
                while icomma < iparen:
                    base, str_options = str_options.split(',', 1)
                    str_options = str_options.strip()
                    icomma = str_options.index(',')
                    iparen = str_options.index('(')
                    options_start.append(base.strip())
                    #print('  icomma=%s iparen=%s' % (icomma, iparen))
                    #print('  options_start=%s' % options_start)

                icomma = str_options.rindex(',')
                iparen = str_options.rindex(')')
                #print('icomma=%s iparen=%s' % (icomma, iparen))
                while icomma > iparen:
                    str_options, end = str_options.rsplit(')', 1)
                    str_options = str_options.strip() + ')'
                    #print('  str_options = %r' % str_options)
                    icomma = str_options.rindex(',')
                    iparen = str_options.rindex(')')
                    options_end.append(end.strip(' ,'))
                    #print('  icomma=%s iparen=%s' % (icomma, iparen))
                    #print('  options_end=%s' % options_end[::-1])

                #print()
                #print('options_start=%s' % options_start)
                #print('options_end=%s' % options_end)
                #print('leftover = %r' % str_options)
                options = options_start + [str_options] + options_end[::-1]

            else:
                options = str_options.split(',')
            param_type = 'STRESS-type'
            key = key.upper()
            #print('options =', options)
            #asdf
        elif line_upper.startswith('BEGIN'):  # begin bulk
            try:
                (key, value) = line_upper.split(' ')
            except:
                msg = 'excepted "BEGIN BULK" found=%r' % line
                raise RuntimeError(msg)
            key = key.upper()
            param_type = 'BEGIN_BULK-type'
            assert key.upper() == key, key
        elif 'PARAM' in line_upper:  # param
            if ',' in line_upper:
                sline = line_upper.split(',')
            elif '\t' in line_upper:
                sline = line_upper.split('\t')
            else:
                raise SyntaxError("trying to parse %r..." % line)

            if len(sline) != 3:
                raise SyntaxError("trying to parse %r..." % line)
            (key, value, options) = sline
            param_type = 'CSV-type'
            assert key.upper() == key, key
        elif ' ' not in line:
            key = line.strip().upper()
            value = line.strip()
            options = None
            param_type = 'KEY-type'
            assert key.upper() == key, key
        else:
            msg = 'generic catch all...line=%r' % line
            key = ''
            value = line
            options = None
            param_type = 'KEY-type'
            assert key.upper() == key, key
        i += 1
        assert key.upper() == key, 'key=%s param_type=%s' % (key, param_type)

        return (i, key, value, options, param_type)

    def finish_subcases(self):
        """
        Removes any unwanted data in the subcase...specifically the SUBCASE
        data member.  Otherwise it will print out after a key like stress.
        """
        for subcase in itervalues(self.subcases):
            subcase.finish_subcase()

    def convert_to_sol_200(self, model):
        """
        Takes a case control deck and changes it from a SOL xxx to a SOL 200

        :param model: BDF()
            the BDF object

        .. todo:: not done...
        """
        analysis = model.rsolmap_toStr[model.sol]
        model.sol = 200

        subcase0 = self.subcases[0]
        subcase0.add_parameter_to_global_subcase('ANALYSIS', analysis)
        #subcase.add_parameter_to_global_subcase('DESSUB', dessub)

    def _add_parameter_to_subcase(self, key, value, options, param_type, isubcase):
        """
        Internal method
        """
        if self.debug:
            a = 'key=%r' % key
            b = 'value=%r' % value
            c = 'options=%r' % options
            d = 'param_type=%r' % param_type
            msg = "_adding isubcase=%s %-12s %-12s %-12s %-12s" % (isubcase, a,
                                                                   b, c, d)
            self.log.debug(msg)

        if key == 'SUBCASE':
            assert value not in self.subcases, 'key=%s value=%s already exists' % (key, value)
            assert isinstance(value, int)
            isubcase = value
            #print("value=", value)
            self.copy_subcase(i_from_subcase=0, i_to_subcase=isubcase,
                              overwrite_subcase=True)
            if self.debug:
                msg = "copied subcase i_from_subcase=%r to i_to_subcase=%r" % (0, isubcase)
                self.log.debug(msg)
        elif isubcase not in self.subcases:  # initialize new subcase
            #self.isubcase += 1 # is handled in the read code
            msg = 'isubcase=%r is not a valid subcase...subcases=%s' % (
                isubcase, str(sorted(self.subcases.keys())))
            raise RuntimeError(msg)

        subcase = self.subcases[isubcase]
        subcase._add_data(key, value, options, param_type)

        #print("\n%s\n" % (self.subcases[isubcase]))
        return isubcase

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        for isubcase, subcase in sorted(iteritems(self.subcases)):
            subcase.cross_reference(model)

    def get_op2_data(self):
        """
        Gets the relevant op2 parameters required for a given subcase

        .. todo:: not done...
        """
        cases = {}
        for isubcase, subcase in sorted(iteritems(self.subcases)):
            if isubcase:
                cases[isubcase] = subcase.getOp2Data(self.sol, subcase.solmap_toValue)
        return cases

    def __repr__(self):
        msg = ''
        subcase0 = self.subcases[0]
        for subcase_id, subcase in sorted(iteritems(self.subcases)):
            msg += subcase.write_subcase(subcase0)
        #if len(self.subcases) == 1:
            #msg += 'BEGIN BULK\n'

        if self.output_lines:
            msg += '\n'.join(self.output_lines) + '\n'
        msg += '\n'.join(self.reject_lines)
        if self.write_begin_bulk:
            msg += ' '.join(self.begin_bulk) + '\n'
        return msg


def verify_card(key, value, options, line):
    if key in ['AUXMODEL', 'BC', 'BCHANGE', 'BCMOVE', 'CAMPBELL', 'CLOAD',
               'CMETHOD', 'CSSCHD', 'DEACTEL', 'DEFORM', 'DESGLB', 'DESSUB',
               'DIVERG', 'DLOAD', 'DRSPAN', 'FMETHOD', 'FREQUENCY', 'GUST',
               'HADAPART', 'LINE', 'LOAD', 'LOADSET', 'MAXLINES', 'MCHSTAT',
               'MFLUID', 'MODES', 'MODTRAK', 'MPC', 'NLHARM',]:
        try:
            value2 = int(value)
        except:
            print('line=%r is invalid; value=%r' % (line, value))
            raise
        assert value2 > 0, 'line=%r is invalid; value=%r must be greater than 0.' % (line, value2)

def verify_card2(key, value, options, line):
    """
    Make sure there are no obvious errors
    """
    # this is purposely made overly strict to catch all the cases
    int_cards = [
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
        '', '', '', '', '', '', '', '', '', '',
        '',
    ]

    # these may only be integers
    #print("key =", key)

    pass_headers = [
        'SUBTITLE', 'TITLE',
        'A2GG', 'M2GG', 'K2GG',
        'K2PP', 'M2PP',
        'K42GG',

        'XMIN', 'XMAX', 'XTITLE', 'XPAPE', 'XPAPER', 'XAXIS', 'XGRID', 'XGRID LINES', 'XLOG',
        'YMIN', 'YMAX', 'YTITLE', 'YPAPE', 'YPAPER', 'YAXIS', 'YGRID', 'YGRID LINES', 'YLOG',
        'XTMIN', 'XTMAX', 'XTGRID', 'XTTITLE', 'XTAXIS', 'XTGRID LINES', 'XTLOG',
        'YTMIN', 'YTMAX', 'YTGRID', 'YTTITLE', 'YTAXIS', 'YTGRID LINES', 'YTLOG',
        'XBMIN', 'XBMAX', 'XBGRID', 'XBAXIS', 'XBGRID LINES', 'XBTITLE', 'XBLOG',
        'YBMIN', 'YBMAX', 'YBGRID', 'YBAXIS', 'YBGRID LINES', 'YBTITLE', 'YBLOG',

        'RIGHT TICS', 'UPPER TICS',
        'TRIGHT TICS',
        'BRIGHT TICS',

        'PLOTTER', 'XYPLOT',

        'PTITLE',
        'HOUTPUT', 'PLOTID', '', '', '', '', '',
        'AXISYMMETRIC', 'CURVELINESYMBOL', 'CURVELINESYMB', 'AECONFIG',
        'B2GG', 'B2PP', 'AESYMXZ', 'TEMP', 'DSAPRT', 'MEFFMASS',
        'MAXMIN', 'RESVEC', 'MODESELECT', 'RIGID', 'TCURVE',
        'SUPER', 'MAXI DEFO', 'P2G',
        'EXTSEOUT', 'FLSTCNT PREFDB', 'AESYMXY',
        'DSYM', '', '', ''
    ]
    if key in ['BCONTACT', 'CURVELINESYMBOL']:
        try:
            value2 = int(value)
        except:
            print('line=%r is invalid; value=%r' % (line, value))
            raise

    # these may only be integers greater than 0
    elif key in int_cards:
        try:
            value2 = int(value)
        except:
            print('line=%r is invalid; value=%r' % (line, value))
            raise
        assert value2 > 0, 'line=%r is invalid; value=%r must be greater than 0.' % (line, value2)

    # these may have a value of all/none/integer, nothing else
    # except commas are allowed
    # 'DISP=ALL', 'DISP=NONE', 'DISP=1', 'DISP=1,2'
    elif key in ['STRESS', 'STRAIN', 'SPCFORCES', 'DISPLACEMENT', 'MPCFORCES', 'SVECTOR',
                 'VELOCITY', 'ACCELERATION', 'FORCE', 'ESE', 'OLOAD', 'SEALL', 'GPFORCE',
                 'GPSTRESS', 'GPSTRAIN', 'FLUX', 'AEROF', 'THERMAL', 'STRFIELD',
                 'NOUTPUT', 'SEDV', 'APRES', 'HTFLOW', 'NLSTRESS', 'GPKE',
                 'SACCELERATION', 'SDISPLACEMENT', 'SEMG', 'HARMONICS', 'PRESSURE', 'VUGRID',
                 'ELSUM', 'SVELOCITY', 'STRFIELD REAL', 'SENSITY', 'MONITOR',
                 'NLLOAD', 'GPSDCON', 'BOUTPUT', '', '', '']:
        if value not in ['ALL', 'NONE']:
            if ',' in value:
                sline = value.split(',')
                for spot in sline:
                    try:
                        value2 = int(spot)
                    except:
                        print('line=%r is invalid; value=%r' % (line, spot))
                        raise
            else:
                try:
                    value2 = int(value)
                except:
                    print('line=%r is invalid; value=%r' % (line, value))
                    raise
                assert value2 > 0, 'line=%r is invalid; value=%r must be greater than 0.' % (line, value2)
    elif key in ['ECHO']:
        #assert value in ['NONE','BOTH','UNSORT','SORT', 'NOSORT', 'PUNCH', ''], 'line=%r is invalid; value=%r.' % (line, value)
        pass
    elif key in ['CSCALE', 'SUBSEQ', 'SYMSEQ', 'DEFORMATION SCALE', '', '']:
        # floats
        pass
    elif 'SET' in key:
        pass

    # weird cards
    elif key in pass_headers:
        pass
    elif key == 'ANALYSIS':
        assert value in ['HEAT', 'ANALYSIS', 'MFREQ', 'STATICS', 'MODES', 'DFREQ',
                         'MTRAN', 'BUCK', 'MCEIG', 'DCEIG', 'SAERO', 'NLSTATIC', 'NLSTAT',
                         'STATIC', 'MTRANS', 'MODE', 'FLUTTER', 'DIVERG', 'NLTRAN',
                         '', '', '', '', ''], 'line=%r is invalid; value=%r' % (line, value)
    elif key == 'AUTOSPC':
        assert value in ['YES'], 'line=%r is invalid; value=%r' % (line, value)
    else:
        raise NotImplementedError('key=%r line=%r' % (key, line))


def _clean_lines(case_control, lines):
    """
    Removes comment characters defined by a *$*.

    :param case_control:  the CaseControlDeck object
    :param lines: the lines to clean.
    """
    lines2 = []
    for line in lines:
        line = line.strip(' \n\r').split('$')[0].rstrip()
        if line:
            lines2.append(line)

    lines3 = []  # TODO: line, comment
    lines_pack = []
    for line in lines2:
        #print(line)
        if len(lines_pack) == 0:
            #print('0--', line)
            lines_pack.append(line)
            if not line.endswith(','):
                #print('next...')
                lines3.append(lines_pack)
                lines_pack = []
        elif line.endswith(','):
            #print('C--', line)
            lines_pack.append(line)
        else:
            if lines_pack[-1][-1] == ',':  # continued
                #print('xx--', line)
                lines_pack.append(line)
                lines3.append(lines_pack)
                #print('pack =', lines_pack)
                lines_pack = []
            else:  # new card
                #print('new--', line)
                lines3.append(lines_pack)
                lines_pack = [line]
    return [''.join(pack) for pack in lines3]

def main():
    lines = [
        'SUBCASE 1',
        '    ACCELERATION(PLOT,PRINT,PHASE) = ALL',
        '    DISPLACEMENT(PLOT,PRINT,PHASE) = ALL',
        '    DLOAD = 32',
        '    M2GG = 111',
        '    SET 88  = 5, 6, 7, 8, 9, 10 THRU 55 EXCEPT 15, 16, 77, 78, 79, '
        '100 THRU 300',
        '    SET 99  = 1 THRU 10',
        '    SET 105 = 1.009, 10.2, 13.4, 14.0, 15.0',
        '    SET 111 = MAAX1,MAAX2',
        '    SET 1001 = 101/T1, 501/T3, 991/R3',
        '    SET = ALL',
        '    SPC = 42',
        '    TSTEPNL = 22',
        '    VELOCITY(PLOT,PRINT,PHASE) = ALL',
        'BEGIN BULK',
    ]
    deck = CaseControlDeck(lines)
    #deck.create_new_subcase(2)
    #deck.add_parameter_to_local_subcase(0, 'SET 2 = 11,12,13,14,15,16,17,18,'
    #   '19,20,21,22,23,24,25,26,'
    #   '1000000000000000000000000000000000000000000000000000000,33')
    str(deck)
    #print('%s\n\n' % deck)

    #deck2 = CaseControlDeck(['ACCELERATION(PLOT,PRINT,PHASE) = ALL',
    #                         'DISPLACEMENT(PLOT,PRINT,PHASE) = ALL',
    #                         'BEGIN BULK'])
    #print('\n\n%s' % deck2)

if __name__ == '__main__':  # pragma: no cover
    main()
