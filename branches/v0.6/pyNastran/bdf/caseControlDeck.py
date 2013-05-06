# pylint: disable=R0904,R0902,C0103
"""
CaseControlDeck parsing and extraction class
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
import copy

from pyNastran.bdf import subcase
from pyNastran.bdf.subcase import Subcase, update_param_name
from pyNastran.utils.log import get_logger


class CaseControlDeckDeprecated(object):
    def __init__(self):
        pass

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
    def __init__(self, lines, log=None):
        """
        :param self:  the CaseControlDeck object
        :param lines: list of lines that represent the case control deck
                ending with BEGIN BULK
        :param log:   a :mod: `logging` object
        """
        # pulls the logger from the BDF object
        self.log = get_logger(log, "debug")
        self.debug = False
        #self.debug = True
        
        #: stores a single copy of 'BEGIN BULK' or 'BEGIN SUPER'
        self.begin_bulk = ['BEGIN', 'BULK']
        
        # allows BEGIN BULK to be turned off
        self.write_begin_bulk = True
        self._begin_count = 0

        
        self.lines = lines
        self.subcases = {0: Subcase(id=0)}
        self._read(self.lines)

    def has_parameter(self, isubcase, param_name):
        """
        Checks to see if a parameter (e.g. STRESS) is defined in a certain
        subcase ID.

        :param self:       the CaseControlDeck object
        :param isubcase:   the subcase ID to check
        :param param_name: the parameter name to look for
        """
        if self.has_subcase(isubcase):
            return self.subcases[isubcase].has_parameter(param_name.upper())

    def get_subcase_parameter(self, isubcase, param_name):
        """
        Get the [value, options] of a subcase's parameter.  For example, for
        STRESS(PLOT,POST)=ALL, param_name=STRESS, value=ALL, options=['PLOT',
        'POST']

        :param self:       the CaseControlDeck object
        :param isubcase:   the subcase ID to check
        :param param_name: the parameter name to get the [value, options] for
        """
        if self.has_subcase(isubcase):
            return self.subcases[isubcase].get_parameter(param_name.upper())
        raise RuntimeError('isubcase=%s does not exist...' % isubcase)

    def has_subcase(self, isubcase):
        """
        Checks to see if a subcase exists.

        :param self:     the CaseControlDeck object
        :param isubcase: the subcase ID
        :type isubcase: int
        :returns val: does_subcase_exist (type = bool)
        """
        if isubcase in self.subcases:
            return True
        return False

    def create_new_subcase(self, isubcase):
        """
        Method create_new_subcase:

        :param isubcase: the subcase ID
        :type isubcase: int
        .. warning:: be careful you dont add data to the global subcase
                     after running this...is this True???
        """
        #print("creating subcase=%s" % isubcase)
        if self.has_subcase(isubcase):
            sys.stderr.write('subcase=%s already exists...skipping\n' %
                             isubcase)
        self.copy_subcase(i_from_subcase=0, i_to_subcase=isubcase,
                          overwrite_subcase=True)
        #self.subcases[isubcase] = Subcase(id=isubcase)

    def delete_subcase(self, isubcase):
        """
        Deletes a subcase.

        :param self:     the CaseControlDeck object
        :param isubcase: the Subcase to delete
        :type isubcase: int
        """
        if not self.has_subcase(isubcase):
            sys.stderr.write('subcase %s doesnt exist...skipping\n' %
                             isubcase)
        del self.subcases[isubcase]

    def copy_subcase(self, i_from_subcase, i_to_subcase, overwrite_subcase=True):
        """
        Overwrites the parameters from one subcase to another.

        :param self:              the CaseControlDeck object
        :param i_from_subcase:    the Subcase to pull the data from
        :param i_to_subcase:      the Subcase to map the data to
        :param overwrite_subcase: NULLs i_to_subcase before copying
                                  i_from_subcase
        """
        #print("copying subcase from=%s to=%s overwrite=%s" % (i_from_subcase, i_to_subcase, overwrite_subcase))
        if not self.has_subcase(i_from_subcase):
            msg = 'iFromSubcase=|%s| does not exist' % i_from_subcase
            raise RuntimeError(msg)
        if overwrite_subcase:
            subcase_from = self.subcases[i_from_subcase]
            subcase_to = copy.deepcopy(subcase_from)
            subcase_to.id = i_to_subcase
            #for key, param in sorted(subcase_from.params.iteritems()):
                #print("going to copy key=%s param=%s" % (key, param))
            self.subcases[i_to_subcase] = subcase_to
        else:
            if not self.has_subcase(i_to_subcase):
                msg = 'i_to_subcase=|%s| does not exist' % i_to_subcase
                raise RuntimeError(msg)
            subcase_to = self.subcases[i_to_subcase]

            for key, param in sorted(subcase_to.iteritems()):
                #print('copying key=%s param=%s' % (key, param))
                if key == 'BEGIN':
                    pass
                subcase_to[key] = copy.deepcopy(param)

    def get_subcase_list(self):
        return sorted(self.subcases.keys())

    def get_local_subcase_list(self):
        key_list = [key for key in self.subcases if key != 0]  # skip the global
        return sorted(key_list)

    def update_solution(self, isubcase, sol):
        """sol = STATICS, FLUTTER, MODAL, etc."""
        self.add_parameter_to_local_subcase(isubcase, 'ANALYSIS %s' % (sol))

    def add_parameter_to_global_subcase(self, param):
        """
        Takes in a single-lined string and adds it to the global Subcase.

        :param self:  the CaseControlDeck object
        :param param: the variable to add
        .. note:: dont worry about overbounding the line
        
        >>> bdf = BDF()
        >>> bdf.read_bdf(bdf_filename)
        >>> bdf.case_control.add_parameter_to_global_subcase('DISP=ALL')
        >>> print bdf.case_control
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

        :param self:     the CaseControlDeck object
        :param isubcase: the subcase ID to add
        :param param:    the variable to add
        .. note::  dont worry about overbounding the line
        
        >>> bdf = BDF()
        >>> bdf.read_bdf(bdf_filename)
        >>> bdf.case_control.add_parameter_to_local_subcase(1, 'DISP=ALL')
        >>> print bdf.case_control
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

        :param self:  the CaseControlDeck object
        :param param: the variable to add
        """
        if '\n' in param or '\r' in param or '\t' in param:
            msg = 'doesnt support embedded endline/tab characters\n'
            msg += '  param = |%r|' % (param)
            raise SyntaxError(msg)
        #self.read([param])
        lines = self._clean_lines([param])
        (j, key, value, options, param_type) = self._parse_entry(lines)
        return (j, key, value, options, param_type)

    def _clean_lines(self, lines):
        """
        Removes comment characters defined by a *$*.

        :param self:  the CaseControlDeck object
        :param lines: the lines to clean.
        """
        lines2 = []
        for line in lines:
            line = line.strip(' \n\r').split('$')[0].rstrip()
            if line:
                lines2.append(line)
        return lines2

    def _read(self, lines):
        """
        Reads the case control deck

        .. note::    supports comment lines
        .. warning:: doesnt check for 72 character width lines, but will
                     follow that when it's written out
        """
        isubcase = 0
        lines = self._clean_lines(lines)
        i = 0
        while i < len(lines):
            line = lines[i]

            lines2 = [line]
            while ',' in lines[i][-1]:
                #print "lines[%s] = %s" %(i,lines[i])
                i += 1
                lines2.append(lines[i])
            (j, key, value, options, paramType) = self._parse_entry(lines2)
            i += 1
            if key == 'BEGIN':
                if 'BULK' not in line and 'SUPER' not in line:
                    raise NotImplementedError('line=%r' % line)
                if self._begin_count == 1:
                    raise NotImplementedError('multiple BEGIN lines are defined...')
                self.begin_bulk = [key, value]
                self._begin_count += 1
                continue
            
            #print("")
            #print("key=|%s| value=|%s| options=|%s| paramType=%s" %(key, value,
            #                                                        options,
            #                                                          paramType))
            isubcase = self._add_parameter_to_subcase(key, value, options,
                                                      paramType, isubcase)

        #print str(self)
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

        :param self:  the CaseControlDeck object
        :param lines: list of lines
        :return: paramName see brief
        :return: value     see brief
        :return: options   see brief
        :return: paramType see brief
        """
        i = 0
        options = []
        value = None
        key = None
        param_type = None

        line = lines[i]
        #print line
        #print "*****lines = ", lines

        equals_count = 0
        for letter in line:
            if letter == '=':
                equals_count += 1
        line_upper = line.upper()

        #print("line_upper = %r" % line)
        if line_upper.startswith('SUBCASE'):
            #print "line = |%r|" %(line)
            line2 = line.replace('=', '')
            sline = line2.split()
            if len(sline) != 2:
                msg = "trying to parse |%s|..." % line
                raise RuntimeError(msg)
            (key, param_type) = sline
            #print("key=|%s| isubcase=|%s|" % (key,isubcase))
            value = int(param_type)
            #self.isubcase = int(isubcase)
            param_type = 'SUBCASE-type'
        elif (line_upper.startswith('LABEL') or
              line_upper.startswith('SUBT') or  # SUBTITLE
              line_upper.startswith('TITL')):   # TITLE
            try:
                eIndex = line.index('=')
            except:
                msg = "cannot find an = sign in LABEL/SUBTITLE/TITLE line\n"
                msg += "line = |%s|" % (line_upper.strip())
                raise RuntimeError(msg)

            key = line[0:eIndex].strip()
            value = line[eIndex + 1:].strip()
            options = []
            param_type = 'STRING-type'
        elif equals_count == 1:  # STRESS
            if '=' in line:
                (key, value) = line.strip().split('=')
            else:
                msg = 'expected item of form "name = value"   line=|%r|' % (
                    line.strip())
                raise RuntimeError(msg)

            key = key.strip()
            value = value.strip()
            if self.debug:
                self.log.debug("key=|%s| value=|%s|" % (key, value))
            param_type = 'STRESS-type'

            if '(' in key:  # comma may be in line - STRESS-type
                #paramType = 'STRESS-type'
                sline = key.strip(')').split('(')
                key = sline[0]
                options = sline[1].split(',')

                # handle TEMPERATURE(INITIAL) and TEMPERATURE(LOAD) cards
                if key == 'TEMPERATURE' or key == 'TEMP':
                    key = 'TEMPERATURE(%s)' % (options[0])
                    options = []
                #print "key=|%s| options=%s" %(key,options)

            elif ' ' in key and ',' in value:  # SET-type
                (key, ID) = key.split()
                key = key + ' ' + ID

                if self.debug:
                    self.log.debug('SET-type key=%s ID=%s' % (key, ID))
                fivalues = value.rstrip(' ,').split(',')  # float/int values

                #: .. todo:: should be more efficient multiline reader...
                # read more lines....
                if line[-1].strip() == ',':
                    i += 1
                    #print "rawSETLine = |%r|" %(lines[i])
                    while 1:
                        if ',' == lines[i].strip()[-1]:
                            fivalues += lines[i][:-1].split(',')
                        else:  # last case
                            fivalues += lines[i].split(',')
                            #print "fivalues last = i=%s |%r|" %(i,lines[i])
                            i += 1
                            break
                        i += 1
                #print "len(fivalues) = ",len(fivalues)
                value = fivalues

                options = ID  # needed a place to put it...
                param_type = 'SET-type'
            elif ',' in value:  # STRESS-type; special TITLE = stuffA,stuffB
                #print 'A ??? line = ',line
                #raise RuntimeError(line)
                pass
            else:  # STRESS-type; TITLE = stuff
                #print 'B ??? line = ',line
                pass
            
            key = update_param_name(key.strip())
            verify_card(key, value, options, line)

        elif line_upper.startswith('BEGIN'):  # begin bulk
            try:
                (key, value) = line_upper.split(' ')
            except:
                msg = 'excepted "BEGIN BULK" found=|%r|' % line
                raise RuntimeError(msg)
            param_type = 'BEGIN_BULK-type'
        elif 'PARAM' in line_upper:  # param
            sline = line.split(',')
            if len(sline) != 3:
                raise SyntaxError("trying to parse |%s|..." % line)
            (key, value, options) = sline
            param_type = 'CSV-type'
        elif ' ' not in line:
            key = line.strip()
            value = line.strip()
            options = None
            param_type = 'KEY-type'
        else:
            msg = 'generic catch all...line=|%r|' % (line)
            key = ''
            value = line
            options = None
            param_type = 'KEY-type'
        i += 1
        return (i, key, value, options, param_type)

    def finish_subcases(self):
        """
        Removes any unwanted data in the subcase...specifically the SUBCASE
        data member.  Otherwise it will print out after a key like stress.

        :param self:  the CaseControlDeck object
        """
        for subcase in self.subcases.itervalues():
            subcase.finish_subcase()

    def convert_to_sol_200(self, model):
        """
        Takes a case control deck and changes it from a SOL xxx to a SOL 200

        :param self: the CaseControlDeck object
        .. todo:: not done...
        """
        analysis = model.rsolmap_toStr[model.sol]
        model.sol = 200

        subcase.add_parameter_to_global_subcase('ANALYSIS', analysis)
        #subcase.add_parameter_to_global_subcase('DESSUB', dessub)

    def _add_parameter_to_subcase(self, key, value, options, param_type,
                                  isubcase):
        """
        Internal method

        self:  the CaseControlDeck object
        """
        if self.debug:
            a = 'key=|%s|' % (key)
            b = 'value=|%s|' % (value)
            c = 'options=|%s|' % (options)
            d = 'param_type=|%s|' % (param_type)
            msg = "_adding isubcase=%s %-12s %-12s %-12s %-12s" % (isubcase, a,
                                                                   b, c, d)
            self.log.debug(msg)

        if key == 'SUBCASE':
            assert value not in self.subcases
            assert isinstance(value, int)
            isubcase = value
            #print "value=", value
            self.copy_subcase(i_from_subcase=0, i_to_subcase=isubcase,
                              overwrite_subcase=True)
            if self.debug:
                msg = "copied subcase iFromSubcase=%s to iToSubcase=%s" % (
                    0, isubcase)
                self.log.debug(msg)
        elif isubcase not in self.subcases:  # initialize new subcase
            #self.isubcase += 1 # is handled in the read code
            msg = 'isubcase=%s is not a valid subcase...' % (isubcase)
            raise RuntimeError(msg)

        subcase = self.subcases[isubcase]
        subcase._add_data(key, value, options, param_type)

        #print "\n%s\n" %(self.subcases[isubcase])
        return isubcase

    def cross_reference(self, model):
        for (isubcase, subcase) in sorted(self.subcases.iteritems()):
            subcase.cross_reference(model)

    def get_op2_data(self):
        """
        Gets the relevant op2 parameters required for a given subcase

        :param self:  the CaseControlDeck object
        .. todo:: not done...
        """
        cases = {}
        for (isubcase, subcase) in sorted(self.subcases.iteritems()):
            if isubcase:
                cases[isubcase] = subcase.getOp2Data(self.sol,
                                                     subcase.solmap_toValue)
        return cases

    def __repr__(self):
        """
        :param self:  the CaseControlDeck object
        """
        msg = ''
        subcase0 = self.subcases[0]
        for subcase in self.subcases.itervalues():
            #print('writing subcase...')
            msg += subcase.write_subcase(subcase0)
        #if len(self.subcases) == 1:
            #msg += 'BEGIN BULK\n'
        if self.write_begin_bulk:
            msg += ' '.join(self.begin_bulk) + '\n'
        return msg

    #def parseParam(self,param):
    #    """
    #    .. warning:: doesnt support comment characters
    #    """
    #    param = param.substr('\n','').substr('\r','') # remove line endings
    #    parse(param)
    #    #param2 = [''.join(param)]
    #    #print 'param2 = ',param2


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

    # these may only be integers
    #print("key =", key)
    if key in ['BCONTACT', 'CURVELINESYMBOL']:
        try:
            value2 = int(value)
        except:
            print('line=%r is invalid; value=%r' % (line, value))
            raise

    # these may only be integers greater than 0
    elif key in ['SPC', 'MPC', 'TRIM', 'FMETHOD', 'METHOD', 'LOAD',
               'SUPORT', 'SUPORT1', 'TEMPERATURE(INITIAL)', 'TEMPERATURE(LOAD)',
               'DLOAD', 'MFLUID', 'CLOAD', 'NLPARM', 'CMETHOD',
               'FREQUENCY', 'TSTEP', 'TSTEPNL', 'SDAMPING', 'DESOBJ',
               'TEMPERATURE(INIT)', 'RANDOM', 'DESSUB', 'ADAPT', 'MAXLINES',
               'TFL','DESGLB', 'SMETHOD', 'DYNRED', 'GUST', 'TEMPERATURE(MATE)',
               'OTIME', 'NONLINEAR', 'AUXM', 'IC', 'BC', 'OUTRCV', 'DIVERG',
               'DATAREC', 'TEMPERATURE(BOTH)', 'DEFORM', 'MODES', 'CASE',
               'SEDR', 'SELG', 'SEFINAL', 'SEKR', 'TEMPERATURE(ESTIMATE)',
               'GPSDCON', 'AUXMODEL',
               'MODTRAK', 'OFREQ', 'DRSPAN', 'OMODES', 'ADACT', 'SERESP', 'STATSUB',
               'CURVESYM', 'ELSDCON', 'CSSCHD', 'NSM', 'TSTRU', 'RANDVAR', ''
               'RGYRO', 'SELR', 'TEMPERATURE(ESTI)', 'RCROSS', 'SERE', 'SEMR',
               '', '', '', '', '', '', '', '', '', '',
               '']:
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
                 'GPSTRESS', 'GPSTRAIN', 'FLUX','AEROF', 'THERMAL', 'STRFIELD',
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
    elif key in ['CSCALE', 'SUBSEQ','SYMSEQ', 'DEFORMATION SCALE', '', '']:
        # floats
        pass
    elif 'SET' in key:
        pass

    # weird cards
    elif key in ['SUBTITLE', 'TITLE',
        'A2GG',  'M2GG', 'K2GG',
        'K2PP', 'M2PP', 
        'K42GG',

        'XMIN', 'XMAX', 'XTITLE','XPAPE', 'XPAPER', 'XAXIS', 'XGRID', 'XGRID LINES', 'XLOG',
        'YMIN', 'YMAX', 'YTITLE','YPAPE', 'YPAPER', 'YAXIS', 'YGRID', 'YGRID LINES', 'YLOG',
        'XTMIN','XTMAX', 'XTGRID', 'XTTITLE', 'XTAXIS', 'XTGRID LINES', 'XTLOG', 
        'YTMIN','YTMAX', 'YTGRID', 'YTTITLE', 'YTAXIS', 'YTGRID LINES', 'YTLOG', 
        'XBMIN', 'XBMAX', 'XBGRID', 'XBAXIS', 'XBGRID LINES', 'XBTITLE', 'XBLOG',
        'YBMIN', 'YBMAX', 'YBGRID', 'YBAXIS', 'YBGRID LINES', 'YBTITLE', 'YBLOG',

         'RIGHT TICS','UPPER TICS',
        'TRIGHT TICS',
        'BRIGHT TICS', 

        'PLOTTER', 'XYPLOT',

        'PTITLE',
        'HOUTPUT', 'PLOTID', '', '', '', '', '', 
        'AXISYMMETRIC', 'CURVELINESYMBOL', 'CURVELINESYMB', 'AECONFIG',
        'B2GG', 'B2PP', 'AESYMXZ', 'TEMP', 'DSAPRT', 'MEFFMASS', 
        'MAXMIN', 'RESVEC',  'MODESELECT', 'RIGID', 'TCURVE',
        'SUPER',  'MAXI DEFO', 'P2G',
        'EXTSEOUT', 'FLSTCNT PREFDB', 'AESYMXY',
        'DSYM', '', '', '']:
        pass
    elif key == 'ANALYSIS':
        assert value in ['HEAT', 'ANALYSIS', 'MFREQ', 'STATICS', 'MODES', 'DFREQ',
            'MTRAN', 'BUCK', 'MCEIG', 'DCEIG', 'SAERO', 'NLSTATIC', 'NLSTAT',
            'STATIC', 'MTRANS', 'MODE', 'FLUTTER', 'DIVERG', 'NLTRAN', 'FLUT', '', '', '', '', ''], 'line=%r is invalid; value=%r' % (line, value)
    elif key == 'AUTOSPC':
        assert value in ['YES'], 'line=%r is invalid; value=%r' % (line, value)
    else:
        raise NotImplementedError('key=%r line=%r' % (key, line))

if __name__ == '__main__':  # pragma: no cover
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
    print('%s\n\n' % deck)

    #deck2 = CaseControlDeck(['ACCELERATION(PLOT,PRINT,PHASE) = ALL',
    #                         'DISPLACEMENT(PLOT,PRINT,PHASE) = ALL',
    #                         'BEGIN BULK'])
    #print('\n\n%s' % deck2)