# pylint: disable=R0904,R0902,C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
import copy

from pyNastran.bdf import subcase
from pyNastran.bdf.subcase import Subcase


class CaseControlDeck(object):
    nlines_max = 10000
    def __init__(self, lines, log=None):
        """
        @param self
          the case control deck object
        @param lines
          list of lines that represent the case control deck ending with
          BEGIN BULK
        @param log
          a logger object
        """
        if log is None:
        #if 1:
            from pyNastran.general.logger import dummyLogger
            word = 'debug'
            loggerObj = dummyLogger()
            log = loggerObj.startLog(word)  # or info
        self.debug = False
        #self.debug = True

        self.log = log
        self.lines = lines
        self.subcases = {0: Subcase(id=0)}
        self._read(self.lines)

    def has_parameter(self, isubcase, param_name):
        """
        Checks to see if a parameter (e.g. STRESS) is defined in a certain
        subcase ID.
        @param self
          the CaseControl object
        @param isubcase
          the subcase ID to check
        @param param_name
          the parameter name to look for
        """
        if self.has_subcase(isubcase):
            return self.subcases[isubcase].has_parameter(param_name.upper())

    def get_subcase_parameter(self, isubcase, param_name):
        if self.has_subcase(isubcase):
            return self.subcases[isubcase].get_parameter(param_name.upper())
        raise RuntimeError('isubcase=%s does not exist...' % isubcase)

    def has_subcase(self, isubcase):
        """
        Checks to see if a subcase exists.
        @param self the case control deck object
        @param isubcase the subcase ID
        @retval does_subcase_exist (type = bool)
        """
        if isubcase in self.subcases:
            return True
        return False

    def create_new_subcase(self, isubcase):
        """
        Method create_new_subcase:
        @warning
         be careful you dont add data to the global subcase after running
         this...is this True???
        """
        if self.has_subcase(isubcase):
            sys.stderr.write('subcase=%s already exists...skipping\n' %
                             (isubcase))
        self.copy_subcase(i_from_subcase=0, i_to_subcase=isubcase,
                          overwrite_subcase=True)
        #self.subcases[iSubcase] = Subcase(id=iSubcase)

    def delete_subcase(self, isubcase):
        if not self.has_subcase(isubcase):
            sys.stderr.write('subcase %s doesnt exist...skipping\n' %
                             (isubcase))
        del self.subcases[isubcase]

    def copy_subcase(self, i_from_subcase, i_to_subcase, overwrite_subcase=True):
        """
        Overwrites the parameters from one subcase to another.
        @param self
          the case control deck object
        @param i_from_subcase
          the subcase to pull the data from
        @param i_to_subcase
          the subcase to map the data to
        @param overwrite_subcase
          NULLs i_to_subcase before copying i_from_subcase
        """
        if not self.has_subcase(i_from_subcase):
            msg = 'iFromSubcase=|%s| does not exist' % (i_from_subcase)
            raise RuntimeError(msg)
        subcase_to = self.subcases[i_from_subcase]
        if overwrite_subcase:
            subcase_to = copy.deepcopy(subcase_to)
            subcase_to.id = i_to_subcase
            self.subcases[i_to_subcase] = subcase_to
        else:
            if not self.has_subcase(i_to_subcase):
                msg = 'i_to_subcase=|%s| does not exist' % (i_to_subcase)
                raise RuntimeError(msg)
            subcase_to = self.subcases[i_to_subcase]
            for key, param in subcase_to.iteritems():
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
        takes in a single-lined string
        @note
            dont worry about overbounding the line
        """
        (j, key, value, options, paramType) = self._parse_data_from_user(param)
        subcase_list = self.get_subcase_list()
        for isubcase in subcase_list:
            self._add_parameter_to_subcase(key, value, options, paramType,
                                           isubcase)

    def add_parameter_to_local_subcase(self, isubcase, param):
        (j, key, value, options, param_type) = self._parse_data_from_user(param)
        self._add_parameter_to_subcase(key, value, options, param_type,
                                       isubcase)

    def _parse_data_from_user(self, param):
        if '\n' in param or '\r' in param or '\t' in param:
            msg = 'doesnt support embedded endline/tab characters\n'
            msg += '  param = |%r|' % (param)
            raise SyntaxError(msg)
        #self.read([param])
        lines = self._clean_lines([param])
        (j, key, value, options, param_type) = self._parse_entry(lines)
        return (j, key, value, options, param_type)

    def _clean_lines(self, lines):
        """removes comment characters $"""
        lines2 = []
        for line in lines:
            line = line.strip(' \n\r').split('$')[0].rstrip()
            if line:
                lines2.append(line)
        return lines2

    def _read(self, lines):
        """
        reads the case control deck
        @note supports comment lines
        @warning
            doesnt check for 72 character width lines, but will follow that
            when it's written out
        """
        iSubcase = 0
        lines = self._clean_lines(lines)
        i = 0
        while i < len(lines):
            line = lines[i]
            #print "rawLine = |%s|" %(line)
            #self.log.debug("rawLine = |%r|" %(line))

            lines2 = [line]
            while ',' in lines[i][-1]:
                #print "lines[%s] = %s" %(i,lines[i])
                i += 1
                lines2.append(lines[i])
                if i > self.nlines_max:
                    msg = 'There are too many lines in case control deck.\n'
                    msg += 'Assuming an infinite loop was found.'
                    raise RuntimeError(msg)
            (j, key, value, options, paramType) = self._parse_entry(lines2)
            i += 1
            #print ""
            #print "key=|%s| value=|%s| options=|%s| paramType=%s" %(key,value, options,paramType)
            iSubcase = self._add_parameter_to_subcase(key, value, options,
                                                      paramType, iSubcase)
            #print "--------------"
            if i == self.nlines_max:
                msg = 'too many lines in Case Control Deck < %i...' %(
                    self.nlines_max)
                raise RuntimeError(msg)
        ###
        #print "done with while loop...\n"

        #print str(self)
        self.finish_subcases()
    ###

    def _parse_entry(self, lines):
        """
        @brief
            internal method for parsing a card of the case control deck

            parses a single case control deck card into 4 sections
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

        @param self  the case control deck object
        @param lines list of lines
        @retval paramName see brief
        @retval value     see brief
        @retval options   see brief
        @retval paramType see brief
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

        if line_upper.startswith('SUBCASE'):
            #print "line = |%r|" %(line)
            line2 = line.replace('=', '')
            sline = line2.split()
            if len(sline) != 2:
                msg = "trying to parse |%s|..." % (line)
                raise RuntimeError(msg)
            (key, param_type) = sline
            #print "key=|%s| iSubcase=|%s|" %(key,iSubcase)
            value = int(param_type)
            #self.iSubcase = int(iSubcase)
            param_type = 'SUBCASE-type'
        elif (line_upper.startswith('LABEL') or
              line_upper.startswith('SUBTITLE') or
              line_upper.startswith('TITLE')):
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

                ## @todo should be more efficient multiline reader...
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
                    ###
                ###
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
            ###
        ### = in line
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
        ###
        i += 1
        #print "done with ",key
        return (i, key, value, options, param_type)

    def finish_subcases(self):
        """
        removes any unwanted data in the subcase...specifically the SUBCASE
        data member.  Otherwise it will print out after a key like stress.
        """
        for subcase in self.subcases.itervalues():
            subcase.finish_subcase()

    def convert_to_sol_200(self, model):
        """
        Takes a case control deck and changes it from a
        @todo not done...
        """
        analysis = model.rsolmap_toStr[model.sol]
        model.sol = 200

        subcase.add_parameter_to_global_subcase('ANALYSIS', analysis)
        #subcase.add_parameter_to_global_subcase('DESSUB', dessub)

    def _add_parameter_to_subcase(self, key, value, options, param_type,
                                  isubcase):
        """internal method"""
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
        for (iSubcase, subcase) in sorted(self.subcases.iteritems()):
            subcase.cross_reference(model)

    def get_op2_data(self):
        """
        returns the relevant op2 parameters required for a given subcase
        """
        cases = {}
        for (iSubcase, subcase) in sorted(self.subcases.iteritems()):
            if iSubcase:
                cases[iSubcase] = subcase.getOp2Data(self.sol,
                                                     subcase.solmap_toValue)
        return cases

    def __repr__(self):
        msg = ''
        subcase0 = self.subcases[0]
        for subcase in self.subcases.itervalues():
            msg += subcase.write_subcase(subcase0)
        if len(self.subcases) == 1:
            msg += 'BEGIN BULK\n'
        return msg

    #def parseParam(self,param):
    #    """
    #    @warning doesnt support comment characters
    #    """
    #    param = param.substr('\n','').substr('\r','') # remove line endings
    #    parse(param)
    #    #param2 = [''.join(param)]
    #    #print 'param2 = ',param2


if __name__ == '__main__':
    test1()
    lines = [
        'SUBCASE 1',
        '    ACCELERATION(PLOT,PRINT,PHASE) = ALL',
        '    DISPLACEMENT(PLOT,PRINT,PHASE) = ALL',
        '    DLOAD = 32',
        '    M2GG = 111',
        '    SET 88  = 5, 6, 7, 8, 9, 10 THRU 55 EXCEPT 15, 16, 77, 78, 79, 100 THRU 300',
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
    deck.create_new_subcase(2)
    deck.add_parameter_to_local_subcase(1, 'SET 2 = 11,12,13,14,15,16,17,18,'
       '19,20,21,22,23,24,25,26,'
       '1000000000000000000000000000000000000000000000000000000,33')
    print(deck + '\n\n')

    deck2 = CaseControlDeck(['ACCELERATION(PLOT,PRINT,PHASE) = ALL',
                             'DISPLACEMENT(PLOT,PRINT,PHASE) = ALL',
                             'BEGIN BULK'])
    print('\n\n' + deck2)
