## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
# pylint: disable=R0904,R0902,C0103

from __future__ import division, print_function
import sys
import copy
from pyNastran.bdf.subcase import Subcase
from pyNastran.bdf.errors import ParamParseError

class CaseControlDeck(object):
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
            log = loggerObj.startLog(word) # or info
        self.debug = False
        #self.debug = True

        self.log = log
        self.lines = lines
        self.subcases = {0:Subcase(id=0)}
        self._read(self.lines)

    def hasParameter(self, iSubcase, paramName):
        """@see has_parameter"""
        return self.has_parameter(iSubcase, paramName)

    def getSubcaseParameter(self, iSubcase, paramName):
        """@see get_subcase_parameter"""
        return self.get_subcase_parameter(iSubcase, paramName)

    def hasSubcase(self, iSubcase):
        """@see has_subcase"""
        return self.has_subcase(iSubcase)

    def createNewSubcase(self, iSubcase):
        """@see create_new_subcase"""
        self.create_new_subcase(iSubcase)

    def deleteSubcase(self, iSubcase):
        """@see delete_subcase"""
        self.delete_subcase(iSubcase)

    def copySubcase(self, iFromSubcase, iToSubcase, overwriteSubcase=True):
        """@see copy_subcase"""
        self.copy_subcase(iFromSubcase, iToSubcase, overwriteSubcase=True)

    def getSubcaseList(self):
        """@see get_subcase_list"""
        return self.get_subcase_list()

    def getLocalSubcaseList(self):
        """@see get_local_subcase_list"""
        self.get_local_subcase_list()

    def updateSolution(self, iSubcase, sol):
        """@see update_solution"""
        self.update_solution(iSubcase, sol)

    def addParameterToGlobalSubcase(self, param):
        """@see add_parameter_to_global_subcase"""
        self.add_parameter_to_global_subcase(param)

    def addParameterToLocalSubcase(self, iSubcase, param):
        """@see add_parameter_to_local_subcase"""
        self.add_parameter_to_local_subcase(iSubcase, param)

#-----------------------
    def has_parameter(self, iSubcase, paramName):
        if self.hasSubcase(iSubcase):
            return self.subcases[iSubcase].hasParameter(paramName.upper())

    def get_subcase_parameter(self, iSubcase, paramName):
        if self.hasSubcase(iSubcase):
            return self.subcases[iSubcase].getParameter(paramName.upper())
        raise RuntimeError('iSubcase=%s does not exist...' %(iSubcase))

    def has_subcase(self, iSubcase):
        """
        Checks to see if a subcase exists.
        @param self the case control deck object
        @param iSubcase the subcase ID
        @retval does_subcase_exist (type = bool)
        """
        if iSubcase in self.subcases:
            return True
        return False

    def create_new_subcase(self, iSubcase):
        """
        @warning
          be careful you dont add data to the global subcase after running
          this...is this True???
        """
        if self.hasSubcase(iSubcase):
            sys.stderr.write('subcase=%s already exists...skipping\n' %(iSubcase))
        self.copy_subcase(iFromSubcase=0, iToSubcase=iSubcase,overwriteSubcase=True)
        #self.subcases[iSubcase] = Subcase(id=iSubcase)

    def delete_subcase(self, iSubcase):
        if not self.hasSubcase(iSubcase):
            sys.stderr.write('subcase %s doesnt exist...skipping\n' %(iSubcase))
        del self.subcases[iSubcase]

    def copy_subcase(self, iFromSubcase, iToSubcase, overwriteSubcase=True):
        """
        overwrites the parameters from one subcase to another
        @param self             the case control deck object
        @param iFromSubcase     the subcase to pull the data from
        @param iToSubcase       the subcase to map the data to
        @param overwriteSubcase NULLs iToSubcase before copying iFromSubcase
        """
        if not self.hasSubcase(iFromSubcase):
            msg = 'iFromSubcase=|%s| does not exist' %(iFromSubcase)
            raise RuntimeError(msg)
        subcaseFrom = self.subcases[iFromSubcase]
        if overwriteSubcase:
            subcaseTo = copy.deepcopy(subcaseFrom)
            subcaseTo.id = iToSubcase
            self.subcases[iToSubcase] = subcaseTo
        else:
            if not self.has_subcase(iToSubcase):
                msg = 'iToSubcase=|%s| does not exist' %(iToSubcase)
                raise RuntimeError(msg)
            subcaseTo = self.subcases[iToSubcase]
            for key,param in subcaseFrom.iteritems():
                subcaseTo[key] = copy.deepcopy(param)
            ###
        ###

    def get_subcase_list(self):
        return sorted(self.subcases.keys())

    def get_local_subcase_list(self):
        keyList = [key for key in self.subcases if key != 0] # skip the global
        return sorted(keyList)

    def update_solution(self, iSubcase, sol):
        """sol = STATICS, FLUTTER, MODAL, etc."""
        self.add_parameter_to_local_subcase(iSubcase, 'ANALYSIS %s' %(sol))

    def add_parameter_to_global_subcase(self, param):
        """
        takes in a single-lined string
        @note
            dont worry about overbounding the line
        """
        (j, key, value, options, paramType) = self._parse_data_from_user(param)
        subcaseList = self.get_subcase_list()
        for iSubcase in subcaseList:
            self._add_parameter_to_subcase(key, value, options, paramType,
                                        iSubcase)

    def add_parameter_to_local_subcase(self, iSubcase, param):
        (j, key, value, options, paramType) = self._parse_data_from_user(param)
        self._addParameterToSubcase(key, value, options, paramType, iSubcase)

    def _parse_data_from_user(self, param):
        if '\n' in param or '\r' in param or '\t' in param:
            msg = 'doesnt support embedded endline/tab characters\n'
            msg += '  param = |%r|' %(param)
            raise SyntaxError(msg)
        #self.read([param])
        lines = self.clean_lines([param])
        (j, key, value, options, paramType) = self._parse_entry(lines)
        return (j, key, value, options, paramType)

    def _clean_lines(self, lines):
        """removes comment characters $"""
        lines2 = []
        for line in lines:
            line = line.strip(' \n\r').split('$')[0].rstrip()
            if line:
                lines2.append(line)
            ###
        ###
        #for line in lines2:
            #print "L2 = ",line
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
                i+=1
                lines2.append(lines[i])
                if i>10000:
                    msg =  'There are too many lines in case control deck.\n'
                    msg += 'Assuming an infinite loop was found.'
                    raise RuntimeError(msg)
            (j, key, value, options, paramType) = self._parse_entry(lines2)
            i+=1
            #print ""
            #print "key=|%s| value=|%s| options=|%s| paramType=%s" %(key,value,options,paramType)
            iSubcase = self._addParameterToSubcase(key, value, options, paramType, iSubcase)
            #print "--------------"
            if i == 10000:
                msg = 'too many lines in Case Control Deck < 10000...'
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
        options   = []
        value     = None
        key       = None
        paramType = None

        line = lines[i]
        #print line
        #print "*****lines = ", lines
        
        equalsCount = 0
        for letter in line:
            if letter == '=':
                equalsCount += 1
        lineUpper = line.upper()

        if lineUpper.startswith('SUBCASE'):
            #print "line = |%r|" %(line)
            line2 = line.replace('=', '')
            sline = line2.split()
            if len(sline) != 2:
                msg = "trying to parse |%s|..." %(line)
                raise RuntimeError(msg)
            (key, iSubcase) = sline
            #print "key=|%s| iSubcase=|%s|" %(key,iSubcase)
            value = int(iSubcase)
            #self.iSubcase = int(iSubcase)
            paramType = 'SUBCASE-type'
        elif (lineUpper.startswith('LABEL') or lineUpper.startswith('SUBTITLE')
              or lineUpper.startswith('TITLE')):
            try:
                eIndex = line.index('=')
            except:
                msg  = "cannot find an = sign in LABEL/SUBTITLE/TITLE line\n"
                msg += "line = |%s|" %(lineUpper.strip())
                raise RuntimeError(msg)

            key   = line[0:eIndex].strip()
            value = line[eIndex+1:].strip()
            options = []
            paramType = 'STRING-type'
        elif equalsCount == 1: # STRESS
            if '=' in line:
                (key, value) = line.strip().split('=')
            else:
                msg = 'expected item of form "name = value"   line=|%r|' %(line.strip())
                raise RuntimeError(msg)

            key = key.strip()
            value = value.strip()
            if self.debug:
                self.log.debug("key=|%s| value=|%s|" %(key, value))
            paramType = 'STRESS-type'

            if '(' in key:  # comma may be in line - STRESS-type
                #paramType = 'STRESS-type'
                sline = key.strip(')').split('(')
                key = sline[0]
                options = sline[1].split(',')

                # handle TEMPERATURE(INITIAL) and TEMPERATURE(LOAD) cards
                if key == 'TEMPERATURE' or key == 'TEMP':
                    key = 'TEMPERATURE(%s)' %(options[0])
                    options = []
                #print "key=|%s| options=%s" %(key,options)

            elif ' ' in key and ',' in value: # SET-type
                (key, ID) = key.split()
                key = key+' '+ID

                if self.debug:
                    self.log.debug('SET-type key=%s ID=%s' %(key, ID))
                fivalues = value.rstrip(' ,').split(',') # float/int values

                ## @todo should be more efficient multiline reader...
                # read more lines....
                if line[-1].strip() == ',':
                    i+=1
                    #print "rawSETLine = |%r|" %(lines[i])
                    while 1:
                        if ','== lines[i].strip()[-1]:
                            fivalues += lines[i][:-1].split(',')
                        else: # last case
                            fivalues += lines[i].split(',')
                            #print "fivalues last = i=%s |%r|" %(i,lines[i])
                            i+=1
                            break
                        i+=1
                    ###
                ###
                #print "len(fivalues) = ",len(fivalues)
                value = fivalues

                options = ID # needed a place to put it...
                paramType = 'SET-type'
            elif ',' in value: # STRESS-type; special TITLE = stuffA,stuffB
                #print 'A ??? line = ',line
                #raise RuntimeError(line)
                pass
            else:  # STRESS-type; TITLE = stuff
                #print 'B ??? line = ',line
                pass
            ###
        ### = in line
        elif lineUpper.startswith('BEGIN'): # begin bulk
            try:
                (key, value) = lineUpper.split(' ')
            except:
                msg = 'excepted "BEGIN BULK" found=|%r|' %(line)
                raise RuntimeError(msg)
            paramType = 'BEGIN_BULK-type'
        elif 'PARAM' in lineUpper: # param
            sline = line.split(',')
            if len(sline) != 3:
                raise ParamParseError("trying to parse |%s|..." %(line))
            (key, value, options) = sline
            ###
            paramType = 'CSV-type'
        elif ' ' not in line:
            key = line.strip()
            value = line.strip()
            options = None
            paramType = 'KEY-type'
        else:
            msg = 'generic catch all...line=|%r|' %(line)
            #print 'C ??? line = ',line
            #raise RuntimeError(msg)
            key = ''
            value = line
            options = None
            paramType = 'KEY-type'
        ###
        i+=1
        #print "done with ",key
        return (i, key, value, options, paramType)

    def finish_subcases(self):
        """
        removes any unwanted data in the subcase...specifically the SUBCASE
        data member.  Otherwise it will print out after a key like stress.
        """
        for (iSubcase, subcase) in sorted(self.subcases.iteritems()):
            subcase.finish_subcase()
        ###
    ###

    def _addParameterToSubcase(self, key, value, options, paramType, iSubcase):
        """internal method"""
        if self.debug:
            a = 'key=|%s|'       %(key)
            b = 'value=|%s|'     %(value)
            c = 'options=|%s|'   %(options)
            d = 'paramType=|%s|' %(paramType)
            msg = "_adding iSubcase=%s %-12s %-12s %-12s %-12s" %(iSubcase, a,
                                                                  b, c, d)
            self.log.debug(msg)

        if key == 'SUBCASE':
            assert value not in self.subcases
            assert isinstance(value, int)
            iSubcase = value
            #print "value=", value
            self.copy_subcase(iFromSubcase=0, iToSubcase=iSubcase,
                             overwriteSubcase=True)
            if self.debug:
                msg = "copied subcase iFromSubcase=%s to iToSubcase=%s" %(0, iSubcase)
                self.log.debug(msg)
        elif iSubcase not in self.subcases: # initialize new subcase
            #self.iSubcase += 1 # is handled in the read code
            msg = 'iSubcase=%s is not a valid subcase...' %(iSubcase)
            raise RuntimeError(msg)

        subcase = self.subcases[iSubcase]
        subcase._add_data(key,value,options,paramType)
        
        #print "\n%s\n" %(self.subcases[iSubcase])
        return iSubcase

    #def __str__(self):
    #    return self.__repr__()

    def crossReference(self, model):
        for (iSubcase, subcase) in sorted(self.subcases.iteritems()):
            subcase.crossReference(model)

    def get_op2_data(self):
        """
        returns the relevant op2 parameters required for a given subcase
        """
        cases = {}
        for (iSubcase, subcase) in sorted(self.subcases.iteritems()):
            if iSubcase != 0:
                cases[iSubcase] = subcase.getOp2Data(self.sol,
                                                     self.solmap_toValue)
        return cases

    def __repr__(self):
        msg = ''
        subcase0 = self.subcases[0]
        for (iSubcase, subcase) in sorted(self.subcases.iteritems()):
            #if iSubcase==0:
            #print("iSubcase = %s" %(iSubcase))
            #print(subcase)
            #print("********")
            msg += subcase.write_subcase(subcase0)
            #msg += str(subcase)
            #else:
            #    msg += subcase.writeSubcase(subcase0)
            #print "\n"
            #break
        if len(self.subcases) == 1:
            msg += 'BEGIN BULK\n'
        #print msg
        return msg
    ###
###

    #def parseParam(self,param):
    #    """
    #    @warning doesnt support comment characters
    #    """
    #    param = param.substr('\n','').substr('\r','') # remove line endings
    #    parse(param)
    #    #param2 = [''.join(param)]
    #    #print 'param2 = ',param2

def test1():
    lines = ['SPC=2',
             'MPC =3',
             'STRESS= ALL',
             'DISPLACEMENT(PLOT,PUNCH) = 8',]
    deck = CaseControlDeck(lines)
    print("has SPC  True  = %s" %(deck.has_parameter(0, 'SPC')))
    print("has sPC  True  = %s" %(deck.has_parameter(0, 'sPC')))
    print("has junk False = %s" %(deck.has_parameter(0, 'JUNK')))
    
    print("getSubcaseParameter(MPC) 3 = ", deck.get_subcase_parameter(0, 'MPC'))
    deck.add_parameter_to_global_subcase('GPFORCE = 7')
    print("")
    print(deck)
    deck.create_new_subcase(1)
    deck.create_new_subcase(2)
    print(deck)
    
    deck.addParameterToLocalSubcase(1, 'STRAIN = 7')
    #print "-----added----"

    out = deck.getSubcaseParameter(0, 'GPFORCE')
    print("getSubcaseParameter(STRAIN) 7 = %s" %(out))

    deck.addParameterToLocalSubcase(1, 'ANALYSIS = SAERO')
    deck.addParameterToLocalSubcase(2, 'ANALYSIS = STATIC')
    print("-----added2----")
    out = deck.getSubcaseParameter(2, 'ANALYSIS')
    print("getSubcaseParameter(ANALYSIS) = %s" %(out))
    

    deck.addParameterToLocalSubcase(1, 'SET 1 = 100')
    deck.addParameterToLocalSubcase(1, 'SET 2 = 200')
    print(deck)
    
if __name__=='__main__':
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
    deck.createNewSubcase(2)
    deck.addParameterToLocalSubcase(1,'SET 2 = 11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,1000000000000000000000000000000000000000000000000000000,33')
    print(deck+'\n\n')

    deck2 = CaseControlDeck(['ACCELERATION(PLOT,PRINT,PHASE) = ALL',
                             'DISPLACEMENT(PLOT,PRINT,PHASE) = ALL',
                             'BEGIN BULK'])
    print('\n\n'+deck2)
