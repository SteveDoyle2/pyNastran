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
import copy
import sys
from .subcase import Subcase
from pyNastran.bdf.errors import *

class CaseControlDeck(object):
    def __init__(self,lines,log=None):
        """
        @param self  the object pointer
        @param lines list of lines that represent the case control deck ending with BEGIN BULK
        @param log   a logger object
        """
        if log is None:
        #if 1:
            from pyNastran.general.logger import dummyLogger
            loggerObj = dummyLogger()
            log = loggerObj.startLog('debug') # or info
        self.debug = False
        #self.debug = True

        self.log = log
        self.lines = lines
        self.subcases = {0:Subcase(id=0)}
        #self.iSubcase = 0
        self._read(self.lines)

    def hasParameter(self,iSubcase,paramName):
        #paramName = self.updateParamName(paramName)
        if self.hasSubcase(iSubcase):
            return self.subcases[iSubcase].hasParameter(paramName.upper())

    def getSubcaseParameter(self,iSubcase,paramName):
        #paramName = self.updateParamName(paramName)
        #print "iSubcase = ",iSubcase,'\n'
        if self.hasSubcase(iSubcase):
            print(str(self.subcases[iSubcase]))
            return self.subcases[iSubcase].getParameter(paramName.upper())
        ###
        raise RuntimeError('iSubcase=%s does not exist...' %(iSubcase))

    def hasSubcase(self,iSubcase):
        if iSubcase in self.subcases:
            return True
        return False

    def createNewSubcase(self,iSubcase):
        """
        @warning be careful you dont add data to the global subcase after running this...is this True???
        """
        if not hasSubcase(iSubcase):
            sys.stderr.write('subcase=%s already exists...skipping')
        self.copySubcase(iFromSubcase=0,iToSubcase=iSubcase,overwriteSubcase=True)
        #self.subcases[iSubcase] = Subcase(id=iSubcase)

    def deleteSubcase(self,iSubcase):
        if not hasSubcase(iSubcase):
            sys.stderr.write('subcase doesnt exist...skipping')
        del self.subcases[iSubcase]

    def copySubcase(self,iFromSubcase,iToSubcase,overwriteSubcase=True):
        """
        overwrites the parameters from one subcase to another
        @param self             the object pointer
        @param iFromSubcase     the subcase to pull the data from
        @param iToSubcase       the subcase to map the data to
        @param overwriteSubcase NULLs iToSubcase before copying iFromSubcase
        """
        if not self.hasSubcase(iFromSubcase):
            raise RuntimeError('iFromSubcase=|%s| does not exist' %(iFromSubcase))
        subcaseFrom = self.subcases[iFromSubcase]
        if overwriteSubcase:
            #print "inside overwrite..."
            self.subcases[iToSubcase] = copy.deepcopy(self.subcases[iFromSubcase])
            self.subcases[iToSubcase].id = iToSubcase
        else:
            if not hasSubcase(iToSubcase):
                raise RuntimeError('iToSubcase=|%s| does not exist' %(iToSubcase))
            for key,param in subcaseFrom.items():
                subcaseTo[key] = copy.deepcopy(param)
            ###
        ###

    def getSubcaseList(self):
        return sorted(self.subcases.keys())

    def getLocalSubcaseList(self):
        keyList = [key for key in self.subcases if key != 0] # dont get the global
        return sorted(keyList)

    def updateSolution(self,iSubcase,sol):
        """sol = STATICS, FLUTTER, MODAL, etc."""
        self.addParameterToLocalSubcase(self,iSubcase,'ANALYSIS %s')

    def addParameterToGlobalSubcase(self,param):
        """
        takes in a single-lined string
        @note
            dont worry about overbounding the line
        """
        (j,key,value,options,paramType) = self._parseDataFromUser(param)
        subcaseList = self.getSubcaseList()
        for iSubcase in subcaseList:
            self._addParameterToSubcase(key,value,options,paramType,iSubcase)

    def addParameterToLocalSubcase(self,iSubcase,param):
        (j,key,value,options,paramType) = self._parseDataFromUser(param)
        self._addParameterToSubcase(key,value,options,paramType,iSubcase)

    def _parseDataFromUser(self,param):
        if '\n' in param or 'r' in param or '\t' in param:
            msg = 'doesnt support embedded endline/tab characters\n'
            msg += '  param = |%r|' %(param)
            raise SyntaxError(msg)
        #self.read([param])
        lines = self.cleanLines([param])
        (j,key,value,options,paramType) = self._parseEntry(lines)
        return (j,key,value,options,paramType)

    def cleanLines(self,lines):
        """removes comment characters $"""
        lines2 = []
        for line in lines:
            line = line.strip(' \n\r').split('$')[0].rstrip()
            if line:
                lines2.append(line)
            ###
        ###
        #for line in lines2:
        #    print "L2 = ",line
        return lines2

    def _read(self,lines):
        """
        reads the case control deck
        @note supports comment lines
        @warning
            doesnt check for 72 character width lines, but will follow that
            when it's written out
        """
        iSubcase = 0
        lines = self.cleanLines(lines)
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
                if i>100:
                    sys.exit('huhh...')
            (j,key,value,options,paramType) = self._parseEntry(lines2)
            #print "i=%s j=%s" %(i,j)
            i+=1
            #print ""
            #print "key=|%s| value=|%s| options=|%s| paramType=%s" %(key,value,options,paramType)
            iSubcase = self._addParameterToSubcase(key,value,options,paramType,iSubcase)
            #print "--------------"
            if i==600:
                raise Exception('too many lines in Case Control Deck <600...')
        ###
        #print "done with while loop...\n"
        
        #print str(self)
        #sys.exit('stopping...')
        self.finishSubcases()
    ###

    def _parseEntry(self,lines):
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

        @param self  the object pointer
        @param lines list of lines
        @retval paramName see below...
        @retval value     see below...
        @retval options   see below...
        @retval paramType see below...
        """
        i = 0
        options   = []
        value     = None
        key       = None
        paramType = None

        line = lines[i]
        #print line
        #print "*****lines = ",lines
        
        equalsCount = 0
        for letter in line:
            if letter=='=':
                equalsCount +=1
        lineUpper = line.upper()

        if lineUpper.startswith('SUBCASE'):
            #print "line = |%r|" %(line)
            line2 = line.replace('=','')
            sline = line2.split()
            if len(sline)!=2:
                raise InvalidSubcaseParseError("trying to parse |%s|..." %(line))
            (key,iSubcase) = sline
            #print "key=|%s| iSubcase=|%s|" %(key,iSubcase)
            value = int(iSubcase)
            #self.iSubcase = int(iSubcase)
            paramType = 'SUBCASE-type'
        elif lineUpper.startswith('LABEL') or lineUpper.startswith('SUBTITLE') or lineUpper.startswith('TITLE'):
            eIndex = line.index('=')
            key   = line[0:eIndex].strip()
            value = line[eIndex+1:].strip()
            options = []
            paramType = 'STRING-type'
        elif equalsCount==1: # STRESS
            if '=' in line:
                (key,value) = line.strip().split('=')
            else:
                msg = 'expected item of form "name = value"   line=|%r|' %(line.strip())
                raise RuntimeError(msg)

            key   = key.strip()
            value = value.strip()
            if self.debug:  print("key=|%s| value=|%s|" %(key,value))
            paramType = 'STRESS-type'

            if '(' in key:  # comma may be in line - STRESS-type
                #paramType = 'STRESS-type'
                sline = key.strip(')').split('(')
                key = sline[0]
                options = sline[1].split(',')

                # handle TEMPERATURE(INITIAL) and TEMPERATURE(LOAD) cards
                if key=='TEMPERATURE' or key=='TEMP':
                    key = 'TEMPERATURE(%s)' %(options[0])
                    options = []
                #print "key=|%s| options=%s" %(key,options)

            elif ' ' in key and ',' in value: # SET-type
                (key,ID) = key.split()
                fivalues = value.rstrip(' ,').split(',') # float/int values

                ## @todo should be more efficient multiline reader...
                # read more lines....
                if line[-1].strip()==',':
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
                #raise Exception(line)
                pass
            else:  # STRESS-type; TITLE = stuff
                #print 'B ??? line = ',line
                pass
            ###
        ### = in line
        elif lineUpper.startswith('BEGIN'): # begin bulk
            try:
                (key,value) = lineUpper.split(' ')
            except:
                msg = 'excepted "BEGIN BULK" found=|%r|' %(line)
                raise RuntimeError(msg)
            paramType = 'BEGIN_BULK-type'
        elif 'PARAM' in lineUpper: # param
            sline = line.split(',')
            if len(sline) != 3:
                raise ParamParseError("trying to parse |%s|..." %(line))
            (key,value,options) = sline
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
        return (i,key,value,options,paramType)

    def finishSubcases(self):
        """
        removes any unwanted data in the subcase...specifically the SUBCASE
        data member.  Otherwise it will print out after a key like stress.
        """
        for (iSubcase,subcase) in sorted(self.subcases.items()):
            subcase.finishSubcase()
        ###
    ###

    def _addParameterToSubcase(self,key,value,options,paramType,iSubcase):
        """internal method"""
        if self.debug:
            a = 'key=|%s|'       %(key)
            b = 'value=|%s|'     %(value)
            c = 'options=|%s|'   %(options)
            d = 'paramType=|%s|' %(paramType)
            print("_adding iSubcase=%s %-12s %-12s %-12s %-12s" %(iSubcase,a,b,c,d))

        if key=='SUBCASE':
            assert value not in self.subcases
            assert isinstance(value,int)
            iSubcase = value
            #print "value=",value
            self.copySubcase(iFromSubcase=0,iToSubcase=iSubcase,overwriteSubcase=True)
            if self.debug:
                print("copied subcase iFromSubcase=%s to iToSubcase=%s" %(0,iSubcase))
        elif iSubcase not in self.subcases: # initialize new subcase
            #self.iSubcase += 1 # is handled in the read code
            raise

        subcase = self.subcases[iSubcase]
        subcase._addData(key,value,options,paramType)
        
        #print "\n%s\n" %(self.subcases[iSubcase])
        return iSubcase

    #def __str__(self):
    #    return self.__repr__()

    def crossReference(self,mesh):
        for (iSubcase,subcase) in sorted(self.subcases.items()):
            subcase.crossReference(mesh)

    def getOp2Data(self):
        """
        returns the relevant op2 parameters required for a given subcase
        """
        cases = {}
        for (iSubcase,subcase) in sorted(self.subcases.items()):
            if iSubcase != 0:
                cases[iSubcase] = subcase.getOp2Data(self.sol,self.solmap_toValue)
            ###
        ###
        return cases

    def __repr__(self):
        msg = ''
        subcase0 = self.subcases[0]
        for (iSubcase,subcase) in sorted(self.subcases.items()):
            #if iSubcase==0:
            msg += subcase.writeSubcase(subcase0)
            #msg += str(subcase)
            #else:
            #    msg += subcase.writeSubcase(subcase0)
            #print "\n"
            #break
        if len(self.subcases)==1:
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

if __name__=='__main__':
    lines = ['SPC=2',
             'MPC =3',
             'STRESS= ALL',
             'DISPLACEMENT(PLOT,PUNCH) = 8',]
    deck = CaseControlDeck(lines)
    print("has SPC  True  = ",deck.hasParameter(0,'SPC'))
    print("has sPC  True  = ",deck.hasParameter(0,'sPC'))
    print("has junk False = ",deck.hasParameter(0,'JUNK'))
    
    print("getSubcaseParameter(MPC) 3 = ",deck.getSubcaseParameter(0,'MPC'))
    deck.addParameterToGlobalSubcase('STRAIN = 7')
    deck.addParameterToLocalSubcase(1,'STRAIN = 7')
    print("-----added----")

    out = deck.getSubcaseParameter(0,'STRAIN')
    print("getSubcaseParameter(STRAIN) 7 = ",out)

    deck.addParameterToLocalSubcase(2,'SOL=200')
    print("-----added2----")
    out = deck.getSubcaseParameter(2,'SOL')
    print("getSubcaseParameter(SOL) 200 = ",out)
    
