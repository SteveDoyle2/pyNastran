from numpy import array
from pyNastran.op2.op2Errors import *
from pyNastran.op2.op2Codes import Op2Codes

class scalarObject(Op2Codes):
    def __init__(self,dataCode,iSubcase):
        self.iSubcase = iSubcase
        self.isTransient = False
        self.dt = None
        self.dataCode = dataCode
        self.applyDataCode()
        self.log.debug(self.codeInformation())

    def isImaginary(self):
        return bool(self.sortBits[1])

    def name(self):
        return self.__class__.__name__

    def writeF06Floats12E(self,vals):
        vals2 = []
        isAllZeros = True
        for v in vals:
            v2 = '%12.5E' %(v)
            if v2==' 0.00000E+00' or v2=='-0.00000E+00':
                v2 = ' 0.0        '
            else:
                isAllZeros = False
            vals2.append(v2)
        return (vals2,isAllZeros)

    def writeF06Floats13E(self,vals):
        vals2 = []
        isAllZeros = True
        for v in vals:
            v2 = '%13E' %(v)
            if v2==' 0.000000E+00' or v2=='-0.000000E+00':
                v2 = ' 0.0         '
            else:
                isAllZeros = False
            vals2.append(v2)
        return (vals2,isAllZeros)

    def applyDataCode(self):
        self.log = self.dataCode['log']
        for key,value in sorted(self.dataCode.items()):
            if key is not 'log':
                self.__setattr__(key,value)
                self.log.debug("  key=%s value=%s" %(key,value))
                #print "  key=%s value=%s" %(key,value)

    def getUnsteadyValue(self):
        name = self.dataCode['name']
        return self.getVar(name)
        
    def getVar(self,name):
        return getattr(self,name)

    def setVar(self,name,value):
        return self.__setattr__(name,value)

    def startDataMember(self,varName,valueName):
        if hasattr(self,varName):
            return True
        elif hasattr(self,valueName):
            self.setVar(varName,[])
            return True
        return False
        
    def appendDataMember(self,varName,valueName):
        """this appends a data member to a variable that may or may not exist"""
        #print "append..."
        hasList = self.startDataMember(varName,valueName)
        if hasList:
            listA = self.getVar(varName)
            if listA is not None:
                #print "has %s" %(varName)
                value = self.getVar(valueName)
                try:
                    n = len(listA)
                except:
                    print "listA = ",listA
                    raise
                listA.append(value)
                assert len(listA)==n+1
            ###
        ###
    def setDataMembers(self):
        for name in self.dataCode['dataNames']:
            #print "name = ",name
            self.appendDataMember(name+'s',name)
        ###

    def printDataMembers(self):
        """
        Prints out the "unique" vals of the case.
        Uses a provided list of dataCode['dataNames'] to set the values for each
        subcase.  Then populates a list of self.name+'s' (by using setattr)
        with the current value.  For example, if the variable name is 'mode', 
        we make self.modes.  Then to extract the values, we build a list of of the
        variables that were set like this and then loop over then to print their values.
        
        This way there is no dependency on one result type having ['mode'] and another 
        result type having ['mode','eigr','eigi'].
        """
        keyVals = []
        for name in self.dataCode['dataNames']:
            vals = getattr(self,name+'s')
            keyVals.append(vals)
            #print "%ss = %s" %(name,vals)
        
        msg = ''
        for name in self.dataCode['dataNames']:
            msg += '%-10s ' %(name)
        msg += '\n'
        
        nModes = len(keyVals[0])
        for i in range(nModes):
            for vals in keyVals:
                msg += '%-10g ' %(vals[i])
            msg += '\n'
        ###
        return msg+'\n'
            
    def printDataMember(word,selfVarName):
        msg = ''
        if self.getVar(selfVarName):
            msg += '%s = %s' %(word,selfVarName)
        raise Exception('do i need this...msg=%s' %(msg))
        return msg

    def recastGridType(self,gridType):
        """converts a gridType integer to a string"""
        if gridType==1:
            gridType = 'G'  # GRID
        elif gridType==2:
            gridType = 'S'  # SPOINT
        elif gridType==7:
            gridType = 'L'  # RIGID POINT ???
        else:
            raise Exception('gridType=%s' %(gridType))
        ###
        return gridType
        
    def updateDt(self,dataCode,dt):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        self.dataCode = dataCode
        self.applyDataCode()
        raise Exception('updateDt not implemented in the %s class' %(self.__class__.__name__))
        #assert dt>=0.
        #print "updating dt...dt=%s" %(dt)
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
        ###

