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
from __future__ import print_function
#from numpy import array
from numpy import angle
from pyNastran.op2.op2Codes import Op2Codes

class baseScalarObject(Op2Codes):
    def __init__(self):
        pass

    def name(self):
        return self.__class__.__name__

    def writeF06(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        msg = ['writeF06 is not implemented in %s\n' %(self.__class__.__name__)]
        return (''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        msg = 'writeF06Transient is not implemented in %s\n' %(self.__class__.__name__)
        return (''.join(msg),pageNum)

    def writeFloats10E(self,vals):
        vals2 = []
        isAllZeros = True
        for v in vals:
            v2 = '%10.3E' %(v)
            if v2==' 0.000E+00' or v2=='-0.000E+00':
                v2 = ' 0.0      '
            else:
                isAllZeros = False
            vals2.append(v2)
        return (vals2,isAllZeros)

    def writeFloats12E(self,vals):
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

    def writeFloats13E(self,vals):
        vals2 = []
        isAllZeros = True
        for v in vals:
            v2 = '%13.6E' %(v)
            if v2==' 0.000000E+00' or v2=='-0.000000E+00':
                v2 = ' 0.0         '
            else:
                isAllZeros = False
            vals2.append(v2)
        return (vals2,isAllZeros)

    def writeImagFloats13E(self,vals,isMagPhase):
        vals2 = []
        isAllZeros = True
        
        if isMagPhase:
            for v in vals:
                v2 = '%13.6E' %(abs(v))
                if v2==' 0.000000E+00' or v2=='-0.000000E+00':
                    v2 = ' 0.0         '
                else:
                    isAllZeros = False
                vals2.append(v2)

            for v in vals:
                v3 = '%13.6E' %(angle(v,deg=True))
                if v3==' 0.000000E+00' or v3=='-0.000000E+00':
                    v3 = ' 0.0         '
                else:
                    isAllZeros = False
                vals2.append(v3)
        else:
            for v in vals:
                v2 = '%13.6E' %(v.real)
                if v2==' 0.000000E+00' or v2=='-0.000000E+00':
                    v2 = ' 0.0         '
                else:
                    isAllZeros = False
                vals2.append(v2)

            for v in vals:
                v3 = '%13.6E' %(v.imag)
                if v3==' 0.000000E+00' or v3=='-0.000000E+00':
                    v3 = ' 0.0         '
                else:
                    isAllZeros = False
                vals2.append(v3)
        return (vals2,isAllZeros)

    def writeFloats8p4F(self,vals):
        vals2 = []
        isAllZeros = True
        for v in vals:
            v2 = '%8.4f' %(v)
            if v2=='  0.0000' or v2==' -0.0000':
                v2 = '  0.0   '
            else:
                isAllZeros = False
            vals2.append(v2)
        return (vals2,isAllZeros)


class scalarObject(baseScalarObject):
    def __init__(self, dataCode, iSubcase):
        assert 'nonlinearFactor' in dataCode, dataCode
        baseScalarObject.__init__(self)
        self.iSubcase = iSubcase
        self.isTransient = False
        self.dt = None
        self.dataCode = dataCode
        self.applyDataCode()
        self.log.debug(self.codeInformation())

    def isImaginary(self):
        return bool(self.sortBits[1])

    def writeMatlabArgs(self,name,iSubcase,f):
        for key,value,in sorted(self.dataCode.iteritems()):
            if key is not 'log':
                if isinstance(value,str):
                    value = "'%s'" %(value)
                    msg = 'fem.%s(%i).%s = %s;\n' %(name,iSubcase,key,value)
                elif isinstance(value,list) and isinstance(value[0],str):
                    msgTemp = "','".join(value)
                    msg = "fem.%s(%i).%s = {'%s'};\n" %(name,iSubcase,key,msgTemp)

                elif value is None:
                    value = "'%s'" %(value)
                else:
                    msg = 'fem.%s(%i).%s = %s;\n' %(name,iSubcase,key,value)
                f.write(msg)
            ###
        ###

    def applyDataCode(self):
        self.log = self.dataCode['log']
        for key,value in sorted(self.dataCode.iteritems()):
            if key is not 'log':
                self.__setattr__(key,value)
                #self.log.debug("  key=%s value=%s" %(key,value))
                #print "  key=%s value=%s" %(key,value)
        #self.log.debug("")

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
                    print("listA = ",listA)
                    raise
                listA.append(value)
                assert len(listA)==n+1
            ###
        ###

    def setDataMembers(self):
        if 'dataNames' not in self.dataCode:
            msg = 'No "transient" variable was set for %s\n' %(self.tableName)
            raise NotImplementedError(msg+self.codeInformation())

        for name in self.dataCode['dataNames']:
            #print "name = ",name
            self.appendDataMember(name+'s',name)
        ###

    def updateDataCode(self,dataCode):
        self.dataCode = dataCode
        self.applyDataCode()
        self.setDataMembers()
        
    def printDataMembers(self):
        """
        Prints out the "unique" vals of the case.
        Uses a provided list of dataCode['dataNames'] to set the values for
        each subcase.  Then populates a list of self.name+'s' (by using 
        setattr) with the current value.  For example, if the variable name is
        'mode', we make self.modes.  Then to extract the values, we build a
        list of of the variables that were set like this and then loop over
        then to print their values.
        
        This way there is no dependency on one result type having ['mode'] and
        another result type having ['mode','eigr','eigi'].
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
        for i in xrange(nModes):
            for vals in keyVals:
                msg += '%-10g ' %(vals[i])
            msg += '\n'
        ###
        return msg+'\n'

    def recastGridType(self,gridType):
        """converts a gridType integer to a string"""
        if gridType==1:
            Type = 'G'  # GRID
        elif gridType==2:
            Type = 'S'  # SPOINT
        elif gridType==7:
            Type = 'L'  # RIGID POINT (e.g. RBE3)
        elif gridType==0:
            Type = 'H'      # SECTOR/HARMONIC/RING POINT
        else:
            raise RuntimeError('gridType=%s' %(gridType))
        ###
        return Type
        
    def updateDt(self,dataCode,dt):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        self.dataCode = dataCode
        self.applyDataCode()
        raise RuntimeError('updateDt not implemented in the %s class' %(self.__class__.__name__))
        #assert dt>=0.
        #print "updating dt...dt=%s" %(dt)
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
        ###

