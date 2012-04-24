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
import sys
import copy
from struct import unpack

# pyNastran
#from pyNastran.op2.resultObjects.ougv1_Objects import (
#    temperatureObject,
#    nonlinearTemperatureObject,
#    fluxObject,nonlinearFluxObject)
from opg_Objects import appliedLoadsObject
from opg_loadVector import loadVectorObject

class OPG(object):
    """Table of element forces"""
    def readTable_OPG(self):
        table3     = self.readTable_OPG_3
        table4Data = self.readOPG1_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OPG()

    def deleteAttributes_OPG(self):
        params = ['lsdvm','mode','eigr','eign','eigi','modeCycle','freq','time','lftsfq','dLoadID','formatCode','numWide','oCode']
        self.deleteAttributes(params)

    def readTable_OPG_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)
        
        #self.printBlock(data)
        
        
        aCode = self.getBlockIntEntry(data,1)
        #print "aCode = ",aCode
        self.parseApproachCode(data)
        #iSubcase = self.getValues(data,'i',4)

        self.addDataParameter(data,'dLoadID',     'i',8,False)   ## dynamic load set ID/random code
        self.addDataParameter(data,'formatCode',  'i',9,False)   ## format code
        self.addDataParameter(data,'numWide',     'i',10,False)  ## number of words per entry in record; @note is this needed for this table ???
        self.addDataParameter(data,'oCode',       'i',11,False)   ## undefined in DMAP...
        self.addDataParameter(data,'thermal',     'i',23,False)  ## thermal flag; 1 for heat ransfer, 0 otherwise

        #print "dLoadID(8)=%s formatCode(9)=%s numWide(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.formatCode,self.numWide,self.oCode,self.thermal)
        
        ## assuming tCode=1
        if self.analysisCode==1:   # statics
            self.addDataParameter(data,'lsdvmn',  'i',5,False)   ## load set number
        elif self.analysisCode==2: # normal modes/buckling (real eigenvalues)
            self.addDataParameter(data,'mode',     'i',5)   ## mode number
            self.addDataParameter(data,'eign',     'f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'modeCycle','f',7,False)   ## mode or cycle @todo confused on the type - F1???
        #elif self.analysisCode==3: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        #elif self.analysisCode==4: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.analysisCode==5:   # frequency
            self.addDataParameter(data,'freq','f',5)   ## frequency

        elif self.analysisCode==6: # transient
            self.addDataParameter(data,'time','f',5)   ## time step
            #print "TIME(5)=%s" %(self.time)
        elif self.analysisCode==7: # pre-buckling
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            #print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.analysisCode==8: # post-buckling
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            self.addDataParameter(data,'eigr',    'f',6,False)   ## real eigenvalue
            #print "LSDVMN(5)=%s  EIGR(6)=%s" %(self.lsdvmn,self.eigr)
        elif self.analysisCode==9: # complex eigenvalues
            self.addDataParameter(data,'mode','i',5)   ## mode number
            self.addDataParameter(data,'eigr','f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'eigi','f',7,False)   ## imaginary eigenvalue
            #print "mode(5)=%s  eigr(6)=%s  eigi(7)=%s" %(self.mode,self.eigr,self.eigi)
        elif self.analysisCode==10: # nonlinear statics
            self.addDataParameter(data,'lftsfq','f',5)   ## load step
            #print "LFTSFQ(5) = %s" %(self.lftsfq)
        elif self.analysisCode==11: # old geometric nonlinear statics
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            #print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.analysisCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            #print "LSDVMN(5)=%s" %(self.lsdvmn)
        else:
            raise InvalidAnalysisCodeError('invalid analysisCode...analysisCode=%s' %(self.analysisCode))
        ###
        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set, Mode number
        
        #print "*iSubcase=%s"%(self.iSubcase)
        #print "analysisCode=%s tableCode=%s thermal=%s" %(self.analysisCode,self.tableCode,self.thermal)
        #print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

    def readOPG1_Data(self):
        #print "self.analysisCode=%s tableCode(1)=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.thermal)
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        self.atfsCode = [self.analysisCode,self.tableCode,self.formatCode,self.sortCode]

        # grid point force balance
        if   tfsCode==[19,1,0]:
            self.readOPG1_Data_table19_format1_sort0()
        #elif tfsCode==[19,1,1]:
        #    self.readOPG1_Data_format1_sort1()
        #elif tfsCode==[19,2,1]:
        #    self.readOPG1_Data_format2_sort1()
        #elif tfsCode==[19,3,0]:
        #    self.readOPG1_Data_format3_sort0()
        #elif tfsCode==[19,3,1]:
        #    self.readOPG1_Data_format3_sort1()
        
        # load vector
        elif tfsCode==[2,1,0]:
            self.readOPG1_Data_table2_format1_sort0()
        #elif tfsCode==[2,1,1]:
        #    self.readOPG1_Data_format1_sort1()
        #elif tfsCode==[2,2,1]:
        #    self.readOPG1_Data_format2_sort1()
        #elif tfsCode==[2,3,0]:
        #    self.readOPG1_Data_format3_sort0()
        #elif tfsCode==[2,3,1]:
        #    self.readOPG1_Data_format3_sort1()

        # Nonlinear force vector
        #elif tfsCode==[12,1,0]:
        #    self.readOPG1_Data_format1_sort0()

        # OGS1- grid point stresses - surface
        #elif tfsCode==[26,1,0]:
        #    self.readOPG1_Data_format1_sort0()

        # OGS1- grid point stresses - volume direct
        #elif tfsCode==[27,1,0]:
        #    self.readOPG1_Data_format1_sort0()

        # OGS1- grid point stresses - principal
        #elif tfsCode==[28,1,0]:
        #    self.readOPG1_Data_format1_sort0()

        # OGS - Grid point stress discontinuities (plane strain)
        #elif tfsCode==[35,1,0]:
        #    self.readOPG1_Data_format1_sort0()

        # OFMPF2M - does this belong here?
        #elif tfsCode==[51,3,3]:
        #    self.readOPG1_Data_format3_sort3()
        
        # OSMPF2M - does this belong here?
        #elif tfsCode==[52,3,3]:
        #    self.readOPG1_Data_format3_sort3()
        
        # OPMPF2M - does this belong here?
        #elif tfsCode==[53,3,3]:
        #    self.readOPG1_Data_format3_sort3()
        
        # OLMPF2M - does this belong here?
        #elif tfsCode==[54,3,3]:
        #    self.readOPG1_Data_format3_sort3()

        # OGMPF2M - does this belong here?
        #elif tfsCode==[55,3,3]:
        #    self.readOPG1_Data_format3_sort3()
        else:
            #raise Exception('bad tableCode/formatCode/sortCode=%s on %s-OPG table' %(self.atfsCode,self.tableName,))
            #print 'bad analysis/tableCode/formatCode/sortCode=%s on %s-OPG table' %(self.atfsCode,self.tableName)
            self.skipOES_Element()
        ###
        #print self.obj

    #def readOPG1_Data_format1_sort1(self):
        #print 'not supported %s-OPG solution...atfsCode=%s' %(self.tableName,self.atfsCode)
        #self.skipOES_Element()

    #def readOPG1_Data_format2_sort1(self):
        #print 'not supported %s-OPG solution...atfsCode=%s' %(self.tableName,self.atfsCode)
        #self.skipOES_Element()

    #def readOPG1_Data_format3_sort0(self):
        #print 'not supported %s-OPG solution...atfsCode=%s' %(self.tableCode,self.atfsCode)
        #self.skipOES_Element()

    #def readOPG1_Data_format3_sort1(self):
        #print 'not supported %s-OPG solution...atfsCode=%s' %(self.tableCode,self.atfsCode)
        #self.skipOES_Element()

    #def readOPG1_Data_format3_sort3(self):
        #print 'not supported %s-OPG solution...atfsCode=%s' %(self.tableCode,self.atfsCode)
        #self.skipOES_Element()

    def readOPG1_Data_table19_format1_sort0(self):
        if self.thermal==0:
            if self.analysisCode==1: # static
                #print "isAppliedLoads"
                self.createTransientObject(self.appliedLoads,appliedLoadsObject)
                self.readOPGForces(self.data,self.obj)
            elif self.analysisCode==6: # transient
                #print "isAppliedLoads"
                self.createTransientObject(self.appliedLoads,appliedLoadsObject)
                self.readOPGForces(self.data,self.obj)
            elif self.analysisCode==10: # nonlinear static
                #print "isAppliedLoads"
                self.createTransientObject(self.appliedLoads,appliedLoadsObject)
                self.readOPGForces(self.data,self.obj)
            #elif self.analysisCode==11: # old nonlinear static
                #self.createTransientObject(self.appliedLoads,appliedLoadsObject)
            else:
                self.skipOES_Element()
                #print 'not supported %s-OPG solution...atfsCode=%s' %(self.tableName,self.atfsCode)
                #raise Exception('not supported OPG solution...')
            ###
        elif self.thermal==1:
            self.skipOES_Element()
        else:
            raise NotImplementedError('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        
        #readCase = True
        #if self.iSubcase in self.expectedTimes and len(self.expectedTimes[self.iSubcase])>0:
        #    readCase = self.updateDtMap()
        
        #if self.obj and readCase:
        #    self.readScalarsOut(debug=False)
        #else:
        #    self.skipOES_Element()
        ###
        

    def readOPG1_Data_table2_format1_sort0(self):
        if self.thermal==0:
            if self.analysisCode==1: # load vector
                self.createTransientObject(self.loadVectors,loadVectorObject)
            elif self.analysisCode==6: # transient
                self.createTransientObject(self.loadVectors,loadVectorObject)
            elif self.analysisCode==7: # pre-buckling
                self.createTransientObject(self.loadVectors,loadVectorObject)
            elif self.analysisCode==10: # nonlinear static
                self.createTransientObject(self.loadVectors,loadVectorObject)
            elif self.analysisCode==11: # old nonlinear static
                self.createTransientObject(self.loadVectors,loadVectorObject)
            else:
                #self.skipOES_Element()
                pass
                #print 'not supported %s-OPG solution...atfsCode=%s' %(self.tableName,self.atfsCode)
                #print self.codeInformation()
                #raise NotImplementedError('bad approach/table/format/sortCode=%s on %s-OPG table' %(self.atfsCode,self.tableName,))
            ###
        elif self.thermal==1:
            self.skipOES_Element()
        else:
            raise NotImplementedError('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        self.readMappedScalarsOut(debug=False) # handles dtMap

    def readOPGForces(self,data,scalarObject):
        deviceCode = self.deviceCode
        #self.printBlock(data[0:self.numWide*4])
        dn = self.numWide*4
        while len(data)>dn:
            #print "len(data) = ",len(data)
            eData = data[0:dn]
            #self.printBlock(data[:dn])
            gridDevice,eid = unpack('ii',data[0:8])
            nodeID = (gridDevice-deviceCode) // 10
            
            source = ''.join(unpack('cccccccc',data[8:16]))
            #self.printBlock(data[16:dn])
            (dx,dy,dz,rx,ry,rz) = unpack('ffffff',data[16:40])
            #print "source = |%s|" %(source)
            #print "len(data[8:40]"
            
            #print "nodeID=%s eid=%s source=|%s| dx=%-4i dy=%-4i dz=%-4i rx=%-4i ry=%-4i rz=%-4i" %(nodeID,eid,source,dx,dy,dz,rx,ry,rz)
            source2 = source.replace('*','').replace('-','').strip()
            assert source2.isalnum(),'source=|%s| contains invalid characters...' %(source)
            
            scalarObject.add(nodeID,eid,source,dx,dy,dz,rx,ry,rz)
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-self.deviceCode) // 10
            #print "grid=%g dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            #print type(scalarObject)
            data = data[dn:]
            #sys.exit('asd')
        ###
        #print "***********"
        #print self.obj
        #sys.exit('check...')

