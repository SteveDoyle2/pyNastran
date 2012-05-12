import sys
import copy
from struct import unpack

# pyNastran
#from pyNastran.op2.resultObjects.ougv1_Objects import (
#    temperatureObject,
#    nonlinearTemperatureObject,
#    fluxObject,nonlinearFluxObject)
from opg_Objects      import appliedLoadsObject
from opg_loadVector   import loadVectorObject,complexLoadVectorObject
from opnl_forceVector import forceVectorObject,complexForceVectorObject

class OPG(object):
    """Table of element forces"""
    def readTable_OPG(self):
        table3     = self.readTable_OPG_3
        table4Data = self.readOPG_Data
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
        if not self.isSort1():
            raise NotImplementedError('sort2...')
        assert self.isThermal()==False,self.thermal
        
        ## assuming tCode=1
        if self.analysisCode==1:   # statics
            self.addDataParameter(data,'lsdvmn',  'i',5,False)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysisCode==2: # normal modes/buckling (real eigenvalues)
            self.addDataParameter(data,'mode',     'i',5)   ## mode number
            self.addDataParameter(data,'eign',     'f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'modeCycle','f',7,False)   ## mode or cycle @todo confused on the type - F1???
            self.applyDataCodeValue('dataNames',['mode','eign','modeCycle'])
        #elif self.analysisCode==3: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        #elif self.analysisCode==4: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.analysisCode==5:   # frequency
            self.addDataParameter(data,'freq','f',5)   ## frequency
            self.applyDataCodeValue('dataNames',['freq'])
        elif self.analysisCode==6: # transient
            self.addDataParameter(data,'time','f',5)   ## time step
            self.applyDataCodeValue('dataNames',['time'])
        elif self.analysisCode==7: # pre-buckling
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode==8: # post-buckling
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            self.addDataParameter(data,'eigr',    'f',6,False)   ## real eigenvalue
            self.applyDataCodeValue('dataNames',['lsdvmn','eigr'])
        elif self.analysisCode==9: # complex eigenvalues
            self.addDataParameter(data,'mode','i',5)   ## mode number
            self.addDataParameter(data,'eigr','f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'eigi','f',7,False)   ## imaginary eigenvalue
            self.applyDataCodeValue('dataNames',['mode','eigr','eigi'])
        elif self.analysisCode==10: # nonlinear statics
            self.addDataParameter(data,'lftsfq','f',5)   ## load step
            self.applyDataCodeValue('dataNames',['lftsfq'])
        elif self.analysisCode==11: # old geometric nonlinear statics
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
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

    def readOPG_Data(self):
        #print "self.analysisCode=%s tableCode(1)=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.thermal)
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        self.atfsCode = [self.analysisCode,self.tableCode,self.formatCode,self.sortCode]

        
        if   self.tableCode==19:
            self.readOPG_Data_table19() # grid point force balance
        elif self.tableCode==2:  # load vector
            self.readOPG_Data_table2()
        elif self.tableCode==12:  # nonlinear force vector
            self.readOPG_Data_table12()

        # Nonlinear force vector
        #elif tfsCode==[12,1,0]:
        #    self.readOPG_Data_format1_sort0()

        # OGS1- grid point stresses - surface
        #elif tfsCode==[26,1,0]:
        #    self.readOPG_Data_format1_sort0()

        # OGS1- grid point stresses - volume direct
        #elif tfsCode==[27,1,0]:
        #    self.readOPG_Data_format1_sort0()

        # OGS1- grid point stresses - principal
        #elif tfsCode==[28,1,0]:
        #    self.readOPG_Data_format1_sort0()

        # OGS - Grid point stress discontinuities (plane strain)
        #elif tfsCode==[35,1,0]:
        #    self.readOPG_Data_format1_sort0()

        # OFMPF2M - does this belong here?
        #elif tfsCode==[51,3,3]:
        #    self.readOPG_Data_format3_sort3()
        
        # OSMPF2M - does this belong here?
        #elif tfsCode==[52,3,3]:
        #    self.readOPG_Data_format3_sort3()
        
        # OPMPF2M - does this belong here?
        #elif tfsCode==[53,3,3]:
        #    self.readOPG_Data_format3_sort3()
        
        # OLMPF2M - does this belong here?
        #elif tfsCode==[54,3,3]:
        #    self.readOPG_Data_format3_sort3()

        # OGMPF2M - does this belong here?
        #elif tfsCode==[55,3,3]:
        #    self.readOPG_Data_format3_sort3()
        else:
            print self.codeInformation()
            #raise Exception('bad tableCode/formatCode/sortCode=%s on %s-OPG table' %(self.atfsCode,self.tableName,))
            #print 'bad analysis/tableCode/formatCode/sortCode=%s on %s-OPG table' %(self.atfsCode,self.tableName)
            self.skipOES_Element()
        ###
        #print self.obj

    def readOPG_Data_table19(self): # Applied Loads
        isSort1 = self.isSort1()
        if self.numWide==8:  # real/random
            if self.thermal==0:
                self.createTransientObject(self.appliedLoads,appliedLoadsObject) # real
            else:
                raise NotImplementedError(self.codeInformation())
            #self.OUG_RealTable()
            self.readOPGForces()
        elif self.numWide==14:  # real/imaginary or mag/phase
            if self.thermal==0:
                self.createTransientObject(self.appliedLoads,complexAppliedLoadsObject) # complex
            else:
                raise NotImplementedError(self.codeInformation())
            #self.OUG_ComplexTable()
            raise NotImplementedError(self.codeInformation())
        else:
            raise NotImplementedError('only numWide=8 or 14 is allowed  numWide=%s' %(self.numWide))
        ###

    def readOPG_Data_table2(self): # Load Vector
        isSort1 = self.isSort1()
        if self.numWide==8:  # real/random
            if self.thermal==0:
                self.createTransientObject(self.loadVectors,loadVectorObject) # real
            else:
                raise NotImplementedError(self.codeInformation())
            self.OUG_RealTable()
        elif self.numWide==14:  # real/imaginary or mag/phase
            if self.thermal==0:
                self.createTransientObject(self.loadVectors,complexLoadVectorObject) # complex
            else:
                raise NotImplementedError(self.codeInformation())
            self.OUG_ComplexTable()
        else:
            raise NotImplementedError('only numWide=8 or 14 is allowed  numWide=%s' %(self.numWide))
        ###

    def readOPG_Data_table12(self): # Nonlinear Force Vector (in progress)
        isSort1 = self.isSort1()
        if self.numWide==8:  # real/random
            if self.thermal==0:
                self.createTransientObject(self.forceVectors,forceVectorObject) # real
            else:
                raise NotImplementedError(self.codeInformation())
            self.OUG_RealTable()
        elif self.numWide==14:  # real/imaginary or mag/phase
            if self.thermal==0:
                self.createTransientObject(self.forceVectors,complexForceVectorObject) # complex
            else:
                raise NotImplementedError(self.codeInformation())
            self.OUG_ComplexTable()
        else:
            raise NotImplementedError('only numWide=8 or 14 is allowed  numWide=%s' %(self.numWide))
        ###

    def readOPGForces(self): ## @todo needs some work...
        dt = self.nonlinearFactor
        (format1,extract) = self.getOUG_FormatStart()
        format1 += 'i'

        #nTotal = self.numWide*4
        nTotal = 40 # same as dn
        while len(data)>dn:
            #print "len(data) = ",len(data)
            eData = self.data[0:dn]
            #self.printBlock(data[:dn])
            (gridDevice,eid) = unpack(format1,data[0:8])
            nodeID = extract(gridDevice,dt)
            
            source = ''.join(unpack('cccccccc',data[8:16]))
            (dx,dy,dz,rx,ry,rz) = unpack('ffffff',data[16:40])
            #print "source = |%s|" %(source)
            
            #print "nodeID=%s eid=%s source=|%s| dx=%-4i dy=%-4i dz=%-4i rx=%-4i ry=%-4i rz=%-4i" %(nodeID,eid,source,dx,dy,dz,rx,ry,rz)
            source2 = source.replace('*','').replace('-','').strip()
            assert source2.isalnum(),'source=|%s| contains invalid characters...' %(source)
            
            self.obj.add(nodeID,eid,source,dx,dy,dz,rx,ry,rz)
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            #print "nodeID=%g dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(nodeID,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            self.data = self.data[dn:]
        ###
        #print "***********"
        #print self.obj

