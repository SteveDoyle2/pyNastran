import sys
import copy
from struct import unpack

# pyNastran
#from pyNastran.op2.resultObjects.ougv1_Objects import (
#    temperatureObject,
#    nonlinearTemperatureObject,
#    fluxObject,nonlinearFluxObject)
from pyNastran.op2.resultObjects.opg_Objects import appliedLoadsObject


class OGP(object):
    """Table of element forces"""
    def readTable_OGP1(self):
        self.tableName = 'OGP'
        table3     = self.readTable_OGP_3
        table4Data = self.readOGP1_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OGP()

    def deleteAttributes_OGP(self):
        params = ['lsdvm','mode','eigr','eign','eigi','modeCycle','freq','time','lftsfq','dLoadID','formatCode','numWide','oCode']
        self.deleteAttributes(params)

    def readTable_OGP_3(self,iTable): # iTable=-3
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

        self.dLoadID    = self.getValues(data,'i',8)  ## dynamic load set ID/random code
        self.formatCode = self.getValues(data,'i',9)  ## format code
        self.numWide    = self.getValues(data,'i',10) ## number of words per entry in record; @note is this needed for this table ???
        self.oCode      = self.getValues(data,'i',11) ## undefined in DMAP...
        self.thermal    = self.getValues(data,'i',23) ## thermal flag; 1 for heat ransfer, 0 otherwise

        self.dataCode = {'analysisCode': self.approachCode,'deviceCode':self.deviceCode,
                         'dLoadID':self.dLoadID,'formatCode':self.formatCode,
                         'numWide': self.numWide,'oCode':self.oCode,
                         'thermal': self.thermal}
        print "dLoadID(8)=%s formatCode(9)=%s numWide(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.formatCode,self.numWide,self.oCode,self.thermal)
        
        ## assuming tCode=1
        if self.approachCode==1:   # statics
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
            self.nonlinearFactor = self.lsdvmn
        elif self.approachCode==2: # normal modes/buckling (real eigenvalues)
            self.mode      = self.getValues(data,'i',5) ## mode number
            self.eign      = self.getValues(data,'f',6) ## real eigenvalue
            self.modeCycle = self.getValues(data,'f',7) ## mode or cycle @todo confused on the type ???
            self.nonlinearFactor = self.mode
        elif self.approachCode==3: # differential stiffness 0
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
            self.nonlinearFactor = self.lsdvmn
        elif self.approachCode==4: # differential stiffness 1
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
            self.nonlinearFactor = self.lsdvmn
        elif self.approachCode==5:   # frequency
            self.freq = self.getValues(data,'f',5) ## frequency
            self.nonlinearFactor = self.freq

        elif self.approachCode==6: # transient
            self.time = self.getValues(data,'f',5) ## time step
            self.nonlinearFactor = self.time
            print "TIME(5)=%s" %(self.time)
        elif self.approachCode==7: # pre-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## load set
            self.nonlinearFactor = self.lsdvmn
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.approachCode==8: # post-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## mode number
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.nonlinearFactor = self.lsdvmn
            print "LSDVMN(5)=%s  EIGR(6)=%s" %(self.lsdvmn,self.eigr)
        elif self.approachCode==9: # complex eigenvalues
            self.mode   = self.getValues(data,'i',5) ## mode
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.eigi   = self.getValues(data,'f',7) ## imaginary eigenvalue
            self.nonlinearFactor = self.mode
            print "LFTSFQ(5)=%s  EIGR(6)=%s  EIGI(7)=%s" %(self.lftsfq,self.eigr,self.eigi)
        elif self.approachCode==10: # nonlinear statics
            self.lftsfq = self.getValues(data,'f',5) ## load step
            self.nonlinearFactor = self.lftsfq
            print "LFTSFQ(5) = %s" %(self.lftsfq)
        elif self.approachCode==11: # old geometric nonlinear statics
            self.lsdvmn = self.getValues(data,'i',5)
            self.nonlinearFactor = self.lsdvmn
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.approachCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.lsdvmn = self.getValues(data,'i',5)
            self.nonlinearFactor = self.lsdvmn
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        else:
            raise RuntimeError('invalid approach code...approachCode=%s' %(self.approachCode))

        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set, Mode number
        
        #print "*iSubcase=%s"%(self.iSubcase)
        #print "approachCode=%s tableCode=%s thermal=%s" %(self.approachCode,self.tableCode,self.thermal)
        #print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

    def readOGP1_Data(self):
        print "self.approachCode=%s tableCode(1)=%s thermal(23)=%g" %(self.approachCode,self.tableCode,self.thermal)
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        self.atfsCode = [self.approachCode,self.tableCode,self.formatCode,self.sortCode]
        print "tfsCode=%s" %(tfsCode)

        # grid point force balance
        if   tfsCode==[19,1,0]:
            self.readOGP1_Data_format1_sort0()
        #elif tfsCode==[19,1,1]:
        #    self.readOGP1_Data_format1_sort1()
        #elif tfsCode==[19,2,1]:
        #    self.readOGP1_Data_format2_sort1()
        #elif tfsCode==[19,3,0]:
        #    self.readOGP1_Data_format3_sort0()
        #elif tfsCode==[19,3,1]:
        #    self.readOGP1_Data_format3_sort1()
        
        # load vector
        #elif tfsCode==[2,1,0]:
        #    self.readOGP1_Data_format1_sort0()
        #elif tfsCode==[2,1,1]:
        #    self.readOGP1_Data_format1_sort1()
        #elif tfsCode==[2,2,1]:
        #    self.readOGP1_Data_format2_sort1()
        #elif tfsCode==[2,3,0]:
        #    self.readOGP1_Data_format3_sort0()
        #elif tfsCode==[2,3,1]:
        #    self.readOGP1_Data_format3_sort1()

        # Nonlinear force vector
        #elif tfsCode==[12,1,0]:
        #    self.readOGP1_Data_format1_sort0()

        # OGS1- grid point stresses - surface
        #elif tfsCode==[26,1,0]:
        #    self.readOGP1_Data_format1_sort0()

        # OGS1- grid point stresses - volume direct
        #elif tfsCode==[27,1,0]:
        #    self.readOGP1_Data_format1_sort0()

        # OGS1- grid point stresses - principal
        #elif tfsCode==[28,1,0]:
        #    self.readOGP1_Data_format1_sort0()

        # OGS - Grid point stress discontinuities (plane strain)
        #elif tfsCode==[35,1,0]:
        #    self.readOGP1_Data_format1_sort0()
        

        # OFMPF2M - does this belong here?
        #elif tfsCode==[51,3,3]:
        #    self.readOGP1_Data_format3_sort3()
        
        # OSMPF2M - does this belong here?
        #elif tfsCode==[52,3,3]:
        #    self.readOGP1_Data_format3_sort3()
        
        # OPMPF2M - does this belong here?
        #elif tfsCode==[53,3,3]:
        #    self.readOGP1_Data_format3_sort3()
        
        # OLMPF2M - does this belong here?
        #elif tfsCode==[54,3,3]:
        #    self.readOGP1_Data_format3_sort3()

        # OGMPF2M - does this belong here?
        #elif tfsCode==[55,3,3]:
        #    self.readOGP1_Data_format3_sort3()
        else:
            #raise Exception('bad tableCode/formatCode/sortCode=%s on OGP table' %(self.atfsCode))
            print 'bad tableCode/formatCode/sortCode=%s on OGP table' %(self.atfsCode)
            self.skipOES_Element(None)
        ###
        #print self.obj
        print "-------finished OGP----------"

    def readOGP1_Data_format1_sort1(self):
        print 'not supported OGP solution...atfsCode=%s' %(self.atfsCode)
        self.skipOES_Element(None)

    def readOGP1_Data_format2_sort1(self):
        print 'not supported OGP solution...atfsCode=%s' %(self.atfsCode)
        self.skipOES_Element(None)

    def readOGP1_Data_format3_sort0(self):
        print 'not supported OGP solution...atfsCode=%s' %(self.atfsCode)
        self.skipOES_Element(None)

    def readOGP1_Data_format3_sort1(self):
        print 'not supported OGP solution...atfsCode=%s' %(self.atfsCode)
        self.skipOES_Element(None)

    def readOGP1_Data_format3_sort3(self):
        print 'not supported OGP solution...atfsCode=%s' %(self.atfsCode)
        self.skipOES_Element(None)

    def readOGP1_Data_format1_sort0(self):
        if self.thermal==0:
            if self.approachCode==1: # displacement
                print "isAppliedLoads"
                self.obj = appliedLoadsObject(self.dataCode,self.iSubcase)
                self.appliedLoads[self.iSubcase] = self.obj
                self.readOGPForces(self.data,self.obj)
            else:
                self.skipOES_Element(None)
                print 'not supported OGP solution...atfsCode=%s' %(self.atfsCode)
                #raise Exception('not supported OGP solution...')
            ###
        elif self.thermal==1:
            self.skipOES_Element(None)
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###

    def readOGPForces(self,data,scalarObject):
        deviceCode = self.deviceCode
        #self.printBlock(data[0:self.numWide*4])
        dn = self.numWide*4
        while len(data)>dn:
            #print "len(data) = ",len(data)
            eData = data[0:dn]
            #self.printBlock(data[:dn])
            gridDevice,eid = unpack('ii',data[0:8])
            nodeID = (gridDevice-deviceCode)/10
            
            source = ''.join(unpack('cccccccc',data[8:16]))
            #self.printBlock(data[16:dn])
            (dx,dy,dz,rx,ry,rz) = unpack('ffffff',data[16:40])
            #print "source = |%s|" %(source)
            #print "len(data[8:40]"
            
            print "nodeID=%s eid=%s source=|%s| dx=%-4i dy=%-4i dz=%-4i rx=%-4i ry=%-4i rz=%-4i" %(nodeID,eid,source,dx,dy,dz,rx,ry,rz)
            scalarObject.add(nodeID,eid,source,dx,dy,dz,rx,ry,rz)
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-self.deviceCode)/10
            #print "grid=%g dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            #print type(scalarObject)
            data = data[dn:]
            #sys.exit('asd')
        ###
        print "***********"
        #print self.obj
        #sys.exit('check...')

