from fortranFile import FortranFile
from op2Codes import Op2Codes
from op2Errors import *
import os
import sys
import struct
from struct import unpack

from geometryTables import GeometryTables
from elementsStressStrain import ElementsStressStrain
from pyNastran.bdf.bdf_helper import getMethods,addMethods,writeMesh
from pyNastran.op2.tables.resultTable import ResultTable

class Op2(getMethods,addMethods,writeMesh, # BDF methods
          FortranFile,Op2Codes,GeometryTables,ElementsStressStrain,ResultTable):

    def bdfInit(self,log=None):
        if log is None:
            from pyNastran.general.logger import dummyLogger
            loggerObj = dummyLogger()
            log = loggerObj.startLog('debug') # or info
        self.log = log

        self.makeOp2Debug = False

        self.sol = None
        self.iSolLine = None
        self.executiveControlLines = []
        self.caseControlDeck = None
        
        self.params  = {}
        self.nodes   = {}
        self.gridSet = None

        self.elements   = {}
        self.properties = {}
        self.materials  = {}

        self.coords = {}


        self.loads   = {}
        self.dareas  = {}
        self.nlparms = {}

        self.flfacts  = {}
        self.aeros    = {}
        self.gusts    = {}
        self.flutters = {}
        self.splines  = {}
        self.caeros   = {}
        self.gravs    = {}

        self.phbdys = {}
        self.convectionProperties = {}
        self.bcs = {}
        self.constraints = {}
        self.suports = {}
        self.spcObject = None
        self.mpcObject = None
        self.dconstrs = {}
        self.desvars = {}
        self.ddvals = {}
        self.rejectCards = []
        self.rejects = []


    def __init__(self,infileName):
        self.bdfInit()

        self.infilename = infileName
        #self.tablesToRead = ['GEOM1','GEOM2','GEOM3','GEOM4','OQG1','OUGV1','OES1X1']  # 'OUGV1','GEOM1','GEOM2'
        #self.tablesToRead = ['GEOM1','GEOM2','OQG1','OUGV1','OES1X1']  # 'OUGV1','GEOM1','GEOM2'
        #self.tablesToRead = ['GEOM1','GEOM2','GEOM3','OQG1','OUGV1','OES1X1']  # 'OUGV1','GEOM1','GEOM2'
        #self.tablesToRead = ['OQG1','OUGV1','OEF1X','OES1X1','OSTR1X','OES1C','OSTR1C','OGPFB1']  # 'OUGV1','GEOM1','GEOM2'
        self.tablesToRead = ['GEOM1','GEOM2','GEOM3','GEOM4',
                             'EPT','MPT','MPTS',
                             'DYNAMIC','DYNAMICS',
                             'DIT',

                             'DESTAB',
                             'OQG1',
                             'OUGV1','OUPV1',
                             'OEF1X','DOEF1',
                             'OPG1','OGPFB1',
                             'OES1X','OES1X1','OSTR1X','OES1C','OSTR1C','OESNLXR','OESNLXD',
                             'ONRGY1',
                             
                             #what is OUPV1
                             ]
                             
        ## GEOM1 & GEOM2 are skippable on simple problems...hmmm

        self.iSubcaseNameMap = {}


        # OUG
        self.displacements = {}           # aCode=1 tCode=1 fCode=1 sortCode=0 thermal=0
        self.temperatures  = {}           # aCode=1 tCode=1 fCode=1 sortCode=0 thermal=1

        self.freqDisplacements = {}       # aCode=5 tCode=1 fCode=3 sortCode=1 thermal=0

        self.nonlinearDisplacements  = {} # aCode=6 tCode=1 fCode=1 sortCode=0 thermal=0
        self.nonlinearTemperatures   = {} # aCode=6 tCode=1 fCode=1 sortCode=0 thermal=1
        self.preBucklingDisplacements= {} # aCode=7 tCode=1 fCode=1 sortCode=0 thermal=0

        self.eigenvectors = {}            # aCode=2 tCode=7 fCode=1 sortCode=1 thermal=0
        self.postBucklingEigenvector = {} # aCode=8 tCode=7 fCode=1 sortCode=1 thermal=0
        self.complexEigenvalues = {}      # aCode=9 tCode=7 fCode=1 sortCode=1 thermal=0

        #self.forces = {}
        #self.fluxes = {}

        # OEF
        
        ## rename to staticLoads/thermalLoads
        self.displacementForces = {}      # aCode=1  tCode=4 fCode=1 sortCode=0 thermal=0
        self.temperatureForces = {}       # aCode=1  tCode=4 fCode=1 sortCode=0 thermal=1

        self.bucklingForces = {}          # aCode=2  tCode=4 fCode=1 sortCode=0 thermal=0

        ## rename to complexEigenvalueLoads ???
        self.complexEigenvalueForces = {} # aCode=9  tCode=4 fCode=2 sortCode=1 thermal=0
        
        ## rename to nonlinearStaticLoads/nonlinearThermalLoads ???
        self.nonlinearForces = {}         # aCode=10 tCode=4 fCode=1 sortCode=0 thermal=0
        self.nonlinearFluxes = {}         # aCode=10 tCode=4 fCode=1 sortCode=0 thermal=1

        # OES
        self.conrodStress   = {}
        self.conrodStrain   = {}

        self.rodStress   = {}
        self.rodStrain   = {}
        self.barStress   = {}
        self.barStrain   = {}
        self.plateStress = {}
        self.plateStrain = {}
        self.solidStress = {}
        self.solidStrain = {}
        self.compositePlateStress = {}
        self.compositePlateStrain = {}


        # OQG
        self.spcForces           = {} # aCode=1  tCode=3 fCode=1 sortCode=0 thermal=0
        self.spcBucklingForces   = {} # aCode=2  tCode=3 fCode=1 sortCode=0 thermal=0
        self.realImagConstraints = {} # aCode=10 tCode=? fCode=1 sortCode=1 thermal=?
        
        # OPG
        self.appliedLoads = {}  # aCode=1 tCode=2 fCode=1 sortCode=0 thermal=0
        
        # OEE
        self.strainEnergy      = {} # aCode=1 tCode=18 fCode=1 sortCode=0
        self.modesStrainEnergy = {} # aCode=2 tCode=18 fCode=1 sortCode=0

    def printResults(self):
        results = [
                   # OUG - Displacements/Velocity/Acceleration/Temperature/Heat Flux/
                   #       SPC Forces
                   self.displacements,self.temperatures,
                   self.eigenvectors,
                   self.nonlinearTemperatures,self.nonlinearDisplacements,
                   self.forces,self.fluxes,
                   
                   # OEF - Applied Forces/Temperatures
                   self.nonlinearForces,self.nonlinearFluxes,
                   self.temperatureForces,
                   
                   # OQG1 - 
                   self.spcForces,
                   
                   # OES - Stress/Strain
                   self.rodStress,self.rodStrain,
                   self.barStress,self.barStrain,
                   self.plateStress,self.plateStrain,
                   self.solidStress,self.solidStrain,
                   self.compositePlateStress,self.compositePlateStrain,
                   
                   # OGP - Applied Force/Moment
                   self.appliedLoads,
                   ]
        
        msg = '---ALL RESULTS---\n'
        for result in results:
            for iSubcase,res in sorted(result.items()):
                msg += 'iSubcase = %s\n' %(iSubcase)
                msg += str(res) + '\n'
            ###
        ###
        return msg
        
    def readTapeCode(self):
        #self.printSection(500)
        #sys.exit('op2-readTapeCode')
        self.readMarkers([3])
        #self.printSection(20)
        ints = self.readIntBlock()
        if self.makeOp2Debug:
            self.op2Debug.write('%s\n' %(str(ints)))
        #print "*ints = ",ints
        self.readMarkers([7])

        word = self.readStringBlock() # Nastran Fort Tape ID Code - 
        #print "word = |%r|" %(word)

        self.readMarkers([2])
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-1])

        #data = self.getData(60)
        #self.printBlock(data)

    def read(self):
        self.op2 = open(self.infilename,'rb')
        
        if self.makeOp2Debug:
            self.op2Debug = open('debug.out','wb')
        
        self.n = self.op2.tell()
        #print "self.n = ",self.n
        self.readTapeCode()
        #sys.exit('end of tape code')

        isAnotherTable = True
        while isAnotherTable:
            print '-'*80
            try:
                tableName = self.readTableName(rewind=True,stopOnFailure=False)
            except EndOfFileError:  # the isAnotherTable method sucks...
                isAnotherTable = False
                print "***ok exit, but it could be better..."
                break
            except InvalidMarkersError:  # the isAnotherTable method sucks...
                isAnotherTable = False
                print "***poor exit, but it worked..."
                #raise
                break
            except:
                raise
            print "tableName = |%r|" %(tableName)
 
            if tableName in self.tablesToRead:
                if tableName=='GEOM1': # nodes,coords,etc.
                    self.readTable_Geom1()
                elif tableName=='GEOM2': # elements
                    self.readTable_Geom2()
                elif tableName=='GEOM3': # static/thermal loads
                    self.readTable_Geom3()
                elif tableName=='GEOM4': # constraints
                    self.readTable_Geom4()

                elif tableName=='EPT':   # element properties
                    self.readTable_EPT()
                elif tableName in ['MPT','MPTS']:  # material properties
                    self.readTable_MPTS()
                elif tableName in ['DYNAMIC','DYNAMICS']:  # dyanmic info
                    self.readTable_DYNAMICS()
                elif  tableName=='DIT':  # tables...TABLED1/TABLEM1/TABLES1/GUST
                    self.readTable_DIT()

                elif tableName=='DESTAB':  # design variable table
                    self.readTable_DesTab()

                elif tableName in ['OPG1','OGPFB1']: # table of applied loads
                    self.readTable_OGP1()

                
                elif tableName in ['OEF1X','DOEF1']:  # applied loads
                    self.readTable_OEF1()
                elif tableName=='OQG1':  # spc forces
                    self.readTable_OQG1()
                elif tableName in ['OUGV1','OUPV1']: # displacements/velocity/acceleration
                    self.readTable_OUG1()
                elif tableName in ['OES1X','OES1X1','OSTR1X','OES1C','OSTR1C','OESNLXR','OESNLXD']: # stress/strain
                    self.readTable_OES1()
                elif tableName in ['ONRGY1']: # energy???
                    self.readTable_OEE1()
                else:
                    raise Exception('unhandled tableName=|%s|' %(tableName))
                #print "---isAnotherTable---"
                (isAnotherTable) = self.hasMoreTables()
                #isAnotherTable = True
            else:
                (isAnotherTable) = self.skipNextTable()
                continue
            #print self.printSection(140)
            print "*** finished tableName = |%r|" %(tableName)
            ###
        ###
        #

        #tableName = self.readTableName(rewind=True)
        #self.readTable_Geom2()
        #self.readTable_Geom3()
        #self.readTable_Geom4()
        #self.readTable_OQG1()
        #self.readTable_OES1X1()
        print "---end of all tables---"

        #self.printSection(4*51+12)
        
    def parseApproachCode(self,data):
        #self.printBlock(data)
        (aCode,tCode,elementType,iSubcase) = unpack('iiii',data[:16])
        self.iSubcase = iSubcase
        self.tableCode = tCode%1000
        self.sortCode = tCode/1000
        self.deviceCode   = aCode%10
        self.approachCode = (aCode-self.deviceCode)/10
        print "aCode(1)=%s analysisCode=%s deviceCode=%s tCode(2)=%s tableCode=%s sortCode=%s elementType(3)=%s iSubcase(4)=%s" %(aCode,self.approachCode,self.deviceCode,tCode,self.tableCode,self.sortCode,elementType,self.iSubcase)
        print self.printTableCode(self.tableCode)
        return (elementType)

    def getValues(self,data,sFormat,iWordStart,iWordStop=None):
        """
        extracts the ith word from the data structure as the provided type
        supports multiple inputs with iWordStop (note this is words, not outputs)
        @warning
            works with nastran syntax, not standard python syntax
            this makes it work with what's documented in the DMAP manual
        """
        if iWordStop==None:
            #print "iWordStart=%s data[%s:%s]" %(iWordStart,iWordStart*4,(iWordStart+1)*4)
            ds = data[(iWordStart-1)*4:iWordStart*4]
            return unpack(sFormat,ds)[0]
            
        #print "type(data) = ",type(data)
        ds = data[(iWordStart-1)*4:(iWordStop-1)*4]
        return unpack(sFormat,ds)
        
    def getValues8(self,data,sFormat,iWordStart,iWordStop=None):
        raise Exception('is this used...')
        if iWordStop==None:
            ds = data[iWordStart*8:(iWordStart+1)*8]
            return unpack(sFormat,ds)[0]
            
        #print "type(data) = ",type(data)
        ds = data[iWordStart*8:iWordStop*8]
        return unpack(sFormat,ds)
        
    def deleteAttributes(self,params):
        params += ['deviceCode','approachCode','tableCode','iSubcase','data','numWide']
        for param in params:
            if hasattr(self,param):
                #print '%s = %s' %(param,getattr(self,param))
                delattr(self,param)

    def createTransientObject(self,storageObj,classObj,dt):
        """@note dt can also be loadStep depending on the class"""
        if self.iSubcase in storageObj:
            self.obj = storageObj[self.iSubcase]
            self.obj.updateDt(dt)
        else:
            self.obj = classObj(self.iSubcase,dt)
        ###

    def getBufferWords(self):
        bufferWords = self.getMarker()
        print "buffMarker = |%s|" %(bufferWords)
        print "bufferWords = ",bufferWords,bufferWords*4
        assert bufferWords >0
        return bufferWords

    def verifyBufferSize(self,bufferWords):
        assert bufferWords>0,self.printSection(220)

    def readTitle(self):
        word = self.readString(384) # titleSubtitleLabel
        Title    = word[0:128].strip()
        Subtitle = word[128:256].strip()
        Label    = word[256:].strip()
        #print "Title    %s |%s|" %(len(Title   ),Title)
        #print "Subtitle %s |%s|" %(len(Subtitle),Subtitle)
        #print "Label    %s |%s|" %(len(Label   ),Label)
        #print "Title    %s |%s|" %(len(Title   ),Title.strip())
        #print "Subtitle %s |%s|" %(len(Subtitle),Subtitle.strip())
        #print "Label    %s |%s|" %(len(Label   ),Label.strip())
        #print "Title    |%s|" %(Title)
        #print "Subtitle |%s|" %(Subtitle)
        #print "Label    |%s|" %(Label)

        self.readHollerith()
        
        self.Title = Title.strip()
        if self.iSubcase not in self.iSubcaseNameMap:
            self.iSubcaseNameMap[self.iSubcase] = [Subtitle,Label]

    def tableInit(self,word):
        msg = '*'*20+word+'*'*20+'\n'
        if self.makeOp2Debug:
            self.op2Debug.write(msg)
