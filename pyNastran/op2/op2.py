from fortranFile import FortranFile
from op2Codes import Op2Codes
from op2Errors import *
import os
import sys
from struct import unpack

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.bdf_helper import getMethods,addMethods,writeMesh
from pyNastran.op2.tables.resultTable import ResultTable
from pyNastran.op2.tables.geom.geometryTables import GeometryTables

class Op2(BDF,
#class Op2(getMethods,addMethods,writeMesh, # BDF methods
          FortranFile,Op2Codes,GeometryTables,ResultTable):

    def __init__(self,op2FileName):
        BDF.__init__(self,debug=True,log=None)
        bdfExtension = '.bdf'
        f06Extension = '.f06'
        (fname,extension) = os.path.splitext(op2FileName)
        
        print "fname=%s ext=%s" %(fname,extension)
        
        self.op2FileName = op2FileName
        self.bdfFileName = fname+bdfExtension
        self.f06FileName = fname+f06Extension
        print "bdfFileName = ",self.bdfFileName
        self.stopCode = False
        self.makeOp2Debug = False
        self.skippedCardsFile = open('skippedCards.out','a')

        self.tablesToRead = ['GEOM1','GEOM2','GEOM3','GEOM4', # nodes/geometry/loads/BCs
                             'EPT','MPT','MPTS', # properties/materials
                             'DYNAMIC','DYNAMICS',
                             'DIT',  # some header table...

                             'DESTAB',                # design variables
                             'OQG1','OQGV1','OQMG1',  # spc/mpc forces
                             
                             'OUGV1',                 # displacements
                             'OUGPSD2','OUGATO2','OUGRMS2','OUGNO2','OUGCRM2',
                             
                             'OEF1X','DOEF1','OEFIT', # applied forces
                             'OGPFB1','OGS1',         # grid point forces/stresses
                             
                             'OES1X','OES1X1','OES1C',      # stress
                             'OSTR1C','OSTR1X',             # strains
                             'OESNLXR','OESNLXD','OESNL1X','OESTRCP', # nonlinear stress
                             'OESCP',                       # cylinder stress???
                             'OESTRCP',                     # cylinder strain???
                             'OESRT',                       # rotational stress?

                             'ONRGY1', # energy
                             'R1TABRG','HISADD',  # SOL 200

                             ## @todo what do these do???
                             'OPG1','OPGV1', # think this is an OUG table...
                             'OPNL1',
                             'OUPV1',
                             'VIEWTB','ERRORN',
                             'OFMPF2M','OSMPF2M','OPMPF2M','OGPMPF2M','OLMPF2M','OPGPSD2',
                             'PCOMPTS',
                             # OMNS
                             'OMM2',
                             
                             # new
                             #'OUGCRM2','OUGNO2','OUGRMS2','OUGATO2','OUGPSD2''OMM2','AGRF','AFRF','AEMONPT','PERF','PMRF','MONITOR','SDF','FOL','STDISP'
                             #'OVGNO2','OVGRMS2','OVGATO2','OVGPSD2'
                             'STDISP','SDF','MONITOR','PMRF','PERF','PFRF','AEMONPT','AFRF','AGRF',
                             #'FOL',
                             ]
                             
        ## GEOM1 & GEOM2 are skippable on simple problems...hmmm

        self.iSubcaseNameMap = {}
        self.tableNames = []


        # OUG - displacement
        self.displacements = {}           # aCode=1 tCode=1 fCode=1 sortCode=0 thermal=0
        self.temperatures  = {}           # aCode=1 ------- ------- sortCode=0 thermal=1
        self.nonlinearDisplacements  = {} # aCode=6 ------- fCode=1 sortCode=0 thermal=0
        self.nonlinearTemperatures   = {} # ------- ------- ------- ---------- thermal=1

        self.eigenvectors = {}            # aCode=2 tCode=7 ------- sortCode=1 thermal=0

        # OUG - velocity
        self.velocities = {}              # aCode=6 tCode=10 fCode=3 sortCode=0 thermal=0

        # OUG - acceleration
        self.accelerations = {}           # aCode=6 tCode=11 fCode=3 sortCode=0 thermal=0

        # OEF
        ## rename to staticLoads/thermalLoads
        self.forces = {}
        self.fluxes = {}
        self.temperatureForces = {}       # aCode=1  tCode=4 fCode=1 sortCode=0 thermal=1
    
        self.modalForces = {}
        
        ## rename to nonlinearStaticLoads/nonlinearThermalLoads ???
        self.nonlinearForces = {}         # aCode=10 tCode=4 fCode=1 sortCode=0 thermal=0
        self.nonlinearFluxes = {}         # aCode=10 tCode=4 fCode=1 sortCode=0 thermal=1

        # OES
        self.celasStress   = {}
        self.celasStrain   = {}

        self.rodStress    = {}
        self.rodStrain    = {}
        self.conrodStress = {}
        self.conrodStrain = {}

        self.barStress  = {}
        self.barStrain  = {}
        self.beamStress = {}
        self.beamStrain = {}

        self.plateStress = {}
        self.plateStrain = {}
        self.solidStress = {}
        self.solidStrain = {}
        self.compositePlateStress = {}
        self.compositePlateStrain = {}
        

        # OQG
        self.spcForces      = {}
        self.modalSPCForces = {}

        self.mpcForces      = {}
        self.modalMPCForces = {}
        
        # OPG
        self.appliedLoads = {}
        
        # OEE
        self.strainEnergy = {}

    def printResults(self):
        results = [
                   # OUG - Displacements/Velocity/Acceleration/Temperature/Heat Flux/
                   #       SPC Forces
                   self.displacements,self.temperatures,
                   self.eigenvectors,
                   self.nonlinearTemperatures,self.nonlinearDisplacements,
                   #self.forces,self.fluxes,
                   
                   # OEF - Applied Forces/Temperatures
                   self.nonlinearForces,self.nonlinearFluxes,
                   self.temperatureForces,
                   
                   # OQG1 - Forces
                   self.spcForces,self.mpcForces,
                   
                   # OES - Stress/Strain
                   self.celasStress,self.celasStrain,
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

    def readOp2(self):
        self.op2 = open(self.op2FileName,'rb')
        
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
                #print "startTell = ",self.op2.tell()
                if tableName=='GEOM1': # nodes,coords,etc.
                    self.readTable_Geom1()
                elif tableName=='GEOM1N':
                    self.readTable_Geom1N()
                elif tableName=='GEOM2': # elements
                    self.readTable_Geom2()
                elif tableName=='GEOM3': # static/thermal loads
                    self.readTable_Geom3()
                elif tableName=='GEOM4': # constraints
                    self.readTable_Geom4()

                elif tableName in ['EPT','EPTS']:  # element properties
                    self.readTable_EPT()
                elif tableName in ['MPT','MPTS']:  # material properties
                    self.readTable_MPTS()
                elif tableName in ['DYNAMIC','DYNAMICS']:  # dyanmic info
                    self.readTable_DYNAMICS()
                elif  tableName=='DIT':  # tables...TABLED1/TABLEM1/TABLES1/GUST
                    self.readTable_DIT()

                elif tableName in ['DESTAB']:  # design variable table
                    self.readTable_DesTab()
                
                elif tableName in ['R1TABRG']: # not done - response table
                    self.readTable_R1TAB()
                    self.isOptimization = True

                elif tableName in ['HISADD']: # not done
                    self.readTable_R1TAB()
                    self.isOptimization = True
                elif tableName in ['ERRORN']: # not done
                    self.readTable_R1TAB()
                elif tableName in ['VIEWTB','STDISP','FOL','OMM2']: # not done
                    self.readTable_R1TAB()

                elif tableName in ['OPG1','OGPFB1','OPNL1','OGS1','OPGV1']: # table of applied loads
                    self.readTable_OGP1()
                elif tableName in ['OFMPF2M','OSMPF2M','OPMPF2M','OGPMPF2M','OLMPF2M','OPGPSD2']: # what are these???
                    self.readTable_OGP1()
                elif tableName in ['SDF']: # no idea
                    self.readTable_OGP1()

                
                elif tableName in ['OEF1X','DOEF1','OEFIT']:  # applied loads
                    self.readTable_OEF1()
                elif tableName in ['OQG1','OQMG1','OQGV1']:  # spc/mpc forces
                    self.readTable_OQG1()

                elif tableName in ['OUGV1','OUPV1']: # displacements/velocity/acceleration
                    self.readTable_OUG1()
                elif tableName in ['OUGPSD2','OUGATO2','OUGRMS2','OUGNO2','OUGCRM2']: # OUG tables???
                    self.readTable_OUG1()

                elif tableName in ['OES1X','OES1X1','OSTR1X','OES1C','OESNLXR','OESNLXD','OESNL1X','OESCP','OESRT']: # stress
                    self.readTable_OES1()
                elif tableName in ['OSTR1X','OSTR1C','OESTRCP']: # strain
                    self.readTable_OES1()

                elif tableName in ['ONRGY1']: # energy
                    self.readTable_OEE1()
                elif tableName in ['PCOMPTS','MONITOR','PMRF','PERF','PFRF','AEMONPT','FOL','AFRF','AGRF',]:
                    self.readTable_PCOMPTS() # 'SDF',
                else:
                    raise Exception('unhandled tableName=|%s|' %(tableName))
                #print "endTell   = ",self.op2.tell()
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
        self.skippedCardsFile.close()

        #self.printSection(4*51+12)
        
    def parseApproachCode(self,data):
        """
        int3 is the 3rd word in table=-3 and may be 
        elementType or something else depending on the table
        """
        #self.printBlock(data)
        (aCode,tCode,int3,iSubcase) = unpack('iiii',data[:16])
        self.iSubcase = iSubcase
        self.tableCode = tCode%1000
        self.sortCode  = tCode/1000
        self.deviceCode   = aCode%10
        self.analysisCode = (aCode-self.deviceCode)/10
        #print "aCode(1)=%s analysisCode=%s deviceCode=%s tCode(2)=%s tableCode=%s sortCode=%s elementType(3)=%s iSubcase(4)=%s" %(aCode,self.analysisCode,self.deviceCode,tCode,self.tableCode,self.sortCode,elementType,self.iSubcase)
        print self.printTableCode(self.tableCode)
        return (int3)

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
        
    def deleteAttributes(self,params):
        params += ['dataCode','deviceCode','analysisCode','tableCode','iSubcase','data','numWide','nonlinearFactor','obj']
        for param in params:
            if hasattr(self,param):
                #print '%s = %s' %(param,getattr(self,param))
                delattr(self,param)

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

        self.readHollerith() # not really a hollerith, just the end of the block (so bufferWords*4)
        
        self.Title = Title.strip()
        if self.iSubcase not in self.iSubcaseNameMap:
            self.iSubcaseNameMap[self.iSubcase] = [Subtitle,Label]

    def tableInit(self,word):
        self.tableName = word.strip()
        self.tableNames.append(word)
        msg = '*'*20+word+'*'*20+'\n'
        if self.makeOp2Debug:
            self.op2Debug.write(msg)

    def getTableNamesFromOP2(self):
       return self.tableNames
