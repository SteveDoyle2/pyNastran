import os
import sys
from struct import unpack

from fortranFile import FortranFile
from op2Codes import Op2Codes
from op2Errors import *

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.bdf_helper import getMethods,addMethods,writeMesh
from pyNastran.op2.tables.resultTable import ResultTable
from pyNastran.op2.tables.geom.geometryTables import GeometryTables

class OP2(BDF,
#class Op2(getMethods,addMethods,writeMesh, # BDF methods
          FortranFile,Op2Codes,GeometryTables,ResultTable):

    def setSubcases(self,iSubcases=[]):
        """
        allows you to read only the subcases in the list of iSubcases
        @param iSubcases list of [subcase1_ID,subcase2_ID]
        @note  the default is all the subcases
        """
        ## stores the set of all subcases that are in the OP2
        self.subcases = set()
        if iSubcases==[]:
            ## stores if the user entered [] for iSubcases
            self.isAllSubcases = True
            self.validSubcases = []
        else:
            ## should all the subcases be read (default=True)
            self.isAllSubcases = False
            ## the set of valid subcases -> set([1,2,3])
            self.validSubcases = set(iSubcases)
        ###
        self.log.debug("setSubcases - iSubcases = %s" %(self.validSubcases))

    def isValidSubcase(self):
        """
        lets the code check whether or not to read a subcase
        """
        if not self.isAllSubcases:
            if self.iSubcase in self.validSubcases:
                return True
            return False
        return True

    def __init__(self,op2FileName,makeGeom=False,debug=True,log=None):
        """
        Initializes the Op2 object
        @param op2FileName the file to be parsed
        @param makeGeom    reads the BDF tables (default=False)
        @param debug       prints data about how the OP2 was parsed (default=False)
        @param log         a logging object to write debug messages to (@see import logging)
        """
        BDF.__init__(self,debug=debug,log=log)
        self.setSubcases() # initializes the variables
        self.log.debug('op2FileName = %s' %(op2FileName))
        bdfExtension = '.bdf'
        f06Extension = '.f06'
        (fname,extension) = os.path.splitext(op2FileName)
        
        #print "fname=%s ext=%s" %(fname,extension)
        
        ## should the BDF tables be parsed
        self.makeGeom = makeGeom
        
        ## the input OP2 filename
        self.op2FileName = op2FileName
        
        ## the expected BDF filename (guessed)
        self.bdfFileName = fname+bdfExtension
        
        ## the expected F06 filename (guessed)
        self.f06FileName = fname+f06Extension
        #print "bdfFileName = ",self.bdfFileName

        ## developer parameter to write the OP2 is ASCII format
        ## to better understand it
        self.makeOp2Debug = False
        
        ## file object containing the skipped cards
        self.skippedCardsFile = open('skippedCards.out','a')

        ## the list of supported tables (dont edit this)
        self.tablesToRead = ['GEOM1','GEOM2','GEOM3','GEOM4', # nodes/geometry/loads/BCs
                             'EPT','MPT','MPTS', # properties/materials
                             'DYNAMIC','DYNAMICS',
                             
                             'DIT',  # some header table...
                            'BGPDT','EQEXIN','EQEXINS','PVT0','CASECC',#'EDOM',
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
                             'OMM2',
                             'OGPWG','EDOM',

                             # new
                             'OUGCRM2','OUGNO2','OUGRMS2','OUGATO2','OUGPSD2','OMM2',
                             'OVGNO2','OVGRMS2','OVGATO2','OVGPSD2',
                             'STDISP','SDF','MONITOR','PMRF','PERF','PFRF','AEMONPT','AFRF','AGRF',
                             'FOL','GEOM1N',
                             
                             'OESNLXR','OESNL1X','OESPSD2','OESNLBR','OESATO2','OESRMS2','OESNO2','OESCRM2',
                             'OAGPSD2','OAGATO2','OAGRMS2','OAGNO2','OAGCRM2',
                             'OVGCRM2',
                             'OPGATO2','OPGPSD2','OPGPSD2','OPGRMS2','OPGNO2','OPGCRM2',
                             'OQGPSD2','OQGATO2','OQGRMS2','OQGNO2','OQGCRM2',
                             ]
        
        ## a dictionary that maps an integer of the subcaseName to the subcaseID
        self.iSubcaseNameMap = {}
        
        ## list of OP2 tables that were read
        ## mainly for debugging
        self.tableNames = []


        ## OUG - displacement
        self.displacements = {}           # aCode=1 tCode=1 fCode=1 sortCode=0 thermal=0

        ## OUG - temperatures
        self.temperatures  = {}            # aCode=1 ------- ------- sortCode=0 thermal=1
        #self.nonlinearDisplacements  = {} # aCode=6 ------- fCode=1 sortCode=0 thermal=0
        #self.nonlinearTemperatures   = {} # ------- ------- ------- ---------- thermal=1

        ## OUG - eigenvectors
        self.eigenvectors = {}            # aCode=2 tCode=7 ------- sortCode=1 thermal=0

        ## OUG - velocity (not done)
        self.velocities = {}              # aCode=6 tCode=10 fCode=3 sortCode=0 thermal=0

        ## OUG - acceleration (not done)
        self.accelerations = {}           # aCode=6 tCode=11 fCode=3 sortCode=0 thermal=0

        # OEF
        # rename to staticLoads/thermalLoads
        #self.forces = {}
        #self.fluxes = {}
        #self.temperatureForces = {}       # aCode=1  tCode=4 fCode=1 sortCode=0 thermal=1

        #self.modalForces = {}
        
        # rename to nonlinearStaticLoads/nonlinearThermalLoads ???
        #self.nonlinearForces = {}         # aCode=10 tCode=4 fCode=1 sortCode=0 thermal=0
        #self.nonlinearFluxes = {}         # aCode=10 tCode=4 fCode=1 sortCode=0 thermal=1

        # OES
        ## OES - CELAS1/CELAS2/CELAS3/CELAS4
        self.celasStress   = {}
        ## OES - CELAS1/CELAS2/CELAS3/CELAS4
        self.celasStrain   = {}
        
        ## OES - isotropic CROD/CONROD/CTUBE
        self.rodStress  = {}
        ## OES - isotropic CROD/CONROD/CTUBE
        self.rodStrain  = {}
        ## OES - isotropic CBAR
        self.barStress  = {}
        ## OES - isotropic CBAR
        self.barStrain  = {}
        ## OES - isotropic CBEAM
        self.beamStress = {}
        ## OES - isotropic CBEAM
        self.beamStrain = {}

        ## OES - isotropic CTRIA3/CQUAD4
        self.plateStress = {}
        ## OES - isotropic CTRIA3/CQUAD4
        self.plateStrain = {}

        ## OES - isotropic CTETRA/CHEXA/CPENTA
        self.solidStress = {}
        ## OES - isotropic CTETRA/CHEXA/CPENTA
        self.solidStrain = {}
        ## OES - composite CTRIA3/CQUAD4
        self.compositePlateStress = {}
        ## OES - composite CTRIA3/CQUAD4
        self.compositePlateStrain = {}
        
        ## OES - CSHEAR
        self.shearStress = {}
        ## OES - CSHEAR
        self.shearStrain = {}
        

        # OQG
        #self.spcForces      = {}
        #self.modalSPCForces = {}

        #self.mpcForces      = {}
        #self.modalMPCForces = {}
        
        ## OPG - summation of loads for each element
        self.appliedLoads = {}
        
        ## OEE - strain energy density
        self.strainEnergy = {}

    def printResults(self):
        results = [
                   # OUG - Displacements/Velocity/Acceleration/Temperature/Heat Flux/
                   #       SPC Forces
                   #self.displacements,self.temperatures,
                   self.eigenvectors,
                   self.velocities,
                   self.accelerations,
                   #self.nonlinearTemperatures,self.nonlinearDisplacements,
                   #self.forces,self.fluxes,
                   
                   # OEF - Applied Forces/Temperatures
                   #self.nonlinearForces,self.nonlinearFluxes,
                   #self.temperatureForces,
                   
                   # OQG1 - Forces
                   #self.spcForces,self.mpcForces,
                   
                   # OGP - Applied Force/Moment
                   self.appliedLoads,

                   # OES - Stress/Strain
                   #self.celasStress,self.celasStrain,
                   #self.rodStress,self.rodStrain,
                   #self.barStress,self.barStrain,
                   #self.beamStress,self.beamStrain,
                   #self.plateStress,self.plateStrain,
                   #self.solidStress,self.solidStrain,
                   #self.compositePlateStress,self.compositePlateStrain,
                   
                   # OEE - Strain Energy
                   self.strainEnergy,
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
        """
        reads the OP2 header
        @todo whats in this table?
        """
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

    def readOP2(self):
        """
        reads the op2 file
        """
        ## the OP2 file object
        self.op2 = open(self.op2FileName,'rb')
        
        if self.makeOp2Debug:
            ## a developer debug file (largely unsupported)
            self.op2Debug = open('debug.out','wb')
        ## the byte position in the OP2
        self.n = self.op2.tell()

        try:
            self.readTapeCode()
        except:
            raise TapeCodeError('when this happens, the analysis failed...check the F06')
        ###

        isAnotherTable = True
        while isAnotherTable:
            self.log.debug('-'*80)
            try:
                tableName = self.readTableName(rewind=True,stopOnFailure=False)
            except EndOfFileError:  # the isAnotherTable method sucks...
                isAnotherTable = False
                self.log.debug("***ok exit, but it could be better...")
                break
            except InvalidMarkersError:  # the isAnotherTable method sucks...
                isAnotherTable = False
                self.log.debug("***poor exit, but it worked...")
                #raise
                break
            except:
                raise
            self.log.debug("tableName = |%r|" %(tableName))
            if tableName==None:
                break
            elif tableName in self.tablesToRead:
                try:
                    #print "startTell = ",self.op2.tell()
                    if tableName=='GEOM1': # nodes,coords,etc.
                        self.readTable_Geom1()
                    #elif tableName=='GEOM1N':
                    #    self.readTable_Geom1N()
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
                    #elif  tableName in ['DIT']:  # tables...TABLED1/TABLEM1/TABLES1/GUST
                    #    self.readTable_DIT()
                    elif tableName in ['VIEWTB','EQEXIN','EQEXINS','OEFIT','GEOM1N','OMM2','OGPWG',]:
                        self.readTable_DUMMY_GEOM(tableName)

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

                    elif tableName in ['OPG1','OGPFB1','OPNL1','OGS1','OPGV1']: # table of applied loads
                        self.readTable_OGP1()


                    elif tableName in ['OEF1X','DOEF1']:  # applied loads
                        self.readTable_OEF1()
                    elif tableName in ['OQG1','OQMG1','OQGV1']:  # spc/mpc forces
                        self.readTable_OQG1()

                    elif tableName in ['OUGV1','OUPV1']: # displacements/velocity/acceleration
                        self.readTable_OUG1()
                    elif tableName in ['OUGPSD2','OUGATO2','OUGRMS2','OUGNO2','OUGCRM2']: # OUG tables???
                        self.readTable_OUG1()

                    elif tableName in ['OES1X','OES1X1','OSTR1X','OES1C','OESNLXD','OESCP','OESRT','OESRMS2','OESNO2','OESCRM2',]: # stress
                        self.readTable_OES1()  # 'OESNLXR','OESNL1X'
                    elif tableName in ['OSTR1X','OSTR1C','OESTRCP']: # strain
                        self.readTable_OES1()

                    elif tableName in ['ONRGY1']: # energy
                        self.readTable_OEE1()

                    elif tableName in ['PCOMPTS']:
                        self.readTable_PCOMPTS() # 'SDF',

                    # not done
                    elif tableName in []:
                        self.readTableB_DUMMY()
                    elif tableName in ['SDF']:
                        self.readTableB_DUMMY()
                    elif tableName in ['MONITOR','PMRF','PERF','PFRF','AEMONPT','FOL','AFRF','AGRF',]:
                        self.readTableB_DUMMY()
                    #elif tableName in []:
                    #    self.readTableB_DUMMY()

                    elif tableName in ['STDISP','FOL','OFMPF2M','OSMPF2M','OPMPF2M','OGPMPF2M','OLMPF2M','DIT','OVGPSD2']:
                        self.readTable_DUMMY_GEOM(tableName)
                    elif tableName in ['OVGATO2','OVGRMS2','OVGNO2']:
                        self.readTable_DUMMY_GEOM(tableName)

                    elif tableName in ['OESNLXR','OESNL1X','OESPSD2','OESNLBR','OESATO2',]:
                        self.readTable_DUMMY_GEOM(tableName)
                    elif tableName in ['OVGCRM2','OAGPSD2','OAGATO2','OAGRMS2','OAGNO2','OAGCRM2','OPGPSD2','OPGPSD2','OPGPSD2','OPGATO2']:
                        self.readTable_DUMMY_GEOM(tableName)
                    elif tableName in ['OPGRMS2','OPGNO2','OPGCRM2','OQGPSD2',]:
                        self.readTable_DUMMY_GEOM(tableName)
                    elif tableName in ['OQGPSD2','OQGATO2','OQGRMS2','OQGNO2','OQGCRM2','PVT0','CASECC','BGPDT','EDOM',]:
                        self.readTable_DUMMY_GEOM(tableName)
                    else:
                        raise InvalidKeywordError('unhandled tableName=|%s|' %(tableName))
                    #print "endTell   = ",self.op2.tell()
                    #print "---isAnotherTable---"
                    (isAnotherTable) = self.hasMoreTables()
                    #isAnotherTable = True
                except EndOfFileError:
                    isAnotherTable = False
                ###
            else:
                (isAnotherTable) = self.skipNextTable()
                continue
            #print self.printSection(140)
            self.log.debug("*** finished tableName = |%r|" %(tableName))
            ###
        ###

        self.log.debug("---end of all tables---")
        self.skippedCardsFile.close()
        
    def parseSortCode(self):
        bits = [0,0,0]
        
        sortCode = self.sortCode
        i=2
        #print "***sortCode = ",self.sortCode
        while sortCode>0:
            value = sortCode%2
            sortCode = (sortCode - value)//2
            bits[i] = value
            #print "    *bit = ",value
            #print "    sortCode = ",sortCode
            i-=1
        #print "sortBits = ",bits
        ## the bytes describe the SORT information
        self.sortBits = bits
        
        self.dataCode['sortBits'] = self.sortBits

    def parseApproachCode(self,data):
        """
        int3 is the 3rd word in table=-3 and may be 
        elementType or something else depending on the table
        """
        #self.printBlock(data)
        (aCode,tCode,int3,iSubcase) = unpack('iiii',data[:16])
        ## the local subcase ID
        self.iSubcase = iSubcase
        self.subcases.add(self.iSubcase) # set notation

        ## the type of result being processed
        self.tableCode = tCode%1000
        ## used to create sortBits
        self.sortCode  = tCode//1000
        ## what type of data was saved from the run; used to parse the approachCode and gridDevice
        self.deviceCode   = aCode%10
        ## what solution was run (e.g. Static/Transient/Modal)
        self.analysisCode = (aCode-self.deviceCode) // 10

        if self.deviceCode==3:
            #sys.stderr.write('The op2 may be inconsistent...\n')
            #sys.stderr.write("  print and plot can cause bad results...if there's a crash, try plot only\n")
            self.deviceCode = 1
            
            self.log.info('The op2 may be inconsistent...')
            self.log.info("  print and plot can cause bad results...if there's a crash, try plot only")
            #pass

        ## dataCode stores the active variables; these pass important
        ## self variables into the result object
        self.dataCode = {'analysisCode': self.analysisCode,
                         'deviceCode'  : self.deviceCode,
                         'tableCode'   : self.tableCode,
                         'sortCode'    : self.sortCode,
                         'dt'          : None,
                         'log'         : self.log,
                         }
        #print "iSubcase = ",self.iSubcase
        self.parseSortCode()

        #print "aCode(1)=%s analysisCode=%s deviceCode=%s tCode(2)=%s tableCode=%s sortCode=%s iSubcase(4)=%s" %(aCode,self.analysisCode,self.deviceCode,tCode,self.tableCode,self.sortCode,self.iSubcase)
        self.log.debug(self.printTableCode(self.tableCode))
        return (int3)

    def getValues(self,data,sFormat,iWordStart,iWordStop=None):
        """
        extracts the ith word from the data structure as the provided type
        supports multiple inputs with iWordStop (note this is words, not outputs)
        @param self the object pointer
        @param data the binary data that is as long as the buffer size
        @param iWordStart the word to start reading from
        @param iWordStop  the word to stop reading on (largely unused)
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
        """
        deletes any parameters before going to the next table to avoid
        messing up data
        """
        params += ['dataCode','deviceCode','analysisCode','tableCode','sortCode','iSubcase',
                   'data','numWide','nonlinearFactor','obj']
        for param in params:
            if hasattr(self,param):
                #print '%s = %s' %(param,getattr(self,param))
                delattr(self,param)

    def getBufferWords(self):
        bufferWords = self.getMarker()
        #print "buffMarker = |%s|" %(bufferWords)
        #print "bufferWords = ",bufferWords,bufferWords*4
        assert bufferWords >0
        return bufferWords

    def verifyBufferSize(self,bufferWords):
        assert bufferWords>0,self.printSection(220)

    def readTitle(self):
        """reads the Title, Subtitle, and Label.
        Puts them in self.iSubcaseNameMap[iSubcase] = [Subtitle,Label]"""
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

        ## the title of the analysis
        self.Title = Title.strip()
        if self.iSubcase not in self.iSubcaseNameMap:
            self.iSubcaseNameMap[self.iSubcase] = [Subtitle,Label]

    def tableInit(self,word):
        """
        starts a new table
        """
        ## the local table name
        self.tableName = word.strip()
        self.tableNames.append(word)
        msg = '*'*20+word+'*'*20+'\n'
        if self.makeOp2Debug:
            self.op2Debug.write(msg)

    def getTableNamesFromOP2(self):
        """
        returns the list of parsed tables
        """
        return self.tableNames

