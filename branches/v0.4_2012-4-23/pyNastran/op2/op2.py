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
import os
import sys
from numpy import array
from struct import unpack

from fortranFile import FortranFile
from op2Codes import Op2Codes
from op2Errors import *

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.tables.resultTable import ResultTable
from pyNastran.op2.tables.geom.geometryTables import GeometryTables
from pyNastran.f06.f06Writer import F06Writer

class OP2(BDF,  # BDF methods
          FortranFile,Op2Codes,GeometryTables,ResultTable,F06Writer):

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

    def setTransientTimes(self,times): ## @todo this name sucks...
        """
        takes a dictionary of list of times in a transient case and
        gets the output closest to those timse
        times = {subcaseID_1: [time1, time2],
                 subcaseID_2: [time3, time4]}
        """
        expectedTimes = {}
        for iSubcase,eTimes in times.items():
            eTimes = list(times)
            eTimes.sort()
            expectedTimes[iSubcase] = array(eTimes)
        ###
        self.expectedTimes = expectedTimes

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
        
        ## BDF Title
        self.Title = ''
        ## limit output DTs
        self.expectedTimes = {}
        #self.expectedTimes = {1:array([0.1,0.12])}
        
        ## file object containing the skipped cards
        self.skippedCardsFile = open('skippedCards.out','a')

        ## the list of supported tables (dont edit this)
        self.tablesToRead = ['GEOM1','GEOM2','GEOM3','GEOM4', # nodes/geometry/loads/BCs
                             'GEOM1S','GEOM2S','GEOM3S','GEOM4S', # nodes/geometry/loads/BCs - superelements
                             'GEOM1N', #???
                             'EPT', 'MPT',  # properties/materials
                             'EPTS','MPTS', # properties/materials - superelements
                             'EDTS',         # ???
                             'DYNAMIC','DYNAMICS',
                             'DIT',  # tables
                             'LAMA',

                            'BGPDT','EQEXIN','EQEXINS','PVT0','CASECC','EDOM',
                             'DESTAB',                      # design variables
                             'OQG1','OQGV1','OQMG1',        # spc/mpc forces
                             
                             'OUGV1',                       # displacements                             
                             'OGPFB1','OGS1',               # grid point forces/stresses
                             
                             'OEF1X', 'DOEF1','OEFIT',      # applied forces
                             'OES1X','OES1X1','OES1C',      # stress
                             'OSTR1C','OSTR1X',             # strains
                             'OESNLXR','OESNLXD','OESNL1X','OESNLBR', # nonlinear stress

                             'OESCP',                       # cylinder stress???
                             'OESTRCP',                     # cylinder strain???
                             'OESRT',                       # rotational stress?

                             'ONRGY1', # energy
                             'ONRGY2', # energy (sort2, unsupported)

                             'R1TABRG','HISADD',  # SOL 200

                              # unsupported frequency results
                              'OAGPSD2', 'OAGATO2', 'OAGRMS2', 'OAGNO2', 'OAGCRM2',
                              'OEFPSD2', 'OEFATO2', 'OEFRMS2', 'OEFNO2', 'OEFCRM2',
                              'OESPSD2', 'OESATO2', 'OESRMS2', 'OESNO2', 'OESCRM2',
                              'OPGPSD2', 'OPGATO2', 'OPGRMS2', 'OPGNO2', 'OPGCRM2',
                              'OQGPSD2', 'OQGATO2', 'OQGRMS2', 'OQGNO2', 'OQGCRM2',
                             'OSTRPSD2','OSTRATO2','OSTRRMS2','OSTRNO2','OSTRCRM2',
                              'OUGPSD2', 'OUGATO2', 'OUGCRM2', 'OUGNO2', 'OUGRMS2',
                              'OVGPSD2', 'OVGATO2', 'OVGRMS2', 'OVGNO2', 'OVGCRM2',
                             
                             ## @todo what do these do???
                             'OPG1','OPGV1',
                             'OPNL1',
                             'OUPV1',
                             'VIEWTB','ERRORN',
                             'OFMPF2M','OSMPF2M','OPMPF2M','OGPMPF2M','OLMPF2M',
                             'PCOMPTS',
                             'OMM2',
                             'EDOM',
                             'STDISP','SDF','MONITOR','AEMONPT',
                             'OGPWG', # grid point weight

                             # new
                             'AFRF',             'AGRF',
                             'PMRF','PERF','PFRF',
                             'FOL',
                             ]
        
        ## a dictionary that maps an integer of the subcaseName to the subcaseID
        self.iSubcaseNameMap = {}
        
        ## list of OP2 tables that were read
        ## mainly for debugging
        self.tableNames = []

        ## ESE
        self.eigenvalues = {}

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
        self.forces = {}
        self.fluxes = {}
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
        
        ## OES - CTRIAX6
        self.ctriaxStress = {}
        self.ctriaxStrain = {}

        ## OES - isotropic CROD/CONROD/CTUBE
        self.rodStress  = {}
        ## OES - isotropic CROD/CONROD/CTUBE
        self.rodStrain  = {}
        self.nonlinearRodStress = {}
        self.nonlinearRodStrain = {}
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
        ## OESNLXR - CTRIA3/CQUAD4
        self.nonlinearPlateStress = {}
        self.nonlinearPlateStrain = {}
        self.hyperelasticPlateStress = {}
        self.hyperelasticPlateStrain = {}

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
        

        # OQG - spc/mpc forces
        self.spcForces      = {}
        #self.modalSPCForces = {}
        self.mpcForces      = {}
        #self.modalMPCForces = {}
        
        ## OGF - grid point forces
        self.gridPointForces = {}

        ## OPG - summation of loads for each element
        self.appliedLoads = {}
        self.loadVectors  = {}
        
        ## OEE - strain energy density
        self.strainEnergy = {}

    def readTapeCode(self):
        """
        reads the OP2 header
        @todo whats in this table?
        """
        #self.printSection(500)
        #sys.exit('op2-readTapeCode')

        if 0:
            marker = 0
            #print self.printSection(200)
            sys.exit('stopping in readTapeCode in op2.py')
            while marker != -1:
                ints = self.readIntBlock()
                marker = ints[0]
                #print "ints1 = ",ints
            #print ""

            while marker != -2:
                ints = self.readIntBlock()
                marker = ints[0]
                #print "ints2 = ",ints
            #print ""
            while marker != -3:
                ints = self.readIntBlock()
                marker = ints[0]
                #print "ints3 = ",ints
            #print ""
            while marker != -4:
                ints = self.readIntBlock()
                marker = ints[0]
                #print "ints4 = ",ints
            #print ""
            #while marker != -1:
                #ints = self.readIntBlock()
                #marker = ints[0]
                #print "ints1 = ",ints
            #print ""
            self.readMarkers([2])

            #print self.printSection(200)
            ints2 = self.readIntBlock()
            ints3 = self.readIntBlock()
            #print "ints2 = ",ints2
            #print "ints3 = ",ints3
            sys.exit('stopping...')

        self.readMarkers([3])
        self.printSection(20)
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

        #self.readTapeCode()
        try:
            self.readTapeCode()
        except:
            msg  = 'When this happens, the analysis failed or the code bombed...check the F06.\n'
            msg += '  If the F06 is OK:\n'
            msg += '      1.  Make sure you used PARAM,POST,-1 in your BDF/DAT\n'
            msg += '      2.  Run the problem on a different Operating System'
            msg += '      3.  Are you running an OP2? :)  fname=%s' %(self.op2FileName)
            raise TapeCodeError(msg)
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
            #print "tableName = |%r|" %(tableName)
            if tableName==None:
                break
            elif tableName in self.tablesToRead:
                self.tableName = tableName
                try:
                    #print "startTell = ",self.op2.tell()
                    if tableName=='GEOM1': # nodes,coords,etc.
                        self.readTable_Geom1()
                    elif tableName=='GEOM1S': # superelements - nodes,coords,etc.
                        self.readTable_Geom1S()
                    elif tableName=='GEOM2S': # superelements - elements
                        self.readTable_Geom2S()
                    elif tableName=='GEOM3S': # superelements - static/thermal loads
                        self.readTable_Geom3S()
                    elif tableName=='GEOM4S': # superelements - constraints
                        self.readTable_Geom4S()

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
                    elif  tableName in ['DIT']:  # tables...TABLED1/TABLEM1/TABLES1/GUST
                        self.readTable_DIT()
                    elif tableName in ['LAMA']: # eigenvalue
                        self.readTable_LAMA()

                    elif tableName in ['VIEWTB','EQEXIN','EQEXINS','OEFIT','GEOM1N','OGPWG',]:
                        self.readTable_DUMMY_GEOM(tableName)
                    elif tableName in ['OMM2']:
                        self.readTable_OMM2()
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

                    elif tableName in ['OPG1','OPNL1','OGS1','OPGV1']: # table of applied loads
                        self.readTable_OPG()
                    elif tableName in ['OGPFB1',]:
                        self.readTable_OGF()


                    elif tableName in ['OEF1X','DOEF1',  'OEFPSD2','OEFATO2','OEFRMS2','OEFNO2','OEFCRM2',]:  # applied loads
                        self.readTable_OEF()
                    elif tableName in ['OQG1','OQMG1','OQGV1']:  # spc/mpc forces
                        self.readTable_OQG()

                    elif tableName in ['OUGV1','OUPV1']: # displacements/velocity/acceleration
                        self.readTable_OUG()
                    elif tableName in ['OUGPSD2','OUGATO2','OUGRMS2','OUGNO2','OUGCRM2']: # OUG tables???
                        self.readTable_OUG()

                    elif tableName in ['OES1X','OES1X1','OSTR1X','OES1C','OESCP','OESRT','OESNLXR','OESNL1X']: # stress
                        self.readTable_OES()  # 
                    elif tableName in ['OSTR1X','OSTR1C',]: # strain
                        self.readTable_OES()
                    elif tableName in ['OESTRCP','OESNLXD','OESNLXR',]: # ??? stress/strain
                        self.readTable_OES()
                    elif tableName in ['OSTRATO2','OSTRPSD2','OESRMS2','OESNO2','OESCRM2','OSTRRMS2','OESRMS2','OSTRNO2','OESCRM2','OSTRCRM2',]: # unhandled
                        self.readTable_OES()

                    #elif tableName in ['OESNLXD',]: # dont use this, testing only
                        #self.readTable_OES() # NLXD
                    #elif tableName in ['OESNLXR',]: # dont use this
                        #self.readTable_OES()  # NLXR
                    
                    elif tableName in ['ONRGY1']: # energy
                        self.readTable_OEE()
                    elif tableName in ['ONRGY2']:
                        self.readTable_OEE()

                    elif tableName in ['PCOMPTS']:
                        self.readTable_PCOMPTS()
                    elif tableName in ['SDF']: # ???
                        self.readTable_SDF()
                    #elif tableName in ['CASECC']:
                        #self.readTable_CASECC()

                    # not done
                    elif tableName in []:
                        self.readTableB_DUMMY()
                    elif tableName in ['MONITOR','PMRF','PERF','PFRF','AEMONPT','FOL','AFRF','AGRF',]:
                        self.readTableB_DUMMY()
                    #elif tableName in []:
                    #    self.readTableB_DUMMY()

                    elif tableName in ['STDISP','FOL','OFMPF2M','OSMPF2M','OPMPF2M','OGPMPF2M','OLMPF2M','OVGPSD2']:
                        self.readTable_DUMMY_GEOM(tableName)
                    elif tableName in ['OVGATO2','OVGRMS2','OVGNO2']:
                        self.readTable_DUMMY_GEOM(tableName)

                    elif tableName in ['OESNLXR','OESNL1X','OESPSD2','OESNLBR','OESATO2',]:
                        self.readTable_DUMMY_GEOM(tableName)
                    elif tableName in ['OVGCRM2','OAGPSD2','OAGATO2','OAGRMS2','OAGNO2','OAGCRM2','OPGPSD2','OPGPSD2','OPGPSD2','OPGATO2']:
                        self.readTable_DUMMY_GEOM(tableName)
                    elif tableName in ['OPGRMS2','OPGNO2','OPGCRM2','OQGPSD2',]:
                        self.readTable_DUMMY_GEOM(tableName)
                    elif tableName in ['OQGPSD2','OQGATO2','OQGRMS2','OQGNO2','OQGCRM2','PVT0','CASECC','EDOM',]:
                        self.readTable_DUMMY_GEOM(tableName)
                    elif tableName in ['BGPDT','BGPDTS','EDTS',]:
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
        """
        sortCode = 0 -> sortBits = [0,0,0]
        sortCode = 1 -> sortBits = [0,0,1]
        sortCode = 2 -> sortBits = [0,1,0]
        sortCode = 3 -> sortBits = [0,1,1]
        etc.
        sortCode = 7 -> sortBits = [1,1,1]

        sortBits[0] = 0 -> isSort1=True  isSort2=False
        sortBits[1] = 0 -> isReal=True   isReal/Imaginary=False
        sortBits[2] = 0 -> isSorted=True isRandom=False
        """
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
        elementType or something else depending on the table type
        """
        (aCode,tCode,int3,iSubcase) = unpack('iiii',data[:16])
        ## the local subcase ID
        self.iSubcase = iSubcase
        #print "iSubcase = ",iSubcase
        self.subcases.add(self.iSubcase) # set notation

        ## the type of result being processed
        self.tableCode = tCode%1000
        ## used to create sortBits
        self.sortCode = tCode//1000
        ## what type of data was saved from the run; used to parse the approachCode and gridDevice
        ## deviceCode defines what options inside a result, STRESS(PLOT,PRINT), are used.
        self.deviceCode = aCode%10
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
        if bufferWords <=0:
            raise BufferError('An invalid buffersize was found...bufferWords=%s tableName=%s section=\n%s' %(bufferWords,self.tableName,self.printSection(200)))
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

    def printResults(self):
        results = [
                   # OUG - Displacements/Velocity/Acceleration/Temperature/Heat Flux/
                   #       SPC Forces
                   self.displacements,self.temperatures,
                   self.eigenvalues,
                   self.eigenvectors,
                   self.velocities,
                   self.accelerations,
                   #self.nonlinearTemperatures,self.nonlinearDisplacements,
                   #self.forces,self.fluxes,
                   
                   # OEF - Applied Forces/Temperatures - ???
                   #self.nonlinearForces,self.nonlinearFluxes,
                   #self.temperatureForces,
                   
                   # OQG1 - Forces
                   #self.spcForces,self.mpcForces,
                   
                   # OGF - Grid Point Forces
                   self.gridPointForces,

                   # OPG - Applied Force/Moment
                   self.appliedLoads,
                   self.loadVectors,

                   # OES - Stress/Strain
                   self.celasStress,self.celasStrain,
                   self.rodStress,self.rodStrain,
                   self.nonlinearRodStress,self.nonlinearRodStrain,

                   self.barStress,self.barStrain,
                   self.beamStress,self.beamStrain,
                   self.plateStress,self.plateStrain,
                   self.solidStress,self.solidStrain,
                   self.compositePlateStress,self.compositePlateStrain,
                   self.ctriaxStress,self.ctriaxStrain, # strain not coded???
                   
                   # OEE - Strain Energy
                   self.strainEnergy,
                   ]
        
        msg = '---ALL RESULTS---\n'
        for result in results:
            for iSubcase,res in sorted(result.items()):
                msg += 'iSubcase = %s\n' %(iSubcase)
                try:
                    msg += str(res) + '\n'
                except:
                    print 'failed on %s' %(res.__class__.__name__)
                    raise
                ###
            ###
        ###
        return msg
        
