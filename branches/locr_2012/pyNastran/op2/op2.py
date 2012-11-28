#pylint: disable=C0103,W0201,W0223,R0901,R0902,R0904
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
from numpy import array
from struct import unpack

from pyNastran.op2.fortranFile import FortranFile
from pyNastran.op2.op2Codes import Op2Codes


from pyNastran.op2.tables.resultTable import ResultTable
from pyNastran.op2.tables.geom.geometryTables import GeometryTables

from pyNastran.bdf.bdf import BDF
from pyNastran.f06.f06Writer import F06Writer
from pyNastran.f06.matlabWriter import MatlabWriter


class OP2(BDF,
          FortranFile, Op2Codes, GeometryTables, ResultTable, F06Writer,
          MatlabWriter):

    def setSubcases(self, iSubcases=None):
        """
        Allows you to read only the subcases in the list of iSubcases
        @param iSubcases
          list of [subcase1_ID,subcase2_ID]  (default=None; all subcases)
        """
        ## stores the set of all subcases that are in the OP2
        self.subcases = set()
        if iSubcases is None or iSubcases == []:
            ## stores if the user entered [] for iSubcases
            self.isAllSubcases = True
            self.validSubcases = []
        else:
            ## should all the subcases be read (default=True)
            self.isAllSubcases = False
            ## the set of valid subcases -> set([1,2,3])
            self.validSubcases = set(iSubcases)
        self.log.debug("setSubcases - iSubcases = %s" % (self.validSubcases))

    def setTransientTimes(self, times):  ## @todo this name sucks...
        """
        takes a dictionary of list of times in a transient case and
        gets the output closest to those timse
        @code
        times = {subcaseID_1: [time1, time2],
                 subcaseID_2: [time3, time4]}
        @endcode
        """
        expectedTimes = {}
        for (iSubcase, eTimes) in times.iteritems():
            eTimes = list(times)
            eTimes.sort()
            expectedTimes[iSubcase] = array(eTimes)
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

    def __init__(self, op2FileName, makeGeom=False, debug=True, log=None):
        """
        Initializes the Op2 object
        @param op2FileName
          the file to be parsed
        @param makeGeom
          reads the BDF tables (default=False)
        @param debug
          prints data about how the OP2 was parsed (default=False)
        @param log
          a logging object to write debug messages to (@see import logging)
        """
        BDF.__init__(self, debug=debug, log=log)
        self.setSubcases()  # initializes the variables
        self.log.debug('op2FileName = %s' % (op2FileName))
        bdfExtension = '.bdf'
        f06Extension = '.f06'
        (fname, extension) = os.path.splitext(op2FileName)
        self.tableName = 'temp'

        ## should the BDF tables be parsed
        self.makeGeom = makeGeom

        ## the input OP2 filename
        self.op2FileName = op2FileName

        ## the expected BDF filename (guessed)
        self.bdfFileName = fname + bdfExtension

        ## the expected F06 filename (guessed)
        self.f06FileName = fname + f06Extension
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
        self.skippedCardsFile = open('skippedCards.out', 'a')

        ## the list of supported tables (dont edit this)
        self.tablesToRead = [
            # nodes/geometry/loads/BCs
            'GEOM1', 'GEOM2', 'GEOM3', 'GEOM4',      # regular
            'GEOM1S', 'GEOM2S', 'GEOM3S', 'GEOM4S',  # superelements

            'GEOM1OLD', 'GEOM1N',  # ???
            'EPT', 'MPT',    # properties/materials
            'EPTS', 'MPTS',  # properties/materials - superelements
            'EDTS',          # ???
            'DYNAMIC', 'DYNAMICS',
            'DIT',                           # tables (e.g. TABLED1)
            'LAMA', 'BLAMA',                 # eigenvalues

            'BGPDT', 'BGPDTS',               # boundary grids???
            'EQEXIN', 'EQEXINS', 'PVT0', 'CASECC', 'EDOM',
            'DESTAB',                        # design variables
            'OQG1', 'OQGV1',                 # spc forces
            'OQMG1',                         # mpc forces

            'OUGV1', 'OUG1',                 # displacements
            'OGPFB1',                        # grid point forces
            'OGS1',                          # grid point stresses

            'OEF1X', 'DOEF1', 'OEFIT',        # element forces
            'OPG1', 'OPGV1', 'OPNL1',          # applied forces
            'OES1', 'OES1X', 'OES1X1', 'OES1C',  # stress
            'OSTR1C', 'OSTR1X',               # strains
            'OESNLXR', 'OESNLXD', 'OESNL1X', 'OESNLBR',  # nonlinear stress

            'OESCP',                       # cylinder stress???
            'OESTRCP',                     # cylinder strain???
            'OESRT',                       # rotational stress?

            'ONRGY1',  # energy
            'ONRGY2',  # energy (sort2, unsupported)

            'R1TABRG', 'HISADD',  # SOL 200

            # unsupported frequency results
            'OAGPSD2', 'OAGATO2', 'OAGRMS2', 'OAGNO2', 'OAGCRM2',
            'OEFPSD2', 'OEFATO2', 'OEFRMS2', 'OEFNO2', 'OEFCRM2',
            'OESPSD2', 'OESATO2', 'OESRMS2', 'OESNO2', 'OESCRM2',
            'OPGPSD2', 'OPGATO2', 'OPGRMS2', 'OPGNO2', 'OPGCRM2',
            'OQGPSD2', 'OQGATO2', 'OQGRMS2', 'OQGNO2', 'OQGCRM2',

            # supported-ish
            'OQMPSD2', 'OQMATO2', 'OQMRMS2', 'OQMNO2', 'OQMCRM2',
            'OSTRPSD2', 'OSTRATO2', 'OSTRRMS2', 'OSTRNO2', 'OSTRCRM2',
            'OUGPSD2', 'OUGATO2', 'OUGCRM2', 'OUGNO2', 'OUGRMS2',
            'OVGPSD2', 'OVGATO2', 'OVGRMS2', 'OVGNO2', 'OVGCRM2',

            ## @todo what do these do???
            'OUPV1',
            'VIEWTB', 'ERRORN',
            'OFMPF2M', 'OSMPF2M', 'OPMPF2M', 'OGPMPF2M', 'OLMPF2M',
            'PCOMPTS',
            'OMM2',
            'EDOM',
            'STDISP', 'SDF', 'MONITOR', 'AEMONPT',
            'OGPWG',  # grid point weight
            'OQP1',
            'OCRUG', 'OCRPG',

            # new
            'AFRF', 'AGRF',
            'PMRF', 'PERF', 'PFRF',
            'FOL',
        ]

        ## a dictionary that maps an integer of the subcaseName to the
        ## subcaseID
        self.iSubcaseNameMap = {}

        ## list of OP2 tables that were read
        ## mainly for debugging
        self.tableNames = []

        self.__objectsInit__()

    def __objectsInit__(self):
        ## ESE
        self.eigenvalues = {}

        ## OUG - displacement
        self.displacements = {}           # tCode=1 thermal=0
        self.displacementsPSD = {}        # random
        self.displacementsATO = {}        # random
        self.displacementsRMS = {}        # random
        self.displacementsCRM = {}        # random
        self.displacementsNO = {}        # random
        self.scaledDisplacements = {}     # tCode=1 thermal=8

        ## OUG - temperatures
        self.temperatures = {}           # tCode=1 thermal=1

        ## OUG - eigenvectors
        self.eigenvectors = {}            # tCode=7 thermal=0

        ## OUG - velocity
        self.velocities = {}              # tCode=10 thermal=0

        ## OUG - acceleration
        self.accelerations = {}           # tCode=11 thermal=0

        # OEF - Forces - tCode=4 thermal=0
        self.rodForces = {}
        self.barForces = {}
        self.bar100Forces = {}
        self.beamForces = {}
        self.bendForces = {}
        self.bushForces = {}
        self.coneAxForces = {}
        self.damperForces = {}
        self.gapForces = {}
        self.plateForces = {}
        self.plateForces2 = {}
        self.shearForces = {}
        self.solidPressureForces = {}
        self.springForces = {}
        self.viscForces = {}

        self.force_VU = {}
        self.force_VU_2D = {}

        #OEF - Fluxes - tCode=4 thermal=1
        self.thermalLoad_CONV = {}
        self.thermalLoad_CHBDY = {}
        self.thermalLoad_1D = {}
        self.thermalLoad_2D_3D = {}
        self.thermalLoad_VU = {}
        self.thermalLoad_VU_3D = {}
        self.thermalLoad_VUBeam = {}
        #self.temperatureForces = {}

        # OES - tCode=5 thermal=0 sCode=0,1 (stress/strain)
        ## OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
        self.celasStress = {}
        ## OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
        self.celasStrain = {}

        ## OES - CTRIAX6
        self.ctriaxStress = {}
        self.ctriaxStrain = {}

        ## OES - isotropic CROD/CONROD/CTUBE stress
        self.rodStress = {}
        ## OES - isotropic CROD/CONROD/CTUBE strain
        self.rodStrain = {}
        ## OES - nonlinear CROD/CONROD/CTUBE stress
        self.nonlinearRodStress = {}
        ## OES - nonlinear CROD/CONROD/CTUBE strain
        self.nonlinearRodStrain = {}
        ## OES - isotropic CBAR stress
        self.barStress = {}
        ## OES - isotropic CBAR strain
        self.barStrain = {}
        ## OES - isotropic CBEAM stress
        self.beamStress = {}
        ## OES - isotropic CBEAM strain
        self.beamStrain = {}
        ## OES - isotropic CBUSH stress
        self.bushStress = {}
        ## OES - isotropic CBUSH strain
        self.bushStrain = {}
         ## OES - isotropic CBUSH1D strain/strain
        self.bush1dStressStrain = {}

        ## OES - isotropic CTRIA3/CQUAD4 stress
        self.plateStress = {}
        ## OES - isotropic CTRIA3/CQUAD4 strain
        self.plateStrain = {}
        ## OESNLXR - CTRIA3/CQUAD4 stress
        self.nonlinearPlateStress = {}
        ## OESNLXR - CTRIA3/CQUAD4 strain
        self.nonlinearPlateStrain = {}
        self.hyperelasticPlateStress = {}
        self.hyperelasticPlateStrain = {}

        ## OES - isotropic CTETRA/CHEXA/CPENTA stress
        self.solidStress = {}
        ## OES - isotropic CTETRA/CHEXA/CPENTA strain
        self.solidStrain = {}
        ## OES - composite CTRIA3/CQUAD4 stress
        self.compositePlateStress = {}
        ## OES - composite CTRIA3/CQUAD4 strain
        self.compositePlateStrain = {}

        ## OES - CSHEAR stress
        self.shearStress = {}
        ## OES - CSHEAR strain
        self.shearStrain = {}

        ## OES - CELAS1 224, CELAS3 225, 
        self.nonlinearSpringStress = {}
        ## OES - GAPNL 86
        self.nonlinearGapStress = {}
        ## OES - CBUSH 226
        self.nolinearBushStress = {}

        # OQG - spc/mpc forces
        self.spcForces = {}  # tCode=3?
        self.mpcForces = {}  # tCode=39
        
        # OQG - thermal forces
        self.thermalGradientAndFlux = {}

        ## OGF - grid point forces
        self.gridPointForces = {}  # tCode=19

        ## OGS1 - grid point stresses
        self.gridPointStresses = {}       # tCode=26
        self.gridPointVolumeStresses = {}  # tCode=27

        ## OPG - summation of loads for each element
        self.loadVectors = {}       # tCode=2  thermal=0
        self.thermalLoadVectors = {}  # tCode=2  thermal=1
        self.appliedLoads = {}       # tCode=19 thermal=0
        self.forceVectors = {}       # tCode=12 thermal=0

        ## OEE - strain energy density
        self.strainEnergy = {}  # tCode=18

    def get_op2_stats(self):
        """
        gets info about the contents of the different attributes of the
        OP2 class
        """
        table_types = [
            ## OUG - displacement
            'displacements',
            'displacementsPSD',
            'displacementsATO',
            'displacementsRMS',
            'displacementsCRM',
            'displacementsNO',
            'scaledDisplacements',

            ## OUG - temperatures
            'temperatures',

            ## OUG - eigenvectors
            'eigenvectors',

            ## OUG - velocity
            'velocities',

            ## OUG - acceleration
            'accelerations',

            # OQG - spc/mpc forces
            'spcForces',
            'mpcForces',
            'thermalGradientAndFlux',

            ## OGF - grid point forces
            'gridPointForces',

            ## OPG - summation of loads for each element
            'loadVectors',
            'thermalLoadVectors',
            'appliedLoads',
            'forceVectors',

            # OES - tCode=5 thermal=0 sCode=0,1 (stress/strain)
            ## OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
            'celasStress',
            ## OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
            'celasStrain',

            ## OES - isotropic CROD/CONROD/CTUBE stress
            'rodStress',
            ## OES - isotropic CROD/CONROD/CTUBE strain
            'rodStrain',
            ## OES - isotropic CBAR stress
            'barStress',
            ## OES - isotropic CBAR strain
            'barStrain',
            ## OES - isotropic CBEAM stress
            'beamStress',
            ## OES - isotropic CBEAM strain
            'beamStrain',

            ## OES - isotropic CTRIA3/CQUAD4 stress
            'plateStress',
            # OES - isotropic CTRIA3/CQUAD4 strain
            'plateStrain',
            ## OES - isotropic CTETRA/CHEXA/CPENTA stress
            'solidStress',
            ## OES - isotropic CTETRA/CHEXA/CPENTA strain
            'solidStrain',

            ## OES - CSHEAR stress
            'shearStress',
            ## OES - CSHEAR strain
            'shearStrain',
            ## OES - CEALS1 224, CELAS3 225
            'nonlinearSpringStress',
            ## OES - GAPNL 86
            'nonlinearGapStress',
            ## OES - CBUSH 226
            'nolinearBushStress',

        ]

        table_types += [
            # LAMA
            'eigenvalues',

            # OEF - Forces - tCode=4 thermal=0
            'rodForces',
            'barForces',
            'bar100Forces',
            'beamForces',
            'bendForces',
            'bushForces',
            'coneAxForces',
            'damperForces',
            'gapForces',
            'plateForces',
            'plateForces2',
            'shearForces',
            'solidPressureForces',
            'springForces',
            'viscForces',

            'force_VU',
            'force_VU_2D',

            #OEF - Fluxes - tCode=4 thermal=1
            'thermalLoad_CONV',
            'thermalLoad_CHBDY',
            'thermalLoad_1D',
            'thermalLoad_2D_3D',
            'thermalLoad_VU',
            'thermalLoad_VU_3D',
            'thermalLoad_VUBeam',
            #self.temperatureForces
        ]
        table_types += [
            ## OES - CTRIAX6
            'ctriaxStress',
            'ctriaxStrain',
            
            'bushStress',
            'bushStrain',
            'bush1dStressStrain',

            ## OES - nonlinear CROD/CONROD/CTUBE stress
            'nonlinearRodStress',
            'nonlinearRodStrain',

            ## OESNLXR - CTRIA3/CQUAD4 stress
            'nonlinearPlateStress',
            'nonlinearPlateStrain',
            'hyperelasticPlateStress',
            'hyperelasticPlateStrain',

            ## OES - composite CTRIA3/CQUAD4 stress
            'compositePlateStress',
            'compositePlateStrain',

            ## OGS1 - grid point stresses
            'gridPointStresses',        # tCode=26
            'gridPointVolumeStresses',  # tCode=27

            ## OEE - strain energy density
            'strainEnergy',  # tCode=18
        ]
        msg = []
        for table_type in table_types:
            table = getattr(self, table_type)
            for isubcase, subcase in sorted(table.iteritems()):
                if hasattr(subcase,'get_stats'):
                    msg.append('op2.%s[%s]\n' % (table_type, isubcase))
                    msg.extend(subcase.get_stats())
                    msg.append('\n')
                else:
                    msg.append('skipping %s op2.%s[%s]\n\n' % (subcase.__class__.__name__, table_type, isubcase))
                    #raise RuntimeError('skipping %s op2.%s[%s]\n\n' % (subcase.__class__.__name__, table_type, isubcase))

        return ''.join(msg)

    def readTapeCode(self):
        """
        Reads the OP2 header.  This table is still very much in development.
        @todo whats in this table?
        """
        #self.printSection(500)
        #sys.exit('op2-readTapeCode')

        self.readMarkers([3])
        #print(self.printSection(20))
        ints = self.readIntBlock()
        if self.makeOp2Debug:
            self.op2Debug.write('%s\n' % (str(ints)))
        #print "*ints = ",ints
        self.readMarkers([7])

        word = self.readStringBlock()  # Nastran Fort Tape ID Code -
        #print "word = |%r|" %(word)

        self.readMarkers([2])
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-1])

        #data = self.getData(60)
        #self.printBlock(data)

    def readTapeCodePost2(self):
        # PVTO
        #self.readMarkers([2])
        #ints = self.readIntBlock()
        #print("*ints = ",ints)
        
        #ints = self.readBlock()
        #print("*ints = ",ints)
        #print(self.printSection(40))

        block = self.readNewBlock()
        #print("len(block) = %s" %(len(block)))
        if len(block)==12:  # post = -1
            form = -1
            (month,day,year) = unpack(b'iii',block)
            print("date = %s/%s/%s" %(month,day,2000+year))

            block = self.readNewBlock()
            #print("len(block) = %s" %(len(block)))
            print(self.printBlock(block))
            
            block = self.readNewBlock()
            #print("len(block) = %s" %(len(block)))
            print(self.printBlock(block))

            self.readMarkers([-1,0])

            block = self.readNewBlock()
            #print("len(block) = %s" %(len(block)))
            tableName, = unpack(b'8s',block)
            print("tableName=%s self.n=%s" %(tableName,self.n))

            marker = self.readMarker()
            print('marker = %s' %(marker))
            #self.readMarkers([-1])
            block = self.readNewBlock()
            print(self.printBlock(block))

            self.readMarkers([-2])
            block = self.readNewBlock()
            print(self.printBlock(block))
            
            marker = self.readMarker()
            print('marker = %s' %(marker))

            #self.readMarkers([-1])
            #self.readPartOfNewTable()  # PVTO

            #self.readNewTable()
            #self.readNewTable()
            print(self.printSection(100))
            
            sys.exit('stoppingD')

        elif len(block)==8: # post = -2
            self.n = 0
            self.op2.seek(self.n)
            tableName, = unpack(b'8s',block)
            print("tableName = |%s|\n" %(tableName))
            #sys.exit('stoppingA')

            block = self.readNewBlock()
            #print(self.printBlock(block))
            self.readMarkers([-1])
            #print("")

            self.readPartOfNewTable()  # PVTO

            block = self.readNewBlock()
            print(self.printBlock(block))
            self.readPartOfNewTable()  # PRTMAXIM=YES AUTOSPC=NO
            
            #block = self.readNewBlock()
            #print(self.printSection(200))
            #return
            
            self.readNewTable()
            self.readNewTable()
            self.readNewTable()
            #print(self.printSection(40))
            sys.exit('stoppingB')

            return

        
        sys.exit()
        
        # -----------------------
        
        foundMoreTables = True
        while foundMoreTables and tableName:
            print("tableName=%s self.n=%s" %(tableName,self.n))
            print("-----------------------------")
            try:
                n = self.n
                tableName = self.readNewTable()
            except:
                self.op2.seek(n)
                foundMoreTables = False

        print("**************")
        self.op2.seek(n)
        print(self.printSection(240))
        sys.exit('stoppingC')
    
    def readNewTable(self):
        self.op2.read(8); self.n+=8
        data = self.op2.read(8); self.n+=12

        if len(data)<8:
            return

        self.op2.seek(self.n)
        tableName, = unpack(b'8s',data)
        print("tableName = |%s|\n" %(tableName))

        self.readMarkers([-1])
        block = self.readNewBlock()
        #print(self.printBlock(block))

        self.readPartOfNewTable()
        return tableName

    def readPartOfNewTable(self):
        n = -2
        keepGoingOnTable = True
        while keepGoingOnTable:
            print("n = %s" %(n))
            try:
                nStar = self.n
                self.readMarkers([n,1,0])
                block = self.readNewBlock()
                #print(self.printBlock(block))
                n -= 1
            except:
                self.n = nStar
                self.op2.seek(nStar)
                keepGoingOnTable = False

    def readNewBlock(self):
        data = self.op2.read(16)
        #print(self.printBlock(data))
        (four,n,four,fourN) = unpack(b'iiii',data)
        #print('n = %s' %(n))
        self.n += 16
        
        if n > 70000:
            asf
        data = self.op2.read(n*4)
        self.n += n*4+4
        #print("self.n = ",self.n)
        self.op2.seek(self.n)
        
        return data
        
    def readOP2(self):
        """
        reads the op2 file
        """
        ## the OP2 file object
        self.op2 = open(self.op2FileName, 'rb')

        if self.makeOp2Debug:
            ## a developer debug file (largely unsupported)
            self.op2Debug = open('debug.out', 'wb')
        ## the byte position in the OP2
        self.n = self.op2.tell()

        try:
            #self.readTapeCodePost2()
            self.readTapeCode()
        except:
            raise
            msg = ('When this happens, the analysis failed or '
                  'the code bombed...check the F06.\n'
                  '  If the F06 is OK:\n'
                  '      1.  Make sure you used PARAM,POST,-1 in your '
                  'BDF/DAT/NAS\n'
                  '      2.  Run the problem on a different Operating System\n'
                  '      3.  Are you running an OP2? :)  \n'
                  'fname=%s' % (self.op2FileName))
            raise RuntimeError("Tape Code Error: %s", msg)

        isAnotherTable = True
        while isAnotherTable:
            self.log.debug('-' * 80)
            try:
                tableName = self.readTableName(
                    rewind=True, stopOnFailure=False)
                print("tableName = %s" % (tableName))
            except EOFError:  # the isAnotherTable method sucks...
                isAnotherTable = False
                self.log.debug("***ok exit, but it could be better...")
                break
            except SyntaxError:  # Invalid Markers, the isAnotherTable method sucks...
                isAnotherTable = False
                self.log.debug("***poor exit, but it worked...")
                #raise
                break
            except:
                raise
            self.log.debug("tableName = |%r|" % (tableName))
            #print("tableName = |%r|" %(tableName))

            if tableName is None:
                break
            else:
                isAnotherTable = self.readTable(tableName)

        self.log.debug("---end of all tables---")
        self.skippedCardsFile.close()

    def readTable(self, tableName):
        if tableName in self.tablesToRead:
            self.tableName = tableName
            self.isRegular = True
            try:
                #print("startTell = %s" %(self.op2.tell()))
                if tableName == 'GEOM1':  # nodes,coords,etc.
                    self.readTable_Geom1()
                elif tableName == 'GEOM1S':  # superelements
                    self.readTable_Geom1S()  #  - nodes,coords,etc.
                elif tableName == 'GEOM2S':  # superelements
                    self.readTable_Geom2S()  #  - elements
                elif tableName == 'GEOM3S':  # superelements
                    self.readTable_Geom3S()  #  - static/thermal loads
                elif tableName == 'GEOM4S':  # superelements
                    self.readTable_Geom4S()  #  - constraints

                #elif tableName=='GEOM1OLD':
                #    self.readTable_Geom1Old()
                #elif tableName=='GEOM1N':
                #    self.readTable_Geom1N()
                elif tableName == 'GEOM2':  # elements
                    self.readTable_Geom2()
                elif tableName == 'GEOM3':  # static/thermal loads
                    self.readTable_Geom3()
                elif tableName == 'GEOM4':  # constraints
                    self.readTable_Geom4()

                elif tableName in ['EPT', 'EPTS']:  # element properties
                    self.readTable_EPT()
                elif tableName in ['MPT', 'MPTS']:  # material properties
                    self.readTable_MPTS()
                elif tableName in ['DYNAMIC', 'DYNAMICS']:  # dyanmic info
                    self.readTable_DYNAMICS()
                elif  tableName in ['DIT']:  # tables...
                    self.readTable_DIT()     # TABLED1/TABLEM1/TABLES1/GUST
                elif tableName in ['LAMA', 'BLAMA']:  # eigenvalue
                    self.readTable_LAMA()

                elif tableName in ['VIEWTB', 'EQEXIN', 'EQEXINS', 'OEFIT',
                                   'GEOM1N', 'OGPWG', 'GEOM1OLD']:
                    self.readTable_DUMMY_GEOM(tableName)
                elif tableName in ['OMM2']:
                    self.readTable_OMM2()
                elif tableName in ['DESTAB']:  # design variable table
                    self.readTable_DesTab()
                elif tableName in ['R1TABRG']:  # not done - response table
                    self.readTable_R1TAB()
                    self.isOptimization = True
                elif tableName in ['HISADD']:  # not done
                    self.readTable_R1TAB()
                    self.isOptimization = True
                elif tableName in ['ERRORN']:  # not done
                    self.readTable_R1TAB()

                elif tableName in ['OPG1', 'OPNL1', 'OGS1', 'OPGV1']:
                    self.readTable_OPG()  # table of applied loads
                elif tableName in ['OGPFB1', ]:
                    self.readTable_OGF()
                elif tableName in ['OCRUG', 'OCRPG']:  # totally guessing...
                    self.readTable_OUG()

                elif tableName in ['OEF1X', 'DOEF1', 'OEFPSD2', 'OEFATO2',
                                   'OEFRMS2', 'OEFNO2', 'OEFCRM2', ]:
                    self.readTable_OEF()  # applied loads
                elif tableName in ['OQG1', 'OQGV1', 'OQP1', ]:  # spc forces
                    self.readTable_OQG()
                elif tableName in ['OQMG1', 'OQMPSD2', 'OQMATO2', 'OQMRMS2',
                                   'OQMNO2', 'OQMCRM2', ]:  # mpc forces
                    #self.readTable_OQG()
                    self.readTable_DUMMY_GEOM(tableName)

                elif tableName in ['OUGV1', 'OUPV1', 'OUG1']:
                    self.readTable_OUG()  # displacements/velocity/acceleration
                elif tableName in ['OUGPSD2', 'OUGATO2', 'OUGRMS2', 'OUGNO2',
                                   'OUGCRM2']:  # OUG tables???
                    self.readTable_OUG2()
                    #self.readTable_OUG()

                elif tableName in ['OES1', 'OES1X', 'OES1X1', 'OSTR1X',
                                   'OES1C', 'OESCP', 'OESRT', 'OESNLXR',
                                   'OESNL1X']:
                    self.readTable_OES()  # stress
                elif tableName in ['OSTR1X', 'OSTR1C', ]:
                    self.readTable_OES()  # strain
                elif tableName in ['OESTRCP', 'OESNLXD', 'OESNLXR', ]:
                    self.readTable_OES()  # ??? stress/strain
                elif tableName in ['OSTRATO2', 'OSTRPSD2', 'OESRMS2', 'OESNO2',
                                   'OESCRM2', 'OSTRRMS2', 'OESRMS2', 'OSTRNO2',
                                   'OESCRM2', 'OSTRCRM2', ]:  # unhandled
                    self.readTable_OES()

                #elif tableName in ['OESNLXD',]: # dont use this, testing only
                    #self.readTable_OES() # NLXD
                #elif tableName in ['OESNLXR',]: # dont use this
                    #self.readTable_OES()  # NLXR

                elif tableName in ['ONRGY1']:  # energy
                    self.readTable_OEE()
                elif tableName in ['ONRGY2']:
                    self.readTable_OEE()

                elif tableName in ['PCOMPTS']:
                    self.readTable_PCOMPTS()
                elif tableName in ['SDF']:  # ???
                    self.readTable_SDF()
                #elif tableName in ['CASECC']:
                    #self.readTable_CASECC()

                # not done
                elif tableName in []:
                    self.readTableB_DUMMY()
                elif tableName in ['MONITOR']:
                    self.readMonitor()
                elif tableName in ['PMRF', 'PERF', 'PFRF',
                                   'AEMONPT', 'FOL', 'AFRF', 'AGRF', ]:
                    self.readTableB_DUMMY()
                #elif tableName in []:
                #    self.readTableB_DUMMY()

                elif tableName in ['STDISP', 'FOL', 'OFMPF2M', 'OSMPF2M',
                                   'OPMPF2M', 'OGPMPF2M', 'OLMPF2M',
                                   'OVGPSD2']:
                    self.readTable_DUMMY_GEOM(tableName)
                elif tableName in ['OVGATO2', 'OVGRMS2', 'OVGNO2']:
                    self.readTable_DUMMY_GEOM(tableName)

                elif tableName in ['OESNLXR', 'OESNL1X', 'OESPSD2', 'OESNLBR',
                                   'OESATO2', ]:
                    self.readTable_DUMMY_GEOM(tableName)
                elif tableName in ['OVGCRM2', 'OAGPSD2', 'OAGATO2', 'OAGRMS2',
                                   'OAGNO2', 'OAGCRM2', 'OPGPSD2', 'OPGPSD2',
                                   'OPGPSD2', 'OPGATO2']:
                    self.readTable_DUMMY_GEOM(tableName)
                elif tableName in ['OPGRMS2', 'OPGNO2', 'OPGCRM2', 'OQGPSD2']:
                    self.readTable_DUMMY_GEOM(tableName)
                elif tableName in ['OQGPSD2', 'OQGATO2', 'OQGRMS2', 'OQGNO2',
                                   'OQGCRM2', 'PVT0', 'CASECC', 'EDOM', ]:
                    self.readTable_DUMMY_GEOM(tableName)
                elif tableName in ['BGPDT', 'BGPDTS', 'EDTS', ]:
                    self.readTable_DUMMY_GEOM(tableName)
                else:
                    msg = 'unhandled tableName=|%s|' % (tableName)
                    raise KeyError(msg)
                #print("endTell   = ",self.op2.tell())
                #print("---isAnotherTable---")
                (isAnotherTable) = self.hasMoreTables()
                #isAnotherTable = True
            except EOFError:
                isAnotherTable = False
        else:
            if tableName not in [None]:
                raise NotImplementedError('%s is not supported' % (tableName))
            (isAnotherTable) = self.skipNextTable()
            #return isAnotherTable
        #print(self.printSection(140))
        self.log.debug("*** finished tableName = |%r|" % (tableName))
        return isAnotherTable

    def readMonitor(self):
        tableName = self.readTableName(rewind=False)  # LAMA
        self.tableInit(tableName)
        print("tableName1 = |%r|" %(tableName))
        print("tableName2 = |%r|" %(self.tableName))

        self.readMarkers([-1, 7], 'MONITOR')
        ints = self.readIntBlock()
        print("*ints = ",ints)

        self.readMarkers([-2, 1, 0], 'MONITOR')
        bufferWords = self.getMarker()
        print("bufferWords = ",bufferWords)

        word = self.readStringBlock()  # MONITOR
        print("word = |%s|" %(word))

        self.readMarkers([-3, 1, 0], 'MONITOR')
        bufferWords = self.getMarker()
        print("bufferWords = ",bufferWords,bufferWords*4)
        
        data = self.op2.read(4)
        data = self.op2.read(bufferWords*4)
        Format = str(bufferWords*4) + 's'
        Format = bytes(Format)
        word, = unpack(Format, data)
        print("word = ",word)
        data = self.op2.read(4)
        self.n += bufferWords*4+8
        

        self.readMarkers([-4, 1, 0], 'MONITOR')
        bufferWords = self.getMarker()
        print("bufferWords = ",bufferWords,bufferWords*4)
        data = self.op2.read(4)
        data = self.op2.read(bufferWords*4)
        Format = str(bufferWords*4) + 's'
        Format = bytes(Format)
        word, = unpack(Format, data)
        print("word = ",word)
        data = self.op2.read(4)
        self.n += bufferWords*4+8

        self.readMarkers([-5, 1, 0], 'MONITOR')

        #word = self.readStringBlock()  # MONITOR
        #print("word = |%s|" %(word))



        #print(self.printSection(200))
        #sys.exit()

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
        bits = [0, 0, 0]

        sortCode = self.sortCode
        i = 2
        #print("***sortCode = ",self.sortCode)
        while sortCode > 0:
            value = sortCode % 2
            sortCode = (sortCode - value) // 2
            bits[i] = value
            #print("    *bit = %s" %(value))
            #print("    sortCode = %s" %(sortCode))
            i -= 1
        #print("sortBits = %s" %(bits))
        ## the bytes describe the SORT information
        self.sortBits = bits

        self.dataCode['sortBits'] = self.sortBits

    def parseApproachCode(self, data):
        """
        int3 is the 3rd word in table=-3 and may be
        elementType or something else depending on the table type
        """
        (aCode, tCode, int3, iSubcase) = unpack(b'iiii', data[:16])
        ## the local subcase ID
        self.iSubcase = iSubcase
        #print("iSubcase = %s" %(iSubcase))
        self.subcases.add(self.iSubcase)  # set notation

        ## the type of result being processed
        self.tableCode = tCode % 1000
        ## used to create sortBits
        self.sortCode = tCode // 1000
        ## what type of data was saved from the run; used to parse the
        ## approachCode and gridDevice.  deviceCode defines what options
        ## inside a result, STRESS(PLOT,PRINT), are used.
        self.deviceCode = aCode % 10
        ## what solution was run (e.g. Static/Transient/Modal)
        self.analysisCode = (aCode - self.deviceCode) // 10

        if self.deviceCode == 3:
            #sys.stderr.write('The op2 may be inconsistent...\n')
            #sys.stderr.write("  print and plot can cause bad results..."
            #                 "if there's a crash, try plot only\n")
            self.deviceCode = 1

            #self.log.info('The op2 may be inconsistent...')
            #self.log.info('  print and plot can cause bad results...'
            #              'if there's a crash, try plot only')

        ## dataCode stores the active variables; these pass important
        ## self variables into the result object
        self.dataCode = {'analysisCode': self.analysisCode,
                         'deviceCode': self.deviceCode,
                         'tableCode': self.tableCode,
                         'sortCode': self.sortCode,
                         'dt': None,
                         'log': self.log,
                         }
        #print("iSubcase = ",self.iSubcase)
        self.parseSortCode()

        #print('aCode(1)=%s analysisCode=%s deviceCode=%s '
        #      'tCode(2)=%s tableCode=%s sortCode=%s iSubcase(4)=%s'
        #      %(aCode, self.analysisCode, self.deviceCode, tCode, 
        #        self.tableCode, self.sortCode, self.iSubcase))
        #self.log.debug(self.printTableCode(self.tableCode))
        return (int3)

    def parseApproachCode2(self, data):
        """
        int3 is the 3rd word in table=-3 and may be
        elementType or something else depending on the table type
        """
        (aCode, tCode, int3, ID) = unpack(b'iiii', data[:16])
        ## the local subcase ID
        self.ID = ID
        #print("iSubcase = %s" %(iSubcase))
        self.subcases.add(self.ID)  # set notation

        ## the type of result being processed
        self.tableCode = tCode % 1000
        ## used to create sortBits
        self.sortCode = tCode // 1000
        ## what type of data was saved from the run; used to parse the
        ## approachCode and gridDevice.  deviceCode defines what options
        ## inside a result, STRESS(PLOT,PRINT), are used.
        self.deviceCode = aCode % 10
        ## what solution was run (e.g. Static/Transient/Modal)
        self.analysisCode = (aCode - self.deviceCode) // 10

        if self.deviceCode == 3:
            #sys.stderr.write('The op2 may be inconsistent...\n')
            #sys.stderr.write("  print and plot can cause bad results..."
            #                 "if there's a crash, try plot only\n")
            self.deviceCode = 1

            #self.log.info('The op2 may be inconsistent...')
            #self.log.info('  print and plot can cause bad results...'
            #              'if there's a crash, try plot only')

        ## dataCode stores the active variables; these pass important
        ## self variables into the result object
        self.dataCode = {'analysisCode': self.analysisCode,
                         'deviceCode': self.deviceCode,
                         'tableCode': self.tableCode,
                         'sortCode': self.sortCode,
                         'dt': None,
                         'log': self.log,
                         }
        #print("iSubcase = ",self.iSubcase)
        self.parseSortCode()

        #print('aCode(1)=%s analysisCode=%s deviceCode=%s '
        #      'tCode(2)=%s tableCode=%s sortCode=%s iSubcase(4)=%s'
        #      %(aCode, self.analysisCode, self.deviceCode, tCode, 
        #        self.tableCode, self.sortCode, self.iSubcase))
        #self.log.debug(self.printTableCode(self.tableCode))
        return (int3)

    def getValues(self, data, iFormat, iWordStart, iWordStop=None):
        """
        extracts the ith word from the data structure as the provided type
        supports multiple inputs with iWordStop (note this is words,
        not outputs)
        
        @param self the object pointer
        @param data the binary data that is as long as the buffer size
        @param iWordStart the word to start reading from
        @param iWordStop  the word to stop reading on (largely unused)
        @warning
            works with nastran syntax, not standard python syntax
            this makes it work with what's documented in the DMAP manual
        """
        if iWordStop is None:
            #print("iWordStart=%s data[%s:%s]" %(iWordStart,iWordStart*4,
            #                                   (iWordStart+1)*4))
            ds = data[(iWordStart - 1) * 4:iWordStart * 4]
            iFormat = bytes(iFormat)
            return unpack(iFormat, ds)[0]
        #print("type(data) = ",type(data))
        ds = data[(iWordStart - 1) * 4:(iWordStop - 1) * 4]
        return unpack(iFormat, ds)

    def deleteAttributes(self, params):
        """
        deletes any parameters before going to the next table to avoid
        messing up data
        """
        params += ['dataCode', 'deviceCode', 'analysisCode', 'tableCode',
                   'sortCode', 'iSubcase', 'data', 'numWide',
                   'nonlinearFactor', 'obj', 'subtitle', 'label']
        for param in params:
            if hasattr(self, param):
                delattr(self, param)

    def getBufferWords(self):
        bufferWords = self.getMarker()
        #print("buffMarker = |%s|" %(bufferWords))
        #print("bufferWords = ",bufferWords,bufferWords*4)
        if bufferWords <= 0:
            raise BufferError('An invalid buffersize was found...'
                              'bufferWords=%s tableName=%s section=\n%s'
                              % (bufferWords, self.tableName,
                                 self.printSection(200)))
        return bufferWords

    def verifyBufferSize(self, bufferWords):
        assert bufferWords > 0, self.printSection(220)

    def readTitle(self):
        """
        reads the Title, Subtitle, and Label.
        Puts them in self.iSubcaseNameMap[iSubcase] = [Subtitle,Label]
        """
        ## the title of the analysis
        word = self.readString(384)  # titleSubtitleLabel
        self.Title = word[0:128].strip()
        ## the subtitle of the subcase
        self.subtitle = word[128:256].strip()
        ## the label of the subcase
        self.label = word[256:328].strip()

        # not really a hollerith, just the end of the block (so bufferWords*4)
        self.readHollerith()

        self.dataCode['subtitle'] = self.subtitle
        self.dataCode['label'] = self.label

        if hasattr(self,'iSubcase'):
            if self.iSubcase not in self.iSubcaseNameMap:
                self.iSubcaseNameMap[self.iSubcase] = [self.subtitle, self.label]
        else:
            pass

    def tableInit(self, word):
        """
        starts a new table
        """
        ## the local table name
        self.tableName = word.strip()
        ## the names of all the tables
        self.tableNames.append(word)
        #msg = '*' * 20 + word + '*' * 20 + '\n'
        #print(msg)

    def getTableNamesFromOP2(self):
        """
        returns the list of parsed tables
        """
        return self.tableNames

    def printResults(self):
        results = [
            # OUG - Displacements/Velocity/Acceleration/Temperature
            self.displacements, self.displacementsPSD,
            self.displacementsATO,
            self.temperatures,
            self.eigenvalues,
            self.eigenvectors,
            self.velocities,
            self.accelerations,

            # OEF - Applied Forces/Temperatures - ???

            # OQG1 - SPC/MPC Forces
            self.spcForces, self.mpcForces,

            # OGF - Grid Point Forces
            self.gridPointForces,

            # OPG - Applied Force/Moment
            self.appliedLoads,
            self.loadVectors, self.thermalLoadVectors,
            self.forceVectors,

            # OES - Stress/Strain
            self.celasStress, self.celasStrain,
            self.rodStress, self.rodStrain,
            self.nonlinearRodStress, self.nonlinearRodStrain,

            self.barStress, self.barStrain,
            self.beamStress, self.beamStrain,
            self.plateStress, self.plateStrain,
            self.solidStress, self.solidStrain,
            self.compositePlateStress, self.compositePlateStrain,
            self.ctriaxStress, self.ctriaxStrain,

            # OEE - Strain Energy
            self.strainEnergy,
        ]

        msg = '---ALL RESULTS---\n'
        for result in results:
            for (iSubcase, res) in sorted(result.iteritems()):
                msg += 'iSubcase = %s\n' % (iSubcase)
                try:
                    msg += str(res) + '\n'
                except:
                    print('failed on %s' % (res.__class__.__name__))
                    raise
        return msg