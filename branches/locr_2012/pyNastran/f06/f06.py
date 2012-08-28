import os
import sys
from itertools import izip

from pyNastran.general.general import printBadPath

#ComplexEigenvalues,strainEnergyDensity,TemperatureGradientObject
from pyNastran.op2.tables.oug.oug_eigenvectors import EigenVectorObject  # ,ComplexEigenVectorObject
from pyNastran.op2.tables.lama_eigenvalues.lama_objects import RealEigenvalues, ComplexEigenvalues

from pyNastran.f06.tables.oes import OES  # OES
from pyNastran.f06.tables.oug import OUG  # OUG
from pyNastran.f06.tables.oqg import OQG  # OUG
from pyNastran.f06.f06_classes import MaxDisplacement  # classes not in op2
from pyNastran.f06.f06Writer import F06Writer



class F06(OES, OUG, OQG, F06Writer):
    def __init__(self, f06FileName, debug=False, log=None):
        """
        Initializes the F06 object
        @param f06FileName
          the file to be parsed
        @param makeGeom
          reads the BDF tables (default=False)
        @param debug
          prints data about how the F06 was parsed (default=False)
        @param log
          a logging object to write debug messages to (@see import logging)
        """
        self.f06FileName = f06FileName
        if not os.path.exists(self.f06FileName):
            msg = 'cant find F06FileName=|%s|\n%s' % (
                self.f06FileName, printBadPath(self.f06FileName))
            raise RuntimeError(msg)
        self.infile = open(self.f06FileName, 'r')
        self.__initAlt__(debug, log)

        self.lineMarkerMap = {
            'R E A L   E I G E N V E C T O R   N O': self.getRealEigenvectors,
        }
        self.markerMap = {
            #'N A S T R A N    F I L E    A N D    S Y S T E M    P A R A M E T E R    E C H O':self.fileSystem,
            #'N A S T R A N    E X E C U T I V E    C O N T R O L    E C H O':self.executiveControl,
            #'C A S E    C O N T R O L    E C H O ':self.caseControl,
            #'M O D E L   S U M M A R Y':self.summary,

            #'E L E M E N T   G E O M E T R Y   T E S T   R E S U L T S   S U M M A R Y'
            'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R': self.getGridWeight,
            #'OLOAD    RESULTANT':self.oload,
            #'MAXIMUM  SPCFORCES':self.getMaxSpcForces,
            #'MAXIMUM  DISPLACEMENTS': self.getMaxDisplacements,
            #'MAXIMUM  APPLIED LOADS': self.getMaxAppliedLoads,
            #'G R I D   P O I N T   S I N G U L A R I T Y   T A B L E': self.gridPointSingularities,


            #------------------------
            #    N O N - D I M E N S I O N A L   S T A B I L I T Y   A N D   C O N T R O L   D E R I V A T I V E   C O E F F I C I E N T S
            #          N O N - D I M E N S I O N A L    H I N G E    M O M E N T    D E R I V A T I V E   C O E F F I C I E N T S
            #                               A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S
            #                              S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S
            #------------------------
            #                                    R O T O R   D Y N A M I C S   S U M M A R Y
            #                              R O T O R   D Y N A M I C S   M A S S   S U M M A R Y
            #                           E I G E N V A L U E  A N A L Y S I S   S U M M A R Y   (COMPLEX LANCZOS METHOD)
            #------------------------

            'R E A L   E I G E N V A L U E S': self.getRealEigenvalues,
            #'C O M P L E X   E I G E N V A L U E   S U M M A R Y':self.getComplexEigenvalues,
            'E L E M E N T   S T R A I N   E N E R G I E S': self.getElementStrainEnergies,
            'D I S P L A C E M E N T   V E C T O R': self.getDisplacement,
            'C O M P L E X   D I S P L A C E M E N T   V E C T O R': self.getComplexDisplacement,
            'F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T': self.getSpcForces,
            'F O R C E S   O F   M U L T I P O I N T   C O N S T R A I N T': self.getMpcForces,
            #'G R I D   P O I N T   F O R C E   B A L A N C E':self.getGridPointForces,

            'S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )': self.getBarStress,
            'S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )': self.getRodStress,

            'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )': self.getTriStress,
            'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )': self.getQuadStress,
            'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN': self.getQuadStressBilinear,
            'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )': self.getQuadCompositeStress,

            'S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )': self.getSolidStressTetra,
            'S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )': self.getSolidStressHexa,
            'S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )': self.getSolidStressPenta,

            'S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )': self.getBarStrain,
            'S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )': self.getRodStrain,

            'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )': self.getQuadStrains,
            'S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )': self.getTriStrain,

            'S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )': self.getSolidStrainTetra,
            'S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )': self.getSolidStrainHexa,
            'S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )': self.getSolidStrainPenta,

            'T E M P E R A T U R E   V E C T O R': self.getTemperatureVector,
            'F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S': self.getTempGradientsFluxes,

            #'* * * END OF JOB * * *': self.end(),
        }
        self.markers = self.markerMap.keys()

    def __initAlt__(self, debug=False, log=None):
        self.i = 0
        self.storedLines = []

        self.displacementsPSD = {}        # random
        self.displacementsATO = {}        # random
        self.displacementsRMS = {}        # random
        self.displacementsCRM = {}        # random
        self.displacementsNO = {}        # random
        self.scaledDisplacements = {}     # tCode=1 thermal=8

        #-----------------------------------
        self.velocities = {}
        self.accelerations = {}
        self.eigenvalues = {}
        self.eigenvectors = {}

        self.thermalLoadVectors = {}
        self.forceVectors = {}
        self.barForces = {}

        self.beamForces = {}
        self.springForces = {}
        self.damperForces = {}
        self.solidPressureForces = {}

        self.rodStrain = {}
        self.nonlinearRodStress = {}
        self.barStrain = {}
        self.plateStrain = {}
        self.nonlinearPlateStrain = {}
        self.compositePlateStrain = {}
        self.solidStrain = {}
        self.beamStrain = {}
        self.ctriaxStrain = {}
        self.hyperelasticPlateStress = {}

        self.rodStress = {}
        self.nonlinearRodStrain = {}
        self.barStress = {}
        self.plateStress = {}
        self.nonlinearPlateStress = {}
        self.compositePlateStress = {}
        self.solidStress = {}
        self.beamStress = {}
        self.ctriaxStress = {}
        self.hyperelasticPlateStrain = {}

        # beam, shear...not done
        self.shearStrain = {}
        self.shearStress = {}

        self.gridPointStresses = {}
        self.gridPointVolumeStresses = {}
        self.gridPointForces = {}
        #-----------------------------------

        self.iSubcases = []
        self.temperatureGrad = {}

        self.iSubcaseNameMap = {}
        self.loadVectors = {}
        self.gridPointForces = {}

        OES.__init__(self)
        OQG.__init__(self)
        OUG.__init__(self)

        ## the TITLE in the Case Control Deck
        self.Title = ''
        self.startLog(log, debug)

    def startLog(self, log=None, debug=False):
        """
        Sets up a dummy logger if one is not provided
        @param self the object pointer
        @param log a python logging object
        @param debug adds debug messages (True/False)
        """
        if log is None:
            from pyNastran.general.logger import dummyLogger
            if debug:
                word = 'debug'
            else:
                word = 'info'
            loggerObj = dummyLogger()
            log = loggerObj.startLog(word)  # or info
        self.log = log

    def getGridPointSingularities(self):  # @todo not done
        """
        @code
                    G R I D   P O I N T   S I N G U L A R I T Y   T A B L E
        POINT    TYPE   FAILED      STIFFNESS       OLD USET           NEW USET
         ID            DIRECTION      RATIO     EXCLUSIVE  UNION   EXCLUSIVE  UNION
          1        G      4         0.00E+00          B        F         SB       S    *
          1        G      5         0.00E+00          B        F         SB       S    *
        @endcode
        """
        pass

    def getMaxSpcForces(self):  # @todo not done
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int, float, float, float, float, float, float])
        #print "max SPC Forces   ",data
        #self.disp[iSubcase] = DisplacementObject(iSubcase,data)
        #print self.disp[iSubcase]

    def getMaxDisplacements(self):  # @todo not done
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int, float, float, float, float, float, float])
        #print "max Displacements",data
        disp = MaxDisplacement(data)
        #print disp.writeF06()
        #self.disp[iSubcase] = DisplacementObject(iSubcase,data)
        #print self.disp[iSubcase]

    def getMaxAppliedLoads(self):  # @todo not done
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int, float, float, float, float, float, float])
        #print "max Applied Loads",data
        #self.disp[iSubcase] = DisplacementObject(iSubcase,data)
        #print self.disp[iSubcase]

    def getGridWeight(self):  # @todo not done
        line = ''
        while 'PAGE' not in line:
            line = self.infile.readline()[1:]
            self.i += 1
        ###

    def readSubcaseNameID(self):
        subcaseName = self.storedLines[-3].strip()
        #print ''.join(self.storedLines)
        self.Title = subcaseName  # 'no title'
        #print "Title = ",self.Title

        #print "subcaseLine = |%r|" %(subcaseName)
        if subcaseName == '':
            iSubcase = 1
        else:
            iSubcase = self.storedLines[-2].strip()[1:]
            if iSubcase == '':  # no subcase specified
                iSubcase = 1
            else:
                iSubcase = int(iSubcase.strip('SUBCASE '))
            ###
            #assert isinstance(iSubcase,int),'iSubcase=|%r|' %(iSubcase)
            #print "subcaseName=%s iSubcase=%s" %(subcaseName,iSubcase)
        ###
        self.iSubcaseNameMap[iSubcase] = [subcaseName,
                                          'SUBCASE %s' % (iSubcase)]
        transient = self.storedLines[-1].strip()
        isSort1 = False
        if transient:
            transWord, transValue = transient.split('=')
            transWord = transWord.strip()
            transValue = float(transValue)
            transient = [transWord, transValue]

            if transWord == 'LOAD STEP':  # nonlinear statics
                analysisCode = 10
            elif transWord == 'TIME STEP':  ## @todo check name
                analysisCode = 6
            elif transWord == 'EIGENVALUE':  # normal modes
                analysisCode = 2
            elif transWord == 'FREQ':  ## @todo check name
                analysisCode = 5
            elif transWord == 'POINT-ID':
                isSort1 = True
                analysisCode = None
            else:
                raise NotImplementedError('transientWord=|%r| is not supported...' % (transWord))
            ###
        else:
            transient = None
            analysisCode = 1

        dt = None
        if transient is not None:
            dt = transient[1]
        return (subcaseName, iSubcase, transient, dt, analysisCode, isSort1)

    def getRealEigenvalues(self):
        """
        @code
                                                   R E A L   E I G E N V A L U E S
         MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED
          NO.       ORDER                                                                       MASS              STIFFNESS
              1         1        6.158494E+07        7.847607E+03        1.248985E+03        1.000000E+00        6.158494E+07
        @endcode
        """
        (subcaseName, iSubcase, transient, dt, analysisCode,
            isSort1) = self.readSubcaseNameID()

        headers = self.skip(2)
        data = self.readTable([int, int, float, float, float, float, float])

        if iSubcase in self.eigenvalues:
            self.eigenvalues[iSubcase].addF06Data(data)
        else:
            self.eigenvalues[iSubcase] = RealEigenvalues(iSubcase)
            self.eigenvalues[iSubcase].addF06Data(data)
        self.iSubcases.append(iSubcase)

    def getComplexEigenvalues(self):
        """
        @code
                               C O M P L E X   E I G E N V A L U E   S U M M A R Y
        ROOT     EXTRACTION                  EIGENVALUE                     FREQUENCY              DAMPING
         NO.        ORDER             (REAL)           (IMAG)                (CYCLES)            COEFFICIENT
             1           6          0.0              6.324555E+01          1.006584E+01          0.0
             2           5          0.0              6.324555E+01          1.006584E+01          0.0
        @endcode
        """
        #(subcaseName,iSubcase,transient,dt,analysisCode,isSort1) = self.readSubcaseNameID()
        iSubcase = 1  # @todo fix this...

        headers = self.skip(2)
        data = self.readTable([int, int, float, float, float, float])

        if iSubcase in self.eigenvalues:
            self.eigenvalues[iSubcase].addF06Data(data)
        else:
            isSort1 = True
            self.eigenvalues[iSubcase] = ComplexEigenvalues(iSubcase)
            self.eigenvalues[iSubcase].addF06Data(data)
        self.iSubcases.append(iSubcase)

    def getRealEigenvectors(self, marker):
        """
        @code
                                                                                                               SUBCASE 1
        EIGENVALUE =  6.158494E+07
            CYCLES =  1.248985E+03         R E A L   E I G E N V E C T O R   N O .          1

        POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
               1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
            2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0

        analysisCode = 2 (Normal modes)
        tableCode    = 7 (Eigenvector)

        deviceCode   = 1 (Print)
        sortCode     = 0 (Sort2,Real,Sorted Results) => sortBits = [0,0,0]
        formatCode   = 1 (Real)
        #sCode        = 0 (Stress)
        numWide      = 8 (???)
        @endcode
        """
        cycle, iMode = marker.strip(
        ).split('R E A L   E I G E N V E C T O R   N O .')
        iMode = int(iMode)

        cycles = cycle.strip().split('=')
        #print smarker
        assert 'CYCLES' == cycles[0].strip(), 'marker=%s' % (marker)
        cycle = float(cycles[1])

        #print "marker = |%s|" %(marker)
        #subcaseName = '???'
        #print self.storedLines
        #iSubcase = self.storedLines[-2].strip()[1:]
        #iSubcase = int(iSubcase.strip('SUBCASE '))
        #print "subcaseName=%s iSubcase=%s" %(subcaseName,iSubcase)
        (subcaseName, iSubcase, transient, dt, analysisCode,
            isSort1) = self.readSubcaseNameID()
        headers = self.skip(2)

        dataCode = {'log': self.log, 'analysisCode': analysisCode,
                    'deviceCode': 1, 'tableCode': 7, 'sortCode': 0,
                    'sortBits': [0, 0, 0], 'numWide': 8, 'formatCode': 1,
                    'mode': iMode, 'eigr': transient[1], 'modeCycle': cycle,
                    'dataNames': ['mode', 'eigr', 'modeCycle'],
                    'name': 'mode', 'tableName': 'OUGV1',
                    'nonlinearFactor': iMode,
                    #'sCode':0,
                    #'elementName':'CBAR','elementType':34,'stressBits':stressBits,
                    }

        dataTypes = [int, str, float, float, float, float, float, float]
        data = self.readTable(dataTypes)

        if iSubcase in self.eigenvectors:
            self.eigenvectors[iSubcase].readF06Data(dataCode, data)
        else:
            isSort1 = True
            self.eigenvectors[iSubcase] = EigenVectorObject(dataCode, isSort1,
                                                            iSubcase, iMode)
            self.eigenvectors[iSubcase].readF06Data(dataCode, data)
        ###

    def getElementStrainEnergies(self):
        """
        @code
        EIGENVALUE = -3.741384E-04
        CYCLES =  3.078479E-03
                                           E L E M E N T   S T R A I N   E N E R G I E S

                ELEMENT-TYPE = QUAD4               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =  -1.188367E-05
                   MODE               1            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =  -1.188367E-05

                                    ELEMENT-ID          STRAIN-ENERGY           PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY
                                             1         -5.410134E-08                -0.0929             -4.328107E-05
                                             2         -3.301516E-09                -0.0057             -2.641213E-06
        @endcode
        """
        iSubcase = 1 ## @todo not correct
        cycles = self.storedLines[-1][1:].strip()
        cycles = float(cycles.split('=')[1])

        eigenvalue = self.storedLines[-2][1:].strip()
        eigenvalue = float(eigenvalue.split('=')[1])
        #print "eigenvalue=%s cycle=%s" %(eigenvalue,cycles)

        eTypeLine = self.skip(2)[1:]
        eType = eTypeLine[30:40]
        totalEnergy1 = eTypeLine[99:114]

        modeLine = self.skip(1)[1:]
        iMode = modeLine[24:40]
        totalEnergy2 = modeLine[99:114]
        #print "eType=%s totalEnergy1=|%s|" %(eType,totalEnergy1)
        #print "iMode=%s totalEnergy2=|%s|" %(iMode,totalEnergy2)
        headers = self.skip(2)

        data = []
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            self.i += 1
            if 'PAGE' in line:
                break
            sline = line.strip().split()
            if sline == []:
                break
            #print sline
            eid = int(sline[0])
            strainEnergy = float(sline[1])
            percentTotal = float(sline[2])
            strainEnergyDensity = float(sline[3])
            out = (eid, strainEnergy, percentTotal, strainEnergyDensity)
            data.append(out)
        ###

        if sline == []:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            self.i += 1
            #print line

        return
        if iSubcase in self.iSubcases:
            self.strainEnergyDensity[iSubcase].readF06Data(data, transient)
        else:
            sed = strainEnergyDensity(data, transient)
            sed.readF06Data(data, transient)
            self.strainEnergyDensity[iSubcase] = sed
        ###

    def getTempGradientsFluxes(self):
        (subcaseName, iSubcase, transient, dt, analysisCode,
            isSort1) = self.readSubcaseNameID()
        #print transient
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readGradientFluxesTable()
        #print data
        return
        if iSubcase in self.temperatureGrad:
            self.temperatureGrad[iSubcase].addData(data)
        else:
            self.temperatureGrad[iSubcase] = TemperatureGradientObject(
                iSubcase, data)
        self.iSubcases.append(iSubcase)

    def readGradientFluxesTable(self):
        data = []
        Format = [int, str, float, float, float, float, float, float]
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            self.i += 1
            if 'PAGE' in line:
                return data
            sline = [line[0:15], line[15:24].strip(), line[24:44], line[44:61], line[61:78], line[78:95], line[95:112], line[112:129]]
            sline = self.parseLineGradientsFluxes(sline, Format)
            data.append(sline)
        return data

    def parseLineGradientsFluxes(self, sline, Format):
        out = []
        for entry, iFormat in izip(sline, Format):
            if entry.strip() is '':
                out.append(0.0)
            else:
                #print "sline=|%r|\n entry=|%r| format=%r" %(sline,entry,iFormat)
                entry2 = iFormat(entry)
                out.append(entry2)
            ###
        return out

    def readTable(self, Format):
        """
        reads displacement, spc/mpc forces
        @param self the object pointer
        @param Format @see parseLine
        """
        sline = True
        data = []
        while sline:
            sline = self.infile.readline()[1:].strip().split()
            self.i += 1
            if 'PAGE' in sline:
                return data
            sline = self.parseLine(sline, Format)
            if sline is None:
                return data
            data.append(sline)
        return data

    def parseLine(self, sline, Format):
        """
        @param self the object pointer
        @param sline list of strings (split line)
        @param Format list of types [int,str,float,float,float] that maps to sline
        """
        out = []
        for entry, iFormat in izip(sline, Format):
            try:
                entry2 = iFormat(entry)
            except:
                #print "sline=|%s|\n entry=|%s| format=%s" %(sline,entry,iFormat)
                return None
            out.append(entry2)
        return out

    def parseLineBlanks(self, sline, Format):
        """allows blanks"""
        out = []
        for entry, iFormat in izip(sline, Format):
            if entry.strip():
                entry2 = iFormat(entry)
            else:
                entry2 = None
                #print "sline=|%s|\n entry=|%s| format=%s" %(sline,entry,iFormat)
            out.append(entry2)
        return out

    def readF06(self):
        """
        Reads the F06 file
        @param self the object pointer
        """
        #print "reading..."
        blank = 0
        while 1:
            #if self.i%1000==0:
                #print "i=%i" %(self.i)
            line = self.infile.readline()
            marker = line[1:].strip()

            #print "marker = %s" %(marker)
            if marker in self.markers:
                blank = 0
                #print "\n*marker = %s" %(marker)
                self.markerMap[marker]()
                self.storedLines = []
                #print "i=%i" %(self.i)
            elif 'R E A L   E I G E N V E C T O R   N O' in marker:
                blank = 0
                #print "\n*marker = %s" %(marker)
                self.lineMarkerMap['R E A L   E I G E N V E C T O R   N O'](
                    marker)
                self.storedLines = []
            elif marker == '':
                blank += 1
                if blank == 20:
                    break
            elif self.isMarker(marker):  # marker with space in it (e.g. Model Summary)
                print("***marker = |%s|" % (marker))

            else:
                blank = 0
            ###
            self.storedLines.append(line)
            self.i += 1
        #print "i=%i" %(self.i)
        self.infile.close()
        self.processF06()

    def processF06(self):
        #data = [self.disp,self.SpcForces,self.stress,self.isoStress,self.barStress,self.solidStress,self.temperature]
        dataPack = [self.solidStress]
        for dataSet in dataPack:
            for key, data in dataSet.iteritems():
                data.processF06Data()

    def isMarker(self, marker):
        """returns True if the word follows the 'N A S T R A N   P A T T E R N'"""
        marker = marker.strip().split('$')[0].strip()

        if len(marker) < 2 or marker == '* * * * * * * * * * * * * * * * * * * *':
            return False
        for i, char in enumerate(marker):
            #print "i=%s i%%2=%s char=%s" %(i,i%2,char)
            if i % 2 == 1 and ' ' is not char:
                return False
            elif i % 2 == 0 and ' ' == char:
                return False
            ###
        ###
        return True

    def skip(self, iskip):
        for i in xrange(iskip - 1):
            self.infile.readline()
        self.i += iskip
        return self.infile.readline()

    def printResults(self):
        msg = ''
        data = [self.displacements, self.spcForces, self.mpcForces, self.temperatures,
                self.eigenvalues, self.eigenvectors,
                self.rodStress, self.rodStrain,
                self.barStress, self.barStrain,
                self.plateStress, self.plateStrain,
                self.compositePlateStress, self.compositePlateStrain,
                ]

        self.iSubcases = list(set(self.iSubcases))
        for iSubcase in self.iSubcases:
            for result in data:
                if iSubcase in result:
                    msg += str(result[iSubcase])
                ###
            ###
        return msg

if __name__ == '__main__':
    f06Name = sys.argv[1]
    model = f06Name.split('.f06')[0]
    f06 = F06(f06Name)
    f06.readF06()

    f06.writeF06(model + 'f06.out')
    f06.printResults()
    print("done...")
