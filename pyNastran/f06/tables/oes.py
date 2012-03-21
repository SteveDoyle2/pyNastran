from pyNastran.op2.tables.oes.oes_rods   import rodStressObject
from pyNastran.op2.tables.oes.oes_bars   import barStressObject
#from pyNastran.op2.tables.oes.oes_beams   import beamStressObject
#from pyNastran.op2.tables.oes.oes_shear   import shearStressObject
from pyNastran.op2.tables.oes.oes_solids import solidStressObject
from pyNastran.op2.tables.oes.oes_plates import plateStressObject
from pyNastran.op2.tables.oes.oes_compositePlates import compositePlateStressObject

#strain...

class OES(object):
    def __init__(self):
        #self.celasStress   = {} # CELASi
        #self.celasStrain   = {}
        
        self.rodStress  = {} # CROD
        self.rodStrain  = {}

        self.barStress  = {} # CBAR
        self.barStrain  = {}

        #self.beamStress = {} # CBEAM
        #self.beamStrain = {}

        self.plateStress = {} # isotropic CTRIA3/CQUAD4
        self.plateStrain = {}

        self.solidStress = {} # CTETRA/CPENTA/CHEXA
        self.solidStrain = {}

        self.compositePlateStress = {}  # composite CTRIA3/CQUAD4
        self.compositePlateStrain = {}
        
        #self.shearStress = {} # CSHEAR
        #self.shearStrain = {}


    def getRodStress(self):
        """
                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
        ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
          ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
             14    2.514247E+04              1.758725E+02                     15    2.443757E+04              2.924619E+01  
        analysisCode = 1 (Statics)
        deviceCode   = 1 (Print)
        tableCode    = 5 (Stress)
        sortCode     = 0 (Sort2,Real,Sorted Results) => sortBits = [0,0,0]
        formatCode   = 1 (Real)
        sCode        = 0 (Stress)
        numWide      = 8 (???)
        """
        (subcaseName,iSubcase,transient,analysisCode,isSort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        
        stressBits = self.makeStressBits()
        dataCode = {'log':self.log,'analysisCode':analysisCode,'deviceCode':1,'tableCode':5,'sortCode':0,
                    'sortBits':[0,0,0],'numWide':8,'sCode':0,'stressBits':stressBits,
                    'formatCode':1,'elementName':'ROD','elementType':1,
                    }

        data = self.readRodStress()
        if iSubcase in self.rodStress:
            self.rodStress[iSubcase].addF06Data(data,transient)
        else:
            self.rodStress[iSubcase] = rodStressObject(dataCode,iSubcase,transient)
            self.rodStress[iSubcase].addF06Data(data,transient)
        self.iSubcases.append(iSubcase)
    
    def readRodStress(self):
        """
                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
        ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
          ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
             14    2.514247E+04              1.758725E+02                     15    2.443757E+04              2.924619E+01  
        """
        data = []
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            sline = [line[0:13],line[13:29],line[29:42],line[42:55],line[55:67],   line[67:78],line[78:94],line[94:107],line[107:120],line[120:131]]
            if 'PAGE' in line:
                break
            print sline
            out = self.parseLineBlanks(sline,[int,float,float,float,float, int,float,float,float,float]) # line 1
            print out
            data.append(out[:5])
            if isinstance(out[5],int):
                data.append(out[5:])
            self.i+=1
            ###
        ###
        #print "--------"
        #for line in data:
            #print line
        #sys.exit()
        return data

    def getBarStress(self):
        """
                                       S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )
        ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
          ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C
             12    0.0            0.0            0.0            0.0            1.020730E+04   1.020730E+04   1.020730E+04 
                   0.0            0.0            0.0            0.0                           1.020730E+04   1.020730E+04 
        analysisCode = 1 (Statics)
        deviceCode   = 1 (Print)
        tableCode    = 5 (Stress)
        sortCode     = 0 (Sort2,Real,Sorted Results) => sortBits = [0,0,0]
        formatCode   = 1 (Real)
        sCode        = 0 (Stress)
        numWide      = 8 (???)
        """
        (subcaseName,iSubcase,transient,analysisCode,isSort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        
        stressBits = self.makeStressBits()
        dataCode = {'log':self.log,'analysisCode':analysisCode,'deviceCode':1,'tableCode':5,'sortCode':0,
                    'sortBits':[0,0,0],'numWide':8,'sCode':0,'stressBits':stressBits,
                    'formatCode':1,'elementName':'CBAR','elementType':34,
                    }

        data = self.readBarStress()
        if iSubcase in self.barStress:
            self.barStress[iSubcase].addF06Data(data,transient)
        else:
            self.barStress[iSubcase] = barStressObject(dataCode,iSubcase,transient)
            self.barStress[iSubcase].addF06Data(data,transient)
        self.iSubcases.append(iSubcase)
    
    def readBarStress(self):
        """
        ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
          ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C
             12    0.0            0.0            0.0            0.0            1.020730E+04   1.020730E+04   1.020730E+04 
                   0.0            0.0            0.0            0.0                           1.020730E+04   1.020730E+04 
        """
        data = []
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            sline = [line[0:12],line[12:27],line[27:42],line[42:57],line[57:64],line[64:86],line[86:101],line[101:116],line[116:131]]
            if 'PAGE' in line:
                break
            #print sline
            out = self.parseLineBlanks(sline,[int,float,float,float,float, float, float,float,float]) # line 1
            out = ['CBAR']+out
            #data.append(sline)
            line = self.infile.readline()[1:].rstrip('\r\n ')
            sline = [line[12:27],line[27:42],line[42:57],line[57:64],line[86:101],line[101:116],line[116:131]]
            #print sline
            out += self.parseLineBlanks(sline,[    float,float,float,float,        float,float,float]) # line 2
            #print "*",out
            data.append(out)
            self.i+=2
            ###
        ###
        #print "--------"
        #for line in data:
            #print line
        #sys.exit()
        return data

    def getQuadCompositeStress(self):
        """
                       S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
        ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
          ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
            181    1   3.18013E+04  5.33449E+05  1.01480E+03   -7.06668E+01  1.90232E+04   89.88  5.33451E+05  3.17993E+04  2.50826E+05
            181    2   1.41820E+05  1.40805E+05  1.25412E+05   -1.06000E+02  2.85348E+04   44.88  2.66726E+05  1.58996E+04  1.25413E+05
        
        elementType = 33 b/c not bilinear
        """
        (subcaseName,iSubcase,transient,analysisCode,isSort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int,int,float,float,float,float,float,float,float,float,float])

        isMaxShear = False  # Von Mises/Max Shear
        if 'SHEAR' in headers.strip():
            isMaxShear = True
        stressBits = self.makeStressBits(isMaxShear=isMaxShear)
        dataCode = {'log':self.log,'analysisCode':analysisCode,'deviceCode':1,'tableCode':5,'sortCode':0,
                    'sortBits':[0,0,0],'numWide':8,'sCode':0,'stressBits':stressBits,
                    'formatCode':1,'elementName':'CQUAD4','elementType':33,
                    }

        if iSubcase in self.compositePlateStress:
            self.compositePlateStress[iSubcase].addF06Data(data,transient,'CQUAD4')
        else:
            self.compositePlateStress[iSubcase] = compositePlateStressObject(dataCode,iSubcase,transient)
            self.compositePlateStress[iSubcase].addF06Data(data,transient,'CQUAD4')
        self.iSubcases.append(iSubcase)
        
    def getTriStress(self):
        """
                                 S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )
        ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 
          ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
              8   -1.250000E-01     -1.303003E+02   1.042750E+04  -1.456123E+02   -89.2100    1.042951E+04   -1.323082E+02   1.049629E+04
                   1.250000E-01     -5.049646E+02   1.005266E+04  -2.132942E+02   -88.8431    1.005697E+04   -5.092719E+02   1.032103E+04
        analysisCode = 1 (Statics)
        deviceCode   = 1 (Print)
        tableCode    = 5 (Stress)
        sortCode     = 0 (Sort2,Real,Sorted Results) => sortBits = [0,0,0]
        formatCode   = 1 (Real)
        sCode        = 0 (Stress)
        numWide      = 8 (???)
        """
        (subcaseName,iSubcase,transient,analysisCode,isSort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        
        isFiberDistance = False
        isMaxShear = False  # Von Mises/Max Shear
        if 'DISTANCE' in headers:
            isFiberDistance = True
        if 'MAX SHEAR' in headers:
            isMaxShear = True
        stressBits = self.makeStressBits(isFiberDistance,isMaxShear)
        dataCode = {'log':self.log,'analysisCode':analysisCode,'deviceCode':1,'tableCode':5,'sortCode':0,
                    'sortBits':[0,0,0],'numWide':8,'sCode':1,'stressBits':stressBits,
                    'formatCode':1,'elementName':'CTRIA3','elementType':74,
                    }

        data = self.readTriStress()
        if iSubcase in self.plateStress:
            self.plateStress[iSubcase].addF06Data(data,transient)
        else:
            self.plateStress[iSubcase] = plateStressObject(dataCode,iSubcase,transient)
            self.plateStress[iSubcase].addF06Data(data,transient)
        self.iSubcases.append(iSubcase)

    def readTriStress(self):
        """
                ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 
                  ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
                      8   -1.250000E-01     -1.303003E+02   1.042750E+04  -1.456123E+02   -89.2100    1.042951E+04   -1.323082E+02   1.049629E+04
                           1.250000E-01     -5.049646E+02   1.005266E+04  -2.132942E+02   -88.8431    1.005697E+04   -5.092719E+02   1.032103E+04
        """
        data = []
        while 1:
            line = self.infile.readline()[1:].strip().split()
            if 'PAGE' in line:
                break
            #print line
            sline = self.parseLine(line,[int,float, float,float,float, float,float,float, float]) # line 1
            sline = ['CTRIA3']+sline
            data.append(sline)
            line = self.infile.readline()[1:].strip().split()
            sline += self.parseLine(line,[    float, float,float,float, float,float,float, float]) # line 2
            data.append(sline)
            self.i+=2
            ###
        ###
        return data

    def getQuadStress(self):
        """
                             S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN

        ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)
          ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       VON MISES
              6    CEN/4  -1.250000E-01  -4.278394E+02  8.021165E+03 -1.550089E+02   -88.9493   8.024007E+03 -4.306823E+02  8.247786E+03
                           1.250000E-01   5.406062E+02  1.201854E+04 -4.174177E+01   -89.7916   1.201869E+04  5.404544E+02  1.175778E+04

                       4  -1.250000E-01  -8.871141E+02  7.576036E+03 -1.550089E+02   -88.9511   7.578874E+03 -8.899523E+02  8.060780E+03
                           1.250000E-01  -8.924081E+01  1.187899E+04 -4.174177E+01   -89.8002   1.187913E+04 -8.938638E+01  1.192408E+04
        """
        (subcaseName,iSubcase,transient,analysisCode,isSort1) = self.readSubcaseNameID()
        headers = self.skip(3)
        #print "headers = %s" %(headers)
        
        isFiberDistance = False
        isMaxShear = False  # Von Mises/Max Shear
        if 'DISTANCE' in headers:
            isFiberDistance = True
        if 'MAX SHEAR' in headers:
            isMaxShear = True
        stressBits = self.makeStressBits(isFiberDistance,isMaxShear)
        dataCode = {'log':self.log,'analysisCode':analysisCode,'deviceCode':1,'tableCode':5,'sortCode':0,
                    'sortBits':[0,0,0],'numWide':8,'sCode':1,'stressBits':stressBits,
                    'formatCode':1,'elementName':'CQUAD4','elementType':144,
                    }

        data = self.readQuadBilinear()
        if iSubcase in self.plateStress:
            self.plateStress[iSubcase].addF06Data(data,transient)
        else:
            self.plateStress[iSubcase] = plateStressObject(dataCode,iSubcase,transient)
            self.plateStress[iSubcase].addF06Data(data,transient)
        self.iSubcases.append(iSubcase)

    def readQuadBilinear(self):
        data = []
        while 1:
            if 1: # CEN/4
                line = self.infile.readline()[1:].strip().split()
                if 'PAGE' in line:
                    return data
                sline = self.parseLine(line,[int,str,float, float,float,float, float,float,float, float]) # line 1
                sline = ['CQUAD4']+sline
                #data.append(sline)
                line = self.infile.readline()[1:].strip().split()
                sline += self.parseLine(line,[        float, float,float,float, float,float,float, float]) # line 2
                data.append(sline)
                line = self.infile.readline() # blank line
                self.i+=3
            ###
            for i in range(4):
                line = self.infile.readline()[1:].strip().split()
                sline = self.parseLine(line,[int,float, float,float,float, float,float,float, float]) # line 1
                #data.append(sline)
                line = self.infile.readline()[1:].strip().split()
                sline += self.parseLine(line,[    float, float,float,float, float,float,float, float]) # line 2
                data.append(sline)
                line = self.infile.readline() # blank line
                self.i+=3
            ###
        ###
        return data

    def getSolidStressHexa(self):
        return self.readSolidStress('CHEXA',8)
    def getSolidStressPenta(self):
        return self.readSolidStress('CPENTA',6)
    def getSolidStressTetra(self):
        return self.readSolidStress('CTETRA',4)

    def readSolidStress(self,eType,n):
        """
        analysisCode = 1 (Statics)
        deviceCode   = 1 (Print)
        tableCode    = 5 (Stress/Strain)
        sortCode     = 0 (Sort2,Real,Sorted Results) => sortBits = [0,0,0]
        formatCode   = 1 (Real)
        sCode        = 0 (Stress)
        numWide      = 8 (???)
        """
        (subcaseName,iSubcase,transient,analysisCode,isSort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        isMaxShear = True
        if 'VON MISES' in headers:
            isMaxShear = False
            
        data = self.read3DStress(eType,n)
        stressBits = self.makeStressBits(isMaxShear=False)
        dataCode = {'log':self.log,'analysisCode':1,'deviceCode':1,'tableCode':5,
                    'sortCode':0,'sortBits':[0,0,0],'numWide':8,'elementName':eType,'formatCode':1,
                    'sCode':0,'stressBits':stressBits}

        if iSubcase in self.solidStress:
            self.solidStress[iSubcase].addF06Data(data,transient)
        else:
            self.solidStress[iSubcase] = solidStressObject(dataCode,iSubcase,transient)
            self.solidStress[iSubcase].addF06Data(data,transient)
        self.iSubcases.append(iSubcase)

    def read3DStress(self,eType,n):
        data = []
        while 1:
            line = self.infile.readline().rstrip('\n\r') #[1:]
                    #              CENTER         X          #          XY             #        A         #
            sline = [line[1:17],line[17:24],line[24:28],line[28:43],line[43:47],line[47:63],line[63:66],line[66:80],  line[80:83],line[83:88],line[88:93],line[93:98],line[99:113],line[113:130]]
            sline = [s.strip() for s in sline]
            if 'PAGE' in line:
                break
            elif '' is not sline[0]:
                sline = [eType]+sline
            data.append(sline)
        ###
        return data

    def makeStressBits(self,isFiberDistance=False,isMaxShear=True):
        print "isMaxShear=%s isFiberDistance=%s" %(isMaxShear,isFiberDistance)
        stressBits = [0,0,0]
        if isMaxShear==False:
            stressBits[0] = 1 # Von Mises
        
        if isFiberDistance:
            stressBits[2] = 1
        return stressBits
        