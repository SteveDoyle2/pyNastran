import os
import sys

from tables.oes import *  # OES
from tables.oug import *  # OUG

class EndOfFileError(Exception):
    pass

class F06Reader(OES,OUG):
    def __init__(self,f06name):
        self.f06name = f06name
        self.i = 0
        self.markerMap = {
          #'N A S T R A N    F I L E    A N D    S Y S T E M    P A R A M E T E R    E C H O':self.fileSystem,
          #'N A S T R A N    E X E C U T I V E    C O N T R O L    E C H O':self.executiveControl,
          #'C A S E    C O N T R O L    E C H O ':self.caseControl,
          #'M O D E L   S U M M A R Y':self.summary,
          
          #'E L E M E N T   G E O M E T R Y   T E S T   R E S U L T S   S U M M A R Y'
          #'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R'
          #'OLOAD    RESULTANT':self.oload,
          'MAXIMUM  SPCFORCES':self.maxSpcForces,
          'MAXIMUM  DISPLACEMENTS': self.maxDisplacements,
          'MAXIMUM  APPLIED LOADS': self.maxAppliedLoads,
          #'G R I D   P O I N T   S I N G U L A R I T Y   T A B L E': self.gridPointSingularities,


#------------------------
#    N O N - D I M E N S I O N A L   S T A B I L I T Y   A N D   C O N T R O L   D E R I V A T I V E   C O E F F I C I E N T S
#          N O N - D I M E N S I O N A L    H I N G E    M O M E N T    D E R I V A T I V E   C O E F F I C I E N T S
#                               A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S
#                              S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S
#------------------------


          'D I S P L A C E M E N T   V E C T O R':self.displacement,
          'F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T':self.spcForces,
          'F O R C E S   O F   M U L T I P O I N T   C O N S T R A I N T': self.mpcForces,
          #'G R I D   P O I N T   F O R C E   B A L A N C E':self.gridPointForces,
          #'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )': self.quadStress,
          #'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )': self.triStress,
          
          'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )': self.quadCompositeStress,


          'S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )':self.barStress,
          'S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )':self.solidStressHexa,
          'S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )':self.solidStressPenta,
          'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN': self.quadStress,
          #'S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )'
          'S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )':self.solidStressTetra,
          'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )': self.triStress,


          'T E M P E R A T U R E   V E C T O R':self.temperatureVector,
          'F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S':self.tempGradientsFluxes,

          #'* * * END OF JOB * * *': self.end(),
         }
        self.markers = self.markerMap.keys()
        self.infile = open(self.f06name,'r')
        self.storedLines = []

        self.disp = {}
        self.SpcForces = {}
        self.iSubcases = []
        self.temperature = {}
        self.temperatureGrad = {}
        OES.__init__(self)
        self.startLog()

    def startLog(self,log=None,debug=False):
        if log is None:
            from pyNastran.general.logger import dummyLogger
            loggerObj = dummyLogger()
            if debug:
                word = 'debug'
            else:
                word = 'info'
            log = loggerObj.startLog(word) # or info
        self.log = log

    def gridPointSingularities(self):
        """
                    G R I D   P O I N T   S I N G U L A R I T Y   T A B L E
        POINT    TYPE   FAILED      STIFFNESS       OLD USET           NEW USET
         ID            DIRECTION      RATIO     EXCLUSIVE  UNION   EXCLUSIVE  UNION
          1        G      4         0.00E+00          B        F         SB       S    *
          1        G      5         0.00E+00          B        F         SB       S    *
        """
        pass

    def maxSpcForces(self):
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int,float,float,float,float,float,float])
        print "max SPC Forces   ",data
        #self.disp[iSubcase] = DisplacementObject(iSubcase,data)
        #print self.disp[iSubcase]

    def maxDisplacements(self):
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int,float,float,float,float,float,float])
        print "max Displacements",data
        #self.disp[iSubcase] = DisplacementObject(iSubcase,data)
        #print self.disp[iSubcase]
    
    def maxAppliedLoads(self):
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int,float,float,float,float,float,float])
        print "max Applied Loads",data
        #self.disp[iSubcase] = DisplacementObject(iSubcase,data)
        #print self.disp[iSubcase]

    def readSubcaseNameID(self):
        subcaseName = self.storedLines[-3].strip()
        iSubcase   = self.storedLines[-2].strip()[1:]
        iSubcase = int(iSubcase.strip('SUBCASE '))
        #print "subcaseName=%s iSubcase=%s" %(subcaseName,iSubcase)

        transient   = self.storedLines[-1].strip()
        if transient:
            transWord,transValue = transient.split('=')
            transWord = transWord.strip()
            transValue = float(transValue)
            transient = [transWord,transValue]
            
            if transWord=='LOAD STEP':
                analysisCode = 10
            elif transWord=='TIME STEP': ## @todo check name
                analysisCode = 6
            elif transWord=='FREQ': ## @todo check name
                analysisCode = 5
            else:
                raise NotImplementedError('transientWord=|%r| is not supported...' %(transWord))
            ###
        else:
            transient = None
            analysisCode = 1

        return (subcaseName,iSubcase,transient,analysisCode)

    def tempGradientsFluxes(self):
        (subcaseName,iSubcase,transient,analysisCode) = self.readSubcaseNameID()
        #print transient
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readGradientFluxesTable()
        #print data
        return
        if iSubcase in self.temperatureGrad:
            self.temperatureGrad[iSubcase].addData(data)
        else:
            self.temperatureGrad[iSubcase] = TemperatureGradientObject(iSubcase,data)
        self.iSubcases.append(iSubcase)
    
    def readGradientFluxesTable(self):
        data = []
        Format = [int,str,float,float,float, float,float,float]
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            if 'PAGE' in line:
                return data
            sline = [line[0:15],line[15:24].strip(),line[24:44],line[44:61],line[61:78],line[78:95],line[95:112],line[112:129]]
            sline = self.parseLineGradientsFluxes(sline,Format)
            data.append(sline)
        return data

    def parseLineGradientsFluxes(self,sline,Format):
        out = []
        for entry,iFormat in zip(sline,Format):
            if entry.strip() is '':
                out.append(0.0)
            else:
                #print "sline=|%r|\n entry=|%r| format=%r" %(sline,entry,iFormat)
                entry2 = iFormat(entry)
                out.append(entry2)
            ###
        return out
    
    def spcForces(self):
        (subcaseName,iSubcase,transient,analysisCode) = self.readSubcaseNameID()
        headers = self.skip(2)
        return
        #print "headers = %s" %(headers)
        data = self.readTable([int,str,float,float,float,float,float,float])

        if iSubcase in self.SpcForces:
            self.SpcForces[iSubcase].addData(data)
        else:
            self.SpcForces[iSubcase] = DisplacementObject(iSubcase,data)
        self.iSubcases.append(iSubcase)
        #print self.SpcForces[iSubcase]
        
    def mpcForces(self):
        (subcaseName,iSubcase,transient,analysisCode) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int,str,float,float,float,float,float,float])

        if iSubcase in self.MpcForces:
            self.MpcForces[iSubcase].addData(data)
        else:
            self.MpcForces[iSubcase] = DisplacementObject(iSubcase,data)
        self.iSubcases.append(iSubcase)
        #print self.SpcForces[iSubcase]

    def readTable(self,Format):
        """reads displacement, spc/mpc forces"""
        sline = True
        data = []
        while sline:
            sline = self.infile.readline()[1:].strip().split()
            self.i+=1
            sline = self.parseLine(sline,Format)
            if sline is None:
                return data
            data.append(sline)
        return data
    
    def parseLine(self,sline,Format):
        out = []
        for entry,iFormat in zip(sline,Format):
            try:
                entry2 = iFormat(entry)
            except:
                #print "sline=|%s|\n entry=|%s| format=%s" %(sline,entry,iFormat)
                return None
            out.append(entry2)
        return out
        
    def ReadF06(self):
        #print "reading..."
        blank = 0
        while 1:
            if self.i%1000==0:
                print "i=%i" %(self.i)
            line = self.infile.readline()
            marker = line[1:].strip()
            
            #print "marker = %s" %(marker)
            if marker in self.markers:
                blank = 0
                print "\n*marker = %s" %(marker)
                self.markerMap[marker]()
                self.storedLines = []
            elif marker=='':
                blank +=1
                if blank==20:
                    break
            elif self.isMarker(marker): # marker with space in it (e.g. Model Summary)
                print "***marker = |%s|" %(marker)
                
            else:
                blank = 0
            ###
            self.storedLines.append(line)
            self.i+=1
        print "i=%i" %(self.i)
        self.infile.close()
        f06.processF06()

    def processF06(self):
        #data = [self.disp,self.SpcForces,self.stress,self.isoStress,self.barStress,self.solidStress,self.temperature]
        dataPack = [self.solidStress]
        for dataSet in dataPack:
            for key,data in dataSet.items():
                data.processF06Data()
        
    def isMarker(self,marker):
        """returns True if the word follows the 'N A S T R A N   P A T T E R N'"""
        marker = marker.strip().split('$')[0].strip()

        if len(marker)<2 or marker=='* * * * * * * * * * * * * * * * * * * *':
            return False
        for i,char in enumerate(marker):
            #print "i=%s i%%2=%s char=%s" %(i,i%2,char)
            if i%2==1 and ' ' is not char:
                return False
            elif i%2==0 and ' '==char:
                return False
            ###
        ###
        return True
            
    def skip(self,iskip):
        for i in range(iskip-1):
            self.infile.readline()
        self.i += iskip
        return self.infile.readline()
    
    def __repr__(self):
        msg = ''
        data = [self.disp,self.SpcForces,self.stress,self.isoStress,self.barStress,self.solidStress,self.temperature]
        data = [self.disp,self.solidStress,self.temperature]
        self.iSubcases = list(set(self.iSubcases))
        for iSubcase in self.iSubcases:
            for result in data:
                if iSubcase in result:
                    msg += str(result[iSubcase])
                ###
            ###
        return msg

if __name__=='__main__':
    f06 = F06Reader('cylinder01.f06')
    #f06 = F06Reader('ssb.f06')
    f06.ReadF06()
    print f06

