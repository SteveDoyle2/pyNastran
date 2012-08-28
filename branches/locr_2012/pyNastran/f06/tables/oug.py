#pylint: disable=E1101,W0612,R0201
from itertools import izip

from pyNastran.op2.tables.oug.oug_displacements import DisplacementObject, ComplexDisplacementObject
from pyNastran.op2.tables.oug.oug_temperatures import TemperatureObject


class OUG(object):
    def __init__(self):
        self.displacements = {}
        self.temperatures = {}

    def getDisplacement(self):
        """
        @code
                                             D I S P L A C E M E N T   V E C T O R

        POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
               1      G      9.663032E-05   0.0           -2.199001E-04   0.0           -9.121119E-05   0.0
               2      G      0.0            0.0            0.0            0.0            0.0            0.0
               3      G      0.0            0.0            0.0            0.0            0.0            0.0

        analysisCode = 1 (Statics)
        deviceCode   = 1 (Print)
        tableCode    = 1 (Displacement)
        sortCode     = 0 (Sort2,Real,Sorted Results) => sortBits = [0,0,0]
        numWide      = 8 (???)
        @endcode
        """
        (subcaseName, iSubcase, transient, dt, analysisCode,
            isSort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        dataCode = {'log': self.log, 'analysisCode': analysisCode,
                    'deviceCode': 1, 'tableCode': 1,
                    'sortCode': 0, 'sortBits': [0, 0, 0], 'numWide': 8,
                    'tableName': 'OUG', 'nonlinearFactor': dt}
        #print "headers = %s" %(headers)
        dataTypes = [int, str, float, float, float, float, float, float]
        data = self.readTable(dataTypes)

        if iSubcase in self.displacements:
            self.displacements[iSubcase].addF06Data(data, transient)
        else:
            isSort1 = True
            disp = DisplacementObject(dataCode, isSort1, iSubcase)
            disp.addF06Data(data, transient)
            self.displacements[iSubcase] = disp
        self.iSubcases.append(iSubcase)

    def getComplexDisplacement(self):
        """
        @code
          BACKWARD WHIRL
                                                                                                                 SUBCASE 2
          POINT-ID =       101
                                           C O M P L E X   D I S P L A C E M E N T   V E C T O R
                                                             (MAGNITUDE/PHASE)

          FREQUENCY   TYPE          T1             T2             T3             R1             R2             R3
        2.000000E+01     G      3.242295E-16   1.630439E-01   1.630439E-01   1.691497E-17   1.362718E-01   1.362718E-01
                                196.0668        90.0000       180.0000        63.4349       180.0000       270.0000

        tableCode    = 1 (Displacement)
        formatCode   = 3 (Magnitude/Phase)
        sortBits     = [0,1,1]  (Sort1,Real/Imaginary,RandomResponse)
        analysisCode = 5 (Frequency)
        sortCode     = 2 (Random Response)
        @endcode
        """
        (subcaseName, iSubcase, transient, dt, analysisCode,
            isSort1) = self.readSubcaseNameID()
        headers = self.skip(3)
        data = []

        dataCode = {'log': self.log, 'analysisCode': 5, 'deviceCode': 1,
                    'tableCode': 1, 'sortCode': 2, 'sortBits': [0, 1, 1],
                    'numWide': 14, 'formatCode': 3, 'tableName': 'OUGV1',
                    'nonlinearFactor': dt,
                    #'mode':iMode,'eigr':transient[1], 'modeCycle':cycle,
                    #'dataNames':['mode', 'eigr', 'modeCycle'],
                    #'name':'mode',
                    #'sCode':0,
                    #'elementName':'CBAR', 'elementType':34, 'stressBits':stressBits,
                    }
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            if 'PAGE' in line:
                break
            sline = line.strip().split()
            line2 = self.infile.readline()[1:].rstrip('\r\n ')
            sline += line2.strip().split()
            out = [float(sline[0]), sline[1].strip(),
                   float(sline[2]), float(sline[3]), float(sline[4]),
                   float(sline[5]), float(sline[6]), float(sline[7]),
                   float(sline[8]), float(sline[9]), float(sline[10]),
                   float(sline[11]), float(sline[12]), float(sline[13]), ]
            #print sline
            #print out
            data.append(out)
            self.i += 2
        ###
        self.i += 1

        if iSubcase in self.displacements:
            self.displacements[iSubcase].addF06Data(data, transient)
        else:
            isSort1 = True
            disp = ComplexDisplacementObject(dataCode, isSort1, iSubcase)
            disp.addF06Data(data, transient)
            self.displacements[iSubcase] = disp
        self.iSubcases.append(iSubcase)

    def getTemperatureVector(self):
        """
        @code
        LOAD STEP =  1.00000E+00
                                              T E M P E R A T U R E   V E C T O R

        POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE
               1      S      1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03
               7      S      1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03
        analysisCode = 1 (Statics)
        deviceCode   = 1 (Print)
        tableCode    = 1 (Displacement/Temperature)
        sortCode     = 0 (Sort2,Real,Sorted Results) => sortBits = [0,0,0]
        formatCode   = 1 (Real)
        sCode        = 0 (Stress)
        numWide      = 8 (???)
        @endcode
        """
        (subcaseName, iSubcase, transient, dt, analysisCode,
            isSort1) = self.readSubcaseNameID()
        #print transient

        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTemperatureTable()

        dataCode = {'log': self.log, 'analysisCode': 1, 'deviceCode': 1,
                    'tableCode': 1, 'sortCode': 0, 'sortBits': [0, 0, 0],
                    'numWide': 8, 'tableName': 'OUG', 'nonlinearFactor': dt,
                    #'thermalCode':1,
                    #'formatCode':1,
                    #'elementName':eType,'sCode':0,'stressBits':stressBits
                    }

        if iSubcase in self.temperatures:
            self.temperatures[iSubcase].addF06Data(data, transient)
        else:
            isSort1 = True
            temp = TemperatureObject(dataCode, isSort1, iSubcase)
            temp.addF06Data(data, transient)
            self.temperatures[iSubcase] = temp
        self.iSubcases.append(iSubcase)

    def readTemperatureTable(self):
        data = []
        Format = [int, str, float, float, float, float, float, float]
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            if 'PAGE' in line:
                return data
            sline = [line[0:15], line[15:22].strip(), line[22:40], line[40:55], line[55:70], line[70:85], line[85:100], line[100:115]]
            sline = self.parseLineTemperature(sline, Format)
            data.append(sline)
        return data

    def parseLineTemperature(self, sline, Format):
        out = []
        for (entry, iFormat) in izip(sline, Format):
            if entry is '':
                return out
            #print "sline=|%r|\n entry=|%r| format=%r" %(sline,entry,iFormat)
            entry2 = iFormat(entry)
            out.append(entry2)
        return out
