from pyNastran.op2.tables.oug.oug_displacements import displacementObject
from pyNastran.op2.tables.oug.oug_temperatures  import temperatureObject

class OUG(object):
    def __init__(self):
        pass

    def displacement(self):
        """
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
        """
        (subcaseName,iSubcase,transient,analysisCode) = self.readSubcaseNameID()
        headers = self.skip(2)
        dataCode = {'log':self.log,'analysisCode':analysisCode,'deviceCode':1,'tableCode':1,
                    'sortCode':0,'sortBits':[0,0,0],'numWide':8}
        #print "headers = %s" %(headers)
        data = self.readTable([int,str,float,float,float,float,float,float])
        if iSubcase in self.disp:
            self.disp[iSubcase].addF06Data(data,transient)
        else:
            #self.disp[iSubcase] = DisplacementObject(iSubcase,data)
            disp = displacementObject(dataCode,iSubcase)
            disp.addF06Data(data,transient)
            self.disp[iSubcase] = disp
        self.iSubcases.append(iSubcase)

    def temperatureVector(self):
        """
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
        """
        (subcaseName,iSubcase,transient,analysisCode) = self.readSubcaseNameID()
        #print transient
        
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTemperatureTable()
        
        dataCode = {'log':self.log,'analysisCode':1,'deviceCode':1,'tableCode':1,'sortCode':0,
                    'sortBits':[0,0,0],'numWide':8,
                    #'formatCode':1,
                    #'elementName':eType,'sCode':0,'stressBits':stressBits
                    }

        if iSubcase in self.temperature:
            self.temperature[iSubcase].addF06Data(data,transient)
        else:
            temp = temperatureObject(dataCode,iSubcase)
            temp.addF06Data(data,transient)
            self.temperature[iSubcase] = temp
        self.iSubcases.append(iSubcase)

    def readTemperatureTable(self):
        data = []
        Format = [int,str,float,float,float,float,float,float]
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            if 'PAGE' in line:
                return data
            sline = [line[0:15],line[15:22].strip(),line[22:40],line[40:55],line[55:70],line[70:85],line[85:100]]
            sline = self.parseLineTemperature(sline,Format)
            data.append(sline)
        return data

    def parseLineTemperature(self,sline,Format):
        out = []
        for entry,iFormat in zip(sline,Format):
            if entry is '':
                return out
            #print "sline=|%r|\n entry=|%r| format=%r" %(sline,entry,iFormat)
            entry2 = iFormat(entry)
            out.append(entry2)
        return out
    
