import os
import sys

class EndOfFileError(Exception):
    pass

class asdfF06Reader(object):
    def __init__(self,f06):
        self.f06 = f06
        #self.f06 = 'quad4_tri3.f06'
        #self.f06 = 'quad_bending.f06'
        self.infile = open(self.f06,'r')
        self.n = 0
        self.elements = {}

    def removeComments(self,lines):
        lines2 = []
        for line in lines:
            if '$' not in line:
                lines2.append(line)
        return lines2

    def readHeader(self):
        self.nHeader = 50
        wordStart = 'N A S T R A N    F I L E    A N D    S Y S T E M    P A R A M E T E R    E C H O'
        wordExec  = 'N A S T R A N    E X E C U T I V E    C O N T R O L    E C H O'
        wordCase  = 'C A S E    C O N T R O L    E C H O '
        wordEnd   = 'M O D E L   S U M M A R Y'

        line = ''
        
        while wordStart not in line:
            line = self.readlines(nLines=1)
        
        linesFileSystem = []
        linesExecControl = []
        linesCaseControl = []

        while wordExec not in line:
            linesFileSystem.append(line)
        while wordCase not in line:
            linesExecControl.append(line)
        while wordEnd not in line:
            linesCaseControl.append(line)
        
        linesFileSystem  = self.removeComments(linesFileSystem)
        linesExecControl = self.removeComments(linesExecControl)
        linesCaseControl = self.removeComments(linesCaseControl)
        
        #for line in linesCaseControl:
        #    print line.strip()

        
    def readDisplacement(self):
        word = 'D I S P L A C E M E N T   V E C T O R'
        typeList=[None,float,float,float]
        table = {}
        i=0
        nHeader=1

        while 1:
            try:
                table = self.readTable(word,table,typeList,nHeader,self.parseTableLine)
            except EndOfFileError:
                break
            #print "i=",i
            i+=1
        return table

    def printTable(self,table):
        for key,value in sorted(table.items()):
            print "  key=%s value=%s" %(key,value)

    def readTable(self,word,table,typeList,nHeader,parseTableLine):
        #nHeader=1
        tableStart = self.readTableHeader(word,typeList,nHeader,parseTableLine)
        if table=={}:
            table=tableStart
        line = self.readlines(nLines=1)
        while('PAGE' not in line):
            C = parseTableLine(line)
            #print C
            tableValues = []
            
            for itype,value in zip(typeList,C[1:]):
                if itype:
                    tableValues.append(itype(value))
            id = int(C[0])
            table[id] = tableValues
            line = self.readlines(nLines=1)
        #print "C=%s" %(C)
        #print "line = |%s|" %(line)
        return table

    def rewind(self):
        self.infile.close()
        self.infile = open(self.f06,'r')
        self.n = 0

    def readSPCForces(self):
        word = 'F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T'
        typeList=[None,float,float,float]
        table = self.readTable(word,typeList,self.parseTableLine)

    def readlines(self,nLines=1,saveAll=False):
        lines = []
        for line in range(nLines):
            line = self.infile.readline()
            if saveAll:
                lines.append(line)
            #print "line[%s]=%s" %(self.n,line.strip())
            ###
            self.n += 1
            #if self.n>347:
            #    raise Exception('too far')
        ###
        if '* * * END OF JOB * * *' in line:
            raise EndOfFileError()
        elif saveAll:
            return lines
        else:
            return line

#    def parseSummaryTableLine(self,line):
#        a,b = line[1:15],line[15:27] # ID, type
#        c,d,e = line[27:40],line[40:56 ],line[56 :72 ]  # T1,T2,T3
#        f,g,h = line[72:86],line[86:101],line[101:110]  # R1,R2,R3
#
#        A = [a,b,c,d,e,f,g,h]
#        A = [a.strip() for a in A]
#        return A

    def parseTriTableLine(self,line):
        a,b   = line[1:10],line[10:31] # ID, Distance
        c,d,e = line[31:46],line[46:60 ],line[60 :76 ]  # ox,oy,Txy
        f,g,h = line[76:88],line[88:103],line[103:120]  # angle,major,minor
        i   = line[120:132]  # VonMises

        A = [a,b, c,d,e, f,g,h, i]
        A = [a.strip() for a in A]
        #print "A = ",A
        #sys.exit('stopping')
        return A

    def parseQuadTableLine(self,line):
        a,b   = line[1:14],line[14:22] # ID, Grid-ID
        c,d,e = line[22:36],line[36:51 ],line[51 :65 ]  # Distance,ox,oy
        f,g,h = line[65:81],line[81:92],line[92:105]  # Txy,angle,major
        i,j   = line[105:120],line[120:135]  # minor,VonMises

        A = [a,b, c,d,e, f,g,h, i,j]
        A = [a.strip() for a in A]
        #print "A = ",A
        #sys.exit('stopping')
        return A

    def parseTableLine(self,line):
        a,b =   line[0:15],line[15:26] # ID, type
        c,d,e = line[26:39],line[39:55 ],line[55 :71 ]  # T1,T2,T3
        f,g,h = line[71:85],line[85:100],line[100:109]  # R1,R2,R3
        
        #print "ID=%s type=%s" %(a,b)
        #print "T1=%s T2=%s T3=%s" %(c,d,e)
        #print "R1=%s R2=%s R3=%s" %(f,g,h)
        A = [a,b,c,d,e,f,g,h]
        A = [a.strip() for a in A]
        return A

    def readTableHeader(self,word,typeList,nHeader,parseTableLine):
        table = {}
        line = ''
        nTypes = len(typeList)
        #print 'word=|%s|' %(word)
        while(word not in line):
            #print "line=|%s|" %(line.strip())
            line = self.readlines(nLines=1)
            #print "word = ",word
            #print "line = ",line
            isWord = word in line
            #print "isWord = ",isWord
            #print ""
        #print "found...%s" %(word)
        A = parseTableLine(line) # Displacement Vector
        #print "A=%s" %(A)
        #print "line = |%s|" %(line)

        line = self.readlines(nLines=1+nHeader)
        B = parseTableLine(line) # All Variable Titles
        #print "B=%s" %(B)
        #print "line = |%s|" %(line)

        Blist = [] # Blist is the reduced set of Variable Titles
        for bi,itype in zip(B[1:],typeList):
            if itype:
                Blist.append(bi)
        table = {'headers':Blist}
        return table

    def readElementTable(self,word,typeList,nNodes,nHeader,nBlank,parseTableLine):
        nHeader-=1
        table = self.readTableHeader(word,typeList,nHeader,parseTableLine)
        #print "table = ",table
        line = self.readlines(nLines=1)
        #return table
        while('JOB' not in line):
            tableNodes = []
            #print "-----------"
            for iNode in range(nNodes):
                layers = []
                layer = []
                #print "line1 = ",line.strip()
                C1 = parseTableLine(line)

                line = self.readlines(nLines=1)
                #print "line2 = ",line.strip()
                C2 = parseTableLine(line)
                #print "C2=",C2
                if iNode==0:
                    id = C1[0]
                    #print "id=",id
                #print C
                for itype,value in zip(typeList,C1[1:]):

                    if itype:
                        #print "C1 itype=%s value=%s" %(itype,value)
                        layer.append(itype(value))
                #print "tvo_1 = ",layer
                #print ""
                layers.append(layer)

                layer = []
                for itype,value in zip(typeList,C2[1:]):
                    if itype:
                        #print "C2 itype=%s value=%s" %(itype,value)
                        layer.append(itype(value))
                layers.append(layer)
                #print "tvo_2 = ",layer,'\n'
                #print "line*** = ",line.strip()

                #print "lineNew = ",line.strip()
                if nBlank>0:
                    line = self.readlines(nLines=nBlank)
                line = self.readlines(nLines=1) # quad
                tableNodes.append(layers)

            #print "lineNew = ",line.strip()
            #print "tvn = ",tableValuesNodes
            #print "id = ",id

            table[id] = tableNodes
            #print "line*** = ",line.strip()
            if 'JOB' in line:
                break
            #sys.exit()
        print "-----------"
        #print "C=%s" %(C)
        #print "line = |%s|" %(line)

        #for key,value in sorted(table.items()):
        #    print "key=%s " %(key)
        #    for row in value:
        #        print "   value=%s" %(row)

        return table


        pass

    def readSummary(self):
        word = 'OLOAD    RESULTANT'
        word = 'SPCFORCE RESULTANT'
        word = 'MAXIMUM  SPCFORCES'
        word = 'MAXIMUM  DISPLACEMENTS'
        word = 'MAXIMUM  APPLIED LOADS'

    def readQuad4Stresses(self):
        word = 'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )'
        typeList = [str]+8*[float]
        typeList = [str,None,float,float,float,float]
        typeList = [None,None,str,str,str]
        nHeader = 3
        nNodes = 5
        nBlank = 1
        table = self.readElementTable(word,typeList,nNodes,nHeader,nBlank,self.parseQuadTableLine)
        #word = 'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )'

        for key,nodes in table.items():
            if key!='headers':
                self.elements[key] = ShellElement(key,nodes)


    def readTri3Stresses(self):
        word = 'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )'
        typeList = 8*[float]
        typeList = [None,None,float,str,str,str]
        typeList = [None,str,str,str] # good
        #typeList = [str,str]
        # [dist,ox,oy,Txy,angle,major,minor,vm]
        nHeader = 2
        nNodes = 1
        nBlank = 0
        table = self.readElementTable(word,typeList,nNodes,nHeader,nBlank,self.parseTriTableLine)

        for key,nodes in table.items():
            if key!='headers':
                self.elements[key] = ShellElement(key,nodes)


        #word = 'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )'

class ShellElement(object):
    def __init__(self,id,nodes):
        self.layers = []
        self.id = id
        print "id = ",id
        for layers in nodes:
            print "layers = ",layers
            self.nLayers = len(layers)
            for iLayer,layer in enumerate(layers):
                #print "iLayer=%s layer=%s" %(iLayer,layer)
                print "layer=",layer
                (ox,oy,Txy) = layer

                #Txy = 'Txy'
                layerObj = Layer(iLayer,ox,oy,Txy)
                self.layers.append(layerObj)
    def __repr__(self):
        out = ''
        for layer in self.layers:
            out += "%s" %(layer)
        return out+'\n'


class Layer(object):
    def __init__(self,iLayer,ox,oy,Txy):
        self.iLayer=iLayer+1
        self.ox = ox
        self.oy = oy
        self.Txy = Txy
    def __repr__(self):
        return "    iLayer=%s ox=%s oy=%s Txy=%s\n" %(self.iLayer,self.ox,self.oy,self.Txy)

    def major(self):
        pass
    def minor(self):
        pass
    def vm(self):
        pass

class Table(dict):
    def __init__(self,tableDict):
#        del tableDict['tableType'],tableDict['headers']
#        self.tableDict = tableDict

        #self._tableType = tableDict['tableType']
        self.headers = tableDict['headers']

        dict.__init__(self)
        self._keys = tableDict.keys()
        for key, value in tableDict.items():
            self[key] = value

    def getSpot(self,headerName):
        #print "headerName = ",headerName
        #loc = self.headers.index(header)
        try:
            loc = self.headers.index(headerName)
        except ValueError:
            msg = 'invalid header.  headerName=|%s| available headers=|%s|' %(headerName,self.headers)
            raise Exception(msg)
        return loc

    def getValue(self,idNum):
        return self[idNum]

    def getValueName(self,idNum,headerName):
        loc = self.getSpot(headerName)
        return self[idNum][loc]

    def getValueSpot(self,idNum,spot):
        """
        Returns the value in the Ith spot in the list of the dictionary
        Requires that you know the index.
        Should use the getSpot function if you're going to use this.
        """
        return self[idNum][spot]

    def printTable(self,table):
        for key,value in sorted(table.items()):
            print "  key=%s value=%s" %(key,value)


class DisplacementObject(object):
    def __init__(self,subcaseID,data):
        #print "*******"
        self.subcaseID = subcaseID
        self.grids = []
        self.gridTypes = []
        self.translations = []
        self.rotations = []
        self.addData(data)
        
    def addData(self,data):
        for line in data:
            (gridID,gridType,t1,t2,t3,t4,t5,t6) = line
            self.grids.append(gridID)
            self.gridTypes.append(gridType)
            self.translations.append([t1,t2,t3])
            self.rotations.append([t4,t5,t6])
        ###
        #print "grids = ",self.grids
    
    def __repr__(self):
        msg  = "SubcaseID = %s\n" %(self.subcaseID)
        msg += '%-8s %8s %10s %10s %10s %10s %10s %10s\n' %('gridID','gridType','t1','t2','t3','t4','t5','t6')
        for (gridID,gridType,translation,rotation) in zip(self.grids,self.gridTypes,self.translations,self.rotations):
            (t1,t2,t3) = translation
            (t4,t5,t6) = rotation
            msg += "%-8i %8s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n" %(gridID,gridType,t1,t2,t3,t4,t5,t6)
        ###
        return msg

class StressObject(object):
    def __init__(self,subcaseID,data):
        #print "*******"
        self.subcaseID = subcaseID
        self.grids = []
        self.gridTypes = []
        self.translations = []
        self.rotations = []
        
        self.data = []
        self.addData(data)
        
    def addData(self,data):
        
        for line in data:
            self.data.append(line)
        return
            #(gridID,gridType,t1,t2,t3,t4,t5,t6) = line
            #self.grids.append(gridID)
            #self.gridTypes.append(gridType)
            #self.translations.append([t1,t2,t3])
            #self.rotations.append([t4,t5,t6])
        ###
        #print "grids = ",self.grids
    
    def __repr__(self):
        msg  = "SubcaseID = %s\n" %(self.subcaseID)
        for line in self.data:
            msg += '%s\n' %(line)
        return msg
            
        msg += '%-8s %8s %10s %10s %10s %10s %10s %10s\n' %('gridID','gridType','t1','t2','t3','t4','t5','t6')
        for (gridID,gridType,translation,rotation) in zip(self.grids,self.gridTypes,self.translations,self.rotations):
            (t1,t2,t3) = translation
            (t4,t5,t6) = rotation
            msg += "%-8i %8s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n" %(gridID,gridType,t1,t2,t3,t4,t5,t6)
        ###
        return msg

class F06Reader(object):
    def __init__(self,f06name):
        self.f06name = f06name
        self.i = 0
        self.markerMap = {
          #'N A S T R A N    F I L E    A N D    S Y S T E M    P A R A M E T E R    E C H O':self.fileSystem,
          #'N A S T R A N    E X E C U T I V E    C O N T R O L    E C H O':self.executiveControl,
          #'C A S E    C O N T R O L    E C H O ':self.caseControl,
          #'M O D E L   S U M M A R Y':self.summary,
#                            E L E M E N T   G E O M E T R Y   T E S T   R E S U L T S   S U M M A R Y
#                           O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R
#0                                                  OLOAD    RESULTANT       
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

          #'* * * END OF JOB * * *': self.end(),
          'MAXIMUM  SPCFORCES':self.maxSpcForces,
          'MAXIMUM  DISPLACEMENTS': self.maxDisplacements,
          'MAXIMUM  APPLIED LOADS': self.maxAppliedLoads,
         }
        self.markers = self.markerMap.keys()
        self.infile = open(self.f06name,'r')
        self.storedLines = []

        self.disp = {}
        self.SpcForces = {}
        self.stress = {}
        self.subcaseIDs = []

    def gridPointSingularities(self):
        """
0                                         G R I D   P O I N T   S I N G U L A R I T Y   T A B L E
0                             POINT    TYPE   FAILED      STIFFNESS       OLD USET           NEW USET
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
        #self.disp[subcaseID] = DisplacementObject(subcaseID,data)
        #print self.disp[subcaseID]

    def maxDisplacements(self):
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int,float,float,float,float,float,float])
        print "max Displacements",data
        #self.disp[subcaseID] = DisplacementObject(subcaseID,data)
        #print self.disp[subcaseID]
    
    def maxAppliedLoads(self):
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int,float,float,float,float,float,float])
        print "max Applied Loads",data
        #self.disp[subcaseID] = DisplacementObject(subcaseID,data)
        #print self.disp[subcaseID]

    def quadCompositeStress(self):
        """
                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
0      181    1   3.18013E+04  5.33449E+05  1.01480E+03   -7.06668E+01  1.90232E+04   89.88  5.33451E+05  3.17993E+04  2.50826E+05
0      181    2   1.41820E+05  1.40805E+05  1.25412E+05   -1.06000E+02  2.85348E+04   44.88  2.66726E+05  1.58996E+04  1.25413E+05
        """
        subcaseName = self.storedLines[-3].strip()
        subcaseID   = self.storedLines[-2].strip()[1:]
        subcaseID = int(subcaseID.strip('SUBCASE '))
        #print "subcaseName=%s subcaseID=%s" %(subcaseName,subcaseID)
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int,int,float,float,float,float,float,float,float,float,float])
        if subcaseID in self.stress:
            self.stress[subcaseID].addData(data)
        else:
            self.stress[subcaseID] = StressObject(subcaseID,data)
        self.subcaseIDs.append(subcaseID)
        #print self.stress[subcaseID]
        

    def displacement(self):
        """
                                             D I S P L A C E M E N T   V E C T O R
 
      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
             1      G      9.663032E-05   0.0           -2.199001E-04   0.0           -9.121119E-05   0.0
             2      G      0.0            0.0            0.0            0.0            0.0            0.0
             3      G      0.0            0.0            0.0            0.0            0.0            0.0
        """
        subcaseName = self.storedLines[-3].strip()
        subcaseID   = self.storedLines[-2].strip()[1:]
        subcaseID = int(subcaseID.strip('SUBCASE '))
        #print "subcaseName=%s subcaseID=%s" %(subcaseName,subcaseID)
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int,str,float,float,float,float,float,float])
        if subcaseID in self.disp:
            self.disp[subcaseID].addData(data)
        else:
            self.disp[subcaseID] = DisplacementObject(subcaseID,data)
        self.subcaseIDs.append(subcaseID)
        #print self.disp[subcaseID]

    def spcForces(self):
        subcaseName = self.storedLines[-3].strip()
        subcaseID   = self.storedLines[-2].strip()[1:]
        subcaseID = int(subcaseID.strip('SUBCASE '))
        #print "subcaseName=%s subcaseID=%s" %(subcaseName,subcaseID)
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int,str,float,float,float,float,float,float])

        if subcaseID in self.SpcForces:
            self.SpcForces[subcaseID].addData(data)
        else:
            self.SpcForces[subcaseID] = DisplacementObject(subcaseID,data)
        self.subcaseIDs.append(subcaseID)
        #print self.SpcForces[subcaseID]
        
    def mpcForces(self):
        subcaseName = self.storedLines[-3].strip()
        subcaseID   = self.storedLines[-2].strip()[1:]
        subcaseID = int(subcaseID.strip('SUBCASE '))
        #print "subcaseName=%s subcaseID=%s" %(subcaseName,subcaseID)
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data = self.readTable([int,str,float,float,float,float,float,float])

        if subcaseID in self.MpcForces:
            self.MpcForces[subcaseID].addData(data)
        else:
            self.MpcForces[subcaseID] = DisplacementObject(subcaseID,data)
        self.subcaseIDs.append(subcaseID)
        #print self.SpcForces[subcaseID]
        
    def readTable(self,Format):
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
        #print "2*************"
        #print dir(self)
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
                #print "\n*marker = %s" %(marker)
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

    def isMarker(self,marker):
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
        data = [self.disp,self.SpcForces,self.stress]
        self.subcaseIDs = list(set(self.subcaseIDs))
        for subcaseID in self.subcaseIDs:
            for result in data:
                if subcaseID in result:
                    msg += str(result[subcaseID])
                ###
            ###
        return msg

if __name__=='__main__':
    f06 = F06Reader('ssb.f06')
    #f06 = F06Reader('quad_bending_comp.f06')
    f06.ReadF06()
    #print f06
    
if 0:
    reader = F06Reader('fem3.f06')
    #reader.readHeader()
    dTable = reader.readDisplacement()
    reader.printTable(dTable)
    
    sys.exit()

    #print "\n-QUAD Stresses-"
    #stable = reader.readQuad4Stresses()
    #reader.rewind()
    #print "\n"

    #print "\n-Tri Stresses-"
    #stable = reader.readTri3Stresses()
    #reader.rewind()

    #for key,element in reader.elements.items():
    #    print "element = ",element
    #    print "id=%s e=%s" %(key,element)
    #print "\n"

    print "-Displacements-"
    dtableTemp = reader.readDisplacement()
    reader.rewind()
    #reader.printTable(dtableTemp)

    #print "dtableTemp['headers'] = ",dtableTemp['headers']
    dTable = Table(dtableTemp)
    print "dTable.headers = ",dTable.headers
    reader.printTable(dTable)

    print "d4  = ",dTable.getValue(4)
    print "dx4 = ",dTable.getValueName(4,'T1')
    print "dz4 = ",dTable.getValueSpot(4,2)

    #print "-SPC Forces-"
    #ftable = reader.readSPCForces()

