import os
import sys

class EndOfFileError(Exception):
    pass

class F06Reader(object):
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


if __name__=='__main__':
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
