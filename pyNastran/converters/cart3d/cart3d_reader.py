import os
import sys
from math import ceil,sqrt
from numpy import array
from pyNastran.general.general import ListPrint

def convertToFloat(svalues):
    values = []
    for value in svalues:
        values.append(float(value))
    return values

#------------------------------------------------------------------

class Cart3DReader(object):
    modelType = 'cart3d'
    isStructured = False
    isOutwardNormals = True

    def __init__(self,log=None,debug=False):
        self.isHalfModel = True
        self.cartType  = None # grid, result
        self.nPoints   = None
        self.nElements = None
        self.infile = None
        self.infilename = None
        self.readHalf = False

        if log is None:
            from pyNastran.general.logger import dummyLogger
            loggerObj = dummyLogger()
            if debug:
                log = loggerObj.startLog('debug') # or info
            else:
                log = loggerObj.startLog('info') # or info
        self.log = log

    def readCart3d(self,infilename):
        """extracts the points, elements, and Cp"""
        self.infilename = infilename
        self.log.info("---starting reading cart3d file...|%s|---" %(self.infilename))
        self.infile = open(self.infilename,'r')
        self.readHeader()
        points = self.getPoints()
        elements = self.getElements(bypass=False)
        regions = self.getRegions(bypass=False)
        loads = self.readResults(0,self.infile)
        #if self.cartType=='results':
            
        #loads = self.getLoads(['Cp'])
        self.log.info("nPoints=%s  nElements=%s" %(self.nPoints,self.nElements))
        self.log.info("---finished reading cart3d file...|%s|---" %(self.infilename))
        return (points,elements,regions,loads)

    def getYZeroNodes(self,nodes):
        yZeroNodes  = {}
        outerNodes  = {}
        isInnerNode = {}
        
        #nodeTypes  = {}
        for (iNode,node) in sorted(nodes.items()):
            if node[1]==0:
                yZeroNodes[iNode] = node
                isInnerNode[iNode] = True
            else:
                outerNodes[iNode] = [node[0],-node[1],node[2]]
                isInnerNode[iNode] = False
            ###
        ###
        return (yZeroNodes,outerNodes,isInnerNode)

    def makeFullModel(self,nodes,elements,regions,loads):
        """assumes a y=0 line of nodes"""
        self.log.info('---starting makeFullModel---')
        
        maxNode     = len(nodes.keys())
        maxElements = len(elements.keys())
        (yZeroNodes,outerNodes,isInnerNode) = self.getYZeroNodes(nodes)

        #nodes2={}
        nodes2Map = {}
        iNode2=1
        for (iNode,node) in sorted(nodes.items()):
            if isInnerNode[iNode]:
                #pass
                #nodes2[iNode] = node  # dont need to set this...
                nodes2Map[iNode] = iNode
            else:
                #testElements2[maxNode+iNode2] = outerNodes[iNode]
                nodes[maxNode+iNode2] = outerNodes[iNode]
                nodes2Map[iNode]      = maxNode+iNode2
                iNode2+=1
            ###
        ###
        
        # based on the conversion (nodes2map) object, renumber each element
        # stick the result in the location for the new element
        #elementsB={}
        #regionsB={}

        nThree = 147120+10 #len(elements)*2   ## todo why is there +10???
        allSet = set([i+1 for i in range(nThree)])
        elementKeysSet = set(elements.keys())
        missingKeys    = list(allSet.difference(elementKeysSet))
        missingKeys.sort()
        #print "missingKeys = ",missingKeys
        #print "nMissing = ",len(missingKeys)
        iMissing = 0
        for (iElement,element) in sorted(elements.items()):
            #(n1,n2,n3) = element
            region = regions[iElement]
            element2 = []
            
            for nID in element:
                nID2 = nodes2Map[nID]
                element2.append(nID2)
            element2.reverse()  # the normal will be flipped otherwise
            
            eMissing = missingKeys[iMissing]
            #print "eMissing = ",eMissing
            elements[eMissing] = element2
            regions[eMissing]  = region
            iMissing += 1
        ###
        
        for i in range(1,len(elements)):
            assert i in elements,'i=%s is not in elements...' %(i)
            #assert i in regions ,'i=%s is not in regions...' %(i)
        self.log.info('---finished makeFullModel---')
        sys.stdout.flush()
        
        return (nodes,elements,regions,loads)

    def makeHalfModel(self,nodes,elements,regions,Cp):
        nNodesStart = len(nodes)

        self.log.info('---starting makeHalfModel---')
        for (iNode,node) in sorted(nodes.items()):
            if node[1]<0:
                #print iNode,'...',node
                del nodes[iNode]
                try:
                    del Cp[iNode]  # assume loads=0
                except:
                    pass
                ###
            ###
        sys.stdout.flush()

        for (iElement,element) in elements.items():
            #print "iElement = |%s|" %(iElement),type(iElement)
            inNodeList = [True,True,True]
            for (i,node) in enumerate(element):
                if node not in nodes:
                    inNodeList[i] = False
                ###
            if False in inNodeList:
                #print iElement,element,inNodeList,type(iElement)
                del elements[iElement]
                del regions[iElement]
            ###
        ###
        assert len(nodes)!=nNodesStart
        self.log.info('---finished makeHalfModel---')
        return (nodes,elements,regions,Cp)

    def makeMirrorModel(self,nodes,elements,regions,loads):
        nNodesStart = len(nodes)

        self.log.info('---starting makeMirrorModel---')
        for (iNode,node) in sorted(nodes.items()):
            if node[1]<0:
                #print iNode,'...',node
                del nodes[iNode]
                #del Cp[iNode]  # assume loads=0
        sys.stdout.flush()

        for (iElement,element) in elements.items():
            #print "iElement = |%s|" %(iElement),type(iElement)
            inNodeList = [True,True,True]
            for (i,node) in enumerate(element):
                if node not in nodes:
                    inNodeList[i] = False
                ###
            if False in inNodeList:
                #print iElement,element,inNodeList,type(iElement)
                del elements[iElement]
                del regions[iElement]
            ###
        ###
        #assert len(nodes)!=nNodesStart
        self.log.info('---finished makeMirrorModel---')
        
        return self.makeFullModel(nodes,elements,regions,loads)
        #return self.renumberMesh(nodes,elements,regions,loads)

    def renumberMesh(self,nodes,elements,regions,loads):
        iNodeCounter    = 1
        iElementCounter = 1
        
        NodeOldToNew    = {}
        ElementOldToNew = {}
        nodes2    = {}
        elements2 = {}
        regions2  = {}
        for (iNode,node) in sorted(nodes.items()):
            NodeOldToNew[iNode] = iNodeCounter
            nodes2[iNodeCounter] = node
            iNodeCounter += 1
        
        
        for (iElement,element) in sorted(elements.items()):
            element2 = []
            for nID in element:
                nID2 = NodeOldToNew[nID]
                element2.append(nID2)
            elements2[iElementCounter] = element2
            regions2[iElementCounter] = regions[iElement]
            iElementCounter += 1
            
        return (nodes,elements,regions,loads)

    def writeOutfile(self,outfilename,points,elements,regions):
        assert len(points)>0
        self.log.info("---writing cart3d file...|%s|---" %(outfilename))
        outfile = open(outfilename,'wb')

        msg = self.writeHeader(f,points,elements)
        msg = self.writePoints(f,points)
        msg = self.writeElements(f,points,elements)
        msg = self.writeRegions(f,regions)
        outfile.close()

    def writeRegions(self,f,regions):
        msg = ''
        for (iElement,region) in sorted(regions.items()):
            msg += "%s\n" %(region)
        f.write(msg)

    def writeElements(self,f,nodes,elements):
        msg = ''
        strElements = len(str(len(elements)))
        format = " %%%ss " %(strElements)
        #format = '%-10s '+format*3 + '\n'
        #format = format*3 + '\n'
        #print "format = ",format

        maxNodes = []
        for (iElement,element) in sorted(elements.items()):
            #self.log.info("element = %s" %(element))
            #msg += format %(iElement,element[0],element[1],element[2])
            msg  += "%s  %s  %s\n" %(element[0],element[1],element[2])
            for nID in element:
                assert nodes.has_key(nID),'nID=%s is missing' %(nID)

            maxNode = max(element)
            maxNodes.append(maxNode)
        self.log.info("maxNodeID = %s" %(max(maxNodes)))
        f.write(msg)

    def writePoints(self,f,points):
        for (iPoint,point) in sorted(points.items()):
            #msg += "%-10s %-15s %-15s %-15s\n" %(iPoint,point[0],point[1],point[2])
            msg  += "%6.6f  %6.6f  %6.6f\n" %(point[0],point[1],point[2])
            if iPoint%1000:
                f.write(msg)
                msg = ''
        f.write(msg)
        
    def writeHeader(self,f,points,elements):
        nPoints = len(points)
        nElements = len(elements)
        msg = "%s %s\n" %(nPoints,nElements)
        f.write(msg)

    def readHeader(self):
        line = self.infile.readline()
        sline = line.strip().split()
        if len(sline)==2:
            self.nPoints,self.nElements = int(sline[0]),int(sline[1])
            self.cartType = 'grid'
        elif len(sline)==3:
            self.nPoints,self.nElements,self.nResults = int(sline[0]),int(sline[1]),int(sline[2])
            self.cartType = 'result'
        else:
            raise Exception('invalid result type')

        if self.readHalf:
            self.nElementsRead = self.nElements/2
            self.nElementsSkip = self.nElements/2
        else:
            self.nElementsRead = self.nElements
            self.nElementsSkip = 0

    def getPoints(self):
        """
        A point is defined by x,y,z and the ID is the location in points.
        """
        points = {}
        for p in range(self.nPoints):
            x,y,z       = self.infile.readline().strip().split()
            points[p+1] = array([float(x),float(y),float(z)])
        
        maxX = self.getMax(points,0)
        maxY = self.getMax(points,1)
        maxZ = self.getMax(points,2)

        minX = self.getMin(points,0)
        minY = self.getMin(points,1)
        minZ = self.getMin(points,2)
        
        self.log.info("X  max=%g min=%g" %(maxX,minX))
        self.log.info("Y  max=%g min=%g" %(maxY,minY))
        self.log.info("Z  max=%g min=%g" %(maxZ,minZ))
        return points

    def getMin(self,points,i):
        minX = points[1][i]
        for key,point in points.items():
            minX = min(minX,point[i])
        return minX

    def getMax(self,points,i):
        maxX = points[1][i]
        for key,point in points.items():
            maxX = max(maxX,point[i])
        return maxX

    def getElements(self,bypass=False):
        """
        An element is defined by n1,n2,n3 and the ID is the location in elements.
        """
        assert bypass==False
        elements = {}

        print "nElementsRead=%s nElementsSkip=%s" %(self.nElementsRead,self.nElementsSkip)
        for e in range(self.nElementsRead):
            element = self.infile.readline().strip().split()
            element = [int(spot) for spot in element] # converts eid strings into ints
            elements[e+1] = element

        for e in range(self.nElementsSkip):
            element = self.infile.readline()

        return elements

    def getRegions(self,bypass=True):
        regions = {}
        if bypass:
            for i in range(self.nElements):
                self.infile.readline()
        else:
            for i in range(self.nElementsRead):
                line = self.infile.readline()
                regions[i+1] = line.strip()
            ###
            for i in range(self.nElementsSkip):
                line = self.infile.readline()
        ###
        return regions
    
    def getLoads(self,outputs):   ## TODO:  outputs isnt used...
        loads = readResults(self,i,self.infile)
        #Cp = [0.1]*self.nElements ## TODO:  Cp is hardcoded to 0.1
        loads = {'Cp':Cp,'rho':[],'U':[],'V':[],'W':[],'pressure':[]}
        raise Exception('DEPRECIATED...')
        return loads

    def readResults(self,i,infile):
        """
        Reads the Cp results.
        Cp = (p - 1/gamma) / (0.5*M_inf*M_inf)
        (pV)^2 = (pu)^2+(pv)^2+(pw)^2
        M^2 = (pV)^2/rho^2

        # ???
        p = (gamma-1)*(e- (rhoU**2+rhoV**2+rhoW**2)/(2.*rho))
        
        # ???
        rho,rhoU,rhoV,rhoW,rhoE 
        """
        Cp = {}
        Mach = {}
        U = {}
        V = {}
        W = {}
        if(self.cartType=='result'):
            nResultLines = int(ceil(self.nResults/5.))-1
            #print "nResultLines = ",nResultLines
            for pointNum in range(self.nPoints):
                #print "pointNum = ", pointNum
                sline = infile.readline().strip().split() # rho rhoU,rhoV,rhoW,pressure/rhoE/E
                i+=1
                for n in range(nResultLines):
                    #print "nResultLines = ",n
                    sline+=infile.readline().strip().split() # Cp
                    #infile.readline()
                    i+=1
                    #print "sline = ",sline
                    #gamma = 1.4
                    #else:
                    #    p=0.
                sline = self.getList(sline)

                #p=0
                cp = sline[0]
                rho = float(sline[1])
                if(rho>abs(0.000001)):
                    rhoU = float(sline[2])
                    rhoV = float(sline[3])
                    rhoW = float(sline[4])
                    #rhoE = float(sline[5])
                    VV = (rhoU)**2+(rhoV)**2+(rhoW)**2/rho**2
                    mach = sqrt(VV)
                    if mach>10:
                        print "nid=%s Cp=%s mach=%s rho=%s rhoU=%s rhoV=%s rhoW=%s" %(pointNum,cp,mach,rho,rhoU,rhoV,rhoW)

                Cp[pointNum+1] = cp
                Mach[pointNum+1] = mach
                U[pointNum+1] = rhoU/rho
                V[pointNum+1] = rhoV/rho
                W[pointNum+1] = rhoW/rho
                #print "pt=%s i=%s Cp=%s p=%s" %(pointNum,i,sline[0],p)
            ###
        ###
        loads = {'Cp':Cp,'Mach':Mach,'U':U,'V':V,'W':W}
        return loads
    ###

    def getList(self,sline):
        """Takes a list of strings and converts them to floats."""
        try:
            sline2 = convertToFloat(sline)
        except ValueError:
            print "sline = %s" %sline
            raise Exception('cannot parse %s' %(sline))
        return sline2

    def getValue(self,sline):
        """Converts a string to a float."""
        value = float(sline[1])
        return value
    ###
    
    def exportToNastran(self,fname,points,elements,regions):
        from pyNastran.bdf.fieldWriter import printCard
        
        f = open(fname,'wb')
        for nid,grid in enumerate(1,points):
            (x,y,z) = grid
            f.write(printCard(['GRID',nid,'',x,y,z]))
        ###

        e = 1e7
        nu = 0.3
        g = ''
        thickness = 0.25
        setRegions = list(set(regions))
        for pidMid in setRegions:
            f.write(printCard(['MAT1',pidMid,e,g,nu]))
            f.write(printCard(['PSHELL',pidMid,pidMid,thickness]))
        ###

        for eid,(nodes,region) in enumerate(1,zip(elements,regions)):
            (n1,n2,n3) = nodes
            f.write(printCard(['CTRIA3',eid,region,n1,n2,n3]))
        ###
        f.close()

#------------------------------------------------------------------

if __name__=='__main__':
    basepath = os.getcwd()
    configpath = os.path.join(basepath,'inputs')
    workpath   = os.path.join(basepath,'outputsFinal')
    #bdfModel   = os.path.join(configpath,'aeroModel.bdf')
    #assert os.path.exists(bdfModel),'|%s| doesnt exist' %(bdfModel)
    os.chdir(workpath)
    print "basepath",basepath


    cart3dGeom  = os.path.join(configpath,'Cart3d_bwb2.i.tri')
    cart3dGeom2 = os.path.join(workpath,'Cart3d_half.i.tri')
    cart3dGeom3 = os.path.join(workpath,'Cart3d_full.i.tri')
    #cart3dGeom4 = os.path.join(workpath,'Cart3d_full2.i.tri')

    cart = Cart3DReader(cart3dGeom)
    (points,elements,regions,loads) = cart.read()
    (points,elements,regions,loads) = cart.makeHalfModel(points,elements,regions,loads)
    cart.writeOutfile(cart3dGeom2,points,elements,regions)

    #cart = Cart3DReader(cart3dGeom)  # bJet.a.tri
    #(cartPoints,elements,regions,loads) = cart.read()
    #points2 = {}
    #for (iPoint,point) in sorted(points.items()):
        #(x,y,z) = point
        #print "p=%s x=%s y=%s z=%s  z2=%s" %(iPoint,x,y,z,z+x/10.)
        #points2[iPoint] = [x,y,z+x/10.]
    #(points,elements,regions,loads) = cart.makeMirrorModel(points2,elements,regions,loads)
    
    cart2 = Cart3DReader(cart3dGeom2)
    (points,elements,regions,loads) = cart2.read()
    
    #makeFullModel
    (points,elements,regions,loads) = cart2.makeFullModel(points,elements,regions,loads)  # half model values going in
    cart2.writeOutfile(cart3dGeom3,points,elements,regions)
    
    #cart3 = Cart3DReader(cart3dGeom2)
    #(points,elements,regions,loads) = cart3.read()
    #(points,elements,regions,loads) = cart3.makeMirrorModel(points,elements,regions,loads)
    

    #print "loads = ",ListPrint(loads),len(loads)
    #cartOutfile = os.path.join(workpath,'bJet.a.tri_test')
    #cart.writeInfile(cartOutfile,cartPoints,elements,regions)
    