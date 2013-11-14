
from math import radians,degrees,cos,sin
from pyNastran.bdf.bdf import BDF
import copy

def makeCircle(r,theta0,theta1):
    """
    Makes a circle from theta0 to theta1 at a radius R.
    """
    theta0 = radians(theta0)
    theta1 = radians(theta1)
    npoints = 10
    dtheta = (theta1-theta0)/npoints

    theta = theta0
    X=[]; Y=[]
    for i in range(npoints+1):
        x = r*cos(theta)
        y = r*sin(theta)
        print "x=%g \ty=%g     \ttheta=%s" %(x,y,degrees(theta))
        theta += dtheta
        X.append(x); Y.append(y)
    return (X,Y)

def readAirfoil(filename):
    infile = open(filename,'rU')
    lines = infile.readlines()
    (nUpper,nLower) = lines[1].split()
    nUpper=int(float(nUpper)); nLower=int(float(nLower))
    print "nUpper=%s nLower=%s" %(nUpper,nLower)

    upperSurface=[]; lowerSurface=[]
    upperLines = lines[3:3+nUpper]
    for line in upperLines:
        x,y = line.split()
        x=float(x); y=float(y)
        upperSurface.append([x,y])

    lowerLines = lines[4+nUpper:4+nUpper+nLower+1]
    for line in lowerLines:
        x,y = line.split()
        x=float(x); y=float(y)
        lowerSurface.append([x,y])

    for u in upperSurface:
        print "%g \t%g" %(u[0],u[1])

#    for u in lowerSurface:
#        print "%g \t%g" %(u[0],u[1])

    return(upperSurface,lowerSurface)


class BdfToP3d(object):
    def __init__(self):
        #self.mesh = None
        pass

    def makeMesh(self,femName,p3d_name, eStart):
        self.mesh = BDF()
        #cardsToInclude = set(['GRID','CQUAD4'])
        #mesh.setCardsToInclude(cardsToInclude)
        self.mesh.read_bdf(femName, xref=False)

        (sides,rSides) = self.makeConnections()

        edgeB = self.getCommonEdgeNumber(1,2)
        edgeC = self.getCommonEdgeNumber(1,51)
        edgeA = self.scaleEdge(edgeC+2)
        edgeD = self.scaleEdge(edgeB+2)
        print "edgeA = ",edgeA

        edge = self.findBoundEdge(eStart)
        #print "edge = ",edge

        edgeNum = self.getEdgeNumber(eStart,edge)
        eidCorner = copy.deepcopy(eStart)

        isBounded2 = True
        eidCorners = []
        savedIDs2  = []
        while isBounded2:
            eid = copy.deepcopy(eidCorner)
            if(eid in eidCorners):
                print "already found eidCorner=%s" %(eid)
                break
            eidCorners.append(eid)
            print "eid = ",eid

            savedIDs = []
            isBounded = True
            while isBounded:
                savedIDs.append(eid)
                (isBounded,leftoverEID,sharedEdge) = self.getBoundEdge(eid,edgeNum)
                #print "leftoverEID = ",leftoverEID
                eid = leftoverEID

            print ""
            print "eLast = ",savedIDs[-1]
            (eidCorner,sharedEdge) = self.turnCorner(eidCorner,edgeNum)
            savedIDs2.append(savedIDs)

            #eid = eidCorner
            #break
        #print "eidCorners = %s" %(eidCorners)

        self.makePlot3D(savedIDs2,eStart,p3d_name,edgeD)
        #print "e"

    def getNode(self,nid):
        node = self.mesh.nodes[nid]
        #spot = self.mesh.nodesmap[nid]
        #node = self.mesh.nodes[spot]
        return node


    def getElement(self,eid):
        element = self.mesh.elements[eid]
        return element

    def getCommonNode(self,edgeA,edgeB):
        setA = set(edgeA)
        setB = set(edgeB)
        intersect = list(setA.intersection(setB))
        nid = intersect[0]
        return nid

    def getEdgeNumberBasedOnNodes(self,eid,node1,node2):
        edge = [node1,node2]
        edge.sort()
        sides = sides[eid]
        edgeNum = sides.index(edge)
        return edgeNum

    def getEdgeOrder(self,eid,edgeA):
        edge1 = edgeA
        edge2 = self.scaleEdge(edgeA+1)
        edge3 = self.scaleEdge(edgeA+2)
        edge4 = self.scaleEdge(edgeA+3)
        return (edge1,edge2,edge3,edge4)

    def makePlot3D(self,elements,eStart,p3dname,edgeA):
        (n0,n1,n2,n3) = self.getEdgeOrder(eStart,edgeA)
        element = self.mesh.elements[eStart]

        print "n0=%s n1=%s" %(n0,n1)
        #node0 = element.get_node_id(n0)

        jmax = len(elements)+1
        imax = len(elements[0])+1
        kmax = 1
        out  = "1\n"
        out += "%s %s %s\n" %(imax,jmax,kmax)

        x=''; y=''; z=''; w=''
        for j,elementsLine in enumerate(elements):
            (x2,y2,z2,w2) = self.writeP3dLine(elementsLine,n0,n1) # 0,1
            x+=x2; y+=y2; z+=z2; w+=w2

        elementsLine = elements[-1]
        (x2,y2,z2,w2) = self.writeP3dLine(elementsLine,n3,n2)
        x+=x2; y+=y2; z+=z2; w+=w2

        #out += "%s\n%s\n%s\n%s\n" %(x,y,z,w)
        out += "%s\n%s\n%s\n" %(x,y,z)

        outfile = open(p3dname,'wb')
        outfile.write(out)
        outfile.close()

    def writeP3dLine(self,elementsLine,n0,n1):
        #mesh = self.mesh
        x=''; y=''; z=''; w=''
        print "eLine = ",elementsLine
        eSpot   = elementsLine[0]
        element = self.getElement(eSpot)
        #print "element = ",element
        print "n0 = ",n0
        node_ids = element.nodeIDs()
        nid0  = node_ids[n0]
        node0 = self.getNode(nid0)
        #print "node0 = ",node0
        print "nid = ",nid0
        x += '%f ' % node0.xyz[0]
        y += '%f ' % node0.xyz[1]
        z += '%f ' % node0.xyz[2]
        w+="%4s " % nid0

        for i,eid in enumerate(elementsLine):
            element = self.getElement(eid)
            node_ids = element.nodeIDs()
            print "n1 = ",n1
            nid1  = node_ids[n1]
            print "nid = ",nid1
            node1 = self.getNode(nid1)
            x += '%f ' % node1.xyz[0]
            y += '%f ' % node1.xyz[1]
            z += '%f ' % node1.xyz[2]
            w+="%4s " % nid1

            msg = 'nid0=%s i=%s nid0+i=%s nid1=%s' %(nid0,i,nid0+i,nid1)
            assert nid0+i+1==nid1,msg

        x+='\n'; y+='\n'; z+='\n'; w+='\n'
        return (x,y,z,w)

    def getCommonEdgeNumber(self,eidA,eidB):
        elementA = self.sides[eidA]

        setA = set(elementA)
        setB = set(self.sides[eidB])
        edgeList = list(setA.intersection(setB))
        edge = edgeList[0]

        edgeNum = elementA.index(edge)
        print "commonEdge = ",edge
        print "commonEdgeNum = ",edgeNum
        return edgeNum



    def turnCorner(self,eid,edgeNum):
        print "turnCorner eid=%s" %(eid)
        dEdgeNums = [1,2,3]
        for dEdgeNum in dEdgeNums:
            isBounded,leftoverEID,sharedEdge = self.getBoundEdge(eid,edgeNum+dEdgeNum)
            print "bound=%s eid=%s" %(isBounded,leftoverEID)
            if(isBounded):
                break
        assert isBounded
        print "leftoverEID = ",leftoverEID
        print "sharedEdge = ",sharedEdge
        return (leftoverEID,sharedEdge)


    def getBoundEdge(self,eid,edgeNum,):
        edgeN = self.getEdgeN(eid,edgeNum)
        eShared = self.rSides[edgeN]
        #print "eShared[%s] = %s" %(eid,eShared)
        nShared = len(eShared)-1

        isBounded = False
        leftoverEID = None
        if nShared>0:
            isBounded = True
            foundIndex = eShared.index(eid)
            leftoverEID = eShared[1-foundIndex]

        return isBounded,leftoverEID,edgeN

    def scaleEdge(self,edgeNum):
        if(edgeNum>=4):
            edgeNum-=4
        return edgeNum

    def getEdgeN(self,eid,edgeNum):
        edgeNum = self.scaleEdge(edgeNum)
        #msg = '%s=edgeNum<4' %(edgeNum)
        #print "edgeNum = ",edgeNum

        #assert edgeNum<4,msg
        edges = self.sides[eid]
        return edges[edgeNum]

    def findBoundEdge(self,eidStart):
        eSides = self.sides[eidStart]
        for side in eSides:
            neighboringElements = self.rSides[side]
            nEdges = len(neighboringElements)-1
            if(nEdges>0):
                break
        print "neighboringElements[%s]=%s" %(side,neighboringElements)
        print "  nEdges = %s" %(nEdges)
        print ""
        assert nEdges>0
        return side

    def getEdgeNumber(self,eid,edge):
        #element = getElement(eid)
        side = self.sides[eid]
        edgeNum = side.index(edge)
        print "side = ",side
        print "edge = ",edge
        print "edgeNum = ",edgeNum
        return edgeNum

    def makeList(self,nodes,i0,i1):
        node0 = nodes[i0]
        node1 = nodes[i1]
        nList = [node0,node1]
        return nList

    def makeTuple(self,nodes,i0,i1):
        nList = self.makeList(nodes,i0,i1)
        nList.sort()
        nList = tuple(nList)
        return nList

    def makeConnections(self):
        sides = {}
        reversedSides = {}
        for eid,element in self.mesh.elements.items():
            #print dir(element)
            nodes = element.nodeIDs()
            #print "element[%s]=%s" %(eid,nodes)
            nList = self.makeTuple(nodes,0,-1)
            sides[eid] = [nList]
            #print "adding nList = ",nList
            reversedSides.setdefault(nList,[]).append(eid)
            #reversedSides[nList].append(eid)
            for i in range(len(nodes)-1):
                nList = self.makeTuple(nodes,i,i+1)
                #print reversedSides
                #print "adding nList = ",nList
                sides[eid].append(nList)
                reversedSides.setdefault(nList,[]).append(eid)

            #print "sides[%s]=%s" %(eid,sides[eid])
            #if(eid>3):
            #    break
        print
        #for side,eid in reversedSides.items():
        #    print "rS[%s]=%s" %(side,eid)
        self.sides = sides
        self.rSides = reversedSides
        return (sides,reversedSides)

def run():
#    (X,Y) = makeCircle(1.,90.,45.)
#    (upper,lower) = readAirfoil('clarky.dat')

    eStart = 1
    mesh = BdfToP3d()
    p3d_name = 'mesh2.p3d'
    mesh.makeMesh('airfoil.bdf', p3d_name, eStart)




if __name__=='__main__':
    run()
