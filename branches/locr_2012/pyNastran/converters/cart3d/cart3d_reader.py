import os
import sys
from math import ceil, sqrt
from itertools import izip
from numpy import array

#from pyNastran.general.general import list_print
from struct import unpack
from pyNastran.op2.fortranFile import FortranFile
from pyNastran.general.general import is_binary


def convertToFloat(svalues):
    values = []
    for value in svalues:
        values.append(float(value))
    return values


class Cart3DAsciiReader(object):
    modelType = 'cart3d'
    isStructured = False
    isOutwardNormals = True

    def __init__(self, log=None, debug=False):
        self.isHalfModel = True
        self.cartType = None  # grid, result
        self.nPoints = None
        self.nElements = None
        self.infile = None
        self.infilename = None
        self.readHalf = False

        if log is None:
            from pyNastran.general.logger import dummyLogger
            if debug:
                word = 'debug'
            else:
                word = 'info'
            loggerObj = dummyLogger()
            log = loggerObj.startLog(word)  # or info
        self.log = log

    def readCart3d(self, infilename):
        """extracts the points, elements, and Cp"""
        self.infilename = infilename
        self.log.info("---starting reading cart3d file...|%s|---" %
                      (self.infilename))
        self.infile = open(self.infilename, 'r')
        self.readHeader()
        points = self.getPoints()
        elements = self.getElements(bypass=False)
        regions = self.getRegions(bypass=False)
        loads = self.readResults(0, self.infile)
        #if self.cartType=='results':

        #loads = self.getLoads(['Cp'])
        self.log.info("nPoints=%s  nElements=%s" % (self.nPoints,
                                                    self.nElements))
        self.log.info("---finished reading cart3d file...|%s|---" %
                      (self.infilename))
        return (points, elements, regions, loads)

    def getYZeroNodes(self, nodes):
        yZeroNodes = {}
        outerNodes = {}
        isInnerNode = {}

        #nodeTypes  = {}
        for (iNode, node) in sorted(nodes.items()):
            if node[1] == 0:
                yZeroNodes[iNode] = node
                isInnerNode[iNode] = True
            else:
                outerNodes[iNode] = [node[0], -node[1], node[2]]
                isInnerNode[iNode] = False

        return (yZeroNodes, outerNodes, isInnerNode)

    def makeFullModel(self, nodes, elements, regions, loads):
        """assumes a y=0 line of nodes"""
        self.log.info('---starting makeFullModel---')

        maxNode = len(nodes.keys())
        maxElements = len(elements.keys())
        (yZeroNodes, outerNodes, isInnerNode) = self.getYZeroNodes(nodes)

        #nodes2={}
        nodes2Map = {}
        iNode2 = 1
        for (iNode, node) in sorted(nodes.items()):
            if isInnerNode[iNode]:
                #pass
                #nodes2[iNode] = node  # dont need to set this...
                nodes2Map[iNode] = iNode
            else:
                #testElements2[maxNode+iNode2] = outerNodes[iNode]
                nodes[maxNode + iNode2] = outerNodes[iNode]
                nodes2Map[iNode] = maxNode + iNode2
                iNode2 += 1

        # based on the conversion (nodes2map) object, renumber each element
        # stick the result in the location for the new element
        #elementsB={}
        #regionsB={}

        nThree = 147120 + 10  # len(elements)*2   ## todo why is there +10???
        allSet = set([i + 1 for i in xrange(nThree)])
        elementKeysSet = set(elements.keys())
        missingKeys = list(allSet.difference(elementKeysSet))
        missingKeys.sort()
        #print "missingKeys = ",missingKeys
        #print "nMissing = ",len(missingKeys)
        iMissing = 0
        for (iElement, element) in sorted(elements.items()):
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
            regions[eMissing] = region
            iMissing += 1

        for i in xrange(1, len(elements)):
            assert i in elements, 'i=%s is not in elements...' % (i)
            #assert i in regions ,'i=%s is not in regions...' %(i)
        self.log.info('---finished makeFullModel---')
        sys.stdout.flush()

        return (nodes, elements, regions, loads)

    def makeHalfModel(self, nodes, elements, regions, Cp):
        nNodesStart = len(nodes)

        self.log.info('---starting makeHalfModel---')
        for (iNode, node) in sorted(nodes.items()):
            if node[1] < 0:
                #print iNode,'...',node
                del nodes[iNode]
                try:
                    del Cp[iNode]  # assume loads=0
                except:
                    pass

        sys.stdout.flush()

        for (iElement, element) in elements.items():
            #print "iElement = |%s|" %(iElement),type(iElement)
            inNodeList = [True, True, True]
            for (i, node) in enumerate(element):
                if node not in nodes:
                    inNodeList[i] = False

            if False in inNodeList:
                #print iElement,element,inNodeList,type(iElement)
                del elements[iElement]
                del regions[iElement]

        assert len(nodes) != nNodesStart
        self.log.info('---finished makeHalfModel---')
        return (nodes, elements, regions, Cp)

    def makeMirrorModel(self, nodes, elements, regions, loads):
        #nNodesStart = len(nodes)

        self.log.info('---starting makeMirrorModel---')
        for (iNode, node) in sorted(nodes.items()):
            if node[1] < 0:
                #print iNode,'...',node
                del nodes[iNode]
                #del Cp[iNode]  # assume loads=0
        sys.stdout.flush()

        for (iElement, element) in elements.items():
            #print "iElement = |%s|" %(iElement),type(iElement)
            inNodeList = [True, True, True]
            for (i, node) in enumerate(element):
                if node not in nodes:
                    inNodeList[i] = False

            if False in inNodeList:
                #print iElement,element,inNodeList,type(iElement)
                del elements[iElement]
                del regions[iElement]

        #assert len(nodes)!=nNodesStart
        self.log.info('---finished makeMirrorModel---')

        return self.makeFullModel(nodes, elements, regions, loads)
        #return self.renumberMesh(nodes,elements,regions,loads)

    def renumberMesh(self, nodes, elements, regions, loads):
        iNodeCounter = 1
        iElementCounter = 1

        NodeOldToNew = {}
        #ElementOldToNew = {}
        nodes2 = {}
        elements2 = {}
        regions2 = {}
        for (iNode, node) in sorted(nodes.items()):
            NodeOldToNew[iNode] = iNodeCounter
            nodes2[iNodeCounter] = node
            iNodeCounter += 1

        for (iElement, element) in sorted(elements.items()):
            element2 = []
            for nID in element:
                nID2 = NodeOldToNew[nID]
                element2.append(nID2)
            elements2[iElementCounter] = element2
            regions2[iElementCounter] = regions[iElement]
            iElementCounter += 1

        return (nodes, elements, regions, loads)

    def writeOutfile(self, outfilename, points, elements, regions):
        assert len(points) > 0
        self.log.info("---writing cart3d file...|%s|---" % (outfilename))
        f = open(outfilename, 'wb')

        self.writeHeader(f, points, elements)
        self.writePoints(f, points)
        self.writeElements(f, points, elements)
        self.writeRegions(f, regions)
        f.close()

    def writeRegions(self, f, regions):
        msg = ''
        for (iElement, region) in sorted(regions.items()):
            msg += "%s\n" % (region)
        f.write(msg)

    def writeElements(self, f, nodes, elements):
        msg = ''
        strElements = len(str(len(elements)))
        Format = " %%%ss " % (strElements)
        #Format = '%-10s '+Format*3 + '\n'
        #Format = Format*3 + '\n'
        #print "Format = ",Format

        maxNodes = []
        for (iElement, element) in sorted(elements.items()):
            #self.log.info("element = %s" %(element))
            #msg += Format %(iElement,element[0],element[1],element[2])
            msg += "%s  %s  %s\n" % (element[0], element[1], element[2])
            for nid in element:
                assert nid in nodes, 'nid=%s is missing' % (nid)

            maxNode = max(element)
            maxNodes.append(maxNode)
        self.log.info("maxNodeID = %s" % (max(maxNodes)))
        f.write(msg)

    def writePoints(self, f, points):
        msg = ''
        for (iPoint, point) in sorted(points.items()):
            #msg += "%-10s %-15s %-15s %-15s\n" %(iPoint,point[0],point[1],point[2])
            msg += "%6.6f  %6.6f  %6.6f\n" % (point[0], point[1], point[2])
            if iPoint % 1000:
                f.write(msg)
                msg = ''
        f.write(msg)

    def writeHeader(self, f, points, elements):
        nPoints = len(points)
        nElements = len(elements)
        msg = "%s %s\n" % (nPoints, nElements)
        f.write(msg)

    def readHeader(self):
        line = self.infile.readline()
        sline = line.strip().split()
        if len(sline) == 2:
            self.nPoints, self.nElements = int(sline[0]), int(sline[1])
            self.cartType = 'grid'
        elif len(sline) == 3:
            self.nPoints, self.nElements, self.nResults = int(
                sline[0]), int(sline[1]), int(sline[2])
            self.cartType = 'result'
        else:
            raise ValueError('invalid result type')

        if self.readHalf:
            self.nElementsRead = self.nElements // 2
            self.nElementsSkip = self.nElements // 2
        else:
            self.nElementsRead = self.nElements
            self.nElementsSkip = 0

    def getPoints(self):
        """
        A point is defined by x,y,z and the ID is the location in points.
        """
        points = {}
        p = 0
        data = []
        while p < self.nPoints:
            data += self.infile.readline().strip().split()
            while len(data) > 2:
                x = data.pop(0)
                y = data.pop(0)
                z = data.pop(0)
                points[p + 1] = array([float(x), float(y), float(z)])
                p += 1

        maxX = self.getMax(points, 0)
        maxY = self.getMax(points, 1)
        maxZ = self.getMax(points, 2)

        minX = self.getMin(points, 0)
        minY = self.getMin(points, 1)
        minZ = self.getMin(points, 2)

        self.log.info("X  max=%g min=%g" % (maxX, minX))
        self.log.info("Y  max=%g min=%g" % (maxY, minY))
        self.log.info("Z  max=%g min=%g" % (maxZ, minZ))
        return points

    def getMin(self, points, i):
        minX = points[1][i]
        for key, point in points.items():
            minX = min(minX, point[i])
        return minX

    def getMax(self, points, i):
        maxX = points[1][i]
        for key, point in points.items():
            maxX = max(maxX, point[i])
        return maxX

    def getElements(self, bypass=False):
        """
        An element is defined by n1,n2,n3 and the ID is the location in elements.
        """
        assert bypass == False
        elements = {}

        print "nElementsRead=%s nElementsSkip=%s" % (
            self.nElementsRead, self.nElementsSkip)

        e = 0
        data = []
        while e < self.nElementsRead:
            data += self.infile.readline().strip().split()
            while len(data) > 2:
                n1 = int(data.pop(0))
                n2 = int(data.pop(0))
                n3 = int(data.pop(0))
                elements[e + 1] = [n1, n2, n3]
                e += 1

        e = 0
        while e < self.nElementsSkip:
            data += self.infile.readline().strip().split()
            while len(data) > 2:
                data.pop()  # popping from the end is faster
                data.pop()
                data.pop()
                e += 1

        return elements

    def getRegions(self, bypass=True):
        regions = {}
        if bypass:
            for i in xrange(self.nElements):
                self.infile.readline()
        else:
            r = 0
            data = []
            while r < self.nElementsRead:
                data += self.infile.readline().strip().split()
                while len(data) > 0:
                    regions[r + 1] = int(data.pop(0))
                    r += 1

            r = 0
            while r < self.nElementsSkip:
                data += self.infile.readline().strip().split()
                while len(data) > 0:
                    data.pop()
                    r += 1

        return regions

    def readResults(self, i, infile):
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
        if(self.cartType == 'result'):
            nResultLines = int(ceil(self.nResults / 5.)) - 1
            #print "nResultLines = ",nResultLines
            for pointNum in xrange(self.nPoints):
                #print "pointNum = ", pointNum
                # rho rhoU,rhoV,rhoW,pressure/rhoE/E
                sline = infile.readline().strip().split()
                i += 1
                for n in xrange(nResultLines):
                    #print "nResultLines = ",n
                    sline += infile.readline().strip().split()  # Cp
                    #infile.readline()
                    i += 1
                    #print "sline = ",sline
                    #gamma = 1.4
                    #else:
                    #    p=0.
                sline = self.getList(sline)

                #p=0
                cp = sline[0]
                rho = float(sline[1])
                if(rho > abs(0.000001)):
                    rhoU = float(sline[2])
                    rhoV = float(sline[3])
                    rhoW = float(sline[4])
                    #rhoE = float(sline[5])
                    VV = (rhoU) ** 2 + (rhoV) ** 2 + (rhoW) ** 2 / rho ** 2
                    mach = sqrt(VV)
                    if mach > 10:
                        print "nid=%s Cp=%s mach=%s rho=%s rhoU=%s rhoV=%s rhoW=%s" % (pointNum, cp, mach, rho, rhoU, rhoV, rhoW)

                Cp[pointNum + 1] = cp
                Mach[pointNum + 1] = mach
                U[pointNum + 1] = rhoU / rho
                V[pointNum + 1] = rhoV / rho
                W[pointNum + 1] = rhoW / rho
                #print "pt=%s i=%s Cp=%s p=%s" %(pointNum,i,sline[0],p)
            ###
        ###
        loads = {}
        if Cp:
            loads['Cp'] = Cp
        if Mach:
            loads['Mach'] = Mach
        if U:
            loads['U'] = U
            loads['V'] = V
            loads['W'] = W
        return loads
    ###

    def getList(self, sline):
        """Takes a list of strings and converts them to floats."""
        try:
            sline2 = convertToFloat(sline)
        except ValueError:
            print("sline = %s" % (sline))
            raise SyntaxError('cannot parse %s' % (sline))
        return sline2

    def getValue(self, sline):
        """Converts a string to a float."""
        value = float(sline[1])
        return value
    ###

    def exportToNastran(self, fname, points, elements, regions):
        from pyNastran.bdf.fieldWriter import printCard

        f = open(fname, 'wb')
        for nid, grid in enumerate(1, points):
            (x, y, z) = grid
            f.write(printCard(['GRID', nid, '', x, y, z]))
        ###

        e = 1e7
        nu = 0.3
        g = ''
        thickness = 0.25
        setRegions = list(set(regions))
        for pidMid in setRegions:
            f.write(printCard(['MAT1', pidMid, e, g, nu]))
            f.write(printCard(['PSHELL', pidMid, pidMid, thickness]))
        ###

        for eid, (nodes, region) in enumerate(1, izip(elements, regions)):
            (n1, n2, n3) = nodes
            f.write(printCard(['CTRIA3', eid, region, n1, n2, n3]))
        ###
        f.close()

#------------------------------------------------------------------
    def getSegments(self, nodes, elements):
        segments = {}  # key=eid,
        lengths = {}
        for eid, e in elements.iteritems():
            a = tuple(sorted([e[0], e[1]]))  # segments of e
            b = tuple(sorted([e[1], e[2]]))
            c = tuple(sorted([e[2], e[0]]))
            print eid, a, b, c

            #a.sort()
            #b.sort()
            #c.sort()
            if a in segments:
                segments[a].append(eid)
            else:
                segments[a] = [eid]
                #print(a)
                #print nodes[1]
                lengths[a] = nodes[a[1]] - nodes[a[0]]

            if b in segments:
                segments[b].append(eid)
            else:
                segments[b] = [eid]
                lengths[b] = nodes[b[1]] - nodes[b[0]]

            if c in segments:
                segments[c].append(eid)
            else:
                segments[c] = [eid]
                lengths[c] = nodes[c[1]] - nodes[c[0]]

        return segments, lengths

    def makeQuads(self, nodes, elements):
        segments, lengths = self.getSegments(nodes, elements)
        for eid, e in elements.iteritems():
            a = tuple(sorted([e[0], e[1]]))  # segments of e
            b = tuple(sorted([e[1], e[2]]))
            c = tuple(sorted([e[2], e[0]]))
            #a.sort()
            #b.sort()
            #c.sort()
            print eid, e
            print segments[a]
            print lengths[a]
            print len(segments[a])

            eidA = self.getSegment(a, eid, segments)
            eidB = self.getSegment(b, eid, segments)
            eidC = self.getSegment(c, eid, segments)
            print "eidA=%s eidB=%s eidC=%s" % (eidA, eidB, eidC)
            if eidA:
                i = 0
                e2 = elements[eidA]
                self.checkQuad(nodes, eid, eidA, e, e2, a, b, c, i)
                del segments[a]
            if eidB:
                i = 1
                e2 = elements[eidB]
                self.checkQuad(nodes, eid, eidB, e, e2, a, b, c, i)
                del segments[b]
            if eidC:
                i = 2
                e2 = elements[eidC]
                self.checkQuad(nodes, eid, eidC, e, e2, a, b, c, i)
                del segments[c]

            print "------"
            #break
        #for segment in segments:
        asdf

    def checkQuad(self, nodes, eid, eidA, e, e2, a, b, c, i):
        """
        @code
        A----B
        | \ e|
        |e2 \|
        C----D
		
        two tests
           1.  folding angle A-B x A-C
           2a. abs(A-C) - abs(B-D)  = 0  (abs to prevent 2L)
           2b. abs(A-B) - abs(C-D)  = 0
        @endcode
        """

        iplus1 = i + 1
        iplus2 = i + 2

        if iplus1 > 2:
            iplus1 -= 3
        if iplus2 > 2:
            iplus2 -= 3
        print (i, iplus1)
        print (iplus1, iplus2)
        print (iplus2, i)
        AD = nodes[e[i]] - nodes[e[iplus1]]
        AB = nodes[e[iplus1]] - nodes[e[iplus2]]
        BD = nodes[e[iplus2]] - nodes[e[i]]

        print AD
        print AB
        print BD
        print e2
        j = e2.index(e[i])

        jplus1 = j + 1
        jplus2 = j + 2
        if jplus1 > 2:
            jplus1 -= 3
        if jplus2 > 2:
            jplus2 -= 3

        print "DA = ", e[j], e[jplus1]
        DA = nodes[e[j]] - nodes[e[jplus1]]
        print DA

        asdf

    def getSegment(self, a, eid, segments):
        if a in segments:
            aElems = segments[a]
            print aElems
            i = aElems.index(eid)
            #print i
            aElems.pop(i)
            #print aElems
            eidA = aElems[0]
            #eidA = elements[a]
            print "eidA = ", eidA
            return eidA
        return None


#------------------------------------------------------------------
class Cart3DBinaryReader(FortranFile, Cart3DAsciiReader):
    modelType = 'cart3d'
    isStructured = False
    isOutwardNormals = True

    def __init__(self, log=None, debug=False):
        Cart3DAsciiReader.__init__(self)
        FortranFile.__init__(self)

        self.isHalfModel = True
        self.cartType = None  # grid, result
        self.nPoints = None
        self.nElements = None
        self.infile = None
        self.infilename = None
        self.readHalf = False
        self.isNodeArray = True

        self.makeOp2Debug = False
        self.n = 0

        if log is None:
            from pyNastran.general.logger import dummyLogger
            loggerObj = dummyLogger()
            if debug:
                log = loggerObj.startLog('debug')  # or info
            else:
                log = loggerObj.startLog('info')  # or info
        self.log = log

    def readCart3d(self, infileName):
        self.infilename = infileName
        self.op2 = open(infileName, 'rb')
        (self.nPoints, self.nElements) = self.readHeader()
        points = self.readNodes(self.nPoints)
        elements = self.readElements(self.nElements)
        regions = self.readRegions(self.nElements)
        loads = {}
        #print "size = ",size
        #print self.printSection2(200,'>')

        #print self.printSection2(200,'>')
        return (points, elements, regions, loads)

    def readHeader(self):
        data = self.op2.read(4)
        size, = unpack('>i', data)
        #print "size = ",size

        data = self.op2.read(size)
        so4 = size // 4  # size over 4
        if so4 == 3:
            (nPoints, nElements, nResults) = unpack('>iii', data)
            print "nPoints=%s nElements=%s nResults=%s" % (
                nPoints, nElements, nResults)
        elif so4 == 2:
            (nPoints, nElements) = unpack('>ii', data)
            print "nPoints=%s nElements=%s" % (nPoints, nElements)
        else:
            raise RuntimeError('in the wrong spot...endian...')
        self.op2.read(8)  # end of first block, start of second block

        return (nPoints, nElements)

    def readNodes(self, nPoints):
        #print "starting Nodes"
        isBuffered = True
        size = nPoints * 12  # 12=3*4 all the points
        Format = '>' + 'f' * 3000  # 3000 floats; 1000 points

        nodes = {}

        np = 1
        while size > 12000:  # 12k = 4 bytes/float*3 floats/point*1000 points
            #print "size = ",size
            n = np + 999
            #print "nStart=%s np=%s" %(n,np)
            data = self.op2.read(4 * 3000)
            nodeXYZs = list(unpack(Format, data))
            while nodeXYZs:
                z = nodeXYZs.pop()
                y = nodeXYZs.pop()
                x = nodeXYZs.pop()
                node = array([x, y, z])
                assert n not in nodes, 'nid=%s in nodes' % (n)
                nodes[n] = node
                n -= 1
                np += 1
            size -= 4 * 3000

        assert size >= 0
        #print "size = ",size
        #while size>0: # 4k is 1000 points
        n = nPoints
        if size > 0:
            #leftover = size-(nPoints-np)*12
            data = self.op2.read(size)
            #print "leftover=%s size//4=%s" %(leftover,size//4)
            Format = '>' + 'f' * (size // 4)

            #print "len(data) = ",len(data)
            nodeXYZs = list(unpack(Format, data))
            while nodeXYZs:
                z = nodeXYZs.pop()
                y = nodeXYZs.pop()
                x = nodeXYZs.pop()
                node = array([x, y, z])
                assert n not in nodes, 'nid=%s in nodes' % (n)
                nodes[n] = node
                n -= 1
                np += 1
            ###
            size = 0

        #for p,point in sorted(nodes.iteritems()):
            #print "%s %s %s" %(tuple(point))

        if isBuffered:
            pass
        else:
            raise RuntimeError('unBuffered')

        for nid in xrange(1, nPoints + 1):
            assert nid in nodes, 'nid=%s not in nodes' % (nid)
        self.op2.read(8)  # end of second block, start of third block
        return nodes

    def readElements(self, nElements):
        self.nElementsRead = nElements
        self.nElementsSkip = 0
        #print "starting Elements"
        isBuffered = True
        size = nElements * 12  # 12=3*4 all the elements
        Format = '>' + 'i' * 3000

        elements = {}

        ne = 0
        while size > 12000:  # 4k is 1000 elements
            #print "size = ",size
            e = ne + 1000
            data = self.op2.read(4 * 3000)
            nodes = list(unpack(Format, data))
            while nodes:
                n3 = nodes.pop()
                n2 = nodes.pop()
                n1 = nodes.pop()
                element = [n1, n2, n3]
                elements[e] = element
                e -= 1
                ne += 1
            size -= 4 * 3000

        assert size >= 0
        #print "size = ",size
        #while size>0: # 4k is 1000 elements
        if size > 0:
            #leftover = size-(nElements-ne)*12
            data = self.op2.read(size)
            #print "leftover=%s size//4=%s" %(leftover,size//4)
            Format = '>' + 'i' * (size // 4)

            #print "len(data) = ",len(data)
            nodes = list(unpack(Format, data))
            e = nElements
            while nodes:
                n3 = nodes.pop()
                n2 = nodes.pop()
                n1 = nodes.pop()
                element = [n1, n2, n3]
                elements[e] = element
                e -= 1
                ne += 1
            ###
            size = 0

        #for p,point in sorted(nodes.iteritems()):
            #print "%s %s %s" %(tuple(point))

        if isBuffered:
            pass
        else:
            raise RuntimeError('unBuffered')

        self.op2.read(8)  # end of third (element) block, start of regions (fourth) block
        #print "finished Elements"
        return elements

    def readRegions(self, nElements):
        #print "starting Regions"
        #isBuffered = True
        size = nElements * 4  # 12=3*4 all the elements
        Format = '>' + 'i' * 3000

        regions = {}

        nr = 0
        while size > 12000:  # 12k is 3000 elements
            #print "size = ",size
            data = self.op2.read(4 * 3000)
            regionData = list(unpack(Format, data))

            r = nr + 3000
            #print "len(regionData) = ",len(regionData)
            while regionData:
                regions[r] = regionData.pop()
                r -= 1
                nr += 1
            size -= 4 * 3000
            #print "size = ",size

        assert size >= 0
        #print "size = ",size

        if size > 0:
            #leftover = size-(nElements-nr)*4
            data = self.op2.read(size)
            #print "leftover=%s size//4=%s" %(leftover,size//4)
            Format = '>' + 'i' * (size // 4)

            #print "len(data) = ",len(data)
            regionData = list(unpack(Format, data))
            r = nElements
            while regionData:
                regions[r] = regionData.pop()
                r -= 1
                nr += 1
            ###
            size = 0

        self.op2.read(4)  # end of regions (fourth) block
        #print "finished Regions"
        return regions

#------------------------------------------------------------------


def genericCart3DReader(infileName, log=None, debug=False):
    print "infileName = ", infileName
    f = open(infileName, 'rb')
    data = f.read(4)
    f.close()

    if is_binary(infileName):

    #eight, = unpack('>i',data)
    #if eight==8:  # is binary
        #print "Binary eight = ",eight
        obj = Cart3DBinaryReader(log, debug)
    else:
        #print "Ascii eight = ",eight
        obj = Cart3DAsciiReader(log, debug)

    return obj


if __name__ == '__main__':
    # binary
    cart3dGeom = os.path.join('models', 'spike.a.tri')
    #outfilename = os.path.join('models','threePlugs2.tri')
    cart = genericCart3DReader(cart3dGeom)
    (points, elements, regions, loads) = cart.readCart3d(cart3dGeom)
    #cart.makeQuads(points,elements)

    cart.writeOutfile(outfilename, points, elements, regions)

    # ascii
    cart3dGeom = os.path.join('models', 'threePlugs.a.tri')
    outfilename = os.path.join('models', 'threePlugs2.a.tri')
    cart2 = genericCart3DReader(cart3dGeom)
    (points, elements, regions, loads) = cart2.readCart3d(cart3dGeom)
    cart2.writeOutfile(outfilename, points, elements, regions)
    #print points

if 0:
    basepath = os.getcwd()
    configpath = os.path.join(basepath, 'inputs')
    workpath = os.path.join(basepath, 'outputsFinal')
    #bdfModel   = os.path.join(configpath,'aeroModel.bdf')
    #assert os.path.exists(bdfModel),'|%s| doesnt exist' %(bdfModel)
    os.chdir(workpath)
    print "basepath", basepath

    cart3dGeom = os.path.join(configpath, 'Cart3d_bwb2.i.tri')
    cart3dGeom2 = os.path.join(workpath, 'Cart3d_half.i.tri')
    cart3dGeom3 = os.path.join(workpath, 'Cart3d_full.i.tri')
    #cart3dGeom4 = os.path.join(workpath,'Cart3d_full2.i.tri')

    cart = Cart3DAsciiReader()
    (points, elements, regions, loads) = cart.readCart3d(cart3dGeom)
    (points, elements, regions, loads) = cart.makeHalfModel(points,
                                                            elements, regions, loads)
    cart.writeOutfile(cart3dGeom2, points, elements, regions)

    #cart = Cart3DAsciiReader(cart3dGeom)  # bJet.a.tri
    #(cartPoints,elements,regions,loads) = cart.read()
    #points2 = {}
    #for (iPoint,point) in sorted(points.items()):
        #(x,y,z) = point
        #print "p=%s x=%s y=%s z=%s  z2=%s" %(iPoint,x,y,z,z+x/10.)
        #points2[iPoint] = [x,y,z+x/10.]
    #(points,elements,regions,loads) = cart.makeMirrorModel(points2,elements,regions,loads)

    cart2 = Cart3DAsciiReader(cart3dGeom2)
    (points, elements, regions, loads) = cart2.read()

    #makeFullModel
    (points, elements, regions, loads) = cart2.makeFullModel(
        points, elements, regions, loads)  # half model values going in
    cart2.writeOutfile(cart3dGeom3, points, elements, regions)

    #cart3 = Cart3DAsciiReader(cart3dGeom2)
    #(points,elements,regions,loads) = cart3.read()
    #(points,elements,regions,loads) = cart3.makeMirrorModel(points,elements,regions,loads)


    #print "loads = ",list_print(loads),len(loads)
    #cartOutfile = os.path.join(workpath,'bJet.a.tri_test')
    #cart.writeInfile(cartOutfile,cartPoints,elements,regions)
