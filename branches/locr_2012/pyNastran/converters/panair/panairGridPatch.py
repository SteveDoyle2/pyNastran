import sys
from numpy import array, zeros, cross, transpose
from numpy.linalg import norm


def sInt(value):
    """
    int represented as a short float
    """
    value = "%f" % (value)
    return value.rstrip('0')


class PanairPatch(object):
    def __init__(self, iNetwork, netName, kt, cpNorm, x, y, z, log):
        self.log = log

        self.iNetwork = iNetwork  # lets it print out in order, does it need a deepcopy???
        self.netName = netName.strip()
        self.log.debug('netName=|%s|' % (netName))
        self.log.debug("****patch.netName=%s" % (self.netName))
        self.kt = kt
        self.cpNorm = cpNorm
        self.matchw = 0
        self.x = x
        self.y = y
        self.z = z
        shape = x.shape
        self.nRows = shape[0]
        self.nCols = shape[1]

        self.log.debug("shape = %s" % (str(shape)))

    def process(self):
        msg = '     network # being processed %3i\n\n' % (self.iNetwork + 1)
        return msg

    def quickSummary(self, cumPts, cumPn):
        msg = ''
        if self.kt == 1:
            src = 1
            dblt = 12
            nlopt1 = 5
            nropt1 = 3
            nlopt2 = 7
            nropt2 = -2
            ipot = 2
        elif self.kt == 5:
            src = 1
            dblt = 12
            nlopt1 = 6
            nropt1 = 9
            nlopt2 = 7
            nropt2 = -2
            ipot = 2
        elif self.kt == 18:
            src = 0
            dblt = 18
            nlopt1 = 0
            nropt1 = 9
            nlopt2 = 15
            nropt2 = 2
            ipot = 2
            if self.matchw == 1.:
                nlopt2 = 6
            self.log.debug("18...matchw = %s" % (self.matchw))
        elif self.kt == 20:
            src = 0
            dblt = 20
            nlopt1 = 0
            nropt1 = 9
            nlopt2 = 6
            nropt2 = 2
            ipot = 2
        else:
            raise NotImplementedError('new kt...kt=%s' % (self.kt))
        pts = self.nPoints()
        pans = self.nPanels()
        #cumPts=33
        #cumPn =50
        msg += ' %-10s %4s %7s %7s %3s ' % (self.netName,
                                            self.iNetwork + 1, self.nRows, self.nCols, self.kt,)
        msg += '%4s %5s %7s %7s %7s %7s ' % (
            src, dblt, nlopt1, nropt1, nlopt2, nropt2)
        msg += '%7s %7s %7s %7s ' % (ipot, pts, pans, self.cpNorm)
        msg += '%7s %7s\n' % (cumPts, cumPn)
        return msg

    def nPanels(self):
        return (self.nRows - 1) * (self.nCols - 1)

    def nPoints(self):
        return (self.nRows) * (self.nCols)

    def getPanelPoints(self, iPanel):
        r = iPanel % (self.nRows - 1)
        c = iPanel / (self.nRows - 1)

        #print "r=%s c=%s" %(r,c)
        p1 = self.getPoint(r, c)
        p2 = self.getPoint(r, c + 1)
        p3 = self.getPoint(r + 1, c + 1)
        p4 = self.getPoint(r + 1, c)
        return (p1, p2, p3, p4)

    def getPanelPointIDs(self, iPanel):
        r = iPanel % (self.nRows - 1)
        c = iPanel / (self.nRows - 1)

        #print "r=%s c=%s" %(r,c)
        p1 = self.getPointID(r, c)
        p2 = self.getPointID(r, c + 1)
        p3 = self.getPointID(r + 1, c + 1)
        p4 = self.getPointID(r + 1, c)
        return (p1, p2, p3, p4)

    def getSubpanelProperties(self, iPanel):
        (p1, p2, p3, p4) = self.getPanelPoints(iPanel)
        p5 = 0.5 * (p1 + p2)
        p6 = 0.5 * (p2 + p3)
        p7 = 0.5 * (p3 + p4)
        p8 = 0.5 * (p4 + p1)
        p9 = 0.25 * (p1 + p2 + p3 + p4)  # centroid

        p10 = 0.5 * (p5 + p6)
        p11 = 0.5 * (p6 + p7)
        p12 = 0.5 * (p7 + p8)
        p13 = 0.5 * (p8 + p5)
        N1 = cross(p10 - p12, p11 - p13)
        n1 = N1 / norm(N1)

        N2 = cross(p5 - p7, p6 - p8)
        n2 = N2 / norm(N2)

    def getPanelProperties(self, iPanel):
        (p1, p2, p3, p4) = self.getPanelPoints(iPanel)
        a = p1 - p3
        b = p2 - p4
        centroid = 0.25 * (p1 + p2 + p3 + p4)

        N = cross(a, b)
        normN = norm(N)
        n = N / normN  # normal vector

        S = 0.5 * normN  # area

        u = (p1 + p2 - p3 - p4) / 2.  # longitudinal
        p = (-p1 + p2 + p3 - p4) / 2.  # transverse

        u = 0.5 * (a + b)  # longitudinal vector in local coordinates
        p = 0.5 * (-a + b)  # transverse vector in local coordinates
        o = cross(n, u)  # normal to both vectors in local coordinates

        diameter = norm(a - b)

        return (S, n, centroid, diameter, u, p, o)

    def getPanelAreaNormal(self, iPanel):
        (p1, p2, p3, p4) = self.getPanelPoints(iPanel)
        a = p1 - p3
        b = p2 - p4

        N = cross(a, b)
        normN = norm(N)
        n = N / normN  # normal vector

        S = 0.5 * normN
        return (S, n)

    def getPanelArea(self, iPanel):
        # iPanel=200
        (p1, p2, p3, p4) = self.getPanelPoints(iPanel)

        a = p1 - p3
        b = p2 - p4
        S = 0.5 * norm(cross(a, b))
        return S

    def getPoint(self, row, col):
        return array([self.x[row][col],
                      self.y[row][col],
                      self.z[row][col]])

    def getPointID(self, row, col):
        return col * self.nRows + row

    def getIPoint(self, iPoint):
        iRow = iPoint / (self.nCols)
        iCol = iPoint % (self.nCols)
        #self.log.debug("iPoint=%s iRow=%s iCol=%s" %(iPoint,iRow,iCol))
        return self.getPoint(iRow, iCol)

    def getEdge(self, edgeNumber):
        """
        gets all the points associated with a given edge
		@code
                edge1
              0  1  2   -> i (row)
        edge4 3  4  5
              6  7  8  edge2
              9  10 11
            |   edge3
            j
        @endcode
        """
        edgeNumber = 2
        #edgeNumber = 4
        #print "edgeNumber=%s" %(edgeNumber)
        if edgeNumber == 1:
            x = self.x[0][:]
            y = self.y[0][:]
            z = self.z[0][:]  # pretty sure edge 1 is the 0th row
            p = [iCol for iCol in xrange(self.nCols)]  # good
        elif edgeNumber == 2:
            self.log.debug("x.shape = %s" % (str(self.x.shape)))
            x = self.x[:][self.nCols - 1]
            y = self.y[:][self.nCols - 1]
            z = self.z[:][self.nCols - 1]  # pretty sure edge 2 is the 0th row
            p = [iCol * (self.nRows) + (
                self.nRows - 1) for iCol in xrange(self.nCols)]
            #p = [iRow*(self.nCols)+(self.nCols-1) for iRow in xrange(self.nRows)]  #
        elif edgeNumber == 3:
            x = self.x[self.nRows - 1][:]
            y = self.y[self.nRows - 1][:]
            z = self.z[self.nRows - 1][:]  # pretty sure edge3 is the last row
            p = [iCol + self.nRows * iCol for iCol in xrange(
                self.nCols)]  # good
            #p = [(self.nCols-1)*(self.nRows)+iRow for iRow in xrange(self.nRows)]
        elif edgeNumber == 4:
            x = self.x[:][0]
            y = self.y[:][0]
            z = self.z[:][0]  # pretty sure edge 2 is the 0th row
            p = [self.nRows * iCol for iCol in xrange(self.nCols)]  # good
        else:
            raise ValueError('invalid edge; edgeNumber=%s' % (edgeNumber))
        self.log.debug("nRows=%s nCols=%s edgeNumber=%s" % (
            self.nRows, self.nCols, edgeNumber))
        #print "nx = ",len(x)
        #print "p = ",p
        #print "x = ",x
        #print "y = ",y
        #print "z = ",z
        p = [iPoint for iPoint in xrange(self.nPoints())]
        for pointID in p:
            #pointID = 2
            p2 = self.getIPoint(pointID)
            #print "point[%s]=%s" %(pointID,p2)

        return (p, x, y, z)

    def getEdges(self):
        nx = 2 * (self.nRows + self.nCols) - 2
        p = zeros(nx)
        x = zeros(nx)
        y = zeros(nx)
        z = zeros(nx)

        i = 0
        for edgeID in xrange(1, 4 + 1):
            (p1, x1, y1, z1) = self.getEdge(edgeID)
            nx1 = len(x1)
            p[i:i + nx1] = p1[0:nx1]
            x[i:i + nx1] = x1[0:nx1]
            y[i:i + nx1] = y1[0:nx1]
            z[i:i + nx1] = z1[0:nx1]
            self.log.debug("-----")
        return (p, x, y, z)

    def getElements(self, pointI):
        panels = []
        #print "pointI=%s" %(pointI)
        for iPanel in xrange(self.nPanels()):
            panel = self.getPanelPointIDs(iPanel)
            panel2 = []

            for p in panel:
                panel2.append(p + pointI)
            #print "panel=%s panel2=%s" %(str(panel),str(panel2))
            panels.append(panel2)
        return panels

    def getPoints(self):
        points = []
        #self.log.debug("size(X) = %s" %( str( self.x.shape ) ))
        #print "size(X) = %s" %( str(X.size())

        for j in xrange(self.nCols):
            for i in xrange(self.nRows):
                point = [self.x[i][j], self.y[i][j], self.z[i][j]]
                points.append(point)

        return points, len(points)

    def writeAsPlot3D(self):
        out = ''
        x = self.x.ravel()  # unravel
        y = self.y.ravel()  # unravel
        z = self.z.ravel()  # unravel

        for xi in x:
            out += "%s " % (xi)
        out += "\n"

        for yi in y:
            out += "%s " % (yi)
        out += "\n"

        for zi in z:
            out += "%s " % (zi)
        out += "\n"
        self.log.debug(out)
        #print x
        #for c in xrange(self.nCols):
        #    nPointsLeft = nFullLines*2+nPartialLines
        #    for r in xrange(0,self.nRows,2)
        return out

    def rotate(self):
        """
        not complete...
        """
        self.x = transpose(self.x)
        self.y = transpose(self.y)
        self.z = transpose(self.z)
        #self.x[0:n][:] = self.x[-n:-1][:] # something like this...

    def __repr__(self):
        """
        $points - body to wing wakes
        =kn                                               cpnorm
        1.
        =kt
        20.
        =nm       nn                                                          netname
        4.        2.                                                          awbw
        """
        #x = self.writeAsPlot3D()

        self.log.debug("*******")
        header = '$points - surface panels\n'
        points = ''

        header += '%-10s%-10s\n' % ('1.', self.cpNorm)  # nNetworks is 1
        header += '%-10s\n' % (sInt(self.kt))
        header += '%-10s%-10s%50s%-10s\n' % (
            sInt(self.nRows), sInt(self.nCols), '', self.netName)

        #nFullLines = nm/2
        #nPartialLines = nm%2
        #nLines = nFullLines+nPartialLines

        nFullLines = self.nRows / 2
        nPartialLines = self.nRows % 2
        nLines = nFullLines + nPartialLines

        for c in xrange(self.nCols):
            nPointsLeft = nFullLines * 2 + nPartialLines
            for r in xrange(0, self.nRows, 2):
                if nPointsLeft > 1:
                    x1 = self.x[r][c]
                    y1 = self.y[r][c]
                    z1 = self.z[r][c]

                    x2 = self.x[r + 1][c]
                    y2 = self.y[r + 1][c]
                    z2 = self.z[r + 1][c]
                    points += self.writePoints([x1, y1, z1], [x2, y2, z2])
                else:
                    x1 = self.x[r][c]
                    y1 = self.y[r][c]
                    z1 = self.z[r][c]
                    points += self.writePoint([x1, y1, z1])
                nPointsLeft -= 2
        return header + points

    def writePoints(self, point1, point2):
        point1 = self.fixPoint(point1)
        point2 = self.fixPoint(point2)

        out = "%-10s" * 6 % (point1[0], point1[1], point1[2],
                             point2[0], point2[1], point2[2])
        return out + '\n'

    def writePoint(self, point1):
        point1 = self.fixPoint(point1)
        out = "%-10s" * 3 % (point1[0], point1[1], point1[2])
        return out + '\n'

    def fixPoint(self, pointIn):
        pointOut = []
        for value in pointIn:
            sValue = '%s' % (value)
            if len(sValue) > 10:
                sValue = sValue[0:9]
            pointOut.append(sValue.rstrip('0'))
            #print "sValue=%s len=%s" %(sValue,len(sValue))
        #print "pointOut = ",pointOut
        return pointOut


class PanairWakePatch(PanairPatch):
    def __init__(self, iNetwork, netName, options, x, y, z, log):
        (kt, cpNorm, matchw, trailedPanel, edgeNumber, xWake, tWake) = options
        PanairPatch.__init__(self, iNetwork, netName, kt, cpNorm, x, y, z, log)

        self.log = log
        self.matchw = matchw
        self.trailedPanel = trailedPanel
        self.edgeNumber = edgeNumber
        self.xWake = xWake
        self.tWake = tWake
        self.log.debug("matchw = %s" % (self.matchw))
        self.log.debug("wake patch")

    def __repr__(self):
        header = '$trailing wakes\n'
        points = ''

        header += '%-10s%-10s\n' % (1., self.cpNorm)  # nNetworks is 1
        header += '%-10s%-10s\n' % (sInt(self.kt), self.matchw)
        #header += '%-10s%-10s%40s%-10s\n' %(self.nRows,self.nCols,'',self.netName)
        header += '%-10s%-10s%-10s%-10s%-30s%-10s\n' % (self.trailedPanel,
                                                        sInt(self.edgeNumber),
                                                        self.xWake,
                                                        sInt(self.tWake),
                                                        ' ',
                                                        self.netName)
        return header


class PanairGridHelper(object):

    def getCases(self, section):
        """
        $cases - no. of solutions
        =nacase
        1.
        """
        self.ncases = int(float(section[1][0:10]))
        self.log.debug("ncases = %s" % (self.ncases))
        return True

    def getMach(self, section):
        """
        $mach number
        =amach
        .6
        """
        self.mach = float(section[1][0:10])
        self.log.debug("mach = %s" % (self.mach))
        return True

    def setMach(self, mach):
        self.mach = mach

    def writeMach(self):
        out = '$mach number\n'
        out += '%-10s' % (self.mach) + '\n'
        return out

    def writeCases(self):
        out = '$cases - number of solutions\n'
        out += '%-10s' % (sInt(self.ncases)) + '\n'
        return out

    def getAlphas(self, section):
        """
        $angles-of-attack
        =alpc
        4.
        =alpha(1) alpha(2)  alpha(3)
        4.        10.       0.
        """
        self.alphas = []
        self.alphaC = float(section[1][0:10])  # alphaCompressibility
        sline = section[2].split()
        self.alphas = [float(slot) for slot in sline]
        self.log.debug("alphaC=%s alphas=%s" % (self.alphaC, self.alphas))
        return True

    def setAlphas(self, alphas, alphaC):
        self.alphaC = alphaC
        self.alphas = alphas
        self.ncases = len(alphas)

    def writeAlphas(self):
        out = '$angles-of-attack\n'
        out += '%-10s' % (self.alphaC) + '\n'
        out += '%-10s' * len(self.alphas) % (tuple(self.alphas)) + '\n'
        return out

    def getBetas(self, section):
        """
        $angles-of-attack
        =alpc
        4.
        =alpha(1) alpha(2)  alpha(3)
        4.        10.       0.
        """
        self.betas = []
        self.betaC = float(section[1][0:10])  # betaCompressibility
        sline = section[2].split()
        self.betas = [float(slot) for slot in sline]
        self.log.debug("betaC=%s betas=%s" % (self.betaC, self.betas))
        return True

    def setBetas(self, betas, betaC):
        self.betaC = betaC
        self.betas = betas
        self.ncases = len(betas)

    def writeBetas(self):
        out = '$yaw\n'
        out += '%-10s' % (self.betaC) + '\n'
        out += '%-10s' * len(self.betas) % (tuple(self.betas)) + '\n'
        return out

    def getReferenceQuantities(self, section):
        """
        $references for accumulated forces and moments
        =xref     yref      zref      nref
        46.       0.        0.
        =sref     bref      cref      dref
        2400.     60.       40.       90.
        """
        self.xref = float(section[1][0:10])  # 0
        self.yref = float(section[1][10:20])
        self.zref = float(section[1][20:30])
       #self.nref = float(section[1][30:40])
        self.sref = float(section[2][0:10])  # 0
        self.bref = float(section[2][10:20])
        self.cref = float(section[2][20:30])
        self.dref = float(section[2][30:40])
        self.log.debug("xref=%s yref=%s zref=%s" % (self.xref,
                                                    self.yref, self.zref))
        self.log.debug("sref=%s bref=%s cref=%s dref=%s " % (
            self.sref, self.bref, self.cref, self.dref))
        return True

    def writeReferenceQuantities(self):
        out = '$references for accumulated forces and moments\n'
        out += '%-10s' * 3 % (self.xref, self.yref, self.zref) + '\n'
        out += '%-10s' * 4 % (self.sref, self.bref, self.cref,
                              self.dref) + '\n'
        return out

    def getEnd(self, section):
        self.isEnd = True
        self.log.debug("end...")
        return True

    def writeEnd(self):
        if self.isEnd:
            return '$end of panair inputs\n '
