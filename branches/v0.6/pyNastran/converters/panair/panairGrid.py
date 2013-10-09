import os
#import sys
#import copy
from itertools import izip, count

from math import ceil, sin, cos, radians
from numpy import array, zeros, ones


from panairGridPatch import PanairPatch, PanairWakePatch
from panairWrite import PanairWrite
from pyNastran.utils.log import get_logger
from pyNastran.utils import print_bad_path

#from pyNastran.utils import list_print

#CL = -Fx*sin(alpha)*cos(beta) + Fy*sin(alpha)*sin(beta) +Fz*cos(alpha)
#CD =  Fx*cos(alpha)*cos(beta) - Fy*cos(alpha)*sin(beta) +Fz*sin(alpha)
#CY =  Fx*sin(beta) +Fy*cos(beta)


def fortran_value(value):
    return "%8.4E" % value


class PanairGrid(object):
    modelType = 'panair'

    def __init__(self, infileName, log=None, debug=True):
        assert os.path.exists(infileName), print_bad_path(infileName)

        self.infileName = infileName
        self.nNetworks = 0
        self.patches = {}
        self.lines = []

        self.alphas = [0.]
        self.ncases = None
        self.betas = [0.]
        self.alphaC = 0.
        self.betaC = 0.

        self.sref = 1.
        self.bref = 1.
        self.cref = 1.
        self.dref = 1.

        self.xref = 0.
        self.yref = 0.
        self.zref = 0.

        self.isings = 0.0
        self.isingp = 0.0
        self.igeomp = 0.0
        self.icontp = 0.0
        self.ibconp = 0.0
        self.iedgep = 0.0
        self.ipraic = 0.0
        self.nexdgn = 0.0
        self.ioutpr = 0.0
        self.ifmcpr = 0.0
        self.icostp = 0.0

        self.isEnd = None

        self.peaSection = ''


        self.mach = 0.0
        self.dataCheck = 2
        self.titleSection = ''
        self.xyzSection = ''
        self.streamlineSection = ''
        self.flowSection = ''
        self.sectionalPropSection = ''
        self.gridSection = ''

        self.msg = ''

        self.log = get_logger(log, 'debug' if debug else 'info')

    def print_file(self):
        msg = ''
        for i, line in enumerate(self.lines):
            msg += "%6s %s" % (i + 1, line)
        msg += '                                      record of input processing\n\n\n'
        return msg

    def get_panel_points(self, iPanel):
        r = iPanel % (self.nRows - 1)
        c = iPanel / (self.nRows - 1)

        #print "r=%s c=%s" %(r,c)
        p1 = self.get_point(r, c)
        p2 = self.get_point(r, c + 1)
        p3 = self.get_point(r + 1, c + 1)
        p4 = self.get_point(r + 1, c)
        return (p1, p2, p3, p4)

    def get_panel_point_IDs(self, iPanel):
        r = iPanel % (self.nRows - 1)
        c = iPanel / (self.nRows - 1)

        #print "r=%s c=%s" %(r,c)
        p1 = self.get_point_ID(r, c)
        p2 = self.get_point_ID(r, c + 1)
        p3 = self.get_point_ID(r + 1, c + 1)
        p4 = self.get_point_ID(r + 1, c)
        return (p1, p2, p3, p4)

    def nPanels(self):
        totalNPanels = 0
        for patchID in xrange(self.nPatches()):
            patch = self.get_patch(patchID)
            if patch.kt == 1:
                totalNPanels += patch.nPanels()
                #self.log.debug("nPanels = %s" % (totalNPanels))
        return totalNPanels

    def nPatches(self):
        return len(self.patches)

    def get_patch(self, ID):
        return self.patches[ID]

    def update_cases(self):
        """
        reduces confusion by only printing cases that will run
        """
        if len(self.alphas) > self.ncases:
            print("self.ncases =", self.ncases)
            self.alphas = self.alphas[:self.ncases]
        if len(self.betas) > self.ncases:
            self.betas = self.alphas[:self.ncases]

    def write_panair(self, outfileName):
        self.update_cases()
        outfile = open(outfileName, 'wb')
        outfile.write(self.titleSection)
        outfile.write(self.write_data_check())
        outfile.write(self.symmetrySection)

        outfile.write(self.write_mach())
        outfile.write(self.write_cases())
        outfile.write(self.write_alphas())
        outfile.write(self.write_reference_quantities())

        #outfile.write(self.alphaSection)
        #outfile.write(self.caseSection)
        for patchName, patch in sorted(self.patches.iteritems()):
            outfile.write(str(patch))

        outfile.write(self.xyzSection)
        outfile.write(self.streamlineSection)
        outfile.write(self.flowSection)
        outfile.write(self.sectionalPropSection)
        outfile.write(self.gridSection)
        outfile.write(self.write_printout())

        outfile.write(self.peaSection)
        outfile.write(self.write_liberalized_abutments())
        outfile.write(self.write_end())

        outfile.close()

    def remove_comments(self, lines):
        lines2 = []
        for line in lines:
            line = line.rstrip().lower()
            if '=' in line:
                #print "line -> |%s|" %(line)
                if '=' is not line[0]:
                    self.log.debug("line[0] -> %s" % (line[0]))
                    line = line.split('=')[0]
                    #self.log.debug("******")
                    lines2.append(line)
                # else: skip
            else:
                lines2.append(line)

        return lines2

    def _read_title(self, section):
        #print "hi"
        #self.title = section[1:]
        self.titleSection = '\n'.join(section) + '\n'
        self.titleLines = section[1:]
        return True

    def _read_data_check(self, section):
        self.dataCheck = int(float(section[1][0:10]))
        #self.dataCheck = 2
        return True

    def _read_symmetry(self, section):
        """
        @code
        $symmetry - xz plane of symmetry
        =misymm   mjsymm
        1.        0.
        @endcode

        @warning
            doesnt consider antisymmetryic
        """
        # doesnt consider antisymmetric
        self.XZsymmetry = int(float(section[1][0:10]))
        self.XYsymmetry = int(float(section[1][10:20]))   # ???

        # doesnt consider antisymmetric        self.XYsymmetry = int(float(section[1][10:20]))
        self.nSymmetryPlanes = self.XZsymmetry + self.XYsymmetry
        self.symmetrySection = '\n'.join(section) + '\n'
        self.log.debug("XZsymmetry=%s XYsymmetry=%s" % (
            self.XZsymmetry, self.XYsymmetry))
        return True

    def split_points(self, lines, nActual, nRemainder):
        """
        reads the points
        """
        points = []
        for n in xrange(nActual):
            #(x1,y1,z1) = lines[n][0 :10].strip(),lines[n][10:20].strip(),lines[n][20:30].strip()
            #print "x1=%s y1=%s z1=%s" %(x1,y1,z1)
            line = lines[n]
            (x1, y1, z1) = float(line[0:10]), float(line[10:20]), float(line[20:30])
            (x2, y2, z2) = float(line[30:40]), float(lines[n][40:50]), float(line[50:60])
            point1 = array([x1, y1, z1])
            point2 = array([x2, y2, z2])
            #point1 = [x1,y1,z1]
            #point2 = [x2,y2,z2]
            #print list_print(point1)
            #print list_print(point2)
            points.append(point1)
            points.append(point2)

        if nRemainder:
            n += 1
            #print "***"
            line = lines[n]
            (x1, y1, z1) = float(line[0:10]), float(line[10:20]), float(line[20:30])
            point1 = array([x1, y1, z1])
            #print list_print(point1)
            points.append(point1)
        #print "points = ",list_print(points)
        return points

    def add_wake_patch(self, netName, options, x, y, z):
        #print "self.nNetworks = ",self.nNetworks
        patch = PanairWakePatch(self.nNetworks, netName, options, x, y, z, self.log)
        self.msg += patch.process()
        #print "get_patch = ",get_patch
        self.patches[patch.iNetwork] = patch  # deepcopy?
        self.nNetworks += 1

    def add_patch(self, netName, kt, cpNorm, x, y, z):
        #print "self.nNetworks = ",self.nNetworks
        patch = PanairPatch(self.nNetworks, netName, kt, cpNorm, x, y, z, self.log)
        self.msg += patch.process()
        #print "get_patch = ",get_patch
        self.patches[patch.iNetwork] = patch  # deepcopy?
        self.nNetworks += 1

    def find_patch_by_name(self, netName):
        names = []
        for patchID, patch in self.patches.iteritems():
            #self.log.debug("patchID=%s" % (patchID))
            #self.log.debug("*get_patch = %s" %(get_patch))
            #self.log.debug("get_patch.netName=%s" % (get_patch.netName))
            if patch.netName == netName:
                return patch
            names.append(patch.netName)
        msg = 'couldnt findPatchbyName name=|%s| names=%s' % (netName, names)
        raise KeyError(msg)

    def _read_points(self, section):
        """
        @code
        $points - wing-body  with composite panels
        =kn                                               cpnorm
        4.                                                2.
        =kt
        1.
        =nm       nn                                                          netname
        11.       3.                                                          winga
        =x(1,1)   y(1,1)    z(1,1)    x(*,*)    y(*,*)    z(*,*)
           69.4737    9.2105    0.0000   63.7818    9.5807    0.7831
        @endcode
        """
        nNetworks = int(float(section[1][0:10]))
        cpNorm = section[1][50:60].strip()
        #    cpNorm    = float(section[1][50:60])

        kt = int(float(section[2][0:10]))
        #if cpNorm:
            #self.log.debug("nNetworks=%s cpNorm=%s" % (nNetworks, cpNorm))
        #else:
            #self.log.debug("nNetworks=%s" % (nNetworks))

        n = 4
        self.msg += '      kn,kt            %i          %i\n' % (nNetworks, kt)

        for iNetwork in xrange(nNetworks):
            #self.log.debug("lines[* %s] = %s" % (n - 1, section[n - 1]))
            nm = int(float(section[n - 1][0:10]))
            nn = int(float(section[n - 1][10:20]))
            netName = section[n - 1][70:80]
            #self.log.debug("kt=%s nm=%s nn=%s netname=%s" % (
                #kt, nm, nn, netName))

            x = zeros([nm, nn])
            y = zeros([nm, nn])
            z = zeros([nm, nn])
            nFullLines = nm / 2
            nPartialLines = nm % 2
            nLines = nFullLines + nPartialLines
            #print "nFullLines=%s nPartialLines=%s nLines=%s" %(nFullLines,nPartialLines,nLines)
            for j in xrange(nn):
                lines = section[n:n + nLines]
                n += nLines
                #print '\n'.join(lines)
                points = self.split_points(lines, nFullLines, nPartialLines)

                for i, point in enumerate(points):
                    x[i][j] = point[0]
                    y[i][j] = point[1]
                    z[i][j] = point[2]

            #print "--X--"
            #print x
            self.add_patch(netName, kt, cpNorm, x, y, z)
            n += 1
        return True

    def _read_circular_section(self, section):
        """
        @code
        $circular sections - nacelle with composite panels
        =kn
        2.
        =kt
        1.
        =nopt                                                                 netname
        0.                                                                    cowlu
        =nm
        20.
        =xs(1)    ri(1)     xs(2)     ri(2)     xs(*)     ri(*)
            2.0000    2.3000    1.5756    2.3000    1.1486    2.3000
            0.7460    2.3030    0.4069    2.3286    0.1624    2.3790
            0.0214    2.4542   -0.0200    2.5485    0.0388    2.6522
            0.2056    2.7554    0.4869    2.8522    0.8883    2.9413
            1.4250    3.0178    2.1188    3.0656    2.9586    3.0658
            3.8551    3.0175    4.6715    2.9439    5.3492    2.8700
            6.0000    2.7842    6.4687    2.7442
        =nn
        5.
        =th(1)    th(2)     th(3)     th(4)     th(5)
        -90.      -45.      0.        45.       90.
        @endcode
        """
        nNetworks = int(float(section[1][0:10]))
        cpNorm = section[1][50:60].strip()

        kt = int(float(section[2][0:10]))
        #if cpNorm:
            #self.log.debug("nNetworks=%s cpNorm=%s" % (nNetworks, cpNorm))
        #else:
            #self.log.debug("nNetworks=%s" % (nNetworks))

        n = 3
        self.msg += '      kn,kt            %i          %i\n' % (nNetworks, kt)

        for iNetwork in xrange(nNetworks):
            #self.log.debug("lines[%s] = %s" % (n, section[n]))
            # 0-noDisplacement; 1-Specify
            nDisplacement = int(float(section[n][0:10]))
            assert nDisplacement in [0, 1], section[n]
            n += 1

            #print("\nsection[n] = ", section[n].strip())
            nXR = int(float(section[n][0:10]))
            assert nXR > 0, section[n]
            #print("nXR = %s" % nXR)
            netName = section[n][70:80]
            #self.log.debug("kt=%s nXR=%s nDisplacement=%s netname=%s" %
                           #(kt, nXR, nDisplacement, netName))

            nFullLines = nXR / 3
            nPartialLines = int(ceil(nXR % 3 / 3.))
            nXR_Lines = nFullLines + nPartialLines
            #print "X,Rad - nFullLines=%s nPartialLines=%s nLines=%s\n" % (
                #nFullLines, nPartialLines, nXR_Lines)
            n += 1

            Xin = []
            R = []
            lines = section[n:n + nXR_Lines]
            n += nXR_Lines

            for line in lines:
                print(line)
                try:
                    (x1, r1) = float(line[0:10]), float(line[10:20])
                    #print("x1=%s r1=%s" % (x1, r1))
                    Xin.append(x1)
                    R.append(r1)
                except ValueError:
                    pass

                try:
                    (x2, r2) = float(line[20:30]), float(line[30:40])
                    #print("x2=%s r2=%s" % (x2, r2))
                    Xin.append(x2)
                    R.append(r2)
                except ValueError:
                    pass

                try:
                    (x3, r3) = float(line[40:50]), float(line[50:60])
                    #print("x3=%s r3=%s" % (x3, r3))
                    Xin.append(x3)
                    R.append(r3)
                except ValueError:
                    pass
            assert nXR == len(Xin), 'nXR=%s Xin=%s' % (nXR, len(Xin))

            #----------------------------------------------------
            #print("section[n] = ", section[n].strip())
            nTheta = int(float(section[n][0:10]))
            assert nTheta > 0, section[n]
            #print("nTheta = %s" % nTheta)
            nFullLines = nTheta / 6
            nPartialLines = int(ceil(nTheta % 6 / 6.))
            nTheta_Lines = nFullLines + nPartialLines
            #print("Theta - nFullLines=%s nPartialLines=%s nLines=%s" % (
            #    nFullLines, nPartialLines, nTheta_Lines))
            n += 1

            lines = section[n:n + nTheta_Lines]
            theta = []
            for line in lines:
                print(line)
                try:
                    (theta1) = float(line[0:10])
                    theta.append(theta1)
                except ValueError:
                    pass
                try:
                    (theta2) = float(line[10:20])
                    theta.append(theta2)
                except ValueError:
                    pass
                try:
                    (theta3) = float(line[20:30])
                    theta.append(theta3)
                except ValueError:
                    pass
                try:
                    (theta4) = float(line[30:40])
                    theta.append(theta4)
                except ValueError:
                    pass
                try:
                    (theta5) = float(line[40:50])
                    theta.append(theta5)
                except ValueError:
                    pass
                try:
                    (theta6) = float(line[50:60])
                    theta.append(theta6)
                except ValueError:
                    pass
            #print("theta = ", theta)
            n += nTheta_Lines

            assert nTheta == len(theta)
            sinThetaR = [sin(radians(thetai)) for thetai in theta]
            cosThetaR = [cos(radians(thetai)) for thetai in theta]

            zi = 0.
            X = zeros([nXR, nTheta])
            Y = zeros([nXR, nTheta])
            Z = zeros([nXR, nTheta])
            #print("Xin=%s \nR=%s" % (Xin, R))
            #print("X.shape = ", X.shape)
            for (i, x, r) in izip(count(), Xin, R):
                for (j, sinTheta, cosTheta) in izip(count(), sinThetaR, cosThetaR):
                    y = r * sinTheta
                    z = r * cosTheta + zi
                    #print("x=%s y=%s z=%s" % (x, y, z))
                    X[i][j] = x
                    Y[i][j] = y
                    Z[i][j] = z

            #for i,point in enumerate(points):
            #    x[i][j] = point[0]
            #    y[i][j] = point[1]
            #    z[i][j] = point[2]

            #print "--X--"
            #print x
            self.add_patch(netName, kt, cpNorm, X, Y, Z)
        return True

    def _read_grid(self, section):
        self.gridSection = '\n'.join(section) + '\n'
        return True

    def _read_xyz(self, section):
        self.xyzSection = '\n'.join(section) + '\n'
        return True

    def _read_streamlines(self, section):
        self.streamlineSection = '\n'.join(section) + '\n'
        return True

    def _read_printout(self, section):
        """
        @code
        isings  igeomp  isingp  icontp  ibconp  iedgep
        ipraic  nexdgn  ioutpr  ifmcpr  icostp
        @endcode
        """
        #self.printoutSection = '\n'.join(section)+'\n'

        #  ROW 1
        self.isings = float(section[1][0:10])
        #self.isings = 5.
        #2.0 output wake singularity grid-singularity strength and gradients at 9 pts/panel
        #3.0 add table of singularity parameter values
        #4.0 add network array of singularity values
        #5.0 add singularity grid for all network

        self.igeomp = float(section[1][10:20])
        self.igeomp = 1.
        # 1.0 outputs panel data - geometry diagnostic data

        self.isingp = float(section[1][20:30])
        #self.isignp = 1.  # ???
        # 1.0 outputs singularity spline data

        self.icontp = float(section[1][30:40])
        self.icontp = 1.
        # 1.0 outputs control pt location, upper surface unit normals, control point
        #     diagnositc data, control pt maps and singularity parameter maps
        # 2.0 add additional control point map for each network

        self.ibconp = float(section[1][40:50])
        #self.ibconp = 0.0
        # 1.0 outputs boundary condition coeffs, diagnositc data,
        #     boundary condition maps, control point maps,
        #     and singularity paramter maps

        self.iedgep = float(section[1][50:60])
        self.iedgep = .0
        # 1.0 outputs edge-matching diagnositic data
        #print "igeomp=|%s| isingp=|%s| icontp=|%s| ibconp=|%s| iedgep=|%s|" %(igeomp,isingp,icontp,ibconp,iedgep)

        #  ROW 2
        self.ipraic = float(section[2][0:10])
        #self.ipraic = 3.
        # xxx inputs control point sequence number for
        #     which AIC values are to be pritned (rarely used)

        self.nexdgn = float(section[2][10:20])
        self.nexdgn = 0.
        # 1.0 outputs edge and extra control point data

        self.ioutpr = float(section[2][20:30])
        self.ioutpr = -1.0
        # -1.0 omits flow parameter output
        #  0.0 outputs 49 flow parameters
        #  1.0 outputs 13 flow parameters

        self.ifmcpr = section[2][30:40].strip()
        self.ifmcpr = 0.
        # 0.0 omits network force and moment output
        # 1.0 outputs network forces and moments summary per column,
        #     per network, and accumulation over all previous networks

        self.icostp = section[2][40:50].strip()
        self.icostp = 0.
        #self.icostp = 1.0
        #print "ipraic=|%s| nexdgn=|%s| ioutpr=|%s| ifmcpr=|%s| ifmcpr=|%s|" %(ipraic,nexdgn,ioutpr,ifmcpr,iedgep)
        # 1.0 outputs detailed job cost for different segments of the
        #     program (rarely used)

        #$printout options
        #=isings   igeomp    isingp    icontp    ibconp    iedgep
        #4.        0.        0.        1.        1.        0.
        #=ipraic   nexdgn    ioutpr    ifmcpr
        #.0        .0        1.        0.                  3.
        return True

    def _read_trailing_wakes(self, section):
        """
        @code
        $trailing wakes from body
        =kn                                               cpnorm
        2.
        =kt       matcw
        18.       1.
        =inat     insd      xwake     twake                                   netname
        bodyl     3.        100.      .0                                      bodylwk
        bodyu     3.        100.      .0                                      bodyuwk
        @endcode
        """
        #return True  # disable wakes
        nNetworks = int(float(section[1][0:10]))
        cpNorm = section[1][50:60].strip()
        if cpNorm:
            cpNorm = float(cpNorm)
        else:
            cpNorm = 0

        #print "section[2] = ",section[2]
        kt = int(float(section[2][0:10]))

        matchw = section[2][10:20].strip()
        if matchw:
            matchw = float(matchw)
        #if cpNorm:
            #print "nNetworks=%s cpNorm=%s" %(nNetworks,cpNorm)
        #else:
            #print "nNetworks=%s" %(nNetworks)
        #print ""
        #matcw = int(float(section[2][10:20]))

        n = 3  # line number
        self.msg += '      kn,kt            %i          %i\n' % (nNetworks, kt)
        #self.log.debug('kt=%s cpNorm=%s matchw=%s' % (kt, cpNorm, matchw))
        for iNetwork in xrange(nNetworks):
            #self.log.debug("lines[* %s] = %s" % (n, section[n]))

            trailedPanel = section[n][0:10].strip()
            edgeNumber = int(float(section[n][10:20]))
            xWake = float(section[n][20:30])  # x distance
            tWake = float(section[n][30:40])  # 0-wake parallel to x axis;  1-wake in direction of compressibility
            netName = section[n][70:80].strip()
            self.log.debug('trailedPanel=%s edgeNumber=%s xWake=%s tWake=%s netName=%s' % (trailedPanel, edgeNumber, xWake, tWake, netName))
            try:
                patch = self.find_patch_by_name(trailedPanel)
            except KeyError:
                self.log.debug('trailedPanel isnt defined...trailedPanel=|%s|' % (trailedPanel))
                raise

            #xPoints = patch.x
            #yPoints = patch.y
            #print "xPoints = ",xPoints
            #print "yPoints = ",yPoints
            (p1, x1, y1, z1) = patch.get_edge(edgeNumber)

            nPoints = len(x1)
            X = zeros([2, nPoints])
            Y = zeros([2, nPoints])
            Z = zeros([2, nPoints])

            X[0][:] = x1
            Y[0][:] = y1
            Z[0][:] = z1

            if tWake == 0.:
                X[1][:] = ones([nPoints]) * xWake
                Y[1][:] = y1
                Z[1][:] = z1
            else:
                #alphaC, betaC
                raise NotImplementedError('tWake isnt supported')
            self.log.debug("--X---")
            self.log.debug(X)
            self.log.debug("--Y---")
            self.log.debug(Y)
            self.log.debug("--Z---")
            self.log.debug(Z)
            #nm = int(float(section[n-1][0 :10]))
            #nn = int(float(section[n-1][10:20]))
            options = [kt, cpNorm, matchw, trailedPanel, edgeNumber, xWake,
                       tWake]
            patch = self.add_wake_patch(netName, options, X, Y, Z)
            print('----------------------------')
            n += 1
        return True

    def _read_partial_edge_abutments(self, section):
        self.peaSection = '\n'.join(section) + '\n'
        return True

    def write_liberalized_abutments(self):
        msg = '$eat - liberalized abutments\n'
        msg += '%10s' * 6 % (self.epsgeo, self.igeoin, self.igeout,
                             self.nwxref, self.triint, self.iabsum) + '\n'
        return msg

    def _read_liberalized_abutments(self, section):
        #self.log.debug("section[1] = %s" % (section[1]))
        self.epsgeo = section[1][0:10].strip()
        if self.epsgeo:
            self.epsgeo = float(self.epsgeo)
        else:
            self.epsgeo = 0.0
        # absolute value of tolerance used to establith equivalent network
        # edge points; minimum panel diagonal times 0.001 <= espgeo <= minimum
        # panel diagonal times 0.03; to not be restricted use negative numbers

        try:
            self.igeoin = float(section[1][10:20])
        except ValueError:
            self.igeoin = -1.
        # prints network geometry before $EAT
        #  0.0 prints geometry
        # -1.0 omits printout

        try:
            self.igeout = float(section[1][20:30])
        except ValueError:
            # prints network geometry after $EAT
            # 0.0 omits printout
            # 1.0 prints geometry
            self.igeout = 0.

        try:
            self.nwxref = float(section[1][30:40])
        except ValueError:
            #self.nwxref = 1.0
            # prints abutment cross-reference; for each network, all abutments and
            # abutment intersections are described in a format similar to that
            # used in the abutment summary and abutment intersection summary output sections
            # 0.0 omits printout
            # 1.0 prints geometry
            self.nwxref = 0.

        try:
            self.triint = float(section[1][40:50])
        except ValueError:
            # optionally checks for intersection of network edges through other panels
            # for final geometry; used only if needed on complex configurations; can
            # increase the cost of a computer run by 1 or 2 times that of a datachek (2.) run
            # 0.0 omits intersection checking
            # 1.0 checks intersection
            self.triint = 0.0

        try:
            self.iabsum = float(section[1][50:60])
        except ValueError:
            #  0.0 abutment summary
            #  1.0 plus intersection summary
            #  2.0 plus complete list of each network and associated abutments
            #  3.0 plus diagnostic data for program developer
            # -1.0 no abutment printout
            self.iabsum = -0.0
        return True

    def _read_sectional_properties(self, section):
        self.sectionalPropSection = '\n'.join(section) + '\n'
        return True

    def _read_flowfield_properties(self, section):
        self.flowSection = '\n'.join(section) + '\n'
        return True

    #def get(self,section):
    #    pass

    #def get(self,section):
    #    pass

    #def get(self,section):
    #    pass

    def group_sections(self, sections, sectionNames):
        #self.Points = []
        #self.Streamlines = []
        #self.Trailing = []

        #self.Title = None
        #self.dataCheck = None
        #self.symmetry = None
        #self.mach = None
        #self.cases = None
        #self.alphas = None

        #self.pea = None
        #self.eat = None

        #self.printout = None
        #self.end = None

        #self.Grid = []

        sectionMap = {
            'tit': self._read_title,
            'ref': self._read_reference_quantities,
            'dat': self._read_data_check,
            'sym': self._read_symmetry,
            'mac': self._read_mach,
            'cas': self._read_cases,

            'ang': self._read_alphas,
            'yaw': self._read_betas,
            'pri': self._read_printout,
            'poi': self._read_points,
            'cir': self._read_circular_section,
            'tra': self._read_trailing_wakes,
            'pea': self._read_partial_edge_abutments,
            'eat': self._read_liberalized_abutments,
            'sec': self._read_sectional_properties,
            'flo': self._read_flowfield_properties,
            'xyz': self._read_xyz,
            'gri': self._read_grid,
            'str': self._read_streamlines,
            'end': self._read_end,
        }
        #validMaps = ['tit','ref','dat','sym','mac','cas','ang','poi','tra','end']
        #validMaps += ['xyz','str','flo','sec','pri','pea','eat','gri']
        #validMaps = ['yaw','poi','tra']
        validMaps = sectionMap.keys()

        for section, sectionName in izip(sections, sectionNames):  # 1st line
            self.msg += '  $%s\n' % (sectionName)
            #print "section = ",len(section)
            #self.log.debug("sectionName=%s" %(sectionName))
            if sectionName in validMaps:
                #self.log.debug("section[0] = %s" % (section[0]))
                functionMap = sectionMap[sectionName]
                ran = functionMap(section)
                assert ran == True, '%s didnt run' % (sectionName)
                #self.log.debug("")
            if sectionName == 'end':
                break
        return sections

    def split_into_sections(self, lines):
        sections = []
        sectionNames = []
        section = None
        sectionName = None
        for line in lines:
            if '$' in line:
                sections.append(section)
                sectionNames.append(sectionName)
                section = []
                sectionName = line[1:4]
                if sectionName == 'end':
                    break
                #print "new section...name=%s line - %s" %(sectionName,line)
            section.append(line)

        # end section
        section.append(line)
        sections.append(section)
        sectionNames.append(sectionName)

        # get rid of first dummy section
        sections = sections[1:]
        sectionNames = sectionNames[1:]

        assert len(sections) == len(sectionNames), "%s %s" % (
            len(sections), len(sectionNames))

        #for section in sections: # 1st line
            #print "section = ",len(section)
            #print "section[0] = ",section[0]
        return (sections, sectionNames)

    def read_panair(self):
        infile = open(self.infileName, 'r')
        self.lines = infile.readlines()
        infile.close()

        lines = self.remove_comments(self.lines)
        (sections, sectionNames) = self.split_into_sections(lines)
        groups = self.group_sections(sections, sectionNames)
        #self.log.debug("nPatches = %s" % (self.nPatches()))
        # split into headings
        #for panel in panels:
        #    points = readPoints()
        #print "self.msg = ",self.msg

    def getPointsElementsRegions(self):
        points = []
        elements = []
        regions = []
        pointI = 0
        for name, panel in sorted(self.patches.iteritems()):
            #panel = self.patches[2]
            (pointsI, pointi) = panel.get_points()
            (elementsI) = panel.get_elements(pointI)
            #print 'panel.iNetwork=%r' % (panel.iNetwork + 1)
            #regions += panel.iNetwork * ones(len((elementsI)), 'int32')
            regions += [panel.iNetwork + 1] * len(elementsI)

            #print("elementsI = ",elementsI)
            points += pointsI
            elements += elementsI
            pointI += pointi
            #break
            #print "name=%r pointi=%s len(AllElements)=%s len(allPoints)=%s" %(name, pointi, len(elements),len(points))

            if pointI != 0 and False:
                if panel.netName in ['bodylwk', 'awbw']:
                    for point in pointsI:
                        print(point)
                    for element in elementsI:
                        (n1, n2, n3, n4) = element
                        print(element)
                        print(points[n1], points[n2], points[n3], points[n4])
                        print("")
        return points, elements, regions

    def _read_cases(self, section):
        """
        $cases - no. of solutions
        =nacase
        1.
        """
        self.ncases = int(float(section[1][0:10]))
        #self.log.debug("ncases = %s" % (self.ncases))
        return True

    def _read_mach(self, section):
        """
        $mach number
        =amach
        .6
        """
        self.mach = float(section[1][0:10])
        #self.log.debug("mach = %s" % (self.mach))
        return True

    def set_mach(self, mach):
        self.mach = mach

    def write_mach(self):
        out = '$mach number\n'
        out += '%-10s\n' % self.mach
        return out

    def write_cases(self):
        out = ''
        if self.ncases is not None:
            out = '$cases - number of solutions\n'
            out += '%-10s' % sInt(self.ncases) + '\n'
        return out

    def _read_alphas(self, section):
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
        #self.log.debug("alphaC=%s alphas=%s" % (self.alphaC, self.alphas))
        return True

    def set_alphas(self, alphas, alphaC):
        self.alphaC = alphaC
        self.alphas = alphas
        self.ncases = len(alphas)

    def write_title(self):
        out = '$title\n'
        out += '%s\n' % self.title.strip()
        return out

    def write_alphas(self):
        out = '$angles-of-attack\n'
        out += '%-s\n' % self.alphaC
        out += '%-10s' * len(self.alphas) % (tuple(self.alphas)) + '\n'
        return out

    def _read_betas(self, section):
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
        #self.log.debug("betaC=%s betas=%s" % (self.betaC, self.betas))
        return True

    def set_betas(self, betas, betaC):
        self.betaC = betaC
        self.betas = betas
        self.ncases = len(betas)

    def write_betas(self):
        out = '$yaw\n'
        out += '%s\n' % self.betaC
        out += '%-10s' * len(self.betas) % (tuple(self.betas)) + '\n'
        return out

    def _read_reference_quantities(self, section):
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
        #self.log.debug("xref=%s yref=%s zref=%s" % (self.xref,
        #                                            self.yref, self.zref))
        #self.log.debug("sref=%s bref=%s cref=%s dref=%s " % (
        #    self.sref, self.bref, self.cref, self.dref))
        return True

    def write_reference_quantities(self):
        out = '$references for accumulated forces and moments\n'
        out += '=%-10s%-10s%s\n' % ('Xref', 'Yref', 'Zref')
        out += '%-10s%-10s%-s\n' % (self.xref, self.yref, self.zref)
        out += '=%-10s%-10s%-10s%s\n' % ('Sref', 'Bref', 'Cref', 'Dref')
        out += '%-10s%-10s%-10s%s\n' % (self.sref, self.bref, self.cref,
                                        self.dref)
        return out

    def _read_end(self, section):
        self.isEnd = True
        #self.log.debug("end...")
        return True

    def write_end(self):
        if self.isEnd:
            return '$end of panair inputs\n '
        return ''

    def write_data_check(self):
        msg = '$datacheck\n'
        msg += '%s.\n' % (self.dataCheck)
        return msg

    def write_printout(self):
        msg = '$printout options\n'
        msg += "%-10s%-10s%-10s%-10s%-10s%-10s\n" % (self.isings, self.igeomp,
                                                     self.isingp, self.icontp, self.ibconp, self.iedgep)
        msg += "%-10s%-10s%-10s%-10s%-10s\n" % (self.ipraic,
                                                self.nexdgn, self.ioutpr, self.ifmcpr, self.icostp)
        return msg


if __name__ == '__main__':
    import sys
    infileName = sys.argv[1]
    outfileName = infileName[:-4] + '_new.inp'
    print "outfile_name = ", outfileName
    grid = PanairGrid(infileName)
    grid.read_panair()
    grid.write_panair(outfileName)
    (points, elements, regions) = grid.getPointsElementsRegions()
    print("\ninfileName=%s" % infileName)
