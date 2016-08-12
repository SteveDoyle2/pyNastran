from __future__ import print_function
from six import iteritems
from six.moves import range
import copy
from math import sin, cos
from numpy import array, radians, dot, zeros


class LaWGS_Panel(object):
    """
    Parameters
    ----------
    rotate : float
        rotates the patch
    translate : float
        translates the patch
    scale : float
        scales the patch
    """
    def __init__(self, key, header, group):
        #print("key=%s \nheader=|%s|" % (key, header))   # ,iSymG
        (ID, nline, npnt, isyml, rx, ry, rz, tx, ty, tz, xScale,
         yScale, zScale, iSymG) = header.strip().split()
        print("ID=%s name=%s imax=%s jmax=%s" % (ID, key, nline, npnt))
        print("Rotate    = <%s,%s,%s>" % (rx, ry, rz))
        print("Translate = <%s,%s,%s>" % (tx, ty, rz))
        print("Scale     = <%s,%s,%s>" % (xScale, yScale, zScale))

        self.name = key
        self.rotate = array([rx, ry, rz], 'd')
        self.translate = array([tx, ty, tz], 'd')
        self.scale = array([xScale, yScale, zScale], 'd')
        npnt = int(npnt)
        nline = int(nline)

        self.nCols = npnt
        self.nRows = nline
        nGroupLines = npnt // 2  # self.nCols
        isSingleLine = False
        if npnt % 2 == 1:  # self.nCols
            isSingleLine = True
            #nGroupLines+=1

        #X = array((nRows,nCols))
        #Y = array((nRows,nCols))
        #Z = array((nRows,nCols))
        points = []

        irow = 0
        #print("*****",group[-1])
        for iline in range(self.nRows):
            for row in range(nGroupLines):
                line = group[irow]
                #print(" line = ",line)
                (x1, y1, z1, x2, y2, z2) = line.strip().split()
                points.append(array([x1, y1, z1], 'd'))
                points.append(array([x2, y2, z2], 'd'))
                #print("points[%s]=%s" %(i,points[i]))
                #print("points[%s]=%s" %(i+1,points[i+1]))
                #print("---")
                irow += 1

            if isSingleLine:
                line = group[irow]
                #print("*line = ",line)
                (x1, y1, z1) = line.strip().split()
                points.append(array([x1, y1, z1], 'd'))
                irow += 1

        #for i,point in enumerate(points):
            #print("point[%s] = %s" %(i,point))

        n = 0
        #for n,point in enumerate(points):
        Points = []
        for i in range(self.nRows):
            points2 = []
            for j in range(self.nCols):
                points2.append(points[n])
                #jj = n%self.nCols
                #ii = n/self.nCols
                #print("n=%-2s i=%-2s j=%-2s ii=%-2s jj=%-2s" %(n,i,j,ii,jj))
                n += 1

            Points.append(points2)
            #print("len(points[%s]) = %s" %(j,len(points2)))

        self.points = Points
        #print("len(self.points) = %s" %(len(self.points)))

    def buildRotationMatrix(self, r):
        """
        Form the rotation matrix used for geometrical transformations
        Taken from NASA TM 85767 defining LaWGS.
        """
        # rotation angles, degrees
        #r = radians([self.phi,self.theta,self.psi])
        #print("self.rotate = ",self.rotate)
        r = [radians(ri) for ri in self.rotate]
        cPhi = cos(r[0])
        sPhi = sin(r[0])
        cTheta = cos(r[1])
        sTheta = sin(r[1])
        cPsi = cos(r[2])
        sPsi = sin(r[2])

        rot = zeros((3, 3), 'd')
        #print(rot)
        rot[0, 0] = cTheta * cPsi
        rot[1, 0] = cTheta * sPsi
        rot[2, 0] = -sTheta

        rot[0, 1] = -sPsi * cPhi + sTheta * cPsi * sPhi
        rot[1, 1] = cPsi * cPhi + sTheta * sPsi * sPhi
        rot[2, 1] = cTheta * sPhi

        rot[0, 2] = sPsi * sPhi + sTheta * cPsi * cPhi
        rot[1, 2] = -cPsi * sPhi + sTheta * sPsi * cPhi
        rot[2, 2] = cTheta * cPhi

        return rot

    def updatePoints(self):
        rot = self.buildRotationMatrix(self.rotate)
        scale = self.scale
        translate = self.translate

        points = self.points
        Points2 = copy.deepcopy(points)
        print("size(points) = (%s,%s)\n" % (len(points), len(points[0])))
        for i in range(self.nRows):
            #points2 = []
            for j in range(self.nCols):
                Points2[i][j] = scale * (dot(rot, points[i][j]) + translate)
        self.Points = Points2

    def get_points(self):
        Points = []
        for i in range(self.nRows):
            #points2 = []
            for j in range(self.nCols):
                Points.append(self.Points[i][j])
        return Points, len(Points)

    def get_elements(self, pointI):
        n = (self.nRows - 1) * (self.nCols - 1)
        #print("name=%s n=%s nout=%s" %(self.name,n))
        elements = []

        for i in range(self.nRows - 1):
            for j in range(self.nCols - 1):
                elements.append(self.getElement(pointI, i, j))
        assert n == len(elements)
        return elements, n

    def getElement(self, pointI, j, i):
        p1 = (j) * self.nCols + (i) + pointI
        p2 = (j) * self.nCols + (i + 1) + pointI
        p3 = (j + 1) * self.nCols + (i + 1) + pointI
        p4 = (j + 1) * self.nCols + (i) + pointI
        return [p1, p2, p3, p4]

    def write_as_plot3d(self, f):
        X = []
        Y = []
        Z = []
        for i in range(self.nRows):
            #points2 = []
            for j in range(self.nCols):
                (x, y, z) = self.Points[i][j]
                X.append(x)
                Y.append(y)
                Z.append(z)

        msg = ''
        for x in X:
            msg += '%s ' % (x)
        f.write(msg + '\n')

        msg = ''
        for y in Y:
            msg += '%s ' % (y)
        f.write(msg + '\n')

        msg = ''
        for z in Z:
            msg += '%s ' % (z)
        f.write(msg + '\n')


class LaWGS(object):
    model_type = 'LaWGS'

    def __init__(self, filename='tmx1242.wgs'):
        self.filename = filename

    def readLaWGS(self):
        lawgs = open(self.filename, 'r')
        lines = lawgs.readlines()

        nlines = len(lines)
        groups = {}
        name = ''
        group = []
        header = ''

        i = 1
        while i < nlines:
            line = lines[i].rstrip()
            if line.strip() == '':
                #print("found blank line, breaking...")
                break
            #print(line)
            if "'" in line:
                #print("if")
                groups[name] = [header, group]
                name = line.strip("' \t\r\n")
                header = lines[i + 1].rstrip()
                group = []
                i += 1
            else:
                group.append(line)
            i += 1

        groups[name] = [header, group]

        del groups['']
        self.panels = {}
        for key, headerGroup in sorted(iteritems(groups)):
            header, group = headerGroup
            #if key=='BODY':
            if 1:
                panel = LaWGS_Panel(key, header, group)
                panel.updatePoints()
                self.panels[key] = panel

    def getPointsElements(self):
        points, elements, regions = self.get_points_elements_regions()
        return points, elements

    def get_points_elements_regions(self):
        points = []
        elements = []
        regions = []
        pointI = 0
        iregion = 0
        for (name, panel) in sorted(iteritems(self.panels)):
            (pointsI, pointi) = panel.get_points()
            (elementsI, n) = panel.get_elements(pointI)
            points += pointsI
            elements += elementsI
            pointI += pointi
            regions += [iregion] * n
            iregion += 1
            #print("name=%s len(AllElements)=%s len(allPoints)=%s" %(name,len(elements),len(points)))

        #for point in points:
            #print(point)
        return points, elements, regions

    def write_as_plot3d(self, p3dname):
        f = open(p3dname, 'wb')

        f.write('%s\n' % (len(self.panels)))
        for (name, panel) in sorted(iteritems(self.panels)):
            f.write('%s %s 1\n' % (panel.nRows, panel.nCols))

        for (name, panel) in sorted(iteritems(self.panels)):
            panel.write_as_plot3d(f)

if __name__ == '__main__':  # pragma: no cover
    lawgs = LaWGS('tmx1242.wgs')
    lawgs.readLaWGS()
