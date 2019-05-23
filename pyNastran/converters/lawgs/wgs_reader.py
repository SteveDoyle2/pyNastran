"""
defines:
 - Panel()
 - read_lawgs(wgs_filename, log=None, debug=False)
 - LaWGS(self, log=None, debug=False)

"""
import copy
from math import sin, cos

import numpy as np
from numpy import array, radians, dot, zeros
from cpylog import get_logger2

class Panel:
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
    def __init__(self, key, header, lines, log):
        #print("key=%s \nheader=|%s|" % (key, header))   # ,iSymG
        (ID, nline, npnt, unused_isyml, rx, ry, rz, tx, ty, tz, xscale,
         yscale, zscale, unused_isymg) = header.strip().split()
        self.log = log
        log.debug("ID=%s name=%s imax=%s jmax=%s" % (ID, key, nline, npnt))
        log.debug("Rotate    = <%s,%s,%s>" % (rx, ry, rz))
        log.debug("Translate = <%s,%s,%s>" % (tx, ty, rz))
        log.debug("Scale     = <%s,%s,%s>" % (xscale, yscale, zscale))

        self.name = key
        self.rotate = array([rx, ry, rz], dtype='float64')
        self.translate = array([tx, ty, tz], dtype='float64')
        self.scale = array([xscale, yscale, zscale], dtype='float64')
        npnt = int(npnt)
        nline = int(nline)

        self.ncols = npnt
        self.nrows = nline
        ngroup_lines = npnt // 2  # self.ncols
        is_single_line = False
        if npnt % 2 == 1:  # self.ncols
            is_single_line = True
            #ngroup_lines += 1

        #X = array((nrows, ncols))
        #Y = array((nrows, ncols))
        #Z = array((nrows, ncols))
        i = 0
        iline = 0
        points = np.zeros((self.nrows * self.ncols, 3), dtype='float64')
        for unused_irow in range(self.nrows):
            for unused_row in range(ngroup_lines):
                line = lines[iline]
                (x1, y1, z1, x2, y2, z2) = line.strip().split()
                points[i, :] = [x1, y1, z1]
                points[i + 1, :] = [x2, y2, z2]
                iline += 1
                i += 2

            if is_single_line:
                line = lines[iline]
                (x1, y1, z1) = line.strip().split()
                points[i, :] = [x1, y1, z1]
                iline += 1
                i += 1

        #for i,point in enumerate(points):
            #print("point[%s] = %s" % (i, point))

        # put in plot3d format
        points = points.reshape(self.nrows, self.ncols, 3)
        self.points = points

    def build_rotation_matrix(self, r):
        """
        Form the rotation matrix used for geometrical transformations
        Taken from NASA TM 85767 defining LaWGS.

        """
        # rotation angles, degrees
        #r = radians([self.phi,self.theta,self.psi])
        #print("self.rotate = ",self.rotate)
        r = [radians(ri) for ri in self.rotate]
        cphi = cos(r[0])
        sphi = sin(r[0])
        ctheta = cos(r[1])
        stheta = sin(r[1])
        cpsi = cos(r[2])
        spsi = sin(r[2])

        rot = zeros((3, 3), dtype='float64')
        #print(rot)
        rot[0, 0] = ctheta * cpsi
        rot[1, 0] = ctheta * spsi
        rot[2, 0] = -stheta

        rot[0, 1] = -spsi * cphi + stheta * cpsi * sphi
        rot[1, 1] = cpsi * cphi + stheta * spsi * sphi
        rot[2, 1] = ctheta * sphi

        rot[0, 2] = spsi * sphi + stheta * cpsi * cphi
        rot[1, 2] = -cpsi * sphi + stheta * spsi * cphi
        rot[2, 2] = ctheta * cphi

        return rot

    def update_points(self):
        rot = self.build_rotation_matrix(self.rotate)
        scale = self.scale
        translate = self.translate

        points = self.points
        points2 = copy.deepcopy(points)
        self.log.debug("size(points) = (%s,%s)\n" % (len(points), len(points[0])))
        for i in range(self.nrows):
            for j in range(self.ncols):
                points2[i][j] = scale * (dot(rot, points[i][j]) + translate)
        self.points = points2

    def get_points(self):
        points = self.points.reshape(self.nrows * self.ncols, 3)
        npoints = self.nrows * self.ncols
        return points.tolist(), npoints

    def get_elements(self, pointI):
        n = (self.nrows - 1) * (self.ncols - 1)
        #print("name=%s n=%s nout=%s" %(self.name,n))
        elements = []

        for i in range(self.nrows - 1):
            for j in range(self.ncols - 1):
                elements.append(self.get_element(pointI, i, j))
        assert n == len(elements)
        return elements, n

    def get_element(self, pointI, j, i):
        p1 = (j) * self.ncols + (i) + pointI
        p2 = (j) * self.ncols + (i + 1) + pointI
        p3 = (j + 1) * self.ncols + (i + 1) + pointI
        p4 = (j + 1) * self.ncols + (i) + pointI
        return [p1, p2, p3, p4]

    def write_as_plot3d(self, p3d_file):
        """writes a plot3d section"""
        X = []
        Y = []
        Z = []
        for j in range(self.ncols):
            for i in range(self.nrows):
                (x, y, z) = self.points[i][j]
                X.append(x)
                Y.append(y)
                Z.append(z)

        msg = ''
        for x in X:
            msg += '%s ' % (x)
        p3d_file.write(msg + '\n')

        msg = ''
        for y in Y:
            msg += '%s ' % (y)
        p3d_file.write(msg + '\n')

        msg = ''
        for z in Z:
            msg += '%s ' % (z)
        p3d_file.write(msg + '\n')

def read_lawgs(wgs_filename, log=None, debug=False):
    """reads an lawgs file"""
    model = LaWGS(log=log, debug=debug)
    model.read_lawgs(wgs_filename)
    return model

class LaWGS:
    """defines a reader for the LaWGS legacy file format"""
    model_type = 'LaWGS'

    def __init__(self, log=None, debug=False):
        """
        Initializes the LaWGS object

        Parameters
        ----------
        debug : bool/None; default=True
            used to set the logger if no logger is passed in
                True:  logs debug/info/error messages
                False: logs info/error messages
                None:  logs error messages
        log : logging module object / None
            if log is set, debug is ignored and uses the
            settings the logging object has

        """
        self.log = get_logger2(log=log, debug=debug)
        self.panels = {}

    def read_lawgs(self, wgs_filename):
        """reads an lawgs file"""
        with open(wgs_filename, 'r') as lawgs_file:
            lines = lawgs_file.readlines()

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
        for key, header_group in sorted(groups.items()):
            header, group = header_group
            #if key=='BODY':
            #if 1:
            panel = Panel(key, header, group, log=self.log)
            panel.update_points()
            self.panels[key] = panel

    def get_points_elements(self):
        points, elements, unused_regions = self.get_points_elements_regions()
        return points, elements

    def get_points_elements_regions(self):
        points = []
        elements = []
        regions = []
        pointI = 0
        iregion = 0
        for (unused_name, panel) in sorted(self.panels.items()):
            (pointsI, pointi) = panel.get_points()
            (elementsI, n) = panel.get_elements(pointI)
            points += pointsI
            elements += elementsI
            pointI += pointi
            regions += [iregion] * n
            iregion += 1
            #print("name=%s len(All_elements)=%s len(all_points)=%s" % (
                #name, len(elements), len(points)))

        #for point in points:
            #print(point)
        return points, elements, regions

    def write_as_plot3d(self, p3dname):
        """writes a plot3d file"""
        with open(p3dname, 'w') as p3d_file:
            p3d_file.write('%s\n' % (len(self.panels)))
            for (unused_name, panel) in sorted(self.panels.items()):
                p3d_file.write('%s %s 1\n' % (panel.nrows, panel.ncols))

            for (unused_name, panel) in sorted(self.panels.items()):
                panel.write_as_plot3d(p3d_file)
