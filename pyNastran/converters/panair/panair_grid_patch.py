from __future__ import print_function
from six.moves import range
from numpy import array, zeros, cross, transpose
from numpy.linalg import norm


def sInt(value):  # sInt #string_float_value
    """
    int represented as a short float
    """
    value = "%f" % value
    return value.rstrip('0')


class PanairPatch(object):
    def __init__(self, inetwork, network_name, kt, cp_norm, x, y, z, log):
        self.log = log

        self.inetwork = inetwork  # lets it print out in order, does it need a deepcopy???
        self.network_name = network_name.strip()
        #self.log.debug('network_name=|%s|' % (network_name))
        #self.log.debug("****patch.network_name=%s" % (self.network_name))
        self.kt = kt
        self.cp_norm = cp_norm
        self.matchw = 0
        self.x = x
        self.y = y
        self.z = z
        shape = x.shape
        self.nrows = shape[0]
        self.ncols = shape[1]

        #self.log.debug("shape = %s" % (str(shape)))

    def is_wake(self):
        print('is_wake %s %s' % (self.network_name, self.kt))
        if self.kt in [18]:
            return True
        return False

    def write_plot3d(self, f, dim):
        """
        ..todo: is the normal defined correctly?
        ..todo: will this load into tecplot
        """
        print("x.shape=%s" % str(self.x.shape))
        if dim == 1:
            data = self.x
        elif dim == 2:
            data = self.y
        elif dim == 3:
            data = self.z
        else:
            raise RuntimeError('dim=1 -> x; dim=2 -> y; dim=3 -> z')
        ni, nj = data.shape
        for j in range(nj):
            msg = ''
            for i in range(ni):
                msg += '%s ' % data[i, j]
            f.write(msg)
        msg += '\n'
        f.write(msg)

    def process(self):
        msg = '     network # being processed %3i\n\n' % (self.inetwork + 1)
        return msg

    def quick_summary(self, cum_pts, cum_pn):
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
            #self.log.debug("18...matchw = %s" % (self.matchw))
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
        msg += ' %-10s %4s %7s %7s %3s ' % (self.network_name,
                                            self.inetwork + 1, self.nrows, self.ncols, self.kt,)
        msg += '%4s %5s %7s %7s %7s %7s ' % (
            src, dblt, nlopt1, nropt1, nlopt2, nropt2)
        msg += '%7s %7s %7s %7s ' % (ipot, pts, pans, self.cp_norm)
        msg += '%7s %7s\n' % (cum_pts, cum_pn)
        return msg

    def nPanels(self):
        return (self.nrows - 1) * (self.ncols - 1)

    def nPoints(self):
        return self.nrows * self.ncols

    def get_panel_points(self, ipanel):
        r = ipanel % (self.nrows - 1)
        c = ipanel // (self.nrows - 1)

        #print "r=%s c=%s" %(r,c)
        p1 = self.get_point(r, c)
        p2 = self.get_point(r, c + 1)
        p3 = self.get_point(r + 1, c + 1)
        p4 = self.get_point(r + 1, c)
        return (p1, p2, p3, p4)

    def get_panel_point_IDs(self, ipanel):
        r = ipanel % (self.nrows - 1)
        c = ipanel // (self.nrows - 1)

        #print "r=%s c=%s" %(r,c)
        p1 = self.get_point_ID(r, c)
        p2 = self.get_point_ID(r, c + 1)
        p3 = self.get_point_ID(r + 1, c + 1)
        p4 = self.get_point_ID(r + 1, c)
        return (p1, p2, p3, p4)

    def get_subpanel_properties(self, ipanel):
        (p1, p2, p3, p4) = self.get_panel_points(ipanel)
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

    def get_panel_properties(self, ipanel):
        (p1, p2, p3, p4) = self.get_panel_points(ipanel)
        a = p1 - p3
        b = p2 - p4
        centroid = 0.25 * (p1 + p2 + p3 + p4)

        N = cross(a, b)
        normN = norm(N)
        n = N / normN  # normal vector

        area = 0.5 * normN  # area

        u = (p1 + p2 - p3 - p4) / 2.  # longitudinal
        p = (-p1 + p2 + p3 - p4) / 2.  # transverse

        u = 0.5 * (a + b)  # longitudinal vector in local coordinates
        p = 0.5 * (-a + b)  # transverse vector in local coordinates
        o = cross(n, u)  # normal to both vectors in local coordinates

        diameter = norm(a - b)

        return (area, n, centroid, diameter, u, p, o)

    def get_panel_area_normal(self, ipanel):
        (p1, p2, p3, p4) = self.get_panel_points(ipanel)
        normal = cross(p1 - p3, p2 - p4)
        normi = norm(normal)
        normal /= normi  # normal vector

        area = 0.5 * normi
        return (area, normal)

    def get_panel_area(self, ipanel):
        # ipanel=200
        (p1, p2, p3, p4) = self.get_panel_points(ipanel)
        area = 0.5 * norm(cross(p1 - p3, p2 - p4))
        return area

    def get_point(self, row, col):
        return array([self.x[row][col],
                      self.y[row][col],
                      self.z[row][col]])

    def get_point_ID(self, row, col):
        return col * self.nrows + row

    def get_ipoint(self, ipoint):
        irow = ipoint // self.ncols
        icol = ipoint % self.ncols
        #self.log.debug("ipoint=%s irow=%s icol=%s" % (ipoint, irow, icol))
        return self.get_point(irow, icol)

    def get_edge(self, edge_number):
        r"""
        gets all the points associated with a given edge
        ::
                 edge1
               0  1  2   -> i (row)
         edge4 3  4  5
               6  7  8  edge2
               9  10 11
             |   edge3
             j
        """
        self.log.debug('---get_edge---')
        #edge_number = 2
        #edge_number = 4
        print("edge_number=%s" % edge_number)
        self.log.debug("x.shape = %s" % (str(self.x.shape)))

        if edge_number == 1:
            pass
            #self.log.debug("self.x[:]\n%54s%s" % ('', self.x[:]))
            #self.log.debug("self.y[:]\n%54s%s" % ('', self.y[:]))
            #self.log.debug("self.z[:]\n%54s%s" % ('', self.z[:]))
            #edge_number = 2

        if edge_number == 1:
            x = self.x[0][:]
            y = self.y[0][:]
            z = self.z[0][:]  # pretty sure edge 1 is the 0th row
            p = [iCol for iCol in range(self.ncols)]  # good
        elif edge_number == 2:
            x = self.x[:][self.ncols - 1]
            y = self.y[:][self.ncols - 1]
            z = self.z[:][self.ncols - 1]  # pretty sure edge 2 is the 0th row
            p = [iCol * (self.nrows) + (self.nrows - 1) for icol in range(self.ncols)]
            #p = [iRow*(self.ncols)+(self.ncols-1) for iRow in range(self.nrows)]  #
        elif edge_number == 3:
            x = self.x[self.nrows - 1][:]
            y = self.y[self.nrows - 1][:]
            z = self.z[self.nrows - 1][:]  # pretty sure edge3 is the last row
            p = [iCol + self.nrows * iCol for iCol in range(self.ncols)]  # good
            #p = [(self.ncols-1)*(self.nrows)+iRow for iRow in range(self.nrows)]
        elif edge_number == 4:
            x = self.x[:][0]
            y = self.y[:][0]
            z = self.z[:][0]  # pretty sure edge 2 is the 0th row
            p = [self.nrows * iCol for iCol in range(self.ncols)]  # good
        else:
            raise ValueError('invalid edge; edge_number=%s' % edge_number)

        #print "x", x

        #self.log.debug("nRows=%s nCols=%s edge_number=%s" % (
            #self.nRows, self.nCols, edge_number))
        #print "nx = ",len(x)
        #print "p = ",p
        #print "x = ",x
        #print "y = ",y
        #print "z = ",z
        p = [ipoint for ipoint in range(self.nPoints())]
        for point_id in p:
            #point_id = 2
            p2 = self.get_ipoint(point_id)
            #print("point[%s]=%s" % (point_id, p2))

        return (p, x, y, z)

    def get_edges(self):
        self.log.debug('get_edges')
        nx = 2 * (self.nrows + self.ncols) - 2
        p = zeros(nx)
        x = zeros(nx)
        y = zeros(nx)
        z = zeros(nx)

        i = 0
        for edge_id in range(1, 4 + 1):
            (p1, x1, y1, z1) = self.get_edge(edge_id)
            nx1 = len(x1)
            p[i:i + nx1] = p1[0:nx1]
            x[i:i + nx1] = x1[0:nx1]
            y[i:i + nx1] = y1[0:nx1]
            z[i:i + nx1] = z1[0:nx1]
            #self.log.debug("-----")
        return (p, x, y, z)

    def get_elements(self, ipoint):
        panels = []
        #print("ipoint=%s" % (ipoint))
        for ipanel in range(self.nPanels()):
            panel = self.get_panel_point_IDs(ipanel)
            panel2 = []

            for p in panel:
                panel2.append(p + ipoint)
            #print("panel=%s panel2=%s" % (str(panel), str(panel2)))
            panels.append(panel2)
        return panels

    def get_points(self):
        points = []
        #self.log.debug("size(X) = %s" %( str( self.x.shape ) ))
        #print "size(X) = %s" %( str(X.size())
        self.log.debug('self.inetwork=%s self.network_name=%r' % (self.inetwork, self.network_name))
        for j in range(self.ncols):
            for i in range(self.nrows):
                point = [self.x[i][j], self.y[i][j], self.z[i][j]]
                points.append(point)

        return points, len(points)

    def write_as_plot3d(self):
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
        #self.log.debug(out)
        #print x
        #for c in range(self.nCols):
        #    nPointsLeft = nFullLines*2+nPartialLines
        #    for r in range(0,self.nrows, 2)
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
        #x = self.write_as_plot3d()

        #self.log.debug("*******")
        header = '$points - surface panels\n'
        points = ''

        header += '%-10s%-10s\n' % ('1.', self.cp_norm)  # nNetworks is 1
        header += '%-10s\n' % sInt(self.kt)
        header += '%-10s%-10s%50s%-10s\n' % (
            sInt(self.nrows), sInt(self.ncols), '', self.network_name)

        #nfull_lines = nm/2
        #npartial_lines = nm%2
        #nlines = nfull_lines + npartial_lines

        nfull_lines = self.nrows / 2
        npartial_lines = self.nrows % 2
        #nlines = nfull_lines + npartial_lines

        for c in range(self.ncols):
            npoints_left = nfull_lines * 2 + npartial_lines
            for r in range(0, self.nrows, 2):
                if npoints_left > 1:
                    x1 = self.x[r][c]
                    y1 = self.y[r][c]
                    z1 = self.z[r][c]

                    x2 = self.x[r + 1][c]
                    y2 = self.y[r + 1][c]
                    z2 = self.z[r + 1][c]
                    points += self.write_points([x1, y1, z1], [x2, y2, z2])
                else:
                    x1 = self.x[r][c]
                    y1 = self.y[r][c]
                    z1 = self.z[r][c]
                    points += self.write_point([x1, y1, z1])
                npoints_left -= 2
        return header + points

    def write_points(self, point1, point2):
        point1 = self.fix_point(point1)
        point2 = self.fix_point(point2)

        out = "%-10s" * 6 % (point1[0], point1[1], point1[2],
                             point2[0], point2[1], point2[2])
        return out + '\n'

    def write_point(self, point1):
        point1 = self.fix_point(point1)
        out = "%-10s" * 3 % (point1[0], point1[1], point1[2])
        return out + '\n'

    def fix_point(self, point_in):
        point_out = []
        for value in point_in:
            str_value = '%s' % (value)
            if len(str_value) > 10:
                str_value = str_value[0:9]
            point_out.append(str_value.rstrip('0'))
            #print("str_value=%s len=%s" % (str_value, len(str_value)))
        #print("point_out = ", point_out)
        return point_out


class PanairWakePatch(PanairPatch):
    def __init__(self, inetwork, network_name, options, x, y, z, log):
        (kt, cp_norm, matchw, trailed_panel, edge_number, xWake, tWake) = options
        PanairPatch.__init__(self, inetwork, network_name, kt, cp_norm, x, y, z, log)

        self.log = log
        self.matchw = matchw
        self.trailed_panel = trailed_panel
        self.edge_number = edge_number
        self.xWake = xWake
        self.tWake = tWake
        #self.log.debug("matchw = %s" % (self.matchw))
        #self.log.debug("wake patch")

    def is_wake(self):
        return True

    def __repr__(self):
        header = '$trailing wakes\n'
        #points = ''

        header += '%-10s%-10s\n' % (1., self.cp_norm)  # nNetworks is 1
        header += '%-10s%-10s\n' % (sInt(self.kt), self.matchw)
        #header += '%-10s%-10s%40s%-10s\n' %(self.nrows,self.nCols,'',self.netName)
        header += '%-10s%-10s%-10s%-10s%-30s%-10s\n' % (self.trailed_panel,
                                                        sInt(self.edge_number),
                                                        self.xWake,
                                                        sInt(self.tWake),
                                                        ' ',
                                                        self.network_name)
        return header
