import numpy as np
from pyNastran.converters.panair.assign_type import integer


def print_float(value):  # sInt #string_float_value
    """int represented as a short float"""
    value = "%f" % value
    return value.rstrip('0')


class PanairPatch:
    """Used for physical surfaces"""
    def __init__(self, inetwork, network_name, kt, cp_norm, xyz, log):
        self.log = log

        self.inetwork = inetwork  # lets it print out in order, does it need a deepcopy???
        self.network_name = network_name.strip()
        #self.log.debug('network_name=%r' % (network_name))
        #self.log.debug("****patch.network_name=%s" % (self.network_name))
        self.kt = kt
        if cp_norm == '':
            cp_norm = 0
        self.cp_norm = integer(cp_norm, 'cp_norm')
        self.matchw = 0
        self.xyz = xyz
        self.nrows = xyz.shape[0]
        self.ncols = xyz.shape[1]

    def is_wake(self):
        #self.log.debug('is_wake %s %s' % (self.network_name, self.kt))
        if self.kt in [18, 20]:
            return True
        return False

    def write_plot3d(self, f, dim):
        """
        ..todo: is the normal defined correctly?
        ..todo: will this load into tecplot
        """
        self.log.debug("xyz.shape=%s" % str(self.xyz.shape))
        try:
            data = self.xyz[:, :, dim - 1]
        except IndexError:
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
            source = 1
            doublet = 12
            nlopt1 = 5
            nropt1 = 3
            nlopt2 = 7
            nropt2 = -2
            ipot = 2
        elif self.kt == 5:
            source = 1
            doublet = 12
            nlopt1 = 6
            nropt1 = 9
            nlopt2 = 7
            nropt2 = -2
            ipot = 2
        elif self.kt == 18:
            source = 0
            doublet = 18
            nlopt1 = 0
            nropt1 = 9
            nlopt2 = 15
            nropt2 = 2
            ipot = 2
            if self.matchw == 1.:
                nlopt2 = 6
            #self.log.debug("18...matchw = %s" % (self.matchw))
        elif self.kt == 20:
            source = 0
            doublet = 20
            nlopt1 = 0
            nropt1 = 9
            nlopt2 = 6
            nropt2 = 2
            ipot = 2
        else:
            raise NotImplementedError('new kt...kt=%s' % (self.kt))
        pts = self.npoints
        pans = self.npanels
        #cumPts=33
        #cumPn =50
        msg += ' %-10s %4s %7s %7s %3s ' % (
            self.network_name, self.inetwork + 1, self.nrows, self.ncols, self.kt)
        msg += '%4s %5s %7s %7s %7s %7s ' % (
            source, doublet, nlopt1, nropt1, nlopt2, nropt2)
        msg += '%7s %7s %7s %7s ' % (ipot, pts, pans, self.cp_norm)
        msg += '%7s %7s\n' % (cum_pts, cum_pn)
        return msg

    @property
    def npanels(self):
        return (self.nrows - 1) * (self.ncols - 1)

    @property
    def npoints(self):
        return self.nrows * self.ncols

    def get_point(self, row, col):
        return self.xyz[row, col]

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
        self.log.info("edge_number=%s" % edge_number)
        self.log.debug("xyz.shape = %s" % (str(self.xyz.shape)))

        if edge_number == 1:
            xyz = self.xyz[0, :, :] # pretty sure edge 1 is the 0th row
            p = [icol for icol in range(self.ncols)]  # good
        elif edge_number == 2:
            xyz = self.xyz[:, self.ncols - 1, :]  # pretty sure edge 2 is the 0th row
            p = [icol * (self.nrows) + (self.nrows - 1) for icol in range(self.ncols)]
            #p = [iRow*(self.ncols)+(self.ncols-1) for iRow in range(self.nrows)]  #
        elif edge_number == 3:
            xyz = self.xyz[self.nrows - 1, :, :]  # pretty sure edge3 is the last row
            p = [icol + self.nrows * icol for icol in range(self.ncols)]  # good
            #p = [(self.ncols-1)*(self.nrows)+iRow for iRow in range(self.nrows)]
        elif edge_number == 4:
            xyz = self.xyz[:, 0, :]  # pretty sure edge 2 is the 0th row
            p = [self.nrows * icol for icol in range(self.ncols)]  # good
        else:
            raise ValueError('invalid edge; edge_number=%s' % edge_number)


        #self.log.debug("nrows=%s ncols=%s edge_number=%s" % (
            #self.nrows, self.ncols, edge_number))
        p = np.arange(self.npoints)
        for point_id in p:
            #point_id = 2
            unused_p2 = self.get_ipoint(point_id)
            #print("point[%s]=%s" % (point_id, unused_p2))

        return (p, xyz)

    def get_edges(self):
        self.log.debug('get_edges')
        nx = 2 * (self.nrows + self.ncols) - 2
        p = np.zeros(nx)
        xyz = np.zeros((nx, 3), dtype='float32')

        i = 0
        for edge_id in range(1, 4 + 1):
            (p1, xyz1) = self.get_edge(edge_id)
            nx1 = xyz1.shape[0]
            p[i:i + nx1] = p1[:nx1]
            xyz[i:i + nx1] = xyz1[:nx1, :]
            #self.log.debug("-----")
        return (p, xyz)

    def get_elements(self, unused_ipoint):
        #print('nrows=%s ncols=%s' % (self.nrows, self.ncols))
        panels = elements_from_quad(self.ncols, self.nrows)
        return panels

    def get_points(self):
        points = []
        #self.log.debug("size(xyz) = %s" %( str( self.xyz.shape ) ))
        self.log.debug('self.inetwork=%s self.network_name=%r' % (self.inetwork, self.network_name))
        for j in range(self.ncols):
            for i in range(self.nrows):
                point = self.xyz[i, j, :]
                points.append(point)
        return points, len(points)

    def write_as_plot3d(self):
        out = ''
        x = self.xyz[:, :, 0].ravel()  # unravel
        y = self.xyz[:, :, 1].ravel()  # unravel
        z = self.xyz[:, :, 2].ravel()  # unravel

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
        #    npoints_left = nfull_lines*2 + npartial_lines
        #    for r in range(0, self.nrows, 2)
        return out

    #def rotate(self):
        #"""
        #not complete...
        #"""
        #self.x = transpose(self.x)
        #self.y = transpose(self.y)
        #self.z = transpose(self.z)
        #self.x[0:n][:] = self.x[-n:-1][:] # something like this...

    def get_header(self):
        header = '$points - surface panels\n'
        header += '%-10s%-10s\n' % ('1.', self.cp_norm)  # nNetworks is 1
        header += '%-10s\n' % print_float(self.kt)
        header += '%-10s%-10s%50s%-10s\n' % (
            print_float(self.nrows), print_float(self.ncols), '', self.network_name)
        return header

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
        points = ''
        header = self.get_header()

        #nfull_lines = nm // 2
        #npartial_lines = nm % 2
        #nlines = nfull_lines + npartial_lines

        nfull_lines = self.nrows // 2
        npartial_lines = self.nrows % 2
        #nlines = nfull_lines + npartial_lines

        for c in range(self.ncols):
            npoints_left = nfull_lines * 2 + npartial_lines
            for r in range(0, self.nrows, 2):
                if npoints_left > 1:
                    x1, y1, z1 = self.xyz[r, c, :]
                    x2, y2, z2 = self.xyz[r + 1, c, :]
                    points += self.write_points([x1, y1, z1], [x2, y2, z2])
                else:
                    x1, y1, z1 = self.xyz[r, c, :]
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
    """Used for explicit wakes"""
    def __init__(self, inetwork, network_name, options, xyz, log):
        (kt, cp_norm, matchw, trailed_panel, edge_number, xwake, twake) = options
        PanairPatch.__init__(self, inetwork, network_name, kt, cp_norm, xyz, log)

        self.log = log
        self.matchw = matchw
        self.trailed_panel = trailed_panel
        self.edge_number = edge_number
        self.xwake = xwake
        self.twake = twake
        #self.log.debug("matchw = %s" % (self.matchw))
        #self.log.debug("wake patch")

    def is_wake(self):
        return True

    def __repr__(self):
        header = '$trailing wakes\n'
        #points = ''

        header += '%-10s%-10s\n' % (1., self.cp_norm)  # nNetworks is 1
        header += '%-10s%-10s\n' % (print_float(self.kt), self.matchw)
        #header += '%-10s%-10s%40s%-10s\n' %(self.nrows,self.nCols,'',self.netName)
        header += '%-10s%-10s%-10s%-10s%-30s%-10s\n' % (self.trailed_panel,
                                                        print_float(self.edge_number),
                                                        self.xwake,
                                                        print_float(self.twake),
                                                        ' ',
                                                        self.network_name)
        return header

def elements_from_quad(nx, ny):
    """creates quad elements for the gui"""
    assert nx > 1
    assert ny > 1

    nelements = (nx - 1) * (ny - 1)
    npoints = nx * ny

    # create a matrix with the point counter
    ipoints = np.arange(npoints, dtype='int32').reshape((nx, ny))

    # move around the CAERO quad and apply ipoints
    elements = np.zeros((nelements, 4), dtype='int32')
    elements[:, 0] = ipoints[:-1, :-1].ravel()  # (i,  j  )
    elements[:, 1] = ipoints[1:, :-1].ravel()   # (i+1,j  )
    elements[:, 2] = ipoints[1:, 1:].ravel()    # (i+1,j+1)
    elements[:, 3] = ipoints[:-1, 1:].ravel()   # (i,j+1  )
    return elements
