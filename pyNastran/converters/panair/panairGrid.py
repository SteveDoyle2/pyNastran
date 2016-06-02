from __future__ import print_function
from six import iteritems
import os
#import sys
#import copy
from six import PY2
from six.moves import zip, range
from itertools import count

from math import ceil, sin, cos, radians
from numpy import array, zeros, ones


from pyNastran.converters.panair.panairGridPatch import PanairPatch, PanairWakePatch, sInt
#from pyNastran.converters.panair.panairWrite import PanairWrite
from pyNastran.utils.log import get_logger
from pyNastran.utils import print_bad_path

#from pyNastran.utils import list_print

#CL = -Fx*sin(alpha)*cos(beta) + Fy*sin(alpha)*sin(beta) +Fz*cos(alpha)
#CD =  Fx*cos(alpha)*cos(beta) - Fy*cos(alpha)*sin(beta) +Fz*sin(alpha)
#CY =  Fx*sin(beta) +Fy*cos(beta)


def fortran_value(value):
    return "%8.4E" % value


class PanairGrid(object):
    model_type = 'panair'

    def __init__(self, log=None, debug=True):
        self.nnetworks = 0
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

        self.epsgeo = 0.0
        self.igeoin = -1.0
        self.igeout = 0.
        self.nwxref = 0.
        self.triint = 0.0
        self.iabsum = -0.0

        self.isEnd = None

        self.pea_section = ''


        self.mach = 0.0
        self.data_check = 2
        self.title_section = ''
        self.xyz_section = ''
        self.streamline_section = ''
        self.flow_section = ''
        self.sectional_prop_section = ''
        self.grid_section = ''
        self.symmetry_section = ''

        self.msg = ''

        self.log = get_logger(log, 'debug' if debug else 'info')

    def write_plot3d(self, p3dname, is_binary=False, is_iblank=False):
        assert is_binary == False
        assert is_iblank == False

        with open(p3dname, 'w') as p3d_file:
            npatches = len(self.patches)
            npatches = 1
            msg = '%i\n' % npatches
            for patch_id, patch in sorted(iteritems(self.patches)):
                #if patchID == 1:
                print("patch_id = %s" % patch_id)
                ni, nj = patch.x.shape
                #nr = patch.nrows
                #nc = patch.ncols
                #msg += '%i %i 1\n' % (nc, nr)
                msg += '%i %i 1\n' % (ni, nj-1)
                break

            p3d_file.write(msg)
            for patch_id, patch in sorted(iteritems(self.patches)):
                #if patch_id == 1:
                patch.write_plot3d(p3d_file, 1) # x
                patch.write_plot3d(p3d_file, 2) # y
                patch.write_plot3d(p3d_file, 3) # z
                p3d_file.write('\n')
                break

    def print_file(self):
        msg = ''
        for i, line in enumerate(self.lines):
            msg += "%6s %s" % (i + 1, line)
        msg += '                                      record of input processing\n\n\n'
        return msg

    #def get_panel_points(self, ipanel):
        #r = ipanel % (self.nrows - 1)
        #c = ipanel / (self.nrows - 1)

        ##print "r=%s c=%s" %(r, c)
        #p1 = self.get_point(r, c)
        #p2 = self.get_point(r, c + 1)
        #p3 = self.get_point(r + 1, c + 1)
        #p4 = self.get_point(r + 1, c)
        #return (p1, p2, p3, p4)

    #def get_panel_point_IDs(self, ipanel):
        #r = ipanel % (self.nrows - 1)
        #c = ipanel / (self.nrows - 1)

        ##print "r=%s c=%s" %(r, c)
        #p1 = self.get_point_ID(r, c)
        #p2 = self.get_point_ID(r, c + 1)
        #p3 = self.get_point_ID(r + 1, c + 1)
        #p4 = self.get_point_ID(r + 1, c)
        #return (p1, p2, p3, p4)

    def nPanels(self):
        total_npanels = 0
        for patch_id in range(self.nPatches()):
            patch = self.get_patch(patch_id)
            if patch.kt == 1:
                total_npanels += patch.nPanels()
                #self.log.debug("nPanels = %s" % (total_npanels))
        return total_npanels

    def nPatches(self):
        return len(self.patches)

    def get_patch(self, patch_id):
        return self.patches[patch_id]

    def update_cases(self):
        """
        reduces confusion by only printing cases that will run
        """
        if len(self.alphas) > self.ncases:
            print("self.ncases =", self.ncases)
            self.alphas = self.alphas[:self.ncases]
        if len(self.betas) > self.ncases:
            self.betas = self.alphas[:self.ncases]

    def write_panair(self, out_filename):
        self.update_cases()
        if PY2:
            wb = 'wb'
        else:
            wb = 'w'

        with open(out_filename, wb) as outfile:
            outfile.write(self.title_section)
            outfile.write(self.write_data_check())
            outfile.write(self.symmetry_section)

            outfile.write(self.write_mach())
            outfile.write(self.write_cases())
            outfile.write(self.write_alphas())
            outfile.write(self.write_reference_quantities())

            #outfile.write(self.alphaSection)
            #outfile.write(self.caseSection)
            for patch_name, patch in sorted(iteritems(self.patches)):
                outfile.write(str(patch))

            outfile.write(self.xyz_section)
            outfile.write(self.streamline_section)
            outfile.write(self.flow_section)
            outfile.write(self.sectional_prop_section)
            outfile.write(self.grid_section)
            outfile.write(self.write_printout())

            outfile.write(self.pea_section)
            outfile.write(self.write_liberalized_abutments())
            outfile.write(self.write_end())

    def remove_comments(self, lines):
        lines2 = []
        for line in lines:
            line = line.rstrip().lower()
            if '=' in line:
                #print "line -> |%s|" %(line)
                if '=' is not line[0]:
                    self.log.debug("line[0] -> %s" % line[0])
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
        self.title_section = '\n'.join(section) + '\n'
        self.title_lines = section[1:]
        return True

    def _read_data_check(self, section):
        self.data_check = int(float(section[1][0:10]))
        #self.data_check = 2
        return True

    def _read_symmetry(self, section):
        """
        .. code-block:: console

          $symmetry - xz plane of symmetry
          =misymm   mjsymm
          1.        0.

        .. warning::

            doesnt consider antisymmetryic
        """
        # doesnt consider antisymmetric
        self.xz_symmetry = int(float(section[1][0:10]))
        self.xy_symmetry = int(float(section[1][10:20]))   # ???

        # doesnt consider antisymmetric        self.xy_symmetry = int(float(section[1][10:20]))
        self.nsymmetry_planes = self.xz_symmetry + self.xy_symmetry
        self.symmetry_section = '\n'.join(section) + '\n'
        self.log.debug("xz_symmetry=%s xy_symmetry=%s" % (
            self.xz_symmetry, self.xy_symmetry))
        return True

    def split_points(self, lines, nactual, nremainder):
        """
        reads the points
        """
        points = []
        for n in range(nactual):
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

        if nremainder:
            n += 1
            #print "***"
            line = lines[n]
            (x1, y1, z1) = float(line[0:10]), float(line[10:20]), float(line[20:30])
            point1 = array([x1, y1, z1])
            #print list_print(point1)
            points.append(point1)
        #print "points = ",list_print(points)
        return points

    def add_wake_patch(self, network_name, options, x, y, z):
        #print "self.nnetworks = ",self.nnetworks
        patch = PanairWakePatch(self.nnetworks, network_name, options, x, y, z, self.log)
        self.msg += patch.process()
        #print "get_patch = ",get_patch
        self.patches[patch.inetwork] = patch  # deepcopy?
        self.nnetworks += 1

    def add_patch(self, network_name, kt, cp_norm, x, y, z):
        #print "self.nnetworks = ",self.nnetworks
        patch = PanairPatch(self.nnetworks, network_name, kt, cp_norm, x, y, z, self.log)
        self.msg += patch.process()
        #print "get_patch = ",get_patch
        self.patches[patch.inetwork] = patch  # deepcopy?
        self.nnetworks += 1

    def find_patch_by_name(self, network_name):
        names = []
        for patch_id, patch in iteritems(self.patches):
            #self.log.debug("patchID=%s" % (patchID))
            #self.log.debug("*get_patch = %s" %(get_patch))
            #self.log.debug("get_patch.network_name=%s" % (get_patch.network_name))
            if patch.network_name == network_name:
                return patch
            names.append(patch.network_name)
        msg = 'couldnt findPatchbyName name=%r names=%s' % (network_name, names)
        raise KeyError(msg)

    def _read_points(self, section):
        """
        .. code-block:: console

          $points - wing-body  with composite panels
          =kn                                               cpnorm
          4.                                                2.
          =kt
          1.
          =nm       nn                                                          netname
          11.       3.                                                          winga
          =x(1,1)   y(1,1)    z(1,1)    x(*,*)    y(*,*)    z(*,*)
             69.4737    9.2105    0.0000   63.7818    9.5807    0.7831
        """
        nnetworks = int(float(section[1][0:10]))
        cp_norm = section[1][50:60].strip()
        #cp_norm = float(section[1][50:60])

        kt = int(float(section[2][0:10]))
        #if cpNorm:
            #self.log.debug("nnetworks=%s cpNorm=%s" % (nnetworks, cpNorm))
        #else:
            #self.log.debug("nnetworks=%s" % (nnetworks))

        n = 4
        self.msg += '      kn,kt            %i          %i\n' % (nnetworks, kt)

        for inetwork in range(nnetworks):
            #self.log.debug("lines[* %s] = %s" % (n - 1, section[n - 1]))
            nm = int(float(section[n - 1][0:10]))
            nn = int(float(section[n - 1][10:20]))
            network_name = section[n - 1][70:80]
            #self.log.debug("kt=%s nm=%s nn=%s netname=%s" % (
                #kt, nm, nn, network_name))

            x = zeros([nm, nn])
            y = zeros([nm, nn])
            z = zeros([nm, nn])
            nfull_lines = nm / 2
            npartial_lines = nm % 2
            nlines = nfull_lines + npartial_lines
            #print "nfull_lines=%s npartial_lines=%s nlines=%s" %(nfull_lines, npartial_lines, nlines)
            for j in range(nn):
                lines = section[n:n + nlines]
                n += nlines
                #print '\n'.join(lines)
                points = self.split_points(lines, nfull_lines, npartial_lines)

                for i, point in enumerate(points):
                    x[i][j] = point[0]
                    y[i][j] = point[1]
                    z[i][j] = point[2]

            #print "--X--"
            #print x
            self.add_patch(network_name, kt, cp_norm, x, y, z)
            n += 1
        return True

    def _read_circular_section(self, section):
        """
        .. code-block:: console
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
        """
        nnetworks = int(float(section[1][0:10]))
        cp_norm = section[1][50:60].strip()

        kt = int(float(section[2][0:10]))
        #if cpNorm:
            #self.log.debug("nnetworks=%s cpNorm=%s" % (nnetworks, cpNorm))
        #else:
            #self.log.debug("nnetworks=%s" % (nnetworks))

        n = 3
        self.msg += '      kn,kt            %i          %i\n' % (nnetworks, kt)

        for inetwork in range(nnetworks):
            #self.log.debug("lines[%s] = %s" % (n, section[n]))
            # 0-noDisplacement; 1-Specify
            ndisplacement = int(float(section[n][0:10]))
            assert ndisplacement in [0, 1], section[n]
            n += 1

            #print("\nsection[n] = ", section[n].strip())
            nxr = int(float(section[n][0:10]))
            assert nxr > 0, section[n]
            #print("nxr = %s" % nxr)
            network_name = section[n][70:80]
            #self.log.debug("kt=%s nxr=%s ndisplacement=%s netname=%s" %
                           #(kt, nxr, ndisplacement, netName))

            nfull_lines = nxr / 3
            npartial_lines = int(ceil(nxr % 3 / 3.))
            nxr_lines = npartial_lines + npartial_lines
            #print "X,Rad - npartial_lines=%s npartial_lines=%s nLines=%s\n" % (
                #npartial_lines, npartial_lines, nxr_lines)
            n += 1

            Xin = []
            R = []
            lines = section[n:n + nxr_lines]
            n += nxr_lines

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
            assert nxr == len(Xin), 'nxr=%s Xin=%s' % (nxr, len(Xin))

            #----------------------------------------------------
            #print("section[n] = ", section[n].strip())
            ntheta = int(float(section[n][0:10]))
            assert ntheta > 0, section[n]
            #print("ntheta = %s" % ntheta)
            nfull_lines = ntheta / 6
            npartial_lines = int(ceil(ntheta % 6 / 6.))
            ntheta_lines = nfull_lines + npartial_lines
            #print("Theta - nFullLines=%s nPartialLines=%s nLines=%s" % (
            #    nFullLines, nPartialLines, ntheta_lines))
            n += 1

            lines = section[n:n + ntheta_lines]
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
            n += ntheta_lines

            assert ntheta == len(theta)
            sin_theta_r = [sin(radians(thetai)) for thetai in theta]
            cos_theta_r = [cos(radians(thetai)) for thetai in theta]

            zi = 0.
            X = zeros([nxr, ntheta])
            Y = zeros([nxr, ntheta])
            Z = zeros([nxr, ntheta])
            #print("Xin=%s \nR=%s" % (Xin, R))
            #print("X.shape = ", X.shape)
            for (i, x, r) in zip(count(), Xin, R):
                for (j, sin_theta, cos_theta) in zip(count(), sin_theta_r, cos_theta_r):
                    y = r * sin_theta
                    z = r * cos_theta + zi
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
            self.add_patch(network_name, kt, cp_norm, X, Y, Z)
        return True

    def _read_grid(self, section):
        self.grid_section = '\n'.join(section) + '\n'
        return True

    def _read_xyz(self, section):
        self.xyz_section = '\n'.join(section) + '\n'
        return True

    def _read_streamlines(self, section):
        self.streamline_section = '\n'.join(section) + '\n'
        return True

    def _read_printout(self, section):
        """
        .. code-block:: console

          isings  igeomp  isingp  icontp  ibconp  iedgep
          ipraic  nexdgn  ioutpr  ifmcpr  icostp
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
        .. code-block:: console

          $trailing wakes from body
          =kn                                               cpnorm
          2.
          =kt       matcw
          18.       1.
          =inat     insd      xwake     twake                                   netname
          bodyl     3.        100.      .0                                      bodylwk
          bodyu     3.        100.      .0                                      bodyuwk
        """
        #return True  # disable wakes
        nnetworks = int(float(section[1][0:10]))
        cp_norm = section[1][50:60].strip()
        if cp_norm:
            cp_norm = float(cp_norm)
        else:
            cp_norm = 0

        #print "section[2] = ",section[2]
        kt = int(float(section[2][0:10]))

        matchw = section[2][10:20].strip()
        if matchw:
            matchw = float(matchw)
        #if cp_norm:
            #print "nnetworks=%s cp_norm=%s" % (nnetworks, cp_norm)
        #else:
            #print "nnetworks=%s" % (nnetworks)
        #print ""
        #matcw = int(float(section[2][10:20]))

        n = 3  # line number
        self.msg += '      kn,kt            %i          %i\n' % (nnetworks, kt)
        #self.log.debug('kt=%s cp_norm=%s matchw=%s' % (kt, cp_norm, matchw))
        for inetwork in range(nnetworks):
            #self.log.debug("lines[* %s] = %s" % (n, section[n]))

            trailed_panel = section[n][0:10].strip()
            edge_number = int(float(section[n][10:20]))
            xWake = float(section[n][20:30])  # x distance
            tWake = float(section[n][30:40])  # 0-wake parallel to x axis;  1-wake in direction of compressibility
            network_name = section[n][70:80].strip()
            self.log.debug('trailed_panel=%s edge_number=%s xWake=%s tWake=%s network_name=%s' % (
                trailed_panel, edge_number, xWake, tWake, network_name))
            try:
                patch = self.find_patch_by_name(trailed_panel)
            except KeyError:
                self.log.debug('trailed_panel isnt defined...trailed_panel=|%s|' % (trailed_panel))
                raise

            #xPoints = patch.x
            #yPoints = patch.y
            #print "xPoints = ",xPoints
            #print "yPoints = ",yPoints
            (p1, x1, y1, z1) = patch.get_edge(edge_number)

            npoints = len(x1)
            X = zeros([2, npoints])
            Y = zeros([2, npoints])
            Z = zeros([2, npoints])

            X[0][:] = x1
            Y[0][:] = y1
            Z[0][:] = z1

            if tWake == 0.:
                X[1][:] = ones([npoints]) * xWake
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
            options = [kt, cp_norm, matchw, trailed_panel, edge_number, xWake,
                       tWake]
            patch = self.add_wake_patch(network_name, options, X, Y, Z)
            print('----------------------------')
            n += 1
        return True

    def _read_partial_edge_abutments(self, section):
        self.pea_section = '\n'.join(section) + '\n'
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
        self.sectional_prop_section = '\n'.join(section) + '\n'
        return True

    def _read_flowfield_properties(self, section):
        self.flow_section = '\n'.join(section) + '\n'
        return True

    #def get(self,section):
    #    pass

    #def get(self,section):
    #    pass

    #def get(self,section):
    #    pass

    def group_sections(self, sections, section_names):
        #self.Points = []
        #self.Streamlines = []
        #self.Trailing = []

        #self.title = None
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

        section_map = {
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
        #valid_maps = ['tit','ref','dat','sym','mac','cas','ang','poi','tra','end']
        #valid_maps += ['xyz','str','flo','sec','pri','pea','eat','gri']
        #valid_maps = ['yaw','poi','tra']
        valid_maps = section_map.keys()

        for section, section_name in zip(sections, section_names):  # 1st line
            self.msg += '  $%s\n' % (section_name)
            #print "section = ",len(section)
            #self.log.debug("section_name=%s" % (section_name))
            if section_name in valid_maps:
                #self.log.debug("section[0] = %s" % (section[0]))
                function_map = section_map[section_name]
                ran = function_map(section)
                assert ran == True, '%s didnt run' % (section_name)
                #self.log.debug("")
            if section_name == 'end':
                break
        return sections

    def split_into_sections(self, lines):
        sections = []
        section_names = []
        section = None
        section_name = None
        for line in lines:
            if '$' in line:
                sections.append(section)
                section_names.append(section_name)
                section = []
                section_name = line[1:4]
                if section_name == 'end':
                    break
                #print "new section...name=%s line - %s" %(section_name, line)
            section.append(line)

        # end section
        section.append(line)
        sections.append(section)
        section_names.append(section_name)

        # get rid of first dummy section
        sections = sections[1:]
        section_names = section_names[1:]

        assert len(sections) == len(section_names), "%s %s" % (
            len(sections), len(section_names))

        #for section in sections: # 1st line
            #print "section = ",len(section)
            #print "section[0] = ",section[0]
        return (sections, section_names)

    def read_panair(self, infilename):
        assert os.path.exists(infilename), print_bad_path(infilename)
        self.infilename = infilename

        with open(self.infilename, 'r') as infile:
            self.lines = infile.readlines()

        lines = self.remove_comments(self.lines)
        (sections, section_names) = self.split_into_sections(lines)
        groups = self.group_sections(sections, section_names)
        #self.log.debug("nPatches = %s" % (self.nPatches()))
        # split into headings
        #for panel in panels:
        #    points = readPoints()
        #print "self.msg = ",self.msg

    def get_points_elements_regions(self, get_wakes=False):
        points = []
        elements = []
        regions = []
        pointI = 0
        for name, panel in sorted(iteritems(self.patches)):
            if not get_wakes:
                if panel.is_wake():
                    continue
            #panel = self.patches[2]
            #self.log.debug("size(X) = %s" %( str( panel.x.shape ) ))
            (pointsI, pointi) = panel.get_points()
            (elementsI) = panel.get_elements(pointI)
            #print 'panel.inetwork=%r' % (panel.inetwork + 1)
            #regions += panel.inetwork * ones(len((elementsI)), 'int32')
            regions += [panel.inetwork + 1] * len(elementsI)

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
        points = array(points, dtype='float32')
        elements = array(elements, dtype='int32')
        regions = array(regions, dtype='int32')
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
        msg += '%s.\n' % (self.data_check)
        return msg

    def write_printout(self):
        msg = '$printout options\n'
        msg += "%-10s%-10s%-10s%-10s%-10s%-10s\n" % (self.isings, self.igeomp,
                                                     self.isingp, self.icontp, self.ibconp, self.iedgep)
        msg += "%-10s%-10s%-10s%-10s%-10s\n" % (self.ipraic,
                                                self.nexdgn, self.ioutpr, self.ifmcpr, self.icostp)
        return msg


def main():
    import sys
    in_filename = sys.argv[1]
    out_filename = in_filename[:-4] + '_new.inp'
    print("outfile_name = %r" % out_filename)
    grid = PanairGrid()
    grid.read_panair(in_filename)
    grid.write_panair(out_filename)
    (points, elements, regions) = grid.get_points_elements_regions()
    print("\nin_filename=%r" % in_filename)


if __name__ == '__main__':  # pragma: no cover
    main()
