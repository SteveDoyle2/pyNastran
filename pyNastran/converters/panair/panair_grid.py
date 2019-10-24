"""
defines:
 - PanairGrid(log=None, debug=False)

"""
import os
from itertools import count
from math import ceil, sin, cos, radians

import numpy as np


from pyNastran.converters.panair.panair_grid_patch import (
    PanairPatch, PanairWakePatch, print_float)
from pyNastran.converters.panair.assign_type import (
    integer, double, integer_or_blank, double_or_blank, fortran_value)
from cpylog import get_logger2
from pyNastran.utils import check_path

#from pyNastran.utils import list_print

#CL = -Fx*sin(alpha)*cos(beta) + Fy*sin(alpha)*sin(beta) +Fz*cos(alpha)
#CD =  Fx*cos(alpha)*cos(beta) - Fy*cos(alpha)*sin(beta) +Fz*sin(alpha)
#CY =  Fx*sin(beta) +Fy*cos(beta)


class PanairGrid:
    """defines the PanairGrid class"""
    model_type = 'panair'

    def __init__(self, log=None, debug=True):
        """
        Initializes the PanairGrid object

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
        self.infilename = None
        self.title_lines = []
        self.xz_symmetry = 0
        self.xy_symmetry = 0
        self.nsymmetry_planes = 0
        self.title = ''

        self.nnetworks = 0
        self.patches = {}
        self.lines = []

        self.alphas = [0.]
        self.ncases = 0
        self.betas = [0.]
        self.alpha_compressibility = 0.
        self.beta_compressibility = 0.

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
        self.is_end = None
        self.pea_section = ''
        self.mach = 0.0
        self.data_check = 2

        self.noff_body_points = 0
        self.nstreamlines = 0

        self.title_section = ''
        self.xyz_section = ''
        self.streamline_section = ''
        self.flow_section = ''
        self.sectional_prop_section = ''
        self.grid_section = ''
        self.symmetry_section = ''

        self.msg = ''
        self.log = get_logger2(log, debug=debug)

    def write_plot3d(self, p3dname, is_binary=False, is_iblank=False):
        assert not is_binary, is_binary
        assert not is_iblank, is_iblank

        with open(p3dname, 'w') as p3d_file:
            npatches = len(self.patches)
            npatches = 1
            msg = '%i\n' % npatches
            for patch_id, patch in sorted(self.patches.items()):
                #if patchID == 1:
                print("patch_id = %s" % patch_id)
                ni, nj = patch.xyz.shape[:2]
                #nr = patch.nrows
                #nc = patch.ncols
                #msg += '%i %i 1\n' % (nc, nr)
                msg += '%i %i 1\n' % (ni, nj-1)
                break

            p3d_file.write(msg)
            for patch_id, patch in sorted(self.patches.items()):
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

    @property
    def npanels(self):
        total_npanels = 0
        for patch_id in range(self.npatches):
            patch = self.get_patch(patch_id)
            if patch.kt == 1:
                total_npanels += patch.npanels
                #self.log.debug("nPanels = %s" % (total_npanels))
        return total_npanels

    @property
    def npatches(self):
        return len(self.patches)

    def get_patch(self, patch_id):
        return self.patches[patch_id]

    def update_cases(self):
        """
        reduces confusion by only printing cases that will run
        """
        assert self.ncases > -1, 'ncases=%s' % self.ncases
        if len(self.alphas) > self.ncases:
            self.alphas = self.alphas[:self.ncases]
        if len(self.betas) > self.ncases:
            self.betas = self.alphas[:self.ncases]

    def write_panair(self, out_filename):
        """writes the panair file"""
        self.update_cases()
        with open(out_filename, 'w') as panair_file:
            panair_file.write(self.title_section)
            panair_file.write(self.write_data_check())
            panair_file.write(self.symmetry_section)

            panair_file.write(self.write_mach())
            panair_file.write(self.write_cases())
            panair_file.write(self.write_alphas())
            panair_file.write(self.write_reference_quantities())

            #panair_file.write(self.alphaSection)
            #panair_file.write(self.caseSection)
            for unused_patch_name, patch in sorted(self.patches.items()):
                panair_file.write(str(patch))

            panair_file.write(self.xyz_section)
            panair_file.write(self.streamline_section)
            panair_file.write(self.flow_section)
            panair_file.write(self.sectional_prop_section)
            panair_file.write(self.grid_section)
            panair_file.write(self.write_printout())

            panair_file.write(self.pea_section)
            panair_file.write(self.write_liberalized_abutments())
            panair_file.write(self.write_end())

    def _read_title(self, section):
        #self.title = section[1:]
        self.title_section = '\n'.join(section) + '\n'
        self.title_lines = section[1:]
        return True

    def _read_data_check(self, section):
        self.data_check = integer(section[1][0:10], 'data_check')
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
        self.xz_symmetry = integer(section[1][0:10], 'xz_symmetry')
        self.xy_symmetry = integer(section[1][10:20], 'xy_symmetry')

        # doesnt consider antisymmetric
        #self.xy_symmetry = integer(section[1][10:20], 'xy_symmetry')
        self.nsymmetry_planes = self.xz_symmetry + self.xy_symmetry
        self.symmetry_section = '\n'.join(section) + '\n'
        self.log.debug("xz_symmetry=%s xy_symmetry=%s" % (
            self.xz_symmetry, self.xy_symmetry))
        return True

    def split_points(self, lines, nfull_lines, npartial_lines):
        """
        reads the points
        """
        ipoint = 0
        nrows = 2 * nfull_lines + npartial_lines
        points = np.zeros((nrows, 3), dtype='float32')

        for n in range(nfull_lines):
            line = lines[n]
            x1 = double(line[0:10], 'x%i' % (ipoint + 1))
            y1 = double(line[10:20], 'y%i' % (ipoint + 1))
            z1 = double(line[20:30], 'z%i' % (ipoint + 1))
            points[ipoint, :] = [x1, y1, z1]
            ipoint += 1

            x2 = double(line[30:40], 'x%i' % (ipoint + 2))
            y2 = double(line[40:50], 'y%i' % (ipoint + 2))
            z2 = double(line[50:60], 'z%i' % (ipoint + 2))
            points[ipoint, :] = [x2, y2, z2]
            ipoint += 1

        if npartial_lines:
            n += 1
            line = lines[n]
            x1 = double(line[0:10], 'x%i' % (ipoint + 1))
            y1 = double(line[10:20], 'y%i' % (ipoint + 1))
            z1 = double(line[20:30], 'z%i' % (ipoint + 1))
            points[ipoint, :] = [x1, y1, z1]

        return points

    def add_wake_patch(self, network_name, options, xyz):
        patch = PanairWakePatch(self.nnetworks, network_name, options, xyz, self.log)
        self.msg += patch.process()
        self.patches[patch.inetwork] = patch  # deepcopy?
        self.nnetworks += 1

    def add_patch(self, network_name, kt, cp_norm, xyz):
        patch = PanairPatch(self.nnetworks, network_name, kt, cp_norm, xyz, self.log)
        self.msg += patch.process()
        self.patches[patch.inetwork] = patch  # deepcopy?
        self.nnetworks += 1

    def find_patch_by_name(self, network_name):
        names = []
        for unused_patch_id, patch in self.patches.items():
            #self.log.debug("patch_id=%s" % (patch_id))
            #self.log.debug("*get_patch = %s" %(get_patch))
            #self.log.debug("get_patch.network_name=%s" % (get_patch.network_name))
            if patch.network_name == network_name:
                return patch
            names.append(patch.network_name)
        msg = "couldn't find_patch_by_name name=%r names=%s" % (network_name, names)
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
        nnetworks = integer(section[1][0:10], 'nnetworks')
        cp_norm = section[1][50:60].strip()
        #cp_norm = double(section[1][50:60], 'cp_norm')

        kt = integer(section[2][0:10], 'kt')
        #if cp_norm:
            #self.log.debug("nnetworks=%s cp_norm=%s" % (nnetworks, cp_norm))
        #else:
            #self.log.debug("nnetworks=%s" % (nnetworks))

        n = 4
        self.msg += '      kn,kt            %i          %i\n' % (nnetworks, kt)

        for unused_inetwork in range(nnetworks):
            #self.log.debug("lines[* %s] = %s" % (n - 1, section[n - 1]))
            nm = integer(section[n - 1][0:10], 'nm')
            nn = integer(section[n - 1][10:20], 'nn')
            network_name = section[n - 1][70:80]
            #self.log.debug("kt=%s nm=%s nn=%s netname=%s" % (
                #kt, nm, nn, network_name))

            xyz = np.zeros([nm, nn, 3], dtype='float32')
            #y = zeros([nm, nn])
            #z = zeros([nm, nn])
            nfull_lines = nm // 2
            npartial_lines = nm % 2
            nlines = nfull_lines + npartial_lines
            #print "nfull_lines=%s npartial_lines=%s nlines=%s" % (
                #nfull_lines, npartial_lines, nlines)
            for j in range(nn):
                lines = section[n:n + nlines]
                n += nlines
                #print '\n'.join(lines)
                points = self.split_points(lines, nfull_lines, npartial_lines)

                for i, point in enumerate(points):
                    xyz[i, j, :] = [point[0], point[1], point[2]]

            self.add_patch(network_name, kt, cp_norm, xyz)
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
              0.0214    2.4542   -0.0200    2.5485    0.0388    2.6522 # 18
              0.2056    2.7554    0.4869    2.8522    0.8883    2.9413
              1.4250    3.0178    2.1188    3.0656    2.9586    3.0658 # 36
              3.8551    3.0175    4.6715    2.9439    5.3492    2.8700 # 42
              6.0000    2.7842    6.4687    2.7442 # 46 ->23
          =nn
          5.
          =th(1)    th(2)     th(3)     th(4)     th(5)
          -90.      -45.      0.        45.       90.

        """
        nnetworks = integer(section[1][0:10], 'nnetworks')
        cp_norm = section[1][50:60].strip()

        kt = integer(section[2][0:10], 'kt')
        if cp_norm:
            self.log.debug("nnetworks=%s cp_norm=%s" % (nnetworks, cp_norm))
        else:
            self.log.debug("nnetworks=%s" % (nnetworks))
        #print("nnetworks=%s cp_norm=%s" % (nnetworks, cp_norm))

        n = 3
        self.msg += '      kn,kt            %i          %i\n' % (nnetworks, kt)

        for unused_inetwork in range(nnetworks):
            #self.log.debug("lines[%s] = %s" % (n, section[n]))
            # 0-noDisplacement; 1-Specify
            ndisplacement = integer(section[n][0:10], 'ndisplacement')
            assert ndisplacement in [0, 1], section[n]
            n += 1

            #print("\nsection[n] = ", section[n].strip())
            nx_nr = integer(section[n][0:10], 'nx_nr')
            assert nx_nr > 0, section[n]
            #print("nxr = %s" % nxr)
            network_name = section[n][70:80]
            #self.log.debug("kt=%s nx_nr=%s ndisplacement=%s network_name=%s" %
                           #(kt, nx_nr, ndisplacement, network_name))

            nfull_lines = nx_nr // 3
            npartial_lines = int(ceil(nx_nr % 3 / 3.))
            nx_nr_lines = nfull_lines + npartial_lines
            #print("X,Rad - nfull_lines=%s npartial_lines=%s nx_nr_lines=%s\n" % (
                #nfull_lines, npartial_lines, nx_nr_lines))
            n += 1

            Xin = []
            R = []
            lines_full = section[n:n + nx_nr_lines]
            n += nx_nr_lines

            assert len(lines_full) == 7, len(lines_full)

            ipoint = 1
            for line in lines_full:
                if len(line[0:60].strip()) > 0:
                    x1 = double(line[0:10], 'x%i' % ipoint)
                    r1 = double(line[10:20], 'r%i' % ipoint)
                    Xin.append(x1)
                    R.append(r1)
                    ipoint += 1
                else:  # pragma: no cover
                    raise RuntimeError(line.rstrip())

                if len(line[20:60].strip()) > 0:
                    x2 = double(line[20:30], 'x%i' % ipoint)
                    r2 = double(line[30:40], 'r%i' % ipoint)
                    Xin.append(x2)
                    R.append(r2)
                    ipoint += 1

                if len(line[40:60].strip()) > 0:
                    x3 = double(line[40:50], 'x%i' % ipoint)
                    r3 = double(line[50:60], 'r%i' % ipoint)
                    Xin.append(x3)
                    R.append(r3)
                    ipoint += 1
            assert nx_nr == len(Xin), 'nx,nr=%s Xin=%s len(Xin)=%s' % (nx_nr, Xin, len(Xin))

            #----------------------------------------------------
            #print("section[n] = ", section[n].strip())
            ntheta = integer(section[n][0:10], 'ntheta')
            assert ntheta > 0, section[n]
            #print("ntheta = %s" % ntheta)
            nfull_lines = ntheta // 6
            npartial_lines = int(ceil(ntheta % 6 / 6.))
            ntheta_lines = nfull_lines + npartial_lines
            #print("Theta - nfull_lines=%s npartial_lines=%s nlines=%s" % (
            #    nfull_lines, npartial_lines, ntheta_lines))
            n += 1

            lines = section[n:n + ntheta_lines]
            theta = []

            ipoint = 1
            for line in lines:
                #print(line)
                if len(line[0:60].rstrip()) > 0:
                    theta1 = double(line[0:10], 'theta%i' % ipoint)
                    theta.append(theta1)
                    ipoint += 1
                else:  # pragma: no cover
                    raise RuntimeError(line.rstrip())

                if len(line[0:60].rstrip()) > 0:
                    theta2 = double(line[10:20], 'theta%i' % ipoint)
                    theta.append(theta2)
                    ipoint += 1

                if len(line[20:60].rstrip()) > 0:
                    theta3 = double(line[20:30], 'theta%i' % ipoint)
                    theta.append(theta3)
                    ipoint += 1

                if len(line[30:60].rstrip()) > 0:
                    theta4 = double(line[30:40], 'theta%i' % ipoint)
                    theta.append(theta4)
                    ipoint += 1

                if len(line[40:60].rstrip()) > 0:
                    theta5 = double(line[40:50], 'theta%i' % ipoint)
                    theta.append(theta5)
                    ipoint += 1

                if len(line[50:60].rstrip()) > 0:
                    theta6 = double(line[50:60], 'theta%i' % ipoint)
                    theta.append(theta6)
                    ipoint += 1
            n += ntheta_lines

            assert ntheta == len(theta)
            sin_theta_r = [sin(radians(thetai)) for thetai in theta]
            cos_theta_r = [cos(radians(thetai)) for thetai in theta]

            zi = 0.
            # TODO: can't we use numpy meshgrid?
            XYZ = np.zeros([nx_nr, ntheta, 3], dtype='float32')
            #print("Xin=%s \nR=%s" % (Xin, R))
            for (i, x, r) in zip(count(), Xin, R):
                for (j, sin_theta, cos_theta) in zip(count(), sin_theta_r, cos_theta_r):
                    y = r * sin_theta
                    z = r * cos_theta + zi
                    #print("x=%s y=%s z=%s" % (x, y, z))
                    XYZ[i, j, :] = [x, y, z]

            #print("--XYZ--")
            #print(xyz)
            self.add_patch(network_name, kt, cp_norm, XYZ)
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
        self.isings = double(section[1][0:10], 'isings')
        #self.isings = 5.
        #2.0 output wake singularity grid-singularity strength and gradients at 9 pts/panel
        #3.0 add table of singularity parameter values
        #4.0 add network array of singularity values
        #5.0 add singularity grid for all network

        self.igeomp = double(section[1][10:20], 'igeomp')
        self.igeomp = 1.
        # 1.0 outputs panel data - geometry diagnostic data

        self.isingp = double(section[1][20:30], 'isingp')
        #self.isignp = 1.  # ???
        # 1.0 outputs singularity spline data

        self.icontp = double(section[1][30:40], 1.0)
        self.icontp = 1.0
        # 1.0 outputs control pt location, upper surface unit normals, control point
        #     diagnositc data, control pt maps and singularity parameter maps
        # 2.0 add additional control point map for each network

        self.ibconp = double(section[1][40:50], 'ibconp')
        #self.ibconp = 0.0
        # 1.0 outputs boundary condition coeffs, diagnositc data,
        #     boundary condition maps, control point maps,
        #     and singularity paramter maps

        self.iedgep = double(section[1][50:60], 'iedgep')
        self.iedgep = 0.0
        # 1.0 outputs edge-matching diagnositic data
        #print("igeomp=%r isingp=%r icontp=%r ibconp=%r iedgep=%r" % (
            #igeomp, isingp, icontp, ibconp, iedgep))

        #  ROW 2
        self.ipraic = double(section[2][0:10], 'ipraic')
        #self.ipraic = 3.
        # xxx inputs control point sequence number for
        #     which AIC values are to be pritned (rarely used)

        self.nexdgn = double(section[2][10:20], 'nexdgn')
        self.nexdgn = 0.
        # 1.0 outputs edge and extra control point data

        self.ioutpr = double(section[2][20:30], 'ioutpr')
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
        #print("ipraic=%r nexdgn=%r ioutpr=%r ifmcpr=%r ifmcpr=%r"
              #% (ipraic, nexdgn, ioutpr, ifmcpr, iedgep))
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
        nnetworks = integer_or_blank(section[1][0:10], 'nnetworks', 0)

        cp_norm = integer_or_blank(section[1][50:60], 'cp_norm', 0)

        #print "section[2] = ",section[2]
        kt = integer(section[2][0:10], 'kt')

        matchw = section[2][10:20].strip()
        if matchw:
            matchw = float(matchw)
        #if cp_norm:
            #print "nnetworks=%s cp_norm=%s" % (nnetworks, cp_norm)
        #else:
            #print "nnetworks=%s" % (nnetworks)
        #print ""
        #matcw = integer(section[2][10:20], 'matcw)

        n = 3  # line number
        self.msg += '      kn,kt            %i          %i\n' % (nnetworks, kt)
        #self.log.debug('kt=%s cp_norm=%s matchw=%s' % (kt, cp_norm, matchw))
        for unused_inetwork in range(nnetworks):
            #self.log.debug("lines[* %s] = %s" % (n, section[n]))

            trailed_panel = section[n][0:10].strip()
            edge_number = integer(section[n][10:20], 'edge_number')
            xwake = double(section[n][20:30], 'xwake')  # x distance

            # 0-wake parallel to x axis
            #1-wake in direction of compressibility
            twake = double(section[n][30:40], 'twake')
            network_name = section[n][70:80].strip()
            self.log.debug('trailed_panel=%s edge_number=%s xwake=%s twake=%s network_name=%s' % (
                trailed_panel, edge_number, xwake, twake, network_name))
            try:
                patch = self.find_patch_by_name(trailed_panel)
            except KeyError:
                self.log.debug('trailed_panel isnt defined...trailed_panel=%r' % (trailed_panel))
                raise

            (unused_p1, xyz1) = patch.get_edge(edge_number)

            npoints = xyz1.shape[0]
            xyz = np.zeros([2, npoints, 3])

            xyz[0, :, :] = xyz1
            if twake == 0.:
                xyz[1, :, 0] = np.ones([npoints]) * xwake
                xyz[1, :, 1] = xyz1[:, 1]
                xyz[1, :, 2] = xyz1[:, 2]
            else:
                #alphaC, betaC
                raise NotImplementedError('twake isnt supported')
            #self.log.debug("--xyz---")
            #self.log.debug(xyz)
            #nm = integer(section[n-1][0 :10], 'nm')
            #nn = integer(section[n-1][10:20], 'nn')
            options = [kt, cp_norm, matchw, trailed_panel, edge_number, xwake,
                       twake]
            patch = self.add_wake_patch(network_name, options, xyz)
            self.log.info('----------------------------')
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

        # absolute value of tolerance used to establith equivalent network
        # edge points; minimum panel diagonal times 0.001 <= espgeo <= minimum
        # panel diagonal times 0.03; to not be restricted use negative numbers
        self.epsgeo = double_or_blank(section[1][0:10].strip(), 'epsgeo', 0.0)

        # prints network geometry before $EAT
        #  0.0 prints geometry
        # -1.0 omits printout
        self.igeoin = double_or_blank(section[1][10:20].strip(), 'igeoin', -1.)

        # prints network geometry after $EAT
        # 0.0 omits printout
        # 1.0 prints geometry
        self.igeout = double_or_blank(section[1][20:30].strip(), 'igeout', 0.)

        # prints abutment cross-reference; for each network, all abutments and
        # abutment intersections are described in a format similar to that
        # used in the abutment summary and abutment intersection summary output sections
        # 0.0 omits printout
        # 1.0 prints geometry
        self.nwxref = double_or_blank(section[1][30:40].strip(), 'nwxref', 0.)

        # optionally checks for intersection of network edges through other panels
        # for final geometry; used only if needed on complex configurations; can
        # increase the cost of a computer run by 1 or 2 times that of a datachek (2.) run
        # 0.0 omits intersection checking
        # 1.0 checks intersection
        self.triint = double_or_blank(section[1][40:50], 'triint', 0.)

        #  0.0 abutment summary
        #  1.0 plus intersection summary
        #  2.0 plus complete list of each network and associated abutments
        #  3.0 plus diagnostic data for program developer
        # -1.0 no abutment printout
        self.iabsum = double_or_blank(section[1][50:60].strip(), 'iabsum', 0.)
        return True

    def _read_sectional_properties(self, section):
        self.sectional_prop_section = '\n'.join(section) + '\n'
        return True

    def _read_flowfield_properties(self, section):
        self.flow_section = '\n'.join(section) + '\n'
        return True

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
            #print "section = ", len(section)
            #self.log.debug("section_name=%s" % (section_name))
            if section_name in valid_maps:
                #self.log.debug("section[0] = %s" % (section[0]))
                function_map = section_map[section_name]
                ran = function_map(section)
                assert ran, '%s didnt run' % (section_name)
                #self.log.debug("")
            if section_name == 'end':
                break
        return sections

    def read_panair(self, infilename):
        """reads a panair input file"""
        check_path(infilename, 'panair_input_Filename')
        self.infilename = infilename

        with open(self.infilename, 'r') as infile:
            self.lines = infile.readlines()

        lines = remove_comments(self.lines, self.log)
        (sections, section_names) = split_into_sections(lines)
        unused_groups = self.group_sections(sections, section_names)
        #self.log.debug("nPatches = %s" % (self.npatches))
        # split into headings
        #for panel in panels:
        #    points = readPoints()
        #print "self.msg = ",self.msg

    def get_points_elements_regions(self, get_wakes=False):
        points = []
        elements = []
        regions = []
        kt = []
        cp_norm = []
        npoints = 0
        for unused_name, panel in sorted(self.patches.items()):
            if not get_wakes:
                if panel.is_wake():
                    continue

            #self.log.debug("size(xyz) = %s" % ( str(panel.xyz.shape) ))
            patch_points, npointsi = panel.get_points()
            patch_elements = panel.get_elements(npoints) + npoints
            npatch_elements = patch_elements.shape[0]
            regions.append(np.ones(npatch_elements, dtype='int32') * (panel.inetwork + 1))
            kt.append(np.ones(npatch_elements, dtype='int32') * panel.kt)
            #print('panel.cp_norm = %r' % panel.cp_norm)
            cp_norm.append(np.ones(npatch_elements, dtype='int32') * panel.cp_norm)

            points += patch_points
            elements.append(patch_elements)
            npoints += npointsi
            #print "name=%r npointsi=%s len(elements)=%s len(points)=%s" %(
                #name, npointsi, len(elements), len(points))

            #if npoints != 0 and False:
                #if panel.network_name in ['bodylwk', 'awbw']:
                    #for point in patch_points:
                        #print(point)
                    #for element in patch_elements:
                        #(n1, n2, n3, n4) = element
                        #print(element)
                        #print(points[n1], points[n2], points[n3], points[n4])
                        #print("")
        points = np.array(points, dtype='float32')
        elements = np.vstack(elements)
        regions = np.hstack(regions)
        kt = np.hstack(kt)
        cp_norm = np.hstack(cp_norm)
        return points, elements, regions, kt, cp_norm

    def _read_cases(self, section):
        """
        $cases - no. of solutions
        =nacase
        1.
        """
        self.ncases = integer(section[1][0:10], 'ncases')
        #self.log.debug("ncases = %s" % (self.ncases))
        return True

    def _read_mach(self, section):
        """
        $mach number
        =amach
        .6
        """
        self.mach = float(section[1][0:10])
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
            out += '%-10s' % print_float(self.ncases) + '\n'
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
        self.alpha_compressibility = float(section[1][0:10])
        sline = section[2].split()
        self.alphas = [float(slot) for slot in sline]
        #self.log.debug("alpha_compressibility=%s alphas=%s" % (
            #self.alpha_compressibility, self.alphas))
        return True

    def set_alphas(self, alphas, alpha_compressibility):
        self.alpha_compressibility = alpha_compressibility
        self.alphas = alphas
        self.ncases = len(alphas)

    def write_title(self):
        out = '$title\n'
        out += '%s\n' % self.title.strip()
        return out

    def write_alphas(self):
        out = '$angles-of-attack\n'
        out += '%-s\n' % self.alpha_compressibility
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
        self.beta_compressibility = float(section[1][0:10])
        sline = section[2].split()
        self.betas = [float(slot) for slot in sline]
        #self.log.debug("beta_compressibility=%s betas=%s" % (
            #self.beta_compressibility, self.betas))
        return True

    def set_betas(self, betas, beta_compressibility):
        self.beta_compressibility = beta_compressibility
        self.betas = betas
        self.ncases = len(betas)

    def write_betas(self):
        out = '$yaw\n'
        out += '%s\n' % self.beta_compressibility
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

    def _read_end(self, unused_section):
        self.is_end = True
        #self.log.debug("end...")
        return True

    def write_end(self):
        if self.is_end:
            return '$end of panair inputs\n '
        return ''

    def write_data_check(self):
        msg = '$datacheck\n'
        msg += '%s.\n' % (self.data_check)
        return msg

    def write_printout(self):
        msg = '$printout options\n'
        msg += "%-10s%-10s%-10s%-10s%-10s%-10s\n" % (
            self.isings, self.igeomp, self.isingp, self.icontp, self.ibconp, self.iedgep)
        msg += "%-10s%-10s%-10s%-10s%-10s\n" % (
            self.ipraic, self.nexdgn, self.ioutpr, self.ifmcpr, self.icostp)
        return msg

    def print_abutments(self):
        msg = ''
        msg += '               SUMMARY OF FACING SURFACES (+:upper, -:lower)\n'
        msg += ' abutment   nw-ident  ntd  knet.edge    nw-ident  ntd  knet.edge\n'
        msg += '        1   winga      12     1.3+      wingwk     18     7.1+\n'
        msg += '            wingwk     18     7.1-      winga      12     1.1+\n'
        msg += '            winga      12     1.1-      winga      12     1.3-\n'
        msg += '        9   1st p-o-s   0    -1.0       bodyl      12     3.4+\n'
        msg += '            bodyl      12     3.4-      1st p-o-s   0    -1.0\n'

        for patch_id, patch in self.patches.items():
            (p1, unused_xyz1) = patch.get_edges()
            self.log.debug("p[%s] = %s" % (patch_id, p1))
        return msg

    def print_options(self):
        msg = (
            '0            options\n'
            '            %i = singularity grid print flag\n'
            '            %i = panel geometry print flag\n'
            '            %i = spline data flag  ( 0 ==> off, nonzero ==> on )\n'
            '            %i = control point information print flag\n'
            '            %i = boundary condition data print flag \n'
            '            %i = edge matching information print flag\n'
            '            %i = index of control point for which aic-s are printed\n'
            '            %i = edge control point flow properties print flag\n'

            # self.ioutpr, self.ifmcpr, self.icostp, fortran_value(self.epsgeo),
            '            %i = output control flag (-1 ==> no surface flow properties, '
            '0 ==> standard output, 1 ==> short form output )\n'
            '            %s = force/moment control flag (-1 ==> no force and moment data, '
            '0 ==> standard output, 1 ==> nw totals only )\n'
            '            %i = print flag for detailed cost information during execution of job\n'
            '            1 = print flag for singularity parameter maps\n'
            '0               abutment processing options\n'
            '   %10s = global edge abutment tolerance specified by user.  '
            'if this value is zero, a default value will be calculated\n'
            '                later.   this default value is taken as:  '
            '.001  * (minimum panel diameter)\n'

            # self.igeoin, self.igeout, self.nwxref,
            '            %i = print flag controlling geometry printout  b e f o r e  '
            'the abutment processing.''  ( nonzero ==> do print )\n'
            '            %i = print flag controlling geometry printout  a f t e r    '
            'the abutment processing.  ( nonzero ==> do print )\n'
            '            %i = network/abutment/abutment-intersection print flag.  '
            '( nonzero ==> generate the cross referenced abutment listing\n'

            # self.triint, 2
            '            %i = control index for panel intersection checking.  '
            '( nonzero ==> do perform the check. )\n'
            '            %i = abutment/abutment-intersection (short listing) print flag '
            '( 0 ==> suppress, nonzero ==> generate usual print )\n'
            ' \n'

            '                force and moment reference parameters\n'
            '   %10s = reference area for force and moment calculations.    (sref)\n'
            '   %10s = rolling moment reference length  (bref)\n'
            '   %10s = pitching moment reference length (cref)\n'
            '   %10s = yawing moment reference length   (dref)\n'
            '   %10s = x - coordinate for the point about which moments will be calculated  (xref)\n'
            '   %10s = y - coordinate for the point about which moments will be calculated  (yref)\n'
            '   %10s = z - coordinate for the point about which moments will be calculated  (zref)\n'
            '            3 = pressure coefficient index (nprcof) '

            '(1=linear, 2=slenderbody, 3=2nd, 4=isentropic)\n'
            '1\n' % (
                self.isings, self.igeomp, self.isingp, self.icontp, self.ibconp,
                self.iedgep, self.ipraic, self.nexdgn,

                self.ioutpr, self.ifmcpr, self.icostp, fortran_value(self.epsgeo),

                self.igeoin, self.igeout, self.nwxref,

                self.triint, 2,

                fortran_value(self.sref), fortran_value(self.bref), fortran_value(self.cref),
                fortran_value(self.dref),
                fortran_value(self.xref), fortran_value(self.yref), fortran_value(self.zref),
            )
        )
        return msg

    def print_out_header(self):
        msg = """\n
         ****************************************************************************************************

          dynamic memory management initialization

          max no. levels         15   max no. arrays        200   maximum scratch storage     900000   total storage provided  900000
             addr(maplev)          0     addr(maplws)          0     addr(scratch storage)          1

        ****************************************************************************************************
          wopen call on unit    1  blocks:    10   status:     0
          wopen call on unit    2  blocks:    10   status:     0
          wopen call on unit    3  blocks:    10   status:     0
        1


        *****************************************************************************
         *                                                                           *
         *                     a502 - pan-air technology program                     *
         *                                                                           *
         *               potential flow about arbitrary configurations               *
         *               version id = ht2 (12 feb 92) boeing ver i00                 *
         *                                                                           *
         *                                   02/12/92                                *
         *                                                                           *
         *                                                                           *
         * simple wing-body with composite panel. (run with a502i)
         * saaris  865-6209  m/s 7c-36
         *                                                                           *
         *                                                                           *
         *                                                                           *
         *****************************************************************************
        1
        0*b*input-da - list of a502 input data cards -\n"""
        return msg

    def print_grid_summary(self):
        msg = ''
        msg2 = ''
        msg3 = ''

        msg += '        1               ***  quick summary of a502 input  ***\n'

        for i, title_line in enumerate(self.title_lines):
            msg += "title%s:%s\n" % (i + 1, title_line)

        msg += (
            '0               processing options\n'
            '            %i = datacheck.   (0=regular run,1=full datacheck,2=short datacheck)\n'
            '            0 = s.p. flag.   (0 ==> no s.p. file (ft09) provided, 1 ==> '
            'local file ft09 with singularity values is provided)\n'
            '            0 = aic flag.    (0 ==>  no aic file (ft04) provided, 1 ==> '
            'local file ft04 with aic-s is provided by the user)\n'
            '            0 = b.l. flag    (0 ==> no boundary layer file requested, 1 ==> '
            'boundary layer data will be written to file ft17)\n'
            '            0 = velocity correction index.  (0 ==> no correction, 1 ==> '
            'mclean correction, 2 ==>  boctor correction)\n'
            '            0 = flow visualization flag.  (nonzero ==> '
            'off-body and streamline processing will be performed)\n'
            '            0 = off-body calculation type. (0 ==> mass flux, nonzero ==> velocity)\n'
            '            0 = streamline calculation type. (0 ==> mass flux, nonzero ==> velocity)\n'

            '           %2i = number of off-body points.\n'
            '           %2i = number of streamlines to be traced.\n'

            '0               case summary\n'
            '           %2i = number of cases\n'
            '     %f = mach number\n'
            '     %f = compressibility axis angle of attack (alpc)\n'
            '     %f = compressibliity axis angle of sideslip (betc)\n' % (
                self.data_check,
                self.noff_body_points, self.nstreamlines,
                self.ncases, self.mach, self.alpha_compressibility, self.beta_compressibility,
            )
        )

        msg3 += '0network id&index   #rows   #cols  kt  src  dblt  nlopt1  nropt1  nlopt2  nropt2    ipot   # pts  # pans  cpnorm  cum pt  cum pn\n'
        msg3 += '---------- -----   -----   -----  --  ---  ----  ------  ------  ------  ------    ----    ----    ----  ------  ------  ------\n'

        total_points = 0
        total_panels = 0
        for patch_id in range(self.npatches):
            patch = self.get_patch(patch_id)
            msg3 += patch.quick_summary(total_points, total_panels)
            total_panels += patch.npanels
            total_points += patch.npoints

        msg2 += '0  case       alpha          beta      mag(f-s-v)\n'
        msg2 += ' ------    ----------    ----------   -----------\n'
        for icase in range(self.ncases):
            alpha = self.alphas[icase]
            beta = self.betas[icase]
            msg2 += '     %2s      %f      %f      1.000000\n' % (
                icase + 1, alpha, beta)

        msg2 += (
            '0               symmetry options\n'
            '            %s = number of planes of symmetry\n'
            '            %s = x-z plane of symmetry flag (0 ==> no symmetry, '
            '1==> flow symmetry, -1 ==> flow antisymmetry)\n'
            '            %s = x-y plane of symmetry flag (0 ==> no symmetry, '
            '1==> flow symmetry, -1 ==> flow antisymmetry)\n'
            '0               configuration summary\n'
            '          %3s = total number of networks read in\n'

            '         %4s = total number of mesh points\n'
            '         %4s = total number of panels\n'
            '' % (
                self.nsymmetry_planes,
                self.xz_symmetry, self.xy_symmetry,
                self.nnetworks,
                total_points, total_panels,
            )
        )

        return msg + msg2 + msg3



def split_into_sections(lines):
    sections = []
    section_names = []
    section = None
    section_name = None
    assert len(lines) > 0, 'lines=%s' % lines
    line = None
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

def remove_comments(lines, log):
    """
    Comment lines in Panair begin with an equal (=) sign.
    The file is easier to parse if we remove them.
    """
    lines2 = []
    for line in lines:
        line = line.rstrip().lower()
        if '=' in line:
            #print "line -> %r" % (line)
            if line[0] != '=':
                log.debug("line[0] -> %s" % line[0])
                line = line.split('=')[0]
                #log.debug("******")
                lines2.append(line)
            # else: skip
        else:
            lines2.append(line)
    return lines2
