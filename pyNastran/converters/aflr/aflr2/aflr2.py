"""
defines:
 - read_bedge(bedge_filename, beta_reverse=179.7, log=None, debug=False)
 - AFLR2(log=None, debug=False)
   - read_bedge(self, bedge_filename, beta_reverse=179.7)
   - write_nastran(self, bdf_filename)
   - write_fixed_points(self, fixed_points_filename)
   - merge_bedge(self, bedge, bedge_filename)
 - export_to_bedge(bedge_filename, nodes, grid_bcs, curves, subcurves, axis=1, log=None)

"""
import os
import sys
from copy import deepcopy

import numpy as np
from numpy import (zeros, array, vstack, hstack, where,
                   arctan2, arccos, sign, isnan, radians, unique)
from numpy.linalg import norm  # type: ignore
from cpylog import get_logger2

from pyNastran.utils import print_bad_path
from pyNastran.bdf.field_writer_8 import print_card_8

def read_bedge(bedge_filename, beta_reverse=179.7, log=None, debug=False):
    """reads a *.bedge file"""
    model = AFLR2(log=log, debug=debug)
    model.read_bedge(bedge_filename, beta_reverse=beta_reverse)
    return model

class AFLR2:
    """defines methods for reading interfacing with AFLR2"""
    def __init__(self, log=None, debug=False):
        """
        Initializes the AFLR2 object

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
        self.debug = debug

        self.nodes = None
        self.bars = None
        self.curves = None
        self.subcurves = None
        self.grid_bc = None
        self.grid_bcs = None
        self.turn_angle = None

    def read_bedge(self, bedge_filename, beta_reverse=179.7):
        """reads a *.bedge file"""
        ext = os.path.splitext(bedge_filename)[1].lower()
        assert ext == '.bedge', print_bad_path(bedge_filename)

        with open(bedge_filename, 'r') as bedge_file:
            data = bedge_file.read().split()

        i = 0
        ncurves = int(data[i])
        self.log.debug('ncurves = %s' % ncurves)
        i += 1

        inode_curve_min = [None] * ncurves
        inode_curve_max = [None] * ncurves

        nsubcurves_per_curve = [None] * ncurves
        icurve = 0
        for icurve in range(ncurves):
            nsubcurvesi = int(data[i])
            nsubcurves_per_curve[icurve] = nsubcurvesi
            i += 1
        del icurve

        self.log.debug('nsubcurves_per_curve = %s' % nsubcurves_per_curve)
        nsubcurves = sum(nsubcurves_per_curve)
        self.log.debug('nsubcurves = %s' % nsubcurves)

        self.log.debug('data[%i] = %s; nnodes[0]\n' % (i, data[i]))

        nnodes_pack = [None] * nsubcurves
        isubcurve = 0
        for isubcurve in range(nsubcurves):
            nnodesi = int(data[i])
            nnodes_pack[isubcurve] = nnodesi
            i += 1
        del isubcurve

        nnodes = sum(nnodes_pack)
        self.log.debug('nnodes_pack = %s' % nnodes_pack)
        self.log.debug('nnodes = %s' % nnodes)


        inode = 0
        isubcurve = 0
        isubcurvei = 0
        isubcurve_to_curve_map = [None] * nsubcurves
        icurve = 0
        for icurve in range(ncurves):
            nsubcurvesi = nsubcurves_per_curve[icurve]
            inode_curve_min[icurve] = inode
            #delta_node_id = 0
            for isubcurvei in range(nsubcurvesi):
                nnodesi = nnodes_pack[isubcurve]
                #max_node_id += nnodesi
                inode += nnodesi
                isubcurve_to_curve_map[isubcurve] = icurve
                isubcurve += 1
            inode_curve_max[icurve] = inode # max_node_id
        inode_curve_min.append(nnodes)
        inode_curve_max.append(nnodes)
        self.log.debug("isubcurve_to_curve_map = %s" % isubcurve_to_curve_map)
        del icurve, isubcurvei

        # grid BC
        self.log.debug('data[%i]=%s; grid_bc[0]\n' % (i, data[i]))
        grid_bc = [None] * nsubcurves
        for isubcurve in range(nsubcurves):
            grid_bci = int(data[i])
            grid_bc[isubcurve] = grid_bci
            i += 1
        del isubcurve
        self.log.debug('grid_bc = %s' % grid_bc)
        #if grid_bc[0] in [2]:
            #raise RuntimeError('fname=%s' % bedge_filename)
        self.log.debug('data[%i] = %s; nid.x\n' % (i, data[i]))
        #=============================================================

        nodes = zeros((nnodes, 3), dtype='float64')

        # just for testing
        #assert data[i] == '4.87406', 'val=%s fname=%s' % (data[i], bedge_filename)
        istart = i
        for inode in range(nnodes):
            try:
                node = [data[i], data[i+1], 0.]
            except IndexError:
                nactual = len(data)
                nexpected = istart + nnodes * 2
                msg = ('error with data; inode=%s; len(data)=%s '
                       'len(data_expected)=%s; nmissing=%s' % (
                           inode, nactual, nexpected, nexpected - nactual))
                raise IndexError(msg)
            #print('node[%i] = %s' % (inode, node))
            nodes[inode, :] = node
            #assert nodes[inode, 0] == 21., node
            i += 2
        self.log.debug('node[0] = %s' % nodes[0, :])
        self.log.debug('node[%i] = %s' % (inode, node))
        #self.log.debug('nodes = %s' % nodes[:, :2])
        del inode


        initial_normal_spacing = zeros(nnodes, dtype='float64')
        if len(data) != i:
            initial_normal_spacingi = data[i]
            if '*' in initial_normal_spacingi:
                self.log.debug('initial_normal_spacingi = %s' % initial_normal_spacingi)
                nnodesi, initial_normal_spacingi = initial_normal_spacingi.split('*')
                nnodesi = int(nnodesi)
                initial_normal_spacing[:nnodesi] = initial_normal_spacingi
                i += 1
            else:
                for inode in range(nnodes):
                    initial_normal_spacingi = data[i]
                    initial_normal_spacing[inode] = initial_normal_spacingi
                    i += 1

        assert len(data) == i, 'len(data)=%s i=%s' % (len(data), i)

        #==================================================
        # build the outputs

        ielement0 = 0
        nelements = nnodes
        curves = zeros(nelements, dtype='int32')
        subcurves = zeros(nelements, dtype='int32')
        grid_bcs = zeros(nelements, dtype='int32')
        self.log.debug('***ncurves = %s' % ncurves)

        for isubcurve in range(nsubcurves):
            nnodesi = nnodes_pack[isubcurve]
            icurve = isubcurve_to_curve_map[isubcurve]
            grid_bci = grid_bc[isubcurve]
            self.log.debug('isubcurve=%s icurve=%s grid_bc=%s' % (isubcurve, icurve, grid_bci))

            curves[ielement0:ielement0+nnodesi] = icurve
            subcurves[ielement0:ielement0+nnodesi] = isubcurve
            grid_bcs[ielement0:ielement0+nnodesi] = grid_bci
            ielement0 += nnodesi
        del icurve

        inode0 = 0
        bars_list = []
        turn_angle_list = []
        for isubcurve in range(nsubcurves):
            nnodesi = nnodes_pack[isubcurve]
            nbars = nnodesi
            inode = inode0 + nnodesi
            inodes = array([inode0 + inodei for inodei in range(0, nnodesi)])

            bars = zeros((nbars, 2), dtype='int32')
            inode_last = inodes[-1]

            icurve = isubcurve_to_curve_map[isubcurve]
            #print('isubcurve=%s icurve=%s inode_last=%s nnodes=%s' % (
                #isubcurve, icurve, inode_last, nnodes))
            #print('inode_curve_max = %s' % (inode_curve_max))
            #print('inode_curve_min = %s' % (inode_curve_min))
            if inode_last + 1 == nnodes:
                inode_last = inode_curve_min[icurve]
                #print('A')
                #inode_last = 0
            elif inode_last + 1 == inode_curve_max[icurve]:
                inode_last = inode_curve_min[icurve]  # inodes[0] ???
                #print('B')
            #elif inode_last == inode_curve_min[icurve + 1]:
            else:
                #print('C')
                inode_last += 1


            #print('inode_last[%i] = %s' % (isubcurve, inode_last))
            bars[:, 0] = inodes
            bars[:-1, 1] = inodes[1:]

            bars[-1, 1] = inode_last
            assert inode_last != nnodes
            bars_list.append(bars)

            n1 = bars[:, 0]
            n2 = bars[:, 1]
            #n1 = n1a
            #n2a = n2
            #print "n1a = ", n1a, len(n1a)
            #print "n2a = ", n2a, len(n2a)

            ibar = list(range(1, nbars))
            ibar.append(0)

            ibar = array(ibar, dtype='int32')
            n1 = n1
            n3 = bars[ibar, 1]
            #print "n1a = ", n1a, len(n1a)
            #print "n2a = ", n2a, len(n2a)
            #print "n1b = ", n1b, len(n1b)
            #print "n2b = ", n2b, len(n2b)

            #if 0:
                #v1 = nodes[n2, :] - nodes[n1, :]
                #v2 = nodes[n3, :] - nodes[n1, :]

                #L1 = norm(v1, axis=1)
                #L2 = norm(v2, axis=1)

                #c = cross(v1, v2)
                #cn = norm(c, axis=1)

                #L1L2 = (L1 * L2)
                #izero = where(L1L2 == 0.0)[0]
                #L1L2[izero] = 1e-10
                #sin_theta = cn / L1L2
                #theta = degrees(arcsin(sin_theta))

            # convention from http://en.wikipedia.org/wiki/Triangle
            c = norm(nodes[n2, :2] - nodes[n1, :2], axis=1) # c
            a = norm(nodes[n3, :2] - nodes[n2, :2], axis=1) # a
            b = norm(nodes[n3, :2] - nodes[n1, :2], axis=1) # b
            assert len(a) == len(n1), 'wrong size...check axis'

            cos_inner = (a**2 + c**2 - b**2)/(2 * a * c)
            beta = arccos(np.clip(cos_inner, -1., 1.))
            inan = where(isnan(beta))[0]
            beta[inan] = 0.0
            i180 = where(abs(beta) > radians(beta_reverse))[0]
            beta[i180] = 0.0
            #print('beta = %s' % beta)
            #print('beta_nonzero = %s' % beta[where(beta != 0)[0]])
            #print('beta_nonzero deg = %s' % degrees(beta[where(beta != 0)[0]]))

            #if 0:
                #for inani in inan:
                    #nids_failed = [n1[inani], n2[inani], n3[inani]]
                    #print('nodes = %s' % nids_failed)
                    #for nid in nids_failed:
                        #print('  nodes[%3i] = %s' % (nid, nodes[nid, :]))

            #centroid = (nodes[n1, :] + nodes[n2, :] + nodes[n3, :]) / 3.
            #xcentroid = centroid[:, 0]
            #ycentroid = centroid[:, 1]

            # maybe shift the node to the centroid of the element to fix sign?
            min_xy = nodes[:, :2].min(axis=0)
            delta_xy = 2.0 * abs(min_xy)
            dx, dy = delta_xy
            self.log.debug('min_xy = %s' % min_xy)
            assert len(min_xy) == 2, min_xy

            # y/x for nodes 1, 2, and 3; find theta
            theta1g = arctan2(nodes[n1, 1] + dy, nodes[n1, 0] + dx)
            theta2g = arctan2(nodes[n2, 1] + dy, nodes[n2, 0] + dx)
            theta3g = arctan2(nodes[n3, 1] + dy, nodes[n3, 0] + dx)
            theta12g = theta2g - theta1g
            theta23g = theta3g - theta2g

            mag_theta = sign(theta12g - theta23g)
            izero = where(beta == 0.)[0]
            mag_theta[izero] = 0.
            #inotzero = where(beta != 0)
            #print('i')
            #izero = where(allclose(mag_theta, 0.))[0]
            #mag_theta[izero] = 1.0
            #print('mag_theta = %s' % mag_theta)
            turn_angle = mag_theta * beta
            #turn_angle = degrees(beta)
            turn_angle_list.append(turn_angle)


            #print('inodes[%i]=%s'  % (icurve, inodes))
            #print('bars[%i]=\n%s\n'  % (isubcurve, bars))
            for bari in bars:
                n1, n2 = bari
                #print(bari, nodes[n1, :], nodes[n2, :])
            #print('nodes[%i]=\n%s\n'  % (isubcurve, nodes[n1a, :2]))
            inode0 += nnodesi
            #break
        del isubcurve
        bars = vstack(bars_list)
        if len(turn_angle_list) == 1:
            pass # turn_angle = turn_angle_list[0]
        else:
            turn_angle = hstack(turn_angle_list)
        #print('nodes = \n%s' % nodes)

        self.nodes = nodes
        self.bars = bars
        self.curves = curves
        self.subcurves = subcurves
        self.grid_bc = grid_bc
        self.grid_bcs = grid_bcs
        self.log.debug('grid_bcs = %s' % unique(grid_bcs))
        self.turn_angle = turn_angle
        sys.stdout.flush()

    def write_nastran(self, bdf_filename):
        """converts the *.bedge to a nastran *.bdf"""
        with open(bdf_filename, 'w') as bdf_file:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')
            nid = 1
            for x, y, z in self.nodes:
                card = ['GRID', nid, None, x, y, z]
                bdf_file.write(print_card_8(card))
                nid += 1
            bdf_file.write('ENDDATA\n')

    def write_fixed_points(self, fixed_points_filename):
        """writes a *.csv file that can be read by pyNastranGUI to show the points"""
        with open(fixed_points_filename, 'w') as point_file:
            for i, node in enumerate(self.nodes):
                x, y, unused_z = node
                #print('bars = \n%s' % self.bars)
                ix, iy = where(i == self.bars)
                #print('ix[%s]=%s iy[%s]=%s' % (i, ix, i, iy))
                inodes_close = self.bars[ix, _flip_value(iy)]
                nodes_close = self.nodes[inodes_close, :]
                #print('p=(%s,%s); close=\n%s' % (x, y, str(nodes_close)))
                #print('p=(%s,%s); iclose=\n%s' % (x, y, str(inodes_close)))

                # (N,3) - (3,) subtraction
                delta = nodes_close - node
                normi = norm(delta, axis=1)
                max_edge_length = normi.max()
                #print(self.bars[ix])
                #print(nodes_close)
                #print('delta\n%s' % delta)
                #print(normi)
                #max_edge_length = max(0.03, max_edge_length)
                point_file.write('%s %s %s\n' % (x, y, max_edge_length))

    def merge_bedge(self, bedge, bedge_filename):
        """merges two bedge models into a single new *.bedge file"""
        bedge = deepcopy(bedge)

        curve1_max = self.curves.max()
        subcurve1_max = self.subcurves.max()

        curve2_max = bedge.curves.max()
        subcurve2_max = bedge.subcurves.max()

        bedge.curves += curve1_max + 1
        bedge.subcurves += subcurve1_max + 1
        self.log.debug('curve1_max=%s curve2_max=%s' % (curve1_max, curve2_max))
        self.log.debug('subcurve1_max=%s subcurve2_max=%s' % (subcurve1_max, subcurve2_max))

        nodes = vstack([self.nodes, bedge.nodes])
        grid_bc = hstack([self.grid_bc, bedge.grid_bc])
        curves = hstack([self.curves, bedge.curves])
        subcurves = hstack([self.subcurves, bedge.subcurves])
        self.log.debug('ugrid_bcs = %s' % unique(grid_bc))

        export_to_bedge(bedge_filename,
                        nodes, grid_bc, curves, subcurves, axis=2, log=self.log)

def _flip_value(lst):
    """flips a 0 to 1 and vice-versa"""
    return [
        0 if val == 1 else 1
        for val in lst
    ]

def export_to_bedge(bedge_filename,
                    nodes, grid_bcs, curves, subcurves, axis=1, log=None):
    """
    Creates a bedge file

    Parameters
    ----------
    bedge_filename : str
        the *.bedge file
    nodes : ???
        ???
    grid_bcs : ???
        ???
        source is model.grid_bc, not model.grid_bcs
    curves : ???
        ???
    subcurves : ???
        ???
    axis : int; default=1
        the axis to remove (nodes in Nx3)
    log : Logger(); default=None
        a required logging object

    """
    log.debug('bedge_filename = %s' % bedge_filename)
    log.debug('grid_bc = %s' % grid_bcs)
    #if bedge_filename == 'farfield.bedge':
        #print(grid_bcs)

    with open(bedge_filename, 'w') as bedge_file:
        # ncurves
        ucurves = unique(curves)
        ncurves = len(ucurves)
        bedge_file.write('%s\n' % ncurves)

        # write the curves/subcurves
        nodes_pack = []

        log.debug('looping over ucurves=%s' % ucurves)
        nsubcurves_list = []

        all_usubcurves = set()
        for ucurve in ucurves:
            i = where(curves == ucurve)[0]

            subcurvesi = subcurves[i]
            log.debug('ucurve=%s i=%s subcurves[i]=%s' % (ucurve, i, subcurvesi))

            usubcurves = unique(subcurvesi)
            nsubcurves = len(usubcurves)
            nsubcurves_list.append(nsubcurves)
            #f.write('  %s\n' % nsubcurves)

            for usubcurve in usubcurves:
                if usubcurve not in all_usubcurves:
                    log.debug(all_usubcurves)
                    # stop???

                all_usubcurves.add(usubcurve)
                j = where(subcurvesi == usubcurve)[0]
                nnodesi = len(j)
                log.debug('usubcurve=%s j=%s subcurves[j]=%s nnodes=%s' % (
                    usubcurve, j, subcurvesi[j], nnodesi))
                #f.write('%s ' % nnodesi)
                nodes_pack.append(nnodesi)
            #f.write('\n')
        del nsubcurves

        for nsubcurvesi in nsubcurves_list:
            bedge_file.write('  %s' % nsubcurvesi)
        bedge_file.write('\n')
        nsubcurves = len(grid_bcs)
        log.debug('nsubcurves = %s?' % nsubcurves)
        if nsubcurves > 30:
            msg = 'Are you sure you are merging model.grid_bc and not model.grid_bcs?\n'
            msg += 'nsubcurves=%s' % (nsubcurves)
            raise RuntimeError(msg)
        #assert nsubcurves == len(nsubcurves_list) #  wrong check...

        for nnodes in nodes_pack:
            bedge_file.write('  %s' % nnodes)
        bedge_file.write('\n')

        # grid bc
        # len(grid_bcs) == nsubcurves, not nnodes ironically enough
        for grid_bc in grid_bcs:
            bedge_file.write('%s ' % grid_bc)
        bedge_file.write('\n')


        # nodes
        #iaxis = where([0, 1, 2] != axis)[0]
        #print(iaxis)
        for node in nodes:
            x, y, z = node
            if axis == 0:
                bedge_file.write('%s %s\n' % (y, z))
            elif axis == 1:
                bedge_file.write('%s %s\n' % (x, z))
            elif axis == 2:
                bedge_file.write('%s %s\n' % (x, y))
            else:
                raise RuntimeError(axis)
        #initial_normal_spacing = '112*1.0e-3'
