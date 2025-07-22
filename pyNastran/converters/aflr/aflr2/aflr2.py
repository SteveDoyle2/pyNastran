"""
defines:
 - read_bedge(bedge_filename, beta_reverse=179.7, log=None, debug=False)
 - AFLR2(log=None, debug=False)
   - read_bedge(self, bedge_filename, beta_reverse=179.7)
   - write_nastran(self, bdf_filename)
   - write_fixed_points(self, fixed_points_filename)
   - merge_bedge(self, bedge, bedge_filename)
 - export_to_bedge(bedge_filename, nodes, grid_bcs, curves, subcurves, axis=1, log=None)

m3.bedge
    4
    4           2           4           4
   60          16          60          16         134           6
   76          22         252           2          66           2
  144           2
    0           0           0           0           1           1
    1           1           1           1           1           1
    1           1
   -6.633581169548380      -4.083969916911601
   ...

ncurves = 4
nsubcurves_per_curve = [4, 2, 4, 4]; n=4
nsubcurves = 14
nnodes_pack = [60, 16, 60, 16, 134, 6, 76, 22, 252, 2, 66, 2, 144, 2]; n=14
grid_bc = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; n=14

"""
import os
import sys
from copy import deepcopy
from typing import Optional

import numpy as np
from numpy.linalg import norm  # type: ignore
from cpylog import SimpleLogger, __version__ as CPYLOG_VERSION
if CPYLOG_VERSION > '1.6.0':
    from cpylog import get_logger
else:  # pragma: no cover
    from cpylog import get_logger2 as get_logger

from pyNastran.utils import print_bad_path, PathLike
from pyNastran.bdf.field_writer_8 import print_card_8


def read_bedge(bedge_filename: PathLike,
               beta_reverse: float=179.7,
               log: Optional[SimpleLogger]=None,
               debug: str | bool=False):
    """reads a *.bedge file"""
    model = AFLR2(log=log, debug=debug)
    model.read_bedge(bedge_filename, beta_reverse=beta_reverse)
    return model


class AFLR2:
    """defines methods for reading interfacing with AFLR2"""
    def __init__(self, log=None, debug: str | bool | None=False):
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
        self.log = get_logger(log, debug)
        self.debug = debug

        self.nodes = np.zeros((0, 3), dtype='float64')
        self.bars = np.zeros((0, 2), dtype='int32')
        self.curves = np.array([], dtype='int32')
        self.subcurves = np.array([], dtype='int32')
        self.grid_bc = np.array([], dtype='int32')
        self.grid_bcs = np.array([], dtype='int32')
        self.turn_angle = 0.0

    def read_bedge(self, bedge_filename: PathLike,
                   beta_reverse: float=179.7) -> None:
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
        nsubcurves_per_curve = -1 * np.zeros(ncurves, dtype='int32')
        #icurve = 0
        for icurve in range(ncurves):
            nsubcurvesi = int(data[i])
            nsubcurves_per_curve[icurve] = nsubcurvesi
            i += 1
        #del icurve

        self.log.debug('nsubcurves_per_curve = %s' % nsubcurves_per_curve)
        nsubcurves = nsubcurves_per_curve.sum()
        self.log.debug('nsubcurves = %s' % nsubcurves)

        self.log.debug('data[%d] = %s; nnodes[0]\n' % (i, data[i]))

        nnodes_pack = -1 * np.zeros(nsubcurves, dtype='int32')
        isubcurve = 0
        for isubcurve in range(nsubcurves):
            nnodesi = int(data[i])
            nnodes_pack[isubcurve] = nnodesi
            i += 1
        del isubcurve

        nnodes = nnodes_pack.sum()
        self.log.debug('nnodes_pack = %s' % nnodes_pack)
        self.log.debug('nnodes = %s' % nnodes)

        inode = 0
        isubcurve = 0
        isubcurvei = 0
        icurve = 0
        isubcurve_to_curve_map = -1 * np.zeros(nsubcurves, dtype='int32')
        for icurve in range(ncurves):
            nsubcurvesi = nsubcurves_per_curve[icurve]
            inode_curve_min[icurve] = inode
            #delta_node_id = 0
            for isubcurvei in range(nsubcurvesi):
                nnodesi: int = nnodes_pack[isubcurve]
                #max_node_id += nnodesi
                inode += nnodesi
                isubcurve_to_curve_map[isubcurve] = icurve
                isubcurve += 1
            inode_curve_max[icurve] = inode  # max_node_id
        inode_curve_min.append(nnodes)
        inode_curve_max.append(nnodes)
        self.log.debug("isubcurve_to_curve_map = %s" % isubcurve_to_curve_map)
        del icurve, isubcurvei

        # grid BC
        self.log.debug('data[%d]=%s; grid_bc[0]\n' % (i, data[i]))
        grid_bc = -1 * np.ones(nsubcurves, dtype='int32')
        for isubcurve in range(nsubcurves):
            grid_bci = int(data[i])
            grid_bc[isubcurve] = grid_bci
            i += 1
        del isubcurve
        self.log.debug('grid_bc = %s' % grid_bc)
        #if grid_bc[0] in [2]:
            #raise RuntimeError('fname=%s' % bedge_filename)
        self.log.debug('data[%d] = %s; nid.x\n' % (i, data[i]))
        #=============================================================

        nodes = np.zeros((nnodes, 3), dtype='float64')

        # just for testing
        #assert data[i] == '4.87406', 'val=%s fname=%s' % (data[i], bedge_filename)
        istart = i
        for inode in range(nnodes):
            try:
                node = [data[i], data[i+1], 0.]
            except IndexError:
                nactual = len(data)
                nexpected = istart + nnodes * 2
                msg = ('error with data; inode=%d; len(data)=%d '
                       'len(data_expected)=%d; nmissing=%d' % (
                           inode, nactual, nexpected, nexpected - nactual))
                raise IndexError(msg)
            #print('node[%d] = %s' % (inode, node))
            nodes[inode, :] = node
            #assert nodes[inode, 0] == 21., node
            i += 2
        self.log.debug('node[0] = %s' % nodes[0, :])
        self.log.debug('node[%d] = %s' % (inode, node))
        #self.log.debug('nodes = %s' % nodes[:, :2])
        del inode

        initial_normal_spacing = np.zeros(nnodes, dtype='float64')
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
        curves = np.zeros(nelements, dtype='int32')
        subcurves = np.zeros(nelements, dtype='int32')
        grid_bcs = np.zeros(nelements, dtype='int32')
        # self.log.debug('***ncurves = %s' % ncurves)

        for isubcurve in range(nsubcurves):
            nnodesi: int = nnodes_pack[isubcurve]
            icurve = isubcurve_to_curve_map[isubcurve]
            grid_bci = grid_bc[isubcurve]
            # self.log.debug('isubcurve=%s icurve=%s grid_bc=%s' % (isubcurve, icurve, grid_bci))

            curves[ielement0:ielement0+nnodesi] = icurve
            subcurves[ielement0:ielement0+nnodesi] = isubcurve
            grid_bcs[ielement0:ielement0+nnodesi] = grid_bci
            ielement0 += nnodesi
        del icurve

        inode0 = 0
        bars_list = []
        turn_angle_list = []
        for isubcurve in range(nsubcurves):
            nnodesi: int = nnodes_pack[isubcurve]
            nbars = nnodesi
            inode = inode0 + nnodesi
            inodes = np.array([inode0 + inodei for inodei in range(0, nnodesi)])

            bars = np.zeros((nbars, 2), dtype='int32')
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

            #print('inode_last[%d] = %s' % (isubcurve, inode_last))
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

            ibar = np.array(ibar, dtype='int32')
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

                #c = np.cross(v1, v2)
                #cn = norm(c, axis=1)

                #L1L2 = (L1 * L2)
                #izero = where(L1L2 == 0.0)[0]
                #L1L2[izero] = 1e-10
                #sin_theta = cn / L1L2
                #theta = degrees(arcsin(sin_theta))

            # convention from http://en.wikipedia.org/wiki/Triangle
            c = norm(nodes[n2, :2] - nodes[n1, :2], axis=1)  # c
            a = norm(nodes[n3, :2] - nodes[n2, :2], axis=1)  # a
            b = norm(nodes[n3, :2] - nodes[n1, :2], axis=1)  # b
            assert len(a) == len(n1), 'wrong size...check axis'

            cos_inner = (a**2 + c**2 - b**2)/(2 * a * c)
            beta = np.arccos(np.clip(cos_inner, -1., 1.))
            inan = np.where(np.isnan(beta))[0]
            beta[inan] = 0.0
            i180 = np.where(abs(beta) > np.radians(beta_reverse))[0]
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
            # self.log.debug('min_xy = %s' % min_xy)
            assert len(min_xy) == 2, min_xy

            # y/x for nodes 1, 2, and 3; find theta
            theta1g = np.arctan2(nodes[n1, 1] + dy, nodes[n1, 0] + dx)
            theta2g = np.arctan2(nodes[n2, 1] + dy, nodes[n2, 0] + dx)
            theta3g = np.arctan2(nodes[n3, 1] + dy, nodes[n3, 0] + dx)
            theta12g = theta2g - theta1g
            theta23g = theta3g - theta2g

            mag_theta = np.sign(theta12g - theta23g)
            izero = np.where(beta == 0.)[0]
            mag_theta[izero] = 0.
            #inotzero = where(beta != 0)
            #print('i')
            #izero = where(allclose(mag_theta, 0.))[0]
            #mag_theta[izero] = 1.0
            #print('mag_theta = %s' % mag_theta)
            turn_angle = mag_theta * beta
            #turn_angle = degrees(beta)
            turn_angle_list.append(turn_angle)

            #print('inodes[%d]=%s'  % (icurve, inodes))
            #print('bars[%d]=\n%s\n'  % (isubcurve, bars))
            for bari in bars:
                n1, n2 = bari
                #print(bari, nodes[n1, :], nodes[n2, :])
            #print('nodes[%d]=\n%s\n'  % (isubcurve, nodes[n1a, :2]))
            inode0 += nnodesi
            #break
        del isubcurve
        bars = np.vstack(bars_list)
        turn_angle = 0.
        if len(turn_angle_list) == 1:
            pass  # turn_angle = turn_angle_list[0]
        else:
            turn_angle = np.hstack(turn_angle_list)
        #print('nodes = \n%s' % nodes)

        self.nodes = nodes
        self.bars = bars
        self.curves = curves
        self.subcurves = subcurves
        self.grid_bc = grid_bc
        self.grid_bcs = grid_bcs
        self.log.debug('grid_bcs = %s' % np.unique(grid_bcs))
        self.turn_angle = turn_angle
        sys.stdout.flush()

    def write_tri(self, tri_filename: PathLike,
                  holes: list[tuple[float, float]] | None=None,
                  circles: list[tuple[float, float, float, int]] | None=None,
                  regions: list[tuple[float, float]] | None=None,
                  curves_to_skip: list[int] | None=None,
                  bc_tags: list[int] | None=None,
                  temp_tags_map: list[int] | None=None,
                  min_angle: float | str=20.0,
                  max_area: float | str=0.05,
                  tri_order: int=1,
                  plot_clear_regions: bool=True,
                  show: bool=False) -> None:  # pragma: no cover
        if bc_tags is None:
            bc_tags = []
        if temp_tags_map is None:
            temp_tags_map = {}

        triangle_input, options = get_triangle_input(
            self, curves_to_skip=curves_to_skip,
            holes=holes, circles=circles, regions=regions,
            min_angle=min_angle, max_area=max_area,
            tri_order=tri_order,
        )
        # -----------------------------------------
        self.log.info(f'options = {options}')
        import triangle as tr
        triangle_output = tr.triangulate(triangle_input, options)

        triangle_output_to_temp_deck(
            tri_filename, triangle_output,
            bc_tags=bc_tags, temp_tags_map=temp_tags_map)

        if plot_clear_regions and 'triangle_attributes' in triangle_output:
            # makes less messy plots
            del triangle_output['triangle_attributes']

        if show and tri_order == 1:
            import matplotlib.pyplot as plt
            tr.compare(plt, triangle_input, triangle_output)
            plt.show()

    def write_esp(self, esp_filename: PathLike,
                  curves_to_skip: Optional[list[int]]=None) -> None:
        """
        LINSEG  straight line segment
        CIRARC  circular arc
        ARC     alternative way of specifying a circular arc
        BEZIER  Bezier curve
        SPLINE  cubic spline

        skbeg       Wing_RX_TE 0          0
            linseg  Wing_TX_TE Wing_Semib 0
            linseg  Wing_TX    Wing_Semib 0
            linseg  Wing_RX    0          0
            linseg  Wing_RX_TE 0          0
        skend

        SKBEG 1.0 2.0 Z
        LINSEG 1.0+L 2.0 Z
        CIRARC 1.0+L-(1-s2)*H 2.0+s2*H Z \
               1.0+L-H        2.0+H    Z
        LINSEG 1.0 2.0+H Z
        LINSEG 1.0 2.0 Z
        SKEND
        """
        if curves_to_skip is None:
            curves_to_skip = []
        ucurves, nodes_pack, nsubcurves_list = get_nnodes_pack(
            self.curves, self.subcurves, self.log)
        # ncurves = len(ucurves)
        # self.log.debug(f'ucurves = {ucurves}')
        # self.log.debug(f'nodes_pack = {nodes_pack}')
        # self.log.debug(f'nsubcurves_list = {nsubcurves_list}')

        # ucurves = [0 1 2 3]
        # nodes_pack = [60, 16, 60, 16, 134, 6, 76, 22, 252, 2, 66, 2, 144, 2]
        # nsubcurves_list = [4, 2, 4, 4]
        # self.log.info(f'curves = {self.curves}')
        # self.log.info(f'subcurves = {self.subcurves}')
        # self.log.info(f'grid_bc = {self.grid_bc}')
        # self.log.info(f'nodes = {self.nodes.shape}')

        all_lines = []
        isubcurve = 0
        inode = 0
        assert len(self.grid_bc) == len(nodes_pack)
        ncurves = len(ucurves)
        nsubcurves_all = sum(nsubcurves_list)
        x0 = y0 = z0 = np.nan
        for icurve, nsubcurves in zip(ucurves, nsubcurves_list):
            word0 = 'skbeg'
            grid_bci = self.grid_bc[isubcurve]

            lines = []
            lines.append(f'# curve {icurve}/{ncurves}\n')
            lines.append(f'# nsubcurves={nsubcurves}\n')
            if icurve in curves_to_skip:
                for isubcurvei in range(nsubcurves):
                    nnodes = nodes_pack[isubcurve]
                    inode += nnodes
                    isubcurve += 1
                continue

            subcurve_word0 = ''
            for isubcurvei in range(nsubcurves):
                lines.append(f'# subcurve={isubcurve}/{nsubcurves_all}; curve={icurve}:{isubcurvei}\n')
                nnodes = nodes_pack[isubcurve]
                lines.append(f'#   nnodes={nnodes}\n')

                if isubcurvei == 0:
                    x0, y0, z0 = self.nodes[inode]

                subcurve_lines = []
                for inodei in range(nnodes):
                    x, y, z = self.nodes[inode]
                    if len(subcurve_lines) == 0:
                        subcurve_lines.append(f'{word0}  {subcurve_word0} {x} {y} 0.\n')
                        word0 = '     '
                    elif len(lines) == 1:
                        subcurve_lines.append(f'       SPLINE {x} {y} 0.3')
                    else:
                        subcurve_lines.append(f'       SPLINE {x} {y} 0.\n')

                    inode += 1
                lines += subcurve_lines
                isubcurve += 1
                subcurve_word0 = 'LINSEG'

            lines.append(f'       SPLINE {x0} {y0} {z0}  # endpt\n')
            lines.append('skend\n')
            all_lines.extend(lines)

        with open(esp_filename, 'w') as esp_file:
            esp_file.writelines(all_lines)

    def write_nastran(self, bdf_filename: PathLike) -> None:
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

    def write_fixed_points(self, fixed_points_filename: PathLike) -> None:
        """writes a *.csv file that can be read by pyNastranGUI to show the points"""
        with open(fixed_points_filename, 'w') as point_file:
            for i, node in enumerate(self.nodes):
                x, y, unused_z = node
                #print('bars = \n%s' % self.bars)
                ix, iy = np.where(i == self.bars)
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

    def merge_bedge(self, bedge, bedge_filename: PathLike) -> None:
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

        nodes = np.vstack([self.nodes, bedge.nodes])
        grid_bc = np.hstack([self.grid_bc, bedge.grid_bc])
        curves = np.hstack([self.curves, bedge.curves])
        subcurves = np.hstack([self.subcurves, bedge.subcurves])
        self.log.debug('ugrid_bcs = %s' % np.unique(grid_bc))

        export_to_bedge(bedge_filename,
                        nodes, grid_bc, curves, subcurves,
                        axis=2, log=self.log)


def _func_circle(R: float, x: float, y: float,
                 N: int) -> tuple[np.ndarray, np.ndarray]:
    i = np.arange(N)
    theta = i * 2 * np.pi / N
    pts = np.stack([np.cos(theta), np.sin(theta)], axis=1) * R
    seg = np.stack([i, i + 1], axis=1)
    seg[-1, 1] = 0
    xy = np.array([x, y])
    pts2 = pts + xy[np.newaxis, :]
    return pts2, seg


def get_triangle_input(model: AFLR2,
                       holes: list[tuple[float, float]] | None=None,
                       circles: list[tuple[float, float, float, int]] | None=None,
                       regions: list[tuple[float, float]] | None=None,
                       curves_to_skip: Optional[list[int]]=None,
                       min_angle: float | str=20.0,
                       max_area: float | str=0.05,
                       tri_order: int=1,
                       ) -> tuple[dict[str, np.ndarray], str]:
    """
    'qpa0.05'

    Triangle needs to be a closed curve and 2d (unlike bedge)
     - triangle has an MIT license

    https://github.com/inducer/meshpy
    https://rufat.be/triangle/examples.html

    Parameters
    ----------
    extra : str
        https://www.cs.cmu.edu/~quake/triangle.switch.html
        -A Assigns a regional attribute to each triangle that
           identifies what segment-bounded region it belongs to.
        -j Jettisons vertices that are not part of the final
           triangulation from the output .node file (including
           duplicate input vertices and vertices ``eaten'' by holes).
        -o2 Generates second-order subparametric elements with six nodes each.

    options : str; default='qpa0.05'
        *p - Triangulates a Planar Straight Line Graph.
        r - Refines a previously generated mesh.
        *q - Quality mesh generation with no angles smaller
            than 20 degrees. An alternate minimum angle may
            be specified after the q.
        *a - Imposes a maximum triangle area constraint.
            A fixed area constraint (that applies to every
            triangle) may be specified after the a, or varying
            areas may be read from the input dictionary.
        c - Encloses the convex hull with segments.
        D - Conforming Delaunay: use this switch if you want all
            triangles in the mesh to be Delaunay, and not just
            constrained Delaunay; or if you want to ensure that
            all Voronoi vertices lie within the triangulation.
        X - Suppresses exact arithmetic.
        S - Specifies the maximum number of added Steiner points.
        i - Uses the incremental algorithm for Delaunay
            triangulation, rather than the divide-and-conquer
            algorithm.
        F - Uses Steven Fortune’s sweepline algorithm for
            Delaunay triangulation, rather than the
            divide-and-conquer algorithm.
        l - Uses only vertical cuts in the divide-and-conquer
            algorithm. By default, Triangle uses alternating
            vertical and horizontal cuts, which usually improve
            the speed except with vertex sets that are small or
            short and wide. This switch is primarily of
            theoretical interest.
        s - Specifies that segments should be forced into the
            triangulation by recursively splitting them at their
            midpoints, rather than by generating a constrained
            Delaunay triangulation. Segment splitting is true to
            Ruppert’s original algorithm, but can create needlessly
            small triangles. This switch is primarily of
            theoretical interest.
        C - Check the consistency of the final mesh. Uses exact
            arithmetic for checking, even if the -X switch is
            used. Useful if you suspect Triangle is buggy.
        n - Return neighbor list in dict key ‘neighbors’
        e - Return edge list in dict key ‘edges’

    """
    assert tri_order in {1, 2}, tri_order
    if circles is None:
        circles = []
    if curves_to_skip is None:
        curves_to_skip = []
    if regions is None:
        regions = []
    if holes is None:
        holes = []

    if not isinstance(min_angle, str):
        min_angle = f'{min_angle:g}'
    if not isinstance(max_area, str):
        max_area = f'{max_area:g}'
    #-----------------------------------------

    ucurves, nodes_pack, nsubcurves_list = get_nnodes_pack(
        model.curves, model.subcurves, model.log)
    isubcurve = 0
    inode = 0
    segments = []
    for icurve, nsubcurves in zip(ucurves, nsubcurves_list):
        if icurve in curves_to_skip:
            for isubcurvei in range(nsubcurves):
                nnodes = nodes_pack[isubcurve]
                inode += nnodes
                isubcurve += 1
            continue

        inode0 = inode
        for isubcurvei in range(nsubcurves):
            nnodes = nodes_pack[isubcurve]
            segmenti = np.zeros((nnodes, 2), dtype='int32')
            i = inode + np.arange(0, nnodes)
            segmenti[:, 0] = i
            segmenti[:, 1] = i + 1
            inode += nnodes
            isubcurve += 1
            segments.append(segmenti)
        segmenti[-1, 1] = inode0
        inode0 = inode

    all_nodes_list = [model.nodes[:, :2]]
    for circle in circles:
        (R, x, y, N) = circle
        pts0, seg0 = _func_circle(R, x, y, N)
        all_nodes_list.append(pts0)
        segments.append(seg0 + inode)
        inode += len(pts0)

    segment = np.vstack(segments)

    all_nodes = np.vstack(all_nodes_list)
    if len(curves_to_skip):
        # filter unused points
        ipoints = np.unique(segment.ravel())
        all_nodes = all_nodes[ipoints, :]
        segment = np.searchsorted(ipoints, segment)

    #all_nodes = self.nodes[ipoints, :2]
    assert all_nodes.shape[1] == 2, all_nodes.shape

    triangle_input = {
        'vertices': all_nodes,
        'segments': segment,
    }
    if len(holes):
        triangle_input['holes'] = holes

    option_flags = [
        'p',  # generate planar line graph (respect the hard edges)
        'j',  # jettison/remove unused points
    ]
    if len(regions):
        triangle_input['regions'] = regions
        option_flags.append('A')

    if min_angle != '':
        option_flags.append(f'q{min_angle}')
    if max_area != '':
        option_flags.append(f'a{max_area}')

    if tri_order == 2:
        option_flags.append('-o2')
    options = ''.join(option_flags)
    # options = 'qpAa0.05'
    # options = 'qpa0.05'
    return triangle_input, options


def triangle_output_to_temp_deck(tri_filename: PathLike,
                                 triangle_output: dict[str, np.ndarray],
                                 bc_tags: list[int] | None=None,
                                 temp_tags_map: list[int] | None=None) -> None:
    #print(triangle_output.keys())

    # Input
    # -----
    # vertices
    # segments
    # holes
    # regions

    # Output
    # ------
    # vertex_markers: useful for BCs?
    # triangles
    # triangle_attributes: generated by A and regiosn
    # segment_markers: useful for BCs; remeshed by triangle
    #print('triangle_output', triangle_output)
    vertices = triangle_output['vertices']
    vertex_markers = triangle_output['vertex_markers'].flatten().tolist()
    segment_markers = triangle_output['segment_markers'].flatten()
    tris = triangle_output['triangles'] + 1
    nvert = len(vertices)
    ntri = len(tris)

    nids = np.arange(1, nvert+1, dtype='int32')
    eids = np.arange(1, ntri+1, dtype='int32')
    if 'triangle_attributes' in triangle_output:
        pids = triangle_output['triangle_attributes'].ravel().astype('int32')
    else:
        pids = np.ones(ntri, dtype='int32')

    # print(f'vertex_markers = {vertex_markers}')
    # print(f'segment_markers = {segment_markers}')
    upids = np.unique(pids)

    segments = triangle_output['segments']
    segment_markers = triangle_output['segment_markers']
    useg_markers = np.unique(segment_markers)
    bc_nodes_list = []
    temp_nodes_list = []
    temp_value_list = []
    for iseg_marker in useg_markers:
        if iseg_marker not in bc_tags and iseg_marker not in temp_tags_map:
            continue

        if iseg_marker in bc_tags:
            assert iseg_marker not in temp_tags_map, f'iseg_marker={iseg_marker} and is in bc_tags and temp_tags_map'
            iseg = np.where(segment_markers == iseg_marker)[0]
            nodes = np.unique(segments[iseg, :].ravel())
            bc_nodes_list.append(nodes + 1)
        elif iseg_marker in temp_tags_map:
            temp = temp_tags_map[iseg_marker]
            iseg = np.where(segment_markers == iseg_marker)[0]
            nodes = np.unique(segments[iseg, :].ravel())
            nnodesi = len(nodes)
            temps = np.full(nnodesi, temp, dtype='float64')
            temp_nodes_list.append(nodes + 1)
            temp_value_list.append(temps)

    spc_id = 101
    load_id = 102
    if len(bc_nodes_list):
        bc_nodes = np.hstack(bc_nodes_list).tolist()
    else:
        bc_nodes = np.array([], dtype='int32').tolist()

    if len(temp_nodes_list):
        temp_nodes = np.hstack(temp_nodes_list).tolist()
        temp_values = np.hstack(temp_value_list).tolist()
    else:
        temp_nodes = np.array([], dtype='int32').tolist()
        temp_values = np.array([], dtype='int32').tolist()

    tri_dim = tris.shape[1]
    with open(tri_filename, 'w') as tri_file:
        tri_file.write('SOL 101\n')
        tri_file.write('CEND\n')
        tri_file.write(f'DISP(PLOT) = ALL\n')
        tri_file.write(f'STRESS(PLOT) = ALL\n')
        tri_file.write(f'STRAIN(PLOT) = ALL\n')
        if len(bc_nodes):
            tri_file.write(f'SUBCASE 1\n')
            tri_file.write(f'  SPC = {spc_id}\n')
        if len(temp_nodes):
            tri_file.write(f'  LOAD = {load_id}\n')
        tri_file.write('BEGIN BULK\n')
        for nid, (x, y) in zip(nids, vertices):
            card = ['GRID', nid, '', x, y, 0.]
            tri_file.write(print_card_8(card))
        if len(bc_nodes):
            card = ['SPC1', spc_id, '123456', ] + bc_nodes
            tri_file.write(print_card_8(card))

        for temp_node, temp_value in zip(temp_nodes, temp_values):
            card = ['TEMP', load_id, temp_node, temp_value]
            tri_file.write(print_card_8(card))

        if tri_dim == 3:
            for eid, pid, (n1, n2, n3) in zip(eids, pids, tris):
                card = ['CTRIA3', eid, pid, n1, n2, n3]
                tri_file.write(print_card_8(card))
        else:
            assert tri_dim == 6, tri_dim
            for eid, pid, (n1, n2, n3, n4, n5, n6) in zip(eids, pids, tris):
                # the output of triangle is weird
                card = ['CTRIA6', eid, pid, n1, n2, n3, n6, n4, n5]
                tri_file.write(print_card_8(card))

        for pid in upids:
            mid = pid
            tri_file.write(print_card_8(['PSHELL', pid, mid, 0.1]))
            tri_file.write(print_card_8(['MAT1', mid, 3.0e7, None, 0.3]))
        tri_file.write('ENDDATA\n')


def _flip_value(lst: list[int]) -> list[int]:
    """flips a 0 to 1 and vice-versa"""
    return [
        0 if val == 1 else 1
        for val in lst
    ]


def get_nnodes_pack(curves: np.ndarray, subcurves: np.ndarray,
                    log: SimpleLogger) -> list:
    # write the curves/subcurves
    ucurves = np.unique(curves)
    # ncurves = len(ucurves)
    nodes_pack = []

    log.debug('looping over ucurves=%s' % ucurves)
    nsubcurves_list = []

    all_usubcurves = set()
    for ucurve in ucurves:
        i = np.where(curves == ucurve)[0]

        subcurvesi = subcurves[i]
        # log.debug('ucurve=%s i=%s subcurves[i]=%s' % (ucurve, i, subcurvesi))

        usubcurves = np.unique(subcurvesi)
        nsubcurves = len(usubcurves)
        nsubcurves_list.append(nsubcurves)
        #f.write('  %s\n' % nsubcurves)

        for usubcurve in usubcurves:
            # if usubcurve not in all_usubcurves:
                # log.debug(f'all_usubcurves = {all_usubcurves}')
                # stop???

            all_usubcurves.add(usubcurve)
            j = np.where(subcurvesi == usubcurve)[0]
            nnodesi = len(j)
            # log.debug('usubcurve=%s j=%s subcurves[j]=%s nnodes=%s' % (
            #     usubcurve, j, subcurvesi[j], nnodesi))
            #f.write('%s ' % nnodesi)
            nodes_pack.append(nnodesi)
        #f.write('\n')
    del nsubcurves

    return ucurves, nodes_pack, nsubcurves_list


def export_to_bedge(bedge_filename: PathLike,
                    nodes, grid_bcs, curves, subcurves,
                    axis: int=1, log: Optional[SimpleLogger]=None):
    """
    Creates a bedge file

    Parameters
    ----------
    bedge_filename : str | Path
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

    ucurves, nodes_pack, nsubcurves_list = get_nnodes_pack(curves, subcurves, log)
    ncurves = len(ucurves)
    with open(bedge_filename, 'w') as bedge_file:
        bedge_file.write('%s\n' % ncurves)

        for nsubcurvesi in nsubcurves_list:
            bedge_file.write('  %s' % nsubcurvesi)
        bedge_file.write('\n')
        nsubcurves = len(grid_bcs)
        log.debug('nsubcurves = %s?' % nsubcurves)
        if nsubcurves > 30:
            msg = 'Are you sure you are merging model.grid_bc and not model.grid_bcs?\n'
            msg += f'nsubcurves={nsubcurves:d}'
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
