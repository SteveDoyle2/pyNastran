from typing import Union
import numpy as np
from numpy import zeros, empty
from cpylog import get_logger2


def convert_to_int(sline):
    return [int(spot) for spot in sline]

def read_avus(avus_filename, log=None, debug=False):
    model = AvusGrid(log=log, debug=debug)
    model.read_avus_grid(avus_filename)
    return model

class AvusGrid:
    def __init__(self, log=None, debug: Union[str, bool, None]=False):
        self.log = get_logger2(log=log, debug=debug)
        self.infilename = None

        self.zones = None
        self.npoints = None
        self.nfaces = None
        self.ndm = None
        self.npatches = []

        self.nodes = empty(shape=0)
        self.tri_elements = empty(shape=0)
        self.quad_elements = empty(shape=0)
        self.tet_elements = empty(shape=0)
        self.hexa_elements = empty(shape=0)

    def _read_points(self, infile, unused_zones, npoints):
        """
        npoints
        x1 y1 z1
        x2 y2 z2
        ...
        """
        npoints_total = sum(npoints)
        points = zeros((npoints_total, 3), dtype='float32')
        for i in range(npoints_total):
            line = infile.readline()
            xyz = line.split()
            points[i, :] = xyz
        return points

    def _read_faces(self, infile, nfaces):
        """
        # reads faces - needed for viewing

        dim n1 n2 n3 ??? ???
        3 16160 7818 18736 70913 141972
        3 7819  7823 7820  78065 -9
        """
        self.log.debug('read_faces; nfaces=%s' % nfaces)
        ipoints = []
        assert len(nfaces) == 1, 'nfaces=%s must be length 1' % nfaces
        for nface in nfaces:
            #faces = np.zeros((nface, 5), dtype='int32')
            iface = 0
            for unused_iface in range(nface):
                line = infile.readline()
                sline_ints = convert_to_int(line.split())

                # I'm pretty sure this just means 3d
                assert len(sline_ints) == 6, sline_ints

                numppf = sline_ints[0] # must be 2 or 3 - dimensionality

                # read face
                ipoint = []
                for ipointi in range(1, numppf + 1):
                    ipoint.append(sline_ints[ipointi])

                # read cell
                # negative - outer boundary element
                # zero     - inner boundary element
                # positive - physical interface between cells
                try:
                    ilow = sline_ints[ipointi + 1]
                    ihigh = sline_ints[ipointi + 2]
                    icell = [ilow, ihigh]
                except IndexError:
                    msg = 'length should be %s; line=%s\n' % (numppf + 2, sline_ints)
                    raise IndexError(msg)

                if min(icell) <= 0:
                    #print(icell)
                    # if it's a boundary cell (inner or outer)
                    ipoints.append(ipoint)
                    #faces[iface, :] = sline_ints[1:]
                    iface += 1
        #print(ipoint)
        #print(sline_ints)
        faces2 = np.array(ipoints, dtype='int32')
        #print('faces2:\n%s' % faces2)
        return faces2

    def _read_header(self, infile):
        line = infile.readline().split()
        ndm, nzones, npatches = convert_to_int(line)
        self.log.debug('read_zones; ndm=%s nzones=%s npatches=%s' % (ndm, nzones, npatches))
        assert nzones >= 1, 'nzones=%s' % nzones
        zones = [z for z in range(nzones)]

        npoints = []
        nfaces = []
        ncells = []
        mxppfs = []
        mxfpcs = []
        for unused_izone in zones:
            sline = infile.readline().split()
            (npoint, nface, ncell, mxppf, mxfpc) = convert_to_int(sline)
            assert npoint > 0, npoint

            npoints.append(npoint)
            nfaces.append(nface)
            ncells.append(ncell)
            mxppfs.append(mxppf)
            mxfpcs.append(mxfpc)
            self.log.debug("  npoint=%s nface=%s ncells=%s mxppf=%s mxfpc=%s" % (
                npoint, nface, ncell, mxppf, mxfpc))

        self.ndm = ndm
        self.npatches = npatches
        self.npoints = npoints
        self.zones = zones
        self.nfaces = nfaces
        return zones, npoints, nfaces, ncells, mxppfs, mxfpcs

    def read_avus_grid(self, avus_filename):
        self.infilename = avus_filename
        self.log.debug('reading %s' % self.infilename)
        with open(self.infilename, 'r') as infile:

            # TODO: what are mxppfs, mxfpcs???
            (zones, npoints, nfaces, unused_ncells,
             unused_mxppfs, unused_mxfpcs) = self._read_header(infile)
            points = self._read_points(infile, zones, npoints)
            faces = self._read_faces(infile, nfaces)

        self.nodes = points
        self.tri_elements = np.asarray(faces, dtype='int32') - 1

    def write_avus(self, avus_filename):
        with open(avus_filename, 'w') as unused_avus_file:
            pass

    @property
    def nnodes(self):
        return self.nodes.shape[0]

    @property
    def nelements(self):
        return self.tri_elements.shape[0]
