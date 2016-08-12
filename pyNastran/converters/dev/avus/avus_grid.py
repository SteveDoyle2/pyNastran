from numpy import array, zeros, empty

def convert_to_int(line):
    sline = line.strip().split()
    return [int(spot) for spot in sline]

def convert_to_float(line):
    sline = line.strip().split()
    return [float(spot) for spot in sline]


class AvusGrid(object):
    def __init__(self, log=None, debug=False):
        self.log = log

        self.infilename = None
        self.infile = None

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

    def _read_points(self, zones, npoints):
        npoints_total = sum(npoints)
        points = zeros((npoints_total, 3), dtype='float32')
        for i in range(npoints_total):
            line = self.infile.readline()
            xyz = line.split()
            points[i, :] = xyz
        return points

    def _read_faces(self, nfaces):
        """
        # reads faces - needed for viewing
        """
        ipoints = []
        for nf in nfaces:
            for nfi in range(nf):
                line = self.infile.readline()
                sline = convert_to_int(line)
                numppf = sline[0] # must be 2 or 3 - dimensionality

                # read face
                ipoint = []
                for ipointi in range(1, numppf + 1):
                    ipoint.append(sline[ipointi])

                # read cell
                # negative - outer boundary element
                # zero     - inner boundary element
                # positve  - physical interface between cells
                ilow = sline[ipointi + 1]
                ihigh = sline[ipointi + 2]
                icell = [ilow, ihigh]

                if min(icell) <= 0:
                    # if it's a boundary cell (inner or outer)
                    ipoints.append(ipoint)
        self.infile.close()
        return ipoints

    def _read_header(self):
        line = self.infile.readline()
        ndm, nzones, npatches = convert_to_int(line)
        zones = [z for z in range(nzones)]

        npoints = []
        nfaces = []
        ncells = []
        mxppfs = []
        mxfpcs = []
        for nz in zones:
            line = self.infile.readline()
            (npoint, nface, ncell, mxppf, mxfpc) = convert_to_int(line)
            sline = line.strip().split()

            npoints.append(npoint)
            nfaces.append(nface)
            ncells.append(ncell)
            mxppfs.append(mxppf)
            mxfpcs.append(mxfpc)
            print("npoint=%s nface=%s ncells=%s mxppf=%s mxfpc=%s" % (npoint, nface, ncell, mxppf, mxfpc))

        self.ndm = ndm
        self.npatches = npatches
        self.npoints = npoints
        self.zones = zones
        self.nfaces = nfaces
        return (zones, npoints, nfaces, ncells, mxppfs, mxfpcs)

    def read_avus_grid(self, avus_filename):
        self.infilename = avus_filename
        self.infile = open(self.infilename, 'r')
        (zones, npoints, nfaces, ncells, mxppfs, mxfpcs) = self._read_header()
        points = self._read_points(zones, npoints)
        faces = self._read_faces(nfaces)

        self.nodes = points
        self.tri_elements = array(faces, dtype='int32') - 1

    @property
    def nnodes(self):
        return self.nodes.shape[0]

    @property
    def nelements(self):
        return self.tri_elements.shape[0]


def main():
    grid = AvusGrid()
    grid.read_avus_grid('avus.grd')

if __name__ == '__main__':
    main()




