from struct import pack, unpack
from numpy import array, transpose, zeros
from pyNastran.utils.log import get_logger

def write_front(model):
    pass

def write_face(model):
    pass

class Usm3dReader(object):
    def __init__(self, log=None, debug=None):
        self.nodes = None
        self.tris = None
        self.tets = None
        self.precision = 'double'
        self.log = get_logger(log, 'debug' if debug else 'info')

    def write_usm3d(self, basename):
        if 1:
            write_usm3d_volume(self, basename)
        else:
            asdf

    def read_bc_format(self, bc_filename):
        """
        0 - Supersonic Inflow
        1 - Reflection plane
        2 - Supersonic Outflow
        3 - Subsonic Outer Boundaries
        4 - Viscous Surfaces
        5 - Inviscid aerodynamic surface
        44 - Blunt base
        55 - Thick Trailing Edges
        n*100+3 - Engine-exhaust (Fan)
        n*100+2 - Engine-exhaust (Jet Core)
        n*100+1 - Engine-intake
        1001 - Special inflow
        1002 - Special Outflow (Fixed Pressure)

        0 - Freestream - Supersonic Inflow (Bounding Box)
        2 - Extrapolation - Supersonic Outflow (Bounding Box)

        1 - Reflection Plane - Tangent Flow - (Symmetry Plane)
        3 - Characteristic Inflow - Subsonic Inflow/Outflow/Sideflow (Bounding Box)

        4 - Inviscid Surface (Physical)
        5 - Viscous Surface (Physical)
        #==============================

        #Thu Dec 19 11:46:03 2013
        #bc.map
        Patch #        BC             Family   #surf   surfIDs         Family
        #----------------------------------------------------------------------
        1              44             44             0            0        Base  -> Blunt base
        2              4              4              0            0        Bay   -> Viscous Surfaces
        3              0              0              0            0        Inflow -> Supersonic Inflow
        4              2              2              0            0        Outflow -> Supersonic Outflow
        5              3              3              0            0        Sideflow -> Characteristic Inflow/Outflow
        7              1              1              0            0        Symmetry -> Reflection plane
        """
        bc = open(bc_filename, 'r')
        lines = bc.readlines()
        lines2 = []
        for line in lines:
            if len(line.strip().split('#')[0]) > 0:
                lines2.append(line)
        lines = lines2

        for line in lines:
            patch_id, bc, family, surf, surf_ids,  = line.split()
        bc.close()

    def read_usm3d(self, basename, dimension_flag):
        cogsg_file = basename + '.cogsg'
        face_file = basename + '.face'
        front_file = basename + '.front'

        return self.read_cogsg(cogsg_file)
        #self.read_front(front_file)
        #self.read_face(face_file)

    def read_cogsg(self, cogsg_file):
        f = open(cogsg_file, 'rb')

        # nelements * 4 * 4 + 32 ???
        dummy = f.read(4)  # 1022848
        dummy_int, = unpack('>i', dummy)
        #assert dummy_int == 1022848, 'dummy_int = %s' % dummy_int

        # file header
        if self.precision=='single':
            Format = '>6if'
            nbytes = 6 * 4 + 4
            self.log.debug("***single precision")
        elif self.precision=='double':
            Format = '>6id'
            nbytes = 6 * 4 + 8
        else:
            raise Exception('invalid precision format')
        data = f.read(nbytes)

        (inew, ne, np, nb, npv, nev, tc) = unpack(Format, data)
        self.header = {
            'dummy'    : dummy_int,
            'inew'     : inew, # dummy int
            'nElements': ne,  # nc;  number of tets
            'nPoints'  : np,  # npo; number of grid points including nbn
            'nBoundPts': nb,  # nbn; number of boundary points including nbc
            'nViscPts' : npv, # npv; number of viscous points (=0 for Euler)
            'nViscElem': nev, # ncv; number of viscous cells (=0 for Euler)
            'tc'       : tc,  # dummy double
                              # nbc
        }
        self.log.info(self.header)

        # nbn nodes
        #
        #del ne, np

        if 1:
            return self._read_cogsg_volume(f)
        #else:
        #----------------------------------------------------------------------
        # elements
        # face elements
        nnodes_per_face = 3
        nfaces = ne

        if nfaces > 0:
            data_length = nnodes_per_face * nfaces
            Format = '>' + 'i' * data_length
            data = f.read(4 * data_length)

            faces = unpack(Format, data)
            faces = array(faces)
            faces = faces.reshape((nfaces, 3))
        else:
            faces = None
        #----------------------------------------------------------------------
        # nodes
        nBoundPts = nb
        nnodes = nBoundPts

        #data_length = nnodes
        if self.precision == 'double':
            data_length = 8 * nnodes
        elif self.precision == 'single':
            data_length = 4 * nnodes
        else:
            raise RuntimeError('precision = %r' % self.precision)

        skip_nodes = False
        if skip_nodes == True:
            t = self.tell()
            f.goto(t + data_length * 3)
            nodes = None
        else:
            if self.precision == 'double':
                Format = '>%sd' % nnodes
                node_array_format = 'float64'
            elif self.precision == 'single':
                Format = '>%sd' % nnodes
                node_array_format = 'float32'
            else:
                raise RuntimeError('precision = %r' % self.precision)

            data = f.read(data_length)
            X = unpack(Format, data)

            data = f.read(data_length)
            Y = unpack(Format, data)

            data = f.read(data_length)
            Z = unpack(Format, data)

            nodes = array([X, Y, Z])

            nodes = zeros((nnodes, 3), node_array_format)
            nodes[:, 0] = X
            nodes[:, 1] = Y
            nodes[:, 2] = Z
            del X, Y, Z


        f.read(nnodes * 3 * 8)  # 3 -> xyz, 8 -> double precision ???

        #----------------------------------------------------------------------
        # elements
        # boundary layer elements
        nnodes_per_tet = 4
        ntets = nev

        if ntets:
            data_length = nnodes_per_tet * ntets
            Format = '>' + 'i' * data_length
            data = f.read(4 * data_length)

            tets = unpack(Format, data)
            tets = array(tets)
            tets = tets.reshape((tets, 4))

        #----------------------------------------------------------------------
        # volume points
        nnodes = npv

        Format = '>%si' % nnodes
        data = f.read(4 * nnodes)

        nodes_vol = unpack(Format, data)
        nodes_vol = array(nodes_vol)
        nodes_vol = nodes_vol.reshape((tets, 3))



    def _read_cogsg_volume(self, f):
        # volume cells
        self.log.debug('tell volume = %s' % f.tell())
        # surface + volume cells ???
        nelements = self.header['nElements']
        Format = '>%si' % nelements

        elements = zeros((nelements, 4), 'int32')

        self.log.debug("fv.tell = %s" % f.tell())
        for i in range(4): #  tets
            data = f.read(4 * nelements)
            elements[:, i] = unpack(Format, data)
            #print "elements[:, %s] =" %i, elements[:, i]
        elements -= 1

        #print 'tell volume2 =', f.tell()
        dummy2 = f.read(4)
        self.log.debug("dummy2 = %s %s" % (unpack('>i', dummy2), unpack('>f', dummy2)))
        dummy_int2, = unpack('>i', dummy2)

        # 32 = dummy_int2 - 4 * nelements * 4
        assert self.header['dummy'] == dummy_int2

        #-----------------------------------
        # nodes
        nnodes = self.header['nPoints']
        Format = '>%sd' % nnodes

        dummy3 = f.read(4)  # nnodes * 3 * 8
        dummy3_int, = unpack('>i', dummy3)
        #assert dummy3_int == 298560
        print "dummy3 = ", unpack('>i', dummy3), unpack('>f', dummy3)

        nodes = zeros((nnodes, 3), 'float64')
        for i in range(3): #  x, y, z
            data = f.read(8 * nnodes)
            assert len(data) == (8 * nnodes)
            nodes[:, i] = unpack(Format, data)

        dummy4 = f.read(4) # nnodes * 3 * 8
        dummy4_int, = unpack('>i', dummy4)
        #print "dummy4 = ", unpack('>i', dummy4), unpack('>f', dummy4)

        assert dummy3_int == dummy4_int
        self.nodes = nodes
        self.tets = elements
        return nodes, elements





def write_usm3d_volume(model, basename):
    cogsg_file = basename + '.cogsg'
    face_file = basename + '.face'
    front_file = basename + '.front'

    write_cogsg_volume(model, cogsg_file)
    #write_front(model, front_file)
    #write_face(model, face_file)

def write_cogsg_volume(model, cogsg_file):
    #n = 0
    self = model
    nnodes, three = self.nodes.shape
    ntets, four = self.tets.shape

    outfile = open(cogsg_file, 'wb')

    # file header

    values = [32 + ntets * 4 * 4,]
    block_size = pack('>i', *values)
    outfile.write(block_size)

    header = self.header
    values_header = [header['inew'],
              header['nElements'],
              header['nPoints'],
              header['nBoundPts'],
              header['nViscPts'],
              header['nViscElem'],
              header['tc'], # d
              ]
    #n = 36
    Format = '>6id'
    model.log.debug("Format = %r" % Format)
    data = pack(Format, *values_header)
    outfile.write(data)
    #print "outfile.tell = ", outfile.tell()
    #outfile.write(block_size)
    #--------------------------------------------------------------------------
    tets = self.tets

    # tet header
    #values = [ntets * 4 * 4]  # n1, n2, n3, n4 -> 4; 4 -> int
    #block_size = pack('>i', *values)
    #outfile.write(block_size)

    # tets
    Format = '>%si' % ntets
    n0 = tets[:, 0] + 1
    n1 = tets[:, 1] + 1
    n2 = tets[:, 2] + 1
    n3 = tets[:, 3] + 1
    #print "n0 =", n0
    outfile.write( pack(Format, *n0) )
    outfile.write( pack(Format, *n1) )
    outfile.write( pack(Format, *n2) )
    outfile.write( pack(Format, *n3) )
    #n += 4 * 8 * ntets
    #print "outfile.tell 2 = ", outfile.tell(), n

    # tet footer
    outfile.write(block_size)
    #--------------------------------------------------------------------------

    # nodes header
    values = [nnodes * 3 * 8]  # xyz -> 3; 8 -> double precision (???)
    block_size = pack('>i', *values)
    outfile.write(block_size)

    # nodes
    #npoints = header['nPoints']
    Format = '>%sd' % nnodes
    nodes = self.nodes
    n0 = nodes[:, 0]
    n1 = nodes[:, 1]
    n2 = nodes[:, 2]
    outfile.write( pack(Format, *n0) )
    outfile.write( pack(Format, *n1) )
    outfile.write( pack(Format, *n2) )

    # nodes footer
    outfile.write(block_size)

    #--------------------------------------------------------------------------
    outfile.close()


def main():
    model = Usm3dReader()
    #basename = 'HSCT_inviscid'
    basename = 'box'
    model.read_usm3d(basename, 3)
    model.write_usm3d(basename + '_2')

    model.read_usm3d(basename + '_2', 3)


if __name__ == '__main__':
    main()