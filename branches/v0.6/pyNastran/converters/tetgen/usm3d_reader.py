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
        self.tets = None
        self.precision = 'double'
        self.log = get_logger(log, 'debug' if debug else 'info')

    def write_usm3d(self, basename):
        if 1:
            write_usm3d_volume(self, basename)
        else:
            asdf

    def read_usm3d(self, basename):
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
        assert dummy_int == 1022848, 'dummy_int = %s' % dummy_int

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

        # nbn nodes
        #
        #del ne, np

        if 1:
            return self._cogsg_volume(f)
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



    def _cogsg_volume(self, f):
        # volume cells

        # surface + volume cells ???
        nelements = self.header['nElements']
        Format = '>%si' % nelements

        elements = zeros((nelements, 4), 'int32')

        for i in range(4): #  tets
            data = f.read(4 * nelements)
            elements[:, i] = unpack(Format, data)
        dummy2 = f.read(4)
        print "dummy2 =", unpack('>i', dummy2), unpack('>f', dummy2)

        dummy_int2, = unpack('>i', dummy2)

        # 32 = dummy_int2 - 4 * nelements * 4
        assert self.header['dummy'] == dummy_int2

        #-----------------------------------
        # nodes
        nnodes = self.header['nPoints']
        Format = '>%sd' % nnodes

        f.read(4)

        nodes = zeros((nnodes, 3), 'float64')
        for i in range(3): #  x, y, z
            data = f.read(8 * nnodes)
            assert len(data) == (8 * nnodes)
            nodes[:, i] = unpack(Format, data)

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
    self = model
    nnodes, three = self.nodes.shape
    ntets, four = self.tets.shape

    outfile = open(cogsg_file, 'wb')

    # file header

    values = [32,]
    block_size = pack('>i', *values)
    outfile.write(block_size)

    header = self.header
    values = [header['inew'],
              header['nElements'],
              header['nPoints'],
              header['nBoundPts'],
              header['nViscPts'],
              header['nViscElem'],
              header['tc'], # d
              ]
    Format = '>6id'
    model.log.debug("Format = %r" % Format)
    data = pack(Format, *values)
    outfile.write(data)
    outfile.write(block_size)
    #--------------------------------------------------------------------------
    tets = self.tets

    # tet header
    values = [ntets * 4 * 4]  # n1, n2, n3, n4 -> 4; 4 -> int
    block_size = pack('>i', *values)
    outfile.write(block_size)

    # tets
    Format = '>%sd' % ntets
    n0 = tets[:, 0]
    n1 = tets[:, 1]
    n2 = tets[:, 2]
    n3 = tets[:, 3]
    outfile.write( pack(Format, *n0) )
    outfile.write( pack(Format, *n1) )
    outfile.write( pack(Format, *n2) )
    outfile.write( pack(Format, *n3) )

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

    # tets header

    # tets

    outfile.close()


def main():
    model = Usm3dReader()
    #basename = 'HSCT_inviscid'
    basename = 'box'
    model.read_usm3d(basename)
    model.write_usm3d(basename + '_2')


if __name__ == '__main__':
    main()