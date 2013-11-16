from struct import pack, unpack
from numpy import array, transpose


def write_front(model):
    pass

def write_face(model):
    pass

class Usm3d(object):
    def __init__(self):
        self.nodes = None
        self.tets = None
        self.precision = 'double'

    def write_usm3d(self, basename):
        write_usm3d(self, basename)

    def read_usm3d(self, basename):
        cogsg_file = basename + '.cogsg'
        face_file = basename + '.face'
        front_file = basename + '.front'

        self.read_cogsg(cogsg_file)
        self.read_front(front_file)
        self.read_face(face_file)

    def read_cogsg(self, cogsg_file):
        f = open(cogsg_file, 'rb')

        # file header
        if(self.precision=='single'):
            Format = '>6if'
            nbytes = 6 * 4 + 4
            self.log.debug("***single precision")
        elif(self.precision=='double'):
            Format = '>6id'
            nbytes = 6 * 4 + 8
        else:
            raise Exception('invalid precision format')
        data = f.read(nbytes)

        (inew, ne, np, nb, npv, nev, tc) = unpack(data, Format)
        self.header = {
            'inew'     : inew,
            'nElements': ne,
            'nPoints'  : np,
            'nBoundPts': nb,
            'nViscPts' : npv,
            'nViscElem': nev,
            'tc'       : tc,
        }

        #----------------------------------------------------------------------
        # elements
        # face elements
        nnodes_per_face = 3
        nfaces = ne

        data_length = nnodes_per_face * nfaces
        Format = '>' + 'i' * data_length
        data = f.read(4 * data_length)

        faces = unpack(data, Format)
        faces = array(faces)
        faces = faces.reshape((nfaces, 3))

        #----------------------------------------------------------------------
        # nodes
        nBoundPts = nb
        nnodes = nBoundPts

        #data_length = nnodes
        if self.precision == 'dobule':
            data_length = 8 * nnodes
        else:
            data_length = 4 * nnodes

        if skip_nodes == True:
            t = self.tell()
            f.goto(t + data_length * 3)
            nodes = None
        else:
            if self.precision == 'dobule':
                Format = '>' + 'd' * nnodes
                node_array_format = 'float64'
            else:
                Format = '>' + 'd' * nnodes
                node_array_format = 'float32'
            data = f.read(data_length)
            X = unpack(data, Format)

            data = f.read(data_length)
            Y = unpack(data, Format)

            data = f.read(data_length)
            Z = unpack(data, Format)

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

        data_length = nnodes_per_tet * ntets
        Format = '>' + 'i' * data_length
        data = f.read(4 * data_length)

        tets = unpack(data, Format)
        tets = array(tets)
        tets = tets.reshape((tets, 4))

        #----------------------------------------------------------------------
        # volume points
        nnodes = npv

        Format = '>' + 'i' * nnodes
        data = f.read(4 * nnodes)

        tets = unpack(data, Format)
        tets = array(tets)
        tets = tets.reshape((tets, 4))









def write_usm3d(model, basename):
    cogsg_file = basename + '.cogsg'
    face_file = basename + '.face'
    front_file = basename + '.front'

    write_cogsg(model, cogsg_file)
    write_front(model, front_file)
    write_face(model, face_file)

def write_cogsg(model, cogsg_file):
    self = model
    nnodes, three = self.nodes.shape
    ntets, four = self.tets.shape

    f = open(cogsg_file, 'wb')

    # file header

    values = [32,]
    block_size = struct.pack('>i', *values)
    f.write(block_size)

    values = [header['inew'],
              header['nElements'],
              header['nPoints'],
              header['nBoundPts'],
              header['nViscPts'],
              header['nViscElem'],
              header['tc'], # d
              ]
    Format = '>6id'
    self.log().debug("Format = %r" % Format)
    data = struct.pack(Format, *values)
    f.write(data)
    f.write(block_size)
    #--------------------------------------------------------------------------

    # nodes header
    block_size = [nnodes * 3 * 8]  # xyz -> 3; 8 -> double precision (???)
    blockSize = struct.pack('>i', *values)
    outfile.write(block_size)

    # nodes
    Format = '>'+'d'*nPoints
    n0 = nodes[:, 0]
    n1 = nodes[:, 1]
    n2 = nodes[:, 2]
    outfile.write( struct.pack(Format, *n0) )
    outfile.write( struct.pack(Format, *n1) )
    outfile.write( struct.pack(Format, *n2) )

    # nodes footer
    outfile.write(block_size)

    #--------------------------------------------------------------------------

    # tets header

    # tets

    f.close()



def main():
    from tetgen_reader import TetgenReader
    m = TetgenReader()
    base = 'tetgen_test_flipped.1'
    m.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=3)
    m.write_nastran(base + '.bdf')

    write_usm3d(base)
