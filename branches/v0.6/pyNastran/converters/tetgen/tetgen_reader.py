from numpy import array, zeros
from pyNastran.utils.log import get_logger


class TetgenReader(object):
    def __init__(self, log=None, debug=False):
        self.log = get_logger(log, 'debug' if debug else 'info')
        self.nodes = None
        self.tri = None

    def read_tetgen(self, node_filename, smesh_filename):
        self.nodes = self.read_nodes(node_filename)
        self.tri = self.read_smesh(smesh_filename)

    def read_smesh(self, smesh_filename):
        f = open(smesh_filename, 'r')
        for i in xrange(6):
            f.readline()

        nelements, zero = f.readline().strip().split() # nelements, 0
        nelements = int(nelements)


        #print "line =", line
        tri_list = []
        for ielement in xrange(nelements):
            sline = f.readline().strip().split()
            nnodes = sline[0]
            element_nodes = sline[1:]
            if nnodes == '3':
                tri_list.append(element_nodes)
            else:
                raise NotImplementedError('nnodes = %s' % nnodes)
            #print element_nodes
        tri = array(tri_list, 'int32') - 1 # subtract 1 so the node ids start at 0
        print tri
        f.close()
        return tri

    def read_nodes(self, node_filename):
        f = open(node_filename, 'r')
        nnodes, three, zero1, zero2 = f.readline().strip().split()
        assert three == '3', three
        assert zero1 == '0', zero1
        assert zero2 == '0', zero2
        nnodes = int(nnodes)
        nodes = zeros((nnodes, 3), 'float64')
        for inode in xrange(nnodes):
            nodes[inode] = f.readline().strip().split()[1:]
        f.close()
        print "nodes =", nodes
        return nodes


def main():
    m = SMesh_Reader()
    m.read_tetgen('tetgen_test_flipped.1.node', 'tetgen_test_flipped.1.smesh')

if __name__ == '__main__':
    main()