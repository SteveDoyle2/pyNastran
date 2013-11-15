from numpy import array, zeros
from pyNastran.utils.log import get_logger


class TetgenReader(object):
    def __init__(self, log=None, debug=False):
        self.log = get_logger(log, 'debug' if debug else 'info')
        self.nodes = None
        self.tri = None
        self.tet = None

    def read_tetgen(self, node_filename, smesh_filename, ele_filename, dimension_flag):
        self.nodes = self.read_nodes(node_filename)
        if dimension_flag == 2:
            self.tri = self.read_smesh(smesh_filename)
        elif dimension_flag == 3:
            self.tet = self.read_ele(ele_filename)
        else:
            raise RuntimeError('dimension_flag = %r and must be 2 or 3.' % dimension_flag)

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

    def read_ele(self, ele_filename):
        print "ele_filename =", ele_filename
        f = open(ele_filename, 'r')
        nelements, four, one = f.readline().strip().split()

        assert four == '4', four
        assert one == '1', one
        nelements = int(nelements)
        print "nelements =", nelements

        tets = zeros((nelements, 4), 'float64')
        for ielement in xrange(nelements):
            # eid n1    n2    n3    n4       flip_flag???
            # 1   13260 15506 16059 16065    -1

            tets[ielement] = f.readline().strip().split()[1:5]
        f.close()
        #print "nodes =", nodes
        return tets - 1
    #self.tet = self.read_ele(ele_filename)


def main():
    m = SMesh_Reader()
    m.read_tetgen('tetgen_test_flipped.1.node', 'tetgen_test_flipped.1.smesh')

if __name__ == '__main__':
    main()