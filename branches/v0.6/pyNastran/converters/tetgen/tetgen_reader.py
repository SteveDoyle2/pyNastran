from numpy import array, zeros
from pyNastran.utils.log import get_logger
from pyNastran.bdf.fieldWriter import print_card


class TetgenReader(object):
    """
    http://www.wias-berlin.de/preprint/1762/wias_preprints_1762.pdf
    """
    def __init__(self, log=None, debug=False):
        self.log = get_logger(log, 'debug' if debug else 'info')
        self.nodes = None
        self.tri = None
        self.tet = None

    def write_nastran(self, bdf_filename):
        f = open(bdf_filename, 'wb')
        msg = 'CEND\n'
        msg += 'BEGIN BULK\n'

        nid = 1
        cid = None
        for (x, y, z) in self.nodes:
            card = ['GRID', nid, cid, x, y, z]
            msg += print_card(card)
            nid += 1
        f.write(msg)

        eid = 1
        mid = 100
        if self.tri is not None:
            pid = 1
            thickness = 0.1
            pshell = ['PSHELL', pid, mid, thickness]
            msg = print_card(pshell)

            for (n0, n1, n2) in (self.tri + 1):
                card = ['CTRIA3', eid, pid, n0, n1, n2]
                msg += print_card(card)
                eid += 1
                if eid % 1000 == 0:
                    f.write(msg)
                    msg = ''
            f.write(msg)

        if self.tet is not None:
            pid = 2
            psolid = ['PSOLID', pid, mid]
            msg = print_card(psolid)
            for (n0, n1, n2, n3) in (self.tet + 1):
                card = ['CTETRA', eid, pid, n0, n1, n2, n3]
                msg += print_card(card)
                eid += 1
                if eid % 1000 == 0:
                    f.write(msg)
                    msg = ''
            f.write(msg)

        E = 1e7
        G = None
        nu = 0.3
        rho = 0.1
        mat1 = ['MAT1', mid, E, G, nu, rho]
        msg = print_card(mat1)
        f.write(msg)
        f.write('ENDDATA\n')

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

    def read_ele(self, ele_filename, form_flag='1'):
        #print "ele_filename =", ele_filename
        f = open(ele_filename, 'r')
        nelements, four, form_flag_enabled = f.readline().strip().split()
        form_flag_enabled = int(form_flag_enabled)

        assert four == '4', four
        nelements = int(nelements)
        #print "nelements =", nelements

        if not form_flag_enabled:
            tets = zeros((nelements, 4), 'int32')
            for ielement in xrange(nelements):
                # eid n1    n2    n3    n4       flip_flag???
                # 1   13260 15506 16059 16065    -1
                tets[ielement] = f.readline().strip().split()[1:]
        else:
            tets = []
            for ielement in xrange(nelements):
                # eid n1    n2    n3    n4       flip_flag???
                # 1   13260 15506 16059 16065    -1
                (n0, n1, n2, n3, flag) = f.readline().strip().split()[1:]
                if flag == form_flag:
                    tets.append( (n0, n1, n2, n3) )
            tets = array(tets, 'int32')

        f.close()
        #print "nodes =", nodes
        return tets - 1
    #self.tet = self.read_ele(ele_filename)


def main():
    import os
    from pyNastran.converters.stl.stl_reader import STLReader

    m1 = STLReader()
    m1.read_stl('tetgen_test.stl')
    m1.flip_normals()
    m1.write_stl('tetgen_test_flipped.stl')
    del m1

    os.system('tetgen.exe -pqcvVqY tetgen_test_flipped.stl')

    m = TetgenReader()
    base = 'tetgen_test_flipped.1'
    m.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=3)
    m.write_nastran(base + '.bdf')

if __name__ == '__main__':
    main()