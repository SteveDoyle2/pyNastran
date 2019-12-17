"""
defines:
 - read_tetgen(base, dimension_flag=2, log=None, debug=False)
 - Tetgen(log=None, debug=False):
   - write_nastran(self, bdf_filename)
   - read_tetgen(self, node_filename, smesh_filename, ele_filename, dimension_flag)
   - read_smesh(self, smesh_filename)
   - read_nodes(self, node_filename)
   - read_ele(self, ele_filename, form_flag='1')

"""
from numpy import array, zeros
from cpylog import get_logger2
from pyNastran.bdf.field_writer_8 import print_card_8


def read_tetgen(base, dimension_flag=2, log=None, debug=False):
    """simplified interface to Tetgen files"""
    model = Tetgen(log=log, debug=debug)
    model.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag)
    return model

class Tetgen:
    """
    http://www.wias-berlin.de/preprint/1762/wias_preprints_1762.pdf
    """
    def __init__(self, log=None, debug=False):
        """
        Initializes the Tetgen object

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
        self.log = get_logger2(log, debug=debug)
        self.nodes = None
        self.tris = None
        self.tets = None

    def write_nastran(self, bdf_filename):
        """writes a nastran bdf"""
        with open(bdf_filename, 'w') as bdf_file:
            msg = 'CEND\n'
            msg += 'BEGIN BULK\n'

            nid = 1
            cid = None
            for (x, y, z) in self.nodes:
                card = ['GRID', nid, cid, x, y, z]
                msg += print_card_8(card)
                nid += 1
            bdf_file.write(msg)

            eid = 1
            mid = 100
            if self.tris is not None:
                pid = 1
                thickness = 0.1
                pshell = ['PSHELL', pid, mid, thickness]
                msg = print_card_8(pshell)

                for n0, n1, n2 in self.tris + 1:
                    card = ['CTRIA3', eid, pid, n0, n1, n2]
                    msg += print_card_8(card)
                    eid += 1
                    if eid % 1000 == 0:
                        bdf_file.write(msg)
                        msg = ''
                bdf_file.write(msg)

            if self.tets is not None:
                pid = 2
                psolid = ['PSOLID', pid, mid]
                msg = print_card_8(psolid)
                for n0, n1, n2, n3 in self.tets + 1:
                    card = ['CTETRA', eid, pid, n0, n1, n2, n3]
                    msg += print_card_8(card)
                    eid += 1
                    if eid % 1000 == 0:
                        bdf_file.write(msg)
                        msg = ''
                bdf_file.write(msg)

            E = 1e7
            G = None
            nu = 0.3
            rho = 0.1
            mat1 = ['MAT1', mid, E, G, nu, rho]
            msg = print_card_8(mat1)
            bdf_file.write(msg)
            bdf_file.write('ENDDATA\n')

    def read_tetgen(self, node_filename, smesh_filename, ele_filename, dimension_flag):
        """reads a tetgen file"""
        self.nodes = read_nodes(node_filename)
        if dimension_flag == 2:
            self.log.debug('reading the *.smesh')
            self.tris = self.read_smesh(smesh_filename)
        elif dimension_flag == 3:
            self.log.debug('reading the *.ele')
            self.tets = read_ele(ele_filename)
        else:
            raise RuntimeError('dimension_flag = %r and must be 2 or 3.' % dimension_flag)

    def read_smesh(self, smesh_filename):
        """reads the *.smesh file"""
        with open(smesh_filename, 'r') as smesh_file:
            lines = smesh_file.readlines()
            lines = clean_lines(lines)
            #iline = 0 # 0  3  0  0 # node list is found in .node file.

            iline = 1
            nelements, unused_zero = lines[iline].split() # nelements, 0
            nelements = int(nelements)
            self.log.debug('nelements = %s' % nelements)

            # facet section
            tri_list = []
            iline += 1
            for unused_ielement in range(nelements):
                sline = lines[iline].split()
                try:
                    nnodes = sline[0]
                except IndexError:
                    print(sline)
                    raise
                element_nodes = sline[1:]
                if nnodes == '3':
                    tri_list.append(element_nodes)
                else:
                    raise NotImplementedError('nnodes = %s' % nnodes)
                iline += 1
            tri = array(tri_list, 'int32') - 1 # subtract 1 so the node ids start at 0
        return tri


def read_nodes(node_filename):
    """reads the *.node file"""
    with open(node_filename, 'r') as node_file:
        nnodes, three, zero1, zero2 = node_file.readline().strip().split()
        assert three == '3', three
        assert zero1 == '0', zero1
        assert zero2 == '0', zero2
        nnodes = int(nnodes)
        nodes = zeros((nnodes, 3), 'float64')
        for inode in range(nnodes):
            nodes[inode] = node_file.readline().strip().split()[1:]
    return nodes

def read_ele(ele_filename, form_flag='1'):
    """reads the *.ele file"""
    #print("ele_filename =", ele_filename)
    with open(ele_filename, 'r') as ele_file:
        nelements, four, form_flag_enabled = ele_file.readline().strip().split()
        form_flag_enabled = int(form_flag_enabled)

        assert four == '4', four
        nelements = int(nelements)
        #print("nelements =", nelements)

        if not form_flag_enabled:
            tets = zeros((nelements, 4), 'int32')
            for ielement in range(nelements):
                # eid n1    n2    n3    n4       flip_flag???
                # 1   13260 15506 16059 16065    -1
                tets[ielement] = ele_file.readline().strip().split()[1:]
        else:
            tets = []
            for ielement in range(nelements):
                # eid n1    n2    n3    n4       flip_flag???
                # 1   13260 15506 16059 16065    -1
                n0, n1, n2, n3, flag = ele_file.readline().strip().split()[1:]
                if flag == form_flag:
                    tets.append((n0, n1, n2, n3))
            tets = array(tets, 'int32')

    #print("nodes =", nodes)
    return tets - 1
    #self.tet = self.read_ele(ele_filename)


def clean_lines(lines):
    """removes blank lines and commented lines"""
    lines2 = []
    for line in lines:
        line2 = line.split('#')[0].strip()
        if line2:
            lines2.append(line2)
    return lines2


def main():  # pragma: no cover
    import os

    #base = 'gear'
    #read_tetgen(base, dimension_flag=2)
    #return

    from pyNastran.converters.stl.stl import STL
    m1 = STL()
    m1.read_stl('tetgen_test.stl')
    m1.flip_normals()
    m1.write_stl('tetgen_test_flipped.stl')
    del m1

    os.system('tetgen.exe -pqcvVqY tetgen_test_flipped.stl')

    m = Tetgen()
    base = 'tetgen_test_flipped.1'
    m.read_tetgen(base + '.node', base + '.smesh', base + '.ele', dimension_flag=3)
    m.write_nastran(base + '.bdf')

if __name__ == '__main__':  # pragma: no cover
    main()
