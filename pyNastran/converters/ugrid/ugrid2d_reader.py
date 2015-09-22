from copy import deepcopy
from numpy import zeros, array
from pyNastran.bdf.field_writer_8 import print_card_8, print_int_card

def ugrid2d_to_nastran_filename(ugrid2d_filename, bdf_filename,
                                axis_order=None,
                                nid_start=1, eid_start=1,
                                cp=2000, pid=2000, mid=2000,
                                punch=False):
    u = UGRID2D_Reader()
    u.read_ugrid(ugrid2d_filename)

    convert_ugrid2d_to_nastran(bdf_filename, u.nodes, u.tris, u.quads,
                               cp, pid, mid,
                               axis_order=axis_order,
                               nid_start=nid_start, eid_start=eid_start, punch=punch)

class UGRID2D_Reader(object):
    def __init__(self, log=None, debug=None):
        pass

    def read_ugrid(self, ugrid_filename):
        f = open(ugrid_filename, 'r').read()
        data = f.split()
        #f.close()

        nnodes, ntrias, nquads, ntets, npyram5, npenta6, nhexas8s = (int(val) for val in data[:7])
        i = 7

        print('nnodes=%s ntrias=%s nquads=%s ntets4=%s npyram5=%s npenta6=%s nhexas8s=%s' % (
            nnodes, ntrias, nquads, ntets, npyram5, npenta6, nhexas8s))


        #nodes = zeros(nnodes * 3, dtype=ndarray_float)
        #tris  = zeros(ntris * 3, dtype='int32')
        #quads = zeros(nquads * 4, dtype='int32')
        #pids = zeros(npids, dtype='int32')

        #tets = zeros(ntets * 4, dtype='int32')
        #penta5s = zeros(npenta5s * 5, dtype='int32')
        #penta6s = zeros(npenta6s * 6, dtype='int32')
        #hexas = zeros(nhexas * 8, dtype='int32')

        # nodes
        iend = i + nnodes * 3
        nodes = array(data[i:iend], dtype='float64')
        nodes = nodes.reshape((nnodes, 3))
        print(nodes[0, :])
        print(nodes[nnodes-1, :])
        print(nodes[-1, :])
        assert nodes[:, 2].max() == 0.0
        assert nodes[:, 2].min() == 0.0
        i = iend

        # tris
        iend = i + ntrias * 3
        tris = array(data[i:iend], dtype='int32')
        tris = tris.reshape((ntrias, 3))
        #print(nodes[0, :])
        #print(nodes[nnodes-1, :])
        #print(nodes[-1, :])
        #assert nodes[:, 2].max() == 0.0
        #assert nodes[:, 2].min() == 0.0
        i = iend

        iend = i + nquads * 4
        quads = array(data[i:iend], dtype='int32')
        quads = quads.reshape((nquads, 4))
        #print(nodes[0, :])
        #print(nodes[nnodes-1, :])
        #print(nodes[-1, :])
        #assert nodes[:, 2].max() == 0.0
        #assert nodes[:, 2].min() == 0.0
        i = iend

        self.nodes = nodes
        self.tris = tris - 1
        self.quads = quads - 1


def convert_ugrid2d_to_nastran(nastran_filename, nodes, tris, quads,
                               cp, pid, mid,
                               axis_order=None,
                               nid_start=1, eid_start=1, punch=True):
    if axis_order is not None:
        nodes = deepcopy(nodes[:, axis_order])

    f = open(nastran_filename, 'wb')
    if not punch:
        f.write('CEND\n')
        f.write('BEGIN BULK\n')

    cp = None
    #card = ['CORD2R', cp, 0] + [0., 0., 0.] + [0., 0., 1.] + [1., 0., 0.]
    #f.write(print_card_8(card))

    nid = nid_start
    for xi, yi, zi in nodes:
        # yes I'm aware...x is the axial distance
        # y/z are the x/y plane
        card = ['GRID', nid, cp, xi, yi, zi]
        f.write(print_card_8(card))
        nid += 1


    t = 0.1
    card = ['PSHELL', pid, mid, t]
    f.write(print_card_8(card))

    E = 3.0E7
    G = None
    nu = 0.3
    card = ['MAT1', mid, E, G, nu]
    f.write(print_card_8(card))

    eid = eid_start
    for n1, n2, n3 in tris + nid_start:
        card = ['CTRIA3', eid, pid, n1, n2, n3]
        f.write(print_int_card(card))
        eid += 1

    for n1, n2, n3, n5 in quads + nid_start:
        card = ['CQUAD4', eid, pid, n1, n2, n3, n5]
        f.write(print_int_card(card))
        eid += 1

    if not punch:
        f.write('ENDDATA\n')
