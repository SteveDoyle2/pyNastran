from copy import deepcopy

from pyNastran.converters.aflr.ugrid.ugrid2d_reader import UGRID2D_Reader
from pyNastran.bdf.field_writer_8 import print_card_8, print_int_card

def ugrid2d_to_nastran_filename(ugrid2d_filename, bdf_filename,
                                axis_order=None,
                                nid_start=1, eid_start=1,
                                cp=2000, pid=2000, mid=2000,
                                punch=False):
    model = UGRID2D_Reader()
    model.read_ugrid(ugrid2d_filename)

    convert_ugrid2d_to_nastran(bdf_filename, model.nodes, model.tris, model.quads,
                               cp=cp, pid=pid, mid=mid,
                               axis_order=axis_order,
                               nid_start=nid_start, eid_start=eid_start, punch=punch)


def convert_ugrid2d_to_nastran(bdf_filename, nodes, tris, quads,
                               cp=2000, pid=2000, mid=2000,
                               axis_order=None,
                               nid_start=1, eid_start=1, punch=True):
    if axis_order is not None:
        nodes = deepcopy(nodes[:, axis_order])

    with open(bdf_filename, 'wb') as bdf_file:
        if not punch:
            bdf_file.write('CEND\n')
            bdf_file.write('BEGIN BULK\n')

        cp = None
        #card = ['CORD2R', cp, 0] + [0., 0., 0.] + [0., 0., 1.] + [1., 0., 0.]
        #f.write(print_card_8(card))

        nid = nid_start
        for xi, yi, zi in nodes:
            # yes I'm aware...x is the axial distance
            # y/z are the x/y plane
            card = ['GRID', nid, cp, xi, yi, zi]
            bdf_file.write(print_card_8(card))
            nid += 1


        t = 0.1
        card = ['PSHELL', pid, mid, t]
        bdf_file.write(print_card_8(card))

        E = 3.0E7
        G = None
        nu = 0.3
        card = ['MAT1', mid, E, G, nu]
        bdf_file.write(print_card_8(card))

        eid = eid_start
        for n1, n2, n3 in tris + nid_start:
            card = ['CTRIA3', eid, pid, n1, n2, n3]
            bdf_file.write(print_int_card(card))
            eid += 1

        for n1, n2, n3, n5 in quads + nid_start:
            card = ['CQUAD4', eid, pid, n1, n2, n3, n5]
            bdf_file.write(print_int_card(card))
            eid += 1

        if not punch:
            bdf_file.write('ENDDATA\n')
