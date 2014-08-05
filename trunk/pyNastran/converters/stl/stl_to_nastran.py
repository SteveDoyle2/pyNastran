from pyNastran.converters.stl.stl_reader import STLReader
from pyNastran.bdf.fieldWriter import print_card

def stl_to_nastran_filename(stl_filename, bdf_filename,
                            nnodes_offset=0, nelements_offset=0, log=None):
    model = STLReader(log=log)
    model.read_stl(stl_filename)

    nid = nnodes_offset + 1
    cid = None
    pid = 100
    mid = 200
    load_id = 10

    nodal_normals = model.get_normals_at_nodes(model.elements)

    bdf = open(bdf_filename, 'wb')
    bdf.write('CEND\n')
    bdf.write('LOAD = %s\n' % load_id)
    bdf.write('BEGIN BULK\n')
    nid2 = 1
    magnitude = 100.
    for (x, y, z) in model.nodes:
        card = ['GRID', nid, cid, x, y, z]
        bdf.write(print_card(card))

        nx, ny, nz = nodal_normals[nid2 - 1]
        card = ['FORCE', load_id, nid, cid, magnitude, nx, ny, nz]
        bdf.write(print_card(card))
        nid += 1
        nid2 += 1

    eid = nelements_offset + 1
    for (n1, n2, n3) in (model.elements + (nnodes_offset + 1)):
        card = ['CTRIA3', eid, pid, n1, n2, n3]
        bdf.write(print_card(card))
        eid += 1

    t = 0.1
    card = ['PSHELL', pid, mid, t]
    bdf.write(print_card(card))

    E = 1e7
    G = None
    nu = 0.3
    card = ['MAT1', mid, E, G, nu]
    bdf.write(print_card(card))

    bdf.write('ENDDATA\n')
    bdf.close()

if __name__ == '__main__':
    import os

    import pyNastran
    root_path = pyNastran.__path__[0]
    print "root_path =",  root_path

    from pyNastran.converters.cart3d.cart3d_to_stl import cart3d_to_stl_filename

    cart3d_filename = os.path.join(root_path, 'converters', 'cart3d', 'threePlugs_bin.tri')
    stl_filename = os.path.join(root_path, 'converters', 'stl', 'threePlugs.stl')
    bdf_filename = os.path.join(root_path, 'converters', 'stl', 'threePlugs.bdf')

    cart3d_to_stl_filename(cart3d_filename, stl_filename)
    stl_to_nastran_filename(stl_filename, bdf_filename)