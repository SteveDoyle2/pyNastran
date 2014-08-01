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
    #bdf.write('LOAD = %s\n' % load_id)
    bdf.write('BEGIN BULK\n')
    for (x, y, z) in model.nodes:
        card = ['GRID', nid, cid, x, y, z]
        bdf.write(print_card(card))

        #nx, ny, nz = nodal_normals[nid - 1]
        #card = ['FORCE', load_id, nid, cid, 100, nx, ny, nz]
        #bdf.write(print_card(card))
        nid += 1

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
    bdf_filename = 'threePlugs.bdf'
    stl_filename = 'threePlugs.stl'
    stl_to_nastran_filename(stl_filename, bdf_filename)