from pathlib import PurePath
from pyNastran.converters.stl.stl import STL
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.utils import PathLike

def stl_to_nastran_filename(stl_filename: PathLike | STL,
                            bdf_filename: PathLike,
                            nnodes_offset: int=0, nelements_offset: int=0,
                            pid: int=100, mid: int=200,
                            size: int=8, is_double: bool=False,
                            log=None) -> None:
    return stl_to_nastran(stl_filename, bdf_filename, nnodes_offset=nnodes_offset,
                          nelements_offset=nelements_offset,
                          pid=pid, mid=mid, size=size,
                          is_double=is_double, log=log)

def stl_to_nastran(stl_filename: PathLike | STL,
                   bdf_filename: PathLike,
                   nnodes_offset: int=0, nelements_offset: int=0,
                   pid: int=100, mid: int=200,
                   size: int=8, is_double: bool=False,
                   log=None) -> None:
    if isinstance(stl_filename, (str, PurePath)):
        model = STL(log=log)
        model.read_stl(stl_filename)
    elif isinstance(stl_filename, STL):
        model = stl_filename
    else:
        raise TypeError('stl_filename must be a string or STL; type=%s' % type(stl_filename))

    nid = nnodes_offset + 1
    cid = None
    unused_load_id = 10

    #nodal_normals = model.get_normals_at_nodes(model.elements)

    if size == 8:
        print_card = print_card_8
    elif size == 16:
        if is_double:
            print_card = print_card_16
        else:
            print_card = print_card_double
    else:  # pragma: no cover
        raise RuntimeError('size=%r' % size)

    with open(bdf_filename, 'w') as bdf_file:
        bdf_file.write('CEND\n')
        #bdf.write('LOAD = %s\n' % load_id)
        bdf_file.write('BEGIN BULK\n')
        nid2 = 1
        unused_magnitude = 100.

        for x, y, z in model.nodes:
            card = ['GRID', nid, cid, x, y, z]
            bdf_file.write(print_card_16(card))

            #nx, ny, nz = nodal_normals[nid2 - 1]
            #card = ['FORCE', load_id, nid, cid, magnitude, nx, ny, nz]
            #bdf.write(print_card_8(card))
            nid += 1
            nid2 += 1

        eid = nelements_offset + 1

        elements = model.elements + (nnodes_offset + 1)
        for (n1, n2, n3) in elements:
            card = ['CTRIA3', eid, pid, n1, n2, n3]
            bdf_file.write(print_card_8(card))
            eid += 1

        t = 0.1
        card = ['PSHELL', pid, mid, t]
        bdf_file.write(print_card_8(card))

        E = 1e7
        G = None
        nu = 0.3
        card = ['MAT1', mid, E, G, nu]
        bdf_file.write(print_card_8(card))
        bdf_file.write('ENDDATA\n')
    return
