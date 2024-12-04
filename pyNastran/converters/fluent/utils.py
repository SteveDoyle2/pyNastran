from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
if TYPE_CHECKING:
    from pyNastran.utils import PathLike
    from pyNastran.converters.fluent.fluent import Fluent


def filter_by_region(model: Fluent,
                     regions_to_remove: list[int],
                     regions_to_include: list[int]) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                                             np.ndarray, np.ndarray]:
    tris = model.tris
    quads = model.quads
    results = model.results
    assert len(model.result_element_id) > 0, model.result_element_id

    tri_regions = tris[:, 1]
    quad_regions = quads[:, 1]

    # we reordered the tris/quads to be continuous to make them easier to add
    iquad = np.searchsorted(model.result_element_id, quads[:, 0])
    itri = np.searchsorted(model.result_element_id, tris[:, 0])
    #-----------------------------
    is_remove = (len(regions_to_remove) == 0)
    is_include = (len(regions_to_include) == 0)
    assert (is_remove and not is_include) or (not is_remove and is_include)
    if regions_to_remove:
        itri_regions = np.logical_and.reduce([(tri_regions != regioni) for regioni in regions_to_remove])
        iquad_regions = np.logical_and.reduce([(quad_regions != regioni) for regioni in regions_to_remove])
    else:
        itri_regions = np.logical_or.reduce([(tri_regions == regioni) for regioni in regions_to_include])
        iquad_regions = np.logical_or.reduce([(quad_regions == regioni) for regioni in regions_to_include])

    quad_results = results[iquad, :][iquad_regions, :]
    tri_results = results[itri, :][itri_regions, :]

    tris = tris[itri_regions, :]
    quads = quads[iquad_regions, :]
    tri_eids = tris[:, 0]
    quad_eids = quads[:, 0]
    element_id = np.unique(np.hstack([quad_eids, tri_eids]))
    return element_id, tris, quads, quad_results, tri_results


def write_daten(daten_filename: PathLike,
                element_id: np.ndarray,
                titles: np.ndarray,
                results: np.ndarray) -> None:
    assert len(titles) > 1, titles
    titles = titles[1:]
    ntitles = len(titles)
    title_str = '# Shell Id, ' + ', '.join(list(titles)) + '\n'
    fmt = '%-8d' + ' %15s' * ntitles + '\n'
    with open(daten_filename, 'w') as daten_file:
        daten_file.write(title_str)
        for eid, res in zip(element_id, results.tolist()):
            daten_file.write(fmt % tuple([eid] + res))

def read_daten(daten_filename: PathLike,
               scale: float=1.0) -> tuple[np.ndarray, np.ndarray]:
    """
    # Shell Id, Cf Components X, Cf Components Y, Cf Components Z, Pressure Coefficient
         1       0.00158262927085      -0.00013470219276       0.00033890815696      -0.09045402854837
   4175456       0.00005536470364      -0.00050913009274       0.00226408640311      -0.44872838534457

    """
    with open(daten_filename, 'r') as vrt_file:
        lines = vrt_file.readlines()

    element_id_list = []
    data_list = []
    for i, line in enumerate(lines[1:]):
        sline = line.strip().split()
        eid, *datai = sline
        element_id_list.append(eid)
        data_list.append(datai)
    element_id = np.array(element_id_list, dtype='int32')

    titles_sline1 = lines[0].split('#')[1].strip().split(',')
    titles_sline2 = [val.strip() for val in titles_sline1]
    titles = np.array(titles_sline2)

    results = np.array(data_list, dtype='float64')
    if scale != 1.0:
        results *= scale

    # we drop the element_id from the check
    assert results.shape[1] == len(titles)-1, f'shape={str(results.shape)}; titles={titles}'
    return element_id, titles, results

def write_vrt(vrt_filename: PathLike, node_id: str, xyz: np.ndarray) -> None:
    title_str = 'PROSTAR_VERTEX\n'
    title_str += '4000         0         0         0         0         0         0         0\n'
    fmt = '%-8d' + ' %15s' * 3 + '\n'
    with open(vrt_filename, 'w') as vrt_file:
        vrt_file.write(title_str)
        for nid, (x, y, z) in zip(node_id, xyz):
            vrt_file.write(fmt % (nid, x, y, z))
    return

def read_vrt(vrt_filename: PathLike) -> tuple[np.ndarray, np.ndarray]:
    """
     PROSTAR_VERTEX
     4000         0         0         0         0         0         0         0
             1     10.86546244027553      0.15709716073078      1.16760397346681
       2118954      8.91578650948743      1.88831278193532      1.27715047467296
    """
    with open(vrt_filename, 'r') as vrt_file:
        lines = vrt_file.readlines()

    node_id_list = []
    xyz_list = []
    for i, line in enumerate(lines[2:]):
        sline = line.strip().split()
        nid, *xyzi = sline
        node_id_list.append(nid)
        xyz_list.append(xyzi)
    node_id = np.array(node_id_list, dtype='int32')
    xyz = np.array(xyz_list, dtype='float64')
    return node_id, xyz

def write_cell(vrt_filename: PathLike, quads: np.ndarray, tris: np.ndarray) -> None:
    title_str = 'PROSTAR_CELL\n'
    title_str += '4000         0         0         0         0         0         0         0\n'
    dim_fmt  = '%-8d %8s %8s %8d %8s\n'
    tri_fmt  = '%-8d %8s %8s %8s\n'
    quad_fmt = '%-8d %8s %8s %8s %8s\n'

    ntri, ntri_col = tris.shape
    nquad, nquad_col = quads.shape
    assert ntri_col == 5, tris.shape
    assert nquad_col == 6, quads.shape
    with open(vrt_filename, 'w') as cel_file:
        # row1 = [eid, dim, pid, 3, 4]
        # row2 = [eid, n1, n2, n3, n4]
        cel_file.write(title_str)
        for eid, pid, n1, n2, n3 in zip(tris[:, 0], tris[:, 1],
                                        tris[:, 2], tris[:, 3], tris[:, 4]):
            # {eid}      {dim}      {pid}          3          4
            # {eid}        407        406        472        473  #  nnodes = dim

            #cel_file.write(f'{eid}          3            {dim}         {pid}         4\n'
            #               f'{eid}          {n1}          {n2}          {n3}\n')
            cel_file.write(dim_fmt % (eid, '3', '3', pid, '4'))
            cel_file.write(tri_fmt % (eid, n1, n2, n3))

        for eid, pid, n1, n2, n3, n4 in zip(quads[:, 0], quads[:, 1],
                                            quads[:, 2], quads[:, 3], quads[:, 4], quads[:, 5]):
            cel_file.write(dim_fmt % (eid, '3', '4', pid, '4'))
            cel_file.write(quad_fmt % (eid, n1, n2, n3, n4))
    return

def read_cell(cell_filename: PathLike) -> tuple[tuple[np.ndarray, np.ndarray],
                                                np.ndarray]:
    """
    PROSTAR_CELL
     4000         0         0         0         0         0         0         0
         1          3          3          3          4
         1          3          4          5
   ...
   4175456          3          3         27          4
   4175456    2118943    2118954    2118941

    """
    with open(cell_filename, 'r') as cell_file:
        lines = cell_file.readlines()

    element_ids_list = []
    quad_list = []
    tri_list = []
    i = 2
    while i < len(lines):
        sline1 = lines[i].strip().split()
        sline2 = lines[i+1].strip().split()

        idi1, *ids1 = sline1
        idi2, *ids2 = sline2
        assert idi1 == idi2, (sline1, sline2)
        assert len(ids1) == 4, ids1

        #790          3          3          3          4
        #790        406        403        472
        #791          3          4          3          4
        #791        407        406        472        473
        #
        #{eid}      {dim}      {pid}          3          4
        #{eid}        407        406        472        473  #  nnodes = dim
        dim = int(ids1[1])
        pid = ids1[2]
        assert ids1[0] == '3', (idi1, ids1)
        assert ids1[1] in '34', (idi1, ids1)
        #assert ids1[2] in '34', (idi1, ids1)
        assert ids1[3] == '4', (idi1, ids1)
        nids1 = len(ids1)
        nids2 = len(ids2)
        assert nids2 == dim, (nids2, dim)
        ids = [idi1, pid] + list(ids2)
        element_ids_list.append(idi1)
        if dim == 3:
            tri_list.append(ids)
        elif dim == 4:
            quad_list.append(ids)
        else:  # pragma: no cover
            raise RuntimeError(dim)
        assert nids1 in {4}, ids1
        assert nids2 in {3, 4}, ids2
        i += 2

    quads = np.array(quad_list, dtype='int32')
    tris = np.array(tri_list, dtype='int32')
    if len(tri_list):
        assert tris.shape[1] == 2+3, (tris.shape, quads.shape)
    if len(quad_list):
        assert quads.shape[1] == 2+4, (tris.shape, quads.shape)

    element_ids = np.array(element_ids_list, dtype='int32')
    return (quads, tris), element_ids
