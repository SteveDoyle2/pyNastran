import os
from pathlib import Path
import numpy as np


class Fluent:
    def __init__(self, log=None):
        self.log = log

    def read_fluent(self, fluent_filename: str):
        base, ext = os.path.splitext(fluent_filename)
        #case = os.path.basename(vrt_filename)
        vrt_filename = base + '.vrt'
        daten_filename = base + '.daten'
        cell_filename = base + '.cel'
        #fld_filename = base + '.fld'
        assert os.path.exists(vrt_filename)
        assert os.path.exists(daten_filename), daten_filename
        assert os.path.exists(cell_filename), cell_filename

        (quads, tris), (element_ids, region, elements_list) = read_cell(cell_filename)
        node, xyz = read_vrt(vrt_filename)
        element_id, pressure, titles, results = read_daten(daten_filename, scale=1.0)

        #tri_centroid, tri_pressure = tri_split(xyz, tris, element_id, results)
        #quad_centroid, quad_pressure = quad_split(xyz, quads, element_id, results)
        self.node_id = node
        self.xyz = xyz
        self.element_id = element_id
        self.results = pressure
        self.quads = quads
        self.tris = tris

        self.element_ids = element_ids
        self.region = region
        self.elements_list = elements_list
        #self.tri_pressure = tri_


def read_daten(daten_filename: Path,
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
    pressure = results[:, -1]
    if scale != 1.0:
        pressure *= scale
    return element_id, pressure, titles, results

def read_vrt(vrt_filename):
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

def read_cell(cell_filename):
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
    regions_list = []
    elements_list = []

    quad_list = []
    tri_list = []
    i = 2
    while i < len(lines):
        sline1 = lines[i].strip().split()
        sline2 = lines[i+1].strip().split()
        #print(sline1)
        #print(sline2)
       # print('---')

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
        regions_list.append(pid)
        elements_list.append(ids2)
        if dim == 3:
            tri_list.append(ids)
        elif dim == 4:
            quad_list.append(ids)
        else:
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
    regions = np.array(regions_list, dtype='int32')
    return (quads, tris), (element_ids, regions, elements_list)

def tri_split(xyz: np.ndarray,
              tris: np.ndarray,
              element_id: np.ndarray,
              results: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # [eid, pid, n1, n2, n3]
    eids = tris[:, 0]
    n1 = tris[:, 2] - 1
    n2 = tris[:, 3] - 1
    n3 = tris[:, 4] - 1
    xyz1 = xyz[n1, :]
    xyz2 = xyz[n2, :]
    xyz3 = xyz[n3, :]
    tri_centroid = (xyz1 + xyz2 + xyz3) / 3
    ieid = np.searchsorted(element_id, eids)
    tri_pressure = results[ieid]
    return tri_centroid, tri_pressure

def quad_split(xyz: np.ndarray,
               quads: np.ndarray,
               element_id: np.ndarray,
               results: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    eids = quads[:, 0]
    n1 = quads[:, 2] - 1
    n2 = quads[:, 3] - 1
    n3 = quads[:, 4] - 1
    n4 = quads[:, 5] - 1
    xyz1 = xyz[n1, :]
    xyz2 = xyz[n2, :]
    xyz3 = xyz[n3, :]
    xyz4 = xyz[n4, :]
    quad_centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4
    ieid = np.searchsorted(element_id, eids)
    quad_pressure = results[ieid]
    return quad_centroid, quad_pressure

def read_fluent(fluent_filename: str, log=None, debug: bool=False) -> Fluent:
    model = Fluent(log=log)
    model.read_fluent(fluent_filename)
    return model
