from collections import defaultdict
import numpy as np
from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.bdf_equivalence import _get_tree

from pyNastran.converters.cart3d.cart3d import read_cart3d
#from pyNastran.converters.fluent.fluent import read_fluent
#from pyNastran.converters.tecplot.tecplot import read_tecplot


def pressure_map(aero_filename: PathLike,
                 nastran_filename: PathLike,
                 structure_eids=np.array([]),
                 structure_sid=0,
                 aero_format: str='cart3d',
                 map_type: str='pressure',
                 scale: float=1.0,
                 output_sid: int=1):
    """
    Parameters
    ----------
    structure_eids: int np.ndarray
        the elemeent ids to map
    structure_sid: int; default=0
        the sid for the element ids to map
    aero_filename: str
        path to input aero filename:
          - cart3d: triq
          - fluent: vrt, cel, daten
    nastran_filename:  str
        path to output nastran filename
    aero_format: str
        cart3d
    map_type:   str
        pressure
    scale : float
        qinf

    """
    idtype = 'int32'
    fdtype = 'float64'
    assert aero_format in {'cart3d', 'fluent', 'tecplot'}, aero_format
    assert map_type in {'pressure', 'force', 'force_moment'}, aero_format
    assert isinstance(scale, float), scale

    if map_type == 'pressure':
        map_location = 'centroid'
    else:
        map_location = 'node'

    if aero_format == 'cart3d':
        assert map_type in {'pressure', 'force', 'force_moment'}, aero_format
        aero_model = read_cart3d(aero_filename)
        aero_elements = aero_model.elements
        aero_xyz_nodal = aero_model.nodes
        xyz1 = aero_xyz_nodal[aero_elements[:, 0], :]
        xyz2 = aero_xyz_nodal[aero_elements[:, 1], :]
        xyz3 = aero_xyz_nodal[aero_elements[:, 2], :]
        aero_xyz_centroid = (xyz1 + xyz2 + xyz3) / 3
        aero_normal = np.cross(xyz2-xyz1, xyz3-xyz1)
        aero_area = np.linalg.norm(aero_normal)
        aero_Cp_centroidal = aero_model.loads['Cp']
        assert len(aero_Cp_centroidal) == len(aero_elements)
        assert len(aero_Cp_centroidal) == len(aero_xyz)
    #elif aero_format == 'tecplot':
    #    aero_model = read_tecplot(aero_filename)
    #elif aero_format == 'fluent':
    #    aero_model = read_fluent(aero_filename)
    else:  # pragma: no cover
        raise RuntimeError(aero_format)

    structure_model = read_bdf(nastran_filename)
    structure_nodes, structure_xyz = get_structure_xyz(structure_model)
    structure_eids, structure_centroids = get_mapped_structure_element_ids(
        structure_model, structure_eids, structure_sid,
        idtype=idtype, fdtype=fdtype)

    tree = _get_tree(aero_xyz_nodal)
    if map_location == 'centroid':
        unued_deq, ieq = tree.query(structure_centroids, k=1)
    elif map_location == 'node':
        unued_deq, ieq = tree.query(structure_xyz, k=4)
    else:  # pragma: no cover
        raise RuntimeError(map_location)

    nnodes = len(aero_xyz_nodal)
    slots = np.where(ieq[:, :] < nnodes)

    # irows: aero?
    # icols: structure?
    irows, icols = slots
    iaero = irows
    istructure = icols

    #mapped_structure_elements = structure_nodes[istructure]

    model2 = BDF(log=structure_model.log)
    if map_type == 'pressure':
        aero_pressure_centroidal = aero_Cp_centroidal[iaero] * scale
        mapped_structure_elements = structure_nodes[istructure]
        for eid, pressure in zip(mapped_structure_elements, aero_pressure_centroidal):
            model2.add_pload2(output_sid, eids=[eid], pressure=pressure)
    elif map_type == 'force':
        aero_pressure_centroidal = aero_Cp_centroidal[iaero] * scale
        aero_force_centroidal = aero_pressure_centroidal * aero_area
        #pressure_nodal = Cp_nodal[iaero] * scale
        mapped_structure_nodes = structure_nodes[istructure]
        #mapped_structure_xyz = structure_xyz[istructure, :]
        forces_temp = defaultdict(array3)
        nforce_total = defaultdict(int)
        for nid, force in zip(mapped_structure_nodes, aero_force_centroidal):
            forces_temp[nid] += force
            nforce_total[nid] += 1
        #forces = {}
        mag = 1.0
        for nid, force in forces_temp.items():
            nforce_totali = nforce_total[nid]
            force2 = force / nforce_totali
            model2.add_force(output_sid, nid, mag, force2, cid=0)
    else:
        raise RuntimeError(map_type)
    return model2


def array3() -> np.ndarray:
    return np.zeros(3, dtype='float64')


def get_element_ids_by_sid(structure_model: BDF,
                           structure_sid: int,
                           idtype: str='int32') -> tuple[np.ndarray]:
    """
    Parameters
    ----------
    structure_sid: int; default
        the sid for the element ids to map

    Returns
    -------
    structure_eids: int np.ndarray
        the element ids to map

    """
    unsupported_loads = set()
    loads = structure_model.loads[structure_sid]
    eids = set()
    for load in loads:
        if load.type == 'PLOAD4':
            eid = load.eid
            eids.add(eid)
        elif load.type == 'PLOAD2':
            eids.update(load.eids)
        else:
            unsupported_loads.add(load.type)

    structure_eids = np.array(list(eids), dtype=idtype)
    structure_eids.sort()
    if len(unsupported_loads):
        structure_model.log.warning(f'unsupported_loads = {unsupported_loads}')
    return structure_eids

def get_structure_xyz(structure_model: BDF) -> tuple[np.ndarray, np.ndarray]:
    (nid_cp_cd, xyz_cid0,
     xyz_cp, unused_icd_transform, unused_icp_transform,
     ) = structure_model.get_xyz_in_coord_array()
    #del xyz_cp, icd_transform, icp_transform
    structure_nodes = nid_cp_cd[:, 0]
    structure_xyz = xyz_cid0
    return structure_nodes, structure_xyz
def get_mapped_structure_element_ids(structure_model: BDF,
                                     structure_eids: np.ndarray,
                                     structure_sid: int = 0,
                                     idtype: str='int32',
                                     fdtype: str='float64') -> tuple[np.ndarray, np.ndarray]:
    """
    Parameters
    ----------
    structure_eids: int np.ndarray
        the element ids to map
    structure_sid: int; default=0
        the sid for the element ids to map

    Returns
    -------
    structure_eids: int np.ndarray
        the element ids to map
    structure_centroids : float np.ndarray
        the associated centroids

    """
    if len(structure_eids) and structure_sid > 0:
        raise RuntimeError('too many choices for element ids')
    elif len(structure_eids) == 0 and structure_sid == 0:
        raise RuntimeError('no elements were passed in')

    if structure_sid > 0:
        structure_eids = get_element_ids_by_sid(
            structure_model, structure_sid, idtype=idtype)
    #else:

    centroids = []
    for eid in structure_eids:
        elem = structure_model.elements[eid]
        assert elem.type in {'CTRIA3', 'CQUAD4', 'CTRIA6', 'CQUAD8', 'CQUAD'}, elem
        centroid = elem.Centroid()
        centroids.append(centroid)
    structure_centroids = np.array(centroids, dtype=fdtype)
    return structure_eids, structure_centroids
