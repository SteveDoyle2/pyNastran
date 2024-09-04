import os
from collections import defaultdict
from typing import cast, Any
import numpy as np

from pyNastran.utils import PathLike, print_bad_path
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.bdf_equivalence import _get_tree

from pyNastran.converters.cart3d.cart3d import read_cart3d, Cart3D
#from pyNastran.converters.fluent.fluent import read_fluent
from pyNastran.converters.tecplot.tecplot import read_tecplot, Tecplot


def get_aero_model(aero_filename: str, aero_format: str,
                   stop_on_failure=True) -> tuple[Any, list[str]]:
    assert os.path.exists(aero_filename), print_bad_path(aero_filename)
    if aero_format == 'Cart3D':
        model = read_cart3d(aero_filename)
        variables = list(model.loads)
    # elif aero_format == 'Fund3D':
    #     return None, []
    # elif aero_format == 'Fluent Vrt':
    #     pass
    # elif aero_format == 'Fluent Press':
    #     pass
    elif aero_format == 'Tecplot':
        model = read_tecplot(aero_filename)
        #print(model.object_stats())
        variables = model.result_variables
    else:  # pragma: no cover
        if stop_on_failure:
            raise NotImplementedError(aero_format)
        else:
            return None, []
    return model, variables

def get_aero_pressure_centroid(aero_model: Cart3D | Tecplot,
                               aero_format: str,
                               map_type: str,
                               variable: str='Cp') -> dict[str, np.ndarray]:
    """
    variable: str
        cart3d: 'Cp' only; variable is ignored
        tecplot: use variable
        fun3d:        ???
        fluent press: ???
        fluent vrt:   ???
    """
    if aero_format == 'cart3d':
        aero_model = cast(Cart3D, aero_model)
        assert map_type in {'pressure', 'force', 'force_moment'}, aero_format
        aero_elements = aero_model.elements
        aero_xyz_nodal = aero_model.nodes
        xyz1 = aero_xyz_nodal[aero_elements[:, 0], :]
        xyz2 = aero_xyz_nodal[aero_elements[:, 1], :]
        xyz3 = aero_xyz_nodal[aero_elements[:, 2], :]
        #aero_xyz_centroid = (xyz1 + xyz2 + xyz3) / 3
        aero_normal = np.cross(xyz2-xyz1, xyz3-xyz1)
        aero_area = np.linalg.norm(aero_normal)
        aero_Cp_centroidal = aero_model.loads['Cp']
        assert len(aero_Cp_centroidal) == len(aero_elements)
        assert len(aero_Cp_centroidal) == len(aero_xyz)
    elif aero_format == 'tecplot':
        aero_model = cast(Tecplot, aero_model)
        #raise RuntimeError(aero_model)
    #elif aero_format == 'fluent':
    #    aero_model = read_fluent(aero_filename)
    else:  # pragma: no cover
        raise RuntimeError(aero_format)
    aero_dict = {
        'xyz_nodal' : aero_xyz_nodal,
        'Cp_centroidal' : aero_Cp_centroidal,
        'area' : aero_area,
    }
    # aero_Cp_centroidal = out_dict['aero_Cp_centroidal']
    # aero_area = out_dict['aero_area']
    return aero_dict

def get_structural_eids_from_csv_load_id(model: BDF,
                                         eids_structure: np.ndarray,
                                         csv_filename: PathLike='',
                                         load_id: int=0,
                                         idtype: str='int32') -> np.ndarray:
    """
    Parameters
    ----------
    eids_structure : np.ndarray
        the direct method to specify eids
    csv_filename : PathLike; default=''
    load_id : int; default=0

    Returns
    -------
    """
    nelements = len(eids_structure)
    # one must be true
    is_nelements = (nelements > 0)
    is_csv = len(csv_filename)
    is_load_id = (load_id > 0)
    # exclusive or / XOR
    assert (is_nelements ^ is_csv) or (is_csv ^ is_load_id), (nelements, csv_filename, load_id)

    assert nelements == '' or csv_filename == '' or load_id == 0, (nelements, csv_filename, load_id)
    assert nelements or load_id > 0, (nelements, csv_filename, load_id)

    if nelements:
        structural_eids_out = eids_structure
    elif csv_filename:
        structural_eids_out = _load_structural_eids_from_csv(csv_filename, idtype=idtype)
    elif load_id > 0:
        # load_id
        structural_eids_out = get_element_ids_by_sid(model, load_id, dtype=idtype)
    else:  # pragma: no cover
        raise RuntimeError('failed to load eids')
    structural_eids_out.sort()
    return structural_eids_out

def _load_structural_eids_from_csv(csv_filename: PathLike, idtype: str='int32'):
    assert os.path.exists(csv_filename), print_bad_path(csv_filename)
    with open(csv_filename, 'r') as csv_file:
        lines = csv_file.readlines()

    structural_eids_set = set([])
    for line in lines:
        line = line.strip().split('#')[0]
        if not line:
            continue
        sline = line.split(',')
        structural_eids_set.update(sline)
    structural_eids_out = np.array(list(structural_eids_set), dtype=idtype)
    return structural_eids_out

def pressure_map(aero_filename: PathLike,
                 nastran_filename: PathLike,
                 eids_structure=np.array([]),
                 eid_csv_filename: PathLike='',
                 eid_load_id: int=10,
                 aero_format: str='cart3d',
                 map_type: str='pressure',
                 scale: float=1.0,
                 pressure_sid: int=1,
                 idtype: str='int32',
                 fdtype: str='float64',
                 pressure_filename: PathLike='pressure.bdf') -> BDF:
    """
    Parameters
    ----------
    aero_filename : PathLike
        the path to the aero model
    nastran_filename : PathLike
        the path to the nastran model
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
    idtype : str; default='int32'
        the type of the integers
    fdtype : str; default='float64'
        the  type of the floats

    """
    assert aero_format in {'cart3d', 'fluent', 'tecplot'}, aero_format
    assert map_type in {'pressure', 'force', 'force_moment'}, aero_format
    assert isinstance(scale, float), scale

    stop_on_failure = True
    structure_model = read_bdf(nastran_filename)
    structure_eids = get_structural_eids_from_csv_load_id(
        structure_model,
        eids_structure, eid_csv_filename, eid_load_id,
        idtype=idtype)
    aero_model, unused_variables = get_aero_model(
        aero_filename, aero_format,
        stop_on_failure=stop_on_failure)

    pressure_model = pressure_map_from_model(
        aero_model, structure_model, structure_eids,
        pressure_sid=pressure_sid,
        #aero_format=aero_format,
        map_type=map_type, scale=scale,
        #idtype=idtype,
        fdtype=fdtype)
    if pressure_filename:
        pressure_model.write_bdf(pressure_filename)
    return pressure_model

def pressure_map_from_model(aero_model,
                            structure_model: BDF,
                            structure_eids: np.ndarray,
                            structure_sid=0,
                            #aero_format: str='cart3d',
                            map_type: str='pressure',
                            scale: float=1.0,
                            pressure_sid: int=1,
                            #idtype: str='int32',
                            fdtype: str='float64', ) -> BDF:

    aero_format = aero_model.__class__.__name__.lower()
    assert aero_format in {'cart3d', 'fluent', 'tecplot'}, aero_format
    assert map_type in {'pressure', 'force', 'force_moment'}, aero_format
    assert isinstance(scale, float), scale

    if map_type == 'pressure':
        map_location = 'centroid'
    else:
        map_location = 'node'
    stop_on_failure = True
    #---------------------------------------------------------------

    aero_dict = get_aero_pressure_centroid(aero_model, aero_format, map_type)
    aero_xyz_nodal = aero_dict['xyz_nodal']
    aero_Cp_centroidal = aero_dict['Cp_centroidal']
    aero_area = aero_dict['area']

    structure_nodes, structure_xyz = get_structure_xyz(structure_model)
    structure_eids, structure_centroids = get_mapped_structure(
        structure_model, structure_eids, fdtype=fdtype)

    tree = _get_tree(aero_xyz_nodal)
    if map_location == 'centroid':
        unused_deq, ieq = tree.query(structure_centroids, k=1)
    elif map_location == 'node':
        unused_deq, ieq = tree.query(structure_xyz, k=4)
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
            model2.add_pload2(pressure_sid, eids=[eid], pressure=pressure)
    elif map_type == 'force_moment':
        aero_pressure_centroidal = aero_Cp_centroidal[iaero] * scale
        aero_force_centroidal = aero_pressure_centroidal * aero_area
        mapped_structure_nodes = structure_nodes[istructure]

        forces_temp = defaultdict(array3)
        nforce_total = defaultdict(int)
        for nid, force in zip(mapped_structure_nodes, aero_force_centroidal):
            forces_temp[nid] += force
            nforce_total[nid] += 1

        #       3
        #      /|\
        #     / | \
        #    /  c  \
        #   / /   \ \
        #  1---------2
        #
        # s: semiperimeter
        a = np.linalg.norm(n2 - n1)
        b = np.linalg.norm(n3 - n2)
        c = np.linalg.norm(n1 - n3)
        s = 0.5 * (a + b + c)
        A = np.sqrt(s * (s-a) * (s-b) * (s-c))
        for nid, force in zip(mapped_structure_nodes, aero_force_centroidal):
            forces_temp[nid] += force * dist / dist_total
            nforce_total[nid] += 1
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
            model2.add_force(pressure_sid, nid, mag, force2, cid=0)
    else:  # pragma: no cover
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
    skip_loads = {
        'FORCE', 'FORCE1', 'FORCE2',
        'MOMENT', 'MOMENT1', 'MOMENT2',
        'PLOAD1', 'GRAV', 'TEMP', 'ACCEL', 'ACCEL1',
    }
    unsupported_loads = set()
    loads = structure_model.loads[structure_sid]
    eids_set = set()
    for load in loads:
        if load.type in skip_loads:
            continue

        if load.type == 'PLOAD4':
            eid = load.eid
            eids_set.add(eid)
        elif load.type == 'PLOAD':
            eids_set.update(load.eids)
        elif load.type == 'PLOAD2':
            eids_set.update(load.eids)
        else:
            unsupported_loads.add(load.type)

    structure_eids = np.array(list(eids_set), dtype=idtype)
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

def get_mapped_structure(structure_model: BDF,
                         structure_eids: np.ndarray,
                         fdtype: str='float64') -> tuple[np.ndarray, np.ndarray]:
    """
    Get element_ids and centroids

    Parameters
    ----------
    structure_eids: int np.ndarray
        the element ids to map

    Returns
    -------
    structure_eids: int np.ndarray
        the element ids to map
    structure_centroids : float np.ndarray
        the associated centroids

    """
    if len(structure_eids) == 0:
        raise RuntimeError('no elements were passed in')

    centroids = []
    for eid in structure_eids:
        elem = structure_model.elements[eid]
        assert elem.type in {'CTRIA3', 'CQUAD4', 'CTRIA6', 'CQUAD8', 'CQUAD'}, elem
        centroid = elem.Centroid()
        centroids.append(centroid)
    structure_centroids = np.array(centroids, dtype=fdtype)
    return structure_eids, structure_centroids
