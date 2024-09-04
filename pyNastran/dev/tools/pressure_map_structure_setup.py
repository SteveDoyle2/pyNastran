import os
import numpy as np

from pyNastran.utils import PathLike, print_bad_path
from pyNastran.bdf.bdf import BDF


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
