import os
from collections import defaultdict
from typing import cast, Any
import numpy as np

from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.bdf_equivalence import _get_tree

from pyNastran.converters.cart3d.cart3d import Cart3D
#from pyNastran.converters.fluent.fluent import read_fluent
from pyNastran.converters.tecplot.tecplot import Tecplot

from pyNastran.dev.tools.pressure_map_structure_setup import (
    get_structural_eids_from_csv_load_id, get_structure_xyz, get_mapped_structure)
from pyNastran.dev.tools.pressure_map_aero_setup import get_aero_model, get_aero_pressure_centroid


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

def pressure_map_from_model(aero_model: Cart3D | Tecplot,
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
    #stop_on_failure = True
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
