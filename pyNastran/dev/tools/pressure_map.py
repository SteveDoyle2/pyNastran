# import os
from itertools import count
from collections import defaultdict
from typing import cast
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
                 method: str='full_model',
                 scale: float=1.0,
                 pressure_sid: int=1,
                 idtype: str='int32',
                 fdtype: str='float64',
                 pressure_filename: PathLike='pressure.bdf',
                 regions_to_include=None,
                 regions_to_remove=None) -> BDF:
    """
    Parameters
    ----------
    aero_filename : PathLike
        the path to the aero model
    nastran_filename : PathLike
        the path to the nastran model
    eid_csv_filename : PathLike
        ???
    eids_structure : int np.ndarray
        the element ids to map
    eid_load_id: int; default=0
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
    if regions_to_include is None:
        regions_to_include = []
    if regions_to_remove is None:
        regions_to_remove = []

    stop_on_failure = True
    structure_model = read_bdf(nastran_filename)
    assert len(structure_model.elements), structure_model.get_bdf_stats()
    assert len(structure_model.nodes), structure_model.get_bdf_stats()
    structure_eids = get_structural_eids_from_csv_load_id(
        structure_model,
        eids_structure, eid_csv_filename, eid_load_id,
        idtype=idtype)
    aero_model, unused_variables = get_aero_model(
        aero_filename, aero_format,
        stop_on_failure=stop_on_failure)

    if method == 'panel_model':
        pressure_model = pressure_map_to_panel_model(
            aero_model, structure_model, structure_eids,
            pressure_sid=pressure_sid,
            # aero_format=aero_format,
            map_type=map_type, scale=scale,
            # idtype=idtype,
            fdtype=fdtype)
    else:
        assert method == 'full_model', method
        pressure_model = pressure_map_to_structure_model(
            aero_model, structure_model, structure_eids,
            pressure_sid=pressure_sid,
            #aero_format=aero_format,
            map_type=map_type, scale=scale,
            #idtype=idtype,
            fdtype=fdtype)
    if pressure_filename:
        pressure_model.write_bdf(pressure_filename)
    return pressure_model


def pressure_map_to_panel_model(aero_model: Cart3D | Tecplot,
                                structure_model: BDF,
                                structure_eids: np.ndarray,
                                structure_sid=0,
                                map_type: str='pressure',
                                scale: float=1.0,
                                pressure_sid: int=1,
                                reference_point=None,
                                sref: float=1.,
                                cref: float=1.,
                                bref: float=1.,
                                #idtype: str='int32',
                                fdtype: str='float64') -> BDF:
    """Maps pressure from an aero model to a structural model. assumes:
     - 3d aero model
     - 2d structure model panel model
    """
    aero_format = aero_model.__class__.__name__.lower()
    assert aero_format in {'cart3d', 'fluent', 'tecplot'}, aero_format
    assert aero_format in {'cart3d', 'fluent', 'tecplot'}, aero_format
    assert map_type in {'pressure', 'force', 'force_moment'}, aero_format
    assert isinstance(scale, float), scale
    if reference_point is None:
        reference_point = np.zeros(3, dtype=fdtype)

    map_location = 'centroid'
    #stop_on_failure = True
    log = aero_model.log
    #---------------------------------------------------------------

    aero_dict = get_aero_pressure_centroid(
        aero_model, aero_format, map_type)
    aero_xyz_nodal = aero_dict['xyz_nodal']
    aero_Cp = aero_dict['Cp_centroidal']
    aero_area = aero_dict['area']
    aero_centroid = aero_dict['centroid']
    aero_normal = aero_dict['normal']
    aero_dr = aero_centroid - reference_point
    aero_force_coeff = (aero_Cp * aero_area)[:, np.newaxis] * aero_normal
    #aero_force = q * sref * aero_force_coeff
    aero_moment_coeff = np.cross(aero_dr, aero_force_coeff)
    total_aero_force_coeff = aero_force_coeff.sum(axis=0)
    total_aero_moment_coeff = aero_moment_coeff.sum(axis=0)
    assert len(total_aero_force_coeff) == 3, total_aero_force_coeff
    log.info(f'reference_point = {reference_point}')
    log.info(f'total_aero_force_coeff = {total_aero_force_coeff}')
    log.info(f'total_aero_moment_coeff = {total_aero_moment_coeff}')

    nelements_structure = len(structure_eids)
    log = structure_model.log
    structure_nodes, structure_xyz = get_structure_xyz(structure_model)
    out = get_mapped_structure(
        structure_model, structure_eids, fdtype=fdtype)
    structure_eids, structure_areas, structure_centroids, structure_normals = out

    abs_normals = np.abs(structure_normals)
    imax = np.argmax(abs_normals, axis=1)
    imax_min = imax.min()
    # log.info(f'imax = {imax.tolist()}')
    assert len(imax) == nelements_structure
    assert imax_min > 0, imax

    iy = np.where(imax == 1)[0]
    iz = np.where(imax == 2)[0]
    eids_y = structure_eids[iy]
    eids_z = structure_eids[iz]
    assert len(eids_y) or len(eids_z), (eids_y, eids_z)

    model2 = BDF(log=structure_model.log)
    ny = len(iy)
    nz = len(iz)
    log.info(f'ny={ny} nz={nz}')
    if nz > 0:
        aero_force_coeff_z = aero_force_coeff[iz, :]
        structure_box_centroid = structure_centroids[iz, :]
        structure_box_moment_center = np.zeros((nz, 3), dtype=fdtype)

        for i, eid, aero_forcei in zip(count(), eids_z, aero_force_coeff):
            elem = structure_model.elements[eid]
            assert elem.type == 'CQUAD4', elem
            inids = np.searchsorted(structure_nodes, elem.nodes)
            quad_xyz = structure_xyz[inids, ]
            quad_xy = quad_xyz[:, [0, 1]]
            p1 = quad_xyz[0, :]
            p4 = quad_xyz[1, :]
            p3 = quad_xyz[2, :]
            p2 = quad_xyz[3, :]
            # print(f'p1={p1}')
            # print(f'p2={p2}')
            # print(f'p3={p3}')
            # print(f'p4={p4}')
            p14 = (p1 + p4) / 2.
            p23 = (p2 + p3) / 2.
            d14 = (p4 - p1)[1:]
            span = np.linalg.norm(d14)
            chord = abs(p23[0] - p14[0])
            x12 = p2[0] - p1[0]
            x43 = p2[0] - p1[0]
            # print(f'x12={x12}')
            assert x12.min() >= 0., x12
            assert x43.min() >= 0., x43
            assert chord > 0, chord
            assert span > 0, span
            # for p1, p2, p3, p4, moment_center in zip(p1, p2, p3, p4, moment_center):
            structure_box_moment_center[i] = p14 + chord/4.

        # for each aero element, get the closest structure box
        assert len(structure_box_centroid) > 0, structure_box_centroid
        assert len(structure_box_centroid) > 0, aero_xyz_nodal
        tree = _get_tree(structure_box_centroid)
        print(structure_box_centroid.shape, aero_xyz_nodal.shape)
        if map_location == 'centroid':
            unused_deq, ieq = tree.query(aero_xyz_nodal, k=1)
        else:
            raise RuntimeError(map_location)

        nnodes = len(structure_box_centroid)
        print(unused_deq, len(unused_deq), nnodes)
        print(ieq, len(ieq))
        slots = np.where(ieq[:, :] < nnodes)

        # irows: aero?
        # icols: structure?
        irows, icols = slots
        iaero = icols
        istructure = irows
        if map_type == 'pressure':
            aero_pressure_centroidal = aero_Cp_centroidal[iaero] * scale
            mapped_structure_elements = structure_nodes[istructure]
            for eid, pressure in zip(mapped_structure_elements, aero_pressure_centroidal):
                model2.add_pload2(pressure_sid, eids=[eid], pressure=pressure)
        elif map_type == 'force_moment':
            map_force_moment_centroid_tri(
                model2,
                aero_area, aero_Cp_centroidal, iaero,
                structure_box_centroid, structure_box_moment_center, istructure,
                pressure_sid=prssure_sid,
                scale=scale,
            )
        else:
            raise RuntimeError(map_type)

    return model2


def map_pressure_centroid_avg(model2: BDF,
                              aero_area, aero_Cp_centroidal, iaero,
                              structure_nodes, istructure,
                              pressure_sid: int=1,
                              scale: float=1.0) -> None:
    """maps the pressure at the centroid of the panel (not the average)"""
    aero_pressure_centroidal = aero_Cp_centroidal[iaero] * scale
    mapped_structure_elements = structure_nodes[istructure]
    mapped_structure_area = structure_area[istructure]

    pressures = []
    for eid, pressure in zip(mapped_structure_elements, aero_pressure_centroidal):
        model2.add_pload2(pressure_sid, eids=[eid], pressure=pressure)
    return

def map_pressure_centroid(model2: BDF,
                          aero_area, aero_Cp_centroidal, iaero,
                          structure_nodes, istructure,
                          pressure_sid: int=1,
                          scale: float=1.0) -> None:
    """maps the pressure at the centroid of the panel (not the average)"""
    aero_pressure_centroidal = aero_Cp_centroidal[iaero] * scale
    mapped_structure_elements = structure_nodes[istructure]
    for eid, pressure in zip(mapped_structure_elements, aero_pressure_centroidal):
        model2.add_pload2(pressure_sid, eids=[eid], pressure=pressure)
    return


def map_force_centroid_tri(model2: BDF,
                           aero_area, aero_Cp_centroidal, iaero,
                           structure_nodes, istructure,
                           pressure_sid: int=1,
                           scale: float=1.0) -> None:
    aero_pressure_centroidal = aero_Cp_centroidal[iaero] * scale
    aero_force_centroidal = aero_pressure_centroidal * aero_area
    # pressure_nodal = Cp_nodal[iaero] * scale
    mapped_structure_nodes = structure_nodes[istructure]
    # mapped_structure_xyz = structure_xyz[istructure, :]
    forces_temp = defaultdict(array3)
    nforce_total = defaultdict(int)
    for nid, force in zip(mapped_structure_nodes, aero_force_centroidal):
        forces_temp[nid] += force
        nforce_total[nid] += 1

    mag = 1.0
    for nid, force in forces_temp.items():
        nforce_totali = nforce_total[nid]
        force2 = force / nforce_totali
        model2.add_force(pressure_sid, nid, mag, force2, cid=0)
    return


def map_force_moment_centroid_tri(model2: BDF,
                                  aero_area, aero_Cp_centroidal, iaero,
                                  structure_nodes, structure_xyz, istructure,
                                  pressure_sid: int=1,
                                  scale: float=1.0) -> None:
    aero_pressure_centroidal = aero_Cp_centroidal[iaero] * scale
    aero_force_centroidal = aero_pressure_centroidal * aero_area
    mapped_structure_nodes = structure_nodes[istructure]
    mapped_structure_xyz = structure_xyz[istructure],

    forces_temp = defaultdict(array3)
    nforce_total = defaultdict(int)
    for nid, xyz_str, force in zip(mapped_structure_nodes, mapped_structure_xyz,
                                   aero_force_centroidal):
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
    return


def pressure_map_to_structure_model(aero_model: Cart3D | Tecplot,
                                    structure_model: BDF,
                                    structure_eids: np.ndarray,
                                    structure_sid=0,
                                    map_type: str='pressure',
                                    scale: float=1.0,
                                    pressure_sid: int=1,
                                    #idtype: str='int32',
                                    fdtype: str='float64', ) -> BDF:
    """Maps pressure from an aero model to a structural model. assumes:
     - 3d aero model
     - 3d structure model
    """
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
    out = get_mapped_structure(
        structure_model, structure_eids, fdtype=fdtype)
    structure_eids, structure_areas, structure_centroids, structure_normals = out

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
        map_pressure_centroid(
            model2,
            aero_area, aero_Cp_centroidal, iaero,
            structure_nodes, istructure,
            pressure_sid=pressure_sid,
            scale=scale,
        )
    elif map_type == 'force_moment':
        map_force_moment_centroid_tri(
            model2,
            aero_area, aero_Cp_centroidal, iaero,
            structure_nodes, istructure,
            pressure_sid=pressure_sid,
            scale=scale,
        )
    elif map_type == 'force':
        map_force_centroid_tri(
            model2,
            aero_area, aero_Cp_centroidal, iaero,
            structure_nodes, istructure,
            pressure_sid=pressure_sid,
            scale=scale,
        )
    else:  # pragma: no cover
        raise RuntimeError(map_type)
    return model2


def array3() -> np.ndarray:
    return np.zeros(3, dtype='float64')
