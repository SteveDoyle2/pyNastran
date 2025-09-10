# import os
from itertools import count, zip_longest
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
                 xyz_units: str='in',
                 pressure_units: str='psi',
                 pressure_sid: int=1,
                 force_sid: int=2,
                 moment_sid: int=3,
                 idtype: str='int32',
                 fdtype: str='float64',
                 pressure_filename: PathLike='pressure.bdf',
                 aero_xyz_scale: float=1.0,
                 qinf: float=1.0,
                 sref: float=1.0,
                 cref: float=1.0,
                 bref: float=1.0,
                 reference_point: np.ndarray | None=None,
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
    qinf : float; default=1.0
        dynamic pressure
    idtype : str; default='int32'
        the type of the integers
    fdtype : str; default='float64'
        the  type of the floats
    sref : float; default=1.0
        reference area
    bref : float; default=1.0
        reference length for Mx and Mz
    cref : float; default=1.0
        reference length for My
    aero_xyz_scale : float; default=1.0
        scales the aero model to the structural units
        1.0 if units are consistent
        39.3701.0 if aero model in meters and structural model in inches
        39,370.1  if aero model in millimeters and structural model in inches
    pressure_sid : int; default=1
        the output pressure load id
    force_sid : int; default=2
        the output force load id
    moment_sid : int; default=3
        the output moment load id
    xyz_units : str; default='in'
        the units for geometry, reference area / lengths
    pressure_units : str; default='psi'
        the units for qinf

    """
    assert aero_format in {'cart3d', 'fluent', 'tecplot'}, aero_format
    assert map_type in {'pressure', 'force', 'force_moment'}, aero_format
    assert isinstance(qinf, float), qinf
    assert isinstance(sref, float), sref
    assert isinstance(bref, float), bref
    assert isinstance(cref, float), cref
    if regions_to_include is None:
        regions_to_include = []
    if regions_to_remove is None:
        regions_to_remove = []
    if reference_point is None:
        reference_point = np.zeros(3, dtype=fdtype)

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
        aero_xyz_scale=aero_xyz_scale,
        stop_on_failure=stop_on_failure,
        xyz_units=xyz_units,
        # regions_to_remove=regions_to_remove,
        # regions_to_include=regions_to_include,
    )
    log = structure_model.log
    log.info(f'qinf ({pressure_units}) = {qinf:.3f}')

    if method == 'panel_model':
        pressure_model = pressure_map_to_panel_model(
            aero_model, structure_model, structure_eids,
            reference_point,
            regions_to_include=regions_to_include,
            regions_to_remove=regions_to_remove,
            pressure_sid=pressure_sid,
            force_sid=force_sid,
            moment_sid=moment_sid,
            # aero_format=aero_format,
            map_type=map_type,
            qinf=qinf,
            sref=sref,
            bref=bref,
            cref=cref,
            xyz_units=xyz_units,
            # idtype=idtype,
            fdtype=fdtype)
    else:
        assert method == 'full_model', method
        pressure_model = pressure_map_to_structure_model(
            aero_model, structure_model, structure_eids,
            reference_point,
            pressure_sid=pressure_sid,
            force_sid=force_sid,
            moment_sid=moment_sid,
            #aero_format=aero_format,
            map_type=map_type,
            qinf=qinf, sref=sref, bref=bref, cref=cref,
            regions_to_include=regions_to_include,
            regions_to_remove=regions_to_remove,
            #idtype=idtype,
            fdtype=fdtype)
    if pressure_filename:
        pressure_model.write_bdf(pressure_filename)
    return pressure_model


def pressure_map_to_panel_model(aero_model: Cart3D | Tecplot,
                                structure_model: BDF,
                                structure_eids: np.ndarray,
                                reference_point: np.ndarray,
                                structure_sid=0,
                                map_type: str='pressure',
                                scale: float=1.0,
                                qinf: float=1.0,
                                pressure_sid: int=1,
                                force_sid: int=2,
                                moment_sid: int=3,
                                sref: float=1.0,
                                cref: float=1.0,
                                bref: float=1.0,
                                xyz_units: str='in',
                                regions_to_include=None,
                                regions_to_remove=None,
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
    assert isinstance(qinf, float), qinf
    if reference_point is None:
        reference_point = np.zeros(3, dtype=fdtype)

    map_location = 'centroid'
    #stop_on_failure = True
    log = aero_model.log
    #---------------------------------------------------------------

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import matplotlib.cm as cm

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_box_aspect((1, 1, 1))

    #---------------------------------------------------------------
    aero_dict = get_aero_pressure_centroid(
        aero_model, aero_format, map_type,
        regions_to_include=regions_to_include,
        regions_to_remove=regions_to_remove,
    )
    #aero_xyz_nodal = aero_dict['xyz_nodal']
    # aero_Cp = aero_dict['Cp_centroid']
    # aero_area = aero_dict['area']
    # aero_centroid = aero_dict['centroid']
    # aero_normal = aero_dict['normal']

    aero_Cp = aero_dict['tri_Cp_centroid']
    aero_area = aero_dict['tri_area']
    aero_centroid = aero_dict['tri_centroid']
    aero_normal = aero_dict['tri_normal']

    aero_tri_nodes = aero_dict['tri_nodes']
    aero_xyz = aero_dict['xyz_nodal']
    aero_node_id = aero_dict['node_id']
    assert len(aero_Cp) == len(aero_centroid), (len(aero_Cp), len(aero_centroid))
    aero_dr = aero_centroid - reference_point

    # F/qinf per panel
    aero_force_coeff_per_q = (aero_Cp * aero_area)[:, np.newaxis] * aero_normal
    # M/qinf per panel
    aero_moment_coeff_per_q = np.cross(aero_dr, aero_force_coeff_per_q)

    # Cp = (p-pinf)/qinf
    # for qinf=0.
    #    Cp = p/qinf; p=qinf*Cp
    # Fi = pi * Ai * ni = qinf * Cpi * Ai * ni
    # F = qinf * sum(Cpi * Ai * ni)
    # CF = F/(qinf*sref) = 1/sref * sum(Cpi * Ai * ni)

    # M = r x F
    # assume Lref = cref
    lref = np.array([bref, cref, bref])
    log.info(f'bref={bref} cref={cref} -> Lref={lref.tolist()}')
    log.info(f'reference_point ({xyz_units}) = {reference_point}')
    # M = (xyz - xyz_ref) x F = dr x F = dr x qinf * sum(Cpi * Ai * ni)
    # M/qinf = dr x sum(Cpi * Ai * ni)
    # CM = M / (q*S*lref) = M / q * 1/(S*lref)
    # CM = dr x sum(Cpi * Ai * ni) * 1/(S*lref) = aero_moment_coeff / (S*lref)
    total_aero_force_coeff = aero_force_coeff_per_q.sum(axis=0) / sref
    total_aero_moment_coeff = aero_moment_coeff_per_q.sum(axis=0) / (sref * lref)
    assert len(total_aero_force_coeff) == 3, total_aero_force_coeff
    log.info(f'total_aero_force_coeff  = F/qS = {total_aero_force_coeff}')
    log.info(f'total_aero_moment_coeff = M/qL = {total_aero_moment_coeff}')

    nelements_structure = len(structure_eids)
    log.info(f'nelements_structure = {nelements_structure}')
    log = structure_model.log
    structure_nodes, structure_xyz = get_structure_xyz(structure_model, xyz_units=xyz_units)
    out = get_mapped_structure(
        structure_model, structure_eids, fdtype=fdtype)
    structure_eids, structure_area, structure_centroid, structure_normal = out

    abs_normals = np.abs(structure_normal)
    imax = np.argmax(abs_normals, axis=1)
    imax_min = imax.min()
    # log.info(f'imax = {imax.tolist()}')
    assert len(imax) == nelements_structure
    assert imax_min > 0, imax

    iy_structure = np.where(imax == 1)[0]
    iz_structure = np.where(imax == 2)[0]
    eids_y = structure_eids[iy_structure]
    eids_z = structure_eids[iz_structure]
    assert len(eids_y) or len(eids_z), (eids_y, eids_z)

    model2 = BDF(log=structure_model.log)
    ny_structure = len(iy_structure)
    nz_structure = len(iz_structure)
    log.info(f'nstructure={len(structure_eids)} = ny_structure={ny_structure} + nz_structure={nz_structure}')
    if nz_structure > 0:
        assert len(aero_centroid) == len(aero_force_coeff_per_q), (len(aero_centroid), len(aero_force_coeff_per_q))
        # aero_centroid_z = aero_centroid[iz_structure, :]
        # aero_Cp_z = aero_Cp[iz]

        structure_eids_z = structure_eids[iz_structure]
        structure_area_z = structure_area[iz_structure]
        structure_box_centroid = structure_centroid[iz_structure, :]
        structure_box_moment_center = np.zeros((nz_structure, 3), dtype=fdtype)

        for i, eid, aero_forcei in zip(count(), eids_z, aero_force_coeff_per_q):
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
            structure_box_moment_center[i, :] = p14 + chord/4.

        # for each aero element, get the closest structure box
        assert len(structure_box_centroid) > 0, structure_box_centroid
        assert len(structure_box_centroid) > 0, aero_xyz_nodal

        structure_box_centroid_panel = structure_box_centroid[:, [0, 1]]
        aero_centroid_panel = aero_centroid[:, [0, 1]]
        # structure_box_moment_center_panel = structure_box_moment_center[:, [0, 1]]

        neidsi = 10
        msg = f'structure_box_centroid[:{neidsi},:]; n={len(structure_box_centroid)}\n'
        for x, y, z in structure_box_centroid[:neidsi, :]:
            msg += f'    [{x:.3f}, {y:.3f}, {z:.3f}]\n'
        log.debug(msg)
        # del structure_box_centroid, structure_box_moment_center
        # del aero_centroid

        tree = _get_tree(structure_box_centroid_panel)
        nnodes_structure = len(structure_box_centroid_panel)
        print(structure_box_centroid_panel.shape, aero_centroid_panel.shape)
        if map_location == 'centroid':
            unused_deq, ieq = tree.query(aero_centroid_panel, k=1)
            slots = np.where(ieq < nnodes_structure)
            iaero = ieq[slots]
        else:
            raise RuntimeError(map_location)

        log.info(f'deq={unused_deq} n={len(unused_deq)} nnodes={nnodes_structure}')
        log.info(f'ieq={ieq} n={len(ieq)} max={ieq.max()}')
        uieq = np.unique(ieq)
        log.info(f'unique ieq={uieq} n={len(uieq)}')

        # irows: aero?
        # icols: structure?
        # iaero = slots[0]
        # print(f'slots = {slots}')
        # print(f'slots[1] = {slots[1]}')
        log.info(f'unique iaero={iaero} n={len(iaero)}')

        istructure = np.arange(len(iz_structure))
        # if map_type == 'pressure':
        #     aero_pressure_centroidal = aero_Cp_centroidal[iaero] * scale
        #     mapped_structure_elements = structure_nodes[istructure]
        #     for eid, pressure in zip(mapped_structure_elements, aero_pressure_centroidal):
        #         model2.add_pload2(pressure_sid, eids=[eid], pressure=pressure)
        if map_type == 'force_moment':
            assert len(aero_area) == len(aero_force_coeff_per_q)
            assert len(structure_eids_z) == len(structure_box_moment_center)
            assert len(structure_box_moment_center) == len(structure_area_z), (len(structure_box_moment_center), len(structure_area_z))

            # n = 100

            eid0 = eids_z[0]
            polygons = []
            for i, eid in enumerate(eids_z):
                if i == 0:
                    continue
                elem0 = structure_model.elements[eid]
                nodes0 = elem0.nodes
                inodes0 = np.searchsorted(structure_nodes, nodes0)
                elem_xyz0 = structure_xyz[inodes0, :]
                polygons.append(elem_xyz0)
            quad = Poly3DCollection(polygons, facecolor='blue', edgecolor='black', linewidth=1, alpha=0.2)
            ax.add_collection3d(quad)

            # the important one...
            polygons2 = []
            elem0 = structure_model.elements[eid0]
            nodes0 = elem0.nodes
            inodes0 = np.searchsorted(structure_nodes, nodes0)
            elem_xyz0 = structure_xyz[inodes0, :]
            polygons2 = [elem_xyz0]
            quad = Poly3DCollection(polygons2, facecolor='red', edgecolor='black', linewidth=1, alpha=0.8)
            ax.add_collection3d(quad)

            #--------------------------

            print(f'iaero = {iaero}')
            show = False
            if 0:
                # plot all aero centroids associated with structure box id=0
                structure_box_id = 0
                iaero0 = np.where(iaero == structure_box_id)[0]
                aero_centroid_plot = aero_centroid[iaero0, :]
                # scatter = ax.scatter(x, y, z, c=color_values, cmap='viridis', s=50)
                ax.scatter3D(aero_centroid_plot[:, 0], aero_centroid_plot[:, 1], aero_centroid_plot[:, 2],
                             c=aero_Cp[iaero0], cmap='viridis')
                show = True

            if 0:
                # plot all aero elements associated with structure box id=0
                structure_box_id = 0
                iaero0 = np.where(iaero == structure_box_id)[0]
                assert len(iaero) == len(aero_tri_nodes), (len(iaero), len(aero_tri_nodes))

                iaero_tri_nodes = np.searchsorted(aero_node_id, aero_tri_nodes[iaero0, :])
                polygon_tris = aero_xyz[iaero_tri_nodes, :]
                print('polygon_tris.shape =', polygon_tris.shape)

                # works - hard to see colors
                # tri = Poly3DCollection(polygon_tris, facecolor='blue', edgecolor='black', linewidth=0, alpha=1.0)

                # works - bit transparent tho
                cmap = cm.get_cmap('viridis')
                scalar_values = aero_Cp[iaero0]
                norm = plt.Normalize(vmin=min(scalar_values), vmax=max(scalar_values))
                facecolor = cmap(norm(scalar_values))
                tri = Poly3DCollection(polygon_tris, linewidth=0, alpha=1.0, facecolor=facecolor)
                ax.add_collection3d(tri)
                show = True

            if 0:
                # plot the whole aero model
                iaero_tri_nodes = np.searchsorted(aero_node_id, aero_tri_nodes)
                polygon_tris = aero_xyz[iaero_tri_nodes, :]
                tri = Poly3DCollection(polygon_tris, facecolor='blue', edgecolor='black', linewidth=1, alpha=0.2)
                ax.add_collection3d(tri)
                show = True

            if show:
                plt.show()

            panel_dim = 'z'
            # assert len(istructure) == len(iaero), (len(istructure), len(iaero))
            structure_eids_z[iaero]
            map_force_moment_centroid_tri(
                model2, panel_dim,
                aero_centroid, aero_force_coeff_per_q,
                structure_eids_z[iaero], structure_box_moment_center[iaero, :], structure_area_z[iaero],
                pressure_sid=pressure_sid,
                qinf=qinf,
            )
        else:  # pragma: no cover
            raise RuntimeError(map_type)

    return model2


def map_pressure_centroid_avg(model2: BDF,
                              aero_area, aero_Cp_centroidal, iaero,
                              structure_nodes, structure_area, istructure,
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
                          qinf: float=1.0) -> None:
    """maps the pressure at the centroid of the panel (not the average)"""
    aero_pressure_centroidal = aero_Cp_centroidal[iaero] * qinf
    mapped_structure_elements = structure_nodes[istructure]
    for eid, pressure in zip(mapped_structure_elements, aero_pressure_centroidal):
        model2.add_pload2(pressure_sid, eids=[eid], pressure=pressure)
    return


def map_force_centroid_tri(model2: BDF,
                           aero_area, aero_Cp_centroidal, iaero,
                           structure_nodes, istructure,
                           pressure_sid: int=1,
                           qinf: float=1.0) -> None:
    aero_pressure_centroidal = aero_Cp_centroidal[iaero] * qinf
    aero_force_centroidal = aero_pressure_centroidal * aero_area
    mapped_structure_nodes = structure_nodes[istructure]
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
                                  panel_dim: str,
                                  aero_centroid, aero_force_per_q_centroidal,
                                  structure_eids: np.ndarray,
                                  structure_xyz: np.ndarray,
                                  structure_area: np.ndarray,
                                  pressure_sid: int=1,
                                  force_sid: int=2,
                                  moment_sid: int=3,
                                  qinf: float=1.0,
                                  fdtype: str='float64') -> None:
    """

    Parameters
    ----------
    model2 : BDF()
        the output model
    panel_dim : str
        y or z
    aero_centroid : (naero, 3) float np.ndarray
        centroids of the aero panels
    aero_force_per_q_centroidal : (naero, 3) float np.ndarray
        F/q for the panel (Cp*area*normal)
    structure_eids : (nstructure,) int np.ndarray
        the element ids
    structure_xyz : (nstructure,3) float np.ndarray
        the centroid / quarter chord locations to sum moments about
    structure_area : (nstructure,) float np.ndarray
        the area of each aero-box to normalize by
    pressure_sid : int; default=1
        the output pressure load id
    force_sid : int; default=2
        the output force load id
    moment_sid : int; default=3
        the output moment load id

    """
    # validation
    panel_dim = panel_dim.lower()
    assert panel_dim in ['y', 'z'], panel_dim
    assert len(structure_eids) == len(structure_xyz), (len(structure_eids), len(structure_xyz))
    assert len(structure_eids) == len(structure_area), (len(structure_eids), len(structure_area))
    assert len(aero_centroid) == len(aero_force_per_q_centroidal), (len(aero_centroid), len(aero_force_per_q_centroidal))

    #---------------------

    aero_force_centroidal = aero_force_per_q_centroidal * qinf
    # aero_moment_centroidal = aero_moment_per_q_centroidal * qinf
    del aero_force_per_q_centroidal, qinf
    # aero_pressure_centroidal = aero_Cp_centroidal[iaero] * qinf
    # aero_force_centroidal = aero_pressure_centroidal * aero_area

    force_temp = {}
    moment_temp = {}
    nforce_total = {}
    print(f'structure_nodes = {structure_eids}')
    assert len(structure_eids) == len(structure_xyz), (len(structure_eids), len(structure_xyz))
    assert len(aero_centroid) == len(aero_force_centroidal), (len(aero_centroid), len(aero_force_centroidal))
    assert len(structure_eids) == len(aero_centroid), (len(structure_eids), len(aero_centroid))

    for eid, scentroid, acentroid, aforce in zip_longest(structure_eids, structure_xyz,
                                                         aero_centroid, aero_force_centroidal):
        if eid not in force_temp:
            force_temp[eid] = np.zeros(3, dtype=fdtype)
            moment_temp[eid] = np.zeros(3, dtype=fdtype)
            nforce_total[eid] = np.zeros(3, dtype=fdtype)

        # print(f'force_temp[{nid}] = {force_temp[nid]}')
        # print(f'aforce = {aforce}')
        force_temp[eid] += aforce

        # print(f'acentroid = {acentroid}')
        # print(f'scentroid = {scentroid}')
        dr = acentroid - scentroid
        # print(f'dr     = {dr}')
        moment_temp[eid] += np.cross(dr, aforce)
        nforce_total[eid] += 1

    # assert len(structure_area) == len(nforce_total), (len(structure_area), len(nforce_total))

    # mapping to a dictionary because not all spots will will be used
    # need to offset it by 1 because indices are 0-based
    structure_area_dict = {eid: area for eid, area in zip(structure_eids, structure_area)}
    mag = 1.0
    for eid, ntotal in nforce_total.items():
        sarea = structure_area_dict[eid]
        force = force_temp[eid]    # / ntotal
        moment = moment_temp[eid]  # / ntotal
        pressure = force / sarea
        if panel_dim == 'z':
            pressurei = pressure[2]
            force[0] = 0.  # fx
            force[1] = 0.  # fy
            moment[0] = 0.  # mx
            moment[2] = 0.  # mz
        elif panel_dim == 'y':
            pressurei = pressure[1]
            force[0] = 0.  # fx
            force[2] = 0.  # fz
            moment[0] = 0.  # mx
            moment[1] = 0.  # my
        else:
            raise RuntimeError(panel_dim)

        model2.add_pload2(pressure_sid, eids=[eid], pressure=pressurei)

        # sid: int, node: int, mag: float, xyz: np.ndarray
        # model2.add_force(force_sid, eid, mag, force)
        # model2.add_moment(moment_sid, eid, mag, moment)
    return


def old_average(n1, n2, n3):
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
    area = np.sqrt(s * (s-a) * (s-b) * (s-c))
    del area


def pressure_map_to_structure_model(aero_model: Cart3D | Tecplot,
                                    structure_model: BDF,
                                    structure_eids: np.ndarray,
                                    structure_sid=0,
                                    map_type: str='pressure',
                                    scale: float=1.0,
                                    sref: float=1.0,
                                    cref: float=1.0,
                                    bref: float=1.0,
                                    qinf: float=1.0,
                                    pressure_sid: int=1,
                                    force_sid: int=2,
                                    moment_sid: int=3,
                                    regions_to_include=None,
                                    regions_to_remove=None,
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

    aero_dict = get_aero_pressure_centroid(
        aero_model, aero_format, map_type,
        regions_to_include=regions_to_include,
        regions_to_remove=regions_to_remove,
    )
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
            qinf=qinf,
        )
    elif map_type == 'force_moment':
        map_force_moment_centroid_tri(
            model2,
            aero_Cp_centroidal, iaero,
            structure_nodes, istructure,
            pressure_sid=pressure_sid,
            force_sid=force_sid,
            moment_sid=moment_sid,
            qinf=qinf,
        )
    elif map_type == 'force':
        map_force_centroid_tri(
            model2,
            aero_area, aero_Cp_centroidal, iaero,
            structure_nodes, istructure,
            pressure_sid=pressure_sid,
            qinf=qinf,
        )
    else:  # pragma: no cover
        raise RuntimeError(map_type)
    return model2


def array3() -> np.ndarray:
    return np.zeros(3, dtype='float64')
