# import os
from itertools import count, zip_longest
from collections import defaultdict
from typing import cast
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.cm as cm

from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.bdf_equivalence import _get_tree

from pyNastran.converters.cart3d.cart3d import Cart3D
#from pyNastran.converters.fluent.fluent import read_fluent
from pyNastran.converters.tecplot.tecplot import Tecplot

from pyNastran.dev.tools.pressure_map_structure_setup import (
    get_structural_eids_from_csv_load_id, get_structure_xyz, get_mapped_structure)
from pyNastran.dev.tools.pressure_map_aero_setup import get_aero_model, get_aero_pressure_centroid


def pressure_filename_to_fa2j(pressure_filename: PathLike,
                              fa2j_filename: PathLike, sid: int=1):
    model = read_bdf(pressure_filename, xref=False)
    loads = model.loads[sid]
    pressures = []
    for load in loads:
        assert load.type == 'PLOAD2', load
        pressures.append(load.pressure)

    reals = np.array(pressures)
    nrows = len(pressures)
    comment = (
        f'FA2J = Cp_q=1\n'
        f'file = {str(pressure_filename)}; sid={sid:d}\n'
        'DMI\tNAME\t0\tFORM\tTIN\tTOUT\t\tnrows\tncol\n'
        f'pressure_range=[{reals.min():g}, {reals.max():g}]')
    model.add_dmi_column('FA2J', tin=1, tout=1, nrows=nrows, reals=reals, comment=comment)
    model.loads = {}
    model.write_bdf(fa2j_filename)


def pressure_filename_to_wkk_diag(cfd_force_moment_bdf_filename1: PathLike | np.ndarray,
                                  cfd_force_moment_bdf_filename2: PathLike | np.ndarray,
                                  nastran_f06_op2_force_moment1: PathLike | np.ndarray,
                                  nastran_f06_op2_force_moment2: PathLike | np.ndarray,
                                  wkk_filename: PathLike,
                                  dtheta_deg: float=1.,
                                  force_sid1: int=2, moment_sid1: int=3,
                                  force_sid2: int=2, moment_sid2: int=3):
    """
    TODO: not done

    # method 1
    dtheta_deg = 2-1 = 1
    case1 = cfd_force_moment1 / nastran_force_moment1
    case2 = cfd_force_moment2 / nastran_force_moment2
    wkk = (case2 - case1) / radians(dtheta_deg)

    # method 2
    CLA_nastran = (nastran_force_moment2 - nastran_force_moment1) / radians(dtheta_deg)
    CLA_cfd     = (cfd_force_moment2 - cfd_force_moment1) / radians(dtheta_deg)
    Wkk = CLA_cfd / CLA_nastran
        = (cfd_force_moment2 - cfd_force_moment1) / (nastran_force_moment2 - nastran_force_moment1)
    """
    force_moment1 = _get_force_moment(cfd_force_moment_bdf_filename1, force_sid1, moment_sid1)
    force_moment2 = _get_force_moment(cfd_force_moment_bdf_filename2, force_sid2, moment_sid2)
    dtheta = np.radians(dtheta_deg)
    force_moment_per_radian = (force_moment2 - force_moment1) / dtheta

    force_per_radian = force_moment_per_radian[:, 0]
    moment_per_radian = force_moment_per_radian[:, 1]
    reals = force_moment_per_radian.flatten()
    nrows = len(reals)
    comment = (
        f'Wkk = (file2 - file1) / dtheta_radians; dtheta_deg={dtheta_deg:g}\n'
        f'cfd_file1 = {str(cfd_force_moment_bdf_filename1)}; force_sid={force_sid1:d} moment_sid={moment_sid1}\n'
        f'cfd_file2 = {str(cfd_force_moment_bdf_filename2)}; force_sid={force_sid2:d} moment_sid={moment_sid2}\n'
        'DMI\tNAME\t0\tFORM\tTIN\tTOUT\t\tnrows\tncol\n'
        f'force_range=[{force_per_radian.min():g}, {force_per_radian.max():g}]\n'
        f'moment_range=[{moment_per_radian.min():g}, {moment_per_radian.max():g}]')

    model = BDF()
    GCj = np.arange(1, nrows+1, dtype='int32')
    GCi = np.arange(1, nrows+1, dtype='int32')
    model.add_dmi('WKK', form='square', tin=1, tout=1,
                  nrows=nrows, ncols=nrows,
                  GCj=GCj, GCi=GCi, Real=reals, comment=comment)
    model.loads = {}
    model.write_bdf(wkk_filename)


def _get_force_moment(pressure_filename: PathLike | np.ndarray,
                      force_sid: int, moment_sid: int) -> np.ndarray:
    """loads the output of pressure_map"""
    if isinstance(pressure_filename, np.ndarray):
        assert pressure_filename.ndim == 2, pressure_filename.shape
        return pressure_filename

    model = read_bdf(pressure_filename, xref=False)
    force_loads = model.loads[force_sid]
    moment_loads = model.loads[moment_sid]
    forces = []
    moments = []
    for load in force_loads:
        assert load.type == 'PLOAD2', load
        forces.append(load.pressure)

    for load in moment_loads:
        assert load.type == 'PLOAD2', load
        moments.append(load.pressure)

    force_moment = np.column_stack([forces, moments])
    return force_moment


def pressure_map(aero_filename: PathLike,
                 nastran_filename: PathLike | BDF,
                 eids_structure=np.array([]),
                 eid_csv_filename: PathLike='',
                 eid_load_id: int=10,
                 aero_format: str='cart3d',
                 map_type: str='pressure',
                 method: str='full_model',
                 xyz_units: str='???',
                 pressure_units: str='???',
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
    map_type: str
        pressure, force, force_moment
    method : str
        full_model, panel_model
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

    +-------------------+----------+-------+---------------+
    | method / map_type | pressure | force | force_moment  |
    +===================+==========+=======+===============+
    | full_model        |    x     |   x   |     TODO      |
    +-------------------+----------+-------+---------------+
    | panel_model       |   NOTE   | NOTE  |      x        |
    +-------------------+----------+-------+---------------+
    *NOTE: panel_model outputs pressure/force/moment as sid=1,2,3. Select what you want.

    TODO: Support vtail pressure box mapping:
      - requires:
        - same aero model
        - different bdfs
        - different aero regions

    TODO: add support for excluding y/z pressure mapping
    TODO: consider aero panel dihedral; need to reduce load by cos(theta); do a dot product
    """
    assert pressure_filename != 'pressure.bdf', pressure_filename
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
    if isinstance(nastran_filename, BDF):
        structure_model = nastran_filename
    else:
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
            pressure_units=pressure_units,
            # idtype=idtype,
            fdtype=fdtype)
    # elif method == 'strip_model':
        # make the strips
        # just run the simple panel_mode code
    elif method == 'full_model':
        pressure_model = pressure_map_to_structure_model(
            aero_model, structure_model, structure_eids,
            reference_point,
            pressure_sid=pressure_sid,
            force_sid=force_sid,
            moment_sid=moment_sid,
            #aero_format=aero_format,
            map_type=map_type,
            qinf=qinf, pressure_units='Pa',
            sref=sref, bref=bref, cref=cref,
            regions_to_include=regions_to_include,
            regions_to_remove=regions_to_remove,
            #idtype=idtype,
            fdtype=fdtype)
    else:
        expected = ['panel_model', 'full_model']
        raise RuntimeError(f'method={method}; expected={list(expected)}')
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
                                pressure_units: str='',
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

    aero_xyz = aero_dict['xyz_nodal']
    aero_node_id = aero_dict['node_id']

    aero_tri_nodes = aero_dict['tri_nodes']
    aero_tri_cp = aero_dict['tri_Cp_centroid']
    aero_tri_area = aero_dict['tri_area']
    aero_tri_centroid = aero_dict['tri_centroid']
    aero_tri_normal = aero_dict['tri_normal']
    assert np.abs(aero_tri_normal).max() <= 1.001, (aero_tri_normal.min(), aero_tri_normal.max())

    aero_quad_nodes = aero_dict['quad_nodes']
    aero_quad_cp = aero_dict['quad_Cp_centroid']
    aero_quad_area = aero_dict['quad_area']
    aero_quad_centroid = aero_dict['quad_centroid']
    aero_quad_normal = aero_dict['quad_normal']
    nquad = len(aero_quad_area)
    if nquad:
        assert np.abs(aero_quad_normal).max() <= 1.001, (aero_quad_normal.min(), aero_quad_normal.max())

    assert len(aero_tri_cp) == len(aero_tri_centroid), (len(aero_tri_cp), len(aero_tri_centroid))
    aero_tri_dr = aero_tri_centroid - reference_point
    aero_quad_dr = aero_quad_centroid - reference_point

    # F/qinf per panel
    aero_tri_force_coeff_per_q = (aero_tri_cp * aero_tri_area)[:, np.newaxis] * aero_tri_normal
    aero_quad_force_coeff_per_q = (aero_quad_cp * aero_quad_area)[:, np.newaxis] * aero_quad_normal

    # M/qinf per panel
    aero_tri_moment_coeff_per_q = np.cross(aero_tri_dr, aero_tri_force_coeff_per_q)
    aero_quad_moment_coeff_per_q = np.cross(aero_quad_dr, aero_quad_force_coeff_per_q)

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
    total_aero_force_coeff = (aero_tri_force_coeff_per_q.sum(axis=0) + aero_quad_force_coeff_per_q.sum(axis=0)) / sref
    total_aero_moment_coeff = (aero_tri_moment_coeff_per_q.sum(axis=0) + aero_quad_moment_coeff_per_q.sum(axis=0)) / (sref * lref)
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

    # find the direction of max |normal| to split:
    #  - z-normal panels (e.g., wing, tail)
    #  - y-normal panels (e.g., rudder)
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

    bdf_model_out = BDF(log=structure_model.log)
    ny_structure = len(iy_structure)
    nz_structure = len(iz_structure)
    log.info(f'nstructure={len(structure_eids)} = ny_structure={ny_structure} + nz_structure={nz_structure}')
    if nz_structure > 0:
        panel_dim = 'z'
        z_structure_sign = np.sign(structure_normal[iz_structure, 2])
        _map_pressure_panel_model(
            structure_model, bdf_model_out,
            map_type, map_location, panel_dim,
            aero_node_id, aero_xyz,
            aero_tri_nodes, aero_tri_area, aero_tri_cp, aero_tri_normal, aero_tri_centroid, aero_tri_force_coeff_per_q,
            aero_quad_nodes, aero_quad_area, aero_quad_cp, aero_quad_normal, aero_quad_centroid, aero_quad_force_coeff_per_q,
            structure_nodes, structure_xyz,
            structure_eids, structure_area, structure_centroid, iz_structure, z_structure_sign,
            qinf=qinf,
            pressure_sid=pressure_sid,
            force_sid=force_sid,
            moment_sid=moment_sid,
            pressure_units=pressure_units,
        )
    if ny_structure > 0:
        panel_dim = 'y'
        y_structure_sign = np.sign(structure_normal[iy_structure, 1])
        _map_pressure_panel_model(
            structure_model, bdf_model_out,
            map_type, map_location, panel_dim,
            aero_node_id, aero_xyz,
            aero_tri_nodes, aero_tri_area, aero_tri_cp, aero_tri_normal, aero_tri_centroid, aero_tri_force_coeff_per_q,
            aero_quad_nodes, aero_quad_area, aero_quad_cp, aero_quad_normal, aero_quad_centroid, aero_quad_force_coeff_per_q,
            structure_nodes, structure_xyz,
            structure_eids, structure_area, structure_centroid, iy_structure, y_structure_sign,
            qinf=qinf,
            pressure_sid=pressure_sid,
            force_sid=force_sid,
            moment_sid=moment_sid,
        )
    return bdf_model_out


def _map_pressure_panel_model(structure_model: BDF,
                              bdf_model_out: BDF,
                              map_type: str, map_location: str, panel_dim: str,
                              # aero nodes
                              aero_node_id, aero_xyz,
                              # aero elements
                              aero_tri_nodes, aero_tri_area, aero_tri_cp, aero_tri_normal, aero_tri_centroid, aero_tri_force_coeff_per_q,
                              aero_quad_nodes, aero_quad_area, aero_quad_cp, aero_quad_normal, aero_quad_centroid, aero_quad_force_coeff_per_q,
                              # structure
                              structure_nodes, structure_xyz,
                              structure_eids, structure_area, structure_centroid, iz_structure, z_structure_sign,
                              qinf: float=1.0, #pressure_units='',
                              pressure_sid: int=1,
                              force_sid: int=2,
                              moment_sid: int=3,
                              pressure_units: str='') -> None:
    eids_z = structure_eids[iz_structure]

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_box_aspect((1, 1, 1))

    nz_structure = len(iz_structure)
    #-------------------------------------------
    log = structure_model.log
    assert len(aero_tri_centroid) == len(aero_tri_force_coeff_per_q), (len(aero_tri_centroid), len(aero_tri_force_coeff_per_q))
    # aero_centroid_z = aero_centroid[iz_structure, :]
    # aero_Cp_z = aero_Cp[iz]

    structure_eids_z = structure_eids[iz_structure]
    structure_area_z = structure_area[iz_structure]
    structure_box_centroid = structure_centroid[iz_structure, :]

    fdtype = structure_area.dtype
    structure_box_moment_center = np.zeros((nz_structure, 3), dtype=fdtype)

    for i, eid, aero_forcei in zip(count(), eids_z, aero_tri_force_coeff_per_q):
        elem = structure_model.elements[eid]
        assert elem.type == 'CQUAD4', elem
        inids = np.searchsorted(structure_nodes, elem.nodes)
        quad_xyz = structure_xyz[inids, ]
        quad_xy = quad_xyz[:, [0, 1]]
        p1 = quad_xyz[0, :]
        p2 = quad_xyz[1, :]
        p3 = quad_xyz[2, :]
        p4 = quad_xyz[3, :]
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
    assert len(aero_xyz) > 0, aero_xyz

    structure_box_centroid_panel = structure_box_centroid[:, [0, 1]]
    aero_tri_centroid_panel = aero_tri_centroid[:, [0, 1]]
    aero_quad_centroid_panel = aero_quad_centroid[:, [0, 1]]
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
    # print(structure_box_centroid_panel.shape, aero_tri_centroid_panel.shape)
    if map_location == 'centroid':
        deq_tri, ieq_tri = tree.query(aero_tri_centroid_panel, k=1)
        deq_quad, ieq_quad = tree.query(aero_quad_centroid_panel, k=1)
        del deq_tri, deq_quad
        slots_tri = np.where(ieq_tri < nnodes_structure)
        slots_quad = np.where(ieq_quad < nnodes_structure)
        iaero_tri = ieq_tri[slots_tri]
        iaero_quad = ieq_quad[slots_quad]
    else:
        raise RuntimeError(map_location)

    # log.info(f'deq={unused_deq} n={len(unused_deq)} nnodes={nnodes_structure}')
    # log.info(f'ieq={ieq} n={len(ieq)} max={ieq.max()}')
    # uieq = np.unique(ieq)
    # log.info(f'unique ieq={uieq} n={len(uieq)}')

    # irows: aero?
    # icols: structure?
    # iaero = slots[0]
    # print(f'slots = {slots}')
    # print(f'slots[1] = {slots[1]}')
    log.info(f'unique iaero_tri={iaero_tri} n={len(iaero_tri)}')
    log.info(f'unique iaero_quad={iaero_quad} n={len(iaero_quad)}')

    # istructure = np.arange(len(iz_structure))
    # if map_type == 'pressure':
    #     aero_pressure_centroidal = aero_Cp_centroidal[iaero] * scale
    #     mapped_structure_elements = structure_nodes[istructure]
    #     for eid, pressure in zip(mapped_structure_elements, aero_pressure_centroidal):
    #         bdf_model_out.add_pload2(pressure_sid, eids=[eid], pressure=pressure)
    if map_type == 'force_moment':
        assert len(aero_tri_area) == len(aero_tri_force_coeff_per_q)
        assert len(structure_eids_z) == len(structure_box_moment_center)
        assert len(structure_box_moment_center) == len(structure_area_z), (len(structure_box_moment_center), len(structure_area_z))

        # plot all aero centroids associated with this structure box id
        #
        # if eid_filter > 0, the single eid Will be
        # plotted so you can debug the loads
        eid_filter = 0
        # eid_filter = 4205006
        #eid_filter = 4_212_002
        # eid_filter = 1_201_035
        istructure_box_id = np.where(eids_z == eid_filter)[0]
        plot_special_quad = (eid_filter > 0)
        plot_base_quads = plot_special_quad
        plot_scatter = False
        plot_aero_model = False
        if plot_base_quads:  # pragma: no cover
            polygons = []
            for i, eid in enumerate(eids_z):
                if i == istructure_box_id:
                    continue
                elem0 = structure_model.elements[eid]
                nodes0 = elem0.nodes
                inodes0 = np.searchsorted(structure_nodes, nodes0)
                elem_xyz0 = structure_xyz[inodes0, :]
                polygons.append(elem_xyz0)
            quad = Poly3DCollection(polygons, facecolor='blue', edgecolor='black', linewidth=1, alpha=0.2)
            ax.add_collection3d(quad)

        if plot_special_quad:  # pragma: no cover
            # make this quad red
            elem0 = structure_model.elements[eid_filter]
            nodes0 = elem0.nodes
            inodes0 = np.searchsorted(structure_nodes, nodes0)
            elem_xyz0 = structure_xyz[inodes0, :]
            polygons2 = [elem_xyz0]
            quad = Poly3DCollection(polygons2, facecolor='red', edgecolor='black', linewidth=1, alpha=0.8)
            ax.add_collection3d(quad)

        #--------------------------

        log.debug(f'iaero_tri = {iaero_tri}')
        show = False

        # if plot_scatter:  # pragma: no cover
        #     # only the first box
        #     # istructure_box_id = 0
        #     iaero0 = np.where(iaero_tri == istructure_box_id)[0]
        #     aero_centroid_plot = aero_tri_centroid[iaero0, :]
        #     # scatter = ax.scatter(x, y, z, c=color_values, cmap='viridis', s=50)
        #     ax.scatter3D(aero_centroid_plot[:, 0], aero_centroid_plot[:, 1], aero_centroid_plot[:, 2],
        #                  c=aero_cp[iaero0], cmap='viridis')
        #     show = True

        if plot_special_quad:  # pragma: no cover
            # plot all aero elements associated with structure box id=0
            iaero0 = np.where(iaero_tri == istructure_box_id)[0]
            assert len(iaero_tri) == len(aero_tri_nodes), (len(iaero_tri), len(aero_tri_nodes))

            iaero_tri_nodes = np.searchsorted(aero_node_id, aero_tri_nodes[iaero0, :])
            polygon_tris = aero_xyz[iaero_tri_nodes, :]

            # works - hard to see colors b/c
            # tri = Poly3DCollection(polygon_tris, facecolor='blue', edgecolor='black', linewidth=0, alpha=1.0)

            # works - bit transparent tho
            cmap = cm.get_cmap('viridis')
            scalar_values = aero_cp[iaero0]
            norm = plt.Normalize(vmin=scalar_values.min(), vmax=scalar_values.max())
            facecolor = cmap(norm(scalar_values))
            tri = Poly3DCollection(polygon_tris, linewidth=0, alpha=1.0, facecolor=facecolor)
            ax.add_collection3d(tri)

            elem0 = structure_model.elements[eid_filter]
            x, y, z = elem0.Centroid()
            areai = elem0.Area()
            delta = scalar_values.max() - scalar_values.min()
            msg = (f'eid={eid_filter}, xyz_centroid=[{x:.3f}, {y:.3f}, {z:.3f}], area_box={areai:.3f}\n'
                   f'Cp: min={scalar_values.min():.3f}; max={scalar_values.max():.3f}; delta={delta:.3f}')
            fig.suptitle(msg)
            show = True

        if plot_aero_model:  # pragma: no cover
            # plot the whole aero model
            if len(aero_tri_nodes):
                iaero_tri_nodes = np.searchsorted(aero_node_id, aero_tri_nodes)
                polygon_tris = aero_xyz[iaero_tri_nodes, :]
                tri = Poly3DCollection(polygon_tris, facecolor='blue', edgecolor='black', linewidth=1, alpha=0.2)
                ax.add_collection3d(tri)

            if len(aero_quad_nodes):
                iaero_quad_nodes = np.searchsorted(aero_node_id, aero_quad_nodes)
                polygon_quads = aero_xyz[iaero_quad_nodes, :]
                quad = Poly3DCollection(polygon_quads, facecolor='blue', edgecolor='black', linewidth=1, alpha=0.2)
                ax.add_collection3d(quad)
            show = True

        if show:  # pragma: no cover
            # Add a color bar which maps values to colors.
            # fig.colorbar(surf, shrink=0.5, aspect=5)
            plt.show()

        # assert len(istructure) == len(iaero), (len(istructure), len(iaero))
        # structure_eids_z[iaero_tri]
        map_panel_force_moment_centroid(
            bdf_model_out, panel_dim,
            aero_tri_area, aero_tri_cp, aero_tri_normal, aero_tri_centroid, aero_tri_force_coeff_per_q,
            aero_quad_area, aero_quad_cp, aero_quad_normal, aero_quad_centroid, aero_quad_force_coeff_per_q,
            structure_eids_z, structure_box_moment_center, structure_area_z, z_structure_sign,
            iaero_tri, iaero_quad,
            qinf,
            pressure_sid=pressure_sid,
            force_sid=force_sid,
            moment_sid=moment_sid,
            pressure_units=pressure_units,
        )
    else:  # pragma: no cover
        raise RuntimeError(map_type)


def map_pressure_centroid_avg(bdf_model_out: BDF,
                              aero_area, aero_cp_centroidal, iaero,
                              structure_nodes, structure_area, istructure,
                              qinf: float,
                              pressure_units: str='',
                              pressure_sid: int=1,
                              scale: float=1.0) -> None:
    """maps the pressure at the centroid of the panel (not the average)"""
    aero_pressure_centroid = aero_cp_centroidal[iaero] * scale
    mapped_structure_elements = structure_nodes[istructure]
    mapped_structure_area = structure_area[istructure]

    comment = f'Pressure; qinf={qinf} {pressure_units} (map_pressure_centroid_avg)'
    for eid, pressure in zip(mapped_structure_elements, aero_pressure_centroid):
        bdf_model_out.add_pload2(pressure_sid, eids=[eid], pressure=pressure, comment=comment)
        comment = ''
    return


def map_pressure_centroid(bdf_model_out: BDF,
                          aero_pressure_centroid,
                          structure_eids,
                          qinf: float,
                          pressure_units: str,
                          pressure_sid: int=1) -> None:
    """maps the pressure at the centroid of the panel (not the average)"""
    assert len(structure_eids) == len(aero_pressure_centroid), (len(structure_eids), len(aero_pressure_centroid))
    comment = f'Pressure; qinf={qinf} {pressure_units} (map_pressure_centroid)'
    for eid, pressure in zip(structure_eids, aero_pressure_centroid):
        bdf_model_out.add_pload2(pressure_sid, eids=[eid], pressure=pressure, comment=comment)
        comment = ''
    return


def map_force_centroid_tri(bdf_model_out: BDF,
                           aero_area, aero_cp_centroid, iaero,
                           structure_nodes, istructure,
                           pressure_sid: int=1,
                           force_sid: int=2,
                           moment_sid: int=3,
                           qinf: float=1.0) -> None:
    aero_pressure_centroid = aero_cp_centroid[iaero] * qinf
    aero_force_centroid = aero_pressure_centroid * aero_area[iaero]
    mapped_structure_nodes = structure_nodes[istructure]

    forces_temp = defaultdict(array3)
    nforce_total = defaultdict(int)
    for nid, force in zip(mapped_structure_nodes, aero_force_centroid):
        forces_temp[nid] += force
        nforce_total[nid] += 1

    mag = 1.0
    for nid, force in forces_temp.items():
        nforce_totali = nforce_total[nid]
        force2 = force / nforce_totali
        bdf_model_out.add_force(pressure_sid, nid, mag, force2, cid=0)
    return


def map_panel_force_moment_centroid(
        bdf_model_out: BDF,
        panel_dim: str,
        aero_tri_area, aero_tri_cp, aero_tri_normal, aero_tri_centroid, aero_tri_force_per_q_centroid,
        aero_quad_area, aero_quad_cp, aero_quad_normal, aero_quad_centroid, aero_quad_force_per_q_centroid,
        structure_eids: np.ndarray,
        structure_xyz: np.ndarray,
        structure_area: np.ndarray,
        structure_z_sign: np.ndarray,
        iaero_tri: np.ndarray,
        iaero_quad: np.ndarray,
        qinf: float, pressure_units: str='',
        pressure_sid: int=1,
        force_sid: int=2,
        moment_sid: int=3,

    ) -> None:
    """
    Maps forces and moments onto the moment reference point,
    which for a panel is the 1/4 chord. We keep quads & tris split
    to help with debugging.

    Parameters
    ----------
    bdf_model_out : BDF()
        the output model
    panel_dim : str
        y or z
    aero_tri_centroid : (naero, 3) float np.ndarray
        centroids of the aero panels
    aero_tri_force_per_q_centroid : (naero, 3) float np.ndarray
        F/q for the panel (Cp*area*normal)
    aero_quad_centroid : (naero, 3) float np.ndarray
        centroids of the aero panels
    aero_quad_force_per_q_centroid : (naero, 3) float np.ndarray
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
    assert isinstance(qinf, float), qinf
    structure_tri_eids = structure_eids[iaero_tri]
    structure_tri_xyz = structure_xyz[iaero_tri]
    structure_tri_area = structure_area[iaero_tri]
    structure_tri_z_sign = structure_z_sign[iaero_tri]

    structure_quad_eids = structure_eids[iaero_quad]
    structure_quad_xyz = structure_xyz[iaero_quad]
    structure_quad_area = structure_area[iaero_quad]
    structure_quad_z_sign = structure_z_sign[iaero_quad]

    # validation
    panel_dim = panel_dim.lower()
    assert panel_dim in ['y', 'z'], panel_dim
    ntri_structure = len(structure_tri_eids)
    nquad_structure = len(structure_quad_eids)
    ntri_aero = len(aero_tri_centroid)
    # nquad_aero = len(aero_quad_centroid)
    assert len(structure_tri_xyz) == ntri_structure, (len(structure_tri_xyz), ntri_structure)
    assert len(structure_tri_area) == ntri_structure, (len(structure_tri_area), ntri_structure)
    assert len(structure_tri_z_sign) == ntri_structure, (structure_tri_z_sign, ntri_structure)

    assert len(structure_quad_xyz) == nquad_structure, (len(structure_quad_xyz), nquad_structure)
    assert len(structure_quad_area) == nquad_structure, (len(structure_quad_area), nquad_structure)
    assert len(structure_quad_z_sign) == nquad_structure, (structure_quad_z_sign, nquad_structure)

    assert len(aero_tri_cp) == ntri_aero, (len(aero_tri_cp), ntri_aero)
    assert len(aero_tri_force_per_q_centroid) == ntri_aero, (len(aero_tri_force_per_q_centroid), ntri_aero)
    assert len(aero_tri_area) == ntri_aero, (len(aero_tri_area), ntri_aero)
    #---------------------

    aero_tri_force_centroid = aero_tri_force_per_q_centroid * qinf
    aero_quad_force_centroid = aero_quad_force_per_q_centroid * qinf
    # aero_moment_centroidal = aero_moment_per_q_centroidal * qinf
    del aero_tri_force_per_q_centroid
    # aero_pressure_centroid = aero_Cp_centroid[iaero] * qinf
    # aero_force_centroid = aero_pressure_centroid * aero_area

    force_tri = {}
    moment_tri = {}
    # print(f'structure_nodes = {structure_tri_eids}')
    assert len(structure_tri_eids) == len(structure_tri_xyz), (len(structure_tri_eids), len(structure_tri_xyz))
    assert len(aero_tri_centroid) == len(aero_tri_force_centroid), (len(aero_tri_centroid), len(aero_tri_force_centroid))
    assert len(structure_tri_eids) == len(aero_tri_centroid), (len(structure_tri_eids), len(aero_tri_centroid))

    eid_filter = 0
    # eid_filter = 4_212_002
    # eid_filter = 1_201_035
    is_eid_filter = (eid_filter > 0)

    structure_eids = np.hstack([structure_tri_eids, structure_quad_eids])
    structure_area = np.hstack([structure_tri_area, structure_quad_area])
    isort = np.argsort(structure_eids)
    structure_eids = structure_eids[isort]
    structure_area = structure_area[isort]
    fdtype = structure_area.dtype
    for eid in np.unique(structure_eids):
        if is_eid_filter and eid != eid_filter:
            continue
        if eid not in force_tri:
            force_tri[eid] = np.zeros(3, dtype=fdtype)
            moment_tri[eid] = np.zeros(3, dtype=fdtype)

    aero_centroid = np.vstack([aero_tri_centroid, aero_quad_centroid])
    stucture_xyz = np.vstack([structure_tri_xyz, structure_quad_xyz])
    aero_force_centroid = np.vstack([aero_tri_force_centroid, aero_quad_force_centroid])
    drs = aero_centroid - stucture_xyz
    aero_moment = np.cross(drs, aero_force_centroid, axis=1)

    aero_force_centroid = aero_moment[isort, :]
    aero_moment = aero_moment[isort, :]
    drs_tri = aero_tri_centroid - structure_tri_xyz
    # drs_quad = aero_tri_centroid - structure_quad_xyz
    aero_tri_moment = np.cross(drs_tri, aero_tri_force_centroid, axis=1)
    # aero_quad_moment = np.cross(drs_quad, aero_quad_force_centroid, axis=1)
    # assert aero_tri_moment.shape == drs_tri.shape, (aero_tri_moment.shape, drs_tri.shape)
    # assert aero_quad_moment.shape == drs_quad.shape, (aero_quad_moment.shape, drs_quad.shape)

    if is_eid_filter:  # pragma: no cover
        # TODO: doesn't show quad loads
        print(f'eid, area, Cp, normal, aforce')

        for eid, scentroid, aarea, acp, anormal, acentroid, aforce, amoment in zip_longest(
                structure_tri_eids, structure_tri_xyz,
                aero_tri_area, aero_tri_cp, aero_tri_normal,
                aero_tri_centroid, aero_tri_force_centroid, aero_tri_moment):
            if eid != eid_filter:
                continue

            # print(f'force_temp[{nid}] = {force_temp[nid]}')
            print(f'{eid}, {aarea:.3f}, {acp:.3f}, {anormal}, {aforce}')
            force_tri[eid] += aforce

            # dr = acentroid - scentroid
            # print(f'acentroid = {acentroid}')
            # print(f'scentroid = {scentroid}')
            # print(f'dr     = {dr}')
            # amoment_expected = np.cross(dr, aforce)
            # assert np.allclose(amoment, amoment_expected), (amoment, amoment_expected)
            moment_tri[eid] += amoment
    else:
        for eid, aforce, amoment in zip_longest(
            structure_eids, aero_force_centroid, aero_moment):
            force_tri[eid] += aforce
            moment_tri[eid] += amoment

    # mapping to a dictionary because not all spots will will be used
    # need to offset it by 1 because indices are 0-based
    structure_area_dict = {eid: area for eid, area in zip(structure_eids, structure_area)}

    comment_p = f'Pressure; qinf={qinf} {pressure_units} (map_panel_force_moment_centroid)'
    comment_f = f'Force; qinf={qinf} {pressure_units} (map_panel_force_moment_centroid)'
    comment_m = f'Moment; qinf={qinf} {pressure_units} (map_panel_force_moment_centroid)'
    for seid, force in force_tri.items():
        sarea = structure_area_dict[seid]
        force = force_tri[seid]
        moment = moment_tri[seid]

        pressure = force / sarea
        if panel_dim == 'z':
            pressurei = pressure[2]
            # force[0] = 0.  # fx
            # force[1] = 0.  # fy
            # moment[0] = 0.  # mx
            # moment[2] = 0.  # mz
            forcei = force[2]
            momenti = moment[1]
        elif panel_dim == 'y':
            pressurei = pressure[1]
            forcei = force[1]
            # force[0] = 0.  # fx
            # force[2] = 0.  # fz
            # moment[0] = 0.  # mx
            # moment[1] = 0.  # my
            momenti = moment[2]
        else:
            raise RuntimeError(panel_dim)

        if is_eid_filter:
            print(f'seid={seid:d} force={force} sarea={sarea:.3f}\n'
                  f'  pressure={pressure} pressurei={pressurei}')
        bdf_model_out.add_pload2(pressure_sid, eids=[seid], pressure=pressurei, comment=comment_p)
        bdf_model_out.add_pload2(force_sid, eids=[seid], pressure=forcei, comment=comment_f)
        bdf_model_out.add_pload2(moment_sid, eids=[seid], pressure=momenti, comment=comment_m)
        comment_p = ''
        comment_f = ''
        comment_m = ''
    return


# def old_average(n1, n2, n3):
#     #       3
#     #      /|\
#     #     / | \
#     #    /  c  \
#     #   / /   \ \
#     #  1---------2
#     #
#     # s: semiperimeter
#     a = np.linalg.norm(n2 - n1)
#     b = np.linalg.norm(n3 - n2)
#     c = np.linalg.norm(n1 - n3)
#     s = 0.5 * (a + b + c)
#     area = np.sqrt(s * (s-a) * (s-b) * (s-c))
#     del area


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
                                    pressure_units='',
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
    aero_centroid = aero_dict['centroid']
    aero_cp_centroid = aero_dict['Cp_centroid']
    aero_area = aero_dict['area']
    naero_elem = len(aero_area)
    naero_node = len(aero_xyz_nodal)

    structure_nodes, structure_xyz = get_structure_xyz(structure_model)
    nstructure_node = len(structure_nodes)
    out = get_mapped_structure(
        structure_model, structure_eids, fdtype=fdtype)
    structure_eids, structure_areas, structure_centroids, structure_normals = out
    nstructure_elem = len(structure_eids)
    assert len(structure_centroids) == len(structure_eids)

    if map_location == 'centroid':
        tree = _get_tree(aero_centroid)
        unused_deq, ieq = tree.query(structure_centroids, k=1)
        slots = np.where(ieq < naero_elem)
        iaero_elem = ieq[slots]
        istructure = None
        # tree = _get_tree(structure_centroids)
        # unused_deq, ieq = tree.query(aero_centroid, k=1)
        # slots = np.where(ieq < nstructure_elem)
        # iaero_elem = None
        # istructure = ieq[slots]

    elif map_location == 'node':
        tree = _get_tree(aero_xyz_nodal)
        unused_deq, ieq = tree.query(structure_xyz, k=4)
        slots = np.where(ieq[:, :] < naero_node)

        # irows: aero?
        # icols: structure?
        irows, icols = slots
        iaero = irows
        istructure = icols
    else:  # pragma: no cover
        raise RuntimeError(map_location)

    #mapped_structure_elements = structure_nodes[istructure]

    bdf_model_out = BDF(log=structure_model.log)
    if map_type == 'pressure':
        aero_pressure_centroid = aero_cp_centroid[iaero_elem] * qinf
        # mapped_structure_elements = structure_nodes[istructure]
        assert nstructure_elem == len(aero_pressure_centroid), (nstructure_elem, len(aero_pressure_centroid))
        map_pressure_centroid(
            bdf_model_out,
            aero_pressure_centroid,
            structure_eids,
            qinf,
            pressure_units,
            pressure_sid=pressure_sid,
        )
    # elif map_type == 'force_moment':
    #     map_force_moment_centroid(
    #         bdf_model_out,
    #         aero_Cp_centroidal, iaero,
    #         structure_nodes, istructure,
    #         pressure_sid=pressure_sid,
    #         force_sid=force_sid,
    #         moment_sid=moment_sid,
    #         qinf=qinf,
    #     )
    elif map_type == 'force':
        map_force_centroid_tri(
            bdf_model_out,
            aero_area, aero_cp_centroid, iaero,
            structure_nodes, istructure,
            pressure_sid=pressure_sid,
            force_sid=force_sid,
            moment_sid=moment_sid,
            qinf=qinf,
        )
    else:  # pragma: no cover
        raise RuntimeError(map_type)
    return bdf_model_out


def array3() -> np.ndarray:
    return np.zeros(3, dtype='float64')
