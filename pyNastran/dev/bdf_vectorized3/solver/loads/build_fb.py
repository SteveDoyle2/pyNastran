import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF
from .pload1 import apply_pload1
from .grav_rforce import apply_rforce, apply_grav

from ..elements.beam import thermal_load_beam
from ..elements.shells import (
    build_pload4_cquad4, build_pload4_ctria3,
    build_thermal_load_cquad4, build_thermal_load_ctria3)

DOF_MAP = dict[tuple[int, int], int]


def build_Fb_from_loadid(model: BDF,
                         dof_map: DOF_MAP,
                         ndof: int,
                         xg: np.ndarray,
                         sset_b: np.ndarray,
                         load_id: int=0,
                         temp_load_id: int=0,
                         fdtype: str='float32'):
    Fb = np.zeros(ndof, dtype=fdtype)
    assert len(xg) == ndof

    if load_id == 0 and temp_load_id == 0:
        return Fb
    log = model.log

    if load_id:
        reduced_loads = model.get_reduced_static_load()

        loads = reduced_loads[load_id]
        for scale, load in loads:
            if load.type == "SLOAD":
                for mag, nid in zip(load.mags, load.nodes):
                    i = dof_map[(nid, 0)]  # TODO: wrong...
                    Fb[i] = mag * scale
            elif load.type == "FORCE":
                # TODO: wrong because it doesn't handle SPOINTs
                fxyz_myz = load.sum_forces_moments()
                fxyz = fxyz_myz[:, :3]
                nids = load.node_id
                log.debug(f"  FORCE nids={nids} Fxyz={fxyz}")
                for fxyzi, nid in zip(fxyz, nids):
                    assert len(fxyzi) == 3
                    fi = dof_map[(nid, 1)]
                    Fb[fi : fi + 3] = fxyzi

            elif load.type == "MOMENT":
                fxyz_myz = load.sum_forces_moments()
                mxyz = fxyz_myz[:, 3:]
                nids = load.node_id
                log.debug(f"  MOMENT nid={nids} Mxyz={mxyz}")
                for mxyzi, nid in zip(mxyz, nids):
                    fi = dof_map[(nid, 4)]
                    assert len(mxyzi) == 3
                    Fb[fi : fi + 3] = mxyzi
            elif load.type == "SPCD":
                for nid, components, enforced in zip(
                    load.nodes, load.components, load.enforced):
                    for component in str(components):
                        dof = int(component)
                        fi = dof_map[(nid, dof)]
                        xg[fi] = enforced
                        Fb[fi] = np.nan
                        sset_b[fi] = True
            elif load.type == "PLOAD1":
                apply_pload1(model, load, scale, Fb, dof_map, log)
            elif load.type == "PLOAD4":
                build_pload4_cquad4(model, Fb, dof_map, load_id)
                build_pload4_ctria3(model, Fb, dof_map, load_id)
            elif load.type == "GRAV":
                apply_grav(model, load, scale, Fb, dof_map, ndof, log)
            elif load.type == "RFORCE":
                apply_rforce(model, load, scale, Fb, dof_map, ndof, log)
            else:
                print(load.get_stats())
                raise NotImplementedError(load)

    # Thermal loads via TEMPERATURE(LOAD) or TEMPERATURE(BOTH)
    if temp_load_id:
        node_temperatures = _get_node_temperatures(model, temp_load_id)
        if node_temperatures:
            log.debug(f"  Thermal load: {len(node_temperatures)} nodes with dT")
            build_thermal_load_cquad4(model, Fb, dof_map, node_temperatures)
            build_thermal_load_ctria3(model, Fb, dof_map, node_temperatures)
            build_thermal_load_beam(model, Fb, dof_map, node_temperatures)
    return Fb


def build_thermal_load_beam(
    model: BDF,
    Fb: np.ndarray,
    dof_map: DOF_MAP,
    node_temperatures: dict[int, float]) -> None:
    """Add beam thermal loads to the global force vector.

    Parameters
    ----------
    model : BDF
        The model.
    Fb : np.ndarray
        Global force vector to accumulate into.
    dof_map : DOF_MAP
        DOF map.
    node_temperatures : dict[int, float]
        Node temperatures.
    """
    for elem in [model.cbar, model.cbeam]:
        if elem.n == 0:
            continue

        area = elem.area()
        xyz1, xyz2 = elem.get_xyz()
        lengths = np.linalg.norm(xyz2 - xyz1, axis=1)
        v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
        e_g_nus = elem.e_g_nu()

        # Get thermal expansion coefficient from MAT1
        mat = model.mat1
        pids = elem.property_id
        for prop in elem.allowed_properties:
            i_lookup, i_all = searchsorted_filter(prop.property_id, pids, msg='')
            if len(i_lookup) == 0:
                continue
            mat_ids = prop.material_id[i_all]
            mat_slice = mat.slice_card_by_material_id(mat_ids)
            alphas = mat_slice.a if hasattr(mat_slice, 'a') else np.zeros(len(mat_ids))
            for idx, ielem in enumerate(i_lookup):
                nid1, nid2 = elem.nodes[ielem]
                T1 = node_temperatures.get(nid1, 0.0)
                T2 = node_temperatures.get(nid2, 0.0)
                dT = 0.5 * (T1 + T2)
                if abs(dT) < 1e-30:
                    continue
                alpha_i = alphas[idx] if alphas[idx] != 0.0 else 0.0
                if alpha_i == 0.0:
                    continue
                e_i = e_g_nus[ielem, 0]
                PG = thermal_load_beam(
                    area[ielem], e_i, alpha_i, lengths[ielem],
                    ihat[ielem], yhat[ielem], zhat[ielem], dT,
                )
                gi1 = dof_map[(nid1, 1)]
                gi2 = dof_map[(nid2, 1)]
                Fb[gi1:gi1 + 6] += PG[:6]
                Fb[gi2:gi2 + 6] += PG[6:]



def _get_node_temperatures(model: BDF,
                           temp_load_id: int) -> dict[int, float]:
    """Build a dict of {node_id: temperature} from TEMP/TEMPD cards.

    Parameters
    ----------
    temp_load_id : int
        Load set ID referencing TEMP/TEMPD cards.

    Returns
    -------
    node_temperatures : dict[int, float]
        Temperature at each grid point.
    """
    node_temperatures: dict[int, float] = {}

    # Default temperature from TEMPD
    default_temp = None
    if model.tempd.n > 0:
        tempd = model.tempd
        idx = np.where(tempd.load_id == temp_load_id)[0]
        if len(idx) > 0:
            default_temp = tempd.temperature[idx[0]]

    # Apply default to all grid points
    if default_temp is not None:
        for nid in model.grid.node_id:
            node_temperatures[nid] = default_temp

    # Override with explicit TEMP cards
    if model.temp.n > 0:
        temp = model.temp
        idx = np.where(temp.load_id == temp_load_id)[0]
        for i in idx:
            inode = temp.inode
            i0, i1 = inode[i]
            for j in range(i0, i1):
                nid = temp.node_id[j]
                t_val = temp.temperature[j]
                node_temperatures[nid] = t_val

    return node_temperatures
