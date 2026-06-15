"""Shell element stress/strain/force/strain_energy recovery for SOL 101.

Computes results at top (z1) and bottom (z2) fibers for CQUAD4 and CTRIA3.
Supports offsets (ZOFFS) which shift the neutral plane relative to grid points.

Output format matches NX Nastran: [fiber_dist, σxx, σyy, τxy, angle, σ_major, σ_minor, σ_vm]
per fiber (top/bottom) per element.
"""
from __future__ import annotations
from itertools import count
from typing import TextIO, TYPE_CHECKING

import numpy as np

from pyNastran.f06.errors import FatalError
from pyNastran.dev.bdf_vectorized3.solver.elements.shells import (
    _dshape_quad4,
    _jacobian,
    _jacobian_inv,
    _membrane_B,
    _bending_B,
    _shape_quad4,
    _get_ABD_for_element,
    _GAUSS_2x2_PTS,
    _GAUSS_2x2_WTS,
    macn2_stiffness,
    _apply_shell_offset,
)


if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    DOF_MAP = dict[tuple[int, int], int]


def _von_mises(sxx: float, syy: float, sxy: float) -> float:
    """Von Mises stress for plane stress."""
    return np.sqrt(sxx**2 - sxx * syy + syy**2 + 3.0 * sxy**2)


def _principal_stresses(sxx: float, syy: float, sxy: float) -> tuple[float, float, float]:
    """Returns (angle_deg, s_major, s_minor) for plane stress."""
    s_avg = 0.5 * (sxx + syy)
    radius = np.sqrt((0.5 * (sxx - syy)) ** 2 + sxy**2)
    s_major = s_avg + radius
    s_minor = s_avg - radius
    if abs(sxx - syy) > 1e-30:
        angle = 0.5 * np.degrees(np.arctan2(2.0 * sxy, sxx - syy))
    else:
        angle = 45.0 if abs(sxy) > 1e-30 else 0.0
    return angle, s_major, s_minor


def recover_shell_stress_cquad4(
    model: BDF,
    dof_map: DOF_MAP,
    xb: np.ndarray) -> dict[int, np.ndarray]:
    """Recover CQUAD4 stresses at centroid, top and bottom fibers.

    Returns
    -------
    results : dict[eid, (2, 8) ndarray]
        Row 0 = bottom (z1), Row 1 = top (z2).
        Columns: [fiber_dist, σxx, σyy, τxy, angle, σ_major, σ_minor, σ_vm]
    """
    cquad4 = model.cquad4
    if cquad4.n == 0:
        return {}

    pshell = model.pshell
    pcomp = model.pcomp
    pcompg = model.pcompg

    grid = model.grid
    nodes = cquad4.nodes
    inid = grid.index(nodes)
    xyz = grid.xyz_cid0()

    p1 = xyz[inid[:, 0], :]
    p2 = xyz[inid[:, 1], :]
    p3 = xyz[inid[:, 2], :]
    p4 = xyz[inid[:, 3], :]

    v13 = p1 - p3
    v24 = p2 - p4
    normal = np.cross(v13, v24)
    ni = np.linalg.norm(normal, axis=1)
    normal /= ni[:, np.newaxis]

    ihat = p2 - p1
    inorm = np.linalg.norm(ihat, axis=1)
    ihat /= inorm[:, np.newaxis]
    jhat = np.cross(normal, ihat, axis=1)
    jhat /= np.linalg.norm(jhat, axis=1)[:, np.newaxis]

    results = {}

    eids = cquad4.element_id
    pids = cquad4.property_id
    zoffsets = cquad4.zoffset

    thetas = get_theta_from_theta_mcid(model, cquad4, normal)

    zoffsets[np.isnan(zoffsets)] = 0.0
    common_pshell_pids = np.intersect1d(pshell.property_id, pids)
    is_pshells = np.array([pid in common_pshell_pids for pid in pids])

    for (i_elem, eid, pid, theta_mat, zoffset, is_pshell) in zip(
            count(), eids, pids, thetas, zoffsets, is_pshells):
        Ti = np.vstack([ihat[i_elem], jhat[i_elem], normal[i_elem]])

        xy = np.zeros((4, 3))
        xy[0] = Ti @ p1[i_elem]
        xy[1] = Ti @ p2[i_elem]
        xy[2] = Ti @ p3[i_elem]
        xy[3] = Ti @ p4[i_elem]
        x_local = xy[:, 0]
        y_local = xy[:, 1]

        A_mat, B_mat, D_mat, Ds_mat, thickness = _get_ABD_for_element(
            model, pid, pshell, pcomp, pcompg)

        # Get z1, z2 (fiber distances from midplane)
        if is_pshell:
            iprop = np.searchsorted(pshell.property_id, pid)
            z1, z2 = pshell.z[iprop, :]
        else:
            z1 = -thickness / 2.0
            z2 = thickness / 2.0

        # _get_ABD_for_element returns A in the MATERIAL frame.
        # For stress recovery in element frame, rotate A by theta/MCID
        # (same rotation used in stiffness assembly).
        #
        # Rotate A from material to element frame
        if abs(theta_mat) > 1e-10:
            c = np.cos(theta_mat)
            s = np.sin(theta_mat)
            T_inv = np.array([
                [c**2, s**2, -2 * s * c],
                [s**2, c**2, 2 * s * c],
                [s * c, -s * c, c**2 - s**2],
            ])
            A_mat = T_inv @ A_mat @ T_inv.T

        # Extract element displacements in local frame
        elem_nodes = nodes[i_elem]
        u_local_20 = np.zeros(20)
        for inode in range(4):
            nid = elem_nodes[inode]
            i_dof = dof_map[(nid, 1)]
            u_g = xb[i_dof : i_dof + 6]
            u_trans = Ti @ u_g[:3]
            u_rot = Ti[0:2, :] @ u_g[3:6]
            u_local_20[5 * inode] = u_trans[0]
            u_local_20[5 * inode + 1] = u_trans[1]
            u_local_20[5 * inode + 2] = u_trans[2]
            u_local_20[5 * inode + 3] = u_rot[0]
            u_local_20[5 * inode + 4] = u_rot[1]

        if abs(zoffset) > 0.0:
            for inode in range(4):
                r = 5 * inode
                u_local_20[r + 0] += zoffset * u_local_20[r + 4]
                u_local_20[r + 1] -= zoffset * u_local_20[r + 3]

        u_mem = np.zeros(8)
        for inode in range(4):
            u_mem[2 * inode] = u_local_20[5 * inode]
            u_mem[2 * inode + 1] = u_local_20[5 * inode + 1]

        u_bend = np.zeros(12)
        for inode in range(4):
            u_bend[3 * inode] = u_local_20[5 * inode + 2]
            u_bend[3 * inode + 1] = u_local_20[5 * inode + 3]
            u_bend[3 * inode + 2] = u_local_20[5 * inode + 4]

        # Strains at centroid (element frame)
        dN_c = _dshape_quad4(0.0, 0.0)
        J_c, det_J_c = _jacobian(dN_c, x_local, y_local)
        J_inv_c = _jacobian_inv(J_c, det_J_c)
        dN_dxy_c = J_inv_c @ dN_c

        Bm_c = _membrane_B(dN_dxy_c[0], dN_dxy_c[1])
        Bb_c = _bending_B(dN_dxy_c[0], dN_dxy_c[1])

        strain_membrane = Bm_c @ u_mem
        curvature = Bb_c @ u_bend

        # Stress in element frame: sigma = Q_elem @ strain_elem
        # A_mat from _get_ABD_for_element is already in element frame (rotated by theta).
        if thickness > 0.0:
            Q_elem = A_mat / thickness
        else:
            Q_elem = A_mat

        stress_data = np.zeros((2, 8))
        for i_fiber, z_fiber in enumerate([z1, z2]):
            strain_fiber = strain_membrane + z_fiber * curvature
            stress_fiber = Q_elem @ strain_fiber
            sxx, syy, sxy = stress_fiber
            angle, s_major, s_minor = _principal_stresses(sxx, syy, sxy)
            svm = _von_mises(sxx, syy, sxy)
            stress_data[i_fiber] = [z_fiber, sxx, syy, sxy, angle, s_major, s_minor, svm]
        results[eid] = stress_data

    return results



def get_theta_from_theta_mcid(model: BDF,
                              elem: CQUAD4 | CTRIA3,
                              normal: np.ndarray) -> np.ndarray:
    thetas = np.radians(elem.theta)
    mcids = elem.mcid
    is_theta = (mcids == -1)
    is_mcid = ~is_theta

    theta_mat = thetas
    if is_mcid.sum():
        # project mcid onto the elements for the mcids
        mcids2 = mcids.copy()
        mcids2[is_theta] = 0
        mcids_ref = model.coord.slice_card_by_id(mcids2)

        i_mcids = mcids_ref.i[is_mcid, :]
        normals = normal[is_mcid, :]
        ihats = ihat[is_mcid,:]

        i_projs = np.einsum('ij,ij->i', i_mcids, normals) * normals
        i_proj_norms = np.linalg.norm(i_proj, axis=1)

        cos_theta = np.einsum('ij,ij->i', ihats, i_projs)

        axb = np.cross(ihats, i_projs, axis=0)
        sin_theta = np.einsum('ij,ij->i', axb, normals)
        theta_mat[is_mcid] = np.arctan2(sin_theta, cos_theta)

        ifailed = np.where(i_proj_norms < 1e-10)
        if len(ifailed):
            eids_mcid = eids[is_mcid]
            eids_failed = eids_mcid[ifailed]
            raise FatalError(f'eids={eids_failed} failed their MCID projection')

    #is_zero = np.isnan(thetas) & (mcids == -1)
    #thetas[is_zero] = 0.0
    if np.any(np.isnan(theta_mat)):
        print('theta_mat', theta_mat)
        raise RuntimeError('incorrect application of theta_mats; found nan')
    return theta_mat

def recover_shell_stress_ctria3(
    model: BDF,
    dof_map: DOF_MAP,
    xb: np.ndarray,
) -> dict[int, np.ndarray]:
    """Recover CTRIA3 stresses at centroid, top and bottom fibers.

    Returns
    -------
    results : dict[eid, (2, 8) ndarray]
        Row 0 = bottom (z1), Row 1 = top (z2).
        Columns: [fiber_dist, σxx, σyy, τxy, angle, σ_major, σ_minor, σ_vm]
    """
    ctria3 = model.ctria3
    if ctria3.n == 0:
        return {}

    pshell = model.pshell
    pcomp = model.pcomp
    pcompg = model.pcompg

    grid = model.grid
    nodes = ctria3.nodes
    inid = grid.index(nodes)
    xyz = grid.xyz_cid0()

    p1 = xyz[inid[:, 0], :]
    p2 = xyz[inid[:, 1], :]
    p3 = xyz[inid[:, 2], :]

    v12 = p2 - p1
    v13 = p3 - p1
    normal = np.cross(v12, v13, axis=1)
    ni = np.linalg.norm(normal, axis=1)
    normal /= ni[:, np.newaxis]

    ihat = v12 / np.linalg.norm(v12, axis=1)[:, np.newaxis]
    jhat = np.cross(normal, ihat, axis=1)
    jhat /= np.linalg.norm(jhat, axis=1)[:, np.newaxis]

    results = {}

    eids = ctria3.element_id
    pids = ctria3.property_id
    zoffsets = ctria3.zoffset

    thetas = get_theta_from_theta_mcid(model, cquad4, normal)

    zoffsets[np.isnan(zoffsets)] = 0.0
    common_pshell_pids = np.intersect1d(pshell.property_id, pids)
    is_pshells = np.array([pid in common_pshell_pids for pid in pids])

    for i_elem, eid, pid, theta, zoffset, is_pshell in zip(
            count(), eids, pids, thetas, zoffsets, is_pshells):

        Ti = np.vstack([ihat[i_elem], jhat[i_elem], normal[i_elem]])

        xy = np.zeros((3, 3))
        xy[0] = Ti @ p1[i_elem]
        xy[1] = Ti @ p2[i_elem]
        xy[2] = Ti @ p3[i_elem]
        x_local = xy[:, 0]
        y_local = xy[:, 1]

        A_mat, B_mat, D_mat, Ds_mat, thickness = _get_ABD_for_element(
            model, pid, pshell, pcomp, pcompg)

        if is_pshell:
            iprop = np.searchsorted(pshell.property_id, pid)
            z1, z2 = pshell.z[iprop, :]
        else:
            z1 = -thickness / 2.0
            z2 = thickness / 2.0

        elem_nodes = nodes[i_elem]
        u_local_15 = np.zeros(15)
        for inode in range(3):
            nid = elem_nodes[inode]
            i_dof = dof_map[(nid, 1)]
            u_g = xb[i_dof : i_dof + 6]
            u_trans = Ti @ u_g[:3]
            u_rot = Ti[0:2, :] @ u_g[3:6]
            u_local_15[5 * inode] = u_trans[0]
            u_local_15[5 * inode + 1] = u_trans[1]
            u_local_15[5 * inode + 2] = u_trans[2]
            u_local_15[5 * inode + 3] = u_rot[0]
            u_local_15[5 * inode + 4] = u_rot[1]

        if abs(zoffset) > 0.0:
            for inode in range(3):
                r = 5 * inode
                u_local_15[r + 0] += zoffset * u_local_15[r + 4]
                u_local_15[r + 1] -= zoffset * u_local_15[r + 3]

        # CST membrane strain (constant)
        area = 0.5 * abs(
            (x_local[1] - x_local[0]) * (y_local[2] - y_local[0])
            - (x_local[2] - x_local[0]) * (y_local[1] - y_local[0])
        )
        dNdx = np.array(
            [y_local[1] - y_local[2], y_local[2] - y_local[0], y_local[0] - y_local[1]]
        ) / (2.0 * area)
        dNdy = np.array(
            [x_local[2] - x_local[1], x_local[0] - x_local[2], x_local[1] - x_local[0]]
        ) / (2.0 * area)

        u_mem = np.zeros(6)
        for inode in range(3):
            u_mem[2 * inode] = u_local_15[5 * inode]
            u_mem[2 * inode + 1] = u_local_15[5 * inode + 1]

        Bm = np.zeros((3, 6))
        Bm[0, 0::2] = dNdx
        Bm[1, 1::2] = dNdy
        Bm[2, 0::2] = dNdy
        Bm[2, 1::2] = dNdx

        strain_membrane = Bm @ u_mem

        # DKT bending curvature at centroid (simplified: use average of nodal rotations)
        # For CST-based approach: curvature from constant strain triangle is zero
        # For DKT: evaluate at centroid (1/3, 1/3, 1/3)
        # Simplified: use linear interpolation of rotations for curvature
        u_bend = np.zeros(9)
        for inode in range(3):
            u_bend[3 * inode] = u_local_15[5 * inode + 2]
            u_bend[3 * inode + 1] = u_local_15[5 * inode + 3]
            u_bend[3 * inode + 2] = u_local_15[5 * inode + 4]

        # Curvature from rotations: kxx = dtheta_y/dx, kyy = -dtheta_x/dy,
        # kxy = dtheta_y/dy - dtheta_x/dx
        theta_x = np.array([u_local_15[5 * i + 3] for i in range(3)])
        theta_y = np.array([u_local_15[5 * i + 4] for i in range(3)])
        kxx = dNdx @ theta_y
        kyy = -(dNdy @ theta_x)
        kxy = dNdy @ theta_y - dNdx @ theta_x
        curvature = np.array([kxx, kyy, kxy])

        if thickness > 0.0:
            Q = A_mat / thickness
        else:
            Q = A_mat

        stress_data = np.zeros((2, 8))
        for i_fiber, z_fiber in enumerate([z1, z2]):
            strain_fiber = strain_membrane + z_fiber * curvature
            stress_fiber = Q @ strain_fiber
            sxx, syy, sxy = stress_fiber
            angle, s_major, s_minor = _principal_stresses(sxx, syy, sxy)
            svm = _von_mises(sxx, syy, sxy)
            stress_data[i_fiber] = [
                z_fiber, sxx, syy, sxy, angle, s_major, s_minor, svm]
        results[eid] = stress_data

    return results


def recover_shell_force_cquad4(
    model: BDF,
    dof_map: DOF_MAP,
    xb: np.ndarray,
) -> dict[int, np.ndarray]:
    """Recover CQUAD4 element forces (stress resultants) at centroid.

    Returns
    -------
    results : dict[eid, (8,) ndarray]
        [Nxx, Nyy, Nxy, Mxx, Myy, Mxy, Qx, Qy]
    """
    stress_results = recover_shell_stress_cquad4(model, dof_map, xb)
    if not stress_results:
        return {}

    pshell = model.pshell
    cquad4 = model.cquad4
    force_results = {}

    eids = list(stress_results)
    ieid = np.searchsorted(cquad4.element_id, eids)
    pids = cquad4.property_id[ieid]
    for (eid, stress_data), pid in zip(stress_results.items(), pids):
        iprop = np.searchsorted(pshell.property_id, pid)
        thickness = float(pshell.t[iprop]) if iprop < pshell.n else 0.1
        z1 = stress_data[0, 0]
        z2 = stress_data[1, 0]

        # Stress resultants from fiber stresses:
        # N = integral(sigma dz) ≈ (sigma_bot + sigma_top)/2 * t
        # M = integral(sigma * z dz) ≈ (sigma_top - sigma_bot)/(z2 - z1) * t^3/12
        s_bot = stress_data[0, 1:4]  # [sxx, syy, sxy] at z1
        s_top = stress_data[1, 1:4]  # [sxx, syy, sxy] at z2

        Nxx = 0.5 * (s_bot[0] + s_top[0]) * thickness
        Nyy = 0.5 * (s_bot[1] + s_top[1]) * thickness
        Nxy = 0.5 * (s_bot[2] + s_top[2]) * thickness

        dz = z2 - z1 if abs(z2 - z1) > 1e-30 else thickness
        Mxx = (s_top[0] - s_bot[0]) / dz * thickness**3 / 12.0
        Myy = (s_top[1] - s_bot[1]) / dz * thickness**3 / 12.0
        Mxy = (s_top[2] - s_bot[2]) / dz * thickness**3 / 12.0

        # Transverse shear forces (Qx, Qy) would require shear strain recovery
        # For now, set to 0 (TODO: compute from transverse shear strains)
        Qx = 0.0
        Qy = 0.0
        force_results[eid] = np.array([Nxx, Nyy, Nxy, Mxx, Myy, Mxy, Qx, Qy])

    return force_results


def recover_shell_strain_energy_cquad4(
    model: BDF,
    dof_map: DOF_MAP,
    xb: np.ndarray) -> dict[int, float]:
    """Recover CQUAD4 element strain energy: SE = 0.5 * u^T @ K @ u.

    Returns
    -------
    results : dict[eid, float]
    """
    cquad4 = model.cquad4
    if cquad4.n == 0:
        return {}

    pshell = model.pshell
    pcomp = model.pcomp
    pcompg = model.pcompg

    grid = model.grid
    nodes = cquad4.nodes
    inid = grid.index(nodes)
    xyz = grid.xyz_cid0()

    p1 = xyz[inid[:, 0], :]
    p2 = xyz[inid[:, 1], :]
    p3 = xyz[inid[:, 2], :]
    p4 = xyz[inid[:, 3], :]

    v13 = p1 - p3
    v24 = p2 - p4
    normal = np.cross(v13, v24)
    ni = np.linalg.norm(normal, axis=1)
    normal /= ni[:, np.newaxis]

    ihat = p2 - p1
    ihat /= np.linalg.norm(ihat, axis=1)[:, np.newaxis]
    jhat = np.cross(normal, ihat, axis=1)
    jhat /= np.linalg.norm(jhat, axis=1)[:, np.newaxis]

    results = {}

    eids = cquad4.element_id
    pids = cquad4.property_id
    zoffsets = cquad4.zoffset
    assert np.abs(zoffsets.min()) >= 0., zoffsets
    apply_zoffset = (np.abs(zoffset) > 0.0)

    for i, eid, pid, zoffset, apply_zoffset in zip(count(), eids, pids, zoffsets, apply_zoffsets):
        Ti = np.vstack([ihat[i], jhat[i], normal[i]])

        xy = np.zeros((4, 3))
        xy[0] = Ti @ p1[i]
        xy[1] = Ti @ p2[i]
        xy[2] = Ti @ p3[i]
        xy[3] = Ti @ p4[i]
        x_local = xy[:, 0]
        y_local = xy[:, 1]

        A_mat, B_mat, D_mat, Ds_mat, thickness = _get_ABD_for_element(
            model, pid, pshell, pcomp, pcompg)

        Ke = macn2_stiffness(x_local, y_local, A_mat, D_mat, Ds_mat, B_mat)

        #if not np.isnan(zoffset) and abs(zoffset) > 0.0:
        if apply_zoffset:
            Ke = _apply_shell_offset(Ke, zoffset)

        # Extract local displacements
        elem_nodes = nodes[i]
        u_local = np.zeros(20)
        for inode in range(4):
            nid = elem_nodes[inode]
            i_dof = dof_map[(nid, 1)]
            u_g = xb[i_dof : i_dof + 6]
            u_trans = Ti @ u_g[:3]
            u_rot = Ti[0:2, :] @ u_g[3:6]
            u_local[5 * inode] = u_trans[0]
            u_local[5 * inode + 1] = u_trans[1]
            u_local[5 * inode + 2] = u_trans[2]
            u_local[5 * inode + 3] = u_rot[0]
            u_local[5 * inode + 4] = u_rot[1]

        se = 0.5 * u_local @ Ke @ u_local
        results[eid] = float(se)

    return results
