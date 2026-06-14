"""Assembly routines for 3D solid elements (CHEXA, CTETRA, CPENTA).

Wires element-level stiffness, mass, and geometric stiffness matrices
into the global DOF assembly.
"""

from __future__ import annotations
from itertools import count
from typing import TYPE_CHECKING

import numpy as np

from .tetra import (
    _solid_B_matrix,
    _tetra4_shape_functions,
    _tetra10_shape_functions,
    _TETRA4_GAUSS,
    _TETRA10_GAUSS,
    tetra4_stiffness,
    tetra4_mass,
    tetra4_geometric_stiffness,
    tetra10_stiffness,
    tetra10_mass,
    tetra10_geometric_stiffness,
)
from .penta import (
    _penta6_shape_functions,
    _penta15_shape_functions,
    _PENTA6_FULL_GAUSS,
    _PENTA15_GAUSS,
    penta6_stiffness,
    penta6_mass,
    penta6_geometric_stiffness,
    penta15_stiffness,
    penta15_mass,
    penta15_geometric_stiffness,
)
from .hexa import (
    _hexa8_shape_functions,
    _hexa20_shape_functions,
    _isotropic_constitutive,
    hexa8_stiffness,
    hexa8_mass,
    hexa8_geometric_stiffness,
    hexa20_stiffness,
    hexa20_mass,
    hexa20_geometric_stiffness,
)

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.solver.build_stiffness import _COOAccumulator
    DOF_MAP = dict[tuple[int, int], int]


def _get_solid_material_D(model: BDF, property_id: int) -> np.ndarray:
    """Get the 6x6 constitutive matrix for a PSOLID property.

    Parameters
    ----------
    model : BDF
        The BDF model.
    property_id : int
        Property ID (must reference a PSOLID card with MAT1).

    Returns
    -------
    D : (6, 6) constitutive matrix
    """
    psolid = model.psolid
    iprop = np.searchsorted(psolid.property_id, property_id)
    mid = int(psolid.material_id[iprop])

    mat1 = model.mat1
    imat = np.searchsorted(mat1.material_id, mid)
    E = float(mat1.E[imat])
    nu = float(mat1.nu[imat])
    return _isotropic_constitutive(E, nu)


def _get_solid_rho(model: BDF, property_id: int) -> float:
    """Get the density for a PSOLID property."""
    psolid = model.psolid
    iprop = np.searchsorted(psolid.property_id, property_id)
    mid = int(psolid.material_id[iprop])

    mat1 = model.mat1
    imat = np.searchsorted(mat1.material_id, mid)
    return float(mat1.rho[imat])


def _solid_dof_indices(dof_map: "DOF_MAP", node_ids: np.ndarray) -> list[int]:
    """Build DOF index list for solid element (3 translational DOFs per node).

    Parameters
    ----------
    dof_map : DOF_MAP
        Maps (nid, dof) to global DOF index.
    node_ids : (nnodes,)
        Node IDs for this element.

    Returns
    -------
    n_ijv : list[int]
        Global DOF indices [T1_n1, T2_n1, T3_n1, T1_n2, ...].
    """
    n_ijv = []
    for nid in node_ids:
        gi = dof_map[(nid, 1)]
        n_ijv.extend([gi, gi + 1, gi + 2])
    return n_ijv


def _get_element_coords(model: BDF, node_ids: np.ndarray) -> np.ndarray:
    """Get nodal coordinates for an element's nodes.

    Parameters
    ----------
    model : BDF
        The BDF model.
    node_ids : (nnodes,)
        Node IDs.

    Returns
    -------
    coords : (nnodes, 3)
        Nodal coordinates in basic frame.
    """
    grid = model.grid
    all_nids = grid.node_id
    xyz_cid0 = grid.xyz_cid0()
    inodes = np.searchsorted(all_nids, node_ids)
    return xyz_cid0[inodes]


def build_kbb_chexa(model: BDF, coo: "_COOAccumulator", dof_map: "DOF_MAP") -> int:
    """Assemble CHEXA element stiffness matrices into the global COO accumulator.

    Parameters
    ----------
    model : BDF
        The BDF model.
    coo : _COOAccumulator
        Sparse COO accumulator for global stiffness.
    dof_map : DOF_MAP
        Maps (nid, dof) to global DOF index.

    Returns
    -------
    nelements : int
        Number of CHEXA elements assembled.
    """
    elem = model.chexa
    if elem.n == 0:
        return 0

    log = model.log
    nelements = 0

    eids = elem.element_id
    pids = elem.property_id
    nodes = elem.nodes
    for i, eid, pid, all_nodes in zip(count(), eids, pids, nodes):
        # Determine if HEXA8 or HEXA20 based on midside nodes
        base_nodes = all_nodes[:8]
        midside_nodes = all_nodes[8:]
        has_midside = midside_nodes.max() > 0 if len(midside_nodes) > 0 else False

        D = _get_solid_material_D(model, pid)

        if has_midside:
            active_nodes = all_nodes[all_nodes > 0]
            if len(active_nodes) == 20:
                coords = _get_element_coords(model, active_nodes)
                Ke = hexa20_stiffness(coords, D)
                n_ijv = _solid_dof_indices(dof_map, active_nodes)
            else:
                # Partial midside nodes - fall back to HEXA8
                coords = _get_element_coords(model, base_nodes)
                Ke = hexa8_stiffness(coords, D)
                n_ijv = _solid_dof_indices(dof_map, base_nodes)
        else:
            coords = _get_element_coords(model, base_nodes)
            Ke = hexa8_stiffness(coords, D)
            n_ijv = _solid_dof_indices(dof_map, base_nodes)

        coo.add_matrix(n_ijv, Ke)
        nelements += 1

    log.debug(f"  assembled {nelements} CHEXA elements")
    return nelements


def build_kbb_ctetra(model: BDF, coo: "_COOAccumulator", dof_map: "DOF_MAP") -> int:
    """Assemble CTETRA element stiffness matrices into the global COO accumulator."""
    elem = model.ctetra
    if elem.n == 0:
        return 0

    log = model.log
    nelements = 0

    eids = elem.element_id
    pids = elem.property_id
    nodes = elem.nodes
    for i, eid, pid, all_nodes in zip(count(), eids, pids, nodes):
        base_nodes = all_nodes[:4]
        midside_nodes = all_nodes[4:]
        has_midside = midside_nodes.max() > 0 if len(midside_nodes) > 0 else False

        D = _get_solid_material_D(model, pid)

        if has_midside:
            active_nodes = all_nodes[all_nodes > 0]
            if len(active_nodes) == 10:
                coords = _get_element_coords(model, active_nodes)
                Ke = tetra10_stiffness(coords, D)
                n_ijv = _solid_dof_indices(dof_map, active_nodes)
            else:
                coords = _get_element_coords(model, base_nodes)
                Ke = tetra4_stiffness(coords, D)
                n_ijv = _solid_dof_indices(dof_map, base_nodes)
        else:
            coords = _get_element_coords(model, base_nodes)
            Ke = tetra4_stiffness(coords, D)
            n_ijv = _solid_dof_indices(dof_map, base_nodes)

        coo.add_matrix(n_ijv, Ke)
        nelements += 1

    log.debug(f"  assembled {nelements} CTETRA elements")
    return nelements


def build_kbb_cpenta(model: BDF, coo: "_COOAccumulator", dof_map: "DOF_MAP") -> int:
    """Assemble CPENTA element stiffness matrices into the global COO accumulator."""
    elem = model.cpenta
    if elem.n == 0:
        return 0

    log = model.log
    nelements = 0

    eids = elem.element_id
    pids = elem.property_id
    nodes = elem.nodes
    for i, eid, pid, all_nodes in zip(count(), eids, pids, nodes):
        base_nodes = all_nodes[:6]
        midside_nodes = all_nodes[6:]
        has_midside = midside_nodes.max() > 0 if len(midside_nodes) > 0 else False

        D = _get_solid_material_D(model, pid)

        if has_midside:
            active_nodes = all_nodes[all_nodes > 0]
            if len(active_nodes) == 15:
                coords = _get_element_coords(model, active_nodes)
                Ke = penta15_stiffness(coords, D)
                n_ijv = _solid_dof_indices(dof_map, active_nodes)
            else:
                coords = _get_element_coords(model, base_nodes)
                Ke = penta6_stiffness(coords, D)
                n_ijv = _solid_dof_indices(dof_map, base_nodes)
        else:
            coords = _get_element_coords(model, base_nodes)
            Ke = penta6_stiffness(coords, D)
            n_ijv = _solid_dof_indices(dof_map, base_nodes)

        coo.add_matrix(n_ijv, Ke)
        nelements += 1

    log.debug(f"  assembled {nelements} CPENTA elements")
    return nelements


def build_mbb_solids(
    model: BDF,
    coo: "_COOAccumulator",
    dof_map: "DOF_MAP",
) -> float:
    """Assemble solid element mass matrices into the global COO accumulator.

    Returns
    -------
    mass_total : float
        Total mass contributed by solid elements.
    """
    mass_total = 0.0
    mass_total += _build_mbb_chexa(model, coo, dof_map)
    mass_total += _build_mbb_ctetra(model, coo, dof_map)
    mass_total += _build_mbb_cpenta(model, coo, dof_map)
    return mass_total


def _build_mbb_chexa(
    model: BDF,
    coo: "_COOAccumulator",
    dof_map: "DOF_MAP",
) -> float:
    """Assemble CHEXA mass matrices."""
    elem = model.chexa
    if elem.n == 0:
        return 0.0

    mass_total = 0.0
    eids = elem.element_id
    pids = elem.property_id
    nodes = elem.nodes
    for i, eid, pid, all_nodes in zip(count(), eids, pids, nodes):
        base_nodes = all_nodes[:8]
        midside_nodes = all_nodes[8:]
        has_midside = midside_nodes.max() > 0 if len(midside_nodes) > 0 else False

        rho = _get_solid_rho(model, pid)
        if rho == 0.0:
            continue

        if has_midside and len(all_nodes[all_nodes > 0]) == 20:
            active_nodes = all_nodes[all_nodes > 0]
            coords = _get_element_coords(model, active_nodes)
            Me = hexa20_mass(coords, rho)
            n_ijv = _solid_dof_indices(dof_map, active_nodes)
        else:
            coords = _get_element_coords(model, base_nodes)
            Me = hexa8_mass(coords, rho)
            n_ijv = _solid_dof_indices(dof_map, base_nodes)

        coo.add_matrix(n_ijv, Me)
        mass_total += Me.trace() / 3.0

    return mass_total


def _build_mbb_ctetra(
    model: BDF,
    coo: "_COOAccumulator",
    dof_map: "DOF_MAP",
) -> float:
    """Assemble CTETRA mass matrices."""
    elem = model.ctetra
    if elem.n == 0:
        return 0.0

    mass_total = 0.0
    eids = elem.element_id
    pids = elem.property_id
    nodes = elem.nodes
    for i, eid, pid, all_nodes in zip(count(), eids, pids, nodes):
        base_nodes = all_nodes[:4]
        midside_nodes = all_nodes[4:]
        has_midside = midside_nodes.max() > 0 if len(midside_nodes) > 0 else False

        rho = _get_solid_rho(model, pid)
        if rho == 0.0:
            continue

        if has_midside and len(all_nodes[all_nodes > 0]) == 10:
            active_nodes = all_nodes[all_nodes > 0]
            coords = _get_element_coords(model, active_nodes)
            Me = tetra10_mass(coords, rho)
            n_ijv = _solid_dof_indices(dof_map, active_nodes)
        else:
            coords = _get_element_coords(model, base_nodes)
            Me = tetra4_mass(coords, rho)
            n_ijv = _solid_dof_indices(dof_map, base_nodes)

        coo.add_matrix(n_ijv, Me)
        mass_total += Me.trace() / 3.0

    return mass_total


def _build_mbb_cpenta(
    model: BDF,
    coo: "_COOAccumulator",
    dof_map: "DOF_MAP",
) -> float:
    """Assemble CPENTA mass matrices."""
    elem = model.cpenta
    if elem.n == 0:
        return 0.0

    mass_total = 0.0
    eids = elem.element_id
    pids = elem.property_id
    nodes = elem.nodes
    for i, eid, pid, all_nodes in zip(count(), eids, pids, nodes):
        base_nodes = all_nodes[:6]
        midside_nodes = all_nodes[6:]
        has_midside = midside_nodes.max() > 0 if len(midside_nodes) > 0 else False

        rho = _get_solid_rho(model, pid)
        if rho == 0.0:
            continue

        if has_midside and len(all_nodes[all_nodes > 0]) == 15:
            active_nodes = all_nodes[all_nodes > 0]
            coords = _get_element_coords(model, active_nodes)
            Me = penta15_mass(coords, rho)
            n_ijv = _solid_dof_indices(dof_map, active_nodes)
        else:
            coords = _get_element_coords(model, base_nodes)
            Me = penta6_mass(coords, rho)
            n_ijv = _solid_dof_indices(dof_map, base_nodes)

        coo.add_matrix(n_ijv, Me)
        mass_total += Me.trace() / 3.0

    return mass_total


def build_KDgg_solids(
    model: BDF,
    KDgg,
    dof_map: "DOF_MAP",
    u_global: np.ndarray,
) -> None:
    """Assemble geometric stiffness for solid elements from preload displacements.

    Parameters
    ----------
    model : BDF
        The BDF model.
    KDgg : dok_matrix
        Global geometric stiffness matrix (modified in place).
    dof_map : DOF_MAP
        Maps (nid, dof) to global DOF index.
    u_global : (ndof,)
        Global displacement vector from preload.
    """
    _build_KDgg_chexa(model, KDgg, dof_map, u_global)
    _build_KDgg_ctetra(model, KDgg, dof_map, u_global)
    _build_KDgg_cpenta(model, KDgg, dof_map, u_global)


def _compute_element_stress(
    model: BDF,
    coords: np.ndarray,
    D: np.ndarray,
    dof_map: "DOF_MAP",
    node_ids: np.ndarray,
    u_global: np.ndarray,
    nnodes: int,
    shape_func,
    gauss_data,
) -> np.ndarray:
    """Compute average element stress from preload displacements.

    Returns
    -------
    stress_avg : (6,) average stress [sxx, syy, szz, sxy, syz, sxz]
    """
    # Extract element displacements
    ue = np.zeros(3 * nnodes)
    for k, nid in enumerate(node_ids):
        gi = dof_map[(int(nid), 1)]
        ue[3 * k : 3 * k + 3] = u_global[gi : gi + 3]

    gauss_pts, weights = gauss_data
    stress_sum = np.zeros(6)
    vol = 0.0

    for gpt, w in zip(gauss_pts, weights):
        N, dNdnat = shape_func(*gpt)
        J = dNdnat @ coords
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        dN_dxyz = Jinv @ dNdnat
        B = _solid_B_matrix(dN_dxyz, nnodes)
        strain = B @ ue
        stress = D @ strain
        stress_sum += stress * detJ * w
        vol += detJ * w

    if vol > 0.0:
        return stress_sum / vol
    return np.zeros(6)


def _build_KDgg_chexa(model: BDF, KDgg, dof_map, u_global):
    """Assemble CHEXA geometric stiffness."""
    elem = model.chexa
    if elem.n == 0:
        return

    k = 1.0 / np.sqrt(3.0)
    gauss_pts_8 = np.array(
        [
            [-k, -k, -k],
            [k, -k, -k],
            [k, k, -k],
            [-k, k, -k],
            [-k, -k, k],
            [k, -k, k],
            [k, k, k],
            [-k, k, k],
        ]
    )
    weights_8 = np.ones(8)

    def shape_8(xi, eta, mu):
        N, dNdxi, dNdeta, dNdmu = _hexa8_shape_functions(xi, eta, mu)
        dNdnat = np.array([dNdxi, dNdeta, dNdmu])
        return N, dNdnat

    eids = elem.element_id
    pids = elem.property_id
    nodes = elem.nodes
    for i, eid, pid, all_nodes in zip(count(), eids, pids, nodes):
        base_nodes = all_nodes[:8]
        midside_nodes = all_nodes[8:]
        has_midside = midside_nodes.max() > 0 if len(midside_nodes) > 0 else False

        D = _get_solid_material_D(model, pid)

        if has_midside and len(all_nodes[all_nodes > 0]) == 20:
            # HEXA20 geometric stiffness
            active_nodes = all_nodes[all_nodes > 0]
            coords = _get_element_coords(model, active_nodes)
            k3 = np.sqrt(3.0 / 5.0)
            pts_1d = np.array([-k3, 0.0, k3])
            wts_1d = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])
            gpts = []
            gwts = []
            for ii in range(3):
                for jj in range(3):
                    for mm in range(3):
                        gpts.append([pts_1d[ii], pts_1d[jj], pts_1d[mm]])
                        gwts.append(wts_1d[ii] * wts_1d[jj] * wts_1d[mm])
            gauss_data = (np.array(gpts), np.array(gwts))

            stress = _compute_element_stress(
                model,
                coords,
                D,
                dof_map,
                active_nodes,
                u_global,
                20,
                _hexa20_shape_functions,
                gauss_data,
            )
            KDe = hexa20_geometric_stiffness(coords, stress)
            n_ijv = _solid_dof_indices(dof_map, active_nodes)
        else:
            coords = _get_element_coords(model, base_nodes)
            gauss_data = (gauss_pts_8, weights_8)
            stress = _compute_element_stress(
                model,
                coords,
                D,
                dof_map,
                base_nodes,
                u_global,
                8,
                shape_8,
                gauss_data,
            )
            KDe = hexa8_geometric_stiffness(coords, stress)
            n_ijv = _solid_dof_indices(dof_map, base_nodes)

        # Scatter into global KDgg
        ndof_e = len(n_ijv)
        for ii in range(ndof_e):
            for jj in range(ndof_e):
                if KDe[ii, jj] != 0.0:
                    KDgg[n_ijv[ii], n_ijv[jj]] += KDe[ii, jj]


def _build_KDgg_ctetra(model: BDF, KDgg, dof_map, u_global):
    """Assemble CTETRA geometric stiffness."""
    elem = model.ctetra
    if elem.n == 0:
        return

    for i in range(elem.n):
        pid = int(elem.property_id[i])
        all_nodes = elem.nodes[i]
        base_nodes = all_nodes[:4]
        midside_nodes = all_nodes[4:]
        has_midside = midside_nodes.max() > 0 if len(midside_nodes) > 0 else False

        D = _get_solid_material_D(model, pid)

        if has_midside and len(all_nodes[all_nodes > 0]) == 10:
            active_nodes = all_nodes[all_nodes > 0]
            coords = _get_element_coords(model, active_nodes)
            stress = _compute_element_stress(
                model,
                coords,
                D,
                dof_map,
                active_nodes,
                u_global,
                10,
                _tetra10_shape_functions,
                _TETRA10_GAUSS,
            )
            KDe = tetra10_geometric_stiffness(coords, stress)
            n_ijv = _solid_dof_indices(dof_map, active_nodes)
        else:
            coords = _get_element_coords(model, base_nodes)
            stress = _compute_element_stress(
                model,
                coords,
                D,
                dof_map,
                base_nodes,
                u_global,
                4,
                _tetra4_shape_functions,
                _TETRA4_GAUSS,
            )
            KDe = tetra4_geometric_stiffness(coords, stress)
            n_ijv = _solid_dof_indices(dof_map, base_nodes)

        ndof_e = len(n_ijv)
        for ii in range(ndof_e):
            for jj in range(ndof_e):
                if KDe[ii, jj] != 0.0:
                    KDgg[n_ijv[ii], n_ijv[jj]] += KDe[ii, jj]


def _build_KDgg_cpenta(model: BDF, KDgg, dof_map, u_global):
    """Assemble CPENTA geometric stiffness."""
    elem = model.cpenta
    if elem.n == 0:
        return

    for i in range(elem.n):
        pid = int(elem.property_id[i])
        all_nodes = elem.nodes[i]
        base_nodes = all_nodes[:6]
        midside_nodes = all_nodes[6:]
        has_midside = midside_nodes.max() > 0 if len(midside_nodes) > 0 else False

        D = _get_solid_material_D(model, pid)

        if has_midside and len(all_nodes[all_nodes > 0]) == 15:
            active_nodes = all_nodes[all_nodes > 0]
            coords = _get_element_coords(model, active_nodes)
            stress = _compute_element_stress(
                model,
                coords,
                D,
                dof_map,
                active_nodes,
                u_global,
                15,
                _penta15_shape_functions,
                _PENTA15_GAUSS,
            )
            KDe = penta15_geometric_stiffness(coords, stress)
            n_ijv = _solid_dof_indices(dof_map, active_nodes)
        else:
            coords = _get_element_coords(model, base_nodes)
            stress = _compute_element_stress(
                model,
                coords,
                D,
                dof_map,
                base_nodes,
                u_global,
                6,
                _penta6_shape_functions,
                _PENTA6_FULL_GAUSS,
            )
            KDe = penta6_geometric_stiffness(coords, stress)
            n_ijv = _solid_dof_indices(dof_map, base_nodes)

        ndof_e = len(n_ijv)
        for ii in range(ndof_e):
            for jj in range(ndof_e):
                if KDe[ii, jj] != 0.0:
                    KDgg[n_ijv[ii], n_ijv[jj]] += KDe[ii, jj]
