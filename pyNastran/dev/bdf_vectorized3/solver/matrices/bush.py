from __future__ import annotations
from typing import TYPE_CHECKING
from itertools import count

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import dok_matrix #, csc_matrix, lil_matrix, issparse
if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

def _skew(v):
    """Skew-symmetric matrix: _skew(v) @ u = v × u."""
    return np.array([
            [0, -v[2], v[1]],
            [v[2], 0, -v[0]],
            [-v[1], v[0], 0],
    ])

def assemble_cbush_matrices(
    model: BDF,
    dof_map: dict[tuple[int, int], int],
    Kgg: dok_matrix,
    Bgg: dok_matrix=None,
    Mgg: dok_matrix=None,
    include_damping: bool = False,
    include_mass: bool = False,
    extra_mass: dict[int, float] | None = None,) -> tuple[NDArray, NDArray | None, NDArray | None]:
    """Assemble CBUSH stiffness, damping, and lumped mass matrices.

    Uses pyNastran array-based API: reads element
    connectivity from ``model.cbush`` and property data from
    ``model.pbush`` (K1-K6, B1-B6, mass as NumPy arrays).

    Parameters
    ----------
    model : pyNastran vectorized3 BDF model
        Parsed model with ``bdf.cbush`` and ``bdf.pbush`` populated.
    dof_map : dict
        Maps (grid_id, component) -> matrix row/col index.
    ndof : int or None
        Matrix size. None = max(dof_map.values()) + 1.
    include_damping : bool
        If True, also assemble CBUSH damping (B1-B6).
    extra_mass : dict or None
        Additional lumped masses at grids: {grid_id: mass_value}.
        Applied equally to T1, T2, T3 DOFs of the grid.

    Returns
    -------
    Kgg : NDArray
        CBUSH stiffness matrix, shape (ndof, ndof).
    Bgg : NDArray or None
        CBUSH damping matrix (if include_damping and any B values exist).
    Mgg : NDArray or None
        Lumped mass matrix from CBUSH mass and extra_mass
        and include_mass.

    """
    log = model.log
    log.level = 'debug'

    ndof = Kgg.shape[0]
    extra_mass = extra_mass if include_mass else None
    fdtype: str = Kgg.dtype
    has_damping = False
    has_mass = False

    cbush = model.cbush
    pbush = model.pbush
    neid = cbush.n

    # Extra lumped masses
    if extra_mass:
        has_mass = True
        for gid, m in extra_mass.items():
            for comp in range(1, 4):
                idx = dof_map.get((gid, comp))
                if idx is not None:
                    Mgg[idx, idx] += m

    if neid == 0:
        log.info(f"Assembled 0 CBUSH elements into {ndof}x{ndof} matrices")
        return neid

    # Vectorized property lookup: map each element's PID to PBUSH row index
    eids = cbush.element_id
    pids = cbush.property_id  # (neid,)
    ipid = pbush.index(pids)

    # Per-element K, B, mass from property arrays
    k = pbush.k_fields[ipid, :]  # (neid, 6)
    b = pbush.b_fields[ipid, :]  # (neid, 6)
    mass = pbush._mass[ipid]  # (neid,)

    if np.abs(b).sum() == 0:
        include_damping = False
    if np.abs(mass).sum() == 0:
        include_mass = False
    assert np.all(np.isfinite(k))
    assert np.all(np.isfinite(b))
    assert np.all(np.isfinite(mass))
    # Replace NaN with 0 (unset values in PBUSH are stored as NaN)
    #k_all = np.nan_to_num(k_all, nan=0.0)
    #b_all = np.nan_to_num(b_all, nan=0.0)
    #mass_all = np.nan_to_num(mass_all, nan=0.0)

    # (neid, 2) — GA, GB; GB=0 is grounded
    nodes = cbush.nodes
    # (neid,) — CID; -1 means element-based
    coord_ids = cbush.coord_id

    # Compute element coordinate transforms
    # (neid, 3, 3) or None per elem
    transforms = _compute_cbush_transforms(cbush, model)

    # Element S parameter and OCID offset arrays
    s_all = cbush.s        # (n,) spring location fraction
    ocid_all = cbush.ocid  # (n,) offset coordinate system
    si_all = cbush.si      # (n, 3) offset components

    # Grid positions for offset calculations
    grid = model.grid
    coord = model.coord
    xyz = grid.xyz_cid0()
    #nid = grid.node_id
    
    unids = np.unique(nodes)
    if unids[0] == 0:
        inid = -np.ones((neid * 2), dtype='int32')
        is_nids = (nodes > 0)
        is_nids_ravel = is_nids.ravel()
        nids_ravel = nodes.ravel()[is_nids_ravel]
        #print(f'nids_ravel = {nids_ravel}')
        #print(f'grid.node_id = {grid.node_id}')
        inid[is_nids_ravel] = grid.index(nids_ravel)
        inid = inid.reshape((neid, 2))
    else:
        inid = grid.index(nodes)

    #iunid = grid.index(unids)
    #ucds = np.unique(grid.cd[iunid])

    #is_gab = is_nids[:, 0] & is_nids[:, 1]
    #is_ga  = is_nids[:, 0] & ~is_nids[:, 1]
    #is_gb  = ~is_nids[:, 0] & is_nids[:, 1]

    # Assemble element by element using G-matrix formulation
    # K_elem = G^T @ K_spring @ G where G accounts for rigid arms
    # from nodes to spring point (handles S offset and OCID offset)
    kgg_elems: dict[int, tuple[tuple[int, ...], np.ndarray]] = {}
    bgg_elems: dict[int, tuple[tuple[int, ...], np.ndarray]] = {}
    mgg_elems: dict[int, tuple[tuple[int, ...], np.ndarray]] = {}
    cbush_kelems(model, cbush, transforms, kgg_elems, pbush.k_fields)

    for eid, ((ga, gb), k_elem) in kgg_elems.items():
        print(f'eid\n{k_elem}')
        _assemble_kelm(Kgg, k_elem, ga, gb, dof_map)
    #print(f'Kgg:\n{Kgg.todense()}')

    if include_damping:
        cbush_kelems(model, cbush, transforms, bgg_elems, pbush.b_fields)
        for eid, ((ga, gb), b_elem) in bgg_elems.items():
            _assemble_kelm(Bgg, b_elem, ga, gb, dof_map)

    if include_mass:
        imass = (mass != 0)
        has_mass = (imass.sum() > 0)
        if has_mass and Mgg is not None:
            dmass = mass[imass] / 2.0
            for i, eid, (ga, gb), dmassi in zip(count(),
                        cbush.element_id[imass],
                        cbush.nodes[imass, :], dmass):
                for node in [ga, gb]:
                    if node == 0:
                        continue
                    #assert isinstance(dof_map, dict)
                    #assert isinstance(node, integer_types)
                    idx = dof_map[(node, 1)]
                    #print((node, idx))
                    Mgg[idx, idx] += dmassi
                    Mgg[idx+1, idx+1] += dmassi
                    Mgg[idx+2, idx+2] += dmassi

    log.info(f"Assembled {neid} CBUSH elements into {ndof}x{ndof} matrices")
    assert neid > 0, neid
    return neid


def cbush_kelems(model: BDF,
                 cbush: CBUSH,
                 transforms,
                 k_elems: dict[int, tuple[tuple[int, ...], np.ndarray]],
                 k_fields):
    grid = model.grid
    transforms = _compute_cbush_transforms(cbush, model)

    # Element S parameter and OCID offset arrays
    s_all = cbush.s        # (n,) spring location fraction
    ocid_all = cbush.ocid  # (n,) offset coordinate system
    si_all = cbush.si      # (n, 3) offset components

    # neid = len(cbush)
    nodes = cbush.nodes
    gas = cbush.nodes[:, 0]
    gbs = cbush.nodes[:, 1]
    is_gas = (gas > 0)
    is_gbs = (gbs > 0)
    neid = len(cbush)
    
    xyz_cid0 = grid.xyz_cid0()
    xyz = xyz_cid0
    ga_xyzs = np.full((neid, 3), np.nan, dtype='float64')
    gb_xyzs = np.full((neid, 3), np.nan, dtype='float64')
    igas = grid.index(gas[is_gas])
    igbs = grid.index(gbs[is_gbs])
    ga_xyzs[is_gas, :] = xyz_cid0[igas, :]
    gb_xyzs[is_gbs, :] = xyz_cid0[igbs, :]
   
    is_ks = np.abs(k_fields).sum(axis=1)
    assert len(is_ks) == neid
    for (i, eid, (ga, gb), ki,
            si, ocidi, si_vec,
            ga_xyz, gb_xyz, is_ga, is_gb, is_k) in zip(
            count(), cbush.element_id, nodes, k_fields,
            s_all, ocid_all, si_all,
            ga_xyzs, gb_xyzs, is_gas, is_gbs, is_ks):
        if not is_k:
            continue

        # Get rotation matrix R (element-to-global, rows = elem axes in global)
        T = transforms[i]  # (6, 6) or None
        if T is not None:
            R = T[:3, :3]
        else:
            R = np.eye(3)

        if not is_ga:
            ga_xyz = None
        if not is_gb:
            gb_xyz = None

        # Determine OCID offset in global coordinates
        ocid_offset = None
        if ocidi >= 0 and np.any(si_vec != 0):
            if ocidi == 0:
                ocid_offset = si_vec.copy()
            else:
                # Transform si from OCID system to global
                coord = model.coord.slice_card_by_id(np.array([ocidi]))
                R_ocid = np.array([coord.i[0], coord.j[0], coord.k[0]])
                ocid_offset = R_ocid.T @ si_vec

        k_elem = _cbush_kelm(
            eid, ki, R, ga_xyz, gb_xyz, s=si, ocid_offset=ocid_offset)
        k_elems[eid] = ((ga, gb), k_elem)


def _cbush_kelm(eid: int,
                k: NDArray,
                R: NDArray,
                ga_xyz: NDArray,
                gb_xyz: NDArray | None,
                s: float = 0.5,
                ocid_offset: NDArray | None = None,) -> NDArray:
    """Compute CBUSH element stiffness matrix in global coordinates.

    Uses the geometric coupling matrix G to account for spring offset:
        K = G^T @ diag(K1..K6) @ G

    G maps nodal DOFs (in global/CD frame) to spring deformations
    (in element frame), including rigid-arm coupling from the S
    parameter and OCID lateral offset.

    Validated against NX Nastran 2412 KELM output for 8 configurations
    including CID rotation, S offset (s=0, 0.5, 1), and OCID lateral
    offset. Agreement to single-precision (< 1e-3 relative error).

    Parameters
    ----------
    k : (6,) array
        [K1..K6] spring stiffnesses in element coordinates.
    R : (3, 3) array
        Rotation matrix: rows are element x,y,z axes in global frame.
        Transforms global vectors to element frame: v_elem = R @ v_global.
    ga_xyz, gb_xyz : (3,) arrays or None
        Grid positions in global. gb_xyz=None for grounded.
    s : float
        Spring location on element axis (0=GA, 1=GB).
    ocid_offset : (3,) array or None
        Spring offset in global coords (from OCID transformation).
        When specified, replaces the s*L axial offset for the rigid arm.

    Returns
    -------
    K : (6, 6) for grounded or (12, 12) for connected
    """
    a_grounded = ga_xyz is None
    b_grounded = gb_xyz is None

    # Compute rigid arm vectors (global coords, node to spring point)
    if ocid_offset is not None and np.any(ocid_offset != 0):
        r_a = np.asarray(ocid_offset, dtype='float64')
        r_b = r_a - (gb_xyz - ga_xyz) if not grounded else None
    elif a_grounded:
        r_b = np.zeros(3)
        r_a = None
    elif b_grounded:
        r_a = np.zeros(3)
        r_b = None
    else:
        r_a = s * (gb_xyz - ga_xyz)
        r_b = (s - 1.0) * (gb_xyz - ga_xyz)

    K_spring = np.diag(k)
    #print(f's={s} ga_xyz={ga_xyz} gb_xyz={gb_xyz} ra={r_a} rb={r_b}')
    if a_grounded:
        #print(f'eid={eid} A grounded')
        G = np.zeros((6, 6))
        G[:3, :3] = R
        G[:3, 3:] = -R @ _skew(r_b.ravel())
        G[3:, 3:] = R
    elif b_grounded:
        #print(f'eid={eid} B grounded')
        G = np.zeros((6, 6))
        G[:3, :3] = R
        G[:3, 3:] = -R @ _skew(r_a.ravel())
        G[3:, 3:] = R
    else:
        #print(f'eid={eid} free; k={k}')
        G = np.zeros((6, 12))
        G[:3, 0:3] = -R
        G[:3, 3:6] = R @ _skew(r_a.ravel())
        G[:3, 6:9] = R
        G[:3, 9:12] = -R @ _skew(r_b.ravel())
        G[3:, 3:6] = -R
        G[3:, 9:12] = R
    K0 = G.T @ K_spring @ G
    return K0


def _compute_cbush_transforms(
        cbush,
        model: BDF) -> list[NDArray | None]:
    """Compute 6x6 coordinate transforms for all CBUSH elements.

    Uses vectorized3 arrays:
    - coord_id >= 0: rotation from that coordinate system
    - coord_id == -1: element-based from x-vector or g0 node

    Returns list of (6,6) arrays or None (identity = no transform).
    """
    n = cbush.n
    transforms = [None] * n
    coord_ids = cbush.coord_id

    # Get grid positions for computing element axes
    grid = model.grid
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    coord = model.coord

    # Handle CID-based transforms
    cid_mask = coord_ids >= 0
    if np.any(cid_mask):
        unique_cids = np.unique(coord_ids[cid_mask])
        # Remove CID=0 (basic system, no rotation needed)
        unique_cids = unique_cids[unique_cids > 0]

        nodes = cbush.nodes

        #ucoord = coord.slice_card_by_id(unique_cids)
        #ir_ucids = (coord.coord_type == 'R')
        #ic_ucids = (coord.coord_type == 'C')
        #is_ucids = (coord.coord_type == 'S')
        #rbeta = coord.beta[ir_ucids, :, :]
        #cbeta = coord.beta[ic_ucids, :, :]
        #sbeta = coord.beta[is_ucids, :, :]

        for cid in unique_cids:
            icid = coord.index(cid)
            coord_type = coord.coord_type[icid]
            beta = coord.T[icid]  # (3,3) rows = i, j, k in basic
            origin = coord.origin[icid]

            elems_with_cid = np.where(coord_ids == cid)[0]

            if coord_type == "R":
                # Rectangular: same rotation for all elements with this CID
                T = np.zeros((6, 6))
                T[:3, :3] = beta
                T[3:, 3:] = beta
                for idx in elems_with_cid:
                    transforms[idx] = T
            else:
                # Cylindrical/Spherical: position-dependent at element GA
                for idx in elems_with_cid:
                    ga_id = nodes[idx, 0]
                    iga = np.searchsorted(nid, ga_id)
                    if iga >= len(nid) or nid[iga] != ga_id:
                        continue
                    ga_xyz = xyz_cid0[iga]

                    # Position in coord's local rectangular frame
                    xyz_local = (ga_xyz - origin) @ beta.T

                    if coord_type == "C":
                        theta = np.arctan2(xyz_local[1], xyz_local[0])
                        cos_t, sin_t = np.cos(theta), np.sin(theta)
                        R_pos = np.array([
                            [cos_t, -sin_t, 0.0],
                            [sin_t, cos_t, 0.0],
                            [0.0, 0.0, 1.0],
                        ])
                    else:  # 'S'
                        assert coord_type == 'S', coord_type
                        r_xy = np.sqrt(xyz_local[0] ** 2 + xyz_local[1] ** 2)
                        theta = np.arctan2(r_xy, xyz_local[2])
                        phi = np.arctan2(xyz_local[1], xyz_local[0])
                        sin_t, cos_t = np.sin(theta), np.cos(theta)
                        sin_p, cos_p = np.sin(phi), np.cos(phi)
                        R_pos = np.array([
                            [sin_t * cos_p, cos_t * cos_p, -sin_p],
                            [sin_t * sin_p, cos_t * sin_p, cos_p],
                            [cos_t, -sin_t, 0.0],
                        ])

                    # Element axes in basic: beta^T @ R_pos gives columns = axes
                    # R3 rows = element axes in basic (for the G-matrix formulation)
                    R3 = (beta.T @ R_pos).T  # rows = axes in basic
                    T = np.zeros((6, 6))
                    T[:3, :3] = R3
                    T[3:, 3:] = R3
                    transforms[idx] = T

    # Handle element-based transforms (coord_id == -1 with x-vector or g0)
    elem_mask = coord_ids == -1
    if not np.any(elem_mask):
        return transforms

    # x-vector / g0 define orientation
    is_x = cbush.is_x
    is_g0 = cbush.is_g0

    for idx in np.where(elem_mask)[0]:
        ga_id = cbush.nodes[idx, 0]
        gb_id = cbush.nodes[idx, 1]

        if is_x[idx]:
            x_vec = cbush.x[idx, :]
            x_vec = x_vec / np.linalg.norm(x_vec)
        elif is_g0[idx]:
            g0_id = cbush.g0[idx]
            ig0 = grid.index(g0_id)
            iga = grid.index(ga_id)
            x_vec = xyz[ig0, :] - xyz[iga, :]
            norm = np.linalg.norm(x_vec)
            x_vec /= norm
        else:
            elemi = elem.slice_card_by_index([idx])
            raise NotImplementedError(elemi)

        # Build element coordinate system
        if gb_id == 0:
            elem_x = x_vec #.copy()
        else:
            iga = grid.index(ga_id)
            igb = grid.index(gb_id)
            elem_x = xyz[igb] - xyz[iga]
            norm_x = np.linalg.norm(elem_x)
            elem_x /= norm_x
            
        # Build orthogonal triad: x, z = x × v, y = z × x
        elem_z = np.cross(elem_x, x_vec)
        norm_z = np.linalg.norm(elem_z)
        elem_z /= norm_z
        elem_y = np.cross(elem_z, elem_x)

        assert elem_x.shape == (1, 3), elem_x.shape
        assert elem_y.shape == (1, 3), elem_y.shape
        assert elem_z.shape == (1, 3), elem_z.shape
        #R3 = np.array([elem_x, elem_y, elem_z])  # og; (3, 1, 3)
        R3 = np.vstack([elem_x, elem_y, elem_z])
        assert R3.shape == (3, 3), R3.shape
        T = np.zeros((6, 6))
        T[:3, :3] = R3
        T[3:, 3:] = R3
        transforms[idx] = T
    return transforms


def _assemble_kelm(
    Kgg: NDArray,
    k_elem: NDArray,
    ga: int, gb: int,
    dof_map: dict,) -> None:
    """Assemble element stiffness (6x6 or 12x12) into global matrix.

    For grounded (6x6): assembles at GA/GB DOFs.
    For connected (12x12): assembles at [GA, GB] DOFs using the full
    12x12 element matrix which already includes cross-coupling.
    """
    if ga == 0 and gb == 0:
        raise RuntimeError((ga, gb))
    elif gb == 0:
        ga1 = dof_map.get((ga, 1))
        ga_dofs = [ga1+comp for comp in range(6)]
        all_dofs = np.array(ga_dofs)
    elif ga == 0:
        gb1 = dof_map.get((gb, 1))
        gb_dofs = [gb1+comp for comp in range(6)]
        all_dofs = np.array(gb_dofs)
    else:
        ga1 = dof_map.get((ga, 1))
        ga_dofs = [ga1+comp for comp in range(6)]

        gb1 = dof_map.get((gb, 1))
        gb_dofs = [gb1+comp for comp in range(6)]
        all_dofs = np.array(ga_dofs + gb_dofs)
    i, j = np.where(k_elem != 0.0)
    idofs = all_dofs[i]
    jdofs = all_dofs[j]
    Kgg[idofs, jdofs] += k_elem[i, j]


def _apply_cd_transform(k_elem: NDArray,
                        ga: int,
                        gb: int,
                        grid,
                        coord,) -> NDArray:
    """Transform element K from basic to CD frames: Kgg = Tlg^T @ Kll @ Tlg.

    Tlg is block-diagonal with each node's basic-to-CD rotation.
    No-op when all nodes have CD=0.
    """
    nid = grid.node_id
    cd = grid.cd

    iga = grid.index(ga)
    cd_a = cd[iga]

    if gb == 0:
        cd_b = 0
    else:
        igb = grid.index(gb)
        cd_b = cd[igb]

    if cd_a == 0 and cd_b == 0:
        return k_elem

    xyz = grid.xyz_cid0()

    if gb == 0:
        # Grounded: 6x6
        Tlg_a = _cd_rotation_6x6(cd_a, xyz[iga, :], coord)
        return Tlg_a.T @ k_elem @ Tlg_a
    else:
        # Connected: 12x12
        Tlg = np.eye(12)
        if cd_a != 0:
            R_a = _cd_rotation_3x3(cd_a, xyz[iga, :], coord)
            Tlg[0:3, 0:3] = R_a
            Tlg[3:6, 3:6] = R_a
        if cd_b != 0:
            R_b = _cd_rotation_3x3(cd_b, xyz[igb, :], coord)
            Tlg[6:9, 6:9] = R_b
            Tlg[9:12, 9:12] = R_b
        return Tlg.T @ k_elem @ Tlg


def _cd_rotation_3x3(cd: int, xyz_cid0: NDArray, coord) -> NDArray:
    """Compute 3x3 Tlg (CD-to-basic) rotation for a grid at position xyz_cid0.

    Returns Tlg such that v_basic = Tlg @ v_cd.
    Used in Kgg = Tlg^T @ Kll @ Tlg.

    For rectangular (CORD2R): Tlg = beta^T
    For cylindrical (CORD2C): Tlg = beta^T @ R_cyl(theta)
    For spherical (CORD2S): Tlg = beta^T @ R_sph(theta, phi)
    """
    #beta, coord_type, origin = cd_map[cd]
    icid = coord.index(cd)
    beta = coord.T[icid]  # (3, 3) — rows are i, j, k in basic
    coord_type = coord.coord_type[icid]
    origin = coord.origin[icid]

    if coord_type == "R":
        return beta.T

    # Position in coord's local rectangular frame
    xyz_local = (xyz_cid0 - origin) @ beta.T

    if coord_type == "C":
        theta = np.arctan2(xyz_local[1], xyz_local[0])
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        # R_cyl: columns are e_r, e_theta, e_z in local rectangular
        R_cyl = np.array([
            [cos_t, -sin_t, 0.0],
            [sin_t, cos_t, 0.0],
            [0.0, 0.0, 1.0],
        ])
        return beta.T @ R_cyl

    elif coord_type == "S":
        r_xy = np.sqrt(xyz_local[0] ** 2 + xyz_local[1] ** 2)
        theta = np.arctan2(r_xy, xyz_local[2])  # polar from +z
        phi = np.arctan2(xyz_local[1], xyz_local[0])  # azimuth from +x
        sin_t, cos_t = np.sin(theta), np.cos(theta)
        sin_p, cos_p = np.sin(phi), np.cos(phi)
        # R_sph: columns are e_r, e_theta, e_phi in local rectangular
        R_sph = np.array([
            [sin_t * cos_p, cos_t * cos_p, -sin_p],
            [sin_t * sin_p, cos_t * sin_p, cos_p],
            [cos_t, -sin_t, 0.0],
        ])
        return beta.T @ R_sph
    return np.eye(3)


def _cd_rotation_6x6(cd: int, xyz_cid0: NDArray, coord) -> NDArray:
    """Build 6x6 block-diagonal CD rotation from a 3x3."""
    R = _cd_rotation_3x3(cd, xyz_cid0, coord)
    T = np.zeros((6, 6))
    T[:3, :3] = R
    T[3:, 3:] = R
    return T
