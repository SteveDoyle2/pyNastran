"""Lumped mass matrix (MGG) assembly from a Nastran BDF (vectorized).

Assembles the global mass matrix in the g-set (all grid DOFs) using the
lumped mass approach (PARAM,COUPMASS,-1 in Nastran). Uses bdf_vectorized3's
array-based card storage for fast, fully-vectorized assembly with no Python
loops over elements.

Coordinate Systems
------------------
Nastran defines two coordinate frames per grid point:

    CP (card field 3): defines the input location of the grid point.
    CD (card field 7): defines the displacement coordinate system —
        i.e. the directions of T1, T2, T3, R1, R2, R3 at that grid.

The mass matrix lives in one of two frames:

    Mbb (basic frame):
        All DOFs expressed in the basic (CID=0) coordinate system.
        This is the frame used internally during assembly — lumped
        translational mass is isotropic so no rotation is needed for
        structural elements, and CONM2 offsets/inertia are rotated
        from their CID to basic.

    Mgg (global/grid frame):
        Each grid's 6 DOFs are in that grid's CD coordinate system.
        This is the standard Nastran output frame for system matrices.
        The transform is: Mgg = T^T @ Mbb @ T, where T is the sparse
        block-diagonal matrix of 6x6 rotation blocks (one per grid).
        For grids with CD=0, T_i = I_6 (no rotation).

    If all grids have CD=0 (common), Mbb = Mgg.

CONM2 CID Handling
~~~~~~~~~~~~~~~~~~
The CONM2 card's offset vector (X1, X2, X3) and inertia tensor are
defined in coordinate system CID (card field 4). Before computing the
6x6 mass matrix, the offset is rotated from CID to basic:

    d_basic = R_cid @ d_cid    (R_cid = 3x3 CID-to-basic rotation)

The inertia tensor is also rotated: I_basic = R @ I_cid @ R^T.

CID = 0: offset already in basic (no rotation needed).
CID = -1: offset is in basic (Nastran convention for absolute).

Lumped Mass Formulation
-----------------------
For an element with total mass m and n_nodes connected grid points:
    m_per_node = m / n_nodes

Each node receives m_per_node on diagonal entries for DOFs 1, 2, 3
(x, y, z translations). Rotational DOFs (4, 5, 6) get zero mass from
distributed elements (only CONM2 can add rotary inertia).

DOF Ordering
~~~~~~~~~~~~
The g-set has 6 DOFs per grid point: [T1, T2, T3, R1, R2, R3].
Grids are ordered by ascending grid ID. The global DOF index for
grid i (0-based), component c (1-based) is: 6 * i + (c - 1).

Output: scipy sparse CSR matrix of shape (n_dof, n_dof) where
n_dof = 6 * n_grids.

See Also
--------
pyNastran.dev.bdf_vectorized3.mesh_utils.inertia_relief :
    Inertia relief using the mass matrix from this module.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray
from scipy import sparse
from cpylog import SimpleLogger

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

log = SimpleLogger(level='warning')

# Element types that contribute no mass (springs, dampers, etc.)
NO_MASS = {
    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
    'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CVISC', 'CGAP',
    'CFAST', 'GENEL',
    'CONM1', 'CONM2',
    'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',
}


def build_mgg_lumped(
    bdf_filename: str | Path | None = None,
    model: 'BDF | None' = None,
    wtmass: float = 1.0,
    output_frame: str = 'basic',
) -> tuple[sparse.csr_matrix, NDArray[np.integer], float]:
    """Assemble the lumped global mass matrix from a BDF (vectorized).

    Parameters
    ----------
    bdf_filename : str, Path, or None
        Path to the Nastran BDF/DAT file. Ignored if model is provided.
    model : BDF or None
        Pre-loaded bdf_vectorized3 BDF object. If None, reads from
        bdf_filename.
    wtmass : float
        Mass unit conversion factor (Nastran PARAM,WTMASS). Default 1.0.
    output_frame : str
        Which coordinate frame for the output matrix:
        - 'basic' (default): Mbb — all DOFs in basic (CID=0) frame.
          CONM2 offsets/inertia are rotated from their CID to basic.
          No CD rotation applied. Use this for analyses that work in
          the basic frame or when all grids have CD=0.
        - 'global': Mgg — each grid's DOFs in its CD frame.
          Applies Mgg = T^T @ Mbb @ T where T is the block-diagonal
          CD rotation matrix. This matches Nastran's standard MGG output.
          For grids with CD=0, no rotation is applied (T_i = I_6).

    Returns
    -------
    M : scipy.sparse.csr_matrix of shape (n_dof, n_dof)
        Global mass matrix. n_dof = 6 * n_grids.
    grid_ids : (n_grids,) int array
        Sorted grid IDs defining DOF ordering.
    total_mass : float
        Total scalar mass (trace of translational DOFs / 3).
        Invariant under coordinate transformation.
    """
    if output_frame not in ('basic', 'global'):
        raise ValueError(
            f"output_frame must be 'basic' or 'global', got {output_frame!r}")

    if model is None:
        from pyNastran.dev.bdf_vectorized3.bdf import read_bdf
        model = read_bdf(str(bdf_filename), xref=True, debug=False)

    grid_ids = model.grid.node_id  # already sorted
    n_grids = len(grid_ids)
    n_dof = 6 * n_grids

    # Accumulator arrays for COO sparse assembly
    all_rows = []
    all_cols = []
    all_vals = []

    # -----------------------------------------------------------------
    # 1. Structural elements: lumped mass to nodes
    # -----------------------------------------------------------------
    _add_structural_elements(model, grid_ids, n_grids, wtmass,
                             all_rows, all_cols, all_vals)

    # -----------------------------------------------------------------
    # 2. CONM2: full 6x6 mass blocks (vectorized, CID-rotated)
    # -----------------------------------------------------------------
    _add_conm2(model, grid_ids, n_grids, wtmass,
               all_rows, all_cols, all_vals)

    # -----------------------------------------------------------------
    # 3. CONM1: direct 6x6 mass matrices
    # -----------------------------------------------------------------
    _add_conm1(model, grid_ids, n_grids, wtmass,
               all_rows, all_cols, all_vals)

    # -----------------------------------------------------------------
    # 4. CMASS2: scalar mass at DOF pairs
    # -----------------------------------------------------------------
    _add_cmass2(model, grid_ids, n_grids, wtmass,
                all_rows, all_cols, all_vals)

    # -----------------------------------------------------------------
    # 5. CMASS1: scalar mass via property lookup
    # -----------------------------------------------------------------
    _add_cmass1(model, grid_ids, n_grids, wtmass,
                all_rows, all_cols, all_vals)

    # -----------------------------------------------------------------
    # Assemble sparse Mbb
    # -----------------------------------------------------------------
    if len(all_rows) == 0:
        log.warning("No mass found in model")
        M = sparse.csr_matrix((n_dof, n_dof), dtype=np.float64)
        return M, grid_ids, 0.0

    rows_arr = np.concatenate(all_rows)
    cols_arr = np.concatenate(all_cols)
    vals_arr = np.concatenate(all_vals)

    M = sparse.coo_matrix(
        (vals_arr, (rows_arr, cols_arr)), shape=(n_dof, n_dof)
    ).tocsr()

    # -----------------------------------------------------------------
    # Transform Mbb -> Mgg if requested
    # -----------------------------------------------------------------
    if output_frame == 'global':
        M = _transform_mbb_to_mgg(M, model, grid_ids, n_grids)

    # Total mass = (sum of T1 + T2 + T3 diagonals) / 3
    diag = M.diagonal()
    total_mass = float(
        (diag[0::6].sum() + diag[1::6].sum() + diag[2::6].sum()) / 3.0
    )
    log.info(f"MGG assembled: {n_dof} DOFs, total mass = {total_mass:.6g}, frame={output_frame}")

    return M, grid_ids, total_mass


# =============================================================================
# COORDINATE TRANSFORMS
# =============================================================================

def _get_cid_rotation(model: 'BDF', cid: int) -> NDArray:
    """Get 3x3 rotation matrix from CID to basic (CID=0).

    Returns I_3 for CID=0 or CID=-1 (both are basic frame).
    """
    if cid == 0 or cid == -1:
        return np.eye(3, dtype=np.float64)
    return model.coord.xyz_to_global_transform[cid]


def _transform_mbb_to_mgg(
    Mbb: sparse.csr_matrix,
    model: 'BDF',
    grid_ids: NDArray,
    n_grids: int,
) -> sparse.csr_matrix:
    """Apply Mgg = T^T @ Mbb @ T for grids with non-zero CD.

    T is a sparse block-diagonal matrix of shape (n_dof, n_dof) where
    each 6x6 diagonal block is:
        T_i = [[R_i, 0  ],
               [0,   R_i]]
    with R_i the 3x3 rotation from basic to CD_i (= transpose of
    CD-to-basic transform, since R is orthogonal).

    For grids with CD=0, T_i = I_6 (no-op). If ALL grids have CD=0,
    returns Mbb unchanged (no allocation).
    """
    cd = model.grid.cd  # (n_grids,) int array of displacement coord IDs
    if np.all(cd == 0):
        return Mbb

    n_dof = 6 * n_grids

    # Vectorized T matrix construction — no Python loop over grids.
    # CD=0 grids get identity blocks; non-zero CID grids get R blocks in batch.
    t_rows_list = []
    t_cols_list = []
    t_vals_list = []

    # CD=0 grids: identity diagonal entries (6 per grid)
    is_cd0 = (cd == 0)
    n_cd0 = int(np.sum(is_cd0))
    if n_cd0 > 0:
        base_cd0 = 6 * np.where(is_cd0)[0]  # (n_cd0,)
        # 6 diagonal entries per grid
        dofs_cd0 = (base_cd0[:, np.newaxis] + np.arange(6)[np.newaxis, :]).ravel()
        t_rows_list.append(dofs_cd0)
        t_cols_list.append(dofs_cd0)
        t_vals_list.append(np.ones(len(dofs_cd0), dtype=np.float64))

    # Non-zero CID grids: scatter R blocks (two 3x3 per grid) in batch per CID
    unique_cds = np.unique(cd[~is_cd0]) if n_cd0 < n_grids else np.array([], dtype=cd.dtype)
    local_ij = np.mgrid[0:3, 0:3]  # (2, 3, 3): local_ij[0]=row offsets, [1]=col offsets
    local_i_flat = local_ij[0].ravel()  # [0,0,0,1,1,1,2,2,2]
    local_j_flat = local_ij[1].ravel()  # [0,1,2,0,1,2,0,1,2]

    for cid_val in unique_cds:
        cid_int = int(cid_val)
        if cid_int == 0:
            continue
        # R: basic-to-local (transpose of coord's to-basic rotation)
        R_to_basic = model.coord.xyz_to_global_transform[cid_int]
        R = R_to_basic.T
        R_flat = R.ravel()  # (9,)

        # Filter out zeros in R to reduce nnz
        nz_mask = R_flat != 0.0
        nz_i = local_i_flat[nz_mask]
        nz_j = local_j_flat[nz_mask]
        nz_v = R_flat[nz_mask]
        n_nz = int(np.sum(nz_mask))

        mask_cid = (cd == cid_val)
        grid_indices = np.where(mask_cid)[0]  # (n_cid,)
        n_cid = len(grid_indices)
        base_dofs = 6 * grid_indices  # (n_cid,)

        # Two 3x3 blocks per grid: offsets 0 and 3
        # rows: base + block_offset + nz_i, cols: base + block_offset + nz_j
        # Shape: (n_cid, 2*n_nz)
        block_offsets = np.array([0, 3], dtype=np.int64)
        # Expand: (n_cid, 2, n_nz)
        rows_block = (base_dofs[:, np.newaxis, np.newaxis] +
                      block_offsets[np.newaxis, :, np.newaxis] +
                      nz_i[np.newaxis, np.newaxis, :])
        cols_block = (base_dofs[:, np.newaxis, np.newaxis] +
                      block_offsets[np.newaxis, :, np.newaxis] +
                      nz_j[np.newaxis, np.newaxis, :])
        vals_block = np.broadcast_to(nz_v[np.newaxis, np.newaxis, :],
                                     (n_cid, 2, n_nz))

        t_rows_list.append(rows_block.ravel())
        t_cols_list.append(cols_block.ravel())
        t_vals_list.append(vals_block.ravel().copy())

    t_rows = np.concatenate(t_rows_list)
    t_cols = np.concatenate(t_cols_list)
    t_vals = np.concatenate(t_vals_list)

    T = sparse.coo_matrix(
        (t_vals, (t_rows, t_cols)), shape=(n_dof, n_dof)
    ).tocsr()

    # Mgg = T^T @ Mbb @ T
    Mgg = T.T @ Mbb @ T
    return Mgg


# =============================================================================
# MASS ASSEMBLY HELPERS
# =============================================================================

def _add_structural_elements(
    model: 'BDF',
    grid_ids: NDArray,
    n_grids: int,
    wtmass: float,
    all_rows: list,
    all_cols: list,
    all_vals: list,
) -> None:
    """Add lumped mass from all structural element cards (vectorized).

    For each card type, computes per-element mass, divides by number of
    nodes, then scatters to T1/T2/T3 diagonals using searchsorted indexing.

    Lumped translational mass (m*I_3) is invariant under rotation, so no
    coordinate transform is needed here regardless of output frame.
    """
    structural_cards = []
    for card in model.element_cards:
        if card.n == 0:
            continue
        if card.type in NO_MASS:
            continue
        structural_cards.append(card)

    for card in structural_cards:
        elem_mass = card.mass()  # (n_elem,)
        if not np.any(elem_mass > 0):
            continue

        nodes = card.nodes  # (n_elem, n_nodes_per_elem)
        if nodes.ndim == 1:
            nodes = nodes.reshape(-1, 1)

        n_elem, n_nodes_per_elem = nodes.shape

        mass_per_node = elem_mass / n_nodes_per_elem  # (n_elem,)

        mask = mass_per_node > 0
        if not np.any(mask):
            continue

        mass_per_node = mass_per_node[mask]
        nodes_masked = nodes[mask]

        node_ids_flat = nodes_masked.ravel()
        mass_flat = np.repeat(mass_per_node, n_nodes_per_elem)

        valid = node_ids_flat > 0
        if not np.all(valid):
            node_ids_flat = node_ids_flat[valid]
            mass_flat = mass_flat[valid]

        if len(node_ids_flat) == 0:
            continue

        grid_idx = np.searchsorted(grid_ids, node_ids_flat)
        base_dof = 6 * grid_idx

        dof_offsets = np.array([0, 1, 2], dtype=np.int64)
        rows = (base_dof[:, np.newaxis] + dof_offsets[np.newaxis, :]).ravel()
        cols = rows  # diagonal: row == col, COO doesn't mutate
        vals = np.repeat(mass_flat * wtmass, 3)

        all_rows.append(rows)
        all_cols.append(cols)
        all_vals.append(vals)


def _add_conm2(
    model: 'BDF',
    grid_ids: NDArray,
    n_grids: int,
    wtmass: float,
    all_rows: list,
    all_cols: list,
    all_vals: list,
) -> None:
    """Add CONM2 mass elements as full 6x6 blocks (vectorized).

    Computes the 6x6 mass matrix for all CONM2 elements simultaneously.
    The offset and inertia are rotated from CID to basic before assembly,
    producing Mbb entries.
    """
    conm2 = model.conm2
    if conm2.n == 0:
        return

    n = conm2.n
    mass = conm2._mass.copy()          # (n,)
    node_id = conm2.node_id.copy()     # (n,)
    offset = conm2.xyz_offset.copy()   # (n, 3)
    inertia = conm2.inertia.copy()     # (n, 6): I11, I21, I22, I31, I32, I33
    coord_id = conm2.coord_id          # (n,) CID for offset/inertia

    # Skip zero-mass entries
    nonzero = mass != 0.0
    if not np.any(nonzero):
        return
    mass = mass[nonzero]
    node_id = node_id[nonzero]
    offset = offset[nonzero]
    inertia = inertia[nonzero]
    coord_id = coord_id[nonzero]
    n = len(mass)

    # -----------------------------------------------------------------
    # Rotate offset and inertia from CID to basic
    # -----------------------------------------------------------------
    unique_cids = np.unique(coord_id)
    needs_rotation = not (len(unique_cids) == 1 and int(unique_cids[0]) in (0, -1))

    if needs_rotation:
        for cid_val in unique_cids:
            cid_int = int(cid_val)
            if cid_int == 0 or cid_int == -1:
                continue
            mask_cid = coord_id == cid_val
            R = _get_cid_rotation(model, cid_int)  # 3x3 CID-to-basic

            # Rotate offset: d_basic = d_cid @ R (row-vector convention)
            offset[mask_cid] = offset[mask_cid] @ R

            # Rotate inertia tensor: I_basic = R @ I_cid @ R^T
            # Unpack symmetric inertia: I11, I21, I22, I31, I32, I33
            inertia_masked = inertia[mask_cid]
            n_masked = inertia_masked.shape[0]
            I_mat = np.zeros((n_masked, 3, 3), dtype=np.float64)
            I_mat[:, 0, 0] = inertia_masked[:, 0]  # I11
            I_mat[:, 1, 0] = inertia_masked[:, 1]  # I21
            I_mat[:, 0, 1] = inertia_masked[:, 1]  # I21 (symmetric)
            I_mat[:, 1, 1] = inertia_masked[:, 2]  # I22
            I_mat[:, 2, 0] = inertia_masked[:, 3]  # I31
            I_mat[:, 0, 2] = inertia_masked[:, 3]  # I31 (symmetric)
            I_mat[:, 2, 1] = inertia_masked[:, 4]  # I32
            I_mat[:, 1, 2] = inertia_masked[:, 4]  # I32 (symmetric)
            I_mat[:, 2, 2] = inertia_masked[:, 5]  # I33

            # I_basic = R^T @ I_cid @ R (row-vector convention: v_basic = v_local @ R)
            # (n_masked, 3, 3) = (3,3)^T @ (n_masked, 3, 3) @ (3,3)
            I_rotated = np.einsum('ji,njk,kl->nil', R, I_mat, R)

            # Repack to (n_masked, 6): I11, I21, I22, I31, I32, I33
            inertia[mask_cid, 0] = I_rotated[:, 0, 0]
            inertia[mask_cid, 1] = I_rotated[:, 1, 0]
            inertia[mask_cid, 2] = I_rotated[:, 1, 1]
            inertia[mask_cid, 3] = I_rotated[:, 2, 0]
            inertia[mask_cid, 4] = I_rotated[:, 2, 1]
            inertia[mask_cid, 5] = I_rotated[:, 2, 2]

    # -----------------------------------------------------------------
    # Build all 6x6 matrices in basic frame
    # -----------------------------------------------------------------
    M_all = np.zeros((n, 6, 6), dtype=np.float64)

    # Translational mass
    M_all[:, 0, 0] = mass
    M_all[:, 1, 1] = mass
    M_all[:, 2, 2] = mass

    # Offset coupling
    x1 = offset[:, 0]
    x2 = offset[:, 1]
    x3 = offset[:, 2]

    # m * S^T (translational rows, rotational columns)
    M_all[:, 0, 4] = mass * x3
    M_all[:, 0, 5] = -mass * x2
    M_all[:, 1, 3] = -mass * x3
    M_all[:, 1, 5] = mass * x1
    M_all[:, 2, 3] = mass * x2
    M_all[:, 2, 4] = -mass * x1

    # m * S (rotational rows, translational columns)
    M_all[:, 3, 1] = -mass * x3
    M_all[:, 3, 2] = mass * x2
    M_all[:, 4, 0] = mass * x3
    M_all[:, 4, 2] = -mass * x1
    M_all[:, 5, 0] = -mass * x2
    M_all[:, 5, 1] = mass * x1

    # Parallel axis theorem: m * S^T @ S
    M_all[:, 3, 3] += mass * (x2**2 + x3**2)
    M_all[:, 4, 4] += mass * (x1**2 + x3**2)
    M_all[:, 5, 5] += mass * (x1**2 + x2**2)
    M_all[:, 3, 4] += mass * (-x1 * x2)
    M_all[:, 4, 3] += mass * (-x1 * x2)
    M_all[:, 3, 5] += mass * (-x1 * x3)
    M_all[:, 5, 3] += mass * (-x1 * x3)
    M_all[:, 4, 5] += mass * (-x2 * x3)
    M_all[:, 5, 4] += mass * (-x2 * x3)

    # Rotary inertia (already in basic after rotation above)
    I11 = inertia[:, 0]
    I21 = inertia[:, 1]
    I22 = inertia[:, 2]
    I31 = inertia[:, 3]
    I32 = inertia[:, 4]
    I33 = inertia[:, 5]

    M_all[:, 3, 3] += I11
    M_all[:, 4, 3] += I21
    M_all[:, 3, 4] += I21
    M_all[:, 4, 4] += I22
    M_all[:, 5, 3] += I31
    M_all[:, 3, 5] += I31
    M_all[:, 5, 4] += I32
    M_all[:, 4, 5] += I32
    M_all[:, 5, 5] += I33

    # -----------------------------------------------------------------
    # Scatter to COO
    # -----------------------------------------------------------------
    grid_idx = np.searchsorted(grid_ids, node_id)
    base_dof = 6 * grid_idx

    local_i, local_j = np.mgrid[0:6, 0:6]
    local_i_flat = local_i.ravel()  # (36,)
    local_j_flat = local_j.ravel()

    rows = (base_dof[:, np.newaxis] + local_i_flat[np.newaxis, :]).ravel()
    cols = (base_dof[:, np.newaxis] + local_j_flat[np.newaxis, :]).ravel()
    vals = (M_all.reshape(n, 36) * wtmass).ravel()

    nonzero_vals = vals != 0.0
    if np.any(nonzero_vals):
        all_rows.append(rows[nonzero_vals])
        all_cols.append(cols[nonzero_vals])
        all_vals.append(vals[nonzero_vals])


def _add_conm1(
    model: 'BDF',
    grid_ids: NDArray,
    n_grids: int,
    wtmass: float,
    all_rows: list,
    all_cols: list,
    all_vals: list,
) -> None:
    """Add CONM1 elements (direct 6x6 mass matrices).

    CONM1 provides its mass matrix directly — no CID rotation is applied
    here (Nastran convention: CONM1 is defined in the grid's displacement
    coordinate system). The Mbb<->Mgg transform handles this via the
    block-diagonal T matrix.
    """
    conm1 = model.conm1
    if conm1.n == 0:
        return

    n = conm1.n
    node_id = conm1.node_id
    M_matrices = conm1._mass  # (n, 6, 6)

    grid_idx = np.searchsorted(grid_ids, node_id)
    base_dof = 6 * grid_idx

    local_i, local_j = np.mgrid[0:6, 0:6]
    local_i_flat = local_i.ravel()
    local_j_flat = local_j.ravel()

    rows = (base_dof[:, np.newaxis] + local_i_flat[np.newaxis, :]).ravel()
    cols = (base_dof[:, np.newaxis] + local_j_flat[np.newaxis, :]).ravel()
    vals = (M_matrices.reshape(n, 36) * wtmass).ravel()

    nonzero_vals = vals != 0.0
    if np.any(nonzero_vals):
        all_rows.append(rows[nonzero_vals])
        all_cols.append(cols[nonzero_vals])
        all_vals.append(vals[nonzero_vals])


def _add_cmass2(
    model: 'BDF',
    grid_ids: NDArray,
    n_grids: int,
    wtmass: float,
    all_rows: list,
    all_cols: list,
    all_vals: list,
) -> None:
    """Add CMASS2 scalar mass elements (vectorized)."""
    cmass2 = model.cmass2
    if cmass2.n == 0:
        return

    mass = cmass2._mass
    nodes = cmass2.nodes
    components = cmass2.components

    nonzero = mass != 0.0
    if not np.any(nonzero):
        return
    mass = mass[nonzero]
    nodes = nodes[nonzero]
    components = components[nonzero]

    n1 = nodes[:, 0]
    c1 = components[:, 0]
    n2 = nodes[:, 1]
    c2 = components[:, 1]

    rows_list = []
    cols_list = []
    vals_list = []

    valid1 = (n1 > 0) & (c1 > 0)
    if np.any(valid1):
        idx1 = np.searchsorted(grid_ids, n1[valid1])
        dof1 = 6 * idx1 + (c1[valid1] - 1)
        rows_list.append(dof1)
        cols_list.append(dof1)
        vals_list.append(mass[valid1] * wtmass)

    valid2 = (n2 > 0) & (c2 > 0)
    if np.any(valid2):
        idx2 = np.searchsorted(grid_ids, n2[valid2])
        dof2 = 6 * idx2 + (c2[valid2] - 1)
        rows_list.append(dof2)
        cols_list.append(dof2)
        vals_list.append(mass[valid2] * wtmass)

        both_valid = valid1 & valid2
        if np.any(both_valid):
            idx1b = np.searchsorted(grid_ids, n1[both_valid])
            dof1b = 6 * idx1b + (c1[both_valid] - 1)
            idx2b = np.searchsorted(grid_ids, n2[both_valid])
            dof2b = 6 * idx2b + (c2[both_valid] - 1)
            m_both = mass[both_valid] * wtmass
            rows_list.append(dof1b)
            cols_list.append(dof2b)
            vals_list.append(-m_both)
            rows_list.append(dof2b)
            cols_list.append(dof1b)
            vals_list.append(-m_both)

    if rows_list:
        all_rows.append(np.concatenate(rows_list))
        all_cols.append(np.concatenate(cols_list))
        all_vals.append(np.concatenate(vals_list))


def _add_cmass1(
    model: 'BDF',
    grid_ids: NDArray,
    n_grids: int,
    wtmass: float,
    all_rows: list,
    all_cols: list,
    all_vals: list,
) -> None:
    """Add CMASS1 scalar mass elements (mass from PMASS property)."""
    cmass1 = model.cmass1
    if cmass1.n == 0:
        return

    mass = cmass1.mass()
    nodes = cmass1.nodes
    components = cmass1.components

    nonzero = mass != 0.0
    if not np.any(nonzero):
        return
    mass = mass[nonzero]
    nodes = nodes[nonzero]
    components = components[nonzero]

    n1 = nodes[:, 0]
    c1 = components[:, 0]
    n2 = nodes[:, 1]
    c2 = components[:, 1]

    rows_list = []
    cols_list = []
    vals_list = []

    valid1 = (n1 > 0) & (c1 > 0)
    if np.any(valid1):
        idx1 = np.searchsorted(grid_ids, n1[valid1])
        dof1 = 6 * idx1 + (c1[valid1] - 1)
        rows_list.append(dof1)
        cols_list.append(dof1)
        vals_list.append(mass[valid1] * wtmass)

    valid2 = (n2 > 0) & (c2 > 0)
    if np.any(valid2):
        idx2 = np.searchsorted(grid_ids, n2[valid2])
        dof2 = 6 * idx2 + (c2[valid2] - 1)
        rows_list.append(dof2)
        cols_list.append(dof2)
        vals_list.append(mass[valid2] * wtmass)

        both_valid = valid1 & valid2
        if np.any(both_valid):
            idx1b = np.searchsorted(grid_ids, n1[both_valid])
            dof1b = 6 * idx1b + (c1[both_valid] - 1)
            idx2b = np.searchsorted(grid_ids, n2[both_valid])
            dof2b = 6 * idx2b + (c2[both_valid] - 1)
            m_both = mass[both_valid] * wtmass
            rows_list.append(dof1b)
            cols_list.append(dof2b)
            vals_list.append(-m_both)
            rows_list.append(dof2b)
            cols_list.append(dof1b)
            vals_list.append(-m_both)

    if rows_list:
        all_rows.append(np.concatenate(rows_list))
        all_cols.append(np.concatenate(cols_list))
        all_vals.append(np.concatenate(vals_list))


# =============================================================================
# GRID POINT WEIGHT
# =============================================================================

def get_grid_point_weight(
    M: sparse.spmatrix,
    grid_ids: NDArray[np.integer],
    node_xyz: NDArray[np.floating],
    ref_point: NDArray[np.floating] | None = None,
) -> NDArray[np.floating]:
    """Compute the 6x6 rigid body mass matrix (grid point weight).

    Equivalent to Nastran's OWEIGHT output. Uses D^T @ M @ D where D
    is the rigid body mode matrix about the reference point.

    The input mass matrix M should be in basic frame (Mbb) for this
    computation to produce correct results.

    Builds D vectorized (no Python loop over grids).

    Parameters
    ----------
    M : sparse matrix, shape (n_dof, n_dof)
        Global mass matrix (should be Mbb for correct OWEIGHT).
    grid_ids : (n_grids,) int array
        Grid IDs in M's DOF ordering.
    node_xyz : (n_grids, 3) float array
        Grid positions in basic coordinate system.
    ref_point : (3,) array or None
        Reference point for mass properties. If None, uses origin [0,0,0].

    Returns
    -------
    M0 : (6, 6) ndarray
        Rigid body mass matrix (grid point weight table).
    """
    if ref_point is None:
        ref_point = np.zeros(3)
    ref_point = np.asarray(ref_point, dtype=np.float64)

    n_grids = len(grid_ids)
    n_dof = 6 * n_grids

    # Arm vectors for all grids
    arms = node_xyz - ref_point[np.newaxis, :]  # (n_grids, 3)
    dx = arms[:, 0]
    dy = arms[:, 1]
    dz = arms[:, 2]

    # Build D matrix (n_dof, 6) vectorized
    D = np.zeros((n_dof, 6), dtype=np.float64)
    base = np.arange(n_grids, dtype=np.int64) * 6

    # Translation modes
    D[base, 0] = 1.0
    D[base + 1, 1] = 1.0
    D[base + 2, 2] = 1.0

    # Rotation modes: cross product contributions
    D[base + 1, 3] = -dz
    D[base + 2, 3] = dy
    D[base, 4] = dz
    D[base + 2, 4] = -dx
    D[base, 5] = -dy
    D[base + 1, 5] = dx

    # Rotational DOFs: identity
    D[base + 3, 3] = 1.0
    D[base + 4, 4] = 1.0
    D[base + 5, 5] = 1.0

    # M0 = D^T @ M @ D
    if sparse.issparse(M):
        MD = M.dot(D)
    else:
        MD = M @ D
    M0 = D.T @ MD

    return M0
