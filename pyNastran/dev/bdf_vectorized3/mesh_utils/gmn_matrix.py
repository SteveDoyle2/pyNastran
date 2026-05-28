"""GMN matrix assembly from RBE3 elements.

The GMN matrix maps independent (n-set) DOFs to dependent (m-set) DOFs:
    u_m = GMN @ u_n

For RBE3, this is a weighted least-squares interpolation:
    - Build rigid-body modes from each independent node relative to the
      dependent (reference) node
    - Apply weighting to the RB modes
    - Solve: GMN_row = (RB^T @ W @ RB)^-1 @ RB^T @ W

In Nastran DMAP, the flow is:
    1. RBE3 cards -> RMG constraint equation matrix
    2. XMCE1 subdmap partitions RMG: GMN = -inv(RMM) * RMN
    3. u_m = GMN @ u_n in displacement set reduction

This module implements both:
    - Per-element GMN rows (for a single RBE3)
    - Global GMN assembly (all RBE3 elements into one sparse matrix)
"""
import numpy as np
from numpy.typing import NDArray
from scipy import sparse


def _expand_dof_string(comp: int | str) -> list[int]:
    """Expand DOF component integer/string to list of 0-based indices.

    E.g. 123 -> [0, 1, 2], 456 -> [3, 4, 5], 123456 -> [0,1,2,3,4,5]
    """
    s = str(int(comp))
    return [int(c) - 1 for c in s]


def _rigid_body_matrix(dx: NDArray) -> NDArray:
    """6x6 rigid-body transformation from dependent to independent node.

    Given offset dx = x_ind - x_dep, the rigid-body relationship is:
        u_ind = M @ u_dep

    where M maps dependent-node motion to independent-node motion.
    Rows = independent DOFs [T1,T2,T3,R1,R2,R3]
    Columns = dependent DOFs [T1,T2,T3,R1,R2,R3]
    """
    x, y, z = dx
    return np.array([
        [1, 0, 0,  0,  z, -y],
        [0, 1, 0, -z,  0,  x],
        [0, 0, 1,  y, -x,  0],
        [0, 0, 0,  1,  0,  0],
        [0, 0, 0,  0,  1,  0],
        [0, 0, 0,  0,  0,  1],
    ], dtype=float)


def compute_rbe3_thermal_load(
    gmn_rows: NDArray,
    xyz_dep: NDArray,
    dep_dofs: list[int],
    ind_nodes_xyz: NDArray,
    ind_dofs_per_node: list[list[int]],
    alpha: float,
    tref: float,
    temperature: float,
) -> NDArray:
    """Compute the thermal load vector contribution from RBE3 ALPHA/TREF.

    In Nastran, the RBE3 ALPHA field does NOT modify the GMN interpolation.
    The dependent displacement is always u_dep = GMN @ u_ind regardless of
    temperature. Instead, alpha/tref produce an equivalent thermal load
    vector that gets applied through the multipoint constraint equations.

    The thermal expansion of the RBE3 offsets generates a "free thermal
    displacement" vector at the independent DOFs:
        u_thermal_k[T1,T2,T3] = alpha * (T - tref) * (x_ind_k - x_dep)
        u_thermal_k[R1,R2,R3] = 0

    This is mapped to the dependent DOFs via GMN:
        u_thermal_dep = GMN @ u_thermal

    The corresponding thermal force is applied as:
        F_thermal_dep = K_dep * u_thermal_dep

    Parameters
    ----------
    gmn_rows : (n_dep_dofs, n_ind_dofs_total) array
        The GMN matrix rows (from compute_rbe3_gmn_row).
    xyz_dep : (3,) array
        Dependent node position.
    dep_dofs : list of int (0-based)
        Dependent DOF indices.
    ind_nodes_xyz : (n_ind, 3) array
        Independent node positions.
    ind_dofs_per_node : list of list[int]
        Active DOFs (0-based) for each independent node.
    alpha : float
        Thermal expansion coefficient of the RBE3 connection.
    tref : float
        Reference temperature (stress-free).
    temperature : float
        Current temperature.

    Returns
    -------
    u_thermal_dep : (n_dep_dofs,) array
        Free thermal displacement at dependent DOFs. Multiply by stiffness
        to get the thermal load contribution.
    """
    dT = temperature - tref
    if abs(alpha) < 1e-30 or abs(dT) < 1e-30:
        return np.zeros(len(dep_dofs))

    # Build thermal expansion vector for independent DOFs
    n_ind_dofs = gmn_rows.shape[1]
    u_thermal = np.zeros(n_ind_dofs)

    col = 0
    for i_node, dofs in enumerate(ind_dofs_per_node):
        dx = ind_nodes_xyz[i_node] - xyz_dep
        for dof in dofs:
            if dof < 3:  # translational DOFs only
                u_thermal[col] = alpha * dT * dx[dof]
            col += 1

    return gmn_rows @ u_thermal


def compute_rbe3_gmn_row(
    xyz_dep: NDArray,
    dep_dofs: list[int],
    ind_nodes_xyz: NDArray,
    ind_dofs_per_node: list[list[int]],
    weights_per_node: NDArray,
    um_nodes_xyz: NDArray | None = None,
    um_dofs_per_node: list[list[int]] | None = None,
) -> NDArray:
    """Compute GMN rows for a single RBE3 element.

    Parameters
    ----------
    xyz_dep : (3,) array
        Position of dependent (reference) node in basic/global coords.
    dep_dofs : list of int (0-based)
        Dependent DOF indices for refgrid, e.g. [0,1,2,3,4,5] for 123456.
    ind_nodes_xyz : (n_ind, 3) array
        Positions of independent nodes in basic/global coords.
    ind_dofs_per_node : list of list[int]
        For each independent node, which DOFs (0-based) are active.
    weights_per_node : (n_ind,) array
        Weight factor for each independent node.
    um_nodes_xyz : (n_um, 3) array or None
        Positions of UM nodes. These are independent nodes whose specified
        DOFs are moved to the dependent (m-set).
    um_dofs_per_node : list of list[int] or None
        For each UM node, which DOFs (0-based) become dependent.

    Returns
    -------
    gmn_rows : (n_dep_dofs_total, n_ind_dofs_remaining) array
        The interpolation matrix rows for this RBE3.
        Rows include refgrid dep_dofs followed by UM DOFs.
        Columns are the remaining independent DOFs (after removing UM DOFs).
    ind_dof_indices : list of (node_idx, dof_idx) tuples
        Mapping from column index to (independent node index, dof).
        Node indices refer to the original ind_nodes_xyz array.
    """
    n_ind_nodes = len(ind_nodes_xyz)

    # Build the full column DOF index map (all independent DOFs)
    all_ind_dof_list = []  # (node_idx, dof_idx) for each column
    for i_node in range(n_ind_nodes):
        for dof in ind_dofs_per_node[i_node]:
            all_ind_dof_list.append((i_node, dof))
    n_all_ind_dofs = len(all_ind_dof_list)

    # Characteristic length for rotational DOF weighting
    deltas = ind_nodes_xyz - xyz_dep[None, :]
    dists = np.linalg.norm(deltas, axis=1)
    Lc = np.mean(dists) if n_ind_nodes > 0 else 1.0
    if Lc < 1e-12:
        Lc = 1.0

    # Build RB matrix and weights for ALL independent DOFs
    # rb shape: (n_all_ind_dofs, 6) — full 6 dep DOFs as columns
    rb = np.zeros((n_all_ind_dofs, 6))
    W = np.zeros(n_all_ind_dofs)

    row = 0
    for i_node in range(n_ind_nodes):
        dx = ind_nodes_xyz[i_node] - xyz_dep
        rb_full = _rigid_body_matrix(dx)  # (6, 6): rows=indep, cols=dep
        w = weights_per_node[i_node]

        for dof in ind_dofs_per_node[i_node]:
            rb[row, :] = rb_full[dof, :]  # full 6 dep columns
            if dof >= 3:
                W[row] = w * Lc * Lc
            else:
                W[row] = w
            row += 1

    # No UM: standard path
    if um_nodes_xyz is None or um_dofs_per_node is None or len(um_nodes_xyz) == 0:
        # Weighted least squares on full 6-DOF system
        rbw = rb.T * W       # (6, n_ind)
        A = rbw @ rb         # (6, 6) normal equation matrix

        rcond = 1e-12
        if np.linalg.matrix_rank(A, tol=rcond * np.max(np.abs(A))) < 6:
            gmn_full, _, _, _ = np.linalg.lstsq(A, rbw, rcond=rcond)
        else:
            gmn_full = np.linalg.solve(A, rbw)  # (6, n_ind)

        gmn_rows = gmn_full[dep_dofs, :]
        return gmn_rows, all_ind_dof_list

    # UM path: UM does NOT change how u_ref is computed from u_ind.
    # The refgrid weighted LS uses ALL listed independent DOFs (including UM ones).
    # UM only changes the DOF set partitioning for the solver — the UM DOFs
    # become m-set (dependent), meaning their values are solved via the coupled
    # constraint + stiffness system rather than independent force balance.
    #
    # For GMN purposes:
    #   - Refgrid rows: SAME as without UM (weighted LS of ALL ind DOFs)
    #   - UM DOFs: cannot be predicted from a simple GMN formula; they depend
    #     on the full system (stiffness + constraint coupling)
    #
    # We return the refgrid GMN rows and mark which columns are UM (m-set).
    # The UM DOF values are determined by the solver, not by GMN alone.

    # Standard weighted LS (same as no-UM case)
    rbw = rb.T * W
    A = rbw @ rb

    rcond = 1e-12
    if np.linalg.matrix_rank(A, tol=rcond * np.max(np.abs(A))) < 6:
        gmn_full, _, _, _ = np.linalg.lstsq(A, rbw, rcond=rcond)
    else:
        gmn_full = np.linalg.solve(A, rbw)  # (6, n_all_ind)

    # Refgrid rows (same as no-UM)
    gmn_rows = gmn_full[dep_dofs, :]

    # Identify which columns are UM DOFs (for caller's information)
    um_node_dof_set = set()
    for i_um, um_xyz in enumerate(um_nodes_xyz):
        for i_node in range(n_ind_nodes):
            if np.allclose(ind_nodes_xyz[i_node], um_xyz, atol=1e-10):
                for dof in um_dofs_per_node[i_um]:
                    um_node_dof_set.add((i_node, dof))
                break

    return gmn_rows, all_ind_dof_list


def compute_rbe3_constraint_matrices(
    xyz_dep: NDArray,
    dep_dofs: list[int],
    ind_nodes_xyz: NDArray,
    ind_dofs_per_node: list[list[int]],
    weights_per_node: NDArray,
) -> tuple[NDArray, NDArray, NDArray]:
    """Compute the RMG constraint matrix and its RMM/RMN partitions.

    The multipoint constraint equation is:
        RMM @ u_m + RMN @ u_n = 0
        => u_m = -inv(RMM) @ RMN @ u_n = GMN @ u_n

    For RBE3 without UM, RMM = I (identity) and RMN = -GMN_full[dep_dofs].

    Parameters
    ----------
    xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights_per_node
        Same as compute_rbe3_gmn_row.

    Returns
    -------
    RMM : (n_dep, n_dep) array
        m-set x m-set block. For standard RBE3, this is the identity.
    RMN : (n_dep, n_ind_total) array
        m-set x n-set block. GMN = -inv(RMM) @ RMN.
    GMN : (n_dep, n_ind_total) array
        The interpolation matrix (same as compute_rbe3_gmn_row output).
    """
    gmn_rows, ind_dof_list = compute_rbe3_gmn_row(
        xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights_per_node,
    )
    n_dep = len(dep_dofs)
    n_ind = gmn_rows.shape[1]

    # For standard RBE3: constraint is u_m - GMN @ u_n = 0
    # So RMM = I, RMN = -GMN
    RMM = np.eye(n_dep)
    RMN = -gmn_rows

    return RMM, RMN, gmn_rows


def assemble_gmn(model, dof_map: dict[tuple[int, int], int] | None = None,
                 ndof: int | None = None) -> tuple[NDArray, dict, dict]:
    """Assemble global GMN matrix from all RBE3 elements in the model.

    Parameters
    ----------
    model : pyNastran BDF (vectorized3)
        Must have RBE3 elements and grid coordinates.
    dof_map : dict mapping (node_id, dof_1based) -> global_index, optional
        If None, builds one from all grids (6 DOF per grid).
    ndof : int, optional
        Total DOFs. If None, determined from dof_map.

    Returns
    -------
    GMN : (n_m, n_n) sparse CSR matrix
        Rows = m-set (dependent) DOFs, columns = all DOFs.
        Only non-zero columns correspond to independent DOFs.
    m_set : dict mapping (node_id, dof_1based) -> m_set_index
        The dependent DOF set.
    dof_map : dict mapping (node_id, dof_1based) -> global_index
        The full DOF map used.
    """
    rbe3 = model.rbe3
    if rbe3.n == 0:
        raise ValueError("Model has no RBE3 elements")

    grid = model.grid
    xyz_cid0 = grid.xyz_cid0()
    node_ids = grid.node_id

    # Build default DOF map if not provided
    if dof_map is None:
        dof_map = {}
        idx = 0
        for nid in node_ids:
            for dof in range(1, 7):
                dof_map[(int(nid), dof)] = idx
                idx += 1
        ndof = idx
    elif ndof is None:
        ndof = max(dof_map.values()) + 1

    # Identify m-set DOFs (dependent DOFs from all RBE3s)
    m_set = {}  # (nid, dof) -> m_set_row_index
    m_row = 0

    # Pre-compute grid index lookup
    nid_to_idx = {int(nid): i for i, nid in enumerate(node_ids)}

    # Collect all GMN contributions
    rows_list = []
    cols_list = []
    vals_list = []

    iweight_all = rbe3.iweight
    igrid_start = 0  # running index into independent_nodes

    for i_elem in range(rbe3.n):
        ref_nid = int(rbe3.ref_grid[i_elem])
        ref_comp = int(rbe3.ref_component[i_elem])
        dep_dofs = _expand_dof_string(ref_comp)

        # Dependent node position
        xyz_dep = xyz_cid0[nid_to_idx[ref_nid]]

        # Register m-set DOFs
        m_rows_this = []
        for dof in dep_dofs:
            key = (ref_nid, dof + 1)
            if key not in m_set:
                m_set[key] = m_row
                m_row += 1
            m_rows_this.append(m_set[key])

        # Gather independent nodes, DOFs, weights for this element
        iw0, iw1 = iweight_all[i_elem]
        nweights = iw1 - iw0

        ind_nodes_xyz_list = []
        ind_dofs_per_node_list = []
        weights_per_node_list = []
        ind_global_dof_map = []  # (node_id, dof_1based) for each column

        for iw in range(nweights):
            weight = float(rbe3.weight[iw0 + iw])
            comp = int(rbe3.independent_dofs[iw0 + iw])
            active_dofs = _expand_dof_string(comp)

            ngrid = int(rbe3.ngrid_per_weight[igrid_start + iw])
            for ig in range(ngrid):
                flat_idx = sum(rbe3.ngrid_per_weight[igrid_start:igrid_start + iw]) + ig
                nid = int(rbe3.independent_nodes[
                    sum(rbe3.ngrid_per_weight[:igrid_start]) + flat_idx
                ])
                xyz_ind = xyz_cid0[nid_to_idx[nid]]
                ind_nodes_xyz_list.append(xyz_ind)
                ind_dofs_per_node_list.append(active_dofs)
                weights_per_node_list.append(weight)
                for dof in active_dofs:
                    ind_global_dof_map.append((nid, dof + 1))

        ind_nodes_xyz = np.array(ind_nodes_xyz_list)
        weights_per_node = np.array(weights_per_node_list)

        # Compute GMN rows for this element
        gmn_rows, ind_dof_list = compute_rbe3_gmn_row(
            xyz_dep, dep_dofs, ind_nodes_xyz,
            ind_dofs_per_node_list, weights_per_node,
        )

        # Place into global matrix
        for i_row, m_row_idx in enumerate(m_rows_this):
            for i_col, (nid_col, dof_col) in enumerate(ind_global_dof_map):
                val = gmn_rows[i_row, i_col]
                if abs(val) > 1e-15:
                    col_idx = dof_map.get((nid_col, dof_col))
                    if col_idx is not None:
                        rows_list.append(m_row_idx)
                        cols_list.append(col_idx)
                        vals_list.append(val)

        # Advance igrid_start
        igrid_start += nweights

    # Build sparse GMN
    n_m = m_row
    GMN = sparse.csr_matrix(
        (vals_list, (rows_list, cols_list)),
        shape=(n_m, ndof),
    )
    return GMN, m_set, dof_map


def assemble_gmn_simple(
    xyz_dep: NDArray,
    dep_comp: int | str,
    ind_nodes_xyz: NDArray,
    ind_comp: int | str,
    weights: NDArray | float = 1.0,
    alpha: float = 0.0,
    tref: float = 0.0,
    temperature: float | None = None,
) -> NDArray | tuple[NDArray, NDArray]:
    """Simple interface: single RBE3, uniform DOFs on all independent nodes.

    Parameters
    ----------
    xyz_dep : (3,) position of dependent node
    dep_comp : DOF component (e.g. 123456)
    ind_nodes_xyz : (n, 3) positions of independent nodes
    ind_comp : DOF component for all independent nodes (e.g. 123)
    weights : scalar or (n,) weights for independent nodes
    alpha : float, default 0.0
        Thermal expansion coefficient.
    tref : float, default 0.0
        Reference temperature (stress-free state).
    temperature : float or None
        Current temperature. If provided and alpha != 0, returns
        (gmn, u_thermal_dep) tuple.

    Returns
    -------
    gmn : (n_dep_dof, n_ind_total_dof) array
        If temperature is None or alpha == 0.
    (gmn, u_thermal_dep) : tuple
        If temperature is provided and alpha != 0.
        u_thermal_dep is (n_dep_dof,) — the free thermal displacement at
        dependent DOFs from the RBE3 offset expansion. Multiply by stiffness
        to get thermal load contribution. Note: in Nastran, u_dep = GMN @ u_ind
        always holds; the thermal contribution enters via the load vector.
    """
    dep_dofs = _expand_dof_string(dep_comp)
    ind_dofs = _expand_dof_string(ind_comp)
    n_ind = len(ind_nodes_xyz)

    if np.isscalar(weights):
        weights = np.full(n_ind, float(weights))
    else:
        weights = np.asarray(weights, dtype=float)

    ind_dofs_per_node = [ind_dofs] * n_ind
    gmn_rows, _ = compute_rbe3_gmn_row(
        xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights
    )

    if temperature is not None and abs(alpha) > 1e-30:
        u_thermal_dep = compute_rbe3_thermal_load(
            gmn_rows, xyz_dep, dep_dofs, ind_nodes_xyz,
            ind_dofs_per_node, alpha, tref, temperature,
        )
        return gmn_rows, u_thermal_dep
    return gmn_rows
