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
    return np.array(
        [
            [1, 0, 0, 0, z, -y],
            [0, 1, 0, -z, 0, x],
            [0, 0, 1, y, -x, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
        ],
        dtype=float,
    )


def _rotation_6x6(R: NDArray) -> NDArray:
    """Build 6x6 block-diagonal rotation from a 3x3 rotation matrix."""
    T = np.zeros((6, 6))
    T[:3, :3] = R
    T[3:, 3:] = R
    return T


def compute_rbe2_gmn_row(
    xyz_ind: NDArray,
    xyz_dep: NDArray,
    dep_dofs: list[int],
    cd_ind: NDArray | None = None,
    cd_dep: NDArray | None = None,
) -> NDArray:
    """Compute GMN rows for one dependent node of an RBE2 element.

    For RBE2, the dependent node moves rigidly with the independent node:
        u_dep = T @ u_ind

    where T is the 6x6 rigid-body transformation and the rows are selected
    based on which DOFs (CM) are constrained.

    When CD frames are specified, the transformation accounts for the
    node-local coordinate systems:
        GMN = R_dep @ T_basic @ R_ind^T

    Parameters
    ----------
    xyz_ind : (3,) array
        Position of independent node in basic coords.
    xyz_dep : (3,) array
        Position of dependent node in basic coords.
    dep_dofs : list of int (0-based)
        Dependent DOF indices from CM, e.g. [0,1,2,3,4,5] for 123456.
    cd_ind : (3, 3) array or None
        Rotation matrix from basic to independent node's CD frame.
        None means CD=0 (identity).
    cd_dep : (3, 3) array or None
        Rotation matrix from basic to dependent node's CD frame.
        None means CD=0 (identity).

    Returns
    -------
    gmn_row : (n_dep_dofs, 6) array
        Rows of the rigid-body transformation for the constrained DOFs.
        Columns correspond to all 6 DOFs of the independent node (in CD frame).
    """
    offset = xyz_dep - xyz_ind
    T = _rigid_body_matrix(offset)

    if cd_ind is not None or cd_dep is not None:
        # T is in basic frame. Transform to CD frames:
        # u_dep_cd = R_dep @ T_basic @ R_ind^T @ u_ind_cd
        if cd_ind is not None:
            R_ind_6 = _rotation_6x6(cd_ind)
            T = T @ R_ind_6.T
        if cd_dep is not None:
            R_dep_6 = _rotation_6x6(cd_dep)
            T = R_dep_6 @ T

    return T[dep_dofs, :]


def compute_rbe1_gmn_row(
    ind_nodes_xyz: NDArray,
    ind_dofs_per_node: list[list[int]],
    dep_nodes_xyz: NDArray,
    dep_dofs_per_node: list[list[int]],
) -> NDArray:
    """Compute GMN rows for an RBE1 element.

    RBE1 defines a rigid body using specific DOFs at the GN (independent)
    nodes. The GM (dependent/UM) DOFs are then interpolated from that
    rigid body. The GN DOFs must collectively provide 6 independent
    equations (or fewer for planar).

    The approach: build a rigid-body mode matrix for the independent DOFs
    relative to a reference point, invert to get the 6 rigid-body
    parameters from the independent DOFs, then compute dependent DOFs
    from those parameters.

    Parameters
    ----------
    ind_nodes_xyz : (n_ind, 3) array
        Positions of independent (GN) nodes.
    ind_dofs_per_node : list of list[int]
        For each GN node, which DOFs (0-based) are independent.
    dep_nodes_xyz : (n_dep, 3) array
        Positions of dependent (GM/UM) nodes.
    dep_dofs_per_node : list of list[int]
        For each GM node, which DOFs (0-based) are dependent.

    Returns
    -------
    gmn : (n_dep_dofs_total, n_ind_dofs_total) array
        Maps independent DOFs to dependent DOFs.
    """
    # Reference point: centroid of independent nodes
    ref_xyz = np.mean(ind_nodes_xyz, axis=0)

    # Build RB matrix for independent DOFs
    # Each row corresponds to one independent DOF
    # Columns are the 6 rigid-body parameters [T1,T2,T3,R1,R2,R3] at ref
    n_ind_dofs = sum(len(d) for d in ind_dofs_per_node)
    rb_ind = np.zeros((n_ind_dofs, 6))

    row = 0
    for i_node, dofs in enumerate(ind_dofs_per_node):
        offset = ind_nodes_xyz[i_node] - ref_xyz
        T = _rigid_body_matrix(offset)
        for dof in dofs:
            rb_ind[row, :] = T[dof, :]
            row += 1

    # Build RB matrix for dependent DOFs
    n_dep_dofs = sum(len(d) for d in dep_dofs_per_node)
    rb_dep = np.zeros((n_dep_dofs, 6))

    row = 0
    for i_node, dofs in enumerate(dep_dofs_per_node):
        offset = dep_nodes_xyz[i_node] - ref_xyz
        T = _rigid_body_matrix(offset)
        for dof in dofs:
            rb_dep[row, :] = T[dof, :]
            row += 1

    # GMN = rb_dep @ inv(rb_ind)
    # rb_ind maps 6 RB params -> ind DOFs, so inv gives RB params from ind DOFs
    # rb_dep maps 6 RB params -> dep DOFs
    rcond = 1e-12
    if n_ind_dofs == 6 and np.linalg.matrix_rank(rb_ind, tol=rcond) == 6:
        rb_inv = np.linalg.inv(rb_ind)
    else:
        rb_inv, _, _, _ = np.linalg.lstsq(rb_ind, np.eye(6), rcond=rcond)

    return rb_dep @ rb_inv


def compute_rbar_gmn_row(
    xyz_a: NDArray,
    xyz_b: NDArray,
    ind_dofs_a: list[int],
    ind_dofs_b: list[int],
    dep_dofs_a: list[int],
    dep_dofs_b: list[int],
) -> NDArray:
    """Compute GMN rows for an RBAR element.

    RBAR connects two nodes (GA, GB) with a rigid bar. CNA/CNB define
    which DOFs are independent at each node, and CMA/CMB define which
    are dependent. Together CNA+CNB should provide 6 independent DOFs
    defining the rigid body.

    Parameters
    ----------
    xyz_a : (3,) array
        Position of node GA.
    xyz_b : (3,) array
        Position of node GB.
    ind_dofs_a : list of int (0-based)
        Independent DOFs at node A (from CNA).
    ind_dofs_b : list of int (0-based)
        Independent DOFs at node B (from CNB).
    dep_dofs_a : list of int (0-based)
        Dependent DOFs at node A (from CMA).
    dep_dofs_b : list of int (0-based)
        Dependent DOFs at node B (from CMB).

    Returns
    -------
    gmn : (n_dep_dofs_total, n_ind_dofs_total) array
        Maps independent DOFs to dependent DOFs.
        Columns ordered as: [ind_dofs_a, ind_dofs_b].
        Rows ordered as: [dep_dofs_a, dep_dofs_b].
    """
    ind_nodes_xyz = []
    ind_dofs_per_node = []
    dep_nodes_xyz = []
    dep_dofs_per_node = []

    if ind_dofs_a:
        ind_nodes_xyz.append(xyz_a)
        ind_dofs_per_node.append(ind_dofs_a)
    if ind_dofs_b:
        ind_nodes_xyz.append(xyz_b)
        ind_dofs_per_node.append(ind_dofs_b)
    if dep_dofs_a:
        dep_nodes_xyz.append(xyz_a)
        dep_dofs_per_node.append(dep_dofs_a)
    if dep_dofs_b:
        dep_nodes_xyz.append(xyz_b)
        dep_dofs_per_node.append(dep_dofs_b)

    ind_nodes_xyz = np.array(ind_nodes_xyz) if ind_nodes_xyz else np.zeros((0, 3))
    dep_nodes_xyz = np.array(dep_nodes_xyz) if dep_nodes_xyz else np.zeros((0, 3))

    return compute_rbe1_gmn_row(
        ind_nodes_xyz,
        ind_dofs_per_node,
        dep_nodes_xyz,
        dep_dofs_per_node,
    )


def compute_rrod_gmn_row(
    xyz_a: NDArray,
    xyz_b: NDArray,
    cma: int,
    cmb: int,
) -> tuple[NDArray, list[tuple[int, int]]]:
    """Compute GMN row for an RROD element.

    RROD is a pin-ended rigid rod that constrains axial elongation to zero:
        e . (u_B - u_A) = 0
    where e is the unit vector from A to B.

    Exactly one of CMA or CMB specifies a single translational DOF as
    dependent. The constraint is rewritten so that DOF equals a linear
    combination of other translational DOFs at both nodes.

    Parameters
    ----------
    xyz_a : (3,) array
        Position of node GA.
    xyz_b : (3,) array
        Position of node GB.
    cma : int
        Dependent DOF component at A (0 if none at A).
    cmb : int
        Dependent DOF component at B (0 if none at B).

    Returns
    -------
    gmn_row : (1, n_ind) array
        Single row mapping independent translational DOFs to the
        dependent DOF.
    ind_dof_map : list of (node_flag, dof_0based) tuples
        node_flag: 0 for node A, 1 for node B.
        Defines which DOF each column corresponds to.
    """
    L_vec = xyz_b - xyz_a
    L = np.linalg.norm(L_vec)
    assert L > 1e-12, "RROD nodes must not be coincident"
    e = L_vec / L  # unit vector A->B

    # Determine which node/DOF is dependent
    if cma > 0:
        dep_node = 0  # A
        dep_dof = int(str(int(cma))[0]) - 1  # 0-based
    else:
        dep_node = 1  # B
        dep_dof = int(str(int(cmb))[0]) - 1  # 0-based

    assert dep_dof < 3, "RROD dependent DOF must be translational (1-3)"

    # Constraint: e . (u_B - u_A) = 0
    # => e[0]*uBx + e[1]*uBy + e[2]*uBz - e[0]*uAx - e[1]*uAy - e[2]*uAz = 0
    #
    # If dep is uA[dep_dof]:
    #   e[dep_dof]*uA[dep_dof] = e[dep_dof]*uB[dep_dof] + sum over other dirs
    #   uA[dep_dof] = (e[0]/e[dep_dof])*uB[0] + ... - (e[j]/e[dep_dof])*uA[j] for j!=dep_dof
    #
    # More cleanly: from e.(uB - uA) = 0:
    #   If dep is at A: e[d]*uA[d] = e.uB - sum_{j!=d} e[j]*uA[j]
    #     => uA[d] = (e/e[d]).uB - sum_{j!=d} (e[j]/e[d])*uA[j]
    #   If dep is at B: e[d]*uB[d] = e.uA - sum_{j!=d} e[j]*uB[j]
    #     => uB[d] = (e/e[d]).uA - sum_{j!=d} (e[j]/e[d])*uB[j]

    ed = e[dep_dof]
    if abs(ed) < 1e-12:
        dof_names = ["T1(x)", "T2(y)", "T3(z)"]
        valid_dofs = [i + 1 for i in range(3) if abs(e[i]) > 1e-12]
        raise ValueError(
            f"RROD: cannot solve for dependent DOF {dof_names[dep_dof]} "
            f"because the rod has near-zero projection in that direction "
            f"(e={e}). Valid dependent DOF choices: {valid_dofs}"
        )

    # Build independent DOF list and coefficients
    ind_dof_map = []  # (node_flag, dof_0based)
    coeffs = []

    if dep_node == 0:
        # Dependent at A: uA[d] = (e/ed).uB - sum_{j!=d} (e[j]/ed)*uA[j]
        # Independent: other T DOFs at A, then all T DOFs at B
        for j in range(3):
            if j != dep_dof:
                ind_dof_map.append((0, j))
                coeffs.append(-e[j] / ed)
        for j in range(3):
            ind_dof_map.append((1, j))
            coeffs.append(e[j] / ed)
    else:
        # Dependent at B: uB[d] = (e/ed).uA - sum_{j!=d} (e[j]/ed)*uB[j]
        # Independent: all T DOFs at A, then other T DOFs at B
        for j in range(3):
            ind_dof_map.append((0, j))
            coeffs.append(e[j] / ed)
        for j in range(3):
            if j != dep_dof:
                ind_dof_map.append((1, j))
                coeffs.append(-e[j] / ed)

    gmn_row = np.array([coeffs])
    return gmn_row, ind_dof_map


def compute_rspline_gmn_row(
    xyz_i1: NDArray,
    xyz_i2: NDArray,
    xyz_dep: NDArray,
    dep_dofs: list[int] | None = None,
) -> NDArray:
    """Compute GMN rows for an RSPLINE element (cubic Hermite spline).

    RSPLINE interpolates a dependent grid's DOFs between two independent
    endpoint grids using cubic Hermite polynomials for transverse directions
    and linear interpolation for axial/torsion.

    Validated against NX Nastran 2412 (CD=0 case, machine-precision match).

    Notes
    -----
    MYSTRAN's RSPLINE_PROC.f90 (lines 759-781) has a bug in the
    basic-to-global coordinate transformation when the two independent
    endpoint grids have *different* CD frames. MYSTRAN groups the
    column-side rotation by row block instead of by column ownership:

        - FR13/FR14 (I2 columns) incorrectly get TR_I1 (should be TR_I2)
        - FR21/FR22 (I1 columns) incorrectly get TR_I2 (should be TR_I1)

    The bug is silent when both independent grids share the same CD frame
    (the common case). This implementation operates in the basic frame
    and does not apply per-grid CD rotations here (those would be handled
    by the global GMN assembly's CD rotation pass).

    Parameters
    ----------
    xyz_i1 : (3,) array
        Position of independent grid 1 (first endpoint) in basic coords.
    xyz_i2 : (3,) array
        Position of independent grid 2 (second endpoint) in basic coords.
    xyz_dep : (3,) array
        Position of the dependent grid in basic coords.
    dep_dofs : list of int (0-based), optional
        Which DOFs of the dependent grid are constrained.
        Default is all 6: [0, 1, 2, 3, 4, 5].

    Returns
    -------
    gmn_rows : (n_dep_dofs, 12) array
        Rows of the spline interpolation matrix in basic frame.
        Columns: [Tx1, Ty1, Tz1, Rx1, Ry1, Rz1, Tx2, Ty2, Tz2, Rx2, Ry2, Rz2]
    """
    if dep_dofs is None:
        dep_dofs = [0, 1, 2, 3, 4, 5]

    v12 = xyz_i2 - xyz_i1
    v1d = xyz_dep - xyz_i1
    L12 = np.linalg.norm(v12)
    L1D = np.linalg.norm(v1d)

    assert L12 > 1e-12, "RSPLINE independent grids must not be coincident"
    zeta = L1D / L12

    # Build element coordinate system (x along I1->I2)
    x_hat = v12 / L12

    # Y-axis: minimum-component perpendicular (same as MYSTRAN)
    abs_x = np.abs(x_hat)
    idx_sorted = np.argsort(abs_x)
    i_min = idx_sorted[0]
    i_mid = idx_sorted[1]
    i_max = idx_sorted[2]
    y_vec = np.zeros(3)
    y_vec[i_mid] = x_hat[i_max]
    y_vec[i_max] = -x_hat[i_mid]
    y_hat = y_vec / np.linalg.norm(y_vec)

    z_hat = np.cross(x_hat, y_hat)
    z_hat /= np.linalg.norm(z_hat)

    # TRSPLINE: rows are element axes in basic coords
    # Transforms basic -> element: v_elem = T @ v_basic
    T = np.array([x_hat, y_hat, z_hat])

    # Build FR matrix in element coordinates
    FR_elem = _rspline_fr_element(zeta, L12)

    # Transform element -> basic (similarity on each 3x3 block)
    FR_basic = np.zeros((6, 12))
    for i_blk in range(2):
        for j_blk in range(4):
            r = slice(i_blk * 3, i_blk * 3 + 3)
            c = slice(j_blk * 3, j_blk * 3 + 3)
            FR_basic[r, c] = T.T @ FR_elem[r, c] @ T

    return FR_basic[dep_dofs, :]


def _rspline_fr_element(zeta: float, L: float) -> NDArray:
    """Build the 6x12 RSPLINE interpolation matrix in element coords.

    Uses cubic Hermite polynomials for transverse (y, z) and linear
    interpolation for axial (x) and torsion (Rx).

    Parameters
    ----------
    zeta : float
        Normalized distance along spline, in [0, 1].
    L : float
        Total spline length (distance between independent endpoints).

    Returns
    -------
    FR : (6, 12) interpolation matrix in element coordinates.
    """
    z = zeta
    z2 = z**2
    z3 = z**3

    FR = np.zeros((6, 12))

    # Hermite basis functions
    H00 = 1.0 - 3.0 * z2 + 2.0 * z3
    H10 = L * (z - 2.0 * z2 + z3)
    H01 = 3.0 * z2 - 2.0 * z3
    H11 = L * (z3 - z2)

    # Derivatives (dH/dx = (1/L) * dH/dz)
    dH00 = (6.0 * z2 - 6.0 * z) / L
    dH10 = 1.0 - 4.0 * z + 3.0 * z2
    dH01 = (6.0 * z - 6.0 * z2) / L
    dH11 = 3.0 * z2 - 2.0 * z

    # Row 0: Tx_dep (axial — linear)
    FR[0, 0] = 1.0 - z  # from Tx1
    FR[0, 6] = z  # from Tx2

    # Row 1: Ty_dep (transverse y — cubic Hermite)
    FR[1, 1] = H00  # from Ty1
    FR[1, 5] = H10  # from Rz1 (slope in y from rotation about z)
    FR[1, 7] = H01  # from Ty2
    FR[1, 11] = H11  # from Rz2

    # Row 2: Tz_dep (transverse z — cubic Hermite, opposite sign coupling)
    FR[2, 2] = H00  # from Tz1
    FR[2, 4] = -H10  # from Ry1 (negative: right-hand rule)
    FR[2, 8] = H01  # from Tz2
    FR[2, 10] = -H11  # from Ry2

    # Row 3: Rx_dep (torsion — linear)
    FR[3, 3] = 1.0 - z  # from Rx1
    FR[3, 9] = z  # from Rx2

    # Row 4: Ry_dep (rotation about y: Ry = -dTz/dx, right-hand rule)
    FR[4, 2] = -dH00  # from Tz1
    FR[4, 4] = dH10  # from Ry1
    FR[4, 8] = -dH01  # from Tz2
    FR[4, 10] = dH11  # from Ry2

    # Row 5: Rz_dep (rotation about z — derivative of Ty curve: Rz = dTy/dx)
    FR[5, 1] = dH00  # from Ty1
    FR[5, 5] = dH10  # from Rz1
    FR[5, 7] = dH01  # from Ty2
    FR[5, 11] = dH11  # from Rz2

    return FR


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
        rbw = rb.T * W  # (6, n_ind)
        A = rbw @ rb  # (6, 6) normal equation matrix

        rcond = 1e-12
        if np.linalg.matrix_rank(A, tol=rcond * np.max(np.abs(A))) < 6:
            gmn_full, _, _, _ = np.linalg.lstsq(A, rbw, rcond=rcond)
        else:
            gmn_full = np.linalg.solve(A, rbw)  # (6, n_ind)

        gmn_rows = gmn_full[dep_dofs, :]
        return gmn_rows, all_ind_dof_list

    # UM path: UM DOFs are moved from n-set to m-set. The RBE3 weighted LS
    # uses ALL listed independent DOFs (including UM ones) to determine the
    # 6 rigid-body parameters. Then:
    #   - Refgrid GMN rows map n-set (non-UM) columns to refgrid dep DOFs,
    #     with UM column contributions folded into the RMM coupling.
    #   - UM DOFs become m-set but have no explicit constraint equation from
    #     the RBE3 — their values are determined by the solver's stiffness
    #     coupling during MPC reduction.
    #
    # For the GMN matrix (which maps n-set → m-set), we:
    #   1. Compute gmn_full (6 × n_all_ind) from the weighted LS
    #   2. Partition columns into n-set and UM
    #   3. The refgrid constraint is: u_ref = G_n @ u_n + G_um @ u_um
    #      Rearranged: [I, -G_um] @ [u_ref; u_um] = G_n @ u_n
    #      So: GMN for refgrid rows = inv([I, -G_um]) @ [G_n; 0] (block form)
    #   4. UM DOFs have no RBE3-derived GMN rows (they need the stiffness)

    # Standard weighted LS (same as no-UM case)
    rbw = rb.T * W
    A = rbw @ rb

    rcond = 1e-12
    if np.linalg.matrix_rank(A, tol=rcond * np.max(np.abs(A))) < 6:
        gmn_full, _, _, _ = np.linalg.lstsq(A, rbw, rcond=rcond)
    else:
        gmn_full = np.linalg.solve(A, rbw)  # (6, n_all_ind)

    # Identify which columns are UM DOFs
    um_node_dof_set = set()
    for i_um, um_xyz in enumerate(um_nodes_xyz):
        for i_node in range(n_ind_nodes):
            if np.allclose(ind_nodes_xyz[i_node], um_xyz, atol=1e-10):
                for dof in um_dofs_per_node[i_um]:
                    um_node_dof_set.add((i_node, dof))
                break

    # Partition columns: n-set vs UM
    n_set_cols = []
    um_cols = []
    n_set_dof_list = []
    for i_col, (node_idx, dof_idx) in enumerate(all_ind_dof_list):
        if (node_idx, dof_idx) in um_node_dof_set:
            um_cols.append(i_col)
        else:
            n_set_cols.append(i_col)
            n_set_dof_list.append((node_idx, dof_idx))

    # G_n: refgrid DOFs from n-set columns, G_um: from UM columns
    G_n = gmn_full[dep_dofs, :][:, n_set_cols]
    G_um = gmn_full[dep_dofs, :][:, um_cols]

    n_ref = len(dep_dofs)
    n_um = len(um_cols)
    n_m = n_ref + n_um

    # Build RMM and solve for GMN
    # Constraint: u_ref - G_um @ u_um = G_n @ u_n
    # RMM = [I_ref, -G_um; 0, I_um], RMN = [-G_n; 0]
    # GMN = -inv(RMM) @ RMN
    RMM = np.eye(n_m)
    RMM[:n_ref, n_ref:] = -G_um
    RMN = np.zeros((n_m, len(n_set_cols)))
    RMN[:n_ref, :] = -G_n

    gmn_rows = -np.linalg.solve(RMM, RMN)
    return gmn_rows, n_set_dof_list


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
        xyz_dep,
        dep_dofs,
        ind_nodes_xyz,
        ind_dofs_per_node,
        weights_per_node,
    )
    n_dep = len(dep_dofs)
    n_ind = gmn_rows.shape[1]

    # For standard RBE3: constraint is u_m - GMN @ u_n = 0
    # So RMM = I, RMN = -GMN
    RMM = np.eye(n_dep)
    RMN = -gmn_rows

    return RMM, RMN, gmn_rows


def _build_cd_rotations(model) -> dict[int, NDArray]:
    """Build 3x3 CD rotation matrices for each node with non-zero CD.

    Returns
    -------
    cd_rotations : dict mapping node_id -> (3, 3) rotation matrix
        Only nodes with CD != 0 are included. The rotation matrix R
        transforms from basic to local: u_local = R @ u_basic.
    """
    grid = model.grid
    cd_values = grid.cd
    node_ids = grid.node_id
    coord = model.coord

    cd_rotations: dict[int, NDArray] = {}
    for i_node, (nid, cd) in enumerate(zip(node_ids, cd_values)):
        if cd == 0:
            continue
        idx = coord.index(cd)
        R = np.array([coord.i[idx], coord.j[idx], coord.k[idx]])
        cd_rotations[int(nid)] = R
    return cd_rotations


def _apply_cd_to_gmn(
    GMN_dense: NDArray,
    m_set: dict[tuple[int, int], int],
    dof_map: dict[tuple[int, int], int],
    cd_rotations: dict[int, NDArray],
) -> NDArray:
    """Apply CD frame rotations to a dense GMN matrix.

    GMN_basic operates on basic-frame DOFs. This transforms it to operate
    on CD-frame DOFs:
        GMN_cd[m_row, n_col] = R_dep @ GMN_basic @ R_ind^T

    Applied per 3-DOF block (T and R separately).
    """
    if not cd_rotations:
        return GMN_dense

    n_m, n_g = GMN_dense.shape

    # Apply column rotations (independent/n-set side): GMN @ R_ind^T
    # Group columns by node, apply R^T to each 3-DOF block
    nid_to_cols: dict[int, list[int]] = {}
    for (nid, dof), col_idx in dof_map.items():
        nid_to_cols.setdefault(nid, [None] * 6)
        nid_to_cols[nid][dof - 1] = col_idx

    for nid, R in cd_rotations.items():
        if nid not in nid_to_cols:
            continue
        cols = nid_to_cols[nid]
        # T block (DOFs 1-3)
        t_cols = [cols[i] for i in range(3) if cols[i] is not None]
        if len(t_cols) == 3:
            GMN_dense[:, t_cols] = GMN_dense[:, t_cols] @ R.T
        # R block (DOFs 4-6)
        r_cols = [cols[i] for i in range(3, 6) if cols[i] is not None]
        if len(r_cols) == 3:
            GMN_dense[:, r_cols] = GMN_dense[:, r_cols] @ R.T

    # Apply row rotations (dependent/m-set side): R_dep @ GMN
    # Group rows by node
    nid_to_rows: dict[int, list[int | None]] = {}
    for (nid, dof), row_idx in m_set.items():
        nid_to_rows.setdefault(nid, [None] * 6)
        nid_to_rows[nid][dof - 1] = row_idx

    for nid, R in cd_rotations.items():
        if nid not in nid_to_rows:
            continue
        rows_map = nid_to_rows[nid]
        # T block (DOFs 1-3)
        t_rows = [rows_map[i] for i in range(3) if rows_map[i] is not None]
        if len(t_rows) == 3:
            GMN_dense[t_rows, :] = R @ GMN_dense[t_rows, :]
        # R block (DOFs 4-6)
        r_rows = [rows_map[i] for i in range(3, 6) if rows_map[i] is not None]
        if len(r_rows) == 3:
            GMN_dense[r_rows, :] = R @ GMN_dense[r_rows, :]

    return GMN_dense


def assemble_gmn(
    model,
    dof_map: dict[tuple[int, int], int] | None = None,
    ndof: int | None = None,
    apply_cd: bool = True,
) -> tuple[NDArray, dict, dict]:
    """Assemble global GMN matrix from all rigid elements in the model.

    Handles RBAR, RBAR1, RBE1, RBE2, RBE3, and RROD elements.

    Parameters
    ----------
    model : pyNastran BDF (vectorized3)
        Must have at least one rigid element and grid coordinates.
    dof_map : dict mapping (node_id, dof_1based) -> global_index, optional
        If None, builds one from all grids (6 DOF per grid).
    ndof : int, optional
        Total DOFs. If None, determined from dof_map.
    apply_cd : bool, default True
        If True and any nodes have non-zero CD, apply coordinate frame
        rotations so GMN operates on CD-frame DOFs (matching Nastran).

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
    n_rbe1 = model.rbe1.n if hasattr(model, "rbe1") else 0
    n_rbe2 = model.rbe2.n if hasattr(model, "rbe2") else 0
    n_rbe3 = model.rbe3.n if hasattr(model, "rbe3") else 0
    n_rbar = model.rbar.n if hasattr(model, "rbar") else 0
    n_rbar1 = model.rbar1.n if hasattr(model, "rbar1") else 0
    n_rrod = model.rrod.n if hasattr(model, "rrod") else 0
    n_total = n_rbe1 + n_rbe2 + n_rbe3 + n_rbar + n_rbar1 + n_rrod
    if n_total == 0:
        raise ValueError("Model has no rigid elements")

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

    m_set: dict[tuple[int, int], int] = {}
    m_row = 0

    nid_to_idx = {int(nid): i for i, nid in enumerate(node_ids)}

    rows_list: list[int] = []
    cols_list: list[int] = []
    vals_list: list[float] = []

    # --- RBE2 contributions ---
    if n_rbe2 > 0:
        rbe2 = model.rbe2
        for i_elem in range(rbe2.n):
            ind_nid = int(rbe2.independent_node[i_elem])
            cm = int(rbe2.independent_dof[i_elem])
            dep_dofs = _expand_dof_string(cm)
            xyz_ind = xyz_cid0[nid_to_idx[ind_nid]]

            idim0, idim1 = rbe2.idim[i_elem]
            dep_nids = rbe2.dependent_nodes[idim0:idim1]

            for dep_nid in dep_nids:
                dep_nid = int(dep_nid)
                xyz_dep = xyz_cid0[nid_to_idx[dep_nid]]

                gmn_row = compute_rbe2_gmn_row(xyz_ind, xyz_dep, dep_dofs)

                m_rows_this = []
                for dof in dep_dofs:
                    key = (dep_nid, dof + 1)
                    if key not in m_set:
                        m_set[key] = m_row
                        m_row += 1
                    m_rows_this.append(m_set[key])

                for i_row, m_row_idx in enumerate(m_rows_this):
                    for j_col in range(6):
                        val = gmn_row[i_row, j_col]
                        if abs(val) > 1e-15:
                            col_idx = dof_map.get((ind_nid, j_col + 1))
                            if col_idx is not None:
                                rows_list.append(m_row_idx)
                                cols_list.append(col_idx)
                                vals_list.append(val)

    # --- RBE1 contributions ---
    if n_rbe1 > 0:
        rbe1 = model.rbe1
        for i_elem in range(rbe1.n):
            # RBE1 naming in vectorized3 storage is swapped:
            # dependent_node/dof = GN nodes (independent in Nastran terms)
            # independent_node/dof = GM/UM nodes (dependent in Nastran terms)
            idep0, idep1 = rbe1.idependent[i_elem]
            iind0, iind1 = rbe1.iindependent[i_elem]

            gn_nids = rbe1.dependent_node[idep0:idep1]
            gn_comps = rbe1.dependent_dof[idep0:idep1]
            gm_nids = rbe1.independent_node[iind0:iind1]
            gm_comps = rbe1.independent_dof[iind0:iind1]

            # GN = independent (define rigid body)
            ind_nodes_xyz_list = []
            ind_dofs_per_node_list = []
            ind_global_dof_map = []
            for gn_nid, gn_comp in zip(gn_nids, gn_comps):
                gn_nid = int(gn_nid)
                dofs = _expand_dof_string(int(gn_comp))
                ind_nodes_xyz_list.append(xyz_cid0[nid_to_idx[gn_nid]])
                ind_dofs_per_node_list.append(dofs)
                for dof in dofs:
                    ind_global_dof_map.append((gn_nid, dof + 1))

            # GM = dependent (UM nodes)
            dep_nodes_xyz_list = []
            dep_dofs_per_node_list = []
            dep_global_dof_map = []
            for gm_nid, gm_comp in zip(gm_nids, gm_comps):
                gm_nid = int(gm_nid)
                dofs = _expand_dof_string(int(gm_comp))
                dep_nodes_xyz_list.append(xyz_cid0[nid_to_idx[gm_nid]])
                dep_dofs_per_node_list.append(dofs)
                for dof in dofs:
                    dep_global_dof_map.append((gm_nid, dof + 1))

            gmn_rows = compute_rbe1_gmn_row(
                np.array(ind_nodes_xyz_list),
                ind_dofs_per_node_list,
                np.array(dep_nodes_xyz_list),
                dep_dofs_per_node_list,
            )

            # Register m-set DOFs and place into global matrix
            m_rows_this = []
            for nid_dep, dof_dep in dep_global_dof_map:
                key = (nid_dep, dof_dep)
                if key not in m_set:
                    m_set[key] = m_row
                    m_row += 1
                m_rows_this.append(m_set[key])

            for i_row, m_row_idx in enumerate(m_rows_this):
                for i_col, (nid_col, dof_col) in enumerate(ind_global_dof_map):
                    val = gmn_rows[i_row, i_col]
                    if abs(val) > 1e-15:
                        col_idx = dof_map.get((nid_col, dof_col))
                        if col_idx is not None:
                            rows_list.append(m_row_idx)
                            cols_list.append(col_idx)
                            vals_list.append(val)

    # --- RBAR contributions ---
    if n_rbar > 0:
        rbar = model.rbar
        for i_elem in range(rbar.n):
            nid_a = int(rbar.nodes[i_elem, 0])
            nid_b = int(rbar.nodes[i_elem, 1])
            cna = int(rbar.independent_dof[i_elem, 0])
            cnb = int(rbar.independent_dof[i_elem, 1])
            cma = int(rbar.dependent_dof[i_elem, 0])
            cmb = int(rbar.dependent_dof[i_elem, 1])

            xyz_a = xyz_cid0[nid_to_idx[nid_a]]
            xyz_b = xyz_cid0[nid_to_idx[nid_b]]

            ind_dofs_a = _expand_dof_string(cna) if cna > 0 else []
            ind_dofs_b = _expand_dof_string(cnb) if cnb > 0 else []
            dep_dofs_a = _expand_dof_string(cma) if cma > 0 else []
            dep_dofs_b = _expand_dof_string(cmb) if cmb > 0 else []

            if not dep_dofs_a and not dep_dofs_b:
                continue

            gmn_rows = compute_rbar_gmn_row(
                xyz_a, xyz_b, ind_dofs_a, ind_dofs_b, dep_dofs_a, dep_dofs_b
            )

            # Build column mapping for independent DOFs
            ind_global_dof_map = []
            for dof in ind_dofs_a:
                ind_global_dof_map.append((nid_a, dof + 1))
            for dof in ind_dofs_b:
                ind_global_dof_map.append((nid_b, dof + 1))

            # Build row mapping for dependent DOFs
            dep_global_dof_map = []
            for dof in dep_dofs_a:
                dep_global_dof_map.append((nid_a, dof + 1))
            for dof in dep_dofs_b:
                dep_global_dof_map.append((nid_b, dof + 1))

            m_rows_this = []
            for nid_dep, dof_dep in dep_global_dof_map:
                key = (nid_dep, dof_dep)
                if key not in m_set:
                    m_set[key] = m_row
                    m_row += 1
                m_rows_this.append(m_set[key])

            for i_row, m_row_idx in enumerate(m_rows_this):
                for i_col, (nid_col, dof_col) in enumerate(ind_global_dof_map):
                    val = gmn_rows[i_row, i_col]
                    if abs(val) > 1e-15:
                        col_idx = dof_map.get((nid_col, dof_col))
                        if col_idx is not None:
                            rows_list.append(m_row_idx)
                            cols_list.append(col_idx)
                            vals_list.append(val)

    # --- RBE3 contributions ---
    if n_rbe3 > 0:
        rbe3 = model.rbe3
        iweight_all = rbe3.iweight
        igrid_start = 0

        for i_elem in range(rbe3.n):
            ref_nid = int(rbe3.ref_grid[i_elem])
            ref_comp = int(rbe3.ref_component[i_elem])
            dep_dofs = _expand_dof_string(ref_comp)

            xyz_dep = xyz_cid0[nid_to_idx[ref_nid]]

            m_rows_this = []
            for dof in dep_dofs:
                key = (ref_nid, dof + 1)
                if key not in m_set:
                    m_set[key] = m_row
                    m_row += 1
                m_rows_this.append(m_set[key])

            iw0, iw1 = iweight_all[i_elem]
            nweights = iw1 - iw0

            ind_nodes_xyz_list = []
            ind_dofs_per_node_list = []
            weights_per_node_list = []
            ind_global_dof_map = []

            for iw in range(nweights):
                weight = float(rbe3.weight[iw0 + iw])
                comp = int(rbe3.independent_dofs[iw0 + iw])
                active_dofs = _expand_dof_string(comp)

                ngrid = int(rbe3.ngrid_per_weight[igrid_start + iw])
                for ig in range(ngrid):
                    flat_idx = sum(rbe3.ngrid_per_weight[igrid_start : igrid_start + iw]) + ig
                    nid = int(
                        rbe3.independent_nodes[sum(rbe3.ngrid_per_weight[:igrid_start]) + flat_idx]
                    )
                    xyz_ind = xyz_cid0[nid_to_idx[nid]]
                    ind_nodes_xyz_list.append(xyz_ind)
                    ind_dofs_per_node_list.append(active_dofs)
                    weights_per_node_list.append(weight)
                    for dof in active_dofs:
                        ind_global_dof_map.append((nid, dof + 1))

            ind_nodes_xyz = np.array(ind_nodes_xyz_list)
            weights_per_node = np.array(weights_per_node_list)

            gmn_rows, _ = compute_rbe3_gmn_row(
                xyz_dep,
                dep_dofs,
                ind_nodes_xyz,
                ind_dofs_per_node_list,
                weights_per_node,
            )

            for i_row, m_row_idx in enumerate(m_rows_this):
                for i_col, (nid_col, dof_col) in enumerate(ind_global_dof_map):
                    val = gmn_rows[i_row, i_col]
                    if abs(val) > 1e-15:
                        col_idx = dof_map.get((nid_col, dof_col))
                        if col_idx is not None:
                            rows_list.append(m_row_idx)
                            cols_list.append(col_idx)
                            vals_list.append(val)

            igrid_start += nweights

    # --- RBAR1 contributions ---
    # RBAR1: GA is always independent (all 6 DOFs), GB has CB as dependent
    if n_rbar1 > 0:
        rbar1 = model.rbar1
        for i_elem in range(rbar1.n):
            nid_a = int(rbar1.nodes[i_elem, 0])
            nid_b = int(rbar1.nodes[i_elem, 1])
            cb = int(rbar1.dependent_dof[i_elem])

            xyz_a = xyz_cid0[nid_to_idx[nid_a]]
            xyz_b = xyz_cid0[nid_to_idx[nid_b]]

            dep_dofs_b = _expand_dof_string(cb) if cb > 0 else []
            if not dep_dofs_b:
                continue

            # RBAR1: same as RBE2 — GA independent, GB dependent
            gmn_row = compute_rbe2_gmn_row(xyz_a, xyz_b, dep_dofs_b)

            m_rows_this = []
            for dof in dep_dofs_b:
                key = (nid_b, dof + 1)
                if key not in m_set:
                    m_set[key] = m_row
                    m_row += 1
                m_rows_this.append(m_set[key])

            for i_row, m_row_idx in enumerate(m_rows_this):
                for j_col in range(6):
                    val = gmn_row[i_row, j_col]
                    if abs(val) > 1e-15:
                        col_idx = dof_map.get((nid_a, j_col + 1))
                        if col_idx is not None:
                            rows_list.append(m_row_idx)
                            cols_list.append(col_idx)
                            vals_list.append(val)

    # --- RROD contributions ---
    if n_rrod > 0:
        rrod = model.rrod
        for i_elem in range(rrod.n):
            nid_a = int(rrod.nodes[i_elem, 0])
            nid_b = int(rrod.nodes[i_elem, 1])
            cma = int(rrod.dependent_dof[i_elem, 0])
            cmb = int(rrod.dependent_dof[i_elem, 1])

            if cma == 0 and cmb == 0:
                continue

            xyz_a = xyz_cid0[nid_to_idx[nid_a]]
            xyz_b = xyz_cid0[nid_to_idx[nid_b]]

            gmn_row, ind_dof_map_local = compute_rrod_gmn_row(xyz_a, xyz_b, cma, cmb)

            # Determine dependent DOF
            if cma > 0:
                dep_nid = nid_a
                dep_dof_0 = int(str(int(cma))[0]) - 1
            else:
                dep_nid = nid_b
                dep_dof_0 = int(str(int(cmb))[0]) - 1

            key = (dep_nid, dep_dof_0 + 1)
            if key not in m_set:
                m_set[key] = m_row
                m_row += 1
            m_row_idx = m_set[key]

            nids_ab = [nid_a, nid_b]
            for i_col, (node_flag, dof_0) in enumerate(ind_dof_map_local):
                val = gmn_row[0, i_col]
                if abs(val) > 1e-15:
                    col_idx = dof_map.get((nids_ab[node_flag], dof_0 + 1))
                    if col_idx is not None:
                        rows_list.append(m_row_idx)
                        cols_list.append(col_idx)
                        vals_list.append(val)

    # Build sparse GMN (in basic frame)
    n_m = m_row
    GMN = sparse.csr_matrix(
        (vals_list, (rows_list, cols_list)),
        shape=(n_m, ndof),
    )

    # Apply CD frame rotations if requested
    if apply_cd:
        cd_rotations = _build_cd_rotations(model)
        if cd_rotations:
            GMN_dense = _apply_cd_to_gmn(GMN.toarray(), m_set, dof_map, cd_rotations)
            GMN = sparse.csr_matrix(GMN_dense)

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
    gmn_rows, _ = compute_rbe3_gmn_row(xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights)

    if temperature is not None and abs(alpha) > 1e-30:
        u_thermal_dep = compute_rbe3_thermal_load(
            gmn_rows,
            xyz_dep,
            dep_dofs,
            ind_nodes_xyz,
            ind_dofs_per_node,
            alpha,
            tref,
            temperature,
        )
        return gmn_rows, u_thermal_dep
    return gmn_rows
