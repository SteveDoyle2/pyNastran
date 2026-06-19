"""Residual vector (RESVEC) computation for modal augmentation.

Computes static correction vectors that capture the response of truncated
modes to applied loads. Orthogonalizes against the retained modal basis
and appends as pseudo-modes.

References
----------
- Dickens, Nakagawa, Wittbrodt (1997). "A Critique of Mode Acceleration and
  Modal Truncation Augmentation Methods for Modal Response Analysis."
  Computers & Structures, 62(6), pp. 985-998.
- Rose (1992). "Using Residual Vectors in MSC/NASTRAN Dynamic Analysis to
  Improve Accuracy." 33rd SDM Conference, AIAA-92-2448.
- Simcenter Nastran Basic Dynamic Analysis User's Guide, Ch. 4 (PARAM,RESVEC).
"""

import numpy as np
from scipy.sparse import issparse, csc_matrix
from scipy.sparse.linalg import splu


def compute_residual_vectors(
    Kaa: np.ndarray | csc_matrix,
    Maa: np.ndarray | csc_matrix,
    phi: np.ndarray,
    eigenvalues: np.ndarray,
    load_vectors: np.ndarray,
    tiny: float = 1.0e-8,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute residual vectors for modal basis augmentation.

    Given a truncated set of eigenvectors (phi) from a normal modes analysis,
    computes residual vectors that capture the static response to specified
    loads NOT already represented by the retained modes.

    Parameters
    ----------
    Kaa : ndarray or sparse matrix, shape (n, n)
        Stiffness matrix at the analysis (a-set) DOFs. Must be non-singular
        (SPCs already applied).
    Maa : ndarray or sparse matrix, shape (n, n)
        Mass matrix at the analysis (a-set) DOFs.
    phi : ndarray, shape (n, p)
        Mass-normalized eigenvectors (columns). Must satisfy phi.T @ Maa @ phi = I.
    eigenvalues : ndarray, shape (p,)
        Eigenvalues (omega_i^2) corresponding to each mode in phi.
    load_vectors : ndarray, shape (n, n_loads)
        Load patterns for which residual vectors are desired. Each column is
        one independent load vector (e.g., from DAREA, unit forces at USET U6
        DOFs, or M*{e_d} for enforced motion directions).
    tiny : float, default 1.0e-8
        Threshold for discarding linearly dependent residual vectors.
        A vector is discarded if its generalized mass after orthogonalization
        is less than tiny.

    Returns
    -------
    psi : ndarray, shape (n, r)
        Mass-normalized residual vectors (columns). r <= n_loads.
    omega_pseudo : ndarray, shape (r,)
        Pseudo-eigenvalues (omega^2) for each residual vector, computed
        from psi.T @ K @ psi. These are typically much larger than the
        highest retained elastic eigenvalue.

    Notes
    -----
    The algorithm follows the Nastran RESVEC implementation:

    1. Solve K * u_static = F for each load pattern (reuses factored K).
    2. Subtract the portion captured by retained modes:
       psi_raw = u_static - Phi * (Phi.T @ F / omega^2)
    3. M-orthogonalize against all retained modes (modified Gram-Schmidt).
    4. M-orthogonalize residual vectors against each other.
    5. Discard vectors with generalized mass < tiny.
    6. Mass-normalize survivors.
    7. Compute pseudo-frequencies from generalized stiffness.

    The returned vectors can be appended to phi to form an augmented modal
    basis that exactly captures the static response to the specified loads,
    regardless of modal truncation level.
    """
    n, p = phi.shape
    n_loads = load_vectors.shape[1]
    assert load_vectors.shape[0] == n

    # Step 1: Factorize K (or reuse if already factored)
    if issparse(Kaa):
        K_factor = splu(csc_matrix(Kaa))
        solve_K = K_factor.solve
    else:
        from scipy.linalg import cho_factor, cho_solve

        K_cho = cho_factor(Kaa)

        def solve_K(b):
            return cho_solve(K_cho, b)

    # Step 2: Compute static displacements and subtract modal content
    # Modal projection: sum_j phi_j * (phi_j^T * F_k) / omega_j^2
    # = Phi @ diag(1/omega^2) @ Phi^T @ F
    inv_omega2 = 1.0 / eigenvalues  # (p,)
    phi_t_F = phi.T @ load_vectors  # (p, n_loads)
    modal_content = phi @ (inv_omega2[:, np.newaxis] * phi_t_F)  # (n, n_loads)

    psi_raw_list = []
    for k in range(n_loads):
        Fk = load_vectors[:, k]
        u_static = solve_K(Fk)
        psi_raw = u_static - modal_content[:, k]
        psi_raw_list.append(psi_raw)

    # Steps 3-7: Orthogonalize, filter, normalize
    psi_accepted = []
    omega2_accepted = []

    # Precompute M @ phi for repeated use in orthogonalization
    if issparse(Maa):
        M_phi = Maa @ phi  # (n, p) sparse @ dense = dense
    else:
        M_phi = Maa @ phi  # (n, p)

    for psi_raw in psi_raw_list:
        psi = psi_raw.copy()

        # Step 3: M-orthogonalize against retained elastic modes
        for j in range(p):
            alpha = psi @ M_phi[:, j]
            psi -= alpha * phi[:, j]

        # Step 4: M-orthogonalize against previously accepted residual vectors
        for i, psi_prev in enumerate(psi_accepted):
            if issparse(Maa):
                M_psi_prev = Maa @ psi_prev
            else:
                M_psi_prev = Maa @ psi_prev
            beta = psi @ M_psi_prev
            psi -= beta * psi_prev

        # Step 5: Check generalized mass (linear independence)
        if issparse(Maa):
            M_psi = Maa @ psi
        else:
            M_psi = Maa @ psi
        gen_mass = psi @ M_psi

        if gen_mass < tiny:
            continue

        # Step 6: Mass-normalize
        psi /= np.sqrt(gen_mass)

        # Step 7: Compute pseudo-eigenvalue (generalized stiffness)
        if issparse(Kaa):
            K_psi = Kaa @ psi
        else:
            K_psi = Kaa @ psi
        omega2_pseudo = psi @ K_psi

        psi_accepted.append(psi)
        omega2_accepted.append(omega2_pseudo)

    if len(psi_accepted) == 0:
        return np.empty((n, 0)), np.empty((0,))

    psi_out = np.column_stack(psi_accepted)
    omega2_out = np.array(omega2_accepted)
    return psi_out, omega2_out


def build_inertia_load_vectors(
    Maa: np.ndarray | csc_matrix,
    ndof: int,
) -> np.ndarray:
    """Build inertia load vectors for enforced-motion residual vectors (RESVINER).

    Computes M * {e_d} for d = 1..6 rigid-body directions (TX, TY, TZ, RX, RY, RZ).
    Each column represents the inertia load from unit rigid-body acceleration
    in one direction.

    Parameters
    ----------
    Maa : ndarray or sparse matrix, shape (n, n)
        Mass matrix at analysis DOFs.
    ndof : int
        Number of DOFs per grid point (6 for 3D grids).

    Returns
    -------
    F_inertia : ndarray, shape (n, 6)
        Inertia load vectors. Column d is M @ {e_d} where {e_d} has 1.0 at
        every DOF in direction d (translation X/Y/Z or rotation X/Y/Z).

    Notes
    -----
    For a model with N grids and 6 DOF/grid, the rigid body vector for
    translation-X has 1.0 at DOFs 0, 6, 12, ... (every 6th DOF starting at 0).
    Rotation vectors would include moment-arm contributions for a proper
    rigid-body rotation, but Nastran's RESVINER uses simplified unit
    acceleration vectors (pure translation or pure rotation at each grid).
    """
    n = Maa.shape[0]
    n_directions = min(ndof, 6)
    F_inertia = np.zeros((n, n_directions), dtype="float64")

    for d in range(n_directions):
        # Build unit acceleration vector: 1.0 at every grid's d-th DOF
        e_d = np.zeros(n, dtype="float64")
        e_d[d::ndof] = 1.0
        result = Maa @ e_d
        if issparse(Maa):
            F_inertia[:, d] = np.asarray(result).ravel()
        else:
            F_inertia[:, d] = result

    return F_inertia


def build_load_vectors_from_subcase(
    model,
    subcase,
    dof_map: dict,
    ndof: int,
    aset: np.ndarray,
) -> np.ndarray:
    """Build applied-force load vectors from a subcase's LOADSET/LSEQ definition.

    Extracts the static load patterns referenced by LOADSET case control and
    assembles them into load vectors at analysis-set DOFs.

    Parameters
    ----------
    model : BDF
        The pyNastran BDF model (parsed).
    subcase : Subcase
        The case control subcase containing LOADSET or LOAD references.
    dof_map : dict
        Mapping of (node_id, component) -> global DOF index.
    ndof : int
        Total number of DOFs in the g-set.
    aset : ndarray of bool, shape (ndof,)
        Boolean mask for analysis-set DOFs (True = in a-set).

    Returns
    -------
    F_loads : ndarray, shape (n_a, n_loads)
        Load vectors at a-set DOFs. Each column is one load pattern.
        Returns empty (n_a, 0) array if no loads are defined.

    Notes
    -----
    The LOADSET/LSEQ mechanism in Nastran works as follows:

    - Case control: LOADSET = sid  (selects LSEQ entries with matching SID)
    - LSEQ entries map EXCITEID -> LID (static load set ID)
    - The static load set (FORCE, MOMENT, PLOAD, etc.) is assembled into
      a force vector, one per unique EXCITEID

    For RESVEC, we need one residual vector per unique spatial load pattern.
    """
    n_a = int(aset.sum())

    # Check for LOADSET -> LSEQ path
    if "LOADSET" in subcase:
        loadset_id, _ = subcase["LOADSET"]
        lseq = model.lseq
        if lseq.n == 0:
            return np.empty((n_a, 0), dtype="float64")

        # Find LSEQ entries matching this LOADSET
        idx = np.where(lseq.lseq_id == loadset_id)[0]
        if len(idx) == 0:
            return np.empty((n_a, 0), dtype="float64")

        load_ids = lseq.load_id[idx]
    elif "LOAD" in subcase:
        load_id, _ = subcase["LOAD"]
        load_ids = np.array([load_id], dtype="int32")
    else:
        return np.empty((n_a, 0), dtype="float64")

    # Build Fg for each load ID using the model's existing load assembly
    load_vectors = []
    reduced_loads = model.get_reduced_static_load()

    for lid in load_ids:
        if lid not in reduced_loads:
            continue
        Fg = np.zeros(ndof, dtype="float64")
        loads = reduced_loads[lid]
        for scale, load in loads:
            if load.type == "FORCE":
                fxyz_myz = load.sum_forces_moments()
                fxyz = fxyz_myz[:, :3]
                nids = load.node_id
                for fxyzi, nid in zip(fxyz, nids):
                    fi = dof_map[(nid, 1)]
                    Fg[fi : fi + 3] += fxyzi * scale
            elif load.type == "MOMENT":
                fxyz_myz = load.sum_forces_moments()
                mxyz = fxyz_myz[:, 3:]
                nids = load.node_id
                for mxyzi, nid in zip(mxyz, nids):
                    fi = dof_map[(nid, 4)]
                    Fg[fi : fi + 3] += mxyzi * scale
            elif load.type == "SLOAD":
                for mag, nid in zip(load.mags, load.nodes):
                    fi = dof_map[(nid, 0)]
                    Fg[fi] += mag * scale

        # Reduce to a-set (partition; MPC/OMIT handled upstream)
        Fa = Fg[aset]
        if np.any(np.abs(Fa) > 0.0):
            load_vectors.append(Fa)

    if len(load_vectors) == 0:
        return np.empty((n_a, 0), dtype="float64")
    return np.column_stack(load_vectors)


def build_uset_u6_load_vectors(
    model,
    dof_map: dict,
    ndof: int,
    aset: np.ndarray,
) -> np.ndarray:
    """Build unit force vectors at DOFs designated by USET U6 entries.

    USET U6 defines DOFs where unit forces are applied for residual vector
    computation, without requiring actual load cards (FORCE, MOMENT, etc.).

    Parameters
    ----------
    model : BDF
        The pyNastran BDF model.
    dof_map : dict
        Mapping of (node_id, component) -> global DOF index.
    ndof : int
        Total number of DOFs in the g-set.
    aset : ndarray of bool, shape (ndof,)
        Boolean mask for analysis-set DOFs.

    Returns
    -------
    F_u6 : ndarray, shape (n_a, n_u6_dofs)
        Unit force vectors at a-set DOFs. Each column is a unit force at one
        U6-designated DOF.

    Notes
    -----
    USET U6 designates DOFs for residual vector computation. Nastran applies
    a unit force at each designated DOF and computes the residual vector.
    This allows residual vectors without defining actual LOADSET/LSEQ loads.

    Example Nastran input::

        USET  U6  100  123  200  123  300  12

    This designates DOFs 1,2,3 at node 100; DOFs 1,2,3 at node 200; DOFs 1,2
    at node 300 — producing 8 unit load vectors and up to 8 residual vectors.
    """
    n_a = int(aset.sum())
    uset = model.uset

    if uset.n == 0:
        return np.empty((n_a, 0), dtype="float64")

    # Find U6 entries
    u6_mask = uset.name == "U6"
    if not np.any(u6_mask):
        return np.empty((n_a, 0), dtype="float64")

    u6_nodes = uset.node_id[u6_mask]
    u6_components = uset.component[u6_mask]

    load_vectors = []
    for nid, comp in zip(u6_nodes, u6_components):
        # comp is an integer encoding (e.g., 123 means DOFs 1, 2, 3)
        for c in str(comp):
            dof = int(c)
            if (nid, dof) not in dof_map:
                continue
            gi = dof_map[(nid, dof)]
            Fg = np.zeros(ndof, dtype="float64")
            Fg[gi] = 1.0
            Fa = Fg[aset]
            if np.any(np.abs(Fa) > 0.0):
                load_vectors.append(Fa)

    if len(load_vectors) == 0:
        return np.empty((n_a, 0), dtype="float64")
    return np.column_stack(load_vectors)


def augment_modes_with_resvec(
    Kaa: np.ndarray | csc_matrix,
    Maa: np.ndarray | csc_matrix,
    phi_a: np.ndarray,
    eigenvalues: np.ndarray,
    model,
    subcase,
    dof_map: dict,
    ndof: int,
    aset: np.ndarray,
    ndof_per_grid: int = 6,
    include_inertia: bool = True,
    tiny: float = 1.0e-8,
) -> tuple[np.ndarray, np.ndarray]:
    """Top-level function: compute and append residual vectors to a modal basis.

    Collects all applicable load patterns (LOADSET/LSEQ, USET U6, inertia)
    and computes residual vectors, returning the augmented modal basis.

    Parameters
    ----------
    Kaa : ndarray or sparse matrix, shape (n_a, n_a)
        Analysis-set stiffness matrix.
    Maa : ndarray or sparse matrix, shape (n_a, n_a)
        Analysis-set mass matrix.
    phi_a : ndarray, shape (n_a, p)
        Mass-normalized eigenvectors at a-set DOFs.
    eigenvalues : ndarray, shape (p,)
        Eigenvalues (omega^2) for the retained modes.
    model : BDF
        The pyNastran BDF model.
    subcase : Subcase
        Case control subcase.
    dof_map : dict
        (node_id, component) -> global DOF index mapping.
    ndof : int
        Total DOFs in g-set.
    aset : ndarray of bool, shape (ndof,)
        Boolean mask for a-set DOFs.
    ndof_per_grid : int, default 6
        DOFs per grid point.
    include_inertia : bool, default True
        If True, include inertia-relief residual vectors (RESVINER) for
        enforced-motion excitation. Adds 6 load vectors (unit acceleration
        in each rigid-body direction).
    tiny : float, default 1.0e-8
        Threshold for discarding dependent vectors.

    Returns
    -------
    phi_augmented : ndarray, shape (n_a, p + r)
        Augmented modal basis (original modes + residual vectors).
    eigenvalues_augmented : ndarray, shape (p + r,)
        Eigenvalues for all modes (original + pseudo-frequencies for residuals).
    """
    n_a = Kaa.shape[0]
    assert phi_a.shape[0] == n_a

    # Collect all load vectors
    all_loads = []

    # Source 1: LOADSET/LSEQ or LOAD
    F_loads = build_load_vectors_from_subcase(model, subcase, dof_map, ndof, aset)
    if F_loads.shape[1] > 0:
        all_loads.append(F_loads)

    # Source 2: USET U6
    F_u6 = build_uset_u6_load_vectors(model, dof_map, ndof, aset)
    if F_u6.shape[1] > 0:
        all_loads.append(F_u6)

    # Source 3: Inertia loads (RESVINER)
    if include_inertia:
        # Build M*{e_d} at a-set level
        Maa_for_inertia = Maa.toarray() if issparse(Maa) else Maa
        n_directions = min(ndof_per_grid, 6)
        F_inertia = np.zeros((n_a, n_directions), dtype="float64")
        for d in range(n_directions):
            e_d = np.zeros(n_a, dtype="float64")
            e_d[d::ndof_per_grid] = 1.0
            F_inertia[:, d] = Maa_for_inertia @ e_d
        all_loads.append(F_inertia)

    if len(all_loads) == 0:
        return phi_a, eigenvalues

    load_vectors = np.hstack(all_loads)

    # Compute residual vectors
    psi, omega2_pseudo = compute_residual_vectors(
        Kaa, Maa, phi_a, eigenvalues, load_vectors, tiny=tiny
    )

    if psi.shape[1] == 0:
        return phi_a, eigenvalues

    # Augment
    phi_augmented = np.hstack([phi_a, psi])
    eigenvalues_augmented = np.concatenate([eigenvalues, omega2_pseudo])
    return phi_augmented, eigenvalues_augmented
