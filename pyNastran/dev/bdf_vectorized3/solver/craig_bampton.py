"""
Craig-Bampton (SOL 31 / GEN CB MODEL) reduction.

Computes the Craig-Bampton reduced stiffness (KXX), mass (MXX),
transformation matrix (PHIXA), and Output Transformation Matrices (OTMs).

References
----------
Craig, R.R. & Bampton, M.C.C. (1968). "Coupling of Substructures for
    Dynamic Analysis", AIAA Journal, Vol. 6, No. 7, pp. 1313-1319.
MYSTRAN Source: LK6/ (SOLVE_DLR, CALC_KRRcb, CALC_MRRcb, CALC_MRN, etc.)
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import factorized, eigsh

from cpylog import SimpleLogger

from pyNastran.op4.op4 import OP4

DOF_MAP = dict[tuple[int, int], int]


def run_craig_bampton(
    Kgg: Any,
    Mgg: Any,
    dof_map: DOF_MAP,
    r_set_dofs: list[tuple[int, int]],
    neigenvalues: int = 20,
    log: Any = None,
) -> dict[str, np.ndarray]:
    """Perform Craig-Bampton model reduction.

    Parameters
    ----------
    Kgg : sparse matrix (ndof, ndof)
        Global stiffness matrix (g-set, after SPC reduction to free set).
    Mgg : sparse matrix (ndof, ndof)
        Global mass matrix (same DOF ordering as Kgg).
    dof_map : DOF_MAP
        Maps (nid, dof) to global DOF index.
    r_set_dofs : list of (nid, dof)
        Boundary (retained) DOFs defined by SUPORT cards.
    neigenvalues : int
        Number of fixed-interface normal modes to compute.
    log : logger, optional
        Logging object.

    Returns
    -------
    result : dict with keys:
        'KXX' : (num_cb, num_cb) CB reduced stiffness matrix
        'MXX' : (num_cb, num_cb) CB reduced mass matrix
        'PHIXA' : (ndofa, num_cb) CB transformation matrix
        'DLR' : (ndofl, ndofr) constraint modes
        'eigenvalues' : (nvec,) fixed-interface eigenvalues (rad/s)^2
        'eigen_vec' : (ndofl, nvec) fixed-interface mode shapes
        'KRRcb' : (ndofr, ndofr) CB boundary stiffness
        'MRRcb' : (ndofr, ndofr) CB boundary mass
        'MRN' : (ndofr, nvec) coupling mass
        'IF_LTM' : (ndofr, num_cb) interface force LTM
        'r_dof_indices' : (ndofr,) R-set DOF indices in the A-set
        'l_dof_indices' : (ndofl,) L-set DOF indices in the A-set
    """
    if log is None:
        log = SimpleLogger(level="debug")

    ndof = Kgg.shape[0]

    # Identify R-set (boundary) and L-set (interior) DOF indices
    r_indices = np.array([dof_map[dof] for dof in r_set_dofs], dtype="int32")
    all_indices = np.arange(ndof, dtype="int32")
    l_indices = np.setdiff1d(all_indices, r_indices)

    ndofr = len(r_indices)
    ndofl = len(l_indices)
    log.info(f"Craig-Bampton: ndofr={ndofr}, ndofl={ndofl}, neig={neigenvalues}")

    # Partition stiffness and mass into R/L sets
    Kgg_csc = csc_matrix(Kgg)
    Mgg_csc = csc_matrix(Mgg)

    KLL = Kgg_csc[np.ix_(l_indices, l_indices)]
    KRL = Kgg_csc[np.ix_(r_indices, l_indices)]
    KRR = Kgg_csc[np.ix_(r_indices, r_indices)]

    MLL = Mgg_csc[np.ix_(l_indices, l_indices)]
    MRL = Mgg_csc[np.ix_(r_indices, l_indices)]
    MRR = Mgg_csc[np.ix_(r_indices, r_indices)]

    # =========================================================================
    # Step 1: Solve for constraint modes (DLR)
    # KLL * DLR = -KRL^T  =>  DLR = -KLL^{-1} * KRL^T
    # =========================================================================
    log.debug("  solving for constraint modes (DLR)...")
    KLR = KRL.T.toarray()  # (ndofl, ndofr)

    solve_KLL = factorized(KLL.tocsc())
    DLR = np.zeros((ndofl, ndofr))
    for j in range(ndofr):
        DLR[:, j] = solve_KLL(-KLR[:, j])

    log.debug(f"  DLR shape: {DLR.shape}")

    # =========================================================================
    # Step 2: Fixed-interface normal modes
    # KLL * phi = lambda * MLL * phi
    # =========================================================================
    log.debug("  solving fixed-interface eigenvalue problem...")
    nvec = min(neigenvalues, ndofl - 1)

    eigenvalues, eigen_vec = eigsh(
        KLL.tocsc(),
        k=nvec,
        M=MLL.tocsc(),
        sigma=0.0,
        which="LM",
    )

    # Sort by ascending eigenvalue
    sort_idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[sort_idx]
    eigen_vec = eigen_vec[:, sort_idx]

    # Normalize to unit generalized mass
    for i in range(nvec):
        gen_mass = eigen_vec[:, i] @ MLL @ eigen_vec[:, i]
        if gen_mass > 0:
            eigen_vec[:, i] /= np.sqrt(gen_mass)

    gen_mass = np.array([eigen_vec[:, i] @ MLL @ eigen_vec[:, i] for i in range(nvec)])

    log.info(f"  eigenvalues (Hz): {np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)}")

    # =========================================================================
    # Step 3: CB Stiffness (KXX)
    # KRRcb = KRR + KRL * DLR
    # KXX = [
    #    [KRRcb, 0              ],
    #    [0,     diag(omega_i^2)]]
    # =========================================================================
    KRR_dense = KRR.toarray()
    KRL_dense = KRL.toarray()
    KRRcb = KRR_dense + KRL_dense @ DLR

    Kee = np.diag(eigenvalues)  # omega_i^2 * m_i, but m_i=1 after normalization

    num_cb = ndofr + nvec
    KXX = np.zeros((num_cb, num_cb))
    KXX[:ndofr, :ndofr] = KRRcb
    KXX[ndofr:, ndofr:] = Kee

    # =========================================================================
    # Step 4: CB Mass (MXX)
    # MRRcb = MRR + MRL*DLR + (MRL*DLR)^T + DLR^T * MLL * DLR
    # MRN = (MRL + DLR^T * MLL) * eigen_vec
    # MXX = [[MRRcb, MRN], [MRN^T, I]]
    # =========================================================================
    MRR_dense = MRR.toarray()
    MRL_dense = MRL.toarray()
    MLL_dense = MLL.toarray()

    MRL_DLR = MRL_dense @ DLR
    MRRcb = MRR_dense + MRL_DLR + MRL_DLR.T + DLR.T @ MLL_dense @ DLR

    MRN = (MRL_dense + DLR.T @ MLL_dense) @ eigen_vec

    Mee = np.eye(nvec)  # generalized mass = 1 after normalization

    MXX = np.zeros((num_cb, num_cb))
    MXX[:ndofr, :ndofr] = MRRcb
    MXX[:ndofr, ndofr:] = MRN
    MXX[ndofr:, :ndofr] = MRN.T
    MXX[ndofr:, ndofr:] = Mee

    # =========================================================================
    # Step 5: Transformation matrix (PHIXA)
    # PHIXA = [[I, 0], [DLR, eigen_vec]]
    # Rows correspond to A-set DOFs, columns to CB DOFs
    # =========================================================================
    ndofa = ndofr + ndofl
    PHIXA = np.zeros((ndofa, num_cb))

    # Map R-set rows and L-set rows into A-set ordering
    # A-set ordering: original DOF indices sorted
    a_to_rl = np.zeros(ndofa, dtype="int32")
    a_is_r = np.zeros(ndofa, dtype=bool)

    # Build mapping from a-set position to R/L position
    r_pos_in_a = np.searchsorted(all_indices[:ndofa], r_indices)
    l_pos_in_a = np.searchsorted(all_indices[:ndofa], l_indices)

    for ia, ri in enumerate(r_indices):
        PHIXA[ri, ia] = 1.0  # Identity for R-set

    for il, li in enumerate(l_indices):
        PHIXA[li, :ndofr] = DLR[il, :]  # Constraint modes
        PHIXA[li, ndofr:] = eigen_vec[il, :]  # Normal modes

    # =========================================================================
    # Step 6: Interface Force LTM
    # IF_LTM = [MRRcb | MRN | KRRcb]
    # Columns: [R-set accel | modal accel | R-set disp]
    # =========================================================================
    num_cb_full = 2 * ndofr + nvec
    IF_LTM = np.zeros((ndofr, num_cb_full))
    IF_LTM[:, :ndofr] = MRRcb  # R-set accelerations
    IF_LTM[:, ndofr : ndofr + nvec] = MRN  # modal accelerations
    IF_LTM[:, ndofr + nvec :] = KRRcb  # R-set displacements

    log.info(f"  Craig-Bampton complete: KXX={KXX.shape}, MXX={MXX.shape}")

    return {
        "KXX": KXX,
        "MXX": MXX,
        "PHIXA": PHIXA,
        "DLR": DLR,
        "eigenvalues": eigenvalues,
        "eigen_vec": eigen_vec,
        "gen_mass": gen_mass,
        "KRRcb": KRRcb,
        "MRRcb": MRRcb,
        "MRN": MRN,
        "IF_LTM": IF_LTM,
        "r_dof_indices": r_indices,
        "l_dof_indices": l_indices,
        "num_cb_dofs": num_cb,
    }


def write_cb_to_op4(
    op4_filename: str,
    cb_result: dict[str, np.ndarray],
    is_binary: bool = True,
    precision: str = "double",) -> None:
    """
    Write Craig-Bampton matrices to Nastran OP4 format.

    Writes the following named matrices:
        KXX    - CB reduced stiffness (num_cb x num_cb)
        MXX    - CB reduced mass (num_cb x num_cb)
        PHIXA  - CB transformation matrix (ndofa x num_cb)
        DLR    - Constraint modes (ndofl x ndofr)
        IF_LTM - Interface force load transformation matrix

    Parameters
    ----------
    op4_filename : str
        Output OP4 file path.
    cb_result : dict
        Result dictionary from run_craig_bampton().
    is_binary : bool
        True for binary OP4, False for ASCII.
    precision : str
        'single', 'double', or 'default'.
    """
    # form: 1=square, 2=rectangular, 6=symmetric
    matrices = {}

    KXX = cb_result["KXX"].astype(np.float64)
    MXX = cb_result["MXX"].astype(np.float64)
    PHIXA = cb_result["PHIXA"].astype(np.float64)
    DLR = cb_result["DLR"].astype(np.float64)
    IF_LTM = cb_result["IF_LTM"].astype(np.float64)

    # KXX and MXX are square symmetric
    matrices["KXX"] = (6, KXX)
    matrices["MXX"] = (6, MXX)
    # PHIXA and DLR are rectangular
    matrices["PHIXA"] = (2, PHIXA)
    matrices["DLR"] = (2, DLR)
    matrices["IF_LTM"] = (2, IF_LTM)

    name_order = ["KXX", "MXX", "PHIXA", "DLR", "IF_LTM"]
    op4 = OP4()
    op4.write_op4(
        op4_filename,
        matrices,
        name_order=name_order,
        precision=precision,
        is_binary=is_binary,)


def write_cb_to_h5(h5_filename: str,
                   cb_result: dict[str, np.ndarray],) -> None:
    """Write Craig-Bampton matrices to HDF5 using PyTables.

    Creates the following datasets under /CRAIG_BAMPTON/:
        KXX        - CB reduced stiffness (num_cb x num_cb)
        MXX        - CB reduced mass (num_cb x num_cb)
        PHIXA      - CB transformation matrix (ndofa x num_cb)
        DLR        - Constraint modes (ndofl x ndofr)
        IF_LTM     - Interface force LTM
        KRRcb      - Boundary stiffness (ndofr x ndofr)
        MRRcb      - Boundary mass (ndofr x ndofr)
        MRN        - Coupling mass (ndofr x nvec)
        eigenvalues - Fixed-interface eigenvalues (rad/s)^2
        gen_mass   - Generalized mass per mode

    Parameters
    ----------
    h5_filename : str
        Output HDF5 file path.
    cb_result : dict
        Result dictionary from run_craig_bampton().
    """
    import tables

    with tables.open_file(h5_filename, mode="w", title="Craig-Bampton Model") as h5:
        cb_group = h5.create_group(
            "/", "CRAIG_BAMPTON", "Craig-Bampton Matrices")

        h5.create_array(
            cb_group, "KXX", cb_result["KXX"].astype(np.float64),
            title="CB Reduced Stiffness")
        h5.create_array(
            cb_group, "MXX", cb_result["MXX"].astype(np.float64),
            title="CB Reduced Mass")
        h5.create_array(
            cb_group,
            "PHIXA",
            cb_result["PHIXA"].astype(np.float64),
            title="CB Transformation Matrix",
        )
        h5.create_array(
            cb_group, "DLR", cb_result["DLR"].astype(np.float64),
            title="Constraint Modes"
        )
        h5.create_array(
            cb_group, "IF_LTM", cb_result["IF_LTM"].astype(np.float64),
            title="Interface Force LTM")
        h5.create_array(
            cb_group, "KRRcb", cb_result["KRRcb"].astype(np.float64),
            title="Boundary Stiffness")
        h5.create_array(
            cb_group, "MRRcb", cb_result["MRRcb"].astype(np.float64),
            title="Boundary Mass")
        h5.create_array(cb_group, "MRN", cb_result["MRN"].astype(np.float64),
            title="Coupling Mass")
        h5.create_array(
            cb_group,
            "eigenvalues",
            cb_result["eigenvalues"].astype(np.float64),
            title="Fixed-Interface Eigenvalues (rad/s)^2",)
        h5.create_array(
            cb_group, "gen_mass", cb_result["gen_mass"].astype(np.float64),
            title="Generalized Mass")

        # Store dimension info as attributes
        cb_group._v_attrs.ndofr = len(cb_result["r_dof_indices"])
        cb_group._v_attrs.ndofl = len(cb_result["l_dof_indices"])
        cb_group._v_attrs.nvec = len(cb_result["eigenvalues"])
        cb_group._v_attrs.num_cb_dofs = cb_result["num_cb_dofs"]

        # DOF index arrays for reconstruction
        h5.create_array(
            cb_group,
            "r_dof_indices",
            cb_result["r_dof_indices"].astype(np.int32),
            title="R-set DOF Indices",)
        h5.create_array(
            cb_group,
            "l_dof_indices",
            cb_result["l_dof_indices"].astype(np.int32),
            title="L-set DOF Indices",)
