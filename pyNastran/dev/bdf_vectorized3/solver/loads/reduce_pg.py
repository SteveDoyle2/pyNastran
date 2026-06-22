"""Reduce the load vector PG from g-set to a-set.

Implements the full Nastran DOF reduction chain for load vectors:

    PG (g-set)
      | MPC elimination: Pn = Gm^T @ Pg_m + Pg_n
    PN (n-set)
      | SPC removal: Pf = Pn_f (partition)
    PF (f-set)
      | OMIT/Guyan: Pa = Pf_a + Kao @ Koo^{-1} @ Pf_o
    PA (a-set)

References
----------
- Simcenter Nastran DMAP Programmer's Guide, "Static Solution Sequence" flow.
- NX Nastran User's Guide, Appendix B: Set Notation System.
"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from scipy.sparse import issparse, csc_matrix
from scipy.sparse.linalg import spsolve
if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def reduce_Pg_to_Pa(
    Pg: np.ndarray,
    ndof: int,
    sset_b: np.ndarray,
    mset_b: np.ndarray | None = None,
    oset_b: np.ndarray | None = None,
    Gm: np.ndarray | None = None,
    Koo: np.ndarray | csc_matrix | None = None,
    Koa: np.ndarray | csc_matrix | None = None,) -> np.ndarray:
    """Reduce load vector from g-set to a-set through MPC, SPC, and OMIT.

    Parameters
    ----------
    Pg : ndarray, shape (ndof,) or (ndof, n_loads)
        Load vector(s) at g-set DOFs.
    ndof : int
        Total number of g-set DOFs.
    sset_b : ndarray of bool, shape (ndof,)
        Boolean mask for SPC DOFs (True = constrained by SPC).
    mset_b : ndarray of bool or None
        Boolean mask for MPC dependent DOFs (True = dependent/eliminated).
        If None, MPC step is skipped.
    oset_b : ndarray of bool or None
        Boolean mask for omitted DOFs (True = interior/condensed out).
        If None, OMIT step is skipped.
    Gm : ndarray or None, shape (n_m, n_n)
        MPC constraint transformation matrix relating dependent DOFs (m-set)
        to independent DOFs (n-set). The dependent displacement is:
            u_m = Gm @ u_n
        If None, MPC step is skipped.
    Koo : ndarray or sparse matrix or None, shape (n_o, n_o)
        Stiffness of the omitted (interior) DOFs. Required if oset_b is given.
    Koa : ndarray or sparse matrix or None, shape (n_o, n_a)
        Coupling stiffness between omitted and analysis DOFs. Required if
        oset_b is given.

    Returns
    -------
    Pa : ndarray, shape (n_a,) or (n_a, n_loads)
        Load vector(s) reduced to the analysis set.

    Notes
    -----
    The three reduction steps:

    **MPC elimination (G → N):**
    The MPC equation is: u_m = Gm @ u_n (dependent = transform * independent).
    Force equivalence (virtual work): Pn = Pg_n + Gm^T @ Pg_m.
    This adds the dependent-DOF forces to the independent DOFs via the
    transpose of the constraint matrix.

    **SPC removal (N → F):**
    Simply partition out the constrained DOFs. Forces at SPC DOFs become
    reaction forces (not part of the free load vector).

    **OMIT/Guyan reduction (F → A):**
    Static condensation of interior DOFs onto the boundary:
        Pa = Pf_a - Kao @ Koo^{-1} @ Pf_o
    The interior load is statically distributed to the boundary through
    the stiffness coupling. Note the NEGATIVE sign — this is because the
    Guyan transformation is:
        u_o = -Koo^{-1} @ Koa @ u_a + Koo^{-1} @ Pf_o
    and the equivalent boundary load from virtual work is:
        Pa = Pf_a + Goa^T @ Pf_o
           = Pf_a + (-Koo^{-1} @ Koa)^T @ Pf_o
           = Pf_a - Kao^T @ Koo^{-T} @ Pf_o
        Goa = -Koo^{-1} @ Koa
    For symmetric K:
        Kao^T = Koa
        Koo^{-T} = Koo^{-1}
        giving:
          Pa = Pf_a - Koa^T @ Koo^{-1} @ Pf_o
    But since Koa = Kao^T (symmetry), this is:
        Pa = Pf_a + Kao @ Koo^{-1} @ Pf_o

    The static condensation of loads uses:
        Pa = Pf_a - Kao @ Koo^{-1} @ Pf_o

    This comes from eliminating u_o from:
        [Kaa Kao] [u_a]   [Pf_a]
        [Koa Koo] [u_o] = [Pf_o]

    Kaa ua + Kao uo = Pfa
    Koa ua + Koo uo = Pfo
 
    Multiply by Koo^-1:
      Koo^-1 Kaa ua + Koo^-1 Kao uo = Koo^-1 Pfa
      Koo^-1 Koa ua +            uo = Koo^-1 Pfo
      uo = Koo^-1 Pfo - Koo^-1 Koa ua
         = Koo^-1 Pfo - Goa ua  ***
      Goa = Koo^-1 Koa  ***
    
    Plug in:
      Kaa ua + Kao (Koo^-1 Pfo - Goa ua) = Pfa
      Kaa ua + Kao (Koo^-1 Pfo - Goa ua) = Pfa
      Kaa ua + Kao Koo^-1 Pfo - Kao Goa ua = Pfa
      Kaa ua - Kao Goa ua = Pfa - Kao Koo^-1 Pfo
        (Kaa - Kao Goa) ua = Pa2
        Kaa' = Kaa - Kao Goa
             = Kaa - Kao Koo^-1 Koa
        Kaa' ua = Pa2

        Pa2 = Pfa - Kao Koo^-1 Pfo  ***
        
        Goa    = Koo^-1 Koa  ***

        ua = (Kaa - Kao Goa       )^-1 Pa2
        ua = (Kaa - Kao Koo^-1 Koa)^-1 Pa2
          

    Row 2: u_o = Koo^{-1} @ (Pf_o - Koa @ u_a)
    Sub into row 1: Kaa @ u_a + Kao @ Koo^{-1} @ (Pf_o - Koa @ u_a) = Pf_a
    Rearrange: (Kaa - Kao @ Koo^{-1} @ Koa) @ u_a = Pf_a - Kao @ Koo^{-1} @ Pf_o

    So:  Pa_reduced = Pf_a - Kao @ Koo^{-1} @ Pf_o
    And: Kaa_reduced = Kaa - Kao @ Koo^{-1} @ Koa
    """
    is_1d = Pg.ndim == 1
    if is_1d:
        Pg = Pg.reshape(-1, 1)
    n_loads = Pg.shape[1]

    # --- Step 1: MPC elimination (G → N) ---
    has_mpc = (
        mset_b is not None and
        Gm is not None and
        np.any(mset_b))
    if has_mpc:
        nset_b = ~mset_b  # n = g - m
        Pg_m = Pg[mset_b, :]  # Forces at dependent DOFs
        Pg_n = Pg[nset_b, :]  # Forces at independent DOFs
        # Virtual work equivalence: Pn = Pg_n + Gm^T @ Pg_m
        Pn = Pg_n + Gm.T @ Pg_m
    else:
        nset_b = np.ones(ndof, dtype="bool")
        Pn = Pg.copy()

    # --- Step 2: SPC removal (N → F) ---
    # Within the n-set, find which DOFs are SPC'd
    if has_mpc:
        # sset_b is defined in g-set; map to n-set indices
        n_indices = np.where(nset_b)[0]
        sset_in_n = np.isin(n_indices, np.where(sset_b)[0])
        fset_in_n = ~sset_in_n
        Pf = Pn[fset_in_n, :]
        # Track which g-set DOFs are in f-set
        fset_b = nset_b.copy()
        fset_b[sset_b] = False
    else:
        fset_b = ~sset_b
        Pf = Pn[~sset_b[nset_b], :]

    # --- Step 3: OMIT reduction (F → A) ---
    has_omit = (
        oset_b is not None and
        np.any(oset_b) and
        Koo is not None and
        Koa is not None)
    if has_omit:
        # oset_b is in g-set; find which f-set DOFs are omitted
        f_indices = np.where(fset_b)[0]
        oset_in_f = np.isin(f_indices, np.where(oset_b)[0])
        aset_in_f = ~oset_in_f

        Pf_a = Pf[aset_in_f, :]
        Pf_o = Pf[oset_in_f, :]

        # Pa = Pf_a - Kao @ Koo^{-1} @ Pf_o
        # Note: Kao = Koa^T for symmetric K
        if issparse(Koo):
            # Solve Koo @ x = Pf_o for each load column
            Koo_inv_Pfo = np.zeros_like(Pf_o)
            for i in range(n_loads):
                Koo_inv_Pfo[:, i] = spsolve(csc_matrix(Koo), Pf_o[:, i])
        else:
            Koo_inv_Pfo = np.linalg.solve(Koo, Pf_o)

        # Kao = Koa^T (stiffness is symmetric)
        Kao = Koa.T
        Pa = Pf_a - Kao @ Koo_inv_Pfo
    else:
        Pa = Pf

    if is_1d:
        return Pa.ravel()
    return Pa


def get_mpc_dependent_dofs(model: BDF) -> tuple[np.ndarray,
                                                np.ndarray | None]:
    """Identify MPC-dependent DOFs and build the Gm transformation matrix.

    Parameters
    ----------
    model : BDF
        The pyNastran BDF model.

    Returns
    -------
    mset_b : ndarray of bool, shape (ndof,)
        Boolean mask identifying dependent (m-set) DOFs.
    Gm : ndarray or None, shape (n_m, n_n)
        MPC constraint transformation. None if no MPCs exist.

    Notes
    -----
    For each MPC equation:  a1*u1 + a2*u2 + ... = 0
    The first DOF (a1*u1) is dependent. Solving for u1:
        u1 = -(a2/a1)*u2 - (a3/a1)*u3 - ...

    So Gm has one row per dependent DOF, with columns corresponding
    to the independent (n-set) DOFs. Entry Gm[i,j] = -a_j/a_1 where
    a_1 is the coefficient of the dependent DOF.
    """
    mpc = model.mpc
    if mpc.n == 0:
        return None, None

    # Build dof_map from the model
    ndof = model.grid.n * 6 + model.spoint.n
    dof_map = {}
    i = 0
    for nid in model.grid.node_id:
        for dof in range(1, 7):
            dof_map[(nid, dof)] = i
            i += 1
    for nid in model.spoint.spoint_id:
        dof_map[(nid, 0)] = i
        i += 1

    # Identify dependent DOFs and build Gm
    dependent_dofs = []
    gm_rows = []

    for mpc_idi, (idim0, idim1) in zip(mpc.mpc_id, mpc.idim):
        coefficients = mpc.coefficients[idim0:idim1]
        components = mpc.components[idim0:idim1]
        nodes = mpc.node_id[idim0:idim1]

        # First entry is dependent
        dep_nid = nodes[0]
        dep_comp = components[0]
        dep_coeff = coefficients[0]
        dep_dof_idx = dof_map[(dep_nid, dep_comp)]
        dependent_dofs.append(dep_dof_idx)

        # Build row: u_dep = sum(-a_j/a_dep * u_j) for j > 0
        row = np.zeros(ndof, dtype="float64")
        for j in range(1, len(nodes)):
            ind_dof_idx = dof_map[(nodes[j], components[j])]
            row[ind_dof_idx] = -coefficients[j] / dep_coeff
        gm_rows.append(row)

    if len(dependent_dofs) == 0:
        return None, None

    mset_b = np.zeros(ndof, dtype="bool")
    mset_b[dependent_dofs] = True

    # Gm maps n-set (independent) to m-set (dependent)
    # Columns correspond to n-set DOFs (all DOFs except dependent)
    nset_b = ~mset_b
    Gm_full = np.array(gm_rows)  # (n_m, ndof)
    Gm = Gm_full[:, nset_b]  # (n_m, n_n) — only independent columns

    return mset_b, Gm


def get_omit_dofs(model, dof_map: dict, ndof: int) -> np.ndarray | None:
    """Identify OMIT (o-set) DOFs from OMIT/OMIT1 cards.

    Parameters
    ----------
    model : BDF
        The pyNastran BDF model.
    dof_map : dict
        (node_id, component) -> global DOF index.
    ndof : int
        Total g-set DOF count.

    Returns
    -------
    oset_b : ndarray of bool or None
        Boolean mask for omitted DOFs. None if no OMIT cards.
    """
    omit = model.omit
    omit1 = model.omit1
    if omit.n == 0 and omit1.n == 0:
        return None

    oset_b = np.zeros(ndof, dtype="bool")
    for nid, comp in zip(omit.node_id, omit.component):
        for c in str(comp):
            dof = int(c)
            if (nid, dof) in dof_map:
                oset_b[dof_map[(nid, dof)]] = True

    assert len(omit1.node_id) == len(omit1.component)
    for nid, comp in zip(omit1.node_id, omit1.component):
        for c in str(comp):
            dof = int(c)
            if (nid, dof) in dof_map:
                oset_b[dof_map[(nid, dof)]] = True
    return oset_b
