import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF
import numpy as np
from scipy.sparse import dok_matrix

from ..elements.shells import (
    build_KDgg_cquad4, build_KDgg_ctria3)
from ..elements.solids import build_KDgg_solids
DOF_MAP = dict[tuple[int, int], int]


def build_KDgg(model: BDF, ndof: int) -> dok_matrix:
    KDgg = dok_matrix((ndof, ndof))
    # no KDgg for crod, ctube, conrod
    build_KDgg_bar_beam(model, KDgg, dof_map, u_global)
    build_KDgg_cquad4(model, KDgg, dof_map, u_global)
    build_KDgg_ctria3(model, KDgg, dof_map, u_global)
    build_KDgg_solids(model, KDgg, dof_map, u_global)
    return KDgg


def build_KDgg_bar_beam(model: BDF,
                        KDbb: dok_matrix,
                        dof_map: DOF_MAP,
                        xb: np.ndarray,
                        fdtype: str = "float64",) -> int:
    """Build geometric stiffness for CBAR and CBEAM from axial forces.

    Only beam elements (CBAR/CBEAM) contribute. Rod elements (CROD/CONROD/CTUBE)
    have no bending DOFs and cannot buckle — they are skipped.

    Parameters
    ----------
    model : BDF
        The model.
    KDbb : dok_matrix or _COOAccumulator
        The geometric stiffness matrix to accumulate into.
    dof_map : DOF_MAP
        DOF map (nid, comp) -> index.
    xb : np.ndarray
        Displacement solution from the linear static step.
    fdtype : str
        Float dtype.

    Returns
    -------
    nelements : int
        Number of elements processed.
    """
    ndof = len(xb)
    coo = _COOAccumulator(ndof)

    total = 0
    for elem in [model.cbar, model.cbeam]:
        if elem.n == 0:
            continue

        area = elem.area()
        inertia = elem.inertia()
        xyz1, xyz2 = elem.get_xyz()
        lengths = np.linalg.norm(xyz2 - xyz1, axis=1)
        v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
        k_arr = elem.k()
        e_g_nus = elem.e_g_nu()

        for (nid1, nid2), areai, inertiai, lengthi, ki, e_g_nu, ihati, jhati, khati in zip(
            elem.nodes, area, inertia, lengths, k_arr, e_g_nus,
            ihat, yhat, zhat):
            i1, i2, i12, j = inertiai
            e, g, nu = e_g_nu
            k1, k2 = ki

            Ke_local = timoshenko_stiffness(areai, e, g, lengthi, i1, i2, j, k1, k2)
            Teb = beam_transform(ihati, jhati, khati)

            gi1 = dof_map[(nid1, 1)]
            gi2 = dof_map[(nid2, 1)]
            q_basic = np.hstack([xb[gi1:gi1 + 6], xb[gi2:gi2 + 6]])
            q_elem = Teb @ q_basic
            Fe = Ke_local @ q_elem
            P = -Fe[0]  # internal axial force (negative = compression)

            KDe = geometric_stiffness(areai, e, g, lengthi, i1, i2, j, k1, k2, P)
            KD_basic = Teb.T @ KDe @ Teb

            n_ijv = [
                gi1, gi1 + 1, gi1 + 2, gi1 + 3, gi1 + 4, gi1 + 5,
                gi2, gi2 + 1, gi2 + 2, gi2 + 3, gi2 + 4, gi2 + 5,
            ]
            coo.add_matrix(n_ijv, KD_basic)
        total += elem.n

    # Merge into KDbb (supports dok_matrix or _COOAccumulator)
    KD_csc = coo.to_csc()
    if hasattr(KDbb, 'add_matrix'):
        # KDbb is a _COOAccumulator
        KDbb._rows.extend(coo._rows)
        KDbb._cols.extend(coo._cols)
        KDbb._vals.extend(coo._vals)
    else:
        # KDbb is a dok_matrix — add CSC entries
        KD_coo = KD_csc.tocoo()
        for r, c, v in zip(KD_coo.row, KD_coo.col, KD_coo.data):
            KDbb[r, c] += v
    return total
