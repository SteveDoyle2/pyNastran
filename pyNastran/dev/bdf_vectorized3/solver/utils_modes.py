from typing import TextIO, Any
import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
from pyNastran.op2.op2 import OP2
from pyNastran.op2.op2_interface.op2_classes import (
    # RealDisplacementArray,
    # RealSPCForcesArray,
    # RealLoadVectorArray,
    # ComplexDisplacementArray,
    # ComplexVelocityArray,
    # ComplexAccelerationArray,
    RealEigenvalues,
)
def todense(Mgg):
    if hasattr(Mgg, "toarray"):
        Mgg_dense = Mgg.toarray()
    else:
        Mgg_dense = np.asarray(Mgg)
    return Mgg_dense


def slice_modal_set(node_gridtype: np.ndarray,
                    phi: np.ndarray,
                    nnode: int, nmode: int,
                    node_set: np.ndarray) -> np.ndarray:
    assert phi.ndim == 2, phi.shape
    assert node_gridtype.shape == (nnode, 2), node_gridtype.shape
    assert phi.shape == (nmode, nnode * 6), (phi.shape, (nmode, nnode * 6))
    if node_set[0] == 0:  # 0=all
        phi = phi.copy().reshape(nmode, nnode, 6)
    else:  # subset
        assert len(np.unique(node_set)), len(node_set)
        # assert phi.shape == (nnode, nmode), phi.shape
        inode = np.searchsorted(node_gridtype[:, 0], node_set)
        assert np.array_equal(node_gridtype[inode, 0], node_set)
        node_gridtype = node_gridtype[inode, :]
        phi = phi.reshape(nmode, nnode, 6)[:, inode, :]
        nnode = len(node_set)
    assert phi.shape == (nmode, nnode, 6), phi.shape
    return node_gridtype, phi, nnode


def get_real_eigenvalue_method(model: BDF,
                               subcase: Subcase) -> [int, str]:
    method_id, options = subcase["METHOD"]
    assert isinstance(method_id, int), method_id
    method = model.methods[method_id]
    if method.type == "EIGRL":
        neigenvalue = method.nd  # nroots
        norm_str = "MASS" if method.norm is None else method.norm
    elif method.type == 'EIGB':
        # C      : None
        # G      : None
        # L1     : 0.0
        # L2     : 100.0
        # method : 'INV'
        # ndn    : 60
        # ndp    : 60
        # nep    : 20
        neigenvalue = method.ndn + method.ndp
        norm_str = method.norm
    else:
        raise RuntimeError(method.get_stats())
    # neigenvalues = 10
    assert isinstance(norm_str, str), norm_str
    assert norm_str in ['MAX', 'MASS', 'POINT'], norm_str
    return neigenvalue, norm_str


def apply_phi_normalization(Mgg: np.ndarray,
                            Kgg: np.ndarray,
                            eigenvalue: np.ndarray,
                            phit: np.ndarray,
                            nmode: int,
                            nnode_g: int,
                            norm_str: str,
                            ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Applies eigenvalue/eigevenctor normalization
    for norm_str:
     - MASS
     - MAX
    """
    # print(phit)
    assert phit.shape == (nmode, nnode_g * 6)
    # print('phit = ', phit)
    # phit6 = phit.reshape(nmode, nnode_g, 6)
    # print(f'norm_max1 = {np.max(np.linalg.norm(phit6,axis=2), axis=1)}')
    # phit *= 2
    # phit6 = phit.reshape(nmode, nnode_g, 6)
    # print(f'norm_max2 = {np.max(np.linalg.norm(phit6,axis=2), axis=1)}')
    if norm_str == "MAX":
        phit6 = phit.reshape(nmode, nnode_g, 6)

        normi = np.linalg.norm(phit6, axis=2)
        assert normi.shape == (nmode, nnode_g), normi.shape
        norm_scale = normi.max(axis=1)
        assert len(norm_scale) == nmode, (len(norm_scale), nmode)
        phit /= norm_scale[:, np.newaxis]
        eigenvalue /= norm_scale**2

    elif norm_str == "MASS":
        phi = phit.T
        Mhh = phit @ Mgg @ phi
        # Khh = phit @ Kgg @ phi
        # print(f'eig(Khh) = {np.diag(Khh)}')

        # print(Mhh)
        mhh_diag = np.diag(Mhh)
        if not np.allclose(mhh_diag, np.ones(nmode)):
            norm_scale = np.sqrt(mhh_diag)
            phit /= norm_scale[:, np.newaxis]
            eigenvalue /= mhh_diag
    else:
        raise RuntimeError(f"norm_str={norm_str!r} and must be [MASS, MAX]")
    phi = phit.T
    Mhh = phit @ Mgg @ phi
    Khh = phit @ Kgg @ phi
    return phit, Mhh, Khh


def save_eigenvalues(op2: OP2,
                     f06_file: TextIO,
                     out: dict[str, Any],
                     subcase: Subcase,
                     title: str,
                     page_stamp: str = '',
                     page_num: int=1) -> int:
    eigenvalue = out["modes_eigenvalue"]
    Mhh = out["modes_Mhh"]
    Khh = out["modes_Khh"]
    nmode = len(eigenvalue)

    eigenvalue_obj = RealEigenvalues(title, "LAMA", nmodes=nmode)

    cycle = np.sqrt(np.abs(eigenvalue)) / (2.0 * np.pi)
    radian = np.sqrt(np.abs(eigenvalue))

    eigenvalue_obj.mode = np.arange(nmode, dtype="int32") + 1
    eigenvalue_obj.extraction_order = np.arange(nmode, dtype="int32") + 1
    eigenvalue_obj.eigenvalues = eigenvalue
    eigenvalue_obj.radians = radian
    eigenvalue_obj.cycles = cycle
    eigenvalue_obj.generalized_mass = np.diag(Mhh)
    eigenvalue_obj.generalized_stiffness = np.diag(Khh)

    op2.eigenvalues[title] = eigenvalue_obj
    # op2.eigenvalues[isubcase] = eigenvalues_obj

    write_eigenvalue_f06 = True
    if write_eigenvalue_f06:
        str(page_num)
        header = []
        page_num = eigenvalue_obj.write_f06(
            f06_file, header=header, page_stamp=page_stamp, page_num=page_num)
        f06_file.write("\n")
    return page_num

def compute_mass_participation(phig: np.ndarray,
                               Mgg: np.ndarray,
                               nnode_g: int, nmode: int,
                               ) -> dict[str, np.ndarray]:
    """Compute modal mass participation factors and effective mass.

    Assumes mass-normalized eigenvectors (phi.T @ M @ phi = I).

    Parameters
    ----------
    phig : ndarray, shape (nmode, nnode_g*6)
        Mass-normalized eigenvectors (rows = modes).
    Mgg : ndarray, shape (nnode_g*6, nnode_g*6)
        Global mass matrix.
    nnode_g : int
        Number of grid points in the g-set.
    nmode : int
        Number of modes.

    Returns
    -------
    dict with keys:
        participation : ndarray, shape (nmode, 6)
            Modal participation factors per direction (Tx, Ty, Tz, Rx, Ry, Rz).
        effective_mass : ndarray, shape (nmode, 6)
            Effective mass per mode per direction.
        effective_mass_ratio : ndarray, shape (nmode, 6)
            Effective mass as fraction of total mass per direction.
        total_mass : ndarray, shape (6,)
            Total mass per direction (diagonal of r.T @ M @ r).
        cumulative_ratio : ndarray, shape (nmode, 6)
            Cumulative effective mass ratio.
    """
    ndof = nnode_g * 6
    assert phig.shape == (nmode, ndof), (phig.shape, (nmode, ndof))

    # Build rigid-body direction vectors:
    #   one per DOF direction (Tx,Ty,Tz,Rx,Ry,Rz)
    # Each column of R selects DOF i at every node
    R = np.zeros((ndof, 6), dtype=phig.dtype)
    for i in range(6):
        R[i::6, i] = 1.0

    # Total mass per direction: r_j.T @ M @ r_j
    Mgg_dense = todense(Mgg)
    total_mass = np.einsum("ij,ik,jk->k", R, Mgg_dense @ R, np.eye(6))
    # Simpler: total_mass[j] = R[:,j].T @ M @ R[:,j]
    total_mass = np.array([R[:, j] @ Mgg_dense @ R[:, j] for j in range(6)])

    # Participation factors: gamma_ij = phi_i @ M @ r_j
    # phig is (nmode, ndof), Mgg is (ndof, ndof), R is (ndof, 6)
    MR = Mgg_dense @ R  # (ndof, 6)
    participation = phig @ MR  # (nmode, 6)

    effective_mass = participation**2

    effective_mass_ratio = np.zeros_like(effective_mass)
    nonzero = total_mass > 0.0
    effective_mass_ratio[:, nonzero] = effective_mass[:, nonzero] / total_mass[nonzero]

    cumulative_ratio = np.cumsum(effective_mass_ratio, axis=0)

    return {
        "participation": participation,
        "effective_mass": effective_mass,
        "effective_mass_ratio": effective_mass_ratio,
        "total_mass": total_mass,
        "cumulative_ratio": cumulative_ratio,
    }
