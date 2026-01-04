import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase


def slice_modal_set(node_gridtype: np.ndarray,
                    phi: np.ndarray,
                    nnode: int, nmode: int,
                    node_set: np.ndarray) -> np.ndarray:
    assert phi.ndim == 2, phi.shape
    assert node_gridtype.shape == (nnode, 2), node_gridtype.shape
    assert phi.shape == (nmode, nnode*6), (phi.shape, (nmode, nnode*6))
    if node_set[0] != 0:  # 0=all
        assert len(np.unique(node_set)), len(node_set)
        # assert phi.shape == (nnode, nmode), phi.shape
        inode = np.searchsorted(node_gridtype[:, 0], node_set)
        assert np.array_equal(node_gridtype[inode, 0], node_set)
        node_gridtype = node_gridtype[inode, :]
        phi = phi.reshape(nmode, nnode, 6)[:, inode, :]
        nnode = len(node_set)
    assert phi.shape == (nmode, nnode, 6)
    return node_gridtype, phi, nnode


def get_real_eigenvalue_method(model: BDF,
                               subcase: Subcase) -> [int, str]:
    method_id, options = subcase['METHOD']
    assert isinstance(method_id, int), method_id
    method = model.methods[method_id]
    if method.type == 'EIGRL':
        neigenvalue = method.nd  # nroots
        norm_str = 'MASS' if method.norm is None else method.norm
    else:
        raise RuntimeError(method)
    #neigenvalues = 10
    return neigenvalue, norm_str


def apply_phi_normalization(Mgg: np.ndarray,
                            Kgg: np.ndarray,
                            eigenvalue: np.ndarray,
                            phit: np.ndarray,
                            nmode: int,
                            nnode_g: int,
                            norm_str: str):
    """
    Applies eigenvalue/eigevenctor normalization
    for norm_str:
     - MASS
     - MAX
    """
    # print(phit)
    assert phit.shape == (nmode, nnode_g*6)
    # print('phit = ', phit)
    # phit6 = phit.reshape(nmode, nnode_g, 6)
    # print(f'norm_max1 = {np.max(np.linalg.norm(phit6,axis=2), axis=1)}')
    # phit *= 2
    # phit6 = phit.reshape(nmode, nnode_g, 6)
    # print(f'norm_max2 = {np.max(np.linalg.norm(phit6,axis=2), axis=1)}')
    if norm_str == 'MAX':
        phit6 = phit.reshape(nmode, nnode_g, 6)

        normi = np.linalg.norm(phit6, axis=2)
        assert normi.shape == (nmode, nnode_g), normi.shape
        norm_scale = normi.max(axis=1)
        assert len(norm_scale) == nmode, (len(norm_scale), nmode)
        phit /= norm_scale[:, np.newaxis]
        eigenvalue /= norm_scale**2

    elif norm_str == 'MASS':
        phi = phit.T
        Mhh = phit @ Mgg @ phi
        # Khh = phit @ Kgg @ phi
        # print(f'eig(Khh) = {np.diag(Khh)}')

        # print(Mhh)
        massh = np.diag(Mhh)
        if not np.allclose(massh, np.ones(nmode)):
            norm_scale = np.sqrt(massh)
            phit /= norm_scale[:, np.newaxis]
            eigenvalue /= massh
    else:
        raise RuntimeError(f'norm_str={norm_str!r} and must be [MASS, MAX]')
    phi = phit.T
    Mhh = phit @ Mgg @ phi
    Khh = phit @ Kgg @ phi
    return phit, Mhh, Khh
