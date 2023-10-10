from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def _nshell_elements(model: BDF) -> tuple[int, str]:
    ncards = 0
    idtype = 'int32'
    for elem in model.shell_element_cards:
        ncards += elem.n
        if elem.element_id.dtype.name == 'int64':
            idtype = 'int64'
    return ncards, idtype


def get_shell_element_coordinate_system(model: BDF) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                                             np.ndarray, np.ndarray]:
    ncards, dtype = _nshell_elements(model)
    if ncards == 0:
        element_id = np.array([], dtype='int32')
        length = np.array([], dtype='int32')
        centroid = np.zeros((0, 3), dtype='float32')
        ielement = np.zeros((0, 3), dtype='float32')
        jelement = np.zeros((0, 3), dtype='float32')
        return element_id, length, centroid, ielement, jelement

    element_id = np.zeros(ncards, dtype=dtype)
    length = np.full(ncards, np.nan, dtype='float32')
    centroid = np.full((ncards, 3), np.nan, dtype='float32')
    ielement = np.full((ncards, 3), np.nan, dtype='float32')
    jelement = np.full((ncards, 3), np.nan, dtype='float32')
    #normal = np.full((n, 3), np.nan, dtype='float32')
    n1 = 0
    for elem in model.shell_element_cards:
        if elem.n == 0:
            continue
        n2 = n1 + elem.n
        dxyzi, centroidi, ielementi, jelementi, normali = elem.element_coordinate_system()
        element_id[n1:n2] = elem.element_id
        length[n1:n2] = dxyzi
        centroid[n1:n2, :] = centroidi
        ielement[n1:n2, :] = ielementi
        jelement[n1:n2, :] = jelementi
        #normal[n1:n2, :] = normali
        n1 = n2
    return element_id, length, centroid, ielement, jelement

def get_shell_material_coordinate_system(model: BDF) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                                              np.ndarray, np.ndarray]:
    ncards, dtype = _nshell_elements(model)
    if ncards == 0:
        element_id = np.array([], dtype='int32')
        length = np.array([], dtype='int32')
        centroid = np.zeros((0, 3), dtype='float32')
        ielement = np.zeros((0, 3), dtype='float32')
        jelement = np.zeros((0, 3), dtype='float32')
        return element_id, length, centroid, ielement, jelement

    element_id = np.zeros(ncards, dtype=dtype)
    length = np.full(ncards, np.nan, dtype='float32')
    centroid = np.full((ncards, 3), np.nan, dtype='float32')
    ielement = np.full((ncards, 3), np.nan, dtype='float32')
    jelement = np.full((ncards, 3), np.nan, dtype='float32')
    #normal = np.full((ncards, 3), np.nan, dtype='float32')
    n1 = 0
    for elem in model.shell_element_cards:
        if elem.n == 0:
            continue
        n2 = n1 + elem.n
        dxyzi, centroidi, ielementi, jelementi, normali = elem.material_coordinate_system()
        eids = elem.element_id
        #assert eids.min() > 0, elem.get_stats()

        neid = len(eids)
        assert neid == elem.n
        element_id[n1:n2] = eids
        length[n1:n2] = dxyzi
        centroid[n1:n2, :] = centroidi
        ielement[n1:n2, :] = ielementi
        jelement[n1:n2, :] = jelementi
        #normal[n1:n2, :] = normali
        n1 = n2
    assert element_id.min() > 0, element_id.min()
    return element_id, length, centroid, ielement, jelement

def material_coordinate_system(element,
                               normal: np.ndarray,
                               xyz1: np.ndarray,
                               xyz2: np.ndarray) -> tuple[np.ndarray,
                                                          np.ndarray]:
    """helper function for material_coordinate_system"""
    #is_theta
    theta = element.theta
    mcid = element.mcid
    itheta = (mcid == -1)
    ntheta = itheta.sum()
    nmcid = len(itheta) - ntheta

    #if element.theta_mcid is None:
        #raise NotImplementedError('theta_mcid=%r' % element.theta_mcid)

    if (ntheta + nmcid == 0):
        raise NotImplementedError((ntheta, nmcid))
    elif ntheta and nmcid:
        raise NotImplementedError((ntheta, nmcid))

    if nmcid:
        #assert element.theta_mcid_ref is not None, f'mcid={element.theta_mcid} not found for\n{element}'
        i = element.theta_mcid_ref.i
        jmat = np.cross(normal, i, axis=1) # k x i
        jnorm = np.linalg.norm(jmat, axis=1)
        try:
            jmat /= jnorm[:, np.newaxis]
        except FloatingPointError:
            raise ValueError(f'Cannot project i-axis onto element normal i={i} normal={normal}\n{element}')
        # we do an extra normalization here because
        # we had to project i onto the elemental plane
        # unlike in the next block
        imat = np.cross(jmat, normal, axis=1)

    if ntheta:
        # rotate by the angle theta
        thetai = theta[itheta]
        #imati = imat[itheta, :, :]
        #jmati = jmat[itheta, :, :]
        normali = normal[itheta, :]
        xyz1i = xyz1[itheta, :]
        xyz2i = xyz2[itheta, :]

        imat, jmat = element_coordinate_system(element, normali, xyz1i, xyz2i)
        if np.all(thetai == 0.):
            pass
        else:
            itheta2 = itheta & (thetai != 0.)
            if itheta2.sum():
                imat, jmat = rotate_by_thetad(theta[itheta2], imat[itheta2, :], jmat[itheta2, :], normal[itheta, :])
    else:
        raise RuntimeError(element.get_stats())
    return imat, jmat


def element_coordinate_system(element,
                              normal: np.ndarray,
                              xyz1: np.ndarray,
                              xyz2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """helper function for material_coordinate_system"""
    imat = xyz2 - xyz1
    normi = np.linalg.norm(imat, axis=1)
    imat /= normi[:, np.newaxis]
    jmat = np.cross(normal, imat, axis=1) # k x i
    jnorm = np.linalg.norm(jmat, axis=1)
    try:
        jmat /= jnorm[:, np.newaxis]
    except FloatingPointError:
        raise ValueError(f'Cannot project i-axis onto element normal i={imat} normal={normal}\n{element}')
    assert xyz1.shape == imat.shape
    assert xyz1.shape == jmat.shape
    return imat, jmat

def rotate_by_thetad(thetad: np.ndarray,
                     imat: np.ndarray,
                     jmat: np.ndarray,
                     normal: np.ndarray):
    theta = np.radians(thetad)
    cos = np.cos(theta)
    sin = np.sin(theta)

    # (n, 3, 3)
    ncoord = len(theta)
    theta_rotation = np.zeros((ncoord, 3, 3), dtype='float64')
    theta_rotation[:, 0, 0] = cos
    theta_rotation[:, 1, 1] = cos
    theta_rotation[:, 0, 1] = sin
    theta_rotation[:, 1, 0] = -sin
    theta_rotation[:, 2, 2] = 1
    #theta_rotation = np.array([
        #[cos, sin, 0.],
        #[-sin, cos, 0.],
        #[0., 0., 1.],
    #], dtype='float64')

    element_axes = np.dstack([imat, jmat, normal])
    assert element_axes.shape == theta_rotation.shape

    rotated_axes = theta_rotation @ element_axes
    assert rotated_axes.shape == theta_rotation.shape
    imat2 = rotated_axes[:, 0, :]
    jmat2 = rotated_axes[:, 1, :]
    return imat2, jmat2
