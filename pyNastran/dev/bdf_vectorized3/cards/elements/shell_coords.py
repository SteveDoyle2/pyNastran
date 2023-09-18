import numpy as np


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
        jmat = np.cross(normal, i) # k x i
        try:
            jmat /= np.linalg.norm(jmat)
        except FloatingPointError:
            raise ValueError(f'Cannot project i-axis onto element normal i={i} normal={normal}\n{element}')
        # we do an extra normalization here because
        # we had to project i onto the elemental plane
        # unlike in the next block
        imat = np.cross(jmat, normal)

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
            imat, jmat = rotate_by_thetad(theta[itheta2], imat[itheta2, :], jmat[itheta2, :], normal[itheta, :])
    else:
        raise RuntimeError(element.theta_mcid)
    return imat, jmat


def element_coordinate_system(element,
                              normal: np.ndarray,
                              xyz1: np.ndarray,
                              xyz2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """helper function for material_coordinate_system"""
    imat = xyz2 - xyz1
    imat /= np.linalg.norm(imat, axis=1)
    jmat = np.cross(normal, imat) # k x i
    try:
        jmat /= np.linalg.norm(jmat)
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
