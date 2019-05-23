"""
Defines general coordinate system related functions including:
 - xyz_to_rtz_array(xyz)
 - xyz_to_rtp_array(xyz)

 - rtz_to_xyz_array(xyz)
 - rtp_to_xyz_array(xyz)

 - rtz_to_rtp_array(xyz)
 - rtp_to_rtz_array(xyz)

 - coords = cylindrical_rotation_matrix(thetar, dtype='float64')

"""
# pylint: disable=C0103
import numpy as np

# xyz to xxx transforms
def xyz_to_rtz_array(xyz):
    """
    xyz to R-theta-z transform::

      y       R
      |     /
      |   /
      | / theta
      *------------x

    .. math:: x = R \cos(\theta)
    .. math:: y = R \sin(\theta)

    Returns
    -------
    xyz : (3,) float ndarray
        the point in the local coordinate system

    """
    xyz = np.atleast_2d(xyz)
    assert len(xyz.shape) == 2, xyz.shape
    x = xyz[:, 0]
    y = xyz[:, 1]
    theta = np.degrees(np.arctan2(y, x))
    R = np.sqrt(x * x + y * y)
    return np.array([R, theta, xyz[:, 2]], dtype=xyz.dtype).T

def xyz_to_rtp_array(xyz):
    """rho-theta-phi to xyz transform"""
    xyz = np.atleast_2d(xyz)
    assert len(xyz.shape) == 2, xyz.shape
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]
    rho = np.sqrt(x * x + y * y + z * z)
    phi = np.degrees(np.arctan2(y, x))

    #i = np.where(rho == 0.0)
    #if len(i):
    theta = np.zeros(len(z), dtype=z.dtype)
    ir = np.where(rho != 0.0)
    theta[ir] = np.degrees(np.arccos(z[ir] / rho[ir]))
    return np.array([rho, theta, phi], dtype=xyz.dtype).T

#---------------------------------------------------------------
# xxx to xyz transforms
def rtz_to_xyz_array(rtz):
    r"""
    R-theta-z to xyz transform::

      y       R
      |     /
      |   /
      | / theta
      *------------x

    .. math:: x = R \cos(\theta)
    .. math:: y = R \sin(\theta)

    Returns
    -------
    xyz : (3,) float ndarray
        the point in the local coordinate system

    """
    rtz = np.atleast_2d(rtz)
    assert len(rtz.shape) == 2, rtz.shape
    R = rtz[:, 0]
    theta = np.radians(rtz[:, 1])
    x = R * np.cos(theta)
    y = R * np.sin(theta)
    xyz = np.array([x, y, rtz[:, 2]], dtype=rtz.dtype).T
    return xyz

def rtp_to_xyz_array(rtp):
    """
    rho-theta-phi to xyz transform

    Returns
    -------
    xyz : (3,) float ndarray
        the x, y, z in the local coordinate system

    """
    rtp = np.atleast_2d(rtp)
    assert len(rtp.shape) == 2, rtp.shape
    R = rtp[:, 0]
    theta = np.radians(rtp[:, 1])
    phi = np.radians(rtp[:, 2])
    x = R * np.sin(theta) * np.cos(phi)
    y = R * np.sin(theta) * np.sin(phi)
    z = R * np.cos(theta)
    return np.array([x, y, z], dtype=rtp.dtype).T

#---------------------------------------------------------------
# rtz/rtp and rtp/rtz transforms
def rtz_to_rtp_array(rtz):
    """R-theta-z to rho-theta-phi transform"""
    rtz = np.atleast_2d(rtz)
    r = rtz[:, 0]
    thetad = rtz[:, 1]
    z = rtz[:, 2]

    rho = (r**2 + z**2)**0.5
    irho0 = np.where(rho > 0.0)[0]
    dtype = thetad.dtype

    # We need to choose a default.  The equation for phi is:
    #    phi = acos(z/sqrt(x^2 + y^2 + z^2))
    #
    # If we let x and y go to 0, we're left with z/z=1 and
    # phi = acos(1) = 0.  The other alternative is to let
    # z -> 0 and x/y be non-zero, but that leaves us with
    # a 90 degree angle, which feels wrong.
    #
    phi = np.full(thetad.shape, 0., dtype=thetad.dtype)

    phi[irho0] = np.degrees(np.arccos(z[irho0] / rho[irho0]))
    return np.array([rho, thetad, phi], dtype=dtype).T

def rtp_to_rtz_array(rtp):
    """rho-theta-phi to R-theta-z transform"""
    rtp = np.atleast_2d(rtp)
    rho = rtp[:, 0]
    thetad = rtp[:, 1]
    phid = rtp[:, 2]
    phi = np.radians(phid)
    r = rho * np.sin(phi)
    z = rho * np.cos(phi)
    return np.array([r, thetad, z], dtype=rtp.dtype).T

def cylindrical_rotation_matrix(thetar, dtype='float64'):
    """
    Creates a series transformation matrices to rotate by some angle theta

    Parameters
    ----------
    thetar : (n, ) float ndarray
        the theta in radians
    dtype : dtype/str
        the type of the output matrix

    Returns
    -------
    rotation : (ntheta, 3, 3)
        the rotation matrices
    """
    theta = np.asarray(thetar, dtype=dtype)
    ntheta = len(theta)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    #rotation = np.array([
        #[cos(theta), -sin(theta), zero],
        #[sin(theta),  cos(theta), zero],
        #[0.zero,       zero,       one],
    #], dtype=dtype)

    # can this be faster?
    # it definitely could look nicer
    rotation = np.zeros((ntheta, 3, 3), dtype=dtype)
    rotation[:, 0, 0] = cos_theta
    rotation[:, 0, 1] = -sin_theta
    rotation[:, 1, 1] = cos_theta
    rotation[:, 1, 0] = sin_theta
    rotation[:, 2, 2] = 1.

    #print('---------')
    #for rot in rotation:
        #print(np.squeeze(rot))
        #print('---------')
    return rotation
