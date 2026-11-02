import numpy as np


def validate_region_rotations(rotations: list[tuple[list[int], np.ndarray, np.ndarray,
                                                     str, float]],
                              region: np.ndarray,
                              ortho_tol: float=1e-6) -> list[tuple[np.ndarray, np.ndarray,
                                                                   np.ndarray, str, float]]:
    """
    Checks the rotate_regions sets (shared by Fluent and Cart3D)

    Parameters
    ----------
    rotations : list[(region_ids, origin, axes, axis, angle)]
        see Fluent.rotate_regions
    region : (nelements,) int ndarray
        the region of every element in the model
    ortho_tol : float; default=1e-6
        tolerance on axes @ axes.T == I and det(axes) == 1

    Returns
    -------
    sets : list[(region_ids, origin, rotation_matrix, axis, angle)]
        rotation_matrix : (3, 3) float ndarray
            the basic-frame rotation; xyz_new = origin + (xyz - origin) @ R.T
    """
    all_regions = np.unique(region)
    sets = []
    used_regions = {}
    for iset, rotation in enumerate(rotations):
        if len(rotation) != 5:
            raise ValueError(f'set {iset}: expected (region_ids, origin, axes, axis, angle); '
                             f'got {len(rotation)} values')
        region_ids, origin, axes, axis, angle = rotation
        region_ids = np.unique(np.asarray(region_ids).ravel())
        origin = np.asarray(origin, dtype='float64')
        axes = np.asarray(axes, dtype='float64')
        if axis not in ('x', 'y', 'z'):
            raise ValueError(f'set {iset}: axis must be x, y, or z; axis={axis!r}')
        if not np.isfinite(angle):
            raise ValueError(f'set {iset}: angle must be finite; angle={angle}')

        missing = np.setdiff1d(region_ids, all_regions)
        if len(missing):
            raise KeyError(f'set {iset}: region ids not found: {missing.tolist()}; '
                           f'regions={all_regions.tolist()}')
        for region_id in region_ids.tolist():
            if region_id in used_regions:
                raise ValueError(f'set {iset}: region {region_id} is also in '
                                 f'set {used_regions[region_id]}')
            used_regions[region_id] = iset
        if origin.shape != (3,):
            raise ValueError(f'set {iset}: origin must be (3,); shape={origin.shape}')
        if axes.shape != (3, 3):
            raise ValueError(f'set {iset}: axes must be (3, 3); shape={axes.shape}')
        ortho_error = np.abs(axes @ axes.T - np.eye(3)).max()
        det = np.linalg.det(axes)
        if ortho_error > ortho_tol or abs(det - 1.) > ortho_tol:
            raise ValueError(f'set {iset}: axes must be orthonormal and right-handed; '
                             f'max|axes @ axes.T - I|={ortho_error:g}; det={det:g}\n'
                             f'axes=\n{axes}')
        rotation_matrix = axes.T @ _axis_rotation_matrix(axis, angle) @ axes
        sets.append((region_ids, origin, rotation_matrix, axis, angle))
    return sets

def _axis_rotation_matrix(axis: str, angle_deg: float) -> np.ndarray:
    """right-hand rotation matrix about x, y, or z (active; v_new = R @ v)"""
    theta = np.radians(angle_deg)
    c, s = np.cos(theta), np.sin(theta)
    if axis == 'x':
        return np.array([[1., 0., 0.], [0., c, -s], [0., s, c]])
    if axis == 'y':
        return np.array([[c, 0., s], [0., 1., 0.], [-s, 0., c]])
    assert axis == 'z', axis
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
