import numpy as np


def delete_bad_shells(model, max_theta=175., max_skew=70., max_aspect_ratio=100.,
                      max_taper_ratio=4.0):
    """
    Removes bad CQUAD4/CTRIA3 elements

    Parameters
    ----------
    model : BDF ()
        this should be checked

    """
    eids_to_delete = get_bad_shells(model, max_theta=max_theta,
                                    max_skew=max_skew, max_aspect_ratio=max_aspect_ratio,
                                    max_taper_ratio=max_taper_ratio)

    for eid in eids_to_delete:
        del model.elements[eid]
    model.log.info('deleted %s bad CTRIA3/CQUAD4s' % len(eids_to_delete))

    model.validate()
    return model


def get_bad_shells(model, max_theta=175., max_skew=70.,
                   max_aspect_ratio=100., max_taper_ratio=4.0):
    """
    Get the bad shell elements

    Parameters
    ----------
    model : BDF()
        the model object
    xyz_cid0 : (N, 3) float ndarray
        the xyz coordinates in cid=0
    nid_map : dict[nid] : index
        nid : int
            the node id
        index : int
            the index of the node id in xyz_cid0
    max_theta : float; default=175.
        the maximum interior angle
    max_skew : float; default=70.
        the maximum skew angle
    max_aspect_ratio : float; default=100.
        the max aspect ratio
    taper_ratio : float; default=4.0
        the taper ratio; applies to CQUAD4s only

    Returns
    -------
    eids_failed : List[int]
        element ids that fail the criteria

    shells with a edge length=0.0 are automatically added

    """
    xyz_cid0 = model.get_xyz_in_coord(cid=0, fdtype='float32')
    max_theta = np.radians(max_theta)
    max_skew = np.radians(max_skew)
    eids_failed = []

    if model.cquad4.n:
        eids_failed = get_bad_cquad4s(
            model, xyz_cid0, max_theta=max_theta, max_skew=max_skew,
            max_aspect_ratio=max_aspect_ratio, max_taper_ratio=max_taper_ratio)
    if model.ctria3.n:
        # taper ratio doesn't apply to triangles
        eids_failed = get_bad_ctria3s(
            model, xyz_cid0, max_theta=max_theta, max_skew=max_skew,
            max_aspect_ratio=max_aspect_ratio)
    return eids_failed

def get_bad_cquad4s(model, xyz_cid0, max_theta=175., max_skew=70.,
                    max_aspect_ratio=100., max_taper_ratio=4.0):
    """Get the bad elements"""
    min_theta_quad = 0.1
    min_theta_quad = np.radians(min_theta_quad)
    eids_failed = np.array([])
    elements = model.cquad4
    nelements = model.cquad4.n
    eids = elements.element_id

    i1, i2, i3, i4 = elements.get_node_indicies()

    p1 = xyz_cid0[i1, :]
    p2 = xyz_cid0[i2, :]
    p3 = xyz_cid0[i3, :]
    p4 = xyz_cid0[i4, :]
    p1 = p1.reshape((nelements, 3))
    p2 = p2.reshape((nelements, 3))
    p3 = p3.reshape((nelements, 3))
    p4 = p4.reshape((nelements, 3))

    v21 = p2 - p1
    v32 = p3 - p2
    v43 = p4 - p3
    v14 = p1 - p4

    v42 = p4 - p2
    v31 = p3 - p1
    p12 = (p1 + p2) / 2.
    p23 = (p2 + p3) / 2.
    p34 = (p3 + p4) / 2.
    p14 = (p4 + p1) / 2.
    normal = np.cross(v31, v42, axis=1)
    #    e3
    # 4-------3
    # |       |
    # |e4     |  e2
    # 1-------2
    #     e1
    e13 = p34 - p12
    e42 = p23 - p14

    #dot_e13_e42 = e13 @ e42
    norm_e13 = np.linalg.norm(e13, axis=1)
    norm_e42 = np.linalg.norm(e42, axis=1)
    assert norm_e13.size == nelements, 'norm_e13.size=%s nelements=%s' % (norm_e13.size, nelements)

    denom = norm_e13 * norm_e42
    cos_skew1 = np.inner(e13, e42) / denom
    cos_skew2 = np.inner(e13, -e42) / denom
    assert cos_skew1.size == nelements, 'cos_skew1.size=%s nelements=%s' % (cos_skew1.size, nelements)

    skew = np.pi / 2. - np.abs(np.arccos(np.clip([cos_skew1, cos_skew2], -1., 1.))).min()
    assert skew.size == nelements, 'skew.size=%s nelements=%s' % (skew.size, nelements)
    ifail = np.where(skew > max_skew)[0]
    if skew > max_skew:
        model.log.debug('fail-skew = %s' % eids[ifail])
        eids_failed = np.union1d(eids_failed, eids[ifail])
        #eids_failed.append(eid)
        #model.log.debug('eid=%s failed max_skew check; skew=%s' % (eid, np.degrees(skew)))
        #continue

    lengths = np.vstack([np.linalg.norm(v21, axis=1),
                         np.linalg.norm(v32, axis=1),
                         np.linalg.norm(v43, axis=1),
                         np.linalg.norm(v14, axis=1)])
    length_min = lengths.min()
    ifail = np.where(length_min == 0.0)[0]
    if len(ifail):
        model.log.debug('fail-length = %s' % eids[ifail])
        eids_failed = np.union1d(eids_failed, eids[ifail])
        #eids_failed.append(eid)
        #model.log.debug('eid=%s failed length_min check; length_min=%s' % (eid, length_min))
        #continue

    aspect_ratio = lengths.max() / length_min
    ifail = np.where(aspect_ratio > max_aspect_ratio)[0]
    if len(ifail):
        model.log.debug('fail-AR = %s' % eids[ifail])
        eids_failed = np.union1d(eids_failed, eids[ifail])
        #eids_failed.append(eid)
        #model.log.debug('eid=%s failed aspect_ratio check; AR=%s' % (eid, aspect_ratio))
        #continue

    area1 = 0.5 * np.linalg.norm(np.cross(-v14, v21, axis=1), axis=1) # v41 x v21
    area2 = 0.5 * np.linalg.norm(np.cross(-v21, v32, axis=1), axis=1) # v12 x v32
    area3 = 0.5 * np.linalg.norm(np.cross(v43, v32, axis=1), axis=1) # v43 x v32
    area4 = 0.5 * np.linalg.norm(np.cross(v14, -v43, axis=1), axis=1) # v14 x v34
    aavg = (area1 + area2 + area3 + area4) / 4.
    assert aavg.size == nelements, 'aavg.size=%s nelements=%s' % (aavg.size, nelements)
    taper_ratio = (abs(area1 - aavg) + abs(area2 - aavg) +
                   abs(area3 - aavg) + abs(area4 - aavg)) / aavg
    ifail = np.where(taper_ratio > max_taper_ratio)[0]
    if len(ifail):
        model.log.debug('fail-taper = %s' % eids[ifail])
        #eids_failed.append(eid)
        #model.log.debug('eid=%s failed taper_ratio check; AR=%s' % (eid, taper_ratioi))
        #continue

    # ixj = k
    # dot the local normal with the normal vector
    # then take the norm of that to determine the angle relative to the normal
    # then take the sign of that to see if we're pointing roughly towards the normal

    # np.sign(np.linalg.norm(np.dot(
    # a x b = ab sin(theta)
    # a x b / ab = sin(theta)
    # sin(theta) < 0. -> normal is flipped
    n2 = np.sign(np.inner(np.cross(v21, v32, axis=1), normal))
    n3 = np.sign(np.inner(np.cross(v32, v43, axis=1), normal))
    n4 = np.sign(np.inner(np.cross(v43, v14, axis=1), normal))
    n1 = np.sign(np.inner(np.cross(v14, v21, axis=1), normal))
    n = np.array([n1, n2, n3, n4])

    theta_additional = np.where(n < 0, 2*np.pi, 0.)

    numer = np.inner(v21, -v14)
    denom = (np.linalg.norm(v21, axis=1) * np.linalg.norm(v14, axis=1))
    assert numer.size == nelements, 'numer.size=%s nelements=%s' % (numer.size, nelements)
    assert denom.size == nelements, 'denom.size=%s nelements=%s' % (denom.size, nelements)

    cos_theta1 = np.inner(v21, -v14) / (np.linalg.norm(v21, axis=1) * np.linalg.norm(v14, axis=1))
    cos_theta2 = np.inner(v32, -v21) / (np.linalg.norm(v32, axis=1) * np.linalg.norm(v21, axis=1))
    cos_theta3 = np.inner(v43, -v32) / (np.linalg.norm(v43, axis=1) * np.linalg.norm(v32, axis=1))
    cos_theta4 = np.inner(v14, -v43) / (np.linalg.norm(v14, axis=1) * np.linalg.norm(v43, axis=1))
    assert cos_theta1.size == nelements, 'cos_theta1.size=%s nelements=%s' % (cos_theta1.size, nelements)
    interior_angle = np.arccos(np.clip(
        [cos_theta1, cos_theta2, cos_theta3, cos_theta4], -1., 1.))
    assert interior_angle.size == 4*nelements, 'interior_angle.size=%s 4*nelements=%s' % (interior_angle.size, 4*nelements)
    theta = n * interior_angle + theta_additional
    assert theta.size == 4*nelements, 'theta.size=%s 4*nelements=%s' % (theta.size, 4*nelements)

    theta_mini = theta.min(axis=0)
    theta_maxi = theta.max(axis=0)
    assert theta_mini.size == nelements, 'theta_mini.size=%s nelements=%s' % (theta_mini.size, nelements)

    ifail = np.where(theta_mini > min_theta_quad)[0]
    if len(ifail):
        model.log.debug('fail-min-theta-quad = %s' % eids[ifail])
        eids_failed = np.union1d(eids_failed, eids[ifail])
        #eids_failed.append(eid)
        #model.log.debug('eid=%s failed min_theta check; theta=%s' % (
            #eid, np.degrees(theta_mini)))
        #continue

    ifail = np.where(theta_maxi > max_theta)[0]
    if len(ifail):
        model.log.debug('fail-max-theta-quad = %s' % eids[ifail])
        eids_failed = np.union1d(eids_failed, eids[ifail])
        #eids_failed.append(eid)
        #model.log.debug('eid=%s failed max_theta check; theta=%s' % (
            #eid, np.degrees(theta_maxi)))
        #continue
    return eids_failed

def get_bad_ctria3s(model, xyz_cid0, max_theta=175., max_skew=70.,
                    max_aspect_ratio=100.):
    """Get the bad elements"""
    min_theta_tri = 0.1
    min_theta_tri = np.radians(min_theta_tri)
    eids_failed = []
    elements = model.ctria3
    nelements = model.ctria3.n
    eids = elements.element_id

    i1, i2, i3 = elements.get_node_indicies()
    p1 = xyz_cid0[i1, :]
    p2 = xyz_cid0[i2, :]
    p3 = xyz_cid0[i3, :]

    v21 = p2 - p1
    v32 = p3 - p2
    v13 = p1 - p3

    lengths = np.linalg.norm([v21, v32, v13], axis=1)

    lengths = np.vstack([
        np.linalg.norm(v21, axis=1),
        np.linalg.norm(v32, axis=1),
        np.linalg.norm(v13, axis=1),
    ])
    #assert len(lengths) == 3, lengths
    length_min = lengths.min()
    ifail = np.where(length_min == 0.0)[0]
    if len(ifail):
        model.log.debug('fail-length = %s' % eids[ifail])
        #eids_failed.append(eid)
        #model.log.debug('eid=%s failed length_min check; length_min=%s' % (eid, length_min))
        #continue

    #     3
    #    / \
    # e3/   \ e2
    #  /    /\
    # /    /  \
    # 1---/----2
    #    e1
    e1 = (p1 + p2) / 2.
    e2 = (p2 + p3) / 2.
    e3 = (p3 + p1) / 2.

    e21 = e2 - e1
    e31 = e3 - e1
    e32 = e3 - e2

    e3_p2 = e3 - p2
    e2_p1 = e2 - p1
    e1_p3 = e1 - p3
    cos_skew1 = np.inner(e2_p1, e31) / (np.linalg.norm(e2_p1, axis=1) * np.linalg.norm(e31, axis=1))
    cos_skew2 = np.inner(e2_p1, -e31) / (np.linalg.norm(e2_p1, axis=1) * np.linalg.norm(e31, axis=1))
    cos_skew3 = np.inner(e3_p2, e21) / (np.linalg.norm(e3_p2, axis=1) * np.linalg.norm(e21, axis=1))
    cos_skew4 = np.inner(e3_p2, -e21) / (np.linalg.norm(e3_p2, axis=1) * np.linalg.norm(e21, axis=1))
    cos_skew5 = np.inner(e1_p3, e32) / (np.linalg.norm(e1_p3, axis=1) * np.linalg.norm(e32, axis=1))
    cos_skew6 = np.inner(e1_p3, -e32) / (np.linalg.norm(e1_p3, axis=1) * np.linalg.norm(e32, axis=1))
    assert cos_skew1.size == nelements, 'cos_skew1.size=%s nelements=%s' % (cos_skew1.size, nelements)
    skew = np.pi / 2. - np.abs(np.arccos(
        np.clip([cos_skew1, cos_skew2, cos_skew3,
                 cos_skew4, cos_skew5, cos_skew6], -1., 1.)
    )).min()
    assert skew.size == nelements, 'skew.size=%s nelements=%s' % (skew.size, nelements)
    ifail = np.where(skew > max_skew)[0]
    if skew > max_skew:
        model.log.debug('fail-skew = %s' % eids[ifail])
        eids_failed = np.union1d(eids_failed, eids[ifail])
        #eids_failed.append(eid)
        #model.log.debug('eid=%s failed max_skew check; skew=%s' % (eid, np.degrees(skew)))
        #continue

    aspect_ratio = lengths.max() / length_min
    ifail = np.where(aspect_ratio > max_aspect_ratio)[0]
    if len(ifail):
        model.log.debug('fail-AR = %s' % eids[ifail])
        eids_failed = np.union1d(eids_failed, eids[ifail])
        #eids_failed.append(eid)
        #model.log.debug('eid=%s failed aspect_ratio check; AR=%s' % (eid, aspect_ratio))
        #continue

    cos_theta1 = np.inner(v21, -v13) / (np.linalg.norm(v21, axis=1) * np.linalg.norm(v13, axis=1))
    cos_theta2 = np.inner(v32, -v21) / (np.linalg.norm(v32, axis=1) * np.linalg.norm(v21, axis=1))
    cos_theta3 = np.inner(v13, -v32) / (np.linalg.norm(v13, axis=1) * np.linalg.norm(v32, axis=1))

    theta = np.arccos(np.clip(
        [cos_theta1, cos_theta2, cos_theta3], -1., 1.))
    theta_mini = theta.min()
    theta_maxi = theta.max()

    assert theta_mini.size == nelements, 'theta_mini.size=%s nelements=%s' % (theta_mini.size, nelements)

    ifail = np.where(theta_mini > min_theta_tri)[0]
    if len(ifail):
        model.log.debug('fail-min-theta-tri = %s' % eids[ifail])
        eids_failed = np.union1d(eids_failed, eids[ifail])
        #eids_failed.append(eid)
        #model.log.debug('eid=%s failed min_theta check; theta=%s' % (
            #eid, np.degrees(theta_mini)))
        #continue

    ifail = np.where(theta_maxi > max_theta)[0]
    if len(ifail):
        model.log.debug('fail-max-theta-tri = %s' % eids[ifail])
        eids_failed = np.union1d(eids_failed, eids[ifail])
        #eids_failed.append(eid)
        #model.log.debug('eid=%s failed max_theta check; theta=%s' % (
            #eid, np.degrees(theta_maxi)))
        #continue
    return eids_failed

#def main():  # pragma: no cover
    #delete_bad_shells()

#if __name__ == '__main__':  # pragma: no cover
    #main()
