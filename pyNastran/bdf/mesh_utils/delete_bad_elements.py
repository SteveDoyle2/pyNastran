# pylint: disable=C0103
"""
defines:
 - model = delete_bad_shells(model, max_theta=175., max_skew=70., max_aspect_ratio=100.,
                              max_taper_ratio=4.0)
 - eids_to_delete = get_bad_shells(model, xyz_cid0, nid_map, max_theta=175., max_skew=70.,
                                   max_aspect_ratio=100., max_taper_ratio=4.0)

"""
import numpy as np

SIDE_MAP = {}
SIDE_MAP['CHEXA'] = {
    1 : [4, 3, 2, 1],
    2 : [1, 2, 6, 5],
    3 : [2, 3, 7, 6],
    4 : [3, 4, 8, 7],
    5 : [4, 1, 5, 8],
    6 : [5, 6, 7, 8],
}
PIOVER2 = np.pi / 2.
PIOVER3 = np.pi / 3.


def delete_bad_shells(model,
                      min_theta=0.1, max_theta=175.,
                      max_skew=70., max_aspect_ratio=100.,
                      max_taper_ratio=4.0):
    """
    Removes bad CQUAD4/CTRIA3 elements

    Parameters
    ----------
    model : BDF ()
        this should be equivalenced
    min_theta : float; default=0.1
        the maximum interior angle (degrees)
    max_theta : float; default=175.
        the maximum interior angle (degrees)
    max_skew : float; default=70.
        the maximum skew angle (degrees)
    max_aspect_ratio : float; default=100.
        the max aspect ratio
    taper_ratio : float; default=4.0
        the taper ratio; applies to CQUAD4s only

    """
    xyz_cid0 = model.get_xyz_in_coord(cid=0, fdtype='float32')
    nid_map = get_node_map(model)
    eids_to_delete = get_bad_shells(model, xyz_cid0, nid_map,
                                    min_theta=min_theta, max_theta=max_theta,
                                    max_skew=max_skew,
                                    max_aspect_ratio=max_aspect_ratio,
                                    max_taper_ratio=max_taper_ratio)

    for eid in eids_to_delete:
        del model.elements[eid]
    model.log.info('deleted %s bad CTRIA3/CQUAD4s' % len(eids_to_delete))

    model.validate()
    return model

def get_node_map(model):
    """gets an nid->inid mapper"""
    nid_map = {}
    for i, nid in enumerate(sorted(model.nodes.keys())):
        nid_map[nid] = i
    return nid_map

def get_bad_shells(model, xyz_cid0, nid_map,
                   min_theta=0.1, max_theta=175., max_skew=70.,
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
    min_theta : float; default=0.1
        the maximum interior angle (degrees)
    max_theta : float; default=175.
        the maximum interior angle (degrees)
    max_skew : float; default=70.
        the maximum skew angle (degrees)
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
    min_theta_quad = min_theta
    min_theta_tri = min_theta
    min_theta_quad = np.radians(min_theta_quad)
    min_theta_tri = np.radians(min_theta_tri)
    max_theta = np.radians(max_theta)
    max_skew = np.radians(max_skew)
    eids_failed = []
    for eid, element in sorted(model.elements.items()):
        if element.type == 'CQUAD4':
            node_ids = element.node_ids
            #pid = element.Pid()
            #for nid in node_ids:
                #if nid is not None:
                    #nid_to_pid_map[nid].append(pid)
            #self.eid_to_nid_map[eid] = node_ids

            n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids]
            p1 = xyz_cid0[n1, :]
            p2 = xyz_cid0[n2, :]
            p3 = xyz_cid0[n3, :]
            p4 = xyz_cid0[n4, :]
            v21 = p2 - p1
            v32 = p3 - p2
            v43 = p4 - p3
            v14 = p1 - p4

            #aspect_ratio = max(p12, p23, p34, p14) / max(p12, p23, p34, p14)
            lengths = np.linalg.norm([v21, v32, v43, v14], axis=1)
            #assert len(lengths) == 3, lengths
            length_min = lengths.min()
            if length_min == 0.0:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed length_min check; length_min=%s' % (eid, length_min))
                continue

            aspect_ratio = lengths.max() / length_min
            if aspect_ratio > max_aspect_ratio:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed aspect_ratio check; AR=%.2f' % (eid, aspect_ratio))
                continue

            p12 = (p1 + p2) / 2.
            p23 = (p2 + p3) / 2.
            p34 = (p3 + p4) / 2.
            p14 = (p4 + p1) / 2.
            normal = np.cross(p3 - p1, p4 - p2)  # v42 x v31
            #    e3
            # 4-------3
            # |       |
            # |e4     |  e2
            # 1-------2
            #     e1
            e13 = p34 - p12
            e42 = p23 - p14

            cos_skew1 = (e13 @ e42) / (np.linalg.norm(e13) * np.linalg.norm(e42))
            cos_skew2 = (e13 @ -e42) / (np.linalg.norm(e13) * np.linalg.norm(e42))
            skew = np.pi / 2. - np.abs(np.arccos(np.clip([cos_skew1, cos_skew2], -1., 1.))).min()
            if skew > max_skew:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed max_skew check; skew=%.2f' % (eid, np.degrees(skew)))
                continue

            area1 = 0.5 * np.linalg.norm(np.cross(-v14, v21)) # v41 x v21
            area2 = 0.5 * np.linalg.norm(np.cross(-v21, v32)) # v12 x v32
            area3 = 0.5 * np.linalg.norm(np.cross(v43, v32)) # v43 x v32
            area4 = 0.5 * np.linalg.norm(np.cross(v14, -v43)) # v14 x v34
            aavg = (area1 + area2 + area3 + area4) / 4.
            taper_ratioi = (abs(area1 - aavg) + abs(area2 - aavg) +
                            abs(area3 - aavg) + abs(area4 - aavg)) / aavg
            if taper_ratioi > max_taper_ratio:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed taper_ratio check; taper=%.2f' % (eid, taper_ratioi))
                continue

            #if 0:
                #areai = 0.5 * np.linalg.norm(normal)
                ## still kind of in development
                ##
                ## the ratio of the ideal area to the actual area
                ## this is an hourglass check
                #areas = [
                    #np.linalg.norm(np.cross(-v14, v21)), # v41 x v21
                    #np.linalg.norm(np.cross(v32, -v21)), # v32 x v12
                    #np.linalg.norm(np.cross(v43, -v32)), # v43 x v23
                    #np.linalg.norm(np.cross(v14, v43)),  # v14 x v43
                #]
                #area_ratioi1 = areai / min(areas)
                #area_ratioi2 = max(areas) / areai
                #area_ratioi = max(area_ratioi1, area_ratioi2)

            # ixj = k
            # dot the local normal with the normal vector
            # then take the norm of that to determine the angle relative to the normal
            # then take the sign of that to see if we're pointing roughly towards the normal

            # np.sign(np.linalg.norm(np.dot(
            # a x b = ab sin(theta)
            # a x b / ab = sin(theta)
            # sin(theta) < 0. -> normal is flipped
            n2 = np.sign(np.cross(v21, v32) @ normal)
            n3 = np.sign(np.cross(v32, v43) @ normal)
            n4 = np.sign(np.cross(v43, v14) @ normal)
            n1 = np.sign(np.cross(v14, v21) @ normal)
            n = np.array([n1, n2, n3, n4])

            theta_additional = np.where(n < 0, 2*np.pi, 0.)

            cos_theta1 = (v21 @ -v14) / (np.linalg.norm(v21) * np.linalg.norm(v14))
            cos_theta2 = (v32 @ -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
            cos_theta3 = (v43 @ -v32) / (np.linalg.norm(v43) * np.linalg.norm(v32))
            cos_theta4 = (v14 @ -v43) / (np.linalg.norm(v14) * np.linalg.norm(v43))
            interior_angle = np.arccos(np.clip(
                [cos_theta1, cos_theta2, cos_theta3, cos_theta4], -1., 1.))
            theta = n * interior_angle + theta_additional

            theta_mini = theta.min()
            theta_maxi = theta.max()

            if theta_mini < min_theta_quad:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed min_theta check; theta=%.2f' % (
                    eid, np.degrees(theta_mini)))
                continue
            if theta_maxi > max_theta:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed max_theta check; theta=%.2f' % (
                    eid, np.degrees(theta_maxi)))
                continue
            #print('eid=%s theta_min=%-5.2f theta_max=%-5.2f '
                  #'skew=%-5.2f AR=%-5.2f taper_ratioi=%.2f' % (
                      #eid,
                      #np.degrees(theta_mini), np.degrees(theta_maxi),
                      #np.degrees(skew), aspect_ratio, taper_ratioi))

        elif element.type == 'CTRIA3':
            node_ids = element.node_ids
            #pid = element.Pid()
            #self.eid_to_nid_map[eid] = node_ids
            #for nid in node_ids:
                #if nid is not None:
                    #nid_to_pid_map[nid].append(pid)

            n1, n2, n3 = [nid_map[nid] for nid in node_ids]
            p1 = xyz_cid0[n1, :]
            p2 = xyz_cid0[n2, :]
            p3 = xyz_cid0[n3, :]

            v21 = p2 - p1
            v32 = p3 - p2
            v13 = p1 - p3
            lengths = np.linalg.norm([v21, v32, v13], axis=1)
            length_min = lengths.min()
            if length_min == 0.0:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed length_min check; length_min=%s' % (
                    eid, length_min))
                continue

            #assert len(lengths) == 3, lengths
            aspect_ratio = lengths.max() / length_min
            if aspect_ratio > max_aspect_ratio:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed aspect_ratio check; AR=%s' % (eid, aspect_ratio))
                continue

            cos_theta1 = (v21 @ -v13) / (np.linalg.norm(v21) * np.linalg.norm(v13))
            cos_theta2 = (v32 @ -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
            cos_theta3 = (v13 @ -v32) / (np.linalg.norm(v13) * np.linalg.norm(v32))

            theta = np.arccos(np.clip(
                [cos_theta1, cos_theta2, cos_theta3], -1., 1.))
            theta_mini = theta.min()
            theta_maxi = theta.max()

            if theta_mini < min_theta_tri:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed min_theta check; theta=%s' % (
                    eid, np.degrees(theta_mini)))
                continue
            if theta_maxi > max_theta:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed max_theta check; theta=%s' % (
                    eid, np.degrees(theta_maxi)))
                continue

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
            cos_skew1 = (e2_p1 @ e31) / (np.linalg.norm(e2_p1) * np.linalg.norm(e31))
            cos_skew2 = (e2_p1 @ -e31) / (np.linalg.norm(e2_p1) * np.linalg.norm(e31))
            cos_skew3 = (e3_p2 @ e21) / (np.linalg.norm(e3_p2) * np.linalg.norm(e21))
            cos_skew4 = (e3_p2 @ -e21) / (np.linalg.norm(e3_p2) * np.linalg.norm(e21))
            cos_skew5 = (e1_p3 @ e32) / (np.linalg.norm(e1_p3) * np.linalg.norm(e32))
            cos_skew6 = (e1_p3 @ -e32) / (np.linalg.norm(e1_p3) * np.linalg.norm(e32))
            skew = np.pi / 2. - np.abs(np.arccos(
                np.clip([cos_skew1, cos_skew2, cos_skew3,
                         cos_skew4, cos_skew5, cos_skew6], -1., 1.)
            )).min()
            if skew > max_skew:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed max_skew check; skew=%s' % (eid, np.degrees(skew)))
                continue

            #print('eid=%s theta_min=%-5.2f theta_max=%-5.2f skew=%-5.2f AR=%-5.2f' % (
                #eid,
                #np.degrees(theta_mini), np.degrees(theta_maxi),
                #np.degrees(skew), aspect_ratio))
    return eids_failed

def element_quality(model, nids=None, xyz_cid0=None, nid_map=None):
    """
    Gets various measures of element quality

    Parameters
    ----------
    model : BDF()
        a cross-referenced model
    nids : (nnodes, ) int ndarray; default=None
        the nodes of the model in sorted order
        includes GRID, SPOINT, & EPOINTs
    xyz_cid0 : (nnodes, 3) float ndarray; default=None
        the associated global xyz locations
    nid_map : Dict[nid]->index; default=None
        a mapper dictionary

    Returns
    -------
    quality : Dict[name] : (nelements, ) float ndarray
        Various quality metrics
        names : min_interior_angle, max_interior_angle, dideal_theta,
                max_skew_angle, max_aspect_ratio,
                area_ratio, taper_ratio, min_edge_length
        values : The result is ``np.nan`` if element type does not define
                 the parameter.  For example, CELAS1 doesn't have an
                 aspect ratio.

    Notes
    -----
     - pulled from nastran_io.py

    """
    if nids is None or xyz_cid0 is None:
        out = model.get_displacement_index_xyz_cp_cd(
            fdtype='float64', idtype='int32', sort_ids=True)
        unused_icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
        nids = nid_cp_cd[:, 0]
        xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
            xyz_cp, nids, icp_transform, cid=0,
            in_place=False)

    if nid_map is None:
        nid_map = {}
        for i, nid in enumerate(nids):
            nid_map[nid] = i

    all_nids = nids
    del nids

    # these normals point inwards
    #      4
    #    / | \
    #   /  |  \
    #  3-------2
    #   \  |   /
    #    \ | /
    #      1
    _ctetra_faces = (
        (0, 1, 2), # (1, 2, 3),
        (0, 3, 1), # (1, 4, 2),
        (0, 3, 2), # (1, 3, 4),
        (1, 3, 2), # (2, 4, 3),
    )

    # these normals point inwards
    #
    #
    #
    #
    #        /4-----3
    #       /       /
    #      /  5    /
    #    /    \   /
    #   /      \ /
    # 1---------2
    _cpyram_faces = (
        (0, 1, 2, 3), # (1, 2, 3, 4),
        (1, 4, 2), # (2, 5, 3),
        (2, 4, 3), # (3, 5, 4),
        (0, 3, 4), # (1, 4, 5),
        (0, 4, 1), # (1, 5, 2),
    )

    # these normals point inwards
    #       /6
    #     /  | \
    #   /    |   \
    # 3\     |     \
    # |  \   /4-----5
    # |    \/       /
    # |   /  \     /
    # |  /    \   /
    # | /      \ /
    # 1---------2
    _cpenta_faces = (
        (0, 2, 1), # (1, 3, 2),
        (3, 4, 5), # (4, 5, 6),

        (0, 1, 4, 3), # (1, 2, 5, 4), # bottom
        (1, 2, 5, 4), # (2, 3, 6, 5), # right
        (0, 3, 5, 2), # (1, 4, 6, 3), # left
    )

    # these normals point inwards
    #      8----7
    #     /|   /|
    #    / |  / |
    #   /  5-/--6
    # 4-----3   /
    # |  /  |  /
    # | /   | /
    # 1-----2
    _chexa_faces = (
        (4, 5, 6, 7), # (5, 6, 7, 8),
        (0, 3, 2, 1), # (1, 4, 3, 2),
        (1, 2, 6, 5), # (2, 3, 7, 6),
        (2, 3, 7, 6), # (3, 4, 8, 7),
        (0, 4, 7, 3), # (1, 5, 8, 4),
        (0, 6, 5, 4), # (1, 7, 6, 5),
    )

    # quality
    nelements = len(model.elements)
    min_interior_angle = np.zeros(nelements, 'float32')
    max_interior_angle = np.zeros(nelements, 'float32')
    dideal_theta = np.zeros(nelements, 'float32')
    max_skew_angle = np.zeros(nelements, 'float32')
    max_warp_angle = np.zeros(nelements, 'float32')
    max_aspect_ratio = np.zeros(nelements, 'float32')
    #area = np.zeros(nelements, 'float32')
    area_ratio = np.zeros(nelements, 'float32')
    taper_ratio = np.zeros(nelements, 'float32')
    min_edge_length = np.zeros(nelements, 'float32')
    #normals = np.full((nelements, 3), np.nan, 'float32')


    #nids_list = []
    ieid = 0
    for unused_eid, elem in sorted(model.elements.items()):
        if ieid % 5000 == 0 and ieid > 0:
            print('  map_elements = %i' % ieid)
        etype = elem.type
        nids = None
        inids = None

        dideal_thetai = np.nan
        min_thetai = np.nan
        max_thetai = np.nan
        #max_thetai = np.nan
        max_skew = np.nan
        max_warp = np.nan
        aspect_ratio = np.nan
        #areai = np.nan
        area_ratioi = np.nan
        taper_ratioi = np.nan
        min_edge_lengthi = np.nan
        #normali = np.nan
        if etype in ['CTRIA3', 'CTRIAR', 'CTRAX3', 'CPLSTN3']:
            nids = elem.nodes
            inids = np.searchsorted(all_nids, nids)
            p1, p2, p3 = xyz_cid0[inids, :]
            out = tri_quality(p1, p2, p3)
            (areai, max_skew, aspect_ratio,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
            #normali = np.cross(p1 - p2, p1 - p3)

        elif etype in ['CQUAD4', 'CQUADR', 'CPLSTN4', 'CQUADX4']:
            nids = elem.nodes
            inids = np.searchsorted(all_nids, nids)
            p1, p2, p3, p4 = xyz_cid0[inids, :]
            out = quad_quality(elem, p1, p2, p3, p4)
            (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warp) = out

        elif etype == 'CTRIA6':
            nids = elem.nodes
            if None in nids:
                inids = np.searchsorted(all_nids, nids[:3])
                nids = nids[:3]
                p1, p2, p3 = xyz_cid0[inids, :]
            else:
                inids = np.searchsorted(all_nids, nids)
                p1, p2, p3, p4, unused_p5, unused_p6 = xyz_cid0[inids, :]
            out = tri_quality(p1, p2, p3)
            (areai, max_skew, aspect_ratio,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out

        elif etype == 'CQUAD8':
            nids = elem.nodes
            if None in nids:
                inids = np.searchsorted(all_nids, nids[:4])
                nids = nids[:4]
                p1, p2, p3, p4 = xyz_cid0[inids, :]
            else:
                inids = np.searchsorted(all_nids, nids)
                p1, p2, p3, p4 = xyz_cid0[inids[:4], :]
            out = quad_quality(elem, p1, p2, p3, p4)
            (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warp) = out
            #normali = np.cross(p1 - p3, p2 - p4)

        elif etype == 'CSHEAR':
            nids = elem.nodes
            inids = np.searchsorted(all_nids, nids)
            p1, p2, p3, p4 = xyz_cid0[inids, :]
            out = quad_quality(elem, p1, p2, p3, p4)
            (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warp) = out

        elif etype == 'CTETRA':
            nids = elem.nodes
            if None in nids:
                nids = nids[:4]
            inids = np.searchsorted(all_nids, nids)
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _ctetra_faces, nids, nid_map, xyz_cid0)

        elif etype == 'CHEXA':
            nids = elem.nodes
            if None in nids:
                nids = nids[:8]
            inids = np.searchsorted(all_nids, nids)
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _chexa_faces, nids, nid_map, xyz_cid0)

        elif etype == 'CPENTA':
            nids = elem.nodes
            if None in nids:
                nids = nids[:6]
            inids = np.searchsorted(all_nids, nids)
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _cpenta_faces, nids, nid_map, xyz_cid0)

        elif etype == 'CPYRAM':
            # TODO: assuming 5
            nids = elem.nodes
            if None in nids:
                nids = nids[:5]
            inids = np.searchsorted(all_nids, nids)
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _cpyram_faces, nids, nid_map, xyz_cid0)
        elif etype in ['CELAS2', 'CELAS4', 'CDAMP4']:
            # these can have empty nodes and have no property
            # CELAS1: 1/2 GRID/SPOINT and pid
            # CELAS2: 1/2 GRID/SPOINT, k, ge, and s
            # CELAS3: 1/2 SPOINT and pid
            # CELAS4: 1/2 SPOINT and k
            continue
            #nids = elem.nodes
            #assert nids[0] != nids[1]
            #if None in nids:
                #assert nids[0] is not None, nids
                #assert nids[1] is None, nids
                #nids = [nids[0]]
            #else:
                #nids = elem.nodes
                #assert nids[0] != nids[1]
            #inids = np.searchsorted(all_nids, nids)
        elif etype in ['CBUSH', 'CBUSH1D', 'CBUSH2D',
                       'CELAS1', 'CELAS3',
                       'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP5',
                       'CFAST', 'CGAP', 'CVISC']:
            #nids = elem.nodes
            #assert nids[0] != nids[1]
            #assert None not in nids, 'nids=%s\n%s' % (nids, elem)
            #inids = np.searchsorted(all_nids, nids)
            continue
        elif etype in ['CBAR', 'CBEAM']:
            nids = elem.nodes
            inids = np.searchsorted(all_nids, nids)
            p1, p2 = xyz_cid0[inids, :]
            min_edge_lengthi = np.linalg.norm(p2 - p1)
        elif etype in ['CROD', 'CTUBE']:
            nids = elem.nodes
            inids = np.searchsorted(all_nids, nids)
            p1, p2 = xyz_cid0[inids, :]
            min_edge_lengthi = np.linalg.norm(p2 - p1)
            #nnodes = 2
            #dim = 1
        elif etype == 'CONROD':
            nids = elem.nodes
            inids = np.searchsorted(all_nids, nids)
            p1, p2 = xyz_cid0[inids, :]
            min_edge_lengthi = np.linalg.norm(p2 - p1)
        #------------------------------
        # rare
        #elif etype == 'CIHEX1':
            #nids = elem.nodes
            #pid = elem.pid
            #cell_type = cell_type_hexa8
            #inids = np.searchsorted(all_nids, nids)
            #min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                #_chexa_faces, nids, nid_map, xyz_cid0)
            #nnodes = 8
            #dim = 3
        elif etype == 'CHBDYE':
            continue
            #self.eid_map[eid] = ieid
            #eid_solid = elem.eid2
            #side = elem.side
            #element_solid = model.elements[eid_solid]

            #mapped_inids = SIDE_MAP[element_solid.type][side]
            #side_inids = [nid - 1 for nid in mapped_inids]
            #nodes = element_solid.node_ids

            ##nnodes = len(side_inids)
            #nids = [nodes[inid] for inid in side_inids]
            #inids = np.searchsorted(all_nids, nids)

            #if len(side_inids) == 4:
                #pass
            #else:
                #msg = 'element_solid:\n%s' % (str(element_solid))
                #msg += 'mapped_inids = %s\n' % mapped_inids
                #msg += 'side_inids = %s\n' % side_inids
                #msg += 'nodes = %s\n' % nodes
                ##msg += 'side_nodes = %s\n' % side_nodes
                #raise NotImplementedError(msg)
        else:
            #raise NotImplementedError(elem)
            nelements -= 1
            continue
        #nids_list.append(nnodes)
        #nids_list.extend(inids)
        #normals[ieid] = normali
        #eids_array[ieid] = eid
        #pids_array[ieid] = pid
        #dim_array[ieid] = dim
        #cell_types_array[ieid] = cell_type
        #cell_offsets_array[ieid] = cell_offset  # I assume the problem is here
        #cell_offset += nnodes + 1
        #eid_map[eid] = ieid

        min_interior_angle[ieid] = min_thetai
        max_interior_angle[ieid] = max_thetai
        dideal_theta[ieid] = dideal_thetai
        max_skew_angle[ieid] = max_skew
        max_warp_angle[ieid] = max_warp
        max_aspect_ratio[ieid] = aspect_ratio
        #area[ieid] = areai
        area_ratio[ieid] = area_ratioi
        taper_ratio[ieid] = taper_ratioi
        min_edge_length[ieid] = min_edge_lengthi
        ieid += 1
    quality = {
        'min_interior_angle' : min_interior_angle,
        'max_interior_angle' : max_interior_angle,
        'dideal_theta' : dideal_theta,
        'max_skew_angle' : max_skew_angle,
        'max_warp_angle' : max_warp_angle,
        'max_aspect_ratio' : max_aspect_ratio,
        #'area' : area,
        'area_ratio' : area_ratio,
        'taper_ratio' : taper_ratio,
        'min_edge_length' : min_edge_length,
    }
    return quality

def tri_quality(p1, p2, p3):
    """gets the quality metrics for a tri"""
    e1 = (p1 + p2) / 2.
    e2 = (p2 + p3) / 2.
    e3 = (p3 + p1) / 2.

    #    3
    #    / \
    # e3/   \ e2
    #  /    /\
    # /    /  \
    # 1---/----2
    #    e1
    e21 = e2 - e1
    e31 = e3 - e1
    e32 = e3 - e2

    e3_p2 = e3 - p2
    e2_p1 = e2 - p1
    e1_p3 = e1 - p3

    v21 = p2 - p1
    v32 = p3 - p2
    v13 = p1 - p3
    length21 = np.linalg.norm(v21)
    length32 = np.linalg.norm(v32)
    length13 = np.linalg.norm(v13)
    min_edge_length = min(length21, length32, length13)
    area = 0.5 * np.linalg.norm(np.cross(v21, v13))

    ne31 = np.linalg.norm(e31)
    ne21 = np.linalg.norm(e21)
    ne32 = np.linalg.norm(e32)
    ne2_p1 = np.linalg.norm(e2_p1)
    ne3_p2 = np.linalg.norm(e3_p2)
    ne1_p3 = np.linalg.norm(e1_p3)
    cos_skew1 = (e2_p1 @ e31) / (ne2_p1 * ne31)
    cos_skew2 = (e2_p1 @ -e31) / (ne2_p1 * ne31)
    cos_skew3 = (e3_p2 @ e21) / (ne3_p2 * ne21)
    cos_skew4 = (e3_p2 @ -e21) / (ne3_p2 * ne21)
    cos_skew5 = (e1_p3 @ e32) / (ne1_p3 * ne32)
    cos_skew6 = (e1_p3 @ -e32) / (ne1_p3 * ne32)
    max_skew = np.pi / 2. - np.abs(np.arccos(np.clip([
        cos_skew1, cos_skew2, cos_skew3,
        cos_skew4, cos_skew5, cos_skew6], -1., 1.))).min()
    lengths = np.linalg.norm([v21, v32, v13], axis=1)
    #assert len(lengths) == 3, lengths
    length_min = lengths.min()
    if length_min == 0.0:
        aspect_ratio = np.nan
        # assume length_min = length21 = nan, so:
        #   cos_theta1 = nan
        #   thetas = [nan, b, c]
        #   min_theta = max_theta = dideal_theta = nan
        min_theta = np.nan
        max_theta = np.nan
        dideal_theta = np.nan
    else:
        aspect_ratio = lengths.max() / length_min

        cos_theta1 = (v21 @ -v13) / (length21 * length13)
        cos_theta2 = (v32 @ -v21) / (length32 * length21)
        cos_theta3 = (v13 @ -v32) / (length13 * length32)
        thetas = np.arccos(np.clip([cos_theta1, cos_theta2, cos_theta3], -1., 1.))
        min_theta = thetas.min()
        max_theta = thetas.max()
        dideal_theta = max(max_theta - PIOVER3, PIOVER3 - min_theta)

    #theta_deg = np.degrees(np.arccos(max_cos_theta))
    #if theta_deg < 60.:
        #print('p1=%s' % xyz_cid0[p1, :])
        #print('p2=%s' % xyz_cid0[p2, :])
        #print('p3=%s' % xyz_cid0[p3, :])
        #print('theta1=%s' % np.degrees(np.arccos(cos_theta1)))
        #print('theta2=%s' % np.degrees(np.arccos(cos_theta2)))
        #print('theta3=%s' % np.degrees(np.arccos(cos_theta3)))
        #print('max_theta=%s' % theta_deg)
        #asdf
    return area, max_skew, aspect_ratio, min_theta, max_theta, dideal_theta, min_edge_length


def quad_quality(element, p1, p2, p3, p4):
    """gets the quality metrics for a quad"""
    v21 = p2 - p1
    v32 = p3 - p2
    v43 = p4 - p3
    v14 = p1 - p4
    length21 = np.linalg.norm(v21)
    length32 = np.linalg.norm(v32)
    length43 = np.linalg.norm(v43)
    length14 = np.linalg.norm(v14)
    min_edge_length = min(length21, length32, length43, length14)

    v42 = p4 - p2
    v31 = p3 - p1
    p12 = (p1 + p2) / 2.
    p23 = (p2 + p3) / 2.
    p34 = (p3 + p4) / 2.
    p14 = (p4 + p1) / 2.
    v31 = p3 - p1
    v42 = p4 - p2
    normal = np.cross(v31, v42)
    area = 0.5 * np.linalg.norm(normal)

    # still kind of in development
    #
    # the ratio of the ideal area to the actual area
    # this is an hourglass check
    areas = [
        np.linalg.norm(np.cross(-v14, v21)), # v41 x v21
        np.linalg.norm(np.cross(v32, -v21)), # v32 x v12
        np.linalg.norm(np.cross(v43, -v32)), # v43 x v23
        np.linalg.norm(np.cross(v14, v43)),  # v14 x v43
    ]
    #
    # for:
    #   area=1; area1=0.5 -> area_ratioi1=2.0; area_ratio=2.0
    #   area=1; area1=2.0 -> area_ratioi2=2.0; area_ratio=2.0
    min_area = min(areas)
    if min_area == 0.:
        print('nan min area; area=%g areas=%s:\n%s' % (area, areas, element))
        #nodes = element.nodes
        #print('  n_%i = %s' % (nodes[0], p1))
        #print('  n_%i = %s' % (nodes[1], p2))
        #print('  n_%i = %s' % (nodes[2], p3))
        #print('  n_%i = %s' % (nodes[3], p4))
        area_ratio = np.nan
        #raise RuntimeError('bad quad...')
    else:
        area_ratioi1 = area / min_area
        area_ratioi2 = max(areas) / area
        area_ratio = max(area_ratioi1, area_ratioi2)

    area1 = 0.5 * areas[0] # v41 x v21
    area2 = 0.5 * np.linalg.norm(np.cross(-v21, v32)) # v12 x v32
    area3 = 0.5 * np.linalg.norm(np.cross(v43, v32)) # v43 x v32
    area4 = 0.5 * np.linalg.norm(np.cross(v14, -v43)) # v14 x v34
    aavg = (area1 + area2 + area3 + area4) / 4.
    taper_ratio = (abs(area1 - aavg) + abs(area2 - aavg) +
                   abs(area3 - aavg) + abs(area4 - aavg)) / aavg

    #    e3
    # 4-------3
    # |       |
    # |e4     |  e2
    # 1-------2
    #     e1
    e13 = p34 - p12
    e42 = p23 - p14
    ne42 = np.linalg.norm(e42)
    ne13 = np.linalg.norm(e13)
    cos_skew1 = (e13 @ e42) / (ne13 * ne42)
    cos_skew2 = (e13 @ -e42) / (ne13 * ne42)
    max_skew = np.pi / 2. - np.abs(np.arccos(
        np.clip([cos_skew1, cos_skew2], -1., 1.))).min()
    #aspect_ratio = max(p12, p23, p34, p14) / max(p12, p23, p34, p14)
    lengths = np.linalg.norm([v21, v32, v43, v14], axis=1)
    #assert len(lengths) == 3, lengths
    aspect_ratio = lengths.max() / lengths.min()

    cos_theta1 = (v21 @ -v14) / (length21 * length14)
    cos_theta2 = (v32 @ -v21) / (length32 * length21)
    cos_theta3 = (v43 @ -v32) / (length43 * length32)
    cos_theta4 = (v14 @ -v43) / (length14 * length43)
    #max_thetai = np.arccos([cos_theta1, cos_theta2, cos_theta3, cos_theta4]).max()

    # dot the local normal with the normal vector
    # then take the norm of that to determine the angle relative to the normal
    # then take the sign of that to see if we're pointing roughly towards the normal
    #
    # np.sign(np.linalg.norm(np.dot(
    # a x b = ab sin(theta)
    # a x b / ab = sin(theta)
    # sin(theta) < 0. -> normal is flipped
    normal2 = np.sign(np.cross(v21, v32) @ normal)
    normal3 = np.sign(np.cross(v32, v43) @ normal)
    normal4 = np.sign(np.cross(v43, v14) @ normal)
    normal1 = np.sign(np.cross(v14, v21) @ normal)
    n = np.array([normal1, normal2, normal3, normal4])
    theta_additional = np.where(n < 0, 2*np.pi, 0.)

    theta = n * np.arccos(np.clip(
        [cos_theta1, cos_theta2, cos_theta3, cos_theta4], -1., 1.)) + theta_additional
    min_theta = theta.min()
    max_theta = theta.max()
    dideal_theta = max(max_theta - PIOVER2, PIOVER2 - min_theta)
    #print('theta_max = ', theta_max)


    # warp angle
    # split the quad and find the normals of each triangl
    # find the angle between the two triangles
    #
    # 4---3
    # | / |
    # |/  |
    # 1---2
    #
    v41 = -v14
    n123 = np.cross(v21, v31)
    n134 = np.cross(v31, v41)
    #v1 o v2 = v1 * v2 cos(theta)
    cos_warp1 = (n123 @ n134) / (np.linalg.norm(n123) * np.linalg.norm(n134))

    # split the quad in the order direction and take the maximum of the two splits
    # 4---3
    # | \ |
    # |  \|
    # 1---2
    n124 = np.cross(v21, v41)
    n234 = np.cross(v32, v42)
    cos_warp2 = (n124 @ n234) / (np.linalg.norm(n124) * np.linalg.norm(n234))

    max_warp = np.abs(np.arccos(
        np.clip([cos_warp1, cos_warp2], -1., 1.))).max()

    out = (area, taper_ratio, area_ratio, max_skew, aspect_ratio,
           min_theta, max_theta, dideal_theta, min_edge_length, max_warp)
    return out

def get_min_max_theta(faces, all_node_ids, nid_map, xyz_cid0):
    """get the min/max thetas for CTETRA, CPENTA, CHEXA, CPYRAM"""
    cos_thetas = []
    ideal_theta = []
    #print('faces =', faces)
    #assert len(faces) > 0, 'faces=%s nids=%s' % (faces, all_node_ids)
    for face in faces:
        if len(face) == 3:
            node_ids = all_node_ids[face[0]], all_node_ids[face[1]], all_node_ids[face[2]]
            n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
            v21 = xyz_cid0[n2, :] - xyz_cid0[n1, :]
            v32 = xyz_cid0[n3, :] - xyz_cid0[n2, :]
            v13 = xyz_cid0[n1, :] - xyz_cid0[n3, :]
            length21 = np.linalg.norm(v21)
            length32 = np.linalg.norm(v32)
            length13 = np.linalg.norm(v13)
            min_edge_length = min(length21, length32, length13)

            cos_theta1 = (v21 @ -v13) / (length21 * length13)
            cos_theta2 = (v32 @ -v21) / (length32 * length21)
            cos_theta3 = (v13 @ -v32) / (length13 * length32)
            cos_thetas.extend([cos_theta1, cos_theta2, cos_theta3])
            ideal_theta.extend([PIOVER3, PIOVER3, PIOVER3])

        elif len(face) == 4:
            try:
                node_ids = (all_node_ids[face[0]], all_node_ids[face[1]],
                            all_node_ids[face[2]], all_node_ids[face[3]])
            except KeyError:
                print(face)
                print(node_ids)
                raise

            n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
            v21 = xyz_cid0[n2, :] - xyz_cid0[n1, :]
            v32 = xyz_cid0[n3, :] - xyz_cid0[n2, :]
            v43 = xyz_cid0[n4, :] - xyz_cid0[n3, :]
            v14 = xyz_cid0[n1, :] - xyz_cid0[n4, :]
            length21 = np.linalg.norm(v21)
            length32 = np.linalg.norm(v32)
            length43 = np.linalg.norm(v43)
            length14 = np.linalg.norm(v14)
            min_edge_length = min(length21, length32, length43, length14)
            cos_theta1 = (v21 @ -v14) / (length21 * length14)
            cos_theta2 = (v32 @ -v21) / (length32 * length21)
            cos_theta3 = (v43 @ -v32) / (length43 * length32)
            cos_theta4 = (v14 @ -v43) / (length14 * length43)
            cos_thetas.extend([cos_theta1, cos_theta2, cos_theta3, cos_theta4])
            ideal_theta.extend([PIOVER2, PIOVER2, PIOVER2, PIOVER2])
        else:
            raise NotImplementedError(face)
    thetas = np.arccos(cos_thetas)
    ideal_theta = np.array(ideal_theta)
    ideal_thetai = max((thetas - ideal_theta).max(), (ideal_theta - thetas).min())

    min_thetai = thetas.min()
    max_thetai = thetas.max()
    return min_thetai, max_thetai, ideal_thetai, min_edge_length
