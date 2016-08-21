from __future__ import print_function
import numpy as np
from six import iteritems


def delete_bad_shells(model, max_theta=175., max_skew=70., max_aspect_ratio=100.,
                      max_taper_ratio=4.0):
    """
    Removes bad CQUAD4/CTRIA3 elements

    Parameters
    ----------
    model : BDF ()
        this should be equivalenced
    """
    xyz_cid0 = model.get_xyz_in_coord(cid=0, dtype='float32')
    nid_map = {}
    for i, (nid, node) in enumerate(sorted(iteritems(model.nodes))):
        #xyz = node.get_position()
        #xyz_cid0[i, :] = xyz
        nid_map[nid] = i
    eids_to_delete = get_bad_shells(model, xyz_cid0, nid_map, max_theta=max_theta,
                                    max_skew=max_skew, max_aspect_ratio=max_aspect_ratio,
                                    max_taper_ratio=max_taper_ratio)

    for eid in eids_to_delete:
        del model.elements[eid]
    model.log.info('deleted %s bad CTRIA3/CQUAD4s' % len(eids_to_delete))

    model.validate()
    return model


def get_bad_shells(model, xyz_cid0, nid_map, max_theta=175., max_skew=70.,
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
    max_theta = np.radians(max_theta)
    max_skew = np.radians(max_skew)
    eids_failed = []
    for eid, element in sorted(iteritems(model.elements)):
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

            v42 = p4 - p2
            v31 = p3 - p1
            p12 = (p1 + p2) / 2.
            p23 = (p2 + p3) / 2.
            p34 = (p3 + p4) / 2.
            p14 = (p4 + p1) / 2.
            normal = np.cross(v31, v42)
            #    e3
            # 4-------3
            # |       |
            # |e4     |  e2
            # 1-------2
            #     e1
            e13 = p34 - p12
            e42 = p23 - p14

            cos_skew1 = np.dot(e13, e42) / (np.linalg.norm(e13) * np.linalg.norm(e42))
            cos_skew2 = np.dot(e13, -e42) / (np.linalg.norm(e13) * np.linalg.norm(e42))
            skew = np.pi / 2. - np.abs(np.arccos(np.clip([cos_skew1, cos_skew2], -1., 1.))).min()
            if skew > max_skew:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed max_skew check; skew=%s' % (eid, np.degrees(skew)))
                continue

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
                model.log.debug('eid=%s failed aspect_ratio check; AR=%s' % (eid, aspect_ratio))
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
                model.log.debug('eid=%s failed taper_ratio check; AR=%s' % (eid, taper_ratioi))
                continue

            # ixj = k
            # dot the local normal with the normal vector
            # then take the norm of that to determine the angle relative to the normal
            # then take the sign of that to see if we're pointing roughly towards the normal

            # np.sign(np.linalg.norm(np.dot(
            # a x b = ab sin(theta)
            # a x b / ab = sin(theta)
            # sin(theta) < 0. -> normal is flipped
            n2 = np.sign(np.dot(np.cross(v21, v32), normal))
            n3 = np.sign(np.dot(np.cross(v32, v43), normal))
            n4 = np.sign(np.dot(np.cross(v43, v14), normal))
            n1 = np.sign(np.dot(np.cross(v14, v21), normal))
            n = np.array([n1, n2, n3, n4])

            theta_additional = np.where(n < 0, 2*np.pi, 0.)

            cos_theta1 = np.dot(v21, -v14) / (np.linalg.norm(v21) * np.linalg.norm(v14))
            cos_theta2 = np.dot(v32, -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
            cos_theta3 = np.dot(v43, -v32) / (np.linalg.norm(v43) * np.linalg.norm(v32))
            cos_theta4 = np.dot(v14, -v43) / (np.linalg.norm(v14) * np.linalg.norm(v43))
            interior_angle = np.arccos(np.clip(
                [cos_theta1, cos_theta2, cos_theta3, cos_theta4], -1., 1.))
            theta = n * interior_angle + theta_additional

            theta_maxi = theta.max()
            if theta_maxi > max_theta:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed max_theta check; theta=%s' % (
                    eid, np.degrees(theta_maxi)))
                continue

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
            e1 = (p1 + p2) / 2.
            e2 = (p2 + p3) / 2.
            e3 = (p3 + p1) / 2.

            #     3
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
            cos_skew1 = np.dot(e2_p1, e31) / (np.linalg.norm(e2_p1) * np.linalg.norm(e31))
            cos_skew2 = np.dot(e2_p1, -e31) / (np.linalg.norm(e2_p1) * np.linalg.norm(e31))
            cos_skew3 = np.dot(e3_p2, e21) / (np.linalg.norm(e3_p2) * np.linalg.norm(e21))
            cos_skew4 = np.dot(e3_p2, -e21) / (np.linalg.norm(e3_p2) * np.linalg.norm(e21))
            cos_skew5 = np.dot(e1_p3, e32) / (np.linalg.norm(e1_p3) * np.linalg.norm(e32))
            cos_skew6 = np.dot(e1_p3, -e32) / (np.linalg.norm(e1_p3) * np.linalg.norm(e32))
            skew = np.pi / 2. - np.abs(np.arccos(
                np.clip([cos_skew1, cos_skew2, cos_skew3,
                         cos_skew4, cos_skew5, cos_skew6], -1., 1.)
            )).min()
            if skew > max_skew:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed max_skew check; skew=%s' % (eid, np.degrees(skew)))
                continue

            lengths = np.linalg.norm([v21, v32, v13], axis=1)
            length_min = lengths.min()
            if length_min == 0.0:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed length_min check; length_min=%s' % (eid, length_min))
                continue

            #assert len(lengths) == 3, lengths
            aspect_ratio = lengths.max() / length_min
            if aspect_ratio > max_aspect_ratio:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed aspect_ratio check; AR=%s' % (eid, aspect_ratio))
                continue

            cos_theta1 = np.dot(v21, -v13) / (np.linalg.norm(v21) * np.linalg.norm(v13))
            cos_theta2 = np.dot(v32, -v21) / (np.linalg.norm(v32) * np.linalg.norm(v21))
            cos_theta3 = np.dot(v13, -v32) / (np.linalg.norm(v13) * np.linalg.norm(v32))

            interior_angle = np.arccos(np.clip(
                [cos_theta1, cos_theta2, cos_theta3], -1., 1.))
            theta = n * interior_angle + theta_additional
            theta_mini = theta.min()
            theta_maxi = theta.max()

            if theta_maxi > max_theta:
                eids_failed.append(eid)
                model.log.debug('eid=%s failed max_theta check; theta=%s' % (
                    eid, np.degrees(theta_maxi)))
                continue
    return eids_failed

#def main():
    #delete_bad_shells()

#if __name__ == '__main__':
    #main()
