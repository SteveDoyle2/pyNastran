from pyNastran.bdf.bdf import BDF
from collections import defaultdict

class MeshTools(BDF):
    def __init__(self, log=None, debug=False):
        BDF.__init__(self, debug=debug, log=log)

    def project_node_to_plane(self, xyz, v, cid):
        """
        http://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
        [[xa-x0, ya-y0, za-z0]] = A * [[t,u,v]]
        """
        raise NotImplementedError('hasnt been validated')
        Ia = xyz
        Ib = Ia + v

        #v1 = cid.i
        #v2 = cid.j
        #v0 = cid.origin
        A = hstack([Ia-Ib, cid.i, cid.j])
        b = xa - cid.origin
        #A = [
        #    [xa-xb, x1-x0, x2-x0],
        #    [ya-yb, y1-y0, y2-y0],
        #    [za-zb, z1-z0, z2-z0]
        #]
        #b = [[xa-x0, ya-y0, za-z0]]
        # Ax=b
        x = solve(A, b)
        t = x[0]
        p_intersect = (Ia + Ib-Ia) * t
        return p_intersect

    def get_elemental_nodeids(self):
        nid_to_eids = defaultdict(list)
        eid_to_nids = {}
        for eid, element in self.elements.iteritems():
            nids = element.nodeIDs()
            eid_to_nids[eid] = nids
            for nid in nids:
                nid_to_eids[nid].append(eid)

        for eid, element in self.rigidElements.iteritems():
            nids = element.nodeIDs()
            #print nids
            eid_to_nids[eid] = nids
            for nid in nids:
                nid_to_eids[nid].append(eid)

        #print('mpc', self.mpcObject)
        for (eid, mpcs) in sorted(self.mpcs.iteritems()):
            for mpc in mpcs:
                #print str(mpc)
                nids = mpc.nodeIDs()
                eid_to_nids[eid] = nids
                for nid in nids:
                    nid_to_eids[nid].append(eid)
        return nid_to_eids, eid_to_nids

    def get_free_nodes(self, eids):
        """
        ===================================
        Example 1 - rods
        ===================================
        1        2         3
        *---11---*----22---*
        nids_to_eids = {
            1 : [11],
            2 : [11, 22],
            3 : [22],
        }
        eid_to_nids = {
            11 : [1, 2],
            22 : [2, 3],
        }
        eids = [11]
        all_used_nids = {1, 2}
        all_used_touching_eids = {11, 22}
        free_nodes = {2}

        ===================================
        Example 2 - quads
        ===================================

        1     3      5     7
        *-----*------*-----*
        |     |      |     |
        | 11  |  22  |  33 |
        |     |      |     |
        *-----*------*-----*
        2     4      6     8

        nids_to_eids = {
            1 : [11],
            2 : [11],
            3 : [11, 22],
            4 : [11, 22],
            5 : [22, 33],
            6 : [22, 33],
            7 : [33],
            8 : [33],
        }
        eid_to_nids = {
            11 : [1, 3, 4, 2],
            22 : [3, 5, 6, 4],
            33 : [5, 7, 8, 6],
        }

        eids = [11]
        used_nids = {1, 2, 3, 4}
        used_touching_eids = {11, 22}
        free_nids = {3, 4}

        eids = [22]
        used_nids = {3, 5, 6, 4}
        used_touching_eids = {11, 22, 33}
        free_nids = {1, 2, 7, 8}
        """
        nid_to_eids, eid_to_nids = self.get_elemental_nodeids()
        used_nids = set([])
        for eid in eids:
            #elem = self.elements[eid]
            e_nids = eid_to_nids[eid]
            #used_nids.add(e_nids)
            used_nids.update(e_nids)

        used_free_eids = set([])
        for nid in used_nids:
            n_eids = nid_to_eids[nid]
            used_free_eids.update(n_eids)

        used_free_nids = set([])
        for eid in used_free_eids:
            e_nids = eid_to_nids[eid]
            used_free_nids.update(e_nids)

        free_nids = used_free_nids.difference(used_nids)
        free_eids = used_free_eids.difference(eids)
        #all_eids = set(self.elements.keys())
        return used_free_eids, used_free_nids, free_eids, free_nids

    def find_volume_from_surface_elements(self):
        """
        untested
        """
        #p0 = array([0., 0., 0.])
        data = ['GRID', 1000000, 0., 0., 0.]
        p0 = GRID(data=data)
        pid = 100 # doesn't matter

        for element in self.elements:
            c = element.Centroid()
            if element.type in ['CQUAD4']:
                (n1, n2, n3) = element.nodes
                data = ['CPYRAM', p0, n1, n2, n3, n4]
                tet = CPYRAMID5(data=data)
                v = tet.Volume()
            elif element.type in ['CTRIA3']:
                (n1, n2, n3) = element.nodes
                data = ['CTETRA', pid, p0, n1, n2, n3]
                tet = CTETRA4(data=data)
                v = tet.Volume()  # this needs to be a signed volume
            else:
                raise RuntimeError('not implemented; %r' % element.type)
            V += v
        return v

    def get_cutting_plane_from_surface_elements(self):
        pass
