from six import iteritems
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
        for eid, element in iteritems(self.elements):
            nids = element.node_ids
            eid_to_nids[eid] = nids
            for nid in nids:
                nid_to_eids[nid].append(eid)

        for eid, element in iteritems(self.rigidElements):
            nids = element.node_ids
            #print(nids)
            eid_to_nids[eid] = nids
            for nid in nids:
                nid_to_eids[nid].append(eid)

        #print('mpc', self.mpcObject)
        for (eid, mpcs) in sorted(iteritems(self.mpcs)):
            for mpc in mpcs:
                #print(str(mpc))
                nids = mpc.node_ids
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
        raise NotImplementedError('hasnt been tested')
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

    def get_cutting_plane_from_surface_elements_x_cone(self):
        raise NotImplementedError('hasnt been implemented')

    def get_cutting_plane_from_surface_elements(self):
        raise NotImplementedError('hasnt been validated')
        cid = 1
        rid = 0
        z = 1.
        y = 2.

        ##====================================================================
        ## setup

        # transform coordinates into cutting plane coordinate system
        # where 0,0,0 is the origin
        # z is normal to the cutting plane
        # TODO: these nodes are probably wrong...
        origin = [x, 0, 0]
        z_axis = [x, 0, z]
        xz_plane = [2*x, 0,0.]
        data = ['CORD2R', cid, rid] + origin + z_axis + xz_plane
        coordA = CORD2R(data)

        # get the location of all the nodes in the model in the
        # local coordinate frame
        g0 = ['GRID', -1, coordA, 0., 0., 0.]
        grid0 = GRID(g0)

        xyz = {-1: grid0.xyz, }
        for grid in self.nodes:
            xyz[grid.nid] = grid.PositionWRT(xyz, coordA)

        ##====================================================================
        ## get the intersecting elements

        # find nodes with z < tol
        nids_tol = []
        for nid, xyz in iteritems(xyz2):
            if xyz[2] < ztol:
                nids_tol.append(nid)

        # find elements associated with nodes
        eids = self.get_x_associated_with_y(self.elements, ['nodes', 'nid'])

        ##====================================================================
        ## create areal elements
        for eid in eids:
            elem = self.elements[eid]
            if element.type == 'CQUAD4':
                #edges = element.get_edges()
                pass
            elif element.type == 'CTRIA3':
                pass
                # find subelement where positive element
                # TODO: do I need this?

                # get common node
                # TODO: do I need this?
                #common_node = union(edges[iedges[0]], edges[iedges[1]])

                #vector = (p_planes[0] + p_planes[1])
                #if common_node
                #for ip in p_planes:

            else:
                raise RuntimeError('not implemented; %r' % element.type)

            edges = element.get_edges()

            # find the two intersection points
            # TODO: what if there is only 1???
            iedges = []
            p_planes = []
            for iedge, edge in enumerate(edges):
                n1, n2 = edge
                if xyz[n2] == 0.: # on the plane
                    p_plane = xyz
                    iedges.append(iedge)
                    p_planes.append(p_plane)
                elif xyz[n1] / xyz[n2] < 0.: # different signs
                    t = xyz[n1] / xyz[n2]
                    p_plane = t * p2 + p1
                    iedges.append(iedge)
                    p_planes.append(p_plane)

            # make a CTRIA3 that connects the 2 intersecting nodes with the origin
            g1 = ['GRID', 1, 0] + p_planes[0]
            g2 = ['GRID', 2, 0] + p_planes[1]

            g1 = GRID(data=g1)
            g2 = GRID(data=g2)

            el = ['CTRIA3', 1, 1, grid0, g1, g2]
            element = CTRIA3(data=e1)

        ##====================================================================
        ## find grouped elements
        ## we can have two independent regions of elements

        # nodal equivalencing

        # find touching elements

        ##====================================================================
        ## sum area based on normals
        ## we can have two independent regions of elements
        A = 0.0
        for group in groups:
            for elem in group:
                    A += A
        return A

    def equivalence_nodes(self, xyz=None, tol=1e-6):
        """
        Collapses nodes within a tolerance
        :param xyz: dictionary of node locations by nodeID
        :param tol: the distance between nodes

        .. warning:: doesn't handle SPOINTs
        """
        if xyz is None:
            #xyz = {}
            xyz_data = zeros((len(self.nodes), 3), dtype='float64')

            i = 0
            for nid, node in sorted(iteritems(self.nodes)):
                point, M = grid.Position(xyz)
                #xyz[nid] = point
                xyz_data[i, :] = point

        from scipy.spatial import KDTree
        tree = KDTree(xyz_data, leafsize=10)

        # find the non-unique nodes
        i = 0
        all_ids = set([])
        for nid, node in sorted(iteritems(self.nodes)):
            point = xyz_data[i]

            # find the nodes at the current location
            nclose_nodes = 5
            while 1:
                dist, i_collapse = tree.query(point, k=nclose_nodes, p=2)
                if len(dist) < nclose_nodes:
                    break
                nclose_nodes += 15

            # if we have 2 or more nodes at this point
            if len(i_collapse) > 1:
                collapse_nodes = set(i_collapse)
                for ID in collapse_nodes:
                    if ID > i:
                        all_ids.add(ID)
            i += 1

        # how do I destroy nodes efficiently...
        # lets destroy nodes
        for ID in all_ids:
            self.nodes[ID] = None
