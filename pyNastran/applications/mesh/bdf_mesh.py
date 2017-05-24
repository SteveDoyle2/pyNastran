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

        for eid, element in iteritems(self.rigid_elements):
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
