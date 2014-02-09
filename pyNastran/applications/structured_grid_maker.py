from numpy import arange, array, where, acos, norm

from pyNastran.bdf.bdf import BDF

def Grid(object):
    def __init__(self, bdf_name):
        self.bdf = BDF()
        self.bdf.read_bdf(bdf_name, xref=True)

        self.xyz = {}

    def build_2d_mesh(self):
        eids = arange(4719, 9878)
        starting_eid = 9876
        starting_element = self.bdf.elements[starting_eid]

        (edges_dict, edge_to_eids_dict) = self.build_edge_dicts()

        # x
        edges = edges_dict[starting_eid]

        # di = [eid1, eid2, eid3, eid4]
        d = [edge_to_eids_dict[edge] for edge in edges]

        # Li = [1, 1, 2, 2]
        L = array([len(di) for di in d])
        print "L =", L
        i = where(L==2)[0]
        assert len(i) == 2, 'eid=%s is not a corner' % starting_eid
        x = i[0]
        y = i[1]

        edge_x = edges[x]
        dx = d[x]
        next_eid_x = dx.pop(starting_eid)


        edge_y = edges[y]
        dy = d[y]
        next_eid_y = dy.pop(starting_eid)

        # x direction
        nx, eids_x = self.walk_edge(starting_eid, edge_x)
        ny, eids_y = self.walk_edge(starting_eid, edge_y)

        starting_nids = starting_element.nodeIDs()
        set_starting_nids = set(starting_nids)
        non_corner_nids = set(edge_x) + set(edge_y)
        corner_nid = set_starting_nids - non_corner_nids
        assert len(corner_nid) == 1, corner_nid

        positions = {}
        for nid, node in self.bdf.nodes:
            positions[nid] = node.Position()

        # size the array
        xyz = array( (nx, ny, 3), 'float32') )
        xyz[0, 0, :] = positions[corner_nid]

        # 7----8
        # | C  |
        # 4----3----6
        # | A  |  B |
        # 1----2----5  ---> x
        #
        # A - starting element
        # 2-3 is edge_x
        # 4-3 is edge_y
        # 1 - corner node and is on the boundary
        # 2 - part of edge_x that's not in edge_y

        # 5-6 is edge_x2
        # A*B = |A||B|cos(theta)
        # 5 - the node that has the smaller angle
        node_2x = list(set(edge_x) - set(edge_y))
        assert len(node_2x) == 1
        node_2x = node_2x[0]
        v12 = positions[node_2x] - positions[corner_nid]
        n12 = norm(v12)

        node_4y = list(set(edge_y) - set(edge_x))  # node 4
        assert len(node_4y) == 1
        node_4y = node_4y[0]
        v14 = positions[node_4y] - positions[corner_nid]
        n14 = norm(v14)

        # find the 0th row
        for icol, eid in enumerate(eids_x):
            # walk_coords(corner_nid)
            element = self.bdf.elements[eid]
            nids = element.nodeIDs()
            nids.pop(corner_nid)

            (n2, n3, n4) = nids
            edge1 = [corner_nid, n2]
            edge2 = [corner_nid, n3]
            edge3 = [corner_nid, n4]
            #edge1.sort()
            #edge2.sort()
            #edge3.sort()
            edges = [tuple(edge1), tuple(edge2), tuple(edge3)]

            A = array[ (positions[nA] - positions[nB]) for (nA, nB) in edges] )
            cosTs = dot(v12, A)/(n12*norm(A))
            thetas = acos(cosTs)

            itheta_min = where(theta==theta.min())[0][0]
            edge_nid = edge_nid_choices[itheta_min]
            xyzi = positions[edge_nid]
            xyz[0, icol, :] = xyzi

        for ix in xrange(nx):
            for iy in xrange(ny):
                xyz[0, icol, :] = xyzi






    def walk_edge(self, starting_eid, edge_x):
        nx = 1
        eids_x = [starting_eid]
        while 1:
            eid = edge_to_eids_dict[edge_x][0]
            edges = edges_dict[eid]
            iedge0 = edges.index(edge_x)

            # go to the opposite edge of the element
            iedge_new = iedge0 + 2

            # correct the edge number
            iedge = iedge_new if iedge_new < 4 else iedge_new - 4
            new_edge = edges[iedge]
            eids_next = edge_to_eids_dict[edge]
            eids_next.pop(eid)
            if len(eids_next) == 0:
                break
            eid = eids_next[0]

            nx += 1
            eids_x.append(eid)
        return nx, eids_x

    def build_edge_dicts(self, eids):
        edges_dict = {}
        edge_to_eids_dict = {}

        for eid, element in self.bdf.elements:
            if isinstance(element, CQUAD4):
                (n1, n2, n3, n4) element.nodeIDs()
                edge1 = [n1, n2]
                edge2 = [n2, n3]
                edge3 = [n3, n4]
                edge4 = [n4, n1]
                edge1.sort()
                edge2.sort()
                edge3.sort()
                edge4.sort()
                edges = [tuple(edge1), tuple(edge2), tuple(edge3), tuple(edge4)]
            else:
                raise RuntimeError('only CQUAD4s are supported')

            edges_dict[eid] = edges
            for edge in edges:
                if edge in edge_to_eids_dict:
                    edge_to_eids_dict[edge].append(eid)
                else:
                    edge_to_eids_dict[edge] = [eid]
        return edges_dict, edge_to_eids_dict