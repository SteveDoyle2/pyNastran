from __future__ import print_function
from numpy import zeros, cross, dot, allclose, sign
from numpy.linalg import norm

from pyNastran.converters.cart3d.cart3d_reader import Cart3D

from scipy.spatial import KDTree

class Intersect(object):
    def __init__(self, nodes, elements, regions):
        self.nodes = nodes
        self.elements = elements
        self.regions = regions

    def intersect_tris(self):
        # ==== Calculate Global Edge Length ====
        elements = self.elements - 1
        nodes = self.nodes

        ne, three = elements.shape
        p1 = nodes[elements[:, 0], :]
        p2 = nodes[elements[:, 1], :]
        p3 = nodes[elements[:, 2], :]
        centroid = (p1 + p2 + p3) / 3.

        a = p2 - p1
        b = p3 - p1
        n = cross(a, b)
        assert len(n) == ne, 'len(n)=%s ne=%s' % (len(n), ne)

        print(n)
        ni = norm(n, axis=1)
        print('n.shape=%s ni.shape=%s' % (n.shape, ni.shape))
        assert len(ni) == ne, 'len(ni)=%s ne=%s' % (len(ni), ne)

        A = 0.5 * ni # area
        print(min(ni))
        assert A.min() > 0, A

        #sys.exit()
        n /= ni[:, None]  # normal vector
        assert len(n) == ne, 'len(n)=%s ne=%s' % (len(n), ne)

        # Global Edge Length
        gel = zeros((ne, 2), dtype='float64')
        gel[:, 0] = norm(a, axis=1)
        gel[:, 1] = norm(b, axis=1)

        gel2 = gel.max(axis=1)
        assert len(gel2) == ne, 'len(gel2)=%s ne=%s' % (len(gel2), ne)

        # single valued "Global Edge Length" (hence the i)
        geli = max(gel2)
        print('global_edge_length = %s' % geli)

        # we increase the search size just cause...
        # we're expecting nice isotropic triangles, but aren't totally
        # relying on it
        geli *= 1.05
        print('global_edge_length_i = %s' % geli)

        # ==== create node -> element map ====
        nid_to_eid_map = [[]] * ne
        for eid, (n1, n2, n3) in enumerate(elements):
            nid_to_eid_map[n1].append(eid)
            nid_to_eid_map[n2].append(eid)
            nid_to_eid_map[n3].append(eid)

        # ==== Create KD Tree of centroids ====
        centroid_tree = KDTree(centroid)

        # ==== Intersect All Mesh Geoms ====
        elements2 = []
        for i, element in enumerate(elements):
            c = centroid[i]
            nodes1 = elements[i]
            snodes = set(nodes1)
            gel2i = gel2[i]

            print('c[%i] = %s' % (i, c))
            pts = centroid_tree.query_ball_point(c, gel2i)
            #print(pts)
            for pt in pts:
                diff = norm(c - centroid[pt])
                nodes2 = elements[pt]
                common_set = snodes.intersection(nodes2)
                if not common_set:
                    print('   c[%i]=%s alt[%i]=%s diff=%s gel2=%s valid=%s' % (i, list(nodes1),
                                                                               pt, list(nodes2),
                                                                               diff,
                                                                               gel2[pt], diff < geli))
                    is_intersection = self.intersect(i, pt, nodes1, nodes2, nodes, n)
            #print(centroid_tree.query(c, k=10))
            #break

    def intersect(self, e1, e2, element1, element2, nodes, n):
        """
        http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/pubs/tritri.pdf
        """
        n2 = n[e2]
        #print("nodes.shape =", nodes.shape)
        pt = nodes[element2[0], :]
        d2 = -dot(n2, pt)  # vo2 - node 0 on element 2
        #dvi = []
        #for i in range(3):
            #ei = element1[i]
            #dvii = dot(n2, nodes[ei, :]) + d2
            #dvi.append(dvii)
        #print("    dvi = %s" % dvi)
        #e1 = elements1
        dvi2 = dot(n2, nodes[element1, :].T) + d2

        sdvi = sign(dvi2)
        sign_range = sdvi.max() - sdvi.min()
        if allclose(dvi2.min(), 0.) or sign_range == 2.:
            print("     element2 = ", element2[0])
            print("     ", pt)
            print("     d2", d2)
            print("     dvi = %s" % dvi2)
            print("     sign_range = %s" % sign_range)
            is_intersection = True
            raise NotImplementedError()
        else:
            is_intersection = False


        #print("    n2=%s" % (n2))
        return is_intersection

    def remove_inner_elements(self):
        pass

def intersect_model(cart3d_filename):
    cart3d = Cart3D()
    cart3d.read_cart3d(cart3d_filename)

    intersect = Intersect(cart3d.points, cart3d.elements, cart3d.regions)
    intersect.intersect_tris()

def main():
    cart3d_filename = 'threePlugs_bin.tri'
    intersect_model(cart3d_filename)

if __name__ == '__main__':
    main()
