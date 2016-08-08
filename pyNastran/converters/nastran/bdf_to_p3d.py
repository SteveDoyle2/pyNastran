from __future__ import print_function
from math import radians, degrees, cos, sin
import copy
from pyNastran.bdf.bdf import BDF

def make_circle(radius, theta0, theta1):
    """
    Makes a circle from theta0 to theta1 at a radius R.
    """
    theta0 = radians(theta0)
    theta1 = radians(theta1)
    npoints = 10
    dtheta = (theta1 - theta0) / npoints

    theta = theta0
    X = []
    Y = []
    for i in range(npoints+1):
        x = radius * cos(theta)
        y = radius * sin(theta)
        print("x=%g \ty=%g     \ttheta=%s" %(x, y, degrees(theta)))
        theta += dtheta
        X.append(x)
        Y.append(y)
    return (X, Y)

def read_airfoil(filename):
    with open(filename, 'rU') as infile:
        lines = infile.readlines()

    (nupper, nlower) = lines[1].split()
    nupper = int(float(nupper))
    nlower = int(float(nlower))
    print("nupper=%s nlower=%s" % (nupper, nlower))

    upper_surface = []
    lower_surface = []
    upper_lines = lines[3:3+nupper]
    for line in upper_lines:
        x, y = line.split()
        x = float(x)
        y = float(y)
        upper_surface.append([x, y])

    lower_lines = lines[4+nupper:4+nupper+nlower+1]
    for line in lower_lines:
        x, y = line.split()
        x = float(x)
        y = float(y)
        lower_surface.append([x, y])

    for u in upper_surface:
        print("%g \t%g" %(u[0], u[1]))

    #for u in lower_surface:
        #print("%g \t%g" %(u[0], u[1]))

    return(upper_surface, lower_surface)


class BdfToP3d(object):
    def __init__(self):
        #self.mesh = None
        pass

    def make_mesh(self, fem_name, p3d_name, estart):
        self.mesh = BDF(debug=False)
        #cards_to_include = set(['GRID','CQUAD4'])
        #mesh.setCardsToInclude(cards_to_include)
        self.mesh.read_bdf(fem_name, xref=False)

        (sides, rsides) = self.make_connections()

        edge_b = self.get_common_edge_number(1, 2)
        edge_c = self.get_common_edge_number(1, 51)
        edge_a = self.scale_edge(edge_c + 2)
        edge_d = self.scale_edge(edge_b + 2)
        print("edge_a = ", edge_a)

        edge = self.find_bound_edge(estart)
        #print("edge = ", edge)

        edge_num = self.get_edge_number(estart, edge)
        eid_corner = copy.deepcopy(estart)

        is_bounded2 = True
        eid_corners = []
        saved_ids2 = []
        while is_bounded2:
            eid = copy.deepcopy(eid_corner)
            if eid in eid_corners:
                print("already found eid_corner=%s" % eid)
                break
            eid_corners.append(eid)
            print("eid = %s" % eid)

            saved_ids = []
            is_bounded = True
            while is_bounded:
                saved_ids.append(eid)
                (is_bounded, leftover_eid, shared_edge) = self.get_bound_edge(eid, edge_num)
                #print("leftover_eid = %s" % leftover_eid)
                eid = leftover_eid

            print("")
            print("eLast = %s" % saved_ids[-1])
            (eid_corner, shared_edge) = self.turn_corner(eid_corner, edge_num)
            saved_ids2.append(saved_ids)

            #eid = eid_corner
            #break
        #print("eid_corners = %s" %(eid_corners))

        self.make_plot3d(saved_ids2, estart, p3d_name, edge_d)
        #print("e")

    def get_node(self, nid):
        node = self.mesh.nodes[nid]
        #spot = self.mesh.nodesmap[nid]
        #node = self.mesh.nodes[spot]
        return node


    def get_element(self, eid):
        element = self.mesh.elements[eid]
        return element

    def get_common_node(self, edge_a, edge_b):
        set_a = set(edge_a)
        set_b = set(edge_b)
        intersect = list(set_a.intersection(set_b))
        nid = intersect[0]
        return nid

    def get_edge_number_by_nodes(self, eid, node1, node2):
        edge = [node1, node2]
        edge.sort()
        sides = sides[eid]
        edge_num = sides.index(edge)
        return edge_num

    def get_edge_order(self, eid, edge_a):
        edge1 = edge_a
        edge2 = self.scale_edge(edge_a+1)
        edge3 = self.scale_edge(edge_a+2)
        edge4 = self.scale_edge(edge_a+3)
        return (edge1, edge2, edge3, edge4)

    def make_plot3d(self, elements, estart, p3dname, edge_a):
        (n0, n1, n2, n3) = self.get_edge_order(estart, edge_a)
        element = self.mesh.elements[estart]

        print('n0=%s n1=%s' % (n0, n1))
        #node0 = element.get_node_id(n0)

        jmax = len(elements) + 1
        imax = len(elements[0]) + 1
        kmax = 1
        out = '1\n'
        out += '%s %s %s\n' % (imax, jmax, kmax)

        x = ''
        y = ''
        z = ''
        w = ''
        for j, elements_line in enumerate(elements):
            (x2, y2, z2, w2) = self.write_p3d_line(elements_line, n0, n1) # 0,1
            x += x2
            y += y2
            z += z2
            w += w2

        elements_line = elements[-1]
        (x2, y2, z2, w2) = self.write_p3d_line(elements_line, n3, n2)
        x += x2
        y += y2
        z += z2
        w += w2

        #out += "%s\n%s\n%s\n%s\n" % (x, y, z, w)
        out += '%s\n%s\n%s\n' % (x, y, z)

        with open(p3dname, 'wb') as outfile:
            outfile.write(out)

    def write_p3d_line(self, elements_line, n0, n1):
        #mesh = self.mesh
        x = ''
        y = ''
        z = ''
        w = ''
        print("eLine = ", elements_line)
        espot = elements_line[0]
        element = self.get_element(espot)
        #print("element = ", element)
        print("n0 = ", n0)
        node_ids = element.node_ids
        nid0 = node_ids[n0]
        node0 = self.get_node(nid0)
        #print("node0 = ", node0)
        print("nid = ", nid0)
        x += '%f ' % node0.xyz[0]
        y += '%f ' % node0.xyz[1]
        z += '%f ' % node0.xyz[2]
        w += "%4s " % nid0

        for i, eid in enumerate(elements_line):
            element = self.get_element(eid)
            node_ids = element.node_ids
            print("n1 = ", n1)
            nid1 = node_ids[n1]
            print("nid = ", nid1)
            node1 = self.get_node(nid1)
            x += '%f ' % node1.xyz[0]
            y += '%f ' % node1.xyz[1]
            z += '%f ' % node1.xyz[2]
            w += "%4s " % nid1

            msg = 'nid0=%s i=%s nid0+i=%s nid1=%s' % (nid0, i, nid0+i, nid1)
            assert nid0 + i + 1 == nid1, msg

        x += '\n'
        y += '\n'
        z += '\n'
        w += '\n'
        return (x, y, z, w)

    def get_common_edge_number(self, eid_a, eid_b):
        element_a = self.sides[eid_a]

        set_a = set(element_a)
        set_b = set(self.sides[eid_b])
        edge_list = list(set_a.intersection(set_b))
        edge = edge_list[0]

        edge_num = element_a.index(edge)
        print("common_edge = ", edge)
        print("common_edge_num = ", edge_num)
        return edge_num



    def turn_corner(self, eid, edge_num):
        print("turn_corner eid=%s" % eid)
        dedge_nums = [1, 2, 3]
        for dedge_num in dedge_nums:
            is_bounded, leftover_eid, shared_edge = self.get_bound_edge(eid, edge_num+dedge_num)
            print("bound=%s eid=%s" % (is_bounded, leftover_eid))
            if is_bounded:
                break
        assert is_bounded
        print("leftover_eid = ", leftover_eid)
        print("shared_edge = ", shared_edge)
        return (leftover_eid, shared_edge)


    def get_bound_edge(self, eid, edge_num):
        edge_n = self.get_edge_n(eid, edge_num)
        eshared = self.reversed_sides[edge_n]
        #print("eshared[%s] = %s" % (eid, eshared))
        nshared = len(eshared) - 1

        is_bounded = False
        leftover_eid = None
        if nshared > 0:
            is_bounded = True
            found_index = eshared.index(eid)
            leftover_eid = eshared[1 - found_index]

        return is_bounded, leftover_eid, edge_n

    def scale_edge(self, edge_num):
        if edge_num >= 4:
            edge_num -= 4
        return edge_num

    def get_edge_n(self, eid, edge_num):
        edge_num = self.scale_edge(edge_num)
        #msg = '%s=edge_num<4' % edge_num
        #print("edge_num = ", edge_num)

        #assert edge_num<4,msg
        edges = self.sides[eid]
        return edges[edge_num]

    def find_bound_edge(self, eid_start):
        esides = self.sides[eid_start]
        for side in esides:
            neighboring_elements = self.reversed_sides[side]
            nedges = len(neighboring_elements)-1
            if nedges > 0:
                break
        print("neighboring_elements[%s]=%s" % (side, neighboring_elements))
        print("  nedges = %s" % nedges)
        print("")
        assert nedges > 0
        return side

    def get_edge_number(self, eid, edge):
        #element = getElement(eid)
        side = self.sides[eid]
        edge_num = side.index(edge)
        print("side = ", side)
        print("edge = ", edge)
        print("edge_num = ", edge_num)
        return edge_num

    def make_list(self, nodes, i0, i1):
        node0 = nodes[i0]
        node1 = nodes[i1]
        nlist = [node0, node1]
        return nlist

    def make_tuple(self, nodes, i0, i1):
        nlist = self.make_list(nodes, i0, i1)
        nlist.sort()
        nlist = tuple(nlist)
        return nlist

    def make_connections(self):
        sides = {}
        reversed_sides = {}
        for eid, element in self.mesh.elements.items():
            #print(dir(element))
            nodes = element.node_ids
            #print("element[%s]=%s" %(eid,nodes))
            nlist = self.make_tuple(nodes, 0, -1)
            sides[eid] = [nlist]
            #print("adding nlist = ",nlist)
            reversed_sides.setdefault(nlist, []).append(eid)
            #reversed_sides[nlist].append(eid)
            for i in range(len(nodes)-1):
                nlist = self.make_tuple(nodes, i, i+1)
                #print(reversed_sides)
                #print("adding nlist = ",nlist)
                sides[eid].append(nlist)
                reversed_sides.setdefault(nlist, []).append(eid)

            #print("sides[%s]=%s" %(eid,sides[eid]))
            #if eid>3:
                #break
        print()
        #for side,eid in reversed_sides.items():
            #print("rS[%s]=%s" % (side, eid))
        self.sides = sides
        self.reversed_sides = reversed_sides
        return sides, reversed_sides

def run():
   #(X, Y) = makeCircle(1., 90., 45.)
   #(upper, lower) = readAirfoil('clarky.dat')

    estart = 1
    mesh = BdfToP3d()
    p3d_name = 'mesh2.p3d'
    mesh.make_mesh('airfoil.bdf', p3d_name, estart)


if __name__ == '__main__':
    run()
