from __future__ import print_function
import os

from six.moves import zip
from numpy import array, cross, dot, abs as npabs
from numpy.linalg import norm, solve #, cond

from pyNastran.applications.cart3d_nastran_fsi.math_functions import (
    is_list_ranged, list_print, shepard_weight)
#from pyNastran.applications.cart3d_nastran_fsi.matTest import fIsNatural
from pyNastran.bdf.field_writer_8 import print_card_8

from pyNastran.utils.log import get_logger

log = get_logger(log=None, level='debug', encoding='utf-8')

#------------------------------------------------------------------

class Tet4(object):
    def __init__(self, p0, p1, p2, p3, ID=None, nodes=None, neighbors=None):
        """
        The internal tets dont need a number or list of nodes.
        The Delauney TETs do!
        """
        if neighbors is None:
            neighbors = []
        self.p0 = array(p0)
        self.p1 = array(p1)
        self.p2 = array(p2)
        self.p3 = array(p3)

        self.ID = ID
        self.nodes = nodes
        self.neighbors = neighbors

        #self.face0 = (self.p0 + self.p1 + self.p2)/3
        #self.face1 = (self.p0 + self.p1 + self.p3)/3
        #self.face2 = (self.p0 + self.p2 + self.p3)/3
        #self.face3 = (self.p1 + self.p2 + self.p3)/3
        #print("face0 = ", self.face0)

        #print("p0 = ", p0)
        #print("p1 = ", p1)
        #print("p2 = ", p2)
        #print("p3 = ", p3)
        self._centroid = self.calculate_centroid()
        self._volume = None

        #vol = self.volume()
        #vol2 = self.testVol(self.p0, self.p1,
        #                    self.p2, self.p3)
        #assert vol == vol2
        #msg = "volume[%s]=%s" % (ID, vol)
        #assert vol>0,msg
        #log.info(msg)


    def __repr__(self):
        nodes = "[%3s, %3s, %3s, %3s]" % (
            self.nodes[0], self.nodes[1], self.nodes[2], self.nodes[3])
        return "ID=%4s nodes=%s p0=%s p1=%s p2=%s p3=%s" % (
            self.ID, nodes, self.p0, self.p1, self.p2, self.p3)

    def calculate_centroid(self):
        """
        centroid = (a+b+c)/4  with point d at the origin
        or more generally...
        centroid = ((a-d)+(b-d)+(c-d))/4 = (a+b+c-3d)/4

        based on http://en.wikipedia.org/wiki/Tetrahedron
        """
        c = (self.p0 + self.p1 + self.p2 - 3 * self.p3) / 4
        return c

    def centroid(self):
        return self._centroid

    def abs_volume(self):
        return npabs(self.volume())

    def test_vol(self, a, b, c, d):
        vol = dot(a-d, cross(b-d, c-d)) / 6.
        return vol

    def volume(self):
        """
        volume = dot(a-d, cross(b-d, c-d))/6
        based on http://en.wikipedia.org/wiki/Tetrahedron
        """
        if self._volume is None:
            self._volume = dot(self.p0-self.p3, cross(self.p1-self.p3, self.p2-self.p3))/6.
        else:
            return self._volume
        return self._volume

    def map_deflections2(self, deflections, aero_node):
        #A =
        pass

    def map_deflections(self, deflections, aero_node):
        """
        determines the mapped x,y,z deflections of the aero node

        Parameters
        ----------
        deflections : ???
            deflections at p0,p1,p2,p3
        aero_node : ???
            node to find the deflections at

        Solves Ax=b  where A is the nodes, x is the abc coeffs,
                     b is the deflection at the tet node points (x,y,z)
        then uses x to find the deflection at an arbitrary point m (aero_node)
        """
        A = array([[1.]+list(self.p0),
                   [1.]+list(self.p1),
                   [1.]+list(self.p2),
                   [1.]+list(self.p3)])
        #print("m1 = ", m)
        (d1, d2, d3, d4) = deflections
        log.info('d1=%s' % d1)
        log.info('d2=%s' % d2)
        log.info('d3=%s' % d3)
        log.info('d4=%s' % d4)
        #log.info("A = \n%s" % (A))
        #condA = cond(A,2)
        #log.info("A cond=%g;  A=\n%s" % (condA, print_matrix(A)))
        #print("ID  = ",self.ID)

        #condMax = 1E6  # max allowable condition number before alternate approach
        if 1:
            nodes = [self.p0, self.p1, self.p2, self.p3]
            weights = shepard_weight(aero_node, nodes)
            du = array([0., 0., 0.])
            for (node, weight) in zip(nodes, weights):
                du += weight * node
        else:
            ux = self.map_deflection(A, d1, d2, d3, d4, 0, aero_node, 'x')
            uy = self.map_deflection(A, d1, d2, d3, d4, 1, aero_node, 'y')
            uz = self.map_deflection(A, d1, d2, d3, d4, 2, aero_node, 'z')
            du = array([ux, uy, uz])
        log.info("vol       = %s" % self.volume())
        log.info("du        = %s" % du)
        return aero_node + du

    def map_deflection(self, A, d1, d2, d3, d4, i, aero_node, dType='x'):
        """
        see mapDeflections...
        based on the posiition matrix, finds the ith deflection component

        Then we ratios the new deflection diMax and compares it to the nearby points to 'test'.
        """
        di = array([d1[i], d2[i], d3[i], d4[i]])
        diMax = max(npabs(di)) # L1 norm
        log.info("d%s = %s" % (dType, list_print(di)))

        abcd = solve(A, di)

        #ui = abcd*aero_node # element-wise multiplication...faster,
        # cant make it work; tried dot & transpose too...
        (a, b, c, d) = abcd
        ui = a + b*aero_node[0] + c*aero_node[1] + d*aero_node[2]
        is_ranged = is_list_ranged(0., abcd, 1.)
        log.info('is_ranged=%s  u%sRatio=%g - a=%g b=%g c=%g d=%g' %(
            is_ranged, dType, npabs(ui/diMax), a, b, c, d))
        return ui

    def is_internal_node(self, pAero):
        r"""
        If there is an internal node, than the sum of the internal volumes without
        an absolute value will equal the volume of the local tet
          *   3
         /  \
        *---*  0 1
         \*/   2
        """
        is_contained, gammas = fIsNatural(pAero, self.p0, self.p1, self.p2, self.p3,
                                          V=self._volume)
        return is_contained, min(gammas)

        #a = Tet4(self.p0, self.p1, self.p2, pAero,ID='a',nodes=[0,1,2,'p'])
        #b = Tet4(self.p0, self.p1, pAero,self.p3, ID='b',nodes=[0,1,'p',3])
        #c = Tet4(pAero,self.p0, self.p2, self.p3, ID='c',nodes=['p',0,2,3])
        #d = Tet4(self.p1,pAero, self.p2, self.p3, ID='d',nodes=[1,'p',2,3])

        #aVol = a.volume()
        #if aVol >= 0:
        #    bVol = b.volume()
        #    if bVol >= 0:
        #        cVol = c.volume()
        #        if cVol >= 0.:
        #            dVol = d.volume()
        #            if dVol >= 0:
        #                return True
        #return False

        #aVol = a.volume()
        #bVol = b.volume()
        #cVol = c.volume()
        #dVol = d.volume()
        #allv = [aVol, bVol, cVol, dVol]
        #minThis = sum(allv) - max(allv)

        #print("vol[%s] = [%g %g %g %g]   %g" %(self.ID, aVol, bVol, cVol, dVol, minThis))
        #minVol = min(aVol, bVol, cVol, dVol)
        #if minVol < 0.:
            #return False, minThis
        #return True, minThis

#------------------------------------------------------------------

class DelauneyReader(object):
    def __init__(self, infilename):
        self.infilename = infilename

    def build_tets(self):
        grids, elements, volume = self.read()

        tets = {}
        for key, element_neighbors in elements.items():
            nodes, neighbors = element_neighbors
            id1, id0, id2, id3 = nodes

            grid0 = grids[id0]
            grid1 = grids[id1]
            grid2 = grids[id2]
            grid3 = grids[id3]

            tet = Tet4(grid0, grid1, grid2, grid3,
                       ID=key, nodes=nodes, neighbors=neighbors)
            tets[key] = tet
            #tets.append(tet)
        return tets, grids, elements, volume

    def distance(self, p1, p2):
        #print("p1=", p1)
        #print("p2=", p2)
        return norm(p1 - p2)

   #def testVolume(self):
       #vol = 0.
       #for tet in tets:
           #vol += tet.volume()
       #print("volume=%s vol=%s" % (volume, -vol/2.))

    def test_mapper(self):
        m = array([0., 0., 0.])
        tets, grids, elements, volume = self.build_tets()
        tet_new = self.find_closest_tet(m, tets)
        self.update_node(m, tets)
        return tet_new

    def update_node(self, m, tets):
        """
        Takes in a new point m (mid tet point) and a list of tets
        from self.buildTets(), finds the tet the node should go in,
        then updates the node based on the deflections at the local
        tet.
        """
        #m = array([21.18, 0.5, -1.87e-6])
        tet_new = self.find_closest_tet(m, tets)

        # def = deflections at p0,p1,p2,p3
        new_aero_node = tet_new.map_deflections(deflections, m)
        return new_aero_node

    def find_closest_tet(self, m, tets):
        """
        Finds the closest tet.  Has the potential to fail, but until then, skipping.
        Solution to failure:  brute force method.

        m = aeroNode
        tets =
        """
        starting_tet = tets[0]
        closest_tet = self.find_closest_tet_recursion(m, starting_tet, tets)
        #print("found tet = ", closestTet.ID)
        #v1 = array([1., 0., 0.])
        return closest_tet


    def find_closest_tet_recursion(self, m, tet0, tets, excluded=None):
        """
        Makes an assumption that there is a direct line from the starting tet
        to the final tet that can be minimized with nearly every subsequent tet.
        At worst, it will get the starting node one tet closer to the solution.

        m = starting point (mid point)
        """
        if excluded is None:
            excluded = []
        #print("working on tet = ",tet0.ID)
        if tet0.is_internal_node(m):
            return tet0

        cent = tet0.centroid()
        dist = m - cent
        #faces = [tet0.face0, tet0.face1, tet0.face2, tet0.face3]
        dists = [9.e9] * 4
        for i, neighbor in enumerate(tet0.neighbors):
            if neighbor > 0:
                dists[i] = self.distance(m, tets[i].centroid())
        #dists[0] = 9.e9
        #print("dists = ", dists)

        min_value = min(dists)
        i = dists.index(min_value)
        tet_new = tets[tet0.neighbors[i]]
        excluded.append(tet0.ID)
        #print("excluding ID=%s\n" % tet0.ID)

        closest_tet = self.find_closest_tet_recursion(m, tet_new, tets, excluded)
        return closest_tet

    def read(self):
        """
        reads all the node points, neighbors, and the IDs tet nodes,
        determines the volume for testing.
        """
        log.info("---reading Delauney Tetrahedralization file...%r" % self.infilename)
        with open(self.infilename, 'r') as infile:
            lines = infile.readlines()

        slines = []
        for line in lines:
            sline = line.strip().split()
            if sline:
                slines.append(sline)

        #log.info('slines[0] = %s' % (slines[0]))
        ngrids, nelements = self.ints(slines[0])
        log.info("ngrid=%s  nelements=%s" % (ngrids, nelements))

        grids = {}
        for sline in slines[1:ngrids+1]:
            node = self.node(sline)
            grids[node[0]] = node[1:]
            #print("id=%s grid=%s" %(node[0], node[1:]))

        max_x = grids[1][0]
        max_y = grids[1][1]
        max_z = grids[1][2]

        min_x = grids[1][0]
        min_y = grids[1][1]
        min_z = grids[1][2]

        for key, grid in grids.items():
            max_x = max(max_x, grid[0])
            max_y = max(max_y, grid[1])
            max_z = max(max_z, grid[2])

            min_x = min(min_x, grid[0])
            min_y = min(min_y, grid[1])
            min_z = min(min_z, grid[2])

        elements = {}
        nelements = 1
        for sline in slines[ngrids+1:]:
            e = self.ints(sline)
            #print(e)
            neighbors = e[1:5]
            #print("e[1:4] = ", element1)
            element = e[5:]
            #print("e[3:] = ", element2)
            #print(len(element1))
            #print(len(element2))
            elements[nelements] = [element, neighbors]
            nelements += 1

            #print("e[%s]=%s" %(e[0], e[1:]))

        log.info("max_x=%s min_x=%s" % (max_x, min_x))
        log.info("max_y=%s min_y=%s" % (max_y, min_y))
        log.info("max_z=%s min_z=%s" % (max_z, min_z))
        dx = max_x - min_x
        dy = max_y - min_y
        dz = max_z - min_z
        volume = dx * dy * dz
        log.info("volume = %s" % volume)

        return grids, elements, volume

    def write_tets_as_nastran(self, nodes, tets, tet_filename='tet.bdf'):
        with open(tet_filename, 'wb') as outfile:
            msg = ''
            mid = 1 # material id
            for nid, node in nodes.items():
                card = ['GRID', ID, 0, node[0], node[1], node[2]]
                msg += print_card_8(card)
            outfile.write(msg)

            msg = ''
            for i, tet in enumerate(tets):
                nodes = tet.nodes
                msg += print_card_8(['CTETRA', tet.ID, mid, nodes[0], nodes[1], nodes[2], nodes[3]])
                if i % 1000 == 0:
                    outfile.write(msg)
                    msg = ''
            outfile.write(msg)
        log.info("finished writing %s" % tet_filename)

    def ints(self, values):
        return [int(value) for value in values]

    def floats(self, values):
        return [float(value) for value in values]

    def node(self, values):
        return [int(values[0])] + self.floats(values[1:])

def run():
    infilename = os.path.join('delauney', 'geometry.morph.in')
    d = DelauneyReader(infilename)
    tets, nodes, elements, volume = d.build_tets()

    m = array([21.18, 0.5, -1.87347e-6])
    d.find_closest_tet(m, tets)

    d.write_tets_as_nastran(nodes, tets)
    #d.read()
    #d.testVolume()

if __name__ == '__main__':
    run()
