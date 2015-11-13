from six import iteritems
from six.moves import zip, range
from numpy import array, dot, norm

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.fieldWriter import print_card
from scipy import spatial


class NastranMesh(BDF):
    def __init__(self):
        BDF.__init__(self)

    def tank_fill(self, mFuel, percent_start=.50, rho_fuel=51.19,
                 tank_elements=None, gravity=None,
                 mass_tol=0.05, nIterMax=10, add_elements=True):
        """
        Fills a single fuel tank in consistent units
        :percent_start: the percentage to start at.
                        percent_start = z/(zMax-zMin); default=0.50
        :mFuel:         mass (or weight) of fuel to add to the tank
        :rhoFuel:       density of fuel (default = 51.19 lb/ft^3 = 6.6 lb/gal)
        :tank_elements: list of elements defining the boundary of the tank
        :gravity:       vector defining the direction of gravity in cid=0; <0,0,-32.2>=default
        :mass_tol:      tolerance on the mass (massError = massTol*mFuel; default=0.05)
        :nIterMax:     the maximum number of iterations (default=10)
        :addElements:  create CONM2 elements if True

        :returns: nodalMasses the masses on each node
        :returns: percentStarts the percentages used in the interpolation
        :returns: massToli = (mFuel-totalMass)/mFuel

        .. note:: massTol should be valid for any tank size
                  (adjust the percent error from 5% if necessary)

        Method:
           1.  Create a new CORD2R coordinate system (cid=-2)
           2.  Rotate the geometry into the gravity coordinate system
           3.  find the z0 (zero fill line) by taking
               z0=percent_start*(zMax-zMin)
           4.  find the nodes below the waterline (z0)
           5.  apply mass to the nodes based on their depth;
               mass = (z-z0)*g*rho
           6.  compare this to the required mass (or weight) for convergence
               test
           7.  interpolate using a spline to get the next point to check

        Requirements:
           1.  No cid=-2 is already being used
           2.  Tank doesnt have to be closed (but it's probably a good thing)
           3.  No requirement on normals
           4.  percentStart < 1.
           5.  massTol < 0.20
        """
        if tank_elements is None:
            tank_elements = []
        assert percent_start < 1.
        assert mass_tol < 0.20
        percent_starts = []
        mass_errors = []
        percent_starts.append(0.)   # x; empty tank
        mass_errors.append(-mFuel)  # y; the mass is too low by mFuel

        # find the vector with the maximum difference with the gravity vector
        if gravity is None:
            gravity = array([0., 0., -32.2])
        #mag_gravity = norm(gravity)  # magnitude
        A = array([1., 0., 0.])  # global x
        B = array([0., 1., 0.])  # global y
        C = array([0., 0., 1.])  # global z

        Adot_gravity = dot(A, gravity)
        Bdot_gravity = dot(B, gravity)
        Cdot_gravity = dot(C, gravity)
        ABC = [A, B, C]
        ABCdot = [Adot_gravity, Bdot_gravity, Cdot_gravity]
        ABCmax = max(ABCdot)
        i = ABCdot.index(ABCmax)

        (bx, by, bz) = ABC[i]
        gx, gy, gz = gravity

        # Create a new CORD2R coordinate system (cid=-2)
        # origin is at [0.,0.,0.] and doesnt matter
        # z axis is the gravity direction (gravity)
        # point on x-z plane is the max unit vector (ABC)
        coordCard = ['CORD2R', -2, 0, 0., 0., 0., gx, gy, gz, bx, by, bz]
        cardName = 'CORD2R'
        self.add_card(self, coordCard, cardName)
        #cardObj = BDF_Card(coordCard)
        #coord = CORD2R(cardObj)

        # convert all points into the gravity frame
        cid = -2
        elementNodeIDs = {}  # CQUAD4, CTRIA3
        node_locations = {}

        for eid in tank_elements:
            elem = self.elements[eid]
            if (elem.type == 'CQUAD4') or (elem.type == 'CTRIA3'):
                nodes = elem.nodes

                elementNodeIDs[eid] = []
                for node in nodes:
                    nid = node.nid
                    elementNodeIDs[eid].append(nid)
                    if nid not in node_locations:
                        p = node.PositionWRT(self, cid)
                        node_locations[nid] = p

        zMax = node_locations[nid][2]
        zMin = zMax
        for nid, node in sorted(iteritems(node_locations)):
            zMax = max(zMax, node[2])
            zMin = min(zMin, node[2])

        # max sure to go into the while looop
        massToli = 2.  # 2*mFuel
        percentFill = percent_start

        nIter = 0
        while massToli > mass_tol and nIter < nIterMax:
            # find the z0 (zero fill line) by taking z0=percentStart*(zMax-zMin)
            z0 = percentFill * (zMax - zMin)

            aboveNodes = set()
            belowNodes = set()
            for nid, node in sorted(iteritems(node_locations)):
                if node[2] >= z0:
                    aboveNodes.add(nid)
                else:
                    belowNodes.add(nid)

            if 0:
                belowElements = set()
                partialElements = set()
                for eid, nodeIDs in sorted(iteritems(elementNodeIDs)):
                    elem = self.elements[eid]

                    isAboveBelow = set()  # True=Above False=Below
                    for nid in nodeIDs:
                        if nid in aboveNodes:
                            isAboveBelow.add(True)
                        else:
                            isAboveBelow.add(False)

                    if   True in isAboveBelow and False not in isAboveBelow:  # all nodes are above
                        pass
                    elif True not in isAboveBelow and False in isAboveBelow:  # all nodes are below
                        belowElements.add(eid)
                    elif True in isAboveBelow and False in isAboveBelow:  # some nodes are above, some below
                        partialElements.add(eid)
                    else:
                        raise RuntimeError('not above, not below, not partial...')

            if 0:
                for eid in belowElements:
                    elem = self.elements[eid]
                    nodeIDs = elementNodeIDs[eid]
                    if elem.type == 'CQUAD4':
                        pass
                    elif elem.type == 'CTRIA3':
                        pass
                    else:
                        raise NotImplementedError()

            # compute the total and elemental masses
            nodalMasses = {}
            totalMass = 0.
            for nid in belowNodes:
                # mass = g*rho*(z0-z)
                # it's (z0-z) b/c magGravity is always positive and z0 is higher than z
                mass = mag_gravity * (z0 - node_locations[nid][2])
                nodalMasses[nid] = mass
                totalMass += mass
            massError = mFuel - totalMass

            percent_starts.append(percent_start)  # x
            mass_errors.append(massError)         # y
            mass_toli = massError / mFuel

            #x=[]; y=[]
            for xi, yi in mass_found:
                x.append(xi)
                y.append(yi)
            i = argsort(x)
            X = array(percent_starts)[i]  # sorted x
            Y = array(mass_toli)[i]       # sorted y

            spline = buildSpline(Y, X)  # reverse interpolation
            yi = 0.  # find 0. mass
            xi = splev(yi, spline)  # the percentFill for 0. mass
            percentFill = xi

            nIter += 1

        if add_elements:
            maxEid = max(self.elements) + 1  # get the next available eid
            for (nid, mass) in sorted(iteritems(nodalMasses)):
                card = ['CONM2', maxEid, nid, 0, mass]
                self.add_card(self, card, 'CONM2')
                maxEid += 1

        del self.coords[cid]
        return masses, X, Y

    def create_plane(self, width, height, nx, ny, eid_start):
        """
        Creates a quadrilateral plane made of CTRIA3s.  This output
        can be intersected with another geometry.

        :param width:     width of plane
        :param height:    height of plane
        :param nx:        number of elements along the width
        :param ny:        number of elements along the height
        :param eid_start: the starting elementID

        ::

          1-----2
          | \ B |
          | A \ |
          4-----3
        """
        assert nx > 0
        assert ny > 0
        dx = width / nx
        dy = height / ny

        nidStart = 200
        n = nidStart
        #x=[]; y=[]
        ij_NMap = {}
        points = {}
        for j in range(ny + 1):
            yi = dy * j
            for i in range(nx + 1):
                xi = dx * i
                points[n] = array([xi, yi, 0.])
                ij_NMap[(i, j)] = n
                n += 1

        elements = []
        for j in range(ny):
            for i in range(nx):
                element = [ij_NMap[(i, j)],
                           ij_NMap[(i + 1, j)],
                           ij_NMap[(i + 1, j + 1)], ]
                elements.append(element)
                element = [ij_NMap[(i, j)],
                           ij_NMap[(i + 1, j + 1)],
                           ij_NMap[(i, j + 1)], ]
                elements.append(element)

        #origin = [width/2,height/2,0.]
        origin = [0., 0., -1.]
        zAxis = array([0., 1., 0.])
        xAxis = array([1., 0., 0.])

        f = open('plane.bdf', 'wb')
        #f.write('SOL 101\n')
        #f.write('CEND\n')
        #f.write('BEGIN BULK\n')

        cid = 200
        rid = 0
        coord = ['CORD2R', cid, rid, origin[0], origin[1], origin[2],
                 zAxis[0], zAxis[1], zAxis[2],
                 xAxis[0], xAxis[1], xAxis[2]]
        f.write(print_card(coord))
        for nid, point in sorted(iteritems(points)):
            (x, y, z) = point
            node = ['GRID', nid, cid, x, y, z]
            f.write(print_card(node))

        pid = 123
        mid = pid
        eidStart = 200
        for eid, element in enumerate(elements):
            (n1, n2, n3) = element
            tri = ['CTRIA3', eid + eidStart, pid, n1, n2, n3]
            f.write(print_card(tri))

        shell = ['PSHELL', pid, mid, 1.0]
        f.write(print_card(shell))
        mat = ['MAT1', mid, 1e7, None, 0.3]
        f.write(print_card(mat))
        #f.write('ENDDATA\n')
        f.close()
        return points, elements

    def intersect(self, eids, maxDist=2.):
        nodeIs = {}
        nodeLocations = {}
        i = 0
        nodes = []
        for nid, node in sorted(iteritems(self.nodes)):
            position = node.get_position()
            nodeLocations[nid] = position
            nodes.append(position)
            nodeIs[i] = nid
            i += 1

        originalElements = {}
        newElements = {}
        newNodes = set()
        oldNodes = set()
        for eid, element in sorted(iteritems(self.elements)):
            if eid in eids:
                newElements[eid] = element
                newNodes = newNodes.union(set(element.node_ids))
            else:
                originalElements[eid] = element
                oldNodes = oldNodes.union(set(element.node_ids))

        #for eid in newElements:
        newNodes = list(newNodes)
        print("newNodes = ", sorted(newNodes))
        #print("newElements.keys = ", newElements.keys())
        #print("originalElements.keys = ", originalElements.keys())

        tree = spatial.KDTree(nodes)
        k = 10

        for nid in newNodes:
            positionLookup = nodeLocations[nid]

            (dists, iLocs) = tree.query(array(positionLookup), k=k)
            #print "iLocs=%s dists=%s" %(iLocs,dists)
            #print "out = ",out

            closeNodes = set()
            for dist, iLoc in zip(dists, iLocs):
                closeNodes.add(nodeIs[iLoc])
                #print "nodeClose[%s] = %s" %(nodeIs[iLoc],nodeLocations[nodeIs[iLoc]])
            intersectionNodes = closeNodes.intersection(oldNodes)
            if intersectionNodes:
                print("intersectionNodes[%s] = %s" % (nid, intersectionNodes))


    # def flip_normals(self, starterEid, eids=None, flipStarter=False):
    #     """
    #     Takes the normals of SHELL elements and flips it to a common direction
    #     This method follows the contour of the body, so assuming
    #     no internal elements, all the normals on the outside will point
    #     outwards (or inwards).
    #
    #     :param starterEid:  the element to copy the normal of
    #     :param eids:        the element IDs to flip to the common direction (default=None -> all)
    #     :param flipStarter: should the staring element be flipped (default=False)
    #
    #     .. todo:: finish method...think i need to build a edge list...
    #               that'd be a lot easier to loop through stuff...
    #     """
    #     raise NotImplementedError()
    #     normals = {}
    #     validNids = set([])
    #     isCorrectNormal = set([])
    #
    #     allEids = eids
    #     if allEids is None:
    #         allEids = self.elements.keys()
    #     setAllEids = set(allEids)
    #
    #     if flipStarter:
    #         elem = self.Element(starterEid)
    #         elem.flipNormal()
    #     normals[starterEid] = elem.Normal()
    #
    #     for eid in allEids:
    #         element = self.elements[eid]
    #         if isinstance(element, ShellElement):
    #             elem = self.Element(starterEid)
    #             normals[starterEid] = elem.Normal()
    #             validNids = validNids.union(set(elem.node_ids))
    #
    #     # clean up the elements that will be considered
    #     elemsToCheck = set([])
    #     nidToEidMap = self.getNodeIDToElementIDsMap()
    #     for (nid, eidsMap) in sorted(iteritems(nidToEidMap)):
    #         if nid not in validNids:  # clean up extra nodes
    #             del nidToEidMap[nid]
    #         else:
    #             eids = list(set(eids))  # do i need this?
    #             for eid in eids:  # clean up ROD/SOLID elements
    #                 eids2 = []
    #                 if eid in setAllEids:
    #                     eids2.append(eid)
    #                 elemsToCheck = elemsToCheck.union(set(eids2))
    #             nidToEidMap[nid] = eids2
    #
    #     # starts with the starter element, loops thru adjacent elements
    #     # and checks to see if the normal is 'close' to the elements
    #     # normal from before
    #     goEid = starterEid
    #
    #     # no recursion to avoid recursion limit
    #     while 1:
    #         elem = self.Element(goEid)
    #         nids = elem.getNodeIDs()
    #         normals = self._get_adjacent_normals(nids, nidToEidMap)
    #         normal = normals[goEid]
    #
    # def _get_adjacent_elements(self, nids, nidToEidMap):
    #     """
    #     .. todo:: doesnt work...
    #     """
    #     raise NotImplementedError()
    #     normals = {}
    #     #for nid in

if __name__ == '__main__':  ## pragma: no cover
    mesh = NastranMesh()
    if 0:
        width = 10.
        height = 20.
        nx = 10
        ny = 10
        eidStart = 10
        mesh.createPlane(width, height, nx, ny, eidStart)
    mesh.read_bdf('combo.bdf')

    # plane eids
    eids = [200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212,
            213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225,
            226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238,
            239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251,
            252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264,
            265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277,
            278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290,
            291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303,
            304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316,
            317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329,
            330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342,
            343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355,
            356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368,
            369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381,
            382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394,
            395, 396, 397, 398, 399, ]
    mesh.intersect(eids)
