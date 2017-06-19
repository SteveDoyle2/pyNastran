from __future__ import print_function
from six import iteritems
from six.moves import zip, range
from numpy import array, dot

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.field_writer import print_card
from scipy import spatial


def tank_fill(model, mfuel, percent_start=0.50, rho_fuel=51.19,
              tank_elements=None, gravity=None,
              mass_tol=0.05, niter_max=10, add_elements=True):
    """
    Fills a single fuel tank in consistent units

    Parameters
    -----------
    model : BDF()
        a BDF object
    percent_start : float
        the percentage to start at
        percent_start = z/(zMax-zMin); default=0.50
    mfuel : float
        mass (or weight) of fuel to add to the tank
    rho_fuel : float
        density of fuel (default = 51.19 lb/ft^3 = 6.6 lb/gal)
    tank_elements : List[int]
        list of elements defining the boundary of the tank
    gravity : float
        vector defining the direction of gravity in cid=0; <0,0,-32.2>=default
    mass_tol : float
        tolerance on the mass (massError = massTol*mfuel; default=0.05)
    niter_max : int; default=10
        the maximum number of iterations
    add_elements : bool; default=True
        create CONM2 elements if True

    Returns
    -------
    nodal_masses : ???
        the masses on each node
    percentStarts : float
        the percentages used in the interpolation
    massToli : float
        massToli = (mfuel - total_mass)/mfuel

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
       4.  percent_start < 1.
       5.  massTol < 0.20
    """
    if tank_elements is None:
        tank_elements = []
    assert percent_start < 1.
    assert mass_tol < 0.20
    percent_starts = []
    mass_errors = []
    percent_starts.append(0.)   # x; empty tank
    mass_errors.append(-mfuel)  # y; the mass is too low by mFuel

    # find the vector with the maximum difference with the gravity vector
    if gravity is None:
        gravity = array([0., 0., -32.2])
    #mag_gravity = norm(gravity)  # magnitude
    A = array([1., 0., 0.])  # global x
    B = array([0., 1., 0.])  # global y
    C = array([0., 0., 1.])  # global z

    a_dot_gravity = dot(A, gravity)
    b_dot_gravity = dot(B, gravity)
    c_dot_gravity = dot(C, gravity)
    abc = [A, B, C]
    abc_dot = [a_dot_gravity, b_dot_gravity, c_dot_gravity]
    abc_max = max(abc_dot)
    i = abc_dot.index(abc_max)

    (bx, by, bz) = abc[i]
    gx, gy, gz = gravity

    # Create a new CORD2R coordinate system (cid=-2)
    # origin is at [0.,0.,0.] and doesnt matter
    # z axis is the gravity direction (gravity)
    # point on x-z plane is the max unit vector (ABC)
    model.add_cord2r(-2, 0,
                     origin=[0., 0., 0.],
                     zaxis=[gx, gy, gz],
                     xzplane=[bx, by, bz])
    #card_obj = BDF_Card(coord_card)
    #coord = CORD2R(card_obj)

    # convert all points into the gravity frame
    cid = -2
    element_node_ids = {}  # CQUAD4, CTRIA3
    node_locations = {}

    for eid in tank_elements:
        elem = model.elements[eid]
        if elem.type in ('CQUAD4', 'CTRIA3'):
            nodes = elem.nodes

            element_node_ids[eid] = []
            for node in nodes:
                nid = node.nid
                element_node_ids[eid].append(nid)
                if nid not in node_locations:
                    p = node.PositionWRT(model, cid)
                    node_locations[nid] = p

    zMax = node_locations[nid][2]
    zMin = zMax
    for nid, node in sorted(iteritems(node_locations)):
        zMax = max(zMax, node[2])
        zMin = min(zMin, node[2])

    # max sure to go into the while looop
    mass_toli = 2.  # 2*mFuel
    percent_fill = percent_start

    niter = 0
    while mass_toli > mass_tol and niter < niter_max:
        # find the z0 (zero fill line) by taking z0=percentStart*(zMax-zMin)
        z0 = percent_fill * (zMax - zMin)

        above_nodes = set()
        below_nodes = set()
        for nid, node in sorted(iteritems(node_locations)):
            if node[2] >= z0:
                above_nodes.add(nid)
            else:
                below_nodes.add(nid)

        if 0:
            below_elements = set()
            partial_elements = set()
            for eid, node_ids in sorted(iteritems(element_node_ids)):
                elem = model.elements[eid]

                is_above_below = set()  # True=Above False=Below
                for nid in node_ids:
                    if nid in above_nodes:
                        is_above_below.add(True)
                    else:
                        is_above_below.add(False)

                if   True in is_above_below and False not in is_above_below:
                    # all nodes are above
                    pass
                elif True not in is_above_below and False in is_above_below:
                    # all nodes are below
                    below_elements.add(eid)
                elif True in is_above_below and False in is_above_below:
                    # some nodes are above, some below
                    partial_elements.add(eid)
                else:
                    raise RuntimeError('not above, not below, not partial...')

        if 0:
            for eid in below_elements:
                elem = model.elements[eid]
                node_ids = element_node_ids[eid]
                if elem.type == 'CQUAD4':
                    pass
                elif elem.type == 'CTRIA3':
                    pass
                else:
                    raise NotImplementedError()

        # compute the total and elemental masses
        nodal_masses = {}
        total_mass = 0.
        for nid in below_nodes:
            # mass = g*rho*(z0-z)
            # it's (z0-z) b/c magGravity is always positive and z0 is higher than z
            mass = mag_gravity * (z0 - node_locations[nid][2])
            nodal_masses[nid] = mass
            total_mass += mass
        mass_error = mfuel - total_mass

        percent_starts.append(percent_start)  # x
        mass_errors.append(mass_error)         # y
        mass_toli = mass_error / mfuel

        #x=[]; y=[]
        for xi, yi in mass_found:
            x.append(xi)
            y.append(yi)
        i = argsort(x)
        X = array(percent_starts)[i]  # sorted x
        Y = array(mass_toli)[i]       # sorted y

        spline = build_spline(Y, X)  # reverse interpolation
        yi = 0.  # find 0. mass
        xi = splev(yi, spline)  # the percent_fill for 0. mass
        percent_fill = xi

        niter += 1

    if add_elements:
        eid = max(model.elements) + 1  # get the next available eid
        for (nid, mass) in sorted(iteritems(nodal_masses)):
            model.add_conm2(eid, nid, 0, mass)
            eid += 1

    del model.coords[cid]
    return masses, X, Y

def create_plane(model, width, height, nx, ny, eid_start):
    r"""
    Creates a quadrilateral plane made of CTRIA3s.  This output
    can be intersected with another geometry.

    Parameters
    ----------
    model : BDF()
        a BDF object
    width : float
        width of plane
    height : float
        height of plane
    nx : int
        number of elements along the width
    ny : int
        number of elements along the height
    eid_start : int
        the starting elementID

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

    nid_start = 200
    n = nid_start
    ij_nmap = {}
    points = {}
    for j in range(ny + 1):
        yi = dy * j
        for i in range(nx + 1):
            xi = dx * i
            points[n] = array([xi, yi, 0.])
            ij_nmap[(i, j)] = n
            n += 1

    elements = []
    for j in range(ny):
        for i in range(nx):
            element = [ij_nmap[(i, j)],
                       ij_nmap[(i + 1, j)],
                       ij_nmap[(i + 1, j + 1)], ]
            elements.append(element)
            element = [ij_nmap[(i, j)],
                       ij_nmap[(i + 1, j + 1)],
                       ij_nmap[(i, j + 1)], ]
            elements.append(element)

    #origin = [width/2,height/2,0.]
    origin = [0., 0., -1.]
    zaxis = array([0., 1., 0.])
    xaxis = array([1., 0., 0.])

    with open('plane.bdf', 'wb') as bdf_file:
        #bdf_file.write('SOL 101\n')
        #bdf_file.write('CEND\n')
        #bdf_file.write('BEGIN BULK\n')

        cid = 200
        rid = 0
        coord = ['CORD2R', cid, rid, origin[0], origin[1], origin[2],
                 zaxis[0], zaxis[1], zaxis[2],
                 xaxis[0], xaxis[1], xaxis[2]]
        bdf_file.write(print_card(coord))
        for nid, point in sorted(iteritems(points)):
            (x, y, z) = point
            node = ['GRID', nid, cid, x, y, z]
            bdf_file.write(print_card(node))

        pid = 123
        mid = pid
        eid_start = 200
        for eid, element in enumerate(elements):
            (n1, n2, n3) = element
            tri = ['CTRIA3', eid + eid_start, pid, n1, n2, n3]
            bdf_file.write(print_card(tri))

        shell = ['PSHELL', pid, mid, 1.0]
        bdf_file.write(print_card(shell))
        mat = ['MAT1', mid, 1e7, None, 0.3]
        bdf_file.write(print_card(mat))
        #bdf_file.write('ENDDATA\n')
    return points, elements

def intersect(model, eids, max_dist=2.):
    nodeIs = {}
    node_locations = {}
    i = 0
    nodes = []
    for nid, node in sorted(iteritems(model.nodes)):
        position = node.get_position()
        node_locations[nid] = position
        nodes.append(position)
        nodeIs[i] = nid
        i += 1

    original_elements = {}
    new_elements = {}
    new_nodes = set()
    old_nodes = set()
    for eid, element in sorted(iteritems(model.elements)):
        if eid in eids:
            new_elements[eid] = element
            new_nodes = new_nodes.union(set(element.node_ids))
        else:
            original_elements[eid] = element
            old_nodes = old_nodes.union(set(element.node_ids))

    #for eid in new_elements:
    new_nodes = list(new_nodes)
    print("new_nodes = ", sorted(new_nodes))
    #print("new_elements.keys = ", new_elements.keys())
    #print("original_elements.keys = ", original_elements.keys())

    tree = spatial.KDTree(nodes)
    k = 10

    for nid in new_nodes:
        position_lookup = node_locations[nid]

        (dists, ilocs) = tree.query(array(position_lookup), k=k)
        #print("ilocs=%s dists=%s" % (ilocs, dists))
        #print("out = ", out)

        close_nodes = set()
        for dist, iloc in zip(dists, ilocs):
            close_nodes.add(nodeIs[iloc])
            #print "nodeClose[%s] = %s" % (nodeIs[iloc], node_locations[nodeIs[iloc]])
        intersection_nodes = close_nodes.intersection(old_nodes)
        if intersection_nodes:
            print("intersection_nodes[%s] = %s" % (nid, intersection_nodes))


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
#         elem.flip_normal()
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

def main():  # pragma: no cover
    model = BDF()
    if 0:
        width = 10.
        height = 20.
        nx = 10
        ny = 10
        eid_start = 10
        create_plane(model, width, height, nx, ny, eid_start)
    model.read_bdf('combo.bdf')

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

if __name__ == '__main__':  # pragma: no cover
    main()
