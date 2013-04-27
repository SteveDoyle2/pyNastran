# pylint: disable=R0904,R0902
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import array, cross

from pyNastran.bdf.cards.loads.staticLoads import Moment, Force
from pyNastran.bdf.cards.elements.shell import ShellElement


class BDFMethodsDeprecated(object):
    def MassProperties(self):
        """
        .. seealso:: mass_properties
        .. deprecated: will be replaced in version 0.7 with mass_properties
        """
        return self.mass_properties()

    def Mass(self):
        """
        .. seealso:: mass
        .. deprecated: will be replaced in version 0.7 with mass
        """
        return self.mass()

    def resolveGrids(self, cid=0):
        """
        .. seealso:: resolve_grids
        .. deprecated: will be replaced in version 0.7 with resolve_grids
        """
        return self.resolve_grids(cid)

    def unresolveGrids(self, femOld):
        """
        .. seealso:: unresolve_grids
        .. deprecated: will be replaced in version 0.7 with unresolve_grids
        """
        return self.unresolve_grids(femOld)

    def sumForces(self):
        """
        .. seealso:: sum_forces
        .. deprecated: will be replaced in version 0.7 with sum_forces
        """
        return self.sum_forces()

    def sumMoments(self, p0):
        """
        .. seealso:: sum_moments
        .. deprecated: will be replaced in version 0.7 with sum_moments
        """
        return self.sum_moments(p0)


class BDFMethods(BDFMethodsDeprecated):
    def __init__(self):
        pass

    def mass_properties(self):
        """
        Caclulates mass properties in the global system about <0,0,0>
        I   = mass*centroid*centroid
        Ixx =  mass*(y^2+z^2)
        Iyz = -mass*y*z
        http://en.wikipedia.org/wiki/Moment_of_inertia#Moment_of_inertia_tensor
        """
                 #Ixx Iyy Izz, Ixy, Ixz Iyz
        I = array([0., 0., 0., 0., 0., 0., ])
        cg = array([0., 0., 0.])
        mass = 0.
        for element in self.elements.itervalues():
            try:
                p = element.Centroid()  # not really coded across the board
                m = element.Mass()
                (x, y, z) = p
                x2 = x * x
                y2 = y * y
                z2 = z * z
                I[0] += m * (y2 + z2)  # Ixx
                I[1] += m * (x2 + z2)  # Iyy
                I[2] += m * (x2 + y2)  # Izz
                I[3] -= m * x * y      # Ixy
                I[4] -= m * x * z      # Ixz
                I[5] -= m * y * z      # Iyz
                mass += m
                cg += m * p
            except:
                self.log().warning("could not get inertia for element"
                                   "...\n%s" % (element))
        cg = cg / mass
        return (mass, cg, I)

    def mass(self):
        """Calculates mass in the global coordinate system"""
        mass = 0.
        for element in self.elements.itervalues():
            m = element.Mass()
            mass += m
        return mass

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
    #             validNids = validNids.union(set(elem.nodeIDs()))
    #
    #     ## clean up the elements that will be considered
    #     elemsToCheck = set([])
    #     nidToEidMap = self.getNodeIDToElementIDsMap()
    #     for (nid, eidsMap) in sorted(nidToEidMap.iteritems()):
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
    #     ## starts with the starter element, loops thru adjacent elements
    #     ## and checks to see if the normal is 'close' to the elements
    #     ## normal from before
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

    def resolve_grids(self, cid=0):
        """
        Puts all nodes in a common coordinate system (mainly for cid testing)

        :param self: the object pointer
        :param cid:  the cid to resolve the nodes to
        .. note:: loses association with previous coordinate systems so to go
                  back requires another fem
        """
        assert cid in self.coords, ('cannot resolve nodes to '
                                    'cid=|%s| b/c it doesnt exist' % (cid))
        for nid, node in sorted(self.nodes.iteritems()):
            p = node.PositionWRT(self, cid)
            node.UpdatePosition(self, p, cid)

    def unresolve_grids(self, fem_old):
        """
        Puts all nodes back to original coordinate system.

        :param self:       the object pointer
        :param fem_old:    the old model that hasnt lost it's connection to
                           the node cids
        .. warning:: hasnt been tested well...
        """
        debug = False
        for (nid, node_old) in fem_old.nodes.iteritems():
            coord = node_old.cp
            (p, matrix) = coord.transformToGlobal(self.xyz, debug=debug)
            p2 = coord.transformToLocal(p, matrix, debug=debug)
            self.nodes[nid].UpdatePosition(self, p2, coord.cid)

    def sum_forces(self):
        """
        Sums applied forces for all load cases.
        Considers FORCE, FORCE1, FORCE2.

        @retval Forces the forces as a numpy array
        .. warning:: not validated
        """
        for (key, loadCase) in self.loads.iteritems():
            F = array([0., 0., 0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load, Force):
                    f = load.mag * load.xyz
                    F += f
            self.log.info("case=%s F=%s\n\n" % (key, F))
        return F

    def sum_moments(self, p0):
        """
        Sums applied forces & moments about a reference point p0
        for all load cases.
        Considers FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2.

        :p0: the reference point
        @retval Moments the moments as a numpy array
        @retval Forces the forces as a numpy array
        ..warning:: not validated
        """
        p = array(p0)
        for (key, loadCase) in self.loads.iteritems():
            M = array([0., 0., 0.])
            F = array([0., 0., 0.])
            #print "loadCase = ",loadCase
            for load in loadCase:
                #print "load = ",load
                if isinstance(load, Force):
                    f = load.mag * load.xyz
                    node = self.Node(load.node)
                    r = node.Position() - p
                    m = cross(r, f)
                    M += m
                    F += f
                elif isinstance(load, Moment):
                    m = load.mag * load.xyz
                    M += m
            print("case=%s F=%s M=%s\n\n" % (key, F, M))
        return (M, F)