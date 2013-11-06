## GNU Lesser General Public License
##
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
##
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
##
## This file is part of pyNastran.
##
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
##
# pylint: disable=R0904,R0902
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import multiprocessing as mp

from numpy import array, cross, zeros, dot


from pyNastran.bdf.cards.loads.staticLoads import Moment, Force
#from pyNastran.bdf.cards.elements.shell import ShellElement


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


def _mass_properties_mass_mp_func(element):
    try:
        p = element.Centroid()
        mass = element.Mass()
    except:
        mass = 0.
        cg = array([0., 0., 0.])
    return mass, p


class BDFMethods(BDFMethodsDeprecated):
    def __init__(self):
        pass

    def mass_properties(self, reference_point=None, sym_axis=None, num_cpus=1):
        """
        Caclulates mass properties in the global system about the reference point.
        :param self: the BDF object
        :param reference_point: an array that defines the origin of the frame.
            default = <0,0,0>.
        :returns mass: the mass of the model
        :returns cg: the cg of the model as an array.
        :returns I: moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz])

        I = mass * centroid * centroid

        .. math:: I_{xx} = m (dy^2 + dz^2)

        .. math:: I_{yz} = -m * dy * dz

        where:
        .. math:: dx = x_{element} - x_{ref}

        .. seealso:: http://en.wikipedia.org/wiki/Moment_of_inertia#Moment_of_inertia_tensor
        """
        if reference_point is None:
            reference_point = array([0., 0., 0.])

        if num_cpus > 1:
            mass, cg, I = self._mass_properties_mp(num_cpus, reference_point=reference_point, sym_axis=sym_axis)
        else:
                     #Ixx Iyy Izz, Ixy, Ixz Iyz
            I = array([0., 0., 0., 0., 0., 0., ])
            cg = array([0., 0., 0.])
            mass = 0.
            # precompute the CG location and make it the reference point
            if reference_point == 'cg':
                for element in self.elements.itervalues():
                    try:
                        p = element.Centroid()
                        m = element.Mass()
                        mass += m
                        cg += m * p
                    except:
                        pass
                if mass == 0.0:
                    return (mass, cg, I)
                reference_point = cg / mass

            cg = array([0., 0., 0.])
            mass = 0.
            for element in self.elements.itervalues():
                try:
                    p = element.Centroid()
                except:
                    if element.type not in ['CBUSH']:
                        print('****', element.type)
                        print(str(element))
                        raise
                    continue

                try:
                    m = element.Mass()
                    (x, y, z) = p - reference_point
                    x2 = x * x
                    y2 = y * y
                    z2 = z * z
                    I[0] += m * (y2 + z2)  # Ixx
                    I[1] += m * (x2 + z2)  # Iyy
                    I[2] += m * (x2 + y2)  # Izz
                    I[3] += m * x * y      # Ixy
                    I[4] += m * x * z      # Ixz
                    I[5] += m * y * z      # Iyz
                    mass += m
                    cg += m * p
                except:
                    if element.type not in ['CBUSH']:
                        raise
                    self.log.warning("could not get the inertia for element"
                                     "...\n%s" % element)
                    raise
            if mass:
                cg = cg / mass

        if sym_axis == None:
            for key, aero in self.aero.iteritems():
                print("aero =", str(aero))
                sym_axis = ''
                if aero.IsSymmetricalXY():
                    sym_axis += 'y'
                if aero.IsSymmetricalXZ():
                    sym_axis += 'z'
                if IsAntiSymmetricalXY():
                    raise NotImplementedError('%s is antisymmetric about the XY plane' % str(aero))
                if IsAntiSymmetricalXZ():
                    raise NotImplementedError('%s is antisymmetric about the XZ plane' % str(aero))

        if None is not sym_axis and 'x' in sym_axis:
            #print("symx")
            I[0] *= 2.0
            I[1] *= 2.0
            I[2] *= 2.0
            I[3] *= 0.0  # Ixy
            I[4] *= 0.0  # Ixz
            I[5] *= 2.0  # Iyz
            cg[0] = 0.0
        if None is not sym_axis and 'y' in sym_axis:
            #print("symy")
            I[0] *= 2.0
            I[1] *= 2.0
            I[2] *= 2.0
            I[3] *= 0.0  # Ixy
            I[4] *= 2.0  # Ixz
            I[5] *= 0.0  # Iyz
            cg[1] = 0.0
        if None is not sym_axis and 'z' in sym_axis:
            #print("symz")
            I[0] *= 2.0
            I[1] *= 2.0
            I[2] *= 2.0
            I[3] *= 2.0  # Ixy
            I[4] *= 0.0  # Ixz
            I[5] *= 0.0  # Iyz
            cg[2] = 0.0
        return (mass, cg, I)


    def _mass_properties_mp(self, num_cpus, reference_point=None, sym_axis=None):
        """
        Caclulates mass properties in the global system about the reference point.
        :param self: the BDF object
        :param reference_point: an array that defines the origin of the frame.
            default = <0,0,0>.
        :returns mass: the mass of the model
        :returns cg: the cg of the model as an array.
        :returns I: moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz])

        I = mass * centroid * centroid

        .. math:: I_{xx} = m (dy^2 + dz^2)
        .. math:: I_{yz} = -m * dy * dz

        where:
        .. math:: dx = x_{element} - x_{ref}

        .. seealso:: http://en.wikipedia.org/wiki/Moment_of_inertia#Moment_of_inertia_tensor
        """
        if num_cpus <= 1:
            raise RuntimeError('num_proc must be > 1; num_cpus=%s' % num_cpus)

        nelements = len(self.elements)
        #-----------------------------------------------------------
        self.log.info("Creating %i-process pool!" % num_cpus)

        pool = mp.Pool(num_cpus)
        result = pool.imap(_mass_properties_mass_mp_func, [(element) for element in self.elements.itervalues() if element.type not in ['CBUSH'] ])
        #result = [_mass_properties_mass_mp_func(element) for element in self.elements.itervalues()]

        mass = zeros((nelements), 'float64')
        xyz = zeros((nelements, 3), 'float64')
        for j, return_values in enumerate(result):
            #self.log.info("%.3f %% Processed"% (j*100./nelements))
            mass[j] = return_values[0]
            xyz[j, :] = return_values[1]
        self.log.info("Shutting down process pool!")
        pool.close()
        pool.join()

        massi = mass.sum()

        #print('xyz.shape =', xyz.shape)
        #print('mass.shape =', mass.shape)

        #cg = (mass * xyz) / massi
        if massi == 0.0:
            cg = array([0., 0., 0.])
            I = array([0., 0., 0., 0., 0., 0., ])
            return massi, cg, I

        cg = dot(mass, xyz) / massi
        #cg = numpy.multiply(mass, xyz) / massi
        #cg = numpy.multiply(xyz, mass) / massi

        #cg = numpy.multiply(xyz, mass).sum(axis=0) / massi
        #cg = numpy.dot(xyz.T, mass).T / massi

        #cg = (xyz * mass).sum(axis=0) / massi  # (186, 1)
        #print('cg =', cg)
        #print('reference_point ', reference_point)
        if reference_point is None:
            x = xyz[:, 0]
            y = xyz[:, 1]
            z = xyz[:, 2]
        elif isinstance(reference_point[0], float):
            x = xyz[:, 0] - reference_point[0]
            y = xyz[:, 1] - reference_point[1]
            z = xyz[:, 2] - reference_point[2]
        elif reference_point in [u'cg', 'cg']:
            x = xyz[:, 0] - cg[0]
            y = xyz[:, 1] - cg[1]
            z = xyz[:, 2] - cg[2]

        x2 = x**2
        y2 = y**2
        z2 = z**2

        #A = y2 + z2
        #print('mass.shape', mass.shape)
        #print('A.shape', A.shape)
        I = array([
           mass * (y2 + z2),  # Ixx
           mass * (x2 + z2),  # Iyy
           mass * (x2 + y2),  # Izz
           mass * (x * y),    # Ixy
           mass * (x * z),    # Ixz
           mass * (y * z),    # Iyz
        ]).sum(axis=1)

        return (massi, cg, I)

        #massi = numpy.sum(mass)
        #cgi = numpy.sum(cg, axis=0)
        #del cg #, mass
        #self.log.info("massi = %s" % massi)
        #------------------------------------------------------------

        #mass = 0.
        #cg = array([0., 0., 0.])
        #I = zeros((nelements, 5), 'float64')
        #result = [_inertia_mp_func(mass[i], element) for element in xrange(self.elements.iterkeys()]
        #for j, return_values in enumerate(result):
            #I[j, :] += return_values[2]
        #I = numpy.sum(I, axis=0)
        #return (massi, cg, I)

    def mass(self):
        """Calculates mass in the global coordinate system"""
        mass = 0.
        for element in self.elements.itervalues():
            try:
                m = element.Mass()
                mass += m
            except:
                self.log.warning("could not get the mass for element"
                                 "...\n%s" % element)
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
    #     # clean up the elements that will be considered
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

    def resolve_grids(self, cid=0):
        """
        Puts all nodes in a common coordinate system (mainly for cid testing)

        :param self: the object pointer
        :param cid:  the cid to resolve the nodes to (default=0)
        .. note:: loses association with previous coordinate systems so to go
                  back requires another fem
        """
        assert cid in self.coords, ('cannot resolve nodes to '
                                    'cid=%r b/c it doesnt exist' % cid)
        for nid, node in sorted(self.nodes.iteritems()):
            p = node.PositionWRT(self, cid)
            node.UpdatePosition(self, p, cid)

    def unresolve_grids(self, model_old):
        """
        Puts all nodes back to original coordinate system.

        :param self:      the object pointer
        :param model_old: the old model that hasnt lost it's connection to
                          the node cids
        .. warning:: hasnt been tested well...
        """
        debug = False
        for (nid, node_old) in model_old.nodes.iteritems():
            coord = node_old.cp
            (p, matrix) = coord.transformToGlobal(self.xyz, debug=debug)
            p2 = coord.transformToLocal(p, matrix, debug=debug)
            self.nodes[nid].UpdatePosition(self, p2, coord.cid)

    def sum_forces(self):
        """
        Sums applied forces for all load cases.
        Considers FORCE, FORCE1, FORCE2.

        :returns Forces: the forces as a numpy array
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
        Sums applied forces & moments about a reference point p0 for all
        load cases.
        Considers FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2.

        :param p0:        the reference point
        :returns Moments: the moments as a numpy array
        :returns Forces:  the forces as a numpy array
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