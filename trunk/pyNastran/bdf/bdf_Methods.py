"""
This file contains additional methods that do not directly relate to the
reading/writing/accessing of BDF data.  Such methods include:
  - Mass
      get the mass of the model
  - Mass Poperties
      get the mass & moment of inertia of the model
  - sumMoments / sum_moments
      find the net force/moment on the model
  - sumForces / sum_forces
      find the net force on the model
  - resolve_grids
      change all nodes to a specific coordinate system
  - unresolve_grids
      puts all nodes back to original coordinate system
"""
# pylint: disable=R0904,R0902
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import multiprocessing as mp

from numpy import array, cross, zeros, dot
from numpy.linalg import det


from pyNastran.bdf.cards.loads.staticLoads import Moment, Force, LOAD
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
        cg = element.Centroid()
        mass = element.Mass()
    except:
        cg = array([0., 0., 0.])
        mass = 0.
    return mass, cg


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
                for pack in [self.elements, self.masses]:
                    for element in pack.itervalues():
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
            for pack in [self.elements, self.masses]:
                for element in pack.itervalues():
                    try:
                        p = element.Centroid()
                    except:
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
                        self.log.warning("could not get the inertia for element\n%s" % element)
                        continue
            if mass:
                cg = cg / mass

        if sym_axis is None:
            for key, aero in self.aero.iteritems():
                sym_axis = ''
                if aero.IsSymmetricalXY():
                    sym_axis += 'y'
                if aero.IsSymmetricalXZ():
                    sym_axis += 'z'
                if aero.IsAntiSymmetricalXY():
                    raise NotImplementedError('%s is antisymmetric about the XY plane' % str(aero))
                if aero.IsAntiSymmetricalXZ():
                    raise NotImplementedError('%s is antisymmetric about the XZ plane' % str(aero))
        if sym_axis is not None:
            self.log.debug('Mass/MOI sym_axis = %r' % sym_axis)

        scale = 1.0
        if 'WTMASS' in self.params:
            scale = self.params['WTMASS'].values[0]
            #print("mass scale =", scale)

        if None is not sym_axis and 'x' in sym_axis:
            mass *= 2.0
            I[0] *= 2.0
            I[1] *= 2.0
            I[2] *= 2.0
            I[3] *= 0.0  # Ixy
            I[4] *= 0.0  # Ixz
            I[5] *= 2.0  # Iyz
            cg[0] = 0.0
        if None is not sym_axis and 'y' in sym_axis:
            mass *= 2.0
            I[0] *= 2.0
            I[1] *= 2.0
            I[2] *= 2.0
            I[3] *= 0.0  # Ixy
            I[4] *= 2.0  # Ixz
            I[5] *= 0.0  # Iyz
            cg[1] = 0.0
        if None is not sym_axis and 'z' in sym_axis:
            mass *= 2.0
            I[0] *= 2.0
            I[1] *= 2.0
            I[2] *= 2.0
            I[3] *= 2.0  # Ixy
            I[4] *= 0.0  # Ixz
            I[5] *= 0.0  # Iyz
            cg[2] = 0.0

        # these are strangely not included in the F06...
        #mass *= scale
        #I *= scale ** 2.0
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

        nelements = len(self.elements) + len(self.masses)
        #-----------------------------------------------------------
        self.log.debug("Creating %i-process pool!" % num_cpus)

        pool = mp.Pool(num_cpus)
        result  = pool.imap(_mass_properties_mass_mp_func, [(element) for element in self.elements.itervalues()
                           if element.type not in ['CBUSH', 'CBUSH1D',
                               'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                               'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
                           ]])
        result2 = pool.imap(_mass_properties_mass_mp_func, [(element) for element in self.masses.itervalues() ])

        mass = zeros((nelements), 'float64')
        xyz = zeros((nelements, 3), 'float64')
        for j, return_values in enumerate(result):
            #self.log.info("%.3f %% Processed"% (j*100./nelements))
            mass[j] = return_values[0]
            xyz[j, :] = return_values[1]
        pool.close()
        pool.join()

        pool = mp.Pool(num_cpus)
        for j2, return_values in enumerate(result2):
            #self.log.info("%.3f %% Processed"% (j*100./nelements))
            mass[j+j2] = return_values[0]
            xyz[j+j2, :] = return_values[1]
        pool.close()
        pool.join()

        massi = mass.sum()

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

        x2 = x ** 2
        y2 = y ** 2
        z2 = z ** 2

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
        for pack in [self.elements, self.masses]:
            for element in pack.itervalues():
                try:
                    m = element.Mass()
                    mass += m
                except:
                    self.log.warning("could not get the mass for element"
                                     "...\n%s" % element)
        return mass

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

    def __gravity_load(loadcase_id):
        gravity_i = model.loads[2][0]
        gi = gravity_i.N * gravity_i.scale
        mass, cg, I = model.mass_properties(reference_point=p0, sym_axis=None, num_cpus=6)

    def sum_forces_moments(self, p0, load_case_id, include_grav=False):
        """
        Sums applied forces & moments about a reference point p0 for all
        load cases.
        Considers:
          - FORCE, FORCE1, FORCE2
          - MOMENT, MOMENT1, MOMENT2
          - PLOAD2, PLOAD4
          - LOAD

        :param p0:
          the reference point
        :returns Forces:
          the forces
        :returns Moments:
          the moments

        :type p0:
          NUMPY.NDARRAY shape=(3,) or integer
        :type Forces:
          NUMPY.NDARRAY shape=(3,)
        :type Moments:
          NUMPY.NDARRAY shape=(3,)

        ..warning:: not validated
        ..todo:: It's super slow for cid != 0.   We can speed this up a lot
                 if we calculate the normal, area, centroid based on
                 precomputed node locations.
        """
        if not isinstance(load_case_id, int):
            raise RuntimeError('load_case_id must be an integer; load_case_id=%r' % load_case_id)
        if isinstance(p0, int):
            p = self.model.nodes[p0].Position()
        else:
            p = array(p0)

        loadCase = self.loads[load_case_id]
        #for (key, loadCase) in self.loads.iteritems():
            #if key != load_case_id:
                #continue

        scale_factors2 = []
        loads2 = []
        for load in loadCase:
            if isinstance(load, LOAD):
                scale_factors, loads = load.getReducedLoads()
                scale_factors2 += scale_factors
                loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)

        F = array([0., 0., 0.])
        M = array([0., 0., 0.])

        xyz = {}
        for nid, node in self.nodes.iteritems():
            xyz[nid] = node.Position()

        unsupported_types = set([])
        for load, scale in zip(loads2, scale_factors2):
            if isinstance(load, Force):  # FORCE, FORCE1, FORCE2
                f = load.mag * load.xyz
                node = self.Node(load.node)
                r = xyz[node.nid] - p
                m = cross(r, f)
                F += f * scale
                M += m * scale
            elif isinstance(load, Moment):  # MOMENT, MOMENT1, MOMENT2
                m = load.mag * load.xyz
                M += m * scale
            elif load.type == 'PLOAD1':
                pass
                #elem = self.elements[load.eid]
                #if elem.type in ['CBAR',]:
                    #pass
            elif load.type == 'PLOAD2':
                p = load.pressures[0] * scale  # there are 4 pressures, but we assume p0
                for eid in load.eids:
                    elem = self.elements[eid]
                    if elem.type in ['CTRIA3',
                                     'CQUAD4', 'CSHEAR']:
                        n = elem.Normal()
                        A = elem.Area()
                        f = p * n * A
                        r = elm.Centroid() - p
                        m = cross(r, f)
                        F += f
                        M += m
                    else:
                        self.log.debug('case=%s etype=%r loadtype=%r not supported' % (load_case_id, elem.type, load.type))
            elif load.type == 'PLOAD4':
                #elem = load.eid
                p = load.pressures[0] * scale  # there are 4 possible pressures, but we assume p0
                assert load.Cid() == 0, 'Cid() = %s' % (load.Cid())
                assert load.sorl == 'SURF', 'sorl = %s' % (load.sorl)
                assert load.ldir == 'NORM', 'ldir = %s' % (load.ldir)

                for elem in load.eids:
                    eid = elem.eid
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',]:
                        # triangles
                        nodes = elem.nodeIDs()
                        n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                        axb = cross(n1 - n2, n1 - n3)
                        nunit = det(axb)
                        A = 0.5 * nunit
                        n = axb / nunit
                        centroid = (n1 + n2 + n3) / 3.

                        n2 = elem.Normal()
                        assert n == n2
                    elif elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        # quads
                        nodes = elem.nodeIDs()
                        n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                        axb = cross(n1 - n3, n2 - n4) # buggy
                        nunit = det(axb)
                        A = 0.5 * nunit
                        n = axb / nunit
                        centroid = (n1 + n2 + n3 + n4) / 4.

                        n2 = elem.Normal()
                        assert n == n2
                    elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
                        A, centroid, normal = elem.getFaceAreaCentroidNormal(load.g34.nid, load.g1.nid)
                    else:
                        self.log.debug('case=%s eid=%s etype=%r loadtype=%r not supported' % (load_case_id, eid, elem.type, load.type))
                        continue
                    r = centroid - p
                    f = p * A * n
                    #load.cid.transformToGlobal()
                    m = cross(r, f)
                    F += f
                    M += m
            elif load.type == 'GRAV':
                if include_grav:  # this will be super slow
                    assert load.Cid() == 0, 'Cid() = %s' % (load.Cid())
                    g = load.N * load.scale * scale
                    for eid, elem in self.elements.iteritems():
                        cenntroid = elem.Centroid()
                        mass = elem.Mass()
                        r = centroid - p
                        f = mass * g
                        m = cross(r, f)
                        F += f
                        M += m
            else:
                # we collect them so we only get one print
                unsupported_types.add(load.type)

        for Type in unsupported_types:
            self.log.debug('case=%s loadtype=%r not supported' % (load_case_id, Type))
        #self.log.info("case=%s F=%s M=%s\n" % (load_case_id, F, M))
        return (F, M)
