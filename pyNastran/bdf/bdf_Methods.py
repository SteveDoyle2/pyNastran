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
from six import iteritems, integer_types
from six.moves import zip, range
import multiprocessing as mp

from numpy import array, cross, zeros, dot
from numpy.linalg import norm


from pyNastran.bdf.cards.loads.staticLoads import Moment, Force, LOAD


def _mass_properties_mass_mp_func(element):
    """helper method for mass properties multiprocessing"""
    try:
        cg = element.Centroid()
        mass = element.Mass()
    except:
        cg = array([0., 0., 0.])
        mass = 0.
    return mass, cg


class BDFMethods(object):
    def __init__(self):
        pass

    def mass_properties(self, element_ids=None, reference_point=None,
                        sym_axis=None, num_cpus=1, scale=None):
        """
        Caclulates mass properties in the global system about the
        reference point.

        Parameters
        ---------------
        self : BDF
            The BDF object.
        element_ids : ndarray, optional
            An array of element ids.
        reference_point : ndarray, optional
            An array that defines the origin of the frame.
            default = <0,0,0>.
        sym_axis : str, optional
            The axis to which the model is symmetric. If AERO cards are used, this can be left blank
            allowed_values = 'x', 'y', 'z', 'xy', 'yz', 'xz', 'xyz'
        scale : float, optional
            The WTMASS scaling value.
            default=None -> PARAM, WTMASS is used
            float > 0.0

        Returns
        ----------
        mass : float
            The mass of the model.
        cg : ndarray
            The cg of the model as an array.
        I : ndarray
            Moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz]).

        I = mass * centroid * centroid

        .. math:: I_{xx} = m (dy^2 + dz^2)

        .. math:: I_{yz} = -m * dy * dz

        where:

        .. math:: dx = x_{element} - x_{ref}

        .. seealso:: http://en.wikipedia.org/wiki/Moment_of_inertia#Moment_of_inertia_tensor

        .. note::
           This doesn't use the mass matrix formulation like Nastran.
           It assumes m*r^2 is the dominant term.
           If you're trying to get the mass of a single element, it
           will be wrong, but for real models will be correct.
        """
        if reference_point is None:
            reference_point = array([0., 0., 0.])

        if element_ids is None:
            elements = self.elements.values()
            masses = self.masses.values()
            nelements = len(self.elements) + len(self.masses)
        else:
            elements = [element for eid, element in self.elements.items() if eid in element_ids]
            masses = [mass for eid, mass in self.masses.items() if eid in element_ids]
            nelements = len(element_ids)

        num_cpus = 1
        if num_cpus > 1:
            mass, cg, I = self._mass_properties_mp(num_cpus, elements, masses,
                                                   nelements,
                                                   reference_point=reference_point)
        else:
            mass, cg, I = self._mass_properties_sp(elements, masses,
                                                   reference_point=reference_point)

        mass, cg, I = self._apply_mass_symmetry(sym_axis, scale, mass, cg, I)
        return (mass, cg, I)

    def _mass_properties_sp(self, elements, masses, reference_point):
        #Ixx Iyy Izz, Ixy, Ixz Iyz
        # precompute the CG location and make it the reference point
        I = array([0., 0., 0., 0., 0., 0., ])
        cg = array([0., 0., 0.])
        if reference_point == 'cg':
            mass = 0.
            for pack in [elements, masses]:
                for element in pack:
                    try:
                        p = element.Centroid()
                        m = element.Mass()
                        mass += m
                        cg += m * p
                    except:
                        pass
            if mass == 0.0:
                return mass, cg, I

            reference_point = cg / mass
        else:
            # reference_point = [0.,0.,0.] or user-defined array
            pass

        mass = 0.
        cg = array([0., 0., 0.])
        for pack in [elements, masses]:
            for element in pack:
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
        return (mass, cg, I)

    def _apply_mass_symmetry(self, sym_axis, scale, mass, cg, I):
        """
        Scales the mass & moement of inertia based on the symmetry axes
        and the PARAM WTMASS card
        """
        if sym_axis is None:
            for key, aero in iteritems(self.aero):
                sym_axis = ''
                if aero.is_symmetric_xy():
                    sym_axis += 'y'
                if aero.is_symmetric_xz():
                    sym_axis += 'z'
                if aero.is_anti_symmetric_xy():
                    raise NotImplementedError('%s is anti-symmetric about the XY plane' % str(aero))
                if aero.is_anti_symmetric_xz():
                    raise NotImplementedError('%s is anti-symmetric about the XZ plane' % str(aero))
        if sym_axis is not None:
            # either we figured sym_axis out from the AERO cards or the user told us
            self.log.debug('Mass/MOI sym_axis = %r' % sym_axis)

        if None is not sym_axis:
            if 'x' in sym_axis:
                mass *= 2.0
                I[0] *= 2.0
                I[1] *= 2.0
                I[2] *= 2.0
                I[3] *= 0.0  # Ixy
                I[4] *= 0.0  # Ixz
                I[5] *= 2.0  # Iyz
                cg[0] = 0.0
            if 'y' in sym_axis:
                mass *= 2.0
                I[0] *= 2.0
                I[1] *= 2.0
                I[2] *= 2.0
                I[3] *= 0.0  # Ixy
                I[4] *= 2.0  # Ixz
                I[5] *= 0.0  # Iyz
                cg[1] = 0.0
            if 'z' in sym_axis:
                mass *= 2.0
                I[0] *= 2.0
                I[1] *= 2.0
                I[2] *= 2.0
                I[3] *= 2.0  # Ixy
                I[4] *= 0.0  # Ixz
                I[5] *= 0.0  # Iyz
                cg[2] = 0.0

        if scale is None and 'WTMASS' in self.params:
            scale = self.params['WTMASS'].values[0]
            self.log.info('WTMASS scale = %r' % scale)
        elif scale is None:
            scale = 1.0
        mass *= scale
        I *= scale
        return (mass, cg, I)


    def _mass_properties_mp(self, num_cpus, elements, masses, nelements,
                            reference_point=None):
        """
        Caclulates mass properties in the global system about the
        reference point.

        :param self: the BDF object
        :param num_cpus: the number of CPUs to use; 2 < num_cpus < 20
        :param reference_point: an array that defines the origin of the frame.
            default = <0,0,0>.
        :returns mass: the mass of the model
        :returns cg: the cg of the model as an array.
        :returns I: moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz])

        .. seealso:: self.mass_properties
        """
        if num_cpus <= 1:
            raise RuntimeError('num_proc must be > 1; num_cpus=%s' % num_cpus)
        if num_cpus > 20:
            # the user probably doesn't want 68,000 CPUs; change it if you want...
            raise RuntimeError('num_proc must be < 20; num_cpus=%s' % num_cpus)

        self.log.debug("Creating %i-process pool!" % num_cpus)
        pool = mp.Pool(num_cpus)
        result = pool.imap(_mass_properties_mass_mp_func,
                           [(element) for element in elements
                            if element.type not in ['CBUSH', 'CBUSH1D',
                                                    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                                                    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
                                                    ]])
        result2 = pool.imap(_mass_properties_mass_mp_func, [(element) for element in masses])

        mass = zeros((nelements), 'float64')
        xyz = zeros((nelements, 3), 'float64')
        i = 0
        for i, return_values in enumerate(result):
            #self.log.info("%.3f %% Processed" % (i*100./nelements))
            mass[i] = return_values[0]
            xyz[i, :] = return_values[1]
        pool.close()
        pool.join()

        pool = mp.Pool(num_cpus)
        for i2, return_values in enumerate(result2):
            mass[i+i2] = return_values[0]
            xyz[i+i2, :] = return_values[1]
        pool.close()
        pool.join()

        massi = mass.sum()

        #cg = (mass * xyz) / massi
        if massi == 0.0:
            cg = array([0., 0., 0.])
            I = array([0., 0., 0., 0., 0., 0., ])
            return massi, cg, I

        cg = dot(mass, xyz) / massi
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

        I = array([
            mass * (y2 + z2),  # Ixx
            mass * (x2 + z2),  # Iyy
            mass * (x2 + y2),  # Izz
            mass * (x * y),    # Ixy
            mass * (x * z),    # Ixz
            mass * (y * z),    # Iyz
        ]).sum(axis=1)

        return (massi, cg, I)

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
        for nid, node in sorted(iteritems(self.nodes)):
            p = node.get_position_wrt(self, cid)
            node.set_position(self, p, cid)

    def unresolve_grids(self, model_old):
        """
        Puts all nodes back to original coordinate system.

        :param self:      the object pointer
        :param model_old: the old model that hasnt lost it's connection to
                          the node cids

        .. warning:: hasnt been tested well...
        """
        debug = False
        for (nid, node_old) in iteritems(model_old.nodes):
            coord = node_old.cp
            (p, matrix) = coord.transformToGlobal(self.xyz, debug=debug)
            p2 = coord.transformToLocal(p, matrix, debug=debug)
            self.nodes[nid].UpdatePosition(self, p2, coord.cid)

    def __gravity_load(self, loadcase_id):
        """
        .. todo::
            1.  resolve the load case
            2.  grab all of the GRAV cards and combine them into one
                GRAV vector
            3.  run mass_properties to get the mass
            4.  multiply by the gravity vector
        """

        gravity_i = self.loads[2][0]  ## .. todo:: hardcoded
        gi = gravity_i.N * gravity_i.scale
        p0 = array([0., 0., 0.])  ## .. todo:: hardcoded
        mass, cg, I = self.mass_properties(reference_point=p0, sym_axis=None,
                                           num_cpus=6)

    def sum_forces_moments_elements(self, p0, loadcase_id, eids, nids,
                                    include_grav=False):
        """
        Sum the forces/moments based on a list of nodes and elements.

        :param eids:  the list of elements to include (e.g. the loads
                      due to a PLOAD4)
        :param nids:  the list of nodes to include (e.g. the loads due
                      to a FORCE card)
        :param p0:    the point to sum moments about
                      type = int
                          sum moments about the specified grid point
                      type = (3, ) ndarray/list (e.g. [10., 20., 30]):
                          the x, y, z location in the global frame
        Nodal Types  : FORCE, FORCE1, FORCE2,
                       MOMENT, MOMENT1, MOMENT2,
                       PLOAD
        Element Types: PLOAD1, PLOAD2, PLOAD4, GRAV

        If you have a CQUAD4 (eid=3) with a PLOAD4 (eid=3) and a FORCE
        card (nid=5) acting on it, you can incldue the PLOAD4, but
        not the FORCE card by using:

        For just pressure:

        .. code-block:: python

          eids = [3]
          nids = []

        For just force:

        .. code-block:: python

          eids = []
          nids = [5]

        or both:

        .. code-block:: python

          eids = [3]
          nids = [5]

        .. note:: If you split the model into sections and sum the loads
                  on each section, you may not get the same result as
                  if you summed the loads on the total model.  This is
                  due to the fact that nodal loads on the boundary are
                  double/triple/etc. counted depending on how many breaks
                  you have.

        .. todo:: not done...
        """
        if not isinstance(loadcase_id, integer_types):
            raise RuntimeError('loadcase_id must be an integer; loadcase_id=%r' % loadcase_id)
        if isinstance(p0, integer_types):
            p = self.nodes[p0].get_position()
        else:
            p = array(p0)

        load_case = self.loads[loadcase_id]
        #for (key, load_case) in iteritems(self.loads):
            #if key != loadcase_id:
                #continue

        scale_factors2 = []
        loads2 = []
        for load in load_case:
            if isinstance(load, LOAD):
                scale_factors, loads = load.get_reduced_loads()
                scale_factors2 += scale_factors
                loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)

        F = array([0., 0., 0.])
        M = array([0., 0., 0.])

        xyz = {}
        for nid, node in iteritems(self.nodes):
            xyz[nid] = node.get_position()

        unsupported_types = set([])
        for load, scale in zip(loads2, scale_factors2):
            if isinstance(load, Force):  # FORCE, FORCE1, FORCE2
                if load.nodeIDs not in nids:
                    continue
                f = load.mag * load.xyz
                node = self.Node(load.node)
                r = xyz[node.nid] - p
                m = cross(r, f)
                F += f * scale
                M += m * scale
            elif isinstance(load, Moment):  # MOMENT, MOMENT1, MOMENT2
                if load.nodeIDs not in nids:
                    continue
                m = load.mag * load.xyz
                M += m * scale

            elif load.type == 'PLOAD':
                nodes = load.nodeIDs()
                nnodes = len(nodes)
                nodesi = 0
                if nnodes == 3:
                    n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                    axb = cross(n1 - n2, n1 - n3)
                    centroid = (n1 + n2 + n3) / 3.

                elif nnodes == 4:
                    n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                    axb = cross(n1 - n3, n2 - n4)
                    centroid = (n1 + n2 + n3 + n4) / 4.
                    if nodes[3] in nids:
                        nodesi += 1
                else:
                    raise RuntimeError('invalid number of nodes on PLOAD card; nodes=%s' % str(nodes))
                if nodes[0] in nids:
                    nodesi += 1
                if nodes[1] in nids:
                    nodesi += 1
                if nodes[2] in nids:
                    nodesi += 1

                nunit = norm(axb)
                A = 0.5 * nunit
                try:
                    n = axb / nunit
                except FloatingPointError:
                    msg = ''
                    for i, nid in enumerate(nodes):
                        msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                    msg += 'a x b = %s\n' % axb
                    msg += 'nunit = %s\n' % nunit
                    raise FloatingPointError(msg)
                r = centroid - p
                f = load.p * A * n * scale
                m = cross(r, f)

                node_scale = nodesi / float(nnodes)
                F += f * node_scale
                M += m * node_scale

            elif load.type == 'PLOAD1':
                #elem = self.elements[load.eid]
                elem = load.eid
                if elem.eid not in eids:
                    continue

                p1 = load.p1 * scale
                p2 = load.p2 * scale
                if elem.type not in ['CBAR', 'CBEAM', 'CBEND']:
                    raise RuntimeError('element.type=%r is not a CBAR, CBEAM, or CBEND' % elem.type)

                nodes = elem.nodeIDs()
                n1, n2 = xyz[nodes[0]], xyz[nodes[1]]
                n1 += elem.wa
                n2 += elem.wb

                deltaL = n2 - n1
                L = norm(deltaL)
                Ldir = deltaL / L
                if load.scale == 'FR':  # x1, x2 are fractional lengths
                    x1 = load.x1
                    x2 = load.x2
                    compute_fx = False
                elif load.scale == 'LE': # x1, x2 are actual lengths
                    x1 = load.x1 / L
                    x2 = load.x2 / L
                elif load.scale == 'LEPR':
                    print('LEPR continue')
                    continue
                    #msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
                    #raise NotImplementedError(msg)
                elif load.scale == 'FRPR':
                    print('FRPR continue')
                    continue
                    #msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
                    #raise NotImplementedError(msg)
                else:
                    msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
                    raise NotImplementedError(msg)

                if x1 != x2:
                    print('x1 != x2 continue')
                    continue

                #print(load)
                v = elem.get_orientation_vector(xyz)
                i = Ldir
                ki = cross(i, v)
                k = ki / norm(ki)
                j = cross(k, i)

                if load.Type in ['FX', 'FY', 'FZ']:
                    #deltaL = n2 - n1
                    r = (1 - x1) * n1 + x1 * n2
                    #print('    r =', r)
                    #print('    n1 =', n1)
                    #print('    n2 =', n2)
                    #print('    x1 =', x1)
                    #print('    1-x1 =', 1-x1)
                    #print('    deltaL =', deltaL)
                    if load.Type == 'FX':
                        if x1 == x2:
                            Fdir = array([1., 0., 0.])
                    elif load.Type == 'FY':
                        if x1 == x2:
                            Fdir = array([0., 1., 0.])
                    elif load.Type == 'FZ':
                        if x1 == x2:
                            Fdir = array([0., 0., 1.])
                    F += p1 * Fdir
                    M += cross(r - p, F)
                elif load.Type in ['MX', 'MY', 'MZ']:
                    if load.Type == 'MX':
                        if x1 == x2:
                            Mdir = array([1., 0., 0.])
                    elif load.Type == 'MY':
                        if x1 == x2:
                            Mdir = array([0., 1., 0.])
                    elif load.Type == 'MZ':
                        if x1 == x2:
                            Mdir = array([0., 0., 1.])
                    M += p1 * Mdir
                elif load.Type in ['FXE', 'FYE', 'FZE']:
                    r = (1 - x1) * n1 + x1 * n2
                    #print('\n    r =', r)
                    #print('    n1 =', n1)
                    #print('    n2 =', n2)
                    #print('    x1 =', x1)
                    #print('    1-x1 =', 1-x1)
                    #print('    i    =', i)
                    #print('    j    =', j)
                    #print('    k    =', k)
                    if load.Type == 'FXE':
                        if x1 == x2:
                            Fdir = i
                    elif load.Type == 'FYE':
                        if x1 == x2:
                            Fdir = j
                    elif load.Type == 'FZE':
                        if x1 == x2:
                            Fdir = k
                    #print('    Fdir =', Fdir, load.Type)
                    try:
                        F += p1 * Fdir
                    except FloatingPointError:
                        msg = 'eid = %s\n' % elem.eid
                        msg += 'i = %s\n' % Ldir
                        msg += 'Fdir = %s\n' % Fdir
                        msg += 'load = \n%s' % str(load)
                        raise FloatingPointError(msg)
                    M += cross(r - p, F)
                    del Fdir

                elif load.Type in ['MXE', 'MYE', 'MZE']:
                    if load.Type == 'MXE':
                        if x1 == x2:
                            Mdir = i
                    elif load.Type == 'MYE':
                        if x1 == x2:
                            Mdir = j
                    elif load.Type == 'MZE':
                        if x1 == x2:
                            Mdir = k
                    try:
                        M += p1 * Mdir
                    except FloatingPointError:
                        msg = 'eid = %s\n' % elem.eid
                        msg += 'Mdir = %s\n' % Mdir
                        msg += 'load = \n%s' % str(load)
                        raise FloatingPointError(msg)
                    del Mdir
                else:
                    raise NotImplementedError('Type=%r is not supported.  Use "FX", "FXE".' % load.Type)

            elif load.type == 'PLOAD2':
                pressure = load.pressures[0] * scale  # there are 4 pressures, but we assume p0
                for eid in load.eids:
                    if eid not in eids:
                        continue
                    elem = self.elements[eid]
                    if elem.type in ['CTRIA3',
                                     'CQUAD4', 'CSHEAR']:
                        n = elem.Normal()
                        area = elem.Area()
                        f = pressure * n * area
                        r = elem.Centroid() - p
                        m = cross(r, f)
                        F += f
                        M += m
                    else:
                        self.log.debug('case=%s etype=%r loadtype=%r not supported' % (loadcase_id, elem.type, load.type))
            elif load.type == 'PLOAD4':
                #elem = load.eid
                pressure = load.pressures[0] * scale  # there are 4 possible pressures, but we assume p0
                assert load.Cid() == 0, 'Cid() = %s' % (load.Cid())
                assert load.sorl == 'SURF', 'sorl = %s' % (load.sorl)
                assert load.ldir == 'NORM', 'ldir = %s' % (load.ldir)
                for elem in load.eids:
                    eid = elem.eid
                    if eid not in eids:
                        continue
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',]:
                        # triangles
                        nodes = elem.nodeIDs()
                        n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                        axb = cross(n1 - n2, n1 - n3)
                        nunit = norm(axb)
                        area = 0.5 * nunit
                        try:
                            n = axb / nunit
                        except FloatingPointError:
                            msg = ''
                            for i, nid in enumerate(nodes):
                                msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                            msg += 'a x b = %s\n' % axb
                            msg += 'nunit = %s\n' % nunit
                            raise FloatingPointError(msg)
                        centroid = (n1 + n2 + n3) / 3.
                    elif elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        # quads
                        nodes = elem.nodeIDs()
                        n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                        axb = cross(n1 - n3, n2 - n4)
                        nunit = norm(axb)
                        area = 0.5 * nunit
                        try:
                            n = axb / nunit
                        except FloatingPointError:
                            msg = ''
                            for i, nid in enumerate(nodes):
                                msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                            msg += 'a x b = %s\n' % axb
                            msg += 'nunit = %s\n' % nunit
                            raise FloatingPointError(msg)

                        centroid = (n1 + n2 + n3 + n4) / 4.
                    elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
                        A, centroid, normal = elem.getFaceAreaCentroidNormal(load.g34.nid, load.g1.nid)
                    else:
                        self.log.debug('case=%s eid=%s etype=%r loadtype=%r not supported' % (loadcase_id, eid, elem.type, load.type))
                        continue
                    r = centroid - p
                    f = pressure * area * n
                    #load.cid.transformToGlobal()
                    m = cross(r, f)
                    F += f
                    M += m
            elif load.type == 'GRAV':
                if include_grav:  # this will be super slow
                    g = load.GravityVector() * scale
                    for eid, elem in iteritems(self.elements):
                        if eid not in eids:
                            continue
                        centroid = elem.Centroid()
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
            self.log.debug('case=%s loadtype=%r not supported' % (loadcase_id, Type))
        #self.log.info("case=%s F=%s M=%s\n" % (loadcase_id, F, M))
        return (F, M)

    def sum_forces_moments(self, p0, loadcase_id, include_grav=False):
        """
        Sums applied forces & moments about a reference point p0 for all
        load cases.
        Considers:
          - FORCE, FORCE1, FORCE2
          - MOMENT, MOMENT1, MOMENT2
          - PLOAD, PLOAD2, PLOAD4
          - LOAD

        :param p0:
          the reference point
        :type p0:
          NUMPY.NDARRAY shape=(3,) or integer (node ID)

        :param loadcase_id:
          the LOAD=ID to analyze
        :type loadcase_id:
          integer

        :param include_grav:
          includes gravity in the summation (not supported)
        :type include_grav:
          bool

        :returns Forces:
          the forces
        :type Forces:
          NUMPY.NDARRAY shape=(3,)

        :returns Moments:
          the moments
        :type Moments:
          NUMPY.NDARRAY shape=(3,)


        .. warning:: not full validated
        .. todo:: It's super slow for cid != 0.   We can speed this up a lot
                 if we calculate the normal, area, centroid based on
                 precomputed node locations.

        Pressure acts in the normal direction per model/real/loads.bdf and loads.f06
        """
        if not isinstance(loadcase_id, integer_types):
            raise RuntimeError('loadcase_id must be an integer; loadcase_id=%r' % loadcase_id)
        if isinstance(p0, integer_types):
            p = self.model.nodes[p0].get_position()
        else:
            p = array(p0)

        load_case = self.loads[loadcase_id]
        #for (key, load_case) in iteritems(self.loads):
            #if key != loadcase_id:
                #continue

        scale_factors2 = []
        loads2 = []
        for load in load_case:
            if isinstance(load, LOAD):
                scale_factors, loads = load.get_reduced_loads()
                scale_factors2 += scale_factors
                loads2 += loads
            else:
                scale_factors2.append(1.)
                loads2.append(load)

        F = array([0., 0., 0.])
        M = array([0., 0., 0.])

        xyz = {}
        for nid, node in iteritems(self.nodes):
            xyz[nid] = node.get_position()

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
            elif load.type == 'PLOAD':
                nodes = load.nodeIDs()
                nnodes = len(nodes)
                if nnodes == 3:
                    n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                    axb = cross(n1 - n2, n1 - n3)
                    centroid = (n1 + n2 + n3) / 3.
                elif nnodes == 4:
                    n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                    axb = cross(n1 - n3, n2 - n4)
                    centroid = (n1 + n2 + n3 + n4) / 4.
                else:
                    msg = 'invalid number of nodes on PLOAD card; nodes=%s' % str(nodes)
                    raise RuntimeError(msg)

                nunit = norm(axb)
                area = 0.5 * nunit
                try:
                    n = axb / nunit
                except FloatingPointError:
                    msg = ''
                    for i, nid in enumerate(nodes):
                        msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                    msg += 'a x b = %s\n' % axb
                    msg += 'nunit = %s\n' % nunit
                    raise FloatingPointError(msg)
                r = centroid - p
                f = load.p * area * n * scale
                m = cross(r, f)

                F += f
                M += m

            elif load.type == 'PLOAD1':
                elem = load.eid

                p1 = load.p1 * scale
                p2 = load.p2 * scale
                if elem.type not in ['CBAR', 'CBEAM', 'CBEND']:
                    raise RuntimeError('element.type=%r is not a CBAR, CBEAM, or CBEND' % elem.type)

                nodes = elem.nodeIDs()
                n1, n2 = xyz[nodes[0]], xyz[nodes[1]]
                n1 += elem.wa
                n2 += elem.wb

                deltaL = n2 - n1
                L = norm(deltaL)
                Ldir = deltaL / L
                if load.scale == 'FR':  # x1, x2 are fractional lengths
                    x1 = load.x1
                    x2 = load.x2
                    #compute_fx = False
                elif load.scale == 'LE': # x1, x2 are actual lengths
                    x1 = load.x1 / L
                    x2 = load.x2 / L
                elif load.scale == 'LEPR':
                    print('LEPR continue')
                    continue
                    #msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
                    #raise NotImplementedError(msg)
                elif load.scale == 'FRPR':
                    print('FRPR continue')
                    continue
                    #msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
                    #raise NotImplementedError(msg)
                else:
                    msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
                    raise NotImplementedError(msg)

                if x1 != x2:
                    print('x1 != x2 continue')
                    continue

                v = elem.get_orientation_vector(xyz)
                i = Ldir
                ki = cross(i, v)
                k = ki / norm(ki)
                j = cross(k, i)

                if load.Type in ['FX', 'FY', 'FZ']:
                    r = (1 - x1) * n1 + x1 * n2
                    if load.Type == 'FX':
                        if x1 == x2:
                            Fdir = array([1., 0., 0.])
                    elif load.Type == 'FY':
                        if x1 == x2:
                            Fdir = array([0., 1., 0.])
                    elif load.Type == 'FZ':
                        if x1 == x2:
                            Fdir = array([0., 0., 1.])
                    F += p1 * Fdir
                    M += cross(r - p, F)
                elif load.Type in ['MX', 'MY', 'MZ']:
                    if load.Type == 'MX':
                        if x1 == x2:
                            Mdir = array([1., 0., 0.])
                    elif load.Type == 'MY':
                        if x1 == x2:
                            Mdir = array([0., 1., 0.])
                    elif load.Type == 'MZ':
                        if x1 == x2:
                            Mdir = array([0., 0., 1.])
                    M += p1 * Mdir
                elif load.Type in ['FXE', 'FYE', 'FZE']:
                    r = (1 - x1) * n1 + x1 * n2
                    if load.Type == 'FXE':
                        if x1 == x2:
                            Fdir = i
                    elif load.Type == 'FYE':
                        if x1 == x2:
                            Fdir = j
                    elif load.Type == 'FZE':
                        if x1 == x2:
                            Fdir = k
                    #print('    Fdir =', Fdir, load.Type)
                    try:
                        F += p1 * Fdir
                    except FloatingPointError:
                        msg = 'eid = %s\n' % elem.eid
                        msg += 'i = %s\n' % Ldir
                        msg += 'Fdir = %s\n' % Fdir
                        msg += 'load = \n%s' % str(load)
                        raise FloatingPointError(msg)
                    M += cross(r - p, F)
                    del Fdir

                elif load.Type in ['MXE', 'MYE', 'MZE']:
                    if load.Type == 'MXE':
                        if x1 == x2:
                            Mdir = i
                    elif load.Type == 'MYE':
                        if x1 == x2:
                            Mdir = j
                    elif load.Type == 'MZE':
                        if x1 == x2:
                            Mdir = k
                    try:
                        M += p1 * Mdir
                    except FloatingPointError:
                        msg = 'eid = %s\n' % elem.eid
                        msg += 'Mdir = %s\n' % Mdir
                        msg += 'load = \n%s' % str(load)
                        raise FloatingPointError(msg)
                    del Mdir
                else:
                    raise NotImplementedError('Type=%r is not supported.  Use "FX", "FXE".' % load.Type)

            elif load.type == 'PLOAD2':
                pressure = load.pressures[0] * scale  # there are 4 pressures, but we assume p0
                for eid in load.eids:
                    elem = self.elements[eid]
                    if elem.type in ['CTRIA3',
                                     'CQUAD4', 'CSHEAR']:
                        n = elem.Normal()
                        area = elem.Area()
                        f = pressure * n * area
                        r = elem.Centroid() - p
                        m = cross(r, f)
                        F += f
                        M += m
                    else:
                        self.log.debug('case=%s etype=%r loadtype=%r not supported' % (loadcase_id, elem.type, load.type))
            elif load.type == 'PLOAD4':
                #elem = load.eid
                pressure = load.pressures[0] * scale  # there are 4 possible pressures, but we assume p0
                assert load.Cid() == 0, 'Cid() = %s' % (load.Cid())
                assert load.sorl == 'SURF', 'sorl = %s' % (load.sorl)
                assert load.ldir == 'NORM', 'ldir = %s' % (load.ldir)
                for elem in load.eids:
                    eid = elem.eid
                    if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',]:
                        # triangles
                        nodes = elem.node_ids
                        n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                        axb = cross(n1 - n2, n1 - n3)
                        nunit = norm(axb)
                        area = 0.5 * nunit
                        try:
                            n = axb / nunit
                        except FloatingPointError:
                            msg = ''
                            for i, nid in enumerate(nodes):
                                msg += 'nid%i=%i node=%s\n' % (i + 1, nid, xyz[nodes[i]])
                            msg += 'a x b = %s\n' % axb
                            msg += 'nunit = %s\n' % nunit
                            raise FloatingPointError(msg)
                        centroid = (n1 + n2 + n3) / 3.
                    elif elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        # quads
                        nodes = elem.node_ids
                        n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                        axb = cross(n1 - n3, n2 - n4)
                        nunit = norm(axb)
                        area = 0.5 * nunit
                        try:
                            n = axb / nunit
                        except FloatingPointError:
                            msg = ''
                            for i, nid in enumerate(nodes):
                                msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                            msg += 'a x b = %s\n' % axb
                            msg += 'nunit = %s\n' % nunit
                            raise FloatingPointError(msg)

                        centroid = (n1 + n2 + n3 + n4) / 4.
                    elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
                        area, centroid, normal = elem.getFaceAreaCentroidNormal(load.g34.nid, load.g1.nid)
                    else:
                        msg = ('case=%s eid=%s etype=%r loadtype=%r not supported'
                               % (loadcase_id, eid, elem.type, load.type))
                        self.log.debug(msg)
                        continue
                    r = centroid - p
                    f = pressure * area * n
                    #load.cid.transformToGlobal()
                    m = cross(r, f)
                    F += f
                    M += m
            elif load.type == 'GRAV':
                if include_grav:  # this will be super slow
                    g = load.GravityVector() * scale
                    for eid, elem in iteritems(self.elements):
                        centroid = elem.Centroid()
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
            self.log.debug('case=%s loadtype=%r not supported' % (loadcase_id, Type))
        return (F, M)

    def MassProperties(self):
        """
        .. seealso:: mass_properties
        """
        self.deprecated('MassProperties()', 'mass_properties()', '0.7')
        return self.mass_properties()

    def Mass(self):
        """
        .. seealso:: mass
        """
        self.deprecated('Mass()', 'mass_properties()[0]', '0.7')
        mass, cg, I = self.mass_properties()
        return mass

    def mass(self, element_ids=None, sym_axis=None):
        """Calculates mass in the global coordinate system"""
        self.deprecated('mass()', 'mass_properties()[0]', '0.7')
        mass, cg, I = self.mass_properties(element_ids=element_ids,
                                           reference_point=None,
                                           sym_axis=sym_axis,
                                           num_cpus=1)
        return mass

    def resolveGrids(self, cid=0):
        """
        .. seealso:: resolve_grids
        """
        self.deprecated('resolveGrids(cid)', 'resolve_grids(cid)', '0.7')
        return self.resolve_grids(cid)

    def unresolveGrids(self, fem_old):
        """
        .. seealso:: unresolve_grids
        """
        self.deprecated('unresolveGrids(fem_old)', 'unresolve_grids(fem_old)', '0.7')
        return self.unresolve_grids(fem_old)