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
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from collections import defaultdict
from copy import deepcopy
from codecs import open

from six import iteritems, string_types, PY2
from six.moves import zip

import numpy as np
from numpy import array, cross, dot, allclose, mean
from numpy.linalg import norm

from pyNastran.utils import integer_types
from pyNastran.bdf.cards.loads.static_loads import LOAD
from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.utils.mathematics import integrate_positive_unit_line


def transform_inertia(mass, xyz_cg, xyz_ref, xyz_ref2, I_ref):
    """
    Transforms mass moment of inertia using parallel-axis theorem.

    Parameters
    ----------
    mass : float
        the mass
    xyz_cg : (3, ) float ndarray
        the CG location
    xyz_ref : (3, ) float ndarray
        the original reference location
    xyz_ref2 : (3, ) float ndarray
        the new reference location
    I_ref : (6, ) float ndarray
        the mass moment of inertias about the original reference point
        [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]

    Returns
    -------
    I_new : (6, ) float ndarray
        the mass moment of inertias about the new reference point
        [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
    """
    xcg, ycg, zcg = xyz_cg
    xref, yref, zref = xyz_ref
    xref2, yref2, zref2 = xyz_ref2

    dx1 = xcg - xref
    dy1 = ycg - yref
    dz1 = zcg - zref

    dx2 = xref2 - xcg
    dy2 = yref2 - ycg
    dz2 = zref2 - zcg
    print('dx1 = <%s, %s, %s>' % (dx1, dy1, dz1))
    print('dx2 = <%s, %s, %s>' % (dx2, dy2, dz2))

    # consistent with mass_properties, not CONM2
    print('I_ref =', I_ref)
    Ixx_ref, Iyy_ref, Izz_ref, Ixy_ref, Ixz_ref, Iyz_ref = I_ref
    Ixx2 = Ixx_ref - mass * (dx1**2 - dx2**2)
    Iyy2 = Iyy_ref - mass * (dy1**2 - dy2**2)
    Izz2 = Izz_ref - mass * (dz1**2 - dz2**2)
    Ixy2 = Ixy_ref - mass * (dx1 * dy1 - dx2 * dy2)
    Ixz2 = Ixz_ref - mass * (dx1 * dz1 - dx2 * dz2)
    Iyz2 = Iyz_ref - mass * (dy1 * dz1 - dy2 * dz2)
    I_new = np.array([Ixx2, Iyy2, Izz2, Ixy2, Ixz2, Iyz2])
    return I_new

class BDFMethods(BDFAttributes):
    """
    Has the following methods:
        mass_properties(element_ids=None, reference_point=None, sym_axis=None,
            scale=None)
        resolve_grids(cid=0)
        unresolve_grids(model_old)
        sum_forces_moments_elements(p0, loadcase_id, eids, nids,
            include_grav=False, xyz_cid0=None)
        sum_forces_moments(p0, loadcase_id, include_grav=False,
            xyz_cid0=None)
    """

    def __init__(self):
        BDFAttributes.__init__(self)

    def get_area_breakdown(self, property_ids=None, sum_bar_area=True):
        """
        gets a breakdown of the area by property region

        TODO: What about CONRODs?
        #'PBRSECT',
        #'PBCOMP',
        #'PBMSECT',
        #'PBEAM3',
        #'PBEND',
        #'PIHEX',
        #'PCOMPS',
        """
        skip_props = [
            'PSOLID', 'PLPLANE', 'PPLANE', 'PELAS',
            'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
            'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
            'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX', 'PLSOLID',
            'PCOMPS', 'PVISC', 'PBCOMP', 'PBEND',
        ]

        pid_eids = self.get_element_ids_dict_with_pids(property_ids)
        pids_to_area = {}
        for pid, eids in iteritems(pid_eids):
            prop = self.properties[pid]
            areas = []
            if prop.type in ['PSHELL', 'PCOMP', 'PSHEAR', 'PCOMPG', ]:
                for eid in eids:
                    elem = self.elements[eid]
                    try:
                        areas.append(elem.Area())
                    except AttributeError:
                        print(prop)
                        print(elem)
                        raise
            elif prop.type in ['PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE']:
                for eid in eids:
                    elem = self.elements[eid]
                    try:
                        if sum_bar_area:
                            areas.append(elem.Area())
                        else:
                            areas = [elem.Area()]
                    except AttributeError:
                        print(prop)
                        print(elem)
                        raise
            elif prop.type in skip_props:
                pass
            else:
                raise NotImplementedError(prop)
            if areas:
                pids_to_area[pid] = sum(areas)
        return pids_to_area

    def get_volume_breakdown(self, property_ids=None):
        """
        gets a breakdown of the volume by property region

        TODO: What about CONRODs?
        #'PBRSECT',
        #'PBCOMP',
        #'PBMSECT',
        #'PBEAM3',
        #'PBEND',
        #'PIHEX',
        """
        pid_eids = self.get_element_ids_dict_with_pids(property_ids)

        pids_to_volume = {}
        for pid, eids in iteritems(pid_eids):
            prop = self.properties[pid]
            volumes = []
            if prop.type in ['PSHELL', 'PSHEAR']:
                t = prop.t
                areas = []
                for eid in eids:
                    elem = self.elements[eid]
                    areas.append(elem.Area())
                volumesi = [area * t for area in areas]
                volumes.extend(volumesi)
            elif prop.type in ['PCOMP', 'PCOMPG',]:
                areas = []
                for eid in eids:
                    elem = self.elements[eid]
                    areas.append(elem.Area())
                t = prop.Thickness()
                volumesi = [area * t for area in areas]
                volumes.extend(volumesi)
            elif prop.type in ['PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE']:
                # what should I do here?
                lengths = []
                for eid in eids:
                    elem = self.elements[eid]
                    length = elem.Length()
                    lengths.append(length)
                area = prop.Area()
                volumesi = [area * length for length in lengths]
                volumes.extend(volumesi)
            elif prop.type in ['PSOLID', 'PCOMPS', 'PLSOLID']:
                for eid in eids:
                    elem = self.elements[eid]
                    volumes.append(elem.Volume())
            elif prop.type in ['PLPLANE', 'PPLANE', 'PELAS',
                               'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
                               'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
                               'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX',
                               'PVISC', 'PBCOMP', 'PBEND']:
                pass
            else:
                raise NotImplementedError(prop)
            if volumes:
                pids_to_volume[pid] = sum(volumes)
        return pids_to_volume

    def get_mass_breakdown(self, property_ids=None):
        """
        gets a breakdown of the mass by property region

        TODO: What about CONRODs, CONM2s?
        #'PBRSECT',
        #'PBCOMP',
        #'PBMSECT',
        #'PBEAM3',
        #'PBEND',
        #'PIHEX',
        #'PCOMPS',
        """
        pid_eids = self.get_element_ids_dict_with_pids(property_ids)

        mass_type_to_mass = {}
        pids_to_mass = {}
        for eid, elem in iteritems(self.masses):
            if elem.type not in mass_type_to_mass:
                mass_type_to_mass[elem.type] = elem.Mass()
            else:
                mass_type_to_mass[elem.type] += elem.Mass()

        for pid, eids in iteritems(pid_eids):
            prop = self.properties[pid]
            masses = []
            if prop.type in ['PSHELL', 'PSHEAR']:
                t = prop.t
                nsm = prop.nsm
                rho = prop.Rho()
                for eid in eids:
                    elem = self.elements[eid]
                    area = elem.Area()
                    masses.append(area * (rho * t + nsm))
            elif prop.type in ['PCOMP', 'PCOMPG']:
                for eid in eids:
                    elem = self.elements[eid]
                    masses.append(elem.Mass())
            elif prop.type in ['PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE']:
                # what should I do here?
                nsm = prop.nsm
                rho = prop.Rho()
                for eid in eids:
                    elem = self.elements[eid]
                    area = prop.Area()
                    length = elem.Length()
                    masses.append(area * (rho * length + nsm))
            elif prop.type in ['PSOLID', 'PCOMPS', 'PLSOLID']:
                rho = prop.Rho()
                for eid in eids:
                    elem = self.elements[eid]
                    masses.append(rho * elem.Volume())
            elif prop.type in ['PLPLANE', 'PPLANE', 'PELAS',
                               'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
                               'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
                               'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX',
                               'PVISC', 'PBCOMP', 'PBEND']:
                pass
            else:
                raise NotImplementedError(prop)
            if masses:
                pids_to_mass[pid] = sum(masses)
        return pids_to_mass, mass_type_to_mass

    def mass_properties(self, element_ids=None, mass_ids=None, reference_point=None,
                        sym_axis=None, scale=None):
        """
        Caclulates mass properties in the global system about the
        reference point.

        Parameters
        ----------
        element_ids : list[int]; (n, ) ndarray, optional
            An array of element ids.
        mass_ids : list[int]; (n, ) ndarray, optional
            An array of mass ids.
        reference_point : ndarray/str/int, optional
            type : ndarray
                An array that defines the origin of the frame.
                default = <0,0,0>.
            type : str
                'cg' is the only allowed string
            type : int
                the node id
        sym_axis : str, optional
            The axis to which the model is symmetric. If AERO cards are used, this can be left blank
            allowed_values = 'no', x', 'y', 'z', 'xy', 'yz', 'xz', 'xyz'
        scale : float, optional
            The WTMASS scaling value.
            default=None -> PARAM, WTMASS is used
            float > 0.0

        Returns
        -------
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

        Example 1
        ---------
        # mass properties of entire structure
        mass, cg, I = model.mass_properties()
        Ixx, Iyy, Izz, Ixy, Ixz, Iyz = I


        Example 2
        ---------
        # mass properties of model based on Property ID
        pids = list(model.pids.keys())
        pid_eids = self.get_element_ids_dict_with_pids(pids)

        for pid, eids in sorted(iteritems(pid_eids)):
            mass, cg, I = model.mass_properties(element_ids=eids)
        """
        if reference_point is None:
            reference_point = array([0., 0., 0.])
        elif isinstance(reference_point, integer_types):
            reference_point = self.nodes[reference_point].get_position()

        elements, masses = self._mass_properties_elements_init(element_ids, mass_ids)
        mass, cg, I = self._mass_properties(elements, masses,
                                            reference_point=reference_point)
        mass, cg, I = self._apply_mass_symmetry(sym_axis, scale, mass, cg, I)
        return (mass, cg, I)

    def _mass_properties_elements_init(self, element_ids, mass_ids):
        """helper method"""
        # if neither element_id nor mass_ids are specified, use everything
        if element_ids is None and mass_ids is None:
            elements = self.elements.values()
            masses = self.masses.values()

        # if either element_id or mass_ids are specified and the other is not, use only the
        # specified ids
        else:
            if element_ids is None:
                elements = []
            else:
                elements = [element for eid, element in self.elements.items() if eid in element_ids]
            if mass_ids is None:
                masses = []
            else:
                masses = [mass for eid, mass in self.masses.items() if eid in mass_ids]
        return elements, masses

    def mass_properties_no_xref(self, element_ids=None, mass_ids=None, reference_point=None,
                                sym_axis=None, scale=None):
        """
        Caclulates mass properties in the global system about the
        reference point.

        Parameters
        ----------
        element_ids : list[int]; (n, ) ndarray, optional
            An array of element ids.
        mass_ids : list[int]; (n, ) ndarray, optional
            An array of mass ids.
        reference_point : ndarray/str/int, optional
            type : ndarray
                An array that defines the origin of the frame.
                default = <0,0,0>.
            type : str
                'cg' is the only allowed string
            type : int
                the node id
        sym_axis : str, optional
            The axis to which the model is symmetric. If AERO cards are used, this can be left blank
            allowed_values = 'no', x', 'y', 'z', 'xy', 'yz', 'xz', 'xyz'
        scale : float, optional
            The WTMASS scaling value.
            default=None -> PARAM, WTMASS is used
            float > 0.0

        Returns
        -------
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

        Example 1
        ---------
        # mass properties of entire structure
        mass, cg, I = model.mass_properties()
        Ixx, Iyy, Izz, Ixy, Ixz, Iyz = I


        Example 2
        ---------
        # mass properties of model based on Property ID
        pids = list(model.pids.keys())
        pid_eids = self.get_element_ids_dict_with_pids(pids)

        for pid, eids in sorted(iteritems(pid_eids)):
            mass, cg, I = model.mass_properties(element_ids=eids)
        """
        if reference_point is None:
            reference_point = array([0., 0., 0.])
        elif isinstance(reference_point, integer_types):
            reference_point = self.nodes[reference_point].get_position()

        elements, masses = self._mass_properties_elements_init(element_ids, mass_ids)
        #nelements = len(elements) + len(masses)

        mass, cg, I = self._mass_properties_no_xref(elements, masses,
                                                    reference_point=reference_point)

        mass, cg, I = self._apply_mass_symmetry(sym_axis, scale, mass, cg, I)
        return (mass, cg, I)

    def _mass_properties_no_xref(self, elements, masses, reference_point):  # pragma: no cover
        """
        Caclulates mass properties in the global system about the
        reference point.

        Parameters
        ----------
        elements : List[int]; ndarray
            the element ids to consider
        masses : List[int]; ndarray
            the mass ids to consider
        reference_point : (3, ) ndarray; default = <0,0,0>.
            an array that defines the origin of the frame.

        Returns
        -------
        mass : float
            the mass of the model
        cg : (3, ) float NDARRAY
            the cg of the model as an array.
        I : (6, ) float NDARRAY
            moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz])

        .. seealso:: self.mass_properties
        """
        #Ixx Iyy Izz, Ixy, Ixz Iyz
        # precompute the CG location and make it the reference point
        I = array([0., 0., 0., 0., 0., 0., ])
        cg = array([0., 0., 0.])
        if isinstance(reference_point, string_types):
            if reference_point == 'cg':
                mass = 0.
                for pack in [elements, masses]:
                    for element in pack:
                        try:
                            p = element.Centroid_no_xref(self)
                            m = element.Mass_no_xref(self)
                            mass += m
                            cg += m * p
                        except:
                            #pass
                            raise
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
                    p = element.Centroid_no_xref(self)
                except:
                    #continue
                    raise

                try:
                    m = element.Mass_no_xref(self)
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
                    # PLPLANE
                    pid_ref = self.Property(element.pid)
                    if pid_ref.type == 'PSHELL':
                        self.log.warning('p=%s reference_point=%s type(reference_point)=%s' % (
                            p, reference_point, type(reference_point)))
                        raise
                    self.log.warning("could not get the inertia for element/property\n%s%s" % (
                        element, element.pid_ref))
                    continue

        if mass:
            cg /= mass
        return (mass, cg, I)

    def _mass_properties(self, elements, masses, reference_point):  # pragma: no cover
        """
        Caclulates mass properties in the global system about the
        reference point.

        Parameters
        ----------
        elements : List[int]; ndarray
            the element ids to consider
        masses : List[int]; ndarray
            the mass ids to consider
        reference_point : (3, ) ndarray; default = <0,0,0>.
            an array that defines the origin of the frame.

        Returns
        -------
        mass : float
            the mass of the model
        cg : (3, ) float NDARRAY
            the cg of the model as an array.
        I : (6, ) float NDARRAY
            moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz])

        .. seealso:: self.mass_properties
        """
        #Ixx Iyy Izz, Ixy, Ixz Iyz
        # precompute the CG location and make it the reference point
        I = array([0., 0., 0., 0., 0., 0., ])
        cg = array([0., 0., 0.])
        if isinstance(reference_point, string_types):
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
                    # PLPLANE
                    if element.pid_ref.type == 'PSHELL':
                        self.log.warning('p=%s reference_point=%s type(reference_point)=%s' % (
                            p, reference_point, type(reference_point)))
                        raise
                    self.log.warning("could not get the inertia for element/property\n%s%s" % (
                        element, element.pid_ref))

                    continue

        if mass:
            cg /= mass
        return (mass, cg, I)

    def _mass_properties_new(self, element_ids=None, mass_ids=None, reference_point=None,
                             sym_axis=None, scale=None, xyz_cid0=None):  # pragma: no cover
        """
        half implemented, not tested, should be faster someday...
        don't use this

        Caclulates mass properties in the global system about the
        reference point.

        Parameters
        ----------
        element_ids : list[int]; (n, ) ndarray, optional
            An array of element ids.
        mass_ids : list[int]; (n, ) ndarray, optional
            An array of mass ids.
        reference_point : ndarray/str/int, optional
            type : ndarray
                An array that defines the origin of the frame.
                default = <0,0,0>.
            type : str
                'cg' is the only allowed string
            type : int
                the node id
        sym_axis : str, optional
            The axis to which the model is symmetric. If AERO cards are used, this can be left blank
            allowed_values = 'no', x', 'y', 'z', 'xy', 'yz', 'xz', 'xyz'
        scale : float, optional
            The WTMASS scaling value.
            default=None -> PARAM, WTMASS is used
            float > 0.0
        xyz_cid0 : dict[nid] : xyz; default=None -> auto-calculate
            mapping of the node id to the global position

        Returns
        -------
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

        Example 1
        ---------
        # mass properties of entire structure
        mass, cg, I = model.mass_properties()
        Ixx, Iyy, Izz, Ixy, Ixz, Iyz = I


        Example 2
        ---------
        # mass properties of model based on Property ID
        pids = list(model.pids.keys())
        pid_eids = self.get_element_ids_dict_with_pids(pids)

        for pid, eids in sorted(iteritems(pid_eids)):
            mass, cg, I = model.mass_properties(element_ids=eids)
        """
        #if reference_point is None:
        reference_point = array([0., 0., 0.])

        if xyz_cid0 is None:
            xyz = {}
            for nid, node in iteritems(self.nodes):
                xyz[nid] = node.get_position()
        else:
            xyz = xyz_cid0

        elements, masses = self._mass_properties_elements_init(element_ids, mass_ids)

        #mass = 0.
        #cg = array([0., 0., 0.])
        #I = array([0., 0., 0., 0., 0., 0., ])
        #if isinstance(reference_point, string_types):
            #if reference_point == 'cg':
                #mass = 0.
                #for pack in [elements, masses]:
                    #for element in pack:
                        #try:
                            #p = element.Centroid()
                            #m = element.Mass()
                            #mass += m
                            #cg += m * p
                        #except:
                            #pass
                #if mass == 0.0:
                    #return mass, cg, I

                #reference_point = cg / mass
            #else:
                ## reference_point = [0.,0.,0.] or user-defined array
                #pass

        mass = 0.
        cg = array([0., 0., 0.])
        I = array([0., 0., 0., 0., 0., 0., ])

        no_mass = [
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', #'CLEAS5',
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
            'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CVISC', 'CGAP', # is this right?
            'CFAST',
            'CRAC2D', 'CRAC3D',

            'CSSCHD', 'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',
            'CBARAO', 'CORD1R', 'CORD2R', 'CORD1C', 'CORD2C', 'CORD1S', 'CORD2S',
            'CORD3G', 'CONV', 'CONVM', 'CSET', 'CSET1', 'CLOAD',
            'CHBDYG', 'CHBDYE', 'CHBDYP',
        ]
        all_eids = np.array(list(self.elements.keys()), dtype='int32')
        all_eids.sort()

        all_mass_ids = np.array(list(self.masses.keys()), dtype='int32')
        all_mass_ids.sort()

        #def _increment_inertia0(centroid, reference_point, m, mass, cg, I):
            #"""helper method"""
            #(x, y, z) = centroid - reference_point
            #mass += m
            #cg += m * centroid
            #return mass

        def _increment_inertia(centroid, reference_point, m, mass, cg, I):
            """helper method"""
            (x, y, z) = centroid - reference_point
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
            cg += m * centroid
            return mass

        def get_sub_eids(all_eids, eids):
            """supports limiting the element/mass ids"""
            eids = np.array(eids)
            ieids = np.searchsorted(all_eids, eids)
            eids2 = eids[all_eids[ieids] == eids]
            return eids2

        etypes_skipped = set([])
        for etype, eids in iteritems(self._type_to_id_map):
            if etype in ['CROD', 'CONROD']:
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    n1, n2 = elem.node_ids
                    length = norm(xyz[n2] - xyz[n1])
                    centroid = (xyz[n1] + xyz[n2]) / 2.
                    mpl = elem.MassPerLength()
                    m = mpl * length
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
            elif etype == 'CTUBE':
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    n1, n2 = elem.node_ids
                    length = norm(xyz[n2] - xyz[n1])
                    centroid = (xyz[n1] + xyz[n2]) / 2.
                    mpl = elem.pid_ref.MassPerLength()
                    m = mpl * length
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
            elif etype == 'CBAR':
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    n1, n2 = elem.node_ids
                    centroid = (xyz[n1] + xyz[n2]) / 2.
                    length = norm(xyz[n2] - xyz[n1])
                    mpl = elem.pid_ref.MassPerLength()
                    m = mpl * length
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
            elif etype == 'CBEAM':
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    prop = elem.pid_ref
                    n1, n2 = elem.node_ids
                    node1 = xyz[n1]
                    node2 = xyz[n2]
                    centroid = (node1 + node2) / 2.
                    length = norm(node2 - node1)
                    #cda = self.nodes[n1].cid_ref
                    #cdb = self.nodes[n2].cid_ref

                    is_failed, wa, wb, _ihat, jhat, khat = elem.get_axes(self)
                    p1 = node1 + wa
                    p2 = node2 + wb
                    if prop.type == 'PBEAM':
                        rho = prop.Rho()
                        mass_per_lengths = []
                        nsm_per_lengths = []
                        for (area, nsm) in zip(prop.A, prop.nsm):
                            mass_per_lengths.append(area * rho)
                            nsm_per_lengths.append(nsm)
                        mass_per_length = integrate_positive_unit_line(prop.xxb, mass_per_lengths)
                        nsm_per_length = integrate_positive_unit_line(prop.xxb, nsm_per_lengths)
                        #nsm = np.mean(prop.nsm)
                        nsm = nsm_per_length * length
                        m = mass_per_length * length
                        nsm_n1 = (p1 + jhat * prop.m1a + khat * prop.m2a)
                        nsm_n2 = (p2 + jhat * prop.m1b + khat * prop.m2b)
                        nsm_centroid = (nsm_n1 + nsm_n2) / 2.
                        #if nsm != 0.:
                            #p1_nsm = p1 + prop.ma
                            #p2_nsm = p2 + prop.mb
                    elif prop.type == 'PBEAML':
                        #mass_per_lengths = self.get_mass_per_lengths()
                        mass_per_length = prop.MassPerLength() # includes simplified nsm
                        m = mass_per_length * length
                        nsm_centroid = np.zeros(3)
                        nsm = prop.nsm[0]  # TODO: simplified
                    elif prop.type == 'PBCOMP':
                        mass_per_length = prop.MassPerLength()
                        m = mass_per_length * length
                        nsm = prop.nsm
                        nsm_n1 = (p1 + jhat * prop.m1 + khat * prop.m2)
                        nsm_n2 = (p2 + jhat * prop.m1 + khat * prop.m2)
                        nsm_centroid = (nsm_n1 + nsm_n2) / 2.
                    else:
                        raise NotImplementedError(prop.type)

                    #mpl = elem.pid_ref.MassPerLength()
                    #m = mpl * length

                    (x, y, z) = centroid - reference_point
                    (xm, ym, zm) = nsm_centroid - reference_point
                    x2 = x * x
                    y2 = y * y
                    z2 = z * z
                    xm2 = xm * xm
                    ym2 = ym * ym
                    zm2 = zm * zm

                    # Ixx, Iyy, Izz, Ixy, Ixz, Iyz
                    I[0] += m * (y2 + z2) + nsm * (ym2 + zm2)
                    I[1] += m * (x2 + z2) + nsm * (xm2 + zm2)
                    I[2] += m * (x2 + y2) + nsm * (xm2 + ym2)
                    I[3] += m * x * y + nsm * xm * ym
                    I[4] += m * x * z + nsm * xm * zm
                    I[5] += m * y * z + nsm * ym * zm
                    mass += m + nsm
                    cg += m * centroid + nsm * nsm_centroid

            elif etype in ['CTRIA3', 'CTRIA6', 'CTRIAR']:
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    n1, n2, n3 = elem.node_ids[:3]
                    prop = elem.pid_ref
                    centroid = (xyz[n1] + xyz[n2] + xyz[n3]) / 3.
                    area = 0.5 * norm(cross(xyz[n1] - xyz[n2], xyz[n1] - xyz[n3]))
                    if prop.type == 'PSHELL':
                        tflag = elem.tflag
                        ti = prop.Thickness()
                        if tflag == 0:
                            # absolute
                            t1 = elem.T1 if elem.T1 else ti
                            t2 = elem.T2 if elem.T2 else ti
                            t3 = elem.T3 if elem.T3 else ti
                        elif tflag == 1:
                            # relative
                            t1 = elem.T1 * ti if elem.T1 else ti
                            t2 = elem.T2 * ti if elem.T2 else ti
                            t3 = elem.T3 * ti if elem.T3 else ti
                        else:
                            raise RuntimeError('tflag=%r' % tflag)
                        assert t1 + t2 + t3 > 0., 't1=%s t2=%s t3=%s' % (t1, t2, t3)
                        t = (t1 + t2 + t3) / 3.

                        # m/A = rho * A * t + nsm
                        #mass_per_area = self.nsm + rho * self.t

                        mpa = prop.nsm + prop.Rho() * t
                        #mpa = elem.pid_ref.MassPerArea()
                        m = mpa * area
                    elif prop.type in ['PCOMP', 'PCOMPG']:
                        # PCOMP, PCOMPG
                        #rho_t = prop.get_rho_t()
                        #nsm = prop.nsm
                        #rho_t = [mat.Rho() * t for (mat, t) in zip(prop.mids_ref, prop.ts)]
                        #mpa = sum(rho_t) + nsm

                        # works for PCOMP
                        # F:\Program Files\Siemens\NXNastran\nxn10p1\nxn10p1\nast\tpl\cqr3compbuck.dat
                        mpa = prop.get_mass_per_area()
                    elif prop.type == 'PLPLANE':
                        continue
                    else:
                        raise NotImplementedError(prop.type)
                    m = area * mpa
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
            elif etype in ['CQUAD4', 'CQUAD8', 'CQUADR']:
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    n1, n2, n3, n4 = elem.node_ids[:4]
                    prop = elem.pid_ref
                    centroid = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
                    area = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n4] - xyz[n2]))

                    if prop.type == 'PSHELL':
                        tflag = elem.tflag
                        ti = prop.Thickness()
                        if tflag == 0:
                            # absolute
                            t1 = elem.T1 if elem.T1 else ti
                            t2 = elem.T2 if elem.T2 else ti
                            t3 = elem.T3 if elem.T3 else ti
                            t4 = elem.T4 if elem.T4 else ti
                        elif tflag == 1:
                            # relative
                            t1 = elem.T1 * ti if elem.T1 else ti
                            t2 = elem.T2 * ti if elem.T2 else ti
                            t3 = elem.T3 * ti if elem.T3 else ti
                            t4 = elem.T4 * ti if elem.T4 else ti
                        else:
                            raise RuntimeError('tflag=%r' % tflag)
                        assert t1 + t2 + t3 + t4 > 0., 't1=%s t2=%s t3=%s t4=%s' % (t1, t2, t3, t4)
                        t = (t1 + t2 + t3 + t4) / 4.

                        # m/A = rho * A * t + nsm
                        #mass_per_area = self.nsm + rho * self.t

                        mpa = prop.nsm + prop.Rho() * t
                        #mpa = elem.pid_ref.MassPerArea()
                        #m = mpa * area
                    elif prop.type in ['PCOMP', 'PCOMPG']:
                        # PCOMP, PCOMPG
                        #rho_t = prop.get_rho_t()
                        #nsm = prop.nsm
                        #rho_t = [mat.Rho() * t for (mat, t) in zip(prop.mids_ref, prop.ts)]
                        #mpa = sum(rho_t) + nsm
                        mpa = prop.get_mass_per_area()
                    elif prop.type == 'PLPLANE':
                        continue
                        #raise NotImplementedError(prop.type)
                    else:
                        raise NotImplementedError(prop.type)
                    m = area * mpa
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
            elif etype == 'CQUAD':
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    n1, n2, n3, n4 = elem.node_ids[:4]
                    prop = elem.pid_ref
                    centroid = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
                    area = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n4] - xyz[n2]))

                    if prop.type == 'PSHELL':
                        t = prop.Thickness()
                        mpa = prop.nsm + prop.Rho() * t
                    elif prop.type in ['PCOMP', 'PCOMPG']:
                        mpa = prop.get_mass_per_area()
                    elif prop.type == 'PLPLANE':
                        continue
                        #raise NotImplementedError(prop.type)
                    else:
                        raise NotImplementedError(prop.type)
                    m = area * mpa
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

            elif etype == 'CSHEAR':
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    n1, n2, n3, n4 = elem.node_ids
                    prop = elem.pid_ref
                    centroid = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
                    area = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n4] - xyz[n2]))
                    mpa = prop.MassPerArea()
                    m = area * mpa
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
            elif etype in ['CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4']:
                eids2 = get_sub_eids(all_mass_ids, eids)
                for eid in eids2:
                    elem = self.masses[eid]
                    m = elem.Mass()
                    centroid = elem.Centroid()
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
            elif etype == 'CTETRA':
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    n1, n2, n3, n4 = elem.node_ids[:4]
                    centroid = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
                    #V = -dot(n1 - n4, cross(n2 - n4, n3 - n4)) / 6.
                    volume = -dot(xyz[n1] - xyz[n4], cross(xyz[n2] - xyz[n4], xyz[n3] - xyz[n4])) / 6.
                    m = elem.Rho() * volume
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
            elif etype == 'CPYRAM':
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    n1, n2, n3, n4, n5 = elem.node_ids[:5]
                    centroid1 = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
                    area1 = 0.5 * norm(cross(xyz[n3]-xyz[n1], xyz[n4]-xyz[n2]))
                    centroid5 = xyz[n5]
                    centroid = (centroid1 + centroid5) / 2.
                    volume = area1 / 3. * norm(centroid1 - centroid5)
                    m = elem.Rho() * volume
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
            elif etype == 'CPENTA':
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    n1, n2, n3, n4, n5, n6 = elem.node_ids[:6]
                    area1 = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n2] - xyz[n1]))
                    area2 = 0.5 * norm(cross(xyz[n6] - xyz[n4], xyz[n5] - xyz[n4]))
                    centroid1 = (xyz[n1] + xyz[n2] + xyz[n3]) / 3.
                    centroid2 = (xyz[n4] + xyz[n5] + xyz[n6]) / 3.
                    centroid = (centroid1 + centroid2) / 2.
                    volume = (area1 + area2) / 2. * norm(centroid1 - centroid2)
                    m = elem.Rho() * volume
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

            elif etype == 'CHEXA':
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    n1, n2, n3, n4, n5, n6, n7, n8 = elem.node_ids[:8]
                    #(A1, c1) = area_centroid(n1, n2, n3, n4)
                    centroid1 = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
                    area1 = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n4] - xyz[n2]))
                    #(A2, c2) = area_centroid(n5, n6, n7, n8)
                    centroid2 = (xyz[n5] + xyz[n6] + xyz[n7] + xyz[n8]) / 4.
                    area2 = 0.5 * norm(cross(xyz[n7] - xyz[n5], xyz[n8] - xyz[n6]))

                    volume = (area1 + area2) / 2. * norm(centroid1 - centroid2)
                    m = elem.Rho() * volume
                    centroid = (centroid1 + centroid2) / 2.
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

            elif etype in no_mass:
                continue
            elif etype == 'CBEND':
                self.log.info('elem.type=%s doesnt have mass' % etype)
                continue
            elif etype.startswith('C'):
                eids2 = get_sub_eids(all_eids, eids)
                for eid in eids2:
                    elem = self.elements[eid]
                    m = elem.Mass()
                    centroid = elem.Centroid()
                    if m > 0.0:
                        self.log.info('elem.type=%s is not supported in new '
                                      'mass properties method' % elem.type)
                        mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
                    elif etype not in etypes_skipped:
                        self.log.info('elem.type=%s doesnt have mass' % elem.type)
                        etypes_skipped.add(etype)

        if mass:
            cg /= mass
        mass, cg, I = self._apply_mass_symmetry(sym_axis, scale, mass, cg, I)
        # Ixx, Iyy, Izz, Ixy, Ixz, Iyz = I
        return mass, cg, I

    def _apply_mass_symmetry(self, sym_axis, scale, mass, cg, I):
        """
        Scales the mass & moement of inertia based on the symmetry axes
        and the PARAM WTMASS card
        """
        if isinstance(sym_axis, string_types):
            sym_axis = [sym_axis]
        elif isinstance(sym_axis, (list, tuple)):
            # basically overwrite the existing values on the AERO/AEROS card
            pass
        else:
            sym_axis = []

            # The symmetry flags on the AERO/AEROS must be the same, so
            # it doesn't matter which we one pick.  However, they might
            # not both be defined.
            if self.aero is not None:
                if self.aero.is_symmetric_xy():
                    sym_axis.append('xy')
                if self.aero.is_symmetric_xz():
                    sym_axis.append('xz')
                if self.aero.is_anti_symmetric_xy():
                    sym_axis.append('xy')
                    #raise NotImplementedError('%s load is anti-symmetric about the XY plane'
                                              #% str(aero))
                if self.aero.is_anti_symmetric_xz():
                    #raise NotImplementedError('%s load is anti-symmetric about the XZ plane'
                                              #% str(aero))
                    sym_axis.append('xz')

            if self.aeros is not None:
                if self.aeros.is_symmetric_xy():
                    sym_axis.append('xy')
                if self.aeros.is_symmetric_xz():
                    sym_axis.append('xz')
                if self.aeros.is_anti_symmetric_xy():
                    sym_axis.append('xy')
                    #raise NotImplementedError('%s load is anti-symmetric about the XY plane'
                                              #% str(aeros))
                if self.aeros.is_anti_symmetric_xz():
                    sym_axis.append('xz')
                    #raise NotImplementedError('%s load is anti-symmetric about the XZ plane' %
                                              #str(aeros))

        sym_axis = list(set(sym_axis))
        short_sym_axis = [sym_axisi.lower() for sym_axisi in sym_axis]
        is_no = 'no' in short_sym_axis
        if is_no and len(short_sym_axis) > 1:
            raise RuntimeError('no can only be used by itself; sym_axis=%s' % (str(sym_axis)))
        for sym_axisi in sym_axis:
            if sym_axisi.lower() not in ['no', 'xy', 'yz', 'xz']:
                msg = 'sym_axis=%r is invalid; sym_axisi=%r; allowed=[no, xy, yz, xz]' % (
                    sym_axis, sym_axisi)
                raise RuntimeError(msg)

        if sym_axis:
            # either we figured sym_axis out from the AERO cards or the user told us
            self.log.debug('Mass/MOI sym_axis = %r' % sym_axis)

            if 'xz' in sym_axis:
                # y intertias are 0
                cg[1] = 0.0
                mass *= 2.0
                I[0] *= 2.0
                I[1] *= 2.0
                I[2] *= 2.0
                I[3] *= 0.0  # Ixy
                I[4] *= 2.0  # Ixz; no y
                I[5] *= 0.0  # Iyz

            if 'xy' in sym_axis:
                # z intertias are 0
                cg[2] = 0.0
                mass *= 2.0
                I[0] *= 2.0
                I[1] *= 2.0
                I[2] *= 2.0
                I[3] *= 2.0  # Ixy; no z
                I[4] *= 0.0  # Ixz
                I[5] *= 0.0  # Iyz

            if 'yz' in sym_axis:
                # x intertias are 0
                cg[0] = 0.0
                mass *= 2.0
                I[0] *= 2.0
                I[1] *= 2.0
                I[2] *= 2.0
                I[3] *= 0.0  # Ixy
                I[4] *= 0.0  # Ixz
                I[5] *= 2.0  # Iyz; no x

        if scale is None and 'WTMASS' in self.params:
            param = self.params['WTMASS']
            #assert isinstance(param, PARAM), 'param=%s' % param
            scale = param.values[0]
            if scale != 1.0:
                self.log.info('WTMASS scale = %r' % scale)
        elif scale is None:
            scale = 1.0
        mass *= scale
        I *= scale
        return (mass, cg, I)

    def resolve_grids(self, cid=0):
        """
        Puts all nodes in a common coordinate system (mainly for cid testing)

        Parameters
        ----------
        cid : int; default=0
            the cid to resolve the nodes to

        .. note::

           loses association with previous coordinate systems so to go
           back requires another FEM
        """
        assert cid in self.coords, ('cannot resolve nodes to '
                                    'cid=%r b/c it doesnt exist' % cid)
        for nid, node in sorted(iteritems(self.nodes)):
            p = node.get_position_wrt(self, cid)
            node.set_position(self, p, cid)

    #def __gravity_load(self, loadcase_id):
        #"""
        #.. todo::
            #1.  resolve the load case
            #2.  grab all of the GRAV cards and combine them into one
                #GRAV vector
            #3.  run mass_properties to get the mass
            #4.  multiply by the gravity vector
        #"""

        #gravity_i = self.loads[2][0]  ## .. todo:: hardcoded
        #gi = gravity_i.N * gravity_i.scale
        #p0 = array([0., 0., 0.])  ## .. todo:: hardcoded
        #mass, cg, I = self.mass_properties(reference_point=p0, sym_axis=None)

    def sum_forces_moments_elements(self, p0, loadcase_id, eids, nids,
                                    include_grav=False, xyz_cid0=None):
        """
        Sum the forces/moments based on a list of nodes and elements.

        Parameters
        ----------
        eids : List[int]
            the list of elements to include (e.g. the loads due to a PLOAD4)
        nids : List[int]
            the list of nodes to include (e.g. the loads due to a FORCE card)
        p0 : int; (3,) ndarray
           the point to sum moments about
           type = int
               sum moments about the specified grid point
           type = (3, ) ndarray/list (e.g. [10., 20., 30]):
               the x, y, z location in the global frame
        loadcase_id : int
            the LOAD=ID to analyze
        include_grav : bool; default=False
            includes gravity in the summation (not supported)
        xyz_cid0 : None / Dict[int] = (3, ) ndarray
            the nodes in the global coordinate system

        Returns
        -------
        Forces : NUMPY.NDARRAY shape=(3,)
            the forces
        Moments : NUMPY.NDARRAY shape=(3,)
            the moments

        Nodal Types  : FORCE, FORCE1, FORCE2,
                       MOMENT, MOMENT1, MOMENT2,
                       PLOAD
        Element Types: PLOAD1, PLOAD2, PLOAD4, GRAV

        If you have a CQUAD4 (eid=3) with a PLOAD4 (sid=3) and a FORCE
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

        if eids is None:
            eids = list(self.element_ids)
        if nids is None:
            nids = list(self.node_ids)

        load_case = self.loads[loadcase_id]
        #for (key, load_case) in iteritems(self.loads):
            #if key != loadcase_id:
                #continue

        scale_factors2 = []
        loads2 = []
        is_grav = False
        for load in load_case:
            if isinstance(load, LOAD):
                scale_factors, loads = load.get_reduced_loads()
                for scale_factor, loadi in zip(scale_factors, loads):
                    if scale_factor == 0.0:
                        continue
                    scale_factors2.append(scale_factor)
                    loads2.append(loadi)
                #scale_factors2 += scale_factors
                #loads2 += loads
            elif load.type in 'GRAV':
                scale_factors2.append(1.)
                loads2.append(load)
                is_grav = True
            else:
                scale_factors2.append(1.)
                loads2.append(load)

        F = array([0., 0., 0.])
        M = array([0., 0., 0.])

        if xyz_cid0 is None:
            xyz = {}
            for nid, node in iteritems(self.nodes):
                xyz[nid] = node.get_position()
        else:
            xyz = xyz_cid0

        unsupported_types = set([])
        for load, scale in zip(loads2, scale_factors2):
            #if load.type not in ['FORCE1']:
                #continue
            #print(load.type)
            if load.type == 'FORCE':
                if load.node_id not in nids:
                    continue
                if load.Cid() != 0:
                    cp = load.cid
                    #from pyNastran.bdf.bdf import CORD2R
                    #cp = CORD2R()
                    f = load.mag * cp.transform_vector_to_global(load.xyz) * scale
                else:
                    f = load.mag * load.xyz * scale

                node = self.Node(load.node_id)
                r = xyz[node.nid] - p
                m = cross(r, f)
                F += f
                M += m

            elif load.type == 'FORCE1':
                not_found_nid = False
                for nid in load.node_ids:
                    if nid not in nids:
                        not_found_nid = True
                        break
                if not_found_nid:
                    continue

                f = load.mag * load.xyz * scale
                node = self.Node(load.node_id)
                r = xyz[node.nid] - p
                m = cross(r, f)
                F += f
                M += m
            elif load.type == 'FORCE2':
                not_found_nid = False
                for nid in load.node_ids:
                    if nid not in nids:
                        not_found_nid = True
                        break
                if not_found_nid:
                    continue

                f = load.mag * load.xyz * scale
                node = self.Node(load.node_id)
                r = xyz[node.nid] - p
                m = cross(r, f)
                F += f
                M += m
            elif load.type == 'MOMENT':
                not_found_nid = False
                for nid in load.node_ids:
                    if nid not in nids:
                        not_found_nid = True
                        break
                if not_found_nid:
                    continue

                if load.Cid() != 0:
                    cp = load.cid
                    m = cp.transform_vector_to_global(load.xyz)
                else:
                    m = load.xyz
                M += load.mag * m * scale
            elif load.type == 'MOMENT1':
                not_found_nid = False
                for nid in load.node_ids:
                    if nid not in nids:
                        not_found_nid = True
                        break
                if not_found_nid:
                    continue
                m = load.mag * load.xyz * scale
                M += m
            elif load.type == 'MOMENT2':
                not_found_nid = False
                for nid in load.node_ids:
                    if nid not in nids:
                        not_found_nid = True
                        break
                if not_found_nid:
                    continue
                m = load.mag * load.xyz * scale
                M += m

            elif load.type == 'PLOAD':
                nodes = load.node_ids
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
                    raise RuntimeError('invalid number of nodes on PLOAD card; '
                                       'nodes=%s' % str(nodes))
                if nodes[0] in nids:
                    nodesi += 1
                if nodes[1] in nids:
                    nodesi += 1
                if nodes[2] in nids:
                    nodesi += 1

                nunit = norm(axb)
                A = 0.5 * nunit
                try:
                    normal = axb / nunit
                except FloatingPointError:
                    msg = ''
                    for i, nid in enumerate(nodes):
                        msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                    msg += 'a x b = %s\n' % axb
                    msg += 'nunit = %s\n' % nunit
                    raise FloatingPointError(msg)
                r = centroid - p
                f = load.pressure * A * normal * scale
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

                nodes = elem.node_ids
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
                    print('PLOAD1 LEPR continue')
                    continue
                    #msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
                    #raise NotImplementedError(msg)
                elif load.scale == 'FRPR':
                    print('PLOAD1 FRPR continue')
                    continue
                    #msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
                    #raise NotImplementedError(msg)
                else:
                    msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
                    raise NotImplementedError(msg)

                if x1 != x2:
                    print('PLOAD1 x1 != x2 continue\n%s%s'% (str(elem), str(load)))
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
                    raise NotImplementedError('Type=%r is not supported.  '
                                              'Use "FX", "FXE".' % load.Type)

            elif load.type == 'PLOAD2':
                pressure = load.pressure * scale
                for eid in load.element_ids:
                    if eid not in eids:
                        continue
                    elem = self.elements[eid]
                    if elem.type in ['CTRIA3', 'CQUAD4', 'CSHEAR']:
                        normal = elem.Normal()
                        area = elem.Area()
                        f = pressure * normal * area
                        r = elem.Centroid() - p
                        m = cross(r, f)
                        F += f
                        M += m
                    else:
                        self.log.debug('case=%s etype=%r loadtype=%r not supported' % (
                            loadcase_id, elem.type, load.type))
            elif load.type == 'PLOAD4':
                assert load.Cid() == 0, 'Cid() = %s' % (load.Cid())
                assert load.line_load_dir == 'NORM', 'line_load_dir = %s' % (load.line_load_dir)
                for elem in load.eids:
                    eid = elem.eid
                    if eid not in eids:
                        continue
                    etype = elem.type
                    if etype in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',]:
                        # triangles
                        nodes = elem.node_ids
                        n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                        axb = cross(n1 - n2, n1 - n3)
                        nunit = norm(axb)
                        area = 0.5 * nunit
                        try:
                            normal = axb / nunit
                        except FloatingPointError:
                            msg = ''
                            for i, nid in enumerate(nodes):
                                msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                            msg += 'a x b = %s\n' % axb
                            msg += 'nunit = %s\n' % nunit
                            raise FloatingPointError(msg)
                        centroid = (n1 + n2 + n3) / 3.
                        nface = 3
                    elif etype in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        # quads
                        nodes = elem.node_ids
                        n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                        axb = cross(n1 - n3, n2 - n4)
                        nunit = norm(axb)
                        area = 0.5 * nunit
                        try:
                            normal = axb / nunit
                        except FloatingPointError:
                            msg = ''
                            for i, nid in enumerate(nodes):
                                msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                            msg += 'a x b = %s\n' % axb
                            msg += 'nunit = %s\n' % nunit
                            raise FloatingPointError(msg)
                        centroid = (n1 + n2 + n3 + n4) / 4.
                        nface = 4

                    elif etype == 'CTETRA':
                        #face = elem.get_face(load.g1.nid, load.g34.nid)
                        face_acn = elem.get_face_area_centroid_normal(load.g1.nid, load.g34.nid)
                        face, area, centroid, normal = face_acn
                        nface = 3

                    elif etype == 'CHEXA':
                        #face = elem.get_face(load.g1.nid, load.g34.nid)
                        face_acn = elem.get_face_area_centroid_normal(load.g1.nid, load.g34.nid)
                        face, area, centroid, normal = face_acn
                        nface = 4

                    elif etype == 'CPENTA':
                        g1 = load.g1.nid
                        if load.g34 is None:
                            #face = elem.get_face(g1)
                            face_acn = elem.get_face_area_centroid_normal(g1)
                            nface = 3
                        else:
                            #face = elem.get_face(load.g1.nid, load.g34.nid)
                            face_acn = elem.get_face_area_centroid_normal(g1, load.g34.nid)
                            nface = 4
                        face, area, centroid, normal = face_acn
                    else:
                        self.log.debug('case=%s eid=%s etype=%r loadtype=%r not supported' % (
                            loadcase_id, eid, etype, load.type))
                        continue
                    r = centroid - p

                    pressures = load.pressures[:nface]
                    assert len(pressures) == nface
                    if min(pressures) != max(pressures):
                        pressure = mean(pressures)
                        msg = ('%s%s\npressure.min=%s != pressure.max=%s using average of %%s; '
                               'load=%s eid=%%s'  % (str(load), str(elem), min(pressures),
                                                     max(pressures), load.sid))
                        #print(msg % (pressure, eid))
                    else:
                        pressure = pressures[0]

                    if  load.surf_or_line == 'SURF':
                        if norm(load.nvector) != 0.0 or load.Cid() != 0:
                            normal = load.nvector / np.linalg.norm(load.nvector)
                            assert load.Cid() == 0, 'cid=%r on a PLOAD4 is not supported\n%s' % (load.Cid(), str(load))
                    else:
                        msg = 'surf_or_line=%r on PLOAD4 is not supported\n%s' % (
                            load.surf_or_line, str(load))
                        raise NotImplementedError(msg)

                    f = pressure * area * normal * scale
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

    def sum_forces_moments(self, p0, loadcase_id, include_grav=False, xyz_cid0=None):
        """
        Sums applied forces & moments about a reference point p0 for all
        load cases.
        Considers:
          - FORCE, FORCE1, FORCE2
          - MOMENT, MOMENT1, MOMENT2
          - PLOAD, PLOAD2, PLOAD4
          - LOAD

        Parameters
        ----------
        p0 : NUMPY.NDARRAY shape=(3,) or integer (node ID)
            the reference point
        loadcase_id : int
            the LOAD=ID to analyze
        include_grav : bool; default=False
            includes gravity in the summation (not supported)
        xyz_cid0 : None / Dict[int] = (3, ) ndarray
            the nodes in the global coordinate system

        Returns
        -------
        Forces : NUMPY.NDARRAY shape=(3,)
            the forces
        Moments : NUMPY.NDARRAY shape=(3,)
            the moments

        .. warning:: not full validated
        .. todo:: It's super slow for cid != 0.   We can speed this up a lot
                 if we calculate the normal, area, centroid based on
                 precomputed node locations.

        Pressure acts in the normal direction per model/real/loads.bdf and loads.f06
        """
        if not isinstance(loadcase_id, integer_types):
            raise RuntimeError('loadcase_id must be an integer; loadcase_id=%r' % loadcase_id)

        cid = 0
        if isinstance(p0, integer_types):
            if cid == 0:
                p = self.nodes[p0].get_position()
            else:
                p = self.nodes[p0].get_position_wrt(self, cid)
        else:
            p = array(p0)

        try:
            load_case = self.loads[loadcase_id]
        except:
            msg = 'load_case=%s is invalid; ' % loadcase_id
            msg += 'load_cases = %s\n' % np.unique(list(self.loads.keys()))
            for subcase_id, subcase in iteritems(self.subcases):
                if 'LOAD' in subcase:
                    load_id = subcase.get_parameter('LOAD')[0]
                    msg += '  SUBCASE %i; LOAD=%s\n' % (subcase_id, load_id)
                else:
                    msg += '  SUBCASE %i has no LOAD\n' % (subcase_id)
                print(msg)
            raise KeyError(msg)
        #for (key, load_case) in iteritems(self.loads):
            #if key != loadcase_id:
                #continue

        scale_factors2 = []
        loads2 = []
        is_grav = False
        for load in load_case:
            if isinstance(load, LOAD):
                scale_factors, loads = load.get_reduced_loads()
                for scale_factor, loadi in zip(scale_factors, loads):
                    if scale_factor == 0.0:
                        continue
                    scale_factors2.append(scale_factor)
                    loads2.append(loadi)
                #scale_factors2 += scale_factors
                #loads2 += loads
            elif load.type in 'GRAV':
                scale_factors2.append(1.)
                loads2.append(load)
                is_grav = True
            else:
                scale_factors2.append(1.)
                loads2.append(load)

        F = array([0., 0., 0.])
        M = array([0., 0., 0.])

        if xyz_cid0 is None:
            xyz = {}
            for nid, node in iteritems(self.nodes):
                xyz[nid] = node.get_position()
        else:
            xyz = xyz_cid0

        unsupported_types = set([])
        for load, scale in zip(loads2, scale_factors2):
            #if load.type not in ['FORCE1']:
                #continue
            if load.type == 'FORCE':
                if load.Cid() != 0:
                    cp = load.cid
                    #from pyNastran.bdf.bdf import CORD2R
                    #cp = CORD2R()
                    f = load.mag * cp.transform_vector_to_global(load.xyz) * scale
                else:
                    f = load.mag * load.xyz * scale

                node = self.Node(load.node_id)
                r = xyz[node.nid] - p
                m = cross(r, f)
                F += f
                M += m
            elif load.type == 'FORCE1':
                f = load.mag * load.xyz * scale
                node = self.Node(load.node_id)
                r = xyz[node.nid] - p
                m = cross(r, f)
                F += f
                M += m
            elif load.type == 'FORCE2':
                f = load.mag * load.xyz * scale
                node = self.Node(load.node_id)
                r = xyz[node.nid] - p
                m = cross(r, f)
                F += f
                M += m
            elif load.type == 'MOMENT':
                if load.Cid() != 0:
                    cp = load.cid
                    #from pyNastran.bdf.bdf import CORD2R
                    #cp = CORD2R()
                    m = load.mag * cp.transform_vector_to_global(load.xyz) * scale
                else:
                    m = load.mag * load.xyz * scale
                M += m
            elif load.type == 'MOMENT1':
                m = load.mag * load.xyz * scale
                M += m
            elif load.type == 'MOMENT2':
                m = load.mag * load.xyz * scale
                M += m

            elif load.type == 'PLOAD':
                nodes = load.node_ids
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
                    normal = axb / nunit
                except FloatingPointError:
                    msg = ''
                    for i, nid in enumerate(nodes):
                        msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                    msg += 'a x b = %s\n' % axb
                    msg += 'nunit = %s\n' % nunit
                    raise FloatingPointError(msg)
                r = centroid - p
                f = load.pressure * area * normal * scale
                m = cross(r, f)

                F += f
                M += m

            elif load.type == 'PLOAD1':
                elem = load.eid

                p1 = load.p1 * scale
                p2 = load.p2 * scale
                if elem.type in ['CBAR', 'CBEAM']:
                    nodes = elem.node_ids
                    n1, n2 = xyz[nodes[0]], xyz[nodes[1]]
                    n1 += elem.wa
                    n2 += elem.wb

                    deltaL = n2 - n1
                    L = norm(deltaL)
                    try:
                        Ldir = deltaL / L
                    except:
                        msg = 'Length=0.0; nid1=%s nid2=%s\n' % (nodes[0], nodes[1])
                        msg += '%s%s' % (str(elem.nodes[0]), str(elem.nodes[1]))
                        raise FloatingPointError(msg)
                    if load.scale == 'FR':  # x1, x2 are fractional lengths
                        x1 = load.x1
                        x2 = load.x2
                        #compute_fx = False
                    elif load.scale == 'LE': # x1, x2 are actual lengths
                        x1 = load.x1 / L
                        x2 = load.x2 / L
                    elif load.scale == 'LEPR':
                        print('PLOAD1 LEPR continue')
                        continue
                        #msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
                        #raise NotImplementedError(msg)
                    elif load.scale == 'FRPR':
                        print('PLOAD1 FRPR continue')
                        continue
                        #msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
                        #raise NotImplementedError(msg)
                    else:
                        msg = 'PLOAD1 scale=%r is not supported.  Use "FR", "LE".' % load.scale
                        raise NotImplementedError(msg)

                    # FY - force in basic coordinate system
                    # FR - fractional;
                    assert x1 <= x2, 'x1=%s x2=%s' % (x1, x2)
                    if x1 != x2:
                        # continue
                        if not load.type in ['FX', 'FY', 'FZ']:
                            print('PLOAD1 x1 != x2 continue; x1=%s x2=%s; scale=%r\n%s%s'% (
                                x1, x2, load.scale, str(elem), str(load)))
                            continue
                        print('check this...PLOAD1 x1 != x2; x1=%s x2=%s; scale=%r\n%s%s'% (
                            x1, x2, load.scale, str(elem), str(load)))

                        # y = (y2-y1)/(x2-x1)*(x-x1) + y1
                        # y = (y2-y1) * (x-x1)/(x2-x1) + y1
                        # y = y2*(x-x1)/(x2-x1) + y1*(1-(x-x1)/(x2-x1))
                        # y = y2 * r + y1 * (1-r)
                        # r = (x-x1)/(x2-x1)
                        #
                        # y = y2 * r + y1 - y1 * r
                        # yi = y2 * ri + y1 * x + y1 * ri
                        # yi = y2 * ri + y1 * (x2-x1) + y1 * ri
                        #
                        # ri = integral(r)
                        # ri = 1/(x2-x1) * (0.5) * (x1-x2)**2
                        #
                        # yi = integral(y)
                        # yi = y2 * ri + y1 * (x2-x1) + y1 * ri
                        # ri = 1./(x2-x1) * (0.5) * (x1-x2)**2
                        # y1 = p1
                        # y2 = p2
                        # yi = y2 * ri + y1 * (x2-x1) + y1 * ri
                        # F = yi
                        if allclose(p1, -p2):
                            Ftotal = p1
                            x = (x1 + x2) / 2.
                        else:
                            Ftotal = L * (x2-x1) * (p1 + p2)/2.
                            Mx = L * p1 * (x2-x1)/2. + L * (p2-p1) * (2./3. * x2 + 1./3. * x1)
                            x = Mx / Ftotal
                        print('L=%s x1=%s x2=%s p1/L=%s p2/L=%s Ftotal=%s Mtotal=%s x=%s' % (
                            L, x1, x2, p1, p2, Ftotal, Mx, x))

                        i = Ldir
                        if load.Type in ['FX', 'FY', 'FZ']:
                            r = (1. - x) * n1 + x * n2
                            # print('r=%s n1=%s n2=%s' % (r, n1, n2))
                            if load.Type == 'FX':
                                Fdir = array([1., 0., 0.])
                            elif load.Type == 'FY':
                                Fdir = array([0., 1., 0.])
                            elif load.Type == 'FZ':
                                Fdir = array([0., 0., 1.])
                            else:
                                raise NotImplementedError('Type=%r is not supported.  '
                                                          'Use "FX", "FY", "FZ".' % load.Type)

                        Fi = Ftotal * Fdir
                        Mi = cross(r - p, Fdir * Ftotal)
                        F += Fi
                        M += Mi
                        print('Fi=%s Mi=%s x=%s' % (Fi, Mi, x))
                    else:
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
                            raise NotImplementedError('Type=%r is not supported.  '
                                                      'Use "FX", "FXE".' % load.Type)
                else:
                    # CBEND
                    raise RuntimeError('element.type=%r is not a CBAR, CBEAM' % elem.type)

            elif load.type == 'PLOAD2':
                pressure = load.pressure * scale
                for eid in load.element_ids:
                    elem = self.elements[eid]
                    if elem.type in ['CTRIA3', 'CQUAD4', 'CSHEAR']:
                        n = elem.Normal()
                        area = elem.Area()
                        f = pressure * n * area
                        r = elem.Centroid() - p
                        m = cross(r, f)
                        F += f
                        M += m
                    else:
                        self.log.debug('case=%s etype=%r loadtype=%r not supported' % (
                            loadcase_id, elem.type, load.type))
            elif load.type == 'PLOAD4':
                assert load.Cid() == 0, 'Cid() = %s' % (load.Cid())
                #assert load.surf_or_line == 'SURF', 'surf_or_line = %r' % (load.surf_or_line)
                assert load.line_load_dir == 'NORM', 'line_load_dir = %s' % (load.line_load_dir)

                for elem in load.eids:
                    eid = elem.eid
                    etype = elem.type
                    if etype in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',]:
                        # triangles
                        nodes = elem.node_ids
                        n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                        axb = cross(n1 - n2, n1 - n3)
                        nunit = norm(axb)
                        area = 0.5 * nunit
                        try:
                            normal = axb / nunit
                        except FloatingPointError:
                            msg = ''
                            for i, nid in enumerate(nodes):
                                msg += 'nid%i=%i node=%s\n' % (i + 1, nid, xyz[nodes[i]])
                            msg += 'a x b = %s\n' % axb
                            msg += 'nunit = %s\n' % nunit
                            raise FloatingPointError(msg)
                        centroid = (n1 + n2 + n3) / 3.
                        nface = 3
                    elif etype in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                        # quads
                        nodes = elem.node_ids
                        n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                        axb = cross(n1 - n3, n2 - n4)
                        nunit = norm(axb)
                        area = 0.5 * nunit
                        try:
                            normal = axb / nunit
                        except FloatingPointError:
                            msg = ''
                            for i, nid in enumerate(nodes):
                                msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
                            msg += 'a x b = %s\n' % axb
                            msg += 'nunit = %s\n' % nunit
                            raise FloatingPointError(msg)

                        centroid = (n1 + n2 + n3 + n4) / 4.
                        nface = 4
                    elif etype == 'CTETRA':
                        #face1 = elem.get_face(load.g1.nid, load.g34.nid)
                        face_acn = elem.get_face_area_centroid_normal(load.g1.nid, load.g34.nid)
                        face, area, centroid, normal = face_acn
                        #assert face == face1
                        nface = 3
                    elif etype == 'CHEXA':
                        #face1 = elem.get_face(load.g34.nid, load.g1.nid)
                        face_acn = elem.get_face_area_centroid_normal(load.g34.nid, load.g1.nid)
                        face, area, centroid, normal = face_acn
                        #assert face == face1
                        nface = 4
                    elif etype == 'CPENTA':
                        g1 = load.g1.nid
                        if load.g34 is None:
                            #face1 = elem.get_face(g1)
                            face_acn = elem.get_face_area_centroid_normal(g1)
                            nface = 3
                        else:
                            #face1 = elem.get_face(g1, load.g34.nid)
                            face_acn = elem.get_face_area_centroid_normal(g1, load.g34.nid)
                            nface = 4
                        face, area, centroid, normal = face_acn
                        #assert face == face1
                    else:
                        msg = ('case=%s eid=%s etype=%r loadtype=%r not supported'
                               % (loadcase_id, eid, etype, load.type))
                        self.log.debug(msg)
                        continue

                    pressures = load.pressures[:nface]
                    assert len(pressures) == nface
                    if min(pressures) != max(pressures):
                        pressure = mean(pressures)
                        msg = ('%s%s\npressure.min=%s != pressure.max=%s using average of %%s; '
                               'load=%s eid=%%s'  % (str(load), str(elem), min(pressures),
                                                     max(pressures), load.sid))

                        #print(msg % (pressure, eid))
                    else:
                        pressure = load.pressures[0]

                    if load.surf_or_line == 'SURF':
                        if norm(load.nvector) != 0.0 or load.Cid() != 0:
                            normal = load.nvector / np.linalg.norm(load.nvector)
                            assert load.Cid() == 0, 'cid=%r on a PLOAD4 is not supported\n%s' % (load.Cid(), str(load))
                    else:
                        msg = 'surf_or_line=%r on PLOAD4 is not supported\n%s' % (
                            load.surf_or_line, str(load))
                        raise NotImplementedError(msg)

                    r = centroid - p
                    f = pressure * area * normal * scale
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

    def get_element_faces(self, element_ids=None, allow_blank_nids=True):
        """
        Gets the elements and faces that are skinned from solid elements.
        This includes internal faces, but not existing shells.

        Parameters
        ----------
        element_ids : List[int] / None
            skin a subset of element faces
            default=None -> all elements
        allow_blank_nids : bool; default=True
            allows for nids to be None

        Returns
        -------
        eid_faces : (int, List[(int, int, ...)])
           value1 : element id
           value2 : face
        """
        if element_ids is None:
            element_ids = self.element_ids

        eid_faces = []
        if allow_blank_nids:
            for eid in element_ids:
                elem = self.elements[eid]
                if elem.type in ['CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM']:
                    faces = elem.faces
                    for face_id, face in iteritems(faces):
                        eid_faces.append((eid, face))
        else:
            for eid in element_ids:
                elem = self.elements[eid]
                if elem.type in ['CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM']:
                    faces = elem.faces
                    for face_id, face in iteritems(faces):
                        if None in face:
                            msg = 'There is a None in the face.\n'
                            msg = 'face_id=%s face=%s\n%s' % (face_id, str(face), str(elem))
                            raise RuntimeError(msg)
                        eid_faces.append((eid, face))
        return eid_faces

    def get_solid_skin_faces(self):
        """
        Gets the elements and faces that are skinned from solid elements.
        This doesn't include internal faces or existing shells.

        Returns
        -------
        eid_set : Dict[tuple(int, int, ...)] = List[int]
           key : sorted face
           value : list of element ids with that face
        face_map : Dict[tuple(int, int, ...)] = List[int]
           key : sorted face
           value : unsorted face
        """
        eid_faces = self.get_element_faces()
        face_set = defaultdict(int)
        eid_set = defaultdict(list)
        face_map = {}
        for eid, face in eid_faces:
            #print(eid, face)
            raw_face = deepcopy(face)
            try:
                face.sort()
            except:
                print('face = %s' % str(face))
                raise
            tface = tuple(face)
            #print(tface)
            face_set[tface] += 1
            eid_set[tface].append(eid)
            face_map[tface] = raw_face

        #print('eid_set:')
        #for tface, eidset in iteritems(eid_set):
            #print(tface, eidset)

        #print('face_set:')
        #for tface, faceset in iteritems(face_set):
            #print(tface, faceset)

        #print('face_map:')
        #for tface, facemap in iteritems(face_map):
            #print(tface, facemap)

        del_faces = []
        for face, face_count in iteritems(face_set):
            if face_count == 2:
                del_faces.append(face)

        for face in del_faces:
            del face_set[face]
            del eid_set[face]
        return eid_set, face_map

    def write_skin_solid_faces(self, skin_filename,
                               write_solids=False, write_shells=True,
                               size=8, is_double=False, encoding=None):
        """
        Writes the skinned elements

        Parameters
        ----------
        skin_filename : str
            the file to write
        write_solids : bool; default=False
            write solid elements that have skinned faces
        write_shells : bool; default=False
            write newly created shell elements
            if there are shells in the model, doesn't write these

        size : int; default=8
            the field width
        is_double : bool; default=False
            double precision flag
        encoding : str; default=None -> system default
            the string encoding
        """
        if(len(self.element_ids) == 0 or len(self.material_ids) == 0 or
           len(self.property_ids) == 0):
            return
        eid_set, face_map = self.get_solid_skin_faces()
        if len(eid_set) == 0:
            return

        eid_set_to_write = set([])
        nid_set_to_write = set([])
        mid_set_to_write = set([])
        if write_solids:
            for face, eids in iteritems(eid_set):
                eid_set_to_write.update(eids)
                for eid in eids:
                    elem = self.elements[eid]
                    pid = elem.Pid()
                    prop = self.properties[pid] # PSOLID
                    mid = prop.Mid()
                    #print(prop)
                    nid_set_to_write.update(elem.node_ids)
                    mid_set_to_write.add(mid)
                    #print('added_mid (a) =', mid)
        elif write_shells:
            for face, eids in iteritems(eid_set):
                eid_set_to_write.update(eids)
                nid_set_to_write.update(face)
                for eid in eids:
                    elem = self.elements[eid]
                    pid = elem.Pid()
                    prop = self.properties[pid] # PSOLID
                    if prop.type in ['PSOLID', 'PLSOLID']:
                        mid = prop.Mid()
                    elif prop.type == 'PCOMPS':
                        mid = prop.mids[0]
                    else:
                        raise NotImplementedError(prop)
                    #except TypeError:
                        #self.log.warning('TypeError: skipping:%s' % prop)
                        #raise
                    #except AttributeError:
                        #self.log.warning('skipping:%s' % prop)
                        #continue
                    mid_set_to_write.add(mid)
                    #print('added eid=%s pid=%s mid=%s (b)' % (eid, pid, mid))
        else:
            raise RuntimeError('write_solids=False write_shells=False')

        eids_to_write = list(eid_set_to_write)
        nids_to_write = list(nid_set_to_write)
        mids_to_write = list(mid_set_to_write)

        #element_ids_to_delete = set(self.element_ids) - eids_to_write

        eid_shell = max(self.elements) + 1
        pid_shell = max(self.properties) + 1
        mid_shell = max(self.materials) + 1
        self._write_skin_solid_faces(skin_filename, face_map,
                                     nids_to_write, eids_to_write, mids_to_write, eid_set,
                                     eid_shell, pid_shell, mid_shell,
                                     write_solids=write_solids, write_shells=write_shells,
                                     size=size, is_double=is_double, encoding=encoding)

    def _write_skin_solid_faces(self, skin_filename, face_map,
                                nids_to_write, eids_to_write, mids_to_write, eid_set,
                                eid_shell, pid_shell, mid_shell,
                                write_solids=False, write_shells=True,
                                size=8, is_double=False, encoding=None):
        """
        helper method for ``write_skin_solid_faces``

        Parameters
        ----------
        skin_filename : str
            the file to write
        face_map : dict[sorted_face] : face
            sorted_face : List[int, int, int] / List[int, int, int, int]
            face : List[int, int, int] / List[int, int, int, int]

        nids_to_write : List[int, int, ...]
            list of node ids to write
        eids_to_write : List[int, int, ...]
            list of element ids to write
        mids_to_write : List[int, int, ...]
            list of material ids to write
        eid_set : ???
            ???

        eid_shell : int
            the next id to use for the shell id
        pid_shell : int
            the next id to use for the shell property
        mid_shell : int
            the next id to use for the shell material

        write_solids : bool; default=False
            write solid elements that have skinned faces
        write_shells : bool; default=True
            write shell elements
        size : int; default=8
            the field width
        is_double : bool; default=False
            double precision flag
        encoding : str; default=None -> system default
            the string encoding
        """
        encoding = self.get_encoding(encoding)
        if PY2:
            wb = 'wb'
        else:
            wb = 'w'
        with open(skin_filename, wb, encoding=encoding) as bdf_file:
            bdf_file.write('$ pyNastran: punch=True\n')
            for nid in sorted(nids_to_write):
                if nid is None:
                    continue
                node = self.nodes[nid]
                bdf_file.write(node.write_card(size=size, is_double=is_double))

            for cid, coord in iteritems(self.coords):
                if cid == 0:
                    continue
                bdf_file.write(coord.write_card(size=size, is_double=is_double))

            if write_solids:
                for eid in sorted(eids_to_write):
                    elem = self.elements[eid]
                    bdf_file.write(elem.write_card(size=size))
                for pid, prop in iteritems(self.properties):
                    bdf_file.write(prop.write_card(size=size, is_double=is_double))
                for mid in sorted(mids_to_write):
                    material = self.materials[mid]
                    bdf_file.write(material.write_card(size=size, is_double=is_double))
                del eid, pid, mid

            if write_shells:
                mids_to_write.sort()
                for imid, mid in enumerate(mids_to_write):
                    card = ['PSHELL', pid_shell + imid, mid_shell + imid, 0.1]
                    try:
                        msg = print_card_8(card)
                    except RuntimeError:
                        msg = print_card_16(card)
                    bdf_file.write(msg)

                    card = ['MAT1', mid_shell + imid, 3.e7, None, 0.3]
                    #bdf_file.write(self.materials[mid].comment)
                    try:
                        msg = print_card_8(card)
                    except RuntimeError:
                        msg = print_card_16(card)
                    bdf_file.write(msg)

                for face, eids in iteritems(eid_set):
                    face_raw = face_map[face]
                    nface = len(face)
                    #assert len(eids) == 1, eids
                    #for eid in sorted(eids):
                        #elem = self.elements[eid]
                        #print(elem)
                        #break

                    #elem = next(itervalues(self.elements)) # old
                    elem = self.elements[eids[0]]
                    #pid = next(iterkeys(self.properties))
                    pid = elem.Pid()
                    prop = self.properties[pid]
                    try:
                        mid = prop.Mid()
                    except AttributeError:
                        continue
                    #print('mids_to_write = %s' % mids_to_write)
                    #print('mids = ', self.materials.keys())
                    imid = mids_to_write.index(mid)

                    if nface == 3:
                        card = ['CTRIA3', eid_shell, pid_shell + imid] + list(face_raw)
                    elif nface == 4:
                        card = ['CQUAD4', eid_shell, pid_shell + imid] + list(face_raw)
                    elif nface == 4:
                        card = ['CQUAD4', eid_shell, pid_shell + imid] + list(face_raw)
                    elif nface == 6:
                        card = ['CTRIA6', eid_shell, pid_shell + imid] + list(face_raw)
                    elif nface == 8:
                        card = ['CQUAD8', eid_shell, pid_shell + imid] + list(face_raw)
                    else:
                        raise NotImplementedError('face=%s len(face)=%s' % (face, nface))
                    try:
                        msg = print_card_8(card)
                    except RuntimeError:
                        msg = print_card_16(card)
                    bdf_file.write(msg)
                    eid_shell += 1

                    #elem = self.elements[eid]
                    #bdf_file.write(elem.write_card(size=size))
                #for pid, prop in iteritems(self.properties):
                    #bdf_file.write(prop.write_card(size=size, is_double=is_double))
            bdf_file.write('ENDDATA\n')
        #if 0:
            #model = self.__class__.__init__()
            #model.read_bdf(skin_filename)
