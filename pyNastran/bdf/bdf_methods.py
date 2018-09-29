"""
This file contains additional methods that do not directly relate to the
reading/writing/accessing of BDF data.  Such methods include:
  - mass_poperties
      get the mass & moment of inertia of the model
  - sum_forces_moments
      find the net force/moment on the model
  - sum_forces_moments_elements
      find the net force/moment on the model for a subset of elements
  - resolve_grids
      change all nodes to a specific coordinate system
  - unresolve_grids
      puts all nodes back to original coordinate system
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from collections import defaultdict
from typing import List, Tuple, Any, Union, Dict

import numpy as np

from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
from pyNastran.bdf.mesh_utils.mass_properties import (
    mass_properties, mass_properties_no_xref, mass_properties_nsm)
from pyNastran.bdf.mesh_utils.loads import sum_forces_moments, sum_forces_moments_elements
from pyNastran.bdf.mesh_utils.skin_solid_elements import write_skin_solid_faces


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

    def get_length_breakdown(self, property_ids=None, stop_if_no_length=True):
        """
        gets a breakdown of the length by property region

        TODO: What about CONRODs?
        """
        #skip_elements = [
            #'CTRIA3', 'CTRIA6', 'CTRIAR',
            #'CQUAD4', 'CQUAD8',
            #'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
            #'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
            #'CBUSH', 'CBUSH1D', 'CBUSH2D',
            #'CRAC2D', 'CRAC3D',
        #]
        skip_props = [
            'PSOLID', 'PLPLANE', 'PPLANE', 'PELAS',
            'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
            'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
            'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX', 'PLSOLID',
            'PCOMPS', 'PVISC',
            'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
            'PBEND',
        ]
        bar_properties = ['PBAR', 'PBARL', 'PBEAM', 'PBEAML',
                          'PROD', 'PTUBE', 'PBRSECT', 'PBMSECT', 'PBCOMP']
        pid_eids = self.get_element_ids_dict_with_pids(
            property_ids, msg=' which is required by get_length_breakdown')
        pids_to_length = {}
        for pid, eids in pid_eids.items():
            prop = self.properties[pid]
            lengths = []
            if prop.type in skip_props:
                continue
            elif prop.type in bar_properties:
                #['CBAR', 'CBEAM', 'CROD', 'CTUBE']:
                # TODO: Do I need to consider the offset on length effects for a CBEAM?
                for eid in eids:
                    elem = self.elements[eid]
                    try:
                        lengths.append(elem.Length())
                    except AttributeError:  # pragma: no cover
                        print(prop)
                        print(elem)
                        raise
            else:  # pragma: no cover
                print('prop\n%s' % prop)
                eid0 = eids[0]
                elem = self.elements[eid0]
                msg = str(prop) + str(elem)
                raise NotImplementedError(msg)

            if lengths:
                pids_to_length[pid] = sum(lengths)

        has_length = len(pids_to_length)

        if not has_length:
            msg = 'No elements with length were found'
            self.log.warning(msg)
            if stop_if_no_length:
                raise RuntimeError(msg)
        return pids_to_length

    def get_area_breakdown(self, property_ids=None, stop_if_no_area=True, sum_bar_area=True):
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

        sum_bar_area : bool; default=True
            sum the areas for CBAR/CBEAM/CROD/CONROD/CTUBE elements
            True : get the area of the model by property id
            False : only get the cross sectional properties
        """
        skip_props = [
            'PSOLID', 'PLPLANE', 'PPLANE', 'PELAS',
            'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
            'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
            'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX', 'PLSOLID',
            'PCOMPS', 'PVISC', 'PBCOMP', 'PBEND',
        ]

        pid_eids = self.get_element_ids_dict_with_pids(
            property_ids, msg=' which is required by get_area_breakdown')
        pids_to_area = {}
        for pid, eids in pid_eids.items():
            prop = self.properties[pid]
            areas = []
            if prop.type in ['PSHELL', 'PCOMP', 'PSHEAR', 'PCOMPG', ]:
                for eid in eids:
                    elem = self.elements[eid]
                    if elem.type in ['CQUADX']:
                        continue
                    try:
                        areas.append(elem.Area())
                    except AttributeError:  # pragma: no cover
                        print(prop)
                        print(elem)
                        raise
            elif prop.type in ['PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE']:
                for eid in eids:
                    elem = self.elements[eids[0]]
                    area = elem.Area()
                    if sum_bar_area:
                        neids = len(eids)
                        areas = [area * neids]
                    else:
                        areas = [area]
            elif prop.type in skip_props:
                pass
            elif prop.type in ['PBRSECT', 'PBMSECT']:
                self.log.warning('skipping:\n%s' % prop)
                continue
            else:  # pragma: no cover
                raise NotImplementedError(prop)
            if areas:
                pids_to_area[pid] = sum(areas)

        has_area = len(pids_to_area)
        if not has_area:
            msg = 'No elements with area were found'
            self.log.warning(msg)
            if stop_if_no_area:
                raise RuntimeError(msg)
        return pids_to_area

    def get_volume_breakdown(self, property_ids=None, stop_if_no_volume=True):
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
        pid_eids = self.get_element_ids_dict_with_pids(
            property_ids, msg=' which is required by get_area_breakdown')

        no_volume = [
            'PLPLANE', 'PPLANE', 'PELAS',
            'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
            'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
            'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX',
            'PVISC', 'PBCOMP', 'PBEND',
        ]
        pids_to_volume = {}
        skipped_eid_pid = set([])
        for pid, eids in pid_eids.items():
            prop = self.properties[pid]
            volumes = []
            if prop.type == 'PSHELL':
                # TODO: doesn't support PSHELL differential thicknesses
                thickness = prop.t
                areas = []
                for eid in eids:
                    elem = self.elements[eid]
                    if elem.type in ['CQUADX']:
                        continue
                    areas.append(elem.Area())
                volumesi = [area * thickness for area in areas]
                volumes.extend(volumesi)
            elif prop.type in ['PCOMP', 'PCOMPG',]:
                areas = []
                for eid in eids:
                    elem = self.elements[eid]
                    areas.append(elem.Area())
                thickness = prop.Thickness()
                volumesi = [area * thickness for area in areas]
                volumes.extend(volumesi)
            elif prop.type in ['PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE']:
                # TODO: Do I need to consider the offset on length effects for a CBEAM?
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
                    if elem.type in ['CTETRA', 'CPENTA', 'CHEXA']:
                        volumes.append(elem.Volume())
                    else:
                        key = (elem.type, prop.type)
                        if key not in skipped_eid_pid:
                            skipped_eid_pid.add(key)
                            self.log.debug('skipping volume %s' % str(key))
            elif prop.type == 'PSHEAR':
                thickness = prop.t
                areas = []
                for eid in eids:
                    elem = self.elements[eid]
                    areas.append(elem.Area())
                volumesi = [area * thickness for area in areas]
                volumes.extend(volumesi)
            elif prop.type in no_volume:
                pass
            elif prop.type in ['PBRSECT', 'PBMSECT']:
                self.log.warning('skipping:\n%s' % prop)
                continue
            else:  # pragma: no cover
                raise NotImplementedError(prop)
            if volumes:
                pids_to_volume[pid] = sum(volumes)

        has_volume = len(pids_to_volume)
        if not has_volume:
            msg = 'No elements with volume were found'
            self.log.warning(msg)
            if stop_if_no_volume:
                raise RuntimeError(msg)
        return pids_to_volume

    def get_mass_breakdown(self, property_ids=None, stop_if_no_mass=True, detailed=False):
        """
        gets a breakdown of the mass by property region

        Parameters
        ----------
        property_ids : List[int] / int
            list of property ID
        stop_if_no_mass : bool; default=True
            prevents crashing if there are no elements
            setting this to False really doesn't make sense for non-DMIG models
        detailed : bool, optional, default : False
            Separates structural and nonstructural mass outputs.

        Returns
        --------
        pids_to_mass : dict {int : float, ...}
            Map from property id to mass (structural mass only if detailed is True).
        pids_to_mass_nonstructural : dict {int : float, ...}, optional
            Map from property id to nonstructural mass (only if detailed is True).
        mass_type_to_mass : dict {str : float, ...}
            Map from mass id to mass for mass elements.

        TODO: What about CONRODs, CONM2s?
        #'PBCOMP',
        #'PBMSECT',
        #'PBEAM3',
        #'PBEND',
        #'PIHEX',
        #'PCOMPS',
        """
        pid_eids = self.get_element_ids_dict_with_pids(
            property_ids, stop_if_no_eids=False, msg=' which is required by get_area_breakdown')

        mass_type_to_mass = {}
        pids_to_mass = {}
        pids_to_mass_nonstructural = {}
        skipped_eid_pid = set([])
        for eid, elem in self.masses.items():
            if elem.type not in mass_type_to_mass:
                mass_type_to_mass[elem.type] = elem.Mass()
            else:
                mass_type_to_mass[elem.type] += elem.Mass()

        properties_to_skip = [
            'PLPLANE', 'PPLANE', 'PELAS',
            'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
            'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
            'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX',
            'PVISC', 'PBCOMP', 'PBEND', ]
        for pid, eids in pid_eids.items():
            prop = self.properties[pid]
            masses = []
            masses_nonstructural = []
            if prop.type == 'PSHELL':
                # TODO: doesn't support PSHELL differential thicknesses
                thickness = prop.t
                nsm = prop.nsm  # per area
                rho = prop.Rho()
                for eid in eids:
                    elem = self.elements[eid]
                    if elem.type in ['CQUADX']:
                        continue
                    area = elem.Area()
                    if detailed:
                        masses.append(area * (rho * thickness))
                        masses_nonstructural.append(area * nsm)
                    else:
                        masses.append(area * (rho * thickness + nsm))
            elif prop.type in ['PCOMP', 'PCOMPG']:
                # TODO: does the PCOMP support differential thickness?
                #       I don't think so...
                for eid in eids:
                    elem = self.elements[eid]
                    if detailed:
                        mass, mass_nonstructural = elem.Mass_breakdown()
                        masses.append(mass)
                        masses_nonstructural.append(mass_nonstructural)
                    else:
                        masses.append(elem.Mass())

            elif prop.type in ['PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE']:
                nsm = prop.nsm # per unit length
                try:
                    rho = prop.Rho()
                except AttributeError:
                    print(prop)
                    raise
                for eid in eids:
                    elem = self.elements[eid]
                    area = prop.Area()
                    length = elem.Length()
                    if detailed:
                        masses.append(length * (rho * area))
                        masses_nonstructural.append(length * nsm)
                    else:
                        masses.append(length * (rho * area + nsm))
            elif prop.type in ['PSOLID', 'PCOMPS', 'PLSOLID']:
                rho = prop.Rho()
                for eid in eids:
                    elem = self.elements[eid]
                    if elem.type in ['CTETRA', 'CPENTA', 'CHEXA']:
                        masses.append(rho * elem.Volume())
                    else:
                        key = (elem.type, prop.type)
                        if key not in skipped_eid_pid:
                            skipped_eid_pid.add(key)
                            self.log.debug('skipping mass %s' % str(key))
            elif prop.type in properties_to_skip:
                pass
            elif prop.type == 'PSHEAR':
                thickness = prop.t
                nsm = prop.nsm # per area
                rho = prop.Rho()
                for eid in eids:
                    elem = self.elements[eid]
                    area = elem.Area()
                    if detailed:
                        masses.append(area * (rho * thickness))
                        masses_nonstructural.append(area * nsm)
                    else:
                        masses.append(area * (rho * thickness + nsm))
            elif prop.type in ['PBRSECT', 'PBMSECT']:
                self.log.warning('skipping:\n%s' % prop)
                continue
            else:  # pragma: no cover
                raise NotImplementedError(prop)
            if masses:
                pids_to_mass[pid] = sum(masses)
            if masses_nonstructural:
                pids_to_mass_nonstructural[pid] = sum(masses_nonstructural)

        has_mass = len(mass_type_to_mass) > 0 or len(pids_to_mass) > 0
        if not has_mass:
            msg = 'No elements with mass were found'
            self.log.warning(msg)
            if stop_if_no_mass:
                raise RuntimeError(msg)
        if detailed:
            return pids_to_mass, pids_to_mass_nonstructural, mass_type_to_mass
        else:
            return pids_to_mass, mass_type_to_mass

    def mass_properties(self, element_ids=None, mass_ids=None,
                        reference_point=None,
                        sym_axis=None, scale=None, inertia_reference='cg'):
        """
        Calculates mass properties in the global system about the
        reference point.

        Parameters
        ----------
        element_ids : list[int]; (n, ) ndarray, optional
            An array of element ids.
        mass_ids : list[int]; (n, ) ndarray, optional
            An array of mass ids.
        reference_point : ndarray/int, optional
            type : ndarray
                An array that defines the origin of the frame.
                default = <0,0,0>.
            type : int
                the node id
        sym_axis : str, optional
            The axis to which the model is symmetric.
            If AERO cards are used, this can be left blank.
            allowed_values = 'no', x', 'y', 'z', 'xy', 'yz', 'xz', 'xyz'
        scale : float, optional
            The WTMASS scaling value.
            default=None -> PARAM, WTMASS is used
            float > 0.0
        inertia_reference : str; default='cg'
            'cg' : inertia is taken about the cg
            'ref' : inertia is about the reference point

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

        Examples
        --------
        Mass properties of entire structure

        >>> mass, cg, I = model.mass_properties()
        >>> Ixx, Iyy, Izz, Ixy, Ixz, Iyz = I

        Mass properties of model based on Property ID

        >>> pids = list(model.pids.keys())
        >>> pid_eids = self.get_element_ids_dict_with_pids(pids)
        >>> for pid, eids in sorted(pid_eids.items()):
        >>>     mass, cg, I = model.mass_properties(element_ids=eids)
        """
        mass, cg, I = mass_properties(
            self,
            element_ids=element_ids, mass_ids=mass_ids,
            reference_point=reference_point,
            sym_axis=sym_axis, scale=scale,
            inertia_reference=inertia_reference)
        return mass, cg, I

    def mass_properties_no_xref(self, element_ids=None, mass_ids=None,
                                reference_point=None,
                                sym_axis=None, scale=None, inertia_reference='cg'):
        """
        Calculates mass properties in the global system about the
        reference point.

        Parameters
        ----------
        element_ids : list[int]; (n, ) ndarray, optional
            An array of element ids.
        mass_ids : list[int]; (n, ) ndarray, optional
            An array of mass ids.
        reference_point : ndarray/int, optional
            type : ndarray
                An array that defines the origin of the frame.
                default = <0,0,0>.
            type : int
                the node id
        sym_axis : str, optional
            The axis to which the model is symmetric.
            If AERO cards are used, this can be left blank.
            allowed_values = 'no', x', 'y', 'z', 'xy', 'yz', 'xz', 'xyz'
        scale : float, optional
            The WTMASS scaling value.
            default=None -> PARAM, WTMASS is used
            float > 0.0
        inertia_reference : str; default='cg'
            'cg' : inertia is taken about the cg
            'ref' : inertia is about the reference point

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

        Examples
        --------
        **mass properties of entire structure**

        >>> mass, cg, I = model.mass_properties()
        >>> Ixx, Iyy, Izz, Ixy, Ixz, Iyz = I


        **mass properties of model based on Property ID**

        >>> pids = list(model.pids.keys())
        >>> pid_eids = self.get_element_ids_dict_with_pids(pids)
        >>> for pid, eids in sorted(pid_eids.items()):
        >>>     mass, cg, I = model.mass_properties(element_ids=eids)
        """
        mass, cg, I = mass_properties_no_xref(
            self, element_ids=element_ids, mass_ids=mass_ids,
            reference_point=reference_point,
            sym_axis=sym_axis, scale=scale,
            inertia_reference=inertia_reference)
        return mass, cg, I

    def mass_properties_nsm(self, element_ids=None, mass_ids=None, nsm_id=None,
                            reference_point=None,
                            sym_axis=None, scale=None, inertia_reference='cg',
                            xyz_cid0_dict=None, debug=False):
        """
        Calculates mass properties in the global system about the
        reference point.  Considers NSM, NSM1, NSML, NSML1.

        Parameters
        ----------
        model : BDF()
            a BDF object
        element_ids : list[int]; (n, ) ndarray, optional
            An array of element ids.
        mass_ids : list[int]; (n, ) ndarray, optional
            An array of mass ids.
        nsm_id : int
            the NSM id to consider
        reference_point : ndarray/int, optional
            type : ndarray
                An array that defines the origin of the frame.
                default = <0,0,0>.
            type : int
                the node id
        sym_axis : str, optional
            The axis to which the model is symmetric.
            If AERO cards are used, this can be left blank.
            allowed_values = 'no', x', 'y', 'z', 'xy', 'yz', 'xz', 'xyz'
        scale : float, optional
            The WTMASS scaling value.
            default=None -> PARAM, WTMASS is used
            float > 0.0
        inertia_reference : str; default='cg'
            'cg' : inertia is taken about the cg
            'ref' : inertia is about the reference point
        xyz_cid0_dict : dict[nid] : xyz; default=None -> auto-calculate
            mapping of the node id to the global position
        debug : bool; default=False
            developer debug; may be removed in the future

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

        Examples
        --------
        **mass properties of entire structure**

        >>> mass, cg, I = model.mass_properties()
        >>> Ixx, Iyy, Izz, Ixy, Ixz, Iyz = I


        **mass properties of model based on Property ID**
        >>> pids = list(model.pids.keys())
        >>> pid_eids = model.get_element_ids_dict_with_pids(pids)
        >>> for pid, eids in sorted(pid_eids.items()):
        >>>     mass, cg, I = mass_properties(model, element_ids=eids)

        Warnings
        --------
         - If eids are requested, but don't exist, no warning is thrown.
           Decide if this is the desired behavior.
         - If the NSMx ALL option is used, the mass from all elements
           will be considered, even if not included in the element set

        """
        mass, cg, I = mass_properties_nsm(
            self, element_ids=element_ids, mass_ids=mass_ids, nsm_id=nsm_id,
            reference_point=reference_point,
            sym_axis=sym_axis, scale=scale,
            inertia_reference=inertia_reference,
            xyz_cid0_dict=xyz_cid0_dict, debug=debug)
        return (mass, cg, I)

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
        # type: (int, int, List[int], List[int], bool, Union[None, Dict[int, np.ndarray]]) -> Tuple[np.ndarray, np.ndarray]
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
        forces : NUMPY.NDARRAY shape=(3,)
            the forces
        moments : NUMPY.NDARRAY shape=(3,)
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

          Notes
        -----
        If you split the model into sections and sum the loads
        on each section, you may not get the same result as
        if you summed the loads on the total model.  This is
        due to the fact that nodal loads on the boundary are
        double/triple/etc. counted depending on how many breaks
        you have.

        .. todo:: not done...
        """
        forces, moments = sum_forces_moments_elements(self, p0, loadcase_id, eids, nids,
                                                      include_grav=include_grav, xyz_cid0=xyz_cid0)
        return forces, moments

    def sum_forces_moments(self, p0, loadcase_id, include_grav=False,
                           xyz_cid0=None):
        # type: (int, int, bool, Union[None, Dict[int, np.ndarray]]) -> Tuple[np.ndarray, np.ndarray]
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
        forces : NUMPY.NDARRAY shape=(3,)
            the forces
        moments : NUMPY.NDARRAY shape=(3,)
            the moments

        .. warning:: not full validated
        .. todo:: It's super slow for cid != 0.   We can speed this up a lot
                 if we calculate the normal, area, centroid based on
                 precomputed node locations.

        Pressure acts in the normal direction per model/real/loads.bdf and loads.f06
        """
        forces, moments = sum_forces_moments(self, p0, loadcase_id,
                                             include_grav=include_grav, xyz_cid0=xyz_cid0)
        return forces, moments

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
                    for face_id, face in faces.items():
                        eid_faces.append((eid, face))
        else:
            for eid in element_ids:
                elem = self.elements[eid]
                if elem.type in ['CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM']:
                    faces = elem.faces
                    for face_id, face in faces.items():
                        if None in face:
                            msg = 'There is a None in the face.\n'
                            msg = 'face_id=%s face=%s\n%s' % (face_id, str(face), str(elem))
                            raise RuntimeError(msg)
                        eid_faces.append((eid, face))
        return eid_faces

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
        return write_skin_solid_faces(
            self, skin_filename,
            write_solids=write_solids, write_shells=write_shells,
            size=size, is_double=is_double, encoding=encoding)

    def update_model_by_desvars(self, xref=True, desvar_values=None):
        """doesn't require cross referenceing"""
        # these are the nominal values of the desvars
        desvar_init = {key : desvar.value
                       for key, desvar in self.desvars.items()}

        # these are the current values of the desvars
        if desvar_values is None:
            desvar_values = {key : min(max(desvar.value + 0.1, desvar.xlb), desvar.xub)
                             for key, desvar in self.desvars.items()}

        # Relates one design variable to one or more other design variables.
        for dlink_id, dlink in self.dlinks.items():
            value = dlink.c0
            desvar = dlink.dependent_desvar #get_stats()
            desvar_ref = self.desvars[desvar]
            for coeff, desvar_idi in zip(dlink.coeffs, dlink.IDv):
                valuei = desvar_values[desvar_idi]
                value += coeff * valuei
            value2 = min(max(value, desvar_ref.xlb), desvar_ref.xub)
            desvar_values[dlink_id] = value2

        # calculates the real delta to be used by DVGRID
        desvar_delta = {key : (desvar_init[key] - desvar_values[key])
                        for key in self.desvars}


        #min(max(self.xinit, self.xlb), self.xub)

        # DVxREL1
        dvxrel2s = {}
        for dvid, dvprel in self.dvprels.items():
            if dvprel.type == 'DVPREL2':
                dvxrel2s[('DVPREL2', dvid)] = dvprel
                continue
            dvprel.update_model(self, desvar_values)

        for dvid, dvmrel in self.dvmrels.items():
            if dvmrel.type == 'DVPREL2':
                dvxrel2s[('DVMREL2', dvid)] = dvmrel
                continue
            dvmrel.update_model(self, desvar_values)

        for dvid, dvcrel in self.dvcrels.items():
            if dvcrel.type == 'DVPREL2':
                dvxrel2s[('DVMREL2', dvid)] = dvcrel
                continue
            dvcrel.update_model(self, desvar_values)

        #+--------+------+-----+-----+-------+----+----+----+
        #|    1   |   2  |  3  |  4  |   5   |  6 |  7 |  8 |
        #+========+======+=====+=====+=======+====+====+====+
        #| DVGRID | DVID | GID | CID | COEFF | N1 | N2 | N3 |
        #+--------+------+-----+-----+-------+----+----+----+

        # grid_i - grid_i0 = sum(coeffj * (x_desvar_j - x0_desvar_j)) * {Nxyz_f}
        dxyzs = defaultdict(list)
        for dvid, dvgrids in self.dvgrids.items():
            for dvgrid in dvgrids:
                dxyz_cid = dvgrid.coeff * desvar_delta[dvid] * dvgrid.dxyz
                dxyzs[(dvgrid.nid, dvgrid.cid)].append(dxyz_cid)

        # TODO: could be vectorized
        for (nid, cid), dxyz in dxyzs.items():
            dxyz2 = np.linalg.norm(dxyz, axis=0)
            assert len(dxyz2) == 3, len(dxyz2)
            grid = self.nodes[nid]
            coord_from = self.coords[cid]
            coord_to = self.coords[grid.cp]
            grid.xyz += coord_from.transform_node_from_local_to_local(
                coord_to, dxyz2)

        if xref:
            for key, dvxrel2 in dvxrel2s.items():
                dvxrel2.update_model(self, desvar_values)
        #self.nid = nid
        #self.cid = cid
        #self.coeff = coeff
        #self.dxyz = np.asarray(dxyz)
        #dvgrid.dxyz
