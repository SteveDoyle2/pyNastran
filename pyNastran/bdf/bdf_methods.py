"""
This file contains additional methods that do not directly relate to the
reading/writing/accessing of BDF data.  Such methods include:
  - resolve_grids
      change all nodes to a specific coordinate system
  - unresolve_grids
      puts all nodes back to original coordinate system

"""
from collections import defaultdict
from typing import Any, Optional

import numpy as np

from pyNastran.bdf.bdf_interface.attributes import BDFAttributes
from pyNastran.bdf.mesh_utils.breakdowns import (
    get_length_breakdown, get_area_breakdown, get_volume_breakdown, get_mass_breakdown)

from pyNastran.nptyping_interface import NDArray3float


class BDFMethods(BDFAttributes):
    """
    Has the following methods:
        resolve_grids(cid=0)
        unresolve_grids(model_old)

    """
    def __init__(self):
        BDFAttributes.__init__(self)

    def get_length_breakdown(self, property_ids=None, stop_if_no_length: bool=True):
        """
        gets a breakdown of the length by property region

        TODO: What about CONRODs?

        """
        return get_length_breakdown(self, property_ids=property_ids, stop_if_no_length=stop_if_no_length)

    def get_area_breakdown(self, property_ids=None,
                           stop_if_no_area: bool=True,
                           sum_bar_area: bool=True) -> dict[int, float]:
        """
        Gets a breakdown of the area by property region

        Parameters
        ----------
        property_ids : list[int] / int
            list of property ID
        stop_if_no_area : bool; default=True
            prevents crashing if there are no elements
        sum_bar_area : bool; default=True
            sum the areas for CBAR/CBEAM/CROD/CONROD/CTUBE elements
            True : get the area of the model by property id (e.g., A=A_pid*Nelements)
            False : only get the cross sectional properties (e.g., A=A_pid)

        Returns
        -------
        pids_to_area : dict[int pid] : float area
            the pid to area dictionary

        TODO: What about CONRODs?
        #'PBRSECT', 'PBCOMP', 'PBMSECT', 'PBEAM3', 'PBEND', 'PCOMPS',

        """
        pids_to_area = get_area_breakdown(self, property_ids=property_ids,
                                          stop_if_no_area=stop_if_no_area,
                                          sum_bar_area=sum_bar_area)
        return pids_to_area

    def get_volume_breakdown(self, property_ids: Optional[int]=None,
                             stop_if_no_volume: bool=True) -> dict[int, float]:
        """
        gets a breakdown of the volume by property region

        TODO: What about CONRODs?
        #'PBRSECT',
        #'PBCOMP',
        #'PBMSECT',
        #'PBEAM3',
        #'PBEND',

        """
        return get_volume_breakdown(self, property_ids=property_ids,
                                    stop_if_no_volume=stop_if_no_volume)

    def get_mass_breakdown(self,
                           property_ids: Optional[int]=None,
                           stop_if_no_mass: bool=True) -> tuple[
                               dict[int, float],
                               dict[str, float]]:
        """
        Gets a breakdown of the mass by property region.

        Parameters
        ----------
        property_ids : list[int] / int
            list of property ID
        stop_if_no_mass : bool; default=True
            prevents crashing if there are no elements
            setting this to False really doesn't make sense for non-DMIG models

        Returns
        -------
        pids_to_mass : dict {int : float, ...}
            Map from property id to mass.
        mass_type_to_mass : dict {str : float, ...}
            Map from mass id to mass for mass elements.
            CONM2s are listed here

        TODO: What about CONRODs?
        #'PBCOMP', 'PBMSECT', 'PBEAM3', 'PBEND', 'PCOMPS',

        """
        pids_to_mass, mass_type_to_mass = get_mass_breakdown(
            self, property_ids=property_ids,
            stop_if_no_mass=stop_if_no_mass, detailed=False)
        return pids_to_mass, mass_type_to_mass

    def get_mass_breakdown_detailed(self,
                                    property_ids: Optional[int]=None,
                                    stop_if_no_mass: bool=True) -> tuple[
                                        dict[int, float],
                                        dict[int, float],
                                        dict[str, float]]:
        """
        Gets a breakdown of the mass by property region.

        Parameters
        ----------
        property_ids : list[int] / int
            list of property ID
        stop_if_no_mass : bool; default=True
            prevents crashing if there are no elements
            setting this to False really doesn't make sense for non-DMIG models

        Returns
        -------
        pids_to_mass : dict {int : float, ...}
            Map from property id to structural mass.
        pids_to_mass_nonstructural : dict {int : float, ...}
            Map from property id to nonstructural mass.
        mass_type_to_mass : dict {str : float, ...}
            Map from mass id to mass for mass elements.
            CONM2s are listed here

        TODO: What about CONRODs?
        #'PBCOMP', 'PBMSECT', 'PBEAM3', 'PBEND', 'PCOMPS',

        """
        # self.deprecated(
        #     'model.get_mass_breakdown(...)',
        #     'from pyNastran.bdf.mesh_utils.breakdowns import get_mass_breakdown' \
        #     'get_mass_breakdown', '1.5')
        pids_to_mass, pids_to_mass_nonstructural, mass_type_to_mass = get_mass_breakdown(
            self, property_ids=property_ids,
            stop_if_no_mass=stop_if_no_mass, detailed=True)
        return pids_to_mass, pids_to_mass_nonstructural, mass_type_to_mass

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
        #mass, cg, I = mass_properties(self, reference_point=p0, sym_axis=None)

    def get_element_faces(self,
                          element_ids: Optional[list[int]]=None,
                          allow_blank_nids: bool=True) -> Any:
        """
        Gets the elements and faces that are skinned from solid elements.
        This includes internal faces, but not existing shells.

        Parameters
        ----------
        element_ids : list[int] / None
            skin a subset of element faces
            default=None -> all elements
        allow_blank_nids : bool; default=True
            allows for nids to be None

        Returns
        -------
        eid_faces : (int, list[(int, int, ...)])
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

    def update_model_by_desvars(self, xref=True, desvar_values=None):
        """doesn't require cross referenceing"""
        # these are the nominal values of the desvars
        desvar_init = {key: desvar.value
                       for key, desvar in self.desvars.items()}

        # these are the current values of the desvars
        if desvar_values is None:
            desvar_values = {key: min(max(desvar.value + 0.1, desvar.xlb), desvar.xub)
                             for key, desvar in self.desvars.items()}

        # Relates one design variable to one or more other design variables.
        for dlink_id, dlink in self.dlinks.items():
            value = dlink.c0
            desvar = dlink.dependent_desvar
            desvar_ref = self.desvars[desvar]
            for coeff, desvar_idi in zip(dlink.coeffs, dlink.independent_desvars):
                valuei = desvar_values[desvar_idi]
                value += coeff * valuei
            value2 = min(max(value, desvar_ref.xlb), desvar_ref.xub)
            desvar_values[dlink_id] = value2

        # calculates the real delta to be used by DVGRID
        desvar_delta = {key: (desvar_init[key] - desvar_values[key])
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
