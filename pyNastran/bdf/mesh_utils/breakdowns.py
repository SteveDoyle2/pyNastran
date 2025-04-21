"""
defines:
 - pid_to_length = get_length_breakdown(
       model, property_ids=None,
       stop_if_no_length=True)
 - pid_to_area = get_area_breakdown(
       model, property_ids=None,
       stop_if_no_area=True, sum_bar_area=True)
 - pid_to_volume = get_volume_breakdown(
       model, property_ids=None,
       stop_if_no_volume=True)
 - pids_to_mass, mass_type_to_mass = get_mass_breakdown(
       model, property_ids=None,
       stop_if_no_mass=True, detailed=False)
 - pids_to_mass, pids_to_mass_nonstructural, mass_type_to_mass = get_mass_breakdown(
       model, property_ids=None,
       stop_if_no_mass=True, detailed=True)

"""
from __future__ import annotations
from typing import Any, TYPE_CHECKING
from collections import defaultdict
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

def get_material_mass_breakdown_table(model: BDF) -> tuple[dict[int, float],
                                                           dict[int, dict[int, float]],
                                                           dict[int, dict[int, float]],
                                                           dict[int, dict[int, float]],
                                                           ]:
    mid_to_mass = defaultdict(float)
    pid_to_mid_to_mass = defaultdict(lambda: defaultdict(float))
    pid_to_mid_to_nsm_mass = defaultdict(lambda: defaultdict(float))
    pid_to_mid_to_rho_mass = defaultdict(lambda: defaultdict(float))
    line_elements = {'CBAR', 'CBEAM', 'CTUBE', 'CROD'}
    solid_elements = {'CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM'}
    shell_elements = {'CQUAD4', 'CTRIA3', 'CQUAD8', 'CQUADR'}
    skip_elements = {'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CGAP', 'CRAC2D', 'CRAC3D'}

    for eid, elem in model.elements.items():
        elem_type = elem.type
        if elem_type in skip_elements:
            continue

        if elem_type in line_elements:
            prop = elem.pid_ref
            length = elem.Length()
            pid = prop.pid
            mid = prop.Mid()
            area = prop.Area()
            nsm = prop.Nsm()
            rho = prop.Rho()
            massi = elem.Mass()
            mid_to_mass[mid] += massi
            pid_to_mid_to_mass[pid][mid] += massi
            pid_to_mid_to_nsm_mass[pid][mid] += length * nsm
            pid_to_mid_to_rho_mass[pid][mid] += length * rho * area
        elif elem_type in solid_elements:
            prop = elem.pid_ref
            #volume = elem.Volume()
            pid = prop.pid
            mid = prop.Mid()
            massi = elem.Mass()
            mid_to_mass[mid] += massi
            pid_to_mid_to_mass[pid][mid] += massi
            pid_to_mid_to_nsm_mass[pid][mid] += 0.
            pid_to_mid_to_rho_mass[pid][mid] += massi
        elif elem_type in {'CONROD'}: # references a MAT directly
            mid = prop.Mid()
            length = elem.Length()
            massi = elem.Mass()
            nsm = elem.Nsm()
            rho = elem.Area()
            rho = elem.Rho()
            mid_to_mass[mid] += massi
            pid_to_mid_to_mass[0][mid] += massi
            pid_to_mid_to_nsm_mass[0][mid] += length * nsm
            pid_to_mid_to_rho_mass[0][mid] += length * rho * area
        elif elem_type in shell_elements:
            prop = elem.pid_ref
            pid = prop.pid
            prop_type = prop.type
            area = elem.Area()
            nsm = prop.Nsm()
            if prop_type in {'PSHELL'}:
                mid = prop.Mid()
                massi = elem.Mass()
                mid_to_mass[mid] += massi
                pid_to_mid_to_mass[pid][mid] += massi
                pid_to_mid_to_nsm_mass[pid][mid] += area * nsm
                pid_to_mid_to_rho_mass[pid][mid] += area * rho * t
            elif prop_type in {'PCOMP', 'PCOMPG'}:
                nlayers = len(prop.mids)
                pid_to_mid_to_nsm_mass[pid][0] += area * nsm
                nsm_n = nsm / nlayers
                for mid, t in zip(prop.mids, prop.thicknesses):
                    mid_ref = model.materials[mid]
                    rho = mid_ref.rho
                    massi = area * (nsm_n + rho * t)
                    mid_to_mass[mid] += massi
                    pid_to_mid_to_rho_mass[pid][mid] += area * rho * t
                    pid_to_mid_to_nsm_mass[pid][mid] += area * nsm_n
                    pid_to_mid_to_mass[pid][mid] += massi
        else:
            raise NotImplementedError(elem)
    return dict(mid_to_mass), dict(pid_to_mid_to_mass), dict(pid_to_mid_to_rho_mass), dict(pid_to_mid_to_nsm_mass)

def get_property_mass_breakdown_table(model: BDF):
    pids_to_thickness = get_thickness_breakdown(model, property_ids=None, stop_if_no_thickness=False)
    pids_to_length = get_length_breakdown(model, property_ids=None, stop_if_no_length=False)
    pids_to_area = get_area_breakdown(model, property_ids=None, sum_bar_area=False, stop_if_no_area=False)
    pids_to_volume = get_volume_breakdown(model, property_ids=None, stop_if_no_volume=False)
    pids_to_mass, mass_type_to_mass = get_mass_breakdown(model, stop_if_no_mass=False, detailed=False)

    pids = set(list(pids_to_length) +
               list(pids_to_area) +
               list(pids_to_thickness) +
               list(pids_to_volume) +
               list(pids_to_mass))
    pids = list(pids)
    pids.sort()
    data = []
    for pid in pids:
        length = pids_to_length.get(pid, 0.)
        area = pids_to_area.get(pid, 0.)
        thickness = pids_to_thickness.get(pid, 0.)
        volume = pids_to_volume.get(pid, 0.)
        mass = pids_to_mass.get(pid, 0.)
        data.append((pid, length, area, thickness, volume, mass))
    return data

def get_length_breakdown(model: BDF, property_ids=None,
                         stop_if_no_length: bool=True):
    """
    Gets a breakdown of the length by property region.

    Returns
    -------
    pids_to_length : dict[int pid] : float length
        the pid to length dictionary

    TODO: What about CONRODs?

    """
    #skip_elements = [
        #'CTRIA3', 'CTRIA6', 'CTRIAR',
        #'CQUAD4', 'CQUAD8',
        #'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        #'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
        #'CBUSH', 'CBUSH1D', 'CBUSH2D',
        #'CRAC2D', 'CRAC3D',
        #'CBEAM3',
    #]
    skip_props = {
        'PSOLID', 'PLPLANE', 'PPLANE', 'PELAS',
        'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
        'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
        'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX', 'PLSOLID',
        'PCOMPS', 'PCOMPLS', 'PVISC',
        'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',

        # Nastran 95
        'PIHEX',

        # lines - should be included
        'PBEND', # 'PBEAM3',

        # acoustic
        'PACABS', 'PAABSF', 'PACBAR', 'PMIC',

        # welds - not sure
        'PWELD',
    }
    bar_properties = {'PBAR', 'PBARL', 'PBEAM', 'PBEAML',
                      'PROD', 'PTUBE', 'PBRSECT', 'PBMSECT', 'PBCOMP',
                      'PBEAM3'}
    pid_eids = model.get_element_ids_dict_with_pids(
        property_ids, stop_if_no_eids=stop_if_no_length,
        msg=' which is required by get_length_breakdown')
    pids_to_length = {}
    for pid, eids in pid_eids.items():
        prop = model.properties[pid]
        lengths = []
        if prop.type in skip_props:
            continue
        elif prop.type in bar_properties:
            #['CBAR', 'CBEAM', 'CROD', 'CTUBE']:
            # TODO: Do I need to consider the offset on length effects for a CBEAM?
            for eid in eids:
                elem = model.elements[eid]
                try:
                    lengths.append(elem.Length())
                except AttributeError:  # pragma: no cover
                    print(prop)
                    print(elem)
                    raise
        else:  # pragma: no cover
            print('prop\n%s' % prop)
            eid0 = eids[0]
            elem = model.elements[eid0]
            msg = str(prop) + str(elem)
            raise NotImplementedError(msg)

        if lengths:
            pids_to_length[pid] = sum(lengths)

    has_length = len(pids_to_length)

    if not has_length:
        msg = 'No elements with length were found'
        model.log.warning(msg)
        if stop_if_no_length:
            raise RuntimeError(msg)
    return pids_to_length

def get_area_breakdown(model: BDF,
                       property_ids=None,
                       stop_if_no_area: bool=True,
                       sum_bar_area: bool=True) -> dict[int, float]:
    """
    Gets a breakdown of the area by property region.

    Parameters
    ----------
    model: BDF
        a BDF object
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
        #'PBRSECT', 'PBCOMP', 'PBMSECT', 'PBEAM3', 'PBEND', 'PIHEX', 'PCOMPS',

    """
    skip_props = {
        'PSOLID', 'PLPLANE', 'PPLANE', 'PELAS',
        'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
        'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
        'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX', 'PLSOLID',
        'PCOMPS', 'PCOMPLS', 'PVISC', 'PBCOMP', 'PBEND',

        # Nastran 95
        'PIHEX',

        # lines - should be included
        'PBEND', # 'PBEAM3',

        # acoustic
        'PACABS', 'PAABSF', 'PACBAR', 'PMIC',

        # weld
        'PWELD',
    }
    bar_properties = {
        'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE', 'PBEAM3'}

    pid_eids = model.get_element_ids_dict_with_pids(
        property_ids, stop_if_no_eids=stop_if_no_area,
        msg=' which is required by get_area_breakdown')
    pids_to_area = {}
    for pid, eids in pid_eids.items():
        prop = model.properties[pid]
        areas = []
        if prop.type in {'PSHELL', 'PCOMP', 'PSHEAR', 'PCOMPG', }:
            for eid in eids:
                elem = model.elements[eid]
                if elem.type in {'CQUADX'}:
                    continue
                try:
                    areas.append(elem.Area())
                except AttributeError:  # pragma: no cover
                    print(prop)
                    print(elem)
                    raise
        elif prop.type in bar_properties:
            for eid in eids:
                elem = model.elements[eids[0]]
                area = elem.Area()
                if sum_bar_area:
                    neids = len(eids)
                    areas = [area * neids]
                else:
                    areas = [area]
        elif prop.type in skip_props:
            pass
        elif prop.type in {'PBRSECT', 'PBMSECT'}:
            model.log.warning('skipping:\n%s' % prop)
            continue
        else:  # pragma: no cover
            raise NotImplementedError(prop)
        if areas:
            pids_to_area[pid] = sum(areas)

    has_area = len(pids_to_area)
    if not has_area:
        msg = 'No elements with area were found'
        model.log.warning(msg)
        if stop_if_no_area:
            raise RuntimeError(msg)
    return pids_to_area


def get_thickness_breakdown(model: BDF,
                            property_ids=None,
                            stop_if_no_thickness: bool=False):
    skip_props = {
        'PBAR', 'PBARL', 'PROD', 'PSOLID', 'PBUSH', 'PBUSH1D',
        'PGAP', 'PRAC2D', 'PRAC3D', 'PWELD'}
    pid_eids = model.get_element_ids_dict_with_pids(
        property_ids, stop_if_no_eids=stop_if_no_thickness,
        msg=' which is required by get_thickness_breakdown')

    pids_to_thickness = {}
    for pid, eids in pid_eids.items():
        if len(eids) == 0:
            continue

        prop = model.properties[pid]
        thickness_ = []
        if prop.type in skip_props:
            continue
        elif prop.type in {'PSHELL', 'PCOMP', 'PCOMPG'}:
            #['CBAR', 'CBEAM', 'CROD', 'CTUBE']:
            # TODO: Do I need to consider the offset on length effects for a CBEAM?
            for eid in eids:
                elem = model.elements[eid]
                try:
                    thickness_.append(elem.Thickness())
                except AttributeError:  # pragma: no cover
                    print(prop)
                    print(elem)
                    raise
        else:  # pragma: no cover
            print(f'prop\n%s' % prop)
            print(f'eids = {eids}\n')
            eid0 = eids[0]
            elem = model.elements[eid0]
            msg = str(prop) + str(elem)
            raise NotImplementedError(msg)

        if thickness_:
            pids_to_thickness[pid] = sum(thickness_) / len(thickness_)

    has_thickness = len(pids_to_thickness)

    if not has_thickness:
        msg = 'No elements with thickness were found'
        model.log.warning(msg)
        if stop_if_no_thickness:
            raise RuntimeError(msg)
    return pids_to_thickness

def get_volume_breakdown(model: BDF, property_ids=None,
                         stop_if_no_volume: bool=True):
    """
    Gets a breakdown of the volume by property region.

    Parameters
    ----------
    model: BDF
        a BDF object
    property_ids : list[int] / int
        list of property ID
    stop_if_no_volume : bool; default=True
        prevents crashing if there are no elements

    TODO: What about CONRODs?
    #'PBRSECT',
    #'PBCOMP',
    #'PBMSECT',
    #'PBEAM3',
    #'PBEND',
    #'PIHEX',

    """
    pid_eids = model.get_element_ids_dict_with_pids(
        property_ids, stop_if_no_eids=stop_if_no_volume,
        msg=' which is required by get_volume_breakdown')

    no_volume = {
        'PLPLANE', 'PPLANE', 'PELAS',
        'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
        'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
        'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX',
        'PVISC', 'PBCOMP', 'PBEND',

        # lines - should be included
        'PBEND', 'PBEAM3',

        # acoustic
        'PACABS', 'PAABSF', 'PACBAR', 'PMIC',

        # weld
        'PWELD',
    }
    bar_properties = {
        'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE', # 'PBEAM3'
    }

    pids_to_volume = {}
    skipped_eid_pid = set()
    for pid, eids in pid_eids.items():
        prop = model.properties[pid]
        volumes = []
        if prop.type == 'PSHELL':
            # TODO: doesn't support PSHELL differential thicknesses
            thickness = prop.t
            areas = []
            for eid in eids:
                elem = model.elements[eid]
                if elem.type in ['CQUADX']:
                    continue
                areas.append(elem.Area())
            volumesi = [area * thickness for area in areas]
            volumes.extend(volumesi)
        elif prop.type in ['PCOMP', 'PCOMPG',]:
            areas = []
            for eid in eids:
                elem = model.elements[eid]
                areas.append(elem.Area())
            thickness = prop.Thickness()
            volumesi = [area * thickness for area in areas]
            volumes.extend(volumesi)
        elif prop.type in bar_properties:
            # TODO: Do I need to consider the offset on length effects for a CBEAM?
            lengths = []
            for eid in eids:
                elem = model.elements[eid]
                length = elem.Length()
                lengths.append(length)
            area = prop.Area()
            volumesi = [area * length for length in lengths]
            volumes.extend(volumesi)
        elif prop.type in ['PBEAM3']:
            for eid in eids:
                elem = model.elements[eid]
                volumei = elem.Volume()
                volumes.append(volumei)
        elif prop.type in ['PSOLID', 'PCOMPS', 'PCOMPLS', 'PLSOLID', 'PIHEX']:
            for eid in eids:
                elem = model.elements[eid]
                if elem.type in ['CTETRA', 'CPENTA', 'CHEXA']:
                    volumes.append(elem.Volume())
                else:
                    key = (elem.type, prop.type)
                    if key not in skipped_eid_pid:
                        skipped_eid_pid.add(key)
                        model.log.debug('skipping volume %s' % str(key))
        elif prop.type == 'PSHEAR':
            thickness = prop.t
            areas = []
            for eid in eids:
                elem = model.elements[eid]
                areas.append(elem.Area())
            volumesi = [area * thickness for area in areas]
            volumes.extend(volumesi)
        elif prop.type in no_volume:
            pass
        elif prop.type in ['PBRSECT', 'PBMSECT']:
            model.log.warning('skipping:\n%s' % prop)
            continue
        else:  # pragma: no cover
            raise NotImplementedError(prop)
        if volumes:
            pids_to_volume[pid] = sum(volumes)

    has_volume = len(pids_to_volume)
    if not has_volume:
        msg = 'No elements with volume were found'
        model.log.warning(msg)
        if stop_if_no_volume:
            raise RuntimeError(msg)
    return pids_to_volume

def get_mass_breakdown(model: BDF,
                       property_ids: list[int]=None,
                       stop_if_no_mass: bool=True,
                       detailed: bool=False) -> Any:
    """
    Gets a breakdown of the mass by property region.

    Parameters
    ----------
    model: BDF
        a BDF object
    property_ids : list[int] / int
        list of property ID
    stop_if_no_mass : bool; default=True
        prevents crashing if there are no elements
        setting this to False really doesn't make sense for non-DMIG models
    detailed : bool, optional, default : False
        Separates structural and nonstructural mass outputs.

    Returns
    -------
    pids_to_mass : dict {int : float, ...}
        Map from property id to mass (structural mass only if detailed is True).
    pids_to_mass_nonstructural : dict {int : float, ...}, optional
        Map from property id to nonstructural mass (only if detailed is True).
    mass_type_to_mass : dict {str : float, ...}
        Map from mass id to mass for mass elements.
        Used for CONM2s

    TODO: What about CONRODs, CONM2s
        #'PBCOMP', 'PBMSECT', 'PBEAM3', 'PBEND', 'PIHEX', 'PCOMPS',

    ..note:: WTMASS is not considered
    """
    pid_eids = model.get_element_ids_dict_with_pids(
        property_ids, stop_if_no_eids=False,
        msg=' which is required by get_mass_breakdown')

    mass_type_to_mass = {}
    pids_to_mass = {}
    pids_to_mass_nonstructural = {}
    skipped_eid_pid = set()
    for eid, elem in model.masses.items():
        if elem.type not in mass_type_to_mass:
            mass_type_to_mass[elem.type] = elem.Mass()
        else:
            mass_type_to_mass[elem.type] += elem.Mass()

    properties_to_skip = {
        'PLPLANE', 'PPLANE', 'PELAS',
        'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
        'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
        'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX',
        'PVISC',

        # lines - should be included
        'PBCOMP', 'PBEND', 'PBEAM3',

        # acoustic
        'PACABS', 'PAABSF', 'PACBAR', 'PMIC',
    }
    bar_properties = {'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE'}
    for pid, eids in pid_eids.items():
        prop = model.properties[pid]
        masses = []
        masses_nonstructural = []
        if prop.type == 'PSHELL':
            # TODO: doesn't support PSHELL differential thicknesses
            thickness = prop.t
            nsm = prop.nsm  # per area
            rho = prop.Rho()
            for eid in eids:
                elem = model.elements[eid]
                if elem.type == 'CQUADX':
                    continue
                area = elem.Area()
                if detailed:
                    masses.append(area * (rho * thickness))
                    masses_nonstructural.append(area * nsm)
                else:
                    masses.append(area * (rho * thickness + nsm))
        elif prop.type in {'PCOMP', 'PCOMPG'}:
            # TODO: does the PCOMP support differential thickness?
            #       I don't think so...
            for eid in eids:
                elem = model.elements[eid]
                if detailed:
                    mass, mass_nonstructural = elem.Mass_breakdown()
                    masses.append(mass)
                    masses_nonstructural.append(mass_nonstructural)
                else:
                    masses.append(elem.Mass())

        elif prop.type in bar_properties:
            nsm = prop.nsm # per unit length
            try:
                rho = prop.Rho()
            except AttributeError:
                print(prop)
                raise
            for eid in eids:
                elem = model.elements[eid]
                area = prop.Area()
                length = elem.Length()
                if detailed:
                    structural_mass_per_length = rho * area
                    masses.append(length * structural_mass_per_length)
                    masses_nonstructural.append(length * nsm)
                else:
                    mass_per_length = rho * area + nsm
                    masses.append(length * mass_per_length)
        #elif prop.type in ['PBEAM3']:
            #for eid in eids:
                #elem = model.elements[eid]
                #massi = elem.Mass()
                #volumes.append(massi)

        elif prop.type in {'PSOLID', 'PCOMPS', 'PCOMPLS', 'PLSOLID', 'PIHEX'}:
            rho = prop.Rho()
            for eid in eids:
                elem = model.elements[eid]
                if elem.type in {'CTETRA', 'CPENTA', 'CHEXA'}:
                    masses.append(rho * elem.Volume())
                else:
                    key = (elem.type, prop.type)
                    if key not in skipped_eid_pid:
                        skipped_eid_pid.add(key)
                        model.log.debug('skipping mass %s' % str(key))
        elif prop.type in properties_to_skip:
            pass
        elif prop.type == 'PSHEAR':
            thickness = prop.t
            nsm = prop.nsm # per area
            rho = prop.Rho()
            for eid in eids:
                elem = model.elements[eid]
                area = elem.Area()
                if detailed:
                    masses.append(area * (rho * thickness))
                    masses_nonstructural.append(area * nsm)
                else:
                    masses.append(area * (rho * thickness + nsm))
        elif prop.type in {'PBRSECT', 'PBMSECT', 'PWELD'}:
            model.log.warning('skipping:\n%s' % prop)
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
        model.log.warning(msg)
        if stop_if_no_mass:
            raise RuntimeError(msg)

    if detailed:
        return pids_to_mass, pids_to_mass_nonstructural, mass_type_to_mass
    return pids_to_mass, mass_type_to_mass
