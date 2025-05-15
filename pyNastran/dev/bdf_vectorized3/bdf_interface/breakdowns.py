from __future__ import annotations
from collections import defaultdict
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

NO_LENGTH = {
    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
    'CVISC', 'CGAP', 'CBUSH',
    'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUAD', 'CSHEAR',
    'CTRIA3', 'CTRIA6', 'CTRIAR',
    'CHEXA', 'CPENTA', 'CPYRAM', 'CTETRA',
    'CPENTCZ', 'CHEXCZ', 'CIFPENT', 'CIFHEX',
    'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4', 'CONM1', 'CONM2',
    'CQUADX4', 'CQUADX8', 'CQUADX', 'CTRIAX', 'CTRIAX6',
    'CPLSTS3', 'CPLSTS4', 'CPLSTS6', 'CPLSTS8',
    'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8',
    'GENEL', 'CAABSF', 'CHACAB', 'CHACBR',
    'CTRAX3', 'CTRAX4', 'CTRAX6',
}
NO_AREA = {
    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
    'CVISC', 'CGAP', 'CBUSH', 'CBUSH1D',
    'CHEXA', 'CPENTA', 'CPYRAM', 'CTETRA',
    'CPENTCZ', 'CHEXCZ', 'CIFPENT', 'CIFHEX',
    'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4', 'CONM1', 'CONM2',
    'CFAST',
    'GENEL',
    'CHACAB', 'CHACBR',
}
NO_VOLUME = {
    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
    'CVISC', 'CGAP', 'CBUSH', 'CBUSH1D',
    'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4', 'CONM1', 'CONM2',
    #'CQUADX4', 'CQUADX8', 'CQUADX', 'CTRIAX', 'CTRIAX6', # not really
    'CFAST',
    'GENEL', 'CAABSF',
}
NO_MASS = {
    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
    'CVISC', 'CGAP', 'CBUSH', 'CBUSH1D',
    #'CQUADX4', 'CQUADX8', 'CQUADX', 'CTRIAX', 'CTRIAX6', # not really
    'CFAST',
    'GENEL', 'CAABSF',
}

def get_mass_breakdown(model: BDF, stop_if_no_mass: bool=True) -> tuple[dict[int, float],
                                                                        dict[str, float]]:
    """
    Returns
    -------
    pids_to_mass : dict {int : float, ...}
        Map from property id to mass.
    mass_type_to_mass : dict {str : float, ...}
        Map from mass id to mass for mass elements.

    """
    pid_to_mass = defaultdict(float)
    mass_type_to_mass = {}
    element_cards = [element for element in model.element_cards
                     if element.n > 0]
    #mass_by_material_type = {}
    for card in element_cards:
        if card.n == 0 or card.type in NO_MASS:
            continue
        mass = card.mass()
        if card.type in {'CONROD', 'CTRIAX6'}:
            pid = -1 # np.full(card.element_id.shape, -1, dtype='int32')
            pid_to_mass[pid] += mass.sum()
            continue
        elif card.type in {'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4', 'CONM1', 'CONM2'}:
            mass_type_to_mass[card.type] = mass.sum()
            continue
        pid = card.property_id
        for pidi, massi in zip(pid, mass):
            pid_to_mass[pidi] += massi
    pid_to_mass = dict(pid_to_mass)
    if stop_if_no_mass and len(pid_to_mass) == 0 and len(pid_to_mass) == 0 and len(mass_type_to_mass) == 0:
        elements_str = ', '.join([str(element.type) for element in element_cards
                                if element.n > 0 ])
        raise RuntimeError(f'No elements with mass found\nelements = [{elements_str}]')

    #if detailed:
        #return pids_to_mass, pids_to_mass_nonstructural, mass_type_to_mass
    #return pids_to_mass, mass_type_to_mass
    return pid_to_mass, mass_type_to_mass

def get_length_breakdown(model: BDF, property_ids=None, stop_if_no_length: bool=True):
    """
    Gets a breakdown of the length by property region.

    Returns
    -------
    pids_to_length : dict[int pid] : float length
        the pid to length dictionary

    TODO: What about CONRODs?

    """
    pids_to_length = {}
    for element in model.element_cards:
        if element.n == 0 or element.type in NO_LENGTH:
            continue
        length = element.length()
        property_id = element.property_id
        for pid, lengthi in zip(property_id, length):
            pids_to_length[pid] += lengthi
    if stop_if_no_length and len(pids_to_length) == 0:
        raise RuntimeError('No elements with length found')
    return dict(pids_to_length)

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
        #'PBRSECT', 'PBCOMP', 'PBMSECT', 'PBEAM3', 'PBEND', 'PCOMPS',

    """
    pids_to_area = defaultdict(float)
    for element in model.element_cards:
        if element.n == 0 or element.type in NO_AREA:
            continue
        area = element.area()
        property_id = element.property_id
        for pid, areai in zip(property_id, area):
            pids_to_area[pid] += areai
    if stop_if_no_area and len(pids_to_area) == 0:
        raise RuntimeError('No elements with area found')
    return dict(pids_to_area)

    #skip_props = {
        #'PSOLID', 'PLPLANE', 'PPLANE', 'PELAS',
        #'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
        #'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
        #'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX', 'PLSOLID',
        #'PCOMPS', 'PCOMPLS', 'PVISC', 'PBCOMP', 'PBEND',

        ## lines - should be included
        #'PBEND', # 'PBEAM3',

        ## acoustic
        #'PACABS', 'PAABSF', 'PACBAR',
    #}
    #bar_properties = {
        #'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE', 'PBEAM3'}

    #pid_eids = model.get_element_ids_dict_with_pids(
        #property_ids, stop_if_no_eids=stop_if_no_area,
        #msg=' which is required by get_area_breakdown')
    #pids_to_area = {}
    #for pid, eids in pid_eids.items():
        #prop = model.properties[pid]
        #areas = []
        #if prop.type in {'PSHELL', 'PCOMP', 'PSHEAR', 'PCOMPG', }:
            #for eid in eids:
                #elem = model.elements[eid]
                #if elem.type in {'CQUADX'}:
                    #continue
                #try:
                    #areas.append(elem.Area())
                #except AttributeError:  # pragma: no cover
                    #print(prop)
                    #print(elem)
                    #raise
        #elif prop.type in bar_properties:
            #for eid in eids:
                #elem = model.elements[eids[0]]
                #area = elem.Area()
                #if sum_bar_area:
                    #neids = len(eids)
                    #areas = [area * neids]
                #else:
                    #areas = [area]
        #elif prop.type in skip_props:
            #pass
        #elif prop.type in {'PBRSECT', 'PBMSECT'}:
            #model.log.warning('skipping:\n%s' % prop)
            #continue
        #else:  # pragma: no cover
            #raise NotImplementedError(prop)
        #if areas:
            #pids_to_area[pid] = sum(areas)

    #has_area = len(pids_to_area)
    #if not has_area:
        #msg = 'No elements with area were found'
        #model.log.warning(msg)
        #if stop_if_no_area:
            #raise RuntimeError(msg)
    #return pids_to_area

def get_volume_breakdown(model: BDF,
                         property_ids=None,
                         stop_if_no_volume=True,
                         ) -> dict[int, float]:
    #pid_eids = model.get_element_ids_dict_with_pids(
        #property_ids, stop_if_no_eids=stop_if_no_volume,
        #msg=' which is required by get_volume_breakdown')

    #no_volume = {
        #'PLPLANE', 'PPLANE', 'PELAS',
        #'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
        #'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
        #'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX',
        #'PVISC', 'PBCOMP', 'PBEND',

        ## lines - should be included
        #'PBEND', 'PBEAM3',


        ## acoustic
        #'PACABS', 'PAABSF', 'PACBAR',
    #}
    bar_properties = {
        'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE', # 'PBEAM3'
    }

    pids_to_volume = defaultdict(float)
    #skipped_eid_pid = set()
    for element in model.element_cards:
        if element.n == 0 or element.type in NO_VOLUME:
            continue
        volume = element.volume()
        property_id = element.property_id
        for pid, volumei in zip(property_id, volume):
            pids_to_volume[pid] += volumei
    if stop_if_no_volume and len(pids_to_volume) == 0:
        raise RuntimeError('No elements with volume found')
    return dict(pids_to_volume)

    #for pid, eids in pid_eids.items():
        #prop = model.properties[pid]
        #volumes = []
        #if prop.type == 'PSHELL':
            ## TODO: doesn't support PSHELL differential thicknesses
            #thickness = prop.t
            #areas = []
            #for eid in eids:
                #elem = model.elements[eid]
                #if elem.type in ['CQUADX']:
                    #continue
                #areas.append(elem.Area())
            #volumesi = [area * thickness for area in areas]
            #volumes.extend(volumesi)
        #elif prop.type in ['PCOMP', 'PCOMPG',]:
            #areas = []
            #for eid in eids:
                #elem = model.elements[eid]
                #areas.append(elem.Area())
            #thickness = prop.Thickness()
            #volumesi = [area * thickness for area in areas]
            #volumes.extend(volumesi)
        #elif prop.type in bar_properties:
            ## TODO: Do I need to consider the offset on length effects for a CBEAM?
            #lengths = []
            #for eid in eids:
                #elem = model.elements[eid]
                #length = elem.Length()
                #lengths.append(length)
            #area = prop.Area()
            #volumesi = [area * length for length in lengths]
            #volumes.extend(volumesi)
        #elif prop.type in ['PBEAM3']:
            #for eid in eids:
                #elem = model.elements[eid]
                #volumei = elem.Volume()
                #volumes.append(volumei)
        #elif prop.type in ['PSOLID', 'PCOMPS', 'PCOMPLS', 'PLSOLID']:
            #for eid in eids:
                #elem = model.elements[eid]
                #if elem.type in ['CTETRA', 'CPENTA', 'CHEXA']:
                    #volumes.append(elem.Volume())
                #else:
                    #key = (elem.type, prop.type)
                    #if key not in skipped_eid_pid:
                        #skipped_eid_pid.add(key)
                        #model.log.debug('skipping volume %s' % str(key))
        #elif prop.type == 'PSHEAR':
            #thickness = prop.t
            #areas = []
            #for eid in eids:
                #elem = model.elements[eid]
                #areas.append(elem.Area())
            #volumesi = [area * thickness for area in areas]
            #volumes.extend(volumesi)
        #elif prop.type in no_volume:
            #pass
        #elif prop.type in ['PBRSECT', 'PBMSECT']:
            #model.log.warning('skipping:\n%s' % prop)
            #continue
        #else:  # pragma: no cover
            #raise NotImplementedError(prop)
        #if volumes:
            #pids_to_volume[pid] = sum(volumes)

    #has_volume = len(pids_to_volume)
    #if not has_volume:
        #msg = 'No elements with volume were found'
        #model.log.warning(msg)
        #if stop_if_no_volume:
            #raise RuntimeError(msg)
    #return pids_to_volume
