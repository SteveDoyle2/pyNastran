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


def get_length_breakdown(model, property_ids=None, stop_if_no_length=True):
    """
    Gets a breakdown of the length by property region.

    Returns
    -------
    pids_to_length : Dict[int pid] : float length
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
        'PCOMPS', 'PVISC',
        'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',

        # Nastran 95
        'PIHEX',

        # lines - should be included
        'PBEND', # 'PBEAM3',

        # acoustic
        'PACABS', 'PAABSF', 'PACBAR',
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

def get_area_breakdown(model, property_ids=None, stop_if_no_area=True, sum_bar_area=True):
    """
    Gets a breakdown of the area by property region.

    Parameters
    ----------
    property_ids : List[int] / int
        list of property ID
    stop_if_no_area : bool; default=True
        prevents crashing if there are no elements
    sum_bar_area : bool; default=True
        sum the areas for CBAR/CBEAM/CROD/CONROD/CTUBE elements
        True : get the area of the model by property id
        False : only get the cross sectional properties

    Returns
    -------
    pids_to_area : Dict[int pid] : float area
        the pid to area dictionary

    TODO: What about CONRODs?
    #'PBRSECT',
    #'PBCOMP',
    #'PBMSECT',
    #'PBEAM3',
    #'PBEND',
    #'PIHEX',
    #'PCOMPS',

    """
    skip_props = {
        'PSOLID', 'PLPLANE', 'PPLANE', 'PELAS',
        'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
        'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
        'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX', 'PLSOLID',
        'PCOMPS', 'PVISC', 'PBCOMP', 'PBEND',

        # Nastran 95
        'PIHEX',

        # lines - should be included
        'PBEND', # 'PBEAM3',

        # acoustic
        'PACABS', 'PAABSF', 'PACBAR',
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
        if prop.type in ['PSHELL', 'PCOMP', 'PSHEAR', 'PCOMPG', ]:
            for eid in eids:
                elem = model.elements[eid]
                if elem.type in ['CQUADX']:
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
        elif prop.type in ['PBRSECT', 'PBMSECT']:
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

def get_volume_breakdown(model, property_ids=None, stop_if_no_volume=True):
    """
    Gets a breakdown of the volume by property region.

    Parameters
    ----------
    property_ids : List[int] / int
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
        'PACABS', 'PAABSF', 'PACBAR',
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
        elif prop.type in ['PSOLID', 'PCOMPS', 'PLSOLID', 'PIHEX']:
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

def get_mass_breakdown(model, property_ids=None, stop_if_no_mass=True, detailed=False):
    """
    Gets a breakdown of the mass by property region.

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
    -------
    pids_to_mass : dict {int : float, ...}
        Map from property id to mass (structural mass only if detailed is True).
    pids_to_mass_nonstructural : dict {int : float, ...}, optional
        Map from property id to nonstructural mass (only if detailed is True).
    mass_type_to_mass : dict {str : float, ...}
        Map from mass id to mass for mass elements.
        Used for CONM2s

    TODO: What about CONRODs?
    #'PBCOMP',
    #'PBMSECT',
    #'PBEAM3',
    #'PBEND',
    #'PIHEX',
    #'PCOMPS',

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
        'PACABS', 'PAABSF', 'PACBAR',
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
                    masses.append(length * (rho * area))
                    masses_nonstructural.append(length * nsm)
                else:
                    masses.append(length * (rho * area + nsm))
        #elif prop.type in ['PBEAM3']:
            #for eid in eids:
                #elem = model.elements[eid]
                #massi = elem.Mass()
                #volumes.append(massi)

        elif prop.type in ['PSOLID', 'PCOMPS', 'PLSOLID', 'PIHEX']:
            rho = prop.Rho()
            for eid in eids:
                elem = model.elements[eid]
                if elem.type in ['CTETRA', 'CPENTA', 'CHEXA']:
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
        elif prop.type in ['PBRSECT', 'PBMSECT']:
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
