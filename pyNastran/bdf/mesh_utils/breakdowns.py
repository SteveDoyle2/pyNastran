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
import numpy as np
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def _build_nid_to_position(model: BDF) -> dict[int, np.ndarray]:
    """Builds a mapping of node id to global xyz position for all grid nodes."""
    nid_to_pos = {}
    for nid, node in model.nodes.items():
        nid_to_pos[nid] = node.get_position()
    return nid_to_pos


def _vectorized_quad4_areas(eids: list[int], model: BDF,
                            nid_to_pos: dict[int, np.ndarray]) -> np.ndarray:
    """Compute areas for a batch of CQUAD4 elements using vectorized numpy."""
    n = len(eids)
    xyz = np.empty((n, 4, 3), dtype='float64')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        nodes = elem.nodes_ref
        xyz[i, 0] = nid_to_pos[nodes[0].nid]
        xyz[i, 1] = nid_to_pos[nodes[1].nid]
        xyz[i, 2] = nid_to_pos[nodes[2].nid]
        xyz[i, 3] = nid_to_pos[nodes[3].nid]
    # area = 0.5 * |diag1 x diag2|
    diag1 = xyz[:, 2] - xyz[:, 0]  # n3 - n1
    diag2 = xyz[:, 3] - xyz[:, 1]  # n4 - n2
    crosses = np.cross(diag1, diag2)
    areas = 0.5 * np.linalg.norm(crosses, axis=1)
    return areas


def _vectorized_tri3_areas(eids: list[int], model: BDF,
                           nid_to_pos: dict[int, np.ndarray]) -> np.ndarray:
    """Compute areas for a batch of CTRIA3 elements using vectorized numpy."""
    n = len(eids)
    xyz = np.empty((n, 3, 3), dtype='float64')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        nodes = elem.nodes_ref
        xyz[i, 0] = nid_to_pos[nodes[0].nid]
        xyz[i, 1] = nid_to_pos[nodes[1].nid]
        xyz[i, 2] = nid_to_pos[nodes[2].nid]
    # area = 0.5 * |(n1-n2) x (n1-n3)|
    a = xyz[:, 0] - xyz[:, 1]
    b = xyz[:, 0] - xyz[:, 2]
    crosses = np.cross(a, b)
    areas = 0.5 * np.linalg.norm(crosses, axis=1)
    return areas


def _vectorized_chexa8_volumes(eids: list[int], model: BDF,
                               nid_to_pos: dict[int, np.ndarray]) -> np.ndarray:
    """Compute volumes for a batch of CHEXA8 elements using vectorized numpy.

    Uses the exact Grandy (1997) decomposition into 3 scalar triple products:
        V = [det3(x6-x0, x1-x0, x2-x5) +
             det3(x6-x0, x4-x0, x5-x7) +
             det3(x6-x0, x3-x0, x7-x2)] / 6
    where det3(a,b,c) = a . (b x c) and x0..x7 are the 8 hex nodes (0-indexed).
    Reference: https://www.osti.gov/servlets/purl/632793/

    Note: ~15% slower than (area_top+area_bot)/2*height per-element, but faster
    when vectorized over batches. Do not use for scalar CHEXA8.Volume().
    """
    n = len(eids)
    xyz = np.empty((n, 8, 3), dtype='float64')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        nodes = elem.nodes_ref
        for j in range(8):
            xyz[i, j] = nid_to_pos[nodes[j].nid]

    x6_x0 = xyz[:, 6] - xyz[:, 0]

    # term 1: det3(x6-x0, x1-x0, x2-x5)
    b1 = xyz[:, 1] - xyz[:, 0]
    c1 = xyz[:, 2] - xyz[:, 5]
    t1 = np.einsum('ij,ij->i', x6_x0, np.cross(b1, c1))

    # term 2: det3(x6-x0, x4-x0, x5-x7)
    b2 = xyz[:, 4] - xyz[:, 0]
    c2 = xyz[:, 5] - xyz[:, 7]
    t2 = np.einsum('ij,ij->i', x6_x0, np.cross(b2, c2))

    # term 3: det3(x6-x0, x3-x0, x7-x2)
    b3 = xyz[:, 3] - xyz[:, 0]
    c3 = xyz[:, 7] - xyz[:, 2]
    t3 = np.einsum('ij,ij->i', x6_x0, np.cross(b3, c3))

    volumes = (t1 + t2 + t3) / 6.0
    return volumes


def _vectorized_ctetra4_volumes(eids: list[int], model: BDF,
                                nid_to_pos: dict[int, np.ndarray]) -> np.ndarray:
    """Compute volumes for a batch of CTETRA4 elements using vectorized numpy.

    V = -(n1-n4) . ((n2-n4) x (n3-n4)) / 6
    """
    n = len(eids)
    xyz = np.empty((n, 4, 3), dtype='float64')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        nodes = elem.nodes_ref
        xyz[i, 0] = nid_to_pos[nodes[0].nid]
        xyz[i, 1] = nid_to_pos[nodes[1].nid]
        xyz[i, 2] = nid_to_pos[nodes[2].nid]
        xyz[i, 3] = nid_to_pos[nodes[3].nid]
    a = xyz[:, 0] - xyz[:, 3]  # n1 - n4
    b = xyz[:, 1] - xyz[:, 3]  # n2 - n4
    c = xyz[:, 2] - xyz[:, 3]  # n3 - n4
    # dot product per row: sum(a * cross(b, c), axis=1)
    volumes = -np.einsum('ij,ij->i', a, np.cross(b, c)) / 6.0
    return volumes


def _vectorized_bar_lengths(eids: list[int], model: BDF,
                            nid_to_pos: dict[int, np.ndarray]) -> np.ndarray:
    """Compute lengths for a batch of CBAR/CBEAM/CROD/CTUBE elements."""
    n = len(eids)
    p1 = np.empty((n, 3), dtype='float64')
    p2 = np.empty((n, 3), dtype='float64')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        p1[i] = nid_to_pos[elem.nodes_ref[0].nid]
        p2[i] = nid_to_pos[elem.nodes_ref[1].nid]
    lengths = np.linalg.norm(p2 - p1, axis=1)
    return lengths


def _split_eids_by_type(eids: list[int], model: BDF) -> dict[str, list[int]]:
    """Group element IDs by element type string."""
    type_to_eids: dict[str, list[int]] = {}
    for eid in eids:
        etype = model.elements[eid].type
        if etype not in type_to_eids:
            type_to_eids[etype] = []
        type_to_eids[etype].append(eid)
    return type_to_eids


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
        elif elem_type in {'CONROD'}:  # references a MAT directly
            mid = elem.Mid()
            length = elem.Length()
            massi = elem.Mass()
            nsm = elem.Nsm()
            area = elem.Area()
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
                rho = elem.Rho()
                t = elem.Thickness()
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
        length: float = pids_to_length.get(pid, 0.)
        area: float = pids_to_area.get(pid, 0.)
        thickness: float = pids_to_thickness.get(pid, 0.)
        volume: float = pids_to_volume.get(pid, 0.)
        mass: float = pids_to_mass.get(pid, 0.)
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

        # lines - should be included
        'PBEND',  # 'PBEAM3',

        # acoustic
        'PACABS', 'PAABSF', 'PACBAR', 'PMIC',

        # welds - not sure
        'PWELD',
        # sometimes included here
        'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',
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
        if prop.type in skip_props:
            continue
        elif prop.type in bar_properties:
            #['CBAR', 'CBEAM', 'CROD', 'CTUBE']
            # TODO: Do I need to consider the offset on length effects for a CBEAM?
            total_length = 0.0
            for eid in eids:
                total_length += model.elements[eid].Length()
            if total_length > 0.0:
                pids_to_length[pid] = total_length
        else:  # pragma: no cover
            print('prop\n%s' % prop)
            eid0 = eids[0]
            elem = model.elements[eid0]
            msg = str(prop) + str(elem)
            raise NotImplementedError(msg)

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
        #'PBRSECT', 'PBCOMP', 'PBMSECT', 'PBEAM3', 'PBEND', 'PCOMPS',

    """
    skip_props = {
        'PSOLID', 'PLPLANE', 'PPLANE', 'PELAS',
        'PDAMP', 'PBUSH', 'PBUSH1D', 'PBUSH2D',
        'PELAST', 'PDAMPT', 'PBUSHT', 'PDAMP5',
        'PFAST', 'PGAP', 'PRAC2D', 'PRAC3D', 'PCONEAX', 'PLSOLID',
        'PCOMPS', 'PCOMPLS', 'PVISC', 'PBCOMP', 'PBEND',

        # lines - should be included
        'PBEND',  # 'PBEAM3',

        # acoustic
        'PACABS', 'PAABSF', 'PACBAR', 'PMIC',

        # weld
        'PWELD',
        # sometimes included here
        'PMASS',
    }
    bar_properties = {
        'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE', 'PBEAM3'}

    pid_eids = model.get_element_ids_dict_with_pids(
        property_ids, stop_if_no_eids=stop_if_no_area,
        msg=' which is required by get_area_breakdown')

    nid_to_pos = _build_nid_to_position(model)
    pids_to_area = {}
    for pid, eids in pid_eids.items():
        prop = model.properties[pid]
        if prop.type in {'PSHELL', 'PCOMP', 'PSHEAR', 'PCOMPG'}:
            type_to_eids = _split_eids_by_type(eids, model)
            total_area = 0.0
            for etype, typed_eids in type_to_eids.items():
                if etype == 'CQUAD4':
                    total_area += _vectorized_quad4_areas(typed_eids, model, nid_to_pos).sum()
                elif etype == 'CTRIA3':
                    total_area += _vectorized_tri3_areas(typed_eids, model, nid_to_pos).sum()
                elif etype == 'CQUADX':
                    continue
                else:
                    for eid in typed_eids:
                    try:
                        total_area += model.elements[eid].Area()
                    except AttributeError:  # pragma: no cover
                        print(prop)
                        print(elem)
                        raise
            if total_area > 0.0:
                pids_to_area[pid] = total_area
        elif prop.type in bar_properties:
            elem = model.elements[eids[0]]
            area = elem.Area()
            if sum_bar_area:
                neids = len(eids)
                pids_to_area[pid] = area * neids
            else:
                pids_to_area[pid] = area
        elif prop.type in skip_props:
            pass
        elif prop.type in {'PBRSECT', 'PBMSECT'}:
            model.log.warning('skipping:\n%s' % prop)
            continue
        else:  # pragma: no cover
            raise NotImplementedError(prop)

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
        'PGAP', 'PRAC2D', 'PRAC3D', 'PWELD',
        # sometimes included here
        'PMASS',
    }
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
        # sometimes included here
        'PMASS',
    }
    bar_properties = {
        'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PTUBE',  # 'PBEAM3'
    }

    nid_to_pos = _build_nid_to_position(model)
    pids_to_volume = {}
    skipped_eid_pid = set()
    for pid, eids in pid_eids.items():
        prop = model.properties[pid]
        if prop.type == 'PSHELL':
            # TODO: does this support PSHELL differential thicknesses?
            thickness = prop.t
            type_to_eids = _split_eids_by_type(eids, model)
            total_area = 0.0
            for etype, typed_eids in type_to_eids.items():
                if etype == 'CQUAD4':
                    total_area += _vectorized_quad4_areas(typed_eids, model, nid_to_pos).sum()
                elif etype == 'CTRIA3':
                    total_area += _vectorized_tri3_areas(typed_eids, model, nid_to_pos).sum()
                elif etype == 'CQUADX':
                    continue
                else:
                    for eid in typed_eids:
                        total_area += model.elements[eid].Area()
            if total_area > 0.0:
                pids_to_volume[pid] = total_area * thickness
        elif prop.type in ['PCOMP', 'PCOMPG']:
            thickness = prop.Thickness()
            type_to_eids = _split_eids_by_type(eids, model)
            total_area = 0.0
            for etype, typed_eids in type_to_eids.items():
                if etype == 'CQUAD4':
                    total_area += _vectorized_quad4_areas(typed_eids, model, nid_to_pos).sum()
                elif etype == 'CTRIA3':
                    total_area += _vectorized_tri3_areas(typed_eids, model, nid_to_pos).sum()
                else:
                    for eid in typed_eids:
                        total_area += model.elements[eid].Area()
            if total_area > 0.0:
                pids_to_volume[pid] = total_area * thickness
        elif prop.type in bar_properties:
            # TODO: Do I need to consider the offset on length effects for a CBEAM?
            #       No? Because the end points should already be offset...they are not
            type_to_eids = _split_eids_by_type(eids, model)
            total_length = 0.0
            for etype, typed_eids in type_to_eids.items():
                if etype in ('CBAR', 'CBEAM', 'CROD', 'CTUBE'):
                    total_length += _vectorized_bar_lengths(typed_eids, model, nid_to_pos).sum()
                else:
                    for eid in typed_eids:
                        total_length += model.elements[eid].Length()
            area = prop.Area()
            if total_length > 0.0:
                pids_to_volume[pid] = area * total_length
        elif prop.type in ['PBEAM3']:
            total_volume = 0.0
            for eid in eids:
                total_volume += model.elements[eid].Volume()
            if total_volume > 0.0:
                pids_to_volume[pid] = total_volume
        elif prop.type in ['PSOLID', 'PCOMPS', 'PCOMPLS', 'PLSOLID']:
            type_to_eids = _split_eids_by_type(eids, model)
            total_volume = 0.0
            for etype, typed_eids in type_to_eids.items():
                if etype == 'CHEXA':
                    total_volume += _vectorized_chexa8_volumes(typed_eids, model, nid_to_pos).sum()
                elif etype == 'CTETRA':
                    total_volume += _vectorized_ctetra4_volumes(typed_eids, model, nid_to_pos).sum()
                elif etype in ('CPENTA', 'CPYRAM'):
                    for eid in typed_eids:
                        total_volume += model.elements[eid].Volume()
                else:
                    key = (etype, prop.type)
                    if key not in skipped_eid_pid:
                        skipped_eid_pid.add(key)
                        model.log.debug('skipping volume %s' % str(key))
            if total_volume > 0.0:
                pids_to_volume[pid] = total_volume
        elif prop.type == 'PSHEAR':
            thickness = prop.t
            type_to_eids = _split_eids_by_type(eids, model)
            total_area = 0.0
            for etype, typed_eids in type_to_eids.items():
                if etype == 'CQUAD4':
                    total_area += _vectorized_quad4_areas(typed_eids, model, nid_to_pos).sum()
                else:
                    for eid in typed_eids:
                        total_area += model.elements[eid].Area()
            if total_area > 0.0:
                pids_to_volume[pid] = total_area * thickness
        elif prop.type in no_volume:
            pass
        elif prop.type in ['PBRSECT', 'PBMSECT']:
            model.log.warning('skipping:\n%s' % prop)
            continue
        else:  # pragma: no cover
            raise NotImplementedError(prop)

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
        #'PBCOMP', 'PBMSECT', 'PBEAM3', 'PBEND', 'PCOMPS',

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

    nid_to_pos = _build_nid_to_position(model)
    for pid, eids in pid_eids.items():
        prop = model.properties[pid]
        if prop.type == 'PSHELL':
            # TODO: doesn't support PSHELL differential thicknesses
            thickness = prop.t
            nsm = prop.nsm  # per area
            rho = prop.Rho()
            type_to_eids = _split_eids_by_type(eids, model)
            total_area = 0.0
            for etype, typed_eids in type_to_eids.items():
                if etype == 'CQUAD4':
                    total_area += _vectorized_quad4_areas(typed_eids, model, nid_to_pos).sum()
                elif etype == 'CTRIA3':
                    total_area += _vectorized_tri3_areas(typed_eids, model, nid_to_pos).sum()
                elif etype == 'CQUADX':
                    continue
                else:
                    for eid in typed_eids:
                        total_area += model.elements[eid].Area()
            if total_area > 0.0:
                if detailed:
                    pids_to_mass[pid] = total_area * (rho * thickness)
                    pids_to_mass_nonstructural[pid] = total_area * nsm
                else:
                    pids_to_mass[pid] = total_area * (rho * thickness + nsm)

        elif prop.type in {'PCOMP', 'PCOMPG'}:
            # PCOMP mass_per_area is constant for all elements under same pid
            mpa = prop.MassPerArea()
            if detailed:
                mpa_struct = prop.MassPerArea_structure()
                nsm = prop.nsm
            type_to_eids = _split_eids_by_type(eids, model)
            total_area = 0.0
            for etype, typed_eids in type_to_eids.items():
                if etype == 'CQUAD4':
                    total_area += _vectorized_quad4_areas(typed_eids, model, nid_to_pos).sum()
                elif etype == 'CTRIA3':
                    total_area += _vectorized_tri3_areas(typed_eids, model, nid_to_pos).sum()
                else:
                    for eid in typed_eids:
                        total_area += model.elements[eid].Area()
            if total_area > 0.0:
                if detailed:
                    pids_to_mass[pid] = total_area * mpa_struct
                    pids_to_mass_nonstructural[pid] = total_area * nsm
                else:
                    pids_to_mass[pid] = total_area * mpa

        elif prop.type in bar_properties:
            nsm = prop.Nsm()  # per unit length (scalar)
            rho = prop.Rho()
            area = prop.Area()
            type_to_eids = _split_eids_by_type(eids, model)
            total_length = 0.0
            for etype, typed_eids in type_to_eids.items():
                if etype in ('CBAR', 'CBEAM', 'CROD', 'CTUBE'):
                    total_length += _vectorized_bar_lengths(typed_eids, model, nid_to_pos).sum()
                else:
                    for eid in typed_eids:
                        total_length += model.elements[eid].Length()
            if total_length > 0.0:
                if detailed:
                    pids_to_mass[pid] = total_length * (rho * area)
                    pids_to_mass_nonstructural[pid] = total_length * nsm
                else:
                    pids_to_mass[pid] = total_length * (rho * area + nsm)

        elif prop.type in {'PSOLID', 'PCOMPS', 'PCOMPLS', 'PLSOLID'}:
            rho = prop.Rho()
            type_to_eids = _split_eids_by_type(eids, model)
            total_volume = 0.0
            for etype, typed_eids in type_to_eids.items():
                if etype == 'CHEXA':
                    total_volume += _vectorized_chexa8_volumes(typed_eids, model, nid_to_pos).sum()
                elif etype == 'CTETRA':
                    total_volume += _vectorized_ctetra4_volumes(typed_eids, model, nid_to_pos).sum()
                elif etype in ('CPENTA', 'CPYRAM'):
                    for eid in typed_eids:
                        total_volume += model.elements[eid].Volume()
                else:
                    key = (etype, prop.type)
                    if key not in skipped_eid_pid:
                        skipped_eid_pid.add(key)
                        model.log.debug('skipping mass %s' % str(key))
            if total_volume > 0.0:
                pids_to_mass[pid] = rho * total_volume

        elif prop.type in properties_to_skip:
            pass
        elif prop.type == 'PSHEAR':
            thickness = prop.t
            nsm = prop.nsm  # per area
            rho = prop.Rho()
            type_to_eids = _split_eids_by_type(eids, model)
            total_area = 0.0
            for etype, typed_eids in type_to_eids.items():
                if etype == 'CQUAD4':
                    total_area += _vectorized_quad4_areas(typed_eids, model, nid_to_pos).sum()
                else:
                    for eid in typed_eids:
                        total_area += model.elements[eid].Area()
            if total_area > 0.0:
                if detailed:
                    pids_to_mass[pid] = total_area * (rho * thickness)
                    pids_to_mass_nonstructural[pid] = total_area * nsm
                else:
                    pids_to_mass[pid] = total_area * (rho * thickness + nsm)
        elif prop.type in {'PBRSECT', 'PBMSECT', 'PWELD'}:
            model.log.warning('skipping:\n%s' % prop)
            continue
        else:  # pragma: no cover
            raise NotImplementedError(prop)

    has_mass = len(mass_type_to_mass) > 0 or len(pids_to_mass) > 0
    if not has_mass:
        msg = 'No elements with mass were found'
        model.log.warning(msg)
        if stop_if_no_mass:
            raise RuntimeError(msg)

    if detailed:
        return pids_to_mass, pids_to_mass_nonstructural, mass_type_to_mass
    return pids_to_mass, mass_type_to_mass
