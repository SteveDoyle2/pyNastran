"""
defines methods to access force/moment/pressure/temperature data:
  - get_forces_moments_array(model, p0, load_case_id,
                            eid_map, nnodes, normals, dependents_nodes,
                            nid_map=None, include_grav=False)
  - get_pressure_array(model, load_case, eids, stop_on_failure=True)
  - get_temperatures_array(model: BDF, load_case_id, nid_map=None, dtype='float32')
  - get_load_arrays(model, subcase_id, nid_map, eid_map, node_ids, normals)

"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.loads.static_loads import update_pload4_vector
from pyNastran.bdf.mesh_utils.loads import _mean_pressure_on_pload4

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def get_forces_moments_array(model: BDF, p0, load_case_id : int,
                             eid_map, nnodes, normals, dependents_nodes,
                             nid_map=None, include_grav=False):
    """
    Gets the forces/moments on the nodes.

    Parameters
    ----------
    p0 : (3, ) float ndarray
        the reference location
    load_case_id : int
        the load id
    nid_map : ???
        ???
    eid_map : Dict[int eid : int index]
        ???
    nnodes : ???
        the number of nodes in nid_map
    normals : (nelements, 3) float ndarray
        the normal vectors for the shells
        what about solids???
    dependents_nodes : ???
        ???
    include_grav : bool; default=False
        is the mass of the elements considered; unused

    Returns
    -------
    temperature_data : tuple(temperature_key, temperatures)
        temperature_key : str
            One of the following:
              TEMPERATURE(MATERIAL)
              TEMPERATURE(INITIAL)
              TEMPERATURE(LOAD)
              TEMPERATURE(BOTH)
        temperatures : (nnodes, 1) float ndarray
            the temperatures
    load_data : tuple(centroidal_pressures, forces, spcd)
        centroidal_pressures : (nelements, 1) float ndarray
            the pressure
        forces : (nnodes, 3) float ndarray
            the pressure
        spcd : (nnodes, 3) float ndarray
            the SPCD load application

    Considers
    FORCE
    PLOAD2 - CTRIA3, CQUAD4, CSHEAR
    PLOAD4 - CTRIA3, CTRIA6, CTRIAR
             CQUAD4, CQUAD8, CQUAD, CQUADR, CSHEAR
             CTETRA, CPENTA, CHEXA
    SPCD

    """
    if nid_map is None:
        nid_map = model.nid_map
    if not any(['FORCE' in model.card_count,
                'PLOAD' in model.card_count, 'PLOAD2' in model.card_count,
                'PLOAD4' in model.card_count, 'SPCD' in model.card_count,
                'SLOAD' in model.card_count]):
        return None, None, None
    assert len(nid_map) == nnodes, f'len(nid_map)={len(nid_map)} nnodes={nnodes}'

    loads, scale_factors = model.get_reduced_loads(
        load_case_id, skip_scale_factor0=True)[:2]

    #eids = sorted(model.elements.keys())
    centroidal_pressures = np.zeros(len(model.elements), dtype='float32')
    nodal_pressures = np.zeros(len(model.node_ids), dtype='float32')

    forces = np.zeros((nnodes, 3), dtype='float32')
    spcd = np.zeros((nnodes, 3), dtype='float32')
    # loop thru scaled loads and plot the pressure
    cards_ignored = set()

    assert normals is not None
    fail_nids = set()
    fail_count = 0
    fail_count_max = 3
    loads_to_skip = ['MOMENT', 'MOMENT1', 'MOMENT2', 'FORCE1', 'TEMP']
    nodes = model.nodes
    for load, scale in zip(loads, scale_factors):
        load_type = load.type
        if load_type in loads_to_skip:
            pass
        elif load_type == 'FORCE':
            scale2 = load.mag * scale  # does this need a magnitude?
            nid = load.node
            if nid in dependents_nodes:
                fail_nids.add(nid)
                fail_count += 1
                if fail_count < fail_count_max:
                    print('    nid=%s is a dependent node and has a FORCE applied\n%s' % (
                        nid, str(load)))
            forces[nid_map[nid]] += load.xyz * scale2

        elif load_type == 'PLOAD':
            pressure = load.pressure * scale
            nnodes = len(load.nodes)
            if nnodes == 4:
                n1, n2, n3, n4 = load.nodes
                xyz1 = nodes[n1].get_position()
                xyz2 = nodes[n2].get_position()
                xyz3 = nodes[n3].get_position()
                xyz4 = nodes[n4].get_position()
                normal_area = np.cross(xyz3 - xyz1, xyz4 - xyz2)  # TODO: not validated
            elif nnodes == 3:
                n1, n2, n3 = load.nodes
                xyz1 = nodes[n1].get_position()
                xyz2 = nodes[n2].get_position()
                xyz3 = nodes[n3].get_position()
                normal_area = np.cross(xyz2 - xyz1, xyz3 - xyz1)  # TODO: not validated
            else:
                model.log.debug('    case=%s nnodes=%r loadtype=%r not supported' % (
                    load_case_id, nnodes, load.type))
                continue
            forcei = pressure * normal_area / nnodes
            for nid in load.nodes:
                forces[nid_map[nid]] += forcei

        elif load_type == 'PLOAD2':
            pressure = load.pressure * scale  # there are 4 pressures, but we assume p0
            for eid in load.eids:
                elem = model.elements[eid]
                if elem.type in ['CTRIA3',
                                 'CQUAD4', 'CSHEAR']:
                    node_ids = elem.node_ids
                    nnodes = len(node_ids)
                    ie = eid_map[eid]
                    normal = normals[ie, :]

                    area = elem.Area()
                    forcei = pressure * normal * area / nnodes
                    # r = elem.Centroid() - p0
                    # m = cross(r, f)
                    for nid in node_ids:
                        if nid in dependents_nodes:
                            fail_nids.add(nid)
                            fail_count += 1
                            if fail_count < fail_count_max:
                                print(f'    nid={nid} is a dependent node and has a '
                                      f'PLOAD2 applied\n{load}')
                        forces[nid_map[nid]] += forcei
                    forces += forcei
                    # F += f
                    # M += m
                else:
                    model.log.debug('    case=%s etype=%r loadtype=%r not supported' % (
                        load_case_id, elem.type, load.type))

        elif load_type == 'PLOAD4':
            # multiple elements
            eids_missing = []
            for elem in load.eids_ref:
                if isinstance(elem, integer_types):
                    # Nastran is NOT OK with missing element ids
                    eids_missing.append(elem)
                    continue
                ie = eid_map[elem.eid]
                normal = normals[ie, :]
                # pressures[eids.index(elem.eid)] += p
                if elem.type in ['CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR']:
                    area = elem.get_area()
                    elem_node_ids = elem.node_ids
                    nface = len(elem_node_ids)

                    if load.surf_or_line == 'SURF':
                        cid = load.Cid()
                        normal = update_pload4_vector(load, normal, cid)
                    else:
                        msg = (f'surf_or_line={load.surf_or_line!r} on PLOAD4 is not supported\n'
                               f'{load}')
                        model.log.debug(msg)
                        continue

                    pressures = load.pressures[:nface]
                    pressure = _mean_pressure_on_pload4(pressures, load, elem)

                    forcei = pressure * area * normal / nface
                    for nid in elem_node_ids:
                        if nid in dependents_nodes:
                            fail_nids.add(nid)
                            fail_count += 1
                            if fail_count < fail_count_max:
                                print('    nid=%s is a dependent node and has a'
                                      ' PLOAD4 applied\n%s' % (nid, str(load)))
                        #forces[nids.index(nid)] += F
                        i = nid_map[nid]
                        try:
                            forces[i, :] += forcei
                        except IndexError:
                            print(f'i = {i}')
                            print('normals.shape = %s' %  str(normals.shape))
                            print('forces.shape = %s' % str(forces.shape))
                            print('normal = ', normal)
                            print('forces[i, :] = ', forces[i, :])
                            raise
                    #nface = 3
                elif elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
                    area = elem.get_area()
                    elem_node_ids = elem.node_ids
                    nface = len(elem_node_ids)

                    if load.surf_or_line == 'SURF':
                        cid = load.Cid()
                        if cid in [0, None] and np.abs(load.nvector).max() == 0.0:
                            # element surface normal
                            pass
                        else:
                            if np.linalg.norm(load.nvector) != 0.0 and cid in [0, None]:
                                normal = load.nvector / np.linalg.norm(load.nvector)
                            else:
                                raise NotImplementedError('cid=%r nvector=%s on a PLOAD4 is not supported\n%s' % (
                                    cid, load.nvector, str(load)))
                    else:  # pragma: no cover
                        msg = (f'surf_or_line={load.surf_or_line} on PLOAD4 is not supported\n'
                               f'{load}')
                        model.log.debug(msg)
                        continue

                    pressures = load.pressures[:nface]
                    pressure = _mean_pressure_on_pload4(pressures, load, elem)

                    forcei = pressure * area * normal / nface

                    for nid in elem_node_ids:
                        if nid in dependents_nodes:
                            fail_nids.add(nid)
                            fail_count += 1
                            if fail_count < fail_count_max:
                                print('    nid=%s is a dependent node and has a'
                                      ' PLOAD4 applied\n%s' % (nid, str(load)))
                        #forces[nids.index(nid)] += F
                        i = nid_map[nid]
                        try:
                            forces[i, :] += forcei
                        except IndexError:
                            print('i = %s' % i)
                            print('normals.shape = %s' %  str(normals.shape))
                            print('forces.shape = %s' % str(forces.shape))
                            print('normal = ', normal)
                            print('forces[i, :] = ', forces[i, :])
                            raise
                        nface = 4
                else:
                    elem_node_ids = elem.node_ids
                    if elem.type == 'CTETRA':
                        #face1 = elem.get_face(load.g1_ref.nid, load.g34_ref.nid)
                        facn = elem.get_face_area_centroid_normal(
                            load.g1_ref.nid, load.g34_ref.nid)
                        face, area, face_centroid, normal = facn
                        #assert face == face1
                        nface = 3
                    elif elem.type == 'CHEXA':
                        #face1 = elem.get_face(load.g34_ref.nid, load.g1_ref.nid)
                        facn = elem.get_face_area_centroid_normal(
                            load.g34_ref.nid, load.g1_ref.nid)
                        face, area, face_centroid, normal = facn
                        #assert face == face1
                        nface = 4
                    elif elem.type == 'CPENTA':
                        g1 = load.g1_ref.nid
                        if load.g34 is None:
                            #face1 = elem.get_face(g1)
                            facn = elem.get_face_area_centroid_normal(g1)
                            face, area, face_centroid, normal = facn
                            nface = 3
                        else:
                            #face1 = elem.get_face(g1, load.g34.nid)
                            facn = elem.get_face_area_centroid_normal(g1, load.g34_ref.nid)
                            face, area, face_centroid, normal = facn
                            nface = 4
                        #assert face == face1
                    #elif elem.type == 'CPYRAM':
                    else:
                        msg = (f'case={load_case_id} eid={eid} etype={elem.type!r} '
                               f'loadtype={load.type!r} not supported')
                        model.log.debug(msg)
                        continue

                    pressures = load.pressures[:nface]
                    assert len(pressures) == nface
                    cid = load.Cid()
                    if load.surf_or_line == 'SURF':
                        pressure = _mean_pressure_on_pload4(pressures, load, elem)
                        load_dir = update_pload4_vector(load, normal, cid)

                        f = pressure * area * load_dir * scale
                    else:
                        msg = (f'surf_or_line={load.surf_or_line!r} on PLOAD4 is not supported\n'
                               '{load}')
                        model.log.debug(msg)
                        continue

                    for inid in face:
                        inidi = nid_map[elem_node_ids[inid]]
                        nodal_pressures[inid] += pressure * scale / nface
                        forces[inidi, :] += f / nface
                    centroidal_pressures[ie] += pressure

                    #r = centroid - p
                    #load.cid.transformToGlobal()
                    #m = cross(r, f)
                    #M += m
            if eids_missing:
                model.log.error('missing PLOAD4 element ids=%s on:\n%s' % (
                    eids_missing, load.rstrip()))

        elif load_type == 'SPCD':
            #self.nodes = [integer(card, 2, 'G1'),]
            #self.constraints = [components_or_blank(card, 3, 'C1', 0)]
            #self.enforced = [double_or_blank(card, 4, 'D1', 0.0)]
            for nid, c1, d1 in zip(load.node_ids, load.components, load.enforced):
                if nid in dependents_nodes:
                    fail_nids.add(nid)
                    fail_count += 1
                    if fail_count < fail_count_max:
                        model.log.warning('    nid=%s is a dependent node and has an'
                                          ' SPCD applied\n%s' % (nid, str(load)))
                c1 = int(c1)
                assert c1 in [1, 2, 3, 4, 5, 6], c1
                if c1 < 4:
                    spcd[nid_map[nid], c1 - 1] = d1
        elif load_type == 'SLOAD':
            for nid, mag in zip(load.nodes, load.mags):
                forces[nid_map[nid]] += np.array([mag, 0., 0.])
        else:
            if load_type not in cards_ignored:
                cards_ignored.add(load_type)
                model.log.warning('  get_forces_moments_array - unsupported '
                                  f'load.type = {load_type}')
    if fail_count:
        fail_nids_list = list(fail_nids)
        fail_nids_list.sort()
        model.log.warning('fail_nids = %s' % np.array(fail_nids_list))
    return centroidal_pressures, forces, spcd

def get_pressure_array(model: BDF, load_case_id, eids, stop_on_failure=True):
    """
    Gets the shell pressures for a load case.

    Parameters
    ----------
    load_case_id : int
        the load case to get the pressure contour for
    eids : (nelements, ) int ndarray
        the element ids in sorted order
    stop_on_failure : bool; default=True
        crashes if the load_case_id doesn't exist

    Returns
    -------
    is_pressure : bool
        is there pressure data
    pressures : (nelements, 1) float ndarray / None
        ndarray : the centroidal pressures
        None : corresponds to is_pressure=False

    """
    if not any(['PLOAD' in model.card_count, 'PLOAD2' in model.card_count,
                'PLOAD4' in model.card_count]):
        return False, None
    cards_ignored = set()
    pressure_loads = ['PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4']

    if not isinstance(load_case_id, integer_types):
        msg = 'load_case_id must be an integer; type=%s, load_case_id=\n%r' % (
            type(load_case_id), load_case_id)
        raise TypeError(msg)

    loads, scale_factors = model.get_reduced_loads(
        load_case_id, stop_on_failure=stop_on_failure)[:2]
    if len(scale_factors) == 0:
        return False, None
    pressures = np.zeros(len(model.elements), dtype='float32')

    shells = {
        'CTRIA3', 'CTRIA6', 'CTRIA', 'CTRIAR',
        'CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR'}
    etypes_skipped = set()
    iload = 0
    nloads = len(loads)
    show_nloads = nloads > 5000
    # loop thru scaled loads and plot the pressure
    for load, scale in zip(loads, scale_factors):
        if show_nloads and iload % 5000 == 0:
            model.log.debug(f'  NastranIOv iload={iload:d}/{nloads:d}')
        if load.type == 'PLOAD4':
            #print(load.object_attributes())
            eids_missing = []
            for elem in load.eids_ref:
                #elem = model.elements[eid]
                if isinstance(elem, integer_types):
                    eids_missing.append(elem)
                    # Nastran is NOT OK with missing element ids
                    continue

                if elem.type in shells:
                    pressure = load.pressures[0] * scale

                    # single element per PLOAD
                    #eid = elem.eid
                    #pressures[eids.index(eid)] = pressure

                    # multiple elements
                    #for elem in load.eids:
                    ie = np.searchsorted(eids, elem.eid)
                    #pressures[ie] += p  # correct; we can't assume model orientation
                    pressures[ie] += pressure

                #elif elem.type in ['CTETRA', 'CHEXA', 'CPENTA']:
                    #A, centroid, normal = elem.get_face_area_centroid_normal(
                        #load.g34_ref.nid, load.g1_ref.nid)
                    #r = centroid - p
                else:
                    etypes_skipped.add(elem.type)
            if eids_missing:
                model.log.error('missing PLOAD4 element ids=%s on:\n%s' % (
                    eids_missing, load.rstrip()))

        elif load.type == 'PLOAD2':
            pressure = load.pressure * scale  # there are 4 pressures, but we assume p0
            for eid in load.eids:
                elem = model.elements[eid]
                ie = np.searchsorted(eids, elem.eid)
                pressures[ie] += pressure

        #elif load.type == 'PLOAD1':
            #pass
        #elif load.type == 'PLOAD':
            # applied to a node, not an element...
            #pressures[ie] = load.pressure * scale
        elif load.type not in pressure_loads:
            continue
        elif load.type in pressure_loads:
            if load.type not in cards_ignored:
                cards_ignored.add(load.type)
                model.log.warning('  get_pressure_array - unsupported '
                                  f'load.type = { load.type}')
        #else:
            #pass
        iload += 1

    if len(etypes_skipped):
        model.log.warning(f'skipping pressure on {list(etypes_skipped)}')
    return True, pressures

def get_temperatures_array(model: BDF, load_case_id, nid_map=None, dtype='float32'):
    """
    Builds the temperature array based on thermal cards.

    Parameters
    ----------
    load_case_id : int
        the load id
    nid_map : ???; default=None -> auto
        ???
    dtype : str; default='float32'
        the type of the temperature array

    Returns
    -------
    is_temperatures : bool
        is there temperature data
    temperatures : (nnodes, ) float ndarray
        the temperatures

    """
    if 'TEMP' not in model.card_count:
        return False, None
    is_temperatures = True

    if nid_map is None:
        nid_map = model.nid_map
    loads, scale_factors = model.get_reduced_loads(load_case_id)[:2]
    tempd = model.tempds[load_case_id].temperature if load_case_id in model.tempds else 0.
    temperatures = np.ones(len(nid_map), dtype=dtype) * tempd

    skip_loads = [
        'FORCE', 'FORCE1', 'FORCE2',
        'MOMENT', 'MOMENT1', 'MOMENT2',
        'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
        'GRAV', 'ACCEL', 'ACCEL1', 'GMLOAD',
        'ACSRCE', 'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',
        'RFORCE', 'RFORCE1', 'SPCD', 'DEFORM',
    ]
    for load, scale in zip(loads, scale_factors):
        if load.type in skip_loads:
            continue

        assert scale == 1.0, str(load)
        if load.type == 'TEMP':
            temps_dict = load.temperatures
            for nid, val in temps_dict.items():
                nidi = nid_map[nid]
                temperatures[nidi] = val
        else:
            model.log.debug(load.rstrip())
    return is_temperatures, temperatures

def get_load_arrays(model: BDF, subcase_id, eid_map, node_ids, normals,
                    nid_map=None, stop_on_failure=True):
    """
    Gets the following load arrays

    Loads include:
     - Temperature
     - Pressure (Centroidal)
     - Forces
     - SPCD

    Parameters
    ----------
    model : BDF()
        the BDF object
    subcase_id : int
        the subcase id
    eid_map : Dict[int eid : int index]
        ???
    node_ids : List[int] / int ndarray
        the node ids in sorted order
    normals : (nelements, 3) ndarray?
        the normal vectors for the shells
        what about solids???

    Returns
    -------
    found_load : bool
        a flag that indicates if load data was found
    found_temperature : bool
        a flag that indicates if temperature data was found
    temperature_data : tuple(temperature_key, temperatures)
        temperature_key : str
            One of the following:
              TEMPERATURE(MATERIAL)
              TEMPERATURE(INITIAL)
              TEMPERATURE(LOAD)
              TEMPERATURE(BOTH)
        temperatures : (nnodes, 1) float ndarray
            the temperatures
    load_data : tuple(centroidal_pressures, forces, spcd)
        centroidal_pressures : (nelements, 1) float ndarray
            the pressure
        forces : (nnodes, 3) float ndarray
            the pressure
        spcd : (nnodes, 3) float ndarray
            the SPCD load application

    """
    subcase = model.subcases[subcase_id]
    if nid_map is None:
        nid_map = model.nid_map
    nnodes = len(node_ids)
    assert len(nid_map) == nnodes, 'len(nid_map)=%s nnodes=%s' % (len(nid_map), nnodes)
    is_loads = False
    is_temperatures = False

    load_keys = (
        'LOAD', 'TEMPERATURE(MATERIAL)', 'TEMPERATURE(INITIAL)',
        'TEMPERATURE(LOAD)', 'TEMPERATURE(BOTH)')
    temperature_keys = (
        'TEMPERATURE(MATERIAL)', 'TEMPERATURE(INITIAL)',
        'TEMPERATURE(LOAD)', 'TEMPERATURE(BOTH)')

    centroidal_pressures = None
    forces = None
    spcd = None
    temperature_key = None
    temperatures = None
    for key in load_keys:
        try:
            load_case_id = subcase.get_parameter(key)[0]
        except KeyError:
            # print('no %s for isubcase=%s' % (key, subcase_id))
            continue
        #model.log.debug('key=%s load_case_id=%s' % (key, load_case_id))

        try:
            unused_load_case = model.get_reduced_loads(
                load_case_id, scale=1.,
                consider_load_combinations=True,
                skip_scale_factor0=False,
                stop_on_failure=True,
                msg='')
        except KeyError:
            msg = f'LOAD={load_case_id} not found'
            if stop_on_failure:
                raise KeyError(msg)
            model.log.warning(msg)
            continue

        if key == 'LOAD':
            p0 = np.array([0., 0., 0.], dtype='float32')
            centroidal_pressures, forces, spcd = get_forces_moments_array(
                model, p0, load_case_id,
                eid_map=eid_map,
                nnodes=nnodes,
                normals=normals,
                dependents_nodes=model.node_ids,
                nid_map=nid_map,
                include_grav=False)
            if centroidal_pressures is not None: # or any of the others
                is_loads = True
        elif key in temperature_keys:
            is_temperatures, temperatures = get_temperatures_array(
                model, load_case_id, nid_map=nid_map)
            temperature_key = key
        else:  # pragma: no cover
            raise NotImplementedError(key)
    temperature_data = (temperature_key, temperatures)
    load_data = (centroidal_pressures, forces, spcd)
    return is_loads, is_temperatures, temperature_data, load_data
