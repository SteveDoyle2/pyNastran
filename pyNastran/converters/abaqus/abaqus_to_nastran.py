"""
https://2021.help.altair.com/2021/hwdesktop/hwx/topics/conversion_between_solvers/abaqus_nastran_conversion_mapping.htm
"""
from __future__ import annotations
import os
import sys
from itertools import count
from collections import defaultdict
from typing import cast, TYPE_CHECKING

import numpy as np

from pyNastran.utils.numpy_utils import integer_types # , float_types
from pyNastran.bdf.bdf import BDF, CaseControlDeck, Subcase, CQUAD4
from pyNastran.bdf.mesh_utils.find_closest_nodes import find_closest_nodes
from pyNastran.converters.abaqus.abaqus import (
    Abaqus, Elements, Step,
    read_abaqus, get_nodes_nnodes_nelements)
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger

chexa_face_map = {
    # take the first and 3rd nodes
    1: (1-1, 3-1),  #face 1: 1-2-3-4
    2: (5-1, 7-1),  #face 2: 5-8-7-6
    3: (1-1, 6-1),  #face 3: 1-5-6-2
    4: (2-1, 7-1),  #face 4: 2-6-7-3
    5: (3-1, 8-1),  #face 5: 3-7-8-4
    6: (4-1, 5-1),  #face 6: 4-8-5-1
    # subtract 1 to make it 0 based
}

ctetra_face_map = {
    1: (1-1, 2-1, 3-1), # Face 1: 1-2-3
    2: (1-1, 4-1, 2-1), # Face 2: 1-4-2
    3: (2-1, 4-1, 3-1), # Face 3: 2-4-3
    4: (3-1, 4-1, 1-1), # Face 4: 3-4-1
}

tri_face_map = {
    #for triangular shell elements:
    1: -1.0,        #Face NEG or 1: in negative normal direction
    2: 1.0,         #Face POS or 2: in positive normal direction
    3: (1-1, 2-1),  #Face 3: 1-2
    4: (2-1, 3-1),  #Face 4: 2-3
    5: (3-1, 1-1),  #Face 5: 3-1
}

quad_face_map = {
    #for quadrilateral shell elements:
    1: -1.0,        #Face NEG or 1: in negative normal direction
    2: 1.0,         #Face POS or 2: in positive normal direction
    3: (1-1, 2-1),  #Face 3: 1-2
    4: (2-1, 3-1),  #Face 4: 2-3
    5: (3-1, 4-1),  #Face 5: 3-4
    6: (4-1, 1-1),  #Face 6: 4-1
}

def _add_part_to_nastran(nastran_model: BDF,
                         elements: Elements, pid: int,
                         nid_offset: int, eid_offset: int) -> int:
    log = nastran_model.log

    log.debug('starting part...')
    element_types = {etype: eids_nids
                     for etype, eids_nids in elements.element_types.items()
                     if eids_nids[0] is not None}

    unstacked_eids = [eids_nodes[0] for eids_nodes in element_types.values()]
    all_eids = np.hstack(unstacked_eids)
    ueids = np.unique(all_eids)
    map_eids = (len(all_eids) != len(ueids))
    abaqus_type_to_etype = {
        'b31h': 'CBEAM',
        's3': 'CTRIA3',
        's3r': 'CTRIA3',
        'cpe3': 'CTRIA3',
        'cpe3r': 'CTRIA3',

        's4': 'CQUAD4',
        's4r': 'CQUAD4',
        'cpe4': 'CQUAD4',
        'cpe4r': 'CQUAD4',

        'c3d4': 'CTETRA',
        'c3d4r': 'CTETRA',
        'c3d10': 'CTETRA',
        'c3d10r': 'CTETRA',
        'c3d10h': 'CTETRA',

        'c3d6': 'CPENTA',
        'c3d15': 'CPENTA',
        'c3d6r': 'CPENTA',
        'c3d15r': 'CPENTA',

        'c3d8': 'CHEXA',
        'c3d20': 'CHEXA',
        'c3d8r': 'CHEXA',
        'c3d20r': 'CHEXA',

        'cohax4': None,
        'coh2d4': None,
        'cax3': None,
        'cax4r': None,
        'r2d2': None,
    }

    #pid = -1
    type_eid_map_to_eid = {}
    for etype, eids_nids in element_types.items():
        eids_, part_nids = eids_nids
        if eids_ is None and part_nids is None:
            continue

        # TODO: well that's annoying...abaqus can have duplicate ids
        #if eid_offset:
        #eids = (eid_offset + 1 - eids_.min()) + eids_
        if map_eids:
            nastran_etype = abaqus_type_to_etype[etype]
            eids = eid_offset + np.arange(len(eids_)) + 1
            for eid1, eid2 in zip(eids_, eids):
                type_eid_map_to_eid[(nastran_etype, eid1)] = eid2
        else:
            eids = eids_


        #print(f'eids[{etype} = {eids}; eid_offset={eid_offset}')
        #log.warning(f'writing etype={etype} eids={eids_}->{eids}')
        elset_name = elements.element_type_to_elset_name[etype]
        comment = elset_name
        if nid_offset > 0:
            # don't use += or it's an inplace operation
            part_nids = part_nids + nid_offset

        if etype == 'r2d2':
            log.warning('skipping r2d2; should this be a RBE1/RBAR?')
            continue

        if etype == 'b31h':
            for eid, nids in zip(eids, part_nids):
                x = None
                g0 = nids[2]
                nids = [nids[0], nids[1]]
                nastran_model.add_cbeam(
                    eid, pid, nids, x, g0, offt='GGG', bit=None,
                    pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment='')
        elif etype in {'s3', 's3r', 'cpe3', 'cpe3r'}:
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_ctria3(eid, pid, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                                         T1=None, T2=None, T3=None, comment=comment)
                comment = ''
        elif etype in {'s4', 's4r', 'cpe4', 'cpe4r'}:
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_cquad4(eid, pid, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                                         T1=None, T2=None, T3=None, T4=None, comment=comment)
                comment = ''
        elif etype in {'s6', 's6r'}:
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_ctria6(eid, pid, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                                         T1=None, T2=None, T3=None, comment=comment)
                comment = ''
        elif etype in {'s8', 's8r'}:
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_cquad8(eid, pid, nids, theta_mcid=0.0, zoffset=0., tflag=0,
                                         T1=None, T2=None, T3=None, T4=None, comment=comment)
                comment = ''
        elif etype in {'c3d4', 'c3d4r'}: # , 'c3d10', 'c3d10r', 'c3d10h'}:
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_ctetra(eid, pid, nids, comment=comment)
                comment = ''
        elif etype in {'c3d10', 'c3d10r', 'c3d10h'}:
            #log.warning('found quadratic_tetra')
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_ctetra(eid, pid, nids, comment=comment)
                comment = ''
        elif etype in {'c3d6', 'c3d6r'}:
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_cpenta(eid, pid, nids, comment=comment)
                comment = ''
        elif etype in {'c3d15', 'c3d15r'}:
            log.warning('found quadratic_penta')
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_cpenta(eid, pid, nids, comment=comment)
                comment = ''
        elif etype in {'c3d8', 'c3d8r'}:
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_chexa(eid, pid, nids, comment=comment)
                comment = ''
        elif etype in {'c3d20', 'c3d20r'}:
            for eid, nids in zip(eids, part_nids):
                # the last 8 nodes are flipped in sets of 2
                (n1, n2, n3, n4,
                 n5, n6, n7, n8,
                 n9, n10, n11, n12,
                 n13, n14, n15, n16,
                 n17, n18, n19, n20) = nids
                nids2 = [
                    n1, n2, n3, n4,
                    n5, n6, n7, n8,
                    n9, n10, n11, n12,

                    n17, n18, n19, n20,
                    n13, n14, n15, n16,]
                nastran_model.add_chexa(eid, pid, nids2, comment=comment)
                comment = ''
        elif etype in {'cohax4', 'coh2d4', 'cax3', 'cax4r'}:
            del eids, part_nids
            log.warning(f'skipping etype={etype!r}')
            continue
        else:
            raise NotImplementedError(etype)
        eid_offset += len(eids)
        del eids, part_nids
        #print(etype, eid_offset)

    #add_lines(grid, nidsi, part.r2d2, nid_offset)

    #add_tris(grid, nidsi, part.cps3, nid_offset)
    #add_tris(grid, nidsi, part.cpe3, nid_offset)

    #add_quads(grid, nidsi, part.cpe4, nid_offset)
    #add_quads(grid, nidsi, part.cpe4r, nid_offset)

    #add_quads(grid, nidsi, part.cps4, nid_offset)
    #add_quads(grid, nidsi, part.cps4r, nid_offset)

    #add_quads(grid, nidsi, part.coh2d4, nid_offset)
    #add_quads(grid, nidsi, part.cohax4, nid_offset)

    #add_tris(grid, nidsi, part.cax3, nid_offset)
    #add_quads(grid, nidsi, part.cax4, nid_offset)
    #add_quads(grid, nidsi, part.cax4r, nid_offset)

    ## solids
    #add_tetras(grid, nidsi, part.c3d10h, nid_offset)
    #add_hexas(grid, nidsi, part.c3d8r, nid_offset)
    assert eid_offset > 0, eid_offset
    return eid_offset


def _create_nastran_nodes_elements(model: Abaqus,
                                   nastran_model: BDF) -> None:
    log = model.log
    nid_offset = 0
    eid_offset = 0

    pid = 1
    if model.nids is not None and len(model.nids):
        nnodesi = model.nodes.shape[0]
        elements = model.elements
        eid_offset = _add_part_to_nastran(
            nastran_model, elements, pid, nid_offset, eid_offset)
        _build_rigid_ties(model, nastran_model, elements, eid_offset)
        nid_offset += nnodesi

    eid_offset = 0
    for unused_part_name, part in model.parts.items():
        log.warning(f'part_name = {unused_part_name!r} eid_offset={eid_offset:d}')
        nnodesi = part.nodes.shape[0]
        #nidsi = part.nids

        elements = part.elements
        eid_offset = _add_part_to_nastran(
            nastran_model, elements, pid, nid_offset, eid_offset)
        #nids.append(nidsi)
        assert eid_offset > 0, eid_offset
        nid_offset += nnodesi
        for shell_section in part.shell_sections:
            log.info('shell')
        for shell_section in part.solid_sections:
            log.info('solid')
    #nids = np.hstack(nids)

    pid = 1
    mid = 1
    cid = 1
    mat_name_to_mid_dict = {}
    pid, mid = _create_solid_properties(
        model, nastran_model, log, pid, mid,
        mat_name_to_mid_dict)
    pid, mid, cid = _create_shell_properties(
        model, nastran_model, log, pid, mid, cid,
        mat_name_to_mid_dict)


def _build_rigid_ties(model: Abaqus,
                      nastran_model: BDF,
                      elements: Elements,
                      eid_offset: int) -> int:
    ties = model.ties
    if len(ties) == 0:
        return eid_offset
    log = model.log
    eid = eid_offset + 1

    log.error(f'Ties found...not supported\n{ties}')
    #elset_name_to_element_type = {
        #value: key for key, value in
        #elements.element_type_to_elset_name.items()}

    all_nodes_xyz = model.nodes
    all_nids = model.nids

    # preallocate to avoid doing some faces multiple times
    face_to_nids_map = {}
    for tie in ties:
        master_surface = model.surfaces[tie.master]
        slave_surface = model.surfaces[tie.slave]
        assert master_surface.surface_type == 'element', master_surface

        for set_name, face in zip(master_surface.set_names, master_surface.faces):
            face_to_nids_map[(master_surface.surface_type, set_name, face)] = []
        for set_name, face in zip(slave_surface.set_names, slave_surface.faces):
            face_to_nids_map[(slave_surface.surface_type, set_name, face)] = []

    _build_face_to_nids(model, elements, face_to_nids_map)

    for tie in ties:
        master_surface = model.surfaces[tie.master]
        slave_surface = model.surfaces[tie.slave]
        assert master_surface.surface_type == 'element', master_surface

        master_nids_list = []
        slave_nids_list = []
        for set_name, face in zip(master_surface.set_names, master_surface.faces):
            mnids = face_to_nids_map[(master_surface.surface_type, set_name, face)]
            master_nids_list.append(mnids)
        for set_name, face in zip(slave_surface.set_names, slave_surface.faces):
            snids = face_to_nids_map[(slave_surface.surface_type, set_name, face)]
            slave_nids_list.append(snids)
        master_nids = np.unique(np.hstack(master_nids_list))
        slave_nids = np.unique(np.hstack(slave_nids_list))

        #face_to_nids_map[(master_surface.surface_type, )]
        imaster = np.searchsorted(all_nids, master_nids)
        islave = np.searchsorted(all_nids, slave_nids)
        master_xyz = all_nodes_xyz[imaster, :]
        slave_xyz = all_nodes_xyz[islave, :]
        # 1 independent, 1 dependent
        try:
            master_nids2 = find_closest_nodes(
                master_xyz, master_nids, slave_xyz,
                neq_max=1, tol=1e-8, msg='')
        except IndexError:
            log.error(f'skipping contact face for tie={tie.name}')
            continue

        cm = '123456'
        #eid = eid_offset + 1
        for i, master_nid, slave_nid in zip(count(), master_nids, slave_nids):
            gn = master_nid
            Gmi = slave_nid
            nastran_model.add_rbe2(
                eid, gn,  # independent
                cm, Gmi,  # dependent
                alpha=0.0, tref=0.0, comment='', validate=False)
            eid += 1
        #islave2 = np.searchsorted(all_nids, slave_nids2)
        #slave_xyz2 = all_nodes_xyz[islave, :]
        #master_nids2 = find_closest_nodes(
            #slave_xyz2, slave_nids2, master_xyz,
            #neq_max=1, tol=None, msg='')
    eid_offset -= 1
    return eid_offset

def _build_face_to_nids(
        model: Abaqus,
        elements: Elements,
        face_to_nids_map: dict[tuple[str, str, str], np.ndarray]) -> None:
    log = model.log

    for set_type, set_name, face in list(face_to_nids_map.keys()):
        face_int = int(face[-1])
        if set_type == 'element':
            eids = model.element_sets[set_name]
            #eids = elset_name_to_element_type[set_name]
        else:
            raise NotImplementedError(set_type)
        nids_list = []
        for etypei, seti in elements.element_type_to_elset_name.items():
            element_eids = getattr(elements, f'{etypei}_eids')
            common_eids = np.intersect1d(element_eids, eids)
            if len(common_eids):
                element_nodes = getattr(elements, etypei)
                ielement = np.searchsorted(element_eids, common_eids)
                if etypei == 'c3d10':
                    #face_int = int(face[-1])
                    in1, in2, in3 = ctetra_face_map[face_int]
                    element_nodes2 = element_nodes[ielement, :][:, [in1, in2, in3]]
                elif etypei == 's6':
                    #face_int = int(face[-1])
                    #if face_int == 4:
                        #log.error(f'etype=s6 and has a face_id={face} for eids={common_eids}?')
                        #continue
                    in1, in2 = tri_face_map[face_int]
                    element_nodes2 = element_nodes[ielement, :][:, [in1, in2]]

                elif etypei in {'cpe3', 'cps3'}:
                    #for triangular plane stress, plane strain and axisymmetric elements:
                    #Face 1: 1-2
                    #Face 2: 2-3
                    #Face 3: 3-1
                    #Face N: in negative normal direction (only for plane stress)
                    #Face P: in positive normal direction (only for plane stress)
                    asdf
                #elif etypei in {'cpe4', 'cps4'}:
                    #asdf
                elif etypei == 's8':
                    in1, in2 = quad_face_map[face_int]
                    element_nodes2 = element_nodes[ielement, :][:, [in1, in2]]
                else:
                    raise RuntimeError(etypei)
                ravel_element_nodes = element_nodes2.ravel()
                nids_list.append(ravel_element_nodes)

        nids = np.hstack(nids_list)
        face_to_nids_map[(set_type, set_name, face)] = nids
        x = 1
    return

def _create_solid_properties(model: Abaqus, nastran_model: BDF,
                             log: SimpleLogger,
                             pid: int, mid: int,
                             mat_name_to_mid_dict: dict[str, int]) -> tuple[int, int]:
    for solid_section in model.solid_sections:
        #print(solid_section)

        mat_name = solid_section.material_name
        mat = model.materials[mat_name]
        element_set_name = solid_section.elset
        log.info(f'element_set_name = {element_set_name}')

        etypes_eids = map_solid_property_ids(model, nastran_model, element_set_name, pid)
        for etype, eids in etypes_eids.items():
            eid = eids[0]
            elem = nastran_model.elements[eid]
            etype_n = (etype, len(elem.nodes))
            if etype_n in {('CTETRA', 10), ('CHEXA', 20), ('CPENTA', 15), ('CPYRAM', 13)}:
                isop = 'FULL'
            else:
                assert etype_n in {('CTETRA', 4), ('CHEXA', 8), ('CPENTA', 6), ('CPYRAM', 5)}, etype_n
                isop = 'REDUCED'
            nastran_model.add_psolid(
                pid, mid, cordm=0,
                integ=None, stress=None, isop=isop,
                fctn='SMECH', comment=element_set_name + f'; etype={etype}')
            _create_material(nastran_model, mat, mid, comment=mat_name+f' for {etype}')
            pid += 1
            mid += 1

        #if mat_name in mat_name_to_mid_dict:
            #midi = mat_name_to_mid_dict[mat_name]
        #else:
            #mat = model.materials[mat_name]
            #mat = _create_material(nastran_model, mat, mid, comment=mat_name)
            #midi = mid
            #delta_mid = 1

    return pid, mid

def build_coord(model: Abaqus,
                nastran_model: BDF,
                cid: int,
                orientation_name: str) -> int:

    orient = model.orientations[orientation_name]
    if orient.axis is None:
        R = np.eye(3)
        comment = orient.name
    else:
        c = np.cos(orient.alpha)
        s = np.sin(orient.alpha)
        if orient.axis == 1:
            R = np.array([
                [1., 0., 0.],
                [0.,  c, -s],
                [0.,  s,  c],
            ])
        elif orient.axis == 2:
            R = np.array([
                [c,  0., -s],
                [0., 1., 0.],
                [-s, 0.,  c],
            ])
        elif orient.axis == 3:
            R = np.array([
                [c,  -s, 0.],
                [s,   c, 0.],
                [0., 0., 1.],
            ])
        else:  # pragma: no cover
            raise NotImplementedError(orient.axis)
        comment = orient.name + f'\n  axis={orient.axis}; alpha={orient.alpha}'

    if orient.system == 'rectangular':
        origin = orient.origin
        i = orient.x_axis - origin
        ijplane = orient.xy_plane - origin
        i /= np.linalg.norm(i)
        ijplane /= np.linalg.norm(ijplane)

        k = np.cross(i, ijplane)
        k /= np.linalg.norm(k)

        j = np.cross(k, i)
        beta = np.vstack([i, j, k])
        beta2 = R @ beta

        xzplane = origin + beta2[0, :]  # xaxis, xzplane
        #yaxis = beta2[:, 1]
        zaxis = origin + beta2[2, :]

        coord = nastran_model.add_cord2r(
            cid, origin, zaxis, xzplane,
            rid=0, setup=True, comment=comment)
        x = 1
    else:
        raise RuntimeError(orient.system)

    return cid

def _create_shell_properties(model: Abaqus, nastran_model: BDF,
                             log: SimpleLogger,
                             pid: int, mid: int, cid: int,
                             mat_name_to_mid_dict: dict[str, int]) -> tuple[int, int]:
    for shell_section in model.shell_sections:
        coord = -1
        #print(shell_section)
        orientation_name = shell_section.orientation
        element_set_name = shell_section.elset
        log.info(f'element_set_name = {element_set_name}')

        shell_comment = element_set_name
        if orientation_name:
            build_coord(model, nastran_model,
                        cid, orientation_name)
            coord = cid
            shell_comment += f'; mcid={cid}'
            cid += 1

        delta_mid = 0
        mat_name = shell_section.material_name
        if isinstance(mat_name, str):
            # PSHELL
            if mat_name in mat_name_to_mid_dict:
                midi = mat_name_to_mid_dict[mat_name]
            else:
                midi = mid
                mat = model.materials[mat_name]
                mat = _create_material(nastran_model, mat, mid, comment=mat_name)
                mat_name_to_mid_dict[mat_name] = midi
                delta_mid = 1

            #*Shell section, Elset=Base, Orientation=OR1, Material=CF, Offset=0
            #1.
            t = shell_section.thickness
            nastran_model.add_pshell(
                pid, mid1=midi, t=t, mid2=midi, twelveIt3=1.0, mid3=None, tst=0.833333,
                nsm=0.0, z1=None, z2=None, mid4=None, comment=shell_comment)
        else:
            mat_names = shell_section.material_name
            #*Shell section, Elset=Internal_Selection-1_Shell_section-1, COMPOSITE
            #0.25,,Steel
            #0.25,,Steel
            umat_names = np.unique(mat_names)
            mids = np.zeros(len(mat_names), dtype='int32')
            imid = 0
            for mat_name in umat_names:
                imidi = (mat_name == mat_name)

                if mat_name in mat_name_to_mid_dict:
                    midi = mat_name_to_mid_dict[mat_name]
                else:
                    midi = mid + imid
                    mat = model.materials[mat_name]
                    unused_material = _create_material(nastran_model, mat, midi, comment=mat_name)
                    mat_name_to_mid_dict[mat_name] = midi
                    imid += 1
                    delta_mid += 1
                mids[imidi] = midi

            thicknesses = shell_section.thickness
            nastran_model.add_pcomp(
                pid, mids, thicknesses, thetas=None, souts=None,
                nsm=0., sb=0., ft=None, tref=0., ge=0.,
                lam=None, z0=None, comment=shell_comment)

        zoffset = shell_section.offset
        map_shell_property_ids(model, nastran_model, element_set_name, pid,
                               coord=coord, zoffset=zoffset)
        pid += 1
        mid += delta_mid
    return pid, mid, cid

def map_solid_property_ids(model: Abaqus, nastran_model: BDF,
                           element_set_name: str, pid: int) -> dict[str, list[int]]:
    log = model.log
    try:
        eids = get_eids_from_recursive_element_set(model, element_set_name)
    except KeyError:
        log.error(f'cant map section with elset={element_set_name}')
        return {}

    if isinstance(eids, str):
        log.debug(f'mapping section with elset={element_set_name} to {eids}')
        element_name_to_type = {value: key for key, value in
                                model.elements.element_type_to_elset_name.items()}
        etype = element_name_to_type[eids]
        eids = getattr(model.elements, f'{etype}_eids')
        #return
    etypes_eids = defaultdict(list)
    for eid in eids:
        elem = nastran_model.elements[eid]
        etypes_eids[elem.type].append(eid)

    for etype, eids in etypes_eids.items():
        for eid in eids:
            elem = nastran_model.elements[eid].pid = pid
        pid += 1
    #return pid
    return etypes_eids

def map_shell_property_ids(model: Abaqus, nastran_model: BDF,
                           element_set_name: str, pid: int,
                           coord: int=-1, zoffset: float=0.0) -> None:
    log = model.log
    try:
        eids = get_eids_from_recursive_element_set(model, element_set_name)
    except KeyError:
        log.error(f'cant map section with elset={element_set_name}')
        return

    if isinstance(eids, str):
        log.debug(f'mapping section with elset={element_set_name} to {eids}')
        element_name_to_type = {value: key for key, value in
                                model.elements.element_type_to_elset_name.items()}
        etype = element_name_to_type[eids]
        eids = getattr(model.elements, f'{etype}_eids')
        #return

    theta_mcid = 0. if coord == -1 else coord
    for eid in eids:
        elem = nastran_model.elements[eid]
        elem.pid = pid
        elem.zoffset = zoffset
        elem.theta_mcid = theta_mcid

def get_eids_from_recursive_element_set(model: Abaqus, element_set_name: str) -> Union[np.ndarray, str]:
    """returns None if we cant find a set of eids"""
    eids = model.element_sets[element_set_name]
    if isinstance(eids, str):
        return eids.lower()
    return eids

def _create_material(nastran_model: BDF, mat, mid: int, comment: str) -> Union[MAT1, MAT8]:
    G = None
    rho = 0.0
    tref = None
    alpha = None

    sections = mat.sections
    assert 'elastic' in sections or 'engineering constants' in sections, sections

    if 'density' in sections:
        densities = sections['density']
        assert len(densities) == 1, sections
        rho = densities[0]
    if 'expansion' in sections:
        tref, alpha = sections['expansion']

    if 'elastic' in sections:
        elastic = sections['elastic']
        E = elastic[0]
        nu = elastic[1]
        material = nastran_model.add_mat1(
            mid, E, G, nu, rho=rho, a=alpha, tref=tref,
            ge=0.0, St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment=comment)
    elif 'engineering constants' in sections:
        e11, e22, e3, nu12, nu13, n23, g12, g13, \
            g23, tref = sections['engineering constants']
        material = nastran_model.add_mat8(
            mid, e11, e22, nu12, g12=g12, g1z=g13, g2z=g23,
            rho=rho, a1=0., a2=0., tref=0.,
            Xt=0., Xc=None, Yt=0., Yc=None, S=0.,
            ge=0., F12=0., strn=0., comment=comment)
    else:
        raise RuntimeError(sections)
    return material


def abaqus_to_nastran_filename(abaqus_inp_filename: str,
                               nastran_filename_out: str,
                               size: int=8,
                               xform: bool=False,
                               encoding: Optional[str]=None,
                               log: Optional[SimpleLogger]=None) -> BDF:
    if isinstance(abaqus_inp_filename, Abaqus):
        model = abaqus_inp_filename
    else:
        model = read_abaqus(abaqus_inp_filename, encoding=encoding, log=log, debug=True)
    log = model.log
    if xform:
        log.error('xform is super buggy and largely untested')

    nnodes, nids, nodes, nelements = get_nodes_nnodes_nelements(
        model, stop_for_no_elements=True)
    assert nnodes > 0, nnodes
    assert nelements > 0, nelements

    nastran_model = BDF(debug=True, log=log, mode='msc')
    comment = '\n'.join(model.heading)
    nastran_model.add_param('POST', -1, comment=comment)
    for nid, xyz in zip(nids, nodes):
        nastran_model.add_grid(nid, xyz)
    _create_nastran_nodes_elements(model, nastran_model)
    #for step in model.steps:
        #print(step)

    _create_nastran_loads(model, nastran_model)
    try:
        str(nastran_model.case_control_deck)
    except AssertionError:
        log.error('No case control deck found...skipping')
        #print(f'{nastran_model.case_control_deck}')
        nastran_model.case_control_deck = None

    _xform_model(nastran_model, xform)
    nastran_model.write_bdf(
        nastran_filename_out, size=size,
        encoding=None,
        #nodes_size=None, elements_size=None,
        #loads_size=None, is_double=False, interspersed=False,
        enddata=True, write_header=True, close=True)

    x = 1
    return nastran_model

def _xform_model(nastran_model: BDF, xform: bool):
    if not xform:
        return
    log = nastran_model.log
    log.error('xform is super buggy and largely untested')
    nastran_model.cross_reference()
    for eid, elem in nastran_model.elements.items():
        if hasattr(elem, 'theta_mcid') and isinstance(elem.theta_mcid, integer_types):
            # mcid is active
            #elem = cast(CQUAD4, elem)
            dxyz2, centroid, imat, jmat, normal1 = elem.material_coordinate_system()
            dxyz, centroid, ielement, jelement, normal2 = elem.element_coordinate_system()
            # angle between two normalized vectors
            # a o b = |a| |b| cos(theta)
            # theta = acos(a o b)
            theta = np.arccos(np.clip(np.dot(imat, ielement), -1, 1))
            theta_deg = np.degrees(theta)
            elem.theta_mcid = theta_deg
            elem.uncross_reference()
    log.error('xform is super buggy and largely untested')

def _write_boundary_as_nastran(model: Abaqus,
                               step: Abaqus,
                               nastran_model: BDF,
                               case_control_deck: CaseControlDeck) -> None:
    if not step.boundaries:
        return
    comment_list = []
    for iboundary, boundary in enumerate(step.boundaries):
        fixed_spcs = defaultdict(list)
        for (nid, dof), value in boundary.nid_dof_to_value.items():
            if isinstance(nid, str) and nid not in comment_list:
                comment_list.append(nid)

            nids = _get_nodes(model, nid)
            #assert isinstance(nid, int), nid
            for nid in nids:
                if value == 0.0:
                    fixed_spcs[dof].append(nid)
                else:
                    raise RuntimeError((nid, dof, value))

        subcase_id = iboundary + 1
        spc_id = iboundary + 1
        subcase = case_control_deck.create_new_subcase(subcase_id)
        subcase.add('SPC', spc_id, [], 'STRESS-type')
        for dof, nids in fixed_spcs.items():
            #  spc_id: int, components: str, nodes
            nastran_model.add_spc1(spc_id, str(dof), nids)
    if comment_list:
        comment = '\n'.join(comment_list)
        nastran_model.spcs[spc_id][0].comment = comment

def _create_nastran_loads(model: Abaqus, nastran_model: BDF):
    log = nastran_model.log
    if nastran_model.case_control_deck is None:
        nastran_model.case_control_deck = CaseControlDeck([], log=nastran_model.log)

    case_control_deck = nastran_model.case_control_deck
    _write_boundary_as_nastran(model, model, nastran_model, case_control_deck)

    output_map = {
        'U': 'DISP',
        'RF': 'SPCFORCE',  # Reaction Force
        'S': 'STRESS',
        'E': 'STRAIN',
        'ENER': 'ESE',
        #'NOD': disables the error estimator; nastran doesn't have one :)
        #'NOE': disables the error estimator; nastran doesn't have one :)
    }
    for istep, step in enumerate(model.steps):
        subcase_id = istep + 1
        load_id = subcase_id
        if subcase_id in nastran_model.subcases:
            subcase = nastran_model.subcases[subcase_id]
        else:
            subcase = case_control_deck.create_new_subcase(subcase_id)

        _write_boundary_as_nastran(model, step, nastran_model, case_control_deck)
        for output in step.node_output + step.element_output:
            output_upper = output.upper()
            if output_upper in {'NOD', 'NOE'}:
                log.warning(f'skipping output request={output}')
                continue

            try:
                base_output = output_map[output_upper]
            except KeyError:
                raise KeyError(f'output={output!r} is not in [u, rf, s, e, ener]')

            request_flags = ['PLOT', 'PRINT']
            if base_output in {'STRESS', 'STRAIN'}:
                request_flags.append('CENTER')
            subcase.add(base_output, 'ALL', request_flags, 'STRESS-type')

        comment_list = []
        _write_distributed_loads(model, step, nastran_model, subcase, load_id, comment_list)
        _write_concentrated_loads(model, step, nastran_model, subcase, load_id, comment_list)
        if comment_list:
            comment = '\n'.join(comment_list)
            nastran_model.loads[load_id][0].comment = comment

def _write_concentrated_loads(model: Abaqus,
                              step: Step,
                              nastran_model: BDF, subcase: Subcase,
                              load_id: int,
                              comment_list: list[str]) -> None:
    def fxyz():
        return np.zeros(3, dtype='float32')

    for cload in step.cloads:
        assert nastran_model.sol is None or nastran_model.sol == 101, nastran_model.sol
        nastran_model.sol = 101
        subcase.add('LOAD', load_id, [], 'STRESS-type')
        #subcase['LOAD'] = load_id
        forces = defaultdict(fxyz)
        moments = defaultdict(fxyz)
        for cloadi in cload:
            nid, dof, mag = cloadi
            if isinstance(nid, str) and nid not in comment_list:
                comment_list.append(nid)
            nids = _get_nodes(model, nid)
            assert dof in [1, 2, 3], cload

            for nid in nids:
                if dof in {1, 2, 3}:
                    forces[nid][dof - 1] = mag
                elif dof in {4, 5, 6}:
                    moments[nid][dof - 4] = mag
                else:
                    raise NotImplementedError(cloadi)

        mag = 1.0
        if len(forces):
            for nid, xyz in forces.items():
                assert isinstance(nid, integer_types), nid
                nastran_model.add_force(load_id, nid, mag, xyz, cid=0, comment='')
        if len(moments):
            for nid, xyz in moments.items():
                assert isinstance(nid, integer_types), nid
                nastran_model.add_moment(load_id, nid, mag, xyz, cid=0, comment='')

        #print(step.cloads)
    #step.cloads

def _write_distributed_loads(model: Abaqus,
                             step: Step,
                             nastran_model: BDF, subcase: Subcase,
                             load_id: int,
                             comment_list: list[str]) -> None:
    log = model.log
    for dload in step.dloads:
        assert nastran_model.sol is None or nastran_model.sol == 101, nastran_model.sol
        nastran_model.sol = 101
        subcase.add('LOAD', load_id, [], 'STRESS-type')
        #pressures = []
        for dloadi in dload:
            eid, tag = dloadi[:2]
            if tag in {'P1', 'P2', 'P3', 'P4', 'P5', 'P6'}:
                face = tag
                assert len(dloadi) == 3, dloadi
                mag = dloadi[-1]
                eids = _get_elements(model, eid)
                pressure = mag

                if isinstance(face, str):
                    assert face[0] == 'P' and len(face) == 2, face

                    #cquad4_cquad8_face_map = {
                        #1: (1-1, 2-1),  # Face 1: 1-2
                        #2: (2-1, 3-1),  # Face 2: 2-3
                        #3: (2-1, 3-1),  # Face 3: 3-4
                        #4: (2-1, 3-1),  # Face 4: 4-1
                    #}
                    pressures = [pressure, None, None, None]
                    face_int = int(face[-1])
                    for eid in eids:
                        element = nastran_model.elements[eid]
                        if element.type == 'CHEXA':
                            i1, i3 = chexa_face_map[face_int]
                            n1 = element.nodes[i1]
                            n3 = element.nodes[i3]
                            nastran_model.add_pload4(
                                load_id, [eid], pressures, g1=n1, g34=n3, cid=None, nvector=None,
                                surf_or_line='SURF', line_load_dir='NORM', comment='')
                        elif element.type == 'CTRIA6':
                            mapped_face = tri_face_map[face_int]
                            if isinstance(mapped_face, float):
                                pressures = [mapped_face*pressure, None, None, None]
                            elif isinstance(mapped_face, tuple):
                                log.error(f'DLOAD uses P{face_int} for eid={eid} {element.type}?  Assuming normal')
                            else:
                                raise RuntimeError(mapped_face)
                            #assert face_int == 1, face_int
                            nastran_model.add_pload4(
                                load_id, [eid], pressures, g1=None, g34=None, cid=None, nvector=None,
                                surf_or_line='SURF', line_load_dir='NORM', comment='')
                        elif element.type == 'CQUAD8':
                            #if face_int != 1:
                            mapped_face = quad_face_map[face_int]
                            if isinstance(mapped_face, float):
                                pressures = [mapped_face*pressure, None, None, None]
                            elif isinstance(mapped_face, tuple):
                                log.error(f'DLOAD uses P{face_int} for eid={eid} {element.type}?  Assuming normal')
                            else:
                                raise RuntimeError(mapped_face)

                            #assert face_int == 1, face_int
                            nastran_model.add_pload4(
                                load_id, [eid], pressures, g1=None, g34=None, cid=None, nvector=None,
                                surf_or_line='SURF', line_load_dir='NORM', comment='')
                        else:
                            raise NotImplementedError(element)
                else:
                    for eid in eids:
                        asfd
                        nastran_model.add_pload2(load_id, pressure, [eid], comment='')
                        #nastran_model.add_pload(load_id, pressure, nodes, comment='')

            elif tag == 'GRAV':
                assert len(dloadi) == 6, dloadi
                mag, gx, gy, gz = dloadi[2:]
                N = [gx, gy, gz]
                nastran_model.add_grav(load_id, mag, N, cid=0, mb=0, comment='')
            else:
                raise RuntimeError(dloadi)

def _get_nodes(model: Abaqus, nid: Union[int, str]) -> list[int]:
    if isinstance(nid, integer_types):
        nids = [nid]
    else:
        nids = model.node_sets[nid.lower()]
    return nids

def _get_elements(model: Abaqus, eid: Union[int, str]) -> list[int]:
    if isinstance(eid, integer_types):
        eids = [eid]
    else:
        eids = model.element_sets[eid.lower()]
    return eids

def cmd_abaqus_to_nastran(argv=None, log: Optional[SimpleLogger]=None,
                          quiet: str=False) -> None:
    """Interface for abaqus_to_nastran"""
    if argv is None:
        argv = sys.argv

    default_encoding = sys.getdefaultencoding()
    other = '[--large] [--encoding ENCODING] [--debug]'
    other += ' [--xform]'
    msg = (
        'Usage:\n'
        f'  abaqus_to_nastran ABAQUS_INP_IN                 {other}\n'
        f'  abaqus_to_nastran ABAQUS_INP_IN NASTRAN_BDF_OUT {other}\n'
        '  abaqus_to_nastran -h | --help\n'
        '  abaqus_to_nastran -v | --version\n'
        '\n'
        'Required Arguments:\n'
        '  ABAQUS_INP_IN       path to abaqus.inp file\n'
        '  NASTRAN_BDF_OUT     path to nastran.bdf file (default=abaqus.bdf)\n'
        '\n'

        'Nastran Options:\n'
        '  --large             writes the data in large field format\n'
        '\n'

        'Abaqus Options:\n' # 'utf8bom' = 'utf-8-sig'
        f'  --encoding ENCODING  Specify the encoding (e.g., latin1, cp1252, utf8, utf-8-sig); default={default_encoding!s}\n'
        f'  --debug              Turns on debugging\n'
        '\n'

        'Info:\n'
        '  -h, --help     show this help message and exit\n'
        "  -v, --version  show program's version number and exit\n"
        '\n'
        'Examples:\n'
        '  # creates model.bdf\n'
        '  abaqus_to_nastran model.inp\n\n'
        '  # creates fem.bdf\n'
        '  abaqus_to_nastran model.inp fem.bdf\n\n'
        '  # creates model.bdf with large field format\n'
        '  abaqus_to_nastran model.inp --large\n\n'
        '  # creates model.bdf an alternate encoding\n'
        '  abaqus_to_nastran model.inp --encoding utf-8-sig\n\n'
        '\n'
    )
    from docopt import docopt
    import pyNastran
    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver, argv=argv[1:])

    encoding = default_encoding
    if data['--encoding']:
        encoding = data['--encoding']
    if not quiet:  # pragma: no cover
        print(data)
    abaqus_inp_filename = data['ABAQUS_INP_IN']
    nastran_filename_out = os.path.splitext(abaqus_inp_filename)[0] + '.bdf'

    if log is None:
        is_debug = data['--debug']
        default_level = 'debug' if is_debug else 'info'
        level = 'warning' if quiet else default_level
        from cpylog import SimpleLogger
        log = SimpleLogger(level=level)

    if data['NASTRAN_BDF_OUT']:
        nastran_filename_out = data['NASTRAN_BDF_OUT']
    else:
        log.info(f"NASTRAN_BDF_OUT wasn't specified; using {nastran_filename_out!r}")

    size = 8
    if data['--large']:
        size = 16

    xform = ['--xform']
    abaqus_to_nastran_filename(
        abaqus_inp_filename, nastran_filename_out,
        encoding=encoding, size=size,
        xform=xform, log=log)

if __name__ == '__main__':
    cmd_abaqus_to_nastran()
