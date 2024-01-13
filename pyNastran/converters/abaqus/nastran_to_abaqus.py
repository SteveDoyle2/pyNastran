from collections import defaultdict
from typing import Any

import numpy as np
from pyNastran.bdf.bdf import read_bdf, BDF
from pyNastran.converters.abaqus.abaqus_cards import Part, Material, Step, ShellSection, SolidSection

from pyNastran.converters.abaqus.abaqus import Abaqus

def nastran_to_abaqus_filename(bdf_filename: str, abaqus_inp_filename: str) -> None:
    nastran_model = read_bdf(bdf_filename)
    log = nastran_model.log
    model = Abaqus(log=None, debug=True)
    name = 'model'
    nids, nodes = _get_nodes(nastran_model, fdtype='float32')

    node_sets = {}
    element_sets = {}

    pid_to_name_map, element_sets_temp, shell_sections, solid_sections = _get_properties(
        nastran_model)
    element_types, element_sets_temp = _get_elements(nastran_model, element_sets_temp)

    for pid, pid_str in pid_to_name_map.items():
        eids = element_sets_temp[pid]
        element_sets[pid_str] = np.array(eids, dtype='int32')
    del pid_to_name_map
    del element_sets_temp
    cloads = _process_constraints(nastran_model, node_sets)

    solid_sections = []
    shell_sections = []
    part = Part(name, nids, nodes, element_types, node_sets, element_sets,
                solid_sections, shell_sections, log=log)
    model.parts = {
        'model': part,
    }

    _get_materials(nastran_model, model)

    name = 'static_step'
    boundaries = []
    node_output = []
    element_output = []
    static_step = Step(name, boundaries, node_output, element_output,
                       cloads=cloads, dloads=[], surfaces=[], is_nlgeom=False)
    model.steps = [static_step]
    model.write(abaqus_inp_filename, is_2d=False)

def _get_nodes(nastran_model: BDF, fdtype: str='float32') -> tuple[np.ndarray, np.ndarray]:
    """creates nodes"""
    nids = []
    nodes = []
    for nid, node in nastran_model.nodes.items():
        xyz = node.get_position()
        nids.append(nid)
        nodes.append(xyz)
    nodes = np.array(nodes, dtype=fdtype)
    return nids, nodes

def _get_properties(nastran_model: BDF) -> tuple[Any, Any, list[ShellSection], list[SolidSection]]:
    pid_to_name_map = {}
    element_sets_temp = {}
    shell_sections = []
    solid_sections = []

    log = nastran_model.log
    for pid, prop in nastran_model.properties.items():
        pid_to_name_map[pid] = f'{prop.type}_{pid}'  # PSHELL_20
        element_sets_temp[pid] = []  # 20
        if prop.type == 'PSHELL':
            mid = prop.mid1
            elset = None
            material_name = f'{prop.mid1_ref.type}_{mid}'
            thickness = prop.t
            shell_section = ShellSection(material_name, elset, thickness, log)
            shell_sections.append(shell_section)
        elif prop.type == 'PSOLID':
            material_name = f'{prop.mid_ref.type}_{mid}'
            elset = None
            thickness = None
            solid_section = SolidSection(material_name, elset, thickness, log)
            solid_sections.append(solid_section)
        else:
            print(prop)
        #elif prop.type == 'PSHELL':
    return pid_to_name_map, element_sets_temp, shell_sections, solid_sections

def _get_elements(nastran_model: BDF, element_sets_temp: dict[str, list[int]]):
    ctria3s = []
    cquad4s = []
    ctetra4s = []
    ctetra10s = []
    chexa8s = []
    chexa20s = []
    element_types = {
        # elements, element_set_name
        'cpe3' : (ctria3s, ''), # CTRIA3
        'cpe4' : (cquad4s, ''),
        'c3d10h' : (ctetra10s, ''),
        'c3d4r' : (ctetra4s, ''),
        'c3d8r' : (chexa8s, ''),
        'c3d20r' : (chexa20s, ''),
    }
    for eid, elem in nastran_model.elements.items():
        pid = elem.pid
        nidsi = elem.nodes
        if elem.type == 'CTRIA3':
            ctria3s.append([eid] + nidsi)
        elif elem.type == 'CQUAD4':
            cquad4s.append([eid] + nidsi)
        elif elem.type == 'CTETRA4':
            ctetra4s.append([eid] + nidsi)
        elif elem.type == 'CTETRA10':
            ctetra10s.append([eid] + nidsi)
        elif elem.type == 'CHEXA8':
            chexa8s.append([eid] + nidsi)
        elif elem.type == 'CHEXA20':
            chexa20s.append([eid] + nidsi)
        else:
            print(elem)
        element_sets_temp[pid].append(eid)
    assert len(element_sets_temp)
    assert len(element_types)

    return element_types, element_sets_temp

def _process_constraints(nastran_model: BDF, node_sets: dict[str, np.ndarray]) -> dict[str, Any]:
    """creates node_sets"""
    all_cloads = []
    spc_dict = defaultdict(list)
    for subcase_id, subcase in nastran_model.subcases.items():
        if subcase_id == 0:
            continue
        spc_id = subcase['SPC'][0]
        load_id = subcase['LOAD'][0]
        spcs = nastran_model.get_reduced_spcs(spc_id, consider_spcadd=True, stop_on_failure=True)
        for spc in spcs:
            if spc.type == 'SPC1':
                name = f'subcase={subcase_id}_spc={spc_id}_comp={spc.components}'
                spc_dict[name].extend(spc.nodes)
                x = 1
            else:
                print(spc)
        _process_loads(nastran_model, subcase_id, load_id, all_cloads)

    for key, mylist in spc_dict.items():
        node_sets[key] = np.array(mylist, dtype='int32')
    return all_cloads

def _process_loads(nastran_model: BDF,
                   subcase_id: int, load_id: int, all_cloads: list[Any]) -> None:
    #if load_id:
        #return
    name = f'subcase={subcase_id}_LOAD={load_id}'
    cloads = []
    loads, scale_factors, is_grav = nastran_model.get_reduced_loads(
        load_id, scale=1., consider_load_combinations=True,
        skip_scale_factor0=False, stop_on_failure=True, msg='')
    for scale, load in zip(scale_factors, loads):
        if load.type == 'FORCE':
            nid = load.node
            if load.cid > 0:
                print(load.get_stats())
            else:
                for i, mag in enumerate(load.scaled_vector):
                    if abs(mag) == 0.0:
                        continue
                    dof = i + 1
                    cload = [nid, dof, scale * mag]
                    cloads.append(cload)
        else:
            print(load)
    if cloads:
        all_cloads.append(cloads)


def _get_materials(nastran_model: BDF, model: Abaqus):
    for mid, mat in nastran_model.materials.items():
        name = f'{mat.type}_mid{mid}'
        density = mat.rho if mat.rho else 0.0
        sections = {
            'elastic': 'cat',
        }
        material = Material(name, sections, density=density,
                            is_elastic=True, ndepvars=None, ndelete=None)
        model.materials[name] = material
