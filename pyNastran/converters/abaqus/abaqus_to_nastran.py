from collections import defaultdict
from typing import List, Dict, Any

import numpy as np
from pyNastran.bdf.bdf import read_bdf, BDF, CaseControlDeck
from pyNastran.converters.abaqus.abaqus import (
    Abaqus, Part, Material, Step,
    ShellSection, SolidSection,
    read_abaqus, get_nodes_nnodes_nelements)

def nastran_to_abaqus_filename(bdf_filename: str, abaqus_inp_filename: str):
    nastran_model = read_bdf(bdf_filename)
    log = nastran_model.log
    model = Abaqus(log=None, debug=True)
    name = 'model'
    nids = []
    nodes = []
    for nid, node in nastran_model.nodes.items():
        xyz = node.get_position()
        nids.append(nid)
        nodes.append(xyz)
    nodes = np.array(nodes, dtype='float32')

    node_sets = {}
    element_sets = {}

    ctria3s = []
    cquad4s = []
    ctetra4s = []
    ctetra10s = []
    chexa8s = []
    chexa20s = []
    element_types = {
        'cpe3' : ctria3s, # CTRIA3
        'cpe4' : cquad4s,
        'c3d10h' : ctetra10s,
        'c3d4r' : ctetra4s,
        'c3d8r' : chexa8s,
        'c3d20r' : chexa20s,
    }
    pid_to_name_map = {}
    element_sets_temp = {}
    shell_sections = []
    for pid, prop in nastran_model.properties.items():
        pid_to_name_map[pid] = f'{prop.type}_{pid}'  # PSHELL_20
        element_sets_temp[pid] = []  # 20
        if prop.type == 'PSHELL':
            mid = prop.mid1
            material_name = f'{prop.mid1_ref.type}_{mid}'
            thickness = prop.t
            shell_section = ShellSection(material_name, thickness, log)
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

    for eid, elem in nastran_model.elements.items():
        pid = elem.pid
        nidsi = elem.nodes
        if elem.type in ['CTRIA3']:
            ctria3s.append([eid] + nidsi)
        elif elem.type in ['CQUAD4']:
            cquad4s.append([eid] + nidsi)
        elif elem.type in ['CTETRA4']:
            ctetra4s.append([eid] + nidsi)
        elif elem.type in ['CTETRA10']:
            ctetra10s.append([eid] + nidsi)
        elif elem.type in ['CHEXA8']:
            chexa8s.append([eid] + nidsi)
        elif elem.type in ['CHEXA20']:
            chexa20s.append([eid] + nidsi)
        else:
            print(elem)
        element_sets_temp[pid].append(eid)

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

    for mid, mat in nastran_model.materials.items():
        name = f'{mat.type}_mid{mid}'
        density = mat.rho if mat.rho else 0.0
        sections = {
            'elastic': 'cat',
        }
        material = Material(name, sections, density=density,
                            is_elastic=True, ndepvars=None, ndelete=None)
        model.materials[name] = material

    name = 'static_step'
    boundaries = []
    outputs = []
    static_step = Step(name, boundaries, outputs, cloads=cloads, is_nlgeom=False)
    model.steps = [static_step]
    model.write(abaqus_inp_filename, is_2d=False)


def _process_constraints(nastran_model: BDF, node_sets: Dict[str, np.ndarray]) -> Dict[str, Any]:
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
                   subcase_id: int, load_id: int, all_cloads: List[Any]):
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

def _add_part_to_nastran(nastran_model: BDF, elements, pid: int, nid_offset: int) -> None:
    for etype, eids_nids in elements.element_types.items():
        eids, part_nids = eids_nids
        if eids is None and part_nids is None:
            continue

        if nid_offset > 0:
            # don't use += or it's an inplace operation
            part_nids = part_nids + nid_offset

        if etype == 'cpe4':
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_cquad4(eid, pid, nids[1:], theta_mcid=0.0, zoffset=0., tflag=0,
                                         T1=None, T2=None, T3=None, T4=None, comment='')
        elif etype == 's8r':
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_cquad8(eid, pid, nids[1:], theta_mcid=0.0, zoffset=0., tflag=0,
                                         T1=None, T2=None, T3=None, T4=None, comment='')
        elif etype == 'c3d4':
            for eid, nids in zip(eids, part_nids):
                nastran_model.add_ctetra(eid, pid, nids[1:], comment='')
        else:
            raise NotImplementedError(etype)

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


def _create_nastran_nodes_elements(model: Abaqus, nastran_model: BDF) -> None:
    log = model.log
    nid_offset = 0

    pid = 1
    if model.nids is not None and len(model.nids):
        nnodesi = model.nodes.shape[0]
        elements = model.elements
        _add_part_to_nastran(nastran_model, elements, pid, nid_offset)
        nid_offset += nnodesi

    for unused_part_name, part in model.parts.items():
        #log.info('part_name = %r' % unused_part_name)
        nnodesi = part.nodes.shape[0]
        #nidsi = part.nids

        elements = part.elements
        _add_part_to_nastran(nastran_model, elements, pid, nid_offset)
        #nids.append(nidsi)

        nid_offset += nnodesi
        for shell_section in part.shell_sections:
            log.info('shell')
        for shell_section in part.shell_sections:
            log.info('solid')
    #nids = np.hstack(nids)

    pid = 1
    mid = 1
    for solid_section in model.solid_sections:
        #print(solid_section)
        mat_name = solid_section.material_name
        mat = model.materials[mat_name]
        element_set = solid_section.elset
        log.warning(f'element_set = {element_set}')
        #print(mat)
        nastran_model.add_psolid(pid, mid, cordm=0,
                                 integ=None, stress=None, isop=None, fctn='SMECH', comment='')
        G = None
        elastic = mat.sections['elastic']
        E = elastic[0]
        nu = elastic[1]
        nastran_model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0,
                               ge=0.0, St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        log.info('shell2')
        mid += 1
    for shell_section in model.shell_sections:
        #print(shell_section)
        mat_name = shell_section.material_name
        mat = model.materials[mat_name]
        #print(mat)
        t = shell_section.thickness
        nastran_model.add_pshell(pid, mid1=mid, t=t, mid2=mid, twelveIt3=1.0, mid3=None, tst=0.833333,
                                 nsm=0.0, z1=None, z2=None, mid4=None, comment='')
        G = None
        elastic = mat.sections['elastic']
        E = elastic[0]
        nu = elastic[1]
        nastran_model.add_mat1(mid, E, G, nu, rho=0.0, a=0.0, tref=0.0,
                               ge=0.0, St=0.0, Sc=0.0, Ss=0.0, mcsid=0, comment='')
        log.info('shell2')
        mid += 1

def abaqus_to_nastran_filename(abaqus_inp_filename: str,
                               nastran_filename_out: str,
                               log=None) -> BDF:
    if isinstance(abaqus_inp_filename, Abaqus):
        model = abaqus_inp_filename
    else:
        model = read_abaqus(abaqus_inp_filename, log=log, debug=True)
    nnodes, nids, nodes, nelements = get_nodes_nnodes_nelements(
        model, stop_for_no_elements=True)
    assert nnodes > 0, nnodes
    assert nelements > 0, nelements

    nastran_model = BDF(debug=True, log=log, mode='msc')
    for nid, xyz in zip(nids, nodes):
        nastran_model.add_grid(nid, xyz)
    _create_nastran_nodes_elements(model, nastran_model)
    #for step in model.steps:
        #print(step)

    _create_nastran_loads(model, nastran_model)
    nastran_model.write_bdf(nastran_filename_out)
    x = 1
    return nastran_model

def _create_nastran_loads(model: Abaqus, nastran_model: BDF):
    def xyz():
        return np.zeros(3, dtype='float32')

    if nastran_model.case_control_deck is None:
        nastran_model.case_control_deck = CaseControlDeck([], log=nastran_model.log)

    for istep, step in enumerate(model.steps):
        subcase_id = istep + 1
        load_id = subcase_id
        if subcase_id in nastran_model.subcases:
            subcase = nastran_model.subcases[subcase_id]
        else:
            subcase = nastran_model.case_control_deck.create_new_subcase(subcase_id)
        for cload in step.cloads:
            subcase.add('LOAD', load_id, [], 'STRESS-type')
            #subcase['LOAD'] = load_id
            forces = defaultdict(xyz)
            moments = defaultdict(xyz)
            for cloadi in cload:
                nid, dof, mag = cloadi
                assert dof in [1, 2, 3], cload
                if dof in {1, 2, 3}:
                    forces[nid][dof - 1] = mag
                elif dof in {4, 5, 6}:
                    moments[nid][dof - 4] = mag
                else:
                    raise NotImplementedError(cloadi)

            mag = 1.0
            if len(forces):
                for nid, xyz in forces.items():
                    nastran_model.add_force(load_id, nid, mag, xyz, cid=0, comment='')
            if len(moments):
                for nid, xyz in moments.items():
                    nastran_model.add_moment(load_id, nid, mag, xyz, cid=0, comment='')


            #print(step.cloads)
        #step.cloads


if __name__ == '__main__':   # pragma: no cover
    nastran_filename = r'C:\NASA\m4\formats\git\pyNastran\models\plate\plate.bdf'
    abaqus_inp_filename = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\converters\abaqus\plate.inp'
    nastran_to_abaqus_filename(nastran_filename, abaqus_inp_filename)

    #model = read_abaqus(abaqus_filename, debug=True)
    nastran_filename_out = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\converters\abaqus\plate2.bdf'
    abaqus_to_nastran_filename(abaqus_inp_filename, nastran_filename_out)
    x = 1
    #nastran_filename = 'junk.bdf'

    #abaqus_to_nastran_filename(abaqus_filename, nastran_filename)
    #nastran_filename = 'junk.bdf'
