from collections import defaultdict
from typing import Dict

import numpy as np
from pyNastran.bdf.bdf import read_bdf, BDF
from pyNastran.converters.abaqus.abaqus import Abaqus, Part, read_abaqus, get_nodes_nnodes_nelements

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
    for pid, prop in nastran_model.properties.items():
        pid_to_name_map[pid] = f'{prop.type}_{pid}'  # PSHELL_20
        element_sets_temp[pid] = []  # 20
        #if prop.type == 'PSHELL':
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
    _process_constraints(nastran_model, node_sets)

    solid_sections = []
    part = Part(name, nids, nodes, element_types, node_sets, element_sets, solid_sections, log)
    model.parts = {
        'model': part,
    }
    model.write(abaqus_inp_filename, is_2d=False)


def _process_constraints(nastran_model: BDF, node_sets: Dict[str, np.ndarray]):
    """creates node_sets"""
    spc_dict = defaultdict(list)
    for subcase_id, subcase in nastran_model.subcases.items():
        if subcase_id == 0:
            continue
        spc_id = subcase['SPC'][0]
        spcs = nastran_model.get_reduced_spcs(spc_id, consider_spcadd=True, stop_on_failure=True)
        for spc in spcs:
            if spc.type == 'SPC1':
                name = f'subcase={subcase_id}_spc={spc_id}_comp={spc.components}'
                spc_dict[name].extend(spc.nodes)
                x = 1
            else:
                print(spc)
                x = 1

    for key, mylist in spc_dict.items():
        node_sets[key] = np.array(mylist, dtype='int32')

def make_nodes_elements(model: Abaqus, nastran_model: BDF) -> None:
    log = model.log
    nid_offset = 0
    for unused_part_name, part in model.parts.items():
        log.info('part_name = %r' % unused_part_name)
        nnodesi = part.nodes.shape[0]
        #nidsi = part.nids

        pid = 1
        for etype, eids_nids in part.element_types.items():
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
            else:
                raise NotImplementedError(etype)
        #nids.append(nidsi)

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

        nid_offset += nnodesi
    #nids = np.hstack(nids)

def abaqus_to_nastran_filename(abaqus_inp_filename: str, nastran_filename_out: str):
    model = read_abaqus(abaqus_inp_filename, log=None, debug=False)
    nnodes, nids, nodes, nelements = get_nodes_nnodes_nelements(model, stop_for_no_elements=True)
    assert nnodes > 0, nnodes
    assert nelements > 0, nelements

    nastran_model = BDF(debug=True, log=None, mode='msc')
    for nid, xyz in zip(nids, nodes):
        grid = nastran_model.add_grid(nid, xyz)
    make_nodes_elements(model, nastran_model)
    nastran_model.write_bdf(nastran_filename_out)
    x = 1


if __name__ == '__main__':
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
