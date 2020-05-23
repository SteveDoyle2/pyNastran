from collections import defaultdict
from typing import Dict

import numpy as np
from pyNastran.bdf.bdf import read_bdf, BDF
from pyNastran.converters.abaqus.abaqus import Abaqus, Part

def nastran_to_abaqus_filename(bdf_filename: str, abaqus_filename: str):
    nastran_model = read_bdf(bdf_filename)
    log = nastran_model.log
    model = Abaqus(log=None, debug=True)
    name = 'model'
    nids = []
    nodes = []
    for nid, node in nastran_model.nodes.items():
        nids.append(nid)
        nodes.append(node.get_position())
    node_sets = {}
    #element_sets = []

    ctria3s = []
    cquad4s = []
    ctetra4s = []
    ctetra10s = []
    chexa8s = []
    chexa20s = []
    element_types = {
        'cpe3' : ctria3s, # CTRIA3
        'cpe4' : cquad4s,
    }
    pid_to_name_map = {}
    element_sets_temp = {}
    element_sets = {}
    for pid, prop in nastran_model.properties.items():
        pid_to_name_map[pid] = f'{prop.type}_{pid}'  # PSHELL_20
        element_sets_temp[pid] = []  # 20
        #if prop.type == 'PSHELL':
        #elif prop.type == 'PSHELL':

    for eid, elem in nastran_model.elements.items():
        pid = elem.pid
        nids = elem.nodes
        if elem.type in ['CTRIA3']:
            ctria3s.append([eid] + nids)
        elif elem.type in ['CQUAD4']:
            cquad4s.append([eid] + nids)
        elif elem.type in ['CTETRA4']:
            ctetra4s.append([eid] + nids)
        elif elem.type in ['CTETRA10']:
            ctetra10s.append([eid] + nids)
        elif elem.type in ['CHEXA8']:
            chexa8s.append([eid] + nids)
        elif elem.type in ['CHEXA20']:
            chexa20s.append([eid] + nids)
        else:
            print(elem)
        element_sets_temp[pid].append(eid)

    for pid, pid_str in pid_to_name_map.items():
        eids = element_sets_temp[pid]
        element_sets[pid_str] = np.array(eids, dtype='int32')
    del pid_to_name_map
    del element_sets_temp
    spc_dict = _process_constraints(nastran_model, node_sets)

    solid_sections = []
    part = Part(name, nids, nodes, element_types, node_sets, element_sets, solid_sections, log)
    model.parts = {
        'model': part,
    }
    model.write(abaqus_filename, is_2d=False)


def _process_constraints(nastran_model: BDF, node_sets: Dict[str, np.ndarray]):
    spc_dict = defaultdict(list)
    for subcase_id, subcase in nastran_model.subcases.items():
        if subcase_id == 0:
            continue
        spc_id = subcase['SPC'][0]
        print('spc =', spc_id)
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

def abaqus_to_nastran_filename(abaqus_filename: str, nastran_filename: str):
    pass

if __name__ == '__main__':
    nastran_filename = r'C:\NASA\m4\formats\git\pyNastran\models\plate\plate.bdf'
    abaqus_filename = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\converters\abaqus\plate.inp'
    nastran_to_abaqus_filename(nastran_filename, abaqus_filename)
    x = 1
    #nastran_filename = 'junk.bdf'

    #abaqus_to_nastran_filename(abaqus_filename, nastran_filename)
    #nastran_filename = 'junk.bdf'
