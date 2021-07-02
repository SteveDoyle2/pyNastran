from typing import Tuple, List, Set, Union
import numpy as np
from pyNastran.bdf.bdf import BDF, PLOAD2

def remove_missing_loads(model: BDF) -> None:
    """
    removes loads where:
     - PLOAD2s missing elements
       - eids=1,THRU,10 with 5-10 missing becomes eids=1,THRU,4 (on separate cards)
    """
    all_eids = set(list(model.elements))
    loads2_dict = {}
    for sid, loads in model.loads.items():
        loads2 = []
        pload2s = []
        for iload, load in enumerate(loads):
            if load.type == 'PLOAD2':
                load.eids_ref = None
                pload2s.append(load)
            else:
                loads2.append(load)

        # fix the PLOAD2s dropping elements
        pload2s_new = []
        for pload2 in pload2s:
            pressure = pload2.pressure
            eids_to_add = [eid for eid in pload2.eids
                           if eid in all_eids]
            for eid in eids_to_add:
                pload2_new = PLOAD2(sid, pressure, [eid])
                pload2s_new.append(pload2_new)
        loads2 = loads2 + pload2s_new
        loads2_dict[sid] = loads2
    model.loads = loads2_dict


def remove_mpc_chain(model: BDF, nodes: List[int]) -> List[int]:
    """deletes any MPC connected to a given set of nodes"""
    mpcs_dict = {}
    nodes_all = []
    for mpc_id, mpcs in model.mpcs.items():
        # identify the indices of the mpcs to delete
        mpcs2 = []
        for impc, mpc in enumerate(mpcs):
            nodes_mpc = mpc.dependent_nodes + mpc.independent_nodes
            for node in nodes_mpc:
                if node in nodes:
                    nodes_all.extend(nodes_mpc)
                    break
            else:
                # these mpcs were not filtered
                mpcs2.append(mpc)
                x = 1

        mpcs_dict[mpc_id] = mpcs2
    model.mpcs = mpcs_dict
    return nodes_all

def remove_rbe_chain(model: BDF, nodes: List[int]) -> Tuple[Set[int], List[int]]:
    """deletes any RBE2, RBE3, RBAR, etc. connected to a given set of nodes"""
    nodes_all = set([])
    eids_to_delete = []
    for eid, rbe in model.rigid_elements.items():
        rbe_nodes = rbe.independent_nodes + rbe.dependent_nodes
        break_flag = False
        for independent_node in rbe.independent_nodes:
            if independent_node in rbe_nodes:
                eids_to_delete.append(eid)
                nodes_all.update(rbe_nodes)
                break_flag = True
                break
        if break_flag:
            break
        for dependent_node in rbe.dependent_nodes:
            if dependent_node in rbe_nodes:
                eids_to_delete.append(eid)
                nodes_all.update(rbe_nodes)
                break
    for eid in eids_to_delete:
        del model.rigid_elements[eid]
    return nodes_all, eids_to_delete


def remove_element_chain(model: BDF, nodes_all: Set[int], element_types: Union[str, Set[str]]) -> Tuple[Set[int], List[int]]:
    if isinstance(element_types, str):
        element_types = set([element_types])
    eids_to_delete = []
    nodes_all2 = set([])
    for eid, elem in model.elements.items():
        if elem.type in element_types:
            for nid in elem.nodes:
                if nid in nodes_all:
                    eids_to_delete.append(eid)
                    nodes_all2.update(elem.nodes)
                    break
    for eid in eids_to_delete:
        del model.elements[eid]
    return nodes_all2, eids_to_delete
