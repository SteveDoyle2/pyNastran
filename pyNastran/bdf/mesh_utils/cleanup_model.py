from copy import deepcopy
from typing import Any
from pyNastran.bdf.bdf import BDF, LOAD, PLOAD2, PLOAD4


def cleanup_model(model: BDF) -> None:
    """
    The model is in in an invalid state (e.g., missing nodes/elements), so
    trace all the following cards and clean them up (e.g., delete SPCs or PLOAD4s
    that reference invalid nodes/elements).

    Fixes:
     - loads:
       - LOAD, PLOAD4
     - constraints
       - SPC

    Assumptions:
     - subcase is correct (you can fix it by hand)
     - elements are correct
    """
    all_nids = set(model.nodes.keys())
    all_eids = set(model.elements.keys())
    _cleanup_loads(model, all_eids, all_nids)
    _cleanup_spcs(model, all_nids)
    model.cross_reference()


def _cleanup_rigid_elements(model: BDF, nids_to_delete: list[int]) -> None:
    rigid_elements2 = {}
    for eid, elem in model.rigid_elements.items():
        if elem.type == 'RBE3':
            Gijs2 = []
            comps2 = []
            weights2 = []
            for weight, comp, Gij in zip(elem.weights, elem.comps, elem.Gijs):
                Gij2 = []
                for nid in Gij:
                    if nid in nids_to_delete:
                        continue
                    Gij2.append(nid)
                if Gij2:
                    Gijs2.append(Gij2)
                    comps2.append(comp)
                    weights2.append(weight)
            if len(weights2) == 0:
                continue
            elem.Gijs = Gijs2
            elem.comps = comps2
            elem.weights = weights2

            # TODO: check the dependent node
        else:
            raise NotImplementedError(elem)
        rigid_elements2[eid] = elem
    model.rigid_elements = rigid_elements2

def _cleanup_pload4(load: PLOAD4,
                    all_eids: list[int],
                    loads2: list[Any]):
    _eids = []
    for eid in load.eids:
        if eid not in all_eids:
            return None
        _eids.append(eid)

    neids = len(_eids)
    if neids == 0:
        return None

    neids_compress = max(_eids) - min(_eids) + 1
    if neids == 1 or neids == neids_compress:
        load.eids = _eids
        loads2.append(load)
    else:
        for eid in _eids:
            load2 = deepcopy(load)
            load2.eids = [eid]
            loads2.append(load2)
    return None

def _cleanup_load(load_combo: LOAD,
                  all_loads2: dict[int, Any],
                  load_combos2: list[LOAD]) -> None:
    assert load_combo.type == 'LOAD', load_combo
    load_ids = []
    scale_factors = []
    for load_id, scale_factor in zip(load_combo.load_ids, load_combo.scale_factors):
        if load_id not in all_loads2:
            return
        load_ids.append(load_id)
        scale_factors.append(scale_factor)
    if len(load_ids) == 0:
        return
    load_combo.load_ids = load_ids
    load_combo.scale_factors = scale_factors
    load_combos2.append(load_combo)

def _cleanup_loads(model: BDF,
                   all_eids: list[int],
                   all_nids: list[int]) -> None:
    """remove loads that aren't referenced"""
    all_loads2 = {}
    for load_id, loads in model.loads.items():
        loads2 = []
        for load in loads:
            load_type = load.type
            if load_type == 'PLOAD4':
                # adding loads handled within the function
                load = _cleanup_pload4(load, all_eids, loads2)
                continue
            elif load_type in {'FORCE', 'MOMENT'}:
                if load.node not in all_nids:
                    continue
            else:
                raise NotImplementedError(load)
            loads2.append(load)
        if loads2:
            all_loads2[load_id] = loads2

    load_combinations2 = {}
    for load_id, load_combos in model.load_combinations.items():
        load_combos2 = []
        for load_combo in load_combos:
            if load_combo.type == 'LOAD':
                _cleanup_load(load_combo, all_loads2, load_combos2)
                continue
            else:
                raise NotImplementedError(load_combo)
            load_combos2.append(load_combo)

        if load_combos2:
            load_combinations2[load_id] = load_combos2
    model.loads = all_loads2
    model.load_combinations = load_combinations2


def _cleanup_spcs(model: BDF, all_nids) -> None:
    """remove spcs that aren't referenced"""
    all_spcs2 = {}
    for spc_id, spcs in model.spcs.items():
        spcs2 = []
        for spc in spcs:
            #if spc.type == 'SPC1':
                #asdf
            if spc.type == 'SPC':
                nodes2 = []
                for nid in spc.nodes:
                    if nid in all_nids:
                        nodes2.append(nid)
                if len(nodes2) == 0:
                    continue
                else:
                    spc.nodes = nodes2
            else:
                raise NotImplementedError(spc)
        if spcs2:
            all_spcs2[spc_id] = spcs2
    model.spcs = all_spcs2

    x = 1

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


def remove_mpc_chain(model: BDF, nodes: list[int]) -> list[int]:
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

def remove_rbe_chain(model: BDF, nodes: list[int]) -> tuple[set[int], list[int]]:
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


def remove_element_chain(model: BDF,
                         nodes_all: set[int],
                         element_types: str | set[str]) -> tuple[set[int], list[int]]:
    if isinstance(element_types, str):
        element_types = {element_types}
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
