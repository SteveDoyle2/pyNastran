from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.optimization import DRESP1

def constrain_solid_stress_from_properties(model: BDF, pids: List[int],
                                          constraint_id: int, dresp_id: int, uallow: float):
    """adds VM stress constraints for solids"""
    skip_properties = [
        'PROD', 'PELAS', 'PDAMP',
        'PBAR', 'PBEAM', 'PBARL', 'PBEAML', 'PBCOMP', 'PBEND',
        'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
    ]
    response_type = 'PSOLID'
    for pid in pids:
        prop = model.properties[pid]
        if prop.type in skip_properties:
            continue
        if prop.type == 'PSOLID':
            label = f'p{pid}'[:8]
            property_type = 'PSOLID'
            atta = 13 # von mises stress
            region = None
            attb = None
            atti = [pid]
            model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)
            model.add_dconstr(constraint_id, dresp_id, uid=uallow)
        else:
            raise NotImplementedError(prop)


def constrain_shell_stress_from_properties(model: BDF, pids: List[int],
                                           constraint_id: int, uallow: float,
                                           dresp_id: int=1):
    """adds VM stress constraints for shells"""
    skip_properties = [
        'PROD', 'PELAS', 'PDAMP',
        'PBAR', 'PBEAM', 'PBARL', 'PBEAML', 'PBCOMP', 'PBEND',
        'PSOLID', 'PLSOLID',
        'PCOMP', 'PCOMPG', 'PSHEAR',
    ]
    response_type = 'STRESS'
    for pid in pids:
        prop = model.properties[pid]
        if prop.type in skip_properties:
            continue
        if prop.type == 'PSHELL':
            property_type = 'PSHELL'
            atta = [9, 17] # von mises upper/lower surface stress
            region = None
            attb = None
            atti = [pid]

            label = f'o{pid}U'[:8]
            model.add_dresp1(dresp_id, label, response_type, property_type, region, atta[0], attb, atti)
            model.add_dconstr(constraint_id, dresp_id, uid=uallow)
            dresp_id += 1

            label = f'o{pid}L'[:8]
            model.add_dresp1(dresp_id, label, response_type, property_type, region, atta[1], attb, atti)
            model.add_dconstr(constraint_id, dresp_id, uid=uallow)
            dresp_id += 1

        #elif prop.type == 'PCOMP':
            #label = 'resp2'
            #property_type = 'PCOMP'
            #layer = 4
            #atta = 9 # von mises upper surface stress
            #region = None
            #attb = layer
            #atti = [pid]
            #DRESP1(dresp_id, label, response_type, property_type, region, atta, attb, atti)
        else:
            raise NotImplementedError(prop)

    DRESP1(dresp_id, label, response_type, property_type, region, atta, attb, atti,
           comment='')

def solid_topology_optimization(model: BDF, pids_to_optimize: List[int], xinit=0.5, uallow=1.e20):
    model.sol = 200
    eids_dict = model.get_element_ids_dict_with_pids(pids_to_optimize)
    pid = max(model.properties) + 1
    model.uncross_reference()

    response_type = 'STRESS'
    dresp_id = 1
    constraint_id = 1
    for pid_old, eids in sorted(eids_dict.items()):
        prop_old = model.properties[pid_old]
        if prop_old.type != 'PSOLID':
            continue
        assert prop_old.type == 'PSOLID', prop_old.get_stats()

        opt_id = 101
        label = f'top{pid}'[:8]
        prop_type = 'PSOLID'
        model.add_topvar(opt_id, label, prop_type, xinit, pid, xlb=0.001, delxv=0.2, power=3)

        label = f'p{pid}'[:8]
        property_type = 'PSOLID'
        atta = 13 # von mises stress
        region = None
        attb = None
        atti = [pid]
        model.add_dresp1(dresp_id, label, response_type, property_type, region, atta, attb, atti)
        model.add_dconstr(constraint_id, dresp_id, uid=uallow)
        dresp_id += 1

        #for eid in eids:
            #elem = model.elements[eid]
            #assert elem.type in ['CTETRA', 'CHEXA', 'CPENTA', 'CPYRAM'], elem.get_stats()
            #elem.pid = pid
            #model.properties[pid] = deepcopy(prop_old)

            #label = f'e{eid}'
            #response_type = 'STRESS'
            #property_type = 'PSHELL'
            #atta = 9 # von mises upper surface stress
            #region = None
            #attb = None
            #atti = [pid]
            #DRESP1(dresp_id, label, response_type, property_type, region, atta, attb, atti)
            #model.add_dresp1(dresp_id, label, response_type, property_type, region,
                             #atta, attb, atti, validate=True, comment='')
            #pid += 1

    #model.properties[pid]
