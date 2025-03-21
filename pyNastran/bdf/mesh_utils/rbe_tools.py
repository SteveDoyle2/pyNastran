import warnings
from collections import defaultdict
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.elements.mass import CONM2
from pyNastran.bdf.cards.elements.rigid import RBE2, RBE3
from pyNastran.bdf.mesh_utils.mass_properties import increment_inertia, transform_inertia

def rbe3_to_rbe2(model: BDF, eids_to_fix: list[int]) -> None:
    rigid_elements2 = {}
    for eid, elem in model.rigid_elements.items():
        if eid not in eids_to_fix:
            warnings.warn(f'skipping eid={eid} because its not an element to fix\n{str(elem)}')
            rigid_elements2[eid] = elem
            continue
        if elem.type not in {'RBE3'}:
            warnings.warn(f'skipping eid={eid} because its not an RBE3\n{str(elem)}')
            rigid_elements2[eid] = elem
            continue

        # ---RBE3---
        #   Gijs   : [[4407195, 4407212, 4407280, 4407288]]
        #   independent_nodes : [4407195, 4407212, 4407280, 4407288]
        #
        #   Cmi    : []
        #   Gmi    : []
        #   Gmi_node_ids : []
        #
        #   comps  : ['123']
        #   dependent_nodes : [10913]
        #   eid    : 11913
        #   nodes  : [4407195, 4407212, 4407280, 4407288, 10913]
        #   ref_grid_id : 10913
        #   refc   : '123456'
        #   refgrid : 10913
        #   weights : [1.0]
        #   wt_cg_groups : [(1.0, '123', [4407195, 4407212, 4407280, 4407288])]
        assert len(elem.Gmi) == 0, elem.Gmi
        assert elem.comps == ['123'], elem.comps
        assert elem.refc == '123456', elem.refc
        assert elem.weights == [1.0], elem.weights
        assert len(elem.comps) == 1
        ind = elem.independent_nodes
        dep = elem.dependent_nodes
        gn = elem.refgrid  # ind; was dep
        cm = '123456'
        Gmi = elem.independent_nodes
        elem = RBE2(eid, gn, cm, Gmi,
                    tref=elem.tref, alpha=elem.alpha,
                    comment=elem.comment.strip())
        rigid_elements2[eid] = elem
    model.rigid_elements = rigid_elements2


def rbe2_to_rbe3(model: BDF, eids_to_fix: list[int]) -> None:
    rigid_elements2 = {}
    for eid, elem in model.rigid_elements.items():
        if eid not in eids_to_fix:
            warnings.warn(f'skipping eid={eid} because its not an element to fix\n{str(elem)}')
            rigid_elements2[eid] = elem
            continue
        if elem.type not in {'RBE2'}:
            warnings.warn(f'skipping eid={eid} because its not an RBE2\n{str(elem)}')
            rigid_elements2[eid] = elem
            continue

        # ---RBE2---
        #   Gmi    : [1, 2]
        #   cm     : '123456'
        #   gn     : 10
        assert len(elem.Gmi) > 0, elem.Gmi
        assert elem.cm == '123456', elem.cm
        # ind = elem.independent_nodes
        # dep = elem.dependent_nodes
        # gn = elem.refgrid  # ind; was dep
        # cm = '123456'
        # Gmi = elem.independent_nodes
        # elem = RBE2(eid, gn, cm, Gmi,
        #             tref=elem.tref, alpha=elem.alpha,
        #             comment=elem.comment.strip())
        refgrid = elem.gn
        refc = elem.cm
        weights = [1.0]
        comps = ['123']
        Gijs = [elem.dependent_nodes]
        elem = RBE3(
            eid, refgrid, refc, weights, comps, Gijs,
            alpha=elem.alpha, tref=elem.tref,
            comment=elem.comment)
        #print(elem)
        rigid_elements2[eid] = elem
    model.rigid_elements = rigid_elements2


def merge_rbe2(model: BDF, rbe_eids_to_fix: list[int]) -> None:
    log = model.log
    eid_to_dependent_nodes = {}
    # eid_to_eids = {}
    dependent_nid_to_eids_map = defaultdict(set)
    nrigid = 0
    for eid, elem in model.rigid_elements.items():
        if eid not in rbe_eids_to_fix:
            log.warning(f'skipping eid={eid} because its not an element to fix\n{str(elem)}')
        if elem.type not in {'RBE2'}:
            log.warning(f'skipping eid={eid} because its not an RBE2\n{str(elem)}')
            continue
        dependent_nodes = elem.dependent_nodes
        eid_to_dependent_nodes[eid] = dependent_nodes
        for nid in dependent_nodes:
            dependent_nid_to_eids_map[nid].add(eid)
        nrigid += 1
    dependent_nid_to_eids_map = dict(dependent_nid_to_eids_map)
    #print(f'dependent_nid_to_eids_map = {dependent_nid_to_eids_map}')

    eid_to_nids_map = get_eid_to_nid_map(eid_to_dependent_nodes, dependent_nid_to_eids_map)
    eid_to_nids_map = get_eid_to_nid_map(eid_to_nids_map, dependent_nid_to_eids_map)
    eid_to_nids_map = get_eid_to_nid_map(eid_to_nids_map, dependent_nid_to_eids_map)
    #print(f'eid_to_nids_map = {eid_to_nids_map}')

    nids_to_eids_map = defaultdict(list)
    for eid, nids in eid_to_nids_map.items():
        nids_list = list(nids)
        nids_list.sort()
        nids_to_eids_map[tuple(nids_list)].append(eid)
    #print(f'nids_to_eids_map = {set(nids_to_eids_map)}')

    mass_nid_to_elem_map: dict[int, CONM2] = {}
    for eid, mass_elem in model.masses.items():
        if mass_elem.type != 'CONM2':
            log.warning(f'skipping eid={eid} because its not a CONM2\n{str(mass_elem)}')
            continue
        #print(elem.get_stats())
        mass_nid_to_elem_map[mass_elem.nid] = mass_elem
        mass_elem.cross_reference(model)

    for dep_nids_tuple, eids in nids_to_eids_map.items():
        eid0 = eids[0]
        elem0: RBE2 = model.rigid_elements[eid0]
        indep_nid0 = elem0.gn

        log.info(f'eids={eids}; dep_nids={dep_nids_tuple}\n')
        if len(eids) == 1:
            log.debug('skipping\n' + str(elem0))
            #log.debug('\n' + str(mass_elem0))
            continue

        mass_elem0 = mass_nid_to_elem_map[indep_nid0]
        mass_node0 = mass_elem0.nid_ref
        log.info('baseline:\n' + str(elem0))
        log.info('\n' + str(mass_elem0))

        mass_comments = []
        reference_point = np.zeros(3, dtype='float64')
        mass_cg = np.zeros(3, dtype='float64')
        mass = 0.

        # [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
        inertia = np.zeros(6, dtype='float64')
        for eid in eids:
            rigid_elemi: RBE2 = model.rigid_elements[eid]
            indep_nid = rigid_elemi.gn
            mass_elem: CONM2 = mass_nid_to_elem_map[indep_nid]
            massi = mass_elem.mass
            mass_node = mass_elem.nid_ref
            centroidi = mass_elem.X + mass_node.get_position()
            log.debug(f'massi={massi}; centroidi={centroidi}')
            I11, I21, I22, I31, I32, I33 = mass_elem.I
            mass = increment_inertia(centroidi, reference_point, massi, mass, mass_cg, inertia)
            log.debug(f'mass_cg={mass_cg}')
            #[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
            inertia[0] += I11  # Ixx
            inertia[1] += I22  # Iyy
            inertia[2] += I33  # Izz
            inertia[3] += I21  # Ixy
            inertia[4] += I31  # Ixz
            inertia[5] += I32  # Iyz
            if mass_elem.comment:
                mass_comments.append(mass_elem.comment)

        # cleanup old elements
        log.warning('removing')
        for eid in eids[1:]:
            rigid_elemi = model.rigid_elements[eid]
            del model.rigid_elements[eid]
            indep_nid = rigid_elemi.gn
            mass_elemi = mass_nid_to_elem_map[indep_nid]
            mass_eid = mass_elemi.eid
            log.warning('\n' + str(rigid_elemi))
            log.warning('\n' + str(mass_elemi))
            del model.masses[mass_eid]

        # update the CONM2 mass
        log.debug(f'mass = {mass}')
        if mass:
            cg = mass_cg / mass
        else:
            cg = mass_cg
        log.debug(f'cg = {cg}')

        inertia2 = transform_inertia(mass, cg, reference_point, cg, inertia)
        Ixx, Iyy, Izz, Ixy, Ixz, Iyz = inertia2
        mass_comments.append(f'Ixx={Ixx:g} Iyy={Iyy:g} Izz={Izz:g}\nIxy={Ixy:g} Ixz={Ixz:g} Iyz={Iyz:g}\n')
        if mass_comments:
            mass_elem0.comment = ''.join(mass_comments[::-1])

        mass_elem0.X = np.zeros(3, dtype='float64')
        mass_elem0.mass = mass
        mass_node0.comment = f'xyz=[{cg[0]:g},{cg[1]:g},{cg[2]:g}]'
        mass_node0.xyz = cg
        # FEMAP Mass Table: 'Ixx, Iyy, Izz, Ixy, Iyz, Izx'
        mass_elem0.I = [Ixx, Ixy, Iyy, Ixz, Iyz, Izz]

        elem0.Gmi = list(dep_nids_tuple)
        log.info('\n' + str(elem0))
        log.info('\n' + str(mass_elem0))
        log.debug('---------------------')
    return

def get_eid_to_nid_map(
        eid_to_nids_map: dict[int, set[int]],
        nid_to_eids_map: dict[int, set[int]],):
    eid_to_nids_map2 = defaultdict(set)
    for eid, nids in eid_to_nids_map.items():
        for nid in nids:
            eid_to_nids_map2[eid].add(nid)
            eids2 = nid_to_eids_map[nid]
            for eid2 in eids2:
                nids2 = eid_to_nids_map[eid2]
                eid_to_nids_map2[eid].update(nids2)
    return dict(eid_to_nids_map2)
