"""
defines some methods for cleaning up a model
 - remove_unassociated_nodes(...)
 - remove_unassociated_properties(...)
 - remove_unused_materials(...)
"""
from __future__ import print_function
from six import iteritems, itervalues
#from six.moves import zip, range


import numpy as np

from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber

def remove_unused(bdf_filename):
    """
    removes unused:
     - nodes
     - properties
     - materials
     - coords
    """
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename, xref=False)

    #nids = model.nodes.keys()
    #cids =
    #nids = set(list(model.nodes.keys()))
    #cids = set(list(model.coords.keys()))
    #pids = set(list(model.properties.keys()))

    nids_used = set([])
    pids_used = set([])
    mids_used = set([])
    cids_used = set([])

    #card_types = list(model.card_count.keys())
    #card_map = model.get_card_ids_by_card_types(
        #card_types=card_types,
        #reset_type_to_slot_map=False,
        #stop_on_missing_card=True)

    #for nid, node in iteritems(model.nodes):
        #cids_used.update([node.Cp(), node.Cd()])

    skip_cards = [
        'ENDDATA', 'PARAM', 'EIGR', 'EIGRL', 'EIGB', 'EIGP', 'EIGC',
        'SPOINT', 'EPOINT', 'PAERO1', 'AEFACT', 'DESVAR', 'AESTAT',
        'AELIST', 'TRIM', 'SET1', 'FREQ', 'FREQ1', 'FREQ2',
        'TSTEP', 'TSTEPNL',
    ]
    load_types = [
        'GRAV', 'RANDPS', 'FORCE', 'FORCE1', 'FORCE2',
        'MOMENT', 'MOMENT1', 'MOMENT2',
    ]

    # could remove some if we look at the rid_trace
    #for cid, coord in iteritems(model.coords):
        #if coord.type in ['CORD1R', 'CORD1C', 'CORD1S']:
            #nids_used.update(node_ids)
        #elif coord.type in ['CORD1R', 'CORD1C', 'CORD1S']:
            #cids_used.update(coord.Rid())
        #else:
            #raise NotImplementedError(coord)

    for card_type, ids in iteritems(model._type_to_id_map):
    #for card_type, ids in iteritems(card_map):
        if card_type in ['CORD1R', 'CORD1C', 'CORD1S']:
            #print(ids)
            for cid in ids:
                coord = model.coords[cid]
                nids_used.update(coord.node_ids)
        elif card_type in ['CORD2R', 'CORD2C', 'CORD2S']:
            #print(ids)
            for cid in ids:
                coord = model.coords[cid]
                cids_used.add(coord.Rid())

        elif card_type in ['MAT1', 'MAT2', 'MAT8', 'MAT9']:
            # todo: MATS1, MATT1, etc.
            pass
        elif card_type in ['CTETRA', 'CPENTA', 'CPYRAM', 'CHEXA']:
            for eid in ids:
                elem = model.elements[eid]
                nids_used.update(elem.node_ids)
                pids_used.add(elem.Pid())
        elif card_type in ['PSOLID']:
            for pid in ids:
                prop = model.properties[pid]
                mids_used.add(prop.Mid())

        elif card_type in ['CONM2']:
            for eid in ids:
                elem = model.masses[eid]
                nids_used.add(elem.Nid())
                cids_used.add(elem.Cid())
                #print(elem.object_attributes())
                #print(elem.object_methods())
                #aaa

        elif card_type in ['CTRIA3', 'CQUAD4']:
            for eid in ids:
                elem = model.elements[eid]
                nids_used.update(elem.node_ids)
                pids_used.add(elem.Pid())
                if isinstance(elem.theta_mcid, int):
                    cids_used.add(elem.theta_mcid)
        elif card_type in ['CROD']:
            for eid in ids:
                elem = model.elements[eid]
                nids_used.update(elem.node_ids)
                pids_used.add(elem.Pid())
        elif card_type in ['PSHELL']:
            for pid in ids:
                prop = model.properties[pid]
                mids = [mid for mid in prop.material_ids if mid is not None]
                mids_used.update(mids)
        elif card_type in ['PCOMP', 'PCOMPG']:
            for pid in ids:
                prop = model.properties[pid]
                mids = prop.Mids()
                mids_used.update(mids)

        elif card_type in ['RBE2']:
            for eid in ids:
                elem = model.rigid_elements[eid]
                #print(elem.object_attributes())
                #print(elem.object_methods())
                nids_used.update(elem.independent_nodes)
                nids_used.update(elem.dependent_nodes)

        elif card_type in ['TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2']:
            pass
        elif card_type in load_types:
            for loads in itervalues(model.loads):
                for load in loads:
                    if load.type == 'FORCE':
                        nids_used.add(load.node_id)
                        cids_used.add(load.Cid())
                    elif load.type == 'GRAV':
                        cids_used.add(load.Cid())
                    elif load.type == 'RANDPS':
                        pass
                    else:
                        raise NotImplementedError(load)

        elif card_type in ['MPCADD', 'MPC']:
            for mpcs in itervalues(model.mpcs):
                for mpc in mpcs:
                    if mpc.type in ['MPCADD']:
                        pass
                    elif mpc.type in ['MPC']:
                        nids_used.update(mpc.node_ids)

        elif card_type in ['SPCADD', 'SPC1', 'SPC']:
            for spcs in itervalues(model.spcs):
                for spc in spcs:
                    if spc.type in ['SPCADD']:
                        pass
                    elif spc.type in ['SPC1', 'SPC']:
                        nids_used.update(spc.node_ids)

        elif card_type in ['TABLED1', 'TABLED2', 'TABLED3', 'TABLED4', 'TABDMP1', 'TABRND1']:
            pass
        elif card_type in ['SUPORT']:
            for suport in model.suport:
                nids_used.update(suport.node_ids)
        elif card_type in ['SUPORT1']:
            for suport1 in itervalues(model.suport1):
                nids_used.update(suport1.node_ids)
        elif card_type in ['GRID']:
            for nid, node in iteritems(model.nodes):
                cids_used.update([node.Cp(), node.Cd()])

        elif card_type in ['CBAR', 'CBEAM']:
            for eid in ids:
                elem = model.elements[eid]
                nids_used.update(elem.node_ids)
                pids_used.add(elem.Pid())
                if elem.g0 is not None:
                    assert isinstance(elem.g0, int), elem.g0
                    nids_used.add(elem.g0)
        elif card_type in ['PBAR', 'PBARL', 'PROD', 'PTUBE']:
            for pid in ids:
                prop = model.properties[pid]
                mids_used.add(prop.Mid())


        elif card_type in ['PBUSH']:
            pass
            #for pid in ids:
                #prop = model.properties[pid]
                #raise RuntimeError(prop)

        elif card_type in ['CBUSH']:
            for eid in ids:
                elem = model.elements[eid]
                nids_used.update(elem.node_ids)
                pids_used.add(elem.Pid())
                if elem.g0 is not None:
                    assert isinstance(elem.g0, int), elem.g0
                    nids_used.add(elem.g0)
                # TODO: cid

        elif card_type == 'AESURF':
            #CID1  | ALID1 | CID2   | ALID2
            for aesurf in itervalues(model.aesurf):
                cids_used.add(aesurf.Cid1())
                cid2 = aesurf.Cid2()
                if cid2 is not None:
                    cids_used.add(cid2)
        elif card_type in ['SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5']:
            pass
            #for spline_id in ids:
                #spline = model.splines[spline_id]
        elif card_type in ['CAERO1']:
            for eid in ids:
                caero = model.caeros[eid]
                # PID, LSPAN, LCHORD
                cids_used.add(caero.Cp())

        elif card_type in skip_cards:
            pass
        elif card_type in ['DCONSTR']:
            pass
        elif card_type == 'DRESP1':
            for dresp_id in ids:
                dresp = model.dresps[dresp_id]
                if dresp.property_type in ['PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM', 'PROD']:
                    pids_used.update(dresp.atti_values())
                elif dresp.property_type is None:
                    if dresp.response_type == 'WEIGHT':
                        pass
                    else:
                        raise NotImplementedError(dresp)
                else:
                    raise NotImplementedError(dresp)
        elif card_type == 'DVPREL1':
            for dvprel_id in ids:
                dvprel = model.dvprels[dvprel_id]
                if dvprel.Type in ['PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM', 'PROD']:
                    pids_used.add(dvprel.Pid())
                else:
                    raise NotImplementedError(dvprel)
        else:
            raise NotImplementedError(card_type)


    #for pid, prop in iteritems(model.properties):
        #prop = model.properties[pid]
        #if prop.type in no_materials:
            #continue
        #elif prop.type == 'PSHELL':
            #mids_used.extend([mid for mid in prop.material_ids if mid is not None])
        #elif prop.type == 'PCONEAX':
            #mids_used.extend([mid for mid in model.Mids() if mid is not None])

        #elif prop.type in prop_mid:
            #mids_used.append(prop.Mid())
        #elif prop.type in ['PCOMP', 'PCOMPG', 'PCOMPS']:
            #mids_used.extend(prop.Mids())

        #elif prop.type == 'PBCOMP':
            #mids_used.append(prop.Mid())
            #mids_used.extend(prop.Mids())
        #else:
            #raise NotImplementedError(prop)

    nids = set(model.nodes.keys())
    pids = set(model.properties.keys())
    cids = set(model.coords.keys())
    mids = set(model.materials.keys())
    nids_to_remove = list(nids - nids_used)
    pids_to_remove = list(pids - pids_used)
    mids_to_remove = list(mids - mids_used)
    cids_to_remove = list(cids - cids_used)
    for nid in nids_to_remove:
        del model.nodes[nid]
    model.log.debug('removed GRIDs %s' % nids_to_remove)

    for cid in cids_to_remove:
        del model.coords[cid]
    model.log.debug('removing coords %s' % cids_to_remove)

    for pid in pids_to_remove:
        del model.properties[pid]
    model.log.debug('removing properties %s' % pids_to_remove)

    for mid in mids_to_remove:
        del model.materials[mid]
    model.log.debug('removing materials %s' % mids_to_remove)



def remove_unassociated_nodes(bdf_filename, bdf_filename_out, renumber=False,
                              size=8, is_double=False):
    """
    Removes nodes from a model that are not referenced

    Parameters
    ----------
    bdf_filename : str
        the path to the bdf input file
    bdf_filename_out : str
        the path to the bdf output file
    renumber : bool
        should the model be renumbered
    size : int; {8, 16}; default=8
        the bdf write precision
    is_double : bool; default=False
        the field precision to write

    .. warning only considers elements
    .. renumber=False is not supported
    """
    model = BDF(debug=False)
    model.read_bdf(bdf_filename, xref=True)

    nids_used = set([])
    for element in itervalues(model.elements):
        nids_used.update(element.node_ids)
    #for element in itervalues(model.masses):
        #nids_used.update(element.node_ids)
    all_nids = set(model.nodes.keys())

    nodes_to_remove = all_nids - nids_used
    for nid in nodes_to_remove:
        del model.nodes[nid]

    if renumber:
        starting_id_dict = {
            'nid' : 1,
            'eid' : 1,
            'pid' : 1,
            'mid' : 1,
        }
        bdf_renumber(model, bdf_filename_out, size=size, is_double=is_double,
                     starting_id_dict=starting_id_dict)
    else:
        model.write_bdf(bdf_filename_out, size=size, is_double=is_double)

def remove_unassociated_properties(model, reset_type_to_slot_map=True):
    """remove_unassociated_properties"""
    pids_used = set()
    #elem_types = ['']
    card_types = list(model.card_count.keys())
    card_map = model.get_card_ids_by_card_types(card_types=card_types,
                                                reset_type_to_slot_map=reset_type_to_slot_map,
                                                stop_on_missing_card=True)
    skip_cards = [
        'GRID', 'PARAM',
        'MAT1', 'MAT2', 'MAT3', 'MAT4', 'MAT5', 'MAT8', 'MAT9', 'MAT10', 'MAT11',
        'CORD2R', 'CORD2C', 'CORD2S', 'CORD1R', 'CORD1C', 'CORD1S',

        'PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM', 'PROD', 'PELAS', 'PBUSH',
        'PBUSH1D', 'PBUSH2D', 'PSOLID', 'PRAC2D', 'PRAC3D',

        'CONM2',
        'RBE2', 'RBE3', 'RSPLINE',

        'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',
        'PAERO1', 'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5',
        'SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5',
        'AEROS', 'TRIM', 'DIVERG',
        'AERO', 'MKAERO1', 'MKAERO2', 'FLFACT', 'FLUTTER', 'GUST',
        'AELIST', 'AESURF', 'AESET1',
        'CONROD',
        'EIGRL', 'EIGB', 'EIGC', 'EIGR',
        'MPC', 'MPCADD', 'SPC1', 'SPCADD', 'SPCAX', 'SPCD',
        'PLOAD4',
        'DCONSTR', 'DESVAR',
        'ENDDATA',
    ]
    for card_type, ids in iteritems(card_map):
        if card_type in ['CTETRA', 'CPENTA', 'CPYRAM', 'CHEXA']:
            for eid in ids:
                elem = model.elements[eid]
                pids_used.add(elem.Pid())
        elif card_type in ['CTRIA3', 'CQUAD4', 'CBAR', 'CBEAM', 'CROD']:
            for eid in ids:
                elem = model.elements[eid]
                pids_used.add(elem.Pid())
        elif card_type in skip_cards:
            pass
        elif card_type == 'DRESP1':
            for dresp_id in ids:
                dresp = model.dresps[dresp_id]
                if dresp.property_type in ['PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM', 'PROD']:
                    pids_used.update(dresp.atti_values())
                elif dresp.property_type is None:
                    pass
                else:
                    raise NotImplementedError(dresp)
        elif card_type == 'DVPREL1':
            for dvprel_id in ids:
                dvprel = model.dvprels[dvprel_id]
                if dvprel.Type in ['PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM', 'PROD']:
                    pids_used.add(dvprel.Pid())
        else:
            raise NotImplementedError(card_type)
    all_pids = model.properties.keys()
    pids_to_remove = np.setdiff1d(all_pids, pids_used)
    for pid in pids_to_remove:
        del model.properties[pid]

def remove_unused_materials(model):
    """
    Removes all unused material cards

    .. warning:: doesn't support many cards
    """
    properties_without_materials = [
        'PELAS', 'PDAMP', 'PBUSH',
        'PELAST', 'PDAMPT', 'PBUSHT',
        'PGAP', 'PBUSH1D', 'PFAST', 'PVISC',
    ]
    prop_mid = [
        'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PSHEAR', 'PSOLID',
        'PROD', 'PRAC2D', 'PRAC3D', 'PLSOLID', 'PLPLANE', 'PPLANE',
        'PTUBE', 'PDAMP5',
    ]
    mids_used = []
    for elem in itervalues(model.elements):
        if elem.type in ['CONROD']:
            mids_used.append(elem.Mid())

    for pid, prop in iteritems(model.properties):
        prop = model.properties[pid]
        if prop.type in properties_without_materials:
            continue
        elif prop.type == 'PSHELL':
            mids_used.extend([mid for mid in prop.material_ids if mid is not None])
        elif prop.type == 'PCONEAX':
            mids_used.extend([mid for mid in model.Mids() if mid is not None])

        elif prop.type in prop_mid:
            mids_used.append(prop.Mid())
        elif prop.type in ['PCOMP', 'PCOMPG', 'PCOMPS']:
            mids_used.extend(prop.Mids())

        elif prop.type == 'PBCOMP':
            mids_used.append(prop.Mid())
            mids_used.extend(prop.Mids())
        else:
            raise NotImplementedError(prop)

    all_mids = set(model.materials.keys())
    for mid in all_mids:
        if mid not in mids_used:
            model.log.debug('removing mid=%s' % mid)
            del model.materials[mid]

    for dvmrel in itervalues(model.dvmrels):
        mids_used.append(dvmrel.Mid())

