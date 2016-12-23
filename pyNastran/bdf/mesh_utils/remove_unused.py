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

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber

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
    no_materials = [
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
        if prop.type in no_materials:
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
