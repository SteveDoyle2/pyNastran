#import os
from typing import cast
from itertools import count
import numpy as np

from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import read_bdf, BDF, SPLINE1, SET1
from pyNastran.bdf.mesh_utils.find_closest_nodes import find_closest_nodes

NEW_CARDS_OLD = [
    'GRID',
    'CORD2R', 'CORD2C', 'CORD2S',
    'CORD1R', 'CORD1C', 'CORD1S',
    'CAERO1', 'PAERO1', 'SET1', 'SET3',
    'CAERO2', 'PAERO2', 'SPLINE1', 'SPLINE2',
    'AERO', 'AEROS',
]
OLD_CARDS = [
    'GRID',
    'CORD2R', 'CORD2C', 'CORD2S',
    'CORD1R', 'CORD1C', 'CORD1S',
    'CAERO1', 'PAERO1', 'SET1', 'SET3',
    'CAERO2', 'PAERO2', 'SPLINE1', 'SPLINE2',
    'FLUTTER', 'MKAERO1', 'MKAERO2', 'AEFACT',
    'FLFACT',
    'AERO', 'AEROS',
]
NEW_CARDS = OLD_CARDS

def map_aero_model(model_old: PathLike | BDF,
                   model_new: PathLike | BDF,
                   bdf_filename_out: PathLike,
                   remove_new_aero_cards: bool=False,
                   include_mapped_grids: bool=False) -> None:
    """
    Parameters
    ----------
    model_old: BDF
        the baseline model
    model_new: BDF
        the model to map the SET points to
    remove_new_aero_cards: bool; default=False
        throw out the aero, aeros, caeros, splines, set1s
        and just use the model_old values as the baseline
    bdf_filename_out: str | Path
        output file

    """
    if not isinstance(model_old, BDF):
        model_old: BDF = read_bdf(model_old, read_cards=OLD_CARDS)

    if not isinstance(model_new, BDF):
        model_new: BDF = read_bdf(model_new, read_cards=NEW_CARDS)

    model_old.log.debug('-'*80)
    model_old.log.debug('Old Model:')
    model_old.log.debug(model_old.get_bdf_stats())
    model_old.log.debug('-'*80)
    model_new.log.debug('New Model:')
    model_new.log.debug(model_new.get_bdf_stats())
    model_old.log.debug('-'*80)

    spline_nids = get_spline_node_ids(model_old)
    xyz_cid0_old = get_xyz_cid0(model_old, spline_nids)

    if remove_new_aero_cards:
        model_new.log.debug('removing aero cards from new model...')
        model_new.caeros = {}
        model_new.splines = {}
        model_new.paeros = {}
        model_new.aero = model_old.aero
        model_new.aeros = model_old.aeros
        model_new.flfacts = model_old.flfacts
        model_new.flutters = model_old.flutters
        model_new.aefacts = model_old.aefacts
        model_new.trims = model_old.trims
        model_new.aesurf = model_old.aesurf
        model_new.aesurfs = model_old.aesurfs
        model_new.aestats = model_old.aestats

        add_methods = model_new._add_methods
        for paero_id, paero in model_old.paeros.items():
            paero.uncross_reference()
            add_methods.add_paero_object(paero)
        for caero_id, caero in model_old.caeros.items():
            caero.uncross_reference()
            add_methods.add_caero_object(caero)
        for spline_id, spline in model_old.splines.items():
            spline.uncross_reference()
            add_methods.add_spline_object(spline)
        for set_id, set1 in model_old.sets.items():
            set1.uncross_reference()
            add_methods.add_set_object(set1)
            #print(set1.rstrip())

    coord_ids_list = []
    for caero in model_new.caeros.values():
        coord_ids_list.append(caero.cp)

    nids_new = np.unique(list(model_new.nodes))
    xyz_cid0_new = get_xyz_cid0(model_new, nids_new)
    spline_nids_new = find_closest_nodes(
        xyz_cid0_new, nids_new,
        xyz_cid0_old)
    assert len(spline_nids_new) == len(spline_nids)

    inid_out = np.searchsorted(nids_new, spline_nids_new)
    xyz_cid0_new_out = xyz_cid0_new[inid_out, :]
    for inid, spline_id_old, spline_id_new in zip(count(), spline_nids, spline_nids_new):
        xyz_oldi = xyz_cid0_old[inid, :]
        xyz_newi = xyz_cid0_new_out[inid, :]
        print(f'{inid}: {spline_id_old} -> {spline_id_new}')
        print(f'    {xyz_oldi} -> {xyz_newi}')
   #aaa
    for spline_id, spline_old in model_old.splines.items():
        spline_new = model_new.splines[spline_id]
        set1_old_id = spline_old.setg
        set1_new_id = spline_new.setg
        assert set1_old_id == set1_new_id
        set1_old: SET1 = model_old.sets[set1_old_id]
        set1_new: SET1 = model_new.sets[set1_new_id]

        inids = np.searchsorted(spline_nids, set1_old.ids)
        nids_newi = spline_nids_new[inids]
        set1_new.ids = nids_newi.tolist()

    model_new.log.debug('Out Model:')
    model_new.log.debug(model_new.get_bdf_stats())
    coord_ids = np.unique(coord_ids_list)
    with open(bdf_filename_out, 'w') as bdf_file:
        model_new._write_aero(bdf_file)
        model_new._write_aero_control(bdf_file)
        model_new._write_static_aero(bdf_file)
        model_new._write_flutter(bdf_file)
        model_new._write_gust(bdf_file)
        model_new._write_sets(bdf_file)

        if include_mapped_grids:
            cid = 0
            mass = 1.0
            for nid in spline_nids_new:
                node = model_new.nodes[nid]
                bdf_file.write(str(node))
                eid = nid
                bdf_file.write(f'CONM2, {eid}, {nid}, {cid}, {mass}\n')
            for cid in coord_ids:
                coord = model_old.coords[cid]
                bdf_file.write(str(coord))

#find_closest_nodes(nodes_xyz: NDArray3float, nids: NDArrayNint,
#                       xyz_compare: NDArray3float,
#                       neq_max: int=1,
#                       tol: Optional[float]=None,
#                       msg: str='') -> NDArrayNint:

def get_spline_node_ids(model: BDF) -> np.ndarray:
    all_nids_list: list[int] = []
    for spline in model.splines.values():
        spline = cast(SPLINE1, spline)
        set1: SET1 = spline.setg_ref
        #print(set1.get_stats())
        nids = set1.ids
        all_nids_list.extend(nids)
    all_nids = np.unique(all_nids_list)
    return all_nids

def get_xyz_cid0(model: BDF, nids: list[int] | np.ndarray) -> np.ndarray:
    xyz_cid0 = np.zeros((len(nids), 3))
    for inid, nid in enumerate(nids):
        node = model.nodes[nid]
        xyz_cid0[inid, :] = node.get_position()
    return xyz_cid0
