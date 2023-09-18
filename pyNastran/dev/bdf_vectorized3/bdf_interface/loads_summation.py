from __future__ import annotations
from collections import defaultdict
from typing import Any, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.cards.loads.static_loads import (
        #LOAD, GRAV,
        SLOAD,
        #PLOAD, PLOAD1, PLOAD2, PLOAD4,
        #FORCE, FORCE1, FORCE2,
        #MOMENT, MOMENT1, MOMENT2,
        #TEMP, TEMPD, DTEMP, DTEMP,
        #SPCD, DEFORM,
        #PLOADX1,
        #RFORCE, RFORCE1,
        Loads as StaticLoad,
        LOAD,
    )


def get_static_loads_by_subcase_id(model: BDF, subcase_ids: list[int]=None) -> dict[int, Any]:
    res = {}
    for subcase_id, subcase in model.subcases:
        loadset_id = 0
        if 'LOADSET' in subcase:
            # Excitation specified by TYPE is applied load
            # --------------------------------------------
            # References DAREA entries as well as static and thermal load set
            # entries specified by the LID and TID fields, respectively, in the
            # selected LSEQ entry corresponding to EXCITEID.

            # Excitation specified by TYPE is enforced motion
            # -----------------------------------------------
            # References SPCD entries specified by the LID field in the selected LSEQ entry
            # corresponding to EXCITEID.
            #
            # If such entries indicate null enforced motion, the program will
            # then assume that the excitation is enforced motion using large mass
            # and will reference static and thermal load set entries specified by
            # the LID and TID fields, respectively, in the selected LSEQ entry
            # corresponding to EXCITEID, just as in the case of applied load excitation.
            loadset_id, unused_options = subcase['LOADSET']
        else:
            # Excitation specified by TYPE is applied load
            # --------------------------------------------
            # There is no LOADSET request in Case Control
            # EXCITEID may reference DAREA, static and thermal load set entries
            #
            # Excitation specified by TYPE is enforced motion
            # -----------------------------------------------
            # - There is no LOADSET request in Case Control
            # EXCITEID will reference SPCD entries. If such entries indicate
            # null enforced motion, the program will then assume that the
            # excitation is enforced motion using large mass and will
            # reference DAREA and static and thermal load set entries just
            # as in the case of applied load excitation.
            pass
    return res

def get_reduced_static_load(model: BDF, load: Optional[StaticLoad]=None) -> dict[int, list[tuple[float, Any]]]:
    """
    if load is None, then use self.load (LOAD)
    otherwise, use self.loadset
    """
    log = model.log
    if load is None:
        load: LOAD = model.load
    #model.dload.get_reduced_loads(filter_zero_scale_factors=False)
    SKIP_LOADS = {
        'LOAD', # 'SLOAD',
        'TEMP', 'TEMPD',
        #'QBDY3',
    }
    ## need a method to get all the load ids...
    filter_load_card_ids = False
    reduced_loads = defaultdict(list)
    if load.n:
        reduced_loads = load.get_reduced_loads(
            remove_missing_loads=False,
            filter_zero_scale_factors=False,
            stop_on_failure=True)
        filter_load_ids = True

    #new_reduced_loads = defaultdict(list)
    for loadi in model.loads:
        log.debug(loadi)
        if loadi.type in SKIP_LOADS or loadi.n == 0:
            continue
        uloads = np.unique(loadi.load_id)
        for uload in uloads:
            if filter_load_card_ids and uload in reduced_loads:
                continue
            i = np.where(loadi.load_id == uload)[0]
            sliced_card = loadi.slice_card_by_index(i)
            #if sliced_card.type == 'PLOAD4':
                #sliced_card.nvector
            scale_value = (1., sliced_card)
            reduced_loads[uload].append(scale_value)

    #for load_id, loads in new_reduced_loads.items():
        #reduced_loads[load_id] = loads
    return dict(reduced_loads)

def sum_forces_moments_elements(model: BDF,
                       p0: np.ndarray,
                       loadcase_id: int,
                       eids: np.ndarray,
                       nids: np.ndarray,
                       cid: int=0,
                       include_grav: bool=False) -> tuple[dict[int, np.ndarray],
                                                          dict[int, np.ndarray],]:
    return sum_forces_moments(model, p0, loadcase_id, cid=cid, include_grav=include_grav)

def sum_forces_moments(model: BDF,
                       p0: np.ndarray,
                       loadcase_id: int=None,
                       cid: int=0,
                       include_grav: bool=False) -> tuple[dict[int, np.ndarray],
                                                          dict[int, np.ndarray],]:
    # loads that aren't really referenceable from the LOAD card (e.g., cant' be scaled?)
    SKIP_LOADS = {
        # not done...
        #'GRAV', # 'SLOAD',
        # not applicible
        'LOAD', 'LSEQ',
        'TEMP', 'TEMPD',
        'LSEQ',
        #'QBDY3',
    }
    skip_summation = {
        'SPCD', 'DEFORM', 'LSEQ',
        # not done
        #'GRAV',
        #'SLOAD',
        # not applicable
        'TEMP', 'TEMPD',
    }
    log = model.log
    #reduced_loads_lseq = self.get_reduced_static_load(self.lseq)
    reduced_loads = model.get_reduced_static_load(model.load)
    loads_by_load_id = {}
    loads_by_subcase_id = {}

    #is_grav_loads = False
    eids = np.array([], dtype='int32')
    mass = np.array([], dtype='float64')
    cg = np.zeros((0, 3), dtype='float64')
    grav_loads = {'GRAV', 'ACCEL', 'ACCEL1'}
    for load_id, reduced_loadsi in sorted(reduced_loads.items()):
        force_moment_global = np.zeros(6, dtype='float64')
        for scale, loads in reduced_loadsi:
            for load in loads:
                if load.type in grav_loads:
                    #is_grav_loads = True
                    #mass_eids, mass = model.mass()
                    #centroid_eids, centroid = model.centroid()
                    eids, mass, cg, inertia = model.inertia()
                    break

    for load_id, reduced_loadsi in sorted(reduced_loads.items()):
        force_moment_global = np.zeros(6, dtype='float64')
        for scale, loads in reduced_loadsi:
            for load in loads:
                #if load.type in SKIP_LOADS or load.n == 0:
                    #continue

                if load.type in skip_summation or load.n == 0:
                    continue

                if not include_grav and load.type in grav_loads:
                    log.warning(f'include_grav={include_grav!r} and load.type={load.type!r}')
                    continue

                if load.type == 'GRAV':
                    #if len(mass) == 1:
                        #ma = (mass * load.N[:]).sum(axis=0)
                    #else:
                    ma = (mass[:, None] * load.N[None, :]).sum(axis=0)
                    force_sumi = np.hstack([ma, np.zeros(ma.shape, dtype='float64')])
                elif load.type in grav_loads:
                    raise NotImplementedError(load.type)
                else:
                    try:
                        force_sumi = load.sum_forces_moments()
                    except:
                        print(load.write())
                        raise
                force_moment_global += scale * force_sumi.sum(axis=0)
        loads_by_load_id[load_id] = force_moment_global

    for subcase_id, subcase in sorted(model.subcases.items()):
        if 'LOAD' not in subcase:
            continue
        load_id, options = subcase['LOAD']
        force_moment_global = loads_by_load_id[load_id]
        log.info(f'subcase={subcase_id} F={force_moment_global[:3]} M={force_moment_global[3:]}')
        loads_by_subcase_id[subcase_id] = force_moment_global

    #load_by_load_id = self.load.get_loads_by_load_id()
    #for load_id, loads in sorted(load_by_load_id.items()):
        #asdf
    #mass = 0.
    #for card in self.elements:
        #if card.n == 0 or card.type in NO_MASS:
            #continue
        #massi = card.mass()
        #if np.any(np.isnan(massi)):
            #self.log.warning(f'{card.type} has nan mass; mass={massi}')
        #mass += massi.sum()
    return loads_by_load_id, loads_by_subcase_id
