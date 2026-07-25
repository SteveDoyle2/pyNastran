"""defines various shell element tests"""
import os
import unittest
from typing import Any

from cpylog import SimpleLogger

import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck
from pyNastran.dev.bdf_vectorized3.bdf_interface.breakdowns import NO_MASS
#from pyNastran.dev.bdf_vectorized3.bdf_interface.mass_properties import mass_properties_nsm
import pyNastran

PKG_PATH = pyNastran.__path__[0]
MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')


def mass_properties_nsm(model: BDF, nsm_id: int, debug: bool=False):
    nsms_dict = model.nsmadd.get_nsms_by_nsm_id()
    #nsmadd = model.nsmadd.slice_card_by_id(nsm_id)
    nsms = nsms_dict[nsm_id]
    shell_pids = []
    beam_pids = []
    bar_pids = []
    rod_pids = []
    conrod_eids = []

    elements_flag = []
    elements = []
    element_values = []
    #print('nsm_id', nsm_id)
    for nsm in nsms:
        #shell_pid_values = []
        _elements = []
        if nsm.type == 'NSM1':
            #nsm_id : array([1000])
            #nsm_type : array(['PSHELL'], dtype='<U6')
            #pid_eid : array([10])
            #value  : array([1.])
            #utypes = np.unique(nsm.nsm_type)
            for nsm_type, (ielement0, ielement1) in zip(nsm.nsm_type, nsm.ielement):
                pid_eid = nsm.pid_eid[ielement0:ielement1]
                if nsm_type in {'PSHELL', 'PCOMP'}:
                    cards = get_cards_by_property_id(pid_eid, model.shell_element_cards)
                    if len(cards) == 0:
                        continue
                    shell_pids.append((nsm_type, cards, 'per', pid_eid, nsm.value))
                elif nsm_type == 'PBARL':
                    cards = get_cards_by_property_id(pid_eid, [model.cbar])
                    if len(cards) == 0:
                        continue
                    bar_pids.append((nsm_type, cards, 'per', pid_eid, nsm.value))
                elif nsm_type == 'PBEAML':
                    cards = get_cards_by_property_id(pid_eid, [model.cbeam])
                    if len(cards) == 0:
                        continue
                    beam_pids.append((nsm_type, cards, 'per', pid_eid, nsm.value))
                elif nsm_type == 'PROD':
                    pid = pid_eid
                    if pid[0] == -1:
                        cards = [model.crod]
                    else:
                        cards = [card for card in [model.crod]
                                 if card.n > 0 and pid in card.property_id]
                    if len(cards) == 0:
                        continue
                    rod_pids.append((nsm_type, cards, 'per', pid_eid, nsm.value))
                elif nsm_type == 'CONROD':
                    cards = get_cards_by_element_id(pid_eid, [model.conrod])
                    if len(cards) == 0:
                        continue
                    conrod_eids.append((nsm_type, cards, 'per', pid_eid, nsm.value))
                elif nsm_type == 'ELEMENT':
                    _elements.append(pid_eid)
                else:
                    raise RuntimeError(nsm_type)
            if len(_elements):
                element = np.hstack(_elements, dtype=nsm.pid_eid.dtype)
                ones = np.ones(len(element), dtype=nsm.value.dtype)
                element_value = ones * nsm.value
                elements_flag.append('smear')
                elements.append(element)
                element_values.append(element_value)
                del _elements
        elif nsm.type == 'NSML':
            #ivalue : array([[0, 1]])
            #nsm_id : array([4000])
            #nsm_type : array(['PSHELL'], dtype='<U7')
            #nvalue : array([1])
            #pid_eid : array([10])
            #value  : array([1.])
            for nsm_type in nsm.nsm_type:
                if nsm_type in {'PSHELL', 'PCOMP'}:
                    cards = get_cards_by_property_id(nsm.pid_eid, model.shell_element_cards)
                    if len(cards) == 0:
                        continue
                    shell_pids.append((nsm_type, cards, 'smear', nsm.pid_eid, nsm.value))
                elif nsm_type in {'PBARL'}:
                    cards = get_cards_by_property_id(nsm.pid_eid, [model.cbar])
                    if len(cards) == 0:
                        continue
                    bar_pids.append((nsm_type, cards, 'smear', nsm.pid_eid, nsm.value))
                elif nsm_type in {'PBEAML'}:
                    cards = get_cards_by_property_id(nsm.pid_eid, [model.cbeam])
                    if len(cards) == 0:
                        continue
                    beam_pids.append((nsm_type, cards, 'smear', nsm.pid_eid, nsm.value))
                elif nsm_type in {'PROD'}:
                    cards = get_cards_by_property_id(nsm.pid_eid, [model.crod])
                    if len(cards) == 0:
                        continue
                    rod_pids.append((nsm_type, cards, 'smear', nsm.pid_eid, nsm.value))
                elif nsm_type in {'CONROD'}:
                    cards = get_cards_by_element_id(nsm.pid_eid, [model.conrod])
                    if len(cards) == 0:
                        continue
                    conrod_eids.append((nsm_type, cards, 'smear', nsm.pid_eid, nsm.value))
                elif nsm_type == 'ELEMENT':
                    _elements.append(nsm.pid_eid)
                else:
                    raise RuntimeError(nsm_type)
            if len(_elements):
                element = np.hstack(_elements, dtype=nsm.pid_eid.dtype)
                ones = np.ones(len(element), dtype=nsm.value.dtype)
                element_value = ones * nsm.value
                elements_flag.append('smear')
                elements.append(element)
                element_values.append(element_value)
                del _elements

        elif nsm.type == 'NSML1':
            distribution_type = 'smear'
            #ivalue : array([[0, 1]])
            #nsm_id : array([4000])
            #nsm_type : array(['PSHELL'], dtype='<U7')
            #pid_eid : array([10])
            #value  : array([1.])

            #nvalue : array([1])
            #assert len(nsm.nvalue) == 1
            #assert nsm.nvalue.max() == 1
            for nsm_type in nsm.nsm_type:
                pid_eid = nsm.pid_eid
                if nsm_type in {'PSHELL', 'PCOMP'}:
                    cards = get_cards_by_property_id(pid_eid, model.shell_element_cards)
                    if len(cards) == 0:
                        continue
                    shell_pids.append((nsm_type, cards, distribution_type, pid_eid, nsm.value))
                elif nsm_type == 'PBARL':
                    cards = get_cards_by_property_id(pid_eid, [model.cbar])
                    if len(cards) == 0:
                        continue
                    bar_pids.append((nsm_type, cards, distribution_type, pid_eid, nsm.value))
                elif nsm_type == 'PBEAML':
                    cards = get_cards_by_property_id(pid_eid, [model.cbeam])
                    if len(cards) == 0:
                        continue
                    beam_pids.append((nsm_type, cards, distribution_type, pid_eid, nsm.value))
                elif nsm_type == 'PROD':
                    cards = get_cards_by_property_id(pid_eid, [model.crod])
                    if len(cards) == 0:
                        continue
                    rod_pids.append((nsm_type, cards, distribution_type, pid_eid, nsm.value))
                elif nsm_type == 'CONROD':
                    cards = get_cards_by_element_id(pid_eid, [model.conrod])
                    if len(cards) == 0:
                        continue
                    conrod_eids.append((nsm_type, cards, distribution_type, pid_eid, nsm.value))

                elif nsm_type == 'ELEMENT':
                    _elements.append(pid_eid)
                else:
                    raise NotImplementedError(nsm_type)
            if len(_elements):
                element = np.hstack(_elements, dtype=nsm.pid_eid.dtype)
                ones = np.ones(len(element), dtype=nsm.value.dtype)
                element_value = ones * nsm.value
                elements_flag.append('per')
                elements.append(element)
                element_values.append(element_value)
                del _elements
        elif nsm.type == 'NSM':
            #nsm_id : array([3000])
            #nsm_type : array(['PSHELL'], dtype='<U7')
            #ivalue : array([[0, 1]])
            #nvalue : array([1])
            #pid_eid : array([10])
            #value  : array([1.])
            distribution_type = 'per'
            for nsm_type, (ivalue0, ivalue1) in zip(nsm.nsm_type, nsm.ivalue):
                pid_eid = nsm.pid_eid
                values = nsm.value[ivalue0 : ivalue1]
                # assert len(values) == 1, values
                if nsm_type == 'PSHELL':
                    cards = get_cards_by_property_id(pid_eid, model.shell_element_cards)
                    if len(cards) == 0:
                        continue
                    shell_pids.append((nsm_type, cards, distribution_type, pid_eid, nsm.value))
                elif nsm_type == 'PBARL':
                    cards = get_cards_by_property_id(pid_eid, [model.cbar])
                    if len(cards) == 0:
                        continue
                    bar_pids.append((nsm_type, cards, distribution_type, pid_eid, nsm.value))
                elif nsm_type == 'PBEAML':
                    cards = get_cards_by_property_id(pid_eid, [model.cbeam])
                    if len(cards) == 0:
                        continue
                    beam_pids.append((nsm_type, cards, distribution_type, pid_eid, nsm.value))
                elif nsm_type == 'PROD':
                    cards = get_cards_by_property_id(pid_eid, [model.crod])
                    if len(cards) == 0:
                        continue
                    rod_pids.append((nsm_type, cards, distribution_type, pid_eid, nsm.value))
                elif nsm_type == 'CONROD':
                    cards = get_cards_by_element_id(pid_eid, [model.conrod])
                    if len(cards) == 0:
                        continue
                    conrod_eids.append((nsm_type, cards, distribution_type, pid_eid, nsm.value))
                elif nsm_type == 'ELEMENT':
                    _elements.append(pid_eid)
                else:
                    raise NotImplementedError(nsm_type)
            if len(_elements):
                element = np.hstack(_elements, dtype=pid_eid.dtype)
                ones = np.ones(len(element), dtype=values.dtype)
                element_value = ones * nsm.value
                elements_flag.append('per')
                elements.append(element)
                element_values.append(element_value)
                del _elements
        else:
            print(nsm.get_stats())
            model.log.warning(f'skipping {nsm.type}')
            asdf

    # print('---------')
    # print(f'nsm_id = {nsm_id}')
    #shell_pids = np.unique(shell_pids)
    mass_list = []
    centroid_list = []

    for eid, value, flag in zip(elements, element_values, elements_flag):
        neid = len(eid)
        eid.sort()
        area_length_type_list = [] # np.full((neid, 5), np.nan, dtype='float64')
        cards = [card for card in model.element_cards if card.n > 0]

        #all_cards = []
        all_cards_eids = []
        for card in cards:
            if card.n == 0 or card.type in NO_MASS:
                continue
            #intersecti = np.intersect1d(card.element_id, eid)
            #if len(intersecti):
                #all_cards_eids.append(element_id)
                #all_cards.append(card)
        #ueids = np.unique(all_cards_eids)

        #for card in all_cards:
            #ieid_save = np.searchsorted(ueids, eid)
            if len(eid) == 1 and eid[0] == -1:
                eid2 = card.element_id
                # print(f'{card.type} base: eid={eid}; eid2={eid2}')
                card2 = card
            else:
                eid2 = [eidi for eidi in card.element_id
                        if eidi in eid]
                # print(f'{card.type} base: eid={eid}; eid2={eid2}')
                if len(eid2) == 0:
                    continue
                card2 = card.slice_card_by_element_id(eid2)
            #ieid = np.searchsorted(card.element_id, eid2)
            # print(f'   element_id={card2.element_id}')
            #print(f'  ieid.max={ieid.max()} neids={len(card.element_id)}')
            #ieid = ieid[ieid.max() < len(card.element_id)]
            #print(f'  ieid updated ieid={ieid}')
            #eids_lookup = (0 <= card.element_id[-1])
            #element_id = card.element_id # [card.element_id <= eid.max()]
            #eid2 = np.array([eidi for eidi in element_id])
            #ieid = ieid[element_id.max() > ieid.max()]
            #print(f'  element_id.max={element_id.max()} ieid.max()={ieid.max()}; '
                  #f'ieid={ieid}')
            #icorrect = (card.element_id[ieid] == eid)
            #ieid = ieid[icorrect]
            #print(f'  icorrect={icorrect}; ieid={ieid}')
            #if len(ieid) == 0:
                #continue

            if hasattr(card, 'length'):
                areai = card.length()
                card_type = 1
            elif hasattr(card, 'area'):
                areai = card.area()
                card_type = 2
            else:
                raise NotImplementedError(card)

            #old_card_types = area_length_type[ieid, 1]
            #if np.any(np.isfinite(old_card_types)):
                #raise RuntimeError(f'old_card_types={old_card_types}')

            centroidi = card.centroid()
            for areaii, centroidii in zip(areai, centroidi):
                area_length_type_list.append((areaii, card_type,
                                              centroidii[0], centroidii[1], centroidii[2]))
            #area_length_type[ieid, 0] = areai
            #area_length_type[ieid, 1] = card_type
            #area_length_type[ieid, 2:5] = centroidi
            del card2, areai, centroidi
        #assert np.isfinite(area_length_type[:, 1].min()), area_length_type[:, 1]
        area_length_type = np.array(area_length_type_list, dtype='float64')
        card_types = area_length_type[:, 1]
        if card_types.min() != card_types.max():
            raise RuntimeError(f'mixed type cards (e.g., shell/bar/beam in nsm={nsm_id}')
        if flag == 'smear':
            massi = area_length_type[:, 0] * value
        else:
            area_total = area_length_type[:, 0].sum()
            massi = area_length_type[:, 0] / area_total * value
        centroidi = area_length_type[:, 2:5]
        mass_list.append(massi)
        centroid_list.append(centroidi)

    _apply_bar_mass(rod_pids, mass_list, centroid_list)
    _apply_bar_mass(bar_pids, mass_list, centroid_list)
    _apply_bar_mass(beam_pids, mass_list, centroid_list)

    _apply_conrod_mass(conrod_eids, mass_list, centroid_list)

    for pid_type, cards, flag, pid, value in shell_pids:
        print(pid_type, cards, flag, pid, value)
        if flag == 'per':
            #total_area = 0.
            for card1 in cards:
                if pid[0] == -1:
                    card2 = card1
                else:
                    ipid = np.where(card1.property_id == pid)[0]
                    card2 = card1.slice_card_by_index(ipid)
                area = card2.area()
                #total_area += area.sum()
                centroid = card2.centroid()
                massi = area * value
                mass_list.append(massi)
                centroid_list.append(centroid)
                del area, card2, centroid, massi
            #print(f'total_area = {total_area}')
            x = 1
        else:
            assert flag == 'smear', flag
            total_area = 0.
            for card1 in cards:
                if pid[0] == -1:
                    card2 = card1
                else:
                    ipid = np.where(card1.property_id == pid)[0]
                    card2 = card1.slice_card_by_index(ipid)
                area = card2.area()
                total_area += area.sum()
            for card1 in cards:
                if pid[0] == -1:
                    card2 = card1
                else:
                    ipid = np.where(card1.property_id == pid)[0]
                    card2 = card1.slice_card_by_index(ipid)
                area = card2.area()
                centroid = card2.centroid()
                massi = area / total_area * value
                mass_list.append(massi)
                centroid_list.append(centroid)

    element_id, massi, centroidi, inertia = model.inertia()
    mass_list.append(massi)
    mass = np.hstack(mass_list)

    centroid_list.append(centroidi)
    centroid = np.vstack(centroid_list)
    #inertia = None
    mass_total = mass.sum()
    #assert np.allclose(mass_total, 8.), mass_total
    return mass_total, centroid, inertia

def get_cards_by_element_id(eid: np.ndarray, elements: list[Any]) -> list[Any]:
    if eid[0] == -1:
        cards = [card for card in elements if card.n > 0]
    else:
        cards = [card for card in elements
                 if card.n > 0 and eid in card.element_id]
    return cards

def get_cards_by_property_id(pid: np.ndarray, elements: list[Any]) -> list[Any]:
    if pid[0] == -1:
        cards = [card for card in elements if card.n > 0]
    else:
        cards = [card for card in elements
                 if card.n > 0 and pid in card.property_id]
    return cards

def _apply_bar_mass(bar_pids, mass_list, centroid_list):
    for pid_type, cards, flag, pid, value in bar_pids:
        if flag == 'per':
            #total_length = 0.
            for card1 in cards:
                if pid[0] == -1:
                    card2 = card1
                else:
                    assert len(pid) == 1
                    ipid = np.where(card1.property_id == pid)[0]
                    print(f'  card1.type={card1.type} property_id={card1.property_id} pid={pid}')
                    card2 = card1.slice_card_by_index(ipid)
                length = card2.length()
                #total_length += length.sum()
                centroid = card2.centroid()
                massi = length * value
                mass_list.append(massi)
                centroid_list.append(centroid)
                del length, card2, centroid, massi
            #print(f'total_length = {total_length}')
            x = 1
        else:
            assert flag == 'smear', flag
            total_length = 0.
            card_lengths = []
            for card1 in cards:
                if pid[0] == -1:
                    card2 = card1
                else:
                    assert len(pid) == 1
                    ipid = np.where(card1.property_id == pid)[0]
                    print(f'  card1.type={card1.type} property_id={card1.property_id} pid={pid}')
                    card2 = card1.slice_card_by_index(ipid)
                length = card2.length()
                total_length += length.sum()
                card_lengths.append((card2, length))

            for card2, length in card_lengths:
                centroid = card2.centroid()
                massi = length / total_length * value
                mass_list.append(massi)
                centroid_list.append(centroid)
                del length, card2, centroid, massi
            print(f'total_length = {total_length}')


def _apply_conrod_mass(conrod_eids, mass_list, centroid_list):
    for eid_type, cards, flag, eid, value in conrod_eids:
        if flag == 'per':
            #total_length = 0.
            for card1 in cards:
                if eid[0] == -1:
                    card2 = card1
                else:
                    assert len(eid) == 1
                    ieid = np.where(card1.element_id == eid)[0]
                    print(f'  card1.type={card1.type} element_id={card1.element_id} eid={eid}')
                    card2 = card1.slice_card_by_index(ieid)
                length = card2.length()
                #total_length += length.sum()
                centroid = card2.centroid()
                massi = length * value
                mass_list.append(massi)
                centroid_list.append(centroid)
                del length, card2, centroid, massi
            #print(f'total_length = {total_length}')
            x = 1
        else:
            total_length = 0.
            card_lengths = []
            for card1 in cards:
                if eid[0] == -1:
                    card2 = card1
                else:
                    assert len(eid) == 1
                    ieid = np.where(card1.element_id == eid)[0]
                    # print(f'  card1.type={card1.type} element_id={card1.element_id} eid={eid}')
                    card2 = card1.slice_card_by_index(ieid)
                length = card2.length()
                total_length += length.sum()
                card_lengths.append((card2, length))

            for card2, length in card_lengths:
                centroid = card2.centroid()
                massi = length * value
                mass_list.append(massi)
                centroid_list.append(centroid)
                del length, card2, centroid, massi
            # print(f'total_length = {total_length}')


class TestNsmV3(unittest.TestCase):
    def test_nsm_1002(self):
        eid_quad = 1
        eid_tri = 2
        pid_pshell = 10
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3
        nids = [1, 2, 3, 4]
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(eid_quad, pid_pshell, nids) # area=1.0
        model.add_ctria3(eid_tri, pid_pshell, nids[:-1]) # area=0.5
        model.add_pshell(pid_pshell, mid1=mid, t=0.1) #, nsm=None)
        model.add_mat1(mid, E, G, nu, rho=0.0)
        model.add_nsm1(1002, 'ELEMENT', 1.0, [eid_quad, eid_tri]) # correct; 1.5
        model.add_nsml1(2000, 'PSHELL', 1.0, pid_pshell, comment='nsml1') # correct; 1.0
        model.add_nsml1(2006, 'PSHELL', 1.0, 'ALL') # correct; 1.0

        model.setup()
        #mass1, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=1002, debug=False)
        #assert np.allclose(mass1.sum(), 1.5), mass1
        #mass1, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=2000, debug=False)
        #assert np.allclose(mass1.sum(), 1.0), mass1
        eids, mass1, unused_cg, unused_I = model.inertia(nsm_id=2006)
        assert np.allclose(mass1.sum(), 1.0), mass1

    def test_nsm_cquad4(self):
        eid_quad = 1
        eid_tri = 2
        eid_conrod = 3
        eid_crod = 4
        eid_pbeaml = 5
        eid_pbarl = 6
        pid_pbeaml = 40
        pid_pshell = 10
        pid_pbeaml = 21
        pid_pbarl = 31
        pid_prod = 41
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3
        nids = [1, 2, 3, 4]
        model = BDF(debug=True)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(eid_quad, pid_pshell, nids) # area=1.0
        model.add_ctria3(eid_tri, pid_pshell, nids[:-1]) # area=0.5
        model.add_conrod(eid_conrod, mid, [1, 2], A=1.0, j=0.0, c=0.0, nsm=0.0, comment='')

        x = [0., 0., 1.]
        g0 = None
        nids_beam = [1, 2]
        model.add_cbar(eid_pbarl, pid_pbarl, nids_beam, x, g0, offt='GGG', pa=0, pb=0,
                       wa=None, wb=None, comment='')
        model.add_cbeam(eid_pbeaml, pid_pbeaml, nids_beam, x, g0, offt='GGG', bit=None,
                        pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment='')
        model.add_crod(eid_crod, pid_prod, [1, 2])
        model.add_prod(pid_prod, mid, A=0.1)
        model.add_pshell(pid_pshell, mid1=mid, t=0.1) #, nsm=None)

        bar_type = 'BAR'
        dims = [1., 2.]
        xxb = [0.]
        model.add_pbarl(pid_pbarl, mid, bar_type, dims, group='MSCBML0', nsm=0., comment='')

        beam_type = 'BAR'
        dims = [[1., 2.]]
        nsm = [0.0]
        model.add_pbeaml(pid_pbeaml, mid, beam_type, xxb, dims, so=None, nsm=nsm,
                         group='MSCBML0', comment='')
        model.add_mat1(mid, E, G, nu, rho=0.0)

        # TODO: these are correct barring incorrect formulas
        model.add_nsm1(1000, 'PSHELL', 1.0, pid_pshell, comment='nsm1') # correct; 1.5; area=1.5 for PSHELL
        model.add_nsm1(1001, 'ELEMENT', 1.0, eid_quad) # correct; 1.0
        model.add_nsm1(1002, 'ELEMENT', 1.0, [eid_quad, eid_tri]) # correct; 1.5
        model.add_nsm1(1003, 'ELEMENT', 1.0, [eid_pbeaml]) # correct; 1.0
        model.add_nsm1(1004, 'ELEMENT', 1.0, eid_pbarl) # correct; 1.0
        model.add_nsm1(1005, 'ELEMENT', 1.0, 'ALL') # crash according to QRG b/c mixed type; 2.5
        model.add_nsm1(1006, 'PSHELL', 1.0, 'ALL') # correct; 1.5
        model.add_nsm1(1007, 'PSHELL', 1.0, [10, 'THRU', 12]) # correct; 1.5
        model.add_nsm1(1008, 'PSHELL', 1.0, [10, 'THRU', 12, 'BY', 2]) # correct; 1.5
        model.add_nsm1(1009, 'PBARL', 1.0, pid_pbarl) # correct; 1.0
        model.add_nsm1(1010, 'PBEAML', 1.0, pid_pbeaml) # correct; 1.0
        model.add_nsm1(1011, 'PROD', 1.0, pid_prod) # correct; 1.0
        model.add_nsm1(1012, 'CONROD', 1.0, eid_conrod) # correct; 1.0

        #model.add_nsml1(sid, nsm_type, value, ids)
        model.add_nsml1(2000, 'PSHELL', 1.0, pid_pshell, comment='nsml1') # correct; 1.0
        model.add_nsml1(2001, 'ELEMENT', 1.0, eid_quad) # correct; 1.0
        model.add_nsml1(2002, 'ELEMENT', 1.0, [eid_quad, eid_tri]) # correct; 1.0
        model.add_nsml1(2003, 'ELEMENT', 1.0, [eid_pbeaml]) # correct; 1.0
        model.add_nsml1(2004, 'ELEMENT', 1.0, eid_pbarl) # correct; 1.0
        model.add_nsml1(2005, 'ELEMENT', 1.0, 'ALL') # crash according to QRG b/c mixed type; 1.0
        model.add_nsml1(2006, 'PSHELL', 1.0, 'ALL') # correct; 1.0
        model.add_nsml1(2007, 'PSHELL', 1.0, [10, 'THRU', 12]) # correct; 1.0
        model.add_nsml1(2008, 'PSHELL', 1.0, [10, 'THRU', 12, 'BY', 2]) # correct; 1.0
        model.add_nsml1(2009, 'PBARL', 1.0, pid_pbarl) # correct; 1.0
        model.add_nsml1(2010, 'PBEAML', 1.0, pid_pbeaml) # correct; 1.0
        model.add_nsml1(2011, 'PROD', 1.0, pid_prod) # correct; 1.0
        model.add_nsml1(2012, 'CONROD', 1.0, eid_conrod) # correct; 1.0

        #model.add_nsml1(2011, 'PSHELL', 1.0, ['1240', 'THRU', '1250', None, None, # correct; 0.0
        #'2567', 'THRU', '2575',
        #'35689', 'THRU', '35700', None, None,
        #'76', 'THRU', '85',])
        #print(model.nsms[2011])

        model.add_nsm(3000, 'PSHELL', pid_pshell, 1.0, comment='nsm') # correct; 1.5
        model.add_nsm(3001, 'ELEMENT', eid_quad, 1.0) # correct; 1.0
        model.add_nsm(3003, 'ELEMENT', [eid_pbeaml], 1.0) # correct; 1.0
        model.add_nsm(3004, 'ELEMENT', eid_pbarl, 1.0) # correct; 1.0
        model.add_nsm(3009, 'PBARL', pid_pbarl, 1.0) # correct; 1.0
        model.add_nsm(3010, 'PBEAML', pid_pbeaml, 1.0) # correct; 1.0
        model.add_nsm(3011, 'PROD', pid_prod, 1.0) # correct; 1.0
        model.add_nsm(3012, 'CONROD', eid_conrod, 1.0) # correct; 1.0

        model.add_nsml(4000, 'PSHELL', pid_pshell, 1.0, comment='nsml') # correct; 1.0
        model.add_nsml(4001, 'ELEMENT', eid_quad, 1.0) # correct; 1.0
        model.add_nsml(4003, 'ELEMENT', [eid_pbeaml], 1.0) # correct; 1.0
        model.add_nsml(4004, 'ELEMENT', eid_pbarl, 1.0) # correct; 1.0
        model.add_nsml(4009, 'PBARL', pid_pbarl, 1.0) # correct; 1.0
        model.add_nsml(4010, 'PBEAML', pid_pbeaml, 1.0) # correct; 1.0
        model.add_nsml(4011, 'PROD', pid_prod, 1.0) # correct; 1.0
        model.add_nsml(4012, 'CONROD', eid_conrod, 1.0) # correct; 1.0

        model.pop_parse_errors()
        model.cross_reference()

        expected_dict = {
            # NSM1
            1000 : 1.5,
            1001 : 1.0,
            1002 : 1.5,
            1003 : 1.0,
            1004 : 1.0,
            1005 : -1.0,  # crash
            1006 : 1.5,
            1007 : 1.5,
            1008 : 1.5,
            1009 : 1.0,
            1010 : 1.0,
            1011 : 1.0,
            1012 : 1.0,

            #model.add_nsml1(sid, nsm_type, value, ids)
            # NSML1
            2000 : 1.0,
            2001 : 1.0,
            2002 : 1.0,
            2003 : 1.0,
            2004 : 1.0,
            2005 : -1.0, # crash
            2006 : 1.0,
            2007 : 1.0,
            2008 : 1.0,
            2009 : 1.0,
            2010 : 1.0,
            2011 : 1.0,
            2012 : 1.0,

            # NSM
            3000 : 1.5,
            3001 : 1.0,
            3003 : 1.0,
            3004 : 1.0,
            3009 : 1.0,
            3010 : 1.0,
            3011 : 1.0,
            3012 : 1.0,

            # NSM1
            4000 : 1.0,
            4001 : 1.0,
            4003 : 1.0,
            4004 : 1.0,
            4009 : 1.0,
            4010 : 1.0,
            4011 : 1.0,
            4012 : 1.0,
        }
        nsm_ids = np.hstack([
            nsm.nsm_id for nsm in model.nsm_cards
            if nsm.n > 0])
        nsm_ids.sort()
        for nsm_id in nsm_ids:
            mass1_expected = expected_dict[nsm_id]
            # assert mass1_expected >= 0, mass1_expected
            if mass1_expected == -1.0:
                with self.assertRaises(RuntimeError):
                    eids, mass1, unused_cg, unused_I = model.inertia(nsm_id=nsm_id)
            else:
                mass1, unused_cg, unused_I = model.inertia_sum(nsm_id=nsm_id)
                if mass1 != mass1_expected:
                    unused_mass2 = model.inertia(nsm_id=nsm_id)[1]
                    raise RuntimeError('nsm_id=%s mass != %s; mass1=%s' % (nsm_id, mass1_expected, mass1))
            #print('mass[%s] = %s' % (nsm_id, mass))
            #print('----------------------------------------------')

        model2 = save_load_deck(model, run_test_bdf=False)

    def test_nsm_prepare(self):
        """tests the NSMADD and all NSM cards using the prepare methods"""
        model = BDF()
        nsm_id = 100
        fields = ['NSM', nsm_id, 'ELEMENT',
                  1, 1.0,
                  2, 2.0,
                  3, 3.0,
                  4, 2.0]
        model.add_card(fields, 'NSM', comment='', is_list=True,
                       has_none=True)
        model.add_card(fields, 'NSML', comment='', is_list=True,
                       has_none=True)

        fields = ['NSM1', nsm_id, 'ELEMENT', 1.0, 1, 2, 3]
        model.add_card(fields, 'NSM1', comment='', is_list=True,
                       has_none=True)
        model.add_card(fields, 'NSML1', comment='', is_list=True,
                       has_none=True)

    def test_nsmadd(self):
        """tests the NSMADD and all NSM cards"""
        eid_quad = 1
        unused_eid_tri = 2
        unused_eid_conrod = 3
        unused_eid_crod = 4
        unused_eid_pbeaml = 5
        unused_eid_pbarl = 6
        unused_pid_pbeaml = 40
        pid_pshell = 10
        unused_pid_pbeaml = 21
        unused_pid_pbarl = 31
        unused_pid_prod = 41
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3
        nids = [1, 2, 3, 4]

        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(eid_quad, pid_pshell, nids) # area=1.0
        model.add_mat1(mid, E, G, nu, rho=0.0)
        model.add_pshell(pid_pshell, mid1=mid, t=0.1) #, nsm=None)

        model.add_nsm1(1000, 'PSHELL', 1.0, pid_pshell, comment='nsm1') # correct; 1.0
        model.add_nsml1(2000, 'PSHELL', 1.0, pid_pshell, comment='nsml1') # correct; 1.0
        model.add_nsml(3000, 'PSHELL', pid_pshell, 1.0, comment='nsml') # correct; 1.0
        model.add_nsml(4000, 'PSHELL', pid_pshell, 1.0, comment='nsml') # correct; 1.0
        model.add_nsmadd(5000, [1000, 2000, 3000, 4000], comment='nsmadd')
        model.add_nsmadd(5000, [1000, 2000, 3000, 4000], comment='nsmadd')
        model.cross_reference()


        #mass, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=1000)
        #self.assertAlmostEqual(mass.sum(), 1.0)
        #mass, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=2000)
        #self.assertAlmostEqual(mass.sum(), 1.0)
        #mass, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=3000)
        #self.assertAlmostEqual(mass.sum(), 1.0)
        #mass, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=4000)
        #self.assertAlmostEqual(mass.sum(), 1.0)

        eids, mass, unused_cg, unused_I = model.inertia(nsm_id=5000)
        self.assertAlmostEqual(mass.sum(), 8.0)
        model2 = save_load_deck(model)
        eids, mass, unused_cg, unused_I = model2.inertia(nsm_id=5000)

    #def test_nsm(self):
        #"""tests a complete nsm example"""
        #bdf_filename = os.path.join(MODEL_PATH, 'nsm', 'nsm.bdf')
        #bdf_filename = os.path.join(MODEL_PATH, 'nsm', 'TEST_NSM_SOL101.bdf')
        #model = read_bdf(bdf_filename)
        #print('    %6s %-9s %s' % ('nsm_id', 'mass', 'nsm'))
        #mass0 = mass_properties_nsm(model, debug=False)[0]
        #for nsm_id in sorted(chain(model.nsms, model.nsmadds)):
            #mass, cg, I = mass_properties_nsm(model, nsm_id=nsm_id, debug=False)
            #print('    %-6s %-9.4g %.4g' % (nsm_id, mass, mass-mass0))

        #area_breakdown = model.get_area_breakdown()
        #for pid in [20000, 20010]:
            #print('pid=%s area=%.3f' % (pid, area_breakdown[pid]))


    def test_nsmadd_short(self):
        """tests the NSMADD and all NSM cards"""
        eid_quad = 1
        pid_pshell = 10
        mid = 100
        E = 3.0e7
        G = None
        nu = 0.3
        nids = [1, 2, 3, 4]

        model = BDF(debug=True)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(eid_quad, pid_pshell, nids) # area=1.0
        model.add_mat1(mid, E, G, nu, rho=0.0)
        model.add_pshell(pid_pshell, mid1=mid, t=0.1) #, nsm=None)

        model.add_nsm1(1000, 'PSHELL', 1.0, pid_pshell, comment='nsm1') # correct; 1.0
        model.add_nsml1(2000, 'PSHELL', 1.0, pid_pshell, comment='nsml1') # correct; 1.0
        model.add_nsml(3000, 'PSHELL', pid_pshell, 1.0, comment='nsml') # correct; 1.0
        model.add_nsml(4000, 'PSHELL', pid_pshell, 1.0, comment='nsml') # correct; 1.0
        model.add_nsmadd(5000, [1000, 2000, 3000, 4000], comment='nsmadd')
        model.add_nsmadd(5000, [1000, 2000, 3000, 4000], comment='nsmadd')
        model.cross_reference()
        # model.pop_xref_errors()

        save_load_deck(model, run_mass_properties=False,
                       run_remove_unused=False, run_convert=False,
                       run_save_load_hdf5=False)

    def test_nsmadd_subset(self):
        model = BDF(debug=True)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4]) # area=1.0

        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 1., 0.])
        model.add_grid(14, [0., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14]) # area=1.0

        mid = 1
        E = 3.0e7
        G = None
        nu = 0.3
        pid_pshell = 1

        t = 0.1
        rho = 1.0
        area = 1.0
        nsm = 0.
        model.add_mat1(mid, E, G, nu, rho=rho)
        model.add_pshell(pid_pshell, mid1=mid, t=t) #, nsm=None)

        base_mass = 2 * area * (rho*t + nsm)
        nsml1_mass = 3.14

        model.add_nsml1(1, 'PSHELL', nsml1_mass, pid_pshell, comment='nsml1') # correct; 1.0
        model.cross_reference()

        expected_mass = base_mass + nsml1_mass
        eids, mass, cg, inertia = model.inertia(nsm_id=1)
        assert np.allclose(mass.sum(), expected_mass), f'mass={mass} expected_mass={expected_mass}'

        expected_mass = base_mass/2 + nsml1_mass/2
        eids, mass, cg, inertia = model.inertia(nsm_id=1, element_id=[1])
        assert np.allclose(mass.sum(), expected_mass), f'mass={mass} expected_mass={expected_mass}'

    def test_nsm1_subset_per_unit_area(self):
        """NSM1 applies nsm_value per unit area; subset should only get
        its own area's contribution.
        """
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1.0

        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [2., 0., 0.])
        model.add_grid(13, [2., 1., 0.])
        model.add_grid(14, [0., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=2.0

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        nsm_per_area = 1.0
        model.add_nsm1(10, 'PSHELL', nsm_per_area, 1)
        model.cross_reference()

        # Full model: both elements, areas 1+2, nsm_per_area=1.0
        # mass = 1*1.0 + 2*1.0 = 3.0
        eids, mass, cg, inertia = model.inertia(nsm_id=10)
        self.assertAlmostEqual(mass.sum(), 3.0)

        # Subset eid=1 only: area=1.0, mass = 1*1.0 = 1.0
        eids, mass, cg, inertia = model.inertia(nsm_id=10, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 1.0)

        # Subset eid=2 only: area=2.0, mass = 2*1.0 = 2.0
        eids, mass, cg, inertia = model.inertia(nsm_id=10, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 2.0)

    def test_nsml_subset_total_mass(self):
        """NSML distributes total mass proportional to area; subset should
        get proportional share based on area fraction.
        """
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1.0

        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [3., 0., 0.])
        model.add_grid(13, [3., 1., 0.])
        model.add_grid(14, [0., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=3.0

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        total_nsm_mass = 4.0
        model.add_nsml(20, 'PSHELL', 1, total_nsm_mass)
        model.cross_reference()

        # Full: total_nsm_mass distributed over area 1+3=4
        # eid=1 gets 4.0 * (1/4) = 1.0
        # eid=2 gets 4.0 * (3/4) = 3.0
        eids, mass, cg, inertia = model.inertia(nsm_id=20)
        self.assertAlmostEqual(mass.sum(), 4.0)

        # Subset eid=1: gets 1.0
        eids, mass, cg, inertia = model.inertia(nsm_id=20, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 1.0)

        # Subset eid=2: gets 3.0
        eids, mass, cg, inertia = model.inertia(nsm_id=20, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 3.0)

    def test_nsm_subset_element_type(self):
        """NSM with nsm_type=ELEMENT applies per-unit-area to specific elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1.0

        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [2., 0., 0.])
        model.add_grid(13, [2., 2., 0.])
        model.add_grid(14, [0., 2., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=4.0

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        # NSM applied to both elements, nsm_per_area=2.0
        model.add_nsm(30, 'ELEMENT', [1, 2], [2.0, 2.0])
        model.cross_reference()

        # Full: eid=1 mass=1*2=2, eid=2 mass=4*2=8, total=10
        eids, mass, cg, inertia = model.inertia(nsm_id=30)
        self.assertAlmostEqual(mass.sum(), 10.0)

        # Subset eid=1: mass=2
        eids, mass, cg, inertia = model.inertia(nsm_id=30, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 2.0)

        # Subset eid=2: mass=8
        eids, mass, cg, inertia = model.inertia(nsm_id=30, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 8.0)

    def test_nsml1_subset_line_elements(self):
        """NSML1 on PROD distributes total mass by length to line elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [4., 0., 0.])

        mid = 1
        pid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_prod(pid, mid, A=1.0)

        model.add_crod(1, pid, [1, 2])  # length=1.0
        model.add_crod(2, pid, [2, 3])  # length=3.0

        total_nsm = 8.0
        model.add_nsml1(40, 'PROD', total_nsm, pid)
        model.cross_reference()

        # Full: distribute 8.0 by length fraction: 1/(1+3)=0.25, 3/(1+3)=0.75
        # eid=1 gets 8*0.25=2, eid=2 gets 8*0.75=6
        eids, mass, cg, inertia = model.inertia(nsm_id=40)
        self.assertAlmostEqual(mass.sum(), 8.0)

        # Subset eid=1: gets 2.0
        eids, mass, cg, inertia = model.inertia(nsm_id=40, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 2.0)

        # Subset eid=2: gets 6.0
        eids, mass, cg, inertia = model.inertia(nsm_id=40, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 6.0)

    def test_nsm1_subset_line_elements(self):
        """NSM1 on PROD applies nsm_value per unit length to line elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [5., 0., 0.])

        mid = 1
        pid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_prod(pid, mid, A=1.0)

        model.add_crod(1, pid, [1, 2])  # length=2.0
        model.add_crod(2, pid, [2, 3])  # length=3.0

        nsm_per_length = 1.5
        model.add_nsm1(50, 'PROD', nsm_per_length, pid)
        model.cross_reference()

        # Full: eid=1 mass=2*1.5=3, eid=2 mass=3*1.5=4.5, total=7.5
        eids, mass, cg, inertia = model.inertia(nsm_id=50)
        self.assertAlmostEqual(mass.sum(), 7.5)

        # Subset eid=1: mass=3
        eids, mass, cg, inertia = model.inertia(nsm_id=50, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 3.0)

        # Subset eid=2: mass=4.5
        eids, mass, cg, inertia = model.inertia(nsm_id=50, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 4.5)

    def test_nsmadd_subset_combined(self):
        """NSMADD combines multiple NSM cards; subset filtering applies to all."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1.0

        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [1., 0., 0.])
        model.add_grid(13, [1., 1., 0.])
        model.add_grid(14, [0., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=1.0

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        # NSM1: 2.0 per unit area on PSHELL 1
        model.add_nsm1(100, 'PSHELL', 2.0, 1)
        # NSML1: 6.0 total distributed by area on PSHELL 1
        model.add_nsml1(200, 'PSHELL', 6.0, 1)
        # NSMADD combining both
        model.add_nsmadd(300, [100, 200])
        model.cross_reference()

        # NSM1: each element gets 2.0*1.0=2.0, total=4.0
        # NSML1: each element gets 6.0*(1/(1+1))=3.0, total=6.0
        # Combined: total=10.0
        eids, mass, cg, inertia = model.inertia(nsm_id=300)
        self.assertAlmostEqual(mass.sum(), 10.0)

        # Subset eid=1: NSM1 gives 2.0, NSML1 gives 3.0, total=5.0
        eids, mass, cg, inertia = model.inertia(nsm_id=300, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 5.0)

    def test_nsml1_subset_unequal_areas(self):
        """NSML1 with 3 elements of different areas validates proportional split."""
        model = BDF(debug=False)
        # eid=1: area=1.0
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])

        # eid=2: area=2.0
        model.add_grid(11, [0., 0., 0.])
        model.add_grid(12, [2., 0., 0.])
        model.add_grid(13, [2., 1., 0.])
        model.add_grid(14, [0., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])

        # eid=3: area=3.0 (triangle, base=6, height=1 -> area=3)
        model.add_grid(21, [0., 0., 0.])
        model.add_grid(22, [6., 0., 0.])
        model.add_grid(23, [0., 1., 0.])
        model.add_ctria3(3, 1, [21, 22, 23])

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        total_nsm = 12.0
        model.add_nsml1(60, 'PSHELL', total_nsm, 1)
        model.cross_reference()

        # total area = 1+2+3 = 6
        # eid=1: 12*(1/6) = 2.0
        # eid=2: 12*(2/6) = 4.0
        # eid=3: 12*(3/6) = 6.0
        eids, mass, cg, inertia = model.inertia(nsm_id=60)
        self.assertAlmostEqual(mass.sum(), 12.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=60, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 2.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=60, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 4.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=60, element_id=[3])
        self.assertAlmostEqual(mass.sum(), 6.0)

        # Two-element subset: eid=1+3 should get 2+6=8
        eids, mass, cg, inertia = model.inertia(nsm_id=60, element_id=[1, 3])
        self.assertAlmostEqual(mass.sum(), 8.0)

    def test_nsm_subset_with_structural_mass(self):
        """NSM subset correctly adds to structural mass subset."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1.0

        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [1., 2., 0.])
        model.add_grid(13, [1., 3., 0.])
        model.add_grid(14, [0., 3., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=1.0

        mid = 1
        rho = 2.0
        t = 0.5
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=rho)
        model.add_pshell(1, mid1=mid, t=t)

        nsm_per_area = 3.0
        model.add_nsm1(70, 'PSHELL', nsm_per_area, 1)
        model.cross_reference()

        # structural mass per element = area * rho * t = 1.0 * 2.0 * 0.5 = 1.0
        # nsm mass per element = area * nsm_per_area = 1.0 * 3.0 = 3.0
        # total per element = 4.0

        # Full model
        eids, mass, cg, inertia = model.inertia(nsm_id=70)
        self.assertAlmostEqual(mass.sum(), 8.0)

        # Subset eid=1
        eids, mass, cg, inertia = model.inertia(nsm_id=70, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 4.0)

        # Subset eid=2
        eids, mass, cg, inertia = model.inertia(nsm_id=70, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 4.0)

    def test_nsm_subset_cg(self):
        """Validate CG is correct when using element_ids subset."""
        model = BDF(debug=False)
        # eid=1 at x=0..1
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # centroid=(0.5, 0.5, 0), area=1

        # eid=2 at x=2..4
        model.add_grid(11, [2., 0., 0.])
        model.add_grid(12, [4., 0., 0.])
        model.add_grid(13, [4., 1., 0.])
        model.add_grid(14, [2., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # centroid=(3, 0.5, 0), area=2

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        nsm_per_area = 1.0
        model.add_nsm1(80, 'PSHELL', nsm_per_area, 1)
        model.cross_reference()

        # Subset eid=1: mass=1.0, centroid=(0.5, 0.5, 0)
        eids, mass, cg, inertia = model.inertia(nsm_id=80, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 1.0)
        assert np.allclose(cg, [0.5, 0.5, 0.0]), f'cg={cg}'

        # Subset eid=2: mass=2.0, centroid=(3.0, 0.5, 0)
        eids, mass, cg, inertia = model.inertia(nsm_id=80, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 2.0)
        assert np.allclose(cg, [3.0, 0.5, 0.0]), f'cg={cg}'

    def test_nsml1_subset_conrod(self):
        """NSML1 on CONROD elements distributes total mass by length."""
        model = BDF(debug=True)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [5., 0., 0.])

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)

        model.add_conrod(1, mid, [1, 2], A=1.0)  # length=2
        model.add_conrod(2, mid, [2, 3], A=1.0)  # length=3

        total_nsm = 10.0
        model.add_nsml1(90, 'CONROD', total_nsm, [1, 2])
        model.cross_reference()

        # total length = 2+3 = 5
        # eid=1: 10*(2/5) = 4.0
        # eid=2: 10*(3/5) = 6.0
        eids, mass, cg, inertia = model.inertia(nsm_id=90)
        self.assertAlmostEqual(mass.sum(), 10.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=0, include_base_mass=False)
        self.assertAlmostEqual(mass.sum(), 0.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=90, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 4.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=90, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 6.0)


    def test_nsml1_subset_pcomp(self):
        """NSML1 on PCOMP distributes total mass by area to composite shells."""
        model = BDF(debug=False)
        # eid=1: area=1.0
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 10, [1, 2, 3, 4])

        # eid=2: area=4.0
        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [2., 2., 0.])
        model.add_grid(13, [2., 4., 0.])
        model.add_grid(14, [0., 4., 0.])
        model.add_cquad4(2, 10, [11, 12, 13, 14])

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pcomp(10, [mid], [0.1], [0.])

        total_nsm = 15.0
        model.add_nsml1(100, 'PCOMP', total_nsm, 10)
        model.cross_reference()

        # total area = 1+4 = 5
        # eid=1: 15*(1/5) = 3.0
        # eid=2: 15*(4/5) = 12.0
        eids, mass, cg, inertia = model.inertia(nsm_id=100)
        self.assertAlmostEqual(mass.sum(), 15.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=100, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 3.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=100, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 12.0)

    def test_nsm1_subset_all_keyword(self):
        """NSM1 with 'ALL' applies nsm to all properties of that type."""
        model = BDF(debug=False)
        # Two PSHELL properties, one element each
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1, pid=1

        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [3., 2., 0.])
        model.add_grid(13, [3., 3., 0.])
        model.add_grid(14, [0., 3., 0.])
        model.add_cquad4(2, 2, [11, 12, 13, 14])  # area=3, pid=2

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)
        model.add_pshell(2, mid1=mid, t=0.1)

        nsm_per_area = 2.0
        model.add_nsm1(110, 'PSHELL', nsm_per_area, 'ALL')
        model.cross_reference()

        # eid=1: 1.0*2.0 = 2.0
        # eid=2: 3.0*2.0 = 6.0
        eids, mass, cg, inertia = model.inertia(nsm_id=110)
        self.assertAlmostEqual(mass.sum(), 8.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=110, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 2.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=110, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 6.0)

    def test_nsml1_subset_all_keyword(self):
        """NSML1 with 'ALL' distributes total mass across all PSHELL elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1, pid=1

        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [3., 2., 0.])
        model.add_grid(13, [3., 3., 0.])
        model.add_grid(14, [0., 3., 0.])
        model.add_cquad4(2, 2, [11, 12, 13, 14])  # area=3, pid=2

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)
        model.add_pshell(2, mid1=mid, t=0.1)

        total_nsm = 20.0
        model.add_nsml1(120, 'PSHELL', total_nsm, 'ALL')
        model.cross_reference()

        # total area = 1+3 = 4
        # eid=1: 20*(1/4) = 5.0
        # eid=2: 20*(3/4) = 15.0
        eids, mass, cg, inertia = model.inertia(nsm_id=120)
        self.assertAlmostEqual(mass.sum(), 20.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=120, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 5.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=120, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 15.0)

    def test_nsm1_subset_cbar(self):
        """NSM1 on PBAR applies per-unit-length mass to CBAR elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [5., 0., 0.])

        mid = 1
        pid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pbar(pid, mid, area=1.0, i1=1.0, i2=1.0)

        model.add_cbar(1, pid, [1, 2], x=[0., 0., 1.], g0=None)  # length=2
        model.add_cbar(2, pid, [2, 3], x=[0., 0., 1.], g0=None)  # length=3

        nsm_per_length = 4.0
        model.add_nsm1(130, 'PBAR', nsm_per_length, pid)
        model.cross_reference()

        # eid=1: 2*4=8, eid=2: 3*4=12, total=20
        eids, mass, cg, inertia = model.inertia(nsm_id=130)
        self.assertAlmostEqual(mass.sum(), 20.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=130, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 8.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=130, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 12.0)

    def test_nsml1_subset_cbar(self):
        """NSML1 on PBAR distributes total mass by length to CBAR elements."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [5., 0., 0.])

        mid = 1
        pid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pbar(pid, mid, area=1.0, i1=1.0, i2=1.0)

        model.add_cbar(1, pid, [1, 2], x=[0., 0., 1.], g0=None)  # length=2
        model.add_cbar(2, pid, [2, 3], x=[0., 0., 1.], g0=None)  # length=3

        total_nsm = 10.0
        model.add_nsml1(140, 'PBAR', total_nsm, pid)
        model.cross_reference()

        # total length=5; eid=1: 10*(2/5)=4, eid=2: 10*(3/5)=6
        eids, mass, cg, inertia = model.inertia(nsm_id=140)
        self.assertAlmostEqual(mass.sum(), 10.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=140, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 4.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=140, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 6.0)

    def test_nsml_subset_element_type_line(self):
        """NSML/ELEMENT assigns total mass to each element directly."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [4., 0., 0.])

        mid = 1
        pid = 10
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_prod(pid, mid, A=1.0)

        model.add_crod(1, pid, [1, 2])  # length=1
        model.add_crod(2, pid, [2, 3])  # length=3

        # NSML/ELEMENT: each element gets its value as total mass
        model.add_nsml(150, 'ELEMENT', [1, 2], [5.0, 7.0])
        model.cross_reference()

        # eid=1 gets 5.0, eid=2 gets 7.0
        eids, mass, cg, inertia = model.inertia(nsm_id=150)
        self.assertAlmostEqual(mass.sum(), 12.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=150, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 5.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=150, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 7.0)

    def test_nsmadd_subset_area_and_line(self):
        """NSMADD with area (shell) and line (rod) NSM cards combined."""
        model = BDF(debug=False)
        # Shell element
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [2., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=2

        # Rod element
        model.add_grid(11, [0., 5., 0.])
        model.add_grid(12, [3., 5., 0.])

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)
        model.add_prod(2, mid, A=1.0)

        model.add_crod(2, 2, [11, 12])  # length=3

        # NSM1 on shell: 5.0 per area -> eid=1 gets 5*2=10
        model.add_nsm1(160, 'PSHELL', 5.0, 1)
        # NSM1 on rod: 2.0 per length -> eid=2 gets 2*3=6
        model.add_nsm1(170, 'PROD', 2.0, 2)
        # NSMADD
        model.add_nsmadd(180, [160, 170])
        model.cross_reference()

        eids, mass, cg, inertia = model.inertia(nsm_id=180)
        self.assertAlmostEqual(mass.sum(), 16.0)

        # Subset shell only
        eids, mass, cg, inertia = model.inertia(nsm_id=180, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 10.0)

        # Subset rod only
        eids, mass, cg, inertia = model.inertia(nsm_id=180, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 6.0)

    def test_nsm_subset_inertia(self):
        """Validate inertia tensor for single-element subset with known geometry."""
        model = BDF(debug=False)
        # Single unit square at origin (centroid at 0.5, 0.5, 0)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1, centroid=(0.5,0.5,0)

        # Second element far away
        model.add_grid(11, [10., 0., 0.])
        model.add_grid(12, [11., 0., 0.])
        model.add_grid(13, [11., 1., 0.])
        model.add_grid(14, [10., 1., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=1, centroid=(10.5,0.5,0)

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        nsm_per_area = 3.0
        model.add_nsm1(190, 'PSHELL', nsm_per_area, 1)
        model.cross_reference()

        # Point mass approximation: mass=3, at (0.5, 0.5, 0), ref=(0,0,0)
        # Ixx = m*(dy^2+dz^2) = 3*(0.25+0) = 0.75
        # Iyy = m*(dx^2+dz^2) = 3*(0.25+0) = 0.75
        # Izz = m*(dx^2+dy^2) = 3*(0.25+0.25) = 1.5
        # Ixy = m*dx*dy = 3*0.5*0.5 = 0.75
        mass, cg, inertia = model.inertia_sum(
            nsm_id=190, element_id=[1],
            reference_point=[0., 0., 0.], inertia_reference='ref')
        self.assertAlmostEqual(mass, 3.0)
        assert np.allclose(cg, [0.5, 0.5, 0.0]), f'cg={cg}'
        # inertia = [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
        self.assertAlmostEqual(inertia[0], 0.75, places=5)
        self.assertAlmostEqual(inertia[1], 0.75, places=5)
        self.assertAlmostEqual(inertia[2], 1.5, places=5)
        self.assertAlmostEqual(inertia[3], 0.75, places=5)
        self.assertAlmostEqual(inertia[4], 0.0, places=5)
        self.assertAlmostEqual(inertia[5], 0.0, places=5)

    def test_nsm1_subset_multiple_pids(self):
        """NSM1 applied to multiple PSHELL pids via THRU; subset filtering."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1, pid=1

        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [2., 2., 0.])
        model.add_grid(13, [2., 3., 0.])
        model.add_grid(14, [0., 3., 0.])
        model.add_cquad4(2, 2, [11, 12, 13, 14])  # area=2, pid=2

        model.add_grid(21, [0., 4., 0.])
        model.add_grid(22, [3., 4., 0.])
        model.add_grid(23, [3., 5., 0.])
        model.add_grid(24, [0., 5., 0.])
        model.add_cquad4(3, 3, [21, 22, 23, 24])  # area=3, pid=3

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)
        model.add_pshell(2, mid1=mid, t=0.1)
        model.add_pshell(3, mid1=mid, t=0.1)

        nsm_per_area = 1.0
        # NSM1 on pids 1 THRU 3
        model.add_nsm1(200, 'PSHELL', nsm_per_area, [1, 'THRU', 3])
        model.cross_reference()

        # eid=1: 1*1=1, eid=2: 2*1=2, eid=3: 3*1=3, total=6
        eids, mass, cg, inertia = model.inertia(nsm_id=200)
        self.assertAlmostEqual(mass.sum(), 6.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=200, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 1.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=200, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 2.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=200, element_id=[3])
        self.assertAlmostEqual(mass.sum(), 3.0)

        # Subset of two elements
        eids, mass, cg, inertia = model.inertia(nsm_id=200, element_id=[1, 3])
        self.assertAlmostEqual(mass.sum(), 4.0)

    def test_nsm_subset_empty_result(self):
        """Subset with element_ids that don't overlap NSM gives zero mass."""
        model = BDF(debug=False)
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [1., 0., 0.])
        model.add_grid(3, [1., 1., 0.])
        model.add_grid(4, [0., 1., 0.])
        model.add_cquad4(1, 1, [1, 2, 3, 4])  # area=1

        model.add_grid(11, [0., 2., 0.])
        model.add_grid(12, [1., 2., 0.])
        model.add_grid(13, [1., 3., 0.])
        model.add_grid(14, [0., 3., 0.])
        model.add_cquad4(2, 1, [11, 12, 13, 14])  # area=1

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        # NSM only on eid=1
        model.add_nsm1(210, 'PSHELL', 5.0, 1)
        model.cross_reference()

        # Full model: only eid=1 has NSM via pid=1
        eids, mass, cg, inertia = model.inertia(nsm_id=210)
        self.assertAlmostEqual(mass.sum(), 10.0)

        # Subset with both elements (both have pid=1)
        eids, mass, cg, inertia = model.inertia(nsm_id=210, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 5.0)

    def test_nsm1_subset_ctria3(self):
        """NSM1 on CTRIA3 elements via PSHELL with subset filtering."""
        model = BDF(debug=False)
        # eid=1: triangle base=2, height=1 -> area=1.0
        model.add_grid(1, [0., 0., 0.])
        model.add_grid(2, [2., 0., 0.])
        model.add_grid(3, [0., 1., 0.])
        model.add_ctria3(1, 1, [1, 2, 3])  # area=1.0

        # eid=2: triangle base=4, height=3 -> area=6.0
        model.add_grid(11, [0., 5., 0.])
        model.add_grid(12, [4., 5., 0.])
        model.add_grid(13, [0., 8., 0.])
        model.add_ctria3(2, 1, [11, 12, 13])  # area=6.0

        mid = 1
        model.add_mat1(mid, 3.0e7, None, 0.3, rho=0.0)
        model.add_pshell(1, mid1=mid, t=0.1)

        nsm_per_area = 2.0
        model.add_nsm1(220, 'PSHELL', nsm_per_area, 1)
        model.cross_reference()

        # eid=1: 1.0*2.0=2.0, eid=2: 6.0*2.0=12.0
        eids, mass, cg, inertia = model.inertia(nsm_id=220)
        self.assertAlmostEqual(mass.sum(), 14.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=220, element_id=[1])
        self.assertAlmostEqual(mass.sum(), 2.0)

        eids, mass, cg, inertia = model.inertia(nsm_id=220, element_id=[2])
        self.assertAlmostEqual(mass.sum(), 12.0)

if __name__ == '__main__':  # pragma: no cover
    unittest.main()
