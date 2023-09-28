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
    print('nsm_id', nsm_id)
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
                assert len(values) == 1, values
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

    print('---------')
    print(f'nsm_id = {nsm_id}')
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
                print(f'{card.type} base: eid={eid}; eid2={eid2}')
                card2 = card
            else:
                eid2 = [eidi for eidi in card.element_id
                        if eidi in eid]
                print(f'{card.type} base: eid={eid}; eid2={eid2}')
                if len(eid2) == 0:
                    continue
                card2 = card.slice_card_by_element_id(eid2)
            #ieid = np.searchsorted(card.element_id, eid2)
            print(f'   element_id={card2.element_id}')
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
                    print(f'  card1.type={card1.type} element_id={card1.element_id} eid={eid}')
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
            print(f'total_length = {total_length}')


class TestNsm(unittest.TestCase):
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
        #assert np.allclose(mass1, 1.5), mass1
        #mass1, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=2000, debug=False)
        #assert np.allclose(mass1, 1.0), mass1
        mass1, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=2006, debug=False)
        assert np.allclose(mass1, 1.0), mass1

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
        model = BDF(debug=False)
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
        model.pop_xref_errors()

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
            if mass1_expected == -1.0:
                with self.assertRaises(RuntimeError):
                    mass1, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=nsm_id, debug=False)
            else:
                mass1, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=nsm_id, debug=False)
                if mass1 != mass1_expected:
                    unused_mass2 = mass_properties_nsm(model, nsm_id=nsm_id, debug=True)[0]
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
        model.pop_xref_errors()


        #mass, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=1000)
        #self.assertAlmostEqual(mass, 1.0)
        #mass, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=2000)
        #self.assertAlmostEqual(mass, 1.0)
        #mass, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=3000)
        #self.assertAlmostEqual(mass, 1.0)
        #mass, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=4000)
        #self.assertAlmostEqual(mass, 1.0)

        mass, unused_cg, unused_I = mass_properties_nsm(model, nsm_id=5000)
        self.assertAlmostEqual(mass, 8.0)
        model2 = save_load_deck(model)
        mass, unused_cg, unused_I = mass_properties_nsm(model2, nsm_id=5000)

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


if __name__ == '__main__':  # pragma: no cover
    unittest.main()
