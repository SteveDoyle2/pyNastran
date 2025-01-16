from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from pyNastran.bdf.cards.properties.bars import _bar_areaL
if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def get_element_table(model: BDF) -> tuple[list[str], list, dict[str, int]]:
    card_types = []
    type_to_n = {}
    log = model.log
    for card in model.element_cards:
        if card.n == 0:
            continue
        card_types.append(card.type)
        type_to_n[card.type] = card.n

    group_type_to_card_type_map = {
        # True: can't collapse?
        # False: can?
        'Spring': (False, model.spring_element_cards),
        'Damper': (False, model.damper_element_cards + [model.cvisc]),
        'Bush': (False, [model.cbush, model.cbush1d, model.cgap]), #  model.cbush2d
        'Bar/Beam': (False, [model.cbar, model.cbeam, model.cbend]),
        'Shell': (False, model.shell_element_cards),
        'Solids': (False, model.solid_element_cards),
    }
    elements = []
    used_cards_set = set(['CONM2'])
    for group, (collapse_group, cards) in group_type_to_card_type_map.items():
        group_card_types = []
        ncards = sum([card.n for card in cards])
        if ncards == 0:
            continue

        ntotal = 0
        for card in cards:
            if card.n == 0 or card.type in used_cards_set:
                continue
            ntotal += card.n
            used_cards_set.add(card.type)
            group_row = (f'{card.type}: n={card.n}', True, card.n, [])
            group_card_types.append(group_row)
        group_tuple = (f'{group}: n={ntotal}', False, -1, group_card_types)

        if len(group_card_types) == 0:
            continue
        if len(group_card_types) == 1:
            #log.debug(f'group_row')
            #log.debug(f'group_row={str(group_row)}')
            elements.append(group_row)
        else:
            elements.append(group_tuple)

    other_cards = []
    other_card_names = []
    for card in model.element_cards:
        if card.n == 0 or card.type in used_cards_set:
            continue
        other_card_names.append(card.type)
        other_cards.append((f'{card.type}: {card.n}', True, card.n, []))

    if len(other_cards):
        elements.append(('Other', True, -1, other_cards))
    model.log.warning(f'other_card_names = {other_card_names}')
    return card_types, elements, type_to_n

def get_property_table(model: BDF) -> tuple[np.ndarray, list, dict[int, tuple[str, str]]]:
    properties = []
    pids = set([])
    pid_to_type_name: dict[int, tuple[str, str]] = {}
    log = model.log
    for card in model.property_cards:
        if card.n == 0:
            continue
        card_type = card.type
        for i, pid in enumerate(card.property_id):
            pids.add(pid)
            comment = card.comment.get(pid, '')
            if comment == '':
                if card_type == 'PBARL':
                    idim0, idim1 = card.idim[i]
                    beam_type =card.Type[i]
                    dim = card.dims[idim0:idim1]
                    area, i1, i2, i12 = _bar_areaL('PBARL', beam_type, dim, card)
                    if i1 is None:
                        i1 = 0.
                        i2 = 0.
                        i12 = 0.
                    comment = (f'Type={beam_type} '
                               f'area={engieering_format_str(area)} '
                               f'I1={engieering_format_str(i1)} '
                               f'I2={engieering_format_str(i2)} '
                               f'I12={engieering_format_str(i12)}')
                elif card_type =='PBEAML':
                    comment = f'Type={card.Type[i]}'
                elif card_type == 'PSHELL':
                    comment = f't={card.t[i]} mid={card.material_id[i, :]}'
                elif card.type == 'PCOMP':
                    ilayer0, ilayer1 = card.ilayer[i]
                    thickness = sum(card.thickness[ilayer0:ilayer1])
                    theta_str = ply_format(card.theta[ilayer0:ilayer1])
                    comment = f't={thickness:g} theta={theta_str}'
                else:
                    log.warning(f'card.type={card.type} is not supported for comments')
            pid_to_type_name[pid] = (card_type, comment)

    pids_array = np.array(list(pids), dtype='int32')
    pids_array.sort()
    pid_to_type_name = {key: value for key, value in sorted(pid_to_type_name.items())}
    for pid, (prop_type, name) in pid_to_type_name.items():
        properties.append((f'{pid}: {prop_type}: {name}', True, pid, []))
    return pids_array, properties, pid_to_type_name

def get_material_table(model: BDF) -> tuple[np.ndarray, list, dict[int, tuple[str, str]]]:
    materials = []
    mids = set([])
    mid_to_type_name: dict[int, tuple[str, str]] = {}
    log = model.log
    for card in model.material_cards:
        if card.n == 0:
            continue
        card_type = card.type
        for i, mid in enumerate(card.material_id):
            mids.add(mid)
            comment = card.comment.get(mid, '')
            if card_type == 'MAT1':
                comment = (f'E={engieering_format_str(card.E[i])} '
                           f'G={engieering_format_str(card.G[i])} '
                           f'nu={card.nu[i]}')
            elif card_type == 'MAT8':
                comment = (f'E11={engieering_format_str(card.E11[i])} '
                           f'E22={engieering_format_str(card.E22[i])} '
                           f'nu12={card.nu12[i]}')
            else:
                log.warning(f'card.type={card_type} is not supported for comments')
            mid_to_type_name[mid] = (card_type, comment)
            if mid == 2:
                assert card_type == 'MAT1', mid_to_type_name[mid]
    mids_array = np.array(list(mids), dtype='int32')
    mids_array.sort()
    mid_to_type_name = {key: value for key, value in sorted(mid_to_type_name.items())}
    for mid, (card_type, name) in mid_to_type_name.items():
        materials.append((f'{mid}: {card_type}: {name}', True, mid, []))
    return mids_array, materials, mid_to_type_name


def ply_format(theta: np.ndarray) -> str:
    thetai = theta.astype('int32')
    if np.allclose(theta, thetai):
        theta = thetai

    ntheta = len(theta)
    if ntheta > 4:
        is_even = (ntheta % 2 == 0)
        ntheta1 = ntheta // 2
        flag = ''
        if is_even:
            theta1 = theta[:ntheta1]
            theta2 = theta[ntheta1:]
            if np.allclose(theta1, theta2[::-1]):
                theta = theta1
                flag = ' Sym'
    theta_str = '/'.join(str(ti) for ti in theta)
    return theta_str + flag

def engieering_format_str(value: float) -> str:
    if value > 1.0:
        if value < 1000.0:
            return f'{value:g}'
        else:
            base, exp_str = f'{value:e}'.split('e')
            exp = int(exp_str)
            exp_remainder3 = exp % 3
            exp3 = exp - exp_remainder3
            base2 = float(base) * exp_remainder3
            return f'{base2}e+{exp3}'
    else:
        if value == 0:
            return '0'
        base, exp_str = f'{value:e}'.split('e')
        exp = int(exp_str)
        exp_remainder3 = -exp % 3
        exp3 = exp - exp_remainder3
        base2 = float(base) * exp_remainder3
        return f'{base2}e-{exp3}'
