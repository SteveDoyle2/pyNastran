from itertools import chain
import numpy as np

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.dmig import NastranMatrix
from pyNastran.bdf.bdf_interface.compare_card_content import compare_card_content
from pyNastran.bdf.mesh_utils.mass_properties import (
    mass_properties, mass_properties_nsm)  # mass_properties_breakdown


def compare(fem1: BDF,
            fem2: BDF,
            xref: bool=True,
            run_mass: bool=True,
            check: bool=True,
            print_stats: bool=True,
            quiet: bool=False) -> list[str]:
    """compares two fem objects"""
    fem1.log.debug('compare')
    diff_cards = compare_card_count(
        fem1, fem2,
        print_stats=print_stats, quiet=quiet)
    if xref and check:
        get_element_stats(fem1, fem2, run_mass=run_mass, quiet=quiet)
        get_matrix_stats(fem1, fem2)
    compare_card_content(fem1, fem2)
    #compare_params(fem1, fem2)
    #print_points(fem1, fem2)
    return diff_cards


#def compare_params(fem1, fem2):
    #raise RuntimeError('is compare_parms used?')
    #compute(fem1.params, fem2.params)

def compare_card_count(fem1: BDF,
                       fem2: BDF,
                       print_stats: bool=False,
                       quiet: bool=False) -> list[str]:
    """Checks that no cards from fem1 are lost when we write fem2"""
    cards1 = fem1.card_count
    cards2 = fem2.card_count
    for key in cards1:
        if key != key.upper():
            raise RuntimeError('Proper capitalization wasnt determined')
    if print_stats and not quiet:
        fem1.log.debug('get_bdf_stats')
        print(fem1.get_bdf_stats())
    else:
        fem1.get_bdf_stats()
    lines = compute_ints(cards1, cards2, fem1, fem2, quiet=quiet)
    return lines


def get_card_difference(fem1: BDF, fem2: BDF,
                        card_name: str, prefix: str='    ') -> str:
    msg = ''
    if not hasattr(fem1, '_slot_to_type_map'):
        return msg
    card_group = ''
    for card_groupi, card_names in fem1._slot_to_type_map.items():
        if card_name in card_names:
            card_group = card_groupi
            break
    assert card_group != '', card_name

    cards1 = getattr(fem1, card_group)
    cards2 = getattr(fem2, card_group)
    scalar_types = ['aero', 'aeros', 'doptprm']
    list_dict_types = ['loads', 'spcs', 'mpcs']
    list_types = ['suport']
    if card_group in scalar_types:
        assert cards1 is None or hasattr(cards1, 'type'), cards1.get_stats()
    elif card_group in list_types:
        # suport
        assert isinstance(cards1, list), card_group
    elif card_group in list_dict_types:
        # loads, spcs, mpcs
        assert isinstance(cards1, dict), card_group
        for sid, cards in cards1.items():
            assert isinstance(cards, list), (card_group, sid)
    else:
        # typical; nodes, elements, rigid_elements, ...
        assert isinstance(cards1, dict), card_group
        keys1 = set(list(cards1))
        keys2 = set(list(cards1))
        missing = list(keys1 - keys2)
        extra = list(keys2 - keys1)
        missing.sort()
        extra.sort()
        # missing_dict = {key: cards1[key] for key in missing}
        # extra_dict = {key: cards2[key] for key in extra}
        if len(missing):
            msg += f'{prefix}missing: {missing}\n'
        if len(missing):
            msg += f'{prefix}extra:   {extra}\n'
            # for key, value in missing_dict.items():
            #     lines.append(f'{prefix}missing:\n')
    return msg


def compute_ints(cards1: dict[str, int],
                 cards2: dict[str, int],
                 fem1: BDF,
                 fem2: BDF,
                 quiet: str=True) -> list[str]:
    """
    computes the difference / ratio / inverse-ratio between
    fem1 and fem2 to verify the number of card are the same:

    Examples
    --------

    name   fem1  fem2  diff  ratio  1/ratio
    ====   ====  ====  ==== ======  =======
    GRID      1     1     1     1.       1.
    *SPOINT  10     1     9    10.      0.1

    The * indicates a change, which may or may not be a problem.

    """
    fem1.log.debug('compute_ints')
    card_keys1 = set(cards1.keys())
    card_keys2 = set(cards2.keys())
    all_keys = card_keys1.union(card_keys2)
    diff_keys1 = list(all_keys.difference(card_keys1))
    diff_keys2 = list(all_keys.difference(card_keys2))
    #fem1.log.debug(f"all_keys = {all_keys}")

    list_keys1 = list(card_keys1)
    list_keys2 = list(card_keys2)
    #fem1.log.debug(f"list_keys1 = {list_keys1}")
    if diff_keys1 or diff_keys2:
        fem1.log.debug(' diff_keys1=%s diff_keys2=%s' % (diff_keys1, diff_keys2))

    for key in sorted(all_keys):
        msg = ''
        value1 = 0
        if key in list_keys1:
            value1 = cards1[key]

        value2 = 0
        if key in list_keys2:
            value2 = cards2[key]

        lost_keys = list(all_keys.difference(card_keys1))
        extra_keys = list(all_keys.difference(card_keys2))

        diff = abs(value1 - value2)
        star = ' '
        if diff and key not in ['INCLUDE']:
            star = '*'
        if key not in fem1.cards_to_read:
            star = '-'

        factor1 = divide(value1, value2)
        factor2 = divide(value2, value1)
        factor_msg = ''
        if not star or factor1 != factor2:
            if factor1 != factor2:
                factor_msg = f'diff={diff} factor1={factor1:g} {factor2:g}\n'
            msg += f'  {star}key={key:<7s} value1={value1:<7d} value2={value2:<7d}{factor_msg}'
            if star:
                msg += get_card_difference(fem1, fem2, key)

        if not quiet and msg:
            msg = msg.rstrip()
            print(msg)
    #return list_keys1 + list_keys2
    return diff_keys1 + diff_keys2


def get_element_stats(fem1: BDF,
                      unused_fem2: BDF,
                      run_mass: bool=True, quiet: bool=False) -> None:
    """verifies that the various element methods work"""
    fem1.log.debug('get_element_stats')
    for (unused_key, loads) in sorted(fem1.loads.items()):
        for load in loads:
            try:
                all_loads = load.get_loads()
                if not isinstance(all_loads, list):
                    raise TypeError('allLoads should return a list...%s'
                                    % (type(all_loads)))
            except Exception:
                raise
                # print("load statistics not available - load.type=%s "
                #       "load.sid=%s" % (load.type, load.sid))

    fem1._verify_bdf()

    if fem1.elements:
        fem1.get_elements_nodes_by_property_type()
    check_mass(fem1, run_mass=run_mass, quiet=quiet)


def check_mass(fem1: BDF, run_mass: bool=True, quiet: bool=False):
    if not run_mass:
        return
    mass1, cg1, inertia1 = mass_properties(fem1, reference_point=None, sym_axis=None)
    mass2, cg2, inertia2 = mass_properties_nsm(fem1, reference_point=None, sym_axis=None)
    #mass3, cg3, inertia3 = mass_properties_breakdown(fem1)[:3]
    if not quiet:
        if fem1.wtmass != 1.0:
            print('weight = %s' % (mass1 / fem1.wtmass))
        print(f'mass = {mass1}')
        print(f'cg   = {cg1}')
        print('Ixx=%s, Iyy=%s, Izz=%s \nIxy=%s, Ixz=%s, Iyz=%s' % tuple(inertia1))
    assert np.allclose(mass1, mass2), f'mass1={mass1} mass2={mass2}'
    assert np.allclose(cg1, cg2), f'mass={mass1}\ncg1={cg1} cg2={cg2}'
    assert np.allclose(inertia1, inertia2, atol=1e-5), f'mass={mass1} cg={cg1}\ninertia1={inertia1}\ninertia2={inertia2}\ndinertia={inertia1-inertia2}'

    for nsm_id in chain(fem1.nsms, fem1.nsmadds):
        mass, unused_cg, unused_inertia = mass_properties_nsm(
            fem1, reference_point=None, sym_axis=None, nsm_id=nsm_id)
        print('nsm_id=%s' % nsm_id)
        print('  mass = %s' % mass)
        print('  cg = %s' % cg1)
        print('  Ixx=%s, Iyy=%s, Izz=%s \n  Ixy=%s, Ixz=%s, Iyz=%s' % tuple(inertia1))

    reference_point = [10., 10., 10.]
    mass1, cg1, inertia1 = mass_properties(fem1, reference_point=reference_point, sym_axis=None)
    mass2, cg2, inertia2 = mass_properties_nsm(fem1, reference_point=reference_point, sym_axis=None)
    assert np.allclose(mass1, mass2), f'reference_point=[10., 10., 10.]; mass1={mass1} mass2={mass2}'
    assert np.allclose(cg1, cg2), f'reference_point=[10., 10., 10.]; mass={mass1} cg1={cg1} cg2={cg2}'
    assert np.allclose(inertia1, inertia2, atol=1e-5), f'reference_point=[10., 10., 10.]; mass={mass1} cg={cg1} inertia1={inertia1} inertia2={inertia2}'


def get_matrix_stats(fem1: BDF, unused_fem2: BDF) -> None:
    """Verifies the dmig.get_matrix() method works."""
    is_matrix = len(fem1.dmig) or len(fem1.dmi) or len(fem1.dmij) or len(fem1.dmik) or len(fem1.dmiji)
    if not is_matrix:
        return
    fem1.log.debug('get_matrix_stats')
    for (unused_key, dmig) in sorted(fem1.dmig.items()):
        try:
            if isinstance(dmig, NastranMatrix):
                dmig.get_matrix()
            else:
                print("statistics not available - "
                      "dmig.type=%s matrix.name=%s" % (dmig.type, dmig.name))
        except Exception:
            print("*stats - dmig.type=%s name=%s  matrix=\n%s"
                  % (dmig.type, dmig.name, str(dmig)))
            raise

    for (unused_key, dmi) in sorted(fem1.dmi.items()):
        try:
            if isinstance(dmi, NastranMatrix):
                dmi.get_matrix()
            else:
                print("statistics not available - "
                      "dmi.type=%s matrix.name=%s" % (dmi.type, dmi.name))
        except Exception:
            print("*stats - dmi.type=%s name=%s  matrix=\n%s"
                  % (dmi.type, dmi.name, str(dmi)))
            raise

    for (unused_key, dmij) in sorted(fem1.dmij.items()):
        try:
            if isinstance(dmij, NastranMatrix):
                dmij.get_matrix()
            else:
                print("statistics not available - "
                      "dmij.type=%s matrix.name=%s" % (dmij.type, dmij.name))
        except Exception:
            print("*stats - dmij.type=%s name=%s  matrix=\n%s"
                  % (dmij.type, dmij.name, str(dmij)))
            raise

    for (unused_key, dmiji) in sorted(fem1.dmiji.items()):
        try:
            if isinstance(dmiji, NastranMatrix):
                dmiji.get_matrix()
            else:
                print("statistics not available - "
                      "dmiji.type=%s matrix.name=%s" % (dmiji.type, dmiji.name))
        except Exception:
            print("*stats - dmiji.type=%s name=%s  matrix=\n%s"
                  % (dmiji.type, dmiji.name, str(dmiji)))
            raise

    for (unused_key, dmik) in sorted(fem1.dmik.items()):
        try:
            if isinstance(dmik, NastranMatrix):
                dmik.get_matrix()
            else:
                print("statistics not available - "
                      "dmik.type=%s matrix.name=%s" % (dmik.type, dmik.name))
        except Exception:
            print("*stats - dmik.type=%s name=%s  matrix=\n%s"
                  % (dmik.type, dmik.name, str(dmik)))
            raise


def divide(value1: int, value2: int) -> float:
    """
    Used to divide the number of cards to check that nothing was lost.
    Handles division by 0 by returning 0, which is the reciprocal.

    """
    if value1 == value2:  # good for 0/0
        return 1.0
    else:
        try:
            div_value = value1 / float(value2)
        except ZeroDivisionError:
            div_value = 0.
    return div_value
