from __future__ import annotations
from itertools import count
from typing import TYPE_CHECKING

from pyNastran.bdf.cards.utils import wipe_empty_fields
from pyNastran.bdf.bdf_interface.assign_type import interpret_value
from pyNastran.bdf.field_writer_8 import print_field_8, print_card_8
from pyNastran.bdf.field_writer_16 import print_field_16
from pyNastran.bdf.mesh_utils.forces_moments import get_temperatures_array
from pyNastran.bdf.mesh_utils.mpc_dependency import (
    get_mpc_node_ids, get_mpc_node_ids_c1, get_mpcs)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.cards.base_card import BaseCard

#def compare_params(fem1, fem2):
    #raise RuntimeError('is compare_parms used?')
    #compute(fem1.params, fem2.params)

def assert_fields(card1: BaseCard, card2: BaseCard) -> None:
    try:
        fields1 = wipe_empty_fields(card1.repr_fields())
        fields2 = wipe_empty_fields(card2.repr_fields())
    except:
        print("card1 = \n%s" % (card1))
        print("card2 = \n%s" % (card2))
        raise

    if len(fields1) != len(fields2):
        msg = ('len(fields1)=%s len(fields2)=%s\n%r\n%r\n%s\n%s'
               % (len(fields1), len(fields2), fields1, fields2,
                  print_card_8(fields1), print_card_8(fields2)))
        raise RuntimeError(msg)

    msg_end = ''
    max_int = 99999999
    for (i, field1, field2) in zip(count(), fields1, fields2):
        if isinstance(field1, int) and field1 > max_int:
            value1a = print_field_16(field1)
            value2a = print_field_16(field2)
        else:
            value1a = print_field_8(field1)
            value2a = print_field_8(field2)
        msg_end += '%-2s: %-8s %-8s\n' % (i, field1, field2)
        if value1a != value2a:
            if isinstance(field1, int) and field1 > max_int:
                value1 = print_field_16(interpret_value(value1a))
                value2 = print_field_16(interpret_value(value2a))
            else:
                value1 = print_field_8(interpret_value(value1a))
                value2 = print_field_8(interpret_value(value2a))

            if value1 != value2:
                msg = 'value1 != value2\n'
                msg += ('card_name=%s ID=%s i=%s field1=%r field2=%r value1=%r '
                        'value2=%r\n%r\n%r\n' % (fields1[0], fields1[1], i,
                                                 field1, field2, value1, value2,
                                                 fields1, fields2))
                raise RuntimeError(msg + msg_end)

def check_length(fem1, fem2, name):
    obj1 = getattr(fem1, name)
    obj2 = getattr(fem2, name)
    if len(obj1) != len(obj2):
        msg = f'len(fem1.{name})={len(obj1)} len(fem2.{name})={len(obj2)}'
        raise AssertionError(msg)


def compare_params(fem1: BDF, fem2: BDF) -> None:
    for key in fem1.params:
        card1 = fem1.params[key]
        card2 = fem2.params[key]
        assert_fields(card1, card2)

def compare_nodes(fem1: BDF, fem2: BDF) -> None:
    for key in fem1.nodes:
        card1 = fem1.nodes[key]
        card2 = fem2.nodes[key]
        assert_fields(card1, card2)

def compare_elements(fem1: BDF, fem2: BDF) -> None:
    for key in fem1.elements:
        card1 = fem1.elements[key]
        card2 = fem2.elements[key]
        assert_fields(card1, card2)

    for key in fem1.rigid_elements:
        card1 = fem1.rigid_elements[key]
        card2 = fem2.rigid_elements[key]
        assert_fields(card1, card2)

    for key in fem1.masses:
        card1 = fem1.masses[key]
        card2 = fem2.masses[key]
        assert_fields(card1, card2)

def compare_properties(fem1: BDF, fem2: BDF) -> None:
    for key in fem1.properties:
        card1 = fem1.properties[key]
        card2 = fem2.properties[key]
        assert_fields(card1, card2)

    for key in fem1.properties_mass:
        card1 = fem1.properties_mass[key]
        card2 = fem2.properties_mass[key]
        assert_fields(card1, card2)

def compare_materials(fem1: BDF, fem2: BDF) -> None:
    for key in fem1.materials:
        card1 = fem1.materials[key]
        card2 = fem2.materials[key]
        assert_fields(card1, card2)

    for key in fem1.creep_materials:
        card1 = fem1.creep_materials[key]
        card2 = fem2.creep_materials[key]
        assert_fields(card1, card2)

def compare_card_content(fem1, fem2):
    check_obj_names = [
        'params', 'nodes', 'spoints', 'epoints', 'points', 'gridb',
        #'elements', 'rigid_elements',
        'nsms', 'nsmadds',
        'properties', 'properties_mass', 'materials', 'creep_materials',
        'loads', 'coords',
        'spcs', 'spcadds', 'spcoffs', 'mpcs', 'mpcadds', 'dareas', 'dphases',
        'nlparms', 'tsteps', 'tstepnls', 'dmig', 'dmiax', 'dmij', 'dmik', 'dmiji', 'dequations',
        'sets', 'asets', 'bsets', 'csets', 'qsets', 'usets',
        'se_sets', 'se_bsets', 'se_csets', 'se_qsets', 'se_usets',
        'tables', 'tables_d', 'tables_m', 'random_tables', 'methods', 'cMethods']
    for name in check_obj_names:
        check_length(fem1, fem2, name)

    compare_params(fem1, fem2)
    compare_nodes(fem1, fem2)
    compare_elements(fem1, fem2)
    compare_properties(fem1, fem2)
    compare_materials(fem1, fem2)

    nid_map = fem2.nid_map
    for key in fem1.loads:
        loads1 = fem1.loads[key]
        loads2 = fem2.loads[key]
        get_temperatures_array(fem1, key, nid_map)
        for (card1, card2) in zip(loads1, loads2):
            assert_fields(card1, card2)

    for key in fem1.frequencies:
        freqs1 = fem1.frequencies[key]
        freqs2 = fem2.frequencies[key]
        assert isinstance(freqs1, list), freqs1
        assert isinstance(freqs2, list), freqs2
        for (card1, card2) in zip(freqs1, freqs2):
            assert_fields(card1, card2)

    for key in fem1.coords:
        card1 = fem1.coords[key]
        card2 = fem2.coords[key]
        assert_fields(card1, card2)

    for spc_id in fem1.spcadds:
        fem1.get_SPCx_node_ids(spc_id)
        fem1.get_SPCx_node_ids_c1(spc_id)
        fem1.get_reduced_spcs(spc_id)
        fem1.get_spcs(spc_id)

    for spc_id in fem1.spcs:
        fem1.get_SPCx_node_ids(spc_id)
        fem1.get_SPCx_node_ids_c1(spc_id)
        fem1.get_reduced_spcs(spc_id)
        fem1.get_spcs(spc_id)
        #card1 = fem1.spcs[key]
        #card2 = fem2.spcs[key]
        #assert_fields(card1, card2)

    for mpc_id in fem1.mpcadds:
        get_mpc_node_ids(fem1, mpc_id, consider_mpcadd=True)
        get_mpc_node_ids_c1(fem1, mpc_id, consider_mpcadd=True)
        fem1.get_reduced_mpcs(mpc_id, consider_mpcadd=True)
        get_mpcs(fem1, mpc_id)

    for mpc_id in fem1.mpcs:
        get_mpc_node_ids(fem1, mpc_id, consider_mpcadd=False)
        get_mpc_node_ids_c1(fem1, mpc_id, consider_mpcadd=False)
        fem1.get_reduced_mpcs(mpc_id, consider_mpcadd=False)
        get_mpcs(fem1, mpc_id)
        #card1 = fem1.mpcs[key]
        #card2 = fem2.mpcs[key]
        #assert_fields(card1, card2)

    for key in fem1.dareas:
        card1 = fem1.dareas[key]
        card2 = fem2.dareas[key]
        assert_fields(card1, card2)

    for key in fem1.tics:
        card1 = fem1.tics[key]
        card2 = fem2.tics[key]
        assert_fields(card1, card2)

    for key in fem1.dphases:
        card1 = fem1.dphases[key]
        card2 = fem2.dphases[key]
        assert_fields(card1, card2)

    for key in fem1.nlparms:
        card1 = fem1.nlparms[key]
        card2 = fem2.nlparms[key]
        assert_fields(card1, card2)

    for key in fem1.tsteps:
        card1 = fem1.tsteps[key]
        card2 = fem2.tsteps[key]
        assert_fields(card1, card2)

    for key in fem1.tstepnls:
        card1 = fem1.tstepnls[key]
        card2 = fem2.tstepnls[key]
        assert_fields(card1, card2)

    for key in fem1.dequations:
        card1 = fem1.dequations[key]
        card2 = fem2.dequations[key]
        msg = 'card1:\n%s\ncard2:\n%s' % (card1.write_card(), card2.write_card())
        assert card1.write_card() == card2.write_card(), msg
        #assert_fields(card1, card2)

    if fem1.dtable:
        card1 = fem1.dtable
        card2 = fem2.dtable
        assert_fields(card1, card2)

    for key in fem1.sets:
        card1 = fem1.sets[key]
        card2 = fem2.sets[key]

        if card1 != card2 and str(card1) != str(card2):
            # and card1.symmetric_difference(card2):
            # TODO: SET1s don't all handle comments...
            msg = 'SETx cards are not the same\n'
            msg += 'card1:\n%s\n' % str(card1)
            msg += 'card2:\n%s\n' % str(card2)
            msg += 'diff = %s' % str(card1.symmetric_difference(card2))
            raise AssertionError(msg)
        #assert_fields(card1, card2)

    dict_groups = [
        #'se_sets',
        #'usets', 'se_usets',
        'tables', 'random_tables',
        'methods', 'cMethods',
    ]
    #list_groups = [
        #'bsets', 'csets', 'qsets',
        #'se_bsets', 'se_csets', 'se_qsets',
    #]
    for name in dict_groups:
        group1 = getattr(fem1, name)
        group2 = getattr(fem2, name)
        for key in group1:
            try:
                card1 = group1[key]
                card2 = group2[key]
            except KeyError:
                msg = 'could not find key=%s for %s' % (key, name)
                raise KeyError(msg)
            assert_fields(card1, card2)

    for key in fem1.se_sets:
        card1 = fem1.se_sets[key]
        card2 = fem2.se_sets[key]
        assert_fields(card1, card2)

    #for key in fem1.tables:
        #card1 = fem1.tables[key]
        #card2 = fem2.tables[key]
        #assert_fields(card1, card2)

    #for key in fem1.random_tables:
        #card1 = fem1.random_tables[key]
        #card2 = fem2.random_tables[key]
        #assert_fields(card1, card2)

    #for key in fem1.methods:
        #card1 = fem1.methods[key]
        #card2 = fem2.methods[key]
        #assert_fields(card1, card2)

    #for key in fem1.cMethods:
        #card1 = fem1.cMethods[key]
        #card2 = fem2.cMethods[key]
        #assert_fields(card1, card2)

    compare_matrices(fem1, fem2)
    compare_optimization_content(fem1, fem2)
    compare_aero_content(fem1, fem2)
    compare_thermal_content(fem1, fem2)


def compare_matrices(fem1, fem2):
    """verifies that the DMIG, DMIJ, DMIJI, and DMIK matrices are the same"""
    for key in fem1.dmig:
        card1 = fem1.dmig[key]
        card2 = fem2.dmig[key]
        assert str(card1) == str(card2)
        #assert_fields(card1, card2)

    #for key in fem1.dmi:
        #card1 = fem1.dmi[key]
        #card2 = fem2.dmi[key]
        #assert str(card1) == str(card2)

    for key in fem1.dmij:
        card1 = fem1.dmij[key]
        card2 = fem2.dmij[key]
        assert str(card1) == str(card2)

    for key in fem1.dmiji:
        card1 = fem1.dmiji[key]
        card2 = fem2.dmiji[key]
        assert str(card1) == str(card2)

    for key in fem1.dmik:
        card1 = fem1.dmik[key]
        card2 = fem2.dmik[key]
        assert str(card1) == str(card2)

    for key in fem1.dmiax:
        card1 = fem1.dmiax[key]
        card2 = fem2.dmiax[key]
        assert str(card1) == str(card2)


def compare_thermal_content(fem1, fem2):
    """compares thermal cards"""
    assert len(fem1.bcs) == len(fem2.bcs)
    assert len(fem1.thermal_materials) == len(fem2.thermal_materials)
    assert len(fem1.phbdys) == len(fem2.phbdys)
    assert len(fem1.convection_properties) == len(fem2.convection_properties)

    for key in fem1.bcs:
        bcs1 = fem1.bcs[key]
        bcs2 = fem2.bcs[key]
        for (card1, card2) in zip(bcs1, bcs2):
            assert_fields(card1, card2)

    for key in fem1.thermal_materials:
        card1 = fem1.thermal_materials[key]
        card2 = fem2.thermal_materials[key]
        assert_fields(card1, card2)

    for key in fem1.phbdys:
        card1 = fem1.phbdys[key]
        card2 = fem2.phbdys[key]
        assert_fields(card1, card2)

    for key in fem1.convection_properties:
        card1 = fem1.convection_properties[key]
        card2 = fem2.convection_properties[key]


def compare_optimization_content(fem1, fem2):
    """compares optimization cards"""
    assert len(fem1.dconstrs) == len(fem2.dconstrs)
    assert len(fem1.desvars) == len(fem2.desvars)
    assert len(fem1.ddvals) == len(fem2.ddvals)
    assert len(fem1.dresps) == len(fem2.dresps)
    assert len(fem1.dvprels) == len(fem2.dvprels)
    assert len(fem1.dvmrels) == len(fem2.dvmrels)
    assert len(fem1.dvcrels) == len(fem2.dvcrels)

    for key in fem1.dconstrs:
        card1 = fem1.dconstrs[key]
        card2 = fem2.dconstrs[key]
        assert len(card1) == len(card2)
        #assert_fields(card1, card2)

    for key in fem1.desvars:
        card1 = fem1.desvars[key]
        card2 = fem2.desvars[key]
        assert_fields(card1, card2)

    for key in fem1.topvar:
        card1 = fem1.topvar[key]
        card2 = fem2.topvar[key]
        assert_fields(card1, card2)

    for key in fem1.ddvals:
        card1 = fem1.ddvals[key]
        card2 = fem2.ddvals[key]
        assert_fields(card1, card2)

    for key in fem1.dresps:
        card1 = fem1.dresps[key]
        card2 = fem2.dresps[key]
        assert_fields(card1, card2)

    for key in fem1.dvcrels:
        card1 = fem1.dvcrels[key]
        card2 = fem2.dvcrels[key]
        assert_fields(card1, card2)

    for key in fem1.dvmrels:
        card1 = fem1.dvmrels[key]
        card2 = fem2.dvmrels[key]
        assert_fields(card1, card2)

    for key in fem1.dvprels:
        card1 = fem1.dvprels[key]
        card2 = fem2.dvprels[key]
        assert_fields(card1, card2)
        card1.get_xinit_lower_upper_bound(fem1)


def compare_aero_content(fem1, fem2):
    """compares aero cards"""
    assert len(fem1.caeros) == len(fem2.caeros)
    assert len(fem1.paeros) == len(fem2.paeros)
    assert (fem1.aero is None) == (fem2.aero is None), 'fem1.aero_is_None=%s fem2.aero_is_None=%s' % (fem1.aero is None, fem2.aero is None)
    assert (fem1.aeros is None) == (fem2.aeros is None), 'fem1.aeros_is_None=%s fem2.aeros_is_None=%s' % (fem1.aeros is None, fem2.aeros is None)
    assert len(fem1.aeparams) == len(fem2.aeparams)
    assert len(fem1.aelinks) == len(fem2.aelinks)
    assert len(fem1.aelists) == len(fem2.aelists)
    assert len(fem1.aesurf) == len(fem2.aesurf)
    assert len(fem1.aesurfs) == len(fem2.aesurfs)
    assert len(fem1.aestats) == len(fem2.aestats)
    assert len(fem1.gusts) == len(fem2.gusts)
    assert len(fem1.flfacts) == len(fem2.flfacts)
    assert len(fem1.flutters) == len(fem2.flutters)
    assert len(fem1.mkaeros) == len(fem2.mkaeros)
    assert len(fem1.splines) == len(fem2.splines)
    assert len(fem1.trims) == len(fem2.trims)

    for key in fem1.caeros:
        card1 = fem1.caeros[key]
        card2 = fem2.caeros[key]
        assert_fields(card1, card2)

    for key in fem1.paeros:
        card1 = fem1.paeros[key]
        card2 = fem2.paeros[key]
        assert_fields(card1, card2)

    if fem1.aero is not None:
        card1 = fem1.aero
        card2 = fem2.aero
        assert_fields(card1, card2)

    if fem1.aeros is not None:
        card1 = fem1.aeros
        card2 = fem2.aeros
        assert_fields(card1, card2)

    for key in fem1.aeparams:
        card1 = fem1.aeparams[key]
        card2 = fem2.aeparams[key]
        assert_fields(card1, card2)

    for key in fem1.aelinks:
        aelinks1 = fem1.aelinks[key]
        aelinks2 = fem2.aelinks[key]
        for (card1, card2) in zip(aelinks1, aelinks2):
            assert_fields(card1, card2)

    for key in fem1.aelists:
        card1 = fem1.aelists[key]
        card2 = fem2.aelists[key]
        assert_fields(card1, card2)

    for key in fem1.aesurf:
        card1 = fem1.aesurf[key]
        card2 = fem2.aesurf[key]
        assert_fields(card1, card2)

    for key in fem1.aesurfs:
        card1 = fem1.aesurfs[key]
        card2 = fem2.aesurfs[key]
        assert_fields(card1, card2)

    for key in fem1.aestats:
        card1 = fem1.aestats[key]
        card2 = fem2.aestats[key]
        assert_fields(card1, card2)

    for key in fem1.gusts:
        card1 = fem1.gusts[key]
        card2 = fem2.gusts[key]
        assert_fields(card1, card2)

    for key in fem1.flfacts:
        card1 = fem1.flfacts[key]
        card2 = fem2.flfacts[key]
        assert_fields(card1, card2)

    for key in fem1.flutters:
        card1 = fem1.flutters[key]
        card2 = fem2.flutters[key]
        assert_fields(card1, card2)

    for i in range(len(fem1.mkaeros)):
        card1 = fem1.mkaeros[i]
        card2 = fem2.mkaeros[i]
        assert_fields(card1, card2)

    for key in fem1.splines:
        card1 = fem1.splines[key]
        card2 = fem2.splines[key]
        assert_fields(card1, card2)

    for key in fem1.trims:
        card1 = fem1.trims[key]
        card2 = fem2.trims[key]
        assert_fields(card1, card2)
