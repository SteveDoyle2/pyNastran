from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import count

from pyNastran.bdf.cards.utils import wipe_empty_fields
from pyNastran.bdf.bdf_interface.assign_type import interpret_value
from pyNastran.bdf.field_writer import print_card
from pyNastran.bdf.field_writer_8 import print_field_8


def assert_fields(card1, card2):
    return
    try:
        fields1 = wipe_empty_fields(card1.repr_fields())
        fields2 = wipe_empty_fields(card2.repr_fields())
    except:
        print("card1 = \n%s" % card1)
        print("card2 = \n%s" % card2)
        raise

    if len(fields1) != len(fields2):
        msg = ('len(fields1)=%s len(fields2)=%s\n%r\n%r\n%s\n%s'
               % (len(fields1), len(fields2), fields1, fields2,
                  print_card(fields1), print_card(fields2)))
        raise RuntimeError(msg)

    for (i, field1, field2) in zip(count(), fields1, fields2):
        value1a = print_field_8(field1)
        value2a = print_field_8(field2)
        if value1a != value2a:
            value1 = print_field_8(interpret_value(value1a))
            value2 = print_field_8(interpret_value(value2a))

            if value1 != value2:
                msg = 'value1 != value2\n'
                msg += ('cardName=%s ID=%s i=%s field1=%r field2=%r value1=%r '
                        'value2=%r\n%r\n%r' % (fields1[0], fields1[1], i,
                                               field1, field2, value1, value2,
                                               fields1, fields2))
                raise RuntimeError(msg)

def check_length(fem1, fem2, name):
    obj1 = getattr(fem1, name)
    obj2 = getattr(fem2, name)
    if not len(obj2) == len(obj2):
        assert len(obj2) == len(obj2), 'len(fem1.%s)=%i len(fem2.%s)=%i' % (name, len(obj2), name, len(obj2))

def compare_card_content(fem1, fem2):
    check_obj_names = [
        'params', 'nodes', 'elements', #'rigid_elements',
        'properties', 'materials', 'creep_materials',
        #'loads',
        'coords',
        'spcs', 'spcoffs', 'mpcs', 'dareas',
        'nlparms', 'tsteps', 'tstepnls', 'dmigs', 'dequations', 'frequencies', 'sets',
        'tables', 'random_tables', 'methods', 'cMethods']
    for name in check_obj_names:
        check_length(fem1, fem2, name)

    for key in fem1.params:
        card1 = fem1.params[key]
        card2 = fem2.params[key]
        assert_fields(card1, card2)

    for key in fem1.nodes:
        card1 = fem1.nodes[key]
        card2 = fem2.nodes[key]
        assert_fields(card1, card2)

    for key in fem1.elements:
        card1 = fem1.elements[key]
        card2 = fem2.elements[key]
        assert_fields(card1, card2)

        assert len(fem1.rbe2) == len(fem2.rbe2), 'len(fem1.rbe2)=%i len(fem2.rbe2)=%i' % (len(fem1.rbe2), len(fem2.rbe2))
    for key in fem1.rbe2:
        card1 = fem1.rbe2[key]
        card2 = fem2.rbe2[key]
        assert_fields(card1, card2)

    for key in fem1.properties:
        card1 = fem1.properties[key]
        card2 = fem2.properties[key]
        assert_fields(card1, card2)

    for mat1, mat2 in zip(fem1.materials, fem2.materials):
        #print("key =", key)
        #card1 = fem1.materials[key]
        #card2 = fem2.materials[key]
        assert_fields(mat1, mat2, mat1.material_id)

    for key in fem1.creep_materials:
        card1 = fem1.creep_materials[key]
        card2 = fem2.creep_materials[key]
        assert_fields(card1, card2)

    #for key in fem1.loads:
    #    loads1 = fem1.loads[key]
    #    loads2 = fem2.loads[key]
    #    for (card1, card2) in zip(loads1, loads2):
    #        assert_fields(card1, card2)

    for key in fem1.coords:
        card1 = fem1.coords[key]
        card2 = fem2.coords[key]
        assert_fields(card1, card2)

    #for key in fem1.spcs:
        #card1 = fem1.spcs[key]
        #card2 = fem2.spcs[key]
        #assert_fields(card1, card2)

    #for key in fem1.spcadds:
        #card1 = fem1.spcadds[key]
        #card2 = fem2.spcadds[key]
        #assert_fields(card1, card2)

    #for key in fem1.mpcs:
        #card1 = fem1.mpcs[key]
        #card2 = fem2.mpcs[key]
        #assert_fields(card1, card2)

    #for key in fem1.mpcadds:
        #card1 = fem1.mpcadds[key]
        #card2 = fem2.mpcadds[key]
        #assert_fields(card1, card2)

    for key in fem1.dareas:
        card1 = fem1.dareas[key]
        card2 = fem2.dareas[key]
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
        card1 = fem1.dequations
        card2 = fem2.dequations
        assert_fields(card1, card2)

    for key in fem1.frequencies:
        card1 = fem1.frequencies[key]
        card2 = fem2.frequencies[key]
        assert_fields(card1, card2)

    for key in fem1.sets:
        card1 = fem1.sets[key]
        card2 = fem2.sets[key]
        assert_fields(card1, card2)

    for key in fem1.tables:
        card1 = fem1.tables[key]
        card2 = fem2.tables[key]
        assert_fields(card1, card2)

    for key in fem1.random_tables:
        card1 = fem1.random_tables[key]
        card2 = fem2.random_tables[key]
        assert_fields(card1, card2, key)

    for key in fem1.methods:
        card1 = fem1.methods[key]
        card2 = fem2.methods[key]
        assert_fields(card1, card2)

    for key in fem1.cMethods:
        card1 = fem1.cMethods[key]
        card2 = fem2.cMethods[key]
        assert_fields(card1, card2)

    compare_matrices(fem1, fem2)
    compare_optimization_content(fem1, fem2)
    compare_aero_content(fem1, fem2)
    compare_thermal_content(fem1, fem2)


def compare_matrices(fem1, fem2):
    for key in fem1.dmigs:
        card1 = fem1.dmigs[key]
        card2 = fem2.dmigs[key]
        assert str(card1) == str(card2)
        #assert_fields(card1, card2)

    #for key in fem1.dmis:
        #card1 = fem1.dmis[key]
        #card2 = fem2.dmis[key]
        #assert str(card1) == str(card2)

    for key in fem1.dmijs:
        card1 = fem1.dmijs[key]
        card2 = fem2.dmijs[key]
        assert str(card1) == str(card2)

    for key in fem1.dmijis:
        card1 = fem1.dmijis[key]
        card2 = fem2.dmijis[key]
        assert str(card1) == str(card2)

    for key in fem1.dmiks:
        card1 = fem1.dmiks[key]
        card2 = fem2.dmiks[key]
        assert str(card1) == str(card2)


def compare_thermal_content(fem1, fem2):
    assert len(fem1.bcs) == len(fem2.bcs)
    assert len(fem1.thermal_materials) == len(fem2.thermal_materials)
    assert len(fem1.phbdys) == len(fem2.phbdys)
    assert len(fem1.convection_properties) == len(fem2.convection_properties)

    for key in fem1.bcs:
        BCs1 = fem1.bcs[key]
        BCs2 = fem2.bcs[key]
        for (card1, card2) in zip(BCs1, BCs2):
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
    assert len(fem1.dconstrs) == len(fem2.dconstrs)
    assert len(fem1.desvars) == len(fem2.desvars)
    assert len(fem1.ddvals) == len(fem2.ddvals)
    assert len(fem1.dresps) == len(fem2.dresps)
    assert len(fem1.dvprels) == len(fem2.dvprels)

    for key in fem1.dconstrs:
        card1 = fem1.dconstrs[key]
        card2 = fem2.dconstrs[key]
        assert_fields(card1, card2)

    for key in fem1.desvars:
        card1 = fem1.desvars[key]
        card2 = fem2.desvars[key]
        assert_fields(card1, card2)

    for key in fem1.ddvals:
        card1 = fem1.ddvals[key]
        card2 = fem2.ddvals[key]
        assert_fields(card1, card2)

    for key in fem1.dresps:
        card1 = fem1.dresps[key]
        card2 = fem2.dresps[key]
        assert_fields(card1, card2)

    for key in fem1.dvprels:
        card1 = fem1.dvprels[key]
        card2 = fem2.dvprels[key]
        assert_fields(card1, card2)


def compare_aero_content(fem1, fem2):
    assert len(fem1.caeros) == len(fem2.caeros)
    assert len(fem1.paeros) == len(fem2.paeros)
    assert len(fem1.aero) == len(fem2.aero)
    assert len(fem1.aeros) == len(fem2.aeros)
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
    #assert len(fem1.trim) == len(fem2.trim)

    for key in fem1.caeros:
        card1 = fem1.caeros[key]
        card2 = fem2.caeros[key]
        assert_fields(card1, card2)

    for key in fem1.paeros:
        card1 = fem1.paeros[key]
        card2 = fem2.paeros[key]
        assert_fields(card1, card2)

    for key in fem1.aero:
        card1 = fem1.aero[key]
        card2 = fem2.aero[key]
        assert_fields(card1, card2)

    for key in fem1.aeros:
        card1 = fem1.aeros[key]
        card2 = fem2.aeros[key]
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

    #for key in fem1.trim:
        #card1 = fem1.trim[key]
        #card2 = fem2.trim[key]
        #assert_fields(card1, card2, key)
