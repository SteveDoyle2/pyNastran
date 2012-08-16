from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import izip, count

from pyNastran.bdf.bdfInterface.BDF_Card import wipeEmptyFields
from pyNastran.bdf.fieldWriter import print_field, printCard


def assert_fields(card1, card2):
    try:
        fields1 = wipeEmptyFields(card1.reprFields())
        fields2 = wipeEmptyFields(card2.reprFields())
    except:
        print("card1 = \n%s" % (card1))
        print("card2 = \n%s" % (card2))
        raise

    if len(fields1) != len(fields2):
        msg = ('len(fields1)=%s len(fields2)=%s\n%r\n%r\n%s\n%s'
               % (len(fields1), len(fields2), fields1, fields2,
                  printCard(fields1), printCard(fields2)))
        raise RuntimeError(msg)

    for (i, field1, field2) in izip(count(), fields1, fields2):
        value1 = print_field(field1)
        value2 = print_field(field2)
        if value1 != value2:
            msg = 'value1 != value2\n'
            msg += ('cardName=%s ID=%s i=%s field1=%r field2=%r value1=%r '
                    'value2=%r\n%r\n%r' % (fields1[0], fields1[1], i,
                                           field1, field2, value1, value2,
                                           fields1, fields2))
            raise RuntimeError(msg)


def compare_card_content(fem1, fem2):
    assert len(fem1.params) == len(fem2.params)
    assert len(fem1.nodes) == len(fem2.nodes)
    assert len(fem1.elements) == len(fem2.elements)
    assert len(fem1.rigidElements) == len(fem2.rigidElements)
    assert len(fem1.properties) == len(fem2.properties)
    assert len(fem1.materials) == len(fem2.materials)
    assert len(fem1.creepMaterials) == len(fem2.creepMaterials)
    assert len(fem1.loads) == len(fem2.loads)
    assert len(fem1.coords) == len(fem2.coords)
    assert len(fem1.spcs) == len(fem2.spcs)
    assert len(fem1.spcadds) == len(fem2.spcadds)
    assert len(fem1.mpcs) == len(fem2.mpcs)
    assert len(fem1.mpcadds) == len(fem2.mpcadds)
    assert len(fem1.dareas) == len(fem2.dareas)
    assert len(fem1.nlparms) == len(fem2.nlparms)
    assert len(fem1.tsteps) == len(fem2.tsteps)
    assert len(fem1.tstepnls) == len(fem2.tstepnls)
    assert len(fem1.dmigs) == len(fem2.dmigs)
    assert len(fem1.dequations) == len(fem2.dequations)
    assert len(fem1.frequencies) == len(fem2.frequencies)
    assert len(fem1.sets) == len(fem2.sets)
    assert len(fem1.setsSuper) == len(fem2.setsSuper)
    assert len(fem1.tables) == len(fem2.tables)
    assert len(fem1.randomTables) == len(fem2.randomTables)
    assert len(fem1.methods) == len(fem2.methods), 'len(fem1.methods)=%s len(fem2.methods)=%s' % (len(fem1.methods), len(fem2.methods))
    assert len(fem1.cMethods) == len(fem2.cMethods), 'len(fem1.cMethods)=%s len(fem2.cMethods)=%s' % (len(fem1.cMethods), len(fem2.cMethods))

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

    for key in fem1.rigidElements:
        card1 = fem1.rigidElements[key]
        card2 = fem2.rigidElements[key]
        assert_fields(card1, card2)

    for key in fem1.properties:
        card1 = fem1.properties[key]
        card2 = fem2.properties[key]
        assert_fields(card1, card2)

    for key in fem1.materials:
        card1 = fem1.materials[key]
        card2 = fem2.materials[key]
        assert_fields(card1, card2)

    for key in fem1.creepMaterials:
        card1 = fem1.creepMaterials[key]
        card2 = fem2.creepMaterials[key]
        assert_fields(card1, card2)

    for key in fem1.loads:
        loads1 = fem1.loads[key]
        loads2 = fem2.loads[key]
        for (card1, card2) in izip(loads1, loads2):
            assert_fields(card1, card2)

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

    for key in fem1.setsSuper:
        card1 = fem1.setsSuper[key]
        card2 = fem2.setsSuper[key]
        assert_fields(card1, card2)

    for key in fem1.tables:
        card1 = fem1.tables[key]
        card2 = fem2.tables[key]
        assert_fields(card1, card2)

    for key in fem1.randomTables:
        card1 = fem1.randomTables[key]
        card2 = fem2.randomTables[key]
        assert_fields(card1, card2)

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
    assert len(fem1.thermalMaterials) == len(fem2.thermalMaterials)
    assert len(fem1.phbdys) == len(fem2.phbdys)
    assert len(fem1.convectionProperties) == len(fem2.convectionProperties)

    for key in fem1.bcs:
        BCs1 = fem1.bcs[key]
        BCs2 = fem2.bcs[key]
        for (card1, card2) in izip(BCs1, BCs2):
            assert_fields(card1, card2)

    for key in fem1.thermalMaterials:
        card1 = fem1.thermalMaterials[key]
        card2 = fem2.thermalMaterials[key]
        assert_fields(card1, card2)

    for key in fem1.phbdys:
        card1 = fem1.phbdys[key]
        card2 = fem2.phbdys[key]
        assert_fields(card1, card2)

    for key in fem1.convectionProperties:
        card1 = fem1.convectionProperties[key]
        card2 = fem2.convectionProperties[key]


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
        for (card1, card2) in izip(aelinks1, aelinks2):
            assert_fields(card1, card2)

    for key in fem1.aelists:
        card1 = fem1.aelists[key]
        card2 = fem2.aelists[key]
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

    for i in xrange(len(fem1.mkaeros)):
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
