from __future__ import print_function
import sys
import traceback
from six import iteritems

def verify_bdf(model, xref):
    #for key, card in sorted(model.params.items()):
        #card._verify(xref)
    for unused_key, card in sorted(iteritems(model.nodes)):
        try:
            card._verify(xref)
        except:
            print(str(card))
            raise

    _verify_dict(model.coords, xref)
    for unused_key, card in sorted(iteritems(model.elements)):
        try:
            card._verify(xref)
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print(repr(traceback.format_exception(exc_type, exc_value,
                                                  exc_traceback)))
            print(str(card))
            raise

    for eid, cbarao in sorted(model.ao_element_flags.items()):
        try:
            assert model.elements[eid].type == 'CBAR', 'CBARAO error: eid=%s is not a CBAR' % eid
        except:
            print(str(cbarao))
            raise

        _verify_dict(model.properties, xref)
        _verify_dict(model.properties_mass, xref)
        _verify_dict(model.materials, xref)

    _verify_dict(model.dresps, xref)
    _verify_dict(model.dvcrels, xref)
    _verify_dict(model.dvmrels, xref)
    _verify_dict(model.dvprels, xref)
    _verify_dict(model.dresps, xref)
    _verify_dict_list(model.dvgrids, xref)

    _verify_dict(model.dresps, xref)
    _verify_dict(model.gusts, xref)
    _verify_dict(model.tics, xref)

def _verify_dict(dict_obj, xref):
    for unused_key, card in sorted(dict_obj.items()):
        try:
            card._verify(xref)
        except:
            print(str(card))
            raise

def _verify_dict_list(dict_list, xref):
    for unused_key, cards in sorted(dict_list.items()):
        for card in cards:
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise

def validate_bdf(model):
    def _print_card(card):
        try:
            return card.write_card(size=8)
        except RuntimeError:
            return ''

    def _validate_dict_list(objects_dict):
        """helper method for validate"""
        ifailed = 0
        nmax_failed = 0
        assert isinstance(objects_dict, dict), type(objects_dict)
        for unused_key, objects in sorted(objects_dict.items()):
            assert isinstance(objects, list), type(objects)
            for obj in objects:
                #print('obj.get_stats =', obj.get_stats())
                #print(obj.rstrip())
                try:
                    obj.validate()
                except(ValueError, AssertionError, RuntimeError, IndexError) as error:
                    #exc_type, exc_value, exc_traceback = sys.exc_info()
                    # format_tb(exc_traceback)  # works; ugly
                    # format_exc(e) # works; short
                    #traceback.format_stack()
                    #print('validate_dict_list')
                    #print('traceback.format_stack()[:-1] = \n', ''.join(traceback.format_stack()[:-1]))
                    #model.log.info('info2')
                    #model.log.error('error2')
                    #msg = (
                        #'\nTraceback (most recent call last):\n' +
                        #''.join(traceback.format_stack()) +
                        ##''.join(traceback.format_tb(exc_traceback)) + #'\n' +
                        #'%s: %s\n' % (exc_type.__name__, exc_value) +
                        #obj.get_stats() +
                        #'----------------------------------------------------------------\n')
                    #model.log.error(msg)
                    model.log.error(('\n' + obj.get_stats() + '\n' + _print_card(obj)).rstrip())
                    ifailed += 1
                    if ifailed > nmax_failed:
                        raise
            if ifailed:
                raise

    def _validate_dict(objects):
        # type : (dict) -> None
        """helper method for validate"""
        assert isinstance(objects, dict), type(objects)
        ifailed = 0
        nmax_failed = 0
        for unused_id, obj in sorted(objects.items()):
            try:
                obj.validate()
            except(ValueError, AssertionError, RuntimeError, IndexError) as error:
                #exc_type, exc_value, exc_traceback = sys.exc_info()
                # format_tb(exc_traceback)  # works; ugly
                # format_exc(e) # works; short
                #traceback.format_stack()
                #msg = (
                    #'\nTraceback (most recent call last):\n' +
                    #''.join(traceback.format_stack()[:-1])
                #)
                #write_error(msg)
                #write_error(''.join(traceback.format_tb(exc_traceback)) + '\n')
                #model.log.error('\n' + obj.rstrip())
                #write_error('----------------------------------------------------------------\n')
                model.log.error(('\n' + obj.get_stats() + '\n' + _print_card(obj)).rstrip())
                _print_card(obj)
                ifailed += 1
                if ifailed > nmax_failed:
                    raise
        if ifailed:
            raise

    def _validate_list(objects):
        # type : (List) -> None
        """helper method for validate"""
        ifailed = 0
        nmax_failed = 0
        assert isinstance(objects, list), type(objects)
        for obj in objects:
            try:
                obj.validate()
            except(ValueError, AssertionError, RuntimeError) as error:
                #exc_type, exc_value, exc_traceback = sys.exc_info()
                # format_tb(exc_traceback)  # works; ugly
                # format_exc(e) # works; short
                #traceback.format_stack()
                #model.log.error(
                    #'\nTraceback (most recent call last):\n' +
                    #''.join(traceback.format_stack()[:-1]) +
                    #''.join(traceback.format_tb(exc_traceback)) + '\n' +
                    #'\n' + obj.rstrip() +
                    #'----------------------------------------------------------------\n')
                model.log.error(('\n' + obj.get_stats() + '\n' + _print_card(obj)).rstrip())
                ifailed += 1
                if ifailed > nmax_failed:
                    raise
        if ifailed:
            raise

    _validate_dict(model.nodes)
    _validate_dict(model.points)
    _validate_dict(model.coords)
    _validate_dict(model.elements)
    _validate_dict(model.properties)
    _validate_dict(model.rigid_elements)
    _validate_dict(model.plotels)
    _validate_dict(model.masses)
    _validate_dict(model.properties_mass)
    #------------------------------------------------
    _validate_dict(model.materials)
    _validate_dict(model.thermal_materials)
    _validate_dict(model.MATS1)
    _validate_dict(model.MATS3)
    _validate_dict(model.MATS8)
    _validate_dict(model.MATT1)
    _validate_dict(model.MATT2)
    _validate_dict(model.MATT3)
    _validate_dict(model.MATT4)
    _validate_dict(model.MATT5)
    _validate_dict(model.MATT8)
    _validate_dict(model.MATT9)
    _validate_dict(model.creep_materials)
    _validate_dict(model.hyperelastic_materials)

    #------------------------------------------------
    _validate_dict_list(model.load_combinations)
    _validate_dict_list(model.loads)
    _validate_dict_list(model.dloads)
    _validate_dict_list(model.dloads)

    #------------------------------------------------
    _validate_dict(model.nlpcis)
    _validate_dict(model.nlparms)
    _validate_dict(model.rotors)
    _validate_dict(model.tsteps)
    _validate_dict(model.tstepnls)
    _validate_dict_list(model.transfer_functions)
    _validate_dict(model.delays)

    #------------------------------------------------
    if model.aeros is not None:
        model.aeros.validate()
    _validate_dict(model.caeros)
    _validate_dict(model.paeros)
    _validate_dict(model.splines)
    _validate_dict(model.aecomps)
    _validate_dict(model.aefacts)
    #_validate_dict(model.panlists)

    _validate_dict_list(model.aelinks)

    _validate_dict(model.aeparams)
    _validate_dict(model.aesurf)
    _validate_dict(model.aesurfs)
    _validate_dict(model.aestats)
    _validate_dict(model.trims)

    _validate_dict(model.divergs)
    _validate_dict(model.csschds)
    _validate_list(model.mkaeros)
    _validate_list(model.monitor_points)

    #------------------------------------------------
    if model.aero is not None:
        model.aero.validate()

    _validate_dict(model.flfacts)
    _validate_dict(model.flutters)
    _validate_dict(model.gusts)

    #------------------------------------------------
    _validate_dict_list(model.bcs)
    _validate_dict(model.phbdys)
    _validate_dict(model.convection_properties)
    _validate_dict(model.tempds)
    #------------------------------------------------
    _validate_dict(model.bcrparas)
    _validate_dict(model.bctadds)
    _validate_dict(model.bctparas)
    _validate_dict(model.bctsets)
    _validate_dict(model.bsurf)
    _validate_dict(model.bsurfs)

    #------------------------------------------------
    _validate_dict(model.suport1)
    _validate_list(model.suport)
    _validate_list(model.se_suport)

    _validate_dict_list(model.spcadds)
    _validate_dict_list(model.spcs)
    _validate_dict_list(model.mpcadds)
    _validate_dict_list(model.mpcs)
    _validate_dict_list(model.spcoffs)

    #------------------------------------------------
    _validate_dict(model.dareas)
    _validate_dict(model.dphases)

    _validate_dict(model.pbusht)
    _validate_dict(model.pdampt)
    _validate_dict(model.pelast)

    _validate_dict_list(model.frequencies)
    #------------------------------------------------
    _validate_dict(model.dmis)
    _validate_dict(model.dmigs)
    _validate_dict(model.dmijs)
    _validate_dict(model.dmijis)
    _validate_dict(model.dmiks)
    #------------------------------------------------
    #model.asets = []
    #model.bsets = []
    #model.csets = []
    #model.qsets = []
    #model.usets = {}

    ##: SExSETy
    #model.se_bsets = []
    #model.se_csets = []
    #model.se_qsets = []
    #model.se_usets = {}
    #model.se_sets = {}

    _validate_dict(model.sets)
    _validate_dict_list(model.usets)

    _validate_list(model.asets)
    _validate_list(model.omits)
    _validate_list(model.bsets)
    _validate_list(model.csets)
    _validate_list(model.qsets)

    _validate_dict(model.se_sets)
    _validate_dict(model.se_usets)

    _validate_list(model.se_bsets)
    _validate_list(model.se_csets)
    _validate_list(model.se_qsets)
    #------------------------------------------------
    _validate_dict(model.tables)
    _validate_dict(model.tables_d)
    _validate_dict(model.tables_m)
    _validate_dict(model.random_tables)
    _validate_dict(model.tables_sdamping)
    _validate_dict(model.sets)

    #------------------------------------------------
    _validate_dict(model.methods)
    _validate_dict(model.cMethods)
    #------------------------------------------------
    _validate_dict(model.dconadds)
    _validate_dict_list(model.dconstrs)

    _validate_dict(model.desvars)
    _validate_dict(model.ddvals)
    _validate_dict(model.dlinks)
    _validate_dict(model.dresps)

    if model.dtable is not None:
        model.dtable.validate()
    if model.doptprm is not None:
        model.doptprm.validate()

    _validate_dict(model.dequations)
    _validate_dict(model.dvprels)
    _validate_dict(model.dvmrels)
    _validate_dict(model.dvcrels)
    for unused_key, dscreen in sorted(model.dscreen.items()):
        dscreen.validate()
    _validate_dict_list(model.dvgrids)
