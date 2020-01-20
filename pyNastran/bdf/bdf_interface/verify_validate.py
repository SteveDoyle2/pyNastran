"""
defines:
 - verify_bdf(model, xref)
 - validate_bdf(model)

"""
from __future__ import annotations
import sys
import traceback
from typing import List, Dict, Tuple, Any, TYPE_CHECKING
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def verify_bdf(model: BDF, xref: bool) -> None:
    #for key, card in sorted(model.params.items()):
        #card._verify(xref)
    for unused_key, card in sorted(model.nodes.items()):
        try:
            card._verify(xref)
        except:
            print(str(card))
            raise

    _verify_dict(model.coords, xref)
    for unused_key, card in sorted(model.elements.items()):
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

    _verify_dict(model.dequations, xref)
    _verify_dict(model.desvars, xref)
    _verify_dict(model.topvar, xref)
    _verify_dict(model.dvcrels, xref)
    _verify_dict(model.dvmrels, xref)
    _verify_dict(model.dvprels, xref)
    _verify_model_dict(model.dresps, model, xref)
    _verify_dict_list(model.dvgrids, xref)

    for unused_id, gust in sorted(model.gusts.items()):
        gust._verify(model, xref)
    _verify_dict(model.tics, xref)
    model.zona.verify(xref)

    for unused_super_id, superelement in model.superelement_models.items():
        verify_bdf(superelement, xref)

def _verify_dict(dict_obj: Dict[Any, Any], xref: bool) -> None:
    """helper for ``verify_bdf``"""
    for unused_key, card in sorted(dict_obj.items()):
        try:
            card._verify(xref)
        except:
            print(str(card))
            raise

def _verify_model_dict(dict_obj: Dict[Any, Any], model: BDF, xref: bool) -> None:
    """helper for ``verify_bdf``"""
    for unused_key, card in sorted(dict_obj.items()):
        try:
            card._verify(model, xref)
        except:
            print(str(card))
            raise

def _verify_dict_list(dict_list: Dict[Any, List[Any]], xref: bool) -> None:
    """helper for ``verify_bdf``"""
    for unused_key, cards in sorted(dict_list.items()):
        for card in cards:
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise

def _print_card(card: Any) -> str:
    """helper for ``_validate_msg``"""
    try:
        return card.write_card(size=8)
    except RuntimeError:
        return ''

def _validate_msg(card_obj: Any) -> str:
    """helper for ``_validate_traceback``"""
    msg = traceback.format_exc()
    try:
        msg += ('\n' + card_obj.get_stats() + '\n' + _print_card(card_obj)).rstrip()
    except KeyError:
        msg += '\n' + card_obj.get_stats().rstrip()
    return msg

def _validate_traceback(model: BDF, obj, unused_error,
                        ifailed: int, nmax_failed: int) -> Tuple[int, Any, Any, Any]:
    """helper method for ``validate_bdf`` to write a traceback"""
    exc_type, exc_value, exc_traceback = sys.exc_info()
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
    model.log.error(_validate_msg(obj))
    ifailed += 1
    if ifailed > nmax_failed:
        # PY3: raise error from None
        raise
    return ifailed, exc_type, exc_value, exc_traceback

def validate_bdf(model: BDF) -> None:
    _validate_dict(model, model.nodes)
    _validate_dict(model, model.points)
    _validate_dict(model, model.coords)
    _validate_dict(model, model.elements)
    _validate_dict(model, model.properties)
    _validate_dict(model, model.rigid_elements)
    _validate_dict(model, model.plotels)
    _validate_dict(model, model.masses)
    _validate_dict(model, model.properties_mass)
    #------------------------------------------------
    _validate_dict(model, model.materials)
    _validate_dict(model, model.thermal_materials)
    _validate_dict(model, model.MATS1)
    _validate_dict(model, model.MATS3)
    _validate_dict(model, model.MATS8)
    _validate_dict(model, model.MATT1)
    _validate_dict(model, model.MATT2)
    _validate_dict(model, model.MATT3)
    _validate_dict(model, model.MATT4)
    _validate_dict(model, model.MATT5)
    _validate_dict(model, model.MATT8)
    _validate_dict(model, model.MATT9)
    _validate_dict(model, model.creep_materials)
    _validate_dict(model, model.hyperelastic_materials)

    #------------------------------------------------
    _validate_dict_list(model, model.load_combinations)
    _validate_dict_list(model, model.loads)
    _validate_dict_list(model, model.dloads)
    _validate_dict_list(model, model.dloads)

    #------------------------------------------------
    _validate_dict(model, model.nlpcis)
    _validate_dict(model, model.nlparms)
    _validate_dict(model, model.rotors)
    _validate_dict(model, model.tsteps)
    _validate_dict(model, model.tstepnls)
    _validate_dict_list(model, model.transfer_functions)
    _validate_dict(model, model.delays)

    #------------------------------------------------
    if model.aeros is not None:
        model.aeros.validate()
    _validate_dict(model, model.caeros)
    _validate_dict(model, model.paeros)
    _validate_dict(model, model.splines)
    _validate_dict(model, model.aecomps)
    _validate_dict(model, model.aefacts)
    #_validate_dict(model, model.panlists)

    _validate_dict_list(model, model.aelinks)

    _validate_dict(model, model.aeparams)
    _validate_dict(model, model.aesurf)
    _validate_dict(model, model.aesurfs)
    _validate_dict(model, model.aestats)
    _validate_dict(model, model.trims)

    _validate_dict(model, model.divergs)
    _validate_dict(model, model.csschds)
    _validate_list(model, model.mkaeros)
    _validate_list(model, model.monitor_points)

    #------------------------------------------------
    if model.aero is not None:
        model.aero.validate()

    _validate_dict(model, model.flfacts)
    _validate_dict(model, model.flutters)
    _validate_dict(model, model.gusts)

    #------------------------------------------------
    _validate_dict_list(model, model.bcs)
    _validate_dict(model, model.phbdys)
    _validate_dict(model, model.convection_properties)
    _validate_dict(model, model.tempds)
    #------------------------------------------------
    _validate_dict(model, model.bcrparas)
    _validate_dict(model, model.bctadds)
    _validate_dict(model, model.bctparas)
    _validate_dict(model, model.bctsets)
    _validate_dict(model, model.bsurf)
    _validate_dict(model, model.bsurfs)

    #------------------------------------------------
    _validate_dict(model, model.suport1)
    _validate_list(model, model.suport)
    _validate_list(model, model.se_suport)

    _validate_dict_list(model, model.spcadds)
    _validate_dict_list(model, model.spcs)
    _validate_dict_list(model, model.mpcadds)
    _validate_dict_list(model, model.mpcs)
    _validate_dict_list(model, model.spcoffs)

    #------------------------------------------------
    _validate_dict(model, model.dareas)
    _validate_dict(model, model.dphases)

    _validate_dict(model, model.pbusht)
    _validate_dict(model, model.pdampt)
    _validate_dict(model, model.pelast)

    _validate_dict_list(model, model.frequencies)
    #------------------------------------------------
    _validate_dict(model, model.dmi)
    _validate_dict(model, model.dmig)
    _validate_dict(model, model.dmij)
    _validate_dict(model, model.dmiji)
    _validate_dict(model, model.dmik)
    _validate_dict(model, model.dmiax)
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

    _validate_dict(model, model.sets)
    _validate_dict_list(model, model.usets)

    _validate_list(model, model.asets)
    _validate_list(model, model.omits)
    _validate_list(model, model.bsets)
    _validate_list(model, model.csets)
    _validate_list(model, model.qsets)

    _validate_dict(model, model.se_sets)
    _validate_dict(model, model.se_usets)

    _validate_list(model, model.se_bsets)
    _validate_list(model, model.se_csets)
    _validate_list(model, model.se_qsets)
    #------------------------------------------------
    _validate_dict(model, model.tables)
    _validate_dict(model, model.tables_d)
    _validate_dict(model, model.tables_m)
    _validate_dict(model, model.random_tables)
    _validate_dict(model, model.tables_sdamping)
    _validate_dict(model, model.sets)

    #------------------------------------------------
    _validate_dict(model, model.methods)
    _validate_dict(model, model.cMethods)
    #------------------------------------------------
    _validate_dict(model, model.dconadds)
    _validate_dict_list(model, model.dconstrs)

    _validate_dict(model, model.desvars)
    _validate_dict(model, model.topvar)
    _validate_dict(model, model.ddvals)
    _validate_dict(model, model.dlinks)
    _validate_dict(model, model.dresps)

    if model.dtable is not None:
        model.dtable.validate()
    if model.doptprm is not None:
        model.doptprm.validate()

    _validate_dict(model, model.dequations)
    _validate_dict(model, model.dvprels)
    _validate_dict(model, model.dvmrels)
    _validate_dict(model, model.dvcrels)
    for unused_key, dscreen in sorted(model.dscreen.items()):
        dscreen.validate()
    _validate_dict_list(model, model.dvgrids)
    model.zona.validate()

    for unused_super_id, superelement in model.superelement_models.items():
        validate_bdf(superelement)

def _validate_dict_list(model: BDF, objects_dict: Dict[Any, Any]) -> None:
    """helper method for validate_bdf"""
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
                ifailed, exc_type, exc_value, exc_traceback = _validate_traceback(
                    model, obj, error, ifailed, nmax_failed)

        if ifailed:
            raise

def _validate_dict(model: BDF, objects: Dict[Any, Any]) -> None:
    """helper method for validate_bdf"""
    assert isinstance(objects, dict), type(objects)
    ifailed = 0
    nmax_failed = 0
    for unused_id, obj in sorted(objects.items()):
        try:
            obj.validate()
        except(ValueError, AssertionError, RuntimeError, IndexError) as error:
            ifailed, exc_type, exc_value, exc_traceback = _validate_traceback(
                model, obj, error, ifailed, nmax_failed)
    if ifailed:
        raise

def _validate_list(model: BDF, objects: List[Any]) -> None:
    """helper method for validate_bdf"""
    ifailed = 0
    nmax_failed = 0
    assert isinstance(objects, list), type(objects)
    for obj in objects:
        try:
            obj.validate()
        except(ValueError, AssertionError, RuntimeError) as error:
            ifailed, exc_type, exc_value, exc_traceback = _validate_traceback(
                model, obj, error, ifailed, nmax_failed)
    if ifailed:
        raise
