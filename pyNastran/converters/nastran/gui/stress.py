# coding: utf-8
from __future__ import annotations
import sys
import traceback
from copy import deepcopy
from functools import wraps
from collections import defaultdict
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.femutils.utils import pivot_table, unique2d

from pyNastran.op2.result_objects.stress_object import _get_nastran_header
from pyNastran.op2.op2_interface.op2_classes import (
    RealBarStressArray, # ComplexBarStressArray,
    #RealSolidStressArray, RealSolidStressArrayNx, ComplexSolidStressArray,
)
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids_nx import RealSolidArrayNx
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidArray

from pyNastran.converters.nastran.gui.result_objects.simple_table_results import SimpleTableResults
from pyNastran.converters.nastran.gui.result_objects.layered_table_results import LayeredTableResults

from pyNastran.converters.nastran.gui.result_objects.plate_stress_results import PlateStrainStressResults2
from pyNastran.converters.nastran.gui.result_objects.solid_stress_results import SolidStrainStressResults2
from pyNastran.converters.nastran.gui.result_objects.composite_stress_results import CompositeStrainStressResults2
from pyNastran.converters.nastran.gui.types import CasesDict, NastranKey, KeysMap, KeyMap

if TYPE_CHECKING: # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.op2.op2 import OP2
    from pyNastran.op2.result_objects.stress_object import CompositeDict


def nocrash_log(func):
    @wraps(func)
    def wrapper(log: SimpleLogger, stop_on_failure: bool,
                cases, *args, **kwargs):
        assert isinstance(cases, dict), case
        ncases = len(cases)
        try:
            icase = func(cases, *args, **kwargs)
        except NotImplementedError:  # pragma: no cover
            raise
        except Exception as e:
            #print('dont crash...')
            traceback.print_exc(file=sys.stdout)
            log.error(str(e))
            if stop_on_failure:
                raise
            ncases2 = len(cases)
            icase += ncases2 - ncases
        return icase
    return wrapper

@nocrash_log
def get_rod_stress_strains(cases: CasesDict,
                           eids: np.ndarray,
                           model: OP2,
                           times: np.ndarray,
                           key: NastranKey,
                           icase: int,
                           form_dict, header_dict,
                           keys_map: KeysMap,
                           log: SimpleLogger,
                           is_stress: bool) -> int:
    """
    helper method for _fill_op2_time_centroidal_stress.
    """
    subcase_id = key[0]
    if is_stress:
        rods = [
            model.crod_stress, model.conrod_stress, model.ctube_stress,
        ]
        word = 'Stress'
    else:
        rods = [
            model.crod_strain, model.conrod_strain, model.ctube_strain,
        ]
        word = 'Strain'

    #titles = []
    rod_cases = []
    rod_ieids = []
    for result in rods:
        if key not in result:
            continue
        case = result[key]

        eidsi = case.element
        i = np.searchsorted(eids, eidsi)
        if len(i) != len(np.unique(i)):
            #print(case.element_node)
            #print('element_name=%s nnodes_per_element=%s' % (case.element_name, nnodes_per_element))
            #print('iplate = %s' % i)
            #print('  eids = %s' % eids)
            #print('  eidsiA = %s' % case.element_node[:, 0])
            #print('  eidsiB = %s' % eidsi)
            msg = 'i%s (rod)=%s is not unique' % (case.element_name, str(i))
            #msg = 'iplate=%s is not unique' % str(i)
            log.warning(msg)
            continue
        if i.max() == len(eids):
            log.error('skipping because lookup is out of range...')
            continue
        #print('i =', i, i.max())
        #print('eids =', eids, len(eids))
        #print('eidsi =', eidsi, eids)
        #print(f'------------adding i={i} for {case.element_name}-----------')
        rod_cases.append(case)
        rod_ieids.append(i)
    if not rod_ieids:
        return icase

    rod_ieids = np.hstack(rod_ieids)
    ieid_max = len(eids)
    #print('ieid_max =', ieid_max)

    case = rod_cases[0]
    case_headers = case.get_headers()
    #print(case_headers)
    if is_stress:
        #sigma = 'σ'
        method_map = {
            'axial' : 'Stress XX',
            'torsion' : 'Shear XY',
            'SMa' : 'MS_axial',
            'SMt' : 'MS_torsion',
            #'omax' : 'σmax',
            #'omin' : 'σmin',
            #'von_mises' : 'σ von Mises',
        }
        data_format = '%.3f'
    else:
        #sigma = 'ϵ'
        method_map = {
            'axial' : 'Strain XX',
            'torsion' : 'Shear XY',
            'SMa' : 'MS_axial',
            'SMt' : 'MS_torsion',
            #'emax' : 'ϵmax',
            #'emin' : 'ϵmin',
            #'von_mises' : 'ϵ von Mises',
        }
        data_format = '%.3e'
    methods = [method_map[headeri] for headeri in case_headers]
    #if 'Mises' in methods:
        #methods.append('Max shear')
    #else:
        #methods.append(f'{sigma} von Mises')

    #if case.is_von_mises:
        #vm_word = 'vonMises'
    #else:
        #vm_word = 'maxShear'

    scalars_array = []
    for case in rod_cases:
        keys_map[key] = KeyMap(case.subtitle, case.label,
                               case.superelement_adaptivity_index,
                               case.pval_step)
        scalars = case.data
        scalars_array.append(scalars)

    if len(scalars_array) == 0:
        return icase

    scalars_array = concatenate_scalars(scalars_array)
    headers = [] # sidebar word
    res = SimpleTableResults(
        subcase_id, headers, rod_ieids, ieid_max, scalars_array, methods,
        data_format=data_format,
        colormap='jet', uname='Rod ' + word)

    icase = add_simple_methods_to_form(icase, cases, key, subcase_id, word, res, case,
                                       form_dict, header_dict, methods,
                                       name='Rod')
    return icase

@nocrash_log
def get_bar_stress_strains(cases: CasesDict,
                           eids: np.ndarray,
                           model: OP2,
                           times: np.ndarray,
                           key: NastranKey,
                           icase: int,
                           form_dict, header_dict,
                           keys_map: KeysMap,
                           log: SimpleLogger,
                           is_stress: bool) -> int:
    """
    helper method for _fill_op2_time_centroidal_stress.
    """
    #print("***stress eids=", eids)
    subcase_id = key[0]
    if is_stress:
        bars = [model.cbar_stress]
        word = 'Stress'
    else:
        bars = [model.cbar_strain]
        word = 'Strain'

    #titles = []
    bar_cases = []
    bar_ieids = []
    for result in bars:
        if key not in result:
            continue
        case = result[key]
        if case.is_complex:
            log.warning(f'skipping {case.class_name}')
            continue

        eidsi = case.element
        i = np.searchsorted(eids, eidsi)
        if len(i) != len(np.unique(i)):
            #print(case.element_node)
            #print('element_name=%s nnodes_per_element=%s' % (case.element_name, nnodes_per_element))
            #print('iplate = %s' % i)
            #print('  eids = %s' % eids)
            #print('  eidsiA = %s' % case.element_node[:, 0])
            #print('  eidsiB = %s' % eidsi)
            msg = 'i%s (bar)=%s is not unique' % (case.element_name, str(i))
            #msg = 'iplate=%s is not unique' % str(i)
            log.warning(msg)
            continue
        if i.max() == len(eids):
            log.error('skipping because lookup is out of range...')
            continue
        #print('i =', i, i.max())
        #print('eids =', eids, len(eids))
        #print('eidsi =', eidsi, eids)
        #print(f'------------adding i={i} for {case.element_name}-----------')
        bar_cases.append(case)
        bar_ieids.append(i)
    if not bar_ieids:
        return icase

    bar_ieids = np.hstack(bar_ieids)
    ieid_max = len(eids)
    #print('ieid_max =', ieid_max)

    case: RealBarStressArray = bar_cases[0]
    #assert isinstance(case, (ComplexBarStressArray, RealBarStressArray)), case
    case_headers = case.get_headers()
    #print(case_headers)

    # real
    #complex:
    # [s1a, s1b, s1c, s1d, axial,
    #  s2a, s2b, s2c, s2d, ]
    if is_stress:
        method_map = {
             's1a' : 'Stress 1A',
             's2a' : 'Stress 2A',
             's3a' : 'Stress 3A',
             's4a' : 'Stress 4A',

             's1b' : 'Stress 1B',
             's2b' : 'Stress 2B',
             's3b' : 'Stress 3B',
             's4b' : 'Stress 4B',

            'axial' : 'Stress XX',
            'MS_tension' : 'MS_tension',
            'MS_compression' : 'MS_compression',
            'smaxa' : 'Max Principal A',
            'smina' : 'Min Principal A',
            'smaxb' : 'Max Principal B',
            'sminb' : 'Min Principal B',
            #'von_mises' : 'σ von Mises',
        }
        data_format = '%.3f'
    else:
        method_map = {
            'e1a' : 'Strain 1A',
            'e2a' : 'Strain 2A',
            'e3a' : 'Strain 3A',
            'e4a' : 'Strain 4A',

            'e1b' : 'Strain 1B',
            'e2b' : 'Strain 2B',
            'e3b' : 'Strain 3B',
            'e4b' : 'Strain 4B',

            'e1c': 'Strain 1C',
            'e1d': 'Strain 1D',
            'e2c': 'Strain 2C',
            'e2d': 'Strain 2D',

           'axial' : 'Strain XX',
           'MS_tension' : 'MS_tension',
           'MS_compression' : 'MS_compression',
           'emaxa' : 'Max Principal A',
           'emina' : 'Min Principal A',
           'emaxb' : 'Max Principal B',
           'eminb' : 'Min Principal B',
            #'emax' : 'ϵmax',
            #'emin' : 'ϵmin',
            #'von_mises' : 'ϵ von Mises',
        }
        data_format = '%.3e'
    methods = [method_map[headeri] for headeri in case_headers]
    #if 'Mises' in methods:
        #methods.append('Max shear')
    #else:
        #methods.append(f'{sigma} von Mises')

    #if case.is_von_mises:
        #vm_word = 'vonMises'
    #else:
        #vm_word = 'maxShear'

    #headersi = case.get_headers()
    #print('headersi =', headersi)

    scalars_array = []
    for case in bar_cases:
        #ntimes, nelements, nresults = case.data.shape
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
        #                                         txy, angle,
        #                                         majorP, minorP, ovm]

        keys_map[key] = KeyMap(case.subtitle, case.label,
                               case.superelement_adaptivity_index,
                               case.pval_step)

        #nnodes_per_element = case.nnodes
        #nelements_nnodes = nnodes_nlayers // 2
        #nelements = nelements_nnodes // nnodes_per_element
        #nlayers = 2
        scalars = case.data
        scalars_array.append(scalars)

    if len(scalars_array) == 0:
        return icase

    scalars_array = concatenate_scalars(scalars_array)

    exclude_tension_margin = False
    exclude_compression_margin = False
    if 'MS_tension' in methods:
        # real
        itension = methods.index('MS_tension')
        icompression = methods.index('MS_compression')
        exclude_tension_margin = np.allclose(np.abs(scalars_array[:, :, itension]).max(), 0.0)
        exclude_compression_margin = np.allclose(np.abs(scalars_array[:, :, icompression]).max(), 0.0)
        iresults = list(range(len(methods)))
        if exclude_compression_margin:
            methods.pop(icompression)
            iresults.pop(icompression)
        if exclude_tension_margin:
            methods.pop(itension)
            iresults.pop(itension)

    if exclude_compression_margin or exclude_compression_margin:
        scalars_array = scalars_array[:, :, iresults]

    #titles = []  # legend title
    headers = [] # sidebar word
    res = SimpleTableResults(
        subcase_id, headers, bar_ieids, ieid_max, scalars_array, methods,
        data_format=data_format,
        colormap='jet', uname='Bar ' + word)

    icase = add_simple_methods_to_form(icase, cases, key, subcase_id, word, res, case,
                                       form_dict, header_dict, methods,
                                       name='Bar')
    return icase

@nocrash_log
def get_beam_stress_strains(cases: CasesDict,
                            eids: np.ndarray,
                            model: OP2,
                            times: np.ndarray,
                            key: NastranKey,
                            icase: int,
                            form_dict, header_dict,
                            keys_map: KeysMap,
                            log: SimpleLogger,
                            is_stress: bool) -> int:
    """
    helper method for _fill_op2_time_centroidal_stress.
    """
    #print("***stress eids=", eids)
    subcase_id = key[0]
    if is_stress:
        beams = [model.cbeam_stress]
        word = 'Stress'
    else:
        beams = [model.cbeam_strain]
        word = 'Strain'

    #titles = []
    beam_cases = []
    beam_ieids = []
    for result in beams:
        if key not in result:
            continue
        case = result[key]

        eidsi = np.unique(case.element_node[:, 0])
        print(eidsi)
        i = np.searchsorted(eids, eidsi)
        if len(i) != len(np.unique(i)):
            #print(case.element_node)
            #print('element_name=%s nnodes_per_element=%s' % (case.element_name, nnodes_per_element))
            #print('iplate = %s' % i)
            #print('  eids = %s' % eids)
            #print('  eidsiA = %s' % case.element_node[:, 0])
            #print('  eidsiB = %s' % eidsi)
            msg = 'i%s (beam)=%s is not unique' % (case.element_name, str(i))
            #msg = 'iplate=%s is not unique' % str(i)
            log.warning(msg)
            continue
        if i.max() == len(eids):
            log.error('skipping because lookup is out of range...')
            continue

        i2 = np.hstack([i, i + 1]).T.flatten()
        #print('i =', i)
        #print('i2 =', i2)
        #aa
        #print('i =', i, i.max())
        #print('eids =', eids, len(eids))
        #print('eidsi =', eidsi, eids)
        #print(f'------------adding i={i} for {case.element_name}-----------')
        beam_cases.append(case)
        beam_ieids.append(i2)
    if not beam_ieids:
        return icase

    beam_ieids = np.hstack(beam_ieids)
    #inid_max = len(nids)
    ieid_max = len(eids)
    #print('ieid_max =', ieid_max)

    case = beam_cases[0]
    case_headers = case.get_headers()
    #print(case_headers)
    if is_stress:
        #sigma = 'σ'
        method_map = {
            'sxc' : 'Stress C',
            'sxd' : 'Stress D',
            'sxe' : 'Stress E',
            'sxf' : 'Stress F',
            #'torsion' : 'τxy',
            'MS_tension' : 'MS_tension',
            'MS_compression' : 'MS_compression',
            'smax' : 'Max Principal',
            'smin' : 'Min Principal',
            #'von_mises' : 'von Mises',
        }
    else:
        #sigma = 'ϵ'
        method_map = {
            'sxc' : 'Strain C',
            'sxd' : 'Strain D',
            'sxe' : 'Strain E',
            'sxf' : 'Strain F',

            'exc' : 'Strain C',
            'exd' : 'Strain D',
            'exe' : 'Strain E',
            'exf' : 'Strain F',
            #'torsion' : 'τxy',
            'MS_tension' : 'MS_tension',
            'MS_compression' : 'MS_compression',
            'smax' : 'Max Principal',
            'smin' : 'Min Principal',
            #'von_mises' : 'ϵ von Mises',
        }
    methods = [method_map[headeri] for headeri in case_headers]
    return icase
    #if 'Mises' in methods:
        #methods.append('Max shear')
    #else:
        #methods.append(f'{sigma} von Mises')

    #if case.is_von_mises:
        #vm_word = 'vonMises'
    #else:
        #vm_word = 'maxShear'

    #headersi = case.get_headers()
    #print('headersi =', headersi)

    scalars_array = []
    for case in beam_cases:
        if case.is_complex:
            log.warning(f'skipping complex Beam {word}')
            continue

        #ntimes, nelements, nresults = case.data.shape
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
        #                                         txy, angle,
        #                                         majorP, minorP, ovm]

        keys_map[key] = KeyMap(case.subtitle, case.label,
                               case.superelement_adaptivity_index,
                               case.pval_step)

        #nnodes_per_element = case.nnodes
        #nelements_nnodes = nnodes_nlayers // 2
        #nelements = nelements_nnodes // nnodes_per_element
        #nlayers = 2
        scalars = case.data
        scalars_array.append(scalars)

    if len(scalars_array) == 0:
        return icase

    scalars_array = concatenate_scalars(scalars_array)

    #titles = []  # legend title
    headers = [] # sidebar word
    res = SimpleTableResults(
        subcase_id, headers, beam_ieids, ieid_max, scalars_array, methods,
        data_formats=None,
        colormap='jet', location='node', uname='Beam ' + word)

    times = case._times
    for itime, dt in enumerate(times):
        #dt = case._times[itime]
        header = _get_nastran_header(case, dt, itime)
        header_dict[(key, itime)] = header

        formi = []
        form = form_dict[(key, itime)]
        form.append(('Beam ' + word, None, formi))
        # formi = form[0][2]

        for imethod, method in enumerate(methods):
            #cases[icase] = (res, (subcase_id, header))
            cases[icase] = (res, (subcase_id, (itime, imethod, header)))
            formi.append((method, icase, []))
            icase += 1
    return icase

@nocrash_log
def get_plate_stress_strains(cases: CasesDict,
                             eids: np.ndarray,
                             model: OP2,
                             times: np.ndarray,
                             key: NastranKey,
                             icase: int,
                             form_dict, header_dict,
                             keys_map: KeysMap,
                             log: SimpleLogger,
                             use_old_sidebar_objects: bool,
                             is_stress: bool,
                             prefix: str='') -> int:
    """
    helper method for _fill_op2_time_centroidal_stress.
    Gets the max/min stress for each layer.

    key : varies
        (1, 5, 1, 0, 0, '', '')
        (isubcase, analysis_code, sort_method, count, ogs,
         superelement_adaptivity_index, pval_step)
    """
    if not use_old_sidebar_objects:  # pragma: no cover
        return icase
    plates, word, subcase_id, analysis_code = _get_plates(model, key, is_stress, prefix)
    word += ' (centroid)'

    #analysis_code = key[1]
    ##print("***stress eids=", eids)
    #subcase_id = key[0]
    #if prefix == 'modal_contribution':
        #results = model.op2_results.modal_contribution
        #preword = 'Modal Contribution '
    #elif prefix == '':
        #results = model
        #preword = ''
    #else:  # pragma: no cover
        #raise NotImplementedError(prefix)

    #if is_stress:
        #plates = [
            #results.ctria3_stress, results.cquad4_stress,
            #results.ctria6_stress, results.cquad8_stress,
            #results.ctriar_stress, results.cquadr_stress, # results.cquad_stress,
        #]
        #word = preword + 'Stress (centroid)'
    #else:
        #plates = [
            #results.ctria3_strain, results.cquad4_strain,
            #results.ctria6_strain, results.cquad8_strain,
            #results.ctriar_strain, results.cquadr_strain, # results.cquad_strain,
        #]
        #word = preword + 'Strain (centroid)'

    #titles = []
    plate_cases = []
    plates_ieids = []
    #print('key =', key)
    for iplate, result in enumerate(plates):
        #if result:
        #print(f'keys[{iplate}] = {result.keys()}')
        if key not in result:
            continue
        case = result[key]

        if analysis_code in [5, 9]:  # complex
            # 5-freq
            # 9-complex eigenvalues
            if not case.is_complex:
                log.info(f'skipping:\n{case}{case.code_information()}')
                continue
        else:
            assert analysis_code in [1, 2, 6, 7, 8, 10], case.code_information()
            # 1-statics
            # 2-modes
            # 6-transient
            # 7-pre buckling
            # 8-post buckling
            # 10-nonlinear statics
            if not case.is_real:
                log.info(f'skipping:\n{case}{case.code_information()}')
                continue

        nnodes_per_element = case.nnodes_per_element
        nlayers_per_element = nnodes_per_element * 2
        eidsi = case.element_node[::nlayers_per_element, 0]  # ::2 is for layer skipping
        #print(case.element_name, eidsi)
        i = np.searchsorted(eids, eidsi)
        if len(i) != len(np.unique(i)):
            #print(case.element_node)
            #print('element_name=%s nnodes_per_element=%s' % (case.element_name, nnodes_per_element))
            #print('iplate = %s' % i)
            #print('  eids = %s' % eids)
            #print('  eidsiA = %s' % case.element_node[:, 0])
            #print('  eidsiB = %s' % eidsi)
            msg = 'i%s (plate)=%s is not unique' % (case.element_name, str(i))
            #msg = 'iplate=%s is not unique' % str(i)
            log.warning(msg)
            continue
        if i.max() == len(eids):
            log.error('skipping because lookup is out of range...')
            continue
        #print('i =', i, i.max())
        #print('eids =', eids, len(eids))
        #print('eidsi =', eidsi, eids)
        #print(f'------------adding i={i} for {case.element_name}-----------')
        plate_cases.append(case)
        plates_ieids.append(i)
    if not plates_ieids:
        return icase

    plates_ieids = np.hstack(plates_ieids)
    ieid_max = len(eids)
    #print('ieid_max =', ieid_max)

    case = plate_cases[0]
    case_headers = case.get_headers()
    #print(case_headers)
    if is_stress:
        method_map = {
            'fiber_curvature' : 'FiberCurvature',
            'fiber_distance' : 'FiberDistance',
            'oxx' : 'σxx',
            'oyy' : 'σyy',
            'txy' : 'τxy',
            'angle' : 'θ',
            'omax' : 'σmax',
            'omin' : 'σmin',
            'von_mises' : 'σ von Mises',
            'max_shear' : 'τmax',
        }
    else:
        method_map = {
            'fiber_curvature' : 'FiberCurvature',
            'fiber_distance' : 'FiberDistance',
            'exx' : 'ϵxx',
            'eyy' : 'ϵyy',
            'exy' : 'ϵxy',
            'angle' : 'θ',
            'emax' : 'ϵmax',
            'emin' : 'ϵmin',
            'evm' : 'ϵ von Mises',
            'von_mises' : 'ϵ von Mises',
            'max_shear' : '𝛾max',
        }
    methods = [method_map[headeri] for headeri in case_headers]
    #if 'Mises' in methods:
        #methods.append('Max shear')
    #else:
        #methods.append(f'{sigma} von Mises')

    #if case.is_von_mises:
        #vm_word = 'vonMises'
    #else:
        #vm_word = 'maxShear'

    #headersi = case.get_headers()
    #print('headersi =', headersi)

    scalars_array = []
    for case in plate_cases:
        #if case.is_complex:
            #print(f"skipping {case.class_name} because it's complex")
            #model.log.warning(f"skipping {case.class_name} because it's complex")
            #continue

        ntimes, nnodes_nlayers, nresults = case.data.shape
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
        #                                         txy, angle,
        #                                         majorP, minorP, ovm]
        keys_map[key] = KeyMap(case.subtitle, case.label,
                               case.superelement_adaptivity_index,
                               case.pval_step)

        nnodes_per_element = case.nnodes_per_element
        nelements_nnodes = nnodes_nlayers // 2
        nelements = nelements_nnodes // nnodes_per_element
        nlayers = 2
        scalars = case.data.reshape(ntimes, nelements, nnodes_per_element, nlayers, nresults)
        scalars_array.append(scalars[:, :, 0, :, :])

    if len(scalars_array) == 0:
        return icase

    element_name = case.element_name
    scalars_array = concatenate_scalars(scalars_array)

    #titles = []  # legend title
    headers = [] # sidebar word
    res = LayeredTableResults(
        subcase_id, headers, plates_ieids, ieid_max, scalars_array, methods,
        data_formats=None,
        colormap='jet', uname='Plate ' + word)

    if 'fiber_curvature' in case_headers:
        layer_names = {
            0 : 'Layer 1 (Mean)',
            1 : 'Layer 2 (Slope)',
        }
    else:
        layer_names = {
            0 : 'Layer 1 (Upper)',
            1 : 'Layer 2 (Lower)',
        }

    times = case._times

    form_names = []
    for itime, dt in enumerate(times):
        #dt = case._times[itime]
        header = _get_nastran_header(case, dt, itime)
        header2 = header.replace(' = ', '=')
        header_dict[(key, itime)] = header

        formi = []
        form = form_dict[(key, itime)]
        form.append(('Plate ' + word, None, formi))
        # formi = form[0][2]
        form_dict[(key, itime)] = form

        form_namesi = []
        for ilayer in range(2):
            layer = layer_names[ilayer]
            form_layeri = []
            formi.append((layer, None, form_layeri))

            form_namesii = []
            for imethod, method in enumerate(methods):
                #cases[icase] = (res, (subcase_id, header))
                cases[icase] = (res, (subcase_id, (itime, ilayer, imethod, header)))
                form_layeri.append((f'{method} ({layer})', icase, []))
                form_name2 = f'{element_name} Plate {word}: {method} ({layer}) {header2}'.rstrip()
                form_namesii.append(form_name2)
                icase += 1
            form_namesi.append(form_namesii)
        form_names.append(form_namesi)
    res.form_names = np.array(form_names)
    return icase

@nocrash_log
def get_plate_stress_strains2(cases: CasesDict,
                              node_id: np.ndarray,
                              element_id: np.ndarray,
                              model: OP2,
                              times: np.ndarray,
                              key: NastranKey,
                              icase: int,
                              form_dict, header_dict,
                              keys_map: KeysMap,
                              eid_to_nid_map: dict[int, list[int]],
                              log: SimpleLogger,
                              use_new_sidebar_objects: bool,
                              use_new_terms: bool,
                              is_stress: bool,
                              prefix: str='') -> int:
    """
    helper method for _fill_op2_time_centroidal_stress.
    Gets the max/min stress for each layer.

    key : varies
        (1, 5, 1, 0, 0, '', '')
        (isubcase, analysis_code, sort_method, count, ogs,
         superelement_adaptivity_index, pval_step)
    """
    if not use_new_sidebar_objects:  # pragma: no cover
        return icase
    plates, word, subcase_id, analysis_code = _get_plates(model, key, is_stress, prefix)

    plate_cases = []
    for iplate, result in enumerate(plates):
        #if result:
        #print(f'keys[{iplate}] = {result.keys()}')
        if key not in result:
            continue
        plate_case = result[key]
        if plate_case.is_complex:
            continue
        eids = np.unique(plate_case.element_node[:, 0])
        common_eids = np.intersect1d(element_id, eids)
        if len(common_eids) == 0:
            continue
        #ieids = np.unique(np.searchsorted(element_id, eids))
        #if element_id[ieids[0]] != eids[0]
        plate_cases.append(plate_case)

    if len(plate_cases) == 0:
        return icase

    #plate_case_headers = plate_case.get_headers()
    is_von_mises = plate_case.is_von_mises
    assert isinstance(is_von_mises, bool), is_von_mises
    von_mises = 7 if is_von_mises else 'von_mises'
    max_shear = 7 if not is_von_mises else 'max_shear'
    #print(case_headers)
    if is_stress:
        # [fiber_dist, 'oxx', 'oyy', 'txy', 'angle', 'omax', 'omin', ovm]
        iresult_to_title_annotation_map = {
            # iresult: (sidebar_label, annotation)
            0 : ('FiberDistance', 'Fiber Distance'),
            1 : ('Stress XX', 'XX'),
            2 : ('Stress YY', 'YY'),
            3 : ('Stress XY', 'XY'),
            4 : ('Theta', 'Theta'),
            5 : ('sMax Principal', 'Max Principal'),
            6 : ('sMin Principal', 'Min Principal'),
            'abs_principal' : ('sAbs Principal', 'Abs Principal'),
            von_mises : ('von Mises', 'von Mises'), # the magnitude is large
            max_shear : ('Max Shear', 'Max Shear'), # the magnitude is large
        }
        word = 'Stress'
    else:
        iresult_to_title_annotation_map = {
            # iresult: (sidebar_label, annotation)
            #'fiber_curvature' : 'FiberCurvature',
            0 : ('FiberDistance', 'Fiber Distance'),
            1 : ('Strain XX', 'XX'),
            2 : ('Strain YY', 'YY'),
            3 : ('Strain XY', 'XY'),
            4 : ('Theta', 'Theta'),
            5 : ('eMax Principal', 'Max Principal'),
            6 : ('eMin Principal', 'Min Principal'),
            'abs_principal' : ('eAbs Principal', 'Abs Principal'),
            von_mises : ('von Mises', 'von Mises'),  # the magnitude is small
            max_shear : ('Max Shear', 'Max Shear'),  # the magnitude is small
        }
        word = 'Strain'

    if not use_new_terms:
        del iresult_to_title_annotation_map[max_shear]
        del iresult_to_title_annotation_map['abs_principal']

    title = f'Plate {word}'
    keys_map[key] = KeyMap(plate_case.subtitle, plate_case.label,
                           plate_case.superelement_adaptivity_index,
                           plate_case.pval_step)

    res = PlateStrainStressResults2(
        subcase_id, model,
        node_id, element_id,
        plate_cases, iresult_to_title_annotation_map, title,
        data_format='%g', is_variable_data_format=False,
        nlabels=None, labelsize=None, ncolors=None, colormap='',
        set_max_min=False,
        is_fiber_distance=plate_case.is_fiber_distance,
        eid_to_nid_map=eid_to_nid_map,
        uname='PlateStressStrainResults2')

    times = plate_case._times
    for itime, dt in enumerate(times):
        #dt = case._times[itime]
        header = _get_nastran_header(plate_case, dt, itime)
        header2 = header.replace(' = ', '=')
        header_dict[(key, itime)] = header

        formi = []
        form = form_dict[(key, itime)]
        form.append(('Plate ' + word, None, formi))
        form_dict[(key, itime)] = form

        for iresult, (sidebar_label, annotation_label) in iresult_to_title_annotation_map.items():
            cases[icase] = (res, (subcase_id, (itime, iresult, header2)))
            formi.append((sidebar_label, icase, []))
            icase += 1
    return icase

def _get_plates(model: OP2,
                key,
                is_stress: bool,
                prefix: str) -> tuple[str, list, int, int]:
    analysis_code = key[1]
    #print("***stress eids=", eids)
    subcase_id = key[0]
    if prefix == 'modal_contribution':
        results = model.op2_results.modal_contribution
        preword = 'Modal Contribution '
    elif prefix == '':
        results = model
        preword = ''
    else:  # pragma: no cover
        raise NotImplementedError(prefix)

    if is_stress:
        plates = [
            results.ctria3_stress, results.cquad4_stress,
            results.ctria6_stress, results.cquad8_stress,
            results.ctriar_stress, results.cquadr_stress, # results.cquad_stress,
        ]
        word = preword + 'Stress'
    else:
        plates = [
            results.ctria3_strain, results.cquad4_strain,
            results.ctria6_strain, results.cquad8_strain,
            results.ctriar_strain, results.cquadr_strain, # results.cquad_strain,
        ]
        word = preword + 'Strain'

    plates2 = [plate for plate in plates if len(plates)]
    return plates2, word, subcase_id, analysis_code

def _stack_composite_results(model: OP2, log: SimpleLogger,
                             is_stress: bool,
                             key=None):
    if is_stress:
        case_map = {
            # element_name
            'CTRIA3' : model.ctria3_composite_stress,
            'CTRIA6' : model.ctria6_composite_stress,
            'CTRIAR' : model.ctriar_composite_stress,
            'CQUAD4' : model.cquad4_composite_stress,
            'CQUAD8' : model.cquad8_composite_stress,
            'CQUADR' : model.cquadr_composite_stress,
        }
    else:
        case_map = {
            # element_name
            'CTRIA3' : model.ctria3_composite_strain,
            'CTRIA6' : model.ctria6_composite_strain,
            'CTRIAR' : model.ctriar_composite_strain,
            'CQUAD4' : model.cquad4_composite_strain,
            'CQUAD8' : model.cquad8_composite_strain,
            'CQUADR' : model.cquadr_composite_strain,
        }

    keys_map = {}
    key_cases = defaultdict(list)
    for res_cases in case_map.values():
        for case_key, case in res_cases.items():
            if key is None:
                continue

            if case.table_name_str in {'OESCP', 'OESTRCP'}:
                log.warning(f'skipping strength ratio {case.table_name_str}')
                continue
            key_cases[key].append(case)
            if (case_key != key or key in keys_map):
                continue
            keys_map[key] = KeyMap(case.subtitle, case.label,
                                   case.superelement_adaptivity_index,
                                   case.pval_step)

    cases2 = {}
    key_cases = dict(key_cases)
    for key_, cases_ in key_cases.items():
        case = deepcopy(cases_[0])
        cases = cases_[1:]
        #nelements = case.nelements
        #ntotal = case.ntotal
        #ntimes = case.ntimes
        #nresults = case.data.shape[2]
        #itotal0 = 0
        #itotal1 = ntotal

        if len(cases):
            nelements_all = case.nelements + sum([casei.nelements for casei in cases])
            ntotal_all = case.ntotal + sum([casei.ntotal for casei in cases])
            #data2 = np.full((ntimes, ntotal, nresults), np.nan, dtype=case.data.dtype)
            #for itime in range(ntimes):
                #data2[itime, itotal0:itotal1, :] = case.data[itime, :, :]

            #element_ids = [case.element_id]
            element_layers = [case.element_layer]
            datas = [case.data]
            eids_list = [np.unique(case.element_layer[:, 0])]
            for casei in cases:
                #itotal0 += case.ntotal
                #itotal1 += case.ntotal
                eidsi = np.unique(casei.element_layer[:, 0])
                element_layers.append(casei.element_layer)
                datas.append(casei.data)
                eids_list.append(eidsi)
                #element_ids.append(case.element_id)
                #data2[itime, itotal0:itotal1, :] = casei.data[itime, :, :]

            eids = np.hstack(eids_list)
            ueids = np.unique(eids)
            if len(eids) != len(ueids):  # pragma: no cover
                msg = f'eids = {eids}\n'
                for casei_ in cases_:
                    msg += str(casei_)
                    msg += f'  eids = {np.unique(casei_.element_layer[:, 0])}\n'
                    msg += f'  element_layer:\n{casei_.element_layer}\n\n'
                    log.error(msg)
                    case_map = {}
                    keys_map2 = {}
                    cases2 = []
                    return case_map, keys_map2, cases2
                raise RuntimeError(msg)
            #  stack on [*nlayers*, eid_layer]
            element_layer2 = np.vstack(element_layers)
            # stack on [itime, *nlayers*, 10]
            data3 = np.hstack(datas)

            #isort1 = np.argsort(element_layer2[:,1])
            #isort2 = np.argsort(element_layer2[isort1,0])
            #indexs = np.arange(len(isort1))
            #isort = indexs[isort1][isort2]
            element_layer4, iisort = get_composite_sort(element_layer2)
            data4 = data3[:, iisort, :]
            assert data3.shape == data4.shape

            #assert np.array_equal(data2, data3), 'wut...'

            case.data = data4
            #case.element_id = np.hstack(element_ids)
            case.element_layer = element_layer4
            case.nelements = nelements_all
            case.ntotal = ntotal_all
            #case.data = data2
        cases2[key] = case

    #keys = keys_map.keys()
    #for key in keys:
        #key_cases =
        #for etype, case in case_map.items():
    return case_map, keys_map, cases2

def get_composite_sort(element_layers: list[np.ndarray]):
    uelement_layers = unique2d(element_layers)
    assert len(element_layers) == len(uelement_layers)

    iisort = np.lexsort((
        element_layers[:, 1],
        element_layers[:, 0],
    ))
    element_layer_stacked = element_layers[iisort]

    uelement_layers_stacked = unique2d(element_layer_stacked)
    assert len(element_layer_stacked) == len(uelement_layers_stacked)
    return element_layer_stacked, iisort

@nocrash_log
def get_composite_plate_stress_strains2(cases: CasesDict,
                                        eids: np.ndarray,
                                        model: OP2,
                                        times: np.ndarray,
                                        key: NastranKey,
                                        icase: int,
                                        form_dict, header_dict,
                                        keys_map: KeysMap,
                                        log: SimpleLogger,
                                        use_new_sidebar_objects: bool,
                                        is_stress: bool=True) -> int:
    if not use_new_sidebar_objects:  # pragma: no cover
        return icase

    case_map, keys_map2, cases2 = _stack_composite_results(model, log, is_stress, key=key)
    subcase_id = key[0]

    if not cases2:
        return icase

    for key, value in keys_map2.items():
        if key not in keys_map:
            keys_map[key] = value

    case = cases2[key]
    case_headers = case.get_headers()
    word, method_map = _composite_method_map(is_stress)
    methods = [method_map[headeri] for headeri in case_headers]

    # verify we don't crash when we try to pivot later
    # why does this happen???
    _datai = case.data[0, :, 0]
    _eids = case.element_layer[:, 0]
    _layer = case.element_layer[:, 1]
    unused_mytable, unused_myrows = pivot_table(_datai, _eids, _layer, shape=1)
    utable = unique2d(case.element_layer)
    assert np.array_equal(case.element_layer, utable)

    if len(case._times) != case.data.shape[0]:
        return icase

    titleii = f'Composite Plate {word}'
    res = CompositeStrainStressResults2(
        subcase_id, model, eids,
        case, method_map, titleii,
        #dim_max=1.0,
        data_format='%g',
        is_variable_data_format=False,
        nlabels=None, labelsize=None, ncolors=None, colormap='',
        set_max_min=False, uname='CompositeStressResults2')

    #form_layers = {'temp': [],}
    form_names = []
    ntimes_max = case.data.shape[0]
    for itime, dt in enumerate(times):
        if itime == ntimes_max:
            break
        #dt = case._times[itime]
        header = _get_nastran_header(case, dt, itime)
        header_dict[(key, itime)] = header

        formi = []
        form = form_dict[(key, itime)]
        form.append((f'Composite Plate {word}', None, formi))
        form_dict[(key, itime)] = form

        for imethod, method in enumerate(methods): # o11, o22, o12, maxprincipal, ...
            cases[icase] = (res, (subcase_id, (itime, imethod, header)))
            formi.append((f'{method}', icase, []))
            form_name2 = f'Composite Plate {word}: {method}'
            form_names.append(form_name2)
            icase += 1
    return icase

@nocrash_log
def get_composite_plate_stress_strains(cases: CasesDict,
                                       eids: np.ndarray,
                                       model: OP2,
                                       times: np.ndarray,
                                       key: NastranKey,
                                       icase: int,
                                       form_dict, header_dict,
                                       keys_map: KeysMap,
                                       composite_data_dict: CompositeDict,
                                       log: SimpleLogger,
                                       use_old_sidebar_objects: bool,
                                       is_stress: bool=True) -> int:
    """
    helper method for _fill_op2_time_centroidal_stress.
    Gets the stress/strain for each layer.

    element_type = 'CTRIA3'
    key = (1, 1, 1, 0, 0, '', '')
    composite_data_dict[element_type][key]
    """
    if not use_old_sidebar_objects:  # pragma: no cover
        return icase
    case_map, keys_map2, cases2 = _stack_composite_results(
        model, log, is_stress, key=key)
    if len(case_map) == 0:
        # we got into an error state
        return icase
    subcase_id = key[0]

    #assert len(cases) == icase
    #icase0 = icase

    composite_plates_ieids = []
    for element_type, composite_data in composite_data_dict.items():
        try:
            element_layer, ueids, data2, vm_word, ntimes, headers = composite_data[key]
        except KeyError:
            print(composite_data)
            raise

        #print(element_type, ueids)
        i = np.searchsorted(eids, ueids)
        if len(i) != len(np.unique(i)):
            log.error(f' duplicate eids for composite {element_type}...'
                      f'i={i} eids={eids} ueids={ueids}')
            continue
        composite_plates_ieids.append(i)
        #for itime2, header in enumerate(headers):
            #header_dict[(key, itime2)] = header
            #asdf

    #  no elements
    if not composite_plates_ieids:
        return icase

    try:
        case_dict = case_map[element_type]
    except KeyError:
        log.warning(f'skipping is_stress={is_stress} element_type={element_type}')
        raise
        #return icase
    case = case_dict[key]

    composite_plates_ieids = np.hstack(composite_plates_ieids)
    ieid_max = len(eids)
    #print('ieid_max =', ieid_max)

    case_headers = case.get_headers()
    #print('case_headers =', case_headers, vm_word)
    word, method_map = _composite_method_map(is_stress)
    methods = [method_map[headeri] for headeri in case_headers]

    #headersi = case.get_headers()
    #print('headersi =', headersi)
    #titles = []

    scalars_array = []
    for element_type, composite_data in composite_data_dict.items():
        unused_element_layer, unused_ueids, data2, unused_vm_word, unused_ntimes, unused_headers = composite_data[key]
        scalars_array.append(data2)
    if len(scalars_array) == 0:
        return icase

    try:
        scalars_array = concatenate_scalars(scalars_array)
    except ValueError:
        log.error('problem concatenating composite plates')
        return icase

    #print('scalars_array.shape =', scalars_array.shape)
    unused_ntimes, unused_nelements, nlayers, unused_nresults = scalars_array.shape

    #titles = []  # legend title
    headers = [] # sidebar word
    uname = f'Composite Plate {word}'
    res = LayeredTableResults(
        subcase_id, headers, composite_plates_ieids,
        ieid_max, scalars_array, methods,
        data_formats=None,
        colormap='jet', uname=uname)
    form_names = res.form_names

    #times = case._times
    for itime, dt in enumerate(times):
        #dt = case._times[itime]
        header = _get_nastran_header(case, dt, itime)
        header_dict[(key, itime)] = header

        formi = []
        form = form_dict[(key, itime)]
        form.append((f'Composite Plate {word}', None, formi))
        # formi = form[0][2]
        form_dict[(key, itime)] = form

        form_layers = {'temp': [],}
        for ilayer in range(nlayers):
            layer_name = f' Layer {ilayer+1}'
            form_layeri = []
            formi.append((layer_name.strip().lstrip(), None, form_layeri))
            #formi.append((layer_name.strip() + ' ' + str(case.element_name), None, form_layeri))
            form_layers[layer_name] = form_layeri

        for imethod, method in enumerate(methods): # o11, o22, o12, maxprincipal, ...
            for ilayer in range(nlayers):
                layer_name = f' Layer {ilayer+1}'
                form_layeri = form_layers[layer_name]
                #cases[icase] = (res, (subcase_id, header))
                #if use_new_sidebar_objects:
                    #cases[icase] = (res2, (subcase_id, (itime, ilayer, imethod, header)))
                    #form_layeri.append((f'{method} ({layer_name})', icase, []))
                    #form_name2 = f'{element_type} Composite Plate2 {word}: {method} ({layer_name})'
                    #form_names.append(form_name2)
                    #icase += 1

                cases[icase] = (res, (subcase_id, (itime, ilayer, imethod, header)))
                form_layeri.append((f'{method} ({layer_name})', icase, []))
                #form_name2 = f'Composite Plate {word}: {method} ({layer_name})'
                #form_names.append(form_name2)
                icase += 1
    #assert len(cases) == icase
    return icase

def _composite_method_map(is_stress: bool,
                          ) -> tuple[str, dict[str, str]]:
    if is_stress:
        word = 'Stress'
        method_map = {
            #'fiber_distance' : 'FiberDist.',
            'o11' : 'Stress 11',
            'o22' : 'Stress 22',
            't12' : 'Stress 12',
            't1z' : 'Stress 1Z',
            't2z' : 'Stress 2Z',
            'angle' : 'Theta',
            'major' : 'Max Principal',
            'minor' : 'Min Principal',
            'max_shear' : 'Max Shear',
            #'von_mises' : 'von Mises',
        }
    else:
        word = 'Strain'
        method_map = {
            #'fiber_distance' : 'FiberDist.',
            'e11' : 'Strain 11',
            'e22' : 'Strain 22',
            'e12' : 'Strain 12',
            'e1z' : 'Strain 1Z',
            'e2z' : 'Strain 2Z',
            'angle' : 'Theta',
            'major' : 'Max Principal',
            'minor' : 'Min Principal',
            'max_shear' : 'Max Shear',
            #'von_mises' : 'von Mises',
        }
    #methods = ['fiber_distance'] + [method_map[headeri] for headeri in case_headers]
    #methods = case_headers
    return word, method_map

@nocrash_log
def get_solid_stress_strains2(cases: CasesDict,
                              node_id: np.ndarray,
                              element_id: np.ndarray,
                              model: OP2,
                              times: np.ndarray,
                              key: NastranKey,
                              icase: int,
                              form_dict, header_dict,
                              keys_map: KeysMap,
                              log: SimpleLogger,
                              use_new_sidebar_objects: bool,
                              use_new_terms: bool,
                              is_stress: bool,
                              prefix: str='') -> int:
    """
    helper method for _fill_op2_time_centroidal_stress.
    Gets the max/min stress for each layer.

    key : varies
        (1, 5, 1, 0, 0, '', '')
        (isubcase, analysis_code, sort_method, count, ogs,
         superelement_adaptivity_index, pval_step)

    """
    if not use_new_sidebar_objects:  # pragma: no cover
        return icase
    solids, word, subcase_id, analysis_code = _get_solids(
        model, key, is_stress, prefix)

    solid_cases = []
    for solid_case in solids:
        if solid_case.is_complex:
            continue
        solid_eids = np.unique(solid_case.element_node[:, 0])
        common_eids = np.intersect1d(element_id, solid_eids)
        if len(common_eids) == 0:
            continue
        solid_cases.append(solid_case)

    if len(solid_cases) == 0:
        return icase

    #solid_case_headers = solid_case.get_headers()
    is_von_mises = solid_case.is_von_mises
    assert isinstance(is_von_mises, bool), is_von_mises
    von_mises = 9 if is_von_mises else 'von_mises'
    max_shear = 9 if not is_von_mises else 'max_shear'

    if is_stress:
        #['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz', 'omax', 'omid', 'omin', von_mises]
        iresult_to_title_annotation_map = {
            # iresult: (sidebar_label, annotation)
            0 : ('Stress XX', 'XX'),
            1 : ('Stress YY', 'YY'),
            2 : ('Stress ZZ', 'ZZ'),
            3 : ('Shear XY', 'XY'),
            4 : ('Shear YZ', 'YZ'),
            5 : ('Shear XZ', 'XZ'),

            6 : ('sMax Principal', 'Max Principal'),
            8 : ('sMin Principal', 'Min Principal'),
            7 : ('sMid Principal', 'Mid Principal'),
            #'abs_principal' : ('sAbs Principal', 'Abs Principal'),
            von_mises : ('von Mises', 'von Mises'), # the magnitude is large
            max_shear : ('Max Shear', 'Max Shear'), # the magnitude is large
        }
        word = 'Stress'
    else:
        iresult_to_title_annotation_map = {
            # iresult: (sidebar_label, annotation)
            0 : ('Strain XX', 'XX'),
            1 : ('Strain YY', 'YY'),
            2 : ('Strain ZZ', 'ZZ'),
            3 : ('Shear XY', 'XY'),
            4 : ('Shear YZ', 'YZ'),
            5 : ('Shear XZ', 'XZ'),

            6 : ('eMax Principal', 'Max Principal'),
            8 : ('eMin Principal', 'Min Principal'),
            7 : ('eMid Principal', 'Mid Principal'),
            von_mises : ('von Mises', 'von Mises'), # the magnitude is small
            max_shear : ('Max Shear', 'Max Shear'), # the magnitude is small
        }
        word = 'Strain'

    if not use_new_terms:
        del iresult_to_title_annotation_map[max_shear]
        #del iresult_to_title_annotation_map['abs_principal']

    title = f'Solid {word}'
    keys_map[key] = KeyMap(solid_case.subtitle, solid_case.label,
                           solid_case.superelement_adaptivity_index,
                           solid_case.pval_step)

    res = SolidStrainStressResults2(
        subcase_id, model,
        node_id, element_id,
        solid_cases, iresult_to_title_annotation_map, title,
        data_format='%g', is_variable_data_format=False,
        nlabels=None, labelsize=None, ncolors=None, colormap='',
        set_max_min=False,
        uname='SolidStressStrainResults2')

    times = solid_case._times
    for itime, dt in enumerate(times):
        #dt = case._times[itime]
        header = _get_nastran_header(solid_case, dt, itime)
        header2 = header.replace(' = ', '=')
        header_dict[(key, itime)] = header

        formi = []
        form = form_dict[(key, itime)]
        form.append(('Solid ' + word, None, formi))
        form_dict[(key, itime)] = form

        for iresult, (sidebar_label, annotation_label) in iresult_to_title_annotation_map.items():
            cases[icase] = (res, (subcase_id, (itime, iresult, header2)))
            formi.append((sidebar_label, icase, []))
            icase += 1
    return icase

def _get_solids(results: OP2,
                key,
                is_stress: bool,
                prefix: str) -> tuple[str, list, int, int]:
    analysis_code = key[1]
    #print("***stress eids=", eids)
    subcase_id = key[0]
    #if prefix == 'modal_contribution':
        #results = model.op2_results.modal_contribution
        #preword = 'Modal Contribution '
    #elif prefix == '':
        #results = model
        #preword = ''
    #else:  # pragma: no cover
        #raise NotImplementedError(prefix)

    if is_stress:
        cards = [
            results.ctetra_stress, results.cpenta_stress, results.chexa_stress, # results.cpyram_stress,
        ]
        word = 'Stress'
    else:
        cards = [
            results.ctetra_strain, results.cpenta_strain, results.chexa_strain, # results.cpyram_strain,
        ]
        word = 'Strain'

    cards2 = []
    for result in cards:
        if key not in result:
            continue
        cards2.append(result[key])
    return cards2, word, subcase_id, analysis_code

@nocrash_log
def get_solid_stress_strains(cases: CasesDict,
                             eids: np.ndarray,
                             model: OP2,
                             times: np.ndarray,
                             key: NastranKey,
                             icase: int,
                             form_dict, header_dict,
                             keys_map: KeysMap,
                             log: SimpleLogger,
                             use_old_sidebar_objects: bool,
                             is_stress: bool) -> int:
    """
    helper method for _fill_op2_time_centroidal_stress.
    """
    if not use_old_sidebar_objects:
        return icase
    #print("***stress eids=", eids)
    subcase_id = key[0]
    results = model
    if is_stress:
        solids = [
            results.ctetra_stress, results.cpenta_stress, results.chexa_stress, # results.cpyram_stress,
        ]
        word = 'Stress (centroid)'
    else:
        solids = [
            results.ctetra_strain, results.cpenta_strain, results.chexa_strain, # results.cpyram_strain,
        ]
        word = 'Strain (centroid)'

    #titles = []
    solid_cases = []
    solid_ieids = []
    for result in solids:
        if key not in result:
            continue
        case = result[key]
        if isinstance(case, RealSolidArrayNx):
            log.info(f'converting {case.class_name}')
            case = case.to_real_solid_array()
            result[key] = case

        #print(case)
        nnodes = case.nnodes_per_element
        all_eidsi = case.element_node[:, 0]
        nall_eidsi = len(all_eidsi)
        nelementsi = nall_eidsi // nnodes
        eidsi = all_eidsi.reshape(nelementsi, nnodes)[:, 0]

        i = np.searchsorted(eids, eidsi)
        if len(i) != len(np.unique(i)):
            #print(case.element_node)
            #print('element_name=%s nnodes_per_element=%s' % (case.element_name, nnodes_per_element))
            #print('iplate = %s' % i)
            #print('  eids = %s' % eids)
            #print('  eidsiA = %s' % case.element_node[:, 0])
            #print('  eidsiB = %s' % eidsi)
            msg = 'i%s (solid)=%s is not unique' % (case.element_name, str(i))
            #msg = 'iplate=%s is not unique' % str(i)
            log.warning(msg)
            continue
        if i.max() == len(eids):
            log.error('skipping because lookup is out of range...')
            continue
        #print('i =', i, i.max())
        #print('eids =', eids, len(eids))
        #print('eidsi =', eidsi, eids)
        #print(f'------------adding i={i} for {case.element_name}-----------')
        solid_element_id = np.unique(case.element_node[:, 0])
        common_eids = np.intersect1d(eids, solid_element_id)
        if len(common_eids) == 0:
            continue

        solid_cases.append(case)
        solid_ieids.append(i)
    if not solid_ieids:
        return icase

    solid_ieids = np.hstack(solid_ieids)
    ieid_max = len(eids)
    #print('ieid_max =', ieid_max)

    case = solid_cases[0]
    case_headers = case.get_headers()
    if is_stress:
        #sigma = 'σ'
        method_map = {
            'oxx' : 'Stress XX',
            'oyy' : 'Stress YY',
            'ozz' : 'Stress ZZ',
            'txy' : 'Shear XY',
            'tyz' : 'Shear YZ',
            'txz' : 'Shear XZ',

            'omax' : 'Max Principal',
            'omin' : 'Min Principal',
            'omid' : 'Mid Principal',
            'von_mises' : 'von Mises',
            'max_shear' : 'Max Shear',
        }
        data_format = '%.3f'
    else:
        #sigma = 'ϵ'
        method_map = {
            'exx' : 'Strain XX',
            'eyy' : 'Strain YY',
            'ezz' : 'Strain ZZ',
            'exy' : 'Shear XY',
            'eyz' : 'Shear YZ',
            'exz' : 'Shear XZ',

            'emax' : 'Max Principal',
            'emin' : 'Min Principal',
            'emid' : 'Mid Principal',
            'von_mises' : 'von Mises',
            'max_shear' : 'Max Shear',
        }
        data_format = '%.3e'
    methods = [method_map[headeri] for headeri in case_headers]
    #if 'Mises' in methods:
        #methods.append('Max shear')
    #else:
        #methods.append(f'{sigma} von Mises')

    #if case.is_von_mises:
        #vm_word = 'vonMises'
    #else:
        #vm_word = 'maxShear'

    #headersi = case.get_headers()
    #print('headersi =', headersi)

    scalars_array = []
    for case in solid_cases:
        #if case.is_complex:
            #log.warning(f'skipping complex Rod {word}')
            #continue

        #ntimes, nelements, nresults = case.data.shape
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
        #                                         txy, angle,
        #                                         majorP, minorP, ovm]

        keys_map[key] = KeyMap(case.subtitle, case.label,
                               case.superelement_adaptivity_index,
                               case.pval_step)

        #nnodes_per_element = case.nnodes
        #nelements_nnodes = nnodes_nlayers // 2
        #nelements = nelements_nnodes // nnodes_per_element
        #nlayers = 2
        nnodes = case.nnodes_per_element
        scalars = case.data
        #print('scalars.shape', scalars.shape)
        ntimes, nall, nresults = scalars.shape
        nelements = nall // nnodes

        if isinstance(case, (RealSolidArray, ComplexSolidArray)):
            scalars_save = scalars.reshape(ntimes, nelements, nnodes, nresults)
            scalars_array.append(scalars_save[:, :, 0, :])  # centroidal stress
        else:
            raise NotImplementedError(case.class_name)

    if len(scalars_array) == 0:
        return icase

    scalars_array = concatenate_scalars(scalars_array)

    #titles = []  # legend title
    headers = [] # sidebar word
    assert scalars_array.shape[2] == len(methods), f'shape={scalars_array.shape}; methods={methods} n={len(methods)}'
    uname = f'Solid {word}'
    res = SimpleTableResults(
        subcase_id, headers, solid_ieids, ieid_max, scalars_array, methods,
        data_format=data_format,
        colormap='jet', uname=uname)
    icase = add_simple_methods_to_form(icase, cases, key, subcase_id, word, res, case,
                                       form_dict, header_dict, methods,
                                       name='Solid')
    return icase

@nocrash_log
def get_spring_stress_strains(cases: CasesDict,
                              eids: np.ndarray,
                              model: OP2,
                              times: np.ndarray,
                              key: NastranKey,
                              icase: int,
                              form_dict, header_dict,
                              keys_map: KeysMap,
                              log: SimpleLogger,
                              is_stress: bool) -> int:
    """
    helper method for _fill_op2_time_centroidal_stress.
    """
    subcase_id = key[0]
    results = model
    if is_stress:
        springs = [
            results.celas1_stress, results.celas2_stress,
            results.celas3_stress, results.celas4_stress]
        word = 'Stress'
    else:
        springs = [
            results.celas1_strain, results.celas2_strain,
            results.celas3_strain, results.celas4_strain]
        word = 'Strain'

    spring_cases = []
    spring_ieids = []
    for result in springs:
        if key not in result:
            continue
        case = result[key]

        eidsi = case.element

        i = np.searchsorted(eids, eidsi)
        if len(i) != len(np.unique(i)):
            #print(case.element_node)
            #print('element_name=%s nnodes_per_element=%s' % (case.element_name, nnodes_per_element))
            #print('iplate = %s' % i)
            #print('  eids = %s' % eids)
            #print('  eidsiA = %s' % case.element_node[:, 0])
            #print('  eidsiB = %s' % eidsi)
            msg = 'i%s (spring)=%s is not unique' % (case.element_name, str(i))
            #msg = 'iplate=%s is not unique' % str(i)
            model.log.warning(msg)
            continue
        if i.max() == len(eids):
            model.log.error('skipping because lookup is out of range...')
            continue
        spring_cases.append(case)
        spring_ieids.append(i)
    if not spring_ieids:
        return icase

    spring_ieids = np.hstack(spring_ieids)
    ieid_max = len(eids)

    case = spring_cases[0]
    case_headers = case.get_headers()
    if is_stress:
        method_map = {
            'spring_stress' : 'Stress XX',
        }
        data_format = '%.3f'
    else:
        method_map = {
            'spring_strain' : 'Strain XX',
        }
        data_format = '%.3e'
    methods = [method_map[headeri] for headeri in case_headers]

    #headersi = case.get_headers()
    #print('headersi =', headersi)

    scalars_array = []
    for case in spring_cases:
        keys_map[key] = KeyMap(case.subtitle, case.label,
                               case.superelement_adaptivity_index,
                               case.pval_step)
        scalars = case.data
        scalars_array.append(scalars)

    if len(scalars_array) == 0:
        return icase

    scalars_array = concatenate_scalars(scalars_array)

    headers = [] # sidebar word
    res = SimpleTableResults(
        subcase_id, headers, spring_ieids, ieid_max, scalars_array, methods,
        data_format=data_format,
        colormap='jet', uname='Spring ' + word)

    icase = add_simple_methods_to_form(icase, cases, key, subcase_id, word, res, case,
                                       form_dict, header_dict, methods,
                                       name='Spring')
    return icase

def add_simple_methods_to_form(icase: int,
                               cases: CasesDict,
                               key: NastranKey,
                               subcase_id: int,
                               word: str,
                               res, case,
                               form_dict, header_dict,
                               methods, name: str) -> int:
    times = case._times
    nmethods = len(methods)
    if nmethods == 1:
        imethod = 0
        method = methods[imethod]
        for itime, dt in enumerate(times):
            header = _get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header

            form = form_dict[(key, itime)]
            form.append((name + ' ' + word, icase, []))
            cases[icase] = (res, (subcase_id, (itime, imethod, header)))
            icase += 1
    else:
        for itime, dt in enumerate(times):
            header = _get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header

            formi = []
            form = form_dict[(key, itime)]
            form.append((name + ' ' + word, None, formi))
            # formi = form[0][2]

            for imethod, method in enumerate(methods):
                #cases[icase] = (res, (subcase_id, header))
                cases[icase] = (res, (subcase_id, (itime, imethod, header)))
                formi.append((method, icase, []))
                icase += 1
    return icase

def concatenate_scalars(scalars_array: list[np.ndarray]) -> np.ndarray:
    if len(scalars_array) == 1:
        scalars_array = scalars_array[0]
    else:
        #print(len(scalars_array))
        #for arrayi in scalars_array:
            #print('   ', arrayi.shape)
        scalars_array = np.concatenate(scalars_array, axis=1)
    return scalars_array
