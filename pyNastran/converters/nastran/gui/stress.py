# encoding: utf-8
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from pyNastran.converters.nastran.gui.results import SimpleTableResults, LayeredTableResults
from pyNastran.op2.result_objects.stress_object import _get_nastran_header

if TYPE_CHECKING: # pragma: no cover
    from pyNastran.op2.op2 import OP2


def get_rod_stress_strains(eids, cases, model: OP2, times, key, icase,
                           form_dict, header_dict, keys_map, is_stress):
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
            model.log.warning(msg)
            continue
        if i.max() == len(eids):
            model.log.error('skipping because lookup is out of range...')
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
            'axial' : 'σxx',
            'torsion' : 'τxy',
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
            'axial' : 'ϵxx',
            'torsion' : 'ϵxy',
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
        keys_map[key] = (case.subtitle, case.label,
                         case.superelement_adaptivity_index, case.pval_step)
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

def get_bar_stress_strains(eids, cases, model: OP2, times, key, icase,
                           form_dict, header_dict, keys_map, is_stress):
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
            model.log.warning(msg)
            continue
        if i.max() == len(eids):
            model.log.error('skipping because lookup is out of range...')
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

    case = bar_cases[0]
    case_headers = case.get_headers()
    #print(case_headers)
    if is_stress:
        #sigma = 'σ'
        method_map = {
             's1a' : 'σxx 1A',
             's2a' : 'σxx 2A',
             's3a' : 'σxx 3A',
             's4a' : 'σxx 4A',

             's1b' : 'σxx 1B',
             's2b' : 'σxx 2B',
             's3b' : 'σxx 3B',
             's4b' : 'σxx 4B',

            'axial' : 'σxx',
            'MS_tension' : 'MS_tension',
            'MS_compression' : 'MS_compression',
            'smaxa' : 'σmax A',
            'smina' : 'σmin A',
            'smaxb' : 'σmax B',
            'sminb' : 'σmin B',
            #'von_mises' : 'σ von Mises',
        }
        data_format = '%.3f'
    else:
        #sigma = 'ϵ'
        method_map = {
            'e1a' : 'ϵxx 1A',
            'e2a' : 'ϵxx 2A',
            'e3a' : 'ϵxx 3A',
            'e4a' : 'ϵxx 4A',

            'e1b' : 'ϵxx 1B',
            'e2b' : 'ϵxx 2B',
            'e3b' : 'ϵxx 3B',
            'e4b' : 'ϵxx 4B',


            'e1c': 'ϵxx 1C',
            'e1d': 'ϵxx 1D',
            'e2c': 'ϵxx 2C',
            'e2d': 'ϵxx 2D',

           'axial' : 'ϵxx',
           'MS_tension' : 'MS_tension',
           'MS_compression' : 'MS_compression',
           'emaxa' : 'ϵmax A',
           'emina' : 'ϵmin A',
           'emaxb' : 'ϵmax B',
           'eminb' : 'ϵmin B',
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

        keys_map[key] = (case.subtitle, case.label,
                         case.superelement_adaptivity_index, case.pval_step)

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

def get_beam_stress_strains(eids, cases, model: OP2, times, key, icase,
                            form_dict, header_dict, keys_map, is_stress):
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
            model.log.warning(msg)
            continue
        if i.max() == len(eids):
            model.log.error('skipping because lookup is out of range...')
            continue

        i2 = np.hstack([i, i + 1]).T.flatten()
        print('i =', i)
        print('i2 =', i2)
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
    print('ieid_max =', ieid_max)

    case = beam_cases[0]
    case_headers = case.get_headers()
    #print(case_headers)
    if is_stress:
        #sigma = 'σ'
        method_map = {
            'sxc' : 'σxx C',
            'sxd' : 'σxx D',
            'sxe' : 'σxx E',
            'sxf' : 'σxx F',
            #'torsion' : 'τxy',
            'MS_tension' : 'MS_tension',
            'MS_compression' : 'MS_compression',
            'smax' : 'σmax',
            'smin' : 'σmin',
            #'von_mises' : 'σ von Mises',
        }
    else:
        #sigma = 'ϵ'
        method_map = {
            'sxc' : 'ϵxx C',
            'sxd' : 'ϵxx D',
            'sxe' : 'ϵxx E',
            'sxf' : 'ϵxx F',

            'exc' : 'ϵxx C',
            'exd' : 'ϵxx D',
            'exe' : 'ϵxx E',
            'exf' : 'ϵxx F',
            #'torsion' : 'τxy',
            'MS_tension' : 'MS_tension',
            'MS_compression' : 'MS_compression',
            'smax' : 'ϵmax',
            'smin' : 'ϵmin',
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
            model.log.warning(f'skipping complex Beam {word}')
            continue

        #ntimes, nelements, nresults = case.data.shape
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
        #                                         txy, angle,
        #                                         majorP, minorP, ovm]

        keys_map[key] = (case.subtitle, case.label,
                         case.superelement_adaptivity_index, case.pval_step)

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

def get_plate_stress_strains(eids, cases, model: OP2, times, key, icase,
                             form_dict, header_dict, keys_map, is_stress,
                             prefix=''):
    """
    helper method for _fill_op2_time_centroidal_stress.
    Gets the max/min stress for each layer.
    """
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
        #print(case)

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
            model.log.warning(msg)
            continue
        if i.max() == len(eids):
            model.log.error('skipping because lookup is out of range...')
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
        #sigma = 'σ'
        method_map = {
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
        #sigma = 'ϵ'
        method_map = {
            'fiber_curvature' : 'FiberCurvature',
            'fiber_distance' : 'FiberDistance',
            'exx' : 'ϵxx',
            'eyy' : 'ϵyy',
            'exy' : 'ϵxy',
            'angle' : 'θ',
            'emax' : 'ϵmax',
            'emin' : 'ϵmin',
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

        keys_map[key] = (case.subtitle, case.label,
                         case.superelement_adaptivity_index, case.pval_step)

        nnodes_per_element = case.nnodes_per_element
        nelements_nnodes = nnodes_nlayers // 2
        nelements = nelements_nnodes // nnodes_per_element
        nlayers = 2
        scalars = case.data.reshape(ntimes, nelements, nnodes_per_element, nlayers, nresults)
        scalars_array.append(scalars[:, :, 0, :, :])

    if len(scalars_array) == 0:
        return icase

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
    for itime, dt in enumerate(times):
        #dt = case._times[itime]
        header = _get_nastran_header(case, dt, itime)
        header_dict[(key, itime)] = header

        formi = []
        form = form_dict[(key, itime)]
        form.append(('Plate ' + word, None, formi))
        # formi = form[0][2]
        form_dict[(key, itime)] = form

        for ilayer in range(2):
            layer = layer_names[ilayer]
            form_layeri = []
            formi.append((layer, None, form_layeri))
            for imethod, method in enumerate(methods):
                #cases[icase] = (res, (subcase_id, header))
                cases[icase] = (res, (subcase_id, (itime, ilayer, imethod, header)))
                form_layeri.append((f'{method} ({layer})', icase, []))
                icase += 1
    return icase

def get_composite_plate_stress_strains(eids, cases, model: OP2, times, key, icase,
                                       form_dict, header_dict, keys_map,
                                       composite_data_dict, log, is_stress=True):
    """
    helper method for _fill_op2_time_centroidal_stress.
    Gets the stress/strain for each layer.
    """
    subcase_id = key[0]

    plates_ieids = []
    for element_type, composite_data in composite_data_dict.items():
        try:
            element_layer, ueids, data2, vm_word, ntimes, headers = composite_data[key]
        except KeyError:
            print(composite_data)
            raise


        #print(element_type, ueids)
        i = np.searchsorted(eids, ueids)
        if len(i) != len(np.unique(i)):
            model.log.error(f' duplicate eids for composite {element_type}...'
                            f'i={i} eids={eids} ueids={ueids}')
            continue
        plates_ieids.append(i)
        #for itime2, header in enumerate(headers):
            #header_dict[(key, itime2)] = header
            #asdf

    if not plates_ieids:
        return icase

    case_map = {
        # is_stress, element_name
        (True, 'CTRIA3') : model.ctria3_composite_stress,
        (False, 'CTRIA3') : model.ctria3_composite_strain,

        (True, 'CQUAD4') : model.cquad4_composite_stress,
        (False, 'CQUAD4') : model.cquad4_composite_strain,

        (True, 'CTRIAR') : model.ctriar_composite_stress,
        (False, 'CTRIAR') : model.ctriar_composite_strain,
        (True, 'CQUADR') : model.cquadr_composite_stress,
        (False, 'CQUADR') : model.cquadr_composite_strain,
    }
    try:
        case_dict = case_map[(is_stress, element_type)]
    except KeyError:
        log.warning(f'skipping is_stress={is_stress} element_type={element_type}')
        return icase
    case = case_dict[key]


    plates_ieids = np.hstack(plates_ieids)
    ieid_max = len(eids)
    #print('ieid_max =', ieid_max)

    case_headers = case.get_headers()
    #print('case_headers =', case_headers, vm_word)
    if is_stress:
        word = 'Stress'
        #sigma = 'σ'
        method_map = {
            'o11' : 'σ11',
            'o22' : 'σ22',
            't12' : 't12',
            't1z' : 'τ1z',
            't2z' : 'τ2z',
            'angle' : 'θ',
            'major' : 'σ major',
            'minor' : 'σ minor',
            'max_shear' : 'MaxShear',
            #'von_mises' : 'σ von Mises',
        }
    else:
        word = 'Strain'
        #sigma = 'ϵ'
        method_map = {
            'e11' : 'ϵ11',
            'e22' : 'ϵ22',
            'e12' : 'ϵ12',
            'e1z' : 'ϵ1z',
            'e2z' : 'ϵ2z',
            'angle' : 'θ',
            'major' : 'ϵ major',
            'minor' : 'ϵ minor',
            'max_shear' : 'MaxShear',
        }
    methods = [method_map[headeri] for headeri in case_headers]
    #methods = case_headers

    #if 'Mises' in methods:
        #methods.append('Max shear')
    #else:
        #methods.append(f'{sigma} von Mises')

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
    res = LayeredTableResults(
        subcase_id, headers, plates_ieids, ieid_max, scalars_array, methods,
        data_formats=None,
        colormap='jet', uname='Composite Plate ' + word)

    #times = case._times
    for itime, dt in enumerate(times):
        #dt = case._times[itime]
        header = _get_nastran_header(case, dt, itime)
        header_dict[(key, itime)] = header

        formi = []
        form = form_dict[(key, itime)]
        form.append(('Composite Plate ' + word, None, formi))
        # formi = form[0][2]
        form_dict[(key, itime)] = form

        form_layers = {}
        for ilayer in range(nlayers):
            layer_name = f' Layer {ilayer+1}'
            form_layeri = []
            formi.append((layer_name, None, form_layeri))
            form_layers[layer_name] = form_layeri

        for imethod, method in enumerate(methods):
            for ilayer in range(nlayers):
                layer_name = f' Layer {ilayer+1}'
                form_layeri = form_layers[layer_name]
                #cases[icase] = (res, (subcase_id, header))
                cases[icase] = (res, (subcase_id, (itime, ilayer, imethod, header)))
                form_layeri.append((f'{method} ({layer_name})', icase, []))
                icase += 1
    return icase

def get_solid_stress_strains(eids, cases, model: OP2, times, key, icase,
                             form_dict, header_dict, keys_map, is_stress):
    """
    helper method for _fill_op2_time_centroidal_stress.
    """
    #print("***stress eids=", eids)
    subcase_id = key[0]
    results = model
    if is_stress:
        solids = [
            results.ctetra_stress, results.cpenta_stress, results.chexa_stress, # results.cpyram_stress,
        ]
        word = 'Stress'
    else:
        solids = [
            results.ctetra_strain, results.cpenta_strain, results.chexa_strain, # results.cpyram_strain,
        ]
        word = 'Strain'

    #titles = []
    solid_cases = []
    solid_ieids = []
    for result in solids:
        if key not in result:
            continue
        case = result[key]

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
            model.log.warning(msg)
            continue
        if i.max() == len(eids):
            model.log.error('skipping because lookup is out of range...')
            continue
        #print('i =', i, i.max())
        #print('eids =', eids, len(eids))
        #print('eidsi =', eidsi, eids)
        #print(f'------------adding i={i} for {case.element_name}-----------')
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
            'oxx' : 'σxx',
            'oyy' : 'σyy',
            'ozz' : 'σzz',
            'txy' : 'τxy',
            'tyz' : 'τyz',
            'txz' : 'τxz',

            'omax' : 'σmax',
            'omin' : 'σmin',
            'omid' : 'σmid',
            'von_mises' : 'σ von Mises',
        }
        data_format = '%.3f'
    else:
        #sigma = 'ϵ'
        method_map = {
            'exx' : 'ϵxx',
            'eyy' : 'ϵyy',
            'ezz' : 'ϵzz',
            'exy' : 'ϵxy',
            'eyz' : 'ϵyz',
            'exz' : 'ϵxz',

            'emax' : 'ϵmax',
            'emin' : 'ϵmin',
            'emid' : 'ϵmid',
            'von_mises' : 'ϵ von Mises',
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
            #model.log.warning(f'skipping complex Rod {word}')
            #continue

        #ntimes, nelements, nresults = case.data.shape
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
        #                                         txy, angle,
        #                                         majorP, minorP, ovm]

        keys_map[key] = (case.subtitle, case.label,
                         case.superelement_adaptivity_index, case.pval_step)

        #nnodes_per_element = case.nnodes
        #nelements_nnodes = nnodes_nlayers // 2
        #nelements = nelements_nnodes // nnodes_per_element
        #nlayers = 2
        nnodes = case.nnodes_per_element
        scalars = case.data
        #print('scalars.shape', scalars.shape)
        ntimes, nall, nresults = scalars.shape
        nelements = nall // nnodes
        scalars_save = scalars.reshape(ntimes, nelements, nnodes, nresults)
        scalars_array.append(scalars_save[:, :, 0, :])

    if len(scalars_array) == 0:
        return icase

    scalars_array = concatenate_scalars(scalars_array)

    #titles = []  # legend title
    headers = [] # sidebar word
    res = SimpleTableResults(
        subcase_id, headers, solid_ieids, ieid_max, scalars_array, methods,
        data_format=data_format,
        colormap='jet', uname='Solid ' + word)

    icase = add_simple_methods_to_form(icase, cases, key, subcase_id, word, res, case,
                                       form_dict, header_dict, methods,
                                       name='Solid')
    return icase

def get_spring_stress_strains(eids, cases, model: OP2, times, key, icase,
                              form_dict, header_dict, keys_map, is_stress):
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
            'spring_stress' : 'σxx',
        }
        data_format = '%.3f'
    else:
        method_map = {
            'spring_strain' : 'ϵxx',
        }
        data_format = '%.3e'
    methods = [method_map[headeri] for headeri in case_headers]

    #headersi = case.get_headers()
    #print('headersi =', headersi)

    scalars_array = []
    for case in spring_cases:
        keys_map[key] = (case.subtitle, case.label,
                         case.superelement_adaptivity_index, case.pval_step)
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

def add_simple_methods_to_form(icase, cases, key, subcase_id, word, res, case,
                               form_dict, header_dict,
                               methods, name):
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

def concatenate_scalars(scalars_array):
    if len(scalars_array) == 1:
        scalars_array = scalars_array[0]
    else:
        scalars_array = np.concatenate(scalars_array, axis=1)
    return scalars_array
