from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from pyNastran.converters.nastran.gui.results import SimpleTableResults # , LayeredTableResults
from .stress import add_simple_methods_to_form, concatenate_scalars

if TYPE_CHECKING: # pragma: no cover
    from pyNastran.op2.op2 import OP2


def get_spring_force(eids, cases, model: OP2, times, key, icase,
                     form_dict, header_dict, keys_map):
    """
    helper method for _fill_op2_time_centroidal_stress.
    """
    subcase_id = key[0]
    results = model

    springs = [
        results.celas1_force, results.celas2_force,
        results.celas3_force, results.celas4_force]
    word = 'Force'

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
    method_map = {
        'spring_force' : 'Fspring',
    }
    data_format = '%.3f'
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

def get_bar_force(eids, cases, model: OP2, times, key, icase,
                  form_dict, header_dict, keys_map):
    """
    helper method for _fill_op2_time_centroidal_stress.
    """
    #print("***stress eids=", eids)
    subcase_id = key[0]
    bars = [model.cbar_force]

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

    word = 'Force'
    method_map = {
        # real
        'bending_moment_a1' : 'My Bending 1',
        'bending_moment_a2' : 'My Bending 2',
        'bending_moment_b1' : 'Mz Bending 1',
        'bending_moment_b2' : 'Mz Bending 2',

        # complex
        'bending_moment_1a' : 'My Bending 1',
        'bending_moment_2a' : 'My Bending 2',
        'bending_moment_1b' : 'Mz Bending 1',
        'bending_moment_2b' : 'Mz Bending 2',
        'shear1' : 'Fy Shear',
        'shear2' : 'Fz Shear',
        'axial' : 'Fx',
        'torque' : 'Torque',
    }
    data_format = '%.3f'
    methods = [method_map[headeri] for headeri in case_headers]

    scalars_array = []
    for case in bar_cases:
        #if case.is_complex:
            #model.log.warning(f'skipping complex Bar {word}')
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
        scalars = case.data
        scalars_array.append(scalars)

    if len(scalars_array) == 0:
        return icase

    scalars_array = concatenate_scalars(scalars_array)

    #include_tension_margin = True
    #include_compression_margin = True

    headers = [] # sidebar word
    res = SimpleTableResults(
        subcase_id, headers, bar_ieids, ieid_max, scalars_array, methods,
        data_format=data_format,
        colormap='jet', uname='Bar ' + word)

    icase = add_simple_methods_to_form(icase, cases, key, subcase_id, word, res, case,
                                       form_dict, header_dict, methods,
                                       name='Bar')
    return icase


def get_plate_force(eids, cases, model: OP2, times, key, icase,
                    form_dict, header_dict, keys_map):
    """
    helper method for _fill_op2_time_centroidal_stress.
    """
    #print("***stress eids=", eids)
    subcase_id = key[0]
    plates = [
        model.ctria3_force, model.ctria6_force, model.ctriar_force,
        model.cquad4_force, model.cquad8_force, model.cquadr_force,
    ]

    plate_cases = []
    plate_ieids = []
    for result in plates:
        if key not in result:
            continue
        case = result[key]

        nnodes = case.nnodes_per_element
        if nnodes > 1:
            # element_node.shape = (20, 2)
            # element type: QUAD144
            #RealPlateBilinearForceArray
            # 8=[mx, my, mxy, bmx, bmy, bmxy, tx, ty]

            all_eids = case.element_node[:, 0]
            nids = case.element_node[:, 1]
            if nids[-1] != 0:
                eidsi = all_eids[::nnodes]
                #print(nids[::nnodes])
                #print('nnodes =', nnodes)
        else:
            eidsi = case.element

        #print(case)
        #print(eids)
        #print(eidsi)
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
        #model.log.info('saving i%s (plate)' % (case.element_name))
        #print(f'------------adding i={i} for {case.element_name}-----------')
        plate_cases.append(case)
        plate_ieids.append(i)
    if not plate_ieids:
        return icase

    plate_ieids = np.hstack(plate_ieids)
    ieid_max = len(eids)
    #print('ieid_max =', ieid_max)

    case = plate_cases[0]
    case_headers = case.get_headers()
    #print(case_headers)

    word = 'Force'
    method_map = {
        # real
        'mx': 'Fx',
        'my': 'Fy',
        'mxy': 'Fxy',
        'bmx': 'Mx',
        'bmy': 'My',
        'bmxy': 'Mxy',
        'tx': 'Vx',
        'ty': 'Vy',

        # complex
        #'bending_moment_1a' : 'M bending 1A',
        #'bending_moment_2a' : 'M bending 2A',
        #'bending_moment_1b' : 'M bending 1B',
        #'bending_moment_2b' : 'M bending 2B',
        #'shear1' : 'F shear 1',
        #'shear2' : 'F shear 2',
        #'axial' : 'F axial',
        #'torque' : 'torque',
        #'s1a' : 'σxx 1A',
        #'e1a' : 'ϵxx 1A',
    }
    data_format = '%.3f'
    methods = [method_map[headeri] for headeri in case_headers]

    scalars_array = []
    for case in plate_cases:
        keys_map[key] = (case.subtitle, case.label,
                         case.superelement_adaptivity_index, case.pval_step)

        nnodes = case.nnodes_per_element
        if nnodes == 1:
            scalars = case.data
        else:
            scalars = case.data[:, ::nnodes, :]
        scalars_array.append(scalars)

    if len(scalars_array) == 0:
        return icase

    scalars_array = concatenate_scalars(scalars_array)

    headers = [] # sidebar word
    res = SimpleTableResults(
        subcase_id, headers, plate_ieids, ieid_max, scalars_array, methods,
        data_format=data_format,
        colormap='jet', uname='Plate ' + word)

    icase = add_simple_methods_to_form(icase, cases, key, subcase_id, word, res, case,
                                       form_dict, header_dict, methods,
                                       name='Plate')
    return icase
