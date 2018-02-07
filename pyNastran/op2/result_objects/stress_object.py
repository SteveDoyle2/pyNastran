from __future__ import print_function
from collections import OrderedDict
from six import iteritems, iterkeys
import numpy as np

#vm_word = get_plate_stress_strain(
    #model, key, is_stress, vm_word, itime,
    #oxx, oyy, txy, max_principal, min_principal, ovm, is_element_on,
    #header_dict, keys_map)

class StressObject(object):
    def __init__(self, model, key, all_eids, is_stress=True):
        #print('--StressObject--')
        self.model = model
        self.vm_word = None
        self.header_dict = OrderedDict()
        self.keys_map = {}
        self.composite_ieids = {}
        self.is_stress = is_stress

        self.composite_data_dict = create_composite_plates(model, key, is_stress, self.keys_map)
        #self.plates_data_dict = create_plates(model, key, is_stress)

        #for key in iterkeys(self.plates_data_dict):
            #(case.element_node, ueids, data2, vm_word, ntimes) = self.plates_data_dict[key]
            #(min_data, max_data) = data2

        for element_type, composite_data in iteritems(self.composite_data_dict):
            for key in composite_data:
                element_layer, ueids, data2, vm_word, ntimes, headers = composite_data[key]
                #all_eids = ueids

                ntimes, neids, nlayers, nresults = data2.shape
                #data2[itime, eid, layer, oxx/oyy/...]
                ieids = np.searchsorted(all_eids, ueids)
                self.composite_ieids[element_type] = ieids
        #if len(self.composite_data_dict) > 0:
            #print(str(self))
            #print('-------------********------------')


    def set_composite_stress_old(self,
                                 key, itime, oxx, oyy, txy, tyz, txz,
                                 max_principal, min_principal, ovm,
                                 is_element_on, header_dict):
        #print("setting...")

        for element_type, composite_data in iteritems(self.composite_data_dict):
            try:
                element_layer, ueids, data2, vm_word, ntimes, headers = composite_data[key]
            except KeyError:
                print(composite_data)
                raise
            #print("vm_word = ", key, vm_word)
            for itime2, header in enumerate(headers):
                header_dict[(key, itime2)] = header

            ntimes, neids, nlayers, nresults = data2.shape
            #data2[itime, eid, layer, oxx/oyy/...]

            ieids = self.composite_ieids[element_type]
            #print('setting...')
            data3 = data2[itime, :, :, :]
            oxx[ieids] = np.nanmax(data3[:, :, 0], axis=1)
            oyy[ieids] = np.nanmax(data3[:, :, 1], axis=1)
            txy[ieids] = np.nanmax(data3[:, :, 2], axis=1)
            txz[ieids] = np.nanmax(data3[:, :, 3], axis=1)
            tyz[ieids] = np.nanmax(data3[:, :, 4], axis=1)
            # angle
            max_principal[ieids] = np.nanmax(data3[:, :, 6], axis=1)
            min_principal[ieids] = np.nanmin(data3[:, :, 7], axis=1)
            ovm[ieids] = np.nanmax(data3[:, :, 8], axis=1)
            if itime == 0:
                is_element_on[ieids] = 1
            #assert oxxi.shape == (ntimes, neids)
        return vm_word

    def set_composite_stress(self,
                             oxx, oyy, txy, tyz, txz,
                             max_principal, min_principal, ovm,
                             ):  # pragma: no cover

        for element_type, composite_data in iteritems(self.composite_data_dict):
            element_layer, ueids, data2, vm_word, ntimes, headers = composite_data[key]

            ntimes, neids, nlayers, nresults = data2.shape
            #data2[itime, eid, layer, oxx/oyy/...]

            ieids = self.composite_ieids[element_type]
            oxx[ieids] = np.nanmax(data2[:, :, :, 0], axis=2)
            oyy[ieids] = np.nanmax(data2[:, :, :, 1], axis=2)
            txy[ieids] = np.nanmax(data2[:, :, :, 2], axis=2)
            txz[ieids] = np.nanmax(data2[:, :, :, 3], axis=2)
            tyz[ieids] = np.nanmax(data2[:, :, :, 4], axis=2)
            # angle
            max_principal[ieids] = np.nanmax(data2[:, :, :, 6], axis=2)
            min_principal[ieids] = np.nanmin(data2[:, :, :, 7], axis=2)
            ovm[ieids] = np.nanmax(data2[:, :, :, 8], axis=2)
            assert oxxi.shape == (ntimes, neids)

    def __repr__(self):
        if self.is_stress:
            msg = 'StressObject:\n'
        else:
            msg = 'StrainObject:\n'
        msg += '    composite_data_dict.keys() = %s\n' % str(list(self.composite_data_dict.keys()))
        msg += '    composite_ieids = %s\n' % self.composite_ieids
        return msg

def create_plates(model, key, is_stress):
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        plates = [
            model.ctria3_stress, model.cquad4_stress,
            model.ctria6_stress, model.cquad8_stress,
            model.ctriar_stress, model.cquadr_stress,
        ]
    else:
        plates = [
            model.ctria3_strain, model.cquad4_strain,
            model.ctria6_strain, model.cquad8_strain,
            model.ctriar_strain, model.cquadr_strain,
        ]

    #[o11, o22, t12, t1z, t2z, angle, major, minor, max_shear]
    isotropic_data_dict = {}
    for obj_dict in plates:
        cases_to_delete = []
        for case_key, case in iteritems(obj_dict):
            if case_key == key:
                #data_dict[element_type] = [ueids, data]
                case_to_delete = case_key
                cases_to_delete.append(case_to_delete)

                if case.is_von_mises:
                    vm_word = 'vonMises'
                else:
                    vm_word = 'maxShear'

                nnodes_per_element = case.nnodes
                nlayers_per_element = nnodes_per_element * 2  # *2 for every other layer
                #eidsi = case.element_node[::nlayers_per_element, 0]  # ::2 is for layer skipping
                ueids = np.unique(case.element_node[:, 0])
                neids = len(ueids)

                #if 0:
                    #i = np.searchsorted(eids, eidsi)
                    #if len(i) != len(np.unique(i)):
                        #print(case.element_node)
                        #print('element_name=%s nnodes_per_element=%s' % (
                            #case.element_name, nnodes_per_element))
                        #print('iplate = %s' % i)
                        #print('eids = %s' % eids)
                        #print('eidsiA = %s' % case.element_node[:, 0])
                        #print('eidsiB = %s' % eidsi)
                        #msg = 'iplate=%s is not unique' % str(i)
                        #raise RuntimeError(msg)

                ntotal = case.data.shape[1]  # (ndt, ntotal, nresults)
                if nlayers_per_element == 1:
                    j = None
                else:
                    j = np.arange(ntotal)[::nlayers_per_element]

                #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
                #                                         txy, angle,
                #                                         majorP, minorP, ovm]

                #print("nlayers_per_element = ", case.element_name, nlayers_per_element)
                max_vals = [1, 2, 3, 5, 7]
                ntimes = case.data.shape[0]
                min_start = case.data[:, j, 6]
                max_start = case.data[:, j, :][:, :, max_vals]
                min_data_list = [min_start]
                max_data_list = [max_start]
                for inode in range(1, nlayers_per_element):
                    datai = case.data[:, j+inode, :]
                    min_data_list.append(case.data[:, j+inode, :][:, :, max_vals])
                    max_data_list.append(case.data[:, j+inode, :])
                if len(min_data_list) == 1:
                    min_data = min_data_list[0]
                    max_data = max_data_list[0]
                else:
                    #print(np.dstack(max_data_list).shape)
                    # (ntimes, neids, nresults)
                    # (1, 4, 13) = np.dstack(max_data_list)
                    min_data = np.amin(np.dstack(min_data_list), axis=0)
                    max_data = np.amax(np.dstack(max_data_list), axis=0)
                #assert min_data.shape == (ntimes, neids, 1), min_data.shape
                #assert max_data.shape == (ntimes, neids, 5), max_data.shape
                del min_data_list, max_data_list

                #eidsi = case.element_node[:, 0]
                #layers = case.element_node[:, 1]
                data2 = (min_data, max_data)
                isotropic_data_dict[key] = [case.element_node, ueids, data2, vm_word, ntimes]
        #for case_key in cases_to_delete:
            #del obj_dict[case_key]
    return isotropic_data_dict

def create_composite_plates(model, key, is_stress, keys_map):
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        cplates = [
            ('CTRIA3', model.ctria3_composite_stress),
            ('CQUAD4', model.cquad4_composite_stress),
            ('CTRIA6', model.ctria6_composite_stress),
            ('CQUAD8', model.cquad8_composite_stress),
            #model.ctriar_composite_stress,
            #model.cquadr_composite_stress,
        ]
    else:
        cplates = [
            ('CTRIA3', model.ctria3_composite_strain),
            ('CQUAD4', model.cquad4_composite_strain),
            ('CTRIA6', model.ctria6_composite_strain),
            ('CQUAD8', model.cquad8_composite_strain),
            #model.ctriar_composite_strain,
            #model.cquadr_composite_strain,
        ]

    #[o11, o22, t12, t1z, t2z, angle, major, minor, max_shear]
    composite_data_dict = {}
    for element_type, obj_dict in cplates:
        cases_to_delete = []
        if len(obj_dict) == 0:
            continue

        composite_data_dict[element_type] = {}
        for case_key, case in iteritems(obj_dict):
            if case_key == key:
                keys_map[key] = (case.subtitle, case.label, case.superelement_adaptivity_index)
                #data_dict[element_type] = [ueids, data]
                case_to_delete = case_key
                cases_to_delete.append(case_to_delete)

                if case.is_von_mises:
                    vm_word = 'vonMises'
                else:
                    vm_word = 'maxShear'
                eidsi = case.element_layer[:, 0]
                layers = case.element_layer[:, 1]
                ntimes = case.data.shape[0]
                data2, ueids = pivot(case.data, eidsi, layers)

                headers = []
                for itime, dt in enumerate(case._times):
                    header = _get_nastran_header(case, dt, itime)
                    headers.append(header)
                composite_data_dict[element_type][key] = [
                    case.element_layer, ueids, data2, vm_word, ntimes, headers]
        if len(composite_data_dict[element_type]) == 0:
            del composite_data_dict[element_type]
        #for case_key in cases_to_delete:
            #del obj_dict[case_key]


    return composite_data_dict


def pivot(data, rows, cols):
    """
    PCOMP: rows=element_ids, cols=layer
    """
    ncount = len(rows)
    icount = np.arange(ncount)
    assert len(data.shape) == 3
    ntimes = data.shape[0]
    nresults = data.shape[-1]

    rows_new, row_pos_new = np.unique(rows, return_inverse=True)
    cols_new, col_pos_new = np.unique(cols, return_inverse=True)
    nrows = len(rows_new)
    ncols = len(cols_new)

    pivot_table = np.full((nrows, ncols), -1, dtype='int32')
    pivot_table[row_pos_new, col_pos_new] = icount
    #print(pivot_table)

    ipivot_row, ipivot_col = np.where(pivot_table != -1)
    data2 = np.full((ntimes, nrows, ncols, nresults), np.nan, dtype=data.dtype)
    data2[:, ipivot_row, ipivot_col, :] = data[:, icount, :]

    return data2, rows_new

def _get_nastran_header(case, dt, itime):
    #if case is None:
        #return None
    try:
        code_name = case.data_code['name']
    except KeyError:
        return 'Static'

    if isinstance(dt, float):
        header = ' %s = %.4E' % (code_name, dt)
    else:
        header = ' %s = %i' % (code_name, dt)

    # cases:
    #   1. lsftsfqs
    #   2. loadIDs, eigrs
    #   3. lsdvmns, eigrs
    #   ???
    if hasattr(case, 'mode_cycle'):
        header += '; freq = %g Hz' % case.mode_cycles[itime]
    elif hasattr(case, 'cycles'):
        header += '; freq = %g Hz' % case.cycles[itime]
    elif hasattr(case, 'eigis'):
        eigi = case.eigis[itime]
        cycle = np.abs(eigi) / (2. * np.pi)
        header += '; freq = %g Hz' % cycle
    elif hasattr(case, 'eigns'):  #  eign is not eigr; it's more like eigi
        eigi = case.eigrs[itime] #  but |eigi| = sqrt(|eign|)
        cycle = np.sqrt(np.abs(eigi)) / (2. * np.pi)
        header += '; freq = %g Hz' % cycle
    elif hasattr(case, 'dt'):
        time = case._times[itime]
        header += '; time = %g sec' % time
    elif hasattr(case, 'lftsfqs') or hasattr(case, 'lsdvmns') or hasattr(case, 'loadIDs'):
        pass
        #raise RuntimeError(header)
    else:
        msg = 'unhandled case; header=%r\n%s' % (header, str(case))
        print(msg)
        #raise RuntimeError(msg)

    return header.strip('; ')

