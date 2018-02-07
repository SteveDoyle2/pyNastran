from __future__ import print_function
from collections import OrderedDict
from six import iteritems
import numpy as np

#vm_word = get_plate_stress_strain(
    #model, key, is_stress, vm_word, itime,
    #oxx, oyy, txy, max_principal, min_principal, ovm, is_element_on,
    #header_dict, keys_map)

class StressObject(object):
    def __init__(self, model, key, is_stress=True):
        self.model = model
        self.vm_word = None
        self.header_dict = OrderedDict()
        self.keys_map = {}

        #vm_word = self.get_composite_plate_stress_strain(
            #model, key, is_stress, vm_word, itime,
            #oxx, oyy, txy, tyz, txz,
            #max_principal, min_principal, ovm, is_element_on,
            #self.header_dict, self.keys_map, eid_map)

        self.composite_data_dict = create_composite_plates(model, key, is_stress)
        self.plates_data_dict = create_plates(model, key, is_stress)
        all_eids = np.arange(100)
        ueids = np.array([1, 2, 10])
        ieids = np.searchsorted(all_eids, ueids)


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
                case_to_delete = [case_key]
                cases_to_delete.append(case_to_delete)

                if case.is_von_mises:
                    vm_word = 'vonMises'
                else:
                    vm_word = 'maxShear'

                nnodes_per_element = case.nnodes
                nlayers_per_element = nnodes_per_element * 2  # *2 for every other layer
                #eidsi = case.element_node[::nlayers_per_element, 0]  # ::2 is for layer skipping
                ueids = np.unique(case.element_node[:, 0])

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
                min_data_list = [case.data[:, j, 6]]
                max_data_list = [case.data[:, j, [1, 2, 3, 5, 7]]]
                for inode in range(1, nlayers_per_element):
                    min_data_list.append(case.data[:, j+inode, 6])
                    max_data_list.append(case.data[:, j+inode, [1, 2, 3, 5, 7]])
                if len(min_data_list) == 1:
                    min_data = min_data_list[0]
                    max_data = max_data_list[0]
                else:
                    min_data = np.amin(np.dstack(min_data_list), axis=0)
                    max_data = np.amax(np.dstack(max_data_list), axis=0)
                del min_data_list, max_data_list

                #eidsi = case.element_node[:, 0]
                #layers = case.element_node[:, 1]
                data2 = (min_data, max_data)
                isotropic_data_dict[key] = [case.element_node, ueids, data2, vm_word]
        for case_key in cases_to_delete:
            del obj_dict[case_key]
    return isotropic_data_dict

def create_composite_plates(model, key, is_stress):
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
    for unused_element_type, obj_dict in iteritems(cplates):
        cases_to_delete = []
        for case_key, case in iteritems(obj_dict):
            if case_key == key:
                #data_dict[element_type] = [ueids, data]
                case_to_delete = [case_key]
                cases_to_delete.append(case_to_delete)

                if case.is_von_mises:
                    vm_word = 'vonMises'
                else:
                    vm_word = 'maxShear'
                eidsi = case.element_layer[:, 0]
                layers = case.element_layer[:, 1]
                data2, ueids = pivot(case.data, eidsi, layers)
                composite_data_dict[key] = [case.element_layer, ueids, data2, vm_word]

        for case_key in cases_to_delete:
            del obj_dict[case_key]
    return composite_data_dict


def pivot(data, rows, cols):
    """
    PCOMP: rows=element_ids, cols=layer
    """
    ncount = len(rows)
    icount = np.arange(ncount)
    assert len(data.shape) == 2
    unused_ntimes = data.shape[0]
    nresults = data.shape[-1]

    rows_new, row_pos_new = np.unique(rows, return_inverse=True)
    cols_new, col_pos_new = np.unique(cols, return_inverse=True)
    nrows = len(rows_new)
    ncols = len(cols_new)

    pivot_table = np.full((nrows, ncols), -1, dtype='int32')
    pivot_table[row_pos_new, col_pos_new] = icount
    #print(pivot_table)

    ipivot_row, ipivot_col = np.where(pivot_table != -1)
    data2 = np.full((nrows, ncols, nresults), np.nan, dtype=data.dtype)
    data2[ipivot_row, ipivot_col, :] = data[icount, :]

    return data2, rows_new
