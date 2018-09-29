from __future__ import print_function
import sys
from collections import OrderedDict
import numpy as np

#vm_word = get_plate_stress_strain(
    #model, key, is_stress, vm_word, itime,
    #oxx, oyy, txy, max_principal, min_principal, ovm, is_element_on,
    #header_dict, keys_map)

IS_TESTING = 'test' in sys.argv[0]
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
        if IS_TESTING:
            self.plates_data_dict = create_plates(model, key, is_stress)

        #for key in self.plates_data_dict.keys():
            #(case.element_node, ueids, data2, vm_word, ntimes) = self.plates_data_dict[key]
            #(min_data, max_data) = data2

        for element_type, composite_data in self.composite_data_dict.items():
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

        for element_type, composite_data in self.composite_data_dict.items():
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

        for element_type, composite_data in self.composite_data_dict.items():
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
        for case_key, case in obj_dict.items():
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
        for case_key, case in obj_dict.items():
            if case_key == key:
                keys_map[key] = (case.subtitle, case.label,
                                 case.superelement_adaptivity_index, case.pval_step)
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


def get_rod_stress_strain(model, key, is_stress, vm_word, itime,
                          oxx, txy,
                          max_principal, min_principal, ovm, is_element_on,
                          eids, header_dict, keys_map):
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        rods = [model.crod_stress, model.conrod_stress, model.ctube_stress,]
    else:
        rods = [model.crod_strain, model.conrod_strain, model.ctube_strain,]

    for result in rods:
        if key not in result:
            continue

        case = result[key]
        if case.is_complex:
            continue
        eidsi = case.element
        i = np.searchsorted(eids, eidsi)
        if len(i) != len(np.unique(i)):
            msg = 'irod=%s is not unique\n' % str(i)
            print('eids = %s\n' % str(list(eids)))
            print('eidsi = %s\n' % str(list(eidsi)))
            raise RuntimeError(msg)

        is_element_on[i] = 1
        dt = case._times[itime]
        header = _get_nastran_header(case, dt, itime)
        header_dict[(key, itime)] = header
        keys_map[key] = (case.subtitle, case.label,
                         case.superelement_adaptivity_index, case.pval_step)

        # data=[1, nnodes, 4] where 4=[axial, SMa, torsion, SMt]
        oxx[i] = case.data[itime, :, 0]
        txy[i] = case.data[itime, :, 2]
        try:
            ovm[i] = np.sqrt(oxx[i]**2 + 3*txy[i]**2) # plane stress
        except FloatingPointError:
            ovm[i] = 0.
            assert np.allclose(oxx[i], 0.)
            assert np.allclose(txy[i], 0.)
        # max_principal[i] = sqrt(oxx[i]**2 + txy[i]**2)
        # min_principal[i] = max_principal[i] - 2 * txy[i]
        # simplification of:
        #   eig(A) = [oxx, txy]
        #            [txy, 0.0]
        # per Equation 7: http://www.soest.hawaii.edu/martel/Courses/GG303/Eigenvectors.pdf
        try:
            max_principal[i] = (oxx[i] + np.sqrt(oxx[i]**2 + 4 * txy[i]**2)) / 2.
            min_principal[i] = (oxx[i] - np.sqrt(oxx[i]**2 + 4 * txy[i]**2)) / 2.
        except FloatingPointError:
            # underflow is a thing...we can promote the dtype from float32 to float64
            # but we'll hold off until there's a real example
            assert np.allclose(oxx[i], 0.)
            assert np.allclose(txy[i], 0.)
            max_principal[i] = 0.
            min_principal[i] = 0.
    del rods
    return vm_word

def get_bar_stress_strain(model, key, is_stress, vm_word, itime,
                          oxx,
                          max_principal, min_principal, ovm, is_element_on,
                          eids, header_dict, keys_map):
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        bars = model.cbar_stress
    else:
        bars = model.cbar_strain

    if key in bars:
        case = bars[key]
        if case.is_complex:
            pass
        else:
            dt = case._times[itime]
            header = _get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header
            keys_map[key] = (case.subtitle, case.label,
                             case.superelement_adaptivity_index, case.pval_step)
            #s1a = case.data[itime, :, 0]
            #s2a = case.data[itime, :, 1]
            #s3a = case.data[itime, :, 2]
            #s4a = case.data[itime, :, 3]

            axial = case.data[itime, :, 4]
            smaxa = case.data[itime, :, 5]
            smina = case.data[itime, :, 6]
            #MSt = case.data[itime, :, 7]

            #s1b = case.data[itime, :, 8]
            #s2b = case.data[itime, :, 9]
            #s3b = case.data[itime, :, 10]
            #s4b = case.data[itime, :, 11]

            smaxb = case.data[itime, :, 12]
            sminb = case.data[itime, :, 13]
            #MSc   = case.data[itime, :, 14]

            eidsi = case.element # [:, 0]

            i = np.searchsorted(eids, eidsi)
            if len(i) != len(np.unique(i)):
                print('ibar = %s' % i)
                print('eids = %s' % eids)
                msg = 'ibar=%s is not unique' % str(i)
                raise RuntimeError(msg)

            is_element_on[i] = 1.
            oxx[i] = axial

            ## TODO :not sure if this block is general for multiple CBAR elements
            samax = np.amax([smaxa, smaxb], axis=0)
            samin = np.amin([smaxa, smaxb], axis=0)
            assert len(samax) == len(i), len(samax)
            assert len(samin) == len(i)
            savm = np.amax(np.abs(
                [smina, sminb,
                 smaxa, smaxb, axial]), axis=0)

            max_principal[i] = samax
            min_principal[i] = samin
            ovm[i] = savm
            del axial, smaxa, smina, smaxb, sminb, eidsi, i, samax, samin, savm
    del bars
    return vm_word

def try_except_return3(func):
    def try_except_func(*args, **kwargs):
        try:
            out = func(*args, **kwargs)
        except(MemoryError, NameError):
            raise
        except(RuntimeError, IndexError):
            if 'test_' in sys.argv[0]:
                raise
            out = args[3]
        return out
    return try_except_func

@try_except_return3
def get_bar100_stress_strain(model, key, is_stress, vm_word, itime,
                             oxx,
                             max_principal, min_principal, ovm, is_element_on,
                             eids, header_dict, keys_map):
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        bars2 = model.cbar_stress_10nodes
    else:
        bars2 = model.cbar_strain_10nodes

    if key in bars2:
        case = bars2[key]
        if case.is_complex:
            pass
        else:
            dt = case._times[itime]
            header = _get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header
            keys_map[key] = (case.subtitle, case.label,
                             case.superelement_adaptivity_index, case.pval_step)

            #  0    1    2    3    4     5     6     7     8
            # [sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]

            eidsi = case.element # [:, 0]
            ueidsi = np.unique(eidsi)
            istart = np.searchsorted(eidsi, ueidsi)
            unused_iend = np.hstack([istart[1:], [len(eidsi)]])
            axial = case.data[itime, :, 5]

            nbars = len(eidsi) // 10
            if nbars * 10 != len(eidsi):
                msg = 'nbars=%s neids=%s; expected neids=10*nbars=%s' % (nbars, len(eidsi), 10*nbars)
                raise RuntimeError(msg)

            axial = case.data[itime, :, 5].reshape(nbars, 10).min(axis=1)
            smax = case.data[itime, :, 6].reshape(nbars, 10).max(axis=1)
            smin = case.data[itime, :, 7].reshape(nbars, 10).min(axis=1)


            i = np.searchsorted(eids, eidsi)
            if len(i) != len(np.unique(i)):
                print('ibar = %s' % i)
                print('eids = %s' % eids)
                msg = 'ibar=%s is not unique' % str(i)
                raise RuntimeError(msg)

            is_element_on[i] = 1.
            oxx[i] = axial

            ## TODO :not sure if this block is general for multiple CBAR elements
            svm = np.amax(np.abs([smin, smin]), axis=0)

            max_principal[i] = smax
            min_principal[i] = smin
            ovm[i] = svm
            #del axial, smaxa, smina, smaxb, sminb, eidsi, i, samax, samin, savm
    del bars2
    return vm_word

def get_beam_stress_strain(model, key, is_stress, vm_word, itime,
                           oxx,
                           max_principal, min_principal, ovm, is_element_on,
                           header_dict, keys_map, eid_map):
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        beams = model.cbeam_stress
    else:
        beams = model.cbeam_strain

    if key in beams:
        case = beams[key]
        if case.is_complex:
            pass
        else:
            eidsi = case.element_node[:, 0]
            ueids = np.unique(eidsi)
            #neids = len(ueids)

            # sxc, sxd, sxe, sxf
            # smax, smin, MSt, MSc
            dt = case._times[itime]
            header = _get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header
            keys_map[key] = (case.subtitle, case.label,
                             case.superelement_adaptivity_index, case.pval_step)
            sxc = case.data[itime, :, 0]
            sxd = case.data[itime, :, 1]
            sxe = case.data[itime, :, 2]
            sxf = case.data[itime, :, 3]
            smax = case.data[itime, :, 4]
            smin = case.data[itime, :, 5]

            imin = np.searchsorted(eidsi, ueids)
            imax = np.searchsorted(eidsi, ueids, side='right')
            #sxxi = smax[imin:imax]
            for eid, imini, imaxi in zip(ueids, imin, imax):
                oxxi = 0.
                smaxi = 0.
                smini = 0.
                eid2 = eid_map[eid]
                is_element_on[eid2] = 1.
                oxxi = max(
                    sxc[imini:imaxi].max(),
                    sxd[imini:imaxi].max(),
                    sxe[imini:imaxi].max(),
                    sxf[imini:imaxi].max(),
                )
                smaxi = smax[imini:imaxi].max()
                smini = smin[imini:imaxi].min()
                ovmi = max(np.abs(smaxi), np.abs(smini))
                oxxi = oxx[eid2]
                max_principal[eid2] = smaxi
                min_principal[eid2] = smini
                ovm[eid2] = ovmi
            del eidsi, ueids, sxc, sxd, sxe, sxf, smax, smin, oxxi, smaxi, smini, ovmi
    del beams
    return vm_word

def get_plate_stress_strain(model, key, is_stress, vm_word, itime,
                            oxx, oyy, txy, max_principal, min_principal, ovm, is_element_on,
                            eids, header_dict, keys_map):
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

    for result in plates:
        if key not in result:
            continue
        case = result[key]
        if case.is_complex:
            continue

        if case.is_von_mises:
            vm_word = 'vonMises'
        else:
            vm_word = 'maxShear'

        nnodes_per_element = case.nnodes
        nlayers_per_element = nnodes_per_element * 2  # *2 for every other layer
        eidsi = case.element_node[::nlayers_per_element, 0]  # ::2 is for layer skipping

        i = np.searchsorted(eids, eidsi)
        if len(i) != len(np.unique(i)):
            print(case.element_node)
            print('element_name=%s nnodes_per_element=%s' % (case.element_name, nnodes_per_element))
            print('iplate = %s' % i)
            print('eids = %s' % eids)
            print('eidsiA = %s' % case.element_node[:, 0])
            print('eidsiB = %s' % eidsi)
            msg = 'iplate=%s is not unique' % str(i)
            raise RuntimeError(msg)
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
        #                                         txy, angle,
        #                                         majorP, minorP, ovm]
        try:
            is_element_on[i] = 1
        except IndexError:
            print(case.element_node)
            print('i = %s' % i)
            print('eids = %s' % eids)
            print('eidsi = %s' % eidsi)
            raise

        ntotal = case.data.shape[1]  # (ndt, ntotal, nresults)
        if nlayers_per_element == 1:
            j = None
        else:
            j = np.arange(ntotal)[::nlayers_per_element]

        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
        #                                         txy, angle,
        #                                         majorP, minorP, ovm]
        dt = case._times[itime]
        header = _get_nastran_header(case, dt, itime)
        header_dict[(key, itime)] = header
        keys_map[key] = (case.subtitle, case.label,
                         case.superelement_adaptivity_index, case.pval_step)
        oxxi = case.data[itime, j, 1]
        oyyi = case.data[itime, j, 2]
        txyi = case.data[itime, j, 3]
        o1i = case.data[itime, j, 5]
        o3i = case.data[itime, j, 6]
        ovmi = case.data[itime, j, 7]

        #print("nlayers_per_element = ", case.element_name, nlayers_per_element)
        for inode in range(1, nlayers_per_element):
            #print('%s - ilayer = %s' % (case.element_name, inode))
            #print(case.data[itime, j + inode, 1])
            #print(case.data[itime, :, 1])
            oxxi = np.amax(np.vstack([oxxi, case.data[itime, j + inode, 1]]), axis=0)
            oyyi = np.amax(np.vstack([oyyi, case.data[itime, j + inode, 2]]), axis=0)
            txyi = np.amax(np.vstack([txyi, case.data[itime, j + inode, 3]]), axis=0)
            o1i = np.amax(np.vstack([o1i, case.data[itime, j + inode, 5]]), axis=0)
            o3i = np.amin(np.vstack([o3i, case.data[itime, j + inode, 6]]), axis=0)
            ovmi = np.amax(np.vstack([ovmi, case.data[itime, j + inode, 7]]), axis=0)
            assert len(oxxi) == len(j)
            #print('-------')

        oxx[i] = oxxi
        oyy[i] = oyyi
        txy[i] = txyi
        max_principal[i] = o1i
        min_principal[i] = o3i
        ovm[i] = ovmi
    return vm_word


def get_solid_stress_strain(model, key, is_stress, vm_word, itime,
                            oxx, oyy, ozz, txy, tyz, txz,
                            max_principal, mid_principal, min_principal, ovm, is_element_on,
                            eids, header_dict, keys_map):
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        solids = [(model.ctetra_stress),
                  (model.cpenta_stress),
                  (model.chexa_stress),]
    else:
        solids = [(model.ctetra_strain),
                  (model.cpenta_strain),
                  (model.chexa_strain),]

    for result in solids:
        if key not in result:
            continue
        case = result[key]
        if case.is_complex:
            continue

        if case.is_von_mises:
            vm_word = 'vonMises'
        else:
            vm_word = 'maxShear'

        nnodes_per_element = case.nnodes
        eidsi = case.element_cid[:, 0]
        ntotal = len(eidsi)  * nnodes_per_element

        i = np.searchsorted(eids, eidsi)
        if len(i) != len(np.unique(i)):
            print('isolid = %s' % str(i))
            print('eids = %s' % eids)
            print('eidsi = %s' % eidsi)
            assert len(i) == len(np.unique(i)), 'isolid=%s is not unique' % str(i)

        is_element_on[i] = 1
        #self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz,
        #                                         txy, tyz, txz,
        #                                         o1, o2, o3, ovm]

        if nnodes_per_element == 1:
            j = None
        else:
            j = np.arange(ntotal)[::nnodes_per_element]
            ueidsi = np.unique(eidsi)
            assert len(j) == len(ueidsi), 'j=%s ueidsi=%s' % (j, ueidsi)

        dt = case._times[itime]
        header = _get_nastran_header(case, dt, itime)
        header_dict[(key, itime)] = header
        keys_map[key] = (case.subtitle, case.label,
                         case.superelement_adaptivity_index, case.pval_step)
        oxxi = case.data[itime, j, 0]
        oyyi = case.data[itime, j, 1]
        ozzi = case.data[itime, j, 2]
        txyi = case.data[itime, j, 3]
        tyzi = case.data[itime, j, 4]
        txzi = case.data[itime, j, 5]
        o1i = case.data[itime, j, 6]
        o2i = case.data[itime, j, 7]
        o3i = case.data[itime, j, 8]
        ovmi = case.data[itime, j, 9]

        for inode in range(1, nnodes_per_element):
            oxxi = np.amax(np.vstack([oxxi, case.data[itime, j + inode, 0]]), axis=0)
            oyyi = np.amax(np.vstack([oyyi, case.data[itime, j + inode, 1]]), axis=0)
            ozzi = np.amax(np.vstack([ozzi, case.data[itime, j + inode, 2]]), axis=0)
            txyi = np.amax(np.vstack([txyi, case.data[itime, j + inode, 3]]), axis=0)
            tyzi = np.amax(np.vstack([tyzi, case.data[itime, j + inode, 4]]), axis=0)
            txzi = np.amax(np.vstack([txzi, case.data[itime, j + inode, 2]]), axis=0)

            o1i = np.amax(np.vstack([o1i, case.data[itime, j + inode, 6]]), axis=0)
            o2i = np.amax(np.vstack([o2i, case.data[itime, j + inode, 7]]), axis=0)
            o3i = np.amin(np.vstack([o3i, case.data[itime, j + inode, 8]]), axis=0)
            ovmi = np.amax(np.vstack([ovmi, case.data[itime, j + inode, 9]]), axis=0)
            assert len(oxxi) == len(j)

        oxx[i] = oxxi
        oyy[i] = oyyi
        ozz[i] = ozzi
        txy[i] = txyi
        tyz[i] = tyzi
        txz[i] = txzi
        max_principal[i] = o1i
        mid_principal[i] = o2i
        min_principal[i] = o3i
        ovm[i] = ovmi
    del solids
    return vm_word
