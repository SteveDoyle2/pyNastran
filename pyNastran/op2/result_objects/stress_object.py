from __future__ import annotations
import sys
from typing import Any, TYPE_CHECKING

import numpy as np
from pyNastran.femutils.utils import pivot_table, abs_min_max

from pyNastran.utils.locale import func_str
from pyNastran.op2.result_objects.op2_objects import (
    real_modes_to_omega_freq, complex_damping_frequency)
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids_nx import RealSolidArrayNx


from pyNastran.gui.gui_objects.types import HeaderDict
from pyNastran.op2.op2_interface.types import KeyMap, KeysMap, NastranKey
if TYPE_CHECKING: # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.op2.op2 import OP2

float_types = (float, np.float32, np.float64)
int_types = (int, np.int32, np.int64)

CompositeResult = tuple[
    np.ndarray, np.ndarray, np.ndarray,
    str, int, str,
]
CompositeDict = dict[str, dict[Any, CompositeResult]]

#vm_word = get_plate_stress_strain(
    #model, key, is_stress, vm_word, itime,
    #oxx, oyy, txy, max_principal, min_principal, ovm, is_element_on,

    #header_dict, keys_map)

class StressObject:
    def __init__(self, model: OP2,
                 key: NastranKey,
                 all_eids: np.ndarray,
                 is_stress: bool=True) -> None:
        #print('--StressObject--')
        self.model = model
        self.vm_word = None
        self.header_dict: dict[Any, str] = {}
        self.keys_map: dict[NastranKey, str] = {}
        self.composite_ieids: dict[str, str] = {}
        self.is_stress = is_stress
        assert is_stress in {True, False}, is_stress

        self.composite_data_dict: CompositeDict = create_composite_plates(
            model, key, is_stress, self.keys_map)
        #self.plates_data_dict = create_plates(model, key, is_stress)

        #for key in self.plates_data_dict.keys():
            #(case.element_node, ueids, data2, vm_word, ntimes) = self.plates_data_dict[key]
            #(min_data, max_data) = data2

        for element_type, composite_data in self.composite_data_dict.items():
            for keyi in composite_data:
                element_layer, ueids, data2, vm_word, unused_ntimes, headers = composite_data[keyi]
                #all_eids = ueids

                ntimes, unused_neids, unused_nlayers, unused_nresults = data2.shape
                #data2[itime, eid, layer, oxx/oyy/...]
                ieids = np.searchsorted(all_eids, ueids)
                self.composite_ieids[element_type] = ieids
        #if len(self.composite_data_dict) > 0:
            #print(str(self))
            #print('-------------********------------')


    def set_composite_stress_old(self,
                                 key: int, itime: int,
                                 oxx: np.ndarray, oyy: np.ndarray,
                                 txy: np.ndarray, tyz: np.nditer, txz: np.ndarray,
                                 max_principal: np.ndarray, min_principal: np.ndarray,
                                 ovm: np.ndarray, is_element_on: np.ndarray,
                                 header_dict: HeaderDict) -> str:
        for element_type, composite_data in self.composite_data_dict.items():
            try:
                element_layer, unused_ueids, data2, vm_word, unused_ntimes, headers = composite_data[key]
            except KeyError:
                print(composite_data)
                raise
            for itime2, header in enumerate(headers):
                header_dict[(key, itime2)] = header

            ntimes, unused_neids, unused_nlayers, unused_nresults = data2.shape
            #data2[itime, eid, layer, oxx/oyy/...]

            ieids = self.composite_ieids[element_type]
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

    def set_composite_stress_by_layer(self,
                                      key: int, itime: int, nelements: int,
                                      header_dict: HeaderDict) -> str:
        """get all the ply stresses/strains"""
        nlayers = 0
        nelements2 = 0
        for element_type, composite_data in self.composite_data_dict.items():
            try:
                unused_element_layer, unused_ueids, data2, vm_word, unused_ntimes, headers = composite_data[key]
            except KeyError:
                print(composite_data)
                raise
            unused_ntimesi, neidsi, nlayersi, unused_nresultsi = data2.shape
            nlayers = max(nlayers, nlayersi)
            nelements2 += neidsi

        shape = (nelements2, nlayers)
        unused_element_ids = np.full(nelements2, np.nan, dtype='float32')
        oxx = np.full(shape, np.nan, dtype='float32')
        oyy = np.full(shape, np.nan, dtype='float32')

        txy = np.full(shape, np.nan, dtype='float32')
        tyz = np.full(shape, np.nan, dtype='float32')
        txz = np.full(shape, np.nan, dtype='float32')

        max_principal = np.full(shape, np.nan, dtype='float32')  # max
        min_principal = np.full(shape, np.nan, dtype='float32')  # min
        ovm = np.full(shape, np.nan, dtype='float32')

        for element_type, composite_data in self.composite_data_dict.items():
            unused_element_layer, unused_ueids, data2, vm_word, unused_ntimes, headers = composite_data[key]

            for itime2, header in enumerate(headers):
                header_dict[(key, itime2)] = header

            unused_ntimes, unused_neids, nlayers, unused_nresults = data2.shape
            #data2[itime, eid, layer, oxx/oyy/...]

            ieids = self.composite_ieids[element_type]
            data3 = data2[itime, :, :, :]

            # neids, nlayers, nresults
            oxx[ieids] = data3[:, :, 0]
            oyy[ieids] = data3[:, :, 1]
            txy[ieids] = data3[:, :, 2]
            txz[ieids] = data3[:, :, 3]
            tyz[ieids] = data3[:, :, 4]
            # angle
            max_principal[ieids] = data3[:, :, 6]
            min_principal[ieids] = data3[:, :, 7]
            ovm[ieids] = data3[:, :, 8]

            #assert oxxi.shape == (ntimes, neids)
        return vm_word

    def set_composite_stress(self,
                             key, oxx, oyy, txy, tyz, txz,
                             max_principal, min_principal, ovm,
                             ):  # pragma: no cover

        for element_type, composite_data in self.composite_data_dict.items():
            unused_element_layer, unused_ueids, data2, unused_vm_word, ntimes, unused_headers = composite_data[key]

            ntimes, neids, unused_nlayers, unused_nresults = data2.shape
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
            assert oxx.shape == (ntimes, neids)

    def __repr__(self) -> str:
        if self.is_stress:
            msg = 'StressObject:\n'
        else:
            msg = 'StrainObject:\n'
        msg += '    composite_data_dict.keys() = %s\n' % str(list(self.composite_data_dict.keys()))
        msg += '    composite_ieids = %s\n' % self.composite_ieids
        return msg

def create_plates(model: Any, key: NastranKey, is_stress: bool) -> dict[str, Any]:
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        stress = model.op2_results.stress
        plates = [
            stress.ctria3_stress, stress.cquad4_stress,
            stress.ctria6_stress, stress.cquad8_stress,
            stress.ctriar_stress, stress.cquadr_stress,
        ]
    else:
        strain = model.op2_results.strain
        plates = [
            strain.ctria3_strain, strain.cquad4_strain,
            strain.ctria6_strain, strain.cquad8_strain,
            strain.ctriar_strain, strain.cquadr_strain,
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
                #neids = len(ueids)

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
                    #datai = case.data[:, j+inode, :]
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

def create_composite_plates(model, key: NastranKey,
                            is_stress: bool,
                            keys_map: KeysMap) -> CompositeDict:
    """helper method for _fill_op2_time_centroidal_stress

    Returns
    -------
    composite_data_dict: dict[key, value]
         value: [o11, o22, t12, t1z, t2z, angle, major, minor, max_shear]

    """
    if is_stress:
        stress = model.op2_results.stress
        cplates_dict = {
            'CTRIA3' : stress.ctria3_composite_stress,
            'CTRIA6' : stress.ctria6_composite_stress,
            'CTRIAR' : stress.ctriar_composite_stress,
            'CQUAD4' : stress.cquad4_composite_stress,
            'CQUAD8' : stress.cquad8_composite_stress,
            'CQUADR' : stress.cquadr_composite_stress,
        }
    else:
        strain = model.op2_results.strain
        cplates_dict = {
            'CTRIA3' : strain.ctria3_composite_strain,
            'CTRIA6' : strain.ctria6_composite_strain,
            'CTRIAR' : strain.ctriar_composite_strain,
            'CQUAD4' : strain.cquad4_composite_strain,
            'CQUAD8' : strain.cquad8_composite_strain,
            'CQUADR' : strain.cquadr_composite_strain,
        }
    composite_data_dict = {}
    ncases = sum([len(res_dict) for res_dict in cplates_dict.values()])
    if ncases == 0:
        return composite_data_dict

    #[o11, o22, t12, t1z, t2z, angle, major, minor, max_shear]
    for element_type, res_dict in cplates_dict.items():
        cases_to_delete = []
        if len(res_dict) == 0:
            continue

        composite_data_dict[element_type] = {}
        for case_key, case in res_dict.items():
            if case_key != key:
                continue
            print(f'adding key={key} to keys_map')
            keys_map[key] = KeyMap(case.subtitle, case.label,
                                   case.superelement_adaptivity_index,
                                   case.pval_step)
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
            data2, ueids = pivot_table(case.data, eidsi, layers)

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


def _get_nastran_header(case: Any,
                        dt: int | float,
                        itime: int) -> str:
    #if case is None:
        #return None
    try:
        code_name = case.data_code['name']
    except KeyError:
        return 'Static'

    if code_name == 'mode' and isinstance(dt, float_types):
        dt = int(dt)
        #print(f'code_name={code_name}; dt={dt}; type={type(dt)}')

    if isinstance(dt, float_types):
        header = ' %s = %.4E' % (code_name, dt)
    elif isinstance(dt, int_types):
        header = ' %s = %i' % (code_name, dt)
    elif dt is None:
        return 'Static'
    else:  # pragma: no cover
        print(case, dt, type(dt))

    # cases:
    #   1. lsftsfqs
    #   2. loadIDs, eigrs
    #   3. lsdvmns, eigrs
    #   ???
    if hasattr(case, 'mode_cycle'):
        header += '; freq = %s Hz' % func_str(case.mode_cycles[itime])
    elif hasattr(case, 'cycles'):
        header += '; freq = %s Hz' % func_str(case.cycles[itime])
    elif hasattr(case, 'eigis'):
        eigi = np.array(case.eigis)
        eigr = np.array(case.eigrs)
        dampings, freqs = complex_damping_frequency(eigr, eigi)
        damping = dampings[itime]
        freq = freqs[itime]
        header += '; freq = %s Hz; g=%s' % (func_str(freq), func_str(damping))
    elif hasattr(case, 'eigns'):
        unused_omegas, freqs = real_modes_to_omega_freq(case.eigns)
        freq = freqs[itime]
        header += '; freq = %s Hz' % func_str(freq)
    elif hasattr(case, 'freqs'):
        header += '; freq = %s Hz' % func_str(case.freqs[itime])

    elif hasattr(case, 'times'):
        time = case._times[itime]
        header += '; time = %s sec' % func_str(time)
    elif hasattr(case, 'dts'):
        time = case._times[itime]
        header += '; time = %s sec' % func_str(time)
        #print(header)

    elif hasattr(case, 'lftsfqs'):
        step = case.lftsfqs[itime]
        header += '; step = %s' % func_str(step)
    elif hasattr(case, 'lsdvmns'):
        step = case.lsdvmns[itime]
        header += '; step = %s' % func_str(step)
    elif hasattr(case, 'loadIDs'):
        step = case.loadIDs[itime]
        header += '; step = %s' % func_str(step)
    elif hasattr(case, 'loadFactors'):
        step = case.loadFactors[itime]
        header += '; step = %s' % func_str(step)
    elif hasattr(case, 'load_steps'):
        step = case.load_steps[itime]
        header += '; step = %s' % func_str(step)
    else:
        msg = 'unhandled case; header=%r\n%s' % (header, str(case))
        print(msg)
        #raise RuntimeError(msg)

    return header.strip('; ')


def get_rod_stress_strain(model: OP2,
                          key: str,
                          is_stress: bool,
                          vm_word: str,
                          itime: int,
                          oxx: np.ndarray, txy: np.ndarray,
                          max_principal: np.ndarray, min_principal: np.ndarray,
                          ovm: np.ndarray, is_element_on: np.ndarray,
                          eids: np.ndarray, header_dict: HeaderDict,
                          keys_map: KeysMap,
                          log: SimpleLogger) -> str:
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        stress = model.op2_results.stress
        rods = [stress.crod_stress, stress.conrod_stress, stress.ctube_stress,]
    else:
        strain = model.op2_results.strain
        rods = [strain.crod_strain, strain.conrod_strain, strain.ctube_strain,]

    for result in rods:
        if key not in result:
            continue

        case = result[key]
        if case.is_complex:
            continue
        eidsi = case.element
        i = np.searchsorted(eids, eidsi)
        if len(i) != len(np.unique(i)):
            msg = 'i%s (rod)=%s is not unique' % (case.element_name, str(i))
            print('  eids = %s\n' % str(list(eids)))
            print('  eidsi = %s\n' % str(list(eidsi)))
            log.warning(msg)
            continue
            #raise RuntimeError(msg)

        try:
            is_element_on[i] = 1
        except IndexError:
            log.warning('missing %ss' % case.element_name)
            continue

        try:
            dt = case._times[itime]
        except IndexError:
            print(case)
            print(f'len(case._times)={len(case._times)} itime={itime}')
            raise
        header = _get_nastran_header(case, dt, itime)
        header_dict[(key, itime)] = header
        keys_map[key] = KeyMap(case.subtitle, case.label,
                               case.superelement_adaptivity_index,
                               case.pval_step)

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

def get_bar_stress_strain(model: OP2, key,
                          is_stress: bool, vm_word: str, itime: int,
                          oxx: np.ndarray,
                          max_principal: np.ndarray, min_principal: np.ndarray,
                          ovm: np.ndarray, is_element_on: np.ndarray,
                          eids: np.ndarray, header_dict: HeaderDict,
                          keys_map: KeysMap,
                          log: SimpleLogger) -> str:
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        stress = model.op2_results.stress
        bars = stress.cbar_stress
    else:
        strain = model.op2_results.strain
        bars = strain.cbar_strain

    if key not in bars:
        return vm_word

    case = bars[key]
    if case.is_complex:
        return vm_word

    dt = case._times[itime]
    header = _get_nastran_header(case, dt, itime)
    header_dict[(key, itime)] = header
    keys_map[key] = KeyMap(case.subtitle, case.label,
                           case.superelement_adaptivity_index,
                           case.pval_step)
    eidsi = case.element # [:, 0]
    common_eids = np.intersect1d(eids, eidsi)
    i = np.searchsorted(eids, common_eids)
    j = np.searchsorted(eidsi, common_eids)

    #s1a = case.data[itime, :, 0]
    #s2a = case.data[itime, :, 1]
    #s3a = case.data[itime, :, 2]
    #s4a = case.data[itime, :, 3]

    axial = case.data[itime, j, 4]
    smaxa = case.data[itime, j, 5]
    smina = case.data[itime, j, 6]
    #MSt = case.data[itime, :, 7]

    #s1b = case.data[itime, :, 8]
    #s2b = case.data[itime, :, 9]
    #s3b = case.data[itime, :, 10]
    #s4b = case.data[itime, :, 11]

    smaxb = case.data[itime, j, 12]
    sminb = case.data[itime, j, 13]
    #MSc   = case.data[itime, :, 14]

    #if len(i) != len(np.unique(i)):
        #print('ibar = %s' % i)
        #print('  eids = %s' % eids)
        #msg = 'ibar=%s is not unique' % str(i)
        #raise RuntimeError(msg)

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
def get_bar100_stress_strain(model: OP2,
                             key,
                             is_stress: bool, vm_word: str, itime: int,
                             oxx,
                             max_principal, min_principal, ovm, is_element_on,
                             eids: np.ndarray,
                             header_dict,
                             keys_map: KeysMap,
                             log: SimpleLogger) -> str:
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        stress = model.op2_results.stress
        bars2 = stress.cbar_stress_10nodes
    else:
        strain = model.op2_results.strain
        bars2 = strain.cbar_strain_10nodes

    if key in bars2:
        case = bars2[key]
        if case.is_complex:
            pass
        else:
            dt = case._times[itime]
            header = _get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header
            keys_map[key] = KeyMap(case.subtitle, case.label,
                                   case.superelement_adaptivity_index,
                                   case.pval_step)

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
                print('  eids = %s' % eids)
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

def get_beam_stress_strain(model: OP2, key, is_stress: bool, vm_word: str, itime: int,
                           oxx,
                           max_principal, min_principal, ovm, is_element_on,
                           header_dict,
                           keys_map: KeysMap,
                           eid_map,
                           log: SimpleLogger) -> str:
    """helper method for _fill_op2_time_centroidal_stress"""
    if is_stress:
        stress = model.op2_results.stress
        beams = stress.cbeam_stress
    else:
        strain = model.op2_results.strain
        beams = strain.cbeam_strain

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
            keys_map[key] = KeyMap(case.subtitle, case.label,
                                   case.superelement_adaptivity_index,
                                   case.pval_step)
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
                try:
                    eid2 = eid_map[eid]
                except KeyError:
                    log.warning('eid=%s is missing' % eid)
                    continue

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
            del eidsi, ueids, sxc, sxd, sxe, sxf, smax, smin, oxxi, smaxi, smini
    del beams
    return vm_word

def get_plate_stress_strain(model: OP2, key,
                            is_stress: bool, vm_word: str, itime: int,
                            oxx, oyy, txy, max_principal, min_principal, ovm, is_element_on,
                            eids: np.ndarray,
                            header_dict,
                            keys_map: KeysMap,
                            log: SimpleLogger) -> str:
    """
    helper method for _fill_op2_time_centroidal_stress.
    Gets the max/min stress across all layers.
    """
    if is_stress:
        stress = model.op2_results.stress
        plates = [
            stress.ctria3_stress, stress.cquad4_stress,
            stress.ctria6_stress, stress.cquad8_stress,
            stress.ctriar_stress, stress.cquadr_stress,
        ]
    else:
        strain = model.op2_results.strain
        plates = [
            strain.ctria3_strain, strain.cquad4_strain,
            strain.ctria6_strain, strain.cquad8_strain,
            strain.ctriar_strain, strain.cquadr_strain,
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
            print('  eids = %s' % eids)
            print('  eidsiA = %s' % case.element_node[:, 0])
            print('  eidsiB = %s' % eidsi)

            msg = 'i%s (plate)=%s is not unique' % (case.element_name, str(i))
            #msg = 'iplate=%s is not unique' % str(i)
            log.warning(msg)
            continue
            #raise RuntimeError(msg)
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
        #                                         txy, angle,
        #                                         majorP, minorP, ovm]
        #print('iplate = ', i)
        #print('  is_element_on = ', is_element_on)
        try:
            is_element_on[i] = 1
        except IndexError:
            #print(case.element_node)
            print('iplate = %s' % i)
            print('  eids = %s' % eids)
            print('  eidsi = %s' % eidsi)
            continue
            #raise

        ntotal = case.data.shape[1]  # (ndt, ntotal, nresults)
        if nlayers_per_element == 1:
            j = None
            nj = 2
        else:
            j = np.arange(ntotal)[::nlayers_per_element]
            nj = len(j)

        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
        #                                         txy, angle,
        #                                         majorP, minorP, ovm]
        dt = case._times[itime]
        header = _get_nastran_header(case, dt, itime)
        header_dict[(key, itime)] = header
        keys_map[key] = KeyMap(case.subtitle, case.label,
                               case.superelement_adaptivity_index,
                               case.pval_step)
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
            oxxi = abs_min_max(np.vstack([oxxi, case.data[itime, j + inode, 1]]), axis=0)
            oyyi = abs_min_max(np.vstack([oyyi, case.data[itime, j + inode, 2]]), axis=0)
            txyi = abs_min_max(np.vstack([txyi, case.data[itime, j + inode, 3]]), axis=0)
            o1i = np.amax(np.vstack([o1i, case.data[itime, j + inode, 5]]), axis=0)
            o3i = np.amin(np.vstack([o3i, case.data[itime, j + inode, 6]]), axis=0)
            ovmi_vstacked = np.vstack([ovmi, case.data[itime, j + inode, 7]])
            if np.any(np.isfinite(ovmi_vstacked)):
                ovmi = np.nanmax(ovmi_vstacked, axis=0)
                #print("---------------")
                #print(f'ovmi:\n{ovmi}')
                #print(f'itime={itime} j={j} inode={inode}\ncase.data[itime, j + inode, 7]:\n{case.data[itime, j + inode, 7]}')
                #print(f'ovmi_vstacked:\n{ovmi_vstacked}')
                #print('type', type(ovmi))
            else:
                if nj is None:
                    nj, unused_nlayers = ovmi.shape
                ovmi = np.full(nj, np.nan, dtype=ovmi.dtype)
                #print("---------------")
                #print(f'ovmi:\n{ovmi}')
                #print(f'itime={itime} j={j} inode={inode}\ncase.data[itime, j + inode, 7]:\n{case.data[itime, j + inode, 7]}')
                #print(f'ovmi_vstacked:\n{ovmi_vstacked}')
                #print('type', type(ovmi))
                #raise

            #print('ovmi', ovmi, ovmi.shape)
            #print('inode', j, j.shape)
            #print('inode', inode)
            #print('ovmi_vstacked.shape', ovmi_vstacked.shape)
            #print('ovmi2', ovmi2)
            #aa
            assert len(oxxi) == len(j)
            #print('-------')

        oxx[i] = oxxi
        oyy[i] = oyyi
        txy[i] = txyi
        max_principal[i] = o1i
        min_principal[i] = o3i
        ovm[i] = ovmi
    return vm_word

def get_solid_stress_strain(model: OP2, key, is_stress: bool, vm_word: str, itime: int,
                            oxx, oyy, ozz, txy, tyz, txz,
                            max_principal, mid_principal, min_principal, ovm, is_element_on,
                            eids: np.ndarray,
                            header_dict,
                            keys_map: KeysMap,
                            log: SimpleLogger) -> str:
    """helper method for _fill_op2_time_centroidal_stress"""
    vm_word0 = vm_word
    if is_stress:
        stress = model.op2_results.stress
        solids = [stress.ctetra_stress,
                  stress.cpenta_stress,
                  stress.chexa_stress,]
    else:
        strain = model.op2_results.strain
        solids = [strain.ctetra_strain,
                  strain.cpenta_strain,
                  strain.chexa_strain,]

    for result in solids:
        if key not in result:
            continue
        case = result[key]
        if case.is_complex:
            continue

        if isinstance(case, RealSolidArrayNx):
            log.info(f'converting {case.class_name}')
            case = case.to_real_solid_array()
            result[key] = case

        if case.is_von_mises:
            vm_word = 'vonMises'
        else:
            vm_word = 'maxShear'

        nnodes_per_element = case.nnodes
        eidsi = case.element_cid[:, 0]

        common_eids = np.intersect1d(eids, eidsi)
        ifilter = np.array([], dtype='int32')
        if len(common_eids) != len(eidsi):
            ifilter = np.searchsorted(eidsi, common_eids)

        if isinstance(case, RealSolidArray):
            out = _get_solid_stress_strain_real(
                eids,
                eidsi, case,
                nnodes_per_element,
                ifilter,
                key, itime,
                header_dict,
                keys_map,
                log,
            )
        else:
            raise NotImplementedError(case.class_name)
        eidsi, oxxi, oyyi, ozzi, txyi, tyzi, txzi, o1i, o2i, o3i, ovmi = out

        common_eids = np.intersect1d(eids, eidsi)
        if len(common_eids) == 0:
            return vm_word0

        i = np.searchsorted(eids, eidsi)
        #j = np.searchsorted(eidsi, common_eids)
        if len(i) != len(np.unique(i)):
            print('isolid = %s' % str(i))
            print('  eids = %s' % eids)
            print('  eidsi = %s' % eidsi)
            if len(i) != len(np.unique(i)):
                msg = 'i%s (solid)=%s is not unique' % (case.element_name, str(i))
                #msg = 'isolid=%s is not unique' % str(i)
                log.warning(msg)
                continue

        is_element_on[i] = 1
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

def _get_solid_stress_strain_real(eids: np.ndarray,
                                  eidsi: np.ndarray,
                                  case: RealSolidArray,
                                  nnodes_per_element: int,
                                  ifilter: np.ndarray,
                                  key: str, itime: int,
                                  header_dict,
                                  keys_map: KeysMap,
                                  log: SimpleLogger,
                                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                             np.ndarray, np.ndarray, np.ndarray,
                                             np.ndarray, np.ndarray, np.ndarray,
                                             np.ndarray, np.ndarray]:
    """
    Gets the worst/extreme solid stress/strain for each element and calls that
    the value.

    For [0, 1, 2, 2, -2.1], the worst stress is -2.1 because it has the
    largest absolute value (2.1), despite the average being positive.
    """
    ntotal = len(eidsi)  * nnodes_per_element

    #self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz,
    #                                         txy, tyz, txz,
    #                                         o1, o2, o3, ovm]

    dt = case._times[itime]
    header = _get_nastran_header(case, dt, itime)
    header_dict[(key, itime)] = header
    keys_map[key] = KeyMap(case.subtitle, case.label,
                           case.superelement_adaptivity_index,
                           case.pval_step)

    ntotal_filtered = ntotal
    nelement_filtered = len(ifilter)
    data = case.data
    eids_out = eidsi
    if nelement_filtered:
        ntotal_filtered = nelement_filtered * nnodes_per_element
        ntime, nelement_nnode, nresult = case.data.shape
        nelement = nelement_nnode // nnodes_per_element
        data2 = case.data.reshape(ntime, nelement, nnodes_per_element, nresult)[:, ifilter, :, :]
        data3 = data2.reshape(ntime, ntotal_filtered, nresult)
        data = data3
        common_eids = np.intersect1d(eids, eidsi)
        #j = np.searchsorted(eidsi, common_eids) * nnodes_per_element
        j = np.arange(ntotal_filtered)[::nnodes_per_element]
        eids_out = common_eids
    else:
        if nnodes_per_element == 1:
            j = None
        else:
            j = np.arange(ntotal_filtered)[::nnodes_per_element]
            ueidsi = np.unique(eidsi)
            if len(j) != len(ueidsi):
                log.warning('solid stress/strain is missing elements; j=%s ueidsi=%s' % (j, ueidsi))

    # stack em up
    oxx_list = [data[itime, j, 0]]
    oyy_list = [data[itime, j, 1]]
    ozz_list = [data[itime, j, 2]]
    txy_list = [data[itime, j, 3]]
    tyz_list = [data[itime, j, 4]]
    txz_list = [data[itime, j, 5]]
    o1_list  = [data[itime, j, 6]] # max
    o2_list  = [data[itime, j, 7]] # min
    o3_list  = [data[itime, j, 8]] # mid
    ovm_list = [data[itime, j, 9]]
    for inode in range(1, nnodes_per_element):
        j_inode = j + inode
        oxxii = case.data[itime, j_inode, 0]
        oyyii = case.data[itime, j_inode, 1]
        ozzii = case.data[itime, j_inode, 2]
        txyii = case.data[itime, j_inode, 3]
        tyzii = case.data[itime, j_inode, 4]
        txzii = case.data[itime, j_inode, 5]
        o1ii  = case.data[itime, j_inode, 6] # max
        o2ii  = case.data[itime, j_inode, 7] # min
        o3ii  = case.data[itime, j_inode, 8] # mid
        ovmii = case.data[itime, j_inode, 9]
        oxx_list.append(oxxii)
        oyy_list.append(oyyii)
        ozz_list.append(ozzii)
        txy_list.append(txyii)
        tyz_list.append(tyzii)
        txz_list.append(txzii)
        o1_list.append(o1ii)
        o2_list.append(o2ii)
        o3_list.append(o3ii)
        ovm_list.append(ovmii)

    # smush them into a vector and get the worst stress
    # we vstack'd, so we use axis=0 (vs axis=1) because it's 2D
    oxxi = abs_min_max(np.vstack(oxx_list), axis=0)
    oyyi = abs_min_max(np.vstack(oyy_list), axis=0)
    ozzi = abs_min_max(np.vstack(ozz_list), axis=0)
    txyi = abs_min_max(np.vstack(txy_list), axis=0)
    tyzi = abs_min_max(np.vstack(tyz_list), axis=0)
    txzi = abs_min_max(np.vstack(txz_list), axis=0)

    o1i = np.amax(np.vstack(o1_list), axis=0) # max
    o2i = abs_min_max(np.vstack(o2_list), axis=0) # mid
    o3i = np.amin(np.vstack(o3_list), axis=0) # min
    ovmi = np.amax(np.vstack(ovm_list), axis=0) # max

    return eids_out, oxxi, oyyi, ozzi, txyi, tyzi, txzi, o1i, o2i, o3i, ovmi

#def _get_solid_stress_strain_real_nx(case: RealSolidArrayNx,
                                     #eidsi,
                                     #nnodes_per_element: int,
                                     #key: str, itime: int,
                                     #header_dict,
                                     #keys_map):

    #ntotal = len(eidsi)  * nnodes_per_element

    ##self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz,
    ##                                         txy, tyz, txz,
    ##                                         o1, o2, o3, ovm]

    #if nnodes_per_element == 1:
        #j = None
    #else:
        #j = np.arange(ntotal)[::nnodes_per_element]
        #ueidsi = np.unique(eidsi)
        #assert len(j) == len(ueidsi), 'j=%s ueidsi=%s' % (j, ueidsi)

    #dt = case._times[itime]
    #header = _get_nastran_header(case, dt, itime)
    #header_dict[(key, itime)] = header
    #keys_map[key] = (case.subtitle, case.label,
                     #case.superelement_adaptivity_index, case.pval_step)
    #oxxi = case.data[itime, j, 0]
    #oyyi = case.data[itime, j, 1]
    #ozzi = case.data[itime, j, 2]
    #txyi = case.data[itime, j, 3]
    #tyzi = case.data[itime, j, 4]
    #txzi = case.data[itime, j, 5]
    #ovmi = case.data[itime, j, 6]

    #nrows = len(oxxi)
    #a_matrix = np.full((nrows, 3, 3), np.nan)
    #a_matrix[:, 0, 0] = oxxi
    #a_matrix[:, 1, 1] = oyyi
    #a_matrix[:, 2, 2] = ozzi

    #a_matrix[:, 0, 1] = a_matrix[:, 1, 0] = txyi
    #a_matrix[:, 0, 2] = a_matrix[:, 2, 0] = txzi
    #a_matrix[:, 1, 2] = a_matrix[:, 2, 1] = tyzi

    #eigs = np.linalg.eigvalsh(a_matrix)  # array = (..., M, M) array

    #o1i = eigs[:, 2]
    #o2i = eigs[:, 1]
    #o3i = eigs[:, 0]

    #for inode in range(nnodes_per_element-1):
        #oxxi = np.amax(np.vstack([oxxi, case.data[itime, j + inode, 0]]), axis=0)
        #oyyi = np.amax(np.vstack([oyyi, case.data[itime, j + inode, 1]]), axis=0)
        #ozzi = np.amax(np.vstack([ozzi, case.data[itime, j + inode, 2]]), axis=0)
        #txyi = np.amax(np.vstack([txyi, case.data[itime, j + inode, 3]]), axis=0)
        #tyzi = np.amax(np.vstack([tyzi, case.data[itime, j + inode, 4]]), axis=0)
        #txzi = np.amax(np.vstack([txzi, case.data[itime, j + inode, 2]]), axis=0)

        ##o1i = np.amax(np.vstack([o1i, case.data[itime, j + inode, 6]]), axis=0)
        ##o2i = np.amax(np.vstack([o2i, case.data[itime, j + inode, 7]]), axis=0)
        ##o3i = np.amin(np.vstack([o3i, case.data[itime, j + inode, 8]]), axis=0)
        ##o1i = np.zeros(ovmi.shape)
        #ovmi = np.amax(np.vstack([ovmi, case.data[itime, j + inode, 6]]), axis=0)
        #assert len(oxxi) == len(j)
    #return oxxi, oyyi, ozzi, txyi, tyzi, txzi, o1i, o2i, o3i, ovmi
