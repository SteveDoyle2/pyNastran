from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.dev.bdf_vectorized3.solver.utils import get_ieids_eids, get_element
from pyNastran.op2.op2_interface.op2_classes import (
    RealStrainEnergyArray,
    RealSpringStrainArray, RealSpringStressArray, RealSpringForceArray,
)
from .utils import save_strain_energy # get_plot_request, 

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdfvectorized3.bdf import BDF

def recover_celas(f06_file, op2,
                  model: BDF, dof_map, isubcase, xg, eids_str,
                  element_name: str,
                  write_f06_strain: bool, write_op2_strain: bool,
                  write_f06_stress: bool, write_op2_stress: bool,
                  write_f06_force: bool, write_op2_force: bool,
                  write_f06_ese: bool, write_op2_ese: bool,
                  fdtype='float32',
                  title: str='', subtitle: str='', label: str='',
                  page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static spring strain"""
    neids, ielas, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids
    xg, nmode = _fix_xg(xg)

    get_strain = write_f06_strain or write_op2_strain
    get_stress = write_f06_stress or write_op2_stress
    get_force = write_f06_force or write_op2_force
    get_strain_energy = write_f06_ese or write_op2_ese

    force = None
    stress = None
    strain = None
    strain_energy = None

    if get_force:
        force = np.full((neids, 1), np.nan, dtype=fdtype)
    if get_stress:
        stress = np.full((neids, 1), np.nan, dtype=fdtype)
    if get_strain:
        strain = np.full((neids, 1), np.nan, dtype=fdtype)
    if get_strain_energy:
        strain_energy = np.full((neids, 1), np.nan, dtype=fdtype)

    if element_name == 'CELAS1':
        #strain_slot = 'celas1_strain'
        #stress_slot = 'celas1_stress'
        #force_slot = 'celas1_force'
        #ese_slot = 'celas1_strain_energy'

        celas1_mapper = {
            # get_strain, get_stress, get_force, get_strain_energy
            # 2^4 = 16
            # ostr, oes, oef,  ese
            (True, True, True, True): get_celas1_ostr_oes_oef_ese,
            (True, True, True, False): get_celas1_ostr_oes_oef,
            (True, True, False, True): get_celas1_ostr_oes_ese,
            (True, True, False, False): get_celas1_ostr_oes,
            (True, False, True, True): get_celas1_ostr_oef_ese,
            (True, False, True, False): get_celas1_ostr_oef,
            (True, False, False, True): get_celas1_ostr_ese,
            (True, False, False, False): get_celas1_ostr,

            #(False, True, True, True): get_celas1_ostr_oes_oef_ese[1:],
            #(False, True, True, False): get_celas1_ostr_oes_oef[1:],
            #(False, True, False, True): get_celas1_ostr_oes_ese[1:],
            #(False, True, False, False): get_celas1_ostr_oes[1:],
            #(False, False, True, True): get_celas1_ostr_oef_ese[1:],
            #(False, False, True, False): get_celas1_ostr_oef[1:],
            #(False, False, False, True): get_celas1_ostr_ese[1:],

            # ostr, oes, oef,  ese
            (False, True, True, True): get_celas1_oes_oef_ese,
            (False, True, True, False): get_celas1_oes_oef,
            (True, True, False, True):  get_celas1_oes_ese,
            (False, True, False, False): get_celas1_oes,
            (False, False, True, True): get_celas1_oef_ese,
            (False, False, True, False): get_celas1_oef,
            (False, False, False, True): get_celas1_ese,
            #(False, False, False, False): passer,
        }
        celeas1_func = celas1_mapper[(get_strain, get_stress, get_force, get_strain_energy)]
        # --------------------------------------------------------------
        # strain
        for unused_key, celeas1_func in celas1_mapper.items():
            celeas1_func(model, ielas, eids, xg, dof_map,
                         strain, stress, force, strain_energy)

    elif element_name == 'CELAS2':
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            si = elem.s
            strain[ieid] = _recover_straini_celas12(xg, dof_map, elem, si)
    elif element_name == 'CELAS3':
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            si = elem.pid_ref.s
            strain[ieid] = _recover_straini_celas34(xg, dof_map, elem, si)
    elif element_name == 'CELAS4':
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            si = 1.0 # TODO: is this right?
            strain[ieid] = _recover_straini_celas34(xg, dof_map, elem, si)
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    case = None
    _save_spring_stress_strain_force(
        op2.op2_results.strain, f06_file, page_num, page_stamp,
        element_name,
        strain, eids,
        isubcase, title, subtitle, label,
        nmode, result_type='strain',
        case=case,
        write_op2=True,
        write_f06=write_f06_strain)

    _save_spring_stress_strain_force(
        op2.op2_results.stress, f06_file, page_num, page_stamp,
        element_name,
        stress, eids,
        isubcase, title, subtitle, label,
        nmode, result_type='stress',
        case=case,
        write_op2=True,
        write_f06=write_f06_stress)

    _save_spring_stress_strain_force(
        op2.op2_results.force, f06_file, page_num, page_stamp,
        element_name,
        force, eids,
        isubcase, title, subtitle, label,
        nmode, result_type='force',
        case=case,
        write_op2=True,
        write_f06=write_f06_force)

    save_strain_energy(
        op2, f06_file, page_num, page_stamp,
        element_name,
        strain_energy, eids,
        isubcase, title, subtitle, label,
        #case=case,
        write_op2=True,
        write_f06=write_f06_ese)
    return neids

def get_celas1_ostr_oes_ese(model, ielas, eids, xg, dof_map,
                            strain, stress, unused_force, strain_energy):
    asdf
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        straini, unused_forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain[ieid] = straini
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)


def get_celas1_ostr_ese(model: BDF, ielas, eids, xg, dof_map,
                        strain, unused_stress, unused_force, strain_energy):
    sasdf
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        straini, unused_forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        strain[ieid] = straini
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_ostr_oef(model: BDF, ielas, eids, xg, dof_map,
                        strain, unused_stress, force, unused_strain_energy):
    asdf
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        straini, forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        force[ieid] = forcei
        strain[ieid] = straini

def get_celas1_oes_ese(model: BDF, ielas, eids, xg, dof_map,
                       unused_strain, stress, unused_force, strain_energy):
    asdf
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_ostr_oes_oef_ese(model: BDF, ielas, eids, xg, dof_map,
                                strain, stress, force, strain_energy):
    asdf
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        straini, forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        force[ieid] = forcei
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain[ieid] = straini
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_ostr_oes_oef(model: BDF, ielas, eids, xg, dof_map,
                            strain, stress, force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        straini, forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        force[ieid] = forcei
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain[ieid] = straini

def get_celas1_ostr_oef_ese(model: BDF, ielas, eids, xg, dof_map,
                            strain, stress, force, unused_strain_energy):
    asdf
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        straini, forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        force[ieid] = forcei
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain[ieid] = straini

def get_celas1_ostr_oes(model: BDF, ielas, eids, xg, dof_map,
                        strain, stress, unused_force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain[ieid] = _recover_straini_celas12(xg, dof_map, elem, si)

def get_celas1_ostr(model: BDF, ielas, eids, xg, dof_map,
                    strain, unused_stress, unused_force, unused_strain_energy):
    asdf
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        si = elem.pid_ref.s
        strain[ieid] = _recover_straini_celas12(xg, dof_map, elem, si)

def get_celas1_oes_oef_ese(model: BDF, ielas, eids, xg, dof_map,
                           unused_strain, stress, force, strain_energy):
    asdf
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        force[ieid] = _recover_forcei_celas12(xg, dof_map, elem, ki)
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_oes_oef(model: BDF, ielas, eids, xg, dof_map,
                       unused_strain, stress, force, unused_strain_energy):
    asdf
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        force[ieid] = _recover_forcei_celas12(xg, dof_map, elem, ki)
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)

def get_celas1_oes(model: BDF, ielas, eids, xg, dof_map,
                   unused_strain, stress, unused_force, unused_strain_energy):
    asdf
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)

def get_celas1_oef_ese(model: BDF, ielas, eids, xg, dof_map,
                       unused_strain, unused_stress, force, strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        force[ieid] = _recover_forcei_celas12(xg, dof_map, elem, ki)
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_oef(model: BDF, ielas, eids, xg, dof_map,
                   unused_strain, unused_stress, force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        force[ieid] = _recover_forcei_celas12(xg, dof_map, elem, ki)

def get_celas1_ese(model: BDF, ielas, eids, xg, dof_map,
                   unused_strain, unused_stress, unused_force, strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        si = elem.pid_ref.s
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def _celas_ks(model: BDF,
              elem,
              element_name: str,
              ieids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if element_name in {'CELAS1', 'CELAS3'}:
        pid = elem.property_id[ieids]
        pelas = model.pelas.slice_card_by_property_id(pid)
        k = pelas.k.ravel()
        s = pelas.s.ravel()
    elif element_name == 'CELAS2':
        k = elem.k[ieids].ravel()
        s = elem.s[ieids].ravel()
    elif element_name == 'CELAS4':  # good
        k = elem.k[ieids].ravel()
        s = 0
        #raise RuntimeError(elem.get_stats())
    else:
        raise RuntimeError(element_name)
    assert len(k) == len(elem)
    return k, s


def _get_recovery_tables(subcase: Subcase) -> tuple[bool, tuple, tuple]:
    is_stress = 'STRESS' in subcase
    is_strain = 'STRAIN' in subcase
    is_force = 'FORCE' in subcase
    is_strain_energy = 'ESE' in subcase
    is_flags = (is_force, is_stress, is_strain, is_strain_energy)

    if not(is_force or is_stress or is_strain or is_strain_energy):
        str_flags = ('NONE', 'NONE', 'NONE', 'NONE')
        return False, is_flags, str_flags

    force_eids_str = 'NONE'
    stress_eids_str = 'NONE'
    strain_eids_str = 'NONE'
    ese_eids_str = 'NONE'
    if is_force:
        force_eids_str = subcase['FORCE'][0]
    if is_stress:
        stress_eids_str = subcase['STRESS'][0]
    if is_strain:
        strain_eids_str = subcase['STRAIN'][0]
    if is_ese:
        ese_eids_str, ese_options = subcase['ESE']
        # tiny = ese_options['TINY']
    str_flags = (force_eids_str, stress_eids_str, strain_eids_str, ese_eids_str)
    return True, is_flags, str_flags


def _recover_celas(model: BDF,
                   subcase: Subcase,
                   dof_map: dict[int, int],
                   xg: np.ndarray,
                   element_name: str, fdtype='float32') -> None:
    """recovers all celas terms"""
    is_result, is_flags, str_flags = _get_recovery_tables(subcase)
    if not is_result:
        return 0, {}
    is_flags = is_force, is_stress, is_strain, is_strain_energy
    force_eids_str, stress_eids_str, strain_eids_str, ese_eids_str = str_flags
    
    elem = get_element(model, element_name, ieids, eids)
    if is_force or is_stress or is_strain_energy:
        k, s = _celas_ks(model, elem, element_name, ieids)
    dx = _spring_dx(
        xg, dof_map,
        elem, element_name, neids, ieids, eids,
        fdtype=fdtype)
    element_id = elem.element_id

    force = None
    stress = None
    strain_energy = None
    obj_name = element_name.lower()
    straini = dx
    forcei = k * dx
    stressi = k * s * dx  # this seems wrong
    strain_energyi = 0.5 * ki * dx ** 2   # 0.5 * F * dx

    force, stress, strain, strain_energy = _slice_results(
        model, element_name, element_id,
        is_force, force, force_eids_str,
        is_stress, stress, stress_eids_str,
        is_strain, strain, strain_eids_str,
        is_ese, strain_energy, ese_eids_str)

    out = {
        ('force', f'{obj_name}_force'): force,
        ('strain', f'{obj_name}_strain'): strain,
        ('stress', f'{obj_name}_stress'): stress,
        ('strain_energy', f'{obj_name}_strain_energy'): strain_energy,
    }
    assert neids > 0
    return neids, out

    
def _slice_results(model: BDF,
                   element_name: str,
                   element_id: np.ndarray,
                   is_force: bool, force, force_eids_str: str,
                   is_stress: bool, stress, stress_eids_str: str,
                   is_strain: bool, strain, strain_eids_str: str,
                   is_ese: bool, strain_energy, ese_eids_str: str,
                   ):
    neid = len(element_id)
    if is_force:
        force_neids, force_ieids, force_eids = get_ieids_eids(
            model, element_name, force_eids_str)
        force = forcei if force_neids == neid else forcei[force_ieids, :]
        force = None if len(force) == 0 else force

    if is_stress:
        stress_neids, stress_ieids, stress_eids = get_ieids_eids(
            model, element_name, stress_eids_str)
        stress = stressi if stress_neids == neid else stressi[stress_ieids, :]
        stress = None if len(stress) == 0 else stress
    if is_strain:
        strain_neids, strain_ieids, strain_eids = get_ieids_eids(
            model, element_name, strain_eids_str)
        strain = straini if strain_neids == neid else straini[strain_ieids, :]
        strain = None if len(strain) == 0 else strain
    if is_ese:
        ese_neids, ese_ieids, ese_eids = get_ieids_eids(
            model, element_name, ese_eids_str)
        strain_energy = esei if ese_neids == neid else esei[ese_ieids, :]
        if strain_energy.sum() == 0.0:
            strain_energy = None
        else:
            strain_energy = None if len(strain_energy) == 0 else strain_energy
    return force, stress, strain, strain_energy


def _fix_xg(xg):
    if xg.ndim == 1:
        ndof = xg.shape
        nmode = 1
    else:
        assert xg.shape == 2, xg.shape
        ndof, nmode = xg.shape
        assert nmode == 1, xg.shape
    return xg, nmode

def _recover_force_celas(f06_file: TextIO, op2: OP2,
                         model: BDF, dof_map: dict[int, int],
                         isubcase: int, xg: np.ndarray, eids_str: str,
                         element_name: str, fdtype='float32',
                         title: str='', subtitle: str='', label: str='',
                         page_num: int=1, page_stamp='PAGE %s',
                         write_f06: bool = True,
                         write_op2: bool = True,
                         case=None) -> None:
    """recovers static/modal spring force"""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    # print(f'eids_str={eids_str} eids={eids} ieids={eids}')
    if not neids:
        return neids
    xg, nmode = _fix_xg(xg)

    model.log.warning(f'xg.shape={xg.shape}')
    elem = get_element(model, element_name, ieids, eids)
    k = _celas_ks(model, elem, element_name, ieids)[0]
    dx = _spring_dx(
        xg, dof_map,
        elem, element_name, neids, ieids, eids,
        fdtype=fdtype)

    force = k[np.newaxis, :, np.newaxis] * dx
    #assert k.shape == dx.shape, (k.shape, dx.shape)

    _save_spring_stress_strain_force(
        op2.op2_results.force, f06_file, page_num, page_stamp,
        element_name,
        force, eids,
        isubcase, title, subtitle, label,
        nmode, result_type='force',
        write_f06=write_f06, write_op2=write_op2,
        case=case)
    assert neids > 0
    return neids

def _recover_stress_celas(f06_file: TextIO, op2: OP2,
                          model: BDF, dof_map: dict[int, int],
                          isubcase: int, xg: np.ndarray, eids_str: str,
                          element_name: str, fdtype='float32',
                          title: str='', subtitle: str='', label: str='',
                          page_num: int=1, page_stamp='PAGE %s',
                          write_f06: bool = True,
                          write_op2: bool = True,
                          case=None) -> None:
    """recovers static spring force"""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    # print(f'eids_str={eids_str} eids={eids} ieids={eids}')
    if not neids:
        return neids
    xg, nmode = _fix_xg(xg)

    elem = get_element(model, element_name, ieids, eids)
    k, s = _celas_ks(model, elem, element_name, ieids)
    dx = _spring_dx(
        xg, dof_map,
        elem, element_name, neids, ieids, eids,
        fdtype=fdtype)

    ks = k * s
    stress = ks[:, np.newaxis, :] * dx  # this seems wrong
    #assert k.shape == dx.shape, (k.shape, dx.shape)

    _save_spring_stress_strain_force(
        op2.op2_results.stress, f06_file, page_num, page_stamp,
        element_name,
        force, eids, write_f06,
        isubcase, title, subtitle, label,
        nmode, result_type='stress',
        write_f06=write_f06, write_op2=write_op2,
        case=case)
    assert neids > 0
    return neids

def _recover_strain_celas(f06_file: TextIO, op2: OP2,
                          model: BDF, dof_map: dict[int, int],
                          isubcase: int, xg: np.ndarray, eids_str: str,
                          element_name: str, fdtype='float32',
                          title: str='', subtitle: str='', label: str='',
                          page_num: int=1, page_stamp='PAGE %s',
                          write_f06: bool = True,
                          write_op2: bool = True,
                          case=None) -> None:
    """recovers static spring force"""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    # print(f'eids_str={eids_str} eids={eids} ieids={eids}')
    if not neids:
        return neids
    xg, nmode = _fix_xg(xg)

    elem = get_element(model, element_name, ieids, eids)
    k, s = _celas_ks(model, elem, element_name, ieids)
    dx = _spring_dx(
        xg, dof_map,
        elem, element_name, neids, ieids, eids,
        fdtype=fdtype)

    strain = dx  # this seems wrong
    assert k.shape == dx.shape, (k.shape, dx.shape)

    _save_spring_stress_strain_force(
        op2.op2_results.strain, f06_file, page_num, page_stamp,
        element_name,
        force, eids, write_f06,
        isubcase, title, subtitle, label,
        nmode, result_type='strain',
        write_f06=write_f06, write_op2=write_op2,
        case=case)
    assert neids > 0
    return neids


def _recover_strain_energy_celas(f06_file, op2,
                                 model: BDF, dof_map, isubcase, xg, eids_str,
                                 element_name: str, fdtype='float32',
                                 title: str='', subtitle: str='', label: str='',
                                 page_num: int=1, page_stamp='PAGE %s',
                                 write_f06: bool = True,
                                 write_op2: bool = True,
                                 case=None) -> None:
    """recovers static spring strain energy"""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids
    xg, nmode = _fix_xg(xg)

    elem = get_element(model, element_name, ieids, eids)
    k = _celas_ks(model, elem, element_name, ieids)[0]
    dx = _spring_dx(
        xg, dof_map,
        elem, element_name, neids, ieids, eids,
        fdtype=fdtype)
    
    # 1/2 * k * x^2
    strain_energy = 0.5 * k[:, np.newaxis] * dx ** 2
    assert np.all(np.isfinite(k)), k
    assert np.all(np.isfinite(strain_energy)), strain_energy

    save_strain_energy(
        op2, f06_file, page_num, page_stamp,
        element_name,
        strain_energy, eids,
        isubcase, title, subtitle, label,
        write_f06=write_f06, write_op2=write_op2,
        case=case)
    return neids


def _spring_dx(xg: np.ndarray,
               dof_map: DOF_MAP,
               elem,
               element_name: str,
               neids: int,
               ieids: np.ndarray,
               eids: np.ndarray,
               fdtype: str='float32'):
    if xg.ndim == 1:
        ndof = len(xg)
        nmode = 1
        xg = xg.reshape(ndof, nmode)
    else:
        assert xg.shape == 2, xg.shape
        ndof, nmode = xg.shape
    #    assert nmode == 1, xg.shape
    #ndof, nmode = xg.shape

    dx = np.full((neids, nmode), np.nan, dtype=fdtype)
    if element_name in {'CELAS1', 'CELAS2'}:
        for ieid, eid, nodes, comps in zip(ieids, eids, elem.nodes, elem.components):
            nid1, nid2 = nodes
            c1, c2 = comps
            i = dof_map[(nid1, c1)]
            j = dof_map[(nid2, c2)]
            dx[ieid, :] = xg[j, :] - xg[i, :]  # TODO: check the sign
    elif element_name in {'CELAS3', 'CELAS4'}:
        for ieid, eid, nodes in zip(ieids, eids, elem.spoints):
            nid1, nid2 = nodes
            i = dof_map[(nid1, 0)]
            j = dof_map[(nid2, 0)]
            dx[ieid, :] = xg[j, :] - xg[i, :]  # TODO: check the sign
    else:  # pragma: no cover
        raise NotImplementedError(element_name)
    assert np.all(np.isfinite(dx)), dx

    #if nmode == 1:
    #    dx = dx.reshape(1, neids, 1)
    return dx


def _save_spring_stress_strain_force(
        results, f06_file, page_num, page_stamp,
        element_name: str,
        strain, eids,
        isubcase: int, title: str, subtitle: str,
        label: str,
        nmode: int,
        result_type: str,
        case=None,
        write_op2: bool=True,
        write_f06: bool=True) -> None:
    if strain is None:
        return
    if case is None:
        case = {}
    data = strain
    #data = strain.reshape(1, *strain.shape)
 
    if result_type == 'stress':
        table_name = 'OES1'
        cls = RealSpringStressArray
    elif result_type == 'strain':
        table_name = 'OSTR1'
        cls = RealSpringStrainArray
    elif result_type == 'force':
        table_name = 'OEF1'
        cls = RealSpringForceArray
    else:  # pragma: no cover
        raise RuntimeError(result_type)

    casetype = case.get('type', 'static')
    if casetype == 'static':
        assert nmode == 1, data.shape
        spring_obj = cls.add_static_case(
            table_name, element_name, eids, data, isubcase,
            is_sort1=True, is_random=False, is_msc=True,
            random_code=0, title=title, subtitle=subtitle, label=label)
    elif casetype == 'modes':
        eigenvalues = case['eigenvalue']
        spring_obj = cls.add_modal_case(
            table_name, element_name, eids, data, isubcase,
            is_sort1=True, is_random=False, is_msc=True,
            random_code=0, title=title, subtitle=subtitle, label=label)
    elif casetype == 'freq':
        frequency = case['frequency']
        spring_obj = cls.add_freq_case(
            table_name, element_name, eids, data, isubcase,
            is_sort1=True, is_random=False, is_msc=True,
            random_code=0, title=title, subtitle=subtitle, label=label)
    else:  # pragma: no cover
        raise RuntimeError(case)

    obj_dict = getattr(results, f'{element_name.lower()}_{result_type}')
    if write_op2:
        obj_dict[isubcase] = spring_obj
    if write_f06:
        spring_obj.write_f06(
            f06_file, header=None, page_stamp=page_stamp,
            page_num=page_num, is_mag_phase=False, is_sort1=True)
