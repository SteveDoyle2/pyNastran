from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.dev.bdf_vectorized3.solver.utils import get_ieids_eids, get_element
from pyNastran.op2.op2_interface.op2_classes import (
    RealStrainEnergyArray,
    RealSpringStrainArray, RealSpringStressArray, RealSpringForceArray,
)
#from .utils import get_plot_request

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
    _save_spring_strain(
        op2, f06_file, page_num, page_stamp,
        element_name,
        strain, eids, write_f06_strain,
        isubcase, title, subtitle, label)
    _save_spring_stress(
        op2, f06_file, page_num, page_stamp,
        element_name,
        stress, eids, write_f06_stress,
        isubcase, title, subtitle, label)
    _save_spring_force(
        op2, f06_file, page_num, page_stamp,
        element_name,
        force, eids, write_f06_force,
        isubcase, title, subtitle, label)
    _save_strain_energy(
        op2, f06_file, page_num, page_stamp,
        element_name,
        strain_energy, eids, write_f06_ese,
        isubcase, title, subtitle, label)
    return neids

def get_celas1_ostr_oes_ese(model, ielas, eids, xg, dof_map,
                            strain, stress, unused_force, strain_energy):
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
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        straini, unused_forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        strain[ieid] = straini
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_ostr_oef(model: BDF, ielas, eids, xg, dof_map,
                        strain, unused_stress, force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        straini, forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        force[ieid] = forcei
        strain[ieid] = straini

def get_celas1_oes_ese(model: BDF, ielas, eids, xg, dof_map,
                       unused_strain, stress, unused_force, strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_ostr_oes_oef_ese(model: BDF, ielas, eids, xg, dof_map,
                                strain, stress, force, strain_energy):
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
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        si = elem.pid_ref.s
        strain[ieid] = _recover_straini_celas12(xg, dof_map, elem, si)

def get_celas1_oes_oef_ese(model: BDF, ielas, eids, xg, dof_map,
                           unused_strain, stress, force, strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        force[ieid] = _recover_forcei_celas12(xg, dof_map, elem, ki)
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_oes_oef(model: BDF, ielas, eids, xg, dof_map,
                       unused_strain, stress, force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        force[ieid] = _recover_forcei_celas12(xg, dof_map, elem, ki)
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)

def get_celas1_oes(model: BDF, ielas, eids, xg, dof_map,
                   unused_strain, stress, unused_force, unused_strain_energy):
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

def _celas_ks(model: BDF, element_name: str,
              elem, ieids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if element_name in {'CELAS1', 'CELAS3'}:
        pid = elem.property_id[ieids]
        pelas = model.pelas.slice_card_by_property_id(pid)
        k = pelas.k.ravel()
        s = pelas.s.ravel()
    else:
        k = elem.k[ieids].ravel()
        s = elem.s[ieids].ravel()
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
    
    if is_force or is_stress or is_strain_energy:
        k, s = _celas_ks(model, element_name, elem, ieids)

    elem, dx = _spring_dx(
        element_name, neids, ieids, eids,
        dtype=dtype)
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


def _recover_force_celas(f06_file: TextIO, op2: OP2,
                         model: BDF, dof_map: dict[int, int],
                         isubcase: int, xg: np.ndarray, eids_str: str,
                         element_name: str, fdtype='float32',
                         title: str='', subtitle: str='', label: str='',
                         page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static spring force"""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    # print(f'eids_str={eids_str} eids={eids} ieids={eids}')
    if not neids:
        return neids

    elem = get_element(model, element_name, ieids, eids)
    k = _celas_ks(model, element_name, elem, ieids)[0]

    elem, dx = _spring_dx(
        element_name, neids, ieids, eids,
        dtype=dtype)
    force = k * dx

    write_f06_force = True
    _save_spring_force(
        op2, f06_file, page_num, page_stamp,
        element_name,
        force, eids, write_f06_force,
        isubcase, title, subtitle, label)
    assert neids > 0
    return neids

def _recover_stress_celas(f06_file: TextIO, op2: OP2,
                          model: BDF, dof_map: dict[int, int],
                          isubcase: int, xg: np.ndarray, eids_str: str,
                          element_name: str, fdtype='float32',
                          title: str='', subtitle: str='', label: str='',
                          page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static spring force"""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    # print(f'eids_str={eids_str} eids={eids} ieids={eids}')
    if not neids:
        return neids

    k, s = _celas_ks(model, element_name, elem, ieids)
    elem, dx = _spring_dx(
        element_name, neids, ieids, eids,
        dtype=dtype)

    stress = k * s * dx  # this seems wrong

    write_f06_stress = True
    _save_spring_stress(
        op2.op2_results.stress, f06_file, page_num, page_stamp,
        element_name,
        force, eids, write_f06_stress,
        isubcase, title, subtitle, label)
    assert neids > 0
    return neids

def _recover_strain_celas(f06_file: TextIO, op2: OP2,
                          model: BDF, dof_map: dict[int, int],
                          isubcase: int, xg: np.ndarray, eids_str: str,
                          element_name: str, fdtype='float32',
                          title: str='', subtitle: str='', label: str='',
                          page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static spring force"""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    # print(f'eids_str={eids_str} eids={eids} ieids={eids}')
    if not neids:
        return neids

    # k, s = _celas_ks(model, element_name, elem, ieids)
    elem, dx = _spring_dx(
        element_name, neids, ieids, eids,
        dtype=dtype)
    strain = dx  # this seems wrong

    write_f06_strain = True
    _save_spring_strain(
        op2.op2_results, f06_file, page_num, page_stamp,
        element_name,
        force, eids, write_f06_stress,
        isubcase, title, subtitle, label)
    assert neids > 0
    return neids


def _recover_strain_energy_celas(f06_file, op2,
                                 model: BDF, dof_map, isubcase, xg, eids_str,
                                 element_name: str, fdtype='float32',
                                 title: str='', subtitle: str='', label: str='',
                                 page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static spring strain energy"""
    neids, ieids, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    write_f06_ese = True
    strain_energies = np.full((neids, 1), np.nan, dtype=fdtype)

    elem = get_element(model, element_name, ieids, eids)
    k = _celas_ks(model, element_name, elem, ieids)[0]
    elem, dx = _spring_dx(
        element_name, neids, ieids, eids,
        dtype=dtype)
    
    # 1/2 * k * x^2
    strain_energy = 0.5 * ki * dx ** 2

    _save_spring_strain_energy(
        op2, f06_file, page_num, page_stamp,
        element_name,
        strain_energies, eids, write_f06_ese,
        isubcase, title, subtitle, label)
    return neids


def _spring_dx(element_name: str,
               neids: int,
               ieids: np.ndarray,
               eids: np.ndarray,
               dtype: str='float32'):

    elem = get_element(model, element_name, ieids, eids)
    dx = np.full((neids, 1), np.nan, dtype=fdtype)
    if element_name in {'CELAS1', 'CELAS2'}:
        for ieid, eid, nodes, comps in zip(ieids, eids, elem.nodes, elem.components):
            nid1, nid2 = elem.nodes
            c1, c2 = elem.c1, elem.c2
            i = dof_map[(nid1, c1)]
            j = dof_map[(nid2, c2)]
            dx[ieid] = xg[j] - xg[i]  # TODO: check the sign
    elif element_name in {'CELAS3', 'CELAS4'}:
        for ieid, eid, nodes in zip(ieids, eids, elem.spoints):
            nid1, nid2 = elem.nodes
            i = dof_map[(nid1, 0)]
            j = dof_map[(nid2, 0)]
            dx[ieid] = xg[j] - xg[i]  # TODO: check the sign
    else:  # pragma: no cover
        raise NotImplementedError(element_name)
    return elem, dx


def _save_spring_force(op2: OP2, f06_file: TextIO, page_num: int, page_stamp: str,
                       element_name: str,
                       force: np.ndarray, eids: np.ndarray,
                       write_f06_force: bool,
                       isubcase: int, title: str, subtitle: str,
                       label: str) -> None:
    if force is None:
        return
    data = force.reshape(1, *force.shape)
    table_name = 'OEF1'
    spring_force_obj = RealSpringForceArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    force = op2.op2_results.force
    obj = getattr(force, f'{element_name.lower()}_force')
    obj[isubcase] = spring_force_obj
    if write_f06_force:
        spring_force_obj.write_f06(
            f06_file, header=None, page_stamp=page_stamp,
            page_num=page_num, is_mag_phase=False, is_sort1=True)


def _save_spring_stress(stress_result,
                        f06_file, page_num, page_stamp,
                        element_name,
                        stress, eids, write_f06_stress: bool,
                        isubcase: int, title: str, subtitle: str,
                        label: str) -> None:
    data = stress.reshape(1, *stress.shape)
    table_name = 'OES1'
    spring_stress_obj = RealSpringStressArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    obj = getattr(stress_result, f'{element_name.lower()}_stress')
    obj[isubcase] = spring_stress_obj
    if write_f06_stress:
        spring_stress_obj.write_f06(
            f06_file, header=None, page_stamp=page_stamp,
            page_num=page_num, is_mag_phase=False, is_sort1=True)


def _save_spring_strain(op2, f06_file, page_num, page_stamp,
                        element_name,
                        strain, eids, write_f06_strain: bool,
                        isubcase: int, title: str, subtitle: str,
                        label: str) -> None:
    if strain is None:
        return
    data = strain.reshape(1, *strain.shape)
    table_name = 'OSTR1'
    spring_strain_obj = RealSpringStrainArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    obj = getattr(strain_result, f'{element_name.lower()}_strain')
    obj[isubcase] = spring_strain_obj
    if write_f06_strain:
        spring_strain_obj.write_f06(
            f06_file, header=None, page_stamp=page_stamp,
            page_num=page_num, is_mag_phase=False, is_sort1=True)

def _save_strain_energy(op2, f06_file, page_num, page_stamp,
                        element_name,
                        strain_energy, eids, write_f06_ese: bool,
                        isubcase: int, title: str, subtitle: str,
                        label: str) -> None:
    if strain_energy is None:
        return
    if strain_energy.sum() == 0.0:
        return
    data = strain_energy.reshape(1, *strain_energy.shape)
    table_name = 'ONRGY1'
    strain_energy_obj = RealStrainEnergyArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    ese = op2.op2_results.strain_energy
    obj = getattr(ese, f'{element_name.lower()}_strain_energy')
    obj[isubcase] = strain_energy_obj
    if write_f06_ese:
        strain_energy_obj.write_f06(
            f06_file, header=None, page_stamp=page_stamp,
            page_num=page_num, is_mag_phase=False, is_sort1=True)
