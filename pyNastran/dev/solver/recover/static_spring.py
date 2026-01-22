from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.dev.solver.utils import get_ieids_eids
from pyNastran.op2.op2_interface.op2_classes import (
    RealStrainEnergyArray,
    RealSpringStrainArray, RealSpringStressArray, RealSpringForceArray,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF # , Subcase

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
    _save_spring_strain_energy(
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


def get_celas1_ostr_ese(model, ielas, eids, xg, dof_map,
                        strain, unused_stress, unused_force, strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        straini, unused_forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        strain[ieid] = straini
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_ostr_oef(model, ielas, eids, xg, dof_map,
                        strain, unused_stress, force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        straini, forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        force[ieid] = forcei
        strain[ieid] = straini

def get_celas1_oes_ese(model, ielas, eids, xg, dof_map,
                       unused_strain, stress, unused_force, strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_ostr_oes_oef_ese(model, ielas, eids, xg, dof_map,
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

def get_celas1_ostr_oes_oef(model, ielas, eids, xg, dof_map,
                            strain, stress, force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        straini, forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        force[ieid] = forcei
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain[ieid] = straini

def get_celas1_ostr_oef_ese(model, ielas, eids, xg, dof_map,
                            strain, stress, force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        straini, forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)
        force[ieid] = forcei
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain[ieid] = straini

def get_celas1_ostr_oes(model, ielas, eids, xg, dof_map,
                        strain, stress, unused_force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain[ieid] = _recover_straini_celas12(xg, dof_map, elem, si)

def get_celas1_ostr(model, ielas, eids, xg, dof_map,
                    strain, unused_stress, unused_force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        si = elem.pid_ref.s
        strain[ieid] = _recover_straini_celas12(xg, dof_map, elem, si)

def get_celas1_oes_oef_ese(model, ielas, eids, xg, dof_map,
                           unused_strain, stress, force, strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        force[ieid] = _recover_forcei_celas12(xg, dof_map, elem, ki)
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_oes_oef(model, ielas, eids, xg, dof_map,
                       unused_strain, stress, force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        force[ieid] = _recover_forcei_celas12(xg, dof_map, elem, ki)
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)

def get_celas1_oes(model, ielas, eids, xg, dof_map,
                   unused_strain, stress, unused_force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)

def get_celas1_oef_ese(model, ielas, eids, xg, dof_map,
                       unused_strain, unused_stress, force, strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        si = elem.pid_ref.s
        force[ieid] = _recover_forcei_celas12(xg, dof_map, elem, ki)
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def get_celas1_oef(model, ielas, eids, xg, dof_map,
                   unused_strain, unused_stress, force, unused_strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        ki = elem.pid_ref.k
        force[ieid] = _recover_forcei_celas12(xg, dof_map, elem, ki)

def get_celas1_ese(model, ielas, eids, xg, dof_map,
                   unused_strain, unused_stress, unused_force, strain_energy):
    for ieid, eid in zip(ielas, eids):
        elem = model.elements[eid]
        si = elem.pid_ref.s
        strain_energy[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, si)

def _save_spring_stress(op2, f06_file, page_num, page_stamp,
                        element_name,
                        stress, eids, write_f06_stress: bool,
                        isubcase: int, title: str, subtitle: str, label: str) -> None:
    data = stress.reshape(1, *stress.shape)
    table_name = 'OES1'
    spring_stress = RealSpringStressArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    stress = op2.op2_results.stress
    if element_name == 'CELAS1':
        stress.celas1_stress[isubcase] = spring_stress
    elif element_name == 'CELAS2':
        stress.celas2_stress[isubcase] = spring_stress
    elif element_name == 'CELAS3':
        stress.celas3_stress[isubcase] = spring_stress
    elif element_name == 'CELAS4':
        stress.celas4_stress[isubcase] = spring_stress
    else:  # pragma: no cover
        raise NotImplementedError(element_name)
    if write_f06_stress:
        spring_stress.write_f06(f06_file, header=None, page_stamp=page_stamp,
                                page_num=page_num, is_mag_phase=False, is_sort1=True)

def _save_spring_strain(op2, f06_file, page_num, page_stamp,
                        element_name,
                        strain, eids, write_f06_strain: bool,
                        isubcase: int, title: str, subtitle: str, label: str) -> None:
    if strain is None:
        return
    data = strain.reshape(1, *strain.shape)
    table_name = 'OSTR1'
    spring_strain = RealSpringStrainArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    strain = op2.op2_results.strain
    if element_name == 'CELAS1':
        strain.celas1_strain[isubcase] = spring_strain
    elif element_name == 'CELAS2':
        strain.celas2_strain[isubcase] = spring_strain
    elif element_name == 'CELAS3':
        strain.celas3_strain[isubcase] = spring_strain
    elif element_name == 'CELAS4':
        strain.celas4_strain[isubcase] = spring_strain
    else:  # pragma: no cover
        raise NotImplementedError(element_name)
    if write_f06_strain:
        spring_strain.write_f06(f06_file, header=None, page_stamp=page_stamp,
                                page_num=page_num, is_mag_phase=False, is_sort1=True)

def _recover_stress_celas(f06_file, op2,
                          model: BDF, dof_map, isubcase, xg, eids_str,
                          element_name: str, fdtype='float32',
                          title: str='', subtitle: str='', label: str='',
                          page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static spring stress"""
    neids, ielas, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    write_f06_stress = True
    stress = np.full((neids, 1), np.nan, dtype=fdtype)
    if element_name == 'CELAS1':
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            pid_ref = elem.pid_ref
            ki = pid_ref.K()
            si = pid_ref.s
            stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
    elif element_name == 'CELAS2':
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            ki = elem.K()
            si = elem.s
            stress[ieid] = _recover_stressi_celas12(xg, dof_map, elem, ki, si)
    elif element_name == 'CELAS3':
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            pid_ref = elem.pid_ref
            ki = pid_ref.K()
            si = pid_ref.s
            stress[ieid] = _recover_stressi_celas34(xg, dof_map, elem, ki, si)
    elif element_name == 'CELAS4':
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            ki = elem.K()
            si = 1.0 # TODO: is this right?
            #si = elem.s
            stress[ieid] = _recover_stressi_celas34(xg, dof_map, elem, ki, si)
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    _save_spring_stress(
        op2, f06_file, page_num, page_stamp,
        element_name,
        stress, eids, write_f06_stress,
        isubcase, title, subtitle, label)
    return neids

def _recover_force_celas(f06_file, op2,
                         model: BDF, dof_map, isubcase, xg, eids_str,
                         element_name: str, fdtype='float32',
                         title: str='', subtitle: str='', label: str='',
                         page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static spring force"""
    neids, ielas, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids
    force = np.full((neids, 1), np.nan, dtype=fdtype)
    if element_name in {'CELAS1', 'CELAS2'}:
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            ki = elem.K()
            forcei = _recover_strain_forcei_celas12(xg, dof_map, elem, ki)[1:]
            force[ieid] = forcei
    elif element_name in {'CELAS3', 'CELAS4'}:
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            ki = elem.K()
            forcei = _recover_forcei_celas34(xg, dof_map, elem, ki)
            force[ieid] = forcei
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    write_f06_force = True
    _save_spring_force(
        op2, f06_file, page_num, page_stamp,
        element_name,
        force, eids, write_f06_force,
        isubcase, title, subtitle, label)
    return neids

def _recover_strain_energy_celas(f06_file, op2,
                                 model: BDF, dof_map, isubcase, xg, eids_str,
                                 element_name: str, fdtype='float32',
                                 title: str='', subtitle: str='', label: str='',
                                 page_num: int=1, page_stamp='PAGE %s') -> None:
    """recovers static spring strain energy"""
    neids, ielas, eids = get_ieids_eids(model, element_name, eids_str)
    if not neids:
        return neids

    write_f06_ese = True
    strain_energies = np.full((neids, 1), np.nan, dtype=fdtype)
    if element_name in {'CELAS1', 'CELAS2'}:
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            ki = elem.K()
            strain_energies[ieid] = _recover_strain_energyi_celas12(xg, dof_map, elem, ki)
    elif element_name in {'CELAS3', 'CELAS4'}:
        for ieid, eid in zip(ielas, eids):
            elem = model.elements[eid]
            ki = elem.K()
            strain_energies[ieid] = _recover_strain_energyi_celas34(xg, dof_map, elem, ki)
    else:  # pragma: no cover
        raise NotImplementedError(element_name)

    _save_spring_strain_energy(
        op2, f06_file, page_num, page_stamp,
        element_name,
        strain_energies, eids, write_f06_ese,
        isubcase, title, subtitle, label)
    return neids


def _recover_strain_forcei_celas12(xg, dof_map, elem, ki: float):
    """get the static spring force"""
    # F = kx
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    strain = xg[j] - xg[i]  # TODO: check the sign
    force = ki * strain
    return strain, force

def _recover_stressi_celas12(xg, dof_map, elem, ki: float, si: float):
    """get the static spring stress"""
    # F = kx
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    strain = xg[j] - xg[i]  # TODO: check the sign
    stress = ki * si * strain
    return stress

def _recover_stressi_celas34(xg, dof_map, elem, ki: float, si: float):
    """get the static spring stress"""
    # F = kx
    nid1, nid2 = elem.nodes
    i = dof_map[(nid1, 0)]
    j = dof_map[(nid2, 0)]
    strain = xg[j] - xg[i]  # TODO: check the sign
    stress = ki * si * strain
    return stress

def _recover_straini_celas12(xg, dof_map, elem, si: float):
    """get the static spring strain"""
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    strain = xg[j] - xg[i]  # TODO: check the sign
    # TODO: why is the strain 0?
    return si * strain

def _recover_straini_celas34(xg, dof_map, elem, si: float):
    """get the static spring strain"""
    nid1, nid2 = elem.nodes
    i = dof_map[(nid1, 0)]
    j = dof_map[(nid2, 0)]
    strain = xg[j] - xg[i]  # TODO: check the sign
    # TODO: why is the strain 0?
    return si * strain

def _recover_forcei_celas12(xg, dof_map, elem, ki: float):
    """get the static spring force"""
    # F = kx
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    strain = xg[j] - xg[i]  # TODO: check the sign
    force = ki * strain
    return force

def _recover_forcei_celas34(xg, dof_map, elem, ki: float):
    """get the static spring force"""
    # F = kx
    nid1, nid2 = elem.nodes
    i = dof_map[(nid1, 0)]
    j = dof_map[(nid2, 0)]
    strain = xg[j] - xg[i]  # TODO: check the sign
    force = ki * strain
    return force

def _recover_strain_energyi_celas12(xg, dof_map, elem, ki: float):
    """get the static spring force"""
    # F = kx
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    strain = xg[j] - xg[i]  # the sign doesn't matter
    force = ki * strain ** 2
    return force


def _recover_strain_energyi_celas34(xg, dof_map, elem, ki: float):
    """get the static spring force"""
    # F = kx
    nid1, nid2 = elem.nodes
    i = dof_map[(nid1, 0)]
    j = dof_map[(nid2, 0)]
    strain = xg[j] - xg[i]  # the sign doesn't matter
    force = ki * strain ** 2
    return force


def _save_spring_force(op2, f06_file, page_num, page_stamp,
                       element_name,
                       force, eids, write_f06_force: bool,
                       isubcase: int, title: str, subtitle: str, label: str) -> None:
    if force is None:
        return
    data = force.reshape(1, *force.shape)
    table_name = 'OEF1'
    spring_force = RealSpringForceArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    force = op2.op2_results.force
    if element_name == 'CELAS1':
        force.celas1_force[isubcase] = spring_force
    elif element_name == 'CELAS2':
        force.celas2_force[isubcase] = spring_force
    elif element_name == 'CELAS3':
        force.celas3_force[isubcase] = spring_force
    elif element_name == 'CELAS4':
        force.celas4_force[isubcase] = spring_force
    else:  # pragma: no cover
        raise NotImplementedError(element_name)
    if write_f06_force:
        spring_force.write_f06(f06_file, header=None, page_stamp=page_stamp,
                               page_num=page_num, is_mag_phase=False, is_sort1=True)

def _save_spring_strain_energy(op2, f06_file, page_num, page_stamp,
                               element_name,
                               strain_energy, eids, write_f06_ese: bool,
                               isubcase: int, title: str, subtitle: str, label: str) -> None:
    if strain_energy is None:
        return
    data = strain_energy.reshape(1, *strain_energy.shape)
    table_name = 'ONRGY1'
    spring_strain_energy = RealStrainEnergyArray.add_static_case(
        table_name, element_name, eids, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)

    strain_energy = op2.op2_results.strain_energy
    if element_name == 'CELAS1':
        strain_energy.celas1_strain_energy[isubcase] = spring_strain_energy
    elif element_name == 'CELAS2':
        strain_energy.celas2_strain_energy[isubcase] = spring_strain_energy
    elif element_name == 'CELAS3':
        strain_energy.celas3_strain_energy[isubcase] = spring_strain_energy
    elif element_name == 'CELAS4':
        strain_energy.celas4_strain_energy[isubcase] = spring_strain_energy
    else:  # pragma: no cover
        raise NotImplementedError(element_name)
    if write_f06_ese:
        spring_strain_energy.write_f06(f06_file, header=None, page_stamp=page_stamp,
                                       page_num=page_num, is_mag_phase=False, is_sort1=True)
