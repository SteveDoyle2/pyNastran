from typing import TextIO
import numpy as np
from pyNastran.bdf.bdf import BDF, Subcase
from .utils import (
    get_f06_op2_pch_set, get_mag_phase_from_options,
    get_ieids_eids)
from pyNastran.op2.tables.oef_forces.oef_force_objects import RealSpringForceArray
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects  import ComplexSpringForceArray


def recover_force_freq(f06_file: TextIO, op2,
                       model: BDF, dof_map,
                       subcase: Subcase,
                       xq: np.ndarray,
                       # phig: np.ndarray,
                       # phib: np.ndarray,
                       # modes, eigns,
                       freqs: np.ndarray,
                       cfdtype: str='complex64',
                       title: str='', subtitle: str='', label: str='',
                       page_num: int=1, page_stamp: str='PAGE %s'):
    """
    recovers the forces from:
     - FORCE = ALL

    """
    # omegas = np.sqrt(np.abs(eigns))
    # freqs = omegas / (2 * np.pi)

    log = model.log
    if 'FORCE' not in subcase:
        return
    write_f06, write_op2, write_pch, options, eids_write = get_f06_op2_pch_set(
        subcase, 'FORCE')
    write_data = np.any((write_f06, write_op2))
    if not write_data:
        log.warning(f'skipping modal force')
        return page_num

    isubcase = subcase.id

    ngrid = len(model.grid)
    nspoint = len(model.spoint)
    # assert nspoint == 0, nspoint
    # nmode = xq.shape[0]
    # xg = phig.reshape(nmode, ngrid, 6)
    nelements = 0
    nelements, page_num = _recover_force_celas(
        f06_file, op2, model, isubcase,
        xq, eids_write, options,
        freqs,
        'CELAS1', cfdtype=cfdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp,
        nelements=nelements,
        write_f06=write_f06, write_op2=write_op2)
    nelements, page_num = _recover_force_celas(
        f06_file, op2, model, isubcase,
        xq, eids_write, options,
        freqs,
        'CELAS2', cfdtype=cfdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp,
        nelements=nelements,
        write_f06=write_f06, write_op2=write_op2)
    nelements, page_num = _recover_force_celas(
        f06_file, op2, model, isubcase,
        xq, eids_write, options,
        freqs,
        'CELAS3', cfdtype=cfdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp,
        nelements=nelements,
        write_f06=write_f06, write_op2=write_op2)
    nelements, page_num = _recover_force_celas(
        f06_file, op2, model, isubcase,
        xq, eids_write, options,
        freqs,
        'CELAS4', cfdtype=cfdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp,
        nelements=nelements,
        write_f06=write_f06, write_op2=write_op2)

    # nelements += _recover_force_rod(
    #     f06_file, op2, model, dof_map, isubcase, xb, eids_write,
    #     'CROD', fdtype=fdtype,
    #     title=title, subtitle=subtitle, label=label,
    #     page_num=page_num, page_stamp=page_stamp)
    # nelements += _recover_force_rod(
    #     f06_file, op2, model, dof_map, isubcase, xb, eids_write,
    #     'CONROD', fdtype=fdtype,
    #     title=title, subtitle=subtitle, label=label,
    #     page_num=page_num, page_stamp=page_stamp)
    # nelements += _recover_force_rod(
    #     f06_file, op2, model, dof_map, isubcase, xb, eids_write,
    #     'CTUBE', fdtype=fdtype,
    #     title=title, subtitle=subtitle, label=label,
    #     page_num=page_num, page_stamp=page_stamp)
    # nelements += _recover_force_cbar(
    #     f06_file, op2, model, dof_map, isubcase, xb, eids_write,
    #     'CBAR', fdtype=fdtype,
    #     title=title, subtitle=subtitle, label=label,
    #     page_num=page_num, page_stamp=page_stamp)
    if nelements == 0:
        log.warning(f'no force output...{model.card_count}; {model.bdf_filename}')
        asdf

def _recover_force_celas(f06_file: TextIO, op2,
                         model: BDF,
                         isubcase: int,
                         xq: np.ndarray,
                         eids_write: np.ndarray,  options: list[str],
                         freqs: np.ndarray,
                         element_name: str, cfdtype='complex64',
                         title: str='', subtitle: str='', label: str='',
                         page_num: int=1, page_stamp='PAGE %s',
                         nelements: int=0,
                         write_f06: bool=False, write_op2: bool=False,
                         ) -> tuple[int, int]:
    """recovers static spring force"""
    neid, elem, ieid, element_id = get_ieids_eids(model, element_name, eids_write)
    nelements += neid
    if not neid:
        return nelements, page_num

    forces = op2.op2_results.force
    slot = getattr(forces, element_name.lower() + '_force')
    modal_obj: RealSpringForceArray = slot[isubcase]

    nfreq = len(freqs)
    # nmode = xq.shape[0]

    # (2, 2, 1)(2052, 2)
    print(modal_obj.data.shape, xq.shape)
    # force = xq @ modal_obj.data
    # f: frequency
    # h: modal
    # o: one
    assert xq.ndim == 2, xq.shape
    assert modal_obj.data.ndim == 3, modal_obj.data.shape
    force = np.einsum('fh,heo->feo', xq, modal_obj.data)
    assert force.shape == (nfreq, neid, 1), force.shape

    # nmode = xq.shape[0]
    table_name = 'OEF1'
    data = force.reshape(nfreq, neid, 1)
    assert data.shape == (nfreq, neid, 1), data.shape

    # table_name: str,
    # element: np.ndarray,
    # data: np.ndarray,
    # isubcase: int,
    # freqs: np.ndarray,
    # element_name: str,
    obj = ComplexSpringForceArray.add_freq_case(
        table_name, element_name,
        element_id, data, isubcase,
        freqs,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)
    slot[isubcase] = obj

    header = ['', '', '']
    if write_f06:
        is_mag_phase = get_mag_phase_from_options(options)
        page_num = obj.write_f06(
            f06_file, header=header,
            page_stamp=page_stamp, page_num=page_num,
            is_mag_phase=is_mag_phase, is_sort1=True)
        f06_file.write('\n')
    else:
        raise RuntimeError('write the f06')

    # if write_op2:
    #     asdf
    # write_f06_force = True
    # _save_spring_force(
    #     op2, f06_file, page_num, page_stamp,
    #     element_name,
    #     force, eids, write_f06_force,
    #     isubcase, title, subtitle, label)
    return nelements, page_num
