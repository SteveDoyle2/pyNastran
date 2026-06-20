from __future__ import annotations
import warnings
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_classes import (
    RealStrainEnergyArray,)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import Subcase
    from pyNastran.op2.op2 import OP2


def fix_xb_shape(xb: np.ndarray):
    if xb.ndim == 1:
        ndof = xb.shape
        nmode = 1
    else:
        assert xb.shape == 2, xb.shape
        ndof, nmode = xb.shape
        assert nmode == 1, xb.shape
    return xb, nmode

def get_plot_request(subcase: Subcase,
                     request: str) -> tuple[str, bool, bool, bool]:
    """
    request = 'SPCFORCES'
    """
    if request not in subcase:
        return '', False, False, True

    value, options = subcase.get_parameter(request)
    write_f06 = True
    write_op2 = False
    if 'PRINT' in options:
        write_f06 = True
    if 'PLOT' in options:
        write_op2 = True
    if not (write_f06 or write_op2):
        write_op2 = True

    write_results = any((write_f06, write_op2))
    quick_return = not write_results
    nids_write = value
    return nids_write, write_f06, write_op2, quick_return


def get_f06_op2_pch_set(subcase: Subcase, key: str,
                        ) -> tuple[bool, bool, bool,
                                   list[str], np.ndarray]:
    """
    write_f06, write_op2, write_pch, options, set_data = get_f06_op2_pch_set(
        subcase, 'DISPLACEMENT')

    set_data = -1 -> don't calculate
    set_data = 0 -> all
    set_data > 0 -> subset
    """
    options = []
    write_pch = False
    write_f06 = False
    write_op2 = False
    set_data = np.array([-1], dtype='int32')
    if key in subcase:
        value, options = subcase[key]
        write_f06 = 'PRINT' in options
        write_op2 = 'PLOT' in options
        write_pch = 'PUNCH' in options
        if value == 'ALL':
            set_data = np.array([0], dtype='int32')
        else:
            set_id = value
            set_value, set_options = subcase[f'SET {set_id}']
            del set_options
            # set_data = set_value_to_set_data(set_value)
            set_data = np.array(set_value, dtype='int32')
            assert isinstance(set_data, np.ndarray), set_value
    return write_f06, write_op2, write_pch, options, set_data


def get_ieids_eids(model: BDF,
                   element_name: str,
                   eids_write: np.ndarray,
                   ) -> tuple[int, Any, np.ndarray, np.ndarray]:
    assert np.array_equal(eids_write, np.unique(eids_write))

    element_obj = getattr(model, element_name.lower())
    if eids_write[0] == 0:
        eids = element_obj.element_id
        assert np.array_equal(eids, np.unique(eids))
        ieids = None
    else:
        eids = np.intersect1d(element.element_id, eids_write)
        assert np.array_equal(eids, np.unique(eids))
        ieids = np.searchsorted(eids, eids_write)
    neid = len(eids)
    return neid, element_obj, ieids, eids


def get_mag_phase_from_options(options: list[str]) -> bool:
    is_mag_phase = False
    if 'PHASE' in options:
        is_mag_phase = True
    return is_mag_phase


def save_strain_energy(
    op2: OP2,
    f06_file, page_num, page_stamp,
    element_name: str,
    strain_energy: np.ndarray,
    eids: np.ndarray,
    isubcase: int,
    title: str, subtitle: str, label: str,
    write_f06: bool=True,
    write_op2: bool=True,
    case=None) -> None:
    """Save strain energy results to OP2 and write F06."""
    if strain_energy is None:
        warnings.warn(f'no strain energy for {element_name}')
    if strain_energy.sum() == 0.0:
        warnings.warn(f'empty strain energy for {element_name}')
        return
    
    if case is None:
        case = {}
     
    assert strain_energy.ndim == 3, strain_energy.shape
    data = strain_energy
    #data = strain_energy.reshape(1, *strain_energy.shape)
    assert np.all(np.isfinite(data)), data
    table_name = 'ONRGY1'
    #try:

    case_type = case.get('type', 'static')
    if case_type == 'static':
        se_obj = RealStrainEnergyArray.add_static_case(
            table_name, element_name, eids, data, isubcase,
            is_sort1=True, is_random=False, is_msc=True,
            random_code=0, title=title, subtitle=subtitle, label=label)
    elif case_type == 'modes':
        eigenvalue = case['eigenvalue']
        se_obj = RealStrainEnergyArray.add_modal_case(
            table_name, element_name, eids, data, isubcase,
            is_sort1=True, is_random=False, is_msc=True,
            random_code=0, title=title, subtitle=subtitle, label=label)
    else:
        raise RuntimeError(case)
    #except (KeyError, NotImplementedError, AssertionError):
    #    return

    ese = op2.op2_results.strain_energy
    obj_dict = getattr(ese, f'{element_name.lower()}_strain_energy')
    if write_op2:
        obj_dict[isubcase] = se_obj
    if write_f06:
        se_obj.write_f06(
            f06_file, header=None, page_stamp=page_stamp,
            page_num=page_num, is_mag_phase=False, is_sort1=True,
        )
