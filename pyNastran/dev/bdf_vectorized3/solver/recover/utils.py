from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import Subcase


def get_plot_request(subcase: Subcase,
                     request: str) -> tuple[str, bool, bool, bool]:
    """
    request = 'SPCFORCES'
    """
    if request not in subcase:
        return '', False, False, True

    value, options = subcase.get_parameter(request)
    # write_f06 = False
    write_f06 = True
    if 'PRINT' in options:
        write_f06 = True
    if 'PLOT' in options:
        write_op2 = True
    if not(write_f06 or write_op2):
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
