from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import Subcase

def get_plot_request(subcase: Subcase, request: str) -> tuple[str, bool, bool, bool]:
    """
    request = 'SPCFORCES'
    """
    if request not in subcase:
        return '', False, False, True

    value, options = subcase.get_parameter(request)
    write_f06 = False
    write_f06 = True
    if 'PRINT' in options:
        write_f06 = True
    if 'PLOT' in options:
        write_op2 = True
    if not(write_f06 or write_op2):
        write_op2 = True
    quick_return = not write_f06 and not write_op2
    nids_write = value
    return nids_write, write_f06, write_op2, quick_return
