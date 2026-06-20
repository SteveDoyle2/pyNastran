from __future__ import annotations
from typing import TextIO, TYPE_CHECKING
import numpy as np

from pyNastran.dev.bdf_vectorized3.solver.utils import (
    get_ieids_eids, get_element)
from .utils import get_plot_request
from .rod import _recover_strain_rod
from .bar import _recover_strain_cbar
from .beam import _recover_strain_cbeam

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
    DOF_MAP = dict[tuple[int, int], int]


def recover_strain_101(
    f06_file: TextIO,
    op2,
    model: BDF,
    dof_map: DOF_MAP,
    subcase: Subcase,
    xb: np.ndarray,
    fdtype: str = "float32",
    title: str = "",
    subtitle: str = "",
    label: str = "",
    page_num: int = 1,
    page_stamp: str = "PAGE %s",
):
    """Recovers element strains from STRAIN = ALL."""
    if "STRAIN" not in subcase:
        return
    eid_str = "ALL"
    unused_eids_write, write_f06, write_op2, quick_return = get_plot_request(
        subcase, "STRAIN")
    if quick_return:
        return page_num
    isubcase = subcase.id

    nelements = 0
    for name in ['CROD', 'CTUBE', 'CONROD']:
        nelements += _recover_strain_rod(
            f06_file,
            op2,
            model,
            dof_map,
            isubcase,
            xb,
            eid_str,
            name,
            fdtype=fdtype,
            title=title,
            subtitle=subtitle,
            label=label,
            page_num=page_num,
            page_stamp=page_stamp,)

    nelements += _recover_strain_cbar(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        eid_str,
        "CBAR",
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,)
    nelements += _recover_strain_cbeam(
        f06_file,
        op2,
        model,
        dof_map,
        isubcase,
        xb,
        eid_str,
        "CBEAM",
        fdtype=fdtype,
        title=title,
        subtitle=subtitle,
        label=label,
        page_num=page_num,
        page_stamp=page_stamp,)

    if nelements == 0:
        model.log.warning(f"no strain output...{model.card_count}; {model.bdf_filename}")
