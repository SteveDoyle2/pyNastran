from __future__ import annotations
from typing import TextIO, TYPE_CHECKING
import numpy as np

from pyNastran.dev.bdf_vectorized3.solver.utils import (
    get_ieids_eids, get_element, lambda1d)
from pyNastran.op2.op2_interface.op2_classes import (
    RealStrainEnergyArray,)
from .static_spring import _recover_strain_energy_celas
from .rod import _recover_strain_energy_rod
from .beam import _recover_strain_energy_beam
from .utils import get_plot_request, save_strain_energy, fix_xb_shape


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
    DOF_MAP = dict[tuple[int, int], int]


def recover_strain_energy_101(f06_file: TextIO, op2,
                              model: BDF, dof_map,
                              subcase: Subcase, xb, fdtype: str='float32',
                              title: str='', subtitle: str='', label: str='',
                              page_num: int=1, page_stamp: str='PAGE %s'):
    """Recovers element strain energy from ESE = ALL."""
    if 'ESE' not in subcase:
        return
    eid_str = 'ALL'
    unused_eids_write, write_f06, write_op2, quick_return = get_plot_request(
        subcase, 'ESE')
    if quick_return:
        return page_num
    isubcase = subcase.id

    nelements = 0
    for name in ['CELAS1', 'CELAS2', 'CELAS3', 'CELAS4']:
        nelements += _recover_strain_energy_celas(
            f06_file, op2, model, dof_map, isubcase, xb, eid_str,
            name, fdtype=fdtype,
            title=title, subtitle=subtitle, label=label,
            page_num=page_num, page_stamp=page_stamp)

    for name in ['CROD', 'CTUBE', 'CONROD']:
        nelements += _recover_strain_energy_rod(
            f06_file, op2, model, dof_map, isubcase, xb, eid_str,
            name, fdtype=fdtype,
            title=title, subtitle=subtitle, label=label,
            page_num=page_num, page_stamp=page_stamp)

    nelements += _recover_strain_energy_beam(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CBAR', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_energy_beam(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CBEAM', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    if nelements == 0:
        model.log.warning(f'no strain energy output...{model.card_count}; {model.bdf_filename}')
