from __future__ import annotations
from typing import TextIO, TYPE_CHECKING
import numpy as np

from pyNastran.dev.solver.utils import lambda1d, get_ieids_eids
from pyNastran.op2.op2_interface.op2_classes import (
    #RealSpringForceArray, RealRodForceArray, RealCBarForceArray,
    RealStrainEnergyArray,
)
from pyNastran.dev.solver.build_stiffness import ke_cbar
from .static_spring import _recover_strain_energy_celas
from .utils import get_plot_request

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, Subcase, CBAR, PBAR, PBARL


def recover_strain_energy_101(f06_file: TextIO, op2,
                              model: BDF, dof_map,
                              subcase: Subcase, xb, fdtype: str='float32',
                              title: str='', subtitle: str='', label: str='',
                              page_num: int=1, page_stamp: str='PAGE %s'):
    """
    recovers the forces from:
     - ESE = ALL

    """
    if 'ESE' not in subcase:
        return
    eid_str = 'ALL'
    unused_eids_write, write_f06, write_op2, quick_return = get_plot_request(
        subcase, 'ESE')
    if quick_return:
        return page_num
    isubcase = subcase.id

    nelements = 0
    nelements += _recover_strain_energy_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS1', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_energy_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS2', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_energy_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS3', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)
    nelements += _recover_strain_energy_celas(
        f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        'CELAS4', fdtype=fdtype,
        title=title, subtitle=subtitle, label=label,
        page_num=page_num, page_stamp=page_stamp)

    #nelements += _recover_force_rod(
        #f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        #'CROD', fdtype=fdtype,
        #title=title, subtitle=subtitle, label=label,
        #page_num=page_num, page_stamp=page_stamp)
    #nelements += _recover_force_rod(
        #f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        #'CONROD', fdtype=fdtype,
        #title=title, subtitle=subtitle, label=label,
        #page_num=page_num, page_stamp=page_stamp)
    #nelements += _recover_force_rod(
        #f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        #'CTUBE', fdtype=fdtype,
        #title=title, subtitle=subtitle, label=label,
        #page_num=page_num, page_stamp=page_stamp)
    #nelements += _recover_force_cbar(
        #f06_file, op2, model, dof_map, isubcase, xb, eid_str,
        #'CBAR', fdtype=fdtype,
        #title=title, subtitle=subtitle, label=label,
        #page_num=page_num, page_stamp=page_stamp)
    if nelements == 0:
        model.log.warning(f'no strain energy output...{model.card_count}; {model.bdf_filename}')
