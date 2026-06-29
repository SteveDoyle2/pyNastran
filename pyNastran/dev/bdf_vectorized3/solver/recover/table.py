from typing import TextIO
import numpy as np
#-------------------------------
import os
import copy
from datetime import date
from typing import TextIO, Any

import numpy as np
#import scipy
#from scipy.linalg import eigh, eig

from cpylog import SimpleLogger

import pyNastran
from pyNastran.nptyping_interface import (
    NDArrayNbool,
    NDArrayNint,
    NDArrayN2int,
    NDArrayNfloat,
    NDArrayNNfloat,
)
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
from pyNastran.op2.op2 import OP2
from pyNastran.op2.op2_interface.op2_classes import (
    RealDisplacementArray,
    RealSPCForcesArray,
    RealLoadVectorArray,
    RealEigenvectorArray,
    RealMPCForcesArray,
    RealGridPointForcesArray,
)
from ..utils import recast_data, get_param, DOF_MAP
from ..utils_statics import save_static_table
from .utils import (
    get_f06_op2_pch_set, get_plot_request)

#-------------------------------
def _save_displacment(
    op2: OP2,
    f06_file: TextIO,
    subcase: Subcase,
    itime: int,
    ntimes: int,
    node_gridtype: NDArrayN2int,
    xg: NDArrayNfloat,
    ngrid: int,
    ndof_per_grid: int,
    title: str = "",
    subtitle: str = "",
    label: str = "",
    fdtype: str = "float32",
    page_num: int = 1,
    page_stamp: str = "PAGE %s",) -> int:
    f06_request_name = "DISPLACEMENT"
    table_name = "OUGV1"
    # self.log.debug(f'xg = {xg}')
    assert page_stamp is not None
    page_stamp % page_num
    page_num = save_static_table(
        f06_file,
        subcase,
        itime,
        ntimes,
        node_gridtype,
        xg,
        RealDisplacementArray,
        f06_request_name,
        table_name,
        op2.displacements,
        ngrid,
        ndof_per_grid,
        title=title,
        subtitle=subtitle,
        label=label,
        fdtype=fdtype,
        page_num=page_num,
        page_stamp=page_stamp,
    )
    return page_num


def _save_spc_forces(
    op2: OP2,
    f06_file: TextIO,
    subcase: Subcase,
    itime: int,
    ntimes: int,
    node_gridtype: NDArrayN2int,
    fspc: NDArrayNfloat,
    ngrid: int,
    ndof_per_grid: int,
    title: str = "",
    subtitle: str = "",
    label: str = "",
    fdtype: str = "float32",
    page_num: int = 1,
    page_stamp: str = "PAGE %s",) -> int:
    f06_request_name = "SPCFORCES"
    table_name = "OQG1"
    # self.log.debug(f'Fg = {Fg}')
    assert page_stamp is not None
    page_stamp % page_num
    page_num = save_static_table(
        f06_file,
        subcase,
        itime,
        ntimes,
        node_gridtype,
        fspc,
        RealSPCForcesArray,
        f06_request_name,
        table_name,
        op2.spc_forces,
        ngrid,
        ndof_per_grid,
        title=title,
        subtitle=subtitle,
        label=label,
        fdtype=fdtype,
        page_num=page_num,
        page_stamp=page_stamp,
    )
    return page_num


def _save_mpc_forces(
    op2: OP2,
    f06_file: TextIO,
    subcase: Subcase,
    itime: int,
    ntimes: int,
    node_gridtype: NDArrayN2int,
    fmpc: NDArrayNfloat,
    ngrid: int,
    ndof_per_grid: int,
    title: str = "",
    subtitle: str = "",
    label: str = "",
    fdtype: str = "float32",
    page_num: int = 1,
    page_stamp: str = "PAGE %s",) -> int:
    """Save MPC forces to F06 and OP2."""

    #page_num = save_static_table(
    #    f06_file,
    #    subcase,
    #    itime,
    #    ntimes,
    #    node_gridtype,
    #    fmpc,
    #    RealMPCForcesArray,
    #    f06_request_name,
    #    table_name,
    #    op2.mpc_forces,
    #    ngrid,
    #    ndof_per_grid,
    #    title=title,
    #    subtitle=subtitle,
    #    label=label,
    #    fdtype=fdtype,
    #    page_num=page_num,
    #    page_stamp=page_stamp,
    #)

    assert page_stamp is not None
    page_stamp % page_num
    f06_request_name = 'MPCFORCES'
    unused_nids_write, write_f06, write_op2, quick_return = get_plot_request(
        subcase, f06_request_name)
    if quick_return:
        return page_num

    idtype2, fdtype2 = recast_data("int32", fdtype)
    isubcase = subcase.id
    nnodes = node_gridtype.shape[0]
    data = np.zeros((ntimes, nnodes, 6), dtype=fdtype2)
    ngrid_dofs = ngrid * ndof_per_grid
    _fgi = fmpc[:ngrid_dofs].reshape(ngrid, ndof_per_grid)
    data[itime, :ngrid, :] = _fgi
    data[itime, ngrid:, 0] = fmpc[ngrid_dofs:]

    # Use OQG1 for data_code creation, then override table_name for write_f06
    table_name = "OQG1"
    mpc_obj = RealMPCForcesArray.add_static_case(
        table_name, node_gridtype, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)
    mpc_obj.table_name = "OQMG1"

    if write_f06:
        page_num = mpc_obj.write_f06(
            f06_file, header=None,
            page_stamp=page_stamp, page_num=page_num,
            is_mag_phase=False, is_sort1=True)
        f06_file.write("\n")
    if write_op2:
        op2.mpc_forces[isubcase] = mpc_obj
    return page_num


def _save_applied_load(
    op2: OP2,
    f06_file: TextIO,
    subcase: Subcase,
    itime: int,
    ntimes: int,
    node_gridtype: NDArrayN2int,
    Fg: NDArrayNfloat,
    ngrid: int,
    ndof_per_grid: int,
    title: str = "",
    subtitle: str = "",
    label: str = "",
    fdtype: str = "float32",
    page_num: int = 1,
    page_stamp: str = "PAGE %s",) -> int:
    f06_request_name = "OLOAD"
    table_name = "OPG1"
    # self.log.debug(f'Fg = {Fg}')
    assert page_stamp is not None
    page_stamp % page_num
    page_num = save_static_table(
        f06_file,
        subcase,
        itime,
        ntimes,
        node_gridtype,
        Fg,
        RealLoadVectorArray,
        f06_request_name,
        table_name,
        op2.load_vectors,
        ngrid,
        ndof_per_grid,
        title=title,
        subtitle=subtitle,
        label=label,
        fdtype=fdtype,
        page_num=page_num,
        page_stamp=page_stamp,
    )
    return page_num
