from typing import TextIO, Any

import numpy as np
from pyNastran.nptyping_interface import (
    NDArrayN2int, NDArrayNfloat)
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
from pyNastran.op2.op2_interface.op2_classes import (
    RealDisplacementArray, RealSPCForcesArray,
    RealLoadVectorArray,)
from .recover.utils import get_plot_request
from .utils import recast_data

def save_static_table(f06_file: TextIO,
                      subcase: Subcase, itime: int, ntimes: int,
                      node_gridtype: NDArrayN2int, Fg: NDArrayNfloat,
                      obj: RealSPCForcesArray,
                      f06_request_name: str,
                      table_name: str, slot: dict[Any, RealSPCForcesArray],
                      ngrid: int, ndof_per_grid: int,
                      title: str = '', subtitle: str = '', label: str = '',
                      idtype: str = 'int32', fdtype: str = 'float32',
                      page_num: int = 1, page_stamp: str = 'PAGE %s') -> int:
    idtype, fdtype = recast_data(idtype, fdtype)
    isubcase = subcase.id
    # self.log.debug(f'saving {f06_request_name} -> {table_name}')
    unused_nids_write, write_f06, write_op2, quick_return = get_plot_request(
        subcase, f06_request_name)
    if quick_return:
        return page_num
    nnodes = node_gridtype.shape[0]
    data = np.zeros((ntimes, nnodes, 6), dtype=fdtype)
    ngrid_dofs = ngrid * ndof_per_grid
    if ndof_per_grid == 6:
        _fgi = Fg[:ngrid_dofs].reshape(ngrid, ndof_per_grid)
        data[itime, :ngrid, :] = _fgi
    else:
        raise NotImplementedError(ndof_per_grid)
    data[itime, ngrid:, 0] = Fg[ngrid_dofs:]

    spc_forces = obj.add_static_case(
        table_name, node_gridtype, data, isubcase,
        is_sort1=True, is_random=False, is_msc=True,
        random_code=0, title=title, subtitle=subtitle, label=label)
    if write_f06:
        page_num = spc_forces.write_f06(
            f06_file, header=None,
            page_stamp=page_stamp, page_num=page_num,
            is_mag_phase=False, is_sort1=True)
        f06_file.write('\n')
    if write_op2:
        slot[isubcase] = spc_forces
    return page_num
