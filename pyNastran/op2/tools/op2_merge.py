import os
import sys
from pathlib import Path
from typing import Optional

from cpylog import SimpleLogger
from pyNastran.utils import PathLike
from pyNastran.op2.op2 import OP2
from pyNastran.op2.tools.solution_combination import _check_op2_file, _read_op2, _add_base_cases, _results

#--------------------------
# multi-combinations

def run_op2_merge(op2_filenames: list[PathLike | OP2],
                  op2_filename_new: Optional[PathLike]=None,
                  exclude_results: Optional[list[str]]=None,
                  include_results: Optional[list[str]]=None,
                  mode: Optional[str]=None,
                  log: Optional[SimpleLogger]=None) -> OP2:
    """
    Create a series of load case combinations based on a single op2

    Parameters
    ----------
    op2_filenames : list[PathLike]
        path to input file
    op2_filename_new : PathLike
        path to output file
    mode : Optional[str]
        the nastran format (msc, nx)
    log : Optional[SimpleLogger]
        a logger object

    """
    if op2_filename_new is None:
        base, ext = os.path.splitext(op2_filenames[0])
        op2_filename_new = base + '_combined.op2'

    # subcases: Optional[list[int]] = None,
    subcases = None
    #------------------------------------------------
    for op2_filename in op2_filenames:
        _check_op2_file(op2_filename)
    print(f'op2_filename_new = {op2_filename_new}')
    #------------------------------------------------
    # TODO: add a subcase check

    op2_filename = op2_filenames[0]
    model = _read_op2(
        op2_filename, mode=mode, log=log,
        exclude_results=exclude_results,
        include_results=include_results,
        subcases=subcases)

    # do the merge
    for i, op2_filenamei in enumerate(op2_filenames[1:]):
        modeli = _read_op2(
            op2_filenamei, mode=mode, log=log,
            exclude_results=exclude_results,
            include_results=include_results,
            subcases=subcases)

        table_res_types = _results(modeli)
        _add_cases(model, modeli, table_res_types)

    model.log.info(f'writing {os.path.abspath(op2_filename_new)}')
    model.write_op2(op2_filename_new)
    model_new = _read_op2(op2_filename_new, mode=None, log=log)
    return model

def _add_cases(model: OP2, model2: OP2, table_res_types) -> None:
    for table_type, res_type in table_res_types:
        slot = model2.get_result(table_type)
        inslot = model.get_result(table_type)
        if len(slot) == 0:
            # model2.log.warning(f'skipping slot = {table_type}')
            continue
        for key, value in slot.items():
            inslot[key] = value
    return

def cmd_line_merge(quiet: bool=False):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('op2_filenames', help='input op2s', nargs='+')
    parser.add_argument('-o', '--out', help='output op2s', default=None)

    args = parser.parse_args()
    if not quiet:  # pragma: no cover
        print(args)
    op2_filenames = args.op2_filenames
    op2_filename_new = args.out
    print(op2_filenames, op2_filename_new)
    run_op2_merge(op2_filenames, op2_filename_new=op2_filename_new)


def main():  # pragma: no cover
    dirname = Path(r'C:\solid_bending\merge')
    op2_filenames = [
        dirname / 'solid_case_axial.op2',
        dirname / 'solid_case_bending.op2',
    ]
    run_op2_merge(op2_filenames, op2_filename_new=None)


if __name__ == '__main__':  # pragma: no cover
    if len(sys.argv) > 1:
        cmd_line_merge()
    else:
        main()
