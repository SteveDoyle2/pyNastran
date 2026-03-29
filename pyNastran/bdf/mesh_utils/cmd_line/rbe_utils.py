from __future__ import annotations
import os
from typing import Optional, Any
import numpy as np
from pyNastran.utils import PathLike, print_bad_path
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF


def load_infile_eids_default(model: BDF,
                             eids_str: Optional[str],
                             infile: Optional[PathLike]) -> np.ndarray:
    log = model.log
    eids_to_fix = np.array([], dtype='int')
    if infile is None and eids_str is not None:
        eids_to_fix = np.array([int(eid) for eid in eids_str.strip(', ').split(',')])
        return eids_to_fix

    if infile:
        model.log.debug(f'infile = {infile!r}')
        eids_to_fix = load_ints_from_defaults(
            model.rigid_elements, infile)

    log.info(f'eids_to_fix = {eids_to_fix}')
    if len(eids_to_fix) == 0:
        eids_to_fix = list(model.rigid_elements)
        if len(eids_to_fix):
            log.warning(f'no eids were specified assuming all rigid elements={eids_to_fix}')
        else:
            # with open(model.bdf_filename, 'r') as bdf_file:
            #     print(bdf_file.read())
            print(model.get_bdf_stats())
            raise RuntimeError('no eids were specified and none were found...')
        log.info(f'eids_to_fix = {eids_to_fix}')
    return np.array(eids_to_fix, dtype='int32')

def load_ints_from_defaults(ids_dict: dict[int, Any],
                            infilename: Optional[str]) -> list[int]:
    if infilename is None:
        eids_to_fix = list(ids_dict)
    else:
        assert os.path.exists(infilename), print_bad_path(infilename)
        eids_to_fix = np.loadtxt(infilename, dtype='int32').flatten().tolist()
    return eids_to_fix
