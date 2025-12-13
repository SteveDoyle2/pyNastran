import os
from typing import Optional, Any
import numpy as np
from pyNastran.utils import print_bad_path


def load_ints_from_defaults(ids_dict: dict[int, Any],
                            infilename: Optional[str]) -> list[int]:
    if infilename is None:
        eids_to_fix = list(ids_dict)
    else:
        assert os.path.exists(infilename), print_bad_path(infilename)
        eids_to_fix = np.loadtxt(infilename, dtype='int32').flatten().tolist()
    return eids_to_fix
