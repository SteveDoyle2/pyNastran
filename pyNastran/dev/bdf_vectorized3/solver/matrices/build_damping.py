import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF
import numpy as np
from scipy.sparse import dok_matrix

DOF_MAP = dict[tuple[int, int], int]


def build_Bbb_cdamp2(model: BDF,
                     Bbb: dok_matrix,
                     dof_map: DOF_MAP) -> int:
    """Fill CDAMP2 damping matrix."""
    elem = model.cdamp2
    if elem.n == 0:
        return 0
    nids1 = elem.nodes[:, 0]
    nids2 = elem.nodes[:, 1]
    c1s = elem.components[:, 0]
    c2s = elem.components[:, 1]
    bs = elem.b
    for nid1, nid2, c1, c2, bi in zip(nids1, nids2, c1s, c2s, bs):
        i = dof_map[(nid1, c1)]
        j = dof_map[(nid2, c2)]
        Bbb[i, i] += bi
        Bbb[j, j] += bi
        Bbb[i, j] -= bi
        Bbb[j, i] -= bi
    return elem.n
