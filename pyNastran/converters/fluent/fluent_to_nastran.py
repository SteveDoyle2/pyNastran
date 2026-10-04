from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF


def fluent_to_nastran(model, bdf_filename: str) -> BDF:
    from pyNastran.bdf.bdf import BDF
    nastran_model = BDF()
    for nid, (x, y, z) in zip(model.node_id, model.xyz):
        nastran_model.add_grid(nid, [x, y, z])

    # row1 = [eid, dim, pid, 3, 4]
    # row2 = [eid, n1, n2, n3, n4]
    tris = model.tris
    quads = model.quads
    all_pids = []
    if len(tris) > 0:
        pids = tris[:, 1]
        all_pids.append(pids)
        for eid, pid, n1, n2, n3 in zip(tris[:, 0], pids,
                                        tris[:, 2], tris[:, 3], tris[:, 4]):
            nastran_model.add_ctria3(eid, pid, [n1, n2, n3])

    if len(quads) > 0:
        pids = tris[:, 1]
        all_pids.append(pids)
        for eid, pid, n1, n2, n3, n4 in zip(quads[:, 0], quads[:, 1],
                                            quads[:, 2], quads[:, 3], quads[:, 4], quads[:, 5]):
            nastran_model.add_cquad4(eid, pid, [n1, n2, n3, n4])
    upids = np.unique(np.hstack(all_pids))

    t = 0.1
    E = 3.0e7
    G = None
    nu = 0.3
    for pid in upids:
        mid = pid
        nastran_model.add_pshell(pid, mid, t)
        nastran_model.add_mat1(mid, E, G, nu)

    if bdf_filename:
        nastran_model.write_bdf(bdf_filename)
    return nastran_model
