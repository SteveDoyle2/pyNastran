import numpy as np
from pyNastran.bdf.cards.aero.utils import points_elements_from_quad_points

def create_structured_cquad4s(model, pid,
                              p1, p2, p3, p4, nx, ny, nid=1, eid=1):
    nid0 = nid
    x = np.linspace(0., 1., nx + 1)
    y = np.linspace(0., 1., ny + 1)
    points, elements = points_elements_from_quad_points(p1, p2, p3, p4, x, y, dtype='int32')
    for point in points:
        model.add_grid(nid, point)
        nid += 1

    for node_ids in elements + nid0:
        model.add_cquad4(eid, pid, node_ids)
        eid += 1
    return nid, eid
