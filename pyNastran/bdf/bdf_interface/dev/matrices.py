# pylint: disable=C0103
"""
defines:
  - make_gpwg(Mgg, reference_point, xyz_cid0, log)

"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def get_Ajj(model: BDF, xyz=None):
    """not finished"""
    if xyz is None:
        xyz = {}
        for nid, node in model.nodes.items():
            xyz[nid] = node.get_position()
    for unused_caero_id, caero in model.caeros.items():
        unused_centroids = caero.get_centroids()

    for unused_spline_id, spline in model.splines.items():
        unused_spline_nodes = spline.spline_nodes

    Ajj = None
    return Ajj
