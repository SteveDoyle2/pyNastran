"""
defines:
 - model = shift(bdf_filename, dxyz, bdf_filename_out=None)

"""
import numpy as np

from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model


def shift(bdf_filename, dxyz, bdf_filename_out=None):
    """shifts the model by some amount"""
    if isinstance(dxyz, list):
        dxyz = np.array(dxyz)
    assert isinstance(dxyz, np.ndarray), dxyz
    print("dxyz = %s" % dxyz)

    model = get_bdf_model(bdf_filename, xref=True, log=None, debug=True)
    for unused_nid, node in model.nodes.items():
        xyz = node.get_position() + dxyz
        node.set_position(model, xyz, cid=0, xref=True)

    for unused_caero_id, caero in model.caeros.items():
        caero.shift(dxyz)

    if bdf_filename_out:
        model.write_bdf(bdf_filename_out)
    return model

def update_nodes(model, nid_cp_cd, xyz_cid0):
    """how does this work for SPOINTs/EPOINTs???"""
    coord = model.coords[0]
    all_node_ids = np.array(list(model.nodes.keys()), dtype=nid_cp_cd.dtype)
    nids = nid_cp_cd[:, 0]
    inids = np.searchsorted(nids, all_node_ids)
    for inid, nid in zip(inids, all_node_ids):
        node = model.nodes[nid]
        xyz = xyz_cid0[inid, :]
        node.xyz = xyz
        node.cp = 0
        node.cp_ref = coord
