import os
from typing import Optional

from pyNastran.bdf.bdf import BDF, read_bdf, print_card_8
from pyNastran.utils import PathLike
import numpy as np

def solid_dof(bdf_filename: BDF | PathLike,
              nid_filename: Optional[PathLike]=None,
              spc_id: int=100) -> tuple[BDF, np.ndarray]:
    """
    Get all solid nodes not associated with other elements

    Parameters
    ----------
    bdf_filename
    nid_filename

    Returns
    -------

    """
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        base, ext = os.path.splitext(bdf_filename)
        nid_filename = base + '.solid_dof_constraint.blk'
        model = BDF()
        model.is_strict_card_parser = False
        model.read_bdf(bdf_filename, xref=False)
    #nids = np.array(list(model.nodes), dtype='int32')

    solid_nids_set = set([])
    associated_nids_set = set([])

    # rigids
    for eid, elem in model.rigid_elements.items():
        nidsi = elem.nodes
        associated_nids_set.update(nidsi)

    # masses
    for eid, elem in model.masses.items():
        if hasattr(elem, 'nodes'):
            nidsi = elem.nodes
        else:
            nidsi = elem.node_ids
        associated_nids_set.update(nidsi)

    for eid, elem in model.elements.items():
        nidsi = elem.nodes
        if elem.type in {'CTETRA', 'CPENTA', 'CPYRAM', 'CHEXA'}:
            solid_nids_set.update(nidsi)
            continue
        associated_nids_set.update(nidsi)

    solid_nids_list = list(solid_nids_set)
    associated_nids_list = list(associated_nids_set)
    if None in associated_nids_list:
        associated_nids_list.remove(None)
    if None in solid_nids_list:
        solid_nids_list.remove(None)

    associated_nids = np.array(associated_nids_list, dtype='int32')
    solid_nids = np.array(solid_nids_list, dtype='int32')

    out_nids = np.setdiff1d(solid_nids, associated_nids)  # A - B
    card = ['SPC1', spc_id, 456] + out_nids.tolist()

    out_nids_str = ' '.join([str(val) for val in out_nids])
    if nid_filename is not None:
        with open(nid_filename, 'w') as nid_file:
            nid_file.write(f'$ add a 456 constraint onto every solid DOF not associated with another element\n')
            nid_file.write(f'$ out_nids = {out_nids_str}\n')
            nid_file.write(f'$ >>> bdf solid_dof {bdf_filename}\n')
            nid_file.write(print_card_8(card))
    return model, out_nids
