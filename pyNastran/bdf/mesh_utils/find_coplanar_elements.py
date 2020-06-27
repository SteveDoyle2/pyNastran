from __future__ import annotations
from typing import List, Union, Optional, TYPE_CHECKING

import numpy as np
from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def find_coplanar_triangles(bdf_filename: Union[BDF, str],
                            eids: Optional[List[int]]=None) -> List[int]:
    """
    Finds coplanar triangles

    Parameters
    ----------
    bdf_filename : BDF/str
        BDF: a model
        str: the path to the bdf input file
    eids : list
        the element ids to consider

    Returns
    -------
    coplanar_eids : List[int]
        the elements that are coplanar

    """
    model = get_bdf_model(bdf_filename, xref=False, log=None, debug=False)
    log = model.log

    if eids is None:
        eids = model.elements.keys()

    i = 0
    eids_removed = []
    neids = len(eids)
    nids = np.zeros((neids, 3), dtype='int32')
    for eid in eids:
        elem = model.elements[eid]
        try:
            nids[i, :] = elem.nodes
        except ValueError:
            eids_removed.append(eid)
            assert len(elem.nodes) != 3, str(elem)
            continue
        i += 1

    if i != neids:
        log.warning(f'removed {neids-i} non-triangles; eids_removed={eids_removed}')
        nids = nids[:i, :]

    #nids = np.array([
        #[10, 20, 30],
        #[20, 30, 10],
        #[10, 30, 20],
    #], dtype='int32')

    # [1, 2, 3]
    # [2, 3, 1]
    # [1, 3, 2]

    #imin = nids.argmin(axis=1)
    #imax = nids.argmax(axis=1)
    imin = nids.min(axis=1)
    imax = nids.max(axis=1)

    #print('imin = %s' % (imin))  # [0, 2, 0]
    #print('imax = %s' % (imax))  # [2, 1, 1]


    imid = []
    for row, imini, imaxi in zip(nids, imin, imax):
        #a = [imini, imaxi]
        #print(row, imini, imaxi)
        a = list(row)
        #a.remove(row[imini])
        #a.remove(row[imaxi])
        #print(a)
        a.remove(imini)
        #print(a)
        a.remove(imaxi)
        #print(a)
        #print('')
        imid.append(a[0])

    #print('imid = %s' % (imid))  # [1, 0, 2]

    nids2 = np.vstack([imin, imid, imax]).T
    aset = set()
    eids_to_remove = set()
    for eid, row in zip(eids, nids2):
        new_row = tuple(list(row))
        if new_row in aset:
            log.debug(f'eid={eid} exists already...')
            eids_to_remove.add(eid)
        else:
            aset.add(new_row)
    return model, eids_to_remove
