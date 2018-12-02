from __future__ import print_function
import numpy as np
from pyNastran.bdf.bdf import read_bdf
from pyNastran.bdf.utils import parse_patran_syntax_dict


def find_coplanar_triangles(bdf_filename, eids):
    """

    finds coplanar triangles

    Parameters
    ----------
    bdf_filename : str
        the path to the bdf input file
    eids : list
        the element ids to consider

    """
    model = read_bdf(bdf_filename, xref=False)

    if eids is None:
        eids = model.elements.keys()
    neids = len(eids)
    nids = np.zeros((neids, 3), dtype='int32')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]
        nids[i, :] = elem.nodes

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
            print('eid=%s exists already...' % eid)
            eids_to_remove.add(eid)
        else:
            aset.add(new_row)
        #print('aset =', aset)
    return model, eids_to_remove

#def main():
    #"""the test case"""
    #eids = parse_patran_syntax_dict(' Element 830:84798')['Element']
    #eids = find_coplanar_elements(bdf_filename, eids)
