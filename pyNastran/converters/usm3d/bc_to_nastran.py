from __future__ import print_function
from itertools import count
import numpy as np

from pyNastran.converters.usm3d.usm3d_reader import Usm3d
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_8 import print_card_8

def bc_to_nastran(cogsg_filename, bc_filename, nastran_filename):
    """
    converts a *.bc file to a *.bdf file
    """
    model = Usm3d(log=None, debug=None)
    nodes, tets = model.read_cogsg(cogsg_filename, stop_after_header=False)
    assert tets.min() == 0, tets.min()
    header, tris, bcs = model.read_bc(bc_filename)
    nids = np.unique(tris.ravel())

    with open(nastran_filename, 'w') as bdf_file:
        bdf_file.write('$ pyNastran : punch=True\n')
        for nid in nids:
            x, y, z = nodes[nid, :]
            bdf_file.write(print_card_16(['GRID', nid + 1, '', x, y, z]))

        for itri, tri, bc in zip(count(), tris + 1, bcs):
            bdf_file.write(
                print_card_8(['CTRIA3', itri + 1, bc] + list(tri))
            )

        mid = 1
        for bc in np.unique(bcs):
            bdf_file.write(print_card_8(['PSHELL', bc, mid, 0.1]))
        bdf_file.write(print_card_8(['MAT1', mid, 3.0e7, None, 0.3]))
