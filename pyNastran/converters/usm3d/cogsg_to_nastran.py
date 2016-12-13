from __future__ import print_function
from itertools import count
from pyNastran.converters.usm3d.usm3d_reader import Usm3d
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_8 import print_card_8

def bc_to_nastran(bc_filename, nastran_filename):
    model = Usm3d(log=None, debug=None)
    nodes, tets = model.read_cogsg(cogsg_file, stop_after_header=False)
    assert tets.min() == 0, tets.min()
    header, tris, bcs = model.read_bc(bc_filename)
    nids = np.unique(tris.ravel())

    with open(nastran_filename, 'w') as bdf_file:
        for nid in nids:
            x, y, z = nodes[nid, :]
            bdf_file.write(print_card_16(['GRID', nid + 1, '', x, y, z]))

        for itri, tri, bc in zip(count(), tris, bcs):
            bdf_file.write(
                print_card_8(['CTRIA3', itri + 1, bc] + list(tri))
            )

        for ibc in np.unique(bcs):
            bdf_file.write(print_card_16(['PSHELL', bc, 0.1]))
