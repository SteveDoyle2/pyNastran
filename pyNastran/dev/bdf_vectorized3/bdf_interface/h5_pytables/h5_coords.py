from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from .utils import get_group_name
if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.cards.coord import COORD
    from tables import Group


def load_h5_coord(model: BDF, group: Group):
    coord_ids = []
    reference_ids = []
    a1 = []
    a2 = []
    a3 = []

    b1 = []
    b2 = []
    b3 = []

    c1 = []
    c2 = []
    c3 = []
    nodes = []
    coord_str = []
    coord_int = []
    domain_id = []
    for h5_element in group._f_iter_nodes():
        class_id = h5_element._c_classid
        if class_id == 'GROUP':
            class_name = get_group_name(h5_element)
            assert class_name in {'TRANSFORMATION'}, class_name
            continue
        elif class_id == 'TABLE':
            pass
        else:
            raise RuntimeError(class_id)
        name = h5_element.name
        #print(f'loading {name}')
        coord_inti = int(name[-2])
        coord_stri = name[-1]
        data = h5_element.read()
        cid = data['CID']
        ncoordsi = len(cid)
        coord_int.append(np.full(ncoordsi, coord_inti, dtype='int32'))
        coord_str.append(np.full(ncoordsi, coord_stri, dtype='|U1'))
        if coord_inti == 1:
            cord1
        elif coord_inti == 2:
            #dtype([('CID', '<i8'), ('RID', '<i8'), ('A1', '<f8'), ('A2', '<f8'), ('A3', '<f8'),
            #       ('B1', '<f8'), ('B2', '<f8'), ('B3', '<f8'), ('C1', '<f8'), ('C2', '<f8'), ('C3', '<f8'),
            #       ('DOMAIN_ID', '<i8')])
            coord_ids.append(cid)
            reference_ids.append(data['RID'])
            nodesi = np.full((ncoordsi, 3), -1, dtype=cid.dtype)
            nodes.append(nodesi)
            a1.append(data['A1'])
            a2.append(data['A2'])
            a3.append(data['A3'])

            b1.append(data['B1'])
            b2.append(data['B2'])
            b3.append(data['B3'])

            c1.append(data['C1'])
            c2.append(data['C2'])
            c3.append(data['C3'])
            domain_id.append(data['DOMAIN_ID'])
        else:
            raise RuntimeError(name)

    if coord_ids:
        coord_id = np.hstack(coord_ids)
        if 0 not in coord_id:
            coord_ids.append([0])
            coord_id = np.hstack(coord_ids)

            coord_str.append(['R'])
            coord_int.append([2])
            reference_ids.append([0])
            nodes.append([-1, -1, -1])
            a1.append([0.])
            a2.append([0.])
            a3.append([0.])

            b1.append([0.])
            b2.append([0.])
            b3.append([1.])

            c1.append([1.])
            c2.append([0.])
            c3.append([0.])
        ncoord = len(coord_id)

        coord = model.coord
        coord_type = np.hstack(coord_str)
        icoord = np.hstack(coord_int)

        nodes = np.vstack(nodes)
        #  np.full((1, 3), -1, dtype='int32')

        ref_coord_id = np.hstack(reference_ids)
        e1 = np.hstack([
            np.hstack(a1).reshape(ncoord, 1),
            np.hstack(a2).reshape(ncoord, 1),
            np.hstack(a3).reshape(ncoord, 1),
        ])
        e2 = np.hstack([
            np.hstack(b1).reshape(ncoord, 1),
            np.hstack(b2).reshape(ncoord, 1),
            np.hstack(b3).reshape(ncoord, 1),
        ])
        e3 = np.hstack([
            np.hstack(c1).reshape(ncoord, 1),
            np.hstack(c2).reshape(ncoord, 1),
            np.hstack(c3).reshape(ncoord, 1),
        ])

        coord_save(coord, coord_id, coord_type, icoord, ref_coord_id, nodes, e1, e2, e3)
        coord.write()
    #x = 1

def coord_save(self: COORD, coord_id, coord_type, icoord, ref_coord_id, nodes, e1, e2, e3):
    ncoord = len(coord_id)
    self.coord_id = coord_id
    self.coord_type = coord_type
    self.is_resolved = (coord_id == 0)
    self.icoord = icoord

    self.ref_coord_id = ref_coord_id
    self.e1 = e1
    self.e2 = e2
    self.e3 = e3
    self.nodes = nodes
    self.origin   = np.full((ncoord, 3), np.nan, dtype='float64')
    self.z_axis   = np.full((ncoord, 3), np.nan, dtype='float64')
    self.xz_plane = np.full((ncoord, 3), np.nan, dtype='float64')

    self.i = np.full((ncoord, 3), np.nan, dtype='float64')
    self.j = np.full((ncoord, 3), np.nan, dtype='float64')
    self.k = np.full((ncoord, 3), np.nan, dtype='float64')
    self.T = np.full((ncoord, 3, 3), np.nan, dtype='float64')
    self.xyz_to_global_transform = {cid: Ti for cid, Ti in zip(coord_id, self.T)}
    self._set_global_coord()

    self.n = ncoord
    self.is_resolved = (coord_id == 0)
    self.sort()
    self.write(size=8)
    #self._save(coord_id, ref_coord_id, e1, e2, e3)
    #self.model.log.warning('adding nan coords')
    #x = 1

