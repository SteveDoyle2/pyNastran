from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from tables import Group


def load_h5_node(model: BDF, input_group: Group):
    for h5_element in input_group._f_iter_nodes():
        name = h5_element.name
        #print(f'loading {name}')
        if name == 'GRID':
            grid = model.grid
            data = h5_element.read()
            #dtype=[('ID', '<i8'), ('CP', '<i8'), ('X', '<f8', (3,)), ('CD', '<i8'),
            # ('PS', '<i8'), ('SEID', '<i8'), ('DOMAIN_ID', '<i8')])
            grid.node_id = data['ID']
            grid.cp = data['CP']
            grid.cd = data['CD']
            grid.xyz = data['X']
            grid.ps = data['PS']
            grid.seid = data['SEID']
            grid.domain_id = data['DOMAIN_ID']
            grid.n = len(grid.node_id)
            grid._xyz_cid0 = np.full((grid.n, 3), np.nan, dtype=grid.xyz.dtype)
            icp0 = np.where(grid.cp == 0)[0]
            grid._xyz_cid0[icp0, :] = grid.xyz[icp0, :]
            grid.write()
        elif name == 'EPOINT':
            epoint = model.epoint
            data = h5_element.read()
            #dtype([('ID', '<i8'), ('DOMAIN_ID', '<i8')])
            epoint.ids = data['ID']
            epoint.domain_id = data['DOMAIN_ID']
            epoint.n = len(epoint.domain_id)
            epoint.write()
        elif name == 'SPOINT':
            spoint = model.spoint
            data = h5_element.read()
            #dtype([('ID', '<i8'), ('DOMAIN_ID', '<i8')])
            spoint.ids = data['ID']
            spoint.domain_id = data['DOMAIN_ID']
            spoint.n = len(spoint.domain_id)
            spoint.write()
        else:
            raise NotImplementedError(name)
        x = 1
    x = 2


