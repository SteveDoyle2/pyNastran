import numpy as np


def chexa_centroid(self) -> np.ndarray:
    xyz = self.model.grid.xyz_cid0()
    nid = self.model.grid.node_id
    nodes = self.base_nodes
    #nelement = nodes.shape[0]

    inode = np.searchsorted(nid, nodes)
    n1 = xyz[inode[:, 0], :]
    n2 = xyz[inode[:, 1], :]
    n3 = xyz[inode[:, 2], :]
    n4 = xyz[inode[:, 3], :]
    n5 = xyz[inode[:, 4], :]
    n6 = xyz[inode[:, 5], :]
    n7 = xyz[inode[:, 6], :]
    n8 = xyz[inode[:, 7], :]
    centroid = (n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8) / 8
    return centroid
