import numpy as np


class Case2D:
    def __init__(self, node_id: np.ndarray, data: np.ndarray):
        nnode = len(node_id)
        self.node_gridtype = np.zeros((nnode, 2), dtype='int32')
        self.node_gridtype[:, 0] = node_id
        self.data = data.reshape(1, nnode, 3)
