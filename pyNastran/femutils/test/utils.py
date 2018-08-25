import numpy as np

def is_array_close(v1, v2):
    """are two arrays close"""
    return np.all(np.isclose(v1, v2))
