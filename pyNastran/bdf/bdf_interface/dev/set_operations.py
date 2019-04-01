"""
defines:
 - array2 = intersect1d_multi(array_set, assume_unique=False)
"""
import numpy as np

def intersect1d_multi(array_set, assume_unique=False):
    """performs multiple boolean operations"""
    if len(array_set) == 2:
        ar1 = np.intersect1d(array_set[0], array_set[1], assume_unique=assume_unique)
    else:
        lens = [a.size for a in array_set]
        isort = np.argsort(lens)
        i0 = isort[0]
        i1 = isort[1]
        if assume_unique:
            ar1 = np.unique(array_set[i0])
            ar2 = np.unique(array_set[i1])
        else:
            ar1 = array_set[i0]
            ar2 = array_set[i1]
        ar1 = np.intersect1d(ar1, ar2, assume_unique=True)
        for i in isort[1:-1]:
            if assume_unique:
                ar2 = array_set[i + 1]
            else:
                ar2 = np.unique(array_set[i + 1])
            ar1 = np.intersect1d(ar1, ar2, assume_unique=True)
    return ar1
