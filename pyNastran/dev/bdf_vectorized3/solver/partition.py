import numpy as np
from pyNastran.nptyping_interface import (
    NDArrayNfloat, NDArrayNNfloat)


def partition_matrix(matrix: np.ndarray,
                     sets: list[tuple[str, np.ndarray]],
                     ) -> dict[tuple[str, str], NDArrayNNfloat]:
    """partitions a matrix"""
    matrices = {}
    for aname, aset in sets:
        for bname, bset in sets:
            matrices[aname + bname] = matrix[aset, :][:, bset]
    return matrices


def partition_vector(vector: np.ndarray,
                     sets: list[tuple[str, np.ndarray]],
                     fdtype: str='float64') -> list[NDArrayNfloat]:  # pragma: no cover
    """partitions a vector"""
    vectors = []
    for unused_aname, aset in sets:
        if len(aset) == 0:
            vectors.append(np.array([], dtype=fdtype))
            continue
        vectori = vector[aset]
        vectors.append(vectori)
    return vectors


def partition_vector2(vector: np.ndarray,
                      sets: list[tuple[str, np.ndarray]],
                      fdtype: str='float64') -> tuple[NDArrayNfloat, NDArrayNfloat]:
    """partitions a vector"""
    assert len(sets) == 2, sets
    #vectors = partition_vector(vector, sets, fdtype=fdtype)
    #return tuple(vectors)
    (unused_name0, set0), (unused_name1, set1) = sets
    vectors = (vector[set0], vector[set1])
    return vectors


def partition_vector3(vector: np.ndarray,
                      sets: list[tuple[str, np.ndarray]],
                      fdtype: str='float64') -> tuple[NDArrayNfloat, NDArrayNfloat, NDArrayNfloat]:
    """partitions a vector"""
    assert len(sets) == 3, sets
    #vectors = partition_vector(vector, sets, fdtype=fdtype)
    (unused_name0, set0), (unused_name1, set1), (unused_name2, set2) = sets
    vectors = (vector[set0], vector[set1], vector[set2])
    return tuple(vectors)
