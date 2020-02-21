try:
    from nptyping import Array
    # vector:
    # [0., 0., 0.]
    NDArray3float = Array[float, 3]

    # vector:
    # [0., 0., 0., 0., 0., 0.]
    NDArray6float = Array[float, 6]

    # vector:
    # [0., 0., 0., ...]
    NDArrayNbool = Array[bool, ...]
    NDArrayNfloat = Array[float, ...]
    NDArrayNint = Array[int, ...]

    # float matrix:
    # [1., 0., 0.]
    # [2., 0., 0.]
    # [3., 0., 0.]
    NDArrayN3float = Array[float, ..., 3]

    NDArrayNNfloat = Array[float, ..., ...]

    # int matrix:
    # [1, 0]
    # [2, 0]
    # [3, 0]
    NDArrayN2int = Array[int, ..., 2]

    # int matrix:
    # [1, 0, 0]
    # [2, 0, 0]
    # [3, 0, 0]
    NDArrayN3int = Array[int, ..., 3]

    # int matrix:
    # [1, 0, 0, 0]
    # [2, 0, 0, 0]
    # [3, 0, 0, 0]
    NDArrayN4int = Array[int, ..., 4]
    NDArrayN6int = Array[int, ..., 6]
    NDArrayN8int = Array[int, ..., 8]
    NDArrayN9int = Array[int, ..., 9]
    NDArrayN10int = Array[int, ..., 10]
    NDArrayN13int = Array[int, ..., 13]
    NDArrayN15int = Array[int, ..., 15]
    NDArrayN20int = Array[int, ..., 20]

except ImportError:
    from typing import Any
    NDArrayNbool = \
        NDArray3float = NDArray6float = \
        NDArrayNint = NDArrayNfloat = \
        NDArrayN3float = \
        NDArrayN2int = NDArrayN3int = NDArrayN4int = NDArrayN6int = NDArrayN8int = NDArrayN9int = \
        NDArrayN10int = NDArrayN13int = NDArrayN15int = NDArrayN20int = \
        NDArrayNNfloat = Any
