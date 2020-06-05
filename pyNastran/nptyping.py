from typing import Any, Optional
try:
    from nptyping import NDArray
    # vector:
    # [0., 0., 0.]
    NDArray3float = NDArray[(3, ), float]
    NDArray33float = NDArray[(3, 3), float]
    test = Optional[NDArray3float]
    del test

    # vector:
    # [0., 0., 0., 0., 0., 0.]
    NDArray6float = NDArray[(6, ), float]

    # vector:
    # [0., 0., 0., ...]
    NDArrayNbool = NDArray[(Any, ), bool]
    NDArrayNfloat = NDArray[(Any, ), float]
    NDArrayNint = NDArray[(Any, ), int]

    # float matrix:
    # [1., 0., 0.]
    # [2., 0., 0.]
    # [3., 0., 0.]
    NDArrayN3float = NDArray[(Any, 3), float]
    NDArrayN4float = NDArray[(Any, 4), float]

    NDArrayNNfloat = NDArray[(Any, Any), float]

    # int matrix:
    # [1, 0]
    # [2, 0]
    # [3, 0]
    NDArrayN2int = NDArray[(Any, 2), int]

    # int matrix:
    # [1, 0, 0]
    # [2, 0, 0]
    # [3, 0, 0]
    NDArrayN3int = NDArray[(Any, 3), int]

    # int matrix:
    # [1, 0, 0, 0]
    # [2, 0, 0, 0]
    # [3, 0, 0, 0]
    NDArrayN4int = NDArray[(Any, 4), int]
    NDArrayN6int = NDArray[(Any, 6), int]
    NDArrayN8int = NDArray[(Any, 8), int]
    NDArrayN9int = NDArray[(Any, 9), int]
    NDArrayN10int = NDArray[(Any, 10), int]
    NDArrayN13int = NDArray[(Any, 13), int]
    NDArrayN15int = NDArray[(Any, 15), int]
    NDArrayN20int = NDArray[(Any, 20), int]

except ImportError:
    NDArrayNbool = \
    NDArray3float = NDArray6float = \
    NDArrayNint = NDArrayNfloat = \
    NDArrayN3float = NDArrayN4float = \
    NDArray33float = \
    NDArrayN2int = NDArrayN3int = NDArrayN4int = NDArrayN6int = NDArrayN8int = NDArrayN9int = \
    NDArrayN10int = NDArrayN13int = NDArrayN15int = NDArrayN20int = \
    NDArrayNNfloat = Any
