def get_nptyping_v2():
    from typing import Optional
    from nptyping import NDArray, Float, Int, Bool, Shape
    # vector:
    # [0., 0., 0.]
    NDArray3float = NDArray[Shape['3'], Float]
    NDArray33float = NDArray[Shape['3, 3'], Float]
    test = Optional[NDArray3float]
    del test

    # vector:
    # [0., 0., 0., 0., 0., 0.]
    NDArray6float = NDArray[Shape['6'], Float]
    NDArray66float = NDArray[Shape['6, 6'], Float]

    # vector:
    # [0., 0., 0., ...]
    NDArrayNbool = NDArray[Shape['Size'], Bool]
    NDArrayNfloat = NDArray[Shape['Size'], Float]
    NDArrayNint = NDArray[Shape['Size'], Int]

    # float matrix:
    # [1., 0., 0.]
    # [2., 0., 0.]
    # [3., 0., 0.]
    NDArrayN3float = NDArray[Shape['Size, 3'], Float]
    NDArrayN4float = NDArray[Shape['Size, 4'], Float]

    NDArrayNNfloat = NDArray[Shape['Size, Size'], Float]

    # int matrix:
    # [1, 0]
    # [2, 0]
    # [3, 0]
    NDArrayN2int = NDArray[Shape['Size, 2'], Int]

    # int matrix:
    # [1, 0, 0]
    # [2, 0, 0]
    # [3, 0, 0]
    NDArrayN3int = NDArray[Shape['Size, 3'], Int]

    # int matrix:
    # [1, 0, 0, 0]
    # [2, 0, 0, 0]
    # [3, 0, 0, 0]
    NDArrayN4int = NDArray[Shape['Size, 4'], Int]
    NDArrayN6int = NDArray[Shape['Size, 6'], Int]
    NDArrayN8int = NDArray[Shape['Size, 8'], Int]
    NDArrayN9int = NDArray[Shape['Size, 9'], Int]
    NDArrayN10int = NDArray[Shape['Size, 10'], Int]
    NDArrayN13int = NDArray[Shape['Size, 13'], Int]
    NDArrayN15int = NDArray[Shape['Size, 15'], Int]
    NDArrayN20int = NDArray[Shape['Size, 20'], Int]

    out = (
        NDArrayNbool,
        NDArray3float, NDArray6float,
        NDArrayNint, NDArrayNfloat,
        NDArray66float,
        NDArrayN3float, NDArrayN4float,
        NDArray33float,
        NDArrayN2int, NDArrayN3int, NDArrayN4int, NDArrayN6int, NDArrayN8int, NDArrayN9int,
        NDArrayN10int, NDArrayN13int, NDArrayN15int, NDArrayN20int,
        NDArrayNNfloat,
    )
    return out

def get_nptyping_v1():
    from typing import Optional, Any
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
    NDArray66float = NDArray[(6, 6), float]

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
    out = (
        NDArrayNbool,
        NDArray3float, NDArray6float,
        NDArrayNint, NDArrayNfloat,
        NDArray66float,
        NDArrayN3float, NDArrayN4float,
        NDArray33float,
        NDArrayN2int, NDArrayN3int, NDArrayN4int, NDArrayN6int, NDArrayN8int, NDArrayN9int,
        NDArrayN10int, NDArrayN13int, NDArrayN15int, NDArrayN20int,
        NDArrayNNfloat,
    )
    return out

def get_nptyping_null():
    from typing import Any
    NDArrayNbool = \
    NDArray3float = NDArray6float = \
    NDArrayNint = NDArrayNfloat = \
    NDArray66float = \
    NDArrayN3float = NDArrayN4float = \
    NDArray33float = \
    NDArrayN2int = NDArrayN3int = NDArrayN4int = NDArrayN6int = NDArrayN8int = NDArrayN9int = \
    NDArrayN10int = NDArrayN13int = NDArrayN15int = NDArrayN20int = \
    NDArrayNNfloat = Any

    out = (
        NDArrayNbool,
        NDArray3float, NDArray6float,
        NDArrayNint, NDArrayNfloat,
        NDArray66float,
        NDArrayN3float, NDArrayN4float,
        NDArray33float,
        NDArrayN2int, NDArrayN3int, NDArrayN4int, NDArrayN6int, NDArrayN8int, NDArrayN9int,
        NDArrayN10int, NDArrayN13int, NDArrayN15int, NDArrayN20int,
        NDArrayNNfloat
    )
    return out

try:
    import nptyping
    if nptyping.__version__ > '2.0':
        out = get_nptyping_v2()
    else:
        out = get_nptyping_v1()

except ImportError:
    out = get_nptyping_null()

(
    NDArrayNbool,
    NDArray3float, NDArray6float,
    NDArrayNint, NDArrayNfloat,
    NDArray66float,
    NDArrayN3float, NDArrayN4float,
    NDArray33float,
    NDArrayN2int, NDArrayN3int, NDArrayN4int, NDArrayN6int, NDArrayN8int, NDArrayN9int,
    NDArrayN10int, NDArrayN13int, NDArrayN15int, NDArrayN20int,
    NDArrayNNfloat,
) = out
