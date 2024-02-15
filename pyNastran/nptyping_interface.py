
def get_numpy_typing():
    """
    https://stackoverflow.com/questions/71109838/numpy-typing-with-specific-shape-and-datatype
    https://stackoverflow.com/questions/66349242/specific-type-annotation-for-numpy-ndarray-using-mypy
    """
    from typing import Optional, Annotated, Literal #, TypeVar
    import numpy as np
    import numpy.typing as npt
    #npt.NDArray[np.complex64]
    #DType = TypeVar("DType", bound=np.generic)
    #Array4 = Annotated[npt.NDArray[DType], Literal[4]]
    #Array3x3 = Annotated[npt.NDArray[DType], Literal[3, 3]]
    #ArrayNxNx3 = Annotated[npt.NDArray[DType], Literal["N", "N", 3]]

    FloatArray = npt.NDArray[np.float64]
    IntArray = npt.NDArray[np.int_]
    BoolArray = npt.NDArray[np.bool_]

    # vector:
    # [0., 0., 0.]
    NDArray3float = Annotated[FloatArray, Literal[3]]
    NDArray33float = Annotated[FloatArray, Literal[3, 3]]
    test = Optional[NDArray3float]
    del test

    # vector:
    # [0., 0., 0., 0., 0., 0.]
    NDArray6float = Annotated[FloatArray, Literal[6]]
    NDArray66float = Annotated[FloatArray, Literal[6, 6]]

    # vector:
    # [0., 0., 0., ...]
    NDArrayNbool = Annotated[BoolArray, Literal['N']]
    NDArrayNfloat = Annotated[FloatArray, Literal['N']]
    NDArrayNint = Annotated[IntArray, Literal['N']]

    # float matrix:
    # [1., 0., 0.]
    # [2., 0., 0.]
    # [3., 0., 0.]
    NDArrayN3float = Annotated[FloatArray, Literal['N', 3]]
    NDArrayN4float = Annotated[FloatArray, Literal['N', 4]]
    NDArrayNNfloat = Annotated[FloatArray, Literal['N', 'N']]

    # int matrix:
    # [1, 0]
    # [2, 0]
    # [3, 0]
    NDArrayN2int = Annotated[IntArray, Literal['N', 2]]

    # int matrix:
    # [1, 0, 0]
    # [2, 0, 0]
    # [3, 0, 0]
    NDArrayN3int = Annotated[IntArray, Literal['N', 3]]

    # int matrix:
    # [1, 0, 0, 0]
    # [2, 0, 0, 0]
    # [3, 0, 0, 0]
    NDArrayN4int = Annotated[IntArray, Literal['N', 4]]
    NDArrayN6int = Annotated[IntArray, Literal['N', 6]]
    NDArrayN8int = Annotated[IntArray, Literal['N', 8]]
    NDArrayN9int = Annotated[IntArray, Literal['N', 9]]
    NDArrayN10int = Annotated[IntArray, Literal['N', 10]]
    NDArrayN13int = Annotated[IntArray, Literal['N', 13]]
    NDArrayN15int = Annotated[IntArray, Literal['N', 15]]
    NDArrayN20int = Annotated[IntArray, Literal['N', 20]]

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




import numpy as np
if np.__version__ > '1.20' and np.__version__ < '2.0':
    out = get_numpy_typing()
    #import nptyping
    #assert nptyping.__version__ > '2.0'
    #out = get_nptyping_v2()
else:
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
