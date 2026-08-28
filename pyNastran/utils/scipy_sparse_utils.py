"""
hack file to fake scipy matrices in order to parse op2s without scipy

Limitations:
 - no hdf5 export
 - no matrix multiplication
 - no op4 writing
 - no trapezoid, quad functions
 - no nodal equivalencing

"""
import warnings

from pyNastran.utils import int_version
from pyNastran.utils.numpy_utils import integer_types
try:
    import scipy
    IS_SCIPY = True
except ImportError as e:
    warnings.warn('Couldnt find scipy...faking')
    IS_SCIPY = False

if IS_SCIPY:
    SCIPY_VERSION = int_version('scipy', scipy.__version__)[:2]

    import scipy.sparse
    # address scipy.sparse.coo -> scipy.sparse._coo
    IS_NEW_SCIPY = (SCIPY_VERSION >= [1, 8])
    IS_OLD_SCIPY = not IS_NEW_SCIPY

    coo_matrix = scipy.sparse.coo_matrix
    csr_matrix = scipy.sparse.csr_matrix
    csc_matrix = scipy.sparse.csc_matrix
    sparse_types = (coo_matrix, csr_matrix, csc_matrix)
else:
    import numpy as np
    class coo_matrix:
        def __init__(self, data_my_indices, shape=None, dtype=None):
            """
            coo_matrix(
                (nrows, ncols), dtype=dtype)
            coo_matrix(
                (real_array, (GCi, GCj)),
                shape=(mrows, ncols), dtype=dtype)
            """
            # print(f'ndata = {len(data_my_indices)}')
            # print(f'data_my_indices = {data_my_indices}')
            if len(data_my_indices) == 2:
                arg1, arg2 = data_my_indices
                if isinstance(arg1, (np.ndarray, list)):
                    if isinstance(arg1, (np.ndarray, list)):
                        arg1 = np.asarray(arg1, dtype=dtype)
                        data = arg1
                        assert isinstance(data, np.ndarray), data

                    if len(arg2) == 2:
                        # coo_matrix(
                        #     (real_array, (GCi, GCj)),
                        #     shape=(mrows, ncols), dtype=dtype)
                        row, col = arg2
                    else:
                        raise NotImplementedError(arg2)
                        #row_col = arg2
                        #nrow, ncol = row_col
                        #assert isinstance(nrow, integer_types), (nrow, type(nrow))
                        #assert isinstance(ncol, integer_types), arg2
                        #shape = (nrow, ncol)
                        #nvalue = nrow * ncol
                        #assert isinstance(nvalue, integer_types), nvalue
                        #values = np.arange(nvalue)
                        #row = values.reshape(shape)
                        #col = values.reshape((ncol, nrow)).T
                    # else:
                    #     raise NotImplementedError(arg2)
                elif isinstance(arg1, integer_types):
                    # data_mat = coo_matrix((nrows, ncols), dtype=dtype)
                    assert isinstance(arg2, integer_types), arg2
                    shape = (arg1, arg2)
                    data = []
                    row = []
                    col = []
                else:  # pragma: no cover
                    raise NotImplementedError((arg1, type(arg1)))
            else:  # pragma: no cover
                raise NotImplementedError(data_my_indices)

            if dtype is None:
                if isinstance(data, np.ndarray):
                    dtype = data.dtype
                else:  # pragma: no cover
                    raise NotImplementedError((data, type(data), dtype))
            #assert dtype is not None, dtype

            self.data = np.asarray(data, dtype=dtype)
            self.row = row
            self.col = col
            assert len(self.row) == len(self.data)
            assert len(self.row) == len(self.col)
            assert shape is not None, shape
            assert len(shape) == 2, shape
            self.shape = tuple(np.asarray(shape, dtype='int32').tolist())
            self.dtype = self.data.dtype
        # @property
        # def shape(self) -> tuple[int, int]:
        #     return self.data.shape
        @property
        def indices(self) -> np.ndarray:
            return np.column_stack([self.row, self.col])
        @property
        def real(self):
            return self.data.real
        @property
        def imag(self):
            return self.data.imag
        def toarray(self):
            mat = np.zeros(self.shape, dtype=self.dtype)
            row_col = (self.row, self.col)
            mat[row_col] = self.data
            return mat
        @property
        def nnz(self) -> int:
            i = np.where(self.data != 0.)[0]
            return len(i)
        def tocsr(self):
            return self
        def tocsc(self):
            return self
        @property
        def has_sorted_indices(self):
            return False
        def sort_indices(self):
            pass
        def __repr__(self):
            # <COOrdinate sparse matrix of dtype 'float32'
            #         with 240 stored elements and shape (30, 20)>
            # shape = tuple(self.shape.tolist())
            msg = (
                f"'<fake.COOrdinate sparse matrix of dtype {self.dtype.name!r}\n"
                f"        with {self.nnz} stored elements and shape {self.shape}>"
            )
            return msg

    csr_matrix = csc_matrix = coo_matrix
    sparse_types = (coo_matrix, csr_matrix, csc_matrix)
