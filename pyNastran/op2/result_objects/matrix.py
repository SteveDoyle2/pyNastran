"""Defines the Matrix class"""
from __future__ import annotations
from typing import Callable, Optional, TextIO, Any, TYPE_CHECKING
from itertools import count
import scipy.sparse
#from scipy.sparse import coo_matrix, csr_matrix  # type: ignore
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.dmig import DMI, dtype_to_tin_tout_str
from pyNastran.bdf.field_writer import print_card_8, print_card_16, print_card_double
from pyNastran.op2.op2_interface.write_utils import export_to_hdf5
from pyNastran.utils import object_attributes, object_methods, object_stats
sparse_types = (scipy.sparse.coo_matrix, scipy.sparse.csr_matrix, scipy.sparse.csc_matrix)
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger


class Matrix:
    """
    Defines a Matrix object that's stored in:
     - op4.matrices
     - op2.matrices

    Attributes
    ----------
    name : str
        the name of the matrix
    data : varies
        dense : np.ndarray
        sparse : coo_matrix
        data is generally initialized by setting the matrix.data attribute externally
    is_matpool : bool
        is this a matpool matrix?  A matpool has (grid, component) values
        similar to a DMIG.  A non-matpool matrix is similar to a DMI.

    """
    def __init__(self, name: str, form: int | str,
                 data: Optional[np.ndarray | scipy.sparse.coo_matrix]=None):
        """
        Initializes a Matrix

        Parameters
        ----------
        name : str
            the name of the matrix
        form : int
            the matrix type
        is_matpool : bool
            is this a matpool matrix?  A matpool has (grid, component) values
            similar to a DMIG.  A non-matpool matrix is similar to a DMI.

        +------+-----------------+
        | Form | Meaning         |
        +======+=================+
        |  1   | Square          |
        |  2   | Rectangular     |
        |  6   | Symmetric       |
        |  9   | Pseudo identity |
        +------+-----------------+

        """
        self.name = name
        self.data = data
        if isinstance(form, integer_types):
            pass
        else:
            form = form_to_int(form)
        self.form: int = form
        self.is_matpool = False

        # only exist for is_matpool = True; automatically set
        self.col_nid = None
        self.col_dof = None
        self.row_nid = None
        self.row_dof = None
        if not isinstance(name, str):  # pragma: no cover
            raise TypeError(f'name={name!r} must be a string; type={type(name)}')

    def set_matpool_data(self,
                         data: np.ndarray,
                         col_nid: np.ndarray, col_dof: np.ndarray,
                         row_nid: np.ndarray, row_dof: np.ndarray) -> None:
        self.is_matpool = True
        self.data = data
        self.col_nid = col_nid
        self.col_dof = col_nid
        self.row_nid = row_nid
        self.row_dof = row_dof

    def symmetric_to_rectangular(self, log: SimpleLogger):
        """enforce symmetry if necessary"""
        if self.shape_str != 'symmetric':
            return
        if not isinstance(self.data, sparse_types):
            return

        matrix = sparse_symmetric_to_rectangular(self.name, self.data, log)
        self.data = matrix

        # matrix is symmetric, but is not stored as symmetric
        self.matrix_shape = 'rectangular'

    @property
    def shape_str(self) -> str:
        """gets the matrix description"""
        if self.form == 0:
            return 'N/A'
        if self.form == 1:
            return 'square'
        elif self.form == 2:
            return 'rectangular'
        elif self.form == 6:
            return 'symmetric'
        elif self.form == 9:
            return 'pseudo-identity'
        else:
            raise RuntimeError(f'form = {self.form!r}')

    def export_to_hdf5(self, group, log: SimpleLogger) -> None:
        """exports the object to HDF5 format"""
        export_to_hdf5(self, group, log)

    def build_dataframe(self) -> None:
        """exports the object to pandas format"""
        import pandas as pd
        matrix = self.data
        if matrix is None:
            return
        if isinstance(matrix, scipy.sparse.coo_matrix):
            data = {'row': matrix.row, 'col': matrix.col, 'data' : matrix.data}
            data_frame = pd.DataFrame(data=data).reindex(columns=['row', 'col', 'data'])
        elif isinstance(matrix, np.ndarray):
            data_frame = pd.DataFrame(data=matrix)
        else:
            raise NotImplementedError(type(matrix))
        self.data_frame = data_frame

    def write(self, mat: TextIO, print_full: bool=True) -> None:
        """writes to the F06"""
        mat.write(str(self) + '\n')

        matrix = self.data
        if self.data is None:
            skip_msg = f'skipping {self.name!r} because data is None\n\n'
            mat.write(skip_msg.encode('ascii'))
            return
        if isinstance(matrix, scipy.sparse.coo_matrix):
            if print_full:
                for row, col, value in zip(matrix.row, matrix.col, matrix.data):
                    mat.write(f'({row:d}, {col:d}) {value}\n')
            else:
                mat.write(str(matrix))
        else:
            mat.write('name=%r; shape=%s; form=%d; Type=%r\n' % (
                self.name, str(self.data.shape).replace('L', ''),
                self.form, self.shape_str))
            if print_full:
                np.savetxt(mat, self.data, fmt='%.18e', delimiter=',')
            #f06.write(str(matrix))
            #print('WARNING: matrix type=%s does not support writing' % type(matrix))
        mat.write('\n\n')

    def get_stats(self, mode: str='public', keys_to_skip=None,
                  filter_properties: bool=False):
        if keys_to_skip is None:
            keys_to_skip = []

        my_keys_to_skip = ['object_methods', 'object_attributes',]
        return object_stats(self, mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                            filter_properties=filter_properties)

    def object_attributes(self, mode: str='public', keys_to_skip=None,
                          filter_properties: bool=False):
        if keys_to_skip is None:
            keys_to_skip = []

        my_keys_to_skip = ['object_methods', 'object_attributes',]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(self, mode: str='public', keys_to_skip=None):
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = []

        my_keys_to_skip = ['object_methods', 'object_attributes',]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def __repr__(self) -> str:
        header = f'Matrix[{self.name!r}];'
        if self.data is None:
            shape = 'data=None; '
            class_name = '<NoneType>'
            dtype = '<NoneType>; '
        else:
            class_name = str(type(self.data)).replace('<class ', '').replace('>', '').replace("'", '') + ';'
            shape = ' shape=%s;' % str(self.data.shape).replace('L', '')
            dtype = '%s;' % self.data.dtype
        msg = '%-18s %-18s type=%-33s dtype=%-10s desc=%s' % (
            header, shape, class_name, dtype, self.shape_str)
        return msg

    @property
    def dtype_str(self) -> str:
        assert isinstance(self.data, (np.ndarray, scipy.sparse.coo_matrix)), type(self.data)
        return self.data.dtype.name

    @property
    def GCi(self) -> np.ndarray:
        return  self.GCi_GCj[0]

    @property
    def GCj(self) -> np.ndarray:
        return self.GCi_GCj[1]

    @property
    def GCi_GCj(self) -> tuple[np.ndarray, np.ndarray]:
        if isinstance(self.data, np.ndarray):
            # same for real/imaginary
            i, j = np.where(self.data != 0.)
        elif isinstance(self.data, scipy.sparse.coo_matrix):
            # same for real/imaginary
            i = self.data.row
            j = self.data.col
        else:
            raise TypeError(self.data)
        return i, j

    def data_i_j(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if isinstance(self.data, np.ndarray):
            i, j = np.where(self.data != 0.)
            data = self.data[i, j]
            #real = data.real[i, j]
            #if is_complex:
                #imag = data.imag[i, j]
        elif isinstance(self.data, scipy.sparse.coo_matrix):
            i = self.data.row
            j = self.data.col
            data = self.data.data
            #real = data.real
            #if is_complex:
                #imag = data.imag
        elif isinstance(self.data, sparse_types):
            coo = self.data.tocoo(copy=False)
            i = coo.row
            j = coo.col
            data = coo.data
        else:  # pragma: no cover
            raise TypeError(self.data)
        assert len(data) == len(i)
        assert len(data) == len(j)
        return data, i, j

    @property
    def Real(self) -> np.ndarray:
        return self.data_i_j()[0].real

    @property
    def Complex(self) -> np.ndarray:
        return self.data_i_j()[0].imag

    def to_gcj_gci_form(self):
        """once in matrix form, we need to transform to GCj, GCi, Real/Complex form"""
        assert isinstance(self.col_nid, np.ndarray), self.col_nid
        assert isinstance(self.col_dof, np.ndarray), self.row_dof
        assert isinstance(self.row_nid, np.ndarray), self.col_nid
        assert isinstance(self.row_dof, np.ndarray), self.row_dof

        data = self.data
        nrows, ncols = data.shape

        # fix columns to consider
        #is_real = self.dtype_str in {'float32', 'float64'}

        #sparse_csr = scipy.sparse.csr_matrix(self.data)
        #sparse_coo = sparse_csr.tocoo()  #  this should be the answer....


        max_cols = (np.abs(data.real) + np.abs(data.imag)).max(axis=1)
        assert ncols == len(max_cols)
        jcols = np.where(max_cols > 0.)[0]
        GCj = []

        ncol = len(jcols)
        col_nidj = self.col_nid[jcols]
        col_dofj = self.col_dof[jcols]
        GCj = np.vstack([col_nidj, col_dofj]).ravel(ncol, 2).T
        print(f'GCj = {GCj}')

        GCi_out = []
        GCj_out = []
        for j, gcj in zip(jcols, GCj):
            print(f'gcj = {gcj}')
            dataj = self.data[:, j]
            abs_dataj = np.abs(dataj.real) + np.abs(dataj.imag)
            irows = np.where(abs_dataj > 0)[0]
            nrow = len(irows)

            row_nidi = self.row_nid[irows]
            row_dofi = self.row_dof[irows]
            GCi = np.vstack([row_nidi, row_dofi]).ravel(nrow, 2).T
            print(f'GCi = {GCi}')
            GCj_out.append([gcj] * nrow)
            GCi_out.append(GCi)
        GCi_out_array = np.array(GCi_out, dtype='int32')
        GCj_out_array = np.array(GCi_out, dtype='int32')
        return GCi_out_array, GCj_out_array

    @property
    def tin(self) -> int:
        tin_str = dtype_to_tin_tout_str(self.data)
        tin_str_to_tin = {
            'float32': 1,
            'float64': 2,
            'complex64': 3,
            'complex128': 4,
        }
        tin = tin_str_to_tin[tin_str]
        return tin

    @property
    def is_real(self) -> bool:
        """
        1-Real, Single Precision
        2=Real, Double Precision
        3=Complex, Single
        4=Complex, Double
        """
        return self.tin in {1, 2}

    @property
    def is_complex(self) -> bool:
        return not self.is_real

    @property
    def is_sparse(self) -> bool:
        return scipy.sparse.issparse(self.data)
    @property
    def is_dense(self) -> bool:
        return not self.is_sparse

    def to_sparse(self, sparse_type: str='coo') -> None:
        if isinstance(self.data, np.ndarray):
            dtype = self.data.dtype
            data_array, i, j = self.data_i_j()
            if sparse_type == 'coo':
                data = scipy.sparse.coo_matrix(
                    (data_array, (i, j)), shape=self.shape, dtype=dtype, copy=False)
            elif sparse_type == 'csr':
                data = scipy.sparse.csr_matrix(
                    (data_array, (i, j)), shape=self.shape, dtype=dtype, copy=False)
            elif sparse_type == 'csc':
                data = scipy.sparse.csc_matrix(
                    (data_array, (i, j)), shape=self.shape, dtype=dtype, copy=False)
            else:  # pragma: no cover
                raise ValueError(f'sparse_type={sparse_type!r}; supports=[coo, csr, csc]')
            self.data = data

        elif isinstance(self.data, sparse_types):
            self.data = self.data.todense()
        else:
            raise TypeError(self.data)

    def to_dense(self) -> None:
        if isinstance(self.data, np.ndarray):
            pass
        elif isinstance(self.data, sparse_types):
            self.data = self.data.todense()
        else:
            raise TypeError(self.data)

    @property
    def shape(self) -> tuple[int, int]:
        return self.data.shape

    def write_dmi(self, size: int=8) -> str:
        """
        DMI Notes:
         - Additional Forms:
           - 3 = Diagonal matrix (elements on the diagonal are stored in a column vector having m rows)
           - 4 = Lower triangular factor
           - 5 = Upper triangular factor
           - 8 = Identity matrix (m = number of rows, n = m)
         - The total number of DMIs and DTIs may not exceed 1000.
         - For symmetric matrices, the entire matrix must be input.
         - Form 7 matrices may not be defined on this entry.
         - Form 3 matrices are converted to Form 6 matrices, which may be used by any module.
         - The DMIG entry is more convenient for matrices with rows and columns that are
           referenced by grid or scalar point degrees-of-freedom.
        """
        #if self.shape_str == 'symmetric':
            #raise ValueError('call symmetric_to_rectangular before writing a DMI')
        msg = '\n$' + '-' * 80
        msg += '\n$ %s Matrix %s\n' % ('DMI', self.name)
        nrows, ncols = self.shape

        tin = self.tin
        tout = tin

        list_fields = ['DMI', self.name, 0, self.form, tin,
                       tout, None, nrows, ncols]
        msg += print_card_8(list_fields)

        func_small: Callable[list[Any]] = print_card_8 if size == 8 else print_card_16
        func: Callable[list[Any]] = func_small if tin in {1, 3} else print_card_double

        tout = self.tin
        nrows, ncols = self.shape
        data_array, GCi, GCj = self.data_i_j()
        if self.is_complex:
            dmi = DMI(self.name, matrix_form=self.form,
                      tin=self.tin, tout=tout,
                      nrows=nrows, ncols=ncols,
                      GCj=GCj, GCi=GCi,
                      Real=data_array.real, Complex=data_array.imag)
            msg += dmi._get_complex_fields(func)
        else:
            dmi = DMI(self.name, matrix_form=self.form,
                      tin=self.tin, tout=tout,
                      nrows=nrows, ncols=ncols,
                      GCj=GCj, GCi=GCi,
                      Real=data_array, Complex=None)
            msg += dmi._get_real_fields(func)
        return msg


def sparse_symmetric_to_rectangular(name: str, matrix: sparse_types,
                                    log: SimpleLogger) -> np.ndarray:
    # get the upper and lower triangular matrices
    upper_tri = scipy.sparse.triu(matrix)
    lower_tri = scipy.sparse.tril(matrix)

    # extracts a [1, 2, 3, ..., n] off the diagonal of the matrix
    # and make it a diagonal matrix
    diagi = scipy.sparse.diags(scipy.sparse.diagional(upper_tri))

    # Check to see which triangle is populated.
    # If they both are, make sure they're equal
    # or average them and throw a warning
    lnnz = (lower_tri - diagi).nnz
    unnz = (upper_tri - diagi).nnz
    assert isinstance(lnnz, int), type(lnnz)
    assert isinstance(unnz, int), type(unnz)

    # both upper and lower triangle are populated
    if lnnz > 0 and unnz > 0:
        upper_tri_t = upper_tri.T
        if lower_tri == upper_tri_t:
            matrix = upper_tri + upper_tri_t - diagi
        else:
            log.warning(
                f'Matrix {name!r} marked as symmetric does not contain '
                'symmetric data.  Data will be symmetrized by averaging.')
            matrix = (matrix + matrix.T) / 2.
    elif lnnz > 0:
        #  lower triangle is populated
        matrix = lower_tri + lower_tri.T - diagi
    elif unnz > 0:
        #  upper triangle is populated
        matrix = upper_tri + upper_tri.T - diagi
    else:
        # matrix is diagonal (or null)
        matrix = diagi
    #data = matrix
    return matrix

def form_to_int(form: str) -> int:
    """gets the matrix description"""
    assert isinstance(form, str), form
    if form == 'N/A':
        return 0
    if form == 'square':
        return 1
    elif form == 'rectangular':
        return 2
    elif form == 'symmetric':
        return 6
    elif form == 'pseudo-identity':
        return 9
    raise RuntimeError(f'form = {form!r}')  # pragma: no cover
