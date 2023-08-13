"""Defines the Matrix class"""
from typing import Union, Optional
from scipy.sparse import coo_matrix, csr_matrix  # type: ignore
import numpy as np
from pyNastran.op2.op2_interface.write_utils import export_to_hdf5
from pyNastran.utils import object_attributes, object_methods


class Matrix:
    """
    Defines a Matrix object that's stored in op2.matrices

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
    def __init__(self, name: str, form: Union[int, str],
                 data: Optional[Union[np.ndarray, coo_matrix]]=None):
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
        if isinstance(form, int):
            self.form = form
        else:
            self.form = form_to_int(form)
        self.is_matpool = False

        # only exist for is_matpool = True; automatically set
        self.col_nid = None
        self.col_dof = None
        self.row_nid = None
        self.row_dof = None
        if not isinstance(name, str):
            raise TypeError(f'name={name!r} must be a string; type={type(name)}')

    def set_matpool_data(self, data: np.ndarray,
                         col_nid: np.ndarray, col_dof: np.ndarray,
                         row_nid: np.ndarray, row_dof: np.ndarray) -> None:
        self.is_matpool = True
        self.data = data
        self.col_nid = col_nid
        self.col_dof = col_nid
        self.row_nid = row_nid
        self.row_dof = row_dof


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

    def export_to_hdf5(self, group, log) -> None:
        """exports the object to HDF5 format"""
        export_to_hdf5(self, group, log)

    def build_dataframe(self):
        """exports the object to pandas format"""
        import pandas as pd
        matrix = self.data
        if matrix is None:
            return
        if isinstance(matrix, coo_matrix):
            data = {'row': matrix.row, 'col': matrix.col, 'data' : matrix.data}
            data_frame = pd.DataFrame(data=data).reindex(columns=['row', 'col', 'data'])
        elif isinstance(matrix, np.ndarray):
            data_frame = pd.DataFrame(data=matrix)
        else:
            raise NotImplementedError(type(matrix))
        self.data_frame = data_frame

    def write(self, mat, print_full: bool=True) -> None:
        """writes to the F06"""
        mat.write(np.compat.asbytes(str(self) + '\n'))

        matrix = self.data
        if self.data is None:
            skip_msg = f'skipping {self.name!r} because data is None\n\n'
            mat.write(skip_msg.encode('ascii'))
            return
        if isinstance(matrix, coo_matrix):
            if print_full:
                for row, col, value in zip(matrix.row, matrix.col, matrix.data):
                    mat.write(np.compat.asbytes("(%d, %d) %s\n" % (row, col, value)))
            else:
                mat.write(str(matrix))
        else:
            mat.write(np.compat.asbytes('name=%r; shape=%s; form=%d; Type=%r\n' % (
                self.name, str(self.data.shape).replace('L', ''),
                self.form, self.shape_str)))
            if print_full:
                np.savetxt(mat, self.data, fmt='%.18e', delimiter=',')
            #f06.write(str(matrix))
            #print('WARNING: matrix type=%s does not support writing' % type(matrix))
        mat.write(np.compat.asbytes('\n\n'))

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
        assert isinstance(self.data, (np.ndarray, coo_matrix)), type(self.data)
        return self.data.dtype.name

    def to_gcj_gci_form(self):
        """once in matrix form, we need to transform to GCj, GCi, Real/Complex form"""
        nrows, ncols = self.data.shape

        # fix columns to consider
        is_real = self.dtype_str in {'float32', 'float64'}

        sparse_csr = csr_matrix(self.data)
        sparse_coo = sparse_csr.tocoo()  #  this should be the answer....


        max_cols = (np.abs(self.data.real) + np.abs(self.data.imag)).max(axis=1)
        assert ncols == len(maxs)
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


    #def write_to_bdf(self, dmig_type: str):
        #list_fields = [dmig_type, ]
        #pass


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
    raise RuntimeError(f'form = {form!r}')

