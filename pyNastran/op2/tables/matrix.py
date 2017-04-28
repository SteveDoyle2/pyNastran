"""Defines the Matrix class"""
from __future__ import print_function
from scipy.sparse import coo_matrix
import numpy as np

#from pyNastran.utils import object_attributes


class Matrix(object):
    """
    Defines a Matrix object

    Attributes
    ----------
    name : str
        the name of the matrix
    data : varies
        dense : np.ndarray
        sparse : coo_matrix
        data is initialized by setting the matrix.data attribute externally
    is_matpool : bool
        is this a matpool matrix
    """
    def __init__(self, name, form, is_matpool=False):
        """
        Initializes a Matrix

        Parameters
        ----------
        name : str
            the name of the matrix
        form : int
            the matrix type
        is_matpool : bool
            is this a matpool matrix

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
        self.data = None
        self.form = form
        self.is_matpool = is_matpool

        # only exist for is_matpool = True
        self.col_nid = None
        self.col_dof = None
        self.row_nid = None
        self.row_dof = None

    @property
    def shape_str(self):
        """gets the matrix description"""
        if self.form == 1:
            return 'square'
        elif self.form == 2:
            return 'rectangular'
        elif self.form == 6:
            return 'symmetric'
        elif self.form == 9:
            return 'pseudo-identity'
        else:
            raise RuntimeError('form = %s' % self.form)

    def write(self, mat, print_full=True):
        """writes to the F06"""
        mat.write(np.compat.asbytes(str(self) + '\n'))

        matrix = self.data
        if isinstance(matrix, coo_matrix):
            if print_full:
                for row, col, value in zip(matrix.row, matrix.col, matrix.data):
                    mat.write(np.compat.asbytes("(%i, %i) %s\n" % (row, col, value)))
            else:
                mat.write(str(matrix))
        else:
            mat.write(np.compat.asbytes('name=%r; shape=%s; form=%i; Type=%r\n' % (
                self.name, str(self.data.shape).replace('L', ''),
                self.form, self.shape_str)))
            if print_full:
                np.savetxt(mat, self.data, fmt='%.18e', delimiter=',')
            #f06.write(str(matrix))
            #print('WARNING: matrix type=%s does not support writing' % type(matrix))
        mat.write(np.compat.asbytes('\n\n'))

    def __repr__(self):
        class_name = str(type(self.data)).replace('<class ', '').replace('>', '').replace("'", '') + ';'
        header = 'Matrix[%r];' % self.name
        shape = ' shape=%s;' % str(self.data.shape).replace('L', '')
        dtype = '%s;' % self.data.dtype
        msg = '%-18s %-18s type=%-33s dtype=%-10s desc=%s' % (
            header, shape, class_name, dtype, self.shape_str)
        return msg

