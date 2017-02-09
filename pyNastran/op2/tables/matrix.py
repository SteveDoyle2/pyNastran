"""Defines the Matrix class"""
from __future__ import print_function
from scipy.sparse import coo_matrix

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
    def __init__(self, name, is_matpool=False):
        """
        Initializes a Matrix

        Parameters
        ----------
        name : str
            the name of the matrix
        is_matpool : bool
            is this a matpool matrix
        """
        self.name = name
        self.data = None
        self.is_matpool = is_matpool

    def write(self, f06, print_full=True):
        """writes to the F06"""
        f06.write(str(self) + '\n')

        matrix = self.data
        if isinstance(matrix, coo_matrix):
            if print_full:
                for row, col, value in zip(matrix.row, matrix.col, matrix.data):
                    f06.write("(%i, %i) %s\n" % (row, col, value))
            else:
                f06.write(str(matrix))
        else:
            f06.write(str(matrix))
            print('WARNING: matrix type=%s does not support writing' % type(matrix))
        f06.write('\n\n')

    def __repr__(self):
        class_name = str(type(self.data)).replace('<class ', '').replace('>', '').replace("'", '') + ';'
        header = 'Matrix[%r];' % self.name
        shape = ' shape=%s;' % str(self.data.shape).replace('L', '')
        msg = '%-18s %-18s type=%-33s dtype=%s' % (
            header, shape, class_name, self.data.dtype)
        return msg

