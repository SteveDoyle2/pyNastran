"""Defines the Matrix class"""
from scipy.sparse import coo_matrix  # type: ignore
import numpy as np
from pyNastran.op2.op2_interface.write_utils import export_to_hdf5
from pyNastran.utils import object_attributes, object_methods
from pyNastran.op2.op2_interface.op2_codes import MSC_ELEMENTS


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
        if not isinstance(name, str):
            raise TypeError('name=%r must be a string; type=%s' % (name, type(name)))

    @property
    def shape_str(self):
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
            raise RuntimeError('form = %r' % self.form)

    def export_to_hdf5(self, group, log):
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

    def write(self, mat, print_full=True):
        """writes to the F06"""
        mat.write(np.compat.asbytes(str(self) + '\n'))

        matrix = self.data
        if self.data is None:
            skip_msg = 'skipping %s because data is None\n\n' % self.name
            mat.write(skip_msg.encode('ascii'))
            return
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

    def object_attributes(self, mode='public', keys_to_skip=None,
                          filter_properties=False):
        if keys_to_skip is None:
            keys_to_skip = []

        my_keys_to_skip = [
            'object_methods', 'object_attributes',
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(self, mode='public', keys_to_skip=None):
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = []

        my_keys_to_skip = [
            'object_methods', 'object_attributes',
        ]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def __repr__(self):
        header = 'Matrix[%r];' % self.name
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


class MatrixDict:
    """storage object for KDICT, MDICT, BDICT, etc. is op2.matdicts"""
    def __init__(self, name):
        self.name = name
        self.element_types = []
        self.numwides = []
        self.numgrids = []
        self.dof_per_grids = []

        self.eids = []
        self.ge = []
        self.address = []
        self.forms = []
        self.sils = []
        self.xforms = []

    def add(self, eltype, numwids, numgrid, dof_per_grid, form,
            eids, ge, address, sil, xform=None):
        """Sets the next set of the KDICT"""
        self.element_types.append(eltype)
        self.numwides.append(numwids)
        self.numgrids.append(numgrid)
        self.dof_per_grids.append(dof_per_grid)
        self.forms.append(form)

        self.eids.append(eids)
        self.ge.append(ge)
        self.address.append(address)
        self.sils.append(sil)
        self.xforms.append(xform)

    #@property
    #def nodes(self):
        #return [sil // 10 for sil in self.sils]

    #@property
    #def dofs(self):
        #return [sil % 10 for sil in self.sils]

    @property
    def nelements(self):
        return sum([len(eids) for eids in self.eids])

    @property
    def element_names(self):
        return [MSC_ELEMENTS[etype] for etype in self.element_types]

    def __repr__(self):
        msg = 'MatrixDict(name=%r, nelements=%s element_types=%s, element_names=[%s])' % (
            self.name, self.nelements, self.element_types, ', '.join(self.element_names))
        return msg
