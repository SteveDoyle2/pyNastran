from pyNastran.bdf.fieldWriter import print_card_8 #, print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from itertools import izip
from numpy import real, imag, ndarray, where, arange
from numpy import matrix as Matrix
from scipy.sparse import coo_matrix


def write_DMIG(f, name, matrix, form,
               row_index_to_component_id, col_index_to_component_id,
               precision='default'):
    """
    ==== ===============
    Form Definition
    ==== ===============
    1.   Square
    2.   Rectangular
    3.   Diagonal
    6.   Symmetric
    8.   Id entity
    9.   Pseudoidentity
    ==== ===============

    ==== =========================
    Type Definition
    ==== =========================
    1.   Real, Single Precision
    2.   Real, Double Precision
    3    Complex, Single Precision
    4.   Complex, Double Precision
    ==== =========================

    :todo: collapse columns
    :todo: support ndarray
    """
    assert isinstance(matrix, coo_matrix), 'type(matrix)=%s' % type(matrix)
    matrix = matrix.todense()
    assert precision in ['default', 'single', 'double'], "precison=%r valid=['default','single','double']"
    ifo = form

    #dtype = get_dtype(Type, precision='default').name
    dtype = matrix.dtype.name

    # get the matrix Type
    if dtype == 'float32':
        precision_type = 1
        size = 8
        is_complex = False
        print_card = print_card_8
    elif dtype == 'float64':
        precision_type = 2
        #size = 16
        is_complex = False
        print_card = print_card_double
    elif dtype == 'complex64':
        precision_type = 3
        #size = 8
        is_complex = True
        print_card = print_card_8
    elif dtype == 'complex128':
        precision_type = 4
        #size = 16
        is_complex = True
        print_card = print_card_double
    else:
        raise RuntimeError('dtype = %r' % dtype)

    tin = precision_type
    tout = precision_type

    # real/imaginary vs magnitude/phase
    polar = 0

    ncols = ''
    form = ifo
    if ifo == 9: # Number of columns in a rectangular matrix. Used only for IFO = 9.
        nrows, ncols = A.shape
    f.write('$ type = %s\n' % type(matrix))
    card = ['DMIG', name, 0, ifo, tin, tout, polar, '', ncols]
    f.write(print_card_8(card))

    index_to_node_id = {}
    index_to_component_id = {}

    if isinstance(matrix, coo_matrix):
        if is_complex:
            reals = real(matrix.data)
            imags = imag(matrix.data)
            for i, j, reali, imagi in izip(matrix.row, matrix.col, reals, imags):
                GJ, CJ = col_index_to_node_component_id[i]
                G1, C1 = row_index_to_node_component_id[i]
                card = ['DMIG', name, GJ, CJ, G1, C1, reali, imagi]
                f.write(print_card(card))
        else:
            for i, j, reali in izip(matrix.row, matrix.col, matrix.data):
                GJ, CJ = col_index_to_node_component_id[i]
                G1, C1 = row_index_to_node_component_id[i]
                card = ['DMIG', name, GJ, CJ, G1, C1, reali]
                f.write(print_card(card))
    elif isinstance(matrix, Matrix):
        nrows, ncols = matrix.shape
        #if name == 'EYE5CD':
            #pass
        if is_complex:
            for icol in xrange(ncols):
                ii = where(matrix[:, icol] != 0.0)[0]
                reals = real(matrix[:, icol])
                imags = imag(matrix[:, icol])
                if ii.shape[1]: #len(ii) > 0:  # Matrix is 2D
                    #index_range = [ii.min(), ii.max()]
                    GJ, CJ = col_index_to_node_component_id[icol]
                    card = ['DMIG', name, GJ, CJ, '']
                    ii_max = ii.max()
                    ii_min = ii.min()
                    if ii_max == ii_min:
                        ii_max += 1
                    for iirow, irow in enumerate(arange(ii_min, ii_max)):
                        G1, C1 = row_index_to_node_component_id[irow]
                        card.append(G1)
                        card.append(C1)
                        card.append(reals[iirow, 0])
                        card.append(imags[iirow, 0])
                    f.write(print_card(card))
        else:
            for icol in xrange(ncols):
                ii = where(matrix[:, icol] != 0.0)[0]
                if ii.shape[1]: #len(ii) > 0:  # Matrix is 2D
                    #index_range = [ii.min(), ii.max()]
                    GJ, CJ = col_index_to_node_component_id[icol]
                    card = ['DMIG', name, GJ, CJ, '']
                    #for i, j, reali in izip(matrix.row, matrix.col, matrix.data):
                    ii_max = ii.max()
                    ii_min = ii.min()
                    if ii_max == ii_min:
                        ii_max += 1
                    for iirow, irow in enumerate(arange(ii_min, ii_max)):
                        G1, C1 = row_index_to_node_component_id[irow]
                        card.append(G1)
                        card.append(C1)
                        card.append(matrix[irow, icol])  # real number
                        card.append('')
                    f.write(print_card(card))
    else:
        raise NotImplementedError('type = %s' % type(matrix))
    f.write('$-----------------------------------------------------------------------\n')
        #I = id_to_node_idj]