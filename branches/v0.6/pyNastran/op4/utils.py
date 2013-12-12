from pyNastran.bdf.fieldWriter import print_card_8, print_card_16
from itertools import izip
from numpy import real, imag, ndarray
from scipy.sparse import coo_matrix

def write_DMIG(f, name, matrix, form, precision='default', isBigMat=False):
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
    """
    assert isinstance(matrix, coo_matrix), 'type(matrix)=%s' % type(matrix)
    #form = 6
    #assert form == 'cat', 'form=%r' % form
    assert precision in ['default', 'single', 'double'], "precison=%r valid=['default','single','double']"
    ifo = form

    #dtype = get_dtype(Type, precision='default').name
    dtype = matrix.dtype.name

    # get the matrix Type
    if dtype == 'float32':
        Type = 1
        size = 8
        is_complex = False
        print_card = print_card_8
    elif dtype == 'float64':
        Type = 2
        #size = 16
        is_complex = False
        print_card = print_card_16
    elif dtype == 'complex64':
        Type = 3
        #size = 8
        is_complex = True
        print_card = print_card_8
    elif dtype == 'complex128':
        Type = 4
        #size = 16
        is_complex = True
        print_card = print_card_16
    else:
        raise RuntimeError('dtype = %r' % dtype)

    tin = Type
    tout = Type

    # real/imaginary vs magnitude/phase
    polar = 0

    ncols = ''
    if ifo == 9: # Number of columns in a rectangular matrix. Used only for IFO = 9.
        nrows, ncols = A.shape
    card = ['DMIG', name, 0, form, precision, ifo, tin, tout, polar, '', ncols]
    f.write(print_card(card))

    index_to_node_id = {}
    index_to_component_id = {}

    for i in matrix.row:
        index_to_node_id[i] = i + 1
    for i in matrix.col:
        index_to_component_id[i] = 1 + i % 6

    row_index_to_node_id = index_to_node_id
    col_index_to_node_id = index_to_node_id

    row_index_to_component_id = index_to_component_id
    col_index_to_component_id = index_to_component_id

    if is_complex:
        reals = real(matrix.data)
        imags = imag(matrix.data)
        for i, j, reali, imagi in izip(matrix.row, matrix.col, reals, imags):
            GJ = col_index_to_node_id[i]
            CJ = col_index_to_node_id[i]
            G1 = row_index_to_node_id[i]
            C1 = row_index_to_node_id[i]
            card = ['DMIG', name, GJ, CJ, G1, C1, reali, imagi]
            f.write(print_card(card))
    else:
        for i, j, reali in izip(matrix.row, matrix.col, matrix.data):
            GJ = col_index_to_node_id[i]
            CJ = col_index_to_node_id[i]
            G1 = row_index_to_node_id[i]
            C1 = row_index_to_node_id[i]
            card = ['DMIG', name, GJ, CJ, G1, C1, reali]
            f.write(print_card(card))
    f.write('$-----------------------------------------------------------------------\n')
        #I = id_to_node_idj]