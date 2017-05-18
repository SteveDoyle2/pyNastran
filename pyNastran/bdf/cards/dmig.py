# pylint: disable=R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from math import sin, cos, radians, atan2, sqrt, degrees
from six import iteritems
from six.moves import zip, range

import numpy as np
from numpy import array, zeros
from scipy.sparse import coo_matrix

from pyNastran.utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, string, parse_components, interpret_value)


class DTI(BaseCard):
    """
    +-----+-------+-----+------+-------+--------+------+-------------+
    |  1  |   2   |  3  |   4  |   5   |    6   |   7  |       8     |
    +=====+=======+=====+======+=======+========+======+=============+
    | DTI | UNITS | "1" | MASS | FORCE | LENGTH | TIME |   STRESS    |
    +-----+-------+-----+------+-------+--------+------+-------------+
    MSC

    +-----+-------+-----+------+-------+--------+------+-------------+
    |  1  |   2   |  3  |   4  |   5   |    6   |   7  |       8     |
    +=====+=======+=====+======+=======+========+======+=============+
    | DTI | UNITS | "1" | MASS | FORCE | LENGTH | TIME | TEMPERATURE |
    +-----+-------+-----+------+-------+--------+------+-------------+
    NX
    """
    type = 'DTI'
    def __init__(self, name, fields, comment=''):
        self.name = name
        self.fields = fields

    @classmethod
    def add_card(cls, card, comment):
        """
        Adds a DTI card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        name = string(card, 1, 'name')
        if name == 'UNITS':
            integer(card, 2, '1')
            mass = string(card, 3, 'mass')
            force = string(card, 3, 'force')
            length = string(card, 3, 'length')
            time = string(card, 3, 'time')
            temp_stress = string(card, 3, 'stress/temperature')
            fields = {
                'mass' : mass,
                'force' : force,
                'length' : length,
                'time' : time,
                'temp_stress' : temp_stress
            }
        else:
            print(card)
            raise NotImplementedError(card)
        return DTI(name, fields, comment=comment)

    def raw_fields(self):
        if self.name ==  'UNITS':
            mass = self.fields['mass']
            force = self.fields['force']
            length = self.fields['length']
            time = self.fields['time']
            temp_stress = self.fields['temp_stress']
            list_fields = ['DTI', self.name, '1', mass, force, length, time, temp_stress]
        else:
            raise NotImplementedError('DTI name=%r' % self.name)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class NastranMatrix(BaseCard):
    """
    Base class for the DMIG, DMIJ, DMIJI, DMIK matrices
    """
    def __init__(self, name, ifo, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment=''):
        """
        Creates a NastranMatrix

        Parameters
        ----------
        name : str
            the name of the matrix
        ifo : int
            matrix shape
            4=Lower Triangular
            5=Upper Triangular
            6=Symmetric
            8=Identity (m=nRows, n=m)
        tin : int
            matrix input precision
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        tout : int
            matrix output precision
            0=same as tin
            1=Real, Single Precision
            2=Real, Double Precision
            3=Complex, Single Precision
            4=Complex, Double Precision
        polar : int; default=0
            Input format of Ai, Bi
            Integer=blank or 0 indicates real, imaginary format
            Integer > 0 indicates amplitude, phase format
        ncols : int
            ???
        GCj  : List[(node, dof)]???
            the jnode, jDOFs
        GCi  : List[(node, dof)]???
            the inode, iDOFs
        Real : List[float]???
            The real values
        Complex : List[float]???; default=None
            The complex values (if the matrix is complex)
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        if Complex is None:
            Complex = []
        self.name = name

        #: 4-Lower Triangular; 5=Upper Triangular; 6=Symmetric; 8=Identity (m=nRows, n=m)
        self.ifo = ifo

        #: 1-Real, Single Precision; 2=Real,Double Precision;
        #  3=Complex, Single; 4=Complex, Double
        self.tin = tin

        #: 0-Set by cell precision
        self.tout = tout

        #: Input format of Ai, Bi. (Integer=blank or 0 indicates real, imaginary format;
        #: Integer > 0 indicates amplitude, phase format.)
        self.polar = polar

        self.ncols = ncols
        self.GCj = GCj
        self.GCi = GCi

        self.Real = Real
        if self.is_complex:
            self.Complex = Complex
        assert isinstance(ifo, integer_types), 'ifo=%r type=%s' % (ifo, type(ifo))
        assert not isinstance(ifo, bool), 'ifo=%r type=%s' % (ifo, type(ifo))

    @classmethod
    def add_card(cls, card, comment=''):
        name = string(card, 1, 'name')
        #zero

        ifo = integer(card, 3, 'ifo')
        tin = integer(card, 4, 'tin')
        tout = integer_or_blank(card, 5, 'tout', 0)
        polar = integer_or_blank(card, 6, 'polar', 0)
        if ifo == 1: # square
            ncols = integer_or_blank(card, 8, 'ifo=%s; ncol' % ifo)
        elif ifo == 6: # symmetric
            ncols = integer_or_blank(card, 8, 'ifo=%s; ncol' % ifo)
        elif ifo in [2, 9]: # rectangular
            ncols = integer(card, 8, 'ifo=%s; ncol' % (ifo))
        else:
            # technically right, but nulling this will fix bad decks
            #self.ncols = blank(card, 8, 'ifo=%s; ncol' % self.ifo)
            raise NotImplementedError('ifo=%s is not supported' % ifo)

        GCj = []
        GCi = []
        Real = []
        Complex = []
        return cls(name, ifo, tin, tout, polar, ncols,
                   GCj, GCi, Real, Complex, comment=comment)

    @property
    def matrix_type(self):
        if not isinstance(self.ifo, integer_types):
            msg = 'ifo must be an integer; ifo=%r type=%s name=%s' % (
                self.ifo, type(self.ifo), self.name)
            raise TypeError(msg)
        if isinstance(self.ifo, bool):
            msg = 'ifo must not be a boolean; ifo=%r type=%s name=%s' % (
                self.ifo, type(self.ifo), self.name)
            raise TypeError(msg)

        if self.ifo == 1:
            matrix_type = 'square'
        elif self.ifo == 6:
            matrix_type = 'symmetric'
        elif self.ifo in [2, 9]:
            matrix_type = 'rectangular'
        else:
            # technically right, but nulling this will fix bad decks
            #self.ncols = blank(card, 8, 'ifo=%s; ncol' % self.ifo)
            raise NotImplementedError('self.ifo=%r is not supported' % self.ifo)
        return matrix_type

    def finalize(self):
        self.GCi = np.asarray(self.GCi)
        self.GCj = np.asarray(self.GCj)
        self.Real = np.asarray(self.Real)
        if self.is_complex:
            self.Complex = np.asarray(self.Complex)

    @property
    def shape(self):
        if self.ifo in [1, 6]: # square, symmetric
            if self.ncols is not None:
                shape = (self.ncols, self.ncols)
            else:
                nrows, ncols = get_row_col_map(
                    self.GCi, self.GCj, self.ifo)[:2]
                shape = (nrows, ncols)
        elif self.ifo in [2, 9]:
            raise NotImplementedError('need to pull the nrows after reading in everything')
            #shape = (self.ncols, self.ncols)
        else:
            raise NotImplementedError('ifo=%s' % self.ifo)
        return shape

    def _add_column_uaccel(self, comment=''):
        raise NotImplementedError('UACCEL')

    def _add_column(self, card, comment=''):
        if comment:
            if hasattr(self, '_comment'):
                self.comment += comment
            else:
                self.comment = comment

        name = string(card, 1, 'name')
        if name == 'UACCEL':
            return  self._add_column_uaccel()

        Gj = integer(card, 2, 'Gj')
        # Cj = integer(card, 3, 'Cj')
        Cj = integer_or_blank(card, 3, 'Cj', 0)
        #Cj = parse_components(card, 3, 'Cj')
        assert 0 <= Cj <= 6, 'C%i must be between [0, 6]; Cj=%s' % (0, Cj)

        nfields = len(card)
        #print("nfields = %i" % nfields)
        #print("card[5:] =", card[5:])
        #print("(nfields - 5) %% 4 = %i" % ((nfields - 5) % 4))

        nloops = (nfields - 5) // 4
        if (nfields - 5) % 4 in [2, 3]:  # real/complex
            nloops += 1
        #assert nfields <= 8,'nfields=%s' % nfields
        #print("nloops = %i" % nloops)
        assert nloops > 0, 'nloops=%s' % nloops

        for i in range(nloops):
            self.GCj.append((Gj, Cj))

        if self.is_complex:
            if self.is_polar:
                for i in range(nloops):
                    n = 5 + 4 * i
                    Gi = integer(card, n, 'Gi')
                    # Ci = integer(card, n + 1, 'Ci')
                    Ci = integer_or_blank(card, n + 1, 'Ci', 0)
                    #Ci = parse_components(card, n + 1, 'Ci')
                    assert 0 <= Ci <= 6, 'C%i must be between [0, 6]; Ci=%s' % (i + 1, Ci)
                    self.GCi.append((Gi, Ci))
                    magi = double(card, n + 2, 'ai')
                    phasei = double(card, n + 3, 'bi')
                    reali = magi * cos(radians(phasei))
                    complexi = magi * sin(radians(phasei))
                    self.Real.append(reali)
                    self.Complex.append(complexi)
            else:
                for i in range(nloops):
                    n = 5 + 4 * i
                    Gi = integer(card, n, 'Gi')
                    # Ci = integer(card, n + 1, 'Ci')
                    Ci = integer_or_blank(card, n + 1, 'Ci', 0)
                    #Ci = parse_components(card, n + 1, 'Ci')
                    assert 0 <= Ci <= 6, 'C%i must be between [0, 6]; Ci=%s' % (i + 1, Ci)
                    self.GCi.append((Gi, Ci))
                    reali = double(card, n + 2, 'real')
                    complexi = double(card, n + 3, 'complex')
                    self.Real.append(reali)
                    self.Complex.append(complexi)
        else:
            # real
            for i in range(nloops):
                n = 5 + 4 * i
                Gi = integer(card, n, 'Gi')
                # Ci = integer(card, n + 1, 'Ci')
                Ci = integer_or_blank(card, n + 1, 'Ci', 0)
                #Ci = parse_components(card, n + 1, 'Ci')
                assert 0 <= Ci <= 6, 'C%i must be between [0, 6]; Ci=%s' % (i + 1, Ci)
                reali = double(card, n + 2, 'real')
                self.GCi.append((Gi, Ci))
                self.Real.append(reali)
                #print("GC=%s,%s real=%s" % (Gi, Ci, reali))

        msg = '(len(GCj)=%s len(GCi)=%s' % (len(self.GCj), len(self.GCi))
        assert len(self.GCj) == len(self.GCi), msg
        #if self.is_complex:
            #self.Complex(double(card, v, 'complex')

    def get_matrix(self, is_sparse=False, apply_symmetry=True):
        """
        Builds the Matrix

        Parameters
        ----------
        is_sparse : bool; default=False
            should the matrix be returned as a sparse matrix.
            Slower for dense matrices.
        apply_symmetry : bool; default=True
            If the matrix is symmetric (ifo=6), returns a symmetric matrix.
            Supported as there are symmetric matrix routines.

        Returns
        -------
        M : numpy.ndarray or scipy.coomatrix
            the matrix
        rows : dict[int] = [int, int]
            dictionary of keys=rowID, values=(Grid,Component) for the matrix
        cols: dict[int] = [int, int]
            dictionary of keys=columnID, values=(Grid,Component) for the matrix

        .. warning:: is_sparse=True WILL fail
        """
        return get_matrix(self, is_sparse=is_sparse, apply_symmetry=apply_symmetry)

    def rename(self, new_name):
        self.name = new_name

    @property
    def is_real(self):
        return not self.is_complex

    @property
    def is_complex(self):
        if self.tin in [1, 2]: # real
            return False
        elif self.tin in [3, 4]: # complex
            return True
        msg = ('Matrix %r must have a value of TIN = [1, 2, 3, 4].\n'
               'TIN defines the type (real, complex) '
               'of the matrix.  TIN=%r.' % (self.name, self.tin))
        raise ValueError(msg)

    @property
    def is_polar(self):
        """
        Used by:
          - DMIG
          - DMIJ
          - DMIJI
          - DMIK

        Not used by:
          - DMI
          - DMIAX
          - DMIG, UACCEL
          - DMIGOUT
          - DMIGROT
        """
        if self.polar == 0: # real, imag
            return False
        elif self.polar == 1: # mag, phase
            return True
        msg = ('Matrix %r must have a value of POLAR = [0, 1].\n'
               'POLAR defines the type (real/imag or mag/phase) complex) '
               'of the matrix.  POLAR=%r.' % (self.name, self.polar))
        raise ValueError(msg)

    @property
    def tin_dtype(self):
        return self._get_dtype(self.tin)

    @property
    def tout_dtype(self):
        return self._get_dtype(self.tout)

    def _get_dtype(self, type_flag):
        if type_flag == 1:
            dtype = 'float32'
        elif type_flag == 2:
            dtype = 'float64'
        elif type_flag == 3:
            dtype = 'complex64'
        elif type_flag == 4:
            dtype = 'complex128'
        elif type_flag == 0:
            if self.is_complex:
                dtype = 'complex128'
            else:
                dtype = 'float64'
        else:
            raise RuntimeError("invalid option for matrix format")
        return dtype

    def __repr__(self):
        return self.write_card(size=8, is_double=False)

    def fill_in_default_components(self, model):
        for i, (Gi, Ci) in enumerate(self.GCi):
            if Ci is None:
                node = model.nodes[Gi]
                if node.type == 'GRID':
                    msg = ('Ci on DMIG card must be 1, 2, 3, 4, 5, or 6; '
                           'Node=%i (GRID); Ci=%s' % (Gi, Ci))
                    raise RuntimeError(msg)
                elif node.type in ['SPOINT', 'EPOINT']:
                    Ci = 0
                else:
                    raise NotImplementedError(node)
                self.GCi[i] = [Gi, Ci]

        for i, (Gj, Cj) in enumerate(self.GCj):
            if Cj is None:
                node = model.nodes[Gj]
                if node.type == 'GRID':
                    msg = ('Cj on DMIG card must be 1, 2, 3, 4, 5, or 6; '
                           'Node=%i (GRID); Cj=%s' % (Gj, Cj))
                    raise RuntimeError(msg)
                elif node.type in ['SPOINT', 'EPOINT']:
                    Cj = 0
                else:
                    raise NotImplementedError(node)
                self.GCj[i] = [Gj, Cj]
        return

    def write_card(self, size=8, is_double=False):
        """
        .. todo:: support double precision
        """
        msg = '\n$' + '-' * 80
        msg += '\n$ %s Matrix %s\n' % (self.type, self.name)
        list_fields = [self.type, self.name, 0, self.ifo, self.tin,
                       self.tout, self.polar, None, self.ncols]
        if size == 8:
            msg += print_card_8(list_fields)
        else:
            msg += print_card_16(list_fields)

        if self.is_complex:
            if self.is_polar:
                for (GCi, GCj, reali, complexi) in zip(self.GCi, self.GCj, self.Real, self.Complex):
                    magi = sqrt(reali**2 + complexi**2)
                    if reali == 0.0:
                        phasei = 0.0
                    else:
                        phasei = degrees(atan2(complexi, reali))
                    list_fields = [self.type, self.name, GCj[0], GCj[1],
                                   None, GCi[0], GCi[1], magi, phasei]
                    if size == 8:
                        msg += print_card_8(list_fields)
                    elif is_double:
                        msg += print_card_double(list_fields)
                    else:
                        msg += print_card_16(list_fields)
            else:
                for (GCi, GCj, reali, complexi) in zip(self.GCi, self.GCj, self.Real, self.Complex):
                    list_fields = [self.type, self.name, GCj[0], GCj[1],
                                   None, GCi[0], GCi[1], reali, complexi]
                    if size == 8:
                        msg += print_card_8(list_fields)
                    elif is_double:
                        msg += print_card_double(list_fields)
                    else:
                        msg += print_card_16(list_fields)
        else:
            for (GCi, GCj, reali) in zip(self.GCi, self.GCj, self.Real):
                list_fields = [self.type, self.name, GCj[0], GCj[1],
                               None, GCi[0], GCi[1], reali, None]
                if size == 8:
                    msg += print_card_8(list_fields)
                elif is_double:
                    msg += print_card_double(list_fields)
                else:
                    msg += print_card_16(list_fields)
        return msg

def get_row_col_map(GCi, GCj, ifo):
    rows = {}
    rows_reversed = {}

    cols = {}
    cols_reversed = {}
    ndim = len(GCi.shape)
    #print('ndim=%s' % ndim)
    #print('GCj=%s' % GCj)
    #print('GCi=%s' % GCi)
    if ndim == 1:
        i = 0
        #nrows = np.unique(GCi)
        #ncols = np.unique(GCj)
        for gci in GCi:
            if gci not in rows:
                rows[gci] = i
                rows_reversed[i] = gci
                i += 1

        if ifo == 6:
            # symmetric
            for gcj in GCj:
                if gcj not in rows:
                    #print('row.gcj = %s' % str(gcj))
                    rows[gcj] = i
                    rows_reversed[i] = gcj
                    i += 1
            cols = rows
            cols_reversed = rows_reversed
        else:
            j = 0
            for gcj in GCj:
                if gcj not in cols:
                    cols[gcj] = j
                    cols_reversed[j] = gcj
                    j += 1
    else:
        #print('i0=%s j0=%s' % (i, j))
        #nrows = len(GCi)
        #ncols = len(GCj)
        #rows_array = np.zeros((nrows, 2), dtype='int32')
        #cols_array = np.zeros((ncols, 2), dtype='int32')
        #for i, (nid, comp) in enumerate(GCi):
            ##print('i=%s nid=%s comp=%s nrows=%s rows_array.shape=%s' % (
                ##i, nid, comp, nrows, str(rows_array.shape)))
            #rows_array[i, :] = [nid, comp]
        #print('rows_array = \n%s' % rows_array)

        #for j, (nid, comp) in enumerate(GCj):
            #cols_array[j, :] = [nid, comp]
        #print('cols_array = \n%s' % cols_array)

        i = 0
        for (nid, comp) in GCi:
            gci = (nid, comp)
            if gci not in rows:
                #print('row.gci = %s' % str(gci))
                rows[gci] = i
                rows_reversed[i] = gci
                i += 1
        if ifo == 6:
            # symmetric
            for (nid, comp) in GCj:
                gcj = (nid, comp)
                if gcj not in rows:
                    #print('row.gcj = %s' % str(gcj))
                    rows[gcj] = i
                    rows_reversed[i] = gcj
                    i += 1
            cols = rows
            cols_reversed = rows_reversed
        else:
            j = 0
            for (nid, comp) in GCj:
                gcj = (nid, comp)
                if gcj not in cols:
                    #print('col.gcj = %s' % str(gcj))
                    cols[gcj] = j
                    cols_reversed[j] = gcj
                    j += 1

    nrows = len(rows)
    ncols = len(cols)
    assert nrows > 0, 'nrows=%s' % nrows
    assert ncols > 0, 'ncols=%s' % ncols
    return nrows, ncols, ndim, rows, cols, rows_reversed, cols_reversed

def _fill_sparse_matrix(self, nrows, ncols):
    """helper method for get_matrix"""
    GCj = array(self.GCj, dtype='int32') - 1
    GCi = array(self.GCi, dtype='int32') - 1
    reals = array(self.Real, dtype='float32')

    # TODO: matrix size:  is this correct?
    nrows = max(GCi) + 1
    ncols = max(GCj) + 1

    dtype = self._get_dtype(self.tin)
    # TODO: no check for symmetry
    # TODO: no check for dtype
    if self.is_complex:
        complexs = array(self.Complex, dtype='float32')
        data = array([reals, complexs]).astype(complex)
    else:
        data = reals

    if self.ifo in [1, 6]:
        nrows = max(nrows, ncols)
        ncols = nrows

    #A = coo_matrix( (entries,(rows,cols)),shape=(nrows,ncols),dtype=dtype) # test
    M = coo_matrix((data, (self.GCi, self.GCj)),
                   shape=(nrows, ncols), dtype=dtype)
    #M = coo_matrix( (data,(self.GCi,self.GCj)),shape=(i,j)) # old
    #M = coo_matrix( (data,(self.GCi,self.GCj)),shape=(nrows,ncols))
    #print(M.todense())
    #print(M)
    return M


def _fill_dense_rectangular_matrix(self, nrows, ncols, ndim, rows, cols, apply_symmetry):
    """helper method for get_matrix"""
    is_sparse = False
    if self.is_complex:
        M = zeros((nrows, ncols), dtype='complex128')
        if self.ifo == 6 and apply_symmetry:  # symmetric
            for (gcj, gci, reali, complexi) in zip(self.GCj, self.GCi,
                                                   self.Real, self.Complex):
                i = rows[(gci[0], gci[1])]
                j = cols[(gcj[0], gcj[1])]
                M[i, j] = complex(reali, complexi)
                M[j, i] = complex(reali, complexi)
        else:
            for (gcj, gci, reali, complexi) in zip(self.GCj, self.GCi,
                                                   self.Real, self.Complex):
                i = rows[(gci[0], gci[1])]
                j = cols[(gcj[0], gcj[1])]
                M[i, j] = complex(reali, complexi)
    else:
        M = zeros((nrows, ncols), dtype='float64')
        if self.ifo == 6 and apply_symmetry:  # symmetric
            try:
                for (gcj, gci, reali) in zip(self.GCj, self.GCi, self.Real):
                    i = rows[(gci[0], gci[1])]
                    j = cols[(gcj[0], gcj[1])]
                    M[i, j] = reali
                    M[j, i] = reali
            except IndexError:
                msg = ('name=%s ndim=%s i=%s j=%s matrix_type=%s '
                       'is_polar=%s ncols=%s M.shape=%s\n' % (
                           self.name, ndim, i, j, self.matrix_type,
                           self.is_polar, self.ncols, M.shape))
                msg += 'Rows:\n'
                for i, row in enumerate(rows):
                    msg += 'i=%s row=%s\n' % (i, row)
                raise RuntimeError(msg)
        else:
            try:
                for (gcj, gci, reali) in zip(self.GCj, self.GCi, self.Real):
                    i = rows[(gci[0], gci[1])]
                    j = cols[(gcj[0], gcj[1])]
                    M[i, j] = reali
            except KeyError:
                msg = ('name=%s ndim=%s gci=%s gcj=%s matrix_type=%s '
                       'is_polar=%s is_sparse=%s ncols=%s M.shape=%s\n\n' % (
                           self.name, ndim, str(gci), str(gcj), self.matrix_type,
                           self.is_polar, is_sparse, self.ncols, M.shape))

                gci2 = (gci[0], gci[1])
                gcj2 = (gcj[0], gcj[1])
                if gci2 in rows:
                    msg += 'gci/row_key=%s found\n' % str(gci2)
                else:
                    msg += 'gci/row_key=%s not found\n' % str(gci2)
                    msg += 'Rows:\n'
                    for i, row in enumerate(rows):
                        msg += '  i=%s row=%s\n' % (i, row)

                if gcj2 in cols:
                    msg += '\ngcj/col_key=%s found\n' % str(gcj2)
                else:
                    msg += '\ngcj/col_key=%s not found\n' % str(gcj2)
                    msg += 'Cols:\n'
                    for j, col in enumerate(cols):
                        msg += '  j=%s row=%s\n' % (j, col)

                msg += '\n'
                print(msg)

                raise KeyError(msg)
            except IndexError:
                msg = ('name=%s ndim=%s i=%s j=%s matrix_type=%s '
                       'is_polar=%s is_sparse=%s ncols=%s M.shape=%s\n' % (
                           self.name, ndim, i, j, self.matrix_type,
                           self.is_polar, is_sparse, self.ncols, M.shape))
                msg += 'Rows:\n'
                for i, row in enumerate(rows):
                    msg += '  i=%s row=%s\n' % (i, row)

                msg += '\nCols:\n'
                for j, row in enumerate(cols):
                    msg += '  j=%s row=%s\n' % (j, col)
                raise RuntimeError(msg)
    return M


def _fill_dense_column_matrix(self, nrows, ncols, ndim, rows, cols, apply_symmetry):
    """helper method for get_matrix"""
    is_sparse = False
    if self.is_complex:
        M = zeros((nrows, ncols), dtype='complex128')
        if self.ifo == 6 and apply_symmetry:  # symmetric
            assert nrows == ncols, 'nrows=%s ncols=%s' % (nrows, ncols)
            for (gcj, gci, reali, complexi) in zip(self.GCj, self.GCi,
                                                   self.Real, self.Complex):
                i = rows[gci]
                j = cols[gcj]
                M[i, j] = complex(reali, complexi)
                M[j, i] = complex(reali, complexi)
        else:
            for (gcj, gci, reali, complexi) in zip(self.GCj, self.GCi,
                                                   self.Real, self.Complex):
                i = rows[gci]
                j = cols[gcj]
    else:
        #print('nrows=%s ncols=%s' % (nrows, ncols))
        M = zeros((nrows, ncols), dtype='float64')
        if self.ifo == 6 and apply_symmetry:  # symmetric
            assert nrows == ncols, 'nrows=%s ncols=%s' % (nrows, ncols)
            for (gcj, gci, reali) in zip(self.GCj, self.GCi, self.Real):
                i = rows[gci]
                j = cols[gcj]
                M[i, j] = reali
                M[j, i] = reali
        else:
            try:
                for (gcj, gci, reali) in zip(self.GCj, self.GCi, self.Real):
                    i = rows[gci]
                    j = cols[gcj]
                    M[i, j] = reali
            except IndexError:
                msg = ('name=%s ndim=%s i=%s j=%s matrix_type=%s '
                       'is_polar=%s is_sparse=%s ncols=%s M.shape=%s\n' % (
                           self.name, ndim, i, j, self.matrix_type,
                           self.is_polar, is_sparse, self.ncols, M.shape))
                msg += 'Rows:\n'
                for i, row in enumerate(rows):
                    msg += '  i=%s row=%s\n' % (i, row)
                raise RuntimeError(msg)
    return M

def get_matrix(self, is_sparse=False, apply_symmetry=True):
    """
    Builds the Matrix

    Parameters
    ----------
    is_sparse : bool
        should the matrix be returned as a sparse matrix (default=True).
        Slower for dense matrices.
    apply_symmetry: bool
        If the matrix is symmetric (ifo=6), returns a symmetric matrix.
        Supported as there are symmetric matrix routines.
        TODO: unused...

    Returns
    -------
    M : ndarray
        the matrix
    rows : Dict[(nid, nid)] = float
        dictionary of keys=rowID,    values=(Grid,Component) for the matrix
    cols : Dict[](int, int)] = float
        dictionary of keys=columnID, values=(Grid,Component) for the matrix

    .. warning:: is_sparse=True WILL fail
    """
    nrows, ncols, ndim, rows, cols, rows_reversed, cols_reversed = get_row_col_map(
        self.GCi, self.GCj, self.ifo)
    #print('rows = ', rows)
    #print('cols = ', cols)
    #print('i=%s j=%s' % (i, j))
    #nrows = len(rows2)
    #ncols = len(cols2)

    #A = ss.lil_matrix((3,3), dtype='d') # double precision
    #rows = []
    #cols = []
    #data = []
    #for i in range(3):
       #for j in range(3):
           #k = float((i+1)*(j+1))
           #rows.append(i)
           #cols.append(j)
           #data.append(k)
           #A[i,j] = k

    #is_sparse = False
    if is_sparse:
        M = _fill_sparse_matrix(self, nrows, ncols)
    else:
        if ndim == 1:
            M = _fill_dense_column_matrix(self, nrows, ncols, ndim, rows, cols, apply_symmetry)
        else:
            M = _fill_dense_rectangular_matrix(self, nrows, ncols, ndim, rows, cols, apply_symmetry)

    #print(M)
    return (M, rows_reversed, cols_reversed)


class DMIG_UACCEL(BaseCard):
    """
    Direct Matrix Input of Enforced Static Acceleration
    Defines rigid body accelerations in the basic coordinate system.

    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    |   1  |   2    |  3  |  4  |  5  |  6  |  7  |   8   |       |
    +======+========+=====+=====+=====+=====+=====+=======+=======+
    | DMIG | UACCEL | "0" | "9" | TIN |     |     |       | NCOL  |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    | DMIG | UACCEL |  L  |     |     |  G1 | C1  |  X1   |       |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    |      |   G2   |  C2 | X2  |     |  G3 | C3  |  X3   |       |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+

    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    | DMIG | UACCEL |  0  |  9  |  1  |     |     |       |   4   |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    | DMIG | UACCEL |  2  |     |     |  2  |  3  | 386.4 |       |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    | DMIG | UACCEL |  3  |     |     |  2  |  4  |  3.0  |       |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    | DMIG | UACCEL |  4  |     |     |  2  |  6  |  1.0  |       |
    +------+--------+-----+-----+-----+-----+-----+-------+-------+
    """
    type = 'DMIG'
    name = 'UACCEL'
    def __init__(self, tin, ncol, load_sequences, comment=''):
        if comment:
            self.comment = comment
        self.tin = tin
        self.ncol = ncol
        self.load_sequences = load_sequences
        #print(str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DMIG,UACCEL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        tin = integer(card, 4, 'tin')
        ncol = integer_or_blank(card, 8, 'ncol')
        return DMIG_UACCEL(tin, ncol, load_sequences={}, comment=comment)

    def _add_column(self, card, comment=''):
        load_seq = integer(card, 2, 'load_seq')

        g1 = integer(card, 5, 'nid1')
        c1 = parse_components(card, 6, 'c1')
        x1 = double(card, 7, 'x1')
        assert len(card) <= 8, 'len=%s card=%s' % (len(card), card)

        gcx = [g1, c1, x1]
        self.load_sequences[load_seq] = [gcx]

    @staticmethod
    def finalize():
        pass

    def raw_fields(self):
        list_fields = [
            'DMI', 'UACCEL', 0, 9, self.tin, None, None, None, self.ncol
        ]
        for lseq, ncx in sorted(iteritems(self.load_sequences)):
            list_fields += [lseq, None, None]
            for ncxi in ncx:
                list_fields += ncxi
           #for (nid, comp, xi) in ncx:
        #print('list_fields= %s' % list_fields)
        self.write_card()
        return list_fields

    def write_card(self, size=8, is_double=False):
        return self.write_card_8()

    def write_card_8(self):
        return self._write_card(print_card_8)

    def _write_card(self, func):
        msg = '\n$' + '-' * 80
        msg += '\n$ DMIG Matrix UACCEL\n'
        list_fields = [
            'DMIG', 'UACCEL', 0, 9, self.tin, None, None, None, self.ncol,
        ]
        msg += print_card_8(list_fields)

        list_fields = ['DMIG', 'UACCEL']
        for lseq, ncx in sorted(iteritems(self.load_sequences)):
            list_fields += [lseq, None, None]
            for ncxi in ncx:
                list_fields += ncxi
        #print('list_fields= %s' % list_fields)
        msg += print_card_8(list_fields)
        #print(msg)
        #if self.is_complex:
            #msg += self._get_complex_fields(func)
        #else:
            #msg += self._get_real_fields(func)
        return msg


class DMIG(NastranMatrix):
    """
    Defines direct input matrices related to grid, extra, and/or scalar points.
    The matrix is defined by a single header entry and one or more column
    entries. A column entry is required for each column with nonzero elements.

    +------+------+----+-----+-----+------+-------+----+------+
    |   1  |  2   | 3  |  4  |  5  |   6  |   7   | 8  |  9   |
    +======+======+====+=====+=====+======+=======+====+======+
    | DMIG | NAME | 0  | IFO | TIN | TOUT | POLAR |    | NCOL |
    +------+------+----+-----+-----+------+-------+----+------+
    | DMIG | NAME | GJ | CJ  |     |  G1  |  C1   | A1 |  B1  |
    +------+------+----+-----+-----+------+-------+----+------+
    |      |  G2  | C2 | A2  |  B2 |      |       |    |      |
    +------+------+----+-----+-----+------+-------+----+------+
    """
    type = 'DMIG'

    def __init__(self, name, ifo, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment=''):
        NastranMatrix.__init__(self, name, ifo, tin, tout, polar, ncols,
                               GCj, GCi, Real, Complex, comment=comment)


class DMIJ(NastranMatrix):
    """
    Direct Matrix Input at js-Set of the Aerodynamic Mesh
    Defines direct input matrices related to collation degrees-of-freedom
    (js-set) of aerodynamic mesh points for CAERO1, CAERO3, CAERO4 and CAERO5
    and for the slender body elements of CAERO2. These include W2GJ, FA2J and
    input pressures and downwashes associated with AEPRESS and AEDW entries.
    The matrix is described by a single header entry and one or more column
    entries. A column entry is required for each column with nonzero elements.
    For entering data for the interference elements of a CAERO2, use DMIJI
    or DMI.
    """
    type = 'DMIJ'

    def __init__(self, name, ifo, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment=''):
        NastranMatrix.__init__(self, name, ifo, tin, tout, polar, ncols,
                               GCj, GCi, Real, Complex, comment=comment)


class DMIJI(NastranMatrix):
    """
    Direct Matrix Input at js-Set of the Interference Body
    Defines direct input matrices related to collation degrees-of-freedom
    (js-set) of aerodynamic mesh points for the interference elements of CAERO2.
    These include W2GJ, FA2J and input pressures and downwashes associated with
    AEPRESS and AEDW entries. The matrix is described by a single header entry
    and one or more column entries. A column entry is required for each column
    with nonzero elements.  For entering data for the slender elements of a
    CAERO2, or a CAERO1, 3, 4 or 5 use DMIJ or DMI.
    """
    type = 'DMIJI'

    def __init__(self, name, ifo, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment=''):
        NastranMatrix.__init__(self, name, ifo, tin, tout, polar, ncols,
                               GCj, GCi, Real, Complex, comment=comment)


class DMIK(NastranMatrix):
    """
    Direct Matrix Input at ks-Set of the Aerodynamic Mesh
    Defines direct input matrices related to physical (displacement)
    degrees-of-freedom (ks-set) of aerodynamic grid points. These include WKK,
    WTFACT and input forces associated with AEFORCE entries. The matrix is
    described by a single header entry and one or more column entries. A column
    entry is required for each column with nonzero elements.
    """
    type = 'DMIK'

    def __init__(self, name, ifo, tin, tout, polar, ncols,
                 GCj, GCi, Real, Complex=None, comment=''):
        NastranMatrix.__init__(self, name, ifo, tin, tout, polar, ncols,
                               GCj, GCi, Real, Complex, comment=comment)


class DMI(NastranMatrix):
    type = 'DMI'

    def __init__(self, name, form, tin, tout, nrows, ncols,
                 GCj, GCi, Real, Complex=None, comment=''):
        """
        +------+-------+------+------+---------+----------+-----------+-----------+------+
        |  1   |   2   |  3   |   4  |    5    |    6     |     7     | 8         |  9   |
        +======+=======+======+======+=========+==========+===========+===========+======+
        | DMI  |  NAME |  0   | FORM |   TIN   |   TOUT   |           |     M     |  N   |
        +------+-------+------+------+---------+----------+-----------+-----------+------+
        | DMI  |  NAME |  J   |  I1  | A(I1,J) |  A(I1,J) | A(I1+1,J) | A(I1+2,J) | etc. |
        +------+-------+------+------+---------+----------+-----------+-----------+------+
        |      |  I2   | etc. |      |         |          |           |           |      |
        +------+-------+------+------+---------+----------+-----------+-----------+------+
        """
        #NastranMatrix.__init__(self, name, ifo, tin, tout, polar, ncols,
                               #GCj, GCi, Real, Complex, comment='')
        if comment:
            self.comment = comment
        if Complex is None:
            Complex = []
        self.name = name
        self.form = form
        self.tin = tin
        self.tout = tout
        self.nrows = nrows
        self.ncols = ncols
        self.GCi = GCi
        self.GCj = GCj
        self.Real = Real
        if self.is_complex:
            self.Complex = Complex

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DMI card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        name = string(card, 1, 'name')
        #zero

        #: Form of the matrix:  1=Square (not symmetric); 2=Rectangular;
        #: 3=Diagonal (m=nRows,n=1);  4=Lower Triangular; 5=Upper Triangular;
        #: 6=Symmetric; 8=Identity (m=nRows, n=m)
        form = integer(card, 3, 'form')

        #: 1-Real, Single Precision; 2=Real,Double Precision;
        #: 3=Complex, Single; 4=Complex, Double
        tin = integer(card, 4, 'tin')

        #: 0-Set by cell precision
        tout = integer_or_blank(card, 5, 'tout', 0)

        nrows = integer(card, 7, 'nrows')
        ncols = integer(card, 8, 'ncols')

        assert len(card) == 9, 'len(DMI card) = %i\ncard=%s' % (len(card), card)

        GCj = []
        GCi = []
        Real = []
        Complex = []
        return DMI(name, form, tin, tout, nrows, ncols,
                   GCj, GCi, Real, Complex, comment=comment)

    def finalize(self):
        self.GCi = np.asarray(self.GCi)
        self.GCj = np.asarray(self.GCj)
        self.Real = np.asarray(self.Real)
        if self.is_complex:
            self.Complex = np.asarray(self.Complex)

    def is_polar(self):
        if self.tin in [1, 2]:
            is_polar = False
        elif self.tin in [3, 4]:
            is_polar = False # TODO: could be wrong...
        else:
            raise NotImplementedError('nrows=%s ncols=%s' % (self.nrows, self.ncols))
        return is_polar

    @property
    def shape(self):
        return (self.nrows, self.ncols)

    @property
    def ifo(self):
        """
        ifo
        #: 4-Lower Triangular; 5=Upper Triangular; 6=Symmetric; 8=Identity (m=nRows, n=m)

        #: Form of the matrix:  1=Square (not symmetric); 2=Rectangular;
        #: 3=Diagonal (m=nRows,n=1);  4=Lower Triangular; 5=Upper Triangular;
        #: 6=Symmetric; 8=Identity (m=nRows, n=m)
        self.form = integer(card, 3, 'form')
        """
        return self.form
        #if self.nrows == self.ncols:
            ## symmetric
            #ifo = 6
        ##elif self.nrows > 1 and self.ncols > 1:
            ##ifo = 2
        #else:
            #raise NotImplementedError('form=%r nrows=%s ncols=%s' % (self.form, self.nrows, self.ncols))
        #return ifo

    def _add_column(self, card, comment=''):
        """
        .. todo:: support comment
        """
        if self.is_complex:
            return self._read_complex(card)
        else:
            return self._read_real(card)

    def _read_real(self, card):
        """reads a real DMI column"""
        # column number
        j = integer(card, 2, 'icol')

        # counter
        i = 0
        fields = [interpret_value(field) for field in card[3:]]

        # Real, starts at A(i1,j), goes to A(i2,j) in a column
        while i < len(fields):
            i1 = fields[i]
            if isinstance(i1, integer_types):
                i += 1
                is_done_reading_floats = False
                while not is_done_reading_floats and i < len(fields):
                    real_value = fields[i]
                    if isinstance(real_value, integer_types):
                        is_done_reading_floats = True
                    elif isinstance(real_value, float):
                        #print('adding j=%s i1=%s val=%s' % (j, i1, real_value))
                        self.GCj.append(j)
                        self.GCi.append(i1)
                        self.Real.append(real_value)
                        i += 1
                        i1 += 1
                    else:
                        real_value = self.Real[-1]
                        end_i = fields[i + 1]
                        for ii in range(i1, end_i + 1):
                            #print('adding j=%s i1=%s val=%s' % (j, ii, real_value))
                            self.GCj.append(j)
                            self.GCi.append(ii)
                            self.Real.append(real_value)
                        i += 1
                        is_done_reading_floats = True

    def _read_complex(self, card):
        """reads a complex DMI column"""
        #msg = 'complex matrices not supported in the DMI reader...'
        #raise NotImplementedError(msg)
        # column number
        j = integer(card, 2, 'icol')
        # counter
        i = 0
        fields = [interpret_value(field) for field in card[3:]]
        # Complex, starts at A(i1,j)+imag*A(i1,j), goes to A(i2,j) in a column
        while i < len(fields):
            i1 = fields[i]
            assert isinstance(i1, int), card
            i += 1
            is_done_reading_floats = False
            while not is_done_reading_floats and i < len(fields):
                value = fields[i]
                #print("i=%s len(fields)=%s value=%s" % (
                    #i, len(fields), value))
                if isinstance(value, integer_types):
                    is_done_reading_floats = True
                elif isinstance(value, float):
                    complex_value = fields[i + 1]
                    assert isinstance(complex_value, float), card
                    self.GCj.append(j)
                    self.GCi.append(i1)
                    self.Real.append(value)
                    self.Complex.append(complex_value)
                    i += 2
                else:
                    raise NotImplementedError()

    def rename(self, new_name):
        self.name = new_name

    @property
    def is_real(self):
        return not self.is_complex

    @property
    def is_complex(self):
        if self.tin in [3, 4]:
            return True
        return False

    def raw_fields(self):
        """
        .. warning:: All the writers are bad because Nastran insists on
                      making columns a single DMI card.  This makes
                      writing a card much harder, so there are a lot of
                      NotImplementedErrors floating about.

                      This is an invalid method, but is not disabled
                      because it's currently needed for checking results
        """
        list_fields = ['DMI', self.name, 0, self.form, self.tin,
                       self.tout, None, self.nrows, self.ncols]

        if self.is_complex:
            for (gci, gcj, reali, imagi) in zip(self.GCi, self.GCj, self.Real, self.Complex):
                list_fields += ['DMI', self.name, gcj, gci, reali, imagi]
        else:
            for (gci, gcj, reali) in zip(self.GCi, self.GCj, self.Real):
                list_fields += ['DMI', self.name, gcj, gci, reali]
        return list_fields

    def write_card_8(self):
        return self._write_card(print_card_8)

    def _get_real_fields(self, func):
        msg = ''
        uGCj = np.unique(self.GCj)
        for gcj in uGCj:
            i = np.where(gcj == self.GCj)[0]
            gcis = self.GCi[i]
            reals = self.Real[i]
            isort = np.argsort(gcis)
            list_fields = ['DMI', self.name, gcj]

            # will always write the first one
            gci_last = -1
            for gci, real in zip(gcis[isort], reals[isort]):
                if gci == gci_last + 1:
                    pass
                else:
                    list_fields.append(gci)
                list_fields.append(real)
                gci_last = gci
            msg += func(list_fields)
        return msg

    def _get_complex_fields(self, func):
        msg = ''
        uGCj = np.unique(self.GCj)
        for gcj in uGCj:
            i = np.where(gcj == self.GCj)[0]
            gcis = self.GCi[i]
            reals = self.Real[i]
            complexs = self.Complex[i]
            isort = np.argsort(gcis)
            list_fields = ['DMI', self.name, gcj]

            # will always write the first one
            gci_last = -1
            #print('gcis=%s \nreals=%s \ncomplexs=%s' % (
                #gcis[isort], reals[isort], complexs[isort]))
            if max(gcis) == min(gcis):
                list_fields += [gcis[0]]
                for reali, complexi in zip(reals, complexs):
                    list_fields.extend([reali, complexi])
                msg += func(list_fields)
            else:
                for gci, real, complexi in zip(gcis[isort], reals[isort], complexs[isort]):
                    if gci != gci_last + 1:
                        pass
                    else:
                        list_fields.append(gci)
                    list_fields.append(real)
                    list_fields.append(complexi)
                    gci_last = gci
                msg += func(list_fields)
        return msg

    def write_card_16(self):
        return self._write_card(print_card_16)

    def write_card_double(self):
        return self._write_card(print_card_16)

    def _write_card(self, func):
        msg = '\n$' + '-' * 80
        msg += '\n$ %s Matrix %s\n' % ('DMI', self.name)
        list_fields = ['DMI', self.name, 0, self.form, self.tin,
                       self.tout, None, self.nrows, self.ncols]
        msg += print_card_8(list_fields)

        if self.is_complex:
            msg += self._get_complex_fields(func)
        else:
            msg += self._get_real_fields(func)
        return msg

    def write_card(self, size=8, is_double=False):
        if size == 8:
            return self.write_card_8()
        elif is_double:
            return self.write_card_double()
        return self.write_card_16()

    def __repr__(self):
        """
        .. todo:: support shortened output format.  There's a stupidly low 1000
                  DMI cap, I assume this is entries and not matrices.
        """
        return self.write_card(size=8, is_double=False)
