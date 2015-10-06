# pylint: disable=R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import integer_types
from six.moves import zip, range
from math import sin, cos, radians, atan2, sqrt, degrees
#from math import (sin,sinh,cos,cosh,tan,tanh,sqrt,atan,atan2,acosh,acos,asin,
#                  asinh,atanh) #,atanh2   # going to be used by DEQATN

from numpy import array, zeros
from scipy.sparse import coo_matrix

from pyNastran.bdf.cards.baseCard import BaseCard

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double

from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, string, interpret_value)


class NastranMatrix(BaseCard):
    """
    Base class for the DMIG, DMIJ, DMIJI, DMIK matrices
    """
    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.name = string(card, 1, 'name')
            #zero

            #: 4-Lower Triangular; 5=Upper Triangular; 6=Symmetric; 8=Identity (m=nRows, n=m)
            self.ifo = integer(card, 3, 'ifo')
            #: 1-Real, Single Precision; 2=Real,Double Precision;
            #  3=Complex, Single; 4=Complex, Double
            self.tin = integer(card, 4, 'tin')
            #: 0-Set by cell precision
            self.tout = integer_or_blank(card, 5, 'tout', 0)

            #: Input format of Ai, Bi. (Integer=blank or 0 indicates real, imaginary format;
            #: Integer > 0 indicates amplitude, phase format.)
            self.polar = integer_or_blank(card, 6, 'polar', 0)
            if self.ifo in [6, 9]:
                self.ncol = integer(card, 8, 'ncol')
            else:  # technically right, but nulling this will fix bad decks
                self.ncol = None
                #self.ncol = blank(card, 8, 'ncol')
        else:
            raise NotImplementedError(data)

        self.GCj = []
        self.GCi = []
        self.Real = []
        if self.is_complex():
            self.Complex = []

    def write_code_aster(self):
        """
        assume set 1 = MAAX1,MAAX2, etc. and 100/n % on each
        """
        # for real combination
        comm = 'K_Mtx_AB=COMB_MATR_ASSE(COMB_R=(\n'
        comm += '    _F(MATR_ASSE = K_Mtx_A,COEF_R = 1.),\n'
        comm += '    _F(MATR_ASSE = K_Mtx_B,COEF_R = 1.)));\n'

        # for complex combination

        comm += "K_Mtx_AB=COMB_MATR_ASSE(COMB_C=(\n"
        comm += "_F(MATR_ASSE=K_Mtx_A,COEF_C=('RI',0.7,0.3,),)\n"
        comm += "_F(MATR_ASSE=K_Mtx_B,COEF_C=('RI',0.7,0.3,),),),);\n"
        comm = 'K_Mtx=ASSE_MATRICE(MATR_ELEM=ElMtx_K,NUME_DDL=%s,);'
        return comm

    def _add_column_uaccel(self, card=None, data=None, comment=''):
        raise NotImplementedError('UACCEL')

    def _add_column(self, card=None, data=None, comment=''):
        if comment:
            if hasattr(self, '_comment'):
                self._comment += comment
            else:
                self._comment = comment

        name = string(card, 1, 'name')
        if name == 'UACCEL':
            return  self._add_column_uaccel(card=None, data=None)

        Gj = integer(card, 2, 'Gj')
        # Cj = integer(card, 3, 'Cj')
        Cj = integer_or_blank(card, 3, 'Cj', 0)
        #Cj = components(card, 3, 'Cj')
        assert 0 <= Cj <= 6, 'C%i must be between [0, 6]; Cj=%s' % (0, Cj)

        nfields = len(card)
        # print("nfields = %i" % nfields)
        # print("card[5:] =", card[5:])
        # print("(nfields - 5) %% 4 = %i" % ((nfields - 5) % 4))

        nloops = (nfields - 5) // 4
        if (nfields - 5) % 4 in [2, 3]:  # real/complex
            nloops += 1
        #assert nfields <= 8,'nfields=%s' % nfields
        # print("nloops = %i" % nloops)
        assert nloops > 0, 'nloops=%s' % nloops

        for i in range(nloops):
            self.GCj.append((Gj, Cj))

        if self.is_complex():
            if self.is_polar():
                for i in range(nloops):
                    n = 5 + 4 * i
                    Gi = integer(card, n, 'Gi')
                    # Ci = integer(card, n + 1, 'Ci')
                    Ci = integer_or_blank(card, n + 1, 'Ci', 0)
                    #Ci = components(card, n + 1, 'Ci')
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
                    #Ci = components(card, n + 1, 'Ci')
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
                #Ci = components(card, n + 1, 'Ci')
                assert 0 <= Ci <= 6, 'C%i must be between [0, 6]; Ci=%s' % (i + 1, Ci)
                reali = double(card, n + 2, 'real')
                self.GCi.append((Gi, Ci))
                self.Real.append(reali)
                #print("GC=%s,%s real=%s" % (Gi, Ci, reali))

        msg = '(len(GCj)=%s len(GCi)=%s' % (len(self.GCj), len(self.GCi))
        assert len(self.GCj) == len(self.GCi), msg
        #if self.is_complex():
            #self.Complex(double(card, v, 'complex')

    def get_matrix(self, is_sparse=False, apply_symmetry=True):
        """
        Builds the Matrix

        Parameters
        ----------
        self : NastranMatrix
            the object pointer
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

    def isComplex(self):
        self.deprecated('isComplex()', 'is_complex()', '0.8')
        return self.is_complex()

    def isReal(self):
        self.deprecated('isReal()', 'is_real()', '0.8')
        return self.is_real()

    def isPolar(self):
        self.deprecated('isPolar()', 'is_polar()', '0.8')
        return self.is_polar()

    def getMatrix(self, isSparse=False, applySymmetry=True):
        self.deprecated('getMatrix()', 'get_matrix()', '0.8')
        return self.get_matrix(is_sparse=isSparse, apply_symmetry=applySymmetry)

    def is_real(self):
        return not self.is_complex()

    def is_complex(self):
        if self.tin in [1, 2]: # real
            return False
        elif self.tin in [3, 4]: # complex
            return True
        msg = ('Matrix %r must have a value of TIN = [1, 2, 3, 4].\n'
               'TIN defines the type (real, complex) of the matrix.  TIN=%r.' % (self.name, self.tin))
        raise ValueError(msg)

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
        msg = 'Matrix %r must have a value of POLAR = [0, 1].\n' % self.name
        msg += 'POLAR defines the type (real/imag or mag/phase) complex) of the matrix.  POLAR=%r.' % self.polar
        raise ValueError(msg)

    def getDType(self, type_flag):
        self.deprecated('getDType()', 'get_dtype()', '0.8')
        return self._get_dtype(type_flag)

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
            if self.is_complex():
                dtype = 'complex128'
            else:
                dtype = 'float64'
        else:
            raise RuntimeError("invalid option for matrix format")
        return dtype

    def __repr__(self):
        return self.write_card(size=8, is_double=False)

    def fill_in_default_C(self, model):
        for i, (Gi, Ci) in enumerate(self.GCi):
            if Ci is None:
                node = model.nodes[Gi]
                if node.type == 'GRID':
                    msg = 'Ci on DMIG card must be 1, 2, 3, 4, 5, or 6; Node=%i (GRID); Ci=%s' % (Gi, Ci)
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
                    msg = 'Cj on DMIG card must be 1, 2, 3, 4, 5, or 6; Node=%i (GRID); Cj=%s' % (Gj, Cj)
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
                       self.tout, self.polar, None, self.ncol]
        if size == 8:
            msg += print_card_8(list_fields)
        else:
            msg += print_card_16(list_fields)

        if self.is_complex():
            if self.is_polar():
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


def get_matrix(self, is_sparse=False, apply_symmetry=True):
    """
    Builds the Matrix

    :param self:     the object pointer
    :param is_sparse: should the matrix be returned as a sparse matrix (default=True).
                      Slower for dense matrices.
    :param apply_symmetry: If the matrix is symmetric (ifo=6), returns a symmetric matrix.
                           Supported as there are symmetric matrix routines.

    :returns M:    the matrix
    :returns rows: dictionary of keys=rowID,    values=(Grid,Component) for the matrix
    :returns cols: dictionary of keys=columnID, values=(Grid,Component) for the matrix
    .. warning:: is_sparse=True WILL fail
    """
    i = 0
    rows = {}
    rows_reversed = {}
    for GCi in self.GCi:
        if GCi not in rows:
            rows[GCi] = i
            rows_reversed[i] = GCi
            i += 1
    #nrows = len(rows2)

    j = 0
    cols = {}
    cols_reversed = {}
    for GCj in self.GCj:
        if GCj not in cols:
            cols[GCj] = j
            cols_reversed[j] = GCj
            j += 1
    #ncols = len(cols2)

    #A = ss.lil_matrix((3,3), dtype='d') # double precision

    #rows=[]; cols=[]; data=[]
    #for i in range(3):
    #    for j in range(3):
    #        k = float((i+1)*(j+1))
    #        rows.append(i)
    #        cols.append(j)
    #        data.append(k)
    #        A[i,j] = k

    #is_sparse = False
    if is_sparse:
        GCj = array(self.GCj, dtype='int32') - 1
        GCi = array(self.GCi, dtype='int32') - 1
        reals = array(self.Real, dtype='float32')

        # TODO: matrix size:  is this correct?
        nrows = max(GCi) + 1
        ncols = max(GCj) + 1

        dtype = self._get_dtype(self.tin)
        # TODO: no check for symmetry
        # TODO: no check for dtype
        if self.is_complex():
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
    else:
        if self.is_complex():
            M = zeros((i, j), dtype='complex128')
            if self.ifo == 6 and apply_symmetry:  # symmetric
                for (gcj, gci, reali, complexi) in zip(self.GCj, self.GCi, self.Real, self.Complex):
                    i = rows[gci]
                    j = cols[gcj]
                    M[i, j] = complex(reali, complexi)
                    M[j, i] = complex(reali, complexi)
            else:
                for (gcj, gci, reali, complexi) in zip(self.GCj, self.GCi, self.Real, self.Complex):
                    i = rows[gci]
                    j = cols[gcj]
                    M[i, j] = complex(reali, complexi)
        else:
            M = zeros((i, j), dtype='float64')
            if self.ifo == 6 and apply_symmetry:  # symmetric
                for (gcj, gci, reali) in zip(self.GCj, self.GCi, self.Real):
                    i = rows[gci]
                    j = cols[gcj]
                    M[i, j] = reali
                    M[j, i] = reali
            else:
                for (gcj, gci, reali) in zip(self.GCj, self.GCi, self.Real):
                    i = rows[gci]
                    j = cols[gcj]
                    M[i, j] = reali
    #print(M)
    return (M, rows_reversed, cols_reversed)



class DMIG(NastranMatrix):
    """
    Defines direct input matrices related to grid, extra, and/or scalar points.
    The matrix is defined by a single header entry and one or more column
    entries. A column entry is required for each column with nonzero elements.
    """
    type = 'DMIG'

    def __init__(self, card=None, data=None, comment=''):
        NastranMatrix.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            assert len(card) <= 9, 'len(DMIG card) = %i' % len(card)
        else:
            raise NotImplementedError(data)


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

    def __init__(self, card=None, data=None, comment=''):
        NastranMatrix.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            assert len(card) <= 9, 'len(DMIJ card) = %i' % len(card)
        else:
            raise NotImplementedError(data)


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

    def __init__(self, card=None, data=None, comment=''):
        NastranMatrix.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            assert len(card) <= 9, 'len(DMIJI card) = %i' % len(card)
        else:
            raise NotImplementedError(data)


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

    def __init__(self, card=None, data=None, comment=''):
        NastranMatrix.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            assert len(card) <= 9, 'len(DMIK card) = %i' % len(card)
        else:
            raise NotImplementedError(data)


class DMI(NastranMatrix):
    type = 'DMI'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.name = string(card, 1, 'name')
            #zero

            #: Form of the matrix:  1=Square (not symmetric); 2=Rectangular;
            #: 3=Diagonal (m=nRows,n=1);  4=Lower Triangular; 5=Upper Triangular;
            #: 6=Symmetric; 8=Identity (m=nRows, n=m)
            self.form = integer(card, 3, 'form')

            #: 1-Real, Single Precision; 2=Real,Double Precision;
            #: 3=Complex, Single; 4=Complex, Double
            self.tin = integer(card, 4, 'tin')

            #: 0-Set by cell precision
            self.tout = integer_or_blank(card, 5, 'tout', 0)

            self.nRows = integer(card, 7, 'nrows')
            self.nCols = integer(card, 8, 'ncols')

            assert len(card) == 9, 'len(DMI card) = %i' % len(card)
        else:
            raise NotImplementedError(data)

        self.GCj = []
        self.GCi = []
        self.Real = []

        if self.is_complex():
            self.Complex = []

    @property
    def ifo(self):
        return self.nRows == self.nCols

    def _add_column(self, card=None, data=None, comment=''):
        """
        .. todo :: support comment
        """
        if not self.is_complex():  # real
            self._read_real(card)

    def _read_real(self, card):
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
                        endI = fields[i + 1]
                        for ii in range(i1, endI + 1):
                            #print('adding j=%s i1=%s val=%s' % (j, ii, real_value))
                            self.GCj.append(j)
                            self.GCi.append(ii)
                            self.Real.append(real_value)
                        i += 1
                        is_done_reading_floats = True

    #def _read_complex(self, card):
        #msg = 'complex matrices not supported in the DMI reader...'
        #raise NotImplementedError(msg)
        ## column number
        #j = integer(card, 2, 'icol')

        ## counter
        #i = 0
        #fields = [interpret_value(field) for field in card[3:]]

        ## Complex, starts at A(i1,j)+imag*A(i1,j), goes to A(i2,j) in a column
        #while i < len(fields):
            #i1 = fields[i]
            #i += 1
            #isDoneReadingFloats = False
            #while not isDoneReadingFloats and i < len(fields):
                ##print("i=%s len(fields)=%s" %(i, len(fields)))
                #realValue = fields[i]
                #if isinstance(floatValue, integer_types):
                    #isDoneReadingFloats = True
                #elif isinstance(realValue, float):
                    #complexValue = fields[i + 1]
                    #self.GCj.append(j)
                    #self.GCi.append(i1)
                    #self.Real.append(realValue)
                    #self.Complex.append(complexValue)
                    #i += 2
                #else:
                      #raise NotImplementedError()

    def rename(self, new_name):
        self.name = new_name

    def is_real(self):
        return not self.is_complex()

    def is_complex(self):
        if self.tin in [3, 4]:
            return True
        return False

    def raw_fields(self):
        list_fields = ['DMI', self.name, 0, self.form, self.tin,
                       self.tout, None, self.nRows, self.nCols]

        if self.is_complex():
            for (gci, gcj, reali, imagi) in zip(self.GCi, self.GCj, self.Real, self.Complex):
                list_fields += ['DMI', self.name, gcj, gci, reali, imagi]
        else:
            for (gci, gcj, reali) in zip(self.GCi, self.GCj, self.Real):
                list_fields += ['DMI', self.name, gcj, gci, reali]
        return list_fields

    def write_card(self, size=8, is_double=False):
        msg = '\n$' + '-' * 80
        msg += '\n$ %s Matrix %s\n' % ('DMI', self.name)
        list_fields = ['DMI', self.name, 0, self.form, self.tin,
                       self.tout, None, self.nRows, self.nCols]
        if size == 8:
            msg += print_card_8(list_fields)
        else:
            msg += print_card_16(list_fields)

        if self.is_complex():
            for (gci, gcj, reali, imagi) in zip(self.GCi, self.GCj, self.Real, self.Complex):
                list_fields = ['DMI', self.name, gcj, gci, reali, imagi]
                if size == 8:
                    msg += print_card_8(list_fields)
                elif is_double:
                    msg += print_card_double(list_fields)
                else:
                    msg += print_card_16(list_fields)
        else:
            for (gci, gcj, reali) in zip(self.GCi, self.GCj, self.Real):
                list_fields = ['DMI', self.name, gcj, gci, reali]
                if size == 8:
                    msg += print_card_8(list_fields)
                elif is_double:
                    msg += print_card_double(list_fields)
                else:
                    msg += print_card_16(list_fields)
        return msg

    def __repr__(self):
        """
        .. todo:: support shortened output format.  There's a stupidly low 1000
                  DMI cap, I assume this is entries and not matrices.
        """
        return self.write_card(size=8, is_double=False)
