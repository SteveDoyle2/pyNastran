# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import izip
from math import log
#from math import (sin,sinh,cos,cosh,tan,tanh,sqrt,atan,atan2,acosh,acos,asin,
#                  asinh,atanh) #,atanh2   # going to be used by DEQATN

from numpy import zeros  # average
from scipy.sparse import coo_matrix

from pyNastran.bdf.cards.baseCard import BaseCard, print_card
from pyNastran.bdf.assign_type import (integer, integer_or_blank,
    double, string, blank, components)

def ssq(*listA):
    """
    sum of squares
    @note used for DEQATN
    """
    out = 0.
    for x in listA:
        out += x * x
    return out


def sum2(*listA):
    """
    sum of listA
    @note used for DEQATN
    """
    return sum(listA)


def mod(x, y):
    """
    x%y
    @note used for DEQATN
    """
    return x % y


def logx(x, y):
    """
    log base x of y
    @note used for DEQATN
    """
    log(y, x)


def dim(x, y):
    """
    @note used for DEQATN
    """
    return x - min(x, y)


def db(p, pref):
    """
    sound pressure in decibels
    would capitalize it, but you wouldnt be able to call the function...
    """
    return 20. * log(p / pref)


class DEQATN(BaseCard):  # needs work...
    type = 'DEQATN'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        newCard = ''
        foundNone = False
        for field in card.card:
            if foundNone is False and field is not None:
                newCard += field + ','
                foundNone = True
            elif foundNone is True and field is not None:
                newCard += field

        #if len(card.card)>1:
            #print "card.card = ",card.card
            #line0 = ','.join(card.card)
        #else:
            #line0 = ''.join(card.card)
        line0 = newCard
        self.eqID = line0[8:16]

        assert len(self.eqID) == 8, 'len(eqID)==%s' % (len(self.eqID))
        eq = line0[16:]
        eq = eq.replace(' ', '').lower()
        (self.name, self.eq) = eq.split('=')
        #print("EQ = %s" %(self.eq))

    def evaluate(self, args):
        #eqLow = self.eq.lower()
        #eval(self.eq)
        pass

    def __repr__(self):
        eq = self.name + '=' + self.eq
        eqLine = eq[0:56]
        eq = eq[56:]
        list_fields = ['DEQATN  ', '%8s' % (self.eqID), eqLine]

        if len(eq):
            eqLine = eq[0:72]
            eq = eq[72:]
            list_fields += ['        ' + eqLine]
        return ''.join(list_fields)


class NastranMatrix(BaseCard):
    """
    Base class for the DMIG, DMIJ, DMIJI, DMIK matrices
    """
    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        self.name = string(card, 1, 'name')
        #zero

        ## 4-Lower Triangular; 5=Upper Triangular; 6=Symmetric; 8=Identity (m=nRows, n=m)
        self.ifo = integer(card, 3, 'ifo')
        ## 1-Real, Single Precision; 2=Real,Double Precision; 3=Complex, Single; 4=Complex, Double
        self.tin = integer(card, 4, 'tin')
        ## 0-Set by cell precision
        self.tout = integer_or_blank(card, 5, 'tout', 0)

        self.polar = integer_or_blank(card, 6, 'polar')
        if self.ifo == 9:
            self.ncol = integer(card, 8, 'ncol')
        else:
            self.ncol = blank(card, 8, 'ncol')

        self.GCj = []
        self.GCi = []
        self.Real = []
        if self.isComplex():
            self.Complex = []

    def writeCodeAster(self):
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

    def addColumn(self, card=None, data=None, comment=''):
        if comment:
            if hasattr(self, '_comment'):
                self._comment += comment
            else:
                self._comment = comment
        Gj = integer(card, 2, 'Gj')
        Cj = components(card, 3, 'Cj')

        nfields = len(card)
        nloops = (nfields - 5) // 4
        minLoops = nloops - 1
        if minLoops <= 0:
            minLoops = 1
        #assert nFields <= 8,'nFields=%s' %(nFields)

        #print "minLoops = ",minLoops
        #print "nloops   = ",nloops
        for i in xrange(minLoops):
            self.GCj.append((Gj, Cj))

        if self.isComplex():
            for i in xrange(minLoops):
                n = 5 + 4 * i
                Gi = integer(card, n, 'Gi')
                Ci = components(card, n + 1, 'Ci')
                self.GCi.append((Gi, Ci))
                self.Real.append(double(card, n + 2, 'real'))
                self.Complex.append(double(card, n + 3, 'complex'))
        else:
            for i in xrange(minLoops):
                n = 5 + 4 * i
                Gi = integer(card, n, 'Gi')
                Ci = components(card, n + 1, 'Ci')
                self.GCi.append((Gi, Ci))
                self.Real.append(double(card, n + 2, 'real'))

        assert len(self.GCj) == len(self.GCi), '(len(GCj)=%s len(GCi)=%s' % (
            len(self.GCj), len(self.GCi))
        #if self.isComplex():
            #self.Complex(double(card, v, 'complex')

    def getMatrix(self, isSparse=False):
        """
        builds the Matrix
        @param self the object pointer
        @param isSparse should the matrix be returned as a sparse matrix (default=True).  Slower for dense matrices.
        @retval M the matrix
        @retval rows dictionary of keys=rowID,    values=(Grid,Component) for the matrix
        @retval cols dictionary of keys=columnID, values=(Grid,Component) for the matrix
        @warning isSparse WILL fail
        """
        i = 0
        rows = {}
        rowsReversed = {}
        for GCi in self.GCi:
            if GCi not in rows:
                rows[GCi] = i
                rowsReversed[i] = GCi
                i += 1
        #nRows = len(rows2)

        j = 0
        cols = {}
        colsReversed = {}
        for GCj in self.GCj:
            if GCj not in cols:
                cols[GCj] = j
                colsReversed[j] = GCj
                j += 1
        #nCols = len(cols2)

        #A = ss.lil_matrix((3,3), dtype='d') # double precision

        #rows=[]; cols=[]; data=[]
        #for i in xrange(3):
        #    for j in xrange(3):
        #        k = float((i+1)*(j+1))
        #        rows.append(i)
        #        cols.append(j)
        #        data.append(k)
        #        A[i,j] = k

        #isSparse = False
        if isSparse:
            data = []
            rows2 = []
            cols2 = []

            nrows = 0
            ncols = 0
            if self.isComplex():
                for (GCj, GCi, reali, complexi) in izip(self.GCj, self.GCi, self.Real, self.Complex):
                    i = rows[GCi]
                    j = cols[GCj]
                    nrows = max(i, nrows)
                    ncols = max(j, ncols)
                    rows2.append(i)
                    cols2.append(j)
                    data.append(complex(reali, complexi))
            else:
                for (GCj, GCi) in izip(self.GCj, self.GCi):
                    i = rows[GCi]
                    j = cols[GCj]
                    nrows = max(i, nrows)
                    ncols = max(j, ncols)
                    rows2.append(i)
                    cols2.append(j)
                data = self.Real
            #print("i=%s j=%s len(rows2)=%s len(cols2)=%s len(data)=%s"
            #    % (i, j, len(self.GCi), len(self.GCj), len(data)))
            # ,dtype=Format
            #print(rows2)

            #print("nrows=%s ncols=%s" %(nrows, ncols))
            if self.ifo in [1, 6]:
                nrows = max(nrows, ncols)
                ncols = nrows
            #print("nrows=%s ncols=%s" %(nrows, ncols))

            dType = self.getDType(self.tin)
            #A = coo_matrix( (entries,(rows,cols)),shape=(nrows,ncols),dtype=dType) # test
            M = coo_matrix((data, (self.GCi, self.GCj)),
                           shape=(nrows, ncols), dtype=dType)
            #M = coo_matrix( (data,(self.GCi,self.GCj)),shape=(i,j)) # old
            #M = coo_matrix( (data,(self.GCi,self.GCj)),shape=(nrows,ncols))
            #print M.todense()
            #print M
        else:
            if self.isComplex():
                M = zeros((i, j), dtype='complex128')
                for (GCj, GCi, reali, complexi) in izip(self.GCj, self.GCi, self.Real, self.Complex):
                    i = rows[GCi]
                    j = cols[GCj]
                    M[i, j] = complex(reali, complexi)
            else:
                M = zeros((i, j), dtype='float64')
                for (GCj, GCi, reali) in izip(self.GCj, self.GCi, self.Real):
                    i = rows[GCi]
                    j = cols[GCj]
                    M[i, j] = reali
        #print(M)
        return (M, rowsReversed, colsReversed)

    def rename(self, newName):
        self.name = newName

    def isReal(self):
        return not self.isComplex()

    def isComplex(self):
        if self.tin in [3, 4]:
            return True
        return False

    def getDType(self, Type):
        if Type == 1:
            dType = 'float32'
        elif Type == 2:
            dType = 'float64'
        elif Type == 3:
            dType = 'complex64'
        elif Type == 4:
            dType = 'complex128'
        elif Type == 0:
            if self.isComplex():
                dType = 'complex128'
            else:
                dType = 'float64'
        else:
            raise RuntimeError("invalid option for matrix format")
        return dType

    def __repr__(self):
        """
        @todo support double precision
        """
        msg = '\n$' + '-' * 80
        msg += '\n$ %s Matrix %s\n' % (self.type, self.name)
        list_fields = [self.type, self.name, 0, self.ifo, self.tin,
                  self.tout, self.polar, None, self.ncol]
        msg += print_card(list_fields)

        if self.isComplex():
            for (GCi, GCj, reali, imagi) in izip(self.GCi, self.GCj, self.Real, self.Complex):
                list_fields = [self.type, self.name, GCj[0], GCj[1],
                          None, GCi[0], GCi[1], reali, imagi]
                msg += print_card(list_fields)
        else:
            for (GCi, GCj, reali) in izip(self.GCi, self.GCj, self.Real):
                list_fields = [self.type, self.name, GCj[0], GCj[1],
                          None, GCi[0], GCi[1], reali, None]
                msg += print_card(list_fields)
        return msg


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


class DMI(BaseCard):
    type = 'DMI'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        self.name = string(card, 1, 'name')
        #zero

        ## Form of the matrix:  1=Square (not symmetric); 2=Rectangular;
        ## 3=Diagonal (m=nRows,n=1);  4=Lower Triangular; 5=Upper Triangular;
        ## 6=Symmetric; 8=Identity (m=nRows, n=m)
        self.form = integer(card, 3, 'form')

        ## 1-Real, Single Precision; 2=Real,Double Precision;
        ## 3=Complex, Single; 4=Complex, Double
        self.tin = integer(card, 4, 'tin')

        ## 0-Set by cell precision
        self.tout = integer_or_blank(card, 5, 'tout', 0)

        self.nRows = integer(card, 7, 'nrows')
        self.nCols = integer(card, 8, 'ncols')

        self.GCj = []
        self.GCi = []
        self.Real = []

        if self.isComplex():
            self.Complex = []

    def addColumn(self, card=None, data=None):

        if not self.isComplex():  # real
            self.readReal(card)

    def readReal(self, card):
        ## column number
        j = integer(card, 2, 'icol')

        # counter
        i = 0
        fields = card[3:]

        # Real, starts at A(i1,j), goes to A(i2,j) in a column
        while i < len(fields):
            i1 = fields[i]
            #print "i1 = ",i1
            if isinstance(i1, int):
                i += 1
                isDoneReadingFloats = False
                while not isDoneReadingFloats and i < len(fields):
                    #print "i=%s len(fields)=%s" %(i,len(fields))
                    realValue = fields[i]
                    if isinstance(realValue, int):
                        isDoneReadingFloats = True
                    elif isinstance(realValue, float):
                        self.GCj.append(j)
                        self.GCi.append(i1)
                        self.Real.append(realValue)
                        #print "i=%s j=%s value=%s" %(i1,j,realValue)
                        i += 1
                    else:
                        #print "*i=%s j=%s value=%s type=%s" %(i1,j,realValue,type(realValue))
                        realValue = self.Real[-1]
                        #print "*i=%s j=%s value=%s" %(i1,j,realValue)
                        endI = fields[i + 1]
                        #print "*i=%s endI=%s j=%s value=%s" %(i1,endI,j,realValue)
                        for ii in xrange(i1, endI + 1):
                            self.GCj.append(j)
                            self.GCi.append(ii)
                            self.Real.append(realValue)

                        #print "i = ",3+i
                        #print 'field i=',fields[i]
                        i += 1
                        isDoneReadingFloats = True

    # def readComplex(self, card):
    #     msg = 'complex matrices not supported in the DMI reader...'
    #     raise NotImplementedError(msg)
    #     ## column number
    #     j = integer(card, 2, 'icol')
    #
    #     # counter
    #     i = 0
    #     fields = card[3:]
    #
    #     # Complex, starts at A(i1,j)+imag*A(i1,j), goes to A(i2,j) in a column
    #     while i < len(fields):
    #         i1 = fields[i]
    #         i += 1
    #         isDoneReadingFloats = False
    #         asdf
    #         while not isDoneReadingFloats and i < len(fields):
    #             #print("i=%s len(fields)=%s" %(i, len(fields)))
    #             realValue = fields[i]
    #             if isinstance(floatValue, int):
    #                 isDoneReadingFloats = True
    #             elif isinstance(realValue, float):
    #                 complexValue = fields[i + 1]
    #                 self.GCj.append(j)
    #                 self.GCi.append(i1)
    #                 self.Real.append(realValue)
    #                 self.Complex.append(complexValue)
    #                 i += 2
    #             else:
    #                 asdf

    def rename(self, newName):
        self.name = newName

    def isReal(self):
        return not self.isComplex()

    def isComplex(self):
        if self.tin in [3, 4]:
            return True
        return False

    def __repr__(self):
        """
        @todo support shortened output format.  There's a stupidly low 1000
        DMI cap, I assume this is entries and not matrices.
        @todo support double precision
        """
        msg = '\n$' + '-' * 80
        msg += '\n$ %s Matrix %s\n' % (self.type, self.name)
        list_fields = [self.type, self.name, 0, self.form, self.tin,
                  self.tout, None, self.nRows, self.nCols]
        msg += self.print_card(list_fields)
        #msg += self.print_card(list_fields,size=16,isD=False)

        if self.isComplex():
            for (GCi, GCj, reali, imagi) in izip(self.GCi, self.GCj, self.Real, self.Complex):
                list_fields = [self.type, self.name, GCj, GCi, reali, imagi]
                msg += self.print_card(list_fields)
        else:
            for (GCi, GCj, reali) in izip(self.GCi, self.GCj, self.Real):
                list_fields = [self.type, self.name, GCj, GCi, reali]
                msg += self.print_card(list_fields)
        return msg
