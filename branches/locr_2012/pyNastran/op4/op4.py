#!/usr/bin/env python
import os
import sys
from struct import pack, unpack
from numpy import (array, zeros, float32, float64, complex64, complex128, 
                  allclose)
from scipy.sparse import coo_matrix
from pyNastran.general.general import is_binary
from pyNastran.general.mathematics import print_matrix, print_annotated_matrix
from pyNastran.op2.fortranFile import FortranFile


class OP4(FortranFile):
#class OP4(object):
    """
    @todo add endian checking
    @todo test on big matrices
    """
    def __init__(self):
        FortranFile.__init__(self)
        self.n = 0
        self.endian = ''

    def readOP4(self, op4Name, matrixNames=None, precision='default'):
        """
        Reads a NASTRAN OUTPUT4 file, regular or sparse, and stores the
        matrices as the output arguments of the function.  The number of matrices
        read is defined by the list matrixNames.  By default, all matrices will
        be read.  The resulting output is a dictionary of matrices that are
        accessed by their name.

		@code
        from pyNastran.op4.op4 import OP4
        op4 = OP4()
        #alternative way to get all the matrices
        matrices = op4.readOP4(op4Name)
        (formA,A) = matrices['A']
        (formB,B) = matrices['B']
        (formC,C) = matrices['C']

        # or to reduce memory usage
        matrices = op4.readOP4(op4Name,matrixNames=['A','B'])
        (formA,A) = matrices['A']
        (formB,B) = matrices['B']

        # or because you only want A
        matrices = op4.readOP4(op4Name,matrixNames='A')
        (formA,A) = matrices['A']
		@endcode

        @param op4Name an OP4 filename.  Type=STRING.
        @param matrixNames list of matrix names (or None); Type=LIST OF STRINGS / NONE.
        @param floatType specifies if the matrices are in single or double precsion
               (values='default','single','double') which means the format will be whatever the file is in

        @retval dictionary of matrices where the key is the name and the value is a matrix:
            Dense Type:  NUMPY.NDARRAY
            Sparse Type: SCIPY.SPARSE.COO_MATRIX

        @note based off the MATLAB code SAVEOP4 developed by ATA-E and later UCSD.
        @note it's strongly recommended that you convert sparse matrices to another
         format before doing math on them.  This is standard with sparse matrices.
        @warning sparse binary is buggy right now        """
        assert precision in ['default', 'single', 'double'], "precison=|%s| valid=['default','single','double']"
        if isinstance(matrixNames, str):
            matrixNames = [matrixNames]

        if not os.path.exists(op4Name):
            raise IOError('cannot find op4FileName=|%s|' % (op4Name))
        if is_binary(op4Name):
            return self.readOP4Binary(op4Name, matrixNames, precision)
        else:
            return self.readOP4Ascii(op4Name, matrixNames, precision)

#--------------------------------------------------------------------------
    def readOP4Ascii(self, op4Name, matrixNames=None, precision='default'):
        """matrixNames must be a list or None, but basically the same"""
        f = open(op4Name, 'r')
        matrices = {}
        name = 'dummyName'

        i = 0
        while name is not None:
            (name, form, matrix) = self.readMatrixAscii(
                f, matrixNames, precision)
            if name is not None:
                if matrixNames is None or name in matrixNames:  # save the matrix
                    matrices[name] = (form, matrix)
            i += 1
        return matrices

    def readMatrixAscii(self, f, matrixNames=None, precision='default'):
        """reads a matrix"""
        line = f.readline().rstrip()
        if line == '':
            f.close()
            return None, None, None
        ncols, nrows, form, Type = line[0:32].split()
        nrows = int(nrows)

        if nrows < 0:  # if less than 0, big
            isBigMat = True
        elif nrows > 0:
            isBigMat = False
        else:
            raise RuntimeError('unknown BIGMAT.  nRows=%s' % (nrows))

        nrows = abs(nrows)
        ncols = int(ncols)
        form = int(form)
        Type = int(Type)
        dType = self.getDType(Type, precision)

        name, size = line[32:].split()
        lineSize = size.split(
            ',')[1].split('E')[1].split('.')[0]  # 3E23.16 to 23
        lineSize = int(lineSize)

        line = f.readline().rstrip()
        (icol, irow, nWords) = line.split()
        irow = int(irow)

        isSparse = False
        if irow == 0:
            isSparse = True

        if Type in [1, 2]:  # real
            (A) = self.readRealAscii(f, nrows, ncols, lineSize,
                                     line, dType, isSparse, isBigMat)
        elif Type in [3, 4]:  # complex
            (A) = self.readComplexAscii(f, nrows, ncols,
                                        lineSize, line, dType, isSparse, isBigMat)
        else:
            raise RuntimeError('invalid matrix type.  Type=%s' % (Type))

        if not(matrixNames is None or name in matrixNames):  # kill the matrix
            A = None
        #print "form=%s name=%s A=\n%s" %(form,name,str(A))
        return (name, form, A)

    def readRealAscii(self, f, nrows, ncols, lineSize, line, dType, isSparse, isBigMat):
        """
        Method readRealAscii:
        @todo possibly split this into readDenseReal and readSparseReal
         to get rid of all the extra isSparse checks.  This would cleanup the
         runLoop condition as well.
        """
        if isSparse:
            rows = []
            cols = []
            entries = []
        else:
            A = zeros((nrows, ncols), dType)  # Initialize a real matrix

        nLoops = 0
        wasBroken = False
        while 1:
            if nLoops > 0 and not wasBroken:
                line = f.readline().rstrip()
            wasBroken = False

            (icol, irow, nWords) = line.split()
            icol = int(icol)

            if icol > ncols:
                break

            irow = int(irow)
            nWords = int(nWords)

            # This loop condition is overly complicated, but the first time
            # it will always execute.
            # Later if there is a sparse continuation line marker of
            # 1 (very large) integer, there will be no scientific notation value.
            # There also may be another sparse marker with 2 values.  These are not large.
            # The scientific check prevents you from getting stuck in an infinite
            # loop b/c no lines are read if there was one float value.
            # The check for 1 (or 2) integers is to prevent the check for 3 integers
            # which starts a new column.  We only want to continue a column.
            runLoop = True
            sline = line.strip().split()
            while (len(sline) == 1 or len(sline) == 2) and 'E' not in line or runLoop:  # next sparse entry
                irow = self.getIRowAscii(f, line, sline,
                                         nWords, irow, isSparse, isBigMat)

                runLoop = False
                i = 0
                iWord = 0
                isDoneReadingRow = False
                while nWords:
                    n = 0
                    line = f.readline().rstrip()
                    nWordsInLine = self.letterCount(line, 'E')
                    if nWordsInLine == 0:
                        wasBroken = True
                        break

                    for i in range(nWordsInLine):
                        word = line[n:n + lineSize]
                        if isSparse:
                            rows.append(irow - 1)
                            cols.append(icol - 1)
                            entries.append(word)
                        else:
                            A[irow - 1, icol - 1] = float(word)
                        n += lineSize
                        irow += 1
                        iWord += 1
                    nWords -= nWordsInLine
                sline = line.strip().split()
                nLoops += 1
            ###
        ###
        f.readline()

        if isSparse:  # Initialize a real matrix
            A = coo_matrix((entries, (rows, cols)), shape=(
                nrows, ncols), dtype=dType)
            #print "type = %s %s" %(type(A),type(A.todense()))
            #A = A.todense()
        return A

    def readComplexAscii(self, f, nrows, ncols, lineSize, line, dType, isSparse, isBigMat):
        """
        Method readComplexAscii:
        @todo possibly split this into readDenseComplex and readSparseComplex
         to get rid of all the extra isSparse checks.  This would cleanup the
         runLoop condition as well.
        """
        if isSparse:
            rows = []
            cols = []
            entries = []
        else:
            A = zeros((nrows, ncols), dType)  # Initialize a complex matrix

        nLoops = 0
        wasBroken = False
        while 1:
            if nLoops > 0 and not wasBroken:
                line = f.readline().rstrip()
            wasBroken = False

            (icol, irow, nWords) = line.split()
            icol = int(icol)

            if icol > ncols:
                break

            irow = int(irow)
            nWords = int(nWords)

            runLoop = True
            sline = line.strip().split()
            while (len(sline) == 1 or len(sline) == 2) and 'E' not in line or runLoop:  # next sparse entry
                irow = self.getIRowAscii(f, line, sline,
                                         nWords, irow, isSparse, isBigMat)
                runLoop = False

                i = 0
                iWord = 0
                isDoneReadingRow = False
                while nWords:
                    n = 0
                    line = f.readline().rstrip()
                    nWordsInLine = self.letterCount(line, 'E')
                    if nWordsInLine == 0:
                        wasBroken = True
                        break

                    for i in range(nWordsInLine):
                        value = float(line[n:n + lineSize])

                        if iWord % 2 == 0:
                            realValue = value
                        else:
                            if isSparse:
                                rows.append(irow - 1)
                                cols.append(icol - 1)
                                entries.append(complex(realValue, value))
                            else:
                                A[irow - 1, icol - 1] = complex(realValue,
                                                                value)
                            irow += 1
                        iWord += 1
                        n += lineSize
                    nWords -= nWordsInLine
                sline = line.strip().split()
                nLoops += 1
            ###
        ###
        if isSparse:  # Initialize a complex matrix
            A = coo_matrix((entries, (rows, cols)), shape=(
                nrows, ncols), dtype=dType)
        f.readline()
        return A

    def getIRowAscii(self, f, line, sline, nWords, irow, isSparse, isBigMat):
        if isSparse:
            #nWords = (nWords-1)//2  ## @todo this cant be right...
            sline = line.strip().split()
            if isBigMat:
                if len(sline) == 2:
                    pass
                else:
                    sline = f.readline().strip().split()
                assert len(sline) == 2, 'sline=%s len(sline)=%s' % (
                    sline, len(sline))
                (idummy, irow) = sline
                irow = int(irow)
            else:
                if len(sline) == 1:
                    IS = int(line)
                else:
                    IS = int(f.readline().strip())
                L = IS // 65536 - 1
                irow = IS - 65536 * (L + 1)
            ###
        ###
        return irow

    def letterCount(self, word, letter):
        """Counts the number of occurrences of a letter in a word/line."""
        n = 0
        for L in word:
            if L == letter:
                n += 1
        return n

#--------------------------------------------------------------------------
    def readOP4Binary(self, op4Name, matrixNames=None, floatType='default'):
        """matrixNames must be a list or None, but basically the same"""
        self.op4 = open(op4Name, 'rb')
        self.op2 = self.op4
        self.makeOp2Debug = False

        # get the endian
        data = self.op4.read(4)
        (recordLengthBig,) = unpack('>' + 'i', data)
        (recordLengthLittle,) = unpack('<' + 'i', data)
        #print
        if recordLengthBig == 24:
            self.endian = '>'
        elif recordLengthLittle == 24:
            self.endian = '<'
        else:
            msg = 'a 4 could not be found as the first word...endian error\n'
            msg += "RL_Big=%s RL_Little=%s" % (
                recordLengthBig, recordLengthLittle)
            raise RuntimeError(msg)
        self.op4.seek(0)

        i = 0
        matrices = {}
        name = 'dummyName'
        while name is not None:
            # checks for the end of the file
            n = self.n
            data1 = self.op4.read(1)
            self.op4.seek(n)
            if len(data1) == 0:
                break

            (name, form, matrix) = self.readMatrixBinary(
                self.op4, floatType, matrixNames)
            #print print_matrix(matrix)
            if name is not None:
                if matrixNames is None or name in matrixNames:  # save the matrix
                    matrices[name] = (form, matrix)

            #print "not f.closed = ",not self.op4.closed,form,name
            # if not self.op4.closed or form is not None:
            #     data = self.op4.read(4); self.n+=4
            #     if len(data)==0:
            #         break
            #     (recordLength,) = unpack(self.endian+'i',data)
            #     print("RL = %s" %(recordLength))
            #     if recordLength==24:
            #         self.n-=4; self.op4.seek(self.n)
            #     else:
            #         data = self.op4.read(4)
            #         if len(data)==0:
            #             break
            #         (recordLength2,) = unpack(self.endian+'i',data)
            #         assert recordLength2==24
            #         self.op4.seek(self.n)
            #
            i += 1
        return matrices

    def getDType(self, Type, precision='default'):
        """reset the type if 'default' not selected"""
        if precision == 'single':
            if Type in [1, 2]:
                dType = 'float32'
            else:
                dType = 'complex64'
        elif precision == 'double':
            if Type in [1, 2]:
                dType = 'float64'
            else:
                dType = 'complex128'
        else:  # default
            if Type == 1:
                dType = 'float32'
            elif Type == 2:
                dType = 'float64'
            elif Type == 3:
                dType = 'complex64'
            else:
                dType = 'complex128'
        return dType

    def readStartMarker(self, f):
        #print '--------------------------------------'
        #print self.printSection(60)
        data = f.read(4)
        self.n += 4
        (recordLength,) = unpack(self.endian + 'i', data)

        recordLength = 16
        data = f.read(recordLength)
        self.n += recordLength
        assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())

        if recordLength == 16:  # b,icol,irow,nWords,
            (a, icol, irow, nWords) = unpack(self.endian + '4i', data)
            #print "a=%s icol=%s irow=%s nWords=%s" %(a,icol,irow,nWords)
        else:
            raise NotImplementedError('recordLength=%s' % (recordLength))
        return (a, icol, irow, nWords)

    def getIRowSmall(self, f, data):
        if len(data) == 0:
            data = f.read(4)
            self.n += 4

        #print "$$$$$$$$ len(data) = ",len(data)
        IS, = unpack('i', data)
        L = IS // 65536 - 1
        irow = IS - 65536 * (L + 1)
        #print "IS=%s L=%s irow=%s" %(IS,L,irow)
        #assert IS>0
        #assert L>0
        return irow, L

    def getIRowBig(self, f, data):
        if len(data) == 0:
            data = f.read(8)
            self.n += 8
        (idummy, irow) = unpack('2i', data)
        #print "idummy=%s irow=%s" %(idummy,irow)
        #assert irow<100
        return (irow, idummy - 1)

    def readMatrixBinary(self, f, floatType, matrixNames=None):
        """reads a matrix"""
        #print self.printSection(60)
        #print "*************************"
        data = f.read(4)
        self.n += 4
        (recordLength,) = unpack(self.endian + 'i', data)
        assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
        #print "RL = %s" %(recordLength)

        #print self.printSection(60)
        if recordLength == 24:
            data = f.read(recordLength)
            self.n += recordLength

            (ncols, nrows, form, Type, name) = unpack(
                self.endian + '4i8s', data)
            #print "nrows=%s ncols=%s form=%s Type=%s name=%s" %(nrows,ncols,form,Type,name)
        else:
            msg = recordLength + self.printBlock(data)
            raise NotImplementedError('recordLength=%s\n%s' % (msg))

        name = name.strip()
        if 0:
            if Type == 1:
                print("Type = Real, Single Precision")
            elif Type == 2:
                print("Type = Real, Double Precision")
            elif Type == 3:
                print("Type = Complex, Single Precision")
            elif Type == 4:
                print("Type = Complex, Double Precision")

        if nrows < 0:  # if less than 0, big
            isBigMat = True
            nrows = abs(nrows)
        elif nrows > 0:
            isBigMat = False
        else:
            raise RuntimeError('unknown BIGMAT.  nRows=%s' % (nrows))

        # jump forward to get if isSparse, then jump back
        nSave = self.n
        (_a, _icol, _irow, _nWords) = self.readStartMarker(f)
        f.seek(nSave)
        self.n = nSave

        (NWV, NBW, d, dType) = self.getMatrixInfo(Type)

        isSparse = False
        if _irow == 0:
            isSparse = True
        #    rows = []
        #    cols = []
        #    entries = []
        #else:
        #    A = zeros((nrows,ncols),dType)
        #A = zeros((nrows,ncols),dType)

        assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
        if Type in [1, 2]:  # real
            A = self.readRealBinary(f, nrows, ncols, Type, isSparse, isBigMat)

        elif Type in [3, 4]:  # complex
            A = self.readComplexBinary(
                f, nrows, ncols, Type, isSparse, isBigMat)
        else:
            raise RuntimeError("Type=%s" % (Type))

        try:
            print_matrix(A.todense())
        except:
            pass

        if d in ['d', 'dd']:
            f.read(8)
            self.n += 8
        elif d in ['f', 'ff']:
            f.read(4)
            self.n += 4
        else:
            raise NotImplementedError(d)
        #f.read(recordLength); self.n+=recordLength
        #self.printSection(10)
        #f.read(4); self.n+=4

        #if isSparse:  # Initialize a real matrix
        #    A = coo_matrix( (entries,(rows,cols)),shape=(nrows,ncols),dtype=dType)
        #print '------end1--------'
        #print self.printSection(60)
        #print '------end2--------'
        return (name, form, A)

    def getMatrixInfo(self, Type):
        if Type == 1:
            NWV = 1  # number words per value
            NBW = 4
            d = 'f'
        elif Type == 2:
            NWV = 2
            NBW = 8
            d = 'd'
        elif Type == 3:
            NWV = 2
            NBW = 8
            d = 'ff'
        elif Type == 4:
            NWV = 4
            NBW = 16
            d = 'dd'
        else:
            raise RuntimeError("Type=%s" % (Type))
        dType = self.getDType(Type)
        return (NWV, NBW, d, dType)

    def readRealBinary(self, f, nrows, ncols, Type, isSparse, isBigMat):
        (NWV, NBW, d, dType) = self.getMatrixInfo(Type)
        if isSparse:
            rows = []
            cols = []
            entries = []
        else:
            A = zeros((nrows, ncols), dType)
        A = zeros((nrows, ncols), dType)

        data = ''
        icol = -1  # dummy value so the loop starts
        while icol < ncols + 1:  # if isDense
            #if recordLength==0:
            assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
            (icol, irow, nWords) = self.getMarkers(f, isSparse, isBigMat)
            L = nWords
            if 0:
                print("N=%s icol=%s irow=%s nWords=%s" % (
                    self.n, icol, irow, nWords))
                print("-----------")

            if icol == ncols + 1:
                #print "breaking***"
                break

            if isSparse:
                if isBigMat:
                    (irow, L) = self.getIRowBig(f, data[:8])
                    data = data[8:]
                else:
                    (irow, L) = self.getIRowSmall(f, data[:4])
                    data = data[4:]

            if L == -1:
                break

            if 0:
                print("next icol")
                print("N=%s icol=%s irow=%s nWords=%s" % (
                    self.n, icol, irow, nWords))

            #if nWords==0 and isBigMat:
                #self.n-=4; f.seek(self.n)
                #break

            recordLength = 4 * nWords
            data = f.read(recordLength)
            self.n += recordLength
            #print "dataFormat=%s RL=%s NNext=%s" %(d,recordLength,self.n)
            #if icol==ncols+1:
                #break

            nValues = nWords // NWV
            #nRead = nWords//4

            #print "recordLength=%s NBW=%s" %(recordLength,NBW)

            nValues2 = L // NWV
            nValues = nValues2
            #nRead2 = (L-1)

            if 0:
                print("inner while real...")
                print("nWords   = %s" % (nWords))
                print("nValues  = %s" % (nValues))
                print("nValues2 = %s" % (nValues2))
                print("NWV      = %s" % (NWV))
                assert nValues == nValues2

            #if nValues==0:
                #assert icol==ncols+1
                #break

            strValues = nValues * d
            if 0:
                print("strValues = %s" % (strValues))
                print("nValues*NBW=%s len(data)=%s" % (
                    nValues * NBW, len(data)))

            i = 0
            while len(data) > 0:
                if i == 0:
                    pass
                else:
                    if isSparse:
                        if isBigMat:
                            (irow, L) = self.getIRowBig(f, data[0:8])
                            data = data[8:]
                        else:
                            (irow, L) = self.getIRowSmall(f, data[0:4])
                            data = data[4:]
                    assert irow > 0
                    nValues2 = L // NWV
                    nValues = nValues2

                valueList = unpack(strValues, data[0:nValues * NBW])
                assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
                #print self.printSection(4)
                if 0:
                    print("valueList = %s" % (valueList))

                #irow-=1
                #icol-=1
                #print "isSparse = ",isSparse
                if isSparse:
                    #cols += [icol]*nValues
                    #rows += [i+irow for i in range(nValues)]
                    for value in valueList:
                        cols.append(icol - 1)
                        rows.append(irow - 1)
                        #if not allclose(value,1):
                            #asdf
                        #print "A[%s,%s] = %s" %(irow-1,icol-1,value)
                        entries.append(value)
                        A[irow - 1, icol - 1] = value
                        irow += 1
                else:
                    for value in valueList:
                        #print "A[%s,%s] = %s" %(irow-1,icol-1,value)
                        A[irow - 1, icol - 1] = value
                        irow += 1
                recordLength -= nValues * NBW
                data = data[nValues * NBW:]
                if 0:
                    print("recordLength=%s NBW=%s len(data)=%s" %
                          (recordLength, NBW, len(data)))
                    #print(A)
                    print("********")  # ,data
                    print(self.printBlock(data))
                i += 1
            #print "-------------------------------"

        if isSparse:  # Initialize a real matrix
            A = coo_matrix((entries, (rows, cols)), shape=(
                nrows, ncols), dtype=dType)

            if isBigMat:
                #f.read(4); self.n+=4
                #raise NotImplementedError()
                f.read(4)
                self.n += 4
                #pass
            else:
                f.read(4)
                self.n += 4
                #pass
        else:
            f.read(4)
            self.n += 4

        return A

    def readComplexBinary(self, f, nrows, ncols, Type, isSparse, isBigMat):
        (NWV, NBW, d, dType) = self.getMatrixInfo(Type)
        if isSparse:
            rows = []
            cols = []
            entries = []
        else:
            A = zeros((nrows, ncols), dType)
        A = zeros((nrows, ncols), dType)

        recordLength = 0
        data = ''
        icol = -1  # dummy value so the loop starts
        while icol < ncols + 1:  # if isDense
            #if recordLength==0:
            assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
            (icol, irow, nWords) = self.getMarkers(f, isSparse, isBigMat)
            if 0:
                print("N=%s icol=%s irow=%s nWords=%s" % (
                    self.n, icol, irow, nWords))
                print("-----------")

            L = nWords

            if icol == ncols + 1:
                break

            if isSparse:
                if isBigMat:
                    (irow, L) = self.getIRowBig(f, data[:8])
                    data = data[8:]
                else:
                    (irow, L) = self.getIRowSmall(f, data[:4])
                    data = data[4:]

            if L == -1:
                break

            if 0:
                print("next icol")
                print("N=%s icol=%s irow=%s nWords=%s" % (
                    self.n, icol, irow, nWords))

            #if nWords==0 and isBigMat:
                #self.n-=4; f.seek(self.n)
                #break

            recordLength = 4 * nWords
            data = f.read(recordLength)
            self.n += recordLength
            #print "dataFormat=%s RL=%s NNext=%s" %(d,recordLength,self.n)
            if icol == ncols + 1:
                continue

            nValues = nWords // NWV
            #nRead = nWords//4
            while recordLength >= NBW:
                if 0:
                    print("inner while...")
                    print("nWords  = %s" % (nWords))
                    print("nValues = %s" % (nValues))
                    print("NWV     = %s" % (NWV))

                #if nValues==0:
                    #assert icol==ncols+1
                    #break
                strValues = nValues * d
                if 0:
                    print("strValues = %s" % (strValues))
                    print("nValues*NBW=%s len(data)=%s" %
                          (nValues * NBW, len(data)))
                valueList = unpack(strValues, data[0:nValues * NBW])
                assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
                #print self.printSection(4)
                #print self.printBlock(data)
                if 0:
                    print("valueList = %s" % (valueList))

                #irow-=1
                #icol-=1
                #print "isSparse = ",isSparse
                irow -= 1
                icol -= 1

                if isSparse:
                    cols += [icol] * nValues
                    rows += [i + irow for i in range(nValues)]
                    for i, value in enumerate(valueList):
                        if i % 2 == 0:
                            realValue = value
                        else:
                            #A[irow,icol] = complex(realValue,value)
                            #print "A[%s,%s] = %s" %(irow,icol,complex(realValue,value))
                            A[irow, icol] = complex(realValue, value)
                            entries.append(complex(realValue, value))
                            irow += 1
                else:
                    for i, value in enumerate(valueList):
                        if i % 2 == 0:
                            realValue = value
                        else:
                            #print "A[%s,%s] = %s" %(irow,icol,complex(realValue,value))
                            A[irow, icol] = complex(realValue, value)
                            irow += 1

                recordLength -= nValues * NBW
                data = data[nValues * NBW:]
                #print "recordLength=%s NBW=%s" %(recordLength,NBW)
                #print print_matrix(A)
                #print "********",data

        if isSparse:  # Initialize a complex matrix
            A = coo_matrix((entries, (rows, cols)), shape=(
                nrows, ncols), dtype=dType)

            if isBigMat:
                #f.read(4); self.n+=4
                #raise NotImplementedError()
                f.read(4)
                self.n += 4
            else:
                f.read(4)
                self.n += 4
                #pass
        else:
            f.read(4)
            self.n += 4
        return A

    def getMarkers(self, f, isSparse, isBigMat):
        if isSparse:
            if isBigMat:
                (a, icol, irow, nWords) = self.readStartMarker(f)
                #(irow) = self.getIRowBig(f)
                nWords -= 2
                #if nWords>1:
                #    nWords -= 2
                #else:
                #    print("nWords0 = %s" %(nWords))
                #    nWords = 0
            else:
                (a, icol, irow, nWords) = self.readStartMarker(f)
                if irow != 0:
                    assert nWords == 1, 'nWords=%s' % (nWords)

                #(irow) = self.getIRowSmall(f)
                nWords -= 1
            ###
        else:
            (a, icol, irow, nWords) = self.readStartMarker(f)
            #print "N=%s a=%s icol=%s irow=%s nWords=%s"%(self.n,a,icol,irow,nWords)
        return (icol, irow, nWords)

#--------------------------------------------------------------------------
    def getTypeNWV(self, A, precision='default'):
        """
        @param A a matrix or entry in a matrix (to save memory)
        @param precision data precision ='default','single','double'
        @retval Type matrix type 1=real,single; 2=real,double; 3=complex,single; 4=complex,double
        @retval NWV Number of Words per Value
        """
        #print A.dtype.type()
        if isinstance(A.dtype.type(), float32):
            NWV = 1
            if precision != 'double':
                Type = 1
            else:
                Type = 2
        elif isinstance(A.dtype.type(), float64):
            NWV = 1
            if precision != 'single':
                Type = 2
            else:
                Type = 1

        # complex
        elif isinstance(A.dtype.type(), complex64):
            NWV = 2
            if precision != 'double':
                Type = 3
            else:
                Type = 4
        elif isinstance(A.dtype.type(), complex128):
            NWV = 2
            if precision != 'single':
                Type = 4
            else:
                Type = 3
        else:
            raise RuntimeError('invalid Type, only float32, float64, complex64, complex128')
        return (Type, NWV)

    def writeMatrixAscii(self, name, matrix, form=2, precision='default'):
        """
        Writes a real/complex matrix

        @param name the name of the matrix
        @param matrix a two-dimensional NUMPY.NDARRAY
        @param form Form is defined as one of the following:
        
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

        Not Supported by all OP4s:
        
        ==== ===============================
        Form Definition
        ==== ===============================
        4.   Lower triangular factor
        5.   Upper triangular factor
        10.  Cholesky factor
        11.  Trapezoidal factor
        13.  Sparse lower triangular factor
        15.  Sparse upper triangular factor
        ==== ===============================

        @note form defaults to 2, but 1 can be easily determined.  Any others must be specified.
        @todo call the actual function for now...not hooked up
        """
        assert isinstance(name, str), name
        assert isinstance(form, int), form

    def writeSparseMatrixAscii(self, f, name, matrix, form=2, isBigMat=False, precision='default', tol=1e-8):
        msg = ''
        assert isinstance(name, str), 'name=%s' % (name)
        #A = matrix.tolil() # list-of-lists sparse matrix
        A = matrix
        #print dir(matrix)
        (Type, NWV) = self.getTypeNWV(A.data[0], precision)
        if Type in [3, 4]:
            complexFactor = 2
        else:
            complexFactor = 1
        (nrows, ncols) = A.shape

        #if nrows==ncols and form==2:
        #    form = 1
        print("Type=%s" % (Type))
        if isBigMat:
            msg += '%8i%8i%8i%8i%-8s1P,3E23.16\n' % (
                ncols, -nrows, form, Type, name)
        else:
            msg += '%8i%8i%8i%8i%-8s1P,3E23.16\n' % (
                ncols, nrows, form, Type, name)

        #print "A.row = ",A.row
        #print "A.col = ",A.col

        cols = {}
        for j in A.col:
            cols[j] = []
        for i, jcol in enumerate(A.col):
            cols[jcol].append(i)
        #print "cols = ",cols

        f.write(msg)
        msg = ''
        for j, col in cols.iteritems():
            print("***********")
            print("j=%s col=%s" % (j, col))
            #col.sort()

            irows = [A.row[jj] for jj in col]
            #print "irows = ",irows
            (packs) = compressColumn(irows)
            print("packs = %s" % (packs))

            nPacks = len(packs)
            nRows = len(irows)
            if isBigMat:
                #L = complexFactor*(2*len(irows))+1
                L = 2 * nPacks * NWV + nRows
                msg = '%8i%8i%8i\n' % (j + 1, 0, L)
            else:
                L = complexFactor * (2 * len(irows))
                msg = '%8i%8i%8i\n' % (j + 1, 0, L + 1)
            f.write(msg)

            for (iPack, ipack) in enumerate(packs):
                msg = ''
                #print "pack = ",pack

                irow = A.row[col[pack[0]]]
                if isBigMat:
                    #L = complexFactor*(2*len(pack))+1
                    #L = (nPacks+1) + nRows*complexFactor
                    L = (len(ipack) + 1) * NWV
                    #if iPack==0:
                        #L+=1

                    #L = complexFactor*(2+nPacks)+1
                    #L = len(pack)+complexFactor*2
                    #msg = '%8i%8i%8i\n' %(j+1,0,L+1)
                    msg += '%8i%8i\n' % (L, irow + 1)
                else:
                    #L = complexFactor*(2*len(pack))
                    #msg = '%8i%8i%8i\n' %(j+1,0,L+1)

                    IS = irow + 65536 * (L + 1) + 1
                    msg += '%8i\n' % (IS)

                i = 0
                valueStr = ''
                print("ipack=%s rowPack=%s" % (ipack, [
                    A.row[p] for p in pack]))
                for p in pack:
                    irow = col[p]
                    val = A.data[irow]
                    irow = A.row[irow]

                    if Type in [1, 2]:
                        if abs(val) > tol:
                            valueStr += '%23.16E' % (val)
                        else:
                            valueStr += ' 0.0000000000000000E+00'
                        if (i + 1) % 3 == 0:
                            msg += valueStr + '\n'
                            #print "adding", valueStr
                            valueStr = ''
                    else:
                        if abs(val.real) > tol:
                            valueStr += '%23.16E' % (val.real)
                        else:
                            valueStr += ' 0.0000000000000000E+00'
                        if (i + 1) % 3 == 0:
                            msg += valueStr + '\n'
                            #print "adding", valueStr
                            valueStr = ''
                        i += 1
                        if abs(val.imag) > tol:
                            valueStr += '%23.16E' % (val.imag)
                        else:
                            valueStr += ' 0.0000000000000000E+00'
                        if (i + 1) % 3 == 0:
                            msg += valueStr + '\n'
                            #print "adding", valueStr
                            valueStr = ''
                    i += 1
                if valueStr:
                    msg += valueStr + '\n'
                f.write(msg)
        f.write('%8i%8i%8i\n' % (ncols + 1, 1, 1))
        f.write(' 1.0000000000000000E+00\n')
        ###

    def writeDenseMatrixBinary(self, name, matrix, form=2, precision='default', tol=1e-15):
        """
        24 is the record length
        """
        msg = ''
        A = matrix
        (Type, NWV) = self.getTypeNWV(A[0, 0], precision)

        (nrows, ncols) = A.shape
        #if nrows==ncols and form==2:
        #    form = 1

        msg += pack(self.endian + '5i8s', 24, ncols, nrows, form,
                    Type, '%-8s' % name)
        for icol in range(ncols):
            (iStart, iEnd) = self.getStartEndRow(A[:, icol], nrows, tol)

            # write the column
            if iStart is not None and iEnd is not None:
                iEnd += 1
                msg += pack(self.endian + '4i', 24, icol +
                            1, iStart + 1, (iEnd - iStart) * NWV)

                if Type in [1, 2]:
                    #nValues = iEnd-iStart+1
                    #msg += pack('d'*nValues,A[iStart:iEnd+1,icol])
                    for i, irow in enumerate(xrange(iStart, iEnd)):
                        msg += pack('d', A[irow, icol])

                else:  # complex
                    for irow in range(iStart, iEnd):
                        msg += pack('dd', A[irow, icol]
                                    .real, A[irow, icol].imag)
                    ###
                ###
            ###
        ###
        msg += pack(self.endian + '4id', 24, ncols + 1, 1, 1, 1.0)

        return msg

    def getStartEndRow(self, A, nrows, tol=1e-15):
        """find the starting and ending points of the matrix"""
        iStart = None
        for irow in range(nrows):
            if abs(A[irow]) > tol:
                iStart = irow
                break
        iEnd = None
        for irow in reversed(range(nrows)):
            if abs(A[irow]) > tol:
                iEnd = irow
                break
        return (iStart, iEnd)

    def writeDenseMatrixAscii(self, name, matrix, form=2, precision='default', tol=1e-15):
        msg = ''
        A = matrix
        (Type, NWV) = self.getTypeNWV(A[0, 0], precision)

        (nrows, ncols) = A.shape
        #if nrows==ncols and form==2:
        #    form = 1

        msg += '%8i%8i%8i%8i%-8s1P,3E23.16\n' % (
            ncols, nrows, form, Type, name)

        for icol in range(ncols):
            valueStr = ''
            (iStart, iEnd) = self.getStartEndRow(A[:, icol], nrows, tol)

            # write the column
            if iStart is not None and iEnd is not None:
                iEnd += 1
                msg += '%8i%8i%8i\n' % (icol + 1, iStart +
                                        1, (iEnd - iStart) * NWV)
                valueStr = ''
                #i=0
                if Type in [1, 2]:
                    #print "iStart=%s iEnd=%s" %(iStart,iEnd)
                    for i, irow in enumerate(xrange(iStart, iEnd)):
                        if abs(A[irow, icol]) > tol:
                            valueStr += '%23.16E' % (A[irow, icol])
                        else:
                            valueStr += ' 0.0000000000000000E+00'
                        if (i + 1) % 3 == 0:
                            msg += valueStr + '\n'
                            #print "adding", valueStr
                            valueStr = ''
                else:
                    i = 0
                    #print "iStart=%s iEnd=%s" %(iStart,iEnd)
                    for irow in range(iStart, iEnd):
                        if abs(A[irow, icol].real) > tol:
                            valueStr += '%23.16E' % (A[irow, icol].real)
                        else:
                            valueStr += ' 0.0000000000000000E+00'
                        if (i + 1) % 3 == 0:
                            msg += valueStr + '\n'
                            valueStr = ''
                        if abs(A[irow, icol].imag) > tol:
                            valueStr += '%23.16E' % (A[irow, icol].imag)
                        else:
                            valueStr += ' 0.0000000000000000E+00'
                        if (i + 2) % 3 == 0:
                            msg += valueStr + '\n'
                            valueStr = ''
                        i += 2
            if valueStr:
                msg += valueStr + '\n'
            ###
        ###
        msg += '%8i%8i%8i\n' % (ncols + 1, 1, 1)
        msg += ' 1.0000000000000000E+00\n'
        return msg


def matrices():
    strings = array([[ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 1   , 0   , 3   , 0   , 5   , 0   , 7   , 0   , 9   , 0   , 11  , 0   , 13  , 0   , 15  , 0   , 17  , 0   , 19  , 0    ],
                     [ 1   , 0   , 3   , 0   , 5   , 0   , 7   , 0   , 9   , 0   , 11  , 0   , 13  , 0   , 15  , 0   , 17  , 0   , 19  , 0    ],
                     [ 1   , 0   , 3   , 0   , 5   , 0   , 7   , 0   , 9   , 0   , 11  , 0   , 13  , 0   , 15  , 0   , 17  , 0   , 19  , 0    ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 1   , 0   , 3   , 0   , 5   , 0   , 7   , 0   , 9   , 0   , 11  , 0   , 13  , 0   , 15  , 0   , 17  , 0   , 19  , 0    ],
                     [ 1   , 0   , 3   , 0   , 5   , 0   , 7   , 0   , 9   , 0   , 11  , 0   , 13  , 0   , 15  , 0   , 17  , 0   , 19  , 0    ],
                     [ 1   , 0   , 3   , 0   , 5   , 0   , 7   , 0   , 9   , 0   , 11  , 0   , 13  , 0   , 15  , 0   , 17  , 0   , 19  , 0    ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18  , 19  , 20   ],
                     [ 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18  , 19  , 20   ],
                     [ 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18  , 19  , 20   ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18  , 19  , 20   ],
                     [ 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18  , 19  , 20   ],
                     [ 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18  , 19  , 20   ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ],
                     [ 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18  , 19  , 20   ],
                     [ 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18  , 19  , 20   ],
                     [ 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13  , 14  , 15  , 16  , 17  , 18  , 19  , 20   ],
                     [ 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0    ]],'f')
    return strings


def compressColumn(col):
    packs = []

    n = 0
    i = 0
    packi = []
    while i < len(col):
        #print "i=%s n=%s col[i]=%s" %(i,n,col[i])
        if col[i] == n + 1:
            #print "i=n=%s" %(i)
            packi.append(i)
            n += 1
        else:
            if packi:
                packs.append(packi)
                #print "pack = ",pack
            packi = [i]
            n = col[i]
        #print "pack = ",pack
        i += 1

    if packi:
        packs.append(packi)
    #print "packs = ",packs
    return (packs)

if __name__ == '__main__':
    #compressColumn([14, 15, 16, 20, 21, 22, 26, 27, 28])
    filenames = [
        #'test/mat_t_dn.op4',
        'test/mat_t_s1.op4',
        #'test/mat_t_s2.op4',
        #'test/mat_b_dn.op4',
        #'test/mat_b_s1.op4',
        #'test/mat_b_s2.op4',
        #'test/b_sample.op4',
        #'binary.op4',
    ]

    #matrixNames = 'EYE10' # identity
    #matrixNames = 'LOW'
    #matrixNames = 'RND1RS' # real,single
    #matrixNames = 'RND1RD' # real,double
    #matrixNames = 'RND1CS' # complex,single
    #matrixNames = 'RND1CD' # complex,double
    #matrixNames = 'STRINGS'
    #matrixNames = 'EYE5CD' # complex identity
    matrixNames = None
    strings = matrices()

    isBigMat = True
    f = open('ascii.op4', 'wb')
    for fname in filenames:
        op4 = OP4()
        op4.endian = '>'
        #if 't' in fname:
        #else:
            #f = open('binary.op4','wb')

        matrices = op4.readOP4(
            fname, matrixNames=matrixNames, precision='default')
        print("keys = %s" % (matrices.keys()))
        #print "#####################################################"
        print("fname=%s" % (fname))
        for name, (form, matrix) in sorted(matrices.items()):
            print("-----------------------------------")
            print("name = |%s|" % (name))
            if isinstance(matrix, coo_matrix):
                print("SPARSE")
                #matrix = matrix.todense()
                #print print_annotated_matrix(matrix)
            else:
                print("DENSE")
                print(print_matrix(matrix))

            #if 't' in fname:
            #f.write(op4.writeDenseMatrixAscii(name,matrix,form=form,precision='default'))
            if isinstance(matrix, coo_matrix):
                op4.writeSparseMatrixAscii(f, name, matrix, form=form,
                                           precision='default', isBigMat=isBigMat)
            else:
                f.write(op4.writeDenseMatrixAscii(name, matrix, 1, 'single'))
                #f.write(op4.writeDenseMatrixBinary(name,matrix,1,'single'))
        #print print_annotated_matrix(matrices['STRINGS'][1]-strings)
    print("-----------------------------")
    print("done")
    print("-----------------------------")
