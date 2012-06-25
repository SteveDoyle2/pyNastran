import os
import sys

from numpy import array,zeros,ceil,sum,sign
from scipy.sparse import coo_matrix
from pyNastran.general.generalMath import printMatrix

class OP4(object):
    def __init__(self):
        pass

    def readOP4(self,op4Name,matrixNames=None,isAscii=True):
        """
        Reads a NASTRAN OUTPUT4 file, regular or sparse, and stores the
        matrices as the output arguments of the function.  The number of matrices
        read is defined by the list matrixNames.  By default, all matrices will
        be read.  The resulting output is a dictionary of matrices that are
        accessed by their name.

        # assume the file contains matrix A and B
        matrices = readOP4(op4Name,isAscii=True,matrixNames=['A','B'])
        A = matrices['A']
        B = matrices['B']

        @param op4Name an OP4 filename.  Type=STRING.
        @param isAscii the OP4 is an ASCII (human-readable file).  Type=BOOL.
        @param matrixNames list of matrix names (or None); Type=LIST OF STRINGS / NONE.
        @retval dictionary of matrices where the key is the name and the value is a matrix:
            Dense Type:  NUMPY.NDARRAY
            Sparse Type: SCIPY.SPARSE.COO_MATRIX

        @note based off the MATLAB code SAVEOP4 developed by ATA-E and later UCSD.
        @note it's strongly recommended that you convert sparse matrices to another
        format before performing math on them.  This is standard with sparse matrices.
        @warning isAscii=False is not supported.
        """
        assert isinstance(isAscii,bool),'isAscii must be a boolean.  isAscii=%s type=%s' %(isAscii,type(isAscii))
        assert isAscii==True,'Only isAscii=True supported; isAscii=%s' %(isAscii)
        if isinstance(matrixNames,str):
            matrixNames = [matrixNames]

        f = open(op4Name,'r')
        matrices = {}
        name = 'dummyName'

        i=0
        while name is not None:

            (name,matrix) = self.readMatrix(f,matrixNames)
            if name is not None:
                if matrixNames is None or name in matrixNames: # save the matrix
                    matrices[name] = matrix
            i+=1
        return matrices

    def readMatrix(self,f,matrixNames=None):
        """reads a matrix"""
        line = f.readline().rstrip()
        if line=='':
            f.close()
            return None,None
        ncols,nrows,form,Type = line[0:32].split()
        nrows = int(nrows)

        if nrows<0: # if less than 0, big
            isBigMat = True
        elif nrows>0:
            isBigMat = False
        else:
            raise RuntimeError('unknown BIGMAT.  nRows=%s' %(nrows))

        nrows = abs(nrows)
        ncols = int(ncols)
        form = int(form)
        Type = int(Type)

        name,size = line[32:].split()
        lineSize = size.split(',')[1].split('E')[1].split('.')[0] # 3E23.16 to 23
        lineSize = int(lineSize)

        line = f.readline().rstrip()
        (icol,irow,nWords) = line.split()
        irow = int(irow)

        isSparse = False
        if irow==0:
            isSparse = True

        if Type in [1,2]: # real
            (A) = self.readReal(f,nrows,ncols,lineSize,line,isSparse,isBigMat)
        elif Type in [3,4]: # complex
            (A) = self.readComplex(f,nrows,ncols,lineSize,line,isSparse,isBigMat)
        else:
            raise RuntimeError('invalid matrix type.  Type=%s' %(Type))

        if not(matrixNames is None or name in matrixNames): # kill the matrix
            A = None
        return name,A


    def readReal(self,f,nrows,ncols,lineSize,line,isSparse,isBigMat):
        """
        @todo possibly split this into readDenseReal and readSparseReal
        to get rid of all the extra isSparse checks.  This would cleanup the
        runLoop condition as well.
        """
        if isSparse:
            rows=[]; cols=[]; data=[]
        else:
            A = zeros((nrows,ncols),'f')       # Initialize a real matrix

        nLoops = 0
        wasBroken=False
        while 1:
            if nLoops>0 and not wasBroken:
                line = f.readline().rstrip()
            wasBroken = False

            (icol,irow,nWords) = line.split()
            icol = int(icol)

            if icol>ncols:
                break

            irow   = int(irow)
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
            while (len(sline)==1 or len(sline)==2) and 'E' not in line or runLoop: # next sparse entry
                irow = self.getIRow(f,line,sline,nWords,irow,isSparse,isBigMat)

                runLoop = False
                i=0
                iWord = 0
                isDoneReadingRow = False
                while nWords:
                    n = 0
                    line = f.readline().rstrip()
                    nWordsInLine = self.letterCount(line,'E')
                    if nWordsInLine==0:
                        wasBroken = True
                        break

                    for i in range(nWordsInLine):
                        word = line[n:n+lineSize]
                        if isSparse:
                            rows.append(irow-1)
                            cols.append(icol-1)
                            data.append(word)
                        else:
                            A[irow-1,icol-1] = float(word)
                        n += lineSize
                        irow+=1
                        iWord +=1
                    nWords-=nWordsInLine
                sline = line.strip().split()
                nLoops+=1
            ###
        ###
        f.readline()

        if isSparse:
            A = coo_matrix( (data,(rows,cols)),shape=(nrows,ncols),dtype='f') # Initialize a real matrix
            #print "type = %s %s" %(type(A),type(A.todense()))
            #A = A.todense()
        return A

    def readComplex(self,f,nrows,ncols,lineSize,line,isSparse,isBigMat):
        """
        @todo possibly split this into readDenseComplex and readSparseComplex
        to get rid of all the extra isSparse checks.  This would cleanup the
        runLoop condition as well.
        """
        if isSparse:
            rows=[]; cols=[]; data=[]
        else:
            A = zeros((nrows,ncols),'complex') # Initialize a complex matrix

        nLoops = 0
        wasBroken=False
        while 1:
            if nLoops>0 and not wasBroken:
                line = f.readline().rstrip()
            wasBroken = False

            (icol,irow,nWords) = line.split()
            icol = int(icol)

            if icol>ncols:
                break

            irow   = int(irow)
            nWords = int(nWords)

            runLoop = True
            sline = line.strip().split()
            while (len(sline)==1 or len(sline)==2) and 'E' not in line or runLoop: # next sparse entry
                irow = self.getIRow(f,line,sline,nWords,irow,isSparse,isBigMat)
                runLoop = False

                i=0
                iWord = 0
                isDoneReadingRow = False
                while nWords:
                    n = 0
                    line = f.readline().rstrip()
                    nWordsInLine = self.letterCount(line,'E')
                    if nWordsInLine==0:
                        wasBroken = True
                        break

                    for i in range(nWordsInLine):
                        word = float(line[n:n+lineSize])

                        if iWord%2==0:
                            realValue = word
                        else:
                            if isSparse:
                                rows.append(irow-1)
                                cols.append(icol-1)
                                data.append(realValue+word*1j)
                            else:
                                A[irow-1,icol-1] = realValue+word*1j
                            irow +=1
                        iWord +=1
                        n += lineSize
                    nWords-=nWordsInLine
                sline = line.strip().split()
                nLoops+=1
            ###
        ###
        if isSparse:
            A = coo_matrix( (data,(rows,cols)),shape=(nrows,ncols),dtype='complex') # Initialize a real matrix
            #print "type = %s %s" %(type(A),type(A.todense()))
            #A = A.todense()
        f.readline()
        return A

    def getIRow(self,f,line,sline,nWords,irow,isSparse,isBigMat):
        if isSparse:
            #nWords = (nWords-1)//2  ## @todo this cant be right...
            sline = line.strip().split()
            if isBigMat:
                if len(sline)==2:
                    pass
                else:
                    sline = f.readline().strip().split()
                assert len(sline)==2,'sline=%s len(sline)=%s' %(sline,len(sline))
                (idummy,irow) = sline
                irow = int(irow)
            else:
                if len(sline)==1:
                    IS = int(line)
                else:
                    IS = int(f.readline().strip())
                L = IS//65536 - 1
                irow = IS - 65536*(L + 1)
            ###
        ###
        return irow

    def letterCount(self,word,letter):
        """Counts the number of occurrences of a letter in a word/line."""
        n=0
        for L in word:
            if L==letter:
                n+=1
        return n

if __name__=='__main__':
    filename = 'mat_t_dn.op4' # works
    #filename = 'mat_t_s1.op4' # works
    #filename = 'mat_t_s2.op4' # works
    matrixNames = 'STRINGS'
    #matrixNames = None
    op4 = OP4()
    matrices = op4.readOP4(filename,matrixNames=matrixNames)
    for name,matrix in sorted(matrices.items()):
        print "name = %s" %(name)
        print matrix
    #print dir(matrices['RND1CS'])

