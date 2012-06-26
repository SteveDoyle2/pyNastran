import os
import sys
from struct import unpack
from numpy import zeros
from scipy.sparse import coo_matrix
from pyNastran.op2.fortranFile import FortranFile
from pyNastran.general.generalMath import printMatrix

class OP4(FortranFile):
    """
    @todo integrate with ASCII reader; get rid of isAscii
    @todo add endian checking
    @todo test on big matrices
    """
    def __init__(self):
        FortranFile.__init__(self)
        self.makeOp2Debug = False # required to make FortranFile work
        self.n = 0
        self.endian = ''

    def readOP4(self,op4Name,matrixNames=None,isAscii=False,floatType='default'):
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
        @param matrixNames list of matrix names (or None); Type=LIST OF STRINGS / NONE.
        @param floatType specifies if the matrices are in single or double precsion
               (values='default','single','double') which means the format will be whatever the file is in
         
        @param isAscii the OP4 is an ASCII (human-readable file).  Type=BOOL.
        @retval dictionary of matrices where the key is the name and the value is a matrix:
            Dense Type:  NUMPY.NDARRAY
            Sparse Type: SCIPY.SPARSE.COO_MATRIX

        @note based off the MATLAB code SAVEOP4 developed by ATA-E and later UCSD.
        @note it's strongly recommended that you convert sparse matrices to another
        format before performing math on them.  This is standard with sparse matrices.
        @warning isAscii=False is not supported.
        """
        assert floatType in ['default','single','double']
        assert isinstance(isAscii,bool),'isAscii must be a boolean.  isAscii=%s type=%s' %(isAscii,type(isAscii))
        assert isAscii==False,'Only isAscii=False supported; isAscii=%s' %(isAscii)
        if isinstance(matrixNames,str):
            matrixNames = [matrixNames]
        self.op2 = open(op4Name,'rb')

        matrices = {}
        name = 'dummyName'

        i=0
        while name is not None:
            # checks for the end of the file
            n = self.n
            data1 = self.op2.read(1)
            self.op2.seek(n)
            if len(data1)==0:
                break

            (name,matrix) = self.readMatrix(self.op2,floatType,matrixNames)
            if name is not None:
                if matrixNames is None or name in matrixNames: # save the matrix
                    matrices[name] = matrix
            i+=1
        return matrices

    def readStartMarker(self,f):
        data = f.read(4); self.n+=4
        (recordLength,) = unpack(self.endian+'i',data)

        recordLength = 16
        data = f.read(recordLength); self.n+=recordLength

        if recordLength==16: # b,icol,irow,nWords,
            (a,icol,irow,nWords) = unpack(self.endian+'4i',data)
        else:
            raise NotImplementedError('recordLength=%s' %(recordLength))
        return (a,icol,irow,nWords)

    def getIRowSmall(self,f):
        data = f.read(4); self.n+=4
        IS, = unpack('i',data)
        L = IS//65536 - 1
        irow = IS - 65536*(L + 1)
        return irow

    def getIRowBig(self,f):
        data = f.read(8); self.n+=8
        (idummy,irow) = unpack('2i',data)
        return irow

    def readMatrix(self,f,floatType,matrixNames=None):
        """reads a matrix"""
        data = f.read(4); self.n+=4
        (recordLength,) = unpack(self.endian+'i',data)
        
        data = f.read(recordLength); self.n+=recordLength
        if recordLength==24:
            (ncols,nrows,form,Type,name) = unpack(self.endian+'4i8s',data)
            #print "N=%s recordLength=%s ncols=%s nrows=%s form=%s Type=%s name=%s" %(self.n,recordLength,ncols,nrows,form,Type,name)
        else:
            raise NotImplementedError('recordLength=%s\n%s' %(recordLength,self.printSection(60)))

        if 0:
            if Type==1:
                print "Type = Real, Single Precision"
            elif Type==2:
                print "Type = Real, Double Precision"
            elif Type==3:
                print "Type = Complex, Single Precision"
            elif Type==4:
                print "Type = Complex, Double Precision"

        if nrows<0: # if less than 0, big
            isBigMat = True
            nrows = abs(nrows)
        elif nrows>0:
            isBigMat = False
        else:
            raise RuntimeError('unknown BIGMAT.  nRows=%s' %(nrows))
        
        if Type==1:
            NWV = 1 # number words per value
            d = 'f'
            dType = 'float32'
        elif Type==2:
            NWV = 2
            d = 'd'
            dType = 'float64'
        elif Type==3:
            NWV = 2
            d = 'ff'
            dType = 'complex64'
        elif Type==4:
            NWV = 4
            d = 'dd'
            dType = 'complex128'
        else:
            raise RuntimeError("Type=%s" %(Type))

        # reset the type if 'default' not selected
        if floatType=='single':
            if Type==[1,2]:
                dType='float32'
            else:
                dType='complex64'
        elif floatType=='double':
            if Type==[1,2]:
                dType='float64'
            else:
                dType='complex128'

        # jump forward to get if isSparse, then jump back
        nSave = self.n
        (_a,_icol,_irow,_nWords) = self.readStartMarker(f)
        f.seek(nSave); self.n=nSave

        isSparse = False
        if _irow==0:
            isSparse = True
            rows = []
            cols = []
            entries = []
        else:
            A = zeros((nrows,ncols),dType)


        icol=-1 # dummy value so the loop starts
        if Type in [1,2]: # real
            while icol<ncols+1: # if isDense
                (icol,irow,nWords) = self.getMarkers(f,isSparse,isBigMat)

                if nWords==0 and isBigMat:
                    self.n-=4; f.seek(self.n)
                    break

                recordLength = 4*nWords
                data = f.read(recordLength); self.n+=recordLength

                if icol==ncols+1:
                    continue

                nValues = nWords//NWV
                if nValues==0:
                    assert icol==ncols+1
                    break

                strValues = nValues*d
                valueList = unpack(strValues,data)
                
                irow-=1
                icol-=1
                if isSparse:
                    cols += [icol]*nValues
                    rows += [i+irow for i in range(nValues)]
                    for value in valueList:
                        entries.append(value)
                        irow+=1
                else:
                    for value in valueList:
                        A[irow,icol] = value
                        irow+=1

        elif Type in [3,4]: # complex
            while icol<ncols+1: # if isDense
                (icol,irow,nWords) = self.getMarkers(f,isSparse,isBigMat)

                if nWords==0 and isBigMat:
                    self.n-=4; f.seek(self.n)
                    break

                recordLength = 4*nWords
                data = f.read(recordLength); self.n+=recordLength

                nValues = nWords//NWV
                if nValues==0:
                    assert icol==ncols+1
                    break

                strValues = nValues*d
                valueList = unpack(strValues,data)

                if icol==ncols+1:
                    continue
                
                irow-=1
                icol-=1

                if isSparse:
                    cols += [icol]*nValues
                    rows += [i+irow for i in range(nValues)]
                    for i,value in enumerate(valueList):
                        if i%2==0:
                            realValue = value
                        else:
                            #A[irow,icol] = complex(realValue,value)
                            entries.append(complex(realValue,value))
                            irow+=1
                else:
                    for i,value in enumerate(valueList):
                        if i%2==0:
                            realValue = value
                        else:
                            A[irow,icol] = complex(realValue,value)
                            irow+=1
        else:
            raise RuntimeError("Type=%s" %(Type))
        #print printMatrix(A)

        if d in ['d','dd']:
            f.read(8); self.n+=8
        elif d in ['f','ff']:
            f.read(4); self.n+=4
        else:
            raise NotImplementedError(d)
        
        if isSparse:
            #print "len(rows)=%s len(cols)=%s len(entries)=%s" %(len(rows),len(cols),len(entries))
            A = coo_matrix( (entries,(rows,cols)),shape=(nrows,ncols),dtype=dType) # Initialize a real matrix
            #print "type = %s %s" %(type(A),type(A.todense()))
            #A = A.todense()
        return (name,A)

    def getMarkers(self,f,isSparse,isBigMat):
        if isSparse:
            if isBigMat:
                (a,icol,irow,nWords) = self.readStartMarker(f)
                (irow) = self.getIRowBig(f)
                if nWords>1:
                    nWords -= 2
                else:
                    nWords = 0
            else:
                (a,icol,irow,nWords) = self.readStartMarker(f)
                if irow!=0:
                    assert nWords==1,'nWords=%s' %(nWords)

                (irow) = self.getIRowSmall(f)
                nWords -= 1
            ###
        else:
            (a,icol,irow,nWords) = self.readStartMarker(f)
        return (icol,irow,nWords)

if __name__=='__main__':
    filename = 'test/mat_b_dn.op4' # works
    filename = 'test/mat_b_s1.op4' # works
    filename = 'test/mat_b_s2.op4'
    #matrixNames = 'RND1RD'
    matrixNames = None
    op4 = OP4()
    matrices = op4.readOP4(filename,matrixNames=matrixNames)
    for name,matrix in sorted(matrices.items()):
        print "name = %s" %(name)
        if isinstance(matrix,coo_matrix):
            matrix = matrix.todense()
        
        print printMatrix(matrix)
    print "-----------------------------"
    print "done"
    print "-----------------------------"
