import os
import sys
from struct import unpack
from numpy import zeros
from scipy.sparse import coo_matrix
from pyNastran.op2.fortranFile import FortranFile
from pyNastran.general.generalMath import printMatrix

class OP4(FortranFile):
    def __init__(self):
        FortranFile.__init__(self)
        self.makeOp2Debug = False
        self.n = 0
        self.endian = ''

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
        self.op2 = open(op4Name,'rb')
        #print self.printSection(80)
        
        print '------------------------'
        #print self.printSection(20)
        #self.op2.read(4); self.n+=4 # record length

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

            (name,matrix) = self.readMatrix(self.op2,matrixNames)
            if name is not None:
                if matrixNames is None or name in matrixNames: # save the matrix
                    matrices[name] = matrix
            i+=1
        return matrices

    def readMatrix(self,f,matrixNames=None):
        """reads a matrix"""
        print "***********starting**************"
        #print self.printSection(80)
        data = f.read(28); self.n+=28
        
        (recordLength,ncols,nrows,form,Type,name) = unpack(self.endian+'5i8s',data)

        print "*N=%s recordLength=%s nrows=%s ncols=%s form=%s Type=%s name=|%s|" %(
               self.n,recordLength,ncols,nrows,form,Type,name)
        assert Type in [1,2,3,4],'Invalid Type.  Type=%s.  1,2=Real; 3,4=Complex; 1,3=Single Precision; 2,4=Double Precision' %(Type)
        assert 0<nrows<40
        assert 0<ncols<40

        if nrows<0: # if less than 0, big
            isBigMat = True
        elif nrows>0:
            isBigMat = False
        else:
            raise RuntimeError('unknown BIGMAT.  nRows=%s' %(nrows))
        #print "isBigMat=%s" %(isBigMat)


        name = name.strip()
        if not name.isalnum():
            sys.exit('matrix name is not an ASCII string...')

        data = f.read(20); self.n+=20
        (recordLength,b,icol,irow,nWords) = unpack(self.endian+'5i',data)
        print "N=%s recordLength=%s b=%s irow=%s icol=%s nWords=%s" %(self.n,recordLength,b,irow,icol,nWords)

        isSparse = False
        if irow==0:
            isSparse = True
        print "isBigMat=%s isSparse=%s" %(isBigMat,isSparse)

        if Type in [1,2]: # real
            (A) = self.readReal(f,name,nrows,ncols,irow,icol,nWords,Type,isSparse,isBigMat)
        elif Type in [3,4]: # complex
            (A) = self.readComplex(f,nrows,ncols,irow,icol,nWords,Type,isSparse,isBigMat)
        else:
            raise RuntimeError('invalid matrix type.  Type=%s' %(Type))

        print "N=%s name=%s is finished\n" %(self.n,name)
        #print dir(A)
        #print A
        print printMatrix(A)

        if not(matrixNames is None or name in matrixNames): # kill the matrix
            A = None
        return name,A


    def readReal(self,f,name,nrows,ncols,irow,icol,nWords,Type,isSparse,isBigMat):
        """
        @todo possibly split this into readDenseReal and readSparseReal
        to get rid of all the extra isSparse checks.  This would cleanup the
        runLoop condition as well.
        """
        print "*********readReal************"
        if isSparse:
            rows=[]; cols=[]; entries=[]
        else:
            A = zeros((nrows,ncols),'f') # Initialize a real matrix

        nLoops = 0
        wasBroken=False
        while 1:
            #if nLoops>0 and not wasBroken:
                #line = f.readline().rstrip()
                #asf
            isRewound = False

            if nLoops>0:
                #print self.printSection(120)
                n = self.n
                icolOld=icol
                irowOld=irow
                print "N = ",self.n
                data = f.read(20); self.n+=20
                (recordLength,b,icol,irow,nWords) = unpack(self.endian+'5i',data)
                print "N=%s iLoop=%s recordLength=%s b=%s irow=%s icol=%s nWords=%s" %(self.n,nLoops,recordLength,b,irow,icol,nWords)

                if icol==0: # rewind
                    self.n = n; f.seek(n)
                    print "N = ",self.n
                    isRewound = True
                    icol=icolOld
                    irow=irowOld
                    
                    print self.printSection(60)
                    #data = f.read(4); self.n+=4
                    data = f.read(16); self.n+=16
                    (b,icol,irow,nWords) = unpack(self.endian+'4i',data)
                    print "****N=%s iLoop=%s b=%s irow=%s icol=%s nWords=%s" %(self.n,nLoops,b,irow,icol,nWords)
                    sys.stdout.flush()
                    #break
                #print "b=%s" %(b)
            ###

            #print "b=%s" %(b)

            irow,isEndSparse = self.getIRow(f,irow,nWords,Type,isSparse,isBigMat)
            if isEndSparse and nWords//2!=0:
                print self.printSection(40)
                break

            if irow>nrows+1:
                sys.exit('what happened irow is too large...nrow=%s irow=%s icol=%s' %(nrow,irow,icol))

            nFloats = nWords//2
            print "nWords=%s nFloats=%s" %(nWords,nFloats)
            #if nWords>2:
                #print self.printSection(40)

            assert nFloats<1000
            assert nWords//4<nrows*ncols

            if nFloats > 0 and not isEndSparse:
                if   Type==1: # real single precision
                    #if isSparse and name != 'STRINGS': ## @todo this is obviously wrong...
                    #    nFloats = nWords//2*2
                    #else:
                    nFloats = nWords
                    print "***nFloats = %s" %(nFloats)
                    data = f.read(nFloats*4); self.n+=nFloats*4
                    packStr = self.endian+'f'*nFloats
                elif Type==2: # real double precision
                    nDoubles = nWords//2
                    print "***nDoubles = %s" %(nDoubles)
                    data = f.read(nDoubles*8); self.n+=nDoubles*8
                    packStr = self.endian+'d'*nDoubles
                else:
                    raise NotImplementedError('Type=%s 3=complex single; 4=complex double' %(Type))
                print "N = %s" %(self.n)
                #print "packStr = ",packStr,len(data)
                valueList = unpack(packStr,data)
                for value in valueList:
                    #print "irow=%s icol=%s nWords=%s nFloats=%s" %(irow,icol,nWords,nFloats)
                    if irow<=nrows:
                        print " A[%s,%s] = %f" %(irow,icol,value)
                        try:
                            if isSparse:
                                rows.append(irow-1)
                                cols.append(icol-1)
                                entries.append(value)
                            else:
                                A[irow-1,icol-1] = value
                            #print "*A[%s,%s] = %s" %(irow,icol,A[irow-1,icol-1])
                        except IndexError:
                            print A.shape
                            raise
                    else:
                        print "*A[%s,%s] = %f" %(irow,icol,value)
                    irow += 1
                #print self.printSection(80)
            else:
                if   Type==1: # real single precision
                    data = f.read(4); self.n+=4
                    dummy, = unpack(self.endian+'f',data)
                    #f.read(4); self.n+=4
                elif Type==2: # real double precision
                    data = f.read(8); self.n+=8
                    dummy, = unpack(self.endian+'d',data)
                    #f.read(8); self.n+=8
                else:
                    raise NotImplementedError('Type=%s 1=real single; 2=real double' %(Type))
                print "***end of matrix; dummy = %s" %(dummy)
                break
            ###
            nLoops += 1
        ###
        print '----------------------------------'
        if not isSparse:
            f.read(4); self.n+=4
        #print self.printSection(100)
        if isSparse:
            A = coo_matrix( (entries,(rows,cols)),shape=(nrows,ncols),dtype='f') # Initialize a real matrix
            #print "type = %s %s" %(type(A),type(A.todense()))
            #A = A.todense()
        return A

    def readComplex(self,f,nrows,ncols,irow,icol,nWords,Type,isSparse,isBigMat):
        print "*********readComplex************"
        if isSparse:
            rows=[]; cols=[]; entries=[]
        else:
            A = zeros((nrows,ncols),'complex') # Initialize a complex matrix

        nLoops = 0
        wasBroken=False
        while 1:
            print "restarting loop..."
            #if nLoops>0 and not wasBroken:
                #line = f.readline().rstrip()
                #asf
            isRewound = False

            if nLoops>0:
                #print self.printSection(120)
                n = self.n
                icolOld=icol
                irowOld=irow
                
                data = f.read(20); self.n+=20
                (recordLength,b,icol,irow,nWords) = unpack(self.endian+'5i',data)
                print "N=%s iLoop=%s recordLength=%s b=%s irow=%s icol=%s nWords=%s" %(self.n,nLoops,recordLength,b,irow,icol,nWords)
                print "~~~~~~~~"

                if icol==0: # rewind
                    self.n = n
                    isRewound = True
                    asdf
                    break
                #print "b=%s" %(b)
            ###


            #runLoop = True
            #while (len(sline)==1 or len(sline)==2) and 'E' not in line or runLoop: # next sparse entry
            (irow,isEndSparse) = self.getIRow(f,irow,nWords,Type,isSparse,isBigMat)

            if isEndSparse and nWords//2!=0:
                print self.printSection(40)
                break

            if irow>nrows+1:
                sys.exit('what happened irow is too large...nrow=%s irow=%s icol=%s' %(nrow,irow,icol))

            nFloats = nWords//2
            print "nWords=%s nFloats=%s" %(nWords,nFloats)
            #if nWords>2:
                #print self.printSection(40)

            assert nFloats<1000
            assert nWords//4<nrows*ncols

            if nFloats > 0 and not isEndSparse:
                if   Type==3: # complex single precision
                    #if nWords%2==1:  # @todo why is this needed???
                        #f.read(4); self.n+=4
                    nFloats = nWords//2*2
                    data = f.read(nFloats*4); self.n+=nFloats*4
                    packStr = self.endian+'f'*nFloats
                elif Type==4: # complex double precision
                    nDoubles = nWords//2
                    data = f.read(nDoubles*8); self.n+=nDoubles*8
                    packStr = self.endian+'d'*nDoubles
                else:
                    raise NotImplementedError('Type=%s 3=complex single; 4=complex double' %(Type))
                #print "packStr = ",packStr,len(data)
                valueList = unpack(packStr,data)
                assert len(valueList)%2==0
                for i,value in enumerate(valueList):
                    #print "irow=%s icol=%s nWords=%s nFloats=%s" %(irow,icol,nWords,nFloats)
                    
                    if i%2==0:
                        #print "A[%s,%s] = %f" %(irow,icol,value)
                        realValue = value
                    else:
                        print "A[%s,%s] = %s %fj" %(irow,icol,realValue,value)
                        if isSparse:
                            rows.append(irow-1)
                            cols.append(icol-1)
                            entries.append(complex(realValue,value))
                        else:
                            A[irow-1,icol-1] = complex(realValue,value)
                        irow += 1

                #print self.printSection(80)
            else:
                if   Type==3: # complex single precision
                    data = f.read(4); self.n+=4
                    dummy, = unpack(self.endian+'f',data)
                    #f.read(4); self.n+=4
                elif Type==4: # complex double precision
                    data = f.read(8); self.n+=8
                    dummy, = unpack(self.endian+'d',data)
                    #f.read(8); self.n+=8
                else:
                    raise NotImplementedError('Type=%s 3=complex single; 4=complex double' %(Type))
                print "***end of matrix; dummy = %s" %(dummy)
                print "broken...."
                break
            ###
            nLoops += 1
        ###
        print '----------------------------------'
        if not isSparse:
            f.read(4); self.n+=4
        #print self.printSection(100)
        if isSparse:
            A = coo_matrix( (entries,(rows,cols)),shape=(nrows,ncols),dtype='complex') # Initialize a real matrix
            #print "type = %s %s" %(type(A),type(A.todense()))
            #A = A.todense()
        return A

    def getIRow(self,f,irow,nWords,Type,isSparse,isBigMat):
        isEndSparse = False
        if isSparse:
            #nWords = (nWords-1)//2  ## @todo this cant be right...
            
            if isBigMat:
                raise NotImplementedError('isBigMat=True')
            #if isBigMat:
                #if len(sline)==2:
                #    pass
                #else:
                #    sline = f.readline().strip().split()
                #assert len(sline)==2,'sline=%s len(sline)=%s' %(sline,len(sline))
                #(idummy,irow) = sline
                #irow = int(irow)
            else:
                print "***************************"
                #print self.printSection(100)
                data = f.read(4); self.n+=4
                IS, = unpack('i',data)
                #IS = int(f.readline().strip())
                L = IS//65536 - 1
                irow = IS - 65536*(L + 1)
                print "IS=%s L=%s irow=%s" %(IS,L,irow)
                if L==-1:
                    isEndSparse = True
                ###
            ###
        ###
        
        return irow,isEndSparse

if __name__=='__main__':
    filename = 'test/mat_b_dn.op4' # works
    filename = 'test/mat_b_s1.op4'
    #filename = 'test/mat_b_s2.op4'
    #matrixNames = 'EYE5CD'
    matrixNames = None
    op4 = OP4()
    matrices = op4.readOP4(filename,matrixNames=matrixNames)
    for name,matrix in sorted(matrices.items()):
        print "name = %s" %(name)
        
        print printMatrix(matrix)
    print "-----------------------------"
    print "done"
    print "-----------------------------"
