import numpy
from numpy import array,cross,allclose,zeros,matrix
from scipy.linalg import solve_banded

from numpy import 
from numpy.linalg import norm, solve

def reduceMatrix(matA,nids):
    """
    takes a list of ids and removes those rows and cols
    """
    nRows = len(nids)
    matB = matrix(zeros((nRows,nRows),'d') )

    for i,irow in enumerate(nids):
        for j,jcol in enumerate(nids):  
            matB[i,j] = matA[irow,jcol]
    return matB

def isListRanged(a,List,b):
    """
    Returns true if a<= x <= b
    or a-x < 0 < b-x
    """
    for x in List:
        if not isFloatRanged(a,x,b):
            return False
        ###
    ###
    return True

def isFloatRanged(a,x,b):
    """
    Returns true if a<= x <= b
    or a-x < 0 < b-x
    """
    if not a<x:
        if not allclose(x,a):
           return False 
        ###
    ###
    if not x<b:
        if not allclose(x,b):
            return False
        ###
    ###
    return True

def printMatrix(A):
    msg = ''
    for row in A:
        msg += ListPrint(row)+'\n'
    ###
    return msg
###

def ListPrint(listA):
    if len(listA)==0:
        return '[]'
    ###

    msg = '['
    for a in listA:
        if isinstance(a,str):
            msg += ' %s,' %(a)
        elif isinstance(a,float):
            msg += ' %-4.2f,' %(a)
        elif isinstance(a,int):
            msg += ' %g,' %(a)
        else:
            try:
                msg += ' %g,' %(a)
            except TypeError:
                print "a = |%s|" %(a)
                raise
            ###
        ###
    ###
    msg = msg[:-1]
    msg += ' ]'
    return msg

def augmentedIdentity(A):
    """
    Creates an Identity Matrix augmented with zeros.
    The location of the extra zeros depends on A.
    """
    (nx,ny) = A.shape
    #(ny,nx) = A.shape
    I = zeros([nx,ny],'d')
    
    for i in range(nx):
        if i==nx or i==ny:
            break
        I[i][i] = 1.
    return I

def solveTridag(A, D):
     # Find the diagonals
     ud = numpy.insert(numpy.diag(A,1), 0, 0) # upper diagonal
     d  = numpy.diag(A) # main diagonal
     ld = numpy.insert(numpy.diag(A,-1), len(d)-1, 0) # lower diagonal
   
     # simplified matrix
     ab = numpy.matrix([
         ud,
         d,
         ld,
     ])
     #print "ab = ",ab
     return solve_banded((1, 1), ab, D,overwrite_ab=True,overwrite_b=True)

def Area(a,b):
    return 0.5*norm(cross(a,b))

def AreaNormal(nodes):
    """
    Returns area,unitNormal
    n = Normal = a x b
    Area   = 1/2 * |a x b|
    V = <v1,v2,v3>
    |V| = sqrt(v1^0.5+v2^0.5+v3^0.5) = norm(V)
    
    Area = 0.5 * |n|
    unitNormal = n/|n|
    """
    (n0,n1,n2) = nodes
    a = n0-n1
    b = n0-n2
    vector = cross(a,b)
    length = norm(vector)
    normal = vector/length
    area = 0.5*length
    if allclose(norm(normal),1.)==False:
        print "a = ",a
        print "b = ",b
        print "normal = ",normal
        print "length = ",length
        sys.exit('check...')
    return (area,normal)

def Triangle_AreaCentroidNormal(nodes):
    """Returns area,centroid,unitNormal"""
    (area,normal) = AreaNormal(nodes)
    centroid = Centroid(*nodes)
    return (area,centroid,normal)

def Normal(a,b):
    """finds the unit normal vector of 2 vectors"""
    vector = cross(a,b)
    length = norm(vector)
    normal = vector/length
    assert allclose(norm(normal),1.)
    return normal

def Centroid(A,B,C):
    """returns the centroid of a triangle"""
    #print "type(A,B,C) = ",type(A),type(C),type(B)
    centroid = (A+B+C)/3.
    return centroid
