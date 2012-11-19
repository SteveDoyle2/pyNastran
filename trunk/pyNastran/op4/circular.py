#!/usr/bin/env python
#
# Test loading and saving op4 files with a circular test:
#   1. create a matrix with having four different attributes
#       a.  square or rectangular
#       b.  sparse or dense      (if sparse, density varies randomly on 0..1)
#       c.  real   or complex
#       d.  text   or binary     (if text, DIGITS varies randomly on 4..26)
#   2. save the matrix to an op4 file, 'junk.op4' in the current directory
#   3. read the matrix back from the op4 file
#   4. compute the error norm of the original matrix made at
#      step 1 with the one read back from file in step 3
#
import numpy        as np
import scipy.sparse as S
import cop4 

nTests      = 200
approx_size =  10 
nErrors     =  0 
max_E       =  0 
filename    = 'junk.op4'

for iSet in range(nTests):
    # Category I:  square or rectangular
    if np.random.random() < 0.5:  #   square
        cat_1 = 'square'
        nRows = np.random.randint(approx_size) + 1
        nCols = nRows
    else:                       #  rectangular
        cat_1 = 'rectangular'
        nRows = np.random.randint(approx_size) + 1
        nCols = np.random.randint(approx_size) + 1

    # Category II:  sparse or dense 
    #    matrix values will be random between -10..10
    if np.random.random() < 0.5:  #  sparse
        cat_2   = 'sparse'
        density = np.random.random()
        M = S.rand(nRows, nCols, density)
        M.data = 20*(0.5 - M.data)
    else:                       #  dense
        M = 2*(0.5 - np.random.rand(nRows, nCols))
        cat_2   = 'dense'
        density = 1

    # Category III:  real or complex
#   if np.random.random() < 0.5:  #  real
    if True:
        cat_3 = 'real'
    else:                       #  complex
        cat_3 = 'complex'
        M = M * (1 + .4j)

    # Category IV:  binary or text
    if np.random.random() < 0.5:  #  binary
#   if False:
        cat_4 = 'binary   '
        digits = 0
    else:                       #  text
        digits = 4 + np.random.randint(22)
        cat_4 = 'text (%2d)' % (digits)

    cop4.Save(filename, M, digits=digits)
    fh = cop4.OP4(filename, 'r')
    L  = fh.Load(nmat=1)

    if cat_2 == 'sparse':
        if type(L) is np.ndarray:
            Diff = M.todense() - L
        else:
            Diff = M.todense() - L.todense()
    else:
        Diff = M - L

    error = np.linalg.norm(Diff)
    I = np.argmin(Diff) / nCols
    J = np.argmin(Diff) % nCols

    print('%3d/%d %-11s %-6s (d=%6.4f) %-6s %-6s %3d x %3d |%8.5e|' % (
            iSet+1, nTests, cat_1, cat_2, density, cat_3, cat_4, 
            nRows, nCols, error))

    max_E = max(error, max_E)

    if digits:
        allowable_error = max(10**(2-digits), 1e-8)
    else:
        allowable_error = 1e-8
    if error > allowable_error:
        nErrors = nErrors + 1
        if cat_2 == 'sparse':
            print('error exceeds %8.5e ([%d,%d]: %e v. %e)' % (
                allowable_error, I, J, M.todense()[I,J], L.todense()[I,J]))
        else:
            print('error exceeds %8.5e ([%d,%d]: %e v. %e)' % (
                allowable_error, I, J, M[I,J], L[I,J]))

print('%d tests, %d errors (worst error %e)' % (nTests, nErrors, max_E))
