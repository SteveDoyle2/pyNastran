## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
#!/usr/bin/env python
#
# Test the libop4 <-> matlab interface by writing random matrices to
# a file, reading the matrices back into new variables, then finding
# the largest difference between the written and read matrices.
#
# The matrices are evenly divided among these categories:
#   1.  square / rectangular
#   2.  sparse / dense           (if sparse, density varies randomly on 0..1)
#   3.  real   / complex
#   4.  text   / binary          (if text, DIGITS varies randomly on 4..26)
#
# Untested with this program:
#   - Multiple matrices per .op4 file.
#   - nSkip
#   - binary single precision
#
# Network delays sometimes cause the file 'junk_<machine name>.op4' to be
# unreadable during the loadop4() step.  Run from a directory mounted on a
# local disk to avoid this.
#
# Albert Danial March 14 2001
#
import numpy        as np
import scipy.sparse as S
import cop4 

nTests      = 200
approx_size =  20
nErrors     =   0 
max_E       =   0 
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
        M = 20*(0.5 - np.random.rand(nRows, nCols))
        cat_2   = 'dense'
        density = 1

    # Category III:  real or complex
    if np.random.random() < 0.5:  #  real
#   if False:
        cat_3 = 'real'
    else:                       #  complex
        cat_3 = 'complex'
        M = M * (1 + .4j)

    # Category IV:  binary or text
    if np.random.random() < 0.5:  #  binary
#   if True:
        cat_4 = 'binary'
        digits = 0
    else:                       #  text
        cat_4 = 'text'
        digits = 4 + np.random.randint(22)

    print('Test %3d/%d  %-11s %-6s (rho=%6.4f) %-7s %-6s  %4d x %4d' % (
            iSet+1, nTests, cat_1, cat_2, density, cat_3, cat_4, nRows, nCols))

    cop4.Save(filename, M, digits=digits)
    fh = cop4.OP4(filename, 'r')
    L = fh.Load(nmat=1)

    error = np.abs(np.max(M - L))
    if cat_2 == 'dense':
        pass
    else:
        error = np.abs(np.max(M - L))
        if type(error) == np.float64:
            # happens when input is 1x1, or error is exactly zero
            pass
        else:
            if error.getnnz() == 0: 
                error = 0.0
            else:
                error = np.max(error.data)
    print('error=%8.5e' % error)
    max_E = max(error, max_E)
#   clear L, M;   % otherwise will leak memory like a sieve

    if digits:
        allowable_error = max(10**(-digits), 1e-5) # min digits is 5
    else:
        allowable_error = 1e-8;
    if error > allowable_error:
        nErrors = nErrors + 1;
        print('error exceeds allowable of %8.5e' % allowable_error);
        break

print('%d tests, %d errors (worst error %e)' % (nTests, nErrors, max_E))
