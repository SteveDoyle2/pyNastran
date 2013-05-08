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

"""
Create a collection of dense binary matrices in the current directory
then see how long it takes to read the values back in using the
pure Python and the Cython OP4 modules.  Plot the timing results.
al.danial@gmail.com March 2013
"""

import time
import numpy as np
import os.path
import sys
from matplotlib import pyplot

import cop4          # Cython
from op4 import OP4  # pure Python

if __name__ == '__main__':
    dimension = range(1000, 15000, 1000) # 14 matrices sized 1000 x 1000 to
                                         # 14000 x 14000.  The largest will
                                         # be 1.5 GB.
    dimension = range(1000, 5000, 1000)
    #print dimension

    File      = []
    cop4_save = []
    cop4_load = []
    pop4_load = []
    for i,n in enumerate(dimension): File.append( 'A_%dx%d.op4' % (n, n) )

    for i,n in enumerate(dimension):
        if not os.path.exists(File[i]):
            print('%s does not exist, will create.' % (File[i]))
            A = np.random.rand(n,n).astype(np.float32)
            start_elapsed = time.time()
            cop4.Save(File[i], A=A)
            end_elapsed = time.time()
            print('wrote %-20s in %8.3f s' % (File[i], end_elapsed-start_elapsed))

    # Method 1: cop4
    for i in range(len(dimension)):
        start_elapsed = time.time()
        B = cop4.Load(File[i])
        end_elapsed = time.time()
        cop4_load.append(end_elapsed - start_elapsed)
        print('cop4 read %-20s in %8.3f s' % (File[i], cop4_load[i]))

    # Method 2: pure Python op4
    for i in range(len(dimension)):
        start_elapsed = time.time()
        op4 = OP4()
        B = op4.readOP4(File[i])
        end_elapsed = time.time()
        pop4_load.append(end_elapsed - start_elapsed)
        print('pop4 read %-20s in %8.3f s' % (File[i], pop4_load[i]))

    print('Dimension     C op4       Python op4')
    print('---------   ---------     ----------')
    for i,n in enumerate(dimension):
        print('%6d     %8.3f       %8.3f' % (n, cop4_load[i], pop4_load[i]))

    pyplot.semilogy(dimension, cop4_load, color='green')
    pyplot.semilogy(dimension, pop4_load, color='blue' )
    pyplot.grid(True)
    pyplot.xlabel('Matrix Dimension')
    pyplot.ylabel('Time [seconds]')
    pyplot.title('Square Matrix Load Performance')
    pyplot.show()