#!/usr/bin/env python

"""
Create a collection of dense binary matrices in the current directory
then see how long it takes to read the values back in using the
pure Python and the Cython OP4 modules.  Plot the timing results.
al.danial@gmail.com March 2013
"""
from __future__ import print_function
import time
import numpy as np
import os
import sys
from matplotlib import pyplot

import cop4          # Cython
from op4 import OP4  # pure Python

def dense_square_test(dimension, is_binary=True):
    cyfiles = []
    pyfiles = []
    cop4_save = []
    cop4_load = []
    pop4_load = []
    for i,n in enumerate(dimension):
        cyfiles.append( 'cop4_%dx%d.op4' % (n, n) )

    if is_cyfiles:
        # create dense matrices using cop4
        for i,n in enumerate(dimension):
            if not os.path.exists(cyfiles[i]):
                print('%s does not exist, will create.' % cyfiles[i])
                A = np.random.rand(n,n).astype(np.float32)
                start_elapsed = time.time()
                mat_name = 'matrix%-3i' % i
                kwargs = {
                    mat_name : A,
                    'digits' : 5,
                }
                #print('%r' % mat_name)
                if is_binary:
                    del kwargs['digits']
                cop4.Save(cyfiles[i], **kwargs)
                #cop4.Save(File[i], A=A, digits=5)  # ascii
                #cop4.Save(File[i], A=A)            # binary
                end_elapsed = time.time()
                print('wrote %-20s in %8.3f s' % (cyfiles[i], end_elapsed-start_elapsed))
        filenames = cyfiles
    else:
        # create dense matrices using pyop4
        for i,n in enumerate(dimension):
            pyfiles.append( 'pyop4_%dx%d.op4' % (n, n) )

        pyop4 = OP4()
        for i,n in enumerate(dimension):
            #if not os.path.exists(File[i]):
                print('%s does not exist, will create.' % pyfiles[i])
                A = np.random.rand(n,n).astype(np.float32)
                #start_elapsed = time.time()
                mat_name = 'mat%i' % i
                #end_elapsed = time.time()
                #print('wrote %-20s in %8.3f s' % (pyfiles[i], end_elapsed-start_elapsed))

                #def write_op4(op4_filename, names, matrices, forms, precisions, is_binary=True):
                names = [mat_name]
                matrices = [A]
                forms = [1]
                #precisions = [']
                pyop4.write_op4(pyfiles[i], names, matrices, forms, precisions='default', is_binary=is_binary)

        filenames = pyfiles
    print("filenames =", filenames)
    if is_cyfiles:
    #if 1:
        # Method 1: cop4
        for i in xrange(len(dimension)):
            start_elapsed = time.time()
            B = cop4.Load(filenames[i], '*')  # dictionary-based
            print(B.keys())
            key0 = B.keys()[0]
            #print(B[key0])
            end_elapsed = time.time()
            del B
            cop4_load.append(end_elapsed - start_elapsed)
            print('cop4 read %-20s in %8.3f s' % (filenames[i], cop4_load[i]))

    # Method 2: pure Python op4
    for i in xrange(len(dimension)):
        start_elapsed = time.time()
        op4 = OP4()
        B = op4.read_op4(filenames[i])
        print("keys = %s" % B.keys())
        key0 = B.keys()[0]
        print("form = %s" % B[key0][0])
        #print(B[key0][1], '\n')

        end_elapsed = time.time()
        del B, op4
        pop4_load.append(end_elapsed - start_elapsed)
        print('pop4 read %-20s in %8.3f s' % (filenames[i], pop4_load[i]))

    print('Dimension     C op4       Python op4')
    print('---------   ---------     ----------')
    for i,n in enumerate(dimension):
        print('%6d     %8.3f       %8.3f' % (n, cop4_load[i], pop4_load[i]))

    if 1:
        pyplot.semilogy(dimension, cop4_load, color='green')
        pyplot.semilogy(dimension, pop4_load, color='blue' )
        pyplot.legend(['cop4', 'op4'], 'upper left')
        pyplot.grid(True)
        pyplot.xlabel('Matrix Dimension')
        pyplot.ylabel('Time [seconds]')
        pyplot.title('Dense Square Matrix Load Performance')

if __name__ == '__main__':
    is_cyfiles = True
    dimension = range(1000, 15000, 1000) # 14 matrices sized 1000 x 1000 to
                                         # 14000 x 14000.  The largest will
                                         # be 1.5 GB.
    dimension = range(1000, 5000, 1000)
    #dimension = range(5, 10, 10)
    #print dimension
    dense_square_test(dimension)

    pyplot.show()