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
from matplotlib import pyplot as plt

import cop4          # Cython
from op4 import OP4  # pure Python

def dense_square_test(dimension, is_binary=True):
    filenames = []
    cop4_save = []
    cop4_load = []
    pop4_load = []

    # generate the file list
    for i,n in enumerate(dimension):
        if is_cyfiles:
            if is_binary:
                filenames.append( 'cop4_%dx%d_binary.op4' % (n, n) )
            else:
                filenames.append( 'cop4_%dx%d_ascii.op4' % (n, n) )
        else:
            if is_binary:
                filenames.append( 'pyop4_%dx%d_binary.op4' % (n, n) )
            else:
                filenames.append( 'pyop4_%dx%d_ascii.op4' % (n, n) )
    print("is_cyfiles=%s is_binary=%s" % (is_cyfiles, is_binary))

    # create dense matrices using cop4
    pyop4 = OP4()
    for i,n in enumerate(dimension):
        if not os.path.exists(filenames[i]):
            print('%s does not exist, will create.' % filenames[i])

            A = np.random.randn(n,n) #.astype(np.float32)
            #print('*****\n', A)
            mat_name = 'matrix%-3i' % i
            if is_cyfiles:
                start_elapsed = time.time()
                kwargs = {
                    mat_name : A,
                    'digits' : 5,
                }
                if is_binary:
                    del kwargs['digits']
                else:
                    print('digits =', kwargs['digits'])
                cop4.Save(filenames[i], **kwargs)
                #cop4.Save(File[i], A=A, digits=5)  # ascii
                #cop4.Save(File[i], A=A)            # binary
                end_elapsed = time.time()
            else:
                start_elapsed = time.time()
                matrices = {'A' : (1, A) } # 1 is the form -> square
                names = None  # name order -> sorted
                pyop4.write_op4(filenames[i], matrices, names, precision='default', is_binary=is_binary)
                end_elapsed = time.time()
            print('wrote %-20s in %8.3f s' % (filenames[i], end_elapsed-start_elapsed))

    print("filenames =", filenames)
    # Method 1: cop4
    for i in xrange(len(dimension)):
        start_elapsed = time.time()
        B = cop4.Load(filenames[i], '*')  # dictionary-based
        end_elapsed = time.time()
        #print(B.keys())
        key0 = B.keys()[0]
        #print(B[key0])
        del B
        cop4_load.append(end_elapsed - start_elapsed)
        print('cop4 read %-20s in %8.3f s' % (filenames[i], cop4_load[i]))

    # Method 2: pure Python op4
    for i in xrange(len(dimension)):
        start_elapsed = time.time()
        B = pyop4.read_op4(filenames[i])
        end_elapsed = time.time()
        #print("keys = %s" % B.keys())
        key0 = B.keys()[0]
        #print("form = %s" % B[key0][0])
        #print(B[key0][1], '\n')

        del B
        pop4_load.append(end_elapsed - start_elapsed)
        print('pop4 read %-20s in %8.3f s' % (filenames[i], pop4_load[i]))

    print('Dimension     C op4       Python op4')
    print('---------   ---------     ----------')
    for i,n in enumerate(dimension):
        print('%6d     %8.3f       %8.3f' % (n, cop4_load[i], pop4_load[i]))

    if make_plot:
        plt.semilogy(dimension, cop4_load, color='green')
        plt.semilogy(dimension, pop4_load, color='blue' )
        plt.legend(['cop4', 'op4'], 'upper left')
        plt.grid(True)
        plt.xlabel('Matrix Dimension')
        plt.ylabel('Time [seconds]')
        plt.title('Dense Square Matrix Load Performance')

if __name__ == '__main__':
    #is_cyfiles = False
    is_cyfiles = True

    make_plot = False
    make_plot = True
    dimension = range(1000, 15000, 1000) # 14 matrices sized 1000 x 1000 to
                                         # 14000 x 14000.  The largest will
                                         # be 1.5 GB.
    dimension = range(1000, 5000, 1000)
    #dimension = range(5, 10, 10)
    #print dimension

    if 1:
        #plt.figure(1)
        #dense_square_test(dimension, is_binary=False)
        #plt.title('ASCII Dense Square Matrix Load Performance')

        plt.figure(2)
        dense_square_test(dimension, is_binary=True)
        plt.title('Binary Dense Square Matrix Load Performance')


    if make_plot:
        plt.show()