#!/usr/bin/env python

import sys
import numpy as np
from op4 import OP4
from pyNastran.op4.op4 import OP4 as pyOP4
from pyNastran.general.mathematics import print_matrix

for in_file in ['test/mat_b_dn.op4',
                'test/mat_b_s1.op4',
                'test/mat_b_s2.op4',
                'test/mat_t_dn.op4',
                'test/mat_t_s1.op4',
                'test/mat_t_s2.op4',
                'test/b_sample.op4',
                ]:
    try:
        op4fh = OP4(in_file, 'r')
    except:
        print('Failed to get header of %s, ignoring.' % (in_file))
        continue

    op4 = pyOP4()
    matrices = op4.readOP4(in_file)
    print('%s\n%s' % ('=' * 61, in_file))
    op4fh.print_header()
    for i in range(op4fh.nmat):
        a = op4fh.Load(nmat=1, skip=i)
        if a is None:
            print('Failed to get %d-th matrix' % i)
        else:
            Name = op4fh.name[i]

            Format, A = matrices[Name]
            error = 0.
            #print type(a),type(A)
            if isinstance(A, np.ndarray):
                pass
            else:  # sparse
                a = a.todense()
                A = A.todense()

            error = abs(a - A).max()

            if error > 0:
                print "Name = |%s|" % (Name)
                print('%s:' % op4fh.name[i])
                print "error[%s] = %s" % (Name, error)

                print "cOP4:"
                print print_matrix(a)

                print "pyOP4:"
                print print_matrix(A)

                print "diff:"
                print print_matrix(a - A)
                print '----------------------------'
                sys.exit('stopping')
